"""Internal helpers for solar_overlay; use build_pfss() and overlay() instead.

Keep this beside solar_overlay.py and the two datamine modules.
Importing this module does not run an observing workflow.
"""
from __future__ import annotations
import copy
import json
import re
import urllib.request
import warnings
from collections import OrderedDict
from contextlib import contextmanager
from functools import partial
from pathlib import Path
import matplotlib
import numpy as np
import astropy.units as u
from astropy.constants import R_sun
from astropy.coordinates import SkyCoord
from astropy.time import Time
from matplotlib.figure import Figure
import sunpy.map
from sunpy.coordinates import SphericalScreen, frames
from sunpy.util.exceptions import SunpyUserWarning
from sunpy.visualization import colormaps as _sunpy_cmaps  # noqa: F401 (registers cmaps)
from sunpy.visualization.colormaps import color_tables as ct
from datamine_calib_helper import AIA_WAVES, _series,_load_hmi, _matching_files,_require_eit_level1
import Multiwb_datamine as _datamine
import datamine_calib_helper as _helper

class _Layer:
	"""One image layer; load on demand, align references, and colour the result."""
	def __init__(self, label, kind, key, times, paths, loader, z, diff, dtype,
		reference_time, event_start, event_end):
		order = sorted(range(len(paths)), key=lambda i: times[i].mjd)
		self.times = [times[i] for i in order]
		self.paths = [paths[i] for i in order]
		self.label, self.kind, self.key = label, kind, key
		self.loader, self.z, self.diff, self.dtype = loader, z, diff, dtype
		self.event = [i for i, t in enumerate(self.times) if event_start <= t <= event_end]
		if not self.event:
			raise ValueError(f'{label}: no images inside the event interval')
		self.base = min(range(len(self.times)),
			key=lambda i: abs((self.times[i] - reference_time).to_value(u.s)))
		self.cache = OrderedDict()
		self.cmap = 'gray' if diff else ('RdBu_r' if kind == 'HMI' else
			(f'sdoaia{key}' if kind == 'AIA' else f'sohoeit{key}' if kind == 'EIT' else 'gray'))

	def index_for(self, time):
		# An external reference can never become a displayed event frame.
		return min(self.event, key=lambda i: abs((self.times[i] - time).to_value(u.s)))

	def processed(self, index, intensity=False):
		if intensity:
			return self.loader(self.paths[index])
		for i in (index, self.base if self.dtype == 'base' else index - 1):
			if not self.diff and i != index:
				continue
			if i < 0:
				raise FileNotFoundError(f'{self.label}: running difference needs an earlier image')
			if i not in self.cache:
				self.cache[i] = self.loader(self.paths[i])
			self.cache.move_to_end(i)
		while len(self.cache) > 3:
			self.cache.popitem(last=False)
		current = self.cache[index]
		if not self.diff:
			return current
		ref = self.cache[self.base if self.dtype == 'base' else index - 1]
		data = np.asarray(current.data, float) - _regrid(ref, current)
		meta = copy.deepcopy(current.meta)
		meta.pop('blank', None)
		meta['diff_ref'], meta['diff_typ'] = ref.date.isot, self.dtype
		meta['diff_dt'] = float((current.date - ref.date).to_value(u.s))
		return sunpy.map.Map(data, meta)

	def render(self, data, intensity=False):
		data = np.asarray(data, float)
		finite = np.isfinite(data)
		values = data[finite]
		if (self.diff and not intensity) or self.kind == 'HMI':
			limit = float(np.percentile(np.abs(values), 99)) if values.size else 1.
			limit = max(limit, np.finfo(float).eps)
			weight = np.clip(0.5 + data / (2 * limit), 0, 1)
		else:
			lo, hi = np.percentile(values, (1, 99.5)) if values.size else (0., 1.)
			weight = np.clip((data - lo) / max(hi - lo, np.finfo(float).eps), 0, 1)
			weight = np.arcsinh(10 * weight) / np.arcsinh(10.)
		weight = np.where(finite, weight, 0.)
		rgb = matplotlib.colormaps[self.cmap](weight)[..., :3]
		# Magnetic polarities and signed differences need equal opacity.
		return rgb, np.ones_like(weight) if (self.diff and not intensity) or self.kind == 'HMI' else weight, finite


@contextmanager
def _quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SunpyUserWarning)
        warnings.simplefilter("ignore", RuntimeWarning)
        yield


def _parse_time(t) -> Time:
    return t if isinstance(t, Time) else Time(t, scale="utc")


def _target_times(sT, eT, dT) -> Time:
    """sT, sT+dT, ... <= eT."""
    t0, t1 = _parse_time(sT), _parse_time(eT)
    if t1 <= t0:
        raise ValueError("eT must be later than sT")
    if not np.isfinite(dT) or dT <= 0:
        raise ValueError("dT must be > 0 (seconds)")
    n = int(np.floor((t1 - t0).to_value(u.s) / dT + 1e-9))
    return t0 + np.arange(n + 1) * float(dT) * u.s


def _get_gong_file(time: Time, cache_dir):
    """Nearest GONG zero-point-corrected synoptic radial map (mrzqs) to `time`."""
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    best, best_dt = None, np.inf
    for dday in (0, -1, 1):
        d = time + dday * u.day
        base = f"https://gong2.nso.edu/oQR/zqs/{d.strftime('%Y%m')}/mrzqs{d.strftime('%y%m%d')}/"
        try:
            html = urllib.request.urlopen(base, timeout=60).read().decode()
        except Exception:
            continue
        for fn in set(re.findall(r"(mrzqs\d{6}t\d{4}c\d{4}_\d{3}\.fits\.gz)", html)):
            ft = Time(f"20{fn[5:7]}-{fn[7:9]}-{fn[9:11]}T{fn[12:14]}:{fn[14:16]}:00")
            if abs((ft - time).to_value(u.s)) < best_dt:
                best, best_dt = (base, fn), abs((ft - time).to_value(u.s))
    if best is None:
        raise RuntimeError("Could not find a GONG synoptic map; pass gong_file=<path> to build_pfss.")
    out = cache / best[1]
    if not out.exists():
        urllib.request.urlretrieve(best[0] + best[1], out)
    return str(out)


def _load_pfss(directory):
	folder = Path(directory).expanduser().resolve()
	metadata = json.loads((folder / 'metadata.json').read_text())
	if metadata.get('format_version') != 1:
		raise ValueError('Unsupported saved PFSS format')
	with np.load(folder / 'field_lines.npz', allow_pickle=False) as saved:
		out = {k: saved[k].copy() for k in ('lon', 'lat', 'r', 'seg', 'pol')}
		out['time'] = Time(str(saved['time']))
		out['rss'] = float(saved['rss'])
	out['outdir'] = str(folder)
	return out


def _nan_sep(x, y, seg):
    br = np.flatnonzero(np.diff(seg) != 0) + 1
    return np.insert(x, br, np.nan), np.insert(y, br, np.nan)


def _pfss_to_fov(tmap, lines):
	"""Return radial continuations of open PFSS lines to this map's FoV edge.

	The source-surface longitude/latitude stays fixed in Carrington coordinates.
	The projected ray is intersected with the actual pixel rectangle, so cropped,
	rotated, C2 and C3 maps need no manually chosen outer radius. Closed loops and
	traces that did not reach the source surface are never lengthened.
	"""
	result = {k: [] for k in ("lon", "lat", "r", "seg", "pol")}
	obs = tmap.observer_coordinate
	hpc_frame = tmap.coordinate_frame
	cx, cy = tmap.world_to_pixel(SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=hpc_frame))
	center = np.array([cx.value, cy.value])
	ny, nx = tmap.data.shape
	rss = float(lines.get("rss", np.max(lines["r"])))
	for seg in np.unique(lines["seg"]):
		idx = np.flatnonzero(lines["seg"] == seg)
		if not len(idx) or lines["pol"][idx[0]] == 0:
			continue
		i = idx[-1] if lines["r"][idx[-1]] >= lines["r"][idx[0]] else idx[0]
		r0 = float(lines["r"][i])
		if not np.isclose(r0, rss, rtol=0.02):
			continue
		lon, lat = float(lines["lon"][i]), float(lines["lat"][i])
		c = SkyCoord(lon * u.deg, lat * u.deg, r0 * R_sun,
			frame=frames.HeliographicCarrington, observer=obs, obstime=hpc_frame.obstime)
		px, py = tmap.world_to_pixel(c.transform_to(hpc_frame))
		direction = np.array([px.value, py.value]) - center
		if not np.all(np.isfinite(direction)) or np.linalg.norm(direction) < 1e-8:
			continue  # An exactly line-of-sight radial has no projected length.
		enter, leave = 1., np.inf
		for pos, delta, limit in zip(center, direction, (nx, ny)):
			if abs(delta) < 1e-12:
				if not -0.5 <= pos <= limit - 0.5:
					leave = -np.inf
				continue
			a, b = sorted(((-0.5 - pos) / delta, (limit - 0.5 - pos) / delta))
			enter, leave = max(enter, a), min(leave, b)
		if not np.isfinite(leave) or leave <= enter:
			continue
		edge = center + leave * direction
		e = tmap.pixel_to_world(edge[0] * u.pix, edge[1] * u.pix)
		# Solve the perspective projection for the radius at the FoV edge.
		n = c.transform_to(frames.Heliocentric(observer=obs, obstime=hpc_frame.obstime)).cartesian.xyz.to_value(R_sun) / r0
		q = np.array([np.tan(e.Tx.to_value(u.rad)),
			np.tan(e.Ty.to_value(u.rad)) / np.cos(e.Tx.to_value(u.rad))])
		j = int(np.argmax(np.abs(n[:2])))
		denom = n[j] + q[j] * n[2]
		if abs(denom) < 1e-14:
			continue
		r_edge = obs.radius.to_value(R_sun) * q[j] / denom
		if not np.isfinite(r_edge) or r_edge <= r0:
			continue  # A far-side ray can asymptote inside the FoV; do not invent an exit.
		r = np.geomspace(r0, r_edge, 128)
		for key, value in (("lon", np.full(len(r), lon)), ("lat", np.full(len(r), lat)),
			("r", r), ("seg", np.full(len(r), seg)), ("pol", np.full(len(r), lines["pol"][i]))):
			result[key].append(value)
	return {k: np.concatenate(v) if v else np.array([]) for k, v in result.items()}


def _closed_loop_limit(value):
	"""Validate the optional cap on the number of plotted closed loops."""
	try:
		value = float(value)
		if np.isnan(value):
			return None
		if np.isfinite(value) and value >= 0 and value.is_integer():
			return int(value)
	except (TypeError, ValueError):
		pass
	raise ValueError('N_Bloop must be np.nan or a non-negative integer.')


def _draw_pfss(ax, tmap, lines, lw=0.6, alpha=0.9,
              colors=None, N_Bloop=np.nan):
    """
    Draw PFSS field lines on an axes whose pixel grid is `tmap` (a sunpy Map).
    Field lines co-rotate with the Sun: they are re-projected for tmap.date.
    Open +/-: red/blue, closed: white. Segments hidden behind the solar disk are removed.
    Open lines always continue radially to the FoV edge: solid = computed PFSS,
    dashed = radial continuation beyond the source surface. Closed loops retain
    their computed geometry. The PFSS source surface itself is not enlarged.
    N_Bloop: np.nan plots all closed loops; an integer keeps at most that many
    visible closed loops, ranked by maximum radius (highest apex first). Zero
    hides closed loops. Open-field lines and continuations are unaffected.
    """
    colors = colors or {1: "#ff4d4d", -1: "#4da6ff", 0: "#ffffff"}
    extension = _pfss_to_fov(tmap, lines)
    extra = np.r_[np.zeros(len(lines["r"]), dtype=bool), np.ones(len(extension["r"]), dtype=bool)]
    lines = {k: np.r_[lines[k], extension[k]] for k in ("lon", "lat", "r", "seg", "pol")}
    obs = tmap.observer_coordinate
    c = SkyCoord(lines["lon"] * u.deg, lines["lat"] * u.deg, lines["r"] * R_sun,
                 frame=frames.HeliographicCarrington, observer=obs, obstime=tmap.coordinate_frame.obstime)
    with _quiet():
        hpc = c.transform_to(tmap.coordinate_frame)
        px, py = tmap.world_to_pixel(hpc)
    px, py = np.asarray(px.value, float), np.asarray(py.value, float)
    rho = np.hypot(hpc.Tx.to_value(u.arcsec), hpc.Ty.to_value(u.arcsec))
    dobs = obs.radius.to_value(u.m)
    rs = tmap.rsun_meters.to_value(u.m)
    behind = (rho < tmap.rsun_obs.to_value(u.arcsec)) & \
             (hpc.distance.to_value(u.m) > np.sqrt(max(dobs ** 2 - rs ** 2, 0)))
    px[behind], py[behind] = np.nan, np.nan
    limit = _closed_loop_limit(N_Bloop)
    if limit is not None:
        ny, nx = tmap.data.shape
        closed = lines['pol'] == 0
        visible = np.isfinite(px) & np.isfinite(py) & (px >= -0.5) & (px <= nx - 0.5) & (py >= -0.5) & (py <= ny - 0.5)
        ids = np.unique(lines['seg'][closed & visible])
        ids = sorted(ids, key=lambda seg: (-np.max(lines['r'][lines['seg'] == seg]), seg))[:limit]
        hide = closed & ~np.isin(lines['seg'], ids)
        px[hide], py[hide] = np.nan, np.nan
    for pol, col in colors.items():
        for is_extra, style in ((False, "-"), (True, "--")):
            s = (lines["pol"] == pol) & (extra == is_extra)
            if s.any():
                x, y = _nan_sep(px[s], py[s], lines["seg"][s])
                ax.plot(x, y, color=col, lw=lw, alpha=alpha, linestyle=style)
    ny, nx = tmap.data.shape
    ax.set_xlim(-0.5, nx - 0.5)
    ax.set_ylim(-0.5, ny - 0.5)


def _target_map(m, npix):
    """Common display grid: the lowest layer's map (down-sampled if it is huge)."""
    if max(m.data.shape) > npix:
        m = m.resample([npix, npix] * u.pix)
    if np.abs(m.rotation_matrix - np.eye(2)).max() > 1e-3:  # e.g. HMI (CROTA2=180)
        m = m.rotate()
    return m


def _regrid(m, tmap):
    """Reproject map `m` onto the pixel grid / WCS of `tmap` (NaN outside m's FOV)."""
    ratio = abs(float((tmap.scale.axis1 / m.scale.axis1).to_value(u.one)))
    f = 1
    while f * 2 <= ratio and m.data.shape[0] % (f * 2) == 0 and m.data.shape[1] % (f * 2) == 0:
        f *= 2
    with _quiet():
        if f > 1:
            m = m.superpixel([f, f] * u.pix, func=np.nanmean)
        with SphericalScreen(m.observer_coordinate, only_off_disk=True):
            out = m.reproject_to(tmap.wcs)
    return out.data


def _composite(parts, alpha):
    """Alpha-over compositing; first layer opaque, others use alpha * intensity weight."""
    h, w = parts[0][0].shape[:2]
    canvas = np.zeros((h, w, 3))
    for i, (rgb, wt, fin) in enumerate(parts):
        al = np.where(fin, 1.0 if i == 0 else alpha * wt, 0.0)[..., None]
        canvas = np.where(fin[..., None], rgb, 0.) * al + canvas * (1 - al)
    return canvas


def _render_frames(layers, frame_times, outdir, start=0, pfss_lines=None, alpha=0.6,
                   npix=1024, dpi=150, figsize=7.0, prefix="overlay", N_Bloop=np.nan, tmap=None, first_frame_intensity=False):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    layers = sorted(layers, key=lambda L: L.z)
    written = []
    previous = None
    template = sunpy.map.Map(tmap) if isinstance(tmap, (str, Path)) else tmap
    for k in range(start, len(frame_times)):
        t = frame_times[k]
        idx = [L.index_for(t) for L in layers]
        if idx[0] == previous:
            continue  # Do not repeat an exposure at a finer requested cadence.
        previous = idx[0]
        t = layers[0].times[idx[0]]
        idx = [L.index_for(t) for L in layers]
        intensity = [first_frame_intensity and not written and L.diff and L.dtype == 'running' for L in layers]
        maps = [L.processed(i, intensity=raw) for L, i, raw in zip(layers, idx, intensity)]
        tmap = _target_map(template if template is not None else maps[0], npix)
        datas = [_regrid(m, tmap) for m in maps] if template is not None else [tmap.data] + [_regrid(m, tmap) for m in maps[1:]]
        parts = [L.render(d, intensity=raw) for L, d, raw in zip(layers, datas, intensity)]
        img = _composite(parts, alpha)

        fig = Figure(figsize=(figsize, figsize), dpi=dpi)
        ax = fig.add_subplot(projection=tmap)
        ax.imshow(img, origin="lower", interpolation="nearest")
        try:
            tmap.draw_limb(axes=ax, color="white", linewidth=0.6, alpha=0.7)
        except Exception:
            pass
        if pfss_lines is not None:
            _draw_pfss(ax, tmap, pfss_lines, N_Bloop=N_Bloop)
        ax.text(0.02, 0.98, t.strftime("%Y-%m-%d %H:%M:%S"), transform=ax.transAxes,
                color="white", va="top", ha="left", fontsize=11,
                bbox=dict(facecolor="black", alpha=0.5, pad=2, lw=0))
        tags = [f"{L.label}" + (" (intensity)" if raw else f" ({L.dtype} diff)" if L.diff else "")
                + f" [{m.date.strftime('%Y-%m-%d %H:%M:%S')} UTC]"
                + (" preview" if "ov_cal" in m.meta else "")
                for L, i, m, raw in zip(layers, idx, maps, intensity)]
        if template is not None:
            tags.append("Fixed WCS: " + template.reference_date.strftime("%Y-%m-%d %H:%M:%S UTC"))
        tags += ["PFSS + radial continuation (dashed)" for _ in [0] if pfss_lines is not None]
        for j, tag in enumerate(tags):
            ax.text(0.98, 0.98 - 0.04 * j, tag, transform=ax.transAxes, color="white",
                    va="top", ha="right", fontsize=8,
                    bbox=dict(facecolor="black", alpha=0.5, pad=1.5, lw=0))
        out = outdir / f"{prefix}_{t.strftime('%Y%m%d_%H%M%S%f')}.png"
        fig.savefig(out, dpi=dpi)
        fig.clear()
        written.append(str(out))
    return written


def _parse_dataset(ds):
    s = ds.strip().upper().replace(" ", "")
    if s == "PFSS":
        return "PFSS", None
    m = re.fullmatch(r"LASCO/(C2|C3)", s)
    if m:
        return "LASCO", m.group(1)
    m = re.fullmatch(r"SDO/AIA(\d+)", s)
    if m and int(m.group(1)) in AIA_WAVES:
        return "AIA", int(m.group(1))
    m = re.fullmatch(r"SOHO/EIT(171|195|284|304)", s)
    if m:
        return "EIT", int(m.group(1))
    if s == "SDO/HMI":
        return "HMI", None
    raise ValueError(f"Unsupported dataset '{ds}'. Use LASCO/C2, LASCO/C3, SDO/AIA<wave>, SDO/HMI, pfss")

