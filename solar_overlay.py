"""
solar_overlay.py -- image overlays, difference frames and PFSS.

Read finished FITS Maps saved by Multiwb_datamine functions in user-selected
subfolders. Public API: build_pfss, overlay and pfss_context_maps. Acquisition and instrument
calibration belong to Multiwb_datamine: get_LASCO_imgs, calib_LASCO,
get_SDO_imgs, get_EIT_imgs. No download is started by overlay unless explicitly
requested with download=True. Importing either module starts no workflow.

All sT/eT values are UTC strings 'YYYY-MM-DDThh:mm:ss'; dT is seconds.
Public API and inputs
---------------------
build_pfss(sT, eT, outdir, rss=2.5, nrho=35, r_seed=1.05,
           n_lon=36, n_lat=18, gong_file=None, cache_dir='./solar_data/gong')
    Required inputs are the event interval and output directory. Computes one
    STATIC 3D PFSS model from the GONG boundary nearest the midpoint. Saves
    model_3d.npz (3D vector field and axes), field_lines.npz, boundary.fits and
    metadata.json. Returns traced-line arrays and the output directory.
    rss: source surface [solar radii]; nrho: radial cells; r_seed: seed radius
    [solar radii], 1 < r_seed < rss; n_lon/n_lat: seed-grid counts. gong_file:
    optional local GONG synoptic FITS; cache_dir: downloaded boundary cache.

pfss_context_maps(sT, mission='SDO', instrument=None, wave=None, *, basedir='./',
                  NOAA_number=None, ...)
    Saves a full-Sun and an AR PFSS map at time sT in basedir (PFSS_fullSun_*.png,
    PFSS_AR<NOAA>_*.png), on an AIA/HMI/EUVI image nearest sT or, with
    instrument=None, on a blank disk with a Stonyhurst grid. The AR map needs
    Moving_structure_frames (noaa_region, pfss_projection) on the Python path.

 overlay(datasets, sT, eT, dT, diffs=None, diff_types=None, outdir='./',
         alpha=0.6, cache_dir='./solar_data', npix=1024, dpi=150, figsize=7.0,
         pfss_kwargs=None, local_files=None,
         N_Bloop=np.nan, need_LASCO_calib=False,
         image_dirs=None, pfss_dir=None, tmap=None, download=False,
         base_diff_DT=None, reference_search=3600, first_frame_intensity=False)
    Returns PNG path list. datasets: ordered names such as ['LASCO/C2',
    'SDO/AIA171', 'SOHO/EIT195', 'SDO/HMI', 'pfss']; choose any image subset.
    image_dirs: one FITS directory, a list of directories, or {dataset: directory
    or list of directories}. Searches recursively and selects each instrument/
    wavelength from FITS headers. No images are downloaded in this mode.
    pfss_dir: saved build_pfss output folder; required with 'pfss' in directory
    mode. tmap: SunPy Map or FITS path defining display FoV/WCS; None uses the
    bottom image layer. No axes or field-line arrays need to be supplied.
    diffs: per-dataset booleans (missing entries False, PFSS always False).
    diff_types: corresponding 'running'/'base' strings; base uses base_diff_DT.
    outdir: PNG destination; alpha: upper-image opacity [0,1]; npix: longest
    output-grid dimension; dpi: resolution; figsize: square size in inches.
    N_Bloop: NaN keeps all closed loops; nonnegative integer keeps the highest
    N visible closed loops (0 = open lines only). Red/blue: outward/inward open
    lines; white: closed. Open lines continue radially to FoV beyond rss,
    dashed to distinguish the continuation from the computed PFSS solution.
    need_LASCO_calib=False: quick looks; warn if returned data are uncalibrated;
    True uses get_LASCO_imgs for downloads (Level-1 -> JP2 -> DN/s), or
    calib_LASCO for local FITS (optional SolarSoft -> Level-1 -> JP2).
    Verified local Level-1 inputs pass through; local calibration can retrieve
    replacements and may attempt the optional IDL backend.
    Acquisition is disabled by default. Set download=True explicitly to call
    Multiwb_datamine when local image inputs are omitted.
    cache_dir stores downloads; pfss_kwargs forwards optional build_pfss
    settings when no pfss_dir is given. local_files: alternative explicit
    {dataset: [FITS paths]} input; missing keys raise unless download=True.

Requirements: pip install 'sunpy[net,map,timeseries]>=8' reproject matplotlib astropy numpy
pfsspy aiapy radiospectra pandas pillow beautifulsoup4. Keep the supplied
Multiwb_datamine.py on the Python import path. Downloaded verified Level-1 products need no IDL; optional local SolarSoft
calibration requires IDL. Keep solar_overlay_helper.py, Multiwb_datamine.py
and datamine_calib_helper.py beside this file.
"""
from __future__ import annotations
import copy
import json
import re
import warnings
from pathlib import Path
import numpy as np
import astropy.units as u
from astropy.constants import R_sun
from astropy.coordinates import SkyCoord
from astropy.time import Time
import sunpy.map


__all__ = ["build_pfss", "overlay"]

def build_pfss(sT, eT, outdir, rss=2.5, nrho=35, r_seed=1.05, n_lon=36, n_lat=18,
               gong_file=None, cache_dir="./solar_data/gong"):
    """Build and save a static 3D PFSS model for the UTC interval sT--eT.

    Required: sT/eT UTC strings and outdir (model output directory). The boundary
    is the GONG synoptic map nearest the interval midpoint. Returns traced-line
    arrays plus model paths; saves field_lines.npz, model_3d.npz, boundary.fits
    and metadata.json. This is one static model, not time-dependent MHD.

    Optional: rss=2.5 source-surface radius [solar radii]; nrho=35 radial cells;
    r_seed=1.05 tracing seed radius [solar radii], between 1 and rss;
    n_lon=36/n_lat=18 seed counts; gong_file=None local GONG synoptic FITS
    instead of download; cache_dir='./solar_data/gong' boundary download cache.
    The saved 3D grid contains Bphi, Btheta, Br [G] and Carrington longitude,
    sine-latitude and radius axes. Beyond rss only radial plotting continuation
    is used; it is not part of the saved PFSS solution.
    """
    import solar_overlay_helper as _overlay_helper
    import pfsspy
    from pfsspy import tracing

    t0, t1 = _overlay_helper._parse_time(sT), _overlay_helper._parse_time(eT)
    if t1 < t0:
        raise ValueError('eT must not precede sT')
    if not 1 < r_seed < rss:
        raise ValueError('Require 1 < r_seed < rss')
    time = t0 + (t1 - t0) / 2
    fn = gong_file or _overlay_helper._get_gong_file(time, cache_dir)
    with _overlay_helper._quiet():
        gmap = sunpy.map.Map(fn)
        data = np.nan_to_num(gmap.data)
        meta = copy.deepcopy(gmap.meta)
        # GONGSynopticMap already converts its native sine-latitude scale.
        # Only fill missing unit keywords required directly by pfsspy.
        meta.setdefault("cunit1", str(gmap.spatial_units.axis1))
        meta.setdefault("cunit2", str(gmap.spatial_units.axis2))
        if str(meta.get("bunit", "")).lower() == "gauss":
            meta["bunit"] = "G"
        _overlay_helper._helper._add_rsun_ref(meta)  # GONG headers lack RSUN_REF; write it explicitly
        gmap = sunpy.map.Map(data - data.mean(), meta)
    pout = pfsspy.pfss(pfsspy.Input(gmap, nrho, rss))

    lon = np.deg2rad(np.linspace(0, 360, n_lon, endpoint=False))
    lat = np.arcsin(np.linspace(-1, 1, n_lat + 2)[1:-1])  # equal-area latitudes
    lo, la = np.meshgrid(lon, lat)
    seeds = SkyCoord(lo.ravel() * u.rad, la.ravel() * u.rad, r_seed * R_sun,
                     frame=pout.coordinate_frame)
    try:
        tracer = tracing.FortranTracer()
    except Exception:
        tracer = tracing.PythonTracer()
    flines = tracer.trace(seeds, pout)

    L = {k: [] for k in ("lon", "lat", "r", "seg", "pol")}
    for i, fl in enumerate(flines):
        c = fl.coords
        if len(c) < 2:
            continue
        L["lon"].append(c.lon.to_value(u.deg))
        L["lat"].append(c.lat.to_value(u.deg))
        L["r"].append(c.radius.to_value(R_sun))
        L["seg"].append(np.full(len(c), i))
        L["pol"].append(np.full(len(c), int(getattr(fl, "polarity", 0))))
    if not L["lon"]:
        raise RuntimeError("PFSS tracing returned no usable field lines; check boundary and seeds.")
    out = {k: np.concatenate(v) for k, v in L.items()}
    out["time"] = time
    out["rss"] = float(rss)
    folder = Path(outdir).expanduser().resolve()
    folder.mkdir(parents=True, exist_ok=True)
    gmap.save(folder / 'boundary.fits', overwrite=True)
    np.savez_compressed(folder / 'field_lines.npz', **{k: out[k] for k in L},
                        time=time.isot, rss=float(rss))
    bg = pout.bg.to_value(u.G)
    np.savez_compressed(folder / 'model_3d.npz', Bphi=bg[..., 0],
                        Btheta=bg[..., 1], Br=bg[..., 2],
                        longitude_deg=(np.rad2deg(pout.grid.pg) + float(gmap.meta['crval1']) - 180) % 360,
                        sin_latitude=pout.grid.sg, radius_rsun=np.exp(pout.grid.rg))
    metadata = dict(format_version=1, sT=t0.isot, eT=t1.isot, model_time=time.isot,
                    boundary_time=gmap.date.isot, boundary_source=str(fn), rss=float(rss),
                    nrho=int(nrho), r_seed=float(r_seed), n_lon=int(n_lon), n_lat=int(n_lat),
                    field_unit='G', components=['Bphi', 'Btheta', 'Br'],
                    coordinate_frame='heliographic_carrington',
                    component_convention='physical spherical: phi=east, theta=south, r=outward',
                    approximation='static midpoint PFSS; radial continuation beyond rss is display only')
    (folder / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    out['outdir'] = str(folder)
    return out


def overlay(datasets, sT, eT, dT, diffs=None, diff_types=None, outdir="./",
            alpha=0.6, cache_dir="./solar_data", npix=1024, dpi=150, figsize=7.0,
            pfss_kwargs=None, local_files=None,
            N_Bloop=np.nan, need_LASCO_calib=False,
            image_dirs=None, pfss_dir=None, tmap=None, download=False,
            base_diff_DT=None, reference_search=3600, first_frame_intensity=False):
    """
    Write overlay PNG frames (timestamp in the top-left corner) to `outdir`.

    download: False by default; True explicitly enables missing-image retrieval
        through Multiwb_datamine. Directory mode never downloads missing inputs.
    image_dirs: FITS directory, list of directories, or dataset-to-directories
        dictionary; recursive header selection, no acquisition in this mode.
    pfss_kwargs: optional build_pfss settings, for example
        {'rss':2.5, 'nrho':35, 'r_seed':1.05, 'n_lon':36, 'n_lat':18}.
        r_seed is the tracing-start radius [solar radii]; n_lon/n_lat count
        seed longitudes/latitudes. None uses defaults; ignored with pfss_dir.
    pfss_dir: directory saved by build_pfss(sT, eT, outdir), reused as-is.
    tmap: target SunPy Map/FITS path specifying fixed display FoV and WCS;
        all images and PFSS are transformed to its actual reference epoch.
        None uses each frame's bottom image and therefore its observation time.
    SOHO/EIT171, SOHO/EIT195, SOHO/EIT284 and SOHO/EIT304 are supported.
    Low-level drawing is private; callers need neither axes nor trace arrays.

    datasets   : list from {'LASCO/C2','LASCO/C3','SDO/AIA94','SDO/AIA171',...,'SDO/HMI','pfss'}
    sT, eT, dT : start, end ("YYYY-MM-DDThh:mm:ss") and cadence [s]
    diffs      : booleans, one per dataset (missing entries -> False; pfss always False)
    diff_types : 'running' | 'base' per dataset ('' where diff is False).
                 'base' uses base_diff_DT (defaults to sT).
    first_frame_intensity: True displays the first running-difference frame as
                 intensity; later frames subtract their preceding selected exposure.
                 No image before sT is needed in this mode. Default False:
                 try a pre-sT image first, falling back to intensity if absent.
    base_diff_DT: UTC reference time, or {dataset: UTC reference time}.
                 Before sT: latest image <= base_diff_DT. After eT: earliest
                 image >= base_diff_DT. Inside [sT,eT]: nearest event image.
    reference_search: maximum search distance from an external reference time,
                 in seconds (default 3600). Missing references raise an error;
                 the event boundary is never substituted.
    alpha      : opacity of upper layers. EUV intensity, HMI and signed differences
                 have uniform opacity; other intensity layers retain brightness weighting.
    local_files: optional dict {dataset: list of FITS paths}, e.g. files already
                 downloaded by Multiwb_datamine.find_coobs(). Missing keys raise unless download=True. Prepared AIA inputs are not calibrated twice.
    N_Bloop    : PFSS only; np.nan = all closed loops, integer = keep the highest
                 N visible closed loops, 0 = open lines only.
    Local LASCO files are first reduced to one product family by the route
                 Level-1 -> Helioviewer -> Level-0.5 (need_LASCO_calib=True) or
                 Level-1 -> Level-0.5 -> Helioviewer (False); see _route_family.
    need_LASCO_calib: False (default) prepares available LASCO as DN/s previews;
                 True uses get_LASCO_imgs for missing inputs (Level-1 -> JP2 ->
                 DN/s), or calib_LASCO for local inputs (optional SolarSoft ->
                 Level-1 -> JP2). This permits remote replacement of local raw
                 inputs; verified local Level-1 inputs pass through unchanged.
                 Local calibration may attempt the optional IDL backend.

    Layer stacking order (bottom -> top): LASCO C3, LASCO C2, HMI, AIA bands (in the order
    given), PFSS field lines. The bottom layer defines the display grid.
    Frames run sT ... eT. Running differences try a predecessor before sT, falling back to intensity if absent.
    reference_search bounds that predecessor search too. With
    first_frame_intensity=True, the first frame is intensity and needs no predecessor.
    Times used for each layer are printed in its label; the top-left time is the
    actual bottom-image time. Event images stay inside [sT, eT]; references may be outside.
    Differences align the reference WCS to the current image before subtraction;
    no differential-rotation correction is applied. Returns the list of PNG paths.
    """
    import solar_overlay_helper as _overlay_helper
    if not isinstance(download, (bool, np.bool_)):
        raise ValueError('download must be True or False')
    if image_dirs is None and local_files is None and not download:
        raise ValueError('Supply image_dirs or local_files containing finished FITS; use download=True explicitly for acquisition via Multiwb_datamine.')
    if not datasets:
        raise ValueError('datasets must contain at least one image dataset')
    if image_dirs is not None:
        if local_files is not None:
            raise ValueError('Use image_dirs or local_files, not both')
        local_files = {}
        for name in datasets:
            kind, key = _overlay_helper._parse_dataset(name)
            if kind == 'PFSS':
                continue
            directories = image_dirs.get(name) if isinstance(image_dirs, dict) else image_dirs
            if directories is None:
                raise ValueError(f'Missing image directory for {name}')
            if isinstance(directories, (str, Path)):
                directories = [directories]
            local_files[name] = _overlay_helper._matching_files(directories, kind, key)
        if any(_overlay_helper._parse_dataset(d)[0] == 'PFSS' for d in datasets) and pfss_dir is None:
            raise ValueError('Directory mode requires pfss_dir for PFSS; call build_pfss first')
    if pfss_dir is not None:
        metadata = json.loads((Path(pfss_dir) / 'metadata.json').read_text())
        if _overlay_helper._parse_time(sT) < Time(metadata['sT']) or _overlay_helper._parse_time(eT) > Time(metadata['eT']):
            warnings.warn('Requested frames extend outside the saved PFSS interval; using the same static model.')
    diffs = list(diffs or [])
    diff_types = list(diff_types or [])
    if not np.isfinite(alpha) or not 0 <= alpha <= 1:
        raise ValueError("alpha must be between 0 and 1")
    if len(diffs) > len(datasets) or len(diff_types) > len(datasets):
        raise ValueError("diffs / diff_types are longer than datasets")
    specs = [_overlay_helper._parse_dataset(d) for d in datasets]
    if any(k == 'PFSS' for k, _ in specs):
        _overlay_helper._closed_loop_limit(N_Bloop)
    if not isinstance(need_LASCO_calib, (bool, np.bool_)):
        raise ValueError('need_LASCO_calib must be True or False.')
    if len({(k, v) for k, v in specs}) != len(specs):
        raise ValueError("Duplicate entries in datasets")
    if all(k == "PFSS" for k, _ in specs):
        raise ValueError("Select at least one image dataset (LASCO, AIA or HMI) to draw on")

    frame_times = _overlay_helper._target_times(sT, eT, dT)
    t_start = frame_times[0]
    t_end = _overlay_helper._parse_time(eT)
    if not np.isfinite(reference_search) or reference_search <= 0:
        raise ValueError('reference_search must be positive seconds')

    layers, want_pfss, any_diff = [], False, False
    uncalibrated_lasco = False
    n_aia = 0
    for i, ((kind, key), name) in enumerate(zip(specs, datasets)):
        if kind == "PFSS":
            want_pfss = True
            continue
        d = bool(diffs[i]) if i < len(diffs) else False
        typ = (diff_types[i] if i < len(diff_types) else "").strip().lower()
        if d and typ not in ("running", "base"):
            raise ValueError(f"diff_types[{i}] must be 'running' or 'base' for '{name}'")
        any_diff |= d
        print(f"[overlay] loading {name} ...")
        files = (local_files or {}).get(name)
        if files is None and not download:
            raise ValueError(f'Missing local FITS for {name}; retrieve them with Multiwb_datamine first.')
        if kind == "LASCO":
            lasco_outdir = Path(cache_dir) / 'prepared' / ('LASCO_' + key)
            if files is not None:
                # One product family by the route Level-1 -> Helioviewer -> Level-0.5
                # (Level-0.5 before Helioviewer when need_LASCO_calib=False).
                files = _overlay_helper._route_family(files, sT, eT, dT, bool(need_LASCO_calib))
            if files is None:
                prepared = _overlay_helper._datamine.get_LASCO_imgs(sT, eT, dT, key,
                    outdir=lasco_outdir, need_calib=need_LASCO_calib)
            else:
                _, selected = _overlay_helper._series('LASCO', key, sT, eT, dT,
                    cache_dir, files=files)
                loaded = [_overlay_helper._helper.load_map(path) for path in selected]
                if need_LASCO_calib and not all(_overlay_helper._helper._lasco_science_ready(m) or
                        m.meta.get('ov_hv', False) for m in loaded):
                    prepared = _overlay_helper._datamine.calib_LASCO(list(selected),
                        instr=key, outdir=lasco_outdir)
                elif need_LASCO_calib:
                    # Level-1 or Helioviewer family: final already; no recalibration.
                    prepared = [_overlay_helper._helper._prepare_lasco(m, key) for m in loaded]
                else:
                    prepared = [_overlay_helper._helper._prepare_lasco(_overlay_helper._helper.load_map(path), key)
                                for path in selected]
                    prepared = _overlay_helper._helper._save_maps(prepared,
                        lasco_outdir / 'maps', 'LASCO_' + key)
                _overlay_helper._helper._lasco_notice(prepared)
            prepared = _overlay_helper._helper._save_maps(prepared, lasco_outdir / 'maps', 'LASCO_' + key)
            uncalibrated_lasco |= any(not _overlay_helper._helper._lasco_science_ready(m) for m in prepared)
            times = [m.date for m in prepared]
            paths = [Path(cache_dir).expanduser().resolve() / 'prepared' / ('LASCO_' + key) /
                     f"maps/LASCO_{key}_{m.date.strftime('%Y%m%dT%H%M%S%f')}.fits" for m in prepared]
            loader = _overlay_helper._helper.load_map
            z = 0 if key == "C3" else 1
        elif kind == "AIA":
            if files is None:
                _overlay_helper._datamine.get_SDO_imgs(sT, eT, dT, instruments=['AIA'],
                             waves=[[key]], nproc=1,
                             outdir=Path(cache_dir) / 'prepared' / f'AIA{key}')
                files = _overlay_helper._matching_files([Path(cache_dir) / 'prepared' / f'AIA{key}'], kind, key)
            times, paths = _overlay_helper._series("AIA", key, sT, eT, dT, cache_dir, files=files)
            for path in paths:
                if _overlay_helper._is_preview_image(_overlay_helper._helper.load_map(path),path):
                    warnings.warn(f'{path}: using a Helioviewer display preview, not calibrated AIA counts.')
                    continue
                meta = sunpy.map.Map(path).meta
                if not (float(meta.get('lvl_num', 0)) >= 1.5 or
                        any(v in str(meta.get('series', '')).lower() for v in ('lev1.5', 'lev15'))):
                    raise ValueError('overlay requires prepared AIA Level-1.5 FITS')
            loader, z = _overlay_helper._helper.load_map, 3 + n_aia
            n_aia += 1
        elif kind == 'EIT':
            if files is None:
                _overlay_helper._datamine.get_EIT_imgs(sT, eT, dT, key, outdir=Path(cache_dir) / f'EIT{key}')
                files = _overlay_helper._matching_files([Path(cache_dir) / f'EIT{key}'], kind, key)
            times, paths = _overlay_helper._series('EIT', key, sT, eT, dT, cache_dir, files=files)
            for path in paths:
                _overlay_helper._require_eit_level1(_overlay_helper._helper.load_map(path), path)
            loader, z = _overlay_helper._helper.load_map, 3 + n_aia
            n_aia += 1
        else:
            if files is None:
                _overlay_helper._datamine.get_SDO_imgs(sT, eT, dT, instruments=['HMI'],
                             waves=[['magnetogram']], nproc=1, outdir=Path(cache_dir) / 'prepared' / 'HMI')
                files = _overlay_helper._matching_files([Path(cache_dir) / 'prepared' / 'HMI'], kind, key)
            times, paths = _overlay_helper._series("HMI", None, sT, eT, dT, cache_dir, files=files)
            loader, z = _overlay_helper._load_hmi, 2
        reference_time = t_start
        if d and not (typ == 'running' and first_frame_intensity):
            requested = base_diff_DT.get(name) if isinstance(base_diff_DT, dict) else base_diff_DT
            if typ == 'running':
                requested = sT
            base = _overlay_helper._parse_time(requested if requested is not None else sT)
            outside = base < t_start or base > t_end
            lo = base - reference_search * u.s if base < t_start else base
            hi = base + reference_search * u.s if base > t_end else base
            if not outside:
                lo, hi = t_start, t_end
            if typ == 'running':
                lo, hi = t_start - reference_search * u.s, t_start - 1e-6 * u.s
            # Inspect original local files too: event cadence selection can omit
            # the closest reference or discard all files outside the event.
            candidates = [(t, p) for t, p in zip(times, paths) if lo <= t <= hi]
            for path in files or []:
                time = sunpy.map.Map(path).date
                if lo <= time <= hi:
                    candidates.append((time, path))
            if not candidates and download:
                refdir = Path(cache_dir) / 'base_reference' / name.replace('/', '_')
                # Query only this bounded window, never base_diff_DT--sT/eT.
                try:
                    if kind == 'LASCO':
                        maps = _overlay_helper._datamine.get_LASCO_imgs(
                            lo.isot, hi.isot, dT, key, outdir=refdir,
                            need_calib=need_LASCO_calib)
                    elif kind == 'EIT':
                        maps = _overlay_helper._datamine.get_EIT_imgs(
                            lo.isot, hi.isot, dT, key, outdir=refdir)
                    else:
                        maps = _overlay_helper._datamine.get_SDO_imgs(
                            lo.isot, hi.isot, dT, instruments=[kind],
                            waves=[[key if kind == 'AIA' else 'magnetogram']],
                            outdir=refdir, nproc=1)[0]
                except (OSError, RuntimeError, ValueError) as exc:
                    if typ != 'running':
                        raise
                    warnings.warn(f'{name}: predecessor retrieval failed: {exc}')
                    maps = []
                for m in maps:
                    m = _overlay_helper._helper.load_map(m)
                    if lo <= m.date <= hi:
                        candidates.append((m.date, m))
            if not candidates and typ == 'running':
                first_frame_intensity = True
                warnings.warn(f'{name}: no pre-sT image found; first frame will show intensity.')
            elif not candidates:
                raise FileNotFoundError(
                    f'{name}: no {typ} reference in {lo.isot} -- {hi.isot}; '
                    'save an eligible image locally or enable download=True.')
            if candidates:
                reference_time, reference = min(candidates,
                    key=lambda item: abs((item[0] - base).to_value(u.s)))
                refmap = _overlay_helper._helper.load_map(reference)
                if _overlay_helper._is_preview_image(refmap,reference):
                    _overlay_helper._helper._hv_tags(refmap.meta)  # canonical Helioviewer tags
                if kind == 'LASCO':
                    if need_LASCO_calib:
                        refmap = _overlay_helper._datamine.calib_LASCO(
                            [refmap], instr=key, outdir=Path(cache_dir) / 'base_reference' / key)[0]
                    else:
                        refmap = _overlay_helper._helper._prepare_lasco(refmap, key)
                elif kind == 'EIT':
                    _overlay_helper._require_eit_level1(refmap, reference)
                elif kind == 'AIA' and not _overlay_helper._is_preview_image(refmap,reference):
                    meta = refmap.meta
                    if not (float(meta.get('lvl_num', 0)) >= 1.5 or
                            any(v in str(meta.get('series', '')).lower() for v in ('lev1.5', 'lev15'))):
                        raise ValueError('Base reference requires prepared AIA Level-1.5 FITS')
                # Calibration can substitute a nearby exposure: enforce the rule again.
                if not lo <= refmap.date <= hi:
                    raise ValueError(f'{name}: calibrated reference falls outside the allowed time window')
                reference_time = refmap.date
                refdir = Path(cache_dir).expanduser().resolve() / 'base_reference' / name.replace('/', '_')
                refdir.mkdir(parents=True, exist_ok=True)
                refpath = refdir / (reference_time.strftime('%Y%m%dT%H%M%S%f') + '.fits')
                refmap.save(refpath, overwrite=True)
                pairs = sorted(dict(zip(paths, times)).items(), key=lambda item: item[1])
                pairs = [(p, t) for p, t in pairs if t != reference_time]
                pairs.append((refpath, reference_time))
                pairs.sort(key=lambda item: item[1])
                paths, times = [p for p, _ in pairs], [t for _, t in pairs]
                print(f'[overlay] {name} {typ} reference: {reference_time.isot}')
        if d and len(paths) < 2 and not (typ == 'running' and first_frame_intensity):
            raise ValueError(f"{name}: difference images require at least two distinct observations.")
        layers.append(_overlay_helper._Layer(name, kind, key, times, paths, loader, z, d, typ if d else "", reference_time, t_start, t_end))

    lines = None
    if want_pfss:
        print("[overlay] building PFSS model ...")
        lines = (_overlay_helper._load_pfss(pfss_dir) if pfss_dir is not None else
                 build_pfss(sT, eT, Path(cache_dir) / 'pfss', **(pfss_kwargs or {})))

    print("[overlay] rendering frames ...")
    out = _overlay_helper._render_frames(layers, frame_times, outdir, start=0,
                         pfss_lines=lines, alpha=alpha, npix=npix, dpi=dpi, figsize=figsize,
                         N_Bloop=N_Bloop, tmap=tmap, first_frame_intensity=first_frame_intensity)
    print(f"[overlay] wrote {len(out)} PNGs to {outdir}")
    if uncalibrated_lasco:
        print('#####################\n LASCO images are not usable for science needs as it is not calibrated\n ##############')
    return out




def pfss_context_maps(sT, mission='SDO', instrument=None, wave=None, *, basedir='./', img_folder=None,
                      pfss_dir='PFSS', PFSS_dT=6., gong_cache='solar_data/gong', rss=2.5, nrho=35,
                      r_seed=1.05, n_lon=36, n_lat=18, NOAA_number=None, SRS_folder='NOAA_SRS',
                      ar_half_lon=12., ar_half_lat=12., AR_seed_spacing_deg=1., trace_longARrooted_loops=True,
                      AR_r_seed=1.01, n_closed=30, n_AR_open=5, n_global_closed=20,
                      AR_closedB_col='lime', AR_openB_col='magenta', closed_field_color='white',
                      npix=1024, dpi=120, maps='both', max_image_dt=3600.):
    """Save a full-Sun and an AR PFSS map for time sT (UTC) in basedir; return their paths.

    mission: 'SDO' (also 'LASCO'; Earth view), 'STEREO-A' or 'STEREO-B' (that spacecraft's view).
    instrument/wave: 'AIA' (e.g. 193), 'HMI' ('magnetogram'/'continuum') or 'EUVI'
    (171/195/284/304) background image nearest sT, downloaded into img_folder if
    missing (SDO: the SDO root, e.g. 'SDO'; STEREO: the STEREO root, e.g. 'STEREO').
    The nearest finished FITS already in img_folder within max_image_dt seconds is
    used first; only if there is none is an image within +-10 min downloaded.
    instrument=None draws the PFSS lines on a blank disk with a Stonyhurst grid.
    An invalid instrument/wave pair raises ValueError; a missing image falls back
    to the blank disk with a warning. On a blank full-Sun disk the AR rectangle
    (+-ar_half_lon x +-ar_half_lat, rotated to sT) is outlined when NOAA_number is set.
    maps: 'both' (default), 'full' (full-Sun map only) or 'AR' (AR map only; needs NOAA_number).
    PFSS: the UTC-anchored snapshot nearest sT in pfss_dir (PFSS_dT hours apart),
    built with build_pfss if missing; global lines drawn as in the frame makers.
    With NOAA_number the full-Sun map also shows the AR lines (closed AR_closedB_col,
    open AR_openB_col) over the global ones. AR map: AR-rooted lines from the shared cache
    (pfss_projection.ar_rooted_lines; needs Moving_structure_frames on the
    Python path), framed on the NOAA SRS position rotated to sT.
    Relative paths are inside basedir. Saves PFSS_fullSun_<label>_<time>.png
    and PFSS_AR<NOAA>_<label>_<time>.png in basedir.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from sunpy.coordinates import frames, get_earth, get_horizons_coord, propagate_with_solar_surface
    import solar_overlay_helper as oh
    import datamine_calib_helper as dh

    root = Path(basedir).expanduser().resolve()
    t = Time(sT)
    mission, inst = str(mission).upper(), (str(instrument).upper() if instrument else None)
    stereo = mission.startswith('STEREO')
    if stereo and mission not in ('STEREO-A', 'STEREO-B'):
        raise ValueError("mission must be 'SDO', 'LASCO', 'STEREO-A' or 'STEREO-B'")
    valid = {'AIA': ('94', '131', '171', '193', '211', '304', '335', '1600', '1700', '4500'),
             'HMI': ('magnetogram', 'continuum'), 'EUVI': ('171', '195', '284', '304'), None: (None,)}
    if inst not in valid or (inst and str(wave).lower() not in valid[inst]) or (inst == 'EUVI') != (stereo and inst is not None):
        raise ValueError(f'Invalid instrument/wave for mission {mission}: {instrument}/{wave}. AIA: {valid["AIA"]}; '
                         f'HMI: {valid["HMI"]} (mission SDO/LASCO); EUVI: {valid["EUVI"]} (mission STEREO-A/B); or instrument=None.')
    maps = str(maps).lower()
    if maps not in ('both', 'full', 'ar'):
        raise ValueError("maps must be 'both', 'full' or 'AR'")
    if maps == 'ar' and not NOAA_number:
        raise ValueError("maps='AR' needs NOAA_number")
    stamp = t.strftime('%Y%m%d_%H%M%S')

    # PFSS snapshot nearest sT (same folder naming and settings check as the frame makers)
    kw = dict(rss=rss, nrho=nrho, r_seed=r_seed, n_lon=n_lon, n_lat=n_lat)
    step = float(PFSS_dT) * 3600
    epoch = Time(np.round(t.unix / step) * step, format='unix')
    folder = root / pfss_dir / epoch.strftime('%Y%m%d_%H%M%S')
    meta = folder / 'metadata.json'
    ok = meta.is_file() and all(json.loads(meta.read_text()).get(k) == v for k, v in kw.items())
    if not ok:
        print('Generating PFSS snapshot:', epoch.isot)
        build_pfss(epoch.isot, epoch.isot, folder, cache_dir=str(root / gong_cache), **kw)
    global_lines = oh._load_pfss(folder)

    # Background image: nearest cached FITS within max_image_dt, else download near sT (+-10 min)
    img = None
    if inst:
        lo, hi = (t - 600 * u.s).isot, (t + 600 * u.s).isot
        src = root / (img_folder or ('STEREO' if inst == 'EUVI' else 'SDO'))
        local = src / (f'STEREO_{mission[-1]}/EUVI/{int(wave)}' if inst == 'EUVI' else f'{inst}/{wave}')
        cached = []
        for f in local.rglob('*.fits'):
            m = re.search(r'(\d{8})[T_](\d{6})', f.name)
            if m and 'raw' not in f.relative_to(local).parts:
                dt_s = abs((Time.strptime(''.join(m.groups()), '%Y%m%d%H%M%S') - t).to_value(u.s))
                if dt_s <= max_image_dt:
                    cached.append((dt_s, f))
        if cached:
            f = min(cached)[1]
            img = dh.load_map(str(f))
            print(f'Using cached {inst} {wave}: {f.name}')
        else:
            try:
                import Multiwb_datamine as dm
                found = (dm.get_STEREO_imgs(lo, hi, 600, mission[-1], ['EUVI'], [[int(wave)]], outdir=str(src), nproc=1)[0]
                         if inst == 'EUVI' else dm.get_SDO_imgs(lo, hi, 600, [inst], [[wave]], outdir=str(src), nproc=1)[0])
                img = min(found, key=lambda m: abs((m.date - t).to_value(u.s))) if found else None
            except Exception as exc:
                warnings.warn(f'No {inst} {wave} image near {t.isot}: {exc}; blank disk used.')
    observer = img.observer_coordinate if img is not None else (
        get_horizons_coord(mission, t) if stereo else get_earth(t))
    rsun = np.degrees(np.arctan((R_sun / observer.radius).decompose().value)) * 3600  # arcsec
    label = ((mission + '_' if stereo else '') + f'{inst}{wave}' if img is not None else (mission if stereo else 'Earth') + '-view')

    def canvas(centre, half):
        """Blank helioprojective map from the observer, centred on centre, +-half arcsec."""
        ref = SkyCoord(centre[0] * u.arcsec, centre[1] * u.arcsec, obstime=t, observer=observer, frame=frames.Helioprojective)
        hdr = sunpy.map.make_fitswcs_header((npix, npix), ref, scale=[2 * half / npix] * 2 * u.arcsec / u.pix)
        return sunpy.map.Map(np.zeros((npix, npix)), hdr)

    def draw(cmap_, lines, colors, nloop, title, path, rect=False, ar=None):
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(projection=cmap_)
        ax.set_facecolor('black')
        if img is not None:
            b = oh._regrid(img, cmap_)
            if inst == 'HMI':
                cm, norm = 'gray', (Normalize(-500, 500, clip=True) if str(wave).lower() == 'magnetogram'
                                    else Normalize(*np.nanpercentile(b, [1, 99.5]), clip=True))
            else:
                cm, norm = dh._image_display_settings(img, inst, wave)
            ax.imshow(b, cmap=cm, norm=norm, origin='lower')
        else:
            ax.imshow(np.full((npix, npix), np.nan), origin='lower')
        cmap_.draw_limb(axes=ax, color='w', lw=.8)
        cmap_.draw_grid(axes=ax, color='w', lw=.4, alpha=.5)
        oh._draw_pfss(ax, cmap_, lines, N_Bloop=nloop, colors=colors)
        if ar is not None:  # AR lines over the global field, frame-maker colours
            oh._draw_pfss(ax, cmap_, ar, N_Bloop=np.nan, colors=ar_colors)
        if rect:  # AR rectangle on the visible disk (Stonyhurst edges at the NOAA epoch, rotated to sT)
            e = np.linspace(-1, 1, 60)
            L = lon + np.r_[e * ar_half_lon, np.full(60, ar_half_lon), -e * ar_half_lon, np.full(60, -ar_half_lon)]
            B = lat + np.r_[np.full(60, -ar_half_lat), e * ar_half_lat, np.full(60, ar_half_lat), -e * ar_half_lat]
            c = SkyCoord(L * u.deg, B * u.deg, frame=frames.HeliographicStonyhurst, obstime=Time(ar_time))
            with propagate_with_solar_surface():
                z = c.transform_to(frames.Heliocentric(observer=observer, obstime=t)).z.value
                x, y = cmap_.wcs.world_to_pixel(c.transform_to(cmap_.coordinate_frame))
            x[z < 0] = np.nan
            ax.plot(x, y, color='yellow', lw=1.2)
            ax.text(np.nanmax(x), np.nanmax(y), f' NOAA {NOAA_number}', color='yellow', fontsize=9)
        ax.set_xlim(0, npix - 1); ax.set_ylim(0, npix - 1)
        ax.coords[0].set_axislabel('Solar X (arcsec)'); ax.coords[1].set_axislabel('Solar Y (arcsec)')
        ax.set_title(title, fontsize=10)
        fig.savefig(path, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        return str(path)

    when = t.strftime('%Y-%m-%d %H:%M UT')
    if NOAA_number:
        try:
            from noaa_region import find_NOAA_location
            from pfss_projection import ar_rooted_lines
        except ImportError as exc:
            raise ImportError('The AR map needs Moving_structure_frames (noaa_region, pfss_projection) on the Python path.') from exc
        lon, lat, ar_time, _ = find_NOAA_location(NOAA_number, t.isot, str(root / SRS_folder), local_file=str(root / SRS_folder))
        sel = ar_rooted_lines(folder, lon, lat, ar_time, ar_half_lon, ar_half_lat, AR_seed_spacing_deg, AR_r_seed,
                              NOAA_number, bool(trace_longARrooted_loops), n_closed, n_AR_open)
    ar_colors = {0: AR_closedB_col, 1: AR_openB_col, -1: AR_openB_col}
    out = []
    if maps in ('both', 'full'):
        out.append(draw(canvas((0, 0), 2.6 * rsun), global_lines, {1: '#ff4d4d', -1: '#4da6ff', 0: closed_field_color},
                        n_global_closed, f'PFSS {epoch.strftime("%Y-%m-%d %H:%M")} on {label} {when}',
                        root / f'PFSS_fullSun_{label}_{stamp}.png', rect=bool(NOAA_number) and img is None,
                        ar=sel if NOAA_number else None))
    if NOAA_number and maps in ('both', 'ar'):
        c = SkyCoord(lon * u.deg, lat * u.deg, frame=frames.HeliographicStonyhurst, obstime=Time(ar_time))
        with propagate_with_solar_surface():
            c = c.transform_to(frames.Helioprojective(observer=observer, obstime=t))
        half = 2 * rsun * np.sin(np.radians(max(ar_half_lon, ar_half_lat))) + .15 * rsun
        out.append(draw(canvas((c.Tx.to_value(u.arcsec), c.Ty.to_value(u.arcsec)), half), sel,
                        ar_colors, np.nan,
                        f'NOAA {NOAA_number} PFSS {epoch.strftime("%Y-%m-%d %H:%M")} on {label} {when}',
                        root / f'PFSS_AR{NOAA_number}_{label}_{stamp}.png'))
    print('Saved PFSS context maps:', *out, sep='\n  ')
    return out
