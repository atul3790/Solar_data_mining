"""Internal retrieval/calibration helpers for Multiwb_datamine and solar_overlay.

Keep this file beside Multiwb_datamine.py and solar_overlay.py. Users should
call the public functions in those modules; this is an implementation module.
Multiwb_datamine imports this module locally inside the APIs that need it.
No downloads, calibration runs or directory changes occur on import.
"""
import copy
import re
import shutil
import urllib.request
import warnings
from pathlib import Path
from contextlib import contextmanager
from functools import partial
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
import sunpy.map
from sunpy.map import all_coordinates_from_map
from sunpy.net import Fido, attrs as a
from sunpy.util.exceptions import SunpyUserWarning

AIA_WAVES = (94, 131, 171, 193, 211, 304, 335, 1600, 1700)


LASCO_FOV = {"C2": (2.4, 6.0), "C3": (3.8, 30.0)}


LASCO_STD = {"C2": ("orange", "clear"), "C3": ("clear", "clear")}


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


def _hdr_time(path) -> Time:
    """Observation time from the first HDU that has DATE-OBS."""
    with fits.open(path) as hl:
        for hdu in hl:
            h = hdu.header
            if "DATE-OBS" in h or "DATE_OBS" in h:
                d = str(h.get("DATE-OBS", h.get("DATE_OBS")))
                if "T" not in d and "TIME-OBS" in h:
                    d = d + "T" + str(h["TIME-OBS"])
                d = d.replace("/", "-").rstrip("Z")
                return Time(d, format="isot", scale="utc")
    raise ValueError(f"No DATE-OBS in {path}")


def _search(t_lo, t_hi, attrs, sample=None):
    q = [a.Time(t_lo, t_hi)] + list(attrs)
    if sample:
        q.append(a.Sample(float(sample) * u.s))
    res = Fido.search(*q)
    tabs = [r for r in res if len(r)]
    if not tabs:
        raise RuntimeError("Fido/VSO returned no records for this query.")
    return max(tabs, key=len)


def _download(rows, cache: Path):
    cache.mkdir(parents=True, exist_ok=True)
    if "Archive URL" in rows.colnames:
        files = []
        for url in rows["Archive URL"]:
            out = cache / str(url).rsplit("/", 1)[-1]
            temp = out.with_suffix(out.suffix + ".part")
            try:
                if not out.exists():
                    with urllib.request.urlopen(str(url), timeout=60) as src, temp.open("wb") as dst:
                        shutil.copyfileobj(src, dst)
                    temp.replace(out)
                files.append(str(out))
            except OSError as exc:
                temp.unlink(missing_ok=True)
                warnings.warn(f"Skipping failed archive download {url}: {exc}")
        return files
    files = Fido.fetch(rows, path=str(cache / "{file}"))
    if getattr(files, "errors", None):
        warnings.warn(f"{len(files.errors)} download(s) failed; they will be skipped.")
    return [str(f) for f in files]


def _lasco_archive_search(t_lo, t_hi, instr):
	"""Fallback LASCO level-0.5 index from the instrument team's daily archive."""
	rows = []
	days = np.arange(np.datetime64(t_lo.utc.isot[:10]),
		np.datetime64(t_hi.utc.isot[:10]) + np.timedelta64(1, "D"))
	for day in days:
		base = f"https://lasco-www.nrl.navy.mil/lz/level_05/{Time(str(day)).strftime('%y%m%d')}/{instr.lower()}/"
		try:
			with urllib.request.urlopen(base + "img_hdr.txt", timeout=60) as src:
				text = src.read().decode()
		except OSError as exc:
			warnings.warn(f"Could not read LASCO archive index for {day}: {exc}")
			continue
		for line in text.splitlines():
			fields = line.split()
			if len(fields) < 4 or not re.fullmatch(r"[23]\d+\.fts", fields[0]):
				continue
			time = Time(fields[1].replace("/", "-") + "T" + fields[2])
			if fields[3].upper() == instr and t_lo <= time <= t_hi:
				rows.append((time, base + fields[0]))
	if not rows:
		raise RuntimeError(f"No LASCO/{instr} records in VSO or the instrument archive.")
	return Table({"Start Time": Time([r[0] for r in rows]), "Archive URL": [r[1] for r in rows]})

def _radio_quicklook_fallback(sT,eT,outdir,downdata=True):
	'''
	Get candidate repository quicklooks for the supplied UTC interval.

	Returns saved image paths; with downdata=False, returns candidate URLs.
	Uses the supplied legacy archive templates and date-label conventions.
	Daily/panel images can extend beyond sT/eT and are not cropped.
	These previews are not numerical or flux-calibrated spectral arrays.
	'''
	from datetime import datetime, timedelta
	from pathlib import Path
	from urllib.request import urlopen
	import warnings

	start=datetime.strptime(sT,'%Y-%m-%dT%H:%M:%S')
	end=datetime.strptime(eT,'%Y-%m-%dT%H:%M:%S')
	if end<=start:
		raise ValueError('eT must be later than sT.')

	root=Path(outdir)/'repository_quicklooks'
	base=('https://www.sws.bom.gov.au/Category/World%20Data%20Centre/'
		'Data%20Display%20and%20Download/Spectrograph/station/')
	jobs={}
	day=start.replace(hour=0,minute=0,second=0)

	while day<end:
		next_day=day+timedelta(days=1)
		lo=max(start,day)
		hi=min(end,next_day)
		tag=day.strftime('%Y%m%d')

		#### Learmonth/Culgoora: legacy Australian date labels ####
		labels=set()
		if lo<day+timedelta(hours=10):
			labels.add(day)
		if hi>day+timedelta(hours=10):
			labels.add(next_day)
		for code,station in (('LEA','learmonth'),('CUL','culgoora')):
			for label in sorted(labels):
				d=label.strftime('%Y%m%d')
				jobs[f'{code}_{d}.gif']=(
					base+station+'/images/'+d[2:4]+'/'+
					d+'spectrograph.gif')

		#### Daily quicklooks ####
		jobs[f'NDA_{tag}.png']=(
			'https://rsdb.obs-nancay.fr/QL/Nda/NdaP1F'+tag+'.png')
		jobs[f'ORF_{tag}.png']=(
			'https://rsdb.obs-nancay.fr/QL/Orfees/orfees'+tag+'.png')
		jobs[f'ART_{tag}.png']=(
			'http://artemis-iv.phys.uoa.gr/QuickLooks/'+
			day.strftime('%Y/%m/')+tag+'_60s.Int01.png')

		#### HI-RAS panels; next-day a panel covers late UTC ####
		panels=[(day,'a'),(day,'b')]
		if hi>day+timedelta(hours=20):
			panels.append((next_day,'a'))
		for label,panel in panels:
			d=label.strftime('%y%m%d')
			jobs[f'HIR_{d}_{panel}.gif']=(
				'https://sunbase.nict.go.jp/solar/denpa/hirasDB/'+
				label.strftime('%Y/%m/')+d+panel+'.gif')

		#### Bleien: all overlapping 15-minute slots, L1 and L2 ####
		slot=lo.replace(minute=(lo.minute//15)*15,second=0,microsecond=0)
		while slot<hi:
			for level in ('L1','L2'):
				name=f'BLEN5M_{slot:%Y%m%d_%H%M%S}_{level}.fit.gz.png'
				jobs[name]=(
					'http://soleil.i4ds.ch/solarradio/qkl/'+
					slot.strftime('%Y/%m/%d/')+name)
			slot+=timedelta(minutes=15)
		day=next_day

	warnings.warn(
		'Falling back to PNG/GIF quicklooks from multiple observatories; '
		'these are previews, not numerical RSTN/e-CALLISTO spectra.')

	if not downdata:
		for url in jobs.values():
			print(url)
		return list(jobs.values())

	root.mkdir(parents=True,exist_ok=True)
	paths=[]
	for name,url in jobs.items():
		path=root/name
		try:
			if not path.is_file():
				with urlopen(url,timeout=30) as response:
					data=response.read()
				if not (
					data.startswith(b'\x89PNG\r\n\x1a\n') or
					data.startswith((b'GIF87a',b'GIF89a'))
				):
					raise ValueError('Response is not a PNG/GIF image.')
				path.write_bytes(data)
			paths.append(str(path))
			print('Quicklook:',path)
		except Exception as exc:
			print(f'Quicklook unavailable: {url}: {exc}')
	return paths
	
def _fetch_nearest(tab, targets: Time, cache: Path, validate=None, max_rounds=4):
    """
    For every target time download the nearest record (de-duplicated). If a
    downloaded file fails `validate(path)` the next-nearest record is tried.
    Returns a time-sorted list of (Time, path).
    """
    st = tab["Start Time"].mjd * 86400.0
    tg = targets.mjd * 86400.0
    order = np.argsort(np.abs(st[None, :] - tg[:, None]), axis=1)
    bad, done = set(), {}

    def _chosen():
        out = []
        for i in range(len(tg)):
            for r in order[i]:
                if int(r) not in bad:
                    out.append(int(r))
                    break
        return out

    for _ in range(max_rounds):
        need = sorted(set(_chosen()) - set(done))
        if not need:
            break
        files = _download(tab[need], cache)
        for f in files:
            try:
                tf = _hdr_time(f)
                r = need[int(np.argmin(np.abs(st[need] - tf.mjd * 86400.0)))]
                if validate is None or validate(f):
                    done[r] = (tf, f)
                else:
                    bad.add(r)
            except (OSError, ValueError, KeyError) as exc:
                warnings.warn(f"Skipping unusable FITS file {f}: {exc}")
        bad |= {r for r in need if r not in done}

    got = {r: done[r] for r in set(_chosen()) if r in done}
    if not got:
        raise RuntimeError("No usable images were downloaded (check network, "
                           "or relax the LASCO filter check with strict_filters=False).")
    return sorted(got.values(), key=lambda x: x[0].mjd)


def _series(kind, key, sT, eT, dT, cache_dir, validate=None, files=None):
    """Download the distinct images nearest to sT, sT+dT, ... eT."""
    targets = _target_times(sT, eT, dT)
    if files is not None:
        # Reuse existing products, including science-ready find_coobs() outputs.
        got = []
        for path in files:
            tf = _hdr_time(path)
            if targets[0] <= tf <= _parse_time(eT) and (validate is None or validate(path)):
                got.append((tf, str(path)))
        if not got:
            raise RuntimeError("No usable local images inside the requested interval.")
        st = np.array([t.mjd for t, _ in got])
        chosen = {int(np.argmin(np.abs(st - t.mjd))) for t in targets}
        return tuple(map(list, zip(*sorted((got[i] for i in chosen), key=lambda x: x[0].mjd))))
    pad = 1800.0 if kind == "LASCO" else 120.0
    t_lo = targets[0] - pad * u.s
    t_hi = _parse_time(eT) + pad * u.s
    sample = None  # Select cadence locally; VSO Sample can discard nearer endpoint images.
    if kind == "LASCO":
        attrs = [a.Source("SOHO"), a.Instrument("LASCO"), a.Detector(key)]
    elif kind == "AIA":
        attrs = [a.Source("SDO"), a.Instrument("AIA"), a.Wavelength(key * u.angstrom)]
    elif kind == "HMI":
        attrs = [a.Source("SDO"), a.Instrument("HMI"), a.Physobs("LOS_magnetic_field")]
    else:
        raise ValueError(kind)
    try:
        tab = _search(t_lo, t_hi, attrs, sample)
    except (RuntimeError, OSError) as exc:
        if kind != "LASCO":
            raise
        warnings.warn(f"VSO LASCO query failed ({exc}); trying the instrument archive.")
        tab = _lasco_archive_search(t_lo, t_hi, key)
    st = tab["Start Time"]
    tab = tab[(st >= targets[0]) & (st <= _parse_time(eT))]
    if not len(tab):
        raise RuntimeError(f"No {kind} images inside {sT} -- {eT}.")
    got = _fetch_nearest(tab, targets, Path(cache_dir) / f"{kind.lower()}{key or ''}", validate)
    return [g[0] for g in got], [g[1] for g in got]


def _lasco_ok(path, instr, strict=True):
    h = fits.getheader(path, 0)
    if str(h.get("DETECTOR", "")).strip().upper() != instr:
        return False
    if h.get("NAXIS1") != 1024 or h.get("NAXIS2") != 1024:
        return False
    if h.get("OV_HV", False):
        return True  # Explicitly marked Helioviewer product; filter fields may be absent.
    if not strict:
        return True
    f_exp, p_exp = LASCO_STD[instr]
    return (str(h.get("FILTER", "")).strip().lower() == f_exp
            and str(h.get("POLAR", "")).strip().lower() == p_exp)


def _lasco_level1_ok(path, instr, strict=True):
	"""Accept scientific Level-1 MSB FITS, never JP2-derived display counts."""
	if not _lasco_ok(path, instr, strict=strict):
		return False
	header = fits.getheader(path)
	if header.get('OV_HV', False):
		return False
	unit = str(header.get('BUNIT', '')).replace(' ', '').lower()
	try:
		level = float(header.get('LEVEL', header.get('LVL_NUM', 0)))
	except (TypeError, ValueError):
		return False
	# NRL's documented older Level-1 headers may omit LEVEL/LVL_NUM.
	processed = 'reduce_level_1' in str(header.get('HISTORY', '')).lower()
	return (level >= 1 or processed) and unit in ('msb', 'bsun', 'b_sun', 'meansolarbrightness')


def _prepare_lasco(m, instr=None):
    """
    LASCO quick-look preparation: exposure normalisation (DN/s) and masking of
    the occulter / outside-FOV pixels.

    NOTE: this is NOT the full SolarSoft reduce_level_1 (no vignetting, stray-light,
    distortion or photometric (MSB) calibration). Use calib_LASCO() for full
    calibration; this helper performs no full radiometric corrections.
    """
    instr = (instr or str(m.meta.get("detector", ""))).upper()
    lo, hi = LASCO_FOV[instr]
    data = m.data.astype(np.float32)
    meta = copy.deepcopy(m.meta)
    meta.pop("blank", None)  # Integer-only FITS keyword is invalid for floating data.
    unit = str(m.meta.get("bunit") or m.unit).replace(" ", "").lower()
    calibrated = unit in ("msb", "bsun", "b_sun", "meansolarbrightness")
    if not calibrated and not meta.get("ov_norm", False) and unit not in ("dn/s", "ct/s"):
        exposure = float(m.meta.get("exptime", 0))
        if not np.isfinite(exposure) or exposure <= 0:
            raise ValueError("LASCO EXPTIME must be finite and positive.")
        data /= exposure
        meta["bunit"], meta["ov_norm"] = "DN/s", True
    if not calibrated:
        warnings.warn("LASCO preparation is exposure normalisation and FOV masking only; "
                      "full Level-1 calibration requires verified MSB products or SolarSoft REDUCE_LEVEL_1.")
    xy = all_coordinates_from_map(m)
    rho = np.hypot(xy.Tx.to_value(u.arcsec), xy.Ty.to_value(u.arcsec)) / m.rsun_obs.to_value(u.arcsec)
    data[(rho < lo) | (rho > hi)] = np.nan
    if not calibrated:
        meta["ov_cal"] = "exposure/FOV only; not full Level 1"
    with _quiet():
        return sunpy.map.Map(data, meta)


def _calib_lasco_backend(images, instr=None, outdir='./solar_data/lasco_level1/',
                sswidl_command='sswidl', timeout=600, overwrite=False, calib_func=None,
                fallback_level1=True, fallback_max_dt=180):
	"""Private SolarSoft backend; remote fallbacks belong to public wrappers."""
	import hashlib
	import json
	import os
	import shlex
	import subprocess
	import tempfile

	if isinstance(images, (list, tuple)):
		return [_calib_lasco_backend(x, instr, outdir, sswidl_command, timeout, overwrite, calib_func, fallback_level1, fallback_max_dt) for x in images]
	if not isinstance(fallback_level1, (bool, np.bool_)):
		raise ValueError('fallback_level1 must be True or False.')
	if not np.isfinite(fallback_max_dt) or fallback_max_dt <= 0:
		raise ValueError('fallback_max_dt must be positive seconds.')
	if not np.isfinite(timeout) or timeout <= 0:
		raise ValueError('timeout must be finite and positive.')
	path = Path(images).expanduser().resolve() if isinstance(images, (str, Path)) else None
	command = shlex.split(sswidl_command) if isinstance(sswidl_command, str) else list(sswidl_command)
	engine = 'callback:' + calib_func.__module__ + '.' + getattr(calib_func, '__qualname__', type(calib_func).__name__) if calib_func else 'SolarSoft:' + repr(command)
	fingerprint = hashlib.sha256((engine + repr((fallback_level1, fallback_max_dt))).encode())
	if path is not None:
		with path.open('rb') as src:
			for block in iter(lambda: src.read(1024 * 1024), b''):
				fingerprint.update(block)
	else:
		if not isinstance(images, sunpy.map.GenericMap):
			raise TypeError('images must be a FITS path, a SunPy Map, or a list of these.')
		fingerprint.update(np.ascontiguousarray(images.data).tobytes())
		fingerprint.update(repr(images.meta).encode())
	fingerprint.update(str(instr).encode())
	cache = Path(outdir).expanduser().resolve() / fingerprint.hexdigest()[:24]
	product = cache / 'calibrated.fits'
	if product.exists() and not overwrite:
		result = sunpy.map.Map(product)
		if not result.meta.get('ov_full', False):
			raise ValueError(f'Invalid calibration cache: {product}')
		return result
	m = sunpy.map.Map(path) if path is not None else images
	detector = str(m.meta.get('detector', '')).strip().upper()
	if detector not in LASCO_FOV or (instr is not None and instr.upper() != detector):
		raise ValueError('Input must be LASCO C2/C3 data matching instr.')
	unit = str(m.meta.get('bunit', '')).replace(' ', '').lower()
	if unit in ('msb', 'bsun', 'b_sun', 'meansolarbrightness') and float(m.meta.get('level', 0)) >= 1:
		return m  # Already calibrated; never apply the corrections twice.
	if m.meta.get('ov_norm', False) or unit in ('dn/s', 'ct/s'):
		raise ValueError('Full calibration requires the original raw LASCO DN image, not a DN/s preview.')
	cache.mkdir(parents=True, exist_ok=True)
	if calib_func is not None:
		result = calib_func(m)
	else:
		try:
			if not command or shutil.which(command[0]) is None:
				raise RuntimeError('Full LASCO calibration requires SolarSoft/IDL with LASCO and LASCO_DATA. '
					'Set sswidl_command to its launcher, or provide lasco_calib_func. '
					'Use need_LASCO_calib=False for quick looks.')
			log = cache / 'sswidl.log'
			with tempfile.TemporaryDirectory(prefix='lasco_', dir=cache) as tmp:
				work = Path(tmp)
				raw = work / 'raw.fts'
				if path is not None and path.suffix.lower() not in ('.gz', '.fz'):
					shutil.copy2(path, raw)
				else:
					m.save(raw, overwrite=True)
				output = work / 'level1'
				output.mkdir()
				# Quote IDL strings; subprocess uses argv with shell=False.
				qraw = str(raw).replace("'", "''")
				qout = str(output).replace("'", "''")
				batch = work / 'calibrate.pro'
				batch.write_text(f"reduce_level_1, '{qraw}', SAVEDIR='{qout}', /NOSCALE\nexit\n")
				env = os.environ.copy()
				env['SSW_INSTR'] = ' '.join(dict.fromkeys(env.get('SSW_INSTR', '').split() + ['lasco']))
				try:
					with log.open('w') as stream:
						proc = subprocess.run(command + [str(batch)], cwd=work, env=env,
							stdout=stream, stderr=subprocess.STDOUT, timeout=timeout, check=False)
				except subprocess.TimeoutExpired as exc:
					raise RuntimeError(f'LASCO calibration timed out; see {log}') from exc
				files = sorted(output.glob('*.fts')) + sorted(output.glob('*.fits'))
				if proc.returncode != 0 or len(files) != 1:
					raise RuntimeError(f'SolarSoft did not produce one calibrated FITS image; see {log}')
				result = sunpy.map.Map(files[0])
				if not _lasco_level1_ok(files[0], detector, strict=False):
					raise ValueError('SolarSoft output is not verified Level-1 MSB FITS.')
				result.meta['level'] = float(result.meta.get('level', result.meta.get('lvl_num', 1)))
		except (RuntimeError, OSError, ValueError):
			raise
	if not isinstance(result, sunpy.map.GenericMap):
		raise TypeError('The LASCO calibration backend must return a SunPy Map.')
	unit = str(result.meta.get('bunit', '')).replace(' ', '').lower()
	if (str(result.meta.get('detector', '')).strip().upper() != detector
		or unit not in ('msb', 'bsun', 'b_sun', 'meansolarbrightness')
		or float(result.meta.get('level', 0)) < 1):
		raise ValueError('Calibration output must have the correct detector, LEVEL>=1 and MSB brightness units.')
	meta = copy.deepcopy(result.meta)
	meta.pop('ov_cal', None)
	meta['ov_full'], meta['ov_calsr'] = True, engine
	result = sunpy.map.Map(result.data, meta)
	staged = cache / 'calibrated.part.fits'
	result.save(staged, overwrite=True)
	staged.replace(product)
	(cache / 'provenance.json').write_text(json.dumps(dict(source=str(path) if path else 'input Map',
		backend=engine, sha256=fingerprint.hexdigest(), product=str(product)), indent=2))
	return result


def _load_lasco(path, instr, calibrate=True, calib_func=None, need_LASCO_calib=False, calib_kwargs=None):
    if need_LASCO_calib:
        from Multiwb_datamine import calib_LASCO
        return calib_LASCO(path, instr=instr, **(calib_kwargs or {}))
    with _quiet():
        m = sunpy.map.Map(path)
    if calib_func is not None:
        return calib_func(m)
    return _prepare_lasco(m, instr) if calibrate else m


def _get_lasco_source(sT, eT, dT, instr, cache_dir, source, outdir=None,
                      files=None, calibrate=False, strict_filters=True, calib_func=None):
    """Internal single-source retrieval; public wrappers own fallback decisions."""
    instr = instr.upper()
    if instr not in LASCO_FOV:
        raise ValueError("instr must be 'C2' or 'C3'")
    if source not in ('level1', 'archive', 'helioviewer'):
        raise ValueError("source must be 'level1', 'archive' or 'helioviewer'")
    if source == 'level1':
        if calib_func is not None:
            raise ValueError('Level-1 FITS are already calibrated; use source="archive" for a raw-data calibration callback.')
        targets = _target_times(sT, eT, dT)
        folder = Path(cache_dir).expanduser().resolve() / ('lasco_level1_' + instr)
        validate = partial(_lasco_level1_ok, instr=instr, strict=strict_filters)
        if files is not None:
            _, paths = _series('LASCO', instr, sT, eT, dT, cache_dir, validate=validate, files=files)
        else:
            tab = _search(targets[0], _parse_time(eT),
                          [a.Source('SOHO'), a.Instrument('LASCO'), a.Detector(instr),
                           a.Provider('SDAC'), a.Level(1)])
            times = tab['Start Time']
            tab = tab[(times >= targets[0]) & (times <= _parse_time(eT))]
            # Paper Table 2: ~4106 KB L1 float FITS versus ~2100 KB L0.5.
            # Broad size sieve only: final FITS metadata must prove Level 1/MSB.
            '''
            if 'Size' in tab.colnames and getattr(tab['Size'], 'unit', None) is not None:
                sizes = u.Quantity(tab['Size']).to_value(u.Mibyte)
                tab = tab[(sizes >= 3.8) & (sizes <= 4.8)]
            '''    
            if not len(tab):
                raise RuntimeError('No LASCO Level-1 FITS candidates in VSO/SDAC for this interval; no raw or JP2 fallback.')
            got = _fetch_nearest(tab, targets, folder,
                                 validate=lambda p: validate(p) and targets[0] <= _hdr_time(p) <= _parse_time(eT))
            paths = [path for _, path in got]
        maps = []
        for path in paths:
            image = sunpy.map.Map(path)
            image.meta['level'] = float(image.meta.get('level', image.meta.get('lvl_num', 1)))
            image.meta['ov_full'] = True
            image.meta['ov_calsr'] = 'Downloaded scientific Level-1 FITS; VSO/SDAC route'
            image.meta.pop('ov_cal', None)
            # Already calibrated: preserve brightness array and WCS unchanged.
            maps.append(image)
        return _save_maps(maps, outdir, 'LASCO_' + instr) if outdir is not None else maps
    if source == 'helioviewer' and files is None:
        if calib_func is not None:
            raise ValueError('Helioviewer JP2 is a visualization route; use raw archive FITS for calibration callbacks.')
        try:
            import hvpy
            from hvpy.datasource import DataSource
        except ImportError as exc:
            raise ImportError('Helioviewer LASCO retrieval needs hvpy; JP2 reading needs glymur and OpenJPEG.') from exc
        targets = _target_times(sT, eT, dT)
        folder = Path(cache_dir).expanduser().resolve() / ('helioviewer_' + instr)
        folder.mkdir(parents=True, exist_ok=True)
        source_id = getattr(DataSource, 'LASCO_' + instr).value
        selected = {}
        for target in targets:
            path = folder / f"LASCO_{instr}_{target.strftime('%Y%m%dT%H%M%S')}.jp2"
            if not path.exists():
                payload = hvpy.getJP2Image(target.datetime, source_id)
                staged = path.with_suffix('.jp2.part')
                staged.write_bytes(payload)
                staged.replace(path)
            image = sunpy.map.Map(path)
            if not _parse_time(sT) <= image.date <= _parse_time(eT):
                continue
            if image.date.isot in selected:
                continue
            meta = copy.deepcopy(image.meta)
            meta['instrume'], meta['detector'] = 'LASCO', instr
            meta['bunit'] = 'ct'  # Display counts, not reconstructed MSB or raw DN.
            meta['ov_hv'], meta['ov_norm'] = True, True
            meta['ov_src'] = 'Helioviewer JP2'
            meta.pop('ov_full', None)
            prepared = sunpy.map.Map(image.data, meta)
            prepared = _prepare_lasco(prepared, instr) if calibrate else prepared
            prepared.meta['ov_cal'] = 'Helioviewer visualization; not photometric FITS'
            selected[image.date.isot] = prepared
        if not selected:
            raise RuntimeError('No Helioviewer LASCO observations inside the requested interval.')
        maps = [selected[t] for t in sorted(selected)]
        warnings.warn('Helioviewer LASCO JP2 is for visualization; quantitative brightness requires scientific calibrated FITS.')
        return _save_maps(maps, outdir, 'LASCO_' + instr) if outdir is not None else maps
    _, paths = _series("LASCO", instr, sT, eT, dT, cache_dir,
                       validate=partial(_lasco_ok, instr=instr, strict=strict_filters), files=files)
    maps = []
    for path in paths:
        image = sunpy.map.Map(path)
        image.meta.setdefault('ov_raw', str(Path(path).expanduser().resolve()))
        maps.append(_prepare_lasco(image, instr) if calibrate and not _lasco_science_ready(image) else image)
    return _save_maps(maps, outdir, "LASCO_" + instr) if outdir else maps


def _lasco_science_ready(image):
	"""Recognize verified Level-1/MSB products without altering their arrays."""
	unit = str(image.meta.get('bunit', '')).replace(' ', '').lower()
	try:
		level = float(image.meta.get('level', image.meta.get('lvl_num', 0)))
	except (TypeError, ValueError):
		level = 0
	processed = 'reduce_level_1' in str(image.meta.get('history', '')).lower()
	return (not image.meta.get('ov_hv', False) and
		unit in ('msb', 'bsun', 'b_sun', 'meansolarbrightness') and
		(level >= 1 or processed))


def _lasco_notice(images):
	if any(not _lasco_science_ready(m) for m in images):
		warnings.warn('LASCO output includes visualization products, not fully calibrated MSB photometry. '
			'Helioviewer JP2 is compressed display data. Level-0.5 DN/s lacks full '
			'vignetting, geometric-distortion and radiometric corrections; exposure normalization alone is insufficient.')


def _load_hmi(path):
    with _quiet():
        return sunpy.map.Map(path)


def _save_maps(maps, outdir, prefix):
	"""Write finished Maps with their complete WCS, retaining the original downloads."""
	folder = Path(outdir).expanduser().resolve()
	folder.mkdir(parents=True, exist_ok=True)
	for m in maps:
		m.save(folder / f"{prefix}_{m.date.strftime('%Y%m%dT%H%M%S%f')}.fits", overwrite=True)
	return maps


def _fits_paths(directory):
	root = Path(directory).expanduser().resolve()
	if root.is_file():
		return [root]
	if not root.is_dir():
		raise FileNotFoundError(root)
	return sorted(p for p in root.rglob('*') if p.is_file() and
		(p.name.lower().endswith(('.fits', '.fts', '.fit', '.fits.gz', '.fts.gz'))))


def _require_eit_level1(m, source):
	level = m.meta.get('data_lev', m.meta.get('level', m.meta.get('lvl_num', m.meta.get('lev_num', 0))))
	if str(level).strip().upper() not in ('L1', '1', '1.0', 'LEVEL 1'):
		raise ValueError(f'{source}: calibrated EIT Level-1 metadata is required')


def _matching_files(directories, kind, key):
	paths = []
	for directory in directories:
		for path in _fits_paths(directory):
			try:
				m = sunpy.map.Map(path)
				ins = str(m.instrument).upper()
				match = (kind in ins or kind in str(m.detector).upper())
				if kind == 'LASCO':
					match &= str(m.detector).upper() == key
				elif kind in ('AIA', 'EIT'):
					match &= int(m.wavelength.to_value(u.angstrom)) == key
				elif kind == 'HMI':
					match &= ('magnetogram' in str(m.measurement).lower() or
						'm_' in str(m.meta.get('content', '')).lower() or
						'magnetic' in str(m.meta.get('physobs', '')).lower())
				if match:
					paths.append(path)
			except (OSError, ValueError, KeyError, TypeError):
				continue
	if not paths:
		raise FileNotFoundError(f'No {kind}/{key} FITS maps in supplied directories')
	return sorted(set(paths))


def _lasco_exposure_id(image):
    """Native exposure identifier from a LASCO archive filename, when present."""
    name = str(image.meta.get('filename', '')).rsplit('/', 1)[-1]
    match = re.fullmatch(r'[23]\d(\d{6})\.(?:fts|fits)', name, re.IGNORECASE)
    return match[1] if match else None
