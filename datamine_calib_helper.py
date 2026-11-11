"""Internal retrieval/calibration helpers for Multiwb_datamine and solar_overlay.

Keep this file beside Multiwb_datamine.py and solar_overlay.py. Users should
call the public functions in those modules; this is an implementation module.
Multiwb_datamine imports this module locally inside the APIs that need it.
No downloads, calibration runs or directory changes occur on import.
"""
import copy
import re
import shutil
import time
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
from astropy.coordinates import SkyCoord
from astropy.visualization import ImageNormalize, AsinhStretch, PowerStretch
from matplotlib.colors import Normalize
from sunpy.visualization import colormaps as _sunpy_cmaps
import sunpy.map
from sunpy.map import all_coordinates_from_map
from sunpy.net import Fido, attrs as a
from sunpy.util.exceptions import SunpyUserWarning

AIA_WAVES = (94, 131, 171, 193, 211, 304, 335, 1600, 1700)


LASCO_FOV = {"C2": (2.4, 6.0), "C3": (3.8, 30.0)}


LASCO_STD = {"C2": ("orange", "clear"), "C3": ("clear", "clear")}


def _image_display_settings(m,instrument,product=None,preview=False):
	"""Instrument colours and native Map stretch; limits from source pixels.

	EUV FITS use a copy of Map.plot_settings['norm'] with 1--99.9 percentile
	limits to retain faint disk emission and bright flare cores. Display JP2
	pixels stay linear (already stretched). Input data/metadata are unchanged.
	"""
	kind=str(instrument).upper()
	data=np.asarray(m.data,float)
	valid=np.isfinite(data)
	if m.mask is not None:valid&=~np.asarray(m.mask,bool)
	values=data[valid]
	if not values.size:raise ValueError('No finite unmasked image pixels.')
	preview=preview or bool(m.meta.get('ov_hv',False))
	if kind in ('AIA','EUVI','EIT'):
		wave=int(float(product)) if product is not None else int(round(m.wavelength.to_value(u.angstrom)))
		cmap=('sdoaia' if kind=='AIA' else 'euvi' if kind=='EUVI' else 'sohoeit')+str(wave)
	else:
		cmap=m.plot_settings.get('cmap','gray')
	if preview:
		lo,hi=(0.,255.) if values.min()>=0 and values.max()<=255 else (float(values.min()),float(values.max()))
		norm=Normalize(vmin=lo,vmax=max(hi,lo+np.finfo(float).eps),clip=True)
	else:
		norm=copy.deepcopy(m.plot_settings.get('norm'))
		if norm is None:
			norm=ImageNormalize(stretch=PowerStretch(.25) if kind=='EUVI' else AsinhStretch(.01))
		lo,hi=np.percentile(values,[1,99.9])
		norm.vmin=float(lo)
		norm.vmax=float(hi) if hi>lo else float(lo)+1.
		norm.clip=True
	return cmap,norm


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
    if Path(path).suffix.lower() in (".jp2", ".j2k"):
        return sunpy.map.Map(path).date
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


DOWNLOAD_RETRIES = 3       # Fido/VSO attempts per record
DOWNLOAD_BACKOFF_S = 10.   # wait grows 10 s, 20 s between attempts


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
    # VSO/SSC answer intermittently with 504s or dropped connections; retry
    # with back-off before giving up on these records (DOWNLOAD_RETRIES tries).
    last = None
    for attempt in range(DOWNLOAD_RETRIES):
        try:
            files = Fido.fetch(rows, path=str(cache / "{file}"))
            if getattr(files, "errors", None) and not len(files):
                raise RuntimeError(f"{len(files.errors)} download error(s): {files.errors[0]}")
            if getattr(files, "errors", None):
                warnings.warn(f"{len(files.errors)} download(s) failed; they will be skipped.")
            return [str(f) for f in files]
        except Exception as exc:
            last = exc
            if attempt + 1 < DOWNLOAD_RETRIES:
                wait = DOWNLOAD_BACKOFF_S * (attempt + 1)
                print(f"Download attempt {attempt + 1}/{DOWNLOAD_RETRIES} failed ({exc}); retrying in {wait:.0f} s")
                time.sleep(wait)
    raise last


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


def _lasco_obstime(m):
	"""Combine split LASCO DATE_OBS/TIME_OBS for WCS, without editing input data.

	SunPy LASCOMap.date combines them, but reference_date can otherwise use
	midnight. Leave other instruments and already complete timestamps alone.
	"""
	if 'LASCO' not in (str(m.meta.get('instrume',''))+' '+str(m.meta.get('detector',''))).upper():
		return m
	date=str(m.meta.get('date-obs',m.meta.get('date_obs',''))).strip()
	if len(date)!=10 or not m.meta.get('time-obs',m.meta.get('time_obs')):
		return m
	meta=copy.deepcopy(m.meta)
	meta['date-obs']=m.date.utc.isot
	if 'date_obs' in meta:meta['date_obs']=meta['date-obs']
	fixed=sunpy.map.Map(m.data,meta,mask=m.mask)
	fixed.plot_settings=copy.deepcopy(m.plot_settings)
	return fixed


def _prepare_lasco(m, instr=None):
    """
    LASCO quick-look preparation: exposure normalisation (DN/s) and masking of
    the occulter / outside-FOV pixels.

    NOTE: this is NOT the full SolarSoft reduce_level_1 (no vignetting, stray-light,
    distortion or photometric (MSB) calibration). Use calib_LASCO() for full
    calibration; this helper performs no full radiometric corrections.
    """
    m = _lasco_obstime(m)
    if _lasco_science_ready(m):
        return m  # Verified Level-1 MSB: geometry and brightness are final.
    instr = (instr or str(m.meta.get("detector", ""))).upper()
    lo, hi = LASCO_FOV[instr]
    data = m.data.astype(np.float32)
    meta = copy.deepcopy(m.meta)
    meta.pop("blank", None)  # Integer-only FITS keyword is invalid for floating data.
    _add_soho_observer(meta, m.date)  # Saved FITS carry SOHO's real position.
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
    if not calibrated and not meta.get("ov_hv", False):  # Helioviewer keeps its canonical tags.
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
        folder = Path(cache_dir).expanduser().resolve() / ('helioviewer_' + instr)
        selected = {}
        # Shared fetcher: one JP2 per cadence bin, cache reused, one inward retry at edges.
        for path in _fill_helioviewer_jp2(sT, eT, dT, [], folder, 'LASCO', instr):
            image = sunpy.map.Map(path)
            if not _parse_time(sT) <= image.date <= _parse_time(eT):
                continue
            if image.date.isot in selected:
                continue
            meta = _hv_tags(copy.deepcopy(image.meta), path)  # Canonical tags; roll 0.
            meta['instrume'], meta['detector'] = 'LASCO', instr
            _add_soho_observer(meta, image.date)
            prepared = sunpy.map.Map(image.data, meta)
            prepared = _prepare_lasco(prepared, instr) if calibrate else prepared
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


def jp2_to_fits(jp2_paths, outdir, skip_times=(), tolerance=1.0):
	"""Convert Helioviewer JP2s to FITS (maps_HV); return the FITS paths.

	Pixels are kept as Helioviewer display values (0-255, BUNIT='ct'); the
	embedded header gives the WCS, with roll 0 because Helioviewer pixels are
	north-up (old value in OV_ROLL). Spacecraft observer keys are retained.
	LASCO also gets the occulter/FOV mask. Tags: OV_HV, OV_SRC, OV_CAL.
	skip_times: observation Times already held as Helioviewer FITS elsewhere.
	Existing outputs are reused. Display products: morphology and motion only.
	"""
	folder = Path(outdir).expanduser().resolve()
	skip = [_parse_time(t) for t in skip_times]
	out = []
	for path in jp2_paths:
		try:
			m = load_map(path)
			if any(abs((m.date - t).to_value(u.s)) <= tolerance for t in skip):
				continue
			meta = copy.deepcopy(m.meta)
			meta.pop('helioviewer', None)  # JP2-only XML block; not FITS-safe.
			det = str(meta.get('detector', '')).strip().upper()
			ins = str(meta.get('instrume', '')).strip().upper()
			obs = re.sub(r'[^A-Z0-9]', '', str(meta.get('obsrvtry', meta.get('telescop', ''))).upper())
			if 'LASCO' in ins + det or det in LASCO_FOV:
				meta['instrume'], meta['detector'] = 'LASCO', det
				name = f'LASCO_{det}'
			elif det in ('EUVI', 'COR1', 'COR2'):
				name = f"STEREO_{'A' if obs.endswith(('A', 'AHEAD')) else 'B'}_{det}"
				if det == 'EUVI':
					name += f'_{int(round(m.wavelength.to_value(u.angstrom)))}'
			elif 'HMI' in ins:
				name = 'HMI_' + ('continuum' if 'continuum' in str(meta.get('content', '')).lower() else 'magnetogram')
			elif 'AIA' in ins:
				name = f'AIA_{int(round(m.wavelength.to_value(u.angstrom)))}'
			else:
				name = re.sub(r'[^A-Za-z0-9]', '', ins or 'IMG')
			_hv_tags(meta, path)  # Same canonical tags as every other JP2 -> FITS route.
			_add_rsun_ref(meta)
			target = folder / f"{name}_HV_{m.date.strftime('%Y%m%dT%H%M%S%f')}.fits"
			if not target.exists():
				folder.mkdir(parents=True, exist_ok=True)  # Only when something is written.
				with _quiet():
					hv = sunpy.map.Map(m.data, meta)
					if name.startswith('LASCO'):
						hv = _prepare_lasco(hv, det)
					hv.save(target)
			out.append(target)
		except Exception as exc:
			warnings.warn(f'Cannot convert {path} to FITS: {exc}')
	return sorted(out)


def _load_hmi(path):
    return load_map(path)


def _hv_north_up(meta):
	"""Helioviewer JP2 pixels are already rotated to solar north.

	Their headers keep the original CROTA, which SunPy ignores only while the
	JP2-only 'helioviewer' key exists; FITS saving drops that key. Zero the
	roll in the header (original kept in OV_ROLL) so the WCS stays truthful.
	"""
	roll = next((meta[k] for k in ('crota2', 'crota', 'crota1') if k in meta), None)
	if roll is not None and 'ov_roll' not in meta:
		meta['ov_roll'] = float(roll)
	for k in ('crota', 'crota1', 'crota2'):
		if k in meta:
			meta[k] = 0.0
	for k, v in (('pc1_1', 1.), ('pc1_2', 0.), ('pc2_1', 0.), ('pc2_2', 1.)):
		if k in meta:
			meta[k] = v
	return meta


HV_TAGS = {'bunit': 'ct', 'ov_hv': True, 'ov_norm': True, 'ov_src': 'Helioviewer JP2',
	'ov_cal': 'Helioviewer visualization; not photometric FITS'}


def _hv_tags(meta, jp2_path=None):
	"""Give a Helioviewer-derived image the one canonical identity/unit tag set.

	Every JP2 -> FITS conversion (module, overlay and flow code) uses this, so
	they share exactly the same OV_HV, OV_SRC, OV_CAL, OV_NORM and BUNIT and
	form one product family. Roll is set to 0 (old value in OV_ROLL).
	"""
	meta.pop('helioviewer', None)  # JP2-only XML block; not FITS-safe.
	for key in ('ov_full', 'ov_calsr'):
		meta.pop(key, None)
	meta.update(HV_TAGS)
	if jp2_path is not None:
		meta['ov_jp2'] = Path(jp2_path).name
	return _hv_north_up(meta)


SOHO_ORBIT_URL = 'https://soho.nascom.nasa.gov/data/ancillary/orbit/predictive/{year}/'
SOHO_ORBIT_CACHE = Path.home() / '.cache' / 'solar_data_mining' / 'soho_orbit'
L1_DISTANCE_KM = 1.5e6
_soho_days = {}


def _soho_orbit_day(day):
	"""(Time array, GSE x/y/z km) of SOHO for one UTC day from the SOHO orbit file.

	Uses the highest-version SO_OR_PRE_YYYYMMDD_Vnn.FITS (10-min samples),
	cached in ~/.cache/solar_data_mining/soho_orbit. Returns None if unavailable.
	"""
	tag = day.strftime('%Y%m%d')
	if tag in _soho_days:
		return _soho_days[tag]
	result = None
	try:
		SOHO_ORBIT_CACHE.mkdir(parents=True, exist_ok=True)
		cached = sorted(SOHO_ORBIT_CACHE.glob(f'SO_OR_PRE_{tag}_V*.FITS'))
		if not cached:
			url = SOHO_ORBIT_URL.format(year=day.strftime('%Y'))
			with urllib.request.urlopen(url, timeout=30) as src:
				names = sorted(set(re.findall(rf'SO_OR_PRE_{tag}_V\d+\.FITS', src.read().decode('utf-8', 'replace'))))
			if names:
				target = SOHO_ORBIT_CACHE / names[-1]
				with urllib.request.urlopen(url + names[-1], timeout=60) as src:
					target.write_bytes(src.read())
				cached = [target]
		if cached:
			with fits.open(cached[-1]) as hl:
				d = hl[1].data
				t = Time([f'{y:04d}:{doy:03d}:00:00:00' for y, doy in zip(d['Year'], d['Day of year'])],
					format='yday', scale='utc') + np.asarray(d['Ellapsed milliseconds of day'], float) * u.ms
				xyz = np.vstack([d['GSE x (km)'], d['GSE y (km)'], d['GSE z (km)']]).astype(float)
			result = (t, xyz)
	except Exception as exc:
		warnings.warn(f'SOHO orbit file for {tag} unavailable: {exc}')
	_soho_days[tag] = result
	return result


def soho_observer(time):
	"""SOHO position as a HeliographicStonyhurst SkyCoord, and its source text.

	From the SOHO orbit file (GSE offset from Earth, interpolated in time);
	if that is unavailable, the L1 approximation: Earth's direction, 1.5e6 km
	sunward. SOHO's own roll/pointing are not involved: they are in CROTA/CRVAL.
	"""
	from sunpy.coordinates import frames as sframes
	time = _parse_time(time)
	samples = [_soho_orbit_day(time + k * u.day) for k in (-1, 0, 1)]
	samples = [s for s in samples if s is not None]
	if samples:
		t = Time(np.concatenate([s[0].jd for s in samples]), format='jd', scale='utc')
		xyz = np.hstack([s[1] for s in samples])
		order = np.argsort(t.jd)
		tj, xyz = t.jd[order], xyz[:, order]
		if tj[0] <= time.jd <= tj[-1]:
			gse = [np.interp(time.jd, tj, axis) for axis in xyz]
			c = SkyCoord(*(gse * u.km), representation_type='cartesian',
				frame=sframes.GeocentricSolarEcliptic(obstime=time))
			return c.transform_to(sframes.HeliographicStonyhurst(obstime=time)), 'SOHO orbit file (GSE)'
	c = SkyCoord(L1_DISTANCE_KM * u.km, 0 * u.km, 0 * u.km, representation_type='cartesian',
		frame=sframes.GeocentricSolarEcliptic(obstime=time))
	return c.transform_to(sframes.HeliographicStonyhurst(obstime=time)), 'L1 approximation (no SOHO orbit file)'


def _is_soho(meta):
	label = ' '.join(str(meta.get(k, '')) for k in ('telescop', 'obsrvtry', 'instrume', 'detector')).upper()
	return 'SOHO' in label or 'LASCO' in label or 'EIT' in label


def _add_soho_observer(meta, time):
	"""Add SOHO HGLN/HGLT/DSUN_OBS and RSUN_REF when a SOHO header lacks them."""
	if not _is_soho(meta) or all(k in meta for k in ('hgln_obs', 'hglt_obs', 'dsun_obs')):
		return meta
	obs, source = soho_observer(time)
	meta.update({'hgln_obs': float(obs.lon.to_value(u.deg)), 'hglt_obs': float(obs.lat.to_value(u.deg)),
		'dsun_obs': float(obs.radius.to_value(u.m)), 'ov_obs': source})
	meta.setdefault('rsun_ref', 695700000.0)  # SunPy's standard photospheric radius, made explicit.
	return meta


def _is_helioviewer(meta):
	"""Helioviewer JP2 or a FITS/difference derived from one."""
	return (bool(meta.get('ov_hv', False)) or 'helioviewer' in str(meta.get('ov_src', '')).lower()
		or 'helioviewer' in meta or any(str(k).lower().startswith('hv_') for k in meta))


RSUN_REF_M = 695700000.  # IAU 2015 nominal solar radius (m); sunpy's default

def _add_rsun_ref(meta):
	"""Write RSUN_REF explicitly when a header lacks it; return True if added.

	SECCHI Level-0.5 (COR1/COR2/HI/EUVI) and some JP2 headers give RSUN
	(arcsec) but not RSUN_REF. Sunpy then silently assumes the nominal
	radius; recording it in the header keeps the value visible and stable.
	The provenance is kept in OV_RSREF.
	"""
	if meta.get('rsun_ref') not in (None, ''):
		return False
	meta['rsun_ref'] = RSUN_REF_M
	meta['ov_rsref'] = 'RSUN_REF absent in source; IAU 2015 nominal radius written'
	return True


def load_map(source):
	"""Load a FITS/JP2 path or Map with a truthful WCS, whatever its origin.

	Single entry point for downloader, overlay and flow code. Archive FITS
	(VSO/SDAC/SSC/JSOC, Level 0.5/1/1.5) already carry the correct roll in
	CROTA/PC and pass unchanged. Helioviewer pixels are always north-up, so a
	Helioviewer-derived image still carrying the spacecraft roll (files saved
	before 2026-11-10) gets CROTA=0 in memory, old value in OV_ROLL. LASCO split
	DATE-OBS/TIME-OBS is merged. SOHO images without an observer position get
	SOHO's real position (orbit file; L1 approximation offline) and RSUN_REF.
	Any other header without RSUN_REF gets the nominal radius written
	(_add_rsun_ref). Consumers then only follow the WCS.
	"""
	with _quiet():
		m = source if isinstance(source, sunpy.map.GenericMap) else sunpy.map.Map(source)
	m = _lasco_obstime(m)
	meta = m.meta
	if _is_soho(meta) and not all(k in meta for k in ('hgln_obs', 'hglt_obs', 'dsun_obs')):
		fixed = _add_soho_observer(copy.deepcopy(meta), m.date)
		with _quiet():
			m = sunpy.map.Map(m.data, fixed)
		meta = m.meta
	if _is_helioviewer(meta) and 'ov_roll' not in meta:
		roll = any(abs(float(meta.get(k, 0) or 0)) > 1e-6 for k in ('crota', 'crota1', 'crota2'))
		pc = abs(float(meta.get('pc1_2', 0) or 0)) > 1e-8 or abs(float(meta.get('pc2_1', 0) or 0)) > 1e-8
		if roll or pc:
			fixed = _hv_north_up(copy.deepcopy(meta))
			with _quiet():
				m = sunpy.map.Map(m.data, fixed)
	if 'rsun_ref' not in m.meta:
		fixed = copy.deepcopy(m.meta)
		_add_rsun_ref(fixed)
		with _quiet():
			m = sunpy.map.Map(m.data, fixed)
	return m


def _save_maps(maps, outdir, prefix):
	"""Write finished Maps with their complete WCS, retaining the original downloads.

	Helioviewer-derived Maps are saved with the rotation SunPy used in memory,
	so reloading the FITS never re-applies the spacecraft roll.
	"""
	folder = Path(outdir).expanduser().resolve()
	folder.mkdir(parents=True, exist_ok=True)
	for m in maps:
		if 'helioviewer' in m.meta and np.allclose(m.rotation_matrix, np.eye(2)):
			_hv_north_up(m.meta)
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

def _stereo_matches(m,spacecraft,instrument,wave,polar):
	"""Check spacecraft, detector, wavelength and COR1 polarizer."""
	obs=str(m.meta.get('obsrvtry',m.observatory)).upper()
	obs=re.sub(r'[^A-Z0-9]','',obs)
	if obs not in ('STEREO'+spacecraft,
		'STEREO'+('AHEAD' if spacecraft=='A' else 'BEHIND')):
		return False
	if str(m.detector).strip().upper()!=instrument:
		return False
	if m.data.ndim!=2:
		return False
	if instrument=='EUVI':
		return np.isclose(m.wavelength.to_value(u.angstrom),wave)
	if instrument=='COR1':
		return np.isclose(float(m.meta.get('polar',np.nan)),polar)
	if instrument=='COR2':
		# Total-brightness (d*c2 files, POLAR=1001) only; never mix in the
		# single-polarizer n*c2 sequences, whose DN are not comparable.
		return np.isclose(float(m.meta.get('polar',np.nan)),1001)
	if instrument in ('HI1','HI2'):
		# Level-0.5 science sequences only (onboard sums, N_IMAGES>1, s*h1 files);
		# single full-resolution or 256-pixel engineering frames are skipped.
		# Level-1 HI files are handled by _hi_level1_maps, never mixed in here.
		return (int(m.meta.get('n_images',1) or 1)>1 and min(m.data.shape)>=512
			and _hi_level(m.meta)=='0.5')
	return True

def _hi_level(meta):
	"""'1' for HI Level-1/2 products (DN/s, MSB, S10 or hi_prep history), else '0.5'."""
	unit=str(meta.get('bunit','')).strip().upper().replace(' ','')
	if unit and unit not in ('DN','ADU','COUNTS','COUNT'):return '1'
	level=str(meta.get('level',meta.get('data_lev',meta.get('lvl_num','')))).strip().upper()
	if level.lstrip('L').startswith(('1','2')):return '1'
	return '1' if 'hi_prep' in str(meta.get('history','')).lower() else '0.5'

def _hi_level1_maps(sT,eT,dT,spacecraft,instrument,root):
	"""Local HI Level-1 FITS under root, one per cadence bin (nearest DATE-OBS).

	HI Level-1 (hi_prep: DN/s, MSB or S10) is distributed only by UKSSDC
	(https://www.ukssdc.ac.uk/solar/stereo/data.html), which needs a registered
	login, so it is not downloaded here. Files you fetch from there and place
	anywhere under root (e.g. root/raw/level1) are used first.
	Returns (maps, covered_bins, n_bins).
	"""
	targets=_target_times(sT,eT,dT)
	start,end=_parse_time(sT),_parse_time(eT)
	tolerance=float(dT)/2
	found={}
	for path in (_fits_paths(root) if root.is_dir() else []):
		try:
			h=fits.getheader(path)
			if str(h.get('DETECTOR','')).strip().upper()!=instrument or _hi_level({k.lower():v for k,v in h.items()})!='1':continue
			m=sunpy.map.Map(path)
			obs=re.sub(r'[^A-Z0-9]','',str(m.meta.get('obsrvtry',m.observatory)).upper())
			if obs not in ('STEREO'+spacecraft,'STEREO'+('AHEAD' if spacecraft=='A' else 'BEHIND')):continue
			if start<=m.date<=end:found.setdefault(float(m.date.mjd),m)
		except Exception as exc:
			warnings.warn(f'Cannot use HI Level-1 candidate {path}: {exc}')
	chosen={}
	for t in targets:
		near=[m for m in found.values() if abs((m.date-t).to_value(u.s))<=tolerance]
		if near:
			m=min(near,key=lambda m:abs((m.date-t).to_value(u.s)))
			chosen[float(m.date.mjd)]=m
	covered=sum(any(abs((m.date-t).to_value(u.s))<=tolerance for m in chosen.values()) for t in targets)
	return [chosen[k] for k in sorted(chosen)],covered,len(targets)

HI_CADENCE={'HI1':2400.,'HI2':7200.}  # nominal science cadence, seconds

def _hi_series(sT,eT,dT,spacecraft,instrument,outdir):
	"""HI route Level-1 -> Helioviewer -> Level-0.5, as for COR2.

	Level-1: local files only (see _hi_level1_maps). Helioviewer: no HI
	source exists there, so the step is reported and skipped. Level-0.5:
	SSC/VSO through _stereo_series. Bins are max(dT, nominal HI cadence).
	Returns the maps of one family: Level-1 if it covers every bin,
	otherwise the family covering more bins (ties: Level-1).
	"""
	step=max(float(dT),HI_CADENCE[instrument])
	root=Path(outdir).expanduser().resolve()/f'STEREO_{spacecraft}'/instrument
	level1,covered,n=_hi_level1_maps(sT,eT,step,spacecraft,instrument,root)
	print(f'STEREO_{spacecraft}/{instrument}: local Level-1 covers {covered}/{n} bins')
	if covered==n:return level1
	print(f'STEREO_{spacecraft}/{instrument}: Helioviewer has no HI source; '
		'last fallback: SSC Level-0.5 DN')
	l05=_stereo_series(sT,eT,step,spacecraft,instrument,None,outdir)
	return level1 if covered>=len(l05) and level1 else l05

def _hv_jp2_bytes(date,source_id,timeout=(20,120),tries=3):
	"""Return one Helioviewer JP2 (bytes) for datetime `date` and `source_id`.

	Same request as hvpy.getJP2Image, but with a network timeout (connect,
	read) in seconds and `tries` attempts: hvpy sets no timeout, so a
	stalled connection would otherwise hang the run forever.
	"""
	import requests
	from hvpy.parameters import getJP2ImageInputParameters
	p=getJP2ImageInputParameters(date=date,sourceId=source_id)
	for k in range(tries):
		try:
			r=requests.get(p.url,params=p.model_dump(),timeout=timeout);r.raise_for_status()
			return r.content
		except (requests.ConnectionError,requests.Timeout) as exc:
			if k==tries-1:raise
			warnings.warn(f'Helioviewer request stalled ({exc}); retry {k+1}/{tries-1}.')
			time.sleep(10*(k+1))

def _fill_helioviewer_jp2(sT,eT,dT,maps,folder,instrument,wave=None,
                         spacecraft='A'):
	"""Save Helioviewer JP2s; return the saved JP2 paths in the window.

	maps: FITS Maps already covering some bins; only the other bins are
	requested. Pass [] to request a JP2 for every cadence bin.
	instrument 'LASCO' takes wave='C2'/'C3'; files are LASCO_<det>_*.jp2,
	other instruments hv<sourceId>_*.jp2. Cached JP2s are reused.
	"""
	targets=_target_times(sT,eT,dT)
	start,end=_parse_time(sT),_parse_time(eT)
	tolerance=float(dT)/2
	fits_times=[m.date for m in maps]
	missing=[t for t in targets if not any(
		abs((obs-t).to_value(u.s))<=tolerance for obs in fits_times)]
	if not missing:
		return []

	try:
		import hvpy
		from hvpy.datasource import DataSource
	except ImportError:
		warnings.warn('Helioviewer fallback requires hvpy.')
		return []

	if instrument=='LASCO':
		source=getattr(DataSource,f'LASCO_{str(wave).upper()}',None)
		source_id=source.value if source is not None else None
	elif instrument=='AIA':
		source=getattr(DataSource,f'AIA_{int(wave)}',None)
		source_id=source.value if source is not None else None
	elif instrument=='HMI':
		source_id={'continuum':18,'magnetogram':19}.get(wave)
	elif instrument=='EUVI':
		source_id={171:20,195:21,284:22,304:23}[int(wave)]
		if spacecraft=='B':
			source_id+=4
	elif instrument in ('COR1','COR2'):
		source_id={('A','COR1'):28,('A','COR2'):29,
		           ('B','COR1'):30,('B','COR2'):31}[spacecraft,instrument]
	else:
		source_id=None
	if source_id is None:
		warnings.warn(f'No Helioviewer fallback configured for {instrument}/{wave}.')
		return []

	folder=Path(folder)
	folder.mkdir(parents=True,exist_ok=True)
	prefix=f'LASCO_{str(wave).upper()}' if instrument=='LASCO' else f'hv{source_id}'
	cached={}
	for path in folder.glob(f'{prefix}_*.jp2'):
		try:
			obs=sunpy.map.Map(path).date
			if start<=obs<=end:
				cached[obs.isot]=(obs,path)
		except Exception as exc:
			warnings.warn(f'Cannot read cached JP2 {path}: {exc}')

	selected={}
	fits_epochs={obs.isot for obs in fits_times}
	for target in missing:
		near=[(obs,path) for obs,path in cached.values()
		      if obs.isot not in fits_epochs
		      and abs((obs-target).to_value(u.s))<=tolerance]
		if near:
			obs,path=min(near,key=lambda item:abs(
				(item[0]-target).to_value(u.s)))
			selected[obs.isot]=path
			continue

		staged=folder/f'.{prefix}_{target.strftime("%Y%m%dT%H%M%S")}.part.jp2'
		request=target
		for attempt in range(2):
			# Helioviewer returns the image nearest the request; if that falls
			# outside the window, retry once nudged half a tolerance inward.
			try:
				staged.write_bytes(_hv_jp2_bytes(request.datetime,source_id))
				obs=sunpy.map.Map(staged).date
				if not start<=obs<=end or abs((obs-target).to_value(u.s))>tolerance:
					step=-1 if obs>target else 1
					request=target+step*tolerance/2*u.s
					if not start<=request<=end:
						break
					continue
				if obs.isot in fits_epochs:
					break
				path=folder/f'{prefix}_{obs.strftime("%Y%m%dT%H%M%S%f")}.jp2'
				if not path.exists():
					staged.replace(path)
				cached[obs.isot]=(obs,path)
				selected[obs.isot]=path
				break
			except Exception as exc:
				warnings.warn(f'Helioviewer {instrument}/{wave} at {target.isot}: {exc}')
				break
			finally:
				staged.unlink(missing_ok=True)
	return [selected[t] for t in sorted(selected)]

def _stereo_series(sT,eT,dT,spacecraft,instrument,wave,outdir,polar=0):
	"""Reuse local FITS; fetch missing cadence bins through Fido/VSO.

	Returns time-sorted unique Maps. Keeps original metadata and pixels.
	No SECCHI_PREP, exposure normalisation or polarization combination.
	"""
	targets=_target_times(sT,eT,dT)
	start,end=_parse_time(sT),_parse_time(eT)
	root=Path(outdir).expanduser().resolve()/f'STEREO_{spacecraft}'/instrument
	if instrument=='EUVI':
		root=root/str(int(wave))
	raw=root/'raw'
	finished=root/'maps'
	raw.mkdir(parents=True,exist_ok=True)
	finished.mkdir(parents=True,exist_ok=True)
	tolerance=float(dT)/2
	cached={}
	selected={}
	messages=[]

	#### Read existing files first ####
	for path in _fits_paths(root):
		try:
			m=sunpy.map.Map(path)
			if not _stereo_matches(m,spacecraft,instrument,wave,polar):
				continue
			if start<=m.date<=end:
				cached.setdefault(float(m.date.mjd),(m,path))
		except Exception as exc:
			messages.append(f'Cannot use {path}: {exc}')

	missing=[]
	for target in targets:
		near=sorted(cached.values(),
			key=lambda item:abs((item[0].date-target).to_value(u.s)))
		if near and abs((near[0][0].date-target).to_value(u.s))<=tolerance:
			m,path=near[0]
			selected[float(m.date.mjd)]=(m,path)
		else:
			missing.append(target)

	#### Search once, only if some bins are missing ####
	tables=[]
	if missing:
		attrs=[a.Source(f'STEREO_{spacecraft}'),
			a.Instrument('SECCHI'),a.Detector(instrument)]
		if instrument=='EUVI':
			attrs.append(a.Wavelength(float(wave)*u.angstrom))
		if instrument!='COR1':
			attrs.append(a.Physobs('intensity'))
		try:
			query=Fido.search(a.Time(start,end),*attrs)
			tables=[tab for tab in query if len(tab)]
		except Exception as exc:
			messages.append(f'STEREO search failed: {exc}')

	#### Try nearest records, then another provider if necessary ####
	attempted=set()
	for target in missing:
		found=None
		for provider,tab in enumerate(tables):
			distance=np.abs((Time(tab['Start Time'])-target).to_value(u.s))
			for index in np.argsort(distance):
				index=int(index)
				if distance[index]>tolerance:
					break
				record=(provider,index)
				if record in attempted:
					continue
				attempted.add(record)
				try:
					paths=_download(tab[index:index+1],raw)
					for path in paths:
						m=sunpy.map.Map(path)
						if not _stereo_matches(m,spacecraft,instrument,wave,polar):
							continue
						if not start<=m.date<=end:
							continue
						cached[float(m.date.mjd)]=(m,Path(path))
				except Exception as exc:
					messages.append(f'{target.isot}: download failed: {exc}')

				near=sorted(cached.values(),
					key=lambda item:abs((item[0].date-target).to_value(u.s)))
				if near and abs((near[0][0].date-target).to_value(u.s))<=tolerance:
					found=near[0]
					break
			if found is not None:
				break
		if found is None:
			message=f'Missing STEREO_{spacecraft}/{instrument}: {target.isot}'
			print(message)
			messages.append(message)
		else:
			m,path=found
			selected[float(m.date.mjd)]=(m,path)

	#### Save selected Maps, preserving spacecraft geometry ####
	maps=[]
	for key in sorted(selected):
		m,path=selected[key]
		name=f'STEREO_{spacecraft}_{instrument}'
		if instrument=='EUVI':
			name+=f'_{int(wave)}'
		elif instrument=='COR1':
			name+=f'_pol{polar:g}'
		destination=finished/f'{name}_{m.date.strftime("%Y%m%d_%H%M%S")}.fits'
		meta=copy.deepcopy(m.meta)
		if _add_rsun_ref(meta):
			with _quiet():
				m=sunpy.map.Map(m.data,meta)
		if not destination.exists():
			m.save(destination)
		maps.append(m)
		messages.append(
			f'{m.date.isot} | input={path} | saved={destination} | '
			f'level={m.processing_level} | unit={m.unit}')

	log=root/'retrieval_log.txt'
	with log.open('a') as stream:
		stream.write(f'\nRequest: {start.isot} -- {end.isot}; dT={dT}s\n')
		stream.write('\n'.join(messages)+'\n')
	return maps
	
def _matching_files(directories, kind, key):
	paths = []
	for directory in directories:
		root = Path(directory).expanduser().resolve()
		candidates = list(_fits_paths(directory))
		if kind in ('AIA', 'HMI') and root.is_dir():
			candidates += sorted(p for p in root.rglob('*')
				if p.suffix.lower() in ('.jp2', '.j2k'))
		for path in candidates:
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
		raise FileNotFoundError(f'No {kind}/{key} maps in supplied directories')
	return sorted(set(paths), key=lambda p: (p.suffix.lower() in ('.jp2', '.j2k'), str(p)))


def _lasco_exposure_id(image):
    """Native exposure identifier from a LASCO archive filename, when present."""
    name = str(image.meta.get('filename', '')).rsplit('/', 1)[-1]
    match = re.fullmatch(r'[23]\d(\d{6})\.(?:fts|fits)', name, re.IGNORECASE)
    return match[1] if match else None
