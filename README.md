# Solar Data Mining

Python tools for mining and cross-matching multi-wavelength solar flare data for event dates of interest. The repository supports retrieval and inspection of radio, EUV, magnetic-field, coronagraph and high-energy solar datasets used in solar flare and space weather studies.

## Purpose

This repository contains research utilities for compiling event-based solar flare datasets and producing multi-instrument image overlays. `Multiwb_datamine.py` retrieves and prepares observations; `solar_overlay.py` combines prepared images, difference images and potential-field source-surface (PFSS) field lines. The remaining scripts support radio analysis, active-region cutouts, observing catalogs and event matching.

A two-page flow chart of how SDO, LASCO and STEREO images are downloaded, which sources are tried, what the fallbacks are and how spacecraft roll is handled is in [`Image_download_and_roll_flowchart.pdf`](Image_download_and_roll_flowchart.pdf).

## Main capabilities

- Query and process event-based multi-wavelength solar observations.
- Retrieve numerical radio spectra and archive quicklooks from RSTN, e-CALLISTO, Learmonth, Culgoora, Wind/WAVES and STEREO/WAVES.
- Download and prepare SDO/AIA, SDO/HMI, SOHO/LASCO, STEREO/SECCHI (EUVI, COR1, COR2, HI1) and official SOHO/EIT Level-1 images, with Helioviewer JP2 back-ups for every cadence slot.
- Load any image (archive FITS, Level-1, Helioviewer JP2 or a FITS made from one) with a truthful, correctly rolled WCS through one function, `load_map`.
- Generate aligned image overlays, running/base differences and PFSS projections from local FITS files or explicitly enabled downloads.
- Download and plot NoRP GHz light curves, polarization and flare summaries; retrieve NoRH images and RHESSI quicklooks.
- Download GOES soft-X-ray flux and plot its time derivative.
- Match MWA observations with HEK events and retrieve observing metadata.
- Download SDO/AIA active-region cutouts and generate MWA observing catalogs.

## Repository contents

| File | Purpose |
|---|---|
| `Multiwb_datamine.py` | Main observation-retrieval module: radio spectra/light curves, SDO/LASCO/STEREO/EIT image retrieval, Helioviewer JP2 back-ups, GOES and RHESSI products, and co-observation searches. Detailed API below. |
| `solar_overlay.py` | Main visualization module: multi-instrument PNG frames, aligned differences and static PFSS models/projections. Detailed API below. |
| `datamine_calib_helper.py` | Archive search, nearest-exposure selection, FITS validation, LASCO calibration/fallback support, Helioviewer JP2 fetching, `load_map` and `jp2_to_fits`. Keep beside the main modules. |
| `solar_overlay_helper.py` | Internal dataset parsing, FITS selection, reprojection, image compositing, frame rendering and PFSS plotting. Also imported by the frame-maker and flow-tracking projects. Keep beside the main modules. |
| `Image_download_and_roll_flowchart.pdf` | Two-page flow chart: download order, fallbacks and roll handling for SDO, LASCO and STEREO. |
| `Matrix_ops.py` | `Reg` regularises matrix axes and marks gaps. Required by the Wind/STEREO numerical-spectrum workflow; the current implementation explicitly uses it for Wind time regularisation. |
| `get_radioDS_fordates.py` | Standalone radio dynamic-spectrum retrieval for selected dates. |
| `downAIAcutouts.py` | Standalone SDO/AIA active-region cutout and light-curve workflow. Inspect its event-table settings and local paths before running. |
| `downNorP_GHzLC_plot.py` | Standalone NoRP GHz light-curve retrieval, plotting and flare-analysis workflow. |
| `MWA_HEK_match.py` | Matches MWA observations with HEK solar-event information. |
| `MWA_obscat_maker.py` | Creates MWA observing-catalog files. |
| `hessi_flare_catalog.csv` | RHESSI flare catalog used by catalog-selection and quicklook utilities. |
| `time_conv.py` | Helper routines for time conversion. |

## Requirements and setup

The scripts are modules rather than an installable package. Clone the repository and run from its directory, or put that directory on your Python import path (`PYTHONPATH`). Keep the four main/helper files together, with the exact case-sensitive filenames listed above.

The supplied code imports NumPy, Astropy, SunPy, Matplotlib, pandas, Pillow, Beautiful Soup, aiapy, radiospectra, ecallisto-ng, SpacePy, multiprocess and certifi. Overlays also require reproject; PFSS generation requires pfsspy; Helioviewer JP2 back-ups require hvpy (and glymur/OpenJPEG to read JP2). SunPy network/map/timeseries extras install supporting dependencies. Excel output may additionally need an Excel writer such as openpyxl. These are dependency names, not a tested version lock.

```bash
python -m pip install 'sunpy[net,map,timeseries]' numpy astropy matplotlib pandas pillow beautifulsoup4 aiapy radiospectra ecallisto-ng spacepy multiprocess certifi reproject pfsspy hvpy glymur openpyxl
```

Use a compatible scientific-Python environment. The overlay helper imports `sunpy.coordinates.SphericalScreen`, so the installed SunPy must provide it. SpacePy's `pycdf` requires a working NASA CDF backend. Downloads depend on archive availability and instrument coverage. Full local LASCO calibration additionally requires a configured SolarSoft/IDL installation and `sswidl` on PATH; downloaded verified Level-1 products do not require IDL.

Times are UTC, normally `YYYY-MM-DDThh:mm:ss`. Image cadence `dT` is in seconds; MWA search duration `DT` is in minutes. Requested cadence selects available exposures and cannot create observations during archive gaps.

## Images: sources, folders and orientation

### Where images come from

| Source | Used for | What it is |
|---|---|---|
| Fido/VSO (JSOC, SDAC, SSC data archive) | AIA, HMI, LASCO, SECCHI EUVI/COR1/COR2/HI1 | Archive FITS with full WCS |
| UKSSDC (registered login; not downloaded by the code) | SECCHI HI Level-1/2 (`hi_prep`, DN/s, MSB, S10) | Place the files under `STEREO_<sc>/HI1/` (e.g. `raw/level1/`); they are then used first |
| LASCO team archive (`img_hdr.txt` index) | LASCO, when the VSO query fails | Level-0.5 FITS |
| Official EIT Level-1 archive | EIT | Calibrated Level-1 FITS |
| Helioviewer API (`hvpy`) | AIA, HMI, LASCO, EUVI, COR1, COR2 | JP2 display images (0-255), already rotated north-up, with the original FITS header embedded |
| Browse JPG/PNG quicklooks (e.g. `stereo-ssc.../browse/`) | Never used as images | No WCS |

### Folder layout

Each image product folder follows the same convention:

| Folder | Holds |
|---|---|
| `raw/` | Original downloads |
| `maps/` | Science FITS only: Level-1, Level-1.5 (AIA), Level-0.5 DN/s (LASCO) or SECCHI Level-0.5 DN |
| `raw/helioviewer*/` | Helioviewer JP2s, one per cadence slot when `helioviewer=True` (the default) |
| `maps_HV/` | FITS copies of the JP2s made by `jp2_to_fits`, on demand |

SDO products live under `outdir/<AIA|HMI>/<wave>/`; STEREO products under `outdir/STEREO_<A|B>/<EUVI|COR1|COR2|HI1>[/<wave>]/`.

### Product route: Level-1 -> Helioviewer -> Level-0.5

When several products exist for the same images, every program prefers them in this order: verified Level-1 (calibrated), then Helioviewer, then Level-0.5. LASCO follows it when `need_calib`/`need_LASCO_calib` is True (False keeps Level-0.5 first). STEREO always follows it; SECCHI Level-1 (`secchi_prep`) is not distributed, so in practice STEREO uses Helioviewer where it exists and Level-0.5 elsewhere. HI1 follows the same route: local Level-1 files (UKSSDC needs a login, so they are never downloaded), then Helioviewer (it has no HI source, so this step is reported and skipped), then SSC Level-0.5 science sums. The flow code uses the first family in this order that covers every cadence slot, otherwise the family covering the most slots (ties by this order), and never mixes families. Detection thresholds are in the units of the chosen family, so changing family changes what is detected.

### Spacecraft roll: one rule

Every saved FITS carries a header that truthfully describes its pixels, and every consumer loads images through `load_map` and simply follows the WCS. Nothing rotates an image because of its source or level.

| Product | Pixels | Header roll |
|---|---|---|
| Archive FITS (LASCO Level-0.5, SECCHI, HMI) | As observed | The spacecraft roll (e.g. LASCO about 7°, STEREO-A about 6.5–6.9°, HMI about 180°) |
| Level-1 / Level-1.5 (LASCO, AIA, EIT) | Rotated to solar north by the instrument team or aiapy | Matches the pixels; used unchanged |
| Helioviewer JP2 or `maps_HV` FITS | Already north-up | 0, with the original roll kept in `OV_ROLL` |

Helioviewer LASCO JP2s keep the original CROTA in their header. SunPy ignores it only while a JP2-only `helioviewer` key is present, and that key is lost when a FITS is written. FITS saved by versions of this module before 2026-11-10 therefore re-apply the roll (about 7° counter-clockwise for C3, 7.7° for C2). `load_map` corrects such files in memory, and `_save_maps`/`jp2_to_fits` now write roll 0.

## Multiwb_datamine.py

This is the retrieval workhorse. Public functions can be imported independently, although module-level imports require the dependencies above. Importing the module does not start an event-processing workflow.

### Image retrieval and preparation

| Function | Inputs and behavior | Returns / saved products |
|---|---|---|
| `get_SDO_imgs(sT,eT,dT,`<br>`instruments=['AIA','HMI'],`<br>`waves=[[94,171,193],['magnetogram']],`<br>`outdir='./SDO',nproc=4,`<br>`helioviewer=True)` | One job per instrument/product. `find_coobs` downloads from JSOC/VSO at cadence `dT`; AIA Level 1 is upgraded to Level 1.5 (`update_pointing` + `register`, DN/s, north-up); HMI is kept as delivered (CROTA2 about 180°). The image nearest each cadence slot is kept. With `helioviewer=True` a JP2 is saved for every slot. | A list of Map lists in input order. FITS in `outdir/<INSTR>/<wave>/`, downloads in `raw/`, JP2s in `raw/helioviewer/` (not returned). |
| `get_LASCO_imgs(sT,eT,dT,instr,`<br>`outdir='./LASCO_1',`<br>`need_calib=False,`<br>`helioviewer=True)` | `instr='C2'` or `'C3'`. `need_calib=True` follows the route **Level-1 -> Helioviewer -> Level-0.5**: verified Level-1/MSB from VSO/SDAC (saved unchanged in `maps/`); if none is usable, Helioviewer images (FITS in `maps_HV/`, canonical tags, roll 0, SOHO observer); only if Helioviewer also fails, archive Level-0.5 DN/s. `need_calib=False`: archive Level-0.5 (VSO, else the LASCO team archive), exposure-normalised to DN/s and FOV-masked, in `maps/`; with `helioviewer=True` a JP2 is also saved for every slot. Only standard 1024² frames (C2 orange, C3 clear) are accepted. | Map list; downloads in `raw/`, JP2s in `raw/helioviewer_<instr>/`. |
| `get_STEREO_imgs(sT,eT,dT,spacecraft='A',`<br>`instruments=['COR2'],waves=None,`<br>`outdir='./STEREO',nproc=4,polar=0,`<br>`helioviewer=True)` | SECCHI Level-0.5 FITS (`secchi_rectify` headers, DN) from the SSC data archive through VSO. Local files are reused; only uncovered slots are searched. EUVI is filtered by wavelength, COR1 by polarizer `polar` (0/120/240), COR2 to total brightness (POLAR 1001, `d*c2` files) only. No `secchi_prep`, exposure normalisation or polarizer combination. `helioviewer=True` saves a JP2 for every slot. **HI1** (`instruments=['HI1']`): route Level-1 -> Helioviewer -> Level-0.5. Local Level-1 HI files are used if they cover every bin; Helioviewer has no HI source; otherwise SSC Level-0.5 science sums (`s*h1*.fts`, 30 summed exposures, 1024 px, DN) are fetched. Bins are max(`dT`, 2400 s), the HI1 cadence. One family is returned (Level-1 if complete, else the one covering more bins). | One Map list per instrument/wave. FITS in `STEREO_<sc>/<instr>[/<wave>]/maps/`, downloads in `raw/`, JP2s in `raw/helioviewer/`, a `retrieval_log.txt`. |
| `calib_LASCO(images,instr=None,`<br>`outdir='./LASCO_calibrated')` | Existing FITS path/Map or list/tuple. Verified Level-1/MSB passes through. Otherwise tries optional SolarSoft, remote Level-1, JP2 and finally an explicitly uncalibrated preview. Remote replacements must meet time/exposure checks. | Map or list matching input form; saved FITS/logs. Preserve original files referenced by `OV_RAW` for local calibration. |
| `get_EIT_imgs(sT,eT,dT,wave,`<br>`outdir='./EIT',archive_url=...)` | `wave=171,195,284,304`. Retrieves official preprocessed Level-1 products and rejects raw EIT substitutions. | Calibrated SunPy Maps and selected FITS. Availability and native cadence limit sampling. |
| `get_helioviewer_jp2(sT,eT,dT,instrument,`<br>`wave=None,spacecraft='A',outdir)` | Saves a Helioviewer JP2 for every cadence slot: `'LASCO'` (`wave='C2'/'C3'`), `'AIA'`, `'HMI'` (`'magnetogram'/'continuum'`), `'EUVI'`, `'COR1'`, `'COR2'`. If the nearest image falls outside the window, the request is retried once nudged inward. Cached JP2s are reused. | List of JP2 paths. |
| `jp2_to_fits(jp2_paths,outdir,`<br>`skip_times=())` | Converts JP2s to FITS (conventionally into `maps_HV/`): Helioviewer display values, roll 0 with the old value in `OV_ROLL`, STEREO observer kept (SOHO observer added for LASCO), LASCO occulter masked. Every JP2-to-FITS route (this function, the LASCO download route and the frame makers) writes the same canonical tags `OV_HV=True`, `OV_SRC='Helioviewer JP2'`, `OV_CAL='Helioviewer visualization; not photometric FITS'`, `OV_NORM=True`, `BUNIT='ct'`, so all Helioviewer FITS form one product family. Existing outputs are reused; `outdir` is created only when something is written. | List of FITS paths. |
| `load_map(source)` | Loads a FITS/JP2 path or Map with a truthful WCS: archive and Level-1 headers unchanged, stale rolls on Helioviewer-derived files set to 0 in memory, LASCO split DATE-OBS/TIME-OBS merged. SOHO images (LASCO, EIT) without an observer position get SOHO's real HGLN/HGLT/DSUN_OBS from the SOHO orbit file (cached in `~/.cache/solar_data_mining/soho_orbit`; L1 approximation offline) plus `RSUN_REF`, recorded in `OV_OBS`. Use it in every consumer. | SunPy Map. |

```python
from pathlib import Path
import Multiwb_datamine as dm

if __name__=='__main__':
    sT,eT='2025-05-30T22:30:00','2025-05-31T03:00:00'
    aia193,hmi=dm.get_SDO_imgs(sT,eT,1800,instruments=['AIA','HMI'],
        waves=[[193],['magnetogram']],outdir='./SDO',nproc=1)
    c3=dm.get_LASCO_imgs(sT,eT,720,'C3',outdir='./LASCO/C3',need_calib=False)
    cor2,=dm.get_STEREO_imgs(sT,eT,900,spacecraft='A',instruments=['COR2'],
        outdir='./STEREO',nproc=1)
    # Helioviewer copies when the science FITS leave gaps:
    hv=dm.jp2_to_fits(sorted(Path('./LASCO/C3/raw/helioviewer_C3').glob('*.jp2')),'./LASCO/C3/maps_HV')
    m=dm.load_map(hv[0])          # always load through load_map
```

**Multiprocessing:** `nproc=1` runs serially. Larger values use `multiprocess.Pool` across instrument/product combinations, not individual frames. The pool is scoped with `with Pool(nproc)`; returned Maps and each worker's scientific-Python imports can still consume substantial memory. Start with serial mode for debugging and choose worker count according to available memory. On macOS, call parallel workflows from a script under `if __name__=='__main__':`. The nested worker relies on `multiprocess` serialization; replacing the import with standard-library `multiprocessing.Pool` is not a drop-in change. Overlay's own SDO acquisition calls explicitly use `nproc=1`.

**Cost of JP2 back-ups:** one JP2 per slot is about 0.3–1 MB for LASCO/STEREO and 2–4 MB for full-disk AIA/HMI. Set `helioviewer=False` to skip them.

### Radio spectra and light curves

| Function | Purpose and outputs |
|---|---|
| `get_norp_LC(sT,eT,outdir='./',`<br>`plot=True)` | Downloads NoRP FITS and derives microwave Stokes-I/V light curves, polarization, peak summaries and flare quantities. Writes event-subfolder data, CSV/Excel/pickle summaries and optional figures. Returns `None`. |
| `get_norh_imgs(sT,eT,outdir='./')` | Retrieves NoRH image products and saves quicklook plots. Returns `None`. |
| `get_RSTN(sT,eT,outdir='./',`<br>`downdata=True,plot=True)` | Attempts numerical Culgoora/Learmonth, e-CALLISTO and RSTN sources. Australian raw retrieval targets intervals overlapping 22–10 UT and uses the start-date filename convention. Archive quicklooks are a fallback when numerical files are unavailable. Returns paths, plots and errors. `downdata=False` reports candidates. |
| `radspec_to_pickle(filname,`<br>`outdir='./')` | Converts SRS/FITS spectrograms or e-CALLISTO DataFrame pickles into saved `.p` data and PNGs; returns both paths. DataFrame inputs retain that format. |
| `plot_culgoora_spec(filename,`<br>`outdir='./')` | Decodes raw `SPECYYMMDD.bz2` records using recorded timestamps and four bands spanning 18–1800 MHz. Saves spectrum/metadata pickle and PNG, returning both paths. Native amplitude codes are not calibrated SFU. |
| `get_DHband_quicklookDS(sT,eT,`<br>`outdir='./')` | Retrieves Wind/STEREO archive images and assembles available quicklook panels. These are visualization products, not numerical spectra. |
| `get_DHband_calibDS(sT,eT,`<br>`outdir='./',`<br>`sats=('sta','stb','wind'),`<br>`Ch_frq=1.0,stats='median',`<br>`vmin=None,vmax=None)` | Downloads multi-day SPDF CDF data; uses STEREO L3 SFU and the supplied Wind Z-channel gain conversion. Produces dynamic spectra and nearest-channel light curves, recording the actual channel frequency. `Ch_frq` is MHz; `vmin/vmax` are log10(SFU) display limits. Returns paths/errors per spacecraft. |

Numerical DH products are organised as `CDFs/`, `DSpickle/` and `PNGs/`, each with `STEREO_A/`, `STEREO_B/` and `Wind/` subdirectories. Existing CDFs/directories are reused.

`Matrix_ops.Reg(DS,R,C,r0,r1,c0,c1,reg_ax='both',stats='min')` takes a 2-D matrix and its row/column coordinates. `reg_ax='row'`, `'col'` or `'both'` selects which axes to regularise; `stats='min'`, `'median'` or `'mean'` selects the spacing statistic. It assigns nearest input samples to the regular grid, leaves sufficiently distant samples as NaNs and returns `(DSn,Rgs,Cgs)` with masked coordinate arrays. This is nearest-sample regridding, not smoothing or flux calibration.

```python
result=dm.get_DHband_calibDS('2014-11-03T06:00:00','2014-11-03T06:16:00',
    outdir='./radio',sats=('sta','wind'),Ch_frq=1.0,stats='median')
```

### Searches and high-energy utilities

| Function | Purpose and outputs |
|---|---|
| `MWAavail(sT,DT=30,outdir='./')` | Queries MWA observations over a duration in minutes, saves an observation-table pickle and returns success/failure. |
| `find_flare(sT,eT,`<br>`flare_cat='hessi_flare_catalog.csv')` | Returns `(IDs,Tab)` for catalog flares fully contained in the interval, including boundary times, sorted by start time. Accepts CSV or a pickled Astropy Table. |
| `find_RHESSI(sT,eT,outdir='./',`<br>`flare_cat=...,baseurl=...)` | Saves the selected flare catalog and matching RHESSI quicklook files, plus available daily orbit panels. Returns `(IDs,Tab)`. |
| `plot_goes_derivative(sT,eT,`<br>`outdir='./',satellite=None,`<br>`savedata=True)` | Retrieves GOES 1–8 Angstrom flux and computes its derivative using actual sample intervals. `satellite=None` selects the highest available satellite number. Saves a two-panel PDF and optional data pickle (`Time`, `Flux`). |
| `find_coobs(sT,eT,instrs,wvs,`<br>`outdir='./',RHESSI_flarecat=...,`<br>`downdata=False,dT=None,`<br>`HMI_product=None,flat_outdir=False)` | Searches matching instrument/wavelength observations; wavelength lists use Astropy units. Optionally downloads and prepares supported products (AIA to Level 1.5). `HMI_product` filters magnetogram/continuum/dopplergram. Returns failed instrument/wavelength searches. |

## solar_overlay.py

This is the image-combination workhorse. Its public functions are `overlay` and `build_pfss`; rendering and projection details live in `solar_overlay_helper.py`. All image layers, difference references and the target map are loaded through `load_map`.

### Overlay frames

```python
from solar_overlay import overlay

frames=overlay(['LASCO/C2','SDO/AIA94'],
    '2011-06-07T06:00:00','2011-06-07T07:15:00',600,
    image_dirs={'LASCO/C2':'./LASCO/C2/maps','SDO/AIA94':'./SDO/AIA/94'},
    diffs=[True,False],diff_types=['running',''],
    first_frame_intensity=True,outdir='./frames')
```

`overlay` returns PNG paths. Supported dataset names include `LASCO/C2`, `LASCO/C3`, `SDO/AIA94` (and other supported AIA bands), `SDO/HMI`, `SOHO/EIT171/195/284/304` and `pfss`.

| Parameter | Meaning |
|---|---|
| `datasets,sT,eT,dT` | Requested layers, UTC interval and cadence in seconds. |
| `image_dirs=None` | A FITS directory, directory list or dataset-to-directories dictionary. Searches recursively and identifies products from headers; AIA/HMI JP2s in the folder are used too. Point LASCO at `maps/` or `maps_HV/` explicitly, so the two product families are not mixed. Directory mode does not acquire missing images. |
| `local_files=None` | Explicit `{dataset:[FITS paths]}` alternative. Missing keys raise unless acquisition is explicitly enabled. |
| `download=False` | Enables missing-image acquisition only when set to `True`; uses `Multiwb_datamine`. |
| `diffs=None,diff_types=None` | Per-layer difference flags and `'running'`/`'base'` choices. Missing flags are false; PFSS is never differenced. A difference between a science FITS and a Helioviewer image is refused and the frame is listed in `skipped_frames.txt`. |
| `first_frame_intensity=False` | Running differences normally search for a pre-start predecessor and fall back to intensity if absent. `True` makes the first frame intensity without requiring a predecessor. |
| `base_diff_DT=None` | Base reference time or per-dataset dictionary; defaults to the event start. Before the event, selects the latest exposure at/before the reference; after the event, earliest at/after; inside the event, nearest exposure. Missing base references raise. |
| `reference_search=3600` | Maximum distance in seconds for external-reference/predecessor searches. |
| `tmap=None` | Target Map/FITS defining fixed FoV and WCS. Otherwise the bottom image provides the grid and reference epoch. |
| `alpha=0.6,npix=1024,dpi=150,`<br>`figsize=7.0` | Upper-layer maximum opacity, longest display-grid dimension, figure resolution and square figure size in inches. |
| `outdir='./',`<br>`cache_dir='./solar_data'` | Output PNG directory and acquisition/reference cache. |
| `need_LASCO_calib=False` | Default uses Level-0.5 DN/s preparation. `True` asks for verified Level-1 first. Verified Level-1 passes through `_prepare_lasco` untouched. |
| `pfss_dir=None,pfss_kwargs=None,`<br>`N_Bloop=np.nan` | Reuse saved PFSS, or supply model-building settings. NaN keeps all closed loops; an integer retains the highest visible loops; zero shows open lines only. |

Layer order is LASCO C3, LASCO C2, HMI, AIA/EIT image layers and PFSS. Labels report actual layer observation times; upper layers use the image nearest the bottom-layer time with no maximum offset, so check the label. Reference images are reprojected to the current WCS before subtraction; differential solar rotation is not corrected. Prepared AIA Level-1.5 inputs are required and are not calibrated twice. EIT must be official Level-1.

### PFSS model and field-line display

```python
from solar_overlay import build_pfss,overlay

model=build_pfss('2011-06-07T06:00:00','2011-06-07T07:15:00',
    './pfss',rss=2.5,nrho=35,r_seed=1.05,n_lon=36,n_lat=18)
frames=overlay(['LASCO/C2','pfss'],
    '2011-06-07T06:00:00','2011-06-07T07:15:00',600,
    image_dirs='./LASCO/C2/maps',pfss_dir='./pfss',N_Bloop=20,
    outdir='./pfss_frames')
```

`build_pfss(sT,eT,outdir,rss=2.5,nrho=35,r_seed=1.05,n_lon=36,n_lat=18,gong_file=None,cache_dir='./solar_data/gong')` uses the GONG synoptic radial boundary nearest the interval midpoint, or a supplied local `gong_file`. It builds one static PFSS solution, not an evolving MHD model. `rss` and `r_seed` are in solar radii, with `1 < r_seed < rss`; `nrho` controls radial resolution and `n_lon/n_lat` control seed counts.

It returns model/trace information including `outdir` and saves `boundary.fits`, `metadata.json`, `field_lines.npz` and `model_3d.npz`. The 3-D archive contains Bphi/Btheta/Br in gauss with longitude, sine-latitude and radial axes. Plotting projects the co-rotating geometry to the target epoch, hides segments behind the disk and extends open lines radially to the display FoV. Red/blue indicate outward/inward open fields; white indicates closed loops. Dashed segments beyond the source surface are display continuations, not computed PFSS field. The PFSS drawing functions in `solar_overlay_helper` (`_draw_pfss`, `_pfss_to_fov`, `_nan_sep`, `_closed_loop_limit`) are imported by the frame-maker and flow-tracking projects, so changes here reach them.

## Companion code: moving-structure and overlay frames

[`Moving_structure_frames`](https://github.com/atul3790/Solar_data_mining/tree/moving-structure-frames/Moving_structure_frames) builds diagnostic movies from these modules. It follows CMEs and other moving structures in LASCO, STEREO COR2 and HI1 images, and makes overlay frames for LASCO, STEREO and PSP/WISPR. Each frame combines:

- bright and dark structure contours with in-plane drift arrows
- projected PFSS field lines from the observing spacecraft's viewpoint
- an EUV or magnetogram disk image
- Fermi-LAT and GOES light curves

It downloads through `get_LASCO_imgs`, `get_STEREO_imgs` (including HI1) and `get_SDO_imgs`, loads every image with `load_map`, and builds PFSS with `solar_overlay.build_pfss`. Keep this repository on `PYTHONPATH` when using it.

## Scientific context

These tools were developed for event-based solar flare and space weather studies, especially for constructing catalogs connecting flare emission, radio bursts, CMEs and energetic particle events. They support Sun-as-a-star analyses and event visualization across radio, EUV, X-ray and coronagraph observations.

## Notes

- Research scripts may contain site-specific paths and event-table column settings. Inspect standalone script inputs before execution; their APIs may differ from the main modules.
- Helioviewer JP2 and `maps_HV` images, DN/s previews and running differences have different scientific uses from verified radiometrically calibrated products. Helioviewer pixels are display values: good for positions, morphology and apparent motion, not for brightness or mass. Check `BUNIT`, `OV_SRC` and timestamps before quantitative analysis.
- Never difference a science FITS against a Helioviewer image; the overlay and flow code refuse such pairs.
- Headers without a reference solar radius get `RSUN_REF = 695700000 m` (IAU 2015 nominal) written explicitly, with an `OV_RSREF` note, instead of relying on SunPy's silent default. This applies to `load_map`, `jp2_to_fits`, saved STEREO maps and PFSS boundaries.
- Fido/VSO downloads retry up to 3 times (10 s, then 20 s back-off) after gateway time-outs or dropped connections.
- Keep the raw downloads and helper modules when sharing or reusing results. Pickle files should only be loaded from trusted sources.
- The 2026-11-10 changes (roll handling, `load_map`, COR2 filter, JP2 back-ups, `jp2_to_fits`) were tested on live downloads and on files from the 2025-05-30/31 event. Full multi-hour frame-maker runs have not been repeated since.
