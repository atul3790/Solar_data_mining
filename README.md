# Solar Data Mining

Python tools for mining and cross-matching multi-wavelength solar flare data for event dates of interest. The repository supports retrieval and inspection of radio, EUV, magnetic-field, coronagraph and high-energy solar datasets used in solar flare and space weather studies.

## Purpose

This repository contains research utilities for compiling event-based solar flare datasets and producing multi-instrument image overlays. `Multiwb_datamine.py` retrieves and prepares observations; `solar_overlay.py` combines prepared images, difference images and potential-field source-surface (PFSS) field lines. The remaining scripts support radio analysis, active-region cutouts, observing catalogs and event matching.

## Main capabilities

- Query and process event-based multi-wavelength solar observations.
- Retrieve numerical radio spectra and archive quicklooks from RSTN, e-CALLISTO, Learmonth, Culgoora, Wind/WAVES and STEREO/WAVES.
- Download and prepare SDO/AIA, SDO/HMI, SOHO/LASCO and official SOHO/EIT Level-1 images.
- Generate aligned image overlays, running/base differences and PFSS projections from local FITS files or explicitly enabled downloads.
- Download and plot NoRP GHz light curves, polarization and flare summaries; retrieve NoRH images and RHESSI quicklooks.
- Download GOES soft-X-ray flux and plot its time derivative.
- Match MWA observations with HEK events and retrieve observing metadata.
- Download SDO/AIA active-region cutouts and generate MWA observing catalogs.

## Repository contents

| File | Purpose |
|---|---|
| `Multiwb_datamine.py` | Main observation-retrieval module: radio spectra/light curves, SDO/LASCO/EIT image preparation, GOES and RHESSI products, and co-observation searches. Detailed API below. |
| `solar_overlay.py` | Main visualization module: multi-instrument PNG frames, aligned differences and static PFSS models/projections. Detailed API below. |
| `datamine_calib_helper.py` | Internal archive search, nearest-exposure selection, FITS validation, image saving and LASCO calibration/fallback support. Keep beside the main modules. |
| `solar_overlay_helper.py` | Internal dataset parsing, FITS selection, reprojection, image compositing, frame rendering and PFSS plotting. Keep beside the main modules. |
| `Matrix_ops.py` | `Reg` regularises matrix axes and marks gaps. Required by the Wind/STEREO numerical-spectrum workflow; the current implementation explicitly uses it for Wind time regularisation. |
| `get_radioDS_fordates.py` | Standalone radio dynamic-spectrum retrieval for selected dates. |
| `downAIAcutouts.py` | Standalone SDO/AIA active-region cutout and light-curve workflow. Inspect its event-table settings and local paths before running. |
| `downNorP_GHzLC_plot.py` | Standalone NoRP GHz light-curve retrieval, plotting and flare-analysis workflow. |
| `MWA_HEK_match.py` | Matches MWA observations with HEK solar-event information. |
| `MWA_obscat_maker.py` | Creates MWA observing-catalog files. |
| `hessi_flare_catalog.csv` | RHESSI flare catalog used by catalog-selection and quicklook utilities. |
| `time_conv.py` | Helper routines for time conversion. |

## Requirements and setup

The scripts are modules rather than an installable package. Clone the repository and run from its directory, or put that directory on your Python import path. Keep the five main/helper files together, with the exact case-sensitive filenames listed above.

The supplied code imports NumPy, Astropy, SunPy, Matplotlib, pandas, Pillow, Beautiful Soup, aiapy, radiospectra, ecallisto-ng, SpacePy, multiprocess and certifi. Overlays also require reproject; PFSS generation requires pfsspy. SunPy network/map/timeseries extras install supporting dependencies. Excel output may additionally need an Excel writer such as openpyxl. These are dependency names, not a tested version lock.

```bash
python -m pip install 'sunpy[net,map,timeseries]' numpy astropy matplotlib pandas pillow beautifulsoup4 aiapy radiospectra ecallisto-ng spacepy multiprocess certifi reproject pfsspy openpyxl
```

Use a compatible scientific-Python environment. The overlay helper imports `sunpy.coordinates.SphericalScreen`, so the installed SunPy must provide it. SpacePy's `pycdf` requires a working NASA CDF backend. Downloads depend on archive availability and instrument coverage. Full local LASCO calibration additionally requires a configured SolarSoft/IDL installation and `sswidl` on PATH; downloaded verified Level-1 products do not require IDL.

Times are UTC, normally `YYYY-MM-DDThh:mm:ss`. Image cadence `dT` is in seconds; MWA search duration `DT` is in minutes. Requested cadence selects available exposures and cannot create observations during archive gaps. Examples below are usage examples, not claims of successful live archive validation.

## Multiwb_datamine.py

This is the retrieval workhorse. Public functions can be imported independently, although module-level imports require the dependencies above. Importing the module does not start an event-processing workflow. The descriptions below follow the active function definitions; some introductory docstring signatures retain older APIs.

### Image retrieval and preparation

| Function | Inputs and behavior | Returns / saved products |
|---|---|---|
| `get_SDO_imgs(sT,eT,dT,instruments=...,waves=...,outdir='./SDO',nproc=4)` | Instrument list and matching list of product lists. AIA products are numeric wavelengths in Angstrom; HMI products are `magnetogram`, `continuum` or `dopplergram` (`None` selects magnetogram). Defaults request AIA 94/171/193 and HMI magnetogram/continuum. Downloads and AIA preparation are delegated to `find_coobs`. | A list of Map lists, one per instrument/product combination in input order. Finished FITS under `outdir/<instrument>/<product>/`; originals under each product's `raw_SDO/`. |
| `get_LASCO_imgs(sT,eT,dT,instr,outdir='./LASCO_1',need_calib=False)` | `instr='C2'` or `'C3'`. Default prepares archive images as exposure-normalised, field-of-view-masked DN/s previews. With `need_calib=True`, tries verified Level-1/MSB FITS, then Helioviewer JP2, then archive DN/s. | Map list; originals in `raw/`, prepared FITS in `maps/`. JP2 is a display product, not MSB photometry. |
| `calib_LASCO(images,instr=None,outdir='./LASCO_calibrated')` | Existing FITS path/Map or list/tuple. Verified Level-1/MSB passes through. Otherwise tries optional SolarSoft, remote Level-1, JP2 and finally an explicitly uncalibrated preview. Remote replacements must meet time/exposure checks. | Map or list matching input form; saved FITS/logs. Preserve original files referenced by `OV_RAW` for local calibration. |
| `get_EIT_imgs(sT,eT,dT,wave,outdir='./EIT',archive_url=...)` | `wave=171,195,284,304`. Retrieves official preprocessed Level-1 products and rejects raw EIT substitutions. Optional archive URL must point to an equivalent archive. | Calibrated SunPy Maps and selected FITS. Availability and native cadence limit sampling. |

```python
import Multiwb_datamine as dm

if __name__=='__main__':
    sT,eT='2011-06-07T06:00:00','2011-06-07T07:15:00'
    groups=dm.get_SDO_imgs(sT,eT,600,instruments=['AIA','HMI'],
        waves=[[94,193],['magnetogram']],outdir='./SDO',nproc=1)
    aia94,aia193,hmi=groups
    c2=dm.get_LASCO_imgs(sT,eT,600,'C2',outdir='./LASCO/C2',need_calib=False)
```

**Multiprocessing:** `nproc=1` runs serially. Larger values use `multiprocess.Pool` across instrument/product combinations, not individual frames. The pool is scoped with `with Pool(nproc)`; returned Maps and each worker's scientific-Python imports can still consume substantial memory. Start with serial mode for debugging and choose worker count according to available memory. On macOS, call parallel workflows from a script under `if __name__=='__main__':`. The nested worker relies on `multiprocess` serialization; replacing the import with standard-library `multiprocessing.Pool` is not a drop-in change. A pool context does not guarantee cleanup if the parent is forcibly killed. Overlay's own SDO acquisition calls explicitly use `nproc=1`.

### Radio spectra and light curves

| Function | Purpose and outputs |
|---|---|
| `get_norp_LC(sT,eT,outdir='./',plot=True)` | Downloads NoRP FITS and derives microwave Stokes-I/V light curves, polarization, peak summaries and flare quantities. Writes event-subfolder data, CSV/Excel/pickle summaries and optional figures. Returns `None`. |
| `get_norh_imgs(sT,eT,outdir='./')` | Retrieves NoRH image products and saves quicklook plots. Returns `None`. |
| `get_RSTN(sT,eT,outdir='./',downdata=True,plot=True)` | Attempts numerical Culgoora/Learmonth, e-CALLISTO and RSTN sources. Australian raw retrieval targets intervals overlapping 22–10 UT and uses the start-date filename convention. Archive quicklooks are a fallback when numerical files are unavailable. Returns paths, plots and errors. `downdata=False` reports candidates. |
| `radspec_to_pickle(filname,outdir='./')` | Converts SRS/FITS spectrograms or e-CALLISTO DataFrame pickles into saved `.p` data and PNGs; returns both paths. DataFrame inputs retain that format. |
| `plot_culgoora_spec(filename,outdir='./')` | Decodes raw `SPECYYMMDD.bz2` records using recorded timestamps and four bands spanning 18–1800 MHz. Saves spectrum/metadata pickle and PNG, returning both paths. Native amplitude codes are not calibrated SFU. |
| `get_DHband_quicklookDS(sT,eT,outdir='./')` | Retrieves Wind/STEREO archive images and assembles available quicklook panels. These are visualization products, not numerical spectra. |
| `get_DHband_calibDS(sT,eT,outdir='./',sats=('sta','stb','wind'),Ch_frq=1.0,stats='median',vmin=None,vmax=None)` | Downloads multi-day SPDF CDF data; uses STEREO L3 SFU and the supplied Wind Z-channel gain conversion. Produces dynamic spectra and nearest-channel light curves, recording the actual channel frequency. `Ch_frq` is MHz; `vmin/vmax` are log10(SFU) display limits. Returns paths/errors per spacecraft. |

Numerical DH products are organised as `CDFs/`, `DSpickle/` and `PNGs/`, each with `STEREO_A/`, `STEREO_B/` and `Wind/` subdirectories. Existing CDFs/directories are reused.

`Matrix_ops.Reg(DS,R,C,r0,r1,c0,c1,reg_ax='both',stats='min')` takes a 2-D matrix and its row/column coordinates. `reg_ax='row'`, `'col'` or `'both'` selects which axes to regularise; `stats='min'`, `'median'` or `'mean'` selects the spacing statistic. It assigns nearest input samples to the regular grid, leaves sufficiently distant samples as NaNs and returns `(DSn,Rgs,Cgs)` with masked coordinate arrays. This is nearest-sample regridding, not smoothing or flux calibration. Use ordered axes with enough samples and nonzero spacing. It accompanies the Wind/STEREO workflow; currently the explicit `Reg` call regularises Wind time rows, while STEREO L3 data follow their own loading path.

```python
result=dm.get_DHband_calibDS('2014-11-03T06:00:00','2014-11-03T06:16:00',
    outdir='./radio',sats=('sta','wind'),Ch_frq=1.0,stats='median')
```

### Searches and high-energy utilities

| Function | Purpose and outputs |
|---|---|
| `MWAavail(sT,DT=30,outdir='./')` | Queries MWA observations over a duration in minutes, saves an observation-table pickle and returns success/failure. |
| `find_flare(sT,eT,flare_cat='hessi_flare_catalog.csv')` | Returns `(IDs,Tab)` for catalog flares fully contained in the interval, including boundary times, sorted by start time. Accepts CSV or a pickled Astropy Table. |
| `find_RHESSI(sT,eT,outdir='./',flare_cat=...,baseurl=...)` | Saves the selected flare catalog and matching RHESSI quicklook files, plus available daily orbit panels. Returns `(IDs,Tab)`; daily panels may extend beyond the event interval. |
| `plot_goes_derivative(sT,eT,outdir='./',satellite=None,savedata=True)` | Retrieves GOES 1–8 Angstrom flux and computes its derivative using actual sample intervals. `satellite=None` selects the highest available satellite number. Saves a two-panel PDF and optional data pickle. |
| `find_coobs(sT,eT,instrs,wvs,outdir='./',RHESSI_flarecat=...,downdata=False,dT=None,HMI_product=None)` | Searches matching instrument/wavelength observations; wavelength lists use Astropy units. Optionally downloads and prepares supported products. `HMI_product` filters magnetogram/continuum/dopplergram; `None` leaves HMI unfiltered here. Returns failed instrument/wavelength searches. |

## solar_overlay.py

This is the image-combination workhorse. Its public functions are `overlay` and `build_pfss`; rendering and projection details live in `solar_overlay_helper.py`.

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
| `image_dirs=None` | A FITS directory, directory list or dataset-to-directories dictionary. Searches recursively and identifies products from headers. Directory mode does not acquire missing images. |
| `local_files=None` | Explicit `{dataset:[FITS paths]}` alternative. Missing keys raise unless acquisition is explicitly enabled. |
| `download=False` | Enables missing-image acquisition only when set to `True`; uses `Multiwb_datamine`. Local raw LASCO calibration with `need_LASCO_calib=True` can independently try remote replacements. |
| `diffs=None,diff_types=None` | Per-layer difference flags and `'running'`/`'base'` choices. Missing flags are false; PFSS is never differenced. |
| `first_frame_intensity=False` | Running differences normally search for a pre-start predecessor and fall back to intensity if absent. `True` makes the first frame intensity without requiring a predecessor. Subsequent frames use the preceding selected exposure. |
| `base_diff_DT=None` | Base reference time or per-dataset dictionary; defaults to the event start. Before the event, selects the latest exposure at/before the reference; after the event, earliest at/after; inside the event, nearest exposure. Missing base references raise. |
| `reference_search=3600` | Maximum distance in seconds for external-reference/predecessor searches. |
| `tmap=None` | Target Map/FITS defining fixed FoV and WCS. Otherwise the bottom image provides the grid and reference epoch. |
| `alpha=0.6,npix=1024,dpi=150,figsize=7.0` | Upper-layer maximum opacity, longest display-grid dimension, figure resolution and square figure size in inches. Upper image opacity also depends on intensity. |
| `outdir='./',cache_dir='./solar_data'` | Output PNG directory and acquisition/reference cache. |
| `need_LASCO_calib=False` | Default uses preview preparation. `True` enables verified Level-1/fallback acquisition or optional local calibration. Inspect product provenance before photometry. |
| `pfss_dir=None,pfss_kwargs=None,N_Bloop=np.nan` | Reuse saved PFSS, or supply model-building settings. Directory mode requires saved PFSS when that layer is requested. NaN keeps all closed loops; an integer retains the highest visible loops; zero shows open lines only. |

Layer order is LASCO C3, LASCO C2, HMI, AIA/EIT image layers and PFSS. AIA bands preserve their supplied order. Labels report actual layer observation times; the corner timestamp is the actual bottom-image time. Cadence selection does not imply perfectly simultaneous exposures. Reference images are reprojected to the current WCS before subtraction; differential solar rotation is not corrected. Prepared AIA Level-1.5 inputs are required and are not calibrated twice. EIT must be official Level-1.

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

It returns model/trace information including `outdir` and saves `boundary.fits`, `metadata.json`, `field_lines.npz` and `model_3d.npz`. The 3-D archive contains Bphi/Btheta/Br in gauss with longitude, sine-latitude and radial axes. Plotting projects the co-rotating geometry to the target epoch, hides segments behind the disk and extends open lines radially to the display FoV. Red/blue indicate outward/inward open fields; white indicates closed loops. Dashed segments beyond the source surface are display continuations, not computed PFSS field.

## Scientific context

These tools were developed for event-based solar flare and space weather studies, especially for constructing catalogs connecting flare emission, radio bursts, CMEs and energetic particle events. They support Sun-as-a-star analyses and event visualization across radio, EUV, X-ray and coronagraph observations.

## Notes

- Research scripts may contain site-specific paths and event-table column settings. Inspect standalone script inputs before execution; their APIs may differ from the main modules.
- Archive quicklooks, JP2 images, DN/s previews and running differences have different scientific uses from verified radiometrically calibrated products. Check units, timestamps and provenance before quantitative analysis.
- Keep the raw downloads and helper modules when sharing or reusing results. Pickle files should only be loaded from trusted sources.
- The supplied main/helper files passed syntax checks for this update. Live archive retrieval, instrument calibration and full macOS process-lifecycle validation were not run.
