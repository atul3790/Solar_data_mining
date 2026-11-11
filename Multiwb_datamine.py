"""
Multiwb_datamine.py -- observation retrieval and instrument preparation.

Public functions
----------------
Image download functions automatically save standard instrument-color PNGs.
FITS/JP2 data, calibration and return values are unchanged by PNG generation.

All event times are UTC. Image cadence dT is in seconds, except get_DHband_DS,
where dT is an event timestamp. outdir/outfold selects the output directory.

MWAavail(sT, DT=30, outfold='./')
    Input: start time sT and duration DT in minutes.
    Output: observation-table pickle in outfold; returns True/False for success.

get_norp_LC(sT, eT, outdir='./', plot=True)
    Input: event start/end times; plot enables light-curve figures.
    Output: downloaded NoRP FITS, flux/polarization CSVs, peak-flux summary and
    optional plots in an event subfolder. Returns None.

get_norh_imgs(sT, eT, outdir='./')
    Input: event start/end times and destination directory.
    Output: downloaded NoRH image products and quick-look plots. Returns None.

find_flare(sT, eT, flare_cat='hessi_flare_catalog.csv')
    Input: event interval and RHESSI catalogue (CSV or pickled Astropy Table).
    Output: (IDs, Tab), the IDs and catalogue rows of fully contained flares.

find_RHESSI(sT, eT, outdir='./', flare_cat='hessi_flare_catalog.csv', baseurl=...)
    Input: interval, RHESSI catalogue and optional archive root.
    Output: downloaded archive products; returns selected flare (IDs, Tab).

find_coobs(sT, eT, instrs, wvs, outdir='./', RHESSI_flarecat=...,
           downdata=False, dT=None, HMI_product=None)
    Input: interval, instrument list and corresponding wavelength lists with
    Astropy units; optional cadence and download flag. HMI_product selects
    'magnetogram', 'continuum' or 'dopplergram'; None leaves it unfiltered.
    Output: search reports and optional downloads, including prepared AIA;
    returns a list of failed instrument/wavelength searches.

plot_goes_derivative(sT, eT, satellite=None, savedata=True, outdir='./')
    Input: interval, optional GOES number and data-saving flag. None chooses
    the highest available satellite number.
    Output: GOES 1-8 Angstrom flux/derivative PDF and optional pickle; None.

radspec_to_pickle(filname, outdir='./')
    Input: RSTN SRS/CALLISTO FITS filename or e-CALLISTO DataFrame pickle.
    Output: .p data and matching spectrogram PNG; returns both paths.
    DataFrame inputs retain their format; SRS/FITS produce a dictionary.

plot_culgoora_spec(filename, outdir='./')
    Input: Culgoora SPECYYMMDD.bz2 raw file and output directory.
    Output: decoded spectrum/metadata .p and spectrogram PNG; returns both
    paths. Uses recorded timestamps and four frequency bands, 18–1800 MHz.
    Amplitudes remain native instrument codes, not calibrated SFU.

get_RSTN(sT, eT, outdir='./', downdata=True, plot=True)
    Input: UTC interval, output directory, download and plotting flags.
    Checks Learmonth/Culgoora directly, e-CALLISTO and all RSTN stations
    through Fido, attempting every numerical source. Australian retrieval
    applies to intervals overlapping 22–10 UT, using sT's date in filenames.
    Output: numerical products, optional .p/PNG products and repository
    quicklooks only when no numerical files are found. Returns a dictionary
    of paths and errors. downdata=False reports candidates without downloads.

get_DHband_quicklookDS(sT, eT, outdir='./')
    Input: UTC interval and output directory.
    Output: Wind/STEREO archive quicklook images and combined plots when
    available. These are visualization products, not numerical spectra.

get_DHband_calibDS(sT, eT, outdir='./', sats=('sta','stb','wind'),
                  Ch_frq=7.0, stats='median', vmin=None, vmax=None)
    Input: UTC interval; spacecraft choices; light-curve frequency Ch_frq
    in MHz; Wind regularisation statistic; optional log10(SFU) color limits.
    Output: dynamic spectra and light curves as PNGs and .p files.
    PNGs/ and DSpickle/ each contain STEREO_A/, STEREO_B/ and Wind/;
    original downloads are retained under CDFs/ with the same subfolders.
    Uses STEREO L3 flux directly and the existing Wind Z-channel gain
    conversion. Saves the actual light-curve channel frequency.
    Returns output paths and errors per spacecraft.
    Requires spacepy and Matrix_ops.Reg; reuses existing directories/CDFs.

get_LASCO_imgs(sT, eT, dT, instr, outdir='./LASCO', need_calib=False, helioviewer=True)
    Input: interval, cadence, detector C2/C3
    need_calib=True: Level-1/MSB (outdir/maps), else Helioviewer (outdir/maps_HV),
    else Level-0.5 DN/s. need_calib=False: archive Level-0.5 DN/s (outdir/maps), and
    with helioviewer=True a JP2 for every bin in raw/helioviewer_<instr>/.
    Output: list of Maps; originals in outdir/raw.

calib_LASCO(images, instr=None, outdir='./LASCO_calibrated')
    Input: existing FITS path/Map or list/tuple; infer C2/C3 if instr is None.
    Tries SolarSoft -> downloaded Level-1 -> JP2 -> uncalibrated DN/s preview.
    Verified MSB passes unchanged. SolarSoft is optional and requires IDL.
    Output: Map or list matching input form, plus saved FITS/logs in outdir.

get_SDO_imgs(sT, eT, dT, instrument='AIA', wave=None, outdir='./SDO')
    Input: interval, cadence, AIA/HMI; AIA wave in Angstrom (94, 131, 171, 193,
    211, 304, 335, 1600, 1700); HMI wave=None selects LOS magnetograms.
    helioviewer=True saves JP2 previews for missing cadence bins under:
	outdir/<instrument>/<wave>/raw_SDO/helioviewer/

    Output: time-sorted Maps and finished FITS in outdir; original downloads
    in outdir/raw_SDO. Uses find_coobs; does not repeat AIA preparation.

get_EIT_imgs(sT, eT, dT, wave, outdir='./EIT', archive_url=...)
    Input: interval, cadence, wave 171/195/284/304 Angstrom and optional
    equivalent official Level-1 archive mirror.
    Output: calibrated Level-1 Maps and selected FITS in outdir; rejects raw EIT.

get_STEREO_imgs(sT,eT,dT,spacecraft='A',instruments=('COR2',),
	waves=None,outdir='./STEREO',nproc=4,polar=0,helioviewer=True)
	spacecraft: 'A' or 'B'.
	instruments: EUVI, COR1, COR2 or HI1; a string or a sequence.
	waves: one entry per instrument, e.g. [[171,195],None,None].
		waves=None selects all four EUVI bands.
	nproc: number of FITS download workers; 1 runs serially.
	polar: COR1 FITS polarizer angle, 0/120/240 degrees.
	COR2 always returns total-brightness images (POLAR=1001, d*c2 files).
	helioviewer: save original JP2 previews for missing cadence bins.
	HI1: route Level-1 -> Helioviewer -> Level-0.5, like COR2. Level-1 HI
		(hi_prep DN/s, MSB or S10) is only at UKSSDC behind a registered
		login, so local Level-1 files under STEREO_<sc>/HI1 (e.g. raw/level1)
		are used first; Helioviewer has no HI source (step reported and
		skipped); otherwise Level-0.5 science sums (s*h1*.fts, N_IMAGES>1)
		come from SSC/VSO into STEREO_<sc>/HI1/maps. Bins are
		max(dT, 2400 s), the HI1 science cadence.
	Returns one time-sorted list of FITS-based Maps per combination.
	JP2s are saved separately and are not returned as Maps.

load_map(source)
	Load any FITS/JP2 path or Map with a truthful WCS. Use it in every
	consumer (overlay, flow tracking). Archive FITS keep their CROTA/PC roll;
	Helioviewer pixels are north-up, so stale rolls on Helioviewer-derived
	files are zeroed (old value in OV_ROLL). Callers then only follow the WCS:
	rotate()/reproject when the WCS says so, never by source or level.
    
get_helioviewer_jp2(sT,eT,dT,instrument,wave=None,spacecraft='A',outdir)
	Save a Helioviewer JP2 for every cadence bin (LASCO/AIA/HMI/EUVI/COR1/COR2).

jp2_to_fits(jp2_paths,outdir,skip_times=())
	Convert JP2s to FITS (conventionally into a maps_HV folder): display
	values, roll 0 (old in OV_ROLL), observer kept, OV_HV/OV_SRC tags.

helioviewer=True in get_SDO_imgs/get_LASCO_imgs/get_STEREO_imgs saves a JP2
for every cadence bin in raw/helioviewer*/ alongside the science FITS.

Sources: FITS via Fido/VSO (SDAC, SSC data archive, JSOC), the LASCO and
EIT instrument archives, and Helioviewer JP2. Browse JPG/PNG quicklooks
(e.g. stereo-ssc .../browse/) carry no WCS and are never used as images.
The image functions save FITS Maps to user-chosen outdir subfolders for reuse by
solar_overlay.overlay(...). AIA from find_coobs and official EIT
Level-1 are already prepared. Local LASCO calibration requires SolarSoft/IDL;
verified downloaded Level-1 products do not. Internal helpers are in the
companion datamine_calib_helper.py; keep it beside this module.
This module does not import solar_overlay or run an observing workflow on import.
"""
__all__ = ['MWAavail', 'get_norp_LC', 'get_norh_imgs', 'find_flare', 'find_RHESSI', 'plot_goes_derivative', 'find_coobs', 'radspec_to_pickle', 'plot_culgoora_spec', 'get_RSTN', 'get_DHband_quicklookDS', 'get_DHband_calibDS', 'get_LASCO_imgs', 'calib_LASCO', 'get_SDO_imgs', 'get_EIT_imgs', 'get_STEREO_imgs', 'get_helioviewer_jp2', 'jp2_to_fits', 'load_map']
import re,urllib,os,shutil
import urllib.request
from urllib.parse import urljoin,urlsplit
import warnings
from pathlib import Path
from astropy.io import fits
from astropy.time import Time
import sunpy.map
from sunpy.net import Fido, attrs as a
import numpy as np,pickle,glob
from PIL import Image
import sunpy.timeseries as ts
from astropy import units as u
from astropy.table import Table
from radiospectra import net
#from radiospectra.sources.callisto import CallistoSpectrogram
from datetime import datetime as dt, timedelta,timezone
from bs4 import BeautifulSoup as bs
from radiospectra.spectrogram import Spectrogram
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from sunpy.visualization import colormaps as _sunpy_cmaps  # Register instrument tables
from aiapy.calibrate import register, update_pointing
from aiapy.calibrate.utils import get_pointing_table
from sunpy.map import Map
import pandas as pd
import ssl
from copy import deepcopy as cpy
import datamine_calib_helper as _helper
from datamine_calib_helper import load_map, jp2_to_fits
import bz2
from ecallisto_ng.plotting.plotting import plot_spectrogram_mpl
from Matrix_ops import Reg
from multiprocess import Pool

import certifi
ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())

def _save_image_previews(images,outdir,instrument=None,product=None):
	"""Save display PNGs with instrument colors; never alter input FITS/JP2 pixels."""
	folder=Path(outdir).expanduser().resolve()
	try:folder.mkdir(parents=True,exist_ok=True)
	except OSError as exc:
		warnings.warn(f'Could not create preview folder: {exc}')
		return []
	written=[]
	for image in images:
		try:
			path=Path(image) if isinstance(image,(str,Path)) else None
			m=sunpy.map.Map(image) if path is not None else image
			kind=str(instrument or m.detector).upper()
			jp2=(path is not None and path.suffix.lower() in ('.jp2','.j2k')) or bool(m.meta.get('ov_hv',False))
			data=np.asarray(m.data,dtype=float)
			if data.ndim!=2:raise ValueError('A two-dimensional image is required.')
			valid=np.isfinite(data)
			if not valid.any():continue
			values=data[valid]
			band=str(product) if product is not None else ''
			if kind in ('AIA','EUVI','EIT'):
				wave=int(round(m.wavelength.to_value(u.angstrom)))
				band=str(wave)
				cmap=plt.get_cmap(('sdoaia' if kind=='AIA' else 'euvi' if kind=='EUVI' else 'sohoeit')+band)
			elif kind=='LASCO':
				band=str(m.detector).upper()
				cmap=plt.get_cmap('soholasco'+band[-1])
			else:
				cmap=m.plot_settings.get('cmap','gray')
			magnetogram=kind=='HMI' and any(tag in (str(product)+' '+str(m.measurement)).lower()
				for tag in ('magnet','los_magnetic','gauss'))
			if magnetogram and not jp2:
				limit=max(float(np.percentile(abs(values),99)),np.finfo(float).eps)
				lo,hi=-limit,limit
			elif jp2 and values.min()>=0 and values.max()<=255:
				lo,hi=0.,255.
			else:
				lo,hi=np.percentile(values,[1,99.5])
			weight=np.clip((data-lo)/max(hi-lo,np.finfo(float).eps),0,1)
			if kind in ('AIA','EUVI','EIT'):
				cmap,norm=_helper._image_display_settings(m,kind,band,jp2)
				weight=norm(np.ma.masked_invalid(data))
			name=re.sub(r'[^A-Za-z0-9_.-]','_',kind+('_'+band if band else ''))
			output=folder/f'{name}_{"JP2" if jp2 else "FITS"}_{m.date.strftime("%Y%m%d_%H%M%S")}.png'
			plt.imsave(output,np.ma.masked_invalid(weight),origin='lower',cmap=cmap,vmin=0,vmax=1)
			written.append(str(output))
		except Exception as exc:
			warnings.warn(f'Could not save image preview: {exc}')
	return written

### MWA obs query for a few date times
def MWAavail(sT,DT=30,outdir='./'):
	'''
	sT is start time in '2014-04-11T00:03:00' format
	DT is observation duration in min eg. 30
	outfold is the final output folder where found observation details are saved
	
	'''
	eT=str(np.datetime64(sT)+np.timedelta64(DT*60,'s'))
	ST='%3A'.join(sT.split(':')[:-1])
	ET='%3A'.join(eT.split(':')[:-1])
	try:
		res=pd.read_html('http://ws.mwatelescope.org/metadata/find/?search=search&html=1&pagesize=100&dataquality=1&notdeleted=on&delete_id=&min_total_tiles=&max_total_tiles=&min_good_tiles=&max_bad_tiles=&projectid=G0002&groupid=&obsname=&creator=&mintime=&maxtime=&mintime_utc='+ST+'&maxtime_utc='+ET+'&minduration=&minra=&maxra=&mindec=&maxdec=&mingal_long=&maxgal_long=&mingal_lat=&maxgal_lat=&minel=&maxel=&minaz=&maxaz=&gridpoint=&minlst=&maxlst=&minsunel=&maxsunel=&minsunpd=&maxsunpd=&mode=&cenchan=&anychan=&freq_res=&int_time=&minfiles=&contigfreq=')
		taBl=res[0]
		pickle.dump(taBl,open(outdir+sT.replace(':','')+'-'+eT.replace(':','')+'_MWAobsIDs.p','wb'))  
		#pickle.dump(taBl,open('MWA_obsIDs/'+sT+'-'+eT+'_MWAobsIDs.p','wb'))
		print('MWA Obs found with ',len(taBl),' IDs saved to ',sT+'-'+eT+'_MWAobsIDs.p')
		return True
	except:
		return False
#################### Download and plot NoRP data 
#### Nobeyama Radio Polarimeter sun-as-a-star data
def get_norp_LC(sT,eT,outdir='./',plot=True):
	'''
	sT is the datetime string of the event 'YYYY-MM-DDThh:mm:ss'.
	eT is the end datetime string of the event 'YYYY-MM-DDThh:mm:ss'.
	outdir is the path to directory where to store the light curve and quick look image data files
	plot =True ensures the data will be plotted for the sT to eT period
	'''
	tfmt=lambda t: np.datetime_as_string(np.datetime64(t,'ms')+np.timedelta64(50,'ms'),unit='ms')[:-2]
	ratio=lambda v,i: np.divide(v,i,out=np.full(np.broadcast(v,i).shape,np.nan),where=np.isfinite(v)&np.isfinite(i)&(np.asarray(i)!=0))
	dtFmt=mdates.DateFormatter('%b-%d %H:%M')
	baseurl='https://solar.nro.nao.ac.jp/norp/fits/'
	sT=np.datetime64(sT)
	eT=np.datetime64(eT)
	if np.isnat(sT) or np.isnat(eT) or eT<=sT:
		raise ValueError('eT must be later than sT; both must be valid UTC times.')
	outdir=os.path.join(os.fspath(outdir),'')
	dates=np.arange(sT.astype('datetime64[D]')-np.timedelta64(1,'D'),
	 eT.astype('datetime64[D]')+np.timedelta64(2,'D'),np.timedelta64(1,'D'))
	chkdts=[str(i).replace('-','/') for i in dates]
	urls=[baseurl+i[:-2]+'norp'+i[2:].replace('/','')+'.fits.gz' for i in chkdts]
	## Expected frequencies
	subdir=str(sT+(eT-sT)/2).replace(':','')+'/'
	Flrtim=str(sT+(eT-sT)/2).replace(':','')
	print('###################################################\n')
	print('Searching for NoRP data in ',str(sT),' - ',str(eT))
	print('Flare time: ',Flrtim)
	print('###################################################\n')

	os.makedirs(outdir+subdir,exist_ok=True)
	#### Check completed outputs before downloads or FITS reads ####
	cache=outdir+subdir+'MaxFlux_'+tfmt(sT).replace(':','')+'-'+tfmt(eT).replace(':','')+'.csv'
	if os.path.isfile(cache):
		if plot==False:
			print('Analysis already completed for the event:',Flrtim)
			return
		Frqs=pd.read_csv(cache,usecols=['Frequency (GHz)'])['Frequency (GHz)'].to_numpy()
		getall=[outdir+subdir+'Flux_'+str(int(np.round(float(fq),0)))+'GHz_'+tfmt(sT).replace(':','').replace('-','')+'-'+tfmt(eT).replace(':','').replace('-','')+'.png' for fq in Frqs]
		if len(getall)>0 and all(os.path.isfile(i) for i in getall):
			print('Analysis and all plots already present! Skipping FITS reads and downloads.')
			return
 ###############################################################
	badurls=[] 
	for url in urls:
	 ####### Checking existing downloaded material and analysis
		if os.path.isfile(outdir+subdir+url.split('/')[-1]):
			print(url,' already analysed!')
			if plot==True and os.path.isfile(cache) and len(getall)>0 and all(os.path.isfile(i) for i in getall):
				print('Plots also present!! So continuing to next url!')
				continue
  ########################################################
		print('Data searching in : ',url,'\nFor flare time: ',Flrtim)
		fil=outdir+subdir+url.split('/')[-1]
		if not os.path.isfile(fil):
			fil=outdir+url.split('/')[-1]
		if not os.path.isfile(fil):
			try:
				urllib.request.urlretrieve(url,fil)
			except:
				print('No data or bad url')
				badurls+=[url]
				continue
		head=fits.getheader(fil,0)
		DelT=float(fits.getheader(fil)['XPOSURE'])
		beg=np.datetime64(head['DATE-BEG'])
		End=np.datetime64(head['DATE-END'])
		if beg>eT: # If Obs begin date time > End date time req. by user then break the loop 
			break
		if End<sT: # If Obs end date time < Start epoch requested then check the next file
			continue 
		if beg<=eT and End>=sT:
			print('Data found at ',url,';   Flare:',Flrtim)
			data=fits.getdata(fil)
			Time=data['Time'][0]  
			tmp=np.array(list(head.keys()))
			tfqs=tmp[[True if 'FREQ' in i else False for i in tmp]]
			Frqs=[str(np.round(float(head[i])*10**-9,1)) for i in tfqs]
			Alltags=['_'+str(int(np.round(float(i),0)))+'GHz' for i in Frqs]
			frq_dict=dict(zip(Frqs,Alltags))

			idlst = dt(1979, 1, 1, hour=0, minute=0, second=0, microsecond=0, tzinfo=timezone.utc)
			times=[dt.fromtimestamp(idlst.timestamp()+i,timezone.utc).replace(tzinfo=None) for i in Time]
			times=np.array(times,dtype='datetime64[us]')
			if not np.any((times>=sT)&(times<=eT)):
				continue
			t1=sT if sT>times[0] else times[0]
			t2=eT if times[-1]>eT else times[-1]

			if os.path.isfile(outdir+subdir+'MaxFlux_'+tfmt(t1).replace(':','')+'-'+tfmt(t2).replace(':','')+'.csv'):
				print('Analysis completed for the event: ',Flrtim)
				if plot==True:
					getall=[outdir+subdir+'Flux'+tag+'_'+tfmt(t1).replace(':','').replace('-','')+'-'+tfmt(t2).replace(':','').replace('-','')+'.png' for tag in Alltags]
					if all(os.path.isfile(i) for i in getall):
						print('All plots are in place for flare within',sT,' - ',eT,' So doing next time slot!! ')
						continue
					print('Not all frequency LC made; recomputing all frequencies to retain the full summary.')
				else:
					continue

			l1,l2=np.where(times>=t1)[0][0],np.where(times<=t2)[0][-1]
			IFlxs=[]
			VFlxs=[]
			VMxt=[]
			nVm=[]
			nIm=[]
			IMxt=[]
			VI=[]
			Imed=[]
			Vmed=[]
			Durs=[]
			TotI=[]
			TotV=[]
			for fq,tag in frq_dict.items():
				Flx_I=data['CalI'+tag][0]*data['Dval'+tag][0]
				Flx_V=data['CalV'+tag][0]*data['Dval'+tag][0]
				Flx_P=np.divide(Flx_V,Flx_I,out=np.full(Flx_I.shape,np.nan),where=np.isfinite(Flx_I)&np.isfinite(Flx_V)&(Flx_I!=0))
    
				## Flux estimation
				Im=np.nanmax(Flx_I[l1:l2+1])
				Vm=np.nanmax(Flx_V[l1:l2+1])
				Imd=np.nanmedian(Flx_I[l1:l2+1])
				Vmd=np.nanmedian(Flx_V[l1:l2+1])
				IFlxs+=[Im]
				VFlxs+=[Vm] 
				VI+=[float(ratio(Vm,Im))]
				Imed+=[Imd]
				Vmed+=[Vmd]
				Ivl=[np.nan,np.nan]
				Iy1,Iy2,y1,y2,Idur,TI,TV=np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan 

				tog_I=True # Toggle STOKES I plotting on 
				tog_V=True # toggle STOKES V plotting on . Will be set False if data is bad.
				if not np.isfinite(Im):
					tog_I=False
				if not np.isfinite(Vm):
					tog_V=False
      
				if tog_I==False and tog_V==False:
					print('Bad data in Freq: ',tag[1:],'\nContinuing with next Freq...!')
					VMxt+=['---']
					IMxt+=['---']
					nVm+=[np.nan]
					nIm+=[np.nan]
					Durs+=[np.nan]
					TotI+=[np.nan]
					TotV+=[np.nan]
					pickle.dump({'Time':times,'Flux_I':Flx_I,'Flux_V':Flx_V,'Flux V/I':Flx_P,'Start time':t1,'End time':t2,'I Flux range':[Iy1,Iy2],'V Flux range':[y1,y2]},open(outdir+subdir+'FluxLC'+tag+'_'+tfmt(t1).replace(':','')+'-'+tfmt(t2).replace(':','')+'.p','wb'))
					continue
     
    # Ylims in plot and max locations.
				if tog_I==True:
					Iy1,Iy2=np.round(np.nanmin(Flx_I[l1:l2+1])*0.9)-1,np.round(Im*1.1)+1
					Ipls=np.where(Flx_I[l1:l2+1]==Im)[0]
					IMxt+=[tfmt(times[l1:l2+1][Ipls[0]]).replace('T',' ')]
					nIm+=[len(Ipls)]
					## Identify flare period
					#Flocs=np.where((Flx_I[l1:l2+1]>Im/2) & (Flx_I[l1:l2+1]>1.5*Imd))[0] # Definition of flare duration : Flux is atleast above 1.5 x median Flux & flux > Half Max
					Flocs=np.where(((Flx_I[l1:l2+1]-Imd)>(Im-Imd)/2) & (Flx_I[l1:l2+1]>1.15*Imd))[0] # Definition of flare duration : Flux is atleast above 1.5 x median Flux & flux > Half Max
					if len(Flocs)>0:
						FI=times[l1:l2+1][Flocs]
						Ivl=[FI[0],FI[-1]] # Will mark vertical period for flare flux estimation
						Idur=(FI[-1]-FI[0])/np.timedelta64(1,'m')
						TI=np.nansum(Flx_I[l1:l2+1][Flocs])*DelT # Estimating total flare I flux
						if tog_V==True:
							TV=np.nansum(Flx_V[l1:l2+1][Flocs])*DelT # Estimating total flare V flux										
				else:
					IMxt+=['---']
					nIm+=[np.nan]
				if tog_V==True:
					y1,y2=np.nanmin(Flx_V[l1:l2+1]),Vm
					pad=max((y2-y1)*0.1,1.)
					y1,y2=y1-pad,y2+pad
					Vpls=np.where(Flx_V[l1:l2+1]==Vm)[0]
					VMxt+=[tfmt(times[l1:l2+1][Vpls[0]]).replace('T',' ')]
					nVm+=[len(Vpls)]
				else:
					VMxt+=['---']
					nVm+=[np.nan]
				TotI+=[TI]
				TotV+=[TV]
         
				### Plotting section
				if plot==True:
					fig=plt.figure(figsize=(8,6),dpi=90)
					if tog_I==True:
						ax=fig.add_subplot(111)
						plt.plot(times,Flx_I,'k',alpha=0.8)
						plt.ylabel('STOKES I Flux density (SFU)',size=16)
						ax.set_ylim([Iy1,Iy2])
						plt.gca().xaxis.set_major_formatter(dtFmt)
						plt.xlabel('Time [UT] +'+str(sT).split('-')[0],size=16)
						plt.xticks(rotation=40,size=14)
						plt.yticks(size=14)

					if tog_I==True and tog_V==True:
						ax1=ax.twinx()
					elif tog_V==True:
						ax1=fig.add_subplot(111)
					if tog_V==True:
						plt.plot(times,Flx_V,'r-',alpha=0.6)
						if tog_I==False:
							plt.xlabel('Time [UT] +'+str(sT).split('-')[0],size=16)      
							plt.gca().xaxis.set_major_formatter(dtFmt)
							plt.xticks(rotation=40,size=14)
						plt.ylabel('STOKES V Flux density (SFU)',size=16,color='r')
						ax1.set_ylim([y1,y2])
						plt.yticks(size=14,color='r')
					for lt in Ivl:
						if np.isfinite(lt):
							plt.axvline(lt,color='b',linewidth=2,linestyle='--')
					plt.xlim([t1,t2])
					plt.title(tag[1:].split('GHz')[0]+' GHz Light curve',size=16)
					plt.tight_layout()
					plt.savefig(outdir+subdir+'Flux'+tag+'_'+tfmt(t1).replace(':','').replace('-','')+'-'+tfmt(t2).replace(':','').replace('-','')+'.png')
					plt.close()
    
				pickle.dump({'Time':times,'Flux_I':Flx_I,'Flux_V':Flx_V,'Flux V/I':Flx_P,'Start time':t1,'End time':t2,'I Flux range':[Iy1,Iy2],'V Flux range':[y1,y2]},open(outdir+subdir+'FluxLC'+tag+'_'+tfmt(t1).replace(':','')+'-'+tfmt(t2).replace(':','')+'.p','wb'))   

			if os.path.abspath(fil)!=os.path.abspath(outdir+subdir+os.path.basename(fil)):
				os.replace(fil,outdir+subdir+os.path.basename(fil))
			IFlxs=np.array(IFlxs)
			VFlxs=np.array(VFlxs)
			DF=pd.DataFrame({'Frequency (GHz)':Frqs,'Median_I_flux (SFU)':Imed,'Median_V_flux (SFU)':Vmed,'Max_I_flux (SFU)':IFlxs,'Max_V_flux (SFU)':VFlxs,'Max Pol':ratio(VFlxs,IFlxs),'Max_to_median_I (SFU)':ratio(IFlxs,np.array(Imed)),'Max_to_median_V (SFU)':ratio(VFlxs,np.array(Vmed)),'Max_I_time (UT)':IMxt,'Max_V_time (UT)':VMxt,'Max_I_repetition':nIm,'Max_V_repetition':nVm,'Flare_Energy_I (SFU s)':TotI,'Flare_Energy_V (SFU s)':TotV,'Mean pol':ratio(np.array(TotV),np.array(TotI))})
			DF=DF.astype(str)
			DF.to_excel(outdir+subdir+'MaxFlux_'+tfmt(t1).replace(':','')+'-'+tfmt(t2).replace(':','')+'.xlsx',index=False)
			DF.to_csv(outdir+subdir+'MaxFlux_'+tfmt(t1).replace(':','')+'-'+tfmt(t2).replace(':','')+'.csv',index=False)
			pickle.dump(DF,open(outdir+subdir+'MaxFlux_'+tfmt(t1).replace(':','')+'-'+tfmt(t2).replace(':','')+'.p','wb'))
	if len(badurls)>0:
		pickle.dump(badurls,open(outdir+subdir+'Bad_urls_'+str(Flrtim)+'.p','wb'))
############### Download NoRH images ####################
def get_norh_imgs(sT,eT,outdir='./'):
	'''
	INPUT
	-----
	sT is the datetime string of the event 'YYYY-MM-DDThh:mm:ss'.
	eT is the end datetime string of the event 'YYYY-MM-DDThh:mm:ss'.
	outdir is where to save the images
		
	'''
	if np.datetime64(eT)<=np.datetime64(sT):
		raise ValueError('eT must be later than sT.')
	outdir=os.path.abspath(outdir)
	os.makedirs(outdir,exist_ok=True)
	cwd=os.getcwd()
	os.chdir(outdir)
	try:
		baseurl='https://solar.nro.nao.ac.jp/norh/fits/3min/'
		if np.datetime64(sT)<np.datetime64(sT.split('T')[0]+'T06:26:05'):
			YR=sT.split('T')[0].replace('-','/')
		else:
			YR=eT.split('T')[0].replace('-','/')
   
		URL=baseurl+YR+'/'
		req=urllib.request.urlopen(URL)
		links=bs(req).find_all('a',href=True)
		dtims=[]
		lnks=[]
		for ln in links:
			if '.gz' not in ln.get('href'):
				continue
			if 'ifs' in ln.get('href'):
				continue 
			lnk=ln.get('href')
			lnks+=[lnk]
			Dtime=lnk.replace('ifa',sT[0:2]).split('.')[0]
  
			dtims+=[np.datetime64(dt.strptime(Dtime,'%Y%m%d_%H%M%S'))]
		dtims=np.array(dtims)
		if len(dtims)==0:
			print('No NoRH images found:',URL)
			return
		order=np.argsort(dtims)
		dtims=dtims[order]
		lnks=[lnks[i] for i in order]
		try:
			p1=np.where(dtims<np.datetime64(sT))[0][-1]
		except:
			p1=0 
		try:
			p2=np.where(dtims>np.datetime64(eT))[0][0]
		except:
			p2=len(lnks)-1
   
		for i in np.arange(p1,p2+1):
			for img in [lnks[i],lnks[i].replace('ifa','ifs')]:
				if os.path.isfile(img) and all(os.path.isfile(img.split('.')[0]+ext) for ext in ('.png','.pdf')):
					continue
    #### Validate cached FITS; retry incomplete or failed downloads ####
				for attempt in range(2):
					try:
						if attempt==0 and os.path.isfile(img):
							data,header=fits.getdata(img,header=True)
						else:
							urllib.request.urlretrieve(URL+img,img+'.part')
							data,header=fits.getdata(img+'.part',header=True)
							os.replace(img+'.part',img)
						break
					except Exception as err:
						if os.path.isfile(img+'.part'):
							os.remove(img+'.part')
						print('Could not read/download NoRH FITS:',img,';',err)
				else:
					continue
				try:
					plt.figure(figsize=(8.5,7),dpi=90)
					polz='STOKES I' if header['POLARIZ']=='r+l' else 'STOKES V'
					dx=header['CDELT1']
					dy=header['CDELT2']
					xref=header.get('CRPIX1',(data.shape[1]+1)/2)
					yref=header.get('CRPIX2',(data.shape[0]+1)/2)
					extent=[(header.get('CRVAL1',0)+(x-xref)*dx)/60 for x in (0.5,data.shape[1]+0.5)]
					extent+=[(header.get('CRVAL2',0)+(y-yref)*dy)/60 for y in (0.5,data.shape[0]+0.5)]
					if polz=='STOKES I':
						data=np.ma.masked_less_equal(data,1000)
						plt.imshow(np.ma.log10(data),aspect='auto',extent=extent,origin='lower',cmap='gist_heat')
						h=plt.colorbar()
						h.set_label(r'$\mathrm{log_{10}T_B\ (K)}$',size=17)
					else:
						vmax=np.nanmax(np.abs(data))
						if not np.isfinite(vmax) or vmax==0:
							vmax=1.
						data1=cpy(data)
						data1[data1!=0]=np.log10(np.abs(data1[data1!=0]))*data1[data1!=0]/np.abs(data1[data1!=0])
						plt.imshow(data1,aspect='auto',extent=extent,origin='lower',cmap='RdBu_r',vmin=-np.log10(vmax),vmax=np.log10(vmax))
						h=plt.colorbar()
						h.set_label(r'$\mathrm{T_B\ (K)}$',size=17)
					h.ax.tick_params(labelsize=16)
					plt.xlabel('X [arcmin]',size=16)
					plt.ylabel('Y [arcmin]',size=16)
					plt.xticks(size=16)
					plt.yticks(size=16)
					plt.title(header['DATE-OBS']+' '+header['TIME-OBS'].split('.')[0]+'  ('+polz+')',size=17)
					plt.tight_layout()
					plt.savefig(img.split('.')[0]+'.png')
					plt.savefig(img.split('.')[0]+'.pdf')
					plt.close()  
					print('Image saved :',outdir+'/'+img.split('.')[0]+'.png')
				except Exception as err:
					print('NoRH plotting failed:',outdir+'/'+img,';',err)
				finally:
					plt.close()
      
	finally:
		os.chdir(cwd)

############### Find the HESSI flare from flare cat for the period #############
def find_flare(sT,eT,flare_cat="hessi_flare_catalog.csv"):
	'''
	Return IDs and catalogue rows for flares fully contained in sT–eT,
	sorted by start time. Boundary times are included.

	sT, eT: UTC strings 'YYYY-MM-DDThh:mm:ss'.
	flare_cat: CSV or pickled Astropy Table.
	'''
	stc=dt.strptime(sT,'%Y-%m-%dT%H:%M:%S')
	etc=dt.strptime(eT,'%Y-%m-%dT%H:%M:%S')
	if etc<=stc:
		raise ValueError('eT must be later than sT.')

	ext=os.path.splitext(os.fspath(flare_cat))[1].lower()
	if ext in ('.p','.pkl','.pickle'):
		with open(flare_cat,'rb') as fil:
			flcat=pickle.load(fil)
	elif ext=='.csv':
		flcat=Table.read(flare_cat,format='ascii.csv')
	else:
		raise ValueError('flare_cat must be a CSV or pickle file.')

	#### Select candidate dates ####
	days=[(stc+timedelta(days=i)).strftime('%d-%b-%Y').lstrip('0') for i in range((etc.date()-stc.date()).days+1)]
	sTab=flcat[np.isin(np.char.lstrip(np.asarray(flcat['Date'],dtype=str),'0'),days)]

	#### Get start and end times, including midnight crossings ####
	tst=[]
	tend=[]
	for row in sTab:
		start=dt.strptime(str(row['Date'])+' '+str(row['Start time']),
		 '%d-%b-%Y %H:%M:%S')
		end=dt.strptime(str(row['Date'])+' '+str(row['End time']),
		 '%d-%b-%Y %H:%M:%S')
		if end<start:
			end+=timedelta(days=1)
		tst.append(start)
		tend.append(end)

	#### Select fully contained flares and sort chronologically ####
	keep=np.asarray([stc<=start<=etc for start in tst],dtype=bool)
	order=np.argsort(np.asarray(tst,dtype=object)[keep])
	Tab=sTab[keep][order]
	IDs=np.asarray(Tab['ID']).copy()

	return IDs,Tab

##### Download RHESSI summary data #####
	
def find_RHESSI(sT,eT,outdir='./',flare_cat="hessi_flare_catalog.csv", baseurl='https://data.darts.isas.jaxa.jp/pub/mirror/rhessi/'):
	'''
	Save the selected flare catalogue and matching RHESSI quicklook files.

	sT, eT: UTC strings 'YYYY-MM-DDThh:mm:ss'.
	outdir: output directory.
	flare_cat: CSV or pickled Astropy Table.
	baseurl: RHESSI mirror root.

	Also download available daily orbit A.png panels.
	These panels cover whole days, rather than precisely sT to eT.
	Returns IDs, Tab. Failed requests are printed.
	'''

	IDs,Tab=find_flare(sT,eT,flare_cat)
	outdir=os.path.abspath(outdir)
	imgdir=os.path.join(outdir,'RHESSI_QLimgs')
	os.makedirs(imgdir,exist_ok=True)

	tag=sT.replace('-','').replace(':','')+'_'+eT.replace('-','').replace(':','')
	Tab.write(os.path.join(outdir,'RHESSI_flares_'+tag+'.csv'),
	 format='ascii.csv',overwrite=True)
	print('Selected flares:',len(IDs))

	stc=dt.strptime(sT,'%Y-%m-%dT%H:%M:%S')
	etc=dt.strptime(eT,'%Y-%m-%dT%H:%M:%S')

	#### Find actual archive filenames for each day ####
	for i in range((etc.date()-stc.date()).days+1):
		day=stc+timedelta(days=i)
		fsl=day.strftime('%Y/%m/%d/')
		dayIDs=[str(v) for v in
		 Tab['ID'][np.char.lstrip(np.asarray(Tab['Date'],dtype=str),'0')==day.strftime('%d-%b-%Y').lstrip('0')]]

		listings=[('orbit',
		 'http://sprg.ssl.berkeley.edu/~shilaire/rhessi/osp/'+fsl)]
		if dayIDs:
			listings.append(('flare',
			 baseurl.rstrip('/')+'/metadata/qlook_image_plot/'+fsl))

		for kind,LNK in listings:
			try:
				with urllib.request.urlopen(LNK,timeout=30) as req:
					links=bs(req.read(),'html.parser').find_all('a',href=True)
			except OSError as err:
				print('Could not read archive:',LNK,'\n',err)
				continue

			seen=set()
			for ln in links:
				url=urljoin(LNK,ln['href'])
				name=os.path.basename(urlsplit(url).path)

				if kind=='orbit':
					match=(name.startswith(day.strftime('%Y%m%d')+'_')
					 and name.endswith('A.png'))
				else:
					match=any(re.search(r'(?<!\d)'+re.escape(ID)+r'(?!\d)',
					 name) for ID in dayIDs)

				if not match or not name or url in seen or ln['href'].endswith('/'):
					continue
				seen.add(url)
				path=os.path.join(imgdir,name)
				if os.path.isfile(path):
					continue

    #### Download to a temporary file before replacing output ####
				try:
					with urllib.request.urlopen(url,timeout=60) as req, \
      open(path+'.part','wb') as fil:
						shutil.copyfileobj(req,fil)
					os.replace(path+'.part',path)
					print('Saved:',path)
				except OSError as err:
					if os.path.exists(path+'.part'):
						os.remove(path+'.part')
					print('Download failed:',url,'\n',err)

	return IDs,Tab

#### Download and plot SXR lightcurve and its derivative

def plot_goes_derivative(sT, eT, outdir='./',satellite=None,savedata=True):
	'''
	Download and plot GOES 1–8 Å flux and its time derivative.
	sT, eT: UTC datetime strings or datetime objects.
	satellite: GOES number; None selects the highest available number.
	savedata: True or False. To save data in the outdir
	'''
	if isinstance(sT, str):
		sT = dt.fromisoformat(sT)
	if isinstance(eT, str):
		eT = dt.fromisoformat(eT)
	if eT <= sT:
		raise ValueError('eT must be later than sT.')

	nametag=sT.isoformat().replace('-','').replace(':','')+'-'+eT.isoformat().replace('-','').replace(':','')
	outdir=os.path.join(os.path.abspath(outdir),'')
	os.makedirs(outdir,exist_ok=True)
	if satellite is None and os.path.isfile(outdir+'GOES_data_deriv_'+nametag+'.pdf') and (not savedata or os.path.isfile(outdir+'GOES_data_deriv_'+nametag+'.p')):
		return
 #### Download data ####
	print('Fetching GOES X-ray data...')
	query = [a.Time(sT, eT), a.Instrument('XRS'), a.Resolution('flx1s')]
	if satellite is not None:
		query.append(a.goes.SatelliteNumber(int(satellite)))
	results = Fido.search(*query)

	if results.file_num == 0:
		raise ValueError('No GOES XRS data found for this interval.')

	records = results['xrs']
	if satellite is None:
		satellite = max(int(v) for v in records['SatelliteNumber'])
		records = records[np.asarray(records['SatelliteNumber'], dtype=int) == satellite]
	print('Using GOES-' + str(satellite))

	downloaded = Fido.fetch(records,path=outdir)
	if downloaded.errors:
		raise RuntimeError('Some GOES downloads failed: ' + str(downloaded.errors))
	if len(downloaded) == 0:
		raise ValueError('No GOES files downloaded.')

 #### Read and restrict to requested interval ####
	goes_ts = ts.TimeSeries(downloaded, source='XRS', concatenate=True)
	goes_ts = goes_ts.truncate(sT, eT)
	data = goes_ts.to_dataframe().sort_index()
	data = data.loc[~data.index.duplicated(keep='first')]

	if len(data) < 2:
		raise ValueError('At least two samples are needed within the interval.')

	times = data.index
	flux_values = data['xrsb'].to_numpy(dtype=float)
	flux_values *= goes_ts.units['xrsb'].to(u.W/u.m**2)
 
	#### Derivative using actual sampling intervals ####
	time_steps = (times[1:] - times[:-1]).total_seconds().to_numpy()
	derivative = np.diff(flux_values) / time_steps
	derivative_times = times[:-1] + (times[1:] - times[:-1]) / 2
 
	if savedata==True:
		pickle.dump({'Time':times,'Flux':flux_values,'(Time,Derivative)':(derivative_times,derivative)},open(outdir+'GOES_data_deriv_'+nametag+'.p','wb'))

 #### Plot ####
	fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
	ax[0].plot(times, np.where(flux_values > 0, flux_values, np.nan),
	 label='GOES-' + str(satellite) + ' 1–8 Å', color='tab:blue')
	ax[0].set_ylabel('Flux (W/m²)')
	ax[0].set_yscale('log')
	ax[0].legend()
	ax[0].grid(True, which='both')

	ax[1].plot(derivative_times, derivative, label='Derivative', color='tab:red')
	ax[1].set_ylabel('Flux derivative ('+r'$\mathrm{Wm^{-2}s^{-1}}$'+')')
	ax[1].set_xlabel('Time (UT) '+pd.Timestamp(times[0]).strftime('%Y-%m-%d'))
	ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
	ax[1].legend()
	ax[1].grid(True)

	fig.tight_layout()
	plt.savefig(outdir+'GOES_data_deriv_'+nametag+'.pdf')
	plt.close()
	return
 

# Example usage:
# plot_goes_derivative("2025-08-10 00:00", "2025-08-10 02:00")

######### Find otherwave coobs
def find_coobs(sT,eT,instrs,wvs,outdir='./',RHESSI_flarecat="hessi_flare_catalog.csv",downdata=False,dT=None,HMI_product=None,flat_outdir=False):
	'''
	sT is the start datetime string of the event 'YYYY-MM-DDThh:mm:ss'.
	eT is the end datetime string of the event 'YYYY-MM-DDThh:mm:ss'.
	instrs = ['EIS', 'AIA'] - list of instruments of choice. or just a string 'EIT'
	wvs=[[],[]] list of list of wavelengths as floats with astropy units for each instrument. leave as [] if all data is needed for an instrument
	outdir = path to directory where in case any data is found will be saved
	flat_outdir: True saves directly in outdir without instrument/wavelength folders;
		intended for single instrument-product calls from get_SDO_imgs.
	dT: optional cadence in seconds.
	HMI_product: 'magnetogram', 'continuum' or 'dopplergram'; None leaves HMI unfiltered.
	These select LOS magnetic field, continuum intensity or LOS velocity respectively.
	The archive names 'LOS_magnetic_field', 'intensity' and 'LOS_velocity' are also accepted.
	HMI_product applies only to HMI, including calls requesting several instruments.
	Example: find_coobs(sT,eT,['HMI'],[],outdir='./HMI/',downdata=True,HMI_product='continuum')
	Instrument-color PNGs are saved automatically under the download folder/PNGs.
	downdata=True will download observations found outdir

	'''
	outdir=os.path.join(os.path.abspath(outdir),'')
	os.makedirs(outdir,exist_ok=True)
	RHESSI_flarecat=os.path.abspath(RHESSI_flarecat)
	if 'str' in str(type(instrs)):
		instrs=[instrs]
	if 'str' in str(type(wvs)):
		wvs=[[wvs]]
	if len(instrs)==1 and len(wvs)>1:
		wvs=[wvs]
	RSTN_list=['ecallisto','holloman', 'learmonth', 'palehua', 'sagamore', 'sagamore hill', 'san vito', 'san-vito']
	bads=[]
	if wvs==[]:
		wvs=[[] for _ in range(len(instrs))] 
	if dT is not None and (not np.isfinite(dT) or dT <= 0):
		raise ValueError('dT must be finite positive seconds')
	query_options=([a.Sample(float(dT)*u.s)] if dT is not None else [])
	if HMI_product is not None:
		hmi_products={'magnetogram':'LOS_magnetic_field','continuum':'intensity',
		 'dopplergram':'LOS_velocity','los_magnetic_field':'LOS_magnetic_field',
		 'intensity':'intensity','los_velocity':'LOS_velocity'}
		if not isinstance(HMI_product,str) or HMI_product.lower() not in hmi_products:
			raise ValueError("HMI_product must be 'magnetogram', 'continuum', 'dopplergram' or their documented archive names.")
		hmi_observable=hmi_products[HMI_product.lower()]
	pointing_table=None
	cwd=os.getcwd()
  
	try:
		for i in range(len(instrs)):
			ins=instrs[i]
			search_options=query_options+([a.Physobs(hmi_observable)] if ins.lower()=='hmi' and HMI_product is not None else [])
			wv=wvs[i]
			if 'rhessi' in ins.lower():
				if downdata==True:
					find_RHESSI(sT,eT,outdir=outdir+ins,flare_cat=RHESSI_flarecat)
				else:
					print('Available RHESSI flares:',find_flare(sT,eT,flare_cat=RHESSI_flarecat)[0])
				continue
			if ins.lower() in RSTN_list:
				print(ins,' is an RSTN network instrument. Please use get_RSTN() from nexttime for more functionality and ease.')
				print('Calling get_RSTN() functionality! ')
				get_RSTN(sT,eT,instrs=[ins],outdir=outdir,downdata=downdata)
				continue
			os.chdir(cwd)
			destination=outdir if flat_outdir else outdir+ins
			os.makedirs(destination,exist_ok=True)
			os.chdir(destination)

			if len(wv)!=0:
				for wvj in wv:
					print('Searching for data: ',ins,'; Wavelength: ',wvj.value,str(wvj.unit))
					try:
						res=Fido.search(a.Time(sT.replace('T',' '),eT.replace('T',' ')),a.Instrument(ins),a.Wavelength(wvj),*search_options)
						print('Data found for ',ins,'\nWavelength: ',wvj)
					except:
						if 'rhessi' in ins.lower():
							res=find_RHESSI(sT,eT,outdir='./',flare_cat=RHESSI_flarecat,baseurl='https://data.darts.isas.jaxa.jp/pub/mirror/rhessi/')
						else:
							bads+=[(ins,wvj)]
						continue
     
					if res.file_num==0:
						print('No data found for ',sT,' - ',eT,' period\n',ins,' at Wavelength: ',wvs[i])
						continue
					print('Provider: ',res.keys())
					prv=res.keys()[0]

					#pickle.dump(Table(res[prv]),open(ins+'_'+str(wvj.value).replace('.','')+str(wvj.unit)+'_search_results.p','wb'))	
					if downdata==True:
						print('Downloading.....',ins,' data for ',wvj)
						print(res[prv])
						wave_folder='.' if flat_outdir else str(wvj.value)+str(wvj.unit)
						os.makedirs(wave_folder,exist_ok=True)
						bdir=os.getcwd()
						os.chdir(wave_folder)
						attempts=[]
						for client in sorted(res.keys(),key=lambda k:str(k).lower()!='vso'):
							table=res[client]
							attempts.append((str(client),table,{}))
							if str(client).lower()=='vso':
								# ROB can be an AIA mirror, even when Provider is JSOC.
								if ins.lower()=='aia':
									attempts.append(('VSO site ROB',table,{'site':'ROB'}))
								if 'Provider' in table.colnames:
									providers=np.asarray([str(p).strip().upper() for p in table['Provider']])
									for provider in sorted(set(providers),key=lambda p:(p!='ROB',p)):
										attempts.append(('VSO provider '+provider,table[providers==provider],{}))
						downloaded_files=[]
						for label,records,options in attempts:
							try:
								print('Trying download route:',label)
								dwnfils=Fido.fetch(records,path='./',**options)
								errors=getattr(dwnfils,'errors',[])
								for error in errors:
									print(error)
								for fil in dwnfils.data:
									if fil not in downloaded_files:
										downloaded_files.append(fil)
								if dwnfils.data and not errors:
									break
							except Exception as exc:
								print('Download failed:',label,exc)
						if not downloaded_files:
							print('No files downloaded from available providers for',ins)
							os.chdir(bdir)
							continue
						
						for fil in downloaded_files:
							if ins.lower()=='aia':
									try:
										aiah=fits.getheader(fil,1)
									except:
										aia,aiah=fits.getdata(fil,header=True)
										del aia
									tflg=0 
									try:
										if aiah['LVL_NUM']<1.5:
											tflg=1
									except:
										if 'lev1_' in fil or 'lev1.fits' in fil:
											tflg=1 
         ### Upgrade AIA images to level 1.5 
									if tflg==1:
										print('AIA data is in LEVEL 1. So upgrading to 1.5!!')
										AMap=Map(fil)
										if pointing_table is None:
											pointing_table=get_pointing_table(time_range=(sT,eT))
										AMap=update_pointing(AMap,pointing_table=pointing_table) # Doing rotation corrections
										AMap15=register(AMap) # Pixels fluxes are now corrected for direction offset issue and pixels rescaled to 0.6" and image y axes are now aligned with solar north-south axis
										Amap15n=AMap15/AMap15.exposure_time # Making map normalised with respect to exposure times to get values in DN/s/pixel

										data_lev1=fits.open(fil,mode='update')
										hdu=1 if len(data_lev1)>1 else 0
										data_lev1[hdu].data=Amap15n.data
										data_lev1[hdu].header.update(Amap15n.fits_header)
										data_lev1[hdu].header['LVL_NUM']=1.5
										data_lev1[hdu].header['BUNIT']='DN/s'
										data_lev1[hdu].header['EXPTIME']=1.00  
										data_lev1[hdu].header['SERIES']='aia.lev1.5'
										data_lev1[hdu].header['SEGMENT']='image_lev1.5.fits'
										data_lev1[hdu].header['PIXLUNIT']='DN/s'
										data_lev1.flush()
										data_lev1.close()
										os.replace(fil,fil.replace('lev1','lev15'))
										del AMap
										del AMap15 

							if ins.lower() in ('aia','hmi','euvi','eit','lasco'):
								preview_file=fil.replace('lev1','lev15') if ins.lower()=='aia' and os.path.isfile(fil.replace('lev1','lev15')) else fil
								_save_image_previews([preview_file],Path.cwd()/'PNGs',ins,HMI_product if ins.lower()=='hmi' else None)
							if ins.lower()=='norh':
								os.replace(fil,fil if fil.endswith('.fits') else fil+'.fits')
								get_norh_imgs(sT,eT,outdir=os.getcwd())
							if ins.lower()=='goes':
								dat=ts.TimeSeries(fil)
								Obs=dat.observatory.replace('-','')
								Tab=dat.to_table()
								tres=str(24/len(Tab)*3600).split('.')[0]+'s'
								rng=dat.time.value[0].replace('-','').replace(' ','T').replace(':','').split('.')[0]+'-'+dat.time.value[-1].replace('-','').replace(' ','T').replace(':','').split('.')[0]
								pickle.dump(Tab,open(Obs+'_'+rng+'_'+tres+'.p','wb'))
								Tab=Tab.to_pandas()
								Tab.to_csv(Obs+'_'+rng+'_'+tres+'.csv',index=False)    
						os.chdir(bdir)
				continue
			else:
				try:
					res=Fido.search(a.Time(sT.replace('T',' '),eT.replace('T',' ')),a.Instrument(ins),*search_options)
					print('Data found for ',ins,'\nWavelength: ',wv)
				except:
					bads+=[ins]
					continue
				if res.file_num==0:
					print('No data found for',ins)
					continue
				prv=res.keys()[0]

				#pickle.dump(res,open(ins+'_search_results.p','wb'))	
				if downdata==True:
					print('Downloading.....',ins,' data for ',wv)
					print(res[prv])
					attempts=[]
					for client in sorted(res.keys(),key=lambda k:str(k).lower()!='vso'):
						table=res[client]
						attempts.append((str(client),table,{}))
						if str(client).lower()=='vso':
							# ROB can be an AIA mirror, even when Provider is JSOC.
							if ins.lower()=='aia':
								attempts.append(('VSO site ROB',table,{'site':'ROB'}))
							if 'Provider' in table.colnames:
								providers=np.asarray([str(p).strip().upper() for p in table['Provider']])
								for provider in sorted(set(providers),key=lambda p:(p!='ROB',p)):
									attempts.append(('VSO provider '+provider,table[providers==provider],{}))
					downloaded_files=[]
					for label,records,options in attempts:
						try:
							print('Trying download route:',label)
							dwnfils=Fido.fetch(records,path='./',**options)
							errors=getattr(dwnfils,'errors',[])
							for error in errors:
								print(error)
							for fil in dwnfils.data:
								if fil not in downloaded_files:
									downloaded_files.append(fil)
							if dwnfils.data and not errors:
								break
						except Exception as exc:
							print('Download failed:',label,exc)
					if not downloaded_files:
						print('No files downloaded from available providers for',ins)
						continue
					
					for fil in downloaded_files:
						if ins.lower()=='aia':
								try:
									aiah=fits.getheader(fil,1)
								except:
									aia,aiah=fits.getdata(fil,header=True)
									del aia
								tflg=0 
								try:
									if aiah['LVL_NUM']<1.5:
										tflg=1
								except:
									if 'lev1_' in fil or 'lev1.fits' in fil:
										tflg=1 
        ### Upgrade AIA images to level 1.5 
								if tflg==1:
									print('AIA data is in LEVEL 1. So upgrading to 1.5!!')
									AMap=Map(fil)
									if pointing_table is None:
										pointing_table=get_pointing_table(time_range=(sT,eT))
									AMap=update_pointing(AMap,pointing_table=pointing_table) # Doing rotation corrections
									AMap15=register(AMap) # Pixels fluxes are now corrected for direction offset issue and pixels rescaled to 0.6" and image y axes are now aligned with solar north-south axis
									Amap15n=AMap15/AMap15.exposure_time # Making map normalised with respect to exposure times to get values in DN/s/pixel

									data_lev1=fits.open(fil,mode='update')
									hdu=1 if len(data_lev1)>1 else 0
									data_lev1[hdu].data=Amap15n.data
									data_lev1[hdu].header.update(Amap15n.fits_header)
									data_lev1[hdu].header['LVL_NUM']=1.5
									data_lev1[hdu].header['BUNIT']='DN/s'
									data_lev1[hdu].header['EXPTIME']=1.00  
									data_lev1[hdu].header['SERIES']='aia.lev1.5'
									data_lev1[hdu].header['SEGMENT']='image_lev1.5.fits'
									data_lev1[hdu].header['PIXLUNIT']='DN/s'
									data_lev1.flush()
									data_lev1.close()
									os.replace(fil,fil.replace('lev1','lev15'))
									del AMap
									del AMap15 

						if ins.lower() in ('aia','hmi','euvi','eit','lasco'):
							preview_file=fil.replace('lev1','lev15') if ins.lower()=='aia' and os.path.isfile(fil.replace('lev1','lev15')) else fil
							_save_image_previews([preview_file],Path.cwd()/'PNGs',ins,HMI_product if ins.lower()=='hmi' else None)
						if ins.lower()=='norh':
							os.replace(fil,fil if fil.endswith('.fits') else fil+'.fits')
							get_norh_imgs(sT,eT,outdir=os.getcwd())
						if ins.lower()=='goes':
							dat=ts.TimeSeries(fil)
							Obs=dat.observatory.replace('-','')
							Tab=dat.to_table()
							tres=str(24/len(Tab)*3600).split('.')[0]+'s'
							rng=dat.time.value[0].replace('-','').replace(' ','T').replace(':','').split('.')[0]+'-'+dat.time.value[-1].replace('-','').replace(' ','T').replace(':','').split('.')[0]
							pickle.dump(Tab,open(Obs+'_'+rng+'_'+tres+'.p','wb'))
							Tab.write(Obs+'_'+rng+'_'+tres+'.csv',format='ascii.csv')    
				else:
					print('Data available: ',res[prv])
	finally:
		os.chdir(cwd)
	return bads
 
########### Convert and plot radio spec #######################	
def radspec_to_pickle(filname,outdir='./'):
	'''
	Read SRS/FITS using Spectrogram; otherwise read an e-CALLISTO
	DataFrame pickle. Save .p and .png; return their paths.
	DataFrame inputs retain their DataFrame format.

	Requires module-level imports: os, pickle, numpy as np,
	pandas as pd, matplotlib.pyplot as plt, astropy.units as u,
	Spectrogram and plot_spectrogram_mpl.
	'''
	outdir=os.path.join(os.path.abspath(outdir),'')
	os.makedirs(outdir,exist_ok=True)

	try:
		spec=Spectrogram(filname)
	except Exception:
		df=pd.read_pickle(filname)
		if not isinstance(df,pd.DataFrame):
			raise TypeError('Pickle must contain an e-CALLISTO DataFrame.')
		if df.empty:
			raise ValueError('DataFrame is empty.')
		specfil=os.path.splitext(os.path.basename(filname))[0]+'.p'
		data=df
		spec=None
	else:
		stt=str(spec.start_time.value).split('.')[0].replace('-','').replace(':','')
		edt=str(spec.end_time.value).split('.')[0].replace('-','').replace(':','')
		specfil=spec.observatory+'_'+'_'.join([stt,edt])+'.p'
		data={'DS':spec.data,'frequencies':spec.frequencies,
			'times':spec.times,'start_time':spec.start_time,
			'end_time':spec.end_time,'observatory':spec.observatory,
			'detector':spec.detector}

	path=outdir+specfil
	if os.path.realpath(path)!=os.path.realpath(filname):
		with open(path,'wb') as fil:
			pickle.dump(data,fil)

	if spec is None:
		fig=plot_spectrogram_mpl(df,cbar_label='Flux density (uncalib. unit)')
		ax=fig.axes[0]
		indices=np.unique(np.linspace(
			0,len(df)-1,min(7,len(df)),dtype=int))
		ax.set_xticks(indices)
		ax.set_xticklabels(pd.DatetimeIndex(df.index)[indices].strftime('%H:%M'),
			rotation=0)
		date=pd.Timestamp(df.index[0]).strftime('%Y-%m-%d')
	else:
		image=spec.plot()
		ax=image.axes
		fig=ax.figure
		indices=np.unique(np.linspace(
			0,len(spec.times)-1,min(7,len(spec.times)),dtype=int))
		ticks=spec.times[indices]
		ax.set_xticks(ax.convert_xunits(ticks))
		ax.set_xticklabels(ticks.utc.strftime('%H:%M'),rotation=0)
		date=spec.start_time.utc.isot[:10]

		# Format frequencies in MHz using the plotted axis unit.
		freq_unit=ax.yaxis.get_units()
		if freq_unit is None:
			freq_unit=spec.frequencies.unit
		ax.yaxis.set_major_formatter(
			plt.FuncFormatter(lambda value,pos:
				f'{(value*u.Unit(freq_unit)).to_value(u.MHz):g}'))

		unit=str(getattr(spec.data,'unit','')) or 'uncalib. unit'
		fig.colorbar(image,ax=ax,label=f'Flux density ({unit})')

	ax.set_xlabel('Time (UT) '+date)
	ax.set_ylabel('Frequency (MHz)')
	try:
		fig.savefig(path+'ng',dpi=200,bbox_inches='tight')
	finally:
		plt.close(fig)
	return path,path+'ng'

### Plot Culgoora data from .bz2 file
def plot_culgoora_spec(filename,outdir='./'):
	'''
	Read Culgoora SPECYYMMDD.bz2 using actual record timestamps.
	Save a spectrum dictionary .p and spectrogram .png.
	Return (pickle_path, plot_path).

	Each record: 8-byte time header, four 8-byte band headers,
	and four sets of 501 amplitudes. Frequencies are in MHz.
	Color limits affect the plot only; saved amplitudes are unchanged.

	Requires module-level imports: os, bz2, pickle, numpy as np,
	pandas as pd, matplotlib.pyplot as plt, matplotlib.dates as mdates.
	'''
	os.makedirs(outdir,exist_ok=True)
	base_name=os.path.splitext(os.path.basename(filename))[0]

	with bz2.open(filename,'rb') as fil:
		raw=fil.read()
	if not raw or len(raw)%2044:
		raise ValueError('File must contain complete 2044-byte Culgoora records.')

	records=np.frombuffer(raw,dtype=np.uint8).reshape(-1,2044)
	headers=records[:,:8]
	band_bytes=records[:,8:40].copy()

	#### Determine byte order from the four frequency bands ####
	expected=np.array([[18,57],[57,180],[180,570],[570,1800]])
	for endian in ('>','<'):
		words=band_bytes.view(endian+'u2').reshape(-1,4,4)
		limits=words[:,:,:2]
		if np.array_equal(limits[0],expected):
			break
	else:
		raise ValueError('Unrecognized Culgoora frequency-band headers.')
	if not np.all(limits==limits[0]):
		raise ValueError('Frequency bands change within this file.')

	#### Read actual UTC timestamps ####
	years=headers[:,0].astype(int)
	years=np.where(years>=70,years+1900,years+2000)
	times=pd.DatetimeIndex(pd.to_datetime({
		'year':years,'month':headers[:,1],'day':headers[:,2],
		'hour':headers[:,3],'minute':headers[:,4],'second':headers[:,5]
	},errors='raise',utc=True))
	order=np.argsort(times.asi8,kind='stable')
	times=times[order]
	if times.has_duplicates:
		raise ValueError('Duplicate timestamps prevent unambiguous plotting.')

	#### Decode amplitudes and construct the frequency array ####
	frequencies=np.concatenate([
		np.linspace(low,high,501,endpoint=False)
		for low,high in limits[0]
	])
	amplitudes=records[:,40:].astype(np.uint16)
	amplitudes=np.where(amplitudes&128,
		256+(amplitudes&127),amplitudes)
	ds_data=amplitudes[order].T

	data={
		'DS':ds_data,'frequencies':frequencies,
		'times':np.asarray((times-times[0]).total_seconds()),
		'start_time':times[0],'end_time':times[-1],
		'observatory':'CULGOORA','detector':'RAW_BOM_SPEC',
		'frequency_unit':'MHz','time_unit':'s',
		'intensity_unit':'decoded amplitude code',
		'record_headers':headers[order],
		'band_header_bytes':band_bytes[order],
		'timestamps':times
	}
	pickle_path=os.path.join(outdir,base_name+'.p')
	with open(pickle_path,'wb') as fil:
		pickle.dump(data,fil)

	#### Build time-bin edges ####
	x=mdates.date2num(times.to_pydatetime())
	if len(x)>1:
		mid=(x[:-1]+x[1:])/2
		xedges=np.r_[x[0]-(mid[0]-x[0]),mid,
			x[-1]+(x[-1]-mid[-1])]
	else:
		xedges=np.array([x[0]-0.5/86400,x[0]+0.5/86400])

	#### Match the default Learmonth figure aspect ratio ####
	fig,ax=plt.subplots()
	vmin,vmax=np.nanpercentile(ds_data,[5,95])
	if vmax<=vmin:
		vmax=vmin+1

	for band,(low,high) in enumerate(limits[0]):
		values=ds_data[band*501:(band+1)*501]
		im=ax.pcolormesh(xedges,np.linspace(low,high,502),values,
			shading='flat',cmap='jet',vmin=vmin,vmax=vmax,
			rasterized=True)

	#### Leave long gaps blank ####
	if len(x)>2:
		gaps=np.diff(x)
		cadence=np.median(gaps)
		for i in np.where(gaps>3*cadence)[0]:
			ax.axvspan(x[i]+cadence/2,x[i+1]-cadence/2,
				color='white',linewidth=0)

	date=times[0].strftime('%Y-%m-%d')
	if times[-1].date()!=times[0].date():
		date+=' — '+times[-1].strftime('%Y-%m-%d')

	ax.set_title('Culgoora spectrograph')
	ax.set_xlabel('Time (UT) '+date)
	ax.set_ylabel('Frequency (MHz)')
	ax.set_yscale('log')
	ax.set_ylim(18,1800)
	ax.set_xlim(xedges[0],xedges[-1])
	ax.set_yticks([18,30,57,100,180,300,570,1000,1800])
	ax.set_yticklabels(['18','30','57','100','180','300','570','1000','1800'])
	ax.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=5,maxticks=9))
	ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
	fig.colorbar(im,ax=ax,label='Flux density (uncalib. unit)')
	fig.tight_layout()

	plot_path=os.path.join(outdir,base_name+'_spectrogram.png')
	try:
		fig.savefig(plot_path,dpi=150,bbox_inches='tight')
	finally:
		plt.close(fig)
	return pickle_path,plot_path


###### Get RSTN spectrograph selected dates ####
def get_RSTN(sT,eT,outdir='./',downdata=True,plot=True):
	'''
	Get metric dynamic spectra and optionally plot downloaded data.

	Check Culgoora and Learmonth directly, then e-CALLISTO and all
	RSTN stations through Fido. Attempt every numerical source,
	even when another source succeeds.

	sT/eT: UTC strings 'YYYY-MM-DDThh:mm:ss'.
	downdata=False: report candidate URLs/search results only.
	plot=True: save plots using radspec_to_pickle/plot_culgoora_spec.

	Australian raw retrieval applies to intervals overlapping 22–10 UT,
	using sT's date in filenames. Download repository quicklooks only
	when no numerical source supplies files.

	Return a dictionary of downloaded paths, plots and errors.
	Requires module-level imports and companion plotting functions.

	'''
	start=dt.strptime(sT,'%Y-%m-%dT%H:%M:%S')
	end=dt.strptime(eT,'%Y-%m-%dT%H:%M:%S')
	if end<=start:
		raise ValueError('eT must be later than sT.')

	root=Path(outdir).expanduser().resolve()
	result={'ecallisto':[],'australian_raw':[],'rstn':[],
		'quicklooks':[],'plots':[],'errors':[]}

	#### Direct Australian raw products: intervals overlapping 22–10 UT ####
	day=start.replace(hour=0,minute=0,second=0)
	if start.hour<10 or start.hour>=22 or end>day+timedelta(hours=22):
		for station,name in (
			('learmonth',f'LM{start:%y%m%d}.srs'),
			('culgoora',f'SPEC{start:%y%m%d}.bz2')
		):
			url=f'https://downloads.sws.bom.gov.au/wdc/wdc_spec/data/{station}/raw/{start:%y}/{name}'
			if not downdata:
				print(url)
				result['australian_raw'].append(url)
				continue

			folder=root/station
			path=folder/'raw'/name
			try:
				path.parent.mkdir(parents=True,exist_ok=True)
				if not path.is_file():
					urllib.request.urlretrieve(url,path)
				if path.stat().st_size==0:
					raise ValueError('Downloaded file is empty.')
				result['australian_raw'].append(str(path))
				print('Raw data:',path)
			except Exception as exc:
				path.unlink(missing_ok=True)
				message=f'{station} raw download failed: {exc}'
				print(message)
				result['errors'].append(message)
				continue

			if plot:
				try:
					if station=='culgoora':
						_,png=plot_culgoora_spec(path,outdir=folder)
					else:
						_,png=radspec_to_pickle(path,outdir=folder)
					result['plots'].append(str(png))
				except Exception as exc:
					message=f'Plot failed for {path}: {exc}'
					print(message)
					result['errors'].append(message)

	#### First: e-CALLISTO, queried one UTC day at a time ####
	try:
		day=start.replace(hour=0,minute=0,second=0)
		while day<end:
			next_day=day+timedelta(days=1)
			lo=max(start,day)
			hi=min(end,next_day)
			try:
				if not downdata:
					urls=get_remote_files_url(
						day,next_day-timedelta(microseconds=1),
						callisto or None)
					for url in urls:
						match=re.search(r'_(\d{8}_\d{6})_',url)
						if match is None:
							continue
						t=dt.strptime(match.group(1),'%Y%m%d_%H%M%S')
						if t<=hi and t+timedelta(minutes=15)>=lo:
							if url not in result['ecallisto']:
								result['ecallisto'].append(url)
								print('e-CALLISTO:',url)
				else:
					data=get_ecallisto_data(
						lo,hi,callisto or None,verbose=True)
					if isinstance(data,pd.DataFrame):
						data={callisto or 'ecallisto':data}
					for station,df in (data or {}).items():
						if df is None or df.empty:
							continue
						df=df.copy()
						df.index=pd.to_datetime(df.index,utc=True)
						df=df.loc[
							(df.index>=pd.Timestamp(lo,tz='UTC')) &
							(df.index<=pd.Timestamp(hi,tz='UTC'))]
						df=df.dropna(how='all').sort_index()
						if df.empty:
							continue

						folder=root/'ecallisto'
						folder.mkdir(parents=True,exist_ok=True)
						name=re.sub(r'[^A-Za-z0-9_.-]','_',str(station))
						path=folder/f'{name}_{lo:%Y%m%dT%H%M%S}_{hi:%Y%m%dT%H%M%S}.p'
						df.to_pickle(path)
						result['ecallisto'].append(str(path))
						print(station,':',df.index.min(),'--',
							df.index.max(),'; channels:',len(df.columns))

						if plot:
							try:
								_,png=radspec_to_pickle(path,outdir=folder)
								result['plots'].append(str(png))
							except Exception as exc:
								message=f'Plot failed for {path}: {exc}'
								print(message)
								result['errors'].append(message)
			except Exception as exc:
				message=f'e-CALLISTO {lo} -- {hi}: {exc}'
				print(message)
				result['errors'].append(message)
			day=next_day
	except Exception as exc:
		message=f'e-CALLISTO unavailable: {exc}'
		print(message)
		result['errors'].append(message)

	#if result['ecallisto']:
	#	return result

	#### Second: all available numerical RSTN files through Fido ####
	try:
		query=Fido.search(a.Time(sT,eT),a.Instrument('RSTN'))
		print(query)
		if query.file_num==0:
			result['errors'].append('RSTN: no numerical files returned.')
		elif not downdata:
			result['rstn'].append(query)
		else:
			folder=root/'RSTN'
			folder.mkdir(parents=True,exist_ok=True)
			downloaded=Fido.fetch(query,path=str(folder/'{file}'))
			result['rstn'].extend(str(p) for p in downloaded)
			for error in getattr(downloaded,'errors',[]):
				result['errors'].append(str(error))

			if plot:
				for filename in downloaded:
					try:
						if Path(filename).name.upper().startswith('SPEC'):
							_,png=plot_culgoora_spec(filename,outdir=folder)
						else:
							_,png=radspec_to_pickle(filename,outdir=folder)
						result['plots'].append(str(png))
					except Exception as exc:
						message=f'Plot failed for {filename}: {exc}'
						print(message)
						result['errors'].append(message)
	except Exception as exc:
		message=f'RSTN/Fido: {type(exc).__name__}: {exc}'
		print(message)
		result['errors'].append(message)

	#if result['rstn'] or result['australian_raw']:
	#	return result

	#### Last: direct repository quicklook images ####
	result['quicklooks']=_helper._radio_quicklook_fallback(
		sT,eT,root,downdata)
	return result

#####################################################
## Get CDAW WIND/WAVES and STEREO A & B DS #########
def get_DHband_quicklookDS(sT,eT,outdir='./'):
	'''
	sT/eT is the datetime string of the event period 'YYYY-MM-DDThh:mm:ss'.
	delta_t: Time in minutes. Dynamic spectrum png will have start time from dT-delta_t. 
	outdir is where plot and daily spectrograph will be saved.
 	Returns combined 3 panel plot of quick look png images from Wind and STEREO spacecraft
	'''
	######## INPUT #############
	WAVES_png_url='https://cdaw.gsfc.nasa.gov/images/wind/waves/' # Webpage to get specific image png for the date of interest. (e.g., YYYY/MM/DD/YYYYMMDD_hhmmss_windwaves.png)
	WAVES_png_url_latest='https://cdaw.gsfc.nasa.gov/images/wind/waves_h2M/' # Webpage where post 2021 data resides for wind waves.
	SWAVES_png_url='https://cdaw.gsfc.nasa.gov/images/stereo/swaves/'
	SWAVES_newpng_url='https://cdaw.gsfc.nasa.gov/images/stereo/swaves_new/'
	##############
	dT=str(np.datetime64(sT)+(np.datetime64(eT)-np.datetime64(sT))/2)
	delta_t=int(np.ceil((np.datetime64(eT)-np.datetime64(sT)).astype('float')/120.))
	outdir=os.path.join(os.fspath(outdir),'')
	if outdir[-1]!='/':
		outdir+='/'
	dtup=dT.replace('-','/').split('T')
	Wu= WAVES_png_url_latest if np.datetime64(dT)>np.datetime64('2021-12-12T23:59:00') else WAVES_png_url
	SWu=SWAVES_png_url if np.datetime64(dT)<np.datetime64('2019-01-01T00:00:00') else SWAVES_newpng_url # If date is before 2019/01/01 png data exisits in old swaves repo

	### Searching for WIND/WAVES############
	tmpl=Wu+dtup[0]
	wv_png=''
	nlinks=[]
	tims=[]
	gtim=np.nan
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	if not os.path.isdir(outdir+'Wind/'):
		os.mkdir(outdir+'Wind/')
 
	try:
		req=urllib.request.urlopen(tmpl)
		links=bs(req,features="lxml").find_all('a',href=True)
		for ln in links:
			if "windwaves" not in ln.get("href"):
				continue
			tims+=[np.datetime64(dt.strptime(dtup[0]+ 'T'+ln.get("href").split('_')[1],"%Y/%m/%dT%H%M%S"))]
			nlinks+=[ln.get("href")]
		del links
		### If there are no data links to WAVES data thats also a bad period!
		if len(nlinks)==0:
			print('\n################\n No WIND data for ',dtup,'\n################')
  ####################################################################	
		else:
			tims=np.array(tims)
			delts=np.abs((tims-np.datetime64(dt.strptime(dtup[0]+'T'+dtup[1],"%Y/%m/%dT%H:%M:%S"))+np.timedelta64(delta_t,'m')).astype('timedelta64[s]').astype(float)/60)
			gtim=np.where(delts==np.min(delts))[0][0]
			if not os.path.isfile(outdir+'Wind/'+dtup[0].replace('/','')+'_'+dtup[1].replace(':','')+'_waves.png'):
				urllib.request.urlretrieve(tmpl+'/'+nlinks[gtim],outdir+'Wind/'+dtup[0].replace('/','')+'_'+dtup[1].replace(':','')+'_waves.png')
			print('\n########## Wind data found for ',dtup)
 
  ##### WAVES PNG file path being saved for concatinating later with Stereo
			wv_png=outdir+'Wind/'+dtup[0].replace('/','')+'_'+dtup[1].replace(':','')+'_waves.png'
  #############################################################
	except:
		print('\n################\n No WIND data for ',dtup,'\n################')

 ##### Checking if the date is before STEREO mission ###
	if np.datetime64(dtup[0].replace('/','-')+'T'+dtup[1]) <np.datetime64('2006-11-01T00:00:00'):
		print('\n--------------------\n Year before Stereo starting giving data!!\n--------------------------\n')
		return

 ########### STEREO data finding #############
	tmpl=SWu+dtup[0]
	flag=0
	timsA=[]
	timsB=[]
	nlinksA=[]
	nlinksB=[]
	gtimA=np.nan
	gtimB=np.nan
	stA_png=''
	stB_png=''
 
	if not os.path.isdir(outdir+'Stereo/'):
		os.mkdir(outdir+'Stereo/')
	try:
		print('Querrying....',tmpl)
		req=urllib.request.urlopen(tmpl)
		links=bs(req,features="lxml").find_all('a',href=True)
		print('Web data extracted!!')
		if SWu==SWAVES_png_url:
			for ln in links:
				if "swaves" not in ln.get("href"):
					continue
				if 'swaves1A' in ln.get("href"):
					timsA+= [np.datetime64(dt.strptime(dtup[0]+'T'+ln.get("href").split('_')[1],"%Y/%m/%dT%H%M%S"))]
					nlinksA+=[ln.get("href")]
					continue
				if 'swaves1B' in ln.get("href"):
					timsB+= [np.datetime64(dt.strptime(dtup[0]+'T'+ln.get("href").split('_')[1],"%Y/%m/%dT%H%M%S"))]
					nlinksB+=[ln.get("href")]
		else:
			for ln in links:
				if "sta_waves" not in ln.get("href"):
					continue
				timsA+= [np.datetime64(dt.strptime(dtup[0]+'T'+ln.get("href").split('_')[1],"%Y/%m/%dT%H%M%S"))]
				nlinksA+=[ln.get("href")]
   
		del links
		### If there are no data links to WAVES data thats also a bad period!
		if len(nlinksA)==0:
			print('\n################\n No STEREO A data for ',dtup,'\n################')
		else:
			print('\n-------------------\nUseful data links found for STEREO A !!!\n---------------------------')
			timsA=np.array(timsA)
			delts=np.abs((timsA-np.datetime64(dt.strptime(dtup[0]+'T'+dtup[1],"%Y/%m/%dT%H:%M:%S"))+np.timedelta64(delta_t,'m')).astype('timedelta64[s]').astype(float)/60)
			gtimA=np.where(delts==np.min(delts))[0][0]
			if not os.path.isfile(outdir+'Stereo/'+dtup[0].replace('/','')+'_'+dtup[1].replace(':','')+'_swavesA.png'):
				urllib.request.urlretrieve(tmpl+'/'+nlinksA[gtimA],outdir+'Stereo/'+dtup[0].replace('/','')+'_'+dtup[1].replace(':','')+'_swavesA.png')
			stA_png=outdir+'Stereo/'+dtup[0].replace('/','')+'_'+dtup[1].replace(':','')+'_swavesA.png'
			print('\n########## Stereo A data found for ',dtup)

		if len(nlinksB)==0:
			print('\n################\n No STEREO B data for ',dtup,'\n################')
		else:
			print('\n-------------------\nUseful data links found for STEREO B !!!\n---------------------------')
			timsB=np.array(timsB)
			delts=np.abs((timsB-np.datetime64(dt.strptime(dtup[0]+'T'+dtup[1],"%Y/%m/%dT%H:%M:%S"))+np.timedelta64(delta_t,'m')).astype('timedelta64[s]').astype(float)/60)
			gtimB=np.where(delts==np.min(delts))[0][0]
			if not os.path.isfile(outdir+'Stereo/'+dtup[0].replace('/','')+'_'+dtup[1].replace(':','')+'_swavesB.png'):
				urllib.request.urlretrieve(tmpl+'/'+nlinksB[gtimB],outdir+'Stereo/'+dtup[0].replace('/','')+'_'+dtup[1].replace(':','')+'_swavesB.png')
			stB_png=outdir+'Stereo/'+dtup[0].replace('/','')+'_'+dtup[1].replace(':','')+'_swavesB.png'
			print('\n########## Stereo B data found for ',dtup,'\n########################')
		if len(nlinksA)==0 and len(nlinksB)==0:
			print('\n########## No Stereo A & B data found for ',dtup,'\n########################')
			flag=1  
	except:
		flag=1

	if flag==1:
		return
 ########### Making combined images ##########
	gds=[True,True,True]
	png_lst=np.array([wv_png,stA_png,stB_png])

	if list(png_lst).count('')>1:
		print('Only 1 png obtained for the date time in list: ',png_lst,'. So not combining pngs!!')
		return
	if not os.path.isdir(outdir+'Combined/'):
		os.mkdir(outdir+'Combined/')
  
	if '' in png_lst:
		gds[np.where(png_lst=='')[0][0]]=False
		png_lst=png_lst[gds]
		print('\n--------------------\nCombining pngs: ',png_lst) 
		im1=Image.open(png_lst[0])
		im2=Image.open(png_lst[1])
		imn=Image.new('RGB',(max(im1.width,im2.width),im1.height+im2.height))
		imn.paste(im1,(0,0))
		imn.paste(im2,(0,im1.height))
		imn.save(outdir+'Combined/'+dtup[0].replace('/','')+'_'+dtup[1].replace(':','')+'_comb.png')  
		del imn
		print('\nCombined images made: ',png_lst,'\n-----------------')
	else:
		im1=Image.open(png_lst[0])
		im2=Image.open(png_lst[1])
		im3=Image.open(png_lst[2])
		print('\n--------------------\nCombining pngs: ',png_lst)
		imn=Image.new('RGB',(max(im1.width,im2.width,im3.width),im1.height+im2.height+im3.height))
		imn.paste(im1,(0,0))
		imn.paste(im2,(0,im1.height))
		imn.paste(im3,(0,im1.height+im2.height))
		imn.save(outdir+'Combined/'+dtup[0].replace('/','')+'_'+dtup[1].replace(':','')+'_comb.png')  
		del imn
		print('\nCombined images made: ',png_lst,'\n-----------------')
 
########### Calibrated DH band DS and light curve for a chosen frequency
def get_DHband_calibDS(sT,eT,outdir='./',sats=('sta','stb','wind'),Ch_frq=1.0, stats='median',vmin=None,vmax=None):
	'''
	Download SPDF CDF spectra for UTC sT--eT and save spectra/light curves.

	sT/eT: UTC strings YYYY-MM-DDThh:mm:ss; any number of days supported.
	outdir: PNGs/, DSpickle/ and CDFs/, each containing STEREO_A/,
		STEREO_B/ and Wind/.
	sats: sta/stb/wind, individually or as a sequence.
	Ch_frq: requested light-curve frequency in MHz; nearest channel used.
	stats: min/median/mean cadence statistic passed to Matrix_ops.Reg.
	vmin/vmax: optional log10(SFU) display limits.

	Return paths and errors for each spacecraft. Pickles use .p.
	Wind RAD1/RAD2 share one frequency axis and color scale.
	STEREO L3 PSD_SFU is used directly; the supplied Wind Z-channel
	gain conversion is preserved. Existing directories/CDFs are reused.
	Plots use 90 dpi, 18-point axis labels and 16-point tick labels.
	Requires the existing module-level imports.
	'''
	from spacepy import pycdf
	start,end=np.datetime64(sT,'us'),np.datetime64(eT,'us')
	if np.isnat(start) or np.isnat(end) or end<=start:
		raise ValueError('eT must be later than sT; both must be valid UTC times.')
	if stats not in ('min','median','mean'):
		raise ValueError('stats must be min, median or mean.')
	if not np.isfinite(Ch_frq) or Ch_frq<=0:
		raise ValueError('Ch_frq must be a positive frequency in MHz.')
	if vmin is not None and vmax is not None and vmax<=vmin:
		raise ValueError('vmax must exceed vmin.')
	if isinstance(sats,str):
		sats=[sats]
	if any(s not in ('sta','stb','wind') for s in sats):
		raise ValueError('sats must contain sta, stb and/or wind.')

	root=Path(outdir).expanduser().resolve()
	names={'sta':'STEREO_A','stb':'STEREO_B','wind':'Wind'}
	for category in ('PNGs','DSpickle','CDFs'):
		for name in names.values():
			(root/category/name).mkdir(parents=True,exist_ok=True)

	dates=np.arange(start.astype('datetime64[D]'),
		end.astype('datetime64[D]')+np.timedelta64(1,'D'))
	tag=str(start).split('.')[0].replace('-','').replace(':','')+'-'+str(end).split('.')[0].replace('-','').replace(':','')
	result={}

	for sat in dict.fromkeys(sats):
		name=names[sat]
		result[name]={'CDFs':[],'PNGs':[],'DSpickle':[],'errors':[]}
		status=result[name]
		parts={}

		for receiver in (('HFR',) if sat!='wind' else ('RAD1','RAD2')):
			chunks=[]
			for day in dates:
				date=str(day)
				if sat!='wind':
					branch='ahead' if sat=='sta' else 'behind'
					url=f'https://spdf.gsfc.nasa.gov/pub/data/stereo/{branch}/l3/waves/hfr-direction-finding/{date[:4]}/'
					prefix=f'{sat}_l3_wav_hfr_{date.replace("-","")}'
				else:
					url=f'https://spdf.gsfc.nasa.gov/pub/data/wind/waves/{receiver.lower()}_l2/{date[:4]}/'
					prefix=f'wi_l2_wav_{receiver.lower()}_{date.replace("-","")}'

				try:
					with urllib.request.urlopen(url,timeout=60) as stream:
						links=bs(stream.read(),'html.parser').find_all('a',href=True)
					files=[]
					for link in links:
						href=link['href']
						match=re.fullmatch(re.escape(prefix)+r'_v(\d+)\.cdf',href)
						if match:
							files.append((int(match[1]),href))
					if not files:
						raise FileNotFoundError('No matching CDF in '+url)

					href=max(files)[1]
					path=root/'CDFs'/name/href
					if not path.is_file():
						temp=path.with_suffix('.cdf.part')
						try:
							with urllib.request.urlopen(urljoin(url,href),timeout=60) as response:
								with temp.open('wb') as fil:
									while True:
										block=response.read(1024*1024)
										if not block:
											break
										fil.write(block)
							if temp.stat().st_size==0:
								raise ValueError('Empty CDF download.')
							temp.replace(path)
						finally:
							temp.unlink(missing_ok=True)
					status['CDFs'].append(str(path))

					with pycdf.CDF(str(path)) as cdf:
						variable=cdf['PSD_SFU' if sat!='wind' else 'PSD_V2_Z']
						values=np.asarray(variable[:],dtype=float)
						fill=variable.attrs.get('FILLVAL')
						if fill is not None:
							values[values==float(fill)]=np.nan
						for bound,op in (('VALIDMIN',np.less),('VALIDMAX',np.greater)):
							if bound in variable.attrs:
								values[op(values,variable.attrs[bound])]=np.nan
						times=np.asarray(cdf['Epoch'][:],dtype='datetime64[us]')
						freq=np.asarray(cdf['FREQUENCY'][:],dtype=float).squeeze()*1e-6

					if freq.ndim!=1 or values.shape!=(len(times),len(freq)):
						raise ValueError('Unexpected time/frequency/data dimensions.')
					good=np.isfinite(freq)&(freq>0)
					freq,values=freq[good],values[:,good]
					values[~np.isfinite(values)|(values<=0)]=np.nan

					if sat=='wind':
						gain=np.full(freq.shape,5.3027412e18 if receiver=='RAD1' else 7.2450814e18)
						if receiver=='RAD2':
							mask=freq>2.5
							gain[mask]=7.2450814e18*np.exp(1.39e-4*(freq[mask]*1000-2500))
						values*=gain
						values,times,freq=Reg(np.ma.masked_invalid(values),times,freq,
							day.astype('datetime64[us]'),(day+np.timedelta64(1,'D')).astype('datetime64[us]'),
							freq[0],freq[-1],'row',stats)
						values=np.asarray(np.ma.filled(values,np.nan))
						times=np.asarray(times,dtype='datetime64[us]')
						freq=np.asarray(freq,dtype=float)

					keep=(times>=start)&(times<=end)&~np.isnat(times)
					if np.any(keep):
						chunks.append((times[keep],freq,values[keep]))
				except Exception as exc:
					message=f'{name} {receiver} {date}: {exc}'
					print(message)
					status['errors'].append(message)

			if not chunks:
				continue
			try:
				freq=chunks[0][1]
				if any(f.shape!=freq.shape or not np.allclose(f,freq) for _,f,_ in chunks):
					raise ValueError('Frequency grid changes across days; cannot concatenate.')
				times=np.concatenate([t for t,_,_ in chunks])
				values=np.concatenate([v for _,_,v in chunks])
				order=np.argsort(times,kind='stable')
				times,values=times[order],values[order]
				_,unique=np.unique(times,return_index=True)
				times,values=times[unique],values[unique]
				order=np.argsort(freq)
				freq,ds=freq[order],values[:,order].T
				if len(times)<2 or len(freq)<2 or not np.any(np.isfinite(ds)&(ds>0)):
					raise ValueError('Insufficient usable spectrum samples.')
				parts[receiver]={'DS':ds,'Frequency (MHz)':freq,'Time (UT)':times,
					'Time (min)':(times-start)/np.timedelta64(1,'m')}
			except Exception as exc:
				status['errors'].append(f'{receiver}: {exc}')

		if not parts:
			continue
		try:
			p=root/'DSpickle'/name/(tag+'_DS.p')
			with p.open('wb') as fil:
				pickle.dump({'Receivers':parts,'Observatory':name,'sT':sT,'eT':eT,
					'Flux unit':'SFU','Source CDFs':status['CDFs']},fil)
			status['DSpickle'].append(str(p))

			logs=np.concatenate([np.log10(v['DS'][np.isfinite(v['DS'])&(v['DS']>0)]) for v in parts.values()])
			low=float(np.percentile(logs,5)) if vmin is None else float(vmin)
			high=float(np.percentile(logs,99)) if vmax is None else float(vmax)
			if high<=low:
				high=low+1

			#### All receivers on one frequency axis ####
			fig,ax=plt.subplots(figsize=(10,6),dpi=90)
			for receiver,values in parts.items():
				ds=values['DS']
				logged=np.full(ds.shape,np.nan)
				mask=np.isfinite(ds)&(ds>0)
				logged[mask]=np.log10(ds[mask])
				image=ax.pcolormesh(
					values['Time (UT)'],values['Frequency (MHz)'],
					np.ma.masked_invalid(logged),shading='auto',
					cmap='jet',vmin=low,vmax=high)

			ax.set_yscale('log')
			ax.set_ylabel('Frequency (MHz)',fontsize=18)
			ax.set_title(name,fontsize=20)
			ax.set_xlim(start,end)
			date=str(start)[:10]
			if str(start)[:10]!=str(end)[:10]:
				date+=' — '+str(end)[:10]
			xlabel='Time (UT) '+date
			ax.set_xlabel(xlabel,fontsize=18)
			ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
			ax.tick_params(axis='both',which='major',labelsize=16)
			cbar=fig.colorbar(image,ax=ax)
			cbar.set_label(r'$\log_{10}$[Flux density (SFU)]',fontsize=18)
			cbar.ax.tick_params(labelsize=15)
			fig.tight_layout()

			p=root/'PNGs'/name/(tag+'_DS.png')
			try:
				fig.savefig(p,dpi=90,bbox_inches='tight')
			finally:
				plt.close(fig)
			status['PNGs'].append(str(p))

			#### Nearest available channel within the requested receiver band ####
			candidates=[(abs(v['Frequency (MHz)'][np.argmin(abs(v['Frequency (MHz)']-Ch_frq))]-Ch_frq),r)
				for r,v in parts.items() if v['Frequency (MHz)'][0]<=Ch_frq<=v['Frequency (MHz)'][-1]]
			if not candidates:
				status['errors'].append(f'No receiver covers {Ch_frq} MHz; no light curve saved.')
				continue

			receiver=min(candidates)[1]
			values=parts[receiver]
			channel=np.argmin(abs(values['Frequency (MHz)']-Ch_frq))
			freq=float(values['Frequency (MHz)'][channel])
			flux=values['DS'][channel]
			p=root/'DSpickle'/name/(tag+'_LC.p')
			with p.open('wb') as fil:
				pickle.dump({'Frequency (MHz)':freq,'Requested frequency (MHz)':Ch_frq,
					'Light curve':flux,'Time (UT)':values['Time (UT)'],'Flux unit':'SFU',
					'Receiver':receiver},fil)
			status['DSpickle'].append(str(p))

			for logarithmic in (False,True):
				fig,ax=plt.subplots(figsize=(10,6),dpi=90)
				ax.plot(values['Time (UT)'],np.where(flux>0,flux,np.nan))
				if logarithmic:
					ax.set_yscale('log')
				ax.set_xlim(start,end)
				ax.set_ylabel('Flux density (SFU)',fontsize=18)
				ax.set_xlabel(xlabel,fontsize=18)
				ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
				ax.set_title(f'{name} {receiver}: {freq:g} MHz',fontsize=20)
				ax.tick_params(axis='both',which='major',labelsize=16)
				fig.tight_layout()
				p=root/'PNGs'/name/(tag+('_LC_log.png' if logarithmic else '_LC.png'))
				try:
					fig.savefig(p,dpi=90,bbox_inches='tight')
				finally:
					plt.close(fig)
				status['PNGs'].append(str(p))
		except Exception as exc:
			message=f'{name} output failed: {exc}'
			print(message)
			status['errors'].append(message)
	return result

###### LASCO / SDO / EIT image acquisition and calibration ####################


# =============================================================================
# LASCO
# =============================================================================


def get_LASCO_imgs(sT, eT, dT, instr, outdir='./LASCO_1', need_calib=False, helioviewer=True):
	"""Return LASCO Maps; originals stay in outdir/raw.

	sT/eT: UTC 'YYYY-MM-DDThh:mm:ss'; dT: cadence seconds; instr: 'C2'/'C3'.
	need_calib=True: route Level-1 -> Helioviewer -> Level-0.5. Verified
		Level-1/MSB FITS from VSO/SDAC (saved unchanged in outdir/maps); if none
		is usable, Helioviewer images (north-up display values, canonical tags,
		outdir/maps_HV); only if Helioviewer also fails, archive Level-0.5 DN/s.
	need_calib=False: archive Level-0.5 FITS, exposure-normalized to DN/s and
		FoV-masked, saved in outdir/maps.
	helioviewer=True (need_calib=False only): also save a Helioviewer JP2 for
		every cadence bin in outdir/raw/helioviewer_<instr>/ (not returned);
		jp2_to_fits(..., outdir/'maps_HV') converts them if needed.
	No SolarSoft is attempted here; calib_LASCO() can try it on Level-0.5 data.
	Roll: every saved FITS has a truthful WCS. Level-0.5/DN/s keep the
	spacecraft CROTA (pixels unrotated); Level-1 is saved unchanged; Helioviewer
	pixels are north-up, so CROTA=0 and the roll is kept in OV_ROLL.
	"""
	if not isinstance(need_calib, (bool, np.bool_)):
		raise ValueError('need_calib must be True or False.')
	instr = instr.upper()
	if instr not in _helper.LASCO_FOV:
		raise ValueError("instr must be 'C2' or 'C3'.")
	_helper._target_times(sT, eT, dT)  # Validate before entering network fallbacks.
	root = Path(outdir).expanduser().resolve()
	if need_calib:
		try:
			maps = _helper._get_lasco_source(sT, eT, dT, instr, root/'raw', 'level1')
			_helper._lasco_notice(maps)
			_save_image_previews(maps,root/'PNGs','LASCO',instr)
			return _helper._save_maps(maps, root/'maps', 'LASCO_' + instr)
		except Exception as exc:
			warnings.warn(f'LASCO {instr} Level-1 retrieval failed: {exc}; trying Helioviewer next.')
		try:
			maps = _helper._get_lasco_source(sT, eT, dT, instr, root/'raw', 'helioviewer', calibrate=True)
			_save_image_previews(maps,root/'PNGs','LASCO',instr)
			maps = _helper._save_maps(maps, root/'maps_HV', 'LASCO_' + instr + '_HV')
			_helper._lasco_notice(maps)
			return maps
		except Exception as exc:
			warnings.warn(f'LASCO {instr} Helioviewer retrieval failed: {exc}; last fallback: archive Level-0.5 DN/s.')
	elif helioviewer:
		try:
			jp2s=_helper._fill_helioviewer_jp2(sT,eT,dT,[],root/'raw'/f'helioviewer_{instr}','LASCO',instr)
			print(f'LASCO/{instr}: {len(jp2s)} Helioviewer JP2s in the window.')
		except Exception as exc:
			warnings.warn(f'LASCO {instr} Helioviewer JP2 retrieval failed: {exc}')
	maps = _helper._get_lasco_source(sT, eT, dT, instr, root/'raw', 'archive', calibrate=True)
	for image in maps:
		if not _helper._lasco_science_ready(image) and not image.meta.get('ov_hv', False):
			image.meta['ov_full'] = False
			image.meta['ov_cal'] = 'Level-0.5 exposure-normalized DN/s; not full Level-1 calibration'
	_save_image_previews(maps,root/'PNGs','LASCO',instr)
	maps = _helper._save_maps(maps, root/'maps', 'LASCO_' + instr)
	_helper._lasco_notice(maps)
	return maps


def get_helioviewer_jp2(sT,eT,dT,instrument,wave=None,spacecraft='A',outdir='./helioviewer'):
	"""Save a Helioviewer JP2 for every cadence bin; return their paths.

	instrument: 'LASCO' (wave='C2'/'C3'), 'AIA' (wave in Angstrom), 'HMI'
	(wave='magnetogram'/'continuum'), 'EUVI' (wave 171/195/284/304), 'COR1'
	or 'COR2' (spacecraft 'A'/'B'). Cached JP2s in outdir are reused.
	JP2 pixels are north-up display values; use jp2_to_fits for FITS copies.
	"""
	return _helper._fill_helioviewer_jp2(sT,eT,dT,[],Path(outdir),str(instrument).upper(),wave,str(spacecraft).upper())


def calib_LASCO(images, instr=None, outdir='./LASCO_calibrated'):
	"""Prepare existing LASCO FITS/Maps: SolarSoft -> Level-1 download -> JP2.

	images: FITS path, SunPy Map, or list/tuple; return corresponding Map/list.
	instr: C2/C3, or None to infer from headers. outdir receives FITS and logs.
	Verified Level-1/MSB inputs pass unchanged. SolarSoft uses the executable
	'sswidl' on PATH and requires LASCO software/current calibration files.
	DN/s products created by this module retain OV_RAW pointing to original DN;
	if originals are absent, skip local calibration and try remote replacements.
	Remote replacements must be within 180 seconds, retaining their actual time;
	Level-1 native exposure IDs must match when available. If all routes fail,
	return an explicitly uncalibrated DN/s/FoV preview. JP2 is never MSB data.
	"""
	if isinstance(images, (list, tuple)):
		return [calib_LASCO(x, instr, outdir) for x in images]
	m = sunpy.map.Map(images) if isinstance(images, (str, Path)) else images
	if not isinstance(m, sunpy.map.GenericMap):
		raise TypeError('images must be FITS paths or SunPy Maps.')
	detector = str(m.meta.get('detector', '')).strip().upper()
	if detector not in _helper.LASCO_FOV or (instr is not None and instr.upper() != detector):
		raise ValueError('Input must be LASCO C2/C3 matching instr.')
	root = Path(outdir).expanduser().resolve()
	if _helper._lasco_science_ready(m):
		return _helper._save_maps([m], root/'maps', 'LASCO_' + detector)[0]
	original = m.meta.get('ov_raw')
	try:
		raw = sunpy.map.Map(original) if original and Path(original).is_file() else m
		if m.meta.get('ov_hv', False):
			raise ValueError('JP2 display data cannot be calibrated as raw DN.')
		result = _helper._calib_lasco_backend(raw, instr=detector, outdir=root/'solarsoft', fallback_level1=False)
		return _helper._save_maps([result], root/'maps', 'LASCO_' + detector)[0]
	except Exception as exc:
		warnings.warn(f'SolarSoft calibration unavailable/failed: {exc}; trying remote products.')
	for source in ('level1', 'helioviewer'):
		try:
			candidates = _helper._get_lasco_source((m.date-180*u.s).isot, (m.date+180*u.s).isot,
			 180., detector, root/'raw', source, strict_filters=False, calibrate=(source == 'helioviewer'))
			if source == 'level1':
				identity = _helper._lasco_exposure_id(m)
				if identity:
					candidates = [x for x in candidates if _helper._lasco_exposure_id(x) == identity]
			if not candidates:
				raise ValueError('No matching LASCO exposure found.')
			result = min(candidates, key=lambda x: abs((x.date-m.date).to_value(u.s)))
			separation = float(abs((result.date-m.date).to_value(u.s)))
			if separation > 180:
				raise ValueError('Replacement observation exceeds 180 seconds.')
			result.meta['ov_origt'], result.meta['ov_tdiff'] = m.date.isot, separation
			warnings.warn(f'Using {source} observation {result.date.isot}; offset {separation:.3f} s.')
			_helper._lasco_notice([result])
			return _helper._save_maps([result], root/'maps', 'LASCO_' + detector)[0]
		except Exception as exc:
			warnings.warn(f'LASCO {source} replacement failed: {exc}.')
	result = _helper._prepare_lasco(m, detector)
	if not m.meta.get('ov_hv', False):
		result.meta['ov_full'] = False
		result.meta['ov_cal'] = 'Level-0.5 DN/s preview; full calibration/retrieval failed'
	_helper._lasco_notice([result])
	return _helper._save_maps([result], root/'maps', 'LASCO_' + detector)[0]

# =============================================================================
# SDO (AIA / HMI)
# =============================================================================

def get_SDO_imgs(sT,eT,dT,instruments=['AIA','HMI'],
	waves=[[94,171,193],['magnetogram']],outdir='./SDO',
	nproc=4,helioviewer=True):
	'''
	Get AIA/HMI FITS Maps in parallel for instrument-wave combinations.

	Raw FITS: outdir/<instrument>/<wave>/raw/
	JP2 previews: outdir/<instrument>/<wave>/raw/helioviewer/
	Finished FITS: outdir/<instrument>/<wave>/

	Return FITS-based Map lists in input combination order.
	JP2s remain separate; AIA preparation is handled by find_coobs.
	Existing nested caches remain searchable. Cached FITS are checked first;
	only missing cadence slots are queried in bounded windows up to +/-60 s.
	A failed slot does not discard other successful downloads.
	'''
	if len(instruments)!=len(waves):
		raise ValueError('instruments and waves must have matching lengths.')
	if int(nproc)!=nproc or nproc<1:
		raise ValueError('nproc must be a positive integer.')
	targets=list(_helper._target_times(sT,eT,dT))
	start,end=Time(sT),Time(eT)
	if targets[-1]<end:targets.append(end)

	jobs=[]
	for instrument,products in zip(instruments,waves):
		instrument=instrument.upper()
		if instrument not in ('AIA','HMI'):
			raise ValueError("instrument must be 'AIA' or 'HMI'")
		for wave in products:
			if instrument=='AIA' and wave not in _helper.AIA_WAVES:
				raise ValueError(f'wave must be one of {_helper.AIA_WAVES}')
			if instrument=='HMI':
				wave='magnetogram' if wave is None else wave
				if wave not in ('magnetogram','dopplergram','continuum'):
					raise ValueError('Invalid HMI product.')
			if (instrument,wave) not in jobs:
				jobs.append((instrument,wave))

	def download(instrument,wave):
		destination=Path(outdir).expanduser().resolve()/instrument/str(wave)
		raw=destination/'raw'
		raw.mkdir(parents=True,exist_ok=True)
		# Reuse cached FITS before requesting anything. Read each file only once.
		cache={};seen={}
		def scan_new_fits():
			for path in _helper._fits_paths(destination):
				path=Path(path).resolve()
				stat=path.stat();signature=(stat.st_mtime_ns,stat.st_size)
				if seen.get(path)==signature:continue
				seen[path]=signature
				try:
					m=sunpy.map.Map(path)
					if m.meta.get('ov_hv',False):continue
					if instrument not in (str(m.instrument)+' '+str(m.detector)).upper():continue
					if instrument=='AIA':
						if not np.isclose(m.wavelength.to_value(u.angstrom),wave):continue
						prepared=float(m.meta.get('lvl_num',0))>=1.5 or any(
							value in str(m.meta.get('series','')).lower() for value in ('lev1.5','lev15'))
						if not prepared:
							warnings.warn(f'AIA FITS is not prepared Level 1.5: {path}')
							continue
					if start<=m.date<=end:cache[path]=m.date
				except Exception as exc:
					warnings.warn(f'Skipping unreadable cached FITS {path}: {exc}')
		scan_new_fits()
		tolerance=float(dT)/2
		padding=min(60.,tolerance)
		for target in targets:
			if any(abs((time-target).to_value(u.s))<=tolerance for time in cache.values()):
				continue
			lo=max(start,target-padding*u.s)
			hi=min(end,target+padding*u.s)
			print(f'Missing {instrument}/{wave}: {lo.isot} -- {hi.isot}')
			cwd=os.getcwd()
			try:
				find_coobs(lo.isot,hi.isot,[instrument],
					[[float(wave)*u.angstrom]] if instrument=='AIA' else [[]],
					outdir=str(raw)+os.sep,downdata=True,dT=float(dT),
					HMI_product=wave if instrument=='HMI' else None,
					flat_outdir=True)
			except Exception as exc:
				warnings.warn(f'{instrument}/{wave} FITS retrieval failed at {target.isot}: {exc}')
			finally:
				os.chdir(cwd)
			scan_new_fits()
		# Keep only real exposures within half a cadence of each requested slot.
		paths=[]
		for target in targets:
			near=[p for p,time in cache.items() if abs((time-target).to_value(u.s))<=tolerance]
			if not near:continue
			path=min(near,key=lambda p:(abs((cache[p]-target).to_value(u.s)),p.parent!=destination,str(p)))
			if path not in paths:paths.append(path)
		paths.sort(key=lambda p:cache[p].mjd)

		maps=[]
		for path in paths:
			try:
				maps.append(sunpy.map.Map(path))
			except Exception as exc:
				warnings.warn(f'Skipping unreadable FITS {path}: {exc}')
		if maps:
			maps=_helper._save_maps(
				maps,destination,instrument+'_'+str(wave))

		_save_image_previews(maps,destination/'PNGs',instrument,wave)
		if helioviewer:
			try:
				jp2s=_helper._fill_helioviewer_jp2(
					sT,eT,dT,[],raw/'helioviewer',instrument,wave)
				_save_image_previews(jp2s,raw/'helioviewer',instrument,wave)
				print(f'{instrument}/{wave}: {len(maps)} FITS Maps; '
					f'{len(jp2s)} JP2 fallback previews.')
			except Exception as exc:
				warnings.warn(f'{instrument}/{wave} JP2 fallback failed: {exc}')
		return maps

	if not jobs:
		return []
	workers=min(int(nproc),len(jobs))
	if workers==1:
		return [download(*job) for job in jobs]
	with Pool(workers) as pool:
		return pool.starmap(download,jobs)
		
def get_EIT_imgs(sT, eT, dT, wave, outdir='./EIT',
  archive_url='https://lasco6.nascom.nasa.gov/pub/eit/l1/'):
	"""Return and save calibrated official EIT Level-1 SunPy Maps.

	sT/eT: UTC interval; dT: cadence seconds; wave: 171, 195, 284 or 304 Angstrom.
	outdir stores the selected FITS. Raw EIT inputs are rejected. archive_url is the
	Level-1 archive root, override only for an equivalent mirror. No SolarSoft
	installation is required to download products already processed by eit_prep.
	Native cadence and archive gaps limit actual sampling; duplicates are merged.
	"""
	if wave not in (171, 195, 284, 304):
		raise ValueError('EIT wave must be 171, 195, 284 or 304 Angstrom')
	targets = _helper._target_times(sT, eT, dT)
	folder = Path(outdir).expanduser().resolve()
	folder.mkdir(parents=True, exist_ok=True)
	rows = []
	day = _helper._parse_time(sT).datetime.date()
	while day <= _helper._parse_time(eT).datetime.date():
		url = archive_url.rstrip('/') + day.strftime('/%Y/%m/%d/')
		try:
			with urllib.request.urlopen(url, timeout=60) as stream:
				index = stream.read().decode()
		except urllib.error.HTTPError as exc:
			if exc.code != 404:
				raise
		else:
			for filename in sorted(set(re.findall(r'SOHO_EIT_\d+_\d{8}T\d{6}_L1\.fits', index))):
				match = re.fullmatch(r'SOHO_EIT_(\d+)_(\d{8}T\d{6})_L1.fits', filename)
				t = Time.strptime(match[2], '%Y%m%dT%H%M%S')
				if int(match[1]) == wave and _helper._parse_time(sT) - 600*u.s <= t <= _helper._parse_time(eT) + 600*u.s:
					rows.append((t, url + filename))
		day += timedelta(days=1)
	if not rows:
		raise FileNotFoundError('No EIT Level-1 images for this wavelength and interval; no raw substitution was made.')
	tab = Table({'Start Time': Time([r[0] for r in rows]), 'Archive URL': [r[1] for r in rows]})
	got = _helper._fetch_nearest(tab, targets, folder,
		validate=lambda p: _helper._parse_time(sT) <= _helper._hdr_time(p) <= _helper._parse_time(eT))
	paths = [p for _, p in got]
	maps = []
	for path in paths:
		m = sunpy.map.Map(path)
		if 'EIT' not in str(m.instrument).upper() or int(m.wavelength.to_value(u.angstrom)) != wave:
			raise ValueError(f'{path}: not the requested EIT wavelength')
		_helper._require_eit_level1(m, path)
		maps.append(m)
	_save_image_previews(maps,folder/'PNGs','EIT',wave)
	return _helper._save_maps(maps, folder, 'EIT' + str(wave))

def _stereo_job(sT,eT,dT,spacecraft,instrument,wave,outdir,polar=0):
	"""One get_STEREO_imgs job: HI route for HI1, SSC/VSO series otherwise."""
	if instrument=='HI1':
		return _helper._hi_series(sT,eT,dT,spacecraft,instrument,outdir)
	return _helper._stereo_series(sT,eT,dT,spacecraft,instrument,wave,outdir,polar)

def get_STEREO_imgs(sT,eT,dT,spacecraft='A',instruments=['COR2'], waves=None,outdir='./STEREO',nproc=4,polar=0,helioviewer=True):
	'''
	Get STEREO FITS images in parallel for instrument-wave combinations.

	spacecraft: 'A' or 'B'.
	instruments: EUVI, COR1, COR2 or HI1; a string or a sequence.
	waves: one entry per instrument, e.g. [[171,195],None,None].
		waves=None selects all four EUVI bands.
	nproc: number of FITS download workers; 1 runs serially.
	polar: COR1 FITS polarizer angle, 0/120/240 degrees.
	COR2: total-brightness only (POLAR=1001); single-polarizer files are skipped.
	helioviewer: save original JP2 previews for missing cadence bins.
	HI1: route Level-1 -> Helioviewer -> Level-0.5, like COR2. Level-1 HI
		(hi_prep DN/s, MSB or S10) is only at UKSSDC behind a registered
		login, so local Level-1 files under STEREO_<sc>/HI1 (e.g. raw/level1)
		are used first; Helioviewer has no HI source (step reported and
		skipped); otherwise Level-0.5 science sums (s*h1*.fts, N_IMAGES>1)
		come from SSC/VSO into STEREO_<sc>/HI1/maps. Bins are
		max(dT, 2400 s), the HI1 science cadence.

	Returns one time-sorted list of FITS-based Maps per combination.
	JP2s are saved separately and are not returned as Maps.
	'''
	_helper._target_times(sT,eT,dT)
	spacecraft=str(spacecraft).strip().upper()
	if spacecraft not in ('A','B'):
		raise ValueError("spacecraft must be 'A' or 'B'.")
	if isinstance(instruments,str):
		instruments=[instruments]
	instruments=[str(i).strip().upper() for i in instruments]
	if any(i not in ('EUVI','COR1','COR2','HI1') for i in instruments):
		raise ValueError("instruments must contain EUVI, COR1, COR2 or HI1.")
	if int(nproc)!=nproc or nproc<1:
		raise ValueError('nproc must be a positive integer.')
	if 'COR1' in instruments and polar not in (0,120,240):
		raise ValueError('COR1 polar must be 0, 120 or 240 degrees.')

	if waves is None:
		waves=[[171,195,284,304] if i=='EUVI' else None
			for i in instruments]
	if len(instruments)!=len(waves):
		raise ValueError('instruments and waves must have matching lengths.')

	jobs=[]
	for instrument,products in zip(instruments,waves):
		if instrument=='EUVI':
			if products is None:
				raise ValueError('Supply EUVI wavelengths, e.g. [171,195].')
			if np.isscalar(products):
				products=[products]
			products=list(products)
			if not products or any(w not in (171,195,284,304) for w in products):
				raise ValueError('EUVI wavelengths must be 171, 195, 284 or 304.')
		else:
			products=[None]
		for wave in products:
			job=(sT,eT,dT,spacecraft,instrument,wave,outdir,polar)
			if job not in jobs:
				jobs.append(job)

	if not jobs:
		return []
	workers=min(int(nproc),len(jobs))
	if workers==1:
		results=[_stereo_job(*job) for job in jobs]
	else:
		with Pool(workers) as pool:
			results=pool.starmap(_stereo_job,jobs)

	for job,maps in zip(jobs,results):
		instrument,wave=job[4:6]
		product=Path(outdir).expanduser().resolve()/f'STEREO_{spacecraft}'/instrument
		if instrument=='EUVI':product=product/str(int(wave))
		_save_image_previews(maps,product/'PNGs',instrument,wave)
	if helioviewer:
		root=Path(outdir).expanduser().resolve()/f'STEREO_{spacecraft}'
		for job,maps in zip(jobs,results):
			instrument,wave=job[4:6]
			if instrument=='HI1':continue  # no Helioviewer HI source
			product=root/instrument
			if instrument=='EUVI':
				product=product/str(int(wave))
			jp2s=_helper._fill_helioviewer_jp2(
				sT,eT,dT,[],product/'raw'/'helioviewer',
				instrument,wave,spacecraft)
			_save_image_previews(jp2s,product/'raw'/'helioviewer',instrument,wave)

	return results