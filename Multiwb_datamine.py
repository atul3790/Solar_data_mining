'''
All observation search modules.

'''
from sunpy.net import Fido, attrs as a
import numpy as np,pickle,glob
from PIL import Image
import sunpy.timeseries as ts
from astropy import units as u
from astropy.table import Table, Row, Column
from radiospectra import net
#from radiospectra.sources.callisto import CallistoSpectrogram
import zipfile as zp
from datetime import datetime as dt, timedelta,timezone
import urllib,os
from bs4 import BeautifulSoup as bs
from radiospectra.spectrogram import Spectrogram
import matplotlib.pyplot as plt
from aiapy.calibrate import register, update_pointing
from sunpy.map import Map
import pandas as pd
import ssl

ssl._create_default_https_context = ssl._create_unverified_context

### MWA obs query for a few date times
def MWAavail(sT,DT=30,outfold='./'):
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
		pickle.dump(taBl,open(outfold+sT.replace(':','')+'-'+eT.replace(':','')+'_MWAobsIDs.p','wb'))		
		#pickle.dump(taBl,open('MWA_obsIDs/'+sT+'-'+eT+'_MWAobsIDs.p','wb'))
		print('MWA Obs found with ',len(taBl),' IDs saved to ',sT+'-'+eT+'_MWAobsIDs.p')
		return True
	except:
		return False
#################### Download and plot NoRP data 
#### Nobeyama Radio Polarimeter sun-as-a-star data
def find_norp(sT,eT,outdir='./',plot=True):
	'''
	sT is the datetime string of the event 'YYYY-MM-DDThh:mm:ss'.
	eT is the end datetime string of the event 'YYYY-MM-DDThh:mm:ss'.
	outdir is the path to directory where to store the light curve and quick look image data files
	plot =True ensures the data will be plotted for the sT to eT period
	'''
	from astropy.io import fits
	import urllib.request
	from matplotlib.dates import DateFormatter
	dtFmt=DateFormatter('%b-%d %H:%M')
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
	cache=outdir+subdir+'MaxFlux_'+str(sT).replace(':','')+'-'+str(eT).replace(':','')+'.csv'
	if os.path.isfile(cache):
		if plot==False:
			print('Analysis already completed for the event:',Flrtim)
			return
		Frqs=pd.read_csv(cache,usecols=['Frequency (GHz)'])['Frequency (GHz)'].to_numpy()
		getall=[outdir+subdir+'Flux_'+str(int(np.round(float(fq),0)))+'GHz_'+str(sT).replace(':','').replace('-','')+'-'+str(eT).replace(':','').replace('-','')+'.png' for fq in Frqs]
		if len(getall)>0 and all(os.path.isfile(i) for i in getall):
			print('Analysis and all plots already present! Skipping FITS reads and downloads.')
			return
	###############################################################
	badurls=[]	
	for url in urls:
		####### Checking existing downloaded material and analysis
		if os.path.isfile(outdir+subdir+url.split('/')[-1]):
			print(url,' already analysed!')
			if plot==True and len(glob.glob(outdir+subdir+'*.png'))>0:
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

			if os.path.isfile(outdir+subdir+'MaxFlux_'+str(t1).replace(':','')+'-'+str(t2).replace(':','')+'.csv'):
				print('Analysis completed for the event: ',Flrtim)
				if plot==True:
					getall=[outdir+subdir+'Flux'+tag+'_'+str(t1).replace(':','').replace('-','')+'-'+str(t2).replace(':','').replace('-','')+'.png' for tag in Alltags]
					if all(os.path.isfile(i) for i in getall):
						print('All plots are in place for flare within',sT,' - ',eT,' So doing next time slot!! ')
						continue
					print('Not all frequency LC made; recomputing all frequencies to retain the full summary.')

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
				VI+=[Vm/Im]
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
					pickle.dump({'Time':times,'Flux_I':Flx_I,'Flux_V':Flx_V,'Flux V/I':Flx_P,'Start time':t1,'End time':t2,'I Flux range':[Iy1,Iy2],'V Flux range':[y1,y2]},open(outdir+subdir+'FluxLC'+tag+'_'+str(t1).replace(':','')+'-'+str(t2).replace(':','')+'.p','wb'))
					continue
					
				# Ylims in plot and max locations.
				if tog_I==True:
					Iy1,Iy2=np.round(np.nanmin(Flx_I[l1:l2+1])*0.9)-1,np.round(Im*1.1)+1
					Ipls=np.where(Flx_I[l1:l2+1]==Im)[0]
					IMxt+=[str(times[l1:l2+1][Ipls[0]]).replace('T',' ')]
					nIm+=[len(Ipls)]
					## Identify flare period
					#Flocs=np.where((Flx_I[l1:l2+1]>Im/2) & (Flx_I[l1:l2+1]>1.5*Imd))[0] # Definition of flare duration : Flux is atleast above 1.5 x median Flux & flux > Half Max
					Flocs=np.where(((Flx_I[l1:l2+1]-Imd)>(Im-Imd)/2) & (Flx_I[l1:l2+1]>1.15*Imd))[0] # Definition of flare duration : Flux is atleast above 1.5 x median Flux & flux > Half Max
					if len(Flocs)>0:
						FI=times[l1:l2+1][Flocs]
						Ivl=[FI[0],FI[-1]]	# Will mark vertical period for flare flux estimation
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
					VMxt+=[str(times[l1:l2+1][Vpls[0]]).replace('T',' ')]
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
					plt.savefig(outdir+subdir+'Flux'+tag+'_'+str(t1).replace(':','').replace('-','')+'-'+str(t2).replace(':','').replace('-','')+'.png')
					plt.close()
				
				pickle.dump({'Time':times,'Flux_I':Flx_I,'Flux_V':Flx_V,'Flux V/I':Flx_P,'Start time':t1,'End time':t2,'I Flux range':[Iy1,Iy2],'V Flux range':[y1,y2]},open(outdir+subdir+'FluxLC'+tag+'_'+str(t1).replace(':','')+'-'+str(t2).replace(':','')+'.p','wb'))			

			if os.path.abspath(fil)!=os.path.abspath(outdir+subdir+os.path.basename(fil)):
				os.replace(fil,outdir+subdir+os.path.basename(fil))
			IFlxs=np.array(IFlxs)
			VFlxs=np.array(VFlxs)
			DF=pd.DataFrame({'Frequency (GHz)':Frqs,'Median_I_flux (SFU)':Imed,'Median_V_flux (SFU)':Vmed,'Max_I_flux (SFU)':IFlxs,'Max_V_flux (SFU)':VFlxs,'Max Pol':VFlxs/IFlxs,'Max_to_median_I (SFU)':IFlxs/np.array(Imed),'Max_to_median_V (SFU)':VFlxs/np.array(Vmed),'Max_I_time (UT)':IMxt,'Max_V_time (UT)':VMxt,'Max_I_repetition':nIm,'Max_V_repetition':nVm,'Flare_Energy_I (SFU s)':TotI,'Flare_Energy_V (SFU s)':TotV,'Mean pol':np.array(TotV)/np.array(TotI)})
			DF=DF.astype(str)
			DF.to_excel(outdir+subdir+'MaxFlux_'+str(t1).replace(':','')+'-'+str(t2).replace(':','')+'.xlsx',index=False)
			DF.to_csv(outdir+subdir+'MaxFlux_'+str(t1).replace(':','')+'-'+str(t2).replace(':','')+'.csv',index=False)
			pickle.dump(DF,open(outdir+subdir+'MaxFlux_'+str(t1).replace(':','')+'-'+str(t2).replace(':','')+'.p','wb'))
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
	from astropy.io import fits
	import urllib.request
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
					edg=data.shape[0]/2+0.5
					x1=np.round(edg*dx/60,1)
					y1=np.round(edg*dy/60,1)
					extent=[-x1,x1,-y1,y1]
					if polz=='STOKES I':
						data=np.ma.masked_less_equal(data,1000)
						plt.imshow(np.ma.log10(data),aspect='auto',extent=extent,origin='lower',cmap='gist_heat')
						h=plt.colorbar()
						h.set_label(r'$\mathrm{log_{10}T_B\ (K)}$',size=17)
					else:
						vmax=np.nanmax(np.abs(data))
						if not np.isfinite(vmax) or vmax==0:
							vmax=1.
						plt.imshow(data,aspect='auto',extent=extent,origin='lower',cmap='RdBu_r',vmin=-vmax,vmax=vmax)
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
	Return IDs and catalogue rows for flares starting after sT
	and ending before eT.

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
	days=[(stc+timedelta(days=i)).strftime('%d-%b-%Y')
		for i in range((etc.date()-stc.date()).days+1)]
	sTab=flcat[np.isin(flcat['Date'],days)]

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
	tst=np.array(tst,dtype=object)
	tend=np.array(tend,dtype=object)
	goods=(tst>stc)&(tend<etc)
	Tab=sTab[goods][np.argsort(tst[goods])]
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
	import re,shutil,urllib.request
	from urllib.parse import urljoin,urlsplit

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
			Tab['ID'][Tab['Date']==day.strftime('%d-%b-%Y')]]

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

######### Find otherwave coobs
def find_coobs(sT,eT,instrs,wvs,outdir='./',RHESSI_flarecat="hessi_flare_catalog.csv",downdata=False):
	'''
	sT is the start datetime string of the event 'YYYY-MM-DDThh:mm:ss'.
	eT is the end datetime string of the event 'YYYY-MM-DDThh:mm:ss'.
	instrs = ['EIS', 'AIA'] - list of instruments of choice. or just a string 'EIT'
	wvs=[[],[]] list of list of wavelengths as floats with astropy units for each instrument. leave as [] if all data is needed for an instrument
	outdir = path to directory where in case any data is found will be saved
	downdata=True will download observations found outdir

	'''
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
	cwd=os.getcwd()
		
	for i in range(len(instrs)):
		ins=instrs[i]
		wv=wvs[i]
		if ins.lower() in RSTN_list:
			print(ins,' is an RSTN network instrument. Please use get_RSTN() from nexttime for more functionality and ease.')
			print('Calling get_RSTN() functionality! ')
			get_RSTN(sT,eT,instrs=[ins],outdir=outdir)
			continue
		os.chdir(cwd)
		if not os.path.isdir(outdir+ins):
			os.mkdir(outdir+ins)
		os.chdir(outdir+ins)

		if len(wv)!=0:
			for wvj in wv:
				print('Searching for data: ',ins,'; Wavelength: ',wvj.value,str(wvj.unit))
				try:
					res=Fido.search(a.Time(sT.replace('T',' '),eT.replace('T',' ')),a.Instrument(ins),a.Wavelength(wvj))
					print('Data found for ',ins,'\nWavelength: ',wvj)
				except:
					if 'rhessi' in ins.lower():
						res=find_RHESSI(sT,eT,outdir='./',flare_cat=RHESSI_flarecat,baseurl='https://data.darts.isas.jaxa.jp/pub/mirror/rhessi/')
					else:
						bads+=[(ins,wvj)]
					continue
					
				if len(res)==0:
					print('No data found for ',sT,' - ',eT,' period\n',ins,' at Wavelength: ',wvs[i])
					continue
				print('Provider: ',res.keys())
				prv=res.keys()[0]

				#pickle.dump(Table(res[prv]),open(ins+'_'+str(wvj.value).replace('.','')+str(wvj.unit)+'_search_results.p','wb'))	
				if downdata==True:
					print('Downloading.....',ins,' data for ',wvj)
					print(res[prv])
					if not os.path.isdir(str(wvj.value)+str(wvj.unit)+'/'):
						os.mkdir(str(wvj.value)+str(wvj.unit)+'/')
					bdir=os.getcwd()	
					os.chdir(str(wvj.value)+str(wvj.unit)+'/')	
					dwnfils=Fido.fetch(res[prv],path='./')
		
					for fil in dwnfils.data:
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
									AMap=update_pointing(AMap) # Doing rotation corrections
									AMap15=register(AMap) # Pixels fluxes are now corrected for direction offset issue and pixels rescaled to 0.6" and image y axes are now aligned with solar north-south axis
									Amap15n=AMap15/AMap15.exposure_time # Making map normalised with respect to exposure times to get values in DN/s/pixel

									data_lev1=fits.open(fil,mode='update')
									data_lev1[1].data=Amap15n.data
									data_lev1[1].header['EXPTIME']=1.00		
									data_lev1.data=AMap.data
									data_lev1[1].header['SERIES']='aia.lev1.5'
									data_lev1[1].header['SEGMENT']='image_lev1.5.fits'
									data_lev1[1].header['PIXLUNIT']='DN/s'
									data_lev1.flush()
									data_lev1.close()
									os.system('mv '+fil+' '+fil.replace('lev1','lev15'))
									del AMap
									del AMap15	

						if ins.lower()=='norh':
							os.system('mv '+fil+' '+fil+'.fits')
							get_norh_imgs(sT,eT,outdir=os.getcwd())
						if ins.lower()=='goes':
							dat=ts.TimeSeries(fil)
							Obs=dat.observatory.replace('-','')
							Tab=dat.to_table()
							tres=str(24/len(Tab)*3600).split('.')[0]+'s'
							rng=dat.time.value[0].replace('-','').replace(' ','T').replace(':','').split('.')[0]+'-'+dat.time.value[0].replace('-','').replace(' ','T').replace(':','').split('.')[0]
							pickle.dump(Tab,open(Obs+'_'+rng+'_'+tres+'.p','wb'))
							Tab=Tab.to_pandas()
							Tab.to_csv(Obs+'_'+rng+'_'+tres+'.csv',index=False)				
					os.chdir(bdir)
			continue
		else:
			try:
				res=Fido.search(a.Time(sT.replace('T',' '),eT.replace('T',' ')),a.Instrument(ins))
				print('Data found for ',ins,'\nWavelength: ',wv)
			except:
				bads+=[ins]
				continue
			prv=res.keys()[0]

			#pickle.dump(res,open(ins+'_search_results.p','wb'))	
			if downdata==True:
				print('Downloading.....',ins,' data for ',wv)
				print(res[prv])
				dwnfils=Fido.fetch(res[prv],path='./')
	
				for fil in dwnfils.data:
					if ins.lower()=='norh':
						os.system('mv '+fil+' '+fil+'.fits')
						get_norh_imgs(sT,eT,outdir=os.getcwd())
					if ins.lower()=='goes':
						dat=ts.TimeSeries(fil)
						Obs=dat.observatory.replace('-','')
						Tab=dat.to_table()
						tres=str(24/len(Tab)*3600).split('.')[0]+'s'
						rng=dat.time.value[0].replace('-','').replace(' ','T').replace(':','').split('.')[0]+'-'+dat.time.value[0].replace('-','').replace(' ','T').replace(':','').split('.')[0]
						pickle.dump(Tab,open(Obs+'_'+rng+'_'+tres+'.p','wb'))
						Tab.write(Obs+'_'+rng+'_'+tres+'.csv',format='ascii.csv')				
			else:
				print('Data available: ',res[prv])
	return bads
	
########### Convert and plot radio spec #######################	
def radspec_to_pickle(filname,outdir='./'):
	'''
	Data of any RSTN telescope downloaded by sunpy or any .srs or callisto .fit standard file format should work.
	'''
	spec=Spectrogram(filname)
	stt=str(spec.start_time.value).split('.')[0].replace('-','').replace(':','')
	edt=str(spec.end_time.value).split('.')[0].replace('-','').replace(':','')
	specfil=spec.observatory+'_'+'_'.join([stt,edt])+'.p'
	pickle.dump({'DS':spec.data,'frequencies':spec.frequencies,'times':spec.times, 'start_time':spec.start_time, 'end_time':spec.end_time,'observatory':spec.observatory,'detector':spec.detector},open(outdir+specfil,'wb'))
	print('File: ',specfil,' written out!!')
	spec.plot()
	plt.savefig(outdir+specfil+'ng')
	plt.close()

#### Download and plot SXR lightcurve and its derivative

def plot_goes_derivative(sT, eT, satellite=None,savedata=True,outdir='./'):
	'''
	Download and plot GOES 1–8 Å flux and its time derivative.
	sT, eT: UTC datetime strings or datetime objects.
	satellite: GOES number; None selects the highest available number.
	savedata: True or False. To save data in the outdir
	'''
	nametag=sT.replace('-','').replace(':','')+'-'+eT.replace('-','').replace(':','')
	if isinstance(sT, str):
		sT = dt.fromisoformat(sT)
	if isinstance(eT, str):
		eT = dt.fromisoformat(eT)
	if eT <= sT:
		raise ValueError('eT must be later than sT.')

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

	downloaded = Fido.fetch(records)
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
		pickle.dump({'Time':times,'Flux':flux_values,'(Time,Derivative)':(derivative_times,derivative)},open(outdir+'GOES_data_deriv_'+nametag+'.p','rb'))

	#### Plot ####
	fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
	ax[0].plot(times, np.where(flux_values > 0, flux_values, np.nan),
		label='GOES-' + str(satellite) + ' 1–8 Å', color='tab:blue')
	ax[0].set_ylabel('Flux (W/m²)')
	ax[0].set_yscale('log')
	ax[0].legend()
	ax[0].grid(True, which='both')

	ax[1].plot(derivative_times, derivative, label='Derivative', color='tab:red')
	ax[1].set_ylabel('dFlux/dt (W/m²/s)')
	ax[1].set_xlabel('Time (UTC)')
	ax[1].legend()
	ax[1].grid(True)

	fig.tight_layout()
	plt.savefig(outdir+'GOES_data_deriv_'+nametag+'.pdf')
	plt.close()
	return
	

# Example usage:
# plot_goes_derivative("2025-08-10 00:00", "2025-08-10 02:00")

###### Get RSTN spectrograph selected dates ####
def get_RSTN(sT,eT,instrs='',callisto='',outdir='./',downdata=True):
	'''
	sT is the datetime string of the event 'YYYY-MM-DDThh:mm:ss'. LEARMONTH daily spectra will be downloaded for the day
	eT is the end datetime string of the event 'YYYY-MM-DDThh:mm:ss'.
	instrs = ['learmonth','sagamore'] - list of instruments of choice. or just a string 'EIT'
		Accepted list:{'eCALLISTO':'ecallisto','Holloman': 'holloman', 'Learmonth': 'learmonth', 'Palehua': 'palehua', 'Sagamore Hill': 'sagamore', 
		'San Vito': 'san-vito', 'holloman': 'Holloman', 'learmonth': 'Learmonth', 'palehua': 'Palehua', 'sagamore': 'Sagamore Hill', 'san-vito': 'San Vito'}

	callisto= Observatory choice for eCallisto data. Optional.
	outdir is where plot and daily spectrograph will be saved.
	downdata == True will download data else just report
	'''
	####### Finding the spectrograph #####
	bads=[]
	cwd=os.getcwd()
	for i in range(len(instrs)):
		ins=instrs[i]
		print('Serching data in ',ins)
		if not os.path.isdir(outdir+ins):
			os.mkdir(outdir+ins)
		os.chdir(outdir+ins)
		if 'callisto' not in ins.lower():
			try:
				res=Fido.search(a.Time(sT.replace('T',' '),eT.replace('T',' ')),a.Instrument('RSTN'),net.Observatory(ins))
			except:
				bads+=[ins]
				continue
			if len(res)==0:
				print('No data found for ',sT,' - ',eT,' period\n',ins)
				continue
			print('Provider: ',res.keys())
			prv=res.keys()
			prv=prv[0]	
				
		else:
			try:
				if callisto=='':
					res=Fido.search(a.Time(sT.replace('T',' '),eT.replace('T',' ')),a.Instrument('ecallisto'))
				else:
					res=Fido.search(a.Time(sT.replace('T',' '),eT.replace('T',' ')),a.Instrument('ecallisto'),net.Observatory(callisto))
			except:
				bads+=[ins]
				continue
					
			if len(res)==0:
				print('No data found for ',sT,' - ',eT,' period\n',ins)
				continue

			prv=res['callisto']['Observatory'].data
			
			print('Provider: ',np.unique(prv))			
			if 'BIR' in prv:
				prv='BIR'
			elif 'OOTY' in prv:
				prv='OOTY'
			else:
				prv=prv[0]
				locs=np.where(prv==prv[0])[0]
				res=res['callisto'][locs]
				
		if downdata==True:
			dats=Fido.fetch(res[prv],path='./')
			for dat in dats.data:
				radspec_to_pickle(dat,outdir='./')
		else:
			print('Available data: ',res[prv])
				
		os.chdir(cwd)			
	return bads
#####################################################
## Get CDAW WIND/WAVES and STEREO A & B DS #########
def get_DHband_DS(dT,outdir='./',delta_t=30):
	'''
	dT is the datetime string of the event 'YYYY-MM-DDThh:mm:ss'.
	delta_t: Time in minutes. Dynamic spectrum png will have start time from dT-delta_t. 
	outdir is where plot and daily spectrograph will be saved.
 	Returns combined 3 panel plot of png images from Wind and STEREO spacecraft
	'''
	######## INPUT #############
	WAVES_png_url='https://cdaw.gsfc.nasa.gov/images/wind/waves/' # Webpage to get specific image png for the date of interest. (e.g., YYYY/MM/DD/YYYYMMDD_hhmmss_windwaves.png)
	WAVES_png_url_latest='https://cdaw.gsfc.nasa.gov/images/wind/waves_h2M/' # Webpage where post 2021 data resides for wind waves.
	SWAVES_png_url='https://cdaw.gsfc.nasa.gov/images/stereo/swaves/'
	SWAVES_newpng_url='https://cdaw.gsfc.nasa.gov/images/stereo/swaves_new/'
	##############
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
			delts=np.abs((tims-np.datetime64(dt.strptime(dtup[0]+'T'+dtup[1],"%Y/%m/%dT%H:%M:%S"))-np.timedelta64(delta_t,'m')).astype('timedelta64[s]').astype(float)/60)
			gtim=np.where(delts==np.min(delts))[0][0]
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
			delts=np.abs((timsA-np.datetime64(dt.strptime(dtup[0]+'T'+dtup[1],"%Y/%m/%dT%H:%M:%S"))-np.timedelta64(delta_t,'m')).astype('timedelta64[s]').astype(float)/60)
			gtimA=np.where(delts==np.min(delts))[0][0]
			urllib.request.urlretrieve(tmpl+'/'+nlinksA[gtimA],outdir+'Stereo/'+dtup[0].replace('/','')+'_'+dtup[1].replace(':','')+'_swavesA.png')
			stA_png=outdir+'Stereo/'+dtup[0].replace('/','')+'_'+dtup[1].replace(':','')+'_swavesA.png'
			print('\n########## Stereo A data found for ',dtup)

		if len(nlinksB)==0:
			print('\n################\n No STEREO B data for ',dtup,'\n################')
		else:
			print('\n-------------------\nUseful data links found for STEREO B !!!\n---------------------------')
			timsB=np.array(timsB)
			delts=np.abs((timsB-np.datetime64(dt.strptime(dtup[0]+'T'+dtup[1],"%Y/%m/%dT%H:%M:%S"))-np.timedelta64(delta_t,'m')).astype('timedelta64[s]').astype(float)/60)
			gtimB=np.where(delts==np.min(delts))[0][0]
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
		imn=Image.new('RGB',(im1.width,im1.height+im2.height))
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
		imn=Image.new('RGB',(im1.width,im1.height+im2.height+im3.height))
		imn.paste(im1,(0,0))
		imn.paste(im2,(0,im1.height))
		imn.paste(im3,(0,im1.height+im2.height))
		imn.save(outdir+'Combined/'+dtup[0].replace('/','')+'_'+dtup[1].replace(':','')+'_comb.png')		
		del imn
		print('\nCombined images made: ',png_lst,'\n-----------------')
	
