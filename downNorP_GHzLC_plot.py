'''
down NoRP
'''

from astropy.io import fits
from datetime import datetime as dt, timezone
from matplotlib.dates import DateFormatter
import pandas as pd
import numpy as np,pickle,glob
#from multiprocessing import Pool
import urllib,os
from matplotlib import pyplot as plt

########### INPUT
basedir='C:\\Users\\amohan\\Documents\\SGRE_workshop\\SGRE_NORP\\'
Table='SGRE_cat.csv'  # Table should ahve columns Date, Time
issgre=True

basedir='C:\\Users\\amohan\\Documents\\SGRE_workshop\\Rest_DHtypeIVs\\'
Table='datetimes.csv'  # Table should ahve columns Date, Time
issgre=False

email='atul.3.7.1990@gmail.com'
###############################################################
'''
nproc=3
'''
#Compn=os.uname().nodename.split('.')[0]
Compn='Dell'

cwd=os.getcwd()
os.chdir(basedir)

tab=pd.read_csv(Table)
plt.ioff()
if issgre==True:
	date,time=tab['EventDate'].to_numpy(),tab['Time'].to_numpy()
else:
	date,time=tab['Date'].to_numpy(),tab['Time'].to_numpy()

sTs=[]
eTs=[]
for i in range(len(date)):
	if len(time[i])<8:
		sTs+=[str(np.datetime64(date[i].replace('/','-')+'T'+time[i]+':00')-np.timedelta64(1,'h'))]
		eTs+=[str(np.datetime64(date[i].replace('/','-')+'T'+time[i]+':00')+np.timedelta64(1,'h'))]
	else:
		sTs+=[str(np.datetime64(date[i].replace('/','-')+'T'+time[i])-np.timedelta64(1,'h'))]
		eTs+=[str(np.datetime64(date[i].replace('/','-')+'T'+time[i])+np.timedelta64(1,'h'))]		

def find_norp(sT,eT,outdir='./',plot=True):
	'''
	sT is the datetime string of the event 'YYYY-MM-DDThh:mm:ss'.
	eT is the end datetime string of the event 'YYYY-MM-DDThh:mm:ss'.
	outdir is the path to directory where to store the light curve and quick look image data files
	plot =True ensures the data will be plotted for the sT to eT period
	'''
	dtFmt=DateFormatter('%b-%d %H:%M')
	baseurl='https://solar.nro.nao.ac.jp/norp/fits/'
	sT=np.datetime64(sT)
	eT=np.datetime64(eT)
	dates=np.arange(sT,eT,np.timedelta64(24,'h'))
	chkdts=[str(dates[0]-np.timedelta64(24,'h')).split('T')[0].replace('-','/')]+[str(i).split('T')[0].replace('-','/') for i in dates]+[str(dates[-1]+np.timedelta64(24,'h')).split('T')[0].replace('-','/')]		
	urls=[baseurl+i[:-2]+'norp'+i[2:].replace('/','')+'.fits.gz' for i in chkdts]
	## Expected frequencies
	subdir=str(sT+(eT-sT)/2).replace(':','')+'/'
	Flrtim=str(sT+(eT-sT)/2).replace(':','')
	print('###################################################\n')
	print('Searching for NoRP data in ',str(sT),' - ',str(eT))
	print('Flare time: ',Flrtim)
	print('###################################################\n')

	if not os.path.isdir(outdir+subdir):
		os.mkdir(outdir+subdir)
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
			os.system('rm -rf '+fil)
			break
		if End<sT: # If Obs end date time < Start epoch requested then check the next file
			os.system('rm -rf '+fil)
			continue	
		if (beg<sT and End<eT) or (beg<sT and End>eT) or (beg>sT and End>eT) or (beg>sT and End<eT):
			print('Data found at ',url,';   Flare:',Flrtim)
			data=fits.getdata(fil)
			Time=data['Time'][0]		
			tmp=np.array(list(head.keys()))
			tfqs=tmp[[True if 'FREQ' in i else False for i in tmp]]
			Frqs=[str(np.round(float(head[i])*10**-9,1)) for i in tfqs]
			Alltags=['_'+str(int(np.round(float(i),0)))+'GHz' for i in Frqs]
			frq_dict=dict(zip(Frqs,Alltags))

			idlst = dt(1979, 1, 1, hour=0, minute=0, second=0, microsecond=0, tzinfo=timezone.utc)
			times=[dt.fromtimestamp(idlst.timestamp()+i,timezone.utc) for i in Time]
			times=np.array(times).astype(np.datetime64)
			t1=sT if sT>times[0] else times[0]
			t2=eT if times[-1]>eT else times[-1]

			if os.path.isfile(outdir+subdir+'MaxFlux_'+str(t1).replace(':','')+'-'+str(t2).replace(':','')+'.csv'):
				rev_d=dict(zip(Alltags,Frqs))
				print('Analysis completed for the event: ',Flrtim)
				## Check if all plots for the period are there
				if plot==True:
					getall=len(glob.glob(outdir+subdir+'Flux*_'+str(t1).replace(':','').replace('-','')+'-'+str(t2).replace(':','').replace('-','')+'.png'))
					if len(getall)==len(Frqs):
						print('All plots are in place for flare within',sT,' - ' ,eT,' So doing next time slot!! ')
						os.system('rm -rf '+fil)
						continue
					elif len(getall)>0:
						print('Not all frequency LC made for ',t1,' - ',t2,' range.')
						btags=['_'+str(int(np.round(float(i),0)))+'GHz' for i in Frqs]
						gtags=[i.split('GHz_')[0].replace('Flux_','')+'GHz_' for i in getall]
						Alltags=list(set(btags)-set(gtags))
						Frqs=[rev_d[i] for i in Alltags]
						frq_dict=dict(zip(Frqs,Alltags))
						print('Undone frequency range: ',frq_dict)
																			 							
			l1,l2=np.where(times<=t1)[0][-1],np.where(times>=t2)[0][0]
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
				Flx_P=Flx_V/Flx_I
				
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
						Idur=(FI[-1]-FI[0]).astype(float)/60*10**-6
						TI=np.nansum(Flx_I[l1:l2+1][Flocs])*DelT # Estimating total flare I flux
						if tog_V==True:
							TV=np.nansum(Flx_V[l1:l2+1][Flocs])*DelT # Estimating total flare V flux										
				else:
					IMxt+=['---']
					nIm+=[np.nan]
				if tog_V==True:
					y1,y2=np.round(np.nanmin(Flx_V[l1:l2+1])*0.9)-1,np.round(Vm*1.1)+1
					Vpls=np.where(Flx_V[l1:l2+1]==Vm)[0]
					VMxt+=[str(times[l1:l2+1][Vpls[0]]).replace('T',' ')]
					nVm+=[len(Vpls)]
				else:
					VMxt+=['---']
					nVm+=[np.nan]
				TotI+=[TI]
				TotV+=[TV]
				Durs+=[Idur]					
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
				
				if tog_V==True:
					Vpls=np.where(Flx_V[l1:l2+1]==Vm)[0]
					VMxt+=[str(times[l1:l2+1][Vpls[0]]).replace('T',' ')]
				pickle.dump({'Time':times,'Flux_I':Flx_I,'Flux_V':Flx_V,'Flux V/I':Flx_P,'Start time':t1,'End time':t2,'I Flux range':[Iy1,Iy2],'V Flux range':[y1,y2]},open(outdir+subdir+'FluxLC'+tag+'_'+str(t1).replace(':','')+'-'+str(t2).replace(':','')+'.p','wb'))			

			os.system('mv '+fil+' '+outdir+subdir)
			IFlxs=np.array(IFlxs)
			VFlxs=np.array(VFlxs)
			DF=pd.DataFrame({'Frequency (GHz)':Frqs,'Median_I_flux (SFU)':Imed,'Median_V_flux (SFU)':Vmed,'Max_I_flux (SFU)':IFlxs,'Max_V_flux (SFU)':VFlxs,'Dur (min)':Durs,'Max Pol':VFlxs/IFlxs,'Max_to_median_I (SFU)':IFlxs/np.array(Imed),'Max_to_median_V (SFU)':VFlxs/np.array(Vmed),'Max_I_time (UT)':IMxt,'Max_V_time (UT)':IMxt,'Max_I_repetition':nIm,'Max_V_repetition':nVm,'Flare_Energy_I (SFU s)':TotI,'Flare_Energy_V (SFU s)':TotV,'Mean pol':np.array(TotV)/np.array(TotI)})
			DF=DF.astype(str)
			DF.to_excel(outdir+subdir+'MaxFlux_'+str(t1).replace(':','')+'-'+str(t2).replace(':','')+'.xlsx',index=False)
			DF.to_csv(outdir+subdir+'MaxFlux_'+str(t1).replace(':','')+'-'+str(t2).replace(':','')+'.csv',index=False)
			pickle.dump(DF,open(outdir+subdir+'MaxFlux_'+str(t1).replace(':','')+'-'+str(t2).replace(':','')+'.p','wb'))
	if len(badurls)>0:
		pickle.dump(badurls,open(outdir+subdir+'Bad_urls_'+str(Flrtim)+'.p','wb'))

for i in range(len(sTs)):
	find_norp(sTs[i],eTs[i])
'''	
p=Pool(nproc)
T1=[sTs[8]]
T2=[eTs[8]]
succ=p.starmap(find_norp,zip(T1,T2))
'''
os.system("echo 'NoRP data downloads done for "+Table+".\nMachine:  "+Compn+"'| mail -s 'Job done' "+email)
os.chdir(cwd)