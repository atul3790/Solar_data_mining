'''
RUN the MWA_obscat_maker.py to gather the observations for each year and make separate tables for each year.

The code takes in tables for each year and group the observations into date periods and find 
HEK reported flares, eruptions (ER) , coronal jets (CJ) and emergent flux (EF) during each period.

A table is written out with each group matched to different HEK events.
'''

from sunpy.net import attrs as a, Fido
import numpy as np,glob,pickle
import pandas as pd,os

##### INPUT #######################################################
basedir='/Users/amohan/Documents/NASA/Data/All_MWA_obs/'
findir='HEK_events/'
#hek_evs=[a.hek.ER,a.hek.CJ,a.hek.EF]
###################################################################

os.chdir(basedir)
os.system('mkdir '+findir)
Tables=sorted(glob.glob('MWAobs*.p'))[:-1]

for Tabn in Tables:
	Tab=pickle.load(open(Tabn,'rb'))
	tims=Tab['Start time (UTC)'].to_numpy().astype('datetime64[m]')
	
	dtims=tims[1:]-tims[:-1]
	
	#### Finding groups of continuous date observations.
	grps=[]
	gst=0
	gse=0
	
	for i in range(len(dtims)):
		ti=dtims[i]
		tnt=0
		if ti>60:
			gse=i
			if gse-gst==0:
				tnt+=1
				if tnt==1:
					grps+=[(tims[gst],(tims[gst]+np.timedelta64(Tab['Duration (sec)'].iloc[gst],'s')).astype('datetime64[m]'))]
					tnt=0
				gst=gse+1
				i+=1
				continue			
			grps+=[(tims[gst],tims[gse])]
			gst=gse+1
			i+=1
	grps+=[(tims[gst],tims[i])]
	Tab=pd.DataFrame(grps,columns=['Start','End'])
	stts=Tab['Start'].to_numpy().astype('datetime64[m]')
	edts=Tab['End'].to_numpy().astype('datetime64[m]')
	FLs=[]
	ERs=[]
	CJs=[]
	EFs=[]
	for i in range(len(Tab)):
		timerange=a.Time(stts[i],edts[i])
		
		## Checking for flares
		res = Fido.search(timerange, a.hek.FL)
		flrs=np.array([])
		if len(res)>0:
			for R in res:
				if len(R)>0:
					tmp=np.unique(np.array(R['fl_goescls']))
					TMP=np.array([j if j !='' else 'NIL' for j in tmp])
					flrs=np.append(flrs,TMP)
				else:
					flrs=np.append(flrs,np.array(['NIL']))	
			FLs+=[','.join(flrs)]
		else:
			FLs+=['NIL']
		### Checking for jets	
		res = Fido.search(timerange, a.hek.CJ)
		if len(res)>0:
			CJs+=['YES']
		else:
			CJs+=['NIL']	

		### Checking for Eruptions	
		res = Fido.search(timerange, a.hek.ER)
		if len(res)>0:
			ERs+=['YES']
		else:
			ERs+=['NIL']	

		### Checking for Eruptions	
		res = Fido.search(timerange, a.hek.EF)
		if len(res)>0:
			EFs+=['YES']
		else:
			EFs+=['NIL']
		print(i)		
	Tab['Flare']=FLs
	Tab['ER']=ERs
	Tab['EF']=EFs
	Tab['CJ']=CJs
	Tab.to_csv(findir+'MWA-HEK_evperiod_'+Tabn.split('.p')[0].split('_')[1]+'.csv',index=False)
