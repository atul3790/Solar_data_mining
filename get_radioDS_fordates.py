'''
Get radio DS or cross match 
'''
from astropy.table import Table,Column,Row
import numpy as np,os,glob
import urllib
############### INPUT ########################
basedir='C:\\Users\\amohan\\Documents\\NASA\\Data\\CDAW_typeIV_database\\'
outdir='Metric_typeIV_counterparts\\metricDS\\'
lea_url='https://www.sws.bom.gov.au/Category/World%20Data%20Centre/Data%20Display%20and%20Download/Spectrograph/station/learmonth/images/'
cul_url='https://www.sws.bom.gov.au/Category/World%20Data%20Centre/Data%20Display%20and%20Download/Spectrograph/station/culgoora/images/'
art_url='http://artemis-iv.phys.uoa.gr/QuickLooks/'
nda_url='https://rsdb.obs-nancay.fr/QL/Nda/'
orfees_url='https://rsdb.obs-nancay.fr/QL/Orfees/'
hiras_url='https://sunbase.nict.go.jp/solar/denpa/hirasDB/'
blen_url='http://soleil.i4ds.ch/solarradio/qkl/'
datetim_fil='DHtypeIV_datetimes.csv'
#############################################

## nda format: NdaP1F20230801.png; lea/cul: 23/20230801spectrograph.gif
## orfees: orfees20131006.png; hiras:2016/02/160201a.gif
## blen: 2004/07/23/BLEN5M_20040723_174500_L2.fit.gz.png & L1

cwd=os.getcwd()
os.chdir(basedir)

Tab=Table.read(datetim_fil,format='ascii.csv')
time_dict={'0-8':['LEA','HIR','CUL'],'8-9.5':['HIR','LEA','CUL'],'8-18':['ART','NDA','ORF'],'14-18':['BLE'],'19-24':['LEA','CUL','HIR']}
url_dict=dict(zip(['LEA','HIR','CUL','ART','NDA','ORF','BLE'],[lea_url,hiras_url,cul_url,art_url,nda_url,orfees_url,blen_url]))

tag=dict(zip(['NDA','ORF'],['NdaP1F','orfees']))
for i in Tab:
	day1=i['Date'].replace('/','')
	day2=np.datetime64(i['Date'].replace('/','-')+'T01:00:00')+np.timedelta64(24,'h')
	day2=str(day2).split('T')[0].replace('-','')
	day=day1
	try:
		os.mkdir(outdir+day1+'\\')
	except:
		print(outdir+day,' already exists.')
	tim=i['Time'].split(':')
	tim=float(tim[0])+float(tim[1])/60.
	obs=[]
	for j in list(time_dict.keys()):
		s,e=j.split('-')
		if tim>float(s) and tim<float(e):
			obs+=time_dict[j]
			print('Getting DS from : ',obs)
	for O in obs:
		day=day1
		Url=url_dict[O]
		if O =='BLE':
			hrs=i['Time'].split(':')[0]
			mins=float(i['Time'].split(':')[1])
			if mins>=0 and mins<15:
				mins='00'
			elif mins>=15 and mins<30:
				mins='15'
			elif mins>=30 and mins<45:
				mins='30'
			else:
				mins='45'
			timtag=hrs+mins+'00'													
			try:
				urllib.request.urlretrieve(Url+day[0:4]+'/'+day[4:6]+'/'+day[6:]+'/BLEN5M_'+day+'_'+timtag+'_L1.fit.gz.png',outdir+day1+'\\'+i['Time'].replace(':','')+O+'_L1.png')
				print('BLEIN L1 data found for ',day)
			except:
				print('No BLEIN data for ',day)			
			try:
				urllib.request.urlretrieve(Url+day[0:4]+'/'+day[4:6]+'/'+day[6:]+'/BLEN5M_'+day+'_'+timtag+'_L1.fit.gz.png',outdir+day1+'\\'+i['Time'].replace(':','')+O+'_L2.png')
				print('BLEIN L2 data found for ',day)
			except:
				print('No BLEIN data for ',day)			
		if O in ['LEA','CUL']:
			try:
				if tim>10:
					day=day2			
				urllib.request.urlretrieve(Url+day[2:4]+'/'+day+'spectrograph.gif',outdir+day1+'\\'+i['Time'].replace(':','')+O+'.png')
				print(O,' data found for ',i['Date'])			
			except:
				print('No ',O,' data for ',i['Date'])
		elif O in ['NDA','ORF']:
			try:		
				urllib.request.urlretrieve(Url+tag[O]+day+'.png',outdir+day1+'\\'+i['Time'].replace(':','')+O+'.png')
				print(O,' data found for ',i['Date'])			
			except:
				print('No ',O,' data for ',i['Date'])
		elif O == 'ART':
			try:
				urllib.request.urlretrieve(Url+day[0:4]+'/'+day[4:6]+'/'+day+'_60s.Int01.png',outdir+day1+'\\'+i['Time'].replace(':','')+O+'.png')
				print(O,' data found for ',i['Date'])			
			except:
				print('No ',O,' data for ',i['Date'])						
		elif O == 'HIR':
			inx='b'
			if tim>20 and tim<24:
				day=day2
				inx='a'	
			if tim<3:
				inx='a'	
			try:
				urllib.request.urlretrieve(Url+day[0:4]+'/'+day[4:6]+'/'+day[2:]+inx+'.gif',outdir+day1+'\\'+i['Time'].replace(':','')+O+'_'+inx+'.png')
				#urllib.request.urlretrieve(Url+day[0:4]+'/'+day[4:6]+'/'+day[2:]+'b.gif',outdir+day+'\\'+i['Time'].replace(':','')+O+'_b.png')
				print(O,' data found for ',i['Date'])			
			except:
				print('No ',O,' data for ',i['Date'])
	
os.chdir(cwd)	