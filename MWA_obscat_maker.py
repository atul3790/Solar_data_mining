'''
	This code helps to obtain MWA observation metadata for each year falling within user-specified start and end year.
	Finally it will build a master catalog for the whole duraiton.
	### NOTE that the end year is also included. So ensure that there is data in the end year as well.
 
	!!!!!! BEWARE !!!!!!!!!!
	MWA server thinks you are a scammer if this code is run more than once within a period
	So if it fails with an HTTP 504 Gateway error its not a bug, just run the code tomorrow.
	!!!!!!!!!!!!!!!!
'''

from datetime import datetime,timedelta
import pandas as pd,glob,pickle,os,numpy as np
from time_conv import utc_gps
from datetime import datetime as dt,timedelta
############################ Input Parameters ########################
basedir='/Users/amohan/Documents/NASA/' # Working directory within which data directories exist
outdir='Data/All_MWA_obs/'# This is the directory where all observation details and catalogues are saved.
pagesize=27000# This is the maximum number of observation metadata returned by the querry at a time.
t0=2018 # Start year as int
tf=2024 # End year as int
buffer_days=4 # the number of days before today, till when the code can search for events.
			# if buffer_days=4. Code will search for events till today - 4 days.
########################################################################
BASEURL = 'http://ws.mwatelescope.org/metadata/' # MWA metadata server

cwd=os.getcwd()
os.chdir(basedir)

if not os.path.isdir(outdir):
	os.mkdir(outdir)
## If already earlier analysis results exist in the repository the code finds the last year 
#for which the metadata is found. Then the start year is reset to the year from which there is no metadata file available
'''
if mode==1:
	all_exist=sorted(glob.glob(outdir+'MWAobs_*.p'))
	yrs=[]
	if len(all_exist)>0:
		for i in all_exist:
			yrs+=[int(i.split('MWAobs_')[1].replace('.p','').split('-')[1])]
		t0=np.max(yrs)+1
	##########################
'''	
DT=t0
cflag=0
while DT<=tf:
	if DT>2021:
		cflag=1

	t1=str(DT)+'-01-01T00:00:00'
	t2=str(DT)+'-12-31T23:59:59'
	if np.datetime64(t2)>np.datetime64('now'):
		t2=str(np.datetime64('now')-np.timedelta64(buffer_days,'D'))
		
	tag=t1.split('-')[0]+'-'+t2.split('-')[0]
	print('Getting data within: ',tag)	
	'''
	minid=utc_gps(t1) # Leapseconds give the utc time + added seconds to correct for leapseconds. subtract it from GPS start date and take number of seconds to get the gps seconds time. 20 seconds are just subtracted so that you don't miss any observation around the date needed.	
	maxid=utc_gps(t2)
	'''
	print ('Starting to get the MWA Obs for the Solar Project....')
	ST='%3A'.join(t1.split(':')[:-1])
	ET='%3A'.join(t2.split(':')[:-1])
	#res=pd.read_html(BASEURL+'find/?search=search&html=1&projectid=G0002&pagesize='+str(pagesize)+'&mode=HW_LFILES&mintime='+str(minid)+'&maxtime='+str(maxid))
	if cflag==0:
		res=pd.read_html(BASEURL+'find/?search=search&html=1&projectid=G0002&dataquality=1&notdeleted=on&pagesize='+str(pagesize)+'&mode=HW_LFILES&mintime_utc='+ST+'&maxtime_utc='+ET+'&page=1')
	else:
		res=pd.read_html(BASEURL+'find/?search=search&html=1&projectid=G0002&dataquality=1&notdeleted=on&pagesize='+str(pagesize)+'&mintime_utc='+ST+'&maxtime_utc='+ET+'&page=1')
	
	print('Found observations between ',tag)
	bads=np.array([cl for cl in range(len(res[0])) if 'sun' not in res[0]['Obs Name'][cl].lower()])
	Tab=res[0].drop(bads)
	Tab.reset_index(drop=True)
	pickle.dump(Tab,open(outdir+'MWAobs_'+tag+'_V1.p','wb'))
	Tab.to_csv(outdir+'MWAobs_'+tag+'_V1.csv',index=False)
	DT+=1

all_Tabs=sorted(glob.glob(outdir+'MWAobs_*.p'))

r1=[]
for i in all_Tabs:
	r1+=i.split('MWAobs_')[1].replace('.p','').split('-')
r1=np.array(r1).astype(int)
tag=str(np.min(r1))+'-'+str(np.max(r1))

Tab=pickle.load(open(all_Tabs[0],'rb'))
for i in np.arange(len(all_Tabs)-1)+1:
	Tab=pd.concat([Tab,pickle.load(open(all_Tabs[i],'rb'))])
Tab.to_csv(outdir+'MWAobs_master_'+tag+'.csv',index=False)

os.chdir(cwd)
