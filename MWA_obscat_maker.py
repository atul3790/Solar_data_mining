'''
	This code helps to obtain MWA observation metadata falling in a user-specified year interval. 
  It then tabulates the observation details with date, duration, observer, and observing mode in a CSV table and 
  a binary pickle file.
	If there exist any previous catalogs for a year range in the output repository, the code will 
	update it to cover the requested full range
	
	!!!!!! BEWARE !!!!!!!!!!
	MWA server thinks you are a scammer if this code is run more than once within a period
	So if it fails with an HTTP 504 Gateway error its not a bug, just run the code tomorrow.
	!!!!!!!!!!!!!!!!
'''

from datetime import datetime,timedelta
import pandas as pd,glob,pickle,os
from time_conv import utc_gps
from datetime import datetime as dt,timedelta
############################ Input Parameters ########################
basedir='/Users/amohan/Documents/NASA/' # Working directory within which data directories exist
outdir='Data/All_MWA_obs/'# This is the directory where all observation details and catalogues are saved.
pagesize=27000# This is the maximum number of observation metadata returned by the querry at a time.
t0=2017 # Start year as int
tf=2024 # End year as int
#########################################################################

BASEURL = 'http://ws.mwatelescope.org/metadata/' # MWA metadata server

cwd=os.getcwd()
os.chdir(basedir)

if not os.path.isdir(outdir):
	os.mkdir(outdir)
## If already earlier analysis results exist in the repository the code finds the last year 
#for which the metadata is found. Then the start year is reset to the year from which there is no metadata file available
	
DT=t0
flag=0
cflag=0
while DT<=tf:
	if DT>2021:
		cflag=1
	t1=str(DT)+'-01-01T00:00:00'
	if DT+3<tf:
		t2=str(DT+3)+'-12-31T23:59:59'
	else:
		t2=str(tf)+'-12-31T23:59:59'
		flag=1
		
	tag=t1.split('-')[0]+'-'+t2.split('-')[0]
	print('Getting data within: ',tag)	
	minid=utc_gps(t1) # Leapseconds give the utc time + added seconds to correct for leapseconds. subtract it from GPS start date and take number of seconds to get the gps seconds time. 20 seconds are just subtracted so that you don't miss any observation around the date needed.	
	maxid=utc_gps(t2)
	
	print ('Starting to get the MWA Obs for the Solar Project....')
	if clfag==0:
		res=pd.read_html(BASEURL+'find/?search=search&html=1&projectid=G0002&pagesize='+str(pagesize)+'&mode=HW_LFILES&mintime='+str(minid)+'&maxtime='+str(maxid))
	else:
		res=pd.read_html(BASEURL+'find/?search=search&html=1&projectid=G0002&pagesize='+str(pagesize)+'&mintime='+str(minid)+'&maxtime='+str(maxid))
	print('Found observations between ',tag)
	
	## Removing calibration observations ##########################
	bads=np.array([cl for cl in range(len(res[0])) if 'sun' not in res[0]['Obs Name'][cl].lower()])
	Tab=res[0].drop(bads)
	Tab.reset_index(drop=True)
	#################################################################
	pickle.dump(Tab,open(outdir+'MWAobs_'+tag+'.p','wb'))
	Tab.to_csv(outdir+'MWAobs_'+tag+'.csv',index=False)	
	if DT+4<=tf:
		DT=DT+4
	else:
		DT=tf	
	if flag==1:
		print('Exit flag reached!! All years requested are covered!!')
		break

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
		
