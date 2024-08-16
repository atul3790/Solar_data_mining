'''
gps utc and viceversa conversion. This module has 3 functions
1. gps_utc : Converts a given GPS time (given as an integer) to UTC
2. utc_gps : Converts a given UTC time string ('YYYY-MM-DDThh:mm:ss') to Gps seconds
3. leap_sec: Estimates leap seconds for a given UTC time stamp ('YYYY-MM-DDThh:mm:ss')
'''
from datetime import datetime as dt,timedelta

import numpy as np


#### Leap seconds computation
def leap_sec(sT):
	'''
	sT: Time in 'YYYY-MM-DDThh:mm:ss' format
	Provides leapseconds for any given UTC time.
	Works only for dates> 1980 Jan 1.
	'''
	Conv=19
	sT=np.datetime64(sT).astype(dt)
	epochs=[dt(1981,7,1),dt(1982,7,1),dt(1983,7,1),dt(1985,7,1),
	dt(1988,1,1),dt(1990,1,1),dt(1991,1,1),dt(1992,7,1),dt(1993,7,1),dt(1994,7,1),
	dt(1996,1,1),dt(1997,7,1),dt(1999,1,1),dt(2006,1,1),dt(2009,1,1),dt(2012,7,1),
	dt(2015,7,1),dt(2017,1,1)]
	
	for epc in epochs:
		if sT>=epc:
			Conv+=1
		else:
			break
	return Conv							

##### GPS to UTC time conversion	
def gps_utc(Gps):
	'''
	Gps: GPS time stamp
	Returns UTC equivalent time.
	Works only for dates> 1980 Jan 1.
	'''
	tst=dt(1980, 1, 6) + timedelta(seconds=Gps - (40 - 19))
	Conv=leap_sec(str(tst).split('.')[0].replace(' ','T'))			
	return dt(1980, 1, 6) + timedelta(seconds=Gps - (Conv - 19))

##### UTC time to GPS conversion	
def utc_gps(sT):
	'''
	sT: Time in 'YYYY-MM-DDThh:mm:ss' format
	Provides leapseconds for any given UTC time.
	Works only for dates> 1980 Jan 1.
	'''	
	Conv=leap_sec(sT)
	sT=np.datetime64(sT).astype(dt)
	return (sT-dt(1980,1,6)).days*86400+(sT-dt(1980,1,6)).seconds+(Conv-19)
