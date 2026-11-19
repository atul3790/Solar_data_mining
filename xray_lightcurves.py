"""
xray_lightcurves.py -- disk-integrated X-ray light-curve retrieval and overlays.

Companion to Multiwb_datamine.py (uses its get_STIX_LC for Solar Orbiter/STIX).
All times are UTC. Functions return data objects / output paths and only write
to the caller-selected folders.

get_GOES_XRS(sT, eT, outdir='./GOES', satellite=None, resolution='avg1m')
	Input: interval, optional GOES number (None = highest available), XRS product
	('avg1m' 1-min averages or 'flx1s' 1-s fluxes).
	Output: outdir/GOES<N>_XRS_<res>_<sT>-<eT>.csv. Returns (DataFrame, csv path, satellite).

plot_STIX_GOES(stix_df, goes_df, outpath, bands=None, stix_rebin=None, goes_channel='xrsb', ...)
	Input: STIX DataFrame from Multiwb_datamine.get_STIX_LC, GOES DataFrame from
	get_GOES_XRS (or any DataFrame with 'datetime' and flux columns in W m^-2).
	Output: figure at outpath (.png/.pdf). Returns outpath.
"""
import os
import numpy as np, pandas as pd
import matplotlib.pyplot as plt, matplotlib.dates as mdates
from astropy import units as u

STIX_COLORS=['#2a78d6','#1baf7a','#eda100','#e34948','#6250d6']	# as in get_STIX_LC
GOES_LABEL={'xrsb':'1–8 Å','xrsa':'0.5–4 Å'}


def _utc_naive(t):
	'''Datetime-like Series/Index -> tz-naive UTC pandas datetimes.'''
	t=pd.to_datetime(t)
	if getattr(t,'dt',None) is not None:
		return t.dt.tz_convert('UTC').dt.tz_localize(None) if t.dt.tz is not None else t
	return t.tz_convert('UTC').tz_localize(None) if t.tz is not None else t


def get_GOES_XRS(sT,eT,outdir='./GOES',satellite=None,resolution='avg1m'):
	'''
	Download GOES XRS fluxes with SunPy Fido (science-ready NOAA L2 product; no recalibration).
	sT,eT: UTC strings/datetimes. satellite: GOES number or None (highest available).
	resolution: 'avg1m' or 'flx1s'.
	Returns (DataFrame[datetime, xrsa, xrsb, xrsa_quality, xrsb_quality] in W m^-2, csv path, satellite).
	'''
	from sunpy.net import Fido, attrs as a
	import sunpy.timeseries as ts
	sT,eT=pd.Timestamp(sT),pd.Timestamp(eT)
	if eT<=sT:raise ValueError('eT must be later than sT.')
	os.makedirs(outdir,exist_ok=True)
	q=[a.Time(sT,eT),a.Instrument('XRS'),a.Resolution(resolution)]
	if satellite is not None:q.append(a.goes.SatelliteNumber(int(satellite)))
	res=Fido.search(*q)
	if res.file_num==0:raise ValueError('No GOES XRS data found for this interval.')
	rec=res['xrs']
	if satellite is None:
		satellite=max(int(v) for v in rec['SatelliteNumber'])
		rec=rec[np.asarray(rec['SatelliteNumber'],dtype=int)==satellite]
	files=Fido.fetch(rec,path=os.path.join(outdir,'raw'))
	if files.errors or len(files)==0:raise RuntimeError('GOES download failed: '+str(files.errors))
	gts=ts.TimeSeries(files,source='XRS',concatenate=True).truncate(sT,eT)
	df=gts.to_dataframe().sort_index();df=df[~df.index.duplicated(keep='first')]
	for c in ('xrsa','xrsb'):df[c]=df[c].to_numpy(float)*gts.units[c].to(u.W/u.m**2)
	keep=[c for c in ('xrsa','xrsb','xrsa_quality','xrsb_quality') if c in df.columns]
	df=df[keep].rename_axis('datetime').reset_index();df['datetime']=_utc_naive(df['datetime'])
	tag=lambda t:t.strftime('%Y%m%dT%H%M%S')
	csv=os.path.join(outdir,f'GOES{satellite}_XRS_{resolution}_{tag(sT)}-{tag(eT)}.csv')
	df.to_csv(csv,index=False)
	print(f'GOES-{satellite} XRS {resolution}: {len(df)} samples -> {csv}')
	return df,csv,satellite


def plot_STIX_GOES(stix_df,goes_df,outpath,bands=None,stix_rebin=None,goes_channel='xrsb',
		goes_label='GOES',stix_ylim=None,goes_ylim=(1e-7,1e-3),xlim=None,time_label='Time (UT)',
		title=None,figsize=(11,5.6),fontsize=15,dpi=200):
	'''
	STIX quick-look counts (left log axis) with GOES XRS flux (right log axis) on a shared time axis.
	stix_df: DataFrame with 'datetime' and '<lo>-<hi> keV' columns (counts/s), e.g. get_STIX_LC output.
	goes_df: DataFrame with 'datetime' and goes_channel in W m^-2.
	bands: list of STIX band column names (None = all keV columns).
	stix_rebin: averaging window in seconds for STIX (None = native 4 s).
	stix_ylim/goes_ylim/xlim: axis limits (None = automatic); xlim as (start, end) strings.
	time_label: x-axis label; say whether STIX times are Earth light-arrival corrected.
	Returns outpath. Inputs are not modified.
	'''
	S=stix_df.copy();S['datetime']=_utc_naive(S['datetime']);S=S.set_index('datetime').sort_index()
	G=goes_df.copy();G['datetime']=_utc_naive(G['datetime']);G=G.set_index('datetime').sort_index()
	bands=bands or [c for c in S.columns if c.endswith('keV')]
	if stix_rebin:S=S[bands].resample(f'{int(stix_rebin)}s').mean().dropna(how='all')
	fig,ax=plt.subplots(figsize=figsize)
	for c,col in zip(bands,STIX_COLORS*2):
		ax.plot(S.index,S[c].where(S[c]>0),color=col,lw=1.1,drawstyle='steps-mid',label='STIX '+c)
	ax.set_yscale('log');ax.set_ylabel(r'STIX counts s$^{-1}$',fontsize=fontsize)
	if stix_ylim:ax.set_ylim(*stix_ylim)
	ax2=ax.twinx();g=G[goes_channel].where(G[goes_channel]>0)
	ax2.plot(G.index,g,color='k',lw=1.8,label=f'{goes_label} {GOES_LABEL.get(goes_channel,goes_channel)}')
	ax2.set_yscale('log');ax2.set_ylabel(r'GOES '+GOES_LABEL.get(goes_channel,'')+r' flux (W m$^{-2}$)',fontsize=fontsize)
	if goes_ylim:ax2.set_ylim(*goes_ylim)
	if goes_channel=='xrsb':
		for lev,lab in zip([1e-6,1e-5,1e-4],'CMX'):
			ax2.axhline(lev,color='0.6',lw=0.6,ls=':')
			ax2.text(0.993,lev*1.12,lab,transform=ax2.get_yaxis_transform(),va='bottom',ha='right',fontsize=fontsize-2,color='0.35',fontweight='bold')
	if xlim:ax.set_xlim(pd.Timestamp(xlim[0]),pd.Timestamp(xlim[1]))
	ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
	ax.set_xlabel(time_label,fontsize=fontsize);ax.tick_params(labelsize=fontsize-2);ax2.tick_params(labelsize=fontsize-2)
	ax.grid(alpha=0.25)
	h1,l1=ax.get_legend_handles_labels();h2,l2=ax2.get_legend_handles_labels()
	ax.legend(h1+h2,l1+l2,fontsize=fontsize-4.5,ncol=3,loc='upper left',framealpha=0.85)
	if title:ax.set_title(title,fontsize=fontsize-1)
	fig.tight_layout();os.makedirs(os.path.dirname(os.path.abspath(outpath)),exist_ok=True)
	fig.savefig(outpath,dpi=dpi);plt.close(fig)
	print('STIX + GOES overlay ->',outpath)
	return outpath
