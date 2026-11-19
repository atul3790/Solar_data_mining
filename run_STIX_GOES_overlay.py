"""
run_STIX_GOES_overlay.py -- task script: Solar Orbiter/STIX quick-look light curves with
GOES XRS overplotted, for one event interval.

Uses Multiwb_datamine.get_STIX_LC (STIX, from ESA SOAR) and xray_lightcurves
(GOES via SunPy Fido + overlay plot). Edit the INPUT block and run:  python run_STIX_GOES_overlay.py
"""

# =========================================================================================
# PARAMETER GUIDE
# -----------------------------------------------------------------------------------------
# ---- Event and paths ----
# sT, eT          : interval start / end, UTC ISO strings ('YYYY-MM-DDTHH:MM:SS').
# module_path     : folder containing Multiwb_datamine.py and xray_lightcurves.py.
# output_folder   : results folder; STIX/ and GOES/ subfolders hold downloads and CSVs,
#                   the overlay figure is written at its top level. Created if missing.
#
# ---- Sampling ----
# stix_rebin      : STIX averaging window [s]. None = native 4 s quick-look bins.
#                   60 matches GOES 1-min averages and suppresses count noise.
# goes_res_option : 1 = 1-min averages ('avg1m'), 2 = 1-s fluxes ('flx1s', large files).
#
# ---- Instrument settings ----
# earth_time      : True shifts STIX times by the file's EAR_TDEL (Sun->Earth minus
#                   Sun->Solar Orbiter light time) so both curves are in Earth
#                   light-arrival time. False keeps STIX spacecraft UTC (then the
#                   two time axes are NOT directly comparable).
# goes_satellite  : GOES number (e.g. 18, 19). None = highest number available.
# goes_channel_option : 1 = 1-8 A long channel ('xrsb', flare class), 2 = 0.5-4 A ('xrsa').
# stix_bands      : list of STIX band names to plot, e.g. ['4-10 keV','10-15 keV'].
#                   None = all five quick-look bands (4-10 ... 50-84 keV).
#
# ---- Output ----
# stix_ylim       : (min, max) STIX counts/s axis; None = automatic.
# goes_ylim       : (min, max) GOES flux axis [W m^-2]; None = automatic. (1e-7,1e-3)
#                   spans A to X class.
# fig_format      : 'png' or 'pdf'.
# out_tag         : plain name tag for the figure; '' = built from sT/eT.
# =========================================================================================

########### INPUT
# ---- Event and paths ----
sT='2025-05-30T07:00:00'	# UTC
eT='2025-05-31T00:00:00'	# UTC
module_path='~/Documents/Codes/Modules/Solar_data_mining-main/'
output_folder='./results/STIX_GOES_20250530/'

# ---- Sampling ----
stix_rebin=60			# s; None = native 4 s
goes_res_option=1		# 1: avg1m, 2: flx1s

# ---- Instrument settings ----
earth_time=True
goes_satellite=None		# None = highest available
goes_channel_option=1	# 1: 1-8 A (xrsb), 2: 0.5-4 A (xrsa)
stix_bands=None			# None = all bands

# ---- Output ----
stix_ylim=(5,400)
goes_ylim=(1e-7,1e-3)
fig_format='png'
out_tag=''
###############################################################


def main():
	import os,sys
	sys.path.insert(0,os.path.expanduser(module_path))
	import Multiwb_datamine as mwd
	import xray_lightcurves as xlc

	goes_res={1:'avg1m',2:'flx1s'}[goes_res_option]
	goes_ch={1:'xrsb',2:'xrsa'}[goes_channel_option]
	out=os.path.abspath(os.path.expanduser(output_folder))
	tag=out_tag or (sT.replace('-','').replace(':','')+'-'+eT.replace('-','').replace(':',''))

	stix_df,stix_csv=mwd.get_STIX_LC(sT,eT,outdir=os.path.join(out,'STIX'),earth_time=earth_time,plot=False)
	goes_df,goes_csv,sat=xlc.get_GOES_XRS(sT,eT,outdir=os.path.join(out,'GOES'),satellite=goes_satellite,resolution=goes_res)

	tlabel='Time (UT, Earth light-arrival)' if earth_time else 'Time (UT; STIX at Solar Orbiter, GOES at Earth)'
	fig=xlc.plot_STIX_GOES(stix_df,goes_df,os.path.join(out,f'STIX_GOES{sat}_{tag}.{fig_format}'),
		bands=stix_bands,stix_rebin=stix_rebin,goes_channel=goes_ch,goes_label=f'GOES-{sat}',
		stix_ylim=stix_ylim,goes_ylim=goes_ylim,xlim=(sT,eT),
		time_label=tlabel+' on '+sT[:10],
		title=f'Solar Orbiter/STIX quick-look and GOES-{sat} XRS')
	return {'stix_csv':stix_csv,'goes_csv':goes_csv,'figure':fig}


if __name__=='__main__':
	main()
