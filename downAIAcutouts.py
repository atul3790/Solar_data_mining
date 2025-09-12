'''
Code downloads calibrates and makes AIA cutout movies of AR patches of given size. Pass a table with Date, stime and ptime for Date, start time and peak time of the flares of interest.
LC_fold: folder with AR avgd light curve
'''

import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np,pickle
from sunpy.net import Fido, attrs as a
#from sunpy.map.sources import sdo
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
from datetime import datetime, timedelta
from tqdm import tqdm
import matplotlib.animation as animation
from aiapy.calibrate import register, update_pointing, correct_degradation, degradation, normalize_exposure
#from aiapy.calibrate.util import get_pointing_table, get_correction_table
from sunpy.map import Map
from aiapy.psf import deconvolve, psf
from sunpy.time import parse_time
from sunpy.coordinates import frames
import pandas as pd
from multiprocess import Pool

####### INPUT #########################
basedir='/gpfsm/dnb34/amohan/DHtypeIII/'
evTab='gC5_TypeIIIFlx_in55lon_200612-202501_refined.csv'
#wavelengths=[1600,304,171,94,193]
wavelengths=[304,94,193]
plotdir='AR_sourceprops_AIA304-94-193/' # Adress witihn basedir where to save plots. LC_fold and mpr_fold will be made by the code within this
LC_fold='AR_LC/' # Avgd AR flux light curve in all bands of choice
mp4_fold='AR_movies/'
date_col='Date' # Column heading for date
time_col='stime' # Column heading for time
time_pk='ptime' # Flare peak time.
Class='Class' # Leave as '' if not needed
src_loc='helio' # Leave as '' if not needed
time_radio_pk='Time_wrt_pflr_min' # If you have a column in table which says around which time with respect to flare peak time in minutes to have best cadence in dynamic cadence fixing. If you want the high cadence AIA data to be around SXR peak time leave this as 0. Else specify the column in the Table from where SXR peak time - chosen time in minutes is listed. Code expect a float which can be -ve or positive if your chosen time is earlier or later to the peak SXR time.
size_arcsec=400 # Size of AR patch in arcsec
#############################################

cwd=os.getcwd()
os.chdir(basedir)
os.system('mkdir '+plotdir)
os.system('mkdir '+plotdir+'temp/')
os.system('mkdir '+plotdir+LC_fold)
os.system('mkdir '+plotdir+LC_fold+'pickle/')
os.system('mkdir '+plotdir+mp4_fold)
LC_fold=plotdir+LC_fold
mp4_fold=plotdir+mp4_fold

Tab=pd.read_csv(evTab)
nproc=os.cpu_count()
####### Functions ##########################################
def download_aia_cutouts(start_time, end_time, peak_time,event_time,lat,lon, size_arcsec, wavelengths,pk_time,mp4_fold):
    """
    Downloads AIA cutouts between start_time and end_time in given wavelengths.
    Returns dict of {wavelength: list of SunPy maps}.
    """
    #print('hi')
    #return 0
    Rad_pktime=peak_time+timedelta(seconds=pk_time)
    os.system('mkdir '+mp4_fold+event_time.strftime("%Y%m%dT%H%M%S")+'/')
    all_maps = {wl: [] for wl in wavelengths}
    LCs={wl:{'Time':[], 'AR flux':[], 'Sun-as-a-star flux':[]} for wl in wavelengths}
    #fluxes = {wl: [] for wl in wavelengths}
    #sunasstar= {wl: [] for wl in wavelengths}
    #times = {wl: [] for wl in wavelengths}
    t_begin = parse_time("2010-03-25T00:00:00")
    now = Time.now()
    time_window = t_begin + np.arange(0, (now - t_begin).to(u.day).value, 7) * u.day
    hgs_coord = SkyCoord(lon, lat,
                     frame=frames.HeliographicStonyhurst,
                     observer='earth',
                     obstime=start_time)
    center_coord = hgs_coord.transform_to(frames.Helioprojective(observer='earth', obstime=start_time))   
    #pointing_table = get_pointing_table(Time(start_time), Time(end_time))
    #correction_table = get_correction_table("JSOC")

    # Query cadence (adjust if needed)
    #cadence = timedelta(seconds=60) # 1 minutes roughly, to limit download volume
    cadence = timedelta(seconds=1800) # 20 min cadence for pre flar start time	    
    time_range = start_time
    #flag=0
    evrise=(peak_time-event_time).seconds
    #num=0
    #temp=[]
    while time_range < end_time+cadence:        
        if (event_time-time_range).seconds<3600:
            cadence=timedelta(seconds=1200) # 20 minute cadence        
        if (event_time-time_range).seconds<1200:
            cadence=timedelta(seconds=180)  # 3 minute cadence
        if time_range>event_time and time_range<peak_time+timedelta(seconds=60):
            cadence=timedelta(seconds=120) # 2 min cadence
        if time_range>peak_time+timedelta(seconds=60):
            cadence=timedelta(seconds=180) # 3 min cadence            
        if time_range>peak_time+timedelta(seconds=evrise):
            cadence=timedelta(seconds=600) # 10 min cadence
        if time_range>peak_time+timedelta(seconds=evrise*1.75):
            cadence=timedelta(seconds=900) # 15 min cadence            
        if time_range>peak_time+timedelta(seconds=4500):			
            cadence=timedelta(seconds=1800) # 30 min cadence
        if time_range>Rad_pktime - timedelta(seconds=int(evrise/2)) and time_range<Rad_pktime + timedelta(seconds=int(evrise/2)):
            cadence=timedelta(seconds=60) # 1 minute cadence
        #num+=1
        print('Time:' , time_range, 'Cadence: ',cadence.seconds)		
        next_time = time_range + cadence
        #temp+=[(time_range,cadence)]
        #time_range=next_time
		
        # Define in Heliographic Stonyhurst (HGS) coordinates

        # Convert to Helioprojective (HPC)
        
        if next_time > end_time+cadence:
            next_time = end_time+cadence
        for wl in wavelengths:
            query = Fido.search(a.Time(time_range, next_time), a.Instrument('AIA'),a.Wavelength(wl * u.angstrom),a.Sample(12 * u.second))
            if len(query) == 0:
                print('Failed query!!')
                continue

            for inx in range(len(query[0])):
                try:
                    files = Fido.fetch(query[0,inx],path=plotdir+'temp/')  # get first match in interval
                    if os.path.isfile(files[0])!=True:
			print('No file downloaded for: ',files[0].split('/')[-1])
                        continue
                    m = Map(files[0])
                    m=update_pointing(m) # Doing rotation corrections
                    m=register(m) # Pixels fluxes are now corrected for direction offset issue and pixels rescaled to 0.6" and image y axes are now aligned with solar north-south axis
                    m= correct_degradation(m) # Time dependent instrument degradation correction
                    m=normalize_exposure(m) # Making map normalised with respect to exposure times to get values in DN/s/pixel
                    # Create submap cutout
                    bottom_left = SkyCoord(center_coord.Tx - size_arcsec/2 * u.arcsec,
                                           center_coord.Ty - size_arcsec/2 * u.arcsec,
                                           frame=m.coordinate_frame)

                    yy, xx = np.indices(m.data.shape)
                    cx, cy = m.reference_pixel.x.value, m.reference_pixel.y.value
                    r_pix = (m.rsun_obs.to('arcsec').value/m.meta['cdelt1'])  # pixels
                    mask = (xx - cx)**2 + (yy - cy)**2 <= r_pix**2

                    arcsec_to_rad = (np.pi/180.0)/3600.0
                    #sun_sr=(16*3600)*arcsec_to_rad**2
                    #sub_sr=size_arcsec*arcsec_to_rad**2
                    #pix_sr = (m.meta['cdelt1']*arcsec_to_rad) * (m.meta['cdelt2']*arcsec_to_rad)

                    fov={"width": size_arcsec * u.arcsec, "height": size_arcsec * u.arcsec}
                    submap = m.submap(bottom_left, **fov)

                    print("Computing integrated AR flux and sun-as-a-star disk averaged flux...")
            
                    # Assuming all wavelengths have maps at roughly the same time,
                    # we will use times from the first wavelength that has data.
                    LCs[wl]['Time'].append(m.date.datetime)
                    LCs[wl]['AR flux'].append(submap.data.sum())
                    LCs[wl]['Sun-as-a-star flux'].append(m.data[mask].sum())
            
                    if time_range==start_time or (peak_time-time_range).seconds<360 or (end_time-time_range).seconds<2000:
                        fig = plt.figure()
                        ax = plt.subplot(projection=submap)
                        submap.plot(axes=ax)
                        fig.savefig(mp4_fold + event_time.strftime("%Y%m%dT%H%M%S") + '/'+time_range.strftime("%Y%m%dT%H%M%S") + '.png')
                        plt.close()
                        all_maps[wl].append(submap)

                    else:
                        os.system('rm -rf '+files[0])
                    #submap.peek()
                    #plt.savefig(mp4_fold+event_time.strftime("%Y%m%dT%H%M%S")+'/'+time_range.strftime("%Y%m%dT%H%M%S")+'.png')
                    #plt.close()
                    #all_maps[wl].append(submap)
                    del m
                    del submap
                    print('Found a file: ',files[0].split('/')[-1])
                except:
                    print(f"Failed fetching {wl} at {time_range}")
                    
        time_range = next_time
        
    pickle.dump(LCs,open(LC_fold+'pickle/'+event_time.strftime("%Y%m%dT%H%M%S")+'_ARflare'+'.p','wb'))
    return all_maps, LCs

def make_mp4_from_maps(maps, clas, src, output_filename, fps=5):
    """
    Create MP4 video from a list of SunPy maps.
    """
    fig = plt.figure(figsize=(6, 6))
    plt.axis('off')
    ims = []

    for m in maps:
        im = plt.imshow(m.data, cmap=m.cmap, origin='lower', animated=True)
        plt.title(clas+'; '+src, size=14)
        ims.append([im])

    ani = animation.ArtistAnimation(fig, ims, interval=1000/fps, blit=True)
    ani.save(output_filename)
    plt.close(fig)

def plot_lightcurves(LCs, start_time,peak_time,end_time, event_time,wavelengths,clas,src,filen):
    """
    Plot integrated flux lightcurves for all wavelengths.
    """
    for tag in ['AR flux','Sun-as-a-star flux']:
        tg=tag.split(' ')[0].replace('-','')
        plt.ioff()
        plt.figure(figsize=(8, 5))
        fig,ax=plt.subplots(figsize=(8, 5))
        flg=1
        for wl in wavelengths:
            try:
                ax.plot(LCs[wl]['Time'], LCs[wl][tag]/np.nanmax(LCs[wl][tag]), label=f'{wl} Å')
            except:
                flg+=1
                pass
        if flg==len(wavelengths):
            plt.savefig(filen+event_time.strftime("%Y%m%dT%H%M%S")+'_'+tg+'_LC.png')
            continue
        plt.xlim([start_time,end_time])
        plt.axvline(peak_time,color='k',linestyle='-')
        plt.axvline(start_time,color='b',linestyle='--')
        plt.axvline(peak_time+timedelta(hours=1),color='b',linestyle='--')
        plt.tight_layout()
        plt.xticks(size=14)
        plt.yticks(size=14)
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        plt.xlabel('Time (UT) +'+event_time.strftime("%Y-%m-%d"), size=14)
        plt.ylabel('Normalised Flux ('+r'$\mathrm{DN\, s^{-1})}$',size=14)
        plt.title(tag+' light curve ('+clas+'; '+src+')',size=14)
        plt.legend(loc='best',prop={'size':14})
        plt.grid(True)        
        plt.tight_layout()
        plt.savefig(filen+event_time.strftime("%Y%m%dT%H%M%S")+'_'+tg+'_LC.png')
        plt.close()
    return
###############################################################################

Tab=Tab[12:].reset_index(drop=True)
p=Pool(10)
#for i in range(len(Tab)):
def Run_fn (i):
    event_time=np.datetime64(Tab[date_col][i]+'T'+Tab[time_col][i]).astype(datetime)
    peak_time=np.datetime64(Tab[date_col][i]+'T'+Tab[time_pk][i]).astype(datetime)    
    end_time=event_time+timedelta(hours=3)
    start_time=event_time-timedelta(hours=2)
    clas=Tab[Class][i]
    src=Tab[src_loc][i]
    sig={'N':1,'S':-1,'E':-1,'W':1}
    lat=int(src[1:3])*sig[src[0].upper()]*u.deg
    lon=int(src[4:])*sig[src[3].upper()]*u.deg
    try:
        pk_time=Tab[time_radio_pk][i]*60 # Converting minutes to seconds
    except:
        pk_time=0
    print("Downloading AIA cutouts...")
    all_maps,LCs = download_aia_cutouts(start_time, end_time, peak_time,event_time,lat,lon, size_arcsec, wavelengths,pk_time,mp4_fold)
    '''	
    print("Processing and saving MP4 videos...")
    for wl in wavelengths:
        if len(all_maps[wl]) == 0:
            print(f"No data for wavelength {wl} Å")
            continue
        output_mp4 = mp4_fold+event_time.strftime("%Y%m%dT%H%M%S")+f'_AIA{wl}A_cutout.mp4'
        make_mp4_from_maps(all_maps[wl], clas, src, output_mp4, fps=2)
        print(f"Saved video: {output_mp4}")
    '''
    #print("Computing integrated flux and plotting lightcurves...")
    #plot_lightcurves(LCs,start_time,peak_time,end_time, event_time,wavelengths,clas,src,LC_fold)

succ=p.map(Run_fn,np.arange(len(Tab[1:2])))

os.chdir(cwd)
