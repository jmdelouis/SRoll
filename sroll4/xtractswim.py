import numpy as np
import healpy as hp
from os import listdir
import os 
from datetime import datetime
import sys
import xarray as xr
import os
import glob
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE' # command line : export HDF5_USE_FILE_LOCKING=FALSE

# ==========================================================================================
#
#                       Convert l2s Swim data to healpix format for SRoll
#
# ==========================================================================================

# Input data, l2s swim cfosat 
path = '/home/ref-cfosat-public/datasets/swi_l2s/'
# path = '/espace/dataref/ref-cfosat-public/datasets/swi_l2s/'
order=['v1.0']

# Reference time : second since 1970
t0=datetime(1970,1,1,0,0,0)

# Ring size
RGSIZE=100000

# Year 
year=[2023]

# Number of days to compute
n_jour=366

# first day 
start_jour = 1

# Convert to db 
convert_to_db = 1

# Healpix size 
nside = 512

# Detector 
angle = 'L2S04'

# Output name
outname='SWIM'+str(int(year[0]))+'_V0_filt'

# Path output 
outpath='/home/datawork-cersat-public/cache/project/deepsee2/sroll_input/'
# outpath='/espace/datawork/datawork-cersat-public/cache/project/deepsee2/sroll_input/'

# ??
time_scale = 1

#============================================================================================
#============================================================================================
# DATA MODEL
#============================================================================================
# 
#
#============================================================================================

# Resolution
print(
    "Approximate resolution at NSIDE {} is {:.2} deg".format(
        nside, hp.nside2resol(nside, arcmin=True) / 60
    )
)

# Number of pixel NPIX of the map
NPIX = hp.nside2npix(nside)
print(
    "Number of pixels NPIX of the map (nside = {}) : {}".format(
        nside,NPIX))

__file__list__={}

def append_file(path,val,dtype,step):
    data=np.zeros([step,RGSIZE],dtype=dtype)
    nn=val.shape[0]//step
    for i in range(step):
        data[i,0:nn]=val[nn*i:nn*(i+1)].astype(dtype)
    if nn*(step)<val.shape[0]:
        data[step-1,nn:nn+val.shape[0]-nn*(step)]=val[nn*(step):].astype(dtype)

    __file__list__[path].write(data.flatten())

def dB2linear(values_ind_db):
    """convert sigma0 from dB to linear"""         
    return np.ma.power(10.,values_ind_db / 10.).data

def linear2dB(values): 
    """convert sigma0 from linear to dB"""
    return 10. * np.ma.log10(values).data

def init_file(path):
    __file__list__[path]=open(path, "wb")

def close_file(path):
    __file__list__[path].close()


init_file(outpath+outname+'_'+angle+'_hit')
init_file(outpath+outname+'_'+angle+'_sig')
init_file(outpath+outname+'_'+angle+'_ptg')
init_file(outpath+outname+'_'+angle+'_ptg_TUPLE_1')
init_file(outpath+outname+'_'+angle+'_ptg_TUPLE_2')
init_file(outpath+outname+'_'+angle+'_ant_azi')
init_file(outpath+outname+'_'+angle+'_azi')
init_file(outpath+outname+'_'+angle+'_time')



# =========================================================================================================
# THE PURPOSE OF THIS CODE IS TO RETRIEVE SWIM DATA AND AVERAGE THE MEASURED VALUE
# OF SCAT IN HEALPIX PIXELS OBSERVED AT THE SAME INCIDENCE AND ORIENTATION DURING AN ORBIT.
# ONCE THIS OPERATION IS COMPLETE, SROLL CAN WORK IN HEALPIX.
# ======================================================================================================
# hit: Number of times a measurement is taken in the pixel
# sig: The data itself
# ptg, ptg_TUPLE_1, ptg_TUPLE_2: Pointing data, respectively colatitude, longitude, and incidence
# ant_azi: Antenna rotation angle
# azi: Measurement azimuth
# time: Date

for y in year: 
    for d in range(n_jour): 

        print('JOUR '+str(d)+'/'+str(n_jour))

        i_ord_fin = 0
        n_ord_max = 0
        i_ord_max = 0

        for i_ord in order :
            if os.path.exists(path+i_ord+'/%d/%03d/'%(y,d+start_jour)):
                filetab=listdir(path+i_ord+'/%d/%03d/'%(y,d+start_jour))
                if len(filetab)>0:
                    if (len(filetab) > 12) :
                        i_ord_fin = i_ord
                        break 
                    elif len(filetab) > n_ord_max :
                        n_ord_max = len(filetab)
                        i_ord_fin = i_ord
        
        print(i_ord_fin)          
        if i_ord_fin==0:
            continue  

        filetab = glob.glob(path+i_ord_fin+'/%d/%03d/'%(y,d+start_jour)+'*/*'+angle+'*.nc',recursive=True)
        filetab=np.sort(filetab)

        num_file = len(filetab)

        for ifile in range(num_file):
            try:

                name_file = os.path.basename(filetab[ifile][:])     #.split('.')[0]
                yy=int(name_file[22:26])
                mm=int(name_file[26:28])
                dd=int(name_file[28:30])
                hh=int(name_file[31:33])
                mi=int(name_file[33:35])
                ss=int(name_file[35:37])
                tt=datetime(yy,mm,dd,hh,mi,ss)

                duration = (tt - t0).total_seconds()                

                # read data
                ds = xr.open_dataset(filetab[ifile][:])
                sigma0 = ds['seg_sigma0_mean'][:,:].data
                inc    = ds['seg_incidence'][:,:].data.flatten()
                lat    = ds['seg_lat'][:,:].data.flatten()
                lon    = ds['seg_lon'][:,:].data.flatten()
                seg_flag    = ds['seg_flag'][:,:].data.flatten()

                # nt = frame, nx = frame_pulse, ny = pulse_slice
                nt,nx = sigma0.shape       
                sigma0 = sigma0.flatten()
                azi = np.repeat(ds['phi_geo'][:].data,nx)

                # data frame relative time of present orbit, dim (frame)
                time = ds['time']-np.datetime64(tt)
                time = time.dt.seconds+time.dt.microseconds*0.000001
                time = time.data
                time = np.repeat(time,nx)

                # Dim of antenna_azimuth : (frame, frame_pulse)
                ant = np.repeat(ds['phi'][:].data,nx)
                ds.close()

                # COnvert lat lon
                lat=(90.0-lat)/180.*np.pi
                lon=2*np.pi-np.fmod(lon/180.*np.pi+2*np.pi,2*np.pi)

                # Detect bad data
                flag_nan = np.isnan(sigma0) | np.isnan(ant)
                flag_sigma = np.bitwise_and(seg_flag.astype(int),2)==2
                flag = np.where( (~flag_nan) & (~flag_sigma) & (sigma0>=0))[0]

                if len(flag)==0:
                    print('No valid data for this file')

                hidx=hp.ang2pix(nside,lat[flag],lon[flag])
                lat=lat[flag]
                lon=lon[flag]
                inc=inc[flag]/180.*np.pi
                ant=ant[flag]/180.*np.pi
                sigma0 = sigma0[flag]

                iazi=((azi[flag]/10).astype('int'))
                iazi[iazi==36]=0
                time_pola=time[flag]+duration #second since 1970

                hidx=hidx+iazi*12*nside**2
                hit = np.bincount(hidx,minlength=36*12*nside**2)
                idx = np.where(hit>0)[0]
                tim = np.bincount(hidx,weights=time_pola,minlength=36*12*nside**2)
                sig = np.bincount(hidx,weights=sigma0,minlength=36*12*nside**2)
                clat =np.bincount(hidx,weights=np.cos(lat),minlength=36*12*nside**2)
                slat = np.bincount(hidx,weights=np.sin(lat),minlength=36*12*nside**2)
                clon = np.bincount(hidx,weights=np.cos(lon),minlength=36*12*nside**2)
                slon = np.bincount(hidx,weights=np.sin(lon),minlength=36*12*nside**2)
                cinc = np.bincount(hidx,weights=np.cos(inc),minlength=36*12*nside**2)
                sinc = np.bincount(hidx,weights=np.sin(inc),minlength=36*12*nside**2)
                cant = np.bincount(hidx,weights=np.cos(ant),minlength=36*12*nside**2)
                sant = np.bincount(hidx,weights=np.sin(ant),minlength=36*12*nside**2)
                cazi = np.bincount(hidx,weights=np.cos(azi[flag]/180.*np.pi),minlength=36*12*nside**2)
                sazi = np.bincount(hidx,weights=np.sin(azi[flag]/180.*np.pi),minlength=36*12*nside**2)
                sig=(sig[idx]/hit[idx])

                # Convert to db 
                if convert_to_db==1:
                    sig=linear2dB(sig)

                lat=np.arctan2(slat[idx],clat[idx])
                lon=np.arctan2(slon[idx],clon[idx])
                azi=(np.arctan2(sazi[idx],cazi[idx]))
                ant=(np.arctan2(sant[idx],cant[idx]))
                inc=np.arctan2(sinc[idx],cinc[idx])

                tim=tim[idx]/hit[idx]  
                hit=hit[idx]
                append_file(outpath+outname+'_'+angle+'_hit',hit,'float32',time_scale)
                append_file(outpath+outname+'_'+angle+'_sig',sig,'float32',time_scale)
                append_file(outpath+outname+'_'+angle+'_ptg',lon,'float64',time_scale)
                append_file(outpath+outname+'_'+angle+'_ptg_TUPLE_1',lat,'float64',time_scale)
                append_file(outpath+outname+'_'+angle+'_ptg_TUPLE_2',inc,'float64',time_scale)
                append_file(outpath+outname+'_'+angle+'_ant_azi',ant,'float32',time_scale)
                append_file(outpath+outname+'_'+angle+'_azi',azi,'float32',time_scale)
                append_file(outpath+outname+'_'+angle+'_time',tim,'float32',time_scale)

                print(filetab[ifile],np.sum(hit>0),RGSIZE*time_scale,tt.ctime(),duration,time.min(),time.max())
                sys.stdout.flush()

            except Exception as err:
                print('PROBLEM WHILE READING ',filetab[ifile])
                print(f"Unexpected {err=}, {type(err)=}")
                append_file(outpath+outname+'_'+angle+'_hit',np.array([]),'float32',time_scale)
                append_file(outpath+outname+'_'+angle+'_sig',np.array([]),'float32',time_scale)
                append_file(outpath+outname+'_'+angle+'_ptg',np.array([]),'float64',time_scale)
                append_file(outpath+outname+'_'+angle+'_ptg_TUPLE_1',np.array([]),'float64',time_scale)
                append_file(outpath+outname+'_'+angle+'_ptg_TUPLE_2',np.array([]),'float64',time_scale)
                append_file(outpath+outname+'_'+angle+'_ant_azi',np.array([]),'float32',time_scale)
                append_file(outpath+outname+'_'+angle+'_azi',np.array([]),'float32',time_scale)
                append_file(outpath+outname+'_'+angle+'_time',np.array([]),'float32',time_scale)
                sys.stdout.flush()
                

close_file(outpath+outname+'_'+angle+'_hit')
close_file(outpath+outname+'_'+angle+'_sig')
close_file(outpath+outname+'_'+angle+'_ptg')
close_file(outpath+outname+'_'+angle+'_ptg_TUPLE_1')
close_file(outpath+outname+'_'+angle+'_ptg_TUPLE_2')
close_file(outpath+outname+'_'+angle+'_ant_azi')
close_file(outpath+outname+'_'+angle+'_azi')
close_file(outpath+outname+'_'+angle+'_time')

print('DONE')

