import numpy as np
import healpy as hp
from os import listdir
import os 
from datetime import datetime
import sys
import xarray as xr
import os

# ==========================================================================================
#
#                       Convert l1b Scat data to healpix format for SRoll
#
# ==========================================================================================


# =============================================================================================
#                                       Define parameters
# =============================================================================================

# Input data, l1b cfosat
path = '/home/datawork-cersat-public/provider/cnes/satellite/l1b/cfosat/scat/sca_l1b___/'

# Available versions ranked in their order of preference
order=['v3.3.1','3.3','3.2','3.0','10.10']

# Reference time : second since 1970
t0=datetime(1970,1,1,0,0,0)

# Ring size (3 000 000 is enough for scat)
RGSIZE=3000000 # nside 512

# Year 
year=[2021]

# Number of days to compute
n_jour=366

# first day (1 for first day of the year)
start_day = 3

# Convert to db 
convert_to_db = 1

# Healpix size 
nside = 512

# Output name
outname='SCAT2021_example_' 

# Output path  
outpath='/home/datawork-cersat-public/cache/project/deepsee2/sroll_input/'

# flag values for each polarisation
value_hh = [2,18]
value_vv = [0,16]

# time scale
time_scale = 1

#============================================================================================
# Define functions
#============================================================================================

# Resolution
print("Approximate resolution at NSIDE {} is {:.2} deg".format(
        nside, hp.nside2resol(nside, arcmin=True) / 60))

# Number of pixel NPIX of the map
NPIX = hp.nside2npix(nside)
print("Number of pixels NPIX of the map (nside = {}) : {}".format(nside,NPIX))


__file__list__={}

def append_file(path,val,dtype,step):
    """ Write data """
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
    """ Init file """
    __file__list__[path]=open(path, "wb")

def close_file(path):
    """ Close file """
    __file__list__[path].close()

for pola in ['hh','vv']:
    init_file(outpath+outname+'_'+pola+'_hit')
    init_file(outpath+outname+'_'+pola+'_sig')
    init_file(outpath+outname+'_'+pola+'_ptg')
    init_file(outpath+outname+'_'+pola+'_ptg_TUPLE_1')
    init_file(outpath+outname+'_'+pola+'_ptg_TUPLE_2')
    init_file(outpath+outname+'_'+pola+'_ant_azi')
    init_file(outpath+outname+'_'+pola+'_azi')
    init_file(outpath+outname+'_'+pola+'_time')
    init_file(outpath+outname+'_'+pola+'_kp')


# =========================================================================================================
# THE PURPOSE OF THIS CODE IS TO RETRIEVE SCAT DATA AND AVERAGE THE MEASURED VALUE
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
        i_ord_fin = 0
        n_ord_max = 0
        i_ord_max = 0

        # To select the good data version
        for i_ord in order :
            if os.path.exists(path+i_ord+'/%d/%03d/'%(y,d+start_day)):
                filetab=listdir(path+i_ord+'/%d/%03d/'%(y,d+start_day))
                if len(filetab)>0:
                    if (len(filetab) > 12) :
                        i_ord_fin = i_ord
                        break 
                    elif len(filetab) > n_ord_max :
                        n_ord_max = len(filetab)
                        i_ord_fin = i_ord
        if i_ord_fin==0:
            continue

        print(i_ord_fin)        
        filetab = listdir(path+i_ord_fin+'/%d/%03d/'%(y,d+start_day))
        filetab=np.sort(filetab)

        print('JOUR '+str(d)+'/'+str(n_jour))
        num_file = len(filetab)
        num_rand = np.random.randint(0,num_file)

        # Loop on files found for day d
        for ifile in range(num_file):
            try:
                yy=int(filetab[ifile][22:26])
                mm=int(filetab[ifile][26:28])
                dd=int(filetab[ifile][28:30])
                hh=int(filetab[ifile][31:33])
                mi=int(filetab[ifile][33:35])
                ss=int(filetab[ifile][35:37])
                tt=datetime(yy,mm,dd,hh,mi,ss)
                duration = (tt - t0).total_seconds()
                name_file = os.path.basename(path+i_ord_fin+'/%d/%03d/%s'%(y,d+start_day,filetab[ifile])).split('.')[0]

                # Read data
                ds = xr.open_dataset(path+i_ord_fin+'/%d/%03d/%s'%(y,d+start_day,filetab[ifile]))
                sigma0 = ds['slice_sigma0'][:,:,:].data
                inc    = ds['slice_incidence'][:,:,:].data.flatten()
                ele    = ds['slice_elevation'][:,:,:].data.flatten()
                lat    = ds['slice_lat'][:,:,:].data.flatten()
                lon    = ds['slice_lon'][:,:,:].data.flatten()
                azi    = ds['slice_azimuth'][:,:,:].data.flatten()
                kpa = ds['slice_kpc_a'][:,:,:].data.flatten()
                kpb = ds['slice_kpc_b'][:,:,:].data.flatten()
                kpc = ds['slice_kpc_c'][:,:,:].data.flatten()
                snr = ds['slice_snr'][:,:,:].data.flatten()


                # nt = frame, nx = frame_pulse, ny = pulse_slice
                nt,nx,ny=sigma0.shape       
                sigma0=sigma0.flatten()

                # noise ratio
                kp = linear2dB(np.sqrt(kpa*sigma0**2+kpb*sigma0+kpc))

                # data frame relative time of present orbit, dim (frame)
                time = ds['orbit_time'][:].data.flatten()
                time = np.repeat(time,nx*ny)

                # Flag used to select polarisation 
                cell_mode_flag = np.repeat(ds['cell_mode_flag'][:,:].data.flatten(),ny)

                # Dim of antenna_azimuth : (frame, frame_pulse)
                ant = np.repeat(ds['antenna_azimuth'][:,:].data.flatten(),ny)
                ds.close()

                # Convert lat to colatitude and lon
                lat=(90.0-lat)/180.*np.pi
                lon=2*np.pi-np.fmod(lon/180.*np.pi+2*np.pi,2*np.pi)

                # flag for bad data
                flag_inc = (inc<48) & (inc>28)
                flag_bad_snr = (snr<=0) & (linear2dB(sigma0)>-20)
                flag_bad_kp = kp<-9
                flag_nan = np.isnan(sigma0)


                # flag hh / flag vv
                for pola in ['hh','vv']:
                    value_pola = eval('value_'+pola)

                    # Select val ok 
                    flag = np.where(np.isin(cell_mode_flag,value_pola) & flag_inc \
                                    & ~flag_bad_snr & ~flag_bad_kp & ~flag_nan)[0]

                    if len(flag)==0:
                        print('No valid data for this file')

                    hidx=hp.ang2pix(nside,lat[flag],lon[flag])
                    lat_pola=lat[flag]
                    lon_pola=lon[flag]
                    inc_pola=inc[flag]/180.*np.pi
                    ant_pola=ant[flag]/180.*np.pi
                    sigma0_pola = sigma0[flag]
                    iazi_pola=((azi[flag]/10).astype('int'))
                    iazi_pola[iazi_pola==36]=0
                    time_pola=time[flag]+duration # second since 1970
                    kp_pola = kp[flag]
                    azi_pola=azi[flag]

                    # detect time outliers
                    val_out = (np.abs(time_pola - np.median(time_pola)) > 3 * 60 * 60)
                    if np.sum(val_out) > 0:
                        print('Warning for file ' + name_file + ' : ' + str(np.sum(val_out)) + ' outliers detected')
                        hidx=hidx[~val_out]
                        lat_pola=lat_pola[~val_out]
                        lon_pola=lon_pola[~val_out]
                        inc_pola=inc_pola[~val_out]
                        ant_pola=ant_pola[~val_out]
                        sigma0_pola = sigma0_pola[~val_out]
                        iazi_pola=iazi_pola[~val_out]
                        time_pola=time_pola[~val_out] # second since 1970
                        kp_pola=kp_pola[~val_out]
                        azi_pola = azi_pola[~val_out]


                    hidx=hidx+iazi_pola*12*nside**2
                    hit = np.bincount(hidx,minlength=36*12*nside**2)
                    idx = np.where(hit>0)[0]
                    tim = np.bincount(hidx,weights=time_pola,minlength=36*12*nside**2)
                    sig = np.bincount(hidx,weights=sigma0_pola,minlength=36*12*nside**2)
                    kp_bin = np.bincount(hidx,weights=kp_pola,minlength=36*12*nside**2)
                    clat =np.bincount(hidx,weights=np.cos(lat_pola),minlength=36*12*nside**2)
                    slat = np.bincount(hidx,weights=np.sin(lat_pola),minlength=36*12*nside**2)
                    clon = np.bincount(hidx,weights=np.cos(lon_pola),minlength=36*12*nside**2)
                    slon = np.bincount(hidx,weights=np.sin(lon_pola),minlength=36*12*nside**2)
                    cinc = np.bincount(hidx,weights=np.cos(inc_pola),minlength=36*12*nside**2)
                    sinc = np.bincount(hidx,weights=np.sin(inc_pola),minlength=36*12*nside**2)
                    cant = np.bincount(hidx,weights=np.cos(ant_pola),minlength=36*12*nside**2)
                    sant = np.bincount(hidx,weights=np.sin(ant_pola),minlength=36*12*nside**2)
                    cazi = np.bincount(hidx,weights=np.cos(azi_pola/180.*np.pi),minlength=36*12*nside**2)
                    sazi = np.bincount(hidx,weights=np.sin(azi_pola/180.*np.pi),minlength=36*12*nside**2)

                    sig=(sig[idx]/hit[idx])
                    kp_bin=(kp_bin[idx]/hit[idx])

                    # Convert to db 
                    if convert_to_db==1:
                        sig=linear2dB(sig)

                    lat_pola=np.arctan2(slat[idx],clat[idx])
                    lon_pola=np.arctan2(slon[idx],clon[idx])
                    azi_pola=(np.arctan2(sazi[idx],cazi[idx]))
                    ant_pola=(np.arctan2(sant[idx],cant[idx]))
                    inc_pola=np.arctan2(sinc[idx],cinc[idx])

                    tim=tim[idx]/hit[idx]  
                    hit=hit[idx]

                    # Write data
                    append_file(outpath+outname+'_'+pola+'_hit',hit,'float32',time_scale)
                    append_file(outpath+outname+'_'+pola+'_sig',sig,'float32',time_scale)
                    append_file(outpath+outname+'_'+pola+'_ptg',lon_pola,'float64',time_scale)
                    append_file(outpath+outname+'_'+pola+'_ptg_TUPLE_1',lat_pola,'float64',time_scale)
                    append_file(outpath+outname+'_'+pola+'_ptg_TUPLE_2',inc_pola,'float64',time_scale)
                    append_file(outpath+outname+'_'+pola+'_ant_azi',ant_pola,'float32',time_scale)
                    append_file(outpath+outname+'_'+pola+'_azi',azi_pola,'float32',time_scale)
                    append_file(outpath+outname+'_'+pola+'_time',tim,'float32',time_scale)
                    append_file(outpath+outname+'_'+pola+'_kp',kp_bin,'float32',time_scale)

                print(filetab[ifile],np.sum(hit>0),RGSIZE*time_scale,tt.ctime(),duration,time.min(),time.max())
                sys.stdout.flush()
            except Exception as err:
                print('PROBLEM WHILE READING ',filetab[ifile])
                print(f"Unexpected {err=}, {type(err)=}")

                for pola in ['hh','vv']:
                    append_file(outpath+outname+'_'+pola+'_hit',np.array([]),'float32',time_scale)
                    append_file(outpath+outname+'_'+pola+'_sig',np.array([]),'float32',time_scale)
                    append_file(outpath+outname+'_'+pola+'_ptg',np.array([]),'float64',time_scale)
                    append_file(outpath+outname+'_'+pola+'_ptg_TUPLE_1',np.array([]),'float64',time_scale)
                    append_file(outpath+outname+'_'+pola+'_ptg_TUPLE_2',np.array([]),'float64',time_scale)
                    append_file(outpath+outname+'_'+pola+'_ant_azi',np.array([]),'float32',time_scale)
                    append_file(outpath+outname+'_'+pola+'_azi',np.array([]),'float32',time_scale)
                    append_file(outpath+outname+'_'+pola+'_time',np.array([]),'float32',time_scale)
                    append_file(outpath+outname+'_'+pola+'_kp',np.array([]),'float32',time_scale)  
                sys.stdout.flush()
                
for pola in ['hh','vv']:
    close_file(outpath+outname+'_'+pola+'_hit')
    close_file(outpath+outname+'_'+pola+'_sig')
    close_file(outpath+outname+'_'+pola+'_ptg')
    close_file(outpath+outname+'_'+pola+'_ptg_TUPLE_1')
    close_file(outpath+outname+'_'+pola+'_ptg_TUPLE_2')
    close_file(outpath+outname+'_'+pola+'_ant_azi')
    close_file(outpath+outname+'_'+pola+'_azi')
    close_file(outpath+outname+'_'+pola+'_time')
    close_file(outpath+outname+'_'+pola+'_kp')

print('DONE')

