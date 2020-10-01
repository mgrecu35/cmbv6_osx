from readLinData import *
from   numpy import *
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colors as col
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator,LogFormatter
from matplotlib.colors import BoundaryNorm
from matplotlib.colors import LogNorm
import pickle
#fname='mc3e_hiwrap_20110520_134702-181402.nc'
#fname='mc3e_hiwrap_20110425_070444-115447_chirp.nc'
fname='IPHEX_HIWRAP_L1B_2014611-195727-2014611-210243_HKu_dist_v01.h5'
fname='IPHEX_EXRAD_L1B_20140523-225808-20140523-235906_nadir_dist_v01.nc'
fname='IPHEX_HIWRAP_L1B_2014523-225723-2014523-235850_HKu_dist_v01.h5'
fname='IPHEX_HIWRAP_L1B_20140523-215716-20140523-230217_HKu_dist_v01.nc'
fname='IPHEX_HIWRAP_L1B_20140523-225723-20140523-235850_HKu_dist_v01.nc'
fname='IPHEX_HIWRAP_L1B_20140612-215738-20140612-225747_HKa_dist_v01.nc'

zka,vka,latKa,lonKa,rKa,altKa,rollKa,tKa=readiphexKu(fname)
from pyproj import Geod
ws_geod=Geod(ellps='WGS84')
d=[]
dist=0
for i in range(latKa.shape[0]-1):
    d.append(dist)
    az1,az2,ddist=ws_geod.inv(lonKa[i],latKa[i],lonKa[i+1],latKa[i+1])
    dist+=ddist/1000.
a=nonzero(abs(rollKa)>2.5)
for i in a[0]:
    zka[i,:]=-99
hKa=altKa.mean()/1000.-rKa/1000.
zkam=ma.array(zka,mask=zka<-10)
vkam=ma.array(vka,mask=zka<-10)
vkam=ma.array(vka,mask=zka!=zka)
zKuL=[]
zKaL=[]
zXL=[]
dL=[]
h1=arange(176)*0.125
fname='IPHEX_HIWRAP_L1B_20140612-215723-20140612-225846_HKu_dist_v01.nc'
zku,vku,lat,lon,r,alt,roll,t=readiphexKu(fname)

fname='IPHEX_EXRAD_L1B_20140612-220807-20140612-230145_nadir_dist_v01.nc'
zx,vx,lat,lon,rx,altx,rollx,tx=readiphexKu(fname)
h=alt.mean()/1000.-r/1000.
hx=altx.mean()/1000.-rx/1000.
dL=[]
for i in range(3000,12000):
    if abs(rollKa[i])<2.0:
        zkaint=interp(h1,hKa[::-1],zka[i,:][::-1])
        iku=argmin(abs(t-tKa[i]))
        zkuint=interp(h1,h[::-1],zku[iku,:][::-1])
        ix=argmin(abs(tx-tKa[i]))
        zxint=interp(h1,hx[::-1],zx[ix,:][::-1])
        zKuL.append(zkuint)
        zXL.append(zxint)
        zKaL.append(zkaint)
        dL.append(d[i])
zKuL=array(zKuL)
zXL=array(zXL)
zKaL=array(zKaL)
dL=array(dL)
import pickle
pickle.dump([zKuL,zKaL,zXL,dL],open('iphex0612.pklz','wb'))
hint=h1
x1=100
x2=400
x1=230
x2=250
for i in range(8):
    plt.figure()
    x1=100+i*60
    x2=100+i*60+60
    plt.subplot(311)
    plt.pcolormesh(dL,hint,zKuL.T,vmin=0,vmax=50,cmap="jet")
    plt.xlim(x1,x2)
    plt.colorbar()
    plt.ylim(0,10)
    plt.subplot(312)
    plt.pcolormesh(dL,hint,(zKaL).T,vmin=0,vmax=40,cmap="jet")
    plt.ylim(0,10)
    plt.xlim(x1,x2)
    plt.colorbar()
    plt.subplot(313)
    plt.pcolormesh(dL,hint,zXL.T,vmin=0,vmax=50,cmap="jet")
    plt.ylim(0,10)
    plt.xlim(x1,x2)
    plt.colorbar()
