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

zku,vku,lat,lon,r,alt,roll,t=readiphexKu(fname)
from pyproj import Geod
ws_geod=Geod(ellps='WGS84')
d=[]
dist=0
for i in range(lat.shape[0]-1):
    d.append(dist)
    az1,az2,ddist=ws_geod.inv(lon[i],lat[i],lon[i+1],lat[i+1])
    dist+=ddist/1000.
a=nonzero(abs(roll)>2.5)
for i in a[0]:
    zku[i,:]=-99
h=alt.mean()/1000.-r/1000.
zkum=ma.array(zku,mask=zku<-10)
vkum=ma.array(vku,mask=zku<-10)
vkum=ma.array(vku,mask=zku!=zku)
nx=400
for i in range(3,7):
    plt.figure()
    plt.subplot(111)
    plt.pcolormesh(d[i*1000+nx:i*1000+nx+1500], h, \
                   zkum[i*1000+nx:(i)*1000+1500+nx].T,cmap='jet',vmin=0,vmax=45)
    plt.title('Ku-band Reflectivity')
    c=plt.colorbar()
    c.ax.set_title('dBZ')
    plt.ylim(0,15)
    plt.savefig('Z_Kuband_2335_2342.png')

fname='IPHEX_HIWRAP_L1B_20140612-215738-20140612-225747_HKa_dist_v01.nc'
fname='IPHEX_HIWRAP_L1B_20140612-215723-20140612-225846_HKu_dist_v01.nc'
zku,vku,lat,lon,r,alt,roll,t=readiphexKu(fname)
h=alt.mean()/1000.-r/1000.
a=nonzero(abs(roll)>2.5)
for i in a[0]:
    zku[i,:]=-99
zkum=ma.array(zku,mask=zku<-10)
vkum=ma.array(vku,mask=zku<-10)
vkum=ma.array(vku,mask=zku!=zku)
nx=400

#IPHEX_HIWRAP_L1B_20140612-225733-20140612-232959_HKu_dist_v01.nc
d=[]
dist=0
for i in range(lat.shape[0]-1):
    d.append(dist)
    az1,az2,ddist=ws_geod.inv(lon[i],lat[i],lon[i+1],lat[i+1])
    dist+=ddist/1000.
for i in range(3,7):
    plt.figure()
    plt.subplot(111)
    plt.pcolormesh(d[i*1000+nx:i*1000+nx+1500], h, \
                   zkum[i*1000+nx:(i)*1000+1500+nx].T,cmap='jet',vmin=0,vmax=45)
    plt.title('Ku-band Reflectivity')
    c=plt.colorbar()
    c.ax.set_title('dBZ')
    plt.ylim(0,15)
    plt.savefig('Z_Kuband_2335_2342.png')
   
stop
cfadv=zeros((80,50),float)
vm=zeros((80),float)
im=zeros((80),float)
for i in range(5000,8000):
    if abs(roll[i])<4:
        for k in range(1000):
            if h[k]>.250:
                i0=int(h[k]/.250)
                if h[k]>5:
                    vt=-3.3993+0.186*zku[i,k]
                    vt=2.6*10.**(0.1*zku[i,k]*0.107)
                else:
                    vt=2.6*10.**(0.1*zku[i,k]*0.107)
                j0=int(vku[i,k]-vt+25)
                
                if i0<80 and i0>0 and j0>=0 and j0<50:
                    cfadv[i0,j0]+=1
                    if(zku[i,k]>5):
                        vm[i0]+=vku[i,k]-vt
                        im[i0]+=1
                    
plt.pcolormesh(-25+arange(50),arange(80)*0.25,cfadv, norm=LogNorm(),\
               cmap='jet')

plt.figure()
for i in range(30):
    plt.plot(vku[6400+i*10,300:1100],h[300:1100])

plt.xlim(-20,20)

plt.figure()
a=nonzero(im>0)
vm[a]/=im[a]
plt.plot(vm,arange(80)*0.25)
plt.xlim(-5,10)
plt.ylim(1,16)
stop
#fname='IPHEX_HIWRAP_L1B_2014611-195739-2014611-210300_HKa_dist_v01.h5'
#zka,latKa,lonKa,rKa,altKa,rollKa,tKa=readKa(fname,250)

h=alt.mean()/1000.-r/1000.
hKa=alt.mean()/1000.-r/1000.
from numpy import *
zkum=ma.array(zku.T,mask=zku.T<0)
zkam=ma.array(zka.T,mask=zka.T<0)
from pyproj import Geod
ws_geod=Geod(ellps='WGS84')
d=[]
d2=[]
dKa=[]
dist=0
distKa=0
it=10
for i in range(it*1000,it*1000+1000):
    d.append(dist)
    az1,az2,ddist=ws_geod.inv(lon[i],lat[i],lon[i+1],lat[i+1])
    dist+=ddist/1000.
    d2.append(distKa)
    az1,az2,ddist=ws_geod.inv(lonKa[i-49],latKa[i-49],lonKa[i-48],latKa[i-48])
    distKa+=ddist/1000.
matplotlib.rcParams.update({'font.size': 12})
plt.suptitle('IPHEX June 11, 2014')
plt.subplot(211)
plt.pcolormesh(d, h, zkum[10000:11000].T,cmap='jet',vmin=0,vmax=50)
plt.ylim(0,12)
plt.ylim(0,15)
plt.ylabel('Height (km)')
#plt.xlabel('Distance (km)')
cbar=plt.colorbar()
plt.subplot(212)
plt.pcolormesh(d2, hKa, zkam[10000-49:11000-49].T,cmap='jet',vmin=0,vmax=50)
plt.ylim(0,12)
plt.ylim(0,15)
plt.ylabel('Height (km)')
plt.xlabel('Distance (km)')
cbar=plt.colorbar()
cbar.ax.set_title('dBZ')
plt.savefig('IPHEXJune112014.png')
stop
for it in range(24):
    d=[]
    dist=0
    for i in range(it*1000,it*1000+1000):
        d.append(dist)
        az1,az2,ddist=ws_geod.inv(lon[i],lat[i],lon[i+1],lat[i+1])
        dist+=ddist/1000.
        #print ddist, dist
    plt.figure()
    plt.subplot(211)
    plt.pcolormesh(d,h,zkum[it*1000:it*1000+1000,:].T,vmin=10,vmax=50,cmap='jet')
    plt.xlim(55,110)
    plt.ylim(0,15)
    plt.ylabel('Height (km)')
    plt.colorbar()
    plt.subplot(212)
    plt.pcolormesh(d,h,zkam[it*1000:it*1000+1000,:].T,vmin=10,vmax=50,cmap='jet')
    plt.xlim(55,110)
    plt.ylim(0,15)
    plt.ylabel('Height (km)')
    plt.xlabel('Distance (km)')
    cbar=plt.colorbar()
    cbar.ax.set_title('dBZ')
    plt.savefig('leg%2.2i_IPHEXKu.png'%it)
