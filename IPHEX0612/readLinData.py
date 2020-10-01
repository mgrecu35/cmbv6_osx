from numpy import *
from scipy.io import netcdf as cdf
from netCDF4 import Dataset
def readKuKa(fname,n1):
    f=cdf.netcdf_file(fname,'r')
    r=f.variables['range'][:]
    zku=f.variables['zku'][:,:]
    zka=f.variables['zka'][:,:]
    lat=f.variables['lat'][:]
    lon=f.variables['lon'][:]
    alt=f.variables['altitude'][:]
    roll=f.variables['roll'][:]
    t=f.variables['timed'][:]
    return zku,zka,lat,lon,r,alt,roll,t

def readKa(fname,n1):
    f=Dataset(fname,'r')
    r=f.variables['rangevec'][:]
    alt=f.variables['height'][:]
    zku=f.variables['stitchedReflectivity'][:,:]
    lat=f.variables['lat'][:]
    lon=f.variables['lon'][:]
    alt=f.variables['height'][:]
    roll=f.variables['roll'][:]
    t=f.variables['timeUTC'][:]
    return zku,lat,lon,r,alt,roll,t

def readiphex(fname):
    f=cdf.netcdf_file(fname,'r')
    r=f.variables['range'][:]
    dist=f.variables['dist'][:]
    zx=f.variables['zx'][:,:]
    zka=f.variables['zka'][:,:]
    zku=f.variables['zku'][:,:]
    
    vx=f.variables['dopx'][:,:]
    vku=f.variables['dopku'][:,:]
    vka=f.variables['dopka'][:,:]
    vw=f.variables['dopw'][:,:]
    zw=f.variables['zw'][:,:]
    zx=zx-0.5
    lat=f.variables['lat'][:]
    lon=f.variables['lon'][:]
    alt=f.variables['altitude'][:]
    return zx,zka,zku,zw,vx,vku,vka,vw,lat,lon,r,alt,dist

def readiphexX(fname):
    f=cdf.netcdf_file(fname,'r')
    r=f.variables['range'][:]
    vDop=f.variables['dopcorr'][:,:]
    zku=f.variables['zku'][:,:]
    lat=f.variables['lat'][:]
    lon=f.variables['lon'][:]
    alt=f.variables['altitude'][:]
    roll=f.variables['roll'][:]
    t=f.variables['timed'][:]
    return zku,vDop,lat,lon,r,alt,roll,t


def readiphexKu(fname):
    f=Dataset(fname,'r')
    r=f.variables['range'][:]
    vDop=f.variables['dopcorr'][:,:]
    zku=f.variables['zku'][:,:]
    lat=f.variables['lat'][:]
    lon=f.variables['lon'][:]
    alt=f.variables['altitude'][:]
    roll=f.variables['roll'][:]
    t=f.variables['timed'][:]
    return zku,vDop,lat,lon,r,alt,roll,t

