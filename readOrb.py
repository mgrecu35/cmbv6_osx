import numpy as np
import time
from datetime import date
from numpy import *

from netCDF4 import Dataset

def readOrb(orb):
    fh=Dataset(orb)
    Lon=fh['NS/Longitude'][:,24]
    Lat=fh['NS/Latitude'][:,24]
    a=nonzero((Lon+85)*(Lon+75)<0)
    b=nonzero((Lat[a]-35)*(Lat[a]-45)<0)
    sfcPrecip=fh['NS/SLV/precipRateNearSurface'][a[0][b],:]
    lon=fh['NS/Longitude'][a[0][b],:]
    lat=fh['NS/Latitude'][a[0][b],:]
    zsfc=fh['NS/PRE/elevation'][a[0][b],:]
    hzero=fh['NS/VER/heightZeroDeg'][a[0][b],:]
    a1=nonzero(zsfc>300)
    bcf=fh['NS/PRE/binClutterFreeBottom'][a[0][b],:]
    brs=fh['NS/PRE/binRealSurface'][a[0][b],:]
    bst=fh['NS/PRE/binStormTop'][a[0][b],:]
    bzd=fh['NS/VER/binZeroDeg'][a[0][b],:]
    zku=fh['NS/PRE/zFactorMeasured'][a[0][b],:,:]
    zka=fh['MS/PRE/zFactorMeasured'][a[0][b],:,:]
    bbPeak=fh['NS/CSF/binBBPeak'][a[0][b],:]
    t=fh['NS/ScanTime/Hour'][a[0][b]]+fh['NS/ScanTime/Minute'][a[0][b]]/60.0
    pType=fh['NS/CSF/typePrecip'][a[0][b],:]
    pType=(pType/1e7).astype(int)
    sfcType=fh['NS/PRE/landSurfaceType'][:,:]
    reliabF=fh['/NS/SRT/reliabFlag'][:,:]
    piaKu=fh['/NS/SRT/pathAtten'][:,:]
    locZAngle=f['/NS/PRE/localZenithAngle'][:,:]
    #print(fh['NS/PRE'])
    return zku,bzd,bst,brs,bcf,zsfc,hzero,lon,lat,fh,t,sfcPrecip,pType,zka,bbPeak,\
        sfcType,reliabF,piaKu,locZAngle
