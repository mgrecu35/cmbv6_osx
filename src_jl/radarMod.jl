module radarMod
export  Δr, rainFract
export Z_Obs,pia_Obs
export invCov
export dNw
export xMean
export nClutF
export iSurf
export gn_srtPIA, gn_relFlag, wpia
function setRadar_dr(dr,n)
    global Δr=dr;
    global rainFract=zeros(n);
    global Z_Obs=zeros(n);
    global invCov=zeros(n,n);
    global xMean=zeros(n);
    global dNw=zeros(n);
end
function setSRTPIA(srtPIAv,relFlagv,wpiav)
    global gn_srtPIA=srtPIAv
    global gn_relFlag=relFlagv
    global wpia=wpiav
end
function setClutFAndSurf(nClutF_Value,iSurf_Value)
    global nClutF=nClutF_Value
    global iSurf=iSurf_Value
end
function setRainFract(fract)
    global rainFract.=fract;
end
function setZObs(z_obs)
    global Z_Obs.=z_obs;
end
function setPIAObs(PIA)
    global pia_Obs=PIA;
end
function setInvCov(mat,xm)
    global invCov.=mat;
    global xMean.=xm;
end
function setdNw(dnwVal)
    global dNw.=dnwVal
end
end
