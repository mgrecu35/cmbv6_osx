function callcmb(nodes,stormType,nmemb,nmu,brs,hzero,bzd,bcf,sfcType,sysdN,
    zku,zka,relFlag,piaSRT,retParam,radarData,radarRet,i,j,nmemb,nmfreq,locZAngle)
    randemiss=Cfloat.(zeros(nmemb*nmfreq*2))
    rms1=Cfloat.(zeros(1))
    rms2=Cfloat.(zeros(1))
    dn=Cfloat.([-0.25])
    nodes[1]=Int32(trunc(bst[i,j,1]/2))
    nodes[3]=Int32(trunc((bzd[i,j])/2))
    nodes[2]=nodes[3]-2
    nodes[4]=nodes[3]+2
    nodes[5]=Int32(trunc(bcf[i,j]/2))
    for k=1:88
        z13obs[k]=Cfloat(10.0*log10((10^(0.1*zku[i,j,2*k])+10^(0.1*zku[i,j,2*k-1]))/2))
    end
    stormType.iSurf=Cint(trunc(brs[i,j]/2))
    stormType.freezH=Cfloat(hzero[i,j]/1000.)
    stormType.nodes=pointer(nodes)
    logdnwf=Array{Cfloat}(undef,9*nmemb)
    xscalev=Cfloat.(zeros(nmemb))
    i0=Cint.([1])
    j0=Cint.([1])
    iit=Cint.([3])
    i300=i%300+1
    j1=Cint.([j])
    i1= Cint.([i300])
    nmu=Int32.([5])
    msFlag=Cint.([1])
    dZms=Cfloat.(zeros(nmemb))
    ccall((:__geophysens_MOD_interpoldnw,"./combAlg"),Cvoid,(Ref{Cint},Ref{Cint},Ref{Cfloat}),
        j1,i1,logdnwf)

    if sfcType[i,j]==0
        wfractPix=Cfloat.([100])
    else
        wfractPix=Cfloat.([1])
    end
    radarRet.logdNw=pointer(logdnwf.+sysdN)
    radarData.z13obs=pointer(z13obs)
    reliabFpoint=Cint.(zeros(1))
    reliabFpoint[1]=Cint(relFlag[i,j])
    ccall((:setrelflag_,"./combAlg"),Cvoid,(Ref{Cint},),reliabFpoint)
    radarData.pia13srt=piaSRT[i,j]
    localZAngle=locZAngle[i,j]
    ccall((:ensradJulia, "./combAlg"),
        Cvoid, (Ref{radarDataType},
        Ref{stormStructType},Ref{retParamType},Ref{Cint},
        Ref{radarRetType},
        Ref{Cint},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},
        Ref{Cint},Ref{Cfloat},Ref{Cfloat},
        Ref{Cfloat},Ref{Cfloat},Ref{Clong},Ref{Cint},
        Ref{Cint},Ref{Cfloat},Ref{Cint}),
        radarData,stormType,retParam,nmu,radarRet,
        iflag,rms1,rms2,dn,iit,
        xscalev,randemiss,localZAngle,wfractPix,ichunk,
        i0,j0,dZms,msFlag)
end
