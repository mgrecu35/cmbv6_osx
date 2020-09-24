function callcmb(nodes,stormType,nmemb,nmu,brs,bst,hzero,bzd,bcf,bbPeak,sfcType,sysdN,
    zku,zka,relFlag,piaSRT,retParam,radarData,radarRet,i,j,nmfreq,locZAngle,nwdm)

    #println(nmemb, "  ",nmfreq, " ",nmemb)
    nodes=Cint.(zeros(5))
    nodes[1]=Int32(trunc(bst[i,j]/2))
    nodes[3]=Int32(trunc((bzd[i,j])/2))
    if bbPeak[i,j]>0
        nodes[3]=Int32(trunc(bbPeak[i,j]/2))
    end
    nodes[2]=nodes[3]-3
    nodes[4]=nodes[3]+2
    nodes[5]=Int32(trunc(bcf[i,j]/2))
    if nodes[1]>nodes[2]
        nodes[1]=nodes[2]-3
    end
    #println(nodes)
    z13obs=Cfloat.(zeros(88))
    for k=1:88
        z13obs[k]=Cfloat(10.0*log10((10^(0.1*zku[i,j,2*k])+10^(0.1*zku[i,j,2*k-1]))/2))
    end
    reliabFpoint=Cint.(zeros(1))
    reliabFpoint[1]=Cint(relFlag[i,j]*0)
    #println(reliabFpoint,size(relFlag))
    ccall((:setrelflag_,"./combAlg"),Cvoid,(Ref{Cint},),reliabFpoint)

    pia13srt=Cfloat(piaSRT[i,j])
    z13true=Cfloat.(zeros(88).-99)
    z35true=Cfloat.(zeros(88).-99)
    z13obs=Cfloat.(z13obs)
    pia13ret=Cfloat(-99)
    pia35ret=Cfloat(-99)
    ic,jc=Cint(1),Cint(1)
    z35Sim=Cfloat.(zeros(88).-99)
    lwc_ret=Cfloat.(zeros(88))
    sysdN1=sysdN-0.1
    log10dNw=Cfloat.(0.00005*randn(88).+sysdN1)
    dr=Cfloat(0.25)
    nodes=Cint.(nodes.+1)
    isurf=Cint(trunc(brs[i,j]/2)+1)
    imu=Cint.(zeros(88).+3)
    ngates=Cint(88)
    nmfreqm=Cint(8)
    hh=Cfloat.((88:-1:1).*0.25)
    itype=Cint(1)
    kext,salb,asym=Cfloat.(zeros(88*8)),Cfloat.(zeros(88*8)),Cfloat.(zeros(88*8))
    rrate,dm=Cfloat.(zeros(88)),Cfloat.(zeros(88))
    imembC=Cint(1)
    hfreez=Cfloat((175-bzd[i,j])*0.125)
    ccall((:fhb11_,"./combAlg"),Cvoid,(Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat}, # z13true,z35true,z13obs,pia13ret,pia35ret,
    Ref{Cint},Ref{Cint}, #ic,jc
    Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat}, #z35Sim,lwc_ret,log10dNw,dr,
    Ref{Cint},Ref{Cint},Ref{Cint},Ref{Cint},Ref{Cint}, #nodes,isurf,imu,ngates,nmfreqm,
    Ref{Cfloat}, #hh
    Ref{Cint},   #itype
    Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat}, #kext,salb,asym,rrate,dm,hfreez,pia13srt,
    Ref{Cint}), #imembC
    z13true,z35true,z13obs,pia13ret,pia35ret,
    ic,jc,
    z35Sim,lwc_ret,log10dNw,dr,
    nodes,isurf,imu,ngates,nmfreqm,
    hh,
    itype,
    kext,salb,asym,rrate,dm,hfreez,pia13srt,
    imembC)
    for i=1:88
        if dm[i]>0.5
            log10dNw[i]-=0.15*(dm[i]-1.5);
            i0dm=Int(trunc((dm[i]-0.5)/0.04));
            i0dm=max(0,i0dm)
            i0dm=min(59,i0dm)
            log10dNw[i]+=0.2*(nwdm[i0dm+1,3]-log10dNw[i]);
        end
    end
    if dm[nodes[5]]>0.5
        log10dNw.-=0.15*(dm[nodes[5]]-1.5);
        i0dm=Int(trunc((dm[nodes[5]]-0.5)/0.04));
        i0dm=max(0,i0dm)
        i0dm=min(59,i0dm)
        log10dNw.+=0.2*(nwdm[i0dm,3].-log10dNw);
        log10dNw=Cfloat.(log10dNw);
    end
    log10dNw.+=0.1
    ccall((:fhb11_,"./combAlg"),Cvoid,(Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat}, # z13true,z35true,z13obs,pia13ret,pia35ret,
    Ref{Cint},Ref{Cint}, #ic,jc
    Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat}, #z35Sim,lwc_ret,log10dNw,dr,
    Ref{Cint},Ref{Cint},Ref{Cint},Ref{Cint},Ref{Cint}, #nodes,isurf,imu,ngates,nmfreqm,
    Ref{Cfloat}, #hh
    Ref{Cint},   #itype
    Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat}, #kext,salb,asym,rrate,dm,hfreez,pia13srt,
    Ref{Cint}), #imembC
    z13true,z35true,z13obs,pia13ret,pia35ret,
    ic,jc,
    z35Sim,lwc_ret,log10dNw,dr,
    nodes,isurf,imu,ngates,nmfreqm,
    hh,
    itype,
    kext,salb,asym,rrate,dm,hfreez,pia13srt,
    imembC)
    #println(lwc_ret)
    #exit(-1)
    #a=findall(rrate.>0)
    #lwc_ret[a].=10 .^lwc_ret[a]
    return rrate[nodes[5]],rrate[nodes[2]], lwc_ret, dm, log10dNw
end