using Statistics
np=pyimport("numpy")
function callcmb(brs,bst,hzero,bzd,bcf,bbPeak,sfcType,sysdN,pType,
    zku,zka,relFlag,piaSRT,nmfreq,nwdm,iBB)

    #println(nmemb, "  ",nmfreq, " ",nmemb)
    nodes=Cint.(zeros(5))
    nodes[1]=Int32(trunc(bst/2))
    nodes[3]=Int32(trunc((bzd)/2))
    if bbPeak>0 && iBB==1
        nodes[3]=Int32(trunc(bbPeak/2))
    end
    nodes[2]=nodes[3]-4
    nodes[4]=nodes[3]+2
    nodes[5]=Int32(trunc(bcf/2))
    if nodes[1]>nodes[2]
        nodes[1]=nodes[2]-3
    end
    #println(nodes)
    z13obs=Cfloat.(zeros(88))
    z35obs=Cfloat.(zeros(88))
    for k=1:88
        z13obs[k]=Cfloat(10.0*log10((10^(0.1*zku[2*k])+10^(0.1*zku[2*k-1]))/2))
        z35obs[k]=Cfloat(10.0*log10((10^(0.1*zka[2*k])+10^(0.1*zka[2*k-1]))/2))-2.0
    end
    reliabFpoint=Cint.(zeros(1))
    reliabFpoint[1]=Cint(relFlag*0)
    #println(reliabFpoint,size(relFlag))
    ccall((:setrelflag_,"./combAlg"),Cvoid,(Ref{Cint},),reliabFpoint)

    pia13srt=Cfloat(piaSRT)
    z13true=Cfloat.(zeros(88).-99)
    z35true=Cfloat.(zeros(88).-99)
    z13obs=Cfloat.(z13obs)
    pia13ret=Cfloat(-99)
    pia35ret=Cfloat(-99)
    ic,jc=Cint(1),Cint(1)
    z35Sim=Cfloat.(zeros(88).-99)
    lwc_ret=Cfloat.(zeros(88))
    kext,salb,asym=Cfloat.(zeros(88*8)),Cfloat.(zeros(88*8)),Cfloat.(zeros(88*8))
    rrate,dm=Cfloat.(zeros(88)),Cfloat.(zeros(88))

    z13true1=Cfloat.(zeros(88).-99)
    z35true1=Cfloat.(zeros(88).-99)
    pia13ret1=Cfloat(-99)
    pia35ret1=Cfloat(-99)
    z35Sim1=Cfloat.(zeros(88).-99)
    lwc_ret1=Cfloat.(zeros(88))
    kext1,salb1,asym1=Cfloat.(zeros(88*8)),Cfloat.(zeros(88*8)),Cfloat.(zeros(88*8))
    rrate1,dm1=Cfloat.(zeros(88)),Cfloat.(zeros(88))

    sysdN1=sysdN-0.1
    log10dNw=Cfloat.(0.00005*randn(88).+sysdN1)
    dr=Cfloat(0.25)
    nodes=Cint.(nodes.+1)
    isurf=Cint(trunc(brs/2)+1)
    imu=Cint.(zeros(88).+3)
    ngates=Cint(88)
    nmfreqm=Cint(8)
    hh=Cfloat.((88:-1:1).*0.25)
    itype=Cint(pType)
    imembC=Cint(3)
    hfreez=Cfloat((175-bzd)*0.125)
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
    if iBB==0
        xs=0.95
    else
        xs=1.0
    end
    println(nodes)
    println(typeof(log10dNw))

    adj=1
    if  adj==1
        for i=1:88
            if dm[i]>0.5 && (i<nodes[2] || i>nodes[4])
                log10dNw[i]-=0.15*xs*(dm[i]-1.5);
                i0dm=Int(trunc((dm[i]-0.5)/0.04));
                i0dm=max(0,i0dm)
                i0dm=min(59,i0dm)
                log10dNw[i]+=0.2*xs*(nwdm[i0dm+1,3]-log10dNw[i]);
            end
        end
        if dm[nodes[5]]>0.5
            log10dNw.-=0.15*xs*(dm[nodes[5]]-1.5);
            i0dm=Int(trunc((dm[nodes[5]]-0.5)/0.04));
            i0dm=max(0,i0dm)
            i0dm=min(59,i0dm)
            log10dNw.+=0.2*xs*(nwdm[i0dm,3].-log10dNw);
            log10dNw=Cfloat.(log10dNw);
        end
        log10dNw.+=0.2
        if iBB==0
            log10dNw.+=0.6
        end
        log10dNw[nodes[1]:nodes[2]].+=0.01;
        for k=nodes[1]:nodes[3]
            if z35obs[k]>10 && z13obs[k]>0
                log10dNw[k]+=0.03*(z35obs[k]-z35Sim[k])
            end
        end
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
    #println(z35Sim)
    #println(lwc_ret)
    #exit(-1)
    #a=findall(rrate.>0)
    #lwc_ret[a].=10 .^lwc_ret[a]
    end
    nmemb=150

    dnw=zeros(88)
    gradZ=zeros(88,88)
    log10dNwref=copy(log10dNw)
    for it=1:2
    for ik=nodes[1]:nodes[5]
        log10dNw=log10dNwref.+dnw
        log10dNw[ik]=log10dNw[ik]+0.1
        log10dNw=Cfloat.(log10dNw)
        ccall((:fhb11_,"./combAlg"),Cvoid,(Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat}, # z13true,z35true,z13obs,pia13ret,pia35ret,
        Ref{Cint},Ref{Cint}, #ic,jc
        Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat}, #z35Sim,lwc_ret,log10dNw,dr,
        Ref{Cint},Ref{Cint},Ref{Cint},Ref{Cint},Ref{Cint}, #nodes,isurf,imu,ngates,nmfreqm,
        Ref{Cfloat}, #hh
        Ref{Cint},   #itype
        Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat}, #kext,salb,asym,rrate,dm,hfreez,pia13srt,
        Ref{Cint}), #imembC
        z13true1,z35true1,z13obs,pia13ret1,pia35ret1,
        ic,jc,
        z35Sim1,lwc_ret1,log10dNw,dr,
        nodes,isurf,imu,ngates,nmfreqm,
        hh,
        itype,
        kext1,salb1,asym1,rrate1,dm1,hfreez,pia13srt,
        imembC)
        gradZ[ik,:]=(z35Sim1-z35Sim)./(0.1)
    end
    a=findall(z13obs.>10)
    b=findall(z35obs[a].>10)
    nobs=size(b)[1]
    invCov=0.25*np.eye(88)
    gradZ1=gradZ[nodes[1]:nodes[5],a[b]]
    dy=gradZ1[:,:]*(z35obs[a[b]]-z35Sim[a[b]])-
    invCov[nodes[1]:nodes[5],nodes[1]:nodes[5]]*(dnw[nodes[1]:nodes[5]]);
    A=gradZ1*transpose(gradZ1)+invCov[nodes[1]:nodes[5],nodes[1]:nodes[5]];
    dx=A\dy;
    dnw[nodes[1]:nodes[5]].=dnw[nodes[1]:nodes[5]].+0.5*dx;
    log10dNw[nodes[1]:nodes[5]].=log10dNwref[nodes[1]:nodes[5]].+dnw[nodes[1]:nodes[5]]

    log10dNw=Cfloat.(log10dNw)
    println(typeof(log10dNw),size(log10dNw))
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
end
    z35Sim=Float64.(copy(z35Sim))


    return rrate[nodes[5]],rrate[nodes[2]], rrate,lwc_ret, dm, log10dNw,z35Sim,z35Sim,z13obs,z35obs,nodes
end
