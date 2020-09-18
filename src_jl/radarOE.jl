using .scatTables

using LinearAlgebra



scipy=pyimport("scipy.ndimage") 


using .radarMod

function simZ(Dm)
    drk=Δr
    fract=rainFract
    nz=size(Dm)[1]
    piaKu=0
    piaKa=0.0
    zSim=zeros(nz).-99.9
    rrateL=zeros(nz)#.-99.9
    noise=dNw
    attKu=0.0
    for i=1:nClutF
        if Z_Obs[i]<5 || Z_Obs[i]!=Z_Obs[i]
            continue
        end
        n1s,n2s=bisection(dmTs,Dm[i])
        n1,n2=bisection(dmT,Dm[i])
        att=(1-fract[i])*attKuTs[n1s]+fract[i]*attKuT[n1]
        att=att*10^noise[i]
        attKa=(1-fract[i])*attKaTs[n1s]+fract[i]*attKaT[n1]
        attKa=attKa*10^noise[i]
        piaKu=piaKu+att*drk
        zSim[i]=log10((1-fract[i])*10^(0.1*zKuTs[n1s])+(fract[i])*10^(0.1*zKuT[n1]))*10.0-piaKu+10*noise[i]
        rrateL[i]=10^noise[i]*((1-fract[i])*rateTs[n1s]+fract[i]*rateT[n1])
        piaKu+=att*drk
        piaKa+=attKa*drk
    end
    piaKu=piaKu+attKu*drk*(iSurf-nClutF)
    return piaKu,rrateL,zSim
end

function simMZ(Dm)
    drk=Δr
    fract=rainFract
    nz=size(Dm)[1]
    piaKu=0
    piaKa=0.0
    piaW=0.0
    zSim=zeros(nz).-99.9
    zSimKa=zeros(nz).-99.9
    zSimW=zeros(nz).-99.9
    rrateL=zeros(nz)#.-99.9
    noise=dNw
    for i=1:nClutF
        if Z_Obs[i]<5 
            continue
        end
        if Z_Obs[i]!=Z_Obs[i]
            continue
        end
        n1s,n2s=bisection(dmTs,Dm[i])
        n1,n2=bisection(dmT,Dm[i])
        att=(1-fract[i])*attKuTs[n1s]+fract[i]*attKuT[n1]
        att=att*10^noise[i]
        attKa=(1-fract[i])*attKaTs[n1s]+fract[i]*attKaT[n1]
        attKa=attKa*10^noise[i]
        piaKu=piaKu+att*drk
        piaKa=piaKa+attKa*drk
        attW=(1-fract[i])*attWTs[n1s]+fract[i]*attWT[n1]
        attW=attW*10^noise[i]
        piaW=piaW+attW*drk
        zSim[i]=log10((1-fract[i])*10^(0.1*zKuTs[n1s])+(fract[i])*10^(0.1*zKuT[n1]))*10.0-piaKu+10*noise[i]
        zSimKa[i]=log10((1-fract[i])*10^(0.1*zKaTs[n1s])+(fract[i])*10^(0.1*zKaT[n1]))*10.0-piaKa+10*noise[i]
        zSimW[i]=log10((1-fract[i])*10^(0.1*zWTs[n1s])+(fract[i])*10^(0.1*zWT[n1]))*10.0-piaW+10*noise[i]

        rrateL[i]=10^noise[i]*((1-fract[i])*rateTs[n1s]+fract[i]*rateT[n1])
        #println(noise[i]," ",fract[i]," ",n1s)
        piaKu+=att*drk
        piaKa+=attKa*drk
        piaW+=attW*drk
    end
    return piaKu,piaKa,piaW,rrateL,zSim,zSimKa,zSimW
end

function simZ_g(Dm)
    drk=Δr
    fract=rainFract
    nz=size(Dm)[1]
    piaKu=0
    piaKa=0
    zSim=zeros(nz).-99.9
    rrateL=zeros(nz).-99.9
    noise=dNw
    dZdDm=zeros(nz+1,nz)
    dpiaKu1=0
    for i=1:nClutF
        if Z_Obs[i]<5 || Z_Obs[i]!=Z_Obs[i]
            continue
        end
        n1s,n2s=bisection(dmTs,Dm[i])
        n1,n2=bisection(dmT,Dm[i])
        if n1==240
            n1=239
        end
        if n1s==240
            n1s=239
        end
       
        att=(1-fract[i])*attKuTs[n1s]+fract[i]*attKuT[n1]
        attKa=(1-fract[i])*attKaTs[n1s]+fract[i]*attKaT[n1]
        n2s1=n2s
        n21=n2
        if n2s==1
            n2s1=2
        end
        if n2==1
            n21=1
        end
        dDm=(1-fract[i])*dmTs[n2s1]+fract[i]*dmT[n21]-
        ((1-fract[i])*dmTs[n1s]+fract[i]*dmT[n1])+1e-3
        att2=(1-fract[i])*attKuTs[n2s1]+fract[i]*attKuT[n21]
        att=att*10^noise[i]
        attKa=attKa*10^noise[i]
        att2=att2*10^noise[i]
        piaKu=piaKu+att*drk
        piaKa=piaKa+attKa*drk
        dpiaKu1=(att2-att)/dDm*Δr
        zSim[i]=log10((1-fract[i])*10^(0.1*zKuTs[n1s])+(fract[i])*10^(0.1*zKuT[n1]))*10.0-piaKu+10*noise[i]
        zSim2=log10((1-fract[i])*10^(0.1*zKuTs[n2s1])+(fract[i])*10^(0.1*zKuT[n21]))*10.0-piaKu+10*noise[i]
        dZdDm1=(zSim2-zSim[i])/dDm
        rrateL[i]=10^noise[i]*((1-fract[i])*rateTs[n1s]+fract[i]*rateT[n1])
        piaKu+=att*drk
        piaKa=piaKa+attKa*drk
        dZdDm[i,i]=dZdDm[i,i]+dZdDm1-dpiaKu1
        if i+1<=nClutF
            dZdDm[i+1:nClutF,i]=dZdDm[i+1:nClutF,i].-2*dpiaKu1
        end
        dZdDm[end,i]=dZdDm[end,i]+2*sqrt(wpia)*dpiaKu1
    end
    dZdDm[end,nClutF]=dZdDm[end,nClutF]+2*sqrt(wpia)*dpiaKu1*(iSurf-nClutF)
    return piaKu,piaKa,rrateL,zSim,dZdDm
end








function f_z(Dm)
    piaKu,rrate,zSim=simZ(Dm)
    n=size(zSim)[1]
    fobj=0.0
    for i=1:nClutF
        if Z_Obs[i]>5
            fobj=fobj+(zSim[i]-Z_Obs[i])^2
        end
    end
    #println(fobj)
    for i=1:nClutF
        for j=1:nClutF
            fobj=fobj+0.1*(Dm[i]-xMean[i])*invCov[i,j]*(Dm[j]-xMean[j])
        end
    end
    if gn_relFlag==1 || gn_relFlag==2
        fobj=fobj+wpia*(piaKu-gn_srtPIA)^2
    end
    return fobj
end

function g_z(gradZ_out,Dm)
    piaKu,piaKa,rrate,zSim,gradZ=simZ_g(Dm)
    n=size(zSim)[1]
    fobj=0.0
    for i=1:nClutF
        if Z_Obs[i]>5 && Z_Obs[i]==Z_Obs[i]
            fobj=fobj+(zSim[i]-Z_Obs[i])^2
        end
    end
    gradZ_out=gradZ
    for i=1:nClutF
        if Z_Obs[i]<5 || Z_Obs[i]!=Z_Obs[i]
            gradZ_out[i,:].=0
            Dm[i]=0.2
        end
    end
    #gradZ_out.=gradZ_out .+ 0.2*invCov*(Dm-xMean)
    for i=1:nClutF
        for j=1:nClutF
            fobj=fobj+0.1*(Dm[i]-xMean[i])*invCov[i,j]*(Dm[j]-xMean[j])
        end
    end
    return gradZ,zSim,piaKu,rrate,piaKa
end


using LinearAlgebra

function gaussNewton(gradZ_out,xsol,i)
    rms_keep=1e8
    xsol_keep=xsol
    eye=Matrix{Float64}(I,88,88)
    rms=1e8
    for it=1:10
        rms=f_z(xsol)
        if rms>rms_keep
            xsol=xsol_keep
            break
        end
        #println("$(it) $i ",rms)
        gradZ,zSim,piaKu,rrate,piaKa=g_z(gradZ_out,xsol)
        #gradZ=transpose(gradZ);
        dY=zeros(89)
        dY[1:88]=Z_Obs-zSim
        #println(size(dZ))
        for k=1:nClutF
            if Z_Obs[k]<5
                dY[k]=0
            end
            if Z_Obs[k]!=Z_Obs[k]
                dY[k]=0
            end
        end
        if gn_relFlag==1 || gn_relFlag==2
            dY[89]=sqrt(wpia)*(gn_srtPIA-piaKu)
        end
        #println(size(gradZ))
        #println(size(dY))
        dy=transpose(gradZ)*(dY)-0.1*invCov*(xsol-xMean);
        A=transpose(gradZ)*gradZ+0.1*invCov;
        A=A+3*eye;
        #println(dy)
        dx=A\dy;
        #println(gradZ)
        #xit(1)
        xsol_keep=copy(xsol)
        rms_keep=rms
        xsol=xsol+0.5*dx;
        a1=findall(xsol.>3.75)
        xsol[a1].=3.75
        a1=findall(xsol.<0.19)
        xsol[a1].=0.19
        if rms<50
            break
        end
    end
    gradZ,zSimf,piaKuf,rretf,piaKaf=g_z(gradZ_out,xsol)
    return  gradZ,zSimf,piaKuf,rretf,piaKaf,rms, xsol
end

gradZ_out=zeros(88,88)





#    zWInt=zWInterp(h1[100:170])
#    radarMod.setZObs(copy(zInt))
    

#gradZ,zSimf,piaKuf,rretf,piaKaf,rms, xsol=gaussNewton(gradZ_out,xsol,i)
#piaKu,piaKa,piaW,rrateL,zSim,zSimKa,zSimW=simMZ(xsol)
