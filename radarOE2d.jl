using PyCall

pickle=pyimport("pickle")

push!(PyVector(pyimport("sys")["path"]), "./")

f=pybuiltin("open")("IMPACTS0207.pklz","rb")
dat=pickle["load"](f)

zKuL,zKaL,ibbL,hZ,r,zWL,xyL,tL=dat

println(size(zKuL))

include("scatTables.jl")
using .scatTables

temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r,
tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r,
tempW,massW,fractionW,bscatW,DeqW,extW,scatW,gW,vfallW,tempW_r,massW_r,bscatW_r,DeqW_r,
extW_r,scatW_r,gW_r,vfallW_r=scatTables.init()

xr=pyimport("xarray")
function saveTables()
    nwTX=xr["DataArray"](nwT)
    pwcTX=xr["DataArray"](pwcT)
    attKaTX=xr["DataArray"](attKaT)
    attKuTX=xr["DataArray"](attKuT)
    zKaTX=xr["DataArray"](zKaT)
    zKuTX=xr["DataArray"](zKuT)
    zKaTsX=xr["DataArray"](zKaTs)
    zKuTsX=xr["DataArray"](zKuTs)
    dmTsX=xr["DataArray"](dmTs)
    dmTX=xr["DataArray"](dmT)
    pwcTsX=xr["DataArray"](pwcTs)
    attKaTsX=xr["DataArray"](attKaTs)
    attKuTsX=xr["DataArray"](attKuTs)
    attWTsX=xr["DataArray"](attWTs)
    attWTX=xr["DataArray"](attWT)
    zWTsX=xr["DataArray"](zWTs)
    zWTX=xr["DataArray"](zWT)
    nwTsX=xr["DataArray"](nwTs)
    rateX=xr["DataArray"](rateT)
    ratesX=xr["DataArray"](rateTs)
    d=Dict("nwT"=>nwTX,"pwcT"=>pwcTX,"attKaT"=>attKaTX,
           "attKuT"=>attKuTX, "dmT"=>dmTX, "zKuT"=>zKuTX, "zKaT"=>zKaTX,
           "nwTs"=>nwTsX,"pwcTs"=>pwcTsX,"attKaTs"=>attKaTsX,
           "attKuTs"=>attKuTsX, "dmTs"=>dmTsX, "zKuTs"=>zKuTsX, "zKaTs"=>zKaTsX,
           "zWTs"=>zWTsX,"attWTs"=>attWTsX,"zWT"=>zWTX,"attWT"=>attWTX, "rate"=>rateX,
           "rateS"=>ratesX)
    tables=xr["Dataset"](d)
    tables["to_netcdf"]("tables_nsv_rho005_dmTs11.nc")
end

iread=0
if iread==1
    ns=9
    scatTables.getDmNwSF(tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
    tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r,
    temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
    temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r,ns,
    tempW,massW,fractionW,bscatW,DeqW,extW,scatW,gW,vfallW,
    tempW_r,massW_r,bscatW_r,DeqW_r,extW_r,scatW_r,gW_r,vfallW_r)
    scatTables.getDmNwR(tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
    tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r,
    temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
    temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r,
    tempW,massW,fractionW,bscatW,DeqW,extW,scatW,gW,vfallW,
    tempW_r,massW_r,bscatW_r,DeqW_r,extW_r,scatW_r,gW_r,vfallW_r)
    saveTables()
    exit(1)
end  
n4=pyimport("netCDF4")
fh=n4["Dataset"]("tables_nsv_rho005_dmTs11.nc")
zKuT=fh["variables"]["get"]("zKuT")["_get"]([0],[240],[1])
zKaT=fh["variables"]["get"]("zKaT")["_get"]([0],[240],[1])
attKuT=fh["variables"]["get"]("attKuT")["_get"]([0],[240],[1])
attKaT=fh["variables"]["get"]("attKaT")["_get"]([0],[240],[1])
dmT=fh["variables"]["get"]("dmT")["_get"]([0],[240],[1])
nwT=fh["variables"]["get"]("nwT")["_get"]([0],[240],[1])
pwcT=fh["variables"]["get"]("pwcT")["_get"]([0],[240],[1])
zKuTs=fh["variables"]["get"]("zKuTs")["_get"]([0],[240],[1])
zKaTs=fh["variables"]["get"]("zKaTs")["_get"]([0],[240],[1])
attKuTs=fh["variables"]["get"]("attKuTs")["_get"]([0],[240],[1])
attKaTs=fh["variables"]["get"]("attKaTs")["_get"]([0],[240],[1])
dmTs=fh["variables"]["get"]("dmTs")["_get"]([0],[240],[1])
nwTs=fh["variables"]["get"]("nwTs")["_get"]([0],[240],[1])
pwcTs=fh["variables"]["get"]("pwcTs")["_get"]([0],[240],[1])
attWTs=fh["variables"]["get"]("attWTs")["_get"]([0],[240],[1])
attWT=fh["variables"]["get"]("attWT")["_get"]([0],[240],[1])
zWTs=fh["variables"]["get"]("zWTs")["_get"]([0],[240],[1])
zWT=fh["variables"]["get"]("zWT")["_get"]([0],[240],[1])
rateTs=fh["variables"]["get"]("rateS")["_get"]([0],[240],[1])
rateT=fh["variables"]["get"]("rate")["_get"]([0],[240],[1])
include("attCorrection.jl")

push!(PyVector(pyimport("sys")["path"]), "./")

np=pyimport("numpy")
h=0.125*((176:-1:100).-100);
include("radarMod.jl")
fractS=Spline1D([0.,1.25,4.5,4.995,7.25,10.],[0,0,0,0.0,0,0],k=1)
fractD=Spline1D([0.,1.25,3.85,6.25,10.],[1.8,1.8,1.7,1.6,0.7])

h1=(176:-1:1)*0.125
fract=fractS(h1[100:170])
dmMV=fractD(h1[100:170])

using PyPlot
using LinearAlgebra



scipy=pyimport("scipy.ndimage") 

#include("radarMod.jl")
using .radarMod
using LinearAlgebra
radarMod.setRadar_dr(0.125,71)
radarMod.setRainFract(fract)
radarMod.setdNw(zeros(71).+0.0)
invCov=Matrix{Float64}(I,71,71)
#
radarMod.setInvCov(invCov,dmMV)
hZ1=(hZ.-r)./1000


function simZ(Dm)
    drk=Δr
    fract=rainFract
    nz=size(Dm)[1]
    piaKu=0
    piaKa=0.0
    zSim=zeros(nz)
    rrateL=zeros(nz)
    noise=dNw
    for i=1:nz
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
    return piaKu,rrateL,zSim
end

function simMZ(Dm)
    drk=Δr
    fract=rainFract
    nz=size(Dm)[1]
    piaKu=0
    piaKa=0.0
    piaW=0.0
    zSim=zeros(nz)
    zSimKa=zeros(nz)
    zSimW=zeros(nz)
    rrateL=zeros(nz)
    noise=dNw
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
        piaKa=piaKa+attKa*drk
        attW=(1-fract[i])*attWTs[n1s]+fract[i]*attWT[n1]
        attW=attW*10^noise[i]
        piaW=piaW+attW*drk
        zSim[i]=log10((1-fract[i])*10^(0.1*zKuTs[n1s])+(fract[i])*10^(0.1*zKuT[n1]))*10.0-piaKu+10*noise[i]
        zSimKa[i]=log10((1-fract[i])*10^(0.1*zKaTs[n1s])+(fract[i])*10^(0.1*zKaT[n1]))*10.0-piaKa+10*noise[i]
        zSimW[i]=log10((1-fract[i])*10^(0.1*zWTs[n1s])+(fract[i])*10^(0.1*zWT[n1]))*10.0-piaW+10*noise[i]

        rrateL[i]=10^noise[i]*((1-fract[i])*rateTs[n1s]+fract[i]*rateT[n1])
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
    zSim=zeros(nz)
    rrateL=zeros(nz)
    noise=dNw
    dZdDm=zeros(nz,nz)
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
            dZdDm[i,i+1:nClutF]=dZdDm[i,i+1:nClutF].-2*dpiaKu1
        end
    end
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
    gradZ_out.=2*gradZ*(zSim-Z_Obs)
    for i=1:nClutF
        if Z_Obs[i]<5 || Z_Obs[i]!=Z_Obs[i]
            gradZ_out[i,:].=0
            Dm[i]=0.2
        end
    end
    gradZ_out.=gradZ_out .+ 0.2*invCov*(Dm-xMean)
    for i=1:nClutF
        for j=1:nCluF
            fobj=fobj+0.1*(Dm[i]-xMean[i])*invCov[i,j]*(Dm[j]-xMean[j])
        end
    end
    return gradZ,zSim,piaKu,rrate,piaKa
    #return fobj
end

#xsol=dmMV

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
        println("$(it) $i ",rms)
        gradZ,zSim,piaKu,rrate,piaKa=g_z(gradZ_out,xsol)
        #gradZ=transpose(gradZ);
        dy=transpose(gradZ)*(Z_Obs-zSim)-invCov*(xsol-xMean);
        A=transpose(gradZ)*gradZ+invCov;
        A=A+3*eye;
        dx=A\dy;
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

gradZ_out=zeros(71,71)

a=findall(zKuL.!=zKuL)
zKuL[a].=-5;
a=findall(zKaL.!=zKaL)
zKaL[a].=-5;
a=findall(zWL.!=zWL)
zWL[a].=-5;

h1d=h1[100:170]




plt=pyimport("matplotlib.pyplot")
zWObs=zeros(200,71)
zWSim=zeros(200,71)
dZ=zeros(71)
dmT=zeros(200,71)
dmTSF=zeros(200,71)
for i=1:200
    dn1=zeros(71).+0.5
    radarMod.setdNw(dn1)
    global dZ
    #global dn1
    zInterp=Spline1D(hZ1[end:-1:160],zKuL[100+i,end:-1:160])
    zKaInterp=Spline1D(hZ1[end:-1:165],zKaL[100+i,end:-1:165])
    zWInterp=Spline1D(hZ1[end:-1:160],zWL[100+i,end:-1:160])
    zInt=zInterp(h1[100:170])
    zKaInt=zKaInterp(h1[100:170])
    zWInt=zWInterp(h1[100:170])
    radarMod.setZObs(copy(zInt))
    
    xsol=dmMV
    gradZ,zSimf,piaKuf,rretf,piaKaf,rms, xsol=gaussNewton(gradZ_out,xsol,i)
    piaKu,piaKa,piaW,rrateL,zSim,zSimKa,zSimW=simMZ(xsol)
    zWObs[i,:]=zWInt
    zWSim[i,:]=zSimW
    dZ=zWObs[i,:]-zWSim[i,:]
    dmTSF[i,:]=xsol
    for k=1:71
        if zWObs[i,k]>0 && dZ[k]>0
            dn1[k]+=0.125*dZ[k]
        end
    end
    radarMod.setdNw(dn1)
    gradZ,zSimf,piaKuf,rretf,piaKaf,rms, xsol=gaussNewton(gradZ_out,xsol,i)
    piaKu,piaKa,piaW,rrateL,zSim,zSimKa,zSimW=simMZ(xsol)
    zWObs[i,:]=zWInt
    zWSim[i,:]=zSimW
    dmT[i,:]=xsol
end
#
h1d=h1[100:170]
plotZ=pyimport("plotZRet")
plotZ["plotret"](zWSim,zWObs,h1d,dmT,dmTSF,tL[101:300])
