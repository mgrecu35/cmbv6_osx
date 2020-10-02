
#pickle=pyimport("pickle")
#fh=pybuiltin("open")("IMPACTS0201.pklz","rb")
#impactsData=pickle["load"](fh)
#pyimport("sys")




module scatTables
using PyCall
include("psdInt.jl")
include("psdInt_2.jl")

pushfirst!(PyVector(pyimport("sys")["path"]), ".")
export zKuT,dmT,nwT,pwcT,attKaT,attKuT,zKaT
export zKaTs,zKuTs,dmTs,pwcTs,attKuTs,attKaTs,nwTs
export zWTs,attWTs, zWT, attWT, rateT,rateTs
nTs=240
nT=240
nwT=zeros(nT)
pwcT=zeros(nT)
attKaT=zeros(nT)
attKuT=zeros(nT)
zKaT=zeros(nT)
zKuT=zeros(nT)
rateT=zeros(nT)
zWT=zeros(nT)
attWT=zeros(nT)
zWTs=zeros(nTs)
attWTs=zeros(nTs)
rateTs=zeros(nTs)
zKaTs=zeros(nTs)
zKuTs=zeros(nTs)
dmTs=zeros(nTs)
pwcTs=zeros(nTs)
attKaTs=zeros(nTs)
attKuTs=zeros(nTs)
nwTs=zeros(nTs)
for i=1:nT
    zKuT[i]=0+(i-1)*0.25
end
for i=1:nTs
    zKuTs[i]=-7+(i-1)*0.25
end
zCoeffs=[ 0.00159248, -0.01061376,  0.33276274]
dmT=zCoeffs[1].*zKuT.*zKuT+zCoeffs[2].*zKuT.+zCoeffs[3]
dmTs.=dmT
ds=maximum(dmT)-minimum(dmT)
for i=1:nT
    dmT[i]=0.75*(dmT[1]+0.75*(dmT[i]-dmT[1]))*1.3
    dmTs[i]=0.75*(dmTs[1]+0.75*(dmTs[i]-dmTs[1]))*1.3
end
for i=1:14
    dmT[i]=dmT[i]+(i-14)*0.01
    dmTs[i]=dmTs[i]+(i-14)*0.01
end
dmTs=dmTs.*1.2
println(dmT[1:14])
println(dmTs[1:14])
#exit(1)
#
function init()
    #global readS
    readS=pyimport("readScattProf")
    fnameIce="ice-self-similar-aggregates_35-GHz_scat.nc"
    fnameIceKu="ice-self-similar-aggregates_13-GHz_scat.nc"
    fnameRain="liquid-water_35-GHz_scat.nc"
    fnameRainKu="liquid-water_13-GHz_scat.nc"
    fnameIceW="ice-self-similar-aggregates_94-GHz_scat.nc"
    fnameRainW="liquid-water_94-GHz_scat.nc"

    fnameIceX="ice-self-similar-aggregates_X_scat.nc"
    fnameRainX="liquid-water_X_scat.nc"

    temp,mass,fraction,bscat,Deq,
    ext,scat,g,vfall=readS["readScatProf"](fnameIce);
    temp_r,mass_r,bscat_r,Deq_r,
    ext_r,scat_r,g_r,vfall_r=readS["readScatProfR"](fnameRain);
    tempKu,massKu,fractionKu,bscatKu,DeqKu,
    extKu,scatKu,gKu,vfallKu=readS["readScatProf"](fnameIceKu);
    tempKu_r,massKu_r,bscatKu_r,DeqKu_r,
    extKu_r,scatKu_r,gKu_r,vfallKu_r=readS["readScatProfR"](fnameRainKu);
    tempW,massW,fractionW,bscatW,DeqW,
    extW,scatW,gW,vfallW=readS["readScatProf"](fnameIceW);
    tempW_r,massW_r,bscatW_r,DeqW_r,
    extW_r,scatW_r,gW_r,vfallW_r=readS["readScatProfR"](fnameRainW);
    tempX,massX,fractionX,bscatX,DeqX,
    extX,scatX,gX,vfallX=readS["readScatProf"](fnameIceX);
    tempX_r,massX_r,bscatX_r,DeqX_r,
    extX_r,scatX_r,gX_r,vfallX_r=readS["readScatProfR"](fnameRainX);
    return temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
    temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r,
    tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
    tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r,
    tempW,massW,fractionW,bscatW,DeqW,extW,scatW,gW,vfallW,tempW_r,massW_r,bscatW_r,DeqW_r,
    extW_r,scatW_r,gW_r,vfallW_r,
    tempX,massX,fractionX,bscatX,DeqX,extX,scatX,gX,vfallX,tempX_r,massX_r,bscatX_r,DeqX_r,
    extX_r,scatX_r,gX_r,vfallX_r
    #exit()
end
export getDmNwR, getDmNwS
function getDmNwR(tempKu::Array{Float32,1},massKu::Array{Float32,2},fractionKu::Array{Float32,1},
    bscatKu::Array{Float32,3},DeqKu::Array{Float32,2},extKu::Array{Float32,3},scatKu::Array{Float32,3},
    gKu::Array{Float32,3},vfallKu::Array{Float32,2},
    tempKu_r::Array{Float32,1},massKu_r::Array{Float32,1},bscatKu_r::Array{Float32,2},
    DeqKu_r::Array{Float32,1},extKu_r::Array{Float32,2},scatKu_r::Array{Float32,2},
    gKu_r::Array{Float32,2},vfallKu_r::Array{Float32,1},
    temp::Array{Float32,1},mass::Array{Float32,2},fraction::Array{Float32,1},
    bscat::Array{Float32,3},Deq::Array{Float32,2},ext::Array{Float32,3},scat::Array{Float32,3},
    g::Array{Float32,3},vfall::Array{Float32,2},
    temp_r::Array{Float32,1},mass_r::Array{Float32,1},bscat_r::Array{Float32,2},
    Deq_r::Array{Float32,1},ext_r::Array{Float32,2},scat_r::Array{Float32,2},
    g_r::Array{Float32,2},vfall_r::Array{Float32,1},
    tempW::Array{Float32,1},massW::Array{Float32,2},fractionW::Array{Float32,1},
    bscatW::Array{Float32,3},DeqW::Array{Float32,2},extW::Array{Float32,3},scatW::Array{Float32,3},
    gW::Array{Float32,3},vfallW::Array{Float32,2},
    tempW_r::Array{Float32,1},massW_r::Array{Float32,1},bscatW_r::Array{Float32,2},
    DeqW_r::Array{Float32,1},extW_r::Array{Float32,2},scatW_r::Array{Float32,2},
    gW_r::Array{Float32,2},vfallW_r::Array{Float32,1})
    nx=size(zKuT)[1]
    println(nx)
    mu=1.0
    dn=0.1
    freqKu=13.8
    wlKu=300/freqKu
    freqKa=35.5
    wlKa=300/freqKa
    zCoeffs=[ 0.00159248, -0.01061376,  0.33276274]
    freqW=94.0
    wlW=300/freqW
    for i=1:nx
        zObs1=zKuT[i]
        dnRet=dn
        for it=1:4
            dn=dnRet
            global dmc
            retZ=get_fZr(zObs1,bscatKu_r,scatKu_r,extKu_r,gKu_r,DeqKu_r,vfallKu_r,wlKu,dn,mu)
            rwc, Z, att,scatInt,gInt, vdop, dm, rrate=retZ
            dmc=0.1*dmT[i]
            dn1=dn*10
            retZ1=get_fZr(zObs1,bscatKu_r,scatKu_r,extKu_r,gKu_r,DeqKu_r,vfallKu_r,wlKu,dn1,mu)
            rwc1, Z1, att1,scatInt1,gInt1, vdop1, dm1, rrate1=retZ1
            gradDm=(dm1-dm)/(log(dn1/dn))
            dnRet=dn*exp(1.0*(dmc-dm)*gradDm/(gradDm*gradDm+0.000001))
        end
        #dnRet=0.1
        retZ=get_fZr(zObs1,bscatKu_r,scatKu_r,extKu_r,gKu_r,DeqKu_r,vfallKu_r,wlKu,dnRet,mu)
        rwc, Z, att,scatInt,gInt, vdop, dm, rrate=retZ
        pwcT[i]=rwc
        attKuT[i]=att*4.343
        nwT[i]=log10(dnRet)
        zKa, attKa,scatIntKa,gIntKa, vdopKa,dmKa,rrate=get_Zr(rwc,bscat_r,scat_r,ext_r,g_r,Deq_r,vfall_r,wlKa,dnRet,mu)
        #println("$Z $(zKu[i]) $(dm) $(dmc)")
        zKaT[i]=zKa
        attKaT[i]=attKa*4.343
        zW, attW,scatIntW,gIntW, vdopW,dmW,rrate=get_Zr(rwc,bscatW_r,scatW_r,extW_r,gW_r,DeqW_r,vfallW_r,wlW,dnRet,mu)
        zWT[i]=zW
        attWT[i]=attW*4.343
        rateT[i]=rrate
        #println("$(zKuT[i]) $(zKaT[i]) $(zWT[i]) $(dm) $(wlW) $(wlKa) $(wlKu)")
        #println("$(zKuT[i]) $(zKaT[i]) $(zWT[i]) $(dm) $(dmc)")
#get_fZr(zKu[i],bscatKu_r,scatKu_r,extKu_r,gKu_r,DeqKu_r,vfallKu_r,wlKu,dn,mu)
    end
    #exit(1)
end
function getDmNwS(tempKu::Array{Float32,1},massKu::Array{Float32,2},fractionKu::Array{Float32,1},
    bscatKu::Array{Float32,3},DeqKu::Array{Float32,2},extKu::Array{Float32,3},scatKu::Array{Float32,3},
    gKu::Array{Float32,3},vfallKu::Array{Float32,2},
    tempKu_r::Array{Float32,1},massKu_r::Array{Float32,1},bscatKu_r::Array{Float32,2},
    DeqKu_r::Array{Float32,1},extKu_r::Array{Float32,2},scatKu_r::Array{Float32,2},
    gKu_r::Array{Float32,2},vfallKu_r::Array{Float32,1},
    temp::Array{Float32,1},mass::Array{Float32,2},fraction::Array{Float32,1},
    bscat::Array{Float32,3},Deq::Array{Float32,2},ext::Array{Float32,3},scat::Array{Float32,3},
    g::Array{Float32,3},vfall::Array{Float32,2},
    temp_r::Array{Float32,1},mass_r::Array{Float32,1},bscat_r::Array{Float32,2},
    Deq_r::Array{Float32,1},ext_r::Array{Float32,2},scat_r::Array{Float32,2},
    g_r::Array{Float32,2},vfall_r::Array{Float32,1})
    nx=size(zKaTs)[1]
    println(nx)
    mu=-1.0
    dn=0.1
    freqKu=13.8
    wlKu=300/freqKu
    freqKa=35.5
    wlKa=300/freqKa
    ns=11
    for i=1:nx
        zObs1=zKaTs[i]
        dnRet=dn
        rwc=0.1
        for it=1:4
            dn=dnRet
            global dmc
            retZ=get_fZs(zObs1,ns,bscat,scat,ext,g,Deq,vfall,wlKa,dn,mu)
            rwc, Z, att,scatInt,gInt, vdop, dm=retZ
            dmc=0.1*dmTs[i]
            dn1=dn*10
            retZ1=get_fZs(zObs1,ns,bscat,scat,ext,g,Deq,vfall,wlKa,dn1,mu)
            rwc1, Z1, att1,scatInt1,gInt1, vdop1, dm1=retZ1
            gradDm=(dm1-dm)/(log(dn1/dn))
            dnRet=dn*exp(1.0*(dmc-dm)*gradDm/(gradDm*gradDm+0.000001))
        end
        retZ=get_fZs(zObs1,ns,bscat,scat,ext,g,Deq,vfall,wlKa,dnRet,mu)
        rwc, Z, att,scatInt,gInt, vdop, dm=retZ
        pwcTs[i]=rwc
        attKaTs[i]=att*4.343
        nwTs[i]=log10(dnRet)
        zKaTs[i]=Z
        zKu, attKu,scatIntKu,gIntKu, vdopKu,dmKu=
        get_Zs(rwc,ns,bscatKu,scatKu,extKu,gKu,DeqKu,vfallKu,wlKu,dnRet,mu)
        #println("$Z $(zKu[i]) $(dm) $(dmc)")

        zKuTs[i]=zKu
        attKuTs[i]=attKu*4.343
    end
end
function getDmNwSF(tempKu::Array{Float32,1},massKu::Array{Float32,2},fractionKu::Array{Float32,1},
    bscatKu::Array{Float32,3},DeqKu::Array{Float32,2},extKu::Array{Float32,3},scatKu::Array{Float32,3},
    gKu::Array{Float32,3},vfallKu::Array{Float32,2},
    tempKu_r::Array{Float32,1},massKu_r::Array{Float32,1},bscatKu_r::Array{Float32,2},
    DeqKu_r::Array{Float32,1},extKu_r::Array{Float32,2},scatKu_r::Array{Float32,2},
    gKu_r::Array{Float32,2},vfallKu_r::Array{Float32,1},
    temp::Array{Float32,1},mass::Array{Float32,2},fraction::Array{Float32,1},
    bscat::Array{Float32,3},Deq::Array{Float32,2},ext::Array{Float32,3},scat::Array{Float32,3},
    g::Array{Float32,3},vfall::Array{Float32,2},
    temp_r::Array{Float32,1},mass_r::Array{Float32,1},bscat_r::Array{Float32,2},
    Deq_r::Array{Float32,1},ext_r::Array{Float32,2},scat_r::Array{Float32,2},
    g_r::Array{Float32,2},vfall_r::Array{Float32,1},ns,
    tempW::Array{Float32,1},massW::Array{Float32,2},fractionW::Array{Float32,1},
    bscatW::Array{Float32,3},DeqW::Array{Float32,2},extW::Array{Float32,3},scatW::Array{Float32,3},
    gW::Array{Float32,3},vfallW::Array{Float32,2},
    tempW_r::Array{Float32,1},massW_r::Array{Float32,1},bscatW_r::Array{Float32,2},
    DeqW_r::Array{Float32,1},extW_r::Array{Float32,2},scatW_r::Array{Float32,2},
    gW_r::Array{Float32,2},vfallW_r::Array{Float32,1})
    nx=size(zKaTs)[1]
    println(nx)
    mu=2.0
    dn=0.1
    freqKu=13.8
    wlKu=300/freqKu
    freqKa=35.5
    wlKa=300/freqKa
    freqW=94.0
    wlW=300/freqW
    for i=1:nx
        zObs1=zKuTs[i]
        dnRet=dn
        rwc=0.1
        rho=0.0093*(0.1*dmTs[i])^(-0.90)
        if rho<0.85
            rho=0.85
        end
        if rho>1
            rho=1
        end
        ns1,ns2=bisection(fractionKu,rho)
        #println("$(ns1) $(rho)")
        for it=1:4
            dn=dnRet
            global dmc
            retZ=get_fZs(zObs1,ns1,bscatKu,scatKu,extKu,gKu,DeqKu,DeqKu_r,vfallKu,vfallKu_r,wlKu,dn,mu)
            rwc, Z, att,scatInt,gInt, vdop, dm, rrate=retZ
            dmc=0.1*dmTs[i]
            dn1=dn*10
            retZ1=get_fZs(zObs1,ns1,bscatKu,scatKu,extKu,gKu,DeqKu,DeqKu_r,vfallKu,vfallKu_r,wlKu,dn1,mu)
            rwc1, Z1, att1,scatInt1,gInt1, vdop1, dm1, rrate1=retZ1
            gradDm=(dm1-dm)/(log(dn1/dn))
            dnRet=dn*exp(1.0*(dmc-dm)*gradDm/(gradDm*gradDm+0.000001))

        end
        println("dmfirst=$(dm) $(rwc) $(dmc)")
        #dnRet=0.1
        retZ=get_fZs(zObs1,ns1,bscatKu,scatKu,extKu,gKu,DeqKu,DeqKu_r,vfallKu,vfallKu_r,wlKu,dnRet,mu)
        rwc, Z, att,scatInt,gInt, vdop, dm,rrate=retZ
        pwcTs[i]=rwc
        attKuTs[i]=att*4.343
        nwTs[i]=log10(dnRet)
        zKuTs[i]=Z
        zKa, attKa,scatIntKa,gIntKa, vdopKa,dmKa,rrate=
        get_rZs(rwc,dmc,ns1,bscat,scat,ext,g,Deq,Deq_r,vfall,vfall_r,wlKa,dnRet,mu)
        zKu1, attKu1,scatIntKu1,gIntKu1, vdopKu1,dmKu1,rrate1=
        get_rZs(rwc,dmc,ns1,bscatKu,scatKu,extKu,gKu,DeqKu,DeqKu_r,vfallKu,vfallKu_r,wlKu,dnRet,mu)
        println("zsims=$(zKu1) $(Z) $(zObs1)")
        zW, attW,scatIntW,gIntW, vdopW,dmW,rrate=
        get_rZs(rwc,dmc,ns1,bscatW,scatW,extW,gW,DeqW,DeqW_r,vfallW,vfallW_r,wlW,dnRet,mu)
        rateTs[i]=rrate
        #println("$Z $(zKu[i]) $(dm) $(dmc)")
        zKaTs[i]=zKa
        attKaTs[i]=attKa*4.343
        zWTs[i]=zW
        attWTs[i]=attW*4.343
        #println("$(zKuTs[i]) $(zKaTs[i]) $(zWTs[i]) $(dmc) $(dmT[i]) ")
        #exit(1)
    end
    #exit(1)
end
end
