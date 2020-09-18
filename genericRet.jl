
using PyCall
using PyPlot
#bpsh -m to run on beorun cluster

fname="../subSets/2A-CS-CONUS.GPM.DPR.V8-20180723.20190402-S172036-E172914.028936.V06A.HDF5"
fCMB="2B-CS-CONUS.GPM.DPRGMI.CORRA2018.20190402-S104530-E105403.028932.V06A.HDF5"
fCMB="2B-CS-CONUS.GPM.DPRGMI.CORRA2018.20190402-S172036-E172914.028936.V06A.HDF5"
fCMB="2B-CS-CONUS.GPM.DPRGMI.CORRA2018.20190402-S172036-E172914.028936.V06A.HDF5"
fCMB="data/2B-CS-CONUS.GPM.DPRGMI.CORRA2018.20190531-S235253-E000132.029858.V06A.HDF5"
fname="../subSets/2A-CS-CONUS.GPM.DPR.V8-20180723.20190531-S235253-E000132.029858.V06A.HDF5"
fname="../subSets/2A.GPM.DPR.V8-20180723.20190531-S225830-E003103.029858.V06A.HDF5"
#node=f['NS/DSD/binNode'][n1:n2,nr1:nr2,:]
#flag=f['/NS/PRE/flagPrecip'][n1:n2,nr1:nr2]
#precipType=f['/NS/CSF/typePrecip'][n1:n2,nr1:nr2]
#sfcType=f['/NS/PRE/landSurfaceType'][n1:n2,nr1:nr2]
#zeroDeg=f['/NS/VER/binZeroDeg'][n1:n2,nr1:nr2]
#binClut=f['/NS/PRE/binClutterFreeBottom'][n1:n2,nr1:nr2]
#binSf=f['/NS/PRE/binRealSurface'][n1:n2,nr1:nr2]
#hZero=f['/NS/VER/heightZeroDeg'][n1:n2,nr1:nr2]

function iterativep(z13obs,fractS,nodes,dn1)
    piaKuH=0
    piaKuR=0
    dpiaKuRdn=0
    #dn1=0.0
    retRate=zeros(88)
    z13c=copy(z13obs)
    println(size(fractS))
    println(size(retRate))
    rHmax=0
    for it=1:5
        piaKuR=0
        dpiaKuRdn=0
        for k=nodes[1]+1:nodes[5]+1
            if z13obs[k]>0
                n1s,n2s=bisection(zKuTs[14:end],z13c[k])
                n1,n2=bisection(zKuT[14:end],z13c[k]-10*dn1)
                ddn=(zKuT[n1]-zKuT[n2+1])/10.0
                piaKuR+=(fractS[k]*10^dn1*attKuT[13+n1]+
                         (1-fractS[k])*attKuTs[13+n1s])*0.25
                #dpiaKuRdn+=2*((fractS[k]*10^(dn1+ddn)*attKuT[13+n2+1]+
                #               (1-fractS[k])*attKuTs[13+n1s])*0.25-
                #              (fractS[k]*10^dn1*attKuT[13+n1]+
                #               (1-fractS[k])*attKuTs[13+n1s])*0.25)/ddn
                z13c[k]=z13obs[k]+piaKuR

                piaKuR+=(fractS[k]*10^dn1*attKuT[13+n1]+
                         (1-fractS[k])*attKuTs[13+n1s])*0.25
                retRate[k]=(fractS[k]*10^dn1*rateT[13+n1]+
                            (1-fractS[k])*rateTs[13+n1s])
                if fractS[k]<0.1 && rHmax<retRate[k]
                    rHmax=retRate[k]
                end
            end
        end
    end
    return piaKuR,z13c,retRate,dpiaKuRdn,rHmax
end
include("src_jl/structDef.jl")
include("src_jl/scatTables.jl")
include("src_jl/radarMod.jl")
include("src_jl/radarOE.jl")
include("src_jl/psdInt_2.jl")
include("src_jl/readInput.jl")
using .radarMod
using LinearAlgebra

radarMod.setRadar_dr(0.25,88)
fract=zeros(88)
radarMod.setRainFract(fract)
radarMod.setdNw(zeros(88).+0.0)
invCov=Matrix{Float64}(I,88,88)

nt=472
nst=0


fnameT="tables_nsv_rho025_2_dmTs11.nc"

using .scatTables
zKuT,zKaT,attKuT,attKaT,dmT,nwT,pwcT,zKuTs,zKaTs,attKuTs,attKaTs,dmTs,
dmTs,nwTs,pwcTs,attWTs,attWT,zWTs,zWT,rateT,rateTs=readTables(fnameT)


np=pyimport("numpy")
dat=[[19,15],
     [17,25],
     [14,33],
     [6,41],
     [5,45],
     [4.75,44],
     [2,30]]

dat=np["array"](dat)

zObs=np["interp"](0:0.25:18.5,dat[end:-1:1,1],dat[end:-1:1,2])

rrate1d=zeros(75)


dr=0.25
piaKu=0.0
hg=0:0.25:18.5

datR=[[19,0.5],
      [17,1.5],
      [14,6.5],
      [6,30],
      [4.25,40],
      [0,44]]
datR=np["array"](datR)
rDat=datR[:,2].*2.5

function simZ(rDat,dat,hg,datR)
    zObs=np["interp"](0:0.25:18.5,dat[end:-1:1,1],dat[end:-1:1,2])
    rInt=np["interp"](0:0.25:18.5,datR[end:-1:1,1],rDat[end:-1:1])
    fnoise=rand(75).+0.5
    fnoiss=rand(75).+0.5
    for k=2:74
        fnoiss[k]=0.7*fnoise[k]+0.15*fnoise[k+1]+0.15*fnoise[k-1]
    end
    rInt=rInt.*fnoise
    piaKu=0.0
    dn=0.5
    zKu=copy(zObs).*0
    for k=75:-1:29
        n1,n2=bisection(rateTs,rInt[k]/10^dn)
        #global piaKu
        #println("$(hg[k]) $(n1) $(piaKu) $(zKuTs[n1])")
        piaKu+=attKuTs[n1]*dr*10^dn
        zKu[k]=zKuTs[n1]-piaKu+10*dn
        piaKu+=attKuTs[n1]*dr*10^dn
    end
    for k=28:-1:19
        f=(k-19)/9.0
        #global piaKu
        n1s,n2s=bisection(rateTs,rInt[k]/10^dn)
        n1,n2=bisection(rateT,rInt[k]/10^dn)
        piaKu+=f*10^dn*attKuTs[n1s]*dr+10^dn*(1-f)*attKuT[n1]*dr
        zKu[k]=10*log10(f*10^(0.1*zKuTs[n1s]+dn)+(1-f)*10^(0.1*zKuT[n1]+dn))-piaKu
        piaKu+=10^dn*(f*attKuTs[n1s]*dr+(1-f)*attKuT[n1]*dr)
    end
    for k=18:-1:1
        #global piaKu
        n1,n2=bisection(rateT,rInt[k]/10^dn)
        piaKu+=attKuT[n1]*dr*10^dn
        zKu[k]=zKuT[n1]-piaKu+10*dn
        piaKu+=attKuT[n1]*dr*10^dn
    end
    return zKu,piaKu,rInt
end
cfadZ=zeros(75,50)
zKuL=[]
rrateL=[]
@time for f in 0.125:0.000125:1.5
    rDat=datR[:,2].*f
    f1=1.5*rand(6).+0.25
    rDat=rDat.*f1
    zKu,piaKu,rrate=simZ(rDat,dat,hg,datR)
    for k=1:75
        if zKu[k]>0 && zKu[k]<49
            ik=Int(trunc(zKu[k]))+1
            if ik>0 && ik<=50
                cfadZ[k,ik]+=1
            end
        end
    end
    push!(zKuL,zKu)
    push!(rrateL,rrate)
    #plot(zKu,hg)
end

fh=pybuiltin("open")("paramZObs.pklz","wb")
pklz=pyimport("pickle")
zKuL=np["array"](zKuL);
rrateL=np["array"](rrateL);
pklz["dump"]([zKuL,rrateL],fh)
fh["close"]()


#for it=1:3
#    piaKu=0.10
#    for k=65:-1:17
#        n1,n2=bisection(zKuTs,zKu[k]-3)
#        piaKu+=10^0.3*attKuTs[n1]*dr
#        zKu[k]=zObs[k]+piaKu
#        piaKu+=10^0.3*attKuTs[n1]*dr
#        rrate1d[k]=10^0.3*rateTs[n1]
#    end
#    n1,n2=bisection(zKuT,zObs[13]+piaKu-6)
#    rrate1d[13]=10^0.6*rateT[n1]
#
#    println(piaKu," ",10^0.6*attKuT[n1])
#end
