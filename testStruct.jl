
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


include("src_jl/structDef.jl")
#include("src_jl/scatTables.jl")
#include("src_jl/radarMod.jl")
#include("src_jl/radarOE.jl")
#include("src_jl/psdInt_2.jl")
include("src_jl/readInput.jl")
#using .radarMod
using LinearAlgebra

#radarMod.setRadar_dr(0.25,88)
fract=zeros(88)
#radarMod.setRainFract(fract)
#radarMod.setdNw(zeros(88).+0.0)
#invCov=Matrix{Float64}(I,88,88)

nt=472
nst=0
sfcPrecip,cmbLat,cmbLon,sfcPrecipMS,piaSRT,relFlag=readcomb(fCMB,nst,nt)
nst=4660

zKu,zKa,pType,binSf,hZero,dprNodes,clutFree,reliabF,piaKu,lat,
lon,piaKa,reliabKaF,binZeroDeg,sfcType=readnetcdf(fname,nst,nt)


#fnameT="tables_nsv_rho085_2_dmTs11.nc"
#makeTables(fnameT)
#exit(1)
fnameT="tables_nsv_rho085_2_dmTs11.nc"

#using .scatTables
#zKuT,zKaT,attKuT,attKaT,dmT,nwT,pwcT,zKuTs,zKaTs,attKuTs,attKaTs,dmTs,
#dmTs,nwTs,pwcTs,attWTs,attWT,zWTs,zWT,rateT,rateTs=readTables(fnameT)

#println(minimum(dmTs))
#println(minimum(dmT))
#println(maximum(rateT))


ccall((:readtablesliang2_,"./combAlg"),Cvoid,(Ref{Int32},Ref{Int32}),
      nmu,nmfreq)
ccall((:cloud_init_,"./combAlg"),Cvoid,(Ref{Cint},),nmfreq)
ccall((:__nbinmod_MOD_init_nbin,"./combAlg"),Cvoid,(),)
sysdN=Cfloat(-0.25) # if stratiform -=0.1
ccall((:initgeophys_,"./combAlg"),Cvoid,(Ref{Cint},Ref{Cint},Ref{Cfloat}),
      nmemb,nmfreq,sysdN)


nt,nr=size(clutFree)
println(size(dprNodes))
piaRetL=[]
rrate2D=zeros(nt,49).-99.9
r1L=[]
r2L=[]
r3L=[]
using Statistics
np=pyimport("numpy")
rrateGNL=[]
rrateHBL=[]
dmL=[]
indL=[]
ind=[[343,14],[343,15],[343,16],[344,13],[344,14],[344,15]]
for i1=1:1 #i=342:345#1:nt
    #for j=13:37
    for ind1 in ind[1:1]
        i,j=ind1
        stormType.rainType=2;
        if Int32(trunc(pType[i,j]/1e7))==2 &&
            dprNodes[i,j,4]<dprNodes[i,j,5] && dprNodes[i,j,1]<binZeroDeg[i,j]-4 &&
            binZeroDeg[i,j]<150 && zKu[i,j,binZeroDeg[i,j]]>10
            nodes[1]=Int32(trunc(dprNodes[i,j,1]/2))
            nodes[3]=Int32(trunc(dprNodes[i,j,3]/2))
            nodes[3]=Int32(trunc((binZeroDeg[i,j])/2))
            nodes[2]=nodes[3]-3
            nodes[4]=nodes[3]+2
            nodes[5]=Int32(trunc(clutFree[i,j]/2))
            if nodes[4]>=nodes[5]
                continue
            end
            if nodes[1]>nodes[2]-3
                nodes[1]=nodes[2]-3
            end
            println(nodes)
            nodesF=copy(nodes).+1
            nodesF[1]=1
            nodesF[5]=88
            #fractInt=Spline1D(nodesF,[0.0,0.0,0.5,1.0,1.0],k=1)
            #fractS=fractInt(1:88)
            #radarMod.setRainFract(fractS)
            #radarMod.setdNw(zeros(88).+0.125)
            #println(fractS)
            for k=1:88
                z13obs[k]=Cfloat(10.0*log10((10^(0.1*zKu[i,j,2*k])+10^(0.1*zKu[i,j,2*k-1]))/2))
            end
            #radarMod.setZObs(copy(z13obs))
            #radarMod.setClutFAndSurf(nodes[5]+1,Int(trunc(binSf[i,j]/2)))
            #continue
            stormType.iSurf=Cint(trunc(binSf[i,j]/2))
            stormType.freezH=Cfloat(hZero[i,j]/1000.)
            stormType.nodes=pointer(nodes)
            logdnwf=Array{Cfloat}(undef,9*nmemb)
            i300=i%300+1
            j1=Cint.([j])
            i1= Cint.([i300])
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
            if relFlag[i,j]==1 || relFlag[i,j]==2
                #radarData.pia13srt=piaKa[i,j-12]/6
                #println(radarData.pia13srt," ",piaKa[i,j-12])
            end
            ccall((:ensradJulia, "./combAlg"),
                  Cvoid, (Ref{radarDataType},
                          Ref{stormStructType},Ref{retParamType},Ref{Cint},
                          Ref{radarRetType},
                          Ref{Cint},Ref{Cfloat},Ref{Cfloat},Ref{Cfloat},
                          Ref{Cint},Ref{Cfloat},Ref{Cfloat},
                          Ref{Cfloat},Ref{Cfloat
                                          },Ref{Clong},Ref{Cint},
                          Ref{Cint},Ref{Cfloat},Ref{Cint}),
                  radarData,stormType,retParam,nmu,radarRet,
                  iflag,rms1,rms2,dn,iit,
                  xscalev,randemiss,localZAngle,wfractPix,ichunk,
                  i0,j0,dZms,msFlag)
            piaRet=[]
            rrateEns=zeros(nmemb,ngates)
            z13Ens=zeros(nmemb,ngates)
            dnEns=zeros(nmemb,ngates)
            dmEns=zeros(nmemb,ngates)
            global Dm=zeros(ngates).+0.19
            rrateM=zeros(ngates)

            for im=1:nmemb
                push!(piaRet,unsafe_load(radarRet.pia13,im))
                for k=1:ngates
                    rrateEns[im,k]=unsafe_load(radarRet.rrate,(im-1)*ngates+k)
                end
                for k=1:ngates
                    z13Ens[im,k]=unsafe_load(radarRet.z13c,(im-1)*ngates+k)
                    dnEns[im,k]=unsafe_load(radarRet.log10dNw,(im-1)*ngates+k)
                    dmEns[im,k]=unsafe_load(radarRet.d0,(im-1)*ngates+k)
                end
            end
            nClutF=nodes[5]
            for k=1:nClutF
                Dm[k]=mean(dmEns[:,k])
                rrateM[k]=mean(rrateEns[:,k])
            end
            rrateM[nClutF+1:end].=rrateM[nClutF]
            xSol=Dm
            push!(dmL,Dm)
            push!(indL,[i,j])
        end
    end
end


np=pyimport("numpy")
piaRetL=np["array"](piaRetL)
#rrateM=np["mean"](np["arrary"](rrateGNL),axis=0)

#plot(rrateM2,(88:-1:1)*0.25)
#gradZ_out=zeros(88,88)
#gradZ,zSim,piaKu,rrate,piaKa=g_z(gradZ_out,Dm)
#piaKu1,piaKa1,piaW1,rrate1,zSimKu1,zSimKa1,zSimW1=simMZ(xSol)
