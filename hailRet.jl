
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
sfcPrecip,cmbLat,cmbLon,sfcPrecipMS,piaSRT,relFlag=readcomb(fCMB,nst,nt)
nst=4660

zKu,zKa,pType,binSf,hZero,dprNodes,clutFree,reliabF,piaKu,lat,
lon,piaKa,reliabKaF,binZeroDeg,sfcType=readnetcdf(fname,nst,nt)


#fnameT="tables_nsv_rho085_2_dmTs11.nc"
#makeTables(fnameT)
#exit(1)
fnameT="tables_nsv_rho025_2_dmTs11.nc"

using .scatTables
zKuT,zKaT,attKuT,attKaT,dmT,nwT,pwcT,zKuTs,zKaTs,attKuTs,attKaTs,dmTs,
dmTs,nwTs,pwcTs,attWTs,attWT,zWTs,zWT,rateT,rateTs=readTables(fnameT)

println(minimum(dmTs))
println(minimum(dmT))
println(maximum(rateT))


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
retRate=zeros(6,88)
retRateL=[]
ic=1
ind=[[343,14],[343,15],[343,16],[344,13],[344,14],[344,15]]
#for i1=4:4 #
    
#    for ind1 in ind[1:6]
for i=1:nt
    for j=13:37
        global ic
        #i,j=ind1
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

            nodesF=copy(nodes).+1
            nodesF[1]=1
            nodesF[5]=88
            fractInt=Spline1D(nodesF,[0.0,0.0,0.5,1.0,1.0],k=1)
            fractS=fractInt(1:88)
            for k=1:88
                z13obs[k]=Cfloat(10.0*log10((10^(0.1*zKu[i,j,2*k])+10^(0.1*zKu[i,j,2*k-1]))/2))
            end

            z13c=copy(z13obs)
            
            dn1=0.20
            piaKuR,z13c,retRate1,dpiaKuRdn,rHmax=iterativep(z13obs,fractS,nodes,dn1)
            #if retRate1[nodes[5]-1]<rHmax
            #    if retRate1[nodes[5]+1]>0
            #        dn1=dn1+0.4*log(rHmax/retRate1[nodes[5]-1])
            #    end
            #end
            #piaKuR,z13c,retRate1,dpiaKuRdn,rHmax=iterativep(z13obs,fractS,nodes,dn1)
            println(" ",piaKuR, " $(dpiaKuRdn) ",piaSRT[i,j])
            push!(retRateL,retRate1)
        end
        ic+=1
    end
end

np=pyimport("numpy")
retRateL=np["array"](retRateL)
