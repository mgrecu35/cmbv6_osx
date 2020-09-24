using PyCall
os=pyimport("os")
sys=pyimport("sys")
glob=pyimport("glob")
sys["path"]["append"]("/Users/mgrecu/ORO/cmbv6x")
include("src_jl/structDef.jl")
include("src_jl/scatTables.jl")
include("src_jl/readInputDPR.jl")
include("src_jl/radarMod.jl")
readOrb=pyimport("readOrb")
fs=glob.glob("DPR-CS/*HDF5")
fnameT="tables_nsv_rho085_2_dmTs11.nc"


#using .scatTables
#zKuT,zKaT,attKuT,attKaT,dmT,nwT,pwcT,zKuTs,zKaTs,attKuTs,attKaTs,dmTs,
#dmTs,nwTs,pwcTs,attWTs,attWT,zWTs,zWT,rateT,rateTs=readTables(fnameT)

#println(minimum(dmTs))
#println(minimum(dmT))
#println(maximum(rateT))
nmu=5
nmfreq=8

ccall((:readtablesliang2_,"./combAlg"),Cvoid,(Ref{Int32},Ref{Int32}),
      nmu,nmfreq)

zKuJ,attKuJ,pwcJ,rrateJ,dmJ=zeros(300),zeros(300),zeros(300),zeros(300),zeros(300)
zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ=zeros(300),zeros(300),zeros(300),zeros(300),zeros(300)
zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ=zeros(300),zeros(300),zeros(300),zeros(300),zeros(300)
nKuJ=Int32.([300])
imu=Int32.([3])

ccall((:get_radarku_,"./combAlg"),Cvoid,
    (Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32),
    Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32),
    Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32},
    Ref{Int32},Ref{Int32})
    zKuJ,attKuJ,pwcJ,rrateJ,dmJ,
    zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
    zKuSJ,attKuSJ,pwcSJ,rrateJS,dmSJ,
    nKuJ,imu)
ccall((:cloud_init_,"./combAlg"),Cvoid,(Ref{Cint},),nmfreq)
ccall((:__nbinmod_MOD_init_nbin,"./combAlg"),Cvoid,(),)
sysdN=Cfloat(-0.25) # if stratiform -=0.1
ccall((:initgeophys_,"./combAlg"),Cvoid,(Ref{Cint},Ref{Cint},Ref{Cfloat}),
      nmemb,nmfreq,sysdN)



include("cmbWrapper.jl")
r1L=[]
r2L=[]
np=pyimport("numpy")
#for ifile=1:100
#    zku,zka,sfcPrecip,zsfc,hzero,bcf,bst,brs,bzd,bbPeak,pType,sfcType,reliabF,piaSRT,locZAngle=readnetcdfj(fs[ifile])
#    print(ifile," ",fs[ifile])
#end
#exit(1)
for ifile=1:2
    println("one")
#zku,bzd,bst,brs,bcf,zsfc,hzero,lon,lat,t,sfcPrecip,pType,zka,
#bbPeak,sfcType,reliabF,piaSRT,locZAngle=readOrb.readOrb(fs[ifile])
    #zku,zka,sfcPrecip,zsfc,hzero,bcf,bst,brs,bzd,bbPeak,pType,
    #sfcType,reliabF,piaSRT,locZAngle=readnetcdfj(fs[ifile])
    zku,zka,sfcPrecip,zsfc,hzero,bcf,bst,brs,bzd,bbPeak,pType,sfcType,
    reliabF,piaSRT,locZAngle=readnetcdfj(fs[ifile])
    println("two")
    nx,ny,nz=size(zku)
    for i=1:nx
        for j=1:ny
            if pType[i,j]==1
                callcmb(nodes,stormType,nmemb,nmu,brs,bst,hzero,bzd,bcf,bbPeak,sfcType,sysdN,
                        zku,zka,reliabF,piaSRT,retParam,radarData,radarRet,i,j,
                        nmfreq,locZAngle)
                sfcRateEns=zeros(nmemb)
                for im=1:nmemb
                    k=Int(trunc(bcf[i,j]/2))+1
                    sfcRateEns[im]=unsafe_load(radarRet.rrate,(im-1)*ngates+k)
                end
                k=Int(trunc(bcf[i,j]/2))
                println("dpr=",sfcPrecip[i,j])
                if sum(sfcRateEns)/nmemb>10
                    exit(-1)
                end
                #println(sum(sfcRateEns)/nmemb, " ",zku[i,j,k])
                if(sfcRateEns[1]>=0.0)
                    push!(r1L,sum(sfcRateEns)/nmemb)
                    push!(r2L,sfcPrecip[i,j])
                end
            end
        end
    end
    println("ratio=",sum(r1L)/sum(r2L))
end
