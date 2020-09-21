using PyCall
os=pyimport("os")
sys=pyimport("sys")
glob=pyimport("glob")
sys["path"]["append"]("/Users/mgrecu/ORO/cmbv6x")
include("src_jl/structDef.jl")
include("src_jl/scatTables.jl")
include("src_jl/readInput.jl")
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
ccall((:cloud_init_,"./combAlg"),Cvoid,(Ref{Cint},),nmfreq)
ccall((:__nbinmod_MOD_init_nbin,"./combAlg"),Cvoid,(),)
sysdN=Cfloat(-0.25) # if stratiform -=0.1
ccall((:initgeophys_,"./combAlg"),Cvoid,(Ref{Cint},Ref{Cint},Ref{Cfloat}),
      nmemb,nmfreq,sysdN)


zku,bzd,bst,brs,bcf,zsfc,hzero,lon,lat,fh,t,sfcPrecip,pType,zka,
bbPeak,sfcType,reliabF,piaSRT,locZAngle=readOrb.readOrb(fs[1])
nx,ny,nz=size(zku)
include("cmbWrapper.jl")
for i=1:nx
    for j=1:ny
        if pType[i,j]==1
            callcmb(nodes,stormType,nmemb,nmu,brs,hzero,bzd,bcf,sfcType,sysdN,
                    zku,zka,reliabF,piaSRT,retParam,radarData,radarRet,i,j,
                    nmemb,nmfreq,locZAngle)
        end
    end
end
