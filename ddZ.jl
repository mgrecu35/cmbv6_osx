using PyCall
os=pyimport("os")
sys=pyimport("sys")
glob=pyimport("glob")
sys["path"]["append"]("/Users/mgrecu/ORO/cmbv6x")
include("src_jl/structDef.jl")
include("src_jl/scatTables.jl")
include("src_jl/readInputDPR.jl")
include("src_jl/readInput.jl")
include("src_jl/radarMod.jl")
include("src_jl/psdInt_2.jl")
readOrb=pyimport("readOrb")
fs=glob.glob("DPR-CS/2A-*.*HDF5")

fnameT="tables_nsv_rho085_2_dmTs11.nc"

#parse(Int,get(fs,0)[40:47])
netcdf=pyimport("netCDF4")

#using .scatTables
#zKuT,zKaT,attKuT,attKaT,dmT,nwT,pwcT,zKuTs,zKaTs,attKuTs,attKaTs,dmTs,
#dmTs,nwTs,pwcTs,attWTs,attWT,zWTs,zWT,rateT,rateTs=readTables(fnameT)

#println(minimum(dmTs))
#println(minimum(dmT))
#println(maximum(rateT))
nmu=5
nmfreq=8

#ccall((:readtablesliang2_,"./combAlg"),Cvoid,(Ref{Int32},Ref{Int32}),
#      nmu,nmfreq)
#include("procTables.jl")
dnr=-0.3
dnbb=-0.3
dns=-0.5
#newTables=get_CMB_Tables(dnr,dnbb,dns)
#dnr_old,dnbb_old,dns_old=-0.0,-0.0,0
#newTables_old=get_CMB_Tables_old(dnr_old,dnbb_old,dns_old)

#zKuJ,attKuJ,pwcJ,rrateJ,dmJ,
#zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
#zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ,nbinsj,
#zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ,dNwJ,dNwBBJ,dNwSJ=newTables

#newTable2=[zKuJ,attKuJ,pwcJ,rrateJ,dmJ,
#zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
#zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ,nbinsj,
#zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ,dNwJ,dNwBBJ,dNwSJ]
np=pyimport("numpy")
nwdm=np.loadtxt("AncData/nw-dm.txt");
nwdm[:,2:3].-=log10(8.0e6);
include("hbprof.jl")
dr=0.25

fs=np.sort(fs)

nc=pyimport("netCDF4")
include("cmbWrapperHB.jl")
#ccall((:cloud_init_,"./combAlg"),Cvoid,(Ref{Cint},),nmfreq)
#ccall((:__nbinmod_MOD_init_nbin,"./combAlg"),Cvoid,(),)
sysdN=Cfloat(-0.25) # if stratiform -=0.1
#ccall((:initgeophys_,"./combAlg"),Cvoid,(Ref{Cint},Ref{Cint},Ref{Cfloat}),
#      nmemb,nmfreq,sysdN)
#nmemb=50


zKuL=[]
zKaL=[]
pickle=pyimport("pickle")
fp=pybuiltin("open")("kMeans.pklz","rb")
kMeans=pickle.load(fp)
zKaCMBL1=[]
zKaCMBL2=[]
zKuCMBL1=[]
nodesL=[]
piaSRTL=[]
piaL=[]
isurfL=[]
sfcPrecipL=[]
nubfL=[]
sfcPrecipL2=[]
datL=[]
for ifile=6:300
    global zKu1
    zku,zka,sfcPrecip,zsfc,hzero,bcf,bst,brs,bzd,bbPeak,pType,sfcType,
    reliabF,piaSRT,locZAngle,piaHB,precip3D,zkuc=readnetcdfj(fs[ifile])
    a=findall(pType.==2)
    n1=size(a)[1]
    nx=size(piaHB)[1]
    for i1=2:n1-1
        i,j=a[i1][1],a[i1][2]
        zKu1=zku[i:i,j,1:168]
        zKu1[zKu1.<0].=0
        node2=bzd[i,j]
        if hzero[i,j]>2000 && hzero[i,j]<4500 && bcf[i,j]>168 && j>13 && j<37 &&
            i>1 && i<nx-1
            iClass=1#kMeans.predict(zKu1)[1]
            nubf=-10*log10(mean(10.0 .^(-0.9*piaHB[i-1:i+1,j-1:j+1])))/
                                mean(piaHB[i-1:i+1,j-1:j+1])/9.0
            if iClass>0 && nubf>0 && node2>0 && piaSRT[i,j]>0
                push!(zKuL,zku[i,j,:])
                push!(zKaL,zka[i,j-12,:])
                iBB=1
                pType1=1
                if zku[i,j,bcf[i,j]]>10 && zku[i,j,bcf[i,j]-2]>10 &&
                    (reliabF[i,j]==1 || reliabF[i,j]==2)

                    #push!(datL,[zku[i,j,bcf[i,j]],(zku[i,j,bbPeak[i,j]]-30)/10,
                    #max(zka[i,j-12,bbPeak[i,j]],0)/10,ddZ,hzero[i,j]/1000.,piaSRT[i,j]/nubf,
                    #(bcf[i,j]-bbPeak[i,j])/(0+brs[i,j]-bbPeak[i,j]),sfcPrecip[i,j]])
                    for k=bzd[i,j]+3:bcf[i,j]
                        if (np.random.random()<0.1 || k==bcf[i,j]) &&
                             zkuc[i,j,k]>10 && zku[i,j,k]>10
                             isurf=0
                             if k==bcf[i,j]
                                 isurf=1
                             else
                                 isurf=0
                             end

                            push!(datL,[zku[i,j,k],(zku[i,j,node2]-30)/10,
                            hzero[i,j]/1000.,piaSRT[i,j],
                            (k-node2)/(0.0+brs[i,j]-node2),
                            precip3D[i,j,k],zkuc[i,j,k]-zku[i,j,k],isurf])
                        end
                    end
                end
            end
        end
    end
    if(size(zKaCMBL1)[1]>15000)
        println(fs[ifile])
        break
    end
end
datL=np.array(datL)
fh=pybuiltin("open")("belowZc.pklz","wb")
pickle.dump(datL,fh)
fh.close()
