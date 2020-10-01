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
fs=glob.glob("DPR-CS/*023465.*HDF5")
fs=glob.glob("DPR-CS/*01573.*HDF5")
fs2=glob.glob("CMB/*001573.*HDF5")
fnameT="tables_nsv_rho085_2_dmTs11.nc"

#parse(Int,get(fs,0)[40:47])
netcdf=pyimport("netCDF4")

using .scatTables
zKuT,zKaT,attKuT,attKaT,dmT,nwT,pwcT,zKuTs,zKaTs,attKuTs,attKaTs,dmTs,
dmTs,nwTs,pwcTs,attWTs,attWT,zWTs,zWT,rateT,rateTs=readTables(fnameT)

#println(minimum(dmTs))
#println(minimum(dmT))
#println(maximum(rateT))
nmu=5
nmfreq=8

ccall((:readtablesliang2_,"./combAlg"),Cvoid,(Ref{Int32},Ref{Int32}),
      nmu,nmfreq)
include("procTables.jl")
dnr=-0.3
dnbb=-0.3
dns=-0.5
newTables=get_CMB_Tables(dnr,dnbb,dns)
dnr_old,dnbb_old,dns_old=-0.0,-0.0,0
newTables_old=get_CMB_Tables_old(dnr_old,dnbb_old,dns_old)

zKuJ,attKuJ,pwcJ,rrateJ,dmJ,
zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ,nbinsj,
zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ,dNw,dNwBB,dNwSJ=newTables
zKuSJ=zKuTs
zKaSJ=zKaTs
attKuSJ=attKuTs
attKaSJ=attKaTs
rrateSJ=rateTs
pwcSJ=log10.(pwcTs)
newTable2=[zKuJ,attKuJ,pwcJ,rrateJ,dmJ,
zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ,nbinsj,
zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ,dNw,dNwBB,dNwSJ]

nwdm=np.loadtxt("AncData/nw-dm.txt");
nwdm[:,2:3].-=log10(8.0e6);
include("hbprof.jl")
dr=0.25
r1L=[]
r2L=[]
r3L=[]
prateBBL=[]
fs=np.sort(fs)
fs2=np.sort(fs2)
zL=[]
zL2=[]
nwdmL=[]
nc=pyimport("netCDF4")
include("cmbWrapperHB.jl")
ccall((:cloud_init_,"./combAlg"),Cvoid,(Ref{Cint},),nmfreq)
ccall((:__nbinmod_MOD_init_nbin,"./combAlg"),Cvoid,(),)
sysdN=Cfloat(-0.25) # if stratiform -=0.1
ccall((:initgeophys_,"./combAlg"),Cvoid,(Ref{Cint},Ref{Cint},Ref{Cfloat}),
      nmemb,nmfreq,sysdN)
#nmemb=50

pickle=pyimport("pickle")
pfile=pybuiltin("open")("IPHEX0612/iphex0612.pklz","rb")
zKuL,zKaL,zXL,dL=pickle.load(pfile)
isfcL=[]
bzdL=[]
rsL=[]
zetaL=[]
ibbL=[]
for i=1:5880
    global zKu1,zKa1,zX1
    zKu1=zKuL[i,end:-1:1]
    zKa1=zKaL[i,end:-1:1]
    zX1=zXL[i,end:-1:1]
    for k=1:176
        if zKa1[k]!=zKa1[k]
            zKa1[k]=-99
        end
        if zKu1[k]!=zKu1[k]
            zKu1[k]=-99
        end
    end
    push!(bzdL,175-30)
    isfc=argmax(zX1)-3
    push!(isfcL,isfc)
    itop=isfc-3
    for k=1:isfc-3
        if zX1[k]>0 && zX1[k+1]>0 && zX1[k+2]>0
            itop=k
            break
        end
    end
    if itop<145
        zeta=np.sum(10 .^(0.1*zX1[145:isfc-2]*0.8))
    else
        zeta=0.001
    end
    ibb=0
    zmax=maximum(zKu1[145-2:145+3])
    zmin=minimum(zKu1[145-2:145+3])
    if zmax-zmin>7
        ibb=1
    end
    push!(zetaL,zeta)
    push!(rsL,itop)
    push!(ibbL,ibb)
end
plt=pyimport("matplotlib.pyplot")
plt.subplot(211)
plt.pcolormesh(zKuL',vmin=0,vmax=45,cmap="jet")
plt.plot(175 .-isfcL)
plt.plot(175 .-bzdL)
plt.plot(175 .-rsL)
plt.ylim(0,100)
plt.xlim(2000,3000)
plt.subplot(212)
plt.plot(ibbL)
plt.xlim(2000,3000)
