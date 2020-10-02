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


ns=9
temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r,
tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r,
tempW,massW,fractionW,bscatW,DeqW,extW,scatW,gW,vfallW,tempW_r,massW_r,bscatW_r,DeqW_r,
extW_r,scatW_r,gW_r,vfallW_r,
tempX,massX,fractionX,bscatX,DeqX,extX,scatX,gX,vfallX,tempX_r,massX_r,bscatX_r,DeqX_r,
extX_r,scatX_r,gX_r,vfallX_r=scatTables.init()
freqKu=13.8
wlKu=300/freqKu
freqX=9.6
wlX=300/freqX
zObs1=10.0
dnRet=0.27
mu=2.0
rwcL=[]
zXT=[]
attXT=[]
for i=-12:0.25:60
    zObs1=i+0.0
    retZ=get_fZr(zObs1,bscatKu_r,scatKu_r,extKu_r,gKu_r,DeqKu_r,vfallKu_r,wlKu,dnRet,mu)
    rwc, Z, att,scatInt,gInt, vdop, dm, rrate=retZ
    zX, attX,scatIntX,gIntX, vdopX,dmX,rrate=get_Zr(rwc,bscatX_r,scatX_r,extX_r,gX_r,DeqX_r,vfallX_r,wlX,dnRet,mu)
    push!(rwcL,rwc)
    push!(zXT,zX)
    push!(attXT,attX*4.343)
end
swcL=[]
dnRet=0.5
zXsT=[]
attXsT=[]
for i=-12:0.25:51
    zObs1=i;
    ns1=9
    retZ=get_fZs(zObs1,ns1,bscatKu,scatKu,extKu,gKu,DeqKu,DeqKu_r,vfallKu,vfallKu_r,wlKu,dnRet,mu)
    rwc, Z, att,scatInt,gInt, vdop, dm,rrate=retZ
    push!(swcL,rwc)
    zX, attX,scatIntX,gIntX, vdopX,dmX,rrateX=
    get_rZs(rwc,dm,ns1,bscatX,scatX,extX,gX,DeqX,DeqX_r,vfallX,vfallX_r,wlX,dnRet,mu);
    push!(zXsT,zX)
    push!(attXsT,attX*4.343)
end

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
zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ=newTables_old


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
include("gaussNewton.jl")
ccall((:cloud_init_,"./combAlg"),Cvoid,(Ref{Cint},),nmfreq)
ccall((:__nbinmod_MOD_init_nbin,"./combAlg"),Cvoid,(),)
sysdN=Cfloat(-0.25) # if stratiform -=0.1
ccall((:initgeophys_,"./combAlg"),Cvoid,(Ref{Cint},Ref{Cint},Ref{Cfloat}),
      nmemb,nmfreq,sysdN)
#nmemb=50
matplotlib=pyimport("matplotlib")
pickle=pyimport("pickle")
pfile=pybuiltin("open")("IPHEX0612/iphex0612.pklz","rb")
zKuL,zKaL,zXL,dL=pickle.load(pfile)
isfcL=[]
bzdL=[]
rsL=[]
zetaL=[]
ibbL=[]
zL1=[]
zL2=[]
rL=[]
sfcD=[]
indL=[]
nodesL=[]
rrateL=[]
pTypeL=[]
#for i=1000:5000
iret=1
if iret==1
for i=1:5000
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
    for k=70:isfc-3
        #global itop
        #println(zX1[k:k+3])
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
    relFlag=0
    piaSRT=-1
    sysdN1=sysdN+0.35
    dnv=(0.00005*randn(88).+sysdN1)
    sfcType=1
    hzero=30*0.125
    bcf=isfc-4
    bzd=143
    bst=itop
    bbPeak=bzd
    nmfreq=8
    brs=isfc
    pType=0
    if bst<140
        if zKu1[148]<40
            pType=1
        else
            pType=2
        end
        #println(pType)
        #bbPeak= np.argmax(zX1[bbPeak-3:bbPeak+5])+bbPeak-3
        bzd=bbPeak
        ret=callcmb(brs,bst,hzero,bzd,bcf,bbPeak,sfcType,sysdN,pType,
        zKu1,zKa1,zX1,relFlag,piaSRT,nmfreq,nwdm,ibb,zXT,attXT,zXsT,attXsT)
        sfcRate,snowRate, rrate,lwc_ret, dm, log10dNw,z35Sim,zKaSimEns1,
        z13obs,z35obs,nodes,zXtrue,zXobs=ret
        push!(zL1,z35Sim)
        push!(zL2,z35obs)
        push!(rrateL,rrate)
        push!(indL,i)
        push!(nodesL,[brs,bst,bzd,bcf,bbPeak])
        push!(pTypeL,pType)
        if z35obs[nodes[5]]>0 && z35Sim[nodes[5]]>0
            push!(sfcD,[z35obs[nodes[5]],z35Sim[nodes[5]],zXtrue[nodes[5]],zXobs[nodes[5]],sfcRate])

        end
    end
end
nodesL=np.array(nodesL)
#plt=pyimport("matplotlib.pyplot")
plt.figure()
plt.subplot(211)
plt.pcolormesh(zKuL',vmin=0,vmax=45,cmap="jet")
plt.plot(175 .-isfcL)
plt.plot(175 .-bzdL)
plt.plot(175 .-rsL)
plt.ylim(0,100)
plt.xlim(2600,3900)
plt.subplot(212)
plt.pcolormesh(zKaL',vmin=0,vmax=40,cmap="jet")
plt.plot(175 .-isfcL)
plt.plot(175 .-bzdL)
plt.plot(175 .-rsL)
plt.ylim(0,100)
plt.xlim(2600,3900)
sfcD=np.array(sfcD)

a=findall(sfcD[:,4].>0)
println(np.corrcoef(sfcD[a,3],sfcD[a,4]))
println(np.corrcoef(sfcD[a,1],sfcD[a,2]))
println(np.mean(sfcD[a,:],axis=0))
rrateL=np.array(rrateL)

figure()
plt.pcolormesh((rrateL.+1e-9)',vmin=0.1,vmax=100,cmap="jet",norm=matplotlib.colors.LogNorm())
#plt.plot(175 .-isfcL)
#plt.plot(175 .-bzdL)
#plt.plot(175 .-rsL)
#plt.plot(nodesL[:,end]./2)
plt.xlim(2200,3200)
plt.ylim(88,40)
plt.ylabel("Range bin")
plt.xlabel("Profile#")
cbar=plt.colorbar()
cbar.ax.set_title("mm/h")
plt.savefig("retrievedRateRate.png")



plt.figure()
zL1=np.array(zL1)
zL2=np.array(zL2)
plt.subplot(211)
plt.pcolormesh((zL1.+0)',vmin=0,vmax=40,cmap="jet")
#plt.plot(nodesL[:,end]./2)
plt.xlim(2200,3200)
plt.ylim(88,40)
plt.ylabel("Range bin")
plt.colorbar()
plt.subplot(212)
plt.pcolormesh(zL2',vmin=0,vmax=40,cmap="jet")
plt.xlim(2200,3200)
plt.ylim(88,40)
plt.ylabel("Range bin")
plt.xlabel("Profile#")
plt.colorbar()
plt.savefig("zKaObsAndSims.png")

plt.figure()
plt.plot(zL1[2800,end:-1:1],(1:88).*2)
plt.plot(zKuL[3643,:],1:176)
plt.plot(zKaL[3643,:],1:176)
plt.xlim(0,50)
plt.ylim(0,100)
end
#[1.0 0.884354; 0.884354 1.0]
#[1.0 0.823746; 0.823746 1.0]
#[21.0718, 21.6027, 26.8847, 30.1095, 2.77357]
