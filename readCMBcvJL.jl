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
fs=glob.glob("DPR-CS/*HDF5")
fs2=glob.glob("CMB/*HDF5")
fnameT="tables_nsv_rho085_2_dmTs11.nc"

#parse(Int,get(fs,0)[40:47])

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
r4L=[]
#exit(0)
profL=[]
profL2=[]
nf=min(200,size(fs2)[1])
f1L=[]
f2L=[]
profCount=Dict()
fL=[["CMB/2B-CS-CONUS.GPM.DPRGMI.CORRA2018.20140608-S184323-E185201.001573.V06A.HDF5" "DPR-CS/2A-CS-CONUS.GPM.DPR.V8-20180723.20140608-S184323-E185201.001573.V06A.HDF5"],
       ["CMB/2B-CS-CONUS.GPM.DPRGMI.CORRA2018.20140403-S231912-E232750.000549.V06A.HDF5" "DPR-CS/2A-CS-CONUS.GPM.DPR.V8-20180723.20140403-S231912-E232750.000549.V06A.HDF5"],
       ["CMB/2B-CS-CONUS.GPM.DPRGMI.CORRA2018.20140528-S215907-E220745.001404.V06A.HDF5" "DPR-CS/2A-CS-CONUS.GPM.DPR.V8-20180723.20140528-S215907-E220745.001404.V06A.HDF5"],
       ["CMB/2B-CS-CONUS.GPM.DPRGMI.CORRA2018.20140420-S180115-E180954.000810.V06A.HDF5" "DPR-CS/2A-CS-CONUS.GPM.DPR.V8-20180723.20140420-S180115-E180954.000810.V06A.HDF5"],
       ["CMB/2B-CS-CONUS.GPM.DPRGMI.CORRA2018.20140611-S174100-E174934.001619.V06A.HDF5" "DPR-CS/2A-CS-CONUS.GPM.DPR.V8-20180723.20140611-S174100-E174934.001619.V06A.HDF5"],
       ["CMB/2B-CS-CONUS.GPM.DPRGMI.CORRA2018.20140429-S064709-E065546.000943.V06A.HDF5" "DPR-CS/2A-CS-CONUS.GPM.DPR.V8-20180723.20140429-S064709-E065546.000943.V06A.HDF5"],
       ["CMB/2B-CS-CONUS.GPM.DPRGMI.CORRA2018.20140510-S033048-E033926.001112.V06A.HDF5" "DPR-CS/2A-CS-CONUS.GPM.DPR.V8-20180723.20140510-S033048-E033926.001112.V06A.HDF5"],
       ["CMB/2B-CS-CONUS.GPM.DPRGMI.CORRA2018.20140603-S195546-E200412.001496.V06A.HDF5" "DPR-CS/2A-CS-CONUS.GPM.DPR.V8-20180723.20140603-S195546-E200412.001496.V06A.HDF5"]]
zProfL=[]
fL2=[]
for ifile=1:60
    push!(fL2,[fs2[ifile],fs[ifile]])
end
fL=fL2
for ifile=1:60
    zku,zka,sfcPrecip,zsfc,hzero,bcf,bst,brs,bzd,bbPeak,pType,sfcType,
    reliabF,piaSRT,locZAngle=readnetcdfj(fL[ifile][2])
    fs_1=fL[ifile][2]
    fs_2=fL[ifile][1]
    fh=nc.Dataset(fL[ifile][1],"r")
    sfcPrecip_CMB=fh.groups["NS"].variables["surfPrecipTotRate"]
    nt1,nr=sfcPrecip_CMB.shape
    sfcPrecip_CMB=sfcPrecip_CMB._get([0,0],[nt1,nr],[1,1])
    nx,ny=size(zku)
    icount=0
    for i=1:nx
        j0=1
        ibreak=0
        for j=13:37
            if pType[i,j]==2 && bzd[i,j]<bcf[i,j]-2 && maximum(zku[i,j,100:bcf[i,j]-1])>30
                if !(fs_2 in f2L)
                    push!(f2L,fs_2)
                    push!(f1L,fs_1)
                    icount=1
                    profCount[fs_2]=[icount,fs_1]
                else
                    icount+=1
                    profCount[fs_2]=[icount,fs_1]
                end
                j0=j
                sysdN1=sysdN+0.75
                dnv=(0.00005*randn(88).+sysdN1)
                piaHB,sfcRRate,
                prateBBB,prateBBT,prateBB,
                zKaSim,zKaObs,nodes,rrate,dmRet,dnRet,z13obs=hb_bb_cv(zku,zka,bcf,bst,brs,bzd,bbPeak,
                    pType,newTables,i,j,dr,dnv)

                println("Flag=$(reliabF[i,j]) pia=$(piaSRT[i,j])")
                iBB=0
                if reliabF[i,j]==-1 || reliabF[i,j]==-2
                    continue
                end
                push!(zProfL,z13obs)
                sfcRainHB,snowRateHB,rrateHB,lwcHB,dmHB,dnHB=callcmb(nodes,stormType,nmemb,nmu,brs,bst,
                    hzero,bzd,bcf,bbPeak,sfcType,sysdN,
                    zku,zka,reliabF,piaSRT,retParam,radarData,radarRet,i,j,nmfreq,locZAngle,nwdm,iBB)
                ibreak=1
                n1=nodes[5]
                if dmHB[n1]>0 && dmRet[n1]>0
                    push!(nwdmL,[dnRet[n1],dmRet[n1],dnHB[n1],dmHB[n1]])
                end
                if nodes[4]+5<nodes[5]
                    push!(profL,rrate[nodes[2]-10:nodes[4]+5])
                    push!(profL2,rrateHB[nodes[2]-10:nodes[4]+5])
                end
                #println(rrateHB[nodes[2]-15:nodes[4]+1])
                #exit(-1)
                push!(prateBBL,[prateBBB,prateBBT,prateBB,snowRateHB])
                #println(nmemb)
                #println(nodes)
                #println(bcf[i,j], " ",nmemb," ",ngates)
                #println(zKaObs[nodes[1]:nodes[5]])
                #println(sfcRainEns)
                #exit(-1)
                push!(r4L,sfcRainHB)
                n1=nodes[2]+1

                #println("$i $j $(zKaSim[n1]) $(zKaObs[n1])")
                if zKaObs[n1]>10 && zKaSim[n1]>10
                    push!(zL,[zKaSim[n1], zKaObs[n1]])
                end
                n1=nodes[5]+1
                if zKaObs[n1]>10 && zKaSim[n1]>10
                    push!(zL2,[zKaSim[n1], zKaObs[n1]])
                end
                if sfcPrecip[i,j]>=0
                    push!(r1L,sfcRRate)
                    push!(r2L,sfcPrecip[i,j])
                    push!(r3L,sfcPrecip_CMB[i,j])
                end
                #return rrateHB
                #exit(-1)
            end
        end

    end
end
prateBBL=np.array(prateBBL)
for i=1:4
    println(np.mean(prateBBL[:,i]))
end
zL=np.array(zL)
zL2=np.array(zL2)
profL=np.array(profL);
profL2=np.array(profL2);
plot(np.mean(profL,axis=0))
plot(np.mean(profL2,axis=0))
figure()
nwdmL=np.array(nwdmL)
scatter(nwdmL[:,2],nwdmL[:,1])
scatter(nwdmL[:,4],nwdmL[:,3])
xlabel("dm")
ylabel("log10(dNw)")
legend(["Scaled Tables","Iterative"])
ik=0
#print("[")
#for k in keys(profCount)
#    if profCount[k][1]>100
#        println("[\"",k, "\" \"", profCount[k][2],"\"],")
#    end
#end
#print("]")
