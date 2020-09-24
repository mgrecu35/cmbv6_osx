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
dns=-0.3
newTables=get_CMB_Tables(dnr,dnbb,dns)
dnr_old,dnbb_old,dns_old=-0.0,-0.0,0
newTables_old=get_CMB_Tables_old(dnr_old,dnbb_old,dns_old)
nwdm=np.loadtxt("AncData/nw-dm.txt");
nwdm[:,2:3].-=log10(8.0e6);
function hb_bb(zku,zka,bcf,bst,brs,bzd,bbPeak,pType,newTables,i,j,dr,dnv)
    zKuJ,attKuJ,pwcJ,rrateJ,dmJ,
    zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
    zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ,nbinsj,
    zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ=newTables
    nodes=[0,0,0,0,0]
    nodes[1]=Int32(trunc(bst[i,j]/2))
    nodes[3]=Int32(trunc((bzd[i,j])/2))
    if bbPeak[i,j]>0
        nodes[3]=Int32(trunc(bbPeak[i,j]/2))
    end
    nodes[2]=nodes[3]-3
    nodes[4]=nodes[3]+2
    nodes[5]=Int32(trunc(bcf[i,j]/2))
    if nodes[1]>nodes[2]
        nodes[1]=nodes[2]-3
    end
    z13obs=zeros(88)
    z35obs=zeros(88)
    for k=1:88
        z13obs[k]=Cfloat(10.0*log10(1e-9+(10^(0.1*zku[i,j,2*k])+10^(0.1*zku[i,j,2*k-1]))/2))
        z35obs[k]=Cfloat(10.0*log10(1e-9+(10^(0.1*zka[i,j-12,2*k])+10^(0.1*zka[i,j-12,2*k-1]))/2))
    end
    z13c=copy(z13obs)
    iSurf=Int32(trunc(brs[i,j]/2))
    zeta1d=zeros(88)
    beta=0.76
    if z13obs[nodes[5]+1]<0
        z13obs[nodes[5]+1]=0
    end
    piamax=maximum(z13obs[nodes[1]+1:nodes[5]+1])+5-z13obs[nodes[5]+1]

    zetaS=0.0
    q=0.2*log(10)
    eps=1
    piaKa=0.0
    zKaSim=zeros(88).-99
    for it=1:2
        zetaS=0.0
        piaKa=0.0
        for i=nodes[1]:nodes[5]
            if i<=nodes[2]
                if z13obs[i+1]>0
                    n1,n2=bisection(zKuSJ,z13c[i+1]-10*dnv[i+1])
                    zetaS+=attKuSJ[n1]*10^dnv[i+1]/10^(z13c[i+1]*0.1*beta)*10^(z13obs[i+1]*0.1*beta)*dr
                end
                zeta1d[i+1]=zetaS
            end
            if i<=nodes[3] && i>nodes[2]
                if z13obs[i+1]>0
                    f=(nodes[3]-i)/(nodes[3]-nodes[2])
                    if z13c[i+1]+0.5<z13c[nodes[3]+1] && z13c[i+1]-0.5>z13c[nodes[2]+1]
                        f=(z13c[nodes[3]+1]-z13c[i+1])/(z13c[nodes[3]+1]-z13c[nodes[2]+1])
                    end
                    #println(size(zKuSJ),size(zKuBBJ))
                    n1,n2=bisection2(10 .^zKuSJ,10 .^zKuBBJ,f,10^(z13c[i+1]-10*dnv[i+1]))
                    attKu=(f*attKuSJ[n1]+(1-f)*attKuBBJ[n1])*10^dnv[i+1]
                    zetaS+=attKu/10^(z13c[i+1]*0.1*beta)*10^(z13obs[i+1]*0.1*beta)*dr
                end
                zeta1d[i+1]=zetaS
            end
            if i<=nodes[4] && i>nodes[3]
                if z13obs[i+1]>0
                    f=(nodes[4]-i)/(nodes[4]-nodes[3])
                    n1,n2=bisection2(zKuBBJ,zKuJ,f,z13c[i+1]-10*dnv[i+1])
                    attKu=(f*attKuBBJ[n1]+(1-f)*attKuJ[n1])*10^dnv[i+1]
                    zetaS+=attKu/10^(z13c[i+1]*0.1*beta)*10^(z13obs[i+1]*0.1*beta)*dr
                end
                zeta1d[i+1]=zetaS
            end
            if i>nodes[4]
                if z13obs[i+1]>0
                    n1,n2=bisection(zKuJ,z13c[i+1]-10*dnv[i+1])
                    zetaS+=attKuJ[n1]*10^dnv[i+1]/10^(z13c[i+1]*0.1*beta)*10^(z13obs[i+1]*0.1*beta)*dr
                end
                zeta1d[i+1]=zetaS
            end
        end
        eps=min(1,(1-10^(-0.1*beta*piamax))/(q*zetaS))
        for i=nodes[1]:nodes[5]
            z13c[i+1]=z13obs[i+1]-10/beta*log10(1-eps*q*beta*zeta1d[i+1])
        end
    end
    retrate=zeros(88)

    if nodes[5]>87
        nodes[5]=87
    end

    for i=nodes[1]:min(nodes[5],nodes[2])
        if z13obs[i+1]>0
            n1,n2=bisection(zKuSJ,z13c[i+1]-10*dnv[i+1])
            retrate[i+1]=rrateSJ[n1]*10^dnv[i+1]
            zKa=zKaSJ[n1]*10^dnv[i+1]
            piaKa+=attKaSJ[n1]*dr*10^dnv[i+1]
            zKaSim[i+1]=zKa-piaKa
            piaKa+=attKaSJ[n1]*dr*10^dnv[i+1]
        end
    end
    for i=nodes[2]+1:min(nodes[5],nodes[3])
        if z13obs[i+1]>0
            f=(nodes[3]-i)/(nodes[3]-nodes[2])
            f=f
            if z13c[i+1]+0.5<z13c[nodes[3]+1] && z13c[i+1]-0.5>z13c[nodes[2]+1]
            end
            f=0.5
            if i==nodes[2]+1
                f=0.5
            end
            n1s,n2s=bisection(zKuSJ,(z13c[i+1]-10*dnv[i+1]))
            n1b,n2b=bisection(zKuBBJ,(z13c[i+1]-10*dnv[i+1]))
            n1,n2=bisection(zKuJ,(z13c[i+1]-10*dnv[i+1]))
            zKa=(f*zKaJ[n1s]+(1-f)*zKaJ[n1b])*10^dnv[i+1]
            attKa=(f*attKaJ[n1s]+(1-f)*attKaJ[n1b])*10^dnv[i+1]
            piaKa+=attKa*dr
            zKaSim[i+1]=zKa-piaKa
            piaKa+=attKa*dr
            pRate=(f*rrateBBJ[n1b]+(1-f)*rrateJ[n1])*10^dnv[i+1]
            if i==nodes[2]+1
                pRate=0.5*pRate+0.5*rrateSJ[n1s]*10^dnv[i+1]
            end
            retrate[i+1]=pRate
        end
    end
    for i=nodes[3]+1:min(nodes[5],nodes[4])
        if z13obs[i+1]>0
            f=(nodes[4]-i+0.0)/(nodes[4]-nodes[3])
            f=0.15
            if i==nodes[3]+2
                f=0.05
            end
            n1b,n2b=bisection(zKuBBJ,z13c[i+1]-10*dnv[i+1])
            n1,n2=bisection(zKuJ,z13c[i+1]-10*dnv[i+1])
            zKa=(f*zKaBBJ[n1b]+(1-f)*zKaJ[n1])*10^dnv[i+1]
            attKa=(f*attKaBBJ[n1b]+(1-f)*attKaJ[n1])*10^dnv[i+1]
            piaKa+=attKa*dr
            zKaSim[i+1]=zKa-piaKa
            piaKa+=attKa*dr
            pRate=(f*rrateBBJ[n1b]+(1-f)*rrateJ[n1])*10^dnv[i+1]
            println("f=$f", " ",rrateBBJ[n1b]," ",rrateJ[n1]," ",pRate)

            retrate[i+1]=pRate
            #retrate[i+1]=rrateJ[n1]*10^dnv[i+1]
        end
    end
    for i=nodes[4]+1:nodes[5]
        if z13obs[i+1]>0
            n1,n2=bisection(zKuJ,z13c[i+1]-10*dnv[i+1])
            zKa=zKaJ[n1]
            attKa=attKaJ[n1]*10^dnv[i+1]
            piaKa+=attKa*dr
            zKaSim[i+1]=zKa-piaKa
            piaKa+=attKa*dr
            retrate[i+1]=rrateJ[n1]*10^dnv[i+1]
        end
    end
    if zetaS!=zetaS
        println(piamax)
        println(zeta1d)
        println(z13obs[nodes[1]+1:nodes[5]+1])
        println(nodes)
        exit(1)
    end
    #println("$(eps) $(beta) $(zetaS)")
    piaHB=-10/beta*log10(1-eps*q*beta*zetaS)
    return piaHB,retrate[nodes[5]+1],retrate[nodes[4]+1],
    retrate[nodes[2]+1],retrate[nodes[3]+1],zKaSim,z35obs,nodes,retrate
end
dr=0.25
r1L=[]
r2L=[]
r3L=[]
prateBBL=[]
fs=np.sort(fs)
fs2=np.sort(fs2)
zL=[]
zL2=[]
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
for ifile=2:20
    zku,zka,sfcPrecip,zsfc,hzero,bcf,bst,brs,bzd,bbPeak,pType,sfcType,
    reliabF,piaSRT,locZAngle=readnetcdfj(get(fs,ifile-1))
    fh=nc.Dataset(get(fs2,ifile-1),"r")
    sfcPrecip_CMB=fh.groups["NS"].variables["surfPrecipTotRate"]
    nt1,nr=sfcPrecip_CMB.shape
    sfcPrecip_CMB=sfcPrecip_CMB._get([0,0],[nt1,nr],[1,1])
    nx,ny=size(zku)
    for i=1:nx
        j0=1
        ibreak=0
        for j=13:37
            if bbPeak[i,j]>0
                j0=j
                sysdN1=sysdN+0.25
                dnv=(0.00005*randn(88).+sysdN1)
                piaHB,sfcRRate,
                prateBBB,prateBBT,prateBB,
                zKaSim,zKaObs,nodes,rrate=hb_bb(zku,zka,bcf,bst,brs,bzd,bbPeak,pType,newTables,i,j,dr,dnv)
                println(rrate[nodes[2]:nodes[4]+1])
                sfcRainHB,snowRateHB,rrateHB=callcmb(nodes,stormType,nmemb,nmu,brs,bst,hzero,bzd,bcf,bbPeak,sfcType,sysdN,
                    zku,zka,reliabF,piaSRT,retParam,radarData,radarRet,i,j,nmfreq,locZAngle,nwdm)
                ibreak=1
                if nodes[4]+5<nodes[5]
                    push!(profL,rrate[nodes[2]-10:nodes[4]+5])
                    push!(profL2,rrateHB[nodes[2]-10:nodes[4]+5])
                end
                println(rrateHB[nodes[2]-15:nodes[4]+1])
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
