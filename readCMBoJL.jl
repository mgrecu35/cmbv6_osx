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
fnameT="tables_nsv_rho085_2_dmTs11.nc"


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
dnr=-0.2
dnbb=-0.2
dns=0.15
newTables=get_CMB_Tables(dnr,dnbb,dns)
function hb_bb(zku,bcf,bst,brs,bzd,bbPeak,pType,newTables,i,j,dr)
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
    for k=1:88
        z13obs[k]=Cfloat(10.0*log10((10^(0.1*zku[i,j,2*k])+10^(0.1*zku[i,j,2*k-1]))/2))
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
                    n1,n2=bisection(zKuSJ,z13c[i+1])
                    zetaS+=attKuSJ[n1]/10^(z13c[i+1]*0.1*beta)*10^(z13obs[i+1]*0.1*beta)*dr
                end
                zeta1d[i+1]=zetaS
            end
            if i<=nodes[3] && i>nodes[2]
                if z13obs[i+1]>0
                    f=(nodes[3]-i)/(nodes[3]-nodes[2])
                    #println(size(zKuSJ),size(zKuBBJ))
                    n1,n2=bisection2(zKuSJ,zKuBBJ,f,z13c[i+1])
                    attKu=f*attKuSJ[n1]+(1-f)*attKuBBJ[n1]
                    zetaS+=attKu/10^(z13c[i+1]*0.1*beta)*10^(z13obs[i+1]*0.1*beta)*dr
                end
                zeta1d[i+1]=zetaS
            end
            if i<=nodes[4] && i>nodes[3]
                if z13obs[i+1]>0
                    f=(nodes[4]-i)/(nodes[4]-nodes[3])
                    n1,n2=bisection2(zKuBBJ,zKuJ,f,z13c[i+1])
                    attKu=f*attKuBBJ[n1]+(1-f)*attKuJ[n1]
                    zetaS+=attKu/10^(z13c[i+1]*0.1*beta)*10^(z13obs[i+1]*0.1*beta)*dr
                end
                zeta1d[i+1]=zetaS
            end
            if i>nodes[4]
                if z13obs[i+1]>0
                    n1,n2=bisection(zKuJ,z13c[i+1])
                    zetaS+=attKuJ[n1]/10^(z13c[i+1]*0.1*beta)*10^(z13obs[i+1]*0.1*beta)*dr
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
            n1,n2=bisection(zKuSJ,z13c[i+1])
            retrate[i+1]=rrateSJ[n1]
            zKa=zKaSJ[n1]
            piaKa+=attKaSJ[n1]*dr
            zKaSim[i+1]=zKa-piaKa
            piaKa+=attKaSJ[n1]*dr
        end
    end
    for i=nodes[2]+1:min(nodes[5],nodes[3])
        if z13obs[i+1]>0
            f=(nodes[3]-i)/(nodes[3]-nodes[2])
            n1,n2=bisection2(zKuSJ,zKuBBJ,f,z13c[i+1])
            zKa=f*zKaSJ[n1]+(1-f)*zKaBBJ[n1]
            attKa=f*attKaSJ[n1]+(1-f)*attKaBBJ[n1]
            piaKa+=attKa*dr
            zKaSim[i+1]=zKa-piaKa
            piaKa+=attKa*dr
            pRate=f*rrateSJ[n1]+(1-f)*rrateBBJ[n1]
            retrate[i+1]=pRate
        end
    end
    for i=nodes[3]+1:min(nodes[5],nodes[4])
        if z13obs[i+1]>0
            f=(nodes[4]-i)/(nodes[4]-nodes[3])
            n1,n2=bisection2(zKuBBJ,zKuJ,f,z13c[i+1])
            zKa=f*zKaBBJ[n1]+(1-f)*zKaJ[n1]
            attKa=f*attKaBBJ[n1]+(1-f)*attKaJ[n1]
            piaKa+=attKa*dr
            zKaSim[i+1]=zKa-piaKa
            piaKa+=attKa*dr
            pRate=f*rrateBBJ[n1]+(1-f)*rrateJ[n1]
            retrate[i+1]=pRate
        end
    end
    for i=nodes[4]+1:nodes[5]
        if z13obs[i]>0
            n1,n2=bisection(zKuJ,z13c[i+1])
            zKa=zKaJ[n1]
            attKa=attKaJ[n1]
            piaKa+=attKa*dr
            zKaSim[i+1]=zKa-piaKa
            piaKa+=attKa*dr
            retrate[i+1]=rrateJ[n1]
        end
    end
    if zetaS!=zetaS
        println(piamax)
        println(zeta1d)
        println(z13obs[nodes[1]+1:nodes[5]+1])
        println(nodes)
        exit(1)
    end
    println("$(eps) $(beta) $(zetaS)")
    piaHB=-10/beta*log10(1-eps*q*beta*zetaS)
    return piaHB,retrate[nodes[5]+1],retrate[nodes[4]+1],retrate[nodes[2]+1],retrate[nodes[3]+1],zKaSim
end
dr=0.25
r1L=[]
r2L=[]
prateBBL=[]
fs=np.sort(fs)
for ifile=1:2
    zku,zka,sfcPrecip,zsfc,hzero,bcf,bst,brs,bzd,bbPeak,pType,sfcType,
    reliabF,piaSRT,locZAngle=readnetcdfj(fs[ifile])
    nx,ny=size(zku)
    for i=1:nx
        j0=1
        for j=1:ny
            if bbPeak[i,j]>0
                j0=j
                piaHB,sfcRRate,
                prateBBB,prateBBT,prateBB,zKaSim=hb_bb(zku,bcf,bst,brs,bzd,bbPeak,pType,newTables,i,j,dr)
                push!(prateBBL,[prateBBB,prateBBT,prateBB])
                println("$i $j $(piaHB) $(sfcRRate)")
                if sfcPrecip[i,j]>0
                    push!(r1L,sfcRRate)
                    push!(r2L,sfcPrecip[i,j])
                end
            end
        end
    end
end
prateBBL=np.array(prateBBL)
for i=1:3
    println(np.mean(prateBBL[:,i]))
end
