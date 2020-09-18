
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

include("src_jl/readInput.jl")
nst=0
nt=472
sfcPrecip,cmbLat,cmbLon,sfcPrecipMS,piaSRT,relFlag=readcomb(fCMB,nst,nt)
nst=4660

zKu,zKa,pType,binSf,hZero,dprNodes,clutFree,reliabF,piaKu,lat,
lon,piaKa,reliabKaF,binZeroDeg,sfcType=readnetcdf(fname,nst,nt)

fnameT="tables_nsv_rho025_2_dmTs11.nc"
include("src_jl/scatTables.jl")
using .scatTables
zKuT,zKaT,attKuT,attKaT,dmT,nwT,pwcT,zKuTs,zKaTs,attKuTs,attKaTs,dmTs,
dmTs,nwTs,pwcTs,attWTs,attWT,zWTs,zWT,rateT,rateTs=readTables(fnameT)

hg=0:.25:18.5
z13obsL=[]
ic=0
nodes=zeros(5)
z13obs=zeros(88)
nodesL=[]
for i=1:nt
    for j=13:37
        global ic
        if Int32(trunc(pType[i,j]/1e7))==2 &&
            dprNodes[i,j,4]<dprNodes[i,j,5] && dprNodes[i,j,1]<binZeroDeg[i,j]-4 &&
            binZeroDeg[i,j]<150 && zKu[i,j,binZeroDeg[i,j]]>10
            nodes[1]=Int32(trunc(dprNodes[i,j,1]/2))
            nodes[3]=Int32(trunc(dprNodes[i,j,3]/2))
            nodes[3]=Int32(trunc((binZeroDeg[i,j])/2))
            nodes[2]=Int32(nodes[3]-2)
            nodes[4]=Int32(nodes[3]+2)
            nodes[5]=Int32(trunc(clutFree[i,j]/2))

            for k=1:88
                z13obs[k]=Cfloat(10.0*log10((10^(0.1*zKu[i,j,2*k])+10^(0.1*zKu[i,j,2*k-1]))/2))
            end
            if z13obs[Int32(nodes[2])]<40
                continue
            end
            push!(z13obsL,z13obs[end:-1:14])
            push!(nodesL,nodes)
        end
        ic+=1
    end
end
using PyPlot
include("src_jl/psdInt_2.jl")
function getR(nodes,z13obs,drate1,drate2,dr)
    dn=0.0
    piaKu=0.0
    rrate=zeros(75)
    α=[0.000127, 0.0002109, 0.0004109,  0.0004109]
    q=0.2*log(10)
    zetaS=0
    β=0.77
    zeta=zeros(75)
    zKuC=zeros(75)
    for k=nodes[1]:-1:nodes[2]
        n1,n2=bisection(zKuTs,z13obs[k]+piaKu-10*dn)
        piaKu+=attKuTs[n1]*dr*10^dn
        rrate[k]=10^dn*rateTs[n1]
        piaKu+=attKuTs[n1]*dr*10^dn
        zetaS+=q*α[1]*10^(0.1*z13obs[k]*0.77)*dr
        zeta[k]=zetaS
    end
    for k=nodes[2]+1:-1:nodes[4]
        f=(k-nodes[4])/(nodes[2]-nodes[4]+0.0)
        rrate[k]=rrate[nodes[2]]+f*drate1
        αm=f*α[2]+(1-f)*α[3]
        zetaS+=q*αm*10^(0.1*z13obs[k]*0.77)*dr
        zeta[k]=zetaS
    end
    for k=nodes[4]+1:-1:nodes[5]
        f=1.0-(k-nodes[5])/(nodes[4]-nodes[5]+0.0)
        rrate[k]=rrate[nodes[4]]+f*drate2
        zetaS+=q*α[4]*10^(0.1*z13obs[k]*0.77)*dr
        zeta[k]=zetaS
    end
    #println(nodes)
    #println(size(z13obs))
    piamax=55.0-z13obs[nodes[5]]
    zetamax=1.0-10^(-piamax/10*β)
    if zetaS>zetamax
        ϵ=0.9999*zetamax/zetaS
        zeta=zeta.*ϵ
    else
        ϵ=1.0
    end
    println(ϵ)
    for k=nodes[1]:-1:nodes[5]
        zKuC[k]=z13obs[k]-10/β*log(1-ϵ*zeta[k])
    end
    for k=nodes[1]:-1:nodes[2]
        n1,n2=bisection(zKuTs,zKuC[k]-10*dn)
        rrate[k]=10^dn*rateTs[n1]
    end
    for k=nodes[4]:-1:nodes[5]
        n1,n2=bisection(zKuT,zKuC[k]-10*dn)
        rrate[k]=10^dn*rateT[n1]
    end
    drate1=rrate[nodes[4]]-rrate[nodes[2]]
    for k=nodes[2]+1:-1:nodes[4]+1
        f=(k-nodes[4])/(nodes[2]-nodes[4]+0.0)
        rrate[k]=f*rrate[nodes[2]]+(1-f)*rrate[nodes[4]]
    end
    #println(-10/β*log(1-ϵ*zetaS))
    return piaKu,rrate,zKuC
end
function simZ(nodes,rInt,dr)
    zKu=zeros(75)
    dn=0.0
    piaKu=0.0
    for k=nodes[1]:-1:nodes[2]
        n1,n2=bisection(rateTs,rInt[k]/10^dn)
        piaKu+=attKuTs[n1]*dr*10^dn
        zKu[k]=zKuTs[n1]-piaKu+10*dn
        piaKu+=attKuTs[n1]*dr*10^dn
    end
    for k=nodes[2]+1:-1:nodes[4]
        f=(k-nodes[4])/(nodes[2]-nodes[4]+0.0)
        n1s,n2s=bisection(rateTs,rInt[k]/10^dn)
        n1,n2=bisection(rateT,rInt[k]/10^dn)
        piaKu+=f*10^dn*attKuTs[n1s]*dr+10^dn*(1-f)*attKuT[n1]*dr
        zKu[k]=10*log10(f*10^(0.1*zKuTs[n1s]+dn)+(1-f)*10^(0.1*zKuT[n1]+dn))-piaKu
        piaKu+=10^dn*(f*attKuTs[n1s]*dr+(1-f)*attKuT[n1]*dr)
    end
    for k=nodes[4]+1:-1:nodes[5]
        n1,n2=bisection(rateT,rInt[k]/10^dn)
        piaKu+=attKuT[n1]*dr*10^dn
        zKu[k]=zKuT[n1]-piaKu+10*dn
        piaKu+=attKuT[n1]*dr*10^dn
    end
    return zKu,piaKu
end
for i=1:4:length(nodesL)
    nodesp=88 .-Int32.(nodesL[i])
    figure()
    dr1=0.0
    dr2=0.0
    dr=0.25
    piaKu,rrate,zKuC=getR(nodesp,z13obsL[i],dr1,dr2,dr)
    zKu,piaKu=simZ(nodesp,rrate,dr)
    plot(zKuC[5:end],hg[5:end])
    plot(zKu[5:end],hg[5:end])
    plot(z13obsL[i][5:end],hg[5:end])
end

np=pyimport("numpy")

