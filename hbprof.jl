function hb_bb(zku,zka,bcf,bst,brs,bzd,bbPeak,pType,newTables,i,j,dr,dnv)
    zKuJ,attKuJ,pwcJ,rrateJ,dmJ,
    zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
    zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ,nbinsj,
    zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ,dNw,dNwBB,dNwSJ=newTables
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
    dm=zeros(88)
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
    piamax=maximum(z13obs[nodes[1]+1:nodes[5]+1])+3-z13obs[nodes[5]+1]

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
    #println("eps=",eps)
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
            #println("f=$f", " ",rrateBBJ[n1b]," ",rrateJ[n1]," ",pRate)

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
            dm[i+1]=dmJ[n1]
            dnv[i+1]+=dNw[n1]
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
    #println("dm=",dm[nodes[5]+1])
    return piaHB,retrate[nodes[5]+1],retrate[nodes[4]+1],
    retrate[nodes[2]+1],retrate[nodes[3]+1],zKaSim,z35obs,nodes,retrate,dm,dnv
end
function hb_bb_cv(zku,zka,bcf,bst,brs,bzd,bbPeak,pType,newTables,i,j,dr,dnv)
    zKuJ,attKuJ,pwcJ,rrateJ,dmJ,
    zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
    zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ,nbinsj,
    zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ,dNw,dNwBB,dNwSJ=newTables
    nodes=[0,0,0,0,0]
    nodes[1]=Int32(trunc(bst[i,j]/2))
    nodes[3]=Int32(trunc((bzd[i,j])/2))
    if bbPeak[i,j]>0
        nodes[3]=Int32(trunc(bzd[i,j]/2))
    end
    nodes[2]=nodes[3]-4
    nodes[4]=nodes[3]+2
    nodes[5]=Int32(trunc(bcf[i,j]/2))
    if nodes[1]>nodes[2]
        nodes[1]=nodes[2]-3
    end
    z13obs=zeros(88)
    z35obs=zeros(88)
    dm=zeros(88)
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
    piamax=maximum(z13obs[nodes[1]+1:nodes[5]+1])+8-z13obs[nodes[5]+1]

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
            if i<=nodes[4] && i>nodes[2]
                if z13obs[i+1]>0
                    f=(nodes[4]-i)/(nodes[4]-nodes[2])
                    if z13c[i+1]+0.5<z13c[nodes[3]+1] && z13c[i+1]-0.5>z13c[nodes[2]+1]
                        f=(z13c[nodes[3]+1]-z13c[i+1])/(z13c[nodes[3]+1]-z13c[nodes[2]+1])
                    end
                    #println(size(zKuSJ),size(zKuBBJ))
                    n1,n2=bisection2(10 .^zKuSJ,10 .^zKuJ,f,10^(z13c[i+1]-10*dnv[i+1]))
                    attKu=(f*attKuSJ[n1]+(1-f)*attKuJ[n1])*10^dnv[i+1]
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
        if eps>1
            dnv.+=0.05*log10(eps)
        end
    end
    #println("eps=",eps,"pia=",piamax)

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
    for i=nodes[2]+1:min(nodes[5],nodes[4])
        if z13obs[i+1]>0
            f=(nodes[4]-i)/(nodes[4]-nodes[2])
            #f=0.0
            n1s,n2s=bisection(zKuSJ,(z13c[i+1]-10*dnv[i+1]))
            n1,n2=bisection(zKuJ,(z13c[i+1]-10*dnv[i+1]))
            zKa=(f*zKaJ[n1]+(1-f)*zKaJ[n1s])*10^dnv[i+1]
            attKa=(f*attKaSJ[n1s]+(1-f)*attKaJ[n1])*10^dnv[i+1]
            piaKa+=attKa*dr
            zKaSim[i+1]=zKa-piaKa
            piaKa+=attKa*dr
            pRate=(f*rrateSJ[n1s]+(1-f)*rrateJ[n1])*10^dnv[i+1]
            retrate[i+1]=pRate
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
            dm[i+1]=dmJ[n1]
            dnv[i+1]+=dNw[n1]
        end
    end
    if zetaS!=zetaS
        println(piamax)
        println(zeta1d)
        println(z13obs[nodes[1]+1:nodes[5]+1])
        println(nodes)
        exit(1)
    end
    #println(nodes)
    #println("$(eps) $(beta) $(zetaS)")
    piaHB=-10/beta*log10(1-eps*q*beta*zetaS)
    #println("dm=",dm[nodes[5]+1])
    return piaHB,retrate[nodes[5]+1],retrate[nodes[4]+1],
    retrate[nodes[2]+1],retrate[nodes[3]+1],zKaSim,z35obs,nodes,retrate,dm,dnv,z13obs
end
