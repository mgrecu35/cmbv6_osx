#using Dierckx
using Interpolations
using SpecialFunctions
include("psdInt.jl")

function bisection(xn,x)
    n=size(xn)[1]
    if x<xn[1]
        return 1,1
    end
    if x>xn[n]
        return n,n
    end
    n1=1
    n2=n
    while n2-n1>1
        n12=Int(trunc((n1+n2)/2))
        if xn[n12]>x
            n2=n12
        else
            n1=n12
        end
    end
    return n1,n2
end

function bisection2(xn1,xn2,f,x)
    n=size(xn1)[1]
    #println(size(xn1)," ",size(xn2))
    if x<f*xn1[1]+(1-f)*xn2[1]
        return 1,1
    end
    if x>f*xn1[n]+(1-f)*xn2[n]
        return n,n
    end
    n1=1
    n2=n
    while n2-n1>1
        n12=Int(trunc((n1+n2)/2))
        if f*xn1[n12]+(1-f)*xn2[n12]>x
            n2=n12
        else
            n1=n12
        end
    end
    return n1,n2
end

function get_Zr(rwc,bscat_r,scat_r,ext_r,g_r,Deq_r,vfall_r,wl,dn,mu)

    wr,Z,att,dn1,nc,vdop,
    scatInt,gInt,dm,rr=nw_lambd_1M(rwc,dn,mu,bscat_r[9,:],ext_r[9,:],
                             scat_r[9,:],g_r[9,:],
                             vfall_r[:],Deq_r[:],wl)

    return Z, att,scatInt,gInt, vdop, dm, 3.6*rr

end



function get_Zm(f,rwc,ns,bscat_r,scat_r,ext_r,g_r,Deq_r,vfall_r,
                bscat,scat,ext,g,Deq,vfall,wl,dn,mu)

    wr,Zr,attr,dnr1,ncr,vdopr,
    scatIntr,gIntr,dm,rr1=nw_lambd_1M(f*rwc,dn,mu,bscat_r[9,:],ext_r[9,:],
                               scat_r[9,:],g_r[9,:],
                               vfall_r[:],Deq_r[:],wl)

    ws,Zs,atts,dns1,ncs,vdops,
    scatInts,gInts,dm,rr2=nw_lambd_1M((1-f)*rwc,dn,mu,bscat[end,ns,:],ext[end,ns,:],
                               scat[end,ns,:],g[end,ns,:],
                               vfall[ns,:],Deq[ns,:],wl)

    Z=10*log10(10^(0.1*Zr)+10^(0.1*Zs))
    att=attr+atts
    scatInt=scatIntr+scatInts
    gInt=(gIntr*scatIntr+gInts*scatInts)/scatInt
    return Z, att,scatInt,gInt,3.6*(f*rr1+(1-f)*rr2)

end


function get_Zs(rwc,ns,bscat,scat,ext,g,Deq,Deq_r,vfall,vfall_r,wl,dn,mu)
    w,Z,att,dn1,nc,vdop,
    scatInt,gInt,dm,rr=nw_lambd_1Ms(rwc,dn,mu,bscat[end,ns,:],ext[end,ns,:],
                                scat[end,ns,:],g[end,ns,:],
                                vfall[ns,:],vfall_r,Deq[ns,:],Deq_r,wl)
    println("$(w),$(rwc)")
    return Z, att, scatInt,gInt,vdop,dm,3.6*rr

end

function get_rZs(rwc,dmc,ns,bscat,scat,ext,g,Deq,Deq_r,vfall,vfall_r,wl,dn,mu)
    rwcS=rwc
    rr=0.0
    Z=-99
    scatInt=-99
    gInt=-99
    dm=-99
    att=-99
    vdop=-99
    w=-99
    for it=1:4
        w,Z,att,dn1,nc,vdop,
        scatInt,gInt,dm,rr=nw_lambd_1Ms(rwcS,dn,mu,bscat[end,ns,:],ext[end,ns,:],
                                        scat[end,ns,:],g[end,ns,:],
                                        vfall[ns,:],vfall_r,Deq[ns,:],Deq_r,wl)
        rwcS1=rwcS+0.02
        w1,Z1,att1,dn11,nc1,vdop1,
        scatInt1,gInt1,dm1,rr1=nw_lambd_1Ms(rwcS1,dn,mu,bscat[end,ns,:],ext[end,ns,:],
                                            scat[end,ns,:],g[end,ns,:],
                                            vfall[ns,:],vfall_r,Deq[ns,:],Deq_r,wl)


        gradW=log(dm1/dm)/(log(rwcS1/rwcS))

        rwcS=rwcS*exp(0.95*(log(dmc)-log(dm))*gradW/(gradW*gradW+0.0001))
        if rwcS<0.000001
            rwcS=0.000001
        end
        println("$(w) $(rwc) $(dm) $(dm1) $(dmc)")
    end
    println("$(rwc),$(rwcS),$(w) $(dm) $(dmc)")
    #exit(1)
    return Z, att, scatInt,gInt,vdop,dm,3.6*rr
end
