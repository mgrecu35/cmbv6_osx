using PyPlot
np=pyimport("numpy")
nwdm=np.loadtxt("AncData/nw-dm.txt")

function get_CMB_Tables(dnr,dnbb,dns)
    zKuJ,attKuJ,pwcJ,rrateJ,dmJ=Float32.(zeros(300)),Float32.(zeros(300)),Float32.(zeros(300)),Float32.(zeros(300)),Float32.(zeros(300))
    zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ=Float32.(zeros(300)),Float32.(zeros(300)),Float32.(zeros(300)),Float32.(zeros(300)),Float32.(zeros(300))
    zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ=Float32.(zeros(300)),Float32.(zeros(300)),Float32.(zeros(300)),Float32.(zeros(300)),Float32.(zeros(300))

    zKaJ,attKaJ,pwcJ_2,rrateJ_2,dmJ_2=Float32.(zeros(300)),Float32.(zeros(300)),Float32.(zeros(300)),Float32.(zeros(300)),Float32.(zeros(300))
    zKaBBJ,attKaBBJ,pwcBBJ_2,rrateBBJ_2,dmBBJ_2=Float32.(zeros(300)),Float32.(zeros(300)),Float32.(zeros(300)),Float32.(zeros(300)),Float32.(zeros(300))
    zKaSJ,attKaSJ,pwcSJ_2,rrateSJ_2,dmSJ_2=Float32.(zeros(300)),Float32.(zeros(300)),Float32.(zeros(300)),Float32.(zeros(300)),Float32.(zeros(300))
    nKuJ=Int32.([300])
    imu=Int32.([3])
    nbinsj=Int32.([0])
    nKaJ=Int32.([300])
    imu=Int32.([3])
    nbinsj_2=Int32.([0])
    ccall((:get_radarku_,"./combAlg"),Cvoid,
        (Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32},
        Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32},
        Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32},
        Ref{Int32},Ref{Int32},Ref{Int32}),
        zKuJ,attKuJ,pwcJ,rrateJ,dmJ,
        zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
        zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ,
        nKuJ,imu,nbinsj)

    ccall((:get_radarka_,"./combAlg"),Cvoid,
        (Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32},
        Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32},
        Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32},Ref{Float32},
        Ref{Int32},Ref{Int32},Ref{Int32}),
        zKaJ,attKaJ,pwcJ_2,rrateJ_2,dmJ_2,
        zKaBBJ,attKaBBJ,pwcBBJ_2,rrateBBJ_2,dmBBJ_2,
        zKaSJ,attKaSJ,pwcSJ_2,rrateSJ_2,dmSJ_2,
        nKaJ,imu,nbinsj_2)

#plot(zKuJ[1:289],dmJ[1:289])
#logdNw=0.5*(dmJ.-2.0);
    res= [0.4195213580163668, -2.8450870544121782];
    nbinsj=nbinsj[1]
#update the rain lookup table
    dNw=0.4*((res[1]*dmJ[1:289].^2 .+ res[2]*dmJ[1:289]).+3).+dnr
    zKuJ=zKuJ[1:289].+10*dNw[1:289]
    rrateJ=rrateJ[1:289].*10 .^(dNw[1:289])
    attKuJ=attKuJ[1:289].*10 .^(dNw[1:289])
    zKaJ=zKaJ[1:289].+10*dNw[1:289]
    attKaJ=attKaJ[1:289].*10 .^(dNw[1:289])
    pwcJ=pwcJ[1:289].+dNw[1:289]
#update the BB lookup table
    dNwBB=0.4*((res[1]*dmBBJ[1:289].^2 .+ res[2]*dmBBJ[1:289]).+3).+dnbb
    zKuBBJ=zKuBBJ[1:289].+10*dNwBB[1:289]
    attKuBBJ=attKuBBJ[1:289].*10 .^(dNwBB[1:289])
    zKaBBJ=zKaBBJ[1:289].+10*dNwBB[1:289]
    attKaBBJ=attKaBBJ[1:289].*10 .^(dNwBB[1:289])
    pwcBBJ=pwcBBJ[1:289].+dNwBB[1:289]
#update the snow lookup table
    dNwS=0.4*((res[1]*dmSJ[1:nbinsj].^2 .+ res[2]*dmSJ[1:nbinsj]).+3).+dns
    zKuSJ=zKuSJ[1:nbinsj].+10*dNwS[1:nbinsj]
    attKuSJ=attKuSJ[1:nbinsj].*10 .^(dNwS[1:nbinsj])
    zKaSJ=zKaSJ[1:nbinsj].+10*dNwS[1:nbinsj]
    attKaSJ=attKaSJ[1:nbinsj].*10 .^(dNwS[1:nbinsj])
    pwcSJ=pwcSJ[1:nbinsj].+dNwS[1:nbinsj]
    #plot(zKuJ[1:289],dmJ[1:289])
    #plot(zKuT[1:240],dmT[1:240])

    newTables=[zKuJ,attKuJ,pwcJ[1:289],rrateJ[1:289],dmJ[1:289],
    zKuBBJ,attKuBBJ,pwcBBJ[1:289],rrateBBJ[1:289],dmBBJ[1:289],
    zKuSJ[1:nbinsj],attKuSJ[1:nbinsj],pwcSJ[1:nbinsj],rrateSJ[1:nbinsj],dmSJ[1:nbinsj],nbinsj,
    zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ]
    return newTables
end
