
push!(PyVector(pyimport("sys")["path"]), "./")
#using .scatTables
readS=pyimport("readScattProf")

function readvar2d(fh,g1,g2,vname,nst,nt)
    var=fh["groups"]["get"](g1)["groups"]["get"](g2)["variables"]["get"](vname)
    nt1,nr=var["shape"]
    var=var["_get"]([nst,0],[nt,nr],[1,1])
    return var
end
function readvar2d_cmb(fh,g1,vname,nst,nt)
    var=fh["groups"]["get"](g1)["variables"]["get"](vname)
    nt1,nr=var["shape"]
    var=var["_get"]([nst,0],[nt,nr],[1,1])
    return var
end
function readvar3d(fh,g1,g2,vname,nst,nt)
    var=fh["groups"]["get"](g1)["groups"]["get"](g2)["variables"]["get"](vname)
    nt1,nr,nz=var["shape"]
    var=var["_get"]([nst,0,0],[nt,nr,nz],[1,1,1])
    return var
end
function readnetcdf(fname,nst,nt)
    println(fname)
    netcdf=pyimport("netCDF4")
    fh=netcdf["Dataset"](fname,"r")
    #println(fh)
    zKu=fh["groups"]["get"]("NS")["groups"]["get"]("PRE")["variables"]["get"]("zFactorMeasured")
    nt1,nr,nz=zKu["shape"]
    zKa=fh["groups"]["get"]("MS")["groups"]["get"]("PRE")["variables"]["get"]("zFactorMeasured")
    nt1,nr,nz=zKu["shape"]
    zKu=zKu["_get"]([nst,0,0],[nt,nr,nz],[1,1,1])
    nt1,nr,nz=zKa["shape"]
    zKa=zKa["_get"]([nst,0,0],[nt,nr,nz],[1,1,1])
    pType=readvar2d(fh,"NS","CSF","typePrecip",nst,nt)
    binSf=readvar2d(fh,"NS","PRE","binRealSurface",nst,nt)
    nodes=readvar3d(fh,"NS","DSD","binNode",nst,nt)
    hZero=readvar2d(fh,"NS","VER","heightZeroDeg",nst,nt)
    clutFree=readvar2d(fh,"NS","PRE","binClutterFreeBottom",nst,nt)
    reliabF=readvar2d(fh,"NS","SRT","reliabFlag",nst,nt)
    reliabKaF=readvar2d(fh,"MS","SRT","reliabFlag",nst,nt)
    piaKu=readvar2d(fh,"NS","SRT","PIAhybrid",nst,nt)
    piaKa=readvar2d(fh,"MS","SRT","PIAhybrid",nst,nt)
    lat=readvar2d_cmb(fh,"NS","Latitude",nst,nt)
    lon=readvar2d_cmb(fh,"NS","Longitude",nst,nt)
    zeroDeg=readvar2d(fh,"NS","VER","binZeroDeg",nst,nt)
    sfcType=readvar2d(fh,"NS","PRE","landSurfaceType",nst,nt)
    fh["close"]()
    return zKu,zKa,pType,binSf,hZero,nodes,clutFree,reliabF,piaKu,lat,lon,piaKa,reliabKaF,zeroDeg,sfcType
end

function readcomb(fname,ns,nt)
    netcdf=pyimport("netCDF4")
    fh=netcdf["Dataset"](fname,"r")
    sfcPrecip=readvar2d_cmb(fh,"NS","surfPrecipTotRate",ns,nt)
    sfcPrecipMS=readvar2d_cmb(fh,"MS","surfPrecipTotRate",ns,nt)
    cmbLat=readvar2d_cmb(fh,"NS","Latitude",ns,nt)
    cmbLon=readvar2d_cmb(fh,"NS","Longitude",ns,nt)
    piaSRT=readvar2d(fh,"NS","Input","piaEffective",ns,nt)
    relFlag=readvar2d(fh,"NS","Input","piaEffectiveReliabFlag",ns,nt)
    fh["close"]
    return sfcPrecip,cmbLat,cmbLon,sfcPrecipMS,piaSRT,relFlag
end


function readTables(fname)
    println(fname)
    n4=pyimport("netCDF4")
    #global rateTs,rateT, attKuT, attKaT, attKuTs, attKaTs
    fh1=n4["Dataset"](fname)#"tables_nsv_rho005_dmTs11.nc"
    zKuT=fh1["variables"]["get"]("zKuT")["_get"]([0],[240],[1])
    zKaT=fh1["variables"]["get"]("zKaT")["_get"]([0],[240],[1])
    attKuT=fh1["variables"]["get"]("attKuT")["_get"]([0],[240],[1])
    attKaT=fh1["variables"]["get"]("attKaT")["_get"]([0],[240],[1])
    dmT=fh1["variables"]["get"]("dmT")["_get"]([0],[240],[1])
    nwT=fh1["variables"]["get"]("nwT")["_get"]([0],[240],[1])
    pwcT=fh1["variables"]["get"]("pwcT")["_get"]([0],[240],[1])
    zKuTs=fh1["variables"]["get"]("zKuTs")["_get"]([0],[240],[1])
    zKaTs=fh1["variables"]["get"]("zKaTs")["_get"]([0],[240],[1])
    attKuTs=fh1["variables"]["get"]("attKuTs")["_get"]([0],[240],[1])
    attKaTs=fh1["variables"]["get"]("attKaTs")["_get"]([0],[240],[1])
    dmTs=fh1["variables"]["get"]("dmTs")["_get"]([0],[240],[1])
    nwTs=fh1["variables"]["get"]("nwTs")["_get"]([0],[240],[1])
    pwcTs=fh1["variables"]["get"]("pwcTs")["_get"]([0],[240],[1])
    attWTs=fh1["variables"]["get"]("attWTs")["_get"]([0],[240],[1])
    attWT=fh1["variables"]["get"]("attWT")["_get"]([0],[240],[1])
    zWTs=fh1["variables"]["get"]("zWTs")["_get"]([0],[240],[1])
    zWT=fh1["variables"]["get"]("zWT")["_get"]([0],[240],[1])
    
    rateTs=fh1["variables"]["get"]("rateS")["_get"]([0],[240],[1])
    rateT=fh1["variables"]["get"]("rate")["_get"]([0],[240],[1])
    fh1["close"]()
    return zKuT,zKaT,attKuT,attKaT,dmT,nwT,pwcT,zKuTs,zKaTs,attKuTs,attKaTs,dmTs,
    dmTs,nwTs,pwcTs,attWTs,attWT,zWTs,zWT,rateTs,rateT
end



xr=pyimport("xarray")
function saveTables(fname)
    nwTX=xr["DataArray"](nwT)
    pwcTX=xr["DataArray"](pwcT)
    attKaTX=xr["DataArray"](attKaT)
    attKuTX=xr["DataArray"](attKuT)
    zKaTX=xr["DataArray"](zKaT)
    zKuTX=xr["DataArray"](zKuT)
    zKaTsX=xr["DataArray"](zKaTs)
    zKuTsX=xr["DataArray"](zKuTs)
    dmTsX=xr["DataArray"](dmTs)
    dmTX=xr["DataArray"](dmT)
    pwcTsX=xr["DataArray"](pwcTs)
    attKaTsX=xr["DataArray"](attKaTs)
    attKuTsX=xr["DataArray"](attKuTs)
    attWTsX=xr["DataArray"](attWTs)
    attWTX=xr["DataArray"](attWT)
    zWTsX=xr["DataArray"](zWTs)
    zWTX=xr["DataArray"](zWT)
    nwTsX=xr["DataArray"](nwTs)
    rateX=xr["DataArray"](rateT)
    ratesX=xr["DataArray"](rateTs)
    d=Dict("nwT"=>nwTX,"pwcT"=>pwcTX,"attKaT"=>attKaTX,
           "attKuT"=>attKuTX, "dmT"=>dmTX, "zKuT"=>zKuTX, "zKaT"=>zKaTX,
           "nwTs"=>nwTsX,"pwcTs"=>pwcTsX,"attKaTs"=>attKaTsX,
           "attKuTs"=>attKuTsX, "dmTs"=>dmTsX, "zKuTs"=>zKuTsX, "zKaTs"=>zKaTsX,
           "zWTs"=>zWTsX,"zWT"=>zWTX,
           "attWTs"=>attWTsX,"attWT"=>attWTX, "rate"=>rateX,
           "rateS"=>ratesX)
    tables=xr["Dataset"](d)
    tables["to_netcdf"](fname)#"tables_nsv_rho02_dmTs11.nc")
end



function makeTables(fname)
    ns=9
    temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
    temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r,
    tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
    tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r,
    tempW,massW,fractionW,bscatW,DeqW,extW,scatW,gW,vfallW,tempW_r,massW_r,bscatW_r,DeqW_r,
    extW_r,scatW_r,gW_r,vfallW_r=scatTables.init()
    scatTables.getDmNwSF(tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
    tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r,
    temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
    temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r,ns,
    tempW,massW,fractionW,bscatW,DeqW,extW,scatW,gW,vfallW,
    tempW_r,massW_r,bscatW_r,DeqW_r,extW_r,scatW_r,gW_r,vfallW_r)
    scatTables.getDmNwR(tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu,
    tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r,
    temp,mass,fraction,bscat,Deq,ext,scat,g,vfall,
    temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r,
    tempW,massW,fractionW,bscatW,DeqW,extW,scatW,gW,vfallW,
                        tempW_r,massW_r,bscatW_r,DeqW_r,extW_r,scatW_r,gW_r,vfallW_r)
    saveTables(fname)

end
