mutable struct stormStructType
    nodes::Ptr{Cint}#{5,Cint};
    iSurf::Cint;
    freezH::Cfloat;
    rainType::Cint;
end

mutable struct radarDataType
  ngates::Cint;     #// number of gates
  z13obs::Ptr{Cfloat};    #// observed Ku-reflectivities
  z35obs::Ptr{Cfloat};    #// observed Ka-reflectivities
  xlong::Cfloat;      #// longitude -not used
  xlat::Cfloat;       #// latitude -not used
  pia13srt::Cfloat;   #// Ku-band SRT PIA
  relPia13srt::Cfloat #// Ku-band SRT PIA
  pia35srt::Cfloat;    #// Ka-band SRT PIA
  relPia35srt::Cfloat;  # // Ka-band SRT PIA
  dr::Cfloat;           # // gate size
  hh::Vector{Cfloat};  # // hh[i] is the height of gate [i]
  hfreez::Cfloat;
  sigmaZeroKu::Cfloat; #//SJM 12/3/2014
  sigmaZeroKa::Cfloat; #//SJM 12/3/2014
end


mutable struct radarRetType
    ngates::Cint;             #f1
    nMemb::Cint;              #f2
    nmfreq::Cint;             #f3
    z13c::Ptr{Cfloat};        #f4
    z35mod0::Ptr{Cfloat};     #f5
    dz35ms::Ptr{Cfloat};      #f6
    z35::Ptr{Cfloat};         #f7
    pwc::Ptr{Cfloat};         #f8
    rrate::Ptr{Cfloat};       #f9
    d0::Ptr{Cfloat};          #f10
    log10dNw::Ptr{Cfloat};    #f11
    tb::Ptr{Cfloat};          #f12
    emTb::Ptr{Cfloat};        #f13
    emis::Ptr{Cfloat};        #f14
    imu::Ptr{Cint};           #f15
    iwc::Ptr{Cint};           #f16
    icc::Ptr{Cint};           #f17
    jcc::Ptr{Cint};           #f18
    sfc_wind::Ptr{Cfloat};    #f19
    sfc_windU::Ptr{Cfloat};   #f20
    sfc_windV::Ptr{Cfloat};   #f21
    pia13::Ptr{Cfloat};       #f22
    pia35::Ptr{Cfloat};       #f23
    simSigmaZeroKu::Ptr{Cfloat}; #f24
    simSigmaZeroKa::Ptr{Cfloat}; #f25
    z35mMean::Ptr{Cfloat};       #f26
    z35mSig::Ptr{Cfloat};        #f27
    pia13mMean::Cfloat           #f28
    pia35mMean::Cfloat;          #f29
    pia13mSig::Cfloat;           #f30
    pia35mSig::Cfloat;           #f31
    ntpw::Cint;
    tpw::Ptr{Cfloat};            #f32
    tpwCldMod::Ptr{Cfloat};      #f33
    logdNw::Ptr{Cfloat};         #f34
end


struct retParamType
    wz::Cfloat;
    w13::Cfloat;
    w35::Cfloat;
    z13thresh::Cfloat;
    z35thresh::Cfloat;
end

#unsafe_load(radarData.z13obs,88) #retrieve value

nodes=Int32.([0, 69, 72, 75, 87])
stormType=stormStructType(pointer(nodes),88,4.2,2)
#using Dierckx
using Interpolations
stormType.nodes=pointer(nodes)
itp=interpolate((Float32.(nodes),),[0,25,32,28,27], Gridded(Linear()))
z13obsA=itp(0:87)
#println(z13obsA)

ngates=Cint(88)
z13obs=Float32.(z13obsA)#zeros(88).-99.9);
z35obs=Float32.(zeros(88).-99.9);
xlong=Cfloat(-99.9);
xlat=Cfloat(-99.9);
pia13srt=-99.9;
rePia13srt=-99.9;
pia35srt=-99.9;
relPia35srt=-99.9;
dr=0.25;
hh=Float32.((88:-1:1).*0.25);
hfrez=4.2;
sigmaZeroKu=-99
sigmaZeroKa=-99

nmemb=Cint(50)
nmfreq=Cint(8)
nNodes=Cint(9)
nc=50
z13c=Float32.(zeros(nmemb*ngates))
z35mod0=Float32.(zeros(nmemb*ngates))
dz35ms=Float32.(zeros(nmemb*ngates))
z35=Float32.(zeros(nmemb*ngates))
pwc=Float32.(zeros(nmemb*ngates))
rrate=Float32.(zeros(nmemb*ngates))
d0=Float32.(zeros(nmemb*ngates))
log10dNw=Float32.(zeros(nmemb*ngates))
logdNw=Float32.(zeros(nmemb*nNodes))
tb=Float32.(zeros(nmemb*2*nmfreq))
emtb=Float32.(zeros(nmemb*2*nmfreq))
emis=Float32.(zeros(nmemb*2*nmfreq))
imu=Cint.(zeros(nmemb).+3)
icc=Cint.((1:nmemb))
jcc=Cint.((1:nmemb))
sfc_wind=Float32.(zeros(nmemb))
sfc_windU=Float32.(zeros(nmemb))
sfc_windV=Float32.(zeros(nmemb))
pia13=Float32.(zeros(nmemb))
pia35=Float32.(zeros(nmemb))
simSigmaZeroKu=Float32.(zeros(nmemb))
simSigmaZeroKa=Float32.(zeros(nmemb))
z35mMean=Float32.(zeros(ngates))
z35mSig=Float32.(zeros(ngates))
tpwCldMod=Float32.(zeros(nc))
tpw=Float32.(zeros(nmemb))
iwc=Cint.(zeros(nc))
nfreq=Cint(8)

pia13mMean=Cfloat(0)           #f28
pia35mMean=Cfloat(0);          #f29
pia13mSig=Cfloat(0);           #f30
pia35mSig=Cfloat(0);
ntpw=50;

radarRet=radarRetType(ngates,nmemb,nmfreq,pointer(z13c),pointer(z35mod0),
                      pointer(dz35ms),pointer(z35),pointer(pwc),pointer(rrate),
                      pointer(d0),pointer(log10dNw),pointer(tb),
                      pointer(emtb),pointer(emis),pointer(imu),pointer(iwc),
                      pointer(icc),
                      pointer(jcc),pointer(sfc_wind),pointer(sfc_windU),
                      pointer(sfc_windV),pointer(pia13),pointer(pia35),
                      pointer(simSigmaZeroKu),pointer(simSigmaZeroKa),
                      pointer(z35mMean),pointer(z35mSig),
                      pia13mMean,pia35mMean,pia13mSig,pia35mSig,
                      ntpw,pointer(tpw),
                      pointer(tpwCldMod), pointer(logdNw))

radarData=radarDataType(ngates,pointer(z13obs),pointer(z35obs),xlong,xlat,
                        pia13srt,rePia13srt,pia35srt,relPia35srt,
                        dr,hh,hfrez,sigmaZeroKu,sigmaZeroKa)

wz=Cfloat(1.0);
w13=Cfloat(1.0);
w35=Cfloat(1.0);
z13thresh=Cfloat(0.0);
z35thresh=Cfloat(0.0);
retParam=retParamType(wz,w13,w35,z13thresh,z35thresh)
#pushfirst!(PyVector(pyimport("sys")."path"), "")
#radarData.z13obs[1]=3.0

nmu=Cint.([5])
iflag=Cint.([1])
rms1=Cfloat.(zeros(1))
rms2=Cfloat.(zeros(1))
nmu=Int32(5)
nmfreq=Int32(8)
dn=Cfloat.([-0.25])
xscalev=Cfloat.(zeros(nmemb))
randemiss=Cfloat.(zeros(nmemb*nmfreq*2))
localZAngle=Cfloat.([-0.1])
wfractPix=Cfloat.([30])
i0=Cint.([1])
j0=Cint.([1])
msFlag=Cint.([1])
iit=Cint.([3])
ichunk=Clong.([1])
dZms=Cfloat.(zeros(nmemb))
