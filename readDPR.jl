using PyCall
push!(PyVector(pyimport("sys")["path"]), "./")
fname="../subSets/2A-CS-CONUS.GPM.DPR.V8-20180723.20190402-S172036-E172914.028936.V06A.HDF5"
readDPR=pyimport("readDPR")
