using Base: String, Bool, Complex
# using Core: Matrix, Vector
using LinearAlgebra
using Dates
using DelimitedFiles
using ProgressMeter
using Statistics

include("CHform.jl")
include("GSmat.jl")
include("MetroHast.jl")
include("ParTemp.jl")
include("extrafunctions.jl")

#General parameters
J = 1.
h = 0.5
L = Int8(1)

avrg = 2

Nmax = 30
AmountN = 3
Nlist = unique!(trunc.(Int8, logrange(8., Float64(Nmax), AmountN))) #convert the logrange to ints (without duplicates)
AmountN = length(Nlist)
# N = Int8(8)
# R = Int64(8)

#parameters for the parallel tempering
Nmc = Int8(7)
Nsteps = UInt64(2)
Tin = logrange(1e-5, 1e-1, Int64(Nmc))
Tfin = logrange(1e-7, 1e-3, Int64(Nmc))

#parameters for the Metropolis algo
Niter = Int64(Nmc * Nsteps)
#temperature
etaIn = 1e-3 
etaFin = 1e-5
eta = 1. ./(etaIn .+ (etaFin - etaIn)*(log.(LinRange(2., 10., Niter)) .-log(2))/(log(10)-log(2)))#step size for the gradient descend in the VQE

#Parallel Tempering
#time averages for each N
avrgtpt = zeros(Float64, AmountN)
errtpt = zeros(Float64, AmountN)
#Energy averages for each N
avrgEpt = zeros(Float64, AmountN)
errEpt = zeros(Float64, AmountN)
#VQE
#time averages for each N
avrgtmh = zeros(Float64, AmountN)
errtmh = zeros(Float64, AmountN)
#Energy averages for each N
avrgEmh = zeros(Float64, AmountN)
errEmh = zeros(Float64, AmountN)
#Memory allocations
memPT = zeros(Float64, AmountN)
memMH = zeros(Float64, AmountN)
#run
for k in 1:AmountN
    # R = Int64(Nlist[k])
    R = Int64(ceil(sqrt(Nlist[k])))
    #temp values for time average
    tmptimePT = zeros(Float64, avrg)
    tmptimeVQE = zeros(Float64, avrg)
    #temp values for energy average
    Ept = zeros(Float64, avrg)
    Emh = zeros(Float64, avrg)

    #computing the analytical energy to check the results
    m = 1:Nlist[k]
    m = m .- floor(Nlist[k]/2) #element-wise substraction
    #computation of k
    k0 = 2*pi .*m ./Nlist[k]
    #Method 2
    E0 = -1 * sum(sqrt.((-h/J .- cos.(k0)).^2 + sin.(k0).^2))
    memPT[k] = @allocated for l in 1:avrg
        Ptp = ParTemp(Nlist[k], R, L, J, h, Nmc, Nsteps, Tin, Tfin)
        tmptimePT[l] = @elapsed Ept[l] = minimum(run(Ptp))
    end
    #storing the variables
    avrgEpt[k] = mean(Ept)
    errEpt[k] = std(Ept)
    avrgtpt[k] = mean(tmptimePT)
    errtpt[k] = std(tmptimePT)

    memMH[k] = @allocated for l in 1:avrg
        nbrStep = 1
        Vtmp = VarMH(Nlist[k], R, L, J, h)
        initVMHtoHF(Vtmp)
        tmptimeVQE[l] = @elapsed step(Vtmp, eta[nbrStep])
        while (Vtmp.E - avrgEpt[l] < 1e-15) & (nbrStep < Niter)
            nbrStep += 1
            tmptimeVQE[l] += @elapsed step(Vtmp, eta[nbrStep])
            # tmptimeVQE[l] += @elapsed Emh[l], params = mhstep(H, circuit, params,  eta[step], Nlist[k])
            if Vtmp.E < E0
                # println(Vmh.Anglelist)
                # println(E0)
                break
            end
        end
    end
    #storing the variables
    avrgEmh[k] = mean(Emh)
    errEmh[k] = std(Emh)
    avrgtmh[k] = mean(tmptimeVQE)
    errtmh[k] = std(tmptimeVQE)
end

#---- Save data ----
#select the name
date = now()
fname = string("/home/rezad/Documents/Arbeit/EPFL/Labal/PDM/scripts/Julia /data/PTvsMH_" , "Nmc=", string(Nmc),"_", "Nsteps=", string(Nsteps),"_","J=", string(J),"_","h=", string(h),"_" , string(date) , ".csv")
#creating the file
touch(fname)
#opening the file in write mode
data = open(fname, "w")
#concatenating the data into a matrix
savedata = hcat(Nlist,avrgEpt, errEpt, avrgtpt, errtpt, memPT, avrgEmh, errEmh, avrgtmh, errtmh, memMH)
#writing it to the file
writedlm(fname, savedata, ",")
close(data)