using Base: String, Bool, Complex
# using Core: Matrix, Vector
using LinearAlgebra
using Dates
using DelimitedFiles
using ProgressMeter

include("CHform.jl")
include("GSmat.jl")
include("MetroHast.jl")
include("ParTemp.jl")
include("extrafunctions.jl")
include("vqe.jl")

#General parameters
J = 1.
h = 0.5
L = Int8(1)

avrg = 3

Nmax = 20
AmountN = 5
Nlist = unique!(trunc.(Int8, logrange(8., Float64(Nmax), AmountN))) #convert the logrange to ints (without duplicates)
AmountN = length(Nlist)
# N = Int8(8)
# R = Int64(8)

#parameters for the parallel tempering
Nmc = Int8(3)
Nsteps = UInt64(2000)
Tin = logrange(1e-5, 1e-1, Int64(Nmc))
Tfin = logrange(1e-7, 1e-3, Int64(Nmc))

#parameters for the VQE
Niter = 100
etaIn = 5e-2
etaFin = 1e-4 #step size
eta = (etaIn .+ (etaFin - etaIn)*(log.(LinRange(2., 10., Niter)) .-log(2))/(log(10)-log(2)))#step size for the gradient descend in the VQE

#Parallel Tempering
#time averages for each N
avrgtpt = zeros(Float64, AmountN)
errtpt = zeros(Float64, AmountN)
#Energy averages for each N
avrgEpt = zeros(Float64, AmountN)
errEpt = zeros(Float64, AmountN)
#VQE
#time averages for each N
avrgtvqe = zeros(Float64, AmountN)
errtvqe = zeros(Float64, AmountN)
#Energy averages for each N
avrgEvqe = zeros(Float64, AmountN)
errEvqe = zeros(Float64, AmountN)
#Memory allocations
memPT = zeros(Float64, AmountN)
memVQE = zeros(Float64, AmountN)
#run
for k in 1:AmountN
    # R = Int64(Nlist[k])
    R = Int64(ceil(sqrt(Nlist[k])))
    #temp values for time average
    tmptimePT = zeros(Float64, avrg)
    tmptimeVQE = zeros(Float64, avrg)
    #temp values for energy average
    Ept = zeros(Float64, avrg)
    Evqe = zeros(Float64, avrg)

    #computing the analytical energy to check the results
    m = 1:Nlist[k]
    m = m .- floor(Nlist[k]/2) #element-wise substraction
    #computation of k
    k0 = 2*pi .*m ./Nlist[k]
    #Method 2
    E0 = -1 * sum(sqrt.((-h/J .- cos.(k0)).^2 + sin.(k0).^2))
    memPT[k] = @allocated for l in 1:avrg
        Ptp = ParTemp(Nlist[k], R, L, J, h, Nmc, Nsteps, Tin, Tfin)
        tmptimePT[l] = @elapsed Ept[l] = minimum(run(Ptp, Nsteps))
    end
    #storing the variables
    avrgEpt[k] = mean(Ept)
    errEpt[k] = std(Ept)
    avrgtpt[k] = mean(tmptimePT)
    errtpt[k] = std(tmptimePT)

    memVQE[k] = @allocated for l in 1:avrg
        H, circuit, params = initvqe(Nlist[k], L, J, h)
        tmptimeVQE[l] = @elapsed Evqe[l], params = vqestep(H, circuit, params, eta[1], Nlist[k])
        step = 1
        while (Evqe[l]-avrgEpt[l] < 1e-15) & (step < Niter-1)
            step += 1
            tmptimeVQE[l] += @elapsed Evqe[l], params = vqestep(H, circuit, params,  eta[step], Nlist[k])
        end
    end
    #storing the variables
    avrgEvqe[k] = mean(Evqe)
    errEvqe[k] = std(Evqe)
    avrgtvqe[k] = mean(tmptimeVQE)
    errtvqe[k] = std(tmptimeVQE)
end

#---- Save data ----
#select the name
date = now()
fname = string("/home/rezad/Documents/Arbeit/EPFL/Labal/PDM/scripts/Julia /data/PTvsVQE_" , "Nmax=", string(Nmax),"_", "Nsteps=", string(Nsteps),"_" , string(date) , ".csv")
#creating the file
touch(fname)
#opening the file in write mode
data = open(fname, "w")
#concatenating the data into a matrix
savedata = hcat(Nlist,avrgEpt, errEpt, avrgtpt, errtpt, memPT, avrgEvqe, errEvqe, avrgtvqe, errtvqe, memVQE)
#writing it to the file
writedlm(fname, savedata, ",")
close(data)