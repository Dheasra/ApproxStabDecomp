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
h = 1.
L = Int8(1)

avrg = 2

Nmax = 12
AmountN = 6
Nlist = unique!(trunc.(Int8, logrange(4., Float64(Nmax), AmountN))) #convert the logrange to ints (without duplicates)
AmountN = length(Nlist)
# N = Int8(8)
R = 4

#parameters for the parallel tempering
Nmc = Int8(7)
Nsteps = UInt64(4000)
Tin = logrange(1e-5, 1e-1, Int64(Nmc))
Tfin = logrange(1e-7, 1e-3, Int64(Nmc))

#parameters for the Metropolis algo
Niter = Int64(Nmc * Nsteps)

#Stabilizers
#Energy averages for each N
avrgEstab = zeros(Float64, AmountN)
errEstab = zeros(Float64, AmountN)
#Compbasis
#Energy averages for each N
avrgEcb = zeros(Float64, AmountN)
errEcb = zeros(Float64, AmountN)
#run
for k in 1:AmountN
    # R = Int64(Nlist[k])
    # R = Int64(ceil(sqrt(Nlist[k])))
    #temp values for energy average
    Estab = zeros(Float64, avrg)
    Ecb = zeros(Float64, avrg)

    #computing the analytical energy to check the results
    m = 1:Nlist[k]
    m = m .- floor(Nlist[k]/2) #element-wise substraction
    #computation of k
    k0 = 2*pi .*m ./Nlist[k]
    #Method 2
    E0 = -1 * sum(sqrt.((-h/J .- cos.(k0)).^2 + sin.(k0).^2))
    #Stabilizer states
    for l in 1:avrg
        Ptp = ParTemp(Nlist[k], R, L, J, h, Nmc, Nsteps, Tin, Tfin)
        Estab[l] = minimum(run(Ptp))
    end
    #storing the variables
    avrgEstab[k] = mean(Estab)
    errEstab[k] = std(Estab)

    for l in 1:avrg
        Ptp = ParTemp(Nlist[k], R, L, J, h, Nmc, Nsteps, Tin, Tfin)
        Ecb[l] = minimum(runCompBasis(Ptp))
    end
    #storing the variables
    avrgEcb[k] = mean(Ecb)
    errEcb[k] = std(Ecb)
end

#---- Save data ----
#select the name
date = now()
fname = string("/home/rezad/Documents/Arbeit/EPFL/Labal/PDM/scripts/Julia /data/CBvStab_" ,"R=", string(R),"_", "Nmc=", string(Nmc),"_", "Nsteps=", string(Nsteps),"_","J=", string(J),"_","h=", string(h),"_" , string(date) , ".csv")
#creating the file
touch(fname)
#opening the file in write mode
data = open(fname, "w")
#concatenating the data into a matrix
savedata = hcat(Nlist,avrgEstab, errEstab, avrgEcb, errEcb)
#writing it to the file
writedlm(fname, savedata, ",")
close(data)