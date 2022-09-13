using Base: String, Bool, Complex
# using Core: Matrix, Vector
using LinearAlgebra
using Dates
using DelimitedFiles
using Statistics

include("CHform.jl")
include("GSmat.jl")
include("MetroHast.jl")
include("ParTemp.jl")

function logrange(x1::Float64, x2::Float64, n::Int64) 
    return [10^y for y in range(log10(x1), log10(x2), length=n)]::Vector{Float64}
end

#This example shows the implementation of the initialisation and first iteration of the ground state optimisation algorithm.
#Every step is explained in order to give a clear explanation of the algorithm.

# dividers = np.array([4,3,2,1,0.5]) #newdividers
# dividers = np.sort(np.append(np.logspace(np.log10(1/8), np.log10(8), num = 10),1)) #moredividers
dividers = logrange(1. /2., 2., 4)
dividers = sort(append!(dividers, [1., 1/3, 3]))

# Rmax = 20
AmountR = 3
# Rlist = unique!(trunc.(Int8, logrange(2, Float64(Rmax), AmountR))) #convert the logrange to ints (without duplicates)
# Rlist = sort(append!([4*i for i in 1:4],[2]))
Rlist = [16]
AmountR = length(Rlist)
println(Rlist)

N = Int8(8)
# N = Int64(ceil(mean(Rlist)))
# println(R)
J = 1.
h = J./dividers
h = [0.7937005259840997, 1.]
L = Int8(1)

#Temperature range
Tin = 0.001 
Tfin = 0.00001


#parameters for the parallel tempering
Nmc = Int8(10)
Nstepsin = UInt64(1000)
Tin = logrange(1e-5, 1e-1, Int64(Nmc))
Tfin = logrange(1e-7, 1e-3, Int64(Nmc))

# E_end = []
E_end = zeros(Float64, length(h), AmountR)
# counter = 0
# minVal = 0
for l in 1:AmountR
    for k in 1:length(h)
        # if h[k] > 1.
        #     Nsteps = UInt64(ceil(3*h[k]*Nstepsin))
        # else
        #     Nsteps = UInt64(Nstepsin)
        # end
        Nsteps = UInt64(ceil((2*h[k]+1)*Nstepsin))
        # Vmh = VarMH(N, R, L, J, h[k] )
        # initVMHtoHF(Vmh)
        Pt = ParTemp(N, Rlist[l], L, J, h[k], Nmc, Nsteps, Tin, Tfin)
        #--- RUN --- 
        println("R = ", Rlist[l]," h = ", h[k])
        Elog = run(Pt)
        # append!(E_end ,[Elog[end]])
        E_end[k, l], minIdx = findmin(Elog)
        println("MinE = ", E_end[k, l])
        println("------------")
        println(Pt.MHlist[minIdx].Anglelist)
    end
end
#concatenating data 
tmpdata = vcat(E_end, transpose(Rlist))
# println(tmpdata)
append!(h, [N])
# println(h)
savedata = hcat(h, tmpdata)
#---- Save data ----
#select the name
date = now()
fname = string("/home/rezad/Documents/Arbeit/EPFL/Labal/PDM/scripts/Julia /data/energy_" , "N=", string(N),"_", "Nstepsin=", string(Nstepsin),"_", "AmountR=", string(AmountR),"_" , string(date) , ".csv")
#creating the file
touch(fname)
#opening the file in write mode
data = open(fname, "w")
# #concatenating the data into a matrix
# savedata = hcat(E_end, h)
#writing it to the file
writedlm(fname, savedata, ",")
close(data)


#améliorations: parallel tempering, quantum natural gradient, todo: optimiser calc gamma et surtout calcul H et P
#définir notion de proche dans hilbert