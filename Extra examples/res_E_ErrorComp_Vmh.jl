using Base: String, Bool, Complex
# using Core: Matrix, Vector
using LinearAlgebra
using Dates
using DelimitedFiles

include("CHform.jl")
include("GSmat.jl")
include("MetroHast.jl")

function logrange(x1::Float64, x2::Float64, n::Int64) 
    return [10^y for y in range(log10(x1), log10(x2), length=n)]::Vector{Float64}
end

#This example shows the implementation of the initialisation and first iteration of the ground state optimisation algorithm.
#Every step is explained in order to give a clear explanation of the algorithm.

# dividers = np.array([4,3,2,1,0.5]) #newdividers
# dividers = np.sort(np.append(np.logspace(np.log10(1/8), np.log10(8), num = 10),1)) #moredividers
dividers = logrange(1. /8., 8., 10)
dividers = sort(append!(dividers, [1.]))

N = Int8(8)
R = Int64(16)
J = 1.
h = J./dividers
h = [1.25992105, 1.25992105, 0.5, 0.5, 0.5]
L = Int8(1)

#Temperature range
Tin = 0.001 
Tfin = 0.00001

# Nsteps = int(np.floor(R*beta*N*np.log2(N)))
Nsteps = UInt64(10000)

# E_end = []
E_end = zeros(Float64, length(h))
# counter = 0
# minVal = 0
for k in 1:length(h)
    Vmh = VarMH(N, R, L, J, h[k] )
    initVMHtoHF(Vmh)
    #--- RUN --- 
    println("tut ", h[k])
    Elog, Tlog = run(Vmh, Nsteps, Tin, Tfin)
    # append!(E_end ,[Elog[end]])
    E_end[k] = Elog[end]
    println("MinE = ", E_end[k])
    println("------------")
end

#---- Save data ----
#select the name
date = now()
fname = string("/home/rezad/Documents/Arbeit/EPFL/Labal/PDM/scripts/Julia /data/energy_" , "R=", string(R),"_", "Nsteps=", string(Nsteps),"_" , string(date) , ".csv")
#creating the file
touch(fname)
#opening the file in write mode
data = open(fname, "w")
#concatenating the data into a matrix
savedata = hcat(E_end, h)
#writing it to the file
writedlm(fname, savedata, ",")
close(data)


#améliorations: parallel tempering, quantum natural gradient, todo: optimiser calc gamma et surtout calcul H et P
#définir notion de proche dans hilbert