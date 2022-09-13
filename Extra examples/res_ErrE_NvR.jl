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
include("extrafunctions.jl")

#List of #qubits
Nmax = 10
AmountN = 2
Nlist = unique!(trunc.(Int8, logrange(9., Float64(Nmax), AmountN))) #convert the logrange to ints (without duplicates)
AmountN = length(Nlist)
println(Nlist)
#List of #stabilizers
Rlist = unique!(trunc.(Int64, logrange(9., Float64(Nmax), AmountN)))
# Rmax = 20
# AmountR = 3
# Nlist = unique!(trunc.(Int8, logrange(9., Float64(Nmax), AmountN))) #convert the logrange to ints (without duplicates)
AmountR = length(Rlist)
# println(Nlist)

# N = Int8(8)
# R = Int64(ceil(mean(Nlist)))
# println(R)
J = 1.
h = 0.5
L = Int8(1)

#parameters for the parallel tempering
NMC = Int8(7)
Nstepsin = UInt64(2)
Tin = logrange(1e-5, 5e-1, Int64(NMC))
Tfin = logrange(1e-7, 1e-2, Int64(NMC))

# E_end = []
E_end = zeros(Float64, AmountR, AmountN)
for l in 1:AmountN
    for o in 1:AmountR
        Nsteps = UInt64(ceil((Rlist[o]+Nlist[l])*Nstepsin))
        # Vmh = VarMH(N, R, L, J, h )
        # initVMHtoHF(Vmh)
        Pt = ParTemp(Nlist[l], Rlist[o], L, J, h, NMC, Nsteps, Tin, Tfin)
        #--- RUN --- 
        println("N = ", Nlist[l]," R = ", Rlist[o])
        Elog = run(Pt)
        # append!(E_end ,[Elog[end]])
        E_end[l, o] = minimum(Elog)
        println("MinE = ", E_end[l, o])
        println("------------")
    end
end
#concatenating data 
tmpdata = vcat(E_end, transpose(Nlist))
# println(tmpdata)
append!(Rlist, [J/h])
# println(h)
savedata = hcat(Rlist, tmpdata)
#---- Save data ----
#select the name
# date = now()
fname = string("/home/rezad/Documents/Arbeit/EPFL/Labal/PDM/scripts/Julia /data/errE_" , "NMC=", string(NMC),"_Nstepsin=", string(Nstepsin),"_", "J=", string(J),"_" , "h=", string(h), ".csv")
#creating the file
touch(fname)
#opening the file in write mode
data = open(fname, "w")
# #concatenating the data into a matrix
# savedata = hcat(E_end, h)
#writing it to the file
writedlm(fname, savedata, ",")
close(data)