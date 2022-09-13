using Base: String, Bool, Complex
# using Core: Matrix, Vector
using LinearAlgebra
using Dates
using DelimitedFiles
using Statistics

include("CHform.jl")
include("GSmat.jl")
# include("MetroHast.jl")
# include("ParTemp.jl")
include("extrafunctions.jl")
include("QNG.jl")

# bim = @elapsed bam = minimum([8, 2,3, 4, 5, 2, 3, 6,1])
# println(bim, bam)


N = Int8(2)
R = Int64(1)
J = 1.
h = 0.5
L = Int8(1)
Ansatz = "Y"
# Ansatz = "XZ"

#Temperature range
Tin = 0.001 
Tfin = 0.00001
Nsteps = 2
etalist = 1 ./(Tin .+ (Tfin - Tin)*(log.(LinRange(2., 10., Nsteps)) .-log(2))/(log(10)-log(2)))

qng = QNG(N, R, L, J, h, Ansatz)
# initQNGtoHF(qng)
# qng.Anglelist[1,1,:] = mod.(rand(Int8, 2*N), 8)
qng.Anglelist[1,1,:] = [6;;; 0;;; 3;;; 0]
# rankd = 0
# while rankd != R
#     qng.Anglelist[1,2,:] = mod.(rand(Int8, 2*N), 8)
#     # qng.Anglelist[1,2,1] = 2
#     global circ2 = genVarCircY(qng, Int8(1), 2)
#     update_circ(qng.HPm.CHlist[2], circ2)
#     # updateHP(qng.HPm, 2, circ2)
#     compHP(qng.HPm)
#     global rankd = rank(qng.HPm.P)
#     println(rankd)
# end
circ1 = genVarCirc(qng, Int8(1), 1, Ansatz)
update_circ(qng.HPm.CHlist[1], circ1)
compHP(qng.HPm)
println(qng.Anglelist)
# println(circ2)
println(qng.HPm.P)
println(rank(qng.HPm.P))
E, alphas = eigen(qng.HPm.H, qng.HPm.P)
println(E)
qng.E, idxMin = findmin(E)
qng.HPm.Phaselist = alphas[:, idxMin]
# println(norm(qng.HPm.Phaselist))
println(qng.E)
# println(qng.HPm.H)
# println(qng.HPm.P)

run(qng, etalist)
println(qng.E)
# println(Vmh.E)
# if (Vmh.E < -9.)
#     print("badabim badaboom")
#     # println(Vmh.E)
#     println(Vmh.Anglelist)
#     println(Vmh.HPm.CHlist)
#     break
# end
# println(Vmh.Anglelist)

# println(Es)
# println(minimum(Es))

# minimumE = Ptp.MHlist[1].E
# for k in Ptp.MHlist
#     print(k.E, " ")
#     if (k.E < minimumE)
#         global minimumE = k.E
#     end
# end
# println(" ")
# println("min = ", minimumE)

# # Computing the analytical energy to check that the result is not abhorrent
# m = 1:Ptp.MHlist[1].HPm.CHlist[1].N
# m = m .- floor(Ptp.MHlist[1].HPm.CHlist[1].N/2) #element-wise substraction
# #computation of k
# k = 2*pi .*m ./Ptp.MHlist[1].HPm.CHlist[1].N
# #Method 2
# E0 = -1 * sum(sqrt.((-Ptp.MHlist[1].HPm.h/Ptp.MHlist[1].HPm.J .- cos.(k)).^2 + sin.(k).^2))

# println("E0 = ", E0)

# for k in Ptp.MHlist
#     println(k.Anglelist)
# end