using Base: String, Bool, Complex
# using Core: Matrix, Vector
using LinearAlgebra

include("CHform.jl")
include("GSmat.jl")
include("MetroHast.jl")

N = Int8(8)
R = Int64(8)
J = 1.
h = 0.5
L = Int8(1)

#Temperature range
Tin = 0.001 
Tfin = 0.00001

# Nsteps = int(np.floor(R*beta*N*np.log2(N)))
Nsteps = UInt64(10000)

# println(digits(6, base = 2, pad = 3))
# println(digits(4, base = 2, pad = 3))

# for i in 1:1
#     Vmh = VarMH(N, R, L, J, h )
#     initVMHtoHF(Vmh)
#     # println(norm(Vmh.HPm.Phaselist))
#     println(Vmh.E)
#     # println(Vmh.HPm.H)
#     # println(Vmh.HPm.P)

#     run(Vmh, Nsteps, Tin, Tfin)
#     println(Vmh.E)
#     if (Vmh.E < -9.)
#         print("badabim badaboom")
#         # println(Vmh.E)
#         println(Vmh.Anglelist)
#         println(Vmh.HPm.CHlist)
#         break
#     end
# end
# # println(Vmh.E)
# # println(Vmh.Anglelist)


# Vmh = VarMH(N, R, L, J, h )
# # Vmh.Anglelist = [7 5 3 1 3 5 2 0;;; 3 4 3 2 3 5 1 3;;; 1 5 3 6 7 3 3 3;;; 1 1 3 6 7 3 7 7;;; 0 6 5 6 1 1 1 7;;; 7 2 1 6 1 1 1 3;;; 0 4 7 6 5 3 7 7;;; 3 2 3 6 1 7 3 7;;; 3 4 2 2 1 2 3 3;;; 3 7 2 5 1 1 3 3;;; 1 5 5 5 0 5 7 1;;; 1 2 5 0 6 5 7 1;;; 5 5 5 7 4 3 3 1;;; 5 0 5 3 4 7 7 1;;; 5 5 1 5 5 7 2 4;;; 1 7 5 1 1 3 4 5]
# Vmh.Anglelist = [2 0 7 0 4 0 5 0;;; 1 5 1 0 0 0 1 0;;; 1 4 4 1 4 0 0 0;;; 6 0 6 3 0 0 0 0;;; 2 0 5 0 0 0 0 0;;; 5 0 3 0 0 0 0 0;;; 6 0 7 0 0 0 0 0;;; 1 0 6 0 0 0 1 0;;; 7 0 3 4 0 0 0 0;;; 5 0 3 0 0 0 0 0;;; 6 4 4 0 4 0 0 0;;; 7 0 1 0 0 0 0 0;;; 4 4 2 0 0 0 0 0;;; 4 2 5 0 0 0 0 0;;; 1 0 0 0 1 0 0 0;;; 5 0 6 0 0 0 0 0]
# # Vmh.Anglelist = [0 3 0 0;;; 6 6 2 0;;; 0 4 4 5;;; 4 6 5 4;;; 4 6 0 4;;; 4 1 1 2;;; 5 3 7 2;;; 4 0 2 3;;; 2 0 2 0;;; 4 5 0 5;;; 1 3 5 1;;; 0 3 4 0;;; 0 6 4 5;;; 7 5 6 6;;; 0 3 0 5;;; 1 3 6 5]
# # Vmh.Anglelist = [5 4 3 4;;; 2 0 6 0;;; 5 4 6 4;;; 6 3 1 0;;; 5 4 0 0;;; 3 2 4 0;;; 4 4 6 4;;; 3 5 3 3;;; 7 0 2 4;;; 7 3 0 5;;; 5 4 7 4;;; 1 3 6 2;;; 3 0 7 0;;; 2 3 7 2;;; 1 0 1 0;;; 7 0 0 4]
# # println(Vmh.Anglelist[1,1,:])
# # Vmh.Anglelist = [4 1 5 1;;; 4 5 0 5;;; 3 4 2 4;;; 2 4 6 5;;; 3 0 0 0;;; 6 0 7 3;;; 3 0 6 0;;; 4 0 3 1;;; 5 0 0 0;;; 3 0 4 0;;; 3 4 2 0;;; 4 0 7 4;;; 5 4 2 4;;; 2 0 2 3;;; 5 4 4 0;;; 5 0 2 0]
# # println(Vmh.Anglelist)
# for l::Int8 in 1:Vmh.L 
#     for r in 1:Vmh.HPm.R
#     # for r in 1:
#         # println(l, r)
#         circ = genVarCirc(Vmh, l, r)
#         # println(circ)
#         updateHP(Vmh.HPm, r,circ, r)
#         # if r == 1
#         #     println(circ[1])
#         #     println(" ")
#         #     println(circ[2])
#         # end
#     end
# end
# println(Vmh.HPm.CHlist[1])
# println(" ")
# println(Vmh.HPm.CHlist[3])
# println(rank(Vmh.HPm.P))
# # println(Vmh.HPm.H)
# # println(" ")
# # println(Vmh.HPm.P)
# # println(" ")
# E, alphas = eigen(Vmh.HPm.H, Vmh.HPm.P) #compute the energy and the coefficients of the decomposition
# println(E)

println(Int64(12^13))

count = 12
r = 13 
julianiquetamère = Int8(1)
while (r > count)
    global julianiquetamère += Int8(1)
    global count += 12^julianiquetamère
    if julianiquetamère >= 2*12
        println(12^julianiquetamère)
        throw(DomainError(julianiquetamère,"Too many stabilizers"))
    end
end
println(julianiquetamère)
#TODO: faire en sorte de stocker les valeurs importantes quand on a un résultat débile pour pouvoir le corriger 

# println(norm(Vmh.HPm.Phaselist))

# a = CHform(Matrix{Int8}(I,N,N), Matrix{Int8}(I,N,N), zeros(Int8,N,N), zeros(Int8,N), zeros(Int8,N), zeros(Int8,N), 1 + 0im, N)
# println(a.F)
# circ = [["H", "CX", "S"],["1", "1 2", "2"]]
# circ = [["X","X"],["1","1"]]
# circa = [["X"],["1"]]
# circb = [["H", "S","H"],["1", "1", "2"]]

# HPm = HPmat(N, R, J, h)
# updateHP(HPm, 1, circa)
# updateHP(HPm, 2, circb)
# println(HPm.H)
# println(HPm.P)
# println(HPm.H[2,1])

# circRztest = [RzGate(Int8(k-1), k) for k::Int8 = 1:8]
# println(circRztest)
# circRxtest = [RxGate(Int8(k-1), k) for k::Int8 = 1:8]
# println(circRxtest)

# a = CHform(N)
# update_circ(a, circa)
# # println(a)

# # # ad = conjugateCHf(a)
# # # println(ad)

# b = CHform(N)
# update_circ(b, circb)
# b = conjugateCHf(b)
# println("=======================================================")

# amp = scalPrdt(b, a, 1., 1.)
# println(amp)

# ad = conjugateCHf(a)
# println(scalPrdt(b, a))

