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

# bim = @elapsed bam = minimum([8, 2,3, 4, 5, 2, 3, 6,1])
# println(bim, bam)


N = Int8(8)
R = Int64(8)
J = 1.
h = 0.5
L = Int8(1)

# points = zeros(Int64, 5)
avrgN = 10

Nmctut = [i for i in 2:10]
Nsteps = UInt64(200)
#Temperature range
bestE = zeros(Float64, length(Nmctut))
bestt = zeros(Float64, length(Nmctut))
for k in 1:length(Nmctut)
    avrgE = zeros(Float64,avrgN)
    avrgt = zeros(Float64,avrgN)
    for l in 1:avrgN
        Tin = logrange(1e-6, 5e-1, Int64(Nmctut[k]))
        Tfin = logrange(1e-8, 1e-2, Int64(Nmctut[k]))
        Ptp = ParTemp(N, R, L, J, h, Int8(Nmctut[k]), Nsteps, Tin, Tfin)

        avrgt[l] = @elapsed avrgE[l] = minimum(run(Ptp))
    end
    bestE[k] = mean(avrgE)
    bestt[k] = mean(avrgt)
end

println(bestE)
println(bestt)

#Temperature range
Tin = 0.001 
Tfin = 0.00001
Emh = zeros(Float64, length(Nmctut))
tmh = zeros(Float64, length(Nmctut))
for i in 1:avrgN
    Vmh = VarMH(N, R, L, J, h )
    initVMHtoHF(Vmh)
    # println(norm(Vmh.HPm.Phaselist))
    # println(Vmh.E)
    # println(Vmh.HPm.H)
    # println(Vmh.HPm.P)

    tmh[i] = @elapsed run(Vmh, Nsteps, Tin, Tfin)
    Emh[i] = Vmh.E
    # println(Vmh.E)
    # if (Vmh.E < -9.)
    #     print("badabim badaboom")
    #     # println(Vmh.E)
    #     println(Vmh.Anglelist)
    #     println(Vmh.HPm.CHlist)
    #     break
    # end
end
println(mean(Emh))
println(mean(tmh))
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