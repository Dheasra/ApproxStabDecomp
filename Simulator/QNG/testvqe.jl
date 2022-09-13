using Base: String, Bool, Complex
# using Core: Matrix, Vector
using LinearAlgebra
# using Dates
# using DelimitedFiles
using ProgressMeter

include("CHform.jl")
include("GSmat.jl")
include("MetroHast.jl")
include("ParTemp.jl")
include("extrafunctions.jl")
include("vqe.jl")

#General parameters
N = Int8(2)
J = 1.
h = 0.5
L = Int8(1)
Niter = 2
etaIn = 5e-2
etaFin = 1e-4 #step size
etalist = (etaIn .+ (etaFin - etaIn)*(log.(LinRange(2., 10., Niter)) .-log(2))/(log(10)-log(2)))

# Ham, circuit, params = initvqe(N, L, J, h)
Ham = hamiltonian(N, J, h)
params = [6, 3, 0, 0]
circuit = ansatz_alt(N, L, params)
@showprogress for i = 1:Niter
    # global Energy, params = vqestep(Ham, circuit, params, etalist[i], N)
    global Energy, params = vqestep(Ham, circuit, params, 1., N)
    # println(params)
    if abs(Energy - (-21.307328078274512)) < 1e-14
        break
        println(i)
    end
end

println("Energy = ", Energy)

# n = Int8(3); depth = Int8(1);
# J = -1; B = -1;
# # n_params = depth*((n-1) + (n-2) + n)
# n_params = depth*2*n
# # print("nparams =", n_params)
# params = zeros(n_params)
# circuit = ansatz(n, depth, params)
# # circuit = dispatch!(variational_circuit(n, depth),:random);

# h = hamiltonian(n, J, B);
# # h = heisenberg(n)

# energies = []

# opt_steps = 10
# eta = 1e-3
# for i = 1:opt_steps
#     # grad = faithful_grad(h, zero_state(n)=>circuit; nshots=100)
#     _, grad = expect'(h, zero_state(Int64(n))=>circuit)
#     energy = expect(h, zero_state(Int64(n))=>circuit)
#     push!(energies, energy)
#     global params -= eta*grad
#     dispatch!(circuit, params)
#     # global circuit = alternate_timedep_ansatz(n, depth, params - eta*grad)
#     # println("Step $i, energy = $(energy)")
# end

# println("Energy = ", energies)
