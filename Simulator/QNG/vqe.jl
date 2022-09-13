# using Yao,YaoExtensions,PyCall,YaoPlots, JSON, Statistics, Yao.AD, Compose
# using Yao,YaoExtensions, Statistics, Yao.AD, Compose
using Base: String, Bool, Complex
# using Core: Matrix, Vector
using LinearAlgebra
using Yao, YaoAPI, Statistics, Compose
# using Yao
# using QuAlgorithmZoo: Adam, update!

# include all the useful functions for VpVQD
# include("pvqd_functions.jl")
# import PyPlot; const plt = PyPlots

function ansatz(n::Int8, depth::Int8, params::Vector{Float64}) #Based on the implementation by D. "Delphouts" Martres
    count = 1
	n = Int64(n)
	circ = chain(n, put(n, i=>I2) for i in 1:n)
	for d in 1:depth
		for i in 1:n
			#push!(circ, Rx(params[count]))
            # push!(circ, chain(n, put(n, i=>H)))
            # println(typeof(out))
			push!(circ, chain(n, put(n, i=>Rx(params[count]))))
			count = count +1
            # push!(circ,chain(n,put(n, i=>Rz(params[count]))))	
			# count = count+1
		end
        for i in 1:n
			push!(circ,chain(n,put(n, i=>Rz(params[count]))))	
			count = count+1
		end

        for i in 1:n
            push!(circ,chain(n,cnot(i,mod(i,n)+1)))
        end
	end
    # print("count= ",count)
	return circ
end

function ansatz_alt(n, depth, params)
    n = Int64(n)
    count = 1
	# circ = chain(n, put(i=> H) for i in 1:n)
    circ = chain(n, put(i=> I2) for i in 1:n)
	
	for d in 1:depth
		for i in 1:n
			#push!(circ, Rx(params[count]))
			push!(circ,chain(n,put(n, i=>Ry(params[count]))))
			count = count +1
		end

        for i in 1:n-1
            push!(circ,chain(n,cnot(i,i+1)))
        end

		for i in 1:n
			push!(circ,chain(n,put(n, i=>Ry(params[count]))))	
			count = count+1
		end
	end
    # print("count= ",count)
	return circ
end

function vqe(Nsteps::UInt64, N::Int8, L::Int8, J::Float64, h::Float64)
    # #initializing relevant dada
    # n_params = L*2*N #number of parameters for the circuit
    # params = zeros(n_params) #Vector of parameters
    # #Constructing the parameterized circuit
    # circuit = ansatz(N, L, params)
    # #Constructing the cyclic Transverse Ising Model
    # H = hamiltonian(N, J, h);
    Ham, circuit, params, n_params = initvqe(N, L, J, h)
    #storing the energy (output)
    energies = []
    eta = 1e-3 #step size
    for i = 1:Nsteps
        # grad = faithful_grad(h, zero_state(n)=>circuit; nshots=100)

        # _, grad = expect'(h, zero_state(n)=>circuit)
        # energy = expect(h, zero_state(n)=>circuit)
        # global params -= eta*grad
        # dispatch!(circuit, params)
        global energy, params = vqestep(Ham, circuit, params, eta, n)
        push!(energies, energy)

        # global circuit = alternate_timedep_ansatz(n, depth, params - eta*grad)
        # println("Step $i, energy = $(energy)")
    end
    return minimum(energies), energies
end

function initvqe(N::Int8, L::Int8, J::Float64, h::Float64)
    #initializing relevant dada
    n_params = L*2*N #number of parameters for the circuit
    params = zeros(n_params) #Vector of parameters
    #Constructing the parameterized circuit
    # circuit = ansatz(N, L, params)
    circuit = ansatz_alt(N, L, params)
    # println("circ = ", circuit)
    #Constructing the cyclic Transverse Ising Model
    Ham = hamiltonian(N, J, h);
    return Ham, circuit, params
    # #storing the energy (output)
    # energies = []
end

function vqestep(Ham, circuit, params, eta::Float64, n::Int8)
    n = Int64(n)
    bimgrad, grad = expect'(Ham, zero_state(n)=>circuit)
    # println("bimgrad = ",bimgrad)
    println("grad = ",grad)
    energy = expect(Ham, zero_state(n)=>circuit)
    params -= eta*grad
    dispatch!(circuit, params)
    return energy, params
end


#---- Functions created by either Delphouts or S. Barisson ---------------
function hamiltonian(n,J,B)
	op_list = []

	for i in 1:n
		push!(op_list,-J*spin_zz(n,i))
	end
	for i in 1:n
		push!(op_list,-B*spin_x(n,i))
	end
	
	obs = sum(op_list)

	return obs
end

function spin_zz(n,site)
	op_list = []

	for i in 1:n
		if i==site || i==mod(site+1,n)
			push!(op_list,Z)
		else 
			push!(op_list,I2)
		end
	end
	
	obs = kron(op_list...)

	return obs
end

function spin_x(n,site)
	op_list = []

	for i in 1:n
		if i==site 
			push!(op_list,X)
		else 
			push!(op_list,I2)
		end
	end

	obs = kron(op_list...)

	return obs
end	

# n = 3; depth = 1;
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
#     _, grad = expect'(h, zero_state(n)=>circuit)
#     energy = expect(h, zero_state(n)=>circuit)
#     push!(energies, energy)
#     global params -= eta*grad
#     dispatch!(circuit, params)
#     # global circuit = alternate_timedep_ansatz(n, depth, params - eta*grad)
#     # println("Step $i, energy = $(energy)")
# end

# # gatecount(circuit)

# # close("all")
# # exact    = last(JSON.parse(open("data/exact/T2.0_dt0.05.dat","r"))["energies"])
# # plt.hlines(exact, 1, opt_steps, linestyle="dashed", color="black")
# # plt.plot(1:opt_steps, energies)
# # plt.xlabel("step")
# # plt.ylabel("energy")
# # plt.gcf()