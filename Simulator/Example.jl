using Base: String, Bool, Complex
using LinearAlgebra

include("CHform.jl")
include("GSmat.jl")
include("MetroHast.jl")
include("ParTemp.jl")
include("extrafunctions.jl")

#This example script simply runs a Metropolis instance and a Parallel tempering instance and prints their best energy.
#More examples can be found in the "Extra examples" folder


N = Int8(8) #Number of qubits

R = Int64(8) #Number of stabilizer states to approximate the Ground state

#Parameters of the model
J = 1.
h = 0.5

#Number of layers of the Ansatz (parametric circuit)
L = Int8(1)

#Number of steps in the optimisation process
Nsteps = UInt64(20)


#---------------Metropolis-Hastings--------------------
#Initial temperature
Tin = 0.001 
#Final temperature
Tfin = 0.00001

#initialise the algorithm
Vmh = VarMH(N, R, L, J, h )
#start at the Hartree-Fock state
initVMHtoHF(Vmh)

#run the algorithm
run(Vmh, Nsteps, Tin, Tfin)

#print the energy
println("Energy (Metropolis): ",Vmh.E)

#---------------Parallel-tempering--------------------
#Number of Markov chains
Nmc = 5

#List of initial temperatures
TinList = logrange(1e-5, 1e-2, Nmc)
#List of final temperatures
TfinList = logrange(1e-8, 1e-4, Nmc)

#initialise the method, every decomposition is initialised to the HF state by default
Ptp = ParTemp(N, R, L, J, h, Int8(Nmc), Nsteps, TinList, TfinList)

#run the algorithm
Energies = run(Ptp)

#print the results
println("Energies (Parallel tempering): ", Energies)

