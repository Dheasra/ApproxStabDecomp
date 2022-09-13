# GENERAL
Ground state optimisation algorithm with approximate state representation using stabilizer state superposition.

This repository contains all the scripts of the final implementation of a ground state optimisation algorithm with approximate state representation.
This implementation serves as the backbone of my [Quentin Pitteloud] Master thesis.
The project was done in conjunction with CERN's Openlab under the direction of Sofia Vallecorsa and EPFL's Laboratory of Theoretical Physics of Nanosystems (LPTN) under the direction of Vincenzo Savona.

This implementation is written in Julia for performance.

The aim and theory behind this implementation is discussed in the "Thesis.pdf" file.


# QUICK START
To get a demo of it running, simply run the "Example.jl" script.
Extra examples of usage are given in the form of scripts used to generate results for my thesis in the "Extra examples" folder.

# PACKAGE REQUIREMENTS
This code uses the following packages:
-Base: String, Bool, Complex
-LinearAlgebra
-Random
-ProgressMeter (not necessary, used in the Metropolis/Parallel tempering "run" functions to keep track of the progress)

Furthermore, this implementation uses multithreading to accelerate the optimisation processes. Do not forget to add the "--threads=auto" argument in your commands!

# QUICK DESCRIPTION - For more info, please refer to the source code where function usage are given as commentaries
-CHform.jl: Quantum simulator of stabilizer states in CH-form representation. Main functions: scalPrdt (Compute scalar products or expectation values of Hamiltonian), update_circ (update the CH-form for a Clifford circuit written as [["gate"],[string(qubit)]])
-GSmat.jl: Implementation of superpositions of CH-forms, along with matrices H and P. Main function: updateHP (Updates the matrices H and P for a Clifford circuit on a chosen element of the decomposition)
-MetroHast.jl: Metropolis-Hastings optimisation algorithm. Main functions: initVMHtoHF (Initialise the decomposition to the state |0^N> and neighbouring states), run (optimise the decomposition)
-ParTemp.jl: Parallel tempering optimisation algorithm. Main function: run (optimise the decomposition)

# CHANGING PHYSICAL MODELS
Currently, only the transverse field Ising model is implemented in this code. It has been hard-coded in the "CHform.jl" script.
To change the physical model, modify the "scalPrdt" function in said script to fit your model.
Note that you'll probably also have to change the initial stabilizer decomposition! This can be done by changing the "initVMHtoHF" function in MetroHast.jl for both the Metropolis and Parallel tempering algorithms.

# Quantum Natural Gradient (QNG)
The QNG has also been implemented. However, there were issues relating to the construction of the quantum gradient tensor (making it singular).
Additionally, the gradient of the energy would always be zero preventing the use of the QNG.
This method is also included here should anyone want to have a look at it themselves.
The implementation consists only of the QNG.jl file. It has been compared to a variational quantum eigensolver implemented in the file vqe.jl. 
An example of its use is included in the file testvqe.jl.

# FINAL WORDS AND ACKNOWLEDGEMENTS
I hope this description gave you a clear idea of what this project is about.
Should you have any question, I can be contacted through Github itself or at quentin.pitteloud@gmail.com.
This was a really interesting project and I am happy that I have had the opportunity to do it.
For that I would express my sincere thanks to Prof. Vincenzo Savona, and to Dr. Sofia Vallecorsa for allowing me to do it in CERN's Openlab.
I am also very grateful to Dr. Fabrizio Minganti, Dr. Riccardo Rota and Dr. Michele Grossi for their help during my semester projects and thesis.
