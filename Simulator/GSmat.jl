using Base: String, Bool, Complex
# using Core: Matrix, Vector
using LinearAlgebra

include("CHform.jl")#TODO: Checker que ce soit la bonne syntaxe

mutable struct HPmat
    H::Matrix{Complex{Float64}}# H_ij = <phi_i|H|phi_j>
    P::Matrix{Complex{Float64}}# P_ij = <phi_i|phi_j>
    J::Float64 #Bond factor (Ising Model)
    h::Float64 #Transverse field intensity (Transverse field Ising Model)
    CHlist::Array{CHform} #List (sum) of the CH-forms of the stabilizer states approximating the ground state
    Phaselist::Array{Complex{Float64}} #Phase coefficients of each of the stabilizer states
    R::Int64 #Number of elements (stab. states) in the decomposition
    #TODO: écrire un constructeur 
    HPmat(N, R, J, h) = new(zeros(Complex,R,R), zeros(Complex,R,R), J::Float64, h::Float64,[CHform(N) for k::Int8 in 1:R], [1/sqrt(R) for k::Int8 in 1:R],R::Int64) #N is the number of qubits
    # CHform(N) = new(Matrix{Int8}(I,N,N), Matrix{Int8}(I,N,N), zeros(Int8,N,N), zeros(Int8,N), zeros(Int8,N), zeros(Int8,N), 1 + 0im, N)
end

#------------------ Computation and update of the matrices H and P------------------------------------------------------------------------

function compHP(HPm::HPmat) #computes the elements of P,H
    for k::Int64 in 1:HPm.R
        update(HPm, k,k)
    end
end

function updateHP(HPm::HPmat,TargetElem::Int64, circ::Vector{Vector{String}}, MaxElem::Int64 = 0) #updates the matrices H and P when the elem^th element of Omega is modified by circ
    update_circ(HPm.CHlist[TargetElem], circ)
    # @views minimises the allocation overhead by not copying splices of arrays into new memory locations
    @views update(HPm, TargetElem, MaxElem) #We must change all the elements of a column, as such the MaxElem is implicitly HPm.Omega.Nsum
end

function update(HPm::HPmat,TargetElem::Int64, MaxElem::Int64 = 0) #updates the matrices for every element corresponding to TargetElem
    #handling the default value as "all column"
    if MaxElem == 0
        MaxElem = HPm.R 
        # MaxElem = 1
    end
    #We conjugate the second vector to obtain its bra form
    tmp = conjugateCHf(HPm.CHlist[TargetElem])
    #compute/update the upper triangular parts of H and P
    Threads.@threads for k in 1:MaxElem
        # println(k)
        HPm.P[TargetElem, k] = scalPrdt(tmp, HPm.CHlist[k]) 
        HPm.H[TargetElem, k] = scalPrdt(tmp, HPm.CHlist[k], HPm.J, HPm.h) 
        # println(TargetElem, " ",k, HPm.P)
        if (k != TargetElem) #adding the off-diagonal lower triangular elements
            HPm.P[k, TargetElem] = conj(HPm.P[TargetElem, k])
            HPm.H[k, TargetElem] = conj(HPm.H[TargetElem, k])
        end
    end
    # #Now compute the lower triangular part of the matrices using the fact that they are Hermitian
    # HPm.P = Hermitian(HPm.P, :L) #by default, the upper triangle part is used as reference and conjugated to the lower triangular part (otherwise use :L as second argument)
    # HPm.H = Hermitian(HPm.H, :L)
end

#------------------ UTILIY FUNCTIONS ------------------------------------------------------------------------------------------------------

function clearElem(HPm::HPmat, TargetStabList::Vector{Int64}) #Resets all stab. states whose indices are in TargetStabList
    for r in TargetStabList
        HPm.CHlist[r] = CHform(HPm.CHlist[1].N)
    end
end