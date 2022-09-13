using Base: String, Bool, Complex
# using Core: Matrix, Vector
using LinearAlgebra
using Yao
# using PyCall

include("CHform.jl")
include("GSmat.jl")
include("MetroHast.jl")

# plt = pyimport("matplotlib.pyplot")

function logrange(x1::Float64, x2::Float64, n::Int64)  #Generates a list of n numbers distributed logarithmically between x1 and x2 
    return [10^y for y in range(log10(x1), log10(x2), length=n)]::Vector{Float64}
end

function translate(circ::Vector{Vector{String}}, N::Int8) #Translate a circuit in list form into a Yao circuit
    N = Int64(N)
    out = chain(N, put(N, i=>I2) for i in 1:N)
    for d in 1:length(circ[1])
        qb = circ[2][d]
        gate = circ[1][d]
        # if circ[1][d] == "H"
            # push!(circ, cnot(c, t))
            # push!(circ, chain(N, put(q => H)))
        if (gate == "S")
            q = parse(Int64, qb)
            # push!(out, chain(N, put(N, q=>S))) #bloody language with too many chromosomes can't recognize the S gate even though it is defined
            push!(out, chain(N, put(N, q=>T)))
            push!(out, chain(N, put(N, q=>T)))
        elseif (gate == "H")
            q = parse(Int64, qb)
            push!(out, chain(N, put(N, q=>H)))
        elseif (gate == "Z")
            q = parse(Int64, qb)
            push!(out, chain(N, put(q=>Z)))
        elseif (gate == "X")
            q = parse(Int64, qb)
            push!(out, chain(N, put(N, q=>X)))
        elseif (gate == "CZ")
            qbidx = findfirst(" ", qb)[1]
            q = parse(Int64,qb[1:qbidx])
            r = parse(Int64,qb[qbidx+1:end])
            push!(out, cz(q, Int64(r)))
        elseif (gate == "CX")
            qbidx = findfirst(" ", qb)[1]
            q = parse(Int64,qb[1:qbidx])
            r = parse(Int64,qb[qbidx+1:end])
            push!(out, cnot(q, Int64(r)))
        # elseif (gate == "Id") #identity gate multiplied by some factor. The qubit is not required and is instead replaced by the factor 
        #     factor = parse(Complex{Float64},qb)
        #     CHf.w *= factor
        end   
    end
    return out
end

function randomCirc(N::Int8, Ngates::Int64) #Generates a random Clifford circuit in list form
    circ::Vector{Vector{String}} = [[],[]]
    gatesSel::Vector{String} = ["H", "S", "X", "Z", "CX", "CZ"]
    for k in 1:Ngates
        gate::String = rand(gatesSel)
        if gate == "CX" || gate == "CZ"
            q = rand(1:N)
            r = rand(1:N)
            while q == r 
                r = rand(1:N) 
            end
            append!(circ[1], [gate])
            append!(circ[2], [string(string(q), " ", string(r))])
        else
            q = rand(1:N)
            append!(circ[1], [gate])
            append!(circ[2], [string(q)])
        end
    end
    return circ
end

function spH(st2:: ArrayReg{2, ComplexF64, Matrix{ComplexF64}}, st1:: ArrayReg{2, ComplexF64, Matrix{ComplexF64}}, N::Int8, J::Float64, h::Float64) #Computes the expectation value of the Transverse Field Ising Hamiltonian in Yao
    # print("tut")
    st2vec = state(st2)
    amp = 0.
    for i in 1:N 
        tmp = deepcopy(st1)
        zzchain = translate([["Z","Z"],[string(i), string(mod(i,N)+1)]], N)
        apply!(tmp, zzchain)
        amp += -J* dot(st2vec, state(tmp))
        tmp = deepcopy(st1)
        xchain = translate([["X"],[string(i)]], N)
        apply!(tmp, xchain)
        amp += -h* dot(st2vec, state(tmp))
    end
    return amp
end