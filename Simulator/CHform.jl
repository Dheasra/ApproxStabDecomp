using Base: String, Bool, Complex
# using Core: Matrix, Vector
using LinearAlgebra
# using Threads

#Code taken from https://stackoverflow.com/questions/64909449/does-for-struct-check-recursively-in-julia-it-seems-not
abstract type Comparable end
import Base.==
==(a::T, b::T) where T <: Comparable =
    getfield.(Ref(a),fieldnames(T)) == getfield.(Ref(b),fieldnames(T))

"CHform of a quantum state"
mutable struct CHform <: Comparable
    F::Matrix{Int8}
    G::Matrix{Int8}
    M::Matrix{Int8}
    y::Vector{Int8}
    v::Vector{Int8}
    s::Vector{Int8}
    w::Complex
    N::Int8
    CHform(N) = new(Matrix{Int8}(I,N,N), Matrix{Int8}(I,N,N), zeros(Int8,N,N), zeros(Int8,N), zeros(Int8,N), zeros(Int8,N), 1 + 0im, N)
end

# --------------------------QUANTUM CIRCUIT EVOLUTION FUNCTIONS---------------------------------------------------------
function update_circ(CHf::CHform, circ::Vector{Vector{String}})
    gates = circ[1]
    qbs = circ[2]
    g = length(gates) #length returns the product of all dimensions
    for i in 1:g
        if gates[i] == "H"
            @views updateH(CHf, qbs[i])
        else
            @views updateCL(CHf, gates[i], qbs[i])
        end
    end
end

function updateCL(CHf::CHform, gate::String, qb::String)
    if (gate == "S")
        q = parse(Int8, qb)
        CHf.M[q,:] .= mod.(CHf.M[q,:] +  CHf.G[q,:],2)
        CHf.y[q] = mod(CHf.y[q] - 1,4)
    elseif (gate == "Z")
        q = parse(Int8, qb)
        CHf.y[q] = mod(CHf.y[q] - 2,4)
    elseif (gate == "X")
        updateH(CHf, qb)
        updateCL(CHf, "Z", qb)
        updateH(CHf, qb)
    elseif (gate == "CZ")
        qbidx = findfirst(" ", qb)[1]
        q = parse(Int8,qb[1:qbidx])
        r = parse(Int8,qb[qbidx+1:end])
        CHf.M[q,:] .= mod.(CHf.M[q,:] + CHf.G[r,:],2)
        CHf.M[r,:] .= mod.(CHf.M[r,:] + CHf.G[q,:],2)
    elseif (gate == "CX")
        qbidx = findfirst(" ", qb)[1]
        q = parse(Int8,qb[1:qbidx])
        r = parse(Int8,qb[qbidx+1:end])
        # #IMPORTANT NOTE: This copy must be done BEFORE the update of M and F otherwise it will not work. 
        #computation of (MF^T)[q][r]
        MFt = dot(deepcopy(CHf.M[q,:]),deepcopy(CHf.F[r,:]))
        # MFt .= transpose(deepcopy(CHf.M[q,:])) * deepcopy(CHf.F[r,:])
        #updating the matrices
        CHf.G[r,:] .= mod.(CHf.G[r,:] + CHf.G[q,:],2)
        CHf.F[q,:] .= mod.(CHf.F[q,:] + CHf.F[r,:],2)
        CHf.M[q,:] .= mod.(CHf.M[q,:] + CHf.M[r,:],2)     
        CHf.y[q] = mod(CHf.y[q] + CHf.y[r] + 2*MFt, 4)
    elseif (gate == "Id") #identity gate multiplied by some factor. The qubit is not required and is instead replaced by the factor 
        factor = parse(Complex{Float64},qb)
        CHf.w *= factor
    end
end

function updateCR(CHf::CHform, gate::String, qb::String)
    if (gate == "S")
        q = parse(Int8, qb)
        CHf.M[:,q] .= mod.(CHf.M[:,q] + CHf.F[:,q],2)
        CHf.y[:] .= mod.(CHf.y[:] - CHf.F[:,q],4)
    elseif (gate == "Z")
        q = parse(Int8, qb)
        CHf.y[:] .= mod.(CHf.y[:] - 2*CHf.F[:,q],4)
    elseif (gate == "CZ") 
        qbidx = findfirst(" ", qb)[1]
        q = parse(Int8,qb[1:qbidx])
        r = parse(Int8,qb[qbidx+1:end])
        CHf.M[:,q] .= mod.(CHf.M[:,q] + CHf.F[:,r],2) 
        CHf.M[:,r] .= mod.(CHf.M[:,r] + CHf.F[:,q],2) 
        CHf.y[:] .= mod.(CHf.y[:] + 2*CHf.F[:,q].*CHf.F[:,r],4)
    elseif (gate == "CX")
        qbidx = findfirst(" ", qb)[1]
        q = parse(Int8,qb[1:qbidx])
        r = parse(Int8,qb[qbidx+1:end])
        CHf.G[:,q] .= mod.(CHf.G[:,q] + CHf.G[:,r],2)
        CHf.F[:,r] .= mod.(CHf.F[:,r] + CHf.F[:,q],2)
        CHf.M[:,q] .= mod.(CHf.M[:,q] + CHf.M[:,r],2)
    end
end

function updateH(CHf::CHform, qb::String)
    q = parse(Int8, qb)
    #construct t,u,v_, a, b
    t = zeros(Int8, CHf.N)
    u = zeros(Int8, CHf.N)
    v_ = zeros(Int8, CHf.N)
    alpha = 0
    beta = 0
    # t,u,v_, alpha,beta = zeros((CHf.N,1), dtype = int), zeros((CHf.N,1), dtype = int), zeros((CHf.N,1), dtype = int), 0, 0
    v_ = Vector{Int8}(mod.(1 .-CHf.v,2))
    t[:] = mod.(CHf.s[:] + CHf.G[q,:] .*CHf.v[:],2)
    u[:] = mod.(CHf.s[:] + CHf.F[q,:] .*v_[:] + CHf.M[q,:] .*CHf.v[:], 2)
    # alpha = Int64(CHf.G[q,:]*( v_[:].*CHf.s[:]))
    alpha = Int64(dot(CHf.G[q,:], v_[:].*CHf.s[:]))
    # beta = Int64(CHf.M[q,:]*(v_[:].*CHf.s[:]) + CHf.F[q,:]*(CHf.v[:].*(CHf.M[q,:] + CHf.s[:])))
    beta = Int64(dot(CHf.M[q,:],v_[:].*CHf.s[:]) + dot(CHf.F[q,:],(CHf.v[:] .*(CHf.M[q,:] + CHf.s[:]))))
    #constuct V0 and V1 
    V0::Vector{Int8} = []
    V1::Vector{Int8} = []
    for k::Int8 in 1:CHf.N
        if t[k] != u[k] #building V0, V1
            if CHf.v[k] == 0 
                append!(V0, [k]) 
            else 
                append!(V1, [k])
            end
        end
    end
    #case t == u
    if isempty(V0) & isempty(V1) #empty lists have logical value False
        CHf.s = deepcopy(t)
        CHf.w *= ((-1)^alpha + (1im^CHf.y[q])*(-1)^beta)/sqrt(2)
    else #case t != u
        #We avoid having to construct Vc as it is constructed in the same direction as the gates are applied
        #i.e. if a CX gate is added in Vc before a CZ gate, the CX gate will be applied to Uc before the CZ gate
        #This way we can get rid of a small loop and a lists-in-a-list at the cost of some clarity in the code
        if ~isempty(V0) #computing the index p of the first qubit where u[p] != t[p]
            p = V0[1]
        else
            p = V1[1]
        end
        #computing some useful values first
        if t[p] == 1
            x = deepcopy(u)
            x[p] = mod.((x[p] + 1),2)
            z = deepcopy(u)
        else 
            x = deepcopy(t)
            z = deepcopy(t)
            z[p] = mod.((z[p] + 1),2)
        end
        delta = Int8(mod(CHf.y[q] + 2*(alpha+beta),4))
        a,b,c = valanalysis(CHf, Int8(x[p]), Int8(CHf.v[p]), delta)
        #"constructing Vc"
        if ~isempty(V0) 
            for k in V0
                if k != p
                    @views updateCR(CHf,"CX",string(string(p) ," " ,string(k)))
                end
            end
            for k in V1
                @views updateCR(CHf,"CZ",string(string(p), " " ,string(k)))
            end
        else #V0 == []
            for k in V1
                if k != p
                    @views updateCR(CHf,"CX",string(string(k), " " ,string(p)))
                end
            end
        end
        #"constructing Wc = VcS^a"
        if a == 1 
            @views updateCR(CHf, "S",string(p))
        end
        #Update Uh
        CHf.v[p] = b
        #update omega
        CHf.w *= (-1)^alpha
        #update s
        CHf.s = x
        CHf.s[p] = c
    end
end

function valanalysis(CHf::CHform, xq::Int8, vq::Int8, delta::Int8) #return the values of a,b,c as required by the algorithm (equation between 53 and 54 of the paper)
    if vq == 1
        if delta == 0
            a = 0
            b = 0
            c = 0
            # a, b, c = 0, 0, 0
        elseif delta == 2
            CHf.w *= (-1)^xq
            a = 0
            b = 0
            c = 1
            # a, b, c = 0, 0, 1
        elseif delta == 1
            a = 1
            b = 1
            c = mod(xq+1,2)
            # a, b, c = 1, 1, (xq+1)%2
            CHf.w *= (1+1im)/sqrt(2)
        elseif delta == 3
            a = 1
            b = 1
            c = xq
            # a, b, c = 1, 1, xq
            CHf.w *= (1-1im)/sqrt(2)
        end
    else 
        if delta == 0
            a = 0
            b = 1
            c = 0
            # a, b, c = 0, 1, 0
        elseif delta == 2
            a = 0
            b = 1
            c = 1
            # a, b, c = 0, 1, 1
        elseif delta == 1
            a = 1
            b = 1
            c = xq
            # a, b, c = 1, 1, xq
        elseif delta == 3
            a = 1
            b = 1
            c = mod(xq+1,2)
            # a, b, c = 1, 1, (xq+1)%2
        end
        #changing the phase
        if xq == 1
            CHf.w *= (1im)^delta 
        end
    end
    return a,b,c
end
# ---------------------------SCALAR PRODUCT FUNCTIONS ------------------------------------------------------------------
function spbase(CHf::CHform, x::Vector{Int8}) #Computes the scalar product <x|CHf> where x is an element of the computational basis and CHf a stabilizer state in CH-form
    #initializing useful variables
    mu = Int8(0)
    vnorm = 0.0 #float64 by default
    u = zeros(Int8, CHf.N)
    t = zeros(Int8, CHf.N)
    #computing these variables (O(N²))
    # u = dot(x, CHf.F)
    u = transpose(x)*CHf.F
    u = mod.(u,2)
    # t = dot(x, CHf.M)
    t = transpose(x)*CHf.M
    t = mod.(t,2)
    #compute mu 
    # yt = transpose(CHf.y)
    b = triangularprdt(x, CHf.F, CHf.M, CHf.N)
    # mu = dot(yt, x) + 2*b
    mu = dot(CHf.y, x) + 2 .*b
    mu = mod(mu, 4)
    #compute |v| 
    vnorm = norm(CHf.v, 1) #norm 1 of v
    #computing the factors inside the probability amplitude
    famp = 1
    fm1 = 1
    for p in 1:CHf.N
        if CHf.v[p] == 1
            fm1 *= (-1)^(u[p]*CHf.s[p])
        elseif u[p] != CHf.s[p]
            famp *= 0 
        end
    end
    return Complex(2^(-vnorm/2) * 1im^mu * fm1 * famp * CHf.w) 
end

function conjugateCHf(CHf::CHform) #conjugate (transpose) of the state represented by Chf
    tmp::CHform = deepcopy(CHf)
    tmp.w = conj(tmp.w)
    #inverting matrices F and G
    tmp.F = invBinMat(tmp.F)
    tmp.G = invBinMat(tmp.G)
    #computing the new M matrix
    tmp.M = mod.((-1) .*(tmp.F*tmp.M*tmp.G),2) #Important: keep the updates of F, G and M in this order because the update of M uses the updated versions of F and G
    #computing the new gamma
    for q::Int8 in 1:tmp.N
        tmp.y[q] = mod(-dot(tmp.F[q,:],CHf.y[:]) - 2 .*croppedTriangularPrdt(q,tmp, CHf),4)
    end
    return tmp
end

function spCHf(ch2::CHform, ch1::CHform) #Corresponds to <phi_2|phi_1>, ch2 must be in bra form
    #conjugating ch2: 
    # ch2conj = ch2.conjugateCH()
    #multiplying the two Uc matrices together
    newch = prdtUc(ch2, ch1)
    #Adding the H gates to newch
    for k::Int8 in 1:ch1.N
        if ch2.v[k] == 1
            updateH(newch, string(k))
        end
    end
    #Now we have a scalar product of type <s|phi'> where |phi'> is represented by newch
    #computing <s|phi'>
    amp = spbase(newch, ch2.s[:])
    # #now multiplying by (omega_2)* since it is not computed in the amplitude method
    amp *= ch2.w
    return amp
end

function scalPrdt(ch2::CHform, ch1::CHform, J::Float64 = NaN, h::Float64 = NaN) #computes either <ch2|ch1> or <ch2|H|ch1> for the transverse field Ising model Hamiltonian, depending on values of J and h
    if isnan(J) || isnan(h) #the values of J and h aren't defined, we simply compute the scalar product <ch2|ch1>
        return spCHf(ch2, ch1)
    else #if J and h are given, then we compute the expectation value of the Transverse Ising model Hamiltonian <ch2|H|ch1>
        amp = 0.0
        for k in 1:ch1.N
            tmp = deepcopy(ch1)
            #apply the Transverse Ising Hamiltonian to ch1, qubit by qubit
            #first apply the Ising part for qubit k
            updateCL(tmp, "Z", string(k))
            updateCL(tmp, "Z", string(mod(k,tmp.N)+1))
            #compute the expectation value
            amp += -J*spCHf(ch2,tmp)
            #undo the change (more efficient than taking a clean copy of ch1 as Z gates take O(1) to be applied)
            updateCL(tmp, "Z", string(k))
            updateCL(tmp, "Z", string(mod(k,tmp.N)+1))
            # tmp = deepcopy(ch1)
            #Now apply the transverse field part, represented by X gates
            updateCL(tmp, "X", string(k))
            amp += -h*spCHf(ch2, tmp)
        end
        return amp
    end
end

# ---------------------------UTILITY FUNCTIONS -------------------------------------------------------------------------
function triangularprdt(x::Vector{Int8},F::Matrix{Int8},M::Matrix{Int8},N::Int8) #implementation of what I call a "triangular product"
    #a "Triangular Product" is defined as (copy/paste in a LaTeX file if it is too unreadable)
    #b = \sum^n_{j=1}\sum^n_{l=1}x_l F_{l,j}(\sum^n_{k=l} x_k M_{k,j})nonzero F,M are the matrices in the CH-form and x is the a n-bit string corresponding to the comp. basis state in the computation of <x|phi> 
    #(see the Amplitude method of the CHform class for more information)
    xFxM = 0
    for l in 1:N 
        #computation of the sum of x_{k}M_{k,j} (elements k = 2,...,n ) 
        # xM = transpose(x[l:end]) *M[l:end,:] #with "l:" element l is included, but with ":l" it is not
        xFxM += dot(x[l] .*F[l,:], transpose(x[l:end]) *M[l:end,:])
    end
    return xFxM
end

function invBinMat(Ain::Matrix{Int8}) #inverts a binary matrix with det=1 using Gaussian elimination
    A = deepcopy(Ain)
    N = Int8(size(A)[1])
    usedIdx = [] #list of indices of all rows that already have been used as pivots
    Ainv = Matrix{Int8}(I,N,N)
    for k in 1:N[1] #sumation of lines to achieve a scrambled identity matrix
        # print('line k = ',k)
        idx = [i for i in 1:N if A[i,k] != 0] #list of indices of 1's in the k^th column of A
        if length(idx) > 1 #checking if the column has more than one non-zero element
            l = Int64(1)
            while idx[l] in usedIdx
                l = l+1
            end
            pvt = idx[l] #Pivot point, must be non-zero for all matrices in our case
            for l in idx
                if l != pvt 
                    A[l,:] .= mod.(A[l,:]+ A[pvt,:],2)
                    Ainv[l,:] .= mod.(Ainv[l,:]+ Ainv[pvt,:],2)
                end
            end
            append!(usedIdx, pvt) #adding the pivot row index
        else
            append!(usedIdx, idx[1]) #adding the pivot row index
        end
    end
    #Now the matrix should ressemble the identity matrix, with scrambled rows
    #swapping lines to recover the identity
    for k in 1:N[1]-1
        if A[k,k] != 1
            # print('k = ',k)
            # print(A[k+1:,k])
            OneIdx = [i for i in k+1:N if A[i,k] != 0][1]#OneIdx should be a single element array
            # OneIdx = k+1+OneIdx[1] #the index in OneIdx is the number of rows below k+1 (included), hence to recover the true row index we have to add k+1
            #swap row k with row OneIdx[1] of A (to keep track of which row have already been swapped)
            tmp = A[k,:]
            A[k,:] = A[OneIdx,:]
            A[OneIdx,:] = tmp 
            #do the same swap for Ainv
            tmp = Ainv[k,:]
            Ainv[k,:] = Ainv[OneIdx,:]
            Ainv[OneIdx,:] = tmp 
        end
    end
    return Ainv 
end

function croppedTriangularPrdt(q::Int8,ch1::CHform, ch2::CHform) #computes the "Triangular product" without diagonal used in the computation of Uc1*Uc2
    ab = Int64(0)
    for j::Int8 in 2:ch1.N
        #computation of the sum of x_{k}M_{k,j} (elements k = 2,...,n )
        ab += ch1.F[q,j] * ((transpose(ch1.F[q,1:j-1])*ch2.M[1:j-1,:] + transpose(ch1.M[q,1:j-1])*ch2.G[1:j-1,:])*ch2.F[j,:])
    end 
    return ab
end

function prdtUc(ch2::CHform, ch1::CHform)#computes F,M,G,y when ch1's matrix Uc has been multiplied to the left by a matrix Uc of another CH-form ch2, i.e. computes Uc2*Uc1
    newch::CHform = deepcopy(ch1) #Computing a scalar product does not change the vector, hence a copy is made
    #computing gamma
    Threads.@threads for q::Int8 in 1:newch.N
        newch.y[q] = mod(ch2.y[q] + transpose(ch2.F[q,:])*ch1.y[:] + 2*croppedTriangularPrdt(q,ch2, ch1),4)
    end
    #computing the matrices 
    newch.F .= mod.(ch2.F*ch1.F,2)
    newch.M .= mod.(ch2.F*ch1.M .+ ch2.M*ch1.G,2)
    newch.G .= mod.(ch2.G*ch1.G,2)
    return newch
end
