using Base: String, Bool, Complex
# using Core: Matrix, Vector
using LinearAlgebra
using Random 
using ProgressMeter

include("CHform.jl")
include("GSmat.jl")
include("MetroHast.jl")

mutable struct QNG
    HPm::HPmat #matrices H and P, coefficient
    E::Float64 #energy
    Anglelist::Array{Int8} #Array (3 dimensional) storing the all angles (for all stabilizers and layers, lines = layers ; columns = stabilizers states; 3rd dimension = Angles)
    L::Int8 #number of layers
    g::Vector{Vector{Matrix{Float64}}} #Quantum Geometric tensor: First vector - stabilizers, second vector - layers, Matrix - parameters
    nL::Vector{Vector{Vector{Complex{Float64}}}} #gradient of L: same arrangement as the QGT
    Ansatz::String
    #Constructor
    QNG(N::Int8, R::Int64, L::Int8,  J::Float64, h::Float64, Ansatz::String) = new(HPmat(N, R, J, h), 0. , zeros(Int8, L, R, 2*N), L::Int8, [[zeros(Float64, 2*N, 2*N) for _i in 1:L] for _ in 1:R],[[zeros(Float64, 2*N) for _i in 1:L] for _ in 1:R], Ansatz)
end

function run(Qng::QNG, Etalist::Vector{Float64})
    #----run the gradient descent----
    for k in 1:length(Etalist)
        step(Qng, Etalist[k])
        println("angles = ",Qng.Anglelist)
        println("g = ",Qng.g)
        println("nL = ",Qng.nL)
    end
    # #Trick to avoid having a singular QGT
    # Qng.g[r][l] += 1e-4*I #I is the identity matrix
    #----compute the energy----
    #First compute H and P 
    compHP(Qng.HPm) #the list of stabilizer states is already initialized to {|phi^L_r>}_r=1,...,R
    #solve the generalized eigenvalue problem
    E, alphas = eigen(Qng.HPm.H, Qng.HPm.P)
    Qng.E, idxMin = findmin(E)
    Qng.HPm.Phaselist = alphas[:, idxMin]
end

function step(Qng::QNG, eta::Float64) #eta is the step size
    #compute the QGT and the gradient of L 
    computeQGTandGradL(Qng)
    #solving the system of equations g*dTheta = -eta*nL, w/ dTheta = theta_{t+1} - theta_t
    for r::Int64 in 1:Qng.HPm.R 
        for l::Int8 in 1:L
            #computing the Moore-Penrose pseudoinverse of the block QGT
            gpinv = pinv(Qng.g[r][l], 1e-2) #the second argument is the relative tolerance, i.e. singular values smaller than this (times the largest sing. val) will not be inverted
            # println(Qng.g[r][l])
            #solving the system of equations
            # dtheta = -gpinv*Qng.nL[r][L]
            dtheta = -Qng.nL[r][L]
            println("dtheta = ",dtheta)
            # for i in 1:length(dtheta)

            #compting theta[t+1]
            Qng.Anglelist[l,r,:] += dtheta
            #note the energy doesn't have to be computed at each step, saving a lot of computation time
        end
    end
end

function computeQGTandGradL(Qng::QNG)
    for r::Int64 in 1:Qng.HPm.R #stabilizers 
        stab_r = CHform(Qng.HPm.CHlist[1].N)
        for l::Int8 in 1:L #layers
            U_l = genVarCirc(Qng, l, r, Qng.Ansatz)
            update_circ(stab_r, U_l)
            #todo
            computeQGTandGradblock(Qng, stab_r, r, l)
        end
        #change the CHlist after computing phi_r[1:L]
        Qng.HPm.CHlist[r] = stab_r
    end
end
function computeQGTandGradblock(Qng::QNG, stab_rl::CHform, r::Int64, l::Int8) #computes a block of the QGT and the gradient of L corresponding to stabilizer r and layer l
    # for r::Int64 in 1:Qng.HPm.R #stabilizers 
    #     stab_r = CHform(Qng.HPm.CHlist[1].N)
    #     for l::Int8 in 1:L #layers
    #         U_l = genVarCirc(Qng, l, r)
    #         update_circ(stab_r, U_l)
    for i in 1:length(Qng.Anglelist[:,:,1]) #à voir pour une optimisation ici
        #contructing the Ki generator of rotation
        # qbi::Int64 = Int64(mod(i,2))
        # Ki::Vector{Vector{String}} = [[],[]]
        # if qbi == 1
        #     Ki = [["Id", "X"],[string(0.5), string(Int64((i+1)/2))]]
        # else
        #     Ki = [["Id", "Z"],[string(0.5), string(Int64(i/2))]]
        # end
        Ki = rotationGen(i, qng.Ansatz)
        #contructing the state Ki|phi_r>
        tmpi = deepcopy(stab_rl)
        update_circ(tmpi, Ki)
        #----compute the gradient of L------
        #computing -i/2 * <phi_r|HKi|phi_r>
        tmpibra =conjugateCHf(tmpi)
        qng.nL[r][l][i] = -1im/2*scalPrdt(tmpibra, tmpi, Qng.HPm.J, Qng.HPm.h)
        print("-------- tut ",qng.nL[r][l][i], " ")
        #add the complex conjugate i/2 * <phi_r|KiH|phi_r>
        qng.nL[r][l][i] += conj(qng.nL[r][l][i])
        # qng.nL[r][l][i] *= -1im/2 
        #----compute the QGT ---------------
        for j in 1:length(Qng.Anglelist[:,:,1])
            #contructing the Ki generator of rotation
            qbj::Int64 = mod(j,2)
            # Kj::Vector{Vector{String}} = [[],[]]
            # if qbj == 1
            #     Kj = [["Id", "X"],[string(0.5), string(Int64((j+1)/2))]]
            # else
            #     Kj = [["Id", "Z"],[string(0.5), string(Int64(j/2))]]
            # end
            Kj = rotationGen(j, qng.Ansatz)
            #contructing the state <phi_r|Kj
            tmpj = deepcopy(stab_rl)
            update_circ(tmpj, Kj)
            tmpj = conjugateCHf(tmpj)
            #computing g_ij = Re(<phi_r|KjKi|phi_r> - <phi_r|Ki|phi_r><phi_r|Kj|phi_r>) i.e. the QGT
            Qng.g[r][l][i,j] = real(scalPrdt(tmpj, tmpi) - scalPrdt(stab_rl, tmpi)*scalPrdt(tmpj, stab_rl))
        end
    end
    #     end
    # end
end

function rotationGen(Param::Int64, AnsatzType::String) #computes a "block" of GradL corresponding to stabilizer r and layer l
    K::Vector{Vector{String}} = [[],[]]
    qb::Int64 = mod(Param,2)
    if AnsatzType == "Y"
        if qb == 1
            K = [["Id", "X", "Z"],[string(0.5*1im), string(Int64((Param+1)/2)),string(Int64((Param+1)/2))]]
        else
            K = [["Id", "X" ,"Z"],[string(0.5*1im), string(Int64(Param/2)),string(Int64(Param/2))]]
        end
    else 
        if qb == 1
            K = [["Id", "X"],[string(0.5), string(Int64((Param+1)/2))]]
        else
            K = [["Id", "Z"],[string(0.5), string(Int64(Param/2))]]
        end
    end
    return K
end

#-------------Initialisation----------------------
function initQNGtoHF(Qng::QNG, recursion::Int8 = Int8(0))
    #updating the CH-form list with the HF state
    # println("stab = init ",Qng.HPm.P)
    VarCirc = genVarCirc(Qng, Qng.L, 1, Qng.Ansatz)
    updateHP(Qng.HPm, 1, VarCirc,1)
    # println("stab = 1 ",Qng.HPm.P)
    # println("=================")
    #generating other HF + small excitation states (the small excitations are generated by changing the angles on Rz gates)
    for r::Int64 in 2:Qng.HPm.R
        tmpHP = genRandHFstab(Qng, r)
        maxTry = 0
        expOffset::Int8 = Int8(0)
        while (rank(tmpHP.P)<r) & (maxTry <= 2*Qng.HPm.CHlist[1].N)
            # println("tut")
            maxTry += 1
            #retry 
            tmpHP = genRandHFstab(Qng, r, expOffset)
            if (maxTry >= 2*Qng.HPm.CHlist[1].N) #there is probably no indep stabilizer obtainable by changing only expOffset angles, try expOffset + 1 angles
                # println("reset")
                maxTry = 0
                expOffset += Int8(1)
            end
            # println(Qng.Anglelist)
        end
        # println("rank = ",rank(tmpHP.P))
        #accept the change
        Qng.HPm = deepcopy(tmpHP)
    end
    # println(rank(Qng.HPm.P))
    # println("=================")
    # println(Qng.HPm.P)
    #computation of the energy and the coefficients of the decomposition
    vals, vecs = eigen(Qng.HPm.H, Qng.HPm.P) 
    # println(vals)
    #find the minimal energy
    Qng.E = 0.0
    # print(typeof(Qng.E))
    Qng.HPm.Phaselist = vecs[:,1] 
    for k in 1:length(vals)
        if (real(vals[k]) < Qng.E) & (abs(imag(vals[k])) < 1e-14)
            Qng.E = real(vals[k])
            Qng.HPm.Phaselist = vecs[:, k] 
        end
    end
    #Error Handling
    #Stop if we can't find a non-NaN value for E
    if recursion >= 10
        throw(DomainError(recursion,"Cannot find an initial state without NaN energy"))
    end
    #If the value is NaN, try again until recursion is too high
    if isnan(Qng.E)
        initVMHtoHF(Qng, recursion + 1)
    end
    #normalize the coefficients in the decomposition
    # nrm = sqrt(dot(Qng.HPm.Phaselist, conj(Qng.HPm.Phaselist)))
    nrm = norm(Qng.HPm.Phaselist) #default norm for the norm function: Euclidean
    Qng.HPm.Phaselist = Qng.HPm.Phaselist/nrm
end

function genVarCirc(Qng::QNG, layer::Int8, stabIdx::Int64, Ansatz::String = "RZ")
    if Ansatz == "Y"
        return genVarCircY(Qng, layer, stabIdx)
    else 
        return genVarCircRZ(Qng, layer, stabIdx)
    end
end

function genVarCircRZ(Qng::QNG, layer::Int8, stabIdx::Int64) #Generates a variational quantum circuit similar to one used in VQE, using angles found in Qng
    circ::Vector{Vector{String}} = [[],[]]
    for k::Int8 in 1:Qng.HPm.CHlist[1].N
        #rotation gates %TODO
        # rxcirc = RxGate(Qng.Anglelist[2*Qng.HPm.CHlist[1].N*(Qng.R*layer + stabIdx)+k], k)
        rxcirc = RxGate(Qng.Anglelist[layer,stabIdx, 2*k-1], k)
        rzcirc = RzGate(Qng.Anglelist[layer,stabIdx, 2*k], k)
        append!(circ[1], append!(rxcirc[1], rzcirc[1]))
        append!(circ[2], append!(rxcirc[2], rzcirc[2]))
    end
    for k::Int8 in 1:Qng.HPm.CHlist[1].N
        #CNOT gates
        append!(circ[1], ["CX"])
        append!(circ[2], [string(string(k), " ", string(mod(k,Qng.HPm.CHlist[1].N)+1))])
    end
    # print(circ)
    return circ
end

function genVarCircY(Qng::QNG, layer::Int8, stabIdx::Int64) #Generates a variational quantum circuit similar to one used in VQE, using angles found in Qng
    circ::Vector{Vector{String}} = [[],[]]
    for k::Int8 in 1:Qng.HPm.CHlist[1].N
        #rotation gates %TODO
        # rxcirc = RxGate(Qng.Anglelist[2*Qng.HPm.CHlist[1].N*(Qng.R*layer + stabIdx)+k], k)
        rycirc = RyGate(Qng.Anglelist[layer,stabIdx, 2*k-1], k)
        # rzcirc = RzGate(Qng.Anglelist[layer,stabIdx, 2*k-1], k)
        append!(circ[1], append!(rycirc[1]))
        append!(circ[2], append!(rycirc[2]))
    end
    for k::Int8 in 1:Qng.HPm.CHlist[1].N
        #CNOT gates
        append!(circ[1], ["CX"])
        append!(circ[2], [string(string(k), " ", string(mod(k,Qng.HPm.CHlist[1].N)+1))])
    end
    for k::Int8 in 1:Qng.HPm.CHlist[1].N
        #rotation gates %TODO
        # rycirc = RyGate(Qng.Anglelist[2*Qng.HPm.CHlist[1].N*(Qng.R*layer + stabIdx)+k], k)
        rycirc = RyGate(Qng.Anglelist[layer,stabIdx, 2*k], k)
        # rzcirc = RzGate(Qng.Anglelist[layer,stabIdx, 2*k], k)
        append!(circ[1], append!(rycirc[1]))
        append!(circ[2], append!(rycirc[2]))
    end
    # print(circ)
    return circ
end

function genRandHFstab(Qng::QNG, r::Int64, expOffset::Int8 = Int8(0), AnsatzType = "XZ")
    # r is the index of the current stabilizer
    tmpHP = deepcopy(Qng.HPm)
    n = Int64(Qng.HPm.CHlist[1].N)
    #computing the amount angles to be changed
    count::Int64 = n
    expnt::Int8 = Int8(1 + expOffset)
    while (r > count)
        expnt += Int8(1)
        count += n ^expnt
        if expnt >= 2*n
            throw(DomainError(expnt,"Too many stabilizers"))
        end
    end
    #changing the required amount of angles to generate a new state
    usedExp::Vector{Int8} = []
    for c in 1:expnt
        #trying some qubit
        idx = Int8(rand(1:2*Qng.HPm.CHlist[1].N))
        while idx in usedExp #checking the qubit has not already been affected
            idx = Int8(rand(1:2*Qng.HPm.CHlist[1].N))
        end
        append!(usedExp, [idx])
        #choosing randomly a new angle
        nwangl = Int8(rand(1:8))
        Qng.Anglelist[Qng.L, r, idx] = nwangl
    end
    #constructing the variational circuit
    # if AnsatzType == "Y"
    #     varCirc = genVarCircY(Qng, Qng.L, r)
    # else
    #     varCirc = genVarCirc(Qng, Qng.L, r)
    # end
    varCirc = genVarCirc(Qng, Qng.L, r, Qng.Ansatz)
    #applying the change
    updateHP(tmpHP, r, varCirc, r)
    return tmpHP
end

function RyGate(exponent::Int8, qb::Int8)
    circ = [[],[]]
    if mod(exponent,8) > 3
        append!(circ[1], ["Id"])
        append!(circ[2], [string(-1)])
    end
    for k in 1:Int8(mod(exponent,4))
        append!(circ[1], ["H","Z"])
        append!(circ[2], [string(qb), string(qb)])
    end
    return circ
end