using Base: String, Bool, Complex
# using Core: Matrix, Vector
using LinearAlgebra
using Random 
using ProgressMeter

include("CHform.jl")
include("GSmat.jl")
include("MetroHast.jl")

mutable struct ParTemp
    MHlist::Vector{VarMH} #list of all the parallel MC 
    BetaList::Matrix{Float64} #matrix containing the lists of temperatures at each step every MC
    Nmc::Int8  #number of Markov Chains
    Nsteps::UInt64 #number of steps in the MCs
    #Constructor
    # ParTemp(N,R,L, J, h, Nmc, TinList, TfinList) = new([VarMH(N,R,L, J, h) for _ in 1:Nmc], Matrix{Float64}([[1 ./(TinList[k] .+ (TfinList[k] - TinList[k])*(log.(LinRange(2., 10., Nsteps)) .-log(2))/(log(10)-log(2)))] for k in length(TinList)]), Nmc)
    # ParTemp(N::Int8,R::Int64,L::Int8, J::Float64, h::Float64, Nmc::Int8, Nsteps::UInt64, TinList::Vector{Float64}, TfinList::Vector{Float64 }) = new([VarMH(N,R,L, J, h) for _ in 1:Nmc], initBetalist(TinList, TfinList, Nmc, Nsteps), Nmc, Nsteps)
    ParTemp(N::Int8,R::Int64,L::Int8, J::Float64, h::Float64, Nmc::Int8, Nsteps::UInt64, TinList::Vector{Float64}, TfinList::Vector{Float64 }) = new(initMHlisttoHF([VarMH(N,R,L, J, h) for _ in 1:Nmc]), initBetalist(TinList, TfinList, Nmc, Nsteps), Nmc, Nsteps)
end

function run(Ptp::ParTemp)
    # Computing the analytical energy to check that the result is not abhorrent
    m = 1:Ptp.MHlist[1].HPm.CHlist[1].N
    m = m .- floor(Ptp.MHlist[1].HPm.CHlist[1].N/2) #element-wise substraction
    #computation of k
    k = 2*pi .*m ./Ptp.MHlist[1].HPm.CHlist[1].N
    #Method 2
    E0 = -1 * sum(sqrt.((-Ptp.MHlist[1].HPm.h/Ptp.MHlist[1].HPm.J .- cos.(k)).^2 + sin.(k).^2))
    @showprogress for i in 1:Ptp.Nsteps
        #do a step 
        step(Ptp, Ptp.BetaList[:,i])
        #if the result is attained sooner than at the end of the Nsteps, break
        if minimum([k.E for k in Ptp.MHlist])-E0 < 1e-14
            break
        end
    end
    return [k.E for k in Ptp.MHlist]
end

function run_log(Ptp::ParTemp)
    # Computing the analytical energy to check that the result is not abhorrent
    m = 1:Ptp.MHlist[1].HPm.CHlist[1].N
    m = m .- floor(Ptp.MHlist[1].HPm.CHlist[1].N/2) #element-wise substraction
    #computation of k
    k = 2*pi .*m ./Ptp.MHlist[1].HPm.CHlist[1].N
    #Method 2
    E0 = -1 * sum(sqrt.((-Ptp.MHlist[1].HPm.h/Ptp.MHlist[1].HPm.J .- cos.(k)).^2 + sin.(k).^2))

    Elog = zeros(Float64, Ptp.Nsteps, Ptp.Nmc)
    Tlog = zeros(Float64, Ptp.Nsteps, Ptp.Nmc)
    @showprogress for i in 1:Ptp.Nsteps
        #do a step 
        step(Ptp, Ptp.BetaList[:,i])
        Elog[i,:] = [k.E for k in Ptp.MHlist]
        Tlog[i,:] = 1. ./Ptp.BetaList[:,i]
        #if the result is attained sooner than at the end of the Nsteps, break
        if minimum([k.E for k in Ptp.MHlist])-E0 < 1e-14
            break
        end
    end
    return Elog, Tlog
end

function step(Ptp::ParTemp, CurrBetalist::Vector{Float64})
    #choose to swap temperatures or not 
    b = Int8(rand(1:30))
    mc1 = Int8(0)
    mc2 = Int8(0)
    swapIdx = zeros(Int8, 2)
    if b == 1 #swap
        
        #select which swap to do
        mc1 = Int8(rand(1:Ptp.Nmc))
        mc2 = Int8(rand(1:Ptp.Nmc))
        while mc1==mc2 & Ptp.Nmc>1
            mc2 = Int8(rand(1:Ptp.Nmc))
        end
        #approve swap 
        p = rand(Float64) #generates a random number in [0,1)
        # P = min(1, exp((Ptp.MHlist[mc1]-Ptp.MHlist[mc2])*(CurrBetalist[mc1] - CurrBetalist[mc2])))
        if p < min(1, exp((Ptp.MHlist[mc1].E-Ptp.MHlist[mc2].E)*(CurrBetalist[mc1] - CurrBetalist[mc2])))
            #storing the swap indices
            swapIdx[1] = mc1
            swapIdx[2] = mc2
            #todo faire le swap
            tmpT = copy(Ptp.BetaList[mc1,:]) 
            Ptp.BetaList[mc1,:] = copy(Ptp.BetaList[mc2,:]) 
            Ptp.BetaList[mc2,:] = tmpT
        else #update the MC normally
            # step(Ptp.MHlist[mc1], CurrBetalist[mc1])
            # step(Ptp.MHlist[mc2], CurrBetalist[mc2])
            Threads.@threads for i in [mc1 mc2]
                step(Ptp.MHlist[i], CurrBetalist[i])
            end
        end
    end
    #In any other cases, the Markov chains are updated as usual
    Threads.@threads for i  in filter(x -> x ∉ swapIdx, 1:Ptp.Nmc) #what this does is creating a range 1:Ptp.Nmc but filters out the elements in swapIdx
        @views step(Ptp.MHlist[i], CurrBetalist[i])
    end
end

#-----------------------Utility function for the constructor----------------------------------------------
function initBetalist(TinList::Vector{Float64}, TfinList::Vector{Float64}, Nmc::Int8, Nsteps::UInt64) #TODO
    out = zeros(Float64, Nmc, Nsteps)
    for k in 1:Nmc
        out[k,:] = 1 ./(TinList[k] .+ (TfinList[k] - TinList[k])*(log.(LinRange(2., 10., Nsteps)) .-log(2))/(log(10)-log(2)))
    end
    return out
end

function initMHlisttoHF(MHlist::Vector{VarMH}) #initialize the list of Markov Chains to the Hartree-Fock state
    # return initVMHtoHF.(deepcopy(MHlist))
    out::Vector{VarMH} = deepcopy(MHlist)
    for k in out
        initVMHtoHF(k)
    end
    return out
end

#----------Run with comp. basis states ---------------------------------
function runCompBasis(Ptp::ParTemp)
    # Computing the analytical energy to check that the result is not abhorrent
    m = 1:Ptp.MHlist[1].HPm.CHlist[1].N
    m = m .- floor(Ptp.MHlist[1].HPm.CHlist[1].N/2) #element-wise substraction
    #computation of k
    k = 2*pi .*m ./Ptp.MHlist[1].HPm.CHlist[1].N
    #Method 2
    E0 = -1 * sum(sqrt.((-Ptp.MHlist[1].HPm.h/Ptp.MHlist[1].HPm.J .- cos.(k)).^2 + sin.(k).^2))
    @showprogress for i in 1:Ptp.Nsteps
        #do a step 
        stepCompBasis(Ptp, Ptp.BetaList[:,i])
        #if the result is attained sooner than at the end of the Nsteps, break
        if minimum([k.E for k in Ptp.MHlist])-E0 < 1e-14
            break
        end
    end
    return [k.E for k in Ptp.MHlist]
end

function stepCompBasis(Ptp::ParTemp, CurrBetalist::Vector{Float64})
    #choose to swap temperatures or not 
    b = Int8(rand(1:30))
    mc1 = Int8(0)
    mc2 = Int8(0)
    swapIdx = zeros(Int8, 2)
    if b == 1 #swap
        #select which swap to do
        mc1 = Int8(rand(1:Ptp.Nmc))
        mc2 = Int8(rand(1:Ptp.Nmc))
        while mc1==mc2 & Ptp.Nmc>1
            mc2 = Int8(rand(1:Ptp.Nmc))
        end
        #approve swap 
        p = rand(Float64) #generates a random number in [0,1)
        # P = min(1, exp((Ptp.MHlist[mc1]-Ptp.MHlist[mc2])*(CurrBetalist[mc1] - CurrBetalist[mc2])))
        if p < min(1, exp((Ptp.MHlist[mc1].E-Ptp.MHlist[mc2].E)*(CurrBetalist[mc1] - CurrBetalist[mc2])))
            #storing the swap indices
            swapIdx[1] = mc1
            swapIdx[2] = mc2
            #todo faire le swap
            tmpT = copy(Ptp.BetaList[mc1,:]) 
            Ptp.BetaList[mc1,:] = copy(Ptp.BetaList[mc2,:]) 
            Ptp.BetaList[mc2,:] = tmpT
        else #update the MC normally
            # step(Ptp.MHlist[mc1], CurrBetalist[mc1])
            # step(Ptp.MHlist[mc2], CurrBetalist[mc2])
            Threads.@threads for i in [mc1 mc2]
                stepCompBasis(Ptp.MHlist[i], CurrBetalist[i])
            end
        end
    end
    #In any other cases, the Markov chains are updated as usual
    Threads.@threads for i  in filter(x -> x ∉ swapIdx, 1:Ptp.Nmc) #what this does is creating a range 1:Ptp.Nmc but filters out the elements in swapIdx
        @views stepCompBasis(Ptp.MHlist[i], CurrBetalist[i])
    end
end