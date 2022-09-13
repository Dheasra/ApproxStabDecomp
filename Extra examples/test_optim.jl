using Base: String, Bool, Complex
# using Core: Matrix, Vector
using LinearAlgebra
using Statistics
using PyCall
using Dates
using DelimitedFiles
using BenchmarkTools

include("CHform.jl")
# include("GSmat.jl")
# include("MetroHast.jl")
include("extrafunctions.jl")

plt = pyimport("matplotlib.pyplot")

println("# threads = ",Threads.nthreads())

Nmax = 40
AmountN = 25
# R = Int64(8)
J = 1.
h = 0.5
# L = Int8(1)
tries = 200

Nlist = unique!(trunc.(Int8, logrange(3., Float64(Nmax), AmountN))) #convert the logrange to ints (without duplicates)

AmountN = length(Nlist)
# #Computing a dummy scalar product to pre-compile the function and have better benchmark results
# ch1tut = CHform(Nlist[1])
# ch2tut = CHform(Nlist[1])
# tut = scalPrdt(ch2tut, ch1tut, J, h)

avrgt = zeros(Float64, AmountN)
errt = zeros(Float64, AmountN)
for l in 1:AmountN
    timelog = zeros(Float64, tries)
    println("N = ", Nlist[l])
    for k in 1:tries
        # println("k = ", k)
        ch1 = CHform(Nlist[l])
        circ1 = randomCirc(Nlist[l], 500)
        update_circ(ch1, circ1)
        
        ch2 = CHform(Nlist[l])
        circ2 = randomCirc(Nlist[l], 500)
        update_circ(ch2, circ2)
        ch2 = conjugateCHf(ch2)

        timelog[k] = @elapsed scalPrdt(ch2, ch1, J, h)
    end
    if l == 1
        avrgt[l] = mean(timelog[2:end])
        errt[l] = std(timelog[2:end])
    else
        avrgt[l] = mean(timelog)
        errt[l] = std(timelog)
    end
end

# avrgt = PyObject(avrgt)
# errt = PyObject(errt)



#---- Save data ----
#select the name
date = now()
fname = string("/home/rezad/Documents/Arbeit/EPFL/Labal/PDM/scripts/Julia /data/TimeH_" , "J=", string(J), " h=", string(h),"_" , string(date) , ".csv")
#creating the file
touch(fname)
#opening the file in write mode
data = open(fname, "w")
#concatenating the data into a matrix
savedata = hcat(Vector{Int8}(Nlist), avrgt, errt)
#writing it to the file
writedlm(fname, savedata, ",")
close(data)


# #Plot using matplotlib

# plt.rc("font", size = 14)
# plt.errorbar(Nlist,avrgt, errt, marker = "x", ls = "", color = "b", zorder = 1) #ls is the acronym of linestyle and "" means no line
# plt.plot(Nlist, avrgt[1]*Nlist.^2/(Nlist[1]^2), color = "b")
# #plotting the fit
# # plt.plot(Rlist[startfit:], pf, zorder = 2) #FITMARK
# # plt.plot(np.exp(Rlist), np.exp(pf))
# #legend and other stuff
# plt.grid(true)
# plt.xscale("log")
# plt.yscale("log")
# plt.legend()
# # else: #FITMARK
# #     plt.legend(["Fit(x) = {:.2e}x³".format(fitcoeffs[0]),"Data"]) #FITMARK
# # #     plt.legend(["Fit(x) = {:.2e}x³+{:.2e}x²+{:.2e}x+{:.2e}".format(fitcoeffs[0],fitcoeffs[1],fitcoeffs[2],fitcoeffs[3]),"Data"]) #FITMARK
# plt.ylabel("Time [s]")
# plt.xlabel("Number of qubits")
# plt.show()
