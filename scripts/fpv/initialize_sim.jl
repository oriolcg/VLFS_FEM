
using DrWatson
@quickactivate "MonolithicFEMVLFS"

include(srcdir("fpv","beam_no_joint_src.jl"))

using WaveSpec
using Plots
using .Constants
using .Jonswap
using .WaveTimeSeries
using DataFrames
using TimerOutputs
using TickTock 

using .beam

# #include(srcdir("beam_no_joint_src.jl"))
# function constant_spectrum(start, finish,threshold) 
#     ω_low = collect(start:0.1:threshold)
#     ω_high = collect(threshold:0.1:finish)[2:end]
#     ω = vcat(ω_low,ω_high)
#     return ω
# end

# #Set the parameters
# h_b = [1]

# #constant wave spectrum 
# ω = constant_spectrum(1,4,3)
# α = vec(repeat([0], length(ω),1))
# A_w = vec(repeat([1], length(ω),1))


function constant_spectrum(start, finish,threshold) 
    ω_low = collect(start:0.02:threshold)
    ω_high = collect(threshold:0.1:finish)[2:end]
    ω = vcat(ω_low,ω_high)
    return ω
end

#Set the parameters
h_b = [1]

#constant wave spectrum 
ω = constant_spectrum(0.2,10,3)
α = vec(repeat([0], length(ω),1))
A_w = vec(repeat([1], length(ω),1))


#Jonswap spectrum 
#= 
Tp = 3.
Hs = 5
d = 30

ω, S, A = jonswap(Hs, Tp; 
    plotflag=false,  nω=121, ωc = 12) #omega cut off should be 10. 

ω = ω[2:end]
α = vec(repeat([0], length(ω),1))
A_w = A[2:end]    
jonswap_plot = plot(ω,A_w, label = false, xlabel = "Frequency [rad/s]",
ylabel = "Amplitude [m]")
jonswap_plot = vline!([1], label = "1 rad/s")
 =#
#Collect all parameters to a dictionaryc 

allparams = Dict(
    "length_beam" => [200],
    "h_b" => h_b,                            
    # "material" => ["eps", "cfrp", "gfrp"],# ["eps", "cfrp", "gfrp", "neopren", "pvc", "hdpe", "steel"], #store youngs modulus and density 
    "material" => ["eps"],
    "phase"=> [α],
    "omega" => [ω],
    "amplitude" => [A_w],
    "mesh_size" => [1.0],
    "numFloat" => [1]
)

#creates all possible combinations
dicts = dict_list(allparams)
tick()


function makesim(d::Dict)
    
    RAO_η, RAO_ηx, η_ϕ, ηₓ_ϕ, RAO_ηxx, RAO_ηxxx = beam.run_beam(d)
    fulld = copy(d)
    fulld["RAO_η"] = RAO_η
    fulld["RAO_ηx"] = RAO_ηx
    fulld["RAO_ηxx"] = RAO_ηxx
    fulld["RAO_ηxxx"] = RAO_ηxxx
    fulld["η_ϕ"] = η_ϕ
    fulld["ηₓ_ϕ"] = ηₓ_ϕ
    #fulld["solution_η"] = solution_η
    return fulld
end


for (i, d) in enumerate(dicts)
    tick()
    result = makesim(d) # is a dict conainting the input and output of simulation
    #safesave("data/sims2011/"*savename(d, "jld2"), result)
    safesave(datadir("fpv_202312", savename(d, "jld2")), result)
    tock()
end 

tock()



