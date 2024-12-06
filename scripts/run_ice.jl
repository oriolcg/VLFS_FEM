using DrWatson
@quickactivate "MonolithicFEMVLFS"

using Revise
using Plots
using LaTeXStrings
using DataFrames
using DataFramesMeta
using CSV


include(srcdir("floating_ice","Step (3).jl"))

# using .Step: Step_params, run_Step  sqrt(eta_im^2+eta_re^2)

resDir = datadir("floating_ice","Step")
mesh_file = projectdir("models","floating_ice-Step05.json")
name = "0.5_varLd_k3"

allparams = Dict(
    "ω" => [0.8],
    "Q" => [1.9],     
)

dicts = dict_list(allparams)

function makesim(d::Dict)

  lname = name*savename(d)
  @show lname

  case = Step.Step_params(
    ω = d["ω"],
    Q = d["Q"],
    Lₜₒₜ = 780,
    # Lb = 520,
    # Ld = 130,
    # xdₒᵤₜ = 650,
    resDir = resDir,
    name = lname,
    mesh_file = mesh_file,
    n_elements = 4, 
    kguess = 0.4,      #k_guesses for Q=1.9 ω=0.8: k=[0.09, 0.29, 0.4]
    μ₀ = 50)

  x,η = Step.run_Step(case)  

end

# # Warmup case
# case_warmup = Step.Step_params(
#   ω = 0.4,
#   Q = 0,
#   Lₜₒₜ = 780, 
#   # Lb = 520,
#   # Ld = 130,
#   # xdₒᵤₜ = 650,
#   resDir = resDir,
#   name = "tmp",
#   mesh_file = projectdir("models","floating_ice-noStep_coarse.json"),
#   n_elements = 1)
# dummy_result = Step.run_Step(case_warmup)

for (i, d) in enumerate(dicts)
  result = makesim(d)  
end 

