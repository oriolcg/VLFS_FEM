using DrWatson
@quickactivate "MonolithicFEMVLFS"

using Plots
using LaTeXStrings
using DataFrames
using DataFramesMeta
using CSV


include(srcdir("floating_ice","Step (2).jl"))

using .Step: Step_params, run_Step

resDir = datadir("floating_ice","Step")
mesh_file = projectdir("models","floating_ice-Step05.json")
name = "Step0.5_"

allparams = Dict(
    "ω" => [0.1,0.2],
    "Q" => 0,     
)

dicts = dict_list(allparams)

function makesim(d::Dict)

  lname = name*savename(d)
  @show lname

  case = Step_params(
    ω = d["ω"],
    Q = d["Q"], 
    Lb = 520,
    Ld = 130,
    xdₒᵤₜ = 650,
    resDir = resDir,
    name = lname,
    mesh_file = mesh_file,
    n_elements = 4)

  x,η = run_Step(case)  

end

for (i, d) in enumerate(dicts)
  result = makesim(d)  
end 

