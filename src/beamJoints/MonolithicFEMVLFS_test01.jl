#module MonolithicFEMVLFS

using DrWatson
@quickactivate "MonolithicFEMVLFS"

using Plots
using LaTeXStrings
using DataFrames
using DataFramesMeta
using CSV

export run_tests

# Include source files
include("Periodic_Beam.jl")
include("Periodic_Beam_FS.jl")
include("Khabakhpasheva_freq_domain.jl")
include("Khabakhpasheva_time_domain.jl")
include("Liu.jl")
include("Yago_freq_domain.jl")
include("Multi_geo_freq_domain.jl")

using .Periodic_Beam: Periodic_Beam_params, run_periodic_beam
using .Periodic_Beam_FS: Periodic_Beam_FS_params, run_periodic_beam_FS
using .Khabakhpasheva_freq_domain: Khabakhpasheva_freq_domain_params, run_Khabakhpasheva_freq_domain
using .Khabakhpasheva_time_domain: Khabakhpasheva_time_domain_params, run_Khabakhpasheva_time_domain
using .Liu: Liu_params, run_Liu
using .Yago_freq_domain: Yago_freq_domain_params, run_Yago_freq_domain
using .MultiGeo_freq_domain: MultiGeo_freq_domain_params, run_MultiGeo_freq_domain

# Extend DrWatson functions
DrWatson.allaccess(c::Periodic_Beam_params) = (:n, :dt, :tf, :orderϕ, :orderη, :k)
DrWatson.default_prefix(c::Periodic_Beam_params) = c.name
DrWatson.allaccess(c::Periodic_Beam_FS_params) = (:n, :dt, :tf, :order, :k)
DrWatson.default_prefix(c::Periodic_Beam_FS_params) = c.name
DrWatson.allaccess(c::Khabakhpasheva_freq_domain_params) = (:nx, :ny, :order, :ξ, :vtk_output)
DrWatson.default_prefix(c::Khabakhpasheva_freq_domain_params) = c.name
DrWatson.allaccess(c::Khabakhpasheva_time_domain_params) = (:nx, :ny, :order, :ξ, :vtk_output)
DrWatson.default_prefix(c::Khabakhpasheva_time_domain_params) = c.name
DrWatson.allaccess(c::Liu_params) = (:ω,)
DrWatson.default_prefix(c::Liu_params) = c.name
DrWatson.allaccess(c::Yago_freq_domain_params) = (:nx, :ny, :nz, :order, :λfactor, :dfactor)
DrWatson.default_prefix(c::Yago_freq_domain_params) = c.name
DrWatson.allaccess(c::MultiGeo_freq_domain_params) = (:mesh_file, :order, :dfactor)
DrWatson.default_prefix(c::MultiGeo_freq_domain_params) = c.name

# Include script files
include("../scripts/5-1-1-periodic-beam-spatial-convergence.jl")
include("../scripts/5-1-2-periodic-beam-time-convergence.jl")
include("../scripts/5-1-3-periodic-beam-energy.jl")
include("../scripts/5-1-4-periodic-beam-free-surface-energy.jl")
include("../scripts/5-2-1-Khabakhpasheva-freq-domain.jl")
include("../scripts/5-2-2-Khabakhpasheva-time-domain.jl")
include("../scripts/5-3-1-Liu.jl")
include("../scripts/5-4-1-Yago.jl")
include("../scripts/5-5-1-MultiGeo.jl")

# Define Execution function
function run_5_1_1(case::Periodic_Beam_params)
  case_name = savename(case)
  println("-------------")
  println("Case: ",case_name)
  e_ϕ, e_η, = run_periodic_beam(case)
  e_ϕ_i = last(e_ϕ)
  e_η_i = last(e_η)
  case_name_suffix = savename(case,"jld2";digits=10)    
  file = datadir("5-1-1-periodic-beam-spatial-convergence", case_name_suffix)
  prefix, data, suffix = DrWatson.parse_savename(case_name_suffix, parsetypes=(Int, Float64))    
  push!(data, "e_ϕ_i"=>e_ϕ_i, "e_η_i"=>e_η_i)
  #@tagsave(file,data)
  return data
end

# Warm-up case
k = 1
H = 1.0
g = 9.81
ω = √(g*k*tanh(k*H))
T = 2π/ω
Δt = T/1
n = 3
order = 2
path = datadir("5-1-1-periodic-beam-spatial-convergence")
case = Periodic_Beam_params(
  name="11Warm-up",
  n=n,
  dt=Δt, 
  tf=T,
  k=k,
  orderϕ=order,
  orderη=order,
  vtk_output=true
)
filename = config -> savename(config; digits=8)
produce_or_load(path,case,run_5_1_1; prefix="", filename = filename)

#end
