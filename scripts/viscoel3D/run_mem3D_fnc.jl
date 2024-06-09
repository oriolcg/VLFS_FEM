using DrWatson
using Parameters
using WaveSpec
using .Constants
@quickactivate "MonolithicFEMVLFS.jl"

# Here you may include files from the source directory
include(srcdir("viscoel3D","mem3D_freq_para_damp_fnc.jl"))

resDir::String = "data/sims_202309/run3D/mem3D_m90_t025_tau10_r1"

# Warm-up run
params = Memb3D.Memb_params(name = resDir)
Memb3D.main(params)

# Production run
@with_kw struct run_params
  name = resDir
  order::Int = 2
  vtk_output::Bool = true

  H0 = 10 #m #still-water depth

  # Wave parameters
  ω = 1.5:0.05:3.5
  T = 2*π./ω
  η₀ = 0.10*ones(length(ω))
  # α = randomPhase1(ω; seed=100)
  # k = dispersionRelAng.(H0, ω; msg=false)

  # Membrane parameters
  Lm = 2*H0 #m
  Wm = Lm  
  mᵨ = 0.9 #mass per unit area of membrane / ρw
  Tᵨ = 0.1/4*g*Lm*Lm #T/ρw
  τ = 0.1#damping coeff


  # Domain 
  nx = 220
  ny = 60
  nz = 10
  mesh_rz = 1.2 #Ratio for Geometric progression of eleSize
  Ld = 6*H0 #damping zone length
  LΩ = 5*Lm + 2*Ld
  WΩ = 6*H0
  x₀ = -Ld
  domain =  (x₀, x₀+LΩ, -WΩ/2, WΩ/2 , -H0, 0.0)
  partition = (nx, ny, nz)
  xdᵢₙ = 0.0
  xdₒₜ = x₀ + LΩ - Ld
  xm₀ = xdᵢₙ + 2*Lm
  xm₁ = xm₀ + Lm
  ym₀ = -Wm/2
  ym₁ = Wm/2

end
params = run_params()
Memb3D.main(params)
