using DrWatson
using Parameters
using WaveSpec
using .Constants
@quickactivate "MonolithicFEMVLFS.jl"

# Here you may include files from the source directory
include(srcdir("viscoel3D","mem3D_freq_para_damp.jl"))

resDir::String = "data/sims_202309/run/mono3D_freq_free"

# Warm-up run
params = Memb3D.Memb_params(
    name = resDir)
Memb3D.main(params)

# Production run
@with_kw struct run_params
  name = resDir
  order::Int = 2
  vtk_output::Bool = true

  H0 = 10 #m #still-water depth

  # Membrane parameters
  Lm = 2*H0 #m
  Wm = Lm  
  mᵨ = 0.9 #mass per unit area of membrane / ρw
  Tᵨ = 0.1/4*g*Lm*Lm #T/ρw
  τ = 0.0#damping coeff
  diriFlag = false

  # Wave parameters
  ω = 1.5#3.45#2.0#2.4
  η₀ = 0.10  

  # Domain 
  nx = 200
  ny = 20
  nz = 10
  mesh_rz = 1.2 #Ratio for Geometric progression of eleSize
  # Ld = 9*H0 #damping zone length
  # LΩ = 5*Lm + 2*Ld
  Ld = 6*H0 #damping zone length
  LΩ = 4*Lm + 2*Ld
  WΩ = 4*H0
  x₀ = -Ld
  domain =  (x₀, x₀+LΩ, -WΩ/2, WΩ/2 , -H0, 0.0)
  partition = (nx, ny, nz)
  xdᵢₙ = 0.0
  xdₒₜ = x₀ + LΩ - Ld
  xm₀ = xdᵢₙ + 1.5*Lm
  xm₁ = xm₀ + Lm
  ym₀ = -Wm/2
  ym₁ = Wm/2
end
params = run_params()
Memb3D.main(params)
