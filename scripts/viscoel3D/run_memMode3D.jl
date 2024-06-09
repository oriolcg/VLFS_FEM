using DrWatson
using Parameters
using WaveSpec
using .Constants
@quickactivate "MonolithicFEMVLFS.jl"

# Here you may include files from the source directory
include(srcdir("viscoel3D","mem3D_modes_para_free_iter.jl"))

lresDir::String = "data/sims_202309/run/mem3D_modes_res_free"

# Warm-up run
mfac = [0.9]
tfac = [0.1]
params = Memb3D_modes.Memb_params(
  resDir = lresDir)
Memb3D_modes.main(params, mfac, tfac)

# Production run
@with_kw struct run_params
  resDir = lresDir
  order::Int = 1
  order2::Int = 2
  vtk_output::Bool = true

  H0 = 10 #m #still-water depth

  # Membrane parameters
  Lm = 2*H0 #m
  Wm = Lm  

  nωₙ = 7
  errLim = 1e-3
  
  # Domain 
  nx = 60
  ny = 60
  nz = 8
  mesh_rz = 1.2 #Ratio for Geometric progression of eleSize
  LΩ = 6*H0 
  WΩ = 6*H0
  x₀ = 0.0
  domain =  (x₀, x₀+LΩ, -WΩ/2, WΩ/2 , -H0, 0.0)
  partition = (nx, ny, nz)
  xm₀ = x₀ + 2*H0
  xm₁ = xm₀ + Lm
  ym₀ = -Wm/2
  ym₁ = Wm/2

end
params = run_params()
Memb3D_modes.main(params, mfac, tfac)
