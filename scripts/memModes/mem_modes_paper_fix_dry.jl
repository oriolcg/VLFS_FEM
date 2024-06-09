#module Membrane_modes

using Revise
using Gridap
using Plots
using DrWatson
using WaveSpec
using .Constants
using LinearAlgebra
using TickTock
using DataFrames
using Printf



function run_case(mfac = 0.9, tfac = 0.1)

  
  caseName = "ten" * @sprintf("%0.2f", tfac) *"_mass" * @sprintf("%0.2f", mfac)
  name::String = "data/sims_202305/mem_modes_res_fix/dry/mem_modes_"*caseName
  order::Int = 2
  vtk_output::Bool = true
  filename = name*"/mem"

  if( isdir(name) )
    return
  end
  mkdir(name)

  ρw = 1025 #kg/m3 water
  H0 = 10 #m #still-water depth

  # Membrane parameters
  Lm = 2*H0 #m
  @show g #defined in .Constants
  mᵨ = mfac  #mass per unit area of membrane / ρw
  Tᵨ = tfac *g*H0*H0 #T/ρw

  # Excitation wave parameters
  ω = 1.0


  # Domain 
  nx = 40
  LΩ = Lm 
  x₀ = 0.0
  domain =  (x₀, x₀+LΩ)
  @show Lm
  @show LΩ
  @show domain
  println()


  # Mesh
  model = CartesianDiscreteModel(domain,nx)


  # Labelling
  labels_Ω = get_face_labeling(model)
  add_tag_from_tags!(labels_Ω,"mem_bnd",[1,2])   # assign the label "surface" to the entity 3,4 and 6 (top corners and top side)
  

  # Triangulations
  Γ = Interior(model) #same as Triangulation()
  
  # Construct the tag for membrane boundary
  Λmb = Boundary(Γ)
  

  writevtk(model, filename*"_model")
  if vtk_output == true
    writevtk(Γ,filename*"_G")
    writevtk(Λmb,filename*"_Lmb")  
  end


  # Measures
  degree = 2*order
  dΓ = Measure(Γ,degree)
  dΛmb = Measure(Λmb,degree)


  # Normals
  @show nΛmb = get_normal_vector(Λmb)


  # Dirichlet Fnc
  gη(x) = ComplexF64(0.0)

  # FE spaces
  reffe = ReferenceFE(lagrangian,Float64,order)
  V_Γ = TestFESpace(Γ, reffe, conformity=:H1, 
    vector_type=Vector{ComplexF64},
    dirichlet_tags=["mem_bnd"]) #diri
  # V_Γη = TestFESpace(Γη, reffe, conformity=:H1, 
  #   vector_type=Vector{ComplexF64},
  #   dirichlet_tags=["mem_bnd"]) #diri
  # V_Γη = TestFESpace(Γη, reffe, conformity=:H1, 
  #   vector_type=Vector{ComplexF64})
  # U_Γ = TrialFESpace(V_Γ)
  U_Γ = TrialFESpace(V_Γ, gη) #diri
  # U_Γη = TrialFESpace(V_Γη)


  # Weak form: Constant matrices
  ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
  m11(η,v) = ∫( mᵨ*v*η )dΓ
  k11(η,v) = ∫( Tᵨ*∇(v)⋅∇(η) )dΓ +  
              ∫(- Tᵨ*v*∇(η)⋅nΛmb )dΛmb #diri

  l1(v) = ∫( 0*v )dΓ
  println("[MSG] Done Weak form")

  # Global matrices: constant matrices
  M11 = get_matrix(AffineFEOperator( m11, l1, U_Γ, V_Γ ))
  K11 = get_matrix(AffineFEOperator( k11, l1, U_Γ, V_Γ ))
  
  println("[MSG] Done Global matrices")
  println(K11 == transpose(K11))

  #xp = range(xm₀, xm₁, size(V,2)+2)

  λ = LinearAlgebra.eigvals(M11\Matrix(K11))
  V = LinearAlgebra.eigvecs(M11\Matrix(K11))
  @show sum(imag.(λ))
  ωₙ = real.(λ)
  ωₙ = ifelse.(ωₙ .<0, 0, ωₙ)
  @show ωₙ = sqrt.(ωₙ)

  println(length(ωₙ))
  println(V[:,1])

  xp = range(0, Lm, nx*2+1)

  data = Dict(
    "xp" => xp,
    "ωₙ" => ωₙ,
    "V" => V  
  )

  wsave(filename*"_modesdata.jld2", data)
end


mfac = [0.1:0.1:1;]
tfac = [0.05:0.05:1.0;]

# mfac = [0.5]
# tfac = [0.5]


for imfac in mfac
  for itfac in tfac
    run_case(imfac, itfac)
  end
end

#end