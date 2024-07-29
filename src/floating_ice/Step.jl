module Step

using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using GridapGmsh
using Parameters
using Roots

export run_Step
export Step_params

@with_kw struct Step_params
  name::String = "Step"
  k::Real = 0.4
end

function run_Step(params::Step_params)
  @unpack name, k = params

  # Fixed parameters
  h_ice = 0.1
  ρ_ice = 917.0
  m = ρ_ice*h_ice
  E = 5.0e9
  ν = 0.33
  I = h_ice^3/12
  EI = E*I/(1-ν^2)
  H₀ = 10.0             
  Lb = 3.0
  Q = 0.0

  # Physics
  g = 9.81
  ρ = 1025
  d₀ = m/ρ
  a₁ = EI/ρ

  # wave properties
  ω = √((EI*k^4 - Q*k^2 + 1) * g*k*tanh(k*H₀))
  λ = 2*π / k                 # wavelength
  @show λ, λ/Lb
  η₀ = 0.01
  ηᵢₙ(x) = η₀*exp(im*k*x[1])
  ϕᵢₙ(x) = -im*(η₀*ω/k)*(cosh(k*(x[2])) / sinh(k*H₀))*exp(im*k*x[1])
  vᵢₙ(x) = (η₀*ω)*(cosh(k*(x[2]+0.075*Lb)) / sinh(k*H₀))*exp(im*k*x[1])
  vzᵢₙ(x) = -im*ω*η₀*exp(im*k*x[1])

  # Numerics constants
  order = 4
  h = Lb/50
  γ = 1.0*order*(order-1)/h
  βₕ = 0.5
  αₕ = -im*ω/g * (1-βₕ)/βₕ

  # Damping [method 5 (added terms dyn BC and kin BC), ramp function shape 1 - Kim(2014)]
  μ₀ = 6.0
  Ld = 4*Lb
  xdₒᵤₜ = 9*Lb
  μ₁ᵢₙ(x) = μ₀*(1.0 - sin(π/2*(x[1])/Ld))
  μ₁ₒᵤₜ(x) = μ₀*(1.0 - cos(π/2*(x[1]-xdₒᵤₜ)/Ld))
  μ₂ᵢₙ(x) = μ₁ᵢₙ(x)*k
  μ₂ₒᵤₜ(x) = μ₁ₒᵤₜ(x)*k
  ηd(x) = μ₂ᵢₙ(x)*ηᵢₙ(x)
  ∇ₙϕd(x) = μ₁ᵢₙ(x)*vzᵢₙ(x)

  # Fluid model
  𝒯_Ω = DiscreteModelFromFile("models/floating_ice_coarse.json")
  println("Model loaded")

  # Triangulations
  Ω = Interior(𝒯_Ω)
  Γ = Boundary(𝒯_Ω,tags=["beam","damping_in","damping_out"])
  Γᵢₙ = Boundary(𝒯_Ω,tags="inlet")
  Γb = Boundary(𝒯_Ω,tags="beam")
  Γd1 = Boundary(𝒯_Ω,tags="damping_in")
  Γd2 = Boundary(𝒯_Ω,tags="damping_out")
  # Γκ = Boundary(𝒯_Ω,tags=["damping_in","damping_out"])
  Λb = Skeleton(Γ)

  filename = "data/VTKOutput/floating_ice/Step/"*name
  writevtk(Ω,filename*"_O_trian")
  writevtk(Γb,filename*"_Gb_trian")
  writevtk(Γd1,filename*"_Gd1_trian")
  writevtk(Γd2,filename*"_Gd2_trian")
  writevtk(Λb,filename*"_Lb_trian")

  # Measures
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓb = Measure(Γb,degree)
  dΓd1 = Measure(Γd1,degree)
  dΓd2 = Measure(Γd2,degree)
  dΓᵢₙ = Measure(Γᵢₙ,degree)
  dΛb = Measure(Λb,degree)

  # Normals
  nΛb = get_normal_vector(Λb)

  # FE spaces
  reffe = ReferenceFE(lagrangian,Float64,order)
  V_Ω = TestFESpace(Ω, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
  # V_Γκ = TestFESpace(Γκ, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
  V_Γη = TestFESpace(Γ, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
  U_Ω = TrialFESpace(V_Ω)
  # U_Γκ = TrialFESpace(V_Γκ)
  U_Γη = TrialFESpace(V_Γη)
  X = MultiFieldFESpace([U_Ω,U_Γη])
  Y = MultiFieldFESpace([V_Ω,V_Γη])

  # Weak form
  ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
  a((ϕ,η),(w,v)) = ∫(  ∇(w)⋅∇(ϕ) )dΩ   +
    ∫(  v*((-ω^2*d₀ + g)*η - im*ω*ϕ) + a₁*Δ(v)*Δ(η) + im*ω*w*η - μ₂ᵢₙ*η*w + μ₁ᵢₙ*∇ₙ(ϕ)*v )dΓd1    +
    ∫(  v*((-ω^2*d₀ + g)*η - im*ω*ϕ) + a₁*Δ(v)*Δ(η) + im*ω*w*η - μ₂ₒᵤₜ*η*w + μ₁ₒᵤₜ*∇ₙ(ϕ)*v )dΓd2   +
    ∫(( v*((-ω^2*d₀ + g)*η - im*ω*ϕ) + a₁*Δ(v)*Δ(η) ) +  im*ω*w*η  )dΓb  +
    ∫(  a₁ * ( - jump(∇(v)⋅nΛb) * mean(Δ(η)) - mean(Δ(v)) * jump(∇(η)⋅nΛb) + γ*( jump(∇(v)⋅nΛb) * jump(∇(η)⋅nΛb) ) ) )dΛb
  l((w,v)) =  ∫( w*vᵢₙ )dΓᵢₙ - ∫( ηd*w - ∇ₙϕd*v )dΓd1


  # # Weak form (bending + tensile force)
  # ## d₀ = m/ρ,  a₁ = EI/ρ,  a2 = Q/ρ, 
  # ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
  # a((ϕ,η),(w,v)) = ∫(  ∇(w)⋅∇(ϕ) )dΩ   +
  #   ∫(  v*((-ω^2*d₀ + g)*η - im*ω*ϕ) + (a₁)*Δ(v)*Δ(η) + Tᵨ*∇(v)⋅∇(η) + im*ω*w*η - μ₂ᵢₙ*η*w + μ₁ᵢₙ*∇ₙ(ϕ)*v )dΓd1    +
  #   ∫(  v*((-ω^2*d₀ + g)*η - im*ω*ϕ) + a₁*Δ(v)*Δ(η) + Tᵨ*∇(v)⋅∇(η) + im*ω*w*η - μ₂ₒᵤₜ*η*w + μ₁ₒᵤₜ*∇ₙ(ϕ)*v )dΓd2   +
  #   ∫(( v*((-ω^2*d₀ + g)*η - im*ω*ϕ) + a₁*Δ(v)*Δ(η) )+ Tᵨ*∇(v)⋅∇(η) +  im*ω*w*η  )dΓb  +
  #   ∫(  a₁ * ( - jump(∇(v)⋅nΛb) * mean(Δ(η)) - mean(Δ(v)) * jump(∇(η)⋅nΛb) + γ*( jump(∇(v)⋅nΛb) * jump(∇(η)⋅nΛb) ) ) )dΛb +
  #  
  # l((w,v)) =  ∫( w*vᵢₙ )dΓᵢₙ - ∫( ηd*w - ∇ₙϕd*v )dΓd1

  op = AffineFEOperator(a,l,X,Y)
  println("Operator created")
  (ϕₕ,ηₕ) = Gridap.solve(op)
  println("Operator solved")

  xy_cp = get_cell_points(get_fe_dof_basis(V_Γη)).cell_phys_point            
  x_cp = [[xy_ij[1] for xy_ij in xy_i] for xy_i in xy_cp]
  η_cdv = get_cell_dof_values(ηₕ)
  p = sortperm(x_cp[1])
  x_cp_sorted = [x_i[p] for x_i in x_cp]
  η_cdv_sorted = [η_i[p] for η_i in η_cdv]
  xs = [(x_i-6*Lb)/Lb for x_i in vcat(x_cp_sorted...)]
  η_rel_xs = [abs(η_i)/η₀ for η_i in vcat(η_cdv_sorted...)]

  # writevtk(Γκ,filename*"_kappa",cellfields=["eta_re"=>real(κₕ),"eta_im"=>imag(κₕ)])
  writevtk(Γ,filename*"_eta",cellfields=["eta_re"=>real(ηₕ),"eta_im"=>imag(ηₕ)])
  writevtk(Ω,filename*"_phi",cellfields=["phi_re"=>real(ϕₕ),"phi_im"=>imag(ϕₕ)])

  return (xs,η_rel_xs)

end

end
