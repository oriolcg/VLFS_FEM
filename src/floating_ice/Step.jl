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
  ω::Real = 0.2
end

function run_Step(params::Step_params)
  @unpack name, ω = params

  # Fixed parameters
  m = 500
  EI = 1.0e10
  H₀ = 60
  Lb = 300.0

  # Physics
  g = 9.81
  ρ = 1025
  d₀ = m/ρ
  a₁ = EI/ρ

  # wave properties
  f(k) = sqrt(g*k*tanh(k*H₀)) - ω
  k = abs(find_zero(f, 0.2))  # wave number
  λ = 2*π / k                 # wavelength
  @show λ, λ/Lb
  η₀ = 0.01
  ηᵢₙ(x) = η₀*exp(im*k*x[1])
  ϕᵢₙ(x) = -im*(η₀*ω/k)*(cosh(k*(x[2]+0.075*Lb)) / sinh(k*H₀))*exp(im*k*x[1])
  vᵢₙ(x) = (η₀*ω)*(cosh(k*(x[2]+0.075*Lb)) / sinh(k*H₀))*exp(im*k*x[1])
  vzᵢₙ(x) = -im*ω*η₀*exp(im*k*x[1])

  # Numerics constants
  order = 4
  h = Lb/50
  γ = 1.0*order*(order-1)/h
  βₕ = 0.5
  αₕ = -im*ω/g * (1-βₕ)/βₕ

  # Damping
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
  𝒯_Ω = GmshDiscreteModel("models/floating_ice_coarse.msh")

  # Triangulations
  Ω = Interior(𝒯_Ω)
  Γ = Boundary(𝒯_Ω,tags=["beam","free_surface","damping_in","damping_out"])
  Γᵢₙ = Boundary(𝒯_Ω,tags="inlet")
  Γb = Boundary(𝒯_Ω,tags="beam")
  Γd1 = Boundary(𝒯_Ω,tags="damping_in")
  Γd2 = Boundary(𝒯_Ω,tags="damping_out")
  Γf = Boundary(𝒯_Ω,tags="free_surface")
  Γκ = Boundary(𝒯_Ω,tags=["free_surface","damping_in","damping_out"])
  Λb = Skeleton(Γb)

  filename = "data/VTKOutput/floating_ice/Step/"*name
  writevtk(Ω,filename*"_O_trian")
  writevtk(Γb,filename*"_Gb_trian")
  writevtk(Γd1,filename*"_Gd1_trian")
  writevtk(Γd2,filename*"_Gd2_trian")
  writevtk(Γf,filename*"_Gf_trian")
  writevtk(Λb,filename*"_Lb_trian")

  # Measures
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓb = Measure(Γb,degree)
  dΓd1 = Measure(Γd1,degree)
  dΓd2 = Measure(Γd2,degree)
  dΓf = Measure(Γf,degree)
  dΓᵢₙ = Measure(Γᵢₙ,degree)
  dΛb = Measure(Λb,degree)

  # Normals
  nΛb = get_normal_vector(Λb)

  # FE spaces
  reffe = ReferenceFE(lagrangian,Float64,order)
  V_Ω = TestFESpace(Ω, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
  V_Γκ = TestFESpace(Γκ, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
  V_Γη = TestFESpace(Γb, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
  U_Ω = TrialFESpace(V_Ω)
  U_Γκ = TrialFESpace(V_Γκ)
  U_Γη = TrialFESpace(V_Γη)
  X = MultiFieldFESpace([U_Ω,U_Γκ,U_Γη])
  Y = MultiFieldFESpace([V_Ω,V_Γκ,V_Γη])

  # Weak form
  ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
  a((ϕ,κ,η),(w,u,v)) = ∫(  ∇(w)⋅∇(ϕ) )dΩ   +
    ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ )dΓf   +
    ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ - μ₂ᵢₙ*κ*w + μ₁ᵢₙ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd1    +
    ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ - μ₂ₒᵤₜ*κ*w + μ₁ₒᵤₜ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd2    +
    ∫(  ( v*((-ω^2*d₀ + g)*η - im*ω*ϕ) + a₁*Δ(v)*Δ(η) ) +  im*ω*w*η  )dΓb  +
    ∫(  a₁ * ( - jump(∇(v)⋅nΛb) * mean(Δ(η)) - mean(Δ(v)) * jump(∇(η)⋅nΛb) + γ*( jump(∇(v)⋅nΛb) * jump(∇(η)⋅nΛb) ) ) )dΛb
  l((w,u,v)) =  ∫( w*vᵢₙ )dΓᵢₙ - ∫( ηd*w - ∇ₙϕd*(u + αₕ*w) )dΓd1

  op = AffineFEOperator(a,l,X,Y)
  (ϕₕ,κₕ,ηₕ) = Gridap.solve(op)

  xy_cp = get_cell_points(get_fe_dof_basis(V_Γη)).cell_phys_point
  x_cp = [[xy_ij[1] for xy_ij in xy_i] for xy_i in xy_cp]
  η_cdv = get_cell_dof_values(ηₕ)
  p = sortperm(x_cp[1])
  x_cp_sorted = [x_i[p] for x_i in x_cp]
  η_cdv_sorted = [η_i[p] for η_i in η_cdv]
  xs = [(x_i-6*Lb)/Lb for x_i in vcat(x_cp_sorted...)]
  η_rel_xs = [abs(η_i)/η₀ for η_i in vcat(η_cdv_sorted...)]

  writevtk(Γκ,filename*"_kappa",cellfields=["eta_re"=>real(κₕ),"eta_im"=>imag(κₕ)])
  writevtk(Γb,filename*"_eta",cellfields=["eta_re"=>real(ηₕ),"eta_im"=>imag(ηₕ)])
  writevtk(Ω,filename*"_phi",cellfields=["phi_re"=>real(ϕₕ),"phi_im"=>imag(ϕₕ)])

  return (xs,η_rel_xs)

end

end
