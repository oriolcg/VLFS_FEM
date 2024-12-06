module Step

using Revise
using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using GridapGmsh
using Parameters
using Roots

export run_Step
export Step_params

@with_kw struct Step_params
  resDir::String = "data/floating_ice/Step/"
  name::String = "Step"
  ω::Real = 0.2
  Q::Float64 = 0.0
  step_ratio::Float64 = 0.5
  mesh_file::String = "floating_ice-Step05.json"
  n_elements::Float64 = 4
  Lₜₒₜ::Float64 = 780
  # Lb::Float64 = 520
  # Ld::Float64 = 130
  # xdₒᵤₜ::Float64 = 650
  kguess::Float64 = 0.5
  μ₀::Float64 = 2.5

end

function run_Step(params::Step_params)

  @unpack name, ω, Q, step_ratio, mesh_file, n_elements, Lₜₒₜ = params
  @unpack resDir = params

  # Fixed parameters
  h_ice = 0.1
  ρ_ice = 917.0
  m = ρ_ice*h_ice
  E = 5.0e9
  ν = 0.33
  I = h_ice^3/12
  EI = E*I/(1-ν^2)
  H₀ = 10.0
  
  # Physics
  g = 9.81
  ρ = 1025
  d₀ = m/(ρ*g)
  a₁ = EI/(ρ*g)
  a₂ = Q*√a₁                 # a₂ => Q = 1.4*√D with D = EI/ρg        

  # # wave properties

  f(k) = √((a₁*k^4 - a₂*k^2 + 1) * g*k*tanh(k*H₀)) - ω   # dispersion relation
  @unpack kguess = params
  @show kguess
  k = abs(find_zero(f, kguess))       # wave number
  λ = 2*π / k                         # wave length

  Ld = 2*λ
  Lb = Lₜₒₜ - 4*λ
  xdₒᵤₜ = Lₜₒₜ - 2*λ

  @show ω, Q, k, λ, λ/Lb, Ld  

  # η₀ = 1
  # ηᵢₙ(x) = η₀*exp(im*k*x[1])
  # ϕᵢₙ(x) = -im*(η₀*ω/k)*(cosh(k*(x[2])) / sinh(k*H₀))*exp(im*k*x[1])
  # vᵢₙ(x) = -1.0*(η₀*ω)*(cosh(k*(x[2])) / sinh(k*H₀))*exp(im*k*x[1])
  # vzᵢₙ(x) = -im*ω*η₀*exp(im*k*x[1])
  # ∇ϕᵢₙ(x) = VectorValue(k*(η₀*ω/k)*(cosh(k*(x[2])) / sinh(k*H₀))*exp(im*k*x[1]), -im*(η₀*ω)*(sinh(k*(x[2])) / sinh(k*H₀))*exp(im*k*x[1]))


  η₀ = 1
  ηᵢₙ(x) = η₀*exp(im*k*x[1])
  # ϕᵢₙ(x) = (-im*η₀*g/ω)*(a₁*k^4 - a₂*k^2 + 1)*(cosh(k*(x[2])) / cosh(k*H₀))*exp(im*k*x[1])
  vᵢₙ(x) = (k*η₀*g/ω)*(a₁*k^4 - a₂*k^2 + 1)*(cosh(k*(x[2])) / cosh(k*H₀))*exp(im*k*x[1])
  vzᵢₙ(x) = (-k*im*η₀*g/ω)*(a₁*k^4 - a₂*k^2 + 1)*(sinh(k*(x[2])) / cosh(k*H₀))*exp(im*k*x[1])
 
  # Numerics constants
  
  order = 4
  h = 1/(n_elements*Lb)
  γ = 100.0*order*(order-1)/h
  βₕ = 0.5
  αₕ = -im*ω/g * (1-βₕ)/βₕ

  # Damping [method 5 (added terms dyn BC and kin BC), ramp function shape 1 - Kim(2014)]
  @unpack μ₀ = params
  @show μ₀

  function μ₁ᵢₙ(x)

      if 0 <= x[1] < Ld
          return μ₀ * (1.0 - sin(π/2*(x[1]) / Ld))
      else
          return 0.0
      end
  end

  function μ₁ₒᵤₜ(x)
    
    if xdₒᵤₜ <= x[1] < Lₜₒₜ
        return μ₀*(1.0 - cos(π/2*(x[1]-xdₒᵤₜ)/Ld))
    else
        return 0.0
    end
  end

  add_term = 0.5*((4*k*g +2) + √((4*k*g +2)^2 - 4))

  μ₂ᵢₙ(x) = μ₁ᵢₙ(x)*k*add_term
  μ₂ₒᵤₜ(x) = μ₁ₒᵤₜ(x)*k*add_term
  ηd(x) = μ₂ᵢₙ(x)*ηᵢₙ(x)
  ∇ₙϕd(x) = μ₁ᵢₙ(x)*vzᵢₙ(x)

  # Fluid model
  𝒯_Ω = DiscreteModelFromFile(mesh_file)
  println("Model loaded")

  # Triangulations
  Ω = Interior(𝒯_Ω)
  Γ = Boundary(𝒯_Ω,tags=["beam","damping_in","damping_out"])
  Γᵢₙ = Boundary(𝒯_Ω,tags="inlet")
  Γb = Boundary(𝒯_Ω,tags=["damping_in","beam","damping_out"])
  # Γd1 = Boundary(𝒯_Ω,tags="damping_in")
  # Γd2 = Boundary(𝒯_Ω,tags="damping_out")
  Λb = Skeleton(Γ)

  
  filename = resDir*name
  writevtk(Ω,filename*"_O_trian.vtu")
  writevtk(Γb,filename*"_Gb_trian.vtu")
  # writevtk(Γd1,filename*"_Gd1_trian.vtu")
  # writevtk(Γd2,filename*"_Gd2_trian.vtu")
  writevtk(Λb,filename*"_Lb_trian.vtu")

  # Measures
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓb = Measure(Γb,degree)
  dΓ = Measure(Γ,degree)
  # dΓd1 = Measure(Γd1,degree)
  # dΓd2 = Measure(Γd2,degree)
  dΓᵢₙ = Measure(Γᵢₙ,degree)
  dΛb = Measure(Λb,degree)

  # Normals
  nΛb = get_normal_vector(Λb)

  # FE spaces
  reffe = ReferenceFE(lagrangian,Float64,order)
  V_Ω = TestFESpace(Ω, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
  V_Γη = TestFESpace(Γ, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
  U_Ω = TrialFESpace(V_Ω)
  U_Γη = TrialFESpace(V_Γη)
  X = MultiFieldFESpace([U_Ω,U_Γη])
  Y = MultiFieldFESpace([V_Ω,V_Γη])


  # Weak form (bending + tensile force)
  ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
  a((ϕ,η),(w,v)) = ∫(  ∇(w)⋅∇(ϕ) )dΩ   +
  ∫(  v*((-ω^2*d₀ + 1)*η - im*ω/g*ϕ) + a₁*Δ(v)*Δ(η) + a₂*∇(v)⋅∇(η) + im*ω*w*η - μ₂ᵢₙ*η*w + (1/g)*(μ₁ᵢₙ*∇ₙ(ϕ)*v) - μ₂ₒᵤₜ*η*w + (1/g)*(μ₁ₒᵤₜ*∇ₙ(ϕ)*v)  )dΓb  +
  ∫(  a₁*( - jump(∇(v)⋅nΛb) * mean(Δ(η)) - mean(Δ(v)) * jump(∇(η)⋅nΛb) + γ*( jump(∇(v)⋅nΛb) * jump(∇(η)⋅nΛb) ) ) )dΛb
  l((w,v)) =  ∫( w*vᵢₙ )dΓᵢₙ - ∫( ηd*w - ∇ₙϕd*v )dΓb


  # solver
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

  # exporting VTK output
  # writevtk(Γκ,filename*"_kappa",cellfields=["eta_re"=>real(κₕ),"eta_im"=>imag(κₕ)])
  writevtk(Γ,filename*"_eta.vtu",cellfields=["eta_re"=>real(ηₕ),"eta_im"=>imag(ηₕ)])
  writevtk(Ω,filename*"_phi.vtu",cellfields=["phi_re"=>real(ϕₕ),"phi_im"=>imag(ϕₕ)])
  writevtk(Γ,filename*"_eta_in.vtu",cellfields=["eta_re"=>x->real(ηᵢₙ(x)),"eta_im"=>x->imag(ηᵢₙ(x))])
  writevtk(Γ,filename*"_vz_in.vtu",cellfields=["vz_re"=>x->real(vzᵢₙ(x)),"vz_im"=>x->imag(vzᵢₙ(x))])

  return (xs,η_rel_xs)

end

end
