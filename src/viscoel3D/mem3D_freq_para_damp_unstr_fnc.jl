module Memb3D

using Revise
using Gridap
using Parameters
using Plots
using DrWatson
using WaveSpec
using .Constants
using TickTock
using DataFrames
using Printf
using GridapGmsh


function main(params)

  function run_freq(ω, η₀)

    tick()

    # Wave parameters
    k = dispersionRelAng(H0, ω; msg = false)
    λ = 2*π/k
    T = 2π/ω
    ηᵢₙ(x) = η₀*exp(im*k*(x[1]-x₀))
    ϕᵢₙ(x) = -im*(η₀*ω/k)*(cosh(k*x[3]) / 
      sinh(k*H0))*exp(im*k*(x[1]-x₀))
    vxᵢₙ(x) = (η₀*ω)*(cosh(k*x[3]) / 
      sinh(k*H0))*exp(im*k*(x[1]-x₀))
    vzfsᵢₙ(x) = -im*ω*η₀*exp(im*k*(x[1]-x₀)) #???
    @show ω, T
    @show λ
    @show η₀
    @show H0, H0/λ

    # Numeric constants
    αₕ = -im*ω/g * (1-βₕ)/βₕ
    @show αₕ
    println()

    # Damping
    μ₀ = 2.5
    μ₁ᵢₙ(x) = μ₀*(1.0 - sin(π/2*(x[1]-x₀)/Ld))
    μ₁ₒᵤₜ(x) = μ₀*(1.0 - cos(π/2*(x[1]-xdₒₜ)/Ld))
    μ₂ᵢₙ(x) = μ₁ᵢₙ(x)*k
    μ₂ₒᵤₜ(x) = μ₁ₒᵤₜ(x)*k
    ηd(x) = μ₂ᵢₙ(x)*ηᵢₙ(x)
    ∇ₙϕd(x) = μ₁ᵢₙ(x)*vzfsᵢₙ(x) #???


    # Weak form
    ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,0.0,1.0)  
    a((ϕ,κ,η),(w,u,v)) =      
      ∫(  ∇(w)⋅∇(ϕ) )dΩ   +
      ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ )dΓfs   +
      ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ 
        - μ₂ᵢₙ*κ*w + μ₁ᵢₙ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd1    +
      # ∫( -w * im * k * ϕ )dΓot +
      ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ 
        - μ₂ₒᵤₜ*κ*w + μ₁ₒᵤₜ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd2    +
      ∫(  v*(g*η - im*ω*ϕ) +  im*ω*w*η
        - mᵨ*v*ω^2*η + Tᵨ*(1-im*ω*τ)*∇(v)⋅∇(η) )dΓm  #+ 
      #∫(- Tᵨ*(1-im*ω*τ)*v*∇(η)⋅nΛmb )dΛmb #diri

    l((w,u,v)) =  ∫( w*vxᵢₙ )dΓin - ∫( ηd*w - ∇ₙϕd*(u + αₕ*w) )dΓd1


    println("[MSG] Formulation complete")

    # Solution
    op = AffineFEOperator(a,l,X,Y)
    (ϕₕ,κₕ,ηₕ) = solve(op)
    
    # xΓκ = get_cell_coordinates(Γκ)  
      
    # Energy flux (Power) calculation
    ηxyz = ∇(ηₕ) 
    ηx = ηxyz ⋅ VectorValue(1.0, 0.0, 0.0)
    ηy = ηxyz ⋅ VectorValue(0.0, 1.0, 0.0)
    Pd = sum(∫( abs(ηx)*abs(ηx) + abs(ηy)*abs(ηy) )dΓm)
    Pd = 0.5*Tᵨ*ρw*τ*ω*ω*Pd

    # Wave energy flux
    kh = k*H0
    wave_n = 0.5*(1 + 2*kh/sinh(2*kh))
    Pin = (0.5*ρw*g*η₀*η₀)*(ω/k)*wave_n * Pin_Wm    

    push!(prbPow, [ω, Pin, Pd])

    println("Power In \t ",Pin," W")    
    println("Power Abs \t ",Pd," W")

    # Function for inlet phase
    κin = interpolate_everywhere(ηᵢₙ, 
      FESpace(Γκ, reffe, conformity=:H1, 
        vector_type=Vector{ComplexF64}))

    κr = κₕ - κin

    # Wave direction
    κₕᵢ = ∇(κₕ)⋅VectorValue(1.0, 0.0, 0.0)
    κₕⱼ = ∇(κₕ)⋅VectorValue(0.0, 1.0, 0.0)
    cx = im*ω*κₕ * κₕᵢ / ( (κₕᵢ*κₕᵢ) + (κₕⱼ*κₕⱼ) )
    cy = im*ω*κₕ * κₕⱼ / ( (κₕᵢ*κₕᵢ) + (κₕⱼ*κₕⱼ) )  
    
    
    # Pressure calculation
    Pres_d = ρw * im * ω * ϕₕ
    FCyl = sum(∫( Pres_d * nΓCyl )dΓcyl)
    @show FCyl
    push!(prbForce, [ω, FCyl[1], FCyl[2], FCyl[3]])

    @printf( dataFile, "%15.6f, %15.6f, %15.6f, %15.6f, %15.6f, %15.6f \n", 
      ω, Pin, Pd, abs(FCyl[1]), abs(FCyl[2]), abs(FCyl[3]))    
    flush(dataFile)
    
    paraFolder = name*"/mem_omg_"*@sprintf("%.2f",ω)
    paraFile = paraFolder*"/mem"

    if vtk_output == true
      if( isdir(paraFolder) )
        rm(paraFolder; recursive=true)
      end
      mkdir(paraFolder)

      writevtk(Ω,paraFile * "_O_sol.vtu",
        cellfields = ["phi_re" => real(ϕₕ),"phi_im" => imag(ϕₕ),
        "phi_abs" => abs(ϕₕ), "phi_ang" => angle∘(ϕₕ)])
      writevtk(Γκ,paraFile * "_Gk_sol.vtu",
        cellfields = ["eta_re" => real(κₕ),"eta_im" => imag(κₕ),
        "eta_abs" => abs(κₕ), "eta_ang" => angle∘(κₕ),
        "etaR_re" => real(κr),"etaR_im" => imag(κr),
        "etaR_abs" => abs(κr), "etaR_ang" => angle∘(κr),
        "ηin_abs" => abs(κin), "ηin_ang" => angle∘(κin),
        "cx" => real(cx), "cy" => real(cy),
        "cx_im" => imag(cx), "cy_im" => imag(cy)])
      writevtk(Γη,paraFile * "_Ge_sol.vtu",
        cellfields = ["eta_re" => real(ηₕ),"eta_im" => imag(ηₕ),
        "eta_abs" => abs(ηₕ), "eta_ang" => angle∘(ηₕ),
        "eta_x_re" => real(ηx), "eta_x_im" => imag(ηx),
        "eta_y_re" => real(ηy), "eta_y_im" => imag(ηy)])
    end    

    tock()
    println("----------")
    println()
    return 0
  end

  @unpack name, order, vtk_output = params
  filename = name*"/mem"

  ρw = 1025 #kg/m3 water
  @unpack H0 = params #m #still-water depth
  @show H0

  # Membrane parameters
  @unpack mᵨ, Tᵨ, τ, Pin_Wm = params
  @show g   #defined in .Constants
  @show mᵨ  #mass per unit area of membrane / ρw
  @show Tᵨ  #T/ρw
  @show τ   #damping coeff
  @show Pin_Wm 
  # Projected Width of membrane for calculation of input power

  # Wave parameters
  @unpack ω, η₀ = params
  @show ω 
  @show η₀  
  @show H0
  println()


  # Domain 
  @unpack mesh_file = params
  @unpack Ld, LΩ, WΩ, x₀ = params
  @unpack xdᵢₙ, xdₒₜ = params
  @show mesh_file  
  @show LΩ, WΩ
  println()

  # Numeric constants
  # h = LΩ / nx
  # γ = 1.0*order*(order-1)/h
  βₕ = 0.5
  # @show h
  @show βₕ  
  println()

  model = DiscreteModelFromFile(mesh_file)  

  println("Num of nodes = ",length(model.grid.node_coordinates))

  # # Triangulations
  Ω = Interior(model) #same as Triangulation()
  Γ = Boundary(model, tags=["FreeSurface","Membrane"]) #same as BoundaryTriangulation()
  Γin = Boundary(model, tags="Inlet")
  Γot = Boundary(model, tags="Outlet")
  Γcyl = Boundary(model, tags="Cylinder")


  # Auxiliar functions
  function is_damping1(xs) # Check if an element is inside the damping zone 1
    n = length(xs)
    x = (1/n)*sum(xs)
    (x₀ <= x[1] <= xdᵢₙ ) * ( x[3] ≈ H0)
  end
  function is_damping2(xs) # Check if an element is inside the damping zone 2
    n = length(xs)
    x = (1/n)*sum(xs)
    (xdₒₜ <= x[1] ) * ( x[3] ≈ H0)
  end

  # Masking and Beam Triangulation
  Γκ = Boundary(model, tags = "FreeSurface")
  xΓκ = get_cell_coordinates(Γκ)
  Γd1_to_Γ_mask = lazy_map(is_damping1, xΓκ)
  Γd2_to_Γ_mask = lazy_map(is_damping2, xΓκ)
  Γd1 = Triangulation(Γκ, findall(Γd1_to_Γ_mask))
  Γd2 = Triangulation(Γκ, findall(Γd2_to_Γ_mask))  
  Γm = Boundary(model, tags = "Membrane")
  Γη = Γm  
  Γfs = Triangulation(Γκ, findall(!, Γd1_to_Γ_mask .| Γd2_to_Γ_mask))


  # # Construct the tag for membrane boundary
  # Λmb = Boundary(Γm)
  # xΛmb = get_cell_coordinates(Λmb)
  # xΛmb_n1 = findall(model.grid_topology.vertex_coordinates .== xΛmb[1])
  # xΛmb_n2 = findall(model.grid_topology.vertex_coordinates .== xΛmb[2])
  # new_entity = num_entities(labels_Ω) + 1
  # labels_Ω.d_to_dface_to_entity[1][xΛmb_n1[1]] = new_entity
  # labels_Ω.d_to_dface_to_entity[1][xΛmb_n2[1]] = new_entity
  # add_tag!(labels_Ω, "mem_bnd", [new_entity])


  writevtk(model, filename*"_model")
  if vtk_output == true
    writevtk(Ω,filename*"_O")
    writevtk(Γ,filename*"_G")
    writevtk(Γm,filename*"_Gm")  
    writevtk(Γd1,filename*"_Gd1")
    writevtk(Γd2,filename*"_Gd2")
    writevtk(Γfs,filename*"_Gfs")
    writevtk(Γcyl,filename*"_Gcyl")
    # writevtk(Λmb,filename*"_Lmb")  
  end


  # Measures
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓm = Measure(Γm,degree)
  dΓd1 = Measure(Γd1,degree)
  dΓd2 = Measure(Γd2,degree)
  dΓfs = Measure(Γfs,degree)
  dΓin = Measure(Γin,degree)
  dΓot = Measure(Γot,degree)
  dΓcyl = Measure(Γcyl,degree)
  # dΛmb = Measure(Λmb,degree)


  # # Normals
  # @show nΛmb = get_normal_vector(Λmb)
  nΓCyl = get_normal_vector(Γcyl)
  @show nΓCyl


  # Dirichlet Fnc
  gη(x) = ComplexF64(0.0)

  # FE spaces
  reffe = ReferenceFE(lagrangian,Float64,order)
  V_Ω = TestFESpace(Ω, reffe, conformity=:H1, 
    vector_type=Vector{ComplexF64})
  V_Γκ = TestFESpace(Γκ, reffe, conformity=:H1, 
    vector_type=Vector{ComplexF64})  
  V_Γη = TestFESpace(Γη, reffe, conformity=:H1, 
      vector_type=Vector{ComplexF64})
  U_Ω = TrialFESpace(V_Ω)
  U_Γκ = TrialFESpace(V_Γκ)
  U_Γη = TrialFESpace(V_Γη)
  X = MultiFieldFESpace([U_Ω,U_Γκ,U_Γη])
  Y = MultiFieldFESpace([V_Ω,V_Γκ,V_Γη])  

  # Power probes
  prbPow = DataFrame(zeros(Float64, 1, 3), :auto)
  prbForce = DataFrame(zeros(ComplexF64, 1, 4), :auto)
  dataFile = open(filename*"_printed_data.txt", "w")

  # Run weak-form for each freq
  run_freq.(ω, η₀)

  prbPow = prbPow[2:end,:]
  prbForce = prbForce[2:end,:]

  k = dispersionRelAng.(H0, ω; msg=false)
  data = Dict(
    "ω" => ω,
    "η₀" => η₀,
    "k" => k,
    "prbPow" => prbPow,
    "prbForce" => prbForce)

  wsave(filename*"_data.jld2", data)
  close(dataFile)
  return (prbPow, prbForce)
end


"""
Memb_params

Parameters for the VIV.jl module.
"""
@with_kw struct Memb_params

  mesh_file::String = "data/sims_202311/northHoyle/mesh/mesh0/meshWorks.msh"
  name::String = "data/sims_202309/run/mono3D_freq_free"
  order::Int = 2
  vtk_output::Bool = true

  H0 = 10 #m #still-water depth

  # Membrane parameters
  mᵨ = 0.005 #970.0*0.005/1025 #mass per unit area of membrane / ρw
  Tᵨ = 48.78 #1e7 stress x 0.005 thick #T/ρw
  τ = 0.0#damping coeff
  Pin_Wm = 12.0 # Projected Width of membrane for calculation of input power

  # Wave parameters
  ω = 2:1:3
  T = 2*π./ω
  η₀ = 0.10*ones(length(ω))
  # α = randomPhase(ω; seed=100)

  # Domain 
  Ld = 20 #damping zone length
  LΩ = 80
  WΩ = 40
  x₀ = -40
  xdᵢₙ = -20
  xdₒₜ = 20

end


end