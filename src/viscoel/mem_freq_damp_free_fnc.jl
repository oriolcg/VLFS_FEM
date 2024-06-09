module Memb2D

using Revise
using Gridap
using Plots
using DrWatson
using WaveSpec
using WaveSpec.Constants
using WaveSpec.Jonswap
using DataFrames:DataFrame
using DataFrames:Matrix
using TickTock
using Parameters


function main(params)

  function run_freq(ω, η₀, α)

    tick()
    # Wave parameters
    k = dispersionRelAng(H0, ω; msg=false)
    λ = 2π/k   
    T = 2π/ω
    ηᵢₙ(x) = η₀*exp(im*k*x[1] + im*α)
    ϕᵢₙ(x) = -im*(η₀*ω/k)*(cosh(k*(H0 + x[2])) / 
      sinh(k*H0))*exp(im*k*x[1] + im*α)
    vxᵢₙ(x) = (η₀*ω)*(cosh(k*(H0 + x[2])) / 
      sinh(k*H0))*exp(im*k*x[1] + im*α)
    vzfsᵢₙ(x) = -im*ω*η₀*exp(im*k*x[1] + im*α) #???
    @show ω, T
    @show λ
    @show η₀
    @show H0, H0/λ

    # Numeric constants
    αₕ = -im*ω/g * (1-βₕ)/βₕ
    @show αₕ
    println()  

    # Damping
    # Ldw = min( 5.0*λ, Ld )
    μ₀ = 2.5#maximum([2.5, 5.24/(ω^0.922)])#2.5
    μ₁ᵢₙ(x) =  μ₀*(1.0 - sin(π/2 * (x[1]-x₀)/Ld ))
    # μ₁ᵢₙ(x) =  μ₀*(1.0 - sin(π/2 * min( (x[1]-x₀)/Ldw, 1.0 ) ))
    μ₂ᵢₙ(x) = μ₁ᵢₙ(x)*k
    # μ₂ₒᵤₜ(x) = μ₁ₒᵤₜ(x)*k
    ηd(x) = μ₂ᵢₙ(x)*ηᵢₙ(x)
    ∇ₙϕd(x) = μ₁ᵢₙ(x)*vzfsᵢₙ(x) #???

    # Weak form
    ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
    a((ϕ,κ,η),(w,u,v)) =      
      ∫(  ∇(w)⋅∇(ϕ) )dΩ   +
      ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ )dΓfs   +
      ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ 
        - μ₂ᵢₙ*κ*w + μ₁ᵢₙ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd1    +
      ∫( -w * im * k * ϕ )dΓot +
      # ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ 
      #   - μ₂ₒᵤₜ*κ*w + μ₁ₒᵤₜ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd2    +
      ∫(  v*(g*η - im*ω*ϕ) +  im*ω*w*η
        - mᵨ*v*ω^2*η + Tᵨ*(1-im*ω*τ)*∇(v)⋅∇(η) )dΓm  #+ 
      #∫(- Tᵨ*(1-im*ω*τ)*v*∇(η)⋅nΛmb )dΛmb

    l((w,u,v)) =  ∫( w*vxᵢₙ )dΓin - ∫( ηd*w - ∇ₙϕd*(u + αₕ*w) )dΓd1


    # Solution
    op = AffineFEOperator(a,l,X,Y)
    (ϕₕ,κₕ,ηₕ) = solve(op)

    # Function for inlet phase
    κin = interpolate_everywhere(ηᵢₙ, 
      FESpace(Γκ, reffe, conformity=:H1, vector_type=Vector{ComplexF64}))
    κr = κₕ - κin

    # Energy flux (Power) calculation
    ηx = ∇(ηₕ)⋅VectorValue(1.0,0.0)
    Pd = sum(∫( abs(ηx)*abs(ηx) )dΓm)
    Pd = 0.5*Tᵨ*ρw*τ*ω*ω*Pd

    # Wave energy flux
    ηrf = abs(κr(Point(prbPowx[1],0.0)))
    ηtr = abs(κₕ(Point(prbPowx[2],0.0)))
    kh = k*H0
    wave_n = 0.5*(1 + 2*kh/sinh(2*kh))
    Pin = (0.5*ρw*g*η₀*η₀)*(ω/k)*wave_n
    Prf = (0.5*ρw*g*ηrf*ηrf)*(ω/k)*wave_n
    Ptr = (0.5*ρw*g*ηtr*ηtr)*(ω/k)*wave_n
    PErr = Pin - Prf - Ptr - Pd
    println("Power In \t ",Pin,"  W/m")
    println("Power Ref \t ",Prf," W/m")
    println("Power Trans \t ",Ptr," W/m")
    println("Power Abs \t ",Pd," W/m")
    println("Error \t ",PErr," W/m")

    # Interpolation on prboes
    prb_κ = zeros(ComplexF64, 1, length(prbxy))
    prb_κ_x = zeros(ComplexF64, 1, length(prbxy))  
    

    prb_κ[prbfs] = κₕ(prbxy[prbfs])
    prb_κ[prbmem] = ηₕ(prbxy[prbmem])

    prb_κ_x[prbfs] = (∇(κₕ)⋅VectorValue(1.0,0.0))(prbxy[prbfs])
    prb_κ_x[prbmem] = (∇(ηₕ)⋅VectorValue(1.0,0.0))(prbxy[prbmem])
  
    push!(prbDa, prb_κ)  
    push!(prbDa_x, prb_κ_x)  

    push!(prbDaΓη, ηₕ(prxΓη))
    push!(prbDaΓκ, κₕ(prxΓκ))

    push!(prbPow, [Pin, Prf, Ptr, Pd, PErr, 0.0])
    
    tock()
    return 0
  end


  @unpack name, order, vtk_output = params
  @show name
  @show order
  @show vtk_output
  filename = name*"/mem"

  @unpack H0, ω, T, η₀, α = params 
  k = dispersionRelAng.(H0, ω; msg=false)
  
  ρw = 1025 #kg/m3 water
  @show H0  #m #still-water depth
  @show ω

  # Peak Wave
  ωₚ, indp = findmax(η₀)
  @show ωₚ = ω[indp]
  kₚ = dispersionRelAng(H0, ωₚ; msg=false)
  println("Peak Wave T, L ", 2*pi/ωₚ, " ", 2*pi/kₚ)


  # Membrane parameters
  @unpack Lm, mᵨ, Tᵨ, τ = params
  @show Lm  #m
  @show g #defined in .Constants
  @show mᵨ #mass per unit area of membrane / ρw
  @show Tᵨ #T/ρw
  @show τ #damping coeff


  # Domain 
  @unpack nx, ny, mesh_ry, Ld, LΩ, x₀ = params
  @unpack domain, partition, xdᵢₙ, xm₀, xm₁ = params
  @show Lm
  @show LΩ, Ld
  @show domain
  @show partition
  @show mesh_ry
  @show (xm₀, xm₁)
  @show isinteger(Lm/LΩ*nx)
  @show LΩ/nx
  @show H0/ny
  println()


  # Numeric constants
  h = LΩ / nx
  γ = 1.0*order*(order-1)/h
  βₕ = 0.5
  # αₕ = -im*ω/g * (1-βₕ)/βₕ
  @show h
  @show βₕ
  # @show αₕ
  println()



  # Mesh
  function f_y(y, r, n, H0; dbgmsg = false)
    # Mesh along depth as a GP
    # Depth is 0 to -H0    
    if(r ≈ 1.0)
      return y  
    else
      a0 = H0 * (r-1) / (r^n - 1)    
      if(dbgmsg)
        ln = 0:n
        ly = -a0 / (r-1) * (r.^ln .- 1)         
        @show hcat( ly, [ 0; ly[1:end-1] - ly[2:end] ] )
      end
      
      if y ≈ 0
        return 0.0
      end
      j = abs(y) / H0 * n  
      return -a0 / (r-1) * (r^j - 1)
    end
  end
  map(x) = VectorValue( x[1], f_y(x[2], mesh_ry, ny, H0; dbgmsg=false) )
  model = CartesianDiscreteModel(domain,partition,map=map)


  # Labelling
  labels_Ω = get_face_labeling(model)
  add_tag_from_tags!(labels_Ω,"surface",[3,4,6])   # assign the label "surface" to the entity 3,4 and 6 (top corners and top side)
  add_tag_from_tags!(labels_Ω,"bottom",[1,2,5])    # assign the label "bottom" to the entity 1,2 and 5 (bottom corners and bottom side)
  add_tag_from_tags!(labels_Ω,"inlet",[7])         # assign the label "inlet" to the entity 7 (left side)
  add_tag_from_tags!(labels_Ω,"outlet",[8])        # assign the label "outlet" to the entity 8 (right side)
  add_tag_from_tags!(labels_Ω, "water", [9])       # assign the label "water" to the entity 9 (interior)


  # Triangulations
  Ω = Interior(model) #same as Triangulation()
  Γ = Boundary(model,tags="surface") #same as BoundaryTriangulation()
  Γin = Boundary(model,tags="inlet")
  Γot = Boundary(model,tags="outlet")


  # Auxiliar functions
  function is_mem(xs) # Check if an element is inside the beam1
    n = length(xs)
    x = (1/n)*sum(xs)
    (xm₀ <= x[1] <= xm₁ ) * ( x[2] ≈ 0.0)
  end
  function is_damping1(xs) # Check if an element is inside the damping zone 1
    n = length(xs)
    x = (1/n)*sum(xs)
    (x₀ <= x[1] <= xdᵢₙ ) * ( x[2] ≈ 0.0)
  end
  # function is_damping2(xs) # Check if an element is inside the damping zone 2
  #   n = length(xs)
  #   x = (1/n)*sum(xs)
  #   (xdₒₜ <= x[1] ) * ( x[2] ≈ 0.0)
  # end

  # Masking and Beam Triangulation
  xΓ = get_cell_coordinates(Γ)
  Γm_to_Γ_mask = lazy_map(is_mem, xΓ)
  Γd1_to_Γ_mask = lazy_map(is_damping1, xΓ)
  #Γd2_to_Γ_mask = lazy_map(is_damping2, xΓ)
  Γm = Triangulation(Γ, findall(Γm_to_Γ_mask))
  Γd1 = Triangulation(Γ, findall(Γd1_to_Γ_mask))
  #Γd2 = Triangulation(Γ, findall(Γd2_to_Γ_mask))
  Γfs = Triangulation(Γ, findall(!, Γm_to_Γ_mask .| 
    Γd1_to_Γ_mask ))# .| Γd2_to_Γ_mask))
  Γη = Triangulation(Γ, findall(Γm_to_Γ_mask))
  Γκ = Triangulation(Γ, findall(!,Γm_to_Γ_mask))


  # Construct the tag for membrane boundary
  Λmb = Boundary(Γm)
  xΛmb = get_cell_coordinates(Λmb)
  xΛmb_n1 = findall(model.grid_topology.vertex_coordinates .== xΛmb[1])
  xΛmb_n2 = findall(model.grid_topology.vertex_coordinates .== xΛmb[2])
  new_entity = num_entities(labels_Ω) + 1
  labels_Ω.d_to_dface_to_entity[1][xΛmb_n1[1]] = new_entity
  labels_Ω.d_to_dface_to_entity[1][xΛmb_n2[1]] = new_entity
  add_tag!(labels_Ω, "mem_bnd", [new_entity])


  writevtk(model, filename*"_model")
  if vtk_output == true
    writevtk(Ω,filename*"_O")
    writevtk(Γ,filename*"_G")
    writevtk(Γm,filename*"_Gm")  
    writevtk(Γd1,filename*"_Gd1")
    # writevtk(Γd2,filename*"_Gd2")
    writevtk(Γfs,filename*"_Gfs")
    writevtk(Λmb,filename*"_Lmb")  
  end


  # Measures
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓm = Measure(Γm,degree)
  dΓd1 = Measure(Γd1,degree)
  # dΓd2 = Measure(Γd2,degree)
  dΓfs = Measure(Γfs,degree)
  dΓin = Measure(Γin,degree)
  dΓot = Measure(Γot,degree)
  dΛmb = Measure(Λmb,degree)


  # Normals
  @show nΛmb = get_normal_vector(Λmb)


  # Dirichlet Fnc
  gη(x) = ComplexF64(0.0)

  # FE spaces
  reffe = ReferenceFE(lagrangian,Float64,order)
  V_Ω = TestFESpace(Ω, reffe, conformity=:H1, 
    vector_type=Vector{ComplexF64})
  V_Γκ = TestFESpace(Γκ, reffe, conformity=:H1, 
    vector_type=Vector{ComplexF64})
  # V_Γη = TestFESpace(Γη, reffe, conformity=:H1, 
  #   vector_type=Vector{ComplexF64},
  #   dirichlet_tags=["mem_bnd"]) #diri
  V_Γη = TestFESpace(Γη, reffe, conformity=:H1, 
    vector_type=Vector{ComplexF64})
  U_Ω = TrialFESpace(V_Ω)
  U_Γκ = TrialFESpace(V_Γκ)
  # U_Γη = TrialFESpace(V_Γη, gη)
  U_Γη = TrialFESpace(V_Γη)
  X = MultiFieldFESpace([U_Ω,U_Γκ,U_Γη])
  Y = MultiFieldFESpace([V_Ω,V_Γκ,V_Γη])

  # Probes
  @unpack prbx, prbPowx = params
  # prbx=[  -20.0, 0.0, 20.0, 40.0, 50.0, 
  #         52.7, 53.7, 55, 60.0, 80.0, 
  #         85.0, 90.0, 95.0, 100.0, 120.0, 
  #         125.0, 140.0, 160.0, 180.0 ]
  # prbPowx=[ 55.0, 125.0 ]
  @show prbxy = Point.(prbx, 0.0)
  prbmem = (xm₀ .<= prbx .<= xm₁ )
  @show prbfs = findall(!,prbmem)
  @show prbmem = findall(prbmem)

  lDa = zeros(ComplexF64, 1, length(prbxy))
  prbDa = DataFrame(lDa, :auto) # for η
  prbDa_x = DataFrame(lDa, :auto) #for ηₓ


  # Storing soln at Γ
  # Difficulties in doing it purely using get_cell_dof_valules()
  # instead going to do it using evaluate
  # maybe a bit slow, but shouldnt matter too much i guess.
  xΓη = get_cell_coordinates(Γη)
  xΓκ = get_cell_coordinates(Γκ)
  prxΓη = [val[1] for val in xΓη]
  tmp = [val[2] for val in xΓη]
  push!(prxΓη,tmp[end])
  prxΓκ = [val[1] for val in xΓκ]
  push!(prxΓκ,prxΓη[1])
  sort!(prxΓκ)
  lDa = zeros(ComplexF64, 1, length(prxΓη))
  prbDaΓη = DataFrame(lDa, :auto)
  lDa = zeros(ComplexF64, 1, length(prxΓκ))
  prbDaΓκ = DataFrame(lDa, :auto)

  prbPow = DataFrame(zeros(Float64, 1, 6), :auto)

  # Run weak-form for each freq
  run_freq.(ω, η₀, α)

  @show prbDa = prbDa[2:end, :]
  prbDa_x = prbDa_x[2:end, :]
  prbDaΓη = prbDaΓη[2:end,:]
  prbDaΓκ = prbDaΓκ[2:end,:]
  prbPow = prbPow[2:end,:]

  k = dispersionRelAng.(H0, ω; msg=false)

  for lprb in 1:length(prbxy)
    plt1 = plot(k*H0, abs.(prbDa[:,lprb]), linewidth=3, 
      xlabel = "kh",
      ylabel = "A (m)",
      title = "Amplitude")  

    plt2 = plot(k*H0, abs.(prbDa_x[:,lprb]), linewidth=3, 
      xlabel = "kh",
      ylabel = "dA/dx",
      title = "Slope Magnitude")
    
    plt3 = plot(k*H0, angle.(prbDa[:,lprb]), linewidth=3, 
      xlabel = "kh",
      ylabel = "α (rad)",
      title = "Phase")  

    plt4 = plot(k*H0, angle.(prbDa_x[:,lprb]), linewidth=3, 
      xlabel = "kh",
      ylabel = "α (rad)",
      title = "Slope Phase")
    
    xloc = prbx[lprb]
    pltAll = plot(plt1, plt2, plt3, plt4, layout=4, dpi=330,
      plot_title = "x = $xloc")

    savefig(pltAll,filename*"_dxPrb_$lprb"*".png")
  end  

  data = Dict("ω" => ω,
              "η₀" => η₀,
              "k" => k,
              "prbxy" => prbxy,
              "prbDa" => prbDa,            
              "prbDa_x" => prbDa_x,
              "prxΓκ" => prxΓκ,
              "prxΓη" => prxΓη,
              "prbDaΓκ" => prbDaΓκ,
              "prbDaΓη" => prbDaΓη,
              "prbPow" => prbPow)

  wsave(filename*"_data.jld2", data)

end


"""
Memb_params

Parameters for the VIV.jl module.
"""
@with_kw struct Memb_params_warmup
  name::String = "data/sims_202306/run/spec_free"
  order::Int = 2
  vtk_output::Bool = true

  H0 = 10 #m #still-water depth

  # Wave parameters
  # ω, S, η₀ = jonswap(0.4, 2.5; 
  #     plotflag=true, plotloc=filename, nω=145)
  # println(ω[1], "\t", ω[2], "\t", ω[end])
  # ω = ω[2:end]
  # S = S[2:end]
  # η₀ = η₀[2:end]
  # ω = [2*π/2.53079486745378, 2*π/2.0]
  # η₀ = [0.25, 0.25]
  ω = 1:1:2
  T = 2*π./ω
  η₀ = 0.10*ones(length(ω))
  α = randomPhase(ω; seed=100)
  # k = dispersionRelAng.(H0, ω; msg=false)

  # Membrane parameters
  Lm = 2*H0 #m
  Wm = Lm  
  mᵨ = 0.9 #mass per unit area of membrane / ρw
  Tᵨ = 0.1/4*g*Lm*Lm #T/ρw
  τ = 0.0#damping coeff


  # Domain 
  nx = 330
  ny = 10
  mesh_ry = 1.2 #Ratio for Geometric progression of eleSize
  Ld = 15*H0 #damping zone length
  LΩ = 18*H0 + Ld #2*Ld
  x₀ = -Ld
  domain =  (x₀, x₀+LΩ, -H0, 0.0)
  partition = (nx, ny)
  xdᵢₙ = 0.0
  xm₀ = xdᵢₙ + 8*H0
  xm₁ = xm₀ + Lm

  # Probes
  prbx=[  -20.0, 0.0, 20.0, 40.0, 50.0, 
          52.7, 53.7, 55, 60.0, 80.0, 
          85.0, 90.0, 95.0, 100.0, 120.0, 
          125.0, 140.0, 160.0, 180.0 ]
  prbPowx=[ 55.0, 125.0 ]

end

@with_kw struct Memb_params
  name::String = "data/sims_202306/run/spec_free"
  order::Int = 2
  vtk_output::Bool = true

  H0 = 10 #m #still-water depth

  # Wave parameters
  # ω, S, η₀ = jonswap(0.4, 2.5; 
  #     plotflag=true, plotloc=filename, nω=145)
  # println(ω[1], "\t", ω[2], "\t", ω[end])
  # ω = ω[2:end]
  # S = S[2:end]
  # η₀ = η₀[2:end]
  # ω = [2*π/2.53079486745378, 2*π/2.0]
  # η₀ = [0.25, 0.25]
  ω = 0.7:0.5:5
  T = 2*π./ω
  η₀ = 0.10*ones(length(ω))
  α = randomPhase(ω; seed=100)
  # k = dispersionRelAng.(H0, ω; msg=false)

  # Membrane parameters
  Lm = 2*H0 #m
  Wm = Lm  
  mᵨ = 0.9 #mass per unit area of membrane / ρw
  Tᵨ = 0.1/4*g*Lm*Lm #T/ρw
  τ = 0.0#damping coeff


  # Domain 
  nx = 1650
  ny = 20
  mesh_ry = 1.2 #Ratio for Geometric progression of eleSize
  Ld = 15*H0 #damping zone length
  LΩ = 18*H0 + Ld #2*Ld
  x₀ = -Ld
  domain =  (x₀, x₀+LΩ, -H0, 0.0)
  partition = (nx, ny)
  xdᵢₙ = 0.0
  xm₀ = xdᵢₙ + 8*H0
  xm₁ = xm₀ + Lm

  # Probes
  prbx=[  -20.0, 0.0, 20.0, 40.0, 50.0, 
          52.7, 53.7, 55, 60.0, 80.0, 
          85.0, 90.0, 95.0, 100.0, 120.0, 
          125.0, 140.0, 160.0, 180.0 ]
  prbPowx=[ 55.0, 125.0 ]

end

end