module Memb3D_modes

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
using Parameters


function run_case(params, mfac = 0.9, tfac = 0.1)
  
  function run_freq(ω)

    k = dispersionRelAng(H0, ω; msg=false)
    @show ω, k
  
    # Weak form: ω dependent
    k22(ϕ,w) = ∫( ∇(w)⋅∇(ϕ) )dΩ +
        ∫( -w * im * k * ϕ )dΓin + ∫( -w * im * k * ϕ )dΓot
  
    c23(κ,w) = ∫( im*ω*w*κ )dΓfs 
  
    c32(ϕ,u) = ∫( -im*ω*u*ϕ )dΓfs 
    
    # Global matrices: ω dependent
    K22 = get_matrix(AffineFEOperator( k22, l2, U_Ω, V_Ω ))
    C23 = get_matrix(AffineFEOperator( c23, l2, U_Γκ, V_Ω ))
    C32 = get_matrix(AffineFEOperator( c32, l3, U_Ω, V_Γκ ))
  
    # Solution
    tick()
    Mϕ = K22 - ( C23 * (Matrix(K33) \ C32) )
    Mhat = C12 * (Mϕ \ C21)
    Mtot = M11 + Mhat
    tock()
  
    # Eigen values
    # λ = LinearAlgebra.eigvals(Mtot\Matrix(K11))
    # V = LinearAlgebra.eigvecs(Mtot\Matrix(K11))  
    # # @show real.(λ[1:nωₙ])
    # # ωₙ = sqrt.(real.(λ))
    Ur, S, Vr = svd(Mtot\Matrix(K11))
    λ = reverse(S)
    V = reverse(Vr, dims=2)
    rλ = real.(λ[1:nωₙ])
    @show rλ
    return(rλ[1:nωₙ], V[:,1:nωₙ])
      
  end

  @unpack resDir, order, order2, vtk_output = params
  caseName = "ten" * @sprintf("%0.2f", tfac) *"_mass" * @sprintf("%0.2f", mfac)
  name::String = resDir*"/mem_modes_"*caseName
  filename = name*"/mem"

  # if( isdir(name) )
  #   return
  # end
  # mkdir(name)

  if( ~isdir(name) )
    mkdir(name)
  end  

  @unpack H0, Lm, Wm,  = params
  ρw = 1025 #kg/m3 water
  @show H0  #m #still-water depth

  # Membrane parameters
  @show Lm #m
  @show Wm #m 
  @show g #defined in .Constants
  mᵨ = mfac  #mass per unit area of membrane / ρw
  Tᵨ = tfac/4*g*Lm*Lm #T/ρw

  # Excitation wave parameters
  ω = 1.0

  # Domain
  @unpack nx, ny, nz, mesh_rz, LΩ, WΩ = params 
  @unpack x₀, domain, partition, xm₀, xm₁, ym₀, ym₁ = params
  @show Lm, Wm
  @show LΩ, WΩ
  @show domain
  @show partition
  @show (xm₀, xm₁) (ym₀, ym₁)
  @show isinteger(Lm/LΩ*nx)
  @show LΩ/nx
  @show H0/nz
  println()
  #@show Ld*k/2/π
  #@show cosh.(k*H0*0.5)./cosh.(k*H0)


  # Mesh
  @unpack dbgmsg = params
  function f_z(z, r, n, H0; dbgmsg = false)
    # Mesh along depth as a GP
    # Depth is 0 to -H0    
    if(r ≈ 1.0)
      return z 
    else
      a0 = H0 * (r-1) / (r^n - 1)    
      if(dbgmsg)
        ln = 0:n
        lz = -a0 / (r-1) * (r.^ln .- 1)         
        @show hcat( lz, [ 0; lz[1:end-1] - lz[2:end] ] )
      end
      
      if z ≈ 0
        return 0.0
      end
      j = abs(z) / H0 * n  
      return -a0 / (r-1) * (r^j - 1)
    end
  end
  map(x) = VectorValue( x[1], x[2], f_z(x[3], 
    mesh_rz, nz, H0; dbgmsg=dbgmsg) )
  model = CartesianDiscreteModel(domain,partition,map=map)


  # Labelling
  labels_Ω = get_face_labeling(model)
  add_tag_from_tags!(labels_Ω,"surface",
    [5,6,7,8,11,12,15,16,22])   
  add_tag_from_tags!(labels_Ω,"bottom",
    [1,2,3,4,9,10,13,14,21])
  add_tag_from_tags!(labels_Ω,"inlet",
    [17,19,25]) 
  add_tag_from_tags!(labels_Ω,"outlet",
    [18,20,26])   
  add_tag_from_tags!(labels_Ω,"swNear",
    [23])   
  add_tag_from_tags!(labels_Ω,"swFar",
    [24])   


  # Triangulations
  Ω = Interior(model) #same as Triangulation()
  Γ = Boundary(model,tags="surface") #same as BoundaryTriangulation()
  Γin = Boundary(model,tags="inlet")
  Γot = Boundary(model,tags="outlet")


  # Auxiliar functions
  function is_mem(xs) # Check if an element is inside the beam1
    n = length(xs)
    x = (1/n)*sum(xs)
    (xm₀ <= x[1] <= xm₁ ) * (ym₀ <= x[2] <= ym₁) * ( x[3] ≈ 0.0)
  end


  # Masking and Beam Triangulation
  xΓ = get_cell_coordinates(Γ)
  Γm_to_Γ_mask = lazy_map(is_mem, xΓ)
  Γm = Triangulation(Γ, findall(Γm_to_Γ_mask))
  Γfs = Triangulation(Γ, findall(!, Γm_to_Γ_mask))
  Γη = Triangulation(Γ, findall(Γm_to_Γ_mask))
  Γκ = Triangulation(Γ, findall(!,Γm_to_Γ_mask))


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
    writevtk(Γfs,filename*"_Gfs")
    # writevtk(Λmb,filename*"_Lmb")  
  end


  # Measures
  degree = 2*order
  degree2 = 2*order2
  dΩ = Measure(Ω,degree)
  dΓm = Measure(Γm,degree2)
  dΓfs = Measure(Γfs,degree2)
  dΓin = Measure(Γin,degree)
  dΓot = Measure(Γot,degree)
  # dΛmb = Measure(Λmb,degree)


  # Normals
  # @show nΛmb = get_normal_vector(Λmb)


  # Dirichlet Fnc
  gη(x) = ComplexF64(0.0)

  # FE spaces
  reffe = ReferenceFE(lagrangian,Float64,order)
  reffe2 = ReferenceFE(lagrangian,Float64,order2)
  V_Ω = TestFESpace(Ω, reffe, conformity=:H1, 
    vector_type=Vector{ComplexF64})
  V_Γκ = TestFESpace(Γκ, reffe2, conformity=:H1, 
    vector_type=Vector{ComplexF64})
  # V_Γη = TestFESpace(Γη, reffe, conformity=:H1, 
  #   vector_type=Vector{ComplexF64},
  #   dirichlet_tags=["mem_bnd"]) #diri
  V_Γη = TestFESpace(Γη, reffe2, conformity=:H1, 
    vector_type=Vector{ComplexF64})
  U_Ω = TrialFESpace(V_Ω)
  U_Γκ = TrialFESpace(V_Γκ)
  # U_Γη = TrialFESpace(V_Γη, gη) #diri
  U_Γη = TrialFESpace(V_Γη)


  # Weak form: Constant matrices
  ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
  m11(η,v) = ∫( mᵨ*v*η )dΓm
  k11(η,v) = ∫( v*g*η + Tᵨ*∇(v)⋅∇(η) )dΓm #+  
              # ∫(- Tᵨ*v*∇(η)⋅nΛmb )dΛmb #diri

  c12(ϕ,v) = ∫( v*ϕ )dΓm

  c21(η,w) = ∫( w*η )dΓm  

  k33(κ,u) = ∫( u*g*κ )dΓfs

  l1(v) = ∫( 0*v )dΓm
  l2(w) = ∫( 0*w )dΩ
  l3(u) = ∫( 0*u )dΓfs 
  println("[MSG] Done Weak form")

  # Global matrices: constant matrices
  M11 = get_matrix(AffineFEOperator( m11, l1, U_Γη, V_Γη ))
  K11 = get_matrix(AffineFEOperator( k11, l1, U_Γη, V_Γη ))
  C12 = get_matrix(AffineFEOperator( c12, l1, U_Ω, V_Γη ))

  C21 = get_matrix(AffineFEOperator( c21, l2, U_Γη, V_Ω ))

  K33 = get_matrix(AffineFEOperator( k33, l3, U_Γκ, V_Γκ ))
  println("[MSG] Done Global matrices")
  println(K11 == transpose(K11))

  #xp = range(xm₀, xm₁, size(V,2)+2)

  @unpack nωₙ, errLim, maxIter = params
  da_ωₙ = zeros(Float64, 1, nωₙ)
  @show ωₙ=zeros(Float64, 1, nωₙ) .+ ω
  da_V = []

  # # For index=1 not looping coz ωₙ[1] = 0.0
  # i = 1
  # ωₙ, V = run_freq(ω)
  # da_ωₙ[i] = ωₙ[1]
  # push!(da_V, V[:,i])

  for i in 1:nωₙ
    # global da_ωₙ, da_V  
    # global ωₙ, ω
    local V, lIter
    lIter = 0    
    Δω = 1
    ω = ωₙ[i]    
    while ((Δω > errLim) && (lIter < maxIter))
      # global ω, ωₙ
      rλ, V = run_freq(ω)
      ωₒ = ω      
      ωᵣ = sqrt(rλ[i])
      # if(i==1)
      #   #ω = 0.2 * ωₙ[i] + 0.8*ω
      #   ω = 0.0
      #   Δω = 0.0
      #   V = V*0.0
      # elseif(i==4)
      #   ω = 0.2 * ωᵣ + 0.8*ωₒ
      #   Δω = abs(ω - ωₒ)/ωₒ
      # else
      #   ω = 0.5 * ωᵣ + 0.5*ωₒ
      #   Δω = abs(ω - ωₒ)/ωₒ
      # end      
      ω = 0.5 * ωᵣ + 0.5*ωₒ
      Δω = abs(ω - ωₒ)/ωₒ
      # @show ωₙ      
      lIter += 1
      @show i, ω, Δω, lIter
    end
    da_ωₙ[i] = ω
    push!(da_V, V[:,i])
  end

  println(da_ωₙ)

  xp = range(xm₀, xm₁, length(da_V[1]))

  data = Dict(
    "xp" => xp,
    "ωₙ" => da_ωₙ,
    "V" => da_V  
  )

  # data = Dict(
  #   "model" => model, 
  #   "Γκ" => Γκ, 
  #   "V_Γκ" => V_Γκ, 
  #   "reffe2" => reffe2, 
  #   "order2" => order2)

  wsave(filename*"_modesdata.jld2", data)
end


function main(params, mfac, tfac)

  # mfac = [0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0 ]
  # tfac = [0.05, 0.10, 0.25, 0.50, 0.75]

  # # mfac = [0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0 ]
  # # tfac = [0.8]

  # mfac = [0.9]
  # tfac = [0.1]

  for imfac in mfac
    for itfac in tfac
      run_case(params, imfac, itfac)
    end
  end
end


"""
Memb_params

Parameters for the VIV.jl module.
"""
@with_kw struct Memb_params
  resDir::String = "data/sims_202309/mem3D_modes_res_free"
  order::Int = 1
  order2::Int = 2
  vtk_output::Bool = true

  H0 = 10 #m #still-water depth

  # Membrane parameters
  Lm = 2*H0 #m
  Wm = Lm  

  nωₙ = 2
  errLim = 1e-1
  maxIter = 15

  dbgmsg = false

  # Domain 
  nx = 30
  ny = 8
  nz = 4
  mesh_rz = 1.2 #Ratio for Geometric progression of eleSize
  LΩ = 6*H0 
  WΩ = 4*H0
  x₀ = 0.0
  domain =  (x₀, x₀+LΩ, -WΩ/2, WΩ/2 , -H0, 0.0)
  partition = (nx, ny, nz)
  xm₀ = x₀ + 2*H0
  xm₁ = xm₀ + Lm
  ym₀ = -Wm/2
  ym₁ = Wm/2

end


end