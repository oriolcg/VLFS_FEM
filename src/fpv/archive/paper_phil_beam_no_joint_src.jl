module beam

using DrWatson
@quickactivate "MonolithicFEMVLFS"

using Gridap
using JLD2
using TickTock
using Plots
using DataFrames
using Roots
using WaveSpec
using Plots
using .Constants
using TimerOutputs

include(srcdir("fpv","materials.jl"))

export run_beam

  function run_beam(params::Dict)
    @unpack length_beam, h_b, material, phase, omega, amplitude,mesh_size = params    
    @unpack numFloat = params
    @show h_b
    mat = material 
    @show α = phase
    @show mesh_size
    @show ω = omega 
    @show A_w = amplitude
    @unpack depth = params
    d = depth


    ## Auxiliary functions

    # Check if an element is inside the beam
    function is_beam(coords) #applies on the cells, not on the nodes
      n = length(coords)
      #println(n)
      x = (1/n)*sum(coords)
      #println(x)
      
      (xb₀ <= x[1] <= xb₁ ) * ( x[2] ≈ d )
    end

    function is_beam_node(coords) 
      n = length(coords)
      #println(n)
      x = (1/n)*sum(coords)
      #println(x)
      
      (xb₀ < x[1] < xb₁ ) * ( x[2] ≈ d )
    end

    # Check if an element is inside the inlet damping zone
    function is_inlet_damping(xs)        
      n = length(xs) 
      x = (1/n)*sum(xs)
      (x[1] <= Ld_in ) * ( x[2] ≈ d )
    end

    # Check if an element is inside the outlet damping zone
    function is_outlet_damping(xs)
      n = length(xs)
      x = (1/n)*sum(xs)
      ( (LΓ - Ld_out) <= x[1] ) * ( x[2] ≈ d )
    end

    # Check if an element is on the beam boundary
    function is_beam_boundary(xs)
      is_on_xb₀ = [x[1]≈xb₀ for x in xs] 
      is_on_xb₁ = [x[1]≈xb₁ for x in xs]
      element_on_xb₀ = minimum(is_on_xb₀)
      element_on_xb₁ = minimum(is_on_xb₁)
      element_on_xb₀ | element_on_xb₁ 
    end

    # Check if an element is a joint
    function is_a_joint(xs)
      is_on_xj = Array{Bool,1}(UndefInitializer(),length(xs))
      is_on_xj .= false
      for xi in xj
        is_on_xj = is_on_xj .| [x[1]≈xi for x in xs] # array of booleans of size the number of points in an element (for points, it will be an array of size 1)
      end
      element_on_xj = minimum(is_on_xj) # Boolean with "true" if at least one entry is true, "false" otherwise.
      element_on_xj 
    end

    #Material properties
    #γ = 1
    ρ_b, E = give_material(mat)      # bending stiffness parameter
    
    ρ_b = ρ_b             # mass density VLFS [kg m^-3]
    #h_b = 1                # height VLFS [m] 
    L = length_beam         # total length VLFS [m]
    

    I = 1/12*h_b^3        # second moment of inertia per unit meter width (I/b) 
    #E = material[1]         # Youngs modulus [N/m^2] for HDPE
    EI_b = E * I     # bending stiffness VLFS per unit meter width [Nm/m]  
    ξ = 0                   # rotational stiffness parameter
    k_r = ξ*EI_b/L          # rotational stiffness connection [Nm]

    ρ_w = 1025                  # water density [kg m^-3]    
    println("Water Depth = ", d)  
    
    
    # Wave
    k = dispersionRelAng.(d, ω; msg=false)
    λ = 2*π./k
    
    # Parameters for simulation
    α1_b = ρ_b*h_b/ρ_w        
    α2_b = EI_b/ρ_w
    massPerArea = α1_b * ρ_w


    # Domain
    Ld_in  = 550  
    L_fs =   50   
    Lb     = L    
    Ld_out = 0 
    LΓ = Ld_in + L_fs + Lb + L_fs + Ld_out 
    @show domain = (0,LΓ,0,d)  

    xᵢₙ = Ld_in       
    xb₀ = Ld_in + L_fs
    xb₁ = xb₀ + Lb
    xd = xb₁ + L_fs
    
    # Floaters
    number_floaters = numFloat #1
    number_joints = number_floaters - 1
    β =1/number_floaters
    xj = Vector{Float64}()
    for i in 1:number_joints
      append!(xj, (xb₀ + β*i*Lb) )
    end  
    @show numFloat
    @show xj
        
    # Resolution
    nx = round(Int, LΓ/mesh_size)
    # nz = 20
    # mesh_rz = 1.15
    nz = 5
    mesh_rz = 1.85
    @show partition = (nx,nz)

    # Mesh
    function f_z(yin, r, n, H0; dbgmsg = false)
      y = yin - H0
      
      # Mesh along depth as a GP
      # Depth is 0 to -H0    
      if(r ≈ 1.0)
        return y + H0
      else
        a0 = H0 * (r-1) / (r^n - 1)    
        if(dbgmsg)
            ln = 0:n
            ly = H0.-a0 / (r-1) * (r.^ln .- 1)         
            @show hcat( ly, [ 0; ly[1:end-1] - ly[2:end] ] )
        end
        
        if y ≈ 0
            return H0
        end
        j = abs(y) / H0 * n  
        return H0-a0 / (r-1) * (r^j - 1)
      end
    end        
    map(x)  = VectorValue(x[1], f_z(x[2], mesh_rz, nz, d; 
      dbgmsg=false) )
    model_Ω = simplexify(CartesianDiscreteModel(domain,partition, map=map))

    labels_Ω = get_face_labeling(model_Ω)

    add_tag_from_tags!(labels_Ω,"surface",[3,4,6])       # assign the label "surface" to the entity 3,4 and 6 (top corners and top side)
    add_tag_from_tags!(labels_Ω,"bottom",[1,2,5])        # assign the label "bottom" to the entity 1,2 and 5 (bottom corners and bottom side)
    add_tag_from_tags!(labels_Ω,"inlet",[7])             # assign the label "inlet" to the entity 7 (left side)
    add_tag_from_tags!(labels_Ω,"outlet",[8])            # assign the label "outlet" to the entity 8 (right side)
    add_tag_from_tags!(labels_Ω, "water", [9])           # assign the label "water" to the entity 9 (interior)

    # Triangulations
    Ω = Interior(model_Ω)
    Γ = Boundary(model_Ω,tags="surface")
    Γin = Boundary(model_Ω,tags="inlet")
    Γout = Boundary(model_Ω, tags = "outlet")

    # Beam Triangulations 
    xΓ = get_cell_coordinates(Γ)                           # returns all the coordinates of the top layer (y = 30)

    grid_dim_0_Γ = Skeleton(Γ) #Edges not included in the Skeleton
    xΓ_dim_0 = get_cell_coordinates(grid_dim_0_Γ)
    Λj_to_Γ_mask = lazy_map(is_a_joint,xΓ_dim_0)
    Λj = Triangulation(grid_dim_0_Γ,Λj_to_Γ_mask)


    Λb_to_Γb_mask = lazy_map(is_beam_node,xΓ_dim_0)
    Λb_to_Γb_mask_without_j = Λb_to_Γb_mask .&& .!(Λj_to_Γ_mask) #remove the joints
    Λb = Triangulation(grid_dim_0_Γ,Λb_to_Γb_mask_without_j)

    Γb_to_Γ_mask = lazy_map(is_beam,xΓ)
    Γb_to_Γ = findall(Γb_to_Γ_mask)                        # results in indices of the beam. findall returns the indices where entry=true
    Γη_to_Γ = findall(Γb_to_Γ_mask)                        #eta is the elevation of  the beam. Is actually equal to Γb_to_Γ
    Γκ_to_Γ = findall(!, Γb_to_Γ_mask)

    Γb  = Triangulation(Γ,Γb_to_Γ)
    Γη  = Triangulation(Γ,Γη_to_Γ)
    Γκ  = Triangulation(Γ,Γκ_to_Γ)

    # Measures
    order = 4
    orderBM = 2
    degree = 2*order
    dΩ = Measure(Ω,degree)
    dΓb = Measure(Γb,degree)
    dΓfs = Measure(Γκ,degree)
    dΓin = Measure(Γin,degree)
    dΛb = Measure(Λb,degree)
    dΛj = Measure(Λj,degree)
    dΓot = Measure(Γout, degree)

    # Normals
    nΛb = get_normal_vector(Λb)
    nΛj = get_normal_vector(Λj)

    # FE spaces
    reffe = ReferenceFE(lagrangian,Float64,order)
    V_Ω = TestFESpace(Ω, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
    V_Γη = TestFESpace(Γη, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
    V_Γk = TestFESpace(Γκ, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
    U_Ω = TrialFESpace(V_Ω)
    U_Γη = TrialFESpace(V_Γη)
    U_Γκ = TrialFESpace(V_Γk)
    X = MultiFieldFESpace([U_Ω,U_Γη, U_Γκ])
    Y = MultiFieldFESpace([V_Ω,V_Γη,V_Γk])
    

    # Probe positions         
    xem_cords = xb₀:1.0:xb₁
    prb_xy = Point.(xem_cords,d)
    prbj_xy = Point.(xj,d)
    prb_wave = Point.( [572.7, 573.7, 575, 725], d )


    # Dataframes
    lDa = zeros(Complex, 1, length(xem_cords))
    da_η = DataFrame(lDa, :auto)
    da_ηx = DataFrame(lDa,:auto)
    da_ηxx = DataFrame(lDa, :auto)
    da_ηxxx = DataFrame(lDa,:auto)
    # da_jump_ηxxx = DataFrame(zeros(Complex, 1, length(xj)),:auto)
    da_wavePrb = DataFrame(zeros(Complex,1,4),:auto)
    

    function run_spectrum(ω,A_w,k,α,λ)

      tick()

      println("--- Freq --- ", ω)

      # Wave input
      ηᵢₙ(x) = A_w*exp(im*k*x[1] + im*α)
      ϕᵢₙ(x) = -im*(A_w*ω/k)*(cosh(k*x[2]) / 
          sinh(k*d))*exp(im*k*x[1] + im*α)
      vxᵢₙ(x) = (A_w*ω)*(cosh(k*x[2]) / 
          sinh(k*d))*exp(im*k*x[1] + im*α)
      vzfsᵢₙ(x) = -im*ω*A_w*exp(im*k*x[1] + im*α) #???

      # Damping
      μ₀ = 2.50
      function μ₁(x)
          lLdin = min(Ld_in, 3*λ)            
          return μ₀ * (1.0 - sin(π/2 * x[1] / lLdin)) * (x[1] < lLdin)
      end
      μ₂(x) = μ₁(x) * k 
      ηd(x) = μ₂(x)*ηᵢₙ(x) * (x[1]<xᵢₙ)
      ∇ₙϕd(x) = μ₁(x)*vzfsᵢₙ(x) * (x[1]<xᵢₙ) #???
      
      # Stabalisation
      h = mesh_size
      γ_m = 1.0*order*(order+1)
      βh_fs = 0.5
      αh_fs = -im*ω/g*(1-βh_fs)/βh_fs

      ## Weak form
      ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
      
      # Philipp
      # a((ϕ,η,κ),(w,v,u)) =  ∫(  ∇(w)⋅∇(ϕ) )dΩ   +   
      #   ∫(  (1im*ω*w*κ)  - μ₁*κ*w - μ₂*ϕ*w/g )dΓfs   +
      #   ∫(  βh_fs*(u + αh_fs*w)*g*κ  +   βh_fs*(-1im*ω)*(u + αh_fs*w)*ϕ )dΓfs   +
      #   ∫(  (-ω^2*α1_b + g)*v*η +  Δ(v)*(α2_b*Δ(η)) +  (-1im*ω*v*ϕ)   +   (1im*ω*w*η)  )dΓb   +
      #   ∫(   - (jump(∇(v)⋅nΛb) * mean(α2_b*Δ(η))) - (mean(Δ(v)) * jump(α2_b*∇(η)⋅nΛb))  + 
      #     γ_m/h*( jump(∇(v)⋅nΛb) * jump(α2_b*∇(η)⋅nΛb))  )dΛb +
      #   ∫(  (1/ρ_w)*(jump(∇(v)⋅nΛj) * k_r * jump(∇(η)⋅nΛj)) )dΛj -
      #   ∫(w * im * k * ϕ )dΓot      
    
      # l((w,v,u)) =            ∫( w*u_in )dΓin - ∫( w*ηd + w*ϕd/g )dΓfs


      # Weak form Beam Corrected
      a((ϕ,η,κ),(w,v,u)) =  
          ∫(  ∇(w)⋅∇(ϕ) )dΩ   +   
          ∫(  βh_fs*(u + αh_fs*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ - 
              μ₂*κ*w + μ₁*∇ₙ(ϕ)*(u + αh_fs*w) )dΓfs +
          ∫(  ( v*((-ω^2*α1_b + g)*η - im*ω*ϕ) + α2_b*Δ(v)*Δ(η) ) +  im*ω*w*η  )dΓb  +
          ∫(  α2_b * ( - jump(∇(v)⋅nΛb) * mean(Δ(η)) - mean(Δ(v)) 
              * jump(∇(η)⋅nΛb) + γ_m/h*( jump(∇(v)⋅nΛb) 
              * jump(∇(η)⋅nΛb) ) ) )dΛb +            
          ∫(  (1/ρ_w)*(jump(∇(v)⋅nΛj) * k_r * jump(∇(η)⋅nΛj)) )dΛj +
          ∫( -w * im * k * ϕ )dΓot        
    
      l((w,v,u)) = ∫( w*vxᵢₙ )dΓin - ∫( ηd*w - ∇ₙϕd*(u + αh_fs*w) )dΓfs


      # # Weak form EmptyTank Corrected
      # a((ϕ,η,κ),(w,v,u)) =  
      #     ∫(  ∇(w)⋅∇(ϕ) )dΩ   +   
      #     ∫(  βh_fs*(u + αh_fs*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ - 
      #         μ₂*κ*w + μ₁*∇ₙ(ϕ)*(u + αh_fs*w) )dΓfs +
      #     ∫(  βh_fs*(v + αh_fs*w)*(g*η - im*ω*ϕ) + im*ω*w*η  )dΓb  +            
      #     ∫( -w * im * k * ϕ )dΓot        
    
      # l((w,v,u)) = ∫( w*vxᵢₙ )dΓin - ∫( ηd*w - ∇ₙϕd*(u + αh_fs*w) )dΓfs
    
    
      ## Sovler
      op = AffineFEOperator(a, l, X, Y)
      ϕh , ηh, κh = Gridap.solve(op)
      

      ## Derivatives
      reffeBM = ReferenceFE(lagrangian,Float64,orderBM)
      V_Γη2 = FESpace(Γη, reffeBM, conformity=:L2, 
          vector_type=Vector{ComplexF64})                
      
      ηhx = ∇(ηh) ⋅ VectorValue(1.0, 0.0)                    
      
      # fTemp = interpolate_everywhere(ηhx, V_Γη2)
      ηhxx = ∇∇(ηh)⊙TensorValue(1.0,0.0,0.0,0.0)

      fTemp = interpolate_everywhere(ηhxx, V_Γη2)
      ηhxxx = ∇(fTemp) ⋅ VectorValue(1.0, 0.0) 

      # jump_ηhxxx = ∫(jump(∇(fTemp) ⋅ nΛj))dΛj
      # @show sum(jump_ηhxxx)       
      
      cache_η = Gridap.Arrays.return_cache(ηh, prb_xy)
      cache_ηx = Gridap.Arrays.return_cache(ηhx, prb_xy)
      cache_ηxx = Gridap.Arrays.return_cache(ηhxx, prb_xy)
      cache_ηxxx = Gridap.Arrays.return_cache(ηhxxx, prb_xy)      

      push!(da_η, evaluate!(cache_η, ηh, prb_xy))    
      push!(da_ηx, evaluate!(cache_ηx, ηhx, prb_xy))
      push!(da_ηxx, evaluate!(cache_ηxx, ηhxx, prb_xy))
      push!(da_ηxxx, evaluate!(cache_ηxxx, ηhxxx, prb_xy))
      # push!(da_jump_ηxxx, jump_ηhxxx.(prbj_xy))
      push!(da_wavePrb, κh.(prb_wave) )

      # println("Max Jump ", maximum(abs.(jump_ηhxxx.(prbj_xy))) )

      # writevtk(Γb, "temp", 
      #   cellfields=["eta"=> real(ηh), "eta_x"=> real(ηhx),
      #     "eta_xx"=> real(ηhxx), "eta_xxx"=> real(ηhxxx)])
      
      tock()
      println()
      # return  da_η, da_ηx, da_ηxx, da_ηxxx, da_jump_ηxxx
      return 0
    end
    
    
    run_spectrum.(omega,amplitude,k,phase,λ)

    
    da_η = da_η[2:end,:]
    da_ηx = da_ηx[2:end,:]
    da_ηxx = da_ηxx[2:end,:]
    da_ηxxx = da_ηxxx[2:end,:]
    # da_jump_ηxxx = da_jump_ηxxx[2:end,:]
    da_wavePrb = da_wavePrb[2:end,:]
    
    
    η_a = Matrix(abs.(da_η))
    η_ϕ = Matrix(angle.(da_η))
    
    ηx_a = Matrix(abs.(da_ηx))
    ηx_ϕ = Matrix(angle.(da_ηx))

    ηxx_a = Matrix(abs.(da_ηxx))

    ηxxx_a = Matrix(abs.(da_ηxxx))

    # ηxxx_jump_a = Matrix(abs.(da_jump_ηxxx))

    
    # ## Empty Tank    
    # function RAO(ηin_a, λin, η_a, ηx_a) 
    
    #   #expand the vectors to the size of length(ω),length(probx)
    #   ηin_a_m = repeat(ηin_a,1,size(η_a)[2])
    #   λin_m = repeat(λin,1,size(η_a)[2])      
      
    #   RAO_η = η_a ./ ηin_a_m
    #   RAO_ηx = ηx_a ./ ((2*π) * ηin_a_m ./ λin_m)
      
    #   return RAO_η, RAO_ηx
    # end
    
    # RAO_η, RAO_ηx = RAO(A_w, λ, η_a, ηx_a)

    # return RAO_η, RAO_ηx, η_ϕ, ηx_ϕ 


    # Not Empty tank    
    function RAO_revised(λin, η_a, ηx_a, ηxx_a, ηxxx_a) 
          
      #expand the vectors to the size of length(ω),length(probx)
      daRef = wload(datadir("fpv_202401","empt",
        "Empt_length_beam=100_mesh_size=1.0.jld2"))
      ηin_a_m = daRef["RAO_η"]
      λin_m = repeat(λin,1,size(η_a)[2])     
      k_m = 2*π ./ λin_m       
      
      RAO_η = η_a ./ ηin_a_m
      RAO_ηx = ηx_a ./ (k_m .* ηin_a_m)

      RAO_ηxx = ηxx_a ./ (k_m .* k_m .* ηin_a_m)
      RAO_ηxxx = ηxxx_a ./ (k_m .* k_m .* k_m .* ηin_a_m)       

      return RAO_η, RAO_ηx, RAO_ηxx, RAO_ηxxx
    end
    
    RAO_η, RAO_ηx, RAO_ηxx, RAO_ηxxx = RAO_revised(λ, η_a, 
      ηx_a, ηxx_a, ηxxx_a)

    return RAO_η, RAO_ηx, η_ϕ, ηx_ϕ, 
      RAO_ηxx, RAO_ηxxx, da_wavePrb,
      EI_b, massPerArea
    
  end

end

