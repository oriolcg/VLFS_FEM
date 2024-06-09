module beam

using DrWatson
@quickactivate "MonolithicFEMVLFS"

using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.FESpaces
using Gridap.Visualization
using JLD2
using TickTock
using Plots
using DataFrames
using Random
using Roots
using WaveSpec
using Plots
using .Constants
using .Jonswap
using .WaveTimeSeries
using TimerOutputs

#using ..mat
include(srcdir("fpv","materials.jl"))

#include(srcdir("materials.jl"))


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
        # xs is a vector of points with the coordinates of a given element of dimension "d"
        n = length(xs) # number of points in each element (in an edge there will be two points)
        x = (1/n)*sum(xs) # take the average (centre of the element)
        (x[1] <= Ld_in ) * ( x[2] ≈ d ) # check if the centre is in the inlet damping zone and on the surface
    end
    # Check if an element is inside the outlet damping zone
    function is_outlet_damping(xs)
        # xs is a vector of points with the coordinates of a given element of dimension "d"
        n = length(xs) # number of points in each element (in an edge there will be two points)
        x = (1/n)*sum(xs) # take the average (centre of the element)
        ( (LΓ - Ld_out) <= x[1] ) * ( x[2] ≈ d ) # check if the centre is in the inlet damping zone and on the surface
    end
    # Check if an element is on the beam boundary
    function is_beam_boundary(xs)
        is_on_xb₀ = [x[1]≈xb₀ for x in xs] # array of booleans of size the number of points in an element (for points, it will be an array of size 1)
        is_on_xb₁ = [x[1]≈xb₁ for x in xs]
        element_on_xb₀ = minimum(is_on_xb₀) # Boolean with "true" if at least one entry is true, "false" otherwise.
        element_on_xb₁ = minimum(is_on_xb₁)
        element_on_xb₀ | element_on_xb₁ # Return "true" if any of the two cases is true
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
    ρ_b, E, γ = give_material(mat)      # bending stiffness parameter
    
    ρ_b = ρ_b             # mass density VLFS [kg m^-3]
    #h_b = 1                # height VLFS [m] 
    L = length_beam         # total length VLFS [m]
    

    I = 1/12*h_b^3 * γ         # second moment of inertia per unit meter width (I/b) 
    #E = material[1]         # Youngs modulus [N/m^2] for HDPE
    EI_b = E * I     # bending stiffness VLFS per unit meter width [Nm/m]  
    ξ = 0                   # rotational stiffness parameter
    k_r = ξ*EI_b/L          # rotational stiffness connection [Nm]

    g = 9.81                    # gravitational constant
    ρ_w = 1025                  # water density [kg m^-3]
    d = 30                      # water depth [m]
    
    
    #Wave
    k = dispersionRelAng.(d, ω; msg=false)
    λ = 2*π./k


    #Set up for Gridap simulation
    # ==================
    
    #parameters for simulation
    α1_b = ρ_b*h_b/ρ_w        
    α2_b = EI_b/ρ_w


    #domain

    Ld_in  = 550  # = 1.5* wavelength omega = 0.1 #250                                              # Inlet damping zone length, depends on the 
    L_fs =   200                                               # Free water surface before and after the structure 
    Lb     = L                                                  # VLFS length
    Ld_out = 0 #250                                                 # Outlet damping zone length
    LΓ = Ld_in+L_fs+Lb+L_fs+Ld_out                                      # length fluid domain inc. damping zones
    @show domain = (0,LΓ,0,d)                                         # fluid domain. 
                                                                # domain is always built in (0.0 x 0.0 y 0.0 z) format

    xᵢₙ = Ld_in                                                  # x-coordinate end frontal damping zone
    xb₀ = Ld_in + 200 #0.5(L_fd-L)                                  # x-coordinate start of VLFS 
    xb₁ = xb₀ + Lb                                               # x-coordinate end of VLFS
    xd = xb₁ + 200 #0.5(L_fd-L)                                       # x-coordinate initial outlet damping zone point
    #parameters for the floaters
    number_floaters = numFloat #1
    number_joints = number_floaters - 1
    β =1/number_floaters
    xj = Vector{Float64}()

    for i in 1:number_joints
    append!(xj, (xb₀ + β*i*Lb) )
    end  
    @show numFloat
    @show xj
    
    #=damping zone
    μ₀ = 10

    μ₁(x::VectorValue) = μ₀* (1.0 - cos(π/2*(x[1]-xd)/Ld_out)) * (x[1] > xd) + μ₀*(1.0 - sin(π/2*x[1]/Ld_in)) * (x[1] < xᵢₙ) 
    μ₂(x::VectorValue) = -(μ₁(x)^2)/4

    ηd = x -> μ₁(x)*η_in(x)*(x[1]<xᵢₙ)
    ϕd = x -> μ₂(x)*ϕ_in(x)*(x[1]<xᵢₙ)
    =#

    #discrete model of the full domain
    
    nx = round(Int, LΓ/mesh_size)
    # nz = 20
    # mesh_rz = 1.15
    nz = 5
    mesh_rz = 1.85
    @show partition = (nx,nz)    #

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
    
    map(x)  = VectorValue(x[1], f_z(x[2], mesh_rz, nz, d; dbgmsg=false) )
    model_Ω = simplexify(CartesianDiscreteModel(domain,partition, map=map))            #CartesianDiscreteModel creates a mesh, 

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
    

    #Postprocessing
    #______Postprocessing______#

    y_cords = [d]
    x_cords = collect(xb₀:1:xb₁) 
    y_cords = repeat(y_cords,length(x_cords),1)
    prb_xy = Point.(x_cords,y_cords)[:,1]

    #prob_x = LinRange(xb₀,xb₁,Integer((xb₁-xb₀)*10)) #creates a probe every 10cm

    lDa = zeros(Complex, 1, length(x_cords))
    solution_η = DataFrame(lDa, :auto)
    solution_ηx = DataFrame(lDa,:auto)
    solution_ηxx = DataFrame(lDa, :auto)
    solution_ηxxx = DataFrame(lDa,:auto)

    #= filename = "F:/Github/Thesis-LFPV/Offshore PV/data/VTKOutput/model_"
    writevtk(Ω,filename*"_O") =#

    function run_spectrum(ω,A_w,k,α,λ)

        println("--- Freq --- ", ω)

        ηᵢₙ(x) = A_w*exp(im*k*x[1] + im*α)
        ϕᵢₙ(x) = -im*(A_w*ω/k)*(cosh(k*x[2]) / 
            sinh(k*d))*exp(im*k*x[1] + im*α)
        vxᵢₙ(x) = (A_w*ω)*(cosh(k*x[2]) / 
            sinh(k*d))*exp(im*k*x[1] + im*α)
        vzfsᵢₙ(x) = -im*ω*A_w*exp(im*k*x[1] + im*α) #???
          
        #damping zone
        #μ₀ = maximum([2.85, 10.26/(ω^0.760)]) #5.61/ω#6/(ω^0.6)        
        # function μ₁(x)
        #     if ω > 0.5
        #         if ( x[1] / (3λ) < 1.0 ) && (x[1] < xᵢₙ)
        #             return μ₀ * (1.0 - sin(π/2 * x[1] / (3λ)))
        #         else
        #             return 0.0
        #         end
        #     else
        #         if x[1] < xᵢₙ
        #             return μ₀ * (1.0 - sin(π/2 * x[1] / Ld_in))
        #         else
        #             return 0.0
        #         end
        #     end
        # end
        μ₀ = 2.50
        function μ₁(x)
            lLdin = min(Ld_in, 3*λ)            
            return μ₀ * (1.0 - sin(π/2 * x[1] / lLdin)) * (x[1] < lLdin)
        end
        # μ₁(x) = μ₀*(1.0 - sin(π/2 * x[1]/Ld_in ))
        μ₂(x) = μ₁(x) * k #Factor 1.9 based on trial and error

        
        ηd(x) = μ₂(x)*ηᵢₙ(x) * (x[1]<xᵢₙ)
        ∇ₙϕd(x) = μ₁(x)*vzfsᵢₙ(x) * (x[1]<xᵢₙ) #???
      
        h = mesh_size
        γ_m = 1.0*order*(order+1)
        βh_fs = 0.5
        αh_fs = -im*ω/g*(1-βh_fs)/βh_fs

        ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
      
      
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
      
      
      
        println("time to set up AffineFEOperator")
        tick()
        op = AffineFEOperator(a, l, X, Y)
        tock()

        println("Time so solve")
        tick()
        ϕh , ηh, κh = Gridap.solve(op)
        tock()
        
        # Derivatives
        reffeBM = ReferenceFE(lagrangian,Float64,orderBM)
        V_Γη2 = FESpace(Γη, reffeBM, conformity=:H1, 
            vector_type=Vector{ComplexF64})                
        
        ηhx = ∇(ηh) ⋅ VectorValue(1.0, 0.0)                    
        
        # fTemp = interpolate_everywhere(ηhx, V_Γη2)
        ηhxx = ∇∇(ηh)⊙TensorValue(1.0,0.0,0.0,0.0)

        fTemp = interpolate_everywhere(ηhxx, V_Γη2)
        ηhxxx = ∇(fTemp) ⋅ VectorValue(1.0, 0.0) 

        jump_ηhxxx = jump(∇(fTemp))

        # tt = ηhx
        # tt2 = ηhxx        
        # tt2b = ηhxx        
        # @show tt(Point(762.0-0.00001,30.0))
        # @show tt(Point(762.0+0.00001,30.0))
        # @show tt2(Point(762.0-0.00001,30.0))
        # @show tt2b(Point(762.0-0.00001,30.0))
        # @show tt2(Point(762.0+0.00001,30.0))
        # @show tt2b(Point(762.0+0.00001,30.0))
        # @show tt2(Point(762.0,30.0))
        # @show tt2b(Point(762.0,30.0))
        # @show jump_ηhxxx(Point(762.0,30.0))
        
        # @show ηhx(Point(761.5, 30.0))
        # @show ηhxx(Point(761.5, 30.0))        
        # @show ηhxxx(Point(761.5, 30.0))        
        
        println("Time to evaluate")
        tick()
        cache_η = Gridap.Arrays.return_cache(ηh,prb_xy)
        cache_ηx = Gridap.Arrays.return_cache(ηhx,prb_xy)
        cache_ηxx = Gridap.Arrays.return_cache(ηhxx,prb_xy)
        cache_ηxxx = Gridap.Arrays.return_cache(ηhxxx,prb_xy)

        push!(solution_η, evaluate!(cache_η, ηh, prb_xy))    
        push!(solution_ηx,evaluate!(cache_ηx, ηhx, prb_xy))
        push!(solution_ηxx,evaluate!(cache_ηxx, ηhxx, prb_xy))
        push!(solution_ηxxx,evaluate!(cache_ηxxx, ηhxxx, prb_xy))
        tock()
        

        return  solution_η, solution_ηx, solution_ηxx, solution_ηxxx
    end
    
    run_spectrum.(omega,amplitude,k,phase,λ)

    solution_η = solution_η[2:end,:]
    solution_ηx = solution_ηx[2:end,:]
    solution_ηxx = solution_ηxx[2:end,:]
    solution_ηxxx = solution_ηxxx[2:end,:]
    

    
    η_amplitude = Matrix(abs.(solution_η))
    η_ϕ = Matrix(angle.(solution_η))
    
    ηx_amplitude = Matrix(abs.(solution_ηx))
    ηx_ϕ = Matrix(angle.(solution_ηx))

    ηxx_amplitude = Matrix(abs.(solution_ηxx))
    ηxx_ϕ = Matrix(angle.(solution_ηxx))

    ηxxx_amplitude = Matrix(abs.(solution_ηxxx))
    ηxxx_ϕ = Matrix(angle.(solution_ηxxx))

    # #A_w tested in an empty water tank
    # A_w_real = [1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0001; 1.0; 1.0001; 1.0001; 1.0; 1.0001; 1.0001; 1.0; 1.0001; 1.0002; 1.0002; 1.0001; 1.0001; 1.0002; 1.0003; 1.0003; 1.0003; 1.0002; 1.0002; 1.0002; 1.0002; 1.0003; 1.0003; 1.0003; 1.0003; 1.0003; 1.0003; 1.0004; 1.0004; 1.0004; 1.0005; 1.0005; 1.0005; 1.0005; 1.0005; 1.0005; 1.0006; 1.0006; 1.0006; 1.0006; 1.0007; 1.0007; 1.0007; 1.0007; 1.0008; 1.0008; 1.0008; 1.0008; 1.0009; 1.0009; 1.0009; 1.001; 1.001; 1.001; 1.001; 1.0011; 1.0011; 1.0012; 1.0012; 1.0012; 1.0012; 1.0013; 1.0013; 1.0014; 1.0014; 1.0014; 1.0015; 1.0015; 1.0016; 1.0016; 1.0017; 1.0017; 1.0017; 1.0018; 1.0018; 1.0019; 1.0019; 1.002; 1.002; 1.0021; 1.0021; 1.0022; 1.0023; 1.0023; 1.0024; 1.0024; 1.0025; 1.0025; 1.0026; 1.0027; 1.0027; 1.0028; 1.0029; 1.0029; 1.003; 1.0031; 1.0031; 1.0032; 1.0033; 1.0033; 1.0034; 1.0035; 1.0036; 1.0037; 1.0037; 1.0038; 1.0039; 1.004; 1.0041; 1.0045; 1.005; 1.0056; 1.0062; 1.0069; 0.9988; 0.9981; 0.9972; 0.9961; 0.9947; 0.993; 0.9911; 0.9889; 0.9864; 0.9836; 0.9805; 0.9771; 0.9734; 0.9695; 0.9652; 0.9607; 0.9559; 0.9509; 0.9458; 0.9446; 0.9436; 0.9426; 0.9417; 0.941; 0.9403; 0.9397; 0.9393; 0.9389; 0.9387; 0.9386; 0.9385; 0.9386; 0.9388; 0.9391; 0.9395; 0.9399; 0.9405; 0.9411; 0.9418; 0.9425; 0.9433; 0.9441; 0.945; 0.9458; 0.9467; 0.9476; 0.9482; 0.9489; 0.9498; 0.9503; 0.9508; 0.951; 0.9511; 0.9512; 0.9509; 0.9506; 0.9497; 0.9489; 0.9473; 0.945; 0.9438; 0.9403; 0.9382; 0.9341; 0.9304]

    
    ## Empty Tank
    # function RAO(amplitude,wavelength, elevation,slope) #ampltiude of the wave, wavelength, elevation and slope of the floater 
    
    #     #expand the vectors to the size of length(ω),length(probx)
    #     amplitude = repeat(amplitude,1,size(elevation)[2])
    #     wavelength = repeat(wavelength,1,size(elevation)[2])
        
        
    #     RAO_η = elevation ./ amplitude
    #     RAO_ηx= slope ./ ((2*π) * amplitude ./ wavelength)
        
    #     return RAO_η, RAO_ηx
    # end
    
    # RAO_η, RAO_ηx= RAO(A_w, λ, η_amplitude, ηₓ_amplitude)

    # return RAO_η, RAO_ηx, η_ϕ, ηₓ_ϕ #, solution_η


    ## Not Empty tank
    function RAO_revised(wavelength, elevation, slope, slopex, slopexx) 
        #ampltiude of the wave, wavelength, elevation and slope of the floater 
    
        #expand the vectors to the size of length(ω),length(probx)
        daRef = wload(datadir("fpv_202312","empt",
            "Empt_length_beam=200_mesh_size=0.25.jld2"))
        amplitude = daRef["RAO_η"]
        # amplitude = repeat(amplitude,1,size(elevation)[2])
        wavelength = repeat(wavelength,1,size(elevation)[2])
        
        
        RAO_η = elevation ./ amplitude
        RAO_ηx= slope ./ ((2*π) * amplitude ./ wavelength)

        RAO_ηxx= slopex ./ ((2*π) * amplitude ./ wavelength) ./ amplitude
        RAO_ηxxx= slopexx ./ ((2*π) * amplitude ./ wavelength) ./ amplitude ./ amplitude
        
        return RAO_η, RAO_ηx, RAO_ηxx, RAO_ηxxx
    end
    
    RAO_η, RAO_ηx, RAO_ηxx, RAO_ηxxx = RAO_revised(λ, η_amplitude, 
        ηx_amplitude, ηxx_amplitude, ηxxx_amplitude)

    return RAO_η, RAO_ηx, η_ϕ, ηx_ϕ, RAO_ηxx, RAO_ηxxx #, solution_η
    
end

#=
ω = [3]
α = [0]
A_w = [1]


#Collect all parameters to a dictionaryc 

allparams = Dict(
    "length_beam" => [200],
    "h_b" => 1,                            
    "γ" => 1,
    "material" => ["hdpe"],#["eps", "cfrp", "gfrp", "neopren", "pvc", "hdpe", "steel"], #store youngs modulus and density 
    "phase"=> [α],
    "omega" => [ω],
    "amplitude" => [A_w] 
)

dicts = dict_list(allparams)
dicts[1]
sim = run_beam(dicts[1])

#plot(ω,A_w,xlims =[0,11], ylims = [0,1.5], title = "constant wave spectrum", label ="ω_start = 1, ω_cutoff = 10")

RAO_η = sim[1]
RAO_ηx = sim[2]
η_phi = sim[3]
ηx_phi = sim[4]

delta_phi = η_phi .- ηx_phi

plot(ω, RAO_η[:,100], title = "RAO_η/RAO_ηx  for center of beam", label = "RAO_η", xlabel = "rad/s")
plot!(ω, RAO_ηx[:,100], label = "RAO_ηx normalized with (4*amplitude/wavelength)", xlabel = "rad/s")


plot(RAO_η[1,:], title = "RAO for the entire beam", label = "omega = 1", xlabel = "m")
plot!(RAO_η[11,:], label = "omega = 2")
plot!(RAO_η[21,:], label = "omega = 3")
 =#
end

