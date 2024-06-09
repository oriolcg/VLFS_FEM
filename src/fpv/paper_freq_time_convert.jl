module Freq_time

using DrWatson
@quickactivate "Offshore PV"

using DataFrames
using Plots
using Statistics
using TimerOutputs
using Interpolations
using TickTock
using CSV
export freq_time_trans
using Gridap



function freq_time_trans(rao,#::DataFrameRow, 
    spectrum; 
    dt = 0.25 ,tEnd =1200) 
    #rao should be a DataFrameRow
    
    
    
    Tp = spectrum["PeakPeriod"]
    Hs = spectrum["Hs"]    
    ωe = spectrum["frequency"] 
    A = spectrum["amplitude"]    
    α = spectrum["phaseshift"]


    ## Create a FESpace for derivatives
    # ---------------------Start---------------------
    Lb = 100
    domain = (0, Lb)
    @show nx = round(Int, Lb/rao["mesh_size"])    

    Γη = CartesianDiscreteModel(domain, nx)
    order = 4
    reffe = ReferenceFE(lagrangian,Float64,order)
    V_Γη = FESpace(Γη, reffe, conformity=:H1, 
      vector_type=Vector{ComplexF64})


    xem_cords = 1:2:Lb
    prb_xy = Point.(xem_cords)
    # xemIndex = 2:2:100
    # ----------------------End----------------------


    ## Inteprolate the η Complex values across all ωe freqs
    # ---------------------Start---------------------
    ηdof_org = rao["ηdof"]
    ω_org = rao["omega"]  #contains the omega with the old step size    

    # @show ω_org
    # @show ωe

    ηdof_intp = zeros(Complex, length(ωe), size(ηdof_org,2))
    
    for i in 1:(size(ηdof_intp,2))
        interp_η = linear_interpolation(ω_org, ηdof_org[:,i])
        ηdof_intp[:,i] = interp_η(ωe)
    end    
    # ----------------------End----------------------    


    ## Scale the response as per spectrum amplitude
    A_mat = repeat(A, 1, size(ηdof_intp,2))    
    α_mat = exp.(im*α)
    α_mat = repeat(α_mat, 1, size(ηdof_intp,2))

    ηdof_act = ηdof_intp .* A_mat .* α_mat    
    
    ti = collect(0:dt:tEnd)    
    nTSteps = length(ti)

    ## Time summation at each point
    # function η_at_t(t, η)
    #     sum(η .* exp.(-im*ωe*t))
    # end
    
    # ηdof_time = zeros(Complex, length(ti), size(ηdof_act,2))    
    # tick()
    # for i in axes(ηdof_act,2)
    #     println(i)
    #     ηdof_time[:,i] = η_at_t.(ti, Ref(ηdof_act[:,i]) )
    # end
    # tock()

    phMat = ti * transpose(ωe)
    phMat = exp.(-im*phMat)
        
    ηdof_time = phMat * ηdof_act
    

    η_real_time = zeros(Float64, size(ηdof_time,1), 
        length(prb_xy))
    ηx_real_time = zeros(Float64, size(ηdof_time,1), 
        length(prb_xy))
    
    tick()
    # Method 4: Error Wrong cache
    # The cache_ηx Only depends on the tree strucutre of ηx
    # Therefore it only has to be built once.
    # And this is the quantity that takes time
    # So as Oriol suggested, only compute this once.
    # ηh = FEFunction(V_Γη, ηdof_time[1,:])
    # ηhx = ∇(ηh) ⋅ VectorValue(1.0) 
    # cache_ηx = Gridap.Arrays.return_cache(ηhx, prb_xy)      

    for i in axes(ηdof_time,1)
        ηh = FEFunction(V_Γη, ηdof_time[i,:])
        ηhx = ∇(ηh) ⋅ VectorValue(1.0) 

        # Method 1: 50 sec
        # η_real_time[i,:] = real.(ηh.(prb_xy))
        # ηx_real_time[i,:] = real.(ηhx.(prb_xy))

        # # Method 2: 16 sec
        # cache_η = Gridap.Arrays.return_cache(ηh, prb_xy)      
        # cache_ηx = Gridap.Arrays.return_cache(ηhx, prb_xy)      
        # η_real_time[i,:] = real.(evaluate!(cache_η, ηh, prb_xy))
        # ηx_real_time[i,:] = real.(evaluate!(cache_ηx, ηhx, prb_xy))        

        # Method 2b: 12 sec        
        cache_ηx = Gridap.Arrays.return_cache(ηhx, prb_xy)      
        ηx_real_time[i,:] = real.(evaluate!(cache_ηx, ηhx, prb_xy))        

        # # Method 3: 
        # η_real_time[i,:] = real.(ηh.free_values[xemIndex])
        # ηx_real_time[i,:] = real.(ηhx.free_values[xemIndex])

        # # Method 4: 
        # cache_ηx = Gridap.Arrays.return_cache(ηhx, prb_xy)              
        # ηx_real_time[i,:] = real.(evaluate!(cache_ηx, ηhx, prb_xy))        
        
    end
    tock()

    θ_time = atan.(ηx_real_time)*180/π #transform to degrees        

    # plt1 = plot()
    # tindex = 1:25:1000
    # for i in tindex
    #     plot!(plt1, xem_cords, η_real_time[i,:],dpi=330)
    # end

    save_df_theta = DataFrame(zeros(nTSteps,size(θ_time,2)+2), :auto)
    rename!(save_df_theta, :x1 => :h_b)
    rename!(save_df_theta, :x2 => :material)
    new_names = Symbol.("theta_x" .* string.(1:2:99))
    rename!(save_df_theta, names(save_df_theta)[3:size(θ_time,2)+2] .=> new_names)    

    save_df_theta[!,"h_b"] .= rao["h_b"]
    save_df_theta[!,"material" ] .= rao["material"]
    save_df_theta[!,3:end] = θ_time

        
    # filename = "time_res_mat_" * rao.material *
    #     "_hb_"*string(rao.h_b) * 
    #     "_numFloat_"*string(rao.numFloat)

    # println(filename)
    
    # CSV.write(datadir("fpv_202403")*"/"
    #     *filename*".csv", 
    #     save_df_theta; delim = ';')

    return save_df_theta
    
end


end
