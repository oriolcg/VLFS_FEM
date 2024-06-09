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




#This file contains a function, which takes the DataframeRow as input and a random spectrum also stored in a Dataframe as an input
#Output is a time solution for θ
#Output is a DataFrame containing all the variations of parameters in the rows, while the different columns represent diffrent input spectrum


function freq_time_trans(rao::DataFrameRow, 
    spectrum; 
    time_steps = 0.25 ,seconds::Int64 =1200) 
    #rao should be a DataFrameRow
    
    #rao is the response amplitude function for a constant wavespectrum A = 1 and intial phaseshift = 0 
    
    Tp = spectrum["PeakPeriod"]
    Hs = spectrum["Hs"]
    ω = spectrum["frequency"]
    k = spectrum["wavenumber"]
    A = spectrum["amplitude"]    
    α = spectrum["phaseshift"]

    xi = ω

    slope_spectrum = k .* A  #Wave Height over wave length

    RAO_η  = rao["RAO_η"] 
    RAO_ηx = rao["RAO_ηx"]
    ϕ_η = rao["η_ϕ"]
    ϕ_ηx = rao["ηx_ϕ"]   

    # Taking only x=1.0, 3.0, 5.0 and so on
    RAO_η  = RAO_η[:, 2:2:end]
    RAO_ηx = RAO_ηx[:, 2:2:end]
    ϕ_η = ϕ_η[:, 2:2:end]
    ϕ_ηx = ϕ_ηx[:, 2:2:end]


    ω_org = rao["omega"]  #contains the omega with the old step size    
    ω_org = repeat(ω_org, 1, size(RAO_ηx)[2])
    #RAOs need to be interpolated 

    #xi = range(0, 2*pi*2, 2401)[freqInd] #Nw 2401 between 0 and 2*pi*2hz -> just take values between 0.2 and 10 rad/s
    RAOη_new = zeros(length(xi) ,size(RAO_η)[2])
    RAOηx_new = zeros(length(xi) ,size(RAO_ηx)[2])
    ϕ_η_new = zeros(length(xi) ,size(ϕ_η)[2])
    ϕ_ηx_new = zeros(length(xi) ,size(ϕ_ηx)[2])

    tick()
    for i in 1:(size(RAOη_new)[2])
        inter_η = linear_interpolation(ω_org[:,i], RAO_η[:,i])
        RAOη_new[:,i] = inter_η(xi)
        inter_ηx = linear_interpolation(ω_org[:,i], RAO_ηx[:,i])
        RAOηx_new[:,i] = inter_ηx(xi)
        inter_ϕ_η = linear_interpolation(ω_org[:,i], ϕ_η[:,i])
        ϕ_η_new[:,i] = inter_ϕ_η(xi) #Philipp mistake!!!
        inter_ϕ_ηx = linear_interpolation(ω_org[:,i], ϕ_ηx[:,i])
        ϕ_ηx_new[:,i] = inter_ϕ_ηx(xi)
    end
    tock()
    #From now on the RAO have more rows, columns represent the postions
    #Overwrite to not change the code
    #RAO_η = RAOη_new
    #RAO_ηx = RAOηx_new

    println(size(RAO_ηx))
    
    A = repeat(A,1,size(RAOη_new)[2])
    slope_spectrum = repeat(slope_spectrum,1,size(RAOηx_new)[2])
    k = repeat(k,1,size(RAOη_new)[2])
    ω = repeat(ω,1,size(RAOη_new)[2])
    α = repeat(α,1,size(RAOη_new)[2])


    println(size(A))
    elevation_amp = A .* RAOη_new 
    elevation_angle = ϕ_η_new .+ α

    slope_amp = slope_spectrum .* RAOηx_new
    slope_angle = ϕ_ηx_new .+ α
    
  
    # #Functions to calculate the sume  Philip incorrect
    # fourier_η(t) =  sum(elevation_amp .* sin.(ω*t.+ elevation_angle),dims=1) #function to sum up all fourier parts for every point for a specific time
    # #fourier_ηᵪ(x) = sum(slope_amp .* sin.(ω*x-k * x .+ slope_angle),dims=1)
    # fourier_ηᵪ(x) = sum(slope_amp .* sin.(ω*x .+ slope_angle),dims=1) #right formulation

    #Functions to calculate the sume 
    fourier_η(t) =  sum(elevation_amp .* 
        cos.(-ω*t.+ elevation_angle), dims=1) #function to sum up all fourier parts for every point for a specific time    
    # fourier_ηᵪ(t) = sum( -slope_amp .* 
    #     sin.(-ω*t.+ elevation_angle), dims=1) 
    fourier_ηᵪ(t) = sum( slope_amp .* 
        cos.(-ω*t.+ slope_angle), dims=1) 

    time_points = collect(0:time_steps:seconds)
    vec(time_points)

    # η = vcat(fourier_η.(time_points)...) #not really needed
    ηₓ= vcat(fourier_ηᵪ.(time_points)...)

    θ_points = atan.(ηₓ)*180/π #transform to degrees
    
    # θ_points = η
    # θ = ones(length(time_points),size(RAOη_new)[2]-1)

    # for i in 1:(size(RAOη_new)[2]-1)
    #     θ[:,i] = (θ_points[:,i] .+ θ_points[:,i+1]) ./ 2 #average the tilt between to points
    # end

    θ = θ_points

    Df_time = DataFrame("seconds" => time_points)
    #Df_η    = DataFrame(η,:auto)
    
    Df_θ    = DataFrame(θ,:auto)
    
    #calculate statistical values

    Df_mean_θ = mean.(eachcol(Df_θ))
    Df_std_θ = [std(Df_θ[!,col], corrected = true) for col in names(Df_θ)]


    θ_time = [Df_time Df_θ]
    #η_time = [Df_time Df_η]

    rao = hcat(DataFrame(rao), [θ_time])
    rename!(rao, :x1 => :theta_time)
    #rao = hcat(DataFrame(rao), [η_time])
    #rename!(rao, :x1 => :η_time)
    rao = hcat(DataFrame(rao), [Df_mean_θ])
    rename!(rao, :x1 => :theta_mean)
    rao = hcat(DataFrame(rao), [Df_std_θ])
    rename!(rao, :x1 => :theta_std)


    rao = select!(rao, Not([:amplitude, :omega, :length_beam, :RAO_ηx, :RAO_η, :η_ϕ, :ηx_ϕ,  :phase, :path ]))
                                     
    println("THETA")
    # 4800 is seconds/dt
    save_df_theta = DataFrame(zeros(4801,size(θ,2)+2), :auto)
    rename!(save_df_theta, :x1 => :h_b)
    rename!(save_df_theta, :x2 => :material)
    new_names = Symbol.("theta_x" .* string.(1:2:99))
    rename!(save_df_theta, names(save_df_theta)[3:size(θ,2)+2] .=> new_names)    

    
    save_df_theta[!,"h_b"] .= rao.h_b[1]
    save_df_theta[!,"material" ] .= rao.material[1]
    save_df_theta[!,3:end] = θ
    
    filepath = "fpv_pp_" * string(Tp) *"_Hs_"*string(Hs)
    
    filename = "/time_res_mat_" * rao.material[1] *
        "_hb_"*string(rao.h_b[1]) * 
        "_numFloat_"*string(rao.numFloat[1])
    
    CSV.write(datadir("fpv_202401")*"/"
        *filepath*filename*".csv", 
        save_df_theta; delim = ';')

    
end


end
