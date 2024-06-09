#module init_freq
using DrWatson
@quickactivate "MonolithicFEMVLFS"

using DataFrames
using WaveSpec
using Plots
using .Constants
using .Jonswap
using .WaveTimeSeries
using TimerOutputs
using JLD2
using CSV
using Statistics
using TickTock


#using .Freq_time
include(srcdir("fpv","freq_time_convert.jl"))


#load the simulation data 

df = collect_results(datadir("fpv_202311"))
@show df

Hs = 5.0
Tp = 11.64
d = 30

#ω, S, A = JonswapComb.comb_jonswap(Hs, Tp;)


ω, S, A = jonswap(Hs, Tp; plotflag=false,  nω=2401, ωc = 2*pi*2)
dw = ω[2]-ω[1]

ω= ω[40:1910]
A = A[40:1910]
S = S[40:1910]
xi = range(0, 2*pi*2, 2401)[40:1910]

#spectrum_df = DataFrame(frequency = ω, amplitude = A, wavenumber = k, phaseshift = α, wavelength = λ)
#only use omega >= 1 as RAO only calculated from 1 - 10.
spectrum_df = DataFrame(Dict("Peak Period" =>  Tp,
                            "significant waveheight" => Hs,
                            "frequency" => ω,
                            "amplitude" => A,
                            "wavenumber" => dispersionRelAng.(d, ω; msg=false),
                            "phaseshift" => randomPhase(ω; seed=100),
                            "wavelength" => 2*π ./ dispersionRelAng.(d, ω; msg=false)))


display(plot(ω,S))


#Initialize a Dataframe. First column is h_b, second:Length, third = material, then 200 columns with mean values, 
#Dataframe was not used in the End
matrix = zeros(size(df)[1],403)
save_df = DataFrame(matrix, :auto)
rename!(save_df, :x1 => :h_b)
rename!(save_df, :x2 => :L_b)
rename!(save_df, :x3 => :material)
new_names_std = Symbol.("std_x" .* string.(1:200))
new_names_mean = Symbol.("mean_x" .* string.(1:200))
rename!(save_df, names(save_df)[4:203] .=> new_names_mean)
rename!(save_df, names(save_df)[204:end] .=> new_names_std)
names(save_df)
save_df.material = string.(save_df.material)
size(save_df[:,4:203])
size(save_df[:,204:end])

global i = 1

for row in eachrow(df)

    tick()
    result = Freq_time.freq_time_trans(row,spectrum_df, Hs, Tp)        
    tock()

end

#filename = "Test file new dw"
#filename = "Time results"*" pp="*string(Tp)*" Hs="*string(Hs) 
#CSV.write("csv data/CSV time results/"*filename *".csv", save_df; delim = ';')


#df_final.theta_time[1].x1

#df_final = DataFrame(freq_time_trans(df[1,:],spectrum_df))
#histogram(df_final.theta_time[1].x1, bins = 100)

#= 
global df_final = DataFrame(freq_time_trans(df[1,:],spectrum_df))

for i in range(2,3) #takes 11 seconds
    df_final = vcat(df_final, freq_time_trans((df[i,:]),spectrum_df))
end
 =#


#end



