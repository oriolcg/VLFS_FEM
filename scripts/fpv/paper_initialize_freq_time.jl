#module init_freq
using DrWatson
@quickactivate "MonolithicFEMVLFS"

using DataFrames
using WaveSpec
using Plots
using .Constants
using .Jonswap
using TimerOutputs
using JLD2
using CSV
using Statistics
using TickTock


# include(srcdir("fpv","paper_phil_freq_time_convert.jl"))
# df = collect_results(datadir("fpv_202401","casePhil","all"))

include(srcdir("fpv","paper_freq_time_convert.jl"))
df = collect_results(datadir("fpv_202402","casePhil","all"))


@show df

Hs = 5.0
Tp = 11.64
d = 30

#freInd = collect(40:1910)
freqInd = collect(40:1140)
#freqInd = [500,1140]

#ω, S, A = JonswapComb.comb_jonswap(Hs, Tp;)


ω, S, A = jonswap(Hs, Tp; plotflag=false,  nω=2401, ωc = 2*pi*2)
dw = ω[2]-ω[1]

ω= ω[freqInd]
A = A[freqInd]
S = S[freqInd]



#spectrum_df = DataFrame(frequency = ω, amplitude = A, wavenumber = k, phaseshift = α, wavelength = λ)
#only use omega >= 1 as RAO only calculated from 1 - 10.
spectrum_dict = Dict(
    "PeakPeriod" =>  Tp,
    "Hs" => Hs,
    "frequency" => ω,
    "amplitude" => A,
    "wavenumber" => dispersionRelAng.(d, ω; msg=false),
    "phaseshift" => randomPhase(ω; seed=100)
)


display(plot(ω,S))


for row in eachrow(df)

    tick()
    result = Freq_time.freq_time_trans(row,
        spectrum_dict)        
    tock()
    println()

end
