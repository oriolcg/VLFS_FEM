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
using Printf


include(srcdir("fpv","paper_freq_time_convert.jl"))


daHs = CSV.read(
  "data/fpv_202403/dataHsTmHourly/sign_wave_height.csv",
  DataFrame )

daTz = CSV.read(
  "data/fpv_202403/dataHsTmHourly/average_period.csv",
  DataFrame )

daDay = CSV.read(
  "data/fpv_202403/dataHsTmHourly/daylight_flag.csv",
  DataFrame )


da = wload("data/fpv_202403/caseHDPE_F100_T0p20/depth=30_h_b=0.2_length_beam=100_material=hdpe_mesh_size=1.0_numFloat=1.jld2")
fileName = datadir("fpv_202403", "caseHDPE_F100_T0p20", "PLA-NF01")*"/"


function rampWindow(t, tRamp)  

  t0 = t[1]
  t1 = t[1]+tRamp
  t2 = t[end]-tRamp
  t3 = t[end]

  rampt1 = ifelse.( t .< t1, 0.5*(1 .- cos.(pi*(t.-t0)/tRamp)), 1 )
  rampt2 = ifelse.( t .> t2, 0.5*(1 .- cos.(pi*(t3.-t)/tRamp)), 1 )

  # for i in axes(ramp,1)
  #   println(t[i], "\t", rampt1[i], "\t", rampt2[i])
  # end

  flag = rampt2 .>= rampt1

  ramp = rampt1.* flag + rampt2 .* (1 .- flag)
  
end


dt = 0.25
tEnd = 3600
rampT = 15 #s each end
t = 0:dt:tEnd
tRamp = rampWindow(t, rampT)


nHours = size(daHs,1)

# for i in 7000:nHours
for i in axes(daHs,1)
  
  println("-x-x-x-x- ")
  @printf("Hour = %i \n", i)
  @printf("Percentage = %0.2f %%\n", i/nHours*100)
  println(fileName)

  Hs = daHs[i,1]
  Tz = daTz[i,1]
  day = daDay[i,2]
  d = da["depth"]  
  @show Hs, Tz, day


  Gamma = 3.3
  # Taking same peakedness para as Alba's previous work

  Tp = Tz*(0.327*exp(-0.315*Gamma)+1.17);
  ωp = 2*pi/Tp


  # I have put cut off freq at 12 rad/s
  # Alba's data has some specs with very small Hs and Tz
  # For example certain specs have Tz = 0.2s
  # Tz = 0.2s => Tp = 0.26s => ωp = 24 rad/s > ωc
  # Our response beyond 6 rad/s is anyways not important
  # So I reset ωp to 10.0 in these cases

  if(ωp > 10.0)
    Tp = 2*pi/10.0
    ω, S, A = jonswap(Hs, Tp, Gamma; 
      plotflag=false,  nω=2401, ωc = 2*pi*2)      
  
  else
    ω, S, A = jonswap(Hs, Tp, Gamma; 
      plotflag=false,  nω=2401, ωc = 2*pi*2)  

  end

  # ω, S, A = jonswap(Hs, Tp, Gamma; 
  #   plotflag=false,  nω=2401, ωc = 2*pi*2)  

  freqInd = findall( 0.2.< ω .<6.0 )
  ω= ω[freqInd]
  A = A[freqInd]
  S = S[freqInd]

  if(day == 0)
    println("--Night---")
    println()
    continue
  end
  println("----Day----")


  #spectrum_df = DataFrame(frequency = ω, amplitude = A, wavenumber = k, phaseshift = α, wavelength = λ)
  #only use omega >= 1 as RAO only calculated from 1 - 10.
  spectrum_dict = Dict(
      "PeakPeriod" =>  Tp,
      "Hs" => Hs,
      "frequency" => ω,
      "amplitude" => A,
      "wavenumber" => dispersionRelAng.(d, ω; msg=false),
      "phaseshift" => randomPhase(ω)
  )

  result = Freq_time.freq_time_trans(da, spectrum_dict;
    dt=dt, tEnd=tEnd)  

  tRampMat = repeat(tRamp, 1, size(result,2)-2)
  vals = Matrix( result[:,3:end] )
  vals = vals .* tRampMat
  vals = round.(vals, digits=4)
  result[:,3:end] = vals 
  
  saveName = fileName * "hour_"*@sprintf("%04i",i)*".csv"

  CSV.write(saveName, result; delim = ';')

  println()

end