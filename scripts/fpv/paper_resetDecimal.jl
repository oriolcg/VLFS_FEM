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



case = "PLA-NF50"

readName = datadir("fpv_202403", "caseHDPE_F100_T0p20", case)*"/"

fileName = datadir("fpv_202403", "caseHDPE_F100_T0p20", "new", case)*"/"


all_files = readdir(readName)

for iName in all_files
  println(iName)

  df = CSV.read(readName*iName, DataFrame)

  vals = Matrix(df[:,3:end])
  vals = round.(vals, digits=4)

  df[:,3:end] = vals
  CSV.write(fileName*iName, df, delim = ";")

end
