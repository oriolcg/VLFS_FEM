using Gridap
using Gridap.Io
using GridapGmsh


# if GridapGmsh.gmsh.isInitialized()
  # GridapGmsh.gmsh.finalize()
# end
# infilename = joinpath(@__DIR__,"floating_ice_coarse.msh")
# outfilename = joinpath(@__DIR__,"floating_ice_coarse.json")
infilename = joinpath(@__DIR__,"floating_ice-step_ratio05_refined.msh")
outfilename = joinpath(@__DIR__,"floating_ice-step_ratio05_refined.json")
model = GmshDiscreteModel(infilename)
to_json_file(model,outfilename)
