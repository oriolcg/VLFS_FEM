module analysis

using Gridap
using DrWatson

d1 = wload("data/sims_202212/BeamMultJoints-freq-method1/dom_data.jld2")
#d2 = wload("data/sims_202212/BeamMultJoints-freq-method2/dom_data.jld2")
d2 = wload("data/sims_202212/BeamMultJoints-freq-method1r2/dom_data.jld2")

v1 = d1["κₕ"]
v2 = d2["κₕ"]

@show norm(v1.cell_dof_values - v2.cell_dof_values)

v1 = d1["ηₕ"]
v2 = d2["ηₕ"]

@show norm(v1.cell_dof_values - v2.cell_dof_values)

v1 = d1["ϕₕ"]
v2 = d2["ϕₕ"]

@show norm(v1.cell_dof_values - v2.cell_dof_values)

end