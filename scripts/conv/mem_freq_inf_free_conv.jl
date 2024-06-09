module BeamMultJoints_freq

using Revise
using Gridap
using Gridap.Adaptivity
using Plots
using DrWatson
using WaveSpec
using .Constants
using .Jonswap
using DataFrames:DataFrame
using DataFrames:Matrix
using TickTock

function run_freq(ω, η₀)

  tick()
  # Wave parameters
  k = dispersionRelAng(H0, ω; msg=false)
  λ = 2π/k   
  T = 2π/ω
  @show ω, T
  @show λ
  @show H0, H0/λ

  vxᵢₙ(x) = (η₀*ω)*(cosh(k*(H0 + x[2])) / 
    sinh(k*H0))*exp(im*k*x[1])

  # Weak form
  ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
  a((ϕ,η),(w,v)) =      
    ∫(  ∇(w)⋅∇(ϕ) )dΩ   +
    ∫(  v*(g*η - im*ω*ϕ) +  im*ω*w*η)dΓm  
      #- mᵨ*v*ω^2*η + Tᵨ*∇(v)⋅∇(η) )dΓm  

  l((w,v)) = ∫( w*vxᵢₙ )dΓin + ∫( w*vxᵢₙ )dΓot


  # Solution
  op = AffineFEOperator(a,l,X,Y)
  (ϕₕ,ηₕ) = solve(op)


  # Interpolation on prboes
  
  push!( prbDaΓη, ηₕ(prxΓη) )

  tock()
  return 0
end


name::String = "data/sims_202305_conv6/spec_free"
order::Int = 1
vtk_output::Bool = true
filename = name*"/mem"

ρw = 1025 #kg/m3 water
H0 = 10 #m #still-water depth


# Membrane parameters
@show Lm = 2π #m
@show g #defined in .Constants
@show mᵨ = 0.09 #mass per unit area of membrane / ρw
@show Tᵨ = 0.09*g #T/ρw

ω = 3.13209194927542
η₀ = 1

# Domain 
rf = 64
nx = 9 * rf
ny = 4 * rf
domain =  (0, Lm, -H0, 0.0)
partition = (nx, ny)
xm₀ = 0
xm₁ = xm₀ + Lm
@show Lm
@show domain
@show partition
@show (xm₀, xm₁)
@show Lm/nx
@show H0/ny
println()


model = CartesianDiscreteModel(domain,partition)#, 
  #isperiodic=(true,false))


# Labelling
labels_Ω = get_face_labeling(model)
add_tag_from_tags!(labels_Ω,"surface",[3,4,6])   # assign the label "surface" to the entity 3,4 and 6 (top corners and top side)
add_tag_from_tags!(labels_Ω,"bottom",[1,2,5])    # assign the label "bottom" to the entity 1,2 and 5 (bottom corners and bottom side)
add_tag_from_tags!(labels_Ω,"inlet",[7])         # assign the label "inlet" to the entity 7 (left side)
add_tag_from_tags!(labels_Ω,"outlet",[8])        # assign the label "outlet" to the entity 8 (right side)


# Triangulations
Ω = Interior(model) #same as Triangulation()
Γ = Boundary(model,tags="surface") #same as BoundaryTriangulation()
Γin = Boundary(model,tags="inlet")
Γot = Boundary(model,tags="outlet")


# Masking and Beam Triangulation
Γm = Γ
Γη = Γ


# Construct the tag for membrane boundary
Λmb = Boundary(Γm)

writevtk(model, filename*"_model")
if vtk_output == true
  writevtk(Ω,filename*"_O")
  writevtk(Γm,filename*"_Gm")  
end


# Measures
degree = 2*order
dΩ = Measure(Ω,degree)
dΓm = Measure(Γm,degree)
dΓin = Measure(Γin,degree)
dΓot = Measure(Γot,degree)


# Normals
@show nΛmb = get_normal_vector(Λmb)


# FE spaces
reffe = ReferenceFE(lagrangian,Float64,order)
V_Ω = TestFESpace(Ω, reffe, conformity=:H1, 
  vector_type=Vector{ComplexF64})
V_Γη = TestFESpace(Γη, reffe, conformity=:H1, 
  vector_type=Vector{ComplexF64})
U_Ω = TrialFESpace(V_Ω)
U_Γη = TrialFESpace(V_Γη)
X = MultiFieldFESpace([U_Ω,U_Γη])
Y = MultiFieldFESpace([V_Ω,V_Γη])


# Storing soln at Γ
prxΓη =  [ Point(i,0.0) for i in range(0.0, Lm, 256*9+1) ]
lDa = zeros(ComplexF64, 1, length(prxΓη))
prbDaΓη = DataFrame(lDa, :auto)


# Run weak-form for each freq
run_freq.(ω, η₀)

prbDaΓη = prbDaΓη[2:end,:]

println(prbDaΓη)

data = Dict(
    "prxΓη" => prxΓη,
    "prbDaΓη" => prbDaΓη)

wsave(filename*"_data.jld2", data)

end