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

function run_freq(ω, η₀, α)

  tick()
  # Wave parameters
  k = dispersionRelAng(H0, ω; msg=false)
  λ = 2π/k   
  T = 2π/ω
  ηᵢₙ(x) = η₀*exp(im*k*x[1] + im*α)
  ϕᵢₙ(x) = -im*(η₀*ω/k)*(cosh(k*(H0 + x[2])) / 
    sinh(k*H0))*exp(im*k*x[1] + im*α)
  vxᵢₙ(x) = (η₀*ω)*(cosh(k*(H0 + x[2])) / 
    sinh(k*H0))*exp(im*k*x[1] + im*α)
  vzfsᵢₙ(x) = -im*ω*η₀*exp(im*k*x[1] + im*α) #???
  @show ω, T
  @show λ
  @show η₀
  @show H0, H0/λ

  # Numeric constants
  αₕ = -im*ω/g * (1-βₕ)/βₕ
  @show αₕ
  println()  

  # Damping
  μ₀ = 2.5#maximum([2.5, 5.24/(ω^0.922)])#2.5
  μ₁ᵢₙ(x) = μ₀*(1.0 - sin(π/2*(x[1]-x₀)/Ld))
  μ₂ᵢₙ(x) = μ₁ᵢₙ(x)*k
  μ₂ₒᵤₜ(x) = μ₁ₒᵤₜ(x)*k
  ηd(x) = μ₂ᵢₙ(x)*ηᵢₙ(x)
  ∇ₙϕd(x) = μ₁ᵢₙ(x)*vzfsᵢₙ(x) #???

  # Weak form
  ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
  a((ϕ,κ,η),(w,u,v)) =      
    ∫(  ∇(w)⋅∇(ϕ) )dΩ   +
    ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ )dΓfs   +
    ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ 
      - μ₂ᵢₙ*κ*w + μ₁ᵢₙ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd1    +
    ∫( -w * im * k * ϕ )dΓot +
    # ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ 
    #   - μ₂ₒᵤₜ*κ*w + μ₁ₒᵤₜ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd2    +
    ∫(  v*(g*η - im*ω*ϕ) +  im*ω*w*η
      - mᵨ*v*ω^2*η + Tᵨ*(1-im*ω*τ)*∇(v)⋅∇(η) )dΓm  #+ 
    #∫(- Tᵨ*(1-im*ω*τ)*v*∇(η)⋅nΛmb )dΛmb

  l((w,u,v)) =  ∫( w*vxᵢₙ )dΓin - ∫( ηd*w - ∇ₙϕd*(u + αₕ*w) )dΓd1


  # Solution
  op = AffineFEOperator(a,l,X,Y)
  (ϕₕ,κₕ,ηₕ) = solve(op)

  # Function for inlet phase
  κin = interpolate_everywhere(ηᵢₙ, 
    FESpace(Γκ, reffe, conformity=:H1, vector_type=Vector{ComplexF64}))
  κr = κₕ - κin

  # Energy flux (Power) calculation
  ηx = ∇(ηₕ)⋅VectorValue(1.0,0.0)
  Pd = sum(∫( abs(ηx)*abs(ηx) )dΓm)
  Pd = 0.5*Tᵨ*ρw*τ*ω*ω*Pd

  # Wave energy flux
  ηrf = abs(κr(Point(prbPowx[1],0.0)))
  ηtr = abs(κₕ(Point(prbPowx[2],0.0)))
  kh = k*H0
  wave_n = 0.5*(1 + 2*kh/sinh(2*kh))
  Pin = (0.5*ρw*g*η₀*η₀)*(ω/k)*wave_n
  Prf = (0.5*ρw*g*ηrf*ηrf)*(ω/k)*wave_n
  Ptr = (0.5*ρw*g*ηtr*ηtr)*(ω/k)*wave_n
  PErr = Pin - Prf - Ptr - Pd
  println("Power In \t ",Pin," W/m")
  println("Power Ref \t ",Prf," W/m")
  println("Power Trans \t ",Ptr," W/m")
  println("Power Abs \t ",Pd," W/m")
  println("Error \t ",PErr," W/m")

  # Interpolation on prboes
  
  push!( prbDaΓη, ηₕ(prxΓη) )

  push!(prbPow, [Pin, Prf, Ptr, Pd, PErr, 0.0])
  
  tock()
  return 0
end


name::String = "data/sims_202305_conv3/spec_free"
order::Int = 2
vtk_output::Bool = true
filename = name*"/mem"

ρw = 1025 #kg/m3 water
H0 = 10 #m #still-water depth

# Wave parameters
# ω, S, η₀ = jonswap(0.4, 2.5; 
#     plotflag=true, plotloc=filename, nω=145)
# println(ω[1], "\t", ω[2], "\t", ω[end])
# ω = ω[2:end]
# S = S[2:end]
# η₀ = η₀[2:end]
# ω = [2*π/2.53079486745378, 2*π/2.0]
# η₀ = [0.25, 0.25]
ω = 2.0:0.1:2.1
T = 2*π./ω
η₀ = 0.25*ones(length(ω))
α = randomPhase(ω; seed=100)
k = dispersionRelAng.(H0, ω; msg=false)


# Peak Wave
ωₚ, indp = findmax(η₀)
@show ωₚ = ω[indp]
kₚ = dispersionRelAng(H0, ωₚ; msg=false)
println("Peak Wave T, L ", 2*pi/ωₚ, " ", 2*pi/kₚ)


# Membrane parameters
@show Lm = 2*H0 #m
@show g #defined in .Constants
@show mᵨ = 0.9 #mass per unit area of membrane / ρw
@show Tᵨ = 0.1*g*H0*H0 #T/ρw
@show τ = 0.0#damping coeff


# Domain 
nx = 170
ny = 4
mesh_ry = 1.0 #Ratio for Geometric progression of eleSize
Ld = 5*H0 #damping zone length
LΩ = 12*H0 + Ld #2*Ld
x₀ = -Ld
domain =  (x₀, x₀+LΩ, -H0, 0.0)
partition = (nx, ny)
xdᵢₙ = 0.0
# xdₒₜ = x₀ + LΩ - Ld
xm₀ = xdᵢₙ + 5*H0
xm₁ = xm₀ + Lm
@show Lm
@show LΩ
@show domain
@show partition
@show (xm₀, xm₁)
@show isinteger(Lm/LΩ*nx)
@show LΩ/nx
@show H0/ny
println()


# Numeric constants
h = LΩ / nx
γ = 1.0*order*(order-1)/h
βₕ = 0.5
# αₕ = -im*ω/g * (1-βₕ)/βₕ
@show h
@show βₕ
# @show αₕ
println()


# Damping
# μ₀ = 2.5
# μ₁ᵢₙ(x) = μ₀*(1.0 - sin(π/2*(x[1]-x₀)/Ld))
# μ₁ₒᵤₜ(x) = μ₀*(1.0 - cos(π/2*(x[1]-xdₒₜ)/Ld))
# μ₂ᵢₙ(x) = μ₁ᵢₙ(x)*kₚ
# μ₂ₒᵤₜ(x) = μ₁ₒᵤₜ(x)*kₚ



# Mesh
function f_y(y, r, n, H0)
  # Mesh along depth as a GP
  # Depth is 0 to -H0
  if(r ≈ 1.0)
    return y  
  else
    a0 = H0 * (r-1) / (r^n - 1)    
    if y ≈ 0
      return 0.0
    end
    j = abs(y) / H0 * n  
    return -a0 / (r-1) * (r^j - 1)
  end
end
map(x) = VectorValue( x[1], f_y(x[2], mesh_ry, ny, H0) )
modelo = CartesianDiscreteModel(domain,partition,map=map)
model = get_model( refine(modelo,32) )


# Labelling
labels_Ω = get_face_labeling(model)
add_tag_from_tags!(labels_Ω,"surface",[3,4,6])   # assign the label "surface" to the entity 3,4 and 6 (top corners and top side)
add_tag_from_tags!(labels_Ω,"bottom",[1,2,5])    # assign the label "bottom" to the entity 1,2 and 5 (bottom corners and bottom side)
add_tag_from_tags!(labels_Ω,"inlet",[7])         # assign the label "inlet" to the entity 7 (left side)
add_tag_from_tags!(labels_Ω,"outlet",[8])        # assign the label "outlet" to the entity 8 (right side)
add_tag_from_tags!(labels_Ω, "water", [9])       # assign the label "water" to the entity 9 (interior)


# Triangulations
Ω = Interior(model) #same as Triangulation()
Γ = Boundary(model,tags="surface") #same as BoundaryTriangulation()
Γin = Boundary(model,tags="inlet")
Γot = Boundary(model,tags="outlet")


# Auxiliar functions
function is_mem(xs) # Check if an element is inside the beam1
  n = length(xs)
  x = (1/n)*sum(xs)
  (xm₀ <= x[1] <= xm₁ ) * ( x[2] ≈ 0.0)
end
function is_damping1(xs) # Check if an element is inside the damping zone 1
  n = length(xs)
  x = (1/n)*sum(xs)
  (x₀ <= x[1] <= xdᵢₙ ) * ( x[2] ≈ 0.0)
end
# function is_damping2(xs) # Check if an element is inside the damping zone 2
#   n = length(xs)
#   x = (1/n)*sum(xs)
#   (xdₒₜ <= x[1] ) * ( x[2] ≈ 0.0)
# end

# Masking and Beam Triangulation
xΓ = get_cell_coordinates(Γ)
Γm_to_Γ_mask = lazy_map(is_mem, xΓ)
Γd1_to_Γ_mask = lazy_map(is_damping1, xΓ)
#Γd2_to_Γ_mask = lazy_map(is_damping2, xΓ)
Γm = Triangulation(Γ, findall(Γm_to_Γ_mask))
Γd1 = Triangulation(Γ, findall(Γd1_to_Γ_mask))
#Γd2 = Triangulation(Γ, findall(Γd2_to_Γ_mask))
Γfs = Triangulation(Γ, findall(!, Γm_to_Γ_mask .| 
  Γd1_to_Γ_mask ))# .| Γd2_to_Γ_mask))
Γη = Triangulation(Γ, findall(Γm_to_Γ_mask))
Γκ = Triangulation(Γ, findall(!,Γm_to_Γ_mask))


# Construct the tag for membrane boundary
Λmb = Boundary(Γm)
xΛmb = get_cell_coordinates(Λmb)
xΛmb_n1 = findall(model.grid_topology.vertex_coordinates .== xΛmb[1])
xΛmb_n2 = findall(model.grid_topology.vertex_coordinates .== xΛmb[2])
new_entity = num_entities(labels_Ω) + 1
labels_Ω.d_to_dface_to_entity[1][xΛmb_n1[1]] = new_entity
labels_Ω.d_to_dface_to_entity[1][xΛmb_n2[1]] = new_entity
add_tag!(labels_Ω, "mem_bnd", [new_entity])


writevtk(model, filename*"_model")
if vtk_output == true
  writevtk(Ω,filename*"_O")
  writevtk(Γ,filename*"_G")
  writevtk(Γm,filename*"_Gm")  
  writevtk(Γd1,filename*"_Gd1")
  # writevtk(Γd2,filename*"_Gd2")
  writevtk(Γfs,filename*"_Gfs")
  writevtk(Λmb,filename*"_Lmb")  
end


# Measures
degree = 2*order
dΩ = Measure(Ω,degree)
dΓm = Measure(Γm,degree)
dΓd1 = Measure(Γd1,degree)
# dΓd2 = Measure(Γd2,degree)
dΓfs = Measure(Γfs,degree)
dΓin = Measure(Γin,degree)
dΓot = Measure(Γot,degree)
dΛmb = Measure(Λmb,degree)


# Normals
@show nΛmb = get_normal_vector(Λmb)


# Dirichlet Fnc
gη(x) = ComplexF64(0.0)

# FE spaces
reffe = ReferenceFE(lagrangian,Float64,order)
V_Ω = TestFESpace(Ω, reffe, conformity=:H1, 
  vector_type=Vector{ComplexF64})
V_Γκ = TestFESpace(Γκ, reffe, conformity=:H1, 
  vector_type=Vector{ComplexF64})
# V_Γη = TestFESpace(Γη, reffe, conformity=:H1, 
#   vector_type=Vector{ComplexF64},
#   dirichlet_tags=["mem_bnd"]) #diri
V_Γη = TestFESpace(Γη, reffe, conformity=:H1, 
  vector_type=Vector{ComplexF64})
U_Ω = TrialFESpace(V_Ω)
U_Γκ = TrialFESpace(V_Γκ)
# U_Γη = TrialFESpace(V_Γη, gη)
U_Γη = TrialFESpace(V_Γη)
X = MultiFieldFESpace([U_Ω,U_Γκ,U_Γη])
Y = MultiFieldFESpace([V_Ω,V_Γκ,V_Γη])

prbPowx=[ 25.0, 95.0 ]


# Storing soln at Γ
prxΓη =  [ Point(i,0.0) for i in range(50.0, 70.0, 128*20+1) ]
lDa = zeros(ComplexF64, 1, length(prxΓη))
prbDaΓη = DataFrame(lDa, :auto)

prbPow = DataFrame(zeros(Float64, 1, 6), :auto)

# Run weak-form for each freq
run_freq.(ω, η₀, α)

prbDaΓη = prbDaΓη[2:end,:]
prbPow = prbPow[2:end,:]

k = dispersionRelAng.(H0, ω; msg=false)

data = Dict("ω" => ω,
            "η₀" => η₀,
            "k" => k,
            "prxΓη" => prxΓη,
            "prbDaΓη" => prbDaΓη,
            "prbPow" => prbPow)

wsave(filename*"_data.jld2", data)

end