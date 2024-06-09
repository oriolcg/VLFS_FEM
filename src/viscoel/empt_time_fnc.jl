module BeamMultJoints_freq

using Revise
using Gridap
using Printf
using Plots
using DrWatson
using WaveSpec
using .Constants
using .WaveTimeSeries
using .Jonswap
using WriteVTK
using LineSearches: BackTracking
using DataFrames:DataFrame
using DataFrames:Matrix


name::String = "data/sims_202301/empt_time_spec"
order::Int = 2
vtk_output::Bool = true
filename = name*"/mem"


H0 = 10 #m #still-water depth


# Wave parameters
ω, S, η₀ = jonswap(0.4, 2.5; 
    plotflag=true, plotloc=filename, nω=145)
println(ω[1], "\t", ω[2], "\t", ω[end])
ω = ω[2:end]
S = S[2:end]
η₀ = η₀[2:end]
α = randomPhase(ω; seed=100)
# ω = [2*π/2.0, 2*π/2.5, 2*π/4]
# η₀ = [0.3, 0.4, 0.1]
# @show α = randomPhase(ω; seed=100)

# Wave parameters
k = dispersionRelAng.(H0, ω; msg=false)
@show λ = 2π/k
@show T = 2π/ω
ηᵢₙ(x,t) = sum( η₀ .* cos.(k*x[1] - ω*t + α) )
ϕᵢₙ(x,t) = (η₀.*ω./k) .* (cosh.(k*(H0 + x[2])) ./ 
  sinh.(k*H0)) .* sin.(k*x[1]-ω*t + α)
vᵢₙ(x,t) = sum( -(η₀.*ω) .* (cosh.(k*(H0 + x[2])) ./ 
  sinh.(k*H0)) .* cos.(k*x[1]-ω*t + α) )
vzᵢₙ(x,t) = sum( ω .* η₀ .* sin.(k*x[1]-ω*t + α) )
ηᵢₙ(t::Real) = x -> ηᵢₙ(x,t)
ϕᵢₙ(t::Real) = x -> ϕᵢₙ(x,t)
vᵢₙ(t::Real) = x -> vᵢₙ(x,t)
vzᵢₙ(t::Real) = x -> vzᵢₙ(x,t)


# Peak Wave
ωₚ, indp = findmax(η₀)
@show ωₚ = ω[indp]
kₚ = dispersionRelAng(H0, ωₚ; msg=false)
Tₚ = 2*pi/ωₚ
println("Peak Wave T, L ", 2*pi/ωₚ, " ", 2*pi/kₚ)


# Domain 
nx = 1200
ny = 20
mesh_ry = 1.1 #Ratio for Geometric progression of eleSize
Ld = 3*H0 #damping zone length
LΩ = 2*Ld + 3*2*H0
x₀ = -Ld
domain =  (x₀, x₀+LΩ, -H0, 0.0)
partition = (nx, ny)
xdᵢₙ = 0.0
xdₒₜ = x₀ + LΩ - Ld
@show LΩ
@show domain
@show partition
@show LΩ/nx
@show H0/ny
println()


# Time stepping
θₜ = 0.5
t₀ = 0.0
Δt = 2.5/40#Tₚ/40
outΔt = 5.0#Tₚ/2
tf = 2.5*80#50*Tₚ
∂uₜ_∂u = θₜ*Δt


# Numeric constants
h = LΩ / nx
γ = 1.0*order*(order-1)/h
βₕ = 0.5
αₕ = ∂uₜ_∂u/g * (1-βₕ)/βₕ
@show h
@show βₕ
@show αₕ
println()


# Damping
μ₀ = 2.5
μ₁ᵢₙ(x::VectorValue) = μ₀*(1.0 - sin(π/2*(x[1]-x₀)/Ld))
μ₁ₒᵤₜ(x::VectorValue) = μ₀*(1.0 - cos(π/2*(x[1]-xdₒₜ)/Ld))
μ₂ᵢₙ(x) = μ₁ᵢₙ(x)*kₚ
μ₂ₒᵤₜ(x) = μ₁ₒᵤₜ(x)*kₚ
ηd(t) = x -> μ₁ᵢₙ(x)*kₚ*ηᵢₙ(x,t)
∇ₙϕd(t) = x -> μ₁ᵢₙ(x)*vzᵢₙ(x,t)


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
model = CartesianDiscreteModel(domain,partition,map=map)


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


# Auxiliar functions
function is_damping1(xs) # Check if an element is inside the damping zone 1
  n = length(xs)
  x = (1/n)*sum(xs)
  (x₀ <= x[1] <= xdᵢₙ ) * ( x[2] ≈ 0.0)
end
function is_damping2(xs) # Check if an element is inside the damping zone 2
  n = length(xs)
  x = (1/n)*sum(xs)
  (xdₒₜ <= x[1] ) * ( x[2] ≈ 0.0)
end

# Masking and Beam Triangulation
xΓ = get_cell_coordinates(Γ)
Γd1_to_Γ_mask = lazy_map(is_damping1, xΓ)
Γd2_to_Γ_mask = lazy_map(is_damping2, xΓ)
Γd1 = Triangulation(Γ, findall(Γd1_to_Γ_mask))
Γd2 = Triangulation(Γ, findall(Γd2_to_Γ_mask))
Γfs = Triangulation(Γ, findall(!, Γd1_to_Γ_mask .| Γd2_to_Γ_mask))
Γκ = Γ


writevtk(model, filename*"_model")
if vtk_output == true
  writevtk(Ω,filename*"_O")
  writevtk(Γ,filename*"_G")
  writevtk(Γd1,filename*"_Gd1")
  writevtk(Γd2,filename*"_Gd2")
  writevtk(Γfs,filename*"_Gfs")
end


# Measures
degree = 2*order
dΩ = Measure(Ω,degree)
dΓd1 = Measure(Γd1,degree)
dΓd2 = Measure(Γd2,degree)
dΓfs = Measure(Γfs,degree)
dΓin = Measure(Γin,degree)


# FE spaces
reffe = ReferenceFE(lagrangian,Float64,order)
V_Ω = TestFESpace(Ω, reffe, conformity=:H1)
V_Γκ = TestFESpace(Γκ, reffe, conformity=:H1)
U_Ω = TransientTrialFESpace(V_Ω)
U_Γκ = TransientTrialFESpace(V_Γκ)
X = TransientMultiFieldFESpace([U_Ω,U_Γκ])
Y = MultiFieldFESpace([V_Ω,V_Γκ])


# Weak form
∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
m((ϕₜ,κₜ),(w,u)) = 
  ∫(  βₕ*(u + αₕ*w)*ϕₜ - w*κₜ )dΓfs +
  ∫(  βₕ*(u + αₕ*w)*ϕₜ - w*κₜ )dΓd1    +
  ∫(  βₕ*(u + αₕ*w)*ϕₜ - w*κₜ )dΓd2   
a((ϕ,κ),(w,u)) =      
  ∫(  ∇(w)⋅∇(ϕ) )dΩ   +
  ∫(  βₕ*(u + αₕ*w)*g*κ )dΓfs   +
  ∫(  βₕ*(u + αₕ*w)*g*κ - μ₂ᵢₙ*κ*w + μ₁ᵢₙ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd1  +
  ∫(  βₕ*(u + αₕ*w)*g*κ - μ₂ₒᵤₜ*κ*w + μ₁ₒᵤₜ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd2  

l(t,(w,u)) =  
  ∫( w * vᵢₙ(t) )dΓin - 
  ∫( ηd(t)*w - ∇ₙϕd(t)*(u + αₕ*w) )dΓd1

# Solution
op = TransientConstantMatrixFEOperator(m, a, l, X, Y)
ls = LUSolver()
ode_solver = ThetaMethod(ls,Δt,θₜ)

# Initial solution
u0 = interpolate_everywhere([0.0,0.0],X(0.0))

@show uht = solve(ode_solver,op,u0,t₀,tf)


# Probes
# prbx=range(-30, 90, 13)
# @show prbxy = Point.(prbx, 0.0)
lDa = zeros(Float64, 1, 2)
prbDa = DataFrame(lDa, :auto) # for η


if vtk_output == true
  pvd_Ω = paraview_collection(filename * "_O_sol", append=false)
  pvd_Γκ = paraview_collection(filename * "_Gk_sol", append=false)
end

if vtk_output == true
  tpr = @sprintf("%5.3f",t₀)                    
  tval = @sprintf("%d",round(Int64,t₀*1000))
  ϕₕ, κₕ = u0
  pvd_Ω[t₀] = createvtk(Ω,
    filename * "_O_sol" * "_$tval.vtu",
    cellfields = ["phi" => ϕₕ])
  pvd_Γκ[t₀] = createvtk(Γκ,
    filename * "_Gk_sol" * "_$tval.vtu",
    cellfields = ["kappa" => κₕ])  
end    


# Execute
@show outMod = round(Int64,outΔt/Δt);

for (uh, t) in uht  
    ϕₕ, κₕ = uh
    tpr = @sprintf("%5.3f",t)                    
    tval = @sprintf("%d",round(Int64,t*1000))    

    lDa[1] = t
    lDa[2] = κₕ.free_values[601] #middle-node
    @show lDa
    push!(prbDa, lDa)


    if(round(t/Δt) % outMod != 0)
      println("Time : $tpr \t", round(t/Δt) / outMod)
      continue
    end

    println("Time : $tpr \t", round(t/Δt) / outMod, "vtk")
    
    if vtk_output == true
      pvd_Ω[t] = createvtk(Ω,
        filename * "_O_sol" * "_$tval.vtu",
        cellfields = ["phi" => ϕₕ])
      pvd_Γκ[t] = createvtk(Γκ,
        filename * "_Gk_sol" * "_$tval.vtu",
        cellfields = ["kappa" => κₕ])      
    end    
end

if vtk_output == true
  vtk_save(pvd_Ω)
  vtk_save(pvd_Γκ)
end


data = Dict("prbDa" => prbDa,
            "ω" => ω,
            "α" => α,
            "η₀" => η₀)

wsave(filename*"_data.jld2", data)
end