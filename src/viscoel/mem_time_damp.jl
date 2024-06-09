module BeamMultJoints_freq

using Revise
using Gridap
using Printf
using Plots
using DrWatson
using WaveSpec
using .Constants
using .WaveTimeSeries
using WriteVTK


name::String = "data/sims_202306/run/mono_time"
order::Int = 2
vtk_output::Bool = true
filename = name*"/mem"

ρw = 1025 #kg/m3 water
H0 = 10 #m #still-water depth

# Membrane parameters
@show Lm = 2*H0 #m
@show g #defined in .Constants
@show mᵨ = 0.9 #mass per unit area of membrane / ρw
@show Tᵨ = 0.1/4*g*Lm*Lm #T/ρw
@show τ = 0.25#damping coeff
diriFlag = false

# Wave parameters
ω = 2.4#3.45#2.0#2.4
η₀ = 0.10
k = dispersionRelAng(H0, ω)
λ = 2*π/k
T = 2π/ω
ph0 = π/2
ηᵢₙ(x,t) = η₀*cos(k*x[1]-ω*t + ph0)
ϕᵢₙ(x,t) = (η₀*ω/k)*(cosh(k*(H0 + x[2])) / sinh(k*H0))*sin(k*x[1]-ω*t + ph0)
vᵢₙ(x,t) = -(η₀*ω)*(cosh(k*(H0 + x[2])) / sinh(k*H0))*cos(k*x[1]-ω*t + ph0)
vzᵢₙ(x,t) = ω*η₀*sin(k*x[1]-ω*t + ph0)
ηᵢₙ(t::Real) = x -> ηᵢₙ(x,t)
ϕᵢₙ(t::Real) = x -> ϕᵢₙ(x,t)
vᵢₙ(t::Real) = x -> vᵢₙ(x,t)
vzᵢₙ(t::Real) = x -> vzᵢₙ(x,t)
@show λ
@show T
@show η₀
@show H0
@show H0/λ
println()


# Domain 
nx = 1000
ny = 20
mesh_ry = 1.2 #Ratio for Geometric progression of eleSize
Ld = Lm #damping zone length
LΩ = 2*Ld + 3*Lm
x₀ = -Ld
domain =  (x₀, x₀+LΩ, -H0, 0.0)
partition = (nx, ny)
xdᵢₙ = 0.0
xdₒₜ = x₀ + LΩ - Ld
xm₀ = xdᵢₙ + Lm
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


# Time stepping
γₜ = 0.5
βₜ = 0.25
t₀ = 0.0
Δt = T/40
outΔt = T/4
tf = 15*T
∂uₜ_∂u = γₜ/(βₜ*Δt)
∂uₜₜ_∂u = 1/(βₜ*Δt^2)


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
μ₂ᵢₙ(x) = μ₁ᵢₙ(x)*k
μ₂ₒᵤₜ(x) = μ₁ₒᵤₜ(x)*k
ηd(t) = x -> μ₂ᵢₙ(x)*ηᵢₙ(x,t)
∇ₙϕd(t) = x -> μ₁ᵢₙ(x)*vzᵢₙ(x,t)


# Mesh
function f_y(y, r, n, H0; dbgmsg = false)
  # Mesh along depth as a GP
  # Depth is 0 to -H0    
  if(r ≈ 1.0)
    return y  
  else
    a0 = H0 * (r-1) / (r^n - 1)    
    if(dbgmsg)
      ln = 0:n
      ly = -a0 / (r-1) * (r.^ln .- 1)         
      @show hcat( ly, [ 0; ly[1:end-1] - ly[2:end] ] )
    end
    
    if y ≈ 0
      return 0.0
    end
    j = abs(y) / H0 * n  
    return -a0 / (r-1) * (r^j - 1)
  end
end
map(x) = VectorValue( x[1], f_y(x[2], mesh_ry, ny, H0; dbgmsg=false) )
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
function is_damping2(xs) # Check if an element is inside the damping zone 2
  n = length(xs)
  x = (1/n)*sum(xs)
  (xdₒₜ <= x[1] ) * ( x[2] ≈ 0.0)
end

# Masking and Beam Triangulation
xΓ = get_cell_coordinates(Γ)
Γm_to_Γ_mask = lazy_map(is_mem, xΓ)
Γd1_to_Γ_mask = lazy_map(is_damping1, xΓ)
Γd2_to_Γ_mask = lazy_map(is_damping2, xΓ)
Γm = Triangulation(Γ, findall(Γm_to_Γ_mask))
Γd1 = Triangulation(Γ, findall(Γd1_to_Γ_mask))
Γd2 = Triangulation(Γ, findall(Γd2_to_Γ_mask))
Γfs = Triangulation(Γ, findall(!, Γm_to_Γ_mask .| 
  Γd1_to_Γ_mask .| Γd2_to_Γ_mask))
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
  writevtk(Γd2,filename*"_Gd2")
  writevtk(Γfs,filename*"_Gfs")
  writevtk(Λmb,filename*"_Lmb")  
end


# Measures
degree = 2*order
dΩ = Measure(Ω,degree)
dΓm = Measure(Γm,degree)
dΓd1 = Measure(Γd1,degree)
dΓd2 = Measure(Γd2,degree)
dΓfs = Measure(Γfs,degree)
dΓin = Measure(Γin,degree)
dΛmb = Measure(Λmb,degree)


# Normals
@show nΛmb = get_normal_vector(Λmb)


# Dirichlet Fnc
gη(x,t) = 0.0
gη(t) = x -> gη(x,t)

# FE spaces
reffe = ReferenceFE(lagrangian,Float64,order)
V_Ω = TestFESpace(Ω, reffe, conformity=:H1)
V_Γκ = TestFESpace(Γκ, reffe, conformity=:H1)
if(diriFlag)
  V_Γη = TestFESpace(Γη, reffe, conformity=:H1,
    dirichlet_tags=["mem_bnd"])
else
  V_Γη = TestFESpace(Γη, reffe, conformity=:H1)
end
U_Ω = TransientTrialFESpace(V_Ω)
U_Γκ = TransientTrialFESpace(V_Γκ)
if(diriFlag)
  U_Γη = TransientTrialFESpace(V_Γη, gη)
else
  U_Γη = TransientTrialFESpace(V_Γη)
end
X = TransientMultiFieldFESpace([U_Ω,U_Γκ,U_Γη])
Y = MultiFieldFESpace([V_Ω,V_Γκ,V_Γη])


# Weak form
∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
if(diriFlag)
  m((ϕₜₜ,κₜₜ,ηₜₜ),(w,u,v)) = ∫( mᵨ*v*ηₜₜ )dΓm
  c((ϕₜ,κₜ,ηₜ),(w,u,v)) = 
    ∫(  βₕ*(u + αₕ*w)*ϕₜ - w*κₜ )dΓfs +
    ∫(  βₕ*(u + αₕ*w)*ϕₜ - w*κₜ )dΓd1    +
    ∫(  βₕ*(u + αₕ*w)*ϕₜ - w*κₜ )dΓd2    +
    ∫(  v*ϕₜ - w*ηₜ + Tᵨ*τ*∇(v)⋅∇(ηₜ) )dΓm +
    ∫(- Tᵨ*τ*v*∇(ηₜ)⋅nΛmb )dΛmb
  a((ϕ,κ,η),(w,u,v)) =      
    ∫(  ∇(w)⋅∇(ϕ) )dΩ   +
    ∫(  βₕ*(u + αₕ*w)*g*κ )dΓfs   +
    ∫(  βₕ*(u + αₕ*w)*g*κ - μ₂ᵢₙ*κ*w + μ₁ᵢₙ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd1  +
    ∫(  βₕ*(u + αₕ*w)*g*κ - μ₂ₒᵤₜ*κ*w + μ₁ₒᵤₜ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd2  +
    ∫(  v*(g*η) + Tᵨ*∇(v)⋅∇(η) )dΓm + 
    ∫(- Tᵨ*v*∇(η)⋅nΛmb )dΛmb

else
  m((ϕₜₜ,κₜₜ,ηₜₜ),(w,u,v)) = ∫( mᵨ*v*ηₜₜ )dΓm
  c((ϕₜ,κₜ,ηₜ),(w,u,v)) = 
    ∫(  βₕ*(u + αₕ*w)*ϕₜ - w*κₜ )dΓfs +
    ∫(  βₕ*(u + αₕ*w)*ϕₜ - w*κₜ )dΓd1    +
    ∫(  βₕ*(u + αₕ*w)*ϕₜ - w*κₜ )dΓd2    +
    ∫(  v*ϕₜ - w*ηₜ + Tᵨ*τ*∇(v)⋅∇(ηₜ) )dΓm #+
    # ∫(- Tᵨ*τ*v*∇(ηₜ)⋅nΛmb )dΛmb
  a((ϕ,κ,η),(w,u,v)) =      
    ∫(  ∇(w)⋅∇(ϕ) )dΩ   +
    ∫(  βₕ*(u + αₕ*w)*g*κ )dΓfs   +
    ∫(  βₕ*(u + αₕ*w)*g*κ - μ₂ᵢₙ*κ*w + μ₁ᵢₙ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd1  +
    ∫(  βₕ*(u + αₕ*w)*g*κ - μ₂ₒᵤₜ*κ*w + μ₁ₒᵤₜ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd2  +
    ∫(  v*(g*η) + Tᵨ*∇(v)⋅∇(η) )dΓm #+ 
    # ∫(- Tᵨ*v*∇(η)⋅nΛmb )dΛmb
end

l(t,(w,u,v)) =  
  ∫( w*vᵢₙ(t) )dΓin - 
  ∫( ηd(t)*w - ∇ₙϕd(t)*(u + αₕ*w) )dΓd1


# Solution
op = TransientConstantMatrixFEOperator(m,c,a,l,X,Y)
ls = LUSolver()
ode_solver = Newmark(ls,Δt,γₜ,βₜ)

# Initial solution
u0 = interpolate_everywhere([0.0,0.0,0.0],X(0.0))
u0t = interpolate_everywhere([0.0,0.0,0.0],X(0.0))
u0tt = interpolate_everywhere([0.0,0.0,0.0],X(0.0))

uht = solve(ode_solver,op,(u0,u0t,u0tt),t₀,tf)

if vtk_output == true
  pvd_Ω = paraview_collection(filename * "_O_sol", append=false)
  pvd_Γκ = paraview_collection(filename * "_Gk_sol", append=false)
  pvd_Γη = paraview_collection(filename * "_Ge_sol", append=false)
end

if vtk_output == true
  tpr = @sprintf("%5.3f",t₀)                    
  tval = @sprintf("%d",round(Int64,t₀*1000))
  ϕₕ, κₕ, ηₕ = u0
  pvd_Ω[t₀] = createvtk(Ω,
    filename * "_O_sol" * "_$tval.vtu",
    cellfields = ["phi" => ϕₕ])
  pvd_Γκ[t₀] = createvtk(Γκ,
    filename * "_Gk_sol" * "_$tval.vtu",
    cellfields = ["kappa" => κₕ])
  pvd_Γη[t₀] = createvtk(Γη,
    filename * "_Ge_sol" * "_$tval.vtu",
    cellfields = ["eta" => ηₕ])
end    


# Execute
@show outMod = round(Int64,outΔt/Δt);

for (uh, t) in uht  
    ϕₕ, κₕ, ηₕ = uh
    tpr = @sprintf("%5.3f",t)                    
    tval = @sprintf("%d",round(Int64,t*1000))

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
      pvd_Γη[t] = createvtk(Γη,
        filename * "_Ge_sol" * "_$tval.vtu",
        cellfields = ["eta" => ηₕ])
    end    
end

if vtk_output == true
  vtk_save(pvd_Ω)
  vtk_save(pvd_Γκ)
  vtk_save(pvd_Γη)
end

end