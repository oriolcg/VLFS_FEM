module BeamMultJoints_freq

using Revise
using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Plots
using DrWatson

name::String = "data/sims_202212/BeamMultJoints-freq-method3"
order::Int = 2
ξ::Float64 = 0.0
vtk_output::Bool = true
filename = name*"/dom"

# Beam parameters
g = 9.81
ρw = 1025 #kg/m3 water
Lb = 10 #m
ρb = 7850 # density kg/m3
EI = 47100.0 
h0 = 0.001 #m
β = (0.2, 0.5, 0.7) # relative location of joints
# β = 0.5
nj = length(β) # num_of_joints
d₀ = ρb*h0/ρw #draft of the floating body
Dᵨ = EI/ρw
kᵣ = ξ*Dᵨ/Lb


# Wave parameters
α = 0.249 # tune wave-length
H0 = 1.1 #m #still-water depth, atleast at inlet
λ = α*Lb #m #wave-length
η₀ = 0.01 #m #wave-amplitude
k = 2π/λ
ω = sqrt(g*k*tanh(k*H0))
T = 2π/ω
ηᵢₙ(x) = η₀*exp(im*k*x[1])
ϕᵢₙ(x) = -im*(η₀*ω/k)*(cosh(k*(H0 + x[2])) / sinh(k*H0))*exp(im*k*x[1])
vxᵢₙ(x) = (η₀*ω)*(cosh(k*(H0 + x[2])) / sinh(k*H0))*exp(im*k*x[1])
vzfsᵢₙ(x) = -im*ω*η₀*exp(im*k*x[1]) #???
@show λ
@show T
@show H0
@show H0/λ


# Domain 
nx = 640
ny = 20
mesh_ry = 1.1 #Ratio for Geometric progression of eleSize
Ld = Lb #damping zone length
LΩ = 2*Ld + 2*Lb
x₀ = -Ld
domain =  (x₀, x₀+LΩ, -H0, 0.0)
partition = (nx, ny)
xdᵢₙ = 0.0
xdₒₜ = x₀ + LΩ - Ld
xb₀ = xdᵢₙ + Lb/2.0
xb₁ = xb₀ + Lb
xbⱼ = xb₀ .+ Lb.*β
@show Lb
@show LΩ
@show domain
@show partition
@show (xb₀, xb₁)
@show xbⱼ
@show Lb/LΩ*nx
@show LΩ/nx
@show H0/ny


# Numeric constants
nx_total = Int(ceil(nx/minimum(β))*ceil(LΩ/Lb))
h = LΩ / nx_total
γ = 1.0*order*(order-1)/h
βₕ = 0.5
αₕ = -im*ω/g * (1-βₕ)/βₕ
@show nx_total
@show h
@show βₕ
@show αₕ


# Damping
μ₀ = 2.5
μ₁ᵢₙ(x) = μ₀*(1.0 - sin(π/2*(x[1]-x₀)/Ld))
μ₁ₒᵤₜ(x) = μ₀*(1.0 - cos(π/2*(x[1]-xdₒₜ)/Ld))
μ₂ᵢₙ(x) = μ₁ᵢₙ(x)*k
μ₂ₒᵤₜ(x) = μ₁ₒᵤₜ(x)*k
ηd(x) = μ₂ᵢₙ(x)*ηᵢₙ(x)
∇ₙϕd(x) = μ₁ᵢₙ(x)*vzfsᵢₙ(x) #???


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
writevtk(model, filename*"_model")


# Triangulations
Ω = Interior(model) #same as Triangulation()
Γ = Boundary(model,tags="surface") #same as BoundaryTriangulation()
Γin = Boundary(model,tags="inlet")


# Auxiliar functions
function is_beam(xs) # Check if an element is inside the beam1
  n = length(xs)
  x = (1/n)*sum(xs)
  (xb₀ <= x[1] <= xb₁ ) * ( x[2] ≈ 0.0)
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
# function is_a_joint(xs) # Check if an element is a joint  
#   is_on_xbⱼ = [x[1]≈xbⱼ && x[2]≈0.0 for x in xs] # array of booleans of size the number of points in an element (for points, it will be an array of size 1)    
#   element_on_xbⱼ = minimum(is_on_xbⱼ) # Boolean with "true" if at least all entries are true, "false" otherwise.
#   element_on_xbⱼ
# end
function is_a_joint(xs) # Check if an element is a joint  
  # @show xs
  is_on_xbⱼ = [x[1]≈y && x[2]≈0.0 for x in xs for y in xbⱼ] 
  element_on_xbⱼ = maximum(is_on_xbⱼ) # Boolean with "true" if at least all entries are true, "false" otherwise.
  element_on_xbⱼ
end

function is_a_beam_node(xs) # Check if an element is a joint  
  #@show xs
  is_on_xb = [(x[2] ≈ 0.0) && (xb₀ < x[1] < xb₁) for x in xs ]
  return maximum(is_on_xb)
end


# Masking and Beam Triangulation
xΓ = get_cell_coordinates(Γ)
Γb_to_Γ_mask = lazy_map(is_beam, xΓ)
Γd1_to_Γ_mask = lazy_map(is_damping1, xΓ)
Γd2_to_Γ_mask = lazy_map(is_damping2, xΓ)
Γb = Triangulation(Γ, findall(Γb_to_Γ_mask))
Γd1 = Triangulation(Γ, findall(Γd1_to_Γ_mask))
Γd2 = Triangulation(Γ, findall(Γd2_to_Γ_mask))
Γfs = Triangulation(Γ, findall(!, Γb_to_Γ_mask .| 
  Γd1_to_Γ_mask .| Γd2_to_Γ_mask))
Γη = Triangulation(Γ, findall(Γb_to_Γ_mask))
Γκ = Triangulation(Γ, findall(!,Γb_to_Γ_mask))


# Construct the mask for the joint (Oriol ref)
Γ_mask_in_Ω_dim_0 = get_face_mask(labels_Ω,"surface",0)
grid_dim_0_Γ = GridPortion(Grid(ReferenceFE{0},model),Γ_mask_in_Ω_dim_0)
xΓ_dim_0 = get_cell_coordinates(grid_dim_0_Γ)
Λj_to_Γ_mask = lazy_map(is_a_joint,xΓ_dim_0)
ΛjOrl = Skeleton(Γ,Λj_to_Γ_mask)

# Construct the mask for the joint
#Note: Skeleton does not include the boundary of the selection
xΓ_skel = get_cell_coordinates(Skeleton(Γ))
Λj_to_Γ_mask = lazy_map(is_a_joint, xΓ_skel)
@show Λj = Triangulation(Skeleton(Γ), findall(Λj_to_Γ_mask))
@show Λj == ΛjOrl
@show get_cell_coordinates(Λj) == get_cell_coordinates(ΛjOrl)
# @show get_cell_coordinates(Λj) 
# @show get_cell_coordinates(ΛjOrl)
# @show get_normal_vector(ΛjOrl).⁺
# @show get_normal_vector(Λj).⁺


# Construct the mask for beam nodes
@show Λb = Skeleton(Γb)
xΓb_skel = get_cell_coordinates(Skeleton(Γ))
Λb_to_Γ_mask = lazy_map(is_a_beam_node, xΓb_skel)
@show Λj_to_Γ_mask
@show Λb_to_Γ_mask
@show Λb = Triangulation(Skeleton(Γ), 
  findall(Λb_to_Γ_mask .&& .!(Λj_to_Γ_mask)))

if vtk_output == true
  writevtk(Ω,filename*"_O")
  writevtk(Γ,filename*"_G")
  writevtk(Γb,filename*"_Gb")  
  writevtk(Γd1,filename*"_Gd1")
  writevtk(Γd1,filename*"_Gd2")
  writevtk(Γfs,filename*"_Gfs")
  writevtk(Λb,filename*"_Lb")  
  writevtk(Λj,filename*"_Lj")
end



# Measures
degree = 2*order
dΩ = Measure(Ω,degree)
dΓb = Measure(Γb,degree)
dΓd1 = Measure(Γd1,degree)
dΓd2 = Measure(Γd2,degree)
dΓfs = Measure(Γfs,degree)
dΓin = Measure(Γin,degree)
dΛb = Measure(Λb,degree)
dΛj = Measure(Λj,degree)


# Normals
nΛb = get_normal_vector(Λb)
nΛj = get_normal_vector(Λj)


# FE spaces
reffe = ReferenceFE(lagrangian,Float64,order)
V_Ω = TestFESpace(Ω, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
V_Γκ = TestFESpace(Γκ, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
V_Γη = TestFESpace(Γη, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
U_Ω = TrialFESpace(V_Ω)
U_Γκ = TrialFESpace(V_Γκ)
U_Γη = TrialFESpace(V_Γη)
X = MultiFieldFESpace([U_Ω,U_Γκ,U_Γη])
Y = MultiFieldFESpace([V_Ω,V_Γκ,V_Γη])


# Weak form
∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
a((ϕ,κ,η),(w,u,v)) =      ∫(  ∇(w)⋅∇(ϕ) )dΩ   +
∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ )dΓfs   +
∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ 
  - μ₂ᵢₙ*κ*w + μ₁ᵢₙ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd1    +
∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ 
  - μ₂ₒᵤₜ*κ*w + μ₁ₒᵤₜ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd2    +
∫(  ( v*((-ω^2*d₀ + g)*η - im*ω*ϕ) 
  + Dᵨ*Δ(v)*Δ(η) ) +  im*ω*w*η  )dΓb  +
∫(  Dᵨ * ( - jump(∇(v)⋅nΛb) * mean(Δ(η)) - mean(Δ(v)) 
  * jump(∇(η)⋅nΛb) + γ*( jump(∇(v)⋅nΛb) * jump(∇(η)⋅nΛb) ) ) )dΛb +
∫(  (jump(∇(v)⋅nΛj) * kᵣ * jump(∇(η)⋅nΛj)) )dΛj

l((w,u,v)) =  ∫( w*vxᵢₙ )dΓin - ∫( ηd*w - ∇ₙϕd*(u + αₕ*w) )dΓd1


# Solution
op = AffineFEOperator(a,l,X,Y)
(ϕₕ,κₕ,ηₕ) = solve(op)

if vtk_output == true
  writevtk(Ω,filename * "_O_sol.vtu",
    cellfields = ["phi_re" => real(ϕₕ),"phi_im" => imag(ϕₕ),
    "phi_abs" => abs(ϕₕ), "phi_ang" => angle∘(ϕₕ)])
  writevtk(Γκ,filename * "_Gk_sol.vtu",
    cellfields = ["eta_re" => real(κₕ),"eta_im" => imag(κₕ),
    "eta_abs" => abs(κₕ), "eta_ang" => angle∘(κₕ)])
  writevtk(Γη,filename * "_Ge_sol.vtu",
    cellfields = ["eta_re" => real(ηₕ),"eta_im" => imag(ηₕ),
    "eta_abs" => abs(ηₕ), "eta_ang" => angle∘(ηₕ)])
end

data = Dict("ϕₕ" => ϕₕ,
            "κₕ" => κₕ,
            "ηₕ" => ηₕ)

wsave(filename*"_data.jld2", data)
end