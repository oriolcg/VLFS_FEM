module Membrane_modes

using Revise
using Gridap
using Plots
using DrWatson
using WaveSpec
using .Constants
using LinearAlgebra
using TickTock

name::String = "data/sims_202302/mem_modes_dry_free"
order::Int = 1
vtk_output::Bool = true
filename = name*"/mem"

ρw = 1025 #kg/m3 water
H0 = 10 #m #still-water depth

# Membrane parameters
Lm = 2*H0 #m
@show g #defined in .Constants
mᵨ = 0.9 #mass per unit area of membrane / ρw
Tᵨ = 0.1*g*H0*H0 #T/ρw
@show τ = 0.5#damping coeff

# Excitation wave parameters
ω = 2.4
k = dispersionRelAng(H0, ω)


# Domain 
nx = 450
ny = 18
mesh_ry = 1.1 #Ratio for Geometric progression of eleSize
Ld = 1.5*H0 #damping zone length
LΩ = 6*H0 + 2*Ld
x₀ = -Ld
domain =  (x₀, x₀+LΩ, -H0, 0.0)
partition = (nx, ny)
xdᵢₙ = 0.0
xdₒₜ = x₀ + LΩ - Ld
xm₀ = xdᵢₙ + 2*H0
xm₁ = xm₀ + Lm
@show Lm
@show LΩ
@show domain
@show partition
@show (xm₀, xm₁)
@show isinteger(Lm/LΩ*nx)
@show LΩ/nx
@show H0/ny
@show Ld*k/2/π
println()


# Damping
μ₀ = 2.5
μ₁ᵢₙ(x) = μ₀*(1.0 - sin(π/2*(x[1]-x₀)/Ld))
μ₁ₒᵤₜ(x) = μ₀*(1.0 - cos(π/2*(x[1]-xdₒₜ)/Ld))
μ₂ᵢₙ(x) = μ₁ᵢₙ(x)*k
μ₂ₒᵤₜ(x) = μ₁ₒᵤₜ(x)*k
ηd(x) = μ₂ᵢₙ(x)*ηᵢₙ(x)
#∇ₙϕd(x) = μ₁ᵢₙ(x)*vzfsᵢₙ(x) #???



# Mesh
function f_y(y, r, n, H0; dbgmsg = false)
  # Mesh along depth as a GP
  # Depth is 0 to -H0    
  if(r ≈ 1.0)
    return y  
  else
    a0 = H0 * (r-1) / (r^n - 1)    
    if(dbgmsg)
      println("[VAL] Smallest dy = ", a0)
      println("[VAL] Largest dy = ", 
        -a0/(r-1)*(r^(n-1)-1) + a0/(r-1)*(r^n-1))
    end
    
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
# U_Γη = TrialFESpace(V_Γη, gη) #diri
U_Γη = TrialFESpace(V_Γη)


# Weak form
∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
m11(η,v) = ∫( mᵨ*v*η )dΓm
k11(η,v) = ∫( Tᵨ*∇(v)⋅∇(η) )dΓm #+  
            #∫(- Tᵨ*v*∇(η)⋅nΛmb )dΛmb #diri


l1(v) = ∫( 0*v )dΓm
println("[MSG] Done Weak form")

# Global matrices
M11 = get_matrix(AffineFEOperator( m11, l1, U_Γη, V_Γη ))
K11 = get_matrix(AffineFEOperator( k11, l1, U_Γη, V_Γη ))
println("[MSG] Done Global matrices")

λ = LinearAlgebra.eigvals(M11\Matrix(K11))
V = LinearAlgebra.eigvecs(M11\Matrix(K11))
@show sum(imag.(λ))
ωₙ = real.(λ)
ωₙ = ifelse.(ωₙ .<0, 0, ωₙ)
@show ωₙ = sqrt.(ωₙ)
xp = range(xm₀, xm₁, size(V,2))

data = Dict(
  "xp" => xp,
  "ωₙ" => ωₙ,
  "V" => V
)

wsave(filename*"_modesdata.jld2", data)

end