module BeamMultJoints_freq

using Revise
using Gridap
using Plots
using DrWatson
using WaveSpec
using .Constants
using .WaveFreqDom


name::String = "data/sims_202303/run/empt_freq"
order::Int = 2
vtk_output::Bool = true
filename = name*"/mem"

H0 = 10 #m #still-water depth


# Wave parameters
@show g #defined in .Constants
ω = 8
k = dispersionRelAng(H0, ω)
λ = 2π/k #21 #13.95 #0.5*Lm #m #wave-length
η₀ = 0.25#ω/g #m #wave-amplitude
T = 2π/ω
ηᵢₙ(x) = η₀*exp(im*k*x[1])
ϕᵢₙ(x) = -im*(η₀*ω/k)*(cosh(k*(H0 + x[2])) / sinh(k*H0))*exp(im*k*x[1])
vxᵢₙ(x) = (η₀*ω)*(cosh(k*(H0 + x[2])) / sinh(k*H0))*exp(im*k*x[1])
vzfsᵢₙ(x) = -im*ω*η₀*exp(im*k*x[1]) #???
@show λ
@show T
@show η₀
@show H0
@show H0/λ
println()


# Domain 
nx = 2600
ny = 20
mesh_ry = 1.1 #Ratio for Geometric progression of eleSize
Ld = 20*H0 #damping zone length
LΩ = Ld + 3*2*H0
x₀ = -Ld
domain =  (x₀, x₀+LΩ, -H0, 0.0)
partition = (nx, ny)
xdᵢₙ = 0.0
@show LΩ
@show domain
@show partition
@show LΩ/nx
@show H0/ny
println()


# Numeric constants
h = LΩ / nx
γ = 1.0*order*(order-1)/h
βₕ = 0.5
αₕ = -im*ω/g * (1-βₕ)/βₕ
@show h
@show βₕ
@show αₕ
println()


# Damping
#μ₀ = 1.5
μ₀ = maximum([1.5, 5.4/(ω^0.760)]) #5.61/ω#6/(ω^0.6)
if(ω>3.5)then
  μ₁ᵢₙ(x) = ifelse( (x[1]-x₀)/(2*λ) < 1.0, 
    μ₀*(1.0 - sin(π/2*(x[1]-x₀)/(2*λ))), 0.0 )
else 
  μ₁ᵢₙ(x) = μ₀*(1.0 - sin(π/2*(x[1]-x₀)/Ld))
end
μ₂ᵢₙ(x) = μ₁ᵢₙ(x)*k
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


# Triangulations
Ω = Interior(model) #same as Triangulation()
Γ = Boundary(model,tags="surface") #same as BoundaryTriangulation()
Γin = Boundary(model,tags="inlet")
Γot = Boundary(model,tags="outlet")


# Auxiliar functions
function is_damping1(xs) # Check if an element is inside the damping zone 1
  n = length(xs)
  x = (1/n)*sum(xs)
  (x₀ <= x[1] <= xdᵢₙ ) * ( x[2] ≈ 0.0)
end


# Masking damping zones
xΓ = get_cell_coordinates(Γ)
Γd1_to_Γ_mask = lazy_map(is_damping1, xΓ)
Γd1 = Triangulation(Γ, findall(Γd1_to_Γ_mask))
Γfs = Triangulation(Γ, findall(!, Γd1_to_Γ_mask ))
Γκ = Γ


writevtk(model, filename*"_model")
if vtk_output == true
  writevtk(Ω,filename*"_O")
  writevtk(Γ,filename*"_G")
  writevtk(Γd1,filename*"_Gd1")
  writevtk(Γfs,filename*"_Gfs")
end


# Measures
degree = 2*order
dΩ = Measure(Ω,degree)
dΓd1 = Measure(Γd1,degree)
dΓfs = Measure(Γfs,degree)
dΓin = Measure(Γin,degree)
dΓot = Measure(Γot,degree)


# FE spaces
reffe = ReferenceFE(lagrangian,Float64,order)
V_Ω = TestFESpace(Ω, reffe, conformity=:H1, 
  vector_type=Vector{ComplexF64})
V_Γκ = TestFESpace(Γκ, reffe, conformity=:H1, 
  vector_type=Vector{ComplexF64})
U_Ω = TrialFESpace(V_Ω)
U_Γκ = TrialFESpace(V_Γκ)
X = MultiFieldFESpace([U_Ω,U_Γκ])
Y = MultiFieldFESpace([V_Ω,V_Γκ])


# Weak form
∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
a((ϕ,κ),(w,u)) =      
  ∫(  ∇(w)⋅∇(ϕ) )dΩ   +
  ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ )dΓfs   +
  ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ 
    - μ₂ᵢₙ*κ*w + μ₁ᵢₙ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd1    +
  ∫( -w * im * k * ϕ )dΓot 
  # ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ 
  #   - μ₂ₒᵤₜ*κ*w + μ₁ₒᵤₜ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd2    

l((w,u)) =  ∫( w*vxᵢₙ )dΓin - ∫( ηd*w - ∇ₙϕd*(u + αₕ*w) )dΓd1


# Solution
op = AffineFEOperator(a,l,X,Y)
(ϕₕ,κₕ) = solve(op)

κₓ = ∇(κₕ) ⋅ VectorValue(1.0,0.0)

if vtk_output == true
  writevtk(Ω,filename * "_O_sol.vtu",
    cellfields = ["phi_re" => real(ϕₕ),"phi_im" => imag(ϕₕ),
    "phi_abs" => abs(ϕₕ), "phi_ang" => angle∘(ϕₕ)])
  writevtk(Γκ,filename * "_Gk_sol.vtu",
    cellfields = ["eta_re" => real(κₕ),"eta_im" => imag(κₕ),
    "eta_abs" => abs(κₕ), "eta_ang" => angle∘(κₕ)] )  
    writevtk(Γκ,filename * "_Gkx_sol.vtu",
    cellfields = ["eta_x_re" => real(κₓ),"eta_x_im" => imag(κₓ),
    "eta_x_abs" => abs(κₓ), "eta_x_ang" => angle∘(κₓ)])  
end

data = Dict("ϕₕ" => ϕₕ,
            "κₕ" => κₕ)

wsave(filename*"_data.jld2", data)

end