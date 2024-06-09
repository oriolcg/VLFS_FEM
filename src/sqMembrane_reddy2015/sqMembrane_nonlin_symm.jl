# Reference
# Reddy (2015) : Reddy, J.N., 2015. An introduction to nonlinear 
# finite element analysis: with applications to heat transfer, 
# fluid mechanics, and solid mechanics, Second edition. ed. 
# Oxford University Press, Oxford.
#
# Example 6.7.3
module sqMembrane

using Revise
using Gridap
using LineSearches: BackTracking
using Printf
using Plots
using DataFrames
using CSV

# Constants
f0 = 1 #constant pressure
a1ux = 0.2
a1uy = 0.2
simΔt = 0.100
outΔt = 0.2
order = 1
probpath = "data/sims_202212/sqMem_nonlin_symm_a1ux"*
  "$a1ux"*"_a1uy"*"$a1uy"*"_mesh1"
probname = probpath*"/sqMem"
if (isdir(probpath))
  rm(probpath, recursive=true)
end
mkdir(probpath)

probesxy = [Point(0.0, 0.0)
            Point(0.5, 0.0)
            Point(0.0, 0.5)];

# Generate Cartesian Domain 2DH 
domX = (0, 1)
domY = (0, 1)
dx = 0.125
dy = dx
domain = (domX[1], domX[2], domY[1], domY[2])
partition = ( Int((domX[2]-domX[1])/dx), 
  Int((domY[2]-domY[1])/dy))
model = CartesianDiscreteModel(domain, partition)
writevtk(model, probname*"_model")


# Label sides of the domain
labels = get_face_labeling(model)


# Define Test Fnc
reffe = ReferenceFE(lagrangian, Float64, order)
Ψ = TestFESpace(model, reffe, 
  conformity=:H1, dirichlet_tags=[2,3,4,6,8])


# Dirichlet BC
g1(x, t::Real) = 0.0
g1(t::Real) = x -> g1(x,t)


# Define Trial Space
U = TransientTrialFESpace(Ψ, g1)


# Define integration space
Ω = Triangulation(model)
dΩ = Measure(Ω, 2*order)
Γ = Boundary(model, tags = [2,3,4,6,8])
dΓ = Measure(Γ, 2*order)

nΓ = get_normal_vector(Γ)

@show a1u = VectorValue(a1ux, a1uy)
a(du) = (1+du⋅a1u)*VectorValue(1.0,0.0) + 
  (1+du⋅a1u)*VectorValue(0.0,1.0)
aSc(du) = 1+du⋅a1u  


# Weak form
resid(t, u, ψ) = 
  ∫( ψ * ∂tt(u) )dΩ +
  ∫( ∇(ψ)⋅ (∇(u) .* a(∇(u))) )dΩ -
  ∫( ψ * f0 )dΩ -
  ∫( ψ * (a(∇(u)) .* ∇(u)) ⋅ nΓ )dΓ

# resid(t, u, ψ) = 
#   ∫( ψ * ∂tt(u) )dΩ +
#   ∫( ∇(ψ)⋅ (∇(u) * aSc(∇(u))) )dΩ -
#   ∫( ψ * f0 )dΩ -
#   ∫( ψ * (aSc(∇(u)) * ∇(u)) ⋅ nΓ )dΓ


op_AD = TransientFEOperator(resid, U, Ψ; order=2)


# Initial Soln
t0 = 0.0
tF = 5.0
u0 = interpolate_everywhere(0.0, U(t0))
u0t = interpolate_everywhere(0.0, U(t0))
u0tt = interpolate_everywhere(0.0, U(t0))


#NL Solver
nls = NLSolver(show_trace=true, method=:newton, 
  linesearch=BackTracking(), iterations=10)
ode_solver = Newmark(nls, simΔt, 0.5, 0.25)        
#ode_solver = GeneralizedAlpha(nls, simΔt, 0.0)    
solnht = solve(ode_solver, op_AD, (u0, u0t, u0tt), t0, tF)

# for (solh, tn) in solnht        
#   println(tn)        
# end

# Probe setup
numP = length(probesxy)   
lDa = zeros(Float64, 1, numP*3+1)    
probeDa = zeros(Float64, 1, numP*3+1)    
probeDa = DataFrame(probeDa, :auto)
probeDa[1,1] = t0;    

createpvd(probname) do pvd
  uh = u0
  tval = @sprintf("%5.3f",t0)                
  println("Time : $tval")
  tval = @sprintf("%d",floor(Int64,t0*1000))                
  pvd[t0] = createvtk(Ω,probname*"_$tval"*".vtu",
      cellfields=["u"=>uh])
end


# Execute
outMod = floor(Int64,outΔt/simΔt);
createpvd(probname, append=true) do pvd    
  cnt=0
  for (solh, t) in solnht                       
    cnt = cnt+1
    uh = solh
    tval = @sprintf("%5.3f",t)                
    println("Time : $tval")
    tval = @sprintf("%d",floor(Int64,t*1000))                

    lDa[1] = t  
    for ip in 1:numP
        ipxy = probesxy[ip];
        ik = (ip-1)*3 + 1
        lDa[ik+1] = ipxy[1]
        lDa[ik+2] = ipxy[2]
        lDa[ik+3] = uh(ipxy)        
    end
    push!(probeDa, lDa)

    if(cnt%outMod != 0) 
        continue
    end
    pvd[t] = createvtk(Ω,probname*"_$tval"*".vtu",
      cellfields=["u"=>uh])
  end
end    

@show probeDa

CSV.write(probname*"_probes.csv", probeDa)

probeDa0 = DataFrame(CSV.File("data/sims_202212/reddy2015_ref6.7.3.csv"))
probeDa_lin = DataFrame(CSV.File("data/sims_202212/sqMem_nonlin_a1ux0.0_a1uy0.0_mesh1_dt100/sqMem_probes.csv"))

labels=["data"
        "lin-dt100"
        "nonlin-symm-dt100"]

plot(probeDa0[:,1], probeDa0[:,2].*1.5.-0.5, 
  seriestype=:scatter, label=labels[1])
plot!(probeDa_lin[:,1], probeDa_lin[:,4], linewidth=3, 
  label=labels[2])
plot!(probeDa[:,1], probeDa[:,4], linewidth=3, label=labels[3], 
  size=(1280,640), legend=:bottomright)
savefig(probname*"_plot.png")

end