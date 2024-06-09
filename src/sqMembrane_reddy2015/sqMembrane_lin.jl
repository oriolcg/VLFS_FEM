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
axx = 1
ayy = 1
f0 = 1 #constant pressure
probpath = "data/sims_202212/sqMem_axx"*
  "$axx"*"_ayy"*"$ayy"*"_mesh1"
probname = probpath*"/sqMem"
if (isdir(probpath))
  rm(probpath, recursive=true)
end
mkdir(probpath)

simΔt = 0.10
outΔt = 0.2
order = 1

probesxy = [Point(0.0, 0.0)
            Point(0.5, 0.0)
            Point(0.0, 0.5)
            Point(-0.5, 0.0)
            Point(0.0, -0.5)];

# Generate Cartesian Domain 2DH 
domX = (-1, 1)
domY = (-1, 1)
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
  conformity=:H1, dirichlet_tags=["boundary"])


# Dirichlet BC
g1(x, t::Real) = 0.0
g1(t::Real) = x -> g1(x,t)


# Define Trial Space
U = TransientTrialFESpace(Ψ, g1)


# Define integration space
Ω = Triangulation(model)
dΩ = Measure(Ω, 2*order)
Γ = Boundary(model, tags = "boundary")
dΓ = Measure(Γ, 2*order)

nΓ = get_normal_vector(Γ)

@show a = VectorValue(axx, ayy)

# Weak form
resid(t, u, ψ) = 
  ∫( ψ * ∂tt(u) )dΩ +
  ∫( ∇(ψ)⋅ (∇(u) .* a) )dΩ -
  ∫( ψ * f0 )dΩ -
  ∫( ψ * (a .* ∇(u)) ⋅nΓ )dΓ


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

plot(probeDa[:,1], probeDa[:,4])
savefig(probname*"_plot.png")

end