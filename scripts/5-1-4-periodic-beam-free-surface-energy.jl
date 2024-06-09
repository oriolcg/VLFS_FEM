function run_5_1_4_periodic_beam_free_surface_energy()

  # Define Execution function
  function run_5_1_4(case::Periodic_Beam_FS_params)
    case_name = savename(case)
    println("-------------")
    println("Case: ",case_name)
    (e_ϕ, e_η, E_kin_f , E_pot_f, E_kin_s, E_ela_s, E_kin_f₀, E_kin_s₀, E_pot_f₀, E_ela_s₀, time) = run_periodic_beam_FS(case)
    e_ϕ_i = last(e_ϕ)
    e_η_i = last(e_η)
    case_name_suffix = savename(case,"jld2";digits=8)
    file = datadir("5-1-4-periodic-beam-free-surface-energy", case_name_suffix)
    prefix, data, suffix = DrWatson.parse_savename(case_name_suffix, parsetypes=(Int, Float64))
    push!(data,
    "E_kin_f"=>E_kin_f,
    "E_pot_f"=>E_pot_f,
    "E_kin_s"=>E_kin_s,
    "E_ela_s"=>E_ela_s,
    "E_kin_f₀"=>E_kin_f₀,
    "E_kin_s₀"=>E_kin_s₀,
    "E_pot_f₀"=>E_pot_f₀,
    "E_ela_s₀"=>E_ela_s₀,
    "time"=>time
    )
    #@tagsave(file,data)
    return data
  end

  # Warm-up case
  k = 1
  H = 1.0
  g = 9.81
  ω = √(g*k*tanh(k*H))
  T = 2π/ω
  Δt = T/1
  n = 3
  order = 2
  path = datadir("5-1-4-periodic-beam-free-surface-energy")
  case = Periodic_Beam_FS_params(
    name="11Warm-up",
    n=n,
    dt=Δt,
    tf=T,
    k=k,
    order=order,
    vtk_output=true
  )
  filename = config -> savename(config; digits=8)
  produce_or_load(path,case,run_5_1_4; filename=filename)

  # parameters
  k  = 15
  ω = √(9.81*k*tanh(k*1.0))
  T = 2π/ω
  tf1 = 10*T

  # Element size Convergence
  Δt = 1.0e-3
  order = 4
  for i in 1:5
    nelem = 2^(i+1)
    case = Periodic_Beam_FS_params(
      name="EnergyEvolution",
      n=nelem,
      dt=Δt,
      tf=tf1,
      k=k,
      order=order
    )
    filename = config -> savename(config; digits=8)
    data, file = produce_or_load(path,case,run_5_1_4; filename=filename)
  end

  # Time step size Convergence
  order = 4
  nelem = 64
  tf2 = 1*T
  for i in 1:5
    Δt = T * 2^(-1.0-i)
    case = Periodic_Beam_FS_params(
      name="EnergyEvolution",
      n=nelem,
      dt=Δt,
      tf=tf2,
      k=k,
      order=order
    )
    filename = config -> savename(config; digits=8)
    data, file = produce_or_load(path,case,run_5_1_4; filename=filename)
  end

  res = collect_results(path)

  println("Ploting Energy evolution")
  plt1 = plot(
    fontsize=12,
    legend=:outerright,
    legendfontsize=10,
    xtickfontsize=10,
    ytickfontsize=10,
  )
  plt2 = plot(
    fontsize=12,
    legend=:topleft,
    legendfontsize=10,
    xtickfontsize=10,
    ytickfontsize=10,
    xaxis=:log,
    yaxis=:log,
  )
  xlabel!(plt1,"Time")
  xlabel!(plt2,"Time step size")
  ylabel!(plt1,"Relative energy error")
  ylabel!(plt2,"Relative energy error")
  styles = [:dash,:dashdot,:dashdotdot,:dot]
  shapes = [:square,:circle,:utriangle,:diamond]
  colors = cgrad(:winter,4,categorical=true)
  factors_plt2 = [10,30,80,80]

  # Plot energy evolution vs elements
  colors = cgrad(:winter,5,categorical=true)
  for i in 2:5
    order = 4
    Δti = 1.0e-3
    nelem = 2^(i+1)
    nx = 2*nelem
    res_elem  = @linq res |> where(:order .== order, :k .== k, :dt .== Δti, :n.==nelem, :tf .== round(tf1,digits=8))
    t = res_elem[!,:time]
    E_tot = res_elem[!,:E_kin_f]+res_elem[!,:E_pot_f]+res_elem[!,:E_kin_s]+res_elem[!,:E_ela_s]
    E_tot₀ = res_elem[!,:E_kin_f₀]+res_elem[!,:E_pot_f₀]+res_elem[!,:E_kin_s₀]+res_elem[!,:E_ela_s₀]
    E_tot = abs.(E_tot[1] .- E_tot₀) ./ E_tot₀
    plot!(plt1,
      t[1],E_tot,
      yaxis=:log,
      color=colors[i],
      label="nx=$nx"
    )
  end

  # Plot energy evolution vs time step size
  Eh = Float64[]
  Δts = Float64[]
  for i in 1:5
    order = 4
    nelem = 64
    Δt = round(T * 2^(-1.0-i),digits=8)
    res_time  = @linq res |> where(:order .== order, :k .== k, :dt .== Δt, :n.==nelem, :tf .== round(tf2,digits=8))
    t = res_time[!,:time]
    E_tot = res_time[!,:E_kin_f]+res_time[!,:E_pot_f]+res_time[!,:E_kin_s]+res_time[!,:E_ela_s]
    E_tot₀ = res_time[!,:E_kin_f₀]+res_time[!,:E_pot_f₀]+res_time[!,:E_kin_s₀]+res_time[!,:E_ela_s₀]
    E_tot = abs.(E_tot[1] .- E_tot₀) ./ E_tot₀
    push!(Eh,E_tot[end])
    push!(Δts,Δt)
  end
  plot!(plt2,
    Δts,Eh,
    shape=shapes[3],
    color=:red,
    style=styles[3],
    linewidth=2,
    label="r=4, n=128"
  )
  plot!(plt2,
    Δts,0.2*Δts.^(2),
    color=:black,
    style=styles[3],
    linewidth=2,
    label=latexstring("dt^{-2}"),
    xticks=(Δts,["T/$(2^(i+1))" for i in 1:length(Δts)])
  )

  # Save
  savefig(plt1,plotsdir("5-1-4-periodic-beam-free-surface-energy","error_mesh_vs_time.png"))
  savefig(plt2,plotsdir("5-1-4-periodic-beam-free-surface-energy","error_convergence_time.png"))

end
