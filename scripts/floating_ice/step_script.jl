function run_tmp_Step()
  case = Step_params(
    k=0.2,
    T=0.5, 
    Lb = 2*62.84,
    Ld = 62.84,
    xdₒᵤₜ = 3*62.84,
    name="Step",mesh_file="floating_ice-step_ratio05.json")
  x,η = run_Step(case)

  # Define execution function
  function run_tmp_step(case::Step_params)
    case_name = savename(case)
    println("-------------")
    println("Case: ",case_name)
    x,η = run_Step(case)
    case_name_suffix = savename(case,"jld2";digits=8)
    file = datadir("floating_ice", case_name_suffix)
    file = datadir("floating_ice", case_name_suffix)
    prefix,data,suffix = DrWatson.parse_savename(case_name_suffix,parsetypes=(Int, Float64))
    push!(data,"x"=>x, "η"=>η)
    save(file,data)
    return data
  end

  # # Case 1: depth_ratio=0.5  Q=0  k=0.4 
  # path = datadir("floating_ice")
  # case = Step_params(
  #   k=0.4,
  #   Lb = 70.0,
  #   Ld_Lb = 2,
  #   xdₒᵤₜ_Lb= 3,
  #   name="Step-mod-05-Q0-k04",mesh_file="floating_ice-step_ratio05.json")
  # @show case
  # data, file = produce_or_load(path,case,run_tmp_step)










  # # Case 2: depth_ratio=0.5  Q=0  k=0.4 
  # case = Step_params(
  #   k=0.4,
  #   Lb = 70.0,
  #   Ld_Lb = 2,
  #   xdₒᵤₜ_Lb= 3,
  #   name="Step05-Q0-k04",mesh_file="floating_ice_modified_step50.json")
  # @show case
  # data, file = produce_or_load(path,case,run_tmp_step)


  # # Case 2: ω=0.8
  # case = Liu_params(ω=0.8,name="omega-08")
  # @show case
  # data, file = produce_or_load(path,case,run_5_3_1)

  # # Gather data
  # res = collect_results(path)]


  # # Reference data
  # Liu_data_04 = CSV.File(datadir("Ref_data/Liu","omega_04.csv");header=false)
  # Liu_data_08 = CSV.File(datadir("Ref_data/Liu","omega_08.csv");header=false)

  # # Plot case 1
  # res1 = @linq res |> where(:ω .== 0.4)
  # xs1 = res1[!,:x][1]
  # η_xs1 = res1[!,:η][1]
  # p = sortperm(xs1)
  # plt1 = plot(xs1[p],η_xs1[p],
  #             xlims=(0,1),
  #             ylims=(0.6,1.4),
  #             lw=2,
  #             label="Monolithic CG/DG",
  #             palette=:rainbow)
  # plot!(plt1,Liu_data_04.Column1,Liu_data_04.Column2,marker="o",line=false,label="Liu et al.")
  # xlabel!("x/L")
  # ylabel!("|η|/η₀")
  # savefig(plt1, plotsdir("5-3-1-Liu","omega_04"))

  # # Plot case 2
  # res1 = @linq res |> where(:ω .== 0.8)
  # xs1 = res1[!,:x][1]
  # η_xs1 = res1[!,:η][1]
  # p = sortperm(xs1)
  # plt1 = plot(xs1[p],η_xs1[p],
  #             xlims=(0,1),
  #             ylims=(0,1.2),
  #             lw=2,
  #             label="Monolithic CG/DG",
  #             palette=:rainbow)
  # plot!(plt1,Liu_data_08.Column1,Liu_data_08.Column2,marker="o",line=false,label="Liu et al.")
  # xlabel!("x/L")
  # ylabel!("|η|/η₀")
  # savefig(plt1, plotsdir("5-3-1-Liu","omega_08"))

end
