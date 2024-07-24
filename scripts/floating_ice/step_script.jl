function run_tmp_Step()

  case = Step_params(ω=0.2,name="Step")
  x,η = run_Step(case)

  # Define execution function
  function run_tmp_step(case::Step_params)
    case_name = savename(case)
    println("-------------")
    println("Case: ",case_name)
    x,η = run_Step(case)
    case_name_suffix = savename(case,"jld2";digits=8)
    file = datadir("floating_ice_output", case_name_suffix)
    file = datadir("floating_ice_output", case_name_suffix)
    prefix,data,suffix = DrWatson.parse_savename(case_name_suffix,parsetypes=(Int, Float64))
    push!(data,"x"=>x, "η"=>η)
    save(file,data)
    return data
  end



  # # Case 1: ω=0.4
  # path = datadir("5-3-1-Liu")
  # case = Liu_params(ω=0.4,name="omega-04")
  # @show case
  # data, file = produce_or_load(path,case,run_5_3_1)

  # # Case 2: ω=0.8
  # case = Liu_params(ω=0.8,name="omega-08")
  # @show case
  # data, file = produce_or_load(path,case,run_5_3_1)

  # # Gather data
  # res = collect_results(path)

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
