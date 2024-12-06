function run_tmp_Step(ω_in)
  case = Step_params(
    ω = ω_in,
    Q = 0, 
    Lb = 520,
    Ld = 130,
    xdₒᵤₜ = 650,
    name="Step0.5-Q0-omega$ω_in",
    mesh_file="floating_ice-Step05.json",
    n_elements = 4)
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




end
