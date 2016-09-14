function model_inhibition(D::DataFrames.DataFrame, M::XPP.Model, c::ConditionDef, name, p, time_in_noc; plotdata = true, xmin = 0, plotwhat = "CycB", normalise = "CycB")

  if plotdata
    #Plot data
    for i in collect(1:length(c.data))
      plotcell(p[name], D, c.data[i], norm = 0.15)
    end
  end

  model_before(M, Dict("katt" => 0, "Inh" => 0), M.sims["Control"].I, time_in_noc);
  model_after(M, c.pars, Dict(), c.time, name, normalise = normalise, t0 = "after");
  colors = Dict(
  1 => "#1D1A26",
  130 => "#293141",
  280 => "#304A5C",
  445 => "#316575",
  640 => "#2F828A",
  870 => "#30A09B"
  )
  plotsim(M, p[name], name, Dict(plotwhat => colors[time_in_noc]), [xmin, c.time], label = "$time_in_noc  min", title = name, leg_title = "Time in noc", ylim = [0,2.5])
  plotsim(M, p[name], "_before", Dict(plotwhat => colors[time_in_noc]), [xmin, c.time], title = name, linestyle = "--", leg_title = "Time in noc", ylim = [0,1.5])
  y = collect(0:0.1:2.5)
  x = zeros(length(y))
  p[name][:plot](x,y, color = "#adadad", linewidth = 0.5)
end

function model_ctrl_noc(M::XPP.Model, conditions::Dict, D::DataFrame, title::AbstractString, celltype)
  M.name = M.name * "_" * title
  checkpoint!(M);
  sp = Dict()
  figure(figsize = (14,6))
    simulations = [
    "Control",
    "Nocodazole",
    ]
  for (i,v) in enumerate(simulations)
    sp[v] = subplot(1,2,i)
  end
  model_condition(D, M, conditions["Control"], "Control", p = sp)
  model_condition(D, M, conditions["Nocodazole"], "Nocodazole", p = sp)
  subplots_adjust(left = 0.08, right = 0.95, top = 0.90, bottom = 0.15, hspace = 0.4, wspace = 0.4)
  fname = M.name * "_main_modelling_figure.pdf"
  return(M)
end

function model_inhibition_cycb(M::XPP.Model, conditions::Dict, D::DataFrame, title::AbstractString)
  model_condition(D, M, conditions["Control"], "Control", plot_data = false)
  println(M.originalState[:init])
  sp = Dict()
  figure(figsize = (12,16))
  simulations = [
    "Nocodazole + 2.5 uM RO3306",
    "2.5 uM not normalised",
    "Nocodazole + 3.0 uM RO3306",
    "3.0 uM not normalised",
    "Nocodazole + 10.0 uM RO3306",
    "10.0 uM not normalised"
    ]
  for (i,v) in enumerate(simulations)
    sp[v] = subplot(3,2,i)
  end
  model_before(M, M.sims["Control"].P, M.sims["Control"].I, 40)
  # model_condition(D, M, conditions["Nocodazole"], "Nocodazole", sp)
  println(M.originalState[:init])
  for t in [1,130,280,445,640,870]
    model_inhibition(D, M, conditions["Nocodazole + 2.5 uM RO3306",], "Nocodazole + 2.5 uM RO3306", sp, t, plotdata = t == 1, xmin = -10, plotwhat = "CycB", normalise = "CycB")
  end

  for t in [1,130,280,445,640,870]
    model_inhibition(D, M, conditions["Nocodazole + 3.0 uM RO3306"], "Nocodazole + 3.0 uM RO3306", sp, t, plotdata = t == 1, xmin = -10, plotwhat = "CycB", normalise = "CycB")
  end

  for t in [1,130,280,445,640,870]
    model_inhibition(D, M, conditions["Nocodazole + 10.0 uM RO3306"], "Nocodazole + 10.0 uM RO3306", sp, t, plotdata = t == 1, xmin = -10, plotwhat = "CycB", normalise = "CycB")
  end

  for t in [1,130,280,445,640,870]
    model_inhibition(D, M, conditions["2.5 uM not normalised"], "2.5 uM not normalised", sp, t, plotdata = false, xmin = -600, plotwhat = "CycB", normalise = false)
  end

  for t in [1,130,280,445,640,870]
    model_inhibition(D, M, conditions["3.0 uM not normalised"], "3.0 uM not normalised", sp, t, plotdata = false, xmin = -150, plotwhat = "CycB", normalise = false)
  end

  for t in [1,130,280,445,640,870]
    model_inhibition(D, M, conditions["10.0 uM not normalised"], "10.0 uM not normalised", sp, t, plotdata = false, xmin = -50, plotwhat = "CycB", normalise = false)
  end
  subplots_adjust(left = 0.08, right = 0.95, top = 0.90, bottom = 0.15, hspace = 0.4, wspace = 0.4)
  return(M)
end

function model_inhibition_sec(M::XPP.Model, conditions::Dict, D::DataFrame, title::AbstractString)
  model_condition(D, M, conditions["Control"], "Control", plot_data = false)
  sp = Dict()
  figure(figsize = (12,16))
  simulations = [
    "Nocodazole + 2.5 uM RO3306",
    "2.5 uM not normalised",
    "Nocodazole + 3.3 uM RO3306",
    "3.3 uM not normalised",
    "Nocodazole + 13.3 uM RO3306",
    "13.3 uM not normalised"
    ]
  for (i,v) in enumerate(simulations)
    sp[v] = subplot(3,2,i)
  end
  model_before(M, M.sims["Control"].P, M.sims["Control"].I, 40)
  # model_condition(D, M, conditions["Nocodazole"], "Nocodazole", sp)
  println(M.originalState[:init])
  for t in [1,130,280,445,640,870]
    model_inhibition(D, M, conditions["Nocodazole + 2.5 uM RO3306*",], "Nocodazole + 2.5 uM RO3306", sp, t, plotdata = t == 1, xmin = -10, plotwhat = "Sec")
  end

  for t in [1,130,280,445,640,870]
    model_inhibition(D, M, conditions["Nocodazole + 3.0 uM RO3306*"], "Nocodazole + 3.3 uM RO3306", sp, t, plotdata = t == 1, xmin = -10, plotwhat = "Sec")
  end

  for t in [1,130,280,445,640,870]
    model_inhibition(D, M, conditions["Nocodazole + 10.0 uM RO3306*"], "Nocodazole + 13.3 uM RO3306", sp, t, plotdata = t == 1, xmin = -10, plotwhat = "Sec")
  end
  for t in [1,130,280,445,640,870]
    model_inhibition(D, M, conditions["2.5 uM not normalised"], "2.5 uM not normalised", sp, t, plotdata = false, xmin = -600, plotwhat = "Sec", normalise = false)
  end

  for t in [1,130,280,445,640,870]
    model_inhibition(D, M, conditions["3.0 uM not normalised"], "3.3 uM not normalised", sp, t, plotdata = false, xmin = -150, plotwhat = "Sec", normalise = false)
  end

  for t in [1,130,280,445,640,870]
    model_inhibition(D, M, conditions["10.0 uM not normalised"], "13.3 uM not normalised", sp, t, plotdata = false, xmin = -50, plotwhat = "Sec", normalise = false)
  end
  subplots_adjust(left = 0.08, right = 0.95, top = 0.90, bottom = 0.15, hspace = 0.4, wspace = 0.4)
  return(M)
end
