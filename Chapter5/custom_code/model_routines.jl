function model_inhibition(D::DataFrames.DataFrame, M::XPP.Model, c::ConditionDef, name, p, time_in_noc; plotdata = true, xmin = 0, plotwhat = "CycB")

  if plotdata
    #Plot data
    for i in collect(1:length(c.data))
      plotcell(p[name], D, c.data[i], norm = 0.15)
    end
  end

  model_before(M, Dict("katt" => 0, "Inh" => 0), M.sims["Control"].I, time_in_noc);
  model_after(M, c.pars, Dict(), c.time, name, normalise = "CycB", t0 = "after");
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

function model_acutenoc(D::DataFrames.DataFrame, M::XPP.Model, c::ConditionDef, name, p, time_att)

  model_before(M, c.pars, M.sims["Control"].I, time_att);
  model_after(M, Dict("katt" => 0, "knoc" => M.pars["knoc_acute"]), Dict(), c.time, name, t0 = "before");

  for (k,i) in M.sims[name].D
    M.sims[name].D[k] = [M.sims["_before"].D[k]; M.sims[name].D[k]]
  end

  colors = Dict(
    12 => "#1D1A26",
    13 => "#293141",
    14 => "#304A5C",
    15 => "#316575",
    16 => "#2F828A",
    17 => "#30A09B",
    18 => "#30A09B",
    24 => "#1D1A26",
    27 => "#293141",
    30 => "#304A5C",
    33 => "#316575",
    36 => "#2F828A",
    39 => "#30A09B",
    42 => "#30A09B"
  )

  colorsukt = Dict(
    12 => "#321A0E",
    13 => "#522D1C",
    14 => "#73412B",
    15 => "#97563B",
    16 => "#BC6D4B",
    17 => "#E2845C",
    18 => "#E2845C",
    24 => "#321A0E",
    27 => "#522D1C",
    30 => "#73412B",
    33 => "#97563B",
    36 => "#BC6D4B",
    39 => "#E2845C",
    42 => "#E2845C",

    )

  # Plot simulation
  plotsim(M, p[name * "-CycB"], name, Dict("CycB" => colors[time_att]), [0, c.time], title = name * ": CycB-Degradation", linestyle = "-", ylabel = "CycB (a.u.)")

  y = collect(0:0.1:1.3)
  x = ones(length(y)) * time_att
  p[name * "-CycB"][:plot](x,y, color = "#adadad", linewidth = 0.5)

  sp = p[name * "-Mad2-localisation"]
  sp[:plot](x,y, color = "#adadad", linewidth = 0.5)
  x = M.sims[name].D["t"]
  y = M.sims[name].D["uKTa"]
  y_m2 = M.sims[name].D["M2loc"]
  sp[:plot](x, y_m2, linewidth = 2, color = colorsukt[time_att], linestyle = "-")
  sp[:set_xlim]([0,c.time])
  sp[:set_ylim](0,1)
  sp[:set_title](name * ": Mad2-Localisation")
  sp[:set_xlabel]("time (min)")
  sp[:set_ylabel]("uKT:Mad2 (a.u.)")
  sp[:legend](loc="best", title = "Nocodazole added after", frameon = false)
end

function main_modelling_figure(M::XPP.Model, conditions::Dict, D::DataFrame, title::AbstractString, celltype)
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
  model_condition(D, M, conditions["Control"], "Control", sp)
  model_condition(D, M, conditions["Nocodazole"], "Nocodazole", sp)
  subplots_adjust(left = 0.08, right = 0.95, top = 0.90, bottom = 0.15, hspace = 0.4, wspace = 0.4)
  fname = M.name * "_main_modelling_figure.pdf"
  return(M)
end

function supplementary_modelling_figure(M::XPP.Model, conditions::Dict, D::DataFrame, title::AbstractString)
  println(M.originalState[:init])
  sp = Dict()
  figure(figsize = (12,16))
  simulations = [
    "Nocodazole + 2.5 uM RO3306",
    "Nocodazole + 3.0 uM RO3306",
    "Nocodazole + 10.0 uM RO3306"
    ]
  for (i,v) in enumerate(simulations)
    sp[v] = subplot(3,1,i)
  end
  model_before(M, M.sims["Control"].P, M.sims["Control"].I, 40)
  # model_condition(D, M, conditions["Nocodazole"], "Nocodazole", sp)
  println(M.originalState[:init])
  for t in [1,130,280,445,640,870]
    model_inhibition(D, M, conditions["Nocodazole + 2.5 uM RO3306",], "Nocodazole + 2.5 uM RO3306", sp, t, plotdata = t == 1, xmin = -10, plotwhat = "CycB")
  end

  for t in [1,130,280,445,640,870]
    model_inhibition(D, M, conditions["Nocodazole + 3.0 uM RO3306"], "Nocodazole + 3.0 uM RO3306", sp, t, plotdata = t == 1, xmin = -10, plotwhat = "CycB")
  end

  for t in [1,130,280,445,640,870]
    model_inhibition(D, M, conditions["Nocodazole + 10.0 uM RO3306"], "Nocodazole + 10.0 uM RO3306", sp, t, plotdata = t == 1, xmin = -10, plotwhat = "CycB")
  end
  subplots_adjust(left = 0.08, right = 0.95, top = 0.95, bottom = 0.08, hspace = 0.4, wspace = 0.4)
  return(M)

  # sp = Dict()
  # figure(figsize = (12,16))
  # simulations = [
  #   "Nocodazole + 2.5 uM RO3306",
  #   "Nocodazole + 3.3 uM RO3306",
  #   "Nocodazole + 13.3 uM RO3306"
  #   ]
  # for (i,v) in enumerate(simulations)
  #   sp[v] = subplot(3,1,i)
  # end
  # model_before(M, M.sims["Control"].P, M.sims["Control"].I, 40)
  # # model_condition(D, M, conditions["Nocodazole"], "Nocodazole", sp)
  # println(M.originalState[:init])
  # for t in [1,130,280,445,640,870]
  #   model_inhibition(D, M, conditions["Nocodazole + 2.5 uM RO3306*",], "Nocodazole + 2.5 uM RO3306", sp, t, plotdata = t == 1, xmin = -10, plotwhat = "Sec")
  # end
  #
  # for t in [1,130,280,445,640,870]
  #   model_inhibition(D, M, conditions["Nocodazole + 3.0 uM RO3306*"], "Nocodazole + 3.3 uM RO3306", sp, t, plotdata = t == 1, xmin = -10, plotwhat = "Sec")
  # end
  #
  # for t in [1,130,280,445,640,870]
  #   model_inhibition(D, M, conditions["Nocodazole + 10.0 uM RO3306*"], "Nocodazole + 13.3 uM RO3306", sp, t, plotdata = t == 1, xmin = -10, plotwhat = "Sec")
  # end
  # subplots_adjust(left = 0.08, right = 0.95, top = 0.95, bottom = 0.08, hspace = 0.4, wspace = 0.4)
  # return(M)
end
