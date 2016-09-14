using XPP
using PyCall
PyDict(pyimport("matplotlib")["rcParams"])["font.family"] = ["Raleway"]
PyDict(pyimport("matplotlib")["rcParams"])["font.size"] = 18.0
using PyPlot
using DataFrames
using IDAT

type ConditionDef
  data::Array{AbstractString,1}
  pars::Dict
  time::Float64
end

#--------------------------------
# PLOTTING FUNCTIONS
#--------------------------------

# Function to plot timelapse data
function plotcell(p, D, id; norm = false, color = "#d2d2d2", label = "", alpha = 1, linewidth = 2)
  C = D[D[:cell_id] .== id, :];
  ana = C[:ana][1]
  x, y = get_xy(C,1)
  if norm != false
    y = (y - norm) ./ (1-norm)
  end
  p[:plot](x,y,color, alpha = alpha, label = label, linewidth = linewidth)
end

# function to plot simulation
function plotsim(M, sp, name, defs, xlim; ylim = [0,1.3], xlabel = "time (min)", label = "", title = "", linestyle = "-", leg_title = "", ylabel = "", leg = false)
  x = M.sims[name].D["t"]
  for (varname,color) in defs
    y = M.sims[name].D[varname]
    if label == "varname"
      sp[:plot](x, y, linewidth = 2, color = color, linestyle = linestyle, label = varname)
    else
      sp[:plot](x, y, linewidth = 2, color = color, linestyle = linestyle, label = label)
    end
  end
  sp[:set_xlim](xlim)
  sp[:set_ylim](ylim)
  sp[:set_title](title)
  sp[:set_xlabel](xlabel)
  sp[:set_ylabel](ylabel)
  if leg == true
    sp[:legend](loc="best", title = leg_title, frameon = false)
  end
end

#--------------------------------
# TWO-PART SIMULATIONS
#--------------------------------

# Simualte part before perturbation (i.e. parameter or initial condition changes)
function model_before(M, pars, ics, time_before)
  restoreModel!(M)
  #Reset inits
  for (k,v) in ics
    M.init[k] = v;
  end
  #Update parameters
  for (k,v) in pars;
    M.pars[k] = v
  end
  simulate!(M, "_before", collect(0:0.01:time_before));
end

# Simualte part following perturbation (i.e. parameter or initial condition changes)
function model_after(M, pars, ics, time_after, title; normalise = false, t0 = "before")
  # Update initials
  for (k,v) in M.init
    M.init[k] = M.sims["_before"].D[k][end]
  end
  # Change initials
  for (k,v) in ics
    M.init[k] = v
  end
  #Update parameters
  for (k,v) in pars;
    M.pars[k] = v
  end

  #Update kscycb
  M.pars["kscycb"] = M.sims["Control"].P["kscycb"] * exp(-M.pars["r_decay"] * M.sims["_before"].D["t"][end])
  M.spec["total"] = time_after
  simulate!(M, title, collect(0:0.01:time_after));

  if normalise != false
    M.sims[title].D[normalise] = M.sims[title].D[normalise] / M.sims["_before"].D[normalise][end];
    M.sims["_before"].D[normalise] = M.sims["_before"].D[normalise] / M.sims["_before"].D[normalise][end];
  end

  if t0 == "before"
    M.sims[title].D["t"] = M.sims[title].D["t"] + M.sims["_before"].D["t"][end];
  elseif t0 == "after"
    M.sims["_before"].D["t"] = M.sims["_before"].D["t"] - M.sims["_before"].D["t"][end];
  end
end

# Function to model an experiment as specified by a ConditionDef-object
function model_condition(D::DataFrames.DataFrame, M::XPP.Model, c::ConditionDef, name; p = false, uKT = 1.0, plot_data = true)
  restoreModel!(M)
  #Plot data
  if plot_data
    for i in collect(1:length(c.data))
      plotcell(p[name], D, c.data[i], norm = 0.15)
    end
    p[name][:set_xlim]([0, c.time])
    p[name][:set_ylim]([0, 2.4])
    p[name][:set_xlabel]("time (min)")
  end

  # Change parameters
  for (k,v) in c.pars
    M.pars[k] = v
  end

  M.init["uKTt"] = uKT;
  M.init["uKTa"] = uKT;

  M.spec["total"] = c.time
  #Simulate
  simulate!(M, name, collect(0:0.01:c.time))

  if p != false
    # Plot simulation
    S = Dict("CycB" => "k", "Sec" => "#75008c", "MCCt" => "#ff6600", "APC" => "#008000", "uKTa" => "#aa0000");
    plotsim(M, p[name], name, S, [0, c.time], label = "varname", title =  name, ylim = [0,2.2], ylabel = "a.u.")
  end
end
