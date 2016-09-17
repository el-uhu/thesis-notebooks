using PyCall
PyDict(pyimport("matplotlib")["rcParams"])["font.family"] = ["Raleway"]
PyDict(pyimport("matplotlib")["rcParams"])["font.size"] = 18.0
using PyPlot
using DataFrames

function do_nothing(x)
  return x
end

function convert_synth_to_time(x)
  # synth = kscycb * exp(-r_decay * t)
  # synth/kscycb = exp(-r_decay * t)
  # log(synth/kscycb) = -r_decay * t
  # t = - log(synth/kscycb)/r_decay
  t = - log(x./0.005)/0.0008
  return(t)
end

function convert_inh_to_activity(x)
  # Cdk1 = 1/(1 + Inh)
  return(1./(1+x) .* 100)
end

function plot_bifurcation(path; xlim = false, ylim = false, xlabel = false, xf = do_nothing, yf = do_nothing, xscale = "linear", sp = false, color = "r")
  names = split(split(path, "_")[1], ".")[1]
  xtitle = split(names, "-")[1]
  ytitle = split(names, "-")[2]
  D = readtable(path, separator = ' ', header = false, names = [symbol(xtitle), symbol(ytitle), symbol(ytitle * "_2"), symbol("stable_1"), symbol("stable_2")])
  if sp == false
    figure()
    sp = subplot(111)
  end
  #get first branch
  i = findfirst(D[:stable_1][2:end], 2) + 1
  x_1 = xf(D[symbol(xtitle)][1:i])
  y_1 = yf(D[symbol(ytitle)][1:i])
  #get second branch
  j = findfirst(D[:stable_1][i+1:end], 1) - 2
  x_2 = xf(D[symbol(xtitle)][i:i+j])
  y_2 = yf(D[symbol(ytitle)][i:i+j])
  #get third branch
  x_3 = xf(D[symbol(xtitle)][i+j+1:end])
  y_3 = yf(D[symbol(ytitle)][i+j+1:end])

  sp[:plot](x_1, y_1, color = color, linewidth = 2)
  sp[:plot](x_2, y_2, color = color, linewidth = 2, linestyle = "--")
  sp[:plot](x_3, y_3, color = color, linewidth = 2)

  if xlabel == false
    xlabel = xtitle
  end
  sp[:set_xlabel](xlabel)
  sp[:set_xscale](xscale)
  sp[:set_ylabel](ytitle)
  if xlim != false
    sp[:set_xlim](xlim)
  end
  if ylim != false
    sp[:set_ylim](ylim)
  end
  return(D)
end

colors = Dict(
"kscycb0005" => "#1D1A26",
"kscycb00045" => "#293141",
"kscycb0004" => "#304A5C",
"kscycb00035" => "#316575",
"kscycb0003" => "#2F828A",
"kscycb000255" => "#30A09B"
)

fig = figure(figsize = [12,6])
sp = subplot(121)
for (n, c) in colors
  path = "uKT-CycB_" * n * ".dat"
  plot_bifurcation(path; xlim = [0,1], ylim = [0,2.25], xf = do_nothing, yf = do_nothing, xscale = "log", sp = sp, color = c)
end

sp = subplot(122)
for (n, c) in colors
  path = "uKT-CycB_" * n * ".dat"
  plot_bifurcation(path; xlim = [0,1], ylim = [0,2.25], xf = do_nothing, yf = do_nothing, xscale = "linear", sp = sp, color = c)
end
# savefig("uKT-CycB.eps")
# close()

fig = figure()
sp = subplot(111)
for (n, c) in colors
  path = "Inh-CycB_" * n * ".dat"
  plot_bifurcation(path; xlim = [0,100], ylim = [0,2.25], xf = convert_inh_to_activity, yf = do_nothing, xscale = "linear", sp = sp, color = c, xlabel = "maximal Cdk1-activity after inhibition (%)")
end
# savefig("Inh-CycB.eps")
# close()

fig = figure(figsize = [12,6])
sp = subplot(121)
plot_bifurcation("kscycb-CycB_upperbranch.dat"; xlim = [0,0.005], ylim = [0,2.25], xf = do_nothing, yf = do_nothing, xscale = "linear", sp = sp, color = "#1D1A26")

sp = subplot( 122)
plot_bifurcation("kscycb-CycB_upperbranch.dat"; xlim = [0,2000], ylim = [0,2.25], xf = convert_synth_to_time, yf = do_nothing, xscale = "linear", sp = sp, color = "#1D1A26", xlabel = "time (min)")
# savefig("time-CycB.eps")
# close()
