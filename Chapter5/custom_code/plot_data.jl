using Interpolations

function main_data_figure(conditions::Dict, D::DataFrame, fname)
  today = Dates.today()
  try
    run(`mkdir $today\_files`)
  end
  try
    run(`mkdir $today\_files/dataplots`)
  end
  sp = Dict()
  figure(figsize = (12,8))
  simulations = [
  "Control",
  "Nocodazole",
  ]
  sp["Control"] = subplot(2,2,1)
  # sp["Control_r"] = subplot(2,2,3)
  sp["Nocodazole"] = subplot(2,2,2)
  # sp["Nocodazole_r"] = subplot(2,2,4)
  sp["Control_r2"] = subplot(2,2,3)
  sp["Nocodazole_r2"] = subplot(2,2,4)


  for (i,v) in enumerate(simulations)
    for n in conditions[v].data
      plottimecourse(sp, v, D, n)
    end
    sp[v][:set_xlim]([0, conditions[v].time])
    sp[v][:set_ylim]([0, 1.6])
    sp[v][:set_xlabel]("time (min)")
    # sp[v * "_r"][:set_xlim]([0, conditions[v].time])
    sp[v * "_r2"][:set_xlim]([0, 1.2])
    if v == "Control"
      # sp[v * "_r"][:set_ylim]([0, 0.25])
      sp[v * "_r2"][:set_ylim]([0, 0.25])
    else
      # sp[v * "_r"][:set_ylim]([0, 0.002])
      sp[v * "_r2"][:set_ylim]([0, 0.25])
    end
    # sp[v * "_r"][:set_xlabel]("time (min)")
    sp[v * "_r2"][:set_xlabel]("cyclin B")
  end
  subplots_adjust(left = 0.08, right = 0.95, top = 0.95, bottom = 0.08, hspace = 0.4, wspace = 0.4)
  savefig("$today\_files/dataplots/$fname")
  close()
end

# Function to plot timelapse data
function plottimecourse(p, v, D, id)
  C = D[D[:cell_id] .== id, :];
  data = readtable(C[1,:data_table])
  l = round(Int,C[1, :imin])
  u = round(Int,C[1, :imax])
  if v == "Nocodazole"
    step = 1
  else
    step = 1
  end
  x = data[:Time_aligned_][l:step:u]
  y = dropna(data[:Normalized])[l:step:u]
  dx = x[2:end] - x[1:end-1]
  dy = y[2:end] - y[1:end-1]
  deriv = - dy ./ dx
  spl = Spline1D(x, [0; deriv]; w=ones(length(x))/maximum(dx), k=3, bc="nearest", s=0)
  r = spl(x)
  # if v == "Nocodazole"
  #   w = 18
  #   rsmooth = zeros(length(r) - w + 1)
  #   for i in 1:w
  #     rsmooth += r[i:end-(w-i)]
  #   end
  #   tail = round(Int, w/2)
  #   rsmooth = [zeros(tail); rsmooth; zeros(tail)]
  #   r = rsmooth ./ w
  # else
  #   rsmooth = (r[1:end-3] + r[2:end-2] + r[3:end-1] + r[4:end])/4
  #   r = [0;0; rsmooth;0;0]
  # end
  ry = convert(Vector{Float64}, [0; r])
  rx = convert(Vector{Float64}, x)
  # ry = interpolate(ry, BSpline(Cubic(Flat())), OnGrid())

  m = round(Int, C[1, :meta]-l)
  a = round(Int, C[1, :ana]-l)

  println(id, "\t", l,"\t", m,"\t", a,"\t", u)

  if v == "Control"
    phases = Dict(
      "pm" => [1,m, "#aa0000"],
      "m" => [m,a, "#ff6600"],
      "a" => [a, length(x), "#008000"]
    )
  elseif v == "Nocodazole"
    phases = Dict(
      "pm" => [1,a, "#aa0000"],
      "a" => [a, length(x), "#008000"]
    )
  end

  phase_data = [k => get_range(x, y, rx, ry, i[1], i[2]) for (k,i) in phases]


  for (k,d) in phase_data
      p[v][:plot](d["x"], d["y"], phases[k][3], linewidth = 2)
      # p[v * "_r"][:plot](d["rx"], d["ry"], phases[k][3], linewidth = 2, alpha = 0.5)
      p[v * "_r2"][:plot](d["y"], d["ry"], phases[k][3], linewidth = 2, alpha = 0.5)
  end
  p[v][:scatter](phase_data["a"]["x"][1], phase_data["a"]["y"][1], color = "k", linewidth = 2)

end

function get_range(x, y, rx, ry, imin, imax)
  return(Dict("x" => x[imin:imax], "y" => y[imin:imax], "rx" => rx[imin:imax], "ry" => ry[imin:imax]))
  # return(Dict("x" => x[imin:imax], "y" => y[imin:imax]))
end
