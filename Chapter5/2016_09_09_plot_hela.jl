include("utilities.jl")
include("plot_data.jl")
include("model_routines.jl")

close("all")

#Load Data
S = Spec("specs.yml")
D = readtable(joinpath(S.path, "map_curated.csv"));
H = D[D[:celltype] .== "HeLa", :]
HR = H[H[:treatment] .== "RO-3306", :]
HN = H[H[:treatment] .== "noc_only", :]
HN = HN[HN[:ana] .!= 0.0, :]

R = D[D[:celltype] .== "RPE1", :]
RR = R[R[:treatment] .== "RO-3306", :]
RN = R[R[:treatment] .== "noc_only", :]
RN = RN[RN[:ana] .!= 0.0, :]

hro = Dict([d => [id for id in HR[HR[:dose_uM] .== d,:][:cell_id]] for d in unique(HR[:dose_uM])])
rro = Dict([d => [id for id in RR[RR[:dose_uM] .== d,:][:cell_id]] for d in unique(RR[:dose_uM])])
hro[0.0] = [id for id in HN[:cell_id]]
rro[0.0] = [id for id in RN[:cell_id]]

ylabel = Dict(
"HeLa" => "normalised securin-mEGFP intensity",
"RPE1" => "normalised cyclin B1-Venus intensity",
)

colors = Dict(
#RPE1
20.0 => "#075251",
10.0 => "#F4806D",
8.0 => "#B56B91",
4.0 => "#536686",
#---
3.0 => "#075251",
2.5 => "#F4806D",
#---
1.5 => "#075251",
0.75 => "#B56B91",
#HeLa
13.3 => "#B56B91",
6.6 => "#075251",
5.0 => "#536686",
#---
3.3 => "#075251",
2.5 => "#F4806D",
#---
1.7 => "#075251",
1.25 => "#F4806D",
0.8 => "#B56B91"
)

function plot_inh(hdose, rdose; xmax = 100, ymax = 1.2, w = 15, h = 6, orient = "side")

  figure(figsize = [w,h])
  if orient == "side"
    sp = Dict(
    "HeLa" => subplot(121),
    "RPE1" => subplot(122)
    )
  elseif orient == "top"
    sp = Dict(
    "HeLa" => subplot(211),
    "RPE1" => subplot(212)
    )
  end

  for k in keys(sp)
   sp[k][:set_xlim]([0,xmax])
   sp[k][:set_ylim]([0,ymax])
   sp[k][:set_xlabel]("time after RO-3306 addition (min)")
   sp[k][:set_ylabel](ylabel[k])
   sp[k][:set_title](k)
  end


  for d in sort([k for k in keys(hdose)])
    for (n,i) in enumerate(hro[d])
      plotcell(sp["HeLa"], D, i; color = hdose[d],
                label = n == 1 ? "$d\ uM" : "", alpha = 0.5)
    end
  end
  sp["HeLa"][:legend]()

  for d in sort([k for k in keys(rdose)])
    for (n,i) in enumerate(rro[d])
      plotcell(sp["RPE1"], D, i; color = rdose[d],
                label = n == 1 ? "$d\ uM" : "", alpha = 0.5)
    end
  end
  sp["RPE1"][:legend]()
end

function plot_noc()
  data = Dict(["HeLa" => hro[0.0], "RPE1" => rro[0.0]])
  for (k,v) in data
    fig = figure(figsize = [15,15])
    sp = Dict()
    for (n,i) in enumerate(v)
      sp[i] = subplot(4,2,ceil(n/4))
      sp[i][:set_xlim]([0,2500])
      sp[i][:set_ylim]([0,1.6])
      sp[i][:set_xlabel]("time (min)")
      sp[i][:set_ylabel]("a.u.")
      plotcell(sp[i], D, i; color = "k", alpha = 0.5)
    end
  end
end

# plot_inh(
#   Dict([
#     0.0 => "#c7c7c7",
#     13.3 => "#B56B91",
#     6.6 => "#F4806D"]),
#   Dict([
#     20.0 => "#075251",
#     10.0 => "#F4806D",
#     8.0 => "#B56B91",
#     4.0 => "#536686",
#     0.0 => "#c7c7c7"]);
#     xmax = 45
# )

# plot_inh(
#   Dict([
#     0.0 => "#c7c7c7",
#     2.5 => "#B56B91",
#     3.3 => "#075251",
#     # 6.6 => "#B56B91",
#     13.3 => "#F4806D"]),
#   Dict([
#     0.0 => "#c7c7c7",
#     2.5 => "#B56B91",
#     3.0 => "#075251",
#     # 4.0 => "#B56B91",
#     10.0 => "#F4806D",
#     # 20.0 => "#F4806D"
#     ]);
#   xmax = 150,
#   orient = "top",
#   w = 15,
#   h = 15
# )

plot_inh(
  Dict([
    0.0 => "#c7c7c7",
    1.25 => "#B56B91",
    # 1.7 => "#075251",
    # 13.3 => "#F4806D"
    ]),
  Dict([
    0.0 => "#c7c7c7",
    0.75 => "#B56B91",
    # 1.5 => "#075251",
    # 10.0 => "#F4806D",
    ]);
  xmax = 2500,
  ymax = 1.6,
  w = 15,
  orient = "top",
  h = 15
)

# plot_noc()
