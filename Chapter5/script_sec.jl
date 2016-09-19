include("custom_code/utilities.jl")
include("custom_code/plot_data.jl")
include("custom_code/model_routines.jl")
close("all")

#Load Data
S = Spec("specs.yml")
D = readtable(joinpath(S.path, "map_curated.csv"));
H = D[D[:celltype] .== "HeLa", :]
HR = H[H[:treatment] .== "RO-3306", :]
HN = H[H[:treatment] .== "noc_only", :]
HN = HN[HN[:ana] .!= 0.0, :]
hctrl = [id for id in H[H[:treatment] .== "untreated", :cell_id]]

R = D[D[:celltype] .== "RPE1", :]
RR = R[R[:treatment] .== "RO-3306", :]
RN = R[R[:treatment] .== "noc_only", :]
RN = RN[RN[:ana] .!= 0.0, :]

hro = Dict([d => [id for id in HR[HR[:dose_uM] .== d,:][:cell_id]] for d in unique(HR[:dose_uM])])
rro = Dict([d => [id for id in RR[RR[:dose_uM] .== d,:][:cell_id]] for d in unique(RR[:dose_uM])])
hro[0.0] = [id for id in HN[:cell_id]]
rro[0.0] = [id for id in RN[:cell_id]]

a = []
for i in rro[0.75]
  if !(contains(i, "x") || contains(i, "!"))
    a = [a;i]
  end
end
rro[0.75] = a

rctrl = [id for id in R[R[:treatment] .== "untreated", :cell_id]]

#Specify model
odes = Dict(
  "CycB" => "synth - (kdcycb + kdcycb_c*APC)*CycB",
  "Sec" => "synths - (kdsec + kdsec_c*APC)*Sec",
  "APC" => "(kdiss + kimcc + kimcc_p31 * 1/(Jimcc^N + Cdk1^N))*APCMCC - (kass + kass_cdk*Cdk1)*APC*MCCfree",
  "MCCt" => "kamcc*M2 * uKTa  - (kimcc + kimcc_p31 * 1/(Jimcc^N + Cdk1^N))*MCCfree - (kimcc + kimcc_p31 * 1/(Jimcc^N + Cdk1^N)) * APCMCC",
  "Inhe" => "a * (Inh - Inhe)",
  "uKTa" => "ka*Cdk1*(uKTt - uKTa) - ki *  uKTa",
  "uKTt"=>"Noc * (1-uKTt) - katt*uKTt",
  "Noc" => "knoc * (Nocmax - Noc) - knoci * Noc"
);
aux = Dict(
  "M2loc" => "uKTa",
);
alg = Dict(
  "synth" => "kscycb * exp(-r_decay * t)",
  "synths" => "kssec * exp(-r_decays * t)",
  "Cdk1" => "CycB/(1 + Inhe)",
  "M2" => "MCCtot-MCCt",
  "APCMCC" => "1-APC",
  "MCCfree" => "MCCt-APCMCC",
);
pars = Dict(
  "Inh" => 0.0,
  "Inh_ro0_75" => 0.4,
  "Inh_ro2_5" => 0.75,
  "Inh_ro3" => 2.2,
  "Inh_ro10" => 5,
  "katt" => 0.0,
  "katt_ctrl" => 0.6,
  "kdcycb" => 0.002,
  "kscycb" => 0.005,
  "kdcycb_c" => 1.5,
  "kdsec" => 0.001,
  "kssec" => 0.0025,
  "kdsec_c" => 2.5,
  "r_decay" => 0.0008,
  "r_decays" => 0.0016,
  "kdiss" => 0.2,
  "kass" => 100,
  "kass_cdk" => 350,
  "kamcc" => 100,
  "kimcc" => 0.1,
  "kimcc_p31" => 0.08,
  "knoc" => 0,
  "knoci" => 0,
  "knoc_acute" => 0.5,
  "Nocmax" => 0.5,
  "ka" => 40,
  "ki" => 40,
  "Jimcc" => 0.09,
  "N" => 2.6,
  "MCCtot" => 2.4,
  "a" => 0.12,
);
inits = Dict(
"Inhe" => 0.0,
"APC" => 0.0,
"CycB" => 1.0,
"Sec" => 1,
"MCCt" => 2.2,
"uKTa" => 1.0,
"uKTt" => 1.0,
"Noc" => 0.0,
);
spec = Dict(
#Settings:
  "meth" => "stiff",
  "total" => 3500,
  "maxstor" => "2000000000",
);
t = Dates.today()
M = Model(odes, inits, pars, name = "$t\_00_RPE1", alg = alg, spec = spec, aux = aux)
checkpoint!(M)

# Define conditions
conditions = Dict(
"Control" => ConditionDef(rctrl, Dict("katt" =>pars["katt_ctrl"]), 30),
"Nocodazole" => ConditionDef(rro[0.0], Dict("katt" => pars["katt"]), 2000.0),
"Nocodazole + 2.5 uM RO3306" => ConditionDef(rro[2.5], Dict("katt" => pars["katt"], "Inh" => pars["Inh_ro2_5"]), 150.0),
"Nocodazole + 3.0 uM RO3306" => ConditionDef(rro[3.0], Dict("katt" => pars["katt"], "Inh" => pars["Inh_ro3"]), 150.0),
"Nocodazole + 10.0 uM RO3306" => ConditionDef(rro[10.0], Dict("katt" => pars["katt"], "Inh" => pars["Inh_ro10"]),
 150.0),
"2.5 uM not normalised" => ConditionDef(rro[2.5], Dict("katt" => pars["katt"], "Inh" => pars["Inh_ro2_5"]), 600.0),
"3.0 uM not normalised" => ConditionDef(rro[3.0], Dict("katt" => pars["katt"], "Inh" => pars["Inh_ro3"]), 150.0),
"10.0 uM not normalised" => ConditionDef(rro[10.0], Dict("katt" => pars["katt"], "Inh" => pars["Inh_ro10"]),
50.0),
"Nocodazole + 2.5 uM RO3306*" => ConditionDef(hro[2.5], Dict("katt" => pars["katt"], "Inh" => pars["Inh_ro2_5"]), 150.0),
"Nocodazole + 3.0 uM RO3306*" => ConditionDef(hro[3.3], Dict("katt" => pars["katt"], "Inh" => pars["Inh_ro3"]), 150.0),
"Nocodazole + 10.0 uM RO3306*" => ConditionDef(hro[13.3], Dict("katt" => pars["katt"], "Inh" => pars["Inh_ro10"]),
150.0),
"Nocodazole + 0.75 uM RO3306" => ConditionDef(rro[0.75], Dict("katt" => pars["katt"], "Inh" => pars["Inh_ro0_75"]), 2000),
"0.75 uM not normalised" => ConditionDef(rro[0.75], Dict("katt" => pars["katt"], "Inh" => pars["Inh_ro0_75"]), 2000),
)

# M = model_ctrl_noc(M, conditions, D, "modelplots", "RPE1")
# M = model_inhibition_cycb(M, conditions, D, "modelplots")
# M = model_inhibition_sec(M, conditions, D, "modelplots")
M = model_prediction(M, conditions, D, "modelplots")
# main_data_figure(conditions, D, "RPE1.pdf")
