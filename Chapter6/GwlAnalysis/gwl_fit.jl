#Load required packages
using XPP
using DataFrames
using PyPlot
using Optim

# Differential equation
odes = Dict(
    "pS" => "(k_pind + k_gwl * Gwl) * (St - pS) - (k_ind + k_b55 * B)*pS",
    "pEt" => "k_egwl * (Et - pEt) * Gwl - (k_diss + k_cat) * pEB",
    "pEB" => "k_ass * B * pE - (k_diss + k_cat) * pEB"
);
#Algebraic equations
alg = Dict(
    "B" => "Bt - pEB",
    "pE" => "pEt - pEB"
)
# Initial condition
init = Dict(
    "pS" => 0.0,
    "pEt" => 0.0,
    "pEB" => 0.0,
);
#Parameters
pars = Dict(
    "k_pind" => 1.0,
    "k_gwl" => 0.0,
    "k_ind" => 0.0,
    "k_b55" => 0.0,
    "k_egwl" => 0.008,
    "k_diss" => 0.001,
    "k_ass" => 8.2551,
    "k_cat" => 0.165,
    "Gwl" => 1.0,
    "Bt" => 0.7303,
    "Et" => 2.3,
    "St" => 1.0,
);
M = Model(odes, init, pars, name = "Phosphorylation-Model", alg = alg);

function simulateConditions(M, D, p; out = "rmsd_only")
    t1 = 45;
    Gwl=1.0;
    Bt=0.7303;
    Et=2.3;
    Gwl_depletion=0.1;
    B55_depletion=0.3099;
    M.pars["k_gwl"] = p[1];
    M.pars["k_b55"] = D[:k_ds][1];
    M.pars["k_ind"] = D[:k_bg][1];
    M.pars["k_pind"] = p[2];

    checkpoint!(M)
    # Simulate Control
    M.pars["Bt"] = Bt
    M.pars["Gwl"] = Gwl
    simulate!(M, "Control_1", 0:0.01:t1)

    # Simulate B55-Depletion
    M.pars["Bt"] = Bt * B55_depletion
    M.pars["Gwl"] = Gwl
    simulate!(M, "B55_1", 0:0.01:t1)

    # Simulate Gwl-Depletion
    M.pars["Gwl"] = Gwl * Gwl_depletion
    M.pars["Bt"] = Bt
    simulate!(M, "Gwl_1", 0:0.01:t1)

    C0 = M.sims["Control_1"].D["pS"][end];
    G0 = M.sims["Gwl_1"].D["pS"][end];
    B0 = M.sims["B55_1"].D["pS"][end];
    sim = [B0/C0, G0/C0];
    C0d = D[:Control_0][1];
    G0d = (1 /D[:CM_Con_H_GWL_L_][1] + D[:CM_Con_L_GWL_H_][1])/2;
    B0d = (1/ D[:CM_Con_H_B55_L_][1] + D[:CM_Con_L_B55_H_][1])/2;
    data = [B0d, G0d]
    rmsd = sqrt(sum((sim - data).^2))/2
    if out == "rmsd_only"
      return(rmsd)
    else
      return(rmsd, sim, data)
    end
end


D = readtable("substrates.csv");

D[:k_gwl] = zeros(size(D)[1]);
D[:k_pind] = zeros(size(D)[1]);
D[:rmsd] = zeros(size(D)[1]);

for i in 1:size(D)[1]
  if D[i, :rmsd] == 0.0
    println(D[i, :Gene_names])
    println("k_b55 = $(D[i, :k_ds])")
    println("k_ind = $(D[i, :k_bg])")
    println((1/D[i,:CM_Con_H_B55_L_][1] + D[i, :CM_Con_L_B55_H_][1])/2);
    println((1/D[i, :CM_Con_H_GWL_L_][1] + D[i, :CM_Con_L_GWL_H_][1])/2);
    println("Datapoint $i of $(size(D)[1])")
    f(p) = simulateConditions(M, D[i,:], p)
    g = DifferentiableFunction(f)
    initial_guess = [0.00001, 0.999];
    lower_bound = [0.0, 0.0];
    upper_bound = [1.0, 1.0];
    r = fminbox(g, initial_guess, lower_bound,  upper_bound, xtol = 1e-5, ftol = 1e-5, grtol= 1e-2, show_trace=true)
    # r = optimize(f, initial_guess);
    D[i, :k_gwl] = r.minimum[1];
    D[i, :k_pind] = r.minimum[2];
    D[i, :rmsd] = r.f_minimum;
    println("\t rmsd: $(r.f_minimum)\t parameters: $(r.minimum)")
    writetable("substrates_gwl.csv", D)
  end
end
