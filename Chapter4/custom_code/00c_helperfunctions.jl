#Calculates attached kinteochores from unattached kinetochores variable
function get_aKT(M, simname, var)
    kttot = M.sims[simname].P["KTtot"]
    if var == "uKTtot"
        V = M.sims[simname].P["Vbulk"] + M.sims[simname].P["Vx"]
    elseif var == "uKT"
        V = M.sims[simname].P["Vx"]
    end
    aKT = ones(length(M.sims[simname].D[var])) * kttot * V - M.sims[simname].D[var]
    aKT = aKT ./ (kttot * V)
    aKT = aKT *100
    return(aKT)
end

#Calculate inhibition factor from residual activity
function calculateInhibition(targetActivity)
    Inh = 1/targetActivity - 1
    return(Inh)
end

#Determine time when a threshold in a given variable is crossed first
function thresholdCrossing(t, V, var, thres; direction = "from above")
    if direction == "from above"
        t_test = t[findfirst([v <= thres for v in V])]
    else
        t_test = t[findfirst([v >= thres for v in V])]
    end
    println(var, " crosses threshold of ", thres, " at t = ", t_test)
    return(t_test)
end

#Import rates_clustered.dat
function parseData(filename)
    f = open(filename)
    ks = ["CycB", "rate", "SD", "SEM"]
    D = Dict("CycB" => Float64[], "rate" => Float64[], "SD" => Float64[], "SEM" => Float64[])
    for l in eachline(f)
        #Remove newline and space at the end of each line
        l = split(string(l), "\n")[1];
        #Split line to get data points
        pts = split(l, "\t");
        #Loop over var-list & append to SimulationData.D-instance based on index
        for i in 1:4
            push!(D[ks[i]], float(pts[i]))
        end
    end
    return(D)
end

#Calculate rates
function rates(M, simname, var)
    V = M.sims[simname].D[var]
    t = M.sims[simname].D["tr"]
    V = V/maximum(V)
    rate = [- calcRate(V, t, i) for i in 4:length(t)]
    return(Dict(var => V[4:end], "rate" => rate))
end

calcRate(D,t, i) = (D[i]- D[i-3]) / (t[i]- t[i-3])
