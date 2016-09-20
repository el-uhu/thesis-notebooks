################################################################################
#                                                                              #
#       MODEL DEFINITIONS - MIRKOVIC ET. AL.                                   #
#                                                                              #
################################################################################

using XPP
export ctrl, tev, colchicine

function ctrl(typ)
    if typ == "stochastic"
        M = fromOdeFile("stochastic.ode")
        checkpoint!(M)
        return(M)
    elseif typ == "stochastic_basic"
        M = fromOdeFile("stochastic.ode")
        M.pars["Jamcc"] = 0
        M.pars["kamcc_cdk"] = 0
        M.pars["kamcc"] = 500
        M.pars["kloc_cdk"] = 0
        M.pars["kloc"] = 3
        checkpoint!(M)
        return(M)
    elseif typ == "stochastic_sac_only"
        M = fromOdeFile("stochastic.ode")
        M.pars["kloc_cdk"] = 0
        M.pars["kloc"] = 3
        checkpoint!(M)
        return(M)
    elseif typ == "deterministic"
        M = fromOdeFile("deterministic.ode")
        checkpoint!(M)
        return(M)
    end
end

function tev(typ)
    if typ == "stochastic"
        M = fromOdeFile("stochastic.ode")
        M.init["AurBc"] = 62/2
        M.init["AurBcla"] = 938/2
        M.pars["krem"] = 10
        M.pars["Stretch"] = 0.15
        checkpoint!(M)
        return(M)
    elseif typ == "stochastic_basic"
        M = fromOdeFile("stochastic.ode")
        M.pars["Jamcc"] = 0
        M.pars["kamcc_cdk"] = 0
        M.pars["kamcc"] = 500
        M.pars["kloc_cdk"] = 0
        M.pars["kloc"] = 3
        M.init["AurBc"] = 62/2
        M.init["AurBcla"] = 938/2
        M.pars["krem"] = 10
        M.pars["Stretch"] = 0.15
        checkpoint!(M)
        return(M)
    elseif typ == "stochastic_sac_only"
        M = fromOdeFile("stochastic.ode")
        M.pars["kloc_cdk"] = 0
        M.pars["kloc"] = 3
        M.init["AurBc"] = 62/2
        M.init["AurBcla"] = 938/2
        M.pars["krem"] = 10
        M.pars["Stretch"] = 0.15
        checkpoint!(M)
        return(M)
    elseif typ == "deterministic"
        M = fromOdeFile("deterministic.ode")
        M.init["AurBc"] = 0.5
        M.pars["krem"] = 10
        M.pars["Stretch"] = 0.15
        checkpoint!(M)
        return(M)
    end
end

function colchicine(typ)
    if typ == "stochastic"
        M = fromOdeFile("stochastic.ode")
        M.pars["katt"] = 0
        checkpoint!(M)
        return(M)
    elseif typ == "deterministic"
        M = fromOdeFile("deterministic.ode")
        M.pars["katt"] = 0
        checkpoint!(M)
        return(M)
    end

end
