#Subfunctions for plotting
plotvar(var, model, aest) = plot(model.sims[simname].D["tr"], model.sims[simname].D[var], linewidth = 2, color = aest, label = var)
plotAKT(var, model, aest) = plot(model.sims[simname].D["tr"], get_aKT(model, simname, var), linewidth = 2, color = aest, label = "aKT")

#Main function for generating a panel
function makePanel(M)
    figure(figsize = (20,9))
    ax1 = subplot2grid([2,7], [1,0], 1, 1)
    title("Control")
    xlim(0,t_ctrl)
    xticks([0,10,20])
    xlabel("time (min)")
    ylim(0,cyc_max)
    yticks([0,2500,5000,7500,10000])
    ylabel("Number of Molecules")

    for var in ["CycB", "APC"]
        plotvar(var, M["Control"], styles[var])
    end

    ax1_2 = ax1[:twinx]()
    xlim(0,t_ctrl)
    ylim(0, mcc_max)
    yticks([0,500,1000,1500,2000], ["", "", "", "", ""])

    plotvar("MCCfree", M["Control"], styles["MCCfree"])



    ax2 = subplot2grid([2,7], [1,1], 1, 2)
    title("Loss of cohesin")
    xlim(0,t_tev)
    xticks([0,10,20,30,40,50])
    xlabel("time (min)")
    ylim(0,cyc_max)
    yticks([0,2500,5000,7500,10000], ["", "", "", "", ""])

    for var in ["CycB", "APC"]
        plotvar(var, M["TEV"], styles[var])
    end

    ax2_2 = ax2[:twinx]()
    xlim(0,t_tev)
    ylim(0, mcc_max)
    yticks([0,500,1000,1500,2000])
    ylabel("Number of Molecules (MCC)")

    plotvar("MCCfree", M["TEV"], styles["MCCfree"])


    ax3 = subplot2grid([2,7], [0,4], 1, 1)
    title("Control")
    xlim(0,t_ctrl)
    xticks([0,10,20])
    xlabel("time (min)")
    ylim(0,aurb_max)
    ylabel("Number of Molecules")
    yticks([0,20, 40, 60, 80])

    for var in ["AurBc", "AurBa"]
        plotvar(var, M["Control"], styles[var])
    end

    ax3_2 = ax3[:twinx]()
    xlim(0,t_ctrl)
    ylim(0,akt_max)
    yticks([0,25,50,75,100], ["", "", "", "", ""])

    plotAKT("uKT", M["Control"], styles["uKT"])


    ax4 = subplot2grid([2,7], [0,5], 1, 2)
    title("Loss of cohesin")
    xlim(0,t_tev)
    xticks([0,10,20,30,40,50])
    xlabel("time (min)")
    ylim(0, aurb_max)
    yticks([0,20, 40, 60, 80], ["", "", "", "", ""])

    for var in ["AurBc", "AurBa"]
        plotvar(var, M["TEV"], styles[var])
    end

    ax4_2 = ax4[:twinx]()
    xlim(0,t_tev)
    ylim(0,akt_max)
    yticks([0,25, 50, 75, 100])
    ylabel("attachment (%)")

    plotAKT("uKT", M["TEV"], styles["uKT"])

    ax5 = subplot2grid([2,7], [1,4], 2, 1)
    xlim(0,t_ctrl)
    xticks([0,10,20])
    xlabel("time (min)")
    ylim(0,aurbt_max)
    ylabel("Number of Molecules")
    yticks([0, 250, 500, 750, 1000])

    for var in ["AurBctot", "AurBatot"]
        plotvar(var, M["Control"], styles[var])
    end

    ax5_2 = ax5[:twinx]()
    xlim(0,t_ctrl)
    ylim(0,akt_max)
    yticks([0,25,50,75,100], ["", "", "", "", ""])

    plotAKT("uKTtot", M["Control"], styles["uKT"])


    ax6 = subplot2grid([2,7], [1,5], 2, 2)
    xlim(0,t_tev)
    xticks([0,10,20,30,40,50])
    xlabel("time (min)")
    ylim(0,aurbt_max)
    yticks([0, 250, 500, 750, 1000], ["", "", "", "", ""])

    for var in ["AurBctot", "AurBatot"]
        plotvar(var, M["TEV"], styles[var])
    end

    ax6_2 = ax6[:twinx]()
    xlim(0,t_tev)
    ylim(0,akt_max)
    yticks([0,25, 50, 75, 100])
    ylabel("attachment (%)")

    plotAKT("uKTtot", M["TEV"], styles["uKT"])

    subplots_adjust(left = 0.07, right = 0.93, hspace = 0.33, wspace = 0.10, top = 0.90, bottom = 0.10)
end
