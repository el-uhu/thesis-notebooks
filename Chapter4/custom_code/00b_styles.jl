################################################################################
#                                                                              #
#       PLOT COLORS AND FONTS                                                  #
#                                                                              #
################################################################################

styles = Dict(
    "CycB" => "#04328b",
    "APC" => "#670a04",
    "MCCfree" => "#288805",
    "AurBc" => "#db0303",
    "AurBctot" => "#db0303",
    "AurBa" => "#37c4a6",
    "AurBatot" => "#37c4a6",
    "uKT" => "#440091",
    "uKTtot" => "#440091",
)

fonts = Dict(
    "subplt_title" => Dict("family" => "Arial", "weight" => "bold", "size" => "20"),
    "subplt_axis" => Dict("family" => "Arial", "weight" => "bold", "size" => "17"),
    "base" => Dict("family" => "Arial", "size" => 17),
)

using PyCall
PyDict(pyimport("matplotlib")["rcParams"])["font.family"] = [fonts["base"]["family"]]
PyDict(pyimport("matplotlib")["rcParams"])["font.size"] = [fonts["base"]["size"]]

global cyc_max =11000
global mcc_max = 2200
global aurb_max = 88
global aurbt_max = 1100
global akt_max = 110
global t_ctrl = 30
global t_tev = 60
