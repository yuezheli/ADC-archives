# date: 7/1/2025 
# author: Yuezhe Li 
# purpose of this code: to track the source of free payload in tissue endothelial cells 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Parameters: @unpack
using DataFrames
using Plots
using CSV

include(@projectroot("Constants.jl"))
include(@projectroot("model/jones_homo.jl") )
include(@projectroot("model/param-pk.jl"))
include(@projectroot("model/init-pk.jl")) 
include(@projectroot("script/helper-infusion-dosing.jl"))
include(@projectroot("script/helper-human-organ-volumes.jl"))
include(@projectroot("script/helper-tissue-mass.jl"))

# set dose and time for simulation 
tspan = (0., hr_per_day*84);      # [hr]
Dose = 3.6  # [mg/kg]

dose_umol = Dose*BW*1E-3/MW_IGG*1E6; # [umol]

# set TMDD in plasma to 0 
p_notmdd = deepcopy(p_base); 
p_notmdd.init_sR = 0; 

# simulation (IV bolus dosing)
u0_tmp = jones_init(Dose*1E3, p_notmdd, BW, V_Plasma);
sol_tmp = solve(ODEProblem(jonesODEs_homo_tumor!, u0_tmp, tspan, p_notmdd), saveat = 0.25, alg = QNDF(autodiff=false), reltol = 1E-18);

total_degPL = sum(sol_tmp.u[end].DegPL.deg_membrane) + sum(sol_tmp.u[end].DegPL.deg_FcRn); # [uM]

organs = ["Lung", "Liver", "Heart", "Muscle", "Skin", "Adipose", "Bone", "Brain", "Kidney", "Long Intestine", "Small Intestine", "Pancreas", "Thymus", "Spleen", "Other"]; 

p_membrane = bar(organs,sol_tmp.u[end].DegPL.deg_membrane/total_degPL*100, xrotation = 45, color = :blue, alpha = 0.8, ylabel = "% of Degraded PL", label = "Membrane-mediated ADC uptake", dpi = 300, size = (400, 400))
p_FcRn = bar(organs,sol_tmp.u[end].DegPL.deg_FcRn/total_degPL*100, xrotation = 45, color = :red, alpha = 0.8, ylabel = "% of Degraded PL", label = "FcRn-mediated ADC uptake", dpi = 300, size = (400, 400))

savefig(p_membrane, @projectroot("deliv/figure/pl/deg-pl-distribution-membrane.png"));
savefig(p_FcRn, @projectroot("deliv/figure/pl/deg-pl-distribution-fcrn.png"));

# save output for visualization using ggplot
free_pl = DataFrame(
    Organ = organs, 
    deg_mem_perc = sol_tmp.u[end].DegPL.deg_membrane/total_degPL*100, 
    deg_fcrn_perc = sol_tmp.u[end].DegPL.deg_FcRn/total_degPL*100
);

CSV.write(@projectroot("data/sim/free-pl-tracking.csv"), free_pl)