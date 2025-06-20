# date: 2/8/2025
# author: Yuezhe Li 
# purpose of this code: to track source of payload 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)

using DifferentialEquations, ComponentArrays
using Plots, DataFrames, LaTeXStrings
using Statistics

include(@projectroot("model/jones_homo.jl") );
include(@projectroot("model/param-pk.jl")); 
include(@projectroot("model/init-pk.jl")); 
include(@projectroot("script/helper.jl"));

p_tdm1 = deepcopy(p_base);
p_tdm1.PS_Score = 6.
p_tdm1.init_sR = 0.004


first_dose = 3.6  # [mg/kg]
tspan = (0., DayToHour*21);      # [hr]

u0_tdm1, _ = jones_init(first_dose*1E3, p_tdm1);
prob_tdm1 = ODEProblem(jonesODEs_homo_tumor!, u0_tdm1, tspan, p_tdm1);

sol_tmp = solve(prob_tdm1, saveat = 1., alg = QNDF(autodiff=false)); 

total_degPL = sum(sol_tmp.u[end].DegPL.deg_membrane) + sum(sol_tmp.u[end].DegPL.deg_FcRn); # [uM]

organs = ["Lung", "Liver", "Heart", "Muscle", "Skin", "Adipose", "Bone", "Brain", "Kidney", "Long Intestine", "Small Intestine", "Pancreas", "Thymus", "Spleen", "Other"]; 

#dep_pl = DataFrame(Organ = [organs;organs], PL = [sol_tmp.u[end].DegPL.deg_membrane/total_degPL*100; sol_tmp.u[end].DegPL.deg_FcRn/total_degPL*100], Source = repeat(["Membrane-mediated", "FcRn-mediated"], inner = 15));

#bar(dep_pl.Organ, dep_pl.PL, group = dep_pl.Source, ylims = [0.0001, 100], yaxis = :log, xrotation = 45, ylabel = "% of Degraded PL")

p_membrane = bar(organs,sol_tmp.u[end].DegPL.deg_membrane/total_degPL*100, xrotation = 45, color = :blue, alpha = 0.8, ylabel = "% of Degraded PL", label = "Membrane-mediated ADC uptake", dpi = 1000, size = (400, 400))
p_FcRn = bar(organs,sol_tmp.u[end].DegPL.deg_FcRn/total_degPL*100, xrotation = 45, color = :red, alpha = 0.8, ylabel = "% of Degraded PL", label = "FcRn-mediated ADC uptake", dpi = 1000, size = (400, 400))

savefig(p_membrane, @projectroot("deliv/figure/payload/deg-pl-distribution-membrane.png"));
savefig(p_FcRn, @projectroot("deliv/figure/payload/deg-pl-distribution-fcrn.png"));
