# author: Yuezhe LI 
# date: 9/28/2023

using Pkg; Pkg.activate("");
# Pkg.instantiate();

using DifferentialEquations 
using Plots
using DataFrames, CSV, DataFramesMeta

obs = CSV.read("data/zalevsky_2010.csv",DataFrame);

include("param.jl");
include("helper.jl");
include("jonesODEs3Monkey.jl");

# mAb IV Dose
Dose_in_mgkg = 4.; # [mg/kg]

BW = 6.2
V_Plasma = 0.187 

# initialization
tspan = (0.0, 90*DayToHour);  # hr
u0, dose_uM = jones_init(15, Dose_in_mgkg, BW, V_Plasma);

p_cyno = deepcopy(p_N_15);
p_cyno.endothelial_scaling_factor = 60.
p_cyno.PS_Score = 5.

sol = solve(ODEProblem(jonesODEs3Monkey!,u0,tspan,p_cyno), alg=AutoTsit5(Rosenbrock23(autodiff = false)), saveat = 1.);
cplasma = [sol.u[i].C_EXG_Plasma for i in 1:length(sol.t)] * MW_EDG; # [ug/L]

pMonkey = plot(legend=:topright, yaxis=:log10, dpi=1000, size=(350,350));
plot!(sol.t./DayToHour, cplasma/1E3, lw = 2, alpha = 0.6, label="endothelial scaling factor = "* string(p_cyno.endothelial_scaling_factor) );
scatter!(obs.day, obs.conc_ug_ml, label = "Bevacizumab C. Macaque", ma = 0.6, linewidth=0);
xlabel!("time (days)");
ylabel!("plasma conc (ug/mL)")

savefig(pMonkey, "img/bevacizumab-cyno.png")
