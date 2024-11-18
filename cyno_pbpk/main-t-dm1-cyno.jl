# author: Yuezhe LI 
# date: 9/27/23
# purpose of this code: to identify endothelial scaling factor for cyno macaque for the PBPK model for mAb distribution (https://pubmed.ncbi.nlm.nih.gov/31464379/)
# ACSINS score upper bound set based on Jain et al., 2017: https://www.pnas.org/doi/10.1073/pnas.1616408114
# monkey and human endothelial density difference between 1.1 and 2.5; Kelley et al., 1984: https://pubmed.ncbi.nlm.nih.gov/6381373/

using DifferentialEquations 
using Plots
using DataFrames, CSV, DataFramesMeta

obs = CSV.read("data/t-dm1-cyno.csv",DataFrame);

include("param.jl");
include("helper.jl");
include("jonesODEs3Monkey.jl");

# mAb IV Dose
Dose_in_mgkg = 10.; # [mg/kg]

BW = 6.2
V_Plasma = 0.187 

# initialization
tspan = (0.0, 42*DayToHour);  # hr
u0, dose_uM = jones_init(15, Dose_in_mgkg, BW, V_Plasma);

p_cyno = deepcopy(p_N_15);
p_cyno.endothelial_scaling_factor = 60.

sol = solve(ODEProblem(jonesODEs3Monkey!,u0,tspan,p_cyno), alg=AutoTsit5(Rosenbrock23(autodiff = false)), saveat = 1.);
cplasma = [sol.u[i].C_EXG_Plasma for i in 1:length(sol.t)] * MW_EDG; # [ug/L]

pMonkey = plot(legend=:topright, yaxis=:log10, dpi=1000, size=(350,350));
plot!(sol.t./DayToHour, cplasma, lw = 2, alpha = 0.6, label="endothelial scaling factor = "* string(p_cyno.endothelial_scaling_factor) );
scatter!(obs.day, obs.adc_conc_ng_mL, label = "T-DM1 C. Macaque", ma = 0.6, linewidth=0);
xticks!([0, 7, 14, 21], ["0", "7", "14", "21"]);
xlabel!("time (days)");
ylabel!("plasma conc (ng/mL)")

savefig(pMonkey, "img/t-dm1-cyno.png")
