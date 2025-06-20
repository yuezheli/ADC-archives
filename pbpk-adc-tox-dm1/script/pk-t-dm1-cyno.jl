# date: 1/3/2025
# author: Yuezhe Li 
# purpose of this code: to generate figure for T-DM1 cyno PK fit 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using DataFrames, CSV

obs = CSV.read(@projectroot("data/t-dm1-cyno.csv"),DataFrame);

using DifferentialEquations 
using Plots

include(@projectroot("model/jones_cyno.jl") );
include(@projectroot("model/init-pk.jl")); 

# ADC IV Dose
Dose_in_mgkg = 10.; # [mg/kg]

BW_cyno = 6.2          # cyno body weight, [kg]
V_Plasma_cyno = 0.187  # [L] 

p_cyno = ComponentArray(
    PS_Score = 6.,              # AC-SINS score, fitted using clincial PK
    PS_kd = -1.0,               # placeholder params; disabled; 
    KD6_WT = 700.0,             # binding affinity between ADC and FcRn, [nM], https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12461
    infusion = 0.               # Related to ADC infusion; dummy variable 
);

tspan = (0.0, 42*DayToHour);  # hr
u0, _ = jones_init(Dose_in_mgkg*1E3, p_cyno, BW_cyno, V_Plasma_cyno);

sol = solve(ODEProblem(jones_cyno!,u0,tspan,p_cyno), alg=AutoTsit5(Rosenbrock23(autodiff = false)), saveat = 1.);
cplasma = [sol.u[i].C_EXG_Plasma for i in 1:length(sol.t)] * MW; # [ug/L]

pMonkey = plot(legend=:topright, yaxis=:log10, dpi=1000, size=(350,350));
plot!(sol.t./DayToHour, cplasma, lw = 2, alpha = 0.6, label = "Simulation from PBPK model");
scatter!(obs.day, obs.adc_conc_ng_mL, label = "T-DM1 C. Macaque", ma = 0.6, linewidth=0);
xticks!([0, 7, 14, 21, 28, 35, 42], ["0", "7", "14", "21", "28", "35", "42"]);
xlabel!("Time (days)");
ylabel!("plasma concentration (ng/mL)");
display(pMonkey)

savefig(pMonkey, "deliv/figure/pk/t-dm1-cyno.png")

