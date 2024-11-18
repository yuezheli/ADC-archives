# author: Yuezhe LI
# date: 9/28/2023
# purpose: to simulate bevacizumab human PK 

using DifferentialEquations 
using Plots
using DataFrames, CSV, DataFramesMeta

obs = CSV.read("data/shin_2020.csv",DataFrame);

include("param.jl");
include("helper.jl");
include("jonesODEs3Homo.jl");

p_homo = deepcopy(p_N_16);
p_homo.init_sR = 0.
p_homo.PS_Score = 5.

BW = 71.; # kg; 
V_Plasma = 3.126
Dose_in_mgkg = 3.; # [mg/kg]

tspan = (0.0, 2016.);  # hr

u0, dose_uM = jones_init(16, Dose_in_mgkg, BW, V_Plasma);

sol = solve(ODEProblem(jonesODE3Homo!,u0,tspan,p_homo), alg=AutoTsit5(Rosenbrock23(autodiff = false)), saveat = 1.);
cplasma = [sol.u[i].C_EXG_Plasma for i in 1:length(sol.t)] * MW_EDG; # [ug/L]

pHomo = plot(legend=:bottomleft, yaxis=:log10, dpi=1000, size=(350,350));
plot!(sol.t, cplasma/1E3, lw = 2, alpha = 0.6, label="sims" );
scatter!(obs.hour, obs.conc_ug_ml, label = "Bevacizumab obs", ma = 0.6, linewidth=0);
yticks!([1, 10, 100], ["1", "10", "100"]);
ylims!(0.6, 100);
xlabel!("time (hour)");
ylabel!("plasma conc (ug/mL)")

savefig(pHomo, "img/bevacizumab-homo.png");
