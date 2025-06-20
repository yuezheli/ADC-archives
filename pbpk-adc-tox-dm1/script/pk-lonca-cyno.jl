# date: 1/21/2025
# author: Yuezhe Li 
# purposd of this code: to generate PK fit for Lonca (using cyno)

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)

# ADCT-402 cyno data; 
# data obtained from Zammarchi et al., 2018; # https://www.sciencedirect.com/science/article/pii/S0006497120324174
using DataFrames
adct_402_cyno = DataFrame(time_hr=[0, 6.03, 25, 47, 69.2, 110.17, 160.44, 210, 327, 496.5], ADC_conc_ngmL=[15848.99, 13332.84, 10111.78, 7668.62, 6677.39, 4891, 4112.9, 3580, 2531.15, 1788.42]); 

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Plots

include(@projectroot("model/jones_cyno.jl") );
include(@projectroot("model/init-pk.jl")); 

# ADC IV Dose
Dose_in_mgkg = 0.6; # [mg/kg]

BW_cyno = 6.2          # cyno body weight, [kg]
V_Plasma_cyno = 0.187  # [L] 

p_cyno = ComponentArray(
    PS_Score = 6.,              # AC-SINS score, tuned for cyno PK
    PS_kd = -1.0,               # placeholder params; disabled; 
    KD6_WT = 700.0,             # binding affinity between ADC and FcRn, [nM], https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12461
    infusion = 0.               # Related to ADC infusion; dummy variable 
);

tspan = (0.0, 42*DayToHour);  # hr
u0, _ = jones_init(Dose_in_mgkg*1E3, p_cyno, BW_cyno, V_Plasma_cyno);

sol = solve(ODEProblem(jones_cyno!,u0,tspan,p_cyno), alg=AutoTsit5(Rosenbrock23(autodiff = false)), saveat = 1.);
cplasma = [sol.u[i].C_EXG_Plasma for i in 1:length(sol.t)] * MW; # [ug/L]

pMonkey = plot(legend=:topright, yaxis=:log10, ylims = [1E2, 1E6], dpi=1000, size=(400,400));
plot!(sol.t, cplasma, lw = 2, alpha = 0.6, label = "Simulation from PBPK model");
scatter!(adct_402_cyno.time_hr, adct_402_cyno.ADC_conc_ngmL, label = "ADCT-402 C. Macaque", ma = 0.6, linewidth=0);
hline!([212.9], linestyle = :dashdot, color = "black", alpha = 0.5, label = "LLOQ of total Ab");
xlabel!("Time (hr)");
ylabel!("plasma conc (ng/mL)");
display(pMonkey)

savefig(pMonkey, "deliv/figure/pk/ADCT-402-cyno.png")


