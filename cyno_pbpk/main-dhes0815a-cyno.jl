# to simulate DHES0815a cyno PK

using DifferentialEquations 
using Plots
using DataFrames, CSV, DataFramesMeta

obs = CSV.read("data/Lewis2024-NHP-PK.csv",DataFrame);

include("param.jl");
include("helper.jl");
include("jonesODEs3Monkey.jl");

# mAb IV Dose
first_dose = [1. 2. 4. 8. 16. 24.]; # [mg/kg]

BW = 6.2
V_Plasma = 0.187  # [L] 

p_cyno = deepcopy(p_N_15);
p_cyno.endothelial_scaling_factor = 60.
p_cyno.PS_Score = 4.5

tspan = (0.0, 42*DayToHour);  # hr

function cynosims(Dose_in_mgkg, BW, V_Plasma, p_cyno = p_cyno)
    u0, dose_umol = jones_init(15, Dose_in_mgkg, BW, V_Plasma);
    dosetimes = [21.]*DayToHour
    condition(u, t, integrator) = t ∈ dosetimes
    affect!(integrator) = integrator.u[1] += dose_umol/V_Plasma
    cb = DiscreteCallback(condition, affect!)
    sol = solve(ODEProblem(jonesODEs3Monkey!,u0,tspan,p_cyno), alg=AutoTsit5(Rosenbrock23(autodiff = false)), saveat = 1., callback = cb, tstops = dosetimes);
    cplasma = [sol.u[i].C_EXG_Plasma for i in 1:length(sol.t)] * MW_EDG/1E3; # [ug/mL]
    tmpdf = DataFrame(t = sol.t/DayToHour, cplasma_ugmL = cplasma); 
    tmpdf.dose .= string(Dose_in_mgkg) * " mg/kg";
    return tmpdf
end

df = DataFrame(t = [], cplasma_ugmL = [], dose = []); 

for init_dose in first_dose
    tmp = cynosims(init_dose, BW, V_Plasma, p_cyno);
    df = vcat(df, tmp)
end

p_pk = plot(xlabel = "Time (day)", ylabel = "ADC concentration (ug/mL)", legend = :outerright, palette = :Set2_6); 
plot!(df.t, df.cplasma_ugmL, group = df.dose, linewidth = 2); 
scatter!(obs.Time_day, obs.Conc_ugmL, group = obs.Dose, markerstrokewidth=0, markersize=6, alpha = 0.8);
xticks!([0, 7, 14, 21], ["0", "7", "14", "21"]);
xlabel!("time (days)");
ylabel!("plasma conc (ug/mL)");
plot!(yaxis = :log, ylims = [1, 1E3]);
plot!(title = "DHES0815a C. Macaque, endothelial scaling factor = "* string(p_cyno.endothelial_scaling_factor), titlefontsize = 6);
p_pk

savefig(p_pk, "img/DHES0815A-cyno.png");
