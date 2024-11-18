# to simulate DHES0815a human PK 

using Pkg; Pkg.activate("");

using DifferentialEquations 
using Plots
using DataFrames, CSV, DataFramesMeta

include("param.jl");
include("helper.jl");
include("jonesODEs3Homo.jl");

# observed data 
pk_lewis =  CSV.read("data/Lewis2024-clinicalPK.csv",DataFrame);

# HER2 param
p_homo = deepcopy(p_N_16);
p_homo.init_sR = 0.004

# DHES0815a params
p_homo.PS_Score = 4.5
p_homo.Kd = 0.8e-3

# dosing of DHES0815a
first_dose = [0.6, 1.2, 2.4, 4.0, 6.0]; # [mg/kg]

BW = 71.; # kg; 
V_Plasma = 3.126

tspan = (0.0, 21*DayToHour);  # hr

function adcsims(Dose_in_mgkg, BW, V_Plasma, p_homo = p_homo)
    u0, _ = jones_init(16, Dose_in_mgkg, BW, V_Plasma);
    sol = solve(ODEProblem(jonesODE3Homo!,u0,tspan,p_homo), alg=AutoTsit5(Rosenbrock23(autodiff = false)), saveat = 1.);
    cplasma = [sol.u[i].C_EXG_Plasma for i in 1:length(sol.t)] * MW_EDG; # [ug/L]
    tmpdf = DataFrame(t = sol.t/DayToHour, cplasma_ugL = cplasma); 
    tmpdf.dose .= string(Dose_in_mgkg) * " mg/kg";
    return tmpdf
end

df = DataFrame(t = [], cplasma_ugL = [], dose = []); 

for init_dose in first_dose
    tmp = adcsims(init_dose, BW, V_Plasma, p_homo);
    df = vcat(df, tmp)
end

p_pk = plot(xlabel = "Time (day)", ylabel = "ADC concentration (ug/mL)", title = "Clinical PK, DHES0815A", legend = :outerright, palette = :Set2_5); 
plot!(df.t, df.cplasma_ugL/1E3, group = df.dose, linewidth = 2); 
scatter!(pk_lewis.Time_day, pk_lewis.Conc_ngmL/1E3, group = pk_lewis.Dose, markerstrokewidth=0, markersize=6, alpha = 0.8);
plot!(yaxis = :log);
p_pk

savefig(p_pk, "img/DHES0815A-homo.png");
