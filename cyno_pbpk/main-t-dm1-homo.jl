# to simulate T-DM1 human data 

using DifferentialEquations 
using Plots
using DataFrames, CSV, DataFramesMeta

pk_girish =  CSV.read("data/girish-2012.csv",DataFrame);

pk_point3 = @rsubset(pk_girish, :Dose .== "0.3mg/kg" ); 
pk_point6 = @rsubset(pk_girish, :Dose .== "0.6mg/kg" ); 
pk_1point2 = @rsubset(pk_girish, :Dose .== "1.2mg/kg" ); 
pk_2point4 = @rsubset(pk_girish, :Dose .== "2.4mg/kg" ); 
pk_3point6 = @rsubset(pk_girish, :Dose .== "3.6mg/kg" ); 
pk_4point8 = @rsubset(pk_girish, :Dose .== "4.8mg/kg" ); 

include("param.jl");
include("helper.jl");
include("jonesODEs3Homo.jl");

p_homo = deepcopy(p_N_16);
p_homo.init_sR = 0.004

first_dose = [0.3, 0.6, 1.2, 2.4, 3.6, 4.8, 1.8]; # [mg/kg]

BW = 71.; # kg; 
V_Plasma = 3.126

tspan = (0.0, 21*DayToHour);  # hr

function tdm1sims(Dose_in_mgkg, BW, V_Plasma, p_homo = p_homo)
    u0, _ = jones_init(16, Dose_in_mgkg, BW, V_Plasma);
    sol = solve(ODEProblem(jonesODE3Homo!,u0,tspan,p_homo), alg=AutoTsit5(Rosenbrock23(autodiff = false)), saveat = 1.);
    cplasma = [sol.u[i].C_EXG_Plasma for i in 1:length(sol.t)] * MW_EDG; # [ug/L]
    return sol.t/DayToHour, cplasma
end

timestamp = []; 
cplasma = [];

for init_dose in first_dose
    tmpt, tmpp = tdm1sims(init_dose, BW, V_Plasma, p_homo);
    append!(timestamp, [tmpt]);
    append!(cplasma, [tmpp]); 
end

p_pk_tdm1 = plot(legend = :topright, ylims = [1E2, 1E6], yaxis = :log, size=(350,350), legendcolumns=3, dpi = 1000);
plot!(timestamp[1], cplasma[1], lw = 2, alpha = 0.6, color = "red", label = "0.3mg/kg");
scatter!(pk_point3.time_day, pk_point3.T_DM1_ugperml*1E3, ma = 0.6, color = "red", label = false, linewidth=0);
plot!(timestamp[2], cplasma[2], lw = 2, alpha = 0.6, color = "seagreen", label = "0.6mg/kg");
scatter!(pk_point6.time_day, pk_point6.T_DM1_ugperml*1E3, ma = 0.6, color = "seagreen", label = false, linewidth=0);
plot!(timestamp[3], cplasma[3], lw = 2, alpha = 0.6, color = "mediumpurple", label = "1.2mg/kg");
scatter!(pk_1point2.time_day, pk_1point2.T_DM1_ugperml*1E3, ma = 0.6, color = "mediumpurple", label = false, linewidth=0);
plot!(timestamp[4], cplasma[4], lw = 2, alpha = 0.6, color = "violetred", label = "2.4mg/kg");
scatter!(pk_2point4.time_day, pk_2point4.T_DM1_ugperml*1E3, ma = 0.6, color = "violetred", label = false, linewidth=0);
plot!(timestamp[5], cplasma[5], lw = 2, alpha = 0.6, color = "blue2", label = "3.6mg/kg");
scatter!(pk_3point6.time_day, pk_3point6.T_DM1_ugperml*1E3, ma = 0.6, color = "blue2", label = false, linewidth=0);
plot!(timestamp[6], cplasma[6], lw = 2, alpha = 0.6, color = "teal", label = "4.8mg/kg");
scatter!(pk_4point8.time_day, pk_4point8.T_DM1_ugperml*1E3, ma = 0.6, color = "teal", label = false, linewidth=0);
xlabel!("time (day)"); ylabel!("plasma T-DM1 conc (ug/L)"); xlims!(0, 21)

savefig(p_pk_tdm1, "img/t-dm1-homo.png");
