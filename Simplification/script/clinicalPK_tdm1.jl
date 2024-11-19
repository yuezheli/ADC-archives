# date: Mar 20, 2024
# author: Yuezhe Li 
# purpose of this script: to fit clinical PK for mPBPK model 

using Pkg; Pkg.activate("../");

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Plots, Statistics, DataFrames, CSV, DataFramesMeta

include("mPBPK_tumor_homo.jl");
include("init.jl");
include("param.jl");

TotalCells = 1E-3/Vc

# observed data 
pk_girish =  CSV.read("../data/girish-2012.csv",DataFrame);

# update params (default params, from T-DM1)
p_base.Rcopies = 1.0E6
p_base.init_sR = 0.004
p_base.thalf_sR_adc = 120.

# fit for T-DM1
p_base.PS_Score = 6.

# simulation 
tspan = (0., DayToHour*21);    # [hr]
dose_mgkg = [0.3, 0.6, 1.2, 2.4, 3.6, 4.8];
df = DataFrame(Time_day = [], ADCconc_ugmL = [], Dose = []); 
for dose in dose_mgkg
    # initialization
    u0_mpbpk, adcdose_umol = pbpk_init(7, TotalCells, 0., dose*1E3, p_base.Rcopies, p_base.init_sR, p_base.k_endo, p_base.k_rec); 
    # simulation 
    sol_mpbpk = solve(ODEProblem(mPBPK_homo_tumor!, u0_mpbpk, tspan, p_base), saveat = 3., alg = QNDF(autodiff=false), reltol = 1e-12);
    tmpdf = DataFrame(Time_day = sol_mpbpk.t/DayToHour, ADCconc_ugmL = [sol_mpbpk.u[i].C_EXG_Plasma for i in 1:length(sol_mpbpk.t)] * MW_EDG / 1E3); 
    tmpdf.Dose .= string(dose) * "mg/kg";
    df = vcat(df, tmpdf); 
end

# plot the data
p_pk = plot(xlabel = "Time (day)", ylabel = "ADC concentration (ug/mL)", title = "Clinical PK, T-DM1", legend = :outerright, palette = :Set2_6); 
plot!(df.Time_day, df.ADCconc_ugmL, group = df.Dose, linewidth = 2); 
scatter!(pk_girish.time_day, pk_girish.T_DM1_ugperml, group = pk_girish.Dose, markerstrokewidth=0, markersize=6, alpha = 0.8);
plot!(yaxis = :log);
p_pk

savefig(p_pk, "../figure/pk/T-DM1.png");
