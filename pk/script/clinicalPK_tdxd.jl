# date: Mar 20, 2024
# author: Yuezhe Li 
# purpose of this script: to fit clinical PK for mPBPK model 

using Pkg; Pkg.activate("../../");

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Plots, Statistics, DataFrames, CSV, DataFramesMeta

include("mPBPK_tumor_homo.jl");
include("init.jl");
include("param.jl");

TotalCells = 1E-3/Vc

# observed data 
pk_doi =  CSV.read("../../data/doi-2017.csv",DataFrame);

# update params (default params, from T-DM1)
p_base.Rcopies = 1.0E6
p_base.init_sR = 0.004
p_base.thalf_sR_adc = 120.

# fit for T-Dxd (same as T-DM1)
p_base.PS_Score = 6.

# other params for Dxd
p_base.DAR = 8.
p_base.k_in = 46.08
p_base.k_out = 32.32
p_base.k_kill_max = 0.201
p_base.ic50_pl = 9.54E-3
p_base.k_PL = 0.82

# simulation 
tspan = (0., DayToHour*66);    # [hr]
dose_mgkg = [0.8, 1.6, 3.2, 5.4, 6.4, 8.0];
df = DataFrame(Time_day = [], ADCconc_ugmL = [], Dose = []); 
for dose in dose_mgkg
    # initialization
    u0_mpbpk, adcdose_umol = pbpk_init(7, TotalCells, 0., dose*1E3, p_base.Rcopies, p_base.init_sR, p_base.k_endo, p_base.k_rec); 
    # second dosing at D21, D42
    dosetimes = [22., 44.]*DayToHour
    condition(u, t, integrator) = t ∈ dosetimes
    affect!(integrator) = integrator.u[1] += adcdose_umol/V_Plasma
    cb = DiscreteCallback(condition, affect!)
    # simulation 
    sol_mpbpk = solve(ODEProblem(mPBPK_homo_tumor!, u0_mpbpk, tspan, p_base), saveat = 3., alg = QNDF(autodiff=false), reltol = 1e-12, callback = cb, tstops = dosetimes);
    tmpdf = DataFrame(Time_day = sol_mpbpk.t/DayToHour, ADCconc_ugmL = [sol_mpbpk.u[i].C_EXG_Plasma for i in 1:length(sol_mpbpk.t)] * MW_EDG / 1E3); 
    tmpdf.Dose .= string(dose) * "mg/kg";
    df = vcat(df, tmpdf); 
end

# plot the data
p_pk = plot(xlabel = "Time (day)", ylabel = "ADC concentration (ug/mL)", title = "Clinical PK, T-Dxd", legend = :outerright, palette = :Set2_6); 
plot!(df.Time_day, df.ADCconc_ugmL, group = df.Dose, linewidth = 2); 
scatter!(pk_doi.time_day, pk_doi.T_DXd_ugperml, group = pk_doi.Dose, markerstrokewidth=0, markersize=6, alpha = 0.8);
plot!(yaxis = :log);
p_pk

savefig(p_pk, "../figure/pk/T-Dxd.png");
