# date: Mar 20, 2024
# author: Yuezhe Li 
# purpose of this script: to fit clinical PK for mPBPK model 

using Pkg; Pkg.activate("..");

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Plots, Statistics, DataFrames, CSV, DataFramesMeta

include("mPBPK_tumor_homo.jl");
include("init.jl");
include("param.jl");

TotalCells = 1E-3/Vc

# observed data 
pk_lewis =  CSV.read("../data/Lewis2024-clinicalPK.csv",DataFrame);

# update params (default params, from T-DM1)
p_base.Rcopies = 1.0E6
p_base.init_sR = 0.004
p_base.thalf_sR_adc = 120.

# fit for DHES0815A
p_base.PS_Score = 4.5

# update params (DHES0815a specific)
p_base.ic50_pl = 21.7E-3
p_base.DAR = 2.
p_base.Kd = 0.8e-3

# update params (from in vitro optimization)
p_base.k_kill_max = 0.036; 
p_base.k_out = 8.24;
p_base.k_PL = 1.;
p_base.tau = 2.;
p_base.k_in = 5.

# simulation 
tspan = (0., DayToHour*21);    # [hr]
dose_mgkg = [0.6, 1.2, 2.4, 4.0, 6.0];
df = DataFrame(Time_day = [], ADCconc_ugmL = [], ADC_lung = [], ADC_liver =[], ADC_skin = [], 
               ints_pl_lung = [], ints_pl_liver = [], ints_pl_skin = [], end_cyto_pl_lung = [], end_cyto_pl_liver = [], end_cyto_pl_skin = [], Dose = []); 
for dose in dose_mgkg
    # initialization
    u0_mpbpk, adcdose_umol = pbpk_init(7, TotalCells, 0., dose*1E3, p_base.Rcopies, p_base.init_sR, p_base.k_endo, p_base.k_rec); 
    # simulation 
    sol_mpbpk = solve(ODEProblem(mPBPK_homo_tumor!, u0_mpbpk, tspan, p_base), saveat = 3., alg = QNDF(autodiff=false), reltol = 1e-12);
    tmpdf = DataFrame(Time_day = sol_mpbpk.t/DayToHour, 
                      ADCconc_ugmL = [sol_mpbpk.u[i].C_EXG_Plasma for i in 1:length(sol_mpbpk.t)] * MW_EDG / 1E3, 
                      ADC_lung = [sol_mpbpk.u[i].C_EXG[1] for i in 1:length(sol_mpbpk.t)] * MW_EDG / 1E3, 
                      ADC_liver = [sol_mpbpk.u[i].C_EXG[2] for i in 1:length(sol_mpbpk.t)] * MW_EDG / 1E3, 
                      ADC_skin = [sol_mpbpk.u[i].C_EXG[3] for i in 1:length(sol_mpbpk.t)] * MW_EDG / 1E3, 
                      ints_pl_lung = [sol_mpbpk.u[i].ints_payload[1] for i in 1:length(sol_mpbpk.t)], 
                      ints_pl_liver = [sol_mpbpk.u[i].ints_payload[2] for i in 1:length(sol_mpbpk.t)], 
                      ints_pl_skin = [sol_mpbpk.u[i].ints_payload[3] for i in 1:length(sol_mpbpk.t)], 
                      end_cyto_pl_lung = [sol_mpbpk.u[i].end_cyto_payload[1] for i in 1:length(sol_mpbpk.t)], 
                      end_cyto_pl_liver = [sol_mpbpk.u[i].end_cyto_payload[2] for i in 1:length(sol_mpbpk.t)], 
                      end_cyto_pl_skin = [sol_mpbpk.u[i].end_cyto_payload[3] for i in 1:length(sol_mpbpk.t)] ); 
    tmpdf.Dose .= string(dose) * " mg/kg";
    df = vcat(df, tmpdf); 
end

# plot plasma PK
p_pk = plot(xlabel = "Time (day)", ylabel = "ADC concentration (ug/mL)", title = "Clinical PK, DHES0815A", legend = :outerright, palette = :Set2_5); 
plot!(df.Time_day, df.ADCconc_ugmL, group = df.Dose, linewidth = 2); 
scatter!(pk_lewis.Time_day, pk_lewis.Conc_ngmL/1E3, group = pk_lewis.Dose, markerstrokewidth=0, markersize=6, alpha = 0.8);
plot!(yaxis = :log);
p_pk

savefig(p_pk, "../figure/pk/DHES0815A.png");

# plot ADC concentration in tissue interstitium 
p_lung_ints_adc = plot(xlabel = "Time (day)", ylabel = "ADC concentration (ug/mL)", title = "Lung interstitium ADC", legend = :topright, palette = :Set2_5); 
plot!(df.Time_day, df.ADC_lung, group = df.Dose, linewidth = 2); 

p_skin_ints_adc = plot(xlabel = "Time (day)", ylabel = "ADC concentration (ug/mL)", title = "Skin interstitium ADC", legend = :topright, palette = :Set2_5); 
plot!(df.Time_day, df.ADC_skin, group = df.Dose, linewidth = 2); 

plot(p_lung_ints_adc, p_skin_ints_adc, ncol = 2)

# plot tissue interstitium payload conc 
p_ints_payload = plot(xlabel = "Time (day)", ylabel = "Payload concentration (uM)", title = "Skin interstitial payload conc", legend = :outerright, palette = :Set2_5); 
plot!(df.Time_day, df.end_cyto_pl_skin, group = df.Dose, linewidth = 2); 
hline!([p_base.ic50_pl], linestyle = :dash, color = :black, label = false);
plot!(ylims = [1E-6, 1E-1], yaxis = :log);
p_ints_payload
