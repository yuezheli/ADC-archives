# date: 7/3/2025 
# author: Yuezhe Li 
# purpose of this code: compare free payload exposure between trastuzumab emtansine (T-DM1) and cantuzumab mertansine 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Parameters: @unpack
using DataFrames
using DataFramesMeta
using Plots
using Plots.PlotMeasures
using CSV

include(@projectroot("Constants.jl"))
include(@projectroot("model/jones_homo.jl") )
include(@projectroot("model/param-pk.jl"))
include(@projectroot("model/init-pk.jl")) 
include(@projectroot("script/helper-infusion-dosing.jl"))
include(@projectroot("script/helper-visualization.jl"))
include(@projectroot("script/helper-process-outcome.jl"))

# parameters for cantuzumab mertansine (from PK fit)
p_cm = deepcopy(p_base);
p_cm.PS_Score = -2      # turn off using ACSIN score to fine tune PK 
p_cm.PS_kd = 0.0001     # tuned for huC242-DM1 PK 
p_cm.Kd = 1E-2          # https://pmc.ncbi.nlm.nih.gov/articles/PMC443115/
p_cm.init_sR = 5E-5     # Cheng et al., 2011; https://pmc.ncbi.nlm.nih.gov/articles/PMC4012263/; assuming soluble MUC1 have molecular weight of 250kDa
p_cm.k_deconj = 0.13    # fitted based on conjugated ADC 
p_cm.CL_PL_plasma = log(2)/3; # https://pubmed.ncbi.nlm.nih.gov/37787918/

# parameters for trastuzumab emtansine (from PK fit)
p_tmd1 = deepcopy(p_base);
p_tmd1.init_sR = 0.004 
p_tmd1.CL_PL_plasma = log(2)/3; # https://pubmed.ncbi.nlm.nih.gov/37787918/


# simulation dose & schedule setup 
tspan = (-0.01, hr_per_day*84);      # [hr]
AddDose_q3w = [0., 21., 42., 63] * hr_per_day  # [hr]
Dose__ = 3.6  # [mg/kg]

sol_cm = InfusionDoses(Dose__, AddDose_q3w, p_cm, infusion_time = 0.5);
sol_tdm1 = InfusionDoses(Dose__, AddDose_q3w, p_tmd1, infusion_time = 0.5);

df_cm = DataFrame(
    time_hr = sol_cm.t, 
    adc_plasma_uM = [ sol_cm.u[i].C_EXG_Plasma for i in 1:length(sol_cm.t) ], 
    pl_liver_uM = [ sol_cm.u[i].end_cyto_payload[2] for i in 1:length(sol_cm.t) ]
); 

df_tdm1 = DataFrame(
    time_hr = sol_tdm1.t, 
    adc_plasma_uM = [ sol_tdm1.u[i].C_EXG_Plasma for i in 1:length(sol_tdm1.t) ], 
    pl_liver_uM = [ sol_tdm1.u[i].end_cyto_payload[2] for i in 1:length(sol_tdm1.t) ]
);

# plot PK 
plt_pk_ = plot(xlabel = "Time (Day)", ylabel = "plasma ADC (nM)", xticks = [0, 21, 42, 63, 84], background_color_legend = nothing, dpi = 300); 
plot!(df_cm.time_hr/hr_per_day, df_cm.adc_plasma_uM * 1E3, label = "Cantuzumab mertansine, 3.6 mg/kg Q3W", lw = 2, alpha = 0.7, linestyle = :dashdot); 
plot!(df_tdm1.time_hr/hr_per_day, df_tdm1.adc_plasma_uM * 1E3, label = "Trastuzumab emtansine, 3.6 mg/kg Q3W", lw = 1.5, alpha = 0.7, linestyle = :solid); 

# plot free DM1 in liver endothelial cells 
plt_pl_ = plot(xlabel = "Time (Day)", ylabel = "Free DM1 in liver endothelial cell (nM)", xticks = [0, 21, 42, 63, 84], background_color_legend = nothing, dpi = 300); 
plot!(df_cm.time_hr/hr_per_day, df_cm.pl_liver_uM * 1E3, label = "Cantuzumab mertansine, 3.6 mg/kg Q3W", lw = 2, alpha = 0.7, linestyle = :dashdot); 
plot!(df_tdm1.time_hr/hr_per_day, df_tdm1.pl_liver_uM * 1E3, label = "Trastuzumab emtansine, 3.6 mg/kg Q3W", lw = 1.5, alpha = 0.7, linestyle = :solid); 

plt_comb = plot(plt_pk_, plt_pl_, ncol = 2, size = (1000, 400), margin = 4mm)

savefig(plt_comb, @projectroot("deliv/figure/pl/dm1-based-adc-comparison.png"));
