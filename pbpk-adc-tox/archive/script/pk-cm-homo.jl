# date: 1/27/2025
# author: Yuezhe Li 
# purpose of this code: to generate simulation of Cantuzumab mertansine at 235 mg/m2 for PK validation

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using DataFrames, CSV, DataFramesMeta
using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Parameters: @unpack
using Plots, Measures

# read observed; https://pubmed.ncbi.nlm.nih.gov/12525512/; The dose-limiting toxicity (DLT) of cantuzumab mertansine was found to be reversible elevations of hepatic transaminases
pk_tolcher =  CSV.read("data/Tolcher-2003.csv",DataFrame); # total ADC, 235 mg/m2

# ADC total & conjugated comparison (dose unknown)
pk_tolcher_total_adc = DataFrame(
    time_d = [0.657, 0.654, 0.749, 1.268, 2.5, 4.5, 14.19], 
    conc_ugmL = [140.452, 130.744, 101.866, 75.238, 61.886, 41.863, 11.885]
);

pk_tolcher_conj_adc = DataFrame(
    time_d = [0.766, 2.3, 4.272, 21], 
    conc_ugmL = [183.6, 36.22, 11.372, 0.031]
);

# simulation 
include(@projectroot("Constants.jl"))
include(@projectroot("model/jones_homo.jl") )
include(@projectroot("model/param-pk.jl"))
include(@projectroot("model/init-pk.jl")) 
include(@projectroot("script/helper-infusion-dosing.jl"))

# DAR value not changed based on https://pubmed.ncbi.nlm.nih.gov/18301896/
# use total ADC PK for fitting PS_kd, then use conjugated ADC 
p_cm = deepcopy(p_base);
p_cm.PS_Score = -2      # turn off using ACSIN score to fine tune PK 
p_cm.PS_kd = 0.0001     # tuned for huC242-DM1 PK 
p_cm.Kd = 1E-2          # https://pmc.ncbi.nlm.nih.gov/articles/PMC443115/
p_cm.init_sR = 5E-5     # Cheng et al., 2011; https://pmc.ncbi.nlm.nih.gov/articles/PMC4012263/; assuming soluble MUC1 have molecular weight of 250kDa
p_cm.k_deconj = 0.13    # fitted based on conjugated ADC 
p_cm.CL_PL_plasma = log(2)/3; # https://pubmed.ncbi.nlm.nih.gov/37787918/

p_cm_totaladc = deepcopy(p_cm);
p_cm_totaladc.k_deconj = 0

tspan = (-0.01, hr_per_day*84);      # [hr]
AddDose_q3w = [0., 21., 42., 63] * hr_per_day  # [hr]

init_dose = 235;  # [mg/m2]

sol_pk_tolcher_2003 = InfusionDoses(init_dose*HT*HT/BW, AddDose_q3w, p_cm, infusion_time = 0.5);
sol_pk_tolcher_totaladc = InfusionDoses(init_dose*HT*HT/BW, AddDose_q3w, p_cm_totaladc, infusion_time = 0.5);

p_pk_cm = plot(legend = :topright, ylims = [1E-4, 2], xlims = [0, 500], yaxis = :log, size=(400,400), dpi = 300);
plot!(sol_pk_tolcher_2003.t, [sol_pk_tolcher_2003.u[i].C_EXG_Plasma for i in 1:length(sol_pk_tolcher_2003.t)], linewidth = 2, alpha = 0.8, color = "red", label = "235mg/m² Q3W"); 
# scatter!(pk_tolcher.time_hour, pk_tolcher.conc_ug_L, ma = 0.6, label = false, linewidth=0);
scatter!(pk_tolcher_conj_adc.time_d * 7, pk_tolcher_conj_adc.conc_ugmL / MW_IGG * 1E3, ma = 0.6, label = false, linewidth=0);
xlabel!("Time (hour)"); ylabel!("Plasma conjugated cantuzumab mertansine (uM)", guidefontsize = 10); 

p_totalpk_cm = plot(legend = :topright, ylims = [1E-4, 2], xlims = [0, 500], yaxis = :log, size=(400,400), dpi = 300);
plot!(sol_pk_tolcher_totaladc.t, [sol_pk_tolcher_totaladc.u[i].C_EXG_Plasma for i in 1:length(sol_pk_tolcher_totaladc.t)], linewidth = 2, alpha = 0.8, color = "red", label = "235mg/m² Q3W"); 
scatter!(pk_tolcher.time_hour, pk_tolcher.conc_ug_L / MW_IGG, ma = 0.6, label = false, linewidth=0);
xlabel!("Time (hour)"); ylabel!("Plasma total cantuzumab mertansine (uM)", guidefontsize = 10); 

savefig(plot(p_pk_cm, p_totalpk_cm, ncol = 2, size = (900, 400), margin = 5mm), @projectroot("deliv/figure/pk/cantuzumab-mertansine-homo.png"));

# plot DM1 (295 mg/m2)
sol_dm1_tolcher_2003 = InfusionDoses(295*HT*HT/BW, AddDose_q3w, p_cm, infusion_time = 0.5);
plasma_dm1 = [sol_dm1_tolcher_2003.u[i].plasma_payload for i in 1:length(sol_dm1_tolcher_2003.t)]; # [uM]

plt_dm1 = plot(xlabel = "Time (d)", ylabel = "DM1 (uM)", dpi = 300, size = (400,400), background_color_legend = nothing);
plot!(sol_dm1_tolcher_2003.t , plasma_dm1, label = "sims, 295mg/m² Q3W");
plot!(xlims = [0, 98], xticks = [0, 24, 48, 72, 96], ylims = [1E-3, 1], yaxis = :log);

savefig(plt_dm1, @projectroot("deliv/figure/pk/cantuzumab-mertansine-homo-dm1.png"));
