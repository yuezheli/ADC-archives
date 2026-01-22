# date: 1/14/2026 
# author: Yuezhe Li 
# purpose of this code: to fit for PK of Cantuzumab mertansine

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using DifferentialEquations, ModelingToolkit
using Plots
using Plots.Measures
using DataFrames
using DataFramesMeta
using SymbolicIndexingInterface  # for parameter_index()
using CSV
using Sundials  # for CVODE_BDF()

#####################

# constant function
include(@projectroot("model/constants.jl"));

# model functions
include(@projectroot("model/mab-pbpk-modular-eye.jl"));
include(@projectroot("model/default-parameters.jl"));

# helper functions
include(@projectroot("script/helper-initiation.jl")); 
include(@projectroot("script/helper-infusion.jl"));
include(@projectroot("script/helper-visualization.jl"));

#####################

# read observed; https://pubmed.ncbi.nlm.nih.gov/12525512/; The dose-limiting toxicity (DLT) of cantuzumab mertansine was found to be reversible elevations of hepatic transaminases
pk_tolcher_total_adc = DataFrame(
    time_d = [0.657, 0.654, 0.749, 1.268, 2.5, 4.5, 14.19, 0.08, 0, 0.57, 1.35, 2.3, 3.86, 21, 0, 21, 0, 21, 0, 21], 
    conc_ugmL = [140.452, 130.744, 101.866, 75.238, 61.886, 41.863, 11.885, 134.99, 234.78, 205.64, 142.07, 111.95, 60.95, 15.16, 118.33, 42.36, 200.37, 10.21, 282.38, 24.36]
);
@transform!(pk_tolcher_total_adc, :tAb_uM = :conc_ugmL*1E3/MW_IGG); 

pk_tolcher_conj_adc = DataFrame(
    time_d = [0.766, 2.3, 4.272, 21, 0, 1.16, 2.08, 3.95, 21, 0, 21, 0, 21, 0, 21, 0.4, 0.25, 1, 2, 4, 21], 
    conc_ugmL = [183.6, 36.22, 11.372, 0.031, 103.75, 58.02, 27.72, 9.9, 0.04, 260.87, 0.03, 247.43, 0.06, 222.73, 0.07, 158.04, 101.1, 58.04, 28.45, 10.43, 0.04]
);
@transform!(pk_tolcher_conj_adc, :ADC_uM = :conc_ugmL*1E3/MW_IGG); 

#####################
## create pbpk model
@time pbpk = create_pbpk();
@time pbpk_simple = mtkcompile(pbpk);

## initial condition (IV dosing)
u0_infusion = pbpk_initial_condition(0, pbpk_simple); 

## simulation time
tspan = (-0.01, hr_per_day*84);      # [hr]
add_dose = [0., 21., 42., 63]   

## parameters 
param_cm = create_base_pbpk_param(8, pbpk_simple, k_PL_ints_clearance = 0.34);  #  https://pubmed.ncbi.nlm.nih.gov/37787918/
param_cm[pbpk_simple.DAR] = 3.5
param_cm[pbpk_simple.k_diff] = 0.14  # as T-DM1, https://pubmed.ncbi.nlm.nih.gov/37787918/
param_cm[pbpk_simple.CL_plasma_PL] =  log(2)/3; # https://pubmed.ncbi.nlm.nih.gov/37787918/
param_cm[pbpk_simple.k_deconj] =  0.13
param_cm[pbpk_simple.VMAX] = 0

param_nodeconj = deepcopy(param_cm);
param_nodeconj[pbpk_simple.k_deconj] = 0

# infusion, 235 mg/m2
adc_infusion_235 = InfusionCallback(235*3/human_WT, pbpk_simple, infusion_d = add_dose, infusion_hr = 0.5);

# simulation of total Ab 
prob_tab_235 = ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_nodeconj), tspan, callback = adc_infusion_235); 
@time sol_tab_235 = solve(prob_tab_235, alg=CVODE_BDF());

# simulation of ADC
prob_adc_235 = ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_cm), tspan, callback = adc_infusion_235); 
@time sol_adc_235 = solve(prob_adc_235, alg=CVODE_BDF());

#####################
# PK, ADC
plt_plasma_adc = plot(xlabel = "Time (d)", ylabel = "ADC (uM)", size = (400, 400), dpi = 300); 
plot!(ylims = [1E-4, 100], yaxis = :log10, xlims = [0, 21], xticks = [0, 7, 14, 21]);
plot!(sol_adc_235.t/hr_per_day, sol_adc_235[pbpk_simple.plasma_exg.C_Plasma], label = "sims", lw = 2); 
plot!(pk_tolcher_conj_adc.time_d, pk_tolcher_conj_adc.ADC_uM, label = "Tolcher et al., 2003", seriestype = :scatter, alpha = 0.6); 

# PK, total Ab 
plt_plasma_tAb = plot(xlabel = "Time (d)", ylabel = "total Ab (uM)", size = (400, 400), dpi = 300); 
plot!(ylims = [1E-4, 100], yaxis = :log10, xlims = [0, 21], xticks = [0, 7, 14, 21]);
plot!(sol_tab_235.t/hr_per_day, sol_tab_235[pbpk_simple.plasma_exg.C_Plasma], label = "sims", lw = 2); 
plot!(pk_tolcher_total_adc.time_d, pk_tolcher_total_adc.tAb_uM, label = "Tolcher et al., 2003", seriestype = :scatter, alpha = 0.6); 

display(plot(plt_plasma_adc, plt_plasma_tAb, layout = (1,2), size = (800, 400)))

# save pk figures  
savefig(plt_plasma_adc, "deliv/figure/pk-cm-adc-plasma.png");
savefig(plt_plasma_tAb, "deliv/figure/pk-cm-tab-plasma.png");

# PK, PL 
plt_plasma_pl = plot(xlabel = "Time (d)", ylabel = "DM1 (uM)", 
                     ylims = [1E-4, 1E2], yaxis = :log10, xlims = [0, 42], xticks = [0, 7, 14, 21, 28, 35, 42],
                     size = (400, 400), dpi = 300); 
plot!(sol_tab_235.t/hr_per_day, sol_tab_235[pbpk_simple.plasma_pl.C_PL_Plasma], label = "sims",  lw = 2); 

savefig(plt_plasma_pl, "deliv/figure/pk-cm-pl-plasma.png");

#####################
# visualization, liver PL 
plt_liver_pl = plot(xlabel = "Time (d)", ylabel = "DM1 (uM)", 
                     ylims = [1E-4, 1E2], yaxis = :log10, xlims = [0, 42], xticks = [0, 7, 14, 21, 28, 35, 42],
                     size = (400, 400), dpi = 300); 
plot!(sol_tab_235.t/hr_per_day, sol_tab_235[pbpk_simple.liver.PL_tissue.C_PL_IntS], label = "liver, interstitium, 235 mg/m2",); 
plot!(sol_tab_235.t/hr_per_day, sol_tab_235[pbpk_simple.liver.PL_tissue.C_PL_endo], label = "liver, endothelial cells, 235 mg/m2");

savefig(plt_liver_pl, "deliv/figure/pk-cm-pl-liver.png");

