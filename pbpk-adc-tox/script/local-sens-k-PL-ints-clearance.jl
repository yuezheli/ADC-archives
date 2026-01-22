# date: 1/14/2026 
# author: Yuezhe Li 
# purpose of this code: to conduct local sensitivity analysis on the rate of PL clearance in tissue interstitium 

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

## create pbpk model
@time pbpk = create_pbpk();
@time pbpk_simple = mtkcompile(pbpk);

## initial condition (IV dosing)
u0_infusion = pbpk_initial_condition(0, pbpk_simple); 

## parameters 
param_tdm1 = create_base_pbpk_param(6, pbpk_simple, k_PL_ints_clearance = 0.34);  #  https://pubmed.ncbi.nlm.nih.gov/37787918/
param_tdm1[pbpk_simple.DAR] = 3.5
param_tdm1[pbpk_simple.k_deconj] =  8.5E-7 * s_per_hr #  https://pubmed.ncbi.nlm.nih.gov/37787918/
param_tdm1[pbpk_simple.k_diff] = 0.14  #  https://pubmed.ncbi.nlm.nih.gov/37787918/
param_tdm1[pbpk_simple.VMAX] = 0.0012
param_tdm1[pbpk_simple.KM] = 0.006
param_tdm1[pbpk_simple.CL_plasma_PL] =  5

## lower k_PL_ints_clearance rate 
param_lowCL = create_base_pbpk_param(6, pbpk_simple, k_PL_ints_clearance = 0.034); 
param_lowCL[pbpk_simple.DAR] = 3.5
param_lowCL[pbpk_simple.k_deconj] =  8.5E-7 * s_per_hr #  https://pubmed.ncbi.nlm.nih.gov/37787918/
param_lowCL[pbpk_simple.k_diff] = 0.14  #  https://pubmed.ncbi.nlm.nih.gov/37787918/
param_lowCL[pbpk_simple.VMAX] = 0.0012
param_lowCL[pbpk_simple.KM] = 0.006
param_lowCL[pbpk_simple.CL_plasma_PL] =  5

## higher k_PL_ints_clearance rate 
param_highCL = create_base_pbpk_param(6, pbpk_simple, k_PL_ints_clearance = 3.4);  
param_highCL[pbpk_simple.DAR] = 3.5
param_highCL[pbpk_simple.k_deconj] =  8.5E-7 * s_per_hr #  https://pubmed.ncbi.nlm.nih.gov/37787918/
param_highCL[pbpk_simple.k_diff] = 0.14  #  https://pubmed.ncbi.nlm.nih.gov/37787918/
param_highCL[pbpk_simple.VMAX] = 0.0012
param_highCL[pbpk_simple.KM] = 0.006
param_highCL[pbpk_simple.CL_plasma_PL] =  5

# simulation time
tspan = (-0.01, hr_per_day*84);      # [hr]
add_dose = [0., 21., 42., 63]   

# infusion, 3.6 mg/kg
adc_infusion_3point6 = InfusionCallback(3.6, pbpk_simple, infusion_d = add_dose, BW = 90, infusion_hr = 0.5);

prob_tdm1_3point6 = ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_tdm1), tspan, callback = adc_infusion_3point6); 
sol_normCL_3point6 = solve(prob_tdm1_3point6, alg=CVODE_BDF());

prob_highCL = remake(prob_tdm1_3point6, p = param_highCL);
prob_lowCL = remake(prob_tdm1_3point6, p = param_lowCL);

sol_highCL_3point6 = solve(prob_highCL, alg=CVODE_BDF());
sol_lowCL_3point6 = solve(prob_lowCL, alg=CVODE_BDF());

#####################
# PK, plasma DM1 
plt_plasma_pl = plot(xlabel = "Time (d)", ylabel = "Plasma DM1 (uM)", size = (400, 400), dpi = 300); 
plot!(ylims = [1E-5, 1E2], yaxis = :log10, xlims = [0, 42], xticks = [0, 7, 14, 21, 28, 35, 42]); 
plot!(sol_highCL_3point6.t/hr_per_day, sol_highCL_3point6[pbpk_simple.plasma_pl.C_PL_Plasma], label = "tissue PL CL = 3.4 hr-1", lw = 2, alpha = 0.5, linestyle = :dash, color = :red); 
plot!(sol_normCL_3point6.t/hr_per_day, sol_normCL_3point6[pbpk_simple.plasma_pl.C_PL_Plasma], label = "tissue PL CL = 0.34 hr-1", lw = 2, alpha = 0.5, linestyle = :dashdot, color = :blue); 
plot!(sol_lowCL_3point6.t/hr_per_day, sol_lowCL_3point6[pbpk_simple.plasma_pl.C_PL_Plasma], label = "tissue PL CL = 0.034 hr-1", lw = 2, alpha = 0.5, color = :green); 

# DM1, liver endothelial cells 
plt_liver_endo_pl = plot(xlabel = "Time (d)", ylabel = "Liver endothelial DM1 (uM)", size = (400, 400), dpi = 300); 
plot!(ylims = [1E-5, 1E2], yaxis = :log10, xlims = [0, 42], xticks = [0, 7, 14, 21, 28, 35, 42]); 
plot!(sol_highCL_3point6.t/hr_per_day, sol_highCL_3point6[pbpk_simple.liver.PL_tissue.C_PL_endo], label = "tissue PL CL = 3.4 hr-1", lw = 2, alpha = 0.5, linestyle = :dash, color = :red); 
plot!(sol_normCL_3point6.t/hr_per_day, sol_normCL_3point6[pbpk_simple.liver.PL_tissue.C_PL_endo], label = "tissue PL CL = 0.34 hr-1", lw = 2, alpha = 0.5, linestyle = :dashdot, color = :blue); 
plot!(sol_lowCL_3point6.t/hr_per_day, sol_lowCL_3point6[pbpk_simple.liver.PL_tissue.C_PL_endo], label = "tissue PL CL = 0.034 hr-1", lw = 2, alpha = 0.5, color = :green); 

# DM1, liver interstitium
plt_liver_ints_pl = plot(xlabel = "Time (d)", ylabel = "Liver endothelial DM1 (uM)", size = (400, 400), dpi = 300); 
plot!(ylims = [1E-5, 1E2], yaxis = :log10, xlims = [0, 42], xticks = [0, 7, 14, 21, 28, 35, 42]); 
plot!(sol_highCL_3point6.t/hr_per_day, sol_highCL_3point6[pbpk_simple.liver.PL_tissue.C_PL_IntS], label = "tissue PL CL = 3.4 hr-1", lw = 2, alpha = 0.5, linestyle = :dash, color = :red); 
plot!(sol_normCL_3point6.t/hr_per_day, sol_normCL_3point6[pbpk_simple.liver.PL_tissue.C_PL_IntS], label = "tissue PL CL = 0.34 hr-1", lw = 2, alpha = 0.5, linestyle = :dashdot, color = :blue); 
plot!(sol_lowCL_3point6.t/hr_per_day, sol_lowCL_3point6[pbpk_simple.liver.PL_tissue.C_PL_IntS], label = "tissue PL CL = 0.034 hr-1", lw = 2, alpha = 0.5, color = :green); 

# save figures 
savefig(plt_plasma_pl, "deliv/figure/local-sensitivity/pl-ints-clearance-tdm1-dm1-plasma.png");
savefig(plt_liver_endo_pl, "deliv/figure/local-sensitivity/pl-ints-clearance-tdm1-dm1-liver-endothelial.png");
savefig(plt_liver_ints_pl, "deliv/figure/local-sensitivity/pl-ints-clearance-tdm1-dm1-liver-interstitium.png");

