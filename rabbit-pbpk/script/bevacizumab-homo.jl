# date: 10/10/2025 
# author: Yuezhe Li 
# purpose of this code: to show fit of bevacizumab (Avastin) in human PK 
# PK data obtained from Shin et al., 2020; https://link.springer.com/article/10.1007/s00280-020-04144-7

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
using Statistics  # for mean()

#####################
# observed human PK; Shin et al., 2020
obs_avastin = DataFrame(
    time_hr = [11.43, 4.3, 26.61, 1.06, 3.25, 13.26, 51.27, 28.27, 50.73, 96.26, 95.85, 2.05, 12.47, 27.55, 47.62, 92.91, 164.2, 168.98, 168.63, 331.05, 338.43, 333.13, 503.25, 505.57, 502.86, 670.41, 672.74, 670.02, 1007.35, 1009.74, 1007.06, 1344.39, 1344.3, 1344.2, 1681.45, 2016.1],
    conc_ug_ml = [75.35, 84.81, 74.75, 70.81, 64.1, 61.54, 60.95, 57.4, 50.3, 47.93, 39.84, 40.43, 45.96, 43.2, 38.86, 31.95, 37.67, 31.76, 24.85, 27.42, 22.88, 18.34, 22.68, 18.54, 14.99, 18.74, 14.6, 11.05, 12.23, 9.27, 6.51, 7.69, 5.92, 3.94, 3.55, 1.97]
)

#####################
# helper functions
include(@projectroot("script/helper-funcs.jl"));
include(@projectroot("script/helper-initiation.jl")); 
include(@projectroot("script/helper-parameters.jl"));
include(@projectroot("script/helper-infusion.jl"));

# post-proressing functions
include(@projectroot("script/helper-compute-tissue-adc.jl"));
include(@projectroot("script/helper-compute-tissue-free-pl.jl"));
include(@projectroot("script/helper-cmax-cavg-adc.jl"));
include(@projectroot("script/helper-cmax-cavg-payload.jl"));

# model functions
include(@projectroot("model/mab-pbpk-modular-eye.jl"));

# constant function
include(@projectroot("model/constants.jl"));

#####################

## create pbpk model
@time pbpk = create_pbpk();
@time pbpk_simple = mtkcompile(pbpk);

## Define the simulation timespan
tspan = (0.0,504);  # hr

## initial condition (IV dosing)
u0_infusion = pbpk_initial_condition(0, pbpk_simple); 

## parameters 
param_avastin = create_base_pbpk_param(5, pbpk_simple); # from previous fitting 

## intravenous infusion, 1 hr, 
iv_3 = InfusionCallback(3, pbpk_simple, infusion_hr = 1.5, BW = 90);  # protocol from Shin et al., Cancer Chemotherapy and Pharmacology, 2020

## simulation 
@time prob_mtk_3 = ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_avastin), tspan, callback = iv_3); 

@time sol_mtk_iv_3 = solve(prob_mtk_3, alg=CVODE_BDF());

#####################

## plasma mAb visualization 
plt_plasma_mAb = plot(ylabel = "plasma ADC Concentration (uM)", ylims = [1E-3, 10],  yaxis = :log10, legend = :topright, dpi = 300, size = (400, 400)); 
plot!(sol_mtk_iv_3, idxs = pbpk_simple.plasma_exg.C_Plasma, label = "PBPK", lw = 2, alpha = 0.9);
plot!(obs_avastin.time_hr, obs_avastin.conc_ug_ml * 1E3/ MW_IGG, label = "Shin 2020", seriestype = :scatter, alpha = 0.6); 
plot!(xlabel = "Time (hr)", xlims = [0, 504], title = "Dose = 3 mg/kg, intravenous infusion over 1.5 hr", titlefontsize = 8);
display(plt_plasma_mAb)

savefig(plt_plasma_mAb, @projectroot("deliv/figure/bevacizumab-homo-plasma-mAb-iv.pdf")) 

