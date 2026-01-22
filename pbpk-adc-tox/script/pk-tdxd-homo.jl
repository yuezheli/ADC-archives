# date: 1/13/2026 
# author: Yuezhe Li 
# purpose of this code: PK fitting of T-Dxd

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

# T-Dxd PK; Doi et al., 2017;  https://pubmed.ncbi.nlm.nih.gov/29037983/

obs_tdxd_pk = DataFrame(
    time_day = [0.16264, 0.55073, 2.60205, 7.88166, 15.20759, 0.16659, 0.65369, 2.70542, 7.88998, 15.22548, 22.46135, 22.50129, 44.50404, 44.63985, 51.77174, 59.10745, 66.44336, 0.75934, 3.00762, 8.09463, 15.33612, 22.57491, 22.60819, 44.42501, 51.88842, 58.93358, 66.46895, 0.66429, 2.71686, 7.9035, 15.14603, 22.68161, 22.61027, 44.53108, 44.55167, 51.79233, 59.1322, 66.47207, 0.57049, 7.81011, 15.24919, 22.49276, 22.61297, 44.53462, 44.65358, 51.69832, 59.13636, 66.4781, 0.57028, 2.72206, 7.90807, 15.24877, 22.59009, 22.51585, 44.53815, 44.65504, 51.79857, 59.13844, 66.48267],
    tdxd_ugmL = [20.02102, 12.89447, 6.50374, 1.96333, 0.16627, 31.85552, 23.18358, 12.27925, 5.21942, 1.36069, 0.29172, 31.85552, 0.3635, 31.0863, 3.28038, 0.87635, 0.23991, 57.27456, 32.64378, 14.57075, 6.04389, 1.82451, 91.12977, 3.36155, 29.60313, 11.69339, 4.85037, 80.6458, 47.10165, 25.56485, 11.98274, 5.09338, 116.36367, 8.72071, 98.06351, 36.88748, 16.06737, 6.99859, 131.49098, 43.77125, 22.07746, 11.69339, 159.89013, 13.21354, 156.02922, 58.69181, 26.19745, 14.21891, 128.31584, 86.78184, 43.77125, 21.02411, 10.86659, 176.31307, 20.02102, 185.14668, 76.79807, 33.45154, 24.34512],
    dose = ["0.8mg/kg", "0.8mg/kg", "0.8mg/kg", "0.8mg/kg", "0.8mg/kg", "1.6mg/kg", "1.6mg/kg", "1.6mg/kg", "1.6mg/kg", "1.6mg/kg", "1.6mg/kg", "1.6mg/kg", "1.6mg/kg", "1.6mg/kg", "1.6mg/kg", "1.6mg/kg", "1.6mg/kg", "3.2mg/kg", "3.2mg/kg", "3.2mg/kg", "3.2mg/kg", "3.2mg/kg", "3.2mg/kg", "3.2mg/kg", "3.2mg/kg", "3.2mg/kg", "3.2mg/kg", "5.4mg/kg", "5.4mg/kg", "5.4mg/kg", "5.4mg/kg", "5.4mg/kg", "5.4mg/kg", "5.4mg/kg", "5.4mg/kg", "5.4mg/kg", "5.4mg/kg", "5.4mg/kg", "6.4mg/kg", "6.4mg/kg", "6.4mg/kg", "6.4mg/kg", "6.4mg/kg", "6.4mg/kg", "6.4mg/kg", "6.4mg/kg", "6.4mg/kg", "6.4mg/kg", "8.0mg/kg", "8.0mg/kg", "8.0mg/kg", "8.0mg/kg", "8.0mg/kg", "8.0mg/kg", "8.0mg/kg", "8.0mg/kg", "8.0mg/kg", "8.0mg/kg", "8.0mg/kg"]
); 

# convert observed data to dictionary 
pk_tdxd = Dict();
for dose_str in unique(obs_tdxd_pk.dose)
    tmp_pk_ = @rsubset(obs_tdxd_pk, :dose == dose_str);
    dose_num = parse(Float64, replace(dose_str, "mg/kg"=>""));
    tmp_df = DataFrame(
        time_d = tmp_pk_.time_day, 
        ADC_uM = tmp_pk_.tdxd_ugmL*1E3/MW_TDXD
    )
    push!(pk_tdxd, dose_num => tmp_df)
end

obs_tdxd_pk_6point4 = @rsubset(obs_tdxd_pk, :dose .== "6.4mg/kg")

# Total antibody PK (6.4 mg/kg)
obs_totalAb_pk = DataFrame(
    time_d = [0.0, 0.77, 2.75, 7.8, 21.75, 21.84, 42.93, 43.03, 49.95, 56.87, 63.9], 
    tAb_ugL = [124262.37, 109077.67, 73777.6, 45748.55, 13554.25, 129779.93, 19187.52, 135542.49, 64762.07, 36816.1, 20039.5]
    ); 
@transform!(obs_totalAb_pk, :tAb_uM = :tAb_ugL/MW_trastuzumab); 


# Dxd PK
obs_dxd_pk = DataFrame(
    time_d = [0.058, 0.2, 0.893, 3, 8, 15, 21.97, 22, 43, 44, 50, 57, 64.2], 
    dxd_ugL = [3.186, 6.56, 4.854, 2.996, 1.542, 0.682, 0.361, 2.2, 0.5, 5.713, 2.24, 0.962, 0.51]
    ); 
@transform!(obs_dxd_pk, :dxd_uM = :dxd_ugL/MW_Dxd); 


#####################

## create pbpk model
@time pbpk = create_pbpk();
@time pbpk_simple = mtkcompile(pbpk);

## initial condition (IV dosing)
u0_infusion = pbpk_initial_condition(0, pbpk_simple); 

## parameters 
### k_PL_ints_clearance from https://pubmed.ncbi.nlm.nih.gov/37787918/
### AC-SINS score tuned based on total Ab PK 
### VMAX and KM tuned based on TMDD across lower doses 
### k_deconj tuned based on Dxd & ADC PK 
### k_diff calculated based on Wood et al., 2025; https://pubmed.ncbi.nlm.nih.gov/40802833/; Peff from https://pmc.ncbi.nlm.nih.gov/articles/PMC8000490/
param_tdxd = create_base_pbpk_param(5, pbpk_simple, k_PL_ints_clearance = 0.34);
param_tdxd[pbpk_simple.DAR] = 8
param_tdxd[pbpk_simple.k_deconj] =  2E-3
param_tdxd[pbpk_simple.k_diff] = 3.0E-6 * (4*pi*cell_radius^2)/(4/3*pi*cell_radius^3) * s_per_hr #  Peff from https://pmc.ncbi.nlm.nih.gov/articles/PMC8000490/
param_tdxd[pbpk_simple.CL_plasma_PL] = 3
param_tdxd[pbpk_simple.VMAX] = 0.0012
param_tdxd[pbpk_simple.KM] = 0.006
param_tdxd[pbpk_simple.Kp_adc_cor_ah] = 1

# simulation time
tspan = (0, hr_per_day*84);      # [hr]
add_dose_doi = [0., 22., 44.]; 

# simulation dose 
sims_dose = [0.8, 1.6, 3.2, 5.4, 6.4, 8.0]; # [mg/kg]

# simulation of T-Dxd
sol_pk_tdxd = Dict(); 
for init_dose in sims_dose
    adc_infusion_ = InfusionCallback(init_dose, pbpk_simple, infusion_d = add_dose_doi, BW = 90, infusion_hr = 0.5);
    # simulation 
    @time prob_mtk_ = ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_tdxd), tspan, callback = adc_infusion_); 
    @time sol_mtk_infusion_ = solve(prob_mtk_, alg=CVODE_BDF());
    # append result 
    push!(sol_pk_tdxd, init_dose => sol_mtk_infusion_);
end

# simulation of total Ab 
param_nodeconj = deepcopy(param_tdxd); 
param_nodeconj[pbpk_simple.k_deconj] =  0
adc_infusion_6point4 = InfusionCallback(6.4, pbpk_simple, infusion_d = add_dose_doi, BW = 90, infusion_hr = 0.5);
prob_tab_6point4 = ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_nodeconj), tspan, callback = adc_infusion_6point4); 
@time sol_tab_6point4 = solve(prob_tab_6point4, alg=CVODE_BDF());

#####################
# PK, ADC
plt_plasma_adc = PlotSimulationPlasma(sol_pk_tdxd, pk_tdxd, pbpk_simple, adc_name = "trastuzumab deruxtecan", colorPALETTE = :batlowKS, xrange = [0, 63], yrange = [1E-4, 1E3], ylog = true, legendcolumnsnum = 2); 

# PK, total Ab, ADC, PL (6.4 mg/kg)
plt_plasma_tAb_dxd = plot(xlabel = "Time (d)", ylabel = "Concentration (uM)", size = (400, 400), dpi = 300); 
plot!(ylims = [1E-4, 1E3], yaxis = :log10, xlims = [0, 42], xticks = [0, 7, 14, 21, 28, 35, 42]);
plot!(sol_tab_6point4.t/hr_per_day, sol_tab_6point4[pbpk_simple.plasma_exg.C_Plasma], label = "total Ab", color = palette(:batlowKS)[2], lw = 2, linestyle = :dot); 
scatter!(obs_totalAb_pk.time_d, obs_totalAb_pk.tAb_uM, label = "Doi et al., 2017; total Ab", markershape=:diamond, color = palette(:batlowKS)[2], alpha = 0.6);
plot!(sol_pk_tdxd[6.4].t/hr_per_day, sol_pk_tdxd[6.4][pbpk_simple.plasma_pl.C_PL_Plasma], label = "Dxd", color = palette(:batlowKS)[2], lw = 2); 
scatter!(obs_dxd_pk.time_d, obs_dxd_pk.dxd_uM, label = "Doi et al., 2017; Dxd", markershape=:star6, color = palette(:batlowKS)[2], alpha = 0.6); 

display(plot(plt_plasma_adc, plt_plasma_tAb_dxd, layout = (1,2), size = (800, 400), margin = 5mm))

# save pk figures  
savefig(plt_plasma_adc, "deliv/figure/pk-tdxd-adc-plasma.png");
savefig(plt_plasma_tAb_dxd, "deliv/figure/pk-tdxd-tab-dxd-plasma.png");

#####################
# visualization, liver PL 
plt_liver_pl = plot(xlabel = "Time (d)", ylabel = "Dxd (uM)", 
                     ylims = [1E-4, 1E2], yaxis = :log10, xlims = [0, 42], xticks = [0, 7, 14, 21, 28, 35, 42],
                     size = (400, 400), dpi = 300); 
plot!(sol_pk_tdxd[8.0].t/hr_per_day, sol_pk_tdxd[8.0][pbpk_simple.liver.PL_tissue.C_PL_IntS], label = "liver, interstitium, 8.0 mg/kg", linestyle = :dashdot, color = :gray39); 
plot!(sol_pk_tdxd[8.0].t/hr_per_day, sol_pk_tdxd[8.0][pbpk_simple.liver.PL_tissue.C_PL_endo], label = "liver, endothelial cells, 8.0 mg/kg", linestyle = :solid, color = :gray39);
plot!(sol_pk_tdxd[5.4].t/hr_per_day, sol_pk_tdxd[5.4][pbpk_simple.liver.PL_tissue.C_PL_IntS], label = "liver, interstitium, 5.4 mg/kg", linestyle = :dashdot, color = :turquoise4); 
plot!(sol_pk_tdxd[5.4].t/hr_per_day, sol_pk_tdxd[5.4][pbpk_simple.liver.PL_tissue.C_PL_endo], label = "liver, endothelial cells, 5.4 mg/kg", linestyle = :solid, color = :turquoise4); 
plot!(sol_pk_tdxd[3.2].t/hr_per_day, sol_pk_tdxd[3.2][pbpk_simple.liver.PL_tissue.C_PL_IntS], label = "liver, interstitium, 3.2 mg/kg", linestyle = :dashdot, color = :orangered1); 
plot!(sol_pk_tdxd[3.2].t/hr_per_day, sol_pk_tdxd[3.2][pbpk_simple.liver.PL_tissue.C_PL_endo], label = "liver, endothelial cells, 3.2 mg/kg", linestyle = :solid, color = :orangered1); 

display(plt_liver_pl)

savefig(plt_liver_pl, "deliv/figure/pk-tdxd-dxd-liver.png");
