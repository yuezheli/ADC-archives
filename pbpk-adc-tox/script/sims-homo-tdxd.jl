# date: 1/15/2026
# author: Yuezhe Li 
# purpose of this code: simulation of T-Dxd (for further analysis in R)
# additional dose was added (DAR-adj of T-DM1 dose)

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
include(@projectroot("script/helper-extract-values.jl"));

#####################

## create pbpk model
@time pbpk = create_pbpk();
@time pbpk_simple = mtkcompile(pbpk);

## initial condition (IV dosing)
u0_infusion = pbpk_initial_condition(0, pbpk_simple); 


## parameters 
param_tdxd = create_base_pbpk_param(5, pbpk_simple, k_PL_ints_clearance = 0.34);
param_tdxd[pbpk_simple.DAR] = 8
param_tdxd[pbpk_simple.k_deconj] =  2E-3
param_tdxd[pbpk_simple.k_diff] = 3.0E-6 * (4*pi*cell_radius^2)/(4/3*pi*cell_radius^3) * s_per_hr #  Peff from https://pmc.ncbi.nlm.nih.gov/articles/PMC8000490/
param_tdxd[pbpk_simple.CL_plasma_PL] = 3
param_tdxd[pbpk_simple.VMAX] = 0.0012
param_tdxd[pbpk_simple.KM] = 0.006
param_tdxd[pbpk_simple.Kp_adc_cor_ah] = 1

# simulation time
tspan = (-0.01, hr_per_day*84);      # [hr]
add_dose = [0., 21., 42., 63]   

# simulation dose 
tmd1_sims_dose = [0.3, 0.6, 1.2, 2.4, 3.6, 4.8]; # [mg/kg]
tdxd_tdm1_dar_adj = tmd1_sims_dose * 3.5/8

tdxd_sims_dose = [0.8, 1.6, 3.2, 5.4, 6.4]; # [mg/kg]; highest dose (8.0 mg/kg) was eliminated 

# simulation, 8.0 mg/kg
adc_8 = InfusionCallback(8, pbpk_simple, infusion_d = add_dose);
sol_mtk_8 = solve(ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_tdxd), tspan, callback = adc_8, infusion_hr = 0.1), alg=CVODE_BDF());
df_adc_ints = extract_adc_tissue_interstitia(sol_mtk_8, pbpk_simple, 8)
df_pl_ints = extract_PL_tissue_interstitia(sol_mtk_8, pbpk_simple, 8)
df_pl_endo = extract_PL_tissue_endothelial_cells(sol_mtk_8, pbpk_simple, 8)

# simulation, lower T-Dxd doses
for tmp_dose in tdxd_sims_dose
    adc__ = InfusionCallback(tmp_dose, pbpk_simple, infusion_d = add_dose);
    sol_mtk__ = solve(ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_tdxd), tspan, callback = adc__, infusion_hr = 0.1), alg=CVODE_BDF());
    tmpdf_adc_ints = extract_adc_tissue_interstitia(sol_mtk__, pbpk_simple, tmp_dose)
    tmpdf_pl_ints = extract_PL_tissue_interstitia(sol_mtk__, pbpk_simple, tmp_dose)
    tmpdf_pl_endo = extract_PL_tissue_endothelial_cells(sol_mtk__, pbpk_simple, tmp_dose)
    df_adc_ints = vcat(df_adc_ints, tmpdf_adc_ints)
    df_pl_ints = vcat(df_pl_ints, tmpdf_pl_ints)
    df_pl_endo = vcat(df_pl_endo, tmpdf_pl_endo)
end

# simulation, T-DM1 DAR-adjusted dose 
for tmp_dose in tdxd_tdm1_dar_adj
    adc__ = InfusionCallback(tmp_dose, pbpk_simple, infusion_d = add_dose);
    sol_mtk__ = solve(ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_tdxd), tspan, callback = adc__, infusion_hr = 0.1), alg=CVODE_BDF());
    tmpdf_adc_ints = extract_adc_tissue_interstitia(sol_mtk__, pbpk_simple, tmp_dose)
    tmpdf_pl_ints = extract_PL_tissue_interstitia(sol_mtk__, pbpk_simple, tmp_dose)
    tmpdf_pl_endo = extract_PL_tissue_endothelial_cells(sol_mtk__, pbpk_simple, tmp_dose)
    df_adc_ints = vcat(df_adc_ints, tmpdf_adc_ints)
    df_pl_ints = vcat(df_pl_ints, tmpdf_pl_ints)
    df_pl_endo = vcat(df_pl_endo, tmpdf_pl_endo)
end

# save outcome 
CSV.write("data/t-dxd-tissue-ints.csv", df_adc_ints)
CSV.write("data/t-dxd-pl-tissue-ints.csv", df_pl_ints)
CSV.write("data/t-dxd-pl-tissue-endo.csv", df_pl_endo)
