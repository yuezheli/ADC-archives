# date: 1/15/2025 
# author: Yuezhe Li 
# purpose of this code: to create simulation for T-DM1 (for further analysis in R)
# additional dose was added 

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
param_tdm1 = create_base_pbpk_param(5, pbpk_simple, k_PL_ints_clearance = 0.34);  #  https://pubmed.ncbi.nlm.nih.gov/37787918/
param_tdm1[pbpk_simple.VMAX] = 0.0012
param_tdm1[pbpk_simple.KM] = 0.006
param_tdm1[pbpk_simple.DAR] = 3.5
param_tdm1[pbpk_simple.k_diff] = 0.14  #  https://pubmed.ncbi.nlm.nih.gov/37787918/
param_tdm1[pbpk_simple.k_deconj] =  6E-3
param_tdm1[pbpk_simple.CL_plasma_PL] =  4
param_tdm1[pbpk_simple.Kp_adc_cor_ah] = 1

# simulation time
tspan = (-0.01, hr_per_day*84);      # [hr]
add_dose = [0., 21., 42., 63]   

# simulation dose 
tmd1_sims_dose = [0.3, 0.6, 1.2, 2.4, 3.6, 4.8]; # [mg/kg]
tdxd_sims_dose = [0.8, 1.6, 3.2, 5.4, 6.4, 8.0]; # [mg/kg]
tdm1_tdxd_dar_adj = tdxd_sims_dose * 8/3.5

# simulation (1.8 mg/kg) (Japanese patients) # https://pmc.ncbi.nlm.nih.gov/articles/PMC3998859/
adc_1point8 = InfusionCallback(1.8, pbpk_simple, infusion_d = add_dose);
sol_mtk_1point8 = solve(ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_tdm1), tspan, callback = adc_1point8, infusion_hr = 0.1), alg=CVODE_BDF());
df_adc_ints = extract_adc_tissue_interstitia(sol_mtk_1point8, pbpk_simple, 1.8)
df_pl_ints = extract_PL_tissue_interstitia(sol_mtk_1point8, pbpk_simple, 1.8)
df_pl_endo = extract_PL_tissue_endothelial_cells(sol_mtk_1point8, pbpk_simple, 1.8)

# simulation, T-DM1 dose
for tmp_dose in tmd1_sims_dose
    adc__ = InfusionCallback(tmp_dose, pbpk_simple, infusion_d = add_dose);
    sol_mtk__ = solve(ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_tdm1), tspan, callback = adc__, infusion_hr = 0.1), alg=CVODE_BDF());
    tmpdf_adc_ints = extract_adc_tissue_interstitia(sol_mtk__, pbpk_simple, tmp_dose)
    tmpdf_pl_ints = extract_PL_tissue_interstitia(sol_mtk__, pbpk_simple, tmp_dose)
    tmpdf_pl_endo = extract_PL_tissue_endothelial_cells(sol_mtk__, pbpk_simple, tmp_dose)
    df_adc_ints = vcat(df_adc_ints, tmpdf_adc_ints)
    df_pl_ints = vcat(df_pl_ints, tmpdf_pl_ints)
    df_pl_endo = vcat(df_pl_endo, tmpdf_pl_endo)
end

# simulation, DAR-adj T-Dxd dose 
for tmp_dose in tdm1_tdxd_dar_adj
    adc__ = InfusionCallback(tmp_dose, pbpk_simple, infusion_d = add_dose);
    sol_mtk__ = solve(ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_tdm1), tspan, callback = adc__, infusion_hr = 0.1), alg=CVODE_BDF());
    tmpdf_adc_ints = extract_adc_tissue_interstitia(sol_mtk__, pbpk_simple, tmp_dose)
    tmpdf_pl_ints = extract_PL_tissue_interstitia(sol_mtk__, pbpk_simple, tmp_dose)
    tmpdf_pl_endo = extract_PL_tissue_endothelial_cells(sol_mtk__, pbpk_simple, tmp_dose)
    df_adc_ints = vcat(df_adc_ints, tmpdf_adc_ints)
    df_pl_ints = vcat(df_pl_ints, tmpdf_pl_ints)
    df_pl_endo = vcat(df_pl_endo, tmpdf_pl_endo)
end

# save outcome 
CSV.write("data/t-dm1-tissue-ints.csv", df_adc_ints)
CSV.write("data/t-dm1-pl-tissue-ints.csv", df_pl_ints)
CSV.write("data/t-dm1-pl-tissue-endo.csv", df_pl_endo)

