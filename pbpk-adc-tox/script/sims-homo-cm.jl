# date: 1/15/2026 
# author: Yuezhe Li 
# purpose of this code: to create simulation for cantuzumab mertansine
# additional dose added 

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
param_cm = create_base_pbpk_param(8, pbpk_simple, k_PL_ints_clearance = 0.34);  #  https://pubmed.ncbi.nlm.nih.gov/37787918/
param_cm[pbpk_simple.DAR] = 3.5
param_cm[pbpk_simple.k_diff] = 0.14  # as T-DM1, https://pubmed.ncbi.nlm.nih.gov/37787918/
param_cm[pbpk_simple.CL_plasma_PL] =  log(2)/3; # https://pubmed.ncbi.nlm.nih.gov/37787918/
param_cm[pbpk_simple.k_deconj] =  0.13
param_cm[pbpk_simple.VMAX] = 0

# simulation time
tspan = (0.0, hr_per_day*84);      # [hr]

# Tolcher et al., 2003; https://ascopubs.org/doi/10.1200/JCO.2003.05.137
# DLT: grade 3 elevation in hepatic transaminase (course 1)
# MTD: 235 mg/m2
add_dose_q3w = [0., 21., 42., 63]   
dose_q3w_2003 = [22, 44, 88, 132, 176, 235, 295]; # [mg/m2]


# additional dosing scheme from Heft et al., 2004
# https://aacrjournals.org/clincancerres/article/10/13/4363/94515/A-Phase-I-Study-of-Cantuzumab-Mertansine
# dosing QW; DLT is liver tox, at 115 mg/m2 QW
dose_qw_2004 = [40, 60, 80, 96, 115, 138]; # [mg/m2]
AddDose_qw_2004 = [0., 7., 14., 21., 28., 35., 42., 49., 56., 63., 70., 77.]  


# additional dose scheme from Rodin et al., 2008 
# https://link.springer.com/article/10.1007/s00280-007-0672-8
# dosing 3 times a week, 3 weeks per 4 weeks (3 weeks on, 1 week off)
dose_qw34_2008 = [30, 45, 60]; # [mg/m2]
AddDose_qw34_2008 = [0, 2, 4, 7, 9, 11, 14, 16, 18, 
                28, 30, 32, 35, 37, 39, 42, 44, 46, 
                56, 58, 60, 63, 65, 67, 70, 72, 74]  


# simulation dose 
tmd1_sims_dose = [0.3, 0.6, 1.2, 2.4, 3.6, 4.8]; # [mg/kg]

# simulation (1.8 mg/kg) (as T-DM1)
adc_1point8 = InfusionCallback(1.8, pbpk_simple, infusion_d = add_dose_q3w);
sol_mtk_1point8 = solve(ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_cm), tspan, callback = adc_1point8, infusion_hr = 0.5), alg=CVODE_BDF());
df_adc_ints = extract_adc_tissue_interstitia(sol_mtk_1point8, pbpk_simple, 1.8)
df_pl_ints = extract_PL_tissue_interstitia(sol_mtk_1point8, pbpk_simple, 1.8)
df_pl_endo = extract_PL_tissue_endothelial_cells(sol_mtk_1point8, pbpk_simple, 1.8)

# simulation (T-DM1 doses)
for tmp_dose in tmd1_sims_dose
    adc__ = InfusionCallback(tmp_dose, pbpk_simple, infusion_d = add_dose_q3w);
    sol_mtk__ = solve(ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_cm), tspan, callback = adc__, infusion_hr = 0.5), alg=CVODE_BDF());
    tmpdf_adc_ints = extract_adc_tissue_interstitia(sol_mtk__, pbpk_simple, tmp_dose)
    tmpdf_pl_ints = extract_PL_tissue_interstitia(sol_mtk__, pbpk_simple, tmp_dose)
    tmpdf_pl_endo = extract_PL_tissue_endothelial_cells(sol_mtk__, pbpk_simple, tmp_dose)
    df_adc_ints = vcat(df_adc_ints, tmpdf_adc_ints)
    df_pl_ints = vcat(df_pl_ints, tmpdf_pl_ints)
    df_pl_endo = vcat(df_pl_endo, tmpdf_pl_endo)
end

# simulation, Q3W 
for tmp_dose in dose_q3w_2003
    adc__ = InfusionCallback(tmp_dose*3/human_WT, pbpk_simple, infusion_d = add_dose_q3w);
    sol_mtk__ = solve(ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_cm), tspan, callback = adc__, infusion_hr = 0.5), alg=CVODE_BDF());
    tmpdf_adc_ints = extract_adc_tissue_interstitia(sol_mtk__, pbpk_simple, tmp_dose)
    tmpdf_pl_ints = extract_PL_tissue_interstitia(sol_mtk__, pbpk_simple, tmp_dose)
    tmpdf_pl_endo = extract_PL_tissue_endothelial_cells(sol_mtk__, pbpk_simple, tmp_dose)
    df_adc_ints = vcat(df_adc_ints, tmpdf_adc_ints)
    df_pl_ints = vcat(df_pl_ints, tmpdf_pl_ints)
    df_pl_endo = vcat(df_pl_endo, tmpdf_pl_endo)
end

# simulation, QW 
for tmp_dose in dose_qw_2004
    adc__ = InfusionCallback(tmp_dose*3/human_WT, pbpk_simple, infusion_d = AddDose_qw_2004);
    sol_mtk__ = solve(ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_cm), tspan, callback = adc__, infusion_hr = 0.5), alg=CVODE_BDF());
    tmpdf_adc_ints = extract_adc_tissue_interstitia(sol_mtk__, pbpk_simple, tmp_dose)
    tmpdf_pl_ints = extract_PL_tissue_interstitia(sol_mtk__, pbpk_simple, tmp_dose)
    tmpdf_pl_endo = extract_PL_tissue_endothelial_cells(sol_mtk__, pbpk_simple, tmp_dose)
    df_adc_ints = vcat(df_adc_ints, tmpdf_adc_ints)
    df_pl_ints = vcat(df_pl_ints, tmpdf_pl_ints)
    df_pl_endo = vcat(df_pl_endo, tmpdf_pl_endo)
end

# simulation, 3 times a week, every 3 out of 4 weeks 
for tmp_dose in dose_qw34_2008
    adc__ = InfusionCallback(tmp_dose*3/human_WT, pbpk_simple, infusion_d = AddDose_qw34_2008);
    sol_mtk__ = solve(ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_cm), tspan, callback = adc__, infusion_hr = 0.5), alg=CVODE_BDF());
    tmpdf_adc_ints = extract_adc_tissue_interstitia(sol_mtk__, pbpk_simple, tmp_dose)
    tmpdf_pl_ints = extract_PL_tissue_interstitia(sol_mtk__, pbpk_simple, tmp_dose)
    tmpdf_pl_endo = extract_PL_tissue_endothelial_cells(sol_mtk__, pbpk_simple, tmp_dose)
    df_adc_ints = vcat(df_adc_ints, tmpdf_adc_ints)
    df_pl_ints = vcat(df_pl_ints, tmpdf_pl_ints)
    df_pl_endo = vcat(df_pl_endo, tmpdf_pl_endo)
end

# visualization 
plot(df_adc_ints.time_hr, df_adc_ints.plasma, group = df_adc_ints.Dose, legend = :outerright, ylims = [1E-5, 1E1], yaxis = :log10)
plot(df_pl_ints.time_hr, df_pl_ints.liver, group = df_pl_ints.Dose, legend = :outerright, ylims = [1E-5, 1E1], yaxis = :log10)
plot(df_pl_endo.time_hr, df_pl_endo.liver, group = df_pl_endo.Dose, legend = :outerright, ylims = [1E-5, 1E1], yaxis = :log10)

# save outcome 
CSV.write("data/CanAg-dm1-tissue-ints.csv", df_adc_ints)
CSV.write("data/CanAg-dm1-pl-tissue-ints.csv", df_pl_ints)
CSV.write("data/CanAg-dm1-pl-tissue-endo.csv", df_pl_endo)
