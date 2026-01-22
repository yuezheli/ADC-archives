# date: 1/20/2026 
# author: Yuezhe Li 
# purpose of this code: to generate PK simulation for ADCT-402

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

# observed PK 
obs_ = CSV.read(@projectroot("data/loncastuximab-tesirine/clinical-data1.csv"), DataFrame, missingstrings=".");
@rsubset!(obs_, :BLQ == 0);
@rsubset!(obs_, :DOSFRQ == "Q3W");

obs_pk_adc = @rsubset(obs_, :DVID == 1);
obs_pk_pl = @rsubset(obs_, :DVID == 2);

# convert observed data to dictionary 
dict_pk_adc = Dict();
for dose_ in [60, 90, 200]
    tmp_pk_ = @rsubset(obs_pk_adc, :ARMCD == dose_);
    tmp_df = DataFrame(
        time_d = tmp_pk_.TIME/hr_per_day, 
        ADC_uM = tmp_pk_.DV/MW_IGG
    )
    dose_mgkg = dose_/1E3
    push!(dict_pk_adc, dose_mgkg => tmp_df)
end

dict_pk_pl = Dict();
for dose_ in [60, 90, 200]
    tmp_pk_ = @rsubset(obs_pk_pl, :ARMCD == dose_);
    tmp_df = DataFrame(
        time_d = tmp_pk_.TIME/hr_per_day, 
        pl_uM = tmp_pk_.DV/MW_PBD 
    )
    dose_mgkg = dose_/1E3
    push!(dict_pk_pl, dose_mgkg => tmp_df)
end

#####################

## create pbpk model
@time pbpk = create_pbpk();
@time pbpk_simple = mtkcompile(pbpk);

## initial condition (IV dosing)
u0_infusion = pbpk_initial_condition(0, pbpk_simple); 

## parameters 
### AC-SINS score, VMAX, KM all from T-Dxd fitting
param_402 = create_base_pbpk_param(6, pbpk_simple, k_PL_ints_clearance = 0.);  # https://pubmed.ncbi.nlm.nih.gov/38406517/
param_402[pbpk_simple.DAR] = 2
param_402[pbpk_simple.k_diff] = 3 * 11E-6 / cell_radius * s_per_hr  # https://pubmed.ncbi.nlm.nih.gov/32787101/
param_402[pbpk_simple.PS_cor_pl] =  11E-6 * 126E-6 * s_per_hr
param_402[pbpk_simple.PS_ret_pl] =  11E-6 * 1363E-6 * s_per_hr
param_402[pbpk_simple.CL_plasma_PL] =  0.2

# simulation time
tspan = (-0.01, hr_per_day*42);      # [hr]
add_dose = [0., 21.]   

# simulation dose 
sims_dose = [0.06, 0.09, 0.2]; # [mg/kg]

# simulation
sol_pk_402 = Dict(); 
for init_dose in sims_dose
    adc_infusion_ = InfusionCallback(init_dose, pbpk_simple, infusion_d = add_dose);
    # simulation 
    @time prob_mtk_ = ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_402), tspan, callback = adc_infusion_, infusion_hr = 0.1); 
    @time sol_mtk_infusion_ = solve(prob_mtk_, alg=CVODE_BDF());
    # append result 
    push!(sol_pk_402, init_dose => sol_mtk_infusion_);
end

#####################
# visualization 
plt_plasma_adc = PlotSimulationPlasma(sol_pk_402, dict_pk_adc, pbpk_simple, 
    adc_name = "loncastuximab tesirine", colorPALETTE = :Set2_3, 
    xrange = [0, 42], yrange = [1E-6, 1], ylog = true, legendcolumnsnum = 2)


plt_plasma_pl = PlotSimulationPlasmaPL(sol_pk_402, dict_pk_pl, pbpk_simple, 
    adc_name = "PBD", colorPALETTE = :Set2_3, 
    xrange = [0, 42], yrange = [1E-6, 1], ylog = true, legendcolumnsnum = 2)


savefig(plt_plasma_adc, "deliv/figure/pk-adct-402-adc-plasma.png");
savefig(plt_plasma_pl, "deliv/figure/pk-adct-402-pbd-plasma.png");
