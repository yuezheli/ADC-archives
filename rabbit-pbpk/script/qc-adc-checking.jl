# date: 11/25/2025 
# author: Yuezhe Li 
# purpose of this code: to track ADC across organs immediately after dosing 
# the goal was to see if the mass balance in the model holds (without accounting for ADC degredation)
# this idea was proposed during client meeting, 11/24/2025 
# ABBV706 was used as an example here bc it does not have TMDD

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
# helper functions
include(@projectroot("script/helper-funcs.jl"));
include(@projectroot("script/helper-initiation.jl")); 
include(@projectroot("script/helper-parameters.jl"));
include(@projectroot("script/helper-infusion.jl"));

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
param_abbv706 = create_base_pbpk_param(7, pbpk_simple);
param_abbv706[pbpk_simple.DAR] = 6
param_abbv706[pbpk_simple.CL_plasma_PL] = 1124/353 * 4; # computed based on CL/V from NONMEM file
# cornea-related parameters tuned in the rabbit model 
param_abbv706[pbpk_simple.PS_cor] = 1.19E-8 * 4E2

## infusion, 
adc_infusion_2point5 = InfusionCallback(2.5, pbpk_simple, BW = 90, MW_EDG = 14.5E4);

## simulation 
@time prob_mtk_2point5 = ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_abbv706), tspan, callback = adc_infusion_2point5);  
@time sol_mtk_infusion_2point5 = solve(prob_mtk_2point5, alg=CVODE_BDF());


#####################

# compute ADC mass in a organ with vasculature, endothelial cells, and interstitium
function ComputeTissueADCMass(C_V, 
    C_VM, C_bound_VM_mem, C_bound_VM, C_bound2_VM, 
    C_E7, C_bound_E7, C_bound2_E7, C_E6a, C_bound_E6a, C_bound2_E6a, C_E7b, C_bound_E7b, C_bound2_E7b, 
    C_ISM, C_bound_ISM_mem, C_bound_ISM, C_bound2_ISM, C_IntS, 
    V_V, V_IntS, Endothelial_Cell_Frac; 
    pino_time = 10.8/60, CL_up_in_nL_per_hour_per_million_cells = 150, Total_Endothelial_Cell = 1.422e9, Scale_Factor = 603.7, 
    tau_VM = 1/60.0, tau_ISM = 1/60.0, E6a_Vol_Pct = 0.33)

    # compute organ volumes 
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    V_endosomal = CL_up_in_nL_per_hour_per_million_cells * pino_time * Organ_Endothelial_Cell * 1e-6 * 1e-9
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_VM = CL_up * tau_VM
    V_E6a = V_endosomal * E6a_Vol_Pct
    E7_Vol_Pct = (1.0 - E6a_Vol_Pct) / 2 # [-]
    E7b_Vol_Pct = (1.0 - E6a_Vol_Pct) / 2 # [-]
    V_E7 = V_endosomal * E7_Vol_Pct
    V_E7b = V_endosomal * E7b_Vol_Pct
    V_ISM = CL_up * tau_ISM

    # compute total ADC mass
    total_ADC_mass =  (
        C_V * V_V 
        + (C_VM + C_bound_VM_mem + C_bound_VM + C_bound2_VM) * V_VM 
        + (C_E7 + C_bound_E7 + C_bound2_E7) * V_E7 
        + (C_E6a + C_bound_E6a + C_bound2_E6a) * V_E6a 
        + (C_E7b + C_bound_E7b + C_bound2_E7b) * V_E7b  
        + (C_ISM + C_bound_ISM_mem + C_bound_ISM + C_bound2_ISM) * V_ISM 
        + C_IntS * V_IntS
    )

    return total_ADC_mass
end

# purpose of this function: to get the name of the variable by string 
function getVariableIndex_byname(var_name, mdl)
    all_sys_vars = (unknowns(mdl))
    all_sys_names = string.(all_sys_vars)
    tmpindex = findfirst(==(var_name), all_sys_names)
    return(all_sys_vars[tmpindex])
end

# purpose of this function: to get the name of the parameter by string 
function getParameterIndex_byname(param_name, mdl)
    all_sys_params = (parameters(mdl))
    all_sys_param_names = string.(all_sys_params)
    tmpindex = findfirst(==(param_name), all_sys_param_names)
    return(all_sys_params[tmpindex])
end


# compute organ total ADC by organ name 
function ComputeOrganADC_massbyname(sol_mtk_, organ_name, param_global, mdl)
    # pull variable names (by organ)
    ind_igg_C_V = getVariableIndex_byname( organ_name*"₊igg_exg₊C_V(t)", mdl ); 
    ind_igg_C_VM = getVariableIndex_byname( organ_name*"₊igg_exg₊C_VM(t)", mdl ); 
    ind_igg_C_bound_VM_mem = getVariableIndex_byname( organ_name*"₊igg_exg₊C_bound_VM_mem(t)", mdl ); 
    ind_igg_C_bound_VM = getVariableIndex_byname( organ_name*"₊igg_exg₊C_bound_VM(t)", mdl ); 
    ind_igg_C_bound2_VM = getVariableIndex_byname( organ_name*"₊igg_exg₊C_bound2_VM(t)", mdl ); 
    ind_igg_C_E7 = getVariableIndex_byname( organ_name*"₊igg_exg₊C_E7(t)", mdl ); 
    ind_igg_C_bound_E7 = getVariableIndex_byname( organ_name*"₊igg_exg₊C_bound_E7(t)", mdl ); 
    ind_igg_C_bound2_E7 = getVariableIndex_byname( organ_name*"₊igg_exg₊C_bound2_E7(t)", mdl ); 
    ind_igg_C_E6a = getVariableIndex_byname( organ_name*"₊igg_exg₊C_E6a(t)", mdl ); 
    ind_igg_C_bound_E6a = getVariableIndex_byname( organ_name*"₊igg_exg₊C_bound_E6a(t)", mdl ); 
    ind_igg_C_bound2_E6a = getVariableIndex_byname( organ_name*"₊igg_exg₊C_bound2_E6a(t)", mdl ); 
    ind_igg_C_E7b = getVariableIndex_byname( organ_name*"₊igg_exg₊C_E7b(t)", mdl ); 
    ind_igg_C_bound_E7b = getVariableIndex_byname( organ_name*"₊igg_exg₊C_bound_E7b(t)", mdl ); 
    ind_igg_C_bound2_E7b = getVariableIndex_byname( organ_name*"₊igg_exg₊C_bound2_E7b(t)", mdl ); 
    ind_igg_C_ISM = getVariableIndex_byname( organ_name*"₊igg_exg₊C_ISM(t)", mdl ); 
    ind_igg_C_bound_ISM_mem = getVariableIndex_byname( organ_name*"₊igg_exg₊C_bound_ISM_mem(t)", mdl ); 
    ind_igg_C_bound_ISM = getVariableIndex_byname( organ_name*"₊igg_exg₊C_bound_ISM(t)", mdl ); 
    ind_igg_C_bound2_ISM = getVariableIndex_byname( organ_name*"₊igg_exg₊C_bound2_ISM(t)", mdl ); 
    ind_igg_C_IntS = getVariableIndex_byname( organ_name*"₊igg_exg₊C_IntS(t)", mdl ); 
    
    # get parameter index 
    ind_V_V = getParameterIndex_byname(organ_name*"₊V_V", mdl)
    ind_V_IntS = getParameterIndex_byname(organ_name*"₊V_IntS", mdl)
    ind_Endothelial_Cell_Frac = getParameterIndex_byname(organ_name*"₊Endothelial_Cell_Frac", mdl)

    organADC = ComputeTissueADCMass(sol_mtk_[ind_igg_C_V], 
    sol_mtk_[ind_igg_C_VM], sol_mtk_[ind_igg_C_bound_VM_mem], sol_mtk_[ind_igg_C_bound_VM], sol_mtk_[ind_igg_C_bound2_VM], 
    sol_mtk_[ind_igg_C_E7], sol_mtk_[ind_igg_C_bound_E7], sol_mtk_[ind_igg_C_bound2_E7], sol_mtk_[ind_igg_C_E6a], 
    sol_mtk_[ind_igg_C_bound_E6a], sol_mtk_[ind_igg_C_bound2_E6a], sol_mtk_[ind_igg_C_E7b], sol_mtk_[ind_igg_C_bound_E7b], sol_mtk_[ind_igg_C_bound2_E7b], 
    sol_mtk_[ind_igg_C_ISM], sol_mtk_[ind_igg_C_bound_ISM_mem], sol_mtk_[ind_igg_C_bound_ISM], sol_mtk_[ind_igg_C_bound2_ISM], sol_mtk_[ind_igg_C_IntS], 
    param_global[ind_V_V], param_global[ind_V_IntS], param_global[ind_Endothelial_Cell_Frac])

    return organADC
end

# compute ADC across all organs and convert to a data frame 
function CollectADCs_across(sol_mtk_, param_global, mdl)
    # compute in the eyes (ICB, retina, choroid)
    # this section was copied from script/helper-compute-tissue-adc.jl
    # this part didn't use the function defined above and called below due to inconsistency in capitalization 
    ADC_Choroid = ComputeTissueADCMass(sol_mtk_[mdl.eye.choroid.igg_exg.C_V], 
    sol_mtk_[mdl.eye.choroid.igg_exg.C_VM], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound_VM_mem], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound_VM], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound2_VM], 
    sol_mtk_[mdl.eye.choroid.igg_exg.C_E7], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound_E7], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound2_E7], sol_mtk_[mdl.eye.choroid.igg_exg.C_E6a], 
    sol_mtk_[mdl.eye.choroid.igg_exg.C_bound_E6a], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound2_E6a], sol_mtk_[mdl.eye.choroid.igg_exg.C_E7b], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound_E7b], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound2_E7b], 
    sol_mtk_[mdl.eye.choroid.igg_exg.C_ISM], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound_ISM_mem], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound_ISM], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound2_ISM], sol_mtk_[mdl.eye.choroid.igg_exg.C_IntS], 
    param_global[mdl.eye.choroid.V_V], param_global[mdl.eye.choroid.V_IntS], param_global[mdl.eye.choroid.Endothelial_Cell_Frac])

    ADC_ICB = ComputeTissueADCMass(sol_mtk_[mdl.eye.icb.igg_exg.C_V], 
    sol_mtk_[mdl.eye.icb.igg_exg.C_VM], sol_mtk_[mdl.eye.icb.igg_exg.C_bound_VM_mem], sol_mtk_[mdl.eye.icb.igg_exg.C_bound_VM], sol_mtk_[mdl.eye.icb.igg_exg.C_bound2_VM], 
    sol_mtk_[mdl.eye.icb.igg_exg.C_E7], sol_mtk_[mdl.eye.icb.igg_exg.C_bound_E7], sol_mtk_[mdl.eye.icb.igg_exg.C_bound2_E7], sol_mtk_[mdl.eye.icb.igg_exg.C_E6a], 
    sol_mtk_[mdl.eye.icb.igg_exg.C_bound_E6a], sol_mtk_[mdl.eye.icb.igg_exg.C_bound2_E6a], sol_mtk_[mdl.eye.icb.igg_exg.C_E7b], sol_mtk_[mdl.eye.icb.igg_exg.C_bound_E7b], sol_mtk_[mdl.eye.icb.igg_exg.C_bound2_E7b], 
    sol_mtk_[mdl.eye.icb.igg_exg.C_ISM], sol_mtk_[mdl.eye.icb.igg_exg.C_bound_ISM_mem], sol_mtk_[mdl.eye.icb.igg_exg.C_bound_ISM], sol_mtk_[mdl.eye.icb.igg_exg.C_bound2_ISM], sol_mtk_[mdl.eye.icb.igg_exg.C_IntS], 
    param_global[mdl.eye.icb.V_V], param_global[mdl.eye.icb.V_IntS], param_global[mdl.eye.icb.Endothelial_Cell_Frac])

    ADC_Retina = ComputeTissueADCMass(sol_mtk_[mdl.eye.retina.igg_exg.C_V], 
    sol_mtk_[mdl.eye.retina.igg_exg.C_VM], sol_mtk_[mdl.eye.retina.igg_exg.C_bound_VM_mem], sol_mtk_[mdl.eye.retina.igg_exg.C_bound_VM], sol_mtk_[mdl.eye.retina.igg_exg.C_bound2_VM], 
    sol_mtk_[mdl.eye.retina.igg_exg.C_E7], sol_mtk_[mdl.eye.retina.igg_exg.C_bound_E7], sol_mtk_[mdl.eye.retina.igg_exg.C_bound2_E7], sol_mtk_[mdl.eye.retina.igg_exg.C_E6a], 
    sol_mtk_[mdl.eye.retina.igg_exg.C_bound_E6a], sol_mtk_[mdl.eye.retina.igg_exg.C_bound2_E6a], sol_mtk_[mdl.eye.retina.igg_exg.C_E7b], sol_mtk_[mdl.eye.retina.igg_exg.C_bound_E7b], sol_mtk_[mdl.eye.retina.igg_exg.C_bound2_E7b], 
    sol_mtk_[mdl.eye.retina.igg_exg.C_ISM], sol_mtk_[mdl.eye.retina.igg_exg.C_bound_ISM_mem], sol_mtk_[mdl.eye.retina.igg_exg.C_bound_ISM], sol_mtk_[mdl.eye.retina.igg_exg.C_bound2_ISM], sol_mtk_[mdl.eye.retina.igg_exg.C_IntS], 
    param_global[mdl.eye.retina.V_V], param_global[mdl.eye.retina.V_IntS], param_global[mdl.eye.retina.Endothelial_Cell_Frac])

    df = DataFrame(
        # time 
        time_hr = sol_mtk_.t, 
        # total tissue ADC payload mass
        lung = ComputeOrganADC_massbyname(sol_mtk_, "lung", param_global, mdl), 
        liver = ComputeOrganADC_massbyname(sol_mtk_, "liver", param_global, mdl), 
        heart = ComputeOrganADC_massbyname(sol_mtk_, "heart", param_global, mdl), 
        muscle = ComputeOrganADC_massbyname(sol_mtk_, "muscle", param_global, mdl), 
        skin = ComputeOrganADC_massbyname(sol_mtk_, "skin", param_global, mdl), 
        adipose = ComputeOrganADC_massbyname(sol_mtk_, "adipose", param_global, mdl), 
        bone = ComputeOrganADC_massbyname(sol_mtk_, "bone", param_global, mdl), 
        brain = ComputeOrganADC_massbyname(sol_mtk_, "brain", param_global, mdl), 
        kidney = ComputeOrganADC_massbyname(sol_mtk_, "kidney", param_global, mdl), 
        sm_int = ComputeOrganADC_massbyname(sol_mtk_, "sm_int", param_global, mdl), 
        la_int = ComputeOrganADC_massbyname(sol_mtk_, "la_int", param_global, mdl), 
        pancreas = ComputeOrganADC_massbyname(sol_mtk_, "pancreas", param_global, mdl), 
        thymus = ComputeOrganADC_massbyname(sol_mtk_, "thymus", param_global, mdl), 
        spleen = ComputeOrganADC_massbyname(sol_mtk_, "spleen", param_global, mdl), 
        other = ComputeOrganADC_massbyname(sol_mtk_, "other", param_global, mdl), 
        marrow = ComputeOrganADC_massbyname(sol_mtk_, "marrow", param_global, mdl), 
        # ADC concentration in the eye 
        choroid = ADC_Choroid, 
        retina = ADC_Retina, 
        icb = ADC_ICB, 
        ah = sol_mtk_[mdl.eye.ah_igg_exg.C_AQ] * param_global[getParameterIndex_byname("V_AQ", mdl)], 
        vh = sol_mtk_[mdl.eye.vh_igg_exg.C_VH] * param_global[getParameterIndex_byname("V_VH", mdl)], 
        cornea = sol_mtk_[mdl.eye.igg_exg.C_COR] * param_global[getParameterIndex_byname("V_Cor", mdl)], 
        # ADC concentration in plasma 
        plasma = sol_mtk_[mdl.plasma_exg.C_Plasma]  * param_global[getParameterIndex_byname("V_Plasma", mdl)], 
        # ADC concentration in LN 
        ln = sol_mtk_[mdl.plasma_exg.C_LN]  * param_global[getParameterIndex_byname("V_LN", mdl)], 
    )

    return df
end


# compute ADC across all organs in the pbpk model 
all_conj_adc = CollectADCs_across(sol_mtk_infusion_2point5, param_abbv706, pbpk_simple); 

# compute total ADC 
tot_conj_adc = [];
for i in 1:size(all_conj_adc)[1]
    tmp_adc = sum(all_conj_adc[i, 2:end])
    append!(tot_conj_adc, tmp_adc)
end

all_conj_adc.tot_adc .= tot_conj_adc;

# the total amount of ADC after infusion ends  (infusion ends at 1 hr)
@rsubset(all_conj_adc, :time_hr >=1)[1, "tot_adc"]

# compute ADC dose 
# dose = 2.5 mg/kg, BW = 90 kg, ADC MW = 14.5E4 g/mol
total_dose_umol = 2.5*90*1E-3/14.5E4*1e6; 

# display 
final_df = DataFrame(
    infused_ADC_umol = total_dose_umol, 
    conj_adc_in_pbpk_after_infusion_ended = @rsubset(all_conj_adc, :time_hr >=1)[1, "tot_adc"]
)

println( final_df )

