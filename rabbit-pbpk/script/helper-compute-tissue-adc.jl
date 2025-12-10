# date: 10/3/2025 
# author: Yuezhe Li 
# purpose of this code: to compute whole tissue ADC concentration 

function ComputeTissueADC(C_V, 
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

    V_Organ = V_V + V_endosomal + V_IntS

    return total_ADC_mass/V_Organ
end

function OrganAvgADC(sol_mtk_, param_global, mdl = pbpk_simple)
    ADC_LI = ComputeTissueADC(sol_mtk_[mdl.la_int.igg_exg.C_V], 
    sol_mtk_[mdl.la_int.igg_exg.C_VM], sol_mtk_[mdl.la_int.igg_exg.C_bound_VM_mem], sol_mtk_[mdl.la_int.igg_exg.C_bound_VM], sol_mtk_[mdl.la_int.igg_exg.C_bound2_VM], 
    sol_mtk_[mdl.la_int.igg_exg.C_E7], sol_mtk_[mdl.la_int.igg_exg.C_bound_E7], sol_mtk_[mdl.la_int.igg_exg.C_bound2_E7], sol_mtk_[mdl.la_int.igg_exg.C_E6a], 
    sol_mtk_[mdl.la_int.igg_exg.C_bound_E6a], sol_mtk_[mdl.la_int.igg_exg.C_bound2_E6a], sol_mtk_[mdl.la_int.igg_exg.C_E7b], sol_mtk_[mdl.la_int.igg_exg.C_bound_E7b], sol_mtk_[mdl.la_int.igg_exg.C_bound2_E7b], 
    sol_mtk_[mdl.la_int.igg_exg.C_ISM], sol_mtk_[mdl.la_int.igg_exg.C_bound_ISM_mem], sol_mtk_[mdl.la_int.igg_exg.C_bound_ISM], sol_mtk_[mdl.la_int.igg_exg.C_bound2_ISM], sol_mtk_[mdl.la_int.igg_exg.C_IntS], 
    param_global[mdl.la_int.V_V], param_global[mdl.la_int.V_IntS], param_global[mdl.la_int.Endothelial_Cell_Frac])

    ADC_SI = ComputeTissueADC(sol_mtk_[mdl.sm_int.igg_exg.C_V], 
    sol_mtk_[mdl.sm_int.igg_exg.C_VM], sol_mtk_[mdl.sm_int.igg_exg.C_bound_VM_mem], sol_mtk_[mdl.sm_int.igg_exg.C_bound_VM], sol_mtk_[mdl.sm_int.igg_exg.C_bound2_VM], 
    sol_mtk_[mdl.sm_int.igg_exg.C_E7], sol_mtk_[mdl.sm_int.igg_exg.C_bound_E7], sol_mtk_[mdl.sm_int.igg_exg.C_bound2_E7], sol_mtk_[mdl.sm_int.igg_exg.C_E6a], 
    sol_mtk_[mdl.sm_int.igg_exg.C_bound_E6a], sol_mtk_[mdl.sm_int.igg_exg.C_bound2_E6a], sol_mtk_[mdl.sm_int.igg_exg.C_E7b], sol_mtk_[mdl.sm_int.igg_exg.C_bound_E7b], sol_mtk_[mdl.sm_int.igg_exg.C_bound2_E7b], 
    sol_mtk_[mdl.sm_int.igg_exg.C_ISM], sol_mtk_[mdl.sm_int.igg_exg.C_bound_ISM_mem], sol_mtk_[mdl.sm_int.igg_exg.C_bound_ISM], sol_mtk_[mdl.sm_int.igg_exg.C_bound2_ISM], sol_mtk_[mdl.sm_int.igg_exg.C_IntS], 
    param_global[mdl.sm_int.V_V], param_global[mdl.sm_int.V_IntS], param_global[mdl.sm_int.Endothelial_Cell_Frac])

    ADC_Lung = ComputeTissueADC(sol_mtk_[mdl.lung.igg_exg.C_V], 
    sol_mtk_[mdl.lung.igg_exg.C_VM], sol_mtk_[mdl.lung.igg_exg.C_bound_VM_mem], sol_mtk_[mdl.lung.igg_exg.C_bound_VM], sol_mtk_[mdl.lung.igg_exg.C_bound2_VM], 
    sol_mtk_[mdl.lung.igg_exg.C_E7], sol_mtk_[mdl.lung.igg_exg.C_bound_E7], sol_mtk_[mdl.lung.igg_exg.C_bound2_E7], sol_mtk_[mdl.lung.igg_exg.C_E6a], 
    sol_mtk_[mdl.lung.igg_exg.C_bound_E6a], sol_mtk_[mdl.lung.igg_exg.C_bound2_E6a], sol_mtk_[mdl.lung.igg_exg.C_E7b], sol_mtk_[mdl.lung.igg_exg.C_bound_E7b], sol_mtk_[mdl.lung.igg_exg.C_bound2_E7b], 
    sol_mtk_[mdl.lung.igg_exg.C_ISM], sol_mtk_[mdl.lung.igg_exg.C_bound_ISM_mem], sol_mtk_[mdl.lung.igg_exg.C_bound_ISM], sol_mtk_[mdl.lung.igg_exg.C_bound2_ISM], sol_mtk_[mdl.lung.igg_exg.C_IntS], 
    param_global[mdl.lung.V_V], param_global[mdl.lung.V_IntS], param_global[mdl.lung.Endothelial_Cell_Frac])

    ADC_Liver = ComputeTissueADC(sol_mtk_[mdl.liver.igg_exg.C_V], 
    sol_mtk_[mdl.liver.igg_exg.C_VM], sol_mtk_[mdl.liver.igg_exg.C_bound_VM_mem], sol_mtk_[mdl.liver.igg_exg.C_bound_VM], sol_mtk_[mdl.liver.igg_exg.C_bound2_VM], 
    sol_mtk_[mdl.liver.igg_exg.C_E7], sol_mtk_[mdl.liver.igg_exg.C_bound_E7], sol_mtk_[mdl.liver.igg_exg.C_bound2_E7], sol_mtk_[mdl.liver.igg_exg.C_E6a], 
    sol_mtk_[mdl.liver.igg_exg.C_bound_E6a], sol_mtk_[mdl.liver.igg_exg.C_bound2_E6a], sol_mtk_[mdl.liver.igg_exg.C_E7b], sol_mtk_[mdl.liver.igg_exg.C_bound_E7b], sol_mtk_[mdl.liver.igg_exg.C_bound2_E7b], 
    sol_mtk_[mdl.liver.igg_exg.C_ISM], sol_mtk_[mdl.liver.igg_exg.C_bound_ISM_mem], sol_mtk_[mdl.liver.igg_exg.C_bound_ISM], sol_mtk_[mdl.liver.igg_exg.C_bound2_ISM], sol_mtk_[mdl.liver.igg_exg.C_IntS], 
    param_global[mdl.liver.V_V], param_global[mdl.liver.V_IntS], param_global[mdl.liver.Endothelial_Cell_Frac])

    ADC_Marrow = ComputeTissueADC(sol_mtk_[mdl.marrow.igg_exg.C_V], 
    sol_mtk_[mdl.marrow.igg_exg.C_VM], sol_mtk_[mdl.marrow.igg_exg.C_bound_VM_mem], sol_mtk_[mdl.marrow.igg_exg.C_bound_VM], sol_mtk_[mdl.marrow.igg_exg.C_bound2_VM], 
    sol_mtk_[mdl.marrow.igg_exg.C_E7], sol_mtk_[mdl.marrow.igg_exg.C_bound_E7], sol_mtk_[mdl.marrow.igg_exg.C_bound2_E7], sol_mtk_[mdl.marrow.igg_exg.C_E6a], 
    sol_mtk_[mdl.marrow.igg_exg.C_bound_E6a], sol_mtk_[mdl.marrow.igg_exg.C_bound2_E6a], sol_mtk_[mdl.marrow.igg_exg.C_E7b], sol_mtk_[mdl.marrow.igg_exg.C_bound_E7b], sol_mtk_[mdl.marrow.igg_exg.C_bound2_E7b], 
    sol_mtk_[mdl.marrow.igg_exg.C_ISM], sol_mtk_[mdl.marrow.igg_exg.C_bound_ISM_mem], sol_mtk_[mdl.marrow.igg_exg.C_bound_ISM], sol_mtk_[mdl.marrow.igg_exg.C_bound2_ISM], sol_mtk_[mdl.marrow.igg_exg.C_IntS], 
    param_global[mdl.marrow.V_V], param_global[mdl.marrow.V_IntS], param_global[mdl.marrow.Endothelial_Cell_Frac])

    ADC_Skin = ComputeTissueADC(sol_mtk_[mdl.skin.igg_exg.C_V], 
    sol_mtk_[mdl.skin.igg_exg.C_VM], sol_mtk_[mdl.skin.igg_exg.C_bound_VM_mem], sol_mtk_[mdl.skin.igg_exg.C_bound_VM], sol_mtk_[mdl.skin.igg_exg.C_bound2_VM], 
    sol_mtk_[mdl.skin.igg_exg.C_E7], sol_mtk_[mdl.skin.igg_exg.C_bound_E7], sol_mtk_[mdl.skin.igg_exg.C_bound2_E7], sol_mtk_[mdl.skin.igg_exg.C_E6a], 
    sol_mtk_[mdl.skin.igg_exg.C_bound_E6a], sol_mtk_[mdl.skin.igg_exg.C_bound2_E6a], sol_mtk_[mdl.skin.igg_exg.C_E7b], sol_mtk_[mdl.skin.igg_exg.C_bound_E7b], sol_mtk_[mdl.skin.igg_exg.C_bound2_E7b], 
    sol_mtk_[mdl.skin.igg_exg.C_ISM], sol_mtk_[mdl.skin.igg_exg.C_bound_ISM_mem], sol_mtk_[mdl.skin.igg_exg.C_bound_ISM], sol_mtk_[mdl.skin.igg_exg.C_bound2_ISM], sol_mtk_[mdl.skin.igg_exg.C_IntS], 
    param_global[mdl.skin.V_V], param_global[mdl.skin.V_IntS], param_global[mdl.skin.Endothelial_Cell_Frac])

    ADC_Choroid = ComputeTissueADC(sol_mtk_[mdl.eye.choroid.igg_exg.C_V], 
    sol_mtk_[mdl.eye.choroid.igg_exg.C_VM], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound_VM_mem], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound_VM], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound2_VM], 
    sol_mtk_[mdl.eye.choroid.igg_exg.C_E7], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound_E7], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound2_E7], sol_mtk_[mdl.eye.choroid.igg_exg.C_E6a], 
    sol_mtk_[mdl.eye.choroid.igg_exg.C_bound_E6a], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound2_E6a], sol_mtk_[mdl.eye.choroid.igg_exg.C_E7b], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound_E7b], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound2_E7b], 
    sol_mtk_[mdl.eye.choroid.igg_exg.C_ISM], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound_ISM_mem], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound_ISM], sol_mtk_[mdl.eye.choroid.igg_exg.C_bound2_ISM], sol_mtk_[mdl.eye.choroid.igg_exg.C_IntS], 
    param_global[mdl.eye.choroid.V_V], param_global[mdl.eye.choroid.V_IntS], param_global[mdl.eye.choroid.Endothelial_Cell_Frac])

    ADC_ICB = ComputeTissueADC(sol_mtk_[mdl.eye.icb.igg_exg.C_V], 
    sol_mtk_[mdl.eye.icb.igg_exg.C_VM], sol_mtk_[mdl.eye.icb.igg_exg.C_bound_VM_mem], sol_mtk_[mdl.eye.icb.igg_exg.C_bound_VM], sol_mtk_[mdl.eye.icb.igg_exg.C_bound2_VM], 
    sol_mtk_[mdl.eye.icb.igg_exg.C_E7], sol_mtk_[mdl.eye.icb.igg_exg.C_bound_E7], sol_mtk_[mdl.eye.icb.igg_exg.C_bound2_E7], sol_mtk_[mdl.eye.icb.igg_exg.C_E6a], 
    sol_mtk_[mdl.eye.icb.igg_exg.C_bound_E6a], sol_mtk_[mdl.eye.icb.igg_exg.C_bound2_E6a], sol_mtk_[mdl.eye.icb.igg_exg.C_E7b], sol_mtk_[mdl.eye.icb.igg_exg.C_bound_E7b], sol_mtk_[mdl.eye.icb.igg_exg.C_bound2_E7b], 
    sol_mtk_[mdl.eye.icb.igg_exg.C_ISM], sol_mtk_[mdl.eye.icb.igg_exg.C_bound_ISM_mem], sol_mtk_[mdl.eye.icb.igg_exg.C_bound_ISM], sol_mtk_[mdl.eye.icb.igg_exg.C_bound2_ISM], sol_mtk_[mdl.eye.icb.igg_exg.C_IntS], 
    param_global[mdl.eye.icb.V_V], param_global[mdl.eye.icb.V_IntS], param_global[mdl.eye.icb.Endothelial_Cell_Frac])

    ADC_Retina = ComputeTissueADC(sol_mtk_[mdl.eye.retina.igg_exg.C_V], 
    sol_mtk_[mdl.eye.retina.igg_exg.C_VM], sol_mtk_[mdl.eye.retina.igg_exg.C_bound_VM_mem], sol_mtk_[mdl.eye.retina.igg_exg.C_bound_VM], sol_mtk_[mdl.eye.retina.igg_exg.C_bound2_VM], 
    sol_mtk_[mdl.eye.retina.igg_exg.C_E7], sol_mtk_[mdl.eye.retina.igg_exg.C_bound_E7], sol_mtk_[mdl.eye.retina.igg_exg.C_bound2_E7], sol_mtk_[mdl.eye.retina.igg_exg.C_E6a], 
    sol_mtk_[mdl.eye.retina.igg_exg.C_bound_E6a], sol_mtk_[mdl.eye.retina.igg_exg.C_bound2_E6a], sol_mtk_[mdl.eye.retina.igg_exg.C_E7b], sol_mtk_[mdl.eye.retina.igg_exg.C_bound_E7b], sol_mtk_[mdl.eye.retina.igg_exg.C_bound2_E7b], 
    sol_mtk_[mdl.eye.retina.igg_exg.C_ISM], sol_mtk_[mdl.eye.retina.igg_exg.C_bound_ISM_mem], sol_mtk_[mdl.eye.retina.igg_exg.C_bound_ISM], sol_mtk_[mdl.eye.retina.igg_exg.C_bound2_ISM], sol_mtk_[mdl.eye.retina.igg_exg.C_IntS], 
    param_global[mdl.eye.retina.V_V], param_global[mdl.eye.retina.V_IntS], param_global[mdl.eye.retina.Endothelial_Cell_Frac])

    df = DataFrame(
        # time 
        time_hr = sol_mtk_.t, 
        # average tissue ADC concentration
        la_int = ADC_LI, 
        sm_int = ADC_SI, 
        lung = ADC_Lung, 
        liver = ADC_Liver, 
        marrow = ADC_Marrow, 
        skin = ADC_Skin, 
        choroid = ADC_Choroid, 
        retina = ADC_Retina, 
        icb = ADC_ICB, 
        # raw free ADC concentration 
        ah = sol_mtk_[mdl.eye.ah_igg_exg.C_AQ], 
        vh = sol_mtk_[mdl.eye.vh_igg_exg.C_VH], 
        cornea = sol_mtk_[mdl.eye.igg_exg.C_COR], 
    )

    return df
end
