# date: 10/3/2025 
# author: Yuezhe Li 
# purpose of this code: to compute total free payload concentration in a organ 

function ComputeTissueFreePayload(C_plasma_PL, C_endo_PL, C_IntS_PL, V_V, V_IntS, Endothelial_Cell_Frac; 
    pino_time = 10.8/60, CL_up_in_nL_per_hour_per_million_cells = 150, Total_Endothelial_Cell = 1.422e9, Scale_Factor = 603.7)

    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    V_endosomal = CL_up_in_nL_per_hour_per_million_cells * pino_time * Organ_Endothelial_Cell * 1e-6 * 1e-9

    total_PL_mass =  C_plasma_PL * V_V + C_endo_PL * V_endosomal + C_IntS_PL * V_IntS

    V_Organ = V_V + V_endosomal + V_IntS

    return total_PL_mass/V_Organ
end

function OrganFreePL(sol_mtk_, param_global, mdl = pbpk_simple)

    PL_LI = ComputeTissueFreePayload(sol_mtk_[mdl.plasma_pl.C_PL_Plasma], sol_mtk_[mdl.la_int.PL_tissue.C_PL_endo], sol_mtk_[mdl.la_int.PL_tissue.C_PL_IntS], 
    param_global[mdl.la_int.V_V], param_global[mdl.la_int.V_IntS], param_global[mdl.la_int.Endothelial_Cell_Frac])

    PL_SI = ComputeTissueFreePayload(sol_mtk_[mdl.plasma_pl.C_PL_Plasma], sol_mtk_[mdl.sm_int.PL_tissue.C_PL_endo], sol_mtk_[mdl.sm_int.PL_tissue.C_PL_IntS], 
    param_global[mdl.sm_int.V_V], param_global[mdl.sm_int.V_IntS], param_global[mdl.sm_int.Endothelial_Cell_Frac])

    PL_Lung = ComputeTissueFreePayload(sol_mtk_[mdl.plasma_pl.C_PL_Plasma], sol_mtk_[mdl.lung.PL_tissue.C_PL_endo], sol_mtk_[mdl.lung.PL_tissue.C_PL_IntS], 
    param_global[mdl.lung.V_V], param_global[mdl.lung.V_IntS], param_global[mdl.lung.Endothelial_Cell_Frac])

    PL_Liver = ComputeTissueFreePayload(sol_mtk_[mdl.plasma_pl.C_PL_Plasma], sol_mtk_[mdl.liver.PL_tissue.C_PL_endo], sol_mtk_[mdl.liver.PL_tissue.C_PL_IntS], 
    param_global[mdl.liver.V_V], param_global[mdl.liver.V_IntS], param_global[mdl.liver.Endothelial_Cell_Frac])

    PL_Marrow = ComputeTissueFreePayload(sol_mtk_[mdl.plasma_pl.C_PL_Plasma], sol_mtk_[mdl.marrow.PL_tissue.C_PL_endo], sol_mtk_[mdl.marrow.PL_tissue.C_PL_IntS], 
    param_global[mdl.marrow.V_V], param_global[mdl.marrow.V_IntS], param_global[mdl.marrow.Endothelial_Cell_Frac])

    PL_Skin = ComputeTissueFreePayload(sol_mtk_[mdl.plasma_pl.C_PL_Plasma], sol_mtk_[mdl.skin.PL_tissue.C_PL_endo], sol_mtk_[mdl.skin.PL_tissue.C_PL_IntS], 
    param_global[mdl.skin.V_V], param_global[mdl.skin.V_IntS], param_global[mdl.skin.Endothelial_Cell_Frac])

    PL_choroid = ComputeTissueFreePayload(sol_mtk_[mdl.plasma_pl.C_PL_Plasma], sol_mtk_[mdl.eye.choroid.PL_tissue.C_PL_endo], sol_mtk_[mdl.eye.choroid.PL_tissue.C_PL_IntS], 
    param_global[mdl.eye.choroid.V_V], param_global[mdl.eye.choroid.V_IntS], param_global[mdl.eye.choroid.Endothelial_Cell_Frac])

    PL_retina = ComputeTissueFreePayload(sol_mtk_[mdl.plasma_pl.C_PL_Plasma], sol_mtk_[mdl.eye.retina.PL_tissue.C_PL_endo], sol_mtk_[mdl.eye.retina.PL_tissue.C_PL_IntS], 
    param_global[mdl.eye.retina.V_V], param_global[mdl.eye.retina.V_IntS], param_global[mdl.eye.retina.Endothelial_Cell_Frac])

    PL_icb = ComputeTissueFreePayload(sol_mtk_[mdl.plasma_pl.C_PL_Plasma], sol_mtk_[mdl.eye.icb.PL_tissue.C_PL_endo], sol_mtk_[mdl.eye.icb.PL_tissue.C_PL_IntS], 
    param_global[mdl.eye.icb.V_V], param_global[mdl.eye.icb.V_IntS], param_global[mdl.eye.icb.Endothelial_Cell_Frac])


    df = DataFrame(
        # time 
        time_hr = sol_mtk_.t, 
        # average tissue free payload concentration
        la_int = PL_LI, 
        sm_int = PL_SI, 
        lung = PL_Lung, 
        liver = PL_Liver, 
        marrow = PL_Marrow, 
        skin = PL_Skin, 
        choroid = PL_choroid, 
        retina = PL_retina, 
        icb = PL_icb, 
        # raw free payload concentration 
        ah = sol_mtk_[mdl.eye.pl_cor_humor.C_PL_AH], 
        vh = sol_mtk_[mdl.eye.pl_cor_humor.C_PL_VH], 
        cornea = sol_mtk_[mdl.eye.pl_cor_humor.C_PL_COR], 
    )

    return df
end
