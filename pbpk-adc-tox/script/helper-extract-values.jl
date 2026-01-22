# date: 1/14/2025 
# author: Yuezhe Li 
# purpose of this code: to extract ADC/ PL concentration in all the organs 

function extract_adc_tissue_interstitia(sol_dict, mdl, dose)
    tmp_df = DataFrame(
        time_hr = sol_dict.t, 
        lung = sol_dict[mdl.lung.igg_exg.C_IntS], 
        liver = sol_dict[mdl.liver.igg_exg.C_IntS], 
        heart = sol_dict[mdl.heart.igg_exg.C_IntS], 
        muscle = sol_dict[mdl.muscle.igg_exg.C_IntS], 
        skin = sol_dict[mdl.skin.igg_exg.C_IntS], 
        adipose = sol_dict[mdl.adipose.igg_exg.C_IntS], 
        bone = sol_dict[mdl.bone.igg_exg.C_IntS], 
        brain = sol_dict[mdl.brain.igg_exg.C_IntS], 
        kidney = sol_dict[mdl.kidney.igg_exg.C_IntS], 
        sm_int = sol_dict[mdl.sm_int.igg_exg.C_IntS], 
        la_int = sol_dict[mdl.la_int.igg_exg.C_IntS], 
        pancreas = sol_dict[mdl.pancreas.igg_exg.C_IntS], 
        thymus = sol_dict[mdl.thymus.igg_exg.C_IntS], 
        spleen = sol_dict[mdl.spleen.igg_exg.C_IntS], 
        other = sol_dict[mdl.other.igg_exg.C_IntS], 
        marrow = sol_dict[mdl.marrow.igg_exg.C_IntS], 
        choroid = sol_dict[mdl.eye.choroid.igg_exg.C_IntS], 
        icb = sol_dict[mdl.eye.icb.igg_exg.C_IntS], 
        retina = sol_dict[mdl.eye.retina.igg_exg.C_IntS], 
        cornea = sol_dict[mdl.eye.igg_exg.C_COR], 
        ah = sol_dict[mdl.eye.ah_igg_exg.C_AQ], 
        vh = sol_dict[mdl.eye.vh_igg_exg.C_VH], 
        plasma = sol_dict[mdl.plasma_exg.C_Plasma], 
        ln = sol_dict[mdl.plasma_exg.C_LN]
        );

    tmp_df.Dose .= dose
    return tmp_df
end


function extract_PL_tissue_interstitia(sol_dict, mdl, dose)
    tmp_df = DataFrame(
        time_hr = sol_dict.t, 
        lung = sol_dict[mdl.lung.PL_tissue.C_PL_IntS], 
        liver = sol_dict[mdl.liver.PL_tissue.C_PL_IntS], 
        heart = sol_dict[mdl.heart.PL_tissue.C_PL_IntS], 
        muscle = sol_dict[mdl.muscle.PL_tissue.C_PL_IntS], 
        skin = sol_dict[mdl.skin.PL_tissue.C_PL_IntS], 
        adipose = sol_dict[mdl.adipose.PL_tissue.C_PL_IntS], 
        bone = sol_dict[mdl.bone.PL_tissue.C_PL_IntS], 
        brain = sol_dict[mdl.brain.PL_tissue.C_PL_IntS], 
        kidney = sol_dict[mdl.kidney.PL_tissue.C_PL_IntS], 
        sm_int = sol_dict[mdl.sm_int.PL_tissue.C_PL_IntS], 
        la_int = sol_dict[mdl.la_int.PL_tissue.C_PL_IntS], 
        pancreas = sol_dict[mdl.pancreas.PL_tissue.C_PL_IntS], 
        thymus = sol_dict[mdl.thymus.PL_tissue.C_PL_IntS], 
        spleen = sol_dict[mdl.spleen.PL_tissue.C_PL_IntS], 
        other = sol_dict[mdl.other.PL_tissue.C_PL_IntS], 
        marrow = sol_dict[mdl.marrow.PL_tissue.C_PL_IntS], 
        choroid = sol_dict[mdl.eye.choroid.PL_tissue.C_PL_IntS], 
        icb = sol_dict[mdl.eye.icb.PL_tissue.C_PL_IntS], 
        retina = sol_dict[mdl.eye.retina.PL_tissue.C_PL_IntS], 
        cornea = sol_dict[mdl.eye.pl_cor_humor.C_PL_COR], 
        ah = sol_dict[mdl.eye.pl_cor_humor.C_PL_AH], 
        vh = sol_dict[mdl.eye.pl_cor_humor.C_PL_VH], 
        plasma = sol_dict[mdl.plasma_pl.C_PL_Plasma], 
        );

    tmp_df.Dose .= dose
    return tmp_df
end


function extract_PL_tissue_endothelial_cells(sol_dict, mdl, dose)
    tmp_df = DataFrame(
        time_hr = sol_dict.t, 
        lung = sol_dict[mdl.lung.PL_tissue.C_PL_endo], 
        liver = sol_dict[mdl.liver.PL_tissue.C_PL_endo], 
        heart = sol_dict[mdl.heart.PL_tissue.C_PL_endo], 
        muscle = sol_dict[mdl.muscle.PL_tissue.C_PL_endo], 
        skin = sol_dict[mdl.skin.PL_tissue.C_PL_endo], 
        adipose = sol_dict[mdl.adipose.PL_tissue.C_PL_endo], 
        bone = sol_dict[mdl.bone.PL_tissue.C_PL_endo], 
        brain = sol_dict[mdl.brain.PL_tissue.C_PL_endo], 
        kidney = sol_dict[mdl.kidney.PL_tissue.C_PL_endo], 
        sm_int = sol_dict[mdl.sm_int.PL_tissue.C_PL_endo], 
        la_int = sol_dict[mdl.la_int.PL_tissue.C_PL_endo], 
        pancreas = sol_dict[mdl.pancreas.PL_tissue.C_PL_endo], 
        thymus = sol_dict[mdl.thymus.PL_tissue.C_PL_endo], 
        spleen = sol_dict[mdl.spleen.PL_tissue.C_PL_endo], 
        other = sol_dict[mdl.other.PL_tissue.C_PL_endo], 
        marrow = sol_dict[mdl.marrow.PL_tissue.C_PL_endo], 
        choroid = sol_dict[mdl.eye.choroid.PL_tissue.C_PL_endo], 
        icb = sol_dict[mdl.eye.icb.PL_tissue.C_PL_endo], 
        retina = sol_dict[mdl.eye.retina.PL_tissue.C_PL_endo], 
        cornea = sol_dict[mdl.eye.pl_cor_humor.C_PL_COR], 
        ah = sol_dict[mdl.eye.pl_cor_humor.C_PL_AH], 
        vh = sol_dict[mdl.eye.pl_cor_humor.C_PL_VH], 
        plasma = sol_dict[mdl.plasma_pl.C_PL_Plasma], 
        );

    tmp_df.Dose .= dose
    return tmp_df
end






