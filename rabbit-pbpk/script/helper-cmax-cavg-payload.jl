# date: 10/6/2025 
# author: Yuezhe Li 
# purpose of this code: postprocessing payload concentration in organ average or interstitium

function cmax_cavg_payload(sol_mtk_infusion_, param_sims, mdl; IncludePlasma = false)
    # compute PL concentration in organ 
    organ_avg_pl = OrganFreePL(sol_mtk_infusion_, param_sims, mdl);

    # Cmax, PL 
    organ_cmax_avg_pl = DataFrame(
        organ = names(select(organ_avg_pl, Not(:time_hr))), 
        Cmax_organ_avg_pl = map(maximum, eachcol(select(organ_avg_pl, Not(:time_hr))))
    );

    # Cavg, PL 
    organ_cavg_avg_pl = DataFrame(
        organ = names(select(organ_avg_pl, Not(:time_hr))), 
        Cavg_organ_avg_pl = map(mean, eachcol(select(organ_avg_pl, Not(:time_hr))))
    );

    # Cmax, PL, tissue interstitium
    organ_cmax_ints_pl = DataFrame(
        organ = ["la_int", "sm_int", "lung", "liver", "marrow", "skin", "choroid", "retina", "icb", "ah", "vh", "cornea"], 
        Cmax_organ_ints_pl = [
            maximum(sol_mtk_infusion_[mdl.la_int.PL_tissue.C_PL_IntS]),
            maximum(sol_mtk_infusion_[mdl.sm_int.PL_tissue.C_PL_IntS]),
            maximum(sol_mtk_infusion_[mdl.lung.PL_tissue.C_PL_IntS]),
            maximum(sol_mtk_infusion_[mdl.liver.PL_tissue.C_PL_IntS]),
            maximum(sol_mtk_infusion_[mdl.marrow.PL_tissue.C_PL_IntS]),
            maximum(sol_mtk_infusion_[mdl.skin.PL_tissue.C_PL_IntS]),
            maximum(sol_mtk_infusion_[mdl.eye.choroid.PL_tissue.C_PL_IntS]),
            maximum(sol_mtk_infusion_[mdl.eye.retina.PL_tissue.C_PL_IntS]),
            maximum(sol_mtk_infusion_[mdl.eye.icb.PL_tissue.C_PL_IntS]),
            maximum(sol_mtk_infusion_[mdl.eye.pl_cor_humor.C_PL_AH]),
            maximum(sol_mtk_infusion_[mdl.eye.pl_cor_humor.C_PL_VH]),
            maximum(sol_mtk_infusion_[mdl.eye.pl_cor_humor.C_PL_COR]),
        ]
    );

    # Cavg, PL, tissue interstitium
    organ_cavg_ints_pl = DataFrame(
        organ = ["la_int", "sm_int", "lung", "liver", "marrow", "skin", "choroid", "retina", "icb", "ah", "vh", "cornea"], 
        Cavg_organ_ints_pl = [
            mean(sol_mtk_infusion_[mdl.la_int.PL_tissue.C_PL_IntS]),
            mean(sol_mtk_infusion_[mdl.sm_int.PL_tissue.C_PL_IntS]),
            mean(sol_mtk_infusion_[mdl.lung.PL_tissue.C_PL_IntS]),
            mean(sol_mtk_infusion_[mdl.liver.PL_tissue.C_PL_IntS]),
            mean(sol_mtk_infusion_[mdl.marrow.PL_tissue.C_PL_IntS]),
            mean(sol_mtk_infusion_[mdl.skin.PL_tissue.C_PL_IntS]),
            mean(sol_mtk_infusion_[mdl.eye.choroid.PL_tissue.C_PL_IntS]),
            mean(sol_mtk_infusion_[mdl.eye.retina.PL_tissue.C_PL_IntS]),
            mean(sol_mtk_infusion_[mdl.eye.icb.PL_tissue.C_PL_IntS]),
            mean(sol_mtk_infusion_[mdl.eye.pl_cor_humor.C_PL_AH]),
            mean(sol_mtk_infusion_[mdl.eye.pl_cor_humor.C_PL_VH]),
            mean(sol_mtk_infusion_[mdl.eye.pl_cor_humor.C_PL_COR]),
        ]
    );

    organ_cmax_pl = leftjoin(organ_cmax_ints_pl, organ_cmax_avg_pl, on=:organ);
    organ_cavg_pl = leftjoin(organ_cavg_avg_pl, organ_cavg_ints_pl, on=:organ);
    organ_pl = leftjoin(organ_cmax_pl, organ_cavg_pl, on=:organ);

    if IncludePlasma
        # add plasma information for reporting only purpse 
        tmp_df_plasma = DataFrame(
            organ = ["plasma"], 
            Cmax_organ_ints_pl = maximum(sol_mtk_infusion_[mdl.plasma_pl.C_PL_Plasma]),
            Cmax_organ_avg_pl = maximum(sol_mtk_infusion_[mdl.plasma_pl.C_PL_Plasma]), 
            Cavg_organ_avg_pl = mean(sol_mtk_infusion_[mdl.plasma_pl.C_PL_Plasma]),
            Cavg_organ_ints_pl = mean(sol_mtk_infusion_[mdl.plasma_pl.C_PL_Plasma]),
        ); 
        full_pl = vcat(organ_pl, tmp_df_plasma)
        return full_pl
    else 
        sort!(organ_pl, :Cmax_organ_ints_pl, rev = true)
        return organ_pl
    end
end
