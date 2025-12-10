# date: 10/6/2025 
# author: Yuezhe Li 
# purpose of this code: to postprocessing ADC concentrations (organ average, intersititum)

function cmax_cavg_adc(sol_mtk_infusion_, param_sims, mdl; IncludePlasma = false)

    # compute ADC, averaged in organ 
    organ_avg_adc = OrganAvgADC(sol_mtk_infusion_, param_sims, mdl);

    # compute Cmax of avg organ ADC
    organ_cmax_avg_adc = DataFrame(
        organ = names(select(organ_avg_adc, Not(:time_hr))), 
        Cmax_organ_avg_adc = map(maximum, eachcol(select(organ_avg_adc, Not(:time_hr))))
    );

    # compute Cavg of avg organ ADC
    organ_cavg_avg_adc = DataFrame(
        organ = names(select(organ_avg_adc, Not(:time_hr))), 
        Cavg_organ_avg_adc = map(mean, eachcol(select(organ_avg_adc, Not(:time_hr))))
    );

    # compute Cmax of ADC interstitium
    organ_cmax_ints_adc = DataFrame(
        organ = ["la_int", "sm_int", "lung", "liver", "marrow", "skin", "choroid", "retina", "icb", "ah", "vh", "cornea"], 
        Cmax_organ_ints_adc = [
            maximum(sol_mtk_infusion_[mdl.la_int.igg_exg.C_IntS]),
            maximum(sol_mtk_infusion_[mdl.sm_int.igg_exg.C_IntS]),
            maximum(sol_mtk_infusion_[mdl.lung.igg_exg.C_IntS]),
            maximum(sol_mtk_infusion_[mdl.liver.igg_exg.C_IntS]),
            maximum(sol_mtk_infusion_[mdl.marrow.igg_exg.C_IntS]),
            maximum(sol_mtk_infusion_[mdl.skin.igg_exg.C_IntS]),
            maximum(sol_mtk_infusion_[mdl.eye.choroid.igg_exg.C_IntS]),
            maximum(sol_mtk_infusion_[mdl.eye.retina.igg_exg.C_IntS]),
            maximum(sol_mtk_infusion_[mdl.eye.icb.igg_exg.C_IntS]),
            maximum(sol_mtk_infusion_[mdl.eye.ah_igg_exg.C_AQ]),
            maximum(sol_mtk_infusion_[mdl.eye.vh_igg_exg.C_VH]),
            maximum(sol_mtk_infusion_[mdl.eye.igg_exg.C_COR]),
        ]
    );

    # compute Cavg of ADC interstitium
    organ_cavg_ints_adc = DataFrame(
        organ = ["la_int", "sm_int", "lung", "liver", "marrow", "skin", "choroid", "retina", "icb", "ah", "vh", "cornea"], 
        Cavg_organ_ints_adc = [
            mean(sol_mtk_infusion_[mdl.la_int.igg_exg.C_IntS]),
            mean(sol_mtk_infusion_[mdl.sm_int.igg_exg.C_IntS]),
            mean(sol_mtk_infusion_[mdl.lung.igg_exg.C_IntS]),
            mean(sol_mtk_infusion_[mdl.liver.igg_exg.C_IntS]),
            mean(sol_mtk_infusion_[mdl.marrow.igg_exg.C_IntS]),
            mean(sol_mtk_infusion_[mdl.skin.igg_exg.C_IntS]),
            mean(sol_mtk_infusion_[mdl.eye.choroid.igg_exg.C_IntS]),
            mean(sol_mtk_infusion_[mdl.eye.retina.igg_exg.C_IntS]),
            mean(sol_mtk_infusion_[mdl.eye.icb.igg_exg.C_IntS]),
            mean(sol_mtk_infusion_[mdl.eye.ah_igg_exg.C_AQ]),
            mean(sol_mtk_infusion_[mdl.eye.vh_igg_exg.C_VH]),
            mean(sol_mtk_infusion_[mdl.eye.igg_exg.C_COR]),
        ]
    );

    # join the data frames 
    organ_cmax_adc = leftjoin(organ_cmax_ints_adc, organ_cmax_avg_adc, on=:organ);
    organ_cavg_adc = leftjoin(organ_cavg_avg_adc, organ_cavg_ints_adc, on=:organ);
    organ_adc = leftjoin(organ_cmax_adc, organ_cavg_adc, on=:organ);

    if IncludePlasma
        # add plasma information for reporting only purpse 
        tmp_df_plasma = DataFrame(
            organ = ["plasma"], 
            Cmax_organ_ints_adc = maximum(sol_mtk_infusion_[mdl.plasma_exg.C_Plasma]),
            Cmax_organ_avg_adc = maximum(sol_mtk_infusion_[mdl.plasma_exg.C_Plasma]),
            Cavg_organ_avg_adc = mean(sol_mtk_infusion_[mdl.plasma_exg.C_Plasma]),
            Cavg_organ_ints_adc = mean(sol_mtk_infusion_[mdl.plasma_exg.C_Plasma]),
        ); 
        full_adc = vcat(organ_adc, tmp_df_plasma)
        sort!(full_adc, :Cmax_organ_ints_adc, rev = true); 
        return full_adc
    else 
        sort!(organ_adc, :Cmax_organ_ints_adc, rev = true); 
        return organ_adc
    end
end
