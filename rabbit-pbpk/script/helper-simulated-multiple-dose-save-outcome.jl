# date: 11/10/2025 
# author: Yuezhe Li
# purpose of this script: to define a helper function to simulate to multiple dose while save the outcome 


function helper_simulated_multiple_dose_save_Cmax(dose_mgkg, parameterset; 
    mdl = pbpk_simple, bw_sims = BW_homo, u0_ = u0_infusion, tspan = tspan, dose_unit = "mg/kg", 
    IncludePlasma = false)

    # dosing 
    adc_infusion__ = InfusionCallback(dose_mgkg, mdl, BW = bw_sims);
    @time prob_mtk__ = ODEProblem(mdl, merge(Dict(u0_), parameterset), tspan, callback = adc_infusion__); 

    # simulation 
    @time sol_mtk_infusion__ = solve(prob_mtk__, alg=CVODE_BDF());

    # Cmax 
    organ_cmax_cvag_adc_pl__ = leftjoin(
    cmax_cavg_adc(sol_mtk_infusion__, parameterset, mdl, IncludePlasma = IncludePlasma), 
    cmax_cavg_payload(sol_mtk_infusion__, parameterset,mdl, IncludePlasma = IncludePlasma), 
    on=:organ);

    organ_cmax_cvag_adc_pl__.Dose .= string(dose_mgkg) *" " * dose_unit;

    return sol_mtk_infusion__, organ_cmax_cvag_adc_pl__
end

