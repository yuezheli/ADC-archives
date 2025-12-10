# date: 10/17/2025 
# author: Yuezhe Li 
# purpose of this code: compute eye avg ADC 

using ProjectRoot
include(@projectroot("script/helper-compute-tissue-free-pl.jl"));

function UveaAvgPL(sol_mtk_, param_global; 
    # volume default are in rabbit volume, [L]
    V_RETINA = 42E-6, V_ICB = 87.8E-6, V_CHOROID = 28.4E-6, 
    include_retina = false,  # retina is generally not included based on uvea definition, but left this option in just in case
    mdl = pbpk_simple)
    # compute avg ADC concentration in different part of the eye 
    eye_pl = @select(OrganFreePL(sol_mtk_, param_global, mdl), :time_hr, :choroid, :retina, :icb, :ah, :vh, :cornea); 
    # compute avg ADC concentration across ICB, choroid, w/o retina 
    if include_retina
        total_uvea_pl = eye_pl.choroid * V_CHOROID .+ eye_pl.retina * V_RETINA .+ eye_pl.icb * V_ICB; # [umol]
        eye_pl.uvea_pl_conc = total_uvea_pl / (V_RETINA + V_ICB + V_CHOROID)
    else
        total_uvea_pl = eye_pl.choroid * V_CHOROID .+ eye_pl.icb * V_ICB; # [umol]
        eye_pl.uvea_pl_conc = total_uvea_pl / (V_ICB + V_CHOROID)
    end
    return eye_pl
end




