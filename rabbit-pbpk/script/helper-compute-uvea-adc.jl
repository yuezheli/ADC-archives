# date: 10/16/2025 
# author: Yuezhe Li 
# purpose of this code: to compute the average ADC concentration in non-humor, not cornea part of the eye (i.e., choroid, iris-ciliary body, retina)

using ProjectRoot
include(@projectroot("script/helper-compute-tissue-adc.jl"));

function UveaAvgADC(sol_mtk_, param_global; 
    # volume default are in rabbit volume, [L]
    V_RETINA = 42E-6, V_ICB = 87.8E-6, V_CHOROID = 28.4E-6, 
    include_retina = false,  # retina is generally not included based on uvea definition, but left this option in just in case
    mdl = pbpk_simple)
    # compute avg ADC concentration in different part of the eye 
    eye_adc = @select(OrganAvgADC(sol_mtk_, param_global, mdl), :time_hr, :choroid, :retina, :icb, :ah, :vh, :cornea); 
    # compute avg ADC concentration across ICB, choroid, w/o retina 
    if include_retina
        total_uvea_adc = eye_adc.choroid * V_CHOROID .+ eye_adc.retina * V_RETINA .+ eye_adc.icb * V_ICB; # [umol]
        eye_adc.uvea_adc_conc = total_uvea_adc / (V_RETINA + V_ICB + V_CHOROID)
    else
        total_uvea_adc = eye_adc.choroid * V_CHOROID .+ eye_adc.icb * V_ICB; # [umol]
        eye_adc.uvea_adc_conc = total_uvea_adc / (V_ICB + V_CHOROID)
    end
    return eye_adc
end

