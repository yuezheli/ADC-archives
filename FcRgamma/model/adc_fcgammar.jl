# date: May 24, 2024
# author: Yuezhe Li 
# purpose of this code: to define a simplified cell model, with ADC uptake into the cell in a fixed rate, fixed degredation rate (through lysosome), and payload-induced killing 

using DifferentialEquations, ComponentArrays
using Parameters: @unpack

const Vc = 4/3 * pi * (5E-6)^3 * 1E3 ; # [L]

function adc_fcgammar!(du,u,p,t) 
    @unpack A_m, A_c, Nc, P_c, P_m = u
    @unpack tdouble, DAR, Emax_Payload, ic50_pl, k_in_FcgR, k_lys, k_in, k_out, k_PL, V_medium, Nmax = p

    flux_ADC_intake = k_in_FcgR * A_m/V_medium
    flux_ADC_lys = k_lys * A_c 

    C_P_c = P_c/(Nc * Vc)
    C_P_m = P_m/V_medium

    Kgrow = log(2)/tdouble * (1 - Nc/Nmax)
    Kkill_eff = Emax_Payload*C_P_c/(ic50_pl + C_P_c)

    du.A_m = -flux_ADC_intake * V_medium
    du.A_c = flux_ADC_intake * V_medium - flux_ADC_lys
    du.P_c = DAR*flux_ADC_lys*(1-Kkill_eff) - k_out * (C_P_c - C_P_m) * (Nc * Vc)
    du.P_m = k_out * (C_P_c - C_P_m) * (Nc * Vc) - k_PL * P_m
    du.Nc = Kgrow * Nc - Kkill_eff * Nc
end

