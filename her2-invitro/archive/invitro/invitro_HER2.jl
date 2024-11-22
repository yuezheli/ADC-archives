# author: Yuezhe Li 
# date: 3/18/24
# purpose of this script: to create an in vitro model for HER2+ cells

const N_av = 6.0221409e+23 # Avogadro constant 
const Vc = 4/3 * pi * (5E-6)^3 * 1E3 ; # [L]
const MW = 15E4;  # molecular weight, g/mol

using ComponentArrays
using Parameters: @unpack

function invitro_model_her2!(du,u,p,t) 
    @unpack Nc_1, Nc_2, Nc_1_neg, Nc_2_neg, R_s, R_e, AR_s, AR_e, P_c, P_m, P_neg_c, A_m = u

    @unpack tdouble, tdouble_neg, tau, DAR, Emax_Payload, ic50_pl, k_PL, Rcopies, k_deg, k_rec, k_endo, Kon, Kd, k_out, k_in, k_lys, V_medium, Nmax = p
    
    Ntot = Nc_1 + Nc_2 
    Ntot_neg = Nc_1_neg + Nc_2_neg 
    C_P_c = P_c/Vc                  # payload conc inside target+ cells
    C_P_neg_c = P_neg_c/ Vc         # payload conc inside target- cells
    C_P_m = P_m/V_medium

    Kgrow = log(2)/tdouble * (1 - (Ntot + Ntot_neg)/Nmax)
    Koff = Kd*Kon
   
    Ksyn = Rcopies/N_av*k_deg*1e6 # unit in umol/hr 

    flux_A_R_s_binding = A_m*R_s*Kon/V_medium
    flux_AR_s_unbinding = AR_s*Koff
    
    flux_R_s_syn = Ksyn
    flux_R_s_degrade = R_s*k_deg
    flux_R_e_rec = R_e*k_rec
    flux_R_s_int = R_s*k_endo
    
    flux_AR_s_int = AR_s*k_endo
    flux_AR_e_rec = AR_e*k_rec
    
    flux_AR_cat = AR_e*k_lys
    flux_P_l_to_P_c = DAR*AR_e*k_lys

    flux_P_c_to_P_m = k_out*(C_P_c - C_P_m)*Vc
    flux_P_neg_c_from_P_m = k_in * (C_P_m - C_P_neg_c)*Vc
    Kkill_eff = Emax_Payload*C_P_c/(ic50_pl + C_P_c)

    du.A_m = flux_AR_s_unbinding * Nc_1 - flux_A_R_s_binding * Nc_1

    du.R_s = flux_R_s_syn - flux_R_s_degrade + flux_R_e_rec - flux_R_s_int + flux_AR_s_unbinding - flux_A_R_s_binding 
    du.R_e = flux_R_s_int - flux_R_e_rec 

    du.AR_s = flux_A_R_s_binding - flux_AR_s_unbinding - flux_AR_s_int + flux_AR_e_rec
    du.AR_e = flux_AR_s_int - flux_AR_e_rec - flux_AR_cat 

    du.Nc_1 = Kgrow*Nc_1 - Kkill_eff*Nc_1
    du.Nc_2 = Kkill_eff*Nc_1 - Nc_2/tau
    
    du.P_c = (1-Kkill_eff)*flux_P_l_to_P_c - flux_P_c_to_P_m 
    
    # warhead in medium
    Du_P_m = flux_P_c_to_P_m * Nc_1 + DAR*1/tau*Nc_2*(AR_e + AR_s) - flux_P_neg_c_from_P_m * Nc_1_neg
    du.P_m = Du_P_m - k_PL * P_m; 

    # receptor-negative cells
    Kgrow_neg = log(2)/ tdouble_neg * (1 - (Ntot + Ntot_neg)/Nmax)
    Kkill_eff_neg = Emax_Payload*C_P_neg_c/(ic50_pl + C_P_neg_c)

    du.Nc_1_neg = Kgrow_neg*Nc_1_neg - Kkill_eff_neg*Nc_1_neg
    du.Nc_2_neg = Kkill_eff_neg*Nc_1_neg - Nc_2_neg/tau
    du.P_neg_c = flux_P_neg_c_from_P_m

end



