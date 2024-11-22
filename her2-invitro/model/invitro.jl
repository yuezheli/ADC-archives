# date: 10/23/2024
# author: Yuezhe Li 
# purpose of this code: to code a generic tumor model 

using Parameters: @unpack

include("adc-constants.jl");

function invitro_model_her2!(du,u,p,t) 
    @unpack Nc_1, Nc_2, Nc_1_neg, Nc_2_neg, R_s, R_e, AR_s, AR_e, P_c_pos, P_m, P_c_neg, A_m, A_e = u # [uM or cell count]

    @unpack tdouble_pos, tdouble_neg, tau, DAR, k_kill_max, ic50_pl, k_PL, k_PL_ex, Rcopies, k_deg, k_rec, k_endo, k_on, Kd, k_out, k_in, k_lys, V_medium, Nmax = p
    
    Ntot = Nc_1 + Nc_2 + Nc_1_neg + Nc_2_neg 
    Vpos = Vc * Nc_1; 
    Vneg = Vc * Nc_1_neg; 
    
    k_grow_pos = log(2)/tdouble_pos * (1 - Ntot/Nmax)
    k_grow_neg = log(2)/tdouble_neg * (1 - Ntot/Nmax)
    k_off = Kd*k_on
    k_taa_syn = Rcopies/N_av*k_deg*1e6 * Nc_1 / V_medium; # [uM/hr]
   
    flux_A_R_s_binding = A_m * R_s * k_on
    flux_AR_s_unbinding = AR_s * k_off

    flux_A_R_e_binding = A_e * R_e * k_on
    flux_AR_e_unbinding = AR_e * k_off
    
    k_kill_pos = k_kill_max*P_c_pos/(ic50_pl + P_c_pos)
    k_kill_neg = k_kill_max*P_c_neg/(ic50_pl + P_c_neg)

    du.A_m = flux_AR_s_unbinding - flux_A_R_s_binding

    du.R_s = k_taa_syn - R_s*k_deg + R_e*k_rec*Vpos/V_medium - R_s*k_endo + flux_AR_s_unbinding - flux_A_R_s_binding - k_kill_pos*R_s
    du.R_e = R_s*k_endo*V_medium/Vpos - R_e*k_rec - k_kill_pos*R_e - flux_A_R_e_binding

    du.AR_s = flux_A_R_s_binding - flux_AR_s_unbinding - AR_s*k_endo + AR_e*k_rec*Vpos/V_medium - k_kill_pos*AR_s
    du.AR_e = AR_s*k_endo*V_medium/Vpos - AR_e*k_rec - k_lys*AR_e - k_kill_pos*AR_e - flux_AR_e_unbinding + flux_A_R_e_binding

    du.A_e = -flux_A_R_e_binding + flux_AR_e_unbinding - k_kill_pos*A_e - k_lys*A_e

    du.Nc_1 = k_grow_pos*Nc_1 - k_kill_pos*Nc_1
    du.Nc_2 = k_kill_pos*Nc_1 - Nc_2/tau
    
    du.P_c_pos = DAR*(k_lys*AR_e + k_lys*A_e) - k_out*(P_c_pos - P_m) - k_kill_pos*P_c_pos - k_PL_ex * P_c_pos
    du.P_m = k_out*(P_c_pos - P_m)*Vpos/V_medium - k_in*(P_m - P_c_neg)*Vneg/V_medium + Nc_2/tau*Vc/V_medium*P_c_pos - k_PL*P_m
    
    # receptor-negative cells
    du.Nc_1_neg = k_grow_neg*Nc_1_neg - k_kill_neg*Nc_1_neg
    du.Nc_2_neg = k_kill_neg*Nc_1_neg - Nc_2_neg/tau
    du.P_c_neg = k_in*(P_m - P_c_neg) - k_kill_neg*P_c_neg

end

