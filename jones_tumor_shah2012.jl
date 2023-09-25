# author: Yuezhe LI
# date: 6/9/23
# purpose of the code: code ODE model for ADC distribution and mechanistic detail of anti-tumor effect in the tumor compartment
# the ADC distribution backbone was adapted from Jones et al., 2019
# https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12461
# original implementation from Metrum Research Group BioPBPK site
# https://github.com/metrumresearchgroup/bioPBPK/tree/main/mAb_Jones2019/julia
# tumor model was adopted from Shah, Haddish-Berhane, and Betts, 2012 and Scheuher et al., 2022 with simplification
# https://pubmed.ncbi.nlm.nih.gov/23151991/
# https://europepmc.org/article/ppr/ppr584999
# in addition, soluble HER in plasma dynamics (also from Scheuher et al., 2022) 
# unit: concentratio [uM]; volume [L]; time [h];

include("params.jl") # load this for number of organs in the model 

using Parameters: @unpack

# Preallocate arrays
N_Organs = 15
const V_endosomal = zeros(N_Organs)
const CL_up = zeros(N_Organs)
const V_VM = zeros(N_Organs)
const V_ISM = zeros(N_Organs)
const V_E7 = zeros(N_Organs)
const V_E6a = zeros(N_Organs)
const V_E7b = zeros(N_Organs)
const Organ_Endothelial_Cell = zeros(N_Organs)
const sigma_V = zeros(N_Organs)
const sigma_IS = zeros(N_Organs)
const Endothelial_Cell_Frac = zeros(N_Organs)
const V_V = zeros(N_Organs)
const V_IntS = zeros(N_Organs)
const V_Organ = zeros(N_Organs)
const PLQ = zeros(N_Organs)
const LF = zeros(N_Organs)

function jones_tumor_shah2012!(du,u,p,t) 
    @unpack C_EXG_Plasma, C_sR_plasma, C_sR_EXG, C_EDG_Plasma, C_EDG_LN, C_EDG, C_EXG_LN, C_EXG, C_FcRn_E6a, C_FcRn_E7, C_FcRn_E7b, C_FcRn_ISM, C_FcRn_VM, 
            tumor, DEGprotein, end_endo_payload, end_cyto_payload, ints_payload = u
    
    @unpack Nc_1, Nc_2, Nc_3, Nc_4, Nc_1_neg, Nc_2_neg, Nc_3_neg, Nc_4_neg, R_s, R_e, AR_s, AR_e, P_c, P_m, P_neg_c, A_m  = tumor
    
    @unpack infusion, tumor_cell_radius, tdouble, tau, PS_Score, PS_kd, KD6_WT, DAR, 
            epsilon, P_ADC, D_ADC, Rcap, Rkrogh, thalf_sR, thalf_sR_adc, init_sR, Kkill_scale, Emax_Payload, IC50_Payload, k_PL, 
            Rcopies, k_deg, k_rec, k_endo, Kon, Kd, k_out, k_in, k_lys = p

    # calculate tumor volume 
    Vc = 4/3 * pi * (tumor_cell_radius*1E-6)^3 * 1E3 ; # [L]
    TumorVolume = (Nc_1 + Nc_2 + Nc_3 + Nc_4 + Nc_1_neg + Nc_2_neg + Nc_3_neg + Nc_4_neg) * Vc / 0.375;
        
    # Organ Indices
    Lung = 1
    Liver = 2
    Heart = 3
    Muscle = 4
    Skin = 5
    Adipose = 6
    Bone = 7
    Brain = 8
    Kidney = 9
    SI = 10
    LI = 11
    Pancreas = 12
    Thymus = 13
    Spleen = 14
    Other = 15
    # Tumor = 16

    # mAb Location Indices
    V = 1
    VM = 2
    E7 = 3
    E6a = 4
    E7b = 5
    ISM = 6
    IntS = 7
    bound_VM = 8
    bound_E7 = 9
    bound_E6a = 10
    bound_E7b = 11
    bound_ISM = 12
    bound2_VM = 13
    bound2_E7 = 14
    bound2_E6a = 15
    bound2_E7b = 16
    bound2_ISM = 17
    bound_VM_mem = 18
    bound_ISM_mem = 19


    # Reflection Coefficients
    sigma_V[Lung] = 0.95
    sigma_V[Liver] = 0.85
    sigma_V[Heart] = 0.95
    sigma_V[Muscle] = 0.95
    sigma_V[Skin] = 0.95
    sigma_V[Adipose] = 0.95
    sigma_V[Bone] = 0.85
    sigma_V[Brain] = 0.99
    sigma_V[Kidney] = 0.9
    sigma_V[SI] = 0.9
    sigma_V[LI] = 0.95
    sigma_V[Pancreas] = 0.9
    sigma_V[Thymus] = 0.9
    sigma_V[Spleen] = 0.85
    sigma_V[Other] = 0.95
    # sigma_V[Tumor] = 0.842 # from Shah and Betts, 2011: https://pubmed.ncbi.nlm.nih.gov/22143261/ 


    # Endothelial Cell Fractions
    Endothelial_Cell_Frac[Lung] = 0.0834
    Endothelial_Cell_Frac[Liver] = 0.1877
    Endothelial_Cell_Frac[Heart] = 0.0011
    Endothelial_Cell_Frac[Muscle] = 0.1928
    Endothelial_Cell_Frac[Skin] = 0.0819
    Endothelial_Cell_Frac[Adipose] = 0.0999
    Endothelial_Cell_Frac[Bone] = 0.1478
    Endothelial_Cell_Frac[Brain] = 0.0115
    Endothelial_Cell_Frac[Kidney] = 0.0157
    Endothelial_Cell_Frac[SI] = 0.0121
    Endothelial_Cell_Frac[LI] = 0.0209
    Endothelial_Cell_Frac[Pancreas] = 0.0001
    Endothelial_Cell_Frac[Thymus] = 0.0001
    Endothelial_Cell_Frac[Spleen] = 0.0499
    Endothelial_Cell_Frac[Other] = 0.0951    
    # Endothelial_Cell_Frac[Tumor] = 0.005    # Lindauer et al., 2017 (placeholder)
    
    # Vascular Volumes [L]
    # based on male BW = 71kg
    V_V[Lung] = 55/1000
    V_V[Liver] = 183/1000
    V_V[Heart] = 13.1/1000
    V_V[Muscle] = 662/1000
    V_V[Skin] = 127/1000
    V_V[Adipose] = 148/1000
    V_V[Bone] = 224/1000
    V_V[Brain] = 31.9/1000
    V_V[Kidney] = 18.2/1000
    V_V[SI] = 6.15/1000
    V_V[LI] = 8.74/1000
    V_V[Pancreas] = 5.7/1000
    V_V[Thymus] = 0.353/1000
    V_V[Spleen] = 26.8/1000
    V_V[Other] = 204/1000     
    # V_V[Tumor] = TumorVolume * 7.0/100.0   # Lindauer et al., 2017

    # Interstitial Volumes [L]
    # based on male BW = 71kg
    V_IntS[Lung] = 300/1000
    V_IntS[Liver] = 429/1000
    V_IntS[Heart] = 48.8/1000
    V_IntS[Muscle] = 3910/1000
    V_IntS[Skin] = 1125/1000
    V_IntS[Adipose] = 2289/1000
    V_IntS[Bone] = 1891/1000
    V_IntS[Brain] = 261/1000
    V_IntS[Kidney] = 49.8/1000
    V_IntS[SI] = 67.1/1000
    V_IntS[LI] = 95.3/1000
    V_IntS[Pancreas] = 18/1000
    V_IntS[Thymus] = 1.09/1000
    V_IntS[Spleen] = 44.3/1000
    V_IntS[Other] = 831/1000      
    # V_IntS[Tumor] = TumorVolume * 0.55 # Lindauer et al., 2017

    # Organ Volumes [L]
    # based on male BW = 71kg
    V_Organ[Lung] = 1000/1000 
    V_Organ[Liver] = 2143/1000
    V_Organ[Heart] = 341/1000
    V_Organ[Muscle] = 30078/1000
    V_Organ[Skin] = 3408/1000
    V_Organ[Adipose] = 13465/1000
    V_Organ[Bone] = 10165/1000
    V_Organ[Brain] = 1450/1000
    V_Organ[Kidney] = 332/1000
    V_Organ[SI] = 385/1000
    V_Organ[LI] = 548/1000
    V_Organ[Pancreas] = 104/1000
    V_Organ[Thymus] = 6.41/1000
    V_Organ[Spleen] = 221/1000
    V_Organ[Other] = 4852/1000    
    # V_Organ[Tumor] = TumorVolume

    # scaled linearly based on the BW using ref. BW = 71 kg
    V_Plasma = 3.126; # [L]
    V_LN     = 0.274; # [L]

    # Blood Flows [L/h]
    # scaled with an allometric exponent of 0.75 based on the BW using ref. BW = 71 kg
    # PLQ[Lung] = 181913/1000 
    PLQ[Liver] = 13210/1000
    PLQ[Heart] = 7752/1000
    PLQ[Muscle] = 33469/1000
    PLQ[Skin] = 11626/1000
    PLQ[Adipose] = 11233/1000
    PLQ[Bone] = 2591/1000
    PLQ[Brain] = 21453/1000
    PLQ[Kidney] = 36402/1000
    PLQ[SI] = 12368/1000
    PLQ[LI] = 12867/1000
    PLQ[Pancreas] = 3056/1000
    PLQ[Thymus] = 353/1000
    PLQ[Spleen] = 6343/1000
    PLQ[Other] = 9190/1000              
    # PLQ[Tumor] = tumorPLQ * TumorVolume      # From https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5270293/
    # rebalance PLQ
    PLQ[Lung] = sum(PLQ[2:end])

    # Lymph Flows [L/h]
    LF_PLQ_frac = 0.002
    LF[Lung]     = PLQ[Lung]*LF_PLQ_frac
    LF[Liver]    = (PLQ[Liver] + PLQ[SI]-LF[SI] + PLQ[LI]-LF[LI] + PLQ[Spleen]-LF[Spleen] + PLQ[Pancreas] - LF[Pancreas])*LF_PLQ_frac 
    LF[Heart]    = PLQ[Heart]*LF_PLQ_frac
    LF[Muscle]   = PLQ[Muscle]*LF_PLQ_frac
    LF[Skin]     = PLQ[Skin]*LF_PLQ_frac
    LF[Adipose]  = PLQ[Adipose]*LF_PLQ_frac
    LF[Bone]     = PLQ[Bone]*LF_PLQ_frac
    LF[Brain]    = PLQ[Brain]*LF_PLQ_frac
    LF[Kidney]   = PLQ[Kidney]*LF_PLQ_frac
    LF[SI]       = PLQ[SI]*LF_PLQ_frac
    LF[LI]       = PLQ[LI]*LF_PLQ_frac
    LF[Pancreas] = PLQ[Pancreas]*LF_PLQ_frac
    LF[Thymus]   = PLQ[Thymus]*LF_PLQ_frac
    LF[Spleen]   = PLQ[Spleen]*LF_PLQ_frac
    LF[Other]    = PLQ[Other]*LF_PLQ_frac
    # LF[Tumor]    = PLQ[Tumor]*0.0029
    L_LymphNode = sum(LF) # [L/h]

    # Setup
    PS_a = 1.8051
    PS_b = 0.2624
    if PS_Score > -1.0
        PS_Kd = 10^(exp(PS_a - PS_b*PS_Score)) 
    else
        PS_Kd = PS_kd
    end
    # k_on_PS = 8.06E+07/1E6 # unit uM.h-1 (this is mouse data)
    k_on_PS = 5.59E8/1E6 
    pino_time = 10.8/60  # min->h
    Scale_Factor = 603.7 # scaling factor for human; if mouse, this number should be 1
    Total_Endothelial_Cell = 1.422e+009 # [number]. Total number of endothelial cells in mouse
    CL_up_in_nL_per_hour_per_million_cells = 150
    @. Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    # Organ_Endothelial_Cell[Tumor] = 0.5/100.0 * TumorVolume[1] / (CL_up_in_nL_per_hour_per_million_cells * pino_time * 1e-6*1e-9) # Based on http://link.springer.com/10.1007/s10928-011-9232-2
    @. CL_up  = CL_up_in_nL_per_hour_per_million_cells*Organ_Endothelial_Cell*1E-15
    k_off_PS = (PS_Kd/1000)*k_on_PS

    tau_VM = 1/60  # [h]  
    tau_ISM = 1/60 # [h]

    E6a_Vol_Pct = 0.33 # [-]
    E7b_Vol_Pct = (1-E6a_Vol_Pct)/2 # [-]
    E7_Vol_Pct = (1-E6a_Vol_Pct)/2 # [-]      
        @. V_VM = CL_up * tau_VM
    @. V_ISM = CL_up * tau_ISM
    @. V_endosomal = CL_up_in_nL_per_hour_per_million_cells * pino_time*Organ_Endothelial_Cell*1e-6*1e-9
    # V_endosomal[Tumor] = 0.5/100.0 * TumorVolume[1]
    @. V_E7 = V_endosomal * E7_Vol_Pct
    @. V_E6a = V_endosomal * E6a_Vol_Pct
    @. V_E7b = V_endosomal * E7b_Vol_Pct

    # kon for EXG    
    k_on_6_EXG  = 8.06E+07/1E6  # [1/(uM*h)]
    k_on_7_EXG  = 1.61E+07/5/1E6  # [1/(uM*h)]
    
    # kon for EDG    
    k_on_6_EDG  = 8.06E+07/1E6 # [1/(uM*h)]
    k_on_7_EDG  = 1.61E+07/5/1E6 # [1/(uM*h)]

    # KD6_WT =  700    # [nM]  [Jones original value]
    # KD7_WT =  154077 # [nM]  [Jones original value]

    KD7_WT = KD6_WT*220 # [nM]  # https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12461

    KD_6_EXG = KD6_WT #  #  [nM]. 
    KD_7_EXG = KD7_WT #   #  [nM]. 
    KD_6_EDG = 700    # [nM]   
    KD_7_EDG = 154077 # [nM]

    on_rate_ratio_1st_vs_2nd_binding = 83.7 # [-]

    k_off_7_EDG = (KD_7_EDG/1000)*k_on_7_EDG # [1/h]
    k_off_6_EDG = (KD_6_EDG/1000)*k_on_6_EDG # [1/h] 
    k_off_7_EXG = (KD_7_EXG/1000)*k_on_7_EXG # [1/h] 
    k_off_6_EXG = (KD_6_EXG/1000)*k_on_6_EXG # [1/h] 

    k_on_6_EDG2  = k_on_6_EDG/on_rate_ratio_1st_vs_2nd_binding # [1/(uM*h)] 
    k_on_7_EDG2  = k_on_7_EDG/on_rate_ratio_1st_vs_2nd_binding # [1/(uM*h)]
    k_on_6_EXG2  = k_on_6_EXG/on_rate_ratio_1st_vs_2nd_binding # [1/(uM*h)]
    k_on_7_EXG2  = k_on_7_EXG/on_rate_ratio_1st_vs_2nd_binding # [1/(uM*h)]

    KD_6_EDG2  = KD_6_EDG*on_rate_ratio_1st_vs_2nd_binding # [nM] 
    KD_7_EDG2  = KD_7_EDG*on_rate_ratio_1st_vs_2nd_binding # [nM] 

    KD_6_EXG2 = on_rate_ratio_1st_vs_2nd_binding*KD6_WT # [nM]
    KD_7_EXG2 = on_rate_ratio_1st_vs_2nd_binding*KD7_WT # [nM]

    k_off_7_EDG2 = (KD_7_EDG2/1000)*k_on_7_EDG2 # [1/h] 
    k_off_6_EDG2 = (KD_6_EDG2/1000)*k_on_6_EDG2 # [1/h] 
    k_off_7_EXG2 = (KD_7_EXG2/1000)*k_on_7_EXG2 # [1/h] 
    k_off_6_EXG2 = (KD_6_EXG2/1000)*k_on_6_EXG2 # [1/h] 

    C_Mem = 18.5 # [uM] 


    on_rate_ratio_1st_vs_2nd_binding = 83.7 # [-]

    kint_PS = 0.0380 # [1/h]. Internalization rate for polyspecificity-membrane bound mAb
    Prob_deg = 0.95 # [-] ~ 95% of free mAb inside endosome will be routed to be degraded
    kdeg_FcRn_Ab = log(2)/11.1

    # Prob_deg = 0.98
    # kdeg_FcRn_Ab = log(2.0)/11.2

    FcRn_recycle_fraction = 0.99

    FR    = 0.715

        # ===================================================================================== 
    #  ---------------------- DIFFERENTIAL VARIABLES/EQUATIONS ----------------------------
    # =====================================================================================
        #  Exogenous IgG 

    #  All organs except Lung and Liver at vascular cmpt (below), exogenous and endogenous IgG, free only

        @views @. du.C_EXG[3:N_Organs,V] = ((PLQ[3:N_Organs]*C_EXG[Lung,V]  #  from Lung 
        - (PLQ[3:N_Organs]- LF[3:N_Organs])*C_EXG[3:N_Organs,V]  #  leave to Main Plasma
        - (1-sigma_V[3:N_Organs])*LF[3:N_Organs]*C_EXG[3:N_Organs,V] #  going via Lymph to Interstitial 
        - CL_up[3:N_Organs]*FR*C_EXG[3:N_Organs,V] #  pinocytosis from vascular to vascular membrane
        + CL_up[3:N_Organs]*FR*C_EXG[3:N_Organs,VM] #  exocytosis from vascular membrane to vascular space 
        )/V_V[3:N_Organs]); 
        
        #  vascular space EDG
        @views @. du.C_EDG[3:N_Organs,V] = ((PLQ[3:N_Organs]*C_EDG[Lung,V] #  from Main Plasma 
        -(PLQ[3:N_Organs]- LF[3:N_Organs])*C_EDG[3:N_Organs,V]  #  leave to Main Plasma
        - (1-sigma_V[3:N_Organs])*LF[3:N_Organs]*C_EDG[3:N_Organs,V] #  going via Lymph to Interstitial 
        - CL_up[3:N_Organs]*FR*C_EDG[3:N_Organs,V] #  pinocytosis from vascular to vascular membrane
        + CL_up[3:N_Organs]*FR*C_EDG[3:N_Organs,VM] #  exocytosis from vascular membrane to vascular space 
        )/V_V[3:N_Organs]); 

    #  Liver, vascular cmpt, exogenous free only
    @views du.C_EXG[Liver,V] =
            ((PLQ[Liver]*C_EXG[Lung,V] #  Inlet distributed from Lung
        + (PLQ[Spleen]-LF[Spleen])*C_EXG[Spleen,V] #  Inlet from Spleen
        + (PLQ[Pancreas]-LF[Pancreas])*C_EXG[Pancreas,V] #  Inlet from Pancrease
        + (PLQ[SI]-LF[SI])*C_EXG[SI,V]  #  Inlet from S.I
        + (PLQ[LI]-LF[LI])*C_EXG[LI,V]   #  Inlet from LF.I
        - C_EXG[Liver,V]*(PLQ[Liver]-LF[Liver] + PLQ[Spleen]-LF[Spleen] 
                                    + PLQ[Pancreas]-LF[Pancreas] + PLQ[SI]-LF[SI]+ PLQ[LI]-LF[LI]) #  Outlet
        - (1-sigma_V[Liver])*LF[Liver]*C_EXG[Liver,V]
        - CL_up[Liver]*FR*C_EXG[Liver,V] #  pinocytosis from vascular to vascular membrane
        + CL_up[Liver]*FR*C_EXG[Liver,VM] #  exocytosis from vascular membrane to vascular space 
        )/V_V[Liver]); 

    #  Lung, vascular cmpt, exogenous free only
    @views du.C_EXG[Lung,V] = (((PLQ[Lung]+LF[Lung])*C_EXG_Plasma[1]
                - PLQ[Lung]*C_EXG[Lung,V]
                -(1-sigma_V[Lung])*LF[Lung]*C_EXG[Lung,V]
                - CL_up[Lung]*FR*C_EXG[Lung,V] #  pinocytosis from vascular to vascular membrane
            + CL_up[Lung]*FR*C_EXG[Lung,VM]#  exocytosis from vascular membrane to vascular space
                        )/V_V[Lung]) # 
                        
    #  Liver, vascular cmpt, endogenous free only
    @views du.C_EDG[Liver,V] =
            ((PLQ[Liver]*C_EDG[Lung,V] #  Inlet distributed from Lung
        + (PLQ[Spleen]-LF[Spleen])*C_EDG[Spleen,V] #  Inlet from Spleen
        + (PLQ[Pancreas]-LF[Pancreas])*C_EDG[Pancreas,V] #  Inlet from Pancrease
        + (PLQ[SI]-LF[SI])*C_EDG[SI,V]  #  Inlet from S.I
        + (PLQ[LI]-LF[LI])*C_EDG[LI,V]   #  Inlet from LF.I
        - C_EDG[Liver,V]*(PLQ[Liver]-LF[Liver] + PLQ[Spleen]-LF[Spleen] 
                                + PLQ[Pancreas]-LF[Pancreas] + PLQ[SI]-LF[SI]+ PLQ[LI]-LF[LI]) #  Outlet
        - (1-sigma_V[Liver])*LF[Liver]*C_EDG[Liver,V]
        - CL_up[Liver]*FR*C_EDG[Liver,V] #  pinocytosis from vascular to vascular membrane
        + CL_up[Liver]*FR*C_EDG[Liver,VM] #  exocytosis from vascular membrane to vascular space 
        )/V_V[Liver]) ; 
        
    #  Lung, vascular cmpt, endogenous free only
    @views du.C_EDG[Lung,V] = (((PLQ[Lung]+LF[Lung])*C_EDG_Plasma[1]
                - PLQ[Lung]*C_EDG[Lung,V]
                -(1-sigma_V[Lung])*LF[Lung]*C_EDG[Lung,V]
                - CL_up[Lung]*FR*C_EDG[Lung,V] #  pinocytosis from vascular to vascular membrane
            + CL_up[Lung]*FR*C_EDG[Lung,VM]#  exocytosis from vascular membrane to vascular space 
                        )/V_V[Lung]) ;
    
    #  Organ: All 
    #  Cmpt:  vascular side membrane, endosome(7.4, 6a, 6b), interstitial side membrane, interstitial space
    #  Species: exogenous IgG (free & bound), endogeneous (free & bound) IgG and FcRn   
    #  vascular side membrane 
    @views @. du.C_EXG[1:N_Organs,VM] = 
     ((CL_up[1:N_Organs]*FR*C_EXG[1:N_Organs,V] #  from vascular space
    - CL_up[1:N_Organs]*FR*C_EXG[1:N_Organs,VM] #  to vascular space
    - CL_up[1:N_Organs]*FR*C_EXG[1:N_Organs,VM] #  to endosomal pH=7.4
    + CL_up[1:N_Organs]*FR*C_EXG[1:N_Organs,E7b] #  from endosomal pH=6
    - k_on_7_EXG*C_EXG[1:N_Organs,VM]*C_FcRn_VM[1:N_Organs]*V_VM[1:N_Organs] #  rxn on 
    + k_off_7_EXG*C_EXG[1:N_Organs,bound_VM]*V_VM[1:N_Organs] #  rxn off
    - k_on_PS*C_EXG[1:N_Organs,VM]*C_Mem*V_VM[1:N_Organs] #  reaction of membrane 
    + k_off_PS*C_EXG[1:N_Organs, bound_VM_mem]*V_VM[1:N_Organs]
        )/V_VM[1:N_Organs]) ;
            
    @views @. du.C_EXG[1:N_Organs, bound_VM_mem] = 
    ((k_on_PS*C_EXG[1:N_Organs,VM]*C_Mem*V_VM[1:N_Organs] #  reaction of membrane 
    - k_off_PS*C_EXG[1:N_Organs, bound_VM_mem]*V_VM[1:N_Organs]
    - kint_PS*C_EXG[1:N_Organs, bound_VM_mem]*V_VM[1:N_Organs]
        )/V_VM[1:N_Organs]) ;
            

    #  endosomal pH=7.4
    @views @. du.C_EXG[1:N_Organs,E7] = 
    ((CL_up[1:N_Organs]*FR*C_EXG[1:N_Organs,VM] #  from vascular membrane
    + CL_up[1:N_Organs]*(1-FR)*C_EXG[1:N_Organs,ISM] #  from is membrane
    - CL_up[1:N_Organs]*C_EXG[1:N_Organs,E7] #  to endosomal pH=6
    - k_on_7_EXG*C_EXG[1:N_Organs,E7]*C_FcRn_E7[1:N_Organs]*V_E7[1:N_Organs] #  rxn on 
    + k_off_7_EXG*C_EXG[1:N_Organs,bound_E7]*V_E7[1:N_Organs] #  rxn off
    + kint_PS*C_EXG[1:N_Organs, bound_VM_mem]*V_VM[1:N_Organs] #  from membrane bound due to polyspecificity
    + kint_PS*C_EXG[1:N_Organs, bound_ISM_mem]*V_ISM[1:N_Organs] #  from membrane bound due to polyspecificity
        )/V_E7[1:N_Organs])  ;

    #  Free EXG mAb in E6a
    @views @. du.C_EXG[1:N_Organs,E6a] = 
    ((CL_up[1:N_Organs]*C_EXG[1:N_Organs,E7] #  from E7
    - CL_up[1:N_Organs]*C_EXG[1:N_Organs,E6a] #  move to E7b or being routed to lysosomal for degradation
    - k_on_6_EXG*C_EXG[1:N_Organs,E6a]*C_FcRn_E6a[1:N_Organs]*V_E6a[1:N_Organs] #  rxn on 
    + k_off_6_EXG*C_EXG[1:N_Organs,bound_E6a]*V_E6a[1:N_Organs] #  rxn off
        )/V_E6a[1:N_Organs]) ;
            
    #  endosomal pH=7.4  b 
    @views @. du.C_EXG[1:N_Organs,E7b] = 
    ((CL_up[1:N_Organs]*(1-Prob_deg)*C_EXG[1:N_Organs,E6a]  #  from 6a
    - CL_up[1:N_Organs]*C_EXG[1:N_Organs,E7b] #  leave from E7b to membranes 
    - k_on_7_EXG*C_EXG[1:N_Organs,E7b]*C_FcRn_E7b[1:N_Organs]*V_E7b[1:N_Organs] #  rxn on 
    + k_off_7_EXG*C_EXG[1:N_Organs,bound_E7b]*V_E7b[1:N_Organs] #  rxn off
        )/V_E7b[1:N_Organs]);
            

    #  IS side mem
    @views @. du.C_EXG[1:N_Organs,ISM] = 
    ((CL_up[1:N_Organs]*(1-FR)*C_EXG[1:N_Organs,IntS] #  from IS space 
    - CL_up[1:N_Organs]*(1-FR)*C_EXG[1:N_Organs,ISM] #  to endosomal pH=7.4
    + CL_up[1:N_Organs]*(1-FR)*C_EXG[1:N_Organs,E7b] #  from endosomal pH=7.4 b
    - CL_up[1:N_Organs]*(1-FR)*C_EXG[1:N_Organs,ISM] #  to interstitial space
    - k_on_7_EXG*C_EXG[1:N_Organs,ISM]*C_FcRn_ISM[1:N_Organs]*V_ISM[1:N_Organs] #  rxn on 
    + k_off_7_EXG*C_EXG[1:N_Organs,bound_ISM]*V_ISM[1:N_Organs] #  rxn off
    - k_on_PS*C_EXG[1:N_Organs,ISM]*C_Mem*V_ISM[1:N_Organs] #  reaction of membrane 
    + k_off_PS*C_EXG[1:N_Organs, bound_ISM_mem]*V_ISM[1:N_Organs]
        )/V_ISM[1:N_Organs]) ;
            
    @views @. du.C_EXG[1:N_Organs, bound_ISM_mem] = ((k_on_PS*C_EXG[1:N_Organs,ISM]*C_Mem*V_ISM[1:N_Organs] #  reaction of membrane 
    - k_off_PS*C_EXG[1:N_Organs, bound_ISM_mem]*V_ISM[1:N_Organs]
    - kint_PS*C_EXG[1:N_Organs, bound_ISM_mem]*V_ISM[1:N_Organs]
        )/V_ISM[1:N_Organs]) ;
            
    #  Interstitial space
    @views @. du.C_EXG[1:N_Organs,IntS] =   
    (((1-sigma_V[1:N_Organs])*LF[1:N_Organs]*C_EXG[1:N_Organs,V] #  going from vascular via Lymph 
    - CL_up[1:N_Organs]*(1-FR)*C_EXG[1:N_Organs,IntS] #  pinocytosis 
    + CL_up[1:N_Organs]*(1-FR)*C_EXG[1:N_Organs,ISM] #  exocytosis
    -(1-sigma_IS[1:N_Organs])*LF[1:N_Organs]*C_EXG[1:N_Organs,IntS]
    )/V_IntS[1:N_Organs]); # (DO I NEED TO UPDATE mAb in TUMOR INTERSTITIAL)

    # # # # # # # #  For FcRn-IgG (bounded)
    @views @. du.C_EXG[1:N_Organs,bound_VM] = 
        ((k_on_7_EXG*C_EXG[1:N_Organs,VM]*C_FcRn_VM[1:N_Organs]*V_VM[1:N_Organs] #  rxn on 
    - k_off_7_EXG*C_EXG[1:N_Organs,bound_VM]*V_VM[1:N_Organs] #  rxn off
    - k_on_7_EXG2*C_EXG[1:N_Organs,bound_VM]*C_FcRn_VM[1:N_Organs]*V_VM[1:N_Organs] #  rxn on 
    + k_off_7_EXG2*C_EXG[1:N_Organs,bound2_VM]*V_VM[1:N_Organs] #  rxn off 
    + CL_up[1:N_Organs]*FR*C_EXG[1:N_Organs,bound_E7b]  #  from Endosomal 6b to Vascular membrane
    - CL_up[1:N_Organs]*FR*C_EXG[1:N_Organs,bound_VM] #  to Endosomal 7
    - kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound_VM]*V_VM[1:N_Organs] #  recover FcRn when IgG-FrRn is destroyed
        )/V_VM[1:N_Organs]);
            
        @views @. du.C_EXG[1:N_Organs,bound2_VM] = 
        ((k_on_7_EXG2*C_EXG[1:N_Organs,bound_VM]*C_FcRn_VM[1:N_Organs]*V_VM[1:N_Organs] #  rxn on 
    - k_off_7_EXG2*C_EXG[1:N_Organs,bound2_VM]*V_VM[1:N_Organs] #  rxn off 
    + CL_up[1:N_Organs]*FR*C_EXG[1:N_Organs,bound2_E7b]  #  from Endosomal 6b to Vascular membrane
    - CL_up[1:N_Organs]*FR*C_EXG[1:N_Organs,bound2_VM] #  to Endosomal 7
    - kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound2_VM]*V_VM[1:N_Organs] #  recover FcRn when IgG-FrRn is destroyed
        )/V_VM[1:N_Organs]);
            

    #  endosomal space 
    @views @. du.C_EXG[1:N_Organs,bound_E7] = 
    ((+ k_on_7_EXG*C_EXG[1:N_Organs,E7]*C_FcRn_E7[1:N_Organs]*V_E7[1:N_Organs] #  rxn on 
    - k_off_7_EXG*C_EXG[1:N_Organs,bound_E7]*V_E7[1:N_Organs] #  rxn off
    - k_on_7_EXG2*C_EXG[1:N_Organs,bound_E7]*C_FcRn_E7[1:N_Organs]*V_E7[1:N_Organs] #  rxn on 
    + k_off_7_EXG2*C_EXG[1:N_Organs,bound2_E7]*V_E7[1:N_Organs] #  rxn off 
    - CL_up[1:N_Organs]*C_EXG[1:N_Organs,bound_E7] #  goes to E6
    + CL_up[1:N_Organs]*FR*C_EXG[1:N_Organs,bound_VM] #  from VM 
    + CL_up[1:N_Organs]*(1-FR)*C_EXG[1:N_Organs,bound_ISM]  #  from ISM
        )/V_E7[1:N_Organs]) ;
            

        @views @. du.C_EXG[1:N_Organs,bound2_E7] = 
        ((k_on_7_EXG2*C_EXG[1:N_Organs,bound_E7]*C_FcRn_E7[1:N_Organs]*V_E7[1:N_Organs] #  rxn on 
    - k_off_7_EXG2*C_EXG[1:N_Organs,bound2_E7]*V_E7[1:N_Organs] #  rxn off 
    - CL_up[1:N_Organs]*C_EXG[1:N_Organs,bound2_E7] #  goes to E6
    + CL_up[1:N_Organs]*FR*C_EXG[1:N_Organs,bound2_VM] #  from VM 
    + CL_up[1:N_Organs]*(1-FR)*C_EXG[1:N_Organs,bound2_ISM]   #  from ISM
        )/V_E7[1:N_Organs]) ;
            
    #  endosomal pH=6 a
    @views @. du.C_EXG[1:N_Organs,bound_E6a] = 
    ((+ k_on_6_EXG*C_EXG[1:N_Organs,E6a]*C_FcRn_E6a[1:N_Organs]*V_E6a[1:N_Organs] #  rxn on 
    - k_off_6_EXG*C_EXG[1:N_Organs,bound_E6a]*V_E6a[1:N_Organs] #  rxn off
    - k_on_6_EXG2*C_EXG[1:N_Organs,bound_E6a]*C_FcRn_E6a[1:N_Organs]*V_E6a[1:N_Organs] #  rxn on 
    + k_off_6_EXG2*C_EXG[1:N_Organs,bound2_E6a]*V_E6a[1:N_Organs] #  rxn off 
    + CL_up[1:N_Organs]*C_EXG[1:N_Organs,bound_E7] #  internalization, from pH7 to pH6a
    - CL_up[1:N_Organs]*C_EXG[1:N_Organs,bound_E6a] #  exocytosis, from E6a to E7b
    )/V_E6a[1:N_Organs]) ;
            
    @views @. du.C_EXG[1:N_Organs,bound2_E6a] = 
    ((+ k_on_6_EXG2*C_EXG[1:N_Organs,bound_E6a]*C_FcRn_E6a[1:N_Organs]*V_E6a[1:N_Organs] #  rxn on 
    - k_off_6_EXG2*C_EXG[1:N_Organs,bound2_E6a]*V_E6a[1:N_Organs] #  rxn off 
    + CL_up[1:N_Organs]*C_EXG[1:N_Organs,bound2_E7] #  internalization, from pH7 to pH6a
    - CL_up[1:N_Organs]*C_EXG[1:N_Organs,bound2_E6a] #  exocytosis, from E6a to E7b
    )/V_E6a[1:N_Organs]);
            
    #  endosomal pH=7.4  b
    @views @. du.C_EXG[1:N_Organs, bound_E7b] = 
    ((+ k_on_7_EXG*C_EXG[1:N_Organs, E7b]*C_FcRn_E7b[1:N_Organs]*V_E7b[1:N_Organs] #  rxn on 
    - k_off_7_EXG*C_EXG[1:N_Organs,bound_E7b]*V_E7b[1:N_Organs] #  rxn off
    - k_on_7_EXG2*C_EXG[1:N_Organs, bound_E7b]*C_FcRn_E7b[1:N_Organs]*V_E7b[1:N_Organs] #  rxn on 
    + k_off_7_EXG2*C_EXG[1:N_Organs,bound2_E7b]*V_E7b[1:N_Organs] #  rxn off 
    + CL_up[1:N_Organs]*C_EXG[1:N_Organs,bound_E6a] #  internalization, from E6a to E7b
    - CL_up[1:N_Organs]*C_EXG[1:N_Organs,bound_E7b] #  exocytosis, from E7b to vascular/interstitial membrane
        )/V_E7b[1:N_Organs]);
            
        @views @. du.C_EXG[1:N_Organs, bound2_E7b] = 
    ((+ k_on_7_EXG2*C_EXG[1:N_Organs, bound_E7b]*C_FcRn_E7b[1:N_Organs]*V_E7b[1:N_Organs] #  rxn on 
    - k_off_7_EXG2*C_EXG[1:N_Organs,bound2_E7b]*V_E7b[1:N_Organs] #  rxn off 
    + CL_up[1:N_Organs]*C_EXG[1:N_Organs,bound2_E6a] #  internalization, from E6a to E7b
    - CL_up[1:N_Organs]*C_EXG[1:N_Organs,bound2_E7b] #  exocytosis, from E7b to vascular/interstitial membrane
        )/V_E7b[1:N_Organs]);
            
    #  Interstitial membrane
    @views @. du.C_EXG[1:N_Organs,bound_ISM] = 
        ((k_on_7_EXG*C_EXG[1:N_Organs,ISM]*C_FcRn_ISM[1:N_Organs]*V_ISM[1:N_Organs] #  rxn on  
    - k_off_7_EXG*C_EXG[1:N_Organs,bound_ISM]*V_ISM[1:N_Organs] #  rxn off
    - k_on_7_EXG2*C_EXG[1:N_Organs,bound_ISM]*C_FcRn_ISM[1:N_Organs]*V_ISM[1:N_Organs] #  rxn on 
    + k_off_7_EXG2*C_EXG[1:N_Organs,bound2_ISM]*V_ISM[1:N_Organs] #  rxn off 
    + CL_up[1:N_Organs]*(1-FR)*C_EXG[1:N_Organs,bound_E7b]  #  from Endosomal 6 to interstitial membrane
    - CL_up[1:N_Organs]*(1-FR)*C_EXG[1:N_Organs,bound_ISM] #  to Endosomal 7
    - kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound_ISM]*V_ISM[1:N_Organs]  #  recover FcRn when IgG-FrRn is destroyed
        )/V_ISM[1:N_Organs]);
            
        @views @. du.C_EXG[1:N_Organs,bound2_ISM] = 
    ((+ k_on_7_EXG2*C_EXG[1:N_Organs,bound_ISM]*C_FcRn_ISM[1:N_Organs]*V_ISM[1:N_Organs] #  rxn on 
    - k_off_7_EXG2*C_EXG[1:N_Organs,bound2_ISM]*V_ISM[1:N_Organs] #  rxn off 
    + CL_up[1:N_Organs]*(1-FR)*C_EXG[1:N_Organs,bound2_E7b]  #  from Endosomal 6 to interstitial membrane
    - CL_up[1:N_Organs]*(1-FR)*C_EXG[1:N_Organs,bound2_ISM] #  to Endosomal 7
    - kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound2_ISM]*V_ISM[1:N_Organs]  #  recover FcRn when IgG-FrRn is destroyed
        )/V_ISM[1:N_Organs]);
    
    #  Free Endogenous IgG 
    
    
    #  vascular side membrane 
    @views @. du.C_EDG[1:N_Organs,VM] =  ((CL_up[1:N_Organs]*FR*C_EDG[1:N_Organs,V] #  from vascular space
    - CL_up[1:N_Organs]*FR*C_EDG[1:N_Organs,VM] #  to vascular space
    - CL_up[1:N_Organs]*FR*C_EDG[1:N_Organs,VM] #  to endosomal pH=7.4
    + CL_up[1:N_Organs]*FR*C_EDG[1:N_Organs,E7b] #  from endosomal pH=7.4 b
    - k_on_7_EDG*C_EDG[1:N_Organs,VM]*C_FcRn_VM[1:N_Organs]*V_VM[1:N_Organs] #  rxn on 
    + k_off_7_EDG*C_EDG[1:N_Organs,bound_VM]*V_VM[1:N_Organs] #  rxn off
    )/V_VM[1:N_Organs]);
    
    #  endosomal pH=7.4
    @views @. du.C_EDG[1:N_Organs,E7] =  ((CL_up[1:N_Organs]*FR*C_EDG[1:N_Organs,VM] #  from vascular membrane
    + CL_up[1:N_Organs]*(1-FR)*C_EDG[1:N_Organs,ISM] #  from is membrane
    - CL_up[1:N_Organs]*C_EDG[1:N_Organs,E7] #  to endosomal pH=6
    - k_on_7_EDG*C_EDG[1:N_Organs,E7]*C_FcRn_E7[1:N_Organs]*V_E7[1:N_Organs] #  rxn on 
    + k_off_7_EDG*C_EDG[1:N_Organs,bound_E7]*V_E7[1:N_Organs] #  rxn off
    )/V_E7[1:N_Organs])# 
    
    #  endosomal pH=6.0 a 
    @views @. du.C_EDG[1:N_Organs,E6a] =  ((CL_up[1:N_Organs]*C_EDG[1:N_Organs,E7] #  from E7
    - CL_up[1:N_Organs]*C_EDG[1:N_Organs,E6a] #  move to E7b or being routed to lysosomal for degradation
    - k_on_6_EDG*C_EDG[1:N_Organs,E6a]*C_FcRn_E6a[1:N_Organs]*V_E6a[1:N_Organs] #  rxn on 
    + k_off_6_EDG*C_EDG[1:N_Organs,bound_E6a]*V_E6a[1:N_Organs] #  rxn off
    )/V_E6a[1:N_Organs]); 
    
    #  endosomal pH=7.4  b 
    @views @. du.C_EDG[1:N_Organs,E7b] =  ((CL_up[1:N_Organs]*(1-Prob_deg)*C_EDG[1:N_Organs,E6a]  #  from 6a
    - CL_up[1:N_Organs]*C_EDG[1:N_Organs,E7b] #  leave from E7b to membranes 
    - k_on_7_EDG*C_EDG[1:N_Organs,E7b]*C_FcRn_E7b[1:N_Organs]*V_E7b[1:N_Organs] #  rxn on 
    + k_off_7_EDG*C_EDG[1:N_Organs,bound_E7b]*V_E7b[1:N_Organs] #  rxn off
    )/V_E7b[1:N_Organs]);   
            
    #  IS side mem
    @views @. du.C_EDG[1:N_Organs,ISM] =  ((CL_up[1:N_Organs]*(1-FR)*C_EDG[1:N_Organs,IntS] #  from IS space 
    - CL_up[1:N_Organs]*(1-FR)*C_EDG[1:N_Organs,ISM] #  to endosomal pH=7.4
    + CL_up[1:N_Organs]*(1-FR)*C_EDG[1:N_Organs,E7b] #  from endosomal pH=7.4 b
    - CL_up[1:N_Organs]*(1-FR)*C_EDG[1:N_Organs,ISM] #  to interstitial space
    - k_on_7_EDG*C_EDG[1:N_Organs,ISM]*C_FcRn_ISM[1:N_Organs]*V_ISM[1:N_Organs] #  rxn on 
    + k_off_7_EDG*C_EDG[1:N_Organs,bound_ISM]*V_ISM[1:N_Organs] #  rxn off
    )/V_ISM[1:N_Organs]);

    #  Interstitial space
    @views @. du.C_EDG[1:N_Organs,IntS] =   
    (((1-sigma_V[1:N_Organs])*LF[1:N_Organs]*C_EDG[1:N_Organs,V] #  going from vascular via Lymph 
    - CL_up[1:N_Organs]*(1-FR)*C_EDG[1:N_Organs,IntS] #  pinocytosis 
    + CL_up[1:N_Organs]*(1-FR)*C_EDG[1:N_Organs,ISM] #  exocytosis
    -(1-sigma_IS[1:N_Organs])*LF[1:N_Organs]*C_EDG[1:N_Organs,IntS];
    )/V_IntS[1:N_Organs]) # 
            
    #  ==============================================================================================
    #  Bound Endogenous IgG 
    
    # # # # # # # #  For FcRn-IgG (bounded)
    @views @. du.C_EDG[1:N_Organs,bound_VM] = 
        ((k_on_7_EDG*C_EDG[1:N_Organs,VM]*C_FcRn_VM[1:N_Organs]*V_VM[1:N_Organs] #  rxn on 
    - k_off_7_EDG*C_EDG[1:N_Organs,bound_VM]*V_VM[1:N_Organs] #  rxn off
    - k_on_7_EDG2*C_EDG[1:N_Organs,bound_VM]*C_FcRn_VM[1:N_Organs]*V_VM[1:N_Organs] #  rxn on 
    + k_off_7_EDG2*C_EDG[1:N_Organs,bound2_VM]*V_VM[1:N_Organs] #  rxn off 
    + CL_up[1:N_Organs]*FR*C_EDG[1:N_Organs,bound_E7b]  #  from Endosomal 6b to Vascular membrane
    - CL_up[1:N_Organs]*FR*C_EDG[1:N_Organs,bound_VM] #  to Endosomal 7
    - kdeg_FcRn_Ab*C_EDG[1:N_Organs,bound_VM]*V_VM[1:N_Organs] #  recover FcRn when IgG-FrRn is destroyed
    )/V_VM[1:N_Organs]) ;
    
    @views @. du.C_EDG[1:N_Organs,bound2_VM] = 
        ((k_on_7_EDG2*C_EDG[1:N_Organs,bound_VM]*C_FcRn_VM[1:N_Organs]*V_VM[1:N_Organs] #  rxn on 
    - k_off_7_EDG2*C_EDG[1:N_Organs,bound2_VM]*V_VM[1:N_Organs] #  rxn off 
    + CL_up[1:N_Organs]*FR*C_EDG[1:N_Organs,bound2_E7b]  #  from Endosomal 6b to Vascular membrane
    - CL_up[1:N_Organs]*FR*C_EDG[1:N_Organs,bound2_VM] #  to Endosomal 7
    - kdeg_FcRn_Ab*C_EDG[1:N_Organs,bound2_VM]*V_VM[1:N_Organs] #  recover FcRn when IgG-FrRn is destroyed
    )/V_VM[1:N_Organs]) ;

    #  endosomal space 
    @views @. du.C_EDG[1:N_Organs,bound_E7] = 
    ((+ k_on_7_EDG*C_EDG[1:N_Organs,E7]*C_FcRn_E7[1:N_Organs]*V_E7[1:N_Organs] #  rxn on 
    - k_off_7_EDG*C_EDG[1:N_Organs,bound_E7]*V_E7[1:N_Organs] #  rxn off
    - k_on_7_EDG2*C_EDG[1:N_Organs,bound_E7]*C_FcRn_E7[1:N_Organs]*V_E7[1:N_Organs] #  rxn on 
    + k_off_7_EDG2*C_EDG[1:N_Organs,bound2_E7]*V_E7[1:N_Organs] #  rxn off 
    - CL_up[1:N_Organs]*C_EDG[1:N_Organs,bound_E7] #  goes to E6
    + CL_up[1:N_Organs]*FR*C_EDG[1:N_Organs,bound_VM] #  from VM 
    + CL_up[1:N_Organs]*(1-FR)*C_EDG[1:N_Organs,bound_ISM];  #  from ISM
    )/V_E7[1:N_Organs]) # 
            
    @views @. du.C_EDG[1:N_Organs,bound2_E7] = 
    ((+ k_on_7_EDG2*C_EDG[1:N_Organs,bound_E7]*C_FcRn_E7[1:N_Organs]*V_E7[1:N_Organs] #  rxn on 
    - k_off_7_EDG2*C_EDG[1:N_Organs,bound2_E7]*V_E7[1:N_Organs] #  rxn off 
    - CL_up[1:N_Organs]*C_EDG[1:N_Organs,bound2_E7] #  goes to E6
    + CL_up[1:N_Organs]*FR*C_EDG[1:N_Organs,bound2_VM] #  from VM 
    + CL_up[1:N_Organs]*(1-FR)*C_EDG[1:N_Organs,bound2_ISM];  #  from ISM
    )/V_E7[1:N_Organs]) # 
    
    #  endosomal pH=6 a
    @views @. du.C_EDG[1:N_Organs,bound_E6a] = 
    ((+ k_on_6_EDG*C_EDG[1:N_Organs,E6a]*C_FcRn_E6a[1:N_Organs]*V_E6a[1:N_Organs] #  rxn on 
    - k_off_6_EDG*C_EDG[1:N_Organs,bound_E6a]*V_E6a[1:N_Organs] #  rxn off
    - k_on_6_EDG2*C_EDG[1:N_Organs,bound_E6a]*C_FcRn_E6a[1:N_Organs]*V_E6a[1:N_Organs] #  rxn on 
    + k_off_6_EDG2*C_EDG[1:N_Organs,bound2_E6a]*V_E6a[1:N_Organs] #  rxn off 
    + CL_up[1:N_Organs]*C_EDG[1:N_Organs,bound_E7] #  internalization, from pH7 to pH6a
    - CL_up[1:N_Organs]*C_EDG[1:N_Organs,bound_E6a] #  exocytosis, from E6a to E7b
    )/V_E6a[1:N_Organs]) #  

    @views @. du.C_EDG[1:N_Organs,bound2_E6a] = 
    ((+ k_on_6_EDG2*C_EDG[1:N_Organs,bound_E6a]*C_FcRn_E6a[1:N_Organs]*V_E6a[1:N_Organs] #  rxn on 
    - k_off_6_EDG2*C_EDG[1:N_Organs,bound2_E6a]*V_E6a[1:N_Organs] #  rxn off 
    + CL_up[1:N_Organs]*C_EDG[1:N_Organs,bound2_E7] #  internalization, from pH7 to pH6a
    - CL_up[1:N_Organs]*C_EDG[1:N_Organs,bound2_E6a] #  exocytosis, from E6a to E7b
    )/V_E6a[1:N_Organs]) ; 

    #  endosomal pH=7.4  b
    @views @. du.C_EDG[1:N_Organs,bound_E7b] = 
    ((+ k_on_7_EDG*C_EDG[1:N_Organs,E7b]*C_FcRn_E7b[1:N_Organs]*V_E7b[1:N_Organs] #  rxn on 
    - k_off_7_EDG*C_EDG[1:N_Organs,bound_E7b]*V_E7b[1:N_Organs] #  rxn off
    - k_on_7_EDG2*C_EDG[1:N_Organs, bound_E7b]*C_FcRn_E7b[1:N_Organs]*V_E7b[1:N_Organs] #  rxn on 
    + k_off_7_EDG2*C_EDG[1:N_Organs,bound2_E7b]*V_E7b[1:N_Organs] #  rxn off 
    + CL_up[1:N_Organs]*C_EDG[1:N_Organs,bound_E6a] #  internalization, from E6a to E7b
    - CL_up[1:N_Organs]*C_EDG[1:N_Organs,bound_E7b] #  exocytosis, from E7b to vascular/interstitial membrane
    )/V_E7b[1:N_Organs]) ; 
    
    @views @. du.C_EDG[1:N_Organs,bound2_E7b] = 
    ((+ k_on_7_EDG2*C_EDG[1:N_Organs, bound_E7b]*C_FcRn_E7b[1:N_Organs]*V_E7b[1:N_Organs] #  rxn on 
    - k_off_7_EDG2*C_EDG[1:N_Organs,bound2_E7b]*V_E7b[1:N_Organs] #  rxn off 
    + CL_up[1:N_Organs]*C_EDG[1:N_Organs,bound2_E6a] #  internalization, from E6a to E7b
    - CL_up[1:N_Organs]*C_EDG[1:N_Organs,bound2_E7b] #  exocytosis, from E7b to vascular/interstitial membrane
    )/V_E7b[1:N_Organs]) ; 
    
    #  Interstitial membrane
    @views @. du.C_EDG[1:N_Organs,bound_ISM] = 
    ((+ k_on_7_EDG*C_EDG[1:N_Organs,ISM]*C_FcRn_ISM[1:N_Organs]*V_ISM[1:N_Organs] #  rxn on  
    - k_off_7_EDG*C_EDG[1:N_Organs,bound_ISM]*V_ISM[1:N_Organs] #  rxn off
    - k_on_7_EDG2*C_EDG[1:N_Organs,bound_ISM]*C_FcRn_ISM[1:N_Organs]*V_ISM[1:N_Organs] #  rxn on 
    + k_off_7_EDG2*C_EDG[1:N_Organs,bound2_ISM]*V_ISM[1:N_Organs] #  rxn off 
    + CL_up[1:N_Organs]*(1-FR)*C_EDG[1:N_Organs,bound_E7b]  #  from Endosomal 6 to interstitial membrane
    - CL_up[1:N_Organs]*(1-FR)*C_EDG[1:N_Organs,bound_ISM] #  to Endosomal 7.4
    - kdeg_FcRn_Ab*C_EDG[1:N_Organs,bound_ISM]*V_ISM[1:N_Organs]  #  recover FcRn when IgG-FrRn is destroyed
    )/V_ISM[1:N_Organs]) ;
    
    @views @. du.C_EDG[1:N_Organs,bound2_ISM] = 
    ((+ k_on_7_EDG2*C_EDG[1:N_Organs,bound_ISM]*C_FcRn_ISM[1:N_Organs]*V_ISM[1:N_Organs] #  rxn on 
    - k_off_7_EDG2*C_EDG[1:N_Organs,bound2_ISM]*V_ISM[1:N_Organs] #  rxn off 
    + CL_up[1:N_Organs]*(1-FR)*C_EDG[1:N_Organs,bound2_E7b]  #  from Endosomal 7.4 b to interstitial membrane
    - CL_up[1:N_Organs]*(1-FR)*C_EDG[1:N_Organs,bound2_ISM] #  to Endosomal 7.4
    - kdeg_FcRn_Ab*C_EDG[1:N_Organs,bound2_ISM]*V_ISM[1:N_Organs]  #  recover FcRn when IgG-FrRn is destroyed
    )/V_ISM[1:N_Organs]) ;
    
    #  ==============================================================================================
    #  Free FcRn in Organs

    #  vascular side membrane 
    @views @. du.C_FcRn_VM[1:N_Organs] = 
    ((+ CL_up[1:N_Organs]*FR*C_FcRn_E7b[1:N_Organs]*(1-FcRn_recycle_fraction)# from endosomal pH=7.4  b
    - CL_up[1:N_Organs]*FR*C_FcRn_VM[1:N_Organs] #  to vascular space
    #  EDG 
    +V_VM[1:N_Organs]*(
    - k_on_7_EDG*C_EDG[1:N_Organs,VM]*C_FcRn_VM[1:N_Organs] # rxn on in vascular side membrane 
    + k_off_7_EDG*C_EDG[1:N_Organs,bound_VM] # rxn off vascular side membrane 
    + kdeg_FcRn_Ab*C_EDG[1:N_Organs,bound_VM]
    - k_on_7_EDG2*C_EDG[1:N_Organs,bound_VM]*C_FcRn_VM[1:N_Organs] # rxn on in vascular side membrane 
    + k_off_7_EDG2*C_EDG[1:N_Organs,bound2_VM] # rxn off vascular side membrane 
    + 2*kdeg_FcRn_Ab*C_EDG[1:N_Organs,bound2_VM] # 
    ########## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    #  EXG 
    - k_on_7_EXG*C_EXG[1:N_Organs,VM]*C_FcRn_VM[1:N_Organs] #  rxn on in vascular side membrane 
    + k_off_7_EXG*C_EXG[1:N_Organs,bound_VM] #  rxn off vascular side membrane 
    + kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound_VM]
    - k_on_7_EXG2*C_EXG[1:N_Organs,bound_VM]*C_FcRn_VM[1:N_Organs] #  rxn on in vascular side membrane 
    + k_off_7_EXG2*C_EXG[1:N_Organs,bound2_VM] #  rxn off vascular side membrane 
    + 2*kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound2_VM] #  
    ))/V_VM[1:N_Organs]);
    
    #  endosomal pH=7.4
    @views @. du.C_FcRn_E7[1:N_Organs] = 
                ((CL_up[1:N_Organs]*FR*C_FcRn_VM[1:N_Organs] #  from vascular membrane
    + CL_up[1:N_Organs]*(1-FR)*C_FcRn_ISM[1:N_Organs] #  from is membrane
    - CL_up[1:N_Organs]*C_FcRn_E7[1:N_Organs] #  to endosomal pH=6
    + V_E7[1:N_Organs]*(
    - k_on_7_EDG*C_EDG[1:N_Organs,E7]*C_FcRn_E7[1:N_Organs] #  rxn on endosomal pH=7.4
    + k_off_7_EDG*C_EDG[1:N_Organs,bound_E7] #  rxn off endosomal pH=7.4.
    - k_on_7_EDG2*C_EDG[1:N_Organs,bound_E7]*C_FcRn_E7[1:N_Organs] #  rxn on endosomal pH=7.4. 
    + k_off_7_EDG2*C_EDG[1:N_Organs,bound2_E7] #  rxn off endosomal pH=7.4. 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    #  EXG
    - k_on_7_EXG*C_EXG[1:N_Organs,E7]*C_FcRn_E7[1:N_Organs] #  rxn on endosomal pH=7.4
    + k_off_7_EXG*C_EXG[1:N_Organs,bound_E7] #  rxn off endosomal pH=7.4
    - k_on_7_EXG2*C_EXG[1:N_Organs,bound_E7]*C_FcRn_E7[1:N_Organs] #  rxn on endosomal pH=7.4. 
    + k_off_7_EXG2*C_EXG[1:N_Organs,bound2_E7] #  rxn off endosomal pH=7.4. 
    ))/V_E7[1:N_Organs])# 
    
    #  endosomal pH=6.0
    @views @. du.C_FcRn_E6a[1:N_Organs] =  ((CL_up[1:N_Organs]*C_FcRn_E7[1:N_Organs] #  from E7
    - CL_up[1:N_Organs]*C_FcRn_E6a[1:N_Organs] #  from E6a to E7b
    + CL_up[1:N_Organs]*C_FcRn_E7b[1:N_Organs]*FcRn_recycle_fraction
    + V_E6a[1:N_Organs]*(
    - k_on_6_EDG*C_EDG[1:N_Organs,E6a]*C_FcRn_E6a[1:N_Organs] #  rxn on endosomal pH=6.0
    + k_off_6_EDG*C_EDG[1:N_Organs,bound_E6a] #  rxn off endosomal pH=6.0
    - k_on_6_EDG2*C_EDG[1:N_Organs,bound_E6a]*C_FcRn_E6a[1:N_Organs] #  rxn on endosomal pH=6.0. 
    + k_off_6_EDG2*C_EDG[1:N_Organs,bound2_E6a] #  rxn off endosomal pH=6.0. 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    #  EXG
    - k_on_6_EXG*C_EXG[1:N_Organs,E6a]*C_FcRn_E6a[1:N_Organs] #  rxn on endosomal pH=6.0
    + k_off_6_EXG*C_EXG[1:N_Organs,bound_E6a] #  rxn off endosomal pH=6.0
    - k_on_6_EXG2*C_EXG[1:N_Organs,bound_E6a]*C_FcRn_E6a[1:N_Organs] #  rxn on endosomal pH=6.0. 
    + k_off_6_EXG2*C_EXG[1:N_Organs,bound2_E6a] #  rxn off endosomal pH=6.0. 
    ))/V_E6a[1:N_Organs]) ; 

    #  endosomal pH=7.4
    @views @. du.C_FcRn_E7b[1:N_Organs] =  ((CL_up[1:N_Organs]*C_FcRn_E6a[1:N_Organs] #  from E6a
    - CL_up[1:N_Organs]*C_FcRn_E7b[1:N_Organs] #  from E7b to membranes
    + V_E7b[1:N_Organs]*(
    - k_on_7_EDG*C_EDG[1:N_Organs,E7b]*C_FcRn_E7b[1:N_Organs] #  rxn on endosomal pH=7.4
    + k_off_7_EDG*C_EDG[1:N_Organs,bound_E7b] #  rxn off endosomal pH=7.4
    - k_on_7_EDG2*C_EDG[1:N_Organs,bound_E7b]*C_FcRn_E7b[1:N_Organs] #  rxn on endosomal pH=7.4. 
    + k_off_7_EDG2*C_EDG[1:N_Organs,bound2_E7b] #  rxn off endosomal pH=7.4. 
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    #  EXG 
    - k_on_7_EXG*C_EXG[1:N_Organs,E7b]*C_FcRn_E7b[1:N_Organs] #  rxn on endosomal pH=7.4
    + k_off_7_EXG*C_EXG[1:N_Organs,bound_E7b] #  rxn off endosomal pH=7.4
    - k_on_7_EXG2*C_EXG[1:N_Organs,bound_E7b]*C_FcRn_E7b[1:N_Organs] #  rxn on endosomal pH=7.4. 
    + k_off_7_EXG2*C_EXG[1:N_Organs,bound2_E7b] #  rxn off endosomal pH=7.4. 
    ))/V_E7b[1:N_Organs]); 
    
    #  IS side mem
    @views @. du.C_FcRn_ISM = 
    ((- CL_up[1:N_Organs]*(1-FR)*C_FcRn_ISM[1:N_Organs] #  from IS membrane to endosomal pH=7.4
    + CL_up[1:N_Organs]*(1-FR)*C_FcRn_E7b[1:N_Organs]*(1-FcRn_recycle_fraction) #  from endosomal pH=7.4  b to IS side mem
    + V_ISM[1:N_Organs]*(
    - k_on_7_EDG*C_EDG[1:N_Organs,ISM]*C_FcRn_ISM[1:N_Organs] #  rxn on IS side mem
    + k_off_7_EDG*C_EDG[1:N_Organs,bound_ISM] #  rxn off IS side mem
    + kdeg_FcRn_Ab*C_EDG[1:N_Organs,bound_ISM]
    - k_on_7_EDG2*C_EDG[1:N_Organs,bound_ISM]*C_FcRn_ISM[1:N_Organs] #  rxn on IS side mem. 
    + k_off_7_EDG2*C_EDG[1:N_Organs,bound2_ISM] #  rxn off IS side mem. 
    + 2*kdeg_FcRn_Ab*C_EDG[1:N_Organs,bound2_ISM] #  
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    #  EXG 
    - k_on_7_EXG*C_EXG[1:N_Organs,ISM]*C_FcRn_ISM[1:N_Organs] #  rxn on IS side mem
    + k_off_7_EXG*C_EXG[1:N_Organs,bound_ISM] #  rxn off IS side mem
    + kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound_ISM]
    - k_on_7_EXG2*C_EXG[1:N_Organs,bound_ISM]*C_FcRn_ISM[1:N_Organs] #  rxn on IS side mem. 
    + k_off_7_EXG2*C_EXG[1:N_Organs,bound2_ISM] #  rxn off IS side mem. 
    + 2*kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound2_ISM] #  
    ))/V_ISM[1:N_Organs]);


    # ; flux of endogenous IgG in Lymph node
    du.C_EDG_LN = 
        (((1-sigma_IS[Heart])*LF[Heart]*C_EDG[Heart, IntS]
    + (1-sigma_IS[Kidney])*LF[Kidney]*C_EDG[Kidney, IntS]
    + (1-sigma_IS[Muscle])*LF[Muscle]*C_EDG[Muscle, IntS]
    + (1-sigma_IS[Skin])*LF[Skin]*C_EDG[Skin, IntS]
    + (1-sigma_IS[Brain])*LF[Brain]*C_EDG[Brain, IntS]
    + (1-sigma_IS[Adipose])*LF[Adipose]*C_EDG[Adipose, IntS]
    + (1-sigma_IS[Thymus])*LF[Thymus]*C_EDG[Thymus, IntS]
    + (1-sigma_IS[Liver])*LF[Liver]*C_EDG[Liver, IntS]
    + (1-sigma_IS[Spleen])*LF[Spleen]*C_EDG[Spleen, IntS]
    + (1-sigma_IS[Pancreas])*LF[Pancreas]*C_EDG[Pancreas, IntS]
    + (1-sigma_IS[SI])*LF[SI]*C_EDG[SI, IntS]
    + (1-sigma_IS[LI])*LF[LI]*C_EDG[LI, IntS]
    + (1-sigma_IS[Bone])*LF[Bone]*C_EDG[Bone, IntS]
    + (1-sigma_IS[Other])*LF[Other]*C_EDG[Other, IntS]
    # + (1-sigma_IS[Tumor])*LF[Tumor]*C_EDG[Tumor, IntS]
    + (1-sigma_IS[Lung])*LF[Lung]*C_EDG[Lung, IntS]
    - L_LymphNode*C_EDG_LN)/V_LN);         
        
    # ; flux of exogenous IgG in Plasma
    w_other_1_EXG_Plasma = 
    ((PLQ[Heart]-LF[Heart])*C_EXG[Heart,V] 
    +(PLQ[Kidney]-LF[Kidney])*C_EXG[Kidney,V]
    +(PLQ[Muscle]-LF[Muscle])*C_EXG[Muscle,V]
    +(PLQ[Skin]-LF[Skin])*C_EXG[Skin,V]
    +(PLQ[Brain]-LF[Brain])*C_EXG[Brain,V]
    +(PLQ[Adipose]-LF[Adipose])*C_EXG[Adipose,V]
    +(PLQ[Thymus]-LF[Thymus])*C_EXG[Thymus,V]
    +(PLQ[Liver]-LF[Liver])*C_EXG[Liver,V]
        +(PLQ[Spleen]-LF[Spleen])*C_EXG[Liver,V]
    +(PLQ[Pancreas]-LF[Pancreas])*C_EXG[Liver,V]
    +(PLQ[SI]-LF[SI])*C_EXG[Liver,V]
    +(PLQ[LI]-LF[LI])*C_EXG[Liver,V]
    +(PLQ[Bone]-LF[Bone])*C_EXG[Bone,V]
    +(PLQ[Other]-LF[Other])*C_EXG[Other,V]
    # +(PLQ[Tumor]-LF[Tumor])*C_EXG[Tumor, V]
    -(PLQ[Lung]+LF[Lung])*C_EXG_Plasma
    +L_LymphNode*C_EXG_LN)
    
    

    # ; flux of exogenous IgG in Lymph node
    du.C_EXG_LN = 
    (((1-sigma_IS[Heart])*LF[Heart]*C_EXG[Heart, IntS]
    + (1-sigma_IS[Kidney])*LF[Kidney]*C_EXG[Kidney, IntS]
    + (1-sigma_IS[Muscle])*LF[Muscle]*C_EXG[Muscle, IntS]
    + (1-sigma_IS[Skin])*LF[Skin]*C_EXG[Skin, IntS]
    + (1-sigma_IS[Brain])*LF[Brain]*C_EXG[Brain, IntS]
    + (1-sigma_IS[Adipose])*LF[Adipose]*C_EXG[Adipose, IntS]
    + (1-sigma_IS[Thymus])*LF[Thymus]*C_EXG[Thymus, IntS]
    + (1-sigma_IS[Liver])*LF[Liver]*C_EXG[Liver, IntS]
    + (1-sigma_IS[Spleen])*LF[Spleen]*C_EXG[Spleen, IntS]
    + (1-sigma_IS[Pancreas])*LF[Pancreas]*C_EXG[Pancreas, IntS]
    + (1-sigma_IS[SI])*LF[SI]*C_EXG[SI, IntS]
    + (1-sigma_IS[LI])*LF[LI]*C_EXG[LI, IntS]
    + (1-sigma_IS[Bone])*LF[Bone]*C_EXG[Bone, IntS]
    + (1-sigma_IS[Other])*LF[Other]*C_EXG[Other, IntS]
    # + (1-sigma_IS[Tumor])*LF[Tumor]*C_EXG[Tumor, IntS]
    + (1-sigma_IS[Lung])*LF[Lung]*C_EXG[Lung, IntS]
    - L_LymphNode*C_EXG_LN
            )/V_LN); 

    du.C_EDG_Plasma = 0.0;

    # ; flux related to soluble HER2 dynamics in plasma 
    kdeg_sR = log(2)/ thalf_sR;
    kdeg_sR_adc = log(2)/ thalf_sR_adc;
    ksyn_sR = kdeg_sR * init_sR; 
    Koff_sR_EXG = Kd*Kon; 

    du.C_sR_plasma = ksyn_sR - kdeg_sR*C_sR_plasma - Kon*C_sR_plasma*C_EXG_Plasma + Koff_sR_EXG*C_sR_EXG;
    du.C_sR_EXG = Kon*C_sR_plasma*C_EXG_Plasma - Koff_sR_EXG*C_sR_EXG - kdeg_sR_adc*C_sR_EXG;


    # tumor dynamics
    V_IntS_Tumor = 0.55 * TumorVolume; # Lindauer et al., 2017; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5270293/

    Kgrow = log(2)/tdouble   
    Koff = Kd*Kon
   
    Ksyn = Rcopies/N_av*k_deg*1e6 # unit in umol/hr

    
    Ntot = Nc_1 + Nc_2 + Nc_3 + Nc_4
    Ntot_neg = Nc_1_neg + Nc_2_neg + Nc_3_neg + Nc_4_neg
    C_P_c = P_c/Vc                  # payload conc inside target+ cells
    C_P_neg_c = P_neg_c/ Vc         # payload conc inside target- cells
    C_P_m = P_m/V_IntS_Tumor

    flux_A_R_s_binding = A_m*R_s*Kon/V_IntS_Tumor
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
    Kkill_eff = Emax_Payload*C_P_c/(IC50_Payload + C_P_c)*Kkill_scale

    Du_A_m = flux_AR_s_unbinding * Ntot - flux_A_R_s_binding * Ntot

    du.tumor.R_s = flux_R_s_syn - flux_R_s_degrade + flux_R_e_rec - flux_R_s_int + flux_AR_s_unbinding - flux_A_R_s_binding 
    du.tumor.R_e = flux_R_s_int - flux_R_e_rec 

    du.tumor.AR_s = flux_A_R_s_binding - flux_AR_s_unbinding - flux_AR_s_int + flux_AR_e_rec
    du.tumor.AR_e = flux_AR_s_int - flux_AR_e_rec - flux_AR_cat 

    du.tumor.Nc_1 = Kgrow*Nc_1 - Kkill_eff*Nc_1
    du.tumor.Nc_2 = Kkill_eff*Nc_1 - Nc_2/tau
    du.tumor.Nc_3 = (Nc_2 - Nc_3)/tau
    du.tumor.Nc_4 = (Nc_3 - Nc_4)/tau
    
    du.tumor.P_c = flux_P_l_to_P_c - flux_P_c_to_P_m
    
    # warhead in tumor interstitial space
    Du_P_m = flux_P_c_to_P_m * Ntot + DAR*1/tau*Nc_4*(AR_e + AR_s) - flux_P_neg_c_from_P_m * Ntot_neg
    du.tumor.P_m = Du_P_m - k_PL * P_m; 

    # receptor-negative cells
    Kgrow_neg = log(2)/ tdouble     
    Kkill_eff_neg = Emax_Payload*C_P_neg_c/(IC50_Payload + C_P_neg_c)*Kkill_scale

    du.tumor.Nc_1_neg = Kgrow_neg*Nc_1_neg - Kkill_eff_neg*Nc_1_neg
    du.tumor.Nc_2_neg = Kkill_eff_neg*Nc_1_neg - Nc_2_neg/tau
    du.tumor.Nc_3_neg = (Nc_2_neg - Nc_3_neg)/tau
    du.tumor.Nc_4_neg = (Nc_3_neg - Nc_4_neg)/tau
    du.tumor.P_neg_c = flux_P_neg_c_from_P_m


    # vascular & surface-based exchange of  tumor
    Rtumor = ( TumorVolume * 1e3/(4/3*pi) ) ^ (1/3) # [cm]
    surface_exchange = 6 * D_ADC/ (Rtumor^2) ; # [h-1]
    vascular_exchange = 2 * P_ADC * Rcap / (Rkrogh^2) ; # [h-1]

    flux_ADC_plasma2tumor = (vascular_exchange + surface_exchange) * (C_EXG_Plasma * epsilon - A_m/V_IntS_Tumor) * V_IntS_Tumor

    du.tumor.A_m = flux_ADC_plasma2tumor + Du_A_m
    # du.C_EXG[Tumor,IntS] = ( du.C_EXG[Tumor,IntS] + Du_A_m/V_IntS[Tumor] ) # ADC in media (µM),  Ab = C_EXG[Tumor,IntS]   

    # update ADC conc in plasma
    du.C_EXG_Plasma = ((w_other_1_EXG_Plasma - flux_ADC_plasma2tumor)/V_Plasma - Kon*C_sR_plasma*C_EXG_Plasma + Koff_sR_EXG*C_sR_EXG) + infusion; 
    

    # degraded mAb/ IgG
    @views @. du.DEGprotein.deg_EXG[1:N_Organs] = CL_up[1:N_Organs]*Prob_deg*C_EXG[1:N_Organs,E6a] + kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound_VM]*V_VM[1:N_Organs] + kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound2_VM]*V_VM[1:N_Organs] + kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound_ISM]*V_ISM[1:N_Organs] + kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound2_ISM]*V_ISM[1:N_Organs]
    @views @. du.DEGprotein.deg_EDG[1:N_Organs] = CL_up[1:N_Organs]*Prob_deg*C_EDG[1:N_Organs,E6a] + kdeg_FcRn_Ab*C_EDG[1:N_Organs,bound_VM]*V_VM[1:N_Organs] + kdeg_FcRn_Ab*C_EDG[1:N_Organs,bound2_VM]*V_VM[1:N_Organs] + kdeg_FcRn_Ab*C_EDG[1:N_Organs,bound_ISM]*V_ISM[1:N_Organs] + kdeg_FcRn_Ab*C_EDG[1:N_Organs,bound2_ISM]*V_ISM[1:N_Organs]
    du.DEGprotein.deg_TumorCellular = flux_AR_cat * Ntot
    du.DEGprotein.deg_plasma = kdeg_sR_adc*C_sR_EXG

    # payload concentration in tissue endothetial cells' endosomes
    @views @. du.end_endo_payload = CL_up[1:N_Organs]*Prob_deg*C_EXG[1:N_Organs,E6a]/V_endosomal[1:N_Organs] * DAR + 
       (kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound_VM]*V_VM[1:N_Organs] + 
        kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound2_VM]*V_VM[1:N_Organs] + 
        kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound_ISM]*V_ISM[1:N_Organs] + 
        kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound2_ISM]*V_ISM[1:N_Organs]
    )/V_endosomal[1:N_Organs]*DAR - k_in*(end_endo_payload - end_cyto_payload);

    # payload concentration in tissue endothetial cells' cytosol
    V_cytosol = V_Organ - V_V - V_IntS - V_endosomal
    @views @. du.end_cyto_payload = (
        k_in*(end_endo_payload - end_cyto_payload)*V_endosomal/V_cytosol[1:N_Organs]   # payload diffuse from endosome into cytosol 
        - k_out*(end_cyto_payload-ints_payload)  # payload from cytosol to tissue interstitium 
        - k_out*end_cyto_payload  # payload diffuse from endothelial cell cytosol to plasma 
        )

    @views @. du.ints_payload = k_out*(end_cyto_payload-ints_payload)*V_cytosol[1:N_Organs]/V_IntS[1:N_Organs] - k_PL*ints_payload
    
end