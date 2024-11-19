# author: Yuezhe LI
# date: 3/14/24
# purpose of the code: to add blood cells to simplify Jones model 
# prior Jones model: see https://ghe.metrumrg.com/yuezhel/ADC-toxicity-2023/blob/main/jones_tumor_homo.jl
# BC compartment: obtained from Shah and Betts, 2012; https://pubmed.ncbi.nlm.nih.gov/22143261/

const N_av = 6.0221409e+23 # Avogadro constant 
using ComponentArrays
using Parameters: @unpack


function mPBPK_bc_homo_tumor!(du,u,p,t) 
    @unpack C_EXG_Plasma, C_EDG_Plasma, C_EDG_LN, C_EDG, C_EXG_LN, C_EXG, C_FcRn_E6a, C_FcRn_E7, C_FcRn_E7b, C_FcRn_ISM, C_FcRn_VM, 
            C_sR_plasma, C_sR_EXG, 
            C_EXG_Plasma_BC, C_BC_EXG, 
            tumor, DEGprotein, end_endo_payload, end_cyto_payload, ints_payload = u
    
    @unpack Nc_1, Nc_2, Nc_3, Nc_4, Nc_1_neg, Nc_2_neg, Nc_3_neg, Nc_4_neg, R_s, R_e, AR_s, AR_e, P_c, P_m, P_neg_c, A_m  = tumor
    
    @unpack infusion, tumor_cell_radius, tdouble, tau, PS_Score, PS_kd, KD6_WT, DAR, 
            epsilon, P_ADC, D_ADC, Rcap, Rkrogh, thalf_sR, thalf_sR_adc, init_sR, Kkill_scale, Emax_Payload, IC50_Payload, k_PL, 
            Rcopies, k_deg, k_rec, k_endo, Kon, Kd, k_out, k_in, k_lys = p
            
    N_Organs = 7
    V_endosomal = zeros(N_Organs)
    V_BC = zeros(N_Organs)
    CL_up = zeros(N_Organs)
    V_VM = zeros(N_Organs)
    V_ISM = zeros(N_Organs)
    V_E7 = zeros(N_Organs)
    V_E6a = zeros(N_Organs)
    V_E7b = zeros(N_Organs)
    Organ_Endothelial_Cell = zeros(N_Organs)
    sigma_V = zeros(N_Organs)
    sigma_IS = ones(N_Organs) * 0.2
    Endothelial_Cell_Frac = zeros(N_Organs)
    V_V = zeros(N_Organs)
    V_IntS = zeros(N_Organs)
    V_Organ = zeros(N_Organs)
    PLQ = zeros(N_Organs)
    BCQ = zeros(N_Organs)
    LF = zeros(N_Organs)

    # calculate tumor volume 
    TumorVolume = (Nc_1 + Nc_2 + Nc_3 + Nc_4 + Nc_1_neg + Nc_2_neg + Nc_3_neg + Nc_4_neg) * Vc / 0.375;
        
    # Organ Indices
    Lung = 1
    Liver = 2
    Skin = 3
    Kidney = 4
    Intestine = 5 # this should be a combination of LI, SI, and pacreas
    Spleen = 6
    Other = 7

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
    sigma_V[Skin] = 0.95
    sigma_V[Kidney] = 0.9
    sigma_V[Intestine] = 0.9 # this should be either 0.9 or 0.95
    sigma_V[Spleen] = 0.85
    sigma_V[Other] = 0.95 

    # Endothelial Cell Fractions
    Endothelial_Cell_Frac[Lung] = 0.0834
    Endothelial_Cell_Frac[Liver] = 0.1877
    Endothelial_Cell_Frac[Skin] = 0.0819
    Endothelial_Cell_Frac[Kidney] = 0.0157
    Endothelial_Cell_Frac[Intestine] = 0.0331 # combined LI, SI, Pancreas
    Endothelial_Cell_Frac[Spleen] = 0.0499
    Endothelial_Cell_Frac[Other] = 0.5483    # combined Heart, Muscle, Adipose, Bone, Brain, Thymus, and Other

    # Vascular Volumes [L]
    # based on male BW = 71kg
    V_V[Lung] = 55/1000
    V_V[Liver] = 183/1000
    V_V[Skin] = 127/1000
    V_V[Kidney] = 18.2/1000
    V_V[Intestine] = (6.15+8.74+5.7)/1000 # combined SI, LI, Pancreas
    V_V[Spleen] = 26.8/1000
    V_V[Other] = (13.1+662+148+224+31.9+0.353+204)/1000     # combine Heart, Muscle, Adipose, Bone, Brain, Thymus, and Other

    # blood cell volume [L]
    # based on male BW = 71kg
    V_BC[Lung] = 45/1000
    V_BC[Liver] = 149/1000
    V_BC[Skin] = 104/1000
    V_BC[Kidney] = 14.9/1000
    V_BC[Intestine] = (5.03+7.15+4.66)/1000 # combined SI, LI, Pancreas
    V_BC[Spleen] = 21.9/1000
    V_BC[Other] = (10.8+541+121+183+26.1+0.288+167)/1000     # combine Heart, Muscle, Adipose, Bone, Brain, Thymus, and Other
    V_BC_plasma = 2558/1000

    # Interstitial Volumes [L]
    # based on male BW = 71kg
    V_IntS[Lung] = 300/1000
    V_IntS[Liver] = 429/1000    
    V_IntS[Skin] = 1125/1000
    V_IntS[Kidney] = 49.8/1000
    V_IntS[Intestine] = (67.1+95.3+18)/1000   # combined SI, LI, Pancreas
    V_IntS[Spleen] = 44.3/1000
    V_IntS[Other] = (48.8+3910+2289+1891+261+1.09+831)/1000      # combine Heart, Muscle, Adipose, Bone, Brain, Thymus, and Other

    # Organ Volumes [L]
    # based on male BW = 71kg
    V_Organ[Lung] = 1000/1000 
    V_Organ[Liver] = 2143/1000
    V_Organ[Skin] = 3408/1000
    V_Organ[Kidney] = 332/1000
    V_Organ[Intestine] = (385+548+104)/1000   # combined SI, LI, Pancreas
    V_Organ[Spleen] = 221/1000
    V_Organ[Other] = (341+30078+13465+10165+1450+6.41+4852)/1000    # combine Heart, Muscle, Adipose, Bone, Brain, Thymus, and Other

    # scaled linearly based on the BW using ref. BW = 71 kg
    V_Plasma = 3.126; # [L]
    V_LN     = 0.274; # [L]

    # Blood Flows [L/h]
    # scaled with an allometric exponent of 0.75 based on the BW using ref. BW = 71 kg
    # PLQ[Lung] = 181913/1000 
    PLQ[Liver] = 13210/1000
    PLQ[Skin] = 11626/1000
    PLQ[Kidney] = 36402/1000
    PLQ[Intestine] = (12368+12867+3056)/1000      # combined SI, LI, Pancreas
    PLQ[Spleen] = 6343/1000
    PLQ[Other] = (7752+33469+11233+2591+21453+353+9190)/1000       # combine Heart, Muscle, Adipose, Bone, Brain, Thymus, and Other
    # rebalance PLQ
    PLQ[Lung] = sum(PLQ[2:end])

    # blood cell flow [L/h]
    # obtained from Shah and Betts, 2012; https://pubmed.ncbi.nlm.nih.gov/22143261/
    # based on male BW = 71kg
    BCQ[Lung] = 148838/1000
    BCQ[Liver] = 10808/1000
    BCQ[Skin] = 9512/1000
    BCQ[Kidney] = 29784/1000
    BCQ[Intestine] = (10120 + 10527 + 2500)/1000  # combined SI, LI, Pancreas
    BCQ[Spleen] = 5189/1000
    BCQ[Other] = (6342+27383+9191+2120+17553+289+4517)/1000    # combine Heart, Muscle, Adipose, Bone, Brain, Thymus, and Other
    BCQ_plasma = 148838/1000

    # Lymph Flows [L/h]
    LF_PLQ_frac = 0.002
    LF[Lung]     = PLQ[Lung]*LF_PLQ_frac
    LF[Liver]    = (PLQ[Liver] + PLQ[Intestine]-LF[Intestine] + PLQ[Spleen]-LF[Spleen])*LF_PLQ_frac 
    LF[Skin]     = PLQ[Skin]*LF_PLQ_frac
    LF[Kidney]   = PLQ[Kidney]*LF_PLQ_frac
    LF[Intestine]= PLQ[Intestine]*LF_PLQ_frac
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

    FcRn_recycle_fraction = 0.99

    FR    = 0.715

        # ===================================================================================== 
    #  ---------------------- DIFFERENTIAL VARIABLES/EQUATIONS ----------------------------
    # =====================================================================================
        
    # blood cell for plasma 
    du.C_EXG_Plasma_BC = ( BCQ[Kidney]*C_BC_EXG[Kidney] + BCQ[Skin]*C_BC_EXG[Skin] + BCQ[Other]*C_BC_EXG[Other] + 
                            (BCQ[Liver]+BCQ[Spleen]+BCQ[Intestine])*C_BC_EXG[Liver] - BCQ[Lung]*C_BC_EXG[Lung] )/ BCQ_plasma; 

    # blood cell for lung 
    # du.C_BC_EXG[Lung] = (BCQ[Lung]*C_BC_EXG[Lung] - sum(BCQ[2:N_Organs])*C_BC_EXG[Lung])/V_BC[Lung]
    du.C_BC_EXG[Lung] = BCQ[Lung]*(C_EXG_Plasma_BC - C_BC_EXG[Lung])/V_BC[Lung]

    # added blood cell (all other organs)
    @views @. du.C_BC_EXG[3:N_Organs] = (BCQ[3:N_Organs] * (C_BC_EXG[Lung] - C_BC_EXG[3:N_Organs]) )/V_BC[3:N_Organs]; 

    # Liver, BC 
    du.C_BC_EXG[Liver] = (BCQ[Liver]*C_BC_EXG[Lung] + BCQ[Spleen]*C_BC_EXG[Spleen] + BCQ[Intestine]*C_BC_EXG[Intestine] 
                         - (BCQ[Liver]+BCQ[Spleen]+BCQ[Intestine])*C_BC_EXG[Liver] )/V_BC[Liver];

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
        + (PLQ[Intestine]-LF[Intestine])*C_EXG[Intestine,V] #  Inlet from LI, SI, Pancreas
        - C_EXG[Liver,V]*(PLQ[Liver]-LF[Liver] + PLQ[Spleen]-LF[Spleen] + PLQ[Intestine]-LF[Intestine]) #  Outlet
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
        + (PLQ[Intestine]-LF[Intestine])*C_EDG[Intestine,V] #  Inlet from LI, SI, Pancreas
        - C_EDG[Liver,V]*(PLQ[Liver]-LF[Liver] + PLQ[Spleen]-LF[Spleen] + PLQ[Intestine]-LF[Intestine]) #  Outlet
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
    )/V_IntS[1:N_Organs]); 

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
    du.C_EDG_LN = (
    ( (1-sigma_IS[Kidney])*LF[Kidney]*C_EDG[Kidney, IntS]
    + (1-sigma_IS[Skin])*LF[Skin]*C_EDG[Skin, IntS]
    + (1-sigma_IS[Liver])*LF[Liver]*C_EDG[Liver, IntS]
    + (1-sigma_IS[Spleen])*LF[Spleen]*C_EDG[Spleen, IntS]
    + (1-sigma_IS[Intestine])*LF[Intestine]*C_EDG[Intestine, IntS]
    + (1-sigma_IS[Other])*LF[Other]*C_EDG[Other, IntS]
    + (1-sigma_IS[Lung])*LF[Lung]*C_EDG[Lung, IntS]
    - L_LymphNode*C_EDG_LN)/V_LN);         
        
    # ; flux of exogenous IgG in Plasma
    w_other_1_EXG_Plasma = 
    (
     (PLQ[Kidney]-LF[Kidney])*C_EXG[Kidney,V]
    +(PLQ[Skin]-LF[Skin])*C_EXG[Skin,V]
    +(PLQ[Liver]-LF[Liver])*C_EXG[Liver,V]
    +(PLQ[Spleen]-LF[Spleen])*C_EXG[Liver,V]
    +(PLQ[Intestine]-LF[Intestine])*C_EXG[Liver,V]
    +(PLQ[Other]-LF[Other])*C_EXG[Other,V]
    -(PLQ[Lung]+LF[Lung])*C_EXG_Plasma
    +L_LymphNode*C_EXG_LN)

    # ; flux of exogenous IgG in Lymph node
    du.C_EXG_LN = 
    ((
      (1-sigma_IS[Kidney])*LF[Kidney]*C_EXG[Kidney, IntS]
    + (1-sigma_IS[Skin])*LF[Skin]*C_EXG[Skin, IntS]
    + (1-sigma_IS[Liver])*LF[Liver]*C_EXG[Liver, IntS]
    + (1-sigma_IS[Spleen])*LF[Spleen]*C_EXG[Spleen, IntS]
    + (1-sigma_IS[Intestine])*LF[Intestine]*C_EXG[Intestine, IntS]
    + (1-sigma_IS[Other])*LF[Other]*C_EXG[Other, IntS]
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


    # Tumor Dynamics
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
    @views @. du.DEGprotein.deg_EXG[1:N_Organs] = CL_up[1:N_Organs]*Prob_deg*C_EXG[1:N_Organs,E6a] + 
                                                  kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound_VM]*V_VM[1:N_Organs] + 
                                                  kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound2_VM]*V_VM[1:N_Organs] + 
                                                  kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound_ISM]*V_ISM[1:N_Organs] + 
                                                  kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound2_ISM]*V_ISM[1:N_Organs]
    @views @. du.DEGprotein.deg_EDG[1:N_Organs] = CL_up[1:N_Organs]*Prob_deg*C_EDG[1:N_Organs,E6a] + 
                                                  kdeg_FcRn_Ab*C_EDG[1:N_Organs,bound_VM]*V_VM[1:N_Organs] + 
                                                  kdeg_FcRn_Ab*C_EDG[1:N_Organs,bound2_VM]*V_VM[1:N_Organs] + 
                                                  kdeg_FcRn_Ab*C_EDG[1:N_Organs,bound_ISM]*V_ISM[1:N_Organs] + 
                                                  kdeg_FcRn_Ab*C_EDG[1:N_Organs,bound2_ISM]*V_ISM[1:N_Organs]
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