# author: Yuezhe Li 
# date: 9/26/2023
# this version of the Jones model for cyno monkey 

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
const V_C = zeros(N_Organs)
const V_Organ = zeros(N_Organs)
const PLQ = zeros(N_Organs)
const LF = zeros(N_Organs)


function jonesODEs3Monkey!(du,u,p,t)

    @unpack C_EXG_Plasma, C_EDG_Plasma, C_EDG_LN, C_EDG, C_EXG_LN, C_EXG, C_FcRn_E6a, C_FcRn_E7, C_FcRn_E7b, C_FcRn_ISM, C_FcRn_VM, DEGprotein = u

    @unpack PS_Score, PS_kd, KD6_WT, endothelial_scaling_factor, infusion = p

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
     # sigma_V = zeros(N_Organs)
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
 
     sigma_IS = 0.2 * ones(15)
 
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
 
     # Vascular Volumes
     V_V[Lung] = 1.96/1000
     V_V[Liver] = 15.9/1000
     V_V[Heart] = 1.09/1000
     V_V[Muscle] = 72/1000
     V_V[Skin] = 25.2/1000
     V_V[Adipose] = 1.7/1000
     V_V[Bone] = 20.9/1000
     V_V[Brain] = 2.07/1000
     V_V[Kidney] = 1.5/1000
     V_V[SI] = 1.57/1000
     V_V[LI] = 2.24/1000
     V_V[Pancreas] = 0.685/1000
     V_V[Thymus] = 0.115/1000
     V_V[Spleen] = 0.72/1000
     V_V[Other] = 5.56/1000
 
     # Interstitial Volumes
     V_IntS[Lung] = 10.7/1000
     V_IntS[Liver] = 37.4/1000
     V_IntS[Heart] = 4.05/1000
     V_IntS[Muscle] = 426/1000
     V_IntS[Skin] = 223/1000
     V_IntS[Adipose] = 26.3/1000
     V_IntS[Bone] = 177/1000
     V_IntS[Brain] = 16.9/1000
     V_IntS[Kidney] = 4.09/1000
     V_IntS[SI] = 17.2/1000
     V_IntS[LI] = 24.4/1000
     V_IntS[Pancreas] = 2.17/1000
     V_IntS[Thymus] = 0.355/1000
     V_IntS[Spleen] = 1.19/1000
     V_IntS[Other] = 22.7/1000
 
     # Cellular Volumes
     V_C[Lung] = 21.2/1000
     V_C[Liver] = 120/1000
     V_C[Heart] = 22.2/1000
     V_C[Muscle] = 2701/1000
     V_C[Skin] = 403/1000
     V_C[Adipose] = 124/1000
     V_C[Bone] = 732/1000
     V_C[Brain] = 72.9/1000
     V_C[Kidney] = 20.3/1000
     V_C[SI] = 78.2/1000
     V_C[LI] = 111/1000
     V_C[Pancreas] = 8.98/1000
     V_C[Thymus] = 1.52/1000
     V_C[Spleen] = 3.42/1000
     V_C[Other] = 98.9/1000
     
     # Organ Volumes
     V_Organ = zeros(N_Organs)
     V_Organ[Lung] = 35.7/1000 
     V_Organ[Liver] = 187/1000
     V_Organ[Heart] = 28.3/1000
     V_Organ[Muscle] = 3273/1000
     V_Organ[Skin] = 674/1000
     V_Organ[Adipose] = 154/1000
     V_Organ[Bone] = 952/1000
     V_Organ[Brain] = 94/1000
     V_Organ[Kidney] = 27.3/1000
     V_Organ[SI] = 98.7/1000
     V_Organ[LI] = 140/1000
     V_Organ[Pancreas] = 12.5/1000
     V_Organ[Thymus] = 2.09/1000
     V_Organ[Spleen] = 5.95/1000
     V_Organ[Other] = 132/1000
 
     V_Plasma = 0.187 
     V_LN     = 25.1/1000
 
     # Blood Flows
     PLQ[Lung] = 22433/1000 
     PLQ[Liver] = 1251/1000
     PLQ[Heart] = 696/1000
     PLQ[Muscle] = 3944/1000
     PLQ[Skin] = 2492/1000
     PLQ[Adipose] = 139/1000
     PLQ[Bone] = 248/1000
     PLQ[Brain] = 1508/1000
     PLQ[Kidney] = 3237/1000
     PLQ[SI] = 3432/1000
     PLQ[LI] = 3571/1000
     PLQ[Pancreas] = 399/1000
     PLQ[Thymus] = 125/1000
     PLQ[Spleen] = 185/1000
     PLQ[Other] = 1206/1000   # this is adjusted based on PLQ balance, not from Shah and Betts, 2011; https://pubmed.ncbi.nlm.nih.gov/22143261/
 
     # Lymph Flows    
     LF[Lung]     = PLQ[Lung]*0.002
     LF[Liver]    = (PLQ[Liver] + PLQ[SI]-LF[SI] + PLQ[LI]-LF[LI] + PLQ[Spleen]-LF[Spleen] + PLQ[Pancreas] - LF[Pancreas])*0.002 
     LF[Heart]    = PLQ[Heart]*0.002
     LF[Muscle]   = PLQ[Muscle]*0.002
     LF[Skin]     = PLQ[Skin]*0.002
     LF[Adipose]  = PLQ[Adipose]*0.002
     LF[Bone]     = PLQ[Bone]*0.002
     LF[Brain]    = PLQ[Brain]*0.002
     LF[Kidney]   = PLQ[Kidney]*0.002
     LF[SI]       = PLQ[SI]*0.002
     LF[LI]       = PLQ[LI]*0.002
     LF[Pancreas] = PLQ[Pancreas]*0.002
     LF[Thymus]   = PLQ[Thymus]*0.002
     LF[Spleen]   = PLQ[Spleen]*0.002
     LF[Other]    = PLQ[Other]*0.002
     L_LymphNode = sum(LF) # [L/h]
 
    # Setup
    PS_a = 1.8051
    PS_b = 0.2624
    if PS_Score > -1.0
        PS_Kd = 10^(exp(PS_a - PS_b*PS_Score)) 
    else
        PS_Kd = PS_kd
    end
     k_on_PS = 7.92E+08/1E6 # unit uM.h-1 (rhesus monkey)
     pino_time = 10.8/60  # min->h
     # Scale_Factor = 603.7  # scales number of endothelial cells to value for human 
     Scale_Factor = endothelial_scaling_factor # (fitted for monkey)
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
    @. V_E7 = V_endosomal * E7_Vol_Pct
    @. V_E6a = V_endosomal * E6a_Vol_Pct
    @. V_E7b = V_endosomal * E7b_Vol_Pct

    # kon for EXG    
    k_on_6_EXG  = 8.06E+07/1E6  # [1/(uM*h)]
    k_on_7_EXG  = 1.61E+07/5/1E6  # [1/(uM*h)]
    
    # kon for EDG    
    k_on_6_EDG  = 8.06E+07/1E6 # [1/(uM*h)]
    k_on_7_EDG  = 1.61E+07/5/1E6 # [1/(uM*h)]

    KD7_WT = KD6_WT*220 # [nM]  

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
+ (1-sigma_IS[Lung])*LF[Lung]*C_EXG[Lung, IntS]
- L_LymphNode*C_EXG_LN
        )/V_LN); 

du.C_EDG_Plasma = 0.0;
du.C_EXG_Plasma = (w_other_1_EXG_Plasma/V_Plasma) + infusion; 

# degraded mAb/ IgG
@views @. du.DEGprotein.deg_EXG[1:N_Organs] = CL_up[1:N_Organs]*Prob_deg*C_EXG[1:N_Organs,E6a] + kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound_VM]*V_VM[1:N_Organs] + kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound2_VM]*V_VM[1:N_Organs] + kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound_ISM]*V_ISM[1:N_Organs] + kdeg_FcRn_Ab*C_EXG[1:N_Organs,bound2_ISM]*V_ISM[1:N_Organs]
@views @. du.DEGprotein.deg_EDG[1:N_Organs] = CL_up[1:N_Organs]*Prob_deg*C_EDG[1:N_Organs,E6a] + kdeg_FcRn_Ab*C_EDG[1:N_Organs,bound_VM]*V_VM[1:N_Organs] + kdeg_FcRn_Ab*C_EDG[1:N_Organs,bound2_VM]*V_VM[1:N_Organs] + kdeg_FcRn_Ab*C_EDG[1:N_Organs,bound_ISM]*V_ISM[1:N_Organs] + kdeg_FcRn_Ab*C_EDG[1:N_Organs,bound2_ISM]*V_ISM[1:N_Organs]          

end

