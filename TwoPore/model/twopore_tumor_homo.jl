# author: Yuezhe Li 
# date: Dec 1, 2023
# purpose of this script: to generate a dummy human model with tumor compartment
# tumor penentration is adopted from Shah, Haddish-Berhane, and Betts, 2012; 
# https://pubmed.ncbi.nlm.nih.gov/23151991/

using Parameters: @unpack 
using ComponentArrays


function lishah_tumor_homo!(du, u, p, t)
    @unpack C_EXG_Plasma, C_EXG_LN, C_EXG, endo_deg, renal_clearance, C_EXG_Tumor, Nc1, Nc2, Nc3, Nc4 = u

    @unpack MW, k_deg, infusion, Vc, P_protein, D_protein, Rcap, Rkrogh, epsilon = p
    
    a_e = 0.0483 * MW^0.386; # a_e in nm, MW in Dalton
    theta = exp(1 - 8.7/(1 + exp(0.028*(-MW/1E3+72.3)))); # MW in Da
    ratio_As_A0s = 0.2352*exp(-0.00008295*MW) + 0.7767*exp(-0.00053095*MW); # MW in Dalton
    ratio_Al_A0l = 0.3429*exp(-0.00012175*MW) + 0.6571*exp(-0.00000421*MW); # MW in Dalton
    sigma_L = 0.000035*MW^0.717; # MW in Dalton
    sigma_S = 1-0.8489*exp(-0.00004*MW); # MW in Dalton
    
    # Organ Indices
    N_Organs = 15
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
    E = 2
    IS = 3

    # Organ Volumes, [L]; homo, 71kg
    global V_Plasma = 3.126 
    global V_LN     = 0.274

    # Organ Volumes, [L]; homo, 71kg
    global V_Organ = zeros(N_Organs)
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

    # vascular volumes (i.e. plasma volume), [L]; homo, 71kg
    global V_V = zeros(N_Organs)
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

    # Interstitial Volumes, [L]; homo, 71kg
    global V_IntS = zeros(N_Organs)
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

    # Endosomal Volumes, [L]; homo, 71kg
    global V_endosomal = zeros(N_Organs)
    V_endosomal[Heart] = 1.71/1000
    V_endosomal[Lung] = 5.0/1000
    V_endosomal[Muscle] = 150.0/1000
    V_endosomal[Skin] = 17.0/1000
    V_endosomal[Adipose] = 67.3/1000
    V_endosomal[Bone] = 50.8/1000
    V_endosomal[Brain] = 7.25/1000
    V_endosomal[Kidney] = 1.66/1000
    V_endosomal[Liver] = 10.7/1000
    V_endosomal[SI] = 1.93/1000
    V_endosomal[LI] = 2.74/1000
    V_endosomal[Pancreas] = 0.518/1000
    V_endosomal[Thymus] = 0.0321/1000
    V_endosomal[Spleen] = 1.11/1000
    V_endosomal[Other] = 24.3/1000

    # Blood Flows, [L/h]; homo, 71kg
    global PLQ = zeros(N_Organs)
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
    # rebalance PLQ
    PLQ[Lung] = sum(PLQ[2:end])

    # Lymph Flows, [L/h];
    LF = zeros(N_Organs)
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
    
    # tissue lymph flow through large and small pores
    alpha_L = 0.042; # Fractional hydraulic conductance of large pores
    alpha_S = 0.942; # Fractional hydraulic conductance of small pores
    r_S = 4.44; # small pore radius, [nm]
    r_L = 22.85; # large pore radius, [nm]
    x_j = 0.38; # derived in appendix
    global J_iso = x_j * LF; 
    LF_L = J_iso .+ alpha_L*LF
    LF_S = -J_iso .+ alpha_S*LF
    X_p = 13197. ; # [nm^3] # see appendix
    X_P_S = X_p * 1/a_e * ratio_As_A0s * (1-alpha_L)/(r_S*r_S);
    X_P_L = X_p * 1/a_e * ratio_Al_A0l * alpha_L/(r_L*r_L);
    global PS_S = X_P_S * LF; 
    global PS_L = X_P_L * LF; 
    # Pe_L, Pe_S could be calculated based on other elements
    Pe_L = unique(round.(LF_L*(1-sigma_L) ./ PS_L, sigdigits=4))[1]
    Pe_S = unique(round.(LF_S*(1-sigma_S) ./ PS_S, sigdigits=4))[1]
    CL_tp_L = PS_L .* ( 1 .- C_EXG[1:N_Organs,IS]./C_EXG[1:N_Organs,V] ) * Pe_L/(exp(Pe_L)-1) .+ LF_L*(1-sigma_L);
    CL_tp_S = PS_S .* ( 1 .- C_EXG[1:N_Organs,IS]./C_EXG[1:N_Organs,V] ) * Pe_S/(exp(Pe_S)-1) .+ LF_S*(1-sigma_S);
    # incorporate conditions when drug conc = 0
    for k in 1:length(CL_tp_L)
        if isnan(CL_tp_L[k])
            CL_tp_L[k] = 0
        end
        if isnan(CL_tp_S[k])
            CL_tp_S[k] = 0
        end
    end
    CLtp = CL_tp_L .+ CL_tp_S

    # endothelial pinocytosis rate, [L/h]; (reported number was the rate of pinocytosis per endosomal space)
    CL_up_per_endosomal = 0.55; # [L/h/L]
    CLup = CL_up_per_endosomal * V_endosomal; 

    # reflection coefficients
    sigma_IS = 0.2 * ones(N_Organs);

    # First order lysosomal degradation rate; mouse, 28g
    # k_deg = 32.2; # [hr-1]

    # GFR; assuming human BSA = 1.9m^2
    GFR = 12/1000*60*1.9/1.73 ; # [L/h]; converted from mL/min/1.73m^2

    # ===================================================================================== # 
    #  ---------------------- DIFFERENTIAL VARIABLES/EQUATIONS ---------------------------- #
    # ===================================================================================== #

    # generaic plasma concentration
    @views @. du.C_EXG[3:N_Organs,V] = ((
        PLQ[3:N_Organs] * C_EXG[Lung, V] 
        - (PLQ[3:N_Organs] - LF[3:N_Organs]) * C_EXG[3:N_Organs, V]  
        - (CLtp[3:N_Organs] + CLup[3:N_Organs]) * C_EXG[3:N_Organs, V] 
    )/V_V[3:N_Organs]); 

    # liver plasma concentration 
    du.C_EXG[Liver, V] = (
        (PLQ[Liver]*C_EXG[Lung,V] + 
        (PLQ[Spleen]-LF[Spleen])*C_EXG[Spleen,V] + (PLQ[Pancreas]-LF[Pancreas])*C_EXG[Pancreas,V] + (PLQ[SI]-LF[SI])*C_EXG[SI,V] + (PLQ[LI]-LF[LI])*C_EXG[LI,V] 
        - C_EXG[Liver,V]*(PLQ[Liver]-LF[Liver] + PLQ[Spleen]-LF[Spleen] + PLQ[Pancreas]-LF[Pancreas] + PLQ[SI]-LF[SI]+ PLQ[LI]-LF[LI]) 
        - (CLtp[Liver] + CLup[Liver])*C_EXG[Liver, V]
        )/V_V[Liver]);

    # lung plasma concentration 
    du.C_EXG[Lung, V] = ((
        (PLQ[Lung]+LF[Lung])*C_EXG_Plasma - PLQ[Lung]*C_EXG[Lung,V] - (CLtp[Lung] + CLup[Lung]) * C_EXG[Lung, V] 
    )/V_V[Lung] );

    # update renal clearance  
    global CL_renal = GFR * theta; 
    du.C_EXG[Kidney, V] = du.C_EXG[Kidney, V] - CL_renal*C_EXG[Kidney, V]/V_V[Kidney]; 

    w_other_1_EXG_Plasma = (
      (PLQ[Heart]-LF[Heart])*C_EXG[Heart,V] 
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

    # endosomal concentration
    @views @. du.C_EXG[1:N_Organs,E] = ((
        CLup[1:N_Organs] * (C_EXG[1:N_Organs, V] + C_EXG[1:N_Organs, IS]) - k_deg * C_EXG[1:N_Organs, E] * V_endosomal[1:N_Organs]
    )/V_endosomal[1:N_Organs]);

    # Interstitial concentration 
    @views @. du.C_EXG[1:N_Organs,IS] = ((
        CLtp[1:N_Organs] * C_EXG[1:N_Organs, V] - (1 .- sigma_IS) * LF[1:N_Organs] * C_EXG[1:N_Organs,IS] 
        - CLup[1:N_Organs] *  C_EXG[1:N_Organs, IS] 
    #    - CLup[1:N_Organs] *  C_EXG[1:N_Organs, V]  # removed based on mass balance
    )/V_IntS[1:N_Organs]);

    # LN concentration 
    du.C_EXG_LN = ((
      (1-sigma_IS[Heart])*LF[Heart]*C_EXG[Heart, IS]
    + (1-sigma_IS[Kidney])*LF[Kidney]*C_EXG[Kidney, IS]
    + (1-sigma_IS[Muscle])*LF[Muscle]*C_EXG[Muscle, IS]
    + (1-sigma_IS[Skin])*LF[Skin]*C_EXG[Skin, IS]
    + (1-sigma_IS[Brain])*LF[Brain]*C_EXG[Brain, IS]
    + (1-sigma_IS[Adipose])*LF[Adipose]*C_EXG[Adipose, IS]
    + (1-sigma_IS[Thymus])*LF[Thymus]*C_EXG[Thymus, IS]
    + (1-sigma_IS[Liver])*LF[Liver]*C_EXG[Liver, IS]
    + (1-sigma_IS[Spleen])*LF[Spleen]*C_EXG[Spleen, IS]
    + (1-sigma_IS[Pancreas])*LF[Pancreas]*C_EXG[Pancreas, IS]
    + (1-sigma_IS[SI])*LF[SI]*C_EXG[SI, IS]
    + (1-sigma_IS[LI])*LF[LI]*C_EXG[LI, IS]
    + (1-sigma_IS[Bone])*LF[Bone]*C_EXG[Bone, IS]
    + (1-sigma_IS[Other])*LF[Other]*C_EXG[Other, IS]
    + (1-sigma_IS[Lung])*LF[Lung]*C_EXG[Lung, IS]
        - L_LymphNode * C_EXG_LN
    )/V_LN);

    
    # tumor dynamicss
    du.Nc1 = 0.
    du.Nc2 = 0.
    du.Nc3 = 0.
    du.Nc4 = 0.
    
    Vtumor = (Nc1 + Nc2 + Nc3 + Nc4) * Vc / 0.375; # tumor volume, [L]
    Rtumor = ((Vtumor*1E3)/(4/3*pi))^(1/3) ; # tumor radius, [cm]
    V_IntS_Tumor = 0.55 * Vtumor
    
    surface_exchange = 6 * D_protein/ (Rtumor^2) ; # [h-1]
    vascular_exchange = 2 * P_protein * Rcap / (Rkrogh^2) ; # [h-1]
    
    flux_ADC_plasma2tumor = (vascular_exchange + surface_exchange) * (C_EXG_Plasma * epsilon - C_EXG_Tumor) * V_IntS_Tumor
    
    du.C_EXG_Tumor = flux_ADC_plasma2tumor/V_IntS_Tumor;
    
    du.C_EXG_Plasma = (w_other_1_EXG_Plasma - flux_ADC_plasma2tumor)/V_Plasma + infusion; 
    
    # tracking degredation and renal clearance 
    du.endo_deg = sum(k_deg * C_EXG[1:N_Organs, E] .* V_endosomal[1:N_Organs]); 
    du.renal_clearance = CL_renal*C_EXG[Kidney, V];

end

