# author: Yuezhe Li 
# date: Nov 16, 2023
# purpose of the script: to build the model described in Li and Shah, 2019. 
# https://pubmed.ncbi.nlm.nih.gov/31028591/

using Parameters: @unpack 
using ComponentArrays

function lishah_mus!(du, u, p, t)
    @unpack C_EXG_Plasma, C_EXG_LN, C_EXG, endo_deg, renal_clearance = u

    # a_e = 0.0483 * MW^0.386; # a_e in nm, MW in Dalton
    # theta = exp(1 - 8.7/(1 + exp(0.028*(-MW+72.3)))); # MW in kDa
    # ratio_As_A0s = 0.2352*exp(-0.00008295*MW) + 0.7767*exp(-0.00053095*MW); # MW in Dalton
    # ratio_Al_A0l = 0.3429*exp(-0.00012175*MW) + 0.6571*exp(-0.00000421*MW); MW in Dalton
    # sigma_L = 0.000035*MW^0.717; # MW in Dalton
    # sigma_S = 1-0.8489*exp(-0.00004*MW); # MW in Dalton
    # Pe_S, Pe_L calculation see below; does not depend on additional input

    @unpack a_e, theta, k_deg, ratio_Al_A0l, ratio_As_A0s, sigma_L, sigma_S, Pe_S, Pe_L, infusion = p

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

    # Organ Volumes, [L]; mouse, 28g
    global V_Plasma = 0.944/1000 
    global V_LN     = 0.113/1000

    # Organ Volumes, [L]; mouse, 28g
    global V_Organ = zeros(N_Organs)
    V_Organ[Lung] = 0.204/1000
    V_Organ[Liver] = 1.93/1000
    V_Organ[Heart] = 0.152/1000
    V_Organ[Muscle] = 11.3/1000
    V_Organ[Skin] = 5.02/1000
    V_Organ[Adipose] = 1.98/1000
    V_Organ[Bone] = 2.82/1000
    V_Organ[Brain] = 0.485/1000
    V_Organ[Kidney] = 0.525/1000
    V_Organ[SI] = 0.728/1000
    V_Organ[LI] = 0.314/1000
    V_Organ[Pancreas] = 0.097/1000
    V_Organ[Thymus] = 0.009/1000
    V_Organ[Spleen] = 0.127/1000
    V_Organ[Other] = 0.465/1000

    # vascular volumes (i.e. plasma volume), [L]; mouse, 28g
    global V_V = zeros(N_Organs)
    V_V[Heart] = 0.00585/1000
    V_V[Lung] = 0.0295/1000
    V_V[Muscle] = 0.249/1000
    V_V[Skin] = 0.188/1000
    V_V[Adipose] = 0.0218/1000
    V_V[Bone] = 0.0621/1000
    V_V[Brain] = 0.0107/1000
    V_V[Kidney] = 0.0289/1000
    V_V[Liver] = 0.164/1000
    V_V[SI] = 0.0116/1000
    V_V[LI] = 0.005/1000
    V_V[Pancreas] = 0.00534/1000
    V_V[Thymus] = 0.0005/1000
    V_V[Spleen] = 0.0154/1000
    V_V[Other] = 0.0195/1000

    # Interstitial Volumes, [L]; mouse, 28g
    global V_IntS = zeros(N_Organs)
    V_IntS[Heart] = 0.0217/1000
    V_IntS[Lung] = 0.0384/1000
    V_IntS[Muscle] = 1.47/1000
    V_IntS[Skin] = 1.66/1000
    V_IntS[Adipose] = 0.337/1000
    V_IntS[Bone] = 0.525/1000
    V_IntS[Brain] = 0.0873/1000
    V_IntS[Kidney] = 0.0788/1000
    V_IntS[Liver] = 0.385/1000
    V_IntS[SI] = 0.127/1000
    V_IntS[LI] = 0.0545/1000
    V_IntS[Pancreas] = 0.0169/1000
    V_IntS[Thymus] = 0.00153/1000
    V_IntS[Spleen] = 0.0254/1000
    V_IntS[Other] = 0.0797/1000

    # Endosomal Volumes, [L]; mouse, 28g
    global V_endosomal = zeros(N_Organs)
    V_endosomal[Heart] = 0.00076/1000
    V_endosomal[Lung] = 0.00102/1000
    V_endosomal[Muscle] = 0.0566/1000
    V_endosomal[Skin] = 0.0251/1000
    V_endosomal[Adipose] = 0.00991/1000
    V_endosomal[Bone] = 0.0141/1000
    V_endosomal[Brain] = 0.00243/1000
    V_endosomal[Kidney] = 0.00263/1000
    V_endosomal[Liver] = 0.00963/1000
    V_endosomal[SI] = 0.00364/1000
    V_endosomal[LI] = 0.00157/1000
    V_endosomal[Pancreas] = 0.000485/1000
    V_endosomal[Thymus] = 0.00005/1000
    V_endosomal[Spleen] = 0.000635/1000
    V_endosomal[Other] = 0.00233/1000

    # Blood Flows, [L/h]; mouse, 28g
    PLQ = zeros(N_Organs)
    PLQ[Heart] = 36.5/1000
    PLQ[Lung] = 373. /1000
    PLQ[Muscle] = 86.1/1000
    PLQ[Skin] = 27.8/1000
    PLQ[Adipose] = 13.4/1000
    PLQ[Bone] = 15.2/1000
    PLQ[Brain] = 11.8/1000
    PLQ[Kidney] = 68.5/1000
    PLQ[Liver] = 10.3/1000
    PLQ[SI] = 58.1/1000
    PLQ[LI] = 17.3/1000
    PLQ[Pancreas] = 6.24/1000
    PLQ[Thymus] = 1.19/1000
    PLQ[Spleen] = 8.18/1000
    PLQ[Other] = 10.9/1000
    PLQ[Lung] = sum(PLQ[2:end]); # this number matched what was provided in the Li and Shah 2019 paper;

    # Lymph Flows, [L/h]; mouse, 28g
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
    L_LymphNode = sum(LF) # [L/h] # note this number is the same as the number listed in Table 1 in Li and Shah, 2019;
    
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
    # Pe_L = unique(round.(LF_L*(1-sigma_L) ./ PS_L, sigdigits=4))[1]
    # Pe_S = unique(round.(LF_S*(1-sigma_S) ./ PS_S, sigdigits=4))[1]
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

    # GFR; mouse, 28g
    GFR = 0.278/1000*60; # [L/h]

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

    du.C_EXG_Plasma = (w_other_1_EXG_Plasma/V_Plasma) + infusion; 

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

    # tracking degredation and renal clearance 
    du.endo_deg = sum(k_deg * C_EXG[1:N_Organs, E] .* V_endosomal[1:N_Organs]); 
    du.renal_clearance = CL_renal*C_EXG[Kidney, V];

end
