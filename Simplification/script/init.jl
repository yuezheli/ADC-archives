# author:Yuezhe Li 
# date: 2/29/2024
# purpose of this script: to set initial condition for PBPK model

function pbpk_init(N_Organs, TotalCells, RposFrac, Dose_in_ugkg, Rcopy, init_sR, k_endo, k_rec)

    InitRpos = TotalCells * RposFrac
    InitRneg = TotalCells * (1-RposFrac)

    MW_EDG = 150000.0; # [Da]
    FcRn_Conc = 44.093; # [uM]  
    EDG_mg_ml = 10.0;  # [mg/mL]
    C_EDG_Plasma0 =  EDG_mg_ml/MW_EDG*1E6; # [uM]
    C_EDG_LN0 = 1E-18*1E6; # [uM]
    C_EDG0 =  1E-18*1E6 * ones(N_Organs,17); # [uM]
    C_EXG_LN0 = 1E-18*1E6; # [uM]
    C_EXG0 =  1E-18*1E6 *ones(N_Organs,19); # [uM]
    C_FcRn_E6a0 =  ones(N_Organs)*FcRn_Conc; # [uM]
    C_FcRn_E7b0 =  ones(N_Organs)*FcRn_Conc; # [uM] 
    C_FcRn_E70 =  ones(N_Organs)*FcRn_Conc; # [uM]      
    C_FcRn_ISM0 =  ones(N_Organs)*FcRn_Conc*1E-4; # [uM] 
    C_FcRn_VM0 =  ones(N_Organs)*FcRn_Conc*1E-4; # [uM] 
    C_EXG_Plasma0 =  Dose_in_ugkg*BW/(V_Plasma)/MW_EDG; # [uM]

    init_tumor = ComponentArray(
        Nc_1 = InitRpos, Nc_2 = 0., Nc_3 = 0., Nc_4 = 0., 
        Nc_1_neg = InitRneg, Nc_2_neg = 0., Nc_3_neg = 0., Nc_4_neg = 0.,
        R_s = Rcopy/N_av*1e6, R_e = k_endo/k_rec*Rcopy/N_av*1e6, AR_s = 0., AR_e = 0., P_c = 0., P_m = 0., P_neg_c = 0., 
        A_m = 0.
        );
    
    init_DEGprotein = ComponentArray(deg_EXG = zeros(N_Organs), deg_EDG = zeros(N_Organs), deg_TumorCellular = 0.0, deg_plasma = 0.0)
    
    u0 = ComponentArray(C_EXG_Plasma = C_EXG_Plasma0, 
                        C_sR_plasma = init_sR, C_sR_EXG = 0., 
                        C_EDG_Plasma = C_EDG_Plasma0, C_EDG_LN = C_EDG_LN0, C_EDG = C_EDG0, C_EXG_LN = C_EXG_LN0,
                        C_EXG = C_EXG0, C_FcRn_E6a = C_FcRn_E6a0, C_FcRn_E7 = C_FcRn_E70, C_FcRn_E7b = C_FcRn_E7b0, 
                        C_FcRn_ISM = C_FcRn_ISM0, C_FcRn_VM = C_FcRn_VM0,tumor = init_tumor, DEGprotein = init_DEGprotein, 
                        end_endo_payload = zeros(N_Organs), end_cyto_payload = zeros(N_Organs), ints_payload = zeros(N_Organs)); 

    return u0, Dose_in_ugkg*BW/MW_EDG ;
end

