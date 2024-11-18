# author: Yuezhe Li 
# date: 9/26/2023

using ComponentArrays
using Parameters: @unpack

include("param.jl") # load this for number of organs in the model 

# create initial value for human status
function jones_init(N_Organs, Dose_in_mgkg, BW, V_Plasma, init_sR = 0.)

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
    C_EXG_Plasma0 =  Dose_in_mgkg*1E3*BW/(V_Plasma)/MW_EDG; # [uM]

    init_DEGprotein = ComponentArray(deg_EXG = zeros(N_Organs), deg_EDG = zeros(N_Organs), deg_plasma = 0.0)


    u0 = ComponentArray(C_EXG_Plasma = C_EXG_Plasma0, 
                            C_sR_plasma = init_sR, C_sR_EXG = 0., 
                            C_EDG_Plasma = C_EDG_Plasma0, C_EDG_LN = C_EDG_LN0, C_EDG = C_EDG0, C_EXG_LN = C_EXG_LN0,
                            C_EXG = C_EXG0, C_FcRn_E6a = C_FcRn_E6a0, C_FcRn_E7 = C_FcRn_E70, C_FcRn_E7b = C_FcRn_E7b0, 
                            C_FcRn_ISM = C_FcRn_ISM0, C_FcRn_VM = C_FcRn_VM0, DEGprotein = init_DEGprotein ); 

    return u0, Dose_in_mgkg*1E3*BW/MW_EDG ;
end

