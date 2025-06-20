# author: Yuezhe LI
# date: 6/8/23
# purpose of the code: 
# 1. to create initial condition 
# 2. to post-process the solution 

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using DataFrames, DataFramesMeta
using Parameters: @unpack

include(@projectroot("model/adc-constants.jl"));
include(@projectroot("model/param-pk.jl")); 

const N_Organs = 15

function jones_init(Dose_in_ugkg, p_base = p_base, bw = BW, V_plasma = V_Plasma)

    FcRn_Conc = 44.093; # [uM]  
    EDG_mg_ml = 10.0;  # [mg/mL]
    C_EDG_Plasma0 =  EDG_mg_ml/MW*1E6; # [uM]
    C_EDG_LN0 = 1E-18*1E6; # [uM]
    C_EDG0 =  1E-18*1E6 * ones(N_Organs,17); # [uM]
    C_EXG_LN0 = 1E-18*1E6; # [uM]
    C_EXG0 =  1E-18*1E6 *ones(N_Organs,19); # [uM]
    C_FcRn_E6a0 =  ones(N_Organs)*FcRn_Conc; # [uM]
    C_FcRn_E7b0 =  ones(N_Organs)*FcRn_Conc; # [uM] 
    C_FcRn_E70 =  ones(N_Organs)*FcRn_Conc; # [uM]      
    C_FcRn_ISM0 =  ones(N_Organs)*FcRn_Conc*1E-4; # [uM] 
    C_FcRn_VM0 =  ones(N_Organs)*FcRn_Conc*1E-4; # [uM] 
    C_EXG_Plasma0 =  Dose_in_ugkg*bw/(V_plasma)/MW; # [uM]

    init_DEGprotein = ComponentArray(deg_EXG = zeros(N_Organs), deg_EDG = zeros(N_Organs), deg_plasma = 0.0)
    init_DegPL = ComponentArray(deg_membrane = zeros(N_Organs), deg_FcRn = zeros(N_Organs))

    u0 = ComponentArray(C_EXG_Plasma = C_EXG_Plasma0, 
        C_sR_plasma = 0., C_sR_EXG = 0., 
        C_EDG_Plasma = C_EDG_Plasma0, C_EDG_LN = C_EDG_LN0, C_EDG = C_EDG0, C_EXG_LN = C_EXG_LN0,
        C_EXG = C_EXG0, C_FcRn_E6a = C_FcRn_E6a0, C_FcRn_E7 = C_FcRn_E70, C_FcRn_E7b = C_FcRn_E7b0, 
        C_FcRn_ISM = C_FcRn_ISM0, C_FcRn_VM = C_FcRn_VM0, tumor_ADC = 0.0, DEGprotein = init_DEGprotein, DegPL = init_DegPL, 
        end_cyto_payload = zeros(N_Organs), ints_payload = zeros(N_Organs)); 
    
    # update initial soluble receptor concentration if existed
    try 
        u0.C_sR_plasma = p_base.init_sR
    catch 
        println("no soluble receptor specified")
    end

    return u0, C_EXG_Plasma0;
end
