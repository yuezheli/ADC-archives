# author: Yuezhe Li 
# date: Nov 16, 2023

using ComponentArrays
using Parameters: @unpack

# nanobody parameter, 13kDa
p_nanobody = ComponentArray(
    a_e = 1.87,
    sigma_S = 0.495, 
    sigma_L = 0.0312, 
    ratio_As_A0s = 8.08E-2, 
    ratio_Al_A0l = 0.693,
    Pe_S = 0.0105, 
    Pe_L = 1.04, 
    theta = 0.667, 
    k_deg = 32.2, # [hr-1]
    infusion = 0.
);

# scFv parameters, 27kDa
p_scFv = ComponentArray(
    a_e = 2.48,
    sigma_S = 0.712, 
    sigma_L = 0.0526, 
    ratio_As_A0s = 2.5E-2, 
    ratio_Al_A0l = 0.599,
    Pe_S = 0.0257, 
    Pe_L = 1.59, 
    theta = 0.403, 
    k_deg = 32.2, # [hr-1]
    infusion = 0.
);

# Fab/ diabody parameters, 50kDa
p_Fab = ComponentArray(
    a_e = 3.15,
    sigma_S = 0.885, 
    sigma_L = 0.0819, 
    ratio_As_A0s = 3.72E-3, 
    ratio_Al_A0l = 0.533,
    Pe_S = 0.0876, 
    Pe_L = 2.15, 
    theta = 0.131, 
    k_deg = 32.2, # [hr-1]
    infusion = 0.
);

# scFv2 parameters, 55kDa
p_scFv2 = ComponentArray(
    a_e = 3.26,
    sigma_S = 0.906, 
    sigma_L = 0.0877, 
    ratio_As_A0s = 2.45E-3, 
    ratio_Al_A0l = 0.522,
    Pe_S = 0.113, 
    Pe_L = 2.27, 
    theta = 0.0986, 
    k_deg = 32.2, # [hr-1]
    infusion = 0.
);

# minibody parameter, 80kDa
p_minibody = ComponentArray(
    a_e = 3.77,
    sigma_S = 0.965, 
    sigma_L = 0.115, 
    ratio_As_A0s = 3.09E-4, 
    ratio_Al_A0l = 0.469,
    Pe_S = 0.381, 
    Pe_L = 2.83, 
    theta = 0.022, 
    k_deg = 32.2, # [hr-1]
    infusion = 0.
);

# Fab2 parameter, 100kDa
p_Fab2 = ComponentArray(
    a_e = 4.11,
    sigma_S = 0.984, 
    sigma_L = 0.135, 
    ratio_As_A0s = 5.87E-5, 
    ratio_Al_A0l = 0.431,
    Pe_S = 0.981, 
    Pe_L = 3.28, 
    theta = 0.00703, 
    k_deg = 32.2, # [hr-1]
    infusion = 0.
);

# IgG parameters, 150kDa
p_IgG = ComponentArray(
    a_e = 4.81,
    sigma_S = 0.998, 
    sigma_L = 0.18, 
    ratio_As_A0s = 9.28E-7, 
    ratio_Al_A0l = 0.349,
    Pe_S = 9.82, 
    Pe_L = 4.48, 
    theta = 0.0011, 
    k_deg = 32.2, # [hr-1]
    infusion = 0.
);

# initial condition 
function lishah_init(N_Organs, Dose_in_mgkg, BW, V_Plasma, MW)
    # MW in Dalton
    C_EXG_Plasma0 = Dose_in_mgkg*1E6*BW/MW/V_Plasma; # [nM]
    u0 = ComponentArray(C_EXG_Plasma = C_EXG_Plasma0, 
                        C_EXG_LN = 0., C_EXG = zeros(N_Organs, 3), 
                        endo_deg = 0., renal_clearance = 0.);
    return u0, Dose_in_mgkg*1E6*BW/MW; 
end


