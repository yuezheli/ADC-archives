# date: 9/24/2025 
# author: Yuezhe Li 
# purpose of this code: to create parameters for pbpk 
## create dictionary with parameter values

using ProjectRoot

# constant function
include(@projectroot("model/constants.jl"));


function get_p_map(model, organ, parameter, value)
    # Use lowercase for case-insensitive matching (e.g., "icb" matches both "eye₊icb₊..." and "sigma_IS_ICB")
    organ_lower = lowercase(organ)
    param_array = filter(var -> occursin(organ_lower, lowercase(string(var))) && occursin(parameter, string(var)), parameters(model))
    p_map = []
    for x in param_array
        push!(p_map, Pair(x, value))
    end
    return vcat(p_map...)
end

function get_p_map_all(model, p_dict)
    param_map = []
    organs = keys(p_dict)

    for x in organs
        parval = p_dict[x]
        organ_map = []
        for j in 1:length(parval)
            param = string(keys(parval)[j])
            value = values(parval)[j]
            param_new = get_p_map(model, x, param, value)
            organ_map = push!(organ_map, param_new)
        end
        organ_map_flat = vcat(organ_map...)
        param_map = push!(param_map, organ_map_flat)
    end
    param_map_flat = vcat(param_map...)

    return param_map_flat
end


function create_base_pbpk_param(PS_Score, mdl; k_PL_ints_clearance = 0,  V_Plasma = 3.126)

    # additional eye things 
    f_cho = 0.81 
    f_ret = 0.05 
    f_icb = 1 - f_cho - f_ret

    ## general
    PS_a = 1.8051 # [-]
    PS_b = 0.2624 # [-]
    # PS_Score = 0.0 # [-]  # polyspecificity score
    PS_Kd = 10^(exp(PS_a - PS_b * PS_Score))
    k_on_PS = 8.06e7 / 1e6 # [1/(uM*h)]  # polyspecificity on-rate
    k_off_PS = (PS_Kd / 1000.0) * k_on_PS
    # V_Plasma = 3.126; # [L]. value for human
    KD6_WT =  700.0    # [nM]  
    KD7_WT =  154077.0 # [nM]  
    # Blood flows
    # PLQ_lung    = 181913 / 1000.0 # [L/h]  # should check for rebalancing 
    PLQ_heart   = 7752 / 1000.0  # [L/h]
    PLQ_muscle  = 33469 / 1000.0 # [L/h]
    PLQ_skin    = 11626 / 1000.0 # [L/h]
    PLQ_adipose = 11233 / 1000.0 # [L/h]
    PLQ_bone    = 2591 / 1000.0  # [L/h]
    PLQ_brain   = 21453 / 1000.0 # [L/h]
    PLQ_kidney  = 36402 / 1000.0 # [L/h]
    PLQ_SI      = 12368 / 1000.0 # [L/h]
    PLQ_LI      = 12867 / 1000.0 # [L/h]
    PLQ_pancreas= 3056 / 1000.0  # [L/h]
    PLQ_thymus  = 353 / 1000.0   # [L/h]
    PLQ_spleen  = 6343 / 1000.0 # [L/h]
    PLQ_Marrow = 8280/1000
    PLQ_Eye     = 0.06 # [L/h]
    PLQ_ICB     = PLQ_Eye * f_icb
    PLQ_Retina  = PLQ_Eye * f_ret
    PLQ_Choroid = PLQ_Eye * f_cho
    PLQ_other   = 9190 / 1000.0 - PLQ_Eye - PLQ_Marrow # [L/h]
    PLQ_liver   = 13210/1000 # [L/h] 
    PLQ_lung = PLQ_heart + PLQ_muscle + PLQ_skin + PLQ_adipose + PLQ_bone + PLQ_brain + PLQ_kidney + PLQ_SI + PLQ_LI + PLQ_pancreas + PLQ_thymus + PLQ_spleen + PLQ_Eye + PLQ_Marrow + PLQ_other + PLQ_liver # rebalance lung PLQ
    ## lymph flows [L/h]
    LF_Lung     = PLQ_lung * 0.002
    LF_Heart    = PLQ_heart * 0.002
    LF_Muscle   = PLQ_muscle * 0.002
    LF_Skin     = PLQ_skin * 0.002
    LF_Adipose  = PLQ_adipose * 0.002
    LF_Bone     = PLQ_bone * 0.002
    LF_Brain    = PLQ_brain * 0.002
    LF_Kidney   = PLQ_kidney * 0.002
    LF_SI       = PLQ_SI * 0.002
    LF_LI       = PLQ_LI * 0.002
    LF_Pancreas = PLQ_pancreas * 0.002
    LF_Thymus   = PLQ_thymus * 0.002
    LF_Spleen   = PLQ_spleen * 0.002
    LF_Other    = PLQ_other * 0.002
    LF_Marrow    = PLQ_Marrow * 0.002
    LF_Liver    = (PLQ_liver + PLQ_SI - LF_SI + PLQ_LI - LF_LI + PLQ_spleen - LF_Spleen + PLQ_pancreas - LF_Pancreas) * 0.002 
    LF_ICB      = PLQ_ICB * 0.002 
    LF_Retina   = PLQ_Retina * 0.002 
    LF_Choroid  = PLQ_Choroid * 0.002
    L_LymphNode = LF_Marrow + LF_ICB + LF_Retina + LF_Choroid + LF_Lung + LF_Heart + LF_Muscle + LF_Skin + LF_Adipose + LF_Bone + LF_Brain + LF_Kidney + LF_SI + LF_LI + LF_Pancreas + LF_Thymus + LF_Spleen + LF_Other + LF_Liver # [L/h]


    p_dict = Dict(
        "lung" => (Endothelial_Cell_Frac = 0.0834, PLQ = PLQ_lung, sigma_V = 0.95, V_IntS = 300/1000, V_V=  55/1000, LF=LF_Lung, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "liver" => (Endothelial_Cell_Frac = 0.1877, PLQ = PLQ_liver, sigma_V = 0.85, V_IntS = 429/1000, V_V = 183/1000, LF=LF_Liver, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "heart" => (Endothelial_Cell_Frac = 0.0011, PLQ = PLQ_heart, sigma_V = 0.95, V_IntS = 48.8/1000, V_V = 13.1/1000, LF=LF_Heart, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "muscle" => (Endothelial_Cell_Frac = 0.1928, PLQ = PLQ_muscle, sigma_V = 0.95, V_IntS = 3910/1000, V_V = 662/1000, LF=LF_Muscle, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "skin" => (Endothelial_Cell_Frac = 0.0819, PLQ = PLQ_skin, sigma_V = 0.95, V_IntS = 1125/1000, V_V = 127/1000, LF=LF_Skin, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "adipose" => (Endothelial_Cell_Frac = 0.0999, PLQ = PLQ_adipose, sigma_V = 0.95, V_IntS = 2289/1000, V_V = 148/1000, LF=LF_Adipose, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "bone" => (Endothelial_Cell_Frac = 0.1478, PLQ = PLQ_bone, sigma_V = 0.85, V_IntS = 1891/1000, V_V = 224/1000, LF=LF_Bone, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "brain" => (Endothelial_Cell_Frac = 0.0115, PLQ = PLQ_brain, sigma_V = 0.99, V_IntS = 261/1000, V_V = 31.9/1000, LF=LF_Brain, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "kidney" => (Endothelial_Cell_Frac = 0.0157, PLQ = PLQ_kidney, sigma_V = 0.9, V_IntS = 49.8/1000, V_V = 18.2/1000, LF=LF_Kidney, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "sm_int" => (Endothelial_Cell_Frac = 0.0121, PLQ = PLQ_SI, sigma_V = 0.9, V_IntS = 67.1/1000, V_V = 6.15/1000, LF=LF_SI, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "la_int" => (Endothelial_Cell_Frac = 0.0209, PLQ = PLQ_LI, sigma_V = 0.95, V_IntS = 95.3/1000, V_V = 8.74/1000, LF=LF_LI, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "pancreas" => (Endothelial_Cell_Frac = 0.0001, PLQ = PLQ_pancreas, sigma_V = 0.9, V_IntS = 18/1000, V_V = 5.7/1000, LF=LF_Pancreas, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "thymus" => (Endothelial_Cell_Frac = 0.0001, PLQ = PLQ_thymus, sigma_V = 0.9, V_IntS = 1.09/1000, V_V = 0.353/1000, LF=LF_Thymus, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "spleen" => (Endothelial_Cell_Frac = 0.0499, PLQ = PLQ_spleen, sigma_V = 0.85, V_IntS = 44.3/1000, V_V = 26.8/1000, LF=LF_Spleen, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "other" => (Endothelial_Cell_Frac = 0.0951 - 1E-7*2 - 1E-8 - 0.0011, PLQ = PLQ_other, sigma_V = 0.95, V_IntS = 831/1000, V_V = 204/1000, LF=LF_Other, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance), # adjusted to remove endothelial cells in the eye & bone marrow
        "marrow" => (Endothelial_Cell_Frac = 0.0011, PLQ = PLQ_Marrow, sigma_V = 0.9, V_IntS = 279/1000, V_V = 150/1000, LF=LF_Marrow, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        # eye components - string matching will find both igg_exg and igg_edg automatically
        "icb" => (Endothelial_Cell_Frac = 1E-7, PLQ = PLQ_ICB, sigma_V = 0.99, V_IntS = 1.38E-4 * 0.34, V_V = 1.38E-4 * 0.23, LF = LF_ICB, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "retina" => (Endothelial_Cell_Frac = 1E-8, PLQ = PLQ_Retina, sigma_V = 0.99, V_IntS = 3.26E-4 * 0.32, V_V = 3.26E-4 * 0.28, LF = LF_Retina, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "choroid" => (Endothelial_Cell_Frac = 1E-7, PLQ = PLQ_Choroid, sigma_V = 0.99, V_IntS = 1.39E-4 * 0.34, V_V = 1.39E-4 * 0.23, LF = LF_Choroid, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
    )

    ## create parameter map

    p_map = Dict(get_p_map_all(mdl, p_dict))
    p_map_global = Dict([
        mdl.CL_up_in_nL_per_hour_per_million_cells => 150.0,
        mdl.Total_Endothelial_Cell => 1.422e9,  # mouse; 1.422e+009 
        mdl.Scale_Factor => 603.7, # scales number of endothelial cells to value for human (1 for mouse)
        mdl.tau_VM => 1/60.0,  # [h] 
        mdl.tau_ISM => 1/60, # [h]
        mdl.FR => 0.715,  # fraction recycled
        mdl.k_on_PS => 8.06e7/1e6,
        mdl.kint_PS => 0.0380,  # [1/h]. Internalization rate for polyspecificity-membrane bound mAb
        mdl.Prob_deg => 0.95, # [-] ~ 95% of free mAb inside endosome will be routed to be degraded
        mdl.E6a_Vol_Pct => 0.33, # [-]
        mdl.pino_time => 10.8/60,
        mdl.kdeg_FcRn_Ab => log(2)/11.1,
        mdl.k_on_6_EXG => 8.06e7 / 1e6,  # [1/(uM*h)]
        mdl.k_on_7_EXG => 1.61e7 / 5 / 1e6,  # [1/(uM*h)]
        mdl.k_on_6_EDG => 8.06e7 / 1e6,  # [1/(uM*h)]
        mdl.k_on_7_EDG => 1.61e7 / 5 / 1e6,  # [1/(uM*h)]
        mdl.KD_6_EDG => KD6_WT,    # [nM]   
        mdl.KD_7_EDG => KD7_WT, # [nM]

        mdl.on_rate_ratio_1st_vs_2nd_binding => 83.7, # [-]
        mdl.FcRn_recycle_fraction => 0.99,

        mdl.C_Mem => 18.5, # [uM]  

        mdl.k_off_PS => k_off_PS,
        mdl.L_LymphNode => L_LymphNode, # [L/h]
        mdl.V_Plasma => V_Plasma, # [L]
        mdl.V_LN     => 0.274,
        mdl.KD_6_EXG => KD6_WT,    # [nM]
        mdl.KD_7_EXG => KD7_WT, # [nM]
        mdl.CLOff => 1.0, 

        # infusion
        mdl.infusion => 0.,

        # added for eyes
        mdl.RC_aq => 0.99, 
        mdl.RC_ret => 0.95, 
        mdl.RC_cho => 0.99,

        mdl.Q_AH => 9.24E-5, 
        mdl.Q_PtA => 1.5E-4, 
        mdl.Q_BF => 9.84E-5, 

        mdl.PS_ret => 8.34E-6, 
        mdl.PS_cho => 3.21E-6, 
        mdl.PS_cor => 3.95E-8, 
        mdl.PS_lens => 0, 

        mdl.V_AQ => 3.1E-4,
        mdl.V_VH => 4E-3, 
        mdl.V_Lens => 2.2E-4, 
        mdl.V_Cor => 7.7E-5, 

        # added for ADC PL
        mdl.DAR => 2, 
        mdl.frac_lys => 1, 
        mdl.k_diff => 3 * 16E-6 / cell_radius * s_per_hr , 
        mdl.CL_plasma_PL => 0,  # PL clearance rate in plasma 
        mdl.Kp_PL_endo_plasma => 1, # partition coefficient, endothelial cells:plasma

        # assumed PL permeability surface area in the eye; computed = permeability * surface area 
        mdl.PS_cor_pl => 16E-6 * 126E-6 * s_per_hr,   # [L/hr]
        mdl.PS_ret_pl => 16E-6 * 1363E-6 * s_per_hr,    # [L/hr]
        mdl.PS_lens_pl => 0,  # Abbvie rat study
        
        # assumed PL partition coefficients in the eye 
        mdl.Kp_ah_pl => 1, 
        mdl.Kp_vh_pl => 1, 

        # assumed ADC partition coefficient between aqueous humor and cornea 
        mdl.Kp_adc_cor_ah => 1, 

        # added nonlinear clearance of ADC in plasma 
        mdl.VMAX => 0, 
        mdl.KM => 1, 
        mdl.k_deconj => 0, 
    ])

    p_map_all = merge(p_map, p_map_global); 

    return p_map_all

end
