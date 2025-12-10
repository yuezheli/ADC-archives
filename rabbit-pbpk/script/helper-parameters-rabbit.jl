# date: 10/1-/2025 
# author: Yuezhe Li 
# purpose of this code: to create parameters for rabbit pbpk 
# source of parameters: Bussing and Shah, 2020; https://pubmed.ncbi.nlm.nih.gov/32876799/

using ProjectRoot

# constant function
include(@projectroot("model/constants.jl"));

function create_base_rabbit_pbpk_param(; PS_Score, Scale_Factor_Rabbit, 
    mdl, k_PL_ints_clearance = 0,  V_Plasma = 0.11)

    ## general (inherited from Jones 2019)
    PS_a = 1.8051 # [-]
    PS_b = 0.2624 # [-]
    # PS_Score = 0.0 # [-]  # polyspecificity score
    PS_Kd = 10^(exp(PS_a - PS_b * PS_Score))
    k_on_PS = 8.06e7 / 1e6 # [1/(uM*h)]  # polyspecificity on-rate
    k_off_PS = (PS_Kd / 1000.0) * k_on_PS
    # V_Plasma = 3.126; # [L]. value for human
    KD6_WT =  700.0    # [nM]  
    KD7_WT =  154077.0 # [nM]  

    # rabbit eye 
    f_cho = 0.86 
    f_ret = 0.00911 
    f_icb = 1 - f_cho - f_ret

    # Blood flows
    PLQ_heart   = 0.768  # [L/h]
    PLQ_kidney  = 3.69 # [L/h]
    PLQ_muscle  = 6.21 # [L/h]
    PLQ_skin    = 1.63 # [L/h]
    PLQ_adipose = 1.28 # [L/h]
    PLQ_bone    = 1.86  # [L/h]
    PLQ_brain   = 0.380 # [L/h]
    PLQ_SI      = 2.21 # [L/h]
    PLQ_LI      = 0.564 # [L/h]
    PLQ_pancreas= 0.144  # [L/h]
    PLQ_thymus  = 0.0224   # [L/h]
    PLQ_spleen  = 0.360 # [L/h]
    PLQ_Eye     = 0.0483 # [L/h]
    PLQ_ICB     = PLQ_Eye * f_icb
    PLQ_Retina  = PLQ_Eye * f_ret
    PLQ_Choroid = PLQ_Eye * f_cho
    PLQ_Marrow  = 0.768            # ASSUMED TO BE SIMILAR TO HEART 
    PLQ_other   = 1.89 - PLQ_Eye - PLQ_Marrow # [L/h]
    PLQ_liver   = 0.163 # [L/h] 
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
        "lung" => (Endothelial_Cell_Frac = 0.0834, PLQ = PLQ_lung, sigma_V = 0.95, V_IntS = 0.0032, V_V=  0.00247, LF=LF_Lung, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "heart" => (Endothelial_Cell_Frac = 0.0011, PLQ = PLQ_heart, sigma_V = 0.95, V_IntS = 0.0006, V_V = 0.000231, LF=LF_Heart, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "kidney" => (Endothelial_Cell_Frac = 0.0157, PLQ = PLQ_kidney, sigma_V = 0.9, V_IntS = 0.003, V_V = 0.00158, LF=LF_Kidney, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "muscle" => (Endothelial_Cell_Frac = 0.1928, PLQ = PLQ_muscle, sigma_V = 0.95, V_IntS = 0.162	, V_V = 0.0351, LF=LF_Muscle, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "brain" => (Endothelial_Cell_Frac = 0.0115, PLQ = PLQ_brain, sigma_V = 0.99, V_IntS = 5.60E-6, V_V = 0.000518, LF=LF_Brain, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "adipose" => (Endothelial_Cell_Frac = 0.0999, PLQ = PLQ_adipose, sigma_V = 0.95, V_IntS = 0.0162, V_V = 0.0012, LF=LF_Adipose, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "thymus" => (Endothelial_Cell_Frac = 0.0001, PLQ = PLQ_thymus, sigma_V = 0.9, V_IntS = 0.000697, V_V = 0.000228, LF=LF_Thymus, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "sm_int" => (Endothelial_Cell_Frac = 0.0121, PLQ = PLQ_SI, sigma_V = 0.9, V_IntS = 0.00564, V_V = 0.00144, LF=LF_SI, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "la_int" => (Endothelial_Cell_Frac = 0.0209, PLQ = PLQ_LI, sigma_V = 0.95, V_IntS = 0.00282, V_V = 0.00072, LF=LF_LI, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "spleen" => (Endothelial_Cell_Frac = 0.0499, PLQ = PLQ_spleen, sigma_V = 0.85, V_IntS = 0.00015, V_V = 0.000282, LF=LF_Spleen, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "pancreas" => (Endothelial_Cell_Frac = 0.0001, PLQ = PLQ_pancreas, sigma_V = 0.9, V_IntS = 0.000432, V_V = 0.000648, LF=LF_Pancreas, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "liver" => (Endothelial_Cell_Frac = 0.1877, PLQ = PLQ_liver, sigma_V = 0.85, V_IntS = 0.0163, V_V = 0.0115, LF=LF_Liver, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "bone" => (Endothelial_Cell_Frac = 0.1478, PLQ = PLQ_bone, sigma_V = 0.85, V_IntS = 0.0310, V_V = 0.0127, LF=LF_Bone, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "skin" => (Endothelial_Cell_Frac = 0.0819, PLQ = PLQ_skin, sigma_V = 0.95, V_IntS = 0.0332, V_V = 0.00209, LF=LF_Skin, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
        "other" => (Endothelial_Cell_Frac = 0.0951 - 0.005*3, PLQ = PLQ_other, sigma_V = 0.95, V_IntS = 0.0422, V_V = 0.0103, LF=LF_Other, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance), # adjusted to remove endothelial cells in the eye & bone marrow 
        "marrow" => (Endothelial_Cell_Frac = 0.0011, PLQ = PLQ_Marrow, sigma_V = 0.9, V_IntS = 279/1000, V_V = 150/1000, LF=LF_Marrow, sigma_IS = 0.2, k_out_ints = 0, k_in_ints = 0, k_PL_ex = k_PL_ints_clearance),
    )

    ## create parameter map

    p_map = Dict(get_p_map_all(mdl, p_dict))
    p_map_global = Dict([
        mdl.Total_Endothelial_Cell => 1.422e9,  # mouse 1.422e+009 
        # mdl.Scale_Factor => 603.7, # scales number of endothelial cells to value for human (1 for mouse)
        mdl.Scale_Factor => Scale_Factor_Rabbit, 

        # below are the same as Jones 2019
        mdl.CL_up_in_nL_per_hour_per_million_cells => 150.0,
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
        mdl.V_LN     => 0.0103,   # [L]
        mdl.KD_6_EXG => KD6_WT,    # [nM]
        mdl.KD_7_EXG => KD7_WT, # [nM]
        mdl.CLOff => 1.0, 

        # infusion 
        mdl.infusion => 0., 

        # eye parts that is parallel to other organs (ICB)
        mdl.eye.icb.Endothelial_Cell_Frac => 0.005, 
        mdl.eye.icb.PLQ => PLQ_ICB, 
        mdl.eye.icb.sigma_V => 0.99,    # Bussing and Shah, 2020; #  https://pubmed.ncbi.nlm.nih.gov/32876799/
        mdl.eye.icb.V_V => 87.8E-6 * 0.23, 
        mdl.eye.icb.V_IntS => 87.8E-6 * 0.4, 
        mdl.eye.icb.LF => LF_ICB, 
        mdl.eye.icb.sigma_IS => 0.2, 
        mdl.eye.icb.k_out_ints => 0, 
        mdl.eye.icb.k_in_ints => 0, 
        mdl.eye.icb.k_PL_ex => k_PL_ints_clearance, 

        # eye parts that is parallel to other organs (retina)
        mdl.eye.retina.Endothelial_Cell_Frac => 0.005, 
        mdl.eye.retina.PLQ => PLQ_Retina, 
        mdl.eye.retina.sigma_V => 0.99,    # Bussing and Shah, 2020; #  https://pubmed.ncbi.nlm.nih.gov/32876799/
        mdl.eye.retina.V_V => 42E-6 * 0.287, 
        mdl.eye.retina.V_IntS => 42E-6 * 0.32, 
        mdl.eye.retina.LF => LF_Retina, 
        mdl.eye.retina.sigma_IS => 0.2, 
        mdl.eye.retina.k_out_ints => 0, 
        mdl.eye.retina.k_in_ints => 0, 
        mdl.eye.retina.k_PL_ex => k_PL_ints_clearance, 

        # eye parts that is parallel to other organs (choroid)
        mdl.eye.choroid.Endothelial_Cell_Frac => 0.005, 
        mdl.eye.choroid.PLQ => PLQ_Choroid, 
        mdl.eye.choroid.sigma_V => 0.99,    # Bussing and Shah, 2020; #  https://pubmed.ncbi.nlm.nih.gov/32876799/
        mdl.eye.choroid.V_V => 28.4E-6 * 0.23, 
        mdl.eye.choroid.V_IntS => 28.4E-6 * 0.4, 
        mdl.eye.choroid.LF => LF_Choroid, 
        mdl.eye.choroid.sigma_IS => 0.2, 
        mdl.eye.choroid.k_out_ints => 0, 
        mdl.eye.choroid.k_in_ints => 0, 
        mdl.eye.choroid.k_PL_ex => k_PL_ints_clearance, 
        
        # added for eyes 
        mdl.RC_aq => 0.99, 
        mdl.RC_ret => 0.95, 
        mdl.RC_cho => 0.99,

        mdl.Q_AH => 0.000212,  # [L/hr]
        mdl.Q_PtA => 0.000234, # [L/hr]
        mdl.Q_BF => 0.00003,   # [L/hr]

        mdl.PS_ret => 2.64E-6, 
        mdl.PS_cho => 2.88E-6, 
        mdl.PS_cor => 1.19E-8, 
        mdl.PS_lens => 0, 

        mdl.V_AQ => 0.306E-3, # [L]
        mdl.V_VH => 1.41E-3,  # [L] 
        mdl.V_Lens => 0.41E-3,  # [L] # https://pmc.ncbi.nlm.nih.gov/articles/PMC2779573/
        mdl.V_Cor => 0.0887E-3, 

        # added for ADC PL
        mdl.DAR => 2, 
        mdl.frac_lys => 1, 
        mdl.k_diff => log(2)/(26.5/min_per_hr),  # DM4; assuming the retention half life was the diffusion half life; See Elahere script 
        mdl.CL_plasma_PL => 0,          # PL clearance rate in plasma 
        mdl.Kp_PL_endo_plasma => 1,     # partition coefficient, endothelial cells:plasma

        # assumed PL permeability surface area in the eye; computed = permeability * surface area 
        mdl.PS_cor_pl => 3.29E-8 * 1.19E-8/3.95E-8,    # [L/hr]; for DM4; see Elahere documentation; scaled based on mAb permeability 
        mdl.PS_ret_pl => 3.56E-7 * 2.64E-6/8.34E-6,    # [L/hr]; for DM4; see Elahere documentation; scaled based on mAb permeability 
        mdl.PS_lens_pl => 0,         # Abbvie rat study
        
        # assumed PL partition coefficients in the eye 
        mdl.Kp_ah_pl => 1, 
        mdl.Kp_vh_pl => 1, 

        # assumed ADC partition coefficient between aqueous humor and cornea 
        mdl.Kp_adc_cor_ah => 1, 

        # added nonlinear clearance of ADC in plasma 
        mdl.VMAX => 0, 
        mdl.KM => 1, 
    ])


    p_map_all = merge(p_map, p_map_global); 

    return p_map_all

end
