# author: Yuezhe LI
# date: 9/26/2023

using ComponentArrays

const DayToHour = 24.0 # convert day to hour
const MW_EDG = 15.0e4; # mAb molecular weight [Da]

# this is for homo (with tumor as the 16th organ)
p_N_16 = ComponentArray(
    PS_Score = 6.,                          # AC-SINS score, [unitless], https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5825195/
    PS_kd = -1.0,                           # placeholder params; disabled; 
    KD6_WT = 700.0,                         # binding affinity between ADC and FcRn, [nM], https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12461
    Kd = 0.314e-3,                          # ADC-receptor binding affinity, [uM],
    Kon = 2.07e3,                           # on-rate for receptor-ADC binding, [uM-1.hr-1], https://www.bio-rad.com/webroot/web/pdf/lsr/literature/Bulletin_6044A.pdf
    # tumor perfusion
    tumorPLQ = 12.7,                        # tumor perfusion, [L/h/L]
    # soluble HER2 params
    thalf_sR = 5.,    # half life of soluble HER2 in plasma [h], https://europepmc.org/article/ppr/ppr584999
    thalf_sR_adc = 120., # half life of sHER2:ADC complex in plasma [h], https://europepmc.org/article/ppr/ppr584999
    init_sR = 0.02,     # calculated from Moreno-Aspitia et al., 2013, [umol], with assumption that sHER2 molecular weight = 1E5 g/mol, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4126833/
    # Related to ADC infusion 
    infusion = 0., 
);

# cyno data
p_N_15 = ComponentArray(
    PS_Score = 6.,                          # AC-SINS score, [unitless], https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5825195/
    PS_kd = -1.0,                           # placeholder params; disabled; 
    KD6_WT = 700.0,                         # binding affinity between ADC and FcRn, [nM], https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12461
    endothelial_scaling_factor = 600.,                       
    # Related to ADC infusion 
    infusion = 0., 
);

