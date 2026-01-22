# author: Yuezhe Li
# date: 6/8/23
# purpose of this script: to host parameters used for simulation 

# base parameters are T-DM1 parameters
p_base = ComponentArray(
    # ADC-related parameters
    Kon = 2.07e3,                           # on-rate for receptor-ADC binding, [uM-1.hr-1], https://www.bio-rad.com/webroot/web/pdf/lsr/literature/Bulletin_6044A.pdf  
    Prob_deg = 0.95,                        # probability of free ADC being trafficked into degredation; https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12461
    PS_Score = 6.,                          # AC-SINS score, [unitless], https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5825195/
    PS_kd = -1.0,                           # placeholder params; disabled; 
    KD6_WT = 700.0,                         # binding affinity between ADC and FcRn, [nM], https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12461
    kint_PS = 0.0380,                       # [1/h]. Internalization rate for polyspecificity-membrane bound mAb
    DAR = 3.5,                              # average drug : antibody ratio, [unitless], https://www.nature.com/articles/s41416-019-0635-y
    frac_lys = 1,                           # fraction of payload that survives endothelial lysosomes to be released to the cytosol
    k_PL = 0,                               # free payload degredation rate, intracellular, [1/hr], 
    CL_PL_plasma = 0,                       # free payload clearance rate, plasma, [1/hr], 
    k_PL_ex = 0.34,                         # free payload degredation rate, tissue intersititium, [1/hr], 
    k_out = 0.14,                           # rate of payload leaving of the cell through diffusion, [1/hr], https://europepmc.org/article/ppr/ppr584999
    Kd = 0.314e-3,                          # ADC-receptor binding affinity, [uM], 
    k_deconj = 8.5E-7 * s_per_hr,           # deconjugation rate between mAb and payload, [1/hr], https://europepmc.org/article/ppr/ppr584999
    # tumor penetration-related params (all dummy for now)
    P_ADC = 334 / hr_per_day,      # rate of permeability [um/h]; https://pubmed.ncbi.nlm.nih.gov/27029797/
    D_ADC = 0.022 / hr_per_day,    # rate of diffusion [cm^2/h]; https://pubmed.ncbi.nlm.nih.gov/27029797/
    Rcap = 8.,                     # radium of tumor blood capillary [um]; Shah et al., 2012 # https://pubmed.ncbi.nlm.nih.gov/23151991/
    Rkrogh = 75. ,                 # radius of average distance between 2 blood vessels [um]; Shah et al., 2012 # https://pubmed.ncbi.nlm.nih.gov/23151991/
    epsilon = 0.24,                # partition coeffficient of ADC between plasma:tumor; Shah et al., 2012 # https://pubmed.ncbi.nlm.nih.gov/23151991/
    TumorVolume_L = 1E-3,          # dummy parameter for tumor volume at time 0
    lesion_num = 1,                # dummy parameters for number of lesions 
    # soluble receptor (sHER2) params
    thalf_sR = 5.,          # half life of soluble HER2 in plasma [h], https://europepmc.org/article/ppr/ppr584999
    thalf_sR_adc = 120.,    # half life of sHER2:ADC complex in plasma [h], https://europepmc.org/article/ppr/ppr584999
    init_sR = 0.004,        # calculated from Moreno-Aspitia et al., 2013, [umol], with assumption that sHER2 molecular weight = 1E5 g/mol, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4126833/
    # ADC infusion rate  
    infusion = 0., 
);

