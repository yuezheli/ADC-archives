# author: Yuezhe Li
# date: 6/8/23
# purpose of this script: to host parameters used for simulation 

# define constants
const DayToHour = 24.0 # convert day to hour
const N_av = 6.0221409e+23 # Avogadro constant

const MW_EDG = 15.0e4; # mAb molecular weight [Da]
const BW = 75.0; # [kg]
const V_Plasma = 3.126; # [L]

using ComponentArrays

p_base = ComponentArray(
    # parameters that are either assumed, or likely to be constants
    Kon = 2.07e3,                           # on-rate for receptor-ADC binding, [uM-1.hr-1], https://www.bio-rad.com/webroot/web/pdf/lsr/literature/Bulletin_6044A.pdf
    Emax_Payload = 0.139,                   # maximum in vitro killing rate, [1/hr]
    tumor_cell_radius = 8,                  # radius of tumor cell, [um]
    k_out = 0.14,                           # rate of payload leaving of the cell through diffusion, [1/hr], https://europepmc.org/article/ppr/ppr584999
    k_in = 0.21,                            # rate of payload entering cell through diffusion, [1/hr], https://europepmc.org/article/ppr/ppr584999
    k_lys = 2.4,                            # rate for ADC to be graded in cells to release payload, [1/hr], computed from endosome/ lysosome half life https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2112524/pdf/jc972508.pdf
    Kkill_scale = 0.1,                      # scaling parameters between in vitro and in vivo anti-tumor effect, [unitless]
    # ADC-related parameters
    PS_Score = 6.,                          # AC-SINS score, [unitless], https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5825195/
    PS_kd = -1.0,                           # placeholder params; disabled; 
    KD6_WT = 700.0,                         # binding affinity between ADC and FcRn, [nM], https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12461
    DAR = 3.5,                              # average drug : antibody ratio, [unitless], https://www.nature.com/articles/s41416-019-0635-y
    k_PL = 0.34,                            # free payload degredation rate, [1/hr], 
    IC50_Payload = 23.8E-3,                 # payload IC50, [uM], https://europepmc.org/article/ppr/ppr584999
    tau = 1.,                               # time delay in tumor killing, [hr], 
    # tumor-related parameters
    tdouble = 40*DayToHour,                 # tummor doubling time, [h]
    # receptor-related parameters (from Scheuher et al., 2022, https://europepmc.org/article/ppr/ppr584999)
    Rcopies = 1.6E6,                        # surface receptor copy number, [molecules], https://pubmed.ncbi.nlm.nih.gov/26766593/
    Kd = 0.314e-3,                          # ADC-receptor binding affinity, [uM], 
    k_deg = 0.4572,                         # HER2 degredation rate, [h-1]
    k_rec = 0.0864,                         # HER2 recycling rate, [h-1]
    k_endo = 0.15,                          # HER2:ADC endocytosis rate, [h-1]
    # tumor penetration-related params
    P_ADC = 50.34 / 24,    # rate of permeability [um/h]; https://europepmc.org/article/ppr/ppr584999
    D_ADC = 1.74E-3 / 24,  # rate of diffusion [cm^2/h]; https://europepmc.org/article/ppr/ppr584999
    Rcap = 8.,             # radium of tumor blood capillary [um]; Shah et al., 2012 # https://pubmed.ncbi.nlm.nih.gov/23151991/
    Rkrogh = 75. ,         # radius of average distance between 2 blood vessels [um]; Shah et al., 2012 # https://pubmed.ncbi.nlm.nih.gov/23151991/
    epsilon = 0.24,        # partition coeffficient of ADC between plasma:tumor; Shah et al., 2012 # https://pubmed.ncbi.nlm.nih.gov/23151991/
    # soluble HER2 params
    thalf_sR = 5.,    # half life of soluble HER2 in plasma [h], https://europepmc.org/article/ppr/ppr584999
    thalf_sR_adc = 120., # half life of sHER2:ADC complex in plasma [h], https://europepmc.org/article/ppr/ppr584999
    init_sR = 0.02,     # calculated from Moreno-Aspitia et al., 2013, [umol], with assumption that sHER2 molecular weight = 1E5 g/mol, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4126833/
    # Related to ADC infusion 
    infusion = 0., 
);


# calculate relevent constants
const Vc = 4/3 * pi * (p_base.tumor_cell_radius*1E-6)^3 * 1E3 ; # [L]

# things to change for T-DM1 and T-Dxd
# k_in, k_out, DAR, k_PL
# for T-DM1, the linker is not cleavable, thus the k_in & k_out should be 0
# check for T-Dxd
# DAR for T-DM1 = 3.5, DAR for T-Dxd = 8
# k_in_PL = 5.95E-5 s-1 (T-DM1), 0.0128 s-1 (T-Dxd)
# k_out_PL = 3.94E-5 s-1 (T-DM1), 0.00898 s-1 (T-Dxd)
# tau could change across cell line, between 1E-3 day to 3.84 day
# comparison between T-DM-1 and brentuximab-vedotin:  
# P_ADC = 334. / 24,   # rate of permeability [um/h]; Shah et al., 2012 # https://pubmed.ncbi.nlm.nih.gov/23151991/
# D_ADC = 0.022 / 24,  # rate of diffusion [cm^2/h]; Shah et al., 2012 # https://pubmed.ncbi.nlm.nih.gov/23151991/