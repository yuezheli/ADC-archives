# date: 10/23/24
# author: Yuezhe Li 
# purpose of this code: to hold T-DM1 in vitro parameter

p_tdm1 = ComponentArray(
    tdouble_pos = 37.4,                     # doubling time, [hr], https://www.cellosaurus.org/CVCL_0033 (SK-BR-3)
    tdouble_neg = 61.,                      # doubling time, [hr], https://www.cellosaurus.org/CVCL_1566 (NCI-H520)
    tau = 1.,                               # time delay in tumor killing, [hr], 
    DAR = 3.5,                              
    ic50_pl = 23.8E-3,                      # payload IC50, [uM], https://link.springer.com/article/10.1007/s10928-023-09884-6
    k_kill_max = 0.014,                     # maximum in vitro killing rate, [1/hr], between 0.02 and 0.04 based on https://link.springer.com/article/10.1007/s10928-023-09884-6
    k_PL = 0.,                              # free payload degredation rate, [1/hr], 
    k_PL_ex = 0.,                           # free payload degredation rate, [1/hr], 
    Rcopies = 1.6E6,                        # surface receptor copy number, [molecules], https://pubmed.ncbi.nlm.nih.gov/26766593/
    k_deg = 0.4572,                         # HER2 degredation rate, [h-1], https://link.springer.com/article/10.1007/s10928-023-09884-6
    k_rec = 0.0864,                         # HER2 recycling rate, [h-1], https://link.springer.com/article/10.1007/s10928-023-09884-6
    k_endo = 0.15,                          # HER2:ADC endocytosis rate, [h-1], https://link.springer.com/article/10.1007/s10928-023-09884-6
    k_on = 2.07E3,                          # on-rate for receptor-ADC binding, [uM-1.hr-1], https://www.bio-rad.com/webroot/web/pdf/lsr/literature/Bulletin_6044A.pdf
    Kd = 0.314E-3,                          # ADC-receptor binding affinity, [uM], https://www.nature.com/articles/s41467-023-44533-z
    k_out = 0.14,                           # rate of payload leaving of the cell through diffusion, [1/hr], 
    k_in = 0.21,                            # rate of payload entering cell through diffusion, [1/hr], 
    k_lys = 2.4,                            # rate for ADC to be graded in cells to release payload, [1/hr], computed from endosome/ lysosome half life https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2112524/pdf/jc972508.pdf
    V_medium = 1E-4,                        # medium volume, https://www.nature.com/articles/s41467-023-44533-z
    Nmax = 1E7,                             # maximum cell per well, https://www.thermofisher.com/us/en/home/references/gibco-cell-culture-basics/cell-culture-protocols/cell-culture-useful-numbers.html?gclid=Cj0KCQjwveK4BhD4ARIsAKy6pMLSiIPwT_HPvpeO6_Ra_oDdqjl1jIZFrRCgDHr3ovuEo5I3r_FUMX8aAvKXEALw_wcB&ef_id=Cj0KCQjwveK4BhD4ARIsAKy6pMLSiIPwT_HPvpeO6_Ra_oDdqjl1jIZFrRCgDHr3ovuEo5I3r_FUMX8aAvKXEALw_wcB:G:s&s_kwcid=AL!3652!3!530416915609!!!g!!!382790548!125487008778&cid=bid_clb_cce_r01_co_cp0000_pjt0000_bid00000_0se_gaw_dy_pur_con&s_kwcid=AL!3652!3!530416915609!!!g!!&gad_source=1
);