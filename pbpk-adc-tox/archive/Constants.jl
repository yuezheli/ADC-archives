# date: 6/25/2025 
# author: Yuezhe Li 
# purpose of this code: to host constance that will be used in modeling & simulation 

const N_av = 6.0221409e+23  # Avogadro constant 
const MW_IGG = 15E4;        # molecular weight for mAb/ ADC [g/mol]
const hr_per_day = 24.0     # conversion number between hour and day
const s_per_hr = 3600.      # conversion between second and hour 
const BW = 75.0;            # body weight for a standard human, [kg]
const HT = 1.85             # height for a standard human, [m]
const V_Plasma = 3.126;     # plasma volume for a standard human, [L]

const MW_MMAE = 718         # molecular weight of MMAE, [Da]; https://pubchem.ncbi.nlm.nih.gov/compound/Monomethyl-Auristatin-E
const MW_DM1 = 738          # molecular weight of DM1, [Da]; https://pubchem.ncbi.nlm.nih.gov/compound/dm-1
const MW_Lys_MCC_DM1 = 1100 # molecular weight of T-DM1 payload metabolite; https://pubmed.ncbi.nlm.nih.gov/37787918/
const MW_Dxd = 493          # molecular weight of T-DM1; https://adc.bocsci.com/product/dxd-cas-1599440-33-1-333861.html
