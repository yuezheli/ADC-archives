# date: 12/29/2025 
# author: Yuezhe Li 
# purpose of this code: to host constants 

# Scaling factors
const s_per_hr = 3600.0 # s/hr
const hr_per_day = 24.0 # hr/day
const min_per_hr = 60.0 # hr/day
const nmol_per_mol = 1E9  # nmol/mol
const mg_per_g = 1E3 # mg/g
const ug_per_g = 1E6 # ug/g
const mL_per_L = 1E3 # ml/L
const um2_per_mm2 = 1E6 # um^2/mm^2
const mm3_per_L = 1E6 # mm^3/L
const N_AV = 6.022E23 # Avogadro's number [#/mol]
const uL_per_L = 1E6 # uL/L
const ng_per_ug = 1E3 # ng/ug 
const um_per_mm = 1E3 # um/mm
const g_cell_per_L = 1E3 # [g/L] density of tumor tissue
const ng_per_g = 1E9  # ng/g
const mm3_per_mL = 1E3 # mm^3/mL

const mouse_WT = 20E-3  # assumes a 20 g mouse (g to kg)
const cyno_WT = 2.6     # assumes a 2.6 kg cyno monkey (kg)
const human_WT = 70.0   # assumes a 70 kg human (kg)
const human_BSA = 1.9

const cell_radius = 5E-4  # [cm] # assuming the cell radius = 5 um

const MW_IGG = 15E4
const MW_TDXD = 15.3E4      # [g/mol] typical ADC molecular weight
const MW_TDM1 = 14.8E4
const MW_trastuzumab = 14.55E4 

const MW_MMAE = 718         # molecular weight of MMAE, [Da]; https://pubchem.ncbi.nlm.nih.gov/compound/Monomethyl-Auristatin-E
const MW_DM1 = 738          # molecular weight of DM1, [Da]; https://pubchem.ncbi.nlm.nih.gov/compound/dm-1
const MW_Lys_MCC_DM1 = 1100 # molecular weight of DM1 payload metabolite; https://pubmed.ncbi.nlm.nih.gov/37787918/
const MW_Dxd = 493          # molecular weight of Dxd; https://adc.bocsci.com/product/dxd-cas-1599440-33-1-333861.html
const MW_PBD = 500          # molecular weight PBD dimer

# Define MRG colors 
MRG_green = RGB(25 / 255, 160 / 255, 146 / 255);
MRG_blue = RGB(41 / 255, 127 / 255, 157 / 255);
MRG_light = RGB(182 /255, 217 / 255, 213 / 255);
MRG_yellow = RGB(255 / 255, 195 / 255, 107 / 255);
MRG_red = RGB(255 / 255, 123 / 255, 118 / 255);

