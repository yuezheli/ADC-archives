# date: 7/2/2025 
# author: Yuezhe Li 
# purpose of this code: to simulate different dosing schemes of Cantuzumab mertansine

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using DataFrames, CSV, DataFramesMeta
using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Parameters: @unpack
using Plots

# simulation 
include(@projectroot("Constants.jl"))
include(@projectroot("model/jones_homo.jl") )
include(@projectroot("model/param-pk.jl"))
include(@projectroot("model/init-pk.jl")) 
include(@projectroot("script/helper-infusion-dosing.jl"))
include(@projectroot("script/helper-process-outcome.jl"))

# DAR value not changed based on https://pubmed.ncbi.nlm.nih.gov/18301896/
p_cm = deepcopy(p_base);
p_cm.PS_Score = -2      # turn off using ACSIN score to fine tune PK 
p_cm.PS_kd = 0.0001     # tuned for huC242-DM1 PK 
p_cm.Kd = 1E-2          # https://pmc.ncbi.nlm.nih.gov/articles/PMC443115/
p_cm.init_sR = 5E-5     # Cheng et al., 2011; https://pmc.ncbi.nlm.nih.gov/articles/PMC4012263/; assuming soluble MUC1 have molecular weight of 250kDa
p_cm.k_deconj = 0.13    # fitted based on conjugated ADC 
p_cm.CL_PL_plasma = log(2)/3; # https://pubmed.ncbi.nlm.nih.gov/37787918/

# simulation based on Tolcher et al., 2003; https://pubmed.ncbi.nlm.nih.gov/12525512/
# MTD 235 mg/m2 Q3W (DLT is severe transaminase elevation)
tspan = (-0.01, hr_per_day*84);      # [hr]
AddDose_q3w_2003 = [0., 21., 42., 63] * hr_per_day  # [hr]
dose_q3w_2003 = [22, 44, 88, 132, 176, 235, 295]; # [mg/m2]

sol_pk_cm_2003 = Dict(); 

for init_dose in dose_q3w_2003
    tmp_sol_pk = InfusionDoses(init_dose*HT*HT/BW, AddDose_q3w_2003, p_cm, infusion_time = 0.5);
    push!(sol_pk_cm_2003, init_dose => tmp_sol_pk);
end

# additional dosing scheme from Heft et al., 2004
# https://aacrjournals.org/clincancerres/article/10/13/4363/94515/A-Phase-I-Study-of-Cantuzumab-Mertansine
# dosing QW; DLT is liver tox, at 115 mg/m2 QW

dose_qw_2004 = [40, 60, 80, 96, 115, 138]; # [mg/m2]
AddDose_qw_2004 = [0., 7., 14., 21., 28., 35., 42., 49., 56., 63., 70., 77.] * hr_per_day  # [hr]

sol_pk_cm_2004 = Dict(); 

for init_dose in dose_qw_2004
    tmp_sol_pk = InfusionDoses(init_dose*HT*HT/BW, AddDose_qw_2004, p_cm, infusion_time = 0.5);
    push!(sol_pk_cm_2004, init_dose => tmp_sol_pk);
end

# additional dose scheme from Rodin et al., 2008 
# https://link.springer.com/article/10.1007/s00280-007-0672-8
# dosing 3 times a week, 3 weeks per 4 weeks (3 weeks on, 1 week off)

dose_qw34_2008 = [30, 45, 60]; # [mg/m2]
AddDose_qw34_2008 = [0, 2, 4, 7, 9, 11, 14, 16, 18, 
                28, 30, 32, 35, 37, 39, 42, 44, 46, 
                56, 58, 60, 63, 65, 67, 70, 72, 74] * hr_per_day  # [hr]

sol_pk_cm_2008 = Dict(); 

for init_dose in dose_qw34_2008
    tmp_sol_pk = InfusionDoses(init_dose*HT*HT/BW, AddDose_qw34_2008, p_cm, infusion_time = 0.5);
    push!(sol_pk_cm_2008, init_dose => tmp_sol_pk);
end

# save simulation outcome 
CSV.write(@projectroot("data/sim/cantuzumab-mertansine-q3w.csv"), ProcessOutcome(sol_pk_cm_2003))
CSV.write(@projectroot("data/sim/cantuzumab-mertansine-qw.csv"), ProcessOutcome(sol_pk_cm_2004))
CSV.write(@projectroot("data/sim/cantuzumab-mertansine-tiw.csv"), ProcessOutcome(sol_pk_cm_2008))
