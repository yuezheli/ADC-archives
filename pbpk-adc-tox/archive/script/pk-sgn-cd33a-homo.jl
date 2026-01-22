# date: 6/30/2025 
# author: Yuezhe Li 
# purpose of this code: to fit for PK of vadastuximab talirine

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Parameters: @unpack
using DataFrames
using DataFramesMeta
using Plots
using CSV

include(@projectroot("Constants.jl"))
include(@projectroot("model/jones_homo.jl") )
include(@projectroot("model/param-pk.jl"))
include(@projectroot("model/init-pk.jl")) 
include(@projectroot("script/helper-infusion-dosing.jl"))
include(@projectroot("script/helper-visualization.jl"))
include(@projectroot("script/helper-process-outcome.jl"))

# observed PK; https://pubmed.ncbi.nlm.nih.gov/29196412/
pk_vt = Dict(
    0.05 => DataFrame(
        time_d = [0, 0.08, 0.25, 1], 
        ADC_ug_L = [364.91, 104, 67.74, 30.89]
    ), 
    0.04 => DataFrame(
        time_d = [0, 0.08, 0.25, 1], 
        ADC_ug_L = [231.92, 69.98, 43, 19.55], 
    ), 
    0.02 => DataFrame(
        time_d = [0, 0.08, 0.25, 1], 
        ADC_ug_L = [137.07, 40.08, 20.32, 7.17], 
    )
); 

for dose__ in collect(keys(pk_vt))
    @transform!(pk_vt[dose__], :ADC_uM = :ADC_ug_L/MW_IGG)
end

p_vt = deepcopy(p_base);
p_vt.k_deconj = 0                       # assumed to be none, since it is unknown
p_vt.DAR = 2                            # https://pubmed.ncbi.nlm.nih.gov/29119409/
p_vt.Kd = 1.2E-3                        # binding affinity between CD33 and vadastuximab, [uM], https://ashpublications.org/blood/article/122/8/1455/32075/SGN-CD33A-a-novel-CD33-targeting-antibody-drug
p_vt.k_out = 4.32                       # https://pubs.acs.org/doi/10.1021/acs.jmedchem.0c00691
p_vt.PS_Score = 30                      # tuned based on PK 
p_vt.kint_PS = 20                       # tuned based on PK 
p_vt.init_sR = 0.012                    # tuned based on PK 
p_vt.thalf_sR_adc = 10                  # tuned based on PK 

tspan = (-0.01, hr_per_day*84);      # [hr]
AddDose_q3w = [0., 21., 42., 63] * hr_per_day  # [hr]

sims_dose = [0.02, 0.04, 0.05]; # [mg/kg]

sol_pk_vt = Dict(); 

for init_dose in sims_dose
    tmp_sol_pk = InfusionDoses(init_dose, AddDose_q3w, p_vt, infusion_time = 0.1);
    push!(sol_pk_vt, init_dose => tmp_sol_pk);
end

plt_pk_vt = PlotSimulationPlasma(sol_pk_vt, pk_vt, adc_name = "vadastuximab talirine", colorPALETTE = :batlowKS, xrange = [0, 1.5], yrange = [1E-5, 1E-2], ylog = true)

savefig(plt_pk_vt, @projectroot("deliv/figure/pk/vadastuximab-talirine-homo.png"));

# PlotSimulationLiverEndoCyto(sol_pk_vt, adc_name = "vadastuximab talirine", colorPALETTE = :batlowKS, yrange = [1E-9, 1E-5], ylog = true)

# additional simulation based on dose used in https://pubmed.ncbi.nlm.nih.gov/29196412/
# in addition, 10 ug/kg was used Q4W for phase III trial; https://clinicaltrials.gov/study/NCT02785900#study-plan

sims_dose_phaseI = [5, 10, 20, 30, 40, 50, 60]*1E-3; # [mg/kg]
sol_pk_vt_q3w = Dict(); 
for init_dose in sims_dose_phaseI
    tmp_sol_pk = InfusionDoses(init_dose, AddDose_q3w, p_vt, infusion_time = 0.1);
    push!(sol_pk_vt_q3w, init_dose => tmp_sol_pk);
end

CSV.write(@projectroot("data/sim/vadastuximab-talirine-q3w.csv"), ProcessOutcome(sol_pk_vt_q3w))

AddDose_q4w = [0., 28., 56] * hr_per_day  # [hr]
sims_dose_phaseIII = 0.01; # [mg/kg]
sol_pk_vt_q4w = Dict(); 
for init_dose in sims_dose_phaseIII
    tmp_sol_pk = InfusionDoses(init_dose, AddDose_q4w, p_vt, infusion_time = 0.1);
    push!(sol_pk_vt_q4w, init_dose => tmp_sol_pk);
end

CSV.write(@projectroot("data/sim/vadastuximab-talirine-q4w.csv"), ProcessOutcome(sol_pk_vt_q4w))
