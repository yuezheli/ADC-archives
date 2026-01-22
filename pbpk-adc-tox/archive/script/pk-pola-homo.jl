# date: 6/30/2025
# author: Yuezhe Li 
# purpose of this code: to fit for PK of polatuzumab vedotin 

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

p_pola = deepcopy(p_base);
p_pola.DAR = 3.5                        # https://www.accessdata.fda.gov/drugsatfda_docs/nda/2019/761121Orig1s000ChemR.pdf
p_pola.k_out = 2.                       # https://bmccancer.biomedcentral.com/articles/10.1186/s12885-019-5336-7
p_pola.Kd = 0.26E-3                     # binding affinity between CD79b and polatuzumab, [uM], https://ashpublications.org/blood/article/114/13/2721/26429/Therapeutic-potential-of-an-anti-CD79b-antibody
p_pola.init_sR = 13.9E3/N_av*157*5E6*1E6   # calculated based on the assumption that there are 157 B cell/ uL blood [https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.b.20547] and, 5 L blood, and 13.9E3 CD79b per B cell [https://pubmed.ncbi.nlm.nih.gov/16531332/]
p_pola.PS_Score = 9              # turned based on PK
p_pola.thalf_sR_adc = 4          # assumed; half life of B cell after polatuzumab vedotin is 4 hr; https://pubmed.ncbi.nlm.nih.gov/17374736/
p_pola.thalf_sR = 24             # assumed
p_pola.k_deconj = 0.006473       # https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12137
p_pola.CL_PL_plasma = 60/42      # https://pmc.ncbi.nlm.nih.gov/articles/PMC8004929

# observed PK; https://pubmed.ncbi.nlm.nih.gov/33029633/
pk_pola = Dict(
    1.8 => DataFrame(
        time_d = [0.08, 0.28, 1, 3, 4, 7, 10, 14, 21], 
        acMMAE_ug_L = [654.68, 468.72, 350.45, 175.08, 152.88, 77.94, 40.63, 27.6, 15.94]
    ), 
    1.0 => DataFrame(
        time_d = [0.07, 0.43, 1, 3, 4, 7, 10, 14, 21], 
        acMMAE_ug_L = [307.17, 230.18, 184.69, 79.82, 62.72, 32.53, 21.57, 9.61, 3.97], 
    )
); 

pk_mmae = Dict(
    1.0 => DataFrame(
        time_d = [0.03, 0.18, 1, 3, 4, 7, 10, 14, 21], 
        MMAE_ug_L = [0.23, 0.44, 0.89, 1.29, 1.25, 0.75, 0.5, 0.28, 0.06]
    ), 
    1.8 => DataFrame(
        time_d = [0.1, 0.25, 1, 3, 4, 7, 10, 14, 21], 
        MMAE_ug_L = [0.28, 0.45, 1.22, 1.73, 1.41, 1.42, 0.65, 0.37, 0.15]
    )
);

for dose__ in collect(keys(pk_pola))
    @transform!(pk_pola[dose__], :ADC_uM = :acMMAE_ug_L/MW_MMAE/p_pola.DAR); # ADC conc = acMMAE/DAR, [uM]
end

tspan = (-0.01, hr_per_day*84);      # [hr]
AddDose_q3w = [0., 21., 42., 63] * hr_per_day  # [hr]

sims_dose = [1.0, 1.8]; # [mg/kg]

sol_pk_pola = Dict(); 
for init_dose in sims_dose
    tmp_sol_pk = InfusionDoses(init_dose, AddDose_q3w, p_pola, infusion_time = 0.5);
    push!(sol_pk_pola, init_dose => tmp_sol_pk);
end

# PK visualization, mAb + ADC
plt_pk_pola = PlotSimulationPlasma(sol_pk_pola, pk_pola, adc_name = "polatuzumab vedotin", colorPALETTE = :batlowKS, xrange = [0, 21], yrange = [1E-3, 1], ylog = true)

savefig(plt_pk_pola, @projectroot("deliv/figure/pk/polatuzumab-vedotin-homo.png"));

# PK visualization, MMAE
plasma_mmae_1 = [sol_pk_pola[1.0].u[i].plasma_payload for i in 1:length(sol_pk_pola[1.0].t)]; # [uM]
plasma_mmae_1point8 = [sol_pk_pola[1.8].u[i].plasma_payload for i in 1:length(sol_pk_pola[1.8].t)]; # [uM]

plt_mmae = plot(xlabel = "Time (d)", ylabel = "MMAE (uM)", dpi = 300, size = (400,400), background_color_legend = nothing);
plot!(sol_pk_pola[1.0].t/hr_per_day , plasma_mmae_1, label = "sims, 1.0 mg/kg Q3W", color = :blue);
plot!(pk_mmae[1].time_d , pk_mmae[1].MMAE_ug_L / MW_MMAE, label = "obs, 1.0 mg/kg Q3W", color = :blue, seriestype = :scatter);
plot!(sol_pk_pola[1.8].t/hr_per_day , plasma_mmae_1point8, label = "sims, 1.8 mg/kg Q3W", color = :red);
plot!(pk_mmae[1.8].time_d , pk_mmae[1.8].MMAE_ug_L / MW_MMAE, label = "obs, 1.8 mg/kg Q3W", color = :red, seriestype = :scatter);
plot!(xlims = [0, 42], xticks = [0, 7, 14, 21, 28, 35, 42], ylims = [1E-4, 5E-2], yaxis = :log);

savefig(plt_mmae, @projectroot("deliv/figure/pk/polatuzumab-vedotin-homo-mmae.png"));

# save outcome for future analysis
CSV.write(@projectroot("data/sim/polatuzumab-vedotin-q3w.csv"), ProcessOutcome(sol_pk_pola))
