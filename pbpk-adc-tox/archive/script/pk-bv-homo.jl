# date: 6/27/2025 
# author: Yuezhe Li 
# purpose of this code: to simulate human PK of brentuximab vedotin, dose at 1.2mg/kg, 1.8mg/kg, 2.7mg/kg.
# https://www.nejm.org/doi/full/10.1056/NEJMoa1002965

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using DataFrames
using DataFramesMeta
using Plots
using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Parameters: @unpack
using CSV

include(@projectroot("Constants.jl"))
include(@projectroot("model/jones_homo.jl") )
include(@projectroot("model/param-pk.jl"))
include(@projectroot("model/init-pk.jl")) 
include(@projectroot("script/helper-infusion-dosing.jl"))
include(@projectroot("script/helper-visualization.jl"))
include(@projectroot("script/helper-process-outcome.jl"))

# observed ADC PK; Younes et al., 2010; https://www.nejm.org/doi/full/10.1056/NEJMoa1002965
pk_bv = Dict(
    1.2 => DataFrame(
        time_d = [0.13, 0.26, 0.26, 1, 3, 7, 14, 20.95, 21, 21.2, 21.26, 22.04, 24, 28, 35, 42], 
        ADC_ug_mL = [20.23, 17.55, 14.18, 6.26, 2.16, 1.27, 0.45, 0.18, 24.16, 18.84, 15.22, 7.75, 3.31, 1.36, 0.58, 0.13]
    ), 
    1.8 => DataFrame(
        time_d = [0.14, 0.27, 1.04, 3., 7, 14, 21.05, 21.14, 21.33, 22.04, 24, 28, 35, 42], 
        ADC_ug_mL = [31, 26.87, 12.3, 6.26, 2.58, 1.14, 0.39, 31, 25.03, 10, 5.63, 2.49, 0.95, 0.49], 
    ), 
    2.7 => DataFrame(
        time_d = [0.2, 0.33, 1.05, 3, 7, 14, 21.05, 21.14, 21.27, 22.05, 24, 28, 35, 42], 
        ADC_ug_mL = [44.18, 35.7, 20.23, 8.62, 4.09, 2.08, 0.72, 50.93, 39.72, 20.96, 9.94, 3.95, 2.32, 1.18], 
    )
); 

# observed plasma MMAE; Younes et al., 2010; https://www.nejm.org/doi/full/10.1056/NEJMoa1002965
pk_mmae = Dict(
    1.2 => DataFrame(
        time_d = [0.25, 0.33, 1, 3, 7, 14, 20.9, 21, 21.2, 22, 24, 28, 35, 42], 
        ADC_ugL = [0.341, 2.316, 5.51, 5.51, 1.673, 0.256, 0.134, 0.566, 0.756, 2.489, 2.774, 1.448, 0.354, 0.07]
    ), 
    1.8 => DataFrame(
        time_d = [0.13, 1, 2, 4, 7, 10, 14, 20.9, 21, 22, 24, 28, 35, 42], 
        ADC_ugL = [0.756, 4.769, 4.437, 3.446, 2.154, 1.085, 0.546, 0.138, 1.448, 3.092, 2.982, 1.798, 0.424, 0.192], 
    ), 
    2.7 => DataFrame(
        time_d = [0.13, 1, 2, 3, 4, 7, 10, 14, 21, 21.25, 22, 24, 28, 35, 42], 
        ADC_ugL = [1.009, 6.601, 7.356, 6.844, 6.601, 4.279, 2.078, 1.448, 0.116, 1.501, 3.092, 4.127, 2.876, 0.812, 0.149], 
    )
); 

for dose__ in collect(keys(pk_bv))
    @transform!(pk_bv[dose__], :ADC_uM = :ADC_ug_mL*1E3/MW_IGG)
end

# simulation
p_bv = deepcopy(p_base);
p_bv.PS_Score = 15          # turned based on PK
p_bv.init_sR = 18.8/88E3    # assuming MW of sCD30 to be 88kDa [https://www.nature.com/articles/bcj201785], concentration to be 18.8 ng/mL [https://pmc.ncbi.nlm.nih.gov/articles/PMC2756127/]
p_bv.thalf_sR = 0.1
p_bv.thalf_sR_adc = 28      # https://pubmed.ncbi.nlm.nih.gov/10881697/
p_bv.k_deconj = 0.0785/hr_per_day   # https://pmc.ncbi.nlm.nih.gov/articles/PMC5574006/         
p_bv.Kd = 0.7E3/MW_IGG      # binding affinity between CD30 and brentuximab, [uM], https://www.sciencedirect.com/science/article/pii/S0006497120443186
p_bv.k_out = 2.             # https://bmccancer.biomedcentral.com/articles/10.1186/s12885-019-5336-7
p_bv.DAR = 4.4              # https://pubmed.ncbi.nlm.nih.gov/23151991/
p_bv.CL_PL_plasma = 60/42   # https://pmc.ncbi.nlm.nih.gov/articles/PMC8004929

tspan = (-0.01, hr_per_day*84);      # [hr]
AddDose_q3w = [0., 21., 42., 63] * hr_per_day  # [hr]

sims_dose = [1.2, 1.8, 2.7]; # [mg/kg]

sol_pk_bv = Dict(); 
for init_dose in sims_dose
    tmp_sol_pk = InfusionDoses(init_dose, AddDose_q3w, p_bv, infusion_time = 0.5);
    push!(sol_pk_bv, init_dose => tmp_sol_pk);
end

## ADC PK visualization
plt_pk_bv = PlotSimulationPlasma(sol_pk_bv, pk_bv, adc_name = "brentuximab vedotin", colorPALETTE = :batlowKS, xrange = [0, 42], yrange = [1E-3, 1], ylog = true);

savefig(plt_pk_bv, @projectroot("deliv/figure/pk/brentuximab-vedotin-homo.png"));

# plasma MMAE PK
plasma_mmae_1point2 = [sol_pk_bv[1.2].u[i].plasma_payload for i in 1:length(sol_pk_bv[1.2].t)]; # [uM]
plasma_mmae_1point8 = [sol_pk_bv[1.8].u[i].plasma_payload for i in 1:length(sol_pk_bv[1.8].t)]; # [uM]
plasma_mmae_2point7 = [sol_pk_bv[2.7].u[i].plasma_payload for i in 1:length(sol_pk_bv[2.7].t)]; # [uM]

plt_mmae = plot(xlabel = "Time (d)", ylabel = "MMAE (uM)", dpi = 300, size = (400,400), background_color_legend = nothing);
plot!(sol_pk_bv[1.2].t/hr_per_day , plasma_mmae_1point2, label = "sims, 1.2 mg/kg Q3W", color = :blue);
plot!(pk_mmae[1.2].time_d , pk_mmae[1.2].ADC_ugL / MW_MMAE, label = "obs, 1.2 mg/kg Q3W", color = :blue, seriestype = :scatter);
plot!(sol_pk_bv[1.8].t/hr_per_day , plasma_mmae_1point8, label = "sims, 1.8 mg/kg Q3W", color = :red);
plot!(pk_mmae[1.8].time_d , pk_mmae[1.8].ADC_ugL / MW_MMAE, label = "obs, 1.8 mg/kg Q3W", color = :red, seriestype = :scatter);
plot!(sol_pk_bv[2.7].t/hr_per_day , plasma_mmae_2point7, label = "sims, 2.7 mg/kg Q3W", color = :green);
plot!(pk_mmae[2.7].time_d , pk_mmae[2.7].ADC_ugL / MW_MMAE, label = "obs, 2.7 mg/kg Q3W", color = :green, seriestype = :scatter);
plot!(xlims = [0, 42], xticks = [0, 7, 14, 21, 28, 35, 42], ylims = [1E-4, 5E-2], yaxis = :log);

savefig(plt_mmae, @projectroot("deliv/figure/pk/brentuximab-vedotin-homo-mmae.png"));

# save simulation outcome for further analysis 
CSV.write(@projectroot("data/sim/brentuximab-vedotin-q3w.csv"), ProcessOutcome(sol_pk_bv))
