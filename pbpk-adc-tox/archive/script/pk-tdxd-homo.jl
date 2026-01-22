# date: 7/1/2025 
# author: Yuezhe Li 
# purpose of this code: to generate PK fit of T-Dxd

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

pk_doi =  CSV.read(@projectroot("data/doi-2017.csv"), DataFrame);  # https://pubmed.ncbi.nlm.nih.gov/29037983/

# convert observed data to dictionary 
pk_tdxd = Dict();
for dose_str in unique(pk_doi.Dose)
    tmp_pk_ = @rsubset(pk_doi, :Dose == dose_str);
    dose_num = parse(Float64, replace(dose_str, "mg/kg"=>""));
    tmp_df = DataFrame(
        time_d = tmp_pk_.time_day, 
        ADC_uM = tmp_pk_.T_DXd_ugperml*1E3/MW_IGG
    )
    push!(pk_tdxd, dose_num => tmp_df)
end

p_dxd = deepcopy(p_base);
p_dxd.DAR = 8.
p_dxd.k_deconj = 1.7E-7 * s_per_hr #  https://pubmed.ncbi.nlm.nih.gov/37787918/
p_dxd.k_out = 32.32                # https://pubmed.ncbi.nlm.nih.gov/37787918/
p_dxd.CL_PL_plasma = 15            # tuned for the clearance slope

tspan = (-0.01, hr_per_day*84);      # [hr]
AddDose_doi = [0., 22., 44.] * hr_per_day  # [hr]

sims_dose = [0.8, 1.6, 3.2, 5.4, 6.4, 8.0]; # [mg/kg]

sol_pk_tdxd = Dict(); 
for init_dose in sims_dose
    tmp_sol_pk = InfusionDoses(init_dose, AddDose_doi, p_dxd, infusion_time = 0.5);
    push!(sol_pk_tdxd, init_dose => tmp_sol_pk);
end

plt_pk_tdxd = PlotSimulationPlasma(sol_pk_tdxd, pk_tdxd, adc_name = "trastuzumab deruxtecan", colorPALETTE = :batlowKS, xrange = [0, 63], yrange = [1E-4, 100], ylog = true, legendcolumnsnum = 2);

savefig(plt_pk_tdxd, @projectroot("deliv/figure/pk/trastuzumab-deruxtecan-homo.png"));

# convert simulation to Q3W 
AddDose_q3w = [0., 21., 42., 63] * hr_per_day  # [hr]
sol_pk_tdxd_q3w = Dict(); 
for init_dose in sims_dose
    tmp_sol_pk = InfusionDoses(init_dose, AddDose_q3w, p_dxd, infusion_time = 0.5);
    push!(sol_pk_tdxd_q3w, init_dose => tmp_sol_pk);
end

CSV.write(@projectroot("data/sim/trastuzumab-deruxtecan-q3w.csv"), ProcessOutcome(sol_pk_tdxd_q3w))

# plot Dxd concentration at 6.4 mg/kg
# note that unit in the original figure was wrong and was corrected based on PK comparison of T-Dxd PK at 6.4 mg/kg 
obs_dxd_plasma_ugL = DataFrame(time_d = [0.058, 0.2, 0.893, 3, 8, 15, 21.97, 22, 43, 44, 50, 57, 64.2], conc = [3.186, 6.56, 4.854, 2.996, 1.542, 0.682, 0.361, 2.2, 0.5, 5.713, 2.24, 0.962, 0.51]); # https://pubmed.ncbi.nlm.nih.gov/29037983/

plasma_dxd = [sol_pk_tdxd[6.4].u[i].plasma_payload for i in 1:length(sol_pk_tdxd[6.4].t)]; # [uM]

plt_dxd = plot(xlabel = "Time (d)", ylabel = "Dxd (uM)", dpi = 300, size = (400,400), background_color_legend = nothing);
plot!(sol_pk_tdxd[6.4].t / hr_per_day, plasma_dxd, label = "sims, 6.4 mg/kg");
plot!(obs_dxd_plasma_ugL.time_d, obs_dxd_plasma_ugL.conc / MW_Dxd, label = "Doi et al., 2017", seriestype = :scatter); 
plot!(xlims = [0, 21], xticks = [0, 7, 14, 21], ylims = [1E-4, 1E2], yaxis = :log);

savefig(plt_dxd, @projectroot("deliv/figure/pk/trastuzumab-deruxtecan-homo-dxd.png"));
