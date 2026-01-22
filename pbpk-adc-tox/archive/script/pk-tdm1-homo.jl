# date: 7/1/2025 
# author: Yuezhe Li 
# purpose of this code: to generate PK fit for T-DM1

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

# read observed 
pk_yamamoto =  CSV.read("data/yamamoto-2015.csv",DataFrame); # https://academic.oup.com/jjco/article-abstract/45/1/12/887438?redirectedFrom=fulltext&login=false
pk_girish =  CSV.read("data/girish-2012.csv",DataFrame);     # https://link.springer.com/article/10.1007/s00280-011-1817-3
pk_tdm1_comb = vcat(pk_yamamoto, pk_girish);

# convert observed data to dictionary 
pk_tdm1 = Dict();
for dose_str in unique(pk_tdm1_comb.Dose)
    tmp_pk_ = @rsubset(pk_tdm1_comb, :Dose == dose_str);
    dose_num = parse(Float64, replace(dose_str, "mg/kg"=>""));
    tmp_df = DataFrame(
        time_d = tmp_pk_.time_day, 
        ADC_uM = tmp_pk_.T_DM1_ugperml*1E3/MW_IGG
    )
    push!(pk_tdm1, dose_num => tmp_df)
end

p_tmd1 = deepcopy(p_base);
p_tmd1.CL_PL_plasma = log(2)/3; # https://pubmed.ncbi.nlm.nih.gov/37787918/

tspan = (-0.01, hr_per_day*84);      # [hr]
AddDose_q3w = [0., 21., 42., 63] * hr_per_day  # [hr]

sims_dose = [0.3, 0.6, 1.2, 2.4, 3.6, 4.8, 1.8]; # [mg/kg]

sol_pk_tmd1 = Dict(); 

for init_dose in sims_dose
    tmp_sol_pk = InfusionDoses(init_dose, AddDose_q3w, p_tmd1, infusion_time = 0.1);
    push!(sol_pk_tmd1, init_dose => tmp_sol_pk);
end

plt_pk_tmd1 = PlotSimulationPlasma(sol_pk_tmd1, pk_tdm1, adc_name = "trastuzumab emtansine", colorPALETTE = :batlowKS, xrange = [0, 42], yrange = [1E-4, 10], ylog = true, legendcolumnsnum = 2)

savefig(plt_pk_tmd1, @projectroot("deliv/figure/pk/trastuzumab-emtansine-homo.png"));

CSV.write(@projectroot("data/sim/trastuzumab-emtansine-q3w.csv"), ProcessOutcome(sol_pk_tmd1))

# plot 3.6 mg/kg group DM1, trastuzumab, and T-DM1 plasma PK 
obs_dm1_plasma_ugL = DataFrame(time_d = [0, 0.005, 7, 14, 21], conc = [0.325, 6.071, 0.599, 0.418, 0.428]); # https://link.springer.com/article/10.1007/s00280-011-1817-3
plasma_dm1 = [sol_pk_tmd1[3.6].u[i].plasma_payload for i in 1:length(sol_pk_tmd1[3.6].t)]; # [uM]

plt_dm1 = plot(xlabel = "Time (d)", ylabel = "DM1 (uM)", dpi = 300, size = (400,400), background_color_legend = nothing);
plot!(sol_pk_tmd1[3.6].t / hr_per_day, plasma_dm1, label = "sims, 3.6 mg/kg");
plot!(obs_dm1_plasma_ugL.time_d, obs_dm1_plasma_ugL.conc / MW_Lys_MCC_DM1, label = "Girish et al., 2012", seriestype = :scatter); 
plot!(xlims = [0, 21], xticks = [0, 7, 14, 21], ylims = [1E-4, 1E-1], yaxis = :log);

savefig(plt_dm1, @projectroot("deliv/figure/pk/trastuzumab-emtansine-homo-dm1.png"));
