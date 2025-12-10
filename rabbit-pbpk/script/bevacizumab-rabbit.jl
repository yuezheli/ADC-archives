# date: 10/10/2025 
# author: Yuezhe Li 
# purpose of this code: to test bevacizumab (Avastin) in rabbit PK 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using DifferentialEquations, ModelingToolkit
using Plots
using Plots.Measures
using DataFrames
using DataFramesMeta
using SymbolicIndexingInterface  # for parameter_index()
using CSV
using Sundials  # for CVODE_BDF()
using Statistics  # for mean()

#####################
# observed PK; digitized from Bakri et al., Ophthalmology, 2007 
# https://pubmed.ncbi.nlm.nih.gov/17467524/

obs_vh = DataFrame(
    time_d = [1, 3, 8, 15, 30], 
    conc_ugmL = [421.95, 343.97, 131.31, 45.29, 4.76]
); 

obs_ah = DataFrame(
    time_d = [1, 3, 8, 15, 30], 
    conc_ugmL = [24.72, 38.91, 11.61, 4.82, 0.62]
);

obs_serum = DataFrame(
    time_d = [1, 3, 8, 15, 30], 
    conc_ugmL = [0.61, 2.43, 3.24, 1.8, 0.4]
);

#####################
# helper functions
include(@projectroot("script/helper-funcs.jl"));
include(@projectroot("script/helper-initiation.jl")); 
include(@projectroot("script/helper-parameters-rabbit.jl"));
include(@projectroot("script/helper-infusion.jl"));

# post-proressing functions
include(@projectroot("script/helper-compute-tissue-adc.jl"));
include(@projectroot("script/helper-compute-tissue-free-pl.jl"));
include(@projectroot("script/helper-cmax-cavg-adc.jl"));
include(@projectroot("script/helper-cmax-cavg-payload.jl"));

# model functions
include(@projectroot("model/mab-pbpk-modular-eye.jl"));

# constant function
include(@projectroot("model/constants.jl"));

#####################

## create pbpk model
@time pbpk = create_pbpk();
@time pbpk_rabbit = mtkcompile(pbpk);

## Define the simulation timespan
tspan = (0.0,504);  # hr

## parameters, based on fitting from human PK
param_avastin = create_base_rabbit_pbpk_param(PS_Score = 5, mdl = pbpk_rabbit, Scale_Factor_Rabbit = 40); 
param_avastin[pbpk_rabbit.kdeg_FcRn_Ab] = 32.2 # Bussing and Shah, 2020; 

## IVT dosing (injection into vitreous humor at time 0), 1.25 mg, as in Bakri et al., 2007; https://pubmed.ncbi.nlm.nih.gov/17467524/
DOSE_IVT = 1.25  # [mg]
V_VH_rabbit = param_avastin[pbpk_rabbit.V_VH]
u0_ivt = pbpk_initial_condition(0, pbpk_rabbit, FcRn_Conc = 49.8); # FcRn_Conc concentration updated based on Bussing ans Shah, 2020; https://pubmed.ncbi.nlm.nih.gov/32876799/
u0_ivt[pbpk_rabbit.eye.vh_igg_exg.C_VH.val] = DOSE_IVT*1E-3/V_VH_rabbit/MW_IGG*1e6; # [uM]

## simulation 
@time prob_mtk_ivt = ODEProblem(pbpk_rabbit, merge(Dict(u0_ivt), param_avastin), tspan);  
@time sol_mtk_ivt = solve(prob_mtk_ivt, alg=CVODE_BDF());

#####################

## plasma mAb visualization 
plt__mAb = plot(ylabel = "mAb Concentration (ug/mL)", legend = :bottom, dpi = 300, size = (400, 400)); 
plot!(ylims = [0.01, 1E3],  yaxis = :log10); 
plot!(sol_mtk_ivt, idxs = pbpk_rabbit.eye.ah_igg_exg.C_AQ * MW_IGG/1E3, label = "AH (pbpk)", lw = 2, alpha = 0.9, color = :blue);
plot!(obs_ah.time_d * hr_per_d, obs_ah.conc_ugmL, seriestype = :scatter, alpha = 0.6, color = :blue, label = "AH (Bakri 2007)");
plot!(sol_mtk_ivt, idxs = pbpk_rabbit.eye.vh_igg_exg.C_VH * MW_IGG/1E3, label = "VH", lw = 2, alpha = 0.9, color = :red);
plot!(obs_vh.time_d * hr_per_d, obs_vh.conc_ugmL, seriestype = :scatter, alpha = 0.6, color = :red, label = "VH (Bakri 2007)");
plot!(sol_mtk_ivt, idxs = pbpk_rabbit.plasma_exg.C_Plasma * MW_IGG/1E3, label = "plasma", lw = 2, alpha = 0.9, color = :green);
plot!(obs_serum.time_d * hr_per_d, obs_serum.conc_ugmL, seriestype = :scatter, alpha = 0.6, color = :green, label = "Plasma (Bakri 2007)");
plot!(xlabel = "Time (hr)", xlims = [0, 504], title = "Dose = 1.25 mg, intravitreal injection", titlefontsize = 8);
display(plt__mAb)

savefig(plt__mAb, @projectroot("deliv/figure/bevacizumab-rabbit-ah-vh-mAb-ivt.pdf")) 


#####################
 
# sensitivity analysis on endothelial scaling factor 
param_1 = create_base_rabbit_pbpk_param(PS_Score = 5, mdl = pbpk_rabbit, Scale_Factor_Rabbit = 1); 
param_10 = create_base_rabbit_pbpk_param(PS_Score = 5, mdl = pbpk_rabbit, Scale_Factor_Rabbit = 10); 
param_100 = create_base_rabbit_pbpk_param(PS_Score = 5, mdl = pbpk_rabbit, Scale_Factor_Rabbit = 100); 

@time sol_ivt_1 = solve(ODEProblem(pbpk_rabbit, merge(Dict(u0_ivt), param_1), tspan), alg=CVODE_BDF());
@time sol_ivt_10 = solve(ODEProblem(pbpk_rabbit, merge(Dict(u0_ivt), param_10), tspan), alg=CVODE_BDF());
@time sol_ivt_100 = solve(ODEProblem(pbpk_rabbit, merge(Dict(u0_ivt), param_100), tspan), alg=CVODE_BDF());

plt__mAb__sens = plot(ylabel = "mAb Concentration (uM)", legend = :outerright, dpi = 300, size = (600, 400)); 
plot!(ylims = [4E-4, 9],  yaxis = :log10, yticks = [1E-3, 1E-2, 1E-1, 1]); 
plot!(sol_ivt_1, idxs = pbpk_rabbit.eye.ah_igg_exg.C_AQ, label = "AH, scaling factor = 1", lw = 2, alpha = 0.6, linestyle = :dash);
plot!(sol_ivt_10, idxs = pbpk_rabbit.eye.ah_igg_exg.C_AQ, label = "AH, scaling factor = 10", lw = 1.5, alpha = 0.7, linestyle = :dash);
plot!(sol_ivt_100, idxs = pbpk_rabbit.eye.ah_igg_exg.C_AQ, label = "AH, scaling factor = 100", lw = 1, alpha = 0.9, linestyle = :dash);
plot!(sol_ivt_1, idxs = pbpk_rabbit.eye.vh_igg_exg.C_VH, label = "VH, scaling factor = 1", lw = 2, alpha = 0.6, linestyle = :dashdotdot);
plot!(sol_ivt_10, idxs = pbpk_rabbit.eye.vh_igg_exg.C_VH, label = "VH, scaling factor = 10", lw = 1.5, alpha = 0.7, linestyle = :dashdotdot);
plot!(sol_ivt_100, idxs = pbpk_rabbit.eye.vh_igg_exg.C_VH, label = "VH, scaling factor = 100", lw = 1, alpha = 0.9, linestyle = :dashdotdot);
plot!(sol_ivt_1, idxs = pbpk_rabbit.plasma_exg.C_Plasma, label = "plasma, scaling factor = 1", lw = 2, alpha = 0.6);
plot!(sol_ivt_10, idxs = pbpk_rabbit.plasma_exg.C_Plasma, label = "plasma, scaling factor = 10", lw = 1.5, alpha = 0.7);
plot!(sol_ivt_100, idxs = pbpk_rabbit.plasma_exg.C_Plasma, label = "plasma, scaling factor = 100", lw = 1, alpha = 0.9);
plot!(xlabel = "Time (hr)", xlims = [0, 504], title = "Dose = 1.25 mg, intravitreal injection", titlefontsize = 8);
display(plt__mAb__sens)

savefig(plt__mAb__sens, @projectroot("deliv/figure/bevacizumab-rabbit-ivt-scaling-factor-sensitivity.pdf")) 
