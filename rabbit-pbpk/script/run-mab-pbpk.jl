using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using DifferentialEquations, ModelingToolkit
using Plots
using Plots.Measures
using DataFrames
using SymbolicIndexingInterface  # for parameter_index()
# using CSV
using Sundials  # for CVODE_BDF()

#####################
# helper functions
include(@projectroot("script/helper-funcs.jl"));
include(@projectroot("script/helper-initiation.jl")); 
include(@projectroot("script/helper-parameters.jl"));
include(@projectroot("script/helper-infusion.jl"));

# model functions
include(@projectroot("model/mab-pbpk-modular-eye.jl"));

#####################

## create pbpk model
@time pbpk = create_pbpk();

## checks
#hierarchy(pbpk)
#equations(pbpk)
#unknowns(pbpk)
#parameters(pbpk)

@time pbpk_simple = mtkcompile(pbpk);

# check
# filter(var -> occursin("C_Plasma(t)) ~", string(var)), equations(pbpk_simple))

######################

#---# assign parameter values #---#
# mAb AC-SINS score
ACSINS_score = 0.0;
# mAb molecular weight
MW_EDG = 150000.0; # Da

p_map_all = create_base_pbpk_param(ACSINS_score, pbpk_simple);

######################

#---# set initial conditions #---#

# mAb IV Dose
Dose_in_mgkg = 5.0; # [mg/kg]

u0_map = pbpk_initial_condition(Dose_in_mgkg, pbpk_simple); 

######################

# Define the simulation timespan
tspan = (0.0,1008.0);  # hr

######################

#---# create ODE problem #---#
u0_p = merge(Dict(u0_map), p_map_all);
@time prob_mtk = ODEProblem(pbpk_simple, u0_p, tspan);  


#####################

#---#  simulate  #---#
@time sol_mtk = solve(prob_mtk, alg=CVODE_BDF());

plt_PL = plot(
    plot(sol_mtk, idxs = pbpk_simple.plasma_pl.C_PL_Plasma, title = "Plasma", ylabel = "Conc (uM)", xlabel = "Time (hr)"),
    plot(sol_mtk, idxs = pbpk_simple.lung.PL_tissue.C_PL_endo, title = "Endothelical PL, Lung", ylabel = "Conc (uM)", xlabel = "Time (hr)"),
    plot(sol_mtk, idxs = pbpk_simple.liver.PL_tissue.C_PL_endo, title = "Endothelical PL, Liver", ylabel = "Conc (uM)", xlabel = "Time (hr)"),
    plot(sol_mtk, idxs = pbpk_simple.marrow.PL_tissue.C_PL_endo, title = "Endothelical PL, Bone Marrow", ylabel = "Conc (uM)", xlabel = "Time (hr)"),
    plot(sol_mtk, idxs = pbpk_simple.eye.icb.PL_tissue.C_PL_endo, title = "Endothelical PL, iris-ciliary body", ylabel = "Conc (uM)", xlabel = "Time (hr)"),
    plot(sol_mtk, idxs = pbpk_simple.eye.retina.PL_tissue.C_PL_endo, title = "Endothelical PL, retina", ylabel = "Conc (uM)", xlabel = "Time (hr)"),
    plot(sol_mtk, idxs = pbpk_simple.eye.choroid.PL_tissue.C_PL_endo, title = "Endothelical PL, choroid", ylabel = "Conc (uM)", xlabel = "Time (hr)"),
    layout = (2,4), margin = 5mm, size = (1600, 800)
)

#####################

#---#  infusion  #---#
u0_infusion = pbpk_initial_condition(0, pbpk_simple);
u0_p_infusion = merge(Dict(u0_infusion), p_map_all)

adc_infusion = InfusionCallback(Dose_in_mgkg, pbpk_simple);

@time prob_mtk = ODEProblem(pbpk_simple, u0_p_infusion, tspan, callback = adc_infusion);  

@time sol_mtk_infusion = solve(prob_mtk, alg=CVODE_BDF());

# plot mAb concentration,  
plt_infusion = plot(
    
    plot(sol_mtk_infusion, idxs = pbpk_simple.marrow.igg_exg.C_IntS, title = "Bone Marrow Interstitium", ylabel = "Conc (uM)", xlabel = "Time (hr)", label = false, titlefontsize = 8),
    plot(sol_mtk_infusion, idxs = pbpk_simple.liver.igg_exg.C_IntS, title = "Liver Interstitium", ylabel = "Conc (uM)", xlabel = "Time (hr)", label = false, titlefontsize = 8),
    plot(sol_mtk_infusion, idxs = pbpk_simple.eye.icb.igg_exg.C_IntS, title = "ICB Interstitium", ylabel = "Conc (uM)", xlabel = "Time (hr)", label = false, titlefontsize = 8),
    plot(sol_mtk_infusion, idxs = pbpk_simple.eye.retina.igg_exg.C_IntS, title = "Retina Interstitium", ylabel = "Conc (uM)", xlabel = "Time (hr)", label = false, titlefontsize = 8),
    plot(sol_mtk_infusion, idxs = pbpk_simple.eye.choroid.igg_exg.C_IntS, title = "Choroid Interstitium", ylabel = "Conc (uM)", xlabel = "Time (hr)", label = false, titlefontsize = 8),

    plot(sol_mtk_infusion, idxs = pbpk_simple.plasma_exg.C_Plasma, title = "Plasma", ylabel = "Conc (uM)", xlabel = "Time (hr)", label = false, titlefontsize = 8),

    plot(sol_mtk_infusion, idxs = pbpk_simple.eye.ah_igg_exg.C_AQ, title = "aqueous humor", ylabel = "Conc (uM)", xlabel = "Time (hr)", label = false, titlefontsize = 8),
    plot(sol_mtk_infusion, idxs = pbpk_simple.eye.igg_exg.C_COR, title = "cornea", ylabel = "Conc (uM)", xlabel = "Time (hr)", label = false, titlefontsize = 8),
    # plot(sol_mtk_infusion, idxs = pbpk_simple.eye.lens_igg_exg.C_LENS, title = "lens", ylabel = "Conc (uM)", xlabel = "Time (hr)", label = false, titlefontsize = 8),
    plot(sol_mtk_infusion, idxs = pbpk_simple.eye.vh_igg_exg.C_VH, title = "vitreous humor", ylabel = "Conc (uM)", xlabel = "Time (hr)", label = false, titlefontsize = 8),

    layout = (3,3), margin = 5mm, size = (800, 800)
)

savefig(plt_infusion, @projectroot("deliv/figure/verification/iv-dosing.pdf"))

# plot(sol_mtk, idxs = (pbpk_simple.lung.PLQ + pbpk_simple.lung.LF) * pbpk_simple.plasma_exg.C_Plasma)

plt_infusion_overlay = plot(ylabel = "Conc (nM)", ylims = [4E-5, 1],  yaxis = :log10, legend = :bottomright); 
plot!(sol_mtk_infusion, idxs = pbpk_simple.plasma_exg.C_Plasma, label = "Plasma", linestyle = :dash, lw = 3);
plot!(sol_mtk_infusion, idxs = pbpk_simple.marrow.igg_exg.C_IntS, label = "bone marrow Interstitium", linestyle = :dashdotdot);
plot!(sol_mtk_infusion, idxs = pbpk_simple.liver.igg_exg.C_IntS, label = "liver Interstitium", linestyle = :dashdotdot);
plot!(sol_mtk_infusion, idxs = pbpk_simple.eye.icb.igg_exg.C_IntS, label = "ICB Interstitium", linestyle = :dashdotdot, lw = 2, alpha = 0.8);
plot!(sol_mtk_infusion, idxs = pbpk_simple.eye.retina.igg_exg.C_IntS, label = "Retina Interstitium", linestyle = :dashdotdot, lw = 2, alpha = 0.8);
plot!(sol_mtk_infusion, idxs = pbpk_simple.eye.choroid.igg_exg.C_IntS, label = "Choroid Interstitium", linestyle = :dashdotdot, lw = 2, alpha = 0.8);
plot!(sol_mtk_infusion, idxs = pbpk_simple.eye.ah_igg_exg.C_AQ, label = "aqueous humor", lw = 2, alpha = 0.7);
plot!(sol_mtk_infusion, idxs = pbpk_simple.eye.igg_exg.C_COR, label = "cornea", lw = 2, alpha = 0.7);
plot!(sol_mtk_infusion, idxs = pbpk_simple.eye.vh_igg_exg.C_VH, label = "vitreous humor", lw = 2, alpha = 0.7);
plot!(xlabel = "Time (hr)");
display(plt_infusion_overlay)

savefig(plt_infusion_overlay, @projectroot("deliv/figure/verification/iv-dosing-overlay.pdf"))

#####################

# IVT dosing 
# add drug into vitreous humor, 1.25 mg 
DOSE_IVT = 1.25
V_VH = 4E-3 # 
u0_ivt = pbpk_initial_condition(0, pbpk_simple);
u0_ivt[pbpk_simple.eye.vh_igg_exg.C_VH.val] = DOSE_IVT*1E-3/V_VH/MW_EDG*1e6; # [uM]
p_tmp = deepcopy(p_map_all); 

tspan2 = (0.0,2800.0);  # hr

@time prob_mtk_ivt = ODEProblem(pbpk_simple, merge(Dict(u0_ivt), p_tmp), tspan2);  
@time sol_mtk_ivt = solve(prob_mtk_ivt, alg=CVODE_BDF());

plt_ivt = plot(
    plot(sol_mtk_ivt, idxs = pbpk_simple.plasma_exg.C_Plasma, title = "Plasma", ylabel = "Conc (uM)", xlabel = "Time (hr)", label = false),
    plot(sol_mtk_ivt, idxs = pbpk_simple.eye.icb.igg_exg.C_IntS, title = "ICB Interstitium", ylabel = "Conc (uM)", xlabel = "Time (hr)", label = false),
    plot(sol_mtk_ivt, idxs = pbpk_simple.eye.retina.igg_exg.C_IntS, title = "Retina Interstitium", ylabel = "Conc (uM)", xlabel = "Time (hr)", label = false),
    plot(sol_mtk_ivt, idxs = pbpk_simple.eye.choroid.igg_exg.C_IntS, title = "Choroid Interstitium", ylabel = "Conc (uM)", xlabel = "Time (hr)", label = false),

    plot(sol_mtk_ivt, idxs = pbpk_simple.eye.ah_igg_exg.C_AQ, title = "aqueous humor", ylabel = "Conc (uM)", xlabel = "Time (hr)", label = false),
    plot(sol_mtk_ivt, idxs = pbpk_simple.eye.igg_exg.C_COR, title = "cornea", ylabel = "Conc (uM)", xlabel = "Time (hr)", label = false),
    # plot(sol_mtk_ivt, idxs = pbpk_simple.eye.lens_igg_exg.C_LENS, title = "lens", ylabel = "Conc (uM)", xlabel = "Time (hr)", label = false),
    plot(sol_mtk_ivt, idxs = pbpk_simple.eye.vh_igg_exg.C_VH, title = "vitreous humor", ylabel = "Conc (uM)", xlabel = "Time (hr)", label = false),

    layout = (2,4), margin = 10mm, size = (2000, 1000)
)

savefig(plt_ivt, @projectroot("deliv/figure/verification/ivt-dosing.pdf"))

plt_ivt_overlay = plot(ylabel = "Conc (nM)", xlabel = "Time (hr)", ylims = [1E-5, 1],  yaxis = :log10); 
plot!(sol_mtk_ivt, idxs = pbpk_simple.plasma_exg.C_Plasma, label = "Plasma");
plot!(sol_mtk_ivt, idxs = pbpk_simple.eye.icb.igg_exg.C_IntS, label = "ICB Interstitium");
plot!(sol_mtk_ivt, idxs = pbpk_simple.eye.retina.igg_exg.C_IntS, label = "Retina Interstitium");
plot!(sol_mtk_ivt, idxs = pbpk_simple.eye.choroid.igg_exg.C_IntS, label = "Choroid Interstitium");
plot(sol_mtk_ivt, idxs = pbpk_simple.eye.ah_igg_exg.C_AQ, label = "aqueous humor");
plot(sol_mtk_ivt, idxs = pbpk_simple.eye.igg_exg.C_COR, label = "cornea");
plot(sol_mtk_ivt, idxs = pbpk_simple.eye.vh_igg_exg.C_VH, label = "vitreous humor");
plt_ivt_overlay



## change visualization to match Naware 2024 
plt_ivt_plasma = plot(sol_mtk_ivt, idxs = pbpk_simple.plasma_exg.C_Plasma * 1E3, title = "Plasma", ylabel = "Conc (nM)", xlabel = "Time (hr)", label = false, ylims = [1E-4, 1E2], yaxis = :log10, xlims = [0, 600], lw = 4)
plt_ivt_vh = plot(sol_mtk_ivt, idxs = pbpk_simple.eye.vh_igg_exg.C_VH * 1E3, title = "vitreous humor", ylabel = "Conc (nM)", xlabel = "Time (hr)", label = false, ylims = [1E-3, 1E5], yticks = [1E-3, 0.1, 10, 1E3, 1E5], yaxis = :log10, xlims = [0, 2800], xticks = [0, 500, 1000, 1500, 2000, 2500], lw = 4)


# add drug into vitreous humor, 1.5 mg 
u0_ivt2 = deepcopy(u0_ivt); 
u0_ivt2[pbpk_simple.eye.vh_igg_exg.C_VH.val] = 1.5*1E-3/V_VH/MW_EDG*1e6; # [uM]

@time prob_mtk_ivt2 = ODEProblem(pbpk_simple, merge(Dict(u0_ivt2), p_tmp), tspan2);  
@time sol_mtk_ivt2 = solve(prob_mtk_ivt2, alg=CVODE_BDF());

plt_ivt_ah = plot(sol_mtk_ivt2, idxs = pbpk_simple.eye.ah_igg_exg.C_AQ * 1E3, title = "aqueous humor", ylabel = "Conc (nM)", xlabel = "Time (hr)", label = false, ylims = [1E-1, 1E4], yticks = [0.1, 10, 1E3], yaxis = :log10, xlims = [0, 2000], xticks = [0, 500, 1000, 1500, 2000], lw = 4)

# change reflection coefficient of AH 
p_tmp2 = deepcopy(p_map_all); 
p_tmp2[pbpk_simple.Q_PtA] = 1.8E-4;

@time prob_mtk_ivt3 = ODEProblem(pbpk_simple, merge(Dict(u0_ivt2), p_tmp2), tspan2);  
@time sol_mtk_ivt3 = solve(prob_mtk_ivt3, alg=CVODE_BDF());

plt_ivt3_ah = plot(sol_mtk_ivt3, idxs = pbpk_simple.eye.ah_igg_exg.C_AQ * 1E3, title = "aqueous humor", ylabel = "Conc (nM)", xlabel = "Time (hr)", label = false, ylims = [1E-1, 1E4], yticks = [0.1, 10, 1E3], yaxis = :log10, xlims = [0, 2000], xticks = [0, 500, 1000, 1500, 2000], lw = 4)
