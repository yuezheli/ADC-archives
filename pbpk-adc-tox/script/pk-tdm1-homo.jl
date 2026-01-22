# date: 1/13/2026 
# author: Yuezhe Li 
# purpose of this code: PK fit of T-DM1

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

#####################

# constant function
include(@projectroot("model/constants.jl"));

# model functions
include(@projectroot("model/mab-pbpk-modular-eye.jl"));
include(@projectroot("model/default-parameters.jl"));

# helper functions
include(@projectroot("script/helper-initiation.jl")); 
include(@projectroot("script/helper-infusion.jl"));
include(@projectroot("script/helper-visualization.jl"));

#####################

# observed data obtained from Girish et al., 2012; # https://link.springer.com/article/10.1007/s00280-011-1817-3
obs_tdm1 = DataFrame(
  time_day = [0.1,0.25,1,2,3,0.1,0.25,1,2,3,7,0.1,0.25,1,2,3,7,0.25,1,2,3,7,10,14,17,21,0.1,0.25,1,2,3,7,10,17,21,0.1,0.25,1,2,3,7,14,17,21], 
  T_DM1_ugperml = [9.66134,7.1674,4.22606,2.54967,1.11527,13.02307,10.35053,7.33393,5.69654,1.37138,0.52264,20.14835,20.14835,11.08888,7.33393,3.208,0.49918,74.61359,55.35313,45.01597,33.39572,17.15592,10.35053,3.68201,2.02644,0.40596,74.61359,66.51853,48.22717,34.96559,27.15909,13.95206,9.44196,2.27305,0.77229,126.54482,118.11885,89.66407,78.12105,65.00812,29.09648,12.43836,7.1674,3.3588],
  Dose_mgkg = [0.3,0.3,0.3,0.3,0.3,0.6,0.6,0.6,0.6,0.6,0.6,1.2,1.2,1.2,1.2,1.2,1.2,2.4,2.4,2.4,2.4,2.4,2.4,2.4,2.4,2.4,3.6,3.6,3.6,3.6,3.6,3.6,3.6,3.6,3.6,4.8,4.8,4.8,4.8,4.8,4.8,4.8,4.8,4.8]
); 
@transform!(obs_tdm1, :tdm1_nM = :T_DM1_ugperml*1E6/15E4); 

# convert observed data to dictionary 
pk_tdm1 = Dict();
for dose_ in unique(obs_tdm1.Dose_mgkg)
    tmp_pk_ = @rsubset(obs_tdm1, :Dose_mgkg == dose_);
    tmp_df = DataFrame(
        time_d = tmp_pk_.time_day, 
        ADC_uM = tmp_pk_.T_DM1_ugperml*1E3/MW_TDM1
    )
    push!(pk_tdm1, dose_ => tmp_df)
end

obs_dm1_pk = DataFrame(
  time_d = [0, 0.005, 7, 14, 21], 
  dm1_ugL = [0.325, 6.071, 0.599, 0.418, 0.428]
  );
@transform!(obs_dm1_pk, :dm1_uM = :dm1_ugL/MW_Lys_MCC_DM1); 

# PK data from Krop et al., 2010; https://pubmed.ncbi.nlm.nih.gov/20421541/ (3.6 mg/kg)
obs_dm1_pk2 = DataFrame(
  time_hr = [1.2, 6, 25.3, 49.4, 73.5, 96.4, 168.7], 
  dm1_ugL = [4.46, 3.32, 1.91, 1.2, 1.01, 0.93, 0.89]
  );
@transform!(obs_dm1_pk2, :time_d = :time_hr/hr_per_day); 
@transform!(obs_dm1_pk2, :dm1_uM = :dm1_ugL/MW_Lys_MCC_DM1); 

obs_tab_pk = DataFrame(
  time_hr = [6.02, 25.3, 49.4, 72.29, 95.18, 168.67, 240.96, 337.35, 408.43, 504.82], 
  tdm1_ugL = [67149.25, 47847.09, 34081.82, 28742.54, 21357.67, 13933.59, 9482.74, 3280.01, 2328.67, 772.19], 
  tAb_ugL = [102405.67, 79393.99, 61527.22, 58885.28, 47608.51, 41738.78, 32236.47, 29425.09, 25799.41, 16805.71]
  ); 
@transform!(obs_tab_pk, :time_d = :time_hr/hr_per_day);
@transform!(obs_tab_pk, :tdm1_uM = :tdm1_ugL/MW_TDM1); 
@transform!(obs_tab_pk, :tAb_uM = :tAb_ugL/MW_TDM1);


#####################

## create pbpk model
@time pbpk = create_pbpk();
@time pbpk_simple = mtkcompile(pbpk);

## initial condition (IV dosing)
u0_infusion = pbpk_initial_condition(0, pbpk_simple); 

## parameters 
### AC-SINS score, VMAX, KM all from T-Dxd fitting
param_tdm1 = create_base_pbpk_param(5, pbpk_simple, k_PL_ints_clearance = 0.34);  #  https://pubmed.ncbi.nlm.nih.gov/37787918/
param_tdm1[pbpk_simple.VMAX] = 0.0012
param_tdm1[pbpk_simple.KM] = 0.006
param_tdm1[pbpk_simple.DAR] = 3.5
param_tdm1[pbpk_simple.k_diff] = 0.14  #  https://pubmed.ncbi.nlm.nih.gov/37787918/
param_tdm1[pbpk_simple.k_deconj] =  6E-3
param_tdm1[pbpk_simple.CL_plasma_PL] =  4
param_tdm1[pbpk_simple.Kp_adc_cor_ah] = 1


# simulation time
tspan = (-0.01, hr_per_day*84);      # [hr]
add_dose = [0., 21., 42., 63]   

# simulation dose 
sims_dose = [0.3, 0.6, 1.2, 2.4, 3.6, 4.8]; # [mg/kg]

# simulation
sol_pk_tmd1 = Dict(); 
for init_dose in sims_dose
    adc_infusion_ = InfusionCallback(init_dose, pbpk_simple, infusion_d = add_dose, BW = 90);
    # simulation 
    @time prob_mtk_ = ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_tdm1), tspan, callback = adc_infusion_, infusion_hr = 0.1); 
    @time sol_mtk_infusion_ = solve(prob_mtk_, alg=CVODE_BDF());
    # append result 
    push!(sol_pk_tmd1, init_dose => sol_mtk_infusion_);
end

# simulation (no deconjugation)
### parameters without deconjugation 
param_tAb = deepcopy(param_tdm1);
param_tAb[pbpk_simple.k_deconj] = 0
### simulation of 3.6 mg/kg
adc_infusion_3point6 = InfusionCallback(6.4, pbpk_simple, infusion_d = add_dose, BW = 90, infusion_hr = 0.5);
prob_tab_3point6 = ODEProblem(pbpk_simple, merge(Dict(u0_infusion), param_tAb), tspan, callback = adc_infusion_3point6); 
@time sol_tab_3point6 = solve(prob_tab_3point6, alg=CVODE_BDF());

#####################
# PK, ADC
plt_plasma_adc = PlotSimulationPlasma(sol_pk_tmd1, pk_tdm1, pbpk_simple, adc_name = "trastuzumab emtansine", colorPALETTE = :batlowKS, xrange = [0, 42], yrange = [1E-4, 10], ylog = true, legendcolumnsnum = 2)

# PK, DM1
plt_plasma_pl = plot(xlabel = "Time (d)", ylabel = "Concentration (uM)", size = (400, 400), dpi = 300); 
plot!(ylims = [1E-4, 1E2], yaxis = :log10, xlims = [0, 21], xticks = [0, 7, 14, 21]);                  
plot!(sol_pk_tmd1[3.6].t/hr_per_day, sol_pk_tmd1[3.6][pbpk_simple.plasma_exg.C_Plasma], label = "sims, T-DM1", color = palette(:batlowKS)[5], lw = 2, linestyle = :dot); 
scatter!(obs_tab_pk.time_d, obs_tab_pk.tdm1_uM, label = "Krop et al., 2010; T-DM1", color = palette(:batlowKS)[5], alpha = 0.6, markershape = :star5);
plot!(sol_tab_3point6.t/hr_per_day, sol_tab_3point6[pbpk_simple.plasma_exg.C_Plasma], label = "sims, total Ab", color = palette(:batlowKS)[5], lw = 2); 
scatter!(obs_tab_pk.time_d, obs_tab_pk.tAb_uM, label = "Krop et al., 2010; total Ab", color = palette(:batlowKS)[5], alpha = 0.6, markershape = :utriangle);
plot!(sol_pk_tmd1[3.6].t/hr_per_day, sol_pk_tmd1[3.6][pbpk_simple.plasma_pl.C_PL_Plasma], label = "sims, DM1", color = palette(:batlowKS)[5], lw = 2, linestyle = :dash); 
scatter!(obs_dm1_pk.time_d, obs_dm1_pk.dm1_uM, label = "Girish et al., 2012; DM1", color = palette(:batlowKS)[5], alpha = 0.6); 
scatter!(obs_dm1_pk2.time_d, obs_dm1_pk2.dm1_uM, label = "Krop et al., 2010; DM1", color = palette(:batlowKS)[5], alpha = 0.6, markershape = :diamond); 

display(plot(plt_plasma_adc, plt_plasma_pl, layout = (1,2), size = (800, 400), margin = 5mm))

# save pk figures  
savefig(plt_plasma_adc, "deliv/figure/pk-tdm1-adc-plasma.png");
savefig(plt_plasma_pl, "deliv/figure/pk-tdm1-dm1-plasma.png");

#####################
# visualization, liver PL 
plt_liver_pl = plot(xlabel = "Time (d)", ylabel = "DM1 (uM)", 
                     ylims = [1E-4, 1E2], yaxis = :log10, xlims = [0, 42], xticks = [0, 7, 14, 21, 28, 35, 42],
                     size = (400, 400), dpi = 300); 
plot!(sol_pk_tmd1[3.6].t/hr_per_day, sol_pk_tmd1[3.6][pbpk_simple.liver.PL_tissue.C_PL_IntS], label = "liver, interstitium, 3.6 mg/kg", linestyle = :dashdot, color = :red); 
plot!(sol_pk_tmd1[3.6].t/hr_per_day, sol_pk_tmd1[3.6][pbpk_simple.liver.PL_tissue.C_PL_endo], label = "liver, endothelial cells, 3.6 mg/kg", linestyle = :solid, color = :red); 
plot!(sol_pk_tmd1[2.4].t/hr_per_day, sol_pk_tmd1[2.4][pbpk_simple.liver.PL_tissue.C_PL_IntS], label = "liver, interstitium, 2.4 mg/kg", linestyle = :dashdot, color = :blue); 
plot!(sol_pk_tmd1[2.4].t/hr_per_day, sol_pk_tmd1[2.4][pbpk_simple.liver.PL_tissue.C_PL_endo], label = "liver, endothelial cells, 2.4 mg/kg", linestyle = :solid, color = :blue); 
plot!(sol_pk_tmd1[1.2].t/hr_per_day, sol_pk_tmd1[1.2][pbpk_simple.liver.PL_tissue.C_PL_IntS], label = "liver, interstitium, 1.2 mg/kg", linestyle = :dashdot, color = :green); 
plot!(sol_pk_tmd1[1.2].t/hr_per_day, sol_pk_tmd1[1.2][pbpk_simple.liver.PL_tissue.C_PL_endo], label = "liver, endothelial cells, 1.2 mg/kg", linestyle = :solid, color = :green); 
display(plt_liver_pl)

savefig(plt_liver_pl, "deliv/figure/pk-tdm1-dm1-liver.png");

