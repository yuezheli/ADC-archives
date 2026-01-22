# date: 7/2/2025
# author: Yuezhe Li 
# purpose of this code: to generate heatmap for Cmax, Cavg, against affinity to FcRn/ membrane 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Parameters: @unpack
using DataFrames
using Statistics
using Plots
using LaTeXStrings
using CSV

include(@projectroot("Constants.jl"))
include(@projectroot("model/jones_homo.jl") )
include(@projectroot("model/param-pk.jl"))
include(@projectroot("model/init-pk.jl")) 
include(@projectroot("script/helper-half-life.jl"))

# use DM1 as payload 
p_0 = deepcopy(p_base);
p_0.PS_Score = -2
p_0.PS_kd = 0.001
p_0.init_sR = 0    # remove TMDD in the system for simplicity 

first_dose = 3.6  # [mg/kg]
tspan = (0., hr_per_day*21);      # [hr]

u0_tmp = jones_init(first_dose*1E3, p_0, BW, V_Plasma);  # IV bolus dose
prob_0 = ODEProblem(jonesODEs_homo_tumor!, u0_tmp, tspan, p_0);

PS_kd_values = 10 .^ range(-3.0, 3, step=0.2); # for mAb with PS_Score = 0, PS_Kd = 1.2E6; with PS_Score = 6, PS_Kd = 18.17; with PS_Score = 10, PS_Kd = 2.7
KD6_WT_values = 10 .^ range(-3.0, 3, step=0.2);

cmax = zeros(length(PS_kd_values), length(KD6_WT_values));
cavg = zeros(length(PS_kd_values), length(KD6_WT_values));
thalf = zeros(length(PS_kd_values), length(KD6_WT_values));

for i in 1:length(PS_kd_values)
    for j in 1:length(KD6_WT_values)
        p_tmp = deepcopy(p_0);
        p_tmp.PS_kd = PS_kd_values[i];
        p_tmp.KD6_WT = KD6_WT_values[j];
        prob_tmp = remake(prob_0, p = p_tmp);
        sol_tmp = solve(prob_tmp, saveat = 1., alg = QNDF(autodiff=false), reltol = 1E-12);
        sol_tmp_plasma = [sol_tmp.u[i].C_EXG_Plasma for i in 1:length(sol_tmp.t)]; # [uM]
        sol_tmp_endo_he = [sol_tmp.u[i].end_cyto_payload[2] for i in 1:length(sol_tmp.t)]; # [uM]
        # display( plot(plot(sol_tmp.t, sol_tmp_plasma), plot(sol_tmp.t, sol_tmp_endo_he), ncol = 2) ) # visual QC
        cmax_tmp = maximum(sol_tmp_endo_he)*1E3  # [nM]
        cavg_tmp = mean(sol_tmp_endo_he)*1E3  # [nM]
        cmax[i,j] = cmax_tmp
        cavg[i,j] = cavg_tmp
        thalf[i,j] = plasma_half_life(sol_tmp.t, sol_tmp_plasma)
    end
end

p_cmax = heatmap(PS_kd_values, KD6_WT_values, cmax, color =:jet1, xlabel=L"K$_D$, membrane ($\mu$M)", ylabel=L"K$_D$, FcRn (nM)", cbar_title=L"C$_{max}$ (nM)", clim=(2,28), axis = :log, size = (400, 400), dpi = 300)
p_cavg = heatmap(PS_kd_values, KD6_WT_values, cavg, color =:jet1, xlabel=L"K$_D$, membrane ($\mu$M)", ylabel=L"K$_D$, FcRn (nM)", cbar_title=L"C$_{avg}$ (nM)", axis = :log, size = (400, 400), dpi = 300)
p_thalf = heatmap(PS_kd_values, KD6_WT_values, thalf/hr_per_day, color =:jet1, xlabel=L"K$_D$, membrane ($\mu$M)", ylabel=L"K$_D$, FcRn (nM)", cbar_title=L"t$_{\frac{1}{2}}$ (Day)", axis = :log, size = (400, 400), dpi = 300)

savefig(p_cmax, @projectroot("deliv/figure/pl/sens-membrane-fcrn-kd-cmax.png"));
savefig(p_cavg, @projectroot("deliv/figure/pl/sens-membrane-fcrn-kd-cavg.png"));
savefig(p_thalf, @projectroot("deliv/figure/pl/sens-membrane-fcrn-kd-thalf.png"));
