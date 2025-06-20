# date: 2/10/2025
# author: Yuezhe Li 
# purpose of this code: to generate heatmap for Cmax, Cavg, against affinity to FcRn/ membrane 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)

using DifferentialEquations, ComponentArrays
using Plots, DataFrames, LaTeXStrings
using Statistics

include(@projectroot("model/jones_homo.jl") );
include(@projectroot("model/param-pk.jl")); 
include(@projectroot("model/init-pk.jl")); 
include(@projectroot("script/helper.jl"));

# use DM1 as payload 
p_0 = deepcopy(p_base);
p_0.PS_Score = -2
p_0.PS_kd = 0.001
p_0.init_sR = 0.004

first_dose = 3.6  # [mg/kg]
tspan = (0., DayToHour*21);      # [hr]

u0_tdm1, _ = jones_init(first_dose*1E3, p_0);
prob_tdm1 = ODEProblem(jonesODEs_homo_tumor!, u0_tdm1, tspan, p_0);

PS_kd_values = 10 .^ range(-3.0, 3, step=0.2);
KD6_WT_values = 10 .^ range(-3.0, 3, step=0.2);

cmax = zeros(length(PS_kd_values), length(KD6_WT_values));
cavg = zeros(length(PS_kd_values), length(KD6_WT_values));

for i in 1:length(PS_kd_values)
    for j in 1:length(KD6_WT_values)
        p_tmp = deepcopy(p_0);
        p_tmp.PS_kd = PS_kd_values[i];
        p_tmp.KD6_WT = KD6_WT_values[j];
        prob_tmp = remake(prob_tdm1, p = p_tmp);
        sol_tmp = tissue_ints_endo(solve(prob_tmp, saveat = 1., alg = QNDF(autodiff=false), reltol = 1E-12));
        cmax_tmp = maximum(sol_tmp.endo_he)*1E3  # [nM]
        cavg_tmp = mean(sol_tmp.endo_he)*1E3  # [nM]
        cmax[i,j] = cmax_tmp
        cavg[i,j] = cavg_tmp
    end
end

p_cmax = heatmap(PS_kd_values, KD6_WT_values, cmax, color =:jet1, xlabel=L"K$_D$, membrane ($\mu$M)", ylabel=L"K$_D$, FcRn (nM)", cbar_title=L"C$_{max}$ (nM)", clim=(2,28), axis = :log, size = (400, 400), dpi = 1000)
p_cavg = heatmap(PS_kd_values, KD6_WT_values, cavg, color =:jet1, xlabel=L"K$_D$, membrane ($\mu$M)", ylabel=L"K$_D$, FcRn (nM)", cbar_title=L"C$_{avg}$ (nM)", axis = :log, size = (400, 400), dpi = 1000)

savefig(p_cmax, @projectroot("deliv/figure/payload/sens-membrane-fcrn-kd-cmax.png"));
savefig(p_cavg, @projectroot("deliv/figure/payload/sens-membrane-fcrn-kd-cavg.png"));
