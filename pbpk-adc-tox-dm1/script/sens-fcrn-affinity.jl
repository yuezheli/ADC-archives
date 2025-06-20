# date:2/4/2025
# author: Yuezhe Li 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)

using DifferentialEquations, ComponentArrays
using Plots, DataFrames, LaTeXStrings

include(@projectroot("model/jones_homo.jl") );
include(@projectroot("model/param-pk.jl")); 
include(@projectroot("model/init-pk.jl")); 
include(@projectroot("script/helper.jl"));

p_tdm1 = deepcopy(p_base);
p_tdm1.PS_Score = 6.
p_tdm1.init_sR = 0.004


first_dose = 3.6  # [mg/kg]
tspan = (0., DayToHour*21);      # [hr]

u0_tdm1, _ = jones_init(first_dose*1E3, p_tdm1);
prob_tdm1 = ODEProblem(jonesODEs_homo_tumor!, u0_tdm1, tspan, p_tdm1);

sol_tdm1 = tissue_ints_endo(solve(prob_tdm1, saveat = 1., alg = QNDF(autodiff=false), reltol = 1E-12));
sol_tdm1.KD6_WT .= 700

KD6_WT_values = [0.1, 1, 10, 100, 1000];

for tmp_kd6_wt in KD6_WT_values
    p_tmp = deepcopy(p_tdm1);
    p_tmp.KD6_WT = tmp_kd6_wt;
    prob_tmp = remake(prob_tdm1, p = p_tmp);
    sol_tmp = tissue_ints_endo(solve(prob_tmp, saveat = 1., alg = QNDF(autodiff=false), reltol = 1E-12));
    sol_tmp.KD6_WT .= tmp_kd6_wt;
    sol_tdm1 = vcat(sol_tdm1, sol_tmp);
end

@rsubset!(sol_tdm1, :time > 0);

p_fcrn_kd_he = plot(dpi = 1000, size = (400, 400), palette = :Spectral_6);
plot!(xlabel = "Time (day)", ylabel = "Liver endothelial cytosol payload (nM)", xticks = [0, 7, 14, 21]);
plot!(sol_tdm1.time/DayToHour, sol_tdm1.endo_he*1E3, group = sol_tdm1.KD6_WT, alpha = 0.9);
plot!(legendtitle = L"K_D, FcRn (nM)", legendtitlefontsize = 8, legendfontsize = 8);
display(p_fcrn_kd_he)

p_fcrn_kd_pk = plot(yaxis = :log10, yticks = [1, 10, 100, 1000], ylims = [1, 1000], dpi = 1000, size = (400, 400), palette = :Spectral_6);
plot!(xlabel = "Time (day)", ylabel = "Plasma ADC (nM)", xticks = [0, 7, 14, 21]);
plot!(sol_tdm1.time/DayToHour, sol_tdm1.adcplasma*1E3, group = sol_tdm1.KD6_WT, alpha = 0.9);
plot!(legendtitle = L"K_D, FcRn (nM)", legendtitlefontsize = 8, legendfontsize = 8, legend = :bottomleft);
display(p_fcrn_kd_pk)

savefig(p_fcrn_kd_pk, @projectroot("deliv/figure/payload/sens-fcrn-kd-pk.png"));
savefig(p_fcrn_kd_he, @projectroot("deliv/figure/payload/sens-fcrn-kd-he.png"));
