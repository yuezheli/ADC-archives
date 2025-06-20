# date: 2/4/2025
# author: Yuezhe Li 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)

using DifferentialEquations, ComponentArrays
using Plots, DataFrames, LaTeXStrings

include(@projectroot("model/jones_homo.jl") );
include(@projectroot("model/param-pk.jl")); 
include(@projectroot("model/init-pk.jl")); 
include(@projectroot("script/helper.jl"));

p_0 = deepcopy(p_base);
p_0.PS_Score = -2
p_0.PS_kd = 0.001
p_0.init_sR = 0.004

first_dose = 3.6  # [mg/kg]
tspan = (0., DayToHour*21);      # [hr]

u0_tdm1, _ = jones_init(first_dose*1E3, p_0);
prob_tdm1 = ODEProblem(jonesODEs_homo_tumor!, u0_tdm1, tspan, p_0);

sol_tdm1 = tissue_ints_endo(solve(prob_tdm1, saveat = 1., alg = QNDF(autodiff=false), reltol = 1E-12));
sol_tdm1.PS_kd .= 0.001; 

# PS_Score = 15 --> PS_Kd = 1.31 --> IMGN853
# PS_Score = 8 --> PS_Kd = 5.56 --> SAR408701
# PS_Score = 6 --> PS_Kd = 18.17 --> T-DM1
PS_kd_values = [0.1, 1., 10, 100];

for tmp_ps_kd in PS_kd_values
    p_tmp = deepcopy(p_0);
    p_tmp.PS_kd = tmp_ps_kd;
    prob_tmp = remake(prob_tdm1, p = p_tmp);
    sol_tmp = tissue_ints_endo(solve(prob_tmp, saveat = 1., alg = QNDF(autodiff=false), reltol = 1E-12));
    sol_tmp.PS_kd .= tmp_ps_kd;
    sol_tdm1 = vcat(sol_tdm1, sol_tmp);
end

@rsubset!(sol_tdm1, :time > 0);

p_membrane_kd_he = plot(dpi = 1000, size = (400, 400), guidefontsize = 10, palette = :Set1_5);
plot!(xlabel = "Time (day)", ylabel = "Liver endothelial cytosol payload (nM)", xticks = [0, 7, 14, 21]);
plot!(sol_tdm1.time/DayToHour, sol_tdm1.endo_he*1E3, group = sol_tdm1.PS_kd, alpha = 0.9);
plot!(legendtitle = L"K$_D$ ($\mu$M)", legendtitlefontsize = 8, legendfontsize = 8);
display(p_membrane_kd_he)

p_membrane_kd_pk = plot(dpi = 1000, size = (400, 400), guidefontsize = 10, palette = :Set1_5, legend = :bottomleft);
plot!(xlabel = "Time (day)", ylabel = "Plasma ADC (nM)", xticks = [0, 7, 14, 21], yaxis = :log10, yticks = [1, 10, 100, 1000], ylims = [1, 1000]);
plot!(sol_tdm1.time/DayToHour, sol_tdm1.adcplasma*1E3, group = sol_tdm1.PS_kd, alpha = 0.9);
plot!(legendtitle = L"K$_D$ ($\mu$M)", legendtitlefontsize = 8, legendfontsize = 8);
display(p_membrane_kd_pk)

savefig(p_membrane_kd_pk, @projectroot("deliv/figure/payload/sens-membrane-kd-pk.png"));
savefig(p_membrane_kd_he, @projectroot("deliv/figure/payload/sens-membrane-kd-he.png"));

