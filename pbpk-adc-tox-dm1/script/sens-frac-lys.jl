# date: 2/5/2025
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
sol_tdm1.frac_lys .= 1;

frac_lys_values = [0.1, 0.3, 0.5, 0.7, 0.9];

for tmp_frac_lys in frac_lys_values
    p_tmp = deepcopy(p_tdm1);
    p_tmp.frac_lys = tmp_frac_lys;
    prob_tmp = remake(prob_tdm1, p = p_tmp);
    sol_tmp = tissue_ints_endo(solve(prob_tmp, saveat = 1., alg = QNDF(autodiff=false), reltol = 1E-12));
    sol_tmp.frac_lys .= tmp_frac_lys;
    sol_tdm1 = vcat(sol_tdm1, sol_tmp);
end

@rsubset!(sol_tdm1, :time > 0);

p_frac_lys_kd_he = plot(dpi = 1000, size = (400, 400), guidefontsize = 10, palette = :Set2_6, ylims = [0.01, 10], yaxis = :log);
plot!(xlabel = "Time (day)", ylabel = "Liver endothelial cytosol payload (nM)", xticks = [0, 7, 14, 21]);
plot!(sol_tdm1.time/DayToHour, sol_tdm1.endo_he*1E3, group = sol_tdm1.frac_lys, alpha = 0.9);
plot!(legendtitle = "Fraction of payload lysosomal release", legendtitlefontsize = 8, legendfontsize = 8, legend = :bottom);
display(p_frac_lys_kd_he)

savefig(p_frac_lys_kd_he, @projectroot("deliv/figure/payload/sens-frac-lys-he.png"));

# calculate Cavg 
using Statistics
mean(@rsubset(sol_tdm1, :frac_lys ==0.1).endo_he)*1E3
mean(@rsubset(sol_tdm1, :frac_lys ==0.3).endo_he)*1E3
mean(@rsubset(sol_tdm1, :frac_lys ==0.5).endo_he)*1E3
mean(@rsubset(sol_tdm1, :frac_lys ==0.7).endo_he)*1E3
mean(@rsubset(sol_tdm1, :frac_lys ==0.9).endo_he)*1E3
mean(@rsubset(sol_tdm1, :frac_lys ==1).endo_he)*1E3

# calculate Cmax 
maximum(@rsubset(sol_tdm1, :frac_lys ==0.1).endo_he)*1E3
maximum(@rsubset(sol_tdm1, :frac_lys ==0.3).endo_he)*1E3
maximum(@rsubset(sol_tdm1, :frac_lys ==0.5).endo_he)*1E3
maximum(@rsubset(sol_tdm1, :frac_lys ==0.7).endo_he)*1E3
maximum(@rsubset(sol_tdm1, :frac_lys ==0.9).endo_he)*1E3
maximum(@rsubset(sol_tdm1, :frac_lys ==1).endo_he)*1E3