# date: 2/4/2025
# author: Yuezhe Li 
# purpose of this code: sensitivity analysis on the payload diffusion rate

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)

using DifferentialEquations, ComponentArrays
using Plots, DataFrames, LaTeXStrings
using Statistics

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

sol_tdm1 = tissue_ints_endo(solve(prob_tdm1, saveat = 1., alg = QNDF(autodiff=false)));
sol_tdm1.Camx_he .= maximum(sol_tdm1.endo_he);
sol_tdm1.Cavg_he .= mean(sol_tdm1.endo_he);
sol_tdm1.k_out .= 0.14;

# below using T-DM1 PK as an example; see what level of payload diffusivity would result in same concentration in tissue endothelial cells and interstitium
k_out_values = [collect(0.1:0.1:3); 2.27; 5; 10; 15; 20; 24; 32.32; 40];

for tmp_k_out in k_out_values
    p_tmp = deepcopy(p_tdm1);
    p_tmp.k_out = tmp_k_out;
    prob_tmp = remake(prob_tdm1, p = p_tmp);
    sol_tmp = tissue_ints_endo(solve(prob_tmp, saveat = 1., alg = QNDF(autodiff=false)));
    sol_tmp.Camx_he .= maximum(sol_tmp.endo_he);
    sol_tmp.Cavg_he .= mean(sol_tmp.endo_he);
    sol_tmp.k_out .= tmp_k_out;
    sol_tdm1 = vcat(sol_tdm1, sol_tmp);
end

@rsubset!(sol_tdm1, :time > 0);

final_pl = @rsubset(sol_tdm1, :time == DayToHour*21)


# Cmax
p_payload_cmax = plot(ylims = [0.01, 20], yticks = [0.01, 0.1, 1, 10], yaxis = :log, dpi = 1000, size = (400, 400));
plot!(final_pl.k_out, final_pl.Camx_he*1E3, seriestype=:scatter, alpha = 0.3, label = false, color = :grey);
scatter!([0.14, 2.27, 24, 32.32], [@rsubset(final_pl, :k_out==0.14).Camx_he[1]*1E3, @rsubset(final_pl, :k_out==2.27).Camx_he[1]*1E3, @rsubset(final_pl, :k_out==24).Camx_he[1]*1E3, @rsubset(final_pl, :k_out==32.32).Camx_he[1]*1E3], alpha = 1, label = false, color = :blue);
annotate!([(0.14, @rsubset(final_pl, :k_out==0.14).Camx_he[1]*1E3, ("DM1 \n", 10, :left, :blue))]);
annotate!([(2.27, @rsubset(final_pl, :k_out==2.27).Camx_he[1]*1E3, ("DM4 \n", 10, :left, :blue))]);
annotate!([(24, @rsubset(final_pl, :k_out==24).Camx_he[1]*1E3, ("PBD \n", 10, :left, :blue))]);
annotate!([(32.32, @rsubset(final_pl, :k_out==32.32).Camx_he[1]*1E3, ("Dxd \n", 10, :left, :blue))]);
plot!(ylabel = "Maximum payload concentration (nM)", xlabel = "payload diffusion rate (1/hr)");
display(p_payload_cmax)

savefig(p_payload_cmax, @projectroot("deliv/figure/payload/sens-pl-diffusion-Cmax.png"));

# Cavg
p_payload_cavg = plot(ylims = [0.01, 20], yticks = [0.01, 0.1, 1, 10], yaxis = :log, dpi = 1000, size = (400, 400));
plot!(final_pl.k_out, final_pl.Cavg_he*1E3, seriestype=:scatter, alpha = 0.3, label = false, color = :grey);
scatter!([0.14, 2.27, 32.32], [@rsubset(final_pl, :k_out==0.14).Cavg_he[1]*1E3, @rsubset(final_pl, :k_out==2.27).Cavg_he[1]*1E3, @rsubset(final_pl, :k_out==32.32).Cavg_he[1]*1E3], alpha = 1, label = false, color = :blue);
scatter!([24], [@rsubset(final_pl, :k_out==24).Cavg_he[1]*1E3], alpha = 1, label = false, color = :blue);
annotate!([(0.14, @rsubset(final_pl, :k_out==0.14).Cavg_he[1]*1E3, ("DM1 \n", 10, :left, :blue))]);
annotate!([(2.27, @rsubset(final_pl, :k_out==2.27).Cavg_he[1]*1E3, ("DM4 \n", 10, :left, :blue))]);
annotate!([(24, @rsubset(final_pl, :k_out==24).Cavg_he[1]*1E3, ("PBD \n", 10, :left, :blue))]);
annotate!([(32.32, @rsubset(final_pl, :k_out==32.32).Cavg_he[1]*1E3, ("Dxd \n", 10, :left, :blue))]);
plot!(ylabel = "Average payload concentration (nM)", xlabel = "payload diffusion rate (1/hr)");
display(p_payload_cavg)

savefig(p_payload_cavg, @projectroot("deliv/figure/payload/sens-pl-diffusion-Cavg.png"));


# payload equilibration

p_ratio = plot(xaxis = :log, xticks = ([0.1, 1, 2, 5, 10, 15, 20, 40], ["0.1", "1", "2", "5", "10", "15", "20", "40"]), dpi = 1000, size = (400, 400));
plot!(final_pl.k_out, final_pl.r_he, seriestype=:scatter, alpha = 0.5, label = false, color = :grey);
plot!(xlabel = L"payload diffusion rate (hr$^{-1}$)", guidefontsize = 10, xguidefont = font("Computer Modern"));
plot!(ylabel = L"$\frac{liver \; endothelial \; cytosol \; payload \; concentration}{liver \; interstitium \; payload \; concentration}$");
plot!([0.14, 2.27, 32.32], [@rsubset(final_pl, :k_out==0.14).r_he[1], @rsubset(final_pl, :k_out==2.27).r_he[1], @rsubset(final_pl, :k_out==32.32).r_he[1]], seriestype=:scatter, alpha = 1, label = false, color = :blue);
scatter!([24], [@rsubset(final_pl, :k_out==24).r_he[1]], alpha = 1, label = false, color = :blue);
annotate!([(0.14, @rsubset(final_pl, :k_out==0.14).r_he[1], ("DM1 \n", 10, :bottom, :blue))]);
annotate!([(2.27, @rsubset(final_pl, :k_out==2.27).r_he[1], ("DM4 \n", 10, :bottom, :blue))]);
annotate!([(24, @rsubset(final_pl, :k_out==24).r_he[1], ("PBD ", 10, :right, :blue))]);
annotate!([(32.32, @rsubset(final_pl, :k_out==32.32).r_he[1], ("Dxd \n", 10, :bottom, :blue))]);
plot!(xrotation = 45);
display(p_ratio)

savefig(p_ratio, @projectroot("deliv/figure/payload/sens-pl-diffusion.png"));
