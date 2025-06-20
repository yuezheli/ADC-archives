# date: 2/11/2025
# author: Yuezhe Li 
# purpose of this code: sensitivity analysis on DAR, payload diffusion rate, and membrane affinity 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)

using DifferentialEquations, ComponentArrays
using DataFrames
using Statistics
using DataFramesMeta
using CSV

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

PS_kd_values = 10 .^ range(-3.0, 3, step=0.5);
DAR_values = [1,2,4,6,8];
k_out_values = 10 .^ range(-3.0, 3, step=0.2);

conc = DataFrame(PS_kd = [], DAR = [], k_out = [], cmax = [], cavg = [], time = [])

for i in 1:length(PS_kd_values)
    for j in 1:length(DAR_values)
        for k in 1:length(k_out_values)
            p_tmp = deepcopy(p_0);
            p_tmp.PS_kd = PS_kd_values[i];
            p_tmp.DAR = DAR_values[j];
            p_tmp.k_out = k_out_values[k];
            prob_tmp = remake(prob_tdm1, p = p_tmp);
            sol_tmp = tissue_ints_endo(solve(prob_tmp, saveat = 1., alg = QNDF(autodiff=false), reltol = 1E-18));
            cmax_tmp = maximum(sol_tmp.endo_he)*1E3  # [nM]
            cavg_tmp = mean(sol_tmp.endo_he)*1E3  # [nM]
            temp_conc = DataFrame(PS_kd = PS_kd_values[i], DAR = DAR_values[j], k_out = k_out_values[k], cmax = cmax_tmp, cavg = cavg_tmp, time = sol_tmp.time[end])
            conc = vcat(conc, temp_conc)
        end
    end
end

# save data 
CSV.write("data/sensitivity/dar-diffusion-membrane-affinity.csv", conc);

# visualize data 
conc = CSV.read("data/sensitivity/dar-diffusion-membrane-affinity.csv", DataFrame);

using Plots
using LaTeXStrings

cmax_dar_1 = zeros(length(PS_kd_values), length(k_out_values));
cmax_dar_8 = zeros(length(PS_kd_values), length(k_out_values));

cavg_dar_1 = zeros(length(PS_kd_values), length(k_out_values));
cavg_dar_8 = zeros(length(PS_kd_values), length(k_out_values));

for i in 1:length(PS_kd_values)
    for k in 1:length(k_out_values)
        cmax_dar_1[i,k] = @rsubset(conc, :DAR==1, :k_out == k_out_values[k], :PS_kd == PS_kd_values[i]).cmax[1]
        cmax_dar_8[i,k] = @rsubset(conc, :DAR==8, :k_out == k_out_values[k], :PS_kd == PS_kd_values[i]).cmax[1]
        cavg_dar_1[i,k] = @rsubset(conc, :DAR==1, :k_out == k_out_values[k], :PS_kd == PS_kd_values[i]).cavg[1]
        cavg_dar_8[i,k] = @rsubset(conc, :DAR==8, :k_out == k_out_values[k], :PS_kd == PS_kd_values[i]).cavg[1]
    end
end

p_cmax_dar_1 = plot(dpi = 1000, size = (400, 400), axis = :log);
contour!(k_out_values,PS_kd_values, log10.(cmax_dar_1), fill=true, levels=20, color=:turbo, clim=(-4,2), lw=0); 
plot!(xlabel = "Payload diffusion (1/hr)", ylabel = L"K$_D$, membrane ($\mu$M)", cbar_title=L"log(C$_{max}$) (nM)");
display(p_cmax_dar_1)

p_cmax_dar_8 = plot(dpi = 1000, size = (400, 400), axis = :log);
contour!(k_out_values,PS_kd_values, log10.(cmax_dar_8), fill=true, levels=20, color=:turbo, clim=(-4,2), lw=0); 
plot!(xlabel = "Payload diffusion (1/hr)", ylabel = L"K$_D$, membrane ($\mu$M)", cbar_title=L"log(C$_{max}$) (nM)");
display(p_cmax_dar_8)

p_cavg_dar_1 = plot(dpi = 1000, size = (400, 400), axis = :log);
contour!(k_out_values,PS_kd_values, log10.(cavg_dar_1), fill=true, levels=20, color=:turbo, clim=(-4,2), lw=0); 
plot!(xlabel = "Payload diffusion (1/hr)", ylabel = L"K$_D$, membrane ($\mu$M)", cbar_title=L"log(C$_{avg}$) (nM)");
display(p_cavg_dar_1)

p_cavg_dar_8 = plot(dpi = 1000, size = (400, 400), axis = :log);
contour!(k_out_values,PS_kd_values, log10.(cavg_dar_8), fill=true, levels=20, color=:turbo, clim=(-4,2), lw=0); 
plot!(xlabel = "Payload diffusion (1/hr)", ylabel = L"K$_D$, membrane ($\mu$M)", cbar_title=L"log(C$_{avg}$) (nM)");
display(p_cavg_dar_8)

savefig(p_cavg_dar_1, @projectroot("deliv/figure/payload/sens-membrane-kd-diffusion-dar-1-cavg.png"));
savefig(p_cavg_dar_8, @projectroot("deliv/figure/payload/sens-membrane-kd-diffusion-dar-8-cavg.png"));
savefig(p_cmax_dar_1, @projectroot("deliv/figure/payload/sens-membrane-kd-diffusion-dar-1-cmax.png"));
savefig(p_cmax_dar_8, @projectroot("deliv/figure/payload/sens-membrane-kd-diffusion-dar-8-cmax.png"));
