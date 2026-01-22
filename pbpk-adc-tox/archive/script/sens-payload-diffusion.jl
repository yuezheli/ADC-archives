# date: 7/3/2025 
# author: Yuezhe Li 
# purpose of this code: to conduct sensitivity analysis for payload diffusion rate 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Parameters: @unpack
using DataFrames
using DataFramesMeta
using Statistics
using Plots
using LaTeXStrings
using CSV

include(@projectroot("Constants.jl"))
include(@projectroot("model/jones_homo.jl") )
include(@projectroot("model/param-pk.jl"))
include(@projectroot("model/init-pk.jl")) 

p_tmd1 = deepcopy(p_base);
p_tmd1.CL_PL_plasma = log(2)/3; # https://pubmed.ncbi.nlm.nih.gov/37787918/

first_dose = 3.6  # [mg/kg]
tspan = (0, hr_per_day*21);      # [hr]

u0_tdm1 = jones_init(first_dose*1E3, p_tmd1, BW, V_Plasma);  # IV bolus dose
prob_0 = ODEProblem(jonesODEs_homo_tumor!, u0_tdm1, tspan, p_tmd1);
sol0 = solve(prob_0, saveat = 1., alg = QNDF(autodiff=false)); 

soldf = DataFrame(
    Cmax = maximum( [sol0.u[i].end_cyto_payload[2] for i in 1:length(sol0.t)] ), 
    Cavg = mean( [sol0.u[i].end_cyto_payload[2] for i in 1:length(sol0.t)] ), 
    k_out = p_tmd1.k_out
)

# below using T-DM1 PK as an example; see what level of payload diffusivity would result in same concentration in tissue endothelial cells and interstitium
k_out_values = [collect(0.1:0.1:3); 4.32; 5; 10; 15; 20; 25; 30; 32.32; 40];

for tmp_k_out in k_out_values
    p_tmp = deepcopy(p_tmd1);
    p_tmp.k_out = tmp_k_out;
    prob_tmp = remake(prob_0, p = p_tmp);
    sol_tmp = solve(prob_tmp, saveat = 1., alg = QNDF(autodiff=false)); 
    tmp_soldf = DataFrame(
            Cmax = maximum( [sol_tmp.u[i].end_cyto_payload[2] for i in 1:length(sol_tmp.t)] ), 
            Cavg = mean( [sol_tmp.u[i].end_cyto_payload[2] for i in 1:length(sol_tmp.t)] ), 
            k_out = p_tmp.k_out
        );
    soldf = vcat(soldf, tmp_soldf)
end


# Cmax
p_payload_cmax = plot(#ylims = [1, 20], yticks = [1, 10], yaxis = :log, 
                        dpi = 300, size = (400, 400));
plot!(soldf.k_out, soldf.Cmax*1E3, seriestype=:scatter, alpha = 0.3, label = false, color = :grey);
scatter!([0.14, 2, 32.32], 
        [@rsubset(soldf, :k_out==0.14).Cmax[1]*1E3, 
         @rsubset(soldf, :k_out==2).Cmax[1]*1E3, 
         @rsubset(soldf, :k_out==32.32).Cmax[1]*1E3], alpha = 1, label = false, color = :blue);
annotate!([(0.14, @rsubset(soldf, :k_out==0.14).Cmax[1]*1E3, ("DM1 \n", 10, :left, :blue))]);
annotate!([(2, @rsubset(soldf, :k_out==2).Cmax[1]*1E3, ("MMAE \n", 10, :left, :blue))]);
# annotate!([(4.32, @rsubset(soldf, :k_out==4.32).Cmax[1]*1E3, ("PBD \n", 10, :left, :blue))]);
annotate!([(32.32, @rsubset(soldf, :k_out==32.32).Cmax[1]*1E3, ("Dxd \n", 10, :left, :blue))]);
plot!(ylabel = "Maximum payload concentration (nM)", xlabel = "Payload diffusion rate (1/hr)");
display(p_payload_cmax)

savefig(p_payload_cmax, @projectroot("deliv/figure/pl/diffusion-sensitivity-analysis.png"));
