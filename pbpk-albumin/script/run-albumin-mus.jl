# date: 2/22/2026 
# author: Yuezhe Li 
# purpose of this code: to test run the two-pore model for albumin that was defined in Liu et al., J Pharmacokinet. Pharmacodyn, 2024
# https://pubmed.ncbi.nlm.nih.gov/38691205/

using ProjectRoot
using Pkg; Pkg.activate(@projectroot);

using ModelingToolkit
using ModelingToolkit: getdefault
using DifferentialEquations
using DataFrames
using Plots
using Plots.PlotMeasures
using SymbolicIndexingInterface
using CSV
using Sundials

# observed PK 
plasma_data = DataFrame(
    time_hr = [0.0, 1.27, 1.27, 2.54, 3.81, 11.44, 11.44, 24.15, 24.15, 24.15, 48.31, 48.31, 47.03, 72.46, 72.46, 71.19, 71.19, 71.19, 96.61, 96.61, 96.61, 95.34, 120.76, 120.76, 120.76, 120.76, 143.64, 143.64, 144.92, 144.92, 144.92],
    conc_nM = [288.85, 211.76, 155.24, 126.22, 92.53, 92.53, 75.23, 75.23, 61.17, 49.73, 55.15, 36.46, 29.64, 32.87, 24.1, 21.73, 19.59, 14.36, 19.59, 15.93, 12.95, 11.68, 15.93, 12.95, 11.68, 9.5, 8.56, 6.96, 5.66, 4.15, 3.04]
); 

# load constants 
include(@projectroot("script/constants.jl"));

# load model 
include(@projectroot("model/twopore_albumin.jl"));
@named pbpk = albumin_pbpk_mtk();
pbpk_sys = structural_simplify(pbpk);

# Set the time span (e.g., 168 hours for 1 week)
tspan = (0.0, 168.0);

# Create and Solve the ODEProblem
# The parameters (pars) are already baked into the system defaults.
prob = ODEProblem(pbpk_sys, [pbpk_sys.C_Plasma_Exo => 1000E-9], tspan);

# Using Rodas5 for stiff PBPK systems is generally recommended for stability
sol = solve(prob, Rodas5(), reltol=1e-6, abstol=1e-9)

# 5. Accessing results
plt_plasma = plot(xlabel="Time (hr)", ylabel="Conc (nM)", yaxis = :log10, yticks = [1E-1, 1, 10, 1E2, 1E3], ylims = [1E-1, 1E3]);
plot!(sol.t, sol[pbpk.C_Plasma_Exo] * 1E9, label = "sims");
plot!(plasma_data.time_hr, plasma_data.conc_nM, seriestype = :scatter, label = "Liu et al., 2024"); 
display(plt_plasma)

savefig(plt_plasma, @projectroot("deliv/figure/plasma-pk-albumin-mus.png"));

