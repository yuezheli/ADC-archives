# author: Yuezhe Li 
# date: 3/1/24
# purpose of this code: to compare the simulation outcome between the Jones model and reduced Jones model (i.e. mPBPK)
# to run this script, set workinng directory in the home directory (check using `pwd()`)

using Pkg; Pkg.activate("");

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Plots, Statistics, DataFrames, CSV, DataFramesMeta

include("model/jones_tumor_homo.jl");
include("model/mPBPK_tumor_homo.jl");
include("script/init.jl");
include("script/param.jl");

TotalCells = 1E-3/Vc

# update params (default params, for T-DM1)
p_base.Rcopies = 1.0E6
p_base.PS_Score = 6.
p_base.init_sR = 0.004
p_base.thalf_sR_adc = 120.

tspan = (0., DayToHour*42);    # [hr]

u0_jones, _ = pbpk_init(15, TotalCells, 0., 3.6E3, p_base.Rcopies, p_base.init_sR, p_base.k_endo, p_base.k_rec); 
u0_mpbpk, _ = pbpk_init(7, TotalCells, 0., 3.6E3, p_base.Rcopies, p_base.init_sR, p_base.k_endo, p_base.k_rec); 

sol_jones = solve(ODEProblem(jonesODEs_homo_tumor!, u0_jones, tspan, p_base), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-12);
sol_mpbpk = solve(ODEProblem(mPBPK_homo_tumor!, u0_mpbpk, tspan, p_base), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-12);

plot_plasma =  plot(ylabel = "Plasma conc (uM)", xlabel = "Time (h)", legend = :bottomleft, yaxis = :log);
plot!(sol_jones.t, sol_jones[:C_EXG_Plasma], label = "sims from Jones et al., 2019", linestyle = :solid, linewidth=2, alpha = 0.8);
plot!(sol_mpbpk.t, sol_mpbpk[:C_EXG_Plasma], label = "sims from mPBPK model", linestyle = :dash, linewidth=3, alpha = 0.8);

plot_ints_skin =  plot(ylabel = "Skin ints conc (uM)", xlabel = "Time (h)", legend = :bottomright, yaxis = :log);
plot!(sol_jones.t, [sol_jones.u[i].C_EXG[5, 7] for i in 1:length(sol_jones)], label = "sims from Jones et al., 2019", linestyle = :solid, linewidth=2, alpha = 0.8);
plot!(sol_mpbpk.t, [sol_mpbpk.u[i].C_EXG[3, 7] for i in 1:length(sol_mpbpk)], label = "sims from mPBPK model", linestyle = :dash, linewidth=3, alpha = 0.8);

plot_ints_liver =  plot(ylabel = "Liver ints conc (uM)", xlabel = "Time (h)", legend = :bottomright, yaxis = :log);
plot!(sol_jones.t, [sol_jones.u[i].C_EXG[2, 7] for i in 1:length(sol_jones)], label = "sims from Jones et al., 2019", linestyle = :solid, linewidth=2, alpha = 0.8);
plot!(sol_mpbpk.t, [sol_mpbpk.u[i].C_EXG[2, 7] for i in 1:length(sol_mpbpk)], label = "sims from mPBPK model", linestyle = :dash, linewidth=3, alpha = 0.8);

plot_ints_kidney =  plot(ylabel = "Kidney ints conc (uM)", xlabel = "Time (h)", legend = :bottomright, yaxis = :log);
plot!(sol_jones.t, [sol_jones.u[i].C_EXG[9, 7] for i in 1:length(sol_jones)], label = "sims from Jones et al., 2019", linestyle = :solid, linewidth=2, alpha = 0.8);
plot!(sol_mpbpk.t, [sol_mpbpk.u[i].C_EXG[4, 7] for i in 1:length(sol_mpbpk)], label = "sims from mPBPK model", linestyle = :dash, linewidth=3, alpha = 0.8);

plot_combined = plot(plot_plasma, plot_ints_skin, plot_ints_liver, plot_ints_kidney, size = (800, 500))

savefig(plot_combined, "deliv/figure/mpbpk_jones_comparison.png");

