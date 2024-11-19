# author: Yuezhe Li 
# date: 3/14/24
# purpose of this code: to compare the simulation outcome between the mPBPK model, with or without blood cell compartment
# to run this script, set workinng directory in the home directory (check using `pwd()`)

using Pkg; Pkg.activate("");

using DifferentialEquations, ComponentArrays
using Plots, Statistics, DataFrames, CSV, DataFramesMeta

include("mPBPK_bc_tumor_homo.jl");
include("mPBPK_tumor_homo.jl");
include("init.jl");
include("param.jl");

TotalCells = 1E-3/Vc

# update params (default params, for T-DM1)
p_base.Rcopies = 1.0E6
p_base.PS_Score = 6.
p_base.init_sR = 0.004
p_base.thalf_sR_adc = 120.

tspan = (0., DayToHour*42);    # [hr]

u0_mpbpk, Dose_umol = pbpk_init(7, TotalCells, 0., 3.6E3, p_base.Rcopies, p_base.init_sR, p_base.k_endo, p_base.k_rec); 
u0_bc_mpbpk = ComponentArray(u0_mpbpk, C_EXG_Plasma_BC = 0., C_BC_EXG = zeros(7));
# ADC 
u0_bc_mpbpk.C_EXG_Plasma = Dose_umol / V_Plasma
u0_bc_mpbpk.C_EXG_Plasma_BC = 0 / V_BC

sol_mpbpk = solve(ODEProblem(mPBPK_homo_tumor!, u0_mpbpk, tspan, p_base), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-12);
sol_bc_mpbpk = solve(ODEProblem(mPBPK_bc_homo_tumor!, u0_bc_mpbpk, tspan, p_base), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-12);

plot_plasma =  plot(ylabel = "Plasma conc (uM)", xlabel = "Time (h)", legend = :bottomleft, yaxis = :log);
plot!(sol_bc_mpbpk.t, sol_bc_mpbpk[:C_EXG_Plasma], label = "sims from model with blood cells", linestyle = :solid, linewidth=2, alpha = 0.8);
plot!(sol_mpbpk.t, sol_mpbpk[:C_EXG_Plasma], label = "sims from mPBPK model w/o blood cells", linestyle = :dash, linewidth=3, alpha = 0.8);

plot_ints_skin =  plot(ylabel = "Skin ints conc (uM)", xlabel = "Time (h)", legend = :bottomright, yaxis = :log);
plot!(sol_bc_mpbpk.t, [sol_bc_mpbpk.u[i].C_EXG[5, 7] for i in 1:length(sol_bc_mpbpk)], label = "sims from model with blood cells", linestyle = :solid, linewidth=2, alpha = 0.8);
plot!(sol_mpbpk.t, [sol_mpbpk.u[i].C_EXG[3, 7] for i in 1:length(sol_mpbpk)], label = "sims from mPBPK model w/o blood cells", linestyle = :dash, linewidth=3, alpha = 0.8);

plot_ints_liver =  plot(ylabel = "Liver ints conc (uM)", xlabel = "Time (h)", legend = :bottomright, yaxis = :log);
plot!(sol_bc_mpbpk.t, [sol_bc_mpbpk.u[i].C_EXG[2, 7] for i in 1:length(sol_bc_mpbpk)], label = "sims from model with blood cells", linestyle = :solid, linewidth=2, alpha = 0.8);
plot!(sol_mpbpk.t, [sol_mpbpk.u[i].C_EXG[2, 7] for i in 1:length(sol_mpbpk)], label = "sims from mPBPK model w/o blood cells", linestyle = :dash, linewidth=3, alpha = 0.8);

plot_ints_kidney =  plot(ylabel = "Kidney ints conc (uM)", xlabel = "Time (h)", legend = :bottomright, yaxis = :log);
plot!(sol_bc_mpbpk.t, [sol_bc_mpbpk.u[i].C_EXG[4, 7] for i in 1:length(sol_bc_mpbpk)], label = "sims from model with blood cells", linestyle = :solid, linewidth=2, alpha = 0.8);
plot!(sol_mpbpk.t, [sol_mpbpk.u[i].C_EXG[4, 7] for i in 1:length(sol_mpbpk)], label = "sims from mPBPK model w/o blood cells", linestyle = :dash, linewidth=3, alpha = 0.8);

plot_combined = plot(plot_plasma, plot_ints_skin, plot_ints_liver, plot_ints_kidney, size = (800, 500))

savefig(plot_combined, "mpbpk_bc_comparison.png");

