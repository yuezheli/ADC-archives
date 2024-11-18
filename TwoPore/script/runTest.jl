# author: Yuezhe Li 
# date: Nov 16, 2023
# purpose of this code: to test run the model 

using Pkg; Pkg.activate("..");
# Pkg.instantiate(); # run this line once when setting up the repo for the first time 

using DifferentialEquations 
using Plots
using DataFrames, CSV

# read in observed data
obs_nanobody = CSV.read("../data/nanobody_Gainkam.csv",DataFrame);
obs2_scFv = CSV.read("../data/scFv_Milenic.csv",DataFrame);
obs3_scFv = CSV.read("../data/scFv_Wu.csv",DataFrame);
obs_Fab = CSV.read("../data/Fab.csv",DataFrame);
obs2_Fab = CSV.read("../data/Fab_Otsuji.csv",DataFrame);
obs3_Fab = CSV.read("../data/Fab_Milenic.csv",DataFrame);
obs_scFv2 = CSV.read("../data/scFv2.csv",DataFrame);
obs2_scFv2 = CSV.read("../data/scFv2_Pavlinkova.csv",DataFrame);
obs_minibody = CSV.read("../data/minibody.csv",DataFrame);
obs_Fab2 = CSV.read("../data/Fab2.csv",DataFrame);
obs_IgG = CSV.read("../data/IgG.csv",DataFrame);

# simulation 
include("helper.jl")
include("../model/twopore_mus.jl");

tspan = (0., 150.);

# simulation for nanobody (13kDa)
u0_nanobody, dose_nano = lishah_init(15, 0.1, 0.02, 0.944/1000, 13.0E3);
sol_nano = solve(ODEProblem(lishah_mus!,u0_nanobody, (0., 2.),p_nanobody), alg=AutoTsit5(Rosenbrock23(autodiff = false)), saveat = 0.01);
cplasma_nano = [sol_nano.u[i].C_EXG_Plasma for i in 1:length(sol_nano.t)]; # [nM]

# simulation for scFv (27kDa)
u0_scFv, dose_scFv = lishah_init(15, 1., 0.02, 0.944/1000, 27.0E3);
sol_scFv = solve(ODEProblem(lishah_mus!,u0_scFv,tspan,p_scFv), alg=AutoTsit5(Rosenbrock23(autodiff = false)), saveat = 0.02);
cplasma_scFv = [sol_scFv.u[i].C_EXG_Plasma for i in 1:length(sol_scFv.t)]; # [nM]

# simulation for Fab (50kDa)
u0_Fab, dose_Fab = lishah_init(15, 1., 0.02, 0.944/1000, 50.0E3);
sol_Fab = solve(ODEProblem(lishah_mus!,u0_Fab,tspan,p_Fab), alg=AutoTsit5(Rosenbrock23(autodiff = false)), saveat = 0.1);
cplasma_Fab = [sol_Fab.u[i].C_EXG_Plasma for i in 1:length(sol_Fab.t)]; # [nM]

# simulation for scFv2 (55kDa)
u0_scFv2, dose_scFv2 = lishah_init(15, 1., 0.02, 0.944/1000, 55.0E3);
sol_scFv2 = solve(ODEProblem(lishah_mus!,u0_scFv2,tspan,p_scFv2), alg=AutoTsit5(Rosenbrock23(autodiff = false)), saveat = 0.1);
cplasma_scFv2 = [sol_scFv2.u[i].C_EXG_Plasma for i in 1:length(sol_scFv2.t)]; # [nM]

# simulation for minibody (80kDa)
u0_minibody, dose_minibody = lishah_init(15, 1., 0.02, 0.944/1000, 80.0E3);
sol_minibody = solve(ODEProblem(lishah_mus!,u0_minibody,tspan,p_minibody), alg=AutoTsit5(Rosenbrock23(autodiff = false)), saveat = 1.);
cplasma_minibody = [sol_minibody.u[i].C_EXG_Plasma for i in 1:length(sol_minibody.t)]; # [nM]

# simulation for Fab2 (100kDa)
u0_Fab2, dose_Fab2 = lishah_init(15, 1., 0.02, 0.944/1000, 100.0E3);
sol_Fab2 = solve(ODEProblem(lishah_mus!,u0_Fab2,tspan,p_Fab2), alg=AutoTsit5(Rosenbrock23(autodiff = false)), saveat = 1.);
cplasma_Fab2 = [sol_Fab2.u[i].C_EXG_Plasma for i in 1:length(sol_Fab2.t)]; # [nM]

# simulation for IgG (150kDa)
u0_IgG, dose_IgG = lishah_init(15, 1., 0.02, 0.944/1000, 150.0E3);
sol_IgG = solve(ODEProblem(lishah_mus!,u0_IgG,tspan,p_IgG), alg=AutoTsit5(Rosenbrock23(autodiff = false)), saveat = 1.);
cplasma_IgG = [sol_IgG.u[i].C_EXG_Plasma for i in 1:length(sol_IgG.t)]; # [nM]

# plot data 
plot_nanobody = plot(title = "nanobody 13kDa", titlefontsize = 8, xlabel = "time (min)", ylabel = "%IA/TBV", xlims = (0., 70)); 
plot!(sol_nano.t*60, cplasma_nano*V_Plasma/dose_nano, label = "sims"); 
plot!(obs_nanobody[:,1], obs_nanobody[:,2]/100, seriestype=:scatter, label="Gainkam et al., 2008");

plot_scFv = plot(title = "scFv 27kDa", titlefontsize = 8, xlabel = "time (hr)", ylabel = "% plasma", xlims = (0, 3)); 
plot!(sol_scFv.t, cplasma_scFv*V_Plasma/dose_scFv*100, label = "sims"); 
plot!(obs2_scFv[:,1]/60, obs2_scFv[:,2], seriestype=:scatter, label="Milenic et al., 1999");
plot!(obs3_scFv[:,1], obs3_scFv[:,2], seriestype=:scatter, label="Wu et al., 1996");

plot_Fab = plot(title = "Fab 50kDa", titlefontsize = 8, xlabel = "time (hr)", ylabel = "normalized conc (L-1)", ylims = (0.01, 1E4)); 
plot!(sol_Fab.t, cplasma_Fab/dose_Fab, label = "sims"); 
plot!(obs_Fab[:,1], obs_Fab[:,2], seriestype=:scatter, label="Li & Shah, 2019");
plot!(yaxis = :log);

plot2_Fab = plot(title = "Fab 50kDa", titlefontsize = 8, xlabel = "time (hr)", ylabel = "% injection dose/ g blood", xlims = (0., 24)); 
plot!(sol_Fab.t, cplasma_Fab/dose_Fab/1000*100, label = "sims"); 
plot!(obs2_Fab[:,1], obs2_Fab[:,2], seriestype=:scatter, label="Otsuji et al., 2003");
plot!(obs3_Fab[:,1]/60, obs3_Fab[:,2]*dose_Fab/(V_Plasma*1000), seriestype=:scatter, label="Milenic et al., 1999");

plot_scFv2 = plot(title = "scFv2 55kDa", titlefontsize = 8, xlabel = "time (hr)", ylabel = "normalized conc (L-1)", ylims = (0.01, 1E4)); 
plot!(sol_scFv2.t, cplasma_scFv2/dose_scFv2, label = "sims"); 
plot!(obs_scFv2[:,1], obs_scFv2[:,2], seriestype=:scatter, label="Li & Shah, 2019");
plot!(yaxis = :log); 

plot2_scFv2 = plot(title = "scFv2 55kDa", titlefontsize = 8, xlabel = "time (hr)", ylabel = "% injected dose", xlims = (0, 30)); 
plot!(sol_scFv2.t, cplasma_scFv2*V_Plasma/dose_scFv2*100, label = "sims"); 
plot!(obs2_scFv2[:,1]/60, obs2_scFv2[:,2], seriestype=:scatter, label="Pavlinkova et al., 1999");

plot_minibody = plot(title = "minibody 80kDa", titlefontsize = 8, xlabel = "time (hr)", ylabel = "normalized conc (L-1)", ylims = (0.01, 1E4)); 
plot!(sol_minibody.t, cplasma_minibody/dose_minibody, label = "sims"); 
plot!(obs_minibody[:,1], obs_minibody[:,2], seriestype=:scatter, label="Li & Shah, 2019");
plot!(yaxis = :log); 

plot_Fab2 = plot(title = "Fab2 100kDa", titlefontsize = 8, xlabel = "time (hr)", ylabel = "normalized conc (L-1)", ylims = (0.01, 1E4)); 
plot!(sol_Fab2.t, cplasma_Fab2/dose_Fab2, label = "sims"); 
plot!(obs_Fab2[:,1], obs_Fab2[:,2], seriestype=:scatter, label="Li & Shah, 2019");
plot!(yaxis = :log); 

plot_IgG = plot(title = "IgG 150kDa", titlefontsize = 8, xlabel = "time (hr)", ylabel = "normalized conc (L-1)", ylims = (0.01, 1E4)); 
plot!(sol_IgG.t, cplasma_IgG/dose_IgG, label = "sims"); 
plot!(obs_IgG[:,1], obs_IgG[:,2], seriestype=:scatter, label="Li & Shah, 2019");
plot!(yaxis = :log); 

# save figures 
savefig(plot_nanobody, "../deliv/figure/nanobody-sims.png");
savefig(plot_scFv, "../deliv/figure/scFv-sims.png");
savefig(plot_Fab, "../deliv/figure/Fab-sims1.png");
savefig(plot2_Fab, "../deliv/figure/Fab-sims2.png");
savefig(plot_minibody, "../deliv/figure/minibody-sims.png");
savefig(plot_Fab2, "../deliv/figure/Fab2-sims.png");
savefig(plot_IgG, "../deliv/figure/IgG-sims.png");

