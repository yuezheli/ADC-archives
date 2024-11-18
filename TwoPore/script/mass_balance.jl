# author: Yuezhe Li 
# date: Nov 17, 2023
# purpose of this script: check mass balance 

using Pkg; Pkg.activate("..");
using DifferentialEquations 
using Plots
using DataFrames, CSV

# simulation 
include("helper.jl")
include("../model/twopore_mus.jl");

tspan = (0., 150.);

p = deepcopy(p_nanobody);

# simulation for nanobody, as an example
u0, dose = lishah_init(15, 1., 0.02, 0.944/1000, 13.0E3);
sol = solve(ODEProblem(lishah_mus!,u0, tspan,p), alg=AutoTsit5(Rosenbrock23(autodiff = false)), saveat = 1.);

system_protein = [];
total_protein = [];

for i in 1:length(sol.t)
    organvasprotein = sum(sol.u[i].C_EXG[:, 1] .*  V_V);
    organendoprotein = sum(sol.u[i].C_EXG[:, 2] .*  V_endosomal);
    organintsprotein = sum(sol.u[i].C_EXG[:, 3] .*  V_IntS);
    tmp_sys_protein = organvasprotein + organendoprotein + organintsprotein + sol.u[i].C_EXG_Plasma*V_Plasma + sol.u[i].C_EXG_LN*V_LN;
    tmp_tot_protein = tmp_sys_protein + sol.u[i].endo_deg + sol.u[i].renal_clearance; 
    append!(system_protein, tmp_sys_protein);
    append!(total_protein, tmp_tot_protein)
end

println(total_protein/dose)

