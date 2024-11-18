# date: 4/2/24
# author: Yuezhe Li 
# purpose of this script: to test a full pbpk model from He et al., 2019; 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6533104/

using Pkg; Pkg.activate("");

using DifferentialEquations, ComponentArrays
using Plots, DataFrames, CSV, DataFramesMeta

include("model/dox_fullpbpk.jl");

# define params 
p_dox = ComponentArray(
    infusion = 0., 
    Kpt_blood = 0.68, 
    Va = 3.126E3,  # [mL];  https://pubmed.ncbi.nlm.nih.gov/22143261/
    Vp = 3.126E3,  # [mL]; https://pubmed.ncbi.nlm.nih.gov/22143261/
    Ve_blood = 2.558E3, # [mL]; https://pubmed.ncbi.nlm.nih.gov/22143261/
    # CLrenal = 2040., # renal clearance, [mL/h]
    # CLhepatic = 22.5E3, # hepatic clearance, [mL/h]
    CLrenal = 0., 
    CLhepatic = 0., 
    PER = 0.000756  # cytoplasmic membrane permeability coefficient [cm/h]
); 

# define initial condition 
N_Organs = 8
u0 = ComponentArray(
    Ca = 0., 
    Cp = 0., 
    Ce_blood = 0., 
    Cet_blood = 0., 
    Ci_org = zeros(N_Organs), 
    Cet_org = zeros(N_Organs), 
    totalclearance = 0. 
);

# dose 
MW = 543.5; # [g/mol]
Dose = 50E-3 * (1.7^2); # [g];
u0.Cp = Dose/MW*1E6/p_dox.Vp; # [umol/mL]

sol = solve(ODEProblem(dox_fullpbpk!, u0, (0., 25), p_dox));
Heart = 6;
Kd = 0.13; 
DNA_conc = 4.57; 
c_int = [sol.u[i].Ci_org[Heart] for i in 1:length(sol.t)]
c_intracellular_free = [0.5 * ( (sol.u[i].Cet_org[Heart] - DNA_conc - Kd ) + ( (sol.u[i].Cet_org[Heart] - DNA_conc - Kd).^2 .+ 4*Kd*sol.u[i].Cet_org[Heart] ).^(0.5) ) for i in 1:length(sol.t)]; 
c_nucleus = [sol.u[i].Cet_org[Heart] for i in 1:length(sol.t)] - c_intracellular_free;

p = plot(title = "heart", titlefontsize = 8, yaxis = :log, ylims = [1E-3, 1E2], legend = :outerright);
plot!(sol.t, sol[:Cp]*1E3, label = "plasma");
plot!(sol.t, c_int*1E3, label = "intersititum");
plot!(sol.t, c_intracellular_free*1E3, label = "intracellular free");
plot!(sol.t, c_nucleus*1E3, label = "nucleus");
plot!(xlabel = "Time (hour)", ylabel = "Conc (uM)"); 
p
