# date: 4/3/24
# author: Yuezhe Li 
# purpose of this script: to troubleshoot the model from He et al., 2019 based on mass balance 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6533104/

using Pkg; Pkg.activate("");

using DifferentialEquations, ComponentArrays
using Plots, DataFrames, CSV, DataFramesMeta

include("model/dox_fullpbpk.jl");

# define params (distribution to the cells turned off)
p_dox = ComponentArray(
    infusion = 0., 
    Kpt_blood = 0.68, 
    Va = 3.126E3,  # [mL];  https://pubmed.ncbi.nlm.nih.gov/22143261/
    Vp = 3.126E3,  # [mL]; https://pubmed.ncbi.nlm.nih.gov/22143261/
    Ve_blood = 2.558E3, # [mL]; https://pubmed.ncbi.nlm.nih.gov/22143261/
    CLrenal = 8190.748, # renal clearance, scaled from mouse, [mL/h]
    CLhepatic = 25.5E3, # hepatic clearance, [mL/h]
    PER = 0.
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
u0.Cp = 1/p_dox.Vp; # [umol/mL]

sol = solve(ODEProblem(dox_fullpbpk!, u0, (0., 25), p_dox));

Lung = 1; Liver = 2; Gut = 3; Spleen = 4; Kidney = 5; Heart = 6; Other = 7; Tumor = 8
Vi_org = zeros(N_Organs); 
Vi_org[Lung] = 234.; Vi_org[Liver] = 275.; Vi_org[Gut] = 287.; Vi_org[Spleen] = 38.;
Vi_org[Kidney] = 42.; Vi_org[Heart] = 44.; Vi_org[Other] = 1969.; Vi_org[Tumor] = 2.;

totalsystem = [ (sol.u[i].Ca*p_dox.Va + sol.u[i].Cp*p_dox.Vp + sol.u[i].Ce_blood*p_dox.Ve_blood + sol.u[i].Cet_blood*p_dox.Ve_blood + 
                 sum(sol.u[i].Ci_org .* Vi_org) + sol.u[i].totalclearance) for i in 1:length(sol.t) ]

plot(sol.t, totalsystem, ylims = [0, 2])
