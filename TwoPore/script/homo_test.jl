# author: Yuezhe Li 
# date: Dec 1, 2023 
# purpose: to test the model with homo physiology param and a dummy tumor component 

using Pkg; Pkg.activate("..");
# Pkg.instantiate(); # run this line once when setting up the repo for the first time 

# load ode
include("../model/twopore_tumor_homo.jl");

using DifferentialEquations 
using Plots
using DataFrames, CSV

p_IgG = ComponentArray(
    MW = 15E4, 
    k_deg = 32.2, # [hr-1]; mouse value; assumed to be a const
    infusion = 0., 
    Vc = 4/3 * pi * (5E-6)^3 * 1E3, # tumor cell volume [L], 
    P_protein = 334/24., # rate of permeability [um/h]; https://pubmed.ncbi.nlm.nih.gov/23151991/
    D_protein = 0.022/24 , # rate of diffusion [cm^2/h]; https://pubmed.ncbi.nlm.nih.gov/23151991/
    Rcap = 8. , # radium of tumor blood capillary [um]; https://pubmed.ncbi.nlm.nih.gov/23151991/
    Rkrogh = 75. , # adius of average distance between 2 blood vessels [um]; https://pubmed.ncbi.nlm.nih.gov/23151991/
    epsilon = 0.24, # partition coeffficient of ADC between plasma:tumor; https://pubmed.ncbi.nlm.nih.gov/23151991/
);

p_BSA = ComponentArray(
    MW = 66.4E3, 
    k_deg = 32.2, # [hr-1]; mouse value; assumed to be a const
    infusion = 0., 
    Vc = 4/3 * pi * (5E-6)^3 * 1E3, # tumor cell volume [L], 
    P_protein = 9.4E-7 * 1E4 * 3600, # rate of permeability [um/h]; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5435480/
    D_protein = 1.0E-6 * 3600, # rate of diffusion [cm^2/h]; https://www.iosrphr.org/papers/vol8-issue3/F0803014146.pdf
    Rcap = 8. , # radium of tumor blood capillary [um]; https://pubmed.ncbi.nlm.nih.gov/23151991/
    Rkrogh = 75. , # adius of average distance between 2 blood vessels [um]; https://pubmed.ncbi.nlm.nih.gov/23151991/
    epsilon = 0.24, # partition coeffficient of ADC between plasma:tumor; https://pubmed.ncbi.nlm.nih.gov/23151991/
);

# define init condition 
Dose_in_mgkg = 1.
dose_nM = Dose_in_mgkg*1E6*71/15E4;
 
u0_IgG = ComponentArray(C_EXG_Plasma = dose_nM/3.126, 
                    C_EXG_LN = 0., C_EXG = zeros(15, 3), endo_deg = 0., renal_clearance = 0., 
                    C_EXG_Tumor = 0., Nc1 = 10000., Nc2 = 0., Nc3 = 0., Nc4 = 0.);

u0_BSA = ComponentArray(C_EXG_Plasma = dose_nM/3.126, 
                    C_EXG_LN = 0., C_EXG = zeros(15, 3), endo_deg = 0., renal_clearance = 0., 
                    C_EXG_Tumor = 0., Nc1 = 10000., Nc2 = 0., Nc3 = 0., Nc4 = 0.);

# solve ODE
sol_IgG = solve(ODEProblem(lishah_tumor_homo!, u0_IgG, (0., 150.), p_IgG), alg=AutoTsit5(Rosenbrock23(autodiff = false)), saveat = 0.5);
cplasma_IgG = [sol_IgG.u[i].C_EXG_Plasma for i in 1:length(sol_IgG.t)]; # [nM]
ctumor_IgG = [sol_IgG.u[i].C_EXG_Tumor for i in 1:length(sol_IgG.t)]; # [nM]

sol_BSA = solve(ODEProblem(lishah_tumor_homo!, u0_BSA, (0., 150.), p_BSA), alg=AutoTsit5(Rosenbrock23(autodiff = false)), saveat = 0.5);
cplasma_BSA = [sol_BSA.u[i].C_EXG_Plasma for i in 1:length(sol_BSA.t)]; # [nM]
ctumor_BSA = [sol_BSA.u[i].C_EXG_Tumor for i in 1:length(sol_BSA.t)]; # [nM]

# visualization
p_conc = plot(xlabel = "Time (h)", ylabel = "concentration (nM)");
plot!(sol_IgG.t, cplasma_IgG, label = "plasma, IgG (150kDa)", color = :red, alpha = 0.7, linestyle = :solid);
plot!(sol_IgG.t, ctumor_IgG, label = "tumor (dummy), IgG (150kDa)", color = :red, alpha = 0.7, linestyle = :dash);
plot!(sol_BSA.t, cplasma_BSA, label = "plasma, BSA (66kDa)", color = :blue, alpha = 0.7, linestyle = :solid);
plot!(sol_BSA.t, ctumor_BSA, label = "tumor (dummy), BSA (66kDa)", color = :blue, alpha = 0.7, linestyle = :dash);
display(p_conc);
