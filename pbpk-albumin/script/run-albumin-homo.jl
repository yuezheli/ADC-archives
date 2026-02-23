# date: 2/23/2026 
# author: Yuezhe Li 
# purpose of this code: to run the model in human 
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
    time_hr = [3.57, 10.1, 28.35, 59.7, 194.41, 265.06, 384.13, 576.52],
    conc_nM = [277.11, 260.39, 198.86, 155.0, 88.33, 73.19, 56.92, 41.52]
);

# load constants 
include(@projectroot("script/constants.jl"));

# load model 
include(@projectroot("model/twopore_albumin.jl"));
@named pbpk = albumin_pbpk_mtk();
pbpk_sys = structural_simplify(pbpk);

# Set the time span 
tspan = (0.0, 800.0);

# update parameters from mouse to human 
function updatevolume_human(mdl)
    p_homo = Dict([
        mdl.BW => 71,  # [unit = u"kg"]
        mdl.V_Plasma => 3.126, 
        mdl.V_LN     => 0.274, 
        # Vascular volumes (i.e. plasma volume), [L]; homo, 71kg
        mdl.V_V_Heart => 13.1/1000, 
        mdl.V_V_Lung => 55/1000, 
        mdl.V_V_Liver => 183/1000, 
        mdl.V_V_Muscle => 662/1000, 
        mdl.V_V_Skin => 127/1000, 
        mdl.V_V_Adipose => 148/1000, 
        mdl.V_V_Bone => 224/1000, 
        mdl.V_V_Brain => 31.9/1000, 
        mdl.V_V_Kidney => 18.2/1000, 
        mdl.V_V_SI => 6.15/1000, 
        mdl.V_V_LI => 8.74/1000, 
        mdl.V_V_Pancreas => 5.7/1000, 
        mdl.V_V_Thymus => 0.353/1000, 
        mdl.V_V_Spleen => 26.8/1000, 
        mdl.V_V_Other => 204/1000, 
        # Interstitial Volumes, [L]; homo, 71kg
        mdl.V_IS_Lung => 300/1000,  
        mdl.V_IS_Liver => 429/1000, 
        mdl.V_IS_Heart => 48.8/1000,
        mdl.V_IS_Muscle => 3910/1000, 
        mdl.V_IS_Skin => 1125/1000, 
        mdl.V_IS_Adipose => 2289/1000, 
        mdl.V_IS_Bone => 1891/1000, 
        mdl.V_IS_Brain => 261/1000, 
        mdl.V_IS_Kidney => 49.8/1000, 
        mdl.V_IS_SI => 67.1/1000, 
        mdl.V_IS_LI => 95.3/1000, 
        mdl.V_IS_Pancreas => 18/1000, 
        mdl.V_IS_Thymus => 1.09/1000, 
        mdl.V_IS_Spleen => 44.3/1000, 
        mdl.V_IS_Other => 831/1000, 
        # Endosomal Volumes, [L]; homo, 71kg
        mdl.V_E_Heart => 1.71/1000, 
        mdl.V_E_Lung => 5.0/1000, 
        mdl.V_E_Muscle => 150.0/1000, 
        mdl.V_E_Skin => 17.0/1000, 
        mdl.V_E_Adipose => 67.3/1000, 
        mdl.V_E_Bone => 50.8/1000, 
        mdl.V_E_Brain => 7.25/1000, 
        mdl.V_E_Kidney => 1.66/1000, 
        mdl.V_E_Liver => 10.7/1000, 
        mdl.V_E_SI => 1.93/1000, 
        mdl.V_E_LI => 2.74/1000, 
        mdl.V_E_Pancreas => 0.518/1000, 
        mdl.V_E_Thymus => 0.0321/1000, 
        mdl.V_E_Spleen => 1.11/1000, 
        mdl.V_E_Other => 24.3/1000, 
        # Blood Flows, [L/h]; homo, 71kg
        mdl.Q_Liver => 13210/1000, 
        mdl.Q_Heart => 7752/1000, 
        mdl.Q_Muscle => 33469/1000, 
        mdl.Q_Skin => 11626/1000, 
        mdl.Q_Adipose => 11233/1000, 
        mdl.Q_Bone => 2591/1000, 
        mdl.Q_Brain => 21453/1000, 
        mdl.Q_Kidney => 36402/1000,
        mdl.Q_SI => 12368/1000, 
        mdl.Q_LI => 12867/1000, 
        mdl.Q_Pancreas => 3056/1000, 
        mdl.Q_Thymus => 353/1000, 
        mdl.Q_Spleen => 6343/1000, 
        mdl.Q_Other => 9190/1000 , 
        # human glomerular filtration rate
        mdl.GFR => 7.2,  # [L/h]
        # others 
        mdl.C_0_endo_Alb => 6E-4, 
        mdl.Ksyn_alb => 8.1E-6, # [M/hr]
        mdl.Spino_alb => 0.185, 
        mdl.kon_FcRn => 2.7E7,  # [1/M/hr]
        mdl.koff_FcRn => 30.24, # [1/hr]
    ]);
    return p_homo
end
p_homo = updatevolume_human(pbpk_sys); 

# Create and Solve the ODEProblem
# The parameters (pars) are already baked into the system defaults.
u0_exg = Dict([pbpk_sys.C_Plasma_Exo => 480E-9]);
prob = ODEProblem(pbpk_sys, merge(p_homo, u0_exg), tspan);

# Using Rodas5 for stiff PBPK systems is generally recommended for stability
sol = solve(prob, Rodas5(), reltol=1e-6, abstol=1e-9)

# visualization 
plt_plasma = plot(xlabel="Time (hr)", ylabel="Conc (nM)", yaxis = :log10, yticks = [1E-1, 1, 10, 1E2, 1E3], ylims = [1E-1, 1E3]);
plot!(sol.t, sol[pbpk.C_Plasma_Exo] * 1E9, label = "sims");
plot!(plasma_data.time_hr, plasma_data.conc_nM, seriestype = :scatter, label = "Liu et al., 2024"); 
display(plt_plasma)

savefig(plt_plasma, @projectroot("deliv/figure/plasma-pk-albumin-homo.png"));


