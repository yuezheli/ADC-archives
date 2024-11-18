# author: Yuezhe Li 
# date: Sep 25, 2023

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Plots, Statistics, DataFrames, CSV, DataFramesMeta

# simulation 
include("params.jl"); 
include("jones_tumor_nodal.jl") 
include("helper.jl"); 

pnodel_tdm1 = deepcopy(p_N_16);
pnodel_tdm1.Rcopies = 1.0E6
pnodel_tdm1.PS_Score = 6.
pnodel_tdm1.init_sR = 0.004
pnodel_tdm1.thalf_sR_adc = 120.

tspan = (0., DayToHour*84);      # [hr]
Dose = 3.6  # [mg/kg]

u_nodal, dose_uM = solid_tumor_init(15, TotalCells, 1.0, Dose*1E3, pnodel_tdm1.Rcopies, pnodel_tdm1.init_sR, pnodel_tdm1.k_endo, pnodel_tdm1.k_rec, true);

sol_nodal = solve(ODEProblem(jones_tumor_nodal!, u_nodal, tspan, pnodel_tdm1), saveat = 3., alg = QNDF(autodiff=false), reltol = 1e-12);

## total ADC in tumor [umol]
tumorm1 = [ (sol_nodal.u[i].tumor.AR_s + sol_nodal.u[i].tumor.AR_e) for i in 1:length(sol_nodal.t)] ;

## degraded mAb mass in tumor cells [umol]
dt1 = [sol_nodal.u[i].DEGprotein.deg_TumorCellular for i in 1:length(sol_nodal.t)];

## total ADC in tissue other than tumor [umol]
tm1 = TissueMass(sol_nodal, V_V, V_VM, V_E7, V_E6a, V_E7b, V_ISM, V_IntS);

## create a table
dist_adc = DataFrame(
    Organ = [ "Lung", "Liver", "Heart", "Muscle", "Skin", "Adipose", "Brain", "Kidney", "Intestin", "Spleen", "Tumor"], 
    frac = [
        (tm1[end, 1] .+ sol_nodal.u[end].DEGprotein.deg_EXG[1])/(dose_uM), 
        (tm1[end, 2] .+ sol_nodal.u[end].DEGprotein.deg_EXG[2])/(dose_uM), 
        (tm1[end, 3] .+ sol_nodal.u[end].DEGprotein.deg_EXG[3])/(dose_uM), 
        (tm1[end, 4] .+ sol_nodal.u[end].DEGprotein.deg_EXG[4])/(dose_uM),
        (tm1[end, 5] .+ sol_nodal.u[end].DEGprotein.deg_EXG[5])/(dose_uM),  
        (tm1[end, 6] .+ sol_nodal.u[end].DEGprotein.deg_EXG[6])/(dose_uM), 
        (tm1[end, 8] .+ sol_nodal.u[end].DEGprotein.deg_EXG[8])/(dose_uM),
        (tm1[end, 9] .+ sol_nodal.u[end].DEGprotein.deg_EXG[9])/(dose_uM),
        (tm1[end, 10] .+ sol_nodal.u[end].DEGprotein.deg_EXG[10])/(dose_uM) + (tm1[end, 11] .+ sol_nodal.u[end].DEGprotein.deg_EXG[11])/(dose_uM),
        (tm1[end, 14] .+ sol_nodal.u[end].DEGprotein.deg_EXG[14])/(dose_uM),
        (tumorm1[end] .+ dt1[end])/(dose_uM), 
    ]
);

sort!(dist_adc, [:frac], rev = true)

p_total_mass_tissue = bar(dist_adc.Organ, dist_adc.frac * 100, label = false,
    xticks =:all, xrotation = 45,
    size = [350, 250], yaxis = :log, 
    xlabel ="", ylabel = "% of the total drug", 
    ylims = (1E-2, 1E1), 
    xtickfontsize=10, 
    dpi = 1000)

savefig(p_total_mass_tissue, "nodal-tumor.png");
