# author: Yuezhe Li 
# date: Sep 25, 2023

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Plots, Statistics, DataFrames, CSV, DataFramesMeta

# simulation 
include("params.jl"); 
include("jones_tumor_shah2012.jl") 
include("helper.jl"); 

p15_tdm1 = deepcopy(p_N_15);
p15_tdm1.Rcopies = 1.0E6
p15_tdm1.PS_Score = 6.
p15_tdm1.init_sR = 0.004
p15_tdm1.thalf_sR_adc = 120.

#---------------------------------- T-DM1 in tissue distribution ----------------------------------#
tspan = (0., DayToHour*84);      # [hr]
Dose = 3.6  # [mg/kg]

u15, dose_uM = solid_tumor_init(N_Organs, TotalCells, 1.0, Dose*1E3, p15_tdm1.Rcopies, p15_tdm1.init_sR, p15_tdm1.k_endo, p15_tdm1.k_rec);

sol15 = solve(ODEProblem(jones_tumor_shah2012!, u15, tspan, p15_tdm1), saveat = 3., alg = QNDF(autodiff=false), reltol = 1e-12);

## total ADC in tumor [umol]
tumorm1 = [ (sol15.u[i].tumor.A_m + sol15.u[i].tumor.AR_s + sol15.u[i].tumor.AR_e) for i in 1:length(sol15.t)] ;

## degraded mAb mass in tumor cells [umol]
dt1 = [sol15.u[i].DEGprotein.deg_TumorCellular for i in 1:length(sol15.t)];

## total ADC in tissue other than tumor [umol]
tm1 = TissueMass(sol15, V_V, V_VM, V_E7, V_E6a, V_E7b, V_ISM, V_IntS);

## create a table
dist_adc = DataFrame(
    Organ = [# "Plasma", 
             "Lung", "Liver", "Heart", "Muscle", "Skin", "Adipose", "Brain", "Kidney", "Intestin", "Spleen", "Tumor"], 
    frac = [
        #(sol15.u[end].C_EXG_Plasma * V_Plasma .+ sol15.u[end].DEGprotein.deg_plasma)/(dose_uM), 
        (tm1[end, 1] .+ sol15.u[end].DEGprotein.deg_EXG[1])/(dose_uM), 
        (tm1[end, 2] .+ sol15.u[end].DEGprotein.deg_EXG[2])/(dose_uM), 
        (tm1[end, 3] .+ sol15.u[end].DEGprotein.deg_EXG[3])/(dose_uM), 
        (tm1[end, 4] .+ sol15.u[end].DEGprotein.deg_EXG[4])/(dose_uM),
        (tm1[end, 5] .+ sol15.u[end].DEGprotein.deg_EXG[5])/(dose_uM),  
        (tm1[end, 6] .+ sol15.u[end].DEGprotein.deg_EXG[6])/(dose_uM), 
        (tm1[end, 8] .+ sol15.u[end].DEGprotein.deg_EXG[8])/(dose_uM),
        (tm1[end, 9] .+ sol15.u[end].DEGprotein.deg_EXG[9])/(dose_uM),
        (tm1[end, 10] .+ sol15.u[end].DEGprotein.deg_EXG[10])/(dose_uM) + (tm1[end, 11] .+ sol15.u[end].DEGprotein.deg_EXG[11])/(dose_uM),
        (tm1[end, 14] .+ sol15.u[end].DEGprotein.deg_EXG[14])/(dose_uM),
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

savefig(p_total_mass_tissue, "solid-tumor-shah2012.png");
