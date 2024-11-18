# author: Yuezhe LI
# date: Sep 25, 2023

using DifferentialEquations, RecursiveArrayTools, DataFrames, CSV, DataFramesMeta
using Plots, Statistics
using ComponentArrays
using PDFmerger

include("jones_tumor_as_another_organ.jl") 
include("params.jl")
include("helper.jl")

p16_tdm1 = deepcopy(p_N_16);
p16_tdm1.Rcopies = 1.0E6
# p16_tdm1.Rcopies = 0.  # for troubleshooting 
p16_tdm1.PS_Score = 6.
p16_tdm1.init_sR = 0.004
# p16_tdm1.init_sR = 0.  # for troubleshooting 
p16_tdm1.thalf_sR_adc = 120.

tspan = (0.0, DayToHour*84);  # hr

Dose = 3.6 ; #[mg/kg]

u16, dose_uM = solid_tumor_init(N_Organs, TotalCells, 1.0, Dose*1E3, p16_tdm1.Rcopies, p16_tdm1.init_sR, p16_tdm1.k_endo, p16_tdm1.k_rec);

sol16 = solve(ODEProblem(jones_tumor_as_another_organ!, u16, tspan, p16_tdm1), saveat = 1.0, alg = QNDF(autodiff=false), reltol = 1e-15);
# cplasma = [sol16.u[i].C_EXG_Plasma for i in 1:length(sol16.t)] * MW_EDG ; # unit: ug/L

## total ADC in tumor [umol]
tumorm1 = [ (sol16.u[i].tumor.AR_s + sol16.u[i].tumor.AR_e) for i in 1:length(sol16.t)] ;

## degraded mAb mass in tumor cells [umol]
dt1 = [sol16.u[i].DEGprotein.deg_TumorCellular for i in 1:length(sol16.t)];

## total ADC in tissue other than tumor [umol]
tm1 = TissueMass(sol16, V_V, V_VM, V_E7, V_E6a, V_E7b, V_ISM, V_IntS);

## track all ADC in the model 
tm1[:, 16] .+= tumorm1;  # add tumor ADC
dadc = EpithelialDegradedProtein(sol16);   # track degraded ADC in organ endothelials
dadc[:, 16] .+= dt1;   # track ADC degraded inside tumor cells
total_adc = hcat(tm1 .+ dadc,                                                                                           # add degraded and existed ADC
                [(sol16.u[i].C_EXG_Plasma * V_Plasma .+ sol16.u[i].DEGprotein.deg_plasma) for i in 1:length(sol16.t)]); # add plasma ADC

## create a table
dist_adc = DataFrame(
    Organ = [# "Plasma",  "Other", 
            "Lung", "Liver", "Heart", "Muscle", "Skin", "Adipose", "Brain", "Kidney", "Intestin", "Spleen","Tumor"], 
    frac = [
        #total_adc[end, 17]/(dose_uM), 
        #(total_adc[end, 7] + total_adc[end, 12] + total_adc[end, 13] + total_adc[end, 15])/(dose_uM),
        total_adc[end, 1]/(dose_uM), 
        total_adc[end, 2]/(dose_uM), 
        total_adc[end, 3]/(dose_uM), 
        total_adc[end, 4]/(dose_uM),
        total_adc[end, 5]/(dose_uM),  
        total_adc[end, 6]/(dose_uM), 
        total_adc[end, 8]/(dose_uM),
        total_adc[end, 9]/(dose_uM),
        (total_adc[end, 10] + total_adc[end, 11])/(dose_uM),
        total_adc[end, 14]/(dose_uM),
        total_adc[end, 16]/(dose_uM), 
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

savefig(p_total_mass_tissue, "solid-tumor-as-an-organ.png");

