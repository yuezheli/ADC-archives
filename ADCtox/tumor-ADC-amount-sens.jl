# Author: Yuezhe Li
# date: 7/31/23
# purpose: sensitivity analysis on the percentage of ADC ended up in tumor 

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Plots, Statistics, DataFrames, CSV, DataFramesMeta

include("jones_tumor_homo.jl") 
include("params.jl"); 
include("helper.jl"); 

# update params (default params, for T-DM1)
p_base.Rcopies = 1.0E6
p_base.PS_Score = 6.
p_base.init_sR = 0.004
p_base.thalf_sR_adc = 120.

TotalCells = 1.13e9              # [cell count]

tspan = (-0.1, DayToHour*84);    # [hr]

# repeated dosing 
AddDose = [0.] * DayToHour        # [hr]
Dose_in_ugkg = ones(size(AddDose)) * 3.6E3  # [ug/kg]
infusion_time = 0.5                         # [hr]
dose_uM = sum(Dose_in_ugkg)*BW/MW_EDG

global cbs0 = [];
if length(AddDose) >= 1
    for i in 1:length(AddDose)
        function affect_infusion_on!(integrator)
            integrator.p.infusion = Dose_in_ugkg[i]*BW/(V_Plasma)/MW_EDG/infusion_time;
        end
        function affect_infusion_off!(integrator)
            integrator.p.infusion = 0.;
        end
        cb01 = PresetTimeCallback(AddDose[i],affect_infusion_on!);
        cb02 = PresetTimeCallback( (AddDose[i]+infusion_time),affect_infusion_off!);
        global cbs0 = push!(cbs0, cb01);
        global cbs0 = push!(cbs0, cb02);
    end
end
cbset0 = CallbackSet(cbs0...);

function dist_adc(p_base, cbset0 = cbset0)
    u0, _ = jones_init(TotalCells, 1.0, 0.0, p_base.Rcopies, p_base.init_sR, p_base.k_endo, p_base.k_rec);
    sol1 = solve(ODEProblem(jonesODEs_homo_tumor!, u0, tspan, p_base), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-12, callback = cbset0);

    tumorm1 = [ (sol1.u[i].tumor.A_m + sol1.u[i].tumor.AR_s + sol1.u[i].tumor.AR_e) for i in 1:length(sol1.t)] ;

    ## degraded mAb mass in tumor cells [umol]
    dt1 = [sol1.u[i].DEGprotein.deg_TumorCellular for i in 1:length(sol1.t)];

    ## total ADC in tissue other than tumor [umol]
    tm1 = TissueMass(sol1, V_V, V_VM, V_E7, V_E6a, V_E7b, V_ISM, V_IntS);

    dist_adc1 = DataFrame(
        Organ = ["Plasma", "Liver", "Muscle", "Lung", "Adipose", "Skin", "Brain", "Kidney",  "Heart", "Tumor"], 
        frac = [
            (sol1.u[end].C_EXG_Plasma * V_Plasma .+ sol1.u[end].DEGprotein.deg_plasma)/(dose_uM*size(AddDose)[1]), 
            (tm1[end, 2] .+ sol1.u[end].DEGprotein.deg_EXG[2])/(dose_uM*size(AddDose)[1]), 
            (tm1[end, 4] .+ sol1.u[end].DEGprotein.deg_EXG[4])/(dose_uM*size(AddDose)[1]), 
            (tm1[end, 1] .+ sol1.u[end].DEGprotein.deg_EXG[1])/(dose_uM*size(AddDose)[1]), 
            (tm1[end, 6] .+ sol1.u[end].DEGprotein.deg_EXG[6])/(dose_uM*size(AddDose)[1]), 
            (tm1[end, 5] .+ sol1.u[end].DEGprotein.deg_EXG[5])/(dose_uM*size(AddDose)[1]), 
            (tm1[end, 8] .+ sol1.u[end].DEGprotein.deg_EXG[8])/(dose_uM*size(AddDose)[1]),
            (tm1[end, 9] .+ sol1.u[end].DEGprotein.deg_EXG[9])/(dose_uM*size(AddDose)[1]),
            (tm1[end, 3] .+ sol1.u[end].DEGprotein.deg_EXG[3])/(dose_uM*size(AddDose)[1]), 
            (tumorm1[end] .+ dt1[end])/(dose_uM*size(AddDose)[1]), 
        ]
    );
    return dist_adc1
end

function efficacy_adc(p_base, cbset0 = cbset0)
    u0, _ = jones_init(TotalCells, 1.0, 0.0, p_base.Rcopies, p_base.init_sR, p_base.k_endo, p_base.k_rec);
    sol1 = solve(ODEProblem(jonesODEs_homo_tumor!, u0, tspan, p_base), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-12, callback = cbset0);

    Rpos_cell = [ (sol1.u[i].tumor.Nc_1 + sol1.u[i].tumor.Nc_2 + sol1.u[i].tumor.Nc_3 + sol1.u[i].tumor.Nc_4) for i in 1:length(sol1.t)] ;
    Rneg_cell = [ (sol1.u[i].tumor.Nc_1_neg + sol1.u[i].tumor.Nc_2_neg + sol1.u[i].tumor.Nc_3_neg + sol1.u[i].tumor.Nc_4_neg) for i in 1:length(sol1.t) ];
    tv = (Rpos_cell .+ Rneg_cell) * Vc;
    tv0 = tv[1];

    return tv/tv0;
end

# soluble binding partner scan 
sHER2_conc = [4.0E-3, 0., 4.0E-4, 0.001, 2.0E-3, 5.0E-3, 8.0E-3, 0.01];
tumor_frac_sHER2_conc = []
for sHER2 in sHER2_conc
    p_tmp = deepcopy(p_base); p_tmp.init_sR = sHER2
    dist_adc_tmp = dist_adc(p_tmp);
    append!(tumor_frac_sHER2_conc, dist_adc_tmp.frac[end]);
end

p_tumor_frac_sHER2_conc = plot(xticks = (sHER2_conc, string.(sHER2_conc)), xrotation = 45, size=(350,300));
scatter!(sHER2_conc, tumor_frac_sHER2_conc * 100, label = false);
xlabel!("soluble receptor concentration in plasma (uM)"); ylabel!("ADC ends up in tumor (%)")

savefig(p_tumor_frac_sHER2_conc, "img/param-scan/sHER2-tumor-ADC-percentage.png");

# recepor expression level scan
HER2copies = [1.0E3, 2.0E3, 5.0E3, 1.0E4, 1.0E5, 1.0E6, 1.0E7, 1.0E8];
tumor_frac_HER2copies = []
ntv_HER2copies = []
for Her2copy in HER2copies
    p_tmp = deepcopy(p_base); p_tmp.Rcopies = Her2copy
    dist_adc_tmp = dist_adc(p_tmp);
    ntv_tmp = efficacy_adc(p_tmp);
    append!(tumor_frac_HER2copies, dist_adc_tmp.frac[end]);
    append!(ntv_HER2copies, ntv_tmp[end]);
end

p_tumor_frac_HER2copies = plot(xticks = (HER2copies, ["1E3", "2E3", "5E3", "1E4", "1E5", "1E6", "1E7", "1E8"]), xaxis = :log10, xrotation = 90, size=(350,250), dpi = 1000);
scatter!(HER2copies, tumor_frac_HER2copies * 100, label = false, markersize = 6);
xlabel!("surface receptor copy number (#)"); ylabel!("% ADC ends up in tumor");
ylims!(0.008, 0.012)

savefig(p_tumor_frac_HER2copies, "img/param-scan/HER2-tumor-ADC-percentage.png");

p_ntv_HER2copies = plot(xticks = (HER2copies, ["1E3", "2E3", "5E3", "1E4", "1E5", "1E6", "1E7", "1E8"]), xaxis = :log10, xrotation = 90, size=(350,300), dpi = 1000);
scatter!(HER2copies, (1 .- ntv_HER2copies) * 100, label = false);
xlabel!("surface receptor copy number (#)"); ylabel!("% tumor volume reduction");
ylims!(10, 30)

savefig(p_ntv_HER2copies, "img/param-scan/HER2-tumor-volume-reduction.png");

# HER2:ADC complex internalization rate scan 
k_c_endocytosis = [0.01, 0.05, 0.1, 0.2, 0.5, 1.5];
tumor_frac_k_endo = []
for k_endo in k_c_endocytosis
    p_tmp = deepcopy(p_base); p_tmp.k_endo = k_endo
    dist_adc_tmp = dist_adc(p_tmp);
    append!(tumor_frac_k_endo, dist_adc_tmp.frac[end]);
end

p_tumor_frac_k_endo = plot(xticks = (k_c_endocytosis, string.(k_c_endocytosis)), xaxis = :log10, xrotation = 90, size=(350,300), ylims = (0.008, 0.012), dpi = 1000);
scatter!(k_c_endocytosis, tumor_frac_k_endo * 100, label = false, markersize = 6);
xlabel!("receptor:ADC endocytosis rate (h-1)"); ylabel!("% ADC ends up in tumor")

savefig(p_tumor_frac_k_endo, "img/param-scan/endocytosis-tumor-ADC-percentage.png");

# HER2:ADC binding affinity
Kd_ADC = [0.3E-5, 0.3E-4, 0.3E-3, 0.3E-2, 0.3];
tumor_frac_Kd_ADC = [];
ntv_Kd_ADC = [];
for Kd in Kd_ADC
    p_tmp = deepcopy(p_base); p_tmp.Kd = Kd
    dist_adc_tmp = dist_adc(p_tmp);
    append!(tumor_frac_Kd_ADC, dist_adc_tmp.frac[end]);
    ntv_tmp = efficacy_adc(p_tmp);
    append!(ntv_Kd_ADC, ntv_tmp[end]);
end

p_tumor_frac_Kd_ADC = plot(xticks = (Kd_ADC, string.(Kd_ADC)), xaxis = :log10, xrotation = 20);
scatter!(Kd_ADC, tumor_frac_Kd_ADC * 100, label = false);
xlabel!("receptor:ADC binding affinity (uM)"); ylabel!("ADC ends up in tumor (%)")

savefig(p_tumor_frac_Kd_ADC, "img/param-scan/Kd-tumor-ADC-percentage.png");

p_ntv_Kd_ADC = plot(xticks = (Kd_ADC, string.(Kd_ADC)), xaxis = :log10, xrotation = 20);
scatter!(Kd_ADC, (1 .- ntv_Kd_ADC) * 100, label = false);
xlabel!("receptor:ADC binding affinity (uM)"); ylabel!("Normalized tumor volume reduction (%)")

savefig(p_ntv_Kd_ADC, "img/param-scan/Kd-target-ADC-tumor-volume-reduction.png");


# tumor perfusion scan
Rkrogh = [1., 5., 10., 25., 50., 75., 100.];
tumor_frac_Rkrogh = []
for Rk in Rkrogh
    p_tmp = deepcopy(p_base); p_tmp.Rkrogh = Rk
    dist_adc_tmp = dist_adc(p_tmp);
    append!(tumor_frac_Rkrogh, dist_adc_tmp.frac[end]);
end

p_tumor_frac_Rkrogh = plot(xticks = (Rkrogh, string.(Rkrogh)), xaxis = :log10, yaxis = :log10, xrotation = 90, size=(350,300), dpi = 1000);
scatter!(Rkrogh, tumor_frac_Rkrogh * 100, label = false, markersize = 6);
xlabel!("dist between blood vessels (um)"); ylabel!("% ADC ends up in tumor")

savefig(p_tumor_frac_Rkrogh, "img/param-scan/Rkrogh-tumor-ADC-percentage.png");


permeability_adc = [1., 2., 4., 8., 12., 16., 20.]; 

tumor_frac_p_adc = []
for p_adc in permeability_adc
    p_tmp = deepcopy(p_base); p_tmp.P_ADC = p_adc
    dist_adc_tmp = dist_adc(p_tmp);
    append!(tumor_frac_p_adc, dist_adc_tmp.frac[end]);
end

p_tumor_frac_p_adc = plot(xticks = (permeability_adc, string.(permeability_adc)), xrotation = 90, size=(350,300), dpi = 1000);
scatter!(permeability_adc, tumor_frac_p_adc * 100, label = false, markersize = 6);
xlabel!("ADC diffusivity (um/h)"); ylabel!("% ADC ends up in tumor")

savefig(p_tumor_frac_p_adc, "img/param-scan/P-ADC-tumor-ADC-percentage.png");

