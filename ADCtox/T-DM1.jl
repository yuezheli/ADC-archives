# author: Yuezhe Li 
# date: 6/27/23
# purpose: to adjust PK for T-DM1 

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Plots, Statistics, DataFrames, CSV, DataFramesMeta

# read observed 
pk_yamamoto =  CSV.read("data/yamamoto-2015.csv",DataFrame);
pk_girish =  CSV.read("data/girish-2012.csv",DataFrame);

pk_point3 = @rsubset(pk_girish, :Dose .== "0.3mg/kg" ); 
pk_point6 = @rsubset(pk_girish, :Dose .== "0.6mg/kg" ); 
pk_1point2 = @rsubset(pk_girish, :Dose .== "1.2mg/kg" ); 
pk_2point4_1 = @rsubset(pk_girish, :Dose .== "2.4mg/kg" ); 
pk_3point6_1 = @rsubset(pk_girish, :Dose .== "3.6mg/kg" ); 
pk_4point8 = @rsubset(pk_girish, :Dose .== "4.8mg/kg" ); 

pk_1point8 = @rsubset(pk_yamamoto, :Dose .== "1.8mg/kg" ); 
pk_2point4_2 = @rsubset(pk_yamamoto, :Dose .== "2.4mg/kg" ); 
pk_3point6_2 = @rsubset(pk_yamamoto, :Dose .== "3.6mg/kg" );

pk_2point4 = vcat(pk_2point4_1, pk_2point4_2);
pk_3point6 = vcat(pk_3point6_1, pk_3point6_2);

# simulation 
include("jones_tumor_homo.jl") 
include("params.jl"); 
include("helper.jl"); 

TotalCells = 1.13e9       # [cell count]; convert to 2mL of initial tumor size

p_tdm1 = deepcopy(p_base);
p_tdm1.Rcopies = 1.0E6
p_tdm1.PS_Score = 6.
p_tdm1.init_sR = 0.004
p_tdm1.thalf_sR_adc = 120.

tspan = (-0.01, DayToHour*42);      # [hr]
AddDose = [0., 21.] * DayToHour  # [hr]

first_dose = [0.3, 0.6, 1.2, 2.4, 3.6, 4.8, 1.8]; # [mg/kg]

#---------------------------------- PK ----------------------------------#

function OneDoseSims(init_dose, p_base = p_base, infusion_time = 0.5)
    u0_tmp, _ = jones_init(TotalCells, 1.0, 0., p_base.Rcopies, p_base.init_sR, p_base.k_endo, p_base.k_rec);

    global cbs0 = [];
    if length(AddDose) >= 1
        for i in 1:length(AddDose)
            function affect_infusion_on!(integrator)
                integrator.p.infusion = init_dose *1E3*BW/(V_Plasma)/MW_EDG/infusion_time;
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

    sol_tmp = solve(ODEProblem(jonesODEs_homo_tumor!, u0_tmp, tspan, p_base), saveat = 3., alg = QNDF(autodiff=false), reltol = 1e-12, callback = cbset0);
    cplasma_tmp = [sol_tmp.u[i].C_EXG_Plasma for i in 1:length(sol_tmp.t)]; # [uM]

    return sol_tmp.t, cplasma_tmp
end

timestamp = []; 
cplasma = [];

for init_dose in first_dose
    tmpt, tmpp = OneDoseSims(init_dose, p_tdm1, 0.5);
    append!(timestamp, [tmpt/DayToHour]);
    append!(cplasma, [tmpp * MW_EDG / 1E3]); # [ug/mL]
end


p_pk_tdm1 = plot(legend = :topright, ylims = [1E-1, 1E3], yaxis = :log, size=(350,350), legendcolumns=3, dpi = 1000);
plot!(timestamp[1], cplasma[1], lw = 2, alpha = 0.6, color = "red", label = "0.3mg/kg");
scatter!(pk_point3.time_day, pk_point3.T_DM1_ugperml, ma = 0.6, color = "red", label = false, linewidth=0);
plot!(timestamp[2], cplasma[2], lw = 2, alpha = 0.6, color = "seagreen", label = "0.6mg/kg");
scatter!(pk_point6.time_day, pk_point6.T_DM1_ugperml, ma = 0.6, color = "seagreen", label = false, linewidth=0);
plot!(timestamp[3], cplasma[3], lw = 2, alpha = 0.6, color = "mediumpurple", label = "1.2mg/kg");
scatter!(pk_1point2.time_day, pk_1point2.T_DM1_ugperml, ma = 0.6, color = "mediumpurple", label = false, linewidth=0);
plot!(timestamp[4], cplasma[4], lw = 2, alpha = 0.6, color = "violetred", label = "2.4mg/kg");
scatter!(pk_2point4.time_day, pk_2point4.T_DM1_ugperml, ma = 0.6, color = "violetred", label = false, linewidth=0);
plot!(timestamp[5], cplasma[5], lw = 2, alpha = 0.6, color = "blue2", label = "3.6mg/kg");
scatter!(pk_3point6.time_day, pk_3point6.T_DM1_ugperml, ma = 0.6, color = "blue2", label = false, linewidth=0);
plot!(timestamp[6], cplasma[6], lw = 2, alpha = 0.6, color = "teal", label = "4.8mg/kg");
scatter!(pk_4point8.time_day, pk_4point8.T_DM1_ugperml, ma = 0.6, color = "teal", label = false, linewidth=0);
xlabel!("time (day)"); ylabel!("plasma T-DM1 conc (ug/mL)"); xlims!(0, 30)


savefig(p_pk_tdm1, "img/t-dm1-pk.png");

#---------------------------------- T-DM1 in tissue interstitium ----------------------------------#
tspan = (0., DayToHour*84);      # [hr]
Dose = 3.6  # [mg/kg]

u0, _ = jones_init(TotalCells, 1.0, 0., p_tdm1.Rcopies, p_tdm1.init_sR, p_tdm1.k_endo, p_tdm1.k_rec);
u0.C_EXG_Plasma = Dose *1E3*BW/(V_Plasma)/MW_EDG; # [uM]

sol0 = solve(ODEProblem(jonesODEs_homo_tumor!, u0, tspan, p_tdm1), saveat = 3., alg = QNDF(autodiff=false), reltol = 1e-12);
adc_ints = InterstitialADC(sol0);

p_ints_adc = plot(xlabel = "time (day)", ylabel = "tissue interstitial \n ADC concentration (uM)", size=(350, 265), dpi = 1000); 
plot!(sol0.t/DayToHour, adc_ints.lung_ints, label = "Lung", lw = 2.5, alpha = 0.7); 
#plot!(sol0.t/DayToHour, adc_ints.liver_ints, label = "Liver"); 
#plot!(sol0.t/DayToHour, adc_ints.heart_ints, label = "Heart"); 
#plot!(sol0.t/DayToHour, adc_ints.muscle_ints, label = "Muscle"); 
plot!(sol0.t/DayToHour, adc_ints.skin_ints, label = "Skin", lw = 2.5, alpha = 0.7); 
#plot!(sol0.t/DayToHour, adc_ints.adipose_ints, label = "Adipose"); 
#plot!(sol0.t/DayToHour, adc_ints.bone_ints, label = "Bone"); 
#plot!(sol0.t/DayToHour, adc_ints.brain_ints, label = "Brain"); 
#plot!(sol0.t/DayToHour, adc_ints.kidney_ints, label = "Kidney"); 
plot!(sol0.t/DayToHour, adc_ints.si_ints, label = "Small intestine", lw = 2.5, alpha = 0.7); 
#plot!(sol0.t/DayToHour, adc_ints.li_ints, label = "Large intestin"); 
#plot!(sol0.t/DayToHour, adc_ints.pancreas_ints, label = "Pancreas"); 
#plot!(sol0.t/DayToHour, adc_ints.thymus_ints, label = "Thymus"); 
#plot!(sol0.t/DayToHour, adc_ints.spleen_ints, label = "Spleen"); 
hline!([8.6E-3], linestyle=:dashdot, label="IC50, T-DM1", color = :black, alpha = 0.5, lw = 2.5);  # https://ascopubs.org/doi/abs/10.1200/JCO.2018.36.4_suppl.256, for HER2-low cells
plot!(legend = :topright)

savefig(p_ints_adc, "img/t-dm1-tissue-ints.png");

#---------------------------------- T-DM1 in tissue distribution ----------------------------------#
dose_uM = Dose *1E3*BW/MW_EDG; # [uM]

## total ADC in tumor [umol]
tumorm1 = [ (sol0.u[i].tumor.A_m + sol0.u[i].tumor.AR_s + sol0.u[i].tumor.AR_e) for i in 1:length(sol0.t)] ;

## degraded mAb mass in tumor cells [umol]
dt1 = [sol0.u[i].DEGprotein.deg_TumorCellular for i in 1:length(sol0.t)];

## total ADC in tissue other than tumor [umol]
tm1 = TissueMass(sol0, V_V, V_VM, V_E7, V_E6a, V_E7b, V_ISM, V_IntS);

## create a table
dist_adc = DataFrame(
    Organ = [#"Plasma", 
            "Lung", "Liver", "Heart", "Muscle", "Skin", "Adipose", "Brain", "Kidney", "Intestin", "Spleen", "Tumor"], 
    frac = [
        # (sol0.u[end].C_EXG_Plasma * V_Plasma .+ sol0.u[end].DEGprotein.deg_plasma)/(dose_uM), 
        (tm1[end, 1] .+ sol0.u[end].DEGprotein.deg_EXG[1])/(dose_uM), 
        (tm1[end, 2] .+ sol0.u[end].DEGprotein.deg_EXG[2])/(dose_uM), 
        (tm1[end, 3] .+ sol0.u[end].DEGprotein.deg_EXG[3])/(dose_uM), 
        (tm1[end, 4] .+ sol0.u[end].DEGprotein.deg_EXG[4])/(dose_uM),
        (tm1[end, 5] .+ sol0.u[end].DEGprotein.deg_EXG[5])/(dose_uM),  
        (tm1[end, 6] .+ sol0.u[end].DEGprotein.deg_EXG[6])/(dose_uM), 
        (tm1[end, 8] .+ sol0.u[end].DEGprotein.deg_EXG[8])/(dose_uM),
        (tm1[end, 9] .+ sol0.u[end].DEGprotein.deg_EXG[9])/(dose_uM),
        (tm1[end, 10] .+ sol0.u[end].DEGprotein.deg_EXG[10])/(dose_uM) + (tm1[end, 11] .+ sol0.u[end].DEGprotein.deg_EXG[11])/(dose_uM),
        (tm1[end, 14] .+ sol0.u[end].DEGprotein.deg_EXG[14])/(dose_uM),
        (tumorm1[end] .+ dt1[end])/(dose_uM), 
    ]
);

sort!(dist_adc, [:frac], rev = true)

p_total_mass_tissue = bar(dist_adc.Organ, dist_adc.frac * 100, label = false,
    xticks =:all, xrotation = 45,
    size = [350, 250], yaxis = :log, 
    xlabel ="", ylabel = "% of the total drug", 
    ylims = (1E-3, 1E1), 
    xtickfontsize=10, 
    dpi = 1000)

savefig(p_total_mass_tissue, "img/t-dm1-total-mass-tissue.png");

#---------------------------------- Payload concentration in tissue interstitium ----------------------------------#

p_payload_dist_dm1 = plot(xlabel = "time (day)", ylabel = "payload concentration (uM)", yaxis = :log, size = (350, 300), dpi = 1000);
plot!(sol0.t/DayToHour, [sol0.u[i].end_cyto_payload[1] for i in 1:length(sol0.t)], label = "lung endothelial", alpha = 0.6, lw = 2, linestyle = :dash, color = :limegreen);
#plot!(sol0.t/DayToHour, [sol0.u[i].ints_payload[1] for i in 1:length(sol0.t)], label = "lung intersititium", alpha = 0.6, lw = 2, linestyle = :dash, color = :green4);
plot!(sol0.t/DayToHour, [sol0.u[i].end_cyto_payload[2] for i in 1:length(sol0.t)], label = "liver endothelial", alpha = 0.6, lw = 2, color = :deepskyblue);
#plot!(sol0.t/DayToHour, [sol0.u[i].ints_payload[2] for i in 1:length(sol0.t)], label = "liver intersititium", alpha = 0.6, lw = 2, color = :dodgerblue);
plot!(sol0.t/DayToHour, [sol0.u[i].end_cyto_payload[4] for i in 1:length(sol0.t)], label = "muscle endothelial", alpha = 0.6, lw = 2, linestyle = :dashdotdot, color = :blue);
#plot!(sol0.t/DayToHour, [sol0.u[i].ints_payload[4] for i in 1:length(sol0.t)], label = "muscle intersititium", alpha = 0.6, lw = 2, linestyle = :dashdotdot, color = :navyblue);
plot!(sol0.t/DayToHour, [sol0.u[i].end_cyto_payload[5] for i in 1:length(sol0.t)], label = "skin endothelial", alpha = 0.6, lw = 2, linestyle = :dashdot, color = :purple2);
#plot!(sol0.t/DayToHour, [sol0.u[i].ints_payload[5] for i in 1:length(sol0.t)], label = "skin intersititium", alpha = 0.6, lw = 2, linestyle = :dashdot, color = :purple4);
plot!(sol0.t/DayToHour, [sol0.u[i].tumor.P_c/Vc for i in 1:length(sol0.t)], label = "tumor cellular payload", alpha = 0.6, lw = 2, color = :black);
hline!([3E-3], linestyle=:dashdot, label="IC50, DM1", color = :black, alpha = 0.5); # Shim, 2020; # https://www.mdpi.com/2218-273X/10/3/360
plot!(ylims = [1E-5, 0.1], xlims = [0, 60])

savefig(p_payload_dist_dm1, "img/t-dm1-tissue-payload-conc.png");


