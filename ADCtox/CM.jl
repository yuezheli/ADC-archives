# Author: Yuezhe Li
# date: 8/3/23
# purpose: to model cantuzumab mertansine PK and payload distribution 
# note: do not use the tumor compartment in this model! The tumor part is not calibrated due to lack of information!

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Plots, Statistics, DataFrames, CSV, DataFramesMeta

include("jones_tumor_homo.jl") 
include("params.jl"); 
include("helper.jl"); 

const HT = 1.85 # [m]

# set dosing params
function adc_dosing(Dose_in_mg_m2, AddDose, infusion_time = 0.5)
    global cbs0 = [];
    if length(AddDose) >= 1
        for i in 1:length(AddDose)
            function affect_infusion_on!(integrator)
                integrator.p.infusion = Dose_in_mg_m2[i]*1E3*HT*HT/(V_Plasma)/MW_EDG/infusion_time;
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
    global cbset0 = CallbackSet(cbs0...);
    return cbset0;
end

# read in PK data
c_dm1_35 = CSV.read("data/Tolcher-2003.csv",DataFrame);

p_cm = deepcopy(p_base); 
# DAR value not changed based on https://pubmed.ncbi.nlm.nih.gov/18301896/
p_cm.PS_Score = -1.0
p_cm.PS_kd = 0.001
p_cm.init_sR = 0.001  

TotalCells = 1e10                # [cell count]

tspan = (-0.1, 500.);    # [hr]

u0, _ = jones_init(TotalCells, 1.0, 0.0, p_cm.Rcopies, p_cm.init_sR, p_cm.k_endo, p_cm.k_rec);

cbset1_cm = adc_dosing([235.], [0.]* DayToHour); 

solcm1 = solve(ODEProblem(jonesODEs_homo_tumor!, u0, tspan, p_cm), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-12, callback = cbset1_cm);

# PK dynamics 
p_pk_cm = plot(legend = :topright, ylims = [1E2, 1E6], yaxis = :log, size=(270,270), dpi = 1000);
plot!(solcm1.t, [solcm1.u[i].C_EXG_Plasma for i in 1:length(solcm1.t)]* MW_EDG, alpha = 0.6, lw = 2, color = "red", label = "235mg/m²");
scatter!(c_dm1_35.time_hour, c_dm1_35.conc_ug_L, ma = 0.6, color = "red", linewidth=0, label = false);
xlabel!("time (hour)"); ylabel!("plasma CM conc (ug/L)")

savefig(p_pk_cm, "img/cantuzumab-mertansine-pk.png");

# payload concentration 
p_payload_dist_1_cm = plot(legend = :bottom, size=(350,300), dpi = 1000);
plot!(solcm1.t/DayToHour, [solcm1.u[i].end_cyto_payload[2] for i in 1:length(solcm1.t)], label = "liver endothelial", alpha = 0.7, lw = 2.2, color = :deepskyblue);
plot!(solcm1.t/DayToHour, [solcm1.u[i].ints_payload[2] for i in 1:length(solcm1.t)], label = "liver intersititium", alpha = 0.7, lw = 2.2, color = :dodgerblue);
plot!(solcm1.t/DayToHour, [solcm1.u[i].end_cyto_payload[5] for i in 1:length(solcm1.t)], label = "skin endothelial", alpha = 0.7, lw = 2.2, linestyle = :dash, color = :purple2);
plot!(solcm1.t/DayToHour, [solcm1.u[i].ints_payload[5] for i in 1:length(solcm1.t)], label = "skin intersititium", alpha = 0.7, lw = 2.2, linestyle = :dash, color = :purple4);
hline!([3E-3], linestyle=:dashdot, label="IC50, DM1", color = :black, alpha = 0.5, lw = 1.5); # Shim, 2020; # https://www.mdpi.com/2218-273X/10/3/360
plot!(xlabel = "time (day)", ylabel = "payload concentration (uM)", ylims = [1E-4, 0.1], yaxis = :log)

savefig(p_payload_dist_1_cm, "img/tissue-payload-cm-235mgm2.png");

# add additional payload analysis, with a focus on liver payload 
# dose obtained from Tolcher et al., 2003; # https://pubmed.ncbi.nlm.nih.gov/12525512/
solcm2 = solve(ODEProblem(jonesODEs_homo_tumor!, u0, tspan, p_cm), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-12, callback = adc_dosing([295.], [0.]* DayToHour); );
solcm3 = solve(ODEProblem(jonesODEs_homo_tumor!, u0, tspan, p_cm), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-12, callback = adc_dosing([176.], [0.]* DayToHour); );
solcm4 = solve(ODEProblem(jonesODEs_homo_tumor!, u0, tspan, p_cm), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-12, callback = adc_dosing([132.], [0.]* DayToHour); );
solcm5 = solve(ODEProblem(jonesODEs_homo_tumor!, u0, tspan, p_cm), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-12, callback = adc_dosing([88.], [0.]* DayToHour); );
solcm6 = solve(ODEProblem(jonesODEs_homo_tumor!, u0, tspan, p_cm), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-12, callback = adc_dosing([44.], [0.]* DayToHour); );
solcm7 = solve(ODEProblem(jonesODEs_homo_tumor!, u0, tspan, p_cm), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-12, callback = adc_dosing([22.], [0.]* DayToHour); );

# Cantuzumab Mertansine in Liver Interstitium
p_payload_dist_2_cm = plot(legend = :top, size=(350,300), legendcolumns=3, dpi = 1000);
plot!(solcm2.t/DayToHour, [solcm2.u[i].ints_payload[2] for i in 1:length(solcm2.t)], label = "295mg/m²", alpha = 0.8, lw = 2.2);
plot!(solcm1.t/DayToHour, [solcm1.u[i].ints_payload[2] for i in 1:length(solcm1.t)], label = "235mg/m²", alpha = 0.8, lw = 2.2);
plot!(solcm3.t/DayToHour, [solcm3.u[i].ints_payload[2] for i in 1:length(solcm3.t)], label = "176mg/m²", alpha = 0.8, lw = 2.2);
plot!(solcm4.t/DayToHour, [solcm4.u[i].ints_payload[2] for i in 1:length(solcm4.t)], label = "132mg/m²", alpha = 0.8, lw = 2.2);
plot!(solcm5.t/DayToHour, [solcm5.u[i].ints_payload[2] for i in 1:length(solcm5.t)], label = "88mg/m²", alpha = 0.8, lw = 2.2);
plot!(solcm6.t/DayToHour, [solcm6.u[i].ints_payload[2] for i in 1:length(solcm6.t)], label = "44mg/m²", alpha = 0.8, lw = 2.2);
plot!(solcm7.t/DayToHour, [solcm7.u[i].ints_payload[2] for i in 1:length(solcm7.t)], label = "22mg/m²", alpha = 0.8, lw = 2.2);
hline!([3E-3], linestyle=:dashdot, label="IC50, DM1", color = :black, alpha = 0.5, lw = 2); # Shim, 2020; # https://www.mdpi.com/2218-273X/10/3/360
plot!(xlabel = "time (day)", ylabel = "payload concentration (uM)", ylims = [1E-4, 1], yaxis = :log)

# Cantuzumab Mertansine in Skin Interstitium
p_payload_dist_3_cm = plot(legend = :outerright, legendtitle = "Dose", legendtitlefontsize = 8);
plot!(solcm2.t/DayToHour, [solcm2.u[i].ints_payload[5] for i in 1:length(solcm2.t)], label = "295mg/m²", alpha = 0.6, lw = 2);
plot!(solcm1.t/DayToHour, [solcm1.u[i].ints_payload[5] for i in 1:length(solcm1.t)], label = "235mg/m²", alpha = 0.6, lw = 2);
plot!(solcm3.t/DayToHour, [solcm3.u[i].ints_payload[5] for i in 1:length(solcm3.t)], label = "176mg/m²", alpha = 0.6, lw = 2);
plot!(solcm4.t/DayToHour, [solcm4.u[i].ints_payload[5] for i in 1:length(solcm4.t)], label = "132mg/m²", alpha = 0.6, lw = 2);
plot!(solcm5.t/DayToHour, [solcm5.u[i].ints_payload[5] for i in 1:length(solcm5.t)], label = "88mg/m²", alpha = 0.6, lw = 2);
plot!(solcm6.t/DayToHour, [solcm6.u[i].ints_payload[5] for i in 1:length(solcm6.t)], label = "44mg/m²", alpha = 0.6, lw = 2);
plot!(solcm7.t/DayToHour, [solcm7.u[i].ints_payload[5] for i in 1:length(solcm7.t)], label = "22mg/m²", alpha = 0.6, lw = 2);
hline!([3E-3], linestyle=:dashdot, label="DM-1 IC50", color = :black, alpha = 0.5); # Shim, 2020; # https://www.mdpi.com/2218-273X/10/3/360
plot!(xlabel = "time (day)", ylabel = "payload concentration (uM)", ylims = [1E-4, 0.1], yaxis = :log)

savefig(p_payload_dist_2_cm, "img/tissue-payload-cm-liver-init.png");
savefig(p_payload_dist_3_cm, "img/tissue-payload-cm-skin-init.png");

