# Author: Yuezhe Li
# date: 8/2/23
# purpose: to model Brentuximab Vedotin PK and payload distribution 

using Pkg; Pkg.activate("");

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Plots, Statistics, DataFrames, CSV, DataFramesMeta

include("jones_tumor_homo.jl") 
include("params.jl"); 
include("helper.jl"); 

# set dosing params
function sgn_35_dosing(Dose_in_ugkg, AddDose, infusion_time = 0.5)
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
    global cbset0 = CallbackSet(cbs0...);
    return cbset0;
end

# read in PK data
sgn_35 = CSV.read("data/Younes-2010.csv",DataFrame);

pk_1point2 = @rsubset(sgn_35, :Dose_mg_kg == 1.2); # used in clinics
pk_1point8 = @rsubset(sgn_35, :Dose_mg_kg == 1.8); # used in clinics
pk_2point7 = @rsubset(sgn_35, :Dose_mg_kg == 2.7); # not used in clinical 

p_bv = deepcopy(p_base); 
# PK parameter for brentuximab vedotin (bv) (tuned based on observed data)
p_bv.PS_Score = 10.
p_bv.init_sR = 0.01  # https://www.nature.com/articles/bcj201785
p_bv.thalf_sR_adc = 120.
# update parameter for brentuximab vedotin (bv)
p_bv.Rcopies = 9.3E5  # https://pubmed.ncbi.nlm.nih.gov/22542449/
p_bv.P_ADC = 334. / 24   # rate of permeability [um/h]; Shah et al., 2012 # https://pubmed.ncbi.nlm.nih.gov/23151991/
p_bv.D_ADC = 0.022 / 24  # rate of diffusion [cm^2/h]; Shah et al., 2012 # https://pubmed.ncbi.nlm.nih.gov/23151991/
p_bv.Kd = 2.0E-3         # binding affinity between CD30 and brentuximab, [uM], Shah et al., 2012 # https://pubmed.ncbi.nlm.nih.gov/23151991/
p_bv.k_endo = 0.61/24    # CD30:ADC internalization rate 
p_bv.k_out = 2.47        # https://bmccancer.biomedcentral.com/articles/10.1186/s12885-019-5336-7
p_bv.k_in = 5.07         # https://bmccancer.biomedcentral.com/articles/10.1186/s12885-019-5336-7
p_bv.IC50_Payload = 0.97E-3 # [uM] # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5024081/
p_bv.DAR = 4.4           # https://pubmed.ncbi.nlm.nih.gov/23151991/

TotalCells = 1.13e9              # [cell count]

tspan = (-0.1, DayToHour*84);    # [hr]

u0_bv, _ = jones_init(TotalCells, 1.0, 0.0, p_bv.Rcopies, p_bv.init_sR, p_bv.k_endo, p_bv.k_rec);

cbset1_bv = sgn_35_dosing([1.2E3, 1.2E3], [0., 21.]* DayToHour); 
cbset2_bv = sgn_35_dosing([1.8E3, 1.8E3], [0., 21.]* DayToHour); 
cbset3_bv = sgn_35_dosing([2.7E3, 2.7E3], [0., 21.]* DayToHour); 

solbv1 = solve(ODEProblem(jonesODEs_homo_tumor!, u0_bv, tspan, p_bv), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-12, callback = cbset1_bv);
solbv2 = solve(ODEProblem(jonesODEs_homo_tumor!, u0_bv, tspan, p_bv), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-12, callback = cbset2_bv);
solbv3 = solve(ODEProblem(jonesODEs_homo_tumor!, u0_bv, tspan, p_bv), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-12, callback = cbset3_bv);

# PK dynamics 
p_pk_bv = plot(legend = :outertop, ylims = [1E-1, 1E2], yaxis = :log, size=(300,300), legendcolumns=3, dpi = 1000);
plot!(solbv1.t/DayToHour, [solbv1.u[i].C_EXG_Plasma for i in 1:length(solbv1.t)]* MW_EDG / 1E3, alpha = 0.6, color = "red", label = "1.2mg/kg", lw=2);
scatter!(pk_1point2.time_day, pk_1point2.ADC_ug_mL, ma = 0.6, color = "red", label = false, linewidth=0);
plot!(solbv2.t/DayToHour, [solbv2.u[i].C_EXG_Plasma for i in 1:length(solbv2.t)]* MW_EDG / 1E3, alpha = 0.6, color = "seagreen", label = "1.8mg/kg", lw=2);
scatter!(pk_1point8.time_day, pk_1point8.ADC_ug_mL, ma = 0.6, color = "seagreen", label = false, linewidth=0);
plot!(solbv3.t/DayToHour, [solbv3.u[i].C_EXG_Plasma for i in 1:length(solbv3.t)]* MW_EDG / 1E3, alpha = 0.6, color = "mediumpurple", label = "2.7mg/kg", lw=2);
scatter!(pk_2point7.time_day, pk_2point7.ADC_ug_mL, ma = 0.6, color = "mediumpurple", label = false, linewidth=0);
xlabel!("time (day)"); ylabel!("plasma BV conc (ug/mL)"); xlims!(0, 64)

savefig(p_pk_bv, "img/brentuximab-vedotin-pk.png");

# payload concentration (1.2mg/kg)
p_payload_dist_1 = plot(legend = :right, size=(350,300), dpi = 1000);
plot!(solbv1.t/DayToHour, [solbv1.u[i].end_cyto_payload[1] for i in 1:length(solbv1.t)], label = "lung endothelial", alpha = 0.6, lw = 2, linestyle = :dash, color = :limegreen);
#plot!(solbv1.t/DayToHour, [solbv1.u[i].ints_payload[1] for i in 1:length(solbv1.t)], label = "lung intersititium", alpha = 0.6, lw = 2, linestyle = :dash, color = :green4);
plot!(solbv1.t/DayToHour, [solbv1.u[i].end_cyto_payload[2] for i in 1:length(solbv1.t)], label = "liver endothelial", alpha = 0.6, lw = 2, color = :deepskyblue);
#plot!(solbv1.t/DayToHour, [solbv1.u[i].ints_payload[2] for i in 1:length(solbv1.t)], label = "liver intersititium", alpha = 0.6, lw = 2, color = :dodgerblue);
plot!(solbv1.t/DayToHour, [solbv1.u[i].end_cyto_payload[5] for i in 1:length(solbv1.t)], label = "skin endothelial", alpha = 0.6, lw = 2, linestyle = :dashdot, color = :purple2);
#plot!(solbv1.t/DayToHour, [solbv1.u[i].ints_payload[5] for i in 1:length(solbv1.t)], label = "skin intersititium", alpha = 0.6, lw = 2, linestyle = :dashdot, color = :purple4);
#plot!(solbv1.t/DayToHour, [solbv1.u[i].tumor.P_c/Vc for i in 1:length(solbv1.t)], label = "tumor cellular payload", alpha = 0.6, lw = 2, color = :black);
hline!([3E-3], linestyle=:dashdot, label="IC50, MMAE", color = :black, alpha = 0.5); # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5024081/
plot!(xlabel = "time (day)", ylabel = "payload concentration (uM)", ylims = [1E-5, 1E-1], xlims = [0, 21], yaxis = :log)

# payload concentration (1.8mg/kg)
p_payload_dist_2 = plot(title = "brentuximab vedotin dose = 1.8mg/kg", titlefontsize = 8, legend = :outerright, dpi = 1000);
plot!(solbv2.t/DayToHour, [solbv2.u[i].end_cyto_payload[1] for i in 1:length(solbv1.t)], label = "lung endothelial", alpha = 0.6, lw = 2, linestyle = :dash, color = :limegreen);
#plot!(solbv2.t/DayToHour, [solbv2.u[i].ints_payload[1] for i in 1:length(solbv1.t)], label = "lung intersititium", alpha = 0.6, lw = 2, linestyle = :dash, color = :green4);
plot!(solbv2.t/DayToHour, [solbv2.u[i].end_cyto_payload[2] for i in 1:length(solbv1.t)], label = "liver endothelial", alpha = 0.6, lw = 2, color = :deepskyblue);
#plot!(solbv2.t/DayToHour, [solbv2.u[i].ints_payload[2] for i in 1:length(solbv1.t)], label = "liver intersititium", alpha = 0.6, lw = 2, color = :dodgerblue);
plot!(solbv2.t/DayToHour, [solbv2.u[i].end_cyto_payload[5] for i in 1:length(solbv1.t)], label = "skin endothelial", alpha = 0.6, lw = 2, linestyle = :dashdot, color = :purple2);
#plot!(solbv2.t/DayToHour, [solbv2.u[i].ints_payload[5] for i in 1:length(solbv1.t)], label = "skin intersititium", alpha = 0.6, lw = 2, linestyle = :dashdot, color = :purple4);
#plot!(solbv2.t/DayToHour, [solbv2.u[i].tumor.P_c/Vc for i in 1:length(solbv1.t)], label = "tumor cellular payload", alpha = 0.6, lw = 2, color = :black);
hline!([0.97E-3], linestyle=:dashdot, label=false, color = :black, alpha = 0.5); # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5024081/
plot!(xlabel = "time (day)", ylabel = "payload concentration (uM)", ylims = [1E-5, 0.1], xlims = [0, 21], yaxis = :log)

savefig(p_payload_dist_1, "img/tissue-payload-bv-1.2mg.png");
savefig(p_payload_dist_2, "img/tissue-payload-bv-1.8mg.png");
