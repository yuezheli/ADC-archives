# author: Yuezhe Li 
# date: 7/11/23
# purpose: to adjust PK for T-DXd 

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Plots, Statistics, DataFrames, CSV, DataFramesMeta

# read observed 
pk_doi =  CSV.read("data/doi-2017.csv",DataFrame);

pk_point8 = @rsubset(pk_doi, :Dose .== "0.8mg/kg" ); 
pk_1point6 = @rsubset(pk_doi, :Dose .== "1.6mg/kg" ); 
pk_3point2 = @rsubset(pk_doi, :Dose .== "3.2mg/kg" ); 
pk_5point4 = @rsubset(pk_doi, :Dose .== "5.4mg/kg" ); 
pk_8 = @rsubset(pk_doi, :Dose .== "8.0mg/kg" ); 


# simulation 
include("jones_tumor_homo.jl") 
include("params.jl"); 
include("helper.jl"); 

TotalCells = 1.13e9       # [cell count]; convert to 2mL of initial tumor size

p_dxd = deepcopy(p_base);
p_dxd.Rcopies = 1.0E6
p_dxd.PS_Score = 6.
p_dxd.init_sR = 0.004
p_dxd.thalf_sR_adc = 120.
p_dxd.DAR = 8.
p_dxd.k_in = 46.08
p_dxd.k_out = 32.32
p_dxd.Emax_Payload = 0.201
p_dxd.IC50_Payload = 9.54E-3
p_dxd.k_PL = 0.82

tspan = (-0.01, DayToHour*63);      # [hr]
AddDose = [0., 22., 44.] * DayToHour  # [hr]

first_dose = [0.8, 1.6, 3.2, 5.4, 8.0]; # [mg/kg]

function MultiDosesSims(init_dose, AddDose, p_base = p_base, infusion_time = 0.5)
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

    #affect!(integrator) = integrator.u[1] += u0_tmp.C_EXG_Plasma;  cb = PresetTimeCallback(AddDose,affect!);

    sol_tmp = solve(ODEProblem(jonesODEs_homo_tumor!, u0_tmp, tspan, p_base), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-12, callback = cbset0);
    cplasma_tmp = [sol_tmp.u[i].C_EXG_Plasma for i in 1:length(sol_tmp.t)]; # [uM]

    return sol_tmp.t, cplasma_tmp
end

timestamp = []; 
cplasma = [];

for init_dose in first_dose
    try 
        tmpt, tmpp = MultiDosesSims(init_dose, AddDose, p_dxd, 0.5);
        append!(timestamp, [tmpt/DayToHour]);
        append!(cplasma, [tmpp * MW_EDG / 1E3]); # [ug/mL]
    catch 
        println("stability issue at dose " *string(init_dose) * " mg/kg");
    end
end

p_pk_tdxd = plot(legend = :top, ylims = [1E-1, 1E4], yaxis = :log, size=(350,350), legendcolumns=3, dpi = 1000);
plot!(timestamp[1], cplasma[1], lw = 2, alpha = 0.6, color = "red", label = "0.8mg/kg");
scatter!(pk_point8.time_day, pk_point8.T_DXd_ugperml, ma = 0.6, color = "red", label = false, linewidth=0);
plot!(timestamp[2], cplasma[2], lw = 2, alpha = 0.6, color = "seagreen", label = "1.6mg/kg");
scatter!(pk_1point6.time_day, pk_1point6.T_DXd_ugperml, ma = 0.6, color = "seagreen", label = false, linewidth=0);
plot!(timestamp[3], cplasma[3], lw = 2, alpha = 0.6, color = "mediumpurple", label = "3.2mg/kg");
scatter!(pk_3point2.time_day, pk_3point2.T_DXd_ugperml, ma = 0.6, color = "mediumpurple", label = false, linewidth=0);
plot!(timestamp[4], cplasma[4], lw = 2, alpha = 0.6, color = "violetred", label = "5.4mg/kg");
scatter!(pk_5point4.time_day, pk_5point4.T_DXd_ugperml, ma = 0.6, color = "violetred", label = false, linewidth=0);
plot!(timestamp[5], cplasma[5], lw = 2, alpha = 0.6, color = "blue2", label = "8mg/kg");
scatter!(pk_8.time_day, pk_8.T_DXd_ugperml, ma = 0.6, color = "blue2", label = false, linewidth=0);
xlabel!("time (day)"); xlims!(0, 63); ylabel!("plasma T-Dxd conc (ug/mL)")

savefig(p_pk_tdxd, "img/t-dxd-pk.png");

#---------------------------------- Payload concentration in tissue interstitium ----------------------------------#
AddDose = [0.] * DayToHour        # [hr]
Dose_in_ugkg = ones(size(AddDose)) * 5.4E3  # [ug/kg]    # https://www.nejm.org/doi/full/10.1056/NEJMoa2203690
infusion_time = 0.5                         # [hr]

u1, _ = jones_init(TotalCells, 1.0, 0.0, p_dxd.Rcopies, p_dxd.init_sR, p_dxd.k_endo, p_dxd.k_rec);

dose_uM = sum(Dose_in_ugkg)*BW/MW_EDG

global cbs1 = [];
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
        global cbs0 = push!(cbs1, cb01);
        global cbs0 = push!(cbs1, cb02);
    end
end
cbset1 = CallbackSet(cbs1...);

sol2 = solve(ODEProblem(jonesODEs_homo_tumor!, u0, tspan, p_dxd), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-12, callback = cbset1);

p_payload_dist_dxd = plot(xlabel = "time (day)", ylabel = "payload concentration (uM)", yaxis = :log, size = (350, 300), dpi = 1000);
plot!(sol2.t/DayToHour, [sol2.u[i].end_cyto_payload[1] for i in 1:length(sol2.t)], label = "lung endothelial", alpha = 0.6, lw = 2, linestyle = :dash, color = :limegreen);
#plot!(sol2.t/DayToHour, [sol2.u[i].ints_payload[1] for i in 1:length(sol2.t)], label = "lung intersititium", alpha = 0.6, lw = 2, linestyle = :dash, color = :green4);
plot!(sol2.t/DayToHour, [sol2.u[i].end_cyto_payload[2] for i in 1:length(sol2.t)], label = "liver endothelial", alpha = 0.6, lw = 2, color = :deepskyblue);
#plot!(sol2.t/DayToHour, [sol2.u[i].ints_payload[2] for i in 1:length(sol2.t)], label = "liver intersititium", alpha = 0.6, lw = 2, color = :dodgerblue);
plot!(sol2.t/DayToHour, [sol2.u[i].end_cyto_payload[4] for i in 1:length(sol2.t)], label = "muscle endothelial", alpha = 0.6, lw = 2, linestyle = :dashdotdot, color = :blue);
#plot!(sol2.t/DayToHour, [sol2.u[i].ints_payload[4] for i in 1:length(sol2.t)], label = "muscle intersititium", alpha = 0.6, lw = 2, linestyle = :dashdotdot, color = :navyblue);
plot!(sol2.t/DayToHour, [sol2.u[i].end_cyto_payload[5] for i in 1:length(sol2.t)], label = "skin endothelial", alpha = 0.6, lw = 2, linestyle = :dashdot, color = :purple2);
#plot!(sol2.t/DayToHour, [sol2.u[i].ints_payload[5] for i in 1:length(sol2.t)], label = "skin intersititium", alpha = 0.6, lw = 2, linestyle = :dashdot, color = :purple4);
plot!(sol2.t/DayToHour, [sol2.u[i].tumor.P_c/Vc for i in 1:length(sol2.t)], label = "tumor cellular payload", alpha = 0.6, lw = 2, color = :black);
hline!([0.357E-3], linestyle=:dashdot, label="IC50, Dxd", color = :black, alpha = 0.5); # https://aacrjournals.org/mct/article/21/4/635/689570/DS-7300a-a-DNA-Topoisomerase-I-Inhibitor-DXd-Based
plot!(ylims = [1E-5, 0.1], xlims = [0, 60])

savefig(p_payload_dist_dxd, "img/t-dxd-tissue-payload-conc.png");

