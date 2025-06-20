# date: 1/15/2025
# author: Yuezhe Li 
# purpose of this script: to house commonly used function for simulation of endothelial payload 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using DataFrames, DataFramesMeta

include(@projectroot("model/adc-constants.jl"));
include(@projectroot("model/jones_homo.jl") );
include(@projectroot("model/param-pk.jl")); 
include(@projectroot("model/init-pk.jl")); 

# repeated, infusion dosing simulation
function InfusionDoses(init_dose, AddDose, p_base = p_base; infusion_time = 0.5, tspan = tspan, Reltol = 1E-12, BW = BW)
    # initialization 
    u0_tmp, _ = jones_init(0., p_base);

    # repeated dosing 
    global cbs0 = [];
    if length(AddDose) >= 1
        for i in 1:length(AddDose)
            function affect_infusion_on!(integrator)
                integrator.p.infusion = init_dose *1E3*BW/(V_Plasma)/MW/infusion_time;
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

    # simulation 
    sol_tmp = solve(ODEProblem(jonesODEs_homo_tumor!, u0_tmp, tspan, p_base), saveat = 1., alg = QNDF(autodiff=false), reltol = Reltol, callback = cbset0);

    
    tmpdf = DataFrame(time = sol_tmp.t)

    # process the endothelial payload concentration 
    tmpdf.pl_lung = [sol_tmp.u[i].end_cyto_payload[1] for i in 1:length(sol_tmp.t)]; # [uM]
    tmpdf.pl_liver = [sol_tmp.u[i].end_cyto_payload[2] for i in 1:length(sol_tmp.t)]; # [uM]
    tmpdf.pl_heart = [sol_tmp.u[i].end_cyto_payload[3] for i in 1:length(sol_tmp.t)]; # [uM]
    tmpdf.pl_muscle = [sol_tmp.u[i].end_cyto_payload[4] for i in 1:length(sol_tmp.t)]; # [uM]
    tmpdf.pl_skin = [sol_tmp.u[i].end_cyto_payload[5] for i in 1:length(sol_tmp.t)]; # [uM]
    tmpdf.pl_adipose = [sol_tmp.u[i].end_cyto_payload[6] for i in 1:length(sol_tmp.t)]; # [uM]
    tmpdf.pl_bone = [sol_tmp.u[i].end_cyto_payload[7] for i in 1:length(sol_tmp.t)]; # [uM]
    tmpdf.pl_brain = [sol_tmp.u[i].end_cyto_payload[8] for i in 1:length(sol_tmp.t)]; # [uM]
    tmpdf.pl_kidney = [sol_tmp.u[i].end_cyto_payload[9] for i in 1:length(sol_tmp.t)]; # [uM]
    tmpdf.pl_si = [sol_tmp.u[i].end_cyto_payload[10] for i in 1:length(sol_tmp.t)]; # [uM]
    tmpdf.pl_li = [sol_tmp.u[i].end_cyto_payload[11] for i in 1:length(sol_tmp.t)]; # [uM]
    tmpdf.pl_pancreas = [sol_tmp.u[i].end_cyto_payload[12] for i in 1:length(sol_tmp.t)]; # [uM]
    tmpdf.pl_thymus = [sol_tmp.u[i].end_cyto_payload[13] for i in 1:length(sol_tmp.t)]; # [uM]
    tmpdf.pl_spleen = [sol_tmp.u[i].end_cyto_payload[14] for i in 1:length(sol_tmp.t)]; # [uM]

    # plasma 
    tmpdf.adcplasma = [sol_tmp.u[i].C_EXG_Plasma for i in 1:length(sol_tmp.t)]; # [uM]

    @rsubset!(tmpdf, :time >=0)

    return tmpdf 
end

# compute Ctrough 
function Ctrough(df_pl, time_d, organ)
    organindex = findall( (x -> x == "pl_"* organ), names(df_pl)) ;
    tmpdf = df_pl[:,[1,organindex[1]]];

    tmp = []
    for timd_d_tmp in time_d
        tmp2 = @rsubset(tmpdf, :time == timd_d_tmp * DayToHour)
        append!(tmp, tmp2[1, 2]);
    end

    return tmp
end

# compute the ratio of payload concentration in tissue interstitium and in tissue endothelial cells 
function tissue_ints_endo(sol_tmp)
    # process the endothelial payload concentration 
    endo_pl_lung = [sol_tmp.u[i].end_cyto_payload[1] for i in 1:length(sol_tmp.t)]; # [uM]
    endo_pl_liver = [sol_tmp.u[i].end_cyto_payload[2] for i in 1:length(sol_tmp.t)]; # [uM]
    endo_pl_heart = [sol_tmp.u[i].end_cyto_payload[3] for i in 1:length(sol_tmp.t)]; # [uM]
    endo_pl_muscle = [sol_tmp.u[i].end_cyto_payload[4] for i in 1:length(sol_tmp.t)]; # [uM]
    endo_pl_skin = [sol_tmp.u[i].end_cyto_payload[5] for i in 1:length(sol_tmp.t)]; # [uM]
    endo_pl_adipose = [sol_tmp.u[i].end_cyto_payload[6] for i in 1:length(sol_tmp.t)]; # [uM]
    endo_pl_bone = [sol_tmp.u[i].end_cyto_payload[7] for i in 1:length(sol_tmp.t)]; # [uM]
    endo_pl_brain = [sol_tmp.u[i].end_cyto_payload[8] for i in 1:length(sol_tmp.t)]; # [uM]
    endo_pl_kidney = [sol_tmp.u[i].end_cyto_payload[9] for i in 1:length(sol_tmp.t)]; # [uM]
    endo_pl_si = [sol_tmp.u[i].end_cyto_payload[10] for i in 1:length(sol_tmp.t)]; # [uM]
    endo_pl_li = [sol_tmp.u[i].end_cyto_payload[11] for i in 1:length(sol_tmp.t)]; # [uM]
    endo_pl_pancreas = [sol_tmp.u[i].end_cyto_payload[12] for i in 1:length(sol_tmp.t)]; # [uM]
    endo_pl_thymus = [sol_tmp.u[i].end_cyto_payload[13] for i in 1:length(sol_tmp.t)]; # [uM]
    endo_pl_spleen = [sol_tmp.u[i].end_cyto_payload[14] for i in 1:length(sol_tmp.t)]; # [uM]

    # process payload concentratio in tissue interstitium 
    ints_pl_lung = [sol_tmp.u[i].ints_payload[1] for i in 1:length(sol_tmp.t)]; # [uM]
    ints_pl_liver = [sol_tmp.u[i].ints_payload[2] for i in 1:length(sol_tmp.t)]; # [uM]
    ints_pl_heart = [sol_tmp.u[i].ints_payload[3] for i in 1:length(sol_tmp.t)]; # [uM]
    ints_pl_muscle = [sol_tmp.u[i].ints_payload[4] for i in 1:length(sol_tmp.t)]; # [uM]
    ints_pl_skin = [sol_tmp.u[i].ints_payload[5] for i in 1:length(sol_tmp.t)]; # [uM]
    ints_pl_adipose = [sol_tmp.u[i].ints_payload[6] for i in 1:length(sol_tmp.t)]; # [uM]
    ints_pl_bone = [sol_tmp.u[i].ints_payload[7] for i in 1:length(sol_tmp.t)]; # [uM]
    ints_pl_brain = [sol_tmp.u[i].ints_payload[8] for i in 1:length(sol_tmp.t)]; # [uM]
    ints_pl_kidney = [sol_tmp.u[i].ints_payload[9] for i in 1:length(sol_tmp.t)]; # [uM]
    ints_pl_si = [sol_tmp.u[i].ints_payload[10] for i in 1:length(sol_tmp.t)]; # [uM]
    ints_pl_li = [sol_tmp.u[i].ints_payload[11] for i in 1:length(sol_tmp.t)]; # [uM]
    ints_pl_pancreas = [sol_tmp.u[i].ints_payload[12] for i in 1:length(sol_tmp.t)]; # [uM]
    ints_pl_thymus = [sol_tmp.u[i].ints_payload[13] for i in 1:length(sol_tmp.t)]; # [uM]
    ints_pl_spleen = [sol_tmp.u[i].ints_payload[14] for i in 1:length(sol_tmp.t)]; # [uM]

    tmpdf = DataFrame(time = sol_tmp.t)
    tmpdf.r_he = endo_pl_liver./ints_pl_liver
    tmpdf.r_skin = endo_pl_skin ./ ints_pl_skin
    tmpdf.r_si = endo_pl_si ./ ints_pl_si
    tmpdf.endo_he = endo_pl_liver
    tmpdf.endo_skin = endo_pl_skin

    tmpdf.adcplasma = [sol_tmp.u[i].C_EXG_Plasma for i in 1:length(sol_tmp.t)]; # [uM]

    return tmpdf
end