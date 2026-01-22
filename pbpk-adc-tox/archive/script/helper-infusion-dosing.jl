# date: 6/25/2025
# author: Yuezhe Li 
# purpose of this code: for IV dosing

function InfusionDoses(init_dose, AddDose, p_base = p_base; infusion_time = 0.5, tspan = tspan, Reltol = 1E-12, BW = BW)
    # initialization 
    u0_tmp = jones_init(0., p_base, BW, V_Plasma);

    # repeated dosing 
    global cbs0 = [];
    if length(AddDose) >= 1
        for i in 1:length(AddDose)
            function affect_infusion_on!(integrator)
                integrator.p.infusion = init_dose *1E3*BW/(V_Plasma)/MW_IGG/infusion_time;
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
    sol_tmp = solve(ODEProblem(jonesODEs_homo_tumor!, u0_tmp, tspan, p_base), saveat = 0.25, alg = QNDF(autodiff=false), reltol = Reltol, callback = cbset0);

    return sol_tmp
end
