# date: Apr 18, 2024
# author: Yuezhe Li 
# purpose of this script: to test PK of sacituzumab govitecan, 
# observed range from Supp Table 1 of Starodub et al., 2016; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4558321/

using Pkg; Pkg.activate("../");

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using DataFrames, DataFramesMeta, CSV
using Plots

# read in data from Santi et al., 2016; https://pubmed.ncbi.nlm.nih.gov/38236523/
pk_santi =  CSV.read("../data/Santi2024.csv",DataFrame);
pk_intact = @rsubset(pk_santi, :type == "intact_Trodelvy");
pk_total = @rsubset(pk_santi, :type == "total_Trodelvy");

include("mPBPK_deconj_homo.jl");
include("init.jl");
include("param.jl");

TotalCells = 1E-3/Vc; # assume tumor volume = 1mL

function infusion_sims(dose_mgkg, infusion_time = 2, sampling_time = 0.5, p_immu_132 = p_immu_132)
    tspan = (0., infusion_time + sampling_time);    # [hr]
    p_tmp = deepcopy(p_immu_132);
    u0_mpbpk, _ = pbpk_init(7, TotalCells, 0., 0, p_tmp.Rcopies, p_tmp.init_sR, p_tmp.k_endo, p_tmp.k_rec); 
    p_tmp.infusion = dose_mgkg*BW*1E3/(V_Plasma)/MW_EDG/infusion_time;
    function affect_infusion_off!(integrator)
        integrator.p.infusion = 0.;
    end
    cb02 = PresetTimeCallback(infusion_time, affect_infusion_off!);
    # simulation 
    sol_mpbpk = solve(ODEProblem(mPBPK_deconj_homo!, u0_mpbpk, tspan, p_tmp), saveat = 0.5, alg = QNDF(autodiff=false), reltol = 1e-12, callback = cb02);
    # return the final plasma conc 
    return sol_mpbpk.u[end].C_EXG_Plasma * MW_EDG * 1E-3  # (ug/mL)
end

# update params 
p_immu_132 = ComponentArray(p_base, k_dec = 0.); 
p_immu_132.k_dec = log(2)/DayToHour; # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4673178/
p_immu_132.DAR = 7.6; # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4558321/
p_immu_132.Rcopies = 3.0E5  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4673178/
p_immu_132.init_sR = 0.  # assumed 
p_immu_132.PS_Score = 8.

p_immu_132.Kd = 0.56E-3         # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4673178/
p_immu_132.ic50_pl = 2.0E-3  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4673178/
p_immu_132.k_kill_max = 0.032
# updated based on dato-dxd internalization data
p_immu_132.k_deg = 6.22
p_immu_132.k_rec = 0.594
p_immu_132.k_endo = 0.575
# updated based on sn38 data 
p_immu_132.k_in = 11.
p_immu_132.tau = 18.16 
p_immu_132.k_in = 27.6

dose_mgkg = [8., 10., 12., 18.]; # [mg/kg]
mean_adc = [141.5, 185.77, 183.3, 258.27]; # [ug/mL]; Supp Table 1, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4558321/
# sd for mean_adc [23.8, 54.14, 71.8, 143.26]
mean_IgG = [193.1, 203.35, 239.2, 409.16]; # conc of hRS7 [ug/mL]; Supp Table 1, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4558321/
# sd for mean_IgG [56.5, 55.72, 70.7, 88,78]
pred_IgG = [];
pred_adc = [];
for i in 1:length(dose_mgkg)
    tmp_p_immu_132 = deepcopy(p_immu_132); 
    tmp_p_immu_132.k_dec = 0.
    append!(pred_IgG, infusion_sims(dose_mgkg[i], 2., 0.5, tmp_p_immu_132)); 
    append!(pred_adc, infusion_sims(dose_mgkg[i], 2., 0.5, p_immu_132)); 
end

pred_adc
pred_IgG


# compare PK from Santi et al., 2024; https://pubmed.ncbi.nlm.nih.gov/38236523/
function double_sims(dose_mgkg, simstime = 7*DayToHour, dosetimes = [7, 21, 28, 42, 49]*DayToHour, p_immu_132 = p_immu_132)
    tspan = (0., simstime);    # [hr]
    p_tmp = deepcopy(p_immu_132);
    # initial bolus IV dose
    u0_mpbpk, _ = pbpk_init(7, TotalCells, 0., 0, p_tmp.Rcopies, p_tmp.init_sR, p_tmp.k_endo, p_tmp.k_rec); 
    u0_mpbpk.C_EXG_Plasma = dose_mgkg*BW*1E3/(V_Plasma)/MW_EDG;
    # define additional dosing 
    condition(u, t, integrator) = t ∈ dosetimes
    affect!(integrator) = integrator.u.C_EXG_Plasma += dose_mgkg*BW*1E3/(V_Plasma)/MW_EDG
    cb = DiscreteCallback(condition, affect!)
    # simulation for ADC 
    sol_mpbpk = solve(ODEProblem(mPBPK_deconj_homo!, u0_mpbpk, tspan, p_tmp), saveat = 3, alg = QNDF(autodiff=false), reltol = 1e-12, callback = cb, tstops = dosetimes);
    c_plasma_adc = [sol_mpbpk.u[i].C_EXG_Plasma for i in 1:length(sol_mpbpk.t)]; # [uM]
    # simulation for total ADC 
    p2_tmp = deepcopy(p_immu_132); p2_tmp.k_dec = 0.; 
    sol2_mpbpk = solve(ODEProblem(mPBPK_deconj_homo!, u0_mpbpk, tspan, p2_tmp), saveat = 3, alg = QNDF(autodiff=false), reltol = 1e-12, callback = cb, tstops = dosetimes);
    c_plasma_total = [sol2_mpbpk.u[i].C_EXG_Plasma for i in 1:length(sol2_mpbpk.t)]; # [uM]
    # return the final plasma conc 
    return sol_mpbpk.t/DayToHour, c_plasma_adc*1E3, c_plasma_total*1E3; # [nM]
end

t, intact, total = double_sims(10., 63*DayToHour);

p_pk = plot(title = "PK", xlabel = "Time (Day)", ylabel = "Conc (nM)");
plot!(t, intact, label = "intact ADC", color = :red);
scatter!(pk_intact.time_day, pk_intact.conc_nM, label = false, color = :red, alpha = 0.7);
plot!(t, total, label = "total ADC", color = :blue);
scatter!(pk_total.time_day, pk_total.conc_nM, label = false, color = :blue, alpha = 0.7);
display(p_pk);

savefig(p_pk, "../figure/pk/immu-132.png");

