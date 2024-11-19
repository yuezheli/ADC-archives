# date: Apr 22, 2024
# author: Yuezhe Li 
# purpose of this script: to fit the PK of Datopotamab deruxtecan
# observed human data were obtained from Tang poster from AACR 2024 (San Diego, CA)

using Pkg; Pkg.activate("../");

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using DataFrames, DataFramesMeta, CSV
using Plots

pk_tang =  CSV.read("../data/Tang2024.csv",DataFrame);


include("mPBPK_deconj_homo.jl");
include("mPBPK_deconj_cyno.jl");
include("init.jl");
include("param.jl");

TotalCells = 1E-3/Vc; # assume tumor volume = 1mL

function single_sims(dose_mgkg, simstime = 21*DayToHour, p_tmp = p_datodxd)
    tspan = (0., simstime);    # [hr]
    # initial bolus IV dose
    u0_mpbpk, _ = pbpk_init(7, TotalCells, 0., 0, p_tmp.Rcopies, p_tmp.init_sR, p_tmp.k_endo, p_tmp.k_rec); 
    u0_mpbpk.C_EXG_Plasma = dose_mgkg*BW*1E3/(V_Plasma)/MW_EDG;
    # simulation for ADC 
    sol_mpbpk = solve(ODEProblem(mPBPK_deconj_homo!, u0_mpbpk, tspan, p_tmp), saveat = 3, alg = QNDF(autodiff=false), reltol = 1e-12);
    c_plasma_adc = [sol_mpbpk.u[i].C_EXG_Plasma for i in 1:length(sol_mpbpk.t)]; # [uM]
    # return the final plasma conc 
    return sol_mpbpk.t/DayToHour, c_plasma_adc*MW_EDG*1E-3 # [ug/mL]
end

# update params 
p_datodxd = ComponentArray(p_base, k_dec = 0.); # assume to have no deconjugation
p_datodxd.DAR = 4; # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10981107/
p_datodxd.Rcopies = 3.0E5  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4673178/
p_datodxd.init_sR = 0.  # assumed 
p_datodxd.PS_Score = 8.4

p_datodxd.Kd = 0.74E-3         # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9398094/
p_datodxd.ic50_pl = 9.54E-3
# updated based on dato-dxd internalization data
p_datodxd.k_deg = 6.22
p_datodxd.k_rec = 0.594
p_datodxd.k_endo = 0.575

t, adc = single_sims(6., 21*DayToHour, p_datodxd);

p_pk_homo = plot(title = "PK, Dato-Dxd (homo)", titlefontsize = 8, xlabel = "Time (Day)", ylabel = "Conc (ug/mL)", yaxis = :log);
plot!(t, adc, label = "sims", color = :red);
scatter!(pk_tang.time_day, pk_tang.conc_ugmL, label = "obs", color = :red, alpha = 0.7);

display(p_pk_homo)

savefig(p_pk_homo, "../figure/pk/dato-dxd.png");

#----------# cyno comparison #----------#
#=
pk_okajima =  CSV.read("../../data/okajima2021-cyno-pk.csv",DataFrame);

# simulation for cyno 
u0_mpbpk, _ = pbpk_init(7, TotalCells, 0., 0, p_datodxd.Rcopies, p_datodxd.init_sR, p_datodxd.k_endo, p_datodxd.k_rec); 
u0_mpbpk.C_EXG_Plasma = 6*6.2*1E3/(0.187)/MW_EDG;
sol_mpbpk = solve(ODEProblem(mPBPK_deconj_cyno!, u0_mpbpk, (0., 21*DayToHour), p_datodxd), saveat = 3, alg = QNDF(autodiff=false), reltol = 1e-12);
c_plasma_adc_cyno = [sol_mpbpk.u[i].C_EXG_Plasma for i in 1:length(sol_mpbpk.t)]*MW_EDG*1E-3; # [ug/mL]

p_pk_cyno = plot(title = "PK, Dato-Dxd (cyno)", titlefontsize = 8, xlabel = "Time (Day)", ylabel = "Conc (ug/mL)", yaxis = :log);
plot!(sol_mpbpk.t/DayToHour, c_plasma_adc_cyno, label = "sims", color = :red);
scatter!(pk_okajima.time_day, pk_okajima.conc_ugmL, label = "obs", color = :red, alpha = 0.7);
display(p_pk_cyno);

savefig(p_pk_cyno, "../figure/pk/dato-dxd-cyno.png");
=#
