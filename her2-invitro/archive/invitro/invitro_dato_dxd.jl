# date: Apr 23, 2024
# author: Yuezhe Li 
# purpose of this script: to fit in vitro data for dato-dxd
# observed data obtained from Okajima et al., 2021; https://aacrjournals.org/mct/article/20/12/2329/675152/Datopotamab-Deruxtecan-a-Novel-TROP2-directed

using Pkg; Pkg.activate("");

using DifferentialEquations, ComponentArrays
using Plots, DataFrames, CSV, DataFramesMeta
using Optimization, OptimizationNLopt, OptimizationBBO

include("invitro_HER2.jl")

# update params for dato-dxd and TROP2 (NCI-N87 cell line) 
p_trop2dxd = ComponentArray(
    tdouble = 30.,                          # doubling time, [hr], https://www.cellosaurus.org/CVCL_1603
    tdouble_neg = 30.,                      # doubling time, [hr], 
    tau = 1.,                               # time delay in tumor killing, [hr], 
    DAR = 4.,                               # average drug : antibody ratio, [unitless], https://aacrjournals.org/mct/article/20/12/2329/675152
    ic50_pl = 9.54E-3,                      # payload IC50, [uM], https://europepmc.org/article/ppr/ppr584999
    Emax_Payload = 0.201,                   # maximum in vitro killing rate, [1/hr], https://europepmc.org/article/ppr/ppr584999
    k_PL = 0.82,                            # free payload degredation rate, [1/hr], https://europepmc.org/article/ppr/ppr584999
    Rcopies = 2.0E5,                        # surface receptor copy number, [molecules], https://aacrjournals.org/mct/article/20/12/2329/675152
    k_deg = 0.4572,                         # receptor degredation rate, [h-1], TO BE OPTIMIZED
    k_rec = 0.01,                           # receptor recycling rate, [h-1], TO BE OPTIMIZED
    k_endo = 0.15,                          # receptor:ADC endocytosis rate, [h-1], TO BE OPTIMIZED
    Kon = 2.07e3,                           # on-rate for receptor-ADC binding, [uM-1.hr-1], https://www.bio-rad.com/webroot/web/pdf/lsr/literature/Bulletin_6044A.pdf
    Kd = 0.74e-3,                           # ADC-receptor binding affinity, [uM], https://aacrjournals.org/mct/article/20/12/2329/675152
    k_out = 32.32,                          # rate of payload leaving of the cell through diffusion, [1/hr], https://jpet.aspetjournals.org/content/374/1/184 and http://xlink.rsc.org/?DOI=C7ME00093F
    k_in = 46.08,                           # rate of payload entering cell through diffusion, [1/hr], https://jpet.aspetjournals.org/content/374/1/184 and http://xlink.rsc.org/?DOI=C7ME00093F
    k_lys = 2.4,                            # rate for ADC to be graded in cells to release payload, [1/hr], computed from endosome/ lysosome half life https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2112524/pdf/jc972508.pdf
    V_medium = 1E-4,                        # medium volume, https://www.nature.com/articles/s41467-023-44533-z
    Nmax = 4E6,                             # maximum cell per well, 
);

# update initial value (based on Okajima et al., 2021) 
u0_trop2dxd = ComponentArray(Nc_1 = 0, Nc_2 = 0, Nc_1_neg = 0, Nc_2_neg = 0, 
                    R_s = p_trop2dxd.Rcopies/N_av*1e6, R_e = p_trop2dxd.k_endo/p_trop2dxd.k_rec*p_trop2dxd.Rcopies/N_av*1e6, 
                    AR_s = 0, AR_e = 0, P_c = 0, P_m = 0, P_neg_c = 0, A_m = 0); 

#======================= internalization optimization =======================# 

function receptor_internalization(p_trop2dxd = p_trop2dxd, u0_trop2dxd = u0_trop2dxd, simstime = 3.)
    p_tmp = deepcopy(p_trop2dxd); 
    # remove toxicity due to the nature of the assay
    p_tmp.Emax_Payload = 0. 
    p_tmp.k_lys = 0.
    # update initial condition to reflect the binding 
    u0_tmp = deepcopy(u0_trop2dxd);
    u0_tmp.R_s = 0. 
    u0_tmp.AR_s = u0_trop2dxd.R_s
    # simulation 
    sol0 = solve(ODEProblem(invitro_model_her2!, u0_tmp, (0., simstime), p_tmp), saveat = [0., 0.25, 0.5, 1., 2., 3.], alg = QNDF(autodiff=false), reltol = 1e-18, abstol = 1E-18);
    AR_s = [sol0.u[i].AR_s for i in 1:length(sol0.t)];
    AR_e = [sol0.u[i].AR_e for i in 1:length(sol0.t)];
    return AR_e ./ (AR_s .+ AR_e) * 100
end

function internalization_loss(p, obs, opt = true, simstime = 3,  p_trop2dxd = p_trop2dxd)
    opt_param = deepcopy(p_trop2dxd);
    opt_param.k_deg = p[1]
    opt_param.k_endo = p[2]
    opt_param.k_rec = p[3]
    # update init condition 
    u0_opt = ComponentArray(Nc_1 = 0, Nc_2 = 0, Nc_3 = 0, Nc_4 = 0, Nc_1_neg = 0, Nc_2_neg = 0, Nc_3_neg = 0, Nc_4_neg = 0, 
                    R_s = opt_param.Rcopies/N_av*1e6, R_e = opt_param.k_endo/opt_param.k_rec*opt_param.Rcopies/N_av*1e6, 
                    AR_s = 0, AR_e = 0, P_c = 0, P_m = 0, P_neg_c = 0, A_m = 0); 
    # simulation 
    sims_internalized_percentage = receptor_internalization(opt_param, u0_opt, simstime)
    # calc diff 
    if opt
        return sum( (sims_internalized_percentage .- obs.percentage).^2 )
    else 
        return sims_internalized_percentage
    end
end

nci_n87_internalization =  CSV.read("../../data/okajima2021-internalization.csv",DataFrame);

f_int = OptimizationFunction(internalization_loss, Optimization.AutoForwardDiff());

p_int_ncin87 = solve(OptimizationProblem(f_int, [0.45, 0.2, 0.01], nci_n87_internalization, ub = [10., 10., 1.], lb = [0., 0., 0.]), NLopt.LN_NELDERMEAD()); # [6.22, 0.575, 0.594]

p_pred_int_nvin87 = internalization_loss(p_int_ncin87, nci_n87_internalization, false);

plot_int = plot(title = "Dato-Dxd internalization", titlefontsize = 8, xlabel = "Time (min)", ylabel = "Internalization rate (%)");
plot!(nci_n87_internalization.time_min, p_pred_int_nvin87, label = "pred, NCI-N87");
plot!(nci_n87_internalization.time_min, nci_n87_internalization.percentage, label = "Okajima et al., 2021", seriestype = :scatter);

display(plot_int);

savefig(plot_int, "figure/invitro/dato-dxd-internalization.png");

#======================= cell-killing optimization =======================# 
p_trop2dxd.k_deg = p_int_ncin87[1]
p_trop2dxd.k_endo = p_int_ncin87[2]
p_trop2dxd.k_rec = p_int_ncin87[3]
p_trop2dxd.k_PL = 0. 

function Dxd_release(receptorcopy, p_trop2dxd = p_trop2dxd, seedingtime = 12., incubationtime = 24., dxd_mw = 494)
    p_tmp = deepcopy(p_trop2dxd);
    p_tmp.Rcopies = receptorcopy;
    p_tmp.V_medium = 2E-3;  # update to account for 6-well plate
    u0 = ComponentArray(Nc_1 = 3E5, Nc_2 = 0, Nc_3 = 0, Nc_4 = 0, Nc_1_neg = 0, Nc_2_neg = 0, Nc_3_neg = 0, Nc_4_neg = 0, 
                    R_s = p_tmp.Rcopies/N_av*1e6, R_e = p_tmp.k_endo/p_tmp.k_rec*p_tmp.Rcopies/N_av*1e6, 
                    AR_s = 0, AR_e = 0, P_c = 0, P_m = 0, P_neg_c = 0, A_m = 0);
    # overnight seeding 
    sol0 = solve(ODEProblem(invitro_model_her2!, u0, (0., seedingtime), p_tmp), saveat = [0., seedingtime], alg = QNDF(autodiff=false));
    # dxd release 
    u02 = deepcopy(sol0.u[end]);
    u02.A_m = p_tmp.V_medium * 0.1; # [umol]
    sol2 = solve(ODEProblem(invitro_model_her2!, u02, (0., incubationtime), p_tmp), saveat = [0., incubationtime], alg = QNDF(autodiff=false), reltol = 1e-18);
    return sol2.u[end].P_m/p_tmp.V_medium * dxd_mw; # [ug/L]
end

function loss_dxd_release(p, obs, opt = true, p_trop2dxd = p_trop2dxd)
    opt_param = deepcopy(p_trop2dxd);
    opt_param.Emax_Payload = p[1]
    # simulation
    sims_dxd = []; 
    for trop2 in obs.TROP2
        append!( sims_dxd, Dxd_release(trop2, opt_param, 12.) );
    end
    if opt
        return sum( (sims_dxd .- obs.dxd_conc_ngmL).^2 )
    else
        return sims_dxd
    end
end

dxd_release =  CSV.read("../../data/okajima2021_dxd_release.csv",DataFrame);

f_dxd = OptimizationFunction(loss_dxd_release, Optimization.AutoForwardDiff());

p_dxd_rel = solve(OptimizationProblem(f_dxd, [0.2], dxd_release, ub = [1.], lb = [0.]), NLopt.LN_NELDERMEAD()); # [1.]

p_pred_dxd_rel = loss_dxd_release(p_dxd_rel, dxd_release, false); 

plot_dxd_rel = plot(title = "Dxd release", titlefontsize = 8, xlabel = "TROP2 copy number (#)", ylabel = "Dxd conc (ng/mL)");
plot!(dxd_release.TROP2, p_pred_dxd_rel, label = "pred, NCI-N87");
plot!(dxd_release.TROP2, dxd_release.dxd_conc_ngmL, label = "Okajima et al., 2021", seriestype = :scatter);

display(plot_dxd_rel);

savefig(plot_dxd_rel, "../figure/dato-dxd-dxd-release.png");

# sensitivity analysis over the Dxd release assay and Emax_Payload
emax_test = [0.1, 0.2, 0.5, 1., 2.];
pred_dxd_rel = [];
for emax in emax_test
    append!(pred_dxd_rel, [loss_dxd_release(p_dxd_rel, dxd_release, false)])
end

plot_dxd_rel_sens = plot(title = "Dxd release", titlefontsize = 8, xlabel = "TROP2 copy number (#)", ylabel = "Dxd conc (ng/mL)");
plot!(dxd_release.TROP2, dxd_release.dxd_conc_ngmL, label = "Okajima et al., 2021", seriestype = :scatter);
plot!(dxd_release.TROP2, pred_dxd_rel[1], label = "pred, Emax = 0.1");
plot!(dxd_release.TROP2, pred_dxd_rel[2], label = "pred, Emax = 0.2");
plot!(dxd_release.TROP2, pred_dxd_rel[3], label = "pred, Emax = 0.5");
plot!(dxd_release.TROP2, pred_dxd_rel[4], label = "pred, Emax = 1.0");
plot!(dxd_release.TROP2, pred_dxd_rel[5], label = "pred, Emax = 2.0");

display(plot_dxd_rel_sens);

savefig(plot_dxd_rel_sens, "../figure/dato-dxd-dxd-release-sensivitity-analysis.png");

# save parameters 
final_param_trop2dxd = DataFrame(
    names = [string.(keys(p_trop2dxd))[i] for i in 1:length(p_trop2dxd)], 
    values = [p_trop2dxd[i] for i in 1:length(p_trop2dxd)]
)

CSV.write("../tab/dadt_dxd_params.csv", final_param_trop2dxd);
