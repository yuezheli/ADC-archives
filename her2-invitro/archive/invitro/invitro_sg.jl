# date: May 16, 2024
# author: Yuezhe Li 
# purpose of this script: to fit in vitro data for sacituzumab govitecan
# observed data obtained from 

using Pkg; Pkg.activate("..");

using DifferentialEquations, ComponentArrays
using Plots, DataFrames, CSV, DataFramesMeta
using Optimization, OptimizationNLopt, OptimizationBBO

include("invitro_HER2.jl")

# update params for SG and TROP2 (NCI-N87 cell line) 
p_trop2sn38 = ComponentArray(
    tdouble = 30.,                          # doubling time, [hr], https://www.cellosaurus.org/CVCL_1603
    tdouble_neg = 30.,                      # doubling time, [hr], 
    tau = 24.,                              # time delay in tumor killing, [hr], https://www.oncotarget.com/article/4318/text/
    DAR = 6.8,                              # average drug : antibody ratio, [unitless], https://www.oncotarget.com/article/4318/text/
    ic50_pl = 2E-3,                         # payload IC50, [uM], https://www.oncotarget.com/article/4318/text/
    Emax_Payload = 0.201,                   # maximum in vitro killing rate, [1/hr], TO BE OPTIMIZED
    k_PL = 0.02,                            # free payload degredation rate, [1/hr], https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5564636/
    Rcopies = 2.0E5,                        # surface receptor copy number, [molecules], https://aacrjournals.org/mct/article/20/12/2329/675152
    k_deg = 6.22,                           # receptor degredation rate, [h-1], OPTIMIZED USING DATO-DXD DATA
    k_rec = 0.594,                          # receptor recycling rate, [h-1], OPTIMIZED USING DATO-DXD DATA
    k_endo = 0.575,                         # receptor:ADC endocytosis rate, [h-1], OPTIMIZED USING DATO-DXD DATA
    Kon = 2.07e3,                           # on-rate for receptor-ADC binding, [uM-1.hr-1], https://www.bio-rad.com/webroot/web/pdf/lsr/literature/Bulletin_6044A.pdf
    Kd = 0.3e-3,                            # ADC-receptor binding affinity, [uM], https://pubs.acs.org/doi/10.1021/acs.bioconjchem.5b00223
    k_out = 32.32,                          # rate of payload leaving of the cell through diffusion, [1/hr], TO BE OPTIMIZED
    k_in = 46.08,                           # rate of payload entering cell through diffusion, [1/hr], TO BE OPTIMIZED
    k_lys = 2.4,                            # rate for ADC to be graded in cells to release payload, [1/hr], computed from endosome/ lysosome half life https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2112524/pdf/jc972508.pdf
    V_medium = 1E-4,                        # medium volume, https://www.nature.com/articles/s41467-023-44533-z
    Nmax = 4E6,                             # maximum cell per well, 
);

# SN38 in vitro cytotoxicity; obtained from Cheng et al., 2022, Figure 7; https://pubmed.ncbi.nlm.nih.gov/36620535/
obs_sn38 = DataFrame(sn38_nM = [1E-3, 1E-2, 0.1, 1, 10, 100., 1E3, 1E4, 1E5], 
                      rel_viability = [83.36,81.75,82.16,49.83,26.36,7.67,0,0,0]); 

p_sn38 = ComponentArray(
    tdouble = 36.6,                         # doubling time, [hr], https://www.cellosaurus.org/CVCL_1258
    tau = 24.,                              # time delay in tumor killing, [hr], https://www.oncotarget.com/article/4318/text/
    ic50_pl = 2E-3,                    # payload IC50, [uM], https://www.oncotarget.com/article/4318/text/
    Emax_Payload = 0.201,                   # maximum in vitro killing rate, [1/hr], TO BE OPTIMIZED
    k_PL = 0.,                              # free payload degredation rate, [1/hr], TO BE OPTIMIZED
    k_in = 46.08,                           # rate of payload entering cell through diffusion, [1/hr], TO BE OPTIMIZED
    V_medium = 1E-4,                        # medium volume, https://www.nature.com/articles/s41467-023-44533-z
    Nmax = 4E6,                             # maximum cell per well, 
); 

u0_sn38 = ComponentArray(Nc_1 = 3E3, Nc_2 = 0, P_c = 0, P_m = 0.);

function invitro_pl!(du,u,p,t) 
    @unpack Nc_1, Nc_2, P_c, P_m = u
    @unpack tdouble, tau, Emax_Payload, ic50_pl, k_PL, k_in, V_medium, Nmax = p
    C_P_c = P_c/Vc                  
    C_P_m = P_m/V_medium
    Kgrow = log(2)/tdouble * (1 - (Nc_1+Nc_2)/Nmax)
    Kkill_eff = Emax_Payload*C_P_c/(ic50_pl + C_P_c)
    du.Nc_1 = Kgrow*Nc_1 - Kkill_eff*Nc_1
    du.Nc_2 = Kkill_eff*Nc_1 - Nc_2/tau
    du.P_c = k_in * (C_P_m - C_P_c) * Vc - Kkill_eff * C_P_c * Vc
    du.P_m = - k_in * (C_P_m - C_P_c) * Vc * Nc_1 - k_PL * P_m
end

function cell_survive(SN38conc, Nseed = 3000, exposuretime = 72., p_sn38 = p_sn38)
    u0 = deepcopy(u0_sn38); 
    u0.Nc_1 = Nseed;
    u0.P_m = p_sn38.V_medium * SN38conc * 1E-3  # [uM]
    sol_adcexposure = solve(ODEProblem(invitro_pl!, u0, (0., exposuretime), p_sn38), saveat = 1.);
    return sol_adcexposure[:Nc_1][end] 
end

# define loss function
function loss_viability(p, obs, opt = true, Nseed = 3000, exposuretime = 72., p_sn38 = p_sn38)
    opt_param = deepcopy(p_sn38);
    opt_param.Emax_Payload = p[1]
    opt_param.k_in = p[2]
    opt_param.tau = p[3]
    Nc_base = cell_survive(0., Nseed, exposuretime, opt_param)
    Nc_dosed = [];
    for sn38conc in obs.sn38_nM
        append!(Nc_dosed, cell_survive(sn38conc, Nseed, exposuretime, opt_param));
    end
    Nc_dosed = Nc_dosed / Nc_base * 100
    diff = sum((Nc_dosed .- obs.rel_viability).^2);
    if opt
        return diff; 
    else 
        return Nc_dosed;
    end
end 

f_tv = OptimizationFunction(loss_viability, Optimization.AutoForwardDiff());
sol_hcc3 = solve(OptimizationProblem(f_tv, [0.01, 5., 24.], obs_sn38, ub = [10., 1E4, 24.], lb = [0., 1., 1.]), NLopt.LN_NELDERMEAD()); 
# sol_hcc3 = [0.032, 11., 18.16]

pred_viability = loss_viability(sol_hcc3, obs_sn38, false); 

# plot simulation data 
plot_sn38 = plot(xlabel = "SN38 conc (nM)", ylabel = "Relative viability (%)");
plot!(obs_sn38.sn38_nM, pred_viability, label = "sims"); 
scatter!(obs_sn38.sn38_nM, obs_sn38.rel_viability, label = "Cheng et al., 2022", alpha = 0.8); 
plot!(xaxis = :log10);
display(plot_sn38); 

savefig(plot_sn38, "../figure/SN38.png");

# update parameters 
p_trop2sn38.Emax_Payload = sol_hcc3[1]
p_trop2sn38.k_in = sol_hcc3[2]
p_trop2sn38.tau = sol_hcc3[3]

# update SN38 uptake by cells based on literature
# in Itoh et al., 2005, 15min, Caco-2 cells uptake of 900 pmol/pg protein; # https://link.springer.com/article/10.1007/s00280-004-0937-4
# the medium SN38 concentration was 25uM 
# assuming 100ug-200ug per 1E6 cells # https://research.fredhutch.org/paulovich/en/typical-protein-yields.html
# then k_in between 27.6 hr-1 and 55.2 hr-1
p_trop2sn38.k_in = 27.6 

# save parameters 
final_param_trop2sn38 = DataFrame(
    names = [string.(keys(p_trop2sn38))[i] for i in 1:length(p_trop2sn38)], 
    values = [p_trop2sn38[i] for i in 1:length(p_trop2sn38)]
)

CSV.write("../tab/sg_params.csv", final_param_trop2sn38);
