# date: May 21, 2024 
# author: Yuezhe Li 
# purpose of this code: to link ADC IC50 with HER2 expression (using T-DM1 as an example)

using Pkg; Pkg.activate("..");

using DifferentialEquations, ComponentArrays
using Plots, DataFrames
using Optimization, OptimizationNLopt

include("invitro_HER2.jl")

p_tdm1 = ComponentArray(
    tdouble = 37.4,                         # doubling time, [hr], https://www.cellosaurus.org/CVCL_0033
    tdouble_neg = 61.,                      # doubling time, [hr], https://www.cellosaurus.org/CVCL_1566
    tau = 1.,                               # time delay in tumor killing, [hr], 
    DAR = 3.5,                              
    ic50_pl = 23.8E-3,                      # payload IC50, [uM], 
    Emax_Payload = 0.0139,                  # maximum in vitro killing rate, [1/hr]
    k_PL = 0.34,                            # free payload degredation rate, [1/hr], 
    Rcopies = 1.6E6,                        # surface receptor copy number, [molecules], https://pubmed.ncbi.nlm.nih.gov/26766593/
    k_deg = 0.4572,                         # HER2 degredation rate, [h-1], https://europepmc.org/article/ppr/ppr584999)
    k_rec = 0.0864,                         # HER2 recycling rate, [h-1], https://europepmc.org/article/ppr/ppr584999)
    k_endo = 0.15,                          # HER2:ADC endocytosis rate, [h-1], https://europepmc.org/article/ppr/ppr584999)
    Kon = 2.07e3,                           # on-rate for receptor-ADC binding, [uM-1.hr-1], https://www.bio-rad.com/webroot/web/pdf/lsr/literature/Bulletin_6044A.pdf
    Kd = 0.314e-3,                          # ADC-receptor binding affinity, [uM], https://www.nature.com/articles/s41467-023-44533-z
    k_out = 0.14,                           # rate of payload leaving of the cell through diffusion, [1/hr], 
    k_in = 0.21,                            # rate of payload entering cell through diffusion, [1/hr], 
    k_lys = 2.4,                            # rate for ADC to be graded in cells to release payload, [1/hr], computed from endosome/ lysosome half life https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2112524/pdf/jc972508.pdf
    V_medium = 1E-4,                        # medium volume, https://www.nature.com/articles/s41467-023-44533-z
    Nmax = 4E6,                             # maximum cell per well, 
);

function cell_survive(ADCconc, Nseed = 6500, adcexposuretime = 48., p_tmp = p_tdm1)
    # update initial value 
    u0_tdm1 = ComponentArray(Nc_1 = Nseed, Nc_2 = 0, Nc_1_neg = 0, Nc_2_neg = 0, AR_s = 0, AR_e = 0, P_c = 0, P_m = 0, P_neg_c = 0, 
        R_s = p_tmp.Rcopies/N_av*1e6, R_e = p_tmp.k_endo/p_tmp.k_rec*p_tmp.Rcopies/N_av*1e6, 
        A_m = p_tmp.V_medium * ADCconc); 
    sol_adcexposure = solve(ODEProblem(invitro_model_her2!, u0_tdm1, (0., adcexposuretime), p_tmp), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-18);
    return sol_adcexposure[:Nc_1][end] 
end

function normalized_cell_survival(ADCconc_list, Rcopies = 1E6, Nseed = 3E5, adcexposuretime = 48., p_tdm1 = p_tdm1)
    p_tmp = deepcopy(p_tdm1); 
    p_tmp.Rcopies = Rcopies; 
    # baseline, no drug 
    N0 = cell_survive(0, Nseed, adcexposuretime, p_tmp); 
    N_list = [];
    for adcconc in ADCconc_list
        append!(N_list, cell_survive(adcconc, Nseed, adcexposuretime, p_tmp) ); 
    end
    return N_list/N0; 
end

tmd1_conc = [1E-8, 1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1, 1E2];
nc = normalized_cell_survival(tmd1_conc, 1E6)

obs = DataFrame(ADCconc = tmd1_conc, normalized_cell_survival = nc); 

# analytical solution: IC50 = -k2 * k1
function sigmoid(x, k1, k2, k3, k4)
    return k3 ./ (1.0 .+ exp.(x./k1 .+ k2)) .+ k4
end

function ic50_loss(p, obs, runOpt = true)
    ADCconc = obs.ADCconc
    normalized_cell_survival = obs.normalized_cell_survival
    log_adc_conc = log10.(ADCconc)
    pred_surv = sigmoid(log_adc_conc, p[1], p[2], p[3], p[4])
    diff = sum( (pred_surv .- normalized_cell_survival).^2 )
    if runOpt
        return diff
    else
        return pred_surv
    end
end

function ic50_fit(obs, p_init = [1., 1., 0.3, 0.1])
    f_ic50 = OptimizationFunction(ic50_loss, Optimization.AutoForwardDiff());
    sol = solve(OptimizationProblem(f_ic50, p_init, obs, ub = [10., 10., 1., 1.], lb = [-10., -10., 0., 0.]), NLopt.LN_NELDERMEAD()); 
    ic50 = 10^( -sol[2] * sol[1] );
    pred = ic50_loss(sol, obs, false)
    return ic50, pred
end

ic50, pred = ic50_fit(obs)

plot(obs.ADCconc, obs.normalized_cell_survival, xaxis = :log10, seriestype=:scatter, label = false, xticks = tmd1_conc);
plot!(obs.ADCconc, pred, label = false)

# HER2 expression scanning 
function HER2_tdm1_ic50(HER2copies = [1E2, 5E2, 1E3, 5E3, 1E4, 1E5, 1E6], tmd1_conc = [1E-8, 1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1, 1E2])
    IC50 = [];
    data = DataFrame(ADCconc = [], normalized_cell_survival = [], type = [], HER2copies = [])
    for her2copies in HER2copies
        nc_tmp = normalized_cell_survival(tmd1_conc, her2copies)
        obs_tmp = DataFrame(ADCconc = tmd1_conc, normalized_cell_survival = nc_tmp, type = "sims"); 
        obs_tmp[!, :HER2copies] .= her2copies;
        ic50_tmp, pred_tmp = ic50_fit(obs_tmp); 
        fitted = DataFrame(ADCconc = tmd1_conc, normalized_cell_survival = pred_tmp, type = "fitted"); 
        fitted[!, :HER2copies] .= her2copies;
        append!(IC50, ic50_tmp);
        data = vcat(data, obs_tmp, fitted); 
    end
    return IC50, data
end

ic50, data = HER2_tdm1_ic50();

sims = filter(:type => type -> type == "sims", data);
fitted = filter(:type => type -> type == "fitted", data);

plt = plot(xaxis = :log10, xticks = tmd1_conc, palette = :Dark2_7, legend = :bottomleft, xlabel = "T-DM1 concentration (uM)", ylabel = "Cell survival"); 
scatter!(sims.ADCconc, sims.normalized_cell_survival, group = sims.HER2copies, ma = 0.6, linewidth=0); 
plot!(fitted.ADCconc, fitted.normalized_cell_survival, group = fitted.HER2copies); 
display(plt);

savefig(plt, "../figure/HER2_TDM1_IC50.png");

