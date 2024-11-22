# date: 11/21/2024
# author: Yuezhe Li 
# purpose of this code: to optimize parameters for ARX788

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)

using DifferentialEquations, ComponentArrays
using Plots, DataFrames
using DataFramesMeta
using Optimization, OptimizationNLopt
using Statistics

# data on ARX788 treated and untreated data; in vitro
# obtained from Shastri et al., 2020; # https://pubmed.ncbi.nlm.nih.gov/32499302/
untreated = DataFrame(
    time_hr = [0., 12., 16., 24.], 
    cell_percentage = [100., 100., 125., 156.]
);

arx788treated = DataFrame(
    time_hr = [0., 12., 16., 24.], 
    cell_percentage = [100., 98., 104., 100.], 
    medium_pl = [0., 0.24, 0.47, 1.01] * 1E-3 / 731.976     # ng -> umol 
);

include(@projectroot("model/invitro.jl"));
include(@projectroot("model/param-tdm1.jl"));

p_tmmaf = deepcopy(p_tdm1);
p_tmmaf.tdouble_pos = 43        # doubling time, [hr], https://www.cellosaurus.org/CVCL_1259
p_tmmaf.tdouble_neg = 40.6      # doubling time, [hr], https://www.cellosaurus.org/CVCL_0419
p_tmmaf.DAR = 2                 # average drug : antibody ratio, [unitless], https://pubmed.ncbi.nlm.nih.gov/32499302/
p_tmmaf.ic50_pl = 0.2E-3        # payload IC50, [uM], https://pubmed.ncbi.nlm.nih.gov/16417259/
p_tmmaf.Rcopies = 1.88E6        # surface receptor copy number, [molecules], calculated based on https://pubmed.ncbi.nlm.nih.gov/26766593/ and https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4622696/
p_tmmaf.V_medium = 2E-3         # 6-well plate, # https://pubmed.ncbi.nlm.nih.gov/32499302/
p_tmmaf.tau = 12                # assumed

# init value
u0_tmmaf = ComponentArray(Nc_1 = 0.3E6, Nc_2 = 0, Nc_1_neg = 0, Nc_2_neg = 0, 
             R_s = p_tmmaf.Rcopies/N_av*1e6, 
             R_e = p_tmmaf.k_endo/p_tmmaf.k_rec*p_tmmaf.Rcopies/N_av*1e6, 
             AR_s = 0, AR_e = 0, P_c_pos = 0, P_m = 0, P_c_neg = 0, A_m = 0, A_e = 0); 

# simulation without drug 
p_tmmaf.tdouble_pos = 34.  # doubling time fitted based on Shastri 2020
sol0 = solve(ODEProblem(invitro_model_her2!, u0_tmmaf, (0., 24.), p_tmmaf), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-18);
# plot(sol0.t, [sol0.u[i].Nc_1 for i in 1:length(sol0.t)]/u0_tmmaf.Nc_1*100, xlabel = "Time (hr)", ylabel = "Normalized Cell Number", label = false, ylim = [60, 180])

# simulation with ARX788 treated
function loss_fit(p, obs, ADCconc_ugmL = 5, opt = true,  p_her2mmaf = p_tmmaf, u0_her2mmaf = u0_tmmaf)
    opt_param = deepcopy(p_her2mmaf);
    opt_param.k_kill_max = exp(p[1])
    opt_param.k_out = exp(p[2])
    opt_param.k_PL = exp(p[3])

    u0_tmp = deepcopy(u0_her2mmaf); 
    u0_tmp.A_m = ADCconc_ugmL / MW * 1E3 * p_her2mmaf.V_medium; # ug/mL -> umol
    sol_opt = solve(ODEProblem(invitro_model_her2!, u0_tmp, (0., 24.), opt_param), saveat = obs.time_hr, alg = QNDF(autodiff=false), reltol = 1e-18);
    diff1 = sum( (([sol_opt.u[i].Nc_1 for i in 1:length(sol_opt.t)]/u0_tmp.Nc_1*100 .- obs.cell_percentage)./obs.cell_percentage).^2 )
    diff3 = sum( (([sol_opt.u[i].P_m for i in 1:length(sol_opt.t)] .- obs.medium_pl)./mean(obs.medium_pl)).^2 )
    diff = diff1/obs.cell_percentage[end] + diff3/obs.medium_pl[end]
    if opt
        return diff; 
    else 
        return [sol_opt.u[i].Nc_1 for i in 1:length(sol_opt.t)]/u0_tmp.Nc_1*100, [sol_opt.u[i].P_m for i in 1:length(sol_opt.t)]
    end
end

# visualization 
function visualization(log_parameters, ADCconc_ugmL = 5.; obs = arx788treated,  p_her2mmaf = p_tmmaf, u0_her2mmaf = u0_tmmaf)

    cellperc, p_m = loss_fit(log_parameters, obs, ADCconc_ugmL, false, p_her2mmaf, u0_her2mmaf); 

    plt1 = plot(xlabel = "Time (hr)", ylabel = "Normalized Cell Number (%)", ylim = [60, 180], legend = :bottomright, dpi = 300); 
    plot!(sol0.t, [sol0.u[i].Nc_1 for i in 1:length(sol0.t)]/u0_tmmaf.Nc_1*100, label = "sims, untreated", color = :blue); 
    plot!(untreated.time_hr, untreated.cell_percentage, seriestype=:scatter, label = "obs, untreated", color = :blue); 
    plot!(arx788treated.time_hr, cellperc, label = "sims, ARX788 treated", color = :red); 
    plot!(arx788treated.time_hr, arx788treated.cell_percentage, seriestype=:scatter, label = "obs, ARX788 treated", color = :red);  

    plt2 = plot(xlabel = "Time (hr)", ylabel = "Medium payload (umol)", legend = :topleft, dpi = 300); 
    plot!(arx788treated.time_hr, p_m, label = "sims, ARX788 treated"); 
    plot!(arx788treated.time_hr, arx788treated.medium_pl, seriestype=:scatter, label = "obs, ARX788 treated"); 

    plotout = plot(plt1, plt2)
    return plotout

end

p_preoptimization  = visualization(log.([0.02, 0.01, 1.]))

#f_tv = OptimizationFunction(loss_fit, Optimization.AutoForwardDiff());
#sol_hcc1954 = solve(OptimizationProblem(f_tv, log.([0.02, 0.01, 1.]), arx788treated  ), NLopt.LN_NELDERMEAD());
#p_postoptimization = visualization(exp.(sol_hcc1954))

savefig(p_preoptimization, @projectroot("deliv/figure/invitro/ARX788_HCC1954.png") );

