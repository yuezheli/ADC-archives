# date: 11/21/2024
# author: Yuezhe Li 
# purpose of this code: to optimize parameters for SYD985, trastuzumab + cleavable linker + duocarmycin; 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)

using DifferentialEquations, ComponentArrays
using Plots, DataFrames
using DataFramesMeta
using Optimization, OptimizationNLopt
using Statistics

# cytotoxicity assay; Fig 2, van der Lee et al., Mol Cancer Ther, 2015; https://pubmed.ncbi.nlm.nih.gov/25589493/
obs_skbr3 = DataFrame(adcconc_ngmL = [0.3, 1.0, 3.32, 10., 31.34, 98.8, 3106.8], 
                      rel_viability = [100., 100., 87.2, 30., 3., 1.6, 1.1]); 

include(@projectroot("model/invitro.jl"));
include(@projectroot("model/param-tdm1.jl"));
                      
p_syd985 = deepcopy(p_tdm1)
p_syd985.DAR = 2.8          # average drug : antibody ratio, [unitless], https://pubmed.ncbi.nlm.nih.gov/25589493/
p_syd985.ic50_pl = 0.08E-3   # payload IC50, [uM], https://pubmed.ncbi.nlm.nih.gov/25589493/

#======================= cytotoxicity =======================# 
function cell_survive(ADCconc, Nseed = 6500, seedtime = 12., adcexposuretime = 48., regrowthtime = 120.,  p_base = p_syd985)
    u0_base = ComponentArray(Nc_1 = Nseed, Nc_2 = 0, Nc_1_neg = 0, Nc_2_neg = 0, 
             R_s = p_base.Rcopies/N_av*1e6*Nseed, 
             R_e = p_base.k_endo/p_base.k_rec*p_base.Rcopies/N_av*1e6*Nseed, 
             AR_s = 0, AR_e = 0, P_c_pos = 0, P_m = 0, P_c_neg = 0, 
             A_m = 0., A_e = 0.); 
    sol_seed = solve(ODEProblem(invitro_model_her2!, u0_base, (0., seedtime), p_base), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-18);
    u0 = deepcopy(sol_seed.u[end]); 
    u0.A_m = ADCconc / MW; # convert ADCconc from ng/mL to uM
    sol_adcexposure = solve(ODEProblem(invitro_model_her2!, u0, (0., adcexposuretime), p_base), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-18);
    u0 = deepcopy(sol_adcexposure.u[end]); 
    u0.A_m = 0. 
    u0.AR_s = 0. 
    u0.P_m = 0.
    if regrowthtime > 0.
        sol_regrowth = solve(ODEProblem(invitro_model_her2!, u0, (0., regrowthtime), p_base), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-18);
        return sol_regrowth.u[end].Nc_1
    else 
        return sol_adcexposure.u[end].Nc_1
    end
end

function cell_viability_visual(p_base = p_syd985)
    Nc_base_skbr3 = cell_survive(0., 6500., 18., 6., 120, p_base);
    Nc_skbr3 = [];
    for adcconc in obs_skbr3.adcconc_ngmL
        append!(Nc_skbr3, cell_survive(adcconc, 6500., 18., 6., 120, p_base));
    end
    Nc_skbr3 = Nc_skbr3 / Nc_base_skbr3 * 100

    plot_syd_skbr3 = plot(xlabel = "ADC conc (ug/L)", ylabel = "Relative viability (%)", dpi = 300);
    plot!(obs_skbr3.adcconc_ngmL, Nc_skbr3, label = "sims"); 
    scatter!(obs_skbr3.adcconc_ngmL, obs_skbr3.rel_viability, label = "van der Lee et al., 2015", alpha = 0.8); 
    plot!(xaxis = :log10);
    
    return plot_syd_skbr3
end

# hand-tuning
p_syd985.k_kill_max = 0.04
p_syd985.k_out = 5.
p_syd985.k_PL = 0.01
p_preoptimzation = cell_viability_visual(p_syd985)

savefig(p_preoptimzation, @projectroot("deliv/figure/invitro/syd985-sk-br-3.png"));

#======================= bystander effect =======================# 
# optimization k_in based on co-culture experiment 
function bystander_cell_survival(ADCconc, Nseed = 1E4, seedtime = 4., adcexposuretime = 48., regrowthtime = 120.,  p_base = p_syd985)
    u0 = ComponentArray(Nc_1 = Nseed/2, Nc_2 = 0, Nc_1_neg = Nseed/2, Nc_2_neg = 0, 
             R_s = p_base.Rcopies/N_av*1e6*Nseed/2, 
             R_e = p_base.k_endo/p_base.k_rec*p_base.Rcopies/N_av*1e6*Nseed/2, 
             AR_s = 0, AR_e = 0, P_c_pos = 0, P_m = 0, P_c_neg = 0, 
             A_m = 0., A_e = 0.); 
    sol_seed = solve(ODEProblem(invitro_model_her2!, u0, (0., seedtime), p_base), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-18);
    u0 = deepcopy(sol_seed.u[end]); 
    u0.A_m = ADCconc / MW; # convert ADCconc from ng/mL to uM
    sol_adcexposure = solve(ODEProblem(invitro_model_her2!, u0, (0., adcexposuretime), p_base), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-18);
    u0 = deepcopy(sol_adcexposure.u[end]); 
    u0.A_m = 0. 
    u0.AR_s = 0. 
    u0.P_m = 0.
    if regrowthtime > 0.
        sol_regrowth = solve(ODEProblem(invitro_model_her2!, u0, (0., regrowthtime), p_base), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-18);
        return (sol_regrowth.u[end].Nc_1_neg + sol_regrowth.u[end].Nc_2_neg)/(sol_regrowth.u[end].Nc_1_neg + sol_regrowth.u[end].Nc_2_neg + sol_regrowth.u[end].Nc_1 + sol_regrowth.u[end].Nc_2)
    else 
        return (sol_adcexposure.u[end].Nc_1_neg + sol_adcexposure.u[end].Nc_2_neg)/(sol_adcexposure.u[end].Nc_1_neg + sol_adcexposure.u[end].Nc_2_neg + sol_adcexposure.u[end].Nc_1 + sol_adcexposure.u[end].Nc_2 ) 
    end
end

# data obtained from Figure 3, van der Lee et al., 2015; # https://aacrjournals.org/mct/article/14/3/692/137131/The-Preclinical-Profile-of-the-Duocarmycin-Based
# bystander cell survival should be ~ 0.2
p_syd985.k_in = 2.
bystander_cell_survival(0.1E3/MW, 1E4, 4, 144., 0., p_syd985)

# save results
final_param_syd985 = DataFrame(
    names = [string.(keys(p_syd985))[i] for i in 1:length(p_syd985)], 
    values = [p_syd985[i] for i in 1:length(p_syd985)]
)
