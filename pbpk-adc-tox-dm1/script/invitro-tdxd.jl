# date: 11/20/2024
# author: Yuezhe Li 
# purpose of this code: to test parameters for in vitro cytotoxicity of T-Dxd

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)

using DifferentialEquations, ComponentArrays
using Plots, DataFrames

# T-Dxd induced cell killing from https://aacrjournals.org/clincancerres/article/22/20/5097/124857/DS-8201a-A-Novel-HER2-Targeting-ADC-with-a-Novel
tdxd_kpl4 = DataFrame(
    adcconc_ugmL = [0.67, 3.68, 16.07, 79.12, 389.8, 2032.31, 9573.16]/1E3, 
    rel_viability = [98.03, 100.75, 73.57, 15.55, 8.46, 7.91, 5.96]/100
);

include(@projectroot("model/invitro.jl"));
include(@projectroot("model/param-invitro-tdm1.jl"));

function cell_survive(ADCconc, Nseed = 1E4, simstime = 48,  p_base = p_tdm1)
    u0_base = ComponentArray(Nc_1 = Nseed, Nc_2 = 0, Nc_1_neg = 0, Nc_2_neg = 0, 
             R_s = p_base.Rcopies/N_av*1e6*Nseed, 
             R_e = p_base.k_endo/p_base.k_rec*p_base.Rcopies/N_av*1e6*Nseed, 
             AR_s = 0, AR_e = 0, P_c_pos = 0, P_m = 0, P_c_neg = 0, 
             A_m = ADCconc / MW * 1E3, A_e = 0.); 
    sol0 = solve(ODEProblem(invitro_model_her2!, u0_base, (0., simstime), p_base), saveat = 1., reltol = 1e-18);
    return sol0.u[end].Nc_1
end 

# update parameters for T-Dxd, based on Scheuher et al., 2023; https://link.springer.com/article/10.1007/s10928-023-09884-6
p_tdxd = deepcopy(p_tdm1);
p_tdxd.DAR = 8.
p_tdxd.k_in = 46.08
p_tdxd.k_out = 32.32
p_tdxd.k_kill_max = 0.201
p_tdxd.ic50_pl = 9.54E-3
# p_tdxd.k_PL_ex = 0.82

# below were tuned
p_tdxd.Rcopies = 1E4;

Nc_tdxd_0 = cell_survive(0., 1E4, 48, p_tdxd);
Nc_tdxd_ = [];
for adcconc in tdxd_kpl4.adcconc_ugmL
    append!(Nc_tdxd_, cell_survive(adcconc, 1E4, 48, p_tdxd));
end
Nc_tdxd_ = Nc_tdxd_ / Nc_tdxd_0 

plot_tdxd_ = plot(xlabel = "ADC conc (ug/mL)", ylabel = "Relative viability (%)", dpi = 300);
plot!(tdxd_kpl4.adcconc_ugmL, Nc_tdxd_*100, label = "sims"); 
scatter!(tdxd_kpl4.adcconc_ugmL, tdxd_kpl4.rel_viability*100, label = "Ogitani et al., 2016", alpha = 0.8); 
plot!(xaxis = :log10);
display(plot_tdxd_);

savefig(plot_tdxd_, "deliv/figure/invitro/tdxd-kpl4.png");
