# date: 10/23/24
# author: Yuezhe Li 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)

using DifferentialEquations, ComponentArrays
using Plots, DataFrames
using DataFramesMeta

# T-DM1 induced cell killing from https://aacrjournals.org/cancerres/article/68/22/9280/542757/Targeting-HER2-Positive-Breast-Cancer-with
tdm1_skbr3 = DataFrame(
    adcconc_ugmL = [0.00149, 0.0047, 0.013, 0.04, 0.114, 0.358, 1.09, 3.25], 
    rel_viability = [99.52, 89., 59.36, 34.64, 32.15, 32.53, 33.1, 32.69]/100
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

# SK-BR3 cell line 
p_skbr3_tdm1 = deepcopy(p_tdm1);
p_skbr3_tdm1.Rcopies = 1E4;
p_skbr3_tdm1.k_kill_max = 0.03;

Nc_skbr3_0 = cell_survive(0., 1E4, 48, p_skbr3_tdm1);
Nc_skbr3 = [];
for adcconc in tdm1_skbr3.adcconc_ugmL
    append!(Nc_skbr3, cell_survive(adcconc, 1E4, 48, p_skbr3_tdm1));
end
Nc_skbr3 = Nc_skbr3 / Nc_skbr3_0 

plot_tdm1_skbr3 = plot(xlabel = "ADC conc (ug/mL)", ylabel = "Relative viability (%)", dpi = 300);
plot!(tdm1_skbr3.adcconc_ugmL, Nc_skbr3*100, label = "sims"); 
scatter!(tdm1_skbr3.adcconc_ugmL, tdm1_skbr3.rel_viability*100, label = "Erickson et al., 2012", alpha = 0.8); 
plot!(xaxis = :log10);
display(plot_tdm1_skbr3);

savefig(plot_tdm1_skbr3, "deliv/figure/invitro/tdm1-sk-br-3.png");

# BT-474EEI cell line 
p_bt474eei_tdm1 = deepcopy(p_tdm1);
p_bt474eei_tdm1.tdouble_pos = 80; # https://www.cellosaurus.org/CVCL_0179
p_bt474eei_tdm1.Rcopies = 0.25E6; # https://aacrjournals.org/mct/article/11/5/1133/91314/The-Effect-of-Different-Linkers-on-Target-Cell
p_bt474eei_tdm1.ic50_pl = 0.1 # tuned

## T-DM1 killing curve; https://aacrjournals.org/mct/article/11/5/1133/91314/The-Effect-of-Different-Linkers-on-Target-Cell
tdm1_bt474eei = DataFrame(
    adcconc_ugmL = [0.000005, 0.00002, 0.00004, 0.00014, 0.00043, 0.001, 0.004, 0.011, 0.033, 0.1], 
    mean_colony_number = [25.9, 27.35, 31, 26, 29, 26.83, 24.37, 2.85, 0.13, 0.778] 
);
@transform!(tdm1_bt474eei, :rel_viability = :mean_colony_number ./ first(:mean_colony_number) );

Nc_bt474eei_0 = cell_survive(0., 1E4, 6*24, p_bt474eei_tdm1);
Nc_bt474eei = [];
for adcconc in tdm1_bt474eei.adcconc_ugmL
    append!(Nc_bt474eei, cell_survive(adcconc, 1E4, 6*24, p_bt474eei_tdm1));
end
Nc_bt474eei = Nc_bt474eei / Nc_bt474eei_0

plot_tdm1_bt474eei = plot(xlabel = "ADC conc (ug/mL)", ylabel = "Relative viability (%)", dpi = 300);
plot!(tdm1_bt474eei.adcconc_ugmL, Nc_bt474eei*100, label = "sims"); 
scatter!(tdm1_bt474eei.adcconc_ugmL, tdm1_bt474eei.rel_viability*100, label = "Erickson et al., 2012", alpha = 0.8); 
plot!(xaxis = :log10);
display(plot_tdm1_bt474eei);

savefig(plot_tdm1_bt474eei, "deliv/figure/invitro/tdm1-bt-474-eei.png");
