# date: 11/21/2024
# author: Yuezhe Li 
# purpose of this code: to optimize in vitro parameters for DHES0815A

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)

using DifferentialEquations, ComponentArrays
using Plots, DataFrames
using DataFramesMeta

# observed data, obtained from Fig 2a, Lewis et al., 2024; https://www.nature.com/articles/s41467-023-44533-z
dhes0815a_skbr3 = DataFrame(
    adcconc_ugmL = [10, 3.33, 1.11, 0.37, 0.123, 0.041, 0.0137, 0.0046, 0.0015, 0.0005, 0.0001], 
    rel_viability = [12, 18, 23, 27, 28, 36, 69, 94, 100, 99, 100]
);

dhes0815a_kpl4 = DataFrame(
    adcconc_ugmL = [100, 33.33, 11.11, 3.7, 1.23, 0.41, 0.13, 0.0457, 0.0152, 0.005, 0.001], 
    rel_viability = [2, 3, 8, 31, 63, 83, 94, 97, 98, 97, 100]
);

include(@projectroot("model/invitro.jl"));
include(@projectroot("model/param-invitro-tdm1.jl"));

p_her2pbd = deepcopy(p_tdm1);
p_her2pbd.tdouble_pos = 37.4        # doubling time, [hr], https://www.cellosaurus.org/CVCL_0033 (SK-BR-3)
p_her2pbd.tdouble_neg = 24.         # doubling time, [hr], https://pubmed.ncbi.nlm.nih.gov/6871841/ (MCF-7)
p_her2pbd.DAR = 2                   # average drug : antibody ratio, [unitless], https://www.nature.com/articles/s41416-019-0635-y
p_her2pbd.ic50_pl = 27.1E-3         # payload IC50, [uM], https://www.nature.com/articles/s41467-023-44533-z
p_her2pbd.Rcopies = 1.1E6           # surface receptor copy number, [molecules], https://aacrjournals.org/mct/article/19/9/1833/92916/ARX788-a-Site-specific-Anti-HER2-Antibody-Drug
p_her2pbd.V_medium = 1E-4           # https://www.nature.com/articles/s41467-023-44533-z
p_her2pbd.Kd = 0.57E-3              # https://www.nature.com/articles/s41467-023-44533-z

#======================= cytotoxicity =======================# 
# cell survival simulation
function cell_survive(ADCconc, Nseed = 1E4, simstime = 48,  p_base = p_tdm1)
    u0_base = ComponentArray(Nc_1 = Nseed, Nc_2 = 0, Nc_1_neg = 0, Nc_2_neg = 0, 
             R_s = p_base.Rcopies/N_av*1e6*Nseed, 
             R_e = p_base.k_endo/p_base.k_rec*p_base.Rcopies/N_av*1e6*Nseed, 
             AR_s = 0, AR_e = 0, P_c_pos = 0, P_m = 0, P_c_neg = 0, 
             A_m = ADCconc / MW * 1E3, A_e = 0.); 
    sol0 = solve(ODEProblem(invitro_model_her2!, u0_base, (0., simstime), p_base), saveat = 1., reltol = 1e-18);
    return sol0.u[end].Nc_1
end 

p_her2pbd.k_kill_max = 0.02  
p_her2pbd.k_out = 2.   

Nc_base_skbr3 = cell_survive(0., 1E4, 72, p_her2pbd);
Nc_skbr3 = [];
for adcconc in dhes0815a_skbr3.adcconc_ugmL
    append!(Nc_skbr3, cell_survive(adcconc, 1E4, 72, p_her2pbd));
end
Nc_skbr3 = Nc_skbr3 / Nc_base_skbr3 * 100

plot_her2pbd_skbr3 = plot(xlabel = "ADC conc (ug/mL)", ylabel = "Relative viability (%)", dpi = 300);
plot!(dhes0815a_skbr3.adcconc_ugmL, Nc_skbr3, label = "sims"); 
scatter!(dhes0815a_skbr3.adcconc_ugmL, dhes0815a_skbr3.rel_viability, label = "Lewis et al., 2024", alpha = 0.8); 
plot!(xaxis = :log10);
plot_her2pbd_skbr3

savefig(plot_her2pbd_skbr3, @projectroot("deliv/figure/invitro/dhes0815a-sk-br-3.png"));

#======================= bystander effect =======================# 
# optimization k_in based on co-culture experiment 

function bystander(ADCconc, Nseed = 1E4, simstime = 120.,  p_base = p_her2pbd)
    u0_base = ComponentArray(Nc_1 = Nseed, Nc_2 = 0, Nc_1_neg = 0, Nc_2_neg = 0, 
             R_s = p_base.Rcopies/N_av*1e6*Nseed, 
             R_e = p_base.k_endo/p_base.k_rec*p_base.Rcopies/N_av*1e6*Nseed, 
             AR_s = 0, AR_e = 0, P_c_pos = 0, P_m = 0, P_c_neg = 0, 
             A_m = ADCconc / MW * 1E3, A_e = 0.); 
    u0_base.A_m = p_base.V_medium * ADCconc / MW * 1E3; # convert ADCconc from ug/mL to umol
    u0_base.Nc_1 = Nseed * 3/4; 
    u0_base.Nc_1_neg = Nseed * 1/4;
    sol0 = solve(ODEProblem(invitro_model_her2!, u0_base, (0., simstime), p_base), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-18);
    return (sol0.u[end].Nc_1_neg + sol0.u[end].Nc_2_neg)/(sol0.u[end].Nc_1 + sol0.u[end].Nc_2 + sol0.u[end].Nc_1_neg + sol0.u[end].Nc_2_neg)
end

# bystander effect reported in Fig 2C; the desired outcome should be somewhere between 0.4 and 0.6 
p_her2pbd.k_in = 10.
bystander(0.4, 1E6, 120., p_her2pbd)

#======================= Save final parameters =======================# 
final_param_her2pbd_skbr3 = DataFrame(
    names = [string.(keys(p_her2pbd))[i] for i in 1:length(p_her2pbd)], 
    values = [p_her2pbd[i] for i in 1:length(p_her2pbd)]
)

using CSV
CSV.write(@projectroot("deliv/tab/dhes0815a_skbr3_params.csv"), final_param_her2pbd_skbr3);
