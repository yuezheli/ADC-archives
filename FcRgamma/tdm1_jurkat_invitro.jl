# author: Yuezhe Li 
# date: May 23, 2024
# purpose of this code: to fit for FcγR intake and cytotoxicity in T cell line (Jurkat) and megakaryoblastic leukemia cell line (MEG01-S)
# data obtained from Aoyama et al., 2022; # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8837541/

using Pkg; Pkg.activate("");

# load ODE model 
include("model/adc_fcgammar.jl"); 

using DifferentialEquations, ComponentArrays
using Plots, Statistics, DataFrames

# T-DM1 uptake by Jurkat cells (FcγR-mediated) (Fig 3a middle panel of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8837541/)
tdm1_uptake = DataFrame(
    tdm1_conc_nM = [5.47, 13., 33.31, 82., 210., 520],
    lumin_ints = [3.97, 3.97, 4.63, 4.63, 5.29, 8.6] 
); 

# T-DM1 cytotoxicity by Jurkat cells (FcγR-mediated) (Fig 3b middle panel of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8837541/)
tdm1_cell_survival = DataFrame(
    tdm1_conc_nM = [0.12, 0.76, 1.94, 5., 12., 32., 78., 201.32, 515], 
    cell_survival = [0.86, 0.85, 0.89, 0.91, 0.89, 0.88, 0.62, 0.24, 0.12]
);

# T-Dxd cytotoxicity by Jurkat cells (FcγR-mediated) (Fig 3b middle right panel of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8837541/)
tdxd_cell_survival = DataFrame(
    tdxd_conc_nM = [0.13, 0.33, 0.83, 2.02, 4.9, 12.48, 30.29, 77.21, 196.84, 477.69], 
    cell_survival = [0.89, 0.87, 0.88, 0.83, 0.86, 0.85, 0.84, 0.79, 0.68, 0.29]
); 

p_tdm1 = ComponentArray(
    tdouble = 20.7,                         # doubling time, [hr], https://www.cellosaurus.org/CVCL_0367
    DAR = 3.5,                              # average drug : antibody ratio, [unitless], https://www.nature.com/articles/s41416-019-0635-y
    IC50_Payload = 4.,                      # payload IC50, [uM], https://pubmed.ncbi.nlm.nih.gov/37799390/
    Emax_Payload = 0.139,                   # maximum in vitro killing rate, [1/hr]
    k_in_FcgR = 0.1,                        # FcγR-mediated ADC uptake
    k_out = 0.14,                           # rate of payload leaving of the cell through diffusion, [1/hr], https://europepmc.org/article/ppr/ppr584999
    k_PL = 0.34,                            # free payload degredation rate, [1/hr], 
    k_in = 0.21,                            # rate of payload entering cell through diffusion, [1/hr], https://europepmc.org/article/ppr/ppr584999
    k_lys = 2.4,                            # rate for ADC to be graded in cells to release payload, [1/hr], computed from endosome/ lysosome half life https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2112524/pdf/jc972508.pdf
    V_medium = 1E-4,                        # medium volume, [L]
    Nmax = 4E5,                             # maximum cell per well, https://www.thermofisher.com/us/en/home/references/gibco-cell-culture-basics/cell-culture-protocols/cell-culture-useful-numbers.html
); 

p_tdxd = deepcopy(p_tdm1);
p_tdxd.DAR = 8.
p_tdxd.IC50_Payload = 0.31; # https://www.medchemexpress.com/Dxd.html
p_tdxd.k_out = 32.32; # https://europepmc.org/article/ppr/ppr584999
p_tdxd.k_in = 46.08; # https://europepmc.org/article/ppr/ppr584999
p_tdxd.k_PL = 0.82; # https://europepmc.org/article/ppr/ppr584999

# T-DM1 intenalization 
function tdm1_internalization_loss(p, obs, opt = true)
    p_tmp = deepcopy(p_tdm1);
    p_tmp.k_in_FcgR = p[1]
    p_tmp.Emax_Payload = p[2]
    p_tmp.k_out = p[3]
    cyto_adc = [];
    for adcconc in obs.tdm1_conc_nM
        u0_tmp = ComponentArray(A_m = adcconc * 1E-3 * p_tmp.V_medium, A_c = 0., Nc = 1E5, P_c = 0., P_m = 0.);
        sol_tmp = solve(ODEProblem(adc_fcgammar!, u0_tmp, (0., 4.), p_tmp), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-18);
        append!(cyto_adc, sol_tmp[:A_c][end])
    end
    if opt
        normalized_adc = (obs.lumin_ints .- obs.lumin_ints[1]) ./ (obs.lumin_ints[end] - obs.lumin_ints[1])
        normalized_sim_adc = (cyto_adc .- cyto_adc[1]) ./ (cyto_adc[end] - cyto_adc[1])
        return sum( (normalized_adc .- normalized_sim_adc).^2 )
    else 
        return (cyto_adc .- cyto_adc[1]) ./ (cyto_adc[end] - cyto_adc[1]) 
    end
end

ints = tdm1_internalization_loss([0.0001, 0.05, 0.5], tdm1_uptake, false);  

plt_int = plot(xlabel = "T-DM1 conc (nM)", ylabel = "Normalized Lumin Intensity", label = :outerright);
plot!(tdm1_uptake.tdm1_conc_nM, (tdm1_uptake.lumin_ints .- tdm1_uptake.lumin_ints[1]) ./ (tdm1_uptake.lumin_ints[end] - tdm1_uptake.lumin_ints[1]), label = "Aoyama et al., 2022", seriestype=:scatter);
plot!(tdm1_uptake.tdm1_conc_nM, ints, label = "fitted");
display(plt_int); 

savefig(plt_int, "figure/tdm1_junkat_internalization.png");

# ADC induced cytotocixity 
function adc_cytotoxicity(p, obs, opt = true, p_ADC = p_tdm1)
    p_tmp = deepcopy(p_ADC);
    p_tmp.k_in_FcgR = p[1]
    p_tmp.Emax_Payload = p[2]
    p_tmp.k_out = p[3]
    # no dose 
    Nc1 = [];
    for adcconc in obs[:,1]
        u0_tmp = ComponentArray(A_m = adcconc * 1E-3 * p_tmp.V_medium, A_c = 0., Nc = 1E5, P_c = 0., P_m = 0.);
        sol_tmp = solve(ODEProblem(adc_fcgammar!, u0_tmp, (0., 72.), p_tmp), saveat = 1., alg = QNDF(autodiff=false), reltol = 1e-18);
        append!(Nc1, sol_tmp[:Nc][end])
    end
    if opt
        return sum( (Nc1 / Nc1[1] .- obs.cell_survival).^2 )
    else 
        return Nc1 / Nc1[1]
    end
end

## T-DM1 induced cytotocixity 
survival = adc_cytotoxicity([0.0001, 0.12, 0.14], tdm1_cell_survival, false, p_tdm1); 

plt_cyt = plot(xlabel = "T-DM1 conc (nM)", ylabel = "Cell survival", label = :outerright, xaxis = :log10, legend = :bottomleft);
plot!(tdm1_cell_survival.tdm1_conc_nM, tdm1_cell_survival.cell_survival, label = "Aoyama et al., 2022", seriestype=:scatter);
plot!(tdm1_cell_survival.tdm1_conc_nM, survival, label = "fitted");
display(plt_cyt); 

savefig(plt_cyt, "figure/tdm1_junkat_toxicity.png");

## T-Dxd induced cytotocixity 
survival_tdxd = adc_cytotoxicity([0.0001, 0.55, 46.08], tdxd_cell_survival, false, p_tdxd); 

plt_cyt_tdxd = plot(xlabel = "T-Dxd conc (nM)", ylabel = "Cell survival", label = :outerright, xaxis = :log10, legend = :bottomleft, ylims = [0.,1]);
plot!(tdxd_cell_survival.tdxd_conc_nM, tdxd_cell_survival.cell_survival, label = "Aoyama et al., 2022", seriestype=:scatter);
plot!(tdxd_cell_survival.tdxd_conc_nM, survival_tdxd, label = "fitted");
display(plt_cyt_tdxd); 

savefig(plt_cyt_tdxd, "figure/tdxd_junkat_toxicity.png");


