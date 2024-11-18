# author: Yuezhe Li 
# date: Nov 17, 2023
# purpose: calculate glomerular sieving coefficient (theta)
# the implmentation of the equation was based on Li and Shah, 2019 (Eq 8), but this equation is probably wrong since higher MW was linked to higher renal clearance. 
# further investigation is required to generate sth like Fig 6.

using Pkg; Pkg.activate("..")

using Plots
using DataFrames
using CSV

# GFR; mouse, 28g
GFR_mus = 0.278; # [mL/min]

# this is the function obtained from Li and Shah, 2019
function gfr_siev_coef(MW_kDa)
    p1 = 1 + exp(-0.028 * (-MW_kDa + 72.3));
    p2 = 1 - 8.7/p1;
    p3 = exp(p2);
    return p3;
end

mw_kda = [13., 27., 50., 55., 80., 100., 150.];
renal_clearance = GFR_mus * [ gfr_siev_coef(tmp_mw) for tmp_mw in mw_kda ];

plot(mw_kda, renal_clearance, ylims =(0.00001, 1), yaxis = :log)

# data obtained from Haraldsson et al., 2008
# https://pubmed.ncbi.nlm.nih.gov/18391170/

haraldsson = DataFrame(
    names = ["nAlbumin", "nAlbumin", "Kappa-dimer", "nHRP", "nHRP", "nMyoglobin", "cMyoglobin", "IgG", "LDH-5"], 
    mol_weight = [66.5, 66.5, 22.5, 44., 44., 17., 17., 150., 140.], # kDa
    sieving_coef = [0.033, 0.026, 0.149, 0.07, 0.11, 0.77, 0.74, 0.0023, 0.0056]
);

# this is adjusted based on known relationship
function gfr_siev_coef2(MW_kDa)
    p1 = 1 + exp(0.028 * (-MW_kDa + 72.3));
    p2 = 1 - 8.7/p1;
    p3 = exp(p2);
    return p3;
end

sievcoef2 = [ gfr_siev_coef2(tmp_mw) for tmp_mw in haraldsson.mol_weight];

p_glo_siev_coef = plot(yaxis = :log, xlabel = "molecular weight (kDa)", ylabel = "glomerular filter coefficient", ylims =(0.0001, 10));
plot!(haraldsson.mol_weight, haraldsson.sieving_coef, label = "Haraldsson et al., 2008", seriestype=:scatter);
plot!(haraldsson.mol_weight, sievcoef2, label = "corrected formula", seriestype=:scatter);

savefig(p_glo_siev_coef, "../deliv/figure/glomerular-sieving-coef-calc.png");

# predict glomerular sieving coefficient (theta) based on the updated fitting 
sievcoef3 =  [ gfr_siev_coef2(tmp_mw) for tmp_mw in mw_kda ];

pred_theta = DataFrame(names = ["nanbody", "scFv", "Fab", "scFv2", "minibody", "Fab2", "IgG"], 
          MW = mw_kda, 
          filtration_coef = sievcoef3, 
          Li2019theta = [0.667, 0.403, 0.131, 0.0986, 0.022, 0.00703, 0.0011], 
          ); 

CSV.write("../deliv/table/glomerular-sieving-coefficient.csv",  pred_theta); 

# visualization of glomerular sieving coef (corrected) and molecular weight 
protein_mw_kda = 10 .^ LinRange(-1.1, 3.1, 50);
theta = [ gfr_siev_coef2(tmp_mw) for tmp_mw in protein_mw_kda ];

p_pred_theta = plot(ylims =(0.0001, 1.1), yaxis = :log, xaxis = :log, xlabel = "molecular weight (kDa)", ylabel = "glomerular sieving coefficient");
plot!(protein_mw_kda, theta, label = false);
#hline!([1], linestyle = :dot, alpha = 0.1);
xticks!([0.1, 1, 10, 100, 1000], ["0.1", "1", "10", "100", "1000"]);
p_pred_theta

savefig(p_pred_theta, "../deliv/figure/pred-glomerular-sieving-coef.png");
