# date: 1/27/2025
# author: Yuezhe Li 
# purpose of this code: to show PK fit of ADCT-402

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using Plots

# ADCT-402 homo data; 
# data obtained from https://www.accessdata.fda.gov/drugsatfda_docs/nda/2021/761196Orig1s000MultidisciplineR.pdf
using DataFrames
obs = DataFrame(
    time_week = [0.05, 0.44, 1.1, 1.87, 2.97, 3.02, 3.24, 3.63, 4.78, 6.04], 
    adc_conc_ugmL = [1.6, 1.02, 0.75, 0.58, 0.44, 2.22, 1.79, 1.36, 0.97, 0.71]
);

include(@projectroot("model/jones_homo.jl") );
include(@projectroot("model/param-pk.jl")); 
include(@projectroot("model/init-pk.jl")); 
include(@projectroot("script/helper.jl")); 

p_adct402 = deepcopy(p_base);
p_adct402.PS_Score = 6         # fitted from cyno
p_adct402.init_sR = 0.77E-3    # https://pmc.ncbi.nlm.nih.gov/articles/PMC3520838/
p_adct402.DAR = 2.3            # https://www.accessdata.fda.gov/drugsatfda_docs/nda/2021/761196Orig1s000MultidisciplineR.pdf (BLA)
p_adct402.k_out = 24           # https://pubs.acs.org/doi/10.1021/acs.jmedchem.0c00691 (computed from PAMPA value reported in Table 1, with the assumption that cell radius = 5um; the range for this rate from different PBDs should be between 2.4 hr-1 and 345.6 hr-1, depends on the pH and modification )
p_adct402.Kd = 665E-6          # https://www.accessdata.fda.gov/drugsatfda_docs/nda/2021/761196Orig1s000MultidisciplineR.pdf (BLA)

tspan = (-0.01, 42*DayToHour);  # hr
dose_q3w = [0, 21]*DayToHour

df_150_q3w = InfusionDoses(0.15, dose_q3w, p_adct402, infusion_time = 0.1, Reltol = 1E-18);
df_150_q3w.dose .= "150 ug/kg";

p_pk_adct402 = plot(xlabel = "Time (day)", ylabel = "ADC concentration (ug/L)", legend = :bottomright, legendcolumns=1, dpi = 1000, size = (400, 400)); 
plot!(df_150_q3w.time/DayToHour, df_150_q3w.adcplasma*MW, linewidth = 2, alpha = 0.8, label = "simulation"); 
scatter!(obs.time_week * 7, obs.adc_conc_ugmL*1E3, markerstrokewidth=0, markersize=6, alpha = 0.5, label="obs, 150 ug/kg");
plot!(yaxis = :log, ylims = [1, 20000], xlims = [0, 21]);
xticks!([0, 7, 14, 21], ["0", "7", "14"]);
display(p_pk_adct402)

savefig(p_pk_adct402, @projectroot("deliv/figure/pk/ADCT-402-homo.png"));

# payload 
pskin = plot(df_150_q3w.time/DayToHour, df_150_q3w.pl_skin, xticks = 0:7:63, title = "skin endothelial payload, ADCT-402")
pliver = plot(df_150_q3w.time/DayToHour, df_150_q3w.pl_liver, xticks = 0:7:63, title = "liver endothelial payload, ADCT-402")

savefig(pskin, @projectroot("deliv/figure/payload/adct-402-q3w-skin.png"));
savefig(pliver, @projectroot("deliv/figure/payload/adct-402-q3w-liver.png"));

df_200_q3w = InfusionDoses(0.2, dose_q3w, p_adct402, infusion_time = 0.1, Reltol = 1E-18);
df_200_q3w.dose .= "200 ug/kg";

df_q3w = vcat(df_150_q3w, df_200_q3w);

using CSV
CSV.write(@projectroot("data/sims/homo-adct402.csv"), df_q3w)
