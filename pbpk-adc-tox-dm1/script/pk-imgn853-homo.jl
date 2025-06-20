# date: 1/30/2025
# author: Yuezhe Li 
# purpose of this code: fit pk of mirvetuximab soravtansine (IMGN853)

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using DataFrames, CSV, DataFramesMeta
pk_moore = CSV.read(@projectroot("data/moore-2017.csv"),DataFrame);

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Plots, Statistics

include(@projectroot("model/jones_homo.jl") );
include(@projectroot("model/param-pk.jl")); 
include(@projectroot("model/init-pk.jl")); 
include(@projectroot("script/helper.jl")); 

p_imgn853 = deepcopy(p_base);
p_imgn853.PS_Score = 15          # fitted
p_imgn853.init_sR = 0.25E-3     # https://pubmed.ncbi.nlm.nih.gov/36402875/
p_imgn853.DAR = 3.5             # https://pmc.ncbi.nlm.nih.gov/articles/PMC10387675/
p_imgn853.k_out = 2.27          # https://pubmed.ncbi.nlm.nih.gov/16618769/
p_imgn853.Kd = 0.01E-3          # https://pubmed.ncbi.nlm.nih.gov/25904506/

q3w_dose = [0.15 0.5 1 2 3.3 5 7]; # [mg/kg]

tspan = (-0.01, 84*DayToHour);  # hr
dose_q3w = [0, 21, 42, 63]*DayToHour

df_ms = DataFrame(time = [], pl_lung = [], pl_liver = [], pl_heart = [], pl_muscle = [], pl_skin = [], pl_adipose = [], pl_bone = [], pl_brain = [], 
                  pl_kidney = [], pl_si = [], pl_li = [], pl_pancreas = [], pl_thymus = [], pl_spleen = [], adcplasma = [], dose = []);

for init_dose in q3w_dose
    # body weight was tuned for PK 
    tmpdf_q3w = InfusionDoses(init_dose, dose_q3w, p_imgn853, infusion_time = 0.1, BW = 120);
    tmpdf_q3w.dose .= string(init_dose) * " mg/kg Q3W"
    df_ms = vcat(df_ms, tmpdf_q3w)
end

# PK 
p_pk_imgn853 = plot(xlabel = "Time (day)", ylabel = "ADC concentration (ug/mL)", legend = :top, legendcolumns=2, palette = :lighttest, dpi = 1000, size = (400, 400)); 
plot!(df_ms.time/DayToHour, df_ms.adcplasma*MW/1E3, group = df_ms.dose, linewidth = 2); 
scatter!(pk_moore.time_d, pk_moore.adc_conc_ugmL, group = pk_moore.dose_mgkg, markerstrokewidth=0, markersize=6, alpha = 0.8, label=false);
plot!(yaxis = :log, ylims = [1, 1E3], yticks = [1, 10, 100, 1000], xlims = [0, 15]);
xticks!([0, 7, 14, 21], ["0", "7", "14", "21"]);
display(p_pk_imgn853)

savefig(p_pk_imgn853, @projectroot("deliv/figure/pk/imgn853-homo.png"));

# payload
df_ms = DataFrame(time = [], pl_lung = [], pl_liver = [], pl_heart = [], pl_muscle = [], pl_skin = [], pl_adipose = [], pl_bone = [], pl_brain = [], 
                  pl_kidney = [], pl_si = [], pl_li = [], pl_pancreas = [], pl_thymus = [], pl_spleen = [], adcplasma = [], dose = []);

q3w_dose_2 = [0.15 0.5 1 2 3.3 5 6 7]; # [mg/kg]

for init_dose in q3w_dose_2
    # body weight was tuned for PK 
    tmpdf_q3w = InfusionDoses(init_dose, dose_q3w, p_imgn853, infusion_time = 0.1);
    tmpdf_q3w.dose .= string(init_dose) * " mg/kg Q3W"
    df_ms = vcat(df_ms, tmpdf_q3w)
end

plot2_pl_he = plot(xlabel = "Time (day)", ylabel = "Payload concentration (uM)", dpi = 1000, size = (400, 400), ylims = [0, 3.5], xticks = 0:7:84); 
plot!(legend = :top, legendcolumns=2, legendfontsize = 8, xlabel = "Time (day)", ylabel = "DM4 conc (nM)");
plot!(df_ms.time/DayToHour, df_ms.pl_liver*1E3, group = df_ms.dose, palette = :RdYlGn_8, alpha = 0.8);
display(plot2_pl_he)

savefig(plot2_pl_he, @projectroot("deliv/figure/payload/imgn853-he-q3w.png"));

# additional dose that matches SAR408701
df_ms_170mgm2 = InfusionDoses(170*HT*HT/BW, dose_q3w, p_imgn853, infusion_time = 0.1);
df_ms_170mgm2.dose .= "170 mg/m2 Q3W";

CSV.write(@projectroot("data/sims/homo-imgn853.csv"), vcat(df_ms, df_ms_170mgm2))
