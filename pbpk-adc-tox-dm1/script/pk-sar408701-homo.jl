# date: 1/6/2025
# author: Yuezhe Li 
# purpose of this code: fit pk of tusamitamab ravtansine -- 1 grade 3 increase in transaminase levels during cycle 1 with 190 mg/m2 Q3W
# https://pmc.ncbi.nlm.nih.gov/articles/PMC10461573/
# PK data obtained from https://pmc.ncbi.nlm.nih.gov/articles/PMC10461573/
# DM4 IC50 between 0.12nM and 23nM, depends on cell line and metabolite form; https://pubs.acs.org/doi/10.1021/jm060319f

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using DataFrames, CSV, DataFramesMeta
pk_tabernero = CSV.read(@projectroot("data/tabernero-2023.csv"),DataFrame);

using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Plots, Statistics

include(@projectroot("model/jones_homo.jl") );
include(@projectroot("model/param-pk.jl")); 
include(@projectroot("model/init-pk.jl")); 
include(@projectroot("script/helper.jl")); 

p_sar408701 = deepcopy(p_base);
p_sar408701.PS_Score = 8
p_sar408701.init_sR = 120/180E3 # https://onlinelibrary.wiley.com/doi/10.1111/cas.12451
p_sar408701.thalf_sR_adc = 105. # https://pubmed.ncbi.nlm.nih.gov/9100477/
p_sar408701.DAR = 3.8           # https://pubmed.ncbi.nlm.nih.gov/33046521/
p_sar408701.k_out = 2.27        # https://pubmed.ncbi.nlm.nih.gov/16618769/
p_sar408701.Kd = 0.017E-3       # https://pubmed.ncbi.nlm.nih.gov/33046521/

q3w_dose = [120 135 150 170 190]; # [mg/m2]
q2w_dose = [5 10 20 40 80 100 120 150]; # [mg/m2]

tspan = (-0.01, 84*DayToHour);  # hr
dose_q3w = [0, 21, 42, 63]*DayToHour
dose_q2w = [0, 14, 28, 42, 56, 70]*DayToHour

df_tr = DataFrame(time = [], pl_lung = [], pl_liver = [], pl_heart = [], pl_muscle = [], pl_skin = [], pl_adipose = [], pl_bone = [], pl_brain = [], 
                  pl_kidney = [], pl_si = [], pl_li = [], pl_pancreas = [], pl_thymus = [], pl_spleen = [], adcplasma = [], dose = [], scheme = []);

for init_dose in q3w_dose
    tmpdf_q3w = InfusionDoses(init_dose*HT*HT/BW, dose_q3w, p_sar408701, infusion_time = 3);
    tmpdf_q3w.dose .= string(init_dose) * " mg/m2 Q3W"
    tmpdf_q3w.scheme .= "Q3W";
    df_tr = vcat(df_tr, tmpdf_q3w)
    
end

for init_dose in q2w_dose
    tmpdf_q2w = InfusionDoses(init_dose*HT*HT/BW, dose_q2w, p_sar408701, infusion_time = 3);
    tmpdf_q2w.dose .= string(init_dose) * " mg/m2 Q2W"
    tmpdf_q2w.scheme .= "Q2W";
    df_tr = vcat(df_tr, tmpdf_q2w)
end

df_q3w = @rsubset(df_tr, :scheme .== "Q3W");

# PK 
p_pk_sar408701 = plot(xlabel = "Time (day)", ylabel = "ADC concentration (ug/mL)", legend = :top, legendcolumns=2, palette = :Set2_5, dpi = 1000, size = (400, 400)); 
plot!(df_q3w.time/DayToHour, df_q3w.adcplasma*MW/1E3, group = df_q3w.dose, linewidth = 2); 
scatter!(pk_tabernero.time_day, pk_tabernero.ADC_conc_ugmL, group = pk_tabernero.dose_mgm2, markerstrokewidth=0, markersize=6, alpha = 0.8, label=false);
plot!(yaxis = :log, ylims = [10, 200], xlims = [0, 15]);
xticks!([0, 7, 14], ["0", "7", "14"]);
display(p_pk_sar408701)

savefig(p_pk_sar408701, @projectroot("deliv/figure/pk/sar408701-homo.png"));

# PAYLOAD 

CSV.write(@projectroot("data/sims/homo-sar408701.csv"), df_tr)

plot_pl_he = plot(xlabel = "Time (day)", ylabel = "Payload concentration (uM)", dpi = 1000, size = (400, 400), ylims = [0, 3.5], xticks = 0:7:84 ); 
plot!(legend = :top, legendcolumns=3, legendfontsize = 6, xlabel = "Time (day)", ylabel = "DM4 conc (nM)");
plot!(df_q3w.time/DayToHour, df_q3w.pl_liver*1E3, group = df_q3w.dose, palette = :Set2_5);
display(plot_pl_he)

savefig(plot_pl_he, @projectroot("deliv/figure/payload/sar408701-q3w-liver.png"));

