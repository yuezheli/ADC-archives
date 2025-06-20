# date: 1/22/2025
# author: Yuezhe Li 
# purpose of this code: to fit homo PK of ADCT-601

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using Plots

# ADCT-601 homo data; 
# data obtained from Van Tine et al., 2024; # https://www.adctmedical.com/wp-content/uploads/2024/04/AACR-2024_Van-Tine_ADCT-601-102-dose-esc-STS_poster_FINAL_crtd-1.pdf
using DataFrames
adct_601_homo_7point5 = DataFrame(time_d=[0, 0, 1.1, 3.8, 20.8, 22.6, 24.3], ADC_conc_ugL=[8139.9, 2187.2, 1000.9, 211.9, 3508, 1364.2, 578.7]); 
adct_601_homo_7point5.Dose_mg .= 7.5; 
adct_601_homo_11 = DataFrame(time_d=[3.8, 0, 0, 0, 1.5, 3.7, 21, 21, 21, 22.4, 24.3, 42], ADC_conc_ugL=[211.9, 1983.6, 1685.3, 1312.8, 527.4, 168.7, 15.4, 1828.3, 1006.1, 584.7, 566.3, 15.6]); 
adct_601_homo_11.Dose_mg .= 11; 
adct_601_homo_13 = DataFrame(time_d=[0,0,0,1.5,21,21,21,22.5,24.2,42], ADC_conc_ugL=[2900.9, 1798.8, 1579, 1001, 15.5, 2060.3, 3368.9, 1379.1, 273.5, 15.8]); 
adct_601_homo_13.Dose_mg .= 13;
adct_601_homo_15 = DataFrame(time_d=[0,0,0,1.7,3.6,21,21,21,22.3,24.4,42], ADC_conc_ugL=[989.8,3029.6,3764.6,1632,601.1,14.9,758.5,3079.3,1588.1,491.7,16]); 
adct_601_homo_15.Dose_mg .= 15;  
adct_601_homo = vcat(adct_601_homo_7point5, adct_601_homo_11, adct_601_homo_13, adct_601_homo_15);

include(@projectroot("model/jones_homo.jl") );
include(@projectroot("model/param-pk.jl")); 
include(@projectroot("model/init-pk.jl")); 
include(@projectroot("script/helper.jl")); 

p_adct601 = deepcopy(p_base);
p_adct601.PS_kd = 0.01        # fitted
p_adct601.PS_Score = -2       # tuned off
p_adct601.KD6_WT = 145        # https://pmc.ncbi.nlm.nih.gov/articles/PMC9377743/
p_adct601.init_sR = 0.6E-3    # https://pmc.ncbi.nlm.nih.gov/articles/PMC6952099/
p_adct601.DAR = 1.9           # https://www.adctherapeutics.com/wp-content/uploads/2020/09/Poster_601_AACR-2018-.pdf
p_adct601.k_out = 24          # https://pubs.acs.org/doi/10.1021/acs.jmedchem.0c00691
p_adct601.Kd = 0.211E-3       # https://pubmed.ncbi.nlm.nih.gov/33046521/

tspan = (-0.01, 42*DayToHour);  # hr
dose_q3w = [0, 21]*DayToHour

first_dose = [7.5, 11, 13, 15]; # [mg]

df_q3w = DataFrame(time = [], pl_lung = [], pl_liver = [], pl_heart = [], pl_muscle = [], pl_skin = [], pl_adipose = [], pl_bone = [], pl_brain = [], 
                  pl_kidney = [], pl_si = [], pl_li = [], pl_pancreas = [], pl_thymus = [], pl_spleen = [], adcplasma = [], dose = []);

for init_dose in first_dose
    tmpdf_q3w = InfusionDoses(init_dose/BW, dose_q3w, p_adct601, infusion_time = 0.5, Reltol = 1E-18);
    tmpdf_q3w.dose .= string(init_dose) * " mg"
    df_q3w = vcat(df_q3w, tmpdf_q3w)
end

# PK 

p_pk_adct601 = plot(xlabel = "Time (day)", ylabel = "ADC concentration (ug/L)", legend = :bottomright, legendcolumns=2, palette = :Set1_4, dpi = 1000, size = (400, 400)); 
plot!(df_q3w.time/DayToHour, df_q3w.adcplasma*MW, group = df_q3w.dose, linewidth = 2, alpha = 0.8); 
scatter!(adct_601_homo.time_d, adct_601_homo.ADC_conc_ugL, group = adct_601_homo.Dose_mg, markerstrokewidth=0, markersize=6, alpha = 0.5, label=false);
plot!(yaxis = :log, ylims = [1, 20000], xlims = [0, 42]);
xticks!([0, 7, 14, 21, 28, 35, 42], ["0", "7", "14", "21", "28", "35", "42"]);
display(p_pk_adct601)

savefig(p_pk_adct601, @projectroot("deliv/figure/pk/ADCT-601-homo.png"));

# payload 
using CSV
CSV.write(@projectroot("data/sims/homo-adct601.csv"), df_q3w)

pskin = plot(df_q3w.time/DayToHour, df_q3w.pl_skin, group = df_q3w.dose, xticks = 0:7:63, title = "skin endothelial payload, ADCT-601")
pliver = plot(df_q3w.time/DayToHour, df_q3w.pl_liver, group = df_q3w.dose, xticks = 0:7:63, title = "liver endothelial payload, ADCT-601")

savefig(pskin, @projectroot("deliv/figure/payload/adct-601-q3w-skin.png"));
savefig(pliver, @projectroot("deliv/figure/payload/adct-601-q3w-liver.png"));
