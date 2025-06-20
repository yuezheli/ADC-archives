# date: 1/24/2025
# author: Yuezhe Li 
# purpose of this code: to fit for PK of MEDI2228 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using Plots

# MEDI2228 homo data; 
# data obtained from https://www.tandfonline.com/doi/full/10.1080/10428194.2024.2373331
using CSV, DataFrames
obs = CSV.read(@projectroot("data/dimopoulos-2024-geomean.csv"),DataFrame);

include(@projectroot("model/jones_homo.jl") );
include(@projectroot("model/param-pk.jl")); 
include(@projectroot("model/init-pk.jl")); 
include(@projectroot("script/helper.jl")); 

p_medi2228 = deepcopy(p_base);
p_medi2228.PS_Score = 8        # fitted
p_medi2228.init_sR = 0.046      # https://pubmed.ncbi.nlm.nih.gov/30315237/
p_medi2228.DAR = 2              # https://www.tandfonline.com/doi/full/10.1080/10428194.2024.2373331
p_medi2228.k_out = 24           # https://pubs.acs.org/doi/10.1021/acs.jmedchem.0c00691
p_medi2228.Kd = 60.7E-3         # https://pubmed.ncbi.nlm.nih.gov/30315237/ # note the binding affinity (Kd) used the value for sBCMA:ADC, which is different from membrane BCMA:ADC


tspan = (-0.01, 42*DayToHour);  # hr
dose_q3w = [0, 21]*DayToHour

first_dose = [0.0125, 0.025, 0.05, 0.1, 0.14, 0.2]; # [mg/kg]

df_q3w = DataFrame(time = [], pl_lung = [], pl_liver = [], pl_heart = [], pl_muscle = [], pl_skin = [], pl_adipose = [], pl_bone = [], pl_brain = [], 
                  pl_kidney = [], pl_si = [], pl_li = [], pl_pancreas = [], pl_thymus = [], pl_spleen = [], adcplasma = [], dose = []);

for init_dose in first_dose
    tmpdf_q3w = InfusionDoses(init_dose, dose_q3w, p_medi2228, infusion_time = 0.001, Reltol = 1E-18);
    tmpdf_q3w.dose .= string(init_dose) * " mg/kg"
    df_q3w = vcat(df_q3w, tmpdf_q3w)
end

# PK 

p_pk_medi2228 = plot(xlabel = "Time (day)", ylabel = "ADC concentration (ug/L)", legend = :bottomright, legendcolumns=2, palette = :Accent_6, dpi = 1000, size = (400, 400)); 
plot!(df_q3w.time/DayToHour, df_q3w.adcplasma*MW, group = df_q3w.dose, linewidth = 2, alpha = 0.8); 
scatter!(obs.time_d, obs.adc_conc_ngmL, group = obs.dose_mgkg, markerstrokewidth=0, markersize=6, alpha = 0.5, label=false);
plot!(yaxis = :log, ylims = [1, 15000], xlims = [0, 42]);
xticks!(0:7:42);
display(p_pk_medi2228)

savefig(p_pk_medi2228, @projectroot("deliv/figure/pk/MEDI2228-homo.png"));

# payload 
CSV.write(@projectroot("data/sims/homo-medi2228.csv"), df_q3w)

pskin = plot(df_q3w.time/DayToHour, df_q3w.pl_skin, group = df_q3w.dose, xticks = 0:7:42, title = "skin endothelial payload, MEDI2228")
pliver = plot(df_q3w.time/DayToHour, df_q3w.pl_liver, group = df_q3w.dose, xticks = 0:7:42, title = "liver endothelial payload, MEDI2228")

savefig(pskin, @projectroot("deliv/figure/payload/medi2228-q3w-skin.png"));
savefig(pliver, @projectroot("deliv/figure/payload/medi2228-q3w-liver.png"));

