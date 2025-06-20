# date: 1/28/2025
# author: Yuezhe Li 
# purpose of this code: to fit PK for PCA062

using ProjectRoot
using Pkg; Pkg.activate(@projectroot)
using DataFrames, CSV, DataFramesMeta

# read observed; https://pubmed.ncbi.nlm.nih.gov/35131875/
pk_duca =  CSV.read("data/duca-2022.csv",DataFrame);


using DifferentialEquations, RecursiveArrayTools, ComponentArrays
using Plots

# simulation 
include(@projectroot("model/jones_homo.jl") );
include(@projectroot("model/param-pk.jl")); 
include(@projectroot("model/init-pk.jl")); 
include(@projectroot("script/helper.jl"));

p_pca062 = deepcopy(p_base);
p_pca062.PS_Score = 8        # fitted
p_pca062.init_sR = 0.       # assumed 

tspan = (-0.01, 350.);      # [hr]

q2w_dose = [0.4, 0.6, 0.9, 1.4, 2.1, 3.6, 4.4, 4.7]; # [mg/kg]

df_q2w = DataFrame(time = [], pl_lung = [], pl_liver = [], pl_heart = [], pl_muscle = [], pl_skin = [], pl_adipose = [], pl_bone = [], pl_brain = [], 
                  pl_kidney = [], pl_si = [], pl_li = [], pl_pancreas = [], pl_thymus = [], pl_spleen = [], adcplasma = [], dose = []);

for init_dose in q2w_dose
    tmpdf = InfusionDoses(init_dose, [0.], p_pca062, infusion_time = 0.2);
    tmpdf.dose .= string(init_dose) * " mg/kg"
    df_q2w = vcat(df_q2w, tmpdf)
end

p_pk_pca062 = plot(xlabel = "Time (hour)", ylabel = "ADC concentration (ug/mL)", legend = :bottomright, legendcolumns=2, palette = :Set2_8, dpi = 1000, size = (400, 400)); 
plot!(df_q2w.time, df_q2w.adcplasma*MW/1E3, group = df_q2w.dose, linewidth = 2, alpha = 0.8); 
scatter!(pk_duca.time_hr, pk_duca.adc_conc_ugmL, group = pk_duca.dose_mgkg, markerstrokewidth=0, markersize=6, alpha = 0.5, label=false);
plot!(yaxis = :log, ylims = [0.01, 1000], yticks = [0.01, 0.1, 1, 10, 100, 1000], xlims = [0, 350]);
display(p_pk_pca062)

savefig(p_pk_pca062, @projectroot("deliv/figure/pk/PCA062-homo.png"));


# payload 
CSV.write(@projectroot("data/sims/homo-pca062.csv"), df_q2w)

pliver = plot(df_q2w.time/DayToHour, df_q2w.pl_liver, group = df_q2w.dose, xticks = 0:7:14, title = "liver endothelial payload, PCA062")

savefig(pliver, @projectroot("deliv/figure/payload/pca062-q2w-liver.png"));
