# date: 7/1/2025 
# author: Yuezhe Li 
# purpose of this code: to return organ volume 

function ReturnHumanOrganVolume(; N_Organs = 15)

    # Preallocate arrays
    V_endosomal = zeros(N_Organs)
    CL_up = zeros(N_Organs)
    V_VM = zeros(N_Organs)
    V_ISM = zeros(N_Organs)
    V_E7 = zeros(N_Organs)
    V_E6a = zeros(N_Organs)
    V_E7b = zeros(N_Organs)
    Organ_Endothelial_Cell = zeros(N_Organs)
    Endothelial_Cell_Frac = zeros(N_Organs)
    V_V = zeros(N_Organs)
    V_IntS = zeros(N_Organs)

    # Organ Indices
    Lung = 1
    Liver = 2
    Heart = 3
    Muscle = 4
    Skin = 5
    Adipose = 6
    Bone = 7
    Brain = 8
    Kidney = 9
    SI = 10
    LI = 11
    Pancreas = 12
    Thymus = 13
    Spleen = 14
    Other = 15
    Marrow = 16


    # Endothelial Cell Fractions
    Endothelial_Cell_Frac[Lung] = 0.0834
    Endothelial_Cell_Frac[Liver] = 0.1877
    Endothelial_Cell_Frac[Heart] = 0.0011
    Endothelial_Cell_Frac[Muscle] = 0.1928
    Endothelial_Cell_Frac[Skin] = 0.0819
    Endothelial_Cell_Frac[Adipose] = 0.0999
    Endothelial_Cell_Frac[Bone] = 0.1478
    Endothelial_Cell_Frac[Brain] = 0.0115
    Endothelial_Cell_Frac[Kidney] = 0.0157
    Endothelial_Cell_Frac[SI] = 0.0121
    Endothelial_Cell_Frac[LI] = 0.0209
    Endothelial_Cell_Frac[Pancreas] = 0.0001
    Endothelial_Cell_Frac[Thymus] = 0.0001
    Endothelial_Cell_Frac[Spleen] = 0.0499
    Endothelial_Cell_Frac[Other] = 0.094    # Modified after adding the BM compartment 0.0951 - 0.0011 = 0.094
    Endothelial_Cell_Frac[Marrow] = 0.0011   # similar to the heart, because it has a similar blood flow rate

    # Vascular Volumes [L]
    # based on male BW = 71kg
    V_V[Lung] = 55/1000
    V_V[Liver] = 183/1000
    V_V[Heart] = 13.1/1000
    V_V[Muscle] = 662/1000
    V_V[Skin] = 127/1000
    V_V[Adipose] = 148/1000
    V_V[Bone] = 224/1000
    V_V[Brain] = 31.9/1000
    V_V[Kidney] = 18.2/1000
    V_V[SI] = 6.15/1000
    V_V[LI] = 8.74/1000
    V_V[Pancreas] = 5.7/1000
    V_V[Thymus] = 0.353/1000
    V_V[Spleen] = 26.8/1000
    V_V[Other] = 54/1000     # modified after adding the BM compartment 204-150 = 54
    V_V[Marrow] = 150/1000    # From Baxter_et_al_1995    

    # Interstitial Volumes [L]
    # based on male BW = 71kg
    V_IntS[Lung] = 300/1000
    V_IntS[Liver] = 429/1000
    V_IntS[Heart] = 48.8/1000
    V_IntS[Muscle] = 3910/1000
    V_IntS[Skin] = 1125/1000
    V_IntS[Adipose] = 2289/1000
    V_IntS[Bone] = 1891/1000
    V_IntS[Brain] = 261/1000
    V_IntS[Kidney] = 49.8/1000
    V_IntS[SI] = 67.1/1000
    V_IntS[LI] = 95.3/1000
    V_IntS[Pancreas] = 18/1000
    V_IntS[Thymus] = 1.09/1000
    V_IntS[Spleen] = 44.3/1000
    V_IntS[Other] = 552/1000      # modified after adding the BM compartment 831-279 = 552
    V_IntS[Marrow] = 279/1000      # From Baxter_et_al_1995    

    pino_time = 10.8/60  # min->h
    Scale_Factor = 603.7 # scaling factor for human; if mouse, this number should be 1
    Total_Endothelial_Cell = 1.422e+009 # [number]. Total number of endothelial cells in mouse
    CL_up_in_nL_per_hour_per_million_cells = 150
    @. Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    @. CL_up  = CL_up_in_nL_per_hour_per_million_cells*Organ_Endothelial_Cell*1E-15

    tau_VM = 1/60  # [h]  
    tau_ISM = 1/60 # [h]

    E6a_Vol_Pct = 0.33 # [-]
    E7b_Vol_Pct = (1-E6a_Vol_Pct)/2 # [-]
    E7_Vol_Pct = (1-E6a_Vol_Pct)/2 # [-]      
        @. V_VM = CL_up * tau_VM
    @. V_ISM = CL_up * tau_ISM
    @. V_endosomal = CL_up_in_nL_per_hour_per_million_cells * pino_time*Organ_Endothelial_Cell*1e-6*1e-9
    @. V_E7 = V_endosomal * E7_Vol_Pct
    @. V_E6a = V_endosomal * E6a_Vol_Pct
    @. V_E7b = V_endosomal * E7b_Vol_Pct

    return V_V, V_VM, V_E7, V_E6a, V_E7b, V_ISM, V_IntS
end


# return total organ volume 
function TotalHumanOrganVolume(; N_Organs = 15)

    V_Organ = zeros(N_Organs)

    # Organ Indices
    Lung = 1
    Liver = 2
    Heart = 3
    Muscle = 4
    Skin = 5
    Adipose = 6
    Bone = 7
    Brain = 8
    Kidney = 9
    SI = 10
    LI = 11
    Pancreas = 12
    Thymus = 13
    Spleen = 14
    Other = 15
    Marrow = 16

    V_Organ[Lung] = 1000/1000 
    V_Organ[Liver] = 2143/1000
    V_Organ[Heart] = 341/1000
    V_Organ[Muscle] = 30078/1000
    V_Organ[Skin] = 3408/1000
    V_Organ[Adipose] = 13465/1000
    V_Organ[Bone] = 10165/1000
    V_Organ[Brain] = 1450/1000
    V_Organ[Kidney] = 332/1000
    V_Organ[SI] = 385/1000
    V_Organ[LI] = 548/1000
    V_Organ[Pancreas] = 104/1000
    V_Organ[Thymus] = 6.41/1000
    V_Organ[Spleen] = 221/1000
    V_Organ[Other] = 3352/1000    # modified after adding the BM compartment 4852-1500 = 3352
    V_Organ[Marrow] = 1500/1000   # From Baxter_et_al_1995  
    
    return V_Organ

end