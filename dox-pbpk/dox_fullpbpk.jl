# date: 4/1/2024
# author: Yuezhe Li 
# purpose of this script: to code up the full PBPK model for doxorubicin
# model published in He et al., 2019; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6533104/

using ComponentArrays
using Parameters: @unpack

function dox_fullpbpk!(du,u,p,t)
    @unpack Ca, Cp, Ce_blood, Cet_blood, Ci_org, Cet_org, totalclearance = u

    @unpack infusion, Kpt_blood, Va, Vp, Ve_blood, CLrenal, CLhepatic, PER = p

    N_Organs = 8
    Lung = 1
    Liver = 2
    Gut = 3
    Spleen = 4
    Kidney = 5
    Heart = 6
    Other = 7
    Tumor = 8

    # DNA concentration, 70kg human [umol/mL]
    CN = zeros(N_Organs)
    CN[Lung] = 23.5
    CN[Liver] = 6.43
    CN[Gut] = 20.5
    CN[Spleen] = 54.1
    CN[Kidney] = 12.5
    CN[Heart] = 4.59
    CN[Other] = 31.9
    CN[Tumor] = 17.

    # surface area of blood cells, 70kg human [cm^2]
    S = zeros(N_Organs)
    S[Lung] = 15.6
    S[Liver] = 28.5
    S[Gut] = 32.5
    S[Spleen] = 2.75
    S[Kidney] = 5.19
    S[Heart] = 6.03
    S[Other] = 226.
    S[Tumor] = 0.2
    S = S*1E5
    S_blood = 25.8E5

    # intracellular volume, 70kg human [mL]
    Ve_org = zeros(N_Organs)
    Ve_org[Lung] = 629.
    Ve_org[Liver] = 1153.
    Ve_org[Gut] = 1315.
    Ve_org[Spleen] = 111.
    Ve_org[Kidney] = 210.
    Ve_org[Heart] = 244.
    Ve_org[Other] = 9142.
    Ve_org[Tumor] = 8.

    # blood flow of organs, 70kg human [mL/h] 
    Q = zeros(N_Organs)
    Q[Lung] = 5600.
    Q[Liver] = 273. 
    Q[Gut] = 1100.
    Q[Spleen] = 77.
    Q[Kidney] = 1240.
    Q[Heart] = 240.
    # Q[Other] = (700. + 750 + 260 + 300) # modified based on https://link.springer.com/article/10.1023/A:1018943613122
    Q[Tumor] = 1.
    Q[Other] = (Q[Lung] - sum(Q[Liver:Heart]) - Q[Tumor])
    Q = Q*60

    # interstitium volume, 70kg human [mL]
    Vi_org = zeros(N_Organs)
    Vi_org[Lung] = 234.
    Vi_org[Liver] = 275.
    Vi_org[Gut] = 287.
    Vi_org[Spleen] = 38.
    Vi_org[Kidney] = 42.
    Vi_org[Heart] = 44.
    Vi_org[Other] = 1969.
    Vi_org[Tumor] = 2.

    # partition coef between interstitial fluid and plasma, 70kg human [-]
    Kpt = ones(N_Organs) * 0.68
    Kpt[Gut] = 0.94
    Kpt[Tumor] = 1.

    # cell/ interstitial fluid partition coefficient (nonspecific protein binding) (optimized in mouse)
    Kp = zeros(N_Organs)
    Kp[Lung] = 3.38
    #Kp[Liver] = 10.9
    Kp[Liver] = 20.
    Kp[Gut] = 1.26
    Kp[Spleen] = 5.74
    #Kp[Kidney] = 5.91
    Kp[Kidney] = 8
    #Kp[Heart] = 5.63
    Kp[Heart] = 15
    Kp[Other] = 13.3
    Kp[Tumor] = 3.6

    # other params 
    Kpp = 2.2836 # partition coefficient associated with pH gradient
    Fr = 0.18
    Fra = 1  # a parameter reflects the joint effect of influx/efflux transporters
    # PER = 0.0756  # cytoplasmic membrane permeability coefficient [cm/h]
    Kd = 0.13 # equilibrium dissociation constant of DNA intercalation [umol/mL]
    # CLrenal = 8190.748 # renal clearance, scaled from mouse, [mL/h]
    # CLhepatic = 25.5E3 # hepatic clearance, [mL/h]
    fup = 0.26

    # nucleus sub-compartments
    Ce_org = zeros(N_Organs); 
    @. Ce_org = 0.5 * ( (Cet_org .- CN .- Kd ) + ( (Cet_org .- CN .- Kd).^2 .+ 4*Kd*Cet_org ).^(0.5) )

    # interstitial subcompartment (modified based on past experience)
    @. du.Ci_org[Gut:Tumor] = ( (Ca .- Ci_org[Gut:Tumor]./Kpt[Gut:Tumor]).*Q[Gut:Tumor] .- PER*S[Gut:Tumor].*Ci_org[Gut:Tumor] .+ PER*S[Gut:Tumor].*Ce_org[Gut:Tumor]./(Kp[Gut:Tumor]*Kpp) ) ./ Vi_org[Gut:Tumor]

    # interstitial fluid in liver (modified based on past experience)
    du.Ci_org[Liver] = ( Ca*Q[Liver] + Ci_org[Gut]/Kpt[Gut]*Q[Gut] + Ci_org[Spleen]/Kpt[Spleen]*Q[Spleen] - Ci_org[Liver]/Kpt[Liver]*(Q[Liver]+Q[Gut]+Q[Spleen])  
                         - PER*S[Liver]*Ci_org[Liver] + PER*S[Liver]*Ce_org[Liver]/(Kp[Liver]*Kpp)
                         - CLhepatic*Ci_org[Liver]) / Vi_org[Liver]

    # interstitial fluid in lung (modified based on the past experience)
    du.Ci_org[Lung] = ( Cp*Q[Lung] - Ci_org[Lung]/Kpt[Lung]*Q[Lung] - PER*S[Lung]*Ci_org[Lung] + PER*S[Lung]*Ce_org[Lung]/(Kp[Lung]*Kpp) ) / Vi_org[Lung]

    # cytoplasmic subcompartment
    @. du.Cet_org = ( PER*S.*Ci_org .- PER*S.*Ce_org./(Kp*Kpp) ) ./ Ve_org

    # DNA-bound data 
    # @. CDNA_bound = Cet_org - Ce_org

    # blood compartment 
    du.Ca = ( Ci_org[Lung]/Kpt[Lung]*Q[Lung] - Ca*sum(Q[Liver:Tumor]) - PER*S_blood*Ca ) / Va
    du.Cp = ( sum(Ci_org[Kidney:Tumor]./Kpt[Kidney:Tumor].*Q[Kidney:Tumor]) + Ci_org[Liver]/Kpt[Liver]*(Q[Liver]+Q[Gut]+Q[Spleen]) - Cp*Q[Lung] - Cp*CLrenal + PER*S_blood*Cet_blood/Kpt_blood ) / Vp + infusion
    du.Ce_blood = 1/Ve_blood * ( PER*S_blood*Ca - PER*S_blood*Cet_blood/Kpt_blood )  

    # dummy variable to track clearance 
    du.totalclearance = Cp*CLrenal + CLhepatic*Ci_org[Liver]

end
