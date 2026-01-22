# date: 1/5/2026 
# author: Yuezhe Li 
# purpose of this code: to set up initial condition of pbpk model 

# all other parameters/initial values are taken from this publication:
# Jones et al., CPT: Pharmacometrics & Systems Pharmacology, 2019
# link:  https://doi.org/10.1002/psp4.12461
FcRn_Conc_homo = 44.093; # [uM]  
EDG_mg_ml = 10.0;  # [mg/mL]
# C_EXG_Plasma0 =  Dose_in_mgkg*BW_homo/(V_Plasma*1000)/MW_EDG*1E6; # [uM]
# C_EDG_Plasma0 =  EDG_mg_ml/MW_EDG*1E6; # [uM]
# C_EDG_LN0 = 1E-18*1E6; # [uM]
# C_EDG0 =  1E-18*1E6 * ones(15,17); # [uM]
# C_EXG_LN0 = 1E-18*1E6; # [uM]
# C_EXG0 =  1E-18*1E6 *ones(15,19); # [uM]
# C_FcRn_E6a0 =  ones(15)*FcRn_Conc; # [uM]
# C_FcRn_E7b0 =  ones(15)*FcRn_Conc; # [uM] 
# C_FcRn_E70 =  ones(15)*FcRn_Conc; # [uM]      
# C_FcRn_ISM0 =  ones(15)*FcRn_Conc*1E-4; # [uM] 
# C_FcRn_VM0 =  ones(15)*FcRn_Conc*1E-4; # [uM] 

function pbpk_initial_condition(Dose_in_mgkg, mdl; 
    EDG_mg_ml = EDG_mg_ml, V_Plasma = 3.126, MW_EDG = MW_IGG, BW = human_WT, FcRn_Conc = FcRn_Conc_homo)
    ## get a u0 of zeros
    u0_map = Dict([var => 0.0 for var in unknowns(mdl)])

    ### plasma (uM) ; dose in plasma_exg
    u0_map[mdl.plasma_exg.C_Plasma.val] = Dose_in_mgkg*BW/(V_Plasma*1000)/MW_EDG*1e6; # [uM]
    u0_map[mdl.plasma_edg.C_Plasma.val] = EDG_mg_ml/MW_EDG*1e6; # [uM]

    ### LN (uM)
    u0_map[mdl.ln_edg.C_LN.val] = 1e-18*1e6; # [uM]

    ### Plasma PL 
    u0_map[mdl.plasma_pl.C_PL_Plasma.val] = 0; # [uM]

    ### all other IgG species (uM)
    subcompartments = [:C_V, 
                    :C_VM, :C_bound_VM, :C_bound2_VM, :C_bound_VM_mem, 
                    :C_E7, :C_bound_E7, :C_bound2_E7, 
                    :C_E6a, :C_bound_E6a, :C_bound2_E6a,
                    :C_E7b, :C_bound_E7b, :C_bound2_E7b,
                    :C_ISM, :C_bound_ISM, :C_bound2_ISM, :C_bound_ISM_mem,
                    :C_IntS]
    # subsystems = nameof.(ModelingToolkit.get_systems(pbpk);)
    for uk in unknowns(mdl)
        if occursin("edg", string(uk)) || occursin("exg", string(uk))
            for sc in subcompartments
                if occursin(String(sc), string(uk))
                    # println(uk)
                    # println(sc)
                    u0_map[uk] = 1e-18*1e6 # [uM]    
                end
            end
        end
    end


    ### FcRn species (uM)
    for uk in unknowns(mdl)
            if occursin("C_FcRn",string(uk), )
                if occursin("VM", string(uk),) || occursin("ISM",string(uk), )
                    u0_map[uk] = FcRn_Conc*1e-4 # [uM]
                else
                    u0_map[uk] = FcRn_Conc # [uM]
                end
            end
    end

    return u0_map
end

