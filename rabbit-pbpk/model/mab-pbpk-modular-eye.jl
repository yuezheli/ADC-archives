## updated for MTK v10.0.1.0
## this script has parameters in the subcompartment model functions ;; so subcomarrment functions can be compiled
## based on mab-pbpk2.jl but define derived params in vars block


Symbolics.option_to_metadata_type(::Val{:symscope}) = SymScope



#  ==============================================================================================
#  ==============================================================================================

#################################################
########### Sub-compartment functions ###########
#################################################

#---------------------------------------#
#--------------# IgG #------------------#
#---------------------------------------#

# IgG (endogenous and exogenous; free and bound)

#  Organ: All 
#  Cmpt:  vascular side membrane, endosome(7.4, 6a, 6b), interstitial side membrane, interstitial space
#  Species: exogenous IgG (free & bound), endogeneous (free & bound) IgG and FcRn  

# vascular space ; all except liver and lung
function V_ALL(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        PLQ, [symscope = ParentScope(LocalScope())]
        V_V, [symscope = ParentScope(LocalScope())]
        sigma_V, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        FR, [symscope = GlobalScope()]
        LF, [symscope = ParentScope(LocalScope())]
    end

    vars = @variables begin
        C_V(t)
        C_VM(t)
        C_V_Lung(t)
    end

    # derived parameters
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15

    CL_up = CL_up * CLOff
    eqs = [
      Dt(C_V) ~ (
        PLQ * C_V_Lung  #  from Lung 
        - (PLQ - LF) * C_V  #  leave to Main Plasma
        - (1.0 - sigma_V) * LF * C_V #  going via Lymph to Interstitial 
        - CL_up * FR * C_V #  pinocytosis from vascular to vascular membrane
        + CL_up * FR * C_VM #  exocytosis from vascular membrane to vascular space 
      ) / V_V
    ]

    System(eqs, t, vars, pars; name=name)
end


# vascular space - liver
function V_LIVER(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        PLQ, [symscope = ParentScope(LocalScope())]
        PLQ_Spleen, [symscope = GlobalScope()]
        PLQ_Pancreas,   [symscope = GlobalScope()]
        PLQ_SI, [symscope = GlobalScope()]
        PLQ_LI, [symscope = GlobalScope()]
        LF,   [symscope = ParentScope(LocalScope())]
        LF_Spleen, [symscope = GlobalScope()]
        LF_Pancreas, [symscope = GlobalScope()]
        LF_SI,  [symscope = GlobalScope()]
        LF_LI, [symscope = GlobalScope()]
        sigma_V, [symscope = ParentScope(LocalScope())]
        V_V, [symscope = ParentScope(LocalScope())]
        FR, [symscope = GlobalScope()]
    end
    
    vars = @variables begin
        C_V(t)
        C_VM(t)
        C_V_Lung(t)
        C_V_Spleen(t)
        C_V_Pancreas(t)
        C_V_SI(t)
        C_V_LI(t)
    end

    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    CL_up = CL_up * CLOff
    eqs = [
        Dt(C_V) ~ (
            PLQ * C_V_Lung #  Inlet distributed from Lung
            + (PLQ_Spleen - LF_Spleen) * C_V_Spleen #  Inlet from Spleen
            + (PLQ_Pancreas - LF_Pancreas) * C_V_Pancreas #  Inlet from Pancreas
            + (PLQ_SI - LF_SI) * C_V_SI  #  Inlet from S.I
            + (PLQ_LI - LF_LI) * C_V_LI   #  Inlet from L.I
            - C_V * (PLQ - LF + PLQ_Spleen - LF_Spleen + PLQ_Pancreas - LF_Pancreas + PLQ_SI - LF_SI + PLQ_LI - LF_LI) #  Outlet
            - (1.0 - sigma_V) * LF * C_V
            - CL_up * FR * C_V #  pinocytosis from vascular to vascular membrane
            + CL_up * FR * C_VM #  exocytosis from vascular membrane to vascular space 
        ) / V_V
    ]

    System(eqs, t, vars, pars; name=name)
end

# vascular space - lung
function V_LUNG(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        V_V, [symscope = ParentScope(LocalScope())]
        LF, [symscope = ParentScope(LocalScope())]
        PLQ, [symscope = ParentScope(LocalScope())]
        sigma_V, [symscope = ParentScope(LocalScope())]
        FR, [symscope = GlobalScope()]
    end
    vars = @variables begin
        C_V(t)
        C_VM(t)
        C_Plasma(t)
    end

    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    CL_up = CL_up * CLOff
    eqs = [
    #  Lung, vascular cmpt, exogenous free only
        Dt(C_V) ~ (
            (PLQ + LF) * C_Plasma
               - PLQ * C_V
                - (1.0 - sigma_V) * LF * C_V
                - CL_up * FR * C_V #  pinocytosis from vascular to vascular membrane
                + CL_up * FR * C_VM #  exocytosis from vascular membrane to vascular space
        ) / V_V
    ]

    System(eqs, t, vars, pars; name=name)
end

################
################

#  Organ: All 
#  Cmpt:  vascular side membrane, endosome(7.4, 6a, 6b), interstitial side membrane, interstitial space
#  Species: exogenous IgG (free & bound), endogeneous (free & bound) IgG and FcRn  

# vascular side membrane - free
function VM(; name, exg)
    exg = Float64.(exg) # ensure exg is a Float64 for calculations
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        tau_VM , [symscope = GlobalScope()]
        k_on_7_EXG, [symscope = GlobalScope()]
        k_on_7_EDG, [symscope = GlobalScope()]
        KD_7_EXG, [symscope = GlobalScope()]
        KD_7_EDG, [symscope = GlobalScope()]
        k_on_PS, [symscope = GlobalScope()]
        k_off_PS, [symscope = GlobalScope()]
        FR, [symscope = GlobalScope()]
        C_Mem, [symscope = GlobalScope()]
    end
            vars = @variables begin
        C_VM(t)
        C_V(t)
        C_E7b(t)
        C_FcRn_VM(t), [symscope = ParentScope(LocalScope())]
        C_bound_VM(t)
        C_bound_VM_mem(t)
    end

    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_VM = CL_up * tau_VM
    k_on_7 = (exg * k_on_7_EXG) + ((1.0 - exg) * k_on_7_EDG)
    KD_7 = (exg * KD_7_EXG) + ((1.0 - exg) * KD_7_EDG)
    k_off_7 = (KD_7 / 1000.0) * k_on_7 # [1.0/h] 

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
        Dt(C_VM) ~ (
            CL_up * FR * C_V                   #  from vascular space
            - CL_up * FR * C_VM                 #  to vascular space
            - CL_up * FR * C_VM                 #  to endosomal pH=7.4
            + CL_up * FR * C_E7b                #  from endosomal pH=6
            - k_on_7 * C_VM * C_FcRn_VM * V_VM  #  rxn on 
            + k_off_7 * C_bound_VM * V_VM       #  rxn off
            - exg * k_on_PS * C_VM * C_Mem * V_VM     #  reaction of membrane ; adding exg flag to switch between EXG and EDG
            + exg * k_off_PS * C_bound_VM_mem * V_VM
        ) / V_VM
    ]


    System(eqs, t, vars, pars; name=name)
end


#---------#

# vascular side membrane - bound
function BOUND_VM_MEM(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space

        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        tau_VM, [symscope = GlobalScope()]
        k_off_PS, [symscope = GlobalScope()]
        kint_PS, [symscope = GlobalScope()]
        k_on_PS, [symscope = GlobalScope()]
        C_Mem, [symscope = GlobalScope()]
    end
        

    vars = @variables begin
        C_bound_VM_mem(t)
        C_VM(t)
    end

    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_VM = CL_up * tau_VM 
    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
        Dt(C_bound_VM_mem) ~ (
            k_on_PS * C_VM * C_Mem * V_VM #  reaction of membrane 
            - k_off_PS * C_bound_VM_mem * V_VM
            - kint_PS * C_bound_VM_mem * V_VM
        ) / V_VM
    ]

    System(eqs, t, vars, pars; name=name)
end

#---------#

# endosomal pH = 7.4
function E7(; name, exg)
    exg = Float64.(exg)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac,      [symscope = ParentScope(LocalScope())]
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        pino_time, [symscope = GlobalScope()]
        tau_VM, [symscope = GlobalScope()]
        E6a_Vol_Pct, [symscope = GlobalScope()]
        tau_ISM, [symscope = GlobalScope()]
        k_on_7_EXG, [symscope = GlobalScope()]
        k_on_7_EDG, [symscope = GlobalScope()]
        KD_7_EXG, [symscope = GlobalScope()]
        KD_7_EDG, [symscope = GlobalScope()]
        kint_PS, [symscope = GlobalScope()]
        FR,     [symscope = GlobalScope()]
    end
            
    vars = @variables begin
        C_E7(t)
        C_VM(t)
        C_ISM(t)
        C_bound_VM_mem(t)
        C_bound_ISM_mem(t)
        C_FcRn_E7(t), [symscope = ParentScope(LocalScope())]
        C_bound_E7(t)
        
    end

    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_endosomal = CL_up_in_nL_per_hour_per_million_cells * pino_time * Organ_Endothelial_Cell * 1e-6 * 1e-9
    V_VM = CL_up * tau_VM
    E7_Vol_Pct = (1.0 - E6a_Vol_Pct) / 2 # [-]
    V_E7 = V_endosomal * E7_Vol_Pct
    V_ISM = CL_up * tau_ISM
    k_on_7 = (exg * k_on_7_EXG) + ((1.0 - exg) * k_on_7_EDG)
    KD_7 = (exg * KD_7_EXG) + ((1.0 - exg) * KD_7_EDG)
    k_off_7 = (KD_7 / 1000.0) * k_on_7 # [1.0/h] 

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
        Dt(C_E7) ~ (
            CL_up * FR * C_VM                              #  from vascular membrane
            + CL_up * (1.0 - FR) * C_ISM                        #  from is membrane
            - CL_up * C_E7                                  #  to endosomal pH=6
            - k_on_7 * C_E7 * C_FcRn_E7 * V_E7 #  rxn on 
            + k_off_7 * C_bound_E7 * V_E7               #  rxn off
            + exg * kint_PS * C_bound_VM_mem * V_VM               #  from membrane bound due to polyspecificity
            + exg * kint_PS * C_bound_ISM_mem * V_ISM             #  from membrane bound due to polyspecificity
        ) / V_E7
    ]

    System(eqs, t, vars, pars; name=name)
end


#---------#

# E6a
function E6A(; name, exg)
    exg = Float64.(exg)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Scale_Factor, [symscope = GlobalScope()]
        pino_time,  [symscope = GlobalScope()]
        E6a_Vol_Pct, [symscope = GlobalScope()]
        k_on_6_EXG, [symscope = GlobalScope()]
        k_on_6_EDG, [symscope = GlobalScope()]
        KD_6_EXG, [symscope = GlobalScope()]
        KD_6_EDG, [symscope = GlobalScope()]
    end
    
        

    vars = @variables begin
        C_E6a(t)
        C_E7(t)
        C_bound_E6a(t)
        C_FcRn_E6a(t), [symscope = ParentScope(LocalScope())]
    end

    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_endosomal = CL_up_in_nL_per_hour_per_million_cells * pino_time * Organ_Endothelial_Cell * 1e-6 * 1e-9
    V_E6a = V_endosomal * E6a_Vol_Pct
    k_on_6 = (exg * k_on_6_EXG) + ((1.0 - exg) * k_on_6_EDG)
    KD_6 = (exg * KD_6_EXG) + ((1.0 - exg) * KD_6_EDG)
    k_off_6 = (KD_6 / 1000.0) * k_on_6 # [1.0/h]

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
       #  Free EXG mAb in E6a
       Dt(C_E6a) ~ (
        CL_up * C_E7 #  from E7
            - CL_up * C_E6a #  move to E7b or being routed to lysosomal for degradation
            - k_on_6 * C_E6a * C_FcRn_E6a * V_E6a #  rxn on 
            + k_off_6 * C_bound_E6a * V_E6a  #  rxn off
        ) / V_E6a
    ]

    System(eqs, t, vars, pars; name=name)
end

#---------#

# endosomal pH = 7.4 b
function E7B(; name, exg)
    exg = Float64.(exg)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        E6a_Vol_Pct, [symscope = GlobalScope()]
        pino_time, [symscope = GlobalScope()]
        Prob_deg, [symscope = GlobalScope()]
        k_on_7_EXG, [symscope = GlobalScope()]
        k_on_7_EDG, [symscope = GlobalScope()]
        KD_7_EXG, [symscope = GlobalScope()]
        KD_7_EDG, [symscope = GlobalScope()]
    end
#    Organ Level Parameters

        vars = @variables begin
        C_E7b(t)
        C_E6a(t)
        C_FcRn_E7b(t), [symscope = ParentScope(LocalScope())]
        C_bound_E7b(t)
    end

    # derived params
    E7b_Vol_Pct = (1.0 - E6a_Vol_Pct) / 2 # [-]
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_endosomal = CL_up_in_nL_per_hour_per_million_cells * pino_time * Organ_Endothelial_Cell * 1e-6 * 1e-9
    V_E7b = V_endosomal * E7b_Vol_Pct
    k_on_7 = (exg * k_on_7_EXG) + ((1.0 - exg) * k_on_7_EDG)
    KD_7 = (exg * KD_7_EXG) + ((1.0 - exg) * KD_7_EDG)
    k_off_7 = (KD_7 / 1000.0) * k_on_7 # [1.0/h]

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [Dt(C_E7b) ~ (
            CL_up * (1.0 - Prob_deg) * C_E6a         #  from 6a
            - CL_up * C_E7b                        #  leave from E7b to membranes 
            - k_on_7 * C_E7b * C_FcRn_E7b * V_E7b  #  rxn on 
            + k_off_7 * C_bound_E7b * V_E7b        #  rxn off
        ) / V_E7b
    ]

    System(eqs, t, vars, pars; name=name)
end

#---------#

# IS side mem - free
function ISM(; name, exg)
    exg = Float64.(exg)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        tau_ISM, [symscope = GlobalScope()]
        k_on_7_EXG, [symscope = GlobalScope()]
        k_on_7_EDG, [symscope = GlobalScope()]
        KD_7_EXG, [symscope = GlobalScope()]
        KD_7_EDG, [symscope = GlobalScope()]
        k_on_PS, [symscope = GlobalScope()]
        k_off_PS, [symscope = GlobalScope()]
        FR, [symscope = GlobalScope()]
        C_Mem, [symscope = GlobalScope()]
    end
    vars = @variables begin
        C_ISM(t)
        C_IntS(t)
        C_E7b(t)
        C_bound_ISM(t)
        C_bound_ISM_mem(t)
        C_FcRn_ISM(t), [symscope = ParentScope(LocalScope())]
    end


    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_ISM = CL_up * tau_ISM
    k_on_7 = (exg * k_on_7_EXG) + ((1.0 - exg) * k_on_7_EDG)
    KD_7 = (exg * KD_7_EXG) + ((1.0 - exg) * KD_7_EDG)
    k_off_7 = (KD_7 / 1000.0) * k_on_7 # [1.0/h]

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [Dt(C_ISM) ~ (
        CL_up * (1.0 - FR) * C_IntS                  #  from IS space 
            - CL_up * (1.0 - FR) * C_ISM                 #  to endosomal pH=7.4
            + CL_up * (1.0 - FR) * C_E7b                 #  from endosomal pH=7.4 b
            - CL_up * (1.0 - FR) * C_ISM                 #  to interstitial space
            - k_on_7 * C_ISM * C_FcRn_ISM * V_ISM  #  rxn on 
            + k_off_7 * C_bound_ISM * V_ISM        #  rxn off
            - exg * k_on_PS * C_ISM * C_Mem * V_ISM          #  reaction of membrane 
            + exg * k_off_PS * C_bound_ISM_mem * V_ISM
        ) / V_ISM
    ]


    System(eqs, t, vars, pars; name=name)
end

#-----------#

# IS side mem - bound
function BOUND_ISM_MEM(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        tau_ISM, [symscope = GlobalScope()]
        k_on_PS, [symscope = GlobalScope()]
        k_off_PS, [symscope = GlobalScope()]
        kint_PS,    [symscope = GlobalScope()]
        C_Mem, [symscope = GlobalScope()]
    end
    vars = @variables begin
        C_bound_ISM_mem(t)
        C_ISM(t)
    end

    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_ISM = CL_up * tau_ISM

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
        Dt(C_bound_ISM_mem) ~ (
        k_on_PS * C_ISM * C_Mem * V_ISM #  reaction of membrane 
            - k_off_PS * C_bound_ISM_mem * V_ISM
            - kint_PS * C_bound_ISM_mem * V_ISM
        ) / V_ISM
    ]

    System(eqs, t, vars, pars; name=name)
end

#-----------#

#  Interstitial space
function INTS(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        sigma_IS, [symscope = ParentScope(LocalScope())]
        PLQ, [symscope = ParentScope(LocalScope())]
        V_IntS, [symscope = ParentScope(LocalScope())]
        sigma_V,    [symscope = ParentScope(LocalScope())]
        FR, [symscope = GlobalScope()]
        LF, [symscope = ParentScope(LocalScope())]
    end
        
    vars = @variables begin
        C_IntS(t)
        C_V(t)
        C_ISM(t)
    end

    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
            Dt(C_IntS) ~ (
            (1.0 - sigma_V) * LF * C_V        #  going from vascular via Lymph 
            - CL_up * (1.0 - FR) * C_IntS     #  pinocytosis 
            + CL_up * (1.0 - FR) * C_ISM      #  exocytosis
            - (1.0 - sigma_IS) * LF * C_IntS
        ) / V_IntS
    ]

    System(eqs, t, vars, pars; name=name)
end

#-----------#
#-----------#

# FcRn-IgG #

##  bounded 
function BOUND_VM(; name, exg)
    exg = Float64.(exg)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        tau_VM, [symscope = GlobalScope()]
        k_on_7_EXG, [symscope = GlobalScope()]
        k_on_7_EDG, [symscope = GlobalScope()]
        KD_7_EXG, [symscope = GlobalScope()]
        KD_7_EDG, [symscope = GlobalScope()]
        kdeg_FcRn_Ab,  [symscope = GlobalScope()]
        on_rate_ratio_1st_vs_2nd_binding, [symscope = GlobalScope()]
        FR, [symscope = GlobalScope()]
    end
        
    vars = @variables begin
        C_bound_VM(t)
        C_VM(t)
        C_FcRn_VM(t), [symscope = ParentScope(LocalScope())]
        C_bound2_VM(t)
        C_bound_E7b(t)
    end

    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_VM = CL_up * tau_VM
    k_on_7 = (exg * k_on_7_EXG) + ((1.0 - exg) * k_on_7_EDG)
    KD_7 = (exg * KD_7_EXG) + ((1.0 - exg) * KD_7_EDG)
    k_off_7 = (KD_7 / 1000.0) * k_on_7 # [1.0/h]
    k_on_7_2 = k_on_7 / on_rate_ratio_1st_vs_2nd_binding # [1.0/(uM*h)]
    KD_7_2 = KD_7 * on_rate_ratio_1st_vs_2nd_binding # [nM] 
    k_off_7_2 = (KD_7_2 / 1000.0) * k_on_7_2 # [1.0/h]

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
            Dt(C_bound_VM) ~ (
            k_on_7 * C_VM * C_FcRn_VM * V_VM              #  rxn on 
            - k_off_7 * C_bound_VM * V_VM                 #  rxn off
            - k_on_7_2 * C_bound_VM * C_FcRn_VM * V_VM    #  rxn on 
            + k_off_7_2 * C_bound2_VM * V_VM              #  rxn off 
            + CL_up * FR * C_bound_E7b                    #  from Endosomal 6b to Vascular membrane
            - CL_up * FR * C_bound_VM                     #  to Endosomal 7
            - kdeg_FcRn_Ab * C_bound_VM * V_VM            #  recover FcRn when IgG-FrRn is destroyed
       ) / V_VM
    ]

    System(eqs, t, vars, pars; name=name)
end

#-----------#

#  bounded2
function BOUND2_VM(; name, exg)
    exg = Float64.(exg)
    @independent_variables t
    Dt = Differential(t)
    
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        tau_VM, [symscope = GlobalScope()]
        k_on_7_EXG, [symscope = GlobalScope()]
        k_on_7_EDG, [symscope = GlobalScope()]
        KD_7_EXG, [symscope = GlobalScope()]
        KD_7_EDG, [symscope = GlobalScope()]
        kdeg_FcRn_Ab, [symscope = GlobalScope()]
        on_rate_ratio_1st_vs_2nd_binding, [symscope = GlobalScope()]
        FR, [symscope = GlobalScope()]
    end
    vars = @variables begin
        C_bound2_VM(t)
        C_bound_VM(t)
        C_FcRn_VM(t), [symscope = ParentScope(LocalScope())]
        C_bound2_E7b(t)
    end

    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_VM = CL_up * tau_VM
    k_on_7 = (exg * k_on_7_EXG) + ((1.0 - exg) * k_on_7_EDG)
    KD_7 = (exg * KD_7_EXG) + ((1.0 - exg) * KD_7_EDG)
    k_on_7_2 = k_on_7 / on_rate_ratio_1st_vs_2nd_binding # [1.0/(uM*h)]
    KD_7_2 = KD_7 * on_rate_ratio_1st_vs_2nd_binding # [nM] 
    k_off_7_2 = (KD_7_2 / 1000.0) * k_on_7_2 # [1.0/h]

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
        Dt(C_bound2_VM) ~ (
            k_on_7_2 * C_bound_VM * C_FcRn_VM * V_VM #  rxn on 
            - k_off_7_2 * C_bound2_VM * V_VM #  rxn off 
            + CL_up * FR * C_bound2_E7b  #  from Endosomal 6b to Vascular membrane
            - CL_up * FR * C_bound2_VM #  to Endosomal 7
            - kdeg_FcRn_Ab * C_bound2_VM * V_VM #  recover FcRn when IgG-FrRn is destroyed
       ) / V_VM
    ]

    System(eqs, t, vars, pars; name=name)
end

#------------#

#  endosomal space 
## bound E7
function BOUND_E7(; name, exg)
    exg = Float64.(exg)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        pino_time, [symscope = GlobalScope()]
        E6a_Vol_Pct,    [symscope = GlobalScope()]
        k_on_7_EXG, [symscope = GlobalScope()]
        k_on_7_EDG, [symscope = GlobalScope()]
        KD_7_EXG, [symscope = GlobalScope()]
        KD_7_EDG, [symscope = GlobalScope()]
        on_rate_ratio_1st_vs_2nd_binding, [symscope = GlobalScope()]
        FR, [symscope = GlobalScope()]
    end
    vars = @variables begin
        C_bound_E7(t)
        C_E7(t)
        C_bound_VM(t)
        C_FcRn_E7(t), [symscope = ParentScope(LocalScope())]
        C_bound2_E7(t)
        C_bound_ISM(t)
    end


    # derived params
    E7_Vol_Pct = (1.0 - E6a_Vol_Pct) / 2 # [-]
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_endosomal = CL_up_in_nL_per_hour_per_million_cells * pino_time * Organ_Endothelial_Cell * 1e-6 * 1e-9
    V_E7 = V_endosomal * E7_Vol_Pct
    k_on_7 = (exg * k_on_7_EXG) + ((1.0 - exg) * k_on_7_EDG)
    KD_7 = (exg * KD_7_EXG) + ((1.0 - exg) * KD_7_EDG)
    k_off_7 = (KD_7 / 1000.0) * k_on_7 # [1.0/h]
    k_on_7_2 = k_on_7 / on_rate_ratio_1st_vs_2nd_binding # [1.0/(uM*h)]
    KD_7_2 = KD_7 * on_rate_ratio_1st_vs_2nd_binding # [nM] 
    k_off_7_2 = (KD_7_2 / 1000.0) * k_on_7_2 # [1.0/h]

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
        Dt(C_bound_E7) ~ (
            k_on_7 * C_E7 * C_FcRn_E7 * V_E7 #  rxn on 
            - k_off_7 * C_bound_E7 * V_E7 #  rxn off
            - k_on_7_2 * C_bound_E7 * C_FcRn_E7 * V_E7 #  rxn on 
            + k_off_7_2 * C_bound2_E7 * V_E7 #  rxn off 
            - CL_up * C_bound_E7 #  goes to E6
            + CL_up * FR * C_bound_VM #  from VM 
            + CL_up * (1.0 - FR) * C_bound_ISM  #  from ISM
         ) / V_E7
    ]

    System(eqs, t, vars, pars; name=name)
end

#-----------#

## bound2 E7
function BOUND2_E7(; name, exg)
    exg = Float64.(exg)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        pino_time, [symscope = GlobalScope()]
        E6a_Vol_Pct, [symscope = GlobalScope()]
        k_on_7_EXG, [symscope = GlobalScope()]
        k_on_7_EDG, [symscope = GlobalScope()]
        KD_7_EXG, [symscope = GlobalScope()]
        KD_7_EDG, [symscope = GlobalScope()]
        on_rate_ratio_1st_vs_2nd_binding,   [symscope = GlobalScope()]
        FR, [symscope = GlobalScope()]
    end
    vars = @variables begin
        C_bound2_E7(t)
        C_bound_E7(t)
        C_bound2_VM(t)
        C_bound2_ISM(t)
        C_FcRn_E7(t), [symscope = ParentScope(LocalScope())]
    end

    # derived params
    E7_Vol_Pct = (1.0 - E6a_Vol_Pct) / 2 # [-]
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_endosomal = CL_up_in_nL_per_hour_per_million_cells * pino_time * Organ_Endothelial_Cell * 1e-6 * 1e-9
    V_E7 = V_endosomal * E7_Vol_Pct
    k_on_7 = (exg * k_on_7_EXG) + ((1.0 - exg) * k_on_7_EDG)
    KD_7 = (exg * KD_7_EXG) + ((1.0 - exg) * KD_7_EDG)
    k_on_7_2 = k_on_7 / on_rate_ratio_1st_vs_2nd_binding # [1.0/(uM*h)]
    KD_7_2 = KD_7 * on_rate_ratio_1st_vs_2nd_binding # [nM] 
    k_off_7_2 = (KD_7_2 / 1000.0) * k_on_7_2 # [1.0/h]

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
        Dt(C_bound2_E7) ~ (
            k_on_7_2 * C_bound_E7 * C_FcRn_E7 * V_E7 #  rxn on 
            - k_off_7_2 * C_bound2_E7 * V_E7 #  rxn off 
            - CL_up * C_bound2_E7 #  goes to E6
            + CL_up * FR * C_bound2_VM #  from VM 
            + CL_up * (1.0 - FR) * C_bound2_ISM   #  from ISM
       ) / V_E7
    ]

    System(eqs, t, vars, pars; name=name)
end

#-----------#

## endosomal pH=6 a
### bound
function BOUND_E6A(; name, exg)
    exg = Float64.(exg)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())] 
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor,   [symscope = GlobalScope()]
        pino_time, [symscope = GlobalScope()]
        E6a_Vol_Pct, [symscope = GlobalScope()]
        k_on_6_EXG, [symscope = GlobalScope()]
        k_on_6_EDG, [symscope = GlobalScope()]
        KD_6_EXG, [symscope = GlobalScope()]
        KD_6_EDG, [symscope = GlobalScope()]
        on_rate_ratio_1st_vs_2nd_binding,   [symscope = GlobalScope()]
    end
            
    vars = @variables begin
        C_bound_E6a(t)
        C_E6a(t)
        C_bound2_E6a(t)
        C_FcRn_E6a(t), [symscope = ParentScope(LocalScope())]
        C_bound_E7(t)
    end

    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_endosomal = CL_up_in_nL_per_hour_per_million_cells * pino_time * Organ_Endothelial_Cell * 1e-6 * 1e-9
    V_E6a = V_endosomal * E6a_Vol_Pct
    k_on_6 = (exg * k_on_6_EXG) + ((1.0 - exg) * k_on_6_EDG)
    KD_6 = (exg * KD_6_EXG) + ((1.0 - exg) * KD_6_EDG)      
    k_off_6 = (KD_6 / 1000.0) * k_on_6 # [1.0/h]
    k_on_6_2 = k_on_6 / on_rate_ratio_1st_vs_2nd_binding # [1.0/(uM*h)]
    KD_6_2 = KD_6 * on_rate_ratio_1st_vs_2nd_binding # [nM] 
    k_off_6_2 = (KD_6_2 / 1000.0) * k_on_6_2 # [1.0/h]

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
        Dt(C_bound_E6a) ~ (
            k_on_6 * C_E6a * C_FcRn_E6a * V_E6a #  rxn on 
            - k_off_6 * C_bound_E6a * V_E6a #  rxn off
            - k_on_6_2 * C_bound_E6a * C_FcRn_E6a * V_E6a #  rxn on 
            + k_off_6_2 * C_bound2_E6a * V_E6a #  rxn off 
            + CL_up * C_bound_E7 #  internalization, from pH7 to pH6a
            - CL_up * C_bound_E6a #  exocytosis, from E6a to E7b
        ) / V_E6a
    ]

    System(eqs, t, vars, pars; name=name)
end


#-----------#

## endosomal pH=6 a
### bound2
function BOUND2_E6A(; name, exg)
    exg = Float64.(exg)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell,     [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        pino_time, [symscope = GlobalScope()]
        E6a_Vol_Pct, [symscope = GlobalScope()]
        k_on_6_EXG, [symscope = GlobalScope()]
        k_on_6_EDG, [symscope = GlobalScope()]
        KD_6_EXG, [symscope = GlobalScope()]
        KD_6_EDG, [symscope = GlobalScope()]
        on_rate_ratio_1st_vs_2nd_binding, [symscope = GlobalScope()]
    end
        
    vars = @variables begin
        C_bound2_E6a(t)
        C_bound_E6a(t)
        C_FcRn_E6a(t), [symscope = ParentScope(LocalScope())]
        C_bound2_E7(t)
    end

    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_endosomal = CL_up_in_nL_per_hour_per_million_cells * pino_time * Organ_Endothelial_Cell * 1e-6 * 1e-9
    V_E6a = V_endosomal * E6a_Vol_Pct
    k_on_6 = (exg * k_on_6_EXG) + ((1.0 - exg) * k_on_6_EDG)
    KD_6 = (exg * KD_6_EXG) + ((1.0 - exg) * KD_6_EDG)
    k_on_6_2 = k_on_6 / on_rate_ratio_1st_vs_2nd_binding # [1.0/(uM*h)]
    KD_6_2 = KD_6 * on_rate_ratio_1st_vs_2nd_binding # [nM] 
    k_off_6_2 = (KD_6_2 / 1000.0) * k_on_6_2 # [1.0/h]

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
        Dt(C_bound2_E6a) ~ (
            k_on_6_2 * C_bound_E6a * C_FcRn_E6a * V_E6a #  rxn on 
            - k_off_6_2 * C_bound2_E6a * V_E6a #  rxn off 
            + CL_up * C_bound2_E7 #  internalization, from pH7 to pH6a
            - CL_up * C_bound2_E6a #  exocytosis, from E6a to E7b
        ) / V_E6a
    ]

    System(eqs, t, vars, pars; name=name)
end

#-----------#

## endosomal pH=7.4 b
### bound
function BOUND_E7B(; name, exg)
    exg = Float64.(exg)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        pino_time, [symscope = GlobalScope()]
        E6a_Vol_Pct, [symscope = GlobalScope()]
        k_on_7_EXG, [symscope = GlobalScope()]
        k_on_7_EDG, [symscope = GlobalScope()]
        KD_7_EXG, [symscope = GlobalScope()]
        KD_7_EDG, [symscope = GlobalScope()]
        on_rate_ratio_1st_vs_2nd_binding, [symscope = GlobalScope()]
    end
            vars = @variables begin
        C_bound_E7b(t)
        C_E7b(t)
        C_FcRn_E7b(t), [symscope = ParentScope(LocalScope())]
        C_bound2_E7b(t)
        C_bound_E6a(t)
    end

    # derived params
    E7b_Vol_Pct = (1.0 - E6a_Vol_Pct) / 2 # [-]

    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_endosomal = CL_up_in_nL_per_hour_per_million_cells * pino_time * Organ_Endothelial_Cell * 1e-6 * 1e-9
    V_E7b = V_endosomal * E7b_Vol_Pct
    k_on_7 = (exg * k_on_7_EXG) + ((1.0 - exg) * k_on_7_EDG)
    KD_7 = (exg * KD_7_EXG) + ((1.0 - exg) * KD_7_EDG)     
    k_off_7 = (KD_7 / 1000.0) * k_on_7 # [1.0/h]
    k_on_7_2 = k_on_7 / on_rate_ratio_1st_vs_2nd_binding # [1.0/(uM*h)]
    KD_7_2 = KD_7 * on_rate_ratio_1st_vs_2nd_binding # [nM] 
    k_off_7_2 = (KD_7_2 / 1000.0) * k_on_7_2 # [1.0/h]

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
        Dt(C_bound_E7b) ~ (
            k_on_7 * C_E7b * C_FcRn_E7b * V_E7b #  rxn on 
            - k_off_7 * C_bound_E7b * V_E7b #  rxn off
            - k_on_7_2 * C_bound_E7b * C_FcRn_E7b * V_E7b #  rxn on 
            + k_off_7_2 * C_bound2_E7b * V_E7b #  rxn off 
            + CL_up * C_bound_E6a #  internalization, from E6a to E7b
            - CL_up * C_bound_E7b #  exocytosis, from E7b to vascular/interstitial membrane
         ) / V_E7b
    ]

    System(eqs, t, vars, pars; name=name)
end

#-----------#

## endosomal pH=7.4 b
### bound2
function BOUND2_E7B(; name, exg)
    exg = Float64.(exg)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        pino_time, [symscope = GlobalScope()]
        E6a_Vol_Pct, [symscope = GlobalScope()]
        k_on_7_EXG, [symscope = GlobalScope()]
        k_on_7_EDG, [symscope = GlobalScope()]
        KD_7_EXG, [symscope = GlobalScope()]
        KD_7_EDG, [symscope = GlobalScope()]
        on_rate_ratio_1st_vs_2nd_binding,   [symscope = GlobalScope()]
    end
            vars = @variables begin
        C_bound2_E7b(t)
        C_bound_E7b(t)
        C_bound2_E6a(t)
        C_FcRn_E7b(t), [symscope = ParentScope(LocalScope())]
    end

    # derived params
    E7b_Vol_Pct = (1.0 - E6a_Vol_Pct) / 2 # [-]
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_endosomal = CL_up_in_nL_per_hour_per_million_cells * pino_time * Organ_Endothelial_Cell * 1e-6 * 1e-9
    V_E7b = V_endosomal * E7b_Vol_Pct
    k_on_7 = (exg * k_on_7_EXG) + ((1.0 - exg) * k_on_7_EDG)
    KD_7 = (exg * KD_7_EXG) + ((1.0 - exg) * KD_7_EDG)       
    k_on_7_2 = k_on_7 / on_rate_ratio_1st_vs_2nd_binding # [1.0/(uM*h)]
    KD_7_2 = KD_7 * on_rate_ratio_1st_vs_2nd_binding # [nM] 
    k_off_7_2 = (KD_7_2 / 1000.0) * k_on_7_2 # [1.0/h]

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
        Dt(C_bound2_E7b) ~ (
            k_on_7_2 * C_bound_E7b * C_FcRn_E7b * V_E7b #  rxn on 
            - k_off_7_2 * C_bound2_E7b * V_E7b #  rxn off 
            + CL_up * C_bound2_E6a #  internalization, from E6a to E7b
            - CL_up * C_bound2_E7b #  exocytosis, from E7b to vascular/interstitial membrane
         ) / V_E7b
    ]

    System(eqs, t, vars, pars; name=name)
end


#-----------#

## Interstitial membrane
### bound
function BOUND_ISM(; name, exg)
    exg = Float64.(exg)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        tau_ISM, [symscope = GlobalScope()]
        k_on_7_EXG, [symscope = GlobalScope()]
        k_on_7_EDG, [symscope = GlobalScope()]
        KD_7_EXG, [symscope = GlobalScope()]
        KD_7_EDG, [symscope = GlobalScope()]
        kdeg_FcRn_Ab, [symscope = GlobalScope()]
        on_rate_ratio_1st_vs_2nd_binding, [symscope = GlobalScope()]
        FR , [symscope = GlobalScope()]
    end
            vars = @variables begin
        C_bound_ISM(t)
        C_ISM(t)
        C_bound2_ISM(t)
        C_bound_E7b(t)
        C_FcRn_ISM(t), [symscope = ParentScope(LocalScope())]
    end

    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_ISM = CL_up * tau_ISM
    k_on_7 = (exg * k_on_7_EXG) + ((1.0 - exg) * k_on_7_EDG)
    KD_7 = (exg * KD_7_EXG) + ((1.0 - exg) * KD_7_EDG)       
    k_off_7 = (KD_7 / 1000.0) * k_on_7 # [1.0/h]
    k_on_7_2 = k_on_7 / on_rate_ratio_1st_vs_2nd_binding # [1.0/(uM*h)]
    KD_7_2 = KD_7 * on_rate_ratio_1st_vs_2nd_binding # [nM] 
    k_off_7_2 = (KD_7_2 / 1000.0) * k_on_7_2 # [1.0/h]

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
        Dt(C_bound_ISM) ~ (
            k_on_7 * C_ISM * C_FcRn_ISM * V_ISM #  rxn on  
            - k_off_7 * C_bound_ISM * V_ISM  #  rxn off
            - k_on_7_2 * C_bound_ISM * C_FcRn_ISM * V_ISM #  rxn on 
            + k_off_7_2 * C_bound2_ISM * V_ISM #  rxn off 
            + CL_up * (1.0 - FR) * C_bound_E7b  #  from Endosomal 6 to interstitial membrane
            - CL_up * (1.0 - FR) * C_bound_ISM #  to Endosomal 7
            - kdeg_FcRn_Ab * C_bound_ISM * V_ISM  #  recover FcRn when IgG-FrRn is destroyed
       ) / V_ISM
    ]

    System(eqs, t, vars, pars; name=name)
end

#-----------#

## Interstitial membrane
### bound2
function BOUND2_ISM(; name, exg)
    exg = Float64.(exg)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        tau_ISM, [symscope = GlobalScope()]
        k_on_7_EXG, [symscope = GlobalScope()]
        k_on_7_EDG, [symscope = GlobalScope()]
        KD_7_EXG,   [symscope = GlobalScope()]
        KD_7_EDG, [symscope = GlobalScope()]
        kdeg_FcRn_Ab, [symscope = GlobalScope()]
        on_rate_ratio_1st_vs_2nd_binding, [symscope = GlobalScope()]
        FR, [symscope = GlobalScope()]
    end
            vars = @variables begin
        C_bound2_ISM(t)
        C_bound_ISM(t)
        C_bound2_E7b(t)
        C_FcRn_ISM(t), [symscope = ParentScope(LocalScope())]
    end

    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_ISM = CL_up * tau_ISM   
    k_on_7 = (exg * k_on_7_EXG) + ((1.0 - exg) * k_on_7_EDG)
    KD_7 = (exg * KD_7_EXG) + ((1.0 - exg) * KD_7_EDG)     
    k_on_7_2 = k_on_7 / on_rate_ratio_1st_vs_2nd_binding # [1.0/(uM*h)]
    KD_7_2 = KD_7 * on_rate_ratio_1st_vs_2nd_binding # [nM] 
    k_off_7_2 = (KD_7_2 / 1000.0) * k_on_7_2 # [1.0/h]

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
        Dt(C_bound2_ISM) ~ (
            k_on_7_2 * C_bound_ISM * C_FcRn_ISM * V_ISM #  rxn on 
            - k_off_7_2 * C_bound2_ISM * V_ISM #  rxn off 
            + CL_up * (1.0 - FR) * C_bound2_E7b  #  from Endosomal 6 to interstitial membrane
            - CL_up * (1.0 - FR) * C_bound2_ISM #  to Endosomal 7
            - kdeg_FcRn_Ab * C_bound2_ISM * V_ISM  #  recover FcRn when IgG-FrRn is destroyed
         ) / V_ISM
    ]

    System(eqs, t, vars, pars; name=name)
end


#  ==============================================================================================
#  ==============================================================================================

#----------------------------------------#
#-----------------# FcRn #---------------#
#----------------------------------------#

#  Free FcRn in Organs

##  vascular side membrane 

function FCRN_VM(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac,  [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        tau_VM, [symscope = GlobalScope()]
        k_on_7_EXG, [symscope = GlobalScope()]
        k_on_7_EDG, [symscope = GlobalScope()]
        KD_7_EXG, [symscope = GlobalScope()]
        KD_7_EDG, [symscope = GlobalScope()] 
        kdeg_FcRn_Ab, [symscope = GlobalScope()]
        on_rate_ratio_1st_vs_2nd_binding, [symscope = GlobalScope()]
        FR, [symscope = GlobalScope()]
        FcRn_recycle_fraction, [symscope = GlobalScope()]
    end
            vars = @variables begin
        C_FcRn_VM(t), [symscope = ParentScope(LocalScope())]
        C_FcRn_E7b(t), [symscope = ParentScope(LocalScope())]
        C_VM_EDG(t)
        C_bound_VM_EDG(t)
        C_bound2_VM_EDG(t)
        C_VM_EXG(t)
        C_bound_VM_EXG(t)
        C_bound2_VM_EXG(t)
    end

    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_VM = CL_up * tau_VM
    k_off_7_EXG = (KD_7_EXG / 1000.0) * k_on_7_EXG # [1.0/h]
    k_off_7_EDG = (KD_7_EDG / 1000.0) * k_on_7_EDG # [1.0/h]
    k_on_7_EDG2 = k_on_7_EDG / on_rate_ratio_1st_vs_2nd_binding # [1.0/(uM*h)]
    k_on_7_EXG2 = k_on_7_EXG / on_rate_ratio_1st_vs_2nd_binding # [1.0/(uM*h)]
    KD_7_EDG2 = KD_7_EDG * on_rate_ratio_1st_vs_2nd_binding # [nM]
    KD_7_EXG2 = KD_7_EXG * on_rate_ratio_1st_vs_2nd_binding # [nM]
    k_off_7_EDG2 = (KD_7_EDG2 / 1000.0) * k_on_7_EDG2 # [1.0/h] 
    k_off_7_EXG2 = (KD_7_EXG2 / 1000.0) * k_on_7_EXG2 # [1.0/h] 


    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
            Dt(C_FcRn_VM) ~ (
              CL_up * FR * C_FcRn_E7b * (1.0 - FcRn_recycle_fraction) # from endosomal pH=7.4  b
            - CL_up * FR * C_FcRn_VM #  to vascular space
            #  EDG 
            + V_VM * (
                - k_on_7_EDG * C_VM_EDG * C_FcRn_VM # rxn on in vascular side membrane 
                + k_off_7_EDG * C_bound_VM_EDG # rxn off vascular side membrane 
                + kdeg_FcRn_Ab * C_bound_VM_EDG
                - k_on_7_EDG2 * C_bound_VM_EDG * C_FcRn_VM # rxn on in vascular side membrane 
                + k_off_7_EDG2 * C_bound2_VM_EDG # rxn off vascular side membrane 
                + 2.0 * kdeg_FcRn_Ab * C_bound2_VM_EDG # 
                ########## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
                #  EXG 
                - k_on_7_EXG * C_VM_EXG * C_FcRn_VM #  rxn on in vascular side membrane 
                + k_off_7_EXG * C_bound_VM_EXG #  rxn off vascular side membrane 
                + kdeg_FcRn_Ab * C_bound_VM_EXG
                - k_on_7_EXG2 * C_bound_VM_EXG * C_FcRn_VM #  rxn on in vascular side membrane 
                + k_off_7_EXG2 * C_bound2_VM_EXG #  rxn off vascular side membrane 
                + 2.0 * kdeg_FcRn_Ab * C_bound2_VM_EXG #  
            )
        ) / V_VM
    ]

    System(eqs, t, vars, pars; name=name)
end

#-----------#

# endosomal pH = 7.4

function FCRN_E7(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        pino_time, [symscope = GlobalScope()]
        k_on_7_EXG, [symscope = GlobalScope()]
        k_on_7_EDG, [symscope = GlobalScope()]
        KD_7_EXG, [symscope = GlobalScope()]
        KD_7_EDG, [symscope = GlobalScope()]
        on_rate_ratio_1st_vs_2nd_binding, [symscope = GlobalScope()]
        FR, [symscope = GlobalScope()]
        E6a_Vol_Pct, [symscope = GlobalScope()]
    end
        

    vars = @variables begin
        C_FcRn_E7(t), [symscope = ParentScope(LocalScope())]
        C_FcRn_VM(t), [symscope = ParentScope(LocalScope())]
        C_FcRn_ISM(t), [symscope = ParentScope(LocalScope())]
        C_E7_EDG(t)
        C_bound_E7_EDG(t)
        C_bound2_E7_EDG(t)
        C_E7_EXG(t)
        C_bound_E7_EXG(t)
        C_bound2_E7_EXG(t)
    end

    # derived params
    E7_Vol_Pct = (1.0 - E6a_Vol_Pct) / 2 # [-]
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_endosomal = CL_up_in_nL_per_hour_per_million_cells * pino_time * Organ_Endothelial_Cell * 1e-6 * 1e-9
    V_E7 = V_endosomal * E7_Vol_Pct
    k_off_7_EXG = (KD_7_EXG / 1000.0) * k_on_7_EXG # [1.0/h]
    k_off_7_EDG = (KD_7_EDG / 1000.0) * k_on_7_EXG # [1.0/h]
    k_on_7_EDG2 = k_on_7_EDG / on_rate_ratio_1st_vs_2nd_binding # [1.0/(uM*h)]
    k_on_7_EXG2 = k_on_7_EXG / on_rate_ratio_1st_vs_2nd_binding # [1.0/(uM*h)]
    KD_7_EDG2 = KD_7_EDG * on_rate_ratio_1st_vs_2nd_binding # [nM]
    KD_7_EXG2 = KD_7_EXG * on_rate_ratio_1st_vs_2nd_binding # [nM]
    k_off_7_EDG2 = (KD_7_EDG2 / 1000.0) * k_on_7_EDG2 # [1.0/h] 
    k_off_7_EXG2 = (KD_7_EXG2 / 1000.0) * k_on_7_EXG2 # [1.0/h] 

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
        Dt(C_FcRn_E7) ~ (
            CL_up * FR * C_FcRn_VM #  from vascular membrane
            + CL_up * (1.0-FR) * C_FcRn_ISM #  from is membrane
            - CL_up * C_FcRn_E7 #  to endosomal pH=6
            + V_E7 * (
                - k_on_7_EDG * C_E7_EDG * C_FcRn_E7 #  rxn on endosomal pH=7.4
                + k_off_7_EDG * C_bound_E7_EDG #  rxn off endosomal pH=7.4.
                - k_on_7_EDG2 * C_bound_E7_EDG * C_FcRn_E7 #  rxn on endosomal pH=7.4. 
                + k_off_7_EDG2 * C_bound2_E7_EDG #  rxn off endosomal pH=7.4. 
        #         # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
        #         #  EXG
                - k_on_7_EXG * C_E7_EXG * C_FcRn_E7 #  rxn on endosomal pH=7.4
                + k_off_7_EXG * C_bound_E7_EXG #  rxn off endosomal pH=7.4
                - k_on_7_EXG2 * C_bound_E7_EXG * C_FcRn_E7 #  rxn on endosomal pH=7.4. 
                + k_off_7_EXG2 * C_bound2_E7_EXG #  rxn off endosomal pH=7.4. 
            )
        ) / V_E7
    ]

    System(eqs, t, vars, pars; name=name)
end

#---------#

# E6a
function FCRN_E6A(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        pino_time, [symscope = GlobalScope()]
        E6a_Vol_Pct, [symscope = GlobalScope()]
        k_on_6_EXG, [symscope = GlobalScope()]
        k_on_6_EDG, [symscope = GlobalScope()]
        KD_6_EXG, [symscope = GlobalScope()]
        KD_6_EDG, [symscope = GlobalScope()]
        on_rate_ratio_1st_vs_2nd_binding, [symscope = GlobalScope()]
        FcRn_recycle_fraction, [symscope = GlobalScope()]
    end
            vars = @variables begin
        C_FcRn_E6a(t), [symscope = ParentScope(LocalScope())]
        C_FcRn_E7(t), [symscope = ParentScope(LocalScope())]
        C_FcRn_E7b(t), [symscope = ParentScope(LocalScope())]
        C_E6a_EDG(t)
        C_bound_E6a_EDG(t)
        C_bound2_E6a_EDG(t)
        C_E6a_EXG(t)
        C_bound_E6a_EXG(t)
        C_bound2_E6a_EXG(t)
    end

    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_endosomal = CL_up_in_nL_per_hour_per_million_cells * pino_time * Organ_Endothelial_Cell * 1e-6 * 1e-9
    V_E6a = V_endosomal * E6a_Vol_Pct
    k_off_6_EXG = (KD_6_EXG / 1000.0) * k_on_6_EXG # [1.0/h]
    k_off_6_EDG = (KD_6_EDG / 1000.0) * k_on_6_EDG # [1.0/h]
    k_on_6_EDG2 = k_on_6_EDG / on_rate_ratio_1st_vs_2nd_binding # [1.0/(uM*h)]
    k_on_6_EXG2 = k_on_6_EXG / on_rate_ratio_1st_vs_2nd_binding # [1.0/(uM*h)]
    KD_6_EDG2 = KD_6_EDG * on_rate_ratio_1st_vs_2nd_binding # [nM]
    KD_6_EXG2 = KD_6_EXG * on_rate_ratio_1st_vs_2nd_binding # [nM]
    k_off_6_EDG2 = (KD_6_EDG2 / 1000.0) * k_on_6_EDG2 # [1.0/h] 
    k_off_6_EXG2 = (KD_6_EXG2 / 1000.0) * k_on_6_EXG2 # [1.0/h]

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
        Dt(C_FcRn_E6a) ~ (
            CL_up * C_FcRn_E7 #  from E7
            - CL_up * C_FcRn_E6a #  from E6a to E7b
            + CL_up * C_FcRn_E7b * FcRn_recycle_fraction
            + V_E6a * (
                - k_on_6_EDG * C_E6a_EDG * C_FcRn_E6a #  rxn on endosomal pH=6.0
                + k_off_6_EDG * C_bound_E6a_EDG #  rxn off endosomal pH=6.0
                - k_on_6_EDG2 * C_bound_E6a_EDG * C_FcRn_E6a #  rxn on endosomal pH=6.0. 
                + k_off_6_EDG2 * C_bound2_E6a_EDG #  rxn off endosomal pH=6.0. 
                # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
                #  EXG
                - k_on_6_EXG * C_E6a_EXG * C_FcRn_E6a #  rxn on endosomal pH=6.0
                + k_off_6_EXG * C_bound_E6a_EXG #  rxn off endosomal pH=6.0
                - k_on_6_EXG2 * C_bound_E6a_EXG * C_FcRn_E6a #  rxn on endosomal pH=6.0. 
                + k_off_6_EXG2 * C_bound2_E6a_EXG #  rxn off endosomal pH=6.0. 
            )
        ) / V_E6a
    ]

    System(eqs, t, vars, pars; name=name)
end

#---------#

# endosomal pH = 7.4 b
function FCRN_E7B(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        pino_time, [symscope = GlobalScope()]
        E6a_Vol_Pct, [symscope = GlobalScope()]
        k_on_7_EXG, [symscope = GlobalScope()]
        k_on_7_EDG, [symscope = GlobalScope()]
        KD_7_EXG, [symscope = GlobalScope()]
        KD_7_EDG,   [symscope = GlobalScope()]
        on_rate_ratio_1st_vs_2nd_binding, [symscope = GlobalScope()]
    end
    vars = @variables begin
        C_FcRn_E7b(t), [symscope = ParentScope(LocalScope())]
        C_FcRn_E6a(t), [symscope = ParentScope(LocalScope())]
        C_E7b_EDG(t)
        C_bound_E7b_EDG(t)
        C_bound2_E7b_EDG(t)
        C_E7b_EXG(t)
        C_bound_E7b_EXG(t)
        C_bound2_E7b_EXG(t)
    end

    # derived params
    E7b_Vol_Pct = (1.0 - E6a_Vol_Pct) / 2 # [-]
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_endosomal = CL_up_in_nL_per_hour_per_million_cells * pino_time * Organ_Endothelial_Cell * 1e-6 * 1e-9
    V_E7b = V_endosomal * E7b_Vol_Pct
    k_off_7_EXG = (KD_7_EXG / 1000.0) * k_on_7_EXG # [1.0/h]
    k_off_7_EDG = (KD_7_EDG / 1000.0) * k_on_7_EDG # [1.0/h]
    k_on_7_EDG2 = k_on_7_EDG / on_rate_ratio_1st_vs_2nd_binding # [1.0/(uM*h)]
    k_on_7_EXG2 = k_on_7_EXG / on_rate_ratio_1st_vs_2nd_binding # [1.0/(uM*h)]
    KD_7_EDG2 = KD_7_EDG * on_rate_ratio_1st_vs_2nd_binding # [nM]
    KD_7_EXG2 = KD_7_EXG * on_rate_ratio_1st_vs_2nd_binding; # [nM]
    k_off_7_EDG2 = (KD_7_EDG2 / 1000.0) * k_on_7_EDG2 # [1.0/h] 
    k_off_7_EXG2 = (KD_7_EXG2 / 1000.0) * k_on_7_EXG2 # [1.0/h]

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
        Dt(C_FcRn_E7b) ~ (
            CL_up * C_FcRn_E6a #  from E6a
            - CL_up * C_FcRn_E7b #  from E7b to membranes
            + V_E7b * (
                - k_on_7_EDG * C_E7b_EDG * C_FcRn_E7b #  rxn on endosomal pH=7.4
                + k_off_7_EDG * C_bound_E7b_EDG #  rxn off endosomal pH=7.4
                - k_on_7_EDG2 * C_bound_E7b_EDG * C_FcRn_E7b #  rxn on endosomal pH=7.4. 
                + k_off_7_EDG2 * C_bound2_E7b_EDG #  rxn off endosomal pH=7.4. 
                # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
                #  EXG 
                - k_on_7_EXG * C_E7b_EXG * C_FcRn_E7b #  rxn on endosomal pH=7.4
                + k_off_7_EXG * C_bound_E7b_EXG #  rxn off endosomal pH=7.4
                - k_on_7_EXG2 * C_bound_E7b_EXG * C_FcRn_E7b #  rxn on endosomal pH=7.4. 
                + k_off_7_EXG2 * C_bound2_E7b_EXG #  rxn off endosomal pH=7.4. 
            )
        ) / V_E7b 
    ]

    System(eqs, t, vars, pars; name=name)
end

#---------#

# IS side mem
function FCRN_ISM(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        tau_ISM, [symscope = GlobalScope()]
        k_on_7_EXG, [symscope = GlobalScope()]
        k_on_7_EDG, [symscope = GlobalScope()]
        KD_7_EXG, [symscope = GlobalScope()]
        KD_7_EDG, [symscope = GlobalScope()]
        kdeg_FcRn_Ab, [symscope = GlobalScope()]
        on_rate_ratio_1st_vs_2nd_binding, [symscope = GlobalScope()]
        FR, [symscope = GlobalScope()]
        FcRn_recycle_fraction, [symscope = GlobalScope()]
    end
            
    vars = @variables begin
        C_FcRn_ISM(t), [symscope = ParentScope(LocalScope())]
        C_FcRn_E7b(t), [symscope = ParentScope(LocalScope())]
        C_ISM_EDG(t)
        C_bound_ISM_EDG(t)
        C_bound2_ISM_EDG(t)
        C_ISM_EXG(t)
        C_bound_ISM_EXG(t)
        C_bound2_ISM_EXG(t)
    end

    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    V_ISM = CL_up * tau_ISM
    k_off_7_EXG = (KD_7_EXG / 1000.0) * k_on_7_EXG # [1.0/h]
    k_off_7_EDG = (KD_7_EDG / 1000.0) * k_on_7_EDG # [1.0/h]
    k_on_7_EDG2 = k_on_7_EDG / on_rate_ratio_1st_vs_2nd_binding # [1.0/(uM*h)]
    k_on_7_EXG2 = k_on_7_EXG / on_rate_ratio_1st_vs_2nd_binding # [1.0/(uM*h)]
    KD_7_EDG2 = KD_7_EDG * on_rate_ratio_1st_vs_2nd_binding # [nM]
    KD_7_EXG2 = KD_7_EXG * on_rate_ratio_1st_vs_2nd_binding # [nM]
    k_off_7_EDG2 = (KD_7_EDG2 / 1000.0) * k_on_7_EDG2 # [1.0/h] 
    k_off_7_EXG2 = (KD_7_EXG2 / 1000.0) * k_on_7_EXG2 # [1.0/h] 

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
      #  IS side mem
      Dt(C_FcRn_ISM) ~ (
            -CL_up * (1.0 - FR) * C_FcRn_ISM #  from IS membrane to endosomal pH=7.4
            + CL_up * (1.0 - FR) * C_FcRn_E7b * (1.0 - FcRn_recycle_fraction) #  from endosomal pH=7.4  b to IS side mem
            + V_ISM * (
                - k_on_7_EDG * C_ISM_EDG * C_FcRn_ISM #  rxn on IS side mem
                + k_off_7_EDG * C_bound_ISM_EDG #  rxn off IS side mem
                + kdeg_FcRn_Ab * C_bound_ISM_EDG
                - k_on_7_EDG2 * C_bound_ISM_EDG * C_FcRn_ISM #  rxn on IS side mem. 
                + k_off_7_EDG2 * C_bound2_ISM_EDG #  rxn off IS side mem. 
                + 2.0 * kdeg_FcRn_Ab * C_bound2_ISM_EDG #  
                # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
                #  EXG 
                - k_on_7_EXG * C_ISM_EXG * C_FcRn_ISM #  rxn on IS side mem
                + k_off_7_EXG * C_bound_ISM_EXG #  rxn off IS side mem
                + kdeg_FcRn_Ab * C_bound_ISM_EXG
                - k_on_7_EXG2 * C_bound_ISM_EXG * C_FcRn_ISM #  rxn on IS side mem. 
                + k_off_7_EXG2 * C_bound2_ISM_EXG #  rxn off IS side mem. 
                + 2.0 * kdeg_FcRn_Ab * C_bound2_ISM_EXG #  
            )
        ) / V_ISM
    ]

    System(eqs, t, vars, pars; name=name)
end


#  ==============================================================================================   
# EYE FUNCTIONS

#-----------------------------------------#
#--------------# iris-ciliary body #-------------#
#-----------------------------------------#
# vascular space
function V_ICB(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        PLQ, [symscope = ParentScope(LocalScope())]
        V_V, [symscope = ParentScope(LocalScope())]
        sigma_V, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        FR, [symscope = GlobalScope()]
        LF, [symscope = ParentScope(LocalScope())]
        # added for iris-ciliary body
        RC_aq, [symscope = GlobalScope()]
        Q_AH, [symscope = GlobalScope()]
        Q_PtA, [symscope = GlobalScope()]
    end

    vars = @variables begin
        C_V(t)
        C_VM(t)
        C_V_Lung(t)
        # added 
        C_IntS(t)
    end

    # derived parameters
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15

    CL_up = CL_up * CLOff
    eqs = [
      Dt(C_V) ~ (
        PLQ * C_V_Lung  #  from Lung 
        - (PLQ - LF) * C_V  #  leave to Main Plasma
        - (1.0 - sigma_V) * LF * C_V #  going via Lymph to Interstitial 
        - CL_up * FR * C_V #  pinocytosis from vascular to vascular membrane
        + CL_up * FR * C_VM #  exocytosis from vascular membrane to vascular space 
        - C_V * Q_PtA * (1 -  RC_aq) 
        + C_IntS * Q_AH
      ) / V_V
    ]

    System(eqs, t, vars, pars; name=name)
end

#  Interstitial space
function INTS_ICB(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        sigma_IS, [symscope = ParentScope(LocalScope())]
        V_IntS, [symscope = ParentScope(LocalScope())]
        sigma_V,    [symscope = ParentScope(LocalScope())]
        FR, [symscope = GlobalScope()]
        LF, [symscope = ParentScope(LocalScope())]
        # added 
        Q_AH, [symscope = GlobalScope()]
    end
        
    vars = @variables begin
        C_IntS(t)
        C_V(t)
        C_ISM(t)
        C_AQ(t) 
    end

    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
            Dt(C_IntS) ~ (
            (1.0 - sigma_V) * LF * C_V        #  going from vascular via Lymph 
            - CL_up * (1.0 - FR) * C_IntS     #  pinocytosis 
            + CL_up * (1.0 - FR) * C_ISM      #  exocytosis
            - (1.0 - sigma_IS) * LF * C_IntS 
            # added 
            + C_AQ * Q_AH  
            - C_IntS * Q_AH     
        ) / V_IntS
    ]

    System(eqs, t, vars, pars; name=name)
end

#-----------------------------------------#
#--------------# retina #-------------#
#-----------------------------------------#
# vascular space is the same as other compartment
#  Interstitial space
function INTS_RETINA(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        sigma_IS, [symscope = ParentScope(LocalScope())]
        V_IntS, [symscope = ParentScope(LocalScope())]
        sigma_V,    [symscope = ParentScope(LocalScope())]
        FR, [symscope = GlobalScope()]
        LF, [symscope = ParentScope(LocalScope())]
        # added 
        PS_ret, [symscope = GlobalScope()]
        PS_cho, [symscope = GlobalScope()]
        RC_ret, [symscope = GlobalScope()]
        RC_cho, [symscope = GlobalScope()]
        Q_BF, [symscope = GlobalScope()]
    end
        
    vars = @variables begin
        C_IntS(t)
        C_V(t)
        C_ISM(t)
        # added 
        C_VH(t) # this was referred as C_vit/ C_Vit in the eqs
        C_IntS_CHO(t) # NEEDS ATTENTION
    end

    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
            Dt(C_IntS) ~ ((
            (1.0 - sigma_V) * LF * C_V        #  going from vascular via Lymph 
            - CL_up * (1.0 - FR) * C_IntS     #  pinocytosis 
            + CL_up * (1.0 - FR) * C_ISM      #  exocytosis
            - (1.0 - sigma_IS) * LF * C_IntS
            # added 
            + C_VH * Q_BF * (1 - RC_ret)  
            - C_IntS * Q_BF * (1 - RC_cho)
        ) / V_IntS 
        + (C_VH - C_IntS) * PS_ret  
        + ( C_IntS_CHO - C_IntS ) * PS_cho
            )
    ]

    System(eqs, t, vars, pars; name=name)
end



#-----------------------------------------#
#--------------# choriod #-------------#
#-----------------------------------------#
# vascular space is the same as other compartment
function V_CHOROID(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        PLQ, [symscope = ParentScope(LocalScope())]
        V_V, [symscope = ParentScope(LocalScope())]
        sigma_V, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        FR, [symscope = GlobalScope()]
        LF, [symscope = ParentScope(LocalScope())]
        # added 
        Q_BF, [symscope = GlobalScope()]
    end

    vars = @variables begin
        C_V(t)
        C_VM(t)
        C_V_Lung(t)
        # added 
        C_IntS(t)
    end

    # derived parameters
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15

    CL_up = CL_up * CLOff
    eqs = [
      Dt(C_V) ~ (
        PLQ * C_V_Lung  #  from Lung 
        - (PLQ - LF) * C_V  #  leave to Main Plasma
        - (1.0 - sigma_V) * LF * C_V #  going via Lymph to Interstitial 
        - CL_up * FR * C_V #  pinocytosis from vascular to vascular membrane
        + CL_up * FR * C_VM #  exocytosis from vascular membrane to vascular space 
        + C_IntS * Q_BF 
      ) / V_V
    ]

    System(eqs, t, vars, pars; name=name)
end

#  Interstitial space
function INTS_CHOROID(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        sigma_IS, [symscope = ParentScope(LocalScope())]
        V_IntS, [symscope = ParentScope(LocalScope())]
        sigma_V,    [symscope = ParentScope(LocalScope())]
        FR, [symscope = GlobalScope()]
        LF, [symscope = ParentScope(LocalScope())]
        # added 
        PS_cho, [symscope = GlobalScope()]
        RC_cho, [symscope = GlobalScope()]
        Q_BF, [symscope = GlobalScope()]
    end
        
    vars = @variables begin
        C_IntS(t)
        C_V(t)
        C_ISM(t)
        # added 
        C_IntS_RET(t) # NEEDS ATTENTION
    end

    # derived params
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15

    CL_up = CL_up * CLOff # apply clearance off from vascular space to interstitial space
    eqs = [
            Dt(C_IntS) ~ ((
            (1.0 - sigma_V) * LF * C_V        #  going from vascular via Lymph 
            - CL_up * (1.0 - FR) * C_IntS     #  pinocytosis 
            + CL_up * (1.0 - FR) * C_ISM      #  exocytosis
            - (1.0 - sigma_IS) * LF * C_IntS
            # added 
            + C_IntS_RET * Q_BF * (1 - RC_cho) 
            - C_IntS * Q_BF
        ) / V_IntS
        + ( C_IntS_RET - C_IntS ) * PS_cho
            )
    ]

    System(eqs, t, vars, pars; name=name)
end

#-----------------------------------------#
#--------------# cornea #-------------#
#-----------------------------------------#
function CORNEA(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        PS_cor, [symscope = GlobalScope()]
        Kp_adc_cor_ah, [symscope = GlobalScope()]
        V_Cor, [symscope = GlobalScope()]
    end
        
    vars = @variables begin
        C_COR(t)
        C_AQ(t)
    end

    eqs = [
            Dt(C_COR) ~ (C_AQ - C_COR/Kp_adc_cor_ah) * PS_cor / V_Cor
    ]

    System(eqs, t, vars, pars; name=name)
end

#-----------------------------------------#
#--------------# lens #-------------#
#-----------------------------------------#
function LENS(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        PS_lens, [symscope = GlobalScope()]
        V_Lens, [symscope = GlobalScope()]
    end
        
    vars = @variables begin
        C_LENS(t)
        C_AQ(t)
        C_VH(t) # this was referred as C_vit/ C_Vit in the eqs
    end

    eqs = [
            Dt(C_LENS) ~ (
                  (C_AQ - C_LENS) * PS_lens
                + (C_VH - C_LENS) * PS_lens
            )/V_Lens
    ]

    System(eqs, t, vars, pars; name=name)
end

#-----------------------------------------#
#--------------# aqueous humor compartment #-------------#
#-----------------------------------------#

function AQUEOUSHUMOR(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        PS_cor, [symscope = GlobalScope()]
        PS_lens, [symscope = GlobalScope()]
        RC_aq, [symscope = GlobalScope()]
        Q_AH, [symscope = GlobalScope()]
        Q_PtA, [symscope = GlobalScope()]
        Q_BF, [symscope = GlobalScope()]
        V_AQ, [symscope = GlobalScope()]
    end
        
    vars = @variables begin
        C_COR(t)
        C_AQ(t)
        C_LENS(t)
        C_VH(t) # this was referred as C_vit/ C_Vit in the eqs
        C_V_ICB(t)  # NEED ATTENTION
    end

    # derived variable 
    Q_VA = Q_AH + Q_BF - Q_PtA

    eqs = [
            Dt(C_AQ) ~ (
                  (C_COR - C_AQ) * PS_cor / V_AQ
                + (C_LENS - C_AQ) * PS_lens / V_AQ
                + (C_VH * Q_VA) / V_AQ 
                + (C_V_ICB * Q_PtA * (1 - RC_aq)) / V_AQ
                - ( C_AQ * (Q_AH + Q_BF + Q_VA)  ) / V_AQ
            )
    ]

    System(eqs, t, vars, pars; name=name)
end

#-----------------------------------------#
#--------------# vitreous humor compartment #-------------#
#-----------------------------------------#

function VITREOUSHUMOR(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        PS_ret, [symscope = GlobalScope()]
        PS_lens, [symscope = GlobalScope()]
        RC_ret, [symscope = GlobalScope()]
        Q_AH, [symscope = GlobalScope()]
        Q_PtA, [symscope = GlobalScope()]
        Q_BF, [symscope = GlobalScope()]
        V_VH, [symscope = GlobalScope()]
    end
        
    vars = @variables begin
        C_VH(t) # this was referred as C_vit/ C_Vit in the eqs
        C_IntS_RET(t)  # NEED ATTENTION
        C_LENS(t)
        C_AQ(t)
    end

    # derived variable 
    Q_VA = Q_AH + Q_BF - Q_PtA

    eqs = [
            Dt(C_VH) ~ (
                  (C_IntS_RET - C_VH) * PS_ret / V_VH
                + (C_LENS - C_VH) * PS_lens / V_VH
                - (C_VH * Q_VA) / V_VH 
                - (C_VH * Q_BF * (1 - RC_ret)) / V_VH
                + ( C_AQ * (Q_BF + Q_VA)  ) / V_VH
            )
    ]

    System(eqs, t, vars, pars; name=name)
end

#  ==============================================================================================   

#-----------------------------------------#
#--------------# lymph node #-------------#
#-----------------------------------------#

# flux of IgG in Lymph node

function LN(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        sigma_IS_Lung, [symscope = GlobalScope()]
        sigma_IS_SI, [symscope = GlobalScope()]
        sigma_IS_LI, [symscope = GlobalScope()]
        sigma_IS_Liver, [symscope = GlobalScope()]
        sigma_IS_Brain, [symscope = GlobalScope()]
        sigma_IS_Other, [symscope = GlobalScope()]
        PLQ_Adipose, [symscope = GlobalScope()]
        sigma_IS_Pancreas, [symscope = GlobalScope()]
        sigma_IS_Muscle, [symscope = GlobalScope()]
        PLQ_SI, [symscope = GlobalScope()]
        PLQ_Heart, [symscope = GlobalScope()]
        sigma_IS_Kidney, [symscope = GlobalScope()]
        V_LN, [symscope = GlobalScope()]
        PLQ_Bone, [symscope = GlobalScope()]
        PLQ_LI, [symscope = GlobalScope()]
        PLQ_Liver, [symscope = GlobalScope()]
        PLQ_Brain, [symscope = GlobalScope()]
        PLQ_Spleen, [symscope = GlobalScope()]
        PLQ_Thymus, [symscope = GlobalScope()]
        PLQ_Lung, [symscope = GlobalScope()]
        sigma_IS_Thymus, [symscope = GlobalScope()]
        PLQ_Skin, [symscope = GlobalScope()]
        sigma_IS_Skin, [symscope = GlobalScope()]
        PLQ_Muscle, [symscope = GlobalScope()]
        PLQ_Pancreas, [symscope = GlobalScope()]
        sigma_IS_Adipose, [symscope = GlobalScope()]
        sigma_IS_Heart, [symscope = GlobalScope()]
        sigma_IS_Spleen, [symscope = GlobalScope()]
        PLQ_Kidney, [symscope = GlobalScope()]
        PLQ_Other, [symscope = GlobalScope()]
        sigma_IS_Bone,  [symscope = GlobalScope()]
        LF_Heart, [symscope = GlobalScope()]
        LF_Kidney, [symscope = GlobalScope()]
        LF_Muscle, [symscope = GlobalScope()]
        LF_Skin, [symscope = GlobalScope()] 
        LF_Brain, [symscope = GlobalScope()]
        LF_Adipose, [symscope = GlobalScope()]
        LF_Thymus, [symscope = GlobalScope()]
        LF_Liver, [symscope = GlobalScope()]
        LF_Spleen, [symscope = GlobalScope()]
        LF_Pancreas, [symscope = GlobalScope()]
        LF_SI, [symscope = GlobalScope()]
        LF_LI, [symscope = GlobalScope()]
        LF_Bone, [symscope = GlobalScope()]
        LF_Other, [symscope = GlobalScope()]
        LF_Lung, [symscope = GlobalScope()]
        L_LymphNode, [symscope = GlobalScope()]
        # added for eye
        sigma_IS_ICB, [symscope = GlobalScope()]
        sigma_IS_Retina, [symscope = GlobalScope()]
        sigma_IS_Choroid, [symscope = GlobalScope()]
        PLQ_ICB, [symscope = GlobalScope()]
        PLQ_Retina, [symscope = GlobalScope()]
        PLQ_Choroid, [symscope = GlobalScope()]
        LF_ICB, [symscope = GlobalScope()]
        LF_Retina, [symscope = GlobalScope()]
        LF_Choroid, [symscope = GlobalScope()]
        # added for bone marrow 
        sigma_IS_Marrow, [symscope = GlobalScope()]
        PLQ_Marrow, [symscope = GlobalScope()]
        LF_Marrow, [symscope = GlobalScope()]
    end
        vars = @variables begin
        C_LN(t)
        C_IntS_Heart(t)
        C_IntS_Kidney(t)
        C_IntS_Muscle(t)
        C_IntS_Skin(t)
        C_IntS_Brain(t)
        C_IntS_Adipose(t)
        C_IntS_Thymus(t)
        C_IntS_Liver(t)
        C_IntS_Spleen(t)
        C_IntS_Pancreas(t)
        C_IntS_SI(t)
        C_IntS_LI(t)
        C_IntS_Bone(t)
        C_IntS_Other(t)
        C_IntS_Lung(t)
        # added, eye 
        C_IntS_ICB(t)
        C_IntS_Retina(t)
        C_IntS_Choroid(t)
        # added, bonr marrow 
        C_IntS_Marrow(t)
    end

    eqs = [
        Dt(C_LN) ~ (
            (1.0 - sigma_IS_Heart) * LF_Heart * C_IntS_Heart
            + (1.0 - sigma_IS_Kidney) * LF_Kidney * C_IntS_Kidney
            + (1.0 - sigma_IS_Muscle) * LF_Muscle * C_IntS_Muscle
            + (1.0 - sigma_IS_Skin) * LF_Skin * C_IntS_Skin
            + (1.0 - sigma_IS_Brain) * LF_Brain * C_IntS_Brain
            + (1.0 - sigma_IS_Adipose) * LF_Adipose * C_IntS_Adipose
            + (1.0 - sigma_IS_Thymus) * LF_Thymus * C_IntS_Thymus
            + (1.0 - sigma_IS_Liver) * LF_Liver * C_IntS_Liver
            + (1.0 - sigma_IS_Spleen) * LF_Spleen * C_IntS_Spleen
            + (1.0 - sigma_IS_Pancreas) * LF_Pancreas * C_IntS_Pancreas
            + (1.0 - sigma_IS_SI) * LF_SI * C_IntS_SI
            + (1.0 - sigma_IS_LI) * LF_LI * C_IntS_LI
            + (1.0 - sigma_IS_Bone) * LF_Bone * C_IntS_Bone
            + (1.0 - sigma_IS_Other) * LF_Other * C_IntS_Other
            + (1.0 - sigma_IS_Lung) * LF_Lung * C_IntS_Lung
            # added 
            + (1.0 - sigma_IS_ICB) * LF_ICB * C_IntS_ICB
            + (1.0 - sigma_IS_Retina) * LF_Retina * C_IntS_Retina
            + (1.0 - sigma_IS_Choroid) * LF_Choroid * C_IntS_Choroid
            + (1.0 - sigma_IS_Marrow) * LF_Marrow * C_IntS_Marrow
            #
            - L_LymphNode * C_LN
        ) / V_LN
    ]

    # listing global parameters here since there is no param block in this function
    System(eqs, t, vars, pars; name=name)
end


#  ==============================================================================================   

#-----------------------------------------#
#---------------# plasma #----------------#
#-----------------------------------------#

function PLASMA(; name, exg)
    exg = Float64.(exg)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        PLQ_Adipose, [symscope = GlobalScope()]
        PLQ_Liver, [symscope = GlobalScope()]
        PLQ_Skin, [symscope = GlobalScope()]
        PLQ_Pancreas, [symscope = GlobalScope()]
        PLQ_Muscle, [symscope = GlobalScope()]
        PLQ_Brain, [symscope = GlobalScope()]
        V_Plasma, [symscope = GlobalScope()]
        PLQ_SI, [symscope = GlobalScope()]
        PLQ_Spleen, [symscope = GlobalScope()]
        PLQ_Heart, [symscope = GlobalScope()]
        PLQ_Thymus, [symscope = GlobalScope()]
        PLQ_Kidney, [symscope = GlobalScope()]
        PLQ_Other, [symscope = GlobalScope()]
        PLQ_Lung, [symscope = GlobalScope()]
        PLQ_Bone, [symscope = GlobalScope()]
        PLQ_LI, [symscope = GlobalScope()] 
        LF_Heart, [symscope = GlobalScope()]
        LF_Kidney, [symscope = GlobalScope()]
        LF_Muscle, [symscope = GlobalScope()]
        LF_Skin, [symscope = GlobalScope()]
        LF_Brain, [symscope = GlobalScope()]
        LF_Adipose, [symscope = GlobalScope()]
        LF_Thymus, [symscope = GlobalScope()]
        LF_Liver, [symscope = GlobalScope()]
        LF_Spleen, [symscope = GlobalScope()]
        LF_Pancreas, [symscope = GlobalScope()]
        LF_SI, [symscope = GlobalScope()]
        LF_LI, [symscope = GlobalScope()] 
        LF_Bone, [symscope = GlobalScope()]
        LF_Other, [symscope = GlobalScope()]
        LF_Lung, [symscope = GlobalScope()]
        L_LymphNode, [symscope = GlobalScope()]
        # added for eye
        PLQ_ICB, [symscope = GlobalScope()]
        PLQ_Retina, [symscope = GlobalScope()]
        PLQ_Choroid, [symscope = GlobalScope()]
        LF_ICB, [symscope = GlobalScope()]
        LF_Retina, [symscope = GlobalScope()]
        LF_Choroid, [symscope = GlobalScope()]
        # added for bone marrow 
        PLQ_Marrow, [symscope = GlobalScope()]
        LF_Marrow, [symscope = GlobalScope()]
        # added for infusion 
        infusion(t), [symscope = GlobalScope()]
        # added for non-linear clearance in plasma 
        VMAX, [symscope = GlobalScope()]
        KM, [symscope = GlobalScope()]
    end
    vars = @variables begin
        C_Plasma(t)
        C_LN(t)
        C_V_Heart(t)
        C_V_Kidney(t)
        C_V_Muscle(t)
        C_V_Skin(t)
        C_V_Brain(t)
        C_V_Adipose(t)
        C_V_Thymus(t)
        C_V_Liver(t)
        C_V_Spleen(t)
        C_V_Pancreas(t)
        C_V_SI(t)
        C_V_LI(t)
        C_V_Bone(t)
        C_V_Other(t)
        C_V_Lung(t)
        # added for eye
        C_V_ICB(t)
        C_V_Retina(t)
        C_V_Choroid(t)
        # added for bone marrow
        C_V_Marrow(t)
    end

    eqs = [
        Dt(C_Plasma) ~ Float64(exg) * (
            (PLQ_Heart - LF_Heart) * C_V_Heart
            + (PLQ_Kidney - LF_Kidney) * C_V_Kidney
            + (PLQ_Muscle - LF_Muscle) * C_V_Muscle
            + (PLQ_Skin - LF_Skin) * C_V_Skin
            + (PLQ_Brain - LF_Brain) * C_V_Brain
            + (PLQ_Adipose - LF_Adipose) * C_V_Adipose
            + (PLQ_Thymus - LF_Thymus) * C_V_Thymus
            + (PLQ_Liver - LF_Liver) * C_V_Liver
            + (PLQ_Spleen - LF_Spleen) * C_V_Liver
            + (PLQ_Pancreas - LF_Pancreas) * C_V_Liver
            + (PLQ_SI - LF_SI) * C_V_Liver
            + (PLQ_LI - LF_LI) * C_V_Liver
            + (PLQ_Bone - LF_Bone) * C_V_Bone
            + (PLQ_Other - LF_Other) * C_V_Other
            # added 
            + (PLQ_ICB - LF_ICB) * C_V_ICB
            + (PLQ_Retina - LF_Retina) * C_V_Retina
            + (PLQ_Choroid - LF_Choroid) * C_V_Choroid
            + (PLQ_Marrow - LF_Marrow) * C_V_Marrow
            # 
            - (PLQ_Lung + LF_Lung) * C_Plasma
            + L_LymphNode * C_LN
        ) / V_Plasma + infusion * Float64(exg) - (VMAX/V_Plasma) / (C_Plasma + KM) * C_Plasma
    ]

    name_symbol = Bool(exg) ? Symbol("plasma_exg") : Symbol("plasma_edg")

    # listing global parameters here since there is no param block in this function
    System(eqs, t, vars, pars; name=name_symbol)
end

#-----------------------------------------#
#--------------# payload #-------------#
#-----------------------------------------#

# payload in endothelial cells 

function PL_organ_endothelial(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        CLOff, [symscope = GlobalScope()] #  clearance off from vascular space to interstitial space
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        pino_time, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        tau_VM , [symscope = GlobalScope()]
        tau_ISM, [symscope = GlobalScope()]
        DAR, [symscope = GlobalScope()]
        Prob_deg, [symscope = GlobalScope()]
        kdeg_FcRn_Ab, [symscope = GlobalScope()]
        frac_lys, [symscope = GlobalScope()]
        k_diff, [symscope = GlobalScope()]
        Kp_PL_endo_plasma, [symscope = GlobalScope()]
        k_out_ints, [symscope = ParentScope(LocalScope())]
        k_in_ints, [symscope = ParentScope(LocalScope())]
    end

    vars = @variables begin
        # cytosol PL concentration, endothelial cells
        C_PL_endo(t)
        # PL concentration, tissue interstitium
        C_PL_IntS(t)
        # PL concentration, plasma 
        C_PL_Plasma(t)  
        # all mAb/ ADC degredation related concentrations 
        C_E6a(t)
        C_bound_VM(t)
        C_bound2_VM(t)
        C_bound_ISM(t)
        C_bound2_ISM(t)        
    end

    # derived parameters
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    CL_up = CL_up_in_nL_per_hour_per_million_cells * Organ_Endothelial_Cell * 1e-15
    CL_up = CL_up * CLOff

    V_VM = CL_up * tau_VM
    V_ISM = CL_up * tau_ISM
    V_endosomal = CL_up_in_nL_per_hour_per_million_cells * pino_time * Organ_Endothelial_Cell * 1e-6 * 1e-9

    eqs = [
      Dt(C_PL_endo) ~ (
        (CL_up * Prob_deg * C_E6a
        + kdeg_FcRn_Ab * C_bound_VM * V_VM
        + kdeg_FcRn_Ab * C_bound2_VM * V_VM
        + kdeg_FcRn_Ab * C_bound_ISM * V_ISM
        + kdeg_FcRn_Ab * C_bound2_ISM * V_ISM    
      ) / V_endosomal * DAR * frac_lys
      - k_out_ints * C_PL_endo
      + k_in_ints * C_PL_IntS
      - k_diff * (C_PL_endo/Kp_PL_endo_plasma - C_PL_IntS) 
      - k_diff * (C_PL_endo/Kp_PL_endo_plasma - C_PL_Plasma)
      )
    ]

    System(eqs, t, vars, pars; name=name)
end


# PL in tissue interstitium 

function PL_tissue_interstitium(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        Endothelial_Cell_Frac, [symscope = ParentScope(LocalScope())]
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]
        pino_time, [symscope = GlobalScope()]
        V_IntS, [symscope = ParentScope(LocalScope())]
        k_diff, [symscope = GlobalScope()]
        Kp_PL_endo_plasma, [symscope = GlobalScope()]
        k_out_ints, [symscope = ParentScope(LocalScope())]
        k_in_ints, [symscope = ParentScope(LocalScope())]
        k_PL_ex, [symscope = ParentScope(LocalScope())] # PL clearance rate in tissue interstitium
    end

    vars = @variables begin
        # cytosol PL concentration, endothelial cells
        C_PL_endo(t)
        # PL concentration, tissue interstitium
        C_PL_IntS(t)    
    end

    # derived parameters
    Organ_Endothelial_Cell = Total_Endothelial_Cell * Endothelial_Cell_Frac * Scale_Factor
    V_endosomal = CL_up_in_nL_per_hour_per_million_cells * pino_time * Organ_Endothelial_Cell * 1e-6 * 1e-9

    flux_PL_endothelial2interstitium = (k_out_ints * C_PL_endo - k_in_ints * C_PL_IntS + k_diff * (C_PL_endo/Kp_PL_endo_plasma - C_PL_IntS) ) * V_endosomal  # unit: umol

    eqs = [
      Dt(C_PL_IntS) ~ (
        flux_PL_endothelial2interstitium / V_IntS - k_PL_ex * C_PL_IntS
      ) 
    ]

    System(eqs, t, vars, pars; name=name)
end


function create_organ_PL()
    # create subcompartment objects
    @independent_variables t
    Dt = Differential(t)
    PL_endo = PL_organ_endothelial(; name = :PL_endo)
    PL_Ints = PL_tissue_interstitium(; name = :PL_Ints)

    ModelingToolkit.extend(System(Equation[], t, name=:PL_tissue), [PL_endo, PL_Ints])
end

# payload in plasma 
function PL_plasma(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        CL_up_in_nL_per_hour_per_million_cells, [symscope = GlobalScope()]
        pino_time, [symscope = GlobalScope()]
        Total_Endothelial_Cell, [symscope = GlobalScope()]
        Scale_Factor, [symscope = GlobalScope()]

        Endothelial_Cell_Frac_Lung, [symscope = GlobalScope()]
        Endothelial_Cell_Frac_Liver, [symscope = GlobalScope()]
        Endothelial_Cell_Frac_Heart, [symscope = GlobalScope()]
        Endothelial_Cell_Frac_Muscle, [symscope = GlobalScope()]
        Endothelial_Cell_Frac_Skin, [symscope = GlobalScope()]
        Endothelial_Cell_Frac_Adipose, [symscope = GlobalScope()]
        Endothelial_Cell_Frac_Bone, [symscope = GlobalScope()]
        Endothelial_Cell_Frac_Brain, [symscope = GlobalScope()]
        Endothelial_Cell_Frac_Kidney, [symscope = GlobalScope()]
        Endothelial_Cell_Frac_SI, [symscope = GlobalScope()]
        Endothelial_Cell_Frac_LI, [symscope = GlobalScope()]
        Endothelial_Cell_Frac_Pancreas, [symscope = GlobalScope()]
        Endothelial_Cell_Frac_Thymus, [symscope = GlobalScope()]
        Endothelial_Cell_Frac_Spleen, [symscope = GlobalScope()]
        Endothelial_Cell_Frac_Other, [symscope = GlobalScope()]
        Endothelial_Cell_Frac_Marrow, [symscope = GlobalScope()]
        Endothelial_Cell_Frac_ICB, [symscope = GlobalScope()]
        Endothelial_Cell_Frac_Retina, [symscope = GlobalScope()]
        Endothelial_Cell_Frac_Choroid, [symscope = GlobalScope()]
        
        k_diff, [symscope = GlobalScope()]
        Kp_PL_endo_plasma, [symscope = GlobalScope()]
        CL_plasma_PL, [symscope = GlobalScope()] # PL clearance rate in plasma
        V_Plasma, [symscope = GlobalScope()] # plasma volume

        Q_PtA, [symscope = GlobalScope()]
        DAR, [symscope = GlobalScope()]
        VMAX, [symscope = GlobalScope()]
        KM, [symscope = GlobalScope()]
    end

    vars = @variables begin
        # cytosol PL concentration, endothelial cells
        C_PL_endo_Lung(t)
        C_PL_endo_Liver(t)
        C_PL_endo_Heart(t)
        C_PL_endo_Muscle(t)
        C_PL_endo_Skin(t)
        C_PL_endo_Adipose(t)
        C_PL_endo_Bone(t)
        C_PL_endo_Brain(t)
        C_PL_endo_Kidney(t)
        C_PL_endo_SI(t)
        C_PL_endo_LI(t)
        C_PL_endo_Pancreas(t)
        C_PL_endo_Thymus(t)
        C_PL_endo_Spleen(t)
        C_PL_endo_Other(t)
        C_PL_endo_Marrow(t)
        C_PL_endo_ICB(t)
        C_PL_endo_Retina(t)
        C_PL_endo_Choroid(t)
        # PL concentration, plasma 
        C_PL_Plasma(t)   
        # ADC concentration, plasma 
        C_Plasma(t) 
    end

    # derived parameters
    V_endosomal_Lung = CL_up_in_nL_per_hour_per_million_cells * pino_time * Total_Endothelial_Cell * Endothelial_Cell_Frac_Lung * Scale_Factor * 1e-15 
    V_endosomal_Liver = CL_up_in_nL_per_hour_per_million_cells * pino_time * Total_Endothelial_Cell * Endothelial_Cell_Frac_Liver * Scale_Factor * 1e-15 
    V_endosomal_Heart = CL_up_in_nL_per_hour_per_million_cells * pino_time * Total_Endothelial_Cell * Endothelial_Cell_Frac_Heart * Scale_Factor * 1e-15 
    V_endosomal_Muscle = CL_up_in_nL_per_hour_per_million_cells * pino_time * Total_Endothelial_Cell * Endothelial_Cell_Frac_Muscle * Scale_Factor * 1e-15 
    V_endosomal_Skin = CL_up_in_nL_per_hour_per_million_cells * pino_time * Total_Endothelial_Cell * Endothelial_Cell_Frac_Skin * Scale_Factor * 1e-15 
    V_endosomal_Adipose = CL_up_in_nL_per_hour_per_million_cells * pino_time * Total_Endothelial_Cell * Endothelial_Cell_Frac_Adipose * Scale_Factor * 1e-15 
    V_endosomal_Bone = CL_up_in_nL_per_hour_per_million_cells * pino_time * Total_Endothelial_Cell * Endothelial_Cell_Frac_Bone * Scale_Factor * 1e-15 
    V_endosomal_Brain = CL_up_in_nL_per_hour_per_million_cells * pino_time * Total_Endothelial_Cell * Endothelial_Cell_Frac_Brain * Scale_Factor * 1e-15 
    V_endosomal_Kidney = CL_up_in_nL_per_hour_per_million_cells * pino_time * Total_Endothelial_Cell * Endothelial_Cell_Frac_Kidney * Scale_Factor * 1e-15 
    V_endosomal_SI = CL_up_in_nL_per_hour_per_million_cells * pino_time * Total_Endothelial_Cell * Endothelial_Cell_Frac_SI * Scale_Factor * 1e-15 
    V_endosomal_LI = CL_up_in_nL_per_hour_per_million_cells * pino_time * Total_Endothelial_Cell * Endothelial_Cell_Frac_LI * Scale_Factor * 1e-15 
    V_endosomal_Pancreas = CL_up_in_nL_per_hour_per_million_cells * pino_time * Total_Endothelial_Cell * Endothelial_Cell_Frac_Pancreas * Scale_Factor * 1e-15 
    V_endosomal_Thymus = CL_up_in_nL_per_hour_per_million_cells * pino_time * Total_Endothelial_Cell * Endothelial_Cell_Frac_Thymus * Scale_Factor * 1e-15 
    V_endosomal_Spleen = CL_up_in_nL_per_hour_per_million_cells * pino_time * Total_Endothelial_Cell * Endothelial_Cell_Frac_Spleen * Scale_Factor * 1e-15 
    V_endosomal_Other = CL_up_in_nL_per_hour_per_million_cells * pino_time * Total_Endothelial_Cell * Endothelial_Cell_Frac_Other * Scale_Factor * 1e-15 
    V_endosomal_Marrow = CL_up_in_nL_per_hour_per_million_cells * pino_time * Total_Endothelial_Cell * Endothelial_Cell_Frac_Marrow * Scale_Factor * 1e-15 
    V_endosomal_ICB = CL_up_in_nL_per_hour_per_million_cells * pino_time * Total_Endothelial_Cell * Endothelial_Cell_Frac_ICB * Scale_Factor * 1e-15 
    V_endosomal_Retina = CL_up_in_nL_per_hour_per_million_cells * pino_time * Total_Endothelial_Cell * Endothelial_Cell_Frac_Retina * Scale_Factor * 1e-15 
    V_endosomal_Choroid = CL_up_in_nL_per_hour_per_million_cells * pino_time * Total_Endothelial_Cell * Endothelial_Cell_Frac_Choroid * Scale_Factor * 1e-15 

    # fluxes, PL diffuse out of endothelial cells 
    flux_diff_lung = k_diff * (C_PL_endo_Lung/Kp_PL_endo_plasma - C_PL_Plasma) * V_endosomal_Lung
    flux_diff_liver = k_diff * (C_PL_endo_Liver/Kp_PL_endo_plasma - C_PL_Plasma) * V_endosomal_Liver
    flux_diff_heart = k_diff * (C_PL_endo_Heart/Kp_PL_endo_plasma - C_PL_Plasma) * V_endosomal_Heart
    flux_diff_muscle = k_diff * (C_PL_endo_Muscle/Kp_PL_endo_plasma - C_PL_Plasma) * V_endosomal_Muscle
    flux_diff_skin = k_diff * (C_PL_endo_Skin/Kp_PL_endo_plasma - C_PL_Plasma) * V_endosomal_Skin
    flux_diff_adipose = k_diff * (C_PL_endo_Adipose/Kp_PL_endo_plasma - C_PL_Plasma) * V_endosomal_Adipose
    flux_diff_bone = k_diff * (C_PL_endo_Bone/Kp_PL_endo_plasma - C_PL_Plasma) * V_endosomal_Bone
    flux_diff_brain = k_diff * (C_PL_endo_Brain/Kp_PL_endo_plasma - C_PL_Plasma) * V_endosomal_Brain
    flux_diff_kidney = k_diff * (C_PL_endo_Kidney/Kp_PL_endo_plasma - C_PL_Plasma) * V_endosomal_Kidney
    flux_diff_si = k_diff * (C_PL_endo_SI/Kp_PL_endo_plasma - C_PL_Plasma) * V_endosomal_SI
    flux_diff_li = k_diff * (C_PL_endo_LI/Kp_PL_endo_plasma - C_PL_Plasma) * V_endosomal_LI
    flux_diff_pancreas = k_diff * (C_PL_endo_Pancreas/Kp_PL_endo_plasma - C_PL_Plasma) * V_endosomal_Pancreas
    flux_diff_thymus = k_diff * (C_PL_endo_Thymus/Kp_PL_endo_plasma - C_PL_Plasma) * V_endosomal_Thymus
    flux_diff_spleen = k_diff * (C_PL_endo_Spleen/Kp_PL_endo_plasma - C_PL_Plasma) * V_endosomal_Spleen
    flux_diff_other = k_diff * (C_PL_endo_Other/Kp_PL_endo_plasma - C_PL_Plasma) * V_endosomal_Other
    flux_diff_marrow = k_diff * (C_PL_endo_Marrow/Kp_PL_endo_plasma - C_PL_Plasma) * V_endosomal_Marrow
    flux_diff_icb = k_diff * (C_PL_endo_ICB/Kp_PL_endo_plasma - C_PL_Plasma) * V_endosomal_ICB
    flux_diff_retina = k_diff * (C_PL_endo_Retina/Kp_PL_endo_plasma - C_PL_Plasma) * V_endosomal_Retina
    flux_diff_choroid = k_diff * (C_PL_endo_Choroid/Kp_PL_endo_plasma - C_PL_Plasma) * V_endosomal_Choroid

    flux_diff = (flux_diff_lung + flux_diff_liver + flux_diff_heart + flux_diff_muscle + flux_diff_skin + flux_diff_adipose + flux_diff_bone + flux_diff_brain
                 + flux_diff_kidney + flux_diff_si + flux_diff_li + flux_diff_pancreas + flux_diff_thymus + flux_diff_spleen + flux_diff_other + 
                 flux_diff_marrow + flux_diff_icb + flux_diff_retina + flux_diff_choroid)

    eqs = [
      Dt(C_PL_Plasma) ~ (
         (flux_diff - C_PL_Plasma * Q_PtA)/V_Plasma 
         + DAR * (VMAX/V_Plasma) / (C_Plasma + KM) * C_Plasma   # from nonspecific ADC clearance
         - CL_plasma_PL * C_PL_Plasma                           # clearance
      ) 
    ]

    System(eqs, t, vars, pars; name=name)

end

# payload in eye 
## assuming the dynamics is comparable to mAb, with the same fluid flow and diffusion to drive PL distribution 
function PL_eye(; name)
    @independent_variables t
    Dt = Differential(t)
    pars = @parameters begin
        # fluid flow 
        Q_AH, [symscope = GlobalScope()]
        Q_PtA, [symscope = GlobalScope()]
        Q_BF, [symscope = GlobalScope()]

        # tissue volume
        V_AQ, [symscope = GlobalScope()]
        V_Cor, [symscope = GlobalScope()]
        V_Lens, [symscope = GlobalScope()]
        V_VH, [symscope = GlobalScope()]

        # permeability surface area 
        PS_cor_pl, [symscope = GlobalScope()]
        PS_lens_pl, [symscope = GlobalScope()]
        PS_ret_pl, [symscope = GlobalScope()]

        # partiton coeffieint (placeholder)
        Kp_ah_pl, [symscope = GlobalScope()]
        Kp_vh_pl, [symscope = GlobalScope()]

    end
    vars = @variables begin
        # PL concentrations in cornea, aqueous humor, vitreous humor, lens
        C_PL_COR(t)
        C_PL_AH(t)
        C_PL_VH(t)
        C_PL_LENS(t)
        # PL concentration in plasma 
        C_PL_Plasma(t)
        # PL concentration in retina interstitium 
        C_PL_IntS_retina(t)
    end

    # derived variable 
    Q_VA = Q_AH + Q_BF - Q_PtA

    
    eqs = [
            Dt(C_PL_VH) ~ (
                  (C_PL_IntS_retina - C_PL_VH) * PS_ret_pl / V_VH
                + (C_PL_LENS - C_PL_VH) * PS_lens_pl / V_VH
                - (C_PL_VH/Kp_vh_pl * Q_VA) / V_VH 
                - (C_PL_VH/Kp_vh_pl * Q_BF) / V_VH
                + ( C_PL_AH/Kp_ah_pl * (Q_BF + Q_VA)  ) / V_VH
            ), 
            Dt(C_PL_COR) ~ (C_PL_AH - C_PL_COR) * PS_cor_pl / V_Cor, 
            Dt(C_PL_AH) ~ (
                  (C_PL_COR - C_PL_AH) * PS_cor_pl / V_AQ
                + (C_PL_LENS - C_PL_AH) * PS_lens_pl / V_AQ
                + (C_PL_VH/Kp_vh_pl * Q_VA) / V_AQ 
                + (C_PL_Plasma * Q_PtA ) / V_AQ
                - ( C_PL_AH/Kp_ah_pl * (Q_AH + Q_BF + Q_VA)  ) / V_AQ
            ), 
            Dt(C_PL_LENS) ~ (
                  (C_PL_AH - C_PL_LENS) * PS_lens_pl
                + (C_PL_VH - C_PL_LENS) * PS_lens_pl
            )/V_Lens
    ]


    System(eqs, t, vars, pars; name=name)

end

#############################################################
#############################################################

##################################################
############ Create organ functions ##############
##################################################

## IgG organ function

function create_organ_igg(organ, exg)
    @independent_variables t
    # dummy variable to track 
    vars = @variables begin
        total_exg_igg(t)
    end
    # create subcompartments objects

    ## All organs except plasma and lymph node (all subcompartments except IgG vascular side)
    vm = VM(; name = :vm, exg = exg)
    e7 = E7(; name = :e7, exg = exg)
    e6a = E6A(; name = :e6a, exg = exg)
    e7b = E7B(; name = :e7b, exg = exg)
    ism = ISM(; name = :ism, exg = exg)
    ints = INTS(; name = :ints)
    bound_vm = BOUND_VM(; name = :bound_vm, exg = exg)
    bound2_vm = BOUND2_VM(; name = :bound2_vm, exg = exg)
    bound_e7 = BOUND_E7(; name = :bound_e7, exg = exg)
    bound2_e7 = BOUND2_E7(; name = :bound2_e7, exg = exg)
    bound_e6a = BOUND_E6A(; name = :bound_e6a, exg = exg)
    bound2_e6a = BOUND2_E6A(; name = :bound2_e6a, exg = exg)
    bound_e7b = BOUND_E7B(; name = :bound_e7b, exg = exg)
    bound2_e7b = BOUND2_E7B(; name = :bound2_e7b, exg = exg)
    bound_ism = BOUND_ISM(; name = :bound_ism, exg = exg)
    bound2_ism = BOUND2_ISM(; name = :bound2_ism, exg = exg)

    if exg
        bound_vm_mem = BOUND_VM_MEM(; name = :bound_vm_mem)
        bound_ism_mem = BOUND_ISM_MEM(; name = :bound_ism_mem)
    end

    ## vascular side subcompartment
    if organ == :liver
        v = V_LIVER(; name = :v)
    elseif organ == :lung
        v = V_LUNG(; name = :v)
    else
        v = V_ALL(; name = :v)
    end

    name_symbol = exg ? Symbol("igg_exg") : Symbol("igg_edg")

    if exg
        ModelingToolkit.extend(System(Equation[];name=name_symbol), [v, vm, bound_vm_mem, e7, e6a, e7b, ism, bound_ism_mem, ints, bound_vm, bound2_vm, bound_e7, bound2_e7, bound_e6a, bound2_e6a, bound_e7b, bound2_e7b, bound_ism, bound2_ism])
    else
        ModelingToolkit.extend(System(Equation[];name=name_symbol), [v, vm, e7, e6a, e7b, ism, ints, bound_vm, bound2_vm, bound_e7, bound2_e7, bound_e6a, bound2_e6a, bound_e7b, bound2_e7b, bound_ism, bound2_ism])
    end
end


## FcRn organ function

function create_organ_fcrn()
    # create subcompartment objects
    @independent_variables t
    Dt = Differential(t)
    fcrn_vm = FCRN_VM(; name = :fcrn_vm)
    fcrn_e7 = FCRN_E7(; name = :fcrn_e7)
    fcrn_e6a = FCRN_E6A(; name = :fcrn_e6a)
    fcrn_e7b = FCRN_E7B(; name = :fcrn_e7b)
    fcrn_ism = FCRN_ISM(; name = :fcrn_ism)

    ModelingToolkit.extend(System(Equation[], t, name=:fcrn), [fcrn_vm, fcrn_e7, fcrn_e6a, fcrn_e7b, fcrn_ism])
end


## join IgG and FcRn
function create_organ(organ)
    @independent_variables t
    ## create organ - IgG
    organ_igg_exg = create_organ_igg(organ, true)
    organ_igg_edg = create_organ_igg(organ, false)

    ## create organ - fcrn
    organ_fcrn = create_organ_fcrn()

    ## create organ PL (added)
    PL_tissue = create_organ_PL()

    connections = [ 

        ## variables ##

        ### igg - fcrn
        #### C_VM
        organ_igg_exg.C_VM ~ organ_fcrn.C_VM_EXG,
        organ_igg_edg.C_VM ~ organ_fcrn.C_VM_EDG,
        #### C_bound_VM
        organ_igg_exg.C_bound_VM ~ organ_fcrn.C_bound_VM_EXG,
        organ_igg_edg.C_bound_VM ~ organ_fcrn.C_bound_VM_EDG,
        #### C_bound2_VM
        organ_igg_exg.C_bound2_VM ~ organ_fcrn.C_bound2_VM_EXG,
        organ_igg_edg.C_bound2_VM ~ organ_fcrn.C_bound2_VM_EDG,
        #### C_E7
        organ_igg_exg.C_E7 ~ organ_fcrn.C_E7_EXG,
        organ_igg_edg.C_E7 ~ organ_fcrn.C_E7_EDG,
        #### C_bound_E7
        organ_igg_exg.C_bound_E7 ~ organ_fcrn.C_bound_E7_EXG,
        organ_igg_edg.C_bound_E7 ~ organ_fcrn.C_bound_E7_EDG,
        ####C_bound2_E7
        organ_igg_exg.C_bound2_E7 ~ organ_fcrn.C_bound2_E7_EXG,
        organ_igg_edg.C_bound2_E7 ~ organ_fcrn.C_bound2_E7_EDG,
        #### C_E6a
        organ_igg_exg.C_E6a ~ organ_fcrn.C_E6a_EXG,
        organ_igg_edg.C_E6a ~ organ_fcrn.C_E6a_EDG,
        #### C_bound_E6a
        organ_igg_exg.C_bound_E6a ~ organ_fcrn.C_bound_E6a_EXG,
        organ_igg_edg.C_bound_E6a ~ organ_fcrn.C_bound_E6a_EDG,
        #### C_bound2_E6a
        organ_igg_exg.C_bound2_E6a ~ organ_fcrn.C_bound2_E6a_EXG,
        organ_igg_edg.C_bound2_E6a ~ organ_fcrn.C_bound2_E6a_EDG,
        #### C_E7b
        organ_igg_exg.C_E7b ~ organ_fcrn.C_E7b_EXG,
        organ_igg_edg.C_E7b ~ organ_fcrn.C_E7b_EDG,
        #### C_bound_E7b
        organ_igg_exg.C_bound_E7b ~ organ_fcrn.C_bound_E7b_EXG,
        organ_igg_edg.C_bound_E7b ~ organ_fcrn.C_bound_E7b_EDG,
        #### C_bound2_E7b
        organ_igg_exg.C_bound2_E7b ~ organ_fcrn.C_bound2_E7b_EXG,
        organ_igg_edg.C_bound2_E7b ~ organ_fcrn.C_bound2_E7b_EDG,
        #### C_ISM
        organ_igg_exg.C_ISM ~ organ_fcrn.C_ISM_EXG,
        organ_igg_edg.C_ISM ~ organ_fcrn.C_ISM_EDG,
        #### C_bound_ISM
        organ_igg_exg.C_bound_ISM ~ organ_fcrn.C_bound_ISM_EXG,
        organ_igg_edg.C_bound_ISM ~ organ_fcrn.C_bound_ISM_EDG,
        #### C_bound2_ISM
        organ_igg_exg.C_bound2_ISM ~ organ_fcrn.C_bound2_ISM_EXG,
        organ_igg_edg.C_bound2_ISM ~ organ_fcrn.C_bound2_ISM_EDG,


        # ### fcrn - igg
        # #### C_FcRn_VM
        organ_fcrn.C_FcRn_VM ~ organ_igg_exg.C_FcRn_VM,
        organ_fcrn.C_FcRn_VM ~ organ_igg_edg.C_FcRn_VM,

        #### C_FcRn_E7
        organ_fcrn.C_FcRn_E7 ~ organ_igg_exg.C_FcRn_E7,
        organ_fcrn.C_FcRn_E7 ~ organ_igg_edg.C_FcRn_E7,

        #### C_FcRn_E6a
        organ_fcrn.C_FcRn_E6a ~ organ_igg_exg.C_FcRn_E6a,
        organ_fcrn.C_FcRn_E6a ~ organ_igg_edg.C_FcRn_E6a,
    
        #### C_FcRn_E7b
        organ_fcrn.C_FcRn_E7b ~ organ_igg_exg.C_FcRn_E7b,
        organ_fcrn.C_FcRn_E7b ~ organ_igg_edg.C_FcRn_E7b,

        #### C_FcRn_ISM
        organ_fcrn.C_FcRn_ISM ~ organ_igg_exg.C_FcRn_ISM,
        organ_fcrn.C_FcRn_ISM ~ organ_igg_edg.C_FcRn_ISM,

        # ### PL 
        PL_tissue.C_E6a ~ organ_igg_exg.C_E6a, 
        PL_tissue.C_bound_VM ~ organ_igg_exg.C_bound_VM, 
        PL_tissue.C_bound2_VM ~ organ_igg_exg.C_bound2_VM, 
        PL_tissue.C_bound_ISM ~ organ_igg_exg.C_bound_ISM, 
        PL_tissue.C_bound2_ISM ~ organ_igg_exg.C_bound2_ISM, 
    ]

    return ModelingToolkit.compose(System(connections, t; name = organ), [organ_igg_exg, organ_igg_edg, organ_fcrn, PL_tissue]; name=organ)
end

##################################################
############ Create eye functions ##############
##################################################

## IgG organ function

function create_icb_igg(exg)    
    v = V_ICB(; name = :v)  # differerent for ICB
    vm = VM(; name = :vm, exg = exg)
    e7 = E7(; name = :e7, exg = exg)
    e6a = E6A(; name = :e6a, exg = exg)
    e7b = E7B(; name = :e7b, exg = exg)
    ism = ISM(; name = :ism, exg = exg)
    ints = INTS_ICB(; name = :ints) # differerent for ICB
    bound_vm = BOUND_VM(; name = :bound_vm, exg = exg)
    bound2_vm = BOUND2_VM(; name = :bound2_vm, exg = exg)
    bound_e7 = BOUND_E7(; name = :bound_e7, exg = exg)
    bound2_e7 = BOUND2_E7(; name = :bound2_e7, exg = exg)
    bound_e6a = BOUND_E6A(; name = :bound_e6a, exg = exg)
    bound2_e6a = BOUND2_E6A(; name = :bound2_e6a, exg = exg)
    bound_e7b = BOUND_E7B(; name = :bound_e7b, exg = exg)
    bound2_e7b = BOUND2_E7B(; name = :bound2_e7b, exg = exg)
    bound_ism = BOUND_ISM(; name = :bound_ism, exg = exg)
    bound2_ism = BOUND2_ISM(; name = :bound2_ism, exg = exg)

    if exg
        bound_vm_mem = BOUND_VM_MEM(; name = :bound_vm_mem)
        bound_ism_mem = BOUND_ISM_MEM(; name = :bound_ism_mem)
    end

    name_symbol = exg ? Symbol("igg_exg") : Symbol("igg_edg")

    if exg
        ModelingToolkit.extend(System(Equation[];name=name_symbol), [v, vm, bound_vm_mem, e7, e6a, e7b, ism, bound_ism_mem, ints, bound_vm, bound2_vm, bound_e7, bound2_e7, bound_e6a, bound2_e6a, bound_e7b, bound2_e7b, bound_ism, bound2_ism])
    else
        ModelingToolkit.extend(System(Equation[];name=name_symbol), [v, vm, e7, e6a, e7b, ism, ints, bound_vm, bound2_vm, bound_e7, bound2_e7, bound_e6a, bound2_e6a, bound_e7b, bound2_e7b, bound_ism, bound2_ism])
    end
end

function create_retina_igg(exg)
    v = V_ALL(; name = :v)  # same for retina
    vm = VM(; name = :vm, exg = exg)
    e7 = E7(; name = :e7, exg = exg)
    e6a = E6A(; name = :e6a, exg = exg)
    e7b = E7B(; name = :e7b, exg = exg)
    ism = ISM(; name = :ism, exg = exg)
    ints = INTS_RETINA(; name = :ints) # differerent for retina 
    bound_vm = BOUND_VM(; name = :bound_vm, exg = exg)
    bound2_vm = BOUND2_VM(; name = :bound2_vm, exg = exg)
    bound_e7 = BOUND_E7(; name = :bound_e7, exg = exg)
    bound2_e7 = BOUND2_E7(; name = :bound2_e7, exg = exg)
    bound_e6a = BOUND_E6A(; name = :bound_e6a, exg = exg)
    bound2_e6a = BOUND2_E6A(; name = :bound2_e6a, exg = exg)
    bound_e7b = BOUND_E7B(; name = :bound_e7b, exg = exg)
    bound2_e7b = BOUND2_E7B(; name = :bound2_e7b, exg = exg)
    bound_ism = BOUND_ISM(; name = :bound_ism, exg = exg)
    bound2_ism = BOUND2_ISM(; name = :bound2_ism, exg = exg)

    if exg
        bound_vm_mem = BOUND_VM_MEM(; name = :bound_vm_mem)
        bound_ism_mem = BOUND_ISM_MEM(; name = :bound_ism_mem)
    end

    name_symbol = exg ? Symbol("igg_exg") : Symbol("igg_edg")

    if exg
        ModelingToolkit.extend(System(Equation[];name=name_symbol), [v, vm, bound_vm_mem, e7, e6a, e7b, ism, bound_ism_mem, ints, bound_vm, bound2_vm, bound_e7, bound2_e7, bound_e6a, bound2_e6a, bound_e7b, bound2_e7b, bound_ism, bound2_ism])
    else
        ModelingToolkit.extend(System(Equation[];name=name_symbol), [v, vm, e7, e6a, e7b, ism, ints, bound_vm, bound2_vm, bound_e7, bound2_e7, bound_e6a, bound2_e6a, bound_e7b, bound2_e7b, bound_ism, bound2_ism])
    end
end

function create_choroid_igg(exg)
    v = V_CHOROID(; name = :v)  # different for choroid
    vm = VM(; name = :vm, exg = exg)
    e7 = E7(; name = :e7, exg = exg)
    e6a = E6A(; name = :e6a, exg = exg)
    e7b = E7B(; name = :e7b, exg = exg)
    ism = ISM(; name = :ism, exg = exg)
    ints = INTS_CHOROID(; name = :ints) # differerent for choroid 
    bound_vm = BOUND_VM(; name = :bound_vm, exg = exg)
    bound2_vm = BOUND2_VM(; name = :bound2_vm, exg = exg)
    bound_e7 = BOUND_E7(; name = :bound_e7, exg = exg)
    bound2_e7 = BOUND2_E7(; name = :bound2_e7, exg = exg)
    bound_e6a = BOUND_E6A(; name = :bound_e6a, exg = exg)
    bound2_e6a = BOUND2_E6A(; name = :bound2_e6a, exg = exg)
    bound_e7b = BOUND_E7B(; name = :bound_e7b, exg = exg)
    bound2_e7b = BOUND2_E7B(; name = :bound2_e7b, exg = exg)
    bound_ism = BOUND_ISM(; name = :bound_ism, exg = exg)
    bound2_ism = BOUND2_ISM(; name = :bound2_ism, exg = exg)

    if exg
        bound_vm_mem = BOUND_VM_MEM(; name = :bound_vm_mem)
        bound_ism_mem = BOUND_ISM_MEM(; name = :bound_ism_mem)
    end

    name_symbol = exg ? Symbol("igg_exg") : Symbol("igg_edg")

    if exg
        ModelingToolkit.extend(System(Equation[];name=name_symbol), [v, vm, bound_vm_mem, e7, e6a, e7b, ism, bound_ism_mem, ints, bound_vm, bound2_vm, bound_e7, bound2_e7, bound_e6a, bound2_e6a, bound_e7b, bound2_e7b, bound_ism, bound2_ism])
    else
        ModelingToolkit.extend(System(Equation[];name=name_symbol), [v, vm, e7, e6a, e7b, ism, ints, bound_vm, bound2_vm, bound_e7, bound2_e7, bound_e6a, bound2_e6a, bound_e7b, bound2_e7b, bound_ism, bound2_ism])
    end
end

function create_eye_parts(organ)
    @independent_variables t
    ## create organ - IgG
    if organ == :icb
        organ_igg_exg = create_icb_igg(true)
        organ_igg_edg = create_icb_igg(false)
    elseif organ == :retina
        organ_igg_exg = create_retina_igg(true)
        organ_igg_edg = create_retina_igg(false)
    else 
        # organ == :choriod
        organ_igg_exg = create_choroid_igg(true)
        organ_igg_edg = create_choroid_igg(false)
    end

    ## create organ - fcrn
    organ_fcrn = create_organ_fcrn()

    ## create organ PL (added)
    PL_tissue = create_organ_PL()

    connections = [ 

        ## variables ##

        ### igg - fcrn
        #### C_VM
        organ_igg_exg.C_VM ~ organ_fcrn.C_VM_EXG,
        organ_igg_edg.C_VM ~ organ_fcrn.C_VM_EDG,
        #### C_bound_VM
        organ_igg_exg.C_bound_VM ~ organ_fcrn.C_bound_VM_EXG,
        organ_igg_edg.C_bound_VM ~ organ_fcrn.C_bound_VM_EDG,
        #### C_bound2_VM
        organ_igg_exg.C_bound2_VM ~ organ_fcrn.C_bound2_VM_EXG,
        organ_igg_edg.C_bound2_VM ~ organ_fcrn.C_bound2_VM_EDG,
        #### C_E7
        organ_igg_exg.C_E7 ~ organ_fcrn.C_E7_EXG,
        organ_igg_edg.C_E7 ~ organ_fcrn.C_E7_EDG,
        #### C_bound_E7
        organ_igg_exg.C_bound_E7 ~ organ_fcrn.C_bound_E7_EXG,
        organ_igg_edg.C_bound_E7 ~ organ_fcrn.C_bound_E7_EDG,
        ####C_bound2_E7
        organ_igg_exg.C_bound2_E7 ~ organ_fcrn.C_bound2_E7_EXG,
        organ_igg_edg.C_bound2_E7 ~ organ_fcrn.C_bound2_E7_EDG,
        #### C_E6a
        organ_igg_exg.C_E6a ~ organ_fcrn.C_E6a_EXG,
        organ_igg_edg.C_E6a ~ organ_fcrn.C_E6a_EDG,
        #### C_bound_E6a
        organ_igg_exg.C_bound_E6a ~ organ_fcrn.C_bound_E6a_EXG,
        organ_igg_edg.C_bound_E6a ~ organ_fcrn.C_bound_E6a_EDG,
        #### C_bound2_E6a
        organ_igg_exg.C_bound2_E6a ~ organ_fcrn.C_bound2_E6a_EXG,
        organ_igg_edg.C_bound2_E6a ~ organ_fcrn.C_bound2_E6a_EDG,
        #### C_E7b
        organ_igg_exg.C_E7b ~ organ_fcrn.C_E7b_EXG,
        organ_igg_edg.C_E7b ~ organ_fcrn.C_E7b_EDG,
        #### C_bound_E7b
        organ_igg_exg.C_bound_E7b ~ organ_fcrn.C_bound_E7b_EXG,
        organ_igg_edg.C_bound_E7b ~ organ_fcrn.C_bound_E7b_EDG,
        #### C_bound2_E7b
        organ_igg_exg.C_bound2_E7b ~ organ_fcrn.C_bound2_E7b_EXG,
        organ_igg_edg.C_bound2_E7b ~ organ_fcrn.C_bound2_E7b_EDG,
        #### C_ISM
        organ_igg_exg.C_ISM ~ organ_fcrn.C_ISM_EXG,
        organ_igg_edg.C_ISM ~ organ_fcrn.C_ISM_EDG,
        #### C_bound_ISM
        organ_igg_exg.C_bound_ISM ~ organ_fcrn.C_bound_ISM_EXG,
        organ_igg_edg.C_bound_ISM ~ organ_fcrn.C_bound_ISM_EDG,
        #### C_bound2_ISM
        organ_igg_exg.C_bound2_ISM ~ organ_fcrn.C_bound2_ISM_EXG,
        organ_igg_edg.C_bound2_ISM ~ organ_fcrn.C_bound2_ISM_EDG,


        # ### fcrn - igg
        # #### C_FcRn_VM
        organ_fcrn.C_FcRn_VM ~ organ_igg_exg.C_FcRn_VM,
        organ_fcrn.C_FcRn_VM ~ organ_igg_edg.C_FcRn_VM,

        #### C_FcRn_E7
        organ_fcrn.C_FcRn_E7 ~ organ_igg_exg.C_FcRn_E7,
        organ_fcrn.C_FcRn_E7 ~ organ_igg_edg.C_FcRn_E7,

        #### C_FcRn_E6a
        organ_fcrn.C_FcRn_E6a ~ organ_igg_exg.C_FcRn_E6a,
        organ_fcrn.C_FcRn_E6a ~ organ_igg_edg.C_FcRn_E6a,
    
        #### C_FcRn_E7b
        organ_fcrn.C_FcRn_E7b ~ organ_igg_exg.C_FcRn_E7b,
        organ_fcrn.C_FcRn_E7b ~ organ_igg_edg.C_FcRn_E7b,

        #### C_FcRn_ISM
        organ_fcrn.C_FcRn_ISM ~ organ_igg_exg.C_FcRn_ISM,
        organ_fcrn.C_FcRn_ISM ~ organ_igg_edg.C_FcRn_ISM,

        # ### PL 
        PL_tissue.C_E6a ~ organ_igg_exg.C_E6a, 
        PL_tissue.C_bound_VM ~ organ_igg_exg.C_bound_VM, 
        PL_tissue.C_bound2_VM ~ organ_igg_exg.C_bound2_VM, 
        PL_tissue.C_bound_ISM ~ organ_igg_exg.C_bound_ISM, 
        PL_tissue.C_bound2_ISM ~ organ_igg_exg.C_bound2_ISM, 

    ]

    return ModelingToolkit.compose(System(connections, t; name = organ), [organ_igg_exg, organ_igg_edg, organ_fcrn, PL_tissue]; name=organ)
end

## 
function create_eye()
    @independent_variables t
    ## create ICB - IgG
        icb = create_eye_parts(:icb)
    ## create retina - IgG
        retina = create_eye_parts(:retina)
    ## create choroid - IgG
        choroid = create_eye_parts(:choroid)
    
    ## create cornea 
        cornea_exg = ModelingToolkit.extend(System(Equation[];name=Symbol("igg_exg")), [CORNEA(; name = :cornea)])
        cornea_edg = ModelingToolkit.extend(System(Equation[];name=Symbol("igg_edg")), [CORNEA(; name = :cornea)])

    ## create lens 
        lens_exg = ModelingToolkit.extend(System(Equation[];name=Symbol("lens_igg_exg")), [LENS(; name = :lens)])
        lens_edg = ModelingToolkit.extend(System(Equation[];name=Symbol("lens_igg_edg")), [LENS(; name = :lens)])

    ## create aqueous humor 
        aqueoushumor_exg = ModelingToolkit.extend(System(Equation[];name=Symbol("ah_igg_exg")), [AQUEOUSHUMOR(; name = :aqueoushumor)])
        aqueoushumor_edg = ModelingToolkit.extend(System(Equation[];name=Symbol("ah_igg_edg")), [AQUEOUSHUMOR(; name = :aqueoushumor)])

    ## create vitreous humor 
        vitreoushumor_exg = ModelingToolkit.extend(System(Equation[];name=Symbol("vh_igg_exg")), [VITREOUSHUMOR(; name = :vitreoushumor)])
        vitreoushumor_edg = ModelingToolkit.extend(System(Equation[];name=Symbol("vh_igg_edg")), [VITREOUSHUMOR(; name = :vitreoushumor)])

    ## create PL in cornea, lens, aqueous humor, vitreous humor
        pl_eye = ModelingToolkit.extend(System(Equation[];name=Symbol("pl_cor_humor")), [PL_eye(; name = :pl_cor_humor)])
        
    connections = [ 

        ## EYE SPECIFIC ##

        #### icb connection to aqueous humor 
        icb.igg_exg.C_AQ ~ aqueoushumor_exg.C_AQ, 
        icb.igg_edg.C_AQ ~ aqueoushumor_edg.C_AQ, 

        #### retina connection to vitreous homor, 
        retina.igg_exg.C_VH ~ vitreoushumor_exg.C_VH, 
        retina.igg_edg.C_VH ~ vitreoushumor_edg.C_VH, 

        #### retina connection to choroid 
        retina.igg_exg.C_IntS_CHO ~ choroid.igg_exg.C_IntS, 
        retina.igg_edg.C_IntS_CHO ~ choroid.igg_edg.C_IntS, 

        #### choroid connection to retina 
        choroid.igg_exg.C_IntS_RET ~ retina.igg_exg.C_IntS, 
        choroid.igg_edg.C_IntS_RET ~ retina.igg_edg.C_IntS, 

        #### cornea to aqueous humor 
        cornea_exg.C_AQ ~ aqueoushumor_exg.C_AQ, 
        cornea_edg.C_AQ ~ aqueoushumor_edg.C_AQ, 

        #### lens to aqueous humor 
        lens_exg.C_AQ ~ aqueoushumor_exg.C_AQ, 
        lens_edg.C_AQ ~ aqueoushumor_edg.C_AQ, 

        #### lens to vitreous humor 
        lens_exg.C_VH ~ vitreoushumor_exg.C_VH, 
        lens_edg.C_VH ~ vitreoushumor_edg.C_VH, 

        #### aqueous humor to ICB 
        aqueoushumor_exg.C_V_ICB ~ icb.igg_exg.C_V, 
        aqueoushumor_edg.C_V_ICB ~ icb.igg_edg.C_V, 

        #### aqueous humor to vitreous humor  
        aqueoushumor_exg.C_VH ~ vitreoushumor_exg.C_VH, 
        aqueoushumor_edg.C_VH ~ vitreoushumor_edg.C_VH, 

        #### aqueous humor to cornea 
        aqueoushumor_exg.C_COR ~ cornea_exg.C_COR, 
        aqueoushumor_edg.C_COR ~ cornea_exg.C_COR,

        #### aqueous humor to lens 
        aqueoushumor_exg.C_LENS ~ lens_exg.C_LENS, 
        aqueoushumor_edg.C_LENS ~ lens_exg.C_LENS,

        #### vitreous humor to retina 
        vitreoushumor_exg.C_IntS_RET ~ retina.igg_exg.C_IntS, 
        vitreoushumor_edg.C_IntS_RET ~ retina.igg_edg.C_IntS, 

        #### vitreous humor to lens 
        vitreoushumor_exg.C_LENS ~ lens_exg.C_LENS, 
        vitreoushumor_edg.C_LENS ~ lens_exg.C_LENS, 
        
        #### vitreous humor to aqueous humor 
        vitreoushumor_exg.C_AQ ~ aqueoushumor_exg.C_AQ, 
        vitreoushumor_edg.C_AQ ~ aqueoushumor_edg.C_AQ, 

        #### PL connections 
        pl_eye.C_PL_IntS_retina ~ retina.PL_tissue.C_PL_IntS, 
        pl_eye.C_PL_Plasma ~ retina.PL_tissue.C_PL_Plasma, 
    ]

    return ModelingToolkit.compose(System(connections, t; name = :eye), [icb, retina, choroid, cornea_exg, cornea_edg, lens_exg, lens_edg, aqueoushumor_exg, aqueoushumor_edg, vitreoushumor_exg, vitreoushumor_edg, pl_eye]; name=:eye)
end



## join IgG-FcRn with plasma and lymph node
function create_pbpk()
    @independent_variables t
    # create plasma and lymph node objects
    plasma_exg = PLASMA(; name = :plasma, exg = true)
    plasma_edg = PLASMA(; name = :plasma, exg = false)
    ln_exg = LN(; name = :ln_exg)
    ln_edg = LN(; name = :ln_edg)

    organs = [:lung, :liver, :heart, :muscle, :skin, :adipose, :bone, :brain, :kidney, :sm_int, :la_int, :pancreas, :thymus, :spleen, :other, :marrow];
    Organs = Dict([organ => create_organ(organ) for organ in organs]);
    Organs[:plasma_exg] = plasma_exg
    Organs[:plasma_edg] = plasma_edg
    Organs[:ln_exg] = ln_exg
    Organs[:ln_edg] = ln_edg

    # create eye 
    Organs[:eye] = create_eye();

    # create plasma PL 
    Organs[:plasma_pl] = PL_plasma(; name = :plasma_pl);

    connections = Equation[]
    # for exg_edg in [:EXG, :EDG]
    connections_exg = [
            # exg - plasma
            Organs[:lung].igg_exg.C_V ~ Organs[:plasma_exg].C_V_Lung,
            Organs[:liver].igg_exg.C_V ~ Organs[:plasma_exg].C_V_Liver,
            Organs[:heart].igg_exg.C_V ~ Organs[:plasma_exg].C_V_Heart,
            Organs[:muscle].igg_exg.C_V ~ Organs[:plasma_exg].C_V_Muscle,
            Organs[:skin].igg_exg.C_V ~ Organs[:plasma_exg].C_V_Skin,
            Organs[:adipose].igg_exg.C_V ~ Organs[:plasma_exg].C_V_Adipose,
            Organs[:bone].igg_exg.C_V ~ Organs[:plasma_exg].C_V_Bone,
            Organs[:brain].igg_exg.C_V ~ Organs[:plasma_exg].C_V_Brain,
            Organs[:kidney].igg_exg.C_V ~ Organs[:plasma_exg].C_V_Kidney,
            Organs[:sm_int].igg_exg.C_V ~ Organs[:plasma_exg].C_V_SI,
            Organs[:la_int].igg_exg.C_V ~ Organs[:plasma_exg].C_V_LI,
            Organs[:pancreas].igg_exg.C_V ~ Organs[:plasma_exg].C_V_Pancreas,
            Organs[:thymus].igg_exg.C_V ~ Organs[:plasma_exg].C_V_Thymus,
            Organs[:spleen].igg_exg.C_V ~ Organs[:plasma_exg].C_V_Spleen,
            Organs[:other].igg_exg.C_V ~ Organs[:plasma_exg].C_V_Other,

            Organs[:marrow].igg_exg.C_V ~ Organs[:plasma_exg].C_V_Marrow,       # added, bone marrow

            Organs[:eye].icb.igg_exg.C_V ~ Organs[:plasma_exg].C_V_ICB,          # added, iris-ciliary body 
            Organs[:eye].retina.igg_exg.C_V ~ Organs[:plasma_exg].C_V_Retina,    # added, retina 
            Organs[:eye].choroid.igg_exg.C_V ~ Organs[:plasma_exg].C_V_Choroid,  # added, choroid

            Organs[:plasma_pl].C_Plasma ~ Organs[:plasma_exg].C_Plasma,         # added nonspecific ADC clearance, ADC -> PL 

            # exg - ln
            Organs[:lung].igg_exg.C_IntS ~ Organs[:ln_exg].C_IntS_Lung,
            Organs[:liver].igg_exg.C_IntS ~ Organs[:ln_exg].C_IntS_Liver,
            Organs[:heart].igg_exg.C_IntS ~ Organs[:ln_exg].C_IntS_Heart,
            Organs[:muscle].igg_exg.C_IntS ~ Organs[:ln_exg].C_IntS_Muscle,
            Organs[:skin].igg_exg.C_IntS ~ Organs[:ln_exg].C_IntS_Skin,
            Organs[:adipose].igg_exg.C_IntS ~ Organs[:ln_exg].C_IntS_Adipose,
            Organs[:bone].igg_exg.C_IntS ~ Organs[:ln_exg].C_IntS_Bone,
            Organs[:brain].igg_exg.C_IntS ~ Organs[:ln_exg].C_IntS_Brain,
            Organs[:kidney].igg_exg.C_IntS ~ Organs[:ln_exg].C_IntS_Kidney,
            Organs[:sm_int].igg_exg.C_IntS ~ Organs[:ln_exg].C_IntS_SI,
            Organs[:la_int].igg_exg.C_IntS ~ Organs[:ln_exg].C_IntS_LI,
            Organs[:pancreas].igg_exg.C_IntS ~ Organs[:ln_exg].C_IntS_Pancreas,
            Organs[:thymus].igg_exg.C_IntS ~ Organs[:ln_exg].C_IntS_Thymus,
            Organs[:spleen].igg_exg.C_IntS ~ Organs[:ln_exg].C_IntS_Spleen,
            Organs[:other].igg_exg.C_IntS ~ Organs[:ln_exg].C_IntS_Other,

            Organs[:marrow].igg_exg.C_IntS ~ Organs[:ln_exg].C_IntS_Marrow,       # added, bone marrow

            Organs[:eye].icb.igg_exg.C_IntS ~ Organs[:ln_exg].C_IntS_ICB,            # added, iris-ciliary body 
            Organs[:eye].retina.igg_exg.C_IntS ~ Organs[:ln_exg].C_IntS_Retina,      # added, retina
            Organs[:eye].choroid.igg_exg.C_IntS ~ Organs[:ln_exg].C_IntS_Choroid,    # added, choroid

            # exg (excluding lung and liver) - lung
            Organs[:heart].igg_exg.C_V_Lung ~ Organs[:lung].igg_exg.C_V,
            Organs[:muscle].igg_exg.C_V_Lung ~ Organs[:lung].igg_exg.C_V,
            Organs[:skin].igg_exg.C_V_Lung ~ Organs[:lung].igg_exg.C_V,
            Organs[:adipose].igg_exg.C_V_Lung ~ Organs[:lung].igg_exg.C_V,
            Organs[:bone].igg_exg.C_V_Lung ~ Organs[:lung].igg_exg.C_V,
            Organs[:brain].igg_exg.C_V_Lung ~ Organs[:lung].igg_exg.C_V,
            Organs[:kidney].igg_exg.C_V_Lung ~ Organs[:lung].igg_exg.C_V,
            Organs[:sm_int].igg_exg.C_V_Lung ~ Organs[:lung].igg_exg.C_V,
            Organs[:la_int].igg_exg.C_V_Lung ~ Organs[:lung].igg_exg.C_V,
            Organs[:pancreas].igg_exg.C_V_Lung ~ Organs[:lung].igg_exg.C_V,
            Organs[:thymus].igg_exg.C_V_Lung ~ Organs[:lung].igg_exg.C_V,
            Organs[:spleen].igg_exg.C_V_Lung ~ Organs[:lung].igg_exg.C_V,
            Organs[:other].igg_exg.C_V_Lung ~ Organs[:lung].igg_exg.C_V,

            Organs[:marrow].igg_exg.C_V_Lung ~ Organs[:lung].igg_exg.C_V,       # added, bone marrow

            Organs[:eye].icb.igg_exg.C_V_Lung ~ Organs[:lung].igg_exg.C_V,          # added, iris-ciliary body 
            Organs[:eye].retina.igg_exg.C_V_Lung ~ Organs[:lung].igg_exg.C_V,       # added, retina 
            Organs[:eye].choroid.igg_exg.C_V_Lung ~ Organs[:lung].igg_exg.C_V,      # added, choroid


            # exg - liver
            Organs[:lung].igg_exg.C_V ~ Organs[:liver].igg_exg.C_V_Lung,
            Organs[:sm_int].igg_exg.C_V ~ Organs[:liver].igg_exg.C_V_SI,
            Organs[:la_int].igg_exg.C_V ~ Organs[:liver].igg_exg.C_V_LI,
            Organs[:pancreas].igg_exg.C_V ~ Organs[:liver].igg_exg.C_V_Pancreas,
            Organs[:spleen].igg_exg.C_V ~ Organs[:liver].igg_exg.C_V_Spleen,

        ]
        push!(connections, connections_exg...)



connections_edg = [
            # edg - plasma
            Organs[:lung].igg_edg.C_V ~ Organs[:plasma_edg].C_V_Lung,
            Organs[:liver].igg_edg.C_V ~ Organs[:plasma_edg].C_V_Liver,
            Organs[:heart].igg_edg.C_V ~ Organs[:plasma_edg].C_V_Heart,
            Organs[:muscle].igg_edg.C_V ~ Organs[:plasma_edg].C_V_Muscle,
            Organs[:skin].igg_edg.C_V ~ Organs[:plasma_edg].C_V_Skin,
            Organs[:adipose].igg_edg.C_V ~ Organs[:plasma_edg].C_V_Adipose,
            Organs[:bone].igg_edg.C_V ~ Organs[:plasma_edg].C_V_Bone,
            Organs[:brain].igg_edg.C_V ~ Organs[:plasma_edg].C_V_Brain,
            Organs[:kidney].igg_edg.C_V ~ Organs[:plasma_edg].C_V_Kidney,
            Organs[:sm_int].igg_edg.C_V ~ Organs[:plasma_edg].C_V_SI,
            Organs[:la_int].igg_edg.C_V ~ Organs[:plasma_edg].C_V_LI,
            Organs[:pancreas].igg_edg.C_V ~ Organs[:plasma_edg].C_V_Pancreas,
            Organs[:thymus].igg_edg.C_V ~ Organs[:plasma_edg].C_V_Thymus,
            Organs[:spleen].igg_edg.C_V ~ Organs[:plasma_edg].C_V_Spleen,
            Organs[:other].igg_edg.C_V ~ Organs[:plasma_edg].C_V_Other,

            Organs[:marrow].igg_edg.C_V ~ Organs[:plasma_edg].C_V_Marrow,       # added, bone marrow

            Organs[:eye].icb.igg_edg.C_V ~ Organs[:plasma_edg].C_V_ICB,          # added, iris-ciliary body 
            Organs[:eye].retina.igg_edg.C_V ~ Organs[:plasma_edg].C_V_Retina,    # added, retina 
            Organs[:eye].choroid.igg_edg.C_V ~ Organs[:plasma_edg].C_V_Choroid,  # added, choroid



            # edg - ln
            Organs[:lung].igg_edg.C_IntS ~ Organs[:ln_edg].C_IntS_Lung,
            Organs[:liver].igg_edg.C_IntS ~ Organs[:ln_edg].C_IntS_Liver,
            Organs[:heart].igg_edg.C_IntS ~ Organs[:ln_edg].C_IntS_Heart,
            Organs[:muscle].igg_edg.C_IntS ~ Organs[:ln_edg].C_IntS_Muscle,
            Organs[:skin].igg_edg.C_IntS ~ Organs[:ln_edg].C_IntS_Skin,
            Organs[:adipose].igg_edg.C_IntS ~ Organs[:ln_edg].C_IntS_Adipose,
            Organs[:bone].igg_edg.C_IntS ~ Organs[:ln_edg].C_IntS_Bone,
            Organs[:brain].igg_edg.C_IntS ~ Organs[:ln_edg].C_IntS_Brain,
            Organs[:kidney].igg_edg.C_IntS ~ Organs[:ln_edg].C_IntS_Kidney,
            Organs[:sm_int].igg_edg.C_IntS ~ Organs[:ln_edg].C_IntS_SI,
            Organs[:la_int].igg_edg.C_IntS ~ Organs[:ln_edg].C_IntS_LI,
            Organs[:pancreas].igg_edg.C_IntS ~ Organs[:ln_edg].C_IntS_Pancreas,
            Organs[:thymus].igg_edg.C_IntS ~ Organs[:ln_edg].C_IntS_Thymus,
            Organs[:spleen].igg_edg.C_IntS ~ Organs[:ln_edg].C_IntS_Spleen,
            Organs[:other].igg_edg.C_IntS ~ Organs[:ln_edg].C_IntS_Other,

            Organs[:marrow].igg_edg.C_IntS ~ Organs[:ln_edg].C_IntS_Marrow,         # added, bone marrow

            Organs[:eye].icb.igg_edg.C_IntS ~ Organs[:ln_edg].C_IntS_ICB,            # added, iris-ciliary body 
            Organs[:eye].retina.igg_edg.C_IntS ~ Organs[:ln_edg].C_IntS_Retina,      # added, retina
            Organs[:eye].choroid.igg_edg.C_IntS ~ Organs[:ln_edg].C_IntS_Choroid,    # added, choroid


            # edg (excluding lung and liver) - lung
            Organs[:heart].igg_edg.C_V_Lung ~ Organs[:lung].igg_edg.C_V,
            Organs[:muscle].igg_edg.C_V_Lung ~ Organs[:lung].igg_edg.C_V,
            Organs[:skin].igg_edg.C_V_Lung ~ Organs[:lung].igg_edg.C_V,
            Organs[:adipose].igg_edg.C_V_Lung ~ Organs[:lung].igg_edg.C_V,
            Organs[:bone].igg_edg.C_V_Lung ~ Organs[:lung].igg_edg.C_V,
            Organs[:brain].igg_edg.C_V_Lung ~ Organs[:lung].igg_edg.C_V,
            Organs[:kidney].igg_edg.C_V_Lung ~ Organs[:lung].igg_edg.C_V,
            Organs[:sm_int].igg_edg.C_V_Lung ~ Organs[:lung].igg_edg.C_V,
            Organs[:la_int].igg_edg.C_V_Lung ~ Organs[:lung].igg_edg.C_V,
            Organs[:pancreas].igg_edg.C_V_Lung ~ Organs[:lung].igg_edg.C_V,
            Organs[:thymus].igg_edg.C_V_Lung ~ Organs[:lung].igg_edg.C_V,
            Organs[:spleen].igg_edg.C_V_Lung ~ Organs[:lung].igg_edg.C_V,
            Organs[:other].igg_edg.C_V_Lung ~ Organs[:lung].igg_edg.C_V,

            Organs[:marrow].igg_edg.C_V_Lung ~ Organs[:lung].igg_edg.C_V,           # added, bone marrow

            Organs[:eye].icb.igg_edg.C_V_Lung ~ Organs[:lung].igg_edg.C_V,          # added, iris-ciliary body 
            Organs[:eye].retina.igg_edg.C_V_Lung ~ Organs[:lung].igg_edg.C_V,       # added, retina 
            Organs[:eye].choroid.igg_edg.C_V_Lung ~ Organs[:lung].igg_edg.C_V,      # added, choroid


            # edg - liver
            Organs[:lung].igg_edg.C_V ~ Organs[:liver].igg_edg.C_V_Lung,
            Organs[:sm_int].igg_edg.C_V ~ Organs[:liver].igg_edg.C_V_SI,
            Organs[:la_int].igg_edg.C_V ~ Organs[:liver].igg_edg.C_V_LI,
            Organs[:pancreas].igg_edg.C_V ~ Organs[:liver].igg_edg.C_V_Pancreas,
            Organs[:spleen].igg_edg.C_V ~ Organs[:liver].igg_edg.C_V_Spleen,
        ]
        push!(connections, connections_edg...)
        # # C_Plasma
        push!(connections, Organs[:lung].igg_exg.C_Plasma ~ Organs[:plasma_exg].C_Plasma)
        push!(connections, Organs[:lung].igg_edg.C_Plasma ~ Organs[:plasma_edg].C_Plasma)

        # # C_LN
        push!(connections, Organs[:ln_exg].C_LN ~ Organs[:plasma_exg].C_LN)
        push!(connections, Organs[:ln_edg].C_LN ~ Organs[:plasma_edg].C_LN)

        # add plasma PL related connections (endosomal concentration)
        push!(connections, Organs[:lung].PL_tissue.C_PL_endo ~ Organs[:plasma_pl].C_PL_endo_Lung)
        push!(connections, Organs[:liver].PL_tissue.C_PL_endo ~ Organs[:plasma_pl].C_PL_endo_Liver)
        push!(connections, Organs[:heart].PL_tissue.C_PL_endo ~ Organs[:plasma_pl].C_PL_endo_Heart)
        push!(connections, Organs[:muscle].PL_tissue.C_PL_endo ~ Organs[:plasma_pl].C_PL_endo_Muscle)
        push!(connections, Organs[:skin].PL_tissue.C_PL_endo ~ Organs[:plasma_pl].C_PL_endo_Skin)
        push!(connections, Organs[:adipose].PL_tissue.C_PL_endo ~ Organs[:plasma_pl].C_PL_endo_Adipose)
        push!(connections, Organs[:bone].PL_tissue.C_PL_endo ~ Organs[:plasma_pl].C_PL_endo_Bone)
        push!(connections, Organs[:brain].PL_tissue.C_PL_endo ~ Organs[:plasma_pl].C_PL_endo_Brain)
        push!(connections, Organs[:kidney].PL_tissue.C_PL_endo ~ Organs[:plasma_pl].C_PL_endo_Kidney)
        push!(connections, Organs[:sm_int].PL_tissue.C_PL_endo ~ Organs[:plasma_pl].C_PL_endo_SI)
        push!(connections, Organs[:la_int].PL_tissue.C_PL_endo ~ Organs[:plasma_pl].C_PL_endo_LI)
        push!(connections, Organs[:pancreas].PL_tissue.C_PL_endo ~ Organs[:plasma_pl].C_PL_endo_Pancreas)
        push!(connections, Organs[:thymus].PL_tissue.C_PL_endo ~ Organs[:plasma_pl].C_PL_endo_Thymus)
        push!(connections, Organs[:spleen].PL_tissue.C_PL_endo ~ Organs[:plasma_pl].C_PL_endo_Spleen)
        push!(connections, Organs[:other].PL_tissue.C_PL_endo ~ Organs[:plasma_pl].C_PL_endo_Other)
        push!(connections, Organs[:marrow].PL_tissue.C_PL_endo ~ Organs[:plasma_pl].C_PL_endo_Marrow)
        push!(connections, Organs[:eye].icb.PL_tissue.C_PL_endo ~ Organs[:plasma_pl].C_PL_endo_ICB)
        push!(connections, Organs[:eye].retina.PL_tissue.C_PL_endo ~ Organs[:plasma_pl].C_PL_endo_Retina)
        push!(connections, Organs[:eye].choroid.PL_tissue.C_PL_endo ~ Organs[:plasma_pl].C_PL_endo_Choroid)

        push!(connections, Organs[:lung].PL_tissue.C_PL_Plasma ~ Organs[:plasma_pl].C_PL_Plasma)
        push!(connections, Organs[:liver].PL_tissue.C_PL_Plasma ~ Organs[:plasma_pl].C_PL_Plasma)
        push!(connections, Organs[:heart].PL_tissue.C_PL_Plasma ~ Organs[:plasma_pl].C_PL_Plasma)
        push!(connections, Organs[:muscle].PL_tissue.C_PL_Plasma ~ Organs[:plasma_pl].C_PL_Plasma)
        push!(connections, Organs[:skin].PL_tissue.C_PL_Plasma ~ Organs[:plasma_pl].C_PL_Plasma)
        push!(connections, Organs[:adipose].PL_tissue.C_PL_Plasma ~ Organs[:plasma_pl].C_PL_Plasma)
        push!(connections, Organs[:bone].PL_tissue.C_PL_Plasma ~ Organs[:plasma_pl].C_PL_Plasma)
        push!(connections, Organs[:brain].PL_tissue.C_PL_Plasma ~ Organs[:plasma_pl].C_PL_Plasma)
        push!(connections, Organs[:kidney].PL_tissue.C_PL_Plasma ~ Organs[:plasma_pl].C_PL_Plasma)
        push!(connections, Organs[:sm_int].PL_tissue.C_PL_Plasma ~ Organs[:plasma_pl].C_PL_Plasma)
        push!(connections, Organs[:la_int].PL_tissue.C_PL_Plasma ~ Organs[:plasma_pl].C_PL_Plasma)
        push!(connections, Organs[:pancreas].PL_tissue.C_PL_Plasma ~ Organs[:plasma_pl].C_PL_Plasma)
        push!(connections, Organs[:thymus].PL_tissue.C_PL_Plasma ~ Organs[:plasma_pl].C_PL_Plasma)
        push!(connections, Organs[:spleen].PL_tissue.C_PL_Plasma ~ Organs[:plasma_pl].C_PL_Plasma)
        push!(connections, Organs[:other].PL_tissue.C_PL_Plasma ~ Organs[:plasma_pl].C_PL_Plasma)
        push!(connections, Organs[:marrow].PL_tissue.C_PL_Plasma ~ Organs[:plasma_pl].C_PL_Plasma)
        push!(connections, Organs[:eye].icb.PL_tissue.C_PL_Plasma ~ Organs[:plasma_pl].C_PL_Plasma)
        push!(connections, Organs[:eye].retina.PL_tissue.C_PL_Plasma ~ Organs[:plasma_pl].C_PL_Plasma)
        push!(connections, Organs[:eye].choroid.PL_tissue.C_PL_Plasma ~ Organs[:plasma_pl].C_PL_Plasma)

        push!(connections, Organs[:plasma_exg].LF_LI ~ Organs[:la_int].LF)
        push!(connections, Organs[:plasma_edg].LF_LI ~ Organs[:la_int].LF)
        push!(connections, Organs[:plasma_exg].LF_Heart ~ Organs[:heart].LF)
        push!(connections, Organs[:plasma_edg].LF_Heart ~ Organs[:heart].LF)
        push!(connections, Organs[:plasma_exg].LF_Kidney ~ Organs[:kidney].LF)
        push!(connections, Organs[:plasma_edg].LF_Kidney ~ Organs[:kidney].LF)
        push!(connections, Organs[:plasma_exg].LF_Muscle ~ Organs[:muscle].LF)
        push!(connections, Organs[:plasma_edg].LF_Muscle ~ Organs[:muscle].LF)
        push!(connections, Organs[:plasma_exg].LF_Skin ~ Organs[:skin].LF)
        push!(connections, Organs[:plasma_edg].LF_Skin ~ Organs[:skin].LF)
        push!(connections, Organs[:plasma_exg].LF_Brain ~ Organs[:brain].LF)
        push!(connections, Organs[:plasma_edg].LF_Brain ~ Organs[:brain].LF)
        push!(connections, Organs[:plasma_exg].LF_Adipose ~ Organs[:adipose].LF)
        push!(connections, Organs[:plasma_edg].LF_Adipose ~ Organs[:adipose].LF)
        push!(connections, Organs[:plasma_exg].LF_Thymus ~ Organs[:thymus].LF)
        push!(connections, Organs[:plasma_edg].LF_Thymus ~ Organs[:thymus].LF)
        push!(connections, Organs[:plasma_exg].LF_Liver ~ Organs[:liver].LF)
        push!(connections, Organs[:plasma_edg].LF_Liver ~ Organs[:liver].LF)
        push!(connections, Organs[:plasma_exg].LF_Spleen ~ Organs[:spleen].LF)
        push!(connections, Organs[:plasma_edg].LF_Spleen ~ Organs[:spleen].LF)
        push!(connections, Organs[:plasma_exg].LF_Pancreas ~ Organs[:pancreas].LF)
        push!(connections, Organs[:plasma_edg].LF_Pancreas ~ Organs[:pancreas].LF)
        push!(connections, Organs[:plasma_exg].LF_SI ~ Organs[:sm_int].LF)
        push!(connections, Organs[:plasma_edg].LF_SI ~ Organs[:sm_int].LF)

        push!(connections, Organs[:plasma_exg].LF_Bone ~ Organs[:bone].LF)
        push!(connections, Organs[:plasma_edg].LF_Bone ~ Organs[:bone].LF)
        push!(connections, Organs[:plasma_exg].LF_Other ~ Organs[:other].LF)
        push!(connections, Organs[:plasma_edg].LF_Other ~ Organs[:other].LF)
        push!(connections, Organs[:plasma_exg].LF_Lung ~ Organs[:lung].LF)
        push!(connections, Organs[:plasma_edg].LF_Lung ~ Organs[:lung].LF)
        # push!(connections, Organs[:plasma_exg].L_LymphNode ~ Organs[:ln_exg].L_LymphNode)
        # push!(connections, Organs[:plasma_edg].L_LymphNode ~ Organs[:ln_edg].L_LymphNode)

        push!(connections, Organs[:plasma_exg].LF_Marrow ~ Organs[:marrow].LF)          # added, bone marrow
        push!(connections, Organs[:plasma_edg].LF_Marrow ~ Organs[:marrow].LF)          # added, bone marrow

        push!(connections, Organs[:plasma_exg].LF_ICB ~ Organs[:eye].icb.LF)            # added, eye connnection, icb
        push!(connections, Organs[:plasma_edg].LF_ICB ~ Organs[:eye].icb.LF)            # added, eye connnection, icb
        push!(connections, Organs[:plasma_exg].LF_Retina ~ Organs[:eye].retina.LF)      # added, eye connnection, retina
        push!(connections, Organs[:plasma_edg].LF_Retina ~ Organs[:eye].retina.LF)      # added, eye connnection, retina
        push!(connections, Organs[:plasma_exg].LF_Choroid ~ Organs[:eye].choroid.LF)    # added, eye connnection, choroid
        push!(connections, Organs[:plasma_edg].LF_Choroid ~ Organs[:eye].choroid.LF)    # added, eye connnection, choroid

        push!(connections, Organs[:plasma_exg].PLQ_Heart ~ Organs[:heart].PLQ)
        push!(connections, Organs[:plasma_edg].PLQ_Heart ~ Organs[:heart].PLQ)
        push!(connections, Organs[:plasma_exg].PLQ_Kidney ~ Organs[:kidney].PLQ)
        push!(connections, Organs[:plasma_edg].PLQ_Kidney ~ Organs[:kidney].PLQ)
        push!(connections, Organs[:plasma_exg].PLQ_Muscle ~ Organs[:muscle].PLQ)
        push!(connections, Organs[:plasma_edg].PLQ_Muscle ~ Organs[:muscle].PLQ)
        push!(connections, Organs[:plasma_exg].PLQ_Skin ~ Organs[:skin].PLQ)
        push!(connections, Organs[:plasma_edg].PLQ_Skin ~ Organs[:skin].PLQ)
        push!(connections, Organs[:plasma_exg].PLQ_Brain ~ Organs[:brain].PLQ)
        push!(connections, Organs[:plasma_edg].PLQ_Brain ~ Organs[:brain].PLQ)
        push!(connections, Organs[:plasma_exg].PLQ_Adipose ~ Organs[:adipose].PLQ)
        push!(connections, Organs[:plasma_edg].PLQ_Adipose ~ Organs[:adipose].PLQ)
        push!(connections, Organs[:plasma_exg].PLQ_Thymus ~ Organs[:thymus].PLQ)
        push!(connections, Organs[:plasma_edg].PLQ_Thymus ~ Organs[:thymus].PLQ)
        push!(connections, Organs[:plasma_exg].PLQ_Liver ~ Organs[:liver].PLQ)
        push!(connections, Organs[:plasma_edg].PLQ_Liver ~ Organs[:liver].PLQ)
        push!(connections, Organs[:plasma_exg].PLQ_Spleen ~ Organs[:spleen].PLQ)
        push!(connections, Organs[:plasma_edg].PLQ_Spleen ~ Organs[:spleen].PLQ)
        push!(connections, Organs[:plasma_exg].PLQ_Pancreas ~ Organs[:pancreas].PLQ)
        push!(connections, Organs[:plasma_edg].PLQ_Pancreas ~ Organs[:pancreas].PLQ)
        push!(connections, Organs[:plasma_exg].PLQ_SI ~ Organs[:sm_int].PLQ)
        push!(connections, Organs[:plasma_edg].PLQ_SI ~ Organs[:sm_int].PLQ)
        push!(connections, Organs[:plasma_exg].PLQ_LI ~ Organs[:la_int].PLQ)
        push!(connections, Organs[:plasma_edg].PLQ_LI ~ Organs[:la_int].PLQ)
        push!(connections, Organs[:plasma_exg].PLQ_Bone ~ Organs[:bone].PLQ)
        push!(connections, Organs[:plasma_edg].PLQ_Bone ~ Organs[:bone].PLQ)
        push!(connections, Organs[:plasma_exg].PLQ_Other ~ Organs[:other].PLQ)
        push!(connections, Organs[:plasma_edg].PLQ_Other ~ Organs[:other].PLQ)
        push!(connections, Organs[:plasma_exg].PLQ_Lung ~ Organs[:lung].PLQ)
        push!(connections, Organs[:plasma_edg].PLQ_Lung ~ Organs[:lung].PLQ)

        push!(connections, Organs[:plasma_exg].PLQ_Marrow ~ Organs[:marrow].PLQ)        # added, bone marrow
        push!(connections, Organs[:plasma_edg].PLQ_Marrow ~ Organs[:marrow].PLQ)        # added, bone marrow


        push!(connections, Organs[:plasma_exg].PLQ_ICB ~ Organs[:eye].icb.PLQ)            # added, eye connnection, icb
        push!(connections, Organs[:plasma_edg].PLQ_ICB ~ Organs[:eye].icb.PLQ)            # added, eye connnection, icb
        push!(connections, Organs[:plasma_exg].PLQ_Retina ~ Organs[:eye].retina.PLQ)      # added, eye connnection, retina
        push!(connections, Organs[:plasma_edg].PLQ_Retina ~ Organs[:eye].retina.PLQ)      # added, eye connnection, retina
        push!(connections, Organs[:plasma_exg].PLQ_Choroid ~ Organs[:eye].choroid.PLQ)    # added, eye connnection, choroid
        push!(connections, Organs[:plasma_edg].PLQ_Choroid ~ Organs[:eye].choroid.PLQ)    # added, eye connnection, choroid


        push!(connections, Organs[:ln_edg].sigma_IS_Heart ~ Organs[:heart].sigma_IS)
        push!(connections, Organs[:ln_exg].sigma_IS_Heart ~ Organs[:heart].sigma_IS)
        push!(connections, Organs[:ln_edg].sigma_IS_Kidney ~ Organs[:kidney].sigma_IS)
        push!(connections, Organs[:ln_exg].sigma_IS_Kidney ~ Organs[:kidney].sigma_IS)
        push!(connections, Organs[:ln_edg].sigma_IS_Muscle ~ Organs[:muscle].sigma_IS)
        push!(connections, Organs[:ln_exg].sigma_IS_Muscle ~ Organs[:muscle].sigma_IS)
        push!(connections, Organs[:ln_edg].sigma_IS_Skin ~ Organs[:skin].sigma_IS)
        push!(connections, Organs[:ln_exg].sigma_IS_Skin ~ Organs[:skin].sigma_IS)
        push!(connections, Organs[:ln_edg].sigma_IS_Brain ~ Organs[:brain].sigma_IS)
        push!(connections, Organs[:ln_exg].sigma_IS_Brain ~ Organs[:brain].sigma_IS)
        push!(connections, Organs[:ln_edg].sigma_IS_Adipose ~ Organs[:adipose].sigma_IS)
        push!(connections, Organs[:ln_exg].sigma_IS_Adipose ~ Organs[:adipose].sigma_IS)
        push!(connections, Organs[:ln_edg].sigma_IS_Thymus ~ Organs[:thymus].sigma_IS)
        push!(connections, Organs[:ln_exg].sigma_IS_Thymus ~ Organs[:thymus].sigma_IS)
        push!(connections, Organs[:ln_edg].sigma_IS_Liver ~ Organs[:liver].sigma_IS)
        push!(connections, Organs[:ln_exg].sigma_IS_Liver ~ Organs[:liver].sigma_IS)
        push!(connections, Organs[:ln_edg].sigma_IS_Spleen ~ Organs[:spleen].sigma_IS)
        push!(connections, Organs[:ln_exg].sigma_IS_Spleen ~ Organs[:spleen].sigma_IS)
        push!(connections, Organs[:ln_edg].sigma_IS_Pancreas ~ Organs[:pancreas].sigma_IS)
        push!(connections, Organs[:ln_exg].sigma_IS_Pancreas ~ Organs[:pancreas].sigma_IS)
        push!(connections, Organs[:ln_edg].sigma_IS_SI ~ Organs[:sm_int].sigma_IS)
        push!(connections, Organs[:ln_exg].sigma_IS_SI ~ Organs[:sm_int].sigma_IS)
        push!(connections, Organs[:ln_edg].sigma_IS_LI ~ Organs[:la_int].sigma_IS)
        push!(connections, Organs[:ln_exg].sigma_IS_LI ~ Organs[:la_int].sigma_IS)
        push!(connections, Organs[:ln_edg].sigma_IS_Bone ~ Organs[:bone].sigma_IS)
        push!(connections, Organs[:ln_exg].sigma_IS_Bone ~ Organs[:bone].sigma_IS)
        push!(connections, Organs[:ln_edg].sigma_IS_Other ~ Organs[:other].sigma_IS)
        push!(connections, Organs[:ln_exg].sigma_IS_Other ~ Organs[:other].sigma_IS)
        push!(connections, Organs[:ln_edg].sigma_IS_Lung ~ Organs[:lung].sigma_IS)
        push!(connections, Organs[:ln_exg].sigma_IS_Lung ~ Organs[:lung].sigma_IS)

        push!(connections, Organs[:ln_edg].sigma_IS_Marrow ~ Organs[:marrow].sigma_IS)          # added, bone marrow
        push!(connections, Organs[:ln_exg].sigma_IS_Marrow ~ Organs[:marrow].sigma_IS)          # added, bone marrow

        push!(connections, Organs[:ln_exg].sigma_IS_ICB ~ Organs[:eye].icb.sigma_IS)            # added, eye connnection, icb
        push!(connections, Organs[:ln_edg].sigma_IS_ICB ~ Organs[:eye].icb.sigma_IS)            # added, eye connnection, icb
        push!(connections, Organs[:ln_exg].sigma_IS_Retina ~ Organs[:eye].retina.sigma_IS)      # added, eye connnection, retina
        push!(connections, Organs[:ln_edg].sigma_IS_Retina ~ Organs[:eye].retina.sigma_IS)      # added, eye connnection, retina
        push!(connections, Organs[:ln_exg].sigma_IS_Choroid ~ Organs[:eye].choroid.sigma_IS)    # added, eye connnection, choroid
        push!(connections, Organs[:ln_edg].sigma_IS_Choroid ~ Organs[:eye].choroid.sigma_IS)    # added, eye connnection, choroid


        
        # added plasma PL related connections (endothelial fraction)
        push!(connections,  Organs[:plasma_pl].Endothelial_Cell_Frac_Lung ~ Organs[:lung].Endothelial_Cell_Frac)
        push!(connections,  Organs[:plasma_pl].Endothelial_Cell_Frac_Liver ~ Organs[:liver].Endothelial_Cell_Frac)
        push!(connections,  Organs[:plasma_pl].Endothelial_Cell_Frac_Heart   ~ Organs[:heart].Endothelial_Cell_Frac)
        push!(connections,  Organs[:plasma_pl].Endothelial_Cell_Frac_Muscle  ~ Organs[:muscle].Endothelial_Cell_Frac)
        push!(connections,  Organs[:plasma_pl].Endothelial_Cell_Frac_Skin    ~ Organs[:skin].Endothelial_Cell_Frac)
        push!(connections,  Organs[:plasma_pl].Endothelial_Cell_Frac_Adipose ~ Organs[:adipose].Endothelial_Cell_Frac)
        push!(connections,  Organs[:plasma_pl].Endothelial_Cell_Frac_Bone    ~ Organs[:bone].Endothelial_Cell_Frac)
        push!(connections,  Organs[:plasma_pl].Endothelial_Cell_Frac_Brain   ~ Organs[:brain].Endothelial_Cell_Frac)
        push!(connections,  Organs[:plasma_pl].Endothelial_Cell_Frac_Kidney  ~ Organs[:kidney].Endothelial_Cell_Frac)
        push!(connections,  Organs[:plasma_pl].Endothelial_Cell_Frac_SI      ~ Organs[:sm_int].Endothelial_Cell_Frac)
        push!(connections,  Organs[:plasma_pl].Endothelial_Cell_Frac_LI      ~ Organs[:la_int].Endothelial_Cell_Frac)
        push!(connections,  Organs[:plasma_pl].Endothelial_Cell_Frac_Pancreas ~ Organs[:pancreas].Endothelial_Cell_Frac)
        push!(connections,  Organs[:plasma_pl].Endothelial_Cell_Frac_Thymus  ~ Organs[:thymus].Endothelial_Cell_Frac)
        push!(connections,  Organs[:plasma_pl].Endothelial_Cell_Frac_Spleen  ~ Organs[:spleen].Endothelial_Cell_Frac)
        push!(connections,  Organs[:plasma_pl].Endothelial_Cell_Frac_Other   ~ Organs[:other].Endothelial_Cell_Frac)
        push!(connections,  Organs[:plasma_pl].Endothelial_Cell_Frac_Marrow  ~ Organs[:marrow].Endothelial_Cell_Frac)
        push!(connections,  Organs[:plasma_pl].Endothelial_Cell_Frac_ICB     ~ Organs[:eye].icb.Endothelial_Cell_Frac)
        push!(connections,  Organs[:plasma_pl].Endothelial_Cell_Frac_Retina  ~ Organs[:eye].retina.Endothelial_Cell_Frac)
        push!(connections,  Organs[:plasma_pl].Endothelial_Cell_Frac_Choroid ~ Organs[:eye].choroid.Endothelial_Cell_Frac)

        # added payload connections 
        push!(connections,  Organs[:eye].pl_cor_humor.C_PL_IntS_retina ~ Organs[:eye].retina.PL_tissue.C_PL_IntS)
        push!(connections,  Organs[:eye].pl_cor_humor.C_PL_Plasma ~ Organs[:plasma_pl].C_PL_Plasma)         


    ModelingToolkit.compose(System(connections, t, name=:pbpk), collect(values(Organs)); name=:pbpk)
end

