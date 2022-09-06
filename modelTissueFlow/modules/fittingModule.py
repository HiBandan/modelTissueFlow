import numpy
from scipy.optimize import curve_fit

# user-import
from modelTissueFlow.modules import inOutTools
from modelTissueFlow.modules import equationModule

def curve_fitting_and_prediction(INPUT,fittingDomain,BOUNDARY_COND_SWITCH,input_parameter,tension_flag,curvature_flag,v_grad_correction_flag,distinguish_myo_switch,fix_PARAMETERS,fit_piv_avg): 
    #################
    # specify-input #
    #################
    s_res,s,piv,total_myoGrad,apical_myoGrad,basal_myoGrad,total_mom_curvGrad,apical_mom_curvGrad,basal_mom_curvGrad,v_grad_correction,friction,viscosity = INPUT   
    myo_dependent_inputs = [tension_flag*apical_myoGrad,tension_flag*basal_myoGrad,curvature_flag*apical_mom_curvGrad,curvature_flag*basal_mom_curvGrad]     
    # distinguish-apical-basal_myo
    apical_myoGrad,basal_myoGrad,apical_mom_curvGrad,basal_mom_curvGrad = distinguish_apical_basal_myo(distinguish_myo_switch,myo_dependent_inputs)
    # turn-on/off-v-grad-correction
    v_grad_correction = -1.0*v_grad_correction_flag*v_grad_correction # NOTE: flip-in-sign 
    ##################
    # all-parameters #
    ##################
    all_PARAMETER_indices = [item for item,_ in enumerate(input_parameter)]
    numAllParameters = len(all_PARAMETER_indices)
    ##################
    # fix-parameters #
    ##################
    fix_PARAMETER_indices = [key for key,_ in fix_PARAMETERS.items()]
    ###################
    # free-parameters #
    ###################
    free_PARAMETER_indices = numpy.sort(list(set(all_PARAMETER_indices) - set(fix_PARAMETER_indices)))
    if free_PARAMETER_indices.size == 0:
        free_PARAMETER_indices = fix_PARAMETER_indices
    free_input_parameter = [input_parameter[indx] for indx in free_PARAMETER_indices]
    free_PARAMETERS_len = [len(p) for p in free_input_parameter] 
    numFreeParameters =len(free_PARAMETERS_len)
    shift_free_PARAMETERS_len = numpy.cumsum(free_PARAMETERS_len)
    ###############
    # constraints #
    ###############
    const_HF = [numpy.PZERO,numpy.PINF]
    const_HV = [numpy.PZERO,numpy.PINF]
    const_l_H = [numpy.PZERO,numpy.PINF]
    const_r_a = [numpy.PZERO,numpy.PINF] 
    const_r_b = [numpy.PZERO,numpy.PINF]
    const_eta = [numpy.PZERO,numpy.PINF] 
    const_gamma = [numpy.PZERO,numpy.PINF] 
    const_param = [const_HF,const_HV,const_l_H,const_r_a,const_r_b,const_gamma,const_eta] 
    const_matrix = numpy.array([p for sub_parameter in [[const_param[indx]] if len(free_input_parameter[count]) == 1 else [const_param[indx] for _ in range(s.shape[0])] for count,indx in enumerate(free_PARAMETER_indices)] for p in sub_parameter]) 
    #################################
    # parameter-spilitting-function #
    #################################
    def split_PARAMETERS_ARRAY_to_PARAMETERS(time_indx,PARAMETERS_ARRAY):
        PARAMETERS = [0]*numAllParameters
        # construct-free-parameters
        free_PARAMETERS = [PARAMETERS_ARRAY[0] if free_PARAMETERS_len[0]==1 else PARAMETERS_ARRAY[0+time_indx]] 
        for shift_count,item in enumerate(shift_free_PARAMETERS_len[:-1]):
            free_PARAMETERS.append(PARAMETERS_ARRAY[item] if free_PARAMETERS_len[shift_count+1]==1 else PARAMETERS_ARRAY[item+time_indx])   
        # update-with-free-parameters
        for free_PARAMETER_counter,item in enumerate(free_PARAMETER_indices):
            PARAMETERS[item] = free_PARAMETERS[free_PARAMETER_counter] 
        # update-with-fix-parameters
        for fixed_key,item in fix_PARAMETERS.items():
            PARAMETERS[fixed_key] = item[time_indx] if isinstance(item,list) else item
        return PARAMETERS
    ##########################
    # curve-fitting-function #
    ##########################
    def curve_fitting_function(apical_myoGrad,basal_myoGrad,apical_mom_curvGrad,basal_mom_curvGrad,v_grad_correction,friction,viscosity,fit_piv_avg):
        def solve_model_equation(arc_len,*PARAMETERS_ARRAY):
            v_sol = []
            PARAMETERS_ARRAY = list(PARAMETERS_ARRAY)
            for time_indx,(v_exp,t_g_apical,t_g_basal,m_c_g_apical,m_c_g_basal,v_g,f,e) in enumerate(zip(piv,apical_myoGrad,basal_myoGrad,apical_mom_curvGrad,basal_mom_curvGrad,v_grad_correction,friction,viscosity)):   
                PARAMETERS = split_PARAMETERS_ARRAY_to_PARAMETERS(time_indx,PARAMETERS_ARRAY)
                BOUNDARY_CONDITIONS = equationModule.boundary_conditions(PARAMETERS,v_exp) 
                s = arc_len[time_indx]
                # solve-equation
                v,_ = equationModule.ODE_solver_with_non_constant_coefficients(s,f,e,m_c_g_apical,m_c_g_basal,t_g_apical,t_g_basal,v_g,BOUNDARY_CONDITIONS[BOUNDARY_COND_SWITCH],PARAMETERS)
                # select-data-only-for-fitting-domain
                s_for_fitting = s[fittingDomain]
                v_for_fitting = v[fittingDomain]
                v_sol.append(inOutTools.area_under_curve(s_for_fitting,v_for_fitting,closed=False)/(s_for_fitting[-1]-s_for_fitting[0])) if fit_piv_avg == 'piv_avg' else v_sol.append(v_for_fitting)  
            return numpy.hstack(v_sol)
        return solve_model_equation 
    ###########################
    # individual/collective ? #
    ###########################
    if len(s.shape) == 1: 
        s_res,s,piv,total_myoGrad,apical_myoGrad,basal_myoGrad,total_mom_curvGrad,apical_mom_curvGrad,basal_mom_curvGrad,v_grad_correction,friction,viscosity = [numpy.array([item]) for item in [s_res,s,piv,total_myoGrad,apical_myoGrad,basal_myoGrad,total_mom_curvGrad,apical_mom_curvGrad,basal_mom_curvGrad,v_grad_correction,friction,viscosity]]
    ###########
    # fitting #
    ###########
    est_PARAMETERS = [p for sub_parameter in free_input_parameter for p in sub_parameter]
    piv_for_fitting = numpy.hstack([inOutTools.area_under_curve(a[fittingDomain],b[fittingDomain],closed=False)/(a[fittingDomain][-1]-a[fittingDomain][0]) for a,b in zip(s,piv)]) if fit_piv_avg == 'piv_avg' else numpy.hstack([p[fittingDomain] for p in piv])
    est_PARAMETERS,_ = curve_fit(curve_fitting_function(apical_myoGrad,basal_myoGrad,apical_mom_curvGrad,basal_mom_curvGrad,v_grad_correction,friction,viscosity,fit_piv_avg),s,piv_for_fitting,p0=est_PARAMETERS,bounds=(const_matrix[:,0],const_matrix[:,1])) 
    ##############
    # prediction #
    ##############
    v_fit = []
    chi_squr =[]
    PARAMETERS_ARRAY = list(est_PARAMETERS)  
    for time_indx,(v_exp,t_g_apical,t_g_basal,m_c_g_apical,m_c_g_basal,v_g,f,e) in enumerate(zip(piv,apical_myoGrad,basal_myoGrad,apical_mom_curvGrad,basal_mom_curvGrad,v_grad_correction,friction,viscosity)):   
        PARAMETERS = split_PARAMETERS_ARRAY_to_PARAMETERS(time_indx,PARAMETERS_ARRAY)
        BOUNDARY_CONDITIONS = equationModule.boundary_conditions(PARAMETERS,v_exp)
        v_sol,_ = equationModule.ODE_solver_with_non_constant_coefficients(s[time_indx],f,e,m_c_g_apical,m_c_g_basal,t_g_apical,t_g_basal,v_g,BOUNDARY_CONDITIONS[BOUNDARY_COND_SWITCH],PARAMETERS)
        v_fit.append(v_sol[fittingDomain]) 
        chi_squr.append(numpy.sum(numpy.square(numpy.array(v_sol[fittingDomain])-numpy.array(v_exp[fittingDomain]))))
    s_fit = [s[fittingDomain] for s in s_res]
    ##########################
    # restructure-parameters #
    ##########################
    # free-parameters 
    free_input_parameter = [PARAMETERS_ARRAY[:shift_free_PARAMETERS_len[0]]]
    for free_PARAMETER_count in range(numFreeParameters-1):
        free_input_parameter.append(PARAMETERS_ARRAY[shift_free_PARAMETERS_len[free_PARAMETER_count]:shift_free_PARAMETERS_len[free_PARAMETER_count+1]]) 
    for free_PARAMETER_count,free_PARAMETER_indx in enumerate(free_PARAMETER_indices):
        input_parameter[free_PARAMETER_indx] = free_input_parameter[free_PARAMETER_count]  
    # fix-parameters 
    for fixed_key,item in fix_PARAMETERS.items():
        input_parameter[fixed_key] = item if isinstance(item,list) else [item]
    return input_parameter,v_fit if len(v_fit) > 1 else v_fit[0],chi_squr,s_fit if len(s_fit) > 1 else s_fit[0]

def distinguish_apical_basal_myo(distinguish_myo_switch,myo_dependent_inputs):
    apical_myoGrad,basal_myoGrad,apical_mom_curvGrad,basal_mom_curvGrad = myo_dependent_inputs
    ############
    # ra != rb #
    ############
    if distinguish_myo_switch:
        # apical: myo/mom-curv-grad
        apical_myoGrad = -1.0*apical_myoGrad # NOTE: flip-in-sign
        apical_mom_curvGrad = 1.0*apical_mom_curvGrad
        # basal: myo/mom-curv-grad
        basal_myoGrad = -1.0*basal_myoGrad # NOTE: flip-in-sign
        basal_mom_curvGrad = -1.0*basal_mom_curvGrad # NOTE: flip-in-sign
    ######################
    # r_eff = ra, rb = 0 #
    ######################
    else:
        # apical: myo/mom-curv-grad
        apical_myoGrad = -1.0*apical_myoGrad -1.0*basal_myoGrad 
        apical_mom_curvGrad = 1.0*apical_mom_curvGrad -1.0*basal_mom_curvGrad
        # basal: myo/mom-curv-grad
        basal_myoGrad = numpy.zeros_like(basal_myoGrad) # NOTE: flip-in-sign
        basal_mom_curvGrad = numpy.zeros_like(apical_mom_curvGrad) # NOTE: flip-in-sign
    return apical_myoGrad,basal_myoGrad,apical_mom_curvGrad,basal_mom_curvGrad

def prediction_by_equation(s,INPUT,friction,viscosity,BOUNDARY_COND_SWITCH,PARAMETERS,tension_flag,curvature_flag,v_grad_correction_flag,distinguish_myo_switch): 
    # extract-input 
    apical_myoGrad,basal_myoGrad,apical_mom_curvGrad,basal_mom_curvGrad,v_grad_correction = INPUT
    myo_dependent_inputs = [tension_flag*apical_myoGrad,tension_flag*basal_myoGrad,curvature_flag*apical_mom_curvGrad,curvature_flag*basal_mom_curvGrad]     
    # distinguish-apical-basal_myo
    apical_myoGrad,basal_myoGrad,apical_mom_curvGrad,basal_mom_curvGrad = distinguish_apical_basal_myo(distinguish_myo_switch,myo_dependent_inputs)    
    # turn-on/off-v-grad-correction    
    v_grad_correction = -1.0*v_grad_correction_flag*v_grad_correction # NOTE: flip-in-sign   
    # boundary-condition
    PARAMETERS = list(PARAMETERS)
    BOUNDARY_CONDITIONS = equationModule.boundary_conditions(PARAMETERS,[1e-6,1e-6])
    # solution
    v_sol,fric = equationModule.ODE_solver_with_non_constant_coefficients(s,friction,viscosity,apical_mom_curvGrad,basal_mom_curvGrad,apical_myoGrad,basal_myoGrad,v_grad_correction,BOUNDARY_CONDITIONS[BOUNDARY_COND_SWITCH],PARAMETERS)
    return v_sol,fric
   
