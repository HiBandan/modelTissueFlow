import numpy

# user-import
from modelTissueFlow.modules import inOutTools

def boundary_conditions(PARAMETERS,piv):
    boundary_conditions  = {'PERIODIC' : [],'LEFT_RIGHT_fixed' : [0.0,0.0],'LEFT_RIGHT_experimental_input' : [piv[0],piv[-1]]}
    return boundary_conditions

def ODE_solver_with_non_constant_coefficients(s,friction,viscosity,apical_mom_curvGrad,basal_mom_curvGrad,apical_myoGrad,basal_myoGrad,v_grad_correction,BOUNDARY_COND,PARAMETERS):
    ###################
    # equation-domain #
    ###################
    N = len(s)
    ds = numpy.average(s[1:] - s[:-1])
    s = numpy.insert(s,[0,len(s)],[s[0],s[-1]],axis=0)
    ##############
    # parameters #
    ##############
    gamma_Delta = PARAMETERS[0]
    eta_Delta = PARAMETERS[1]
    inv_sqr_L_H = PARAMETERS[2]
    r_a = PARAMETERS[3]
    r_b = PARAMETERS[4]
    a_gamma = PARAMETERS[5]
    a_eta = PARAMETERS[6]
    # rescaling-parameters
    inv_sqr_L_H = inv_sqr_L_H*(a_gamma/a_eta)
    r_a = r_a/a_eta
    r_b = r_b/a_eta
    #########################
    # equation-coefficients #
    #########################
    # viscosity
    A = 1.0 + viscosity*eta_Delta
    A = numpy.insert(A,[0,len(A)],[A[0],A[-1]],axis=0)
    # nothing
    B = numpy.ones_like(s)*0.0
    B = numpy.insert(B,[0,len(B)],[B[0],B[-1]],axis=0)
    # friction
    friction_0 = 1.0 + friction*gamma_Delta
    C = -1.0*friction_0*inv_sqr_L_H
    C = numpy.insert(C,[0,len(C)],[C[0],C[-1]],axis=0)
    # active-tension/curvature
    D = r_a*apical_mom_curvGrad + r_a*apical_myoGrad + r_b*basal_mom_curvGrad + r_b*basal_myoGrad + v_grad_correction  
    D = numpy.insert(D,[0,len(D)],[D[0],D[-1]],axis=0)
    ###################
    # equation-output #
    ###################
    y = None
    # coefficient-matrix 
    coeff_Mat = numpy.zeros((N + 2,N + 2))
    for indx in range(1,N + 1):
        coeff_Mat[indx][indx-1]  = A[indx]/(ds*ds) - B[indx]/(2.0*ds)
        coeff_Mat[indx][indx]  = C[indx] - 2.0*A[indx]/(ds*ds)
        coeff_Mat[indx][indx+1]  = A[indx]/(ds*ds) + B[indx]/(2.0*ds)
    inhomo_Vec = numpy.zeros(N + 2)
    for indx in range(1,N + 1):
        inhomo_Vec[indx] = D[indx]
    #########################
    # non-periodic-boundary #
    #########################
    if BOUNDARY_COND: 
        inhomo_Vec[2] = inhomo_Vec[2] - coeff_Mat[2][1]*BOUNDARY_COND[0]
        inhomo_Vec[N-1] = inhomo_Vec[N-1] - coeff_Mat[N-1][N]*BOUNDARY_COND[1]
        coeff_Mat = numpy.delete(coeff_Mat, [0,1,N,N+1], 0)
        coeff_Mat = numpy.delete(coeff_Mat, [0,1,N,N+1], 1)
        inhomo_Vec = numpy.delete(inhomo_Vec, [0,1,N,N+1], 0)
        # solve-equation
        y = numpy.linalg.solve(numpy.array(coeff_Mat),numpy.array(inhomo_Vec))
        inOutTools.error(not numpy.allclose(numpy.dot(coeff_Mat,y),inhomo_Vec),inOutTools.get_linenumber(),'equationsModule -> no solution of the equation !')
        y = numpy.insert(y,[0,len(y)],[BOUNDARY_COND[0],BOUNDARY_COND[1]],axis=0)# virtual-node-entry
    #####################
    # periodic-boundary #
    #####################
    else:
        coeff_Mat[2][N] = coeff_Mat[2][1]
        coeff_Mat[N][2] = coeff_Mat[N][N+1]*(1 + coeff_Mat[1][2]/coeff_Mat[1][0])
        coeff_Mat[N][N-1] = coeff_Mat[N][N-1] + coeff_Mat[N][N+1]
        coeff_Mat[N][N] = coeff_Mat[N][N] + (coeff_Mat[N][N+1]*coeff_Mat[1][1])/coeff_Mat[1][0]
        inhomo_Vec[N] = inhomo_Vec[N] + (coeff_Mat[N][N+1]*inhomo_Vec[1])/coeff_Mat[1][0]
        coeff_Mat = numpy.delete(coeff_Mat, [0,1,N+1], 0)
        coeff_Mat = numpy.delete(coeff_Mat, [0,1,N+1], 1)
        inhomo_Vec = numpy.delete(inhomo_Vec, [0,1,N+1], 0)
        # solve-equation
        y = numpy.linalg.solve(numpy.array(coeff_Mat),numpy.array(inhomo_Vec))
        inOutTools.error(not numpy.allclose(numpy.dot(coeff_Mat,y),inhomo_Vec),inOutTools.get_linenumber(),'equationsModule -> no solution of the equation !')
        y = numpy.insert(y,0,y[-1],axis=0)
    return y,friction_0    
