#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 13 19:13:34 2022

@author: bandan
"""

#############################################################################################
#                                       1D-tissue-flow                                      #
############################################################################################# 
import analysis_movies
from modelTissueFlow import tissueFlow1D,inOutTools

#******************************************#
#           INITIALIZE-TISSUE-FLOW         #              
#******************************************#
tF = tissueFlow1D(parameters=inOutTools.read_parameters_from_file('./analysis_parameters.txt'),path='.')
gen_ID_List = list(analysis_movies.GenoTypes.keys())
for gen_counter,gen_ID in enumerate(gen_ID_List): 
    #****************************#
    #           GENOTYPE         #              
    #****************************#
    GENOTYPE = inOutTools.load_data_via_pickle('./inputData/' +'output_' + gen_ID)   
    print('genotype:',GENOTYPE.ID)
    
    print(GENOTYPE.time_List)
    #****************************#
    #           FITTING          #              
    #****************************#
    tF.run_FITTING(GENOTYPE,inOutTools.read_parameters_from_file('./fitting_parameters.txt'),hypothesis='curvature',changingParameter='lh',parameters_FIXED={},BOUNDARY_COND_SWITCH = 'PERIODIC',spatial_fitting_domain=[],fit_piv_Type = 'piv',figFormat='png')
    
    


