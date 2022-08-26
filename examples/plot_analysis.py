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
from modelTissueFlow import inOutTools

########################################
#               PLOTTING               #
########################################
# analysed-data-path
gen_ID_List = list(analysis_movies.GenoTypes.keys())
for gen_counter,gen_ID in enumerate(gen_ID_List): 
    #****************************#
    #           GENOTYPE         #              
    #****************************#
    GENOTYPE = inOutTools.load_data_via_pickle('./inputData/' +'output_' + gen_ID) 
    print('genotype:',GENOTYPE.ID)
    GENOTYPE.view_piv('./movies' +'/'+ gen_ID,analysis_movies.figFormat)
    #****************************#
    #           EMBRYOS          #              
    #****************************#
    EMBRYOS = GENOTYPE.ANIMALS
    for embryo_counter,EMBRYO in enumerate(EMBRYOS[0:1]):
        print('embryo:',EMBRYO.ID)
        EMBRYO.view_frames('./movies' +'/'+gen_ID,inOutTools.read_parameters_from_file('./plot_analysis_parameters.txt'),figFormat='png')
