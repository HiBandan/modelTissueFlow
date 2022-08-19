#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 13 19:13:34 2022

@author: bandan
"""

#############################################################################################
#                                       1D-tissue-flow                                      #
############################################################################################# 
import movie_read
from modelTissueFlow import inOutTools

########################################
#               PLOTTING               #
########################################
# analysed-data-path
inOut_path = movie_read.userPath  if inOutTools.operating_system() == 'Windows' else movie_read.userPath
gen_ID_List = list(movie_read.GenoTypes.keys())
parameters = inOutTools.read_parameters_from_file(inOutTools.current_program_path() + '/plot_parameters.txt')
for gen_counter,gen_ID in enumerate(gen_ID_List): 
    #****************************#
    #           GENOTYPE         #              
    #****************************#
    GENOTYPE = inOutTools.load_data_via_pickle(inOut_path +'/'+ gen_ID +'/output_' + gen_ID) 
    print('genotype:',GENOTYPE.ID)
    GENOTYPE.view_piv(inOut_path +'/'+ gen_ID,movie_read.figFormat)
    #****************************#
    #           EMBRYOS          #              
    #****************************#
    EMBRYOS = GENOTYPE.ANIMALS
    for embryo_counter,EMBRYO in enumerate(EMBRYOS[:]):
        print('embryo:',EMBRYO.ID)
        EMBRYO.view_frames(inOut_path +'/'+gen_ID,parameters,figFormat='png')
