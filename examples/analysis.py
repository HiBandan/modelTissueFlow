#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 13 19:13:34 2022

@author: bandan
"""

#############################################################################################
#                                       1D-tissue-flow                                      #
############################################################################################# 
#import movie_read

import movie_read
from modelTissueFlow import tissueFlow1D,inOutTools

# color-list
colorsList = inOutTools.colors 
colorMaps = inOutTools.colorMaps

########################################
#               ANALYSIS               #
########################################
# locate-movies
inOut_path = movie_read.userPath  if inOutTools.operating_system() == 'Windows' else movie_read.userPath
# initialize-tissue
gen_ID_List = list(movie_read.GenoTypes.keys())
parameters = inOutTools.read_parameters_from_file(inOutTools.current_program_path() + '/analysis_parameters.txt')
tissueFlow_genotype_List = dict(zip(gen_ID_List,[tissueFlow1D(parameters,inOut_path+'/'+gen_ID,movie_read.embryoTypes_frameIndex_Map,movie_read.GenoTypes[gen_ID]) for gen_ID in gen_ID_List]))
#****************************#
#           GENOTYPE         #              
#****************************#
for gen_counter,gen_ID in enumerate(tissueFlow_genotype_List.keys()): 
    print(gen_ID)
    # load-initialized-tissue
    tF = tissueFlow_genotype_List[gen_ID]
    # run-GENOTYPE/embryo-analysis
    GENOTYPE = tF.analysis_SYSTEM(sys_ID=gen_ID) # ANIMAL = tF.analysis_ANIMAL(sys_ID=gen_ID,animal_ID=animal_ID) 
    # save-GENOTYPE-information 
    fileName = tF.inOut_path + '/output_' + gen_ID
    inOutTools.dump_data_via_pickle(fileName,GENOTYPE)
