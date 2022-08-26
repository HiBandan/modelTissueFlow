#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 13 19:13:34 2022

@author: bandan
"""

#############################################################################################
#                                       1D-tissue-flow                                      #
############################################################################################# 

import time
import analysis_movies
from modelTissueFlow import tissueFlow1D,inOutTools

# color-list
colorsList = inOutTools.colors 
colorMaps = inOutTools.colorMaps

########################################
#               ANALYSIS               #
########################################
start_time = time.time()
# initialize-tissue
gen_ID_List = list(analysis_movies.GenoTypes.keys())
tissueFlow_genotype_List = dict(zip(gen_ID_List,[tissueFlow1D(inOutTools.read_parameters_from_file('./analysis_parameters.txt'),'.',analysis_movies.embryoTypes_frameIndex_Map,analysis_movies.GenoTypes[gen_ID]) for gen_ID in gen_ID_List]))
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
    fileName = './inputData/output_' + gen_ID
    inOutTools.dump_data_via_pickle(fileName,GENOTYPE)
print("--- %s minutes ---" % ((time.time() - start_time)/60))