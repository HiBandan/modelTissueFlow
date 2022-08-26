#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 13 19:13:34 2022

@author: bandan
"""

#############################################################################################
#                                       1D-tissue-flow                                      #
############################################################################################# 

from modelTissueFlow import tissueFlow1D,inOutTools

#******************************************#
#           INITIALIZE-TISSUE-FLOW         #              
#******************************************#
tF = tissueFlow1D(parameters=inOutTools.read_parameters_from_file('./analysis_parameters.txt'),path='.')
#******************************#
#           MODELLING          #              
#******************************#
tF.run_MODEL(inOutTools.read_parameters_from_file('./model_parameters.txt'),figFormat='png')


