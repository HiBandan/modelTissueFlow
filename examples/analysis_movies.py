#*******************#
# system-to-analyze #
#*******************#
GenoTypes = { 
                #'SYSTEM': transition-cutoff
                'wild_type': 0.2, # > 0.0
                #'ets': 0.2 # > 0.0
            }

movie_path = './movies'
data_path = './data'

#####################################################################
#                               all-systems                         #
#####################################################################
embryoTypes_frameIndex_Map = {    
                                #####################
                                # enrty-information #
                                #####################
                                # 'SYSTEM':
                                # {
                                     # 'ANIMAL' :  [[A,B],[C,[D,E]]]  
                                     ## A: location-of-posterior-pole, which-quadrant-(Q) ?             
                                     ## B: orientation-of-epithelium, clock-wise(cw)-or-anti-clock-wise(ccw) ?     
                                     ## C: intial-frame, cellular-front-passing-nucleous     
                                     ## D: first-segmented-frame
                                     ## E: last-segmented-frame  
                                # }
                                
                                  
                                # wild-type 
                                'wild_type':
                                { 
                                    'movie_1': [['Q-2','D-cw'],[44,[44,53]]],
                                    'movie_2': [['Q-3','D-cw'],[47,[47,56]]]
                                },     
                            }
    
            
    
