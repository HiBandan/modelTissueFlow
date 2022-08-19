GenoTypes = { 
                #'SYSTEM': transition-cutoff
                'wild_type': 0.2 # > 0.0
            }

figFormat = 'png'

userPath = './movies' 

# embryo-types
embryoTypes_frameIndex_Map = {    
                                #******************************************************************************************#
                                # 'SYSTEM':
                                # {
                                     # 'ANIMAL' :  [[A,B],[C,[D,E]]]  
                                     ## A: location-of-posterior-pole, which-quadrant-(Q) ?             
                                     ## B: orientation-of-epithelium, clock-wise(cw)-or-anti-clock-wise(ccw) ?     
                                     ## C: intial-frame, cellular-front-passing-nucleous     
                                     ## D: first-segmented-frame
                                     ## E: last-segmented-frame  
                                # }
                                #******************************************************************************************#
                                  
                                ######
                                # WT #
                                ######
                                'wild_type':
                                { 
                                    'movie_1': [['Q-2','D-cw'],[44,[44,53]]],
                                    'movie_2': [['Q-3','D-cw'],[47,[47,56]]]
                                },     
                            }
    
            
    
