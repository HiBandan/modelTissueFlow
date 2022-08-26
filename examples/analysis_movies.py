GenoTypes = { 
                #'SYSTEM': transition-cutoff
                'wild_type': 0.2, # > 0.0
                #'ets': 0.2 # > 0.0
            }

figFormat = 'png'

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
                                
                                'ets':
                                { 
                                   '210616_ets_GAP43mSc_sqhGFP_1_3_--': [['Q-1','D-cw'],[8,[8,45]]], # [['Q-1','D-cw'],[8,[8,45]]]
                                   '210526_ets_GAP43mSc_sqhGFP_1_2': [['Q-1','D-ccw'],[36,[36,75]]], # [['Q-1','D-ccw'],[36,[36,75]]]
                                   '210616_ets_GAP43mSc_sqhGFP_2_1_--': [['Q-1','D-cw'],[17,[17,56]]], # [['Q-1','D-cw'],[17,[17,56]]]
                                   '210617_ets_GAP43mSc_sqhGFP_1_2': [['Q-3','D-cw'],[57,[57,96]]], # [[['Q-3','D-cw'],[57,[57,96]]]
                                   '210726_ets_GAP43mSc_sqhGFP_1_4_--': [['Q-3','D-ccw'],[1,[1,40]]], # [['Q-3','D-ccw'],[1,[1,40]]]
                                   '210604_ets_GAP43mSc_sqhGFP_1_2': [['Q-3','D-ccw'],[18,[18,55]]] # [['Q-3','D-ccw'],[18,[18,55]]],  # reference-for-paper
                                },
                            }
    
            
    
