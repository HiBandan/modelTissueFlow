import time
import analysis_movies
from modelTissueFlow import tissueFlow1D

# start-time
start_time = time.time()
# initialize-tissue-flow
tF = tissueFlow1D()
# loop-over-GENOTYPES
for gen_counter,gen_ID in enumerate(analysis_movies.GenoTypes.keys()):   
    # fitting      
    tF.run_FITTING(expDataPath=analysis_movies.data_path+'/output_'+ gen_ID,parametersFileName='./fit_parameters.txt',hypothesis='curvature',changingParameter='lh',parameters_FIXED={},BOUNDARY_COND_SWITCH = 'PERIODIC',spatial_fitting_domain=[],fit_piv_Type = 'piv',figFormat='pdf')
# end-time    
print("--- %s minutes ---" % ((time.time() - start_time)/60))    
    


