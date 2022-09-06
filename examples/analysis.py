import time
import analysis_movies
from modelTissueFlow import Genotype1D

# start-time
start_time = time.time()
# loop-over-GENOTYPES
for gen_ID,transition_cutOff_val in analysis_movies.GenoTypes.items():
    # initialize
    GENOTYPE = Genotype1D('./analysis_parameters.txt',gen_ID,analysis_movies.embryoTypes_frameIndex_Map[gen_ID],transition_cutOff_val,analysis_movies.movie_path)
    # analysis
    GENOTYPE.analysis_GENOTYPE(outPutDir=analysis_movies.data_path)
# end-time    
print("--- %s minutes ---" % ((time.time() - start_time)/60))    
    
