import time
import analysis_movies
from modelTissueFlow import tissueFlow1D

# start-time
start_time = time.time()
# initialize-tissue
initialized_GENOTYPE_List = dict(zip(list(analysis_movies.GenoTypes.keys()),[tissueFlow1D('./analysis_parameters.txt',analysis_movies.movie_path,analysis_movies.embryoTypes_frameIndex_Map,analysis_movies.GenoTypes[gen_ID]) for gen_ID in analysis_movies.GenoTypes.keys()]))
# loop-over-GENOTYPES
for gen_ID,GENOTYPE in initialized_GENOTYPE_List.items(): 
    # specific-genotype
    GENOTYPE.analysis_SYSTEM(sys_ID=gen_ID,outPutDir=analysis_movies.data_path)
# end-time    
print("--- %s minutes ---" % ((time.time() - start_time)/60))