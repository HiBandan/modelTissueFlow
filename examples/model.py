import time
import analysis_movies
from modelTissueFlow import Genotype1D

# start-time
start_time = time.time()
# initialize-tissue-flow
GENOTYPE = Genotype1D()
# modelling
GENOTYPE.run_MODEL(outputPath=analysis_movies.movie_path,parametersFileName='./modelling_parameters.txt',figFormat='png')
# end-time    
print("--- %s minutes ---" % ((time.time() - start_time)/60))    
    


