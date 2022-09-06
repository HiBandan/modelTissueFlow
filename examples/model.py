import time
import analysis_movies
from modelTissueFlow import tissueFlow1D

# start-time
start_time = time.time()
# initialize-tissue-flow
tF = tissueFlow1D()
# modelling
tF.run_MODEL(outputPath=analysis_movies.movie_path,parametersFileName='./modelling_parameters.txt',figFormat='pdf')
# end-time    
print("--- %s minutes ---" % ((time.time() - start_time)/60))    
    


