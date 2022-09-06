import time
import analysis_movies
from modelTissueFlow import inOutTools

# start-time
start_time = time.time()
# loop-over-GENOTYPES
for gen_counter,gen_ID in enumerate(analysis_movies.GenoTypes.keys()): 
    # load-genotype-data
    GENOTYPE = inOutTools.load_data_via_pickle(inPutPath=analysis_movies.data_path+'/output_'+gen_ID) 
    # view-genotype(embryo-averaged)-data
    GENOTYPE.view_piv(figFormat='pdf')
    # loop-over-EMBRYOS
    for embryo_counter,EMBRYO in enumerate(GENOTYPE.ANIMALS):
        # view-individual-embryo-data
        EMBRYO.view_frames(parametersFileName='./plot_parameters.txt',figFormat='pdf')
# end-time    
print("--- %s minutes ---" % ((time.time() - start_time)/60))