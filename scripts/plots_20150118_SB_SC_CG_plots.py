import os
base_dir = 'C:\Users\Ben\Documents\GitHub\expression_broad_data'
os.chdir(base_dir) 
from core import expression_plots 
from core import io_library 
#from IPython.core.debugger import Tracer
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt 
import seaborn as sns


def main():
    #Species List: S. Cerevisiae, K. Lactis, S. Bayanus, C. Glabrata, N. Castelli
    #species_list = ['Saccharomyces cerevisiae', 'Kluyveromyces lactis', 'Candida glabrata', 'Saccharomyces castellii' , 'Saccharomyces bayanus']
    species_list = ['Saccharomyces bayanus']
    #filename of combined raw and perturbation expression data
    #fname_out_bases = ['SCer', 'KLac', 'CGla', 'SCas','SBay']
    fname_out_bases = ['SBay']
    
    pct_expression = (0.90,1)
    pct_stability = (0.80,1)
    
    sort_order = ['mean_gene_expression',0]
    #sort_order is either ['stability_score',1] which means sort stability scores in ascending order or
    #['mean_gene_expression',0] which means sort mean_expression scores in descending order
    N_plotted = 20
    
    fig_fname_suffix = '_highexp_stable_by_exp.png'
    
    # Optional: Parse original data file and make CSV for raw expression and perturbations
    # Only need to do this the first time you use the data set.  Subsequent runs can read data directly from the CSV. 
    parse_data = False
    
    if parse_data: 
        io_library.make_data_tables(species_list,fname_out_bases, base_dir)
 
    for jj in range(len(species_list)):
        
        #Load raw expression data
        fname = os.path.normpath(base_dir + "\microarray_data\\raw_exp\\"  + fname_out_bases[jj] + '_raw_exp.csv')
        raw_exp = pd.read_csv(fname, index_col = 'orf_name')
        print fname + ' raw expression dataset loaded'
        #extract mean expression only
        mean_gene_expression = raw_exp['Mean']
        mean_gene_expression.name = 'mean_gene_expression'
        #Run mean extraction histogram plot
        
        
        #Load data for microarrays
        fname = os.path.normpath(base_dir + "\microarray_data\\GSE36253_Growth\\"  + fname_out_bases[jj] + '_growth.csv' )
        growth_exp = pd.read_csv(fname,header = [0,1,2], index_col = [0,1])
        print fname + ' growth microarray dataset loaded'
        
        #group by conditions and take mean
        growth_replicate_groups = growth_exp.groupby(axis = 1, level = 'conditions')
        growth_exp_avg = growth_replicate_groups.aggregate(np.mean)
        
        if species_list[jj] != 'Saccharomyces bayanus':
            #There is no stress dataset for S. bayanus
            fname = os.path.normpath(base_dir + "\microarray_data\\GSE38478_Stress\\"  + fname_out_bases[jj] + '_stress.csv' )
            stress_exp = pd.read_csv(fname,header = [0,1,2], index_col = [0,1])
            print fname + ' stress microarray dataset loaded'
            
            #group by condition and take mean
            stress_replicate_groups = stress_exp.groupby(axis = 1, level = 'conditions')
            stress_exp_avg = stress_replicate_groups.aggregate(np.mean)
            
            #combine growth and stress average expression datasets. 
            if False in stress_exp_avg.index==growth_exp_avg.index:
                print "Error: ID mismatch between condition data. Species = {}".format(species)
                break
            condition_arrays = pd.concat([growth_exp_avg,stress_exp_avg], axis = 1)
        else:
            condition_arrays = growth_exp_avg
        
        #Run Promoter decision plots

  
        Nconds = len(condition_arrays.columns)
        weight_vector = [1 for kk in range(Nconds)]
        #Should I make this a unit vector?
    
        fig_fname =  os.path.normpath("C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\Bioinformatics\Microarrays\\20160215_promoter_Pick_CG_SB_NCAS\\"  + fname_out_bases[jj] + fig_fname_suffix )
        
        fig, ax = expression_plots.promoter_choice_plot(mean_gene_expression, condition_arrays, pct_expression, pct_stability, weight_vector, sort_order = sort_order, N_plotted = N_plotted)
        
        plt.yticks(rotation = 0)
        plt.xticks(rotation = 90) 
        plt.show()
        
        fig.savefig(fig_fname,format = 'png')
        
        print fig_fname + ' plotted and saved'
   
        

if __name__=="__main__":
	main()