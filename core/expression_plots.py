import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt 
import matplotlib as mpl
from IPython.core.debugger import Tracer
 
def stability_calc(condition_values, weight_vector): 
    weighted_values = np.multiply(condition_values,weight_vector)
    stability_score = np.linalg.norm(weighted_values)
        
    return stability_score

def promoter_choice_plot(mean_gene_expression, condition_arrays, pct_expression, pct_stability, weight_vector, sort_order = ['stability_score',1], N_plotted = 20):
    #sort_order is either ['stability_score',1] which means sort stability scores in ascending order or
    #['mean_gene_expression',0] which means sort mean_expression scores in descending order
    
    #Calculate stability score
    #set NaNs to zero for calculation of stability score, but save original array
    condition_arrays_na0 = condition_arrays.fillna(0)
    stability_scores = condition_arrays_na0.apply(stability_calc, axis =  1, args = tuple([weight_vector]) )
    stability_scores.name = 'stability_score'
    #Leave NaNs for visualization. 
    
    #Add on mean_gene_expression and stability_score column
    mean_gene_expression.index = condition_arrays.index
    condition_arrays = pd.concat([mean_gene_expression, stability_scores, condition_arrays], axis = 1)
    
    #convert percentile to number of items
    Ngenes = condition_arrays.index.size
    
    expression_indices = [int(np.ceil(Ngenes*(1-percentile))) for percentile in pct_expression]
    stability_indices = [int(np.ceil(Ngenes*(1-percentile))) for percentile in pct_stability]
    
    #Flag data with desired mean expression
    condition_arrays = condition_arrays.sort_values('mean_gene_expression',ascending = 0)
    condition_arrays['top_N_expression'] = False
    condition_arrays['top_N_expression'][expression_indices[1]:expression_indices[0]] = True
    
    #Flag data with desired stability score
    condition_arrays = condition_arrays.sort_values('stability_score',ascending = 1)
    condition_arrays['top_N_stability'] = False
    condition_arrays['top_N_stability'][stability_indices[1]:stability_indices[0]] = True
    
    #extract data with both desired stability and mean scores. 
    desired_genes = condition_arrays.loc[(condition_arrays['top_N_expression'] == True) & (condition_arrays['top_N_stability'] == True)]
    
    desired_genes = desired_genes.sort_values(sort_order[0],ascending = sort_order[1])
    N_plotted = min([len(desired_genes.index),N_plotted])
    plotted_genes = desired_genes.iloc[0:N_plotted,0:-2]
    #plt.pcolor(desired_genes.iloc[:,0:-2])
        
    cmap = mpl.cm.RdBu_r
    cmap.set_bad('k',1.0)
    fig, ax = plt.subplots()
    sns.heatmap(plotted_genes, cmap = cmap, ax = ax)

    return fig, ax
    
    