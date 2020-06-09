#from core import io_library 
#import io_library
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt 
import matplotlib as mpl
from IPython.core.debugger import Tracer
import plotly.graph_objs as pygo
 
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

    return fig, ax, condition_arrays
    
def promoter_ortholog_plot(SC_genenames, species, native_orfs_empirical, condition_arrays):

    #Makes a dictionary to look up orfs by gene names.  This won't include all orfs - those without names had NaN in column 4 so 
    #are presumably left out. 
    SC_orfs_lookup, SC_genename_lookup = io_library.read_SGD_features()
    
    #For a given species read in the ortholog file, make a dictionary
    orth_lookup = io_library.read_orth_lookup_table('Saccharomyces cerevisiae', species)
    
    #For the given list of SC genes find the orthologs
    #This handles a situation where one SC orf has two orthologs. 
    native_orf_orthologs = []
    for genename in SC_genenames: 
        orf = SC_orfs_lookup[genename]
        orf_orthologs = orth_lookup[orf]
        for orf_ortholog in orf_orthologs: 
            native_orf_orthologs.append([orf_ortholog, genename])
    
    orth_lookup_rev = io_library.read_orth_lookup_table(species, 'Saccharomyces cerevisiae')
    
    for native_orf in native_orfs_empirical:
        orf_orthologs_rev = orth_lookup_rev[native_orf]
        for orf_ortholog_rev in orf_orthologs_rev:
            try: 
                SC_genename = SC_genename_lookup[orf_ortholog_rev]
            except KeyError as e:
                print('No S.Cerevisiae ortholog for : ' + native_orf)
                SC_genename = orf_ortholog_rev
            native_orf_orthologs.append([native_orf, SC_genename])

    native_orfs = [orth[0] for orth in native_orf_orthologs]
    SC_genenames = [orth[1] for orth in native_orf_orthologs]
    #Tracer()()
    native_orfs, SC_genenames = remove_duplicate_orfs(native_orfs, SC_genenames)
    
    
    display_order = [jj for jj in range(len(native_orfs))]
    
    native_orf_ortholog_dict = {}
    native_orf_ortholog_dict['display_order'] = display_order
    native_orf_ortholog_dict['SC_genename'] = SC_genenames
    print('reloaded 1011')
    #Tracer()()
    native_orf_ortholog_df = pd.DataFrame(native_orf_ortholog_dict, index = native_orfs)
    
    #Make a dataframe from the list by selecting correct data from condition arrays
    #sorting is required to extract data with slice 
    condition_arrays.sort_index(inplace = True)
    plotted_genes = condition_arrays.loc[(slice(None),[ortholog[0] for ortholog in native_orf_orthologs]),: ]  
    #Get rid of ID index
    #Tracer()()
    plotted_genes.index = plotted_genes.index.droplevel(0)
    
    #Add new columns for each ortholog
    print(species)   
    
    plotted_genes = pd.concat([native_orf_ortholog_df, plotted_genes], axis = 1)
    
    #Make new multiindex with the SC_genename first
    plotted_genes.set_index(['SC_genename',plotted_genes.index], inplace = True)
    plotted_genes.index.names = ['SC_genename', 'native_orf']
    
    #Sort by N_display
    plotted_genes.sort_values('display_order', inplace = True)
    
    #remove expression and stability flags as well as order columns for plotting
    plotted_genes = plotted_genes.iloc[:,1:-2]
    
    #Visualize with heatmap
    cmap = mpl.cm.RdBu_r
    cmap.set_bad('k',1.0)
    fig, ax = plt.subplots()
    sns.heatmap(plotted_genes, cmap = cmap, ax = ax)
            
    return fig, ax
    
def remove_duplicate_orfs(native_orfs, SC_genenames):
    
    remove_ind = []
    for orf in list(set(native_orfs)):
        count = 0
        for ii, jj in enumerate(native_orfs):
            if jj == orf:
                count = count + 1
                if count >1:
                    remove_ind.append(ii)
                    print('duplicate orf mapping for {}, removed {} '.format(orf,SC_genenames[ii]))
                
    remove_ind = sorted(remove_ind, reverse=True)
    for ind in remove_ind:
        SC_genenames.pop(ind)
        native_orfs.pop(ind)
        
    return native_orfs, SC_genenames
    
def ROC_plot(single_expression_dataset, target_gene_set , threshold_type, quantile_range, ax , plot_params):
    #dataset is a series with heading as the name of the condition and fold change values for each orf. 
    #
    #threshold type is either 'activated' or 'repressed'
    #
    #quantile_range is a numpy array that indicates which points will be used in the ROC plot. 'activated'
    #quantiles are typically between 0.5 and 1 and 'repressed' quantiles are between 0 and 0.5.
    #Example activated quantile_range: quantile_range = np.linspace(0.45,0.995,num =  100)
    #Example repressed quantile_range: quantile_range = np.linspace(0.005,0.55,num =  100)
    #
    #plot_params is a set of arguments that determines various properties of the plot
    #eg: plot_params = {'color':'b', 'linewidth': 1.5}
    
    N_genes_total = len(single_expression_dataset.index)
    N_true = len(set(target_gene_set))
    N_false = N_genes_total-N_true
    FPR = []
    TPR = []
    for quantile in quantile_range: 
        if threshold_type == 'activated':
            classified_genes = set(single_expression_dataset[single_expression_dataset > single_expression_dataset.quantile(quantile)].index)
        elif threshold_type =='repressed':
            classified_genes = set(single_expression_dataset[single_expression_dataset < single_expression_dataset.quantile(quantile)].index)
        FPR.append(float(len(classified_genes - set(target_gene_set)))/float(N_false))
        TPR.append(float(len(classified_genes & set(target_gene_set)))/float(N_true))        
    
    ax.plot(FPR,TPR, **plot_params)
    ax.set_xlabel('False Positive Rate', fontsize = 14)
    ax.set_ylabel('True Positive Rate', fontsize = 14)
    
    return TPR, FPR
    
def multi_scatter_plot(expression_data, conditions, xlimit = [], xticks = [], ylimit = [], yticks = []): 
    fig, axarr = plt.subplots(len(conditions), len(conditions), sharex = True, sharey = True)
    for jj in range(len(conditions)): 
        condition_data_J = expression_data[conditions[jj]]
        for kk in range(len(conditions)):
            ax = axarr[jj,kk]
            condition_data_K = expression_data[conditions[kk]]
            ax.scatter(condition_data_K, condition_data_J)
            if len(xlimit) != 0:
                ax.set_xlim(xlimit[0],xlimit[1])
            if len(xticks) !=0:
                ax.set_xticks(xticks)
            if len(ylimit) != 0 : 
                ax.set_ylim(ylimit[0], ylimit[1])
            if len(yticks) != 0:
                ax.set_yticks(yticks)
            if kk == 0:
                ax.set_ylabel(conditions[jj], rotation = 45)
            if jj == 0:
                ax.set_title(conditions[kk])
    return fig, axarr
    
def lfc_padj_plot_with_lines(x_data,y_data,hover_text,lines):
    #lines is a dictionary with the line name as key and the item as
    #[(x1,y1),(x2,y2)]
   

    data = []

    trace = pygo.Scatter(
                x = x_data, 
                y =  y_data,
                text = hover_text,
                mode = 'markers',
                marker = {'opacity': 0.5}, 
                          #'color': 'rgba'+str(cmap(NN/10))}, 
                name = 'genes'
            )
        
    data.append(trace)
    
    for line_name, line in lines.items():

        x1,y1 = line[0]
        x2,y2 = line[1]
        
        x = np.array([x1,x2])
        y = np.array([y1,y2])

        trace2 = pygo.Scatter(
            x = x,
            y = y,
            mode = 'lines',
            marker = {'color': 'black',
                      'size': 5},
            name = line_name
        )
    
        data.append(trace2)
    
    layout = pygo.Layout(
        xaxis= {
            #"range":[-2, 20],
            "title":'LFC'
        },
        yaxis= {
            #"range":[-2, 20],
            "title":'-log10(padj)'
        },
        showlegend=False
    )


    fig = pygo.Figure(data=data, layout = layout) 
    
    return fig

def pval_min_line(line_coords, ymin, x_data): 
    x1,y1 = line_coords[0]
    x2,y2 = line_coords[1]
    
    if y2-y1<0: #Line has negative slope - looking at right side of axis
        x_at_ymin1 = x1
        x_at_ymin2 = max(x_data)
        
    elif y2-y1 >0:   #Line has positive slope - looking at left side of axis
        x_at_ymin1 = x2
        x_at_ymin2 = min(x_data)
    
    p1 = (x_at_ymin1, ymin)
    p2 = (x_at_ymin2, ymin)
    
    return p1,p2