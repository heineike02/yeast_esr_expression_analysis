import os
base_dir = os.path.normpath(os.path.dirname(os.getcwd()))
import sys
sys.path.append(base_dir + '/core')
import io_library
from IPython.core.debugger import Tracer
from IPython.core.debugger import Tracer
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import matplotlib as mpl


# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #


POSITIVE_CUTOFF = 'none'
NEGATIVE_CUTOFF = 'none'


# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #

def get_genes(positive_cutoff, negative_cutoff):
    fname = base_dir + '/Susan_UpAndDownRegulatedGenes4fold_MyAndOShea.xlsx'
    f = open(fname)
    data_high = pd.read_excel(f, sheetname = 0)
    sorted_data_high = data_high.sort_values('My foldchange', ascending = False)
    sorted_data_high.reset_index(drop = True, inplace = True)
    susan_foldchange = sorted_data_high['My foldchange']

    index = -1
    test_value = susan_foldchange[0]

    if positive_cutoff == 'none':
    	index = len(susan_foldchange)

    else:
	    while test_value > positive_cutoff:
	        test_value = susan_foldchange[index + 1]
	        index += 1

    high_gene_names = sorted_data_high['Genes'][0 : index]
    high_gene_change = sorted_data_high['My foldchange'][0 : index]
    high_gene_oshea = sorted_data_high['Oshea foldchange'][0 : index]

    high_gene_data = {'Gene' : high_gene_names, 'Log_Change' : high_gene_change, 'Oshea' : high_gene_oshea}
    high_genes = pd.DataFrame(high_gene_data)
    f.close()

    f = open(fname)
    data_low = pd.read_excel(f, sheetname = 1)
    sorted_data_low = data_low.sort_values('My foldchange')
    sorted_data_low.reset_index(drop = True, inplace = True)
    susan_foldchange = sorted_data_low['My foldchange']

    index = -1
    test_value = susan_foldchange[0]

    if negative_cutoff == 'none':
    	index = len(susan_foldchange)

    else:
	    while test_value < negative_cutoff:
	        test_value = susan_foldchange[index + 1]
	        index += 1
        
    low_gene_names = sorted_data_low['Genes'][0 : index]
    low_gene_change = sorted_data_low['My foldchange'][0 : index]
    low_gene_oshea = sorted_data_low['Oshea foldchange'][0 : index]
    low_gene_data = {'Gene' : low_gene_names, 'Log_Change' : low_gene_change, 'Oshea' : low_gene_oshea}
    low_genes = pd.DataFrame(low_gene_data)
    

    concat_genes = pd.concat([high_genes, low_genes])
    return high_genes, low_genes, concat_genes


# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #

def create_corr_plot(high_genes, low_genes, concat_genes):
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	#  Upregulated (Susan) genes
	high_genes_susan = high_genes['Log_Change'].values.tolist()
	high_genes_oshea = high_genes['Oshea'].values.tolist()	
	m, b = np.polyfit(high_genes_susan, high_genes_oshea, 1)
	plt.plot(high_genes_susan, high_genes_oshea, '*')
	xlims = ax.get_xlim()
	x = np.linspace(xlims[0], xlims[1], 100)
	high_best_fit = m * x + b
	corr_coef = np.corrcoef(high_genes_susan, high_genes_oshea)[1,0]
	print corr_coef
	plt.plot(x, high_best_fit)
	plt.title('Susan Upregulated Genes')
	plt.xlabel('Susan Logfold Change')
	plt.ylabel('Oshea Logfold Change')
	ax.text(9,-1, '$r = 0.61144709392$')
	plt.show()


	#  Downregulated (Susan) genes
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	low_genes_susan = low_genes['Log_Change'].values.tolist()
	low_genes_oshea = low_genes['Oshea'].values.tolist()
	m, b = np.polyfit(low_genes_susan, low_genes_oshea, 1)	
	plt.plot(low_genes_susan, low_genes_oshea, '*')
	xlims = ax.get_xlim()
	x = np.linspace(xlims[0], xlims[1], 100)
	low_best_fit = m * x + b
	corr_coef = np.corrcoef(low_genes_susan, low_genes_oshea)[1,0]
	ax.text(-8,-5.5, '$r = 0.157239640506$')
	print corr_coef
	plt.plot(x, low_best_fit)
	plt.title('Susan Downregulated Genes')
	plt.xlabel('Susan Logfold Change')
	plt.ylabel('Oshea Logfold Change')
	plt.show()

# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #



high_genes, low_genes, concat_genes = get_genes(POSITIVE_CUTOFF, NEGATIVE_CUTOFF)
create_corr_plot(high_genes, low_genes, concat_genes)


