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
import seaborn as sns


# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #


# Parameters for choosing up and down regulated genes in S.Cer pka inhibition
POSITIVE_LOG_CUTOFF = 6.5
NEGATIVE_LOG_CUTOFF = -6.0


# Create ortholog lookup tables for each species
# species_list = ['Kluyveromyces lactis', 'Saccharomyces castellii', 
                 # 'Candida glabrata', 'Saccharomyces bayanus', 'Saccharomyces cerevisiae']


species_list = ['Saccharomyces cerevisiae', 'Saccharomyces bayanus', 
				'Candida glabrata', 'Saccharomyces castellii',
				'Kluyveromyces lactis']
# Refernce dictionary for creating file names
species_name_dict = {'Saccharomyces cerevisiae' : 'SCer',
                    'Kluyveromyces lactis': 'KLac', 
                    'Candida glabrata' : 'CGla', 
                    'Saccharomyces castellii' : 'SCas', 
                    'Saccharomyces bayanus' : 'SBay'}

# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #

def get_genes_susan(positive_cutoff, negative_cutoff):
    fname = base_dir + '/Susan_UpAndDownRegulatedGenes4fold_MyAndOShea.xlsx'
    f = open(fname)
    data_high = pd.read_excel(f, sheetname = 0)
    sorted_data_high = data_high.sort_values('My foldchange', ascending = False)
    sorted_data_high.reset_index(drop = True, inplace = True)
    susan_foldchange = sorted_data_high['My foldchange']

    index = -1
    test_value = susan_foldchange[0]

    while test_value > positive_cutoff:
        test_value = susan_foldchange[index + 1]
        index += 1

    high_gene_names = sorted_data_high['Genes'][0 : index]
    high_gene_change = sorted_data_high['My foldchange'][0 : index]
    high_gene_data = {'Gene' : high_gene_names, 'Log_Change' : high_gene_change}
    high_genes = pd.DataFrame(high_gene_data)
    f.close()

    f = open(fname)
    data_low = pd.read_excel(f, sheetname = 1)
    sorted_data_low = data_low.sort_values('My foldchange')
    sorted_data_low.reset_index(drop = True, inplace = True)
    susan_foldchange = sorted_data_low['My foldchange']

    index = -1
    test_value = susan_foldchange[0]

    while test_value < negative_cutoff:
        test_value = susan_foldchange[index + 1]
        index += 1
        
    low_gene_names = sorted_data_low['Genes'][0 : index]
    low_gene_change = sorted_data_low['My foldchange'][0 : index]
    low_gene_data = {'Gene' : low_gene_names, 'Log_Change' : low_gene_change}
    low_genes = pd.DataFrame(low_gene_data)
    
    return high_genes, low_genes

# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #

def concat_SC_genes(high_genes, low_genes):
	# Grab ortholog data from tab file and create lookup tables
	SC_orfs_lookup, SC_genename_lookup = io_library.read_SGD_features()

	# Gather all the ORF names for the high genes
	high_orfs = []
	for gene_name in high_genes['Gene']:    
		if gene_name in SC_genename_lookup:
			high_orfs.append(gene_name)
		else:
			high_orfs.append(SC_orfs_lookup[gene_name])


	# Append the orf names to the existing high_gene data
	high_gene_data = {'Gene' : high_genes['Gene'], 'Log_Change' : high_genes['Log_Change'],
	                 'ORF' : high_orfs}
	high_genes = pd.DataFrame(high_gene_data)

	# Repeat with the low genes
	low_orfs = []
	for gene_name in low_genes['Gene']:    
	    if gene_name in SC_genename_lookup:
	        low_orfs.append(gene_name)
	    else:
	        low_orfs.append(SC_orfs_lookup[gene_name])

	# Append the orf names to the existing high_gene data
	low_gene_data = {'Gene' : low_genes['Gene'], 'Log_Change' : low_genes['Log_Change'],
	                 'ORF' : low_orfs}
	low_genes = pd.DataFrame(low_gene_data)

	SC_gene_data = {'Gene' : high_genes['Gene'] + low_genes['Gene'], 
	                'Log_Change' : high_genes['Log_Change'] + low_genes['Log_Change'],
	                'ORF' : low_orfs + high_orfs}


	SC_genes = pd.concat([high_genes, low_genes])

	return SC_genes


# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #

def compile_species_data():

	species_data = {}
	SC_genes_repeats = {}
	conditions = []
	for species in species_list:
	    SC_genes_repeats[species] = {}
	    species_data[species] = {}
	    species_data[species]['condition_names'] = []
	    orth_dict = io_library.read_orth_lookup_table('Saccharomyces cerevisiae', species)

	    fname = os.path.normpath(base_dir + 
	                             '/scripts/expression_broad_data_datafiles/microarray_data/GSE36253_Growth/' + 
	                             species_name_dict[species] + '_growth.csv' )
	    growth_exp = pd.read_csv(fname,header = [0,1,2], index_col = [0,1])
	    growth_replicate_groups = growth_exp.groupby(axis = 1, level = 'conditions')
	    growth_exp_avg = growth_replicate_groups.aggregate(np.mean)



	    if species != 'Saccharomyces bayanus':
	        #There is no stress dataset for S. bayanus
	        fname = os.path.normpath(base_dir + 
	                                 '/scripts/expression_broad_data_datafiles/microarray_data/GSE38478_Stress/' + 
	                                 species_name_dict[species] + '_stress.csv' )
	        stress_exp = pd.read_csv(fname,header = [0,1,2], index_col = [0,1])
	        # Group by condition and take mean
	        stress_replicate_groups = stress_exp.groupby(axis = 1, level = 'conditions')
	        stress_exp_avg = stress_replicate_groups.aggregate(np.mean)

	        if False in stress_exp_avg.index==growth_exp_avg.index:
	            print "Error: ID mismatch between condition data. Species = {}".format(species)
	            break
	        condition_arrays = pd.concat([growth_exp_avg,stress_exp_avg], axis = 1)	
	    else:
	        condition_arrays = growth_exp_avg

	    num_conditions = 0
	    for condition in condition_arrays:
	        species_data[species]['condition_names'].append(condition)
	        if condition not in conditions:
	            conditions.append(condition)
	        species_data[species][condition] = condition_arrays[condition]
	        num_conditions += 1

	    species_data[species]['num_conditions'] = num_conditions

	    for condition_index, condition in enumerate(condition_arrays):
	        species_data[species][condition_index] = []

	    species_data[species]['orth_genes'] = []
	    for gene in SC_genes['ORF']:
	        SC_genes_repeats[species][gene] = 1
	        SC_genes_repeats[gene] = 1
	        if species != 'Saccharomyces cerevisiae':
	            try: 
	                orth_gene_list = orth_dict[gene]
	            except KeyError:
	                orth_gene_list = ['NONE']

	        else:
	            orth_gene_list = [gene]


	        if SC_genes_repeats[species][gene] < len(orth_gene_list):
	            SC_genes_repeats[species][gene] = len(orth_gene_list)
	        for orth_gene in orth_gene_list:
	            species_data[species]['orth_genes'].append(orth_gene)

	            for condition_index, condition in enumerate(condition_arrays):
	                condition_data = species_data[species][condition]

	                if orth_gene != 'NONE':
	                    species_data[species][condition_index].append(species_data[species][condition].loc[(slice(None), orth_gene)].values[0])
	                else:
	                    species_data[species][condition_index].append(float('nan'))

	return species_data, conditions, SC_genes_repeats

# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #

def create_plotting_data(SC_genes, SC_genes_repeats):
    SC_gene_plot_data = SC_genes['Log_Change'].values.tolist()    

    num_columns = 1
    num_rows = 0
    SC_genes_labels = []
    for species in species_data:
        num_columns += species_data[species]['num_conditions']

        for gene_index, gene in enumerate(SC_genes['ORF']):
            if SC_genes_repeats[gene] < SC_genes_repeats[species][gene]:
                SC_genes_repeats[gene] = SC_genes_repeats[species][gene]

    num_rows = 0
    row_index = 0
    for gene_index, gene in enumerate(SC_genes['ORF']):
        num_rows += SC_genes_repeats[gene]
        # Switch this append call depending on whether you want the column labels
        # to be the standardized ORF name or the common gene name
        # SC_genes_labels.append(gene)
        SC_genes_labels.append(SC_genes['Gene'].values[gene_index])
        # SC_genes_labels.append(str(SC_genes['Gene'].values[gene_index]) + ' ' + str(SC_gene_plot_data[gene_index]))
        row_index += 1
        for i in range(SC_genes_repeats[gene] - 1):
            SC_gene_plot_data.insert(row_index, SC_gene_plot_data[row_index - 1])
            SC_genes_labels.append(' ')    
            row_index += 1

      #   for i in range(SC_genes_repeats[gene] - 1):
	    	# # print '------------------------------'
	    	# # print SC_genes['Gene'].values[gene_index]
	    	# # print SC_genes_repeats[gene] - 1
	    	# # print i
	    	# # print '------------------------------'
	     #    SC_gene_plot_data.insert(gene_index + 1, SC_gene_plot_data[gene_index])
	     #    # SC_gene_plot_data.insert(0, SC_gene_plot_data[gene_index])
	     #    # SC_genes_labels.append(' ')
	     #    SC_genes_labels.append(str(SC_genes['Gene'].values[gene_index]) + ' ' + str(SC_gene_plot_data[gene_index]))


    plotting_data = np.zeros((num_rows, num_columns))
    plotting_data[:, 0] = SC_gene_plot_data

    return plotting_data, SC_genes_labels, num_columns

# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #

def plot_data(plotting_data, species_data):
	condition_labels = ['Susan']
	species_label_indices = []
	col_index = 1
	col_offset = 0
	for species_index, species in enumerate(species_list):

	    species_label_indices.append(col_index)
	    
	    repeat_diff = {}
	    
	    for gene_index, gene in enumerate(SC_genes['ORF']):
	        repeat_diff[gene] = SC_genes_repeats[gene] - SC_genes_repeats[species][gene]

	    for condition_index in range(species_data[species]['num_conditions']):
	        row_index = 0
	        condition_plotting_data = species_data[species][condition_index]
	        for gene_index, gene in enumerate(SC_genes['ORF']):
	            row_index += 1
	            for i in range(repeat_diff[gene]):
	                condition_plotting_data.insert(row_index, float('nan'))
	                row_index += 1

	        condition_labels.append(species_data[species]['condition_names'][condition_index])
	        plotting_data[:, col_index] = condition_plotting_data


	        col_index += 1
	cmap = mpl.cm.RdBu_r
	cmap.set_bad('k',1.0)
	fig, ax = plt.subplots()
	ax = sns.heatmap(plotting_data, cmap = cmap, ax = ax, yticklabels = SC_genes_labels, xticklabels = condition_labels)
	ylocs, ylabels = plt.yticks()
	plt.setp(ylabels, rotation = 0)
	xlocs, xlabels = plt.xticks()
	plt.setp(xlabels, rotation = 90)
	species_labels = []
	for species in species_list:
	    species_labels.append(species_name_dict[species])
	newax = ax.twiny()
	for i in range(len(species_label_indices)):
	    species_label_indices[i] = species_label_indices[i] / float(num_columns)

	newax.set_xticklabels(species_labels)
	newax.set_xticks(species_label_indices)
	plt.tight_layout()
	plt.show()

# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #


# Now we have a list of the genes in S.Cer whose log changes are high and low respectively in the
# pka experiment
high_genes, low_genes = get_genes_susan(POSITIVE_LOG_CUTOFF, NEGATIVE_LOG_CUTOFF)
SC_genes = concat_SC_genes(high_genes, low_genes)
species_data, conditions, SC_genes_repeats = compile_species_data()
plotting_data, SC_genes_labels, num_columns = create_plotting_data(SC_genes, SC_genes_repeats)
plot_data(plotting_data, species_data)












