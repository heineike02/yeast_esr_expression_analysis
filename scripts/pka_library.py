# Need to import matplotlib first and set the backend before importing anything else (specifically before pyplot)
# Especially needs to before the imports for core files that may themselves import 
# pyplot without setting the mpl backend first
import matplotlib as mpl
mpl.use('Qt4Agg') # Need this (at least on Mac OSX to get window to display properly)

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
import seaborn as sns
# import seaborn as sns


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


species_abb_dict = { 'SCer' : 'Saccharomyces cerevisiae',
                     'KLac' : 'Kluyveromyces lactis', 
                     'CGla' : 'Candida glabrata', 
                     'SCas' : 'Saccharomyces castellii', 
                     'SBay' : 'Saccharomyces bayanus'}
# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #

def collapse_values_list(values_list):

    values_max = np.zeros(len(values_list[0]))
    values_min = np.zeros(len(values_list[0]))
    values_avg = np.zeros(len(values_list[0]))
    for values in values_list:
        values_max = np.maximum(values_max, values)
        values_min = np.minimum(values_max, values)

    values_avg = np.mean(values_list, axis = 0)

    values_max = values_max.tolist()
    values_min = values_min.tolist()
    values_avg = values_avg.tolist()

    return values_avg, values_max, values_min





# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #


def compile_total_data():

    total_data = {}
    condition_list = []

    # --------------------------------------------------------------- #
    #   # First need to get the single columns from the Regev Data #
    # --------------------------------------------------------------- #
    # This doesn't include Susan's NMPP data
    for species in species_list:
        total_data[species] = {}
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

        for condition in condition_arrays:
            total_data[species][condition] = {}
            condition_key = species + ':' + condition
            condition_list.append(condition_key)
            total_data[species][condition]['Values'] = condition_arrays[condition].values.tolist()
            total_data[species][condition]['Genes'] = condition_arrays[condition].index.get_level_values('orf_name').values.tolist()


    # --------------------------------------------------------------- #
    #   # Add Susan's Data #
    # --------------------------------------------------------------- #
    ORF_gene_table, gene_ORF_table = io_library.read_SGD_features()
    fname = base_dir + '/Susan_UpAndDownRegulatedGenes4fold_MyAndOShea.xlsx'
    f = open(fname)
    data_high = pd.read_excel(f, sheetname = 0)
    high_gene_names = data_high['Genes'].values.tolist()
    high_gene_susan = data_high['My foldchange'].values.tolist()
    high_gene_oshea = data_high['Oshea foldchange'].values.tolist()
    high_gene_ORFs = []
    for gene in high_gene_names:
        if gene in ORF_gene_table:
            high_gene_ORFs.append(ORF_gene_table[gene])
        else:
            high_gene_ORFs.append(gene)


    f = open(fname)
    data_low = pd.read_excel(f, sheetname = 1)
    low_gene_names = data_low['Genes'].values.tolist()
    low_gene_susan = data_low['My foldchange'].values.tolist()
    low_gene_oshea = data_low['Oshea foldchange'].values.tolist()
    low_gene_ORFs = []
    for gene in low_gene_names:
        if gene in ORF_gene_table:
            low_gene_ORFs.append(ORF_gene_table[gene])
        else:
            low_gene_ORFs.append(gene)

    condition_list.append('Saccharomyces cerevisiae:susan')
    total_data['Saccharomyces cerevisiae']['susan'] = {}
    total_data['Saccharomyces cerevisiae']['susan']['Values'] = high_gene_susan + low_gene_susan
    # total_data['Saccharomyces cerevisiae']['susan']['Genes'] = high_gene_names + low_gene_names
    total_data['Saccharomyces cerevisiae']['susan']['Genes'] = high_gene_ORFs + low_gene_ORFs

    condition_list.append('Saccharomyces cerevisiae:oshea')
    total_data['Saccharomyces cerevisiae']['oshea'] = {}
    total_data['Saccharomyces cerevisiae']['oshea']['Values'] = high_gene_oshea + low_gene_oshea
    # total_data['Saccharomyces cerevisiae']['oshea']['Genes'] = high_gene_names + low_gene_names
    total_data['Saccharomyces cerevisiae']['oshea']['Genes'] = high_gene_ORFs + low_gene_ORFs




    # ----------------------------------------------------------------------- #
    #   # Get average, max, and min column values from common conditions #
    #   E.g. four timepoints (columns) in NaCl condition for each species
    # ----------------------------------------------------------------------- #

    # First, find the common conditions based on naming convention
    # Look for the last '_' in the condition name. If the text before
    # the last '_' is the same in a group of columns for a given species,
    # those columns should all be timepoints in the same condition. E.g.
    # NaCl_015, NaCl_030, etc

    # Could incorporate this into the for loop above, but kept 
    # separate for now to make cleaner      

    for species in species_list:
        conditions_to_append = {}
        prev_condition = ''
        test_name_prev = prev_condition.split('_')
        test_name_prev = test_name_prev[len(test_name_prev) - 2]
        temp_values_list = []

        for condition in total_data[species]:
            test_name_current = condition.split('_')
            test_name_current = test_name_current[len(test_name_current) - 2]
            if test_name_prev == test_name_current:
                temp_values_list.append(total_data[species][condition]['Values'])
            else: 

                if len(temp_values_list) > 0:
                    values_avg, values_max, values_min = collapse_values_list(temp_values_list)
                    genes = total_data[species][prev_condition]['Genes']
                    condition_name = test_name_prev
                    
                    condition_name_avg = condition_name + '_avg'
                    condition_name_max = condition_name + '_max'
                    condition_name_min = condition_name + '_min'
                    condition_list.append(species + ':' + condition_name_avg)
                    condition_list.append(species + ':' + condition_name_max)
                    condition_list.append(species + ':' + condition_name_min)


                    conditions_to_append[condition_name_avg] = {'Genes' : genes, 'Values' : values_avg}
                    conditions_to_append[condition_name_max] = {'Genes' : genes, 'Values' : values_max}
                    conditions_to_append[condition_name_min] = {'Genes' : genes, 'Values' : values_min}


                temp_values_list = [total_data[species][condition]['Values']]
                prev_condition = condition
                test_name_prev = test_name_current

        for condition in conditions_to_append:
            total_data[species][condition] = conditions_to_append[condition]


    return total_data, condition_list
# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #


def data_by_key(condition_key, total_data):
    species = condition_key.split(':')[0]
    condition = condition_key.split(':')[1]
    return total_data[species][condition]


# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #

def get_seed_data(activated_number, repressed_number, gene_list, total_data, condition_key, mode = 'top'):
    # Mode can be 'top', 'all', or 'name'
    high_genes = []
    low_genes = []
    high_values = []
    low_values = []
    gene_value_map = {}
    value_gene_map = {}
    ORF_map = {}
    ORFs = []

    data = data_by_key(condition_key, total_data)

    data_frame = pd.DataFrame({'Genes' : data['Genes'], 'Values' : data['Values']})
    data_frame_high = data_frame.sort_values('Values', ascending = False)
    data_frame_low = data_frame.sort_values('Values', ascending = True)

    if mode == 'top':
        activated_index = activated_number
        repressed_index = repressed_number
    else:
        activated_index = len(data_frame_high.values.tolist())
        repressed_index = len(data_frame_low.values.tolist())

    high_data = data_frame_high.values.tolist()[0:activated_index]
    low_data = data_frame_low.values.tolist()[0:repressed_index]

    for data_pair in high_data:
        if mode != 'name' or data_pair[0] in gene_list:
            high_genes.append(data_pair[0])
            high_values.append(data_pair[1])
            gene_value_map[data_pair[0]] = data_pair[1]
            value_gene_map[data_pair[1]] = data_pair[0]

    for data_pair in low_data:
        if mode != 'name' or data_pair[0] in gene_list:
            low_genes.append(data_pair[0])
            low_values.append(data_pair[1])
            gene_value_map[data_pair[0]] = data_pair[1]
            value_gene_map[data_pair[1]] = data_pair[0]

    genes = high_genes + low_genes
    values = high_values + low_values


    if get_species(condition_key) == 'Saccharomyces cerevisiae':
        ORF_gene_table, gene_ORF_table = io_library.read_SGD_features()
        for gene in genes:
            if gene in ORF_gene_table:
                ORF_map[gene] = ORF_gene_table[gene]
                ORFs.append(ORF_gene_table[gene])
            else:
                ORF_map[gene] = gene
                ORFs.append(gene)
    else:
        for gene in genes:
            ORF_map[gene] = gene
            ORFs.append(gene)

    seed_data = {}
    seed_data['Genes'] = genes
    seed_data['Values'] = values
    seed_data['High Genes'] = high_genes
    seed_data['Low Genes'] = low_genes
    seed_data['High Values'] = high_values
    seed_data['Low Values'] = low_values
    seed_data['Condition Key'] = condition_key
    seed_data['Gene Map'] = gene_value_map
    seed_data['Value Map'] = value_gene_map
    seed_data['ORF Map'] = ORF_map
    seed_data['ORFs'] = ORFs

    return seed_data
# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #

def create_condition_key(species_arg, condition):
    species = ''
    if species_arg in species_abb_dict:
        species = species_abb_dict[species_arg]
    elif species in species_list:
        species = species_arg
    else:
        raise KeyError('Invalid Species Name : ' + species_arg)

    condition_key = species + ':' + condition

    return condition_key

# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #
def get_species(condition_key):
    return condition_key.split(':')[0]

def get_condition(condition_key):
    return condition_key.split(':')[1]
# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #

def create_cross_data(seed_data, total_data, condition_list):

    cross_data = {}

    seed_genes = seed_data['Genes']
    seed_genes_ORFs = seed_data['ORFs']
    seed_values = seed_data['Values']
    seed_condition = seed_data['Condition Key']
    seed_species = get_species(seed_condition)

    for condition in condition_list:
        if condition != seed_condition:
            cross_data[condition] = {}
            cross_data[condition]['Genes'] = []
            cross_data[condition]['Values'] = []
            cross_data[condition]['Gene Map'] = {}
            cross_data[condition]['Value Map'] = {}
            cross_data[condition]['Seed Gene Map'] = {}

            
            current_species = get_species(condition)
            if seed_species != current_species:
                orth_table = io_library.read_orth_lookup_table(seed_species, current_species)
                
                for seed_gene_ORF in seed_genes_ORFs:

                    cross_data[condition]['Seed Gene Map'][seed_gene_ORF] = []

                    if seed_gene_ORF not in orth_table:
                        cross_data[condition]['Genes'].append('NONE')
                        # cross_data[condition]['Values'].append('NONE')
                        cross_data[condition]['Seed Gene Map'][seed_gene_ORF].append('NONE')
                    elif orth_table[seed_gene_ORF] == ['NONE']:
                        cross_data[condition]['Genes'].append('NONE')
                        # cross_data[condition]['Values'].append('NONE')
                        cross_data[condition]['Seed Gene Map'][seed_gene_ORF].append('NONE')
                    else:
                        for gene in orth_table[seed_gene_ORF]:
                            cross_data[condition]['Genes'].append(gene)
                            cross_data[condition]['Seed Gene Map'][seed_gene_ORF].append(gene)
            else:
                cross_data[condition]['Genes'] = seed_genes_ORFs
                for gene in seed_genes_ORFs:
                    cross_data[condition]['Seed Gene Map'][gene] = [gene]



    for condition in cross_data:
        search_species = get_species(condition)
        search_condition = get_condition(condition)
        
        # if search_species != seed_species:
        for search_gene in cross_data[condition]['Genes']:
            if search_gene == 'NONE':
                cross_data[condition]['Values'].append(float('nan'))
                cross_data[condition]['Gene Map'][search_gene] = float('nan')
            elif total_data[search_species][search_condition]['Genes'].count(search_gene) == 0:
                cross_data[condition]['Values'].append(float('nan'))
                cross_data[condition]['Gene Map'][search_gene] = float('nan')
            else:
                gene_index = total_data[search_species][search_condition]['Genes'].index(search_gene)
                gene_value = total_data[search_species][search_condition]['Values'][gene_index]
                cross_data[condition]['Values'].append(gene_value)
                cross_data[condition]['Gene Map'][search_gene] = gene_value
                cross_data[condition]['Value Map'][gene_value] = search_gene

        # else:
        #     for search_gene in cross_data[condition]['Genes']:
        #         gene_index = total_data[search_species][search_condition]['Genes'].index(search_gene)
        #         gene_value = total_data[search_species][search_condition]['Values'][gene_index]
        #         cross_data[condition]['Values'].append(gene_value)
        #         cross_data[condition]['Gene Map'][search_gene] = gene_value
        #         cross_data[condition]['Value Map'][gene_value] = search_gene

    
    return cross_data


# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #

def sort_conditions(conditions):

    sorted_conditions = list(conditions)
    sorted_conditions.sort()
    return sorted_conditions



# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #

def create_plotting_data(seed_data, cross_data):

    repeats = {}
    plotting_columns = []
    plotting_labels_x = []

    seed_genes = seed_data['ORFs']

    for seed_gene in seed_genes:
        repeats[seed_gene] = 1

    plotting_seed_values = seed_data['Values']
    plotting_labels_y = seed_data['Genes']

    num_columns = 1
    for condition in cross_data:
        num_columns += 1
        for seed_gene in seed_genes:
            repeats[seed_gene] = max(repeats[seed_gene], len(cross_data[condition]['Seed Gene Map'][seed_gene]))

    for seed_gene in seed_genes:
        for i in range(repeats[seed_gene] - 1):
            gene_index = plotting_labels_y.index(seed_gene)
            gene_value = plotting_seed_values[gene_index]
            plotting_seed_values.insert(gene_index, gene_value)
            plotting_labels_y.insert(gene_index, ' ')

    num_rows = len(plotting_seed_values)
    plotting_data = np.zeros((num_rows, num_columns))
    plotting_data[:, 0] = plotting_seed_values


    if 'susan' in seed_data['Condition Key'] or 'oshea' in seed_data['Condition Key']:
        plotting_labels_y = []
        ORF_gene_table, gene_ORF_table = io_library.read_SGD_features()
        for seed_gene in seed_data['Genes']:
            if seed_gene != ' ':
                plotting_labels_y.append(gene_ORF_table[seed_gene])
            else:
                plotting_labels_y.append(' ')

    plotting_labels_x = sort_conditions(cross_data.keys())

    # for column_counter, condition in enumerate(cross_data):
    for column_counter, condition in enumerate(plotting_labels_x):
        condition_genes = cross_data[condition]['Genes']
        condition_values = []
        # plotting_labels_x.append(condition)

        for seed_gene in seed_genes:


            for gene in cross_data[condition]['Seed Gene Map'][seed_gene]:
                if gene != 'NONE':
                    # print cross_data[condition]['Gene Map']
                    condition_values.append(cross_data[condition]['Gene Map'][gene])
                else:
                    condition_values.append(float('nan'))
            for i in range(repeats[seed_gene] - len(cross_data[condition]['Seed Gene Map'][seed_gene])):
                condition_values.append(float('nan'))

        plotting_data[:, column_counter] = condition_values





    return plotting_data, plotting_labels_x, plotting_labels_y


# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #

def plot_data(plotting_data, plotting_labels_x, plotting_labels_y):


    species_label_indices = np.ones(len(species_list))
    species_label_indices *= float('nan')
    # plotting_data_search = plotting_data[ : , 1: ]

    for species_index, species in enumerate(species_list):
        for condition_index, condition in enumerate(plotting_labels_x):
            if species in condition:
                species_label_indices[species_index] = condition_index + 1
                break

    for label_index, condition in enumerate(plotting_labels_x):
        plotting_labels_x[label_index] = condition[0:5]



    num_columns = len(plotting_labels_y) + 1
    cmap = mpl.cm.RdBu_r
    cmap.set_bad('k',1.0)
    fig, ax = plt.subplots()
    ax = sns.heatmap(plotting_data, cmap = cmap, ax = ax, yticklabels = plotting_labels_y, xticklabels = plotting_labels_x)
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
    plt.get_current_fig_manager().window.raise_()
    plt.show()


# --------------------------------------------------------------------------------- #
# #### MAIN #####
# --------------------------------------------------------------------------------- #

total_data, condition_list = compile_total_data()
condition_key = create_condition_key('SCer', 'susan')

seed_gene_list = []
seed_data = get_seed_data(10, 10, seed_gene_list, total_data, condition_key)

cross_data = create_cross_data(seed_data, total_data, condition_list)
plotting_data, plotting_labels_x, plotting_labels_y = create_plotting_data(seed_data, cross_data)
plot_data(plotting_data, plotting_labels_x, plotting_labels_y)


# plotting_labels_x.sort()
# for element in plotting_labels_x:
    # print element



