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
    fname = base_dir + '/Susan_UpAndDownRegulatedGenes4fold_MyAndOShea.xlsx'
    f = open(fname)
    data_high = pd.read_excel(f, sheetname = 0)
    high_gene_names = data_high['Genes'].values.tolist()
    high_gene_susan = data_high['My foldchange'].values.tolist()
    high_gene_oshea = data_high['Oshea foldchange'].values.tolist()


    f = open(fname)
    data_low = pd.read_excel(f, sheetname = 1)
    low_gene_names = data_low['Genes'].values.tolist()
    low_gene_susan = data_low['My foldchange'].values.tolist()
    low_gene_oshea = data_low['Oshea foldchange'].values.tolist()

    condition_list.append('Saccharomyces cerevisiae:susan')
    total_data['Saccharomyces cerevisiae']['susan'] = {}
    total_data['Saccharomyces cerevisiae']['susan']['Values'] = high_gene_susan + low_gene_susan
    total_data['Saccharomyces cerevisiae']['susan']['Genes'] = high_gene_names + low_gene_names

    condition_list.append('Saccharomyces cerevisiae:oshea')
    total_data['Saccharomyces cerevisiae']['oshea'] = {}
    total_data['Saccharomyces cerevisiae']['oshea']['Values'] = high_gene_oshea + low_gene_oshea
    total_data['Saccharomyces cerevisiae']['oshea']['Genes'] = high_gene_names + low_gene_names




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

def get_seed_data(activated_number, repressed_number, total_data, condition_key):
    high_genes = []
    low_genes = []
    high_values = []
    low_values = []
    gene_value_map = {}
    value_gene_map = {}
    ORF_map = {}

    data = data_by_key(condition_key, total_data)

    data_frame = pd.DataFrame({'Genes' : data['Genes'], 'Values' : data['Values']})
    data_frame_high = data_frame.sort_values('Values', ascending = False)
    data_frame_low = data_frame.sort_values('Values', ascending = True)

    high_data = data_frame_high.values.tolist()[0:activated_number]
    low_data = data_frame_low.values.tolist()[0:repressed_number]

    for data_pair in high_data:
        high_genes.append(data_pair[0])
        high_values.append(data_pair[1])
        gene_value_map[data_pair[0]] = data_pair[1]
        value_gene_map[data_pair[1]] = data_pair[0]

    for data_pair in low_data:
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
            else:
                ORF_map[gene] = gene
    else:
        for gene in genes:
            ORF_map[gene] = gene

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

def create_cross_data(seed_data, total_data, condition_list, selection = 'All'):

    cross_data = {}
    condition_key = seed_data['Condition Key']
    seed_species = get_species(condition_key)
    
    if selection == 'High':
        genes_selector = 'High Genes'
    elif selection == 'Low':
        genes_selector = 'Low Genes'
    else:
        genes_selector = 'Genes'

    for condition in condition_list:
        print condition


        if condition != condition_key:
            condition_name = get_condition(condition)
            cross_data[condition] = {}
            species = get_species(condition)
            # orth_table_forward = io_library.read_orth_lookup_table(seed_species, species)
            orth_table_backward = io_library.read_orth_lookup_table(species, seed_species)
            genes_list = total_data[species][condition_name]['Genes']
            values_list = total_data[species][condition_name]['Values']

            print len(genes_list)
            print len(values_list)

            indices_to_remove = []
            # print condition
            # print genes_list
            for gene_index, gene in enumerate(genes_list):
                # print species
                # print gene
                if gene not in orth_table_backward:
                    indices_to_remove.append(gene_index)
                    # continue
                # elif orth_table_backward[gene] not in seed_data[genes_selector]:
                elif 'NONE' in orth_table_backward[gene]:
                    # continue
                    indices_to_remove.append(gene_index)

                else:
                    continue

                #     in_seed_data = False
                #     for lookup_gene in orth_table_backward[gene]:
                #         if lookup_gene in seed_data['ORF Map']:
                #             in_seed_data = True
                # # elif orth_table_backward[gene] not in seed_data['ORF Map']:
                #     if not in_seed_data:
                #         indices_to_remove.append(gene_index)


            # print indices_to_remove
            for index_counter, index in enumerate(indices_to_remove):
                genes_list.pop(index - index_counter)
                values_list.pop(index - index_counter)
                # print 'genes list = ' + str(len(genes_list))
                # print 'values list = ' + str(len(values_list))


            cross_data[condition]['Genes'] = genes_list
            cross_data[condition]['Values'] = values_list

            
            print len(genes_list)
            print len(values_list)

    return cross_data




# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #

def create_plotting_data(seed_data, cross_data, selection = 'All'):

    if selection == 'High':
        genes_selector = 'High Genes'
        values_selector = 'High Values'
    elif selection == 'Low':
        genes_selector = 'Low Genes'
        values_selector = 'Low Values'
    else:
        genes_selector = 'Genes'
        values_selector = 'Values'

    gene_repeats = np.ones(len(seed_data[genes_selector]))
    repeat_numbers = {}
    seed_species = get_species(seed_data['Condition Key'])

    for condition in cross_data:
        species = get_species(condition)
        orth_table = io_library.read_orth_lookup_table(seed_species, species)

        # print orth_table
        for gene_index, gene in enumerate(seed_data[genes_selector]):
            lookup_name = seed_data['ORF Map'][gene]
            if lookup_name in orth_table:
                gene_repeats[gene_index] = max(gene_repeats[gene_index], len(orth_table[lookup_name]))

    for condition in cross_data:
        # print condition
        # print cross_data[condition]
        repeat_numbers[condition] = np.ones(len(cross_data[condition]['Genes']))

        for gene_index, gene in enumerate(seed_data[genes_selector]):
            lookup_name = seed_data['ORF Map'][gene]
            if lookup_name in orth_table:
                # print len(repeat_numbers[condition])
                repeat_numbers[condition][gene_index] += gene_repeats - len(orth_table[gene])
            else:
                repeat_numbers[condition][gene_index] = -1



    # num_columns = 





    return species

# --------------------------------------------------------------------------------- #
# #### MAIN #####
# --------------------------------------------------------------------------------- #

# high_genes, low_genes = get_seed_genes(10, 10, SC_source = 'oshea')

# orth_table_test = io_library.read_orth_lookup_table('Candida glabrata', 'Saccharomyces castellii')


total_data, condition_list = compile_total_data()

condition_key = create_condition_key('SCer', 'susan')
seed_data = get_seed_data(10, 10, total_data, condition_key)

print total_data['Saccharomyces cerevisiae']['hydrogen peroxide_avg']

# cross_data = create_cross_data(seed_data, total_data, condition_list)
# condition_key = cross_data.keys()[1]
# print cross_data[condition_key]



# plotting_data = create_plotting_data(seed_data, cross_data)
# print cross_data.keys()
# print len(cross_data)






