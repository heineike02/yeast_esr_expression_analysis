# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
import regev_library
import numpy as np
import os
base_dir = os.path.normpath(os.path.dirname(os.getcwd()))
import sys
sys.path.append(base_dir + '/core')
import io_library
import matplotlib as mpl
import seaborn as sns
import matplotlib.pyplot as plt
import params
import scipy
# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
# 
# This library got used a lot more than the single score one since it invovles calculating the cross-species correlation
# scores and manipulating them. The two libraries have a similar layout though and similar funciton names
# 
# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
def is_uniform(test_list):
    """
    Auxiliary function for checking if a list is made up of identical values

    Args;
        test_list: List of values of test

    Returns:
        is_uniform: True if the list is made up of identical values.
                    False otherwise
    """
    is_uniform = True

    test_value = test_list[0]

    for value in test_list:
        if value != test_value:
            is_uniform = False

    return is_uniform


# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------

def get_scores(total_data, condition_list):
    """
    Get the (total) score data corresponding all the genes in total_data (and the conditions in condition_list)
    
    Args:
        total_data: A total data struct 
        condition_list: The list of conditions corresponding to the total data struct

    Returns:
        score_data: Cross-species score data for all genes in the total data struct
    """
    score_data = {}
    gene_list_susan = []

    # gene_list = total_data[]
    # print total_data.keys()

    for species_counter, current_species in enumerate(total_data):
        condition_temp = total_data[current_species].keys()[0]
        if condition_temp == 'susan' or condition_temp == 'oshea':
            condition_temp = total_data[current_species].keys()[1]
        if condition_temp == 'susan' or condition_temp == 'oshea':
            condition_temp = total_data[current_species].keys()[2]
        gene_list = total_data[current_species][condition_temp]['Genes']

        for condition in total_data[current_species]:
            if gene_list != total_data[current_species][condition]['Genes']:
                for gene in total_data[current_species][condition]['Genes']:
                    # if species == 'Saccharomyces cerevisiae':
                    #     ORF_gene_table, gene_ORF_table = io_library.read_SGD_features()

                    if gene not in gene_list and gene not in gene_list_susan:
                        # print current_species, condition
                        gene_list_susan.append(gene)
        
        condition_key = current_species + ':' + condition
        print '--------------------------------------------------'
        print condition_key
        print '--------------------------------------------------'

        condition_list = regev_library.sort_conditions(condition_list)

        # ----------------------------------------------------
        # index_cutoff = 25
        index_cutoff = 'all'
        if index_cutoff != 'all':
            gene_list_temp = gene_list[0:index_cutoff]
        else:
            gene_list_temp = gene_list
        # ----------------------------------------------------


        seed_data_current = regev_library.get_seed_data(1000, 1000, gene_list_temp, total_data, condition_key, mode = 'name')
        cross_data_current = regev_library.create_cross_data(seed_data_current, total_data, condition_list)


        current_species_values = {}

        for gene_counter, gene in enumerate(gene_list_temp):
            score_data[gene] = {}
            score_data[gene][current_species] = {}
            print '------------------'
            print 'Species ' + str(species_counter + 1) + ' of ' + str(len(total_data))
            print 'Gene ' + str(gene_counter + 1) + ' of ' + str(len(gene_list_temp))
            print gene
            print '------------------'
            current_species_values[condition] = seed_data_current['Gene Map'][gene]
            for test_species in total_data:
                test_species_values = {}
                for test_condition in total_data[test_species]:
                    test_condition_key = test_species + ':' + test_condition
                    if test_condition_key != condition_key:
                        if test_species == current_species:
                            current_species_values[test_condition] = cross_data_current[test_condition_key]['Gene Map'][gene]
                            test_species_values[test_condition] = [cross_data_current[test_condition_key]['Gene Map'][gene]]
                        else:
                            test_species_values[test_condition] = []
                            for orth_gene in cross_data_current[test_condition_key]['Seed Gene Map'][gene]:
                                # print test_condition_key, orth_gene, cross_data_current[test_condition_key]['Gene Map'][orth_gene]

                                if orth_gene in cross_data_current[test_condition_key]['Gene Map']:
                                    test_species_values[test_condition].append(cross_data_current[test_condition_key]['Gene Map'][orth_gene])
                                else:
                                    test_species_values[test_condition].append(float('nan'))
                # print ' '
                # print ' '
                # print test_species
                # print test_species_values
                # print ' '

                current_common_values = []
                test_common_values = []

                for current_condition in current_species_values:
                    for test_condition in test_species_values:
                        if current_condition == test_condition:
                        # if current_condition in test_condition or test_condition in current_condition:
                            current_common_values.append(current_species_values[current_condition])
                            best_distance = float('inf')
                            best_value = float('nan')

                            for test_value in test_species_values[test_condition]:
                                if np.isnan(test_value):
                                    continue
                                elif abs(test_value - current_species_values[current_condition]) < best_distance:
                                    best_distance = abs(test_value - current_species_values[current_condition])
                                    best_value = test_value
                            test_common_values.append(best_value)


                nan_indices = []
                for index, value in enumerate(test_common_values):
                    if np.isnan(value) or np.isnan(current_common_values[index]):
                        nan_indices.append(index)

                for counter, nan_index in enumerate(nan_indices):
                    del test_common_values[nan_index - counter]
                    del current_common_values[nan_index - counter]

                # print ' '
                # print current_common_values
                # print ' '
                # print test_common_values
                # print ' '
                # print len(current_common_values) == len(test_common_values)
                if len(current_common_values) > 1 and is_uniform(current_common_values):
                    current_common_values[0] = current_common_values[0] + 1e-8
                if len(test_common_values) > 1 and is_uniform(test_common_values):
                    test_common_values[0] = test_common_values[0] + 1e-8


                if len(current_common_values) == 0:
                    score_data[gene][current_species][test_species] = 0.0
                elif len(current_common_values) == 1:
                    # I don't like this metric of setting the single-point coeficients to zero,
                    # but I can't think of anything else for now
                    score_data[gene][current_species][test_species] = 0.0
                elif current_common_values == test_common_values:
                    score_data[gene][current_species][test_species] = 1.0
                else:
                    coeff, p_val = scipy.stats.pearsonr(current_common_values, test_common_values)
                    score_data[gene][current_species][test_species] = abs(coeff)
                    
                if np.isnan(score_data[gene][current_species][test_species]):
                    print '*************************************************************'
                    print current_common_values
                    print test_common_values
                    print '*************************************************************'

            
            # print ' '
            # print ' '
            # print 'HOME SPECIES'
            # print current_species
            # print current_species_values
            # print ' '

    # Add in the hanging genes from susan/oshea data
    current_species = 'Saccharomyces cerevisiae'
    condition = 'susan'
    condition_key = current_species + ':' + condition
    # if species == 'Saccharomyces cerevisiae':
    ORF_gene_table, gene_ORF_table = io_library.read_SGD_features()
    print '--------------------------------------------------'
    print condition_key
    print '--------------------------------------------------'
    seed_data_current = regev_library.get_seed_data(1000, 1000, gene_list_susan, total_data, condition_key, mode = 'name')
    cross_data_current = regev_library.create_cross_data(seed_data_current, total_data, condition_list)
    current_species_values = {}
    for gene_counter, gene in enumerate(gene_list_susan):
        score_data[gene] = {}
        score_data[gene][current_species] = {}
        print '------------------'
        # print 'Species ' + str(species_counter + 1) + ' of ' + str(len(total_data))
        print 'Adding Susan Genes'
        print 'Gene ' + str(gene_counter + 1) + ' of ' + str(len(gene_list_susan))
        print gene
        print '------------------'
        current_species_values[condition] = seed_data_current['Gene Map'][gene]
        for test_species in total_data:
            test_species_values = {}
            for test_condition in total_data[test_species]:
                test_condition_key = test_species + ':' + test_condition
                if test_condition_key != condition_key:
                    if test_species == current_species:
                        current_species_values[test_condition] = cross_data_current[test_condition_key]['Gene Map'][gene]
                    else:
                        test_species_values[test_condition] = []
                        for orth_gene in cross_data_current[test_condition_key]['Seed Gene Map'][gene]:
                            # print test_condition_key, orth_gene, cross_data_current[test_condition_key]['Gene Map'][orth_gene]

                            if orth_gene in cross_data_current[test_condition_key]['Gene Map']:
                                test_species_values[test_condition].append(cross_data_current[test_condition_key]['Gene Map'][orth_gene])
                            else:
                                test_species_values[test_condition].append(float('nan'))
            # print ' '
            # print ' '
            # print test_species
            # print test_species_values
            # print ' '

            current_common_values = []
            test_common_values = []

            for current_condition in current_species_values:
                for test_condition in test_species_values:
                    if current_condition == test_condition:
                    # if current_condition in test_condition or test_condition in current_condition:
                        current_common_values.append(current_species_values[current_condition])
                        best_distance = float('inf')
                        best_value = float('nan')

                        for test_value in test_species_values[test_condition]:
                            if np.isnan(test_value):
                                continue
                            elif abs(test_value - current_species_values[current_condition]) < best_distance:
                                best_distance = abs(test_value - current_species_values[current_condition])
                                best_value = test_value
                        test_common_values.append(best_value)


            nan_indices = []
            for index, value in enumerate(test_common_values):
                if np.isnan(value) or np.isnan(current_common_values[index]):
                    nan_indices.append(index)


            for counter, nan_index in enumerate(nan_indices):
                del test_common_values[nan_index - counter]
                del current_common_values[nan_index - counter]
                # print test_common_values
                # print current_common_values

            # print ' '
            # print current_common_values
            # print ' '
            # print test_common_values
            # print ' '
            # print len(current_common_values) == len(test_common_values)


            if len(current_common_values) > 1 and is_uniform(current_common_values):
                current_common_values[0] = current_common_values[0] + 1e-8
            if len(test_common_values) > 1 and is_uniform(test_common_values):
                test_common_values[0] = test_common_values[0] + 1e-8

            if len(current_common_values) == 0:
                if test_species == 'Saccharomyces cerevisiae': # for the hanging genes
                    score_data[gene][current_species][test_species] = 1.0
                else:
                    score_data[gene][current_species][test_species] = 0.0
            elif len(current_common_values) == 1:
                # I don't like this metric of setting the single-point coeficients to zero,
                # but I can't think of anything else for now
                if test_species == 'Saccharomyces cerevisiae': # for the hanging genes
                    score_data[gene][current_species][test_species] = 1.0
                else:
                    score_data[gene][current_species][test_species] = 0.0
            elif current_common_values == test_common_values:
                score_data[gene][current_species][test_species] = 1.0
            else:
                coeff, p_val = scipy.stats.pearsonr(current_common_values, test_common_values)
                score_data[gene][current_species][test_species] = abs(coeff)

            if np.isnan(score_data[gene][current_species][test_species]):
                print '*************************************************************'
                print current_common_values
                print test_common_values
                print '*************************************************************'


    np.save('./score_data_' + str(index_cutoff), score_data)
    return score_data

# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
def plot_scores(gene_list, score_data):
    """
    Plots the cross-species score data pertaining to a given set of genes

    Args:
        gene_list: The set of genes of interest
        score_data: Total score_data structure to grab plotting data from

    Returns:
        None
    """

    num_columns = len(params.species_list)
    num_rows = len(gene_list)
    plot_matrix = np.zeros((num_rows, num_columns))

    for row_index, gene in enumerate(gene_list):
        home_species = score_data[gene].keys()[0]
        for column_index, species in enumerate(params.species_list):
            plot_matrix[row_index, column_index] = score_data[gene][home_species][species]


    plotting_labels_x = params.species_list
    # plotting_labels_y = gene_list
    plotting_labels_y = [' ']

    # cmap = mpl.cm.RdBu_r
    cmap = mpl.cm.Blues
    cmap.set_bad('k',1.0)
    fig, ax = plt.subplots()
    ax = sns.heatmap(plot_matrix, cmap = cmap, ax = ax, yticklabels = plotting_labels_y, xticklabels = plotting_labels_x)
    ylocs, ylabels = plt.yticks()
    plt.setp(ylabels, rotation = 0)
    xlocs, xlabels = plt.xticks()
    plt.setp(xlabels, rotation = 90)
    # title = plt.title('Seed Condition ' + seed_data['Condition Key'])
    # title = plt.title('Seed Condition ' + params.condition_key + '_' + str(params.HIGH_NUM) + '_' + str(params.LOW_NUM))
    # title.set_y(1.05)
    # plt.subplots_adjust(top=0.86)
    # title = 'Seed Condition ' + seed_data['Condition Key']
    # plt.text(0.5, 2, title, horizontalalignment = 'center', fontsize = 12)
    plt.tight_layout()
    plt.get_current_fig_manager().window.raise_()
    plt.show()

# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------

def get_worst_genes(score_data, number, test_species, root_species):
    """
    Given the total score data and two species, a test and root species, 
    finds the set of worst genes (as determined by number argument) in terms
    of correlation score between the two
    
    Args:
        score_data: Total score data structure
        number: The number of worst genes to find
        test_species: The species in question (the correlation scores will come from this
                      species)
        root_species: The 'home' species (the gene list will come from
                      this species)
    
    Returns:
        gene_list_worst: A list of the worst correlated genes between the two species. This list 
                         will be the gene names from the root_species
    """
    gene_list = []

    for gene in score_data:
        # print '------------------------------------------------------------------------------'
        # print score_data[gene][test_species][root_species]
        # print score_data[gene][root_species][test_species]
        if score_data[gene].keys() == [test_species]:
            # print gene
            # print score_data[gene][test_species][root_species]
            gene_list.append((gene, score_data[gene][test_species][root_species]))
        # print '------------------------------------------------------------------------------'

    # Sort the gene list by values
    gene_list.sort(key = lambda tup : tup[1])

    gene_list_worst = []
    index = 0
    while len(gene_list_worst) < number:
        if gene_list[index][1] != 0.0:
            gene_list_worst.append(gene_list[index])
        index += 1


    return gene_list_worst

# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------


def get_susan_rank(gene_tuple_list, species, total_data):
    """
    Finds the pka rank of the genes in the gene_tuple_list

    Args:
        gene_tuple_list: List of tuples with format
                                (gene, value)
                         for which we want a correpsonding 
                         pka activation/inhibition rank
        species: The species from which the gene_tuple_list comes from
        total_data: Total data struct

    Returns: List of rank tuples of form 
                    (susan rank, oshea rank)
             Each rank is a measure of how much the gene's ortholog in SCer changes in 
             pka inhibition (in terms of absolute value of logfold chnage of gene expression)
    """
    rank_list = []
    # print '---------------------------------------------------------------------'
    # print species

    sorted_susan_tuples = total_data['Saccharomyces cerevisiae']['susan']['Tuples']
    sorted_susan_tuples.sort(key = lambda tup : abs(tup[1]), reverse = True)

    sorted_oshea_tuples = total_data['Saccharomyces cerevisiae']['oshea']['Tuples']
    sorted_oshea_tuples.sort(key = lambda tup : abs(tup[1]), reverse = True)

    orth_table = io_library.read_orth_lookup_table(species, 'Saccharomyces cerevisiae')
    # for gene_tuple in gene_tuple_list:
    for gene_tuple in gene_tuple_list:
        gene = gene_tuple[0]
        value = gene_tuple[1]
        orth = orth_table[gene][0]
        susan_rank = -1
        found_susan = False
        oshea_rank = -1
        found_oshea = False
        # print orth

        for susan_index, susan_tuple in enumerate(sorted_susan_tuples):
            susan_rank += 1
            # print susan_tuple[0]
            if susan_tuple[0] == orth:
                found_susan = True
                break
        # print susan_rank, found_susan
        # if not found_susan:
            # susan_rank = -1

        for oshea_index, oshea_tuple in enumerate(sorted_oshea_tuples):
            oshea_rank += 1
            # print susan_tuple[0]
            if oshea_tuple[0] == orth:
                found_oshea = True
                break
        # print oshea_rank, found_oshea
        # if not found_oshea:
            # oshea_rank = -1

        rank_list.append((susan_rank, oshea_rank))

    # print '---------------------------------------------------------------------'


    return rank_list

# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------