# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #

def get_genes_oshea(SC_genes):

	SC_genes_oshea = SC_genes['Oshea'].values
	return SC_genes_oshea



# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #

def compile_species_data(species_list, SC_genes):

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
    SC_gene_plot_data_susan = SC_genes['Log_Change'].values.tolist()
    SC_gene_plot_data_oshea = SC_genes['Oshea'].values.tolist()

    num_columns = 2
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
        row_index += 1
        for i in range(SC_genes_repeats[gene] - 1):
            SC_gene_plot_data_susan.insert(row_index, SC_gene_plot_data_susan[row_index - 1])
            SC_gene_plot_data_oshea.insert(row_index, SC_gene_plot_data_oshea[row_index - 1])
            SC_genes_labels.append(' ')    
            row_index += 1



    plotting_data = np.zeros((num_rows, num_columns))
    plotting_data[:, 0] = SC_gene_plot_data_susan
    plotting_data[:, 1] = SC_gene_plot_data_oshea

    return plotting_data, SC_genes_labels, num_columns

# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #

def plot_data(plotting_data, species_data):
	condition_labels = ['Susan', 'Oshea']
	species_label_indices = []
	col_index = 2
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



# def concat_SC_genes(high_genes, low_genes):
#     # Grab ortholog data from tab file and create lookup tables
#     SC_orfs_lookup, SC_genename_lookup = io_library.read_SGD_features()

#     # Gather all the ORF names for the high genes
#     high_orfs = []
#     for gene_name in high_genes['Gene']:    
#         if gene_name in SC_genename_lookup:
#             high_orfs.append(gene_name)
#         else:
#             high_orfs.append(SC_orfs_lookup[gene_name])


#     # Append the orf names to the existing high_gene data
#     high_gene_data = {'Gene' : high_genes['Gene'], 'Log_Change' : high_genes['Log_Change'],
#                      'ORF' : high_orfs, 'Oshea' : high_genes['Oshea']}
#     high_genes = pd.DataFrame(high_gene_data)

#     # Repeat with the low genes
#     low_orfs = []
#     for gene_name in low_genes['Gene']:    
#         if gene_name in SC_genename_lookup:
#             low_orfs.append(gene_name)
#         else:
#             low_orfs.append(SC_orfs_lookup[gene_name])

#     # Append the orf names to the existing high_gene data
#     low_gene_data = {'Gene' : low_genes['Gene'], 'Log_Change' : low_genes['Log_Change'],
#                      'ORF' : low_orfs, 'Oshea' : low_genes['Oshea']}
#     low_genes = pd.DataFrame(low_gene_data)

#     SC_gene_data = {'Gene' : high_genes['Gene'] + low_genes['Gene'], 
#                     'Log_Change' : high_genes['Log_Change'] + low_genes['Log_Change'],
#                     'ORF' : low_orfs + high_orfs, 'Oshea' : high_genes['Oshea'] + low_genes['Oshea']}


#     SC_genes = pd.concat([high_genes, low_genes])

#     return SC_genes











# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------


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

    seed_genes = seed_data[genes_selector]
    print seed_genes


    for index, condition in enumerate(condition_list):
        if condition != condition_key:
            condition_name = get_condition(condition)
            cross_data[condition] = {}
            species = get_species(condition)
            orth_table_forward = io_library.read_orth_lookup_table(seed_species, species)
            orth_table_backward = io_library.read_orth_lookup_table(species, seed_species)

            seed_genes_orth = []

            for gene in seed_genes:
                ORF_map = seed_data['ORF Map']
                if ORF_map[gene] in orth_table_forward:
                    for orth in orth_table_forward[ORF_map[gene]]:
                        seed_genes_orth.append(orth)
                else:
                    seed_genes_orth.append('NONE')
                    # seed_genes_orth.append(orth_table_forward[ORF_map[gene]])
            # print seed_genes_orth
            # print len(seed_genes_orth)
            
            # print len(total_data[species][condition_name]['Genes']) == len(total_data[species][condition_name]['Values'])
            # print len(total_data['Saccharomyces cerevisiae']['hydrogen peroxide_avg']['Genes']) == len(total_data['Saccharomyces cerevisiae']['hydrogen peroxide_avg']['Values'])

            genes_list = total_data[species][condition_name]['Genes']
            values_list = total_data[species][condition_name]['Values']
            # print condition
            # print species
            # print condition_name
            # print len(genes_list) == len(values_list)

            # genes_list = []
            # values_list = []

            # if len(genes_list) != len(values_list):
            #     print '-----------------------------------'
            #     print condition
            #     print condition_name
            #     print '-----------------------------------'


            # print len(genes_list)
            # print len(values_list)

            indices_to_remove = []
            # print condition
            # print genes_list


            for gene_index, gene in enumerate(genes_list):
                # print species
                # print gene
                # if gene not in orth_table_backward:
                    # indices_to_remove.append(gene_index)
                    # continue
                # elif orth_table_backward[gene] not in seed_data[genes_selector]:
                # elif 'NONE' in orth_table_backward[gene]:
                    # continue
                    # indices_to_remove.append(gene_index)

                # elif gene not in seed_genes_orth:
                if gene not in seed_genes_orth:
                    indices_to_remove.append(gene_index)
                else:
                    continue

            # indices_to_remove = [1]

                #     in_seed_data = False
                #     for lookup_gene in orth_table_backward[gene]:
                #         if lookup_gene in seed_data['ORF Map']:
                #             in_seed_data = True
                # # elif orth_table_backward[gene] not in seed_data['ORF Map']:
                #     if not in_seed_data:
                #         indices_to_remove.append(gene_index)


            # print indices_to_remove

            # indices_to_remove = []


            for index_counter, index in enumerate(indices_to_remove):
                # genes_list.pop(index - index_counter)
                # values_list.pop(index - index_counter)
                # del genes_list[0]
                # del values_list[0]
                genes_list = genes_list[0 : index - index_counter] + genes_list[index - index_counter + 1 : ]
                values_list = values_list[0 : index - index_counter] + values_list[index - index_counter + 1 : ]
                # genes_list = genes_list[1:]
                # values_list = values_list[1:]



                # if len(genes_list) != len(values_list):
                    # print condition + ' ' + str(len(genes_list) == len(values_list)) + ' INDEX ' + str(index_counter)
            




                # genes_list.pop(0)
                # values_list.pop(0)
                # i = 1

            # print 'genes list = ' + str(len(genes_list))
            # print 'values list = ' + str(len(values_list))


            print condition + ' ' + str(len(genes_list))
            cross_data[condition]['Genes'] = genes_list
            cross_data[condition]['Values'] = values_list

            # print len(genes_list) == len(values_list)

            

            
            # print len(genes_list)
            # print len(values_list)
    # print '--------------------------------------------------------'
    # for condition in cross_data:
        # print condition + ' ' + str(len(cross_data[condition]['Genes']))

    return cross_data





# --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------

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
