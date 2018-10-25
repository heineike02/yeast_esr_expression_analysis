# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
import regev_library
import params
import single_score_library
import species_score_library
import numpy as np
import os
import matplotlib as mpl
import seaborn as sns
import matplotlib.pyplot as plt
base_dir = os.path.normpath(os.path.dirname(os.getcwd()))
import sys
sys.path.append(base_dir + '/core')
import io_library

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# This is a mix of scripts that I was running to try to get results. They're generally broken up into the individual results I was trying to generate
# Some of it may be useful to look at for examples on what to do in the future
# 
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create the pka (or single condition) comparison score matrix
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# seed_data, cross_data, plotting_data, plotting_labels_x, plotting_labels_y = single_score_library.check_and_create(params.HIGH_NUM, params.LOW_NUM, params.condition_key, params.ext)

# score_data = single_score_library.get_scores(seed_data, cross_data)

# # temp_list = [x / 17.6635664592 for x in score_data[regev_library.create_condition_key('SCer', 'oshea')]['Values'] ]
# # print temp_list

# score_matrix = single_score_library.create_score_matrix(plotting_data, score_data, seed_data)

# # print score_matrix[ : , 1] - temp_list
# single_score_library.plot_score_matrix(score_matrix, plotting_labels_x, plotting_labels_y, seed_data)

# # regev_library.plot_data(seed_data, plotting_data, plotting_labels_x, plotting_labels_y)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Creat the species (Scer) comparison score matrix
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# total_data, condition_list = regev_library.compile_total_data()
# # # seed_data, cross_data, plotting_data, plotting_labels_x, plotting_labels_y = single_score_library.check_and_create(params.HIGH_NUM, params.LOW_NUM, params.condition_key, params.ext)
# # score_data = species_score_library.get_scores(total_data, condition_list)


# cwd = os.getcwd()

# score_data = np.load(cwd + '/stored_data/score_data_all.npy')[()]

# # print score_data.keys()

# gene_list = total_data['Saccharomyces cerevisiae']['susan']['Genes'][0:300]
# species_score_library.plot_scores(gene_list, score_data)



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Get the least correlated genes in each species, check their pka inhibition rank (absolute value of logfold change in susan and oshea)
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# CHECK_GENES = 300
# cwd = os.getcwd()
# score_data = np.load(cwd + '/stored_data/score_data_all.npy')[()]
# total_data = np.load(cwd + '/stored_data/total_data/total_data.npy')[()]
# seed_data, cross_data, plotting_data, plotting_labels_x, plotting_labels_y = single_score_library.check_and_create(params.HIGH_NUM, params.LOW_NUM, params.condition_key, params.ext)

# worst_genes = {}
# ranks = {}
# # temp_species_list = ['Kluyveromyces lactis']
# # for species in temp_species_list:
# # print total_data['Saccharomyces cerevisiae']['oshea']['Value Map']


# species_list_test = ['Saccharomyces bayanus', 
# 				'Candida glabrata', 'Saccharomyces castellii',
# 				'Kluyveromyces lactis']

# # sorted_susan_tuples = total_data['Saccharomyces cerevisiae']['susan']['Tuples']
# # sorted_susan_tuples.sort(key = lambda tup : abs(tup[1]), reverse = True)
# # print sorted_susan_tuples

# for species in species_list_test:
# 	worst_genes[species] = species_score_library.get_worst_genes(score_data, CHECK_GENES, species, 'Saccharomyces cerevisiae')
# 	ranks[species] = species_score_library.get_susan_rank(worst_genes[species], species, total_data)

# 	# Don't need to do this, but doing it for visualization. If the rank is low enough for either dataset, hard code it to 0
# 	# for rank_index, rank_tuple in enumerate(ranks[species]):
# 	# 	if rank_tuple[0] < CHECK_GENES:
# 	# 		ranks[species][rank_index] = (0, rank_tuple[1])
# 	# 	if rank_tuple[1] < CHECK_GENES:
# 	# 		ranks[species][rank_index] = (rank_tuple[0], 0)
# 	# 	if rank_tuple[1] > CHECK_GENES and rank_tuple[0] > CHECK_GENES:
# 	# 		ranks[species][rank_index] = (max(max(ranks[species])), max(max(ranks[species])))




# # ------------------------------------------------------------------------------------
# #  Plotting
# # ------------------------------------------------------------------------------------

# num_columns = len(species_list_test) * 2
# num_rows = CHECK_GENES
# plot_matrix = np.zeros((num_rows, num_columns))

# # for row_index, gene in enumerate(gene_list):
#     # home_species = score_data[gene].keys()[0]
#     # for column_index, species in enumerate(params.species_list):
#         # plot_matrix[row_index, column_index] = score_data[gene][home_species][species]


# column_index = 0
# for species in species_list_test:
# 	for row_index, rank_tuple in enumerate(ranks[species]):
# 		for tuple_index in range(2):
# 			plot_matrix[row_index, column_index + tuple_index] = rank_tuple[tuple_index]
# 	column_index += 2

# plotting_labels_x = []
# for species in species_list_test:
# 	for i in range(2):
# 		if i == 0:
# 			# plotting_labels_x.append(species + '_susan')
# 			plotting_labels_x.append(params.species_name_dict[species] + '_susan')
# 		elif i == 1:
# 			# plotting_labels_x.append(species + '_oshea')
# 			plotting_labels_x.append(params.species_name_dict[species] + '_oshea')
# plotting_labels_y = [' ']

# # plotting_labels_x = params.species_list
# # plotting_labels_y = gene_list

# # cmap = mpl.cm.RdBu_r
# cmap = mpl.cm.Blues
# cmap.set_bad('k',1.0)
# fig, ax = plt.subplots()
# ax = sns.heatmap(plot_matrix, cmap = cmap, ax = ax, xticklabels = plotting_labels_x, yticklabels = plotting_labels_y)
# # ylocs, ylabels = plt.yticks()
# # plt.setp(ylabels, rotation = 0)
# # xlocs, xlabels = plt.xticks()
# # plt.setp(xlabels, rotation = 90)


# plt.tight_layout()
# plt.get_current_fig_manager().window.raise_()
# plt.show()


# # ------------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------------

# # print worst_genes



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Get the top genes from Susan's (or Oshea's) data and check their scores across species
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
CHECK_GENES = 300
pka_condition = 'susan'
cwd = os.getcwd()
score_data = np.load(cwd + '/stored_data/score_data_all.npy')[()]
total_data = np.load(cwd + '/stored_data/total_data/total_data.npy')[()]
seed_data, cross_data, plotting_data, plotting_labels_x, plotting_labels_y = single_score_library.check_and_create(params.HIGH_NUM, params.LOW_NUM, params.condition_key, params.ext)

high_genes = []
scores = {}
home_scores = {}
# temp_species_list = ['Kluyveromyces lactis']
# for species in temp_species_list:
# print total_data['Saccharomyces cerevisiae']['oshea']['Value Map']


species_list_test = ['Saccharomyces bayanus', 
				'Candida glabrata', 'Saccharomyces castellii',
				'Kluyveromyces lactis']

sorted_pka_tuples = total_data['Saccharomyces cerevisiae'][pka_condition]['Tuples']
sorted_pka_tuples.sort(key = lambda tup : abs(tup[1]), reverse = True)
# print sorted_pka_tuples
sorted_pka_tuples = sorted_pka_tuples[0 : CHECK_GENES]

for species in species_list_test:
	scores[species] = []
	home_scores[species] = []
	orth_table_forward = io_library.read_orth_lookup_table('Saccharomyces cerevisiae', species)
	orth_table_backward = io_library.read_orth_lookup_table(species, 'Saccharomyces cerevisiae')
	for gene, value in sorted_pka_tuples:
		# best_val = float('-inf')
		best_val = float('-inf')
		best_orth = 'NONE'
		if gene in orth_table_forward:
			for orth in orth_table_forward[gene]:
				if orth not in score_data:
					best_val = float('-inf')
					# best_val = float('nan')
					best_orth = 'NONE'
				elif orth == 'NONE':
					best_val = float('-inf')
					# best_val = float('nan')
					best_orth = 'NONE'
				elif score_data[orth][species]['Saccharomyces cerevisiae'] > best_val:
					best_val = score_data[orth][species]['Saccharomyces cerevisiae']
					best_orth = orth
					# if orth not in orth_table_backward:
					# 	print 'HERE!'

		scores[species].append((best_orth, best_val))
		home_scores[species].append((gene, best_val))



		

	# Don't need to do this, but doing it for visualization. If the rank is low enough for either dataset, hard code it to 0
	# for rank_index, rank_tuple in enumerate(ranks[species]):
	# 	if rank_tuple[0] < CHECK_GENES:
	# 		ranks[species][rank_index] = (0, rank_tuple[1])
	# 	if rank_tuple[1] < CHECK_GENES:
	# 		ranks[species][rank_index] = (rank_tuple[0], 0)
	# 	if rank_tuple[1] > CHECK_GENES and rank_tuple[0] > CHECK_GENES:
	# 		ranks[species][rank_index] = (max(max(ranks[species])), max(max(ranks[species])))




# ------------------------------------------------------------------------------------
#  Plotting
# ------------------------------------------------------------------------------------

# num_columns = len(species_list_test)
# num_rows = CHECK_GENES
# plot_matrix = np.zeros((num_rows, num_columns))

# # for row_index, gene in enumerate(gene_list):
#     # home_species = score_data[gene].keys()[0]
#     # for column_index, species in enumerate(params.species_list):
#         # plot_matrix[row_index, column_index] = score_data[gene][home_species][species]


# column_index = 0
# for species in species_list_test:
# 	for row_index, rank_tuple in enumerate(scores[species]):
# 		if np.isinf(rank_tuple[1]):
# 			plot_matrix[row_index, column_index] = float('nan')
# 		else:	
# 			plot_matrix[row_index, column_index] = rank_tuple[1]
# 	column_index += 1



# score_matrix = np.copy(plot_matrix)
# # sort matrix by the K Lactis column
# score_matrix = score_matrix[score_matrix[: , 3].argsort()]


# plotting_labels_x = []
# for species in species_list_test:
# 			plotting_labels_x.append(species)


# # plotting_labels_y = [' ']

# # plotting_labels_y = np.zeros(CHECK_GENES)
# plotting_labels_y = []
# # for gene, value in sorted_pka_tuples:
# 	# plotting_labels_y.append(gene) 

# # plotting_labels_x = params.species_list
# # plotting_labels_y = gene_list

# # cmap = mpl.cm.RdBu_r
# cmap = mpl.cm.Blues
# cmap.set_bad('k',1.0)
# fig, ax = plt.subplots()
# ax = sns.heatmap(score_matrix, cmap = cmap, ax = ax, xticklabels = plotting_labels_x, yticklabels = plotting_labels_y)
# # ylocs, ylabels = plt.yticks()
# # plt.setp(ylabels, rotation = 0)
# # xlocs, xlabels = plt.xticks()
# # plt.setp(xlabels, rotation = 90)


# plt.tight_layout()
# plt.get_current_fig_manager().window.raise_()
# plt.show()


# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------
# Sort the K Lac values and get genes
# ------------------------------------------------------------------------------------

klac_scores = np.copy(scores['Kluyveromyces lactis'])
klac_home_scores = np.copy(home_scores['Kluyveromyces lactis'])


klac_scores_sorted = sorted(klac_scores, key = lambda tup : tup[1])
klac_home_scores_sorted = sorted(klac_home_scores, key = lambda tup : tup[1])

orth_table = io_library.read_orth_lookup_table('Kluyveromyces lactis', 'Saccharomyces cerevisiae')

del_indices = []
for index, tup in enumerate(klac_scores_sorted):
	# print '---------------------------'
	# print tup
	if tup[0] == 'NONE':
		del_indices.append(index)
		# if not np.isinf(float(tup[1])):
			# del_indices.append(index)
		# else:
			# klac_scores_sorted[index] = (tup[0], 0.0)

		if np.isinf(float(tup[1])):
			klac_home_scores_sorted[index] = (klac_home_scores_sorted[index][0], 0.0)


for counter, index in enumerate(del_indices):
	del klac_scores_sorted[index - counter]


num_genes = len(klac_scores_sorted)
# for index, tup in enumerate(klac_home_scores_sorted):
	# print tup[0], tup[1], index

klac_home_scores_lowest = klac_home_scores_sorted[0:30]
# for index, tup in enumerate(klac_home_scores_lowest):
	# print tup[0]

# for gene, value in klac_scores_sorted:
# 	print gene, value, orth_table[gene]

# print num_genes

# print klac_scores_sorted[0:90]
# print '---------------------------------------------'
# print klac_scores_sorted[90:180]
# print '---------------------------------------------'
# print klac_scores_sorted[180:270]

low_score_genes = []
middle_score_genes = []
high_score_genes = []

for gene, value in klac_scores_sorted[0:90]:
	low_score_genes.append(gene)
for gene, value in klac_scores_sorted[90:180]:
	middle_score_genes.append(gene)
for gene, value in klac_scores_sorted[180:270]:
	high_score_genes.append(gene)


# orth_table_temp = io_library.read_orth_lookup_table('Saccharomyces cerevisiae', 'Kluyveromyces lactis')
# print orth_table_temp['YBR285W']



# for gene in low_score_genes:
# 	for orth in orth_table[gene]:
# 		if orth in total_data['Saccharomyces cerevisiae'][pka_condition]['Genes']:
# 			print orth



# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------





# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
