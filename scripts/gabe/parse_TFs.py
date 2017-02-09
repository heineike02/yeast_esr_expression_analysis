# ----------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------
import sys
sys.path.append('/Users/greder/Dropbox/UCSF/El-Samad/expression_broad_data/scripts')
sys.path.append('/Users/greder/Dropbox/UCSF/El-Samad/expression_broad_data/core')
import regev_library
import os
import csv
import scipy.stats
import params
import scipy.cluster.hierarchy as cluster
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
import matplotlib as mpl
cwd = os.getcwd()
cwd = cwd + '/stored_data/text_data'
# ----------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------
# 
# This script is more of a single workflow than cross_species_GO, though at the end it started piling up. Again, this code might be
# pretty reusable. This is the workflow for getting the TF p-vals. The text and csv files are stored in stored_data/text_files. New data
# can be generated or regenereated from downloading the yeastract csv files. 
# 
# ----------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------
# Change SELECTION to either run analysis on the scer pka genes whose orthologs have low correlation
# scores (but the orthologs exist) or the scer pka genes that have no klac orthologs
# ----------------------------------------------------------------------------------------------------

# SELECTION = 'LOW_CORR'
SELECTION = 'NO_ORTH'

# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------------------
# 
# Goes through the yeastract csv files and compiles the TF data. 
# 
# ----------------------------------------------------------------------------------------------------
if SELECTION == 'LOW_CORR':
	fname = os.path.normpath(cwd+'/TFs_klac_DNA_and_expression.csv')
else:
	fname = os.path.normpath(cwd+'/TFs_klac_no_orth_AND.csv')


counts = {}
genes_Klac_DNA_And_Expression = {}
TFs_Klac_DNA_And_Expression = {}
csvfile = open(fname, 'rb')
reader = csv.reader(csvfile, delimiter = '\t')
num_klac_genes = 0
for row in reader:
	TF = row[2]
	gene = row[1]
	# print row
	if TF in counts:
		counts[TF] += 1
	else:
		counts[TF] = 1

	if TF in TFs_Klac_DNA_And_Expression:
		TFs_Klac_DNA_And_Expression[TF].append(gene)
	else:
		TFs_Klac_DNA_And_Expression[TF] = [gene]

	if gene in genes_Klac_DNA_And_Expression:
		genes_Klac_DNA_And_Expression[gene].append(TF)
	else:
		num_klac_genes += 1
		genes_Klac_DNA_And_Expression[gene] = [TF]


Klac_TF_DNA_And_Expression_tuples = []
for TF in counts:
	Klac_TF_DNA_And_Expression_tuples.append((TF, counts[TF]))

# ----------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------
if SELECTION == 'LOW_CORR':
	fname = os.path.normpath(cwd+'/TFs_klac_DNA_only.csv')
else:
	fname = os.path.normpath(cwd+'/TFs_klac_no_orth_DNA.csv')


counts = {}
genes_Klac_DNA_Only = {}
TFs_Klac_DNA_Only = {}
csvfile = open(fname, 'rb')
reader = csv.reader(csvfile, delimiter = '\t')
num_klac_genes = 0
for row in reader:
	TF = row[2]
	gene = row[1]
	# print row
	if TF in counts:
		counts[TF] += 1
	else:
		counts[TF] = 1

	if TF in TFs_Klac_DNA_Only:
		TFs_Klac_DNA_Only[TF].append(gene)
	else:
		TFs_Klac_DNA_Only[TF] = [gene]

	if gene in genes_Klac_DNA_Only:
		genes_Klac_DNA_Only[gene].append(TF)
	else:
		num_klac_genes += 1
		genes_Klac_DNA_Only[gene] = [TF]


Klac_TF_DNA_Only_tuples = []
for TF in counts:
	Klac_TF_DNA_Only_tuples.append((TF, counts[TF]))

# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
if SELECTION == 'LOW_CORR':
	fname = os.path.normpath(cwd+'/TFs_klac_DNA_plus_Expression.csv')
else:
	fname = os.path.normpath(cwd+'/TFs_klac_no_orth_PLUS.csv')

counts = {}
genes_Klac_DNA_Plus_Expression = {}
TFs_Klac_DNA_Plus_Expression = {}
csvfile = open(fname, 'rb')
reader = csv.reader(csvfile, delimiter = '\t')
num_klac_genes = 0
for row in reader:
	TF = row[2]
	gene = row[1]
	# print row
	if TF in counts:
		counts[TF] += 1
	else:
		counts[TF] = 1

	if TF in TFs_Klac_DNA_Plus_Expression:
		TFs_Klac_DNA_Plus_Expression[TF].append(gene)
	else:
		TFs_Klac_DNA_Plus_Expression[TF] = [gene]

	if gene in genes_Klac_DNA_Plus_Expression:
		genes_Klac_DNA_Plus_Expression[gene].append(TF)
	else:
		num_klac_genes += 1
		genes_Klac_DNA_Plus_Expression[gene] = [TF]


Klac_TF_DNA_Plus_Expression_tuples = []
for TF in counts:
	Klac_TF_DNA_Plus_Expression_tuples.append((TF, counts[TF]))

# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
if SELECTION == 'LOW_CORR':
	fname = os.path.normpath(cwd+'/TFs_klac_expression_only.csv')
else:
	fname = os.path.normpath(cwd+'/TFs_klac_no_orth_EXPRESSION.csv')


counts = {}
genes_Klac_Expression_Only = {}
TFs_Klac_Expression_Only = {}
csvfile = open(fname, 'rb')
reader = csv.reader(csvfile, delimiter = '\t')
num_klac_genes = 0
for row in reader:
	TF = row[2]
	gene = row[1]
	# print row
	if TF in counts:
		counts[TF] += 1
	else:
		counts[TF] = 1

	if TF in TFs_Klac_Expression_Only:
		TFs_Klac_Expression_Only[TF].append(gene)
	else:
		TFs_Klac_Expression_Only[TF] = [gene]

	if gene in genes_Klac_Expression_Only:
		genes_Klac_Expression_Only[gene].append(TF)
	else:
		num_klac_genes += 1
		genes_Klac_Expression_Only[gene] = [TF]


Klac_TF_Expression_Only_tuples = []
for TF in counts:
	Klac_TF_Expression_Only_tuples.append((TF, counts[TF]))

# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
fname = os.path.normpath(cwd+'/TFs_pka.csv')
counts = {}
genes_pka = {}
TFs_pka = {}
csvfile = open(fname, 'rb')
reader = csv.reader(csvfile, delimiter = '\t')
num_pka_genes = 0
for row in reader:
	TF = row[2]
	gene = row[1]
	if TF in counts:
		counts[TF] += 1
	else:
		counts[TF] = 1

	if TF in TFs_pka:
		TFs_pka[TF].append(gene)
	else:
		TFs_pka[TF] = [gene]

	if gene in genes_pka:
		genes_pka[gene].append(TF)
	else:
		num_pka_genes += 1
		genes_pka[gene] = [TF]

pka_TF_tuples = []
for TF in counts:
	pka_TF_tuples.append((TF, counts[TF]))

# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
	
# print num_pka_genes
# print len(Klac_TF_tuples)
# print len(pka_TF_tuples)

pka_dict = {}
for tup in pka_TF_tuples:
	TF = tup[0]
	if TF not in pka_dict:
		pka_dict[TF] = tup[1]



TF_dict = {}
TF_search_list = []
for tup in Klac_TF_DNA_Plus_Expression_tuples:
	TF = tup[0]
	p = pka_dict[TF] / float(num_pka_genes)
	DNA_Plus_Expression_observation = tup[1]
	TF_dict[TF] = {}
	TF_dict[TF]['p'] = p
	TF_dict[TF]['pka'] = pka_dict[TF]
	TF_dict[TF]['PLUS'] = DNA_Plus_Expression_observation

for tup in Klac_TF_DNA_And_Expression_tuples:
	TF = tup[0]
	TF_dict[TF]['AND'] = tup[1]


for tup in Klac_TF_DNA_Only_tuples:
	TF = tup[0]
	TF_dict[TF]['DNA'] = tup[1]

for tup in Klac_TF_Expression_Only_tuples:
	TF = tup[0]
	TF_dict[TF]['EXPRESSION'] = tup[1]

TF_tuples = []
for TF in TF_dict:
	p = TF_dict[TF]['p']
	PLUS_observation = TF_dict[TF]['PLUS']
	PLUS_genes = TFs_Klac_DNA_Plus_Expression[TF]
	PLUS_p = scipy.stats.binom_test(PLUS_observation, n = float(num_klac_genes), p = p) 

	if 'AND' in TF_dict[TF].keys():
		AND_observation = TF_dict[TF]['AND']
		AND_genes = TFs_Klac_DNA_And_Expression[TF]
		AND_p = scipy.stats.binom_test(AND_observation, n = float(num_klac_genes), p = p)

	else:
		AND_observation = 0
		AND_genes = []
		AND_p = scipy.stats.binom_test(AND_observation, n = float(num_klac_genes), p = p)

	if 'DNA' in TF_dict[TF].keys():
		DNA_observation = TF_dict[TF]['DNA']
		DNA_genes = TFs_Klac_DNA_Only[TF]
		DNA_p = scipy.stats.binom_test(DNA_observation, n = float(num_klac_genes), p = p)
	else:
		DNA_observation = 0
		DNA_genes = []
		DNA_p = scipy.stats.binom_test(DNA_observation, n = float(num_klac_genes), p = p)

	if 'EXPRESSION' in TF_dict[TF].keys():
		EXPRESSION_observation = TF_dict[TF]['EXPRESSION']
		EXPRESSION_genes = TFs_Klac_Expression_Only[TF]
		EXPRESSION_p = scipy.stats.binom_test(EXPRESSION_observation, n = float(num_klac_genes), p = p)
	else:
		EXPRESSION_observation = 0
		EXPRESSION_genes = []
		EXPRESSION_p = scipy.stats.binom_test(EXPRESSION_observation, n = float(num_klac_genes), p = p)


	tup = (TF, TF_dict[TF]['pka'], p, PLUS_observation, PLUS_p, AND_observation, AND_p, DNA_observation, DNA_p, EXPRESSION_observation, EXPRESSION_p, PLUS_genes)
	TF_tuples.append(tup)

TF_index = 0
pka_count_index = 1
pka_p_index = 2
PLUS_observation_index = 3
PLUS_p_index = 4
AND_observation_index = 5
AND_p_index = 6
DNA_observation_index = 7
DNA_p_index = 8
EXPRESSION_observation_index = 9
EXPRESSION_p_index = 10
PLUS_genes_index = 11

TF_tuples_sorted = sorted(TF_tuples, key = lambda tup : tup[PLUS_p_index])

# for tup in TF_tuples_sorted:
	# print tup[TF_index], tup[PLUS_p_index], tup[AND_p_index], tup[DNA_p_index], tup[EXPRESSION_p_index], tup[pka_count_index] / float(num_pka_genes), tup[PLUS_observation_index] / float(num_klac_genes)

top_TFs = []
for TF_tuple in TF_tuples_sorted:
	if max(TF_tuple[PLUS_p_index], TF_tuple[AND_p_index], TF_tuple[DNA_p_index], TF_tuple[EXPRESSION_p_index]) < 0.01:
		top_TFs.append(TF_tuple)

# for TF_tuple in top_TFs:
	# print TF_tuple[TF_index], TF_tuple[2], TF_tuple[3] / float(num_klac_genes), TF_tuple[4]




# ----------------------------------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------------------------------
if SELECTION == 'LOW_CORR':
	fname = os.path.normpath(cwd+'/Klac_low_genes.txt')
else:
	fname = os.path.normpath(cwd+'/klac_no_orth.txt')

gene_list = []
csvfile = open(fname, 'rb')
reader = csv.reader(csvfile, delimiter = '\n')
for row in reader:
	gene = row[0]
	gene_list.append(gene)


# gene_list = genes_Klac_DNA_Plus_Expression.keys()
# print gene_list

total_data, condition_list = regev_library.compile_total_data()
gene_list_susan = total_data['Saccharomyces cerevisiae']['susan']['Genes']

condition_list_klac = []
for condition in condition_list:
	# if params.get_species(condition) == 'Kluyveromyces lactis' or params.get_species(condition) == 'Saccharomyces cerevisiae':
	if params.get_species(condition) == 'Saccharomyces cerevisiae':
		condition_list_klac.append(condition)

# for condition in condition_list:
	# print condition
# condition_key = 'Kluyveromyces lactis:Depletion_avg'
condition_key = 'Saccharomyces cerevisiae:susan'
data = regev_library.data_by_key(condition_key, total_data)
seed_data = regev_library.get_seed_data(0, 0, gene_list, total_data, condition_key, mode = 'name')
# print len(seed_data['Genes'])
# print len(gene_list)
# print seed_data

cross_data = regev_library.create_cross_data(seed_data, total_data, condition_list_klac)
plotting_data, plotting_labels_x, plotting_labels_y = regev_library.create_plotting_data(seed_data, cross_data)

# regev_library.plot_data(seed_data, plotting_data, plotting_labels_x, plotting_labels_y)


# mask = np.zeros((plotting_data.shape), dtype = bool)
for row_index in range(len(plotting_data[: , 0])):
	for col_index in range(len(plotting_data[0, : ])):
		if np.isnan(plotting_data[row_index, col_index]):
			# plotting_data[row_index, col_index] = 0.0
			plotting_data[row_index, col_index] = np.average(plotting_data[ row_index , : ])
			# mask[row_index, col_index] = True






# ----------------------------------------------------------------------------------------------------------------
# 
# Hierarchical clustering and storing of the data using seaborn function.
# This is an example of creating TF and clustering matrices for using the regev_library
# plot_data_TF function
# 
# ----------------------------------------------------------------------------------------------------------------

# fig, ax = plt.subplots(1, 1, sharex=True)
# sns.despine(left = True)
clustergrid = sns.clustermap(plotting_data, col_cluster = False, metric = 'chebyshev',
							 xticklabels = plotting_labels_x, yticklabels = plotting_labels_y)
# cluster_matrix = clustergrid.data
# cluster_matrix = cluster_matrix.values.tolist()
# cluster_matrix = np.asarray(cluster_matrix)
ax = clustergrid.ax_heatmap
new_order_gene_list = []
# print len(cluster_matrix.values.tolist()[101])

for index, label in enumerate(ax.get_yticklabels()):
	gene = label.get_text()
	new_order_gene_list.append(gene)

cmap = mpl.cm.RdBu_r
cmap.set_bad('k',1.0)

new_indices = clustergrid.dendrogram_row.reordered_ind
cluster_matrix = np.copy(plotting_data)
for row_index in range(len(cluster_matrix[ : , 0])):
	cluster_matrix[ row_index , : ] = plotting_data[new_indices[row_index]]
# ax = sns.heatmap(cluster_matrix, cmap = cmap, ax = ax, yticklabels = plotting_labels_y, xticklabels = plotting_labels_x)




plotting_labels_y = [' ']



TF_matrix = np.zeros((len(new_order_gene_list), len(top_TFs)))
TF_labels = []
for column_index, TF_tuple in enumerate(top_TFs):
	# print TF_tuple[TF_index]
	TF_labels.append(TF_tuple[TF_index])
	for row_index, gene in enumerate(new_order_gene_list):
		if gene in TF_tuple[PLUS_genes_index]:
			TF_matrix[row_index, column_index] = 1



# print TF_matrix

regev_library.plot_data_TF(seed_data, cluster_matrix, TF_matrix, plotting_labels_x, plotting_labels_y, TF_labels)






# # print num_pka_genes
# # print test_dict
# test_dict_klac = {}
# for tup in Klac_TF_DNA_And_Expression_tuples:
# 	klac_gene = tup[0]
# 	p = test_dict[klac_gene] / float(num_pka_genes)
# 	observation = tup[1]
# 	test_dict_klac[klac_gene] = (klac_gene, observation, p)






# binom_list = []
# for klac_gene in test_dict_klac:
# 	tup = test_dict_klac[klac_gene]
# 	gene = tup[0]
# 	observation = tup[1]
# 	p = tup[2]
# 	# print '------------------'
# 	# print observation
# 	# print p
# 	binom_val = scipy.stats.binom_test(observation, n = float(num_klac_genes), p = p)
# 	# print '------------------'
# 	binom_list.append((gene, binom_val))


# binom_list_sorted = sorted(binom_list, key = lambda tup : tup[1])

# for tup in binom_list_sorted:
# 	print tup


# print counts


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------