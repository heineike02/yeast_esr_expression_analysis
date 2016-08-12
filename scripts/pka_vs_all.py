import os
base_dir = os.path.normpath(os.path.dirname(os.getcwd()))
import sys
sys.path.append(base_dir + '/core')
import expression_plots
import io_library
from IPython.core.debugger import Tracer
from IPython.core.debugger import Tracer
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt 
import matplotlib as mpl
import seaborn as sns
import time
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt

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




# Parameters for choosing up and down regulated genes in S.Cer pka inhibition
POSITIVE_LOG_CUTOFF = 8.0
NEGATIVE_LOG_CUTOFF = -6.0

# Now we have a list of the genes in S.Cer whose log changes are high and low respectively in the
# pka experiment
high_genes, low_genes = get_genes_susan(POSITIVE_LOG_CUTOFF, NEGATIVE_LOG_CUTOFF)


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

      
# print low_genes
# SC_orfs_lookup['ARSI-1']

SC_gene_data = {'Gene' : high_genes['Gene'] + low_genes['Gene'], 
                'Log_Change' : high_genes['Log_Change'] + low_genes['Log_Change'],
                'ORF' : low_orfs + high_orfs}


SC_genes = pd.concat([high_genes, low_genes])


# Create ortholog lookup tables for each species


species_list = ['Kluyveromyces lactis', 'Candida glabrata', 
                'Saccharomyces castellii' , 'Saccharomyces bayanus', 'Saccharomyces cerevisiae']

# species_list = ['Kluyveromyces lactis', 'Candida glabrata', 'Saccharomyces bayanus']
# species_list = ['Kluyveromyces lactis', 'Saccharomyces bayanus']

# species_list = ['Candida glabrata']
# species_list = ['Saccharomyces bayanus']
# species_list = ['Saccharomyces castellii']
# species_list = ['Kluyveromyces lactis']


# Refernce dictionary for creating file names
species_name_dict = {'Saccharomyces cerevisiae' : 'SCer',
                    'Kluyveromyces lactis': 'KLac', 
                    'Candida glabrata' : 'CGla', 
                    'Saccharomyces castellii' : 'SCas', 
                    'Saccharomyces bayanus' : 'SBay'}


species_data = {}
SC_genes_repeats = {}
conditions = []
for species in species_list:
    SC_genes_repeats[species] = {}
    species_data[species] = {}
    species_data[species]['condition_names'] = []
    orth_dict = io_library.read_orth_lookup_table('Saccharomyces cerevisiae', species)
#     print "Species: " + species
#     print orth_dict['YPL223C']

    fname = os.path.normpath(base_dir + 
                             '/scripts/expression_broad_data_datafiles/microarray_data/GSE36253_Growth/' + 
                             species_name_dict[species] + '_growth.csv' )
    growth_exp = pd.read_csv(fname,header = [0,1,2], index_col = [0,1])
#     print fname + ' growth microarray dataset loaded'
    #group by conditions and take mean
    growth_replicate_groups = growth_exp.groupby(axis = 1, level = 'conditions')
    growth_exp_avg = growth_replicate_groups.aggregate(np.mean)



    if species != 'Saccharomyces bayanus':
        #There is no stress dataset for S. bayanus
        # Gabe - 7/11/16
        fname = os.path.normpath(base_dir + 
                                 '/scripts/expression_broad_data_datafiles/microarray_data/GSE38478_Stress/' + 
                                 species_name_dict[species] + '_stress.csv' )
        stress_exp = pd.read_csv(fname,header = [0,1,2], index_col = [0,1])
        # Group by condition and take mean
        stress_replicate_groups = stress_exp.groupby(axis = 1, level = 'conditions')
        stress_exp_avg = stress_replicate_groups.aggregate(np.mean)
#         print fname + ' stress microarray dataset loaded'
#         print stress_exp

        #combine growth and stress average expression datasets. 
        if False in stress_exp_avg.index==growth_exp_avg.index:
            print "Error: ID mismatch between condition data. Species = {}".format(species)
            break
        condition_arrays = pd.concat([growth_exp_avg,stress_exp_avg], axis = 1)	
    else:
        condition_arrays = growth_exp_avg



    # print condition_arrays.index

    # print condition_arrays['PS/LOG']
    # print condition_arrays.loc[(slice(None),'KLLA0C01672g'),: ]

    # plotted_genes = condition_arrays.loc[(slice(None),[ortholog[0] for ortholog in native_orf_orthologs]),: ]  
        
    # print condition_arrays['PS/LOG']

    num_conditions = 0
    for condition in condition_arrays:
        species_data[species]['condition_names'].append(condition)
        if condition not in conditions:
            conditions.append(condition)
        species_data[species][condition] = condition_arrays[condition]
        num_conditions += 1

    species_data[species]['num_conditions'] = num_conditions


    


    # species_data[species]['orth'] = orth_dict
    for condition_index, condition in enumerate(condition_arrays):
        species_data[species][condition_index] = []

    species_data[species]['orth_genes'] = []
    for gene in SC_genes['ORF']:
        SC_genes_repeats[species][gene] = 1
        SC_genes_repeats[gene] = 1
        # print '-----------------------------------'
        # print 'S.Cer Gene Name: ' + gene
        if species != 'Saccharomyces cerevisiae':
            try: 
                orth_gene_list = orth_dict[gene]
            except KeyError:
                orth_gene_list = ['NONE']
                # continue

        else:
            orth_gene_list = [gene]


        if SC_genes_repeats[species][gene] < len(orth_gene_list):
        # if SC_genes_repeats[species][gene] < len(orth_dict[gene]):
            SC_genes_repeats[species][gene] = len(orth_gene_list)
            # SC_genes_repeats[species][gene] = len(orth_dict[gene])
        for orth_gene in orth_gene_list:
            # print orth_gene
            # if orth_gene != 'NONE':
                # SC_genes_final.append(gene)
            # print orth_gene
            species_data[species]['orth_genes'].append(orth_gene)

            for condition_index, condition in enumerate(condition_arrays):
                condition_data = species_data[species][condition]

                # print '--------' + condition  + '-------------'
                # print condition_arrays.loc[(slice(None), orth_gene),:][condition][0]
                if orth_gene != 'NONE':
                    species_data[species][condition_index].append(species_data[species][condition].loc[(slice(None), orth_gene)].values[0])
                    # print condition_data.loc[(slice(None), orth_gene)].values[0]
                else:
                    species_data[species][condition_index].append(float('nan'))
                    # species_data[species][condition_index].append(0.0)
                    # print 'nan'
                # print species_data[species][condition].loc[(slice(None), orth_dict[gene]),:]
        # print '-----------------------------------'




num_conditions = len(conditions)
SC_gene_plot_data = SC_genes['Log_Change'].values.tolist()

# data_index = 0
# for gene_index, gene in enumerate(SC_genes['ORF']):
#     # num_rows += SC_genes_repeats[gene]
#     data_index += 1
#     for repeat_index in range(SC_genes_repeats[gene]):
#         num_rows += 1
#         if repeat_index == 0:
#             SC_genes_labels.append(gene)
#         else:
#             SC_genes_labels.append(' ')
#             SC_gene_plot_data.insert(data_index, SC_gene_plot_data[data_index - 1])
    


num_columns = 1
num_rows = 0
SC_genes_labels = []
for species in species_data:
    # if len(species_data[species]['orth_genes']) > num_rows:
    #     num_rows = len(species_data[species]['orth_genes'])
    num_columns += species_data[species]['num_conditions']

    for gene_index, gene in enumerate(SC_genes['ORF']):
        # for repeat_index, repeat_number in range(SC_genes_repeats[gene]):
        if SC_genes_repeats[gene] < SC_genes_repeats[species][gene]:
            SC_genes_repeats[gene] = SC_genes_repeats[species][gene]

num_rows = 0
for gene_index, gene in enumerate(SC_genes['ORF']):
    num_rows += SC_genes_repeats[gene]
    SC_genes_labels.append(gene)
    for i in range(SC_genes_repeats[gene] - 1):
        SC_gene_plot_data.insert(gene_index + 1, SC_gene_plot_data[gene_index])
        SC_genes_labels.append(' ')




plotting_data = np.zeros((num_rows, num_columns))

plotting_data[:, 0] = SC_gene_plot_data
# print num_rows
# print num_columns

# for i in range(len(SC_gene_plot_data)):
#     print SC_genes_labels[i], SC_gene_plot_data[i]



# print plotting_data

# for species in species_list:


# print plotting_data

condition_labels = ['Susan']
species_label_indices = []
col_index = 1
col_offset = 0
for species_index, species in enumerate(species_list):



    # print '---------------------------------------------------'
    # print species
    # print '---------------------------------------------------'




    # condition_labels.append(species)
    # species_label_indices.append(col_index + species_index)
    species_label_indices.append(col_index)
    
    repeat_diff = {}
    
    for gene_index, gene in enumerate(SC_genes['ORF']):
        repeat_diff[gene] = SC_genes_repeats[gene] - SC_genes_repeats[species][gene]
        # if SC_genes_repeats[gene] > SC_genes_repeats[species][gene]:
            # print species_data[species][condition].loc[(slice(None), orth_gene)].values[0]


    # for condition_index, condition in enumerate(conditions):
    
    for condition_index in range(species_data[species]['num_conditions']):
        row_index = 0
        condition_plotting_data = species_data[species][condition_index]
        for gene_index, gene in enumerate(SC_genes['ORF']):
            row_index += 1
            for i in range(repeat_diff[gene]):
                # condition_plotting_data.append(float('nan'))



                condition_plotting_data.insert(row_index, float('nan'))
                row_index += 1

                # if species == 'Saccharomyces cerevisiae':
                #     print gene
                # if gene == 'YDL083C' and species == 'Saccharomyces cerevisiae':
                    # print 'HEREHERE'
                    # print condition_plotting_data

        # print '------------------------------------'
        # print species_data[species][condition].loc[15293]

        # condition_labels.append(condition)

        # print species_data[species]['condition_names'][condition_index]
        condition_labels.append(species_data[species]['condition_names'][condition_index])
        # condition_labels.append(conditions[col_index - 1])

        # print '------------------------------------------------'
        # print len(species_data[species][condition_index])
        # print species_data[species][condition_index]
        # print '------------------------------------------------'

        # plotting_data[:, species_index + condition_index + 1] = species_data[species][condition_index]

        # plotting_data[:, species_index + condition_index + 1] = condition_plotting_data

        # plotting_data[:, col_offset + condition_index + 1] = condition_plotting_data
        plotting_data[:, col_index] = condition_plotting_data


        col_index += 1



        # print species_data[species][condition][10001]['KLLA0A11440g']
        # for measurement in species_data[species][condition]:




        

# for row_index, gene in enumerate(SC_genes_final):
    # for condition_index, condition in enumerate(conditions):
        # for species_index, species in enumerate(species_list):
            # orth_dict = species_data[species]['orth']
            # print species_data[species][condition].loc[(slice(None), orth_dict[gene]),:]
    

# for condition_index, condition in enumerate(conditions):
    # for species_index, species in enumerate(species_list):
        # print species_data[species][condition][:]
        # plotting_data[]


cmap = mpl.cm.RdBu_r
cmap.set_bad('k',1.0)
# cmap.set_bad('k')
fig, ax = plt.subplots()

# condition_labels = []
ax = sns.heatmap(plotting_data, cmap = cmap, ax = ax, yticklabels = SC_genes_labels, xticklabels = condition_labels)
ylocs, ylabels = plt.yticks()
plt.setp(ylabels, rotation = 0)



xlocs, xlabels = plt.xticks()
plt.setp(xlabels, rotation = 90)


# xticks = range(num_columns + len(species_list))
# xticks = range(num_columns)

# ax.set_xticks(xticks, minor = True)

# va = np.zeros(num_columns + len(species_list))
# va = np.zeros(num_columns)

species_labels = []
for species in species_list:
    species_labels.append(species_name_dict[species])

newax = ax.twiny()
# newax.grid(b=False)

# species_label_indices = [0, 0.1,0.2,0.3,0.4]
for i in range(len(species_label_indices)):
    species_label_indices[i] = species_label_indices[i] / float(num_columns)


newax.set_xticklabels(species_labels)
newax.set_xticks(species_label_indices)

# newax.set_frame_on(True)
# newax.patch.set_visible(False)
# for sp in newax.spines.itervalues():
#     sp.set_visible(False)
# newax.spines["bottom"].set_visible(True)





# ax2 = fig.add_axes((0.1,0.1,0.8,0.0))
# ax2.yaxis.set_visible(False)
# # new_tick_locations = np.array([1,6,10,12])
# new_tick_locations = [xlocs[i] for i in species_label_indices]

# # def tick_function(X):

# print species_labels
# ax2.set_xticks(new_tick_locations)
# ax2.set_xticklabels(species_labels)


# newax.set_xticks(new_tick_locations)
# newax.set_xticklabels(tick_function(new_tick_locations))
# newax.set_xlabel(r"Modified x-axis: $1/(1+X)$")





# fig.subplots_adjust(bottom=0.8)
# newax.set_frame_on(True)
# newax.patch.set_visible(False)
# newax.xaxis.set_ticks_position('bottom')
# newax.xaxis.set_label_position('bottom')
# newax.spines['bottom'].set_position(('outward', 40))

# newax.set_xlabel('blah')

# newax.set_xticklabels(species_labels)



# for i in range(len(va)):
#     if i in species_label_indices:
#         va[i] -= 0.1



# for t, y in zip(ax.get_xticklabels(), va ):
#     t.set_y(y)

# print ax.get_xticklabels()


plt.tight_layout()
plt.show()




