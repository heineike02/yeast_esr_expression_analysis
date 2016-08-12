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
import pka_library
import scipy


# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #

POSITIVE_LOG_CUTOFF = 5.5
NEGATIVE_LOG_CUTOFF = -5.5


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




# --------------------------------------------------------------------------------- #

# --------------------------------------------------------------------------------- #



orth_dict = io_library.read_orth_lookup_table('Saccharomyces cerevisiae', 'Candida glabrata')
high_genes, low_genes = pka_library.get_genes_susan(POSITIVE_LOG_CUTOFF, NEGATIVE_LOG_CUTOFF)
SC_genes = pka_library.concat_SC_genes(high_genes, low_genes)
species_data, conditions, SC_genes_repeats = pka_library.compile_species_data(species_list, SC_genes)
gla_data = species_data['Candida glabrata']



high_genes, low_genes = pka_library.get_genes_susan(POSITIVE_LOG_CUTOFF, -1000.0)
SC_genes_high = pka_library.concat_SC_genes(high_genes, low_genes)
species_data_high, conditions_high, SC_genes_repeats_high = pka_library.compile_species_data(species_list, SC_genes_high)
gla_data_high = species_data_high['Candida glabrata']


high_genes, low_genes = pka_library.get_genes_susan(1000, NEGATIVE_LOG_CUTOFF)
SC_genes_low = pka_library.concat_SC_genes(high_genes, low_genes)
species_data_low, conditions_low, SC_genes_repeats_low = pka_library.compile_species_data(species_list, SC_genes_low)
gla_data_low = species_data_low['Candida glabrata']




missing_total_orthologs = 0
for SC_gene in orth_dict:
	CG_orth = orth_dict[SC_gene][0]
	if CG_orth == 'NONE':
		missing_total_orthologs += 1	


missing_top_orthologs = 0
for gene in gla_data['orth_genes']:
	if gene == 'NONE':
		missing_top_orthologs += 1

missing_high_orthologs = 0
for gene in gla_data_high['orth_genes']:
	if gene == 'NONE':
		missing_high_orthologs += 1


missing_low_orthologs = 0
for gene in gla_data_low['orth_genes']:
	if gene == 'NONE':
		missing_low_orthologs += 1



# total_rand_vect = np.zeros(len(orth_dict))
# total_rand_vect[0:missing_total_orthologs] = 1


# SAMPLES = 10000
# missing_high_values = []
# missing_low_values = []
# missing_both_values = []

# for i in range(SAMPLES):
# 	missing_high_counter = 0
# 	missing_low_counter = 0
# 	missing_both_counter = 0 

# 	for total_counter in range(len(SC_genes)):
# 		missing_both_counter += np.random.choice(total_rand_vect)
# 	missing_both_values.append(float(missing_both_counter) / len(SC_genes))

# 	for total_counter in range(len(SC_genes_high)):
# 		missing_high_counter += np.random.choice(total_rand_vect)
# 	missing_high_values.append(float(missing_high_counter) / len(SC_genes_high))

# 	for total_counter in range(len(SC_genes_low)):
# 		missing_low_counter += np.random.choice(total_rand_vect)
# 	missing_low_values.append(float(missing_low_counter) / len(SC_genes_low))

# expected_high_value = np.average(missing_high_values)
# expected_low_value = np.average(missing_low_values)
# expected_both_value = np.average(missing_both_values)

# std_high = np.std(missing_high_values)
# std_low = np.std(missing_low_values)
# std_both = np.std(missing_both_values)


# print expected_both_value
# print expected_low_value
# print expected_high_value


# total_observation = float(missing_total_orthologs) / len(orth_dict) 
# both_observation = float(missing_top_orthologs) / len(SC_genes)
# high_observation = float(missing_high_orthologs) / len(SC_genes_high)
# low_observation = float(missing_low_orthologs) / len(SC_genes_low)


# scipy.stats.binom_test(missing_top_orthologs, n=len(SC_genes), p=float(missing_top_orthologs) / len(SC_genes))

print '------------------------------------------------------------------------------------------------------------'
print '------------------------------------------------------------------------------------------------------------'
print 'Number of Missing C.Gla orthologs of total S.Cer pka knockout list: ' + str(missing_total_orthologs)
print 'Total Number of activated/repressed S.Cer genes: ' + str(len(orth_dict))
print 'Fraction of missing C.Gla orthologs in top activated/repressed S.Cer genes: ' + str(float(missing_total_orthologs) / len(orth_dict))
print 'p-value (two-sided binomial test): ' + str(scipy.stats.binom_test(missing_total_orthologs, n=len(orth_dict), p=0.173129251701))
print ' '
print ' '
print ' '
print 'Number of Missing C.Gla orthologs of most activated/repressed S.Cer genes: ' + str(missing_top_orthologs)
print 'Total Number of top activated/repressed S.Cer genes: ' + str(len(SC_genes))
print 'Fraction of missing C.Gla orthologs in top activated/repressed S.Cer genes: ' + str(float(missing_top_orthologs) / len(SC_genes))
# print 'Expected Fraction of missing C.Gla orthologs in top activated/repressed S.Cer genes: ' + str(expected_both_value)
print 'p-value (two-sided binomial test): ' + str(scipy.stats.binom_test(missing_top_orthologs, n=len(SC_genes), p=0.173129251701))
print ' '
print ' '
print ' '
print 'Number of Missing C.Gla orthologs of most activated S.Cer genes: ' + str(missing_high_orthologs)
print 'Total Number of top activated S.Cer genes: ' + str(len(SC_genes_high))
print 'Fraction of missing C.Gla orthologs in top activated S.Cer genes: ' + str(float(missing_high_orthologs) / len(SC_genes_high))
# print 'Expected Fraction of missing C.Gla orthologs in top activated S.Cer genes: ' + str(expected_high_value)
print 'p-value (two-sided binomial test): ' + str(scipy.stats.binom_test(missing_high_orthologs, n=len(SC_genes_high), p=0.173129251701))
print ' '
print ' '
print ' '
print 'Number of Missing C.Gla orthologs of most repressed S.Cer genes: ' + str(missing_low_orthologs)
print 'Total Number of top repressed S.Cer genes: ' + str(len(SC_genes_low))
print 'Fraction of missing C.Gla orthologs in top repressed S.Cer genes: ' + str(float(missing_low_orthologs) / len(SC_genes_low))
# print 'Expected Fraction of missing C.Gla orthologs in top repressed S.Cer genes: ' + str(expected_low_value)
print 'p-value (two-sided binomial test): ' + str(scipy.stats.binom_test(missing_low_orthologs, n=len(SC_genes_low), p=0.173129251701))
print '------------------------------------------------------------------------------------------------------------'
print '------------------------------------------------------------------------------------------------------------'




