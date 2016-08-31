# import regev_library
# import numpy as np
# import os


# cwd = os.getcwd()

# total_data = np.load(cwd + '/stored_data/total_data/total_data.npy')[()]

# susan_list = total_data['Saccharomyces cerevisiae']['susan']['Genes']

# for gene in susan_list:
# 	print gene


import sys
import os
import csv
import scipy.stats
import single_score_library
import params
import numpy as np


cwd = os.getcwd()
total_data = np.load(cwd + '/stored_data/total_data/total_data.npy')[()]

susan_genes = total_data['Saccharomyces cerevisiae']['susan']['Genes']
high_genes = []
low_genes = []

for index, gene in enumerate(susan_genes):
	if total_data['Saccharomyces cerevisiae']['susan']['Values'][index] > 0:
		high_genes.append(susan_genes[index])
	elif total_data['Saccharomyces cerevisiae']['susan']['Values'][index] < 0:
		low_genes.append(susan_genes[index])


klac_genes = []

fname = os.path.normpath(cwd+'/Klac_low_genes.txt')
counts = {}
genes = {}
csvfile = open(fname, 'rb')
reader = csv.reader(csvfile, delimiter = '\n')

for row in reader:
	klac_genes.append(row[0])

klac_high = []
klac_low = []
klac_none = []
for gene in klac_genes:
	if gene in high_genes:
		klac_high.append(gene)
	elif gene in low_genes:
		klac_low.append(gene)
	else:
		print gene
		klac_none.append(gene)

print len(klac_high)


