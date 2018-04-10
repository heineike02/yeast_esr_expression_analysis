# -*- coding: utf-8 -*-
#from IPython.core.debugger import Tracer
import os 
import pandas as pd
import numpy as np
import re
import math
import scipy.stats as stats
from collections import Counter
import subprocess
print('I am importing io_library')
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from ete3 import Tree
import requests
from lxml import etree    #parses xml output

#Indicate operating environment and import core modules
location_input = input("what computer are you on? a = Bens, b = gpucluster, c = other   ")
location_dict = {'a': "C:\\Users\\heine\\github\\expression_broad_data", 'b': "/home/heineike/github/expression_broad_data",'c':'you need to add your location to the location_dict'}
base_dir = location_dict[location_input]
print("base directory is " + base_dir)
data_processing_dir = base_dir + os.sep + os.path.normpath("expression_data") + os.sep
print("data processing dir is " + data_processing_dir )

def tryfloatconvert(value, default):
    try:
        return float(value)
    except ValueError:
        return default

def parse_raw_exp(species):
    #Parses raw expression data for a given species.  Ouputs a dataframe. 
    #org_dict = {'Kluyveromyces lactis': ('GSE22198_family.soft','SystematicName','KLac',KL_soft_alt_dict), 'Saccharomyces cerevisiae': ('GSE22204_family.soft','ACCESSION_STRING','SCer',SC_soft_alt_dict)}
    raw_exp_datasets = {'Kluyveromyces lactis': 'GSE22198_family.soft', 'Saccharomyces cerevisiae': 'GSE22204_family.soft', 'Candida glabrata':'GSE22194_family.soft', 'Naumovozyma castellii' : 'GSE22200_family.soft', 'Saccharomyces bayanus' : 'GSE22205_family.soft', 'Saccharomyces mikatae': 'GSE22201_family.soft', 'Lachancea waltii': 'GSE22199_family.soft', 'Saccharomyces paradoxus': 'GSE22193_family.soft', 'Lacancea kluyverii': 'GSE22202_family.soft', 'Debaryomyces hansenii': 'GSE22191_family.soft'}
    #raw_exp_fnames = {'Kluyveromyces lactis': 'KLac_raw.csv', 'Saccharomyces cerevisiae': 'SCer_raw.csv'}
    orf_header_names = {spec : 'SystematicName' for spec in ['Kluyveromyces lactis', 'Candida glabrata', 'Naumovozyma castellii', 'Saccharomyces bayanus', 'Saccharomyces mikatae', 'Lachancea waltii', 'Saccharomyces paradoxus','Debaryomyces hansenii']}
    orf_header_names['Saccharomyces cerevisiae']='ACCESSION_STRING'
    orf_header_names['Lacancea kluyverii'] = 'Gene'
    
    #data cleaning notes: 
    #Kluyveromyces Lactis: No expression value listed - made NA
    #gene, line
    #11892, 7334
    #14495, 9937
    #13495, 14311
    #14052, 14868
    #10964, 17154
    #12500, 18690
    #Did not need to clean data for other species - made an if loop to catch it. 
    
    #KLac had four samples    
    input_fname = os.path.normpath(data_processing_dir + 'regev_data/raw_exp/' + raw_exp_datasets[species])
    with open(input_fname) as f:
        # identify samples in the dataset
        # Skips text before the line that says beginning of the interesting block:
        sample_datasets = []
        for line in f:
            if line.split()[0]== '!Series_sample_id':
                break
        
        sample_datasets.append(line.split()[2])
        #Tracer()()
        for line in f: 
            if line.split()[0] != '!Series_sample_id':
                break
            sample_datasets.append(line.split()[2])
        
        # get list of orfs with index numbers
        
        # Find line that starts table listing gene names and index numbers
        for line in f: 
            if line.split()[0] == '!platform_table_begin':
                break
        
        #Find index of header that will identify orf
        line = next(f)
        linesp = line.split()
        orf_header_ind = linesp.index(orf_header_names[species])
        
        data_dict = {}
        platform_dict = {}
        for line in f: 
            linesp = line.split('\t')
            if linesp[0] == "!platform_table_end\n":
                #include new line because this splits on tab.  
                break 
            
            orf_ind = linesp[0]
            #S. Cerevisiae orf names are formatted differently than other species. 
            if species == 'Saccharomyces cerevisiae':                       
                orf_name = linesp[orf_header_ind].split('|')[1]
            else: 
                orf_name = linesp[orf_header_ind]
                            
            platform_dict[orf_ind] = orf_name  
        
        data_dict['orf_name'] = platform_dict
        
        #builds a dictionary with primary key and experimental value for each id number. 
        #adds these as part of the data_dict.  
        for jj in range(len(sample_datasets)):
            sample_data = {}
            #finds sample name: 
            for line in f: 
                linesp = line.split()
                if linesp[0] == '^SAMPLE':
                    break 
            sample_name = linesp[2]
            
            #builds dictionary of id numbers and experimental values.     
            for line in f: 
                linesp = line.split()
                if linesp[0] == '!sample_table_begin':
                    break 
            
            #skip header line
            next(f)
            for line in f: 
                linesp = line.split()
                if linesp[0] == '!sample_table_end':
                    break
                if len(linesp) == 1:
                    linesp.append('nan')
                    print('Species: {}, Line {} was blank, made nan'.format(species, linesp[0]))
                sample_data[linesp[0]] = tryfloatconvert(linesp[1],np.nan)
                
                
            data_dict[sample_name]=sample_data
        
        raw_exp = pd.DataFrame(data_dict) 
        #make orf index and ID a field
    
    #raw_exp[sample_datasets]= raw_exp[sample_datasets].applymap(lambda x: tryfloatconvert(x,np.nan))
    raw_exp['Mean'] = raw_exp[sample_datasets].mean(axis = 'columns')
    ids = raw_exp.index
    raw_exp['ID']=ids
    orfs = raw_exp['orf_name']
    raw_exp = raw_exp.set_index('orf_name')
    orf_lookup = pd.Series(orfs.values, index = ids.values ) 
    orf_lookup = orf_lookup.sort_index()
    return raw_exp, orf_lookup                    
    
    
def parse_micro_data(species, exptype, orf_lookup): 
    #load raw data for a given species/experiment type pair.  Third argument is an orf lookup table
    #which is a from parse_raw_exp
    exptype_prefixes = {'Growth': '36253', 'Stress': '38478'}
    platform_prefixes = {'Kluyveromyces lactis': '10499', 'Saccharomyces cerevisiae': '9294', 'Candida glabrata':'10497', 'Naumovozyma castellii' : '10501', 'Saccharomyces bayanus' : '10505', 'Saccharomyces mikatae': '10502', 'Lachancea waltii': '10500', 'Saccharomyces paradoxus': '10496', 'Lacancea kluyverii': '10503', 'Debaryomyces hansenii': '15298', 'Vanderwaltozyma polyspora': '15297' }
    
    input_fname = os.path.normpath(data_processing_dir + '/regev_data/GSE' + exptype_prefixes[exptype] + '_' + exptype + '/GSE'+ exptype_prefixes[exptype] + '-GPL' + platform_prefixes[species] + '_series_matrix.txt')
    with open(input_fname) as f:
        if exptype == 'Growth': 
            
            #Make condition list and replicate list
            for line in f:
                if line == '\n':
                    #skips line if it is just a new line
                    pass
                elif line.split()[0]== '!Sample_title':
                    condition_titles = line
                    conditions = [re.split("[\-\[\]]",title)[1] for title in condition_titles.split('\t')[1:]]
                    replicates = [re.split("[\-\[\]]",title)[2] for title in condition_titles.split('\t')[1:]]
                    break
        
        elif exptype == 'Stress': 
            
            #Make replicate list 
            for line in f:
                if line == '\n':
                    #skips line if it is just a new line
                    pass
                elif line.split()[0]== '!Sample_source_name_ch1':
                    replicate_titles = line
                    replicates = [re.split("[\"_]",title)[2] for title in replicate_titles.split('\t')[1:]]
                    break
                    
            #Make condition list by concatenating stress and time information for each microarray
            for line in f:
                if line.split()[0]== '!Sample_characteristics_ch1':
                    condition_stress_titles = line
                    condition_stresses = [re.split("[\"\:]",title)[2].strip() for title in condition_stress_titles.split('\t')[1:]]
                    break
            for line in f:
                if line.split()[0]== '!Sample_characteristics_ch1':
                    condition_time_titles = line
                    condition_times = [re.split("[\"\:]",title)[2].strip() for title in condition_time_titles.split('\t')[1:]]
                    break
            conditions = ['{}_{:03d}'.format(tup[0],int(tup[1])) for tup in zip(condition_stresses,condition_times)]

            
                
        #scroll to beginning of data table and extract ids for each experiment
        for line in f:
            if line.split()[0]== '"ID_REF"':
                condition_ids = line.split()[1:]
                condition_ids = [condition_id.strip('"') for condition_id in condition_ids]
                break
        
        #store data for each experiment
        expdata_ungrouped = {}    
        for line in f: 
            linesp = line.split()
            if linesp[0] in ['"DarkCorner"','!series_matrix_table_end' ]:
                #"DarkCorner" is the beginning in a set of controls for the SCer chip.  
                # The data table for the other species ends at the line !series_matrix_table_end
                break
            chip_location_ID = linesp[0].strip('"')
            chip_location_values = linesp[1:]
            chip_location_values = [tryfloatconvert(value,np.nan) for value in chip_location_values]
            expdata_ungrouped[chip_location_ID] = chip_location_values
        
        #build dataframe from expdata_ungrouped
        expdata = pd.DataFrame.from_dict(expdata_ungrouped, orient = 'index')
        #sort ids alphanumerically (not sure whether this is necessary)
        expdata.sort_index(axis=0, inplace=True)
        #Sets a multi index for conditions and replicates
        col_mindex_arrays = [conditions,replicates,condition_ids]
        col_mindex_tuples = list(zip(*col_mindex_arrays))
        col_mindex = pd.MultiIndex.from_tuples(col_mindex_tuples, names=['conditions', 'replicates','array_ids'])
        expdata.columns = col_mindex
        #sort by conditions
        expdata_sorted = expdata.sort_index(level = 'conditions',axis = 'columns')
        #add in index of gene names (maybe have multi-index or just replace id numbers)
        if False in (expdata_sorted.index == orf_lookup.index):
            print("Error: ID mismatch between experiment data and orf lookup table. Species = {}, Experiment Type = {}".format(species, exptype))
            return
        
        #if no error, continue here
        print("All ID's match between experiment data and orf lookup table. Species = {}, Experiment Type = {}".format(species, exptype))
        
        mindex_arrays = [np.array(orf_lookup.index),np.array(orf_lookup)]
        mindex_tuples = list(zip(*mindex_arrays))
        expdata_sorted.index = pd.MultiIndex.from_tuples(mindex_tuples, names=['ID','orf_name'])
        

          
        

    
    return expdata_sorted
    

def make_data_tables(species_list,fname_out_bases, base_dir):
    
    for jj,spec in enumerate(species_list):  
    
        if spec != 'Vanderwaltozyma polyspora':     #Vpol doesn't have raw expression data
            #Generate raw expression data
            raw_exp, orf_lookup = parse_raw_exp(spec)

            #save raw expression data to a csv file
            fname = os.path.normpath(data_processing_dir + "regev_data/raw_exp/"  + fname_out_bases[jj] + '_raw_exp.csv')
            raw_exp.to_csv(fname)
            print(fname + ' saved')
        elif spec == 'Vanderwaltozyma polyspora':    #still need an orf lookup.  perhaps should do them all with these tables. 
            array_table = pd.read_table(data_processing_dir+os.path.normpath('regev_data/Vpol_array_table.txt'), dtype = 'str')
            orf_lookup = array_table.loc[:,['ID','ORF']].set_index('ID').squeeze()
        
        #Generate data for microarrays
        growth_exp = parse_micro_data(spec,'Growth',orf_lookup)
        fname = os.path.normpath(data_processing_dir + "regev_data/GSE36253_Growth/"  + fname_out_bases[jj] + '_growth.csv' )
        growth_exp.to_csv(fname)
        print(fname + ' saved')
        
        if not(spec in {'Saccharomyces bayanus','Saccharomyces mikatae','Saccharomyces paradoxus', 'Lacancea kluyverii','Debaryomyces hansenii','Vanderwaltozyma polyspora'}):
            #There is no stress dataset for these species
            stress_exp = parse_micro_data(species_list[jj],'Stress',orf_lookup)
            fname = os.path.normpath(data_processing_dir + "regev_data/GSE38478_Stress/"  + fname_out_bases[jj] + '_stress.csv' )
            stress_exp.to_csv(fname)
            print(fname + ' saved')
    
    return 

def combine_growth_stress_datasets(species):
    #species can be SCer, CGla, SCas, KLac.  No stress dataset for SBay 
    fname = os.path.normpath(data_processing_dir + "regev_data/GSE36253_Growth/"  + species + '_growth.csv' )
    growth_exp = pd.read_csv(fname,header = [0,1,2], index_col = [0,1])
    print(fname + ' growth microarray dataset loaded')

    #group by conditions and take mean
    growth_replicate_groups = growth_exp.groupby(axis = 1, level = 'conditions')
    growth_exp_avg = growth_replicate_groups.aggregate(np.mean)

    fname = os.path.normpath(data_processing_dir + "regev_data/GSE38478_Stress/"  + species + '_stress.csv' )
    stress_exp = pd.read_csv(fname,header = [0,1,2], index_col = [0,1])
    print(fname + ' stress microarray dataset loaded')

    #group by condition and take mean
    stress_replicate_groups = stress_exp.groupby(axis = 1, level = 'conditions')
    stress_exp_avg = stress_replicate_groups.aggregate(np.mean)

    #combine growth and stress average expression datasets. 
    if False in stress_exp_avg.index==growth_exp_avg.index:
        print("Error: ID mismatch between condition data. Species = {}".format(species))
    growth_stress_data = pd.concat([growth_exp_avg,stress_exp_avg], axis = 1)

    #gets rid of ID index
    growth_stress_data.reset_index(level=0, inplace=True)
    fname_out = os.path.normpath(data_processing_dir + 'regev_data/' + species + '_growth_stress.csv')  
    growth_stress_data.to_csv(fname_out)
    print('combined dataset saved as ' + fname_out )

    return growth_stress_data
    
def read_SGD_features():
    
    #Read in orf/name file and make it a dictionary
    # Gabe 7/12/16
    # SC_features_fname = os.path.normpath(data_processing_dir + "\ortholog_files\\SGD_features.tab")
    SC_features_fname = os.path.normpath(data_processing_dir + "/ortholog_files/SGD_features.tab")

    SC_features = pd.read_csv(SC_features_fname, sep = '\t', header=None)
    SC_orfs = SC_features.groupby(1).get_group('ORF')
    
    #Makes a dictionary to look up orfs by gene names.  This won't include all orfs - those without names had NaN in column 4 so 
    #are presumably left out. 
    SC_orfs_lookup = dict(zip(SC_orfs[4], SC_orfs[3]))
    SC_genename_lookup = dict(zip(SC_orfs[3], SC_orfs[4]))
    SC_features_lookup = dict(zip(SC_orfs[3], SC_orfs[15]))
       
    return SC_orfs_lookup, SC_genename_lookup, SC_features_lookup

def get_sgd_description(sc_genename_list):
    SGD_features = read_SGD_features()
    description_dict = SGD_features[2]

    description_list = []
    for gene in sc_genename_list:
        try: 
            description = description_dict[gene]
        except KeyError: 
            print('description lookup failed for ' + gene)
            description = "Lookup Failed"
        description_list.append(description)
    
    return description_list


def read_orth_lookup_table(species1, species2, orth_dir):
    #ensure orth_dir has a separator at the end of it.
    #For a given species read in the ortholog file, make a dictionary
    #formerly used the full name of the species, but 2/20/2018 switched to using four letter abbreviations. 
    #orth_file_abbrev = {'Kluyveromyces lactis': 'Klac', 'Saccharomyces cerevisiae': 'Scer', 'Candida glabrata':'Cgla', 'Saccharomyces castellii' : 'Scas', 'Saccharomyces bayanus' : 'Sbay'}
    orth_fname = species1 + "-" + species2 + "-orthologs.txt"
    orth_fname = os.path.normpath(orth_dir + orth_fname)
    
    #There are some orfs that have multiple orthologs - in that case both will be used
    with open(orth_fname) as f:
        orth_lookup = {}
        for line in f:
            linesp = line.split()
            orth_lookup[linesp[0]]= linesp[1:]

    return orth_lookup  
    
def get_gasch_ESR_list(act_rep):
    #For a file from the Gasch 2000 supplement, read in the data
    fname = os.path.normpath(data_processing_dir + "/gasch_data/gasch_fig3_" + act_rep + "_ESR.txt")
    
    with open(fname) as f:
        out = []
        #Skips first two lines
        next(f)
        next(f)
        for line in f:
            linesp = line.split()
            out.append(linesp[0])

    return out

def all_gasch_conditions(fname):
    #returns a list of all the conditions in the Gacsch data set
    #
    #two main options for fname are 
    #  gasch_complete_dataset.txt
    #  gasch_fig1_all_conditions.txt
        
    fname_full = os.path.normpath(data_processing_dir + "/gasch_data/" + fname)
    with open(fname_full) as f:
        conditions = next(f).split("\t")
    #remove newline from final condition
    conditions[-1] = conditions[-1].strip('\n')
    #remove items that are not conditions
    conditions = [condition for condition in conditions if not(condition in {'GID','UID','NAME','GWEIGHT',''})]
    return conditions

def read_gasch_data(conditions,fname):
    #For selected conditions from the Gasch 2000 supplement, extract the data as 
    #a dataframe with a row for each gene and each column a condition
    #possible conditions are given by the function all_gasch_conditions()
    #
    #You must choose a data file to extract the data from.  
    #The data in gasch_complete_dataset.txt is not "zero-transformed" per the paper to reflect expression v.s. no condition
    #usually it is expression v.s. a reference pool of all similar conditions. 
    #
    #The data in i.e. gasch_fig1_all_conditions.txt is zero transformed (per the figures) except for the
    #carbon sources which are v.s. a reference pool. 
    
    if fname == "gasch_complete_dataset.txt":
        gene_ind = 0
    else:      # This is for when we use the data from the figure.  It has an extra E-weight row - not sure what that is for.  
        gene_ind = 1
    fname_full = os.path.normpath(data_processing_dir + "/gasch_data/" + fname)
    
    with open(fname_full) as f:
        header = next(f).split("\t")
        #remove newline from final condition in header
        header[-1] = header[-1].strip('\n')
        condition_inds = [header.index(condition) for condition in conditions]
        #Skips first two lines
        next(f)
        exp_data = []
        gene_names = []
        for line in f:
            exp_values = line.split("\t")
            exp_data.append([tryfloatconvert(exp_values[condition_ind],None) for condition_ind in condition_inds])
            gene_names.append(exp_values[gene_ind])
        
    
    exp_df = pd.DataFrame(data=exp_data, index=gene_names, columns=conditions) 
    
    return exp_df

def parse_data_series_matrix_SC(desired_conditions, data_dir, GEO_accession):
    #Extracts data for desired conditions from a series matrix file for S.Cerevisiae
    #Extract dictionary for each gene from family.soft file
    #Desired condition is a list of tuples, 
    #   the first entry is the name you want to have for your columns (need not match file)
    #   the second entry is the array designator (should match file)
    
    #Extract dictionary for each gene from family.soft file
    soft_fname = os.path.normpath(data_dir + GEO_accession + '_family.soft' )
    with open(soft_fname) as f:
        for line in f: 
            if line.split()[0] == '!platform_table_begin':
                break

        #Find index of header that will identify orf
        line = f.next()
        linesp = line.split()
        orf_header_ind = linesp.index('ORF')
        #skip first three lines to get to data
        platform_dict = {}
        for line in f: 
            linesp = line.split('\t')
            if linesp[0] == "!platform_table_end\n":
                #include new line because this splits on tab.  
                break 

            if linesp[0][0:6] == 'A_06_P':
                orf_ind = linesp[0]
                orf_name = linesp[orf_header_ind]
                platform_dict[orf_ind] = orf_name  

    
    
    #Extracts data for desired conditions from a series matrix file for S.Cerevisiae
    series_fname = os.path.normpath(data_dir + GEO_accession + '_series_matrix.txt')
    #GSM1423542: nmpp1 treatment (20 min) - ACY142 +nmpp1 / ACY142

    with open(series_fname) as f:
        #Find line that starts table listing gene names and index numbers
        for line in f: 
            if line.split('\t')[0] == '!series_matrix_table_begin\n':
                break

        exp_line = f.next()
        exp_list = [item.strip('"') for item in exp_line.split()]

        #extract data for only the listed conditions

        desired_arrays = [condition[1] for condition in desired_conditions]
        array_inds = [exp_list.index(array) for array in desired_arrays]

        orf_list = []
        exp_value_list = []

        for line in f:
            linesp = line.split('\t') 
            if linesp[0] == '!series_matrix_table_end\n':
                break       

            if linesp[0].strip('"')[0:6] == 'A_06_P':
                orf_ind = linesp[0].strip('"')
                orf = platform_dict[orf_ind]
                orf_list.append(orf)
                exp_values = [float(linesp[array_ind]) for array_ind in array_inds]
                exp_value_list.append(exp_values)

    #Make dataframe with orf values as index. 
    data = pd.DataFrame(exp_value_list, index = orf_list, columns = [condition[0] for condition in desired_conditions])

    return data

def load_oshea_NMPP1_data(): 
    #Import SCer NMPP1 microarray data from O'shea microarrays
    title_list = '"treatment: 0 min: No 1-NM-PP1"	"treatment: 10 min: 3 µM 1-NM-PP1"	"treatment: 20 min: 3 µM 1-NM-PP1"	"treatment: 30 min: 3 µM 1-NM-PP1"	"treatment: 40 min: 3 µM 1-NM-PP1"	"treatment: 50 min: No 1-NM-PP1"	"treatment: 60 min: No 1-NM-PP1"	"treatment: 70 min: No 1-NM-PP1"	"treatment: 0 min: No 1-NM-PP1"	"treatment: 10 min: 3 µM 1-NM-PP1"	"treatment: 20 min: 3 µM 1-NM-PP1"	"treatment: 30 min: No 1-NM-PP1"	"treatment: 40 min: No 1-NM-PP1"	"treatment: 50 min: No 1-NM-PP1"	"treatment: 60 min: No 1-NM-PP1"	"treatment: 70 min: No 1-NM-PP1"	"treatment: 0 min: No 1-NM-PP1"	"treatment: 10 min: 120 nM 1-NM-PP1"	"treatment: 20 min: 120 nM 1-NM-PP1"	"treatment: 30 min: No 1-NM-PP1"	"treatment: 40 min: No 1-NM-PP1"	"treatment: 50 min: No 1-NM-PP1"	"treatment: 60 min: No 1-NM-PP1"	"treatment: 70 min: No 1-NM-PP1"	"treatment: 0 min: No 1-NM-PP1"	"treatment: 5 min: 750 nM 1-NM-PP1"	"treatment: 10 min: No 1-NM-PP1"	"treatment: 15 min: 750 nM 1-NM-PP1"	"treatment: 20 min: No 1-NM-PP1"	"treatment: 25 min: 750 nM 1-NM-PP1"	"treatment: 30 min: No 1-NM-PP1"	"treatment: 35 min_chip_1: 750 nM 1-NM-PP1"	"treatment: 35 min_chip 2: 750 nM 1-NM-PP1"	"treatment: 40 min: No 1-NM-PP1"	"treatment: 45 min: 750 nM 1-NM-PP1"	"treatment: 50 min: No 1-NM-PP1"	"treatment: 55 min: 750 nM 1-NM-PP1"	"treatment: 60 min: No 1-NM-PP1"	"treatment: 65 min: No 1-NM-PP1"	"treatment: 70 min: No 1-NM-PP1"	"treatment: 0 min: No 1-NM-PP1"	"treatment: 5 min: 750 nM 1-NM-PP1"	"treatment: 10 min: No 1-NM-PP1"	"treatment: 15 min: No 1-NM-PP1"	"treatment: 20 min: No 1-NM-PP1"	"treatment: 25 min: 750 nM 1-NM-PP1"	"treatment: 30 min_chip_1: No 1-NM-PP1"	"treatment: 35 min: No 1-NM-PP1"	"treatment: 30 min_chip_2: No 1-NM-PP1"	"treatment: 40 min: No 1-NM-PP1"	"treatment: 45 min: 750 nM 1-NM-PP1"	"treatment: 50 min: No 1-NM-PP1"	"treatment: 55 min: No 1-NM-PP1"	"treatment: 60 min: No 1-NM-PP1"	"treatment: 65 min: No 1-NM-PP1"	"treatment: 70 min: No 1-NM-PP1"'
    title_list = title_list.split('\t')
    title_list = [item.strip("\"") for item in title_list]
    id_ref_list = '"GSM812516"	"GSM812517"	"GSM812518"	"GSM812519"	"GSM812520"	"GSM812521"	"GSM812522"	"GSM812523"	"GSM812524"	"GSM812525"	"GSM812526"	"GSM812527"	"GSM812528"	"GSM812529"	"GSM812530"	"GSM812531"	"GSM812532"	"GSM812533"	"GSM812534"	"GSM812535"	"GSM812536"	"GSM812537"	"GSM812538"	"GSM812539"	"GSM812540"	"GSM812541"	"GSM812542"	"GSM812543"	"GSM812544"	"GSM812545"	"GSM812546"	"GSM812547"	"GSM812548"	"GSM812549"	"GSM812550"	"GSM812551"	"GSM812552"	"GSM812553"	"GSM812554"	"GSM812555"	"GSM812556"	"GSM812557"	"GSM812558"	"GSM812559"	"GSM812560"	"GSM812561"	"GSM812562"	"GSM812563"	"GSM812564"	"GSM812565"	"GSM812566"	"GSM812567"	"GSM812568"	"GSM812569"	"GSM812570"	"GSM812571"'
    id_ref_list = id_ref_list.split('\t')
    id_ref_list = [item.strip("\"") for item in id_ref_list]
    
    ind_minus = 0
    #[0,8,16,24,40] are indices for all 0 min no NMPP1 but those other sets are hybridized to a different combination of conditions. 
    id_ref_minus = [id_ref_list[ind_minus]]
    
    #add in titles for for 30,40 min +3uM NMPP1 microarrays
    titles_plus = ['treatment: 30 min: 3 µM 1-NM-PP1','treatment: 40 min: 3 µM 1-NM-PP1']
    inds_plus = [title_list.index(title)  for title in titles_plus]
    ids_ref_plus = [id_ref_list[jj] for jj in inds_plus]
    
    oshea_microarray_name_minus = ["no 1-NMPP1, 0min"]
    oshea_microarray_names_plus = ["3uM 1-NMPP1, 30min", "3uM 1-NMPP1, 40min"]
    oshea_microarray_names = oshea_microarray_name_minus+oshea_microarray_names_plus
    
    desired_conditions = zip(oshea_microarray_names, id_ref_minus+ids_ref_plus)
    
    oshea_exp_data_dir = data_processing_dir + '\GSE32703_NMPP1_SC\\' 
    GEO_accession = 'GSE32703'
    oshea_SC_PKA_data = parse_data_series_matrix_SC(desired_conditions, oshea_exp_data_dir, GEO_accession)
    
    #average controls and + NMPP1 conditions. 
    #using the other controls is silly because they are normalized to a pool of a different experiment. 
    #oshea_SC_PKA_data['minus_NMPP1'] = oshea_SC_PKA_data[oshea_microarray_names_minus].mean(axis = 1)
    oshea_SC_PKA_data['plus_NMPP1'] = oshea_SC_PKA_data[oshea_microarray_names_plus].mean(axis = 1)
    oshea_SC_PKA_data['SC_PKA(AS)+1NMPP1'] = oshea_SC_PKA_data.loc[:,'plus_NMPP1'].sub(oshea_SC_PKA_data.loc[:,'no 1-NMPP1, 0min'])
    
    #Add in SC common name
    SC_common_names = SC_common_name_lookup(oshea_SC_PKA_data.index)
        
    oshea_SC_PKA_data['SC_common_name'] = SC_common_names
    
    oshea_SC_PKA_data_summary = oshea_SC_PKA_data[['SC_PKA(AS)+1NMPP1', 'SC_common_name']]
    
    return oshea_SC_PKA_data_summary

def load_solis_NMPP1_data(column_to_use): 
    #Load Solis 2016 PKA inhibition data
    fname_solis_SC_PKA_data = os.path.normpath(data_processing_dir + '\SCer_NMPP1_RNA_Seq\solis_2016.xlsx')
    solis_SC_all_data = pd.read_excel(fname_solis_SC_PKA_data, header = 3)
    
    #The data set has a lot of duplicated indices - remove duplicates. 
    #Make a dictionary by going down the list - this will leave in a bunch of NaNs for duplicates. 
    solis_inhib_unique = {}
    columns = ['Inhib','X466.1NM']
    column_dict = {'Inhib':'HS-AA','X466.1NM':'WT'}
    for item in solis_SC_all_data.iterrows():
        gene = item[1]['locus']
        exp_val = [item[1][column] for column in columns]
        solis_inhib_unique[gene] = exp_val
    
    repeat_counter = Counter(solis_SC_all_data['locus'])
    for gene in repeat_counter.keys():
        if repeat_counter[gene] > 1:
            dupe_rows = solis_SC_all_data[solis_SC_all_data['locus']==gene]
            solis_inhib_unique[gene] = [np.nanmean(dupe_rows[column]) for column in columns]
    
    solis_SC = pd.DataFrame.from_dict(solis_inhib_unique, orient = 'index')
    solis_SC.columns = [column_dict[column] for column in columns]
    
    #Add an Sc Common Name
    SC_common_names = SC_common_name_lookup(solis_SC.index)
        
    solis_SC['SC_common_name'] = SC_common_names
    
    solis_PKA_data_summary = solis_SC[[column_to_use, 'SC_common_name']]
    solis_PKA_data_summary.columns = ['SC_PKA(AS)+1NMPP1', 'SC_common_name']
    
    return solis_PKA_data_summary   


def write_YGOB_orth_lookup_table(species1, species2, base_dir, all_ortholog_file):
    #for each position in species 1
    #Assign orthologs or 'NONE' to each column from Position 2 and 3
    fname = os.path.normpath(base_dir + all_ortholog_file)
    orth_positions = {'Kluyveromyces lactis': [15], 'Saccharomyces cerevisiae' : [11,21]}
    #YGOB_Pillars.txt order of species: 
    #    0    V. polyspora Position 1
    #    1    T. phaffii Position 1
    #    2    T. blattae Position 1
    #    3    N. dairenensis Position 1
    #    4    N. castellii Position 1
    #    5    K. naganishii Position 1
    #    6    K. africana Position 1
    #    7    C. glabrata Position 1
    #    8    S. bayanus var. uvarum Position 1
    #    9    S. kudriavzevii Position 1
    #    10   S. mikatae Position 1
    #    11   S. cerevisiae Position 1
    #    12   Ancestral Gene Order
    #    13   Z. rouxii
    #    14   T. delbrueckii
    #    15   K. lactis
    #    16   E. gossypii 
    #    17   E. cymbalariae
    #    18   L. kluyveri
    #    19   L. thermotolerans
    #    20   L. waltii
    #    21   S. cerevisiae Position 2
    #    22   S. mikatae Position 2
    #    23   S. kudriavzevii Position 2
    #    24   S. bayanus var. uvarum Position 2
    #    25   C. glabrata Position 2
    #    26   K. africana Position 2
    #    27   K. naganishii Position 2
    #    28   N. castellii Position 2
    #    29   N. dairenensis Position 2
    #    30   T. blattae Position 2
    #    31   T. phaffii Position 2
    #    32   V. polyspora Position 2
    
    species1_columns = orth_positions[species1]
    species2_columns = orth_positions[species2]
    
    with open(fname) as f:
        orth_lookup = []
        for line in f:
            linesp = line.split()
            for column1 in species1_columns: 
                if linesp[column1]!= '---':
                    species1_gene = [linesp[column1]]
                    species2_genes = []
                    for column2 in species2_columns: 
                        if linesp[column2]!='---': 
                            species2_genes.append(linesp[column2])
                    if len(species2_genes) == 0:
                        species2_genes.append('NONE')
                    orth_lookup.append(species1_gene + species2_genes)
                    
    orth_file_abbrev = {'Kluyveromyces lactis': 'Klac', 'Saccharomyces cerevisiae': 'Scer', 'Candida glabrata':'Cgla', 'Saccharomyces castellii' : 'Scas', 'Saccharomyces bayanus' : 'Sbay'}
    orth_lookup_outputfname = os.path.normpath(base_dir + '\expression_data\ortholog_files_YGOB\\' + orth_file_abbrev[species1] + "-" + orth_file_abbrev[species2] + "-orthologs.txt"  )
    orth_lookup_outputfile = open(orth_lookup_outputfname, 'w')
    
    for gene in orth_lookup:
        line = '\t'.join(gene)+'\n'
        orth_lookup_outputfile.write(line)
    orth_lookup_outputfile.close()
    
    return orth_lookup 

def load_YGOB_annotations(species, base_dir, species_tab_file):
    fname = os.path.normpath(base_dir + species_tab_file)
    
    with open(fname) as f:
        annotation_lookup = {}
        for line in f:
            linesp = line.split('\t')
            gene = linesp[0]
            annotation = linesp[8]
            annotation_lookup[gene] = annotation
    
    return annotation_lookup

def SC_common_name_lookup(gene_list):
    #SC Common Name lookup
    #Input is a list of orfs, output is a list of common names
    
    
    SC_orfs_lookup, SC_common_name_lookup, SC_features_lookup = read_SGD_features()
    SC_common_names = []
    
    
    for gene in gene_list: 
        try:  
            SC_common_name = SC_common_name_lookup[gene]
            if isinstance(SC_common_name,float):
                if math.isnan(SC_common_name):
                    SC_common_names.append(gene)
                else: 
                    print('float but not nan uh oh!')
            else: 
                SC_common_names.append(SC_common_name_lookup[gene])
        except KeyError: 
            SC_common_names.append(gene)
    
    return SC_common_names

def SC_orf_lookup_by_name(name_list):
    #Input is a list of common names, output is a list of orfs
    #would be nice if this worked with a single string as well as a list 
    #e.g. 'ERV14' as well as ['ERV14']
    
    SC_orfs_lookup, SC_common_name_lookup, SC_features_lookup = read_SGD_features()
    sc_genenames = []
    
    for name in name_list: 
        try:  
            sc_genename = SC_orfs_lookup[name]
        except KeyError: 
            print("S.Cer orf for " + name + "not found")
            sc_genename = name
        
        sc_genenames.append(sc_genename)
    
    return sc_genenames

def get_species_from_gene(gene, species_info):
    #extracts species for a given gene name: 
    #loop through species that have a prefix
    
    gene_spec = ''    
    for row in species_info[~species_info['ygob_gene_prefix'].isna()].loc[:,['ygob_gene_prefix','abbreviation']].itertuples():
        prefix = row.ygob_gene_prefix
        spec = row.abbreviation
        #if prefix is not a string, then it should be an np.nan
        if type(prefix) == str:
            found = re.findall(prefix,gene)
            if len(found)==1: 
                if gene_spec != '':
                    print('problem - gene had more than one prefix match')
                gene_spec = spec
            elif len(found)>1: 
                print('problem - gene name had more than one match to a prefix')
    if gene_spec == '':
        print(gene)
        #didn't match with a prefix - check Egos and Scer
        first_letter = gene[0]
        if first_letter == 'Y':
            gene_spec = 'Scer'
        elif first_letter == 'A':
            gene_spec = 'Egos'
        else:
            print('problem - no match with prefix or with Egos and Scer. ' + gene)
    
    return gene_spec

def extract_promoter_sequence(gene, prom_length):
    #extracts promoter sequence from NCBI for a given gene. 
    
    species_info = pd.read_csv(os.path.normpath(data_processing_dir + 'promoter_phylogenies/ygob_species_info.csv'))
       
    spec = get_species_from_gene(gene, species_info)
    
    #Check if the species is in NCBI
    ncbi_specs = species_info[species_info['in_ncbi']]['abbreviation']
    full_species_dict = dict(zip(species_info['abbreviation'],species_info['ncbi_name']))
    
    if spec in set(ncbi_specs):
        
        strand_inds = {"plus":"1","minus":"2"}    
        
        #converts gene name from ygob name to name in NCBI
        if spec == 'Vpol':
            gene = gene.replace('.','p')
        if spec == 'Klac':
            gene = gene[0:5] + '_' + gene[5:]
        if spec == 'Egos': 
            gene = 'AGOS_' + gene
            
        genus, spec_full = full_species_dict[spec].split(' ')
        gene_search_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=(' + gene + '%5BGene%20Name%5D)%20AND%20' + genus + '%20' + spec_full + '%5BOrganism%5D&retmode=json'
        gene_search_response = requests.get(gene_search_url)
        
        if gene_search_response.ok:
            #verify there is only one search result
            search_count = gene_search_response.json()['esearchresult']['count']
            if int(search_count) != 1: 
                raise ValueError('More or less than one search result for ' + gene + ' count = ' + search_count)
            entrez_gene_id = gene_search_response.json()['esearchresult']['idlist'][0]
        else:
            raise ValueError('Gene search response for ' + gene + ' not ok')
    
        # Extract ID fetch data
        gene_data_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=' + entrez_gene_id + '&retmode=xml'
        gene_data_response = requests.get(gene_data_url) 
        if gene_data_response.ok:
            gene_data_xml = etree.fromstring(gene_data_response.content)
            entrezgene_element = gene_data_xml.getchildren()[0]
            locus_element = entrezgene_element.find('Entrezgene_locus')
            #check that Entrezgene_locus is not None
            if locus_element == None: 
                raise ValueError('Entrezgene_locus not found ' + gene)
            locus_gene_commentary_list = locus_element.findall("Gene-commentary") 
            if len(locus_gene_commentary_list) == 1: 
                locus_gene_commentary = locus_gene_commentary_list[0]
            else: 
                print('NCBI has more than one gene commentary elements in the Entrezgene Locus looking for type 1 = genomic')
                #Keep only the gene commentary that has type = 1 
                locus_gene_commentary_shortlist = [locus_gene_commentary for locus_gene_commentary in locus_gene_commentary_list if locus_gene_commentary.find("Gene-commentary_type").text=='1']
                if len(locus_gene_commentary_shortlist) == 1: 
                    locus_gene_commentary = locus_gene_commentary_shortlist[0]
                else: 
                    print('Error - more than one gene commentary listed for genomic dna.  gene = ' + gene)
            if locus_gene_commentary.find("Gene-commentary_type").text != '1':
                print('Error - gene commentary type is not 1: genomic dna.  gene = ' + gene)
            seq_id = locus_gene_commentary.find("Gene-commentary_accession").text
            seq_ver = locus_gene_commentary.find("Gene-commentary_version").text
            #it seems that the "Seq-interval_strand" element only exists if it is plus so set default to plus
            seq_strand = "plus"
            seq_start = seq_stop = None
            for element in locus_gene_commentary.find("Gene-commentary_seqs").getchildren()[0].iter():
                if element.tag == "Seq-interval_from": 
                    seq_start = str(int(element.text)+1)  #query seems one off from the data in the gene commentary
                if element.tag == "Seq-interval_to":
                    seq_stop = str(int(element.text)+1)   #query seems one off from the data in the gene commentary
                if element.tag == "Seq-interval_strand":
                    seq_strand = element.getchildren()[0].get("value")
        else:
                raise ValueError('Gene Id fetch response for ' + gene + ' not ok.  ID = ' + entrez_gene_id)
    
        # extract chromosome sequence number, coordinates and strand.  Save into database. 
    
        # Add the operational promoter length to the correct side. 
        if seq_strand == "plus": 
            seq_start_new = str(int(seq_start)-prom_length)
            seq_stop_new = seq_stop
        elif seq_strand == "minus":
            seq_start_new = seq_start
            seq_stop_new = str(int(seq_stop)+prom_length)
    
        sequence_query = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=" + seq_id + "." + seq_ver + "&rettype=fasta&seq_start=" + seq_start_new + "&seq_stop=" + seq_stop_new + "&strand=" + strand_inds[seq_strand] 
        sequence_data_response = requests.get(sequence_query) 
        sequence_data = ''.join(sequence_data_response.text.split('\n')[1:-2])
        #should probably have check here to ensure that I don't parse the sequence output in a bad way - for now assuming it is second line through
        #the second to last line per the example. 
    
        # Add all protein sequences plus promoter sequen
        promoter_output = (seq_id + "." + seq_ver,[seq_start,seq_stop],seq_strand,sequence_data[0:prom_length],sequence_data[prom_length:])
    
    else: 
        print(gene + ' not in NCBI - no method to get promoter yet')
        promoter_output = np.nan
        
    return spec, promoter_output

def write_promoter_file(promoter_database, gene_list,fname):
    #Changed name from write_ame_promoter_file to write_promoter_file 10JAN18    
    with open(fname,'w') as f: 
        for gene in gene_list:
            try: 
                row = promoter_database.loc[gene,]
                header_line = '>' + row.name + ' ' + str(len(row.prom_seq)) + 'bp_upstream\n'
                seq_line = row['prom_seq'] + '\n'
                f.write(header_line)
                f.write(seq_line)
            except KeyError: 
                print(gene + " not in promoter data set.")
    return

def build_motif_dict(fname = data_processing_dir + os.path.normpath('motifs/JASPAR_CORE_2016_fungi.meme')): 
    motif_dict = {}
    with open(fname,'r') as f: 
        for line in f: 
            if line != '\n':
                if line.split()[0]=='MOTIF':
                    motif_dict[line.split()[1]]=line.split()[2]
    return motif_dict

def read_ame_output(fname):
    #reads in ame program output file from meme suits
    #an older version of Meme didn't list the name of the motif, so needed a motif_dict, but not necessary anymore. 
    #assumes first 12 lines are not data
    motif_id = []
    motif_name = []
    motif_consensus = []
    pval = []
    pval_corr = []
    with open(fname,'r') as f: 
        for jj in range(1,6):
            next(f)
        #extract the name of the motif file from the command
        command = next(f)
        motif_file = command.split()[-1].split('/')[-1]
        for jj in range(7,13):
            next(f)
        for line in f:
            if motif_file == 'JASPAR2018_CORE_fungi_redundant.meme':
                motif_id.append(line.split()[5])
                motif_name.append(line.split()[6])
                motif_consensus.append(line.split()[7].strip('()'))
            else: 
                motif_id.append(line.split()[5])
                motif_name.append(line.split()[5])
                motif_consensus.append(line.split()[6].strip('()'))
            #Add right-tailed p-value to list
            pval.append(float(line.split()[-9]))
            #Add right-tailed corrected p-value to list
            pval_corr.append(float(line.split()[-2]))


    ame_dict = {"motif_id": motif_id, "motif_name": motif_name, "motif_consensus": motif_consensus, "pval":pval, "pval_corr": pval_corr}
    ame_data = pd.DataFrame.from_dict(ame_dict)
    
    return ame_data

def run_ame_analysis(spec, target_gene_list, control_gene_list,target_fname_prefix, control_fname_prefix, motif, 
                     promoter_dir = {'KL': data_processing_dir+'/kl_promoters/' , 'SC': data_processing_dir + '/sc_promoters/'},
                     promoter_fname = {'KL': 'kl_promoters.pkl', 'SC': 'sc_promoters.pkl'},
                     ame_scoring = 'avg',
                     ame_method = 'ranksum',
                     ame_pvalue_threshold = '0.05' ):
    #runs ame program from meme software for given target gene lit, control gene list and set of motifs. 
    #extract promoters
    #spec: Species - either 'KL' or 'SC'
    #target gene list 
    #
    #ame_scoring.  Chose 'avg' as default - the average odds score for the sequence.  Since all the sequences are the same length, this seems like a good choice.  'total_hits' is also an option if you believe the number of hits is important, 
    #ame_method.  chose 'ranksum' as default - this is the default for the website.  'fisher' is also used, but this requires setting a threshold (default on the website is 1) for the score of each hit to be called a true positive.  
    
    promoters = pd.read_pickle(os.path.normpath(promoter_dir[spec] + promoter_fname[spec])) 
    target_promoters = promoters.loc[target_gene_list,]
    control_promoters = promoters.loc[control_gene_list,]  #for control using just promoters which have orthologs - but should I use all promoters? 


    fname_prefixes = {'target':target_fname_prefix, 'control': control_fname_prefix}
    promoter_lists = {'target':target_promoters, 'control': control_promoters}

    for gene_set in ['target','control']:
        fname = promoter_dir[spec] + 'promoter_sets/' + fname_prefixes[gene_set] + '_promoters.fasta'
        with open(fname,'w') as f:
            for row in promoter_lists[gene_set].itertuples():
                if isinstance(row.prom_seq,str):
                    header_line = '>' + row.Index + ' 700bp_upstream\n'
                    seq_line = row.prom_seq + '\n'
                    f.write(header_line)
                    f.write(seq_line)
                else:
                    print(row.Index + ' in ' + gene_set + ' set promoter value is not a string : ' + str(row.prom_seq)+ '. Skipping for promoter file.')


    #Use subprocess to run meme commands: 

    #ame --verbose 1 --oc . --control all_kl_promoters.fasta --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 0.05 mito_promoters_kl.fasta db/JASPAR/JASPAR_CORE_2016_fungi.meme
    motif_db = os.path.normpath(data_processing_dir + '/motifs/'+ motif['fname'])
    target_sequences = os.path.normpath(promoter_dir[spec] + 'promoter_sets/' + fname_prefixes['target'] + '_promoters.fasta')
    control_sequences = os.path.normpath(promoter_dir[spec] + 'promoter_sets/' + fname_prefixes['control'] + '_promoters.fasta')
    output_dir = os.path.normpath(promoter_dir[spec] + 'ame_output')
    file_prefix = target_fname_prefix + '_vs_' + control_fname_prefix + '_motif_' + motif['name'] + '_pVal_' + ame_pvalue_threshold


    ame_command = [ "/home/kieran/meme/bin/ame",
                  "--verbose", "2",
                  "--oc", output_dir,
                  "--control", control_sequences,
                  "--bgformat", "1", 
                  "--scoring", ame_scoring,
                  "--method", ame_method, 
                  "--pvalue-report-threshold", ame_pvalue_threshold, 
                  target_sequences,
                  motif_db]

    ame_output = subprocess.run(ame_command,stdout = subprocess.PIPE) 

    print("ame output return code = " + str(ame_output.returncode))

    #change file prefix
    for fname in ["ame.txt","ame.html"]:
        os.rename(output_dir + os.sep + fname, output_dir + os.sep + file_prefix + fname)
    
    return


#def run_meme_analysis(input_promoters_fname, background_promoters_fname):
#
#    pspgen_command = [ "/home/kieran/meme/bin/ame",
#                   "--verbose", "2",
#                   "--oc", output_dir,
#                   "--control", control_sequences,
#                   "--bgformat", "1", 
#                   "--scoring", ame_scoring,
#                   "--method", ame_method, 
#                   "--pvalue-report-threshold", ame_pvalue_threshold, 
#                   target_sequences,
#                   motif_db]

def merge_overlap_column(series_a, series_b):
    #for two string columns in an outer merge that should have identical entries except where 
    #there are np.nan values, outputs a list that has just one value and replaces np.nan with 
    #value that is present in either series a or series b. 
    #set consensus of x equal to None when not present instead of np.nan

    series_a = series_a.where(series_a.notnull(), other=None)
    series_b = series_b.where(series_b.notnull(), other=None)

    out_list = []
    for a,b in zip(series_a,series_b):
        if a==b: #np.isnan(b):
            out = a
        elif a:
            out = a
        elif b: 
            out = b
        else: 
            print("a doesn't match b a: " + a + " b: " + b)
        out_list.append(out)

    return out_list
    
def make_foldchange_subsets(kl_sc_PKA_data, pthreshold_KL, pthreshold_SC): 
    ##would be nice to have option to make thresholds on either padj or log2FoldChange
    #Highlight hits that are statistically significant for K.Lactis and S. Cerevisiae and breaks them down into activated and 
    #repressed for each group. 
    #output is a dictionary which links a gene subset name to a list of genes using the S.Cer orf name. 
    
    kl_sc_PKA_data_klsig = kl_sc_PKA_data[kl_sc_PKA_data['padj_KL'] < pthreshold_KL]
    kl_sc_PKA_data_scsig = kl_sc_PKA_data[kl_sc_PKA_data['padj_SC'] < pthreshold_SC]
    
    klscsig_ind = set(kl_sc_PKA_data_klsig.index)&set(kl_sc_PKA_data_scsig.index)
    klsig_scunsig_ind = set(kl_sc_PKA_data_klsig.index)-klscsig_ind
    scsig_klunsig_ind = set(kl_sc_PKA_data_scsig.index)-klscsig_ind
    unsig_ind = set(kl_sc_PKA_data.index)-(set(kl_sc_PKA_data_klsig.index)|set(kl_sc_PKA_data_scsig.index))
    
    print("At an adjusted pvalue threshold of {:.2E} for K.Lac and {:.2E} for S.Cer, there are"
           " {:d} genes significant for both species, {:d} genes significant for KL only,"
           " {:d} genes significant for SC only, and {:d} unsignificant genes".format(
                pthreshold_KL,
                pthreshold_SC,
                len(klscsig_ind), 
                len(klsig_scunsig_ind),
                len(scsig_klunsig_ind),
                len(unsig_ind)))
    
    kl_sc_PKA_data_klsig_scunsig = kl_sc_PKA_data.loc[list(klsig_scunsig_ind)]
    kl_sc_PKA_data_scsig_klunsig = kl_sc_PKA_data.loc[list(scsig_klunsig_ind)]
    kl_sc_PKA_data_klscsig = kl_sc_PKA_data.loc[list(klscsig_ind)]
    kl_sc_PKA_data_unsig = kl_sc_PKA_data.loc[list(unsig_ind)]
    
    #Significant for both, but up in SC and down in kl. 
    upsc_downkl = list(kl_sc_PKA_data_klscsig[(kl_sc_PKA_data_klscsig['log2FoldChange_SC']>0)&(kl_sc_PKA_data_klscsig['log2FoldChange_KL']<0)]['sc_genename'])
    print("Genes that are up in S.Cer and Down in K.Lac")
    print(upsc_downkl)
    #HEF3 {'YNL014W'}
    
    upkl_downsc = list(kl_sc_PKA_data_klscsig[(kl_sc_PKA_data_klscsig['log2FoldChange_SC']<0)&(kl_sc_PKA_data_klscsig['log2FoldChange_KL']>0)]['sc_genename'])
    print("Genes that are up in K.Lac and Down in S.Cer")
    print(upkl_downsc)
    
    #Add in the one gene that is significant for both but repressed in SC and activated in KL. 
    klsig_act_genes = list(kl_sc_PKA_data_klsig_scunsig[kl_sc_PKA_data_klsig_scunsig['log2FoldChange_KL']>0]['sc_genename'])
    for gene in upkl_downsc: 
        klsig_act_genes.append(gene)
    
    #Add in the one gene that is significant for both but activated in SC and repressed in KL.
    klsig_rep_genes = list(kl_sc_PKA_data_klsig_scunsig[kl_sc_PKA_data_klsig_scunsig['log2FoldChange_KL']<0]['sc_genename'])
    for gene in upsc_downkl: 
        klsig_rep_genes.append(gene)
    
    #Add in the one gene that is significant for both but activated in SC and repressed in KL.
    scsig_act_genes = list(kl_sc_PKA_data_scsig_klunsig[kl_sc_PKA_data_scsig_klunsig['log2FoldChange_SC']>0]['sc_genename'])
    for gene in upsc_downkl: 
        scsig_act_genes.append(gene)
    
    #Add in the one gene that is significant for both but repressed in SC and activated in KL. 
    scsig_rep_genes = list(kl_sc_PKA_data_scsig_klunsig[kl_sc_PKA_data_scsig_klunsig['log2FoldChange_SC']<0]['sc_genename'])
    for gene in upkl_downsc: 
        scsig_rep_genes.append(gene)
    
    gene_set_dict = {'klsig_act': klsig_act_genes,
                    'klsig_rep': klsig_rep_genes,
                    'scsig_act': scsig_act_genes,
                    'scsig_rep': scsig_rep_genes,
                    'klscsig_act': list(kl_sc_PKA_data_klscsig[(kl_sc_PKA_data_klscsig['log2FoldChange_SC']>0) & (kl_sc_PKA_data_klscsig['log2FoldChange_KL']>0)]['sc_genename']),
                    'klscsig_rep': list(kl_sc_PKA_data_klscsig[(kl_sc_PKA_data_klscsig['log2FoldChange_SC']<0) & (kl_sc_PKA_data_klscsig['log2FoldChange_KL']<0)]['sc_genename']),
                    'klscunsig': list(kl_sc_PKA_data_unsig['sc_genename'])
                    }
    return gene_set_dict

def kl_genename_convert(df): 
    #probably don't really need this one - just list one below. 
    #Input is dataframe with GFF version of kl genename as the index. 
    #Ouput is dataframe with standard version of kl genename as the index. 

    kl_genename = kl_genename_convert(df.index)
    df['kl_genename'] = kl_genename
    df.set_index('kl_genename',inplace = True)

    return df

def kl_genename_convert_list(kl_genes): 
    kl_genename = []
    for gene in kl_genes: 
        if gene[0:5]=='KLLA0':
            new_gene = gene.split('_')[0]+gene.split('_')[1]
        else: 
            new_gene = gene
        kl_genename.append(new_gene)

    return kl_genename
    
def load_goslim_data(GO_aspect):
    #The three GO_aspect values are: 
    #C = cellular_component
    #F = molecular_function
    #P = biological_process
    go_slims = pd.read_table(data_processing_dir + os.path.normpath('/go_terms/go_slim_mapping.tab'),header = None)
    go_slims.columns = ['sc_genename','sc_common_name','sgd_ID','GO_aspect','GO_term','GO_term_ID','feature_type']
    
    go_slims_aspect = go_slims[go_slims['GO_aspect']==GO_aspect]
    go_term_list = list(set(go_slims_aspect['GO_term']))
    
    return go_slims_aspect, go_term_list

def go_terms_for_genelist(gene_set_list, go_slims_aspect, go_term_list): 
    #for a given gene list provides a dataframe listing genes in that list for each go term. 
    go_term_data = []

    go_term_index = []

    for term in go_term_list: 
        term_genes = list(go_slims_aspect[go_slims_aspect['GO_term']==term]['sc_genename'])
        if len(term_genes)> len(set(term_genes)):
            print("Duplicate Term: " + term)
        subset_genes_in_goterm =  set(gene_set_list) & set(term_genes)
        N_subset_genes_in_goterm = len(subset_genes_in_goterm)
        N_genes_in_goterm = len(term_genes)
        if N_subset_genes_in_goterm >0:
            subset_genes_in_goterm_commonname = SC_common_name_lookup(subset_genes_in_goterm)
            go_term_data.append((N_subset_genes_in_goterm,
                                subset_genes_in_goterm,
                                subset_genes_in_goterm_commonname,
                                N_genes_in_goterm))
            go_term_index.append(term)

    go_term_df = pd.DataFrame(go_term_data, index = go_term_index,columns = ['N subset genes in goterm',
                                                                            'genes',
                                                                            'genes common name',
                                                                            'N genes in goterm'])
                                                                            
    return go_term_df

def go_term_enrichment(gene_set_list, background_genes, go_term_list, go_slims_aspect): 
    #input is gene_set_list, background_genes list, go_term_list, and go_slims_aspect dataframe with data for each go term.  
    #output is a dataframe with enrichment columns and full list of enriched genes. 
       
    missing_genes = ['Scer_YGOB_YDR134C', 'Scer_RDN18-1', 
                     'Scer_YGOB_SDC25', 'Scer_RDN25-1', 
                     'Scer_RDN58-1', 'Scer_YGOB_Anc_7.495', 
                     'Scer_YGOB_ADL119W', 'Scer_RDN5-1']
    
    
    bg_filtered_genes = set(background_genes) & set(missing_genes)
    if len(bg_filtered_genes)>0:
        print("The following missing genes are filtered out of the background set:")
        print(bg_filtered_genes)
    background_genes_filtered = list(set(background_genes)-set(missing_genes))
    
    
    #filter out missing genes 
    gene_set_filtered_genes = set(gene_set_list) & set(missing_genes)
    if len(gene_set_filtered_genes)>0:
        print("The following missing genes are filtered out of the input set:")
        print(gene_set_filtered_genes)
        gene_set_list = list(set(gene_set_list)-set(missing_genes))
    go_term_data = []
    go_term_list_filtered = go_term_list.copy()
    for term in go_term_list: 
        term_genes = list(go_slims_aspect[go_slims_aspect['GO_term']==term]['sc_genename'])
        if len(term_genes)> len(set(term_genes)):
            print("Duplicate Term for term " + term)
        subset_genes_in_goterm =  set(gene_set_list) & set(term_genes)
        N_subset_genes_in_goterm = len(subset_genes_in_goterm)
        N_subset_genes_notin_goterm = len(gene_set_list)-N_subset_genes_in_goterm
        N_bg_genes_in_goterm = len(set(background_genes_filtered) & set(term_genes))
        if N_bg_genes_in_goterm >0:
            N_bg_genes_notin_goterm = len(background_genes_filtered)-N_bg_genes_in_goterm
            oddsratio, pvalue = stats.fisher_exact([[N_subset_genes_in_goterm, N_bg_genes_in_goterm], [N_subset_genes_notin_goterm, N_bg_genes_notin_goterm]],alternative = 'greater')
            subset_genes_in_goterm_commonname = SC_common_name_lookup(subset_genes_in_goterm)
            go_term_data.append((N_subset_genes_in_goterm,
                                 len(gene_set_list),
                                 float(N_subset_genes_in_goterm)/float(len(gene_set_list)),
                                 float(N_bg_genes_in_goterm)/float(len(background_genes_filtered)),
                                 float(N_subset_genes_in_goterm)/float(N_bg_genes_in_goterm),
                                 pvalue,
                                 subset_genes_in_goterm,
                                 subset_genes_in_goterm_commonname,
                                 oddsratio))
        else:
            print(term + " term removed: no background genes")
            go_term_list_filtered.remove(term) 
    
    go_term_enrichment = pd.DataFrame(go_term_data, index = go_term_list_filtered,columns = ['N subset genes in goterm',
                                                                            'N genes in subset',
                                                                            'pct goterm in subset',
                                                                            'pct go term in background',
                                                                            'pct of go terms genes in subset',
                                                                            'pvalue',
                                                                            'genes',
                                                                            'genes common name',
                                                                            'oddsratio'])
    
    go_term_enrichment.sort_values(by = 'pvalue', inplace = True)
    
    return go_term_enrichment

def go_terms_by_gene(sc_genename_list): 
    GO_aspect = 'P'
    go_slims_aspect, go_term_list = load_goslim_data(GO_aspect)

    go_term_list = []
    for gene in sc_genename_list:
        go_term_list.append(", ".join(list(go_slims_aspect[go_slims_aspect['sc_genename'] == gene].loc[:,'GO_term'])))

    return go_term_list

def make_meme_promoter_files(gene_list, fname_prefix, spec_comparision_data): 

    #Read in the KL promoter database.  
    kl_promoters = pd.read_pickle(base_dir + os.sep + os.path.join("expression_data","kl_promoters","kl_promoters.pkl"))

    #Read in the SC promoter database.  
    sc_promoters = pd.read_pickle(base_dir + os.sep + os.path.join("expression_data","sc_promoters","sc_promoters.pkl"))

    # Make a subset of kl promoters 
    gene_subset_kl_orth = list(set(spec_comparision_data[spec_comparision_data['SC_common_name'].isin(gene_list)]['kl_genename']))
    #including set here removes duplicates

    kl_promoters_subset = kl_promoters.loc[gene_subset_kl_orth,]

    promoter_fname_kl = data_processing_dir + os.path.join("kl_promoters","promoter_sets", fname_prefix + "_kl.fasta")
    prom_seq_column = 4
    with open(promoter_fname_kl,'w') as f: 
        for row in kl_promoters_subset.itertuples():
            header_line = '>' + row[0] + ' 700bp_upstream\n'
            seq_line = row[prom_seq_column] + '\n'
            f.write(header_line)
            f.write(seq_line)

    # Make a subset of sc promoters for mitochondrial translation.  
    gene_subset_sc = list(set(spec_comparision_data[spec_comparision_data['SC_common_name'].isin(gene_list)]['sc_genename']))
    sc_promoters_subset = sc_promoters.loc[gene_subset_sc,]

    promoter_fname_sc = data_processing_dir + os.path.join("sc_promoters","promoter_sets", fname_prefix + "_sc.fasta")
    prom_seq_column = 2
    with open(promoter_fname_sc,'w') as f: 
        for row in sc_promoters_subset.itertuples():
            header_line = '>' + row[0] + ' 700bp_upstream\n'
            seq_line = row[prom_seq_column] + '\n'
            f.write(header_line)
            f.write(seq_line)
    
    return promoter_fname_kl, promoter_fname_sc


def run_fimo_command(promoter_fname, thresh, fname_prefix, output_dir, 
                     motif_fname = data_processing_dir + os.path.normpath('motifs/JASPAR_CORE_2016_fungi.meme'), 
                     motif = "All"):
    #Default is the Jaspar core 2016 database with all motifs
    #should eventually make motif command accept common name, look up number in dictionary, and accept more than one. 
    #test with this motif: MA0398.1
    #motif_dict = io_library.build_motif_dict(fname)
    #motif_dict['MA0371.1']
    
    if motif == "All":
        motif_arg = []
    else:
        motif_arg = ["--motif",motif]

    fimo_command = ([ "/home/kieran/meme/bin/fimo",
                      "--oc", output_dir,
                      "--verbosity", "1",
                      "--thresh", str(thresh)] +
                     motif_arg + 
                     [ motif_fname,
                       promoter_fname]
                   )

    fimo_output = subprocess.run(fimo_command,stdout = subprocess.PIPE) 

    print("fimo output return code = " + str(fimo_output.returncode))

    #change file prefix and delete output other than .txt file

    files_to_remove = [output_dir + os.sep + fname for fname in ['cisml.xml', 'fimo.html', 'fimo.xml', 'fimo.gff']]
    for file_to_remove in files_to_remove:
        os.remove(file_to_remove)


    fimo_fname_out = output_dir + os.sep + fname_prefix + '_fimo.txt'
    os.rename(output_dir + os.sep + 'fimo.txt', fimo_fname_out)
    
    return fimo_fname_out

def nucleotide_frequency(species):
    promoters = pd.read_pickle(base_dir + os.sep + os.path.join("expression_data",
                                                                species.lower()+ "_promoters",species.lower() + "_promoters.pkl"))
    actg_count = np.array([0,0,0,0])
    alphabet = 'ACGT'

    for gene, row in promoters.iterrows():
        prom_seq = row.loc['prom_seq']
        actg_count_gene = np.array([prom_seq.count(letter) for letter in alphabet])
        actg_count = actg_count + actg_count_gene

    nucleotide_frequency_out = actg_count/sum(actg_count)
    
    print("nucleotide frequency for " + alphabet + " = " + str(nucleotide_frequency_out))
    
    return nucleotide_frequency_out
    
def hex_color_dictionary(input_list):
    cmap = plt.get_cmap('hsv')
    N = len(input_list)
    cmap_discrete = [cmap(nn/N) for nn in range(0,N)]
    #cmap=plt.get_cmap('tab10')
    hex_colors = []

    for rgb_color in cmap_discrete:
        int_color = tuple(int(rgb_element*255) for rgb_element in rgb_color)
        hex_colors.append('#{:02x}{:02x}{:02x}'.format(*int_color))

    input_color_dict = dict(zip(input_list,hex_colors))
    
    return input_color_dict

def exact_promoter_scan_genelist(gene_list, motif_dict, promoter_database, output_format = 'counts', sequence_context = 0): 
    #finds nonoverlapping exact matches forward and backward for motifs. 
    #input:  motif dictionary, promoter data structure, gene list from dataframe (genes must be primary key of promoter data structure)
    #output_format: 'counts'(default) or 'full'.  If 'full' is selected each entry is a list of tuples containing the location of the 
    #               found motif, and the sequence (with sequence context) of the motif 
    #sequence_context: default 0.  How many bases on either side of the motif to display in the output. 
    #output: dataframe with primary key as gene list.  columns are presence of scanned motifs - either counts or location data. 

    output_data_frame = pd.DataFrame(index = gene_list)

    for motif_name, motif in motif_dict.items():
        L_motif = len(motif)
        output_motif = []
        for gene in gene_list: 

            #Make a sequence object
            promoter_database.loc[gene]
            prom_seq = Seq(promoter_database.loc[gene]['prom_seq'],alphabet=generic_dna)

            #fwd search
            prom_seq_fwd = str(prom_seq)
            L_prom = len(prom_seq_fwd)
            motif_sites_fwd = [m.start() for m in re.finditer(motif, prom_seq_fwd)]
            #to find overlapping motifs
            #[m.start() fr m in re.finditer('(?=' + motif + ')', prom_seq_rev)]

            #reverse search
            #rev search
            prom_seq_rev = str(prom_seq.reverse_complement())
            #location is in reference to the fwd seq
            motif_sites_rev = [len(prom_seq_rev) - m.start() - L_motif for m in re.finditer(motif, prom_seq_rev)]

            #This funtion doesn't do anything to locate the sites. 
            if output_format == 'count':  
                n_motifs_gene = len(motif_sites_fwd) + len(motif_sites_rev)
                output_motif.append(n_motifs_gene)
            elif output_format == 'full': 
                loc_hitseq_gene = []
                for motif_site in motif_sites_fwd: 
                    loc = L_prom-1-motif_site
                    hitseq = prom_seq_fwd[motif_site-sequence_context:motif_site+L_motif+sequence_context]       
                    loc_hitseq_gene.append((loc,hitseq))
                for motif_site in motif_sites_rev: 
                    loc = L_prom - motif_site+L_motif
                    hitseq = prom_seq_rev[(L_prom-motif_site-L_motif-sequence_context):(L_prom-motif_site+sequence_context)]
                    loc_hitseq_gene.append((loc,hitseq))
                
                if len(loc_hitseq_gene)==0: 
                    output_motif.append(None)
                else:
                    output_motif.append(loc_hitseq_gene)
            else: 
                print("choose output_format: 'count' or 'full' ")


        output_data_frame[motif_name] = output_motif    
    
    return output_data_frame


def exact_promoter_scan(motif, prom_seq, output_format='count', sequence_context=0):
    #this function is for a single sequence - same name used to be for one that searched a list of genes 
    # now called exact_promoter_scan_genelist
    # 
    # scans a sequence (biopython object) for a motif in both forward and backward directions 
    # If output format is 'count', then just lists the number of hits. 
    # If output format is 'full', gives the location, direction (relative to the gene), and
    # sequence context (number of bases on either side of the hit)output_format
    

    L_motif = len(motif)

    #fwd search
    prom_seq_fwd = str(prom_seq)
    L_prom = len(prom_seq_fwd)
    motif_sites_fwd = [m.start() for m in re.finditer(motif, prom_seq_fwd)]
    #to find overlapping motifs
    #[m.start() fr m in re.finditer('(?=' + motif + ')', prom_seq_rev)]

    #reverse search
    #rev search
    prom_seq_rev = str(prom_seq.reverse_complement())
    #location is in reference to the fwd seq
    motif_sites_rev = [len(prom_seq_rev) - m.start() - L_motif for m in re.finditer(motif, prom_seq_rev)]

    if output_format == 'count':  
        n_motifs_gene = len(motif_sites_fwd) + len(motif_sites_rev)
        output = n_motifs_gene
    elif output_format == 'full': 
        loc_hitseq_gene = []
        for motif_site in motif_sites_fwd: 
            loc = L_prom-1-motif_site
            hitseq = prom_seq_fwd[motif_site-sequence_context:motif_site+L_motif+sequence_context]  
            hitdir = 'fwd'
            loc_hitseq_gene.append((loc,hitdir,hitseq))
        for motif_site in motif_sites_rev: 
            loc = L_prom - motif_site+L_motif
            hitseq = prom_seq_rev[(L_prom-motif_site-L_motif-sequence_context):(L_prom-motif_site+sequence_context)]
            hitdir = 'rev'
            loc_hitseq_gene.append((loc,hitdir,hitseq))

        if len(loc_hitseq_gene)==0: 
            output = None
        else:
            output = loc_hitseq_gene
    else: 
        print("choose output_format: 'count' or 'full' ")
        
    return output



def threshold_sign(x, high_threshold, low_threshold): 
    if x > high_threshold:
        x_sign = 1
    elif x< low_threshold:
        x_sign = -1
    else: 
        x_sign = 0
        
    return x_sign

def threshold_group_SC(a,b,high_threshold, low_threshold):
    #input: a, b, high_threshold, low_threshold
    #output: group - either {'up_up', 'up_flat','up_down','flat_flat', 'down_flat','down_down'}

    a_sign = threshold_sign(a, high_threshold, low_threshold)
    b_sign = threshold_sign(b, high_threshold, low_threshold)

    sign_sum = a_sign+b_sign

    if sign_sum ==2:
        group = 'up_up'
    elif sign_sum == -2: 
        group = 'down_down'
    elif sign_sum == 1: 
        group = 'up_flat'
    elif sign_sum == -1: 
        group = 'down_flat'
    elif sign_sum == 0:
        if a_sign-b_sign == 0:
            group = 'flat_flat'
        else: 
            group = 'up_down'
    else: 
        print('Error - did not assign group')

    return group

def threshold_group_SC_series(A,B, high_threshold, low_threshold ):
    
    AB = pd.concat([A,B], axis = 1)

    group = []
    for ind,row in AB.iterrows():
        group.append(threshold_group_SC(row[0],row[1],high_threshold, low_threshold))

    group = pd.Series(group, index = AB.index)

    return group

def threshold_group(a,high_threshold, low_threshold):
    #change to this more general name from threshold_group_KL
    #input: a, b, high_threshold, low_threshold
    #output: group - either up, flat, or down

    a_sign = threshold_sign(a, high_threshold, low_threshold)

    if a_sign ==1:
        group = 'up'
    elif a_sign == 0: 
        group = 'flat'
    elif a_sign == -1: 
        group = 'down'
    else: 
        print('Error - did not assign group')

    return group

def threshold_group_series(A, high_threshold, low_threshold ):
    
    group = []
    for ind,item in A.iteritems():
        group.append(threshold_group(item,high_threshold, low_threshold))

    group = pd.Series(group, index = A.index)

    return group

def build_phylo_conversion_sc():
    #builds file for name conversion from phylome_db tree files. 
    #that file can be read in with pd.read_table()
    
    fname_gene_conversion = data_processing_dir + os.path.normpath("phylogenetic_trees/all_id_conversion.txt")
    fname_gene_conversion_sc = data_processing_dir + os.path.normpath("phylogenetic_trees/sc_id_conversion.txt")
    n_start = 2531
    n_stop = 24766
    
    with open(fname_gene_conversion) as gc_all:
        with open(fname_gene_conversion_sc,"w") as gc_sc:
            gene_name_options = []
            gene_name = 'init'
            for i, line in enumerate(gc_all):
                if i >= n_start-1:
                    #Get first line - gene_name
                    newline_gene_name = line.split()[0]
                    #print('gene_name: ' + gene_name)
                    #print('newline_gene_name:' + newline_gene_name)
                    if newline_gene_name == gene_name: 
                        #if new gene name matches last one, add option to list
                        gene_name_options.append(line.split())
                    else: 
                        #if new gene name doesn't match, 
                        
                        #Process last set of options (except on first run)
                        if len(gene_name_options)>0:
                            phylome_name = gene_name
                            gene_name_option_dict = {item[1]:item[2] for item in gene_name_options}
                            if "ensembl" in set(gene_name_option_dict.keys()):
                                ref_name = gene_name_option_dict["ensembl"]
                            elif "genename" in set(gene_name_option_dict.keys()):
                                ref_name = gene_name_option_dict["genename"]
                                print('no ensemble id for ' + gene_name + ', using genename: ' + ref_name)
                            else: 
                                print('no ensemble id or genename for ' + gene_name + ', using phylome DB id')
                                ref_name = phylome_name
                            gc_sc.write(phylome_name + '\t' + ref_name + '\n') 
                        
                        #Reset gene name and reset gene name options
                        gene_name = newline_gene_name
                        #print('next gene: ' + newline_gene_name)
                        gene_name_options = []                  
                        gene_name_options.append(line.split())
                    if i >= n_stop-1:
                        #Process final set of options
                        phylome_name = gene_name
                        gene_name_option_dict = {item[1]:item[2] for item in gene_name_options}
                        if "ensembl" in set(gene_name_option_dict.keys()):
                            ref_name = gene_name_option_dict["ensembl"]
                        elif "genename" in set(gene_name_option_dict.keys()):
                            ref_name = gene_name_option_dict["genename"]
                            print('no ensemble id for ' + gene_name + ', using genename: ' + ref_name)
                        else: 
                            print('no ensemble id or genename for ' + gene_name + ', using phylome DB id')
                            ref_name = phylome_name
                        gc_sc.write(phylome_name + '\t' + ref_name + '\n') 
                        print('Done')
                        break

    return

def get_phylomedb_name(sc_genename_list):
    sc_phylomedb_convert = pd.read_table(data_processing_dir + os.path.normpath("phylogenetic_trees/sc_id_conversion.txt"), header=None)
    sc_phylomedb_convert.columns = ["phylomedb_name","sc_genename"]

    phylomedb_list = []

    for gene in sc_genename_list:
        phylomedb_matches = sc_phylomedb_convert[sc_phylomedb_convert['sc_genename']==gene]['phylomedb_name']
        if len(phylomedb_matches) == 1:
            phylomedb_name = phylomedb_matches.values[0]
        elif len(phylomedb_matches) == 0:
            print('No match for ' + gene + ' in phylome namespace')
            phylomedb_name = None
        phylomedb_list.append(phylomedb_name)
    
    return phylomedb_list

def build_phylo_conversion_kl():
    #builds file for name conversion from phylome_db tree files for K.Lactis genes. 
    #that file can be read in with pd.read_table()
    fname_gene_conversion = data_processing_dir + os.path.normpath("phylogenetic_trees/all_id_conversion.txt")
    fname_gene_conversion_kl = data_processing_dir + os.path.normpath("phylogenetic_trees/kl_id_conversion.txt")
    n_start = 88796
    n_stop = 98798
    
    with open(fname_gene_conversion) as gc_all:
        with open(fname_gene_conversion_kl,"w") as gc_kl:
            gene_name_options = []
            gene_name = 'init'
            for i, line in enumerate(gc_all):
                if i >= n_start-1:
                    #Get first line - gene_name
                    newline_gene_name = line.split()[0]
                    #print('gene_name: ' + gene_name)
                    #print('newline_gene_name:' + newline_gene_name)
                    if newline_gene_name == gene_name: 
                        #if new gene name matches last one, add option to list
                        gene_name_options.append(line.split())
                    else: 
                        #if new gene name doesn't match, 
    
                        #Process last set of options (except on first run)
                        if len(gene_name_options)>0:
                            phylome_name = gene_name
                            gene_name_option_dict = {item[1]:item[2] for item in gene_name_options}
                            if "Genolevures" in set(gene_name_option_dict.keys()):
                                ref_name = gene_name_option_dict["Genolevures"]
                            elif "genename" in set(gene_name_option_dict.keys()):
                                ref_name = gene_name_option_dict["genename"]
                                print('no Genolevures id for ' + gene_name + ', using genename: ' + ref_name)
                            else: 
                                print('no Genolevures id or genename for ' + gene_name + ', using phylome DB id')
                                ref_name = phylome_name
                            gc_kl.write(phylome_name + '\t' + ref_name + '\n') 
    
                        #Reset gene name and reset gene name options
                        gene_name = newline_gene_name
                        #print('next gene: ' + newline_gene_name)
                        gene_name_options = []                  
                        gene_name_options.append(line.split())
                    if i >= n_stop-1:
                        #Process final set of options
                        phylome_name = gene_name
                        gene_name_option_dict = {item[1]:item[2] for item in gene_name_options}
                        if "Genolevures" in set(gene_name_option_dict.keys()):
                            ref_name = gene_name_option_dict["Genolevures"]
                        elif "genename" in set(gene_name_option_dict.keys()):
                            ref_name = gene_name_option_dict["genename"]
                            print('no Genolevures id for ' + gene_name + ', using genename: ' + ref_name)
                        else: 
                            print('no Genolevures id or genename for ' + gene_name + ', using phylome DB id')
                            ref_name = phylome_name
                        gc_kl.write(phylome_name + '\t' + ref_name + '\n') 
                        print('Done')
                        break
    return


def assign_topologies(trees, WGH_branch_seed, KLE_branch_seed, ZT_branch_seed, outgroups, verbose = False):
    topology_dict = {WGH_branch_seed: "C", KLE_branch_seed: "A", ZT_branch_seed:"B"}
    
    topologies = []
    for gene,row in trees.iterrows():
        #print(gene)
        tree_string = row['tree']
        #tree_string = trees.loc[gene]["tree"]
        #"(Phy000D0CO_YEAST:0.0377129,Phy000NQ6C_SACBA:0.0653462,((Phy004FKRY_NAUDC:0.0912444,Phy000NS92_SACCA:0.115243)0.99985:0.0765755,((((Phy004F4O4_588726:0.253337,(Phy004F7XM_1071379:0.214507,(Phy000JNPU_VANPO:0.153757,Phy004FCCE_TETPH:0.171821)0.99985:0.111541)0.998964:0.0332526)0.980573:0.0259834,(((Phy000JRAI_KLUWA:0.065763,Phy002423J_LACTH:0.0853684)0.99985:0.118867,(Phy000NWPH_SACKL:0.118849,(Phy0000C6I_ASHGO:0.163073,Phy0008PPH_KLULA:0.187974)0.999239:0.0291304)0.99985:0.036076)0.99985:0.0367641,(Phy004FVSC_TORDC:0.14757,Phy00244CV_ZYGRO:0.280471)0.99985:0.0455249)0.99985:0.0234786)0.952493:0.0167291,(Phy004FGQ7_KAZAF:0.157457,(Phy00042WG_CANGA:0.179752,(Phy003M12E_51660:0.219913,Phy003LT0W_51914:0.298372)0.99985:0.0535978)0.99985:0.0516007)0.953575:0.0196913)0.972812:0.0235719,((Phy000NM4H_SACBA:0.0797609,Phy000CVUB_YEAST:0.111285)0.99985:0.241968,(((Phy004G2AN_HANAN:0.607512,((Phy0043HAA_DEKBR:0.613043,(Phy0002J8L_CANAL:0.261718,(Phy000M23V_PICST:0.182477,Phy0005LWI_DEBHA:0.179178)0.99985:0.0567132)0.99985:0.131607)0.99985:0.0941789,(Phy000EWFJ_YARLI:0.646307,((Phy004D9WH_45786:0.295069,Phy004DAUY_45786:0.786718)0.99985:0.125919,(Phy000D148_SCHPO:0.930355,((Phy004D9WR_45786:0.198548,((Phy004G3IF_HANAN:0.174478,(Phy0000C5V_ASHGO:0.212532,(((Phy004F8NL_1071379:0.141367,(Phy004FD9R_TETPH:0.141431,Phy000JLBS_VANPO:0.116295)0.99985:0.0650686)0.99985:0.0356093,((((Phy004F1VH_588726:0.164592,(Phy000NUHY_SACCA:0.0565971,Phy004FI75_NAUDC:0.0592237)0.99985:0.032963)0.967128:0.0146652,(Phy003M2DA_51660:0.0718156,Phy003LSK0_51914:0.123436)0.99985:0.0352864)0.993589:0.0159514,(Phy000CY5Z_YEAST:0.0752616,(Phy00042VW_CANGA:0.0982544,(Phy004FH3D_KAZAF:0.0433095,Phy004FGT3_KAZAF:0.0291784)0.99985:0.0655249)0.998957:0.0272474)0.991005:0.0227564)0.99985:0.0179636,(Phy004FYCC_TORDC:0.0799393,Phy00244D9_ZYGRO:0.119371)0.994601:0.0197164)0.99985:0.0232753)0.99985:0.060257,(Phy0023ZJ5_SACKL:0.0790576,(Phy0008M3B_KLULA:0.284,(Phy002433O_LACTH:0.034923,Phy000JPRE_KLUWA:0.0533176)0.99985:0.0827682)0.99985:0.0294055)0.934264:0.0162554)0.99985:0.0331771)0.99985:0.0953258)0.99985:0.104003,(Phy0043EA7_DEKBR:0.329095,(Phy0002IUY_CANAL:0.113784,(Phy0005HV2_DEBHA:0.168247,Phy000M1Q7_PICST:0.0681382)0.2608:0.0373838)0.99985:0.142364)0.949433:0.0296014)0.99985:0.058341)0.99985:0.439321,((Phy000D147_SCHPO:0.337686,(Phy0043DGZ_DEKBR:0.179192,((Phy000NXUS_SACKL:0.0119718,(Phy000JQ10_KLUWA:0.00356342,Phy00240LJ_LACTH:0.0173874)0.99985:0.0420066)0.973509:0.0124443,(Phy000PHQX_KLULA:0.0837091,((Phy004FW12_TORDC:0.0270097,Phy00245C9_ZYGRO:0.0799646)0.99985:0.0293122,((Phy003LVON_51914:0.00488701,Phy003LVW9_51914:0.0214526)0.99985:0.049813,((Phy004FB62_TETPH:0.0450363,Phy000JO2T_VANPO:0.0157048)0.99985:0.0284698,((Phy004F55K_588726:0.043957,((Phy003ET24_SACBA:9.8e-08,Phy003ESPQ_SACBA:0.00517011)0.99985:0.0245326,(Phy000CZIO_YEAST:0.00368032,Phy000CW40_YEAST:0.00684232)0.99985:0.0147427)0.99985:0.0233092)0.99985:0.0129149,(((Phy004FL32_NAUDC:8.89e-08,Phy004FJ8W_NAUDC:0.0210664)0.99985:0.018709,(Phy000NR38_SACCA:0.00441582,Phy000NRWI_SACCA:0.0338758)0.99985:0.0200783)0.848919:0.00637556,(Phy000HL5M_CANGA:0.0380492,(Phy004FFGC_KAZAF:0.00245352,Phy004FHDM_KAZAF:0.00791846)0.99985:0.0383082)0.99985:0.0206773)0.929977:0.010802)0.99985:0.0134442)0.98884:0.0114619)0.997822:0.00886798)0.99985:0.0274939)0.99985:0.0182694)0.99985:0.146112)0.99985:0.27075)0.99985:0.616952,((Phy000NXBP_SACKL:0.0276937,((Phy0000DBL_ASHGO:0.0543494,(Phy000JPNA_KLUWA:0.0151433,Phy00243RQ_LACTH:0.00665339)0.99985:0.0350212)0.99985:0.02568,(Phy0008O77_KLULA:0.0552353,((Phy004FCXO_TETPH:0.0382825,(Phy000JO2F_VANPO:0.030324,Phy000JND8_VANPO:0.0667573)0.99985:0.0551445)0.99985:0.0188317,(Phy004F80P_1071379:0.0787646,((Phy004F1YX_588726:0.0731204,(Phy004FFWL_KAZAF:0.0212898,((Phy000CYMB_YEAST:0.0105865,Phy000NLRW_YEAST:0.0212602)0.99985:0.0105349,(Phy003ESNK_SACBA:0.0190025,Phy000NLWS_SACBA:0.0299389)0.99985:0.00986892)0.99985:0.0227824)0.97486:0.00873849)0.970244:0.00473995,((Phy000457A_CANGA:0.0545661,(Phy004FKDA_NAUDC:0.0379504,Phy000NR20_SACCA:0.0526524)0.953743:0.014551)0.994092:0.0148149,((Phy003LUN2_51914:0.0422204,Phy003M1LP_51660:0.036874)0.995219:0.011347,(Phy004FVH1_TORDC:0.0373219,Phy00244HG_ZYGRO:0.100079)0.99985:0.0412169)0.921431:0.00903586)0.998632:0.0123235)0.99985:0.0134243)0.985352:0.00847894)0.99985:0.0739419)0.99985:0.0195691)0.998527:0.0148528)0.99985:0.0749685,(Phy004G2XL_HANAN:0.14119,(((Phy000EWZ3_YARLI:0.0811337,(Phy000EUZU_YARLI:0.0243446,(Phy003G9YQ_YARLI:0.0484158,Phy000EXEP_YARLI:0.0328251)0.995874:0.0170054)0.99985:0.0245836)0.99985:0.0600747,(Phy004D8CO_45786:0.151996,(Phy004DBR1_45786:0.0452565,Phy004DANJ_45786:0.0687919)0.99985:0.028464)0.99985:0.0358608)0.99985:0.0642867,((Phy0002HWT_CANAL:0.0628017,(Phy0005M0G_DEBHA:0.102678,Phy000M2EQ_PICST:0.0437075)0.99985:0.0338248)0.99985:0.0940611,((Phy0002LM8_CANAL:0.0839271,Phy0005M8P_DEBHA:0.180701)0.99985:0.0500089,(((Phy001SWG1_SCHPO:0.0352369,Phy001SWG2_SCHPO:0.063444)0.99985:0.199495,(Phy0043E6V_DEKBR:0.103961,Phy0043H82_DEKBR:0.133998)0.99985:0.110704)0.99985:0.042867,(Phy004G3G7_HANAN:0.0982853,((Phy00240WA_LACTH:0.0451633,Phy000JR1W_KLUWA:0.0681361)0.99985:0.0587049,((Phy003FKGG_ASHGO:0.0811873,(Phy000NYEU_SACKL:0.0536055,Phy0008N6P_KLULA:0.0923059)0.77818:0.00986821)0.998707:0.0219914,(Phy004FVHE_TORDC:0.063156,(((Phy000NU01_SACCA:0.0890922,(((Phy004F8XZ_1071379:0.132533,(Phy003LSW3_51914:0.107247,Phy003M3EY_51660:0.0677824)0.99985:0.0412051)0.858533:0.00671357,(Phy004FKS7_NAUDC:0.0665065,(Phy000JLXL_VANPO:0.0625492,(Phy004FCBV_TETPH:0.145474,Phy000JNJG_VANPO:0.0837134)0.99985:0.034268)0.99985:0.0430787)0.998837:0.0208369)0.99985:0.0196212,(Phy000457W_CANGA:0.0711854,(Phy004F2YU_588726:0.0564193,(Phy003ET83_SACBA:0.0412039,Phy000CX4H_YEAST:0.0312003)0.99985:0.0837373)0.947318:0.0223552)0.941307:0.00984725)0.99985:0.011665)0.99985:0.0192451,(Phy004FH1A_KAZAF:0.0907382,(Phy000CVMU_YEAST:0.0230885,Phy000NM19_SACBA:0.0378338)0.99985:0.077243)0.993676:0.0240291)0.99985:0.0203553,(Phy00244KO_ZYGRO:0.227036,(Phy000NT8F_SACCA:0.0706802,Phy004FLUW_NAUDC:0.102082)0.99985:0.0400108)0.951593:0.0156376)0.821965:0.0142286)0.99985:0.0555068)0.995297:0.0215114)0.99985:0.0617542)0.99985:0.0766563)0.99985:0.0391873)0.99985:0.0512995)0.868152:0.00804878)0.986886:0.0281588)0.895443:0.0194272)0.99985:0.201861)0.987571:0.16679)0.99985:3.18348)0.99985:0.571344)0.871722:0.0680403)0.99985:0.176123)0.99985:0.0800853)0.99985:0.182843,(Phy004F7V6_1071379:1.13573,(Phy003M25Z_51660:1.38478,Phy003LWAK_51914:1.88325)0.99985:0.607834)0.985457:0.149056)0.99985:0.0802184,((Phy004FDRT_TETPH:0.413917,Phy000JMLE_VANPO:0.43574)0.99985:0.102264,(Phy004FLX8_NAUDC:0.23091,Phy000NSTN_SACCA:0.23412)0.99985:0.102889)0.337602:0.0337093)0.99985:0.0474337)0.924651:0.0456035)0.99985:0.0507224)0.99985:0.102691);"
    
        tree = Tree(tree_string)
    
        #Extract only leaves that are are WGH_branch_seed, KLE_branch_seed, ZT_branch seed 
        leaf_names = tree.get_leaf_names()
        leaf_subset = [name for name in leaf_names if name.split('_')[1] in {WGH_branch_seed,KLE_branch_seed,ZT_branch_seed}]
    
        #print(leaf_subset)
    
        #Check for 1-2 WGH proteins and only 1 each of KL and ST branch
        leaf_subset_spec = [name.split("_")[1] for name in leaf_subset]
        leaf_subset_spec_counter = Counter(leaf_subset_spec)
        
        for branch_seed in [KLE_branch_seed,ZT_branch_seed]:
            if (leaf_subset_spec_counter[branch_seed]>1):
                #Throw out extra WGH branch gene if children of common ancestor between it and current gene are in an outgroup.  
    
                #Prune down tree to include spec of interest and outgroups
                leaf_subset_outgroups = [name for name in leaf_names if name.split('_')[1] in set([WGH_branch_seed,KLE_branch_seed,ZT_branch_seed] + outgroups)]
                tree.prune(leaf_subset_outgroups) 
                #print(tree)
                #for each of the branch genes. 
                branch_genes = [leaf for ind,leaf in enumerate(leaf_subset) if leaf_subset_spec[ind]==branch_seed]
                for branch_gene in branch_genes:
                    ancestor = tree.get_common_ancestor([gene,branch_gene])
                    ancestor_leaves = ancestor.get_leaf_names()
                    ancestor_leaves_spec = [name.split("_")[1] for name in ancestor_leaves]
    
                    #if the leaves of the parent node include an outgroup, throw out the gene
                    if len(set(ancestor_leaves_spec) & set(outgroups)) > 0:
                        if verbose: 
                            print('Threw out ' + branch_gene + ' because outgroup gene between it and seed')
                        leaf_subset.remove(branch_gene)    
                        leaf_subset_spec = [name.split("_")[1] for name in leaf_subset]
                        leaf_subset_spec_counter = Counter(leaf_subset_spec)
    
                    if  (leaf_subset_spec_counter[branch_seed]==0):
                        if verbose: 
                            print('Alert threw out both ' + branch_seed + ' genes for seed: ' + gene)
    
        if (leaf_subset_spec_counter[KLE_branch_seed]==1) & (leaf_subset_spec_counter[WGH_branch_seed] > 0) & (leaf_subset_spec_counter[ZT_branch_seed]==1):
    
            #if 2 or more S.Cer proteins, remove any other paralogs from the leaf subset
            if leaf_subset_spec_counter[WGH_branch_seed]>1:
                if leaf_subset_spec_counter[WGH_branch_seed]>2:
                    if verbose: 
                        print(gene + " has more than 2 paralogs")
                        print(leaf_subset)
                not_scer = [spec != WGH_branch_seed for spec in leaf_subset_spec]
                is_gene = [leaf == gene for leaf in leaf_subset]
                is_gene_or_not_scer = [ig or ns for ig,ns in zip(is_gene, not_scer)]
                leaf_subset = [gene for gene,cond in zip(leaf_subset,is_gene_or_not_scer) if cond]
    
            tree.prune(leaf_subset)  
            children = [child.name for child in tree.get_children()]
            if len(children)==2:
                children.remove('')
                distant_child = children[0]
                distant_child_spec = distant_child.split("_")[1]
                topology = topology_dict[distant_child_spec]
            elif len(children)==3:
                if verbose: 
                    print("three children for " + gene + ", assume topology D, leaf subset =")
                    print(leaf_subset)
                topology = "D"
            else: 
                if verbose: print("Something must be wrong: " + gene)
        elif (leaf_subset_spec_counter[KLE_branch_seed]==1) & (leaf_subset_spec_counter[WGH_branch_seed] > 0) & (leaf_subset_spec_counter[ZT_branch_seed]==0):
            if verbose: 
                print(gene + " missing " + ZT_branch_seed + " protein")
                print(leaf_subset)
            topology = "missing " + ZT_branch_seed + " protein"
        elif (leaf_subset_spec_counter[KLE_branch_seed]==0) & (leaf_subset_spec_counter[WGH_branch_seed] > 0) & (leaf_subset_spec_counter[ZT_branch_seed]==1):
            if verbose: 
                print(gene + " missing " + KLE_branch_seed + " protein")
                print(leaf_subset)
            topology = "missing " + KLE_branch_seed + " protein"
        elif (leaf_subset_spec_counter['KLULA']==0) & (leaf_subset_spec_counter[WGH_branch_seed] > 0) & (leaf_subset_spec_counter[ZT_branch_seed]==0):
            if verbose: 
                print(gene + " missing " + ZT_branch_seed + " and " + KLE_branch_seed + " protein")
                print(leaf_subset)
            topology = "missing " + ZT_branch_seed + " and " + KLE_branch_seed + " protein"
        else: 
            if verbose: 
                print("non-standard phylogeny for " + gene + ", leaf subset =")
                print(leaf_subset)
            topology = "under construction"
    
    
        topologies.append(topology)
    
    return(topologies)

def assign_single_topology(gene,leaf_subset,tree, KLE_branch_seed, WGH_branch_seed, ZT_branch_seed, verbose):
    topology_dict = {WGH_branch_seed: "C", KLE_branch_seed: "A", ZT_branch_seed:"B"}
    leaf_subset_spec = [name.split("_")[1] for name in leaf_subset]
    leaf_subset_spec_counter = Counter(leaf_subset_spec)

    #if (leaf_subset_spec_counter[KLE_branch_seed]==1) & (leaf_subset_spec_counter[WGH_branch_seed] in {1,2}) & (leaf_subset_spec_counter[ZT_branch_seed]==1):
    if (leaf_subset_spec_counter[KLE_branch_seed]==1) & (leaf_subset_spec_counter[ZT_branch_seed]==1):
        #If there is 1 KLE, 1 ZT and 1 or 2 WGH, 
        #Check for topology for low and high
        #prune tree to three leaves 
        tree_three_leaf = tree.copy()
        
        tree_three_leaf.prune(leaf_subset)  
        children = [child.name for child in tree_three_leaf.get_children()]
        if len(children)==2:
            #remove '' because that is an internal node.
            children.remove('')  
            distant_child = children[0]
            distant_child_spec = distant_child.split("_")[1]
            #topology is assigned based on which of the three leaves are most distant. 
            topology_out = topology_dict[distant_child_spec]
        elif len(children)==3:
            if verbose: 
                print("three children for " + gene + ", assume topology D, leaf subset =")
                print(leaf_subset)
            topology_out = "D"
        else: 
            if verbose: 
                print("Something must be wrong: " + gene) 
    elif (leaf_subset_spec_counter[KLE_branch_seed]==0) & (leaf_subset_spec_counter[ZT_branch_seed]==0):
        if verbose: 
            print(gene + " missing " + ZT_branch_seed + " and " + KLE_branch_seed + " protein")
        topology_out = "missing " + ZT_branch_seed + " and " + KLE_branch_seed + " protein"
    elif leaf_subset_spec_counter[ZT_branch_seed]==0:
        if verbose: 
            print(gene + " missing " + ZT_branch_seed + " protein")
        topology_out = "missing " + ZT_branch_seed + " protein"
    elif (leaf_subset_spec_counter[KLE_branch_seed]==0):
        if verbose: 
            print(gene + " missing " + KLE_branch_seed + " protein")
        topology_out = "missing " + KLE_branch_seed + " protein"
    else: 
        if verbose: 
            print("non-standard phylogeny for " + gene + ", leaf subset =")
            print(leaf_subset)
        topology_out = "under construction"
    
    return topology_out


def assign_joint_topology(genes, leaf_subset, tree, single_topologies, paralog_labels, sc_genes, KLE_branch_seed, WGH_branch_seed, ZT_branch_seed, verbose):
    
    leaf_subset_spec = [name.split("_")[1] for name in leaf_subset]
    leaf_subset_spec_counter = Counter(leaf_subset_spec)
    
    branch_timing_check_dict = {"A": [ZT_branch_seed], "B": [KLE_branch_seed], "C": [ZT_branch_seed, KLE_branch_seed]}
    #In an N-N2 topology a late brancher is more closely related to the other species designated by this dictionary
    #than its paralog is.    

    #Make sure topology of either case is not a pathological case. 
    pathological_cases = {"missing " + ZT_branch_seed + " protein", 
                         "missing " + KLE_branch_seed + " protein",
                         "missing " + ZT_branch_seed + " and " + KLE_branch_seed + " protein",
                         "under construction",
                         None, 
                        "D"}

    topology_updated = single_topologies

    if len( {single_topologies[label] for label in paralog_labels}  & pathological_cases) > 0:
        joint_topology_out = 'pathological'
    elif single_topologies[paralog_labels[0]] != single_topologies[paralog_labels[1]]:
        joint_topology_list = list(single_topologies.values()) 
        joint_topology_list.sort()
        joint_topology_out = '-'.join(joint_topology_list)
    else: 
        sister_check = len(tree.get_common_ancestor(list(genes.values())).get_leaf_names()) 
        if sister_check == 2:    #If two WGH genes are most closely related to one another it is a 1 topology 
            joint_topology_out = single_topologies[paralog_labels[0]] + '-' + single_topologies[paralog_labels[1]] + '1'
            #both single topologies should be the same, but use different labels anyway.  If a see something like A-B1 I will know 
            #that is an issue. 
        elif sister_check in {3,4}: #2 topologies
            joint_topology_out = single_topologies[paralog_labels[0]] + '-' + single_topologies[paralog_labels[1]] + '2'

            #determine which paralog branched early or late.  late brancher is sister with other species.
            shared_topology = single_topologies[paralog_labels[0]]
            if shared_topology != single_topologies[paralog_labels[1]]:
                print('bad logic in joint topology assignment ' + str(sc_genes))
            for label in paralog_labels:
                #use early_branch_check_dict to look up the branch/branches that are closest to the WGH genes. 
                branch_timing_check_list = [leaf for ind, leaf in enumerate(leaf_subset) if leaf_subset_spec[ind] in branch_timing_check_dict[shared_topology]]
                #add in the gene we are testing
                branch_timing_check_list = branch_timing_check_list + [genes[label]]
                #late brancher if SC gene is more closely related to other group than its paralog
                if len(tree.get_common_ancestor(branch_timing_check_list).get_leaf_names()) == len(branch_timing_check_list):
                    #the length of the list of leaves is the same as the check list so this gene is more closely related than 
                    #its paralog and it is a late brancher
                    topology_updated[label] = single_topologies[label] + '_late'
                elif len(tree.get_common_ancestor(branch_timing_check_list).get_leaf_names()) == len(branch_timing_check_list)+1: 
                    #the length of the list of leaves is one longer than the check list so it must now include the other paralog. 
                    #this gene is therefore the early brancher. 
                    topology_updated[label] = single_topologies[label] + '_early'
                else:
                    print("problem :early branch check broken " + str(branch_timing_check_list))
        else:
            print('something went wrong ' + sc_genes[paralog_labels[0]] + ' ' + sc_genes[paralog_labels[1]] + ' sister check = ' + str(sister_check))
            print(tree)
     
    return joint_topology_out, topology_updated

def joint_topology_type(joint_topology): 
    joint_topologies_gene_conv = {'A-A1','B-B1','C-C1','C-C2'}
    joint_topologies_diff_lineages = {'A-A2','B-B2','A-C','A-B','B-C'}
    
    if joint_topology in joint_topologies_gene_conv:
        joint_topology_type = 'gene_conv'
    elif joint_topology in joint_topologies_diff_lineages:
        joint_topology_type = 'diff_lineages'
    else:
        print("Joint topology not in defined set: {}".format(joint_topology))
        joint_topology_type = None
    
    return joint_topology_type

def consistant_joint_topologies(ortholog_dataset, paralog_labels, resolution):
    #ortholog_dataset must have joint_topology_A_tree and joint_topology_B_tree columns for paralog labels (A,B) 
    ind_list = []
    for ind,row in ortholog_dataset.iterrows(): 
        joint_topology_A_tree = row['joint_topologies_'+paralog_labels[0]+'_tree']
        joint_topology_B_tree = row['joint_topologies_'+paralog_labels[1]+'_tree']
        #Make sure none of the joint topologies are pathological or None
        if len({None, 'pathological', 'C-outgroup-?'} - {joint_topology_A_tree,joint_topology_B_tree})==3:
            #check to see if both joint topologies agree
            if resolution == 'exact':
                if joint_topology_A_tree==joint_topology_B_tree:
                    ind_list.append(ind)
                    #print("{} and {} have matching joint topologies: {} and {} respectively".format(row['SC_common_name_low'],row['SC_common_name_high'],joint_topology['low'],joint_topology['high']))           
            elif resolution == 'type': 
                joint_topology_types = {'A':joint_topology_type(joint_topology_A_tree), 
                                       'B':joint_topology_type(joint_topology_B_tree)}
                if joint_topology_types['A']==joint_topology_types['B']:
                    ind_list.append(ind)
                    #print("{} and {} have similar joint topologies: {} and {} respectively".format(row['SC_common_name_low'],row['SC_common_name_high'],joint_topology['low'],joint_topology['high']))

    return(ind_list)

def get_other_paralogs_from_dataframe(genes, dataframe): 
    #using an input list of genes that have paralogs, and a dataframe that contains its paralogs, outputs a dataframe that
    #has just the other paralogs that are not listed. 
    #the genes must be listed as the SC_common_name, and the dataframe must have sc_genename and kl_genename fields. 
    N_genes_in = len(genes)
    dataframe_genes = dataframe[dataframe['SC_common_name'].isin(genes)]
    sc_genenames = dataframe_genes['sc_genename']
    kl_genenames = dataframe_genes['kl_genename']
    dataframe_all = dataframe[dataframe['kl_genename'].isin(kl_genenames)]
    sc_genenames_all = dataframe_all['sc_genename']
    sc_genenames_paralogs = list(set(sc_genenames_all)-set(sc_genenames))
    N_genes_out = len(sc_genenames_paralogs)
    if N_genes_in != N_genes_out: 
        print("number of genes in not equal to number of genes out - maybe one of the input genes doesn't have a paralog, or it is not contained in the dataframe") 
    dataframe_paralogs = dataframe_all[dataframe_all['sc_genename'].isin(sc_genenames_paralogs)]
    
    return(dataframe_paralogs)
    
    
    return 

