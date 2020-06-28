# -*- coding: utf-8 -*-
#from IPython.core.debugger import Tracer
print('Importing yeast_esr_exp.  If autoreload, may need to reset base_dir and data_processing dir \n  yeast_esr_exp.base_dir=base_dir \n yeast_esr_exp.data_processing_dir = data_processing_dir')
import os 
import pandas as pd
import numpy as np
import re
import math
import scipy.stats as stats
import scipy.spatial.distance as spd
from collections import Counter
import subprocess
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
import gffutils  
#from ete3 import Tree
#ete3 is not officially supported on windows, and so must be loaded via pip: 
# pip install -U https://github.com/etetoolkit/ete/archive/qt5.zip
# ref: https://groups.google.com/forum/#!topic/etetoolkit/6NblSBPij4o
#20181031: got this error message: twisted 18.7.0 requires PyHamcrest>=1.9.0, which is not installed.

# import requests   #Maybe I can remove these routines or do them in another file
# from lxml import etree    #parses xml output
from itertools import product
#import pickle  #Try not to use pickle
from itertools import chain

base_dir = ''
data_processing_dir = ''
#These need to be set after importing the module based on file structure 
#set in std_libraries.py
#I could probably do it automatically with relative paths. 

spec_lookup = {'Klac' : 'Kluyveromyces lactis', 'Scer': 'Saccharomyces cerevisiae', 
 'Cgla' : 'Candida glabrata' , 'Ncas': 'Naumovozyma castellii', 
 'Sbay' : 'Saccharomyces bayanus', 'Smik': 'Saccharomyces mikatae',
 'Lwal' : 'Lachancea waltii', 'Spar' : 'Saccharomyces paradoxus', 
 'Lklu' : 'Lachancea kluyverii', 'Dhan': 'Debaryomyces hansenii', 
 'Calb' : 'Candida albicans', 'Ylip': 'Yarrowia lipolytica', 
 'Sjap' : 'Schizosaccharomyces japonicus' , 'Spom' : 'Schizosaccharomyces pombe' }


## External Data
###SGD Data

def read_SGD_features():
    
    #Read in orf/name file and make it a dictionary
    # Gabe 7/12/16
    # SC_features_fname = os.path.normpath(data_processing_dir + "\ortholog_files\\SGD_features.tab")
    SC_features_fname = os.path.normpath(data_processing_dir + "/ortholog_files_regev/SGD_features.tab")

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

### Other external data

def make_platform_dict(spec):
    #For a given species, uses a gene expression dataset from the Tsankov et al Raw Expression data to build a map from 
    #gene ID to gene name.  Stores the output as a .csv
    
    id_source_data_inds = {'Klac': 'GSE22198', 'Scer': 'GSE22204', 
                            'Cgla':'GSE22194', 'Ncas' : 'GSE22200', 
                            'Sbay' : 'GSE22205', 'Suva' : 'GSE22205',  #Spar uses the same gene names as sbayanus - in this study they thought of it as a variant of bayanus.    
                            'Smik': 'GSE22201', 
                            'Lwal': 'GSE22199', 'Spar': 'GSE22193', 
                            'Lklu': 'GSE22202', 'Dhan': 'GSE22191', 
                            'Calb': 'GSE22197', 'Ylip': 'GSE22192'}
    orf_header_names = {spec_ohn : 'SystematicName' for spec_ohn in ['Klac', 'Cgla', 'Ncas', 'Sbay', 'Suva', 'Smik', 'Lwal', 'Spar','Dhan', 'Calb', 'Ylip']}
    orf_header_names['Scer']='ACCESSION_STRING'
    orf_header_names['Lklu'] = 'Gene'
    for spec_ohn in ['Vpol','Sjap','Spom']: 
        orf_header_names[spec_ohn] = 'ORF'
    
    orfmap_fname_out = os.path.normpath(data_processing_dir + 'regev_data/platform_ids/' + spec + '_ids.csv')
    
    if spec in ['Vpol','Sjap','Spom']: 
        orfmap_fname = os.path.normpath(data_processing_dir + 'regev_data/platform_ids/' + spec + '_ids.txt')
        platform_df = pd.read_table(orfmap_fname, comment='#', index_col=0)
        #Switch name for 'Vpol' to match orthogroups
        if spec == 'Vpol':
            vpol_rename = []
            for orf_name in platform_df.loc[:,orf_header_names[spec]]:
                kpol_ind = orf_name.split('Kpol')[1]
                orf_name_new = 'Kpol_' + kpol_ind.split('p')[0] + '.' + kpol_ind.split('p')[1]
                vpol_rename.append(orf_name_new)
            platform_df['YGOB_orf'] = vpol_rename
            platform_dict =  dict(platform_df.loc[:,'YGOB_orf'])
        else: 
            platform_dict =  dict(platform_df.loc[:,orf_header_names[spec]])
    
    else: 
    
        orfmap_fname = os.path.normpath(data_processing_dir + 'regev_data/platform_ids/' + id_source_data_inds[spec] + '_family.soft')

        with open(orfmap_fname) as f:
            # Find line that starts table listing gene names and index numbers
            for line in f: 
                if line.split()[0] == '!platform_table_begin':
                    break

            #Find index of header that will identify orf
            line = next(f)
            linesp = line.split()
            orf_header_ind = linesp.index(orf_header_names[spec])

            platform_dict = {}
            for line in f: 
                linesp = line.split('\t')
                if linesp[0] == "!platform_table_end\n":
                    #include new line because this splits on tab.  
                    break 

                orf_ind = linesp[0]
                #S. Cerevisiae orf names are formatted differently than other species. 
                if spec == 'Scer':                       
                    orf_name = linesp[orf_header_ind].split('|')[1]
                else: 
                    orf_name = linesp[orf_header_ind]
                
                #converts index to integer except in SCer where it is a string
                if spec == 'Scer': 
                    platform_dict[orf_ind] = orf_name  
                else:
                    platform_dict[int(orf_ind)] = orf_name

    
    pd.Series(platform_dict).to_csv(orfmap_fname_out)
    print(orfmap_fname_out + ' saved')
    
    return platform_dict

def parse_raw_exp(spec, platform_dict, save_file=True):
    #for a given species abbreviation and optional filename returns a dataframe with expression data from 
    #Tsankov et al 2010.   
    #Returns a dictionary which lists the platform mapping (from microarray to gene) and a dataframe of 
    #expression values, platform
    raw_exp_data_inds = {'Klac': 'GSE22198', 'Scer': 'GSE22204', 
                        'Cgla':'GSE22194', 'Ncas' : 'GSE22200', 
                        'Sbay' : 'GSE22205', 'Smik': 'GSE22201', 
                        'Lwal': 'GSE22199', 'Spar': 'GSE22193', 
                        'Lklu': 'GSE22202', 'Dhan': 'GSE22191', 
                        'Calb': 'GSE22197', 'Ylip': 'GSE22192'}
    
    #Load expression data
    raw_exp_fname = os.path.normpath(data_processing_dir + 'regev_data/raw_exp/' + raw_exp_data_inds[spec] + '_series_matrix.txt')
    raw_exp = pd.read_table(raw_exp_fname, comment = '!')

    #Make column for orf names
    orf_names = [platform_dict[platform_id] for platform_id in raw_exp['ID_REF']]
    raw_exp['orf_name'] = orf_names
    raw_exp.drop(columns='ID_REF', inplace = True)
    #Take median of all data that has more than one spot per orf 
    #Tsankov et al used median rather than mean so will stay with that. 
    grouped = raw_exp.groupby('orf_name')
    raw_exp_med = grouped.median()

    #Quantile normalize data across replicates and then calculate median of the normalized columns. 
    raw_exp_med_qnorm = quantileNormalize(raw_exp_med)
    raw_exp_med_qnorm['med_qnorm'] = raw_exp_med_qnorm.median(axis = 'columns')
    raw_exp_med_qnorm['std_qnorm'] = raw_exp_med_qnorm.std(axis = 'columns')
    raw_exp_out = raw_exp_med.merge(raw_exp_med_qnorm, left_index =True, right_index = True, how = 'outer', suffixes = ('','_qnorm'))

    ## Didn't seem to have a problem with nan's now that I am using pd.read_table.  If it comes up this might help: 
    # raw_exp[sample_datasets]= raw_exp[sample_datasets].applymap(lambda x: tryfloatconvert(x,np.nan))

    if save_file: 
        fname = data_processing_dir + os.path.normpath('regev_data/raw_exp/' + spec + '_raw_exp.csv')
        raw_exp_out.to_csv(fname)
        print(fname + ' saved')
    return raw_exp_out                    
    
def parse_micro_data(spec, exptype, platform_dict): 
    #load raw data for a given species/experiment type pair.  
    #
    #Inputs: 
    #    spec: species four letter abbreviation
    #    exptype: "Growth" or "Stress"
    #    platform_dict:  dictionary linking orf_name to spot id on the microarray.  comes from make_platform_dict
    #
    #Output: 
    #    expdata_sorted: dataframe with multindex with ID_REF and orf_name as levels for rows, and 
    #        columns for all the data with a multiindex 
    #        that contains condition, replicate name, and the array id for the data. 
    #
    
    
    exptype_prefixes = {'Growth': '36253', 'Stress': '38478'}
    platform_prefixes = {'Klac': '10499', 'Scer': '9294',
            'Cgla':'10497', 'Ncas' : '10501', 
            'Sbay' : '10505', 'Suva': '10505',  #Sbay and Suva in the same dataset - require special handling
            'Smik': '10502', 
            'Lwal': '10500', 'Spar': '10496', 
            'Lklu': '10503', 'Dhan': '15298', 
            'Vpol': '15297', 'Calb': '10498', 
            'Ylip': '10495', 'Sjap': '15299',
            'Spom': '15300' }

    input_fname = os.path.normpath(data_processing_dir + '/regev_data/GSE' + exptype_prefixes[exptype] + '_' + exptype + '/GSE'+ exptype_prefixes[exptype] + '-GPL' + platform_prefixes[spec] + '_series_matrix.txt')

    expdata = pd.read_table(input_fname, comment = '!', index_col=0)

    #if species is Scer drops everything after the index "DarkCorner" - those are a set of controls built into the microarray that we don't have data for and don't use. 
    if spec == 'Scer':
        expdata.drop(index=expdata.loc["DarkCorner":].index, inplace=True)


    # Get list of conditions and repicates 
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

    array_ids = list(expdata.columns)

    col_mindex_arrays = [conditions,replicates, array_ids]
    col_mindex_tuples = list(zip(*col_mindex_arrays))
    col_mindex = pd.MultiIndex.from_tuples(col_mindex_tuples, names=['conditions', 'replicates','array_ids'])
    expdata.columns = col_mindex
    #sort by conditions
    expdata_sorted = expdata.sort_index(level = 'conditions',axis = 'columns')
    #add in index of gene names 
    if not(list(expdata_sorted.index)==list(platform_dict.keys())):
        print("Error: ID mismatch between experiment data and orf lookup table. Species = {}, Experiment Type = {}".format(spec, exptype))
        return

    #if no error, continue here
    print("All ID's match between experiment data and orf lookup table. Species = {}, Experiment Type = {}".format(spec, exptype))

    orf_name = [platform_dict[spot_id] for spot_id in expdata_sorted.index]
    expdata_sorted['orf_name'] = orf_name
    
    #Add orf_name to index
    expdata_sorted.reset_index(inplace=True)
    expdata_sorted.set_index(['ID_REF', 'orf_name'], inplace=True)
    
    #Sbay and Suva are combined on the same microarray - this splits them apart. 
    if spec in {'Sbay', 'Suva'}: 
        if spec =='Sbay':
            conds_to_keep = [condition for condition in list(expdata_sorted.columns.levels[0]) if condition[0:6]!='sbayuv']
            expdata_sorted = expdata_sorted.loc[:,conds_to_keep]
        elif spec =='Suva': 
            conds_to_keep = [condition for condition in list(expdata_sorted.columns.levels[0]) if condition[0:6]=='sbayuv']
            new_cond_names = {condition: condition.split('_')[1] for condition in conds_to_keep}
            expdata_sorted = expdata_sorted.loc[:,conds_to_keep]
            expdata_sorted.rename(columns=new_cond_names, inplace=True)

    return expdata_sorted
 
def make_data_tables(spec_list):
    
    for spec in spec_list:  
        
        platform_dict = make_platform_dict(spec)
        if not(spec in {'Vpol', 'Sjap', 'Spom', 'Suva'}):     #Vpol, Sjap, Spom, and Suva don't have raw expression data
            #Generate raw expression data
            raw_exp = parse_raw_exp(spec, platform_dict, save_file=True)
        
        #Generate data for microarrays
        #files are saved within parse_raw_data - maybe it would be more consistent to do that in parse_micro_data as well. 
        #print(platform_dict)
        growth_exp = parse_micro_data(spec,'Growth',platform_dict)
        fname = os.path.normpath(data_processing_dir + "regev_data/GSE36253_Growth/"  + spec + '_growth.csv' )
        growth_exp.to_csv(fname)
        print(fname + ' saved')
        
        if (spec in {'Scer', 'Cgla', 'Calb', 'Klac', 'Lwal', 'Ncas', 'Sjap','Spom'}):
            #Stress data is only in these species
            stress_exp = parse_micro_data(spec,'Stress',platform_dict)
            fname = os.path.normpath(data_processing_dir + "regev_data/GSE38478_Stress/"  + spec + '_stress.csv' )
            stress_exp.to_csv(fname)
            print(fname + ' saved')
    
    return 

def combine_growth_stress_datasets(spec):
    #This function takes median of replicates and combines growth and stress data into single dataframe and also saves that as a .csv. 
    #Also takes median of spots that have the same gene ID.  Output has index which is gene names and column names for each condition. 

    fname = os.path.normpath(data_processing_dir + "regev_data/GSE36253_Growth/"  + spec + '_growth.csv' )
    growth_exp = pd.read_csv(fname,header = [0,1,2], index_col = [0,1])
    print(fname + ' growth microarray dataset loaded')

    #group by conditions and take median
    growth_replicate_groups = growth_exp.groupby(axis = 1, level = 'conditions')
    growth_exp_med = growth_replicate_groups.aggregate(np.median)

    if (spec in {'Scer', 'Cgla', 'Calb', 'Klac', 'Lwal', 'Ncas', 'Sjap','Spom'}):
        fname = os.path.normpath(data_processing_dir + "regev_data/GSE38478_Stress/"  + spec + '_stress.csv' )
        stress_exp = pd.read_csv(fname,header = [0,1,2], index_col = [0,1])
        print(fname + ' stress microarray dataset loaded')

        #group by condition and take median
        stress_replicate_groups = stress_exp.groupby(axis = 1, level = 'conditions')
        stress_exp_med = stress_replicate_groups.aggregate(np.median)

        #combine growth and stress average expression datasets. 
        if False in stress_exp_med.index==growth_exp_med.index:
            print("Error: ID mismatch between condition data. Species = {}".format(spec))
        growth_stress_data = pd.concat([growth_exp_med,stress_exp_med], axis = 1)
        
    else: 
        growth_stress_data = growth_exp_med
        
    #gets rid of ID index and takes median of any spots that represent the same gene
    growth_stress_data_gene_ind = growth_stress_data.groupby('orf_name').median()
    if len(growth_stress_data_gene_ind) < len(growth_stress_data) :
        print(spec + ' has some duplicate spots')

    fname_out = os.path.normpath(data_processing_dir + 'regev_data/' + spec + '_growth_stress.csv')  
    growth_stress_data_gene_ind.to_csv(fname_out)
    print('combined dataset saved as ' + fname_out )

    return growth_stress_data

def regev_ohnolog_expression_data_SSD_combine(exp_df, spec_sets, spec_conditions, combine_method = 'mean'):
    #Takes a datframe which consists gene expression data for various conditions from
    #Thompson et al and Roy et al listed by their relationship to pairs of S.Cer WGH 
    #paralogs. Converts data for other duplications to single floating point numbers. 
    #
    #exp_df: input dataframe.  For consistent naming 
    #
    #spec_conditions: dictionary listing desired conditions for each species. 
    #
    #spec_sets: Dictionary linking name of species sets to a list of species to use.  Of the form: 
    #
    #spec_sets = {'Post WGH low' : [], 
    #         'Post WGH high' : [], 
    #         'Pre WGH' : []} 
    #
    # combine method:  Right now I only have mean
    # mean: set the value to the mean fold change between all duplicate orthologs. 
    #
    #other good options: 
        # A: set the value to the mean fold change between all
        # B: exclude that gene [start here]
        # C: Pick the one that is most highly expressed in raw data
        # D: Pick the one with highest fold change
        # E: Pick the one with lowest fold change


    levels = {'Post WGH low': 'low', 
              'Post WGH high': 'high', 
              'Pre WGH' : ''}

    expression_data = {}

    for row in exp_df.iterrows(): 

        high_gene_common_name = row[1]['SC_common_name_high'] 
        low_gene_common_name = row[1]['SC_common_name_low']
        
        expression_data_row = []
        
        #If there are no orthologs, set the value to np.Nan 
        for spec_set_name, spec_set in spec_sets.items():
            level = levels[spec_set_name]
            if level == '': 
                level_sep = ''
            else: 
                level_sep = '_'

            for spec in spec_set: 
                for condition in spec_conditions[spec]: 
                    vals = row[1][spec + '_' + condition + level_sep + level]
                    if isinstance(vals,list):   #if data isn't a list, it should be an np.nan
                        Nval = sum(np.logical_not(np.isnan(vals)))
                        if Nval > 0:  #if at least one item has a value
                            vals_clean = [val for val in vals if not(np.isnan(val))]
                            if combine_method == 'mean': #Change this to alter behavior for multiple orthologs
                                val_est = np.mean(vals_clean)  
                            expression_data_row.append(val_est)
                        else:  #Both paralogs had no value
                            expression_data_row.append(np.nan)
                    else: 
                        assert np.isnan(vals), ("Value " + str(vals) + " for " + spec + '_' + condition + ' is not a list, but not a nan')
                        expression_data_row.append(vals)


        expression_data[high_gene_common_name + '_' + low_gene_common_name] = expression_data_row

    columns = []
    for spec_set_name, spec_set in spec_sets.items():
        level = levels[spec_set_name]
        if level == '': 
            level_sep = ''
        else: 
            level_sep = '_'
        for spec in spec_set: 
            spec_columns = [spec + '_' + condition + level_sep + level for condition in list(spec_conditions[spec])]
            columns = columns + spec_columns
    expression_data_df = pd.DataFrame.from_dict(expression_data, orient = 'index', columns = columns)

    return expression_data_df

def load_regev_data_gois(ohnologs_goi, sort_column, self_spec, spec_order_post_WGH, spec_order_pre_WGH):
    #Load data for all species for a given set of gois
    #ohnologs_goi is indexed by YGOB ancestor, and has genename_low and genename_high referring to the level 
    #of expression in the sort_column (e.g. log2FoldChange or PKAest)

    ## SMIK data often seems to be missing

    #Makes dataframe with data for all experiments

    #Must use ortholog mapping for regev data. 
    orth_dir = data_processing_dir + 'ortholog_files_regev' + os.sep 

    columns_to_keep=[]
    for level in ['low', 'high']:
        for column_base in ['genename', sort_column]:
            columns_to_keep.append(column_base + '_' + level)

    ohnologs_goi_array = ohnologs_goi.loc[:,columns_to_keep]
    ohnologs_goi_array.rename(columns = {'genename_'+level : self_spec+'_genename_'+level for level in ['low','high']}, inplace=True)

    spec_conditions = {}

    #adding data just for self_spec

    #load data
    fname_array_data = os.path.normpath(data_processing_dir + 'regev_data/' + self_spec + '_growth_stress_norm.csv')  
    spec_data = pd.read_csv(fname_array_data, index_col=0)
    conditions = spec_data.columns
    spec_conditions[self_spec] = conditions

    for level in ['low','high']:
        condition_data_dict = {condition: [] for condition in conditions}
        for gene in ohnologs_goi_array[self_spec + '_genename_' + level ]:
            for condition in conditions:
                cond_data = []
                try: 
                    cond_data.append(spec_data.loc[gene,condition])
                except KeyError:
                    cond_data.append(np.nan)
                    print('Mismatch between goi index and expression data index ' + self_spec + ' : ' + gene + ' ' + condition)
                condition_data_dict[condition].append(cond_data)

        for condition in conditions: 
            ohnologs_goi_array[self_spec + '_' + condition + '_' + level] = condition_data_dict[condition]



    #post WGH:


    for spec in spec_order_post_WGH: 

        #load ortholog mapping
        orth_lookup = read_orth_lookup_table(self_spec, spec, orth_dir)

        #load data
        fname_array_data = os.path.normpath(data_processing_dir + 'regev_data/' + spec + '_growth_stress_norm.csv')  
        spec_data = pd.read_csv(fname_array_data, index_col=0)
        conditions = spec_data.columns
        spec_conditions[spec] = conditions

        for level in ['low','high']: 
            N_orth_list = []
            orth_name_list = []
            condition_data_dict = {condition: [] for condition in conditions}
            for gene in ohnologs_goi_array[self_spec + '_genename_' + level ]:
                try: 
                    if orth_lookup[gene][0]=='NONE':
                        N_orth_list.append(0)
                        orth_name_list.append('NONE_in_orthogroups')
                        for condition in conditions: 
                            condition_data_dict[condition].append(np.nan)
                        print(gene + ' has NONE listed in ortholog file for ' + spec)
                    else: 
                        orth_names = orth_lookup[gene]
                        N_orth_list.append(len(orth_names))
                        orth_name_list.append(orth_names)
                        for condition in conditions:
                            cond_data = []
                            for spec_gene in orth_names: 
                                try: 
                                    cond_data.append(spec_data.loc[spec_gene,condition])
                                except KeyError:
                                    cond_data.append(np.nan)
                                    print('Mismatch between ortholog file and expression data index ' + spec + ' : ' + spec_gene + ' ' + condition)
                            condition_data_dict[condition].append(cond_data)
                except KeyError:
                    N_orth_list.append(0)
                    orth_name_list.append('NONE_not_in_orthogroups')
                    for condition in conditions: 
                        condition_data_dict[condition].append(np.nan)
                    print(gene + ' not present in orthogroup file for ' + spec)

            ohnologs_goi_array[spec + '_N_' + level] = N_orth_list
            ohnologs_goi_array[spec + '_genename_' + level] = orth_name_list
            for condition in conditions: 
                ohnologs_goi_array[spec + '_' + condition + '_' + level] = condition_data_dict[condition]



    #pre WGH:


    for spec in spec_order_pre_WGH: 

        #load ortholog mapping
        orth_lookup = read_orth_lookup_table(self_spec, spec, orth_dir)

        #load data
        fname_array_data = os.path.normpath(data_processing_dir + 'regev_data/' + spec + '_growth_stress_norm.csv')  
        spec_data = pd.read_csv(fname_array_data, index_col=0)
        conditions = spec_data.columns
        spec_conditions[spec] = conditions

        N_orth_list = []
        orth_name_list = []
        condition_data_dict = {condition: [] for condition in conditions}

        for row in ohnologs_goi_array.loc[:,[self_spec + '_genename_low',self_spec + '_genename_high']].iterrows():

            gene_low = row[1][self_spec+'_genename_low']
            gene_high = row[1][self_spec+'_genename_high']

            try:         
                if orth_lookup[gene_low] == orth_lookup[gene_high]: 
                    gene_to_test = gene_low
                    if orth_lookup[gene_to_test][0]=='NONE':
                        N_orth_list.append(0)
                        orth_name_list.append('NONE_in_orthogroups')
                        for condition in conditions: 
                            condition_data_dict[condition].append(np.nan)
                        print(gene_to_test + ' has NONE listed in ortholog file for ' + spec)
                    else: 
                        orth_names = orth_lookup[gene_to_test]
                        N_orth_list.append(len(orth_names))
                        orth_name_list.append(orth_names)
                        for condition in conditions:
                            cond_data = []
                            for spec_gene in orth_names:
                                try: 
                                    cond_data.append(spec_data.loc[spec_gene,condition])
                                except KeyError:
                                    cond_data.append(np.nan)
                                    print('Mismatch between ortholog file and expression data index ' + spec + ' : ' + spec_gene + ' ' + condition)
                            condition_data_dict[condition].append(cond_data)

                else: 
                    raise LookupError('high and low paralogs do not have same ortholog in ' + spec + ': ' + gene_low + ' and ' + gene_high )
            except KeyError: 
                print('either ' + gene_low + ' or ' + gene_high + 'is not in index for ortholog table in ' + spec)    


                gene_to_test = 'NONE'
                if gene_low in orth_lookup.keys():
                    print(gene_low + ' is in index for lookup table')
                    gene_to_test = gene_low
                elif gene_high in orth_lookup.keys():
                    print(gene_high + ' is in index for lookup table')
                    gene_to_test = gene_high

                if gene_to_test != 'NONE':
                    print(gene_to_test + ' did have an ortholog in ' + spec)
                    if orth_lookup[gene_to_test][0]=='NONE':
                        N_orth_list.append(0)
                        orth_name_list.append('NONE_in_orthogroups')
                        for condition in conditions: 
                            condition_data_dict[condition].append(np.nan)
                        print(gene_to_test + ' has NONE listed in ortholog file for ' + spec)
                    else: 
                        orth_names = orth_lookup[gene_to_test]
                        N_orth_list.append(len(orth_names))
                        orth_name_list.append(orth_names)
                        for condition in conditions:
                            cond_data = []
                            for spec_gene in orth_names:
                                try: 
                                    cond_data.append(spec_data.loc[spec_gene,condition])
                                except KeyError:
                                    cond_data.append(np.nan)
                                    print('Mismatch between ortholog file and expression data index ' + spec + ' : ' + spec_gene + ' ' + condition)
                            condition_data_dict[condition].append(cond_data)
                else:
                    raise RuntimeError('Neither gene is index for ortholog table, but it was not caught in the first section')


        ohnologs_goi_array[spec + '_N'] = N_orth_list
        ohnologs_goi_array[spec + '_genename'] = orth_name_list
        for condition in conditions: 
            ohnologs_goi_array[spec + '_' + condition] = condition_data_dict[condition]


    return ohnologs_goi_array, spec_conditions

def sort_conservaton_by_col_avg(expression_data_df, spec_order, columns_to_combine, induced_thresh ):
    #for a given expression_data_df (e.g. from load_regev_data_gois), uses conditions in columns_to_combine 
    #for each species in spec order to sort rows by induction above a certain threshold (induced_thresh)
    #in the species at the end of the spec_order (which is presumed to be ordered phylogenetically)
    #Make a column which is average fold change for these columns
    for spec in ['Scer'] + spec_order: 
        columns_to_combine_spec = [column for column in expression_data_df.columns if 
         (column.split('_')[0]==spec) and 
         (column.split('_')[1] in columns_to_combine) and
         (column.split('_')[2]=='high')]
        expression_data_df[spec + '_PKAest_high'] = expression_data_df.loc[:,columns_to_combine_spec].apply(np.mean, axis = 1)
        expression_data_df[spec + '_PKAestInd_high'] = (expression_data_df[spec + '_PKAest_high']>induced_thresh)

    expression_data_df_sorted = expression_data_df.sort_values(by=[spec + '_PKAestInd_high' for spec in reversed(spec_order)], ascending=False)
    return expression_data_df_sorted

def sort_regev_stress_conditions(all_conds, spec_sets, expression_data_df):
    #Sorts out conditions in a combined expression_data_df (see load_regev_data_gois) 
    #so that all are in the order all_conds when present.  The different species sets are given by 
    #spec sets which is a dictionary of the form: 
    #
    #spec_sets = {'Post WGH low' : [], 
    #         'Post WGH high' : [], 
    #         'Pre WGH' : []} 

    levels = {'Post WGH low': 'low', 
              'Post WGH high': 'high', 
              'Pre WGH' : ''}

    new_col_order = []
    for spec_set_name in ['Post WGH high', 'Post WGH low', 'Pre WGH']:
        spec_set = spec_sets[spec_set_name]
        level = levels[spec_set_name]
        if level == '': 
            level_sep = ''
        else: 
            level_sep = '_'

        for spec in spec_set:    

            if level_sep == '':
                columns_spec = [column for column in expression_data_df.columns if 
                 (column.split('_')[0]==spec)]
            else: 
                columns_spec = [column for column in expression_data_df.columns if 
                 (column.split('_')[0]==spec) and 
                 (column.split('_')[-1]==level)]

            columns_spec_sorted = []
            for cond in all_conds: 
                cond_spec = spec + '_' + cond + level_sep + level
                if cond_spec in columns_spec: 
                    columns_spec_sorted.append(cond_spec)

            new_col_order = new_col_order + columns_spec_sorted

    return new_col_order

def de_stress_gois(spec, columns_to_combine, goi_criteria):
    #Gets differential expression set for a given species based on data from growth/stress microarryas from 2013 Broad study
    #
    #Inputs: 
    #    spec: species of interest 4 letter abbreviation (e.g. Scer)
    #    columns_to_combine:  Columns from growth_stress_norm.csv data that will be combined to form
    #                         PKAest, an estimate of induction under PKA inhibition.  It will not
    #                         Throw an error if those columns do not exist in the dataset (should probably add a warning)
    #    goi_criteria:  dictionary of goi criteria with the following fields: 
    #        min_high_lfc:  Minimum LFC for the high activation paralog
    #        lfc_diff:  Difference in LFC to call a differentially expressed paralog
    #        max_low_lfc: Maximum LFC that low activation paralog can have.  
    #
    #Required data:  Requires a normalized data file: <spec>_growth_stress_norm.csv
    #
    #Outputs: 
    #    ohnologs_goi: dataframe whose intex in the YGOB ancestor and with 'low' and 'high' paralogs
    #                  based on the value of PKAest.  Only contains genes that meet crietria for differential expression
    #    ohnologs_sorted: superset of ohnologs_goi, but contains all ohnologs from a given species
    #    
    #
    #Get all WGH paralogs pairs for a given species.  

    ohnologs = get_WGH_pairs_by_spec(spec)
    
    #For species that have different names in YGOB (for ohnolog names) than in expression data (regev names), 
    #makes genename_gene1 and genename_gene2 the ohnolog name
    if spec in {'Ncas', 'Smik'}: 
        print('Species is ' + spec + ', translating gene names from YGOB to regev')
        
        orth_dir_regev = data_processing_dir + "ortholog_files_regev" + os.sep

        YGOB_regev = read_orth_lookup_table(spec + '_YGOB',spec + '_regev', orth_dir_regev)

        ohnologs.rename(columns={'genename_gene1':'genename_YGOB_gene1', 'genename_gene2': 'genename_YGOB_gene2'}, inplace=True)

        for geneN in ['gene1', 'gene2']: 
            genename_regev_list = []
            for genename_YGOB in ohnologs['genename_YGOB_' + geneN]:
                genenames_regev = YGOB_regev[genename_YGOB]
                if len(genenames_regev)==1:
                    if genenames_regev[0]=='NONE':
                        print(genename_YGOB + ' has no match in regev dataset')
                    genename_regev_list.append(genenames_regev[0])
                else: 
                    print(genename_YGOB + ' has multiple matches in regev dataset: ')
                    print(genenames_regev)
                    genename_regev_list.append('MULTIPLE')

            ohnologs['genename_' + geneN] = genename_regev_list

    
    
    #load data
    fname_array_data = os.path.normpath(data_processing_dir + 'regev_data/' + spec + '_growth_stress_norm.csv')  
    spec_data = pd.read_csv(fname_array_data, index_col=0)
    conditions = spec_data.columns
    
    #threshold for induction using the average of the combined columns that make up the PKA inhibition estimate.  
    #All data is normalized already
    
    #Make a column which is average fold change for these columns
    columns_to_combine_spec = [column for column in conditions if column in columns_to_combine]
    spec_data['PKAest'] = spec_data.loc[:,columns_to_combine_spec].apply(np.mean, axis = 1)


    #Merge in expression data for all those paralogs, sort on average/individual conditions
    sort_column = 'PKAest'
    
    ohnologs_sorted = join_ohnologs_and_sort(spec_data, ohnologs, sort_column)

    ohnologs_goi = ohnologs_sorted[(ohnologs_sorted[sort_column + '_high']>goi_criteria['min_high_lfc']) &
                                          ((ohnologs_sorted[sort_column + '_high']-ohnologs_sorted[sort_column +'_low'])>goi_criteria['lfc_diff']) &
                                          (ohnologs_sorted[sort_column + '_low'] < goi_criteria['max_low_lfc'])]

    print(str(len(ohnologs_goi)) + ' differentially expressed genes identified')
    
    return ohnologs_goi, ohnologs_sorted, goi_criteria


def load_tsankov_data_gois(ohnologs_goi, sort_column, seed_spec, spec_order_post_WGH, spec_order_pre_WGH):
    #Load data for all species for a given set of gois
    #ohnologs_goi is indexed by YGOB ancestor, and has genename_low and genename_high referring to the level 
    #of expression in the sort_column (e.g. log2FoldChange or PKAest)
    #The seed_spec is assumed to be a post WGH species
    #the spec order for post WGH does not include the seed_spec

    orth_dir = data_processing_dir + 'ortholog_files_regev' + os.sep 

    columns_to_keep=[]
    for level in ['low', 'high']:
        for column_base in ['genename', sort_column]:
            columns_to_keep.append(column_base + '_' + level)

    ohnologs_goi_array = ohnologs_goi.loc[:,columns_to_keep]
    ohnologs_goi_array.rename(columns = {'genename_'+level : seed_spec+'_genename_'+level for level in ['low','high']}, inplace=True)

    #adding data just for seed_spec

    #load expression data
    fname_array_data = data_processing_dir + os.path.normpath('regev_data/raw_exp/' + seed_spec + '_raw_exp_norm.csv')  
    raw_exp_tsankov_norm = pd.read_csv(fname_array_data, index_col=0, header=None)
    spec_data = raw_exp_tsankov_norm.to_dict()[1]

    for level in ['low','high']: 
        raw_exp_tsankov_norm_level = []
        for gene in ohnologs_goi_array[seed_spec + '_genename_' + level ]: 
            try: 
                raw_exp_tsankov_norm_level.append([spec_data[gene]]) 
            except KeyError:
                print(gene + ' has no data in tsankov experiment for ' + seed_spec)
                raw_exp_tsankov_norm_level.append('NONE_not_in_tsankov_index')
        ohnologs_goi_array[seed_spec + '_raw_exp_tsankov_' + level] = raw_exp_tsankov_norm_level


    for spec in spec_order_post_WGH: 
        #load ortholog mapping
        seed_spec_orth_lookup = read_orth_lookup_table(seed_spec, spec, orth_dir)

        #Load expression data
        fname_array_data = data_processing_dir + os.path.normpath('regev_data/raw_exp/' + spec + '_raw_exp_norm.csv')  
        raw_exp_tsankov_norm = pd.read_csv(fname_array_data, index_col=0, header=None)
        spec_data = raw_exp_tsankov_norm.to_dict()[1]

        for level in ['low','high']: 
            raw_exp_tsankov_norm_level = []
            for gene in ohnologs_goi_array[seed_spec + '_genename_' + level ]: 
                try:
                    spec_genes = seed_spec_orth_lookup[gene]
                except KeyError:
                    print(gene + ' is not in orthogroup index ' + spec)
                    raw_exp_tsankov_norm_level.append('NONE_not_in_orthogroup_index')
                    continue
                if spec_genes[0]=='NONE':
                    raw_exp_tsankov_norm_level.append('NONE_no_ortholog')    
                else: 
                    raw_exp_list = []
                    for spec_gene in spec_genes: 
                        try:
                            raw_exp_val = spec_data[spec_gene]
                            raw_exp_list.append(raw_exp_val)
                        except KeyError: 
                            print(spec_gene + ' has no entry in Tsankov dataset for ' + spec + '. Ortholog of ' + gene)
                    if raw_exp_list==[]:
                        print('Neither of the orthologs for ' + gene + ' have an entry in Tsankov dataset')
                        raw_exp_tsankov_norm_level.append('NONE_not_in_tsankov_index')
                    else: 
                        raw_exp_tsankov_norm_level.append(raw_exp_list)
            ohnologs_goi_array[spec + '_raw_exp_tsankov_' + level] = list(raw_exp_tsankov_norm_level)


    for spec in spec_order_pre_WGH: 
        #load ortholog mapping
        seed_spec_orth_lookup = read_orth_lookup_table(seed_spec, spec, orth_dir)

        #Load expression data
        fname_array_data = data_processing_dir + os.path.normpath('regev_data/raw_exp/' + spec + '_raw_exp_norm.csv')  
        raw_exp_tsankov_norm = pd.read_csv(fname_array_data, index_col=0, header=None)
        spec_data = raw_exp_tsankov_norm.to_dict()[1]

        raw_exp_tsankov_norm_spec = []

        #this is currently set up without assuming that the ortholog of post WGH species with a similar ancestor will map 
        #to the same pre-WGH ortholog.  If the ortholog files are from pillars this is unneccessarily complex. 
        for row in ohnologs_goi_array.loc[:,[seed_spec + '_genename_low',seed_spec + '_genename_high']].iterrows():
            gene_low = row[1][seed_spec + '_genename_low']
            gene_high = row[1][seed_spec + '_genename_high']

            gene_to_test=None
            try:         
                spec_genes_low = seed_spec_orth_lookup[gene_low]
                gene_to_test=gene_low
            except KeyError: 
                print('low gene' + gene_low + ' is not in orthogroup index ' + spec)
                spec_genes_low = None
            try: 
                spec_genes_high = seed_spec_orth_lookup[gene_high]
                gene_to_test=gene_high
            except KeyError: 
                print('high gene' + gene_high + ' is not in orthogroup index ' + spec)            
                spec_genes_low = None

            if gene_to_test==None: 
                print('Niether high or low gene in orthogroup index ' + gene_low + ' '  + gene_high + ' for ' + spec)
                raw_exp_tsankov_norm_spec.append('NONE_not_in_orthogroup_index')
                continue

            if spec_genes_low == spec_genes_high:   #this should only be true if both of the genes have the same set of orthologs
                assert gene_to_test ==gene_high, 'If both genes have same orthologs gene_high will always be the gene to test'
                assert spec_genes_low!=None, "If both genes are None, then neither should be in the orthogroup index"
                spec_genes = seed_spec_orth_lookup[gene_to_test]
                if spec_genes[0]=='NONE':   
                    raw_exp_tsankov_norm_spec.append('NONE_no_ortholog')
                else: 
                    raw_exp_list = []
                    for spec_gene in spec_genes: 
                        try:
                            raw_exp_val = spec_data[spec_gene]
                            raw_exp_list.append(raw_exp_val)
                        except KeyError: 
                            print(spec_gene + ' has no entry in Tsankov dataset for ' + spec + '. Ortholog of ' + gene)
                    if raw_exp_list==[]:
                        print('Neither of the orthologs for ' + gene + ' have an entry in Tsankov dataset')
                        raw_exp_tsankov_norm_spec.append('NONE_not_in_tsankov_index')
                    else: 
                        raw_exp_tsankov_norm_spec.append(raw_exp_list)
            elif (spec_genes_low==None) or (spec_genes_high==None):  #Case for when one of the seed_spec genes does not have an ortholog but the other one does
                spec_genes = seed_spec_orth_lookup[gene_to_test]
                if spec_genes[0]=='NONE':   
                    raw_exp_tsankov_norm_spec.append('NONE_no_ortholog')
                else: 
                    raw_exp_list = []
                    for spec_gene in spec_genes: 
                        try:
                            raw_exp_val = spec_data[spec_gene]
                            raw_exp_list.append(raw_exp_val)
                        except KeyError: 
                            print(spec_gene + ' has no entry in Tsankov dataset for ' + spec + '. Ortholog of ' + gene)
                    if raw_exp_list==[]:
                        print('Neither of the orthologs for ' + gene + ' have an entry in Tsankov dataset')
                        raw_exp_tsankov_norm_spec.append('NONE_not_in_tsankov_index')
                    else: 
                        raw_exp_tsankov_norm_spec.append(raw_exp_list) 
            else: 
                raise ValueError('Orthologs of high and low paralogs of ' + gene_low + ' and ' + gene_high + ' do not match for ' + spec + '. ' + str(spec_genes_low) + str(spec_genes_high))

        ohnologs_goi_array[spec + '_raw_exp_tsankov'] = raw_exp_tsankov_norm_spec

    return ohnologs_goi_array

def tsankov_ohnolog_expression_data_SSD_combine(exp_df, spec_sets, combine_method = 'mean'):
    #Takes a datframe which consists of raw expression data for various conditions from
    #Tsankov et al listed by their relationship to pairs of WGH
    #paralogs. Converts data for other duplications to single floating point numbers. 
    #
    #exp_df: input dataframe.  For consistent naming 
    #
    #spec_sets: Dictionary linking name of species sets to a list of species to use.  Of the form: 
    #
    #spec_sets = {'Post WGH low' : [], 
    #         'Post WGH high' : [], 
    #         'Pre WGH' : []} 
    #
    # combine method:  Right now I only have mean
    # mean: set the value to the mean fold change between all duplicate orthologs. 
    #
    #other good options: 
        # A: set the value to the mean fold change between all
        # B: exclude that gene [start here]
        # C: Pick the one that is most highly expressed in raw data
        # D: Pick the one with highest fold change
        # E: Pick the one with lowest fold change

    levels = {'Post WGH low': 'low', 
                  'Post WGH high': 'high', 
                  'Pre WGH' : ''}

    raw_expression_data = {}

    for row in exp_df.iterrows(): 
        high_gene_common_name = row[1]['SC_common_name_high'] 
        low_gene_common_name = row[1]['SC_common_name_low']

        raw_expression_data_row = []

        for spec_set_name in ['Post WGH high', 'Post WGH low', 'Pre WGH']: 
            spec_set = spec_sets[spec_set_name]
            level = levels[spec_set_name]
            if level == '': 
                level_sep = ''
            else: 
                level_sep = '_'

            for spec in spec_set: 
                vals = row[1][spec + '_raw_exp_tsankov' + level_sep + level]
                if isinstance(vals,list):   #if data isn't a list, it should be an np.nan
                    #unlike with regev data I don't include none values for raw exp if orthologs exist. 
                    if combine_method == 'mean': #Change this to alter behavior for multiple orthologs
                        val_est = np.mean(vals)  
                        raw_expression_data_row.append(val_est)
                    else:  
                        raise ValueError('Bad Combine Method: ' + combine_method)
                else: 
                    assert isinstance(vals, str), ("Value " + str(vals) + " for " + spec + ' is not a list, but not a string')
                    raw_expression_data_row.append(np.nan)

        raw_expression_data[high_gene_common_name + '_' + low_gene_common_name] = raw_expression_data_row

    columns = []
    for spec_set_name in ['Post WGH high', 'Post WGH low', 'Pre WGH']: 
        spec_set = spec_sets[spec_set_name]
        level = levels[spec_set_name]
        if level == '': 
            level_sep = ''
        else: 
            level_sep = '_'
        for spec in spec_set: 
            columns.append(spec + level_sep + level)

    raw_expression_data_df = pd.DataFrame.from_dict(raw_expression_data, orient = 'index', columns = columns) 


    return raw_expression_data_df    

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
    title_list = '"treatment: 0 min: No 1-NM-PP1"   "treatment: 10 min: 3 µM 1-NM-PP1"  "treatment: 20 min: 3 µM 1-NM-PP1"  "treatment: 30 min: 3 µM 1-NM-PP1"  "treatment: 40 min: 3 µM 1-NM-PP1"  "treatment: 50 min: No 1-NM-PP1"    "treatment: 60 min: No 1-NM-PP1"    "treatment: 70 min: No 1-NM-PP1"    "treatment: 0 min: No 1-NM-PP1" "treatment: 10 min: 3 µM 1-NM-PP1"  "treatment: 20 min: 3 µM 1-NM-PP1"  "treatment: 30 min: No 1-NM-PP1"    "treatment: 40 min: No 1-NM-PP1"    "treatment: 50 min: No 1-NM-PP1"    "treatment: 60 min: No 1-NM-PP1"    "treatment: 70 min: No 1-NM-PP1"    "treatment: 0 min: No 1-NM-PP1" "treatment: 10 min: 120 nM 1-NM-PP1"    "treatment: 20 min: 120 nM 1-NM-PP1"    "treatment: 30 min: No 1-NM-PP1"    "treatment: 40 min: No 1-NM-PP1"    "treatment: 50 min: No 1-NM-PP1"    "treatment: 60 min: No 1-NM-PP1"    "treatment: 70 min: No 1-NM-PP1"    "treatment: 0 min: No 1-NM-PP1" "treatment: 5 min: 750 nM 1-NM-PP1" "treatment: 10 min: No 1-NM-PP1"    "treatment: 15 min: 750 nM 1-NM-PP1"    "treatment: 20 min: No 1-NM-PP1"    "treatment: 25 min: 750 nM 1-NM-PP1"    "treatment: 30 min: No 1-NM-PP1"    "treatment: 35 min_chip_1: 750 nM 1-NM-PP1" "treatment: 35 min_chip 2: 750 nM 1-NM-PP1" "treatment: 40 min: No 1-NM-PP1"    "treatment: 45 min: 750 nM 1-NM-PP1"    "treatment: 50 min: No 1-NM-PP1"    "treatment: 55 min: 750 nM 1-NM-PP1"    "treatment: 60 min: No 1-NM-PP1"    "treatment: 65 min: No 1-NM-PP1"    "treatment: 70 min: No 1-NM-PP1"    "treatment: 0 min: No 1-NM-PP1" "treatment: 5 min: 750 nM 1-NM-PP1" "treatment: 10 min: No 1-NM-PP1"    "treatment: 15 min: No 1-NM-PP1"    "treatment: 20 min: No 1-NM-PP1"    "treatment: 25 min: 750 nM 1-NM-PP1"    "treatment: 30 min_chip_1: No 1-NM-PP1" "treatment: 35 min: No 1-NM-PP1"    "treatment: 30 min_chip_2: No 1-NM-PP1" "treatment: 40 min: No 1-NM-PP1"    "treatment: 45 min: 750 nM 1-NM-PP1"    "treatment: 50 min: No 1-NM-PP1"    "treatment: 55 min: No 1-NM-PP1"    "treatment: 60 min: No 1-NM-PP1"    "treatment: 65 min: No 1-NM-PP1"    "treatment: 70 min: No 1-NM-PP1"'
    title_list = title_list.split('\t')
    title_list = [item.strip("\"") for item in title_list]
    id_ref_list = '"GSM812516"  "GSM812517" "GSM812518" "GSM812519" "GSM812520" "GSM812521" "GSM812522" "GSM812523" "GSM812524" "GSM812525" "GSM812526" "GSM812527" "GSM812528" "GSM812529" "GSM812530" "GSM812531" "GSM812532" "GSM812533" "GSM812534" "GSM812535" "GSM812536" "GSM812537" "GSM812538" "GSM812539" "GSM812540" "GSM812541" "GSM812542" "GSM812543" "GSM812544" "GSM812545" "GSM812546" "GSM812547" "GSM812548" "GSM812549" "GSM812550" "GSM812551" "GSM812552" "GSM812553" "GSM812554" "GSM812555" "GSM812556" "GSM812557" "GSM812558" "GSM812559" "GSM812560" "GSM812561" "GSM812562" "GSM812563" "GSM812564" "GSM812565" "GSM812566" "GSM812567" "GSM812568" "GSM812569" "GSM812570" "GSM812571"'
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

def parse_data_series_matrix_kemmeren(desired_conditions, series_fname, microarray_family_fname): 
    #Extracts data as a pandas dataframefor desired conditions from a series matrix file similar to the kemmeren gene deletion dataset
    #
    #Desired condition is a dictionary in which: 
    #   the key is the name you want to have for your columns (need not match file)
    #   the value is the array designator (should match file)
    #
    # series data is the series_matrix.txt file from NCBI for the data series. 
    #
    #microarray_family name is the name of the file with the table of microarray ids. 
    #

    microarray_lookup = get_microarray_lookup(microarray_family_fname)
    ids = set(microarray_lookup.keys())

    with open(series_fname) as f:
        #Find line that starts table listing gene names and index numbers
        for line in f: 
            if line.split('\t')[0] == '!series_matrix_table_begin\n':
                break

        exp_line = next(f)
        exp_list = [item.strip('"') for item in exp_line.split()]

        #extract data for only the listed conditions
        array_inds = {condition: exp_list.index(array) for condition,array in desired_conditions.items()}

        id_list = []
        orf_list = []
        exp_value_list = []

        for line in f:
            linesp = line.split('\t') 
            if linesp[0] == '!series_matrix_table_end\n':
                break       

            spot_id = int(linesp[0])
            if spot_id in ids:
                id_list.append(spot_id)
                orf = microarray_lookup[spot_id]
                orf_list.append(orf)
                exp_values = [float(linesp[array_ind]) for array_ind in array_inds.values()]
                exp_value_list.append(exp_values)

    #Make dataframe with orf values as index. 
    data = pd.DataFrame(exp_value_list, index = id_list, columns = array_inds.keys())

    data['sc_genename']=orf_list
    #take mean of columns with duplicated values. 
    grouped = data.groupby('sc_genename')
    data_mean = grouped.mean()

    return data_mean

def get_microarray_lookup(fname,nrows_skip=17):
    #given the filename of a full micrarray table from geo, make a lookup table for the genes using the
    # SYSTEMATIC_NAME column
    microarray_info = pd.read_table(fname, skiprows=nrows_skip)
    #only keeps rows for 'gene' reporter group (throws out controls) 
    microarray_info = microarray_info[microarray_info['REPORTER_GROUP']=='gene']
    microarray_lookup = dict(zip(microarray_info['ID'], microarray_info['SYSTEMATIC_NAME']))
    
    return microarray_lookup

## Id conversions
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

def SC_common_name_lookup_KL_generate_dict():
    #Builds database for SC_common_name_lookup_KL.  
    #Only need to run once. The database is stored as a csv file after that
    
    #Load ortholog data to add SCer common names for the KL genes.  For adds both orthologs separated by
    # an _.  
    # For some reason this takes a while.  SC_common_name_lookup should probably be optimized
    orth_dir = os.path.normpath(data_processing_dir + 'ortholog_files_YGOB') + os.sep
    orth_lookup = read_orth_lookup_table('Klac', 'Scer', orth_dir)

    sc_common_name_labels = {}

    for kl_gene, orth_list in orth_lookup.items():
        #print(kl_gene)
        label = '_'.join([SC_common_name_lookup([sc_orf])[0] for sc_orf in orth_list])
        sc_common_name_labels[kl_gene] = label
        #print(kl_gene + ':' + label)
        
    fname = os.path.normpath(data_processing_dir + 'ortholog_files_YGOB/kl_sc_common_name.csv')
    pd.Series(sc_common_name_labels).to_csv(fname)
    
    print("sc_common_name_labels dict saved at " + fname)
    #pickle.dump(sc_common_name_labels, open(os.path.normpath(data_processing_dir + 'ortholog_files_YGOB/kl_dict_sc_common_name.pkl'),'wb'))
  
    return sc_common_name_labels

def SC_common_name_lookup_KL(kl_genename_list): 
    #Return a list of SC common names.  If there is a paralog the paralog names are separated by an underscore.  
    #If there is no ortholog in SC, returns the KL name.  

    #Load labels for K.Lac genes
    
    #previous version used pickle
    #sc_common_name_labels = pickle.load(open(os.path.normpath(data_processing_dir + 'ortholog_files_YGOB/kl_dict_sc_common_name.pkl'),mode='rb'))
    
    #Load the csv and save as a dictionary
    fname = os.path.normpath(data_processing_dir + 'ortholog_files_YGOB/kl_sc_common_name.csv')
    sc_common_name_labels = pd.read_csv(fname, index_col = 0)['0'].to_dict()
    
    sc_common_name_label_list = []
    for kl_gene in kl_genename_list:
        try: 
            sc_common_name_label = sc_common_name_labels[kl_gene]
            if sc_common_name_label == 'NONE':
                sc_common_name_label = kl_gene

        except KeyError:
            sc_common_name_label = kl_gene
        sc_common_name_label_list.append(sc_common_name_label)

    return sc_common_name_label_list 

def SC_common_name_columns_ohnologs(ohnologs_sorted, spec, orth_dir):
    #Adds SC_common_name_high, low and high_low columns to a file with ohnologs from a species in YGOB
    #Input:
    #    ohnologs_sorted: dataframe with ohnologs sorted into low and high based on some criteria.  Need
    #        to have columns "genename_low" and "genename_high"
    #    spec: four letter abbreviation of species.
    #    orth_dir:  directory of ortholog file
    #         YGOB: data_processing_dir + "ortholog_files_YGOB" + os.sep
    #         regev: data_processing_dir + "ortholog_files_regev" + os.sep
    #         y1000: data_processing_dir + "ortholog_files_y1000" + os.sep

    if spec=='Scer': 
        for level in ['low', 'high']:
            ohnologs_sorted['SC_common_name_' + level] = SC_common_name_lookup(ohnologs_sorted['genename_' + level])
    else: 
        spec_scer_lookup = read_orth_lookup_table(spec, 'Scer', orth_dir)

        for level in ['low', 'high']: 
            SC_common_names = []
            for spec_gene in list(ohnologs_sorted['genename_' + level]): 
                scer_genes = spec_scer_lookup[spec_gene]
                if scer_genes == ['NONE']:
                    scer_genes = [spec_gene]
                scer_common_names = SC_common_name_lookup(scer_genes)
                #if there are two orthologs for one SC gene, they are joined with an underscore here
                SC_common_names.append('_'.join(scer_common_names))
            ohnologs_sorted['SC_common_name_' + level] = SC_common_names    
    
    ohnologs_sorted['SC_common_name_high_low']=ohnologs_sorted['SC_common_name_high'] + '_' + ohnologs_sorted['SC_common_name_low']
    
    return ohnologs_sorted

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

# def kl_genename_convert(df): 
    # #probably don't really need this one - just list one below. 
    # #Input is dataframe with GFF version of kl genename as the index. 
    # #Ouput is dataframe with standard version of kl genename as the index. 

    # kl_genename = kl_genename_convert_list(df.index)
    # df['kl_genename'] = kl_genename
    # df.set_index('kl_genename',inplace = True)

    # return df

def kl_genename_convert_list(kl_genes): 
    kl_genename = []
    for gene in kl_genes: 
        if gene[0:5]=='KLLA0':
            new_gene = gene.split('_')[0]+gene.split('_')[1]
        else: 
            new_gene = gene
        kl_genename.append(new_gene)

    return kl_genename
    
## Dealing with Orthologs
def write_YGOB_orth_lookup_table(spec1, spec2):
    ##!! Post WGH to post WGH is incorrect here - can't just use Pillars!!
    #Makes an orhtolog lookup table for all species in the YGOB Pillars data. 
    #spec1 the index, spec2 is the lookup.  multiple orthologs separated by tabs, NONE if there are none.  
    #Note that this cannot lookup Small Scale Duplications after the WGH because that data is not part of the pillars 
    #SSDs in pillars are just listed on their own - often without ancestors.  

    #older version required an input file - assuming we are always using YGOB_pillars_bmh.txt
    #older version didn't have all the cases and also returned a list instead of a dictionary

    fname = os.path.normpath(data_processing_dir + "ortholog_files_YGOB/YGOB_Pillars_bmh.txt")

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

    #position in pillars file of ancestor
    #anc_pos = 12

    #for regev lab expression data, Suva and Sbay have same gene names. 
    orth_positions_post_WGH = {'Vpol': [0,32], 'Tpha': [1,31], 'Tbla': [2,30], 'Ndai': [3,29], 
                      'Ncas': [4,28], 'Knag': [5,27], 'Kafr': [6,26], 'Cgla': [7,25],
                      'Suva': [8,24], 'Sbay':[8,24], 'Skud': [9,23], 'Smik': [10,22], 'Scer' : [11,21]}
                      

    orth_positions_pre_WGH = {'Zrou':13, 'Tdel':14, 'Klac':15, 'Egos':16, 'Ecym': 17, 'Lklu':18, 'Lthe':19, 'Lwal':20}

    #Case1: spec1 post WGH, spec2 pre WGH

    if ((spec1 in orth_positions_post_WGH.keys()) & (spec2 in orth_positions_pre_WGH.keys())): 
        spec1_columns = orth_positions_post_WGH[spec1]
        spec2_column = orth_positions_pre_WGH[spec2]

        with open(fname) as f:
            orth_lookup = {}
            for line in f:
                linesp = line.split()
                for spec1_column in spec1_columns: 
                    spec1_gene = linesp[spec1_column]
                    if spec1_gene!= '---':
                        spec2_gene = linesp[spec2_column]
                        if spec2_gene!='---': 
                            orth_lookup[spec1_gene] = spec2_gene
                        else: 
                            orth_lookup[spec1_gene] = 'NONE'
        
    #Case2: spec1 pre WGH, spec2 pre WGH

    if ((spec1 in orth_positions_pre_WGH.keys()) & (spec2 in orth_positions_pre_WGH.keys())): 
        spec1_column = orth_positions_pre_WGH[spec1]
        spec2_column = orth_positions_pre_WGH[spec2]

        with open(fname) as f:
            orth_lookup = {}
            for line in f:
                linesp = line.split()
                spec1_gene = linesp[spec1_column]
                if spec1_gene!= '---':
                    spec2_gene = linesp[spec2_column]
                    if spec2_gene!='---': 
                        orth_lookup[spec1_gene] = spec2_gene
                    else: 
                        orth_lookup[spec1_gene] = 'NONE'

        


    #Case 3: spec1 pre WGH, spec2 post WGH
    #In this case there will often be more than one ortholog per pre WGH gene

    if ((spec1 in orth_positions_pre_WGH.keys()) & (spec2 in orth_positions_post_WGH.keys())): 
        spec1_column = orth_positions_pre_WGH[spec1]
        spec2_columns = orth_positions_post_WGH[spec2]

        with open(fname) as f:
            orth_lookup = {}
            for line in f:
                linesp = line.split()
                spec1_gene = linesp[spec1_column]
                if spec1_gene!= '---':
                    spec2_genes = []
                    for spec2_column in spec2_columns: 
                        spec2_gene = linesp[spec2_column]
                        if spec2_gene!='---':
                            spec2_genes.append(spec2_gene)
                    if len(spec2_genes)==0:
                        orth_lookup[spec1_gene] = 'NONE'
                    else: 
                        orth_lookup[spec1_gene] = spec2_genes


    #Case 4: spec1 pre WGH, spec2 pre WGH
    #use tracks in pillars file to match up WGH genes.  If there are two paralogs, this will match tracks. 
    #If there is one paralog, and there are two orthologs, it will list both.  
    #If there is only one ortholog for one gene, it will list the other one regardless of track.  

    if ((spec1 in orth_positions_post_WGH.keys()) & (spec2 in orth_positions_post_WGH.keys())): 
        print(spec1 + ' and ' + spec2 + ' are both post WGH species and the pillars file does not assign syntenic orthologs.  Use write_YGOB_orth_lookup_table_goi')
        return
    #         spec1_columns = orth_positions_post_WGH[spec1]
    #         spec2_columns = orth_positions_post_WGH[spec2]

    #         with open(fname) as f:
    #             orth_lookup = {}
    #             for line in f:
    #                 linesp = line.split()
    #                 spec1_genes = []
    #                 for spec1_column in spec1_columns: 
    #                     spec1_gene = linesp[spec1_column]
    #                     if spec1_gene!= '---':
    #                         spec1_genes.append(spec1_gene)

    #                 if len(spec1_genes)==2: 
    #                     for jj,spec1_gene in enumerate(spec1_genes): 
    #                         spec2_column = spec2_columns[jj]
    #                         spec2_gene = linesp[spec2_column]
    #                         if spec2_gene!='---':
    #                             orth_lookup[spec1_gene] = [spec2_gene]
    #                         else: 
    #                             orth_lookup[spec1_gene] = ['NONE']

    #                 elif len(spec1_genes)==1: 
    #                     spec1_gene = spec1_genes[0]
    #                     spec2_genes = []
    #                     for spec2_column in spec2_columns: 
    #                         spec2_gene = linesp[spec2_column]
    #                         if spec2_gene!='---':
    #                             spec2_genes.append(spec2_gene)
    #                     if len(spec2_genes)==0:
    #                         orth_lookup[spec1_gene] = 'NONE'
    #                     else: 
    #                         orth_lookup[spec1_gene] = spec2_genes





    orth_lookup_outputfname = os.path.normpath(data_processing_dir + 'ortholog_files_YGOB/' + spec1 + "-" + spec2 + "-orthologs.txt"  )

    print(orth_lookup_outputfname)
    print_ortholog_file(orth_lookup_outputfname, orth_lookup)
                    
    return orth_lookup 

def get_WGH_pairs_by_spec(spec):
    #for a given post WGH species, return dataframe with Ancestor as key and pos1/pos2 orthologs. 

    fname = os.path.normpath(data_processing_dir + "ortholog_files_YGOB/YGOB_Pillars_bmh.txt")

    #position in pillars file of ancestor
    anc_pos = 12

    orth_positions = {'Vpol': [0,32], 'Tpha': [1,31], 'Tbla': [2,30], 'Ndai': [3,29], 
                      'Ncas': [4,28], 'Knag': [5,27], 'Kafr': [6,26], 'Cgla': [7,25],
                      'Suva': [8,24], 'Skud': [9,23], 'Smik': [10,22], 'Scer' : [11,21]}
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

    p1_col = orth_positions[spec][0]
    p2_col = orth_positions[spec][1]
    with open(fname) as f:
        paralog_list = []
        for line in f:
            linesp = line.split()
            if ((linesp[p1_col]!= '---') & (linesp[p2_col]!= '---')):
                paralog_list.append([linesp[anc_pos],linesp[p1_col],linesp[p2_col]])

    ohnologs = pd.DataFrame(paralog_list)
    ohnologs.columns = ['anc', 'genename_gene1', 'genename_gene2']
    ohnologs.set_index('anc', inplace=True)

    return ohnologs

def write_YGOB_orth_lookup_table_goi(spec1, spec2):
    #This is necessary because pillars does not specify syntenic assignment for post-WGH genes
    #uses data generated from YGOB in 20190708_DEpka_post_WGH_YGOB_synteny_assign

    #Makes an orhtolog lookup table for all species in the YGOB Pillars data for just goi
    #spec1 the index, spec2 is the lookup.  multiple orthologs separated by tabs, NONE if there are none.  

    #Note that this cannot lookup Small Scale Duplications after the WGH because that data is not part of the pillars 


    fname = os.path.normpath(data_processing_dir + "ortholog_files_YGOB/YGOB_pillars_goi.txt")

    #YGOB_Pillars.txt order of species: 
    #    0    V. polyspora Position B
    #    1    T. phaffii Position B
    #    2    T. blattae Position B
    #    3    N. dairenensis Position B
    #    4    N. castellii Position B
    #    5    K. naganishii Position B
    #    6    K. africana Position B
    #    7    C. glabrata Position B
    #    8    S. bayanus var. uvarum Position B
    #    9    S. kudriavzevii Position B
    #    10   S. mikatae Position B
    #    11   S. cerevisiae Position B
    #    12   L. waltii
    #    13   L. thermotolerans
    #    14   L. kluyveri
    #    15   E. cymbalariae
    #    16   E. gossypii 
    #    17   K. lactis
    #    18   T. delbrueckii
    #    19   Z. rouxii
    #    20   Ancestral Gene Order
    #    21   S. cerevisiae Position A
    #    22   S. mikatae Position A
    #    23   S. kudriavzevii Position A
    #    24   S. bayanus var. uvarum Position A
    #    25   C. glabrata Position A
    #    26   K. africana Position A
    #    27   K. naganishii Position A
    #    28   N. castellii Position A
    #    29   N. dairenensis Position A
    #    30   T. blattae Position A
    #    31   T. phaffii Position A
    #    32   V. polyspora Position A

    #position in pillars file of ancestor
    #anc_pos = 20

    #for regev lab expression data, Suva and Sbay have same gene names. 
    orth_positions_post_WGH = {'Vpol': [0,32], 'Tpha': [1,31], 'Tbla': [2,30], 'Ndai': [3,29], 
                      'Ncas': [4,28], 'Knag': [5,27], 'Kafr': [6,26], 'Cgla': [7,25],
                      'Suva': [8,24], 'Sbay':[8,24], 'Skud': [9,23], 'Smik': [10,22], 'Scer' : [11,21]}


    orth_positions_pre_WGH = {'Zrou':19, 'Tdel':18, 'Klac':17, 'Egos':16, 'Ecym': 15, 'Lklu':14, 'Lthe':13, 'Lwal':12}

    #Case1: spec1 post WGH, spec2 pre WGH

    if ((spec1 in orth_positions_post_WGH.keys()) & (spec2 in orth_positions_pre_WGH.keys())): 
        spec1_columns = orth_positions_post_WGH[spec1]
        spec2_column = orth_positions_pre_WGH[spec2]

        with open(fname) as f:
            orth_lookup = {}
            for line in f:
                linesp = line.split()
                for spec1_column in spec1_columns: 
                    spec1_gene = linesp[spec1_column]
                    if spec1_gene!= '---':
                        spec2_gene = linesp[spec2_column]
                        if spec2_gene!='---': 
                            orth_lookup[spec1_gene] = spec2_gene
                        else: 
                            orth_lookup[spec1_gene] = 'NONE'

    #Case2: spec1 pre WGH, spec2 pre WGH

    if ((spec1 in orth_positions_pre_WGH.keys()) & (spec2 in orth_positions_pre_WGH.keys())): 
        spec1_column = orth_positions_pre_WGH[spec1]
        spec2_column = orth_positions_pre_WGH[spec2]

        with open(fname) as f:
            orth_lookup = {}
            for line in f:
                linesp = line.split()
                spec1_gene = linesp[spec1_column]
                if spec1_gene!= '---':
                    spec2_gene = linesp[spec2_column]
                    if spec2_gene!='---': 
                        orth_lookup[spec1_gene] = spec2_gene
                    else: 
                        orth_lookup[spec1_gene] = 'NONE'




    #Case 3: spec1 pre WGH, spec2 post WGH
    #In this case there will often be more than one ortholog per pre WGH gene

    if ((spec1 in orth_positions_pre_WGH.keys()) & (spec2 in orth_positions_post_WGH.keys())): 
        spec1_column = orth_positions_pre_WGH[spec1]
        spec2_columns = orth_positions_post_WGH[spec2]

        with open(fname) as f:
            orth_lookup = {}
            for line in f:
                linesp = line.split()
                spec1_gene = linesp[spec1_column]
                if spec1_gene!= '---':
                    spec2_genes = []
                    for spec2_column in spec2_columns: 
                        spec2_gene = linesp[spec2_column]
                        if spec2_gene!='---':
                            spec2_genes.append(spec2_gene)
                    if len(spec2_genes)==0:
                        orth_lookup[spec1_gene] = 'NONE'
                    else: 
                        orth_lookup[spec1_gene] = spec2_genes


    #Case 4: spec1 post WGH, spec2 post WGH
    #use tracks in pillars file to match up WGH genes.  If there are two paralogs, this will match tracks. 
    #If there is one paralog, and there are two orthologs, it will list both.  
    #If there is only one ortholog for one gene, it will list the other one regardless of track.  

    if ((spec1 in orth_positions_post_WGH.keys()) & (spec2 in orth_positions_post_WGH.keys())): 
        spec1_columns = orth_positions_post_WGH[spec1]
        spec2_columns = orth_positions_post_WGH[spec2]

        with open(fname) as f:
            orth_lookup = {}
            for line in f:
                linesp = line.split()
                spec1_genes = []
                for spec1_column in spec1_columns: 
                    spec1_gene = linesp[spec1_column]
                    if spec1_gene!= '---':
                        spec1_genes.append(spec1_gene)

                if len(spec1_genes)==2: 
                    for jj,spec1_gene in enumerate(spec1_genes): 
                        spec2_column = spec2_columns[jj]
                        spec2_gene = linesp[spec2_column]
                        if spec2_gene!='---':
                            orth_lookup[spec1_gene] = [spec2_gene]
                        else: 
                            orth_lookup[spec1_gene] = ['NONE']

                elif len(spec1_genes)==1: 
                    spec1_gene = spec1_genes[0]
                    spec2_genes = []
                    for spec2_column in spec2_columns: 
                        spec2_gene = linesp[spec2_column]
                        if spec2_gene!='---':
                            spec2_genes.append(spec2_gene)
                    if len(spec2_genes)==0:
                        orth_lookup[spec1_gene] = 'NONE'
                    else: 
                        orth_lookup[spec1_gene] = spec2_genes





    orth_lookup_outputfname = os.path.normpath(data_processing_dir + 'ortholog_files_YGOB/' + spec1 + "-" + spec2 + "-goi-orthologs.txt"  )

    print(orth_lookup_outputfname)
    print_ortholog_file(orth_lookup_outputfname, orth_lookup)

    return orth_lookup 

def make_orth_file_YGOB_regev(spec1,spec2):
    #makes an ortholog file that goes from 
    #spec1_spec2-orthologs.txt (from YGOB dir) -> spec2_YGOB-spec2_regev-orthologs.txt (from regev dir)

    #Load Spec1-Spec2 from YGOB
    orth_dir_YGOB =  data_processing_dir + "ortholog_files_YGOB" + os.sep
    spec1_spec2_YGOB = read_orth_lookup_table(spec1,spec2, orth_dir_YGOB)

    #Load spec2_YGOB-spec2_Regev
    orth_dir_regev = data_processing_dir + "ortholog_files_regev" + os.sep
    spec2_YGOB_regev = read_orth_lookup_table(spec2 + '_YGOB',spec2 + '_regev', orth_dir_regev)

    spec1_spec2_regev = {}
    for spec1_gene_YGOB, spec2_genes_YGOB in spec1_spec2_YGOB.items(): 
        if len(spec2_genes_YGOB)==1: 
            spec2_gene_YGOB= spec2_genes_YGOB[0]
            if spec2_gene_YGOB == 'NONE':
                spec1_spec2_regev[spec1_gene_YGOB] = ['NONE']
            else: 
                spec1_spec2_regev[spec1_gene_YGOB] = spec2_YGOB_regev[spec2_gene_YGOB]
        else:  
            spec2_genes_regev = []
            for spec2_gene_YGOB in spec2_genes_YGOB: 
                spec2_genes_regev_jj = spec2_YGOB_regev[spec2_gene_YGOB]
                if not(spec2_genes_regev_jj[0]=='NONE'):
                    spec2_genes_regev.append(spec2_genes_regev_jj)
            if len(spec2_genes_regev)==0:
                spec1_spec2_regev[spec1_gene_YGOB] = ['NONE']
            else: 
                spec1_spec2_regev[spec1_gene_YGOB] = list(chain.from_iterable(spec2_genes_regev))
    
    orth_lookup_outputfname = os.path.normpath(data_processing_dir + 'ortholog_files_regev/' + spec1 + "-" + spec2 + "-orthologs.txt"  )
    
    print('saving ' + orth_lookup_outputfname)
    print_ortholog_file(orth_lookup_outputfname, spec1_spec2_regev)
    
    return spec1_spec2_regev

def make_orth_file_regev_YGOB(spec1, spec2):
    #makes an ortholog file that goes from 
    #spec1_regev-spec1_YGOB-orthologs.txt (from regev dir) -> spec1_spec2-orthologs.txt (from YGOB dir)

    #Load spec1_regev-spec1_YGOB from regev dir
    orth_dir_regev = data_processing_dir + "ortholog_files_regev" + os.sep
    spec1_regev_YGOB = read_orth_lookup_table(spec1 + '_regev',spec1 + '_YGOB', orth_dir_regev)

    #Load Spec1-Spec2 from YGOB
    orth_dir_YGOB =  data_processing_dir + "ortholog_files_YGOB" + os.sep
    spec1_spec2_YGOB = read_orth_lookup_table(spec1,spec2, orth_dir_YGOB)

    spec1_spec2_regev_YGOB = {}
    for spec1_gene_regev, spec1_genes_YGOB in spec1_regev_YGOB.items(): 
        if len(spec1_genes_YGOB)==1: 
            spec1_gene_YGOB = spec1_genes_YGOB[0]
            if spec1_gene_YGOB == 'NONE':
                spec1_spec2_regev_YGOB[spec1_gene_regev] = ['NONE']
            else: 
                spec1_spec2_regev_YGOB[spec1_gene_regev] = spec1_spec2_YGOB[spec1_gene_YGOB]
        else:  
            spec2_genes_YGOB = []
            for spec1_gene_YGOB in spec1_genes_YGOB: 
                spec2_genes_YGOB_jj = spec1_spec2_YGOB[spec1_gene_YGOB]
                if not(spec2_genes_YGOB_jj[0]=='NONE'):
                    spec2_genes_YGOB.append(spec2_genes_YGOB_jj)
            if len(spec2_genes_YGOB)==0:
                spec1_spec2_regev_YGOB[spec1_gene_regev] = ['NONE']
            else: 
                spec1_spec2_regev_YGOB[spec1_gene_regev] = list(chain.from_iterable(spec2_genes_YGOB))

    orth_lookup_outputfname = os.path.normpath(data_processing_dir + 'ortholog_files_regev/' + spec1 + "-" + spec2 + "-orthologs.txt"  )

    print('saving ' + orth_lookup_outputfname)
    print_ortholog_file(orth_lookup_outputfname, spec1_spec2_regev_YGOB)

    return spec1_spec2_regev_YGOB

def pairwise_hits_to_orth_dict(pairwise_score, high_thresh, low_thresh, diff_thresh):
    #Given a dictionary from a species genename list to a list of options with scores for hits, 
    #makes a dictionary 1:many from one to the other using a high threshold, low threshold and difference in threshold to drop
    #for throwing out possibilities. 
    
    YGOB_regev = {}
    for genename_YGOB, options in pairwise_score.items(): 
        genename_max = options.idxmax()
        max_val = options[genename_max]
        if max_val < low_thresh: 
            print('No Match ' + genename_YGOB)
            YGOB_regev[genename_YGOB] = ['NONE']
        else: 
            assert max_val > high_thresh , 'Max val not greater than high threshold'
            options_rest = options.drop(genename_max)
            genename_nextmax = options_rest.idxmax()
            nextmax_val = options_rest[genename_nextmax]
            if nextmax_val > high_thresh: 
                diff = max_val-nextmax_val
                print('More than one match {}, diff = {:.2f}'.format(genename_YGOB, diff))
                genename_max_list = [genename_max]

                #Keep adding regev genenames to the list until they are past the threshold
                while diff<diff_thresh:
                    genename_max_list.append(genename_nextmax)
                    nextmax_val_old = nextmax_val
                    genename_nextmax_old = genename_nextmax
                    if len(options_rest)==1: 
                        print("diff thresh never met " + genename_YGOB)
                        genename_nextmax = options_rest.index[0]
                        genename_max_list.append(genename_nextmax)
                        diff = diff_thresh + 2
                    else: 
                        options_rest = options_rest.drop(genename_nextmax_old)
                        genename_nextmax = options_rest.idxmax()
                        nextmax_val = options_rest[genename_nextmax]
                        diff = nextmax_val_old-nextmax_val
                YGOB_regev[genename_YGOB] = genename_max_list
            else: 
                YGOB_regev[genename_YGOB] = [genename_max]

    #     maxvals.append(maxvals)

    return YGOB_regev

def load_YGOB_annotations(species, species_tab_file):
    #previously base_dir was the middle input
    fname = os.path.normpath(base_dir + species_tab_file)
    
    with open(fname) as f:
        annotation_lookup = {}
        for line in f:
            linesp = line.split('\t')
            gene = linesp[0]
            annotation = linesp[8]
            annotation_lookup[gene] = annotation
    
    return annotation_lookup

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

def read_orth_lookup_table(species1, species2, orth_dir):
    #ensure orth_dir has a separator at the end of it.
    #examples: 
    #YGOB: data_processing_dir + "ortholog_files_YGOB" + os.sep
    #regev: data_processing_dir + "ortholog_files_regev" + os.sep
    #y1000: data_processing_dir + "ortholog_files_y1000" + os.sep
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

def join_ohnologs_and_sort(data_to_add, ohnologs, sort_column):
    #Adds data_to_add, a dataframe indexed on genenames to a dataframe of ohnologs
    #from YGOB_pillars (i.e. from get_WGH_pairs_by_spec) indexed on ancestor, and with column names 
    #genename_gene1 and genename_gene2.  
    #sort_column is a column in data_to_add and output paralogs are sorted into low and high sets based on that
    #Example Sort column: 'log2FoldChange
    # 
    # For S. Cer, we have the ohnologs saved in a file, and need to reindex and rename prior to using it.  
    # 
    # ohnologs = pd.read_csv(data_processing_dir + os.path.normpath("ortholog_files_YGOB/ohnologs.csv"), index_col=0)
    # ohnologs.set_index('Ancestor', inplace=True)
    # ohnologs.rename(columns={"Gene 1": "genename_gene1",
    #                          "Gene 2": "genename_gene2"}, inplace=True)
    # 
    # After running for S.Cer, to add SC_common_name, 
    #
    # for level in ['low','high']: 
    #     ohnologs_expression_sorted['SC_common_name_' + level] = yeast_esr_exp.SC_common_name_lookup(ohnologs_expression_sorted['genename_' + level])
    # ohnologs_expression_sorted.drop(columns=['Gene Name 1', 'Gene Name 2'], inplace=True)
    
    # # Rename Gene Name 1 and Gene Name 2 columns
    # ohnologs_expression.rename(columns = {'Gene 1' : 'genename_gene1',
    #                                       'Gene 2' : 'genename_gene2'}, inplace = True)
    
    ohnologs_expression_gene1 = pd.merge(ohnologs, data_to_add, how = 'inner', left_on = 'genename_gene1', right_index = True)
    ohnologs_expression = pd.merge(ohnologs_expression_gene1, data_to_add, how = 'inner', left_on = 'genename_gene2', right_index = True, suffixes = ['_gene1', '_gene2'])

    ## This seems important if you are not merging on the index
    ## Drop Gene 1 and Gene 2 columns
    #ohnologs_expression.drop(['Gene 1', 'Gene 2'], axis = 1, inplace = True)
    
    new_columns = {}
    for level in ['low','high']:
        new_columns['genename_' + level] = []
        for column_name in data_to_add.columns:
            new_columns[column_name + '_' + level] = []

    for index, row in ohnologs_expression.iterrows():
        #Decide if gene1 or gene2 is low expression
        if row[sort_column + '_gene1']<row[sort_column + '_gene2']:
            low_gene = 'gene1'
            high_gene = 'gene2'
        elif  row[sort_column + '_gene1']>row[sort_column + '_gene2']:
            low_gene = 'gene2'
            high_gene = 'gene1'
        else:
            print('problems with {} and {}'.format(row['genename_gene1'],row['genename_gene2']))

        level_dict = {'low':low_gene, 'high': high_gene}

        for level, level_gene in level_dict.items():
            new_columns['genename_' + level].append(row['genename_'+level_gene])
            for column_name in data_to_add.columns:
                new_columns[column_name + '_' + level].append(row[column_name + '_' + level_gene])

    ohnologs_expression_sorted = ohnologs_expression.copy()

    for level in ['low','high']:
        ohnologs_expression_sorted['genename_' + level] = new_columns['genename_' + level]
        for column_name in data_to_add.columns:
            ohnologs_expression_sorted[column_name + '_' + level] = new_columns[column_name + '_' + level] 

    columns_to_drop = []

    for gene in ['gene1','gene2']:
        columns_to_drop.append('genename_' + gene)
        for column_name in data_to_add.columns:
            columns_to_drop.append(column_name + '_' + gene)

    ohnologs_expression_sorted.drop(columns=columns_to_drop, inplace=True)

    return ohnologs_expression_sorted

def join_ohnologs_and_sort_filter(data_to_add, ohnologs, sort_column, filter_column, filter_value, filter_direction):
    #Adds data_to_add, a dataframe indexed on genenames to a dataframe of ohnologs, sorting based on sort column, after filtering based on filter column
    #
    #from YGOB_pillars (i.e. from get_WGH_pairs_by_spec) indexed on ancestor, and with column names 
    #genename_gene1 and genename_gene2.  
    #sort_column is a column in data_to_add and output paralogs are sorted into low and high sets based on that
    #Example Sort column: 'log2FoldChange
    #
    #filter_column is the column to use to filter for significance prior to sorting, and you filter based on filter direction and value. 
    #This is for cases such as RFS1/PST2 where RFS1 has a lower LFC than PST2, its ohnolog, but PST2 does not meet the p-value threshold)
    #
    #In this routine each column is first filtered and 0 assigned if the filter is not passed.  
    #
    # For S. Cer, we have the ohnologs saved in a file, and need to rindex and rename prior to using it.  
    # 
    # ohnologs = pd.read_csv(data_processing_dir + os.path.normpath("ortholog_files_YGOB/ohnologs.csv"), index_col=0)
    # ohnologs.set_index('Ancestor', inplace=True)
    # ohnologs.rename(columns={"Gene 1": "genename_gene1",
    #                          "Gene 2": "genename_gene2"}, inplace=True)
    # 
    # After running for S.Cer, to add SC_common_name, 
    #
    # for level in ['low','high']: 
    #     ohnologs_expression_sorted['SC_common_name_' + level] = yeast_esr_exp.SC_common_name_lookup(ohnologs_expression_sorted['genename_' + level])
    # ohnologs_expression_sorted.drop(columns=['Gene Name 1', 'Gene Name 2'], inplace=True)

    # # Rename Gene Name 1 and Gene Name 2 columns
    # ohnologs_expression.rename(columns = {'Gene 1' : 'genename_gene1',
    #                                       'Gene 2' : 'genename_gene2'}, inplace = True)

    ohnologs_expression_gene1 = pd.merge(ohnologs, data_to_add, how = 'inner', left_on = 'genename_gene1', right_index = True)
    ohnologs_expression = pd.merge(ohnologs_expression_gene1, data_to_add, how = 'inner', left_on = 'genename_gene2', right_index = True, suffixes = ['_gene1', '_gene2'])

    
    for gene_label in ['gene1','gene2']:
    
        ohnologs_expression[sort_column + '_' + gene_label + '_filt'] = 0
        if filter_direction=='less':
            pass_inds = ohnologs_expression[filter_column + '_' + gene_label]< filter_value
        elif filter_direction == 'more':
            pass_inds = ohnologs_expression[filter_column + '_' + gene_label]> filter_value
        else:
            raise ValueError('filter direction should be less or more')

        ohnologs_expression.loc[pass_inds, sort_column + '_' + gene_label + '_filt'] = ohnologs_expression.loc[pass_inds, sort_column + '_' + gene_label]

    
    #Identify when gene1 is higher: 
    gene1_higher = ohnologs_expression[sort_column + '_gene1_filt']>ohnologs_expression[sort_column + '_gene2_filt']

    #Identify when gene2 is higher: 
    gene2_higher = ohnologs_expression[sort_column + '_gene1_filt']<ohnologs_expression[sort_column + '_gene2_filt']

    #For genes in which at least one passes the filter, assigne low and high genes

    ohnologs_expression['low_gene'] = 'blank'
    ohnologs_expression['high_gene'] = 'blank'
    ohnologs_expression.loc[gene1_higher,'low_gene']='gene2'
    ohnologs_expression.loc[gene1_higher,'high_gene']='gene1'
    ohnologs_expression.loc[gene2_higher,'low_gene']='gene1'
    ohnologs_expression.loc[gene2_higher,'high_gene']='gene2'

    #When neither pass the filter, just use LFC
    #identify when neither pass the filter: 
    low_sig = ohnologs_expression[sort_column + '_gene1_filt']==ohnologs_expression[sort_column + '_gene2_filt']
    gene1_higher_low_sig = low_sig & (ohnologs_expression[sort_column + '_gene1']>ohnologs_expression[sort_column + '_gene2']) 
    gene2_higher_low_sig = low_sig & (ohnologs_expression[sort_column + '_gene1']<ohnologs_expression[sort_column + '_gene2']) 

    ohnologs_expression.loc[gene1_higher_low_sig,'low_gene']='gene2'
    ohnologs_expression.loc[gene1_higher_low_sig,'high_gene']='gene1'
    ohnologs_expression.loc[gene2_higher_low_sig,'low_gene']='gene1'
    ohnologs_expression.loc[gene2_higher_low_sig,'high_gene']='gene2'

    #Swap names based on which is higher

    #add new columns for low/high data
    for level in ['low','high']:
        for column_name in list(data_to_add.columns) + ['genename']:
            ohnologs_expression[column_name + '_' + level]='blank'

    #add data to new columns based on which gene is low and which is high

    opposites = {'gene1':'gene2','gene2':'gene1'}

    for low_gene in ['gene1','gene2']: 
        high_gene = opposites[low_gene]
        level = 'low'
        for column_name in list(data_to_add.columns) + ['genename']: 
            column_name_low_gene = column_name + '_' + low_gene
            column_name_low = column_name + '_low'
            column_name_high_gene = column_name + '_' + high_gene
            column_name_high = column_name + '_high'
            ohnologs_expression.loc[(ohnologs_expression['low_gene']==low_gene),column_name_low] =  ohnologs_expression.loc[(ohnologs_expression['low_gene']==low_gene),column_name_low_gene]
            ohnologs_expression.loc[(ohnologs_expression['low_gene']==low_gene),column_name_high] =  ohnologs_expression.loc[(ohnologs_expression['low_gene']==low_gene),column_name_high_gene]

    #delete all columns with gene1 and gene2

    columns_to_drop = ['log2FoldChange_gene1_filt', 'log2FoldChange_gene2_filt', 'low_gene', 'high_gene']

    for gene_label in ['gene1','gene2']:
        columns_to_drop.append('genename_' + gene_label)
        for column_name in list(data_to_add.columns) + ['genename']:
            columns_to_drop.append(column_name + '_' + gene_label)

    ohnologs_expression_sorted = ohnologs_expression.drop(columns=columns_to_drop)

        
    return ohnologs_expression_sorted
        
def print_ortholog_file(orth_lookup_outputfname, orth_lookup):
    with open(orth_lookup_outputfname, 'w') as fw:
        for spec1_gene, spec2_genes in orth_lookup.items():
            if isinstance(spec2_genes,str): 
                #if this is the case spec2_genes is just one gene
                line = spec1_gene + '\t' + spec2_genes + '\n'
                fw.write(line)
            if isinstance(spec2_genes,list): 
                line = spec1_gene + '\t' + '\t'.join(spec2_genes) + '\n'
                fw.write(line)
    
    return

## Grouping genes

def build_pka_act_rep_target_sets():
    #Using genes that were identified as activated and repressed in S. cer and K. lac, builds sets that combine those 
    #For functional enrichment with go terms and enrichmenet calculation for ohnologs.  
    #Also builds background sets.  Background for SC is all genes that have DEseq data, background set for KL is all genes 
    #that have orthologs in SC, and have data in both SC and KL.  
    #
    #Uses target setfiles created in these scripts:
    #20181031_klscpka_r1g1_m24_SC_analysis and
    #20181031_klscpka_r1g1_m24_KL_analysis
    #
    #Scer-Klac ortholog file saved in ortholog_files_YGOB 
    #
    #Deseq input from  
    #20181017_klscpka_m24_r1g1_deseq2.rmd

    target_set_files = {'SC': {'pkainh_act':'20200608_pkainh_act_SC.csv',
                               'pkainh_rep':'20200608_pkainh_rep_SC.csv'},
                        'KL': {'pkainh_act':'20200608_pkainh_act_KL.csv',
                               'pkainh_rep':'20200608_pkainh_rep_KL.csv'}
                        }

    target_sets = {}

    for spec,spec_sets in target_set_files.items():
        target_sets_spec = {}

        for set_name, fname in spec_sets.items():
            target_set_df = pd.read_csv(os.path.normpath(data_processing_dir + 'kl_sc_PKA_as_m24_r1g1_20181017/' + fname ), index_col=0)
            target_set = set(target_set_df.index)
            target_sets_spec[set_name] = target_set

        target_sets[spec] = target_sets_spec

    # kl_orthologs = pd.read_pickle(data_processing_dir + "ortholog_files_YGOB/kl_orthologs.pkl")
    kl_orthologs = read_orth_lookup_table('Klac', 'Scer', data_processing_dir + os.sep + 'ortholog_files_YGOB' + os.sep)

    #Make different gene sets to look up orthologs (using SC_genenames as a key).   For SC will have all the genes, for KL will only have
    #genes with SC orthologs. 

    gene_set_names = ['kl_only_act', 'kl_only_rep','sc_only_act','sc_only_rep','klsc_act','klsc_rep','sc_act','sc_rep','kl_act','kl_rep']

    #For kl_act we include the names of all SC orthologs of genes activated   
    #If there is a kl_act gene that has two orthologs in S.cer, we include both of them.  

    kl_sets = {}
    for setname in ['pkainh_act','pkainh_rep']:
        kl_set = []
        for kl_gene in list(target_sets['KL'][setname]): 
            if not(kl_gene in kl_orthologs.keys()):
                print(kl_gene + ' from ' + setname + ' not in kl ortholog index')
            else: 
                sc_orth = kl_orthologs[kl_gene]
                if sc_orth != ['NONE']:
                    kl_set = kl_set + sc_orth
        kl_sets[setname] = set(kl_set)

    kl_act = kl_sets['pkainh_act']
    kl_rep = kl_sets['pkainh_rep']


    #kl_act = set(kl_orthologs.loc[kl_orthologs['kl_genename'].isin(target_sets['KL']['pkainh_act'],'sc_genename'])
    #kl_rep = set(kl_orthologs.loc[kl_orthologs['kl_genename'].isin(target_sets['KL']['pkainh_rep']),'sc_genename'])


    sc_act = target_sets['SC']['pkainh_act']
    sc_rep = target_sets['SC']['pkainh_rep']

    #Load DEseq data for SCer PKA AS -/+ NMPP1 
    #pkainh_deseq_SC = pd.read_csv(os.path.normpath(data_processing_dir + '\\kl_sc_PKA_as_m24_r1g1_20181017\\20181017_deseq_SC_AS_WT_nmpp1.csv'), index_col=0)
    pkainh_deseq_SC = pd.read_csv(os.path.normpath(data_processing_dir + '\\kl_sc_PKA_as_m24_r1g1_20181017\\20200603_deseq_SC_AS_WT_nmpp1.csv'), index_col=0)

    pkainh_deseq_SC['SC_common_name'] = SC_common_name_lookup(pkainh_deseq_SC.index)

    #Load DEseq data for KLac PKA AS -/+ NMPP1
    #Load DEseq data for SCer PKA AS -/+ NMPP1 
    #pkainh_deseq_KL = pd.read_csv(os.path.normpath(data_processing_dir + '\\kl_sc_PKA_as_m24_r1g1_20181017\\20181017_deseq_KL_AS_WT_nmpp1.csv'), index_col=0)
    pkainh_deseq_KL = pd.read_csv(os.path.normpath(data_processing_dir + '\\kl_sc_PKA_as_m24_r1g1_20181017\\20200603_deseq_KL_AS_WT_nmpp1.csv'), index_col=0)
    kl_genes_w_data = set(kl_genename_convert_list(list(pkainh_deseq_KL.index)))

    gene_sets = { 'kl_only_act' : kl_act - sc_act, 
                  'kl_only_rep' : kl_rep - sc_rep, 
                  'sc_only_act' : sc_act - kl_act, 
                  'sc_only_rep' : sc_rep - kl_rep, 
                  'klsc_act' : kl_act & sc_act, 
                  'klsc_rep' : kl_rep & sc_rep,
                  'sc_act':sc_act,
                  'sc_rep':sc_rep,
                  'kl_act':kl_act,
                  'kl_rep':kl_rep,
                  'sc_no_change':set(pkainh_deseq_SC.index)-(sc_act | sc_rep),
                  'kl_no_change':kl_genes_w_data-(kl_act | kl_rep) 
                  }


    #Prepare background set


    background_map = {'kl_only_act' : 'KL_orth',
                     'kl_only_rep' : 'KL_orth', 
                     'sc_only_act' : 'SC',
                     'sc_only_rep': 'SC',
                     'sc_no_change': 'SC',
                     'klsc_act': 'KL_orth',
                     'klsc_rep': 'KL_orth',
                     'sc_act':'SC',
                     'sc_rep':'SC',
                     'kl_act':'KL',
                     'kl_rep':'KL', 
                     'sc_no_change':'SC',
                     'kl_no_change':'KL'
                     }


    #Only keeps KL gene if it has data in the KL DEpka set
    all_kl_orths = []
    #5188 genes with KL orthologs
    #of those, only 4953 have data in KL
    #of those, only 4802 have data for the SC ortholog as well.  
    for kl_genename, kl_orth in kl_orthologs.items():
        if kl_genename in kl_genes_w_data:
            all_kl_orths = all_kl_orths + kl_orth
    all_kl_orths_set = set(all_kl_orths)-{'NONE'}

    background_genes = {'SC' : set(pkainh_deseq_SC.index), 
                        'KL': kl_genes_w_data,
                        'KL_orth' : all_kl_orths_set & set(pkainh_deseq_SC.index)
                       }
    return gene_sets, background_map, background_genes

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

def get_gis1_rph1_sets(act_threshold, inh_threshold, kl_sc_PKA_data_subset, condition = 'log'): 
    #Use cutoff to assign activated, inhibited, and repressed for gis1, rph1 and gis1rph1 for each condition in the westholm dataset
    #input is
    #   threshold levels 
    #   kl_sc_PKA_data_subset - must have the columns for the data: ['log_g-wt','log_r-wt', 'log_gr-wt']
    # condition can be 'log', 'PDS' or '3d' for 3days stationary phase
    
    rg_columns = [condition + '_g-wt',condition + '_r-wt', condition + '_gr-wt']
    for column in rg_columns: 
        kl_sc_PKA_data_subset[column+'_label'] = threshold_group_series(kl_sc_PKA_data_subset[column], act_threshold, inh_threshold )

    strains = ['g','r','gr']
    gr_sets = {}
    exp_profiles = tuple(product(['up','flat','down'],['up','flat','down'],['up','flat','down']))
    for g_label, r_label, gr_label in exp_profiles:
        gr_sets['g-'+g_label + '_r-'+r_label + '_gr-'+gr_label] = list(kl_sc_PKA_data_subset[(kl_sc_PKA_data_subset[condition + '_g-wt_label']==g_label) &
                                                                      (kl_sc_PKA_data_subset[condition + '_r-wt_label']==r_label) &
                                                                      (kl_sc_PKA_data_subset[condition + '_gr-wt_label']==gr_label)]['sc_genename'])


    return kl_sc_PKA_data_subset, gr_sets

## Go analysis
def load_goslim_data(GO_aspect, go_slim_fname = 'go_slim_mapping_20181204.tab'):
    #The three GO_aspect values are: 
    #C = cellular_component
    #F = molecular_function
    #P = biological_process
    go_slims = pd.read_table(data_processing_dir + os.path.normpath('go_terms/' + go_slim_fname),header = None)
    go_slims.columns = ['sc_genename','SC_common_name','sgd_ID','GO_aspect','GO_term','GO_term_ID','feature_type']
    
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
        subset_genes_in_goterm =   list(set(gene_set_list) & set(term_genes))
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

## Promoter analysis

class PromHits(dict):
    def __init__(self, prom_counts_subset, motif_dict, prefix = ''):
        #7/29/19 Changed format of motif dict so that it included length of the promoter to check as well as the sequence
        if prefix != '': 
            prefix = prefix + '_'
        self.motif_dict = motif_dict   #records the motif_dict input as metadata
        self['total']= len(prom_counts_subset)
        self.prom_counts = prom_counts_subset
        self['hits'] = {}
        self['hits_vec'] = {}
        self['pct'] = {}
        self['avg'] = {}
        for motif_name in motif_dict.keys():
            #hits_motif = sum(prom_counts_subset[prefix + motif_name + '_count']>0)
            L_prom = motif_dict[motif_name][1]
            hits_vec = []
            
            for hit_list in prom_counts_subset[prefix + motif_name + '_full_features']:
                N_hits = 0
                for (loc, strand, context) in hit_list: 
                    if loc<L_prom:
                        N_hits = N_hits + 1
                hits_vec.append(N_hits)
            
            hits_motif = sum(np.array(hits_vec)>0)          
            self['hits'][motif_name] = hits_motif
            
            hits_vec_series = pd.Series(index = prom_counts_subset.index, data= hits_vec)
            self['hits_vec'][motif_name]=hits_vec_series           
            
            if self['total']==0:
                self['pct'][motif_name] = np.nan
            else: 
                self['pct'][motif_name] = hits_motif/self['total']
            self['avg'][motif_name] = prom_counts_subset[prefix + motif_name + '_count'].sum()/self['total']
    
    def STRE_TATA_combined(self,ranges): 
        #Quantify combined STRE and TATA presence
        
        hits = []

        for row in self.prom_counts.iterrows():
            motif_conditions = {}
            for motif in ['STRE','TATA']:
                motif_conditions[motif] = False
                if row[1][motif + '_count']>0:
                    for (loc,strand,motif_found) in row[1][motif+'_full_features']:
                        if loc<ranges[motif]:
                            motif_conditions[motif]=True

            if motif_conditions['STRE'] & motif_conditions['TATA']:
                hits.append(row[0])
                
        self['hits']['STRE_TATA'] = len(hits)
        self['STRE_TATA_hits'] = hits
        if self['total']==0:
            self['pct']['STRE_TATA'] = np.nan
        else: 
            self['pct']['STRE_TATA'] = len(hits)/self['total']
                   
class PromComparison(dict):
    #Object that gathers data on motif hits between different sets of promoters.
    def __init__(self, promsets):
        self['promsets'] = promsets
    
    def generate_pval(self, promhits1, promhits2, set2type, motif_dict):
        #set2type can be 'all' or 'another_set'
        self.PromHits1 = promhits1
        self.PromHits2 = promhits2
        self.set2type = set2type
        N1_total = promhits1['total']
        self['pval'] = {}
        for motif_name in motif_dict.keys():
            N1_hits = promhits1['hits'][motif_name]
            if set2type=='all': 
                N_all_hits = promhits2['hits'][motif_name]
                N_all_total = promhits2['total'] 
                N2_hits = N_all_hits-N1_hits
                N2_total = N_all_total - N1_total
            elif set2type == 'another_set': 
                N2_hits = promhits2['hits'][motif_name]
                N2_total = promhits2['total'] 

            oddsratio, pvalue = stats.fisher_exact([[N1_hits, N2_hits],
                                                [N1_total, N2_total]], 
                                                  alternative = 'two-sided')
            self['pval'][motif_name] = pvalue
            
    #Would be nice to have a test to see if the average number of STREs is 
    #different - would need to look at distributions first.  Some discrete 
    #distribution with lots of zeros.  
    

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


def ygob_promoter_extract(spec, L_prom): 
    #extracts promoters from a YGOB genome.tab file and deposits them 
    #in a file.  Also deposits short promoters in a file.  
    #only works on gpucluster right now because that's where the genomes are
    
    ygob_genome_dir = "/home/heineike/genomes/YGOB/"

    ygob_fname_spec = ygob_genome_dir + spec + '_genome.tab'

    features = pd.read_table(ygob_fname_spec, index_col=0, header=None)

    features.columns = ['strand', 'start', 'end', 'ygob on/off', 'chrm', 'short_name', 'coords', 'notes']
    #  - NAME - unambiguous name used to identify the gene
    #  - ORIENTATION - 0 = Crick Strand, 1 = Watson Strand
    #  - START COORDINATE - 5' start of the co-ordinate range of the feature
    #  - STOP COORDINATE - 3' end of the co-ordinate range of the feature.  I call it 'end' to match typical GTFs
    #  - ON/OFF - whether feature is displayed or not in YGOB
    #  - CHROMOSOME/CONTIG/SCAFFOLD NUMBER - identifying number of source sequence
    #  - SHORT NAME - the shorter name that will appear in the gene box on screen in YGOB
    #  - COORDINATES - complete gene co-ordinates with intron/exon annotation and complement tag if appropriate
    #  - NOTES - Gen

    genome_fname = ygob_genome_dir + spec + '_sequence.fsa'

    #given a features file, extract the promoters for all the genes: 
    #would be easy to extend to extracting for a subset of genes, but might be better to do that from the 
    #all promoter file anyway

    #genes = ['NCAS0A00120', 'NCAS0A00150', 'NCAS0J02180']

    #probably should separate out genes by chromosome

    #get list of chromosomes
    chrm_list = list(SeqIO.index(genome_fname, "fasta"))

    #make dict of chromosome name in genome file to chromosome name in features table 
    chrm_dict = {chrm: int(chrm.split('_')[-1]) for chrm in chrm_list}

    strand_dict = {0:'-', 1:'+'}

    seq_records = SeqIO.parse(genome_fname, "fasta")

    promoter_dir = data_processing_dir + os.path.normpath('promoter_phylogeny/promoter_sets/' + spec)
    #os.mkdir(promoter_dir)
    promoters_fname = os.path.normpath(promoter_dir + '/all_promoters_' + str(L_prom) + '.fasta')

    with open(promoters_fname, 'w') as f: 

        short_promoters = {}
        for chrm_seq in seq_records:
            chrm = chrm_seq.id
            chrm_features = chrm_dict[chrm]

            features_chrm = features[features['chrm']==chrm_features]
            #genes_chrm = list(set(genes) & set(features_chrm.index))

            #for each gene, extract promoter sequence
            for feature in features_chrm.iterrows():
                gene = feature[0]
                strand = strand_dict[feature[1].strand]
                start = feature[1].start
                end = feature[1].end
                #Adjust coordinates to get L_prom "promoter" sequences
                if strand == '-': 
                    prom_end = end
                    prom_start = prom_end + L_prom   #should do min of this and the total length of the scaffold, 
                elif strand == '+': 
                    prom_end = start - 1
                    prom_start = max(0,prom_end - L_prom)


                L_scaffold = len(chrm_seq)

                if strand == '-': 
                    if prom_start > L_scaffold: 
                        print('promoter region extends past the scaffold, spec = ' + spec + 'Gene = ' + gene + ', L_prom = ' + str(L_prom))
                        prom_start = L_scaffold
                    if prom_end > L_scaffold: 
                        print('scaffold ends at the end of the gene, spec = ' + spec + ' Gene = ' + gene)
                        prom_end = L_scaffold

                    promoter = chrm_seq.seq[prom_end:prom_start].reverse_complement()
                elif strand == '+': 
                    promoter = chrm_seq.seq[prom_start:prom_end]

                if abs(prom_end-prom_start)<L_prom:
                    short_promoters[gene] = abs(prom_end-prom_start)

                #do not add promoter if it has L=0 
                if abs(prom_end-prom_start)>0: 
                    f.write('>' + gene + ' scaffold=' + chrm + " strand=" + strand + " start=" + str(prom_start) + ' end=' + str(prom_end) +  ' L=' + str(abs(prom_end-prom_start)) + '\n')
                    f.write(str(promoter.upper()) + '\n')  #I wonder why some of the bases were in lower case

    
    print('All promoters for ' + spec + ' saved in ' + promoters_fname)
    short_promoters_fname = os.path.normpath(promoter_dir + '/short_promoters.fasta')
    with open(short_promoters_fname,'w') as f: 
        for gene, L_short_prom in short_promoters.items(): 
            f.write(gene + '\t' + str(L_short_prom) + '\n')     
    
    print('short promoters saved in ' + short_promoters_fname)

    return 


def ygob_AA_extract(spec):
    #Makes fasta file of AA sequences for YGOB species
    print(spec)
    
    ygob_genome_dir = "/home/heineike/genomes/YGOB/"

    ygob_fname_spec = ygob_genome_dir + spec + '_genome.tab'

    features = pd.read_table(ygob_fname_spec, index_col=0, header=None)

    features.columns = ['strand', 'start', 'end', 'ygob on/off', 'chrm', 'short_name', 'coords', 'notes']
    #  - NAME - unambiguous name used to identify the gene
    #  - ORIENTATION - 0 = Crick Strand, 1 = Watson Strand
    #  - START COORDINATE - 5' start of the co-ordinate range of the feature
    #  - STOP COORDINATE - 3' end of the co-ordinate range of the feature.  I call it 'end' to match typical GTFs
    #  - ON/OFF - whether feature is displayed or not in YGOB
    #  - CHROMOSOME/CONTIG/SCAFFOLD NUMBER - identifying number of source sequence
    #  - SHORT NAME - the shorter name that will appear in the gene box on screen in YGOB
    #  - COORDINATES - complete gene co-ordinates with intron/exon annotation and complement tag if appropriate
    #  - NOTES - Gen

    genome_fname = ygob_genome_dir + spec + '_sequence.fsa'

    #get list of chromosomes
    chrm_list = list(SeqIO.index(genome_fname, "fasta"))

    #make dict of chromosome name in genome file to chromosome name in features table 
    chrm_dict = {chrm: int(chrm.split('_')[-1]) for chrm in chrm_list}

    strand_dict = {0:'-', 1:'+'}

    seq_records = SeqIO.parse(genome_fname, "fasta")

    AA_dir = data_processing_dir + os.path.normpath('ortholog_files_YGOB/AA_lists/')
    # #os.mkdir(promoter_dir)
    AA_fname = AA_dir + os.sep + os.path.normpath(spec + '_AA_list.fasta')

    with open(AA_fname, 'w') as f: 

        for chrm_seq in seq_records:
            chrm = chrm_seq.id
            print('scaffold: ' + chrm)
            chrm_features = chrm_dict[chrm]

            features_chrm = features[features['chrm']==chrm_features]
            
            #genes_chrm = list(set(genes) & set(features_chrm.index))

            #for each gene, extract AA sequence
            for feature in features_chrm.iterrows():
                gene = feature[0]
                strand = strand_dict[feature[1].strand]
                start = feature[1].start
                end = feature[1].end

                if strand == '-': 
                    CDS = chrm_seq.seq[start-1:end].reverse_complement()

                elif strand == '+': 
                    CDS = chrm_seq.seq[start-1:end]

                AA_seq = str(CDS.translate())

                f.write('>' + gene + ' scaffold=' + chrm + " strand=" + strand + " start=" + str(start) + ' end=' + str(end) + '\n')
                f.write(AA_seq + '\n')  #I wonder why some of the bases were in lower case


    print('All AA seqs for ' + spec + ' saved in ' + AA_fname)

    return

   
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

def build_motif_dict(fname):
    #removed default - it should be:  = data_processing_dir + os.path.normpath('motifs/JASPAR_CORE_2016_fungi.meme')
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
                     ame_evalue_threshold = '10' ):
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

    motif_db = os.path.normpath(data_processing_dir + '/motifs/'+ motif['fname'])
    target_sequences = os.path.normpath(promoter_dir[spec] + 'promoter_sets/' + fname_prefixes['target'] + '_promoters.fasta')
    control_sequences = os.path.normpath(promoter_dir[spec] + 'promoter_sets/' + fname_prefixes['control'] + '_promoters.fasta')
    output_dir = os.path.normpath(promoter_dir[spec] + 'ame_output')
    file_prefix = target_fname_prefix + '_vs_' + control_fname_prefix + '_motif_' + motif['name'] + '_EVal_' + ame_evalue_threshold + '_' 


    #  Version 5 command
    ame_command = [ "/home/heineike/meme/bin/ame",
                  "--verbose", "2",
                  "--oc", output_dir,
                  "--control", control_sequences,
                  "--scoring", ame_scoring,
                  "--method", ame_method, 
                  "--evalue-report-threshold", ame_evalue_threshold, 
                  target_sequences,
                  motif_db]

    #  Version 4.0.12 command
    #    ame_command = [ "/home/kieran/meme/bin/ame",
    #                  "--verbose", "2",
    #                  "--oc", output_dir,
    #                  "--control", control_sequences,
    #                  "--bgformat", "1", 
    #                  "--scoring", ame_scoring,
    #                  "--method", ame_method, 
    #                  "--pvalue-report-threshold", ame_pvalue_threshold, 
    #                  target_sequences,
    #                  motif_db]
    
    #print(' '.join(ame_command))    

    ame_output = subprocess.run(ame_command,stdout = subprocess.PIPE) 

    print("ame output return code = " + str(ame_output.returncode))

    #change file prefix
    for fname in ["ame.tsv","ame.html"]:
        os.rename(output_dir + os.sep + fname, output_dir + os.sep + file_prefix + fname)
    
    return

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
    #motif_dict = yeast_esr_exp.build_motif_dict(fname)
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

def exact_motif_pattern_len(motif): 
    #Finds the length of patterns such as 'TATA[AT]A[AT]' by counting everything within the brackets as just one letter
    split_left = motif.split('[')

    L = len(split_left[0])

    if len(split_left)>1:
        for section in split_left[1:]:
            L = L + 1 + len(section.split(']')[1])
            
    return L

def exact_promoter_scan(motif, prom_seq, output_format='count', sequence_context=0):
    #Finds a motif in a sequence object
    # 
    #Inputs: 
    #    motif:  A string which is the motif we wish to search (used in re.finditer) ex: 'CCCCT'. 
    #    prom_seq: biopython Seq object.  It will be scanned both forward and backward.  
    #        output_format: if 'count', then just lists the number of hits, if 'full', gives the location, 
    #        direction (relative to the gene), and sequence context (number of bases on either side of the hit)
    #
    #Output: depending on output format is either a number of found motifs, or a list of tuples with one tuple for each hit.  
    #    The tuple contains the location (with refrence to the forward sequence), 
    #    fwd or rev, and the match with surrounding context bases. 
    #
    #Revision notes: this function is for a single sequence - same name used to be for one that searched a list of genes 
    # now that can be done with exact_promoter_scan_genelist or exact_promoter_scan_from_fasta
    #
    # For some reason I had str(prom_seq) rather than str(prom_seq.seq) until 121718
    #
    ## output of full changed to be an empty list if there are no hits.  Before it was None. 
    
    L_motif = exact_motif_pattern_len(motif)

    #fwd search
    prom_seq_fwd = str(prom_seq.seq)
    L_prom = len(prom_seq_fwd)
    motif_sites_fwd = [m.start() for m in re.finditer(motif, prom_seq_fwd)]
    #to find overlapping motifs
    #[m.start() fr m in re.finditer('(?=' + motif + ')', prom_seq_rev)]

    #reverse search
    #rev search
    prom_seq_rev = str(prom_seq.reverse_complement().seq)
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
        
        output = loc_hitseq_gene
        
#         if len(loc_hitseq_gene)==0: 
#             output = None
#         else:
#             output = loc_hitseq_gene
    else: 
        print("choose output_format: 'count' or 'full' ")
        
    return output

def exact_promoter_scan_from_fasta(promoters_fname, motif_dict, output_format = 'counts', sequence_context = 0, seq_key_func = lambda seq : seq.id): 
    #For a fasta file of promoters (with Id and sequence for each line) and a dictionary of exact motifs finds counts, location, and sequence context of found motifs.
    #
    #Input: 
    #    promoters_fname: filename of promoter fasta file
    #    motif_dict: Dictionary of exact motifs to be used e.g. {'STRE': ('CCCCT',L_prom)} 
    #    output_format: 'counts'(default) or 'full'.  If 'full' is selected will also output a list of tuples containing the location of the 
    #               found motif, and the sequence (with sequence context) of the motif 
    #    sequence_context:  number of bases surrounding the found motif to print in full output mode (default is 0)
    #    #deprecated argument, now L_prom is included for each motif   L_prom: desired promoter length - will shorten the length of the promoter if the input from the fasta is longer than this.  If the promoter is shorter than the given length then the full promoter will be used.  The default is None which uses the whole sequence.  
    #
    #Output: Dataframe with each row being the id of the gene and columns for each motif with counts and, if 'full' output mode is selected, 
    #        a full_features column with all features for each hit (a list of tuples), and None if there are no hits

    record_iterator = SeqIO.parse(promoters_fname, "fasta")                                                         

    output_dict = {}
    columns = []

    for motif_name, (motif,L_prom) in motif_dict.items():
        if output_format=='count': 
            columns.append(motif_name + '_count')
        elif output_format=='full': 
            columns = columns + [motif_name + '_' + feature for feature in ['count', 'full_features']]   
        else: 
            print('Invalid output_format : ' + output_format)

    for seq in record_iterator: 
        gene_id = seq_key_func(seq)
        L_seq = len(seq)
        output_row = []
        for motif_name, (motif,L_prom) in motif_dict.items():
            if L_prom == None:
                seq_cropped=seq
            else:
                if L_seq > L_prom: 
                    seq_cropped = seq[(L_seq-L_prom):]
                else: 
                    seq_cropped = seq
            if output_format=='count':
                prom_hits = exact_promoter_scan(motif, seq_cropped, output_format=output_format, sequence_context = sequence_context)
                output_row.append(prom_hits)   #append counts
            elif output_format=='full': 
                prom_hits = exact_promoter_scan(motif, seq_cropped, output_format=output_format, sequence_context = sequence_context)
                #add output in the order of columns
                counts = len(prom_hits)
                output_row.append(counts)     #append count
                output_row.append(prom_hits)  #append full_features

        output_dict[gene_id] = output_row

    output = pd.DataFrame.from_dict(output_dict, orient="index", columns = columns)

    return output

def exact_promoter_scan_genelist(gene_list, motif_dict, promoter_database, output_format = 'counts', sequence_context = 0): 
    #finds nonoverlapping exact matches forward and backward for motifs. 
    #input:  motif_dict,  Dictionary of exact motifs to be used e.g. {'STRE': 'CCCCT'} 
    #        promoter_database 
    #        gene_list from dataframe (genes must be primary key of promoter data structure)
    #        output_format: 'counts'(default) or 'full'.  If 'full' is selected each entry is a list of tuples containing the location of the 
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

def exact_promoter_scan_genelist_dict(gene_list, promoter_dict, motif_dict, output_format = 'counts', sequence_context = 0): 
    #For a dictionary of all promoters and a list of genes, find counts, location and sequence context of desired motifs.  
    #
    #Input: 
    #    gene_list: List of genes to analyze from the promoter dict
    #    promoter_dict: dictionary linking gene name to promoter sequence (a biopython sequence object)
    #    motif_dict: Dictionary of exact motifs to be used e.g. {'STRE': 'CCCCT'} 
    #    output_format: 'counts'(default) or 'full'.  If 'full' is selected will also output a list of tuples containing the location of the 
    #               found motif, and the sequence (with sequence context) of the motif 
    #    sequence_context:  number of bases surrounding the found motif to print in full output mode (default is 0)
    #    L_prom: desired promoter length - will shorten the length of the promoter if the input from the fasta is longer than this.  If the promoter is shorter than the given length then the full promoter will be used.  The default is None which uses the whole sequence.  
    #
    #Output: Dataframe with each row being the id of the gene and columns for each motif with counts and, if 'full' output mode is selected, 
    #        a full_features column with all features for each hit (a list of tuples), and None if there are no hits

    output_dict = {}
    columns = []

    for motif_name, (motif, L_prom) in motif_dict.items():
        if output_format=='count': 
            columns.append(motif_name + '_count')
        elif output_format=='full': 
            columns = columns + [motif_name + '_' + feature for feature in ['count', 'full_features']]   
        else: 
            print('Invalid output_format : ' + output_format)

    for gene in gene_list:
        if gene in promoter_dict.keys():
            seq = promoter_dict[gene]
            L_seq = len(seq)
            output_row = []
            for motif_name, (motif, L_prom) in motif_dict.items():
                #crop seq to appropriate length
                if L_prom == None:
                    seq_cropped=seq
                else:
                    if L_seq > L_prom: 
                        seq_cropped = seq[(L_seq-L_prom):]
                    else: 
                        seq_cropped = seq
                
                if output_format=='count': 
                    prom_hits = exact_promoter_scan(motif, seq_cropped, output_format=output_format, sequence_context = sequence_context)
                    output_row.append(prom_hits)   #append counts
                elif output_format=='full': 
                    prom_hits = exact_promoter_scan(motif, seq_cropped, output_format=output_format, sequence_context = sequence_context)
                    #add output in the order of columns
                    counts = len(prom_hits)
                    output_row.append(counts)     #append count
                    output_row.append(prom_hits)  #append full_features

            output_dict[gene] = output_row
        else:
            print(gene + ' not in promoter_dict')

    output = pd.DataFrame.from_dict(output_dict, orient="index", columns = columns)
    
    return output

def exact_promoter_scan_genelist_orths_dict(gene_list_orths, rule, promoter_dict, motif_dict, output_format = 'counts', sequence_context = 0, L_prom = None): 

    #For a dictionary of all promoters and a list of genes, find counts, location and sequence context of desired motifs.  
    #
    #Input: 
    #    gene_list_orths: List of orthologous genes to analyze from the promoter dict.  This is a list
    #                     of lists because there may be more than one ortholog for a given gene.  
    #    rule:  the rule for choosing how to deconflict when there is more than one ortholog for a given gene. 
    #           right now the only rule is to pick the one with the max STRE.  If there are two with the 
    #           same number of STREs, it just picks the one with the lowest index.
    #    promoter_dict: dictionary linking gene name to promoter sequence (a biopython sequence object)
    #    motif_dict: Dictionary of exact motifs to be used e.g. {'STRE': 'CCCCT'} 
    #    output_format: 'counts'(default) or 'full'.  If 'full' is selected will also output a list of tuples containing the location of the 
    #               found motif, and the sequence (with sequence context) of the motif 
    #    sequence_context:  number of bases surrounding the found motif to print in full output mode (default is 0)
    #    L_prom: desired promoter length - will shorten the length of the promoter if the input from the fasta is longer than this.  If the promoter is shorter than the given length then the full promoter will be used.  The default is None which uses the whole sequence.  
    #
    #Output: 
    #
    #output_df = Dataframe with each row being the id of the gene and columns for each motif with counts and, if 'full' output mode is selected, 
    #        a full_features column with all features for each hit (a list of tuples).
    #
    #chosen_orths = list that has the same order as the input ortholog list that has values of either "NONE", or 
    #               the gene that was chosen to represent the ortholog by the given rule

    output_dict = {}
    columns = []

    for motif_name, motif in motif_dict.items():
        if output_format=='count': 
            columns.append(motif_name + '_count')
        elif output_format=='full': 
            columns = columns + [motif_name + '_' + feature for feature in ['count', 'full_features']]   
        else: 
            print('Invalid output_format : ' + output_format)
    
    if rule=='max_STRE':
        stre_ind = [jj for jj, column in enumerate(columns) if column=='STRE_count'][0]
        
    chosen_orths = []
    for genes in gene_list_orths:
        if len(genes)==0: 
            chosen_orths.append('NONE')
        else:   #1 or more genes in orthologs
            for gene in genes:
                if gene in promoter_dict.keys():
                    seq = promoter_dict[gene]
                    L_seq = len(seq)
                    if L_prom == None:
                        seq_cropped=seq
                    else:
                        if L_seq > L_prom: 
                            seq_cropped = seq[(L_seq-L_prom):]
                        else: 
                            seq_cropped = seq
                    output_row = []
                    for motif_name, motif in motif_dict.items():
                        if output_format=='count': 
                            prom_hits = exact_promoter_scan(motif, seq_cropped, output_format=output_format, sequence_context = sequence_context)
                            output_row.append(prom_hits)   #append counts
                        elif output_format=='full': 
                            prom_hits = exact_promoter_scan(motif, seq_cropped, output_format=output_format, sequence_context = sequence_context)
                            #add output in the order of columns
                            counts = len(prom_hits)
                            output_row.append(counts)     #append count
                            output_row.append(prom_hits)  #append full_features

                    output_dict[gene] = output_row
                else:
                    print(gene + ' not in promoter_dict')
            if len(genes)==1:
                chosen_orths.append(gene)
            else:
                if rule=='max_STRE':
                    stre_counts = []
                    for gene in genes:
                        stre_counts.append(output_dict[gene][stre_ind])
                    chosen_orths.append(genes[np.argmax(stre_counts)])  #Chooses the one with the most STRES, if there are two, chooses the first
                else: 
                    raise ValueError('No rule for choosing orths')


    output_df = pd.DataFrame.from_dict(output_dict, orient="index", columns = columns)
    
    return output_df, chosen_orths

def motif_scan_YGOB_specs(ohnologs_goi, seed_spec, spec_order_pre_WGH, spec_order_post_WGH, motif_dict, L_prom=700,output_format='full', sequence_context=2 ):
    #searches for STRE within a set of ohnologs
    #
    #inputs:
    #    ohnologs_goi:  Dataframe of ohnologs that we want to search
    #    seed_spec:     Species which defined the ohnologs_goi set.  this species genename should be in genename_low/high columns
    #    spec_order_pre/post_WGH: list of pre/post_WGH species (four letter abbrev) that we want to search
    #    motif_dict: dictionary of motifs to search
    #    L_prom, output_format, sequence_context:  motif search paramters we pass
    #
    #outputs: 
    #    ohnologs_goi:  Adds columns with counts and results for motif searches for each included species
    #    motif_calc:  returns dictionary with pvalues and counts for presence of motif for each input motif
    #
    #A file of promoters for the set of gois is also saved during the routine in the promoter_dir

    spec_sets = {'Post WGH low' : spec_order_post_WGH, 
                 'Post WGH high' : spec_order_post_WGH, 
                 'Pre WGH' : spec_order_pre_WGH} 

    levels = {'Post WGH low': 'low', 
              'Post WGH high': 'high', 
              'Pre WGH' : ''}

    sc_kl_abbrev_lookup = {'Scer': 'sc', 'Klac': 'kl'}

    motif_calcs = {}  #already set up in S.Cer/K.Lac routine above
    for spec_set_name, spec_set in spec_sets.items():
        print(spec_set_name)
        level = levels[spec_set_name]
        if level == '': 
            level_sep = ''
        else: 
            level_sep = '_'
        motif_calcs[spec_set_name] = {}  #already set up in S.Cer/K.Lac routine above
        for spec in spec_set:
            print(spec)
            if spec in {'Klac', 'Scer'}:
                promoter_dir = data_processing_dir + os.path.normpath(sc_kl_abbrev_lookup[spec] + '_promoters/promoter_sets')
                all_promoters_fname = os.path.normpath(promoter_dir + '/all_' + sc_kl_abbrev_lookup[spec] + '_promoters.fasta')
            else: 
                promoter_dir = data_processing_dir + os.path.normpath('promoter_phylogeny/promoter_sets/' + spec)
                all_promoters_fname = os.path.normpath(promoter_dir + '/all_promoters_' + str(L_prom) + '.fasta')
            all_promoters = SeqIO.to_dict(SeqIO.parse(all_promoters_fname, "fasta"))  

            #scan all promoters for motifs
            all_promoters_scan = exact_promoter_scan_from_fasta(all_promoters_fname, 
                                                      motif_dict, 
                                                      output_format = 'full', 
                                                      sequence_context = 2, 
                                                      L_prom = None)
            motif_calcs_spec = {}
            motif_calcs_spec['all'] = {'total':len(all_promoters_scan)}
            for motif_name in motif_dict.keys():
                motif_calcs_spec['all'][motif_name]={'hits': sum(all_promoters_scan[motif_name + '_count']>0)}
                motif_calcs_spec['all'][motif_name]['pct'] = (motif_calcs_spec['all'][motif_name]['hits'])/(motif_calcs_spec['all']['total'])

            #     print('N ' + motif_name + ' in promoters of ' + spec + ' : ' + str(N_all_hits))
            #     print('N total promoters for ' + spec + ' : ' + str(N_all_total))

            genes = []


            #for each level, scan for motifs and add a column to the ohnologs_goi dataset.  Also saves a
            #promoter set file for gois, and outputs enrichment calculations

            if spec == seed_spec: 
                genes = ohnologs_goi['genename' + level_sep + level]
            elif spec_set_name == 'Pre WGH' :
                #Get goi orthologs for this species
                #load ortholog mapping
                orth_lookup = read_orth_lookup_table(seed_spec, spec, data_processing_dir + os.sep + "ortholog_files_YGOB" + os.sep)

                #for pre_wgh species checks to see if ortholog matches
                for row in ohnologs_goi.loc[:,['genename_low','genename_high']].iterrows(): 
                    orth_low = orth_lookup[row[1].genename_low][0]
                    orth_high = orth_lookup[row[1].genename_high][0]
                    if ((len(orth_lookup[row[1].genename_low])==2)|(len(orth_lookup[row[1].genename_low])==2)):
                        raise ValueError('more than one ortholog for a given ' + seed_spec + ' gene :' + row[1].genename_low + ' -> ' + orth_lookup[row[1].genename_low] +
                                         ', ' + row[1].genename_high + ' -> ' + orth_lookup[row[1].genename_low])
                    if orth_low != orth_high:
                        raise ValueError('Orthologs of WGH paralogs do not match : ' + row[1].genename_low + ' -> ' + orth_low +
                                         ', ' + row[1].genename_high + ' -> ' + orth_high) # + '. Used high ortholog')
                        #genes.append(orth_high)
                    else:
                        genes.append(orth_high)
            else: 
                #Get goi orthologs for this species
                #load ortholog mapping
                orth_lookup = read_orth_lookup_table(seed_spec, spec, data_processing_dir + os.sep + "ortholog_files_YGOB" + os.sep)

                for seed_gene in ohnologs_goi['genename' + level_sep + level]:
                    gene = orth_lookup[seed_gene]
                    if len(gene)==2:
                        print(gene + ' has more than one ortholog for ' + seed_gene)
                    genes.append(gene[0])

            ohnologs_goi[spec + '_genename' + level_sep + level] = genes

            goi_promoters_fname = promoter_dir + os.sep + 'depka_' + seed_spec + level_sep + level + '_' + str(L_prom) + '.fasta'

            seqs = []
            motif_scan_results = {motif_name:[] for motif_name in motif_dict.keys()}
            motif_scan_counts = {motif_name:[] for motif_name in motif_dict.keys()}
            for gene in genes:
                if gene != 'NONE':
                    try: 
                        seq = all_promoters[gene]
                        seqs.append(seq)
                        for motif_name, motif in motif_dict.items():
                            motif_scan_result = exact_promoter_scan(motif, seq, output_format=output_format, sequence_context=2)
                            motif_scan_results[motif_name].append(motif_scan_result)
                            motif_scan_counts[motif_name].append(len(motif_scan_result))
                    except KeyError:
                        print(gene + ' not present in promoter dict for ' + spec)
                        for motif_name, motif in motif_dict.items():
                            motif_scan_results[motif_name].append('NONE')
                            motif_scan_counts[motif_name].append(np.nan)
                else:
                    for motif_name, motif in motif_dict.items():
                        motif_scan_results[motif_name].append('NONE')
                        motif_scan_counts[motif_name].append(np.nan)

            #outputs a filename for a given promoter set
            with open(goi_promoters_fname, 'w') as f: 
                SeqIO.write(seqs, f, 'fasta')


            #Add column to dataframe and make enrichment calculations. 
            orthologs_present = [gene for gene in genes if gene!='NONE']
            motif_calcs_spec['goi'] = {'total': len(orthologs_present)}
            N2_total = motif_calcs_spec['all']['total'] - motif_calcs_spec['goi']['total']
            for motif_name in motif_dict.keys(): 
                ohnologs_goi[spec + '_' + motif_name + '_count' + level_sep + level] = motif_scan_counts[motif_name]
                ohnologs_goi[spec + '_' + motif_name + '_result'+ level_sep + level] = motif_scan_results[motif_name]

                motif_gt0 = [n_motif for n_motif in motif_scan_counts[motif_name] if ((not(np.isnan(n_motif))) & (n_motif>0))]
                motif_calcs_spec['goi'][motif_name] = {'hits':len(motif_gt0)}
                motif_calcs_spec['goi'][motif_name]['pct'] = (motif_calcs_spec['goi'][motif_name]['hits'])/(motif_calcs_spec['goi']['total'])
                N2_hits = motif_calcs_spec['all'][motif_name]['hits'] - motif_calcs_spec['goi'][motif_name]['hits']

                oddsratio, pvalue = stats.fisher_exact([
                                                        [motif_calcs_spec['goi'][motif_name]['hits'], N2_hits],
                                                        [motif_calcs_spec['goi']['total'], N2_total]
                                                       ], 
                                                       alternative = 'two-sided')        
                motif_calcs_spec['goi'][motif_name]['pval'] = pvalue

            motif_calcs[spec_set_name][spec] = motif_calcs_spec

    return ohnologs_goi, motif_calcs

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

def generate_dubious_orf_list(db_created, savefile): 
    #Generates a series which is a list of dubious orfs from the S. cerevisiae genome file. 
    #If running this for the first time, db_created will need to be False.  Otherwise it should be True
    #the database should be at: "genomes/scer_20181114/saccharomyces_cerevisiae_R64-2-2_20170117.db"
    #If you want to save the file again, savefile should be True
    
    gtf_fname = data_processing_dir + os.path.normpath("genomes/scer_20181114/saccharomyces_cerevisiae_R64-2-2_20170117.gff")
    db_fname = data_processing_dir + os.path.normpath("genomes/scer_20181114/saccharomyces_cerevisiae_R64-2-2_20170117.db")
 
    if db_created: 
        db = gffutils.FeatureDB(db_fname)
    else: 
        db = gffutils.create_db(gtf_fname, db_fname)

    sc_features = pd.read_sql('select * from features;', db.conn)

    sc_genes = sc_features[sc_features['featuretype']=='gene']

    dubious_orfs = []
    for row in sc_genes.iterrows():
        attributes = row[1]['attributes'].strip("{}").split(",")
        for attribute in attributes: 
            if attribute.split(":")[0]=='"orf_classification"':
                if attribute.split(":")[1]=='["Dubious"]':
                    dubious_orfs.append(row[1]['id'])
    print("There are {} dubious orfs in the saccharomyces_cerevisiae_R64-2-2_20170117.gff".format(len(dubious_orfs)))
    
    dubious_orfs = pd.Series(dubious_orfs)
    
    if savefile: 
        dubious_orfs.to_csv(data_processing_dir + os.path.normpath("genomes/scer_20181114/dubious_orfs.csv"))
    
    return dubious_orfs

def wsl_filename_convert(fname):
    #Convert a filename on a windows system to the corresponding Windows Subsystem for Linux filename
    wsl_fname = '/' + '/'.join(['mnt','c'] + fname.split('\\')[1:])
    return wsl_fname

## Topology determination
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

#Helper functions
def tryfloatconvert(value, default):
    try:
        return float(value)
    except ValueError:
        return default

def quantileNormalize(df_input):
    #From: https://github.com/ShawnLYU/Quantile_Normalize
    #computes quantile normalization on a dataframe composed of floating point numbers
    #For instance the columns could all be replicates of a microarray experiment. 
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df
    
def correlation_nan_filt(v1, v2): 
    
    v1_nan = np.isnan(v1)
    v2_nan = np.isnan(v2)

    v1_v2_nan = v1_nan+v2_nan

    v1_filt = [item for jj, item in enumerate(v1) if ~v1_v2_nan[jj]]
    v2_filt = [item for jj, item in enumerate(v2) if ~v1_v2_nan[jj]]

    corr_dist = spd.correlation(v1_filt, v2_filt)
    
    #assign 1 
    if np.isnan(corr_dist):
        corr_dist=1.0  
        print("yeast_esr_exp.correlation_nan_filt.  Warning: encountered vector of all 0's.  Assigned correlation distance of 1.")
    
    return corr_dist
