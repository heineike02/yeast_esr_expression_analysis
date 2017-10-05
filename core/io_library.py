# -*- coding: utf-8 -*-
#from IPython.core.debugger import Tracer
import os 
import pandas as pd
import numpy as np
import re
import math
from collections import Counter
import subprocess
print('I am importing io_library')
#data_dir is a global variable where all data files are stored. 
# Gabe 7/12/16
#data_dir = os.path.normpath("C:\\Users\Ben\Documents\GitHub\expression_broad_data\expression_data")		
#data_dir = os.path.normpath(os.path.dirname(os.getcwd()) + '/scripts/expression_broad_data_datafiles/microarray_data/')
#data_dir = '/home/heineike/github/expression_broad_data/expression_data'
data_dir = os.path.normpath("C:\\Users\\heine\\github\\expression_broad_data\\expression_data")		

print(data_dir)

def tryfloatconvert(value, default):
    try:
        return float(value)
    except ValueError:
        return default

def parse_raw_exp(species):
    #Parses raw expression data for a given species.  Ouputs a dataframe and saves a .csv file. 
    #org_dict = {'Kluyveromyces lactis': ('GSE22198_family.soft','SystematicName','KLac',KL_soft_alt_dict), 'Saccharomyces cerevisiae': ('GSE22204_family.soft','ACCESSION_STRING','SCer',SC_soft_alt_dict)}
    raw_exp_datasets = {'Kluyveromyces lactis': 'GSE22198_family.soft', 'Saccharomyces cerevisiae': 'GSE22204_family.soft', 'Candida glabrata':'GSE22194_family.soft', 'Saccharomyces castellii' : 'GSE22200_family.soft', 'Saccharomyces bayanus' : 'GSE22205_family.soft'}
    #raw_exp_fnames = {'Kluyveromyces lactis': 'KLac_raw.csv', 'Saccharomyces cerevisiae': 'SCer_raw.csv'}
    orf_header_names = {'Kluyveromyces lactis': 'SystematicName', 'Saccharomyces cerevisiae': 'ACCESSION_STRING', 'Candida glabrata':'SystematicName', 'Saccharomyces castellii' : 'SystematicName', 'Saccharomyces bayanus' : 'SystematicName'}
    
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
    input_fname = os.path.normpath(data_dir + '/raw_exp/' + raw_exp_datasets[species])
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
        line = f.next()
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
            f.next()
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
    platform_prefixes = {'Kluyveromyces lactis': '10499', 'Saccharomyces cerevisiae': '9294', 'Candida glabrata':'10497', 'Saccharomyces castellii' : '10501', 'Saccharomyces bayanus' : '10505'}
    
    input_fname = os.path.normpath(data_dir + '/GSE' + exptype_prefixes[exptype] + '_' + exptype + '/GSE'+ exptype_prefixes[exptype] + '-GPL' + platform_prefixes[species] + '_series_matrix.txt')
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
        #sort ids alphanumerically
        expdata = expdata.sort()
        #Sets a multi index for conditions and replicates
        col_mindex_arrays = [conditions,replicates,condition_ids]
        col_mindex_tuples = list(zip(*col_mindex_arrays))
        col_mindex = pd.MultiIndex.from_tuples(col_mindex_tuples, names=['conditions', 'replicates','array_ids'])
        expdata.columns = col_mindex
        #sort by conditions
        expdata_sorted = expdata.sortlevel('conditions',axis = 'columns')
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
    
    for jj in range(len(species_list)): 
    
        #Generate raw expression data
        raw_exp, orf_lookup = parse_raw_exp(species_list[jj])
        
        #save raw expression data to a csv file
        # Gabe 7/12/16
        # fname = os.path.normpath(base_dir + "\microarray_data\\raw_exp\\"  + fname_out_bases[jj] + '_raw_exp.csv')
        fname = os.path.normpath(base_dir + "/expression_data/raw_exp/"  + fname_out_bases[jj] + '_raw_exp.csv')
        raw_exp.to_csv(fname)
        print(fname + ' saved')
        
        #Generate data for microarrays
        growth_exp = parse_micro_data(species_list[jj],'Growth',orf_lookup)
        # Gabe 7/12/16
        # fname = os.path.normpath(base_dir + "\microarray_data\\GSE36253_Growth\\"  + fname_out_bases[jj] + '_growth.csv' )
        fname = os.path.normpath(base_dir + "/expression_data/GSE36253_Growth/"  + fname_out_bases[jj] + '_growth.csv' )
        growth_exp.to_csv(fname)
        print(fname + ' saved')
        
        if species_list[jj] != 'Saccharomyces bayanus':
            #There is no stress dataset for S. bayanus
            stress_exp = parse_micro_data(species_list[jj],'Stress',orf_lookup)
            # Gabe 7/12/16
            # fname = os.path.normpath(base_dir + "\microarray_data\\GSE38478_Stress\\"  + fname_out_bases[jj] + '_stress.csv' )
            fname = os.path.normpath(base_dir + "/expression_data/GSE38478_Stress/"  + fname_out_bases[jj] + '_stress.csv' )
            stress_exp.to_csv(fname)
            print(fname + ' saved')
    
    return 
    
def read_SGD_features():
    
    #Read in orf/name file and make it a dictionary
    # Gabe 7/12/16
    # SC_features_fname = os.path.normpath(data_dir + "\ortholog_files\\SGD_features.tab")
    SC_features_fname = os.path.normpath(data_dir + "/ortholog_files/SGD_features.tab")

    SC_features = pd.read_csv(SC_features_fname, sep = '\t', header=None)
    SC_orfs = SC_features.groupby(1).get_group('ORF')
    
    #Makes a dictionary to look up orfs by gene names.  This won't include all orfs - those without names had NaN in column 4 so 
    #are presumably left out. 
    SC_orfs_lookup = dict(zip(SC_orfs[4], SC_orfs[3]))
    SC_genename_lookup = dict(zip(SC_orfs[3], SC_orfs[4]))
    SC_features_lookup = dict(zip(SC_orfs[3], SC_orfs[15]))
       
    return SC_orfs_lookup, SC_genename_lookup, SC_features_lookup

def read_orth_lookup_table(species1, species2, orth_dir):
    #For a given species read in the ortholog file, make a dictionary
    orth_file_abbrev = {'Kluyveromyces lactis': 'Klac', 'Saccharomyces cerevisiae': 'Scer', 'Candida glabrata':'Cgla', 'Saccharomyces castellii' : 'Scas', 'Saccharomyces bayanus' : 'Sbay'}
    orth_fname = orth_file_abbrev[species1] + "-" + orth_file_abbrev[species2] + "-orthologs.txt"
    orth_fname = os.path.normpath(orth_dir + os.sep + orth_fname)
    
    #There are some orfs that have multiple orthologs - in that case both will be used
    with open(orth_fname) as f:
        orth_lookup = {}
        for line in f:
            linesp = line.split()
            orth_lookup[linesp[0]]= linesp[1:]

    return orth_lookup  
    
def get_gasch_ESR_list(act_rep):
    #For a file from the Gasch 2000 supplement, read in the data
    fname = os.path.normpath(data_dir + "/gasch_data/gasch_fig3_" + act_rep + "_ESR.txt")
    
    with open(fname) as f:
        out = []
        #Skips first two lines
        next(f)
        next(f)
        for line in f:
            linesp = line.split()
            out.append(linesp[0])

    return out
    
def read_gasch_data(conditions,fname):
    #For selected conditions from the Gasch 2000 supplement, extract the data as 
    #a dataframe with a row for each gene and each column a condition
    
    if fname == "gasch_complete_dataset.txt":
        gene_ind = 0
    else: 
        gene_ind = 1
    fname_full = os.path.normpath(data_dir + "/gasch_data/" + fname)
    
    with open(fname_full) as f:
        header = next(f).split("\t")
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
    series_fname = GEO_accession + '_series_matrix.txt'
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

    #Find line that starts table listing gene names and index numbers

    series_fname = os.path.normpath(data_dir + series_fname)
    #GSM1423542: nmpp1 treatment (20 min) - ACY142 +nmpp1 / ACY142

    with open(series_fname) as f:
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
    
    oshea_exp_data_dir = data_dir + '\GSE32703_NMPP1_SC\\' 
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
    fname_solis_SC_PKA_data = os.path.normpath(data_dir + '\SCer_NMPP1_RNA_Seq\solis_2016.xlsx')
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

def build_motif_dict(fname): 
    motif_dict = {}
    with open(fname,'r') as f: 
        for line in f: 
            if line != '\n':
                if line.split()[0]=='MOTIF':
                    motif_dict[line.split()[1]]=line.split()[2]
    return motif_dict

def write_ame_promoter_file(promoter_database, gene_list,fname):
        
    with open(fname,'w') as f: 
        for gene in gene_list:
            try: 
                row = promoter_database.loc[gene,]
                header_line = '>' + row.name + ' 700bp_upstream\n'
                seq_line = row['prom_seq'] + '\n'
                f.write(header_line)
                f.write(seq_line)
            except KeyError: 
                print(gene + " not in promoter data set.")
    
    return

def read_ame_output(fname, motif_dict):
    #reads in ame program output file from meme suits
    #assumes first 13 lines are not data
    motif_name_list = []
    pval_list = []
    with open(fname,'r') as f: 
        for jj in range(1,13):
            f.next()
        for line in f:
            #Add motif name
            motif_id = line.split()[5]
            motif_name_list.append(motif_dict[motif_id])
            #Add two tailed corrected p-value to list (why two tailed? what is U-value? What is corrected p-Value?)
            pval_list.append(line.split()[-1][0:-1])
    
    
    ame_dict = {"motif_name": motif_name_list, "pval":pval_list}
    ame_data = pd.DataFrame.from_dict(ame_dict)
    
    return ame_data

def run_ame_analysis(spec, target_gene_list, control_gene_list,target_fname_prefix, control_fname_prefix, motif, 
                     promoter_dir = {'KL': data_dir+'/kl_promoters/' , 'SC': data_dir + '/sc_promoters/'},
                     promoter_fname = {'KL': 'kl_promoters.pkl', 'SC': 'sc_promoters.pkl'},
                     ame_scoring = 'totalhits',
                     ame_method = 'fisher',
                     ame_pvalue_threshold = '0.05' ):
    #runs ame program from meme software for given target gene lit, control gene list and set of motifs. 
    #extract promoters
    promoters = pd.read_pickle(promoter_dir[spec] + promoter_fname[spec]) 
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
    motif_db = promoter_dir[spec]+ motif['fname']
    target_sequences = promoter_dir[spec] + 'promoter_sets/' + fname_prefixes['target'] + '_promoters.fasta'
    control_sequences = promoter_dir[spec] + 'promoter_sets/' + fname_prefixes['control'] + '_promoters.fasta'
    output_dir = promoter_dir[spec] + 'ame_output'
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




def make_foldchange_subsets(kl_sc_PKA_data, pthreshold_KL, pthreshold_SC): 
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
    #Input is dataframe with GFF version of kl genename as the index. 
    #Ouput is dataframe with standard version of kl genename as the index. 

    kl_genename = []
    for gene in df.index: 
        if gene[0:5]=='KLLA0':
            new_gene = gene.split('_')[0]+gene.split('_')[1]
        else: 
            new_gene = gene
        kl_genename.append(new_gene)
    
    df['kl_genename'] = kl_genename
    df.set_index('kl_genename',inplace = True)

    return df
    
def load_goslim_data(GO_aspect):
    #The three GO_aspect values are: 
    #C = cellular_component
    #F = molecular_function
    #P = biological_process
    go_slims = pd.read_table(data_dir + os.path.normpath('/go_terms/go_slim_mapping.tab'),header = None)
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