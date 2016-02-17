from IPython.core.debugger import Tracer
import os 
import pandas as pd
import numpy as np
import re

#data_dir is a global variable where all data files are stored. 
data_dir = os.path.normpath("C:\Users\Ben\Documents\GitHub\expression_broad_data\microarray_data")		

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
                    print 'Species: {}, Line {} was blank, made nan'.format(species, linesp[0])
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
                    conditions = [re.split("[\-\[\]]",title)[1] for title in condition_titles.split('\t')[1:-1]]
                    replicates = [re.split("[\-\[\]]",title)[2] for title in condition_titles.split('\t')[1:-1]]
                    break
        
        elif exptype == 'Stress': 
            
            #Make replicate list 
            for line in f:
                if line == '\n':
                    #skips line if it is just a new line
                    pass
                elif line.split()[0]== '!Sample_source_name_ch1':
                    replicate_titles = line
                    replicates = [re.split("[\"_]",title)[2] for title in replicate_titles.split('\t')[1:-1]]
                    break
                    
            #Make condition list by concatenating stress and time information for each microarray
            for line in f:
                if line.split()[0]== '!Sample_characteristics_ch1':
                    condition_stress_titles = line
                    condition_stresses = [re.split("[\"\:]",title)[2].strip() for title in condition_stress_titles.split('\t')[1:-1]]
                    break
            for line in f:
                if line.split()[0]== '!Sample_characteristics_ch1':
                    condition_time_titles = line
                    condition_times = [re.split("[\"\:]",title)[2].strip() for title in condition_time_titles.split('\t')[1:-1]]
                    break
            conditions = ['{}_{:03d}'.format(tup[0],int(tup[1])) for tup in zip(condition_stresses,condition_times)]

            
                
        #scroll to beginning of data table and extract ids for each experiment
        for line in f:
            if line.split()[0]== '"ID_REF"':
                condition_ids = line.split()[1:-1]
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
            chip_location_values = linesp[1:-1]
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
            print "Error: ID mismatch between experiment data and orf lookup table. Species = {}, Experiment Type = {}".format(species, exptype)
            return
        
        #if no error, continue here
        print "All ID's match between experiment data and orf lookup table. Species = {}, Experiment Type = {}".format(species, exptype)
        
        mindex_arrays = [np.array(orf_lookup.index),np.array(orf_lookup)]
        mindex_tuples = list(zip(*mindex_arrays))
        expdata_sorted.index = pd.MultiIndex.from_tuples(mindex_tuples, names=['ID','orf_name'])
        

          
        

    
    return expdata_sorted
    

def make_data_tables(species_list,fname_out_bases, base_dir):
    
    for jj in range(len(species_list)): 
    
        #Generate raw expression data
        raw_exp, orf_lookup = parse_raw_exp(species_list[jj])
        
        #save raw expression data to a csv file
        fname = os.path.normpath(base_dir + "\microarray_data\\raw_exp\\"  + fname_out_bases[jj] + '_raw_exp.csv')
        raw_exp.to_csv(fname)
        print fname + ' saved'
        
        #Generate data for microarrays
        growth_exp = parse_micro_data(species_list[jj],'Growth',orf_lookup)
        fname = os.path.normpath(base_dir + "\microarray_data\\GSE36253_Growth\\"  + fname_out_bases[jj] + '_growth.csv' )
        growth_exp.to_csv(fname)
        print fname + ' saved'
        
        if species_list[jj] != 'Saccharomyces bayanus':
            #There is no stress dataset for S. bayanus
            stress_exp = parse_micro_data(species_list[jj],'Stress',orf_lookup)
            fname = os.path.normpath(base_dir + "\microarray_data\\GSE38478_Stress\\"  + fname_out_bases[jj] + '_stress.csv' )
            stress_exp.to_csv(fname)
            print fname + ' saved'
    
    return 
        
        
  
