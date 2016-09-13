from IPython.core.debugger import Tracer
import os 
import pandas as pd
import numpy as np
import re

#data_dir is a global variable where all data files are stored. 
# Gabe 7/12/16
data_dir = os.path.normpath("C:\Users\Ben\Documents\GitHub\expression_broad_data\microarray_data")		
#data_dir = os.path.normpath(os.path.dirname(os.getcwd()) + '/scripts/expression_broad_data_datafiles/microarray_data/')

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
        # Gabe 7/12/16
        # fname = os.path.normpath(base_dir + "\microarray_data\\raw_exp\\"  + fname_out_bases[jj] + '_raw_exp.csv')
        fname = os.path.normpath(base_dir + "/microarray_data/raw_exp/"  + fname_out_bases[jj] + '_raw_exp.csv')
        raw_exp.to_csv(fname)
        print fname + ' saved'
        
        #Generate data for microarrays
        growth_exp = parse_micro_data(species_list[jj],'Growth',orf_lookup)
        # Gabe 7/12/16
        # fname = os.path.normpath(base_dir + "\microarray_data\\GSE36253_Growth\\"  + fname_out_bases[jj] + '_growth.csv' )
        fname = os.path.normpath(base_dir + "/microarray_data/GSE36253_Growth/"  + fname_out_bases[jj] + '_growth.csv' )
        growth_exp.to_csv(fname)
        print fname + ' saved'
        
        if species_list[jj] != 'Saccharomyces bayanus':
            #There is no stress dataset for S. bayanus
            stress_exp = parse_micro_data(species_list[jj],'Stress',orf_lookup)
            # Gabe 7/12/16
            # fname = os.path.normpath(base_dir + "\microarray_data\\GSE38478_Stress\\"  + fname_out_bases[jj] + '_stress.csv' )
            fname = os.path.normpath(base_dir + "/microarray_data/GSE38478_Stress/"  + fname_out_bases[jj] + '_stress.csv' )
            stress_exp.to_csv(fname)
            print fname + ' saved'
    
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
       
    return SC_orfs_lookup, SC_genename_lookup

def read_orth_lookup_table(species1, species2):
    #For a given species read in the ortholog file, make a dictionary
    orth_file_abbrev = {'Kluyveromyces lactis': 'Klac', 'Saccharomyces cerevisiae': 'Scer', 'Candida glabrata':'Cgla', 'Saccharomyces castellii' : 'Scas', 'Saccharomyces bayanus' : 'Sbay'}
    orth_fname = orth_file_abbrev[species1] + "-" + orth_file_abbrev[species2] + "-orthologs.txt"
    # Gabe 7/12/16
    # orth_fname = os.path.normpath(data_dir + "\ortholog_files\\" + orth_fname)
    #orth_fname = os.path.normpath(data_dir + "/ortholog_files/" + orth_fname)
    orth_fname = os.path.normpath(data_dir + "/ortholog_files_YGOB/" + orth_fname)
    
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
    
def read_gasch_data(conditions):
    #For selected conditions from the Gasch 2000 supplement, extract the data as 
    #a dataframe with a row for each gene and each column a condition
    fname = os.path.normpath(data_dir + "/gasch_data/gasch_fig1_all_conditions.txt")
    
    with open(fname) as f:
        header = next(f).split("\t")
        condition_inds = [header.index(condition) for condition in conditions]
        #Skips first two lines
        next(f)
        exp_data = []
        gene_names = []
        for line in f:
            exp_values = line.split("\t")
            exp_data.append([tryfloatconvert(exp_values[condition_ind],None) for condition_ind in condition_inds])
            gene_names.append(exp_values[1])
        exp_df = pd.DataFrame(data=exp_data, index=gene_names, columns=conditions) 
    
    return exp_df

def parse_data_osheaNMPP1(desired_conditions): 
    #load raw data for NMPP1 experiments from Oshea Paper
    
    #Extract dictionary for each gene from family.soft file
    oshea_pka_dir = 'C:\Users\Ben\Documents\GitHub\expression_broad_data\microarray_data\GSE32703_NMPP1_SC'
    soft_fname = os.path.normpath(oshea_pka_dir + '/GSE32703_family.soft')
    with open(soft_fname) as f:
        for line in f: 
            if line.split()[0] == '!platform_table_begin':
                break
            
        #Find index of header that will identify orf
        line = f.next()
        linesp = line.split()
        orf_header_ind = linesp.index('ORF')
        #skip first three lines to get to data
        f.next()
        f.next()
        f.next()
        #data_dict = {}
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
    
    soft_fname = os.path.normpath(oshea_pka_dir + '/GSE32703_series_matrix.txt')
    #GSM812516: No NMPP1 0min
    #GSM812520: 40 min: 3uM 1-NMPP1
    #Want the ratio of 40min to 0min

    with open(soft_fname) as f:
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
    orth_lookup_outputfname = os.path.normpath(base_dir + '\microarray_data\ortholog_files_YGOB\\' + orth_file_abbrev[species1] + "-" + orth_file_abbrev[species2] + "-orthologs.txt"  )
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