import os
import matplotlib.pyplot as plt
import numpy as np
import copy 
from scipy import optimize
from IPython.core.debugger import Tracer

class Orf:  
  #Class for an ORF in the K.Lactis Genome 
  def __init__(self, orf, no):
    assert type(orf) in [str, unicode], 'orf name should be a string'
    #assert type(sequence) in [str, unicode], 'sequence should be a string'
    self.orf = orf
    self.no = no
    self.sequence_pm2000 = ''
    self.sequence = ''
    self.exp_raw = []
    self.exp_mean = []
    self.exp_var = []
    self.conditions = {}
    self.variability = 0.0
    # Lag, Log, Late log (LL), Diauxic shift (DS), Post shift (PS), and plateau (P)
    #Make a funciton that calculates mean and variance
  
  def exp_calcs(self):
    self.exp_mean = np.mean(self.exp_raw)
    self.exp_var = np.var(self.exp_raw)
  
  def var_calcs(self, cond_order, weights, na_penalty): 
      #weights is a vector of constants for each condition that adjusts how important each condition is
      #na_penalty is what is used in the calculation if there are NA values.  This will be set to the mean value for all non-NA values for each condition in the entire dataset.

      #Need to make sure it is in the right order - dictionary is unordered...
      var_vec = []
      for jj in range(len(cond_order)): 
          cond_val = self.conditions[cond_order[jj]]
          if cond_val == 'NA':
              var_vec.append((weights[jj]*na_penalty[jj])**2)
          else:
              var_vec.append((weights[jj]*cond_val)**2)
                  
      self.variability = sum(var_vec)
                      
  def __repr__(self):
    # Overload the __repr__ operator to make printing simpler.
    return 'Orf %s Number %s ' %(self.orf ,self.no)


class OrfDict(dict):
    
    #Class for holding multiple orfs
    def __init__(self): 
        dict.__init__(self)
        self.expression_experiments = []
   
    def add_orf(self, id, name, sequence):
        if id in self: 
            print 'Orf %s already in dictionary, ignoring' %id
            return
        self[id] = Orf(id, name, sequence)

def orf_cond_mean(orf, conditions_sorted): 
    new_orf = copy.deepcopy(orf)
    new_orf.conditions = {}
    old_cond = ''
    unique_conds = []
    for cond_full in conditions_sorted:
        cond = cond_full.split(',')[0]
        if cond != old_cond: 
            unique_conds.append(cond)
            new_orf.conditions[cond] = []
        try: 
            cond_val = orf.conditions[cond_full]
        except KeyError: 
            print 'KeyError gene %s does not have data for condition %s' %(orf.orf, cond_full)
            cond_val = 'NA'
        #First store all non-NA values for each new unique condition in a list
        if cond_val != 'NA':
            new_orf.conditions[cond].append(cond_val)
        old_cond = cond
     
    #now replace all lists with values
    for cond in unique_conds: 
        cond_list = new_orf.conditions[cond]
        if len(cond_list)>0: 
            new_orf.conditions[cond] = np.mean(cond_list)
        else: 
            new_orf.conditions[cond] = 'NA'
    
    return new_orf, unique_conds

def mean_NAs(list_in): 
    list_len = 0.0
    list_sum = 0.0
    for item in list_in: 
        if item != 'NA': 
            list_len = list_len + 1.0
            list_sum = list_sum + item
    list_mean = list_sum/list_len
    
    return list_mean

def read_orf_dict(fname,orf_header_name,organism):
    #Note: For gene 11892 (line 7334), gene 14495 (line 9937), gene 13495 (Line 14311), gene 14052 (line 14868), gene 10964 (line 17154), gene 12500 (line 18690) there was no expression value, so I made it NA.  
    lines = open(fname,'r').readlines()
    
    #Find names of experiments in the file
    jj = 0
    line = lines[jj]
    word1 = line.split()[0]
    while word1 != '!Series_sample_id': 
        jj = jj+1
        line = lines[jj]
        if len(line.split())>0:
            word1 = line.split()[0]  
        else:
            word1 = 'no line' 
    
    experiments = []
    
    while word1 == '!Series_sample_id':
        experiment = line.split()[2]
        experiments.append(experiment)
        jj = jj+1
        line = lines[jj]
        if len(line.split())>0:
            word1 = line.split()[0]  
        else:
            word1 = 'no line' 

    print experiments
    
    #Find Line that starts listing names of genes
     
    line = lines[jj]
    word1 = line.split()[0]
    while word1 != 'ID':
        jj = jj+1
        line = lines[jj]
        if len(line.split())>0:
            word1 = line.split()[0]  
        else:
            word1 = 'no line' 
    
    #Find index of header that will identify orf
    linesp = line.split()
    orf_ind = linesp.index(orf_header_name)
    
    #Make dictionary of numbers to gene names
    jj = jj + 1
    line = lines[jj]
    linesp = line.split('\t')
    word1 = linesp[0]  
    gene_dict = {}
    orf_dict = OrfDict()
    orf_dict.expression_experiments = experiments
    
    #Could get this from file - line starts with ^SAMPLE 
    
    while word1!= '!platform_table_end\n':   
        gene_no = linesp[0]
        if organism == 'Kluyveromyces lactis':
            gene_name = linesp[orf_ind]
        elif organism == 'Saccharomyces cerevisiae':
            gene_name = linesp[orf_ind].split('|')[1]
        gene_dict[gene_no] = gene_name  
        orf_dict[gene_name] = Orf(gene_name,gene_no)                                      
        jj = jj + 1
        line = lines[jj]
        linesp = line.split('\t')
        word1 = linesp[0]

    #Make list of genes with experimental data from each experiment
    for kk in range(len(experiments)):
        #Find Line that starts listing exp values of genes
        while word1 != 'ID_REF':
            jj = jj+1
            line = lines[jj]
            if len(line.split())>0:
                word1 = line.split()[0]  
            else:
                word1 = 'no line' 
        
        jj = jj+1
        line = lines[jj]
        linesp = line.split()
        word1 = linesp[0]
        
        #Make dictionary of Orfs
        while word1!= '!sample_table_end':
            gene_no = linesp[0]
            gene_name = gene_dict[gene_no]
            orf = orf_dict[gene_name]
            if linesp[1] != 'NA':
                orf.exp_raw.append(float(linesp[1]))
            jj = jj+1
            line = lines[jj]
            linesp = line.split()
            word1 = linesp[0]
    
    return orf_dict, gene_dict
  
      
def read_micro_data(orf_dict, fname, gene_dict, organism,orf_header_name,org_soft_alt_dict,exp_type):
    #Note: This makes changes to the orfs in orf_dict
    #exp_type is either 'Growth' or 'HeatOx' and determines how experimental conditions are parsed
    with open(fname,'r') as f: 
        jj = 0
        
        line = f.readline()
        
        #Verify that ID numbers are the same for all genes
        
        #advance to line where platform is mentioned
        while line != ('!Platform_organism = ' + organism + '\n'):
            line = f.readline()
        
        #advance to line where table begins
        while line != '!platform_table_begin\n':
            line = f.readline()
        
        line = f.readline()
        line = line.split()
        #Find index of header that will identify orf
        orf_ind = line.index(orf_header_name) 
        line = f.readline()
        line = line.split('\t')
                
        #advance until platform table ends
        
        
        while not('!platform_table_end' in line[0]):
            if organism == 'Kluyveromyces lactis':
                gene_name = line[orf_ind]
            elif organism == 'Saccharomyces cerevisiae':
                try: 
                    gene_name = line[orf_ind].split('|')[1]
                except IndexError:
                    gene_name = org_soft_alt_dict[line[0]] 
            if orf_dict[gene_name].no != line[0]:
                print 'Error in Line:'
                print gene_name
                print orf_dict[gene_name].no
                print line
            line = f.readline()
            line = line.split('\t')
        #Note: several of the oligos on the S.C. Chip correspond to the same ORF. We only use one for these calculations.
        print 'Gene expression chip same as varying condition chip'    
        
        exp_conds = []
        
        while line != '':
            #advance to line where sample is mentioned    
            while  line != '!Sample_organism_ch1 = ' + organism + '\n':
            #while  line != '!Sample_organism_ch1 = Candida glabrata\n':
                line = f.readline()
                if line == '':
                    print 'Got To End'
                    f.close()
                    return exp_conds
            
            #for jj in range(40):
            #    line = f.readline()
            #    print line
                
                    
            if exp_type == 'Growth':
                while not('Sample_description =' in line):
                    line = f.readline()
                linesp = line.split()
                exp_cond = linesp[2]
            
                line = f.readline()
                #exp_cond = exp_cond + ', bio rep ' + linesp[4] + ' tech rep ' + linesp[7] 
                exp_cond = exp_cond + ', ' + line
                exp_conds.append(exp_cond)
                print exp_cond
            elif exp_type == 'HeatOx':
                while not('!Sample_characteristics_ch1' in line):
                    line = f.readline()
                #First line is in the format: '!Sample_characteristics_ch1 = treatment: hydrogen peroxide'
                linesp = line.replace('\n','').split(': ')
                if linesp[-1] == 'hydrogen peroxide':
                    exp_cond = 'H2O2'
                elif linesp[-1] == 'heat shock':
                    exp_cond = 'HS'
                elif linesp[-1] == 'NaCl': 
                    exp_cond = 'NaCl'      
                #second line is of the form !Sample_characteristics_ch1 = time point (minutes): 5
                line = f.readline()
                linesp = line.replace('\n','').split(': ')
                    
                exp_cond = exp_cond + ' T' + linesp[-1]
            
                while not('!Sample_description' in line):
                    line = f.readline()
                exp_cond = exp_cond + ', ' + line 
                exp_conds.append(exp_cond)
                print exp_cond
                                    
            while line != '!sample_table_begin\n': 
                line = f.readline()
            
            #if exp_cond == 'NaCl T5, !Sample_description = Biological replicate B, technical replicate 1\n':
            #    Tracer()()
            
            #advance one more line to avoid header
            line = f.readline()
            line = f.readline()
            linesp = line.split()
            
            while line != '!sample_table_end\n': 
                try: 
                    orf = orf_dict[gene_dict[linesp[0]]]
                    if orf.orf == 'YGR192C':
                        print '!!!!!!!!' + exp_cond
                    if len(linesp) == 2:
                        orf.conditions[exp_cond] = float(linesp[1])
                    elif len(linesp) == 1:
                        orf.conditions[exp_cond] = 'NA'
                        print 'No value, gene # %s, experiment %s' %(linesp[0], exp_cond)
                except KeyError: 
                    print 'No gene %s' %linesp[0]               
                line = f.readline()
                linesp = line.split()
                


        
def heatmap_plot(orfs, conds, row_labels, column_labels):
    
    column_labels = ['Exp'] + column_labels
    data_mat = np.zeros((len(orfs),len(conds)+1))
    
    for jj in range(len(orfs)): 
        orf = orfs[jj] 
        #First add raw expression
        data_mat[jj,0] = orf.exp_mean
        #Then plot log fold change    
        for kk in range(1,len(conds)+1):
            cond = conds[kk-1]
            val = orf.conditions[cond]
            if val == 'NA': 
                print 'NA changed to -7.0: orf %s, condition %s' %(orf.orf,cond)
                data_mat[jj,kk] = -7.0
            else: 
                data_mat[jj,kk] = val
            
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(data_mat, cmap=plt.cm.Spectral)
    
            
    # put the major ticks at the middle of each cell
    # x axis is the columns, y axis is the rows
    ax.set_xticks(np.arange(data_mat.shape[1])+0.5, minor=False)
    ax.set_yticks(np.arange(data_mat.shape[0])+0.5, minor=False)
    
    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    
    ax.set_xticklabels(column_labels, minor=False)
    ax.set_yticklabels(row_labels, minor=False)
    cbar = plt.colorbar(heatmap)
    plt.xticks(rotation=75)
    plt.show()
    
    return fig, ax
    
            
def main():
        #Add heat and ox stress datasets
        #Display the same for S.C. 
        
        #os.chdir(os.path.normpath('C://Users/Ben/Documents/GitHub/bioinformatics/bioinformatics_pycode/'))
       
        #Calculate values for S. Cerevisiae
	#SC soft file alterations: 
	#line 550 A_06_P1460 bmh|XX1
	#line 1577 A_06_P2487 bmh|XX2
        #line 1910 A_06_P2820 bmh|XX3
        #line 3059 A_06_P3969 bmh|XX4
        #line 3064 A_06_P3974 bmh|XX5
        #line 3127 A_06_P4037 bmh|XX6
        #line 4771 A_06_P5681 bmh|XX7
        #added 'sample_table_end\n!' to line 6437 
        #
        #Changes to 'GSE36253_family.soft'
        #added '!platform_table_end to line 6803 because last few lines of SC table were not genes (probably microarray controls)
        #
        #Changes to 'GSE38478_family.soft'
        #added '!platform_table_end to line 6566 because last few lines of SC table were not (probably microarray controls)
        #added A_06_P3413 in line 72374, 66068, 59762, 53456 because the line was blank.  They were all TDH3
        #Very strange that data for TDH3, a very popular promoter seems to have been removed from the NaCl experiments
        
        
        SC_soft_alt_dict = {'A_06_P1460': 'XX1', 'A_06_P2487':'XX2','A_06_P2820': 'XX3','A_06_P3969':'XX4', 'A_06_P3974': 'XX5', 'A_06_P4037': 'XX6', 'A_06_P5681' : 'XX7'}
        #Perhaps make this automatic for KL?
        KL_soft_alt_dict = {}
                        
        org_dict = {'Kluyveromyces lactis': ('GSE22198_family.soft','SystematicName','KLac',KL_soft_alt_dict), 'Saccharomyces cerevisiae': ('GSE22204_family.soft','ACCESSION_STRING','SCer',SC_soft_alt_dict)}
        #org_dict has the following values for each organism: 
        # File name for the gene dictionary
        # Tag used in gene dictionary file to identify primary key number
        # Abbreviation for species used in figures
        #Calculate values for K. Lactis	
        organism = 'Kluyveromyces lactis'
        #organism = 'Saccharomyces cerevisiae'	
        data_dir = os.path.normpath("C:/Users/Ben/Google Drive/UCSF/ElSamad_Lab/PKA/Bioinformatics/From_Desktop/bioinformatics/data/microarray_data")		
        fname = os.path.normpath(data_dir + '/' +org_dict[organism][0]) 
        orf_header_name = org_dict[organism][1]
        org_soft_alt_dict = org_dict[organism][3]
        orf_dict, gene_dict = read_orf_dict(fname,orf_header_name,organism)
        
        org_abbrev = org_dict[organism][2]
        

              	      

        #YIL018W, YOR063W, YPL131W, YLR448W
        
        
        for key in orf_dict.keys():
            orf = orf_dict[key]
            orf.exp_calcs()
               
        exp_mean = [orf_dict[key].exp_mean for key in orf_dict.keys()]
        
        #Plot expression histogram
        fig1, ax1 = plt.subplots()
        ax1.hist(exp_mean, bins = 20)  
        ax1.set_title('Histogram of ' + org_abbrev + ' mRNA expression values')
        ax1.set_xlabel('Mean Exp mRNA/gDNA')
        ax1.set_ylabel('Count')    
        plt.show()
                                                         
        #Add data from microarray experiments to gene dictionary
        fname_growth = os.path.normpath(data_dir + '/' + 'GSE36253_family.soft') 
               
        exp_conds = read_micro_data(orf_dict, fname_growth, gene_dict,organism,orf_header_name,org_soft_alt_dict,'Growth')
        
        #The other dataset has h2O2 and heat stress treatments.
        fname_heat_ox = os.path.normpath(data_dir + '/' + 'GSE38478_family.soft')
        #exp_conds_heat_ox = read_micro_data_heat_ox(orf_dict, fname_heat_ox, gene_dict)
        exp_conds_heat_ox = read_micro_data(orf_dict, fname_heat_ox, gene_dict,organism,orf_header_name,org_soft_alt_dict,'HeatOx')
           
        exp_conds = exp_conds + exp_conds_heat_ox
        cond_order = {'LAG': 'A', 'LL': 'B', 'DS': 'C', 'PS': 'D','PLAT':'E', 'H2O2 T5':'F','H2O2 T15':'G','H2O2 T30': 'H','H2O2 T60':'I','HS T5' : 'J', 'HS T15' : 'K', 'HS T30' : 'L', 'HS T45' : 'M' , 'HS T60' : 'N' , 'NaCl T5' : 'O', 'NaCl T15': 'P', 'NaCl T30' : 'Q', 'NaCl T60' : 'R'}
        #Sort experimental conditions
        exp_conds_sorted = sorted(exp_conds, key = lambda cond: cond_order[cond.split(',')[0]])
        #Sort genes by top expression
        top_orfs = [orf_dict[key] for key in orf_dict.keys() if orf_dict[key].exp_mean > 2.8] 
        top_orfs_sorted = sorted(top_orfs, key = lambda orf: orf.exp_mean)
        
        
        #for orf in top_orfs_sorted:
        #    orf.cond_calcs(exp_conds_sorted)
                        
        column_labels = [cond.split(',')[0] for cond in exp_conds_sorted]
        row_labels = [orf.orf for orf in top_orfs_sorted]
        
        #Heatmap of genes with high expression over different experimental conditions                
        fig2, ax2 = heatmap_plot(top_orfs_sorted, exp_conds_sorted, row_labels, column_labels)
                    
        #Replace data value with average value from all replicates.  
        orf_dict_mean = OrfDict()
        #Tracer()()
        for key in orf_dict.keys():
            orf = orf_dict[key]
            orf_mean = orf_cond_mean(orf,exp_conds_sorted)[0]
            orf_dict_mean[key] = orf_mean
        
        #plot average values for all replicates        
        top_orfs_sorted_mean = [orf_dict_mean[orf.orf] for orf in top_orfs_sorted]
        exp_conds_sorted_mean = orf_cond_mean(top_orfs_sorted[0],exp_conds_sorted)[1]
               
        row_labels_mean = [orf.orf for orf in top_orfs_sorted_mean]
        fig3, ax3 = heatmap_plot(top_orfs_sorted_mean, exp_conds_sorted_mean, row_labels_mean, exp_conds_sorted_mean)
        
        #Calculate variability for each orf: 
        #Calculate average values for each summary condition (to use as NA penalties in variability calculation)
        cond_means = []
        for cond in exp_conds_sorted_mean: 
            cond_all_vals = [orf_dict_mean[key].conditions[cond] for key in orf_dict_mean.keys()]
            cond_means.append(mean_NAs(cond_all_vals))
        
        #Set weights for each mean condition: ['LAG', 'LL', 'DS', 'PS', 'PLAT']
        #SC includes NaCl
        cond_weights = [1.0, 1.0, 1.0, 1.0, 1.0,0.25,0.25,0.25,0.25,0.2,0.2,0.2,0.2,0.2,0.25,0.25,0.25,0.25]
                          
        #Calculate variability score for all orfs
        var_vec = []
        for key in orf_dict_mean.keys():
            orf = orf_dict_mean[key]
            orf.var_calcs(exp_conds_sorted_mean , cond_weights, cond_means)
            var_vec.append(orf.variability)
            
        #Plot histogram of variability
        fig4, ax4 = plt.subplots()
        ax4.hist(var_vec, bins = 20)  
        ax4.set_title('Histogram of KL expression variability scores')
        ax4.set_xlabel('Variability')
        ax4.set_ylabel('Count')    
        plt.show()
        
        #Pick out top 1000 variable genes and sort by expression level
        orfs_by_var = sorted(orf_dict_mean.keys(), key = lambda orf_name: orf_dict_mean[orf_name].variability)
        top_orfs_by_var_names = orfs_by_var[0:2000]    
        top_orfs_by_var = [orf_dict_mean[key] for key in top_orfs_by_var_names]
        
        top_orfs_by_var_exp = sorted(top_orfs_by_var, key = lambda orf: orf.exp_mean, reverse = True)
        top_orfs_by_var_exp = top_orfs_by_var_exp[0:20]
        top_orfs_by_var_exp_names = [orf.orf for orf in top_orfs_by_var_exp]
                        
        fig5, ax5 = heatmap_plot(top_orfs_by_var_exp, exp_conds_sorted_mean, top_orfs_by_var_exp_names, exp_conds_sorted_mean)
                                    
        #Pick a high, medium and low promoter

        print top_orfs_by_var_exp_names
        
        
        #Plot of genes of interest over different experimental conditions.   
        #goi_names = ['ADH1', 'VPS4', 'NOP7', 'CBF5', 'TEF1','TDH3', 'Strong1 ~SC.RPS21A', 'Strong2 ~SC.RPS25A', 'Strong3 ~RPL37A/B','LowVar_high ~GCN4','LowVar_med1 ~BGL2','LowVar_med2 ~ COX8','LowVar_low1 ~COX6','LowVarLow2 ~GLC7','LowVarLow3 ~YNL0190W' ,'KLLA0D17864g', 'KLLA0A01584g', 'KLLA0A08580g']
        #goi_list = ['KLLA0F21010g','KLLA0B10846g','KLLA0D15268g','KLLA0D04796g','KLLA0B08998g','KLLA0A11858g','KLLA0A05700g','KLLA0B06182g','KLLA0C01870g','KLLA0D14113g','KLLA0F03036g','KLLA0B03311g','KLLA0F21890g','KLLA0F12496g','KLLA0D08624g', 'KLLA0D17864g', 'KLLA0A01584g', 'KLLA0A08580g']
        #goi_names = [ 'KL.MCM1','KL.RAP1','KL.RPS24','KL.RPS25','KLLA0D16027g', 'KLLA0F16511g', 'KLLA0D06941g', 'KLLA0B04686g']
        #goi_list = ['KLLA0F26807g','KLLA0D19294g','KLLA0C07755g','KLLA0B06182g','KLLA0D16027g', 'KLLA0F16511g', 'KLLA0D06941g', 'KLLA0B04686g']
             
        #goi_names = ['ADH1','NOP7','RPS25A','RPS25B','MSN2','TPK1']
        #goi_list = ['YOL086C', 'YGR103W','YGR027C','YLR333C','YMR037C','YJL164C']
        
        #SC PKA Pathway
        #goi_names = ['TPK1','TPK2','TPK3', 'BCY1','RAS2','RSR1','RAS1','CYR1','PDE1','PDE2','CDC25','BUD5','IRA2','IRA1','BUD2']
        #goi_list = ['YJL164C', 'YPL203W', 'YKL166C', 'YIL033C', 'YNL098C','YGR152C','YOR101W','YJL005W','YGL248W','YOR360C','YLR310C','YCR038C','YOL081W','YBR140C','YKL092C']
        
        #KL PKA Pathway
        #goi_names = ['TPK1','TPK2','BCY1','RAS2','RSR1','CYR1','PDE1','PDE2','CDC25','BUD5','IRA2','BUD2']
        #goi_list = ['KLLA0D03190g','KLLA0B07205g','KLLA0E04070g','KLLA0C13387g','KLLA0C12001g','KLLA0F15708g','KLLA0D11528g','KLLA0A03619g','KLLA0D09306g','KLLA0C03410g','KLLA0E17963g','KLLA0B06457g']

        
        #SC TFs
        #goi_names = ['MSN2','MSN4','DOT6','TOD6','STB3','CRZ1']
        #goi_list = ['YMR037C','YKL062W','YER088C','YBL054W','YDR169C','YNL027W']
        
        #KL TFs
        goi_names = ['MSN2','DOT6','STB3','CRZ1']
        goi_list = ['KLLA0F26961g','KLLA0E20031g','KLLA0D17270g','KLLA0E08679g']
        
        
        #goi_names = ['CDC37','SSK2','SSK22','STE11','STE50']
        #goi_list = []
        
        
        
        goi_orf_list = [orf_dict[key] for key in goi_list]
        
        #TDH3 - 'YGR192C',
        
        fig6, ax6 = heatmap_plot(goi_orf_list, exp_conds_sorted, goi_names, column_labels)
        
        goi_orf_list_mean = [orf_dict_mean[key] for key in goi_list]        
        fig7, ax7 = heatmap_plot(goi_orf_list_mean, exp_conds_sorted_mean, goi_names, exp_conds_sorted_mean)

        
        
if __name__=="__main__":
	main()