import os
base_dir = 'C:\Users\Ben\Documents\GitHub\expression_broad_data'
os.chdir(base_dir) 
from core import expression_plots 
from core import io_library 
from IPython.core.debugger import Tracer
#import numpy as np
#import pandas as pd
#import re
#import matplotlib.pyplot as plt 
#import seaborn as sns

def main():
    ESR_genes = io_library.read_gasch_data('gasch_fig3_ESR.txt')
    print(ESR_genes)
    
if __name__=="__main__":
    main()