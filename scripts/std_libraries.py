#to keep my scripts consistent, I am adding this to the top of all scripts with %load std_libraries.py
import sys

#Indicate operating environment and import core modules
location_input = input("what computer are you on? a = Ben's laptop, b = gpucluster, c = Ben's desktop, d = other")
location_dict = {'a': "C:\\Users\\BMH_work\\github\\expression_broad_data", 'b': "/home/heineike/github/expression_broad_data",
                 'c': "C:\\Users\\Ben\\Documents\\GitHub\\expression_broad_data", 'd':'you need to add your location to the location_dict'}
figsave_dict = {'a': "C:\\Users\\BMH_work\\Google Drive\\UCSF\\ElSamad_Lab\\PKA\\Manuscript\\" , 
                'b': "/home/heineike/scratch/",
                'c': "C:\\Users\\Ben\\Google Drive\\UCSF\\ElSamad_Lab\\PKA\\Manuscript\\", 
                'd': 'you need to add your location to the figsave dict'}

base_dir = location_dict[location_input]
figsave_dir = figsave_dict[location_input]

print("base directory is " + base_dir)

if sys.path[-1] != base_dir:
    sys.path.append(base_dir)
    print("Added " + base_dir + " to path: " )
    print(sys.path)

import os

print("Importing io_library and setting base_dir and data_processing_dir")
#from core import expression_plots 
import io_library 
io_library.base_dir = base_dir
data_processing_dir = base_dir + os.sep + os.path.normpath("expression_data") + os.sep
io_library.data_processing_dir = data_processing_dir

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Patch, Rectangle  #Circle, Wedge, 
from matplotlib.collections import PatchCollection
from matplotlib.gridspec import GridSpec
#import matplotlib.colormap as cm
from matplotlib_venn import venn2
from matplotlib_venn import venn3
import squarify  #Makes treegraphs. 
#for my windows10 laptop I had to install this package using pip rather than anaconda.  
import seaborn as sns; sns.set(style="ticks", color_codes=True)  #not sure why I set color codes on ticks to be true
#from sklearn import linear_model
import pickle
import subprocess
#import networkx as nx
import scipy.stats as stats
import statsmodels.api as sm
import scipy.spatial.distance as spd
#import statsmodels.graphics.gofplots as stats_graph
import scipy.cluster.hierarchy as sch
import FisherExact


from Bio import SeqIO
from Bio import pairwise2
# from Bio import SeqFeature as sf
# from Bio.Alphabet import generic_dna
# from Bio.Seq import Seq
from Bio import Phylo
import gffutils

# import re

from collections import Counter, OrderedDict
# import scipy.stats as stats
from itertools import chain
#from itertools import product
#this only works if you are online
online_input = input("are you online? Yes/No ")
if online_input == "Yes": 
    import plotly.plotly as py
    import plotly.graph_objs as pygo
    import plotly.tools as pytools
    import plotly.io as pio
    py.sign_in('heineike02_student','9dMTMZgJMgUP0YX0P5mQ')
    #py.sign_in('heineike02', 'APPjKrtARaN2ZgUYIkqr')
    
# # for phylogenetic trees: 
from ete3 import Tree, SeqMotifFace, TreeStyle, add_face_to_node  #the last three are for visualization
# # In order to view ete3 created trees on the gpucluster, you need to use a virtual X server:
# from pyvirtualdisplay import Display
# display = Display(visible=False, size=(1024, 768), color_depth=24)
# display.start()
#ete3 is not officially supported on windows, and so must be loaded via pip: 
# pip install -U https://github.com/etetoolkit/ete/archive/qt5.zip
# ref: https://groups.google.com/forum/#!topic/etetoolkit/6NblSBPij4o

#for scraping internet data (e.g. ncbi)
import requests
from bs4 import BeautifulSoup
#from lxml import etree    #parses xml output

