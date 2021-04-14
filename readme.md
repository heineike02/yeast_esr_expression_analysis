# Yeast Environmental Stress Response (ESR) analysis tools

This code provides the scripts used to analyse gene expression data from [1].  This includes analysis of RNA-seq data under PKA inhibition gathered for the paper as well as analysis of gene expression data across multiple species from [2],[3], and [4].  

## Getting Started

The main functions are contained in the file yeast_esr_exp.py

Scripts to generate analysis and figures are in the scripts folder
The scripts used for ref [1] are listed in PKA_evolution_fig_notes.txt
Also in the scripts folder is notebook_setup.py which you can load at the beginning of a script with: 

%load notebook_setup.py

It runs the std_libraries.py file which loads all necessary libraries for the provided scripts. 

Data should be in the folder expression_data which can be obtained from figshare.  

The code was primarily run on a windows environment using conda (see environment.yml) but should work with Linux as well. 

## References

[1] Heineike, B., and El-Samad, H. (2021). Paralogs in the PKA regulon traveled different evolutionary routes to divergent expression in budding yeast. Front. Fungal Biol. 2. doi:10.3389/ffunb.2021.642336.

[2] Thompson, D. A., Roy, S., Chan, M., Styczynsky, M. P., Pfiffner, J., French, C., et al. (2013). Evolutionary principles of modular gene regulation in yeasts. Elife 2, e00603. doi:10.7554/eLife.00603.

[3] Roy, S., Wapinski, I., Pfiffner, J., French, C., Socha, A., Konieczka, J., et al. (2013). Arboretum: reconstruction and analysis of the evolutionary history of condition-specific transcriptional modules. Genome Res. 23, 1039–1050. doi:10.1101/gr.146233.112.

[4] Tsankov, A. M., Thompson, D. A., Socha, A., Regev, A., and Rando, O. J. (2010). The role of nucleosome positioning in the evolution of gene regulation. PLoS Biol. 8, e1000414. doi:10.1371/journal.pbio.1000414.


## Authors

* **Benjamin Heineike** [heineike02](https://github.com/heineike02)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
