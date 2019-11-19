# Yeast Environmental Stress Response (ESR) analysis tools

This code provides the scripts used to analyse gene expression data related to PKA inhibition and Environmental Stress Response in budding yeast species for in [1]

## Getting Started

The main functions are contained in the file yeast_esr_exp.py

Scripts to generate analysis and figures are in the scripts folder
The scripts used for ref [1] are listed in PKA_evolution_fig_notes.txt
Also in the scripts folder is notebook_setup.py which you can load at the beginning of a script with: 

%load notebook_setup.py

It runs the std_libraries.py file which loads all necessary libraries for the provided scripts. 


Data should be in the folder expression_data which is tracked using git-lfs.  The folder itself can be obtained from figshare [???]


## References

[1] "Paralogs in the PKA regulon traveled different evolutionary routes to divergent expression in budding yeast." by Ben Heineike and Hana El-Samad.


## Authors

* **Benjamin Heineike** [heineike02](https://github.com/heineike02)


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
