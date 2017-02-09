# -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
import matplotlib as mpl
mpl.use('Qt4Agg') # Need this (at least on Mac OSX to get window to display properly)
# import regev_library
# -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# 
# The point of this file was to have a constant parameter file where I could change values that
# show up everywhere in the code from one spot. Seemed to work pretty well overall. Changing
# the species_list here, for example, should change the compiled data and plots across the board.
# When writing new code, you'll want to reference these variables too in order to keep things
# consistent. 
# 
# -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# species_list = ['Saccharomyces cerevisiae', 'Kluyveromyces lactis']


species_list = ['Saccharomyces cerevisiae', 'Saccharomyces bayanus', 
				'Candida glabrata', 'Saccharomyces castellii',
				'Kluyveromyces lactis']



# Refernce dictionary for creating file names
species_name_dict = {'Saccharomyces cerevisiae' : 'SCer',
                    'Kluyveromyces lactis': 'KLac', 
                    'Candida glabrata' : 'CGla', 
                    'Saccharomyces castellii' : 'SCas', 
                    'Saccharomyces bayanus' : 'SBay'}


species_abb_dict = { 'SCer' : 'Saccharomyces cerevisiae',
                     'KLac' : 'Kluyveromyces lactis', 
                     'CGla' : 'Candida glabrata', 
                     'SCas' : 'Saccharomyces castellii', 
                     'SBay' : 'Saccharomyces bayanus'}


# -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------

def create_condition_key(species_arg, condition):
    species = ''
    if species_arg in species_abb_dict:
        species = species_abb_dict[species_arg]
    elif species in species_list:
        species = species_arg
    else:
        raise KeyError('Invalid Species Name : ' + species_arg)

    condition_key = species + ':' + condition

    return condition_key


def get_species(condition_key):
    return condition_key.split(':')[0]

def get_condition(condition_key):
    return condition_key.split(':')[1]

# -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------

HIGH_NUM = 15
LOW_NUM = 15
condition_key = create_condition_key('SCer', 'susan')
condition_ext = species_name_dict[get_species(condition_key)] + '-' + get_condition(condition_key).replace('/', '-')
# condition_key_ext = condition_key.replace(':', '-')
cross_ext = ''
for index, species in enumerate(species_list):
	cross_ext = cross_ext + species_name_dict[species]
	if index != len(species_list) - 1:
		cross_ext = cross_ext + '_'

ext = str(HIGH_NUM) + '_' + str(LOW_NUM) + '_' + condition_ext + '_' + cross_ext

# -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------


# Parameters for choosing up and down regulated genes in S.Cer pka inhibition
# POSITIVE_LOG_CUTOFF = 6.5
# NEGATIVE_LOG_CUTOFF = -6.0


# Create ortholog lookup tables for each species
# species_list = ['Kluyveromyces lactis', 'Saccharomyces castellii', 
#                  'Candida glabrata', 'Saccharomyces bayanus', 'Saccharomyces cerevisiae']


# species_list = ['Saccharomyces cerevisiae', 'Saccharomyces bayanus', 
# 				'Candida glabrata', 'Saccharomyces castellii',
# 				'Kluyveromyces lactis']