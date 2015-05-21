
# coding: utf-8

# This notebook provides a demonstration of the spatial metabolomics pipeline applied to a single molecule from a database. During deployment the process is repeated for all 10K database entries. 
# 
# The pipeline can be broadly split into two stages:
# 
# 1. Predict spectral pattern of a molecule
# 2. Search data and score results
# 

# In[1]:

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import csv
get_ipython().magic(u'matplotlib inline')
import sys
from IPython.display import display, clear_output
import os
sys.path.append('/Users/palmer/Documents/python_codebase/')
from pyIMS.hdf5.inMemoryIMS_hdf5 import inMemoryIMS_hdf5
from pyMS.pyisocalc import pyisocalc
import time


## Predict Spatial Pattern

# An entry from the database is parsed and the sumformula is extracted. We take this and the database key so we can link back later.
# 
# We also define a set of 'adducts' these are molecules that can become attached during the mass spectrometry process and so change the mass of the molecule. It's possible that some or all of the adducts are detected in the data but there's no way of knowing a priori which ones.

# In[2]:

def read_kegg_compounds(filename):
    sum_formulae = {}
    with open(filename) as filein_db:
        data_in = csv.reader(filein_db)
        for line in data_in:
            if line == []: #catch empty line(s)
                continue
            db_id,name,sf,mw = line
            if float(mw) > 0: # useful filter for incomplete entries, not robust
                if sf in sum_formulae:
                    sum_formulae[sf]['name'].append(name)
                    sum_formulae[sf]['db_id'].append(db_id)
                else:
                    sum_formulae[sf] = {}
                    sum_formulae[sf]['name'] = [name]
                    sum_formulae[sf]['db_id'] = [db_id]
                    sum_formulae[sf]['mw'] = mw
        return sum_formulae
    
def read_hmdb_compounds(filename):
    sum_formulae = {}
    with open(filename) as filein_db:
        data_in = csv.reader(filein_db,delimiter='\t')
        for line in data_in:
            if line == []: #catch empty line(s)
                continue
            db_id,mw,sf,name = line
            if mw=='':
                mw = '0'
            if float(mw) > 0: # useful filter for incomplete entries, not robust
                if sf in sum_formulae:
                    sum_formulae[sf]['name'].append(name)
                    sum_formulae[sf]['db_id'].append(db_id)
                else:
                    sum_formulae[sf] = {}
                    sum_formulae[sf]['name'] = [name]
                    sum_formulae[sf]['db_id'] = [db_id]
                    sum_formulae[sf]['mw'] = mw
    return sum_formulae
            
#sum_formulae = read_kegg_compounds('/Users/palmer/Copy/kegg_compounds.csv')
sum_formulae = read_hmdb_compounds('/Users/palmer/Copy/hmdb_database.tsv')
#    sum_formula = 'C23H45NO4' #Glycocholic acid: http://89.223.39.196:2347/substance/00138
adducts = ('H','Na','K')
charge = 1

def generate_isotope_patterns(sum_formulae,save_dir):
	# We simulate a mass spectrum for each sum formula/adduct combination. This generates a set of isotope patterns (see http://www.mi.fu-berlin.de/wiki/pub/ABI/QuantProtP4/isotope-distribution.pdf) which can provide additional informaiton on the molecule detected. This gives us a list of m/z centres for the molecule

	# In[3]:

	mz_list={}
	for sum_formula in sum_formulae:
	    for adduct in adducts:
		isotope_ms = pyisocalc.isodist(sum_formula+adduct,plot=False,sigma=0.01,charges=charge,resolution=200000.0,do_centroid=True)
		if not sum_formula in mz_list:
		    mz_list[sum_formula] = {}
		mz_list[sum_formula][adduct] = isotope_ms.get_spectrum(source='centroids')
	print 'done'

def parse_input
	# Parse data
	from pyIMS.hdf5.inMemoryIMS_hdf5 import inMemoryIMS_hdf5
	from pyIMS.image_measures import level_sets_measure
	filename_in =  '/Users/palmer/Documents/tmp_data/14037s1_Spheroids24h_DHBSub_centroids_IMS.hdf5' #using a temporary hdf5 based format
	IMS_dataset=inMemoryIMS_hdf5(filename_in)

	## Search data and score results

	# We parse a dataset and then search it one adduct at a time using the following sequence:
	#     1. generate images for the peaks
	#     2. score this image using a measure of spatial chaos (checks if image is likley noise)
	#     3. score the distributions (should correlate with the monoisotopic)
	#     4. score the isotope ratio (average image intensity, should match predicted intensities)
	#     5. use scores to filter whether molecule is present 

	# In[4]:

	# Parse data
	from pyIMS.hdf5.inMemoryIMS_hdf5 import inMemoryIMS_hdf5
	from pyIMS.image_measures import level_sets_measure
	filename_in =  '/Users/palmer/Documents/tmp_data/14037s1_Spheroids24h_DHBSub_centroids_IMS.hdf5' #using a temporary hdf5 based format
	IMS_dataset=inMemoryIMS_hdf5(filename_in)


	# In[5]:
def run_search()
	ppm = 5.; #parts per million -  a measure of how accuracte the mass spectrometer is
	nlevels = 30 # parameter for measure of chaos
	q=99
	tstart=time.time()
	measure_value_score={}
	iso_correlation_score = {}
	iso_ratio_score = {}

	for sum_formula in sum_formulae:
	    for adduct in adducts:
		# 1. Geneate ion images
		mz_list[sum_formula][adduct][0] #get centroid mz values
		#tol = mz_list*ppm/1e6 # turn ppm accuracy into a mass range
		ion_datacube = IMS_dataset.get_ion_image(mz_list[sum_formula][adduct][0],ppm) #for each spectrum, sum the intensity of all peaks within tol of mz_list
		for xic in ion_datacube.xic:
		    xic_q = np.percentile(xic,q)
		    xic[xic>xic_q]=xic_q
		    xic[tissue_mask]==0
		# 2. Spatial Chaos 
		if not sum_formula in measure_value_score:
		    measure_value_score[sum_formula] = {}

		if np.sum(ion_datacube.xic_to_image(0)) == 0:
		    measure_value_score[sum_formula][adduct] = 0 # this is now coded into measure_of_chaos
		else:
		    measure_value_score[sum_formula][adduct] = 1-level_sets_measure.measure_of_chaos(ion_datacube.xic_to_image(0),nlevels,interp=True)[0]
		    if measure_value_score[sum_formula][adduct] == 1:
		         measure_value_score[sum_formula][adduct] = 0
		    clear_output(wait=True)
		    
		# 3. Score correlation with monoiso
		notnull_monoiso = ion_datacube.xic[0] > 0 # only compare pixels with values in the monoisotopic (otherwise high correlation for large empty areas)
		if not sum_formula in iso_correlation_score:
		    iso_correlation_score[sum_formula] = {}
		if len(mz_list[sum_formula][adduct][1]) > 1:
		    ## select just positive pixels
		       # iso_correlation_score[sum_formula][adduct] = np.average(
		       #     [np.corrcoef(ion_datacube.xic[0][notnull_monoiso],ion_datacube.xic[ii][notnull_monoiso])[0,1] 
		       #          for ii in range(1,len(mz_list[sum_formula][adduct][1]))] 
		       #             ,weights=mz_list[sum_formula][adduct][1][1:])
		    ## use all pixels (i.e. tissuemask - mask not coded so not editable :S)
		    iso_correlation_score[sum_formula][adduct] = np.average(
		            np.corrcoef(ion_datacube.xic)[1:,0],weights=mz_list[sum_formula][adduct][1][1:]
		             ) # slightly faster to compute all correlations and pull the elements needed
		else: # only one isotope peak, so correlation doesn't make sense
		    iso_correlation_score[sum_formula][adduct] = 1
	    
		# 4. Score isotope ratio
		if not sum_formula in iso_ratio_score:
		    iso_ratio_score[sum_formula] = {}
		isotope_intensity = mz_list[sum_formula][adduct][1] 
		image_intensities = [sum(ion_datacube.xic[ii]) for ii in range(0,len(mz_list[sum_formula][adduct][0]))]
		iso_ratio_score[sum_formula][adduct] = 1-np.mean(abs( isotope_intensity/np.linalg.norm(isotope_intensity) - image_intensities/np.linalg.norm(image_intensities)))    
	print 'Elapsed: {:5.2f} seconds'.format(time.time() - tstart)
		

def score_results()
	measure_tol = 0.99#heuristic tolerances for filter
	iso_corr_tol = 0.4
	iso_ratio_tol = 0.4
	pass_formula = []
	for sum_formula in sum_formulae:
	    for adduct in adducts:
		if measure_value_score[sum_formula][adduct] > measure_tol and iso_correlation_score[sum_formula][adduct] > iso_corr_tol and iso_ratio_score[sum_formula][adduct] > iso_ratio_tol:
		    pass_formula.append('{} {}'.format(sum_formula,adduct))

def output_results()
	# Save the processing results
	filename_out = '{}_results.txt'.format(os.path.dirname(filename_in))


def spatial_metabolomics(ini_file)
	get_variables()
	get_isotopes()
	run_search()
	score_results()
	output_results()


