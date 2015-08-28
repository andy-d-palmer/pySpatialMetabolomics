import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.image as mpimg
import csv
import sys
import os
sys.path.append('/Users/palmer/Documents/python_codebase/')
from pyIMS.hdf5.inMemoryIMS_hdf5 import inMemoryIMS_hdf5
from pyMS.pyisocalc import pyisocalc

def get_variables(json_filename):
    import json
    config = json.loads(open(json_filename).read())
    return config

def generate_isotope_patterns(config):
    from pySpatialMetabolomics.parse_databases import parse_databases
    import pickle
    ### We simulate a mass spectrum for each sum formula/adduct combination. This generates a set of isotope patterns (see http://www.mi.fu-berlin.de/wiki/pub/ABI/QuantProtP4/isotope-distribution.pdf) which can provide additional informaiton on the molecule detected. This gives us a list of m/z centres for the molecule
    def calcualte_isotope_patterns(sum_formulae,adducts='',isocalc_sig=0.01,isocalc_resolution = 200000.,isocalc_do_centroid = True, charge='1'):
        ### Generate a mz list of peak centroids for each sum formula with the given adduct
        # todo - parse sum formula and adduct properly so that subtractions (losses) can be utilised (this code already exists somewhere)
        mz_list={}
        for n,sum_formula in enumerate(sum_formulae):   
            isotope_ms = pyisocalc.isodist(sum_formula+adduct,plot=False,sigma=isocalc_sig,charges=charge,resolution=isocalc_resolution,do_centroid=isocalc_do_centroid)
            if not sum_formula in mz_list:
                 mz_list[sum_formula] = {}
            mz_list[sum_formula][adduct] = isotope_ms.get_spectrum(source='centroids')
        return mz_list
    # Extract variables from config dict
    db_filename = config['file_inputs']['database_file']
    db_dump_folder = config['file_inputs']['database_load_folder']  
    isocalc_sig = config['isotope_generation']['isocalc_sig']  
    isocalc_resolution = config['isotope_generation']['isocalc_resolution']  
    if len(config['isotope_generation']['charge']) > 1:
        print 'Warning: only first charge state currently accepted'
    charge = int('{}{}'.format(config['isotope_generation']['charge'][0]['polarity'], config['isotope_generation']['charge'][0]['n_charges'])) #currently only supports first charge!!
    adducts=[a['adduct'] for a in config['isotope_generation']['adducts']]
  
    # Read in molecules
    sum_formulae = parse_databases.read_generic_csv(db_filename) 
    # Check if already genrated and load if possible, otherwise calculate fresh   
    db_name =  os.path.splitext(os.path.basename(db_filename))[0] 
    mz_list={}
    for adduct in adducts:
        load_file = '{}/{}_{}_{}_{}.dbasedump'.format(db_dump_folder,db_name,adduct,isocalc_sig,isocalc_resolution)
        if os.path.isfile(load_file):
            mz_list_tmp = pickle.load(open(load_file,'r'))
        else:
            mz_list_tmp = calcualte_isotope_patterns(sum_formulae,adducts=(adduct,),isocalc_sig=isocalc_sig,isocalc_resolution=isocalc_resolution,charge=charge)
            if db_dump_folder != "":
                pickle.dump(mz_list_tmp,open(load_file,'w'))
        # add patterns to total list
        for sum_formula in mz_list_tmp:
            if not sum_formula in mz_list:
                mz_list[sum_formula] = {}
            mz_list[sum_formula][adduct] = mz_list_tmp[sum_formula][adduct]
    print 'all isotope patterns generated and loaded'
    return sum_formulae,adducts,mz_list
   
def hot_spot_removal(xics,q):
    for xic in xics:
        xic_q = np.percentile(xic,q)
        xic[xic>xic_q]=xic_q
    return xics
def run_search(config):
    ### Runs the main pipeline
    # Get sum formula and predicted m/z peaks for molecules in database
    sum_formulae,adducts,mz_list = generate_isotope_patterns(config)
    # Parse dataset
    from pyIMS.hdf5.inMemoryIMS_hdf5 import inMemoryIMS_hdf5
    from pyIMS.image_measures import level_sets_measure,isotope_image_correlation,isotope_pattern_match
    IMS_dataset=inMemoryIMS_hdf5(config['file_inputs']['data_file'])
        
    ppm = config['image_generation']['ppm'] #parts per million -  a measure of how accuracte the mass spectrometer is
    nlevels = config['image_generation']['nlevels']  # parameter for measure of chaos
    q=config['image_generation']['q'] 
    measure_value_score={}
    iso_correlation_score = {}
    iso_ratio_score = {}
    
    for sum_formula in sum_formulae:
        for adduct in adducts:
            # 1. Geneate ion images
            ion_datacube = IMS_dataset.get_ion_image(mz_list[sum_formula][adduct][0],ppm) #for each spectrum, sum the intensity of all peaks within tol of mz_list
            ion_datacube.xic=hot_spot_removal(ion_datacube.xic,q)

            # 2. Spatial Chaos 
            if not sum_formula in measure_value_score:
                measure_value_score[sum_formula] = {}
            measure_value_score[sum_formula][adduct] = 1-level_sets_measure.measure_of_chaos(ion_datacube.xic_to_image(0),nlevels,interp=False)[0]
            if measure_value_score[sum_formula][adduct] == 1:
                measure_value_score[sum_formula][adduct] = 0
            # only compare pixels with values in the monoisotopic (otherwise high correlation for large empty areas)
            notnull_monoiso = ion_datacube.xic[0] > 0 
            #for xic in ion_datacube.xic:
            #    xic = xic[notnull_monoiso]
            # 3. Score correlation with monoiso
            if not sum_formula in iso_correlation_score:
                iso_correlation_score[sum_formula] = {}
            if len(mz_list[sum_formula][adduct][1]) > 1:
                iso_correlation_score[sum_formula][adduct] = isotope_image_correlation.isotope_image_correlation(ion_datacube.xic,weights=mz_list[sum_formula][adduct][1][1:])
            else: # only one isotope peak, so correlation doesn't make sense
                iso_correlation_score[sum_formula][adduct] = 1
        
            # 4. Score isotope ratio
            if not sum_formula in iso_ratio_score:
                iso_ratio_score[sum_formula] = {}
            iso_ratio_score[sum_formula][adduct]  = isotope_pattern_match.isotope_pattern_match(ion_datacube.xic,mz_list[sum_formula][adduct][1])
    return measure_value_score,iso_correlation_score,iso_ratio_score
        
def score_results(config,measure_value_score, iso_correlation_score, iso_ratio_score):
    def check_pass(pass_thresh,pass_val):
        tf = []
        for v,t in zip(pass_val,pass_thresh):
            tf.append(v>t)
        if all(tf):
            return True
        else:
            return False
    measure_tol = config['results_thresholds']['measure_tol']
    iso_corr_tol = config['results_thresholds']['iso_corr_tol']
    iso_ratio_tol = config['results_thresholds']['iso_ratio_tol']
    sum_formulae,adducts,mz_list = generate_isotope_patterns(config)
    pass_formula = []
    for sum_formula in sum_formulae:
        for adduct in adducts:
            if check_pass((measure_tol,iso_corr_tol,iso_ratio_tol),(measure_value_score[sum_formula][adduct],iso_correlation_score[sum_formula][adduct],iso_ratio_score[sum_formula][adduct])):
                    pass_formula.append('{} {}'.format(sum_formula,adduct))
    return pass_formula

def output_results(config,measure_value_score,iso_correlation_score,iso_ratio_score):
    import os
    filename_in = config['file_inputs']['data_file']
    output_dir = config['file_inputs']['results_folder']
    measure_tol = config['results_thresholds']['measure_tol']
    iso_corr_tol = config['results_thresholds']['iso_corr_tol']
    iso_ratio_tol = config['results_thresholds']['iso_ratio_tol']
    sum_formulae,adducts,mz_list = generate_isotope_patterns(config)     
    # Save the processing results
    if os.path.isdir(output_dir)==False:
        os.mkdir(output_dir)
    filename_out = '{}{}{}_full_results.txt'.format(output_dir,os.sep,os.path.splitext(os.path.basename(filename_in))[0])
    with open(filename_out,'w') as f_out:
        f_out.write('sf,adduct,mz,moc,spec,spat,pass\n'.format())
        for sum_formula in sum_formulae:
            for adduct in adducts:
                moc_pass =  measure_value_score[sum_formula][adduct] > measure_tol and iso_correlation_score[sum_formula][adduct] > iso_corr_tol and iso_ratio_score[sum_formula][adduct] > iso_ratio_tol
                f_out.write('{},{},{},{},{},{},{}\n'.format(
                        sum_formula,
                        adduct,
                        mz_list[sum_formula][adduct][0][0],
                        measure_value_score[sum_formula][adduct],
                        iso_correlation_score[sum_formula][adduct],
                        iso_ratio_score[sum_formula][adduct],
                        moc_pass)) 

    filename_out = '{}{}{}_pass_results.txt'.format(output_dir,os.sep,os.path.splitext(os.path.basename(filename_in))[0])
    with open(filename_out,'w') as f_out:
        f_out.write('ID,sf,adduct,mz,moc,spec,spat\n'.format())
        for sum_formula in sum_formulae:
            for adduct in adducts:
                if measure_value_score[sum_formula][adduct] > measure_tol and iso_correlation_score[sum_formula][adduct] > iso_corr_tol and iso_ratio_score[sum_formula][adduct] > iso_ratio_tol:
                    f_out.write('{},{},{},{},{},{},{}\n'.format(
                        sum_formulae[sum_formula]['db_id'],
                        sum_formula,adduct,
                        mz_list[sum_formula][adduct][0][0],
                        measure_value_score[sum_formula][adduct],
                        iso_correlation_score[sum_formula][adduct],
                        iso_ratio_score[sum_formula][adduct]))


def run_pipeline(JSON_config_file):
    config = get_variables(JSON_config_file)
    measure_value_score, iso_correlation_score, iso_ratio_score = run_search(config)
    pass_list = score_results(config,measure_value_score, iso_correlation_score, iso_ratio_score)
    output_results(config,measure_value_score,iso_correlation_score,iso_ratio_score)


