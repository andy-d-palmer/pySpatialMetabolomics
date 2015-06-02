
def main_pipeline(sum_formulae,adducts,):
    ppm = 3.; #parts per million -  a measure of how accuracte the mass spectrometer is
    nlevels = 30 # parameter for measure of chaos
    q=99
    tstart=time.time()
    measure_value_score={}
    iso_correlation_score = {}
    iso_ratio_score = {}

    for n,sum_formula in enumerate(sum_formulae):
        for adduct in adducts:
            # 1. Geneate ion images
            #mz_list[sum_formula][adduct][0] #get centroid mz values
            #tol = mz_list*ppm/1e6 # turn ppm accuracy into a mass range
            ion_datacube = IMS_dataset.get_ion_image(mz_list[sum_formula][adduct][0],ppm) #for each spectrum, sum the intensity of all peaks within tol of mz_list
            for xic in ion_datacube.xic:
                xic_q = np.percentile(xic,q)
                xic[xic>xic_q]=xic_q

            # 2. Spatial Chaos 
            if not sum_formula in measure_value_score:
                measure_value_score[sum_formula] = {}

            if np.sum(ion_datacube.xic_to_image(0)) == 0:
                measure_value_score[sum_formula][adduct] = 0 # this is now coded into measure_of_chaos
            else:
                measure_value_score[sum_formula][adduct] = 1-level_sets_measure.measure_of_chaos(ion_datacube.xic_to_image(0),nlevels,interp=False)[0]
                if measure_value_score[sum_formula][adduct] == 1:
                     measure_value_score[sum_formula][adduct] = 0
                clear_output(wait=True)

            # 3. Score correlation with monoiso
            if not sum_formula in iso_correlation_score:
                iso_correlation_score[sum_formula] = {}
            if len(mz_list[sum_formula][adduct][1]) > 1:
                ## select just positive pixels
                   # iso_correlation_score[sum_formula][adduct] = np.average(
                   #     [np.corrcoef(ion_datacube.xic[0][notnull_monoiso],ion_datacube.xic[ii][notnull_monoiso])[0,1] 
                   #          for ii in range(1,len(mz_list[sum_formula][adduct][1]))] 
                   #             ,weights=mz_list[sum_formula][adduct][1][1:])
                ## use all pixels (i.e. tissuemask - mask not coded so not editable :S)

                #iso_correlation = []
                #for ii in enumerate(mz_list[sum_formula][adduct][1])
                #   ios_correlation.append(np.corrcoef(ion_datacube.xic)[0],np.corrcoef(ion_datacube.xic)[ii])
                iso_correlation = np.corrcoef(ion_datacube.xic)[1:,0]
                iso_correlation[np.isnan(iso_correlation)] = 0 # when alll values are the same (e.g. zeros) then correlation is undefined
                iso_correlation_score[sum_formula][adduct] = np.average(
                        iso_correlation,weights=mz_list[sum_formula][adduct][1][1:]
                         ) # slightly faster to compute all correlations and pull the elements needed
            else: # only one isotope peak, so correlation doesn't make sense
                iso_correlation_score[sum_formula][adduct] = 1

            # 4. Score isotope ratio
            if not sum_formula in iso_ratio_score:
                iso_ratio_score[sum_formula] = {}
            isotope_intensity = mz_list[sum_formula][adduct][1] 
            image_intensities = [sum(ion_datacube.xic[ii]) for ii in range(0,len(mz_list[sum_formula][adduct][0]))]
            iso_ratio_score[sum_formula][adduct] = 1-np.mean(abs( isotope_intensity/np.linalg.norm(isotope_intensity) - image_intensities/np.linalg.norm(image_intensities)))    
        if np.mod(n,10)==0:
            clear_output(wait=True)
            sys.stdout.flush()
            print '{} {:2.3f}\% complete\r'.format(sum_formula,100*float(n)/len(sum_formulae),end="\r")
    print 'Elapsed: {:5.2f} seconds'.format(time.time() - tstart)

