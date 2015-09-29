# Option to read back in saved data
def load_results(filename_in):
    with open(filename_in,'r') as f_in:
        f_in.readline()
        measure_value_score={}
        iso_correlation_score={}
        iso_ratio_score={}
        moc_pass={}
        for line in f_in.readlines():
            (sum_formula, 
            adduct, 
            this_mw, 
            this_measure_value, 
            this_iso_correlation , 
            this_iso_ratio , 
            this_pass )= line=line.split(',')
            """adduct = line[1]
                    adduct,
                    measure_value_score[sum_formula][adduct],
                    iso_correlation_score[sum_formula][adduct],
                    iso_ratio_score[sum_formula][adduct],
                    moc_pass)) """
            if not sum_formula in measure_value_score:
                measure_value_score[sum_formula]={}
                iso_correlation_score[sum_formula]={}
                iso_ratio_score[sum_formula]={}
                moc_pass[sum_formula]={}                
            measure_value_score[sum_formula][adduct]=float(this_measure_value)
            iso_correlation_score[sum_formula][adduct]=float(this_iso_correlation)
            iso_ratio_score[sum_formula][adduct]=float(this_iso_ratio)
            moc_pass[sum_formula][adduct]=this_pass
            
        return (measure_value_score,
                    iso_correlation_score,
                    iso_ratio_score,
                    moc_pass)


def check_pass(pass_thresh,pass_val):
    tf = []
    for v,t in zip(pass_val,pass_thresh):
        tf.append(v>t)
    if all(tf):
        return True
    else:
        return False


def plot_images(ion_datacube,iso_spect,iso_max,ax=[],q_val=99):
    import matplotlib.pyplot as plt
    if ax==[]:
        fig = plt.figure(figsize=(20,15),dpi=300)
        ax = [   plt.subplot2grid((2, 4), (0, 0)),
             plt.subplot2grid((2, 4), (0, 1)),
             plt.subplot2grid((2, 4), (0, 2)),
             plt.subplot2grid((2, 4), (0, 3)),
             plt.subplot2grid((2, 4), (1, 0), colspan=4, rowspan=1)
         ]
    for a in ax:
        a.cla()
    # plot images
    for ii in range(0,iso_max):
        im = ion_datacube.xic_to_image(ii)
        # hot-spot removal
        notnull=im>0 
        if np.sum(notnull==False)==np.shape(im)[0]*np.shape(im)[1]:
            im=im
        else:
            im_q = np.percentile(im[notnull],q_val)
            im_rep =  im>im_q       
            im[im_rep] = im_q 

        ax[ii].imshow(im,cmap=c_map)
        ax[ii].set_title('m/z: {:3.4f}'.format(mz_list[sum_formula][adduct][0][ii]))
    # plot spectrum
    notnull=ion_datacube.xic_to_image(0)>0
    data_spect = [np.sum(ion_datacube.xic_to_image(ii)) for ii in range(0,iso_max)]
    data_spect = data_spect / np.linalg.norm(data_spect)
    iso_spect = iso_spect/np.linalg.norm(iso_spect)

    markerline, stemlines, baseline = ax[4].stem( mz_list[sum_formula][adduct][0][0:iso_max],iso_spect,'g')
    plt.setp(stemlines, linewidth=2, color='g')     # set stems  colors
    plt.setp(markerline, 'markerfacecolor', 'g','markeredgecolor','g')    # make points 

    markerline, stemlines, baseline = ax[4].stem( mz_list[sum_formula][adduct][0][0:iso_max],data_spect,'r')
    plt.setp(stemlines, linewidth=2, color='r')     # set stems colors
    plt.setp(markerline, 'markerfacecolor', 'r','markeredgecolor','r')    # make points 
    
    #plot proxy artist
    proxies=[]
    h, = plt.plot(mz_list[sum_formula][adduct][0][0],[0],'-g')
    proxies.append(h)
    h, = plt.plot(mz_list[sum_formula][adduct][0][0],[0],'-r')
    proxies.append(h)

    
    ax[4].legend(proxies,('predicted pattern','data pattern'), numpoints=1)
    return fig,ax
    plt.setp(markerline, 'markerfacecolor', 'g','markeredgecolor','g')    # make points

    markerline, stemlines, baseline = ax[4].stem( mz_list[0:iso_max],data_spect,'r')
    plt.setp(stemlines, linewidth=2, color='r')     # set stems colors
    plt.setp(markerline, 'markerfacecolor', 'r','markeredgecolor','r')    # make points

    #plot proxy artist
    proxies=[]
    h, = plt.plot(mz_list[0],[0],'-g')
    proxies.append(h)
    h, = plt.plot(mz_list[0],[0],'-r')
    proxies.append(h)


    ax[4].legend(proxies,('predicted pattern','data pattern'), numpoints=1)
    return fig,ax
