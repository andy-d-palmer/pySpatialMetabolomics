__author__ = 'palmer'

def select_passes(score_data,t_holds):
    import numpy as np
    moc_t,spec_t,spat_t = t_holds
    pass_tf = np.all([score_data['moc'] > moc_t, score_data['spat'] > spat_t, score_data['spec'] > spec_t],axis=0)
    pass_data = score_data[pass_tf]
    return pass_data

def select_passes_linear(score_data,t_hold):
    import numpy as np
    pass_tf = np.asarray(score_data['moc'] + score_data['spat'] + score_data['spec']) > t_hold
    pass_data = score_data[pass_tf]
    return pass_data

def select_passes_mult(score_data,t_hold):
    import numpy as np
    pass_tf = np.asarray(score_data['moc'] * score_data['spat'] * score_data['spec']) > t_hold
    pass_data = score_data[pass_tf]
    return pass_data

def count_adducts(pass_data,adducts):
    ## number of false hits by adduct
    fd_count = pass_data['adduct'].value_counts() # count adduct occurances
    set_diff = set(adducts) - set(fd_count.keys())
    for s in set_diff: #add zero occurance back in
        fd_count[s] = 0
    return fd_count.loc[[f in adducts for f in fd_count.keys()]]

def calc_fdr_adducts(fd_count,adducts,plausible_adducts):
    import numpy as np
    implausible_adducts = [a for a in adducts if a not in plausible_adducts]
    implausible_list = fd_count[implausible_adducts]
    implausible_hits = sum(implausible_list)
    mean_implausible = float(implausible_hits)/len(implausible_list)
    mean_implausible = np.percentile(implausible_list,50)
    plausible_hits = np.sum(fd_count[list(plausible_adducts)])
    if plausible_hits == 0:
        return 1,plausible_hits,implausible_hits
    fdr = (len(plausible_adducts)*mean_implausible) / plausible_hits
    return fdr,plausible_hits,implausible_hits
