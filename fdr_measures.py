__author__ = 'palmer'

def calc_fdr_df(target_df,decoy_df,col='mult',ascending=False):
    import numpy as np
    import pandas as pd
    if len(target_df) != len(decoy_df):
        raise TypeError('target should be same length as decoy')
    score_vect = pd.concat((target_df,decoy_df),ignore_index=True)
    score_vect=score_vect.iloc[np.random.permutation(len(score_vect))]
    score_vect=score_vect.sort(columns=col,ascending=ascending)
    decoy_hits = np.cumsum(score_vect.index.values > len(target_df),dtype=float) #decoy hits have upper index after concat
    target_hits = np.cumsum(score_vect.index.values < len(target_df),dtype=float) #decoy hits have upper index after concat

    #decoy_hits = np.cumsum(score_tf,dtype=float)
    #target_hits = np.cumsum(score_tf==False,dtype=float)
    fdr_curve = decoy_hits/target_hits
    fdr_curve = np.asarray([d/t for d,t in zip(decoy_hits,target_hits)])
    fdr_curve[target_hits==0]=0
    return fdr_curve,target_hits,score_vect

def get_decoy_df(im_adducts,results_data,sf):
    import numpy as np
    import bisect
    import pandas as pd
    # decoy selection
    nDecoy=len(sf)
    decoy_adducts = np.random.choice(im_adducts,size=nDecoy,replace=True)
    dc_df=[]
    for s,da in zip(sf,decoy_adducts):
        dc_df.append(results_data[da].iloc[bisect.bisect_left(sf,s)])
        dc_df = pd.DataFrame.from_records(dc_df)
    return dc_df

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

def select_passes_l2(score_data,t_hold):
    import numpy as np
    pass_tf = np.sqrt(np.asarray(score_data['moc']**2 + score_data['spat']**2 + score_data['spec'])) > t_hold
    pass_data = score_data[pass_tf]
    return pass_data

def count_adducts(pass_data,adducts):
    ## number of false hits by adduct
    fd_count = pass_data['adduct'].value_counts() # count adduct occurances
    set_diff = set(adducts) - set(fd_count.keys())
    for s in set_diff: #add zero occurance back in
        fd_count[s] = 0
    return fd_count.loc[[f in adducts for f in fd_count.keys()]]

def calc_fdr_adducts(fd_count,adducts,plausible_adducts,average='mean'):
    import numpy as np
    implausible_adducts = [a for a in adducts if a not in plausible_adducts]
    implausible_list = fd_count[implausible_adducts]
    implausible_hits = sum(implausible_list)
    if average == "mean":
        mean_implausible = float(implausible_hits)/len(implausible_list)
    elif average == 'median':
        mean_implausible = np.percentile(implausible_list,50)
    else:
        raise ValueError('average type not recognised: {}'.format(average))
    plausible_hits = np.sum(fd_count[list(plausible_adducts)])
    if plausible_hits == 0:
        return 1,plausible_hits,implausible_hits
    fdr = (len(plausible_adducts)*mean_implausible) / plausible_hits
    return fdr,plausible_hits,implausible_hits
