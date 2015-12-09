__author__ = 'palmer'

def calc_fdr_df(target_df,decoy_df,col='mult',ascending=False):
    import numpy as np
    import pandas as pd
    if len(target_df) != len(decoy_df):
        raise TypeError('target should be same length as decoy {} {}'.format(len(target_df),len(decoy_df)))
    score_vect = pd.concat((target_df,decoy_df),ignore_index=True)
    score_vect=score_vect.iloc[np.random.permutation(len(score_vect))]
    score_vect=score_vect.sort(columns=col,ascending=ascending)
    decoy_hits = np.cumsum(score_vect.index.values > len(target_df),dtype=float) #decoy hits have upper index after concat
    target_hits = np.cumsum(score_vect.index.values < len(target_df),dtype=float) #decoy hits have upper index after concat

    #decoy_hits = np.cumsum(score_tf,dtype=float)
    #target_hits = np.cumsum(score_tf==False,dtype=float)
    #fdr_curve = decoy_hits/target_hits
    fdr_curve = np.asarray([d/t for d,t in zip(decoy_hits,target_hits)])
    fdr_curve[target_hits==0]=0
    return fdr_curve,target_hits,score_vect[col]

def get_decoy_df(decoy_adducts,results_data,sf):
    import numpy as np
    import bisect
    import pandas as pd
    # decoy selection
    nDecoy=len(sf)
    decoy_adducts = np.random.choice(decoy_adducts,size=nDecoy,replace=True)
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


def find_crossing(m,fdr_target):
    import numpy as np
    mm = np.sign(m[0:-1]-fdr_target) * np.sign(m[1:]-fdr_target) #looking for zero crossings
    mm_idx = np.where(mm < 0)
    if len(mm_idx[0]) == 0:
        t=-1
        return t
    t = mm_idx[0][-1]
    return t


def score_msm(score_data_df):
    return score_data_df['moc'] * score_data_df['spat'] * score_data_df['spec']

import pandas as pd
import numpy as np
from collections import OrderedDict

class decoy_adducts():
#### Some class for dealing with FDR results ####
    def __init__(self,fname,target_adducts,decoy_adducts):
        self.target_adducts = target_adducts
        self.decoy_adducts = decoy_adducts
        # read in raw score file and calculate metabolite signal match
        with open(fname) as f_in:
            self.score_data_df = pd.read_csv(f_in, quotechar ='"').fillna(0)
        self.score_data_df["msm"] = score_msm(self.score_data_df)
        self.score_data_df.sort("sf") #should this be here?
        # store some data info
        self.sf_l = np.unique(self.score_data_df["sf"])
        self.n_sf = len(self.sf_l)

    def decoy_adducts_get_pass_list(self,fdr_target,n_reps,col='msm',return_decoy=False):
        # Get MSM threshold @ target fdr
        # Return molecules with higher MSM value
        msm_vals = self.get_msm_threshold(fdr_target,n_reps)
        pass_list={}
        for a in self.target_adducts:
            target_df = self.score_data_df.ix[self.score_data_df["adduct"]==a]
            pass_list[a] = target_df.ix[target_df[col]>msm_vals[a]][['sf','msm']].values
        if not return_decoy:
            return pass_list
        else:
            pass_list_decoy={}
            for a in self.target_adducts:
                decoy_df = self.score_data_df.ix[~self.score_data_df["adduct"].isin(self.target_adducts)]
                pass_list_decoy[a] = decoy_df.ix[decoy_df[col]>msm_vals[a]][['sf','adduct','msm']].values
            return pass_list,pass_list_decoy

    def get_msm_threshold(self, fdr_target, n_reps=10, col='msm'):
        # Repeatedly calcualte FDR curves
        #   Find target crossing point -> find correspdoning msm score
        # return average score per adduct
        msm_vals = {}
        for a in self.target_adducts:
            msm_vals[a]=[]
            fdr_curves,target_hits,score_vects =self.get_fdr_curve(a,n_reps,col)
            for n in range(n_reps):
                crossing_idx = find_crossing(fdr_curves[n],fdr_target)
                if crossing_idx >-1:
                    msm_vals[a].append(score_vects[n].iloc[crossing_idx])

        # calculate average
        for a in self.target_adducts:
            msm_vals[a] = np.median(msm_vals[a])
        return msm_vals

    def get_fdr_curve(self,adduct,n_reps=10,col='msm'):
        # for a particular adduct, calcualte n_reps fdr curves
        target_df = self.score_data_df.ix[self.score_data_df["adduct"]==adduct]
        data_reps = len(self.score_data_df)/self.n_sf - len(self.target_adducts)
        col_vector_decoy = self.score_data_df.ix[self.score_data_df['adduct'].isin(self.decoy_adducts)][col].values
        col_vector_decoy = col_vector_decoy.reshape((self.n_sf,data_reps))
        _ = [np.random.shuffle(i) for i in col_vector_decoy] #shuffle the values in each row
        fdr_curves = []
        target_hits= []
        score_vects= []
        for n in range(n_reps):
            col_vector=col_vector_decoy[:,n]
            decoy_df = pd.DataFrame({"sf":self.sf_l,col:col_vector})
            fdr_curve,target_hit,score_vect = calc_fdr_df(target_df,decoy_df,col=col,ascending=False)
            fdr_curves.append(fdr_curve)
            target_hits.append(target_hit)
            score_vects.append(score_vect)
        return fdr_curves,target_hits,score_vects

