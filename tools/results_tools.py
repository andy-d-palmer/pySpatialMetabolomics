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
