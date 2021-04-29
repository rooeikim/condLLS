#!/usr/bin/python

import pandas as pd
import numpy as np
import sys

inputfile= sys.argv[1]

bfdata = pd.read_table(inputfile,header=0,index_col=0)

def quantileNormalize(df_input):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic[col] = df[col].sort_values(na_position='first').values
    sorted_df = pd.DataFrame(dic)
    #rank = sorted_df.mean(axis = 1).tolist()
    rank = sorted_df.median(axis = 1).tolist()
    #sort
    for col in df:
        # compute percentile rank [0,1] for each score in column 
        t = df[col].rank( pct=True, method='max' ).values
        # replace percentile values in column with quantile normalized score
        # retrieve q_norm score using calling rank with percentile value
        df[col] = [ np.nanpercentile( rank, i*100 ) if ~np.isnan(i) else np.nan for i in t ]
    return df


bfdata_qtnorm = quantileNormalize(bfdata)

bfdata_qtnorm.to_csv(inputfile+".qtnorm",sep='\t')
