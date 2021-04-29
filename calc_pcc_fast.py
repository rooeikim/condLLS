#!/usr/bin/python
import sys, getopt
import pandas as pd
import numpy as np
import scipy.stats as stats
import os.path
from math import * 


def func_pcc(x,y,pre_sqsumroot_x,pre_sqsumroot_y):
    return np.matmul(x.transpose(),y)/(pre_sqsumroot_x*pre_sqsumroot_y)

def calc_pearsonr_fast(essdata_df,outputfile):
    # custom pearsonr
    if essdata_df.isnull().values.any() == True:
        print("There is NaN in the table -> Run slow")
        return calc_pearsonr_slow(essdata_df)
    essdata_submean = dict()
    for g in essdata_df.index:
        essdata_submean[g] = (essdata_df.loc[g] - essdata_df.loc[g].mean()).to_numpy()
    essdata_submean_sqsumroot = dict()
    for g in essdata_submean.keys():
        essdata_submean_sqsumroot[g] = (essdata_submean[g]**2).sum()**0.5

    #edges = dict()
    with open(outputfile,"w") as fout:
        for g1 in essdata_submean.keys():
            for g2 in essdata_submean.keys():
                if g1<g2:
                    pcc = func_pcc(essdata_submean[g1],essdata_submean[g2],essdata_submean_sqsumroot[g1],essdata_submean_sqsumroot[g2])
                    #edges[(g1,g2)] = pcc
                    fout.write(f"{g1}\t{g2}\t{pcc}\n")
    #return edges

inputfile = sys.argv[1]
outputfile = inputfile+".pcc"

bfdata = pd.read_csv(inputfile,sep="\t",index_col=0,header=0)

edges_all = calc_pearsonr_fast(bfdata,outputfile)




