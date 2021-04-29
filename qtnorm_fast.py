#!/usr/bin/python

import sys
import pandas as pd


tablefile = sys.argv[1]
outputfile = tablefile+".qtnorm"

bfdata = pd.read_csv(tablefile,sep="\t",index_col=0,header=0)

rank_mean = bfdata.stack().groupby(bfdata.rank(method='first').stack().astype(int)).mean()


bfdata.rank(method='min').stack().astype(int).map(rank_mean).unstack().to_csv(outputfile,sep="\t")
