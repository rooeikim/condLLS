#!/usr/bin/python
# calculate conditional LLS based on filtered gene list
# automatically cut-off network by GAM and LLS 1
# usage) ./thisscript.py [network] [gspath] [genelist]

import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, getopt
import pandas as pd
import numpy as np
import scipy.stats as stats
import seaborn as sns
import os.path
import click
from pygam import LinearGAM
from math import * 

sns.set_style('white')



def load_network_tab(filepath, statname = 'PCC', max_pair=1000000):
    if os.path.isfile(filepath) == False:
        print(f"File {filepath} is not existed!")
        sys.exit(1)
    try:
        if max_pair <= 0:
            network = pd.read_csv(filepath,index_col=(0,1),header=None,sep="\t")[statname]
        else:
            network = pd.read_csv(filepath,index_col=(0,1),header=None,sep="\t", nrows=max_pair)[statname]
    except:
        print(f"Fail to load {filepath}. Check file format")
        sys.exit(1)
    if len(network) == 0:
        print(f"File is empty")
        sys.exit(1)
    return network

def main():

    networkpath = sys.argv[1]#"table.txt"
    gsdatapath = sys.argv[2]#"gmtfile.txt" # [name] [des] [gene1] [g2] [g3] ....
    genelistfile = sys.argv[3]#"genelist.txt"

    # lls_cutoff -> threshold that a GAM model meets selected LLS
    if len(sys.argv) > 4:
        lls_cutoff = int(sys.argv[4])
    else:
        lls_cutoff = 0
    genelist = list()
    with open(genelistfile,'r') as fp:
        for line in fp:
            if line[0] == "#":
                continue
            genelist.append(line.rstrip())

    # load precalculated edges
    edges_all = load_network_tab(networkpath, statname=2) # row : gene, column : sample, no NA values

    # load gold standard
    gspairs = load_gold_standard_pairs(gsdatapath, genelist) 

    # sort edges
    edges_all = pd.Series(edges_all).sort_values(ascending=False).rename(index = "PCC")
                
    # benchmark
    benchdata = run_benchmark_condLLS(edges_all, gspairs )
    benchdata.to_csv(f"{networkpath}.bench.csv")

    # cutoff - Conditional LLS - GAM 
    XX, predictedy, pair_threshold = run_pyGAM(benchdata,outputplot=f"{networkpath}.benchgraph.pdf",lls_threshold=lls_cutoff)
    edges_all_cutoff = edges_all[:int(pair_threshold)]
    edges_all_cutoff.to_csv(f"{networkpath}.cutoff.csv")


def load_gold_standard_pairs(gmtfile,genespace,maxgenes=50):
    networkgenes = set(genespace)
    network = dict()
    with open(gmtfile,"r") as fp:
        for line in fp:
            linearray = line.rstrip().split("\t")
            termid = linearray.pop(0)
            desc = linearray.pop(0)
            if len(linearray)>maxgenes:  # skip large terms
                continue
            genes = list(networkgenes&set(linearray))
            for i in range(len(genes)):
                if genes[i] not in network:
                        network[genes[i]] = dict()
                for j in range(i+1,len(genes)):
                    if genes[j] not in network:
                        network[genes[j]] = dict()
                    network[genes[i]][genes[j]] = True
                    network[genes[j]][genes[i]] = True
    return network


def run_benchmark_condLLS(target_network, gold_standard_pairs, already_sorted = True, ascending = False, max_pair = 1000000,binsize=1000):

    # calculate positive, negative pairs

    pos = 0
    for g in gold_standard_pairs:
        pos+=len(gold_standard_pairs[g])
    pos /= 2
    neg = ( len(gold_standard_pairs) * (len(gold_standard_pairs)-1) / 2 ) - pos

    # calculate pr curve

    TP = 0
    FP = 0
    FN = pos
    binTP = 0
    binFP = 0

    count = 0

    #print("%f"%(float(pos)/float(neg)))

    randomprecision = float(pos)/float(pos+neg)
    priorprob = float(pos)/float(neg)
    #print(priorprob)
    #print randomprecision
    already = set()
    benchdata = dict()
    binstats=list()

    if already_sorted == False:
        sorted_keys = sorted(target_network,key=target_network.get,reverse = not ascending)
    else:
        sorted_keys = target_network.keys()


    for pair in sorted_keys:
        if count>=max_pair:
            break

        already.add(pair[0])
        already.add(pair[1])


        binstats.append(target_network[pair])
        count+=1
        if pair[0] in gold_standard_pairs and pair[0] in gold_standard_pairs:
            if pair[0] in gold_standard_pairs and pair[1] in gold_standard_pairs[pair[0]]:
                TP += 1
                FN -= 1
                binTP+= 1
            else:
                FP += 1
                binFP += 1

        if count%binsize == 0:
            precision = float(TP)/float(TP+FP)
            recall = float(TP)/float(TP+FN)
            oddratio = precision/randomprecision
            if FP!=0:
                lls = log((float(TP)/float(FP)) / priorprob) / log(2)
            else:
                lls = 20.0
            if binTP!=0:
                if binFP!=0:
                    binlls = log((float(binTP)/float(binFP)) / priorprob) / log(2)
                else:
                    binlls = 20.0
            else:
                binlls = -5.0
            
            # save result
            benchdata[count] = {"CondLLS":lls,"TP":TP,"FP":FP,"MeanBinStatistics":np.mean(binstats),"GeneCoverage":len(already),"BinCondLLS":binlls}
            #print("%d\t%.3f\t%d\t%d\t%d\t%f"%(count,lls,TP,FP,len(already),binlls))


            # initialize
            binTP=0
            binFP=0
            binstats=list()


    if count%binsize != 0:
            
        precision = float(TP)/float(TP+FP)
        recall = float(TP)/float(TP+FN)
        oddratio = precision/randomprecision
        lls = log((float(TP)/float(FP)) / priorprob) / log(2)

        if binTP!=0:
            if binFP!=0:
                binlls = log((float(binTP)/float(binFP)) / priorprob) / log(2)
            else:
                binlls = 20.0
        else:
            binlls = -5.0
        benchdata[count] = {"CondLLS":lls,"TP":TP,"FP":FP,"MeanBinStatistics":np.mean(binstats),"GeneCoverage":len(already),"BinCondLLS":binlls}
        #print("%d\t%.3f\t%d\t%d\t%d\t%f"%(count,lls,TP,FP,len(already),binlls))
    return pd.DataFrame().from_dict(benchdata).T


    
def run_pyGAM(benchdata,outputplot="GAM_bench.pdf",lls_threshold=0.0,min_pair_n = 1000):
    
    
    X = np.asarray([np.asarray([c]) for c in benchdata.index])
    y = benchdata['BinCondLLS'].values
    gam = LinearGAM(n_splines=10).gridsearch(X, y)

    gam.summary()
    plt.figure()
    XX = gam.generate_X_grid(term=0, n=500)
    plt.plot(XX, gam.predict(XX), 'r--')
    plt.plot(XX, gam.prediction_intervals(XX, width=.95), color='b', ls='--')
    plt.scatter(X, y, facecolor='gray', edgecolors='none')
    plt.xlabel("Rank of Pairs")
    plt.ylabel("Conditional LLS")
    plt.savefig(f"{outputplot}",format='pdf')

    predictedy = gam.predict(XX)
    pair_threshold = min_pair_n
    for i in range(len(XX)):
        if predictedy[i] > lls_threshold:
            pair_threshold = XX[i]



    print("# of significant pairs = %f"%pair_threshold)

    return XX, predictedy, pair_threshold

if __name__ == "__main__":
    main()

