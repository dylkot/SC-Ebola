import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, fisher_exact, ranksums

def load_geneset(gmtfn, genes=None, minsize=0):
    '''
    Load genesets stored in gmt format (e.g. as provided by msigdb)

    gmtfn : str
        path to gmt file

    genes : list, optional
        only include genes in this input
    
    minsize : int, optional
        minimum geneset size to keep
        
    Returns
    -------
    
    gsets : dict
        gene_set_name : set of genes
        
    allsetgenes : set
        set of genes found in all genesets combined
    '''
    allsetgenes = set()
    if genes is not None:
        genes = set(genes)
    gsets = {}
    effect_min_size = minsize+2 # account for gset name cols
    with open(gmtfn) as F:
        for l in F.readlines():
            words = [x for x in l.rstrip().split('\t')]
            gsetname = words[0]
            setgenes = words[2:]
            if genes is not None:
                setgenes = set(setgenes).intersection(genes)
            else:
                setgenes = set(setgenes[2:])

            if len(setgenes) >= effect_min_size:
                gsets[gsetname] = setgenes
                allsetgenes = allsetgenes.union(setgenes)

    return(gsets, allsetgenes)



def fishertestbygep(gsets, signatures):
    cols = []
    for sig in signatures.columns:
        cols.append((sig, 'OR'))
        cols.append((sig, 'P'))
    
    res = pd.DataFrame(index=list(gsets.keys()), columns=pd.MultiIndex.from_tuples(cols))
    total = res.shape[0]
    N = signatures.shape[0]
    for sig in signatures.columns:
        print(sig)
        siggenes = set(signatures.index[signatures[sig]])
        for (count,gs) in enumerate(res.index):
            (OR,P, table) = run_fisher_exact(siggenes, gsets[gs], N)
            res.at[gs,(sig,'OR')]=OR
            res.at[gs,(sig,'P')]=P
    return(res)


def run_fisher_exact(siggenes, setgenes, num_total):
    numinter = len(setgenes.intersection(siggenes))
    numunion = len(setgenes.union(siggenes))

    table = [[numinter, len(siggenes) - numinter],
                    [len(setgenes)-numinter, num_total-numunion]]
        
    (OR,P) = fisher_exact(table, alternative='two-sided')    
    return(OR, P, table)


def ranksumtestbygep(gsets, signatures):
    cols = []
    for sig in signatures.columns:
        cols.append((sig, 'H'))
        cols.append((sig, 'P'))
    
    res = pd.DataFrame(index=list(gsets.keys()), columns=pd.MultiIndex.from_tuples(cols))
    total = res.shape[0]
    for (count,gs) in enumerate(res.index):
        if (count % 100) == 0:
            print('%d out of %d' % (count, total))
        ind = signatures.index.isin(gsets[gs])
        for sig in signatures.columns:
            x1 = signatures.loc[ind, sig].dropna()
            x2 = signatures.loc[~ind, sig].dropna()
            (H,P) = ranksums(x1, x2)
            res.at[gs,(sig,'H')]=H
            res.at[gs,(sig,'P')]=P
    return(res)
