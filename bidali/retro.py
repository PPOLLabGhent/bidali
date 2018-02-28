#!/usr/bin/env python3
#R exposed top routines
import pickle, re, sys, operator, LSD
import pandas as pd, numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt
from rpy2.rinterface import RRuntimeError
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
#Activate automatic pandas/r conversion
from rpy2.robjects import pandas2ri
pandas2ri.activate()

# General importr's
base = importr('base')
stats = importr('stats')
utils = importr('utils')

# Gene-set enrichment analysis
def get_gcindices_r(countsGeneLabels,correctBackground=False,remove_empty=True):
    """
    >>> indices_r = get_gcindices(counts.index,correctBackground=False) # doctest: +SKIP
    """
    limma = importr('limma')
    gc = LSD.get_msigdb6()

    if correctBackground:
        gc = {gsc:{gs:[g for g in gc[gsc][gs] if g in countsGeneLabels] for gs in gc[gsc]} for gsc in gc}
        
    countsGeneLabels_r = ro.StrVector(countsGeneLabels)        
    gc_indices_r = {gsc:limma.ids2indices(ro.ListVector(gc[gsc]),countsGeneLabels_r,remove_empty=remove_empty)
                    for gsc in gc}
    
    return gc_indices_r

def romer(counts_r,index_r,design_r,contrast_r):
    """
    counts_r -> should be log-expression values (voom transformed RNAseq counts also good)
    index_r -> collection of geneses (ListVector of IntVector's)
    design_r -> design matrix
    contrast_r -> coef or contrast vector
      e.g. contrasts_r.rx(True,"celllineIMR32.shairpinsh25+celllineIMR32.shairpinsh27"

    >>> romer(counts_r,index_r,design_r,contrast_r) # doctest: +SKIP
    """
    limma = importr('limma')
    rr = limma.romer(counts_r,index=index_r,design=design_r,contrast=contrast_r)
    rr_py = pd.DataFrame(pandas2ri.ri2py(rr),columns=base.colnames(rr),index=base.rownames(rr))
    rr_py['Up_adj'] = stats.p_adjust(rr.rx(True,'Up'),method='fdr')
    rr_py['Down_adj'] = stats.p_adjust(rr.rx(True,'Down'),method='fdr')
    rr_py['min_p_adj'] = rr_py[['Up_adj','Down_adj']].T.min()
    rr_py.sort_values('min_p_adj',inplace=True)
    return rr_py

def romer_fullMSigDB(counts_r,countsGeneLabels,**kwargs):
    """
    >>> romer_fullMSigDB(counts_r,countsGeneLabels) # doctest: +SKIP
    """
    gc = get_gcindices_r(countsGeneLabels)
    return OrderedDict([
        (gsc,romer(counts_r,index_r=gc[gsc],**kwargs))
        for gsc in sorted(gc)
    ])

def mroast(counts_r,index_r,design_r,contrast_r,set_statistic="mean"):
    """
    >>> mroast(counts_r,index_r,design_r,contrast_r) # doctest: +SKIP
    """
    limma = importr('limma')
    return base.data_frame(limma.mroast(counts_r,index_r,design=design_r,contrast=contrast_r,set_statistic=set_statistic))

def genesets2indices_r(genesets,countsGeneLabels):
    """
    genesets should be a dictionary of genesets
    >>> from .tests.test_retro import testGenesets, testCountsGenelabels
    >>> print(genesets2indices_r(testGenesets,testCountsGenelabels))
    [[1]]
    [1] 1 2
    <BLANKLINE>
    [[2]]
    [1] 2 3
    <BLANKLINE>
    <BLANKLINE>
    """
    limma = importr('limma')    
    genesets_r = ro.ListVector(genesets)
    return limma.ids2indices(genesets, countsGeneLabels)

def barcodeplot(fit_r,contrast,indices,geneset,geneset2=None,pngname=None,width=512,height=512,**kwargs):
    """
    >>> barcodeplot(fit_r,3,indices['C5'],'GO_GLUTATHIONE_TRANSFERASE_ACTIVITY',pngname=f.name) # doctest: +SKIP
    """
    limma = importr('limma')
    if pngname:
        grdevices = importr('grDevices')
        grdevices.png(file=pngname, width=width, height=height)
    limma.barcodeplot(fit_r.rx2('t').rx(True,contrast),indices.rx2(geneset),
                      index2 = indices.rx2(geneset2) if geneset2 else ro.NULL,**kwargs)
    if pngname: grdevices.dev_off()

def rankSumProbsMSigDB(ranks,universe,adjpmethod='fdr'):
    from bidali.fegnome import RankSumResult,rankSumProbability
    gc = LSD.get_msigdb6()
    results = OrderedDict()
    for gsc in sorted(gc):
        results[gsc] = pd.DataFrame({gs:rankSumProbability(ranks,gc[gsc][gs]) for gs in gc[gsc]}).T
        results[gsc].columns = RankSumResult._fields
        results[gsc]['fepa'] = stats.p_adjust(results[gsc].fe_p,method=adjpmethod)
        results[gsc]['respa'] = stats.p_adjust(results[gsc].probability,method=adjpmethod)
        results[gsc].sort_values(['fepa','respa'],inplace=True)

    return results
