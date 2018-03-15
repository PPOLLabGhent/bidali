#!/usr/bin/env python3
#R exposed top routines
import pickle, re, sys, operator
from bidali import LSD
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

# Differential expression analysis
def prepareDesign(metadata,design,reflevels=None,RReturnOnly=True):
    """
    Takes pandas metadata table and design string
    to return R design object

    e.g. design = '~0+treatment+batch'
    """
    metacols_r = {
        col: ro.r.relevel(
            ro.r.factor(metadata[col]),
            ref=reflevels[col] if reflevels else None
        )
        for col in metadata
    }
    metadata_r = ro.r['data.frame'](**metacols_r)
    design_r = ro.r['model.matrix'](ro.r.formula(design),data=metadata_r)
    design_r.colnames = ro.StrVector(
        [c.replace(':','.') for c in ro.r.colnames(design_r)]
    ) # resolves naming issue when using interaction factors in design
    #design = pd.DataFrame(np.array(design_r),columns=design_r.colnames,index=design_r.rownames)
    #design.doxnox*2 + (-design.shairpinTBX2sh25 - design.shairpinTBX2sh27)*3*design.doxnox
    if RReturnOnly:
        return design_r
    else:
        design_p = pd.DataFrame(pandas2ri.ri2py(design_r), index = metadata.index, columns = design_r.colnames)
        return design_r, design_p
    
def prepareContrasts(design_r,contrasts, RReturnOnly = True):
    """
    Expects R design_r and list of contrasts that contain design_r column names
    Returns contrasts_r and also pd.DataFrame of contrasts_r if not RReturnOnly

    e.g. contrasts = [
      'treatmentSHC002_dox - treatmentTBX2sh25_dox','treatmentSHC002_dox - treatmentTBX2sh27_dox',
      'treatmentSHC002_nox - treatmentTBX2sh25_nox','treatmentSHC002_nox - treatmentTBX2sh27_nox'
    ]
    """
    # import R limma package
    limma = importr('limma')
    contrasts_r = limma.makeContrasts(*contrasts,levels=design_r)
    if RReturnOnly:
        return contrasts_r
    else:
        contrasts_p = pd.DataFrame(
            pandas2ri.ri2py(contrasts_r),
            columns=ro.r.colnames(contrasts_r),
            index=ro.r.rownames(contrasts_r)
        )
        return contrasts_r, contrasts_p

def DEA(counts,design_r,contrasts=None):
    """
    contrasts needs to be a dictionary with string keys and values to be used as design contrasts or 
    a list of int values to be used as design coefficients

    Returns results and pd.DataFrame of normalized counts
    """
    if isinstance(contrasts, list):
        coefs = True
        if not all(isinstance(k, int) for k in contrasts):
            raise ValueError('coefficient list should be all integers')
        coefnames = list(ro.r.colnames(design_r))
        contrasts = OrderedDict([(coefnames[c-1],c) for c in contrasts])
    else:
        coefs = False
        if not (all(isinstance(k, str) for k in contrasts) and
                all(isinstance(v, str) for v in contrasts.values())):
            raise ValueError('contrast dict should be all string pairs')
    
    # import R limma package
    limma = importr('limma')
    # tranform counts with voom
    voomedCounts_r = limma.voom(counts,design=design_r,plot=True,normalize="quantile")
    fit_r = limma.lmFit(voomedCounts_r,design_r)
    fit_r = limma.eBayes(fit_r)
    coefficients_r = fit_r.rx2('coefficients') #fit_r$coefficients
    if coefs:
        fit_contrasts_r = fit_r
    else:
        contrasts_r, contrasts_p = prepareContrasts(design_r, contrasts.values(), RReturnOnly = False)
        fit_contrasts_r = limma.contrasts_fit(fit_r, contrasts_r)
        fit_contrasts_r = limma.eBayes(fit_contrasts_r)
    print(ro.r.summary(fit_contrasts_r))

    #Full results
    results = OrderedDict()
    for res in contrasts:
        result_r = limma.topTable(
            fit_contrasts_r,coef=contrasts[res],n=len(counts)
        )
        results[res] = pandas2ri.ri2py(result_r)
        results[res].index = ro.r.rownames(result_r)
        #results[res]['gene_label'] = results[res].index.map(lambda x: counts.index[int(x)-1])
        print('# sig',res,'->',(results[res]['adj.P.Val']<=0.05).sum())
        
    return results, pd.DataFrame(
        pandas2ri.ri2py(voomedCounts_r.rx2('E')),
        columns=counts.columns,index=counts.index
    )
#limma volcanoplot and plotMDS need to be done either in qconsole or direct R environment

# Gene-set enrichment analysis
## limma based
def get_gcindices_r(countsGeneLabels,correctBackground=False,remove_empty=True):
    """
    >>> indices_r = get_gcindices(counts.index,correctBackground=False) # doctest: +SKIP

    TODO overlapping with genesets2indices_r => refactor code
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

def genesets2indices_r(genesets, geneLabels, remove_empty=True):
    """
    genesets should be a dictionary of geneset lists, and geneLabels a list of gene labels
    >>> from .tests.test_retro import testGenesets, testCountsGenelabels
    >>> print(genesets2indices_r(testGenesets,testCountsGenelabels))
    [[1]]
    [1] 1 2
    <BLANKLINE>
    [[2]]
    [1] 2 3
    <BLANKLINE>
    <BLANKLINE>

    TODO test function for gene level conversion between python and R
    """
    limma = importr('limma')
    genesets_r = ro.ListVector({gs:ro.StrVector(genesets[gs]) for gs in genesets})
    return limma.ids2indices(genesets_r, geneLabels, remove_empty = remove_empty)

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

## fgsea based
def fgsea(genesets_r, ranks, minSize = 15, maxSize = 500, nperm = 10000):
    """
    Ranks should be pd.Series with as index all genes in geneLabels used for generating genesets_r
    with e.g. genesets2indices_r

    TODO test function for gene level conversion between python and R
    """
    fgsea = importr('fgsea')
    fgseaRes_r = fgsea.fgsea(
        pathways = genesets_r, stats = pd.Series(list(ranks),index=range(1,len(ranks)+1)),
        minSize = minSize, maxSize = maxSize, nperm = nperm
    )
    fgseaRes = pandas2ri.ri2py(fgseaRes_r.rx(True,-8))
    fgseaRes['leadingEdge'] = [
        [ranks.index[int(gn)-1] for gn in fgseaRes_r.rx2(i,'leadingEdge')]
        for i in range(1,len(fgseaRes_r.rownames)+1)
    ]
    return fgseaRes.sort_values('padj')

## rank sums
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
