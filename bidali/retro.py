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
def prepareDesign(metadata,design,reflevels=None):
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
    print(design_r.colnames)
    return design_r

def prepareContrasts(design_r,contrasts):
    """
    Expects R design_r and list of contrasts that contain design_r column names
    Returns contrasts_r and pd.DataFrame of contrasts_r

    e.g. contrasts = [
      'treatmentSHC002_dox - treatmentTBX2sh25_dox','treatmentSHC002_dox - treatmentTBX2sh27_dox',
      'treatmentSHC002_nox - treatmentTBX2sh25_nox','treatmentSHC002_nox - treatmentTBX2sh27_nox'
    ]
    """
    # import R limma package
    limma = importr('limma')
    contrasts_r = limma.makeContrasts(*contrasts,levels=design_r)
    return contrasts_r, pd.DataFrame(
        pandas2ri.ri2py(contrasts_r),
        columns=ro.r.colnames(contrasts_r),
        index=ro.r.rownames(contrasts_r)
    )
    
def DEA(counts,design_r,contrasts=None,coefs=None):
    """
    coefs can be given instead of contrasts, when the contrasts match the design parameters directly
    Returns results and pd.DataFrame of normalized counts
    """
    # import R limma package
    limma = importr('limma')
    # tranform counts with voom
    voomedCounts_r = limma.voom(counts,design=design_r,plot=True,normalize="quantile")
    fit_r = limma.lmFit(voomedCounts_r,design_r)
    fit_r = limma.eBayes(fit_r)
    coefficients_r = fit_r.rx2('coefficients') #fit_r$coefficients
    if not (contrasts or coefs):
        return coefficients_r
    elif coefs:
        fit_contrasts_r = fit_r
        contrasts = coefs
    else:
        fit_contrasts_r = limma.contrasts_fit(fit_r,contrasts)
        fit_contrasts_r = limma.eBayes(fit_contrasts_r)
        print(ro.r.summary(fit_contrasts_r))

    #Full results
    results = {}
    for res in contrasts:
        result = limma.topTable(fit_contrasts_r,coef=res,n=len(counts))
        results[res] = result
        results[res]['gene_label'] = results[res].index.map(lambda x: mf[x])
        results[res]['gene_label_signed'] = results[res].logFC.map(
            lambda x: '+' if x > 0 else '-')+results[res].gene_label
        print('# sig',res,'->',(results[res]['adj.P.Val']<=0.05).sum())
        
    return results, pd.DataFrame(
        pandas2ri.ri2py(voomedCounts_r.rx2('E')),
        columns=counts.columns,index=counts.index
    )
#limma volcanoplot and plotMDS need to be done either in qconsole or direct R environment

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
