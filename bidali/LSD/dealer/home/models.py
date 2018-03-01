#!/usr/bin/env python
from bidali import LSD
import pandas as pd, numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects as ro

#Activate automatic pandas/r conversion
from rpy2.robjects import pandas2ri
pandas2ri.activate()

@LSD.storeDatasetLocally
def get_THMYCN():
    """
    Source: ~/Dropbiz/Lab/z_archive/Datasets/2015_Hyperplasia_ABC/uniData.tx
    Source: ~/Dropbiz/Lab/z_archive/Datasets/2015_Hyperplasia_ABC/SampleAnnotation.txt
    """
    #TH-MYCN mouse model incorporation
    importr('limma')
    importr('preprocessCore')
    thannot = pd.read_table(LSD.datadir+'2015_Hyperplasia_ABC/SampleAnnotation.txt')
    thannot.group = thannot.group.apply(lambda x: x.replace('.','_'))
    thdata = pd.read_table(LSD.datadir+'2015_Hyperplasia_ABC/uniData.txt',sep=' ')
    thdata.columns = thannot.group + (['_r'+str(i+1) for i in range(4)]*3*2)
    genotype = thannot.genotype
    age = thannot.age
    ro.globalenv["thdata"] = thdata
    ro.globalenv["genotype"] = ro.r.relevel(ro.FactorVector(genotype),ref='WT')
    ro.globalenv["age"] = ro.IntVector(age)
    ro.r('thdatanorm <- normalize.quantiles(as.matrix(thdata))')
    ro.r('colnames(thdatanorm) <- colnames(thdata)')
    ro.r('rownames(thdatanorm) <- rownames(thdata)')
    ro.r('designSamples = model.matrix(formula(~genotype:age))')
    ro.r('colnames(designSamples)')
    ro.r('fit <- lmFit(thdatanorm,designSamples)')
    ro.r('fit <- eBayes(fit)')
    thcoeffs = ro.r('data.frame(fit$coefficients)')
    thcoeffs['lineardiff'] = thcoeffs['genotypeTG.age']-thcoeffs['genotypeWT.age']

    return LSD.Dataset(exprData=thdata,
                       metadata=thannot,
                       lmCoeffs=thcoeffs)
