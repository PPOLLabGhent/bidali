# -*- coding: utf-8 -*-
"""Expression analysis module

Defines the Expan object for all your expression analysis needs.
Relies under the hood on the retro module and R packages.
"""
import pandas as pd
from bidali import retro
from collections import OrderedDict

class Expan:
    """
    This class is the starting point for a bidali expression analysis workflow.
    """

    def __init__(self,
                 counts, metadata, export = None,
                 annotations = None, annotatedOnly = True,
                 counts_kwargs={'index_col':0},
                 metadata_kwargs={'index_col':0}, **kwargs):
        """Expression analysis class

        Expects the filename for the counts and metadata table.
        Common counts and metadata pd.read_table arguments can 
        be provided as extra key word arguments. Different arguments
        for pd.read_table need to be specified in the counts_kwargs
        and metadata_kwargs dictionaries
        """
        self.counts = pd.read_table(counts,**counts_kwargs)
        self.metadata = pd.read_table(counts,**metadata_kwargs)
        if annotations: self.annotations = annotations
        else:
            from bidali.LSD.dealer.external.ensembl import get_biomart
            biomart = get_biomart()
            biomart = biomart[~biomart['Gene stable ID'].duplicated()]
            biomart = biomart.set_index('Gene stable ID')
            self.annotations = biomart
        if annotatedOnly:
            self.counts = self.counts[self.counts.index.isin(self.annotations.index)]
        if export:
            self.export(export)

    def designator(self):
        pass

    def exdif(self, design, reflevels, contrasts, coefs = True, countfilter=1, quantro = False):
        """Differential expression method

        Runs differential expression workflow

        e.g. reflevels = {'treatment':'SHC002_nox','rep':'rep1'}

        e.g. contrasts = {'c1': 1, 'c2': 2} if coefs == True 
        else contrasts = {
          'c1': 'treatmentSHC002_dox - treatmentTBX2sh25_dox','treatmentSHC002_dox - treatmentTBX2sh27_dox',
          'c2': 'treatmentSHC002_nox - treatmentTBX2sh25_nox','treatmentSHC002_nox - treatmentTBX2sh27_nox'
        }
        """
        self.counts_fltd = (
            self.counts[self.counts.T.sum() > countfilter * len(self.counts.columns)]
            if countfilter else self.counts
        )
        self.design = design
        design_r = retro.prepareDesign(self.metadata,self.design,reflevels)
        if coefs:
            self.results, self.counts_norm = retro.DEA(self.counts,self.design_r,coefs=contrasts.values())
        else:
            contrasts_r, self.contrasts = retro.prepareContrasts(design_r,contrasts.values())
            self.results, self.counts_norm = retro.DEA(self.counts,self.design_r,contrasts=contrasts.values())

    def export(location):
        """Export expression analysis object
        
        All relevant tables and figures are saved in a zipfolder at the specified location.
        """
        pass

    @class_method
    def import(location):
        """Import tables from a previous Expan

        Only counts, metadata and annotation are imported.
        """
        pass
