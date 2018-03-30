# -*- coding: utf-8 -*-
"""Expression analysis module

Defines the Expan object for all your expression analysis needs.
Relies under the hood on the retro module and R packages.

Example:
    >>> from bidali.expan import Expan
    >>> expan = Expan(                                                                                                      
    ...   counts = '/home/christophe/Dropbiz/Basecamp/BRIP1/2016/2016_008_DDR_RNA sequencing/raw data-results/2015_TMPYP4_SVH_RSEMcounts_pp.csv',
    ...   metadata = '/home/christophe/Dropbiz/Basecamp/BRIP1/2016/2016_008_DDR_RNA sequencing/raw data-results/2015_TMPYP4_metadata.csv'
    ... )
    >>> expan.designator(design = '~batch+treatment+cellline', reflevels = {'batch':'seq1','cellline':'IMR32','treatment':'control'})
    >>> expan.exdif(contrasts = [4, 5], countfilter = 1)
    >>> brip1 = expan['BRIP1']
    >>> brip1.plotCounts('cellline','treatment')
"""
import pandas as pd
import numpy as np
from bidali import retro
from collections import OrderedDict

class Expan:
    """Expression analysis class

    This class is the starting point for a bidali expression analysis workflow.
    """

    def __init__(self,
                 counts, metadata, export = None,
                 annotations = None, annotatedOnly = True,
                 counts_kwargs={}, metadata_kwargs={}, **kwargs):
        """Expression analysis class

        Args:
            counts (str): The filename for the counts table.
            metadata (str): The filename for the metadata table.
            export (str): Export location. Not implemented yet!
            annotations (pd.DataFrame): Annotations for the counts
                gene identifiers. Default loads biomart table.
            annotatedOnly (bool): Only analyze genes that have annotation.
            counts_kwargs (dict):
            metadata_kwargs (dict):
                Common counts and metadata pd.read_table arguments can 
                be provided as extra key word arguments. Different arguments
                for pd.read_table need to be specified in the counts_kwargs
                and metadata_kwargs dictionaries. The default index col is 0
                and the default sep is ','. 
        """
        # Include common kwargs in counts and metadata specific kwargs
        default_kwargs = {'index_col':0,'sep':','}
        default_kwargs.update(kwargs)
        counts_kwargs.update(default_kwargs)
        metadata_kwargs.update(default_kwargs)
        # Prepare tables
        self.counts = pd.read_table(counts,**counts_kwargs)
        ## Prepare metadata
        if isinstance(metadata,str):
            self.metadata = pd.read_table(metadata,**metadata_kwargs)
        else: # prepare metadata from counts columns
            if 'rename_cols' in metadata:
                rename_cols = metadata.pop('rename_cols')
                self.counts.columns = [rename_cols(c) for c in self.counts.columns]
            self.metadata = pd.DataFrame({
                m:[metadata[m](c) for c in self.counts.columns] for m in metadata
            }, index = self.counts.columns)
        if len(self.metadata) != sum(self.counts.columns.isin(self.metadata.index)):
            raise Exception('count column names and metatada row names need to match!')
        # Set column order of counts, to row order of metadata
        self.counts = self.counts[self.metadata.index]
        self._robjects = {} #dict to save private R objects
        if annotations is not None: self.annotations = annotations
        else:
            from bidali.LSD.dealer.external.ensembl import get_biomart
            biomart = get_biomart()
            biomart = biomart[~biomart['Gene stable ID'].duplicated()]
            biomart = biomart.set_index('Gene stable ID')
            self.annotations = biomart
        if annotatedOnly:
            self.counts = self.counts[self.counts.index.isin(self.annotations.index)]
        if export:
            self.exporter(export)

    def __getitem__(self,key):
        """
        First assumes a gene name, and queries the different tables
        for all gene-related information
        """
        if key in self.annotations.index:
            geneid = key
        else:
            geneid = self.annotations.index[self.annotations['Gene name'] == key]
        return GeneResult(
            name = key,
            counts = self.counts.loc[geneid],
            counts_norm = self.counts_norm.loc[geneid],
            exdif = OrderedDict([(k,self.results[k].loc[geneid]) for k in self.results]),
            parent = self
        )

    def __delitem__(self,sample):
        """
        Deletes a sample from counts and metadata
        """
        self.counts = self.counts.drop(sample,axis=1)
        self.metadata = self.metadata.drop(sample)
        

    def designator(self, design, reflevels):
        """Prepare experimental design

        Args:
            design (str): The design of the experiment, e.g. '~treatment+batch'
            reflevels (dict): The reference level for each factor of the metadata

        Example:
            reflevels = {'treatment':'SHC002_nox','rep':'rep1'}
        """
        self.design = design
        self._robjects['design'], self.designmatrix = retro.prepareDesign(
            self.metadata[self.metadata.columns[self.metadata.columns.isin(reflevels)]],
            self.design, reflevels, RReturnOnly = False
        )

    def exdif(self, contrasts, countfilter=1, quantro = False):
        """Differential expression method

        Runs differential expression workflow

        Args:
            contrasts (list or dict): If design coefficients suffice provide
                a list of the design column numbers of interest, else provide
                a dict with contrast name and contrast constructor string.
            countfilter (int): Mininum average number of reads per sample.
                Genes that do not meet this requirement are filtered.
            quantro (bool): Include quantro analysis. Not implemented yet!

        Examples:
            >>> contrasts = [1,2]
            >>> contrasts = {
            ...   'c1': 'treatmentSHC002_dox - treatmentTBX2sh25_dox','treatmentSHC002_dox - treatmentTBX2sh27_dox',
            ...   'c2': 'treatmentSHC002_nox - treatmentTBX2sh25_nox','treatmentSHC002_nox - treatmentTBX2sh27_nox'
            ...  }
        """
        self.contrasts = contrasts
        self.counts_fltd = (
            self.counts[self.counts.T.sum() > countfilter * len(self.counts.columns)]
            if countfilter else self.counts
        )
        self.results, self.counts_norm = retro.DEA(self.counts_fltd, self._robjects['design'], self.contrasts)
        for k in self.results:
             self.results[k]['Gene name'] = self.annotations['Gene name']

    def gsea(self, genesets = None, exclude_results = {}):
        """Geneset enrichment analysis

        Expan.exdif should have already been executed. The gene exdif results will be deduplicated.
        As such, genes that have multiple locations in the genome will only maintain there best result
        in the geneset enrichment analysis. This creates a small bias for such genes.

        Args:
            genesets (dict): Should be a dict with string keys designating the geneset names,
                and gene list values containing all the geneset members. If not provided all
                MSigDB geneset collections are used.
            exclude_results (set): A set of exdif results for which the gsea analysis is not required.
        Returns:
            gsea result
        """
        from bidali.retro import genesets2indices_r, fgsea
        if genesets:
            gscollections = {'custom': genesets}
        else:
            from bidali import LSD
            gscollections = LSD.get_msigdb6()
        # deduplicate gene exdif results
        results_deduplicated = {
            r:self.results[r][~self.results[r]['Gene name'].duplicated()]
            for r in self.results if r not in exclude_results
        }
        self.gsea_results = {}
        for r in results_deduplicated:
            geneLabels = results_deduplicated[r]['Gene name'].as_matrix()
            gscollections_r = {
                gsc:genesets2indices_r(gscollections[gsc], geneLabels)
                for gsc in gscollections
            }
            self.gsea_results[r] = {
                gsc:fgsea(
                    gscollections_r[gsc],
                    results_deduplicated[r].set_index('Gene name').t,
                    minSize = 15
                )
                for gsc in gscollections   
            }

    def exporter(self, location):
        """Export expression analysis object
        
        All relevant tables and figures are saved in a zipfolder at the specified location.
        """
        pass

    @staticmethod
    def importer(location):
        """Import tables from a previous Expan

        Only counts, metadata and annotation are imported.
        """
        pass

    @staticmethod
    def convert_excell_to_csv(filename,sheet_name=0):
        """
        Convenience function to e.g. convert an metadata excell table to a csv equivalent

        Args:
            filename (str): Excell filename.
            sheet_name (int or str): Either sheet number or sheet name.

        Returns:
            filename of converted table
        """
        newFilename = filename[:filename+rindex('.')]+'.csv'
        pd.read_excel(filename,sheet_name=sheet_name).to_csv(newFilename)
        return newFilename

class GeneResult:
    """
    Collects gene specific results.
    Methods for plotting, calculating knockdown
    """
    def __init__(self, name, counts, counts_norm, exdif, parent):
        """
        Differential expression gene result
        """
        self.name = name
        self.counts = counts
        self.counts_norm = counts_norm
        self.exdif = exdif
        self.parent = parent

    def __repr__(self):
        return '<GeneResult {}>'.format(self.name)

    def __str__(self):
        return '{}:\nCounts:\n{}\n\nNormalised:\n{}\n\nResults:\n{}\n'.format(
            self.name,
            self.counts.T,
            self.counts_norm.T,
            '\n'.join(['{} =>\n{}\n'.format(k,self.exdif[k].T) for k in self.exdif])
        )

    def calcKD(self, groupingCol, knockdownCol, knockdownControl = 'control'):
        """Calculate knockdown

        Args:
            groupingCol (str): Grouping column name in analysis metadata.
            knockdownCol (str): KD column name in analysis metadata.
            knockdownControl (str): Within knockdownCol the value of the control condition

        Returns:
            float

        Note:
            * Missing values should be provided as np.NaN
        """
        treatments = {t for t in self.parent.metadata[knockdownCol] if t != knockdownControl}
        kd = self.counts.groupby(self.parent.metadata[groupingCol]
        ).aggregate(
            {t:lambda x,t=t: (x.filter(regex=t).get_values()[0] if x.filter(regex=t).any() else np.nan)/
             x.filter(regex=knockdownControl).get_values()[0]
             for t in treatments})
        return kd

    def plotCounts(self, x, hue, normalisedCounts = True, **kwargs):
        """
        Plot counts for a gene
        """
        from .visualizations import plotGeneCounts
        df = (self.counts_norm if normalisedCounts else self.counts).T.join(self.parent.metadata)
        df.columns = ['counts'] + list(df.columns[1:])
        return plotGeneCounts(df, x=x, hue=hue, **kwargs)

class GenesetResult(GeneResult):
    """
    Inherits from GeneResult.
    Offers extra methods for geneset related functionality
    """
    def __repr__(self):
        return '<GenesetResult {}>'.format(self.name)
