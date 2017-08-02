#!/usr/bin/env python
import LSD, gzip
import pandas as pd, numpy as np
from os.path import expanduser, exists
from itertools import count
from LSD import storeDatasetLocally, datadir, Dataset

#TODO => NRC segmentation data from ~/Dropbiz/Lab/z_archive/AFW/project CONEXIC/NRC data/01. JISTIC

@storeDatasetLocally
def get_NRC(datadir=datadir+'NRC_data_AFW/'):
    metadata = pd.read_excel(datadir+'20111216_NRC_samples.xlsx',skiprows=4)
    metadata.pop('Unnamed: 0')
    metadata.index = metadata.pop('NRC_ID')
    exprdata = pd.read_table(datadir+'20170214_untransformedR2_expression.txt',skiprows=48,
                         names=open(datadir+'20170214_untransformedR2_expression.txt').readline().split())
    exprdata.pop('probeset')
    exprdata.index = exprdata.pop('#H:hugo')
    #exprdata = pd.read_table(gzip.open(datadir+'core_exon_nrc283_zi.tsv.gz','rt'),skiprows=8)
    #exprdata.index = exprdata.pop('#H:probeset')
    #probeset2genes =  exprdata.pop('hugo')
    aCGH = pd.read_table(datadir+'20170214_untransformedR2_aCGH.txt',skiprows=44,
                         names=open(datadir+'20170214_untransformedR2_aCGH.txt').readline().split())
    aCGH.pop('probeset')
    aCGH.index = aCGH.pop('#H:hugo')
    return Dataset(**locals())

#TODO: FischerData needs to be transformed to use the RNA seq counts in GSE62564

@storeDatasetLocally
def get_FischerData():
    filteredOn = { #For reference, to know how dataset has been filtered
        'minimalSurvivalLastFollowup': 365*5
        }
    
    #Metadata
    metadata = pd.read_table(gzip.open(datadir+
                                       "R2_grabbed_data/Fischer498/metadata_src/GSE49710_series_matrix.txt.gz",'rt',
                                       encoding="UTF-8"),
                             skiprows=47,skipfooter=44799-66,engine='python',header=None)
    metadata.index = metadata[0].apply(lambda x: x.replace('!',''))
    del metadata[0]
    metadata = metadata.T
    i = count()
    metadata.columns = [c.replace('ch1',str(next(i))) if c.startswith('Sample_char') else c for c in metadata.columns]
    del i, metadata['Sample_source_name_ch1'], metadata['Sample_status'], metadata['Sample_organism_ch1']
    metadata.columns = [metadata[c][metadata.first_valid_index()].split(':')[0].replace(' ','_')
                        if c.startswith('Sample_char') else c for c in metadata.columns]
    metadata = metadata.applymap(lambda x: x.split(': ')[1] if ': ' in x else x)
    
    metadatasurv = pd.read_table(gzip.open(datadir+"R2_grabbed_data/Fischer498/metadata_src/GSE62564_series_matrix.txt.gz",'rt',encoding="UTF-8"),
                                 skiprows=51,skipfooter=102-74,engine='python')
    metadatasurv.index = metadatasurv['!Sample_geo_accession'].apply(lambda x: x.replace('!',''))
    del metadatasurv['!Sample_geo_accession']
    metadatasurv = metadatasurv.T
    metadatasurv.columns = range(len(metadatasurv.columns))
    metadatasurv = metadatasurv[list(range(7,21))]
    metadatasurv.columns = [v.split(':')[0].replace(' ','_') for v in metadatasurv.ix[metadatasurv.first_valid_index()]]
    metadatasurv = metadatasurv.applymap(lambda x: x.split(': ')[1] if not x is np.nan else x)
    metadatasurv.index = metadata.index #Both sample sets are sorted the same way, but different gse names
    assert sum(metadatasurv.age == metadata.age_at_diagnosis) == len(metadata)
    metadata['overall_survival'] = metadatasurv.os_day.apply(int)
    metadata['eventfree_survival'] = metadatasurv.efs_day.apply(int)
    del metadatasurv
    metadata.index = metadata.Sample_geo_accession.apply(lambda x: x.lower())
    metadata.Sample_title = metadata.Sample_title.apply(lambda x: x[5:].replace(' patient ',''))
    metadata.death_from_disease = metadata.death_from_disease == '1'
    metadata.progression = metadata.progression == '1'

    #Expression data
    exprdata = pd.read_table(datadir+'R2_grabbed_data/Fischer498/GSE49710_R2.txt')
    exprdata.index = exprdata.pop('#H:hugo')
    del exprdata['probeset']

    #aCGH
    aCGH = pd.read_table(datadir+'R2_grabbed_data/Fischer498/SEQC_aCGH/SEQC_aCGH_all_146.txt')
    geosearch = metadata[['Sample_title','Sample_geo_accession']].copy()
    geosearch.Sample_geo_accession = geosearch.index
    geosearch.index = geosearch.Sample_title
    aCGH.Sample = aCGH.Sample.apply(lambda x: geosearch.ix[x].Sample_geo_accession)
    del geosearch
    aCGH['log2ratio'] = (aCGH.CN/2).apply(np.log2)
    #Convert coordinates to hg38
    lo = LSD.get_lift19to38()
    aCGH['Start38'] = aCGH.T.apply(lambda x: lo.convert_coordinate(x.Chromosome,x.Start)).apply(lambda x: x[0][1] if x else np.nan)
    aCGH['End38'] = aCGH.T.apply(lambda x: lo.convert_coordinate(x.Chromosome,x.End)).apply(lambda x: x[0][1] if x else np.nan)
    del lo, aCGH['Start'], aCGH['End']
    aCGH = aCGH.dropna().copy()
    #Assign genes to regions
    genannot = LSD.get_ensemblGeneannot()
    aCGH['genes'] = aCGH.T.apply(lambda x: {f.attributes['gene_name'][0] for f in genannot.region('{}:{}-{}'
                                                    .format(x.Chromosome[3:],int(x.Start38),int(x.End38)),featuretype='gene')})
    aCGH['nrGenes'] = aCGH.genes.apply(len)
    del genannot
    #To set cut offs look at hist => aCGH.log2ratio.hist(ax=ax,bins='auto')
    aCGH['annotation'] = aCGH.log2ratio.apply(lambda x: 'gain' if x > 0.3 else ('loss' if x < -0.3 else 'normal'))

    # Filter patients whom according to metadata survived, but had last follow up before 5 years
    metadata = metadata[metadata.death_from_disease |
        (~metadata.death_from_disease &
        (metadata.overall_survival > filteredOn['minimalSurvivalLastFollowup']))]
    exprdata = exprdata[metadata.index]
    aCGH = aCGH[aCGH.Sample.isin(metadata.index)]
    
    return Dataset(**locals())

#TODO => FD sequencing data
#@storeDatasetLocally
def get_SequencedFischerData():
    """
    For the sequence data, the Trente TUC files were chosen, following
    'cufflinks on RefSeq genes and reported as log2(1 + FPKM)' see TUC readme

    Reference: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49711

    Source: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE49711&format=file&file=GSE49711%5FSEQC%5FNB%5FTUC%5FG%5Flog2%2Etxt%2Egz
    """
    filteredOn = { #For reference, to know how dataset has been filtered
        'minimalSurvivalLastFollowup': 365*5
        }
    print('Filtering settings:', filteredOn)

    #Metadata
    metadata = pd.read_table(gzip.open(datadir+
                                       "R2_grabbed_data/Fischer498/metadata_src/GSE49710_series_matrix.txt.gz",'rt',
                                       encoding="UTF-8"),
                             skiprows=47,skipfooter=44799-66,engine='python',header=None)
    metadata.index = metadata[0].apply(lambda x: x.replace('!',''))
    del metadata[0]
    metadata = metadata.T
    i = count()
    metadata.columns = [c.replace('ch1',str(next(i))) if c.startswith('Sample_char') else c for c in metadata.columns]
    del i, metadata['Sample_source_name_ch1'], metadata['Sample_status'], metadata['Sample_organism_ch1']
    metadata.columns = [metadata[c][metadata.first_valid_index()].split(':')[0].replace(' ','_')
                        if c.startswith('Sample_char') else c for c in metadata.columns]
    metadata = metadata.applymap(lambda x: x.split(': ')[1] if ': ' in x else x)
    
    metadatasurv = pd.read_table(gzip.open(datadir+"R2_grabbed_data/Fischer498/metadata_src/GSE62564_series_matrix.txt.gz",'rt',encoding="UTF-8"),
                                 skiprows=51,skipfooter=102-74,engine='python')
    metadatasurv.index = metadatasurv['!Sample_geo_accession'].apply(lambda x: x.replace('!',''))
    del metadatasurv['!Sample_geo_accession']
    metadatasurv = metadatasurv.T
    metadatasurv.columns = range(len(metadatasurv.columns))
    metadatasurv = metadatasurv[list(range(7,21))]
    metadatasurv.columns = [v.split(':')[0].replace(' ','_') for v in metadatasurv.ix[metadatasurv.first_valid_index()]]
    metadatasurv = metadatasurv.applymap(lambda x: x.split(': ')[1] if not x is np.nan else x)
    metadatasurv.index = metadata.index #Both sample sets are sorted the same way, but different gse names
    assert sum(metadatasurv.age == metadata.age_at_diagnosis) == len(metadata)
    metadata['overall_survival'] = metadatasurv.os_day.apply(int)
    metadata['eventfree_survival'] = metadatasurv.efs_day.apply(int)
    del metadatasurv
    metadata.Sample_geo_accession = metadata.Sample_geo_accession.apply(lambda x: x.lower())
    metadata.set_index("Sample_geo_accession",inplace=True)
    metadata.Sample_title = metadata.Sample_title.apply(lambda x: x[5:].replace(' patient ',''))
    metadata.death_from_disease = metadata.death_from_disease == '1'
    metadata.progression = metadata.progression == '1'

    #Expression data
    ## Array expression
    exprdata_A = pd.read_table(datadir+'R2_grabbed_data/Fischer498/GSE49710_R2.txt')
    exprdata_A.index = exprdata_A.pop('#H:hugo')
    del exprdata_A['probeset']
    
    ## RNAseq gene expression
    exprdata_G = pd.read_table(
        gzip.open(
            datadir+'R2_grabbed_data/Fischer498/sequencedData/GSE49711_SEQC_NB_TUC_G_log2.txt.gz',
            'rt',encoding="UTF-8"
        ),
        index_col = '00gene_id'
    )
    exprdata_G.columns = [
        metadata.reset_index().set_index('Sample_title').ix[c.split('_')[1]].Sample_geo_accession
        for c in exprdata_G.columns
    ]

    ## RNAseq transcript level expression
    exprdata_T = pd.read_table(
        gzip.open(
            datadir+'R2_grabbed_data/Fischer498/sequencedData/GSE49711_SEQC_NB_TUC_T_log2.txt.gz',
            'rt',encoding="UTF-8"
        ),
        index_col = '00transcript_id'
    )
    exprdata_T.columns = exprdata_G.columns

    ## RNAseq junction level expression
    exprdata_J = pd.read_table(
        gzip.open(
            datadir+'R2_grabbed_data/Fischer498/sequencedData/GSE49711_SEQC_NB_TUC_J_log2.txt.gz',
            'rt',encoding="UTF-8"
        ),
        index_col = 'sample_ID'
    )
    exprdata_J.columns = exprdata_G.columns

    # Default expression dataset -> exprdata_G
    exprdata = exprdata_G
    
    #aCGH
    aCGH = pd.read_table(datadir+'R2_grabbed_data/Fischer498/SEQC_aCGH/SEQC_aCGH_all_146.txt')
    geosearch = metadata[['Sample_title']].copy()
    geosearch.reset_index(inplace=True)
    geosearch.set_index('Sample_title',inplace=True)
    aCGH.Sample = aCGH.Sample.apply(lambda x: geosearch.ix[x].Sample_geo_accession)
    del geosearch
    aCGH['log2ratio'] = (aCGH.CN/2).apply(np.log2)
    #Convert coordinates to hg38
    lo = LSD.get_lift19to38()
    aCGH['Start38'] = aCGH.T.apply(lambda x: lo.convert_coordinate(x.Chromosome,x.Start)).apply(lambda x: x[0][1] if x else np.nan)
    aCGH['End38'] = aCGH.T.apply(lambda x: lo.convert_coordinate(x.Chromosome,x.End)).apply(lambda x: x[0][1] if x else np.nan)
    del lo, aCGH['Start'], aCGH['End']
    aCGH = aCGH.dropna().copy()
    #Assign genes to regions
    genannot = LSD.get_ensemblGeneannot()
    aCGH['genes'] = aCGH.T.apply(lambda x: {f.attributes['gene_name'][0] for f in genannot.region('{}:{}-{}'
                                                    .format(x.Chromosome[3:],int(x.Start38),int(x.End38)),featuretype='gene')})
    aCGH['nrGenes'] = aCGH.genes.apply(len)
    del genannot
    #To set cut offs look at hist => aCGH.log2ratio.hist(ax=ax,bins='auto')
    aCGH['annotation'] = aCGH.log2ratio.apply(lambda x: 'gain' if x > 0.3 else ('loss' if x < -0.3 else 'normal'))

    # Filter patients whom according to metadata survived, but had last follow up before 5 years
    metadata = metadata[metadata.death_from_disease |
        (~metadata.death_from_disease &
        (metadata.overall_survival > filteredOn['minimalSurvivalLastFollowup']))]
    exprdata_A = exprdata_A[metadata.index]
    exprdata_G = exprdata_G[metadata.index]
    exprdata_T = exprdata_T[metadata.index]
    exprdata_J = exprdata_J[metadata.index]
    aCGH = aCGH[aCGH.Sample.isin(metadata.index)]
    
    return Dataset(**locals())
