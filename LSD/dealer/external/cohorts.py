from LSD import storeDatasetLocally, datadir, get_lift19to38, get_ensemblGeneannot, Dataset
from os.path import expanduser, exists
import gzip
import pandas as pd, numpy as np
from itertools import count

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
    lo = get_lift19to38()
    aCGH['Start38'] = aCGH.T.apply(lambda x: lo.convert_coordinate(x.Chromosome,x.Start)).apply(lambda x: x[0][1] if x else np.nan)
    aCGH['End38'] = aCGH.T.apply(lambda x: lo.convert_coordinate(x.Chromosome,x.End)).apply(lambda x: x[0][1] if x else np.nan)
    del lo, aCGH['Start'], aCGH['End']
    aCGH = aCGH.dropna().copy()
    #Assign genes to regions
    genannot = get_ensemblGeneannot()
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
    exprdata = [metadata.index]
    aCGH = aCGH[aCGH.Sample.isin(metadata.index)]
    
    return Dataset(**locals())
