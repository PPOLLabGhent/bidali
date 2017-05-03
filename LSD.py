#!/usr/bin/env python
# LSD: Lab Speleman Datasets
# Module that preprocesses commonly used datasets at Speleman Lab
# and makes them available for use in python.

import gzip, pickle
from zipfile import ZipFile
from io import TextIOWrapper
import pandas as pd, numpy as np
from os.path import expanduser, exists
from itertools import count

## Defaults
datadir = expanduser('~/Dropbox (speleman lab)/Lab/z_archive/Datasets/')
processedDataStorage = expanduser('~/Data/LSDpy/')

## Utility functions
def getLSDataset(name,**kwargs):
    try: return name(**kwargs)
    except TypeError:
        return globals()[name](**kwargs)

def listLSDatasets():
    return [i for i in sorted(globals()) if i.startswith('get_')]

class Dataset:
    """
    A Dataset object is a collection of datasets, accessible as
    attributes from the Dataset object.

    to_R pushes the sub datasets to R making them available in the
    global namespace
    """
    def __init__(self,**kwargs):
        self.__datasets__ = set(kwargs)
        for kw in kwargs:
            self.__setattr__(kw,kwargs[kw])

    def to_R(self):
        pass


def retrieveSources(dataset_getfunction):
    """
    A dataset_getfunction function that contains 'Source:' lines
    in the docstring, can be decorated with this function.
    If a source is not locally available, it will be downloaded
    and added to the processedDataStorage location.

    A source line has to be formatted accordingly:
    Source: [filename] url

    If filename is not provided, the last part of the url (after last '/')
    is taken as filename.
    """
    import inspect, requests

    def wrapper(*args, **kwargs):
        try:
            return dataset_getfunction(*args, **kwargs)
        except FileNotFoundError:
            for docline in inspect.getdoc(dataset_getfunction).split('\n'):
                if docline.startswith('Source:'):
                    docline = docline.split()
                    if len(docline) == 2:
                        url = docline[1]
                        filename = url[url.rindex('/')+1:]
                    elif len(docline) == 3:
                        url = docline[2]
                        filename = docline[1]
                if not exists(processedDataStorage+filename):
                    r = requests.get(url)
                    with open(processedDataStorage+filename,'wb') as f:
                        f.write(r.content)
            try: return dataset_getfunction(*args, **kwargs)
            except FileNotFoundError:
                print('Either not all source files are documented correctly in docstring,',
                      'or there is a source file unrelated issue')
                raise

    return wrapper

def storeDatasetLocally(dataset_getfunction):
    """
    Can be used as a decorator for 'get_dataset' functions.
    It will check if a processed dataset is locally available,
    and if so, load that one instead of processing from the source
    files.

    Should only be used for functions that do not process the data
    differently depending on the 'get_dataset' function arguments.
    """
    import inspect, hashlib
    
    def wrapper(*args, **kwargs):
        #Check if data was already processed
        functionSource = inspect.getsource(dataset_getfunction)
        hashvalue = hashlib.md5(functionSource.encode()).hexdigest()
        datastorage = '{}{}_{}.pickle'.format(processedDataStorage,
                                              dataset_getfunction.__name__.replace('get_',''),
                                              hashvalue)
        if exists(datastorage):
            dataset = pickle.load(open(datastorage,'rb'))
            return dataset
        else:
            dataset = dataset_getfunction(*args, **kwargs)
            try: pickle.dump(dataset,open(datastorage,'wb'))
            except FileNotFoundError:
                print('Not possible to store dataset locally. Create',processedDataStorage,
                      'if you want to avoid reprocessing dataset on every call.')
            return dataset

    return wrapper

## Datasets
### References/annotations
def get_ensembl(onlyGeneLabeled=True,onlyInChromosomes=None):
    ensembl = pd.read_table(datadir+'Genomes/Ensembl/Biomart/idmapping_extended.txt',
                        header=None,index_col=0,names=('egid','etid','start','stop','gcC','chr','strand','TSS','typeg','typet','gene_label','entrez'))
    ensembl = ensembl[~ensembl.index.duplicated()]
    if onlyGeneLabeled: ensembl = ensembl[ensembl.gene_label.isnull().apply(lambda x: not(x))]
    ensembl.chr = ensembl.chr.apply(lambda x: 'chr'+x.replace('MT','M'))
    if onlyInChromosomes: ensembl = ensembl[ensembl.chr.isin(onlyInChromosomes)]
    return ensembl

@retrieveSources
def get_ensemblGeneannot():
    """
    Info: http://www.ensembl.org/info/data/ftp/index.html
    Source: ftp://ftp.ensembl.org/pub/release-88/gtf/homo_sapiens/Homo_sapiens.GRCh38.88.gtf.gz
    """
    import gffutils
    try: db = gffutils.FeatureDB(processedDataStorage+'Homo_sapiens.GRCh38.86.sqlite3')
    except ValueError:
        if not exists(processedDataStorage+'Homo_sapiens.GRCh38.86.gtf.gz'):
            raise FileNotFoundError
        db = gffutils.create_db(processedDataStorage+'Homo_sapiens.GRCh38.86.gtf.gz',
                                processedDataStorage+'Homo_sapiens.GRCh38.86.sqlite3',
                                disable_infer_genes=True,disable_infer_transcripts=True)
    return db
    
def get_entrez():
    """
    Source: ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/RefSeqGene/gene_RefSeqGene
    """
    entrez = pd.read_table('Dropbiz/Lab/z_archive/Datasets/Genomes/Entrez/gene_RefSeqGene', index_col='GeneID')
    return entrez

def get_liftover(frm=19,to=38):
    """
    Info: http://hgdownload.cse.ucsc.edu/downloads.html
    Source: http://hgdownload.cse.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz
    """
    from pyliftover import LiftOver
    liftoverfile = 'hg{}ToHg{}.over.chain.gz'.format(frm,to)
    try: lo = LiftOver(processedDataStorage+liftoverfile)
    except FileNotFoundError:
        import requests
        liftoverurl = 'http://hgdownload.cse.ucsc.edu/gbdb/hg{}/liftOver/{}'.format(frm,liftoverfile)
        r = requests.get(liftoverurl)
        with open(processedDataStorage+liftoverfile,'wb') as file:
            file.write(r.content)
        lo = LiftOver(processedDataStorage+liftoverfile)
    return lo

def get_lift19to38():
    return get_liftover(frm=19,to=38)

@storeDatasetLocally
def get_proteinNetworks():
    """
    Sources: 
      https://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.4.147/BIOGRID-ALL-3.4.147.tab2.zip
      http://string-db.org/download/protein.links.v10/9606.protein.links.v10.txt.gz
      http://string-db.org/mapping_files/entrez_mappings/entrez_gene_id.vs.string.v10.28042015.tsv
    """
    import networkx as nx

    #Biogrid
    with ZipFile(datadir+'ProteinNetworks/BIOGRID-ALL-3.4.147.tab2.zip') as biogridzip:
        ds = pd.read_table(TextIOWrapper(biogridzip.open('BIOGRID-ALL-3.4.147.tab2.txt','r')),low_memory=False)
    ds = ds[ds['Organism Interactor A'] == 9606]
    Gbio = nx.Graph()
    ds.T.apply(lambda x: Gbio.add_edge(x['Official Symbol Interactor A'],x['Official Symbol Interactor B']))

    #String-DB
    stringdb = pd.read_table(gzip.open(datadir+'ProteinNetworks/9606.protein.links.v10.txt.gz','rt'),sep=' ')
    stringids = pd.read_table(datadir+'ProteinNetworks/entrez_gene_id.vs.string.v10.28042015.tsv',index_col='STRING_Locus_ID')
    entrez = get_entrez()
    stringids = stringids[stringids['#Entrez_Gene_ID'].isin(entrez.index)]
    stringdb = stringdb[stringdb.combined_score > 400] #study stringdb.combined_score.hist(bins='auto') to set threshold
    stringdb = stringdb[stringdb.protein1.isin(stringids.index) & stringdb.protein2.isin(stringids.index)]
    stringdb.protein1 = stringdb.protein1.apply(lambda x: stringids.ix[x]['#Entrez_Gene_ID'])
    stringdb.protein2 = stringdb.protein2.apply(lambda x: stringids.ix[x]['#Entrez_Gene_ID'])
    stringdb.protein1 = stringdb.protein1.apply(lambda x: entrez.ix[x].Symbol)
    stringdb.protein2 = stringdb.protein2.apply(lambda x: entrez.ix[x].Symbol)
    Gstring = nx.Graph()
    stringdb.T.apply(lambda x: Gstring.add_edge(x.protein1,x.protein2))
    
    return Dataset(biogridnx = Gbio, biogrid = ds,
                   stringnx = Gstring, string = stringdb)
    
@storeDatasetLocally
def get_centromeres():
    """
    Source: from R bioconductor GWASTools: data(centromeres.hg38)
    """
    from seqanalysis import loadHumanGenome

    # Centromere positions
    centromereshg38="""1      1 122026460  125184587
2      2  92188146   94090557
3      3  90772459   93655574
4      4  49708101   51743951
5      5  46485901   50059807
6      6  58553889   59829934
7      7  58169654   60828234
8      8  44033745   45877265
9      9  43236168   45518558
10    10  39686683   41593521
11    11  51078349   54425074
12    12  34769408   37185252
13    13  16000001   18051248
14    14  16000001   18173523
15    15  17000001   19725254
16    16  36311159   38280682
17    17  22813680   26885980
18    18  15460900   20861206
19    19  24498981   27190874
20    20  26436233   30038348
21    21  10864561   12915808
22    22  12954789   15054318
X      X  58605580   62412542
Y      Y  10316945   10544039"""

    centromereshg38 = pd.DataFrame([c.split()[-3:] for c in centromereshg38.split('\n')],
                                   columns= "chrom left_base right_base".split())
    centromereshg38.index = centromereshg38.chrom.apply(lambda x: 'chr'+x)
    centromereshg38['left_base'] = centromereshg38.pop('left_base').apply(int)
    centromereshg38['right_base']=centromereshg38.pop('right_base').apply(int)

    genome = loadHumanGenome()
    centromereshg38['len'] = centromereshg38.apply(lambda x: len(genome.chromosomes[x.name]),axis=1)
    centromereshg38['qlen'] = centromereshg38.len - centromereshg38.right_base
    centromereshg38['chr_weight'] = centromereshg38.len/centromereshg38.len.max()
    centromereshg38['q_weight'] = centromereshg38.qlen/centromereshg38[['left_base','qlen']].max().max()
    centromereshg38['p_weight'] = centromereshg38.left_base/centromereshg38[['left_base','qlen']].max().max()
    del genome

    ensembl = get_ensembl()
    ensembl = ensembl[ensembl.chr.isin(centromereshg38.index)]
    ensembl['chrarm'] = ensembl.apply(lambda x: 'p' if x.stop < centromereshg38.ix[x.chr].left_base else
                                  ('q' if x.start > centromereshg38.ix[x.chr].right_base else 'pq'),axis=1)
    centromereshg38['chr_genes'] = ensembl.groupby('chr').size()
    centromereshg38['p_genes'] = ensembl[ensembl.chrarm=='p'].groupby('chr').size()
    centromereshg38['q_genes'] = ensembl[ensembl.chrarm=='q'].groupby('chr').size()
    centromereshg38['chr_gweight'] = centromereshg38.chr_genes/centromereshg38.chr_genes.max()
    centromereshg38['q_gweight'] = centromereshg38.q_genes/centromereshg38[['p_genes','q_genes']].max().max()
    centromereshg38['p_gweight'] = centromereshg38.p_genes/centromereshg38[['p_genes','q_genes']].max().max()
    del ensembl

    return centromereshg38

### Cohorts
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
    metadata['overall_survival'] = metadatasurv.os_day
    metadata['eventfree_survival'] = metadatasurv.efs_day
    del metadatasurv
    metadata.index = metadata.Sample_geo_accession.apply(lambda x: x.lower())
    metadata.Sample_title = metadata.Sample_title.apply(lambda x: x[5:].replace(' patient ',''))
    metadata.overall_survival = metadata.overall_survival.apply(int)
    metadata.eventfree_survival = metadata.eventfree_survival.apply(int)
    metadata.death_from_disease = metadata.death_from_disease == '1'
    metadata.progression = metadata.progression == '1'

    #Expression data
    exprdata = pd.read_table(datadir+'R2_grabbed_data/Fischer498/GSE49710_R2.txt')
    exprdata.index = exprdata.pop('#H:hugo')
    del exprdata['probeset']

    #aCGH
    from pyliftover import LiftOver
    import gffutils
    lo = get_lift19to38()
    genannot = gffutils.FeatureDB('/home/christophe/Data/Genomes/Homo_sapiens.GRCh38.86.sqlite3')
    aCGH = pd.read_table(datadir+'R2_grabbed_data/Fischer498/SEQC_aCGH/SEQC_aCGH_all_146.txt')
    geosearch = metadata[['Sample_title','Sample_geo_accession']].copy()
    geosearch.Sample_geo_accession = geosearch.index
    geosearch.index = geosearch.Sample_title
    aCGH.Sample = aCGH.Sample.apply(lambda x: geosearch.ix[x].Sample_geo_accession)
    del geosearch
    aCGH['log2ratio'] = (aCGH.CN/2).apply(np.log2)
    #Convert coordinates to hg38
    aCGH['Start38'] = aCGH.T.apply(lambda x: lo.convert_coordinate(x.Chromosome,x.Start)).apply(lambda x: x[0][1] if x else np.nan)
    aCGH['End38'] = aCGH.T.apply(lambda x: lo.convert_coordinate(x.Chromosome,x.End)).apply(lambda x: x[0][1] if x else np.nan)
    del aCGH['Start'], aCGH['End']
    aCGH = aCGH.dropna()
    #Assign genes to regions
    aCGH['genes'] = aCGH.T.apply(lambda x: {f.attributes['gene_name'][0] for f in genannot.region('{}:{}-{}'
                                                .format(x.Chromosome,int(x.Start38),int(x.End38))) if f.featuretype == 'gene'})
    aCGH['nrGenes'] = aCGH.genes.apply(len)
    #To set cut offs look at hist => aCGH.log2ratio.hist(ax=ax,bins='auto')
    aCGH['annotation'] = aCGH.log2ratio.apply(lambda x: 'gain' if x > 0.3 else ('loss' if x < -0.3 else 'normal'))
    
    return Dataset(**locals())
