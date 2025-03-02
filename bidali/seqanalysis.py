#!/bin/env python3
#CVN set of functions for sequence analysis
#General imports
import os
import numpy as np
import pandas as pd
from lostdata.processing import storeDatasetLocally, retrieveSources, processedDataStorage
from lostdata.dealer.ensembl import get_ensemblGeneannot as get_ensembl

def recomplement(dna):
    """
    Returns the complement of a DNA motive for regex searching
    """
    return dna.translate(recomplement.dict)[::-1]
try: recomplement.dict=str.maketrans('ACTGactg()[]{}','TGACtgac)(][}{')
except AttributeError: #py2#
    import string
    recomplement.dict=string.maketrans('ACTGactg()[]{}','TGACtgac)(][}{')

class DNA:
    def __init__(self,name,sequence):
        self.name = name
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, key):
        return self.sequence[key]

    def __str__(self):
        return self.name
    
    def __repr__(self):
        return ('{} ({}...)'.format(self.name,self.sequence[:10]))

class DNAregion:
    def __init__(self,dna,region):
        self.dna = dna
        self.region = region

    def __str__(self):
        return self.dna[region]
    
    def __repr__(self):
        return self.dna[region]

class Genome:
    """
    A simple class to represent a genome.
    """
    def __init__(self,species,assembly,chromosomes):
        self.species = species
        self.assembly = assembly
        self.chromosomes = {}
        for ch in chromosomes:
            self.addChromosome(ch.name,ch)

    def addChromosome(self,name,sequence):
        self.chromosomes[name] = sequence

    def __len__(self):
        return len(self.chromosomes)

    def __getitem__(self, key):
        if key[1].step and key[1].step < 0:
            return recomplement(self.chromosomes[key[0]][
                key[1].start-1:key[1].stop:abs(key[1].step)])
        else:
            return self.chromosomes[key[0]][key[1].start-1:key[1].stop:key[1].step]

    def windowSlider(self,windowSize=100,overlapping=True):
        for ch in sorted(self.chromosomes):
            for i in range(0,len(self.chromosomes[ch]),1 if overlapping else windowSize):
                yield (ch,slice(i,i+windowSize))

@retrieveSources
def loadHumanGenome():
    """
    Loads the GRCh38 human genome

    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.2.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.3.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.4.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.5.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.6.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.8.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.9.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.10.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.12.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.14.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.15.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.16.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.17.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.18.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
    Source: ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz

    """
    from glob import glob
    import gzip
    files = glob(os.path.join(processedDataStorage, 'Homo_sapiens.GRCh38.dna.chromosome.*'))
    if not files: raise FileNotFoundError #TODO raise when amount of files is not as expected
    chromosomes = []
    for f in files:
        with gzip.open(f,mode='rt') as fh:
            f = fh.readlines()
        chromosomes.append(
            DNA(
                'chr'+f.pop(0).strip()[1:].split()[0],
                ''.join([l.strip() for l in f])
            )
        )
    return Genome('human','GRCh38',chromosomes)

@storeDatasetLocally
def get_centromeres():
    """
    Source: from R bioconductor GWASTools: data(centromeres.hg38)
    """
    from bidali.seqanalysis import loadHumanGenome

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
    ensembl = pd.DataFrame(
        [(g.id, 'chr'+g.chrom, g.start, g.stop, g.strand, g.attributes['gene_name'][0])
            for g in ensembl.features_of_type('gene')],
        columns=('id','chr','start','stop','strand','name')
    )
    ensembl = ensembl[ensembl.chr.isin(centromereshg38.index)]
    ensembl['chrarm'] = ensembl.apply(lambda x: 'p' if x.stop < centromereshg38.loc[x.chr].left_base else
                                  ('q' if x.start > centromereshg38.loc[x.chr].right_base else 'pq'),axis=1)
    centromereshg38['chr_genes'] = ensembl.groupby('chr').size()
    centromereshg38['p_genes'] = ensembl[ensembl.chrarm=='p'].groupby('chr').size()
    centromereshg38['q_genes'] = ensembl[ensembl.chrarm=='q'].groupby('chr').size()
    centromereshg38['chr_gweight'] = centromereshg38.chr_genes/centromereshg38.chr_genes.max()
    centromereshg38['q_gweight'] = centromereshg38.q_genes/centromereshg38[['p_genes','q_genes']].max().max()
    centromereshg38['p_gweight'] = centromereshg38.p_genes/centromereshg38[['p_genes','q_genes']].max().max()
    del ensembl

    return centromereshg38


def get_centromeres_hg19():
    """
    Source: hg19 centromere positions.
    """
    from bidali.seqanalysis import loadHumanGenome
    
    # Centromere positions for hg19 (as per UCSC Genome Browser)
    centromereshg19 = """1      1  121535434  122227999
2      2  92280923   93235083
3      3  90624140   93645123
4      4  49653707   51907379
5      5  45461453   49233429
6      6  58550267   59672769
7      7  58499891   60930441
8      8  43861953   45937043
9      9  43347223   45372555
10    10  39552494   41455798
11    11  51074379   54259268
12    12  34842886   37310189
13    13  16000001   18046601
14    14  16000001   18174892
15    15  17000001   19778853
16    16  36280611   38347607
17    17  22802406   26792423
18    18  15472766   20853106
19    19  24464388   27053569
20    20  26473426   29924481
21    21  10864396   12911260
22    22  12945960   15062510"""
    
    # Parse centromere data into a DataFrame
    centromereshg19 = pd.DataFrame([c.split()[-3:] for c in centromereshg19.split('\n')],
                                   columns="chrom left_base right_base".split())
    centromereshg19.index = centromereshg19.chrom.apply(lambda x: 'chr'+x)
    centromereshg19['left_base'] = centromereshg19.pop('left_base').apply(int)
    centromereshg19['right_base'] = centromereshg19.pop('right_base').apply(int)
    
    # Load human genome to get chromosome lengths
    genome = loadHumanGenome()
    
    # Calculate additional centromere properties
    centromereshg19['len'] = centromereshg19.apply(lambda x: len(genome.chromosomes[x.name]), axis=1)
    centromereshg19['qlen'] = centromereshg19.len - centromereshg19.right_base
    centromereshg19['chr_weight'] = centromereshg19.len / centromereshg19.len.max()
    centromereshg19['q_weight'] = centromereshg19.qlen / centromereshg19[['left_base', 'qlen']].max().max()
    centromereshg19['p_weight'] = centromereshg19.left_base / centromereshg19[['left_base', 'qlen']].max().max()
    
    # Clean up genome to free memory
    del genome
    
    # Load Ensembl gene data
    ensembl = get_ensembl()
    ensembl = pd.DataFrame(
        [(g.id, 'chr' + g.chrom, g.start, g.stop, g.strand, g.attributes['gene_name'][0])
         for g in ensembl.features_of_type('gene')],
        columns=('id', 'chr', 'start', 'stop', 'strand', 'name')
    )
    
    # Filter Ensembl genes that intersect with centromere data
    ensembl = ensembl[ensembl.chr.isin(centromereshg19.index)]
    
    # Classify genes by chromosomal arm (p, q, or pq)
    ensembl['chrarm'] = ensembl.apply(lambda x: 'p' if x.stop < centromereshg19.loc[x.chr].left_base else
                                      ('q' if x.start > centromereshg19.loc[x.chr].right_base else 'pq'), axis=1)
    
    # Calculate number of genes per chromosome and arm
    centromereshg19['chr_genes'] = ensembl.groupby('chr').size()
    centromereshg19['p_genes'] = ensembl[ensembl.chrarm == 'p'].groupby('chr').size()
    centromereshg19['q_genes'] = ensembl[ensembl.chrarm == 'q'].groupby('chr').size()
    
    # Weighting based on the number of genes
    centromereshg19['chr_gweight'] = centromereshg19.chr_genes / centromereshg19.chr_genes.max()
    centromereshg19['q_gweight'] = centromereshg19.q_genes / centromereshg19[['p_genes', 'q_genes']].max().max()
    centromereshg19['p_gweight'] = centromereshg19.p_genes / centromereshg19[['p_genes', 'q_genes']].max().max()
    
    # Clean up Ensembl to free memory
    del ensembl
    
    return centromereshg19


class PFM():
    def __init__(self,pfmfile):
        pfm = open(pfmfile).readlines()
        self.name = pfm.pop(0).strip().split().pop()
        pfm = [l.replace('[','').replace(']','').strip().split() for l in pfm]
        self.nucleotides = np.array([l.pop(0) for l in pfm])
        self.pfm = np.array([[int(i) for i in l] for l in pfm])
        self.pfm_norm = (self.pfm/self.pfm.sum(0)).round(3)

    def setMinOfMax(self,minimal=0.30):
        self.includeInMotive = (self.pfm/self.pfm.max(0))>minimal

    def generateMotive(self,minimal=0.30):
        self.setMinOfMax(minimal)
        return (r'['+
                r']['.join([''.join(self.nucleotides[self.includeInMotive[:,i]]) for i in range(self.includeInMotive.shape[1])])+
                r']')

# Gene annotation functions
def literatureLinkSearch(term,referenceTerm='quadruplex'):
    """
    Returns all pubmed ids for pubs containing both term and referenceTerm
    """
    import requests
    import xml.etree.ElementTree as ET
    try: r = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',
                          params={'db':'pubmed','term':'{}[Title/Abstract] AND {}[Title/Abstract]'.format(referenceTerm,term)})
    except Exception:
        return np.nan
    root = ET.fromstring(r.text)
    return [i.text for i in root.findall('IdList/Id')]

# Signature enrichment test
def calcSignature(counts,up,down=None,annotation=None,ax=None,xax=None,grouping=None,g4sig=False,
                  medianLines=True,cortests=True,**kwargs):
    """
    - signature (up,down) should be provided as (a) pandas.Series
    - grouping should be of the form (grouping_func,markers_dict) -> e.g. (lambda x: x.split('_')[0],{'SH-SY5Y':'s','IMR32':'o'})
            grouping_func should provide as result all markers_dict keys
    - xax if used should be pd.Series with index corresponding to counts.columns
    """
    import pandas as pd
    if type(up) == str:
        if not 'signatures' in dir(testSignature):
            testSignature.signatures = pd.read_table('/home/christophe/Dropbiz/Lab/z_archive/Datasets/dnupGenesetsOfInterest.gmt',
                                                     index_col=0,header=None)
            testSignature.sigdescriptions = testSignature.signatures.pop(1)
        up = testSignature.signatures.ix[up].dropna()
    if annotation is not None:
        up = annotation[annotation.gene_label.isin(up)].index
    ranks = counts.rank()
    rankSum = ranks[ranks.index.isin(up)].sum()
    if g4sig:
        rankSum = calcGlobalG4sig(ranks,ensembl,rank=False) if type(g4sig) == bool else g4sig
    xax = xax.sort_values() if xax is not None else pd.Series(range(len(rankSum)),rankSum.index)
    if ax:
        rankRel = rankSum/rankSum.max()
        if grouping and down is None:
            for gname,grp in rankRel.groupby(grouping[0]):
                ax.scatter([xax.ix[i] for i in grp.index],grp,c='r',s=30,marker=grouping[1][gname])
        else: ax.scatter([xax.ix[i] for i in rankRel.index],rankRel,c='r',label='up ({})'.format(len(up)),**kwargs)
        ax.set_xticks(xax)
        ax.set_xticklabels(xax.index,rotation=-30,ha='left')
        padding = (max(xax)-min(xax))*0.01
        ax.set_xlim((min(xax)-padding,max(xax)+padding))
        if medianLines:
            ax.axvline(xax.median())
            ax.axhline(rankRel.median())
        if cortests:
            from scipy.stats import spearmanr,pearsonr,fisher_exact
            spearesults = spearmanr([xax.ix[i] for i in rankRel.index],rankRel)
            ax.text(.99,.99,
                    'Spearman: {:.2f} ({:.4f})\nPearson: {:.2f} ({:.4f})\nFisherET: {:.2f} ({:.4f})'.format(
                        *spearesults,
                        *pearsonr([xax.ix[i] for i in rankRel.index],rankRel),
                        *fisher_exact([[sum((xax<xax.median())&(rankRel>rankRel.median())),
                                        sum((xax>xax.median())&(rankRel>rankRel.median()))],
                                       [sum((xax<xax.median())&(rankRel<rankRel.median())),
                                        sum((xax>xax.median())&(rankRel<rankRel.median()))]],
                                      alternative = 'less' if spearesults[0]<0 else 'greater')
                    ),
                    transform=ax.transAxes,va='top',ha='right')
    if down is not None:
        if type(down) == str:
            down = testSignature.signatures.ix[down].dropna()
        if annotation is not None:
            down = annotation[annotation.gene_label.isin(down)].index
        ranksDown = counts.rank(ascending=False)
        rankDownSum = ranksDown[ranksDown.index.isin(down)].sum()
        rankTotalSum = pd.DataFrame({'up':rankSum,'down':rankDownSum}).mean(axis=1)
        if ax:
            ax.scatter([xax.ix[i] for i in rankDownSum.index],rankDownSum/rankDownSum.max(),c='b',label='down ({})'.format(len(down)))
            if grouping:
                for gname,grp in (rankTotalSum/rankTotalSum.max()).groupby(grouping[0]):
                    ax.scatter([xax.ix[i] for i in grp.index],grp,c='g',s=40,marker=grouping[1][gname])
            else: ax.scatter([xax.ix[i] for i in rankTotalSum.index],rankTotalSum/rankTotalSum.max(),c='g',s=40,marker='s',label='mean')
            ax.legend()
        return rankTotalSum/rankTotalSum.max()
    else:
        return rankSum/rankSum.max()

def calcGlobalG4sig(countRanks,geneG4annotation,colG4='G4s',rank=True):
    """
    If rank=True, assume counts were given and rank them
    Scores like this -> sum(log(geneRank**G4s)) == sum(G4s*log(geneRank))
    """
    countRanks = countRanks[countRanks.index.isin(geneG4annotation.index)]
    if rank: countRanks = countRanks.rank()
    #countRanks.apply(lambda x: x.apply(lambda y,x=x: y**int(ensembl.ix[x.name].G4s)),axis=1).applymap(np.log).sum()
    return countRanks.applymap(
        np.log).apply(
            lambda x: x.apply(
                lambda y,x=x: int(geneG4annotation.ix[x.name][colG4])*y),axis=1).sum()    
    
    
# In main section below, some G4 'nifty' programming -> should be moved to dedicated G4 research script
if __name__ == '__main__':
    genome = loadHumanGenome()
    import re
    import matplotlib.pylab as plt
    q4m = r'.{1,7}'.join((r'G{3,7}' for g in range(2)))
    q4m_AT = q4m.replace('.','[AT]')
    q4motif = re.compile(q4m, re.IGNORECASE)
    q4motif_complement = re.compile(q4m.replace('G','C'), re.IGNORECASE)
    
    allowG4bulges = False
    if allowG4bulges:
        q4motif,q4motif_complement = (re.compile(r'|'.join([q4motif.pattern[:i*12]+
                                                            bulge+q4motif.pattern[i*12+6:]
                                                            for bulge in (r'GG.G',r'G.GG')
                                                            for i in range(4)]) , re.IGNORECASE)
                                      for q4motif in (q4motif,q4motif_complement))
    
    import random
    windowSize=1000
    binSize=50
    includeRandomizedWindows = False
    ATcontent = []
    G4content = []
    rG4content = []
    for s in genome.windowSlider(windowSize=windowSize,overlapping=False):
        ATcontent.append(genome[s].upper().count('A')+genome[s].upper().count('T'))
        G4content.append(len(q4motif.findall(genome[s]))+len(q4motif_complement.findall(genome[s])))
        if includeRandomizedWindows:
            seq = list(genome[s])
            random.shuffle(seq)
            seq = ''.join(seq)
            rG4content.append(len(q4motif.findall(seq))+len(q4motif_complement.findall(seq)))
            
    G4perATcontent = [[] for i in range(windowSize+1)]
    for a,g in zip(ATcontent,G4content):
        G4perATcontent[a].append(g)
    if includeRandomizedWindows:
        rG4perATcontent = [[] for i in range(windowSize+1)]
        for a,g in zip(ATcontent,rG4content):
            rG4perATcontent[a].append(g)
        
    #G4perATcontent_filtered = [l for l in G4perATcontent if len(l) > 30]
    #plt.violinplot(G4perATcontent_filtered[::15])
    G4perATbin=[]
    ticks = [100*(bin+(bin+binSize))/(2*windowSize) for bin in range(binSize,windowSize,binSize)]
    for bin in range(binSize,windowSize,binSize):
        l = []
        for i in G4perATcontent[bin:bin+binSize]: l+=i
        G4perATbin.append(l if l else [0])
    if includeRandomizedWindows:
        rG4perATbin=[]
        for bin in range(binSize,windowSize,binSize):
            l = []
            for i in rG4perATcontent[bin:bin+binSize]: l+=i
            rG4perATbin.append(l if l else [0])
    
    #plt.violinplot(G4perATbin,positions=ticks)
    plt.violinplot(G4perATbin)
    plt.title('G4s per AT content in {} bp windows'.format(windowSize))
    plt.xlabel('AT content (%)')
    plt.ylabel('# G4 (min GG.G tracts)')
    plt.xticks(range(1,len(ticks)+1),ticks)
    for i in range(len(G4perATbin)):
        g4, counts = np.unique(G4perATbin[i],return_counts=True)
        plt.scatter([i+1]*g4.shape[0],g4,s=(1+30*np.log10(counts)),c=counts,alpha=0.7)
    
    plt.show(block=False)    
