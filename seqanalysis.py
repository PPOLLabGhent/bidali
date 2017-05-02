#!/bin/env python3
#CVN set of functions for sequence analysis
#General imports
import numpy as np

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

def loadHumanGenome():
    from glob import glob
    files = glob('/home/christophe/Data/Genomes/chroms/chr??.fa')
    files += glob('/home/christophe/Data/Genomes/chroms/chr?.fa')
    chromosomes = []
    for f in files:
        f = open(f).readlines()
        chromosomes.append(DNA(f.pop(0).strip()[1:],''.join([l.strip() for l in f])))
    return Genome('human','GRCh38',chromosomes)

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
