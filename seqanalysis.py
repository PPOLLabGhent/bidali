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
