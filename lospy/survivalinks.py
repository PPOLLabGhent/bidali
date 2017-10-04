#!/usr/bin/env python
# Functions for linking survival to genes

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt

# Calculate survival impact for each gene
def geneImpactSurvival(gene,expressions,metadata,groupingByQuantile=0.5,grouping=None,filter=None,
                       metacensorcol="overall_survival",metaDFDcol="death_from_disease",
                       plot=False,rounding=2):
    """
    Returns None if no meaningful groups
    """
    if filter is not None:
        expressions = expressions[expressions.columns[filter]]
        metadata = metadata[filter]
    if grouping is None:
        cutoff = expressions.ix[gene].quantile(groupingByQuantile) if groupingByQuantile else expressions.ix[gene].mean()
        groupHigh = expressions.ix[gene] > cutoff
    else: groupHigh = grouping

    kmf = KaplanMeierFitter()
    
    try: kmf.fit(metadata[metacensorcol][~groupHigh], metadata[metaDFDcol][~groupHigh], label='low')
    except ValueError: return None
    lastlow = float(kmf.survival_function_.ix[kmf.survival_function_.last_valid_index()])
    if plot:
        if type(plot) is bool:
            ax = kmf.plot()
        else: ax = kmf.plot(ax=plot)

    try: kmf.fit(metadata[metacensorcol][groupHigh], metadata[metaDFDcol][groupHigh], label='high')
    except ValueError: return None
    lasthigh = float(kmf.survival_function_.ix[kmf.survival_function_.last_valid_index()])
    if plot: kmf.plot(ax=ax)

    results = logrank_test(metadata[metacensorcol][groupHigh], metadata[metacensorcol][~groupHigh],
                           metadata[metaDFDcol][groupHigh], metadata[metaDFDcol][~groupHigh], alpha=.99)
    #results.print_summary()
    if not rounding:
        result = (lastlow-lasthigh,results.p_value)
    else:
        result = (round(lastlow-lasthigh,rounding),round(results.p_value,rounding))
    if plot: return result,ax
    else: return result

def geneCombinationImpactSurvival(genes,expressions,metadata,groupingByQuantile=0.5,
                       metacensorcol="overall_survival",metaDFDcol="death_from_disease",
                       plot=False,rounding=2):
    from itertools import combinations
    cutoff = {gene:expressions.ix[gene].quantile(groupingByQuantile) if groupingByQuantile else expressions.ix[gene].mean() for gene in genes}
    groupHigh = [expressions.ix[gene] > cutoff[gene] for gene in genes]
    kmf = KaplanMeierFitter()

    genis = [(i,s) for i in range(len(genes)) for s in ('H','L')]
    lastvalues = {}
    
    for combi in set(combinations(['H','L']*len(genes),len(genes))):
        selection = groupHigh[0] if combi[0] == 'H' else ~groupHigh[0]
        for gsel in zip(range(1,len(genes)),combi[1:]):
            selection = selection & (groupHigh[gsel[0]] if gsel[1] == 'H' else ~groupHigh[gsel[0]])
        if sum(selection) == 0: continue
        kmf.fit(metadata[metacensorcol][selection], metadata[metaDFDcol][selection], label=''.join(combi))
        lastvalues[combi] = (sum(selection),float(kmf.survival_function_.ix[kmf.survival_function_.last_valid_index()]))
        try: kmf.plot(ax=ax)
        except NameError: ax = kmf.plot()

    ax.set_title(':'.join(genes))
    return lastvalues

    #results = logrank_test(metadata[metacensorcol][groupHigh], metadata[metacensorcol][~groupHigh],
    #                       metadata[metaDFDcol][groupHigh], metadata[metaDFDcol][~groupHigh], alpha=.99)
    #results.print_summary()
    #if not rounding:
    #    return (lastlow-lasthigh,results.p_value)
    #else:
    #    return (round(lastlow-lasthigh,rounding),round(results.p_value,rounding))

def subsetsImpactSurvival(subsets,metadata,metacensorcol="overall_survival",
                          metaDFDcol="death_from_disease",plot=False,title=None,rounding=2):
    """
    subsets is a dictionary,
    e.g.: subsets={'cluster {}'.format(i):metadata.index.isin(fitrue.columns[kmeans.labels_==i]) for i in range(4)}
    """
    kmf = KaplanMeierFitter()

    lastvalues = {}
    for subset in subsets:
        kmf.fit(metadata[metacensorcol][subsets[subset]], metadata[metaDFDcol][subsets[subset]], label=subset)
        lastvalues[subset] = (sum(subsets[subset]),float(kmf.survival_function_.ix[kmf.survival_function_.last_valid_index()]))
        try: kmf.plot(ax=ax)
        except NameError: ax = kmf.plot()

    if title: ax.set_title(title)
    return lastvalues

# 2-way survival curves
def twoGeneSurvivalPlot(g1,g2,expressionData,metadata,ax=None):
    if not ax: f,ax = plt.subplots()
    g1med = expressionData.ix[g1].median()
    g2med = expressionData.ix[g2].median()
    gHH = expressionData.columns[(expressionData.ix[g1]>g1med)&(expressionData.ix[g2]>g2med)] #High g1, high g2
    gHL = expressionData.columns[(expressionData.ix[g1]>g1med)&(expressionData.ix[g2]<=g2med)]
    gLH = expressionData.columns[(expressionData.ix[g1]<=g1med)&(expressionData.ix[g2]>g2med)]
    gLL = expressionData.columns[(expressionData.ix[g1]<=g1med)&(expressionData.ix[g2]<=g2med)]
    kmf = KaplanMeierFitter()
    for g,a in zip((gHH,gHL,gLH,gLL),('HH','HL','LH','LL')):
        kmf.fit(metadata.lastfollowup[metadata.index.isin(g)],
                ~metadata.survived[metadata.index.isin(g)], label='{}[{}]'.format(a,len(g)))
        kmf.plot(ax=ax)
    ax.set_title('{} and {}'.format(g1,g2))    

def twoGeneSurvivalPlots(expressionData,metadata,genes=('BRIP1','TK1','BIRC5','EME1','WNT3')):
    from itertools import combinations, count
    numberOfPlots = len(list(combinations(range(len(genes)),2)))
    prows = int(np.ceil(numberOfPlots**0.5))
    pcols = int(np.floor(numberOfPlots**0.5))
    f,axes = plt.subplots(prows,pcols,sharex=True)
    c = count(0)
    for g1 in genes:
        for g2 in genes[genes.index(g1)+1:]:
            ax = axes.flatten()[next(c)]
            twoGeneSurvivalPlot(g1,g2,expressionData=expressionData,metadata=metadata,ax=ax)
