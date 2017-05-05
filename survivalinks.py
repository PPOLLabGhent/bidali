#!/usr/bin/env python
# Functions for linking survival to genes

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt

#Calculate survival impact for each gene
def geneImpactSurvival(gene,expressions,metadata,groupingByQuantile=0.5,grouping=None,filter=None,
                       metacensorcol="overall_survival",metaDFDcol="death_from_disease",
                       plot=False,rounding=2):
    if filter is not None:
        expressions = expressions[expressions.columns[filter]]
        metadata = metadata[filter]
    if grouping is None:
        cutoff = expressions.ix[gene].quantile(groupingByQuantile) if groupingByQuantile else expressions.ix[gene].mean()
        groupHigh = expressions.ix[gene] > cutoff
    else: groupHigh = grouping

    kmf = KaplanMeierFitter()
    
    kmf.fit(metadata[metacensorcol][~groupHigh], metadata[metaDFDcol][~groupHigh], label='low')
    lastlow = float(kmf.survival_function_.ix[kmf.survival_function_.last_valid_index()])
    if plot: ax = kmf.plot()

    kmf.fit(metadata[metacensorcol][groupHigh], metadata[metaDFDcol][groupHigh], label='high')
    lasthigh = float(kmf.survival_function_.ix[kmf.survival_function_.last_valid_index()])
    if plot: kmf.plot(ax=ax)

    results = logrank_test(metadata[metacensorcol][groupHigh], metadata[metacensorcol][~groupHigh],
                           metadata[metaDFDcol][groupHigh], metadata[metaDFDcol][~groupHigh], alpha=.99)
    #results.print_summary()
    if not rounding:
        return (lastlow-lasthigh,results.p_value)
    else:
        return (round(lastlow-lasthigh,rounding),round(results.p_value,rounding))

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
