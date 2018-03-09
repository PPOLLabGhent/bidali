#!/usr/bin/env python
# Plotting imports
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
import seaborn as sns
# General configuration
sns.set(style="whitegrid", palette="pastel", color_codes=True)
# Set matplotlib to interactive mode when executed from interactive shell
import os
if 'ps1' in vars(os.sys):
    plt.ion()
    fig, ax = plt.subplots() #TODO for the moments prevents crashes, but very hacky solution

# Data processing imports
import numpy as np, pandas as pd
import networkx as nx
from inspect import getmembers
from unittest.mock import MagicMock as Mock
from collections import OrderedDict
from bidali import LSD
import itertools

#TODO http://bokeh.pydata.org/en/latest/docs/gallery.html#gallery

# PLotting functions
def plotGeneCounts(data,x,hue,y='counts',dodge=True,jitter=False,ax=None,**kwargs):
        """
        Plot counts for a gene.
        The default counts columns x is set as 'counts'
        
        e.g. x = 'cellline'
             hue = 'treatment'
        """
        if ax is None: fig, ax = plt.subplots()
        ax = sns.stripplot(
            x=x, y=y, hue=hue, data=data,
            dodge=dodge, jitter=jitter, ax=ax, **kwargs
        )
        return ax

def drawGeneEnvNetwork(gene,interactome='string',addNeighborEdges=True,node_color='r',layout='random_layout'):
    """
    interactome -> one of the networks in LSD.get_proteinNetworks
    >> drawGeneEnvNetwork('BRIP1') # doctest: +ELLIPSIS
    <matplotlib.figure.Figure ...>
    """
    nws = LSD.get_proteinNetworks()
    nw = nws.__getattribute__(interactome+'nx')
    G = nx.Graph()
    neighbors = nw.neighbors(gene)
    [G.add_edge('BRIP1',n) for n in neighbors]
    if addNeighborEdges:
        for n in neighbors:
            nn = nw.neighbors(n)
            for nnn in nn:
                if nnn in neighbors: G.add_edge(n,nnn)
    fig,ax = plt.subplots()
    nx.draw_networkx(G,pos=drawGeneEnvNetwork.layouts[layout](G),
                     with_labels=True,
                     node_size=40,node_color=node_color,
                     ax=ax)
    ax.axis('off')
    return fig
drawGeneEnvNetwork.layouts = dict(getmembers(nx.layout))

def drawCNAcircos(cnaPositions,cnaTotal=False,chrRange=None,sortPositions=True,
                  genePositions=None,geneAnnotations=False,
                  startAngle=0,color='r',wedgebgshade='0.9',genecolor='k',ax=None):
    """
    color: either one color or a list of cnaPositions size

    Example BRIP1 on 17q
    >>> drawCNAcircos([(36094885,83084062),(59577514,83084062)],cnaTotal=10,chrRange=(26885980,83257441),
    ... genePositions={'BRIP1':61863521}) # doctest: +ELLIPSIS
    <matplotlib.figure.Figure ...>
    """
    if sortPositions:
        cnaPositions = sorted(cnaPositions,key=lambda x: max(x)-min(x),reverse=True)
    if not cnaTotal: cnaTotal = len(cnaPositions)
    wedgeDegrees = 360/cnaTotal
    maxChr,minChr = max(chrRange),min(chrRange)
    r = maxChr-minChr

    if ax: fig = ax.get_figure()
    else: fig,ax = plt.subplots()
    
    for cna in cnaPositions:
        if wedgebgshade: ax.add_patch(ptch.Wedge((0,0),r,startAngle,startAngle+wedgeDegrees,fc=wedgebgshade))
        wedge = ptch.Wedge((0,0),max(cna)-minChr,startAngle,startAngle+wedgeDegrees,width=max(cna)-min(cna),fc=color)#,ec=color)
        startAngle+=wedgeDegrees
        ax.add_patch(wedge)
    # Add circle edge for chr zoom
    ax.add_patch(ptch.Circle((0,0), radius=r, ec='k', fc='none'))
    # Add genes
    for g in genePositions:
        ax.add_patch(ptch.Circle((0,0), radius=genePositions[g]-minChr, ec=genecolor, fc='none'))
        if geneAnnotations == True: raise NotImplementedError

    ax.axis('off')
    ax.set_xlim((-r,r))
    ax.set_ylim((-r,r))
    return fig

def curvedHeatPlot(dataframe,columns,topDisplayed=10,cellwidth=.2,cellheight=.1,cmap='hot_r',
                   headingTextSize=14,curveLabels=True,filename=None):
    """
    >>> curvedHeatPlot(Mock(),['col1','col2']) # doctest: +ELLIPSIS
    <matplotlib.figure.Figure ...>
    """
    from itertools import count
    cmap = plt.get_cmap(cmap)
    topDisplayed_l_pos,topDisplayed_r_pos = topDisplayed,len(dataframe)-topDisplayed
    def curvedHeat(x,iterposition,ax,columns=columns):
        position = next(iterposition)
        if position < topDisplayed_l_pos:
            for c in columns:
                ax.add_patch(ptch.Rectangle((-1+cellwidth*columns.index(c),(topDisplayed_l_pos - position - 1)*cellheight),
                                            width=cellwidth,height=cellheight,color=cmap(x[c])))
                if position == 0:
                    ax.annotate('{:.2f}'.format(x[c]),
                                (-1+cellwidth*columns.index(c)+cellwidth/2,((topDisplayed_l_pos - position - 1)*cellheight)+.05),
                                ha='center',va='center')
            ax.annotate(x.name,(-1,((topDisplayed_l_pos - position - 1)*cellheight)+.05),ha='right',va='center')
        elif position >= topDisplayed_r_pos:
            for c in columns:
                ax.add_patch(ptch.Rectangle((1-cellwidth-(cellwidth*columns.index(c)),(position-topDisplayed_r_pos)*cellheight),
                                            width=cellwidth,height=cellheight,color=cmap(x[c])))
                if position == 57:
                    ax.annotate('{:.2f}'.format(x[c]),
                                (1-cellwidth*columns.index(c)-cellwidth/2,((position-topDisplayed_r_pos)*cellheight)+.05),
                                ha='center',va='center')
            ax.annotate(x.name,(1,((position-topDisplayed_r_pos)*cellheight)+.05),ha='left',va='center')
        else:
            startAngle = 180+(180*(position - topDisplayed)/(topDisplayed_r_pos - topDisplayed))
            endAngle = 180+(180*(1+position - topDisplayed)/(topDisplayed_r_pos - topDisplayed))
            for c in columns:
                ax.add_patch(ptch.Wedge((0,0),1-(cellwidth*columns.index(c)),startAngle,endAngle,cellwidth,color=cmap(x[c])))
            if curveLabels: ax.annotate(x.name,(np.pi*(startAngle+endAngle)/360,1),xycoords='polar',
                                        ha='right' if position < len(dataframe)/2 else 'left',va='top',
                                        rotation=(startAngle+endAngle)/2 + (180 if position < len(dataframe)/2 else 0))
    figheatmap,ax = plt.subplots(figsize=(6,4))
    ax.set_xlim((-1.2,1.2))
    ax.set_ylim((-1.2,(topDisplayed+2)*cellheight))
    for c in columns:
        ax.annotate(c,(-1+cellwidth*columns.index(c)+cellwidth/2,topDisplayed*cellheight+.05),ha='center',size=headingTextSize)
        ax.annotate(c,(1-cellwidth*columns.index(c)-cellwidth/2,topDisplayed*cellheight+.05),ha='center',size=headingTextSize)

    ax.axis('off')
    c = count(0)
    dataframe.T.apply(curvedHeat,args=(c,ax))
    
    if filename: figheatmap.savefig(filename,transparent=True)
    return figheatmap

def dosageViolin(gene,dataset,ax=None,cntype='gain',risksToPlot=3):
    """
    >>> dosageViolin('BRIP1',dataset=Mock()) # doctest: +SKIP
    """
    ds = dataset
    genedosage = pd.DataFrame({'expression': ds.exprdata.ix[gene],
                               'cna': ds.geneCNA.ix[gene],
    }).dropna()
    genedosage['risk_status'] = (genedosage.T.apply(lambda x: 'high' if ds.metadata.ix[x.name].high_risk=='1' else 'low')
                                 if risksToPlot == 2 else
                                 genedosage.T.apply(lambda x: 'lowrisk' if ds.metadata.ix[x.name].high_risk == '0' else
                                                   'highrisk_amp' if ds.metadata.ix[x.name].mycn_status == '1' else 'highrisk_sc')
    )
    
    genedosage = genedosage[genedosage.cna!=('loss' if cntype=='gain' else 'gain')]
    groupsizes = genedosage.groupby(['cna','risk_status']).size()
    #print(groupsizes)
        
    #Plots
    viofig,violax = (ax.get_figure(),ax) if ax else plt.subplots()
    try: violax.set_title('{} ({})'.format(gene,ds.name))
    except AttributeError: violax.set_title(gene)
    
    if risksToPlot == 2:
        order=['normal',cntype]
        hue_order = ['low','high']
        combinations = list(itertools.product(order,hue_order))
        sns.violinplot(x='cna', y='expression', hue='risk_status', data=genedosage, split=True,
                       inner='points', palette={'high':'b','low':'y'},order=order, hue_order=hue_order, ax=violax)
    elif risksToPlot == 3:
        order = ['lowrisk','highrisk_sc','highrisk_amp']
        hue_order = ['normal',cntype]
        combinations = list(itertools.product(order,hue_order))
        sns.violinplot(x='risk_status', y='expression', hue='cna', data=genedosage, split=True, inner='points',
                       palette={'normal':'y',cntype:'b' if cntype=='loss' else 'r'}, order=order, hue_order=hue_order, ax=violax)
    sns.despine(left=True)
    
    #Annotate sizes
    for c,i in zip(combinations,range(len(combinations))):
        violax.annotate(str(groupsizes[c[::-1]]) if c[::-1] in groupsizes else '0',(int(i/2)+(0.2 if i%2 else -0.2),.5),
                        xycoords=('data','axes fraction'),ha='center')
    
    return viofig

def draw_cellcycle(annotations,totalCompartments=100,phases=OrderedDict([('G1',0),('S',45),('G2',65),('M',85)]),ax=None,innerAnnotRing=.5):
    """
    annotations => pd.Series with genes in index and peak time values

    >>> ax = draw_cellcycle()
    """
    xytransform = lambda r,t: (r*np.cos(2*np.pi*(360-t)/360),r*np.sin(2*np.pi*(360-t)/360))
    if ax:
        fig = ax.get_figure()
    else: fig,ax = plt.subplots()
    padding = .1
    innerRadius = .8
    ax.axis('off')
    ax.set_xlim((-1-padding,1+padding))
    ax.set_ylim((-1-padding,1+padding))
    ax.add_patch(ptch.Circle(xy=(0,0),radius=1,fill=False,lw=2))
    ax.add_patch(ptch.Circle(xy=(0,0),radius=innerRadius,fill=False,lw=2))
    compartments = np.arange(0,360,360/totalCompartments) - 90
    phaseStart = list(phases.values())
    for i,p in enumerate(phases):
        ax.add_patch(ptch.PathPatch(ptch.Path([xytransform(1,compartments[phases[p]]),
                                               xytransform(innerRadius,compartments[phases[p]])]),color='k',lw=5))
        anAngle = np.mean((compartments[phaseStart[i]],compartments[phaseStart[i+1]] if i+1 < len(phases) else 360-90))
        ax.annotate(p,xytransform(np.mean((1,innerRadius)),anAngle),va='center',ha='center')
    for p,grp in annotations.groupby(by=annotations):
        ax.add_patch(ptch.PathPatch(ptch.Path([xytransform(1,compartments[p]),
                                               xytransform(innerRadius,compartments[p])]),color='r',lw=3))
        ax.annotate(','.join(grp.index),xytransform((1+padding) if len(grp) == 1 else innerAnnotRing,compartments[p]),va='center',ha='center',rotation=180-compartments[p])
    return ax
