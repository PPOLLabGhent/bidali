#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
import networkx as nx
from inspect import getmembers

def drawGeneEnvNetwork(gene,interactome='string',addNeighborEdges=True,node_color='r',layout='random_layout'):
    """
    interactome -> one of the networks in LSD.get_proteinNetworks
    """
    import LSD
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
    ... genePositions={'BRIP1':61863521})
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
