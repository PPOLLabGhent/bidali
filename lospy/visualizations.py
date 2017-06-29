#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
import numpy as np
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

def curvedHeatPlot(dataframe,columns,topDisplayed=10,cellwidth=.2,cellheight=.1,cmap='hot_r',
                   headingTextSize=14,curveLabels=True,filename=None):
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
