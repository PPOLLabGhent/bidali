#!/usr/bin/env python
import matplotlib.pyplot as plt
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
