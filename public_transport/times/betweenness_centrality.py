import numpy as np
import pandas as pd
import networkx as nx

filename = "transport_n_IJ.txt"

edgelist = pd.read_csv(filename, sep=' ', header=None, names=['i', 'j', 'value'])
edgelist['distance'] = 1/edgelist['value']
graph = nx.from_pandas_edgelist(edgelist, 'i', 'j', ['value', 'distance'], create_using=nx.DiGraph)

centrality = nx.betweenness_centrality(graph, weight='distance')


centrality_df = pd.DataFrame.from_dict(centrality, orient='index', columns=['centrality']).sort_index()
centrality_df.to_csv('betweenness_centrality.txt', sep=' ', header=False, index=True)