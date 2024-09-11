import numpy as np
import pandas as pd
import networkx as nx

filename = "transport_n_IJ.txt"

edgelist = pd.read_csv(filename, sep=' ', header=None, names=['i', 'j', 'value'])
graph = nx.from_pandas_edgelist(edgelist, 'i', 'j', ['value'], create_using=nx.DiGraph)

centrality = nx.eigenvector_centrality(graph, weight='value')

centrality_df = pd.DataFrame.from_dict(centrality, orient='index', columns=['centrality']).sort_index()
centrality_df.to_csv('eigenvector_centrality.txt', sep=' ', header=False, index=True)