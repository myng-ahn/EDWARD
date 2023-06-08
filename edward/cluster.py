import numpy as np
from igraph import Graph
from sknetwork.clustering import Louvain
from sklearn.neighbors import NearestNeighbors

def knn(data, k):
    """ Perform k nearest neighbors clustering
    """
    k = int(k)
    neighbors = NearestNeighbors(n_neighbors=k, algorithm='auto').fit(data)
    neighbors = neighbors.kneighbors_graph(data).toarray()
    return neighbors

def leiden(data, k, obj="modularity", res=0.5, n_iter=-1):
    """ TODO: documentation
    """
    neighbors = knn(data, k)
    nbrs_graph = Graph.Adjacency(neighbors.tolist())
    nbrs_graph = nbrs_graph.as_undirected(mode="collapse")
    partition = nbrs_graph.community_leiden(objective_function=obj, 
                                            resolution=res, 
                                            n_iterations=n_iter)
    return partition.membership

def louvain(data, k, res=0.5, n_iter=-1):
    """ TODO: documentation
    """
    neighbors = knn(data, k)
    louvain = Louvain(resolution=res, n_aggregations=n_iter)
    labels = louvain.fit_predict(neighbors)
    return labels