import numpy as np
from igraph import Graph
from sknetwork.clustering import Louvain
from sklearn.neighbors import NearestNeighbors

def knn(data, k):
    """ Perform k nearest neighbors clustering using sklearn implementation of 
    knn clustering.

    Parameters
    ----------
    data : ndarray
        The data on which to perform knn clustering. Rows should represent samples
        or variants and columns should represent genes or variants.
    k : int

    Returns
    -------
    ndarray
        A sparse ndarray where the i, j element represents the weight between
        node i and node j.
    """
    k = int(k)
    # fit and transform data using sklearn clustering with given k
    neighbors = NearestNeighbors(n_neighbors=k, algorithm='auto').fit(data)
    neighbors = neighbors.kneighbors_graph(data).toarray()
    return neighbors

def leiden(data, k, obj="modularity", res=0.5, n_iter=-1):
    """ Perform leiden clustering using igraph implementation. 

    Parameters
    ----------
    data : ndarray
        The data on which to perform leiden clustering
    k : int
    obj: str
        "CPM" or "modularity
    res: float
        Resolution parameters. Smaller resolution generates smaller clusters.
    n_iter: int
        Number of iterations to run leiden clustering. A negative number will run 
        the algorithm until the objective function no longer improves. 

    Returns
    -------
    array:
        A nx1 list where for each of the n samples there is a number corresponding
        to which cluster it belongs to.
    """
    neighbors = knn(data, k) # perform knn with given value of k
    nbrs_graph = Graph.Adjacency(neighbors.tolist())
    nbrs_graph = nbrs_graph.as_undirected(mode="collapse")
    partition = nbrs_graph.community_leiden(objective_function=obj, 
                                            resolution=res, 
                                            n_iterations=n_iter)
    return partition.membership

def louvain(data, k, res=0.5, n_iter=-1):
    """ Perform louvain clustering using igraph implementation. 

    Parameters
    ----------
    data : ndarray
        The data on which to perform leiden clustering
    k : int
    res: float
        Resolution parameters. Smaller resolution generates smaller clusters.
    n_iter: int
        Number of iterations to run leiden clustering. A negative number will run 
        the algorithm until the objective function no longer improves. 

    Returns
    -------
    array:
        A nx1 list where for each of the n samples there is a number corresponding
        to which cluster it belongs to.
    """
    neighbors = knn(data, k) # perform knn with given value of k
    louvain = Louvain(resolution=res, n_aggregations=n_iter)
    labels = louvain.fit_predict(neighbors)
    return labels