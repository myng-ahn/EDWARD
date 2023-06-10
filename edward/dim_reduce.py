import numpy as np
import umap as ump
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE

def pca(data, nComponents):
    """ Perform Principal Component Analysis

    Calculates principal components (PC) of data and apply dimensionality 
    reduction using the top nComponents PCs.

    Parameters
    ----------
    data : ndarray
        The data on which to perform PCA. Rows represent observations and 
        columns represent features.
    nComponents : int > 0
        Number of PCs to use for dimensionality reduction.

    Returns
    -------
    ndarray
        PCA dimensionality reduced data.

    Raises
    ------
    ValueError
        if nComponents is greater than the number of features in original data.
    """
    
    if nComponents > data.shape[1]:
        raise ValueError("nComponents cannot exceed number of columns in data.") 

    # perform mean centering
    mean = np.mean(data, axis=0)
    norm_data = data - mean
    
    # compute eigenvector, eigenvalue of covariance matrix
    covar = np.cov(norm_data, rowvar=False)
    covar_eigvals, covar_eigvecs = np.linalg.eig(covar)
    
    # sort eigenvalues
    idxs = np.argsort(covar_eigvals)[::-1]
    sorted_eigvals = covar_eigvals[idxs]
    sorted_eigvecs = covar_eigvecs[:,idxs]

    # project onto PCA space and keep only nComponents
    pca_proj = np.dot(norm_data, sorted_eigvecs)
    pca_proj = pca_proj[:, 0:nComponents]

    return pca_proj, sorted_eigvals, sorted_eigvecs

def tsne(data, perplexity):
    """ Perform t-SNE dimensionality reduction.

    Performs automatic PCA of data followed by sklearn implementation of 
    t-SNE dimensionality reduction. 

    Parameters
    ----------
    data : ndarray
        The data on which to perform t-SNE. Rows represent observations and 
        columns represent features.
    perplexity : int
        Related to number of nearby points considered as nearest neighbors.

    Returns
    -------
    ndarray
        t-SNE reduced data.

    Raises
    ------
    ValueError
        if perplexity is greater than number of samples in data.
    """
    if perplexity > data.shape[0]:
        raise ValueError("Perplexity cannot exceed number of samples in data.")

    # run PCA with 50 PCs
    scaled_data = StandardScaler().fit_transform(data)
    num_pca_components = min(50, data.shape[1])
    pca_transformed, __, __ = pca(scaled_data, num_pca_components)

    # fit and transform data using sklearn implementation of t-SNE
    tsne = TSNE(n_components=2, perplexity=perplexity)
    # pca_transformed = np.array(pca_transformed, dtype='float32')
    tsne_transformed = tsne.fit_transform(pca_transformed)
    return tsne_transformed

def umap(data):
    """ Perform UMAP dimensionality reduction.

    Performs automatic PCA of data followed by UMAP dimensionality reduction
    using external UMAP package.

    Parameters
    ----------
    data : ndarray
        The data on which to perform UMAP reduction. Rows represent observations and 
        columns represent features.

    Returns
    -------
    ndarray
        UMAP dimensionality reduced data.
    """
    # run PCA with 50 PCs
    scaled_data = StandardScaler().fit_transform(data)
    num_pca_components = min(50, data.shape[1])
    pca_transformed, __, __ = pca(scaled_data, num_pca_components)

    # fit and transform data using external UMAP implementation with default params
    um = ump.UMAP()
    umap_transformed = um.fit_transform(pca_transformed)
    return umap_transformed