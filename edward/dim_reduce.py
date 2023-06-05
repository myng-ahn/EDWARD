import numpy as np
import umap
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
        columns represents features.
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

def tsne(data, perplexity=30):
    scaled_data = StandardScaler().fit_transform(data)
    num_pca_components = min(50, data.shape[1])
    pca_transformed, __, __ = pca(scaled_data, num_pca_components)
    tsne = TSNE(n_components=2, perplexity=perplexity, learning_rate='auto')
    tsne_transformed = tsne.fit_transform(pca_transformed)
    return tsne_transformed

def umap(data):
    scaled_data = StandardScaler().fit_transform(data)
    num_pca_components = min(50, data.shape[1])
    pca_transformed, __, __ = pca(scaled_data, num_pca_components)
    umap = umap.UMAP()
    umap_transformed = umap.fit_transform(pca_transformed)
    return umap_transformed