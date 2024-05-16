import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import kendalltau, pearsonr, spearmanr

def scran_computeSumFactors(dataframe):
    pass


def calculate_spearman_distance(x):
    # Convert the input to a NumPy array
    x = np.array(x)

    # Calculate the pairwise Spearman's rank correlation distances
    # between the rows of the input matrix
    distances = pdist(x.T, metric='correlation')

    # Convert the condensed distance vector to a square distance matrix
    distance_matrix = squareform(distances)

    return distance_matrix

def cor(x, y=None, use="everything", method="pearson"):
    '''entirely chatgpt creation not tested yet'''
    # Convert data to numpy arrays
    x = np.asarray(x)
    if y is not None:
        y = np.asarray(y)
    
    # Define handling of missing values
    na_methods = ["all.obs", "complete.obs", "pairwise.complete.obs", "everything", "na.or.complete"]
    na_method_idx = na_methods.index(use) if use in na_methods else None
    if na_method_idx is None:
        raise ValueError("Invalid 'use' argument")
    
    # Select correlation method
    if method not in ["pearson", "kendall", "spearman"]:
        raise ValueError("Invalid correlation method")
    
    # Calculate correlation based on method and missing value handling
    if method == "pearson":
        if y is None:
            return pearsonr(x)
        else:
            return pearsonr(x, y)
    elif na_method_idx in [1, 4]:
        if y is None:
            return kendalltau(x)
        else:
            return kendalltau(x, y)
    else:
        if y is None:
            return spearmanr(x)
        else:
            return spearmanr(x, y)
