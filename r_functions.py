import numpy as np
from scipy.spatial.distance import pdist, squareform

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