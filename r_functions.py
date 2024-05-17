import numpy as np
import pandas as pd

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

def cor(
        x,
        y=None,
        use='everything',
        method='pearson'
):
    '''not entirely sure this works as intended yet'''
    possible_use_values = [
        'all.obs',
        'complete.obs',
        'pairwise.complete.obs',
        'everything',
        'na.or.complete',
    ]
    possible_method_values = ['pearson', 'kendall', 'spearman']
    if use not in possible_use_values:
        raise ValueError(f'invalid "use" argument, must be one of {possible_use_values}')
    if method not in possible_method_values:
        raise ValueError(f'invalid "method" argument, must be one of {possible_method_values}')
    
    def rank(u):
        '''implementation only needed if methods == "kendall" or "spearman"'''
        pass

    if method == 'pearson':
        corr = np.corrcoef(x.T, y.T, rowvar=False)[:x.T.shape[1], x.T.shape[1]:]
        corr = pd.DataFrame(data=corr, index=x.T.columns, columns=y.T.columns)
        #corr = np.corrcoef(x, rowvar=False)
    elif method in ['kendall', 'spearman']:
        print("Methods kendall and spearman are not yet implemented (and likely won't be)")
    
    return corr