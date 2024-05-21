import pandas as pd
import numpy as np
import scipy
from scipy.sparse import csr_matrix
import r_functions as r
import RaceID_utils

def downsample(dataframe, n, dsn):
    pass

class SCseq():

    def __init__(
            self,
            expdata
    ):
        super().__init__()        

        if expdata.shape[0] < 2:
            raise ValueError('expdata must have more than one row')
        if expdata.shape[1] < 2:
            raise ValueError('expdata must have more than one column')
        if expdata.isna().any().any():
            raise ValueError("NA values are not allowed in expdata")
        if (expdata < 0).any().any():
            raise ValueError("Negative values are not allowed in expdata")

        self.expdata = expdata
        #self.expdata.indices = expdata.columns
        #self.expdata.indptr = expdata.index
        '''The raw expression data matrix with cells as columns and genes as rows in sparse matrix format.'''
        self.cluster = {}
        '''dict storing information on the initial clustering step of the RaceID3 algorithm'''
        self.background = {}
        '''dict storing the polynomial fit for the background model of gene expressoin variability computed by RaceID3, 
        which is used for outlier identification'''
        self.filterpar = {}
        '''dict containing the parameters used for cell and gene filterung'''



    def filterdata(
            self,
            mintotal=3000,
            minexpr=5,
            minnumber=5,
            LBatch=None,
            knn=10,
            CGenes=[],
            FGenes=[],
            ccor=0.4,
            bmode='RaceID',
            verbose=True,
    ):
        '''
        This function allows filtering of genes and cells to be used in the RaceID3 analysis. 
        not yet -> It also can perform batch effect correction using an internal method or a recently published alternative \code{mnnCorrect} from the \pkg{batchelor} package.

        Args:
            mintotal (int): minimum total transcript number required. Cells with less than \code{mintotal} transcripts are filtered out. Default is 3000.
        '''
        if not isinstance(mintotal, (int, float)):
            raise ValueError('mintotal has to be a numeric value')
        if mintotal <= 0:
            raise ValueError('mintotal has to be a positive number')
        if not isinstance(minexpr, (int, float)):
            raise ValueError('mintotal has to be a numeric value')
        if minexpr < 0:
            raise ValueError('mintotal has to be a non-negative number')
        if not isinstance(minnumber, int):
            raise ValueError('minnumber has to be a non-negative integer value')
        if minnumber < 0:
            raise ValueError('minnumber has to be a non-negative integer value')
        if not bmode in ['RaceID', 'mnnCorrect']:
            raise ValueError('bmode has to be one of RaceID, mnnCorrect')
        
        self.dimRed = []

        #total transcript counts
        counts = self.expdata.sum()
        #filtering cells
        f = counts >= mintotal
        self.counts = counts[f]

        # Filtering of genes and normalization
        g = (self.expdata.loc[:, f].apply(lambda x: x >= minexpr, axis=1).sum(axis=1) >= minnumber)
        self.ndata = self.expdata.loc[:,f].div(counts[f])
        genes = self.ndata.index[g]
        genes = genes[~genes.isin(FGenes)]

        ## batch effect correctoin by discarding batch effect signature genes
        #bG = None
        bl = None
        #if LBatch is not None and len(LBatch) > 1 and bmode == 'RaceID':
        #    x = self.ndata.loc[genes, :].values * min(self.counts) + 0.1
        #    d = r.calculate_spearman_distance(x)

        # discard genes correlating to genes in CGenes
        # this section isnt fully tested to be working correctly yet
        if CGenes:
            CGenes = [gene for gene in CGenes if gene in genes]
            h = None
            if len(CGenes) > 0:
                if len(CGenes) == 1:
                    pass
                else:
                    genes_df = self.ndata[self.ndata.index.isin(genes)]
                    CGenes_df = self.ndata[self.ndata.index.isin(CGenes)]
                    k = r.cor(genes_df, CGenes_df)
                h = np.apply_along_axis(lambda x: np.max(np.abs(x)), 1, k) < ccor
                h[np.isnan(h)] = True
            if len(h) > 0:
                genes = [gene for i, gene in enumerate(genes) if h[i]]


        self.genes = genes
        self.filterpar = {
            'mintotal': mintotal,
            'minexpr': minexpr,
            'minnumber': minnumber,
            'CGenes': CGenes,
            'FGenes': FGenes,
            'BGenes': bl,
            'bmode': bmode,
        } 

        # compute polynomial fit of background model used for outlier identification from non-normalized data
        bg = RaceID_utils.fitbackground(self.expdata.loc[self.genes, :][self.ndata.columns]) 
        self.background['vfit'] = bg['fit']

        # compute genes with variability above background level from normalized data for feature selection
        bg = RaceID_utils.fitbackground(self.getfdata())
        self.cluster['features'] = bg['n']

        # batch correction by batchelor::mnnCorrect after filtering (skipping this for now)
        if LBatch is not None and len(LBatch) > 1 and bmode == 'mnnCorrect':
            pass

        #self.ndata = self.expdata.loc[:, self.expdata.apply(lambda col: col.sum(skipna=True), axis=0) >= mintotal]

        #x = self.ndata
        #self.fdata = x[x.apply(lambda row: row.sum(skipna=True), axis=1) >= minnumber]
        #x = self.fdata
        #self.fdata = x[x.apply(lambda row: row.max(skipna=True), axis=1) < maxexpr]

        #if downsample:
        #    pass

        #if FGenes != None:
        #    pass

        #if CGenes != None:
        #    pass
    def getfdata(self, g=None, n=None):
        fgenes = self.genes if g is None else [gene for gene in self.ndata.index if gene in g]
        n = list(self.counts.index) if n is None else [name for name in self.counts.index if name in n]
        return self.ndata.loc[fgenes, n] * self.counts.min() + 0.1

    def getExpData(self, genes=None):
        expdata = self.expdata
        if genes == None:
            genes = self.genes
        return expdata.loc[genes, :][self.ndata.columns]



