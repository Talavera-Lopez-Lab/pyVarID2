import pandas as pd
import r_functions

def downsample(dataframe, n, dsn):
    pass

class SCseq():

    def __init__(
            self,
            expdata
    ):
        super().__init__()
        self.expdata = expdata
        '''The raw expression data matrix with cells as columns and genes as rows in sparse matrix format.'''
        self.ndata = expdata
        '''Filtered data with expression normalized to one for each cell.'''
        self.fdata = expdata
        self.filterpar = {}
        '''dict containing the parameters used for cell and gene filterung'''

    def filterdata(
            self,
            mintotal=3000,
            minexpr=5,
            minnumber=1,
            maxexpr=float('inf'),
            downsample=False,
            sfn=False,
            hkn=False,
            dsn=1,
            rseed=17000,
            CGenes=None,
            FGenes=None,
            ccor=0.4,
    ):

        self.filterpar = {
            'mintotal': mintotal,
            'minexpr': minexpr,
            'minnumber': minnumber,
            'maxexpr': maxexpr,
            'downsample': downsample,
            'sfn': sfn,
            'dsn': dsn,
            'CGenes': CGenes,
            'FGenes': FGenes
        } 
        self.ndata = self.expdata.loc[:, self.expdata.apply(lambda col: col.sum(skipna=True), axis=0) >= mintotal]

        if downsample:
            pass

        if sfn:
            d = r_functions.scran_computeSumFactors(self.ndata)
            self.ndata = ((self.ndata.T/d).T)+0.1

        if hkn:
            pass

        if not (downsample or sfn or hkn):
            pass

        x = self.ndata
        self.fdata = x[x.apply(lambda row: row.sum(skipna=True), axis=1) >= minnumber]
        x = self.fdata
        self.fdata = x[x.apply(lambda row: row.max(skipna=True), axis=1) < maxexpr]

        if downsample:
            pass

        if FGenes != None:
            pass

        if CGenes != None:
            pass

