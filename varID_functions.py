import pandas as pd
import numpy as np
import multiprocessing
import r_functions as r
import statsmodels.api as sm

def getFormula(b=None, reg_var=None):
    if reg_var is None:
        model_formula = 'beta'
    else:
        model_formula = 'beta + ' + ' + '.join(reg_var.columns)
    
    if b is not None:
        model_formula = f'( {model_formula} ): b + b + 0'
    
    model_formula = f'x ~ {model_formula}'
    
    return model_formula


def getRegData(k, b=None, reg_var=None):
    # Create a DataFrame with 'beta' column set to the values of k
    reg_data = pd.DataFrame({'beta': k})

    # If reg_var is not None, add its columns to reg_data based on the index of k
    if reg_var is not None:
        reg_var_subset = reg_var.loc[k.index]  # Subset reg_var to match the index of k
        reg_data = reg_data.join(reg_var_subset)
    
    # If b is not None, add the 'b' column to reg_data based on the index of k
    if b is not None:
        reg_data['b'] = b.loc[k.index]
    
    return reg_data

def compResiduals(
        expData,
        batch=None,
        regVar=None,
        span=0.75,
        no_cores=1,
        n_genes=None,
        seed=12345
):
    k = np.log(expData.sum(axis=0))
    modelFormula = getFormula(batch, regVar)
    regData = getRegData(k, batch, regVar)
    # Generate random Poisson-distributed values
    x = np.random.poisson(5, regData.shape[0])
    # Add the new 'x' column to regData
    regData['x'] = x
    # Define the model formula
    model_formula = 'x ~ beta + var1 + var2'
    # Fit the GLM model with a Poisson family
    model = sm.GLM.from_formula(model_formula, data=regData, family=sm.families.Poisson()).fit()
    # Get the names of the coefficients
    reg_names = model.params.index.tolist()
    # Compute the row means of expData
    mx = expData.mean(axis=1)
    # Get the row names (index) of expData
    genes = expData.index.tolist()

def compResiduals0():
    pass

def pruneKnn(
        expData,
        distM=None,
        large=True,
        regNB=True,
        bmethod=None,
        batch=None,
        regVar=None,
        offsetModel=True,
        thetaML=False,
        theta=10,
        ngenes=2000,
        span=0.75,
        pcaComp=None,
        tol=1e-5,
        algorithm='kd_tree',
        metric='pearson',
        genes=None,
        knn=25,
        do_prune=True,
        alpha=1,
        nb=3,
        no_cores=None,
        FSelect=False,
        pca_scale=False,
        ps=1,
        seed=12345

):
    if type(expData) == pd.core.frame.DataFrame:
        if ps < 0:
            raise ValueError('Pseudocount needs to be greater or equal to 0')
        rs = np.sum(expData > 0, axis=1)
        cs = np.sum(expData > 0, axis=0)
        expData.loc[(rs > 0)][cs[cs > 0].index]
        genes = expData.index if genes is None else genes
        batch = batch[expData.columns] if batch != None else batch
    if batch == None and bmethod == None:
        if bmethod == 'harmony':
            hflag = True
            hbatch = batch
            batch = None
    expData = expData.loc[genes]
    colS = cs[cs > 0]
    Xpca = None
    bg = None
    no_cores = max(no_cores, multiprocessing.cpu_count()) if no_cores == None else min(no_cores, multiprocessing.cpu_count())

    if large:
        distM = None
        if regNB:
            if offsetModel:
                #regdata = compresiduals0()