import pandas as pd
import numpy as np
import multiprocessing
import r_functions as r
import statsmodels.api as sm
from MASS import theta_ml
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.genmod.families import Poisson

# Assuming the theta_ml function is already defined
# def theta_ml(y, mu, n=None, weights=None, limit=10, eps=np.finfo(float).eps**0.25, trace=False):
#     # Implementation of theta_ml

def nbRegr(x, model_formula, reg_data, reg_names):
    reg_data['x'] = x
    fit_theta = False
    fit = False

    try:
        rg = GLM.from_formula(model_formula, reg_data, family=Poisson()).fit()
        fit = True
    except Exception as e:
        fit = False

    if fit:
        try:
            rg.theta = float(theta_ml(y=x, mu=rg.fittedvalues, limit=50))
            fit_theta = True
        except Exception as e:
            rg.theta = np.nan

    if fit:
        co = rg.params.to_dict()
        co['theta'] = rg.theta
        return co
    else:
        co = {name: np.nan for name in reg_names}
        co['theta'] = np.nan
        return co


def getFormula(b=None, reg_var=None):
    if reg_var is None:
        model_formula = 'beta'
    else:
        model_formula = 'beta + ' + ' + '.join(reg_var.columns)
    
    if b is not None:
        model_formula = f'( {model_formula} ): b + b + 0'
    
    model_formula = f'x ~ {model_formula}'
    
    return model_formula

def getFormula0(b=None, reg_var=None):
    if b is None and reg_var is None:
        model_formula = 'x ~ 0'
    else:
        if reg_var is not None:
            model_formula = ' + '.join(reg_var.columns)
            if b is not None:
                model_formula = f"( {model_formula} ): b + b"
            elif b is not None:
                model_formula = 'b'
            model_formula = f"x ~ {model_formula} + 0"
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
                pass #regdata = compResiduals0()
            else:
                pass #regData = compResiduals()
            z = regData.pearsonRes
        else:
            regData = None
            z = np.log((expData/colS*min(colS)) + 0.1)

        if FSelect:
            bg = fitBackVar(expData.loc[genes])
            backModel = bg.fit
            genes = bg.genes
            expData = expData.loc[genes]

        z = z.loc[genes]