import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

def fitbackground(x, mthr=-1):

    def uvar(x, fit):
        params = fit.params
        bse = fit.bse
        log2_x = np.log2(x)
        intercept = params.iloc[0]
        intercept_err = bse.iloc[0]
        ml = params.iloc[1]
        ml_err = bse.iloc[1]
        Iml2 = params.iloc[2]
        Iml2_err = bse.iloc[2]
        return 2**(intercept + intercept_err + log2_x * (ml + ml_err) + (Iml2 + Iml2_err) * log2_x ** 2)

    def calculate_vln(v, m, fit):
        # Apply the uvar function to each element in m
        predicted_variances = np.array([uvar(mean, fit) for mean in m])
        # Calculate log2 of v and the predicted variances
        vln = np.log2(v) - np.log2(predicted_variances)
        return vln

    m = np.mean(x, axis=1)
    v = np.var(x, axis=1)
    ml = np.log2(m)
    vl = np.log2(v)
    f = (ml > -np.inf) & (vl > -np.inf)
    ml = ml[f]
    vl = vl[f]
    mm = -8
    while True:
        data = pd.DataFrame({'ml': ml, 'vl': vl})
        fit = smf.ols(formula='vl ~ml + I(ml**2)', data=data).fit()
        if fit.params.iloc[2] >0 or mm >=mthr:
            break
        mm += 0.5
        f = ml > mm
        ml = ml[f]
        vl = vl[f]
    vln = calculate_vln(v, m, fit)
    n = vln[vln>0].index.tolist()
    return {'fit': fit, 'n': n}