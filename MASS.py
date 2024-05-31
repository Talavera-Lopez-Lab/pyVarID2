import numpy as np
import scipy.special as sp

def theta_ml(y, mu, n=None, weights=None, limit=10, eps=np.finfo(float).eps**0.25, trace=False):
    def score(n, th, mu, y, w):
        return np.sum(w * (sp.digamma(th + y) - sp.digamma(th) + np.log(th) + 1 - np.log(th + mu) - (y + th) / (mu + th)))

    def info(n, th, mu, y, w):
        return np.sum(w * (-sp.polygamma(1, th + y) + sp.polygamma(1, th) - 1/th + 2/(mu + th) - (y + th) / (mu + th)**2))
    
    if hasattr(y, 'fittedvalues'):
        mu = y.fittedvalues
        y = y.endog if y.endog is not None else mu + y.resid_response
    
    if weights is None:
        weights = np.ones_like(y)
    
    if n is None:
        n = np.sum(weights)
    
    t0 = n / np.sum(weights * (y / mu - 1)**2)
    it = 0
    del_val = 1
    
    if trace:
        print(f"theta_ml: iter {it} 'theta = {t0}'")
    
    while (it := it + 1) < limit and abs(del_val) > eps:
        t0 = abs(t0)
        del_val = score(n, t0, mu, y, weights) / (i := info(n, t0, mu, y, weights))
        t0 += del_val
        if trace:
            print(f"theta_ml: iter {it} theta = {t0}")
    
    if t0 < 0:
        t0 = 0
        print("estimate truncated at zero")
        t0 = 0
    
    if it == limit:
        print("iteration limit reached")
    
    se = np.sqrt(1 / i)
    return t0, se
