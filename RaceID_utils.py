import numpy as np
from scipy.stats import linregress

def fitbackground(x, mthr=-1):
    '''this one needs some work still not doing as in the original R source code:
    fitbackground <- function(x,mthr=-1){
    x <- as.matrix(x)
    m <- rowMeans(x)
    v <- rowVars(x)
  
    ml <- log2(m)
    vl <- log2(v)
    f <- ml > -Inf & vl > -Inf
    ml <- ml[f]
    vl <- vl[f]
    mm <- -8
    repeat{
        fit <- lm(vl ~ ml + I(ml^2)) 
        if( coef(fit)[3] >= 0 | mm >= mthr){
            break
        }
        mm <- mm + .5
        f <- ml > mm
        ml <- ml[f]
        vl <- vl[f]
    }
    
    vln <- log2(v)  - log2(sapply(m,FUN=uvar,fit=fit))
    n <- names(m)[vln>0]
    return(list(fit=fit,n=n))
}
'''
    x = np.array(x)
    m = np.mean(x, axis=1)
    v = np.var(x, axis=1)

    ml = np.log2(m)
    vl = np.log2(v)
    f = (ml > -np.inf) & (vl > -np.inf)
    ml = ml[f]
    vl = vl[f]
    mm = -8
    while True:
        fit = linregress(ml, vl)
        if fit[3] >= 0 or mm >= mthr:
            break
        mm += 0.5
        f = ml > mm
        ml = ml[f]
        vl = vl[f]

    def uvar(x, fit):
        return 2 ** (fit[1] + np.log2(x) * fit[0] + fit[2] * (np.log2(x) ** 2))

    vln = np.log2(v) - np.log2(np.apply_along_axis(uvar, 0, m, fit=fit))
    n = [name for name, val in zip(np.arange(len(m)), vln) if val > 0]
    return {'fit': fit, 'n': n}
