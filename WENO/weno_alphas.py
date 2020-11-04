import numpy as np
import pandas as pd
import itertools

###############################################################################
# Alpha
###############################################################################
def weno_alphas(irs_1, irs_2 , irs_3):
    """
    Parameters
    ----------
    irs_1 : IS(1,n) Regularity index of the scheme.
    irs_2 : IS(2,n) Regularity index of the scheme.
    irs_3 : IS(3,n) Regularity index of the scheme.

    Returns
    -------
    alpha_1 : First stencil's weight.
    alpha_2 : Second stencil's weight.
    alpha_3 : Thrid stencil's weight.
    
    Wehere the sum of alphas is equal to one.

    """
    C = np.array([(3/16), (10/16), (3/16)])
    
    for i, j, k in zip(irs_1, irs_2, irs_3):
        bi_1 = bi_2 = bi_3 = beta_j = np.array([], float)
        b1 = b2 = b3 = b_j = 0
        for c in range(len(C)):
            if c == 0:
                b1 = C[0]/(((10**-6) + irs_1[i])**2)
                bi_1 = np.append(bi_1, b1)
            elif c == 1:
                b2 = C[1]/(((10**-6) + irs_2[j])**2)
                bi_2 = np.append(bi_2, b2)
            else:
                b3 = C[2]/(((10**-6) + irs_3[k])**2)
                bi_3 = np.append(bi_3, b3)
                bi_1 = pd.DataFrame(bi_1)
                bi_2 = pd.DataFrame(bi_2)
                bi_3 = pd.DataFrame(bi_3)
        
        for i, j, k in zip(bi_1, bi_2, bi_3):
            b_j = bi_1[i] + bi_2[j] + bi_3[k]
            beta_j = np.append(beta_j, b_j)
            beta_j = pd.DataFrame(beta_j)
            
            alpha_1 = bi_1[0] / beta_j[0]
            alpha_1 = pd.DataFrame(alpha_1)
            alpha_2 = bi_2[0] / beta_j[0]
            alpha_2 = pd.DataFrame(alpha_2)
            alpha_3 = bi_3[0] / beta_j[0]
            alpha_3 = pd.DataFrame(alpha_3)
            
    return alpha_1, alpha_2, alpha_3