import numpy as np
import pandas as pd
import itertools

###############################################################################
# WENO_4_y
###############################################################################
def weno_4_y(y, y_3_1_odd_points_df, y_2_2_odd_points_df, y_1_3_odd_points_df, alpha_1, alpha_2, alpha_3):
    """
    
    Parameters
    ----------
    y : y values of the coordenates.
    y_3_1_odd_points_df : Odd values of the Lagrange interpolatory schema S(3,1).
    y_2_2_odd_points_df : Odd values of the Lagrange interpolatory schema S(3,1).
    y_1_3_odd_points_df : Odd values of the Lagrange interpolatory schema S(3,1).
    alpha_1 : First stencil's weight.
    alpha_2 : Second stencil's weight.
    alpha_3 : Thrid stencil's weight.

    Returns
    -------
    weno_4_y : y values of the WENO 4 subdivision schema.

    """
    
    weno_odd = np.array([], float)
    weno = 0
    for i, j, k in zip(y_3_1_odd_points_df, y_2_2_odd_points_df, y_1_3_odd_points_df):  
        for a1, a2, a3 in zip(alpha_1, alpha_2, alpha_3):
            weno = alpha_1[a1]*y_3_1_odd_points_df[i] + alpha_2[a2]*y_2_2_odd_points_df[j] + alpha_3[a3]*y_1_3_odd_points_df[k]
            weno_odd = np.append(weno_odd, weno)
        
    weno_y = [None]*(len(y)+len(weno_odd))
    weno_y[::2] = y
    weno_y[1::2] = weno_odd
    weno_4_y = np.array(weno_y)
    
    return weno_4_y
