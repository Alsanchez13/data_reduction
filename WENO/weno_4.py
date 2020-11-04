import numpy as np
import pandas as pd
import itertools

import WENO
from WENO import weno_is
from WENO import weno_alphas
from WENO import weno_lagrange_schemas
from WENO import weno_4_y

def weno_4(x,y,final_scale):
    """


    Parameters
    ----------
    x : x values of the coordenates.
    y : y values of the coordenates.
    final_scale : Final scale to reach.

    Raises
    ------
    Exception
        final_scale should be greater than initial_scale..

    Returns
    -------
    weno_4 : Coordenates of the WENO 4 subdivision schema at the level of discretization indicated. Data frame type.

    """
    
    x = np.array(x,float)
    y = np.array(y,float)

    final_scale = int(final_scale)
    initial_scale = int(np.log2(len(y)))
    
    if final_scale <= initial_scale:
        raise Exception('final_scale should be greater than initial_scale. The value of initial_scale is: {}'.format(initial_scale))
    
    else:

        for s in range(initial_scale, final_scale):
            irs_1 = irs_2 = irs_3 = np.array([], float)
            alpha_1 = alpha_2 = alpha_3 = np.array([], float)
            y_3_1_odd_points_df = y_2_2_odd_points_df = y_1_3_odd_points_df = np.array([], float)
            weno_4_y_ = np.array([], float)
            xplt = yplt = np.array([], float)
            yp = 0
            nfinal = (2**(s + 1)+1)
            
            #IS
            irs_1 = weno_is.is_1(y)
            irs_2 = weno_is.is_2(y)
            irs_3 = weno_is.is_3(y)
            
            #Alphas
            alpha_1 = weno_alphas.weno_alphas(irs_1, irs_2, irs_3)[0]
            alpha_2 = weno_alphas.weno_alphas(irs_1, irs_2, irs_3)[1]
            alpha_3 = weno_alphas.weno_alphas(irs_1, irs_2, irs_3)[2]
            
            #Schemes
            y_3_1_odd_points_df = weno_lagrange_schemas.s_3_1_odds(y)[0]
            y_2_2_odd_points_df = weno_lagrange_schemas.s_2_2_odds(y)[0]
            y_1_3_odd_points_df = weno_lagrange_schemas.s_1_3_odds(y)[0]
            
            #WENO 4 y odd points
            weno_4_y_ = weno_4_y.weno_4_y(y, y_3_1_odd_points_df, y_2_2_odd_points_df, y_1_3_odd_points_df, alpha_1, alpha_2, alpha_3)
            
            #X and Y
            y = np.array(weno_4_y_)
            x = pd.DataFrame(x)
            xm = x.rolling(2,center=True).mean()
            xplt = pd.concat((x, xm), axis=1).stack().sort_values(ascending=True).reset_index(drop=True).to_numpy()
            x = np.array(xplt)
        x
        y
        x_w = np.array(x).reshape(nfinal,1)
        y_w = np.array(y).reshape(nfinal,1)
        weno_4 = pd.DataFrame(np.column_stack((x_w,y_w)), columns=['x','y'])
        
    return weno_4
