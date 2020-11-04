import numpy as np
import pandas as pd
import itertools
from scipy.stats.mstats import hmean
from PPH.PPH_regular import pph_interpolation_sign


def pph_interpolation(x,y,final_scale):
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
    pph_interpolation : Coordenates of the PPH 4 points interpolatory subdivision scheme at the level of discretization indicated over a regular grid.
    To build the values x=1 and x=n-1, Lagrange interpolatory subdivision scheme S(1,3) and S(3,1) respectively are applied.
    Data Frame type.

    """
    x = np.array(x, float)
    y = np.array(y, float)
    
    final_scale = int(final_scale)
    initial_scale = int(np.log2(len(x)))
    
    if final_scale <= initial_scale:
        raise Exception('final_scale should be greater than initial_scale. The value of initial_scale is: {}'.format(initial_scale))
    
    else:

        for scale in range(initial_scale, final_scale):
            xplt = np.array([], float)
            yplt = np.array([], float)
            yp = sign = 0
            nfinal = (2**(scale+1)+1)
            
            for i in range(len(y)):
                if i == 0:
                    yp = (5/16)*y[int(i)] + (15/16)*y[int(i + 1)] - (5/16)*y[int(i + 2)] + (1/16)*y[int(i + 3)]
                    yplt = np.append(yplt, y[i])
                    yplt = np.append(yplt, yp)
                elif 0 < i < len(y)-2:
                    sign = ((pph_interpolation_sign.sign_pph(y[int(i - 1)] - 2*y[int(i)] + y[int(i + 1)]) + pph_interpolation_sign.sign_pph(y[int(i)] - 2*y[int(i + 1)] + y[int(i + 2)])) / 2)
                    yp = (y[int(i)] + y[int(i + 1)])/2 - (1/8)*sign*hmean(np.array([abs(y[int(i - 1)] - 2*y[int(i)] + y[int(i + 1)]), abs(y[int(i)] - 2*y[int(i + 1)] +  y[int(i + 2)])], float))
                    yplt = np.append(yplt, y[i])
                    yplt = np.append(yplt, yp)
                elif i == len(y)-2:
                    yp = (1/16)*y[int(i - 2)] - (5/16)*y[int(i - 1)] + (15/16)*y[int(i)] + (5/16)*y[int(i + 1)]
                    yplt = np.append(yplt, y[i])
                    yplt = np.append(yplt, yp)
                else:
                    yplt = np.append(yplt, y[i])

            for i,j in zip(x,y):
                x = pd.DataFrame(x)
                xm = x.rolling(2,center=True).mean()
                xplt = pd.concat([x, xm], axis=1).stack().sort_values(ascending=True).reset_index(drop=True).to_numpy().reshape(nfinal,1)
            x = xplt
            yplt = np.array([yplt], float)
            y = yplt.reshape(nfinal,1)
            pph_interpolation = pd.DataFrame(np.column_stack((x,y)), columns=['x','y'])
        return pph_interpolation
