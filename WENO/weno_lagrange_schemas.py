import numpy as np
import pandas as pd
import itertools

###############################################################################
# S(3,1)_2n+1
###############################################################################
def s_3_1_odds(y):
    """

    Parameters
    ----------
    y : y values of the coordenates.

    Returns
    -------
    y_3_1_odd_points_df : Odd values of the Lagrange interpolatory schema S(3,1). Data frame type.
    y_3_1_df : y values of the Lagrange interpolatory schema S(3,1). Data frame type.

    """
    
    yplt = np.array([], float)
    yp = 0
    for i in range(len(y)):
        if i == 0:
            yp = (5/16)*y[int(i)] + (15/16)*y[int(i + 1)] - (5/16)*y[int(i + 2)] + (1/16)*y[int(i + 3)]
            yplt = np.append(yplt, y[i])
            yplt = np.append(yplt, yp)
        elif i == 1:
            yp = -(1/16)*y[int(i - 1)] + (9/16)*y[int(i)] + (9/16)*y[int(i + 1)] - (1/16)*y[int(i + 2)]
            yplt = np.append(yplt, y[i])
            yplt = np.append(yplt, yp)
        elif 1 < i < len(y)-1:
            yp = (1/16)*y[int(i - 2)] - (5/16)*y[int(i - 1)] + (15/16)*y[int(i)] + (5/16)*y[int(i + 1)]
            yplt = np.append(yplt, y[i])
            yplt = np.append(yplt, yp)
        else:
            yplt = np.append(yplt, y[i])
            y_3_1 = np.array(yplt).transpose()
            y_3_1_odd_points = y_3_1[1::2]
            y_3_1_df = pd.DataFrame(y_3_1)
            y_3_1_odd_points_df = pd.DataFrame(y_3_1_odd_points)
    return y_3_1_odd_points_df, y_3_1_df

###############################################################################
# S(2,2)_2n+1
###############################################################################
def s_2_2_odds(y):
    """

    Parameters
    ----------
    y : y values of the coordenates.

    Returns
    -------
    y_2_2_odd_points_df : Odd values of the Lagrange interpolatory schema S(2,2). Data frame type.
    y_2_2_df : y values of the Lagrange interpolatory schema S(2,2). Data frame type.

    """
    
    yplt = np.array([], float)
    yp = 0
    for i in range(len(y)):
        if i == 0:
            yp = (5/16)*y[int(i)] + (15/16)*y[int(i + 1)] - (5/16)*y[int(i + 2)] + (1/16)*y[int(i + 3)]
            yplt = np.append(yplt, y[i])
            yplt = np.append(yplt, yp)
        elif 0 < i < len(y)-2:
            yp = -(1/16)*y[int(i - 1)] + (9/16)*y[int(i)] + (9/16)*y[int(i + 1)] - (1/16)*y[int(i + 2)]
            yplt = np.append(yplt, y[i])
            yplt = np.append(yplt, yp)
        elif i == len(y)-2:
            yp = (1/16)*y[int(i - 2)] - (5/16)*y[int(i - 1)] + (15/16)*y[int(i)] + (5/16)*y[int(i + 1)]
            yplt = np.append(yplt, y[i])
            yplt = np.append(yplt, yp)
        else:
            yplt = np.append(yplt, y[i])
            y_2_2 = np.array(yplt).transpose()
            y_2_2_odd_points = y_2_2[1::2]
            y_2_2_df = pd.DataFrame(y_2_2)
            y_2_2_odd_points_df = pd.DataFrame(y_2_2_odd_points)
    return y_2_2_odd_points_df, y_2_2_df

###############################################################################
# S(1,3)_2n+1
###############################################################################
def s_1_3_odds(y):
    """

    Parameters
    ----------
    y : y values of the coordenates.

    Returns
    -------
    y_1_3_odd_points_df : Odd values of the Lagrange interpolatory schema S(1,3). Data frame type.
    y_1_3_df : y values of the Lagrange interpolatory schema S(1,3). Data frame type.

    """
    
    yplt = np.array([], float)
    yp = 0
    for i in range(len(y)):
        if i < len(y)-3:
            yp = (5/16)*y[int(i)] + (15/16)*y[int(i + 1)] - (5/16)*y[int(i + 2)] + (1/16)*y[int(i + 3)]
            yplt = np.append(yplt, y[i])
            yplt = np.append(yplt, yp)
        elif i == len(y)-3:
            yp = -(1/16)*y[int(i - 1)] + (9/16)*y[int(i)] + (9/16)*y[int(i + 1)] - (1/16)*y[int(i + 2)]
            yplt = np.append(yplt, y[i])
            yplt = np.append(yplt, yp)
        elif i == len(y)-2:
            yp = (1/16)*y[int(i - 2)] - (5/16)*y[int(i - 1)] + (15/16)*y[int(i)] + (5/16)*y[int(i + 1)]
            yplt = np.append(yplt, y[i])
            yplt = np.append(yplt, yp)
        else:
            yplt = np.append(yplt, y[i])
            y_1_3 = np.array(yplt).transpose()
            y_1_3_odd_points = y_1_3[1::2]
            y_1_3_df = pd.DataFrame(y_1_3)
            y_1_3_odd_points_df = pd.DataFrame(y_1_3_odd_points)   
    return y_1_3_odd_points_df, y_1_3_df
