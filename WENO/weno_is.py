import numpy as np
import pandas as pd
import itertools

###############################################################################
# IS(1,n)
###############################################################################
def is_1(y):
    """

    Parameters
    ----------
    y : y values of the coordenates.

    Returns
    -------
    irs_1 : TYPE
        DESCRIPTION.

    """
    
    initial_scale = int(np.log2(len(y) - 1))
    irs = np.array([], float)
    ir = 0      
    for i in range(len(y)):
        if i < 2:
            irs = np.append(irs, y[i])
        elif 2 <= i < len(y)-1:
            ir = 2**(2*initial_scale - 1)*((y[int(i - 2)] - 2*y[int(i - 1)] + y[int(i)])**2 
                                   + (y[int(i - 1)] - 2*y[int(i)] + y[int(i + 1)])**2 
                                   + 2*(y[int(i - 2)] - 3*y[int(i - 1)] + 3*y[int(i)]
                                        - y[int(i + 1)])**2)
            irs = irs = np.append(irs, ir)
        # elif i == len(y)-1:
        #     irs = np.append(irs, y[i])
        irs_1 = pd.DataFrame(irs)        
    return irs_1

###############################################################################
# IS(2,n)
###############################################################################
def is_2(y):
    """

    Parameters
    ----------
    y : y values of the coordenates.

    Returns
    -------
    irs_2 : TYPE
        DESCRIPTION.

    """
    
    initial_scale = int(np.log2(len(y) - 1))
    irs = np.array([], float)
    ir = 0 
    for i in range(len(y)):
        if i == 0:
            irs = np.append(irs, y[i])
        elif 0 < i < len(y)-2:
            ir = 2**(2*initial_scale - 1)*((y[int(i - 1)] - 2*y[int(i)] + y[int(i + 1)])**2 
                                   + (y[int(i)] - 2*y[int(i + 1)] + y[int(i + 2)])**2 
                                   + 2*(y[int(i - 1)] - 3*y[int(i)] + 3*y[int(i + 1)]
                                        - y[int(i + 2)])**2)
            irs = irs = np.append(irs, ir)
        elif i == len(y)-2:
            irs = np.append(irs, y[i])
        irs_2 = pd.DataFrame(irs)     
    return irs_2

###############################################################################
# IS(3,n)
###############################################################################
def is_3(y):
    """

    Parameters
    ----------
    y : y values of the coordenates.
    
    Returns
    -------
    irs_3 : TYPE
        DESCRIPTION.

    """
    
    initial_scale = int(np.log2(len(y) - 1))
    irs = np.array([], float)
    ir = 0 
    for i in range(len(y)):
        if i < len(y)-3:
            ir = 2**(2*initial_scale - 1)*((y[int(i)] - 2*y[int(i + 1)] + y[int(i + 2)])**2 
                                   + (y[int(i + 1)] - 2*y[int(i + 2)] + y[int(i + 3)])**2 
                                   + 2*(y[int(i)] - 3*y[int(i + 1)] + 3*y[int(i + 2)]
                                        - y[int(i + 3)])**2)
            irs = irs = np.append(irs, ir)
        elif len(y)-3 <= i < len(y)-1:
            irs = np.append(irs, y[i])
        irs_3 = pd.DataFrame(irs)
    return irs_3
