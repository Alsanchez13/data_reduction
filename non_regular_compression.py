# -*- coding: utf-8 -*-
"""
Created on Sun Aug 30 15:37:50 2020

@author: alvar
"""
# Libraries
import numpy as np
import pandas as pd
import math

###############################################################################
# Compression
###############################################################################
def non_regular_compression(x, y, compression):
    '''

    Parameters
    ----------
    x : x values of the coordenates.
    y : y values of the coordenates.
    compression : Percentage of the origianl dataset to stock. The choices are: 50%, 25%, 12.5%, 6.3%, 3.2% and 1.6%.

    Raises
    ------
    Exception
        compression might be a value among: 0.50, 0.25, 0.125, 0.063, 0.032, 0.016.

    Returns
    -------
    compressed_data : Coordenates (x,y) corresponding to the dataset compressed.The values associeted to errors are stocked from 9 points onwards. At every loop the value associated to maximal error is stocked. The operator performing the prediction is the adaptation of PPH to non-regular grids.
    iterations : Number of iterations needed to reconstrunct the curves up to the original data size.

    '''
    
    x_original = np.array(x, float)
    y_original = np.array(y, float)
    
    compression = float(compression)
    initial_points = int(len(x_original))
    compressed_points = int(math.ceil(compression * initial_points))
    
    if compression not in [0.50, 0.25, 0.125, 0.063, 0.032, 0.016]:
        raise Exception('compression might be a value among: 0.50, 0.25, 0.125, 0.063, 0.032, 0.016. The value of compression is: {}'.format(compression))
        
    else:
        # Decimation operartor. Takes even values of the dataset introduced until the number of items is below the compressed points wished.
        x = x_original[::]
        y = y_original[::]
        iterations = 0
        while len(x) >= 10:
            x = x[::2]
            y = y[::2]
            iterations += 1 
        x = np.append(x, x_original[1])
        y = np.append(y, y_original[1])
        x = np.append(x, x_original[-2])
        y = np.append(y, y_original[-2])
        x, y = zip(*sorted(zip(x, y)))
        x_decimation = np.array(x[::])
        y_decimation = np.array(y[::])
        pph_x = pph_ni = x_errors = y_errors = np.array([], float)
        
        for itr in range(compressed_points - len(x_decimation)):
            x = np.array(x_decimation[::])
            y = np.array(y_decimation[::])
            pph_x = pph_ni = np.array([], float)
            pph_x = np.append(pph_x, x[0])
            pph_ni = np.append(pph_ni, y[0])
            
            for i in range(1,len(y_decimation)):
                h0 = h1 = h2 =dj = dj1 = wj = wj1 = v = 0
                if i == 0:
                    x = np.append(x, x_original[1])
                    y = np.append(y, y_original[1])
                    x = np.append(x, x_original[-2])
                    y = np.append(y, y_original[-2])
                    x, y = zip(*sorted(zip(x, y)))
                    x = np.array(x[::])
                    y = np.array(y[::])
                    pph_x = np.append(pph_x, x[0])
                    pph_ni = np.append(pph_ni, y[0])
                   
                elif 0 < i < len(y)-2:
                    # Non regular distances
                    h0 = x[i] - x[i-1]
                    h1 = x[i+1] - x[i]
                    h2 = x[i+2] - x[i+1]
                   
                    #Divided differences
                    dj = (y[i - 1])/(h0 *(h0 + h1)) -  (y[i])/(h0*h1) + (y[i + 1])/(h1*(h0 + h1))
                    dj1 = ((y[i])/(h1*(h1 + h2))) - ((y[i+1])/(h1*h2)) + ((y[i + 2])/(h2*(h1 + h2)))
                   
                    #Weights
                    wj = (h1 + 2*h2)/(2*(h0 + h1 + h2))
                    wj1 = 1.0 - wj
                    
                    # Weighted Harmonic mean over non-uniform grid
                    if dj*dj1 > 0:
                        v = ((dj*dj1)/(wj*dj1 + wj1*dj))
                    else:
                        v = 0.0
                        
                    # X evaluated in a value of the original x, placed in the middle position of the interval [i, i + n] where n = 2^iterations.    
                    X = np.array([], float)
                    xm = 0    
                    for k in range(len(x_original)-2):
                        if x[i] < x_original[k] < x[i+1]:
                            xm = x_original[k]
                            X = np.append(X, xm)
                
                    # Coeffcients + PPH
                    pph_ni = np.append(pph_ni, y[i])
                    pph = a0_1 = a1_1 = a2_1 = a3_1 = a0_2 = a1_2 = a2_2 = a3_2 = 0
                    for p in range(len(X)):
                        if abs(dj) <= abs(dj1):
                            # Case 1
                            a0_1 = ((y[i] + y[i + 1])/(2)) - (((h1**2)/(4))*v)
                            a1_1 = (-y[i] + y[i + 1])/(h1) + (((h1**2)/(4*h0 + 2*h1))*(dj - v))
                            a2_1 = v
                            a3_1 = (-(2)/(2*h0 + h1))*(dj - v)
                            pph = a0_1 + (a1_1 * (X[p] - ((x[i] + x[i + 1])/2))) + (a2_1 * ((X[p] - ((x[i] + x[i + 1])/2))**2)) + (a3_1 * ((X[p] - ((x[i] + x[i + 1])/2))**3))
                            pph_ni = np.append(pph_ni, pph)
                        else:
                            # Case 2
                            a0_2 = ((y[i] + y[i + 1])/(2)) - (((h1**2)/(4))*v)
                            a1_2 = ((-y[i] + y[i + 1])/(h1)) + (((h1**2)/(2*h1 + 4*h2))*(-dj1 + v))
                            a2_2 = v
                            a3_2 = -(((2)/(h1) + 2*h2)*(-dj1 + v))
                            pph = a0_2 + (a1_2 * (X[p] - ((x[i] + x[i])/2))) + (a2_2 * ((X[p] - ((x[i] + x[i + 1])/2))**2)) + (a3_2 * ((X[p] - ((x[i] + x[i + 1])/2))**3))
                            pph_ni = np.append(pph_ni, pph)
                    # Joining the x
                    pph_x = np.append(pph_x, x[i])
                    pph_x = np.append(pph_x, X)    
                    
                elif i == len(y)-1:
                    pph_x = np.append(pph_x, x[i])
                    pph_ni = np.append(pph_ni, y[i])
                else:
                    pph_x = np.append(pph_x, x[i])
                    pph_ni = np.append(pph_ni, y[i])

            # Aboslute error values. Difference between the original values and the output of the Prediction operator.
            error = np.array([], float)
            e = 0
            for r in range(len(y_original)):
                e = y_original[r] - pph_ni[r]
                e = abs(e)
                error = np.append(error, e)
           
            ### Option 2: Saving the maximal error per loop. Joining the y values and the error to save the maximal error.
            x_e = y_e = 0
            for e in range(len(error)):
                if error[e] == max(error):       
                    y_e = y_original[e]
                    y_errors = np.append(y_errors, y_e)
                    x_e = x_original[e]
                    x_errors = np.append(x_errors, x_e)
                 
            # Saving the maximal error
            x_decimation = np.append(x_decimation, x_e)
            y_decimation = np.append(y_decimation, y_e)
            x_decimation, y_decimation = zip(*sorted(zip(x_decimation, y_decimation)))
        # Dataset with the error joined.
        x_decimation = np.array(x_decimation[::])
        y_decimation = np.array(y_decimation[::]) 
    # Final dataset.
    iterations = int(iterations - 1)
    compressed_data = pd.DataFrame(np.column_stack((x_decimation, y_decimation)), columns=['x','y'])
        
    return compressed_data, iterations
###############################################################################
# Reconstruction
###############################################################################
def non_regular_reconstrunction_original_values(x, y, x_original):
    '''

    Parameters
    ----------
    x : x values of the coordenates.
    y : y values of the coordenates.
    x_original : x values corresponding to the original data-set before compression.

    Returns
    -------
    non_regular_reconstructed : (x,y) values corresponding to the reconstructed curve. The position on the x axis for the new points corresponds to the original x values. Reconstruction based in PPH scheme on non-regular grids.

    '''
    
    x = np.array(x, float)
    y = np.array(y, float)
    x_original = np.array(x_original, float)
    pph_x = pph_ni = np.array([], float)
    pph_x = np.append(pph_x, x[0])
    pph_ni = np.append(pph_ni, y[0])
    
    for i in range(1,len(y)):
    
        h0 = h1 = h2 =dj = dj1 = wj = wj1 = v = 0
        if i == 0:
            pph_x = np.append(pph_x, x[0])
            pph_ni = np.append(pph_ni, y[0])
           
        elif 0 < i < len(y)-2:
            # Non regular distances
            h0 = x[i] - x[i-1]
            h1 = x[i+1] - x[i]
            h2 = x[i+2] - x[i+1]
           
            #Divided differences
            dj = (y[i - 1])/(h0 *(h0 + h1)) -  (y[i])/(h0*h1) + (y[i + 1])/(h1*(h0 + h1))
            dj1 = ((y[i])/(h1*(h1 + h2))) - ((y[i+1])/(h1*h2)) + ((y[i + 2])/(h2*(h1 + h2)))
    
            #Weights
            wj = (h1 + 2*h2)/(2*(h0 + h1 + h2))
            wj1 = 1.0 - wj
            
            # Weighted Harmonic mean over non-uniform grid
            if dj*dj1 > 0:
                v = ((dj*dj1)/(wj*dj1 + wj1*dj))
            else:
                v = 0.0
           
            # X evaluated in a value of the original x, placed in the middle position of the interval [i, i + n] where n = 2^iterations.    
            X = np.array([], float)
            xm = 0    
            for k in range(len(x_original)-2):
                if x[i] < x_original[k] < x[i+1]:
                    xm = x_original[k]
                    X = np.append(X, xm)
        
            # Coeffcients + PPH
            pph_ni = np.append(pph_ni, y[i])
            pph = a0_1 = a1_1 = a2_1 = a3_1 = a0_2 = a1_2 = a2_2 = a3_2 = 0
            for p in range(len(X)):
                if abs(dj) <= abs(dj1):
                    # Case 1
                    a0_1 = ((y[i] + y[i + 1])/(2)) - (((h1**2)/(4))*v)
                    a1_1 = (-y[i] + y[i + 1])/(h1) + (((h1**2)/(4*h0 + 2*h1))*(dj - v))
                    a2_1 = v
                    a3_1 = (-(2)/(2*h0 + h1))*(dj - v)
                    pph = a0_1 + (a1_1 * (X[p] - ((x[i] + x[i + 1])/2))) + (a2_1 * ((X[p] - ((x[i] + x[i + 1])/2))**2)) + (a3_1 * ((X[p] - ((x[i] + x[i + 1])/2))**3))
                    pph_ni = np.append(pph_ni, pph)
                else:
                    # Case 2
                    a0_2 = ((y[i] + y[i + 1])/(2)) - (((h1**2)/(4))*v)
                    a1_2 = ((-y[i] + y[i + 1])/(h1)) + (((h1**2)/(2*h1 + 4*h2))*(-dj1 + v))
                    a2_2 = v
                    a3_2 = -(((2)/(h1) + 2*h2)*(-dj1 + v))
                    pph = a0_2 + (a1_2 * (X[p] - ((x[i] + x[i])/2))) + (a2_2 * ((X[p] - ((x[i] + x[i + 1])/2))**2)) + (a3_2 * ((X[p] - ((x[i] + x[i + 1])/2))**3))
                    pph_ni = np.append(pph_ni, pph)
            # Joining the x
            pph_x = np.append(pph_x, x[i])
            pph_x = np.append(pph_x, X)    
            
        elif i == len(y)-1:
            pph_x = np.append(pph_x, x[i])
            pph_ni = np.append(pph_ni, y[i])
        else:
            pph_x = np.append(pph_x, x[i])
            pph_ni = np.append(pph_ni, y[i])      
    non_regular_reconstructed = pd.DataFrame(np.column_stack((pph_x, pph_ni)), columns=['x','y'])
    
    return non_regular_reconstructed
    
    
