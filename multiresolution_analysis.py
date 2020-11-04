# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 14:05:46 2020

@author: alvar
"""
import numpy as np
import pandas as pd
import itertools
from PPH.PPH_regular import pph_interpolation as pph
from PPH.PPH_regular import pph_interpolation_sign
from scipy.stats.mstats import hmean

###############################################################################
# Compression
###############################################################################
def multiresolution_analysis_compression_number_errors(x, y, number_points):
    '''
    
    Parameters
    ----------
    x : x values of the coordenates.
    y : y values of the coordenates.
    number_points : Number of points to preserve.
    
    Raises
    ------
    Exception
        number_points might be grater than 17.

    Returns
    -------
    ma_compressed_dataset : Coordenates x,y corresponding to scale '4' of 17 points on a regular grid.
    error_save : Coordenates x,y and error values for the 'n' (number_points - 17) biggest error values.
    initial_resolution : Scale of the original dataset.

    '''
    x = np.array(x, float)
    y = np.array(y, float)
    
    error_x = error_y = error_value = np.array([],float)
    number_points = int(number_points)
    final_resolution = int(4)
    initial_resolution = int(np.log2(len(x)))
    
    if number_points < 17:
        raise Exception('number_points might be grater than 17. The value of number_points is: {}'.format(number_points))
        
    else:
        
        for scale in reversed(range(final_resolution, initial_resolution)):
            
            # Decimation operartor. Takes the even values of the dataset introduced, it corresponds to the scale 'j-1'.
            for i, j in zip(x,y):
                x_decimation = x[::2]
                y_decimation = y[::2]
                ma_decimation = pd.DataFrame(np.column_stack((x_decimation,y_decimation)), columns=['x','y'])
            
            # Prediction operator of the values corresponding to the scale 'j-1'.
            ma_prediction_pph = pph.pph_interpolation(ma_decimation['x'], ma_decimation['y'], final_scale = (scale + 1))
            ma_prediction_pph_y = pd.DataFrame(ma_prediction_pph['y'])
            
            #Error values. Difference between the original values and the output of the Prediction operator.
            y = pd.DataFrame(y)
            ma_error = np.array([], float)
            ma_e = 0
    
            for i, j in zip(y, ma_prediction_pph_y):
                ma_e = abs(y[i]) - abs(ma_prediction_pph_y[j])
                ma_error = np.append(ma_error, ma_e)
            
            # y when error > 0
            y = np.array(y, float).reshape(ma_error.shape)
            ma_error_threshold = np.array([], float)
            ma_error_t = 0
            
            for e in range(len(ma_error)):
                if abs(ma_error[e]) != 0:
                    ma_error_t = y[e]
                    ma_error_threshold = np.append(ma_error_threshold, ma_error_t)
                elif abs(ma_error[e]) == 0:
                    ma_error_t = float('NaN')
                    ma_error_threshold = np.append(ma_error_threshold, ma_error_t)
                ma_error_threshold
                
            ma_error_df = pd.DataFrame(np.column_stack((x, y, ma_error_threshold, ma_error)), columns=['x','y','error', 'error_value'])
    
            # x,y and error when error > 0
            ma_x_y_errors = np.array([], float)
            ma_x_y_errors = ma_error_df.dropna(inplace=False)
            error_point_x = np.array(ma_x_y_errors['x'], float)
            error_point_y = np.array(ma_x_y_errors['error'], float)
            error_point = np.array(ma_x_y_errors['error_value'], float)
            
            for e in range(len(error_point_y)):
                error_x = np.append(error_x, error_point_x[e])
                error_y = np.append(error_y, error_point_y[e])
                error_value = np.append(error_value, abs(error_point[e]))
            
            ### New datset of coordinates x,y corresponding to scale 'j-1'
            ma_compressed = ma_error_df[::2].reset_index(drop=True)
            x = np.array(ma_compressed['x'], float)
            y = np.array(ma_compressed['y'], float)
            ma_compressed_dataset = pd.DataFrame(np.column_stack((x, y)), columns=['x','y'])
            
        error_save = pd.DataFrame(np.column_stack((error_x, error_y, error_value)), columns=['error_x','error_y', 'error_value']).sort_values(by='error_value', ascending=False)
        error_save = error_save.head(int(number_points - 17))
        
        return ma_compressed_dataset, error_save, initial_resolution

###############################################################################
# Reconstruction
###############################################################################
def multiresolution_analysis_reconstruction(x, y, error_x, error_y, final_scale):
    '''

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    error_x : TYPE
        DESCRIPTION.
    error_y : TYPE
        DESCRIPTION.
    final_scale : 

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    x = np.array(x, float)
    y = np.array(y, float)
    error_x = np.array(error_x, float)
    error_y = np.array(error_y, float)
    
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
            
            for i, item in enumerate(x):
                
                for e in range(len(error_x)):
                    if item == error_x[e]:
                        y[i] = error_y[e]
    
            ma_reconstructed = pd.DataFrame(np.column_stack((x,y)), columns=['x','y'])
        return ma_reconstructed

    

