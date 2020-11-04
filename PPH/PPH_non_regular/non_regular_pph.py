# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 10:42:06 2020

@author: alvar
"""
# Libraries
import numpy as np
import pandas as pd

# Function
def non_regular_pph(x, y, iterations):
    '''

    Parameters
    ----------
    x : x values of the coordenates.
    y : y values of the coordenates.
    iterations : number of iterations to perform.

    Raises
    ------
    Exception
        'iterations should be greater than one'.

    Returns
    -------
    non_regular_pph : coordenates (x,y) of the prediction in a Data Frame type.
    
    Coordenates (x,y) of the PPH 4 points interpolatory subdivision scheme on a non-regular grid at the level of discretization indicated.
    Data Frame type.

    '''
    
    x = np.array(x, float)
    y = np.array(y, float)
    
    iterations = int(iterations)
    
    
    if iterations <= 0:
        raise Exception('iterations should be greater than one. The value of iterations is: {}'.format(iterations))
    
    else:
        
        for it in range(1, (iterations + 1)):
            
            # Non regular distances
            H = np.array([],float)
            h = 0
            for i in range(len(x)):
                if i <= 0:
                    H = np.append(H, 1.)
                else: 
                    h = x[i] - x[i - 1]  
                    H = np.append(H, h)
                    
            #Divided differences
            Dj = np.array([], float)
            dj = 0
            for i in range(len(y)):
                if 0 < i < len(y)-2:
                    dj = ((y[int(i - 1)])/(H[int(i)]*(H[int(i)] + H[int(i + 1)]))) -  ((y[int(i)])/(H[int(i)]*H[int(i + 1)])) + ((y[int(i + 1)])/(H[int(i + 1)]*(H[int(i)] + H[int(i + 1)])))
                    Dj = np.append(Dj, dj)   
            
            Dj1 = np.array([], float)
            dj1 = 0
            for i in range(len(y)):
                if 0 < i < len(y)-2:
                    dj1 = ((y[int(i)])/(H[int(i + 1)]*(H[int(i + 1)] + H[int(i + 2)]))) -  ((y[int(i + 1)])/(H[int(i + 1)]*H[int(i + 2)])) + ((y[int(i + 2)])/(H[int(i + 2)]*(H[int(i + 1)] + H[int(i + 2)])))
                    Dj1 = np.append(Dj1, dj1) 
            
            #Weights
            wj = np.array([], float)
            w = 0
            for i in range(len(H)):
                if 0 < i < len(H)-2:
                    w = (H[int(i + 1)] + 2*H[int(i + 2)])/(2*(H[int(i)] + H[int(i + 1)] + H[int(i + 2)]))
                    wj = np.append(wj, w)    
                for j in range(len(wj)):
                    wj1 = 1.0 - wj
            
            # Weighted Harmonic mean over non-uniform grid
            V = np.array([], float)
            v = 0
            for j in range(len(Dj1)):
                if Dj[j]*Dj1[j] > 0:
                    v = ((Dj[j]*Dj1[j])/(wj[j]*Dj1[j] + wj1[j]*Dj[j]))
                    V = np.append(V, v)
                else:
                    v = 0.0
                    V = np.append(V, v)
                    
            # X new values to predict, placed in the middle position between original x values.
            X = np.array([], float)
            xm = 0
            for j in range(len(x)):
                if 0 < j < len(x)-2:    
                    xm = (x[j] + x[j + 1])/2.0
                    X = np.append(X, xm)
                    
            # Coeffcients + PPH
            pph_ni = np.array([], float)
            pph_ni = np.append(pph_ni, y[0])
            yp = (1/2)*y[int(0)] + (1/2)*y[int(1)]
            pph_ni = np.append(pph_ni, yp)
            pph = a0_1 = a1_1 = a2_1 = a3_1 = a0_2 = a1_2 = a2_2 = a3_2 = 0
            for i in range(len(y)):
                if i < len(V):
                    if abs(Dj[i]) <= abs(Dj1[i]):
                        # Case 1
                        a0_1 = ((y[int(i + 1)] + y[int(i + 2)])/(2)) - (((H[int(i + 2)]**2)/(4))*V[i])
                        a1_1 = (-y[int(i + 1)] + y[int(i + 2)])/(H[int(i + 2)]) + (((H[int(i + 2)]**2)/(4*H[int(i + 1)] + 2*H[int(i + 2)]))*(Dj[i] - V[i]))
                        a2_1 = V[i]
                        a3_1 = (-(2)/(2*H[int(i + 1)] + H[int(i + 2)]))*(Dj[i] - V[i])
                        pph = a0_1 + (a1_1 * (X[int(i)] - ((x[int(i + 1)] + x[int(i + 2)])/2))) + (a2_1 * ((X[int(i)] - ((x[int(i + 1)] + x[int(i + 2)])/2))**2)) + (a3_1 * ((X[int(i)] - ((x[int(i + 1)] + x[int(i + 2)])/2))**3))
                        pph_ni = np.append(pph_ni, y[int(i + 1)])  
                        pph_ni = np.append(pph_ni, pph)  
                        
                    else:
                        # Case 2
                        a0_2 = ((y[int(i + 1)] + y[int(i + 2)])/(2)) - (((H[int(i + 2)]**2)/(4))*V[i])
                        a1_2 = ((-y[int(i + 1)] + y[int(i + 2)])/(H[int(i + 2)])) + (((H[int(i + 2)]**2)/(2*H[int(i + 2)] + 4*H[int(i + 3)]))*(-Dj1[i] + V[i]))
                        a2_2 = V[i]
                        a3_2 = -(((2)/(H[int(i + 2)]) + 2*H[int(i + 3)])*(-Dj1[i] + V[i]))
                        pph = a0_2 + (a1_2 * (X[int(i)] - ((x[int(i + 1)] + x[int(i + 2)])/2))) + (a2_2 * ((X[int(i)] - ((x[int(i + 1)] + x[int(i + 2)])/2))**2)) + (a3_2 * ((X[int(i)] - ((x[int(i + 1)] + x[int(i + 2)])/2))**3))
                        pph_ni = np.append(pph_ni, y[int(i + 1)])
                        pph_ni = np.append(pph_ni, pph)
            pph_ni = np.append(pph_ni, y[int(len(y)-2)])
            yp = (1/2)*y[int(len(y)-2)] + (1/2)*y[int(len(y)-1)]
            pph_ni = np.append(pph_ni, yp)
            pph_ni = np.append(pph_ni, y[int(len(y)-1)])
            
            # Joining the x
            x_pph = np.array([], float)
            x_pph = np.append(x_pph, x[0])
            xm = (x[0] + x[1])/2.0
            x_pph = np.append(x_pph, xm)
            for j in range(len(X)):
                if j < len(X):
                    x_pph = np.append(x_pph, x[int(j + 1)])
                    x_pph = np.append(x_pph, X[j])
            x_pph = np.append(x_pph, x[int(len(x)-2)])
            xm = (x[int(len(x)-2)] + x[int(len(x)-1)])/2.0
            x_pph = np.append(x_pph, xm)
            x_pph = np.append(x_pph, x[int(len(x)-1)])
            
            # Returning the x and y
            x = x_pph
            y = pph_ni
            non_regular_pph = pd.DataFrame(np.column_stack((x_pph, pph_ni)), columns=['x','y'])
        return non_regular_pph
