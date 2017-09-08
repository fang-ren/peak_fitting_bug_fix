"""
Created on Feb 8, 2017

@author: fangren
"""
'''
edited by Jae
I am gonig to set the script up so that 
it churns through an entire experiment from 
the January beamtime, via stupid, stupid
copy and paste

Re-evaluating the script for May reveals that a memory leak prevents this from being easily automated.  It seems to benefit
from being able to stop in between calls.

'''

import os.path
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import basinhopping
from scipy.interpolate import interp1d
from numpy.polynomial.chebyshev import chebval
global Q
global intensity



def set_window(windows_size, Q_original, intensity_original):
    global Q
    global intensity
    Q = []
    intensity = []
    for i in np.arange(0, len(Q_original), windows_size):
        Q.append(Q_original[i])
        intensity.append(intensity_original[i])
    Q = np.array(Q)
    intensity = np.array(intensity)
    return Q, intensity

def file_index(index):
    if len(str(index)) == 1:
        return '000' + str(index)
    elif len(str(index)) == 2:
        return '00' + str(index)
    elif len(str(index)) == 3:
        return '0' + str(index)
    elif len(str(index)) == 4:
        return str(index)

def func(x, *params):
    params = params[0]
    y = chebval(x, params[:4])
    E = params[4]
    y = y + E /x
    return y

def object_func(*params):
    params = params[0]
    J = 0
    fit = chebval(Q, params[:4])
    E = params[4]
    fit = fit + E/Q
    for i in range(len(intensity)):
        if Q[i] < 1:
            if intensity[i] < fit[i]:
                J = J + (intensity[i]-fit[i])**4
            elif intensity[i] >= fit[i]:
                J = J + (intensity[i]-fit[i])**4
        else:
            if intensity[i] < fit[i]:
                J = J + (intensity[i] - fit[i]) ** 4
            elif intensity[i] >= fit[i]:
                J = J + (intensity[i] - fit[i]) ** 2
    return J

def background_subtraction(path, basename, start_index, end_index):
    save_path = path + '\\subtracted\\'
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    index = start_index
    total_num_scan = end_index
    while (index <= total_num_scan):
        file_name = os.path.join(path, basename + file_index(index) + '_1D.csv')
        print 'processing', file_name
        if os.path.exists(file_name):
            spectrum = np.genfromtxt(file_name, delimiter= ',')
            intensity_original = spectrum[:, 1][30:-70]
            Q_original =  spectrum[:, 0][30:-70]

            # create a sparse set of data for fitting, window size: 20
            Q, intensity = set_window(20, Q_original, intensity_original)

            # fit  the result according the the defined object_func with the initialized parameters and bounds
            x0 = [1, 1, 1, 1, 1]
            result = basinhopping(object_func, x0)
            bckgrd_sparse = func(Q, result.x)
            f = interp1d(Q, bckgrd_sparse, kind='cubic', bounds_error=False)
            bckgrd = f(Q_original)

            # save the plot
            plt.subplot(211)
            plt.plot(Q_original, intensity_original)
            plt.plot(Q, intensity, 'o')
            plt.plot(Q, bckgrd_sparse, 'o', c = 'r')
            plt.xlim(0.8, 5.8)
            plt.subplot(212)
            plt.plot(Q_original, intensity_original-bckgrd)
            plt.plot(Q, [0] * len(Q), 'r--')
            plt.xlim(0.8, 5.8)
            plt.savefig(os.path.join(save_path, basename + file_index(index) + '_bckgrd_subtracted.png'))
            plt.close('all')
            np.savetxt(os.path.join(save_path, basename + file_index(index) + '_bckgrd_subtracted.csv'), np.concatenate(([Q_original], [intensity_original - bckgrd])).T, delimiter= ',')
        index += 1
