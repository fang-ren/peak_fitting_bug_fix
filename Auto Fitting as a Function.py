# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 15:37:32 2017
This program takes the functions that Fang and Travis built (with some of 
my feedback) and puts them into a single call program that fits and extracts
the FWHM from all subfolders in a specified directory.

The way that this program is supposed to work is as follows:
    
    1.) Import all of the needed packages (including some we cobbled together)
    2.) Look in the user set path for folders 
        a.) create a list of those directories and get its length
    3.) That is used to fed up a for loop that iterates through all sub-directories
    one by one and:
        a.) performs the background subtraction using Fang's function
        b.) performs the peak fitting on the bckg subtracted files
        c.) extracts the max FWHM and min FWHM from those fits and 
        places them in a single file with a structure similar to combiview

Current known problems include:
    a ValueError during peak fitting with an unknown origin

                    


@author: jrh8
"""
import scipy
from peak_detection_matlab import peakdet
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from os.path import basename
from scipy.optimize import curve_fit
import pandas as pd
from bckgrd_subtract_functions import background_subtraction
from peak_fitting_auto_with_FWHM_extraction_error_handling_as_functions import XRDfitting
from peak_fitting_auto_with_FWHM_extraction_error_handling_as_functions import extract_max_FWHM
from peak_fitting_auto_with_FWHM_extraction_error_handling_as_functions import extract_min_FWHM
#%%
path = 'test_data\\'

background_subtraction(path, 'SampleB2_7_400C_24x24_t30_', 389, 390)

path2 = path + '\\subtracted\\'
XRDfitting(path2)
#%%
extract_max_FWHM(path2, 'SampleB2_7_400C_24x24_t30_')
extract_min_FWHM(path2, 'SampleB2_7_400C_24x24_t30_')