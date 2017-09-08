
from peak_detection_matlab import peakdet
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from os.path import basename
from scipy.optimize import curve_fit
import pandas as pd
import scipy
from scipy import signal
#%%
'''


This is the program that Fang and JrHS use to do our fitting

It has been updated on the back end to automatically extract the maximum and 
minimum FWHM from the data set and push them into separate files.

Fang has written a short script that cycles through criterions.  I wrote 
a sub-routine that uses the Average Intensity of the pattern to determine
the appropriate starting criterion and increment.



'''

def func(x, *params):
    """
    create a Lorentzian fitted curve according to params
    """
    y = np.zeros_like(x)
    for i in range(0, len(params), 4):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        n = params[i+3]
        n = np.clip(n, 0,1) 
        y = y + n * amp * np.exp(-4 * np.log(2) * ((x - ctr) / wid) ** 2) + (1 - n) * amp * wid ** 2 / 4 / (
        (x - ctr) ** 2 + wid ** 2 / 4)
    return y

#path = 'C:\\Users\\jrh8\\Documents\\Documents\\Projects\\Metallic Glasses\\Experiment B2 Data\\May2017\\SampleB2_24_200C\\Processed\\subtracted'

def XRDfitting(path):

    save_path = path + '\\peak_fitting_GLS\\'
    
    if not os.path.exists(save_path):
            os.mkdir(save_path)
    
    
    # parameters for 8, 9, 10, 17, 18
    # guess = [2.05, 25, 0.233, 0.5, 3.1, 5000, 0.33, 0.5, 3.5, 23, 0.15, 0.5]
    # high = [2.1, 30, 0.3, 1, 3.4, 10000, 1, 1, 3.55, 30, 0.2, 1]
    # low = [2.0, 20, 0.1, 0, 2.1, 0, 0, 0, 3.45, 10, 0.1, 0]
    
    # parameters for 14
    # guess = [3.9, 100, 0.3, 0.5, 3.1, 5000, 0.33, 0.5]
    # high = [4.0, 175, 0.5, 1, 3.4, 10000, 1, 1]
    # low = [3.8, 75, 0.1, 0, 2.1, 0, 0, 0]
    criterion_list = pd.DataFrame()
    
    for filename in glob.glob(os.path.join(path, '*.csv')):
        if basename(filename)[-5] == 'd' or 1:
            print 'processing', filename
            data = np.genfromtxt(filename, delimiter = ',', skip_header = 1)
            Qlist = data[:,0][:647]
            IntAve = data[:,1][:647]
            #IntAve = scipy.signal.savgol_filter(IntAve, 31, 2, mode ='mirror')
            Average_intensity = np.mean(IntAve)
            max_intensity = np.max(IntAve)
            criterion = Average_intensity /5
            criterion_increment = criterion * 3
            ''' While fitting our data we noticed that samples too close to the
            background were failing during Peak Fitting. This script uses the
            average intensity to remove samples with too low intensity. It 
            does not solve all ValueError Problems'''
            if abs(Average_intensity) < 7:
                    print "Why won't this raise?"
                    nofit = np.matrix('0,0,0,0')
                    np.savetxt(save_path + basename(filename)[:-4] + '_peak_analysis_GLS.csv',nofit, delimiter =",")
            else:
                while 1:
                    try:
                        maxs, mins = peakdet(IntAve, criterion)
                        peaks = maxs[:,0].astype(int)
                        guess = []
                        low = []
                        high = []
                        print peaks
                        for peak in peaks:
                            guess.append(Qlist[peak])
                            low.append(Qlist[peak]-Qlist[peak]*.15)
                            high.append(Qlist[peak]+Qlist[peak]*.15)
        
                            guess.append(IntAve[peak])
                            low.append(0)
                            high.append(IntAve[peak]+10)
                            '''the first append is to append the guess and bounds 
                            for width, the second append is to add guess and bounds 
                            for the n (different components for peak shape)
                            '''
                            guess.append(0.2)
                            low.append(0)
                            high.append(1.2)
        
                            guess.append(0.5)
                            low.append(0)
                            high.append(1.2)
                        
                        popt, pcov = curve_fit(func, Qlist, IntAve, p0=guess, bounds = (low, high))
                        #popt, pcov = curve_fit(func, Qlist, IntAve, p0=guess)
                        fit = func(Qlist, *popt)
                        plt.figure(1)
                        plt.plot(Qlist, IntAve)
                        plt.plot( Qlist, fit, 'r--')
                        plt.plot(Qlist[peaks], IntAve[peaks], 'o')
        
                        for i in range(0, len(popt), 4):
                            ctr1 = popt[i]
                            amp1 = popt[i+1]
                            wid1 = popt[i+2]
                            n1 = popt[i+3]
                            curve1 = n1 * amp1 * np.exp( -4 * np.log(2) * ((Qlist - ctr1)/wid1)**2) + (1-n1) * amp1 * wid1**2 / 4 / ((Qlist-ctr1)**2 + wid1**2 / 4)
                            plt.plot(Qlist, curve1)
        
                        plt.savefig(save_path + basename(filename)[:-4] + '_peak_analysis_GLS')
                        plt.close()
        
                        popt = np.reshape(popt, (popt.size/4, 4))
                        np.savetxt(save_path + basename(filename)[:-4] + '_peak_analysis_GLS.csv', popt, delimiter=",")
                        break
                    except RuntimeError:
                        print "Failed to fit", filename
                        criterion = criterion + criterion_increment
                        print "Increase criterion by " + str(criterion_increment) +" . The criterion is now", criterion
                        
    
    
                    np.savetxt(save_path + basename(filename)[:-4] + '_peak_analysis_GLS.csv', popt, delimiter=",")
                 
                

'''
The next part of this program scrapes the created fit files for the maximum FWHM and places them in an easily
parsable file for subsequent harvesting.


sample_name = "SampleB2_24_200C_24x24_t30_"
'''

#%%
def extract_max_FWHM(path, sample_name):
    path = path + 'peak_fitting_GLS'
    sample_name2 = "_bckgrd_subtracted_peak_analysis_GLS"
    os.chdir( path )
    sample_n = 441
    x0 = 0 
    filename2 = sample_name + '0001' + sample_name2+'.csv'
    draft = pd.read_csv(filename2, sep = ',', header = -1)
    '''draft = draft[["'Id:", "Unnamed: 1"]]'''
    
    draft7 = ''
    
    
    for i in range(sample_n):
        x = x0 + i+1
        def file_index(x):
            if len(str(x)) == 1:
                return '000' + str(x)
            elif len(str(x)) == 2:
                return '00' + str(x)
            elif len(str(x)) == 3:
                return '0' + str(x)
            elif len(str(x)) == 4:
                return str(x)
        sample_name1 = sample_name + file_index(x) + sample_name2 + ".csv"
        print(sample_name1)
        draft1 = pd.read_csv(sample_name1, sep = ',', header = -1)
        print(draft1)
        '''draft1 = draft1[["Unnamed: 1"]]'''
        draft2 = draft1[0] >1.7
        draft3 = draft1[0] <4
        draft4 = draft1[draft2&draft3]
        draft5 = draft4[2]
        draft6 = draft5.max(0)
        draft = pd.concat([draft,draft1], axis = 0)
        draft7 = draft7 + '\n' + str(draft6)
        print(draft)
        combiviewmax = 'Combiview' + sample_name + '.txt'
        draft.to_csv(combiviewmax, sep = '\t')
        maxfwhmfile = 'Maximum FWHM' + sample_name + '.txt'
        np.savetxt(maxfwhmfile, [draft7], fmt='%s')
      
'''
The next part of this program scrapes the created fit files for the minimum FWHM and puts them into an easily parasable file 
for subsequent harvesting 

''' 
#%%
def extract_min_FWHM(path, sample_name):
    path = path +'peak_fitting_GLS'
    sample_name2 = "_bckgrd_subtracted_peak_analysis_GLS"
    os.chdir( path )
    sample_n = 441
    x0 = 0 
    filename2 = sample_name + '0001' + sample_name2+'.csv'
    draft = pd.read_csv(filename2, sep = ',', header = -1)
    '''draft = draft[["'Id:", "Unnamed: 1"]]'''
    
    draft7 = ''
    
    
    for i in range(sample_n):
        x = x0 + i+1
        def file_index(x):
            if len(str(x)) == 1:
                return '000' + str(x)
            elif len(str(x)) == 2:
                return '00' + str(x)
            elif len(str(x)) == 3:
                return '0' + str(x)
            elif len(str(x)) == 4:
                return str(x)
        sample_name1 = sample_name + file_index(x) + sample_name2 + ".csv"
        print(sample_name1)
        draft1 = pd.read_csv(sample_name1, sep = ',', header = -1)
        print(draft1)
        '''draft1 = draft1[["Unnamed: 1"]]'''
        draft2 = draft1[0] >1.7
        draft3 = draft1[0] <4
        draft4 = draft1[draft2&draft3]
        draft5 = draft4[2]
        draft6 = draft5.min(0)
        draft = pd.concat([draft,draft1], axis = 0)
        draft7 = draft7 + '\n' + str(draft6)
        print(draft)
        minfwhmfile = 'Minimum FWHM' + sample_name + '.txt'
        np.savetxt(minfwhmfile, [draft7], fmt='%s')
  
   

