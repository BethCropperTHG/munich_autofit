## autosignal.py
## Signal processing functions for automatic fitting of transfer reaction spectra
## Ben Cropper 2019

import scipy as sp
import numpy as np
from scipy.signal import find_peaks_cwt
from scipy.signal import gaussian
import matplotlib.pyplot as plt
import random

import adminfunctions as ad


def min_spectrum_split(x, y, fstgo, base_split = None):
    if base_split == None:
        while True:    
            try:        
                base_split = int(input('How many pieces should the spectrum be separated into?\n'))
                if base_split < 0: raise ValueError
                break
            except ValueError:
                print('Invalid number of sections entered')

    x_split_list = np.array_split(x,base_split)
    y_split_list = np.array_split(y,base_split)

    min_x_list = []

    for i, partx in enumerate(x_split_list):
        
        min_index = np.where(y_split_list[i] == min(y_split_list[i]))[0]        
        #sometimes there will be more than one channel with the minimum y. As to not bias this algortithm, I'll pick this randomly.        
        if len(min_index) > 1:
            if fstgo: print("The number of minima to choose from section ", i, " is:", len(min_index))
            min_index = random.choice(min_index)
        min_index = int(min_index)
        min_x = partx[min_index]
        min_x_list.append(min_x)

    min_indices = np.intersect1d(x, min_x_list, return_indices = True)[1]
    #print(min_indices)

    min_y_list = y[min_indices]    

    x_split_list = np.split(x, min_indices)
    y_split_list = np.split(y, min_indices)

    return(min_indices, min_x_list, min_y_list, x_split_list, y_split_list)

def min_region_finder(xlist, ylist, no_regions, fstgo):
    len_xlist = []
    for x in xlist:
        len_xlist.append(len(x))
    len_xlist = np.array(len_xlist)
    total_length = np.sum(len_xlist)

    fitting_regions = []

    for x, y in zip(xlist, ylist):
        no_individual_regions = np.around((len(x)/total_length) * no_regions)
        if no_individual_regions != 0:        
            split = min_spectrum_split(x,y, fstgo, base_split = no_individual_regions)
        else:
            split = [[],[],[],[x],[y]]        
        regions_list_x = split[3]
        regions_list_y = split[4]
        for x_region, y_region in zip(regions_list_x, regions_list_y):
            fitting_regions.append([x_region, y_region])
    return(fitting_regions)   


def contaminant_clipper(x, yinit, ybase, split = False):
    
    '''
    clips bits out of the spectrum after prompting the user to do so
    '''    
    yclipfit = yinit #clip things out of this yclip variable so yinit doesn't change
    yclipfind = ybase

    if split:
        xlist, ylist = [[],[]]
        borders = [0, len(x)]

    while True:
        
        #set a buffer such that it can be reinstated to the last one if the user changes their minddef contaminant_c_split(x, y)
        yclipfitbuffer = yclipfit
        yclipfindbuffer = yclipfind                

        #plot what the current spectrum is
        clip_spectrum = ad.spectrum_plotter(x, yclipfind)
        clip_spectrum.show()        
        
        #get the user to input a valid set of bounds in which to delete the contaminant
        while True:    
            lb = input('input a lower bound on the contaminant in channels')
            ub = input('input an upper bound on the contaminant in channels')
            if np.chararray.isdigit(ub + lb): #has to all be digits
                lb = int(lb)
                ub = int(ub)
                if (lb > min(x) and ub > lb and ub < max(x)): #lower and upper bounds have to be in the spectrum, upper bound has to be bigger than lower
                    break             
            print('Invalid bounds, try again')

        #get which indices x is less than the lower bound and greater than the upper bound set
        x21 = np.where(x < lb)[0]
        x22 = np.where(x > ub)[0]    
        
        #get y arrays above and below the contaminant
        y21fit = yclipfit[x21]
        y22fit = yclipfit[x22]

        y21find = yclipfind[x21]
        y22find = yclipfind[x22]

        #how many channels does the contaminant have
        length_difference  = len(x) - len(x21) - len(x22)

        #get an array of 0s to replace the contaminant ys with
        y20 = np.zeros(length_difference)

        #now to stick together the arrays
        #need the ys below the contaminant, then the zeroes where the contaminant was, and then the ys above the contaminant
        ynewfit = np.append(y21fit, y20)
        ynewfit = np.append(ynewfit, y22fit)

        ynewfind = np.append(y21find, y20)
        ynewfind = np.append(ynewfind, y22find)


        #now replace the old ys with the new one
        yclipfit = ynewfit
        yclipfind = ynewfind
        
        #now plot that and ask if they're happy with it
        clip_spectrum = ad.spectrum_plotter(x, yclipfind)
        clip_spectrum.show() 
        
        confirmclip = ad.y_or_n('Are you happy with this removal?')

        #what to do now? if they aren't happy, ask if they want another go. if they think there aren't any contaminants, return to yinit
        #if there still are, then they can have another go
        if confirmclip == 'n':

            cancel = ad.y_or_n('Are you sure there are any contaminants?')
            if cancel == 'n':
                yclipfit = yinit
                yclipfind = ybase
                if split: xlist, ylist = [[x],[yinit]]
                break
            if cancel == 'y':
                yclipfit = yclipfitbuffer
                yclipfind = yclipfindbuffer
                continue
            
        else: #all good, so can continue
            if split:
                borders.append(int(max(np.intersect1d(x, x[x21], return_indices = True)[1])))
                borders.append(int(min(np.intersect1d(x, x[x22], return_indices = True)[1])))
            else:            
                pass

        #do they want to get rid of more contaminants? continue loop if they do, don't if they don't
        moreclip = ad.y_or_n('Would you like to remove any more contaminants?')
        if moreclip == 'y':
            continue
        else:
            break

    if split:
        borders = np.sort(borders)
        for i, border in enumerate(borders):
            if ((i % 2 == 0) and (i is not len(borders))):
                ylist.append(yinit[border:borders[i+1]])
                xlist.append(x[border:borders[i+1]])
                

        return(xlist, ylist, yclipfind)

    return yclipfit, yclipfind

def peak_finder(peak_regions, w, fwhm):
    '''
    Tries to find the peak positions for each region, first by applying a low-pass filter
    then  a continuous wavelet transform, varying wavelet width over the plausible peak widths and finding ridges of local maxima in the cwt matrix.
    (see the scipy docs for a better description of how this works. It is very interesting!)
    '''
    regions_positions = []
    peak_positions = []
    for x2, y in peak_regions:
    
        #now do the cwt. the widths here are the widths of the wavelet you convolve the signal width
        #A peak which persists when you convolve the signal with wavelets of different widths will be considered a peak
        reg_ix = find_peaks_cwt(y, widths = np.arange( w/8 ,fwhm))
        try:        
            reg_positions = x2[reg_ix]
            #this one finds the absolute peak positions on the spectrum and lists them.
            peak_positions += reg_positions.tolist()
        except IndexError:
            reg_positions = [] 
      
        #this one sorts them into regions
        regions_positions.append(reg_positions)
    return peak_positions, regions_positions

