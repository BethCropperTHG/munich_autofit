## adminfunctions.py
## Functions for admin purposes (reading files, getting yes or no answers etc) for automatic fitting of transfer reaction spectra
## Ben Cropper 2019

import os
import numpy as np
import matplotlib.pyplot as plt

def file_reader(handle, delete_0 = True):
    '''
    Reads a file given a handle. The types of files that it reads are 2 column data - histograms with slightly weird binning so you can't just use 1 column and then 1,2,3,4, etc for the bin data.
    I've done it so it can read things that are comma separated and space separated. Again it's a bit weird because it's specific to things created with a 'dump' function that I have for converting
    ROOT histograms to comma separated text files, and a function in RadWare (https://radware.phy.ornl.gov/) (https://web.archive.org/web/20180612114401/https://radware.phy.ornl.gov/) which converts from 
    .spe files used in RadWare to ascii .asc text files. Thus, this function might need to be tweaked if not being used for these specific formats.
    '''
    
    print('reading ', handle)
    
    infile = open(handle,'r')
    xvals=[];ymeas=[]
    
    #need this conditional because the asciis are comma separated and txts are space separated 
    if handle[-1] == 'c': #for .asc extension
        #this bit reads in the files
        while True:
            line = infile.readline()
            if not line: break
        
            items = line.split()
            #one column is x, one is y
            xvals.append(float(items[0][:-1]))#the -1 to get rid of the comma at the end
            ymeas.append(float(items[1]))

    else: #for 2 column space separated text files
        while True:
            line = infile.readline()
            if not line: break
        
            items = line.split()
            xvals.append(float(items[0]))
            ymeas.append(float(items[1]))

    infile.close()

    #x and y should be numpy arrays so we can do more stuff with them
    ymeas = np.array(ymeas)
    xvals = np.array(xvals)

    if delete_0:
        #gets rid of the 0 channel dump. Usually in these things the first channel contains all of the dead time counts, which are to be ignored.
        ymeas = np.delete(ymeas, [0])
        xvals = np.delete(xvals, [0])

    print('successfully read ', handle)

    return xvals, ymeas

def y_or_n(question):
    '''
    simple function for yes/no inputs. resistant to most typos, except if the user doesn't input anything. A simple spell, but not really that unbreakable. Convenient though.
    '''
    while True: #ininite loop that breaks whenever the user inputs either 'y' or 'n'
        yorn = input(question + 'y/n')
        if yorn == 'y': break
        elif yorn == 'n' : break 
        else: print('invalid input, please try again')
        
    return yorn

def listfiles(directory):
    #set an index, since it's going to print out the index before each file so the user can
    #select which file to access more easily
    fileindex = 0
    #simple loop through the directory
    for roots, dirs, filenames in os.walk(directory):
        for f in filenames:
            fileindex = fileindex + 1
            #need to print the name of each file out so the user can select which one to open 
            print( str(fileindex) + ': ' + f)
    #returns the total number of files so the program can detect whether the
    #selected file number is valid (ie falls in between 0 and the number of files)
    return fileindex

def inputfileselector(fileno): #fileno is the number of files in the directory
    #Initialize the desired quantity
    fileselect = None
    while fileselect == None:
        #prompt the user to input their selection
        fileselect = input('Which file would you like to open?\n')
        try:
            
            #failure test
            #in this case, it checks to see if the input matches a file index
            
            #create a list of the possible file indices
            fileindices = list(range(1, fileno + 1))
            #filetest doesn't do annything, but if there's an error in the assignment, it
            #triggers the error
            filetest = fileindices.index(int(fileselect))

        except:
            fileselect = None
    return int(fileselect)

def openfilefromlist(fileselect, directory):
    '''
    given an integer fileselect, it will return the filename of the fileselect-th file in the directory. Combine with inputfileselector.
    '''
    for roots, dirs, filenames in os.walk(directory):
        f = filenames[fileselect-1]
    return f

def spectrum_plotter(x, y, ticks = 100):
    '''
    Plots a spectrum of the right size and with the right axes, so you don't end up with tiny spectra or channel ticks of '0, 500, 1000, 1500, 2000' or something equally as useless.
    '''
    fig = plt.figure(figsize = (15,7)) # short and wide fig size
    ax = fig.add_subplot(111)
    ax.plot(x, y)
    ax.plot(x, np.zeros(len(x)))
    #ax.set_xticks(np.arange(0,max(x), ticks)) #Labels every 100 channels

    return fig

def template_finder(peak_regions, peak_positions, reg_pos, fstgo, peaks_figure = None, pos = None, A_arr = None):
    '''
    Function for getting the user to select a template region from a set of peak regions. specifically it will ask for the best peak, and then it will fit the region with that best peak in it
    '''
    if fstgo:
        #print out all of the peak positions, assigning them a number
        for i,p in enumerate(peak_positions):
            print( '[' + str(i) + ']: ' + str(p))

    #define this - will be the number assigned to the peak which will be selected as the template
    #or none if it is not a valid number
    template_select = None

    while template_select == None:
    
        try: #this try block will try and cast the user's input to an int, will keep pestering the user if they don't input one
            if(fstgo):
                peaks_figure.show()            
                template_select = int(input('Which peak is the most suitable as a template to fix parameters to The region containing it will be what is fitted.'))
                i = 0
                for x2, y in peak_regions:
                    if len(np.intersect1d(peak_positions[template_select], x2)) == 1: #if the peak position is in the peak region, it selects that region
                        template = [x2, y]
                        template_pos = reg_pos[i] #returns which region the teplate is
                        if A_arr is not None:
                            template_a = A_arr[i]

                        peak_pos = np.intersect1d(peak_positions[template_select], x2)
                    i = i + 1
            else:            
                i = 0
                for x2, y in peak_regions:
                    print(pos, min(x2),min(y))
                    if len(np.intersect1d(pos, x2)) == 1: #if the peak position is in the peak region, it selects that region
                        template = [x2, y]
                        template_pos = reg_pos[i] #returns which region the teplate is
                        if A_arr is not None:
                            template_a = A_arr[i]
                        peak_pos = np.intersect1d(pos, x2)
                    i = i + 1
                template_select = 'fish'
        except:
            template_select = None #resets the loop


    if A_arr is None: return template, template_pos, peak_pos
    else: return template, template_pos, peak_pos, template_a

def printfit(fit, x, y, endstop = True):
    '''
    Prints out the fit parameters for each fit similarly to gf3(RadWare), and plots out the fit, along with the constituent peaks and original data, on one figure.
    '''

    #unpacks all the needed variables - in this case the fit (ymod), the fit parameters (p1) and covariance matrix (p1cov), the yields (yieldarr) and uncertainties (yerrarr)
    ymod, p1, p1cov, yieldarr, yerrarr = fit[0], fit[1], fit[2], fit[3], fit[4]
    #these parameters are the same for all of the peaks so should be printed separately.
    off, sig, r, beta = p1[-4], p1[-3], p1[-2], p1[-1]

    #print(p1)
    
    #need the number of peaks, I want to loop over all of the peaks to print their parameters.
    #-3 is to get rid of the metaparameters, divide by 2 because it stores the amplitude and position of each peak     
    nopeaks = int((len(p1) - 4)/2)
    #calculate chi-squared. poisson errors so sigma^2 = y
    chi2 = np.sum((y - ymod)**2/(y+1))
    #if parameters are fixed, the number of degrees of freedom changes: dof = number of data points - number of free parameters     
    if p1cov[-1][-1] == 0:
        ndof = len(y) - len(p1) + 2
    else:
        ndof = len(y) - len(p1)
    #Print the general parameters of the fit - the range of the fit, the number of channels, and the shape parameters that stay constant across the fit
    print('Fit complete. Fitted channels ', min(x), ' to ', max(x), ', with ', nopeaks, ' peak(s).') 
    print('Width = ', paramprint(sig, np.sqrt(p1cov[-3][-3])), '\nR = ', paramprint(r, np.sqrt(p1cov[-2][-2])), '\nBeta = ', paramprint(beta, np.sqrt(p1cov[-1][-1])), '\nBackground Offset', paramprint(off, np.sqrt(p1cov[-4][-4])) , '\nReduced chi-squared = ', np.round(chi2/ndof, decimals = 2)) 

    #now to print the positions and yields
    print('    Positions           Yields              Heights')
    #loop over the number of peaks
    for i in range(nopeaks):
        #get the right thing from each array of relevant parameters
        Y = yieldarr[i]
        sY = yerrarr[i]
        pos = p1[2*i + 1]
        spos = np.sqrt(p1cov[2*i + 1][2*i + 1]) #errors come from the square root of the diagonals of the covariance matrix
        height = p1[2*i]
        sheight = np.sqrt(p1cov[2*i][2*i])

        #need to line up the columns - will make each one 30 characters long
        posspaceslen = 20 - len(paramprint(pos, spos))
        yspaceslen = 20 - len(paramprint(Y, sY))
        posspaces, yspaces = '',''

        for j in range(posspaceslen):
            posspaces += ' '

        for j in range(yspaceslen):
            yspaces += ' '

        print('[' + str(i+1) + '] ' + paramprint(pos, spos) + posspaces + paramprint(Y, sY) + yspaces + paramprint(height, sheight))
    #show the plot generated by the fitting function of the fit and constituent peaks
    if fit[5] is not None:
        fit[5].show()
        #slightly clumsy way of not just printing everything at once - making the user input something to move it on. Annoyingly, pyplot's blocking behaviour doesn't apply in the object-oriented interface
        if endstop: input('Press enter to continue')
        plt.close(fit[5])    

    
def paramprint(param, error):
    '''
    This prints out a parameter with its error, rounding both to the preision of the uncertainty
    '''
    try:
        errormag =  int(np.floor(np.log10(error))) #sets the precision of the uncertainty on the error. Rounds down on the order of magnitude
    except: #for errors of zero    
        errormag = 0

    if np.isnan(error):
        printout = str(param) + ' +- ' + str(error)
        return(printout)

    if errormag > -1: #round to the error if there's a decimal point
        param = int(np.around(param))
        error = int(np.around(error))
    else: #round to the nearest integer if the error is greater than 1
        param = np.around(param, -errormag)
        error = np.around(error, -errormag)

    
    printout = str(param) + ' +- ' + str(error)
    return(printout)

def find_nearest(val, arr):
    'find element of numpy array arr which is closest to value val'
    return (np.abs(arr - val)).argmin()

