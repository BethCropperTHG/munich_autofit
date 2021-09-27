## autofit.py
## Automatic fitting routine for transfer reaction spectra
## Ben Cropper 2019
'''
####################################### Important Hyperparameter Setting #############################################
'''

'''The number of iterations of fitting that this will do before it works out which peaks will be in the final fit.
For example, if this number is set to 25, it will fit the spectrum 25 times with random selections of fitting region 
boundaries, and the peaks it finds consistently will be added to the final list for one final fit.'''
totalnfits = 25

'''The significance a peak requires to be counted as a peak. For example, if this number is set to 3,
any peak with 99 counts in it will need to have an uncertainty on this yield of 33 or less, or it will be
discarded and the region refitted.'''
significance = 2.5 

'''The consistency a peak requires to be counted as a peak in the final fit. A decimal between 0 and 1,
this number is the fraction of iterations this peak needs to be found during by the fitting. For example,
if this number is set to 0.5 and totalnfits is set to 20, and a peak is found at a particular position in
fewer than 10 of these iterations, it will not be counted as a true peak and discounted.''' 
consistency = 0.5

'''The significance between two peak positions required for peaks in different iterations in order to be classified
as the same peak. For example, if a peak has a a position of 900+-0 channels in one spectrum, the distance criterion
is set to 3, and the peak has a position of 903 in another iteration, the error on this position has to be greater
than 3 in order for them to be grouped together to be judged by the consistency criterion after the iterations are done.'''
distance = 4

'''This parameter is the distance the initial estimate of the position in the fitting is away from where the peak finding
finds the position to be. For example, if this is set to 5 and a peak is found at channel 900, the initial estimate
is set to 895 channels. This is to make sure some updating of the inverse Hessian for the position and height
components, making for better uncertainty estimation.'''
pos_offset = 8


'''
################################################################################################################
'''


import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
from scipy.signal import find_peaks_cwt, find_peaks
import glob
import os
import copy
import pickle

import adminfunctions as ad
import fittingfunctions as fit
import signalfunctions as sig

fstgo = True
#plt.rc('text',usetex = True)
'''
####################################### Read File and Show Spectrum #############################################
'''

directory = os.getcwd() + '/spectra'

fileno = ad.listfiles(directory)
fileselect = ad.inputfileselector(fileno)
f = ad.openfilefromlist(fileselect, directory)

x, yinit = ad.file_reader(directory + '/' + f)


spe = ad.spectrum_plotter(x, yinit)
spe.show()



'''
###########################################Contaminant Deletion##################################################
'''

cont_del =  ad.y_or_n("Would you like to remove any contaminants?")

if cont_del == 'y':
    xlist, ylist, yclip  = sig.contaminant_clipper(x,yinit,yinit, split = True)
else:
    xlist, ylist = [[x], [yinit]]
    yclip = yinit

for x2, y2 in zip(xlist, ylist):
    plt.plot(x2, y2)
plt.show()

 

'''
##########################################Threshold set and Peak Region Detection###########################################################
'''
rlist, betalist = [],[]

fitlistlist = []

for i in range(totalnfits):  
    while True:
        if fstgo:
            spe2 = ad.spectrum_plotter(x, yclip, 5)
            spe2.show()
    
            FWHM = float(input('What is the approximate FWHM of your peaks in channels'))
            width = FWHM/(2 * np.sqrt(2*np.log(2)))
    
            len_regions = float(input('Approximately how many channels long would you like your fitting segments to be?'))
            no_regions = len(x)/len_regions     
    
        peak_regions = sig.min_region_finder(xlist, ylist, no_regions, fstgo)
        i = 0
        none_indices = []
        for x2, y2 in peak_regions:
            if len(x2) == 0:
                none_indices.append(i)         
            if fstgo: plt.plot(x2, y2)
            i = i + 1
        if fstgo: plt.show()
        peak_regions = np.delete(peak_regions, none_indices, axis = 0)
    
    
        if fstgo:
            recontaminant = ad.y_or_n('Would you like to re-clip for any more contaminants?')
        
            if recontaminant == 'y':
                yclipfit, yclipfind  = sig.contaminant_clipper(x,yclipfit, yclipfind)
                continue
            else:
                pass
    
            final_regions = ad.y_or_n('Would you like to re-tune the width and number of peaks?')
    
            if final_regions == 'y':
                continue
            else:
                break
        else: break
    




    '''
    ###########################################Fitting##################################
    '''
    if fstgo:
        bg = ad.y_or_n('Would you like to have a local linear background added to these fits?')

        if bg == 'y':
            bg = True
        else:
            bg = False

    peak_positions, region_positions = sig.peak_finder(peak_regions, width, FWHM)

    pos = np.zeros(len(x))
    pos[np.intersect1d(x, peak_positions, return_indices = True)[1]] = 1.2 * max(yclip)

    if fstgo:
        fig = plt.figure(figsize = (15,7))
        axis = fig.add_subplot(111)
        axis.set_xticks(np.arange(0,max(x), 100))
    
        axis.plot(x,yinit, color = 'xkcd:light grey')
        for x2,y in peak_regions:
            axis.plot(x2,y)
    
        axis.plot(x,pos, linewidth = 1)

        plt.show()

        template, template_pos, peak_pos = ad.template_finder(peak_regions, peak_positions, region_positions, fstgo, peaks_figure = fig)
    else:
        template, template_pos, peak_pos = ad.template_finder(peak_regions, peak_positions, region_positions, fstgo, pos = peak_pos)



    template_fit = fit.fit(template[0], template[1], template_pos, 0.1, FWHM, rbfix = False, background = bg, fig = False)
    ad.printfit(template_fit, template[0], template[1])

    rfix, betafix = template_fit[1][-2],template_fit[1][-1]
    rlist.append(rfix)
    betalist.append(betafix)

    fitlist = []
    for i, region in enumerate(peak_regions):
        xreg = region[0]
        yreg = region[1]
        if len(region_positions[i]) == 0: continue     
        ft = fit.fit(xreg, yreg, region_positions[i], 0.1, FWHM, r = rfix, beta = betafix, background = bg, rbfix = False, fig = False)
        #ft = fit.fit(xreg, yreg, region_positions[i], width, FWHM, rbfix = False, background = bg, fig = True)
        ad.printfit(ft, xreg, yreg)

        fitlist.append(ft)

    '''
    ############################################2nd pass#############################
    '''
    muarrnew = []
    aarrnew = []
    while True:
        oldlen = len(aarrnew)    
        muarrnew = []
        aarrnew = []

        for ft in fitlist:
            for i, yiel in enumerate(ft[3]):
                if not np.isnan(ft[4][i]):
                    #if ft[9]:
                    if yiel > significance * ft[4][i]:
                        muarrnew.append(ft[1][2 * i + 1])
                        aarrnew.append(ft[1][2 * i])
                    

    
        newlen = len(aarrnew)
    
        if oldlen == newlen: break
        fitlist = []    
        muarrnew = np.array(muarrnew)
        aarrnew = np.array(aarrnew)
    
        for i, region in enumerate(peak_regions):
            xreg = region[0]
            yreg = region[1]
            


            reg_pos = muarrnew[(muarrnew > min(xreg)) & (muarrnew < max(xreg))]
            a_arr = aarrnew[(muarrnew > min(xreg)) & (muarrnew < max(xreg))]
    

            if len(region_positions[i]) == 0: continue     
            #ft = fit.fit(xreg, yreg, reg_pos, width, FWHM, r = rfix, beta = betafix, background = bg, Aarr = a_arr, fig = False, rbfix = True)
            ft = fit.fit(xreg, yreg, reg_pos, 0.1, FWHM, rbfix = False, background = bg, Aarr = a_arr, fig = False)
            ad.printfit(ft, xreg, yreg)

            fitlist.append(ft)
     
    fitlistlist.append(fitlist)    
    fstgo = False

#now want to plot pos vs yield for every peak in this list.
#print(fitlistlist)

afin = []
mfin = []
smfin = []
yfin = []
nlist = []

for j, fitlist in enumerate(fitlistlist):
    for ft in fitlist:
        for i, yiel in enumerate(ft[3]):
            mfin.append(ft[1][2 * i + 1])
            afin.append(ft[1][2 * i])
            smfin.append(np.sqrt(ft[2][2*i + 1][2*i + 1]))
            yfin.append(yiel)
            nlist.append(j)

plt.errorbar(nlist, mfin, smfin, linestyle = "")
plt.show(block = False)


poslist = []
sposlist = []
alist = []
for mu, smu, a,y in zip(mfin, smfin, afin, yfin):
    z = 10000
    for pos, spos, h in zip(poslist,sposlist, alist):
        p,sp = np.average(pos, weights = 1/np.array(spos)**2, returned = True)
        sp = sp ** (-0.5)
        s = np.sqrt(sp**2 + smu**2)
        z = np.abs(mu - p)/s
        if np.isnan(p): z = 10000
 
        if z < distance:
            pos.append(mu)
            spos.append(smu)
            h.append(a)
            break
        '''
        if y > 1000:
            pos.append(mu)
            spos.append(smu)
            h.append(a)
            break
        '''
    if z < distance: continue
    
    poslist.append([mu])
    sposlist.append([smu])
    alist.append([a])

for i, pos in enumerate(poslist):
    plt.errorbar(np.full(len(pos), i), pos, sposlist[i], linestyle = "" )

plt.show(block = False)



muarrnew = []
aarrnew = []
#print(alist)
for pos, a in zip(poslist, alist):
    if len(pos) > consistency * totalnfits:
        muarrnew.append(np.average(pos))
        aarrnew.append(np.average(a))

muarrnew = np.array(muarrnew)
aarrnew = np.array(aarrnew) 

#print(muarrnew, aarrnew)
peak_regions = sig.min_region_finder(xlist, ylist, no_regions, fstgo)
i = 0
none_indices = []


for x2, y2 in peak_regions:
    if len(x2) == 0:
        none_indices.append(i)         
        i = i + 1
    peak_regions = np.delete(peak_regions, none_indices, axis = 0)



reg_pos = []
a_arr = []

for i, region in enumerate(peak_regions):
    xreg = np.array(region[0])
    yreg = np.array(region[1])
       
    this_reg_pos = muarrnew[(muarrnew > min(xreg)) & (muarrnew < max(xreg))]
    this_a_arr = aarrnew[(muarrnew > min(xreg)) & (muarrnew < max(xreg))]

    pos_order = np.argsort(this_reg_pos)
    this_reg_pos = this_reg_pos[pos_order]
    this_a_arr = this_a_arr[pos_order]        
    
    reg_pos.append(this_reg_pos)
    a_arr.append(this_a_arr)
  

rfix, betafix = np.average(rlist), np.average(betalist)


fitlist = []


for i, region in enumerate(peak_regions):
    xreg = np.array(region[0])
    yreg = np.array(region[1])

    if len(reg_pos[i]) == 0: continue
    xreg = xreg[xreg < max(reg_pos[i]) + 2 * FWHM]
    yreg = yreg[:len(xreg)]

    reg_pos[i][reg_pos[i] - pos_offset < min(xreg)] = min(xreg) + pos_offset     
    rbfind = fit.fit(xreg, yreg, np.round(reg_pos[i]) - pos_offset + 0.5, 0.1, FWHM, rbfix = False, background = bg, posfix = False, r = 50, beta = FWHM/2, fig = False)
    rfix, betafix = rbfind[1][-2],rbfind[1][-1]
    ft = fit.fit(xreg, yreg, np.round(reg_pos[i]) - pos_offset + 0.5, 0.1, FWHM, r = rfix, beta = betafix, background = bg, posfix = False, rbfix = True)

    for err in ft[4]:
        if np.isnan(err):
            print("No peak found with this initial position. Trying again with initialisation closer to peak")
            ft = fit.fit(xreg, yreg, np.round(reg_pos[i]) - 2.5, 0.1, FWHM, r = rfix, beta = betafix, background = bg, posfix = False, rbfix = True)

    #ad.printfit(rbfind, xreg, yreg)
    ad.printfit(ft, xreg, yreg)



    fitlist.append(ft)



'''
##########################################Plot Fits and Save######################
'''

ysub = copy.deepcopy(yclip)

fitplot = plt.figure(figsize = (15,7))
a = fitplot.add_subplot(111)
a.plot(x, yinit, color = 'xkcd:grey')

ofilename = str(f[:-4]) + '_fit.txt'
fil = open(ofilename, 'w+')
topline = 'POSITION sPOSITION AREA sAREA WIDTH sWIDTH R sR BETA sBETA\n'
fil.write(topline)

yieldlist = []
yerrlist = []
mufitlist = []

lines = []

for fit in fitlist:
    #plot whole fit    
    a.plot(fit[7], fit[0], linewidth = 3)
    
    #plot background
    if fit[2][-4][-4] != 0:
        plt.plot(fit[7], np.full(len(fit[7]), fit[1][-4]))

    for y_individual in fit[6]:
        a.plot(fit[7], y_individual)

    #plot residuals
    j = 0
    ix = np.intersect1d(x, fit[7], return_indices = True)[1]
    for i in ix:
        ysub[i] = ysub[i] - fit[0][j]
        j = j + 1
    #save peaks to file
    for i, peak in enumerate(fit[3]):
        line = (fit[1][2 * i + 1], np.sqrt(fit[2][2 * i + 1][2 * i + 1]), peak, fit[4][i], fit[1][-3], np.sqrt(fit[2][-3][-3]), fit[1][-2],np.sqrt(fit[2][-2][-2]), fit[1][-1],np.sqrt(fit[2][-1][-1]))
        lines.append(line)
        
        yieldlist.append(peak)
        yerrlist.append(fit[4][i])
        mufitlist.append(fit[1][2 * i + 1])

linetype = [('POSITION',float), ('sPOSITION',float), ('AREA',float), ('sAREA',float), ('WIDTH',float), ('sWIDTH',float), ('R',float), ('sR',float), ('BETA',float), ('sBETA',float)]
lines = np.array(lines,dtype = linetype)
lines = np.sort(lines, order = 'POSITION')#[::-1]

for line in lines:
    linestring = ''
    for element in line:
        linestring += str(element) + ' '
    linestring = linestring[:-1] + '\n'
    fil.write(linestring)

fil.close() 
try:
    yieldlist = np.array(yieldlist)
    yerrlist = np.array(yerrlist)
    mufitlist = np.array(mufitlist)

    yielderrorplot = plt.figure(figsize = (16,9))
    yieldax = yielderrorplot.add_subplot(221)
    yieldax2 = yielderrorplot.add_subplot(222)
    yieldax3 = yielderrorplot.add_subplot(223)
    yieldax4 = yielderrorplot.add_subplot(224)

    yieldax.errorbar(mufitlist, yieldlist, yerrlist, fmt = '.')
    yieldax.set_xlabel('Position', size = 'xx-large')
    yieldax.set_ylabel('Yield', size = 'xx-large')

    yieldax2.plot(mufitlist, yerrlist, 'o', color = 'r', label = 'Error')
    yieldax2.plot(mufitlist, np.sqrt(yieldlist) , 'o', color = 'b', label = 'Scale (sqrt(yields))')
    yieldax2.set_xlabel('Position', size = 'xx-large')
    yieldax2.set_ylabel('Yield Error', size = 'xx-large')
    yieldax2.legend()

    yieldax3.plot(mufitlist, yerrlist/yieldlist, 'o', color = 'g')
    yieldax3.set_xlabel('Position', size = 'xx-large')
    yieldax3.set_ylabel('Yield Fractional Error', size = 'xx-large')
    
    yieldax4.plot(mufitlist,yerrlist/np.sqrt(yieldlist), 'o', color = 'xkcd:purple')
#    yieldax4.set_xlabel(r'$\sigma/\sqrt{N}$', size = 'xx-large')
    yieldax4.set_ylabel('sigma/root(N)', size = 'xx-large')
    yieldax4.set_xlabel('Position', size = 'xx-large')
    yieldax4.plot(np.arange(len(x)),np.full(len(x),2))

except:
    print('Error in plotting errors in the yields.')

    print('\nThe array of yields is:')
    print(yieldlist)
    print('The length of this array is: ',len(yieldlist))  

    print('\nThe array of yield errors is:')
    print(yerrlist)
    print('The length of this array is: ',len(yerrlist)) 

    print('\nThe array of peak positions is:')
    print(mufitlist)
    print('The length of this array is: ',len(yieldlist)) 
    plt.plot(np.arange(len(x)),np.full(len(x),2))
    
pos = np.zeros(len(x))
pos_ix = []
for mu in muarrnew:
    pos_ix.append(ad.find_nearest(mu, x))
pos[pos_ix] = max(yclip) * 1.2

a.plot(x, pos, linewidth = 1)

a.plot(x, ysub)
f2 = str(f[:-4]) + '_fig.pkl'
print(os.getcwd(), f2)
pickle.dump(fitplot, open(f2, 'wb'))
plt.show()

print(FWHM, len_regions, bg)
if cont_del == 'y':
    for slice in xlist:
        print(min(slice), max(slice)) 
    
