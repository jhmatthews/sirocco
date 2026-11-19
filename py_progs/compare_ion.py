#!/usr/bin/env python 

'''
This routine is used for diagnosing isses with spectra produced during ionization cycles

Command line usage (if any):

    usage: Compare.py root1 root2

Description:  

    This routine compares two spectral files created in
    ionization cycles, specifically the log_spec_tot files

Primary routines:

    compare_ion

Notes:

    As run from the command line, two WCreated spectra
    are compared.  

    One has somewhat more flexibility if run from a 
    Jupyter notebook
                                       
History:

211003 ksl Coding begun

'''

import sys
from astropy.io import ascii
import os
from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np

from scipy.signal.windows import boxcar
from scipy.signal import convolve

def xsmooth(flux,smooth=21):
    '''
    boxcar smooth the flux
    '''
    if (smooth)>1:
        q=convolve(flux,boxcar(smooth)/float(smooth),mode='same')
        return(q)
    else:
        return(flux)



def edge(ymin=1e45,ymax=2e45):
    speed_of_light=3e18
    nu=3e18/912
    plt.plot([nu,nu],[ymin,ymax],'k')
    plt.text(nu,ymax,r'Ly',ha='center',va='bottom')


def compare_ion(one='one/xmatrix_one',two='two/yjump',col='WCreated',outfile=''):
    
    label1=one
    label2=two

    if outfile=='':
        word1=one.split('/')
        word2=two.split('/')
        outfile='diff_%s_%s.png' % (word1[len(word1)-1],word2[len(word2)-1])


    
    try:
        xone=ascii.read(one+'.log_spec_tot')
    except:
        print('Could not find ',one)
        return
    
    try:
        xtwo=ascii.read(two+'.log_spec_tot')
    except:
        print('Could not find ',two)
        return

    plt.semilogx(xone['Freq.'],xsmooth(xone['Freq.']*xone[col]),label=label1)
    plt.semilogx(xtwo['Freq.'],xsmooth(xtwo['Freq.']*xtwo[col]),label=label2)
    plt.semilogx(xone['Freq.'],(xsmooth(xone['Freq.']*xone[col])-xsmooth(xtwo['Freq.']*xtwo[col])),label='%s-%s' % (label1,label2))
# edge()
    plt.xlim(1e15,1e17)
# plt.ylim(1e39,1e46)
    plt.legend()
    plt.savefig(outfile)    

    


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>2:
        compare_ion(sys.argv[1],sys.argv[2])
    else:
        print ('usage: Compare.py root1 root2')

