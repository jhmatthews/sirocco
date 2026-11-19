#!/usr/bin/env python 

'''
This is a general purpose routine for plotting variables  written out by windsave2table for 1d models


Command line usage (if any):

    usage::

        plot_wind_1d.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

191102 ksl
    Coding begun.  This routine is largely parallel to the routine
    plot_wind.py (which is intended for 2d models

'''

import sys
from astropy.io import ascii
import numpy
import matplotlib.pyplot as pylab


def get_data(filename='fiducial_agn_master.txt', var='t_r',grid='ij',inwind='',scale='guess',zmin=-1e50,zmax=1e50):
    '''
    This routine reads and scales the data from a single variable that is read from an ascii table 
    representation of one or more of the parameters in a windsave file

    '''

    try:
        data=ascii.read(filename)
        for one in data.colnames:
            if one=='j':
                print('Error: This appears to be a file for a 2d wind')
                raise ValueError
    except IOError:
        print('Error: file %s does not appear to exist' % filename)
        raise ValueError

    # data.info()





    title=var


    ii=data['i']


    iimax=numpy.max(ii)
    if grid.count('i'):
        x=numpy.array(data['i'])
        xlabel='i'
    elif grid=='log':
        x=numpy.array(data['x'])
        xmin=1e50
        while i<len(x):
            if x[i]>1 and x[i]<xmin:
                xmin= x[i]
            i+=1
        xlogmin=numpy.log10(xmin/10)
        x=numpy.select([x>1],[numpy.log10(x)],default=xlogmin)
        xlabel='log(x)'
    else:
        x=numpy.array(data['x'])
        xlabel='x'

    xvar=numpy.array(data[var])
    xwind=numpy.array(data['inwind'])
    # now mask out the appropriate data

    if inwind=='all':
        xwind.fill(0) # Keep everything
    elif inwind.count('partial'):
        xwind=numpy.select([xwind<0],[xwind],default=0)

    # Now make the mask regardless
    mask=numpy.ma.make_mask(xwind)


    xvar=numpy.ma.array(xvar,mask=mask)


    vmin=numpy.ma.min(xvar)
    vmax=numpy.ma.max(xvar)



    if vmin<zmin:
        vmin=zmin
    if vmax>zmax:
        vmax=vmax



    if scale=='guess':
        if (vmax-vmin)>1e4:
            scale='log'

    # Need this to go from 
    if scale=='log':
        if vmin<vmax/1e10:
            vmin=vmax/1e10
        xvar=numpy.select([xvar>vmin],[xvar],default=vmin)
        xvar=numpy.ma.log10(xvar)
        vmax=numpy.log10(vmax)
        vmin=numpy.log10(vmin)
        title='log '+title

    # Reapply the mask; need to make bad color work
    xvar=numpy.ma.array(xvar,mask=mask)
    ylabel=var



    return x,xvar,title,xlabel,ylabel


def just_plot(x,xvar,root,title,xlabel,ylabel,fig_no=1):
    '''
    This routine simply is produces a plot of a variable from
    that has been printed to an astropy table with a routine like
    windsave2table.    This function is simply a plotting routine


    '''

    pylab.figure(fig_no,(6,6))
    pylab.clf()
    pylab.title(title,size=16)
    pylab.xlabel(xlabel,size=16)
    pylab.ylabel(ylabel,size=16)
    pylab.plot(x,xvar,'-')
    pylab.tight_layout()

    pylab.draw()
    title=title.replace(' ','_')
    if root=='':
        filename=title+'.png'
    else:
        filename='%s_%s.png' %(root,title)
    pylab.savefig(filename)
    return filename





def doit(filename='7MsolBigGapEXT.0.master.txt',var='t_r',grid='ij',inwind='',scale='guess',zmin=-1e50,zmax=1e50,
        plot_dir='',root=''):

    '''
    Plot a single variable from an astropy table (normally created with windsave2table, with various
    options

    Arguments:
        var:
            is the variable to plot
        grid:
            can be ij, log, or anything else.  If ij then the plot will be in grid coordinates, if log
            the plot will be in on a log scale in physical coordiantes.  If anything else, the plot will be
            on a linear scale in physical coordiantes
        scale:
            indicates how the variable should be plotted.  guess tells the routine to make a sensible choice
            linear implies the scale should be linear and log implies a log scale should be used
        zmin:
            overide the max and mimimum in the array (assuming these limits are with the range of
            the variable)
        zmax:
            overide the max and mimimum in the array (assuming these limits are with the range of
            the variable)

    Description:

    Notes:

    History:


    '''

    if root=='':
        root=filename.split('.')
        root=root[0]

    try:
        x,xvar,title,xlabel,ylabel=get_data(filename,var,grid,inwind,scale,zmin,zmax)
    except:
        return


    if plot_dir!='':
        root='%s/%s' % (plot_dir,root)
    plotfile=just_plot(x,xvar,root,title,xlabel,ylabel)
    

    return plotfile


    





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        doit(sys.argv[1],sys.argv[2])
    else:
        print ('usage: plot_wind_1d.py filename var grid')
