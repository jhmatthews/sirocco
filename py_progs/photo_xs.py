#!/usr/bin/env python
'''
University of Southampton -- JM -- March 2015

Synopsis:

    various utilities for processing photoionization cross-sections

Usage:

    to tabulate and save VFKY data::

        photo_xs.py tab photo_fkvy.data output_filename

    to plot all vfky xsections::

        photo_xs.py plotv photo_fkvy.data

'''

#import pylab as p
from pylab import *
import py_read_output as r
import numpy as np
import os, sys
import time

HEV = 4.13620e-15
N_VERNER_TAB = 100        # number of points in each tabulation


def sigma_phot(vfky, freq):

    '''
    sigma_phot uses Verner et al.'s interpolation formulae for the photoionization crossection
    to calculate the bound free (or photoionization) optical depth.  The data must
    have been into the photoionization structures xphot with get_atomic_data and
    the densities of individual ions must have been calculated previously.
    '''

    ix,iz = np.indices(freq.shape)

    # get arrays in right shape for frequency array
    y0 = vfky.y0[ix]
    f0 = vfky.f0[ix]
    y1 = vfky.y1[ix]
    pp = vfky.p[ix]
    ya = vfky.ya[ix]
    yw = vfky.yw[ix]
    sigma = vfky.sigma[ix]


    x = freq / f0 - y0
    y = sqrt (x * x + y1 * y1)

    # This was line fixed by CK in 1998 Jul */
    f1 = (x - 1.0) * (x - 1.0) + yw * yw

    #     f2 = pow (y, 0.5 * x_ptr->p - 5.5);
    f2 = exp ((0.5 * pp - 5.5) * log (y))

    #     f3 = pow (1.0 + sqrt (y / x_ptr->ya), -x_ptr->p);
    f3 = exp ((-pp) * log ((1.0 + sqrt (y / ya))))

    xsection = sigma * f1 * f2 * f3    # the photoinization xsection

    return xsection




class Photo(object):
    '''This is a general class for photoionization data'''

    def __init__(self):
        self.type = None
        self.z = None
        self.ion = None
        self.islp = None
        self.l = None
        self.E0 = None
        self.np = None
        self.energy = None
        self.XS = None
        self.fname = None

    def read_topbase_file(self, filename, mode = "Top"):
        '''
        read in XS info from Topbase XS data in Python format

        INPUT:

            filename (string):
                atomic data filename e.g. topbase_h1_phot.py

        OUTPUT:

            top (topbase class instance):
                topbase class instance containing information for this filename

        '''
        self.fname = filename
        self.type = "Topbase"
        self.tabulated = True

        # read in summary records
        z,ion,islp,l, E0, num_records = sum_records = np.loadtxt(filename,
                dtype={'names': ('Z', 'ion', 'islp', 'l', 'E0', 'np'),
                'formats': ('i4', 'i4', 'i4', 'i4', 'float', 'i4')},
                            comments='Phot%s ' % mode , delimiter=None, converters=None,
                            skiprows=0, usecols=(1,2,3,4,5,6), unpack=True, ndmin=0)

        # then read the actual cross sections
        energy, XS = np.loadtxt(filename, dtype='float',
                        comments='Phot%sS' % mode, delimiter=None, converters=None,
                        skiprows=0, usecols=(1,2), unpack=True, ndmin=0)

        nline = 0

        for i in range(len(z)):

            n_p = int(num_records[i])
            nmax = nline + n_p

            self.z[i] = z[i]
            self.ion[i] = ion[i]
            self.islp[i] = islp[i]
            self.l[i]= l[i]
            self.E0[i] = E0[i]
            self.np[i] = num_records[i]
            self.energy[i] = energy[nline:nmax]
            self.XS[i] = XS[nline:nmax]

            nline = nmax

        return 0



    def read_vfky_file(self, filename):
        '''
        read in XS info from Verner XS data in Python format

        INPUT:

            filename (string):
                atomic data filename e.g. `topbase_h1_phot.py`

        OUTPUT:

            populates members of class such as z, f0, etc.
        '''

        self.fname = filename

        # read in vfky data
        self.z, nelectrons, et, emax, self.E0, sigma, self.ya, self.p, self.yw, self.y0, self.y1 = np.loadtxt(filename, dtype={'names': ('Z', 'ne', 'et', 'emax','e0', 'sigma', 'ya', 'p', 'yw', 'y0', 'y1'), 'formats': ('i4', 'i4', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float')}, comments='#', delimiter=None, converters=None, skiprows=0, usecols=(1,2,3,4,5,6,7,8,9,10,11), unpack=True, ndmin=0)

        self.ion = self.z - nelectrons + 1
        self.ft = et / HEV
        self.fmax = emax / HEV
        self.f0 = self.E0 / HEV
        self.sigma = sigma * 1e-18

        self.type = "Vfky"

        return 0

    def write_file(self, output_filename):
        '''
        write a class out to a tabulated file
        '''

        if self.tabulated == False or self.type == None:
            print("Error: Must have read and tabulated data before writing to file! run tabulate_vfky()")
            return 0

        if self.type == "Topbase":
            prefix = "PhotTop"
        elif self.type == "Vfky":
            prefix = "PhotVfky"

        f = open(output_filename, "w")
        for i in range(len(self.z)):

            f.write("%sS %i %i %i %i %8.4e %i\n" % (prefix, self.z[i], self.ion[i], self.islp[i], self.l[i], self.E0[i], self.np[i]))

            for j in range(len(self.energy[i])):

                f.write("%s %8.4e %8.4e\n" % (prefix, self.energy[i][j], self.XS[i][j]) )

        return 0


    def tabulate_vfky(self):
        '''
        read in XS info from Topbase XS data in Python format

        INPUT:

            self

        :OUTPUT:

            top (topbase class instance):
                topbase class instance containing information for this filename
        '''

        # We need to start our tabulation just a tiny way up from from the threshold, otherwise it is equal to zero.
        very_small = 1e-6

        f1 = self.ft * (1 + very_small)
        f2 = self.fmax * (1 + very_small)

        lf1 = np.log(f1)
        lf2 = np.log(f2)

        dlogf=(lf2-lf1)/(N_VERNER_TAB)

        # create an index array to help with quick manipulation of tabulation
        ion_index, freq_index = np.indices( (len(f1), N_VERNER_TAB ) )

        # this next line does the actual frequency tabulation- quite a lot of numpy trickery going on here
        freq = np.exp(lf1[ion_index] + (dlogf[ion_index]*freq_index))

        # energy in EV
        self.energy = HEV * freq

        # calculate cross_sections
        self.XS = sigma_phot(self, freq)

        self.islp = np.ones(len(self.z)) * -999        # islp not relevant to Verner
        self.l = np.ones(len(self.z)) * -999        # islp not relevant to Verner
        self.np = np.ones(len(self.z)) * N_VERNER_TAB    # all np values the same...not really necessary but more flexible

        self.tabulated = True

        return 0

    def plot_all(self):
        '''plot all Xsections - makes a lot of plots!!'''

        if self.tabulated == False or self.type == None:
            print("Error: Must have read and tabulated data before plotting! run tabulate_vfky()")
            return 0

        #rcParams["text.usetex"] = "True"
        for i in range(len(self.z)):

            plot(self.energy[i]/HEV, self.XS[i], linewidth = 2)
            xlabel(r"nu", fontsize = 16)
            ylabel(r"sigma", fontsize = 16)
            loglog()

            if self.type == "Vfky":
                title("Vfky Z=%i, ion=%i" % (self.z[i], self.ion[i]))
                savefig("xs_vfky_z%i_i%i.png" % (self.z[i], self.ion[i]))

            clf()

        return 0

    def associate_levels(self, levels):

        ground_state_select = (levels.et == 0.0)


        for i in range(len(self.z)):

            z_select = (levels.z == self.z[i])
            ion_select = (levels.ion == self.ion[i])
            select = z_select * ion_select * ground_state_select

            islp = levels.islp[select]
            nlev = levels.l[select]
            z = levels.z[select]
            et = levels.et[select]
            ion = levels.ion[select]

            if len(z) == 1:
                print("Match: %i %i %i %i %8.4e" % (z[0], ion[0], islp[0], nlev[0], et[0]))
                self.islp[i] = float(islp)
                self.l[i] = float(nlev)
            elif len(z) < 1:
                print("No matches for ion %i %i. Not writing level info." % (self.z[i], self.ion[i]))
            else:
                print("Multiple matches for ion %i %i. Not writing level info." % (self.z[i], self.ion[i]))

        return 0



class TopBaseLevel(object):

    def __init__(self, fname):

        self.z, self.ion, self.islp, self.l, self.E, self.et = sum_records = np.loadtxt(fname,
                dtype={'names': ('Z', 'ion', 'islp', 'l', 'E0', 'np'),
                'formats': ('i4', 'i4', 'i4', 'i4', 'float', 'float')},
                            comments='#', delimiter=None, converters=None,
                            skiprows=0, usecols=(1,2,3,4,5,6), unpack=True, ndmin=0)









# Next lines permit one to run the routine from the command line with various options -- see docstring
if __name__ == "__main__":

    if len(sys.argv)<3:
        print(__doc__)
        sys.exit(1)

    mode = sys.argv[1]

    if mode == "tab":
        input_fname = sys.argv[2]
        output_fname = sys.argv[3]

        v = Photo()

        v.read_vfky_file(input_fname)
        v.tabulate_vfky()

        fnames = ["topbase_levels_he.py", "topbase_levels_h.py", "topbase_levels_cno.py", "topbase_levels_fe.py"]

        for f in fnames:
            lev = TopBaseLevel("data/atomic77/%s" % f)
            v.associate_levels(lev)

        v.write_file(output_fname)

        print("Tabulated and saved output in %s" % output_fname)

    elif mode == "plotv":
        input_fname = sys.argv[2]

        v = Photo()

        v.read_vfky_file(input_fname)
        v.tabulate_vfky()
        v.plot_all()

        print("Plotted all tabulated XS")

    else:
        print(__doc__)
        sys.exit(1)
