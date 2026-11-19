#!/usr/bin/env python 
'''
This is a set of classes for use with the radiative transfer code python

Usage:
	

Arguments:

Notes: 
	Not all fo these classes are necessarily used by the plotting routines
	These are slightly less relevant now we've started using astropy.io.ascii
'''




class spectotclass:
    '''This is a class for storing any values read from a Python spec_tot file'''	
    def __init__(self, fre, wave, create, emi, cen, dis, win, sca, hit):
        self.freq = fre
        self.wavelength = wave
        self.created = create
        self.emitted = emi
        self.censrc = cen
        self.disk = dis
        self.wind = win
        self.scattered = sca
        self.hitSurf = hit



class specclass:
    '''This is a class for storing any values read from a Python spec file'''	
    def __init__(self, fre, wave, prod, emi, cen, dis, win, sca, hit, spe):
        self.freq = fre
        self.wavelength = wave
        self.created = prod
        self.emitted = emi
        self.censrc = cen
        self.disk = dis
        self.wind = win
        self.scattered = sca
        self.hitSurf = hit
        self.spec = spe



# line class: analogous to line ptr in python. contains freq, oscillator strength, 
class line:
	'''This is a class analogous to line ptr in python'''
	def __init__(self, _z, _ion, _wavelength, _freq, _osc, _g_l, _g_u, _ll, _lu):
		self.z = _z
		self.ion = _ion
		self.wavelength = _wavelength
		self.freq = _freq
		self.osc = _osc
		self.g_l = _g_l
		self.g_u = _g_u
		self.ll = _ll
		self.lu = _lu
		
		
		
# line class: analogous to line ptr in python. contains freq, oscillator strength, 
class level:
	'''Stores information from a Python levels file'''
	def __init__(self, _z, _ion, _lvl, _ionpot, _E, _g, _rad_rate, _bracks, _nnstring):
		self.z = _z
		self.ion = _ion
		self.lvl = _lvl
		self.ionpot = _ionpot
		self.E = _E
		self.g = _g
		self.rad_rate = _rad_rate
		self.bracks = _bracks
		self.nnstring = _nnstring
		
class chianti_level:
	'''
	Stores chianti level data from file format .elvlc
	See http://www.chiantidatabase.org/cug.pdf Page 8 for description
	'''
	def __init__(self, _index,  _config, _notation, _spin, _l, _l_symbol, _j, _multiplicity, _E_obs, _E_obs2, _E_th, _E_th2, _n):
		self.index = _index
		self.config = _config
		self.notation = _notation
		self.spin = _spin
		self.l = _l
		self.l_symbol = _l_symbol 
		self.j = _j
		self.multiplicity = _multiplicity
		self.E_obs = _E_obs
		self.E_obs2 = _E_obs2
		self.E_th = _E_th
		self.E_th2 = _E_th2
		self.n = _n
		
class chianti_rad:
	'''
	Stores radiative chianti information from wgfa file.
	Contains the wavelengths, gf and A values of the transitions and the indices initial and
	final level corresponding to the indices of the levels as given in the .elvlc file
	See http://www.chiantidatabase.org/cug.pdf Page 9 for description
	'''
	def __init__(self, _ll, _lu, _wave, _freq, _osc, _A, _note_low, _note_up):
		self.ll = _ll
		self.lu = _lu
		self.wave = _wave
		self.freq = _freq
		self.osc = _osc
		self.A = _A
		self.note_low = _note_low
		self.note_up = _note_up
		
		
		
		
		



