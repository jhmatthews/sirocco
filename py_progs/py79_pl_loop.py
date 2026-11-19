#!/usr/bin/env python


'''
This routine carries out a series of thin shell python simulations.

The wind mdot is set to produce a hydrogen density of 1e7 - to change
this one had to change the wind_mdot parameter. The loop is carried out
over the 2-10kev luminosity of the central source. A luminosity of
1e21 give a mainly neutral plasma (IP=5e-10) whilst 1e37 is almost 
totally ionized (IP=5e6). The code strips out ion densities and constructs 
data files with the cloudy type ionization parameter as the first 
line.
If one wants a constant temperature file, then one needs to hack python
to stop the temperature changing
We currently write out H,He,C,N,O and Fe. We also output heating and cooling
mechanisms.


Command line usage (if any):

	usage: python_pl_loop *optional file suffix*

Description:  

Primary routines:

Notes:
									   
History:

081214 nsh Commented

'''


if __name__ == "__main__":		# allows one to run from command line without running automatically with write_docs.py

	import sys, subprocess
	import numpy as np

	#Use an optional suffix so file will be py_hydrogen_*suffix*.dat, 
	#If nothing is supplied, default is PL.
	
	alpha=-0.9
	npoints=101 			#this is the parameter which says how many runs to perform - we do 10 per order of magnitude in lum
	python_ver="py"   	#This is where we define the version of python to use
	atomic="data/standard78"
	python_opts=" "
	ncycles=20
	nphot=1000000
	t_e=10000
	nh=1e7
	lum_start=1e25
	nprocs=0
	
	
	
	
	
	
	wind_mdot=4.72694719429145E-20*(nh/1e7)


	if  len(sys.argv) > 1:
		run_name=sys.argv[1]
		inp=open(run_name+".param")
		for line in inp.readlines():
			data=line.split()
			if len(data)!=2:
				print ("Improperly structured param file - each line should be a param name and value seperated by a space")
			else:
				print(data[0],data[1])
				if data[0]=='alpha': 
					alpha=data[1]
				if data[0]=='nprocs': 
					nprocs=int(data[1])
				if data[0]=='npoints':
					npoints=int(data[1])
				if data[0]=='t_e':
					t_e=float(data[1])
				if data[0]=='python_ver':
					python_ver=data[1]
				if data[0]=='atomic':
					atomic=data[1]
				if data[0]=='python_opts':
					python_opts=data[1]
				if data[0]=='ncycles':
					ncycles=int(data[1])
				if data[0]=='nphot':
					nphot=data[1]
				if data[0]=='nh':
					nh=data[1]
				if data[0]=='lum_start':
					lum_start=float(data[1])
	else:
		run_name='PL'


		
	

	#Open empty files to contain data

	Hout=open('py_hydrogen_'+run_name+'.dat','w')
	Hout.write('IP H1 H2\n')
	Heout=open('py_helium_'+run_name+'.dat','w')
	Heout.write('IP He1 He2 He3\n')
	Cout=open('py_carbon_'+run_name+'.dat','w')
	Cout.write('IP C1 C2 C3 C4 C5 C6 C7\n')
	Nout=open('py_nitrogen_'+run_name+'.dat','w')
	Nout.write('IP N1 N2 N3 N4 N5 N6 N7 N8\n')
	Oout=open('py_oxygen_'+run_name+'.dat','w')
	Oout.write('IP O1 O2 O3 O4 O5 O6 O7 O8 O9\n')
	Feout=open('py_iron_'+run_name+'.dat','w')
	Feout.write('IP')
	for i in range(27):
			Feout.write(" Fe"+str(i+1))
	Feout.write('\n')
	Tout=open('py_temperature_'+run_name+'.dat','w')
	Tout.write("Lum Sim_ip Cl_ip T_e N_e heat cool\n")
	Heat=open('py_heat_'+run_name+'.dat','w')
	Heat.write('IP total photo ff compton ind_comp lines auger\n')
	Cool=open('py_cool_'+run_name+'.dat','w')
	Cool.write('IP total recomb ff compton DR DI lines Adiabatic\n')



	for i in range(npoints):
		lum=10**((float(i)+np.log10(lum_start)*10.0)/10.0)   #The 210 means the first luminosity is 21.0
		print('Starting cycle '+str(i+1)+' of '+str(npoints))
		print('Lum= '+str(lum))
		inp =open('input.pf','w')
		inp.write("Wind_type() 9\n")
		inp.write("Atomic_data "+atomic+"\n")
		inp.write("write_atomicdata      0\n")
		inp.write("photons_per_cycle "+str(nphot)+"\n")
		inp.write("Ionization_cycles "+str(ncycles)+"\n")
		inp.write("spectrum_cycles                                   0\n")
		inp.write("Coord.system()                   0\n")
		inp.write("Wind.dim.in.x_or_r.direction                     4\n")
		inp.write("Wind_ionization()                   9\n")
		inp.write("Line_transfer()                    3\n")
		inp.write("Thermal_balance_options 		0\n")
		inp.write("System_type()                   2\n")
		inp.write("disk.type()             							0\n")
		inp.write("Star_radiation(y=1)                              0\n")
		inp.write("Disk_radiation(y=1)                              0\n")
		inp.write("Wind_radiation(y=1)                              0\n")
		inp.write("QSO_BH_radiation(y=1)                            1\n")
		inp.write("Rad_type_for_star(0=bb,1=models)_to_make_wind                   0\n")
		inp.write("Rad_type_for_agn()_to_make_wind                   4\n")
		inp.write("mstar(msol)                                     0.8\n")
		inp.write("rstar(cm)                                     1e10\n")
		inp.write("tstar                                      1000000\n")
		inp.write("disk.type()                   0\n")
		inp.write("lum_agn(ergs/s) "+str(lum)+"\n")
		inp.write("agn_power_law_index  "+str(alpha)+"\n")
		inp.write("low_energy_break(ev)				0.136\n")
		inp.write("high_energy_break(ev)			20000\n")
		inp.write("Torus(0=no,1=yes) 				 0\n")
		inp.write("wind.radmax(cm)                   1.00000000001e11\n")
		inp.write("wind.t.init       "+str(t_e)+"\n")
		inp.write("shell_wind_mdot(msol/yr)  "+str(wind_mdot)+"\n")
		inp.write("shell.wind.radmin(cm)                       1e11\n")
		inp.write("shell.wind_v_at_rmin(cm)                    1.00000\n")
		inp.write("shell.wind.v_at_rmax(cm)                    1.000010\n")
		inp.write("shell.wind.acceleration_exponent                   1\n")
		inp.write("wind.filling_factor(1=smooth,<1=clumped) (1)       1\n")
		inp.write("spec.type(flambda(1),fnu(2),basic(other)                    2\n")
		inp.write("reverb.type()									0\n")
		inp.write("Extra.diagnostics(0=no)                           0\n")
		inp.write("Use.standard.care.factors(1=yes)                    1\n")
		inp.write("Photon.sampling.approach()                   5\n")
		inp.write("Num.of.frequency.bands                           3\n")
		inp.write("Lowest_energy_to_be_considered(eV)                   0.0001\n")
		inp.write("Highest_energy_to_be_considered(eV)           100000000   \n")
		inp.write("Band.boundary(eV)                             2000\n")
		inp.write("Band.boundary(eV)                            10000\n")
		inp.write("Band.minimum_fraction)                         0.3\n")
		inp.write("Band.minimum_fraction)                         0.4\n")
		inp.write("Band.minimum_fraction)                         0.3\n")
		inp.close()
		if nprocs==0:
			cmd="time "+python_ver+" "+python_opts+" input > output"
		else:
			cmd="time mpirun -n "+str(nprocs)+" "+python_ver+" "+python_opts+" input > output"
		print(cmd)
		subprocess.check_call(cmd,shell=True)	   #This is where we define the version of python to use
		subprocess.check_call("tail -n 60 output  | grep OUTPUT > temp",shell=True)#Strip the last 60 lines from the output
		inp=open('temp','r')
		for line in inp.readlines():
			data=line.split()
			if (data[1]=='Lum_agn='):       #This marks the start of the data we want. Field 14 is the cloudy ionization parameter
				Hout.write(data[14])
				Heout.write(data[14])
				Cout.write(data[14])
				Nout.write(data[14])
				Oout.write(data[14])
				Feout.write(data[14])
				Heat.write(data[14])
				Cool.write(data[14])
				Tout.write(data[2]+' '+data[12]+' '+data[14]+' '+data[4]+' '+data[8]+' ')
			if (data[1]=='Absorbed_flux(ergs-1cm-3)'):
				Heat.write(' '+data[2]+' '+data[4]+' '+data[6]+' '+data[8]+' '+data[10]+' '+data[12]+' '+data[14]+'\n')
				Tout.write(data[2]+' ')
			if (data[1]=='Wind_cooling(ergs-1cm-3)'):
				Cool.write(' '+data[2]+' '+data[4]+' '+data[6]+' '+data[8]+' '+data[10]+' '+data[12]+' '+data[16]+' '+data[14]+'\n')
				Tout.write(data[2]+'\n')
			if (data[1]=='H'):
				for j in range(2):
					Hout.write(' '+data[j+2])
				Hout.write('\n')
			if (data[1]=='He'):
				for j in range(3):
					Heout.write(' '+data[j+2])
				Heout.write('\n')
			if (data[1]=='C'):
				for j in range(7):
					Cout.write(' '+data[j+2])
				Cout.write('\n')
			if (data[1]=='N'):
				for j in range(8):
					Nout.write(' '+data[j+2])
				Nout.write('\n')
			if (data[1]=='O'):
				for j in range(9):
					Oout.write(' '+data[j+2])
				Oout.write('\n')
			if (data[1]=='Fe'):
				for j in range(27):
					Feout.write(' '+data[j+2])
				Feout.write('\n')
	#Flush the output files so one can see progress and if the loop crashes all is not lost
		Hout.flush()
		Heout.flush()
		Cout.flush()
		Nout.flush()
		Oout.flush()
		Feout.flush()
		Tout.flush()
		Heat.flush()
		Cool.flush()	
		print('Finished cycle '+str(i+1)+' of '+str(npoints))
	#Close the files.
	Hout.close()
	Heout.close()
	Cout.close()
	Nout.close()
	Oout.close()
	Feout.close()
	Tout.close()
	Heat.close()
	Cool.close()

