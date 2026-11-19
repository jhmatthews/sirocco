#!/usr/bin/env python 

'''
This routine runs a set of models.

Each model is run in single processor mode, so the parallelism here is job parallel.


Command line usage (if any):

    usage::

        run_many.py -jobs x -np y filename

    `-njobs`
        is the number of jobs to run simultaneously
    `-np`
        is the number of processors per job
    `filename`
        is a list of the .pf files to run


Description:  

Primary routines:

    doit  - oversees multiprocesing of the jobs
    run_one - runs each job
    steer - parses the command line

Notes:

    At present, there is no checking that the various .pf files
    do not drop into interactive mode, or that the routine stalls.
    This really should be added.   Also, we don't currently gather
    information on whether the runs succeeded.
                                       
History:

    190319 ksl Coding begun

'''

import sys
from astropy.io import ascii
import numpy
import time
import multiprocessing
# import date
import subprocess



def read_file(filename,char=''):
    '''
    Read a file and split it into words, eliminating comments
    
    char is an optional parameter used as the delimiter for
    splitting lines into words.  Otherwise white space is
    assumed.

    History:
    
    110729    ksl    Added optional delimiters
    141209    ksl    Reinstalled in my standard startup script so there was flexibility to read any ascii file
    '''

    try:
        f=open(filename,'r')
        xlines=f.readlines()
        f.close()
    except IOError :
        print ("The file %s does not exist" % filename)
        return []   
    
    lines=[]
    
    i=0
    while i<len(xlines):
        z=xlines[i].strip()
        if char=='':
            z=z.split()
        else:
            z=z.split(char)
        if len(z)>0:
            if z[0][0]!='#':
                lines=lines+[z]
        i=i+1
    return lines


def run_one(filename,command='py'):
    '''
    Run and verify a single model
    '''

    xcommand='%s %s' % (command,filename) 

    print(xcommand)

    proc=subprocess.Popen(xcommand,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout,stderr=proc.communicate()

def get_no_jobs(jobs):
    '''
    Check how many jobs are running
    '''
    njobs=0
    for one in jobs:
        if one.is_alive():
            njobs=njobs+1
    return njobs


def doit(filename='all.ls',nj=2,np=1,outputfile='out.txt'):
    '''
    Do something magnificent

    Description:

    Notes:

    History:


    '''

    print('List =',filename)


    command='nice 5 py -p '
    if np>1:
        command='nice 5 mpirun -np %d py -p ' % np


    try:
        # master=read_file(filename)
        
        master=ascii.read(filename,format='no_header')
        master.rename_column('col1','File')
    except:
        print('Error: Could not read %s' % filename)
        return

    print(master)


    # for one in master:
    #     print(one['File'])


    #     run_one(one['File'])

    
    jobs=[]
    for one in master:
        p=multiprocessing.Process(target=run_one,args=(one['File'],command))
        jobs.append(p)

    i=0
    while i<nj and i<len(jobs):
        print('!Starting %s' % master['File'][i])
        one=jobs[i]
        one.start()
        i=i+1

    njobs=get_no_jobs(jobs)
    while i<len(jobs):
        time.sleep(10)
        njobs=get_no_jobs(jobs)
        print('Running %d jobs,including job %d (%s) of %d total' % (njobs,i,master['File'][i-1],len(master)))
        if njobs<nj:
            print('Starting: ',master['File'][i])
            one=jobs[i]
            one.start()
            i=i+1

    for one in jobs:
        one.join()
    # p.join()
    print('Completed multiprocessing')


    return


def steer(argv):
    '''
    Process the command line
    '''

    njobs=4
    np=1
    models=''

    i=0
    while i < len(argv):
        if argv[i]=='-h':
            print(__doc__)
            return
        elif argv[i]=='-np':
            i=i+1
            np=int(argv[i])
        elif argv[i]=='-jobs':
            i+=1
            njobs=int(argv[i])
        elif argv[i][0]=='-':
            print('Error: Unknown switch ---  %s' % argv[i])
            return
        else:
            models=argv[i]
        i+=1

    if models=='':
        print('No models to process without a list of models')
        return

    doit(models,njobs,np)







# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        # doit(sys.argv[1])
        steer(sys.argv)
    else:
        print (__doc__)
