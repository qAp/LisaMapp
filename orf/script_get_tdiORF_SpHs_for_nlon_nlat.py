#!/usr/bin/env python
import os
import glob
import sys
import optparse
import cPickle as cpkl
import numpy as np
import synthlisa
import myLISAmodule as mlisar

"""
Define LISA here
For synthlisa.EccentricInclined(), lisa = { 'type':'EccentricInclined' , 'initial conditions':( eta0 , xi0 , sw , lisa_t0 ) } , where
eta0 --- true anomaly of LISA's guiding centre
xi0  --- rotational phase of the LISA array
sw --- sw < 0 swaps spacecraft so LISA becomes 1->3->2->1 viewed from above
t0 --- the time at which LISA has the above
"""
eta0 , xi0 , sw , lisa_t0 = 0 , 0 , 1 , 0
lisa = { 'type':'EccentricInclined' , 'initial conditions':( eta0 , xi0 , sw , lisa_t0 ) }

"Define the sky"
nlon , nlat = 100 , 51 ; sky = { 'nlon':nlon , 'nlat':nlat }

"Time at which to evaluate LISA parameters"
t0 = 100.

"Specify which tdiORF"
tdiI , tdiJ = ( 'Michelson' , 'G0' , '1' , '1' ) , ( 'Michelson' , 'G0' , '1' , '1' )

"Which frequencies?"
f0 , df , Nf = 1e-5 , 1e-5 , int( 1e5 - 1 )

"lmax"
lmax = 20

"Where to save?"
orfdir = 'data_nlon_%d_nlat_%d/' % ( nlon , nlat )


setup = { 'lisa':lisa , 'sky':sky , 't0':t0 , 
          'tdiI':tdiI , 'tdiJ':tdiJ ,
          'f0':f0 , 'df':df , 'Nf':Nf ,
          'lmax':lmax ,
          'orfdir':orfdir }

setupname = 'setup_tdiORF_SpHs_for_nlon_%d_nlat_%d.pkl' % ( nlon , nlat )
file = open( setupname , 'wb' ) ; cpkl.dump( setup , file , -1 ) ; file.close()

submitname = 'get_tdiORF_SpHs_for_nlon_%d_nlat_%d.sub' % ( nlon , nlat )
file = open( submitname , 'w' )
file.writelines( [ '#!/bin/bash\n' ,
                   '#PBS -N %s\n' % submitname ,
                   '#PBS -q compute\n' ,
                   '#PBS -j oe\n' ,
                   '#PBS -l nodes=2:ppn=2\n' ,
                   '#PBS -l walltime=5:00:00\n' ,
                   'cd $PBS_O_WORKDIR\n' ,
                   '\n' ,
                   './x_tdiORF_SpHs.py %s\n' % setupname ] ) ; file.close()


os.system( 'qsub %s' % submitname )
