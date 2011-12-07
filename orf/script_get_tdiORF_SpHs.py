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
nlon , nlat = 40 , 21 ; sky = { 'nlon':nlon , 'nlat':nlat }

"Time at which to evaluate LISA parameters"
t0 = 100.

"Specify which tdiORF"
tdiI , tdiJ = ( 'Michelson' , 'G0' , '1' , '1' ) , ( 'Michelson' , 'G0' , '1' , '1' )

"Which frequencies?"
f0 , df , Nf = 1e-5 , 1e-5 , int( 1e5 - 1 )
print Nf
"lmax"
lmax = 4

"Where to save?"
orfdir = 'tdiI_%s_%s_%s_tdiJ_%s_%s_%s_lmax_%d_f0_%.1e_df_%.1e_Nf_%d/data_nlon_%d_nlat_%d/' % ( tdiI[:-1] + tdiJ[:-1] + ( lmax , f0 , df , Nf , nlon , nlat ) )


setup = { 'lisa':lisa , 'sky':sky , 't0':t0 , 
          'tdiI':tdiI , 'tdiJ':tdiJ ,
          'f0':f0 , 'df':df , 'Nf':Nf ,
          'lmax':lmax ,
          'orfdir':orfdir
          }

setupname = 'setup_tdiORF_SpHs.pkl'
file = open( 'setup_tdiORF_SpHs.pkl' , 'wb' ) ; cpkl.dump( setup , file , -1 ) ; file.close()



os.system( './x_tdiORF_SpHs.py %s' % 'setup_tdiORF_SpHs.pkl' )
