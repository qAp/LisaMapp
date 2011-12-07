#!/usr/bin/env python
import os
import glob
import sys
import cPickle as cpkl
import numpy as np
#from LISAresponse import *
import synthlisa
import myLISAmodule as mlisar



"""
Define LISA here
eta0 --- true anomaly of LISA's guiding centre
xi0  --- rotational phase of the LISA array
sw --- sw < 0 swaps spacecraft so LISA becomes 1->3->2->1 viewed from above
t0 --- the time at which LISA has the above
"""
eta0 , xi0 , sw , t0 = 0 , 0 , 1 , 0
lisa = synthlisa.EccentricInclined( eta0 , xi0 , sw , t0 )

"Define the sky"
nlon , nlat = 40 , 21
sky  = mlisar.mySpharmt( nlon , nlat )

"Which first & last days?"
days = [ 180 ]

"Specify the function over the sky"
whichfunc = 'tdiORF'
tdiI , tdiJ = ( 'Michelson' , 'G2' , 'A' , '1' ) , ( 'Michelson' , 'G2' , 'E' , '1' )

"Which frequencies?"
f0 , df , Nf = 1e-5 , 1e-5 , 99999
f = f0 + df * np.arange( Nf )

"lmax"
lmax = nlat - 1

"Where to save?"
#orfdir = 'tdiI_%s_%s_%s_tdiJ_%s_%s_%s_lmax_%d_f0_%.1e_df_%.1e_Nf_%.1e/data/' % ( tdiI[:-1] + tdiJ[:-1] + ( lmax , f0 , df , Nf ) )
orfdir = 'test/'
if orfdir not in glob.glob( orfdir ) :
    os.system( 'mkdir -p %s' % orfdir )



dinsecs = 24*60.**2

lisky = mlisar.LISA_in_the_Sky( lisa , sky )

for day in days :
    print '>>>>>> Day %3d <<<<<<' % day
    t = ( day - 0.5 )*dinsecs

    SpHreal , SpHimag = lisky.get_SpHs( t = t , *( lmax , whichfunc , tdiI , tdiJ , f ) )

    orfdic = {'LISAData':
              {'Type': 'EccentricInclined' ,
               'Parameters':
               {'eta0': eta0 , 'xi0': xi0 , 'sw': sw , 't0': t0} } ,
              'OrfMultipleMoments':
              { 'ntrunc': lmax , 'f': f , 'Antenna': tdiI[2]+tdiJ[2] ,
                'real': SpHreal , 'imag': SpHimag }
              }


    filename = 'orf-d%03d.pkl' % day
    filepath = orfdir + filename
    file = open( filepath , 'wb' )
    cpkl.dump( orfdic , file , -1 )
    file.close()

    print " Orf multiple moments saved in %s " % filepath
    

                                
