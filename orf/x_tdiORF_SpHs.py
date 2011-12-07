#!/usr/bin/env python
import os
import glob
import sys
import optparse
import cPickle as cpkl
import numpy as np
import synthlisa
import myLISAmodule as mlisar


parser = optparse.OptionParser( 'usage: ./x_tdiORF_SpHs.py SETUP.pkl' )
( options , args ) = parser.parse_args()
if len( args ) < 1 :
    parser.error( 'You must specify SETUP.pkl containing input parameters!' )
file = open( args[0] , 'rb' ) ; setup = cpkl.load( file ) ; file.close()


if setup['lisa']['type'] == 'EccentricInclined' :
    eta0 , xi0 , sw , lisa_t0 = setup['lisa']['initial conditions']
    lisa = synthlisa.EccentricInclined( eta0 , xi0 , sw , lisa_t0 )

sky  = mlisar.mySpharmt( setup['sky']['nlon'] , setup['sky']['nlat'] )
t0 = setup['t0']
tdiI , tdiJ = setup['tdiI'] , setup['tdiJ']
f0 , df , Nf = setup['f0'] , setup['df'] , setup['Nf']
lmax = setup['lmax']
orfdir = setup['orfdir']

orfpath = orfdir + 'orf_t0_%f.pkl' % t0
if orfpath not in glob.glob( orfpath ) :
    f = f0 + df * np.arange( Nf )
    lisky = mlisar.LISA_in_the_Sky( lisa , sky )
#    SpHreal , SpHimag = lisky.get_SpHs( t = t0 , *( lmax , 'tdiORF' , tdiI , tdiJ , f ) )
    print lmax , tdiI , tdiJ , f , t0
    SpHreal , SpHimag = lisky.get_SpHs( lmax , 'tdiORF' , tdiI , tdiJ , f , t = t0 )
    print SpHreal
    orfdict = {'OrfMultipleMoments':
               { 'ntrunc': lmax , 'f': f , 'Antenna': tdiI[2]+tdiJ[2] ,
                 'real': SpHreal , 'imag': SpHimag } }
    if orfdir not in glob.glob( orfdir ) :
        os.system( 'mkdir -p %s' % orfdir )
    file = open( orfpath , 'wb' ) ; cpkl.dump( orfdict , file , -1 ) ; file.close()
    print " Orf multipole moments saved in %s " % orfpath
else :
    print "Orf multipole moments already saved in %s.  Nothing to do." % orfpath
