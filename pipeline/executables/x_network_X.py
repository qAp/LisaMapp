#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import AnisotropySearch as AS
import optparse

usage = """
$prog [ XPATH,... ] NETWORKXPATH\n
XPATH --- path to file of Xs to be summed(can have more than one)
NETWORKXPATH --- path to file in which to save network X
"""
parser = optparse.OptionParser( usage = usage )

( options , args ) = parser.parse_args()

if len( args ) < 2 :
    parser.error( 'You must specify at least an XPATH and a NETWORKXPATH! ( See help: x_network_X.py -h )' )
else :
    Xpaths = args[ : -1 ] ; netXpath = args[ -1 ]

for i , Xpath in enumerate( Xpaths ) :
    file = open( Xpath , 'rb' ) ; X = cpkl.load( file ) ; file.close()
    if i == 0 :
        netXdata = np.copy( X.xlm )
        lmax = X.ntrunc
    else :
        netXdata += np.copy( X.xlm )

netX = AS.xlmSkyMap( xlm = netXdata )

netXdir = os.path.dirname( netXpath )
if netXdir not in glob.glob( netXdir ) :
    os.system( 'mkdir -p %s' % netXdir )

file = open( netXpath , 'wb') ; cpkl.dump( netX , file , -1 ) ; file.close()
