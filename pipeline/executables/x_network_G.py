#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import AnisotropySearch as AS
import optparse

usage = """
$prog [ GPATH,... ] NETWORKGPATH\n
GPATH --- path to file of Gs to be summed(can have more than one)
NETWORKGPATH --- path to file in which to save network G
"""
parser = optparse.OptionParser( usage = usage )

( options , args ) = parser.parse_args()

if len( args ) < 2 :
    parser.error( 'You must specify at least an GPATH and a NETWORKPATH! ( See help: x_network_G.py -h )' )
else :
    Gpaths = args[ : -1 ] ; netGpath = args[ -1 ]

for i , Gpath in enumerate( Gpaths ) :
    file = open( Gpath , 'rb' ) ; Gdict = cpkl.load( file ) ; file.close()
    if i == 0 :
        netGdata = np.copy( Gdict['G'].data )
        lmax = Gdict['ntrunc']
    else :
        netGdata += np.copy( Gdict['G'].data )

netG = AS.Coarsable( netGdata )
Gdict = {} ; Gdict['G'] = netG ; Gdict['ntrunc'] = lmax

netGdir = os.path.dirname( netGpath )
if netGdir not in glob.glob( netGdir ) :
    os.system( 'mkdir -p %s' % netGdir )
file = open( netGpath , 'wb') ; cpkl.dump( Gdict , file , -1 ) ; file.close()
