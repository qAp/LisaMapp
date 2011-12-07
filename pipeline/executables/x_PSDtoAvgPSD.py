#!/usr/bin/env python
import os
import glob
import sys
import AnisotropySearch as AS
import cPickle as cpkl
from optparse import OptionParser



parser = OptionParser( "usage: ./x_PSDtoAvgPSD.py PSDDIR AVGPSDDIR" )
parser.add_option( '-d' , '--day' , action='append' , dest='days' , type='int' , nargs=1 ,
                   help = 'Days on which to calculate average PSDs.' )

( options , args ) = parser.parse_args()

if len( args ) < 2 :
    parser.error( "You must specify PSDDIR, directory containing PSDs, and AVGPSDDIR, directory to save average PSDs in!" )

psddir , avgpsddir = args


if avgpsddir not in glob.glob( avgpsddir ) :
    os.system( 'mkdir %s' % avgpsddir )

for day in options.days :

    print "++++ DAY %d ++++" % day

    dayl , dayr = day - 1 , day + 1
    psdpathl , psdpathr = ( psddir + 'd%03d.pkl' % dayl ) , ( psddir + 'd%03d.pkl' % dayr )

    if psdpathl not in glob.glob( psdpathl ) :
        print 'Left PSD not available on day %d. Nothing to do.' % day
        continue
    if psdpathr not in glob.glob( psdpathr ) :
        print 'Right PSD not available on day %d. Nothing to do.' % day
        continue

    print "PSDs of adjacent days available, taking the average...",
    AS.avg_psdPKL_avg( day , psdpathl , psdpathr , avgpsddir )
    print "done"

