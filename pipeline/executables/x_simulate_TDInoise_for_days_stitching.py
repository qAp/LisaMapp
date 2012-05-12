#!/usr/bin/env python
import os
import sys
import glob
import optparse
import cPickle as cpkl
import numpy as np
import myLISAmodule as mlisar
import AnisotropySearch as AS
import myUsefuls as mufls


def get_comatrix( f , tditype , tdigen , whichtdis ) :
    """
    Returns the covariance matrix of a set of TDI observables
    INPUT:
    f --- frequencies (numpy array)
    tditype --- type of TDI observables ('Michelson','Sagnac')
    tdigen --- generation of TDI observables ('G0','G1','Gm','G2')
    whichtdis --- set of TDI observables ('123','AET')
    OUTPUT:
    comatrix --- covariance matrix ( Nvar x Nvar x Nf numpy array)
    """
    if whichtdis == '123' :
        pairs = [ ('1','1') , ('1','2') , ('1','3') ,
                  ('2','1') , ('2','2') , ('2','3') ,
                  ('3','1') , ('3','2') , ('3','3') ]
        cIJs = [ mlisar.get_tdiNSD( tditype , tdigen , pair[0] , pair[1] , f ) for pair in pairs ]
    elif whichtdis == 'AET' :
        pairs = [ ('A','A') , ('A','E') , ('A','T') ,
                  ('E','A') , ('E','E') , ('E','T') ,
                  ('T','A') , ('T','E') , ('T','T') ]
        cIJs = [ mlisar.get_tdiNSD( tditype , tdigen , pair[0] , pair[1] , f ) for pair in pairs ] 
    comatrix = np.array( [ [ cIJs[0] , cIJs[1] , cIJs[2] ] ,
                           [ cIJs[3] , cIJs[4] , cIJs[5] ] ,
                           [ cIJs[6] , cIJs[7] , cIJs[8] ] ] )
    return comatrix



usage = """
%prog TSDIR\n
TSDIR --- directory in which to save the time-series files
"""
parser = optparse.OptionParser( usage=usage )
parser.add_option( '-d' , '--day' , action='append' , dest='days' , type='int' , nargs=1 , 
                   help='Days for which to simulate noise time-series' )
parser.add_option( '--Tseg' , action='store' , dest='Tseg' , type='float' , nargs=1 ,
                   help='Duration of stitched segments making up the output time-series' )
parser.add_option( '--stime' , action='store' , dest='stime' , type='float' , nargs=1 , default='1.0' ,
                   help='sampling time of time-series' )
parser.add_option( '--seed' , action='store' , dest='seed' , type='string' , nargs=1 , default='random' ,
                   help='seed for random number generator' )
parser.add_option( '--tditype' , action='store' , dest='tditype' , type='string' , nargs=1 ,
                   default='Michelson' , help="Type of TDI observable ('Michelson','Sagnac')" )
parser.add_option( '--tdigen' , action='store' , dest='tdigen' , type='string' , nargs=1 , default='G0' ,
                   help="Generation of TDI observable ('G0','G1','Gm','G2')" )
parser.add_option( '--whichtdis' , action='store' , dest='whichtdis' , type='string' , nargs=1 ,
                   default='1' , help="Ordinary ('123') or optimal ('AET') TDI observables?" )



( options , args ) = parser.parse_args()

if len( args ) < 1 :
    parser.error( 'You must specify TSDIR!' )
else :
    tsdir = args[ 0 ]

dayinsecs = 86400. 
Nvar = 3

N = int( dayinsecs / options.stime )
n = int( options.Tseg / options.stime )

if N % n == 0 :
    Nseg = int( N / n ) * 2
else :
    Nseg = ( int( N / n ) + 1 ) * 2

duration = n * options.stime

freqdict = mufls.get_freqs_from_duration_and_stime( options.stime , duration )

comatrix = get_comatrix( freqdict['f'] , options.tditype , options.tdigen , options.whichtdis )

for day in options.days :
    print 'Day %d' % day
    t0 = ( day - 1 )*dayinsecs
    t = t0 + options.stime * np.arange( N )
    Npd_before_today = (day - 1) * Nseg * 2 * Nvar * freqdict['Nf'] 

    s = 1 ; tails = [ None ]*Nvar ; TSs = [ [] for v in range( Nvar ) ]
    while s <= Nseg :
        print 'Segments (%d|%d)' % ( s , s+1 )
        Npd_before_segl = Npd_before_today + (s-1)*2*Nvar*freqdict['Nf']
        Npd_before_segr = Npd_before_today + s*2*Nvar*freqdict['Nf']
        tl , tsl = mufls.get_noise_freq_domain_CovarMatrix(comatrix,freqdict['df'],
                                                           t0,freqdict['parityN'],
                                                           options.seed,
                                                           Npd_before_segl )
        tr , tsr = mufls.get_noise_freq_domain_CovarMatrix(comatrix,freqdict['df'],
                                                           t0,freqdict['parityN'],
                                                           options.seed,
                                                           Npd_before_segr )

        for v in range( Nvar ) :
            print 'v = ' , v
            if s == 1 :
                print 'tails[v]' , tails[v]
            else :
                print 'tails[v][:3]' , tails[v][:3]
            ts , tail = mufls.window_and_join( tsl[v] , tsr[v] , tails[v] )
            TSs[ v ] +=  list( ts )  ; tails[ v ] = np.copy( tail )
            print 'tails[v][:3]' , tails[v][:3]
        s += 2

    TSs = np.array( TSs )[ : , :N ]
    
    tscale = { 'Cadence1':options.stime , 'Offset1':t0 }
    tsdict = { 't':AS.Coarsable( t , **tscale ) , 
               '1':AS.Coarsable( TSs[0] , **tscale ) ,
               '2':AS.Coarsable( TSs[1] , **tscale ) ,
               '3':AS.Coarsable( TSs[2] , **tscale ) }    
    if tsdir == '' :
        pass
    elif tsdir not in glob.glob( tsdir ) :
        os.system( 'mkdir -p %s' % tsdir )
    print 'saving time-series to disk...'
    file = open( tsdir+'/d%03d.pkl' % day , 'wb' )
    cpkl.dump( tsdict , file , -1 ) ; file.close()    


