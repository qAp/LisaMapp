#!/usr/bin/env python
import os
import sys
import numpy as np
import cPickle as cpkl
import glob 
import AnisotropySearch as AS
import optparse 

"""
This simulate the time-series for A,E and T for the source and days described in SETUP.pkl
"""

usage = """
%prog PPATH ORFDIR TSDIR\n
PPATH --- path to the skymap whose signal is to be simulated.\n
ORFDIR --- directory to the overlap reductions' multipole moments(effective the GW response to the skymap).\n
TSDIR --- path to the file in which to save the signal time-series.\n
"""
parser = optparse.OptionParser( usage = usage )

parser.add_option( '--GWslope' , action='store' , dest='GWslope' , type='int' , default=0 , nargs=1 ,
                   help='Gravitational spectral slope' )

parser.add_option( '--tditype' , action='store' , dest='tditype' , type='string' , default='Michelson' , nargs=1 ,
                   help='Type of TDI observable' )

parser.add_option( '--tdigen' , action='store' , dest='tdigen' , type='string' , default='G0' , nargs=1 ,
                   help='Generation of TDI observable' )

parser.add_option( '--whichtdi' , action='store' , dest='whichtdi' , type='choice' , default='optimal' , choices=('ordinary','optimal') , nargs=1 ,
                   help="Which TDI observables to simulate: 'ordinary' or 'optimal'" )

parser.add_option( '--lmax' , action='store' , dest='lmax' , type='int' , default=15 , nargs=1 , 
                   help='l up to which to convolve P and ORF' )

parser.add_option( '--stime' , action='store' , dest='stime' , type='float' , default=1.0 , nargs=1 ,
                   help='Sample time of the time-series' )

parser.add_option( '--f0' , action='store' , dest='f0' , type='float' , default=0.0001 , nargs=1 , 
                   help='Lower frequency of ORF' )

parser.add_option( '--df' , action='store' , dest='df' , type='float' , default=0.0001 , nargs=1 ,
                   help='Frequency resolution of ORF' )

parser.add_option( '--Nf' , action='store' , dest='Nf' , type='int' , default=4999 , nargs=1 ,
                   help='Number of frequencies of ORF' )

parser.add_option( '-d' , '--day' , action='append' , dest='days' , type='int' , nargs=1 ,
                   help='Days to simulate signals for' )

parser.add_option( '--seed' , action='store' , dest='seed' , type='string' , nargs=1 , default='None' ,
                   help='Seed for the random number generator.' )

( options , args ) = parser.parse_args()

if len( args ) < 3 :
    parser.error( "You must specify PPATH, ORFDIR AND TSDIR!  Type './x_P_to_TS.py -h' for help " )

Ppath , orfdir , tsdir = args[ :3 ]

#fNyq = 1. / ( 2 * options.stime )
#if ( ( options.Nf + 1 ) * options.df ) != fNyq :
#    raise InputError , 'The Nyquist frequency of ORF does not match the desired Nyquist frequency.  Cannot simulate.'

if options.whichtdi == 'ordinary' :
    IJlist = [ '11' , '12' , '13' , '22' , '23' , '33' ]
elif options.whichtdi == 'optimal' :
    IJlist = [ 'AA' , 'AE' , 'AT' , 'EE' , 'ET' , 'TT' ]

XYlist = [ '11' , '12' , '13' , '22' , '23' , '33' ]
slist  = [ 's1' , 's2' , 's3' ]
duration = 86400.

file = open( Ppath , 'rb' ) ; skymap = cpkl.load( file ) ; file.close

for day in options.days :

    print "" ; print "==== Day %d ====" % day ; print ""

    T0 = ( day - 1 ) * duration

    orfpaths = [ orfdir + 'tdiI_%s_tdiJ_%s_lmax_%d_f0_%f_df_%f_Nf_%d_g00_9/data/orf_d%03d.pkl' 
                 % ( IJ[0] , IJ[1] , options.lmax , options.f0 , options.df , options.Nf , day ) for IJ in IJlist ]
    orfs     = [ AS.OrfMultipleMoments( orfpath ) for orfpath in orfpaths ]

    cspecs = [ AS.Convolve( orf , skymap , options.GWslope ) for orf in orfs ]
    cspec_dict = dict( zip( XYlist , cspecs ) )

    """>>>>> Sort the covariance matrix elements into PSDs and CSDs and save them to disk """
    psddict = { 'f':orf.f , 'AA':cspec_dict['11'] , 'EE':cspec_dict['22'] , 'TT':cspec_dict['33'] }
    csddict = { 'f':orf.f , 'AE':cspec_dict['12'] , 'AT':cspec_dict['13'] , 'ET':cspec_dict['23'] }
    psddir , csddir = 'psd_injected/' , 'csd_injected/'
    if psddir not in glob.glob( psddir ) :
        os.system( 'mkdir %s' % psddir )
    file = open( psddir + 'd%03d.pkl' % day , 'wb' ) ; cpkl.dump( psddict , file , -1 ) ; file.close()
    if csddir not in glob.glob( csddir ) :
        os.system( 'mkdir %s' % csddir )
    file = open( csddir + 'd%03d.pkl' % day , 'wb' ) ; cpkl.dump( csddict , file , -1 ) ; file.close()
    """ <<<<< """

    """ Get from orfs the frequency resolution, in order to get Tseg """
    options.df   = orfs[0].f.Cadence1
    Tseg =  1. / options.df

    if ( duration / Tseg ) % 1 == 0 :
        Nseg = 2 * int( duration / Tseg )
    else:
        Nseg = 2 * ( int(  duration / Tseg ) + 1 )
        
    Ajoins , Ejoins , Tjoins = [] , [] , []
    leftover = None
        
    for i , seg in enumerate( range( 0 , Nseg , 2 ) ) :
        t0l = T0 + Tseg/2 * seg
        t0r = T0 + Tseg/2 * ( seg + 1 )
        print ""
        print "~~~~ Segment pair number: %d ~~~~" % ( i + 1 )
        print "Initial times:" , t0l , ',' , t0r
        print ""
        if options.seed == 'None' :
            shortftl = AS.CSpectra_to_ShortTermFT( cspec_dict , seed=None ) ; shortftr = AS.CSpectra_to_ShortTermFT( cspec_dict , seed=None )
        else :
            seedl = int( options.seed ) + ( day - 1 )*Nseg + seg ; seedr = int( options.seed ) + ( day - 1 )*Nseg + ( seg + 1 )
            print 'seeds: left, right' , seedl , seedr
            shortftl = AS.CSpectra_to_ShortTermFT( cspec_dict , seed=seedl ) ; shortftr = AS.CSpectra_to_ShortTermFT( cspec_dict , seed=seedr )
        stsl = AS.InverseFT( shortftl , Neven=True , tOffset=t0l ) ; stsr = AS.InverseFT( shortftr , Neven=True , tOffset=t0r )
        if i == 0 :
            Offset  = stsl.s1.Offset1 ; Cadence = stsl.s1.Cadence1
        joinrems = [ AS.WindowAndJoin( getattr( stsl , s ) , getattr( stsr , s ) , leftover = getattr( leftover , s , None ) ) for s in slist ]
        joins = [ joinrem[0] for joinrem in joinrems ]  
        rems  = [ joinrem[1] for joinrem in joinrems ]
        rem_dict = dict( zip( slist , rems ) )
        leftover = AS.sTimeSeries( **rem_dict )
        oldandnews = zip( [ Ajoins , Ejoins , Tjoins ] , joins )
        [ oldandnew[0].append( oldandnew[1].data ) for oldandnew in oldandnews ]
        
    tsdatas = [ np.concatenate( tuple( Ijoins ) ) for Ijoins in [ Ajoins , Ejoins , Tjoins ]  ]

    numberofsamplesinaday = int( duration / Cadence )
    
    onedaytsdatas = [ tsdata[ : numberofsamplesinaday ] for tsdata in tsdatas ]
    
    tscoarsables = [ AS.Coarsable( onedaytsdata , Offset1=Offset , Cadence1=Cadence ) for onedaytsdata in onedaytsdatas ]
    
    tsdict = dict( zip( [ '1' , '2' , '3' ] , tscoarsables ) )
    tsdict['t'] = AS.Coarsable( Offset + Cadence * np.arange( numberofsamplesinaday ) , Offset1=Offset , Cadence1=Cadence )

    if tsdir not in glob.glob( tsdir ) :
        os.system( 'mkdir -p %s' % tsdir )
    file = open( tsdir + '/d%03d.pkl' % day , 'wb' ) ; cpkl.dump( tsdict , file , -1 ) ; file.close()


file = open( tsdir + '/x_simulate_signal_FIN.pkl' , 'rb' ) ; FIN = cpkl.load( file ) ; file.close()
FIN += [ 1 ] ; file = open( tsdir + '/x_simulate_signal_FIN.pkl' , 'wb' ) ; cpkl.dump( FIN , file , -1 ) ; file.close()


                                                            

            
