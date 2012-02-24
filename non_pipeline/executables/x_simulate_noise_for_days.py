#!/usr/bin/env python
import os
import sys
import glob
import optparse

usage = """
%prog TSDIR\n
TSDIR --- name of directory in which to save the noise time-series files
"""

parser = optparse.OptionParser( usage = usage )

parser.add_option( '-d' , '--day' , action='append' , dest='days' , type='int' , nargs=1 , 
                   help='Days for which to simulate noise time-series' )

parser.add_option( '--stime' , action='store' , dest='stime' , type='string' , nargs=1 , default='1.' ,
                   help='Sampling time [s]' )

parser.add_option( '--seed' , action='store' , dest='seed' , type='string' , nargs=1 , default='None' ,
                   help='Seed for the random number generator in numpy.random' )

( options , args ) = parser.parse_args()


if len( args ) < 1 :
    parser.error( "You must specify TSDIR!  Type './x_simulate_noise_for_days.py -h' for help." )
else :
    tsdir = args[0]

dayinsecs = 86400.
#Work out Nf for getting the number of random numbers drawn
N = np.round( dayinsecs / stime )
if N % 2 == 0 :
    Nf = int( N/2 - 1 ) 
else :
    Nf = int( ( N-1 ) / 2 )

if options.days == None :
    print 'No days are selected.  Nothing to do.' ; sys.exit()

for day in options.days :
    inittime = ( day - 1 )*dayinsecs
    if options.seed = 'random' :
        os.system( ( './x_simulate_stationary_noise.py --inittime %f --stime %s --duration %f ' + tsdir + '/d%03d.pkl\n' ) % ( inittime , options.stime , dayinsecs , day ) )
    else :
        np.random.seed( options.seed )
        #draw preceeding random numbers first
        np.random.standard_normal( (day - 1)*Nf*2*Nts ) ;
        os.system( ( './x_simulate_stationary_noise.py --inittime %f --stime %s --duration %f ' + tsdir + '/d%03d.pkl\n' ) % ( inittime , options.stime , dayinsecs , day ) )


