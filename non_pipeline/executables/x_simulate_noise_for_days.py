#!/usr/bin/env python
import os
import sys
import glob
import optparse
import cPickle as cpkl
import numpy as np


usage = """
%prog TSDIR\n
TSDIR --- name of directory in which to save the noise time-series files
Comment:\n
This executable needs to be used together with noise-simulating executable, like x_simulate_stationary_noise.py
"""

parser = optparse.OptionParser( usage = usage )

parser.add_option( '-d' , '--day' , action='append' , dest='days' , type='int' , nargs=1 , 
                   help='Days for which to simulate noise time-series' )

parser.add_option( '--stime' , action='store' , dest='stime' , type='string' , nargs=1 , default='1.' ,
                   help='Sampling time [s]' )

parser.add_option( '--seed' , action='store' , dest='seed' , type='string' , nargs=1 , default='random' ,
                   help='Seed for the random number generator in numpy.random' )

( options , args ) = parser.parse_args()


if len( args ) < 1 :
    parser.error( "You must specify TSDIR!  Type './x_simulate_noise_for_days.py -h' for help." )
else :
    tsdir = args[0]

#Number of time-series output by x_simulate_stationary_noise.py (<--Need to look at)
Nts = 3

dayinsecs = 86400.
#Work out Nf for getting the number of random numbers drawn
N = np.round( dayinsecs / float( options.stime ) )
if N % 2 == 0 :
    Nf = int( N/2 - 1 ) 
else :
    Nf = int( ( N-1 ) / 2 )

if options.days == None :
    print 'No days are selected.  Nothing to do.' ; sys.exit()

for day in options.days :
    print 'day %d' % day
    inittime = ( day - 1 )*dayinsecs
    if options.seed == 'random' :
        os.system( ( './x_simulate_stationary_noise.py --inittime %f --stime %s --duration %f --seed %s ' + tsdir + '/d%03d.pkl\n' ) % ( inittime , options.stime , dayinsecs , options.seed , day ) )
    else :
        #Number of random numbers drawn before this day, starting from day 1
        N_previous_draws = (day - 1)*Nf*2*Nts
        os.system( ( './x_simulate_stationary_noise.py --inittime %f --stime %s --duration %f --seed %s --N_previous_draws %d ' + tsdir + '/d%03d.pkl\n' ) % ( inittime , options.stime , dayinsecs , options.seed , N_previous_draws , day ) )



