#!/usr/bin/env python
import os
import sys
import glob
import re
import cPickle as cpkl
import numpy as np
import myUsefuls as mufls
import AnisotropySearch as AS
import optparse


def get_comatrix( f ) :
    """ Define the covariance matrix here!!! """
    Nf = f.shape[ 0 ] 
    p12 , p13 , p23 = np.zeros( (Nf,) ) , np.zeros( (Nf,) ) , np.zeros( (Nf,) )
    p11 , p22 , p33 = 1. * np.ones( (Nf,) ) , 1. * np.ones( (Nf,) ) , 1. * np.ones( (Nf,) )
    comatrix = np.array( [ [ p11 , p12 , p13 ] ,
                           [ np.conj(p12) , p22 , p23 ] ,
                           [ np.conj(p13) , np.conj(p23) , p33 ] ] )
    return comatrix



usage = """
%prog TSPATH\n
This simulates noise time-series from the covariance matrix defined above here.
TSPATH --- path to file in which to save the time-series
"""
parser = optparse.OptionParser( usage = usage )

parser.add_option( '--stime' , action='store' , dest='stime' , type='string' , nargs=1 , default='1.' ,
                   help='Sampling time [s]' )

parser.add_option( '--duration' , action='store' , dest='duration' , type='string' , nargs=1 , default='86400.' ,
                   help='Duration [s]' )

parser.add_option( '--inittime' , action='store' , dest='inittime' , type='string' , nargs=1 , default='0.' ,
                   help='Initial time [s]' )

parser.add_option( '--seed' , action='store' , dest='seed' , type='string' , nargs=1 , default='random' , help='Seed for random number generation' )

parser.add_option( '--N_previous_draws' , action='store' , dest='N_previous_draws' , type='string' , nargs=1 , default='0' , help='Number of draws to discard first' )


( options , args ) = parser.parse_args()
if len( args ) < 1 :
    parser.error( "You must specify a TSPATH!  Type './x_simulate_defined_stationary_noise.py' " )
else :
    tspath = args[ 0 ] ; tsdir = os.path.dirname( tspath )

stime , duration , inittime , N_previous_draws = float(options.stime) , float(options.duration) , float(options.inittime) , int( options.N_previous_draws )

N = np.round( duration / stime )
if N % 2 == 0 :
    Nf = int( N/2 - 1 ) ; parityN = 'Even'
else :
    Nf = int( ( N-1 ) / 2 ) ; parityN = 'Odd'

df = 1 / (N*stime)
f = df * np.arange( 1 , Nf + 1 )

comatrix = get_comatrix( f )

t , n = mufls.get_noise_freq_domain_CovarMatrix( comatrix , df , inittime , parityN , options.seed , N_previous_draws ) 

tscale = { 'Cadence1':stime , 'Offset1':inittime }
tsdict = { 't':AS.Coarsable( t , **tscale ) , 
           '1':AS.Coarsable( n[0] , **tscale ) , '2':AS.Coarsable( n[1] , **tscale ) , '3':AS.Coarsable( n[2] , **tscale ) }

if tsdir not in glob.glob( tsdir ) :
    os.system( 'mkdir %s' % tsdir )
print 'saving time-series to disk...'
file = open( tspath , 'wb' ) ; cpkl.dump( tsdict , file , -1 ) ; file.close()

