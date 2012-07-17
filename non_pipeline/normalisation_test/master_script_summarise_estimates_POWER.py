#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import myLISAmodule as mlisar
import scipy.integrate as  integrate
#import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab


reali_nums = range( 1 , 25+1 )


"Expectation values (as calculated by expectation_values/script.py)"
P00_ana = 1 / 90.
sigP00_ana = 0.00052552385787
sigP00w_ana = 0.000394832846923
sigP00opt_ana = 2.75828491643e-05
sigP00optw_ana = 2.07233500415e-05


"Load up the theoretical estimates from analysis"
results = []
for r in reali_nums :
    print 'Working on realisation %d' % r
    workdir = 'analysis_seed_%d_stime_0.5_hanning_ORFsim/' % r

    Ppath = workdir + '/GW_slope_0/AE/P/P_lmax_0.pkl'
    stdPpath = workdir + '/GW_slope_0/AE/stdP/stdP_lmax_0_strong.pkl'
    stdPwpath = workdir + '/GW_slope_0/AE/stdP/stdP_lmax_0.pkl'
    sumpath = workdir + '/GW_slope_0/AE/optimals/summary_lmax_0.pkl'
    
    file = open( Ppath , 'rb' ) ; P = cpkl.load( file ) ; file.close()
    file = open( stdPpath , 'rb' ) ; stdP = cpkl.load( file ) ; file.close()
    file = open( stdPwpath , 'rb' ) ; stdPw = cpkl.load( file ) ; file.close()
    file = open( sumpath , 'rb' ) ; summary = cpkl.load( file ) ; file.close()

    P00s = summary['PPs'] 
    sigP00s = summary['stdPPs'] 
    sigP00sw = summary['stdPPws'] 
    P00opt = summary['P_opt'] 
    sigP00opt = summary['stdP_opt'] 
    sigP00optw = summary['stdPw_opt'] 

    results += [( r , np.std( P00s ) , np.mean( sigP00s ) , np.mean( sigP00sw ) ,
                  P00opt , ( P00opt - P00_ana ) , sigP00opt , sigP00optw ,
                  sigP00opt_ana )]
    
    print 'P00 estimates from the first 10 measurements(days)'
    print P00s[ :10 ]
    print 'theoretical sigmas from the first 10 measurements(days)'
    print sigP00s[ :10 ]
    print ''
    print 'empirical sigma |' , np.std( P00s )
    print 'sample mean of theoretical sigmas |' , np.mean( sigP00s )
    print 'sample mean of theoretical sigmas(weak signal limit) |' , np.mean( sigP00sw )
    print ''
    print 'optimal P00 estimate |' , P00opt
    print 'optimal P00 estimate - true P00 |' , P00opt - P00_ana
    print 'theoretical optimal sigma |' , sigP00opt
    print 'theoretical optimal sigma (weak signal limit) |' , sigP00optw
    print 'analytical optimal sigma |' , sigP00opt_ana
    print ''

lines = [ ( '%20d' + 8*'%20.10f' + '\n' ) % result for result in results ]
file = open( 'test_summary.dat' , 'w' )
print >> file , ( ( 9*'%20s' + '\n' ) % ( 'realisation' , 'S_emp' , 's.m. of S_the' ,
                                          's.m. of S_the(w)' , 'O_opt' ,
                                          'O_opt - O_true' , 'S_opt_theo' ,
                                          'S_opt_theo(w)' , 'S_opt_ana' ) )
file.writelines( lines )
file.close()
    



