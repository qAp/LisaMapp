#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np


seeds = range( 1 , 25+1 )
days = range( 2 , 364+1 )
avgsddir = 'cIJ_seedsavg_daysavg'


Nseeds = len( seeds )
Ndays = len( days )
first_iter = True
for day in days :
    print 'Day %d' % day
    for seed in seeds :
        csdpath = 'analysis_seed_%d_stime_0.5_hanning_ORFsim/GW_slope_0/AE/cIJ/d%03d.pkl' % ( seed , day )
        file = open( csdpath , 'rb' )
        csddict = cpkl.load( file )
        file.close()
        if first_iter :
            P12 = csddict['AE'].data / (Nseeds*Ndays)
            P11 = csddict['AA'].data / (Nseeds*Ndays)
            P22 = csddict['EE'].data / (Nseeds*Ndays)
            g12 = csddict['gAE'].data[0]
            f = csddict['f'].data
        else :
            P12 += csddict['AE'].data / (Nseeds*Ndays)
            P11 += csddict['AA'].data / (Nseeds*Ndays)
            P22 += csddict['EE'].data / (Nseeds*Ndays)

        first_iter = False

#    sddict = {'f':f,'11':P11,'12':P12,'13':P13,'22':P22,'23':P23,'33':P33}
#    sdpath = avgsddir + '/d%03d.cpkl' % day
#    file = open( sd_path , 'wb' )
#    cpkl.dump( sddict , file , -1 )
#    file.close()
sdpath = avgsddir + '/avgsd.dat' 
if avgsddir not in glob.glob( avgsddir ) :
    os.system( 'mkdir -p %s' % avgsddir )
file = open( sdpath , 'w' )
print >> file , '\n'
print >> file , '#CSD'
print >> file , ('#'+'%19s'+8*'%20s') % ('freq[Hz]','Re{g12}','Im{g12}','Re{P12}','Im{P12}','Re{P11}',
                                         'Im{P11}','Re{P22}','Im{P22}')
for k in range( f.shape[0] ) :
    print >> file , (9*'%20.10e') % ( f[k] ,np.real(g12[k]),np.imag(g12[k]), np.real(P12[k]),np.imag(P12[k]),
                                      np.real(P11[k]),np.imag(P11[k]),
                                      np.real(P22[k]),np.imag(P22[k]) )
file.close()
    





    
        

