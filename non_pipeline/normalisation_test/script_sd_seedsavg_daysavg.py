#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np


seeds = range( 1 , 25+1 )
days = range( 1 , 365+1 )
avgsddir = 'sd_seedsavg_daysavg'


Nseeds = len( seeds )
Ndays = len( days )
first_iter = True
for day in days :
    print 'Day %d' % day
    for seed in seeds :
        psdpath = 'analysis_seed_%d_stime_0.5_hanning_ORFsim/psd/d%03d.pkl' % ( seed , day )
        file = open( psdpath , 'rb' )
        psddict = cpkl.load( file )
        file.close()
        if first_iter :
            P11 = psddict['AA'].data / (Nseeds*Ndays)
            P22 = psddict['EE'].data / (Nseeds*Ndays)
            P33 = psddict['TT'].data / (Nseeds*Ndays)
            f = psddict['f'].data
        else :
            P11 += psddict['AA'].data / (Nseeds*Ndays)
            P22 += psddict['EE'].data / (Nseeds*Ndays)
            P33 += psddict['TT'].data / (Nseeds*Ndays)

        csdpath = 'analysis_seed_%d_stime_0.5_hanning_ORFsim/csd/d%03d.pkl' % ( seed , day )
        file = open( csdpath , 'rb' )
        csddict = cpkl.load( file )
        file.close()
        if first_iter :
            P12 = csddict['AE'].data / (Nseeds*Ndays)
            P13 = csddict['AT'].data / (Nseeds*Ndays)
            P23 = csddict['ET'].data / (Nseeds*Ndays)
        else :
            P12 += csddict['AE'].data / (Nseeds*Ndays)
            P13 += csddict['AT'].data / (Nseeds*Ndays)
            P23 += csddict['ET'].data / (Nseeds*Ndays)

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
print >> file , '#PSD'
print >> file , ('#'+'%19s'+3*'%20s') % ('freq[Hz]','P11','P22','P33')
for k in range( f.shape[0] ) :
    print >> file , (4*'%20.10f') % (f[k],P11[k],P22[k],P33[k])
print >> file , '\n'
print >> file , '#CSD'
print >> file , ('#'+'%19s'+6*'%20s') % ('freq[Hz]','Re{P12}','Im{P12}','Re{P13}',
                                         'Im{P13}','Re{P23}','Im{P23}')
for k in range( f.shape[0] ) :
    print >> file , (7*'%20.10f') % ( f[k] , np.real(P12[k]),np.imag(P12[k]),
                                      np.real(P13[k]),np.imag(P13[k]),
                                      np.real(P23[k]),np.imag(P23[k]) )
file.close()
    





    
        

