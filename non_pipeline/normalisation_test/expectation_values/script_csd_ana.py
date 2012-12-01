#!/usr/bin/env python
import os
import sys
import glob
import numpy as np
import cPickle as cpkl
import synthlisa
import myLISAmodule as mlisar
import AnisotropySearch as AS


days = range( 1 , 365+1 )
csddir = 'csd_ana'
Ppath = ''
norm_factor = ''  #If 
orfdir = ''
GWSpectralSlope, H0 = 0, 1.0
IJs = [ 'AE' , 'AT' , 'ET' ]







file = open(Ppath, 'rb'); skymap = cpkl.load(file); file.close()

for day in days :
    csddict = {}
    for IJ in IJs :
        orf = AS.OrfMultipleMoments(orfdir + '/' + IJ + '/d%03d.pkl' % day)
        csddict[IJ] = AS.Convolve(orf, skymap, GWSpectralslope, H0)
        csddict['f'] = orf.f

    if csddir not in glob.glob( csddir ):
        os.system('mkdir %s' % csddir)
    file = open(csddir + '/d%03d.pkl' % day, 'wb'); cpkl.dump( csddict , file , -1 ); file.close()



    



