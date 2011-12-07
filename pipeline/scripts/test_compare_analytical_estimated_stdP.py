#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np


IJs = [ 'AE' , 'AT'  , 'ET' , 'AE_AT' , 'AE_ET' , 'AT_ET' , 'AE_AT_ET' ]
xdir = 'GW_slope_0_test1/'
edir = 'GW_slope_0_estimated2/'



for IJ in IJs :
    print 'Comparing standard deviations for %s' % IJ
    xpath = xdir + '%s/stdP/stdP.pkl' % IJ
    file = open( xpath , 'rb' ) ; xstdP = cpkl.load( file ) ; file.close()
    epath = edir + '%s/stdP/stdP.pkl' % IJ
    file = open( epath , 'rb' ) ; estdP = cpkl.load( file ) ; file.close()
    print 'expected:' , xstdP 
    print 'estimated:' , estdP * 1e-40

    
    





