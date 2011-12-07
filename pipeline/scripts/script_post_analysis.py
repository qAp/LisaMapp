#!/usr/bin/env python
import os
import sys
import cPickle as cpkl


#Where are the executables for this? (Note we're on jackdt now.)
execdir = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/pipeline-running-scripts/new_executables/'

workdir = os.getcwd() + '/'

file = open( 'setup.pkl' , 'rb' )
setup = cpkl.load( file )
file.close()


for IJ in setup['postproc']['antennas'] :
    os.chdir( IJ )

    print "Post analysis: %s" % IJ
    
    os.system( 'cp %s .' % ( execdir+'x_postprocess.py' ) )

    os.system( './x_postprocess.py %s' % ( workdir+'setup.pkl' ) )

    os.chdir( workdir )
    
    
