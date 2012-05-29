#!/usr/bin/env python
import os
import sys
import glob
import time
import cPickle as cpkl
from optparse import OptionParser
import subprocess

parser = OptionParser( 'usage:x_analyse.py SETUP.pkl' )
parser.add_option( '--do_psd' , action='store_true' , help='Estimate PSDs' )
parser.add_option( '--do_csd' , action='store_true' , help='Estimate CSDs' )
parser.add_option( '--do_avgpsd' , action='store_true' , help='Calculate average PSDs' )
parser.add_option( '--do_X' , action='store_true' , help='Esimtate X, the dirty map' )
parser.add_option( '--do_G' , action='store_true' , help='Estimate G, the Fisher Matrix' )
parser.add_option( '--do_S' , action='store_true' , help='Estimate S, the strong signal bias matrix' )
parser.add_option( '--plot_ts' , action='store_true' , help='Plot time-series' )


( options , args ) = parser.parse_args()


if len( args ) < 1 :
    parser.error( "You must specify a SETUP file of parameters!" )

setupname = args[0]
file = open( setupname , 'rb' ) ; setup = cpkl.load( file ) ; file.close()

execdir = setup['execdir'] 
tsdir = setup['tsdir'] 
workdir = os.getcwd() + '/' + setup['workdir']

scale_ts = setup['scale_ts']


if workdir not in glob.glob( workdir ) :
    os.system( 'mkdir -p %s' % workdir )


if options.do_psd :
    print 'Submitting job for PSD estimation...'
    os.chdir( workdir )
    submitname = 'x_pklTStoPSD'
    os.system( 'cp %s .' % ( execdir + '/' + submitname + '.py' ) )
    psddir = workdir + '/psd/'
    days = setup['psd']['days']
    segduration = setup['psd']['segduration']
    file = open( submitname + '.sub' , 'w' )
    file.writelines( [ '#!/bin/bash\n' ,
                       '#PBS -N %s\n' % submitname ,
                       '#PBS -o %s.out\n' % submitname ,
                       '#PBS -j oe\n' ,
                       '#PBS -q compute\n' ,
                       '#PBS -l nodes=1:ppn=1\n' ,
                       '#PBS -l walltime=10:00:00\n' , '\n' ,
                       'cd $PBS_O_WORKDIR\n' , '\n' ,
                       ( './'+submitname+'.py ' + '-d%d '*len(days) + '--segduration %f --scale_ts %f %s %s\n' )
                       % tuple( days + [ segduration , scale_ts , tsdir , psddir ] ) ,
                       '\n' ,
                       "echo 'psd done'" ] )
    file.close()  
    p = subprocess.Popen( 'qsub %s.sub' % submitname , shell=True ,
                          stdout=subprocess.PIPE )
    psd_jobid = p.communicate()[0].rstrip()
    file = open( submitname+'_jobid' , 'w' )
    print >> file , psd_jobid
    file.close()
    print 'done'
    os.chdir( workdir )


if options.do_csd :
    print 'Submitting job for CSD estimation...'
    os.chdir( workdir )
    submitname = 'x_pklTStoCSD'
    os.system( 'cp %s .' % ( execdir + '/' + submitname + '.py' ) )
    csddir = workdir + '/csd/'
    days = setup['csd']['days']
    segduration = setup['csd']['segduration']
    file = open( submitname + '.sub' , 'w' )
    file.writelines( [ '#!/bin/bash\n' ,
                       '#PBS -N %s\n' % submitname ,
                       '#PBS -o %s.out\n' % submitname ,
                       '#PBS -j oe\n' ,
                       '#PBS -q compute\n' ,
                       '#PBS -l nodes=1:ppn=1\n' ,
                       '#PBS -l walltime=10:00:00\n' , '\n' ,
                       'cd $PBS_O_WORKDIR\n' , '\n' ,
                       ( './'+submitname+'.py ' + '-d%d '*len(days) + '--segduration %f --scale_ts %f %s %s\n' )
                       % tuple( days + [ segduration , scale_ts , tsdir , psddir ] ) ,
                       '\n' ,
                       "echo 'psd done'" ] )
    file.close()  
    p = subprocess.Popen( 'qsub %s.sub' % submitname , shell=True ,
                          stdout=subprocess.PIPE )
    csd_jobid = p.communicate()[0].rstrip()
    file = open( submitname+'_jobid' , 'w' )
    print >> file , csd_jobid
    file.close()
    print 'done'
    os.chdir( workdir )


if options.do_avgpsd :
    print 'Submitting job for averaging PSDs...'
    os.chdir( workdir )
    submitname = 'x_PSDtoAvgPSD'
    os.system( 'cp %s .' % ( execdir + '/' + submitname + '.py' ) )    
    psddir = workdir + '/psd/'
    avgpsddir = workdir + '/avgpsd/'
    days = setup['avgpsd']['days']
    file = open( submitname+'.sub' , 'w' )
    file.writelines( [ '#!/bin/bash\n' ,
                       '#PBS -N %s\n' % submitname ,
                       '#PBS -o %s.out\n' % submitname ,
                       '#PBS -j oe\n' ,
                       '#PBS -q compute\n' ,
                       '#PBS -l nodes=1:ppn=1\n' ,
                       '#PBS -l walltime=10:00:00\n' , '\n' ,
                       'cd $PBS_O_WORKDIR\n' , '\n' ,
                       ('./x_PSDtoAvgPSD.py ' + '-d%d '*len(days) + '%s %s\n' )
                       % tuple( days + [ psddir , avgpsddir ] ) , '\n' ,
                       "echo 'avgpsd done'" ] )
    file.close()
    if options.do_psd :
        afterok = psd_jobid
    else :
        afterok = ''
    p = subprocess.Popen( 'qsub -W depend=afterok:%s %s.sub' % ( afterok , submitname ) ,
                          shell=True , stdout=subprocess.PIPE )
    avgpsd_jobid = p.communicate()[0].rstrip()
    file = open( submitname+'_jobid' , 'w' )
    print >> file , avgpsd_jobid
    file.close()
    print 'done'
    os.chdir( workdir )


if options.do_X :
    print 'Submitting job/jobs for estimating the dirty map X^{IJ}_{\alpha}...'
    os.chdir( workdir )
    submitname = 'x_pklX'
    os.system( 'cp %s .' % ( execdir + '/' + submitname + '.py' ) )
    if options.do_avgpsd :
        afterok = avgpsd_jobid
    else :
        afterok = ''
    psddir = workdir + '/avgpsd/'
    orfdir = setup['X']['orfdir']
    days = setup['X']['days']
    for slope in setup['X']['GWslopes'] :
        for IJ in setup['X']['IJs'] :
            print ( 'GW spectral slope = %d' % slope ) , ( 'IJ = %s' % IJ )
            orfIJdir = orfdir + '/%s' % IJ
            Xdir = workdir + '/GW_slope_%d/%s/X/' % ( slope , IJ )
            cIJdir = workdir + '/GW_slope_%d/%s/cIJ/' % ( slope , IJ )
            jobname = submitname + '_slope_%d_IJ_%s' % ( slope , IJ )
            file = open( jobname+'.sub' , 'w' )
            file.writelines( [ '#!/bin/bash\n' ,
                               '#PBS -N %s\n' % jobname ,
                               '#PBS -o %s.out\n' % jobname ,
                               '#PBS -q compute\n' ,
                               '#PBS -j oe\n' ,
                               '#PBS -l nodes=1:ppn=1\n' ,
                               '#PBS -l walltime=5:00:00\n' ,
                               'cd $PBS_O_WORKDIR\n' ,
                               '\n' ,
                               ( './x_pklX.py ' + '-d%d '*len(days) +
                                 '--GWslope %d --scale_ts %f --flow %f --fhigh %f --lmax %d --window %s \
                                 %s %s %s %s %s\n')
                               % tuple( days + [ slope , setup['scale_ts'] , setup['X']['flow'] ,
                                                 setup['X']['fhigh'] , setup['X']['lmax'] ,
                                                 setup['X']['window'] ,
                                                 tsdir , orfIJdir , psddir , Xdir , cIJdir ] ) ,
                               'echo done' ] )
            file.close()
            p = subprocess.Popen( 'qsub -W depend=afterok:%s %s' % ( afterok , jobname+'.sub' ) ,
                                  shell=True , stdout=subprocess.PIPE )
            x_pklX_jobid = p.communicate()[0].rstrip() #this variable is overwritten for multiple IJ or GW slopes
            file = open( jobname+'_jobid' , 'w' )
            print >> file , x_pklX_jobid
            file.close()
    print 'done'
    os.chdir( workdir )



if options.do_G :
    print 'Submitting job/jobs for estimating Fisher matrix \Gamma^{IJ}_{\alpha\beta} ...'
    os.chdir( workdir )
    submitname = 'x_G'
    os.system( 'cp %s .' % ( execdir + '/' + submitname + '.py' ) )
    if options.do_avgpsd :
        afterok = avgpsd_jobid
    else :
        afterok = ''
    psddir = workdir + '/avgpsd/'
    orfdir = setup['G']['orfdir']
    days = setup['G']['days']
    for slope in setup['G']['GWslopes'] :
        for IJ in setup['G']['IJs'] :
            print ( 'GW spectral slope = %d' % slope ) , ( 'IJ = %s' % IJ )
            orfIJdir = orfdir + '/%s' % IJ
            Gdir = workdir + '/GW_slope_%d/%s/G/' % ( slope , IJ )
            cIIdir = workdir + '/GW_slope_%d/%s/cII/' % ( slope , IJ )            
            jobname = submitname+'_slope_%d_IJ_%s' % ( slope , IJ )
            file = open( jobname+'.sub' , 'w' )
            file.writelines( [ '#!/bin/bash\n' ,
                               '#PBS -N %s\n' % jobname ,
                               '#PBS -o %s.out\n' % jobname ,
                               '#PBS -q compute\n' ,
                               '#PBS -j oe\n' ,
                               '#PBS -l nodes=1:ppn=1\n' ,
                               '#PBS -l walltime=5:00:00\n' ,
                               'cd $PBS_O_WORKDIR\n' ,
                               '\n' ,
                               ( './x_G.py ' + '-d%d '*len(days) +
                                 '--GWslope %d --flow %f --fhigh %f --lmax %d --window %s %s %s %s %s %s\n')
                               % tuple( days + [ slope , setup['G']['flow'] , setup['G']['fhigh'] ,
                                                 setup['G']['lmax'] , setup['X']['window'] ,
                                                 tsdir , orfIJdir , psddir , Gdir , cIIdir ] ) ,
                               'echo done' ] )
            file.close()
            p = subprocess.Popen( 'qsub -W depend=afterok:%s %s' % ( afterok , jobname+'.sub' ) ,
                                  shell=True , stdout=subprocess.PIPE )
            x_G_jobid = p.communicate()[0].rstrip() #this variable is overwritten for multiple IJ or GW slopes
            file = open( jobname+'_jobid' , 'w' )
            print >> file , x_G_jobid
            file.close()
    print 'done'
    os.chdir( workdir )


if options.do_S :
    print 'Submitting job for estimating strong-signal-bias matrix \Psi^{IJ}_{\alpha\beta} ...'
    os.chdir( workdir )
    submitname = 'x_S'
    os.system( 'cp %s .' % ( execdir + '/'+submitname+'.py' ) )
    if options.do_avgpsd :
        afterok = avgpsd_jobid
    else :
        afterok = ''
    csddir = setup['S']['csddir']
    psddir = workdir + '/avgpsd/'
    orfdir = setup['S']['orfdir']
    days = setup['S']['days']
    for slope in setup['S']['GWslopes'] :
        for IJ in setup['S']['IJs'] :
            print ( 'GW spectral slope = %d' % slope ) , ( 'IJ = %s' % IJ )
            orfIJdir = orfdir + '/%s' % IJ
            Spath = workdir + '/GW_slope_%d/%s/S/S.pkl' % ( slope , IJ )            
            jobname = submitname+'_slope_%d_IJ_%s' % ( slope , IJ )
            file = open( jobname+'.sub' , 'w' )
            file.writelines( [ '#!/bin/bash\n' ,
                               '#PBS -N %s\n' % jobname ,
                               '#PBS -o %s.out\n' % jobname ,
                               '#PBS -q compute\n' ,
                               '#PBS -j oe\n' ,
                               '#PBS -l nodes=1:ppn=1\n' ,
                               '#PBS -l walltime=5:00:00\n' ,
                               'cd $PBS_O_WORKDIR\n' , '\n' ,
                               ( './x_S.py ' + '-d%d '*len(days) +
                                 '--GWslope %d --flow %f --fhigh %f --lmax %d --window %s %s %s %s %s %s\n' )
                               % tuple( days + [ slope , setup['S']['flow'] , setup['S']['fhigh'] ,
                                                 setup['S']['lmax'] , setup['X']['window'] ,
                                                 tsdir , csddir , orfIJdir , psddir , Spath ] ) ,
                               'echo done' ] )
            file.close()
            p = subprocess.Popen( 'qsub -W depend=afterok:%s %s' % ( afterok , jobname+'.sub' ) ,
                                  shell=True , stdout=subprocess.PIPE )
            x_S_jobid = p.communicate()[0].rstrip() #this variable is overwritten for multiple IJ or GW slopes
            file = open( jobname+'_jobid' , 'w' )
            print >> file , x_S_jobid
            file.close()
    print 'done'
    os.chdir( workdir )

            



        


