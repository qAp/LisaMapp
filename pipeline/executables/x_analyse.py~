#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
from optparse import OptionParser

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
    print 'Estimating PSDs...'
    os.chdir( workdir )
    psddir = workdir + '/psd/'
    days = setup['psd']['days']
    segduration = setup['psd']['segduration']
    os.system( 'cp %s .' % ( execdir + 'x_pklTStoPSD.py' ) )
    os.system( ( './x_pklTStoPSD.py ' + '-d%d '*len(days) +
                 '--segduration %f --scale_ts %f %s %s' ) % tuple( days + [ segduration , scale_ts , tsdir , psddir ] ) )
    print 'done'
    os.chdir( workdir )


if options.do_csd :
    print 'Estimating CSDs...'
    os.chdir( workdir )
    csddir = workdir + '/csd/'
    days = setup['csd']['days']
    segduration = setup['csd']['segduration']
    os.system( 'cp %s .' % ( execdir + 'x_pklTStoCSD.py' ) )
    os.system(
        ( './x_pklTStoCSD.py ' + '-d%d '*len(days) + '--segduration %f --scale_ts %f %s %s' )
        % tuple( days + [ segduration , scale_ts , tsdir , csddir ] ) )
    print 'done'
    os.chdir( workdir )


if options.do_avgpsd :
    print 'Averaging PSDs...'
    os.chdir( workdir )
    psddir = workdir + '/psd/'
    avgpsddir = workdir + '/avgpsd/'
    days = setup['avgpsd']['days']
    os.system( 'cp %s .' % ( execdir + 'x_PSDtoAvgPSD.py' ) )
    os.system( ('./x_PSDtoAvgPSD.py ' + '-d%d '*len(days) + '%s %s' )
               % tuple( days + [ psddir , avgpsddir ] ) )
    print 'done'
    os.chdir( workdir )


if options.do_X :
    print 'Calculating the dirty map...'
    os.chdir( workdir )
    os.system( 'cp %s .' % ( execdir + '/x_pklX.py' ) )
    psddir = workdir + '/avgpsd/'
    orfdir = setup['X']['orfdir']
    days = setup['X']['days']
    for slope in setup['X']['GWslopes'] :
        for IJ in setup['X']['IJs'] :
            print ( 'H(f) = f^{%d}' % slope ) , ( '%s' % IJ )
            orfIJdir = orfdir + '/tdiI_%s_tdiJ_%s_lmax_0_f0_0.000000_df_0.000174_Nf_5761/data/' % tuple( IJ )
            Xdir = workdir + '/GW_slope_%d/%s/X/' % ( slope , IJ )
            cIJdir = workdir + '/GW_slope_%d/%s/cIJ/' % ( slope , IJ )
#            os.system(
#                ( './x_pklX.py ' + '-d%d '*len(days) + '--GWslope %d --scale_ts %f --flow %f --fhigh %f --lmax %d %s %s %s %s')
#                % tuple( days + [ slope , setup['scale_ts'] , setup['X']['flow'] , setup['X']['fhigh'] , setup['X']['lmax'] ,
#                                  tsdir , orfIJdir , psddir , Xdir ] ) )
            submitname = 'x_pklX_slope_%d_IJ_%s.sub' % ( slope , IJ )
            file = open( submitname , 'w' )
            file.writelines( [ '#!/bin/bash\n' ,
                               '#PBS -N %s\n' %submitname ,
                               '#PBS -q compute\n' ,
                               '#PBS -j oe\n' ,
                               '#PBS -l nodes=1:ppn=1\n' ,
                               '#PBS -l walltime=5:00:00\n' ,
                               'cd $PBS_O_WORKDIR\n' ,
                               '\n' ,
                               ( './x_pklX.py ' + '-d%d '*len(days) + '--GWslope %d --scale_ts %f --flow %f --fhigh %f --lmax %d --window %s %s %s %s %s %s')
                               % tuple( days + [ slope , setup['scale_ts'] , setup['X']['flow'] , setup['X']['fhigh'] ,
                                                 setup['X']['lmax'] , setup['X']['window'] ,
                                                 tsdir , orfIJdir , psddir , Xdir , cIJdir ] ) ] ) ; file.close()
            print 'Submitting job' ; os.system( 'qsub %s' % submitname ) ; print 'done'
    os.chdir( workdir )



if options.do_G :
    print 'Calculating the Fisher matrix...'
    os.chdir( workdir )
    os.system( 'cp %s .' % ( execdir + '/x_G.py' ) )
    psddir = workdir + '/avgpsd/'
    orfdir = setup['G']['orfdir']
    days = setup['G']['days']
    for slope in setup['G']['GWslopes'] :
        for IJ in setup['G']['IJs'] :
            print ( 'H(f) = f^{%d}' % slope ) , ( '%s' % IJ )
            orfIJdir = orfdir + '/tdiI_%s_tdiJ_%s_lmax_0_f0_0.000000_df_0.000174_Nf_5761/data/' % tuple( IJ )
            Gdir = workdir + '/GW_slope_%d/%s/G/' % ( slope , IJ )
            cIIdir = workdir + '/GW_slope_%d/%s/cII/' % ( slope , IJ )            
#            os.system(
#                ( './x_G.py ' + '-d%d '*len(days) + '--GWslope %d --flow %f --fhigh %f --lmax %d %s %s %s %s')
#                % tuple( days + [ slope , setup['G']['flow'] , setup['G']['fhigh'] , setup['G']['lmax'] ,
#                                  tsdir , orfIJdir , psddir , Gdir ] ) )
            submitname = 'x_G_slope_%d_IJ_%s.sub' % ( slope , IJ )
            file = open( submitname , 'w' )
            file.writelines( [ '#!/bin/bash\n' ,
                               '#PBS -N %s\n' %submitname ,
                               '#PBS -q compute\n' ,
                               '#PBS -j oe\n' ,
                               '#PBS -l nodes=1:ppn=1\n' ,
                               '#PBS -l walltime=5:00:00\n' ,
                               'cd $PBS_O_WORKDIR\n' ,
                               '\n' ,
                               ( './x_G.py ' + '-d%d '*len(days) + '--GWslope %d --flow %f --fhigh %f --lmax %d %s %s %s %s %s')
                               % tuple( days + [ slope , setup['G']['flow'] , setup['G']['fhigh'] , setup['G']['lmax'] ,
                                                 tsdir , orfIJdir , psddir , Gdir , cIIdir ] ) ] ) ; file.close()
            print 'Submitting job' ; os.system( 'qsub %s' % submitname ) ; print 'done'
    os.chdir( workdir )




if options.do_S :
    print 'Calculating the bias matrix of the covariance of the Plm in the strong-signal limit...'
    os.chdir( workdir )
    os.system( 'cp %s .' % ( execdir + '/x_S.py' ) )
    csddir = setup['S']['csddir']
    psddir = workdir + '/avgpsd/'
    orfdir = setup['S']['orfdir']
    days = setup['S']['days']
    for slope in setup['S']['GWslopes'] :
        for IJ in setup['S']['IJs'] :
            print 'GWslope = %d , IJ = %s' % ( slope , IJ )
            orfIJdir = orfdir + '/tdiI_%s_tdiJ_%s_lmax_0_f0_0.000000_df_0.000174_Nf_5761/data/' % tuple( IJ )
            Spath = workdir + '/GW_slope_%d/%s/S/S.pkl' % ( slope , IJ )            
#            os.system( ( './x_S.py ' + '-d%d '*len(days) + '--GWslope %d --flow %f --fhigh %f --lmax %d %s %s %s %s %s' )
#                       % tuple( days + [ slope , setup['S']['flow'] , setup['S']['fhigh'] , setup['S']['lmax'] ,
#                                         tsdir , csddir , orfIJdir , psddir , Spath ] ) )
            submitname = 'x_S_slope_%d_IJ_%s.sub' % ( slope , IJ )
            file = open( submitname , 'w' )
            file.writelines( [ '#!/bin/bash\n' , '#PBS -N %s\n' %submitname , '#PBS -q compute\n' ,
                               '#PBS -j oe\n' , '#PBS -l nodes=1:ppn=1\n' , '#PBS -l walltime=5:00:00\n' ,
                               'cd $PBS_O_WORKDIR\n' , '\n' ,
                               ( './x_S.py ' + '-d%d '*len(days) + '--GWslope %d --flow %f --fhigh %f --lmax %d %s %s %s %s %s' )
                               % tuple( days + [ slope , setup['S']['flow'] , setup['S']['fhigh'] , setup['S']['lmax'] ,
                                                 tsdir , csddir , orfIJdir , psddir , Spath ] ) ] ) ; file.close()
            print 'Submitting job' ; os.system( 'qsub %s' % submitname ) ; print 'done'


    os.chdir( workdir )

            


if options.plot_ts :
    print 'Plotting time-series...'
    os.chdir( workdir )
    os.system( 'cp %s .' % ( execdir + '/plot_ts.py' ) )
    for day in setup['plot_ts']['days'] :
        tspath = tsdir + '/d%03d.pkl' % day
        figpath = workdir + '/figures_ts/d%03d.png' % day
        os.system( './plot_ts.py --scale_ts %f %s %s' % ( scale_ts , tspath , figpath ) )
    print 'done'
    os.chdir( workdir )
        


