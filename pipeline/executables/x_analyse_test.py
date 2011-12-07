#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
from optparse import OptionParser

parser = OptionParser( 'usage:x_analyse_test.py SETUP.pkl' )
parser.add_option( '--do_X' , action='store_true' , help='Esimtate X, the dirty map' )
parser.add_option( '--do_G' , action='store_true' , help='Estimate G, the Fisher Matrix' )
parser.add_option( '--do_S' , action='store_true' , help='Estimate Psi, the bias matrix of the covariance of the Plms.' )
parser.add_option( '--do_sigma_avg' , action='store_true' , help='Calculate sky-average of sigma, the std of Plm' )


( options , args ) = parser.parse_args()


if len( args ) < 1 :
    parser.error( "You must specify a SETUP file of parameters!" )

setupname = args[0]
file = open( setupname , 'rb' ) ; setup = cpkl.load( file ) ; file.close()

execdir = setup['execdir'] 
tsdir = setup['tsdir'] 
workdir = os.getcwd() + setup['workdir'] 

scale_ts = setup['scale_ts']

if workdir not in glob.glob( workdir ) :
    os.system( 'mkdir -p %s' % workdir )



""" Set the directories containing the analytical CSDs and the analytical PSDs manually """
psddir = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/runs/EccenctricInclined_eta0_0_xi0_0_sw_1_t0_0/signal/spherical_harmonics/Y_l0_m0_GWSlope_0/TDI/Michelson/G2/AET/simulation/stime_1.0_N_86400_for_days/psd_injected/'

csddir = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/runs/EccenctricInclined_eta0_0_xi0_0_sw_1_t0_0/signal/spherical_harmonics/Y_l0_m0_GWSlope_0/TDI/Michelson/G2/AET/analysis/stime_1.0_N_86400_for_days_for_N_realisations_all_estimated/r50/csd'




if options.do_X :
    print 'Calculating the dirty map...'
    os.chdir( workdir )
    os.system( 'cp %s .' % ( execdir + '/x_pklX_test.py' ) )
    orfdir = setup['X']['orfdir']
    days = setup['X']['days']
    for slope in setup['X']['GWslopes'] :
        for IJ in setup['X']['IJs'] :
            print ( 'H(f) = f^{%d}' % slope ) , ( '%s' % IJ )
            orfIJdir = orfdir + '/tdiI_Michelson_G2_%s_tdiJ_Michelson_G2_%s_lmax_20_f0_0.000100_df_0.000100_Nf_4999/data_nlon_120_nlat_61/' % tuple( IJ )
            Xdir = workdir + '/GW_slope_%d/%s/X/' % ( slope , IJ )
            cIJdir = workdir + '/GW_slope_%d/%s/cIJ/' % ( slope , IJ )
#            os.system(
#                ( './x_pklX_test.py ' + '-d%d '*len(days) + '--GWslope %d --scale_ts %f --flow %f --fhigh %f --lmax %d %s %s %s %s')
#                % tuple( days + [ slope , setup['scale_ts'] , setup['X']['flow'] , setup['X']['fhigh'] , setup['X']['lmax'] ,
#                                  csddir , orfIJdir , psddir , Xdir ] ) )
            submitname = 'x_pklX_test_slope_%d_IJ_%s.sub' % ( slope , IJ )
            file = open( submitname , 'w' )
            file.writelines( [ '#!/bin/bash\n' ,
                               '#PBS -N %s\n' %submitname ,
                               '#PBS -q compute\n' ,
                               '#PBS -j oe\n' ,
                               '#PBS -l nodes=1:ppn=1\n' ,
                               '#PBS -l walltime=5:00:00\n' ,
                               'cd $PBS_O_WORKDIR\n' ,
                               '\n' ,
                               ( './x_pklX_test.py ' + '-d%d '*len(days) + '--GWslope %d --scale_ts %f --flow %f --fhigh %f --lmax %d --window %s %s %s %s %s %s %s')
                               % tuple( days + [ slope , setup['scale_ts'] , setup['X']['flow'] , setup['X']['fhigh'] ,
                                                 setup['X']['lmax'] , setup['X']['window'] , 
                                                 tsdir , csddir , orfIJdir , psddir , Xdir , cIJdir ] ) ] ) ; file.close()
            print 'Submitting job' ; os.system( 'qsub %s' % submitname ) ; print 'done'
    os.chdir( workdir )



if options.do_G :
    print 'Calculating the Fisher matrix...'
    os.chdir( workdir )
    os.system( 'cp %s .' % ( execdir + '/x_G.py' ) )
    orfdir = setup['G']['orfdir']
    days = setup['G']['days']
    for slope in setup['G']['GWslopes'] :
        for IJ in setup['G']['IJs'] :
            print ( 'H(f) = f^{%d}' % slope ) , ( '%s' % IJ )
            orfIJdir = orfdir + '/tdiI_Michelson_G2_%s_tdiJ_Michelson_G2_%s_lmax_20_f0_0.000100_df_0.000100_Nf_4999/data_nlon_120_nlat_61/' % tuple( IJ )
            Gdir = workdir + '/GW_slope_%d/%s/G/' % ( slope , IJ )
#            os.system(
#                ( './x_G_test.py ' + '-d%d '*len(days) + '--GWslope %d --flow %f --fhigh %f --lmax %d %s %s %s %s')
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
                               ( './x_G.py ' + '-d%d '*len(days) + '--GWslope %d --flow %f --fhigh %f --lmax %d %s %s %s %s')
                               % tuple( days + [ slope , setup['G']['flow'] , setup['G']['fhigh'] , setup['G']['lmax'] ,
                                                 tsdir , orfIJdir , psddir , Gdir ] ) ] ) ; file.close()
            print 'Submitting job' ; os.system( 'qsub %s' % submitname ) ; print 'done'
    os.chdir( workdir )



if options.do_S :
    print 'Calculating the bias matrix of the covariance of the Plm in the strong-signal limit...'
    os.chdir( workdir )
    os.system( 'cp %s .' % ( execdir + '/x_S.py' ) )
    orfdir = setup['S']['orfdir']
    days = setup['S']['days']
    for slope in setup['S']['GWslopes'] :
        for IJ in setup['S']['IJs'] :
            print 'GWslope = %d , IJ = %s' % ( slope , IJ )
            orfIJdir = orfdir + '/tdiI_Michelson_G2_%s_tdiJ_Michelson_G2_%s_lmax_20_f0_0.000100_df_0.000100_Nf_4999/data_nlon_120_nlat_61/' % tuple( IJ )
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

    


            
if options.do_sigma_avg :
    print 'Calculating sky-average of standard deviation of Plm...'
    os.chdir( workdir )
    os.system( 'cp %s .' % ( execdir + 'x_sigma_avg.py' ) )
    for slope in setup['sigma_avg']['GWslopes'] :
        for IJ in setup['sigma_avg']['IJs'] :
            print ( 'H(f) = f^{%d}' % slope ) , ( '%s' % IJ )
            Gdir = workdir + 'GW_slope_%d/%s/G/' % ( slope , IJ )
            avgsigdir = workdir + 'GW_slope_%d/%s/sigma_avg/' % ( slope , IJ )

#            os.system( ( './x_sigma_avg.py --regMethod %d --regCutoff %f --lmax %d --nlon %d --nlat %d --plot %s %s' ) %
#                       ( setup['sigma_avg']['regMethod'] , setup['sigma_avg']['regCutoff'] , setup['sigma_avg']['lmax'] ,
#                         setup['sigma_avg']['nlon'] , setup['sigma_avg']['nlat'] , Gdir , avgsigdir ) )
            
            submitname = 'x_sigma_avg_slope_%d_IJ_%s.sub' % ( slope , IJ )
            file = open( submitname , 'w' )
            file.writelines( [ '#!/bin/bash\n' ,
                               '#PBS -N %s\n' % submitname ,
                               '#PBS -q compute\n' , 
                               '#PBS -j oe\n'
                               '#PBS -l nodes=1:ppn=1\n' ,
                               '#PBS -l walltime=1:00:00\n' ,
                               'cd $PBS_O_WORKDIR\n' ,
                               '\n' ,
                               ( './x_sigma_avg.py --regMethod %d --regCutoff %f --lmax %d --nlon %d --nlat %d %s %s' ) %
                               ( setup['sigma_avg']['regMethod'] , setup['sigma_avg']['regCutoff'] , setup['sigma_avg']['lmax'] ,
                                 setup['sigma_avg']['nlon'] , setup['sigma_avg']['nlat'] , Gdir , avgsigdir )
                                ] )
            file.close()
            print 'Submitting job'
            os.system( 'qsub %s' % submitname )
            print 'done'
    os.chdir( workdir )


        


