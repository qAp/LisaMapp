#!/usr/bin/env python
import os
import sys
import glob
import time
import cPickle as cpkl
from optparse import OptionParser

parser = OptionParser( '%prog SETUP.pkl' )
parser.add_option( '--do_P' , action='store_true' , help = 'Calculate clean map P' )
parser.add_option( '--do_PP' , action='store_true' , help = 'Calculate daily clean maps, PPs' )
parser.add_option( '--do_network_G' , action='store_true' , help = 'Calculate G for networks' )
parser.add_option( '--do_network_X' , action='store_true' , help = 'Calculate X for networks' )
parser.add_option( '--do_sigma_avg' , action='store_true' , help = 'Calculate sky-average of sigma, the std of Plm' )
parser.add_option( '--do_sigmamap' , action='store_true' , help = 'Calculate the sigma map, the standard deviation of the clean map in pixel basis)' )
parser.add_option( '--do_stdP' , action='store_true' , help='Calculate the standard deviations of the clean map Plm' )
parser.add_option( '--do_stdPP' , action='store_true' , help='Calculate the standard deviations of the daily clean maps PPlms' )
parser.add_option( '--do_SNRmap' , action='store_true' , help='Calculate the SNR map' )
parser.add_option( '--do_temp_optimals' , action='store_true' , help='Run the x_temp_optimals.py executable' )
parser.add_option( '--plot_X' , action='store_true' , help = 'Plot the dirty map X.' )
parser.add_option( '--plot_P' , action='store_true' , help = 'Plot the dirty map P.' )
parser.add_option( '--plot_sigmamap' , action='store_true' , help = 'Plot sigma map' )
parser.add_option( '--plot_SNRmap' , action='store_true' , help='Plot SNR map' )
parser.add_option( '--plot_sigma_avg' , action='store_true' , help = 'Plot sky-average of sigma, the std of Plm, versus lmax' )
parser.add_option( '--plot_singular_values' , action='store_true' , help = "Plot Fisher matrix's singular values in descending order." )
parser.add_option( '--plot_stdP' , action='store_true' , help='Plot standard deviation of the clean map, Plm' )


( options , args ) = parser.parse_args()

if len( args ) < 1 :
    parser.error( "You must specify a SETUP file of parameters!" )

setupname = args[0]
file = open( setupname , 'rb' ) ; setup = cpkl.load( file ) ; file.close()

execdir = setup['execdir'] 
workdir = os.getcwd() + '/' + setup['workdir']





if options.do_network_G :
    print 'Calculating G for a network...'
    os.chdir( workdir ) ; os.system( 'cp %s .' % ( execdir + '/x_network_G.py' ) )
    for network in setup['network_G']['networks'] :
        print 'Network for GWslope = %d, between %s' % ( network['GWslope'] , network['IJs'] )
        netGpath = workdir + '/GW_slope_%d/%s/G/G.pkl' % ( network['GWslope'] , '_'.join( network['IJs'] ) ) 
        Gpaths = [ ( workdir + '/GW_slope_%d/%s/G/G.pkl' % ( network['GWslope'] , IJ ) ) for IJ in network['IJs'] ]
        os.system( './x_network_G.py %s %s' % ( ' '.join( Gpaths ) , netGpath ) )
    print 'done'
    os.chdir( workdir )

if options.do_network_X :
    print 'Calculating X for a network...'
    os.chdir( workdir ) ; os.system( 'cp %s .' % ( execdir + '/x_network_X.py' ) )
    for network in setup['network_X']['networks'] :
        print 'Network for GWslope = %d, between %s' % ( network['GWslope'] , network['IJs'] )
        netXpath = workdir + '/GW_slope_%d/%s/X/X.pkl' % ( network['GWslope'] , '_'.join( network['IJs'] ) ) 
        Xpaths = [ ( workdir + '/GW_slope_%d/%s/X/X.pkl' % ( network['GWslope'] , IJ ) ) for IJ in network['IJs'] ]
        os.system( './x_network_X.py %s %s' % ( ' '.join( Xpaths ) , netXpath ) ) ; print 'done'


    os.chdir( workdir )





if options.do_PP :
    os.chdir( workdir )
    print 'Calculating the daily clean map...'
    os.system( 'cp %s .' % ( execdir + '/x_P.py' ) )
    for slope in setup['PP']['GWslopes'] :
        for IJ in setup['PP']['IJs'] :
            print 'GWslope = %d , IJ = %s' % ( slope , IJ )
            commands = []
            for day in setup['PP']['days'] :
#                print 'Day %d' % day
                GGpath = workdir + '/GW_slope_%d/%s/GG/GG_d%03d.pkl' % ( slope , IJ , day )
                XXpath = workdir + '/GW_slope_%d/%s/XX/XX_d%03d.pkl' % ( slope , IJ , day )
                PPpath = workdir + '/GW_slope_%d/%s/PP/PP_d%03d_lmax_%d.pkl' % ( slope , IJ , day , setup['PP']['lmax'] )
                if GGpath not in glob.glob( GGpath ) :
                    continue
                if XXpath not in glob.glob( XXpath ) :
                    continue
                commands += [  './x_P.py --regMethod %d --regCutoff %s --lmax %d --mapnorm %f --nlon %d --nlat %d --N_keptSV %s %s %s\n'
                               % ( setup['PP']['regMethod'] , setup['PP']['regCutoff'] , setup['PP']['lmax'] , setup['PP']['mapnorm'] , setup['PP']['nlon'] , setup['PP']['nlat'] , GGpath , XXpath , PPpath ) ]
            jobname = 'do_PP_slope_%d_IJ_%s' % ( slope , IJ )
            file = open( '%s.sub' % jobname , 'w' )
            file.writelines( [ '#!/bin/bash\n' ,
                               '#PBS -N %s\n' % jobname ,
                               '#PBS -o %s.out\n' % jobname ,
                               '#PBS -j oe\n' ,
                               '#PBS -q compute\n' ,
                               '#PBS -l nodes=1:ppn=1\n' ,
                               '#PBS -l walltime=10:00:00\n' , '\n' ,
                               'cd $PBS_O_WORKDIR\n' , '\n' ] + commands +
                             [ "echo 'do_PP done'" ] ) ; file.close()
            file = open( '%s.out' % jobname , 'w' ) ; file.write( 'dummpy output file' ) ; file.close()
            print 'Submitting job for %s...' % jobname ; os.system( 'qsub %s.sub' % jobname ) ; print 'done'
    os.chdir( workdir )


            


if options.do_P :
    os.chdir( workdir )
    print 'Calculating the clean map...'
    os.system( 'cp %s .' % ( execdir + '/x_P.py' ) )
    for slope in setup['P']['GWslopes'] :
        for IJ in setup['P']['IJs'] :
            print 'GWslope = %d , IJ = %s' % ( slope , IJ )
            Gdir = workdir + '/GW_slope_%d/%s/G/' % ( slope , IJ ) ; Gpath = Gdir + '/G.pkl'
            Xdir = workdir + '/GW_slope_%d/%s/X/' % ( slope , IJ ) ; Xpath = Xdir + '/X.pkl'
            Pdir = workdir + '/GW_slope_%d/%s/P/' % ( slope , IJ ) ; Ppath = Pdir + '/P_lmax_%d.pkl' % setup['P']['lmax']

            jobname = 'do_P_slope_%d_IJ_%s' % ( slope , IJ )
            file = open( '%s.sub' % jobname , 'w' )
            file.writelines( [ '#!/bin/bash\n' ,
                               '#PBS -N %s\n' % jobname ,
                               '#PBS -o %s.out\n' % jobname ,
                               '#PBS -j oe\n' ,
                               '#PBS -q compute\n' ,
                               '#PBS -l nodes=1:ppn=1\n' ,
                               '#PBS -l walltime=10:00:00\n' , '\n' ,
                               'cd $PBS_O_WORKDIR\n' , '\n' ,
                               './x_P.py --regMethod %d --regCutoff %s --lmax %d --mapnorm %f --nlon %d --nlat %d --N_keptSV %s %s %s\n'
                               % ( setup['P']['regMethod'] , setup['P']['regCutoff'] , setup['P']['lmax'] , setup['P']['mapnorm'] ,
                         setup['P']['nlon'] , setup['P']['nlat'] , Gpath , Xpath , Ppath ) ,
                               "echo 'do_P done'" ] ) ; file.close()
            file = open( '%s.out' % jobname , 'w' ) ; file.write( 'dummy output file' ) ; file.close()
            print 'Submitting job for do_P...' ; os.system( 'qsub %s.sub' % jobname ) ; print 'done.'
#            os.system( './x_P.py --regMethod %d --regCutoff %s --lmax %d --mapnorm %f --nlon %d --nlat %d --N_keptSV %s %s %s' %
#                       ( setup['P']['regMethod'] , setup['P']['regCutoff'] , setup['P']['lmax'] , setup['P']['mapnorm'] ,
#                         setup['P']['nlon'] , setup['P']['nlat'] , Gpath , Xpath , Ppath ) ) ; print 'done'
    os.chdir( workdir )





if options.do_sigma_avg :
    print 'Calculating sky-average of standard deviation of Plm...'
    os.chdir( workdir )
    os.system( 'cp %s .' % ( execdir + '/x_sigma_avg.py' ) )
    for slope in setup['sigma_avg']['GWslopes'] :
        for IJ in setup['sigma_avg']['IJs'] :
            print ( 'H(f) = f^{%d}' % slope ) , ( '%s' % IJ )
            Gdir = workdir + '/GW_slope_%d/%s/G/' % ( slope , IJ )
            avgsigdir = workdir + '/GW_slope_%d/%s/sigma_avg/' % ( slope , IJ )
            avgsigpath = avgsigdir + '/sigma_avg.pkl'
            submitname = 'x_sigma_avg_slope_%d_IJ_%s.sub' % ( slope , IJ )
            file = open( submitname , 'w' )
            file.writelines( [ '#!/bin/bash\n' , '#PBS -N %s\n' % submitname , '#PBS -q compute\n' , '#PBS -j oe\n' ,
                               '#PBS -l nodes=2:ppn=1\n' , '#PBS -l walltime=3:00:00\n' , 'cd $PBS_O_WORKDIR\n' , '\n' ,
                               ( './x_sigma_avg.py --regMethod %d --regCutoff %s --lmax %d --nlon %d --nlat %d %s %s' ) %
                               ( setup['sigma_avg']['regMethod'] , setup['sigma_avg']['regCutoff'] , setup['sigma_avg']['lmax'] ,
                                 setup['sigma_avg']['nlon'] , setup['sigma_avg']['nlat'] , Gdir , avgsigdir ) ] ) ; file.close()
            print 'Submitting job' ; os.system( 'qsub %s' % submitname ) ; print 'done'


    os.chdir( workdir )





if options.do_sigmamap :
    print 'Calculating sigma map...'
    os.chdir( workdir ) ; os.system( 'cp %s .' % ( execdir + '/x_sigmamap.py' ) )
    for slope in setup['sigmamap']['GWslopes'] :
        for IJ in setup['sigmamap']['IJs'] :
            print 'GWslope = %d , IJ = %s' % ( slope , IJ )
            Gdir = workdir + '/GW_slope_%d/%s/G/' % ( slope , IJ ) ; Gpath = Gdir + '/G.pkl'
            sigmamapdir = workdir + '/GW_slope_%d/%s/sigmamap/' % ( slope , IJ ) ; sigmamappath = sigmamapdir + '/S.pkl'
            submitname = 'x_sigmamap_slope_%d_IJ_%s.sub' % ( slope , IJ )
            file = open( submitname , 'w' )
            file.writelines( [ '#!/bin/bash\n' , '#PBS -N %s\n' % submitname , '#PBS -q compute\n' , '#PBS -j oe\n' ,
                               '#PBS -l nodes=1:ppn=1\n' , '#PBS -l walltime=3:00:00\n' , 'cd $PBS_O_WORKDIR\n' , '\n' ,
                               './x_sigmamap.py --regMethod %d --regCutoff %s --nlon %d --nlat %d --lmax %d --mapnorm %s %s %s' %
                               ( setup['sigmamap']['regMethod'] , setup['sigmamap']['regCutoff'] , setup['sigmamap']['nlon'] , setup['sigmamap']['nlat'] ,
                                 setup['sigmamap']['lmax'] , setup['sigmamap']['mapnorm'] , Gpath , sigmamappath ) ] ) ; file.close()
            print 'Submitting job' ; os.system( 'qsub %s' % submitname ) ; print 'done'


    os.chdir( workdir )






if options.do_stdPP :
    os.chdir( workdir )
    print 'Calculating standard deviation of daily clean maps PPlm...'
    os.system( 'cp %s .' % ( execdir + '/x_stdP.py' ) )
    for slope in setup['stdPP']['GWslopes'] :
        for IJ in setup['stdPP']['IJs'] :
            command = []
            for day in setup['stdPP']['days'] :
                GGpath = workdir + '/GW_slope_%d/%s/GG/GG_d%03d.pkl' % ( slope , IJ , day )
                if GGpath not in glob.glob( GGpath ) :
                    continue
                if setup['stdPP']['signal limit'] == 'weak' :
                    stdPPpath = workdir + '/GW_slope_%d/%s/stdPP/stdPP_d%03d_lmax_%d.pkl' % ( slope , IJ , day , setup['stdPP']['lmax'] )
                    command += [ './x_stdP.py --regMethod %d --regCutoff %s --lmax %d %s %s\n' %
                               ( setup['stdPP']['regMethod'] , setup['stdPP']['regCutoff'] , setup['stdPP']['lmax'] , GGpath , stdPPpath ) ]
                elif setup['stdPP']['signal limit'] == 'strong' :
                    SSpath = workdir + '/GW_slope_%d/%s/SS/SS_d%03d.pkl' % ( slope , IJ , day )
                    if SSpath not in glob.glob( SSpath ) :
                        continue
                    stdPPpath = workdir + '/GW_slope_%d/%s/stdPP/stdPP_d%03d_lmax_%d_strong.pkl' % ( slope , IJ , day , setup['stdPP']['lmax'] )           
                    command += [ './x_stdP.py --regMethod %d --regCutoff %s --lmax %d --strong_signal --Spath %s %s %s\n' %
                               ( setup['stdPP']['regMethod'] , setup['stdPP']['regCutoff'] , setup['stdPP']['lmax'] , SSpath , GGpath , stdPPpath ) ]
                elif setup['stdPP']['signal limit'] == 'both' :
                    SSpath = workdir + '/GW_slope_%d/%s/SS/SS_d%03d.pkl' % ( slope , IJ , day )
                    if SSpath not in glob.glob( SSpath ) :
                        continue
                    stdPPpath = workdir + '/GW_slope_%d/%s/stdPP/stdPP_d%03d_lmax_%d.pkl' % ( slope , IJ , day , setup['stdPP']['lmax'] )
                    command_weak = [ './x_stdP.py --regMethod %d --regCutoff %s --lmax %d --Spath %s %s %s\n' %
                               ( setup['stdPP']['regMethod'] , setup['stdPP']['regCutoff'] , setup['stdPP']['lmax'] , SSpath , GGpath , stdPPpath ) ]
                    stdPPpath = workdir + '/GW_slope_%d/%s/stdPP/stdPP_d%03d_lmax_%d_strong.pkl' % ( slope , IJ , day , setup['stdPP']['lmax'] )           
                    command_strong = [ './x_stdP.py --regMethod %d --regCutoff %s --lmax %d --strong_signal --Spath %s %s %s\n' %
                               ( setup['stdPP']['regMethod'] , setup['stdPP']['regCutoff'] , setup['stdPP']['lmax'] , SSpath , GGpath , stdPPpath ) ]
                    command += ( command_weak + command_strong )
            jobname = 'do_stdPP_slope_%d_IJ_%s' % ( slope , IJ )
            file = open( '%s.sub' % jobname , 'w' )
            file.writelines( [ '#!/bin/bash\n' ,
                               '#PBS -N %s\n' % jobname ,
                               '#PBS -o %s.out\n' % jobname ,
                               '#PBS -j oe\n' ,
                               '#PBS -q compute\n' ,
                               '#PBS -l nodes=1:ppn=1\n' ,
                               '#PBS -l walltime=10:00:00\n' , '\n' ,
                               'cd $PBS_O_WORKDIR\n' , '\n' ] + command + [ "echo 'do_stdPP done'" ]
                             ) ; file.close()
            file = open( '%s.out' % jobname , 'w' ) ; file.write( 'dummy output file' ) ; file.close()
            print 'Submitting job for %s...' % jobname ; os.system( 'qsub %s.sub' % jobname ) ; print 'done.'
    os.chdir( workdir )






if options.do_stdP :
    os.chdir( workdir )
    print 'Calculating standard deviation of clean map Plm...'
    os.system( 'cp %s .' % ( execdir + '/x_stdP.py' ) )
    for slope in setup['stdP']['GWslopes'] :
        for IJ in setup['stdP']['IJs'] :
            print 'GWslope = %d , IJ = %s' % ( slope , IJ )
            Gpath = workdir + '/GW_slope_%d/%s/G/G.pkl' % ( slope , IJ )
            if Gpath not in glob.glob( Gpath ) :
                print 'G not found at %s. Skip...' % Gpath ; continue
            if setup['stdP']['signal limit'] == 'weak' :
                stdPpath = workdir + '/GW_slope_%d/%s/stdP/stdP_lmax_%d.pkl' % ( slope , IJ , setup['stdP']['lmax'] )
                command = [ './x_stdP.py --regMethod %d --regCutoff %s --lmax %d %s %s\n'
                            % ( setup['stdP']['regMethod'] , setup['stdP']['regCutoff'] , setup['stdP']['lmax'] , Gpath , stdPpath ) ]
            elif setup['stdP']['signal limit'] == 'strong' :
                stdPpath = workdir + '/GW_slope_%d/%s/stdP/stdP_lmax_%d_strong.pkl' % ( slope , IJ , setup['stdP']['lmax'] )
                command = [ './x_stdP.py --regMethod %d --regCutoff %s --lmax %d --strong_signal --Spath %s %s %s\n'
                            % ( setup['stdP']['regMethod'] , setup['stdP']['regCutoff'] , setup['stdP']['lmax'] , setup['stdP']['Spath'] , Gpath , stdPpath ) ]
            elif setup['stdP']['signal limit'] == 'both' :
                stdPpath_weak = workdir + '/GW_slope_%d/%s/stdP/stdP_lmax_%d.pkl' % ( slope , IJ , setup['stdP']['lmax'] )
                command_weak = [ './x_stdP.py --regMethod %d --regCutoff %s --lmax %d %s %s\n'
                                 % ( setup['stdP']['regMethod'] , setup['stdP']['regCutoff'] , setup['stdP']['lmax'] , Gpath , stdPpath_weak ) ]
                stdPpath_strong = workdir + '/GW_slope_%d/%s/stdP/stdP_lmax_%d_strong.pkl' % ( slope , IJ , setup['stdP']['lmax'] )
                command_strong = [ './x_stdP.py --regMethod %d --regCutoff %s --lmax %d --strong_signal --Spath %s %s %s\n'
                                   % ( setup['stdP']['regMethod'] , setup['stdP']['regCutoff'] , setup['stdP']['lmax'] , setup['stdP']['Spath'] , Gpath , stdPpath_strong ) ] 
                command = command_weak + command_strong
            jobname = 'do_stdP_slope_%d_IJ_%s' % ( slope , IJ )
            file = open( '%s.sub' % jobname , 'w' )
            file.writelines( [ '#!/bin/bash\n' ,
                               '#PBS -N %s\n' % jobname ,
                               '#PBS -o %s.out\n' % jobname ,
                               '#PBS -j oe\n' ,
                               '#PBS -q compute\n' ,
                               '#PBS -l nodes=1:ppn=1\n' ,
                               '#PBS -l walltime=10:00:00\n' , '\n' ,
                               'cd $PBS_O_WORKDIR\n' , '\n' ] + command + [ "echo 'do_stdP done'" ]
                             ) ; file.close()
            file = open( '%s.out' % jobname , 'w' ) ; file.write( 'dummy output file' ) ; file.close()
            print 'Submitting job for %s...' % jobname ; os.system( 'qsub %s.sub' % jobname ) ; print 'done.'
    os.chdir( workdir )





if options.do_SNRmap :
    print 'Calculating SNR map...'
    os.chdir( workdir ) ; os.system( 'cp %s .' % ( execdir + '/x_SNRmap.py' ) )
    for slope in setup['SNRmap']['GWslopes'] :
        for IJ in setup['SNRmap']['IJs'] :
            print 'GWslope = %d , IJ = %s' % ( slope , IJ )
            Xpath = workdir + '/GW_slope_%d/%s/X/X.pkl' % ( slope , IJ )
            Gpath = workdir + '/GW_slope_%d/%s/G/G.pkl' % ( slope , IJ )
            SNRmappath = workdir + '/GW_slope_%d/%s/SNRmap/SNRmap.pkl' % ( slope , IJ )
            submitname = 'x_SNRmap_slope_%d_IJ_%s.sub' % ( slope , IJ )
            file = open( submitname , 'w' )
            file.writelines( [ '#!/bin/bash\n' , '#PBS -N %s\n' % submitname , '#PBS -q compute\n' , '#PBS -j oe\n' ,
                               '#PBS -l nodes=1:ppn=1\n' , '#PBS -l walltime=3:00:00\n' , 'cd $PBS_O_WORKDIR\n' , '\n' ,
                               './x_SNRmap.py --regMethod %d --regCutoff %s --nlon %d --nlat %d --lmax %d --mapnorm %s %s %s %s' %
                               ( setup['SNRmap']['regMethod'] , setup['SNRmap']['regCutoff'] ,
                                 setup['SNRmap']['nlon'] , setup['SNRmap']['nlat'] , setup['SNRmap']['lmax'] ,
                                 setup['SNRmap']['mapnorm'] , Gpath , Xpath , SNRmappath ) ] ) ; file.close()
            print 'Submitting job' ; os.system( 'qsub %s' % submitname ) ; print 'done'


    os.chdir( workdir )


    


if options.plot_X :
    print 'Plotting the dirty map X...' 
    os.chdir( workdir ) ; os.system( 'cp %s .' % ( execdir + '/plot_X.py' ) )
    for slope in setup['plot_X']['GWslopes'] :
        for IJ in setup['plot_X']['IJs'] :
            print 'GWslope = %d , IJ = %s' % ( slope , IJ )
            Xdir = workdir + '/GW_slope_%d/%s/X/' % ( slope , IJ ) ; Xpath = Xdir + '/X.pkl'
            figdir = workdir + '/GW_slope_%d/%s/figures_X/' % ( slope , IJ ) ; figpath = figdir + '/X_real.png'
            os.system( './plot_X.py --lmax %d --mapnorm %s --nlon %d --nlat %d --plot_imag %s %s' %
                       ( setup['plot_X']['lmax'] , setup['plot_X']['mapnorm'] , setup['plot_X']['nlon'] , setup['plot_X']['nlat'] ,
                         Xpath , figpath ) )
            
    print 'done' ; os.chdir( workdir )

if options.plot_P :
    print 'Plotting the clean map P...' 
    os.chdir( workdir ) ; os.system( 'cp %s .' % ( execdir + '/plot_P.py' ) )
    for slope in setup['plot_P']['GWslopes'] :
        for IJ in setup['plot_P']['IJs'] :
            print 'GWslope = %d , IJ = %s' % ( slope , IJ )
            Pdir = workdir + '/GW_slope_%d/%s/P/' % ( slope , IJ ) ; Ppath = Pdir + '/P_lmax_%d.pkl' % setup['plot_P']['lmax'] 
            figdir = workdir + '/GW_slope_%d/%s/figures_P/' % ( slope , IJ ) ; figpath = figdir + '/P_real_lmax_%d.png' % setup['plot_P']['lmax'] 
            os.system( './plot_P.py --lmax %d --mapnorm %s --nlon %d --nlat %d --plot_imag %s %s' %
                       ( setup['plot_P']['lmax'] , setup['plot_P']['mapnorm'] , setup['plot_P']['nlon'] , setup['plot_P']['nlat'] ,
                         Ppath , figpath ) )

            
    print 'done' ;  os.chdir( workdir )




if options.plot_sigmamap :
    print 'Plotting sigmamap...'
    os.chdir( workdir ) ; os.system( 'cp %s .' % ( execdir + '/plot_sigmamap.py' ) )
    for slope in setup['plot_sigmamap']['GWslopes'] :
        for IJ in setup['plot_sigmamap']['IJs'] :
            print 'GWslope = %d , IJ = %s' % ( slope , IJ )
            sigmamappath = workdir + '/GW_slope_%d/%s/sigmamap/S.pkl' % ( slope , IJ )
            figpath = workdir + '/GW_slope_%d/%s/figures_sigmamap/S_real.png' % ( slope , IJ )
            os.system( './plot_sigmamap.py --mapnorm %s --lmax %d --nlon %d --nlat %d %s %s' %
                       ( setup['plot_sigmamap']['mapnorm'] , setup['plot_sigmamap']['lmax'] ,
                         setup['plot_sigmamap']['nlon'] , setup['plot_sigmamap']['nlat'] , sigmamappath , figpath ) )


    print 'done' ; os.chdir( workdir )




if options.plot_sigma_avg :
    print 'Plotting sky-average of standard deviation of Plm...'
    os.chdir( workdir ) ; os.system( 'cp %s .' % ( execdir + '/plot_sigma_avg.py' ) )
    for slope in setup['plot_sigma_avg']['GWslopes'] :
        for IJ in setup['plot_sigma_avg']['IJs'] :
            print ( 'H(f) = f^{%d}' % slope ) , ( '%s' % IJ )
            avgsigdir = workdir + '/GW_slope_%d/%s/sigma_avg/' % ( slope , IJ )    
            avgsigfigdir = workdir + '/GW_slope_%d/%s/figures_sigma_avg/' % ( slope , IJ )
            avgsigpath = avgsigdir + '/sigma_avg.pkl'
            if avgsigpath not in glob.glob( avgsigpath ) :
                print 'sigma_avg.pkl not in %s. Nothing to plot.' % avgsigdir
                continue

            print 'sigma_avg.pkl in %s. Plotting average sigma versus lmax...' % avgsigdir
            os.system( './plot_sigma_avg.py %s %s' % ( avgsigdir , avgsigfigdir ) )

    print 'done' ; os.chdir( workdir )




if options.plot_stdP :
    print 'Plotting standard deviations of the clean map multipole moments, Plm,...'
    os.chdir( workdir ) ; os.system( 'cp %s .' % ( execdir + '/plot_stdP.py' ) )
    for slope in setup['plot_stdP']['GWslopes'] :
        for IJ in setup['plot_stdP']['IJs'] :
            print 'GWslope = %d , IJ = %s' % ( slope , IJ )
            stdPpath = workdir + '/GW_slope_%d/%s/stdP/stdP.pkl' % ( slope , IJ )
            figpath = workdir + '/GW_slope_%d/%s/figures_stdP/stdP.png' % ( slope , IJ )
            os.system( './plot_stdP.py --lmax %d --norm %s %s %s' % ( setup['plot_stdP']['lmax'] , setup['plot_stdP']['norm'] , stdPpath , figpath ) )


    print 'done' ; os.chdir( workdir )



if options.plot_SNRmap :
    print 'Plotting SNR map ...'
    os.chdir( workdir ) ; os.system( 'cp %s .' % ( execdir + '/plot_SNRmap.py' ) )
    for slope in setup['plot_SNRmap']['GWslopes'] :
        for IJ in setup['plot_SNRmap']['IJs'] :
            print 'GWslope = %d , IJ = %s' % ( slope , IJ )
            SNRmappath = workdir + '/GW_slope_%d/%s/SNRmap/SNRmap.pkl' % ( slope , IJ )
            figpath = workdir + '/GW_slope_%d/%s/figures_SNRmap/SNRmap_real.png' % ( slope , IJ )
            os.system( './plot_SNRmap.py --mapnorm %s %s %s' % ( setup['plot_SNRmap']['mapnorm'] , SNRmappath , figpath ) )


    print 'done' ; os.chdir( workdir )




if options.plot_singular_values :
    print 'Plotting singular values of Fisher matrix...'
    os.chdir( workdir ) ; os.system( 'cp %s .' % ( execdir + '/plot_singular_values.py' ) )
    for slope in setup['plot_singular_values']['GWslopes'] :
        for IJ in setup['plot_singular_values']['IJs'] :
            print 'GWslope = %d , IJ = %s' % ( slope , IJ )
            fishdir = workdir + '/GW_slope_%d/%s/G/' % ( slope , IJ )
            figdir = workdir + '/GW_slope_%d/%s/figures_singular_values/' % ( slope , IJ )
            fishpath = fishdir + '/G.pkl'
            figpath = figdir + '/singular_values.png'  
            if fishpath not in glob.glob( fishpath ) :
                print 'G.pkl not in %s.  Nothing to plot.' % fishdir
                continue

            print 'G.pkl in %s.  Plotting singular values...' % fishdir
            os.system( './plot_singular_values.py --lmax %d %s %s' % ( setup['plot_singular_values']['lmax'] , fishpath , figpath ) )


    print 'done' ; os.chdir( workdir )




if options.do_temp_optimals :
    os.chdir( workdir )
    for slope in setup['temp_optimals']['GWslopes'] :
        for IJ in setup['temp_optimals']['IJs'] :
            jobname_do_PP = 'do_PP_slope_%d_IJ_%s' % ( slope , IJ )
            jobname_do_stdPP = 'do_stdPP_slope_%d_IJ_%s' % ( slope , IJ )
            jobname_x_S = 'x_S_slope_%d_IJ_%s' % ( slope , IJ )
            print 'Waiting for do_PP, do_stdPP and x_S to finish for lope_%d IJ_%s...' % ( slope , IJ )
            while True :
                files = [ open( '%s.out' % jobname , 'r' ) for jobname
                          in [ jobname_do_PP , jobname_do_stdPP , jobname_x_S ] ]
                dones = [ file.readlines()[-1] for file in files ]
                [ file.close() for file in files ]
                if dones == [ 'do_PP done\n' , 'do_stdPP done\n' , 'done\n' ] :
                    print 'done' ; break
                else :
                    time.sleep( 5 ) ; continue
            print 'Calculating optimal clean map and standard deviation from daily estimates...'                    
            os.system( 'cp %s .' % ( execdir + '/x_temp_optimals.py' ) )
            PPdir = workdir + '/GW_slope_%d/%s/PP' % ( slope , IJ )
            stdPPdir = workdir + '/GW_slope_%d/%s/stdPP' % ( slope , IJ )
            SSdir = workdir + '/GW_slope_%d/%s/SS' % ( slope , IJ )
            sumdir = workdir + '/GW_slope_%d/%s/optimals/summary_lmax_%d.pkl' % ( slope , IJ , setup['temp_optimals']['lmax'] )
            jobname = 'do_optimals_slope_%d_IJ_%s' % ( slope , IJ )
            file = open( '%s.sub' % jobname , 'w' )
            file.writelines( [ '#!/bin/bash\n' ,
                               '#PBS -N %s\n' % jobname ,
                               '#PBS -o %s.out\n' % jobname , '#PBS -j oe\n' ,
                               '#PBS -q compute\n' ,
                               '#PBS -l nodes=1:ppn=1\n' , '#PBS -l walltime=10:00:00\n' , '\n' , 
                               'cd $PBS_O_WORKDIR\n' , '\n' ,
                               ( './x_temp_optimals.py ' + '-d%d '*len(setup['temp_optimals']['days'])
                                 + '--lmax %d %s %s %s %s\n' ) %
                               tuple( setup['temp_optimals']['days'] + [ setup['temp_optimals']['lmax'] ,
                                          PPdir , stdPPdir , SSdir , sumdir ] ) ,
                               '\n' , "echo 'do_temp_optimals done'" ] ) ; file.close()
            file = open( '%s.out' % jobname , 'w' ) ; file.write( 'dummy output file' ) ; file.close()
            print 'Submitting job for %s...' % jobname ; os.system( 'qsub %s.sub' % jobname ) ; print 'done.'
    os.chdir( workdir )


