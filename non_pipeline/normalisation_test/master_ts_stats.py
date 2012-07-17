#!/usr/bin/env python
import os
import sys
import glob
import time
import cPickle as cpkl
import numpy as np
import subprocess 


seeds = range( 1 , 25 + 1 )
ana_dir = 'time_series_statistical_properties'


workdir = os.getcwd()

first_seed = True
for seed in seeds :
    print '------------------------------------------------------------------'
    print '|                         Working on seed %d                        |' % seed
    print '------------------------------------------------------------------'

    print " Simulate signal (non-pipeline, freq-domain , from SpHs ) | sigsim "
    sigsim_dir = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/labs/sample_signals/LISA_signal/gIJ00day1_stime8_N10800_fordays_xxx_stime8_n500/'
    os.chdir( sigsim_dir )
    file = open( 'script_simulate_AETnoise_from_arbitrary_SpH_for_days_stitching.py' , 'r' ) ; lines = file.readlines() ; file.close()
    if seed == 'random' :
        lines[ 11 ] = "seed = %s\n" % seed
    else :
        lines[ 11 ] = "seed = '%d'\n" % seed
    file = open( 'script_simulate_AETnoise_from_arbitrary_SpH_for_days_stitching.py' , 'w' )
    file.writelines( lines ) ; file.close()
    os.system( 'python script_simulate_AETnoise_from_arbitrary_SpH_for_days_stitching.py' )
    os.chdir( workdir )

    print " Simulate noise (pipeline, freq-domain , stitch) | noisim "
    noisim_dir = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/labs/sample_noises/uncorr_white_noise_fd/stime_8_N_10800_for_days_xxx_stime_8_n_2000/'
    os.chdir( noisim_dir )
    file = open( 'script_simulate_noise_for_days_stitching.py' , 'r' ) ; lines = file.readlines() ; file.close()
    if seed == 'random' :
        lines[ 11 ] = "seed = %s\n" % seed
    else :
        lines[ 11 ] = "seed = '%d'\n" % ( seed + 1000000 )
    file = open( 'script_simulate_noise_for_days_stitching.py' , 'w' )
    file.writelines( lines ) ; file.close()
    os.system( 'python script_simulate_noise_for_days_stitching.py' )
    os.chdir( workdir )


    print " Merge signal and noise | merge "
    merge_dir = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/labs/sample_outputs/stime_8_N_10800_for_days/'
    os.chdir( merge_dir )
    file = open( sigsim_dir + '/signal_jobids' , 'r' )
    siglines = file.readlines() ; file.close()
    file = open( noisim_dir + '/noise_jobids' , 'r' )
    noilines = file.readlines() ; file.close()
    file = open( 'output_afteroks' , 'w' )
    file.writelines( siglines + noilines ) ; file.close()
    os.system( 'python script_setup_merge_tss.py' )
    os.system( './x_merge_tss.py setup_merge_tss.pkl' )
    os.chdir( workdir )

    print 'Waiting for merging of noise and signal...'
    while True :
        file = open( merge_dir + '/output_jobids' , 'r' )
        outlines = file.readlines() ; file.close()
        jobids = [ outline.rstrip() for outline in outlines ]
        jobxits = []
        for jobid in jobids :
            p = subprocess.Popen( 'qstat -f %s' % jobid ,
                                  shell=True , stdout=subprocess.PIPE ,
                                  stderr=subprocess.PIPE )
            stdout , stderr = p.communicate()
            if len(stdout)==0 and stderr==('qstat: Unknown Job Id %s' % jobid) :
                jobxits += [1]
            elif '\n' in stdout :
                qstats = [ line.strip() for line in stdout.split( '\n' ) ]
                if 'exit_status = 0' in qstats :
                    print 'found the exit_status=0 line!'
                    jobxits += [1]
                else :
                    jobxits += [0]
            else :
                jobxits += [0]
        print 'jobxits' , jobxits
        if jobxits.count(1) == len(jobids) :
            print 'Merging has finished.'
            break
        else :
            time.sleep( 20 )
            continue
    

    print " Calculate sample mean and covariance of y_{1} and y_{2} | ana_dir "
    if first_seed :
        stats = np.zeros( ( 6 , 365 ) )
        first_seed = False
    tspaths = glob.glob( merge_dir + '/data/d*.pkl' )
    days = [ int( os.path.basename( tspath ).rstrip('.pkl').lstrip('d').lstrip('0') ) for tspath in tspaths ]
    for d in range( len( days ) ) :
        file = open( tspaths[ d ] , 'rb' ) ; tsdict = cpkl.load( file ) ; file.close()
        s1 = tsdict[ '1' ].data ; s2 = tsdict[ '2' ].data
        mu_s1 = np.mean( s1 ) ; mu_s2 = np.mean( s2 )
        covmatrix = np.cov( np.array( [ s1 , s2 ] ) )
        sig_11  , sig_12 , sig_21 , sig_22 = covmatrix[0,0] , covmatrix[0,1] , covmatrix[1,0] , covmatrix[1,1]
        stats[ : , days[d]-1 ] += np.array( [ mu_s1 , mu_s2 , sig_11 , sig_12 , sig_21 , sig_22 ] )

    print " Clean up data-simulation directories from this realisation "
    print 'clean signal'
    os.chdir( sigsim_dir )
    os.system( 'rm -r data *.out *.sub' )
    os.chdir( workdir )
    print 'clean noise'
    os.chdir( noisim_dir )
    os.system( 'rm -r data *.out *.sub' )
    os.chdir( workdir )
    print 'clean output'
    os.chdir( merge_dir )
    os.system( 'rm -r data *.out *.sub' )
    os.chdir( workdir )

stats *= 1. / len( seeds )

if ana_dir not in glob.glob( ana_dir ) :
    os.system( 'mkdir -p %s' % ana_dir )
file = open( ana_dir + '/stats.pkl' , 'wb' )
cpkl.dump( stats , file , -1 ) ; file.close()


