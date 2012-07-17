#!/usr/bin/env python
import os
import sys
import glob
import time
import cPickle as cpkl
import subprocess


seeds = range( 1 , 25+1 )



workdir = os.getcwd()

for seed in seeds :
    print '------------------------------------------------------------------'
    print '|                         Testing seed %d                        |' % seed
    print '------------------------------------------------------------------'

    print " Simulate signal (non-pipeline, freq-domain , from SpHs , stitch) | sigsim "
    sigsim_dir = workdir + '/signal'
    os.chdir( sigsim_dir )
    file = open( 'script_simulate_AETnoise_from_arbitrary_SpH_for_days_stitching.py' , 'r' )
    lines = file.readlines()
    file.close()
    if seed == 'random' :
        lines[ 11 ] = "seed = %s\n" % seed
    else :
        lines[ 11 ] = "seed = '%d'\n" % seed
    file = open( 'script_simulate_AETnoise_from_arbitrary_SpH_for_days_stitching.py' , 'w' )
    file.writelines( lines )
    file.close()
    os.system( 'python script_simulate_AETnoise_from_arbitrary_SpH_for_days_stitching.py' )
    os.chdir( workdir )
            

    print " Simulate noise (pipeline, freq-domain, stitch) | noisim "
    noisim_dir = workdir + '/noise/'
    os.chdir( noisim_dir )
    file = open( 'script_simulate_TDInoise_for_days_stitching.py' , 'r' )
    lines = file.readlines()
    file.close()
    if seed == 'random' :
        lines[ 11 ] = "seed = %s\n" % seed
    else :
       lines[ 11 ] = "seed = '%d'\n" % ( seed + 1000000 )
    file = open( 'script_simulate_TDInoise_for_days_stitching.py' , 'w' )
    file.writelines( lines )
    file.close()
    os.system( 'python script_simulate_TDInoise_for_days_stitching.py' )
    os.chdir( workdir )


    print " Merge signal and noise | merge "
    merge_dir = workdir + '/output'
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


    print " Analysis | main "
    main_dir = 'analysis_seed_%d_stime_0.5_hanning_ORFsim' % seed
    os.chdir( workdir )
    os.system( 'python script_setup_analysis.py' )
    file = open( 'setup_analysis.pkl' , 'rb' ) ; setup_main = cpkl.load( file ) ; file.close()
    setup_main['tsdir'] = merge_dir + '/data/'
    setup_main['workdir'] = main_dir
    file = open( 'setup_analysis.pkl' , 'wb' ) ; cpkl.dump( setup_main , file , -1 ) ; file.close()
    os.system( './x_analyse.py --do_psd --do_csd --do_avgpsd --do_X --do_G --do_S setup_analysis.pkl' )
    os.chdir( workdir ) 


    print " Post-analysis | post "
    post_dir = main_dir #'analysis_seed_%d_stime_0.5_hanning' % seed
    os.chdir( workdir )
    os.system( 'python script_setup_post_analysis.py' )
    file = open( 'setup_post_analysis.pkl' , 'rb' ) ; setup_post = cpkl.load( file ) ; file.close()
    setup_post['workdir'] = post_dir
    file = open( 'setup_post_analysis.pkl' , 'wb' ) ; cpkl.dump( setup_post , file , -1 ) ; file.close()
    os.system( './x_postanalyse.py --do_PP --do_P --do_stdPP --do_stdP --do_temp_optimals setup_post_analysis.pkl' )
    os.chdir( workdir )

    print 'Waiting for post-analysis jobs to finish...'
    jobidpaths = glob.glob( post_dir+'/do_PP_*jobid' ) + glob.glob( post_dir+'/do_P_*_jobid' ) + glob.glob(
        post_dir+'/do_stdPP_*jobid' ) + glob.glob( post_dir+'/do_stdP_*jobid' ) + glob.glob(
        post_dir+'/do_optimals_*jobid' )
    jobids = []
    for jobidpath in jobidpaths :
        file = open( jobidpath , 'r' )
        jobids += [ file.readlines()[0].rstrip() ]
        file.close()
    while True :
        jobxits = []
        for jobid in jobids :
            p = subprocess.Popen( 'qstat -f %s' % jobid ,
                                  shell=True , stdout=subprocess.PIPE , stderr=subprocess.PIPE )
            stdout , stderr = p.communicate()
            if len(stdout)==0 and stderr==('qstat: Unknown Job Id %s\n' % jobid) :
                jobxits += [1]
            elif 'exit_status = 0' in stdout :
                jobxits += [1]
            else :
                jobxits += [0]
        print 'jobxits' , jobxits
        if jobxits.count(1) == len(jobids) :
            print 'Post-analysis has finished.'
            break
        else :
            time.sleep( 20 )
            continue

    
#    print "Clean up data-simulation directories from this realisation "
#    os.chdir( sigsim_dir )
#    os.system( 'rm -r data/ *.sub *.out *jobids' )
#    os.chdir( workdir )
#    os.chdir( noisim_dir )
#    os.system( 'rm -r data/ *.sub *.out *jobids' )
#    os.chdir( workdir )
#    os.chdir( merge_dir )
#    os.system( 'rm -r data/ *.sub *.out *jobids' )
#    os.chdir( workdir )
