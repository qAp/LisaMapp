#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import AnisotropySearch as AS
import PostProcess as PP
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser


parser = OptionParser( "usage: ./x_postprocess.py SETUP.pkl" )
( options , args ) = parser.parse_args()
if len( args ) < 1 :
    parser.error( 'You must specify a SETUP file of parameters!' )
setupname = args[0]
file = open( setupname , 'rb' ) ; setup = cpkl.load( file ) ; file.close()


# what to do ?
proj_X = setup['postproc']['proj_X']
proj_P , save_P =  setup['postproc']['proj_P'] , setup['postproc']['save_P']
proj_sigmamap , save_sigmamap = setup['postproc']['proj_Sg'] , setup['postproc']['save_Sg']
proj_SNRmap , save_SNRmap = setup['postproc']['proj_SNR'] , setup['postproc']['save_SNR']
proj_stdP = setup['postproc']['proj_stdP']

#sourcename = setup['inputdata']['sourcename']
Xdir = setup['maxlike']['Xdir']
Gdir = setup['maxlike']['Gdir']
#regularisation 
lmax , regMethod , regCutoff = setup['postproc']['lmax'] , setup['postproc']['regMethod'] , setup['postproc']['regCutoff']
#projection 
nlat , nlon  = setup['postproc']['nlat'] , setup['postproc']['nlon']





if True not in [ proj_X ,
                 proj_P , save_P , 
                 proj_sigmamap , save_sigmamap , 
                 proj_SNRmap , save_SNRmap ,
                 proj_stdP ] :
    print "Nothing to do."

if proj_X :
    print "Plotting X"
    mapnorm = 1.
    Xpath = Xdir + 'X.pkl'
    file = open( Xpath , 'rb' ) ; map_X = cpkl.load( file ) ; file.close()
    Xlm = AS.get_lmax_subset_from( map_X.xlm , lmax )
    Xmap = AS.xlmSkyMap( xlm = mapnorm * Xlm )
    Xmap.xlm_to_plm() ; Xmap.xlm_to_qlm()
    Xmap.create_sky( nlon=nlon , nlat=nlat )
    Xmap.plm_to_P() ; Xmap.qlm_to_Q()
    projdir = 'post_process/figures/X/'
    if projdir not in glob.glob( projdir ) :
        os.system( 'mkdir -p %s' % projdir )
    PP.project_SkyMap( Xmap , Ppath = projdir + 'Xreal.png' )
    


if True in [ proj_P , save_P ,
             proj_sigmamap , save_sigmamap ,
             proj_SNRmap , save_SNRmap ,
             proj_stdP ] :
    print "Calculating regularised inverse of Fisher matrix"
    Gpath = Gdir + 'G.pkl'
    fish = AS.FisherMatrix( Gpath , lmax=lmax )
    fish.regularise( regMethod=regMethod , regCutoff=regCutoff )
    regfishinv = fish.reginvert()

    if regMethod == None :
        regdir = 'post_process/nonreg/'
    else :
        regdir = 'post_process/reg%d_%.5f/' % ( regMethod , regCutoff )
    


if True in [ proj_P , save_P ,
             proj_SNRmap , save_SNRmap ] :
    print "Calculating clean map"
    mapnorm = 1.
    Xpath = Xdir + 'X.pkl'
    file = open( Xpath , 'rb' ) ; map_X = cpkl.load( file ) ; file.close()
    Xlm = AS.get_lmax_subset_from( map_X.xlm , lmax )
    Pdata = np.dot( regfishinv , Xlm )
    Pmap = AS.xlmSkyMap( xlm = mapnorm * Pdata )
    Pmap.xlm_to_plm() ; Pmap.xlm_to_qlm()
    Pmap.create_sky( nlon=nlon , nlat=nlat )
    Pmap.plm_to_P() ; Pmap.qlm_to_Q() ; Pmap.PQ_to_X()

    if proj_P :
        print "Plotting P"
        projdir = regdir + 'figures/P/'
        if projdir not in glob.glob( projdir ) :
            os.system( 'mkdir -p %s' % projdir )
        PP.project_SkyMap( Pmap , Ppath = projdir + 'Preal.png' )

    if save_P :
        print "Saving P to disk"
        datadir = regdir + 'data/P/' 
        if datadir not in glob.glob( datadir ) :
            os.system( 'mkdir -p %s' % datadir )
        datapath = datadir + 'P.pkl'
        file = open( datapath , 'wb' )
        cpkl.dump( Pmap , file , -1 )
        file.close()



if True in [ proj_sigmamap , save_sigmamap ,
             proj_SNRmap , save_SNRmap ,
             proj_stdP ] :
    print "Calculating the covariance of P_lm"
    covarm = np.dot( regfishinv , np.dot( fish.fish , regfishinv ) )


if True in [ proj_sigmamap , save_sigmamap ,
             proj_SNRmap , save_SNRmap ] :
    print "Calculating sigma map"
    sigmaP , lats , lons = AS.getSigmaMap( covarm , nlat , nlon )
    mapdata = np.reshape( sigmaP , ( nlat , nlon ) )
    sigmamap = AS.XSkyMap( X = mapdata )
    sigmamap.X_to_P() ; sigmamap.X_to_Q()

    if proj_sigmamap :
        print "Plotting sigma map"
        projdir = regdir + 'figures/sigmamap/'
        if projdir not in glob.glob( projdir ) :
            os.system( 'mkdir -p %s' % projdir )
        PP.project_SkyMap( sigmamap , Ppath = projdir + 'sigmaPreal.png' )

        rel_var = ( np.max(sigmamap.P) - np.min(sigmamap.P) ) / np.mean( sigmamap.P )
        N_keptSV = fish.N_keptSV
        infopath = projdir + 'sigmaPinfo.txt'
        file = open( infopath , 'w' )
        file.write( 'lmax = %d\n' % lmax )
        file.write( 'Relative variation = %f\n' % rel_var )
        file.write( 'regCutoff = %f\n' % regCutoff )
        file.write( 'Number of singular values retained = %d\n' % N_keptSV )
        file.close()

    if save_sigmamap :
        print "Saving sigma map"
        datadir = regdir + 'data/sigmamap/'
        if datadir not in glob.glob( datadir ) :
            os.system( 'mkdir -p %s' % datadir )
        datapath = datadir + 'sigmamap.pkl'
        file = open( datapath , 'wb' )
        cpkl.dump( sigmamap , file , -1 )
        file.close()

    
if proj_stdP :
    print "Calculating standard deviations of Plm"
    variance = np.diag( covarm )
    indxpn = AS.getMLvec( lmax )
    for i , ml in enumerate( indxpn ) :
        if variance[i] < 0 :
            print "Warning: negative variance for (m,l)=(%d,%d)" % ml
            print "Taking the absolute value"
            variance[i] = np.abs( variance[i] )
    std = np.sqrt( variance )
    std_repakd_list = []
    for l in range( lmax+1 ) :
        for m in range( -l , l+1 ) :
            std_repakd_list += [ std[ indxpn.index( (m,l) ) ] ]
    std_repakd = np.array( std_repakd_list )

    print "Plotting standard deviations of Plm"
    plt.figure()
    ax = plt.subplot(111)
    ax.semilogy( std_repakd )
    dmls = [ 2*l+1 for l in range( lmax+1 ) ]
    xticksat = np.array( [ sum(dmls[:ii]) for ii in range(len(dmls)) ] )
    xticksare = np.array( range( lmax+1 ) )
    plt.xticks( xticksat , xticksare , rotation=45 )
    ax.grid( b='on' )
    plt.xlim( ( 0 , std_repakd.shape[0] ) )
    plt.ylabel( 'standard deviation of $P_{lm}$' )
    plt.xlabel( 'l ((l,m) packed as $ l^{2} + l + m + 1 $)' )
    projdir = regdir + 'figures/stdSpH/'
    if projdir not in glob.glob( projdir ) :
        os.system( 'mkdir -p %s' % projdir )
    plt.savefig( projdir + 'stdSpH.png' )
    


if proj_SNRmap or save_SNRmap :
    print "Calculating SNR map"
    XSNR = Pmap.P / sigmamap.P 
    SNRmap = AS.XSkyMap( X=XSNR )
    SNRmap.X_to_P() ; SNRmap.X_to_Q()
    SNRmap.ntrunc_equals( lmax )
    SNRmap.P_to_plm() ; SNRmap.Q_to_qlm() , SNRmap.plmqlm_to_xlm()

    if proj_SNRmap :
        print "Plotting SNR map"
        projdir = regdir + 'figures/SNRmap/'
        if projdir not in glob.glob( projdir ) :
            os.system( 'mkdir %s' % projdir )
        PP.project_SkyMap( SNRmap , Ppath = projdir + 'SNRreal.png' )

    if save_SNRmap :
        print 'Saving SNR map'
        datadir = regdir + 'data/SNRmap/'
        if datadir not in glob.glob( datadir ) :
            os.system( 'mkdir -p %s' % datadir )
        file = open( datadir + 'SNRmap.pkl' , 'wb' )
        cpkl.dump( SNRmap , file , -1 )
        file.close()

        
        

