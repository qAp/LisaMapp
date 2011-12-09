#!/usr/bin/env python
import os
import sys
import glob
import cPickle as cpkl
import numpy as np
import optparse
import AnisotropySearch as AS

usage = """
%prog PPDIR STDPPDIR\n
PPDIR -- directory containing the PPs (daily clean maps)\n
STDPPDIR -- directory containing the standard deviations of the PPs\n
SSDIR -- directory containing the daily strong signal covariance bias matrices\n
SUMPATH -- file to save annual summaries in
"""
parser = optparse.OptionParser( usage = usage )

parser.add_option( '-d' , '--day' , action='append' , dest='days' , nargs=1 , type='int' , 
                   help='Days whose clean maps and standard deviations to be included' )

parser.add_option( '--lmax' , action='store' , dest='lmax' , nargs=1 , type='int' , default=0 ,
                   help='Maximum l for which the PPs and stdPPs were computed.' )

( options , args ) = parser.parse_args()

if len( args ) < 4 :
    parser.error( "You must specify PPDIR, STDPPDIR, SSDIR and SUMPATH! Type ./x_Popt.py -h for help." )

PPdir , stdPPdir , SSdir , sumpath = args[ :4 ]

days , PPs , stdPPs , stdPPws , SSs = [] , [] , [] , [] , []
for day in options.days :
    print 'Day %d' % day
    PPpath = PPdir + '/PP_d%03d_lmax_%d.pkl' % ( day , options.lmax )
    stdPPpath = stdPPdir + '/stdPP_d%03d_lmax_%d.pkl' % ( day , options.lmax )
    stdPPwpath = stdPPdir + '/stdPP_d%03d_lmax_%d_weak.pkl' % ( day , options.lmax )
    SSpath = SSdir + '/SS_d%03d.pkl' % day
    if PPpath not in glob.glob( PPpath ) :
        print 'PP not found at %s. Skip' % PPpath ; continue
    if stdPPpath not in glob.glob( stdPPpath ) :
        print 'stdPP not found at %s.  Skip' % stdPPpath ; continue
    if stdPPwpath not in glob.glob( stdPPwpath ) :
        print 'stdPP_weak not found at %s.  Skip' % stdPPwpath ; continue
    if SSpath not in glob.glob( SSpath ) :
        print 'stdPP not found at %s.  Skip' % SSpath ; continue
    file = open( PPpath , 'rb' ) ; PP = cpkl.load( file ) ; file.close()
    file = open( stdPPpath , 'rb' ) ; stdPP = cpkl.load( file ) ; file.close()
    file = open( stdPPwpath , 'rb' ) ; stdPPw = cpkl.load( file ) ; file.close()
    file = open( SSpath , 'rb' ) ; SS = cpkl.load( file ) ; file.close()
    days += [ day ] ; PPs += [ np.real( PP.plm[0] ) ] ; stdPPs += [ np.real( stdPP[0] ) ] ; stdPPws += [ np.real( stdPPw[0] ) ] ; SSs += [ SS['S'].data[0] ]

days , PPs , stdPPs , stdPPws , SSs = np.array( days ) , np.array( PPs ) , np.array( stdPPs ) , np.array( stdPPws ) , np.array( SSs )
P_opt = np.sum( PPs / stdPPs**2 ) / np.sum( 1. / stdPPs**2 )
stdP_opt = np.sqrt( 1 / np.sum( 1. / stdPPs**2 ) )
stdPw_opt = np.sqrt( 1 / np.sum( 1. / stdPPws**2 ) )

summary = { 'days':days , 'PPs':PPs , 'stdPPs':stdPPs , 'stdPPws':stdPPws , 
            'P_opt':P_opt , 'stdP_opt':stdP_opt , 'stdPw_opt':stdPw_opt }

sumdir = os.path.dirname( sumpath )
if sumdir not in glob.glob( sumdir ) :
    os.system( 'mkdir -p %s' % sumdir )

file = open( sumpath , 'wb' ) ; cpkl.dump( summary , file , -1 ) ; file.close()


        
