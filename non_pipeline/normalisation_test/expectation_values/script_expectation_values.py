#!/usr/bin/env python
import os
import sys
import glob
import numpy as np
import cPickle as cpkl
import myLISAmodule as mlisar
import AnisotropySearch as AS
import scipy.integrate as integrate
import synthlisa

P00 = 1/90.
tditype , tdigen = 'Michelson' , 'G2'
stime = 8.0


fNyq = 1 / (2*stime)

"""
\gamma^{IJ}_{00}(f)s
"""

def get_gAA( f ) :
    try:
        g00 = 288. * np.ones( f.shape )
    except TypeError :
        g00 = 288.
    return g00

def get_gAE( f ) :
    try:
        g00 = 288. * np.ones( f.shape )
    except TypeError :
        g00 = 288.
    return g00

def get_gAT( f ) :
    try:
        g00 = 288. * np.ones( f.shape )
    except TypeError :
        g00 = 288.
    return g00

def get_gEA( f ) :
    try:
        g00 = 288. * np.ones( f.shape )
    except TypeError :
        g00 = 288.
    return g00

def get_gEE( f ) :
    try:
        g00 = 288. * np.ones( f.shape )
    except TypeError :
        g00 = 288.
    return g00

def get_gET( f ) :
    try:
        g00 = 288. * np.ones( f.shape )
    except TypeError :
        g00 = 288.
    return g00

def get_gTA( f ) :
    try:
        g00 = 288. * np.ones( f.shape )
    except TypeError :
        g00 = 288.
    return g00

def get_gTE( f ) :
    try:
        g00 = 288. * np.ones( f.shape )
    except TypeError :
        g00 = 288.
    return g00

def get_gTT( f ) :
    try:
        g00 = 288. * np.ones( f.shape )
    except TypeError :
        g00 = 288.
    return g00


"""
Power spectral densities
"""
def get_PhAhA( f ) :
    return P00 * get_gAA(f)
def get_PhAhE( f ) :
    return P00 * get_gAE(f)
def get_PhAhT( f ) :
    return P00 * get_gAT(f)
def get_PhEhA( f ) :
    return P00 * get_gEA(f)
def get_PhEhE( f ) :
    return P00 * get_gEE(f)
def get_PhEhT( f ) :
    return P00 * get_gET(f)
def get_PhThA( f ) :
    return P00 * get_gTA(f)
def get_PhThE( f ) :
    return P00 * get_gTE(f)
def get_PhThT( f ) :
    return P00 * get_gTT(f)

def get_PnAnA( f ) :
    return 16e39 * mlisar.get_tdiNSD( tditype , tdigen , 'A' , 'A' , f )
def get_PnAnE( f ) :
    return 16e39 * mlisar.get_tdiNSD( tditype , tdigen , 'A' , 'E' , f )
def get_PnAnT( f ) :
    return 16e39 * mlisar.get_tdiNSD( tditype , tdigen , 'A' , 'T' , f )
def get_PnEnA( f ) :
    return 16e39 * mlisar.get_tdiNSD( tditype , tdigen , 'E' , 'A' , f )
def get_PnEnE( f ) :
    return 16e39 * mlisar.get_tdiNSD( tditype , tdigen , 'E' , 'E' , f )
def get_PnEnT( f ) :
    return 16e39 * mlisar.get_tdiNSD( tditype , tdigen , 'E' , 'T' , f )
def get_PnTnA( f ) :
    return 16e39 * mlisar.get_tdiNSD( tditype , tdigen , 'T' , 'A' , f )
def get_PnTnE( f ) :
    return 16e39 * mlisar.get_tdiNSD( tditype , tdigen , 'T' , 'E' , f )
def get_PnTnT( f ) :
    return 16e39 * mlisar.get_tdiNSD( tditype , tdigen , 'T' , 'T' , f )

def get_PAA( f ) :
    return get_PhAhA(f) + get_PnAnA(f)
def get_PAE( f ) :
    return get_PhAhE(f) + get_PnAnE(f)
def get_PAT( f ) :
    return get_PhAhT(f) + get_PnAnT(f)
def get_PEA( f ) :
    return get_PhEhA(f) + get_PnEnA(f)
def get_PEE( f ) :
    return get_PhEhE(f) + get_PnEnE(f)
def get_PET( f ) :
    return get_PhEhT(f) + get_PnEnT(f)
def get_PTA( f ) :
    return get_PhThA(f) + get_PnTnA(f)
def get_PTE( f ) :
    return get_PhThE(f) + get_PnTnE(f)
def get_PTT( f ) :
    return get_PhThT(f) + get_PnTnT(f)

"""
Write some PSDs and CSDs to disk
"""
n = 2000
df = 1 / ( n*stime )
Nf = n / 2 - 1
f = df + df * np.arange( Nf )

PhAhA = get_PhAhA( f )
PhAhE = get_PhAhE( f )
PhAhT = get_PhAhT( f )
PhEhA = get_PhEhA( f )
PhEhE = get_PhEhE( f )
PhEhT = get_PhEhT( f )
PhThA = get_PhThA( f )
PhThE = get_PhThE( f )
PhThT = get_PhThT( f )

PnAnA = get_PnAnA( f )
PnAnE = get_PnAnE( f )
PnAnT = get_PnAnT( f )
PnEnA = get_PnEnA( f )
PnEnE = get_PnEnE( f )
PnEnT = get_PnEnT( f )
PnTnA = get_PnTnA( f )
PnTnE = get_PnTnE( f )
PnTnT = get_PnTnT( f )

PAA = PnAnA + PhAhA
PAE = PnAnE + PhAhE
PAT = PnAnT + PhAhT
PEA = PnEnA + PhEhA
PEE = PnEnE + PhEhE
PET = PnEnT + PhEhT
PTA = PnTnA + PhThA
PTE = PnTnE + PhThE
PTT = PnTnT + PhThT

tableAA = np.transpose( np.array( [ f , np.real(PhAhA),np.imag(PhAhA) ,
                                    np.real(PnAnA),np.imag(PnAnA) , np.real(PAA),np.imag(PAA) ] ) )
tableAE = np.transpose( np.array( [ f , np.real(PhAhE),np.imag(PhAhE) ,
                                    np.real(PnAnE),np.imag(PnAnE) , np.real(PAE),np.imag(PAE) ] ) )
tableAT = np.transpose( np.array( [ f , np.real(PhAhT),np.imag(PhAhT) ,
                                    np.real(PnAnT),np.imag(PnAnT) , np.real(PAT),np.imag(PAT) ] ) )
tableEA = np.transpose( np.array( [ f , np.real(PhEhA),np.imag(PhEhA) ,
                                    np.real(PnEnA),np.imag(PnEnA) , np.real(PEA),np.imag(PEA) ] ) )
tableEE = np.transpose( np.array( [ f , np.real(PhEhE),np.imag(PhEhE) ,
                                    np.real(PnEnE),np.imag(PnEnE) , np.real(PEE),np.imag(PEE) ] ) )
tableET = np.transpose( np.array( [ f , np.real(PhEhT),np.imag(PhEhT) ,
                                    np.real(PnEnT),np.imag(PnEnT) , np.real(PET),np.imag(PET) ] ) )
tableTA = np.transpose( np.array( [ f , np.real(PhThA),np.imag(PhThA) ,
                                    np.real(PnTnA),np.imag(PnTnA) , np.real(PTA),np.imag(PTA) ] ) )
tableTE = np.transpose( np.array( [ f , np.real(PhThE),np.imag(PhThE) ,
                                    np.real(PnTnE),np.imag(PnTnE) , np.real(PTE),np.imag(PTE) ] ) )
tableTT = np.transpose( np.array( [ f , np.real(PhThT),np.imag(PhThT) ,
                                    np.real(PnTnT),np.imag(PnTnT) , np.real(PTT),np.imag(PTT) ] ) )

file = open( 'PSD.dat' , 'w' )
file.write('#PAA\n')
file.write( ('#'+'%19s'+6*'%20s')
            % ('freq[Hz]', 'Re{PhIhJ}','Im{PhIhJ}','Re{PnInJ}','Im{PnInJ}','Re{PIJ}','Im{PIJ}\n') )
np.savetxt( file , tableAA , fmt='%20.10e' )
file.write('\n\n')
file.write('#PAE\n')
file.write( ('#'+'%19s'+6*'%20s')
            % ('freq[Hz]', 'Re{PhIhJ}','Im{PhIhJ}','Re{PnInJ}','Im{PnInJ}','Re{PIJ}','Im{PIJ}\n') )
np.savetxt( file , tableAE , fmt='%20.10e' )
file.write('\n\n')
file.write('#PAT\n')
file.write( ('#'+'%19s'+6*'%20s')
            % ('freq[Hz]', 'Re{PhIhJ}','Im{PhIhJ}','Re{PnInJ}','Im{PnInJ}','Re{PIJ}','Im{PIJ}\n') )
np.savetxt( file , tableAT , fmt='%20.10e' )
file.write('\n\n')
file.write('#PEA\n')
file.write( ('#'+'%19s'+6*'%20s')
            % ('freq[Hz]', 'Re{PhIhJ}','Im{PhIhJ}','Re{PnInJ}','Im{PnInJ}','Re{PIJ}','Im{PIJ}\n') )
np.savetxt( file , tableEA , fmt='%20.10e' )
file.write('\n\n')
file.write('#PEE\n')
file.write( ('#'+'%19s'+6*'%20s')
            % ('freq[Hz]', 'Re{PhIhJ}','Im{PhIhJ}','Re{PnInJ}','Im{PnInJ}','Re{PIJ}','Im{PIJ}\n') )
np.savetxt( file , tableEE , fmt='%20.10e' )
file.write('\n\n')
file.write('#PET\n')
file.write( ('#'+'%19s'+6*'%20s')
            % ('freq[Hz]', 'Re{PhIhJ}','Im{PhIhJ}','Re{PnInJ}','Im{PnInJ}','Re{PIJ}','Im{PIJ}\n') )
np.savetxt( file , tableET , fmt='%20.10e' )
file.write('\n\n')
file.write('#PTA\n')
file.write( ('#'+'%19s'+6*'%20s')
            % ('freq[Hz]', 'Re{PhIhJ}','Im{PhIhJ}','Re{PnInJ}','Im{PnInJ}','Re{PIJ}','Im{PIJ}\n') )
np.savetxt( file , tableTA , fmt='%20.10e' )
file.write('\n\n')
file.write('#PTE\n')
file.write( ('#'+'%19s'+6*'%20s')
            % ('freq[Hz]', 'Re{PhIhJ}','Im{PhIhJ}','Re{PnInJ}','Im{PnInJ}','Re{PIJ}','Im{PIJ}\n') )
np.savetxt( file , tableTE , fmt='%20.10e' )
file.write('\n\n')
file.write('#PTT\n')
file.write( ('#'+'%19s'+6*'%20s')
            % ('freq[Hz]', 'Re{PhIhJ}','Im{PhIhJ}','Re{PnInJ}','Im{PnInJ}','Re{PIJ}','Im{PIJ}\n') )
np.savetxt( file , tableTT , fmt='%20.10e' )
file.close()




"""
Write some g00(f)s to disk (what follows below will using the same set of g00(f) values)
"""
#n = 500
#df = 1 / ( n*stime )
#Nf = n / 2 - 1
#f = df + df * np.arange( Nf )

#Or get coarse frequencies from an analysis run
file = open( '../analysis_seed_1_stime_0.5_hanning_ORFsim/GW_slope_0/AE/cIJ/d002.pkl' , 'rb' )
cIJdict = cpkl.load( file )
file.close()
df = cIJdict['f'].Cadence1
Nf = cIJdict['f'].data.shape[0]
f = cIJdict['f'].data

gAApath = 'gAA_f0_%.6f_df_%.6f_Nf_%d.dat' % ( df , df , Nf ) 
if gAApath not in glob.glob( gAApath ) :
    gAA = get_gAA( f )
    file = open( gAApath , 'w' )
    print >> file , ( '#'+'%19s' + 2*'%20s' ) % ( 'freq[Hz]' , 'Re{gAA}' , 'Im{gAA}' )
    for k in range( f.shape[0] ) :
        print >> file , ( 3*'%20.10f' ) % ( f[k] , np.real( gAA[k] ) , np.imag( gAA[k] ) )
    file.close()
else :
    file = open( gAApath , 'r' )
    lines = file.readlines()
    file.close()
    gAAlist = []
    for k in range( 1 , Nf+1 ) :
        realpart = float( lines[k].strip().split()[-2] )
        imagpart = float( lines[k].strip().split()[-1] )
        gAAlist += [ realpart + 1j*imagpart ]
    gAA = np.array( gAAlist )

gAEpath = 'gAE_f0_%.6f_df_%.6f_Nf_%d.dat' % ( df , df , Nf ) 
if gAEpath not in glob.glob( gAEpath ) :
    gAE = get_gAE( f )
    file = open( gAEpath , 'w' )
    print >> file , ( '#'+'%19s' + 2*'%20s' ) % ( 'freq[Hz]' , 'Re{gAE}' , 'Im{gAE}' )
    for k in range( f.shape[0] ) :
        print >> file , ( 3*'%20.10f' ) % ( f[k] , np.real( gAE[k] ) , np.imag( gAE[k] ) )
    file.close()
else :
    file = open( gAEpath , 'r' )
    lines = file.readlines()
    file.close()
    gAElist = []
    for k in range( 1 , Nf+1 ) :
        realpart = float( lines[k].strip().split()[-2] )
        imagpart = float( lines[k].strip().split()[-1] )
        gAElist += [ realpart + 1j*imagpart ]
    gAE = np.array( gAElist )

gEApath = 'gEA_f0_%.6f_df_%.6f_Nf_%d.dat' % ( df , df , Nf ) 
if gEApath not in glob.glob( gEApath ) :
    gEA = get_gEA( f )
    file = open( gEApath , 'w' )
    print >> file , ( '#'+'%19s' + 2*'%20s' ) % ( 'freq[Hz]' , 'Re{gEA}' , 'Im{gEA}' )
    for k in range( f.shape[0] ) :
        print >> file , ( 3*'%20.10f' ) % ( f[k] , np.real( gEA[k] ) , np.imag( gEA[k] ) )
    file.close()
else :
    file = open( gEApath , 'r' )
    lines = file.readlines()
    file.close()
    gEAlist = []
    for k in range( 1 , Nf+1 ) :
        realpart = float( lines[k].strip().split()[-2] )
        imagpart = float( lines[k].strip().split()[-1] )
        gEAlist += [ realpart + 1j*imagpart ]
    gEA = np.array( gEAlist )

gEEpath = 'gEE_f0_%.6f_df_%.6f_Nf_%d.dat' % ( df , df , Nf ) 
if gEEpath not in glob.glob( gEEpath ) :
    gEE = get_gEE( f )
    file = open( gEEpath , 'w' )
    print >> file , ( '#'+'%19s' + 2*'%20s' ) % ( 'freq[Hz]' , 'Re{gEE}' , 'Im{gEE}' )
    for k in range( f.shape[0] ) :
        print >> file , ( 3*'%20.10f' ) % ( f[k] , np.real( gEE[k] ) , np.imag( gEE[k] ) )
    file.close()
else :
    file = open( gEEpath , 'r' )
    lines = file.readlines()
    file.close()
    gEElist = []
    for k in range( 1 , Nf+1 ) :
        realpart = float( lines[k].strip().split()[-2] )
        imagpart = float( lines[k].strip().split()[-1] )
        gEElist += [ realpart + 1j*imagpart ]
    gEE = np.array( gEElist )    



PhAhA = P00 * gAA
PhAhE = P00 * gAE
PhEhA = P00 * gEA
PhEhE = P00 * gEE

PnAnA = get_PnAnA( f )
PnAnE = get_PnAnE( f )
PnEnA = get_PnEnA( f )
PnEnE = get_PnEnE( f )

PAA = PnAnA + PhAhA
PAE = PnAnE + PhAhE
PEA = PnEnA + PhEhA
PEE = PnEnE + PhEhE
#write fcoarse,g00,PAE,PAA and PEE to disk so they can be compared with g00,CAE,PAA,PEE from analysis pipeline
datatosave = np.transpose( np.array( [ f , np.real(gAE),np.imag(gAE) ,
                                       np.real(PAE),np.imag(PAE) , np.real(PAA),np.imag(PAA) ,
                                       np.real(PEE),np.imag(PEE) ] ) )
cIJ_ana_dir = 'cIJ_ana'
os.system( 'mkdir %s' % cIJ_ana_dir )
file = open( cIJ_ana_dir+'/anasd.dat' , 'w' )
file.write( ('#'+'%19s'+8*'%20.10s' + '\n') % ( 'freq[Hz]', 'Re{gAE}','Im{gAE}',
                                                'Re{PAE}','Im{PAE}','Re{PAA}','Im{PAA}', 'Re{PEE}','Im{PEE}' ) )
np.savetxt( file , datatosave , fmt='%20.10e' )
file.close()



"""
Analytical variances 
"""
sAA = df * np.sum( PAA )
sAE = df * np.sum( np.real( PAE ) )
sEA = df * np.sum( np.real( PEA ) )
sEE = df * np.sum( PEE )
print 'sAA' , sAA
print 'sAE' , sAE
print 'sEA' , sEA
print 'sEE' , sEE



"""
Analytical Sigmas
"""
N = 10800
flow , fhigh = 0.0005 , 0.0620
Ndays = 363

window = np.hanning( N )
winfactor = ( np.sum( window**2 ) / N )**2 / ( np.sum( window**4 ) / N )


ilow = int( np.round( ( flow - f[0] ) / df ) )
ihigh = int( np.round( ( fhigh - f[0] ) / df ) )
if ilow != 0 :
    print f[ ilow-1 ] ,
else :
    print 'none'  ,
print f[ ilow ]  , f[ ilow+1 ]
print f[ ihigh-1 ] , f[ ihigh ] ,
if ihigh != f.shape[0]-1 :
    print f[ ihigh+1 ]
else :
    print 'none'
print 'ilow , ihigh' , ilow , ihigh
print 'norm' , winfactor*N*stime*df

#print 'frequency integral in Fisher matrix: value' , df * np.sum(
#    np.abs( gAE[ilow:ihigh+1] )**2 / ( PAA[ilow:ihigh+1] * PEE[ilow:ihigh+1] ) )
#print '2*winfactor*N*stime * df ' , 2*winfactor*N*stime * df
print 'np.sum( np.abs( gAE[ilow:ihigh+1] )**2 )' , np.sum( np.abs( gAE[ilow:ihigh+1] )**2 )
#print 'np.sum( 1. / ( PAA[ilow:ihigh+1] * PEE[ilow:ihigh+1] ) )' , np.sum( 1. / ( PAA[ilow:ihigh+1] * PEE[ilow:ihigh+1] ) )
#print 'frequency integral in strong signal bias matrix: value' , df * np.sum(
#    np.abs( gAE[ilow:ihigh+1] )**4 / ( PAA[ilow:ihigh+1] * PEE[ilow:ihigh+1] )**2 )
#print '2*winfactor*N*stime * df ' , 2*winfactor*N*stime * df
print 'np.sum( np.abs( gAE[ilow:ihigh+1] )**4 )' , np.sum( np.abs( gAE[ilow:ihigh+1] )**4 )
#print 'np.sum( 1 / ( PAA[ilow:ihigh+1] * PEE[ilow:ihigh+1] )**2 )' , np.sum( 1 / ( PAA[ilow:ihigh+1] * PEE[ilow:ihigh+1] )**2 )


Gana = 2*winfactor*N*stime * df * np.sum( np.abs( gAE[ilow:ihigh+1] )**2 / ( PAA[ilow:ihigh+1] * PEE[ilow:ihigh+1] ) )

Sana = 2*winfactor*N*stime*P00**2 * df * np.sum( np.abs( gAE[ilow:ihigh+1] )**4 / ( PAA[ilow:ihigh+1] * PEE[ilow:ihigh+1] )**2 )

var_P00_weak_ana = 1 / Gana
var_P00_strong_bias_ana = Sana / Gana**2
var_P00_ana = var_P00_weak_ana + var_P00_strong_bias_ana
var_P00_opt_ana = var_P00_ana / Ndays

print 'Analytical optimal variance = ' , var_P00_opt_ana
print 'Analytical sigma = ' , np.sqrt( var_P00_ana )
print 'Analytical optimal sigma = ' , np.sqrt( var_P00_opt_ana )

print 'Analytical sigma (weak-signal limit) = ' , np.sqrt( var_P00_weak_ana )
print 'Analytical optimal sigma (weak-signal limit) = ' , np.sqrt( var_P00_weak_ana / Ndays )






