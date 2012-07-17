#!/usr/bin/env python
import os
import sys
import glob
import numpy as np


sd_seedsavg_daysavg_cIJ_path = '..//sd_seedsavg_daysavg_cIJ/avgsd.dat'
sd_seedsavg_daysavg_cII_path = '..//sd_seedsavg_daysavg_cII/avgsd.dat'
cIJana_path = '..//expectation_values/cIJ_ana/anasd.dat'
P00 = 1 / 90.



file = open( sd_seedsavg_daysavg_cII_path , 'r' )
lines = file.readlines()
file.close()
cIItable = []
for line in lines[4:] :
    cIItable += [ [ float(entry) for entry in line.split() ] ]

cIItable = np.array( cIItable )
f_cII = cIItable[:,0]
g00_cII = cIItable[:,1] + 1j*cIItable[:,2]
P11_cII = cIItable[:,3] #+ 1j*cIItable[:,4]
P22_cII = cIItable[:,5] #+ 1j*cIItable[:,6]



file = open( sd_seedsavg_daysavg_cIJ_path , 'r' )
lines = file.readlines()
file.close()

cIJtable = []
for line in lines[4:] :
    cIJtable += [ [float(entry) for entry in line.split()] ]

cIJtable = np.array( cIJtable )
f_cIJ = cIJtable[:,0]
g00_cIJ = cIJtable[:,1] + 1j*cIJtable[:,2]
P12_cIJ = cIJtable[:,3] + 1j*cIJtable[:,4]
P11_cIJ = cIJtable[:,5] # + 1j*cIJtable[:,6]
P22_cIJ = cIJtable[:,7] # + 1j*cIJtable[:,8]



print 'f_cII same as f_cIJ:' , np.allclose( f_cII , f_cIJ )
print 'P11_cII same as P11_cIJ:' , np.allclose( P11_cII , P11_cIJ )
print 'P22_cII same as P22_cIJ:' , np.allclose( P22_cII , P22_cIJ )
print 'g00_cII same as g00_cIJ:' , np.allclose( g00_cII , g00_cIJ )


"""
Max. Estimators using above loaded spectral denstities and glms
"""
N = 10800
stime = 8.0
flow , fhigh = 0.0005 , 0.0620
Ndays = 363

window = np.hanning( N )
winfactor = ( np.sum( window**2 ) / N )**2 / ( np.sum( window**4 ) / N )

f = np.copy( f_cII )
df = f[1] - f[0]
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
    
#Estimate G00 using sd_seedsavg_daysavg_cII
#print 'frequency integral in Fisher matrix: value' , df * np.sum(
#    np.abs( g00_cII[ilow:ihigh+1] )**2 / ( P11_cII[ilow:ihigh+1] * P22_cII[ilow:ihigh+1] ) )
#print '2*winfactor*N*stime * df' , 2*winfactor*N*stime * df
print 'np.sum( np.abs( g00_cII[ilow:ihigh+1] )**2 )' , np.sum( np.abs( g00_cII[ilow:ihigh+1] )**2 )
#print 'np.sum( 1. / ( P11_cII[ilow:ihigh+1] * P22_cII[ilow:ihigh+1] ) )' , np.sum( 1. / ( P11_cII[ilow:ihigh+1] * P22_cII[ilow:ihigh+1] ) )
Gana = 2*winfactor*N*stime * df * np.sum( np.abs( g00_cII[ilow:ihigh+1] )**2 / ( P11_cII[ilow:ihigh+1] * P22_cII[ilow:ihigh+1] ) )

#Estimate S00 using sd_seedsavg_daysavg_cII 
#print 'frequency integral in strong signal bias matrix: value' , df * np.sum(
#    np.abs( g00_cII[ilow:ihigh+1] )**4 / ( P11_cII[ilow:ihigh+1] * P22_cII[ilow:ihigh+1] )**2 )
#print '2*winfactor*N*stime * df ' , 2*winfactor*N*stime * df
print 'np.sum( np.abs( g00_cIJ[ilow:ihigh+1] )**4 )' , np.sum( np.abs( g00_cII[ilow:ihigh+1] )**4 )
#print 'np.sum( 1 / ( P11_cII[ilow:ihigh+1] * P22_cII[ilow:ihigh+1] )**2 )' , np.sum( 1 / ( P11_cII[ilow:ihigh+1] * P22_cII[ilow:ihigh+1] )**2 )

Sana = 2*winfactor*N*stime*P00**2 * df * np.sum( np.abs( g00_cII[ilow:ihigh+1] )**4 /
                                                 ( P11_cII[ilow:ihigh+1] * P22_cII[ilow:ihigh+1] )**2 )

#Estimate X00 usindg sd_seedsavg_daysavg_cIJ
Xana = 2*winfactor*N*stime * df * ( np.sum( np.real( np.conj(g00_cIJ[ilow:ihigh+1])*P12_cIJ[ilow:ihigh+1] )
                                            / ( P11_cIJ[ilow:ihigh+1] * P22_cIJ[ilow:ihigh+1] ) ) )

#Estimate P00
P00est = Xana / Gana

var_P00_weak_ana = 1 / Gana
var_P00_strong_bias_ana = Sana / Gana**2
var_P00_ana = var_P00_weak_ana + var_P00_strong_bias_ana
var_P00_opt_ana = var_P00_ana / Ndays

print 'P00est' , P00est
print '(P00est - P00)' , (P00est - P00)
print 'Analytical optimal variance = ' , var_P00_opt_ana
print 'Analytical sigma = ' , np.sqrt( var_P00_ana )
print 'Analytical optimal sigma = ' , np.sqrt( var_P00_opt_ana )

print 'Analytical sigma (weak-signal limit) = ' , np.sqrt( var_P00_weak_ana )
print 'Analytical optimal sigma (weak-signal limit) = ' , np.sqrt( var_P00_weak_ana / Ndays )






#Spectral densities residuals
file = open( cIJana_path , 'r' )
lines = file.readlines()
file.close()
cIJanatable = []
for line in lines[1:] :
    cIJanatable += [ [ float(entry) for entry in line.split() ] ]

cIJanatable = np.array( cIJanatable )

f_ana = cIJanatable[:,0]
g00_ana = cIJanatable[:,1] + 1j * cIJanatable[:,2]
P12_ana = cIJanatable[:,3] + 1j * cIJanatable[:,4]
P11_ana = cIJanatable[:,5] + 1j * cIJanatable[:,6]
P22_ana = cIJanatable[:,7] + 1j * cIJanatable[:,8]

resi_f = f_cIJ - f_ana
resi_g00 = (g00_cIJ - g00_ana) / g00_ana
resi_P12 = (P12_cIJ - P12_ana) / P12_ana
resi_P11 = (P11_cIJ - P11_ana) / P11_ana
resi_P22 = (P22_cIJ - P22_ana) / P22_ana

resi_table = np.transpose( np.array( [ f_ana , np.real(resi_g00) , np.imag(resi_g00) ,
                                       np.real(resi_P12) , np.imag(resi_P12) ,
                                       np.real(resi_P11) , np.imag(resi_P11) , np.real(resi_P22) , np.imag(resi_P22) ] ) )
file = open( 'sd_residuals.dat' , 'w' )
print '#Spectral densities residuals'
print >> file ,  ( '#'+'%19s'+8*'%20s' ) % ('freq[Hz]','Re{g00}','Im{g00}','Re{P12}','Im{P12}','Re{P11}','Im{P11}','Re{P22}','Im{P22}')
np.savetxt( file , resi_table , fmt='%20.10e' )
file.close()


