#!/usr/bin/env python
import os
import sys
import glob
import numpy as np

Nb = 2

days = range( 1 , 10+1 )
stime = 0.5
seed = 100
GWSpectralSlope = 0
lmax = 0
compute_ORF_SpHs = True
Ppath = '/gpfs1/JC0311443/workhere/stochasGW/Mapp/skymaps/library/sphericalhar\
monics/Y_l0_m0_x1_over_90_lmax_0/Y_l0_m0.pkl'
orfdir = ''
tsdir = ''


#Divide the days into batches so they can run in parallel
if days == None :
    print 'No days are selected.  Nothing to do.' ; sys.exit()
else :
    Ndays = len( days )
    if Ndays % Nb == 0 :
        Ndb = Ndays / Nb
        days_batches = [ days[ b*Ndb : (b+1)*Ndb  ]  for b in range( Nb ) ]
    elif Ndays % Nb > 0 :
        Ndb = Ndays / ( Nb-1 )
        Ndbl = Ndays % ( Nb-1 )
        days_batches = [ days[ b*Ndb : (b+1)*Ndb ] for b in range( Nb-1 ) ] + [ days[ (Nb-1)*Ndb : ] ]
    else :
        raise Exception , "Both the number of days and number of batches have to be postivie integer!"

