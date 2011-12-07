#!/usr/bin/env python


import sys
import numpy
import cPickle
import LISAresponse 
import Utilities4 as U4



mapdir = ''


l =                   #of Ylm
m =                   #of Ylm
Amplitude =           #Real number only please!

nlon = 360
nlat = 181

ntrunc = 20
skymap = U4.xlmSkyMap( ntrunc=ntrunc )
if m % 2 == 0 :
    if m == 0:
        skymap.alter_ml( False , (m,l,Amplitude) )
    else:
        skymap.alter_ml( False , (m,l,Amplitude/2.) , (-m,l,Amplitude/2.) )
else:
    skymap.alter_ml( False , (m,l,Amplitude/2.) , (-m,l,-Amplitude/2.) )
skymap.xlm_to_plm()
skymap.xlm_to_qlm()
skymap.create_sky( nlon=nlon , nlat=nlat )
skymap.plm_to_P()
skymap.qlm_to_Q()
skymap.PQ_to_X()


mapname = 'Y_l%d-m%d.pkl' % ( l , m )
mappath = mapdir + mapname
file = open( mappath , 'wb' )
cPickle.dump( skymap , file , -1 )
file.close()


U4.project_SkyMap_Mollweide( skymap )
