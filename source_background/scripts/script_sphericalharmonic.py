#!/usr/bin/env python
import os
import sys
import cPickle as cpkl
import AnisotropySearch as AS
import PostProcess as PP



mapdir = './'


l = 0                  #of Ylm
m = 0                  #of Ylm
Amplitude = 1e-38          #Real number only please!

nlon = 360
nlat = 181

ntrunc = 20
skymap = AS.xlmSkyMap( ntrunc=ntrunc )
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


mapname = 'Y_l%d_m%d.pkl' % ( l , m )
mappath = mapdir + mapname
file = open( mappath , 'wb' )
cpkl.dump( skymap , file , -1 )
file.close()


PP.project_SkyMap( skymap )
