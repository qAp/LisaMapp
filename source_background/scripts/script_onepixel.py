#!/usr/bin/env python


import sys
import numpy
import cPickle
import LISAresponse
import Utilities4 as U4


mapdir  = ''

lon0 =             #0 ~ 359 only
lat0 =             #90 ~ -90 only
Amplitude =        #Real number only please!

ntrunc = 20

nlon = 360
nlat = 181
sky = LISAresponse.mySpharmt( nlon , nlat )
X = numpy.zeros( ( nlat , nlon ) , dtype='complex' )
X[ list( sky.lats ).index( lat0 ) , list( sky.lons ).index( lon0 ) ] = Amplitude

skymap = U4.XSkyMap( X=X )
skymap.X_to_P()
skymap.X_to_Q()
skymap.ntrunc_equals( ntrunc )
skymap.P_to_plm()
skymap.Q_to_qlm()
skymap.plmqlm_to_xlm()

mapname = 'pixel_lon%d-lat%d.pkl' % ( lon0 , lat0 )
mappath = mapdir + mapname
file = open( mappath , 'wb' )
cPickle.dump( skymap , file , -1 )
file.close()

U4.project_SkyMap_Mollweide( skymap )


