#!/usr/bin/env python


import sys
import numpy
import cPickle
import LISAresponse
import Utilities4 as U4


sourcename = ''  #Name your custom function something.

mapdir = ''


nlon =
nlat = 

ntrunc = 20

sky = LISAresponse.mySpharmt( nlon , nlat )
X = numpy.zeros( ( nlat , nlon ) )
for i,lat in enumerate( sky.lats ):
    theta = ( 90 - lat ) / 180 * numpy.pi
    for k,lon in enumerate( sky.lons ):
        lon = lon / 180 * numpy.pi
        X[ i,k ] = #Enter your function here in terms of theta and lon

skymap = U4.XSkyMap( X=X )
skymap.X_to_P()
skymap.X_to_Q()
skymap.ntrunc_equals( ntrunc )
skymap.P_to_plm()
skymap.Q_to_qlm()
skymap.plmqlm_to_xlm()

mapname = sourcename + '.pkl'
mappath = mapdir + mapname
file = open( mappath , 'wb' )
cPickle.dump( skymap , file , -1 )
file.close()


U4.project_SkyMap_Mollweide( skymap )

