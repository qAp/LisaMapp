
import cmath
import numpy as np
import scipy as sp
import matplotlib
import pylab as pl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap , addcyclic




def project_SkyMap_Mollweide( skymap , Pnorm=1 , Qnorm=1
                              , Ppath=None , Qpath=None ):

    Pw , lonsw = addcyclic( skymap.P , skymap.sky.lons )
    Qw , lonsw = addcyclic( skymap.Q , skymap.sky.lons )
    meshlon , meshlat = sp.meshgrid( lonsw , skymap.sky.lats )

    projection = Basemap( projection='moll' , lon_0=180 , resolution='c' )
#    projection = Basemap( projection='ortho' , lat_0=60 , lon_0=180 , resolution='c' )
    projection.drawmapboundary()
    x , y = projection( meshlon , meshlat )

    if Ppath == None :
        pass
    else :
        fig = plt.figure( figsize = (8,8) )
        ax = fig.add_axes( [ 0.05 , 0.15 , 0.8 , 0.8 ] )
        projection.contourf( x , y , Pnorm*Pw , 30 )
        projection.drawmeridians( np.arange(0,360,30) )
        projection.drawparallels( np.arange(-90,90,30) , labels=[1,0,0,0] )
        pos = ax.get_position()
        l , b , w , h = pos.bounds
        cax = plt.axes( [ l+w+0.03 , b , 0.04 , h ] )
        plt.colorbar( cax=cax , orientation='vertical' , format='%.1e' )
        pl.savefig( Ppath )

    if Qpath == None :
        pass
    else :
        fig = plt.figure( figsize = (8,8) )
        ax = fig.add_axes( [ 0.05 , 0.15 , 0.8 , 0.8 ] )
        projection.contourf( x , y , Qnorm*Qw , 30 )
        projection.drawmeridians( np.arange(0,360,30) )
        projection.drawparallels( np.arange(-90,90,30) , labels=[1,0,0,0] )
        pos = ax.get_position()
        l , b , w , h = pos.bounds
        cax = plt.axes( [ l+w+0.03 , b , 0.04 , h ] )
        plt.colorbar( cax=cax , orientation='vertical' )
        pl.savefig( Qpath )
    return


                                        

def project_SkyMap( skymap , Pnorm=1 , Qnorm=1 ,
                    Ppath = None , Qpath = None ) :

    Lons , Lats = np.meshgrid( skymap.sky.lons , skymap.sky.lats )

    if Ppath == None :
        pass
    else :
        plt.figure()

        plt.contourf( Lons , Lats , Pnorm * skymap.P )
        plt.grid( b='on' ) ; plt.xlim( (0,360) ) ; plt.ylim( (-90,90) )
        plt.xticks( np.arange( 0 , 360+30 , 30 ) )
        plt.yticks( np.arange( -90 , 90+30 , 30 ) )
        plt.xlabel( 'longitude[degrees]' ) ; plt.ylabel( 'latitude[degrees]' )
        plt.title( 'real part' )
        plt.colorbar( format='%.1e' )
        plt.savefig( Ppath )

    if Qpath == None :
        pass
    else :
        plt.figure()
        plt.contourf( Lons , Lats , Qnorm * skymap.Q )
        plt.grid( b='on' ) ; plt.xlim( (0,360) ) ; plt.ylim( (-90,90) )
        plt.xticks( np.arange( 0 , 360+30 , 30 ) )
        plt.yticks( np.arange( -90 , 90+30 , 30 ) )
        plt.xlabel( 'longitude[degrees]' ) ; plt.ylabel( 'latitude[degrees]' )
        plt.title( 'imaginary part' )
        plt.colorbar( format='%.1e' )
        plt.savefig( Qpath )
    return


