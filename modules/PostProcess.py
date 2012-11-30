import os
import glob
import cmath
import numpy as np
import scipy as sp
#import matplotlib
#import pylab as pl
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap , addcyclic




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


                                        

def SkyMap_to_datfile(skymap, datpath):
    """For each pixel in SKYMAP, its ecliptic longitude, latitude, and real and imaginary parts of its value
    are saved in a .dat file at DATPATH
    INPUT:
    skymap --- SkyMap object (from AnsiotropySearch.py)
    datpath --- path to save .dat file in
    """
    longrid , latgrid = sp.meshgrid( skymap.sky.lons , skymap.sky.lats )
    file = open(datpath, 'w')
    for k in range(skymap.sky.lons.shape[0]) :
        datatable = np.transpose( np.array([ longrid[:,k] , latgrid[:,k] ,
                                             skymap.P[:,k] , skymap.Q[:,k] ]) )
        np.savetxt( file , datatable , fmt='%20.10e' )
        file.write('\n')
    file.close()


def project_SkyMap(datpath,  Ppath = '', Qpath = '',
                   zlabel = 'amplitude'):
    if not (Ppath and Qpath): return
    lines = ['#!/usr/local/gnuplot -persist\n',
             'set terminal postscript enhanced color\n',
             'set pm3d\n',
             'set contour base\n', 
             'set view map\n', #45, 60
             'unset surface\n',
             'set grid\n',
             '\n',
             'set xtics 0, 30, 360\n', 'set ytics -90, 30, 90\n',
             '\n',
             'set xrange [0:360]\n', 'set yrange [-90:90]\n',
             '\n',
             'set xlabel "elon [degree]"\n', 'set ylabel "elat [degree]"\n', 'set zlabel "%s"\n' % zlabel,
             '\n']
    if Ppath:
        lines += ['set output "%s"\n' % Ppath,
                  'set title "real part"\n',
                  'splot "%s" using 1:2:3 with lines notitle\n' % datpath]
        Pdir = os.path.dirname(Ppath)
        if Pdir not in glob.glob(Pdir): os.system('mkdir -p %s' % Pdir)

    if Qpath:
        lines += ['set output "%s"\n' % Qpath,
                  'set title "imaginary part"\n',
                  'splot "%s" using 1:2:4 with lines notitle\n' % datpath]
        Qdir = os.path.dirname(Qpath)
        if Qdir not in glob.glob(Qdir): os.system('mkdir -p %s' % Qdir)
    
    file = open('project_SkyMap.plt', 'w'); file.writelines(lines); file.close()
    os.system('gnuplot project_SkyMap.plt')
    os.system('rm project_SkyMap.plt')
    return




