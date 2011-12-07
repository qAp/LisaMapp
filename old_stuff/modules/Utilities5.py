
import numpy as np
import AnisotropySearch as AS
import scipy.integrate as Integrate
import scipy.special as scisp
import scipy.misc.common as scimcm




def GpretdiMichelson_to_G1Michelson( FT ) :

    L = 16.6782

    f = FT.s1.Offset1 + FT.s1.Cadence1 * np.arange( FT.s1.data.shape[0] )
    fdict = { 'Offset1' : FT.s1.Offset1 , 'Cadence1' : FT.s1.Cadence1 }

    Xdata , Ydata , Zdata = FT.s1.data , FT.s2.data , FT.s3.data

    e2 = np.exp( - 1j * 4 * np.pi * L * f )

    G1Xdata = ( e2 - 1 ) * Xdata
    G1Ydata = ( e2 - 1 ) * Ydata
    G1Zdata = ( e2 - 1 ) * Zdata

    G1X = AS.Coarsable( G1Xdata , **fdict )
    G1Y = AS.Coarsable( G1Ydata , **fdict )
    G1Z = AS.Coarsable( G1Zdata , **fdict )

    G1Michelson = AS.ShortTermFT( G1X , G1Y , G1Z )

    return G1Michelson




def G1Michelson_to_G1Sagnac( FT ) :

    L = 16.6782

    f = FT.s1.Offset1 + FT.s1.Cadence1 * np.arange( FT.s1.data.shape[0] )
    fdict = { 'Offset1' : FT.s1.Offset1 , 'Cadence1' : FT.s1.Cadence1 }

    Xdata , Ydata , Zdata = FT.s1.data , FT.s2.data , FT.s3.data
    e1 = np.exp( - 1j * 2 * np.pi * L * f )
    e2 = np.exp( - 1j * 4 * np.pi * L * f )
    e3 = np.exp( - 1j * 6 * np.pi * L * f )

    A = 1 - e1 + e2
    B = 1 - e1 - e2 + e3

    alphadata = ( A*Xdata + e1 * ( Ydata + Zdata ) ) / B
    betadata  = ( A*Ydata + e1 * ( Zdata + Xdata ) ) / B
    gammadata = ( A*Zdata + e1 * ( Xdata + Ydata ) ) / B


    alpha = AS.Coarsable( alphadata , **fdict )
    beta  = AS.Coarsable( betadata  , **fdict )
    gamma = AS.Coarsable( gammadata , **fdict )

    SagnacFT = AS.ShortTermFT( alpha , beta , gamma )
    
    return SagnacFT


def G1Sagnac_to_GmSagnac( FT ) :

    L = 16.6782

    f = FT.s1.Offset1 + FT.s1.Cadence1 * np.arange( FT.s1.data.shape[0] )
    fdict = { 'Offset1' : FT.s1.Offset1 , 'Cadence1' : FT.s1.Cadence1 }

    G1adata , G1bdata , G1cdata = FT.s1.data , FT.s2.data , FT.s3.data

    e3 = np.exp( - 1j * 6 * np.pi * L * f )

    Gmadata = ( e3 - 1 ) * G1adata
    Gmbdata = ( e3 - 1 ) * G1bdata
    Gmcdata = ( e3 - 1 ) * G1cdata
    
    Gma = AS.Coarsable( Gmadata , **fdict )
    Gmb = AS.Coarsable( Gmbdata , **fdict )
    Gmc = AS.Coarsable( Gmcdata , **fdict )
    
    GmSagnac = AS.ShortTermFT( Gma , Gmb , Gmc )
    
    return GmSagnac


def FracPhase_to_FracFreq( FT ) :
    
    f = FT.s1.Offset1 + FT.s1.Cadence1 * np.arange( FT.s1.data.shape[0] )
    fdict = { 'Offset1' : FT.s1.Offset1 , 'Cadence1' : FT.s1.Cadence1 }

    phase1data , phase2data , phase3data = FT.s1.data , FT.s2.data , FT.s3.data

    a = 1 / ( 1j * 2 * np.pi * f )

    freq1data = phase1data * a
    freq2data = phase2data * a
    freq3data = phase3data * a

    freq1 = AS.Coarsable( freq1data , **fdict )
    freq2 = AS.Coarsable( freq2data , **fdict )
    freq3 = AS.Coarsable( freq3data , **fdict )

    freqFT = AS.ShortTermFT( freq1 , freq2 , freq3 )

    return freqFT


#Some rotation transformation functions here


def SphericalToCartesian( sphericaltuple ) :
    r , phi , the = sphericaltuple
    x = r * np.sin( the ) * np.cos( phi )
    y = r * np.sin( the ) * np.sin( phi )
    z = r * np.cos( the )
    return x , y , z

def CartesianToSpherical( cartesiantuple ) :
    x , y , z = cartesiantuple
    r   = np.sqrt( x**2 + y**2 + z**2 )

    if r == 0 :
        ( r , phi , the ) = ( 0 , 0 , 0 )
        return r , phi , the

    the = np.arccos( z / r )
    arctan = np.arctan2( y , x )
    if arctan < 0 :
        phi = 2*np.pi + arctan
    else:
        phi = arctan
#    if x == 0 and y == 0 :
#        phi = 0
#    if y >= 0 and x > 0 :
#        phi = np.arctan( abs( y / x ) )
#    elif y >= 0 and x < 0 :
#        phi = np.pi - np.arctan( abs( y / x ) )
#    elif y <= 0 and x < 0 :
#        phi = np.pi + np.arctan( abs( y / x ) )
#    elif y <= 0 and x > 0 :
#        phi = 2*np.pi - np.arctan( abs( y / x ) )
#    elif y > 0 and x == 0 :
#        phi = np.pi / 2
#    elif y < 0 and x == 0 :
#        phi = 3 * np.pi / 2
    return r , phi , the




def RotationMatrix( angle , axis=1 ) :

    s , c = np.sin( angle ) , np.cos( angle )

    if axis == 1 :
        rmatrix = np.array( [ [ 1 , 0 , 0 ] ,
                                 [ 0 , c , -s ] ,
                                 [ 0 , s , c ] ] )
    elif axis == 2 :
        rmatrix = np.array( [ [ c , 0 , s ] ,
                                 [ 0 , 1 , 0 ] ,
                                 [ -s , 0 , c ] ] )
    elif axis == 3 :
        rmatrix = np.array( [ [ c , -s , 0 ] ,
                                 [ s , c , 0 ] ,
                                 [ 0 , 0 , 1 ] ] )
    else :
        raise ValueError , "Axis has to be 1 , 2 , or 3"
    return rmatrix
        


def RotationMatrix_323( phi , the , psi ) :
    rphi = RotationMatrix( -phi , axis=3 )
    rthe = RotationMatrix( -the , axis=2 )
    rpsi = RotationMatrix( -psi , axis=3 )
    rmatrix = np.dot( rpsi , np.dot( rthe , rphi ) )
    return rmatrix
    


def TransformPixels_Rot323( pixels , rphi , rthe , rpsi ) :
    """
    pixels -- numpy array. rows evenly span 180 degrees
                           columns even span 360 degrees
    rphi -- angle, about z-axis, of 1st rotation
    rthe -- angle, about y-axis, of 2nd rotation
    rpsi -- angle, about z-axis, of 3rd rotation
    """
    p = np.copy( pixels )
    rmatrix = RotationMatrix_323( rphi , rthe , rpsi )

    nlat , nlon = p.shape
    dlat , dlon = np.pi / ( nlat - 1 ) , 2 * np.pi / nlon

    colats = dlat * np.arange( nlat )
    lons   = dlon * np.arange( nlon )

    newp = np.zeros( ( nlat , nlon ) )
    for i,colat in enumerate( colats ) :
        for j,lon in enumerate( lons ) :
            kp = SphericalToCartesian( ( 1 , lon , colat ) ) 
            k = tuple( np.dot( np.array( kp ) , rmatrix ) )
            rnew , lonnew , colatnew = CartesianToSpherical( k )
            inew = int( np.around( colatnew / dlat ) )
            jnew = int( np.around( lonnew / dlon ) )
            if jnew == nlon :
                jnew = 0
            newp[ inew , jnew ] += np.copy( p[ i , j ] )
    return newp



def TransformPixels_Rot323_b( pixels , rphi , rthe , rpsi ) :
    """
    pixels -- numpy array. rows evenly span 180 degrees
                           columns even span 360 degrees
    rphi -- angle, about z-axis, of 1st rotation
    rthe -- angle, about y-axis, of 2nd rotation
    rpsi -- angle, about z-axis, of 3rd rotation
    """

    rotphi = RotationMatrix( rphi , axis=3 )
    rotthe = RotationMatrix( rthe , axis=2 )
    rotpsi = RotationMatrix( rpsi , axis=3 )
    rmatrix = np.dot( rotphi , np.dot( rotthe , rotpsi ) )

    oldp = np.copy( pixels )

    nlat , nlon = oldp.shape
    dlat , dlon = np.pi / ( nlat - 1 ) , 2 * np.pi / nlon
    colats = dlat * np.arange( nlat )
    lons   = dlon * np.arange( nlon )

    newp = np.zeros( ( nlat , nlon ) )
    for i , the in enumerate( colats ) :
        for j , phi in enumerate( lons ) :
            newk = SphericalToCartesian( ( 1 , phi , the ) )
            oldk = tuple( np.dot( np.array( newk ) , rmatrix ) )
            oldr , oldphi , oldthe = CartesianToSpherical( oldk )
            oldi = int( np.around( oldthe / dlat ) )
            oldj = int( np.around( oldphi / dlon ) )
            if oldj == nlon :
                oldj = 0
            newp[ i , j ] = oldp[ oldi , oldj ]

    return newp


def WignerD( the , l , m , n ) :
    """
    Wigner D matrices:  d_{mn}^{l}(the)
    the -- angle in radians
    l -- degree (integer)
    m -- order  (integer)
    n -- order  (integer)
    """
    if m >= n :
        factor = (-1)**(l-m) * np.sqrt( ( scimcm.factorial( l+m ) * scimcm.factorial( l-m ) ) /
                                           ( scimcm.factorial( l+n ) * scimcm.factorial( l-n ) ) )
        a = m + n
        b = m - n
        print a
        print b
        return factor * np.cos( the/2. )**a * ( -np.sin( the/2. ) )**b * scisp.jacobi( l-m , a , b )( -np.cos( the ) )

    else :
        return (-1)**(m-n) * WignerD( the , l , n , m )
    


def TransformMultipoles_Rot323( xlm , rphi , rthe , rpsi ) :

    xlm = np.copy( xlm )

    ntrunc = np.sqrt( xlm.shape[0] ) - 1
    if ( ntrunc % 1 ) != 0 :
        print "The length of the entered numpy array does not correspond to an integer ntrunc."
        raise ValueError

    indxpn = AS.getMLvec( ntrunc )

    ylm_rphi = np.array( [ np.exp( - 1j * ml[0] * rphi ) for ml in indxpn ] ) * xlm

    ylm_rthe = np.zeros( ylm_rphi.shape , dtype=ylm_rphi.dtype )
    for i , ml in enumerate( indxpn ) :
        m , l = ml
        nls   = [ ( n , l ) for n in range( -l , l + 1 , 1 ) ]
        ylns  = [ ylm_rphi[ indxpn.index( nl ) ] for nl in nls ]
        dlnms = [ WignerD( rthe , l , nl[0] , m ) for nl in nls ]
        ylm_rthe[ i ] = np.dot( np.array( dlnms ) , np.array( ylns ) )

    ylm_rpsi = np.array( [ np.exp( - 1j * ml[0] * rpsi ) for ml in indxpn ] ) * ylm_rthe
    return ylm_rpsi



#help from Teviet

def Equatorial_to_Galactic( delta , alpha ) :
    """
    Transforms from equatorial to galactic coordinates
    INPUT:
    delta --- declination [radians]
    alpha --- right ascension (RA) [radians]
    OUTPUT:
    b --- galactic latitude [radians]
    l --- galactic longitude [radians]
    """
    aNGP = np.radians( 192.8594813 )
    dNGP = np.radians( 27.1282511 )
    lascend = np.radians( 33 )

    arg1 = np.cos(delta)*np.cos(dNGP)*np.cos(alpha-aNGP) + np.sin(delta)*np.sin(dNGP)
    arg2 = np.sin(delta)*np.cos(dNGP) - np.cos(delta)*np.cos(alpha-aNGP)*np.sin(dNGP)
    arg3 = np.cos(delta)*np.sin(alpha-aNGP)

    b = np.arcsin( arg1 )
    arctan = np.arctan2( arg2 , arg3 )
    if arctan < 0 :
        O_O = 2*np.pi + arctan
    else:
        O_O = arctan
    l = np.mod( O_O + lascend , 2*np.pi )
    return b , l


def Galactic_to_Equatorial( b , l ) :
    """
    Transforms from galactic coordinates to equatorial coordinates.
    INPUT:
    b --- galactic latitude [radians]
    l --- galactic longitude [radians]
    OUTPUT:
    delta --- declination [radians]
    alpha --- right ascension (RA) [radians]
    """
    aNGP = np.radians( 192.8594813 )
    dNGP = np.radians( 27.1282511 )
    lascend = np.radians( 33 )

    arg1 = np.cos(b)*np.cos(dNGP)*np.sin(l-lascend) + np.sin(b)*np.sin(dNGP)
    arg2 = np.cos(b)*np.cos(l-lascend)
    arg3 = np.sin(b)*np.cos(dNGP) - np.cos(b)*np.sin(l-lascend)*np.sin(dNGP)

    delta = np.arcsin( arg1 )
    arctan = np.arctan2( arg2 , arg3 )
    if arctan < 0 :
        O_O = 2*np.pi + arctan
    else:
        O_O = arctan
    alpha = np.mod( O_O + aNGP , 2*np.pi )
    return delta , alpha


def Equatorial_to_Ecliptic( delta , alpha ) :
    """
    Transforms from equatorial coordinates to ecliptic coordinates
    INPUT:
    delta --- declination [radians]
    alpha --- right ascension (RA) [radians]
    OUTPUT:
    beta --- eclipitc latitude [radians]
    lamb --- ecliptic longitude [radians]
    """
    epsi = np.radians( 23.4392911 )

    arg1 = np.sin(delta)*np.cos(epsi) - np.cos(delta)*np.sin(alpha)*np.sin(epsi)
    arg2 = np.cos(delta)*np.sin(alpha)*np.cos(epsi) + np.sin(delta)*np.sin(epsi)
    arg3 = np.cos(delta)*np.cos(alpha)

    beta = np.arcsin( arg1 )
    arctan = np.arctan2( arg2 , arg3 )
    if arctan < 0 :
        lamb = 2*np.pi + arctan
    else:
        lamb = arctan
    return beta , lamb


def Ecliptic_to_Equatorial( beta , lamb ) :
    """
    Transforms from ecliptic to equatorial coordinates
    INPUT:
    beta --- ecliptic latitude [radians]
    lamb --- ecliptic longitude [radians]
    OUTPUT:
    delta --- declination [radians]
    alpha --- right ascension (RA) [radians]
    """
    epsi = np.radians( 23.4392911 )
    
    arg1 = np.cos(beta)*np.sin(lamb)*np.sin(epsi) + np.sin(beta)*np.cos(epsi)
    arg2 = np.cos(beta)*np.sin(lamb)*np.cos(epsi) - np.sin(beta)*np.sin(epsi)
    arg3 = np.cos(beta)*np.cos(lamb)

    delta = np.arcsin( arg1 )
    arctan = np.arctan2( arg2 , arg3 )
    if arctan < 0 :
        alpha = 2*np.pi + arctan
    else:
        alpha = arctan
    return delta , alpha


def getCoordinates_on_Grid( nlat=180 , nlon=360 , transform="" ) :
    """
    Transform coordinates of points on a grid to new coordinates.
    INPUT:
    nlat --- number of latitudes on the grid
    nlon --- number of longitudes on the grid
    transform --- [string] the coordinate transformation from the grid.
    OUTPUT:
    lats_on_grid --- new latitudes on the grid
    lons_on_grid --- new longitudes on the grid
    """
    dlat , dlon = np.pi / ( nlat - 1 ) , 2*np.pi / nlon
    lats = np.pi/2 - dlat*np.arange( nlat )
    lons = dlon*np.arange( nlon )

    lats_on_grid = np.empty( ( nlat , nlon ) )
    lons_on_grid = np.empty( ( nlat , nlon ) )
    for ii,lat in enumerate( lats ) :
        for jj,lon in enumerate( lons ) :

            if transform == 'Equatorial_to_Galactic' :
                coords = Equatorial_to_Galactic( lat , lon )
            elif transform == 'Galactic_to_Equatorial' :
                coords = Galactic_to_Equatorial( lat , lon )
            elif transform == 'Equatorial_to_Ecliptic' :
                coords = Equatorial_to_Ecliptic( lat , lon )
            elif transform == 'Ecliptic_to_Equatorial' :
                coords = Ecliptic_to_Equatorial( lat , lon )
            else:
                print "Coordintates transformation not recognised."
                raise ValueError

            lats_on_grid[ ii , jj ] = coords[0]
            lons_on_grid[ ii , jj ] = coords[1]
    return lats_on_grid , lons_on_grid



def bilinear_sum( F , Dx , Dy , dx , dy ) :
    """
    Returns f(x,y) from bilinear interpolation of the surrounding 4 grid points:
    { (xi,yi),(xi,yj),(xj,yi) and (xj,yj) }
    INPUT:
    F --- matrix np.array( [ [f(xi,yi),f(xi,yj)] , [f(xj,yi),f(xj,yj)] ] )
    Dx --- difference between two adjacent x grid points, xj-xi.
    Dy --- difference between two adjacent y grid points, yj-yi.
    dx --- x - xi
    dy --- y - yi
    OUTPUT:
    f(x,y) interpolated from f(xi,yi), f(xi,yj), f(xj,yi) and f(xj,yj)
    """
    xvec = np.array( [ (Dx-dx) , dx ] )
    yvec = np.array( [ (Dy-dy) , dy ] )
    return np.dot( xvec  , np.dot( F , yvec ) ) / ( Dx*Dy )



def getValues_on_Grid( vmatrix , latsgrid , lonsgrid ) :
    """
    Get value(from vmatrix) given the coordinates(given by latsgrid and lonsgrid) at each point in the grid
    INPUT:
    vmatrix  --- [array] values in a grid on the coordinates to transform from
    latsgrid --- [array] latitudes in the coordinates to transform from, in a grid on the coordinates to transform to
    lonsgrid --- [array] longitudes in the coordinates to transform from, in a grid on the coordinates to transform to
    OUTPUT:
    valuesgrid --- values in the coordinates to transform from, in a grid on the coodinates to transform to
    """
    if latsgrid.shape != lonsgrid.shape :
        print "latitude array and longitude array must have the same dimensions."
        raise Exception

    nlat , nlon = latsgrid.shape

    nlat0 , nlon0 = vmatrix.shape
    dlat0 , dlon0 = np.pi / (nlat0-1) , 2*np.pi / nlon0
    lat0s = np.pi/2 - dlat0 * np.arange( nlat0 )
    lon0s = dlon0 * np.arange( nlon0 )
    
    valuesgrid = np.zeros( ( nlat , nlon ) )
    for k in range( nlat ) :
        for l in range( nlon ) :

            lat , lon = latsgrid[ k , l ] , lonsgrid[ k , l ]
            j = int(  np.floor( ( lon - lon0s[0] ) / dlon0 )  )            
            i = int(  np.ceil( ( lat0s[0] - lat ) / dlat0 )  )
            
            dx , dy = ( lon-lon0s[j] ) , ( lat-lat0s[i] )
            Dx , Dy = dlon0 , dlat0
            
            if j+1 > nlon0 - 1 :
                jOO = 0
            else:
                jOO = j

            F = np.array( [ [ vmatrix[i,jOO] , vmatrix[i-1,jOO] ] ,
                               [ vmatrix[i,jOO+1] , vmatrix[i-1,jOO+1] ] ] )

            valuesgrid[ k , l ] = bilinear_sum( F , Dx , Dy , dx , dy )
    return valuesgrid


"""
GALACTIC MODELS
"""


def get_SunCentred_Number_Density( thing , r , t , p ) :
    """
    returns the number density of thing at (r,t,p) in sun-centred galatic spherical polar coordinates
    INPUT:
    thing --- string containing name of thing whose volume density is to be returned.
    r --- radius of spherical polar coordinates
    t --- polar angle of spherical polar coordinates [radians]
    p --- asimuthal angle of spherical polar coordinates [radians]
    OUTPUT :
    density --- number density [ L^{-3} ]
    """

    RG = 8 #distance between the Sun and the galactic centre [kpc]

    if thing == 'sech disc' :
        """
        Disc: sech ( as given in Matt's note )
        N --- total number of stars
        R0 --- scale radius of the galactic disc [kpc]
        z0 --- scale height of the galaxy [kpc]
        RG --- distance between the Sun and the galactic centre [kpc]
        """
        N , R0 , z0 = 10**7 , 2.5 , 0.5 
        sqrtarg = np.sqrt( r**2*np.sin(t)**2 + 2*RG*r*np.sin(t)*np.cos(p) + RG**2 )
        density = N * np.exp( - sqrtarg / R0 ) / ( 4*np.pi*R0*z0 * sqrtarg * np.cosh( r*np.cos(t)/z0 )**2 )
    elif thing == 'exponential disc' :
        """
        Disc: exponential
        rho0 --- normalisation factor
        R0 --- scale length
        z0 --- scale height
        """
        rho0 , R0 , z0 , RG = 10**7 , 2.5 , 0.5 , 8 #what normlatisation factor? 20110322
        sqrtarg = np.sqrt( r**2*np.sin(t)**2 + 2*RG*r*np.sin(t)*np.cos(p) + RG**2 )
        density = rho0 * np.exp( - np.abs( r*np.cos(t) ) / z0 ) * np.exp( - sqrtarg / R0 )
    elif thing == 'power-law spheroid' :
        """
        Spheroid: power-law formulation
        rho --- normalisation factor
        n --- power index
        a0 --- core radius
        """
        rho0 , a0 , n = 10**7 , 0.1 , 2
        sqrtarg = np.sqrt( r**2 + 2*RG*r*np.sin(t)*np.cos(p) + RG**2 )
        density = rho0 / ( a0**n + sqrtarg**n )
    elif thing == 'de Vaucouleurs spheroid' :
        """
        Spheroid: de Vaucouleurs
        rho --- normalisation factor
        Re --- effective/half-light radius
        """
        rho0 , Re = 10**7 , 0.1
        sqrtarg = np.sqrt( r**2 + 2*RG*r*np.sin(t)*np.cos(p) + RG**2 )
        density = rho0 * np.exp( -7.669 * ( sqrtarg / Re )**.25 ) * Re**.875 / sqrtarg**.875
    else :
        raise Exception , 'You have inputed an unidentified thing.'
    return density

def get_SunCentred_Number_Per_Omega( thing , t , p ) :
    """
    returns the number per unit solid angle of thing at (t,p) in sun-centred galactic spherical polar coordinates
    INPUT :
    thing --- string containing name of thing whose number per unit solid angle is to be returned.
    t --- polar angle of spherical polar coordinates [radians]
    p --- asimuthal angle of spherical polar coordinates [radians]
    OUTPUT :
    integral --- number per unit solid angle
    error --- error associated with the integration 
    """
    integral , error = Integrate.quad(
        lambda r : r**2 * get_SunCentred_Number_Density( thing , r , t , p ) , 0 , 375 )
    return integral , error

    
