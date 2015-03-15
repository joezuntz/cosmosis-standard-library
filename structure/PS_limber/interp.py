import numpy as np

def lagrange(x, xv, yv, n=3, check_bounds=True):
    """ n-point lagrange interpolation of the curve (xv,yv) at point x. """
    assert( n > 1 ) #behavior not yet defined for n <= 1.

    if check_bounds and ((x < xv[0]) or (x > xv[-1])):
        raise ValueError("x out of bounds. xlo, x, xhi = (%2.2e, %2.2e, %2.2e)" % (xv[0], x, xv[-1]))

    dx    = int( np.floor(0.5*n) )
    ixmin = min( max( np.searchsorted( xv, x, side='left' ) - dx, 0 ), len(xv) - n )
    idxs  = np.arange(0,n) + ixmin

    iv    = np.arange(0, n)
    xs    = xv[idxs]
    ys    = yv[idxs]

    lx = 0.0
    for i, tx, ty in zip(iv, xs, ys):
        lx += 1. * ty * np.prod( (x - xs)[iv != i] ) / np.prod( (tx - xs)[iv != i] )

    return lx

def linterp2d(x, y, xv, yv, f, check_bounds=True):
    """ bilinear 2d interpolation of the curve f(xv,yv) at point (x,y) """
    
    if check_bounds and ((x < xv[0]) or (x > xv[-1])):
        raise ValueError("x out of bounds. xlo, x, xhi = (%2.2e, %2.2e, %2.2e)" % (xv[0], x, xv[-1]))

    if check_bounds and ((y < yv[0]) or (y > yv[-1])):
        raise ValueError("y out of bounds. ylo, y, yhi = (%2.2e, %2.2e, %2.2e)" % (yv[0], y, yv[-1]))

    ix = min( max( np.searchsorted( xv, x, side='left' ) - 1, 0 ), len(xv) - 2 )
    iy = min( max( np.searchsorted( yv, y, side='left' ) - 1, 0 ), len(yv) - 2 )

    # bilinear interpolation of mat
    den = (xv[ix+1] - xv[ix]) * (yv[iy+1] - yv[iy])
    ret = ( f[  ix,  iy] * (xv[ix+1] - x) * (yv[iy+1] - y) +
            f[ix+1,  iy] * (x - xv[ix  ]) * (yv[iy+1] - y) +
            f[  ix,iy+1] * (xv[ix+1] - x) * (y - yv[iy  ]) +
            f[ix+1,iy+1] * (x - xv[ix  ]) * (y - yv[iy  ]) ) / den 

    return ret


