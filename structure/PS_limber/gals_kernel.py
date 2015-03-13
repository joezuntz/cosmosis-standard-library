'''

Galaxies PS

Use it as a template for specific galaxies

Part of a series of external utils that creates kernels for Limber integrals. This one is for Galaxies in general.\

You want to return a spline function W(l,chi,z) with l multipole chi comiving distance z redsfhit which is what is needed for limber.

EVERYTHING IS IN h UNITS

'''





class kern():

    def __init__(self, zdist, dndzfun, omm, h0, b=1.):
        '''
        Galaxies KERNEL (h units):

        Args:

            zdist: redshift distribution of the spline
            dndzfun: galaxies redshift distribution
            omm: Omega matter
            h0:hubble constant
            b: Galaxies bias


        Return:

            kern().w_lxz: kernel for limber integral


       '''

        wb = np.zeros(np.size(zdist))
        # use the z's from the P(k,z) array
        zmax = zdist[np.size(zdist) - 1]
        zmin = zdist[0]
        zmax = zdist[np.size(zdist) - 1]
        self.h0 = h0
        self.b = b
        self.omm = omm
        self.zmin = zmin
        self.zmax = zmax
        self.dndzfun = dndzfun
        self.norm = scipy.integrate.quad(dndzfun, self.zmin, self.zmax)[0]

    def w_lxz(self, l, x, z):
        '''
        Galaxies KERNEL (h units):

        w = dN/dz * b(z)

        with
       '''

        return dndzfun(z) / self.norm * b
