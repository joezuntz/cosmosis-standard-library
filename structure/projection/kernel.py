import numpy as np
import scipy.interpolate as interp

class KernelSpline(object):
    def __init__(self, x, y, clip=1.e-6, min_xmin=0.1, ymin=1.e-12, norm=True):
        self.x = x
        assert np.all(y>=0.)
        self.y = y

        #Find out if we've got zeros or very small values padding
        #the kernel - remove those that are unnecessary.
        quick_norm = np.sum(np.diff(x) * y[:-1])
        too_small = self.y/quick_norm < ymin
        good_inds, = np.where(~too_small)
        start, end = good_inds[0], good_inds[-1]
        if start>0:
            start-=1
            self.y[start]=0.
        if end<len(x)-1:
            end+=1
            self.y[end]=0.

        self.x, self.y = self.x[start:end+1], self.y[start:end+1]
        self.xmin, self.xmax = self.x[0], self.x[-1]

        #Find range of indices at start and end of kernel that are too small
        self.spline = interp.InterpolatedUnivariateSpline(self.x, self.y)
        self.norm = self.spline.integral(self.x[0], self.x[-1])
        if norm==True:
            self.y /= self.norm
            self.norm = 1.
            self.spline = interp.InterpolatedUnivariateSpline(self.x, self.y)

        #Compute a mimumim and maximum chi below and above which
        #a fraction clip of the integral is contained
        cumsum_y = (np.cumsum(self.y)-self.y[0])/(np.sum(self.y)-self.y[0]) #cumsum_y goes from 0 to 1
        pos_diff = np.zeros(len(cumsum_y), dtype=bool)
        pos_diff[0] = True
        pos_diff[1:] = np.diff(cumsum_y)>0.
        cumsum_of_x_spline = interp.InterpolatedUnivariateSpline(cumsum_y[pos_diff], self.x[pos_diff])
        self.xmin_clipped = cumsum_of_x_spline(clip)
        self.xmax_clipped = cumsum_of_x_spline(1-clip)

        #Also compute a mean and width
        dx = np.diff(self.x)
        x_mid = 0.5*(self.x[:-1]+self.x[1:])
        y_mid = 0.5*(self.y[:-1]+self.y[1:])
        self.mean = np.sum(dx * x_mid * y_mid) / np.sum(dx * y_mid)
        var = np.sum(dx * (x_mid - self.mean)**2 * y_mid) / np.sum(dx*y_mid)
        self.sigma = np.sqrt(var)
        #print("Setup kernel with mean(chi), sigma(chi) = %.3f, %.3f"%(self.mean, self.sigma))

    def __call__(self, x, fill_value=0.):
        if fill_value is not None:
            in_range = (x>self.xmin)*(x<self.xmax)
            return np.where(in_range, self.spline(x), fill_value)
        else:
            return self.spline(x)

class TomoNzKernel(object):
    def __init__(self, z, nzs, norm=True):
        self.z = z
        self.nzs = {}
        self.nbin = 0
        for i,nz in enumerate(nzs):
            self.nbin += 1
            if norm:
                nz_spline = interp.InterpolatedUnivariateSpline(self.z, nz)
                norm = nz_spline.integral(self.z[0], self.z[-1])
                nz = nz/norm
            self.nzs[i+1] = nz

        #These get set later
        self.nchi_splines = {} #n(chi) splines
        self.wchi_splines = {} #W(chi) splines

    @classmethod
    def from_block(cls, block, section_name, norm=True):
        nbin = block[section_name, "nbin"]
        z = block[section_name, "z"]
        nzs = []
        for i in range(1, nbin+1):
            nz = block[section_name, "bin_%d"%i]
            nzs.append(nz)
        return cls(z, nzs, norm=norm)

    def set_nofchi_splines(self, chi_of_z, dchidz, clip=1.e-6):
        chi = chi_of_z(self.z)
        dchidz = dchidz(self.z)
        self.nchi_splines = {}
        for i in range(1, self.nbin+1):
            self.nchi_splines[i] = KernelSpline(chi, self.nzs[i]/dchidz, clip=clip) 

    def set_wofchi_splines(self, chi_of_z, dchidz, a_of_chi, clip=1.e-6, nchi=5000):
        if len(self.nchi_splines) == 0:
            self.set_nofchi_splines(chi_of_z, dchidz, clip=clip)
        for i in range(1, self.nbin+1):
            nchi_spline = self.nchi_splines[i]
            chi_vals = np.linspace(0., nchi_spline.xmax, nchi)
            a_vals = a_of_chi(chi_vals)
            w_of_chi = np.zeros_like(chi_vals)
            for j,chi in enumerate(chi_vals):
                #integral \int_{chi}^{chi_max} dchi' n(chi') * (chi' - chi)/chi'
                chi_start = max(chi, nchi_spline.xmin)
                integrand = np.where(nchi_spline.x>chi, nchi_spline.y * (nchi_spline.x - chi)/nchi_spline.x, 0.)
                integrand_spline = interp.InterpolatedUnivariateSpline(nchi_spline.x, integrand)
                w_of_chi[j] = chi * integrand_spline.integral(chi_start, nchi_spline.xmax) / a_vals[j]
            self.wchi_splines[i] = KernelSpline(chi_vals, w_of_chi, norm=False)

    def get_kernel_spline(self, kernel_type, bin_number):
        if kernel_type == "N":
            return self.nchi_splines[bin_number]
        elif kernel_type == "W":
            return self.wchi_splines[bin_number]