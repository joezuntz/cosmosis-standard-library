#python code for P(k) to C(l) calculations including
#Limber and exact ("non-Limber"). If this is too slow
#should be straightforward to convert to c.
from __future__ import print_function
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
from scipy.interpolate import interp1d, RectBivariateSpline
from scipy.integrate import quad
import sys
from LOG_HT import fft_log
import time

inv_sqrt2pi = 1./np.sqrt(2*np.pi)

def gaussian(x, mu, sigma):
    d = x-mu
    inv_sigma = 1./sigma
    return np.exp(-0.5*d*d*inv_sigma*inv_sigma) * inv_sqrt2pi * inv_sigma

def resample_power(k_in, k_out, P_in):
    #If k grid not the same, interpolate nl onto linear grid
    on_output_grid=True
    if len(k_in)!=len(k_out):
        on_output_grid = False
    elif not np.allclose(k_in, k_out, rtol=1.e-8):
        on_output_grid = False
    if not on_output_grid:
        P_out = np.zeros((P_in.shape[0],len(k_out)))
        for i in range(P_in.shape[0]):
            P_out[i] = interp1d(np.log(k_in), P_in[i], bounds_error=False, kind='cubic',
                                      fill_value=(P_in[i,0], P_in[i,-1]))(np.log(k_out))
    return P_out

class Kernel(object):
    """Class for storing n(z) kernel"""
    def __init__(self, x, n, clip=1.e-5):
        xmin, xmax = x.min(), x.max()
        self.spline = IUS(x, n)
        norm = self.spline.integral(x.min(), x.max())
        n /= norm
        if clip > 0.:
            clipped_xmin = xmin
            clipped_xmax = xmax
            dx = (x[1]-x[0])/10.
            area_min, area_max=0., 0.
            while area_min<clip:
                area_min = self.spline.integral(xmin, clipped_xmin)
                clipped_xmin += dx
            while area_max<clip:
                area_max = self.spline.integral(clipped_xmax, xmax)
                clipped_xmax -= dx
            xmin, xmax = clipped_xmin, clipped_xmax
            num_x = np.ceil((xmax-xmin)/dx).astype(int)
            x = np.linspace(xmin, xmax, num_x)
            n = self.spline(x)
            clipped_norm = self.spline.integral(xmin, xmax)
            n /= clipped_norm
        self.xmin, self.xmax = xmin, xmax
        self.x = x
        self.n = n
        self.spline = interp1d(x, n, bounds_error=False, 
            fill_value=0., kind='cubic')

    def get_chi_kernel(self, chi_of_z):
        chi_vals = chi_of_z(self.x)
        spl = IUS(chi_vals, self.n)
        norm = spl.integral(chi_vals[0],chi_vals[-1])

        kernel_interp = interp1d(chi_vals, self.n/norm, bounds_error=False, 
            kind='cubic', fill_value=0.)
        kernel_interp.chi_vals = chi_vals
        kernel_interp.mean = np.sum(chi_vals * self.n) / np.sum(chi_vals)
        kernel_interp.var = np.sum((chi_vals - kernel_interp.mean)**2 * self.n) / np.sum(chi_vals)
        kernel_interp.norm = norm
        return kernel_interp

def limber_integrand(chi, ell, kernel1_interp, kernel2_interp, pk_interp_logk):
    if chi<1.e-9:
        return 0.
    kernel1 = kernel1_interp(chi)
    if kernel2_interp is kernel1_interp:
        kernel2 = kernel1
    else:
        kernel2 = kernel2_interp(chi)
    nu = ell+0.5
    k = nu/chi
    pk = pk_interp_logk(chi, np.log(k))
    return kernel1 * kernel2 * pk / chi / chi

def get_cl_limber(ells, chimin, chimax, kernel1_interp, kernel2_interp, pk_interp_logk, rel_tol=1.e-4, abs_tol=0.):
    c_ells, c_ell_errs = np.zeros_like(ells), np.zeros_like(ells)
    for i,ell in enumerate(ells):
        c_ell, c_ell_err = quad(limber_integrand, chimin, chimax, 
                             args = (ell, kernel1_interp, kernel2_interp, pk_interp_logk),
                             epsrel=rel_tol, epsabs=abs_tol, limit=500)
        c_ells[i], c_ell_errs[i] = c_ell, c_ell_err
    return c_ells, c_ell_errs

def nearest_power_of_2(x):
    #Find nearest, greater power of 2 to x. 
    return 1<<(x-1).bit_length()

def get_dlogchi(dchi, chimax):
    #Since our chi values need to be log-spaced (for the fftlog), this dchi
    #will correspond to the last chi increment. Hence our dlogchi is given by
    #dlogchi = log(chi_N) - log(chi_{N-1}) = log(chi_N / chi_{N-1})
    #We can substitute chi_N = chimax and chi_{N-1} = chi_N - dchi to get
    #dlogchi = log(chimax / (chimax - dchi))
    dlogchi = np.log(chimax / (chimax-dchi))
    return dlogchi

def get_cl_exact(ells, chimin, chimax, dlogchi, kernel1_interp, kernel2_interp, 
    pk0_interp_logk, growth_interp, chi_pad_upper=2., chi_pad_lower=2., 
    verbose=True):
    """full integral is \int_0^\inf k^{-1} dk P(k,0) I_1(k) I_2(k)
    where I_1(k) = \int_0^{\inf} k dr_1 F_1(r_1) r^{-0.5} D(r_1) J_{l+0.5}(kr_1),
    and F_1(r_1) is the radial kernel for tracer 1.
    We want to use a fftlog for the I_1(k) calculation, so write it in the form 
    I_1(k) = \int_0^{\inf} k dr_1 f_1(r_1) J_{mu}(kr_1).
    So f_1(r_1_) = F_1(r_1) D(r_1) r^{-0.5}
    We actually do the integral in log(k), so calculate \int_0^\inf dlogk P(k,0) I_1(k) I_2(k).
    """
    q=0
    assert chimin>0.
    log_chimin, log_chimax = np.log(chimin), np.log(chimax)
    if verbose:
        print("padding chi values by e^%.2f/%.2f at lower/upper ends"%(chi_pad_lower,chi_pad_upper))
    log_chimin_padded, log_chimax_padded = log_chimin-chi_pad_lower, log_chimax+chi_pad_upper
    nchi_orig = np.ceil((log_chimax-log_chimin)/dlogchi).astype(int)
    nchi = nearest_power_of_2(nchi_orig) #use nchi that is a power of 2 for fast fft.
    log_chi_vals = np.linspace(log_chimin_padded, log_chimax_padded, nchi)
    chi_vals = np.exp(log_chi_vals)
    if verbose:
        print("chimin padded, chimax padded, nchi padded:", chi_vals[0], chi_vals[-1], len(chi_vals))
    growth_vals = growth_interp(chi_vals)
    kernel1_vals = kernel1_interp(chi_vals)
    auto=False
    if kernel2_interp is kernel1_interp:
        kernel2_vals = kernel1_vals
        auto=True
    else:
        kernel2_vals = kernel2_interp(chi_vals)

    cell = np.zeros_like(ells)
    for i_ell, ell in enumerate(ells):

        f1_vals = kernel1_vals * growth_vals * np.power(chi_vals, -0.5)
        k_vals, I_1 = fft_log(chi_vals, f1_vals, q, ell+0.5)
        if auto:
            I_2 = I_1
        else:
            f2_vals = kernel2_vals * growth_vals * np.power(chi_vals, -0.5)
            _, I_2 = fft_log(chi_vals, f2_vals, q, ell+0.5)
        logk_vals = np.log(k_vals)
        pk_vals = pk0_interp_logk(logk_vals)
        #Now we can compute the full integral \int_0^{\inf} k dk P(k,0) I_1(k) I_2(k)
        #We are values logspaced in k, so calculate as \int_0^{inf} dlog(k) P(k,0) I_1(k) I_2(k)
        integrand_vals = pk_vals * I_1 * I_2
        #Spline and integrate the integrand.
        integrand_interp = IUS(logk_vals, integrand_vals)
        integral = integrand_interp.integral(logk_vals.min(), logk_vals.max())
        cell[i_ell] = integral
    return cell

def test():
    #Calculate C(l)s for a Gaussian redshift windows 
    #(of a couple of different widths) for ells:
    ells = np.linspace(0, 1000, 100)
    mean_z = 1.
    sigma_zs = [0.05, 0.1, 0.3]
    sig_over_dchi = 50
    chi_pad_lower, chi_pad_upper = 3., 3.

    from os.path import join as pjoin
    from numpy import loadtxt
    import pylab

    #For the Limber calculation we need
    #- A chi(z) spline
    #- A 2d non-linear matter power spectrum spline
    #For the exact calculation we need in addition
    #- A P_lin(z=0,k) spline
    #- D(chi) growth spline
    #- A 2d P_nl-P_lin spline for the nonlinear part (for which we'll use Limber)
    test_output_dir = "test/test_output_limber_test"
    dists = pjoin(test_output_dir,"distances")
    #Read in bg z and chi.
    z_bg, chi_bg = loadtxt(pjoin(dists, "z.txt")), loadtxt(pjoin(dists, "d_m.txt"))
    chi_of_z_spline = interp1d(z_bg, chi_bg, kind='cubic', bounds_error=False,
        fill_value=0.)

    #Read in linear and non-linear matter P(k)
    plin_dir = pjoin(test_output_dir, "matter_power_lin")
    z_lin, k_lin, Pk_lin = (loadtxt(pjoin(plin_dir, "z.txt")),
                            loadtxt(pjoin(plin_dir, "k_h.txt")),
                            loadtxt(pjoin(plin_dir, "p_k.txt")))
    #make linear P(k,z=0) and growth splines
    pk0_interp_logk = interp1d(np.log(k_lin), Pk_lin[0], bounds_error=False, fill_value=0.)
    k_growth=1.e-3
    growth_ind=np.argmin(np.abs(k_lin-k_growth))
    growth_vals=np.sqrt(np.divide(Pk_lin[:,growth_ind], Pk_lin[0,growth_ind], 
                    out=np.zeros_like(Pk_lin[:,growth_ind]), where=Pk_lin[:,growth_ind]!=0.))
    growth_interp=interp1d(chi_of_z_spline(z_lin), growth_vals, kind='cubic',
        bounds_error=False, fill_value=0.)

    pnl_dir = pjoin(test_output_dir, "matter_power_nl")
    z_nl, k_nl, Pk_nl = (loadtxt(pjoin(pnl_dir, "z.txt")),
                         loadtxt(pjoin(pnl_dir, "k_h.txt")),
                         loadtxt(pjoin(pnl_dir, "p_k.txt")))   
    #make sure linear and nonlinear P(k)s on same z grid
    assert np.allclose(z_nl, z_lin)
    Pk_nl = resample_power(k_nl, k_lin, Pk_nl)
    k_nl = k_lin

    #make P_nl and P_nl - P_lin 2d splines
    P_lin_from_growth = np.outer(growth_vals**2, Pk_lin[0])
    Pk_nl_sub = Pk_nl - P_lin_from_growth
    pk_nl_interp_logk = RectBivariateSpline(chi_of_z_spline(z_nl), np.log(k_nl), Pk_nl)
    pk_sublin_interp_logk = RectBivariateSpline(chi_of_z_spline(z_nl), np.log(k_nl), Pk_nl_sub)

    #Also make a P_lin 2d spline for testing
    pk_lin_interp_logk = RectBivariateSpline(chi_of_z_spline(z_lin), np.log(k_lin), 
                                              P_lin_from_growth)

    #Do calculation for Gaussian windows of a few different widths
    #make a plot showing fractional difference as a funcion of ell.
    fig=pylab.figure()
    ax = fig.add_subplot(111)
    colors=[u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']

    for i,sigma_z in enumerate(sigma_zs):
        dz, zmax = 0.001, 2.00
        z_vals = np.arange(0.,zmax+dz,dz)
        nz_vals = gaussian(z_vals, mean_z, sigma_z)
        nz_kernel = Kernel(z_vals, nz_vals)
        #Make chi kernel
        chi_kernel = nz_kernel.get_chi_kernel(chi_of_z_spline)

        #Now we can do the calculation
        chimin, chimax = chi_kernel.chi_vals.min(), chi_kernel.chi_vals.max()
        cl_limber,_ = get_cl_limber(ells, chimin, chimax, chi_kernel, chi_kernel, pk_nl_interp_logk)
         
        #Exact calculation. For this we need to choose a dchi (and dlogchi)
        dchi = np.sqrt(chi_kernel.var) / 50 #choose a small dchi w.r.t to the width of the kernel.
        dlogchi = get_dlogchi(dchi, chi_kernel.chi_vals.max())

        cl_sublin,_ = get_cl_limber(ells, chimin, chimax, chi_kernel, chi_kernel, pk_sublin_interp_logk)
        cl_lin_exact = get_cl_exact(ells, chimin, chimax, dlogchi, chi_kernel, chi_kernel,
                                    pk0_interp_logk, growth_interp, chi_pad_lower=chi_pad_lower,
                                    chi_pad_upper=chi_pad_upper)
        #cl_lin_limber,_ = get_cl_limber(ells, chimin, chimax, chi_kernel, chi_kernel, pk_lin_interp_logk)
        cl_exact = cl_sublin + cl_lin_exact


        frac_diff_vals = cl_exact/cl_limber - 1
        pos = frac_diff_vals>0.
        ax.plot(ells[pos], frac_diff_vals[pos], ".", label=r"$\sigma_z=%.2f$"%sigma_z,
            color = colors[i] )
        ax.plot(ells[~pos], -frac_diff_vals[~pos], ".",
            mec = colors[i], mfc='none' )

    ax.set_xlabel(r"$l$")
    ax.set_ylabel(r"$C_{\mathrm{exact}}(l)/C_{\mathrm{limber}}(l) - 1$")
    ax.set_yscale('log')
    ax.legend()
    pylab.show()
    fig.savefig("cl_ratio.png")
    return

if __name__=="__main__":
    test()




