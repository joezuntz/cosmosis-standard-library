#Library for calculating IA contributions up to delta**4
#Includes 'linear alignment' and 'tidal torque' contributions
#contact Niall MacCrann with any issues.
import scipy.interpolate as interp
import numpy as np
import scipy.integrate

pi=np.pi

ACC_LOW = {"log_alpha_coarse_min":-4., "log_alpha_coarse_max":4., 
            "n_alpha_coarse":50, "alpha_fine_min":0.98,
            "alpha_fine_max":1.02, "n_alpha_fine":50,
            "mu_fine_min":0.99, "n_mu_coarse":50,
            "n_mu_fine":50}

ACC_MED = {"log_alpha_coarse_min":-6., "log_alpha_coarse_max":6., 
            "n_alpha_coarse":500, "alpha_fine_min":0.98,
            "alpha_fine_max":1.02, "n_alpha_fine":200,
            "mu_fine_min":0.99, "n_mu_coarse":200,
            "n_mu_fine":200}

ACC_HI = {"log_alpha_coarse_min":-5., "log_alpha_coarse_max":5., 
            "n_alpha_coarse":1000, "alpha_fine_min":0.98,
            "alpha_fine_max":1.02, "n_alpha_fine":200,
            "mu_fine_min":0.99, "n_mu_coarse":200,
            "n_mu_fine":200}


def no_smooth(k):
    return np.ones_like(k)

class Del2kInterpFunc(object):
    """A class for interpolation of the dimensionless power spectrum,
    del^2(k) = k^3 P(k,z) / 2pi^2"""
    def __init__(self,ks,Pk):
        self.ks=ks
        self.Pk=Pk
        self.nk=len(ks)
        self.logks=np.log(ks)
        self.del2 = self.ks**3 * Pk / 2 / pi**2
        self.interp_func=interp.interp1d(self.logks, np.log(self.del2), bounds_error=False, 
                                        fill_value=np.nan, kind='cubic')  
        self.sigma2_S = self.calc_sigma2_S()

    def __call__(self, k):
        del2 = np.atleast_1d(np.exp(self.interp_func(np.log(k))))
        del2[np.isnan(del2)] = 0.
        if len(del2)==1:
            return del2[0]
        else:
            return del2

    def P(self, k):
        return 2*pi**2 * self(k) / k**3

    def sigma2_r(self,R,uk=5):
        return scipy.integrate.quad(self.sigint, self.logks[0], self.logks[-1], args=(R), epsrel=1e-7, limit=10000, epsabs=0.)[0]

    def calc_sigma2_S(self):
        return scipy.integrate.quad(self.sigma2_S_int, self.logks[0], self.logks[-1], 
            epsrel=1e-7, limit=10000, epsabs=0.)[0]

    def calc_sigma2n_S(self, n, S=no_smooth):
        return scipy.integrate.quad(self.sigma2n_S_integrand, self.logks[0], self.logks[-1], 
            args=(n,S), epsrel=1e-7, limit=10000, epsabs=0.)[0]

    def sigma2n_S_integrand(self, lnk, n, S_k):
        """n is the power to which P(k) is raised,
           S is a smoothing function"""    
        k = np.exp(lnk)    
        return 4 * pi * k**3 * (S_k(k) * self.P(k))**n / (2*pi)**3

    def sigma2_S_int(self, lnk):
        k=np.exp(lnk)
        return 2 * pi**2 * (self(k))**2 / np.exp(lnk)**3

    def sigint(self,lnk, R):
        k=np.exp(lnk)
        x=k*R
        w=3*(-x*np.cos(x)+np.sin(x))/x**3
        return w**2 * self(k)

def trapz(x,y):
    return 0.5*np.sum(np.diff(x)*(y[:-1]+y[1:]), axis=0)

def trapz_2d(x,y):
    x_diff=np.diff(x)
    _,diff_grid = np.meshgrid(np.empty(y.shape[1]), x_diff)
    #print 'diff_grid',diff_grid
    return 0.5 * np.sum(diff_grid * (y[:-1] + y[1:]), axis=0)

def grow_Pk(P0,zs,Dzs,n=2):
    """Get P(k,z) from P(k,0) by rescaling using the growth factor"""
    P_out=np.zeros((len(zs), len(P0)))
    for i,(z,Dz) in enumerate(zip(zs,Dzs)):
        P_out[i] = P0*(Dz**n)
    return P_out


def alpha_mu_integral(k_lins, del2interp, kernel_func, Pq=1, Pq2=0, 
            log_alpha_coarse_min=-5.,log_alpha_coarse_max=5.,n_alpha_coarse=200,
            alpha_fine_min=0.98,alpha_fine_max=1.02,n_alpha_fine=200,
            mu_fine_min=0.99,n_mu_coarse=200, n_mu_fine=200,
            smooth_q=no_smooth, smooth_q2=no_smooth, absP_kernel=False,
            smooth_q_power=1, smooth_q2_power=1):
    """This performs integrals of the type in eqn. 19 of 
    Mackey+ 01 (astro-ph/0106364). It doesn't include the prefactor.
    The 'kernel_func' argument expects a function of alpha and mu depending
    on the IA term being calculated. For example of the tidal torque II term
    use the hehe function.
    The integral is performed in two steps in both alpha and mu - with
    finer binning when either variable is close to 1 (since here the denominator
    goes to zero eqn. 19). 
    The function checks that the specified ranges are such that both 
    can't be exactly 1.

    @param k_lins           The values of k at which to evaluate the integral
    @param del2interp       A Del2kInterpFunc object
    @param kernel_func      An IA kernel function
    @param smooth           Optional function to smooth functions of k. [default: No smoothing]
    @param Pq2              integral includes del2(q2)
    @absP_kernel            Set to true if the kernel assumes the integral is over P(k) 
                            rather than del2
    @returns the integral obvs.
    """
    mu_min,mu_max=-1.,1.
    alphas=np.sort(np.concatenate((np.logspace(log_alpha_coarse_min,log_alpha_coarse_max,n_alpha_coarse),
                                   np.linspace(alpha_fine_min,alpha_fine_max,n_alpha_fine))))
    mus=np.concatenate((np.linspace(mu_min,mu_fine_min,n_mu_coarse),
                        np.linspace(mu_fine_min,mu_max,n_mu_fine)))
    #Do some checks on integral ranges. 
    #Check if alpha and mu both get very close to 1.
    if np.any(np.isclose(alphas,1,atol=1e-9,rtol=0.)):
        raise ValueError("""alpha gets rather close to 1....
                        choose different alpha values...""")
    #Check mus are sensible
    assert (mu_fine_min>=-1.) and (mu_fine_min<=1.)
    assert np.abs(alpha_fine_min)>np.abs(10**log_alpha_coarse_min)
    assert np.abs(alpha_fine_max)<np.abs(10**log_alpha_coarse_max)
    alpha_integrand=np.zeros((len(alphas),len(k_lins)))
    for ialpha,alpha in enumerate(alphas):
        q = alpha*k_lins
        mu_integrand=np.zeros((len(mus),len(k_lins)))
        for imu,mu in enumerate(mus):
            kernel = kernel_func(alpha,mu)
            q2 = np.sqrt(1 + alpha**2 - 2*alpha*mu)*k_lins
            if Pq2>0: 
                mu_integrand[imu] = kernel * smooth_q2(q2)**smooth_q2_power
            else:
                del2q2 = del2interp(q2)
                mu_integrand[imu] = del2q2**Pq2 * kernel * smooth_q2(q2)**smooth_q2_power
                if absP_kernel:
                    mu_integrand[imu] *= (2*pi**2 / q2**3)**Pq2

        mu_integral = trapz_2d(mus, mu_integrand)
        if Pq>0:
            alpha_integrand[ialpha] = (smooth_q(q)**smooth_q_power * del2interp(q)**Pq
                                         * mu_integral)
            if absP_kernel:            
                alpha_integrand[ialpha] *= (2*pi**2 / q**3)**Pq
        else:
            alpha_integrand[ialpha] = smooth_q(q)**smooth_q_power * mu_integral

    alpha_integral = trapz_2d(alphas, alpha_integrand)      
    return alpha_integral / (2*pi)**2

def alpha_integral(k_lins, del2interp, kernel_func, 
            log_alpha_coarse_min=-5.,log_alpha_coarse_max=5.,n_alpha_coarse=300,
            alpha_fine_min=0.98,alpha_fine_max=1.02,n_alpha_fine=300,
            smooth=no_smooth, absP_kernel=False):
    """
    @param k_lins           The values of k at which to evaluate the integral
    @param del2interp       A Del2kInterpFunc object
    @param kernel_func      An IA kernel function
    @param smooth           Optional function to smooth functions of k. [default: No smoothing]
    @absP_kernel            Set to true if the kernel assumes the integral is over P(k) 
                            rather than del2
    @returns the integral obvs.
    """
    alphas=np.sort(np.concatenate((np.logspace(log_alpha_coarse_min,log_alpha_coarse_max,n_alpha_coarse),
                                   np.linspace(alpha_fine_min,alpha_fine_max,n_alpha_fine))))

    #Do some checks on integral ranges. 
    #Check if alpha and mu both get very close to 1.
    if np.any(np.isclose(alphas,1,atol=1e-9,rtol=0.)):
        raise ValueError("""alpha gets rather close to 1....
                        choose different alpha values...""")
    alpha_integrand=np.zeros((len(alphas),len(k_lins)))
    for ialpha,alpha in enumerate(alphas):
        q = alpha*k_lins
        alpha_integrand[ialpha] = (smooth(q, del2interp(q))
                                     * kernel_func(alpha))
        if absP_kernel:
            alpha_integrand[ialpha] *= (2*pi**2 / q**3)

    alpha_integral = trapz_2d(alphas, alpha_integrand)      
    return alpha_integral / (2*pi)**2



def alpha_mu_integral_mackey(k_lins, del2interp, kernel_func, log_alpha_coarse_min=-5.,log_alpha_coarse_max=5.,n_alpha_coarse=300,
            alpha_fine_min=0.98,alpha_fine_max=1.02,n_alpha_fine=300,mu_fine_min=0.99,n_mu_coarse=300,
            n_mu_fine=300,smooth=no_smooth):
    """This performs integrals of the type in eqn. 19 of 
    Mackey+ 01 (astro-ph/0106364). It doesn't include the prefactor.
    The 'kernel_func' argument expects a function of alpha and mu depending
    on the IA term being calculated. For example of the tidal torque II term
    use the hehe function.
    The integral is performed in two steps in both alpha and mu - with
    finer binning when either variable is close to 1 (since here the denominator
    goes to zero eqn. 19). 
    The function checks that the specified ranges are such that both 
    can't be exactly 1.

    @param k_lins           The values of k at which to evaluate the integral
    @param del2interp       A Del2kInterpFunc object
    @param kernel_func      An IA kernel function
    @param smooth           Optional function to smooth functions of k. [default: No smoothing]
    @returns the integral obvs.
    """
    mu_min,mu_max=-1.,1.
    alphas=np.sort(np.concatenate((np.logspace(log_alpha_coarse_min,log_alpha_coarse_max,n_alpha_coarse),
                                   np.linspace(alpha_fine_min,alpha_fine_max,n_alpha_fine))))
    mus=np.concatenate((np.linspace(mu_min,mu_fine_min,n_mu_coarse),
                        np.linspace(mu_fine_min,mu_max,n_mu_fine)))
    #Do some checks on integral ranges. 
    #Check if alpha and mu both get very close to 1.
    if np.any(np.isclose(alphas,1,atol=1e-9,rtol=0.)):
        raise ValueError("""alpha gets rather close to 1....
                        choose different alpha values...""")
    #Check mus are sensible
    assert (mu_fine_min>=-1.) and (mu_fine_min<=1.)
    assert np.abs(alpha_fine_min)>np.abs(10**log_alpha_coarse_min)
    assert np.abs(alpha_fine_max)<np.abs(10**log_alpha_coarse_max)
    alpha_integrand=np.zeros((len(alphas),len(k_lins)))
    for ialpha,alpha in enumerate(alphas):
        mu_integrand=np.zeros((len(mus),len(k_lins)))
        for imu,mu in enumerate(mus):
            polyee = kernel_func(alpha,mu)
            k_diff_ang=np.sqrt(1 + alpha**2 - 2*alpha*mu)
            ks_mu = k_lins * k_diff_ang
            del2 = del2interp(ks_mu)
            mu_integrand[imu] = (smooth(ks_mu,del2)
                                    * polyee) / k_diff_ang**7
        mu_integral = trapz_2d(mus, mu_integrand)
        alpha_integrand[ialpha] = ((smooth(k_lins, del2interp(alpha * k_lins)) 
                                    / alpha) * mu_integral)
    alpha_integral = trapz_2d(alphas, alpha_integrand)      
    return alpha_integral
