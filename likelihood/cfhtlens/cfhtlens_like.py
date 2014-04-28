import numpy as np
import os
import scipy.interpolate as interpolate

n_z_bin = 6
n_theta = 5
n_z_pair = n_z_bin*(n_z_bin+1)/2
n_tot = n_z_pair*n_theta*2

dirname = os.path.split(__file__)[0]
DEFAULT_COVMAT = os.path.join(dirname,'covariance_matrix.dat')
DEFAULT_DATA = os.path.join(dirname,'cfhtlens_xipm_6bin.dat')


def load_covariance(filename, plus_only, cut_low_theta):
    data = np.loadtxt(filename).T
    assert data.shape[0] == data.shape[1] == n_tot
    if plus_only:
        #remove all except top-left of matrix (xiplus-xiplus)
        n_plus = n_tot/2
        data = data[:n_plus,:n_plus]
    if cut_low_theta:
        # delete every 5th row and column
        to_cut = slice(0,data.shape[0],n_theta)
        #rows:
        data = np.delete(data, to_cut, axis=0)
        #columns:
        data = np.delete(data, to_cut, axis=1)
    return data

def load_data_vector(filename, plus_only, cut_low_theta):
    data = np.loadtxt(filename).T
    if cut_low_theta:
        #cut 0th row and 5th row, the theta=smallest one
        data = np.delete(data, n_theta, axis=1)
        data = np.delete(data, 0, axis=1)
        xi = data[1:]
    if plus_only:
        #cut out first half of data
        data = data[:,:data.shape[1]/2]
        theta = data[0]
        xi = data[1:]
    else:
        theta = data[0,:data.shape[1]/2]
        xi_plus = data[1:,:data.shape[1]/2]
        xi_minus = data[1:,data.shape[1]/2:]
        xi_plus_vector = np.concatenate(xi_plus)
        xi_minus_vector = np.concatenate(xi_minus)
        xi = np.concatenate((xi_plus_vector,xi_minus_vector))

    
    return theta, xi.flatten()

def theta_lims_from_bin_mids(theta):
    #Assume theta mids are equally spaced in log space, and that
    #log10(theta_mid) = 0.5 * (log10(theta_lower)+log10(theta_upper))
    #delta_log_theta is difference in log space between bin mids - this
    #will also be the difference between bin limits in log space
    delta_log_theta=np.log10(theta[1])-np.log10(theta[0])                                                           
    log10_theta_lims=np.zeros(len(theta)+1)

    for i in range(len(theta)+1):
        log10_theta_lims[i]=np.log10(theta[0])+(i-0.5)*delta_log_theta
    return 10**log10_theta_lims

class CFHTLensLikelihood(object):
    def __init__(self, covmat_filename, data_filename,plus_only,cut_low_theta):
        self.plus_only=plus_only
        self.cut_low_theta=cut_low_theta
        theta, self.data = load_data_vector(data_filename, self.plus_only, self.cut_low_theta)
        #Convert thetas from arcminutes to radians?? to match theory thetas
        self.theta = np.radians(theta/60.)
        #print 'theta bin mids',self.theta
        covmat = load_covariance(covmat_filename, self.plus_only, self.cut_low_theta)
        weight = np.linalg.inv(covmat)
        n_mu = 1656
        if self.plus_only:
            if self.cut_low_theta: p=84
            else: p=105
        else:
            if self.cut_low_theta: p=168
            else: p=210

        alpha = (n_mu - p - 2.0) / (n_mu - 1.0)
        print "Need to check the Anderson Hartlap when cutting matrix - cut first?"
        print "xi_plus only?",self.plus_only
        print "Cut low thetas?",self.cut_low_theta
        weight *= alpha
        self.weight = weight


    def interpolate_to_bin(self, theta, xi):
        return np.interp(self.theta, theta, xi)
    
    """Don't need this now! We have the correct theta values to interpolate at"""
    def integrate_theta_across_bins(self, theta_theory, xi_theory):
        #First get theta_bin_lims from thetas...(or self eventually)
        theta_bin_limits=theta_lims_from_bin_mids(self.theta)
        #Check there's the right number
        assert len(theta_bin_limits)==len(self.theta)+1
        #To get predicted xi for a bin, need to integrate xi * weight function
        #between lower and upper limit of bin (theta_lo,theta_hi)
        #i.e. xi_bin = \int_{theta_lo}^{theta_hi} W(theta) * xi(theta) /  [\int_{theta_lo}^{theta_hi} W(theta)]
        #Geometrically we expect W(theta) = theta^2 (is there some extra effect due to clustering?)
        #First interpolate theory xi - we will use this for all the integrals
        #In fact interpolate xi*theta^2 since this is what we're integrating...use scipy.interplate.splrep (cubic spline)
        integrand_interp = interpolate.splrep(theta_theory, xi_theory*theta_theory**2)
        #Normalization is just \int dtheta theta^2 = (1/3)*(theta_hi^3 - theta_lo^3)
        norms = (theta_bin_limits[1:]**3 - theta_bin_limits[:-1]**3)/3. 
        xi_binned=[]
        for i in range(len(self.theta)):
            theta_lo,theta_hi=theta_bin_limits[i],theta_bin_limits[i+1]
            #Use scipy.interpolate.splint to integrate spline integral
            integrand = interpolate.splint(theta_lo,theta_hi,integrand_interp)
            xi_binned.append(integrand/norms[i])
        return np.array(xi_binned)

    def __call__(self, xi_data):

        if self.plus_only:
            #xi_data is a dictionary containing all the bin pair data vectors
            xi_vector = []
            #loop through the bins loading the theory data
            for i in xrange(1, n_z_bin+1):
                for j in xrange(i, n_z_bin+1):
                    #try both orderings for flexibility
                    try:
                        theta_theory, xi_theory = xi_data[(i,j)]
                    except KeyError:
                        theta_theory, xi_theory = xi_data[(j,i)]
                    theta_theory=theta_theory[:len(theta_theory)/2]
                    xi_theory = xi_theory[:len(theta_theory)]
                    #interpolate to data values - this is correct if self.theta values are mean for each bin
                    xi_binned = self.interpolate_to_bin(theta_theory, xi_theory)
                    #build up longer vector of data
                    xi_vector.append(xi_binned)
            #flatten to single vector
            xi_vector = np.concatenate(xi_vector)
            assert xi_vector.shape == self.data.shape
            #get chi2 and return log like
            delta = xi_vector - self.data
            chi2 = np.dot(delta, np.dot(self.weight, delta))
            return -chi2/2.0

        else:
            #xi_data is a dictionary containing all the bin pair data vectors
            xi_plus_vector = []
            xi_minus_vector = []
            #loop through the bins loading the theory data
            for i in xrange(1, n_z_bin+1):
                for j in xrange(i, n_z_bin+1):
                    #try both orderings for flexibility
                    try:
                        theta_theory, xi_theory = xi_data[(i,j)]
                    except KeyError:
                        theta_theory, xi_theory = xi_data[(j,i)]
                    theta_theory=theta_theory[:len(theta_theory)/2]
                    xi_plus_theory = xi_theory[:len(theta_theory)]
                    xi_minus_theory = xi_theory[len(theta_theory):]
                    #print xi_plus_theory,xi_minus_theory
                    #interpolate to data values - this is correct if self.theta values are mean for each binp the
                    xi_plus_binned = self.interpolate_to_bin(theta_theory, xi_plus_theory) 
                    xi_minus_binned = self.interpolate_to_bin(theta_theory, xi_minus_theory)
                    '''
                    pylab.subplot(121)
                    pylab.loglog(theta_theory,xi_plus_theory)
                    pylab.loglog(self.theta,xi_plus_binned,'o')
                    pylab.subplot(122)
                    pylab.loglog(theta_theory,xi_minus_theory)
                    pylab.loglog(self.theta,xi_minus_binned,'o')
                    pylab.show()
                    '''

                    #build up longer vector of data
                    xi_plus_vector.append(xi_plus_binned)
                    xi_minus_vector.append(xi_minus_binned)

            #flatten to single vector
            xi_plus_vector = np.concatenate(xi_plus_vector)
            xi_minus_vector = np.concatenate(xi_minus_vector)
            #print 'len(xi_plus_vector)', len(xi_plus_vector)
            #concatenate xi_plus_vector and xi_minus_vector to (hopefully) match form of data vector 
            xi_vector=np.concatenate((xi_plus_vector,xi_minus_vector))
            #print 'len(xi_vector)', len(xi_vector)
            assert xi_vector.shape == self.data.shape
            #get chi2 and return log like
            delta = xi_vector - self.data
            chi2 = np.dot(delta, np.dot(self.weight, delta))
            return -chi2/2.0


# ordering is (1,1) (1,2) (1,3) (1,4) (1,5) (1,6) (2,2) (2,3) ... (5,5), (5,6), (6,6)
'''
import pylab
covmat = load_covariance(DEFAULT_COVMAT,False,False)
theta, xi = load_data_vector(DEFAULT_DATA,False,False)
theta = theta[:len(theta)/2]
n_theta_theory=len(theta)
xi_plus,xi_minus = xi[:n_theta_theory],xi[n_theta_theory:]
sigma = covmat.diagonal()**0.5
sigma_plus,sigma_minus = sigma[:n_theta_theory],sigma[n_theta_theory:]
# print theta
# print xi[:4]
pylab.errorbar(theta, xi_plus[0:5], sigma_plus[0:5], fmt='.')
pylab.errorbar(theta, xi_minus[0:5], sigma_minus[0:5], fmt='.')
# pylab.gca().set_yscale('log',nonposy='clip')
# pylab.gca().set_xscale('log')
pylab.show()
'''