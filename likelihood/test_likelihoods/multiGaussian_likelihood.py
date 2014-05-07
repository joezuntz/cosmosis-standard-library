import os
import pydesglue
import numpy as np
import os.path
try:
    import pylab
except ImportError:
    pylab=None


class Multigaussian(object):
    def __init__(self,ndim):
        self.ndim = ndim
        if os.path.isfile('./means10D.out'):
            self.means=np.loadtxt('./means10D.out')
            #print "loading saved means"
        # else:
        #     self.means = np.random.rand(ndim)
        #     np.savetxt('./means10D.out', self.means)
        
        if os.path.isfile('./covariance10D.out'):
            cov=np.loadtxt('./covariance10D.out')
            #print "loading saved covariance matrix"

        # else:
        #     cov = 0.5-np.random.rand(ndim**2).reshape((ndim, ndim))
        #     cov = np.triu(cov)
        #     cov += cov.T - np.diag(cov.diagonal())
        #     cov = np.dot(cov,cov)
        #     np.savetxt('covariance10D.out', cov)

        self.icov = np.linalg.inv(cov)
        self.det=np.linalg.det(cov)
        self.norm = 0.5*np.log((2.0*np.pi)**self.ndim * self.det)

    def lnprob(self,x):

        diff = x-self.means

        return -np.dot(diff,np.dot(self.icov,diff))/2.0 - self.norm

    def __call__(self,x):
        return lnprob(self,x)


def execute(handle):
    try:
        package = pydesglue.DesDataPackage.from_fits_handle(handle)
        section = pydesglue.section_names.cosmological_parameters
        x=[]
        x.append(package.get_param(section,"MG1"))
        x.append(package.get_param(section,"MG2"))
        x.append(package.get_param(section,"MG3"))
        x.append(package.get_param(section,"MG4"))
        x.append(package.get_param(section,"MG5"))
        x.append(package.get_param(section,"MG6"))
        x.append(package.get_param(section,"MG7"))
        x.append(package.get_param(section,"MG8"))
        x.append(package.get_param(section,"MG9"))
        x.append(package.get_param(section,"MG10"))
        mg = Multigaussian(len(x))
        #compute likelihood
        like =mg.lnprob(x)
    
        section = pydesglue.section_names.likelihoods
        package.set_param(section,"MULTIGAUSSIAN_LIKE",like)
        package.write_to_fits_handle(handle)
    
    except KeyboardInterrupt:
        #If the user presses Ctrl-C we respect that
        raise KeyboardInterrupt
    except Exception as E:
        #But any other kind of error is a problem
        print "There was an error calculating the MultiGaussian likelihood: "
        return 1
        
        return 0

