import scipy.interpolate as interp
import numpy as np
import scipy.integrate
import del4
import matplotlib
import pylab
import inspect
import kernels

pi=np.pi

def DEE_term(k_lins, del2_interp, acc=del4.ACC_LOW):
        kernel = kernels.feF2he_D
        return (k_lins**3*del4.alpha_mu_integral(k_lins, del2_interp, kernel, Pq2=True, absP_kernel=True, **acc),
                2 * -4. * del2_interp.sigma2_S/45)

def DBB_term(k_lins, del2_interp, acc=del4.ACC_LOW):
        kernel = kernels.fbF2hb_D
        return (k_lins**3*del4.alpha_mu_integral(k_lins, del2_interp, kernel, Pq2=True, absP_kernel=True, **acc),
                2 * -4. * del2_interp.sigma2_S/45)

def get_all_terms(z_lin, k_lins, P_lin, z_growth, Dz, acc="low", plot=False):

    if acc=="low":
        acc_dict=del4.ACC_LOW
    elif acc=="medium":
        acc_dict=del4.ACC_MED
    elif acc=="high":
        acc_dict=del4.ACC_HI
    else:
        raise KeyError("invalid accuracy setting")

    p_lin_z0=P_lin[0]
    #We need the growth function to rescale the power spectrum at the end...get growth at same z-values as P_lin
    assert np.isclose(z_growth[0],0.)
    assert ((z_growth[-1]>z_lin[-1]) or np.isclose(z_growth[-1],z_lin[-1]))
    if not np.allclose(z_lin, z_growth):
        Dz=interp.interp1d(z_growth,Dz)(z_lin)

    DEE,DEE_const = DEE_term(k_lins, del2_interp, acc=acc_dict)
    DBB,DBB_const = DBB_term(k_lins, del2_interp, acc=acc_dict)
    #print "EE_const, BB_const, DEE_const, DBB_const"
    #print EE_const, BB_const, DEE_const, DBB_const
    if plot:
        import pylab
        fig1=pylab.figure()
        ax=fig1.add_subplot(111)
        ax.plot(k_lins, EE, 'b', label=r'$II_{EE}(k)$')
        ax.plot(k_lins, BB, 'g', label=r'$II_{BB}(k)$')
        ax.plot(k_lins, A, 'c',label=r'$A(k)$')
        ax.plot(k_lins, -A, 'c--')
        ax.plot(k_lins, B_solved, 'gold',label=r'$B(k)$')
        ax.plot(k_lins, -B_solved,'--',color='gold')
        ax.plot(k_lins, DEE, 'm',label=r'$D_{EE}(k)$')
        ax.plot(k_lins, -DEE,'m--')  
        ax.plot(k_lins, DBB, 'r',label=r'$D_{BB}(k)$')
        ax.plot(k_lins, -DBB,'r--')     
        ax.plot(k_lins, EE+EE_const, 'b')
        ax.plot(k_lins, -(EE+EE_const), 'b--')
        ax.plot(k_lins, BB+BB_const, 'g')
        ax.plot(k_lins, -(BB+BB_const), 'g--')
        ax.plot(k_lins, DEE+DEE_const, 'm')
        ax.plot(k_lins, -(DEE+DEE_const),'m--')  
        ax.plot(k_lins, DBB+DBB_const, 'r')
        ax.plot(k_lins, -(DBB+DBB_const),'r--')               
        #ax.plot(k_lins, mackey_EE, 'm-',label='mackey')
        ax.plot(k_lins, p_lin_z0, 'k-', label=r'$P_{\mathrm{lin}}(k,z=0)$')
        ax.set_xlabel(r'$k$')
        ax.set_xlim([1.e-4,20])
        ax.set_ylim([1e-2,1e6])
        ax.set_ylabel(r'$P(k,z=0)$')
        ax.set_xscale('log')
        ax.set_yscale('log')
        pylab.legend(loc=1,ncol=2)
        pylab.tight_layout()
        pylab.savefig('TT.pdf')
        pylab.show()

    specs_2d = {"DEE":DEE+DEE_const,"DBB":DBB+DBB_const}
    specs_3d = {}
    for key,val in specs_2d.iteritems():
        specs_3d[key]=del4.grow_Pk(spec, z_lin, Dz, n=4)
    return specs_3d


