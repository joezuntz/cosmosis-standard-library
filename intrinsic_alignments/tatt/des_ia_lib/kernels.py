#IA 
def X1(alpha,mu):
    return alpha**2 - 2*alpha*mu + 1
def X2(alpha,mu):
    return -3*(alpha**2+1)*mu**2 + alpha**2 + 3*alpha*mu**3 + alpha*mu + 1
def X3(alpha,mu):
    return alpha*(10*mu**2 -3) - 7*mu
def X4(alpha,mu):
    return -10*alpha + 7*(1 + alpha**2)*mu - 4*alpha*mu**2
def X5(alpha,mu):
    a = alpha**3*(-38*mu**5 + 4*mu**3 + 2*mu)
    b = alpha**2*(19*mu**6 + 44*mu**4 - 17*mu**2 + 2)
    c = alpha*(-34*mu**5 - 4*mu**3 + 6*mu)
    return a + b + c
def X6(alpha,mu):
    return 19*mu**4 - 14*mu**2 + 3
def X7(alpha,mu):
    return (1-mu**2)*(alpha-mu)**2*(1-2*alpha*mu)**2
def X8(alpha,mu):
    a = -4 + 12*mu**2 + 4*alpha*mu*(1 - 9*mu**2)
    b = alpha**2*(-3 + 10*mu**2 + 41*mu**4)
    c = alpha**3*mu*(9 - 22*mu**2 - 19*mu**4)
    return a + b + c

#Tidal Torque kernels
def hehe(alpha, mu):
    """This is \dphi/2pi h_e^2(q,k-q)/ 2pi.
    """
    return 2*alpha**2*(X5(alpha,mu) + (1+alpha**4)*X6(alpha,mu))/72/X1(alpha,mu)**2

def hbhb(alpha, mu):
    """This is \dphi/2pi h_b^2(q,k-q)/ 2pi.
    """
    return 2*alpha**2*X7(alpha,mu)/18/X1(alpha,mu)**2

#GI
def F2hE_A(alpha,mu):
    """This is \dphi/2pi F_2(q,k-q) * h_e(q,k-q) / 2pi"""
    return 2*alpha*X2(alpha,mu)*X3(alpha,mu)/84/X1(alpha,mu)**2

def F2hE_B(alpha,mu):
    """This is \dphi/2pi h_e(q-k,-q) * F_2(q,-k) / 2pi"""
    return 2*alpha**2*(X2(alpha,mu)*X4(alpha,mu)/84/alpha/X1(alpha,mu) - 29./630)

def F2he_B_mu_solved_all_alpha(alpha):
    a = 2*alpha*(225 - 600*alpha**2 + 1198*alpha**4 - 600*alpha**6 + 225*alpha**8)
    b = 225*(alpha**2 - 1)**4*(alpha**2 + 1)*np.log(np.abs((alpha-1)/(alpha+1)))
    c = 20160*alpha**5
    z_alpha = (a + b)/c
    return 2*alpha**2 * (z_alpha - 29./315)

def F2he_B_mu_solved_large_alpha(alpha):
    return -16./147 + 353./4704/alpha**2

def F2he_B_mu_solved(alpha, alpha_thresh=100):
    if alpha>alpha_thresh:
        return F2he_B_mu_solved_large_alpha(alpha)
    else:
        return F2he_B_mu_solved_all_alpha(alpha)

#Mixed IA kernels
def feF2he_D(alpha,mu):
    """This \int dphi/2pi f_e(k-q) F_2(q,k-q) h_e(q,k-q)"""
    return 2*alpha**2*(X8(alpha,mu)+alpha**4*X6(alpha,mu))/24/X1(alpha,mu)**2

def fbF2hb_D(alpha,mu):
    """This \int dphi/2pi f_b(k-q) F_2(q,k-q) h_b(q,k-q)"""
    return 2*alpha**2*(alpha*(mu-alpha)*(alpha*mu-1)*
            (2*alpha*mu-1)*(mu**2-1))/3/X1(alpha,mu)**2

#Tidal alignment kernels
#GI
def P_2(mu):
    return 0.5*(3*mu*mu - 1)

def F2q2q(alpha, mu):
    """with 2*alpha**2 P_2 prefactor"""
    return 2*alpha**2 * P_2(mu) * (7*mu + 3*alpha - 10*alpha*mu**2)/14/alpha/(1 + alpha**2 - 2*alpha*mu)

def F2minuskq2(alpha,mu):
    """with 2*alpha**2 P_2 prefactor"""
    return 2*alpha**2 * P_2(mu) * alpha**2*(3 + 7*alpha*mu - 10*mu**2)/(1 + alpha**2 - 2*alpha*mu)/14

#II
def fEfE_A(alpha,mu):
    """with alpha**2 prefactor"""
    return alpha**2 * (3- 14*mu**2 + 19*mu**4)/8.

def fEfE_B(alpha,mu):
    a = 8*alpha*mu*(3*mu**2 - 1)
    b = 4*(3*mu**2 - 1)
    c = alpha**2 * (3 - 14*mu**2 + 19*mu**4)
    d = 8*(1 + alpha**2 - 2*alpha*mu)
    return (a+b+c)/d


#old
def g_ee(alpha, mu):
    """the 'g_ee' polynomial from Mackey+01
    """
    a= 3-14*mu**2+19*mu**4
    b= 3*mu-2*mu**3-17*mu**5
    c= 2-17*mu**2+44*mu**4+19*mu**6
    d= mu+2*mu**3-19*mu**5
    return 0.125 * ((1+alpha**4)*a + 2*alpha*b + (alpha**2) * c + 2*(alpha**3)*d)
