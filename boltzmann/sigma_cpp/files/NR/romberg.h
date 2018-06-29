# include "quadrature.h"
# include "nr3.h"



template <class T>
double simpinteg (T &func, Doub kmin, Doub kmax, Doub eps=1e-6){
			double k,wkR,I=0,f,dk=(kmax-kmin)*eps;
			for(k=kmin;k<=kmax;k=k+dk){
                            			    f=func(k);
                            			    I=I+f*dk;
                            			    }

			return I;
			}

template <class T>
Doub qromb(T &func, Doub a, Doub b, const Doub eps=1.0e-10) {
	const Int JMAX=40, JMAXP=JMAX+1, K=5;//printf("a=%lf b=%lf\n",a,b);
	VecDoub s(JMAX),h(JMAXP);
	Poly_interp polint(h,s,K);
	h[0]=1.0;
	Trapzd<T> t(func,a,b);//printf("sucsess3.1\n");
	for (Int j=1;j<=JMAX;j++) {
		s[j-1]=t.next();
		if (j >= K) {
			Doub ss=polint.rawinterp(j-K,0.0);//printf("qromb %lf %lf %d %d %lf\n",a,b,j,j-K,ss);
			if (abs(polint.dy) <= eps*abs(ss)) return ss;
		}
		h[j]=0.25*h[j-1];
	}
	throw("Too many steps in routine qromb");
}
template<class T>
Doub qromo(Midpnt<T> &q, const Doub eps=3.0e-9) {
	const Int JMAX=40, JMAXP=JMAX+1, K=5;
	VecDoub h(JMAXP),s(JMAX);
	Poly_interp polint(h,s,K);
	h[0]=1.0;
	for (Int j=1;j<=JMAX;j++) {
		s[j-1]=q.next();
		if (j >= K) {
			Doub ss=polint.rawinterp(j-K,0.0);
			if (abs(polint.dy) <= eps*abs(ss)) return ss;
		}
		h[j]=h[j-1]/9.0;
	}
	throw("Too many steps in routine qromo");
}
