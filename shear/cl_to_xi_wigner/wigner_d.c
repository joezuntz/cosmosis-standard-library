// wigner_d.c -- compute the Wigner d-matrix recursively in angular momentum
// ----
// author: Nicolas Tessore <nicolas.tessore@manchester.ac.uk>
// date: 31 Aug 2018
// ----
// uses SSE intrinsics by default if SSE3 is detected; disabled with -DNO_SSE

#include <stdlib.h>
#include <math.h>

#ifndef NO_SSE
#ifdef __SSE3__
#include <x86intrin.h>
#define USE_SSE
#endif
#endif

static inline int binom(int n, int k)
{
    int b, i;
    
    if(k > n/2)
        k = n-k;
    
    b = 1;
    for(i = 1; i <= k; ++i, --n)
    {
        if(n%i == 0)
            b = (n/i)*b;
        else if(b%i == 0)
            b = (b/i)*n;
        else
            b = (n*b)/i;
    }
    return b;
}

void legendre_p(double x, int n0, int n1, double* p)
{
    int n;
    double p0, p1 = 1, p2 = x;
    if(0 >= n0)
        *(p++) = 1;
    if(1 >= n0)
        *(p++) = x;
    for(n = 2; n <= n1; ++n)
    {
        p0 = p1;
        p1 = p2;
        p2 = ((2*n-1)*x*p1 - (n-1)*p0)/n;
        if(n >= n0)
            *(p++) = p2;
    }
}

void wigner_d(int l0, int l1, int n, int m, double theta, double* d)
{
    double d0, u, v, x;
    int l, lp, a, b, c;
    
#ifndef USE_SSE
    double d1, d2, j;
#else
    __m128d o, j, z, r, s, t;
#endif
    
    if(n == 0 && m == 0)
    {
        legendre_p(cos(theta), l0, l1, d);
        return;
    }
    
    if(abs(n) > abs(m))
    {
        if(n > 0)
            lp = n, a = n - m, b = n + m, c = n - m;
        else
            lp = -n, a = m - n, b = -n - m, c = 0;
    }
    else
    {
        if(m > 0)
            lp = m, a = m - n, b = n + m, c = 0;
        else
            lp = -m, a = n - m, b = -n - m, c = n - m;
    }
    
    u = sin(0.5*theta);
    v = cos(0.5*theta);
    x = v*v - u*u;
    
    d0 = (1 - 2*(c&1))*sqrt(binom(a+b, a))*pow(u, a)*pow(v, b);
    
#ifndef USE_SSE
    d1 = 0;
    j = n*m;
#else
    t = _mm_set_pd(0, d0);
    o = _mm_set1_pd(1);
    j = _mm_set_pd(n, m);
    z = _mm_set_pd(1, -n*m);
#endif
    
    for(l = l0; l < lp; ++l)
        *(d++) = 0;
    if(l == lp)
        *(d++) = d0;
    for(l = lp+1; l <= l1; ++l)
    {
#ifndef USE_SSE
        u = (1.-1./(l-n))*(1.-1./(l+n));
        v = (1.-1./(l-m))*(1.-1./(l+m));
        
        d2 = d1;
        d1 = d0;
        d0 = (l*x-j/(l-1))*sqrt((1-u)*(1-v))*d1 - (1.+1./(l-1))*sqrt(u*v)*d2;
        
        if(l >= l0)
            *(d++) = d0;
#else
        r = _mm_set1_pd(l);
        s = _mm_sub_pd(o, _mm_div_pd(o, _mm_add_pd(r, j)));
        r = _mm_sub_pd(o, _mm_div_pd(o, _mm_sub_pd(r, j)));
        s = _mm_mul_pd(r, s);
        r = _mm_sub_pd(o, s);
        s = _mm_mul_pd(_mm_unpacklo_pd(r, s), _mm_unpackhi_pd(r, s));
        
        r = _mm_add_pd(_mm_set_pd(1, l*x), _mm_div_pd(z, _mm_set1_pd(l-1)));
        r = _mm_mul_pd(_mm_mul_pd(r, _mm_sqrt_pd(s)), t);
        r = _mm_hsub_pd(r, r);
        t = _mm_unpacklo_pd(r, t);
        
        if(l >= l0)
            _mm_store_sd(d++, t);
#endif
    }
}
