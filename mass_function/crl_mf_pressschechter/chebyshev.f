      SUBROUTINE chebft(a,b,c,n,func)
      INTEGER n,NMAX
      REAL a,b,c(n),func,PI
      EXTERNAL func
      PARAMETER (NMAX=50, PI=3.141592653589793d0)
      INTEGER j,k
      REAL bma,bpa,fac,y,f(NMAX)
      DOUBLE PRECISION sum
      bma=0.5*(b-a)
      bpa=0.5*(b+a)
      do 11 k=1,n
        y=cos(PI*(k-0.5)/n)
        f(k)=func(y*bma+bpa)
11    continue
      fac=2./n
      do 13 j=1,n
        sum=0.d0
        do 12 k=1,n
          sum=sum+f(k)*cos((PI*(j-1))*((k-0.5d0)/n))
12      continue
        c(j)=fac*sum
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software D041&0(9p#3.
      SUBROUTINE chder(a,b,c,cder,n)
      INTEGER n
      REAL a,b,c(n),cder(n)
      INTEGER j
      REAL con
      cder(n)=0.
      cder(n-1)=2*(n-1)*c(n)
      if(n.ge.3)then
        do 11 j=n-2,1,-1
          cder(j)=cder(j+2)+2*j*c(j+1)
11      continue
      endif
      con=2./(b-a)
      do 12 j=1,n
        cder(j)=cder(j)*con
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software D041&0(9p#3.
      FUNCTION chebev(a,b,c,m,x)
      INTEGER m
      REAL chebev,a,b,x,c(m)
      INTEGER j
      REAL d,dd,sv,y,y2
      if ((x-a)*(x-b).gt.0.) then
         write(*,*) a,b,x
         print*, 'x not in range in chebev'
      endif
      d=0.
      dd=0.
      y=(2.*x-a-b)/(b-a)
      y2=2.*y
      do 11 j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
11    continue
      chebev=y*d-dd+0.5*c(1)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software D041&0(9p#3.
