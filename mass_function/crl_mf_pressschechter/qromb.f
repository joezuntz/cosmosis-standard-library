      SUBROUTINE qromb(func,a,b,ss,q) ! one parameter
      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,q,EPS
      EXTERNAL func
      PARAMETER (EPS=1.d-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint,trapzd
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      h(1)=1d0
      do 11 j=1,JMAX
        call trapzd(func,a,b,s(j),j,q)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0d0,ss,dss)
          if (dabs(dss).le.EPS*dabs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=0.25d0*h(j)
11    continue
      print*, 'too many steps in qromb'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software D041&0(9p#3.
      SUBROUTINE trapzd(func,a,b,s,n,q) ! one paramete
      INTEGER n
      REAL*8 a,b,s,func,q
      EXTERNAL func
      INTEGER it,j
      REAL*8 del,sum,tnm,x
      if (n.eq.1) then
        s=0.5d0*(b-a)*(func(a,q)+func(b,q))
      else
        it=2**(n-2)
        tnm=dble(it)
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0d0
        do 11 j=1,it
          sum=sum+func(x,q)
          x=x+del
11      continue
        s=0.5d0*(s+(b-a)*sum/tnm)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software D041&0(9p#3.
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=50)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=dabs(x-xa(1))
      do 11 i=1,n
        dift=dabs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0d0)print*, 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software D041&0(9p#3.
