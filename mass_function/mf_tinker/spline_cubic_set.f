      SUBROUTINE spline_cubic_set(n,t,y,ibcbeg,ybcbeg,ibcend,ybcend,ypp)
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!  Modified:
!    06 June 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
      implicit none

      integer n

      double precision a1(n)
      double precision a2(n)
      double precision a3(n)
      double precision a4(n)
      double precision a5(n)
      double precision b(n)
      integer i
      integer ibcbeg
      integer ibcend
      double precision t(n)
      double precision y(n)
      double precision ybcbeg
      double precision ybcend
      double precision ypp(n)
c
c  Check.
c
      if ( n .le. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
        write ( *, '(a)' ) '  The number of knots must be at least 2.'
        write ( *, '(a,i8)' ) '  The input value of N = ', n
        stop
      end if

      do i = 1, n - 1
        if ( t(i+1) .le. t(i) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
          write ( *, '(a)' ) 
     &      '  The knots must be strictly increasing, but'
          write ( *, '(a,i8,a,g14.6)' ) '  T(',  i,') = ', t(i)
          write ( *, '(a,i8,a,g14.6)' ) '  T(',i+1,') = ', t(i+1)
          stop
        end if
      end do
c
c  Zero out the matrix.
c
      do i = 1, n
        a1(i) = 0.0D+00
        a2(i) = 0.0D+00
        a3(i) = 0.0D+00
        a4(i) = 0.0D+00
        a5(i) = 0.0D+00
      end do
c
c  Set the first equation.
c
      if ( ibcbeg .eq. 0 ) then
        b(1) = 0.0D+00
        a3(1) = 1.0D+00
        a4(1) = -1.0D+00
      else if ( ibcbeg .eq. 1 ) then
        b(1) = ( y(2) - y(1) ) / ( t(2) - t(1) ) - ybcbeg
        a3(1) = ( t(2) - t(1) ) / 3.0D+00
        a4(1) = ( t(2) - t(1) ) / 6.0D+00
      else if ( ibcbeg .eq. 2 ) then
        b(1) = ybcbeg
        a3(1) = 1.0D+00
        a4(1) = 0.0D+00
      else if ( ibcbeg .eq. 3 ) then
        b(1) = 0.0D+00
        a3(1) = - ( t(3) - t(2) )
        a4(1) =   ( t(3)        - t(1) )
        a5(1) = - (        t(2) - t(1) )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
        write ( *, '(a)' ) 
     &    '  The boundary flag IBCBEG must be 0, 1, 2, or 3.'
        write ( *, '(a,i8)' ) '  The input value is IBCBEG = ', ibcbeg
        stop
      end if
c
c  Set the intermediate equations.
c
      do i = 2, n - 1
        b(i) = ( y(i+1) - y(i) ) / ( t(i+1) - t(i) ) 
     &       - ( y(i) - y(i-1) ) / ( t(i) - t(i-1) )
        a2(i) = ( t(i+1) - t(i)   ) / 6.0D+00
        a3(i) = ( t(i+1) - t(i-1) ) / 3.0D+00
        a4(i) = ( t(i)   - t(i-1) ) / 6.0D+00
      end do
c
c  Set the last equation.
c
      if ( ibcend .eq. 0 ) then
        b(n) = 0.0D+00
        a2(n) = -1.0D+00
        a3(n) = 1.0D+00
      else if ( ibcend .eq. 1 ) then
        b(n) = ybcend - ( y(n) - y(n-1) ) / ( t(n) - t(n-1) )
        a2(n) = ( t(n) - t(n-1) ) / 6.0D+00
        a3(n) = ( t(n) - t(n-1) ) / 3.0D+00
      else if ( ibcend .eq. 2 ) then
        b(n) = ybcend
        a2(n) = 0.0D+00
        a3(n) = 1.0D+00
      else if ( ibcend .eq. 3 ) then
        b(n) = 0.0D+00
        a1(n) = - ( t(n) - t(n-1) )
        a2(n) =   ( t(n)          - t(n-2) )
        a3(n) = - (        t(n-1) - t(n-2) )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
        write ( *, '(a)' ) 
     &    '  The boundary flag IBCEND must be 0, 1, 2, or 3.'
        write ( *, '(a,i8)' ) '  The input value is IBCEND = ', ibcend
        stop
      end if
c
c  Special case:
c    N = 2, IBCBEG = IBCEND = 0.
c
      if ( n .eq. 2 .and. ibcbeg .eq. 0 .and. ibcend .eq. 0 ) then

        ypp(1) = 0.0D+00
        ypp(2) = 0.0D+00
c
c  Solve the linear system.
c
      else

        call penta ( n, a1, a2, a3, a4, a5, b, ypp )

      end if

      return
      END 
      SUBROUTINE penta ( n, a1, a2, a3, a4, a5, b, x )

c*********************************************************************72
c
cc PENTA solves a pentadiagonal system of linear equations.
c
c  Discussion:
c
c    The matrix A is pentadiagonal.  It is entirely zero, except for
c    the main diagaonal, and the two immediate sub- and super-diagonals.
c
c    The entries of Row I are stored as:
c
c      A(I,I-2) -> A1(I)
c      A(I,I-1) -> A2(I)
c      A(I,I)   -> A3(I)
c      A(I,I+1) -> A4(I)
c      A(I,I-2) -> A5(I)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 June 2013
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Cheney, Kincaid,
c    Numerical Mathematics and Computing,
c    1985, pages 233-236.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A1(N), A2(N), A3(N), A4(N), A5(N), the
c nonzero
c    elements of the matrix.  Note that the data in A2, A3 and A4
c    is overwritten by this routine during the solution process.
c
c    Input, double precision B(N), the right hand side of the linear
c system.
c
c    Output, double precision X(N), the solution of the linear system.
c
      implicit none

      integer n

      double precision a1(n)
      double precision a2(n)
      double precision a3(n)
      double precision a4(n)
      double precision a5(n)
      double precision b(n)
      integer i
      double precision x(n)
      double precision xmult

      do i = 2, n - 1
        xmult = a2(i) / a3(i-1)
        a3(i) = a3(i) - xmult * a4(i-1)
        a4(i) = a4(i) - xmult * a5(i-1)
        b(i) = b(i) - xmult * b(i-1)
        xmult = a1(i+1) / a3(i-1)
        a2(i+1) = a2(i+1) - xmult * a4(i-1)
        a3(i+1) = a3(i+1) - xmult * a5(i-1)
        b(i+1) = b(i+1) - xmult * b(i-1)
      end do

      xmult = a2(n) / a3(n-1)
      a3(n) = a3(n) - xmult * a4(n-1)
      x(n) = ( b(n) - xmult * b(n-1)) / a3(n)
      x(n-1) = ( b(n-1) - a4(n-1) * x(n)) / a3(n-1)
      do i = n - 2, 1, -1
        x(i) = ( b(i)- a4(i) * x(i+1) - a5(i) * x(i+2) ) / a3(i)
      end do

      return
      END
      SUBROUTINE spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, 
     &  yppval )

c*********************************************************************72
c
cc SPLINE_CUBIC_VAL evaluates a piecewise cubic spline at a point.
c
c  Discussion:
c
c    SPLINE_CUBIC_SET must have already been called to define the
c    values of YPP.
c
c    For any point T in the interval T(IVAL), T(IVAL+1), the form of
c    the spline is
c
c      SPL(T) = A
c             + B * ( T - T(IVAL) )
c             + C * ( T - T(IVAL) )^2
c             + D * ( T - T(IVAL) )^3
c
c    Here:
c      A = Y(IVAL)
c      B = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
c        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
c      C = YPP(IVAL) / 2
c      D = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 November 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carl deBoor,
c    A Practical Guide to Splines,
c    Springer, 2001,
c    ISBN: 0387953663.
c
c  Parameters:
c
c    Input, integer N, the number of data values.
c
c    Input, double precision T(N), the knot values.
c
c    Input, double precision Y(N), the data values at the knots.
c
c    Input, double precision YPP(N), the second derivatives of the
c    spline at the knots.
c
c    Input, double precision TVAL, a point, typically between T(1) and
c    T(N), at which the spline is to be evalulated.  If TVAL lies
c utside
c    this range, extrapolation is used.
c
c    Output, double precision YVAL, YPVAL, YPPVAL, the value of the
c spline, and
c    its first two derivatives at TVAL.
c
      implicit none

      integer n

      double precision dt
      double precision h
      integer left
      integer right
      double precision t(n)
      double precision tval
      double precision y(n)
      double precision ypp(n)
      double precision yppval
      double precision ypval
      double precision yval
c
c  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
c  Values below T(1) or above T(N) use extrapolation.
c
      call r8vec_bracket ( n, t, tval, left, right )
c
c  Evaluate the polynomial.
c
      dt = tval - t(left)
      h = t(right) - t(left)

      yval = y(left) 
     &     + dt * ( ( y(right) - y(left) ) / h 
     &            - ( ypp(right) / 6.0D+00 + ypp(left) / 3.0D+00 ) * h 
     &     + dt * ( 0.5D+00 * ypp(left) 
     &     + dt * ( ( ypp(right) - ypp(left) ) / ( 6.0D+00 * h ) ) ) )

      ypval = ( y(right) - y(left) ) / h 
     &     - ( ypp(right) / 6.0D+00 + ypp(left) / 3.0D+00 ) * h 
     &     + dt * ( ypp(left) 
     &     + dt * ( 0.5D+00 * ( ypp(right) - ypp(left) ) / h ) )

      yppval = ypp(left) + dt * ( ypp(right) - ypp(left) ) / h

      return
      END
      SUBROUTINE r8vec_bracket ( n, x, xval, left, right )

c*********************************************************************72
c
cc R8VEC_BRACKET searches a sorted array for successive brackets of a
cvalue.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    If the values in the vector are thought of as defining intervals
c    on the real line, then this routine searches for the interval
c    nearest to or containing the given value.
c
c  Modified:
c
c    24 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, length of input array.
c
c    Input, double precision X(N), an array that has been sorted into
c    ascending order.
c
c    Input, double precision XVAL, a value to be bracketed.
c
c    Output, integer LEFT, RIGHT, the results of the search.
c    Either:
c      XVAL < X(1), when LEFT = 1, RIGHT = 2;
c      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
c    or
c      X(LEFT) <= XVAL <= X(RIGHT).
c
      implicit none

      integer n

      integer i
      integer left
      integer right
      double precision x(n)
      double precision xval

      do i = 2, n - 1

        if ( xval .lt. x(i) ) then
          left = i - 1
          right = i
          return
        end if

       end do

      left = n - 1
      right = n

      return
      END
