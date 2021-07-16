c  To get r1mach, mail netlib
c       send r1mach from core
      real function seroin(ax,bx,f,tol)
      real ax,bx,f,tol
c
c      a zero of the function  f(x)  is computed in the interval ax,bx .
c
c  input..
c
c  ax     left endpoint of initial interval
c  bx     right endpoint of initial interval
c  f      function subprogram which evaluates f(x) for any x in
c         the interval  ax,bx
c  tol    desired length of the interval of uncertainty of the
c         final result (.ge.0.)
c
c  output..
c
c  zeroin abscissa approximating a zero of  f  in the interval ax,bx
c
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  this is checked, and an error message is printed if this is not
c  satisfied.   seroin  returns a zero  x  in the given interval
c  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is
c  the  relative machine precision defined as the smallest representable
c  number such that  1.+macheps .gt. 1.
c      this function subprogram is a slightly  modified  translation  of
c  the algol 60 procedure  zero  given in  richard brent, algorithms for
c  minimization without derivatives, prentice-hall, inc. (1973).
c
      real  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
      real   abs, r1mach
   10 eps = r1mach(4)
      tol1 = eps+1.0e0
c
      a=ax
      b=bx
      fa=f(a)
      fb=f(b)
c     check that f(ax) and f(bx) have different signs
      if (fa .eq.0.0e0 .or. fb .eq. 0.0e0) go to 20
      if (fa * (fb/ abs(fb)) .le. 0.0e0) go to 20
         write(6,2500)
2500     format(1x,'f(ax) and f(bx) do not have different signs,',
     1             ' seroin is aborting')
         return
   20 c=a
      fc=fa
      d=b-a
      e=d
   30 if ( abs(fc).ge. abs(fb)) go to 40
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
   40 tol1=2.0e0*eps* abs(b)+0.5e0*tol
      xm = 0.5e0*(c-b)
      if (( abs(xm).le.tol1).or.(fb.eq.0.0e0)) go to 150
c
c see if a bisection is forced
c
      if (( abs(e).ge.tol1).and.( abs(fa).gt. abs(fb))) go to 50
      d=xm
      e=d
      go to 110
   50 s=fb/fa
      if (a.ne.c) go to 60
c
c linear interpolation
c
      p=2.0e0*xm*s
      q=1.0e0-s
      go to 70
c
c inverse quadratic interpolation
c
   60 q=fa/fc
      r=fb/fc
      p=s*(2.0e0*xm*q*(q-r)-(b-a)*(r-1.0e0))
      q=(q-1.0e0)*(r-1.0e0)*(s-1.0e0)
   70 if (p.le.0.0e0) go to 80
      q=-q
      go to 90
   80 p=-p
   90 s=e
      e=d
      if (((2.0e0*p).ge.(3.0e0*xm*q- abs(tol1*q))).or.(p.ge.
     * abs(0.5e0*s*q))) go to 100
      d=p/q
      go to 110
  100 d=xm
      e=d
  110 a=b
      fa=fb
      if ( abs(d).le.tol1) go to 120
      b=b+d
      go to 140
  120 if (xm.le.0.0e0) go to 130
      b=b+tol1
      go to 140
  130 b=b-tol1
  140 fb=f(b)
      if ((fb*(fc/ abs(fc))).gt.0.0e0) go to 20
      go to 30
  150 seroin=b
      return
      end
