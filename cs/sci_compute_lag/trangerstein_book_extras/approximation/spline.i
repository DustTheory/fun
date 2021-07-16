      integer nspline
c     parameter (nspline=20)
      parameter (nspline=5)
      double precision
     &  potential_min_spline,potential_max_spline,dpotential_spline,
     &  potential_spline(0:nspline),
     &  ak1_spline(0:nspline),ak1_derivative_spline(0:nspline),
     &  bk1_spline(0:nspline),bk1_derivative_spline(0:nspline)
      common/spline_common/
     &  potential_min_spline,potential_max_spline,dpotential_spline,
     &  potential_spline,
     &  ak1_spline,ak1_derivative_spline,
     &  bk1_spline,bk1_derivative_spline
