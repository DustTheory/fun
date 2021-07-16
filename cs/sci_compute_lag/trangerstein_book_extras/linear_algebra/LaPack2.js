function LaPack2() {
}
//************************************************************************
LaPack2.dgbbrd = function( vect, m, n, ncc, kl, ku, AB, ldab,
d, e, q, ldq, PT, ldpt, C, ldc, work, info) {
  throw new Error("not programmed: band matrix");
}
LaPack2.zgbbrd = function( vect, m, n, ncc, kl, ku, AB, ldab,
d, e, q, ldq, PT, ldpt, C, ldc, work, info) {
  throw new Error("not programmed: complex band matrix");
}
//************************************************************************
LaPack2.dgbcon = function( norm, n, kl, ku, AB, ldab, ipiv,
anorm, rcond, work, iwork, info ) {
  throw new Error("not programmed: band matrix");
}
LaPack2.zgbcon = function( norm, n, kl, ku, AB, ldab, ipiv,
anorm, rcond, work, iwork, info ) {
  throw new Error("not programmed: complex band matrix");
}
//************************************************************************
LaPack2.dgebd2 = function( m, n, A, lda, d, e, tauq, taup,
work, info, ioffa, ioffd, ioffe, iofftauq, iofftaup, ioffwork) {
  info.setValue( 0 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -4 );
  if ( info.getValue() < 0 ) {
    Blas2.xerbla( 'dgebd2', -info.getValue() );
    return;
  }
  var alpha = new NumberReference();
  var tau = new NumberReference();
  if ( m >= n ) {
    for ( var i = 1; i <= n; i ++ ) {
      alpha.setValue( A[ ioffa + i - 1 + ( i - 1 ) * lda ] );
      tau.setValue( tauq[ iofftauq + i - 1 ] );
      LaPack1.dlarfg( m - i + 1, alpha, A, 1, tau,
        ioffa + Math.min( i + 1, m ) - 1 + ( i - 1 ) * lda );
      A[ ioffa + i - 1 + ( i - 1 ) * lda ] = alpha.getValue();
      tauq[ iofftauq + i - 1 ] = tau.getValue();
      d[ ioffd + i - 1 ] = A[ ioffa + i - 1 + ( i - 1 ) * lda ];
      A[ ioffa + i - 1 + ( i - 1 ) * lda ] = 1.;
      if ( i < n ) {
        LaPack1.dlarf( 'Left', m - i + 1, n - i, A, 1,
          tauq[ iofftauq + i - 1 ], A, lda, work,
          ioffa + i - 1 + ( i - 1 ) * lda,
          ioffa + i - 1 + i * lda, ioffwork );
      }
      A[ ioffa + i - 1 + ( i - 1 ) * lda ] = d[ ioffd + i - 1 ];
      if ( i < n ) {
        alpha.setValue( A[ ioffa + i - 1 + i * lda ] );
        tau.setValue( taup[ iofftaup + i - 1 ] );
        LaPack1.dlarfg( n - i, alpha, A, lda, tau,
          ioffa + i - 1 + ( Math.min( i + 2, n ) - 1 ) * lda );
        A[ ioffa + i - 1 + i * lda ] = alpha.getValue();
        taup[ iofftaup + i - 1 ] = tau.getValue();
        e[ ioffe + i - 1 ] = A[ ioffa + i - 1 + i * lda ];
        A[ ioffa + i - 1 + i * lda ] = 1.;
        LaPack1.dlarf( 'Right', m - i, n - i, A, lda,
          taup[ iofftaup + i - 1 ], A, lda, work,
          ioffa + i - 1 + i * lda, ioffa + i + i * lda, ioffwork );
        A[ ioffa + i - 1 + i * lda ] = e[ ioffe + i - 1 ];
      } else taup[ iofftaup + i - 1 ] = 0.;
    }
  } else {
    for ( i = 1; i <= m; i ++ ) {
      alpha.setValue( A[ ioffa + i - 1 + ( i - 1 ) * lda ] );
      tau.setValue( taup[ iofftaup + i - 1 ] );
      LaPack1.dlarfg( n - i + 1, alpha, A, lda, tau,
        ioffa + i - 1 + ( Math.min( i + 1, n ) - 1 ) * lda );
      A[ ioffa + i - 1 + ( i - 1 ) * lda ] = alpha.getValue();
      taup[ iofftaup + i - 1 ] = tau.getValue();
      d[ ioffd + i - 1 ] = A[ ioffa + i - 1 + ( i - 1 ) * lda ];
      A[ ioffa + i - 1 + ( i - 1 ) * lda ] = 1.;
      if ( i < m ) {
        LaPack1.dlarf( 'Right', m - i, n - i + 1, A, lda,
          taup[ iofftaup + i - 1 ], A, lda, work,
          ioffa + i - 1 + ( i - 1 ) * lda, ioffa + i + ( i - 1 ) * lda,
          ioffwork );
      }
      A[ ioffa + i - 1 + ( i - 1 ) * lda ] = d[ ioffd + i - 1 ];
      if ( i < m ) {
        alpha.setValue( A[ ioffa + i + ( i - 1 ) * lda ] );
        tau.setValue( tauq[ iofftauq + i - 1 ] );
        LaPack1.dlarfg( m - i, alpha, A, 1, tau,
          ioffa + Math.min( i + 2, m ) - 1 + ( i - 1 ) * lda );
        A[ ioffa + i + ( i - 1 ) * lda ] = alpha.getValue();
        tauq[ iofftauq + i - 1 ] = tau.getValue();
        e[ ioffe + i - 1 ] = A[ ioffa + i + ( i - 1 ) * lda ];
        A[ ioffa + i + ( i - 1 ) * lda ] = 1.;
        LaPack1.dlarf( 'Left', m - i, n - i, A, 1,
          tauq[ iofftauq + i - 1 ], A, lda, work,
          ioffa + i + ( i - 1 ) * lda, ioffa + i + i * lda, ioffwork );
        A[ ioffa + i + ( i - 1 ) * lda ] = e[ ioffe + i - 1 ];
      } else tauq[ iofftauq + i - 1 ] = 0.;
    }
  }
}
LaPack2.zgebd2 = function( m, n, A, lda, d, e, tauq, taup,
work, info) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dgecon = function( norm, n, A, lda, anorm, rcond, work,
iwork, info, ioffa, ioffwork, ioffiwork ) {
  var isave = new Array( 3 );
  info.setValue( 0 );
  var onenrm = ( norm.charAt(0).toUpperCase() == '1' ) ||
    ( norm.charAt(0).toUpperCase() == 'O' );
  if ( ! onenrm && norm.charAt(0).toUpperCase() != 'I' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -4 );
  else if ( anorm < 0. ) info.setValue( -5 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgecon', -info.getValue() );
    return;
  }
  rcond.setValue( 0. );
  if ( n == 0 ) {
    rcond.setValue( 1. );
    return;
  } else if ( anorm == 0. ) return;
  var smlnum = LaPack0.dlamch( 'Safe minimum' );
  var ainvnm = new NumberReference();
  var normin = 'N';
  var kase1 = ( onenrm ? 1 : 2 );
  var kase = new IntReference( 0 );
  var sl = new NumberReference();
  var su = new NumberReference();
  while ( true ) { // 10
    LaPack0.dlacn2( n, work, work, iwork, ainvnm, kase, isave,
      ioffwork + n, ioffwork, ioffiwork, 0 );
    if ( kase.getValue() != 0 ) {
      if ( kase.getValue() == kase1 ) {
        LaPack1.dlatrs( 'Lower', 'No transpose', 'Unit', normin, n, A,
          lda, work, sl, work, info, ioffa, ioffwork, ioffwork + 2*n );
        LaPack1.dlatrs( 'Upper', 'No transpose', 'Non-unit', normin, n,
          A, lda, work, su, work, info, ioffa, ioffwork,
          ioffwork + 3*n );
      } else {
        LaPack1.dlatrs( 'Upper', 'Transpose', 'Non-unit', normin, n,
          A, lda, work, su, work, info, ioffa, ioffwork,
          ioffwork + 3*n );
        LaPack1.dlatrs( 'Lower', 'Transpose', 'Unit', normin, n, A,
          lda, work, sl, work, info, ioffa, ioffwork, ioffwork + 2*n );
      }
      var scale = sl.getValue() * su.getValue();
      normin = 'Y';
      if ( scale != 1. ) {
        var ix = Blas1.idamax( n, work, 1, ioffwork );
        if ( scale != 1. ) {
          ix = Blas1.idamax( n, work, 1, ioffwork );
          if ( scale < Math.abs( work[ ioffwork + ix - 1 ] ) * smlnum
          || scale == 0. ) {
            return;
          }
          LaPack1.drscl( n, scale, work, 1, ioffwork );
        }
      }
    } else break;
  }
  if ( ainvnm.getValue() != 0. ) {
    rcond.setValue( ( 1. / ainvnm.getValue() ) / anorm );
  }
}
LaPack2.zgecon = function( norm, n, A, lda, anorm, rcond, work,
iwork, info ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dgehd2 = function( n, ilo, ihi, A, lda, tau, work,
info, ioffa, iofftau, ioffwork ) {
  info.setValue( 0 );
  if ( n < 0 ) info.setValue( -1 );
  else if ( ilo < 1 || ilo > Math.max( 1, n ) ) info.setValue( -2 );
  else if ( ihi < Math.min( ilo, n ) || ihi > n ) info.setValue( -3 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -5 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgehd2', - info.getValue() );
    return;
  }
  for ( var i = ilo; i <= ihi - 1; i ++ ) {
    var alpha =
      new NumberReference( A[ ioffa + i + ( i - 1 ) * lda ] );
    var tauref =
      new NumberReference( tau[ iofftau + i - 1 ] );
    LaPack1.dlarfg( ihi - i, alpha, A, 1, tauref,
      ioffa + Math.min( i + 2, n ) - 1 + ( i - 1 ) * lda );
    A[ ioffa + i + ( i - 1 ) * lda ] = alpha.getValue();
    tau[ iofftau + i - 1 ] = tauref.getValue();
    var aii = A[ ioffa + i + ( i - 1 ) * lda ];
    A[ ioffa + i + ( i - 1 ) * lda ] = 1.;
    LaPack1.dlarf( 'Right', ihi, ihi - i, A, 1, tau[ iofftau + i - 1 ],
      A, lda, work, ioffa + i + ( i - 1 ) * lda,
      ioffa + i * lda, ioffwork );
    LaPack1.dlarf( 'Left', ihi - i, n - i, A, 1,
      tau[ iofftau + i - 1 ], A, lda, work,
      ioffa + i + ( i - 1 ) * lda, ioffa + i + i * lda, ioffwork );
    A[ ioffa + i + ( i - 1 ) * lda ] = aii;
  }
}
LaPack2.zgehd2 = function( n, ilo, ihi, A, lda, tau, work,
info ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dgelq2 = function( m, n, A, lda, tau, work, info,
ioffa, iofftau, ioffwork ) {
  info.setValue( 0 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -4 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgelq2', - info.getValue() );
    return;
  }
  var k = Math.min( m, n );
  for ( var i = 1; i <= k; i ++ ) {
    var alpha =
      new NumberReference( A[ ioffa + i - 1 + ( i - 1 ) * lda ] );
    var tauref =
      new NumberReference( tau[ iofftau + i - 1 ] );
    LaPack1.dlarfg( n - i + 1, alpha, A, lda, tauref,
      ioffa + i - 1 + ( Math.min( i + 1, n ) - 1 ) * lda );
    A[ ioffa + i - 1 + ( i - 1 ) * lda ] = alpha.getValue();
    tau[ iofftau + i - 1 ] = tauref.getValue();
    if ( i < m ) {
      var aii = A[ ioffa + i - 1 + ( i - 1 ) * lda ];
      A[ ioffa + i - 1 + ( i - 1 ) * lda ] = 1.;
      LaPack1.dlarf( 'Right', m - i, n - i + 1, A, lda,
        tau[ iofftau + i - 1 ], A, lda, work,
        ioffa + i - 1 + ( i - 1 ) * lda,
        ioffa + i + ( i - 1 ) * lda, ioffwork );
      A[ ioffa + i - 1 + ( i - 1 ) * lda ] = aii;
    }
  }
}
LaPack2.zgelq2 = function( m, n, A, lda, tau, work, info ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dgeql2 = function( m, n, A, lda, tau, work, info,
ioffa, iofftau, ioffwork ) {
  info.setValue( 0 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -4 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgeql2', - info.getValue() );
    return;
  }
  var k = Math.min( m, n );
  for ( var i = k; i >= 1; i -- ) {
    var alpha = new NumberReference(
      A[ ioffa + m - k + i - 1 + ( n - k + i - 1 ) * lda ] );
    var tauref =
      new NumberReference( tau[ iofftau + i - 1 ] );
    LaPack1.dlarfg( m - k + i, alpha, A, 1, tauref,
      ioffa + ( n - k + i - 1 ) * lda );
    A[ ioffa + m - k + i - 1 + ( n - k + i - 1 ) * lda ] = alpha.getValue();
    tau[ iofftau + i - 1 ] = tauref.getValue();
    var aii =
      A[ ioffa + m - k + i - 1 + ( n - k + i - 1 ) * lda ];
    A[ ioffa + m - k + i - 1 + ( n - k + i - 1 ) * lda ] = 1.;
    LaPack1.dlarf( 'Left', m - k + i, n - k + i - 1, A, 1,
      tau[ iofftau + i - 1 ], A, lda, work,
      ioffa + ( n - k + i - 1 ) * lda, ioffa, ioffwork );
    A[ ioffa + m - k + i - 1 + ( n - k + i - 1 ) * lda ] = aii;
  }
}
LaPack2.zgeql2 = function( m, n, A, lda, tau, work, info ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dgeqr2 = function( m, n, A, lda, tau, work, info,
ioffa, iofftau, ioffwork ) {
  info.setValue( 0 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -4 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgeqr2', - info.getValue() );
    return;
  }
  var k = Math.min( m, n );
  for ( var i = 1; i <= k; i ++ ) {
    var alpha =
      new NumberReference( A[ ioffa + i - 1 + ( i - 1 ) * lda ] );
    var tauref =
      new NumberReference( tau[ iofftau + i - 1 ] );
    LaPack1.dlarfg( m - i + 1, alpha, A, 1, tauref,
      ioffa + Math.min( i + 1 , m ) - 1 + ( i - 1 ) * lda );
    A[ ioffa + i - 1 + ( i - 1 ) * lda ] = alpha.getValue();
    tau[ iofftau + i - 1 ] = tauref.getValue();
    if ( i < n ) {
      var aii = A[ ioffa + i - 1 + ( i - 1 ) * lda ];
      A[ ioffa + i - 1 + ( i - 1 ) * lda ] = 1.;
      LaPack1.dlarf( 'Left', m - i + 1, n - i, A, 1,
        tau[ iofftau + i - 1 ], A, lda, work,
        ioffa + i - 1 + ( i - 1 ) * lda, ioffa + i - 1 + i * lda,
        ioffwork );
      A[ ioffa + i - 1 + ( i - 1 ) * lda ] = aii;
    }
  }
}
LaPack2.zgeqr2 = function( m, n, A, lda, tau, work, info ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dgeqr2p = function( m, n, A, lda, tau, work, info,
ioffa, iofftau, ioffwork ) {
  throw new Error("not tested: Fortran version missing");
  info.setValue( 0 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -4 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgeqr2p', - info.getValue() );
    return;
  }
  var k = Math.min( m, n );
  for ( var i = 1; i <= k; i ++ ) {
    var alpha =
      new NumberReference( A[ ioffa + i - 1 + ( i - 1 ) * lda ] );
    var tauref =
      new NumberReference( tau[ iofftau + i - 1 ] );
    LaPack1.dlarfgp( m - i + 1, alpha, A, 1, tauref,
      ioffa + Math.min( i + 1, m ) - 1 + ( i - 1 ) * lda );
    A[ ioffa + i - 1 + ( i - 1 ) * lda ] = alpha.getValue();
    tau[ iofftau + i - 1 ] = tauref.getValue();
    if ( i < n ) {
      var aii = A[ ioffa + i - 1 + ( i - 1 ) * lda ];
      A[ ioffa + i - 1 + ( i - 1 ) * lda ] = 1.;
      LaPack1.dlarf( 'Left', m - i + 1, n - i, A, 1,
        tau[ iofftau + i - 1 ], A, lda, work,
        ioffa + i - 1 + ( i - 1 ) * lda, ioffa + i - 1 + i * lda,
        ioffwork );
      A[ ioffa + i - 1 + ( i - 1 ) * lda ] = aii;
    }
  }
}
//************************************************************************
LaPack2.dgeqrt2 = function( m, n, A, lda, T, ldt, info, ioffa,
iofft ) {
  throw new Error("not programmed: compact WY storage");
}
//************************************************************************
LaPack2.dgeqrt3 = function( m, n, A, lda, T, ldt, info, ioffa,
iofft ) {
  throw new Error("not programmed: compact WY storage");
}
//************************************************************************
LaPack2.dgerfs = function( trans, n, nrhs, A, lda, AF, ldaf,
ipiv, B, ldb, X, ldx, ferr, berr, work, iwork, info, ioffa, ioffaf,
ioffipiv, ioffb, ioffx, ioffferr, ioffberr, ioffwork, ioffiwork) {
  var itmax = 5;
  var isave = new Array( 3 );
  info.setValue( 0 );
  var notran = ( trans.charAt(0).toUpperCase() == 'N' );
  if ( ! notran && trans.charAt(0).toUpperCase() != 'T'
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( nrhs < 0 ) info.setValue( -3 );
  if ( lda < Math.max( 1, n ) ) info.setValue( -5 );
  if ( ldaf < Math.max( 1, n ) ) info.setValue( -7 );
  if ( ldb < Math.max( 1, n ) ) info.setValue( -10 );
  if ( ldx < Math.max( 1, n ) ) info.setValue( -12 );
  if ( info.getValue() != 0) {
    Blas2.xerbla( 'dgerfs', -info.getValue() );
    return;
  }
  if ( n == 0 || nrhs == 0 ) {
    for ( var j = 1; j <= nrhs; j ++ ) {
      ferr[ ioffferr + j - 1 ] = 0.;
      berr[ ioffberr + j - 1 ] = 0.;
    }
    return;
  }
  var transt = ( notran ? 'T' : 'N' );
  var nz = n + 1;
  var eps = LaPack0.dlamch( 'Epsilon' );
  var safmin = LaPack0.dlamch( 'Safe minimum' );
  var safe1 = nz * safmin;
  var safe2 = safe1 / eps;
  for ( j = 1; j <= nrhs; j ++ ) {
    var count = 1;
    var lstres = 3.;
    while ( true ) { // 20
      Blas1.dcopy( n, B, 1, work, 1, ioffb + ( j - 1 ) * ldb,
        ioffwork + n );
      Blas2.dgemv( trans, n, n, -1., A, lda, X, 1, 1., work, 1, ioffa,
        ioffx + ( j - 1 ) * ldx, ioffwork + n );
      for ( var i = 1; i <= n; i ++ ) {
        work[ ioffwork + i - 1 ] =
          Math.abs( B[ ioffb + i - 1 + ( j - 1 ) * ldb ] );
      }
      if ( notran ) {
        for ( var k = 1; k <= n; k ++ ) {
          var xk =
            Math.abs( X[ ioffx + k - 1 + ( j - 1 ) * ldx ] );
          for( i = 1; i <= n; i ++ ) {
            work[ ioffwork + i - 1 ] +=
              Math.abs( A[ ioffa + i - 1 + ( k - 1 ) * lda ] ) * xk;
          }
        }
      } else {
        for ( k = 1; k <= n; k ++ ) {
          var s = 0.;
          for ( i = 1; i <= n; i ++ ) {
            s += Math.abs( A[ ioffa + i - 1 + ( k  - 1 ) * lda ] )
              * Math.abs( X[ ioffx + i - 1 + ( j - 1 ) * ldx ] );
          }
          work[ ioffwork + k - 1 ] += s;
        }
      }
      s = 0.;
      for ( i = 1; i <= n; i ++ ) {
        if ( work[ ioffwork + i - 1 ] > safe2 ) {
          s = Math.max( s, Math.abs( work[ ioffwork + n + i - 1 ] )
            / work[ ioffwork + i - 1 ] );
        } else {
          s = Math.max( s,
            ( Math.abs( work[ ioffwork + n + i - 1 ] ) + safe1 )
            / ( work[ ioffwork + i - 1 ] + safe1 ) );
        }
      }
      berr[ ioffberr + j - 1 ] = s;
      if ( berr[ ioffberr + j - 1 ] > eps &&
      2. * berr[ ioffberr + j - 1 ] <= lstres
      && count <= itmax ) {
        LaPack1.dgetrs( trans, n, 1, AF, ldaf, ipiv, work, n,
          info, ioffaf, ioffipiv, ioffwork + n );
        Blas1.daxpy( n, 1., work, 1, X, 1, ioffwork + n,
          ioffx + ( j - 1 ) * ldx );
        lstres = berr[ ioffberr + j - 1 ];
        count ++;
      } else break;
    }
    for ( i = 1; i <= n; i ++ ) {
      if ( work[ ioffwork + i - 1 ] > safe2 ) {
        work[ ioffwork + i - 1 ] =
          Math.abs( work[ ioffwork + n + i - 1 ] )
          + nz * eps * work[ ioffwork + i - 1 ] ;
      } else {
        work[ ioffwork + i - 1 ] =
          Math.abs( work[ ioffwork + n + i - 1 ] )
          + nz * eps * work[ ioffwork + i - 1 ] + safe1;
      }
    }
    var kase = new IntReference( 0 );
    while ( true ) { // 100
      var est =
        new NumberReference( ferr[ ioffferr + j - 1 ] );
      LaPack0.dlacn2( n, work, work, iwork, est,
        kase, isave, ioffwork + 2 * n, ioffwork + n, ioffiwork, 0 );
      ferr[ ioffferr + j - 1 ] = est.getValue();
      if ( kase.getValue() != 0 ) {
        if ( kase.getValue() == 1 ) {
          LaPack1.dgetrs( transt, n, 1, AF, ldaf, ipiv, work, n, info,
            ioffaf, ioffipiv, ioffwork + n );
          for ( i = 1; i <= n; i ++ ) {
            work[ ioffwork + n + i - 1 ] *= work[ ioffwork + i - 1 ];
          }
        } else {
          for ( i = 1; i <= n; i ++ ) {
            work[ ioffwork + n + i - 1 ] *= work[ ioffwork + i - 1 ];
          }
          LaPack1.dgetrs( transt, n, 1, AF, ldaf, ipiv, work, n, info,
            ioffaf, ioffipiv, ioffwork + n );
        }
      } else break;
    }
    lstres = 0.;
    for ( i = 1; i <= n; i ++ ) {
      lstres = Math.max( lstres,
        Math.abs( X[ ioffx + i - 1 + ( j - 1 ) * ldx ] ) );
    }
    if ( lstres != 0. ) ferr[ ioffferr + j - 1 ] /= lstres;
  }
}
LaPack2.zgerfs = function( trans, n, nrhs, A, lda, AF, ldaf,
ipiv, B, ldb, X, ldx, ferr, berr, work, iwork, info ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dgerq2 = function( m, n, A, lda, tau, work, info,
ioffa, iofftau, ioffwork ) {
  info.setValue( 0 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -4 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgerq2', - info.getValue() );
    return;
  }
  var k = Math.min( m, n );
  for ( var i = k; i >= 1; i -- ) {
    var alpha = new NumberReference(
      A[ ioffa + m - k + i - 1 + ( n - k + i - 1 ) * lda ] );
    var tauref =
      new NumberReference( tau[ iofftau + i - 1 ] );
    LaPack1.dlarfg( n - k + i, alpha, A, lda, tauref,
      ioffa + m - k + i - 1 );
    A[ ioffa + m - k + i - 1 + ( n - k + i - 1 ) * lda ] = alpha.getValue();
    tau[ iofftau + i - 1 ] = tauref.getValue();
    var aii =
      A[ ioffa + m - k + i - 1 + ( n - k + i - 1 ) * lda ];
    A[ ioffa + m - k + i - 1 + ( n - k + i - 1 ) * lda ] = 1.;
    LaPack1.dlarf( 'Right', m - k + i - 1, n - k + i, A, lda,
      tau[ iofftau + i - 1 ], A, lda, work, ioffa + m - k + i - 1,
      ioffa, ioffwork );
    A[ ioffa + m - k + i - 1 + ( n - k + i - 1 ) * lda ] = aii;
  }
}
LaPack2.zgerq2 = function( m, n, A, lda, tau, work, info ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dgetrf = function( m, n, A, lda, ipiv, info, ioffa,
ioffipiv ) {
  info.setValue( 0 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -4 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgetrf', -info.getValue() );
    return;
  }
  if ( m == 0 || n == 0 ) return;
  var nb = LaPack0.ilaenv( 1, 'dgetrf', ' ', m, n, -1, -1 );
  if ( nb <= 1 || nb >= Math.min( m, n ) ) {
    LaPack1.dgetf2( m, n, A, lda, ipiv, info, ioffa, ioffipiv );
  } else {
    for ( var j = 1; j <= Math.min( m, n ); j += nb ) {
      var jb = Math.min( Math.min( m , n ) - j + 1, nb );
      var iinfo = new NumberReference( 0 );
      LaPack1.dgetf2( m - j + 1, jb, A, lda, ipiv, iinfo,
        ioffa + j - 1 + ( j - 1 ) * lda, ioffipiv + j - 1 ); 
      if ( info.getValue() == 0 && iinfo.getValue() > 0 ) {
        info.setValue( iinfo.getValue() + j - 1 );
      }
      for ( var i = j; i <= Math.min( m, j + jb - 1 ); i ++ ) {
        ipiv[ ioffipiv + i - 1 ] += j - 1;
      }
      LaPack0.dlaswp( j - 1, A, lda, j, j + jb - 1, ipiv, 1, ioffa,
        ioffipiv );
      if ( j + jb <= n ) {
        LaPack0.dlaswp( n - j - jb + 1, A, lda, j, j + jb - 1, ipiv, 1,
          ioffa + ( j + jb - 1 ) * lda, ioffipiv );
        Blas3.dtrsm( 'Left', 'Lower', 'No transpose', 'Unit', jb,
          n - j - jb + 1, 1., A, lda, A, lda,
          ioffa + j - 1 + ( j - 1 ) * lda,
          ioffa + j - 1 + ( j + jb - 1 ) * lda );
        if ( j + jb <= m ) {
          Blas3.dgemm( 'No transpose', 'No transpose', m - j - jb + 1,
            n - j - jb + 1, jb, -1., A, lda, A, lda, 1., A, lda,
            ioffa + j + jb - 1 + ( j - 1 ) * lda,
            ioffa + j - 1 + ( j + jb - 1 ) * lda,
            ioffa + j + jb - 1 + ( j + jb - 1 ) * lda );
        }
      }
    }
  }
}
LaPack2.zgetrf = function( m, n, A, lda, ipiv, info ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dgghrd = function( compq, compz, n, ilo, ihi, A, lda,
B, ldb, Q, ldq, Z, ldz, info, ioffa, ioffb, ioffq, ioffz ) {
  throw new Error("not tested: generalized eigenvalue problem");
  if ( compq.charAt(0).toUpperCase() == 'N' ) {
    var ilq = false;
    var icompq = 1;
  } else if ( compq.charAt(0).toUpperCase() == 'V' ) {
    ilq = true;
    icompq = 2;
  } else if ( compq.charAt(0).toUpperCase() == 'I' ) {
    ilq = true;
    icompq = 3;
  } else icompq = 0;
  if ( compz.charAt(0).toUpperCase() == 'N' ) {
    var ilz = false;
    var icompz = 1;
  } else if ( compz.charAt(0).toUpperCase() == 'V' ) {
    ilz = true;
    icompz = 2;
  } else if ( compz.charAt(0).toUpperCase() == 'I' ) {
    ilz = true;
    icompz = 3;
  } else icompz = 0;
  info.setValue( 0 );
  if ( icompq <= 0 ) info.setValue( -1 );
  else if ( icompz <= 0 ) info.setValue( -2 );
  else if ( n < 0 ) info.setValue( -3 );
  else if ( ilo < 1 ) info.setValue( -4 );
  else if ( ihi > n || ihi < ilo - 1 ) info.setValue( -5 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -7 );
  else if ( ldb < Math.max( 1, n ) ) info.setValue( -9 );
  else if ( ( ilq && ldq < n ) || ldq < 1 ) info.setValue( -11 );
  else if ( ( ilz && ldz < n ) || ldz < 1 ) info.setValue( -13 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgghrd', - info.getValue() );
    return;
  }
  if ( icompq == 3 ) {
    LaPack0.dlaset( 'Full', n, n, 0., 1., Q, ldq, ioffq );
  }
  if ( icompz == 3 ) {
    LaPack0.dlaset( 'Full', n, n, 0., 1., Z, ldz, ioffz );
  }
  if ( n <= 1 ) return;
  for ( var jcol = 1; jcol <= n - 1; jcol ++ ) {
    for ( var jrow = jcol + 1; jrow <= n; jrow ++ ) {
      B[ ioffb + jrow - 1 + ( jcol - 1 ) * ldb ] = 0.;
    }
  }
  for ( jcol = ilo; jcol <= ihi - 2; jcol ++ ) {
    for ( jrow = ihi; jrow >= jcol + 2; jrow -- ) {
      var temp = A[ ioffa + jrow - 2 + ( jcol - 1 ) * lda ];
      var c = new NumberReference();
      var s = new NumberReference();
      var r = new NumberReference(
        A[ ioffa + jrow - 2 + ( jcol - 1 ) * lda ] );
      LaPack1.dlartg( temp, A[ ioffa + jrow - 1 + ( jcol - 1 ) * lda ],
        c, s, r );
      A[ ioffa + jrow - 2 + ( jcol - 1 ) * lda ] = r.getValue();
      A[ ioffa + jrow - 1 + ( jcol - 1 ) * lda ] = 0.;
      Blas1.drot( n - jcol, A, lda, A, lda, c.getValue(), s.getValue(),
        ioffa + jrow - 2 + jcol * lda, ioffa + jrow - 1 + jcol * lda );
      Blas1.drot( n + 2 - jrow, B, ldb, B, ldb, c.getValue(), s.getValue(),
        ioffb + jrow - 2 + ( jrow - 2 ) * ldb,
        ioffb + jrow - 1 + ( jrow - 2 ) * ldb );
      if ( ilq ) {
        Blas1.drot( n, Q, 1, Q, 1, c.getValue(), s.getValue(),
          ioffq + ( jrow - 2 ) * ldq, ioffq + ( jrow - 1 ) * ldq );
      }
      temp = B[ ioffb + jrow - 1 + ( jrow - 1 ) * ldb ];
      r.setValue( B[ ioffb + jrow - 1 + ( jrow - 1 ) * ldb ] );
      LaPack1.dlartg( temp, B[ ioffb + jrow - 1 + ( jrow - 2 ) * ldb ],
        c, s, r );
      B[ ioffb + jrow - 1 + ( jrow - 1 ) * ldb ] = r.getValue();
      B[ ioffb + jrow - 1 + ( jrow - 2 ) * ldb ] = 0.;
      Blas1.drot( ihi, A, 1, A, 1, c.getValue(), s.getValue(),
        ioffa + ( jrow - 1 ) * lda, ioffa + ( jrow - 2 ) * lda );
      Blas1.drot( jrow - 1, B, 1, B, 1, c.getValue(), s.getValue(),
        ioffb + ( jrow - 1 ) * ldb, ioffb + ( jrow - 2 ) * ldb );
      if ( ilz ) {
        Blas1.drot( n, Z, 1, Z, 1, c.getValue(), s.getValue(),
          ioffz + ( jrow - 1 ) * ldz, ioffz + ( jrow - 2 ) * ldz );
      }
    }
  }
}
LaPack2.zgghrd = function( compq, compz, n, ilo, ihi, A, lda,
B, ldb, Q, ldq, Z, ldz, info, ioffa, ioffb, ioffq, ioffz ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dgsvj0 = function( jobv, m, n, A, lda, d, sva, mv, V,
ldv, eps, sfmin, tol, nsweep, work, lwork, info, ioffa, ioffd, ioffsva,
ioffv, ioffwork ) {
  throw new Error("not tested: complex input");
  var fastr = new Array( 5 );
  var applv = ( jobv.charAt(0).toUpperCase() == 'A' );
  var rsvec = ( jobv.charAt(0).toUpperCase() == 'V' );
  if ( ! ( rsvec || applv || jobv.charAt(0).toUpperCase() == 'N') ) {
    info.setValue( -1 );
  } else if ( m < 0 ) info.setValue( -2 );
  else if ( n < 0 || n > m ) info.setValue( -3 );
  else if ( lda < m ) info.setValue( -5 );
  else if ( ( rsvec || applv ) && mv < 0 ) info.setValue( -8 );
  else if ( ( rsvec && ldv < n ) || ( applv && ldv < mv ) ) {
    info.setValue( -10 );
  } else if ( tol <= eps ) info.setValue( -13 );
  else if ( nsweep < 0 ) info.setValue( -14 );
  else if ( lwork < m ) info.setValue( -16 );
  else info.setValue( 0 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgsvj0', - info.getValue() );
    return;
  }
  if ( rsvec ) var mvl = n;
  else if ( applv ) mvl = mv;
  rsvec = rsvec || applv;
  var rooteps = Math.sqrt( eps );
  var rootsfmin = Math.sqrt( sfmin );
  var small = sfmin / eps;
  var big = 1. / sfmin;
  var rootbig = 1. / rootsfmin;
  var bigtheta = 1. / rooteps;
  var roottol = Math.sqrt( tol );
  var emptsw = ( n * ( n - 1 ) ) / 2;
  var notrot = 0;
  fastr[ 0 ] = 0.;
  var swband = 0;
  var kbl = Math.min( 8, n );
  var nbl = n / kbl;
  if ( nbl * kbl != n ) nbl ++;
  var blskip = ( kbl * kbl ) + 1;
  var rowskip = Math.min( 5, kbl );
  var lkahead = 1;
//    swband = 0;
  var pskipped = 0;
  var goto1994 = false;
  for ( var i = 1; i <= nsweep; i ++ ) { // ==> 1993
    var mxaapq = 0.;
    var mxsinj = 0.;
    var iswrot = 0;
    notrot = 0;
    pskipped = 0;
    for ( var ibr = 1; ibr <= nbl; ibr ++ ) { // ==> 2000
      var igl = ( ibr - 1 ) * kbl + 1;
      for ( var ir1 = 0; ir1 <= Math.min( lkahead, nbl - ibr );
      ir1 ++ ) { // ==> 1002
        igl += ir1 * kbl;
        for ( var p = igl; p <= Math.min( igl + kbl - 1, n - 1 );
        p ++ ) { // ==> 2001
          var q = Blas1.idamax( n - p + 1, sva, 1,
            ioffsva + p - 1 ) + p - 1;
          if ( p != q ) {
            Blas1.dswap( m, A, 1, A, 1, ioffa + ( p - 1 ) * lda,
              ioffa + ( q - 1 ) * lda );
            if ( rsvec ) {
              Blas1.dswap( mvl, V, 1, V, 1, ioffv + ( p - 1 ) * ldv,
                ioffv + ( q - 1 ) * ldv );
            }
            var temp1 = sva[ ioffsva + p - 1 ];
            sva[ ioffsva + p - 1 ] = sva[ ioffsva + q - 1 ];
            sva[ ioffsva + q - 1 ] = temp1;
            temp1 = d[ ioffd + p - 1 ];
            d[ ioffd + p - 1 ] = d[ ioffd + q - 1 ];
            d[ ioffd + q - 1 ] = temp1;
          }
          if ( ir1 == 0 ) {
            if ( sva[ ioffsva + p - 1 ] == rootbig &&
            sva[ ioffsva + p - 1 ] > rootsfmin ) {
              sva[ ioffsva + p - 1 ] =
                Blas1.dnrm2( m, A, 1, ioffa + ( p - 1 ) * lda )
                * d[ ioffd + p - 1 ];
            } else {
              var temp1ref = new NumberReference( 0. );
              var aappref = new NumberReference( 1. );
              LaPack0.dlassq( m, A, 1, temp1ref, aappref,
                ioffa + ( p - 1 ) * lda );
              temp1 = temp1ref.getValue();
              var aapp = aappref.getValue();
              sva[ ioffsva + p - 1 ] = temp1 * Math.sqrt( aapp )
                * d[ ioffd + p - 1 ];
            }
            aapp = sva[ ioffsva + p - 1 ];
          } else aapp = sva[ ioffsva + p - 1 ];
          if ( aapp > 0. ) {
            pskipped = 0;
            for ( q = p + 1; q <= Math.min( igl + kbl - 1, n ); q ++ )
            { // 2002
              var aaqq = sva[ ioffsva + q - 1 ];
              if ( aaqq > 0. ) {
                var aapp0 = aapp;
                if ( aaqq >= 1. ) {
                  var rotok = ( small * aapp <= aaqq );
                  if ( aapp < big / aaqq ) {
                    var aapq =
                      ( Blas1.ddot( m, A, 1, A, 1,
                        ioffa + ( p - 1 ) * lda,
                        ioffa + ( q - 1 ) * lda )
                      * d[ ioffd + p - 1 ] * d[ ioffd + q - 1 ]
                      / aaqq ) / aapp;
                  } else {
                    Blas1.dcopy( m, A, 1, work, 1,
                      ioffa + ( p - 1 ) * lda, ioffwork );
                    var ierr = new IntReference();
                    LaPack1.dlascl( 'G', 0, 0, aapp,
                      d[ ioffd + p - 1 ], m, 1, work, lda, ierr,
                      ioffwork );
                    aapq = Blas1.ddot( m, work, 1, A, 1, ioffwork,
                      ioffa + ( q - 1 ) * lda ) * d[ ioffd + q - 1 ]
                      / aaqq;
                  }
                } else {
                  rotok = ( aapp <= aaqq / small );
                  if ( aapp > small / aaqq ) {
                    aapq = ( Blas1.ddot( m, A, 1, A, 1,
                        ioffa + ( p - 1 ) * lda,
                        ioffa + ( q - 1 ) * lda )
                      * d[ ioffd + p - 1 ] * d[ ioffd + q - 1 ]
                      / aaqq ) / aapp;
                  } else {
                    Blas1.dcopy( m, A, 1, work, 1,
                      ioffa + ( q - 1 ) * lda, ioffwork );
                    LaPack1.dlascl( 'G', 0, 0, aaqq,
                      d[ ioffd + q - 1 ], m, 1, work, lda, ierr, 
                      ioffwork );
                    aapq = Blas1.ddot( m, work, 1, A, 1, ioffwork,
                      ioffa + ( p - 1 ) * lda ) * d[ ioffd + p - 1 ]
                      / aapp;
                  }
                }
                mxaapq = Math.max( mxaapq, Math.abs( aapq ) );
                if ( Math.abs( aapq ) > tol ) {
                  if ( ir1 == 0 ) {
                    notrot = 0;
                    pskipped = 0;
                    iswrot ++;
                  }
                  if ( rotok ) {
                    var aqoap = aaqq / aapp;
                    var apoaq = aapp / aaqq;
                    var theta =
                      -0.5 * Math.abs( aqoap - apoaq ) / aapq;
                    if ( Math.abs( theta ) > bigtheta ) {
                      var t = 0.5 / theta;
                      fastr[ 2 ] =
                        t * d[ ioffd + p - 1 ] / d[ ioffd + q - 1 ];
                      fastr[ 3 ] =
                        - t * d[ ioffd + q - 1 ] / d[ ioffd + p - 1 ];
                      Blas1.drotm( m, A, 1, A, 1, fastr,
                        ioffa + ( p - 1 ) * lda,
                        ioffa + ( q - 1 ) * lda, 0 );
                      if ( rsvec ) {
                        Blas1.drotm( mvl, V, 1, V, 1, fastr,
                          ioffv + ( p - 1 ) * ldv,
                          ioffv + ( q - 1 ) * ldv, 0 );
                      }
                      sva[ ioffsva + q - 1 ] = aaqq * Math.sqrt(
                        Math.max( 0., 1. + t * apoaq * aapq ) );
                      aapp *= Math.sqrt(
                        Math.max( 0., 1. - t * aqoap * aapq ) );
                      mxsinj = Math.max( mxsinj, Math.abs( t ) );
                    } else {
                      var thsign = ( aapq >= 0. ? -1. : 1. );
                      t = 1. / ( theta
                        + thsign * Math.sqrt( 1. + theta * theta ) );
                      var cs = Math.sqrt( 1. / ( 1. + t * t ) );
                      var sn = t * cs;
                      mxsinj = Math.max( mxsinj, Math.abs( sn ) );
                      sva[ ioffsva + q - 1 ] = aaqq * Math.sqrt(
                        Math.max( 0., 1. + t * apoaq * aapq ) );
                      aapp *= Math.sqrt( Math.max( 0.,
                        1. - t * aqoap * aapq ) );
                      apoaq = d[ ioffd + p - 1 ] / d[ ioffd + q - 1 ];
                      aqoap = d[ ioffd + q - 1 ] / d[ ioffd + p - 1 ];
                      if ( d[ ioffd + p - 1 ] >= 1. ) {
                        if ( d[ ioffd + q - 1 ] >= 1. ) {
                          fastr[ 2 ] = t * apoaq;
                          fastr[ 3 ] = - t * aqoap;
                          d[ ioffd + p - 1 ] *= cs;
                          d[ ioffd + q - 1 ] *= cs;
                          Blas1.drotm( m, A, 1, A, 1, fastr,
                            ioffa + ( p - 1 ) * lda,
                            ioffa + ( q - 1 ) * lda, 0 );
                          if ( rsvec ) {
                            Blas1.drotm( mvl, V, 1, V, 1, fastr,
                             ioffv + ( p - 1 ) * ldv,
                             ioffv + ( q - 1 ) * ldv, 0 );
                          }
                        } else {
                          Blas1.daxpy( m, - t * aqoap, A, 1, A, 1,
                            ioffa + ( q - 1 ) * lda,
                            ioffa + ( p - 1 ) * lda );
                          Blas1.daxpy( m, cs * sn * apoaq, A, 1, A, 1,
                            ioffa + ( p - 1 ) * lda,
                            ioffa + ( q - 1 ) * lda );
                          d[ ioffd + p - 1 ] *= cs;
                          d[ ioffd + q - 1 ] /= cs;
                          if ( rsvec ) {
                            Blas1.daxpy( mvl, -t * aqoap, V, 1, V, 1,
                              ioffv + ( q - 1 ) * ldv,
                              ioffv + ( p - 1 ) * ldv );
                            Blas1.daxpy( mvl, cs * sn * apoaq, V, 1,
                              V, 1, ioffv + ( p - 1 ) * ldv,
                              ioffv + ( q - 1 ) * ldv );
                          }
                        }
                      } else {
                        if ( d[ ioffd + q - 1 ] >= 1. ) {
                          Blas1.daxpy( m, t * apoaq, A, 1, A, 1,
                            ioffa + ( p - 1 ) * lda,
                            ioffa + ( q - 1 ) * lda );
                          Blas1.daxpy( m, - cs * sn * aqoap, A, 1,
                            A, 1, ioffa + ( q - 1 ) * lda,
                            ioffa + ( p - 1 ) * lda );
                          d[ ioffd + p - 1 ] /= cs;
                          d[ ioffd + q - 1 ] *= cs;
                          if ( rsvec ) {
                            Blas1.daxpy( mvl, t * apoaq, V, 1, V, 1,
                              ioffv + ( p - 1 ) * ldv,
                              ioffv + ( q - 1 ) * ldv );
                            Blas1.daxpy( mvl, - cs * sn * aqoap, V, 1,
                              V, 1, ioffv + ( q - 1 ) * ldv,
                              ioffv + ( p - 1 ) * ldv );
                          }
                        } else {
                          if ( d[ ioffd + p - 1 ] >=
                          d[ ioffd + q - 1 ] ) {
                            Blas1.daxpy( m, - t * aqoap, A, 1, A, 1,
                              ioffa + ( q - 1 ) * lda,
                              ioffa + ( p - 1 ) * lda );
                            Blas1.daxpy( m, cs * sn * apoaq, A, 1,
                              A, 1, ioffa + ( p - 1 ) * lda,
                              ioffa + ( q - 1 ) * lda );
                            d[ ioffd + p - 1 ] *= cs;
                            d[ ioffd + q - 1 ] /= cs;
                            if ( rsvec ) {
                              Blas1.daxpy( mvl, - t * aqoap, V, 1,
                                V, 1, ioffv + ( q - 1 ) * ldv,
                                ioffv + ( p - 1 ) * ldv );
                              Blas1.daxpy( mvl, cs * sn * apoaq, V, 1,
                                V, 1, ioffv + ( p - 1 ) * ldv,
                                ioffv + ( q - 1 ) * ldv );
                            }
                          } else {
                            Blas1.daxpy( m, t * apoaq, A, 1, A, 1,
                              ioffa + ( p - 1 ) * lda,
                              ioffa + ( q - 1 ) * lda );
                            Blas1.daxpy( m, - cs * sn * aqoap, A, 1,
                              A, 1, ioffa + ( q - 1 ) * lda,
                              ioffa + ( p - 1 ) * lda );
                            d[ ioffd + p - 1 ] /= cs;
                            d[ ioffd + q - 1 ] *= cs;
                            if ( rsvec ) {
                              Blas1.daxpy( mvl, t * apoaq, V, 1,
                                V, 1, ioffv + ( p - 1 ) * ldv,
                                ioffv + ( q - 1 ) * ldv );
                              Blas1.daxpy( mvl, - cs * sn * aqoap,
                                V, 1, V, 1, ioffv + ( q - 1 ) * ldv,
                                ioffv + ( p - 1 ) * ldv );
                            }
                          }
                        }
                      }
                    }
                  } else {
                    Blas1.dcopy( m, A, 1, work, 1,
                      ioffa + ( p - 1 ) * lda, ioffwork );
                    LaPack1.dlascl( 'G', 0, 0, aapp, 1., m, 1, work,
                      lda, ierr, ioffwork );
                    LaPack1.dlascl( 'G', 0, 0, aaqq, 1., m, 1, A,
                      lda, ierr, ioffa + ( q - 1 ) * lda );
                    temp1 = - aapq * d[ ioffd + p - 1 ]
                      / d[ ioffd + q - 1 ];
                    Blas1.daxpy( m, temp1, work, 1, A, 1, ioffwork,
                      ioffa + ( q - 1 ) * lda );
                    LaPack1.dlascl( 'G', 0, 0, 1., aaqq, m, 1, A, lda,
                      ierr, ioffa + ( q - 1 ) * lda );
                    sva[ ioffsva + q - 1 ] = aaqq
                      * Math.sqrt( Math.max( 0., 1. - aapq * aapq ) );
                    mxsinj = Math.max( mxsinj, sfmin );
                  } // end if rotok 
                  if ( Math.pow( sva[ ioffsva + q - 1 ] / aaqq, 2 )
                  <= rooteps ) {
                    if ( aaqq < rootbig && aaqq > rootsfmin ) {
                      sva[ ioffsva + q - 1 ] =
                        Blas1.dnrm2( m, A, 1, ioffa + ( q - 1 ) * lda )
                        * d[ ioffd + q - 1 ];
                    } else {
                      tref.setValue( 0. );
                      var aaqqref =
                        new NumberReference( 1. );
                      LaPack0.dlassq( m, A, 1, tref, aaqqref,
                        ioffa + ( q - 1 ) * lda );
                      t = tref.getValue();
                      aaqq = aaqqref.getValue();
                      sva[ ioffsva + q - 1 ] = t * Math.sqrt( aaqq )
                        * d[ ioffd + q - 1 ];
                    }
                  }
                  if ( aapp / aapp0 <= rooteps ) {
                    if ( aapp < rootbig && aapp > rootsfmin ) {
                      aapp =
                        Blas1.dnrm2( m, A, 1, ioffa + ( p - 1 ) * lda )
                        * d[ ioffd + p - 1 ];
                    } else {
                      tref.setValue( 0. );
                      aappref.setValue( 1. );
                      LaPack0.dlassq( m, A, 1, tref, aappref,
                        ioffa + ( p - 1 ) * lda );
                      t = tref.getValue();
                      aapp = aappref.getValue();
                      aapp = t * Math.sqrt( aapp )
                        * d[ ioffd + p - 1 ];
                    }
                    sva[ ioffsva + p - 1 ] = aapp;
                  }
                } else {
                  if ( ir1 == 0 ) notrot ++;
                  pskipped ++;
                }
              } else {
                if ( ir1 == 0 ) notrot ++;
                pskipped ++;
              }
              if ( i <= swband  && pskipped > rowskip ) {
                if ( ir1 == 0 ) aapp = - aapp;
                notrot = 0;
                break; // goto 2103
              }
            } // 2002
            sva[ ioffsva + p - 1 ] = aapp; // 2103
          } else {
            sva[ ioffsva + p - 1 ] = aapp;
            if ( ir1 == 0 && aapp == 0. ) {
              notrot += Math.min( igl + kbl - 1, n ) - p;
            }
          }
        }  // 2001
      } // 1002
      igl = ( ibr - 1 ) * kbl + 1;
      var goto2011 = false;
      for ( var jbc = ibr + 1; jbc <= nbl; jbc ++ ) { // ==> 2010
        var jgl = ( jbc - 1 ) * kbl + 1;
        var ijblsk = 0;
        for ( p = igl; p <= Math.min( igl + kbl - 1, n ); p ++ )
        { // ==> 2100
          aapp = sva[ ioffsva + p - 1 ];
          if ( aapp > 0. ) {
            pskipped = 0;
            for ( q = jgl; q <= Math.min( jgl + kbl - 1, n ); q ++ )
            { // ==> 2200
              aaqq = sva[ ioffsva + q - 1 ];
              if ( aaqq > 0. ) {
                aapp0 = aapp;
                if ( aaqq >= 1. ) {
                  rotok = ( aapp >= aaqq ? small * aapp <= aaqq :
                    small * aaqq <= aapp );
                  if ( aapp < big / aaqq ) {
                    aapq = ( Blas1.ddot( m, A, 1, A, 1,
                      ioffa + ( p - 1 ) * lda,
                      ioffa + ( q - 1 ) * lda ) * d[ ioffd + p - 1 ]
                      * d[ ioffd + q - 1 ] / aaqq ) / aapp;
                  } else {
                    Blas1.dcopy( m, A, 1, work, 1,
                      ioffa + ( p - 1 ) * lda, ioffwork );
                    LaPack1.dlascl( 'G', 0, 0, aapp, d[ ioffd + p - 1],
                      m, 1, work, lda, ierr, ioffwork );
                    aapq = Blas1.ddot( m, work, 1, A, 1, ioffwork,
                      ioffa + ( q - 1 ) * lda ) * d[ ioffd + q - 1 ]
                      / aaqq;
                  }
                } else {
                  rotok = ( aapp >= aaqq ? aapp <= aaqq / small :
                    aaqq <= aapp / small );
                  if ( aapp > small / aaqq ) {
                    aapq = ( Blas1.ddot( m, A, 1, A, 1,
                      ioffa + ( p - 1 ) * lda,
                      ioffa + ( q - 1 ) * lda ) * d[ ioffd + p - 1 ]
                      * d[ ioffd + q - 1 ] / aaqq ) / aapp;
                  } else {
                    Blas1.dcopy( m, A, 1, work, 1,
                      ioffa + ( q - 1 ) * lda, ioffwork );
                    LaPack1.dlascl( 'G', 0, 0, aaqq, d[ ioffd + q - 1],
                      m, 1, work, lda, ierr, ioffwork );
                    aapq = Blas1.ddot( m, work, 1, A, 1, ioffwork,
                      ioffa + ( p - 1 ) * lda ) * d[ ioffd + p - 1 ]
                      / aapp;
                  }
                }
                mxaapq = Math.max( mxaapq, Math.abs( aapq ) );
                if ( Math.abs( aapq ) > tol ) {
                  notrot = 0;
                  pskipped = 0;
                  iswrot ++;
                  if ( rotok ) {
                    aqoap = aaqq / aapp;
                    apoaq = aapp / aaqq;
                    theta = - 0.5 * Math.abs( aqoap - apoaq ) / aapq;
                    if ( aaqq > aapp0 ) theta = - theta;
                    if ( Math.abs( theta ) > bigtheta ) {
                      t = 0.5 / theta;
                      fastr[ 2 ] = t * d[ ioffd + p - 1 ]
                        / d[ ioffd + q - 1 ];
                      fastr[ 3 ] = - t * d[ ioffd + q - 1 ]
                        / d[ ioffd + p - 1 ];
                      Blas1.drotm( m, A, 1, A, 1, fastr,
                        ioffa + ( p - 1 ) * lda,
                        ioffa + ( q - 1 ) * lda, 0 );
                      if ( rsvec ) {
                        Blas1.drotm( mvl, V, 1, V, 1, fastr,
                          ioffv + ( p - 1 ) * ldv,
                          ioffv + ( q - 1 ) * ldv, 0 );
                      }
                      sva[ ioffsva + q - 1 ] = aaqq * Math.sqrt(
                        Math.max( 0., 1. + t * apoaq * aapq ) );
                      aapp *= Math.sqrt( Math.max( 0.,
                        1. - t * aqoap * aapq ) );
                      mxsinj = Math.max( mxsinj, Math.abs( t ) );
                    } else {
                      thsign = ( aapq >= 0. ? -1. : 1. );
                      if ( aaqq > aapp0 ) thsign = - thsign;
                      t = 1. / ( theta
                        + thsign * Math.sqrt( 1. + theta * theta ) );
                      cs = Math.sqrt( 1. / ( 1. + t *  t ) );
                      sn = t * cs;
                      mxsinj = Math.max( mxsinj, Math.abs( sn ) );
                      sva[ ioffsva + q - 1 ] = aaqq * Math.sqrt(
                        Math.max( 0., 1. + t * apoaq * aapq ) );
                      aapp *= Math.sqrt( Math.max( 0.,
                        1. - t * aqoap * aapq ) );
                      apoaq = d[ ioffd + p - 1 ] / d[ ioffd + q - 1 ];
                      aqoap = d[ ioffd + q - 1 ] / d[ ioffd + p - 1 ];
                      if ( d[ ioffd + p - 1 ] >= 1. ) {
                        if ( d[ ioffd + q - 1 ] >= 1. ) {
                          fastr[ 2 ] = t * apoaq;
                          fastr[ 3 ] = -t * aqoap;
                          d[ ioffd + p - 1 ] *= cs;
                          d[ ioffd + q - 1 ] *= cs;
                          Blas1.drotm( m, A, 1, A, 1, fastr,
                            ioffa + ( p - 1 ) * lda,
                            ioffa + ( q - 1 ) * lda, 0 );
                          if ( rsvec ) {
                            Blas1.drotm( mvl, V, 1, V, 1, fastr,
                              ioffv + ( p - 1 ) * ldv,
                              ioffv + ( q - 1 ) * ldv, 0 );
                          }
                        } else {
                          Blas1.daxpy( m, - t * aqoap, A, 1, A, 1,
                            ioffa + ( q - 1 ) * lda,
                            ioffa + ( p - 1 ) * lda );
                          Blas1.daxpy( m, cs * sn * apoaq, A, 1, A, 1,
                            ioffa + ( p - 1 ) * lda,
                            ioffa + ( q - 1 ) * lda );
                          if ( rsvec ) {
                            Blas1.daxpy( mvl, - t * aqoap, V, 1, V, 1,
                              ioffv + ( q - 1 ) * ldv,
                              ioffv + ( p - 1 ) * ldv );
                            Blas1.daxpy( mvl, cs * sn * apoaq, V, 1,
                              V, 1, ioffv + ( p - 1 ) * ldv,
                              ioffv + ( q - 1 ) * ldv );
                          }
                          d[ ioffd + p - 1 ] *= cs;
                          d[ ioffd + q - 1 ] /= cs;
                        }
                      } else {
                        if ( d[ ioffd + q - 1 ] >= 1. ) {
                          Blas1.daxpy( m, t * apoaq, A, 1, A, 1,
                            ioffa + ( p - 1 ) * lda,
                            ioffa + ( q - 1 ) * lda );
                          Blas1.daxpy( m, -cs * sn * aqoap, A, 1, A, 1,
                            ioffa + ( q - 1 ) * lda,
                            ioffa + ( p - 1 ) * lda );
                          if ( rsvec ) {
                            Blas1.daxpy( mvl, t * apoaq, V, 1, V, 1,
                              ioffv + ( p - 1 ) * ldv,
                              ioffv + ( q - 1 ) * ldv );
                            Blas1.daxpy( mvl, -cs * sn * aqoap, V, 1,
                              V, 1, ioffv + ( q - 1 ) * ldv,
                              ioffv + ( p - 1 ) * ldv );
                          }
                          d[ ioffd + p - 1 ] /= cs;
                          d[ ioffd + q - 1 ] *= cs;
                        } else {
                          if ( d[ ioffd + p - 1 ] >=
                          d[ ioffd + q - 1 ] ) {
                            Blas1.daxpy( m, - t * aqoap, A, 1, A, 1,
                              ioffa + ( q - 1 ) * lda,
                              ioffa + ( p - 1 ) * lda );
                            Blas1.daxpy( m, cs * sn * aqoap, A, 1,
                              A, 1, ioffa + ( p - 1 ) * lda,
                              ioffa + ( q - 1 ) * lda );
                            d[ ioffd + p - 1 ] *= cs;
                            d[ ioffd + q - 1 ] /= cs;
                            if ( rsvec ) {
                              Blas1.daxpy( mvl, -t * apoaq, V, 1, V, 1,
                                ioffv + ( q - 1 ) * ldv,
                                ioffv + ( p - 1 ) * ldv );
                              Blas1.daxpy( mvl, cs * sn * aqoap, V, 1,
                                V, 1, ioffv + ( p - 1 ) * ldv,
                                ioffv + ( q - 1 ) * ldv );
                            }
                          } else {
                            Blas1.daxpy( m, t * apoaq, A, 1, A, 1,
                              ioffa + ( p - 1 ) * lda,
                              ioffa + ( q - 1 ) * lda );
                            Blas1.daxpy( m, - cs * sn * aqoap, A, 1,
                              A, 1, ioffa + ( q - 1 ) * lda,
                              ioffa + ( p - 1 ) * lda );
                            d[ ioffd + p - 1 ] /= cs;
                            d[ ioffd + q - 1 ] *= cs;
                            if ( rsvec ) {
                              Blas1.daxpy( mvl, t * apoaq, V, 1, V, 1,
                                ioffv + ( p - 1 ) * ldv,
                                ioffv + ( q - 1 ) * ldv );
                              Blas1.daxpy( mvl, -cs * sn * aqoap, V, 1,
                                V, 1, ioffv + ( q - 1 ) * ldv,
                                ioffv + ( p - 1 ) * ldv );
                            }
                          }
                        }
                      }
                    }
                  } else {
                    if ( aapp > aaqq ) {
                      Blas1.dcopy( m, A, 1, work, 1,
                        ioffa + ( p - 1 ) * lda, ioffwork );
                      LaPack1.dlascl( 'G', 0, 0, aapp, 1., m, 1, work,
                        lda, ierr, ioffwork );
                      LaPack1.dlascl( 'G', 0, 0, aaqq, 1., m, 1, A,
                        lda, ierr, ioffa + ( q - 1 ) * lda );
                      temp1 = - aapq * d[ ioffd + p - 1 ]
                        / d[ ioffd + q - 1 ];
                      Blas1.daxpy( m, temp1, work, 1, A, 1,
                        ioffwork, ioffa + ( q - 1 ) * lda );
                      LaPack1.dlascl( 'G', 0, 0, 1., aaqq, m, 1,
                        A, lda, ierr, ioffa + ( q - 1 ) * lda );
                      sva[ ioffsva + q - 1 ] = aaqq * Math.sqrt(
                        Math.max( 0., 1. - aapq * aapq ) );
                      mxsinj = Math.max( mxsinj, sfmin );
                    } else {
                      Blas1.dcopy( m, A, 1, work, 1,
                        ioffa + ( q - 1 ) * lda, ioffwork );
                      LaPack1.dlascl( 'G', 0, 0, aaqq, 1., m, 1, work,
                        lda, ierr, ioffwork );
                      LaPack1.dlascl( 'G', 0, 0, aapp, 1., m, 1, A,
                        lda, ierr, ioffa + ( p - 1 ) * lda );
                      temp1 = - aapq * d[ ioffd + q - 1 ]
                        / d[ ioffd + p - 1 ];
                      Blas1.daxpy( m, temp1, work, 1, A, 1,
                        ioffwork, ioffa + ( p - 1 ) * lda );
                      LaPack1.dlascl( 'G', 0, 0, 1., aapp, m, 1,
                        A, lda, ierr, ioffa + ( p - 1 ) * lda );
                      sva[ ioffsva + p - 1 ] = aapp * Math.sqrt(
                        Math.max( 0., 1. - aapq * aapq ) );
                      mxsinj = Math.max( mxsinj, sfmin );
                    }
                  } // rotok
                  if ( Math.pow( sva[ ioffsva + q - 1 ] / aaqq, 2 )
                  <= rooteps ) {
                    if ( aaqq < rootbig && aaqq > rootsfmin ) {
                      sva[ ioffsva + q - 1 ] =
                        Blas1.dnrm2( m, A, 1, ioffa + ( q - 1 ) * lda )
                        * d[ ioffd + q - 1 ];
                    } else {
                      tref.setValue( 0. );
                      aaqqref.setValue( 1. );
                      LaPack0.dlassq( m, A, 1, tref, aaqqref,
                        ioffa + ( q - 1 ) * lda );
                      t = tref.getValue();
                      aaqq = aaqqref.getValue();
                      sva[ ioffsva + q - 1 ] = t
                        * Math.sqrt( aaqq ) * d[ ioffd + q - 1 ];
                    }
                  }
                  if ( Math.pow( aapp / aapp0, 2 ) <= rooteps ) {
                    if ( aapp < rootbig && aapp > rootsfmin ) {
                      aapp =
                        Blas1.dnrm2( m, A, 1, ioffa + ( p - 1 ) * lda )
                        * d[ ioffd + p - 1 ];
                    } else {
                      tref.setValue( 0. );
                      aappref.setValue( 1. );
                      LaPack0.dlassq( m, A, 1, tref, aappref,
                        ioffa + ( p - 1 ) * lda );
                      t = tref.getValue();
                      aapp = aappref.getValue();
                      aapp = t * Math.sqrt( aapp )
                        * d[ ioffd + p - 1 ];
                    }
                    sva[ ioffsva + p - 1 ] = aapp;
                  }
                } else {
                  notrot ++;
                  pskipped ++;
                  ijblsk ++;
                } // abs( aapq ) > tol
              } else {
                notrot ++;
                pskipped ++;
                ijblsk ++;
              } // aaqq > 0
              if ( i <= swband && ijblsk >= blskip ) {
                sva[ ioffsva + p - 1 ] = aapp;
                notrot = 0;
                goto2011 = true;
                break;
              }
              if ( i <= swband && pskipped > rowskip ) {
                aapp = - aapp;
                notrot = 0;
                break; // ==> 2203
              }
            } // 2200
            if ( goto2011 ) break;
            sva[ ioffsva + p - 1 ] = aapp; // 2203
          } else {
            if ( aapp == 0. ) {
              notrot += Math.min( jgl + kbl - 1, n ) - jgl + 1;
            }
            if ( aapp < 0. ) notrot = 0;
          } // aapp > 0
        } // 2100
        if ( goto2011 ) break;
      }  // 2010
      for ( p = igl; p <= Math.min( igl + kbl - 1, n ); p ++ ) {
        sva[ ioffsva + p - 1 ] = Math.abs( sva[ ioffsva + p - 1 ] );
      } // 2012
    } // 2000
    if ( sva[ ioffsva + n - 1 ] < rootbig &&
    sva[ ioffsva + n - 1 ] > rootsfmin ) {
      sva[ ioffsva + n - 1 ] =
        Blas1.dnrm2( m, A, 1, ioffa + ( n - 1 ) * lda )
        * d[ ioffd + n - 1 ];
    } else {
      var tref = new NumberReference( 0. );
      aappref.setValue( 1. );
      LaPack0.dlassq( m, A, 1, tref, aappref,
        ioffa + ( n - 1 ) * lda );
      t = tref.getValue();
      aapp = aappref.getValue();
      sva[ ioffsva + n - 1 ] = t * Math.sqrt( aapp )
        * d[ ioffd + n - 1 ];
    }
    if ( i < swband && ( mxaapq <= roottol || iswrot <= n ) ) {
      swband = i;
    }
    if ( i > swband + 1 && mxaapq < Number( n ) * tol &&
    Number( n ) * mxaapq * mxsinj < tol ) {
      goto1994 = true;
      break;
    }
    if ( notrot >= emptsw ) {
      goto1994 = true;
      break;
    }
  } // 1993
  if ( ! goto1994 ) info.setValue( nsweep - 1 );
  else info.setValue( 0 ); // 1994
  for ( p = 1; p <= n - 1; p ++ ) { // 1995
    q = Blas1.idamax( n - p + 1, sva, 1, ioffsva + p - 1 ) + p - 1;
    if ( p != q ) {
      temp1 = sva[ ioffsva + p - 1 ];
      sva[ ioffsva + p - 1 ] = sva[ ioffsva + q - 1 ];
      sva[ ioffsva + q - 1 ] = temp1;
      temp1 = d[ ioffd + p - 1 ];
      d[ ioffd + p - 1 ] = d[ ioffd + q - 1 ];
      d[ ioffd + q - 1 ] = temp1;
      Blas1.dswap( m, A, 1, A, 1, ioffa + ( p - 1 ) * lda,
        ioffa + ( q - 1 ) * lda );
      if ( rsvec ) {
        Blas1.dswap( mvl, V, 1, V, 1, ioffv + ( p - 1 ) * ldv,
          ioffv + ( q - 1 ) * ldv );
      }
    }
  } // 5991
}
//************************************************************************
LaPack2.dgsvj1 = function( jobv, m, n, n1, A, lda, d, sva, mv,
V, ldv, eps, sfmin, tol, nsweep, work, lwork, info, ioffa, ioffd,
ioffsva, ioffv, ioffwork ) {
  throw new Error("not tested: complex input");
  var fastr = new Array( 5 );
  var applv = ( jobv.charAt(0).toUpperCase() == 'A' );
  var rsvec = ( jobv.charAt(0).toUpperCase() == 'V' );
  if ( ! ( rsvec || applv || jobv.charAt(0).toUpperCase() == 'N' ) ) {
    info.setValue( -1 );
  } else if ( m < 0 ) info.setValue( -2 );
  else if ( n < 0 || n > m ) info.setValue( -3 );
  else if ( n1 < 0 ) info.setValue( -4 );
  else if ( lda < m ) info.setValue( -6 );
  else if ( ( rsvec || applv ) && mv < 0 ) info.setValue( -9 );
  else if ( ( rsvec && ldv < n ) || ( applv && ldv < mv ) ) {
    info.setValue( -11 );
  } else if ( tol <= eps ) info.setValue( -14 );
  else if ( nsweep < 0 ) info.setValue( -15 );
  else if ( lwork < m ) info.setValue( -17 );
  else info.setValue( 0 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgsvj1', - info.getValue() );
    return;
  }
  if ( rsvec ) var mvl = n;
  else if ( applv ) mvl = mv;
  rsvec = rsvec || applv;
  var rooteps = Math.sqrt( eps );
  var rootsfmin = Math.sqrt( sfmin );
  var small = sfmin / eps;
  var big = 1. / sfmin;
  var rootbig = 1. / rootsfmin;
  var large = big / Math.sqrt( Number( m * n ) );
  var bigtheta = 1. / rooteps;
  var roottol = Math.sqrt( tol );
  var emptsw = n1 * ( n - n1 );
  var notrot = 0;
  fastr[ 0 ] = 0.;
  var kbl = Math.min( 8, n );
  var nblr = n1 / kbl;
  if ( nblr * kbl != n1 ) nblr ++;
  var nblc = ( n - n1 ) / kbl;
  if ( nblc * kbl != n - n1 ) nblc ++;
  var blskip = kbl * kbl + 1;
  var rowskip = Math.min( 5, kbl );
  var swband = 0;
  var goto1994 = false;
  for ( var i = 1; i <= nsweep; i ++ ) { // ==> 1993
    var mxaapq = 0.;
    var mxsinj = 0.;
    var iswrot = 0;
    notrot = 0;
    var pskipped = 0;
    for ( var ibr = 1; ibr <= nblr; ibr ++ ) { // ==> 2000
      var igl = ( ibr - 1 ) * kbl + 1;
      var goto2011 = false;
      for ( var jbc = 1; jbc <= nblc; jbc ++ ) { // ==> 2010
        var jgl = n1 + ( jbc - 1 ) * kbl + 1;
        var ijblsk = 0;
        for ( var p = igl; p <= Math.min( igl + kbl - 1, n1 );
        p ++ ) { // ==> 2100
          var aapp = sva[ ioffsva + p - 1 ];
          if ( aapp > 0. ) {
            pskipped = 0;
            for ( var q = jgl; q <= Math.min( jgl + kbl - 1, n );
            q ++ ) { // ==> 2200
              var aaqq = sva[ ioffsva + q - 1 ];
              if ( aaqq > 0. ) {
                var aapp0= aapp;
                if ( aaqq >= 1. ) {
                  var rotok = ( aapp >= aaqq ?
                    small * aapp <= aaqq : small * aaqq <= aapp );
                  if ( aapp < big / aaqq ) {
                    var aapq =
                      ( Blas1.ddot( m, A, 1, A, 1,
                        ioffa + ( p - 1 ) * lda,
                        ioffa + ( q - 1 ) * lda )
                      * d[ ioffd + p - 1 ] * d[ ioffd + q - 1 ]
                      / aaqq ) / aapp;
                  } else {
                    Blas1.dcopy( m, A, 1, work, 1,
                      ioffa + ( p - 1 ) * lda, ioffwork );
                    var ierr = new IntReference();
                    LaPack1.dlascl( 'G', 0, 0, aapp,
                      d[ ioffd + p - 1 ], m, 1, work, lda, ierr,
                      ioffwork );
                    aapq = Blas1.ddot( m, work, 1, A, 1, ioffwork,
                      ioffa + ( q - 1 ) * lda ) * d[ ioffd + q - 1 ]
                      / aaqq;
                  }
                } else {
                  rotok = ( aapp >= aaqq ? aapp <= aaqq / small :
                    aaqq <= aapp / small );
                  if ( aapp > small / aaqq ) {
                    aapq = ( Blas1.ddot( m, A, 1, A, 1,
                      ioffa + ( p - 1 ) * lda,
                      ioffa + ( q - 1 ) * lda )
                      * d[ ioffd + p - 1 ] * d[ ioffd + q - 1 ]
                      / aaqq ) / aapp;
                  } else {
                    Blas1.dcopy( m, A, 1, work, 1,
                      ioffa + ( q - 1 ) * lda, ioffwork );
                    LaPack1.dlascl( 'G', 0, 0, aaqq,
                      d[ ioffd + q - 1 ], m, 1, work, lda, ierr,
                      ioffwork );
                    aapq = Blas1.ddot( m, work, 1, A, 1, ioffwork,
                      ioffa + ( p - 1 ) * lda ) * d[ ioffd + p - 1 ]
                      / aapp;
                  }
                } // aaqq >= 1
                mxaapq = Math.max( mxaapq, Math.abs( aapq ) );
                if ( Math.abs( aapq ) > tol ) {
                  notrot = 0;
                  pskipped = 0;
                  iswrot ++;
                  if ( rotok ) {
                    var aqoap = aaqq / aapp; 
                    var apoaq = aapp / aaqq; 
                    var theta =
                      - 0.5 * Math.abs( aqoap - apoaq ) / aapq;
                    if ( aaqq > aapp0 ) theta = - theta;
                    if ( Math.abs( theta ) > bigtheta ) {
                      var t = 0.5 / theta;
                      fastr[ 2 ] = t * d[ ioffd + p - 1 ]
                        / d[ ioffd + q - 1 ];
                      fastr[ 3 ] = - t * d[ ioffd + q - 1 ]
                        / d[ ioffd + p - 1 ];
                      Blas1.drotm( m, A, 1, A, 1, fastr,
                        ioffa + ( p - 1 ) * lda,
                        ioffa + ( q - 1 ) * lda, 0 );
                      if ( rsvec ) {
                        Blas1.drotm( mvl, V, 1, V, 1, fastr,
                          ioffv + ( p - 1 ) * ldv,
                          ioffv + ( q - 1 ) * ldv, 0 );
                      }
                      sva[ ioffsva + q - 1 ] = aaqq * Math.sqrt(
                        Math.max( 0., 1. + t * apoaq * aapq ) );
                      aapp *= Math.sqrt( Math.max( 0.,
                        1. - t * aqoap * aapq ) );
                      mxsinj = Math.max( mxsinj, Math.abs( t ) );
                    } else {
                      var thsign = ( aapq >= 0. ? -1. : 1. );
                      if ( aaqq > aapp0 ) thsign = - thsign;
                      t = 1. / ( theta
                        + thsign * Math.sqrt( 1. + theta * theta ) );
                      var cs = Math.sqrt( 1. / ( 1. + t * t ) );
                      var sn = t * cs;
                      mxsinj = Math.max( mxsinj, Math.abs( sn ) );
                      sva[ ioffsva + q - 1 ] = aaqq * Math.sqrt(
                        Math.max( 0., 1. + t * apoaq * aapq ) );
                      aapp *= Math.sqrt( Math.max( 0.,
                        1. - t * aqoap * aapq ) );
                      apoaq = d[ ioffd + p - 1 ] / d[ ioffd + q - 1 ];
                      aqoap = d[ ioffd + q - 1 ] / d[ ioffd + p - 1 ];
                      if ( d[ ioffd + p - 1 ] >= 1. ) {
                        if ( d[ ioffd + q - 1 ] >= 1. ) {
                          fastr[ 2 ] = t * apoaq;
                          fastr[ 3 ] = - t * aqoap;
                          d[ ioffd + p - 1 ] *= cs;
                          d[ ioffd + q - 1 ] *= cs;
                          Blas1.drotm( m, A, 1, A, 1, fastr,
                            ioffa + ( p - 1 ) * lda,
                            ioffa + ( q - 1 ) * lda, 0 );
                          if ( rsvec ) {
                            Blas1.drotm( mvl, V, 1, V, 1, fastr,
                              ioffv + ( p - 1 ) * ldv,
                              ioffv + ( q - 1 ) * ldv, 0 );
                          }
                        } else {
                          Blas1.daxpy( m, - t * aqoap, A, 1, A, 1,
                            ioffa + ( q - 1 ) * lda,
                            ioffa + ( p - 1 ) * lda );
                          Blas1.daxpy( m, cs * sn * apoaq, A, 1, A, 1,
                            ioffa + ( p - 1 ) * lda,
                            ioffa + ( q - 1 ) * lda );
                          if ( rsvec ) {
                            Blas1.daxpy( mvl, -t * aqoap, V, 1, V, 1,
                              ioffv + ( q - 1 ) * ldv,
                              ioffv + ( p - 1 ) * ldv );
                            Blas1.daxpy( mvl, cs * sn * apoaq, V, 1,
                              V, 1, ioffv + ( p - 1 ) * ldv,
                              ioffv + ( q - 1 ) * ldv );
                          }
                          d[ ioffd + p - 1 ] *= cs;
                          d[ ioffd + q - 1 ] /= cs;
                        }
                      } else {
                        if ( d[ ioffd + q - 1 ] >= 1. ) {
                          Blas1.daxpy( m, t * apoaq, A, 1, A, 1,
                            ioffa + ( p - 1 ) * lda,
                            ioffa + ( q - 1 ) * lda );
                          Blas1.daxpy( m, -cs * sn * aqoap, A, 1, A, 1,
                            ioffa + ( q - 1 ) * lda,
                            ioffa + ( p - 1 ) * lda );
                          if ( rsvec ) {
                            Blas1.daxpy( mvl, t * apoaq, V, 1, V, 1,
                              ioffv + ( p - 1 ) * ldv,
                              ioffv + ( q - 1 ) * ldv );
                            Blas1.daxpy( mvl, -cs * sn * aqoap, V, 1,
                              V, 1, ioffv + ( q - 1 ) * ldv,
                              ioffv + ( p - 1 ) * ldv );
                          }
                          d[ ioffd + p - 1 ] /= cs;
                          d[ ioffd + q - 1 ] *= cs;
                        } else {
                          if ( d[ ioffd + p - 1 ] >=
                          d[ ioffd + q - 1 ] ) {
                            Blas1.daxpy( m, - t * aqoap, A, 1, A, 1,
                              ioffa + ( q - 1 ) * lda,
                              ioffa + ( p - 1 ) * lda );
                            Blas1.daxpy( m, cs * sn * apoaq, A, 1,
                              A, 1, ioffa + ( p - 1 ) * lda,
                              ioffa + ( q - 1 ) * lda );
                            d[ ioffd + p - 1 ] *= cs;
                            d[ ioffd + q - 1 ] /= cs;
                            if ( rsvec ) {
                              Blas1.daxpy( mvl, -t * aqoap, V, 1, V, 1,
                                ioffv + ( q - 1 ) * ldv,
                                ioffv + ( p - 1 ) * ldv );
                              Blas1.daxpy( mvl, cs * sn * apoaq, V, 1,
                                V, 1, ioffv + ( p - 1 ) * ldv,
                                ioffv + ( q - 1 ) * ldv );
                            }
                          } else {
                            Blas1.daxpy( m, t * apoaq, A, 1, A, 1,
                              ioffa + ( p - 1 ) * lda,
                              ioffa + ( q - 1 ) * lda );
                            Blas1.daxpy( m, -cs * sn * aqoap, A, 1,
                              A, 1, ioffa + ( q - 1 ) * lda,
                              ioffa + ( p - 1 ) * lda );
                            d[ ioffd + p - 1 ] /= cs;
                            d[ ioffd + q - 1 ] *= cs;
                            if ( rsvec ) {
                              Blas1.daxpy( mvl, t * apoaq, V, 1, V, 1,
                                ioffv + ( p - 1 ) * ldv,
                                ioffv + ( q - 1 ) * ldv );
                              Blas1.daxpy( mvl, -cs * sn * aqoap, V, 1,
                                V, 1, ioffv + ( q - 1 ) * ldv,
                                ioffv + ( p - 1 ) * ldv );
                            } // rsvec
                          } // d_p >= d_q
                        } // d_q >= 1
                      } // d[ p ] >= 1
                    } // | theta | > bigtheta
                  } else { // rotok
                    if ( aapp > aaqq ) {
                      Blas1.dcopy( m, A, 1, work, 1,
                        ioffa + ( p - 1 ) * lda, ioffwork );
                      LaPack1.dlascl( 'G', 0, 0, aapp, 1., m, 1, work,
                        lda, ierr, ioffwork );
                      LaPack1.dlascl( 'G', 0, 0, aaqq, 1., m, 1, A,
                        lda, ierr, ioffa + ( q - 1 ) * lda );
                      var temp1 = - aapq * d[ ioffd + p - 1 ]
                        / d[ ioffd + q - 1 ];
                      Blas1.daxpy( m, temp1, work, 1, A, 1, ioffwork,
                        ioffa + ( q - 1 ) * lda );
                      LaPack1.dlascl( 'G', 0, 0, 1., aaqq, m, 1,
                        A, lda, ierr, ioffa + ( q - 1 ) * lda );
                      sva[ ioffsva + q - 1 ] = aaqq * Math.sqrt(
                        Math.max( 0., 1. - aapq * aapq ) );
                      mxsinj = Math.max( mxsinj, sfmin );
                    } else {
                      Blas1.dcopy( m, A, 1, work, 1,
                        ioffa + ( q - 1 ) * lda, ioffwork );
                      LaPack1.dlascl( 'G', 0, 0, aaqq, 1., m, 1, work,
                        lda, ierr, ioffwork );
                      LaPack1.dlascl( 'G', 0, 0, aapp, 1., m, 1, A,
                        lda, ierr, ioffa + ( p - 1 ) * lda );
                      temp1 = - aapq * d[ ioffd + q - 1 ]
                        / d[ ioffd + p - 1 ];
                      Blas1.daxpy( m, temp1, work, 1, A, 1, ioffwork,
                        ioffa + ( p - 1 ) * lda );
                      LaPack1.dlascl( 'G', 0, 0, 1., aapp, m, 1,
                        A, lda, ierr, ioffa + ( p - 1 ) * lda );
                      sva[ ioffsva + p - 1 ] = aapp * Math.sqrt(
                        Math.max( 0., 1. - aapq * aapq ) );
                      mxsinj = Math.max( mxsinj, sfmin );
                    }
                  } // rotok
                  if ( Math.pow( sva[ ioffsva + q - 1 ] / aaqq, 2 )
                  <= rooteps ) {
                    if ( aaqq < rootbig && aaqq > rootsfmin ) {
                      sva[ ioffsva + q - 1 ] =
                        Blas1.dnrm2( m, A, 1, ioffa + ( q - 1 ) * lda )
                        * d[ ioffd + q - 1 ];
                    } else {
                      var tref =
                        new NumberReference( 0. );
                      var aaqqref = 
                        new NumberReference( 1. );
                      LaPack0.dlassq( m, A, 1, tref, aaqqref,
                        ioffa + ( q - 1 ) * lda );
                      t = tref.getValue();
                      aaqq = aaqqref.getValue();
                      sva[ ioffsva + q - 1 ] = t * Math.sqrt( aaqq )
                        * d[ ioffd + q - 1 ];
                    }
                  }
                  if ( Math.pow( aapp / aapp0, 2 ) <= rooteps ) {
                    if ( aapp < rootbig && aapp > rootsfmin ) {
                      aapp =
                        Blas1.dnrm2( m, A, 1, ioffa + ( p - 1 ) * lda )
                        * d[ ioffd + p - 1 ];
                    } else {
                      tref.setValue( 0. );
                      var aappref =
                        new NumberReference( 1. );
                      LaPack0.dlassq( m, A, 1, tref, aappref,
                        ioffa + ( p - 1 ) * lda );
                      t = tref.getValue();
                      aapp = aappref.getValue()
                      aapp = t * Math.sqrt( aapp )
                        * d[ ioffd + p - 1 ];
                    }
                    sva[ ioffsva + p - 1 ] = aapp;
                  }
                } else {
                  notrot ++;
                  pskipped ++;
                  ijblsk ++;
                } // | aapq | > tol
              } else {
                notrot ++;
                pskipped ++;
                ijblsk ++;
              } // aaqq > 0
              if ( i <= swband && ijblsk >= blskip ) {
                sva[ ioffsva + p - 1 ] = aapp;
                notrot = 0;
                goto2011 = true;
                break;
              }
              if ( i <= swband && pskipped > rowskip ) {
                aapp = - aapp;
                notrot = 0;
                break;
              }
            } // 2200
            if ( goto2011 ) break; // 2203
            sva[ ioffsva + p - 1 ] = aapp;
          } else {
            if ( aapp == 0. ) {
              notrot += Math.min( jgl + kbl - 1, n ) - jgl + 1;
            }
            if ( aapp < 0. ) notrot = 0;
          } // aapp > 0
        } // 2100
        if ( goto2011 ) break;
      } // 2010
      for ( p = igl; p <= Math.min( igl + kbl - 1, n ); p ++ ) {
        sva[ ioffsva + p - 1 ] = Math.abs( sva[ ioffsva + p - 1 ] );
      }
    } // 2000
    if ( sva[ ioffsva + n - 1 ] < rootbig &&
    sva[ ioffsva + n - 1 ] > rootsfmin ) {
      sva[ ioffsva + n - 1 ] =
        Blas1.dnrm2( m, A, 1, ioffa + ( n - 1 ) * lda )
        * d[ ioffd + n - 1 ];
    } else {
      tref.setValue( 0. );
      aappref.setValue( 1. );
      LaPack0.dlassq( m, A, 1, tref, aappref,
        ioffa + ( n - 1 ) * lda );
      t = tref.getValue();
      aapp = aappref.getValue();
      sva[ ioffsva + n - 1 ] =
        t * Math.sqrt( aapp ) * d[ ioffd + n - 1 ];
    }
    if ( i < swband && ( mxaapq <= roottol || iswrot <= n ) ) {
      swband = i;
    }
    if ( i > swband + 1 && mxaapq < Number( n ) * tol &&
    Number( n ) * mxaapq * mxsinj < tol ) {
      goto1994 = true;
      break;
    }
    if ( notrot >= emptsw ) {
      goto1994 = true;
      break;
    }
  } // 1993
  if ( ! goto1994 ) info.setValue( nsweep - 1 );
  else info.setValue( 0 );
  for ( p = 1; p <= n - 1; p ++ ) { // 5991
    q = Blas1.idamax( n - p + 1, sva, 1, ioffsva + p - 1 ) + p - 1;
    if ( p != q ) {
      temp1 = sva[ ioffsva + p - 1 ];
      sva[ ioffsva + p - 1 ] = sva[ ioffsva + q - 1 ];
      sva[ ioffsva + q - 1 ] = temp1;
      temp1 = d[ ioffd + p - 1 ];
      d[ ioffd + p - 1 ] = d[ ioffd + q - 1 ];
      d[ ioffd + q - 1 ] = temp1;
      Blas1.dswap( m, A, 1, A, 1, ioffa + ( p - 1 ) * lda,
        ioffa + ( q - 1 ) * lda );
      if ( rsvec ) {
        Blas1.dswap( mvl, V, 1, V, 1, ioffv + ( p - 1 ) * ldv,
          ioffv + ( q - 1 ) * ldv );
      }
    }
  } // 5991
}
//************************************************************************
LaPack2.dgtcon = function( norm, n, dl, d, du, du2, ipiv,
anorm, rcond, work, iwork, info, ioffdl, ioffd, ioffdu, ioffdu2,
ioffipiv, ioffwork, ioffiwork) {
  var isave = new Array( 3 );
  info.setValue( 0 );
  var onenrm = ( norm.charAt(0).toUpperCase() == '1' ||
    norm.charAt(0).toUpperCase() == 'O' );
  if ( ! onenrm && norm.charAt(0).toUpperCase() != 'I' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( anorm < 0. ) info.setValue( -8 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgtcon', - info.getValue() );
    return;
  }
  rcond.setValue( 0. );
  if ( n == 0 ) {
    rcond.setValue( 1. );
    return;
  } else if ( anorm == 0. ) return;
  for ( var i = 1; i <= n; i ++ ) {
    if ( d[ ioffd + i - 1 ] == 0. ) return;
  }
  var ainvnm = new NumberReference( 0. );
  var kase1 = ( onenrm ? 1 : 2 );
  var kase = new IntReference( 0 );
  while ( true ) {
    LaPack0.dlacn2( n, work, work, iwork, ainvnm, kase, isave,
      ioffwork + n, ioffwork, ioffiwork, 0 );
    if ( kase.getValue() == 0 ) break;
    if ( kase.getValue() == kase1 ) {
      LaPack1.dgttrs( 'No transpose', n, 1, dl, d, du, du2, ipiv, work,
        n, info, ioffdl, ioffd, ioffdu, ioffdu2, ioffipiv, ioffwork );
    } else {
      LaPack1.dgttrs( 'transpose', n, 1, dl, d, du, du2, ipiv, work,
        n, info, ioffdl, ioffd, ioffdu, ioffdu2, ioffipiv, ioffwork );
    }
  }
  if ( ainvnm.getValue() != 0. ) {
    rcond.setValue( ( 1. / ainvnm.getValue() ) / anorm );
  }
}
LaPack2.zgtcon = function( norm, n, dl, d, du, du2, ipiv,
anorm, rcond, work, iwork, info ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dgtrfs = function( trans, n, nrhs, dl, d, du, dlf, df,
duf, du2, ipiv, B, ldb, X, ldx, ferr, berr, work, iwork, info, ioffdl,
ioffd, ioffdu, ioffdlf, ioffdf, ioffduf, ioffdu2, ioffipiv, ioffb, ioffx,
ioffferr, ioffberr, ioffwork, ioffiwork ) {
  var itmax = 5;
  var isave = new Array( 3 );
  info.setValue( 0 );
  var notran = ( trans.charAt(0).toUpperCase() == 'N' );
  if ( ! notran && trans.charAt(0).toUpperCase() != 'T'
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( nrhs < 0 ) info.setValue( -3 );
  else if ( ldb < Math.max( 1, n ) ) info.setValue( -13 );
  else if ( ldx < Math.max( 1, n ) ) info.setValue( -15 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgtrfs', - info.getValue() );
    return;
  }
  if ( n == 0 || nrhs == 0 ) {
    for ( var j = 1; j <= nrhs; j ++ ) {
      ferr[ ioffferr + j - 1 ] = 0.;
      berr[ ioffberr + j - 1 ] = 0.;
    }
    return;
  }
  if ( notran ) {
    var transn = 'N';
    var transt = 'T';
  } else {
    transn = 'T';
    transt = 'N';
  }
  var nz = 4;
  var eps = LaPack0.dlamch( 'Epsilon' );
  var safmin = LaPack0.dlamch( 'Safe minimum' );
  var safe1 = nz * safmin;
  var safe2 = safe1 / eps;
  for ( j = 1; j<= nrhs; j ++ ) {
    var count = 1;
    var lstres = 3.;
    while ( true ) { // 20
      Blas1.dcopy( n, B, 1, work, 1, ioffb + ( j - 1 ) * ldb,
        ioffwork + n );
      LaPack0.dlagtm( trans, n, 1, -1., dl, d, du, X, ldx, 1., work, n,
        ioffdl, ioffd, ioffdu, ioffx + ( j - 1 ) * ldx, ioffwork + n );
      if ( notran ) {
        if ( n == 1 ) {
          work[ ioffwork ] = Math.abs( B[ ioffb + ( j - 1 ) * ldb ] )
            + Math.abs( d[ ioffd ] * X[ ioffx + ( j - 1 ) * ldx ] );
        } else {
          work[ ioffwork ] = Math.abs( B[ ioffb + ( j - 1 ) * ldb ] )
            + Math.abs( d[ ioffd ] * X[ ioffx + ( j - 1 ) * ldx ] );
            + Math.abs( du[ ioffdu ]
            * X[ ioffx + 1 + ( j - 1 ) * ldx ] );
          for ( var i = 2; i <= n - 1; i ++ ) {
            work[ ioffwork + i - 1 ] =
              Math.abs( B[ ioffb + i - 1 + ( j - 1 ) * ldb ] )
              + Math.abs( dl[ ioffdl + i - 2 ]
              * X[ ioffx + i - 2 + ( j - 1 ) * ldx ] )
              + Math.abs( d[ ioffd + i - 1 ]
              * X[ ioffx + i - 1 + ( j - 1 ) * ldx ] )
              + Math.abs( du[ ioffdu + i - 1 ]
              * X[ ioffx + i + ( j - 1 ) * ldx ] );
          }
          work[ ioffwork + n - 1 ] =
            Math.abs( B[ ioffb + n - 1 + ( j - 1 ) * ldb ] )
            + Math.abs( dl[ ioffdl + n - 2 ]
            * X[ ioffx + n - 2 + ( j - 1 ) * ldx ] )
            + Math.abs( d[ ioffd + n - 1 ]
            * X[ ioffx + n - 1 + ( j - 1 ) * ldx ] );
        }
      } else {
        if ( n == 1 ) {
          work[ ioffwork ] = Math.abs( B[ ioffb + ( j - 1 ) * ldb ] )
            + Math.abs( d[ ioffd ] * X[ ioffx + ( j - 1 ) * ldx ] );
        } else {
          work[ ioffwork ] = Math.abs( B[ ioffb + ( j - 1 ) * ldb ] )
            + Math.abs( d[ ioffd ] * X[ ioffx + ( j - 1 ) * ldx ] );
            + Math.abs( dl[ ioffdl ]
            * X[ ioffx + 1 + ( j - 1 ) * ldx ] );
          for ( i = 2; i <= n - 1; i ++ ) {
            work[ ioffwork + i - 1 ] =
              Math.abs( B[ ioffb + i - 1 + ( j - 1 ) * ldb ] )
              + Math.abs( du[ ioffdu + i - 2 ]
              * X[ ioffx + i - 2 + ( j - 1 ) * ldx ] )
              + Math.abs( d[ ioffd + i - 1 ]
              * X[ ioffx + i - 1 + ( j - 1 ) * ldx ] )
              + Math.abs( dl[ ioffdl + i - 1 ]
              * X[ ioffx + i + ( j - 1 ) * ldx ] );
          }
          work[ ioffwork + n - 1 ] =
            Math.abs( B[ ioffb + n - 1 + ( j - 1 ) * ldb ] )
            + Math.abs( du[ ioffdu + n - 2 ]
            * X[ ioffx + n - 2 + ( j - 1 ) * ldx ] )
            + Math.abs( d[ ioffd + n - 1 ]
            * X[ ioffx + n - 1 + ( j - 1 ) * ldx ] );
        }
      }
      var s = 0.;
      for ( i = 1; i <= n; i ++ ) {
        if ( work[ ioffwork + i - 1 ] > safe2 ) {
          s = Math.max( s,
            Math.abs( work[ ioffwork + n + i - 1 ] )
            / work[ ioffwork + i - 1 ] );
        } else {
          s = Math.max( s,
              ( Math.abs( work[ ioffwork + n + i - 1 ] ) + safe1 )
              / ( work[ ioffwork + i - 1 ] + safe1 ) );
        }
      }
      berr[ ioffberr + j - 1 ] = s;
      if ( berr[ ioffberr + j - 1 ] > eps &&
      2. * berr[ ioffberr + j - 1 ] <= lstres &&
      count <= itmax ) {
        LaPack1.dgttrs( trans, n, 1, dlf, df, duf, du2, ipiv, work, n,
          info, ioffdlf, ioffdf, ioffduf, ioffdu2, ioffipiv,
          ioffwork + n );
        Blas1.daxpy( n, 1., work, 1, X, 1, ioffwork + n,
          ioffx + (j - 1 ) * ldx );
        lstres = berr[ ioffberr + j - 1 ];
        count ++;
      } else break;
    }
    for ( i = 1; i <= n; i ++ ) {
      if ( work[ ioffwork + i - 1 ] > safe2 ) {
        work[ ioffwork + i - 1 ] =
          Math.abs( work[ ioffwork + n + i - 1 ] )
          + nz * eps * work[ ioffwork + i - 1 ];
      } else {
        work[ ioffwork + i - 1 ] =
          Math.abs( work[ ioffwork + n + i - 1 ] )
          + nz * eps * work[ ioffwork + i - 1 ] + safe1;
      }
    }
    var kase = new IntReference( 0 );
    while ( true ) { // 70
      var est =
        new NumberReference( ferr[ ioffferr + j - 1 ] );
      LaPack0.dlacn2( n, work, work, iwork, est, kase, isave,
        ioffwork + 2*n, ioffwork + n, ioffiwork, 0 );
      ferr[ ioffferr + j - 1 ] = est.getValue();
      if ( kase.getValue() == 0 ) break;
      if ( kase.getValue() == 1 ) {
        LaPack1.dgttrs( transt, n, 1, dlf, df, duf, du2, ipiv, work,
          n, info, ioffdlf, ioffdf, ioffduf, ioffdu2, ioffipiv,
          ioffwork + n );
        for ( i = 1; i <= n; i ++ ) {
          work[ ioffwork + n + i - 1 ] *= work[ ioffwork + i - 1 ];
        }
      } else {
        for ( i = 1; i <= n; i ++ ) {
          work[ ioffwork + n + i - 1 ] *= work[ ioffwork + i - 1 ];
        }
        LaPack1.dgttrs( transn, n, 1, dlf, df, duf, du2, ipiv, work,
          n, info, ioffdlf, ioffdf, ioffduf, ioffdu2, ioffipiv,
          ioffwork + n );
      }
    }
    lstres = 0.;
    for ( i = 1; i <= n; i ++ ) {
      lstres = Math.max( lstres,
        Math.abs( X[ ioffx + i - 1 + ( j - 1 ) * ldx ] ) );
    }
    if ( lstres != 0. ) ferr[ ioffferr + j - 1 ] /= lstres;
  }
}
LaPack2.zgtrfs = function( trans, n, nrhs, dl, d, du, dlf, df,
duf, du2, ipiv, B, ldb, X, ldx, ferr, berr, work, iwork, info ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dhgeqz = function( job, compq, compz, n, ilo, ihi, A,
lda, B, ldb, alphar, alphai, beta, Q, ldq, Z, ldz, work, lwork, info ) {
  throw new Error("not programmed: generalized eigenvalue problem");
}
//************************************************************************
LaPack2.dla_gercond = function( trans, n, A, lda, AF, ldaf,
ipiv, cmode, c, info, work, iwork, ioffa, ioffaf, ioffipiv, ioffc,
ioffwork, ioffiwork ) {
  throw new Error("not tested");
  var isave = new Array( 3 );
  info.setValue( 0 );
  var notrans = ( trans.charAt(0).toUpperCase() == 'N' );
  if ( ! notrans && trans.charAt(0).toUpperCase() != 'T'
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -4 );
  else if ( ldaf < Math.max( 1, n ) ) info.setValue( -6 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dla_gercond', - info.getValue() );
    return 0.;
  }
  if ( n == 0 ) return 1.;
  if ( notrans ) {
    for ( var i = 1; i <= n; i ++ ) {
      var tmp = 0.;
      if ( cmode == 1 ) {
        for ( var j = 1; j <= n; j ++ ) {
          tmp += Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ]
            * c[ ioffc + j - 1 ] );
        }
      } else if ( cmode == 0 ) {
        for ( j = 1; j <= n; j ++ ) {
          tmp += Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          tmp += Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ]
            / c[ ioffc + j - 1 ] );
        }
      }
      work[ ioffwork + 2 * n + i - 1 ] = tmp;
    }
  } else {
    for ( i = 1; i <= n; i ++ ) {
      tmp = 0.;
      if ( cmode == 1 ) {
        for ( j = 1; j <= n; j ++ ) {
          tmp += Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ]
            * c[ ioffc + j - 1 ] );
        }
      } else if ( cmode == 0 ) {
        for ( j = 1; j <= n; j ++ ) {
          tmp += Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ] );
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          tmp += Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ]
            / c[ ioffc + j - 1 ] );
        }
      }
      work[ ioffwork + 2 * n + i - 1 ] = tmp;
    }
  }
  var ainvnm = new NumberReference( 0. );
  var kase = new IntReference( 0 );
  while ( true ) { // 10
    LaPack0.dlacn2( n, work, work, iwork, ainvnm, kase, isave,
      ioffwork + n, ioffwork, ioffiwork, 0 );
    if ( kase.getValue() != 0 ) {
      if ( kase.getValue() == 2 ) {
        for ( i = 1; i <= n; i ++ ) {
          work[ ioffwork + i - 1 ] *= work[ ioffwork + 2 * n + i - 1 ];
        }
        if ( notrans ) {
          LaPack1.dgetrs( 'No transpose', n, 1, AF, ldaf, ipiv, work,
            n, info, ioffaf, ioffipiv, ioffwork );
        } else {
          LaPack1.dgetrs( 'Transpose', n, 1, AF, ldaf, ipiv, work,
            n, info, ioffaf, ioffipiv, ioffwork );
        }
        if ( cmode == 1 ) {
          for ( i = 1; i <= n; i ++ ) {
            work[ ioffwork + i - 1 ] /= c[ ioffc + i - 1 ];
          }
        } else if ( cmode == -1 ) {
          for ( i = 1; i <= n; i ++ ) {
            work[ ioffwork + i - 1 ] *= c[ ioffc + i - 1 ];
          }
        }
      } else {
        if ( cmode == 1 ) {
          for ( i = 1; i <= n; i ++ ) {
            work[ ioffwork + i - 1 ] /= c[ ioffc + i - 1 ];
          }
        } else if ( cmode == -1 ) {
          for ( i = 1; i <= n; i ++ ) {
            work[ ioffwork + i - 1 ] *= c[ ioffc + i - 1 ];
          }
        }
        if ( notrans ) {
          LaPack1.dgetrs( 'Transpose', n, 1, AF, ldaf, ipiv, work,
            n, info, ioffaf, ioffipiv, ioffwork );
        } else {
          LaPack1.dgetrs( 'No transpose', n, 1, AF, ldaf, ipiv, work,
            n, info, ioffaf, ioffipiv, ioffwork );
        }
        for ( i = 1; i <= n; i ++ ) {
          work[ ioffwork + i - 1 ] *= work[ ioffwork + 2 * n + i - 1 ];
        }
      }
    } else break;
  }
  if ( ainvnm.getValue() != 0. ) return 1. / ainvnm.getValue();
  else return 0.;
}
//************************************************************************
LaPack2.dlabrd = function( m, n, nb, A, lda, d, e, tauq, taup,
X, ldx, Y, ldy, ioffa, ioffd, ioffe, iofftauq, iofftaup, ioffx, ioffy) {
  throw new Error("not tested");
  if ( m <= 0 || n <= 0 ) return;
  if ( m >= n ) {
    for ( var i = 1; i <= nb; i ++ ) {
      Blas2.dgemv( 'No transpose', m - i + 1, i - 1, -1., A, lda,
        Y, ldy, 1., A, 1, ioffa + i - 1, ioffy + i - 1,
        ioffa + i - 1 + ( i - 1 ) * lda );
      Blas2.dgemv( 'No transpose', m - i + 1, i - 1, -1., X, ldx,
        A, 1, 1., A, 1, ioffx + i - 1, ioffa + ( i - 1 ) * lda,
        ioffa + i - 1 + ( i - 1 ) * lda );
      var alpha =
        new NumberReference( A[ ioffa + i - 1 + ( i - 1 ) * lda ] );
      var tau =
        new NumberReference( tauq[ iofftauq + i - 1 ] );
      LaPack1.dlarfg( m - i + 1, alpha, A, 1, tau,
        ioffa + Math.min( i + 1, m ) - 1 + ( i - 1 ) * lda );
      A[ ioffa + i - 1 + ( i - 1 ) * lda ] = alpha.getValue();
      tauq[ iofftauq + i - 1 ] = tau.getValue();
      d[ ioffd + i - 1 ] = A[ ioffa + i - 1 + ( i - 1 ) * lda ];
      if ( i < n ) {
        A[ ioffa + i - 1 + ( i - 1 ) * lda ]= 1.;
        Blas2.dgemv( 'Transpose', m - i + 1, n - i, 1., A, lda, A, 1,
          0., Y, 1, ioffa + i - 1 + i * lda,
          ioffa + i - 1 + ( i - 1 ) * lda,
          ioffy + i + ( i - 1 ) * ldy );
        Blas2.dgemv( 'Transpose', m - i + 1, i - 1, 1., A, lda, A, 1,
          0., Y, 1, ioffa + i - 1, ioffa + i - 1 + ( i - 1 ) * lda,
          ioffy + ( i - 1 ) * ldy );
        Blas2.dgemv( 'No transpose', n - i, i - 1, -1., Y, ldy, Y, 1,
          1., Y, 1, ioffy + i, ioffy + ( i - 1 ) * ldy,
          ioffy + i + ( i - 1 ) * ldy );
        Blas2.dgemv( 'Transpose', m - i + 1, i - 1, 1., X, ldx, A, 1,
          0., Y, 1, ioffx + i - 1, ioffa + i - 1 + ( i - 1 ) * lda,
          ioffy + ( i - 1 ) * ldy );
        Blas2.dgemv( 'Transpose', i - 1, n - i, -1., A, lda, Y, 1,
          1., Y, 1, ioffa + i * lda, ioffy + ( i - 1 ) * ldy,
          ioffy + i + ( i - 1 ) * ldy );
        Blas1.dscal( n - i, tauq[ iofftauq + i- 1 ], Y, 1,
          ioffy + i + ( i - 1 ) * ldy );
        Blas2.dgemv( 'No transpose', n - i, i, -1., Y, ldy, A, lda,
          1., A, lda, ioffy + i, ioffa + ( i - 1 ) * lda,
          ioffa + i - 1 + i * lda );
        Blas2.dgemv( 'Transpose', i - 1, n - i, -1., A, lda, X, ldx,
          1., A, lda, ioffa + i * lda, ioffx + i - 1,
          ioffa + i - 1 + i * lda );
        alpha.setValue( A[ ioffa + i - 1 + i * lda ] );
        tau.setValue( taup[ iofftaup + i - 1 ] );
        LaPack1.dlarfg( n - i, alpha, A, lda, tau,
          ioffa + i - 1 + ( Math.min( i + 2, n ) - 1 ) * lda );
        A[ ioffa + i - 1 + i * lda ] = alpha.getValue();
        taup[ iofftaup + i - 1 ] = tau.getValue();
        e[ ioffe + i - 1 ] = A[ ioffa + i - 1 + i * lda ];
        A[ ioffa + i - 1 + i * lda ] = 1.;
        Blas2.dgemv( 'No transpose', m - i, n - i, 1., A, lda, A, lda,
          0., X, 1, ioffa + i + i * lda, ioffa + i - 1 + i * lda,
          ioffx + i + ( i - 1 ) * ldx );
        Blas2.dgemv( 'Transpose', n - i, i, 1., Y, ldy, A, lda,
          0., X, 1, ioffy + i, ioffa + i - 1 + i * lda,
          ioffx + ( i - 1 ) * ldx );
        Blas2.dgemv( 'No transpose', m - i, i, -1., A, lda, X, 1,
          1., X, 1, ioffa + i, ioffx + ( i - 1 ) * ldx,
          ioffx + i + ( i - 1 ) * ldx );
        Blas2.dgemv( 'No transpose', i - 1, n - i, 1., A, lda, A, lda,
          0., X, 1, ioffa + i * lda, ioffa + i - 1 + i * lda,
          ioffx + ( i - 1 ) * ldx );
        Blas2.dgemv( 'No transpose', m - i, i - 1, -1., X, ldx, X, 1,
          1., X, 1, ioffx + i * ldx, ioffx + ( i - 1 ) * ldx,
          ioffx + i + ( i - 1 ) * ldx );
        Blas1.dscal( m - i, taup[ iofftaup + i - 1 ], X, 1,
          ioffx + i + ( i - 1 ) * ldx );
      }
    }
  } else {
    for ( i = 1; i <= nb; i ++ ) {
      Blas2.dgemv( 'No transpose', n - i + 1, i - 1, -1., Y, ldy,
        A, lda, 1., A, lda, ioffy + i - 1, ioffa + i - 1,
        ioffa + i - 1 + ( i - 1 ) * lda );
      Blas2.dgemv( 'Transpose', i - 1, n - i + 1, -1., A, lda,
        X, ldx, 1., A, lda, ioffa + ( i - 1 ) * lda,
        ioffx + i - 1, ioffa + i - 1 + ( i - 1 ) * lda );
      alpha.setValue( A[ ioffa + i - 1 + ( i - 1 ) * lda ] );
      tau.setValue( taup[ iofftaup + i - 1 ] );
      LaPack1.dlarfg( n - i + 1, alpha, A, lda, tau,
        ioffa + i - 1 + ( Math.min( i + 1, n ) - 1 ) * lda );
      A[ ioffa + i - 1 + ( i - 1 ) * lda ] = alpha.getValue();
      taup[ iofftaup + i - 1 ] = tau.getValue();
      d[ ioffd + i - 1 ] = A[ ioffa + i - 1 + ( i - 1 ) * lda ];
      if ( i < m ) {
        A[ ioffa + i - 1 + ( i - 1 ) * lda ]= 1.;
        Blas2.dgemv( 'No transpose', m - i, n - i + 1, 1., A, lda,
          A, lda, 0., X, 1, ioffa + i + ( i - 1 ) * lda,
          ioffa + i - 1 + ( i - 1 ) * lda,
          ioffx + i + ( i - 1 ) * ldx );
        Blas2.dgemv( 'Transpose', n - i + 1, i - 1, 1., Y, ldy, A, lda,
          0., X, 1, ioffy + i - 1, ioffa + i - 1 + ( i - 1 ) * lda,
          ioffx + ( i - 1 ) * ldx );
        Blas2.dgemv( 'No transpose', m - i, i - 1, -1., A, lda, X, 1,
          1., X, 1, ioffa + i, ioffx + ( i - 1 ) * ldx,
          ioffx + i + ( i - 1 ) * ldx );
        Blas2.dgemv( 'No transpose', i - 1, n - i + 1, 1., A, lda,
          A, lda, 0., X, 1, ioffa + ( i - 1 ) * lda,
          ioffa + i - 1 + ( i - 1 ) * lda, ioffx + ( i - 1 ) * ldx );
        Blas2.dgemv( 'No transpose', m - i, i - 1, -1., X, ldx, X, 1,
          1., X, 1, ioffx + i, ioffx + ( i - 1 ) * ldx,
          ioffx + i + ( i - 1 ) * ldx );
        Blas1.dscal( m - i, taup[ iofftaup + i - 1 ], X, 1,
          ioffx + i + ( i - 1 ) * ldx );
        Blas2.dgemv( 'No transpose', m - i, i - 1, -1., A, lda, Y, ldy,
          1., A, 1, ioffa + i, ioffy + i - 1,
          ioffa + i + ( i - 1 ) * lda );
        Blas2.dgemv( 'No transpose', m - i, i, -1., X, ldx, A, 1,
          1., A, 1, ioffx + i, ioffa + ( i - 1 ) * lda,
          ioffa + i + ( i - 1 ) * lda );
        alpha.setValue( A[ ioffa + i + ( i - 1 ) * lda ] );
        tau.setValue( tauq[ iofftauq + i - 1 ] );
        LaPack1.dlarfg( m - i, alpha, A, 1, tau,
          ioffa + Math.min( i + 2, m ) - 1 + ( i - 1 ) * lda );
        A[ ioffa + i + ( i - 1 ) * lda ] = alpha.getValue();
        tauq[ iofftauq + i - 1 ] = tau.getValue();
        e[ ioffe + i - 1 ] = A[ ioffa + i + ( i - 1 ) * lda ];
        A[ ioffa + i + ( i - 1 ) * lda ] = 1.;
        Blas2.dgemv( 'Transpose', m - i, n - i, 1., A, lda, A, 1,
          0., Y, 1, ioffa + i + i * lda, ioffa + i + ( i - 1 ) * lda,
          ioffy + i + ( i - 1 ) * ldy );
        Blas2.dgemv( 'Transpose', m - i, i - 1, 1., A, lda, A, 1,
          0., Y, 1, ioffa + i, ioffa + i + ( i - 1 ) * lda,
          ioffy + ( i - 1 ) * ldy );
        Blas2.dgemv( 'No transpose', n - i, i - 1, -1., Y, ldy, Y, 1,
          1., Y, 1, ioffy + i, ioffy + ( i - 1 ) * ldy,
          ioffy + i + ( i - 1 ) * ldy );
        Blas2.dgemv( 'Transpose', m - i, i, 1., X, ldx, A, 1,
          0., Y, 1, ioffx + i, ioffa + i + ( i - 1 ) * lda,
          ioffy + ( i - 1 ) * ldy );
        Blas2.dgemv( 'Transpose', i, n - i, -1., A, lda, Y, 1,
          1., Y, 1, ioffa + i * lda, ioffy + ( i - 1 ) * ldy,
          ioffy + i + ( i - 1 ) * ldy );
        Blas1.dscal( n - i, tauq[ iofftauq + i - 1 ], Y, 1,
          ioffy + i + ( i - 1 ) * ldy );
      }
    }
  }
}
//************************************************************************
LaPack2.dlaed4 = function( n, i, d, z, delta, rho, dlam, info,
ioffd, ioffz, ioffdelta ) {
//document.getElementById("debug_textarea").value +=
//  "entering dlaed4, n,i,rho = " + n + " " + i + " " + rho + "\n";
//document.getElementById("debug_textarea").value += "d = ";
//for ( var jj = 0; jj < n; jj ++ ) {
//  var di = d[ ioffd + jj ];
//  document.getElementById("debug_textarea").value += di + " ";
//}
//document.getElementById("debug_textarea").value += "\n";
//document.getElementById("debug_textarea").value += "z = ";
//for ( var jj = 0; jj < n; jj ++ ) {
//  var zi = z[ ioffz + jj ];
//  document.getElementById("debug_textarea").value += zi + " ";
//}
//document.getElementById("debug_textarea").value += "\n";

  var zz = new Array(3);
  var maxit = 30;
  info.setValue( 0 );
  if ( n == 1 ) {
    dlam.setValue( d[ ioffd ] + rho * z[ ioffz ] * z[ ioffz ] );
    delta[ ioffdelta ] = 1.;
//  document.getElementById("debug_textarea").value +=
//    "n==1, dlam = " + dlam.getValue() + "\n";
    return;
  }
  if ( n == 2 ) {
    LaPack0.dlaed5( i, d, z, delta, rho, dlam, ioffd, ioffz, ioffdelta );
//  document.getElementById("debug_textarea").value +=
//    "n==2, dlam = " + dlam.getValue() + "\n";
//  document.getElementById("debug_textarea").value += "delta = ";
//  for ( var jj = 0; jj < n; jj ++ ) {
//    var deltai = delta[ ioffdelta + jj ];
//    document.getElementById("debug_textarea").value += deltai + " ";
//  }
//  document.getElementById("debug_textarea").value += "\n";
    return;
  }
  var eps = LaPack0.dlamch( 'Epsilon' );
  var rhoinv = 1. / rho;
  if ( i == n ) {
//  document.getElementById("debug_textarea").value += "i==n\n";
    var ii = n - 1;
    var niter = 1;
    var midpt = rho / 2.;

//  document.getElementById("debug_textarea").value +=
//    "ii,niter,midpt = " + ii + " " + niter + " " + midpt + "\n";

    for ( var j = 1; j <= n; j ++ ) { // loop 10
      delta[ ioffdelta + j - 1 ] =
        ( d[ ioffd + j - 1 ] - d[ ioffd + i - 1 ] ) - midpt;
    }

//  document.getElementById("debug_textarea").value += "L10 :delta = ";
//  for ( var jj = 0; jj < n; jj ++ ) {
//    var deltai = delta[ ioffdelta + jj ];
//    document.getElementById("debug_textarea").value += deltai + " ";
//  }
//  document.getElementById("debug_textarea").value += "\n";

    var psi = 0.;
    for ( j = 1; j <= n - 2; j ++ ) { // loop 20
      psi += z[ ioffz + j - 1 ] * z[ ioffz + j - 1 ]
        / delta[ ioffdelta + j - 1 ];
    }
    var c = rhoinv + psi;
    var w = c + z[ ioffz + ii - 1 ] * z[ ioffz + ii - 1 ]
      / delta[ ioffdelta + ii - 1 ]
      + z[ ioffz + n - 1 ] * z[ ioffz + n - 1 ]
      / delta [ ioffdelta + n - 1 ];
//  document.getElementById("debug_textarea").value +=
//    "L20 : psi,c,w = " + psi + " " + c + " " + w + "\n";
    if ( w <= 0. ) {
      var temp = z[ ioffz + n - 2 ] * z[ ioffz + n - 2 ]
        / ( d[ ioffd + n - 1 ] - d[ ioffd + n - 2 ] + rho )
        + z[ ioffz + n - 1 ] * z[ ioffz + n - 1 ] / rho;
      if ( c <= temp ) var tau = rho;
      else {
        var del = d[ ioffd + n - 1 ] - d[ ioffd + n - 2 ];
        var a = - c * del
          + z[ ioffz + n - 2 ] * z[ ioffz + n - 2 ]
          + z[ ioffz + n - 1 ] * z[ ioffz + n - 1 ];
        var b = z[ ioffz + n - 1 ] * z[ ioffz + n - 1 ] * del;
        tau = ( a < 0. ?
          2. * b / ( Math.sqrt( a * a + 4. * b * c ) - a ) :
          ( a + Math.sqrt( a * a + 4. * b * c ) ) / ( 2. * c ) );
      }
      var dltlb = midpt;
      var dltub = rho;
//    document.getElementById("debug_textarea").value +=
//      "w < 0 : temp,tau,del,a,b,dltlb,dltub = " + temp + " " + tau
//      + " " + del + " " + a + " " + b + " " + dltlb + " " + dltub
//      + "\n";
    } else {
      del = d[ ioffd + n - 1 ] - d[ ioffd + n - 2 ];
      a = - c * del + z[ ioffz + n - 2 ] * z[ ioffz + n - 2 ]
        + z[ ioffz + n - 1 ] * z[ ioffz + n - 1 ];
      b = z[ ioffz + n - 1 ] * z[ ioffz + n - 1 ] * del;
      tau = ( a < 0. ? 2. * b / ( Math.sqrt( a * a + 4. * b * c ) - a )
        : ( a + Math.sqrt( a * a + 4. * b * c ) ) / ( 2. * c ) );
      dltlb = 0.;
      dltub = midpt;
//    document.getElementById("debug_textarea").value +=
//      "w >= 0 : del,a,b,tau,dltlb,dltub = " + del + " " + a + " " + b
//      + " " + tau + " " + dltlb + " " + dltub + "\n"
    }
    for ( j = 1; j <= n; j ++ ) { // loop 30
      delta[ ioffdelta + j - 1 ] =
        ( d[ ioffd + j - 1 ] - d[ ioffd + i - 1 ] ) - tau;
    }

//  document.getElementById("debug_textarea").value += "L30 : delta = ";
//  for ( var jj = 0; jj < n; jj ++ ) {
//    var deltai = delta[ ioffdelta + jj ];
//    document.getElementById("debug_textarea").value += deltai + " ";
//  }
//  document.getElementById("debug_textarea").value += "\n";

    var dpsi = 0.;
    psi = 0.;
    var erretm = 0.;
    for ( j = 1; j <= ii; j ++ ) { // loop 40
      temp = z[ ioffz + j - 1 ] / delta[ ioffdelta + j - 1 ];
      psi += z[ ioffz + j - 1 ] * temp;
      dpsi += temp * temp;
      erretm += psi;
    }
    erretm = Math.abs( erretm );

//  document.getElementById("debug_textarea").value +=
//    "L40 : psi,dpsi,erretm = " + psi + " " + dpsi + " " + erretm + "\n";

    temp = z[ ioffz + n - 1 ] / delta[ ioffdelta + n - 1 ];
    var phi = z[ ioffz + n - 1 ] * temp;
    var dphi = temp * temp;
    erretm = 8. * ( - phi - psi ) + erretm - phi + rhoinv
      + Math.abs( tau ) * ( dpsi + dphi );
    w = rhoinv + phi + psi;

//  document.getElementById("debug_textarea").value +=
//    "temp,phi,dphi,erretm,w = " + temp + " " + phi + " " + dphi + " "
//    + erretm + " " + w + "\n";

    if ( Math.abs( w ) <= eps * erretm ) {
      dlam.setValue( d[ ioffd + i - 1 ] + tau );
//    document.getElementById("debug_textarea").value +=
//      "| w | <= eps * erretm\n";
      return;
    }
    if ( w <= 0. ) dltlb = Math.max( dltlb, tau);
    else dltub = Math.min( dltub, tau );
    niter ++;
    c = w - delta[ ioffdelta + n - 2 ] * dpsi
      - delta[ ioffdelta + n - 1 ] * dphi;
    a = ( delta[ ioffdelta + n - 2 ]
      + delta[ ioffdelta + n - 1 ] ) * w
      - delta[ ioffdelta + n - 2 ] * delta[ ioffdelta + n - 1 ]
      * ( dpsi + dphi );
    b = delta[ ioffdelta + n - 2 ] * delta[ ioffdelta + n - 1 ] * w;
    if ( c < 0. ) c = Math.abs( c );
    if ( c == 0. ) var eta = dltub - tau;
    else if ( a >= 0. ) {
      eta = ( a + Math.sqrt( Math.abs( a * a - 4. * b * c ) ) )
        / ( 2. * c );
    } else {
      eta = 2. * b
        / ( a - Math.sqrt( Math.abs( a * a - 4. * b * c ) ) );
    }

//  document.getElementById("debug_textarea").value +=
//    "dltlb,dltub,niter,c,a,b,eta = " + dltlb + " " + dltub + " "
//    + niter + " " + c + " " + a + " " + b + " " + eta + "\n";

    if ( w * eta > 0. ) eta = -w / ( dpsi + dphi );
    temp = tau + eta;
    if ( temp > dltub || temp < dltlb ) {
      eta = ( w < 0. ? ( dltub - tau ) / 2. : ( dltlb - tau ) / 2. );
    }
    for ( j = 1; j <= n; j ++ ) delta[ ioffdelta + j - 1 ] -= eta; // L50
    tau += eta;

//  document.getElementById("debug_textarea").value +=
//    "L50 : eta,temp,tau = " + eta + " " + temp + " " + tau + "\n";
//  document.getElementById("debug_textarea").value += "delta = ";
//  for ( var jj = 0; jj < n; jj ++ ) {
//    var deltai = delta[ ioffdelta + jj ];
//    document.getElementById("debug_textarea").value += deltai + " ";
//  }
//  document.getElementById("debug_textarea").value += "\n";

    dpsi = 0.
    psi = 0.;
    erretm = 0.;
    for ( j = 1; j <= ii; j ++ ) { // loop 60
      temp = z[ ioffz + j - 1 ] / delta[ ioffdelta + j - 1 ];
      psi += z[ ioffz + j - 1 ] * temp;
      dpsi += temp * temp;
      erretm += psi;
    }
    erretm = Math.abs( erretm );

//  document.getElementById("debug_textarea").value +=
//    "L60 : dpsi,psi,erretm = " + dpsi + " " + psi + " " + erretm + "\n";

    temp = z[ ioffz + n - 1 ] / delta[ ioffdelta + n - 1 ];
    phi = z[ ioffz + n - 1 ] * temp;
    dphi = temp * temp;
    erretm += 8. * ( -phi - psi ) - phi + rhoinv
      + Math.abs( tau ) * ( dpsi + dphi );
    w = rhoinv + phi + psi;
    var iter = niter + 1;

//  document.getElementById("debug_textarea").value +=
//    "temp,phi,dphi,erretm,w,iter = " + temp + " " + phi + " "
//    + dphi + " " + erretm + " " + w + " " + iter + "\n";

    for ( niter = iter; niter <= maxit; niter ++ ) {
//    document.getElementById("debug_textarea").value +=
//      "\n\nniter = " + niter + "\n\n";
      if ( Math.abs( w ) <= eps * erretm ) {
        dlam.setValue( d[ ioffd + i - 1 ] + tau );
//      document.getElementById("debug_textarea").value +=
//        "| w | <= eps * erretm: dlam = " + dlam.getValue() + "\n";
//      document.getElementById("debug_textarea").value += "delta = ";
//      for ( var jj = 0; jj < n; jj ++ ) {
//        var deltai = delta[ ioffdelta + jj ];
//        document.getElementById("debug_textarea").value += deltai + " ";
//      }
//      document.getElementById("debug_textarea").value += "\n";
        return;
      }
      if ( w <= 0. ) dltlb = Math.max( dltlb, tau );
      else dltub = Math.min( dltub, tau );
      c = w - delta[ ioffdelta + n - 2 ] * dpsi
        - delta[ ioffdelta + n - 1 ] * dphi;
      a = ( delta[ ioffdelta + n - 2 ]
        + delta[ ioffdelta + n - 1 ] ) * w
        - delta[ ioffdelta + n - 2 ] * delta[ ioffdelta + n - 1 ]
        * ( dpsi + dphi );
      b = delta[ ioffdelta + n - 2 ] * delta[ ioffdelta + n - 1 ] * w;
      if ( a >= 0. ) {
        eta = ( a + Math.sqrt( Math.abs( a * a - 4. * b * c ) ) )
          / ( 2. * c );
      } else {
        eta = 2. * b
          / ( a - Math.sqrt( Math.abs( a * a - 4. * b * c ) ) );
      }

//    document.getElementById("debug_textarea").value +=
//      "dltlb,dltub,c,a,b,eta = " + dltlb + " " + dltub + " " + c + " "
//      + a + " " + b + " " + eta + "\n";

      if ( w * eta > 0. ) eta = -w / ( dpsi + dphi );
      temp = tau + eta;
      if ( temp > dltub || temp < dltlb ) {
        eta = ( w < 0. ? ( dltub - tau ) / 2. : ( dltlb - tau ) / 2. );
      }
      for ( j = 1; j <= n; j ++ ) delta[ ioffdelta + j - 1 ] -= eta; //L70
      tau += eta;

//    document.getElementById("debug_textarea").value +=
//      "L70 : eta,temp = " + eta + " " + temp + "\n";
//    document.getElementById("debug_textarea").value += "delta = ";
//    for ( var jj = 0; jj < n; jj ++ ) {
//      var deltai = delta[ ioffdelta + jj ];
//      document.getElementById("debug_textarea").value += deltai + " ";
//    }
//    document.getElementById("debug_textarea").value += "\n";

      dpsi = 0.;
      psi = 0.;
      erretm = 0.;
      for ( j = 1; j <= ii; j ++ ) { // loop 80
        temp = z[ ioffz + j - 1 ] / delta[ ioffdelta + j - 1 ];
        psi += z[ ioffz + j - 1 ] * temp;
        dpsi += temp * temp;
        erretm += psi;
      }
      erretm = Math.abs( erretm );

//    document.getElementById("debug_textarea").value +=
//      "L80 : dpsi,psi,erretm = " + dpsi + " " + psi + " " + erretm
//      + "\n";

      temp = z[ ioffz + n - 1 ] / delta[ ioffdelta + n - 1 ];
      phi = z[ ioffz + n - 1 ] * temp;
      dphi = temp * temp;
      erretm += 8. * ( -phi - psi ) - phi + rhoinv
        + Math.abs( tau ) * ( dpsi + dphi );
      w = rhoinv + phi + psi;

//    document.getElementById("debug_textarea").value +=
//      "temp,phi,dphi,erretm,w = " + temp + " " + phi + " " + dphi + " "
//      + erretm + " " + w + "\n";

    }
    info.setValue( 1 );
    dlam.setValue( d[ ioffd + i - 1 ] + tau );

//  document.getElementById("debug_textarea").value +=
//    "dlam = " + dlam.getValue() + "\n";
//  document.getElementById("debug_textarea").value += "delta = ";
//  for ( var jj = 0; jj < n; jj ++ ) {
//    var deltai = delta[ ioffdelta + jj ];
//    document.getElementById("debug_textarea").value += deltai + " ";
//  }
//  document.getElementById("debug_textarea").value += "\n";

    return;
  } else { // i < n
//  document.getElementById("debug_textarea").value += "i<n\n";
    niter = 1;
    var ip1 = i + 1;
    del = d[ ioffd + ip1 - 1 ] - d[ ioffd + i - 1 ];
    midpt = del / 2.;
//  document.getElementById("debug_textarea").value +=
//    "niter,ip1,del,midpt = " + niter + " " + ip1 + " " + del + " "
//    + midpt + "\n";
    for ( j = 1; j <= n; j ++ ) { // loop 100
      delta[ ioffdelta + j - 1 ] =
        ( d[ ioffd + j - 1 ] - d[ ioffd + i- 1 ] ) - midpt;
    }

//  document.getElementById("debug_textarea").value += "L100 : delta = ";
//  for ( var jj = 0; jj < n; jj ++ ) {
//    var deltai = delta[ ioffdelta + jj ];
//    document.getElementById("debug_textarea").value += deltai + " ";
//  }
//  document.getElementById("debug_textarea").value += "\n";

    psi = 0.;
    for ( j = 1; j <= i - 1; j ++ ) { // loop 110
      psi += z[ ioffz + j - 1 ] * z[ ioffz + j - 1 ]
        / delta[ ioffdelta + j - 1 ];
    }
    phi = 0.;
    for ( j = n; j >= i + 2; j -- ) { // loop 120
      phi += z[ ioffz + j - 1 ] * z[ ioffz + j - 1 ]
        / delta[ ioffdelta + j - 1 ];
    }
    c = rhoinv + psi + phi;
    w = c + z[ ioffz + i - 1 ] * z[ ioffz + i - 1 ]
      / delta[ ioffdelta + i - 1 ]
      + z[ ioffz + ip1 - 1 ] * z[ ioffz + ip1 - 1 ]
        / delta [ ioffdelta + ip1 - 1 ];

//  document.getElementById("debug_textarea").value += 
//    "L120 : psi,phi,c,w = " + psi + " " + phi + " " + c + " " + w
//    + "\n";

    if ( w > 0. ) {
//    document.getElementById("debug_textarea").value += "w>0\n"; 
      var orgati = true;
      a = c * del + z[ ioffz + i - 1 ] * z[ ioffz + i - 1 ]
        + z[ ioffz + ip1 - 1 ] * z[ ioffz + ip1 - 1 ];
      b = z[ ioffz + i - 1 ] * z[ ioffz + i - 1 ] * del;
      tau = ( a > 0. ?  tau = 2. * b
          / ( a + Math.sqrt( Math.abs( a * a - 4. * b * c ) ) ) :
        ( a - Math.sqrt( Math.abs( a * a - 4. * b * c ) ) )
          / ( 2. * c ) );
      dltlb = 0.;
      dltub = midpt;
    } else {
//    document.getElementById("debug_textarea").value += "w<=0\n"; 
      orgati = false;
      a = c * del - z[ ioffz + i - 1 ] * z[ ioffz + i - 1 ]
        - z[ ioffz + ip1 - 1 ] * z[ ioffz + ip1 - 1 ];
      b = z[ ioffz + ip1 - 1 ] * z[ ioffz + ip1 - 1 ] * del;
      tau = ( a < 0. ?  2. * b
          / ( a - Math.sqrt( Math.abs( a * a + 4. * b * c ) ) ) :
        - ( a + Math.sqrt( Math.abs( a * a + 4. * b * c ) ) )
          / ( 2. * c ) );
      dltlb = - midpt;
      dltub = 0.;
    }
//  document.getElementById("debug_textarea").value +=
//    "orgati,a,b,tau,dltlb,dltub = " + orgati + " " + a + " " + b + " "
//    + tau + " " + dltlb + " " + dltub + "\n";
    if ( orgati ) {
      for ( j = 1; j <= n; j ++ ) { // loop 130
        delta[ ioffdelta + j - 1 ] =
          ( d[ ioffd + j - 1 ] - d[ ioffd + i - 1 ] ) - tau;
      }
    } else {
      for ( j = 1; j <= n; j ++ ) { // loop 140
        delta[ ioffdelta + j - 1 ] =
          ( d[ ioffd + j - 1 ] - d[ ioffd + ip1 - 1 ] ) - tau;
      }
    }

//  document.getElementById("debug_textarea").value += "L140 : delta = ";
//  for ( var jj = 0; jj < n; jj ++ ) {
//    var deltai = delta[ ioffdelta + jj ];
//    document.getElementById("debug_textarea").value += deltai + " ";
//  }
//  document.getElementById("debug_textarea").value += "\n";

    ii = ( orgati ? i : i + 1 );
    var iim1 = ii - 1;
    var iip1 = ii + 1;
    dpsi = 0.;
    psi = 0.;
    erretm = 0.;
    for ( j = 1; j <= iim1; j ++ ) {
      temp = z[ ioffz + j - 1 ] / delta[ ioffdelta + j - 1 ];
      psi += z[ ioffz + j - 1 ] * temp;
      dpsi += temp * temp;
      erretm += psi;
    }
    erretm = Math.abs( erretm );

//  document.getElementById("debug_textarea").value +=
//    "ii,iim1,iip1,dpsi,psi,erretm = " + ii + " " + iim1 + " " + iip1
//    + " " + dpsi + " " + psi + " " + erretm + "\n";

    dphi = 0.;
    phi = 0.;
    for ( j = n; j >= iip1; j -- ) { // loop 160
      temp = z[ ioffz + j - 1 ] / delta[ ioffdelta + j - 1 ];
      phi += z[ ioffz + j - 1 ] * temp;
      dphi += temp * temp;
      erretm += phi;
    }
    w = rhoinv + phi + psi;

//  document.getElementById("debug_textarea").value +=
//    "L160 : dphi,phi,erretm,w = " + dphi + " " + phi + " " + erretm
//    + " " + w + "\n";

    var swtch3 = false;
    if ( orgati ) {
      if ( w < 0. ) swtch3 = true;
    } else {
      if ( w > 0. ) swtch3 = true;
    }
    if ( ii == 1 || ii == n ) swtch3 = false;
    temp = z[ ioffz + ii - 1 ] / delta[ ioffdelta + ii - 1 ];
    var dw = dpsi + dphi + temp * temp;
    temp *= z[ ioffz + ii - 1 ];
    w += temp;
    erretm += 8. * ( phi - psi ) + 2. * rhoinv
      + 3. * Math.abs( temp ) + Math.abs( tau ) * dw;

//  document.getElementById("debug_textarea").value +=
//    "swtch3,temp,dw,w,erretm = " + swtch3 + " " + temp + " " + dw
//    + " " + w + " " + erretm + "\n";

    if ( Math.abs( w ) <= eps * erretm ) {
      dlam.setValue( ( orgati ? d[ ioffd + i - 1 ] + tau
        : d[ ioffd + ip1 - 1 ] + tau ) );
//    document.getElementById("debug_textarea").value +=
//      "| w | <= eps * erretm : dlam = " + dlam.getValue() + "\n";
      return;
    }
    if ( w <= 0. ) dltlb = Math.max( dltlb, tau );
    else dltub = Math.min( dltub, tau );
    niter ++;

//  document.getElementById("debug_textarea").value +=
//    "dltlb,dltub,niter = " + dltlb + " " + dltub + " " + niter + "\n";

    if ( ! swtch3 ) {
      if ( orgati ) {
        c = w - delta[ ioffdelta + ip1 - 1 ] * dw
          - ( d[ ioffd + i - 1 ] - d[ ioffd + ip1 - 1 ] )
          * Math.pow( z[ ioffz + i - 1 ]
          / delta[ ioffdelta + i - 1 ] , 2 );
      } else {
        c = w - delta[ ioffdelta + i - 1 ] * dw
          - ( d[ ioffd + ip1 - 1 ] - d[ ioffd + i - 1 ] )
          * Math.pow( z[ ioffz + ip1 - 1 ]
          / delta[ ioffdelta + ip1 - 1 ] , 2 );
      }
      a = ( delta[ ioffdelta + i - 1 ]
        + delta[ ioffdelta + ip1 - 1 ] ) * w
        - delta[ ioffdelta + i - 1 ] * delta[ ioffdelta + ip1 - 1 ]
        * dw;
      b = delta[ ioffdelta + i - 1 ] * delta[ ioffdelta + ip1 - 1 ]
        * w;
      if ( c == 0. ) {
        if ( a == 0. ) {
          if ( orgati ) {
            a = z[ ioffz + i - 1 ] * z[ ioffz + i - 1 ]
              + delta[ ioffdelta + ip1 - 1 ]
              * delta[ ioffdelta + ip1 - 1 ] * ( dpsi + dphi );
          } else {
            a = z[ ioffz + ip1 - 1 ] * z[ ioffz + ip1 - 1 ]
              + delta[ ioffdelta + i - 1 ]
              * delta[ ioffdelta + i - 1 ] * ( dpsi + dphi );
          }
        }
        eta = b / a;
      } else if ( a <= 0. ) {
        eta = ( a - Math.sqrt( Math.abs( a * a - 4. * b * c ) ) )
          / ( 2. * c );
      } else {
        eta = 2. * b
          / ( a + Math.sqrt( Math.abs( a * a - 4. * b * c ) ) );
      }

//    document.getElementById("debug_textarea").value +=
//      "! swtch3, c,a,b,eta = " + c + " " + a + " " + b + " " + eta
//      + "\n";

    } else {
      temp = rhoinv + psi + phi;
      if ( orgati ) {
        var temp1 = z[ ioffz + iim1 - 1 ]
          / delta[ ioffdelta + iim1 - 1 ];
        temp1 = temp1 * temp1;
        c = temp - delta[ ioffdelta + iip1 - 1 ] * ( dpsi + dphi )
          - ( d[ ioffd + iim1 - 1 ] - d[ ioffd + iip1 - 1 ] ) * temp1;
        zz[ 0 ] = z[ ioffz + iim1 - 1 ] * z[ ioffz + iim1 - 1 ];
        zz[ 2 ] = delta[ ioffdelta + iip1 - 1 ]
          * delta[ ioffdelta + iip1 - 1 ]
          * ( ( dpsi - temp1 ) + dphi );
      } else {
        temp1 = z[ ioffz + iip1 - 1 ] / delta[ ioffdelta + iip1 - 1 ];
        temp1 = temp1 * temp1;
        c = temp - delta[ ioffdelta + iim1 - 1 ] * ( dpsi + dphi )
          - ( d[ ioffd + iip1 - 1 ] - d[ ioffd + iim1 - 1 ] ) * temp1;
        zz[ 0 ] = delta[ ioffdelta + iim1 - 1 ]
          * delta[ ioffdelta + iim1 - 1 ]
          * ( dpsi + ( dphi - temp1 ) );
        zz[ 2 ] = z[ ioffz + iip1 - 1 ] * z[ ioffz + iip1 - 1 ]
      }
      zz[ 1 ] = z[ ioffz + ii - 1 ] * z[ ioffz + ii - 1 ];
      var etaref = new NumberReference( eta );
      LaPack1.dlaed6( niter, orgati, c, delta, zz, w, etaref, info,
        ioffdelta + iim1 - 1, 0 );
      eta = etaref.getValue();

//    document.getElementById("debug_textarea").value +=
//      "swtch3, temp,temp1,c,eta = " + temp + " " + temp1 + " " + c
//      + " " + eta + "\n";
//    document.getElementById("debug_textarea").value += 
//      "zz = " + zz[ 0 ] + " " + zz[ 1 ] + " " + zz[ 2 ] + "\n";

      if ( info.getValue() == 0 ) return;
    }
    if ( w * eta >= 0. ) eta = - w / dw;
    temp = tau + eta;
    if ( temp > dltub || temp < dltlb ) {
      eta = ( w < 0. ? ( dltub - tau ) / 2. : ( dltlb - tau ) / 2. );
    }
    var prew = w;

//  document.getElementById("debug_textarea").value +=
//    "eta,temp,prew = " + eta + " " + temp + " " + prew + "\n";

    for ( j = 1; j <= n; j ++ ) delta[ ioffdelta + j - 1 ] -= eta; // L180

//  document.getElementById("debug_textarea").value += "L180 : delta = ";
//  for ( var jj = 0; jj < n; jj ++ ) {
//    var deltai = delta[ ioffdelta + jj ];
//    document.getElementById("debug_textarea").value += deltai + " ";
//  }
//  document.getElementById("debug_textarea").value += "\n";

    dpsi = 0.;
    psi = 0.;
    erretm = 0.;
    for ( j = 1; j <= iim1; j ++ ) { // loop 190
      temp = z[ ioffz + j - 1 ] / delta[ ioffdelta + j - 1 ];
      psi += z[ ioffz + j - 1 ] * temp;
      dpsi += temp * temp;
      erretm += psi;
    }
    erretm = Math.abs( erretm );

//  document.getElementById("debug_textarea").value +=
//    "L190 : dpsi,psi,erretm = " + dpsi + " " + psi + " " + erretm
//    + "\n";

    dphi = 0.;
    phi = 0.;
    for ( j = n; j >= iip1; j -- ) { // loop 200
      temp = z[ ioffz + j - 1 ] / delta[ ioffdelta + j - 1 ];
      phi += z[ ioffz + j - 1 ] * temp;
      dphi += temp * temp;
      erretm += phi;
    }

//  document.getElementById("debug_textarea").value +=
//    "L200 : dphi,phi,erretm = " + dphi + " " + phi + " " + erretm
//    + "\n";

    temp = z[ ioffz + ii - 1 ] / delta[ ioffdelta + ii - 1 ];
    dw = dpsi + dphi + temp * temp;
    temp = z[ ioffz + ii - 1 ] * temp;
    w = rhoinv + phi + psi + temp;
    erretm += 8. * ( phi - psi ) + 2. * rhoinv
      + 3. * Math.abs( temp ) + Math.abs( tau + eta ) * dw;
    var swtch = false;
    if ( orgati ) {
      if ( - w > Math.abs( prew ) / 10. ) swtch = true;
    } else {
      if ( w > Math.abs( prew ) / 10. ) swtch = true;
    }
    tau += eta;
    iter = niter + 1;

//  document.getElementById("debug_textarea").value +=
//    "temp,dw,erretm,swtch,tau,iter = " + temp + " " + dw + " "
//    + erretm + " " + swtch + " " + tau + " " + iter + "\n";

    for ( niter = iter; niter <= maxit; niter ++ ) {
//    document.getElementById("debug_textarea").value +=
//      "\n\nniter = " + niter + "\n\n";
      if ( Math.abs( w ) <= eps * erretm ) {
        dlam.setValue( ( orgati ? d[ ioffd + i - 1 ] + tau
          : d[ ioffd + ip1 - 1 ] + tau ) );
//      document.getElementById("debug_textarea").value +=
//        "| w | <= eps * erretm: dlam = " + dlam.getValue() + "\n";
//      document.getElementById("debug_textarea").value += "delta = ";
//      for ( var jj = 0; jj < n; jj ++ ) {
//        var deltai = delta[ ioffdelta + jj ];
//        document.getElementById("debug_textarea").value += deltai + " ";
//      }
//      document.getElementById("debug_textarea").value += "\n";
        return;
      }
      if ( w <= 0. ) dltlb = Math.max( dltlb, tau );
      else dltub = Math.min( dltub, tau );

//    document.getElementById("debug_textarea").value +=
//      "dltlb,dltub = " + dltlb + " " + dltub + "\n";

      if ( ! swtch3 ) {
//      document.getElementById("debug_textarea").value += "! swtch3\n";
        if ( ! swtch ) {
          if ( orgati ) {
            c = w - delta[ ioffdelta + ip1 - 1 ] * dw
              - ( d[ ioffd + i - 1 ] - d[ ioffd + ip1 - 1 ] )
              * Math.pow( z[ ioffz + i - 1 ]
              / delta[ ioffdelta + i - 1 ] , 2 );
          } else {
            c = w - delta[ ioffdelta + i - 1 ] * dw
              - ( d[ ioffd + ip1 - 1 ] - d[ ioffd + i - 1 ] )
              * Math.pow( z[ ioffz + ip1 - 1 ]
              / delta[ ioffdelta + ip1 - 1 ] , 2 );
          }
        } else {
          temp = z[ ioffz + ii - 1 ] / delta[ ioffdelta + ii - 1 ];
          if ( orgati ) dpsi += temp * temp;
          else dphi += temp * temp;
          c = w - delta[ ioffdelta + i - 1 ] * dpsi
            - delta[ ioffdelta + ip1 - 1 ] * dphi;
        }
        a = ( delta[ ioffdelta + i - 1 ]
          + delta[ ioffdelta + ip1 - 1 ] ) * w
          - delta[ ioffdelta + i - 1 ] * delta[ ioffdelta + ip1 - 1 ]
          * dw;
        b = delta[ ioffdelta + i - 1 ] * delta[ ioffdelta + ip1 - 1 ]
          * w;

//      document.getElementById("debug_textarea").value +=
//        "c,a,b = " + c + " " + a + " " + b + "\n";

        if ( c == 0 ) {
          if ( a == 0. ) {
            if ( ! swtch ) {
              a = ( orgati ?  z[ ioffz + i - 1 ] * z[ ioffz + i - 1 ]
                + delta[ ioffdelta + ip1 - 1 ]
                * delta[ ioffdelta + ip1 - 1 ] * ( dpsi + dphi )
                : z[ ioffz + ip1 - 1 ] * z[ ioffz + ip1 - 1 ]
                + delta[ ioffdelta + i - 1 ]
                * delta[ ioffdelta + i - 1 ] * ( dpsi + dphi ) );
            } else {
              a = delta[ ioffdelta + i - 1 ]
                * delta[ ioffdelta + i - 1 ] * dpsi
                + delta[ ioffdelta + ip1 - 1 ]
                * delta[ ioffdelta + ip1 - 1 ] * dphi;
            }
          }
          eta = b / a;
        } else if ( a <= 0. ) {
          eta = ( a - Math.sqrt( Math.abs( a * a - 4. * b * c ) ) )
            / ( 2. * c );
        } else {
          eta = 2. * b
            / ( a + Math.sqrt( Math.abs( a * a - 4. * b * c ) ) );
        }

//      document.getElementById("debug_textarea").value +=
//        "a,eta = " + a + " " + eta + "\n";

      } else {
//      document.getElementById("debug_textarea").value += "swtch3\n";
        temp = rhoinv + psi + phi;
        if ( swtch ) {
          c = temp - delta[ ioffdelta + iim1 - 1 ] * dpsi
            - delta[ ioffdelta + iip1 - 1 ] * dphi;
          zz[ 0 ] = delta[ ioffdelta + iim1 - 1 ]
            * delta[ ioffdelta + iim1 - 1 ] * dpsi;
          zz[ 2 ] = delta[ ioffdelta + iip1 - 1 ]
            * delta[ ioffdelta + iip1 - 1 ] * dphi;
        } else {
          if ( orgati ) {
            temp1 = z[ ioffz + iim1 - 1 ]
              / delta[ ioffdelta + iim1 - 1 ];
            temp1 = temp1 * temp1;
            c = temp - delta[ ioffdelta + iip1 - 1 ] * ( dpsi + dphi )
              - ( d[ ioffd + iim1 - 1 ]
              - d[ ioffd + iip1 - 1 ] ) * temp1;
            zz[ 0 ] = delta[ ioffdelta + iim1 - 1 ]
              * z[ ioffz + iim1 - 1 ];
            zz[ 2 ] = delta[ ioffdelta + iip1 - 1 ]
              * delta[ ioffdelta + iip1 - 1 ]
              * ( ( dpsi - temp1 ) + dphi );
          } else {
            temp1 = z[ ioffz + iip1 - 1 ]
              / delta[ ioffdelta + iip1 - 1 ];
            temp1 = temp1 * temp1;
            c = temp - delta[ ioffdelta + iim1 - 1 ] * ( dpsi + dphi )
              - ( d[ ioffd + iip1 - 1 ]
              - d[ ioffd + iim1 - 1 ] ) * temp1;
            zz[ 0 ] = delta[ ioffdelta + iim1 - 1 ]
              * delta[ ioffdelta + iim1 - 1 ]
              * ( dpsi + ( dphi - temp1 ) );
            zz[ 2 ] = z[ ioffz + iip1 - 1 ] * z[ ioffz + iip1 - 1 ];
          }
        }
        etaref.setValue( eta );
        LaPack1.dlaed6( niter, orgati, c, delta, zz, w, etaref, info,
          ioffdelta + iim1 - 1, 0 );
        eta = etaref.getValue();

//      document.getElementById("debug_textarea").value +=
//        "after dlaed6, temp,c,eta = " + temp + " " + c + " " + eta
//        + "\n";
//      document.getElementById("debug_textarea").value +=
//        "zz = " + zz[ 0 ] + " " + zz[ 1 ] + " " + zz[ 2 ] + "\n";

        if ( info.getValue() != 0 ) return;
      }
      if ( w * eta >= 0. ) eta = - w / dw;
      temp = tau + eta;
      if ( temp > dltub || temp < dltlb ) {
        eta = ( w < 0. ? ( dltub - tau ) / 2. : ( dltlb - tau ) / 2. );
      }

//    document.getElementById("debug_textarea").value +=
//      "eta,temp = " + eta + " " + temp + "\n";

      for ( j = 1; j <= n; j ++ ) delta[ ioffdelta + j - 1 ] -= eta;//L210

//    document.getElementById("debug_textarea").value +=
//      "L210 : delta = ";
//    for ( var jj = 0; jj < n; jj ++ ) {
//      var deltai = delta[ ioffdelta + jj ];
//      document.getElementById("debug_textarea").value += deltai + " ";
//    }
//    document.getElementById("debug_textarea").value += "\n";

      tau += eta;
      prew = w;
      dpsi = 0.;
      psi = 0.;
      erretm = 0.;
      for ( j = 1; j <= iim1; j ++ ) { // loop 220
        temp = z[ ioffz + j - 1 ] / delta[ ioffdelta + j - 1 ];
        psi += z[ ioffz + j - 1 ] * temp;
        dpsi += temp * temp;
        erretm += psi;
      }
      erretm = Math.abs( erretm );

//    document.getElementById("debug_textarea").value +=
//      "L220 : tau,prew,dpsi,psi,erretm = " + tau + " " + prew + " "
//      + dpsi + " " + psi + " " + erretm + "\n";

      dphi = 0.;
      phi = 0.;
      for ( j = n; j >= iip1; j -- ) { // loop 230
        temp = z[ ioffz + j - 1 ] / delta[ ioffdelta + j - 1 ];
        phi += z[ ioffz + j - 1 ] * temp;
        dphi += temp * temp;
        erretm += phi;
      }

//    document.getElementById("debug_textarea").value +=
//      "L230 : dphi,phi,erretm = " + dphi + " " + phi + " " + erretm
//      + "\n";

      temp = z[ ioffz + ii - 1 ] / delta[ ioffdelta + ii - 1 ];
      dw = dpsi + dphi + temp * temp;
      temp *= z[ ioffz + ii - 1 ];
      w = rhoinv + phi + psi + temp;
      erretm += 8. * ( phi - psi ) + 2. * rhoinv
        + 3. * Math.abs( temp ) + Math.abs( tau ) * dw;
      if ( w * prew > 0. && Math.abs( w ) > Math.abs( prew ) / 10. ) {
        swtch = ! swtch;
      }

//    document.getElementById("debug_textarea").value +=
//      "temp,dw,w,erretm,swtch = " + temp + " " + dw + " " + w + " "
//      + erretm + " " + swtch + "\n";

    }
    info.setValue( 1 );
    dlam.setValue( ( orgati ? d[ ioffd + i - 1 ] + tau
      : d[ ioffd + ip1 - 1 ] + tau ) );
  }
//document.getElementById("debug_textarea").value +=
//  "leaving dlaed4, dlam = " + dlam.getValue() + "\n";
//document.getElementById("debug_textarea").value += "delta = ";
//for ( var jj = 0; jj < n; jj ++ ) {
//  var deltai = delta[ ioffdelta + jj ];
//  document.getElementById("debug_textarea").value += deltai + " ";
//}
//document.getElementById("debug_textarea").value += "\n";
}
//************************************************************************
LaPack2.dlaein = function( rightv, noinit, n,
H, ldh, wr, wi, vr, vi, B, ldb, work, eps3, smlnum, bignum, info, ioffh,
ioffvr, ioffvi, ioffb, ioffwork ) {
  info.setValue( 0 );
  var rootn = Math.sqrt( Number( n ) );
  var growto = 0.1 / rootn;
  var nrmsml = Math.max( 1., eps3 * rootn ) * smlnum;
  for ( var j = 1; j <= n; j ++ ) {
    for ( var i = 1; i <= j - 1; i ++ ) {
      B[ ioffb + i - 1 + ( j - 1 ) * ldb ] =
        H[ ioffh + i - 1 + ( j - 1 ) * ldh ];
    }
    B[ ioffb + j - 1 + ( j - 1 ) * ldb ] =
      H[ ioffh + j - 1 + ( j - 1 ) * ldh ] - wr;
  }
  var scale = new NumberReference();
  var ierr = new IntReference();
  if ( wi == 0. ) {
    if ( noinit ) {
      for ( i = 1; i <= n; i ++ ) vr [ ioffvr + i - 1 ] = eps3;
    } else {
      var vnorm = Blas1.dnrm2( n, vr, 1, ioffvr );
      Blas1.dscal( n, ( eps3 * rootn ) / Math.max( vnorm, nrmsml ),
        vr, 1, ioffvr );
    }
    if ( rightv ) {
      for ( i = 1; i <= n - 1; i ++ ) {
        var ei = H[ ioffh + i + ( i - 1 ) * ldh ];
        if ( Math.abs( B[ ioffb + i - 1 + ( i - 1 ) * ldb ] )
        < Math.abs( ei ) ) {
          var x = B[ ioffb + i - 1 + ( i - 1 ) * ldb ] / ei;
          B[ ioffb + i - 1 + ( i - 1 ) * ldb ] = ei;
          for ( j = i + 1; j <= n; j ++ ) {
            var temp = B[ ioffb + i + ( j - 1 ) * ldb ];
            B[ ioffb + i + ( j - 1 ) * ldb ] =
              B[ ioffb + i - 1 + ( j - 1 ) * ldb ] - x * temp;
            B[ ioffb + i - 1 + ( j - 1 ) * ldb ] = temp;
          }
        } else {
          if ( B[ ioffb + i - 1 + ( i - 1 ) * ldb ] == 0. ) {
            B[ ioffb + i - 1 + ( i - 1 ) * ldb ] = eps3;
          }
          x = ei / B[ ioffb + i - 1 + ( i - 1 ) * ldb ];
          if ( x != 0. ) {
            for ( j = i + 1; j <= n; j ++ ) {
              B[ ioffb + i + ( j - 1 ) * ldb ] -=
                x * B[ ioffb + i - 1 + ( j - 1 ) * ldb ];
            }
          }
        }
      }
      if ( B[ ioffb + n - 1 + ( n - 1 ) * ldb ] == 0. ) {
        B[ ioffb + n - 1 + ( n - 1 ) * ldb ] = eps3;
      }
      var trans = 'N';
    } else {
      for ( j = n; j >= 2; j -- ) {
        var ej = H[ ioffh + j - 1 + ( j - 2 ) * ldh ];
        if ( Math.abs( B[ ioffb + j - 1 + ( j - 1 ) * ldb ] )
        < Math.abs( ej ) ) {
          x = B[ ioffb + j - 1 + ( j - 1 ) * ldb ] / ej;
          B[ ioffb + j - 1 + ( j - 1 ) * ldb ] = ej;
          for ( i = 1; i <= j - 1; i ++ ) {
            temp = B[ ioffb + i - 1 + ( j - 2 ) * ldb ];
            B[ ioffb + i - 1 + ( j - 2 ) * ldb ] =
              B[ ioffb + i - 1 + ( j - 1 ) * ldb ] - x * temp;
            B[ ioffb + i - 1 + ( j - 1 ) * ldb ] = temp;
          }
        } else {
          if ( B[ ioffb + j - 1 + ( j - 1 ) * ldb ] == 0. ) {
            B[ ioffb + j - 1 + ( j - 1 ) * ldb ] = eps3;
          }
          x = ej / B[ ioffb + j - 1 + ( j - 1 ) * ldb ];
          if ( x != 0. ) {
            for ( i = 1; i <= j - 1; i ++ ) {
              B[ ioffb + i - 1 + ( j - 2 ) * ldb ] -=
                x * B[ ioffb + i - 1 + ( j - 1 ) * ldb ];
            }
          }
        }
      }
      if ( B[ ioffb ] == 0. ) B[ ioffb ] = eps3;
      trans = 'T';
    }
    var normin = 'N';
    var goto120 = false;
    for ( var its = 1; its <= n; its ++ ) {
      LaPack1.dlatrs( 'Upper', trans, 'Nonunit', normin, n, B, ldb,
        vr, scale, work, ierr, ioffb, ioffvr, ioffwork );
      normin = 'Y';
      vnorm = Blas1.dasum( n, vr, 1, ioffvr );
      if ( vnorm >= growto * scale.getValue() ) {
        goto120 = true;
        break;
      }
      temp = eps3 / ( rootn + 1. );
      vr[ ioffvr ] = eps3;
      for ( i = 2; i <= n; i ++ ) vr[ ioffvr + i - 1 ] = temp;
      vr[ ioffvr + n - its ] -= eps3 * rootn;
    } // 110
    if ( ! goto120 ) info.setValue( 1 );
    i = Blas1.idamax( n, vr, 1, ioffvr );
    Blas1.dscal( n, 1. / Math.abs( vr[ ioffvr + i - 1 ] ), vr, 1,
      ioffvr );
  } else { // wi != 0
    if ( noinit ) {
      for ( i = 1; i <= n; i ++ ) {
        vr[ ioffvr + i - 1 ] = eps3;
        vi[ ioffvi + i - 1 ] = 0.;
      }
    } else {
      var norm = LaPack0.dlapy2(
        Blas1.dnrm2( n, vr, 1, ioffvr ),
        Blas1.dnrm2( n, vi, 1, ioffvi ) );
      var rec = ( eps3 * rootn ) / Math.max( norm, nrmsml );
      Blas1.dscal( n, rec, vr, 1, ioffvr );
      Blas1.dscal( n, rec, vi, 1, ioffvi );
    }
    if ( rightv ) {
      B[ ioffb + 1 ] = - wi;
      for ( i = 2; i <= n; i ++ ) B[ ioffb + i ] = 0.;
      for ( i = 1; i <= n - 1; i ++ ) {
        var absbii =
          LaPack0.dlapy2( B[ ioffb + i - 1 + ( i - 1 ) * ldb ],
            B[ ioffb + i + ( i - 1 ) * ldb ] );
        ei = H[ ioffh + i + ( i - 1 ) * ldh ];
        if ( absbii < Math.abs( ei ) ) {
          var xr = B[ ioffb + i - 1 + ( i - 1 ) * ldb ] / ei;
          var xi = B[ ioffb + i + ( i - 1 ) * ldb ] / ei;
          B[ ioffb + i - 1 + ( i - 1 ) * ldb ] = ei;
          B[ ioffb + i + ( i - 1 ) * ldb ] = 0.;
          for ( j = i + 1; j <= n; j ++ ) {
            temp = B[ ioffb + i + ( j - 1 ) * ldb ];
            B[ ioffb + i + ( j - 1 ) * ldb ] =
              B[ ioffb + i - 1 + ( j - 1 ) * ldb ] - xr * temp;
            B[ ioffb + j + i * ldb ] =
              B[ ioffb + j + ( i - 1 ) * ldb ] - xi * temp;
            B[ ioffb + i - 1 + ( j - 1 ) * ldb ] = temp;
            B[ ioffb + j + ( i - 1 ) * ldb ] = 0.;
          }
          B[ ioffb + i + 1 + ( i - 1 ) * ldb ] = - wi;
          B[ ioffb + i + i * ldb ] -= xi * wi;
          B[ ioffb + i + 1 + i * ldb ] += xr * wi;
        } else {
          if ( absbii == 0. ) {
            B[ ioffb + i - 1 + ( i - 1 ) * ldb ] = eps3;
            B[ ioffb + i + ( i - 1 ) * ldb ] = 0.;
            absbii = eps3;
          }
          ei = ( ei / absbii ) / absbii;
          xr = B[ ioffb + i - 1 + ( i - 1 ) * ldb ] * ei;
          xi = - B[ ioffb + i + ( i - 1 ) * ldb ] * ei;
          for ( j = i + 1; j <= n; j ++ ) {
            B[ ioffb + i + ( j - 1 ) * ldb ] +=
              - xr * B[ ioffb + i - 1 + ( j - 1 ) * ldb ]
              + xi * B[ ioffb + j + ( i - 1 ) * ldb ];
            B[ ioffb + j + i * ldb ] =
              - xr * B[ ioffb + j + ( i - 1 ) * ldb ]
              - xi * B[ ioffb + i - 1 + ( j - 1 ) * ldb ];
          }
          B[ ioffb + i + 1 + i * ldb ] -= wi;
        }
        work[ ioffwork + i - 1 ] =
          Blas1.dasum( n - i, B, ldb, ioffb + i - 1 + i * ldb )
          + Blas1.dasum( n - i, B, 1,
            ioffb + i + 1 + ( i - 1 ) * ldb );
      } // 170
      if ( B[ ioffb + n - 1 + ( n - 1 ) * ldb ] == 0. &&
      B[ ioffb + n + ( n - 1 ) * ldb ] == 0. ) {
        B[ ioffb + n - 1 + ( n - 1 ) * ldb ] = eps3;
      }
      work[ ioffwork + n - 1 ] = 0.;
      var i1 = n;
      var i2 = 1;
      var i3 = -1;
    } else { // not rightv
      B[ ioffb + n + ( n - 1 ) * ldb ] = wi;
      for ( j = 1; j <= n - 1; j ++ ) {
        B[ ioffb + n + ( j - 1 ) * ldb ] = 0.;
      }
      for ( j = n; j >= 2; j -- ) {
        ej = H[ ioffh + j - 1 + ( j - 2 ) * ldh ];
        var absbjj = LaPack0.dlapy2(
          B[ ioffb + j - 1 + ( j - 1 ) * ldb ],
          B[ ioffb + j + ( j - 1 ) * ldb ] );
        if ( absbjj < Math.abs( ej ) ) {
          xr = B[ ioffb + j - 1 + ( j - 1 ) * ldb ] / ej;
          xi = B[ ioffb + j + ( j - 1 ) * ldb ] / ej;
          B[ ioffb + j - 1 + ( j - 1 ) * ldb ] = ej;
          B[ ioffb + j + ( j - 1 ) * ldb ] = 0.;
          for ( i = 1; i <= j - 1; i ++ ) {
            temp = B[ ioffb + i - 1 + ( j - 2 ) * ldb ];
            B[ ioffb + i - 1 + ( j - 2 ) * ldb ] =
              B[ ioffb + i - 1 + ( j - 1 ) * ldb ] - xr * temp;
            B[ ioffb + j - 1 + ( i - 1 ) * ldb ] =
              B[ ioffb + j + ( i - 1 ) * ldb ] - xi * temp;
            B[ ioffb + i - 1 + ( j - 1 ) * ldb ] = temp;
            B[ ioffb + j + ( i - 1 ) * ldb ] = 0.;
          }
          B[ ioffb + j + ( j - 2 ) * ldb ] = wi;
          B[ ioffb + j - 2 + ( j - 2 ) * ldb ] += xi * wi;
          B[ ioffb + j - 1 + ( j - 2 ) * ldb ] -= xr * wi;
        } else {
          if ( absbjj == 0. ) {
            B[ ioffb + j - 1 + ( j - 1 ) * ldb ] = eps3;
            B[ ioffb + j + ( j - 1 ) * ldb ] = 0.;
            absbjj = eps3;
          }
          ej = ( ej / absbjj ) / absbjj;
          xr = B[ ioffb + j - 1 + ( j - 1 ) * ldb ] * ej;
          xi = - B[ ioffb + j + ( j - 1 ) * ldb ] * ej;
          for ( i = 1; i <= j - 1; i ++ ) {
            var bij = B[ ioffb + i - 1 + ( j - 1 ) * ldb ];
            var bjp1i = B[ ioffb + j + ( i - 1 ) * ldb ];
            var bijm1 = B[ ioffb + i - 1 + ( j - 2 ) * ldb ];
            var jm1 = j - 1;
            var jp1 = j + 1;
            B[ ioffb + i - 1 + ( j - 2 ) * ldb ] +=
              - xr * B[ ioffb + i - 1 + ( j - 1 ) * ldb ]
              + xi * B[ ioffb + j + ( i - 1 ) * ldb ];
            bijm1 = B[ ioffb + i - 1 + ( j - 2 ) * ldb ];

            var bji = B[ ioffb + j - 1 + ( i - 1 ) * ldb ];
            bjp1i = B[ ioffb + j + ( i - 1 ) * ldb ];
            bij = B[ ioffb + i - 1 + ( j - 1 ) * ldb ];
            B[ ioffb + j - 1 + ( i - 1 ) * ldb ] =
             - xr * B[ ioffb + j + ( i - 1 ) * ldb ]
             - xi * B[ ioffb + i - 1  + ( j - 1 ) * ldb ];
            bji = B[ ioffb + j - 1 + ( i - 1 ) * ldb ];
          }
          jm1=j-1;
          var bjjm1 = B[ ioffb + j - 1 + ( j - 2 ) * ldb ];
          B[ ioffb + j - 1 + ( j - 2 ) * ldb ] += wi;
          bjjm1 = B[ ioffb + j - 1 + ( j - 2 ) * ldb ];
        }
        work[ ioffwork + j - 1 ] =
          Blas1.dasum( j - 1, B, 1, ioffb + ( j - 1 ) * ldb )
          + Blas1.dasum( j - 1, B, ldb, ioffb + j );
      } // 210
      if ( B[ ioffb ] == 0. && B[ ioffb + 1 ] == 0. ) {
        B[ ioffb ] = eps3;
      }
      work[ ioffwork ] = 0.;
      i1 = 1;
      i2 = n;
      i3 = 1;
    }
    var goto280 = false;
    for ( its = 1; its <= n; its ++ ) {
      scale.setValue( 1. );
      var vmax = 1.;
      var vcrit = bignum;
      for ( i = i1; ( i3 == 1 && i <= i2 ) || ( i3 == -1 && i >= i2 );
      i += i3 ) {
        if ( work[ ioffwork + i - 1 ] > vcrit ) {
          rec = 1. / vmax;
          Blas1.dscal( n, rec, vr, 1, ioffvr );
          Blas1.dscal( n, rec, vi, 1, ioffvi );
          scale.setValue( scale.getValue() * rec );
          vmax = 1.;
          vcrit = bignum;
        }
        xr = vr[ ioffvr + i - 1 ];
        xi = vi[ ioffvi + i - 1 ];
        if ( rightv )  {
          for ( j = i + 1; j <= n; j ++ ) {
            xr += - B[ ioffb + i - 1 + ( j - 1 ) * ldb ]
              * vr[ ioffvr + j - 1 ]
              + B[ ioffb + j + ( i - 1 ) * ldb ]
                * vi[ ioffvi + j - 1 ];
            xi += - B[ ioffb + i - 1 + ( j - 1 ) * ldb ]
              * vi[ ioffvi + j - 1 ]
              - B[ ioffb + j + ( i - 1 ) * ldb ]
                * vr[ ioffvr + j - 1 ];
          }
        } else {
          for ( j = 1; j <= i - 1; j ++ ) {
            xr += - B[ ioffb + j - 1 + ( i - 1 ) * ldb ]
              * vr[ ioffvr + j - 1 ]
              + B[ ioffb + i + ( j - 1 ) * ldb ]
                * vi[ ioffvi + j - 1 ];
            xi += - B[ ioffb + j - 1 + ( i - 1 ) * ldb ]
              * vi[ ioffvi + j - 1 ]
              - B[ ioffb + i + ( j - 1 ) * ldb ]
                * vr[ ioffvr + j - 1 ];
          }
        }
        var w = Math.abs( B[ ioffb + i - 1 + ( i - 1 ) * ldb ] )
          + Math.abs( B[ ioffb + i + ( i - 1 ) * ldb ] );
        if ( w > smlnum ) {
          if ( w < 1. ) {
            var w1 = Math.abs( xr ) + Math.abs( xi );
            if ( w1 > w * bignum ) {
              rec = 1. / w1;
              Blas1.dscal( n, rec, vr, 1, ioffvr );
              Blas1.dscal( n, rec, vi, 1, ioffvi );
              xr = vr[ ioffvr + i - 1 ];
              xi = vi[ ioffvi + i - 1 ];
              scale.setValue( scale.getValue() * rec );
              vmax *= rec;
            }
          }
          var p = new NumberReference();
          var q = new NumberReference();
          var k=ioffb + i - 1 + ( i - 1 ) * ldb;
          var kk=ioffb + i + ( i - 1 ) * ldb;
          LaPack0.dladiv( xr, xi, B[ ioffb + i - 1 + ( i - 1 ) * ldb ],
            B[ ioffb + i + ( i - 1 ) * ldb ], p, q );
          vr[ ioffvr + i - 1 ] = p.getValue();
          vi[ ioffvi + i - 1 ] = q.getValue();
          k=ioffvr + i - 1;
          kk=ioffvi + i - 1;
          vmax = Math.max( Math.abs( vr[ ioffvr + i - 1 ] )
            + Math.abs( vi[ ioffvi + i - 1 ] ), vmax );
          vcrit = bignum / vmax;
        } else {
          for ( j = 1; j <= n; j ++ ) {
            vr[ ioffvr + j - 1 ] = 0.;
            vi[ ioffvi + j - 1 ] = 0.;
          }
          vr[ ioffvr + i - 1 ] = 1.;
          vi[ ioffvi + i - 1 ] = 1.;
          scale.setValue( 0. );
          vmax = 1.;
          vcrit = bignum;
        }
      } // 250
      vnorm = Blas1.dasum( n, vr, 1, ioffvr )
        + Blas1.dasum( n, vi, 1, ioffvi );
      if ( vnorm >= growto * scale.getValue() ) {
        goto280 = true;
        break;
      }
      var y = eps3 / ( rootn + 1. );
      vr[ ioffvr ] = eps3;
      vi[ ioffvi ] = 0.;
      for ( i = 2; i <= n; i ++ ) {
        vr[ ioffvr + i - 1 ] = y;
        vi[ ioffvi + i - 1 ] = 0.;
      }
      vr[ ioffvr + n - its ] -= eps3 * rootn;
    } // 270
    if ( ! goto280 ) info.setValue( 1 );
    vnorm = 0.;
    for ( i = 1; i <= n; i ++ ) {
      vnorm = Math.max( vnorm, Math.abs( vr[ ioffvr + i - 1 ] )
        + Math.abs( vi[ ioffvi + i - 1 ] ) );
    }
    Blas1.dscal( n, 1. / vnorm, vr, 1, ioffvr );
    Blas1.dscal( n, 1. / vnorm, vi, 1, ioffvi );
  }
}
LaPack2.zlaein = function( rightv, noinit, n,
H, ldh, w, v, B, ldb, rwork, eps3, smlnum, info ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dlags2 = function( upper, a1, a2, a3, b1, b2,
b3, csu, snu, csv, snv, csq, snq ) {
  throw new Error("not tested: generalized eigenvalue problem");
  var s1 = new NumberReference();
  var s2 = new NumberReference();
  var snr = new NumberReference();
  var csr = new NumberReference();
  var snl = new NumberReference();
  var csl = new NumberReference();
  var r = new NumberReference();
  if ( upper ) {
    var a = a1 * b3;
    var d = a3 * b1;
    var b = a2 * b1 - a1 * b2;
    LaPack1.dlasv2( a, b, d, s1, s2, snr, csr, snl, csl );
    if ( Math.abs( csl.getValue() ) >= Math.abs( snl.getValue() ) ||
    Math.abs( csr.getValue() ) >= Math.abs( snr.getValue() ) ) {
      var ua11r = csl.getValue() * a1;
      var ua12 = csl.getValue() * a2 + snl.getValue() * a3;
      var vb11r = csr.getValue() * b1;
      var vb12 = csr.getValue() * b2 + snr.getValue() * b3;
      var aua12 = Math.abs( csl.getValue() ) * Math.abs( a2 )
        + Math.abs( snl.getValue() ) * Math.abs( a3 );
      var avb12 = Math.abs( csr.getValue() ) * Math.abs( b2 )
        + Math.abs( snr.getValue() ) * Math.abs( b3 );
      if ( ( Math.abs( ua11r ) + Math.abs( ua12 ) ) != 0. ) {
        if ( aua12 / ( Math.abs( ua11r ) + Math.abs( ua12 ) ) 
        <= avb12 / ( Math.abs( vb11r) + Math.abs( vb12 ) ) ) {
          LaPack1.dlartg( - ua11r, ua12, csq, snq, r );
        } else LaPack1.dlartg( - vb11r, vb12, csq, snq, r );
      } else LaPack1.dlartg( - vb11r, vb12, csq, snq, r );
      csu.setValue( csl.getValue() );
      snu.setValue( - snl.getValue() );
      csv.setValue( csr.getValue() );
      snv.setValue( - snr.getValue() );
    } else {
      var ua21 = - snl.getValue() * a1;
      var ua22 = - snl.getValue() * a2 + csl.getValue() * a3;
      var vb21 = - snr.getValue() * b1;
      var vb22 = - snr.getValue() * b2 + csr.getValue() * b3;
      var aua22 = Math.abs( snl.getValue() ) * Math.abs( a2 )
        + Math.abs( csl.getValue() ) * Math.abs( a3 );
      var avb22 = Math.abs( snr.getValue() ) * Math.abs( b2 )
        + Math.abs( csr.getValue() ) * Math.abs( b3 );
      if ( ( Math.abs( ua21 ) + Math.abs( ua22 ) ) != 0. ) {
        if ( aua22 / ( Math.abs( ua21 ) + Math.abs( ua22 ) )
        <= avb22 / ( Math.abs( vb21 ) + Math.abs( vb22 ) ) ) {
          LaPack1.dlartg( - ua21, ua22, csq, snq, r );
        } else LaPack1.dlartg( - vb21, vb22, csq, snq, r );
      } else LaPack1.dlartg( - vb21, vb22, csq, snq, r );
      csu.setValue( snl.getValue() );
      snu.setValue( csl.getValue() );
      csv.setValue( snr.getValue() );
      snv.setValue( csr.getValue() );
    }
  } else {
    a = a1 * b3;
    d = a3 * b1;
    var c = a2 * b3 - a3 * b2;
    LaPack1.dlasv2( a, c, d, s1, s2, snr, csr, snl, csl );
    if ( Math.abs( csr.getValue() ) >= Math.abs( snr.getValue() ) ||
    Math.abs( csl.getValue() ) >= Math.abs( snl.getValue() ) ) {
      ua21 = - snr.getValue() * a1 + csr.getValue() * a2;
      var ua22r = csr.getValue() * a3;
      vb21 = - snl.getValue() * b1 + csl.getValue() * b2;
      var vb22r = csl.getValue() * b3;
      var aua21 = Math.abs( snr.getValue() ) * Math.abs( a1 )
        + Math.abs( csr.getValue() ) * Math.abs( a2 );
      var avb21 = Math.abs( snl.getValue() ) * Math.abs( b1 )
        + Math.abs( csl.getValue() ) * Math.abs( b2 );
      if ( ( Math.abs( ua21 ) + Math.abs( ua22r ) ) != 0. ) {
        if ( aua21 / ( Math.abs( ua21 ) + Math.abs( ua22r ) )
        <= avb21 / ( Math.abs( vb21 ) + Math.abs( vb22r ) ) ) {
          LaPack1.dlartg( ua22r, ua21, csq, snq, r );
        } else LaPack1.dlartg( vb22r, vb21, csq, snq, r );
      } else LaPack1.dlartg( vb22r, vb21, csq, snq, r );
      csu.setValue( csr.getValue() );
      snu.setValue( - snr.getValue() );
      csv.setValue( csl.getValue() );
      snv.setValue( - snl.getValue() );
    } else {
      var ua11 = csr.getValue() * a1 + snr.getValue() * a2;
      ua12 = snr.getValue() * a3;
      var vb11 = csl.getValue() * b1 + snl.getValue() * b2;
      vb12 = snl.getValue() * b3;
      var aua11 = Math.abs( csr.getValue() ) * Math.abs( a1 )
        + Math.abs( snr.getValue() ) * Math.abs( a2 );
      var avb11 = Math.abs( csl.getValue() ) * Math.abs( b1 )
        + Math.abs( snl.getValue() ) * Math.abs( b2 );
      if ( ( Math.abs( ua11 ) + Math.abs( ua12 ) ) != 0. ) {
        if ( aua11 / ( Math.abs( ua11 ) + Math.abs( ua12 ) ) 
        <= avb11 / ( Math.abs( vb11 ) + Math.abs( vb12 ) ) ) {
          LaPack1.dlartg( ua12, ua11, csq, snq, r );
        } else LaPack1.dlartg( vb12, vb11, csq, snq, r );
      } else LaPack1.dlartg( vb12, vb11, csq, snq, r );
      csu.setValue( snr.getValue() );
      snu.setValue( csr.getValue() );
      shqr
      v.setValue( snl.getValue() );
      snv.setValue( csl.getValue() );
    }
  }
}
LaPack2.zlags2 = function( upper, a1, a2, a3, b1, b2,
b3, csu, snu, csv, snv, csq, snq ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dlagv2 = function( A, lda, B, ldb, alphar, alphai,
beta, csl, snl, csr, snr, ioffa, ioffb, ioffalphar, ioffalphai,
ioffbeta
) {
  throw new Error("not programmed: generalized eigenvalue problem");
}
//************************************************************************
LaPack2.dlahqr = function( wantt, wantz, n, ilo, ihi, H, ldh,
wr, wi, iloz, ihiz, Z, ldz, info, ioffh, ioffwr, ioffwi, ioffz ) {
  throw new Error("not tested: complicated input");
  var itmax = 30;
  var dat1 = 0.75;
  var dat2 = -0.4375;
  var v = new Array( 3 );
  info.setValue( 0 );
  if ( n == 0 ) return;
  if ( ilo == ihi ) {
    wr[ ioffwr + ilo - 1 ] = H[ ioffh + ilo - 1 + ( ilo - 1 ) * ldh ];
    wi[ ioffwi + ilo - 1 ] = 0.;
    return;
  }
  for ( var j = ilo; j <= ihi - 3; j ++ ) {
    H[ ioffh + j + 1 + ( j - 1 ) * ldh ] = 0.;
    H[ ioffh + j + 2 + ( j - 1 ) * ldh ] = 0.;
  }
  if ( ilo <= ihi - 2 ) H[ ioffh + ihi - 1 + ( ihi - 3 ) * ldh ] = 0.;
  var nh = ihi - ilo + 1;
  var nz = ihiz - iloz + 1;
  var safmin =
    new NumberReference( LaPack0.dlamch( 'Safe minimum' ) );
  var safmax =
    new NumberReference( 1. / safmin.getValue() );
  LaPack0.dlabad( safmin, safmax );
  var ulp = LaPack0.dlamch( 'Precision' );
  var smlnum = safmin.getValue() * ( Number( nh ) / ulp );
  if ( wantt ) {
    var i1 = 1;
    var i2 = n;
  }
  var i = ihi;
  var goto150 = false;
  while ( true ) { // 20
    var l = ilo;
    if ( i < ilo ) return;
    for ( var its = 0; its <= itmax; its ++ ) {
      for ( var k = i; k <= l + 1; k -- ) {
        if ( Math.abs( H[ ioffh + k - 1 + ( k - 2 ) * ldh ] )
        <= smlnum ) {
          break;
        }
        var tst =
          Math.abs( H[ ioffh + k - 2 + ( k - 2 ) * ldh ] )
          + Math.abs( H[ ioffh + k - 1 + ( k - 1 ) * ldh ] );
        if ( tst == 0. ) {
          if ( k - 2 >= ilo ) {
            tst += Math.abs( H[ ioffh + k - 2 + ( k - 3 ) * ldh ] );
          }
          if ( k + 1 <= ihi ) {
            tst += Math.abs( H[ ioffh + k + ( k - 1 ) * ldh ] );
          }
        }
        if ( Math.abs( H[ ioffh + k - 1 + ( k - 2 ) * ldh ] )
        <= ulp * tst ) {
          var ab = Math.max(
            Math.abs( H[ ioffh + k - 1 + ( k - 2 ) * ldh ] ),
            Math.abs( H[ ioffh + k - 2 + ( k - 1 ) * ldh ] ) );
          var ba = Math.min(
            Math.abs( H[ ioffh + k - 1 + ( k - 2 ) * ldh ] ),
            Math.abs( H[ ioffh + k - 2 + ( k - 1 ) * ldh ] ) );
          var aa = Math.max(
            Math.abs( H[ ioffh + k - 1 + ( k - 1 ) * ldh ] ),
            Math.abs( H[ ioffh + k - 2 + ( k - 2 ) * ldh ]
                    - H[ ioffh + k - 1 + ( k - 1 ) * ldh ] ) );
          var bb = Math.min(
            Math.abs( H[ ioffh + k - 1 + ( k - 1 ) * ldh ] ),
            Math.abs( H[ ioffh + k - 2 + ( k - 2 ) * ldh ]
                    - H[ ioffh + k - 1 + ( k - 1 ) * ldh ] ) );
          var s = aa + bb;
          if ( ba * ( ab / s ) <=
          Math.max( smlnum, ulp * ( bb * ( aa / s ) ) ) ) {
            break;
          }
        }
      } // 30
      l = k;
      if ( l > ilo ) H[ ioffh + l - 1 + ( l - 2 ) * ldh ] = 0.;
      if ( l >= i - 1 ) {
        goto150 = true;
        break;
      }
      if ( ! wantt ) {
        i1 = l;
        i2 = i;
      }
      if ( its == 10 ) {
        s = Math.abs( H[ ioffh + l + ( l - 1 ) * ldh ] )
          + Math.abs( H[ ioffh + l + 1 + l * ldh ] );
        var h11 =
          dat1 * s + H[ ioffh + l - 1 + ( l - 1 ) * ldh];
        var h12 = dat2 * s;
        var h21 = s;
        var h22 = h11;
      } else if ( its == 20 ) {
        s = Math.abs( H[ ioffh + i + ( i - 1 ) * ldh ] )
          + Math.abs( H[ ioffh + i - 2 + ( i - 3 ) * ldh ] );
        h11 = dat1 * s + H[ ioffh + i - 1 + ( i - 1 ) * ldh];
        h12 = dat2 * s;
        h21 = s;
        h22 = h11;
      } else {
        h11 = H[ ioffh + i - 2 + ( i - 2 ) * ldh ];
        h21 = H[ ioffh + i - 1 + ( i - 2 ) * ldh ];
        h12 = H[ ioffh + i - 2 + ( i - 1 ) * ldh ];
        h22 = H[ ioffh + i - 1 + ( i - 1 ) * ldh ];
      }
      s = Math.abs( h11 ) + Math.abs( h12 )
        + Math.abs( h21 ) + Math.abs( h22 );
      if ( s == 0. ) {
        var rt1r = 0;
        var rt1i = 0;
        var rt2r = 0;
        var rt2i = 0;
      } else {
        h11 /= s;
        h21 /= s;
        h12 /= s;
        h22 /= s;
        var tr = ( h11 + h22 ) / 2.;
        var det = ( h11 - tr ) * ( h22 - tr ) - h12 * h21;
        var rtdisc = Math.sqrt( Math.abs( det ) );
        if ( det >= 0. ) {
          rt1r = tr * s;
          rt2r = rt1r;
          rt1i = rtdisc * s;
          rt2i = - rt1i;
        } else {
          rt1r = tr + rtdisc;
          rt2r = tr - rtdisc;
          if ( Math.abs( rt1r - h22 ) <= Math.abs( rt2r - h22 ) ) {
            rt1r *= s;
            rt2r = rt1r;
          } else {
            rt2r *= s;
            rt1r = rt2r;
          }
          rt1i=0.;
          rt2i=0.;
        }
      }
      for ( var m = i - 2; m >= l; m -- ) {
        var h21s = H[ ioffh + m + ( m - 1 ) * ldh ];
        s = Math.abs( H[ ioffh + m - 1 + ( m - 1 ) * ldh ] - rt2r )
          + Math.abs( rt2i ) + Math.abs( h21s );
        h21s = H[ ioffh + m + ( m - 1 ) * ldh ] / s;
        v[ 0 ] = h21s * H[ ioffh + m - 1 + m * ldh ]
          + ( H[ ioffh + m - 1 + ( m - 1 ) * ldh ] - rt1r )
          * ( ( H[ ioffh + m - 1 + ( m - 1 ) * ldh ] - rt2r ) / s )
          - rt1i * ( rt2i / s );
        v[ 1 ] = h21s * ( H[ ioffh + m - 1 + ( m - 1 ) * ldh ]
          + H[ ioffh + m + m * ldh ] - rt1r - rt2r );
        v[ 2 ] = h21s * H[ ioffh + m + 1 + m * ldh ];
        s = Math.abs( v[ 0 ] ) + Math.abs( v[ 1 ] )
          + Math.abs( v[ 2 ] );
        v[ 0 ] /= s;
        v[ 1 ] /= s;
        v[ 2 ] /= s;
        if ( m == l ) break;
        if ( Math.abs( H[ ioffh + m - 1 + ( m - 2 ) * ldh ] )
        * ( Math.abs( v[ 1 ] ) + Math.abs( v[ 2 ] ) ) <=
        ulp * Math.abs( v[ 0 ] )
        * ( Math.abs( H[ ioffh + m - 2 + ( m - 2 ) * ldh ] )
        + Math.abs( H[ ioffh + m - 1 + ( m - 1 ) * ldh ] )
        + Math.abs( H[ ioffh + m + m * ldh ] ) ) ) break;
      } // 50
      for ( k == m; k <= i - 1; k ++ ) {
        var nr = Math.min( 3, i - k + 1 );
        if ( k > m ) {
          Blas1.dcopy( nr, H, 1, v, 1,
            ioffh + k - 1 + ( k - 2 ) * ldh, 0 );
        }
        var alpha = new NumberReference( v[ 0 ] );
        var t1 = new NumberReference();
        LaPack1.dlarfg( nr, alpha, v, 1, t1, 1 );
        v[ 0 ] = alpha.getValue();
        if ( k > m ) {
          H[ ioffh + k - 1 + ( k - 2 ) * ldh ] = v[ 0 ];
          H[ ioffh + k + ( k - 2 ) * ldh ] = 0.;
          if ( k < i - 1 ) H[ ioffh + k + 1 + ( k - 2 ) * ldh ] = 0.;
        } else if ( m > l ) {
          H[ ioffh + k - 1 + ( k - 2 ) * ldh ] *= 1. - t1.getValue();
        }
        var v2 = v[ 1 ];
        var t2 = t1.getValue() * v2;
        if ( nr == 3 ) {
          var v3 = v[ 2 ];
          var t3 = t1.getValue() * v3;
          for ( j = k; j <= i2; j ++ ) {
            var sum = H[ ioffh + k - 1 + ( j - 1 ) * ldh ]
              + v2 * H[ ioffh + k + ( j - 1 ) * ldh ]
              + v3 * H[ ioffh + k + 1 + ( j - 1 ) * ldh ];
            H[ ioffh + k - 1 + ( j - 1 ) * ldh ] -= sum * t1.getValue();
            H[ ioffh + k + ( j - 1 ) * ldh ] -= sum * t2;
            H[ ioffh + k + 1 + ( j - 1 ) * ldh ] -= sum * t3;
          }
          for ( j = i1; j <= Math.min( k + 3, i ); j ++ ) {
            sum = H[ ioffh + j - 1 + ( k - 1 ) * ldh ]
              + v2 * H[ ioffh + j - 1 + k * ldh ]
              + v3 * H[ ioffh + j - 1 + ( k + 1 ) * ldh ];
            H[ ioffh + j - 1 + ( k - 1 ) * ldh ] -= sum * t1.getValue();
            H[ ioffh + j - 1 + k * ldh ] -= sum * t2;
            H[ ioffh + j - 1 + ( k + 1 ) * ldh ] -= sum * t3;
          }
          if ( wantz ) {
            for ( j = iloz; j <= ihiz; j ++ ) {
              sum = Z[ ioffz + j - 1 + ( k - 1 ) * ldz ]
                + v2 * Z[ ioffz + j - 1 + k * ldz ]
                + v3 * Z[ ioffz + j - 1 + ( k + 1 ) * ldz ]
              Z[ ioffz + j - 1 + ( k - 1 ) * ldz ] -= sum * t1.getValue();
              Z[ ioffz + j - 1 + k * ldz ] -= sum * t2;
              Z[ ioffz + j - 1 + ( k + 1 ) * ldz ] -= sum * t3;
            }
          }
        } else if ( nr == 2 ) {
          for ( j = k; j <= i2; j ++ ) {
            sum = H[ ioffh + k - 1 + ( j - 1 ) * ldh ]
              + v2 * H[ ioffh + k + ( j - 1 ) * ldh ];
            H[ ioffh + k - 1 + ( j - 1 ) * ldh ] -= sum * t1.getValue();
            H[ ioffh + k + ( j - 1 ) * ldh ] -= sum * t2;
          }
          for ( j = i1; j <= i; j ++ ) {
            sum = H[ ioffh + j - 1 + ( k - 1 ) * ldh ]
              + v2 * H[ ioffh + j - 1 + k * ldh ];
            H[ ioffh + j - 1 + ( k - 1 ) * ldh ] -= sum * t1.getValue();
            H[ ioffh + j - 1 + k * ldh ] -= sum * t2;
          }
          if ( wantz ) {
            for ( j = iloz; j <= ihiz; j ++ ) {
              sum = Z[ ioffz + j - 1 + ( k - 1 ) * ldz ]
                + v2 * Z[ ioffz + j - 1 + k * ldz ];
              Z[ ioffz + j - 1 + ( k - 1 ) * ldz ] -= sum * t1.getValue();
              Z[ ioffz + j - 1 + k * ldz ] -= sum * t2;
            }
          }
        }
      } // 130
    } // 140
    if ( ! goto150 ) {
      info.setValue( i );
      return;
    }
    if ( l == i ) {
      wr[ ioffwr + i - 1 ] = H[ ioffh + i - 1 + ( i - 1 ) * ldh ];
      wi[ ioffwi + i - 1 ] = 0.;
    } else if ( l == i - 1 ) {
      var a =
        new NumberReference( H[ ioffh + i - 2 + ( i - 2 ) * ldh ] );
      var b =
        new NumberReference( H[ ioffh + i - 2 + ( i - 1 ) * ldh ] );
      var c =
        new NumberReference( H[ ioffh + i - 1 + ( i - 2 ) * ldh ] );
      var d =
        new NumberReference( H[ ioffh + i - 1 + ( i - 1 ) * ldh ] );
      var rt1rref = new NumberReference();
      var rt1iref = new NumberReference();
      var rt2rref = new NumberReference();
      var rt2iref = new NumberReference();
      var cs = new NumberReference();
      var sn = new NumberReference();
      LaPack1.dlanv2( a, b, c, d, rt1rref, rt1iref, rt2rref, rt2iref,
        cs, sn );
      H[ ioffh + i - 2 + ( i - 2 ) * ldh ] = a.getValue();
      H[ ioffh + i - 2 + ( i - 1 ) * ldh ] = b.getValue();
      H[ ioffh + i - 1 + ( i - 2 ) * ldh ] = c.getValue();
      H[ ioffh + i - 1 + ( i - 1 ) * ldh ] = d.getValue();
      wr[ ioffwr + i - 2 ] = rt1rref.getValue();
      wi[ ioffwi + i - 2 ] = rt1iref.getValue();
      wr[ ioffwr + i - 1 ] = rt2rref.getValue();
      wi[ ioffwi + i - 1 ] = rt2iref.getValue();
      if ( wantt ) {
        if ( i2 > i ) {
          Blas1.drot( i2 - i, H, ldh, H, ldh, cs.getValue(), sn.getValue(),
            ioffh + i - 2 + i * ldh, ioffh + i - 1 + i * ldh );
        }
        Blas1.drot( i - i1 - 1, H, 1, H, 1, cs.getValue(), sn.getValue(),
          ioffh + i1 - 1 + ( i - 2 ) * ldh,
          ioffh + i1 - 1 + ( i - 1 ) * ldh );
      }
      if ( wantz ) {
        Blas1.drot( nz, Z, 1, Z, 1, cs.getValue(), sn.getValue(),
          ioffz + iloz - 1 + ( i - 2 ) * ldz,
          ioffz + iloz - 1 + ( i - 1 ) * ldz );
      }
    }
    i = l - 1;
  }
}
LaPack2.zlahqr = function( wantt, wantz, n, ilo, ihi, H, ldh,
w, iloz, ihiz, Z, ldz, info ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dlahr2 = function( n, k, nb, A, lda, tau, T, ldt, Y,
ldy, ioffa, iofftau, iofft, ioffy ) {
  throw new Error("not tested");
  if ( n <= 1 ) return;
  for ( var i = 1; i <= nb; i ++ ) {
    if ( i > 1 ) {
      Blas2.dgemv( 'No transpose', n - k, i - 1, -1., Y, ldy, A, lda,
        1., A, 1, ioffy + k, ioffa + k + i - 2,
        ioffa + k + ( i - 1 ) * lda );
      Blas1.dcopy( i - 1, A, 1, T, 1, ioffa + k + ( i - 1 ) * lda,
        iofft + ( nb - 1 ) * ldt );
      Blas2.dtrmv( 'Lower', 'Transpose', 'Unit', i - 1, A, lda, T, 1,
        ioffa + k, iofft + ( nb - 1 ) * ldt );
      Blas2.dgemv( 'Transpose', n - k - i + 1, i - 1, 1., A, lda, A, 1,
        1., T, 1, ioffa + k + i - 1,
        ioffa + k + i - 1 + ( i - 1 ) * lda,
        iofft + ( nb - 1 ) * ldt );
      Blas2.dtrmv( 'Upper', 'Transpose', 'Non-unit', i - 1, T, ldt,
        T, 1, iofft, iofft + ( nb - 1 ) * ldt );
      Blas2.dgemv( 'No transpose', n - k - i + 1, i - 1, -1., A, lda,
        T, 1, 1., A, 1, ioffa + k + i - 1, iofft + ( nb - 1 ) * ldt,
        ioffa + k + i - 1 + ( i - 1 ) * lda );
      Blas2.dtrmv( 'Lower', 'No transpose', 'Unit', i - 1, A, lda,
        T, 1, ioffa + k, iofft + ( nb - 1 ) * ldt );
      Blas1.daxpy( i - 1, -1., T, 1, A, 1, iofft + ( nb - 1 ) * ldt,
        ioffa + k + ( i - 1 ) * lda );
      A[ ioffa + k + i - 2 + ( i - 2 ) * lda ] = ei;
    }
    var alpha = new
      NumberReference( A[ ioffa + k + i - 1 + ( i - 1 ) * lda ] );
    var tauref =
      new NumberReference( tau[ iofftau + i - 1 ] );
    LaPack1.dlarfg( n - k - i + 1, alpha, A, 1, tauref,
      ioffa + Math.min( k + i + 1, n ) - 1 + ( i - 1 ) * lda );
    A[ ioffa + k + i - 1 + ( i - 1 ) * lda ] = alpha.getValue();
    tau[ iofftau + i - 1 ] = tauref.getValue();
    var ei = A[ ioffa + k + i - 1 + ( i - 1 ) * lda ];
    A[ ioffa + k + i - 1 + ( i - 1 ) * lda ] = 1.;
    Blas2.dgemv( 'No transpose', n - k, n - k - i + 1, 1., A, lda,
      A, 1, 0., Y, 1, ioffa + k + i * lda,
      ioffa + k + i - 1 + ( i - 1 ) * lda,
      ioffy + k + ( i - 1 ) * ldy );
    Blas2.dgemv( 'Transpose', n - k - i + 1, i - 1, 1., A, lda,
      A, 1, 0., T, 1, ioffa + k + i - 1,
      ioffa + k + i - 1 + ( i - 1 ) * lda,
      iofft + ( i - 1 ) * lda );
    Blas2.dgemv( 'No transpose', n - k, i - 1, -1., Y, ldy, T, 1, 1.,
      Y, 1, ioffy + k, iofft + ( i- 1 ) * ldt,
      ioffy + k + ( i - 1 ) * ldy );
    Blas1.dscal( n - k, tau[ iofftau + i - 1 ], Y, 1,
      ioffy + k + ( i - 1 ) * ldy );
    Blas1.dscal( i - 1, - tau[ iofftau + i -1 ], T, 1,
      iofft + ( i - 1 ) * ldt );
    Blas2.dtrmv( 'Upper', 'No transpose', 'Non-unit', i - 1, T, ldt,
      T, 1, iofft, iofft + ( i - 1 ) * ldt );
    T[ iofft + i - 1 + ( i - 1 ) * ldt ] = tau[ iofftau + i - 1 ];
  }
  A[ ioffa + k + nb - 1 + ( nb - 1 ) * lda ] = ei;
  LaPack0.dlacpy( 'All', k, nb, A, lda, Y, ldy, ioffa + lda, ioffy );
  Blas3.dtrmm( 'Right', 'Lower', 'No transpose', 'Unit', k, nb, 1.,
    A, lda, Y, ldy, ioffa + k, ioffy );
  if ( n > k + nb ) {
    Blas3.dgemm( 'No transpose', 'No transpose', k, nb, n - k - nb,
      1., A, lda, A, lda, 1., Y, ldy, ioffa + ( nb + 1 ) * lda,
      ioffa + k + nb, ioffy );
  }
  Blas3.dtrmm( 'Right', 'Upper', 'No transpose', 'Non-unit', k, nb,
    1., T, ldt, Y, ldy, iofft, ioffy );
}
//************************************************************************
//  will be deprecated in a future release
LaPack2.dlahrd = function( n, k, nb, A, lda, tau, T, ldt, Y,
ldy, ioffa, iofftau, iofft, ioffy ) {
  throw new Error("not tested: no longer used");
  if ( n <= 1 ) return;
  var ei = Number.POSITIVE_INFINITY;
  for ( var i = 1; i <= nb; i ++ ) {
    if ( i > 1 ) {
      Blas2.dgemv( 'No transpose', n, i - 1, -1., Y, ldy, A, lda, 1.,
        A, 1, ioffy, ioffa + k + i - 2, ioffa + ( i - 1 ) * lda ); 
      Blas1.dcopy( i - 1, A, 1, T, 1, ioffa + k + ( i - 1 ) * lda,
        iofft + ( nb - 1 ) * ldt );
      Blas2.dtrmv( 'Lower', 'Transpose', 'Unit', i - 1, A, lda, T, 1,
        ioffa + k, iofft + ( nb - 1 ) * ldt );
      Blas2.dgemv( 'Transpose', n - k - i + 1, i - 1, 1., A, lda, A, 1,
        1., T, 1, ioffa + k + i - 1,
        ioffa + k + i - 1 + ( i - 1 ) * lda,
        iofft + ( nb - 1 ) * ldt );
      Blas2.dtrmv( 'Upper', 'Transpose', 'Non-unit', i - 1, T, ldt,
        T, 1, iofft, iofft + ( nb - 1 ) * ldt );
      Blas2.dgemv( 'No transpose', n - k - i + 1, i - 1, 1., A, lda,
        T, 1, 1., A, 1, ioffa + k + i - 1, iofft + ( nb - 1 ) * ldt,
        ioffa + k + i - 1 + ( i - 1 ) * lda );
      Blas2.dtrmv( 'Lower', 'No transpose', 'Unit', i - 1, A, lda,
        T, 1, ioffa + k, iofft + ( nb - 1 ) * ldt );
      Blas1.daxpy( i - 1, -1, T, 1, A, 1, iofft + ( nb - 1 ) * ldt,
        ioffa + k + ( i - 1 ) * lda );
      A[ ioffa + k + i - 2 + ( i - 2 ) * lda ] = ei;
    }
    var alpha =
      new NumberReference( A[ ioffa + k + i - 1 + ( i - 1 ) * lda ] );
    var taui =
      new NumberReference( tau[ iofftau + i - 1 ] );
    LaPack1.dlarfg( n - k - i + 1, alpha, A, 1, taui,
      ioffa + Math.min( k + k + 1 , n ) - 1 + ( i - 1 ) * lda );
    A[ ioffa + k + i - 1 + ( i - 1 ) * lda ] = alpha.getValue();
    tau[ iofftau + i - 1 ] = taui.getValue();
    ei = A[ ioffa + k + i - 1 + ( i - 1 ) * lda ];
    A[ ioffa + k + i - 1 + ( i - 1 ) * lda ] = 1.;
    Blas2.dgemv( 'No transpose', n, n - k - i + 1, 1., A, lda, A, 1,
      0., Y, 1, ioffa + i * lda, ioffa + k + i - 1 + ( i - 1 ) * lda,
      ioffy + ( i - 1 ) * ldy ); 
    Blas2.dgemv( 'Transpose', n - k - i + 1, i - 1, 1., A, lda, A, 1,
      0., T, 1, ioffa + k + i - 1, ioffa + k + i - 1 + ( i - 1 ) * lda,
      iofft + ( i - 1 ) * ldt ); 
    Blas2.dgemv( 'No transpose', n, i - 1, -1., Y, ldy, T, 1,
      1., Y, 1, ioffy, iofft + ( i - 1 ) * ldt,
      ioffy + ( i - 1 ) * ldy ); 
    Blas1.dscal( n, tau[ iofftau + i - 1 ], Y, 1,
      ioffy + ( i - 1 ) * ldy );
    Blas1.dscal( i - 1, - tau[ iofftau + i - 1 ], T, 1,
      iofft + ( i - 1 ) * ldt );
    Blas2.dtrmv( 'Upper', 'No transpose', 'Non-unit', i - 1, T, ldt,
      T, 1, iofft, iofft + ( i - 1 ) * ldt );
    T[ iofft + i - 1 + ( i - 1 ) * ldt ] = tau[ iofftau + i - 1 ];
  }
  A[ ioffa + k + nb - 1 + ( nb - 1 ) * lda ] = ei;
}
LaPack2.zlahrd = function( n, k, nb, A, lda, tau, T, ldt, Y,
ldy ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dlals0 = function( icompq, nl, nr, sqre, nrhs, B, ldb,
Bx, ldbx, perm, givptr, givcol, ldgcol, givnum, ldgnum, poles, difl,
difr, Z, k, c, s, work, info, ioffb, ioffbx, ioffperm, ioffgivnum,
ioffpoles, ioffdifl, ioffdifr, ioffz, ioffwork ) {
  throw new Error("not tested");
  info.setValue( 0 );
  if ( icompq < 0 || icompq > 1 ) info.setValue( -1 );
  else if ( nl < 1 ) info.setValue( -2 );
  else if ( nr < 1 ) info.setValue( -3 );
  else if ( sqre < 0 || sqre > 1 ) info.setValue( -4 );
  var n = nl + nr + 1;
  if ( nrhs < 1 ) info.setValue( -5 );
  else if ( ldb < n ) info.setValue( -7 );
  else if ( ldbx < n ) info.setValue( -9 );
  else if ( givptr < 0 ) info.setValue( -11 );
  else if ( ldgcol < n ) info.setValue( -13 );
  else if ( ldgnum < n ) info.setValue( -15 );
  else if ( k < 1 ) info.setValue( -20 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlals0', -info.getValue() );
    return;
  }
  var m = n + sqre;
  var nlp1 = nl + 1;
  if ( icompq == 0 ) {
    for ( var i = 1; i <= givptr; i ++ ) {
      Blas1.drot( nrhs, B, ldb, B, ldb,
        givnum[ ioffgivnum + i - 1 + ldgnum ],
        givnum[ ioffgivnum + i - 1 ],
        ioffb + givcol[ ioffgivcol + i - 1 + ldgcol ] - 1,
        ioffb + givcol[ ioffgivcol + i - 1 ] - 1,
        ioffgivnum + i - 1 + ldgnum, ioffgivnum + i - 1 );
    }
    Blas1.dcopy( nrhs, B, ldb, Bx, ldbx, ioffb + nlp1 - 1,
      ioffbx );
    for ( i = 2; i <= n; i ++ ) {
      Blas1.dcopy( nrhs, B, ldb, Bx, ldbx,
        ioffb + perm[ ioffperm + i - 1 ] - 1, ioffbx + i - 1 );
    }
    if ( k == 1 ) {
      Blas1.dcopy( nrhs, Bx, ldbx, B, ldb, ioffbx, ioffb );
      if ( z[ ioffz ] < 0. ) {
        Blas1.dscal( nrhs, -1., B, ldb, ioffb );
      }
    } else {
      for ( j = 1; j <= k; j ++ ) {
        var diflj = difl[ ioffdifl + j - 1 ];
        var dj = poles[ ioffpoles + j - 1 ];
        var dsigj = - poles[ ioffpoles + j - 1 + ldgnum ];
        if ( j < k ) {
          difrj = - difr[ ioffdifr + j - 1 ];
          var dsigjp = - poles[ ioffpoles + j + ldgnum ];
        }
        if ( z[ ioffz + j - 1 ] == 0. ||
        poles[ ioffpoles + j - 1 + ldgnum ] == 0. ) {
          work[ ioffwork + j - 1 ] = 0.;
        } else {
          work[ ioffwork + j - 1 ] =
            - poles[ ioffpoles + j - 1 + ldgnum ]
            * z[ ioffz + j - 1 ] / diflj
            / ( poles[ ioffpoles + j - 1 + ldgnum ] + dj );
        }
        for ( i = 1; i <= j - 1; i ++ ) {
          if ( z[ ioffz + i - 1 ] == 0. ||
          poles[ ioffpoles + i - 1 + ldgnum ] == 0. ) {
            work[ ioffwork + i - 1 ] = 0.;
          } else {
            work[ ioffwork + i - 1 ] =
              poles[ ioffpoles + i - 1 + ldgnum ]
              * z[ ioffz + i - 1 ]
              / ( LaPack0.dlamc3( poles[ ioffpoles + i - 1 + ldgnum ],
              dsigj ) - diflj )
              / ( poles[ ioffpoles + i - 1 + ldgnum ] + dj );
          }
        }
        for ( i = j + 1; i <= k; i ++ ) {
          if ( z[ ioffz + i - 1 ] == 0. ||
          poles[ ioffpoles + i - 1 + ldgnum ] == 0. ) {
            work[ ioffwork + i - 1 ] = 0.;
          } else {
            work[ ioffwork + i - 1 ] =
              poles[ ioffpoles + i - 1 + ldgnum ]
              * z[ ioffz + i - 1 ]
              / ( LaPack0.dlamc3( poles[ ioffpoles + i - 1 + ldgnum ],
              dsigjp ) + difrj )
              / ( poles[ ioffpoles + i - 1 + ldgnum ] + dj );
          }
        }
        work[ ioffwork ] = -1.;
        var temp = Blas1.dnrm2( k, work, 1, ioffwork );
        Blas3.dgemv( 'T', k, nrhs, 1., Bx, ldbx, work, 1, 0., B, ldb,
          ioffbx, ioffwork, ioffb + j - 1 );
        LaPack1.dlascl( 'G', 0, 0, temp, 1., 1, nrhs, B, ldb, info,
          ioffb + j - 1 );
      } // 50
    }
    if ( k < Math.max( m, n ) ) {
      LaPack0.dlacpy( 'A', n - k, nrhs, Bx, ldbx, B, ldb,
        ioffbx + k, ioffb + k );
    }
  } else {
    if ( k == 1 ) {
      Blas1.dcopy( nrhs, B, ldb, Bx, ldbx, ioffb, ioffbx );
    } else {
      for ( j = 1; j <= k; j ++ ) {
        dsigj = poles[ ioffpoles + j - 1 + ldgnum ];
        if ( z[ ioffz + j - 1 ] == 0. ) {
          work[ ioffwork + j - 1 ] = 0.;
        } else {
          work[ ioffwork + j - 1 ] =
            - z[ ioffz + j - 1 ] / difl[ ioffdifl + j - 1 ]
            / ( dsigj + poles[ ioffpoles + j - 1 ] )
            / difr[ ioffdifr + j - 1 + ldgnum ];
        }
        for ( i = 1; i <= j - 1; i ++ ) {
          if ( z[ ioffz + j - 1 ] == 0. ) {
            work[ ioffwork + i - 1 ] = 0.;
          } else {
            work[ ioffwork + i - 1 ] = z[ ioffz + j - 1 ]
              / ( LaPack0.dlamc3( dsigj,
              - poles[ ioffpoles + i + ldgnum ] )
              - difr[ ioffdifr + i - 1 ] )
              / ( dsigj + poles[ ioffpoles + i - 1 ] )
              / difr[ ioffdifr + i - 1 + ldgnum ];
          }
        }
        for ( i = j + 1; i <= k; i ++ ) {
          if ( z[ ioffz + j - 1 ] == 0. ) {
            work[ ioffwork + i - 1 ] = 0.;
          } else {
            work[ ioffwork + i - 1 ] = z[ ioffz + j - 1 ]
              / ( LaPack0.dlamc3( dsigj,
              - poles[ ioffpoles + i - 1 + ldgnum ] )
              - difl[ ioffdifl + i - 1 ] )
              / ( dsigj + poles[ ioffpoles + i - 1 ] )
              / difr[ ioffdifr + i - 1 + ldgnum ];
          }
        }
        Blas3.dgemv( 'T', k, nrhs, 1., B, ldb, work, 1, 0., Bx, ldbx,
          ioffb, ioffwork, ioffbx + j - 1 );
      } // 80
      if ( sqre == 1 ) {
        Blas1.dcopy( nrhs, B, ldb, Bx, ldbx, ioffb + m - 1,
          ioffbx + m - 1 );
        Blas1.drot( nrhs, Bx, ldbx, Bx, ldbx, c, s, ioffbx,
          ioffbx + m - 1 );
      }
      if ( k < Math.max( m, n ) ) {
        LaPack0.dlacpy( 'A', n - k, nrhs, B, ldb, Bx, ldbx,
          ioffb + k, ioffbx + k );
      }
      Blas1.dcopy( nrhs, Bx, ldbx, B, ldb, ioffbx, ioffb + nlp1 - 1 );
      if ( sqre == 1 ) {
        Blas1.dcopy( nrhs, Bx, ldbx, B, ldb, ioffbx + m - 1,
          ioffb + m - 1 );
      }
      for ( i = 2; i <= n; i ++ ) {
        Blas1.dcopy( nrhs, Bx, ldbx, B, ldb, ioffbx + i - 1,
          ioffb + perm[ ioffperm + i - 1 ] - 1 );
      }
      for ( i = givptr; i >= 1; i -- ) {
        Blas1.drot( nrhs, B, ldb, B, ldb,
          givnum[ ioffgivnum + i - 1 + ldgnum ],
          - givnum[ ioffgivnum + i - 1 ],
          ioffb + givcol[ ioffgivcol + i - 1 + ldgcol ] - 1,
          ioffb + givcol[ ioffgivcol + i - 1 ] - 1 );
      }
    }
  }
}
//************************************************************************
LaPack2.dlapll = function( n, x, incx, y, incy, ssmin, ioffx,
ioffy ) {
  if ( n <= 1 ) {
    ssmin.setValue( 0. );
    return;
  }
  var alpha = new NumberReference( x[ ioffx ] );
  var tau = new NumberReference();
  LaPack1.dlarfg( n, alpha, x, incx, tau, ioffx + incx );
  x[ ioffx ] = alpha.getValue();
  var a11 = x[ ioffx ];
  x[ ioffx ] = 1.;
  var c =
    - tau.getValue() * Blas1.ddot( n, x, incx, y, incy, ioffx, ioffy );
  Blas1.daxpy( n, c, x, incx, y, incy, ioffx, ioffy );
  alpha.setValue( y[ ioffy + incy ] );
  LaPack1.dlarfg( n - 1, alpha, y, incy, tau, ioffy + 2 * incy );
  y[ ioffy + incy ] = alpha.getValue();
  var a12 = y[ ioffy ];
  var a22 = y[ ioffy + incy ];
  var ssmax = new NumberReference();
  LaPack0.dlas2( a11, a12, a22, ssmin, ssmax );
}
LaPack2.zlapll = function( n, x, incx, y, incy, ssmin ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dlaqp2 = function( m, n, offset, A, lda, jpvt, tau,
vn1, vn2, work, ioffa, ioffjpvt, iofftau, ioffvn1, ioffvn2, ioffwork ) {
  throw new Error("not tested");
  var mn = Math.min( m - offset, n );
  var tol3z = Math.sqrt( LaPack0.dlamch( 'Epsilon' ) );
  for ( var i = 1; i <= mn; i ++ ) {
    var offpi = offset + i;
    var pvt = ( i - 1 )
      + Blas1.idamax( n - i + 1, vn1, 1, ioffvn1 + i - 1 );
    if ( pvt != i ) {
      Blas1.dswap( m, A, 1, A, 1, ioffa + ( pvt - 1 ) * lda,
        ioffa + ( i - 1 ) * lda );
      var itemp = jpvt[ ioffjpvt + pvt - 1 ];
      jpvt[ ioffjpvt + pvt - 1 ] = jpvt[ ioffjpvt + i - 1 ];
      jpvt[ ioffjpvt + i - 1 ] = itemp;
      vn1[ ioffvn1 + pvt - 1 ] = vn1[ ioffvn1 + i - 1 ];
      vn2[ ioffvn2 + pvt - 1 ] = vn2[ ioffvn2 + i - 1 ];
    }
    if ( offpi < m ) {
      var alpha = new
        NumberReference( A[ ioffa + offpi - 1 + ( i - 1 ) * lda ] );
      var tauref =
        new NumberReference( tau[ iofftau + i - 1 ] );
      LaPack1.dlarfg( m - offpi + 1, alpha, A, 1, tauref,
        ioffa + offpi + ( i - 1 ) * lda );
      A[ ioffa + offpi - 1 + ( i - 1 ) * lda ] = alpha.getValue();
      tau[ iofftau + i - 1 ] = tauref.getValue();
    } else {
      alpha.setValue( A[ ioffa + m - 1 + ( i - 1 ) * lda ] );
      tauref.setValue( tau[ iofftau + i - 1 ] );
      LaPack1.dlarfg( 1, alpha, A, 1, tauref,
        ioffa + m - 1 + ( i - 1 ) * lda );
      A[ ioffa + m - 1 + ( i - 1 ) * lda ] = alpha.getValue();
      tau[ iofftau + i - 1 ] = tauref.getValue();
    }
    if ( i <= n ) {
      var aii = A[ ioffa + offpi - 1 + ( i - 1 ) * lda ];
      A[ ioffa + offpi - 1 + ( i - 1 ) * lda ] = 1.;
      LaPack1.dlarf( 'Left', m - offpi + 1, n - i, A, 1,
        tau[ iofftau + i - 1 ], A, lda, work,
        ioffa + offpi - 1 + ( i - 1 ) * lda,
        ioffa + offpi - 1 + i * lda, ioffwork );
      A[ ioffa + offpi - 1 + ( i - 1 ) * lda ] = aii;
    }
    for ( var j = i + 1; j <= n; j ++ ) {
      if ( vn1[ ioffvn1 + j - 1 ] != 0. ) {
        var temp = 1.  - Math.pow(
          Math.abs( A[ ioffa + offpi - 1 + ( j - 1 ) * lda ] )
          / vn1[ ioffvn1 + j - 1 ], 2 );
        temp = Math.max( temp, 0. );
        var temp2 = temp * Math.pow(
          vn1[ ioffvn1 + j - 1 ] / vn2[ ioffvn2 + j - 1 ], 2 );
        if ( temp2 <= tol3z ) {
          if ( offpi < m ) {
            vn1[ ioffvn1 + j - 1 ] = Blas1.dnrm2( m - offpi, A, 1,
              ioffa + offpi + ( j - 1 ) * lda );
            vn2[ ioffvn2 + j - 1 ] = vn1[ ioffvn1 + j - 1 ];
          } else {
            vn1[ ioffvn1 + j - 1 ] = 0.;
            vn2[ ioffvn2 + j - 1 ] = 0.;
          }
        } else vn1[ ioffvn1 + j - 1 ] *= Math.sqrt( temp );
      }
    }
  }
}
//************************************************************************
LaPack2.dlaqps = function( m, n, offset, nb, kb, A, lda, jpvt,
tau, vn1, vn2, auxv, F, ldf, ioffa, ioffjpvt, iofftau, ioffvn1, ioffvn2,
ioffauxv, iofff ) {
  throw new Error("not tested");
  var lastrk = Math.min( m, n + offset );
  var lsticc = 0;
  var k = 0;
  var tol3z = Math.sqrt( LaPack0.dlamch( 'Epsilon' ) );
  while ( true ) { // 10
    if ( k < nb && lsticc == 0 ) {
      k ++;
      var rk = offset + k;
      var pvt = ( k - 1 )
        + Blas1.idamax( n - k + 1, vn1, 1, ioffvn1 + k - 1 );
      if ( pvt != k ) {
        Blas1.dswap( m, A, 1, A, 1, ioffa + ( pvt - 1 ) * lda,
          ioffa + ( k - 1 ) * lda );
        Blas1.dswap( k - 1, F, ldf, F, ldf, iofff + pvt - 1,
          iofff + k - 1 );
        var itemp = jpvt[ ioffjpvt + pvt - 1 ];
        jpvt[ ioffjpvt + pvt - 1 ] = jpvt[ ioffjpvt + k - 1 ];
        jpvt[ ioffjpvt + k - 1 ] = itemp;
        vn1[ ioffvn1 + pvt - 1 ] = vn1[ ioffvn1 + k - 1 ];
        vn2[ ioffvn2 + pvt - 1 ] = vn2[ ioffvn2 + k - 1 ];
      }
      if ( k > 1 ) {
        Blas2.dgemv( 'No transpose', m - rk + 1, k - 1, -1., A, lda,
          F, ldf, 1., A, 1, ioffa + rk - 1, iofff + k - 1,
          ioffa + rk - 1 + ( k - 1 ) * lda );
      }
      var alpha =
        new NumberReference( A[ ioffa + rk - 1 + ( k - 1 ) * lda ] );
      var tauref =
        new NumberReference( tau[ iofftau + k - 1 ] );
      if ( rk < m ) {
        LaPack1.dlarfg( m - rk + 1, alpha, A, 1, tauref,
          ioffa + rk + ( k - 1 ) * lda );
      } else {
        LaPack1.dlarfg( 1, alpha, A, 1, tauref,
          ioffa + rk - 1 + ( k - 1 ) * lda );
      }
      A[ ioffa + rk - 1 + ( k - 1 ) * lda ] = alpha.getValue();
      tau[ iofftau + k - 1 ] = tauref.getValue();
      var akk = A[ ioffa + rk - 1 + ( k - 1 ) * lda ];
      A[ ioffa + rk - 1 + ( k - 1 ) * lda ] = 1.;
      if ( k < n ) {
        Blas2.dgemv( 'Transpose', m - rk + 1, n - k,
          tau[ iofftau + k - 1 ], A, lda, A, 1, 0., F, 1,
          ioffa + rk - 1 + k * lda, ioffa + rk - 1 + ( k - 1 ) * lda,
          iofff + k + ( k - 1 ) * ldf );
      }
      for ( var j = 1; j <= k; j ++ ) {
        F[ iofff + j - 1 + ( k - 1 ) * ldf ] = 0.;
      }
      if ( k > 1 ) {
        Blas2.dgemv( 'Transpose', m - rk + 1, k - 1,
          - tau[ iofftau + k - 1 ], A, lda, A, 1, 0., auxv, 1,
          ioffa + rk - 1, ioffa + rk - 1 + ( k - 1 ) * lda,
          ioffauxv );
        Blas2.dgemv( 'No transpose', n, k - 1, 1., F, ldf, auxv, 1, 1.,
          F, 1, iofff, ioffauxv, iofff + ( k - 1 ) * ldf );
      }
      if ( k < n ) {
        Blas2.dgemv( 'No transpose', n - k, -1., F, ldf, A, lda, 1.,
          A, lda, iofff + k, ioffa + rk - 1,
          ioffa + rk - 1 + k * lda );
      }
      if ( rk < lastrk ) {
        for ( j = k + 1; j <= n; j ++ ) {
          if ( vn1[ ioffvn1 + j - 1 ] != 0. ) {
            var temp = 
              Math.abs( A[ ioffa + rk - 1 + ( j - 1 ) * lda ] )
              / vn1[ ioffvn1 + j - 1 ];
            temp = Math.max( 0., ( 1. + temp ) * ( 1. - temp ) );
            var temp2 = temp * Math.pow( vn1[ ioffvn1 + j - 1 ]
              / vn2[ ioffvn2 + j - 1 ], 2 );
            if ( temp2 <= tol3z ) {
              vn2[ ioffvn2 + j - 1 ] = Number( lsticc );
              lsticc = j;
            } else {
              vn1[ ioffvn1 + j - 1 ] *= Math.sqrt( temp );
            }
          }
        }
      }
      A[ ioffa + rk - 1 + ( k - 1 ) * lda ] = akk;
    } else break;
  }
  kb.setValue( k );
  rk = offset + kb.getValue();
  if ( kb.getValue() < Math.min( n, m - offset ) ) {
    Blas3.dgemm( 'No transpose', 'Transpose', m - rk, n - kb.getValue(),
      kb.getValue(), -1., A, lda, F, ldf, 1., A, lda, ioffa + rk,
      iofff + kb.getValue(), ioffa + rk + kb.getValue() * lda );
  }
  while ( lsticc > 0 ) {
    itemp = Math.round( vn2[ ioffvn2 + lsticc - 1 ] );
    vn1[ ioffvn1 + lsticc - 1 ] =
      Blas1.dnrm2( m - rk, A, 1, ioffa + rk + ( lsticc - 1 ) * lda );
    vn2[ ioffvn2 + lsticc - 1 ] = vn1[ ioffvn1 + lsticc - 1 ];
    lsticc = itemp;
  }
}
//************************************************************************
LaPack2.dlaqr5 = function( wantt, wantz, kacc22, n, ktop, kbot,
nshfts, sr, si, H, ldh, iloz, ihiz, Z, ldz, V, ldv, U, ldu, nv, WV, ldwv,
nh, WH, ldwh, ioffsr, ioffsi, ioffh, ioffz, ioffv, ioffu, ioffwv,
ioffwh ) {
  var vt = new Array( 3 );
  if ( nshfts < 2 ) return;
  if ( ktop >= kbot ) return;
  for ( var i = 1; i <= nshfts - 2; i += 2 ) {
    if ( si[ ioffsi + i - 1 ] != - si[ ioffsi + i ] ) {
      var swap = sr[ ioffsr + i - 1 ];
      sr[ ioffsr + i - 1 ] = sr[ ioffsr + i ];
      sr[ ioffsr + i ] = sr[ ioffsr + i + 1 ];
      sr[ ioffsr + i + 1 ] = swap;
      swap = si[ ioffsi + i - 1 ];
      si[ ioffsi + i - 1 ] = si[ ioffsi + i ];
      si[ ioffsi + i ] = si[ ioffsi + i + 1 ];
      si[ ioffsi + i + 1 ] = swap;
    }
  } // 10
  var ns = nshfts - nshfts % 2;
  var safmin = 
    new NumberReference( LaPack0.dlamch( 'Safe minimum' ) );
  var safmax = 
    new NumberReference( 1. / safmin.getValue() );
  LaPack0.dlabad( safmin, safmax );
  var ulp = LaPack0.dlamch( 'Precision' );
  var smlnum = safmin.getValue() * ( Number( n ) / ulp );
  var accum = ( kacc22 == 1 ) || ( kacc22 == 2 );
  var blk22 = ( ns > 2 ) && ( kacc22 == 2 );
  if ( ktop + 2 <= kbot ) {
    H[ ioffh + ktop + 1 + ( ktop - 1 ) * ldh ] = 0.;
  }
  var nbmps = ns / 2;
  var kdu = 6 * nbmps - 3;
  for ( var incol = 3 * ( 1 - nbmps ) + ktop - 1;
  incol <= kbot - 2; incol += 3 * nbmps - 2 ) { // ==> 220
    var ndcol = incol + kdu;
    if ( accum ) {
      LaPack0.dlaset( 'All', kdu, kdu, 0., 1., U, ldu, ioffu );
    }
    for ( var krcol = incol;
    krcol <= Math.min( incol + 3 * nbmps - 3, kbot - 2 ); krcol ++ )
    { // ==> 150
      var mtop =
        Math.max( 1, ( ( ktop - 1 ) - krcol + 2 ) / 3 + 1 );
      var mbot = Math.min( nbmps, ( kbot - krcol ) / 3 );
      var m22 = mbot + 1;
      var bmp22 = ( mbot < nbmps ) &&
        ( ( krcol + 3 * ( m22 - 1 ) ) == ( kbot - 2 ) );
      for ( var m = mtop; m <= mbot; m ++ ) { // ==> 20
        var k = krcol + 3 * ( m - 1 );
        if ( k == ktop - 1 ) {
          LaPack0.dlaqr1( 3, H, ldh, sr[ ioffsr + 2 * m - 2],
            si[ ioffsi + 2 * m - 2 ], sr[ ioffsr + 2 * m - 1 ],
            si[ ioffsi + 2 * m - 1 ], V,
            ioffh + ktop - 1 + ( ktop - 1 ) * ldh,
            ioffv + ( m - 1 ) * ldv );
          var alpha =
            new NumberReference( V[ ioffv + ( m - 1 ) * ldv ] );
          var tau =
            new NumberReference( V[ ioffv + ( m - 1 ) * ldv ] );
          LaPack1.dlarfg( 3, alpha, V, 1, tau,
            ioffv + 1 + ( m - 1 ) * ldv );
          V[ ioffv + ( m - 1 ) * ldv ] = tau.getValue();
        } else {
          var beta =
            new NumberReference( H[ ioffh + k + ( k - 1 ) * ldh ] );
          V[ ioffv + 1 + ( m - 1 ) * ldv ] =
            H[ ioffh + k + 1 + ( k - 1 ) * ldh ];
          V[ ioffv + 2 + ( m - 1 ) * ldv ] =
            H[ ioffh + k + 2 + ( k - 1 ) * ldh ];
          tau.setValue( V[ ioffv + ( m - 1 ) * ldv ] );
          LaPack1.dlarfg( 3, beta, V, 1, tau,
            ioffv + 1 + ( m - 1 ) * ldv );
          V[ ioffv + ( m - 1 ) * ldv ] = tau.getValue();
          if ( H[ ioffh + k + 2 + ( k - 1 ) * ldh ] != 0. ||
          H[ ioffh + k + 2 + k * ldh ] != 0. ||
          H[ ioffh + k + 2 + ( k + 1 ) * ldh ] != 0. ) {
            H[ ioffh + k + ( k - 1 ) * ldh ] = beta.getValue();
            H[ ioffh + k + 1 + ( k - 1 ) * ldh ] = 0.;
            H[ ioffh + k + 2 + ( k - 1 ) * ldh ] = 0.;
          } else {
            LaPack0.dlaqr1( 3, H, ldh, sr[ ioffsr + 2 * m - 2],
              si[ ioffsi + 2 * m - 2 ], sr[ ioffsr + 2 * m - 1 ],
              si[ ioffsi + 2 * m - 1 ], vt, ioffh + k + k * ldh, 0 );
            alpha.setValue( vt[ 0 ] );
            tau.setValue( vt[ 0 ] );
            LaPack1.dlarfg( 3, alpha, vt, 1, tau, 1 );
            vt[ 0 ] = tau.getValue();
            var refsum = vt[ 0 ]
              * ( H[ ioffh + k + ( k - 1 ) * ldh ] + vt[ 1 ] *
              H[ ioffh + k + 1 + ( k - 1 ) * ldh ] );
            if ( Math.abs( H[ ioffh + k + 1 + ( k - 1 ) * ldh ]
            - refsum * vt[ 1 ] ) + Math.abs( refsum * vt[ 2 ] ) >
            ulp * ( Math.abs( H[ ioffh + k - 1 + ( k - 1 ) * ldh ] )
            + Math.abs( H[ ioffh + k + k * ldh ] )
            + Math.abs( H[ ioffh + k + 1 + ( k + 1 ) * ldh ] ) ) ) {
              H[ ioffh + k + ( k - 1 ) * ldh ] = beta.getValue();
              H[ ioffh + k + 1 + ( k - 1 ) * ldh ] = 0.;
              H[ ioffh + k + 2 + ( k - 1 ) * ldh ] = 0.;
            } else {
              H[ ioffh + k + ( k - 1 ) * ldh ] -= refsum;
              H[ ioffh + k + 1 + ( k - 1 ) * ldh ] = 0.;
              H[ ioffh + k + 2 + ( k - 1 ) * ldh ] = 0.;
              V[ ioffv + ( m - 1 ) * ldv ] = vt[ 0 ];
              V[ ioffv + 1 + ( m - 1 ) * ldv ] = vt[ 1 ];
              V[ ioffv + 2 + ( m - 1 ) * ldv ] = vt[ 2 ];
            }
          }
        }
      } // 20
      k = krcol + 3 * ( m22 - 1 );
      if ( bmp22 ) {
        if ( k == ktop - 1 ) {
          LaPack0.dlaqr1( 2, H, ldh, sr[ ioffsr + 2 * m22 - 2 ],
            si[ ioffsi + 2 * m22 - 2 ], sr[ ioffsr + 2 * m22 - 1 ],
            si[ ioffsi + 2 * m22 - 1 ], V, ioffh + k + k * ldh,
            ioffv + ( m22 - 1 ) * ldv );
          beta.setValue( V[ ioffv + ( m22 - 1 ) * ldv ] );
          tau.setValue( V[ ioffv + ( m22 - 1 ) * ldv ] );
          LaPack1.dlarfg( 2, beta, V, 1, tau,
            ioffv + 1 + ( m22 - 1 ) * ldv );
          V[ ioffv + ( m22 - 1 ) * ldv ] = tau.getValue();
        } else {
          beta = H[ ioffh + k + ( k - 1 ) * ldh ];
          V[ ioffv + 1 + ( m22 - 1 ) * ldv ] =
            H[ ioffh + k + 1 + ( k - 1 ) * ldh ];
          tau.setValue( V[ ioffv + ( m22 - 1 ) * ldv ] );
          LaPack1.dlarfg( 2, beta, V, 1, tau,
            ioffv + 1 + ( m22 - 1 ) * ldv );
          V[ ioffv + ( m22 - 1 ) * ldv ] = tau.getValue();
          H[ ioffh + k + ( k - 1 ) * ldh ] = beta.getValue();
          H[ ioffh + k + 1 + ( k - 1 ) * ldh ] = 0.;
        }
      }
      if ( accum ) var jbot = Math.min( ndcol, kbot );
      else if ( wantt ) jbot = n;
      else jbot = kbot;
      for ( var j = Math.max( ktop, krcol ); j <= jbot; j ++ ) {
        var mend = Math.min( mbot, ( j - krcol + 2 ) / 3 );
        for ( m = mtop; m <= mend; m ++ ) {
          k = krcol + 3 * ( m - 1 );
          refsum = V[ ioffv + ( m - 1 ) * ldv ]
            * ( H[ ioffh + k + ( j - 1 ) * ldh ]
            + V[ ioffv + 1 + ( m - 1 ) * ldv ]
            * H[ ioffh + k + 1 + ( j - 1 ) * ldh ]
            + V[ ioffv + 2 + ( m - 1 ) * ldv ]
            * H[ ioffh + k + 2 + ( j - 1 ) * ldh ] );
          H[ ioffh + k + ( j - 1 ) * ldh ] -= refsum;
          H[ ioffh + k + 1 + ( j - 1 ) * ldh ] -=
            refsum * V[ ioffv + 1 + ( m - 1 ) * ldv ];
          H[ ioffh + k + 2 + ( j - 1 ) * ldh ] -=
            refsum * V[ ioffv + 2 + ( m - 1 ) * ldv ];
        } // 30
      } // 40
      if ( bmp22 ) {
        k = krcol + 3 * ( m22 - 1 );
        for ( j = Math.max( k + 1, ktop ); j <= jbot; j ++ ) {
          refsum = V[ ioffv + ( m22 - 1 ) * ldv ]
            * ( H[ ioffv + k + ( j - 1 ) * ldh ]
            + V[ ioffv + 1 + ( m22 -1 ) * ldv ]
            * H[ ioffh + k + 1 + ( j - 1 ) * ldh ] );
          H[ ioffh + k + ( j - 1 ) * ldh ] -= refsum;
          H[ ioffh + k + 1 + ( j - 1 ) * ldh ] -=
            refsum * V[ ioffv + 1 + ( m22 - 1 ) * ldv ];
        } // 50
      }
      if ( accum ) var jtop = Math.max( ktop, incol );
      else if ( wantt ) jtop = 1;
      else jtop = ktop;
      for ( m = mtop; m <= mbot; m ++ ) {
        if ( V[ ioffv + ( m - 1 ) * ldv ] != 0. ) {
          k = krcol + 3 * ( m - 1 );
          for ( j = jtop; j <= Math.min( kbot, k + 3 ); j ++ ) {
            refsum = V[ ioffv + ( m - 1 ) * ldv ]
              * ( H[ ioffh + j - 1 + k * ldh ]
              + V[ ioffv + 1 + ( m - 1 ) * ldv ]
              * H[ ioffh + j - 1 + ( k + 1 ) * ldh ]
              + V[ ioffv + 2 + ( m - 1 ) * ldv ]
              * H[ ioffv + j - 1 + ( k + 2 ) * ldh ] );
            H[ ioffh + j - 1 + k * ldh ] -= refsum;
            H[ ioffh + j - 1 + ( k + 1 ) * ldh ] -=
              refsum * V[ ioffv + 1 + ( m - 1 ) * ldv ];
            H[ ioffh + j - 1 + ( k + 2 ) * ldh ] -=
              refsum * V[ ioffv + 2 + ( m - 1 ) * ldv ];
          }  // 60
          if ( accum ) {
            var kms = k - incol;
            for ( j = Math.max( 1, ktop - incol ); j <= kdu; j ++ ) {
              refsum = V[ ioffv + ( m - 1 ) * ldv ]
                * ( U[ ioffu + j - 1 + kms * ldu ]
                + V[ ioffv + 1 + ( m - 1 ) * ldv ]
                * U[ ioffu + j - 1 + ( kms + 1 ) * ldu ]
                + V[ ioffv + 2 + ( m - 1 ) * ldv ]
                * U[ ioffu + j - 1 + ( kms + 2 ) * ldu ] );
              U[ ioffu + j - 1 + kms * ldu ] -= refsum;
              U[ ioffu + j - 1 + ( kms + 1 ) * ldu ] -=
                refsum * V[ ioffv + 1 + ( m - 1 ) * ldv ];
              U[ ioffu + j - 1 + ( kms + 2 ) * ldu ] -=
                refsum * V[ ioffv + 2 + ( m - 1 ) * ldv ];
            } // 70
          } else if ( wantz ) {
            for ( j = iloz; j <= ihiz; j ++ ) {
              refsum = V[ ioffv + ( m - 1 ) * ldv ]
                * ( Z[ ioffu + j - 1 + k * ldz ]
                + V[ ioffv + 1 + ( m - 1 ) * ldv ]
                * Z[ ioffz + j - 1 + ( k + 1 ) * ldz ]
                + V[ ioffv + 2 + ( m - 1 ) * ldv ]
                * Z[ ioffz + j - 1 + ( k + 2 ) * ldz ] );
              Z[ ioffz + j - 1 + k * ldz ] -= refsum;
              Z[ ioffz + j - 1 + ( k + 1 ) * ldz ] -=
                refsum * V[ ioffv + 1 + ( m - 1 ) * ldv ];
              Z[ ioffz + j - 1 + ( k + 2 ) * ldz ] -=
                refsum * V[ ioffv + 2 + ( m - 1 ) * ldv ];
            } // 80
          }
        }
      } // 90
      k = krcol + 3 * ( m22 - 1 );
      if ( bmp22 ) {
        if ( V[ ioffv + ( m22 - 1 ) * ldv ] != 0. ) {
          for ( j = jtop; j <= Math.min( kbot, k + 3 ); j ++ ) {
            refsum = V[ ioffv + ( m22 - 1 ) * ldv ]
              * ( H[ ioffh + j - 1 + k * ldh ]
              + V[ ioffv + 1 + ( m22 - 1 ) * ldv ]
              * H[ ioffh + j - 1 + ( k + 1 ) * ldh ] );
            H[ ioffh + j - 1 + k * ldh ] -= refsum;
            H[ ioffh + j - 1 + ( k + 1 ) * ldh ] -=
              refsum * V[ ioffv + 1 + ( m22 - 1 ) * ldv ];
          } // 100
          if ( accum ) {
            kms = k - incol;
            for ( j = Math.max( 1, ktop - incol ); j <= kdu; j ++ ) {
              refsum = V[ ioffv + ( m22 - 1 ) * ldv ]
                * ( U[ ioffu + j - 1 + kms * ldu ]
                + V[ ioffv + 1 + ( m22 - 1 ) * ldv ]
                * U[ ioffu + j - 1 + ( kms + 1 ) * ldu ] );
              U[ ioffu + j - 1 + kms * ldu ] -= refsum;
              U[ ioffu + j - 1 + ( kms + 1 ) * ldu ] -=
                refsum * V[ ioffv + 1 + ( m22 - 1 ) * ldv ];
            } // 110
          } else if ( wantz ) {
            for ( j = iloz; j <= ihiz; j ++ ) {
              refsum = V[ ioffv + ( m22 - 1 ) * ldv ]
                *( Z[ ioffz + j - 1 + k * ldz ]
                + V[ ioffv + 1 + ( m22 - 1 ) * ldv ]
                * Z[ ioffz + j - 1 + ( k + 1 ) * ldz ] );
              Z[ ioffz + j - 1 + k * ldz ] -= refsum;
              Z[ ioffz + j - 1 + ( k + 1 ) * ldz ] -=
                refsum * V[ ioffv + 1 + ( m22 - 1 ) * ldv ];
            } // 120
          }
        }
      }
      var mstart = mtop;
      if ( krcol + 3 * ( mstart - 1 ) < ktop ) mstart ++;
      mend = mbot;
      if ( bmp22 ) mend ++;
      if ( krcol == kbot - 2 ) mend ++;
      for ( m = mstart; m <= mend; m ++ ) {
        k = Math.min( kbot - 1, krcol + 3 * ( m - 1 ) );
        if ( H[ ioffh + k + ( k - 1 ) * ldh ] != 0. ) {
          var tst1 =
            Math.abs( H[ ioffh + k - 1 + ( k - 1 ) * ldh ] )
            + Math.abs( H[ ioffh + k + k * ldh ] );
          if ( tst1 == 0. ) {
            if ( k >= ktop + 1 ) {
              tst1 += Math.abs( H[ ioffh + k - 1 + ( k - 2 ) * ldh ] );
            }
            if ( k >= ktop + 2 ) {
              tst1 += Math.abs( H[ ioffh + k - 1 + ( k - 3 ) * ldh ] );
            }
            if ( k >= ktop + 3 ) {
              tst1 += Math.abs( H[ ioffh + k - 1 + ( k - 4 ) * ldh ] );
            }
            if ( k <= kbot - 2 ) {
              tst1 += Math.abs( H[ ioffh + k + 1 + k * ldh ] );
            }
            if ( k <= kbot - 3 ) {
              tst1 += Math.abs( H[ ioffh + k + 2 + k * ldh ] );
            }
            if ( k <= kbot - 4 ) {
              tst1 += Math.abs( H[ ioffh + k + 3 + k * ldh ] );
            }
          }
          if ( Math.abs( H[ ioffh + k + ( k - 1 ) * ldh ] )
          <= Math.max( smlnum, ulp * tst1 ) ) {
            var h12 = Math.max(
              Math.abs( H[ ioffh + k + ( k - 1 ) * ldh ] ),
              Math.abs( H[ ioffh + k - 1 + k * ldh ] ) );
            var h21 = Math.min(
              Math.abs( H[ ioffh + k + ( k - 1 ) * ldh ] ),
              Math.abs( H[ ioffh + k - 1 + k * ldh ] ) );
            var h11 = Math.max(
              Math.abs( H[ ioffh + k + k * ldh ] ),
              Math.abs( H[ ioffh + k - 1 + ( k - 1 ) * ldh ]
              - H[ ioffh + k + k * ldh ] ) );
            var h22 = Math.min(
              Math.abs( H[ ioffh + k + k * ldh ] ),
              Math.abs( H[ ioffh + k - 1 + ( k - 1 ) * ldh ]
              - H[ ioffh + k + k * ldh ] ) );
            var scl = h11 + h12;
            var tst2 = h22 * ( h11 / scl );
            if ( tst2 == 0. ||
            h21 * ( h12 / scl ) <= Math.max( smlnum, ulp * tst2 ) ) {
              H[ ioffh + k + ( k - 1 ) * ldh ] = 0.;
            }
          }
        }
      } // 130
      mend = Math.min( nbmps, ( kbot - krcol - 1 ) / 3 );
      for ( m = mtop; m <= mend; m ++ ) {
        k = krcol + 3 * ( m - 1 );
        refsum = V[ ioffv + ( m - 1 ) * ldv ]
          * V[ ioffv + 2 + ( m - 1 ) * ldv ]
          * H[ ioffh + k + 3 + ( k + 2 ) * ldh ];
        H[ ioffh + k + 3 + k * ldh ] -= refsum;
        H[ ioffh + k + 3 + ( k + 1 ) * ldh ] -=
          refsum * V[ ioffv + 1 + ( m - 1 ) * ldv ]; 
        H[ ioffh + k + 3 + ( k + 2 ) * ldh ] -=
          refsum * V[ ioffv + 2 + ( m - 1 ) * ldv ]; 
      } // 140
    } // 150
    if ( accum ) {
      if ( wantt ) {
        jtop = 1;
        jbot = n;
      } else {
        jtop = ktop;
        jbot = kbot;
      }
      if ( ! blk22 || incol < ktop || ndcol > kbot || ns <= 2 ) {
        var k1 = Math.max( 1, ktop - incol );
        var nu = ( kdu - Math.max( 0, ndcol - kbot ) ) - k1 + 1;
        for ( var jcol = Math.min( ndcol, kbot ) + 1; jcol <= jbot;
        jcol += nh ) {
          var jlen = Math.min( nh, jbot - jcol + 1 );
          Blas3.dgemm( 'C', 'N', nu, jlen, nu, 1., U, ldu, H, ldh, 0.,
            WH, ldwh, ioffu + k1 - 1 + ( k1 - 1 ) * ldu,
            ioffh + incol + k1 - 1 + ( jcol - 1 ) * ldh, ioffwh );
          LaPack0.dlacpy( 'All', nu, jlen, WH, ldwh, H, ldh,
            ioffwh, ioffh + incol + k1 - 1 + ( jcol - 1 ) * ldh );
        } // 160
        for ( var jrow = jtop; jrow <= Math.max( ktop, incol ) - 1;
        jrow += nv ) {
          jlen = Math.min( nv, Math.max( ktop, incol ) - jrow );
          Blas3.dgemm( 'N', 'N', jlen, nu, nu, 1., H, ldh, U, ldu, 0.,
            WV, ldwv, ioffh + jrow - 1 + ( incol + k1 - 1 ) * ldh,
            ioffu + k1 - 1 + ( k1 - 1 ) * ldu, ioffwv );
          LaPack0.dlacpy( 'All', jlen, nu, WV, ldwv, H, ldh,
            ioffwv, ioffh + jrow - 1 + ( incol + k1 - 1 ) * ldh );
        } // 170
        if ( wantz ) {
          for ( jrow = iloz; jrow <= ihiz; jrow += nv ) {
            jlen = Math.min( nv, ihiz - jrow + 1 );
            Blas3.dgemm( 'N', 'N', jlen, nu, nu, 1., Z, ldz, U, ldu,
              0., WV, ldwv,
              ioffz + jrow - 1 + ( incol + k1 - 1 ) * ldz,
              ioffu + k1 - 1 + ( k1 - 1 ) * ldu, ioffwv );
            LaPack0.dlacpy( 'All', jlen, nu, WV, ldwv, Z, ldz,
              ioffwv, ioffz + jrow - 1 + ( incol + k1 - 1 ) * ldz );
          } // 180 
        }
      } else {
        var i2 = ( kdu + 1 ) / 2;
        var i4 = kdu;
        var j2 = i4 - i2;
        var j4 = kdu;
        var kzs = ( j4 - j2 ) - ( ns + 1 );
        var knz = ns + 1;
        for ( jcol = Math.min( ndcol, kbot ) + 1; jcol <= jbot;
        jcol += nh ) {
          jlen = Math.min( nh, jbot - jcol + 1 );
          LaPack0.dlacpy( 'All', knz, jlen, H, ldh, WH, ldwh,
            ioffh + incol + j2 + ( jcol - 1 ) * ldh, ioffwh + kzs );
          LaPack0.dlaset( 'All', kzs, jlen, 0., 0., WH, ldwh, ioffwh );
          Blas3.dtrmm( 'L', 'U', 'C', 'N', knz, jlen, 1., U, ldu,
            WH, ldwh, ioffu + j2 + kzs * ldu, ioffwh + kzs );
          Blas3.dgemm( 'C', 'N', i2, jlen, j2, 1., U, ldu, H, ldh, 1.,
            WH, ldwh, ioffu, ioffh + incol + ( jcol - 1 ) * ldh,
            ioffwh );
          LaPack0.dlacpy( 'All', j2, jlen, H, ldh, WH, ldwh,
            ioffh + incol + ( jcol -1 ) * ldh, ioffwh + i2 );
          Blas3.dtrmm( 'L', 'L', 'C', 'N', j2, jlen, 1., U, ldu,
            WH, ldwh, ioffu + i2 * ldu, ioffwh + i2 );
          Blas3.dgemm( 'C', 'N', i4 - i2, jlen, j4 - j2, 1., U, ldu,
            H, ldh, 1., WH, ldwh, ioffu + j2 + i2 * ldu,
            ioffh + incol + j2 + ( jcol - 1 ) * ldh, ioffwh + i2 );
          LaPack0.dlacpy( 'All', kdu, jlen, WH, ldwh, H, ldh,
            ioffwh, ioffh + incol + ( jcol - 1 ) * ldh );
        } // 190
        for ( jrow = jtop; jrow <= Math.max( incol, ktop ) - 1;
        jrow += nv ) {
          jlen = Math.min( nv, Math.max( incol, ktop ) - jrow );
          LaPack0.dlacpy( 'All', jlen, knz, H, ldh, WV, ldwv,
            ioffh + jrow - 1 + ( incol + j2 ) * ldh,
            ioffwv + kzs * ldwv );
          LaPack0.dlaset( 'All', jlen, kzs, 0., 0., WV, ldwv, ioffwv );
          Blas3.dtrmm( 'R', 'U', 'N', 'N', jlen, knz, 1., U, ldu,
            WV, ldwv, ioffu + j2 + kzs * ldu, ioffwv + kzs * ldwv );
          Blas3.dgemm( 'N', 'N', jlen, i2, j2, 1., H, ldh, U, ldu, 1.,
            WV, ldwv, ioffh + jrow - 1 + incol * ldh, ioffu, ioffwv );
          LaPack0.dlacpy( 'All', jlen, j2, H, ldh, WV, ldwv,
            ioffh + jrow - 1 + incol * ldh, ioffwv + i2 * ldwv );
          Blas3.dtrmm( 'R', 'L', 'N', 'N', jlen, i4 - i2, 1., U, ldu,
            WV, ldwv, ioffu + i2 * ldu, ioffwv + i2 * ldwv );
          Blas3.dgemm( 'N', 'N', jlen, i4 - i2, j4 - j2, 1., H, ldh,
            U, ldu, 1., WV, ldwv,
            ioffh + jrow - 1 + ( incol + j2 ) * ldh,
            ioffu + j2 + i2 * ldu, ioffwv + i2 * ldwv );
          LaPack0.dlacpy( 'All', jlen, kdu, WV, ldwv, H, ldh,
            ioffwv, ioffh + jrow - 1 + incol * ldh );
        } // 200
        if ( wantz ) {
          for ( jrow = iloz; jrow <= ihiz; jrow += nv ) {
            jlen = Math.min( nv, ihiz - jrow + 1 );
            LaPack0.dlacpy( 'All', jlen, knz, Z, ldz, WV, ldwv,
              ioffz + jrow - 1 + ( incol + j2 ) * ldz,
              ioffwv + kzs * ldwv );
            LaPack0.dlaset( 'All', jlen, kzs, 0., 0., WV, ldwv,
              ioffwv );
            Blas3.dtrmm( 'R', 'U', 'N', 'N', jlen, knz, 1., U, ldu,
              WV, ldwv, ioffu + j2 + kzs * ldu, ioffwv + kzs * ldwv );
            Blas3.dgemm( 'N', 'N', jlen, i2, j2, 1., Z, ldz, U, ldu,
              1., WV, ldwv, ioffz + jrow - 1 + incol * ldz, ioffu,
              ioffwv );
            LaPack0.dlacpy( 'All', jlen, j2, Z, ldz, WV, ldwv,
              ioffz + jrow - 1 + incol * ldz, ioffwv + i2 * ldwv );
            Blas3.dtrmm( 'R', 'L', 'N', 'N', jlen, i4 - i2, 1., U, ldu,
              WV, ldwv, ioffu + i2 * ldu, ioffwv + i2 * ldwv );
            Blas3.dgemm( 'N', 'N', jlen, i4 - i2, j4 - j2, 1., Z, ldz,
              U, ldu, 1., WV, ldwv,
              ioffz + jrow - 1 + ( incol + j2 ) * ldz,
              ioffu + j2 + i2 * ldu, ioffwv + i2 * ldwv );
            LaPack0.dlacpy( 'All', jlen, kdu, WV, ldwv, Z, ldz,
              ioffwv, ioffz + jrow - 1 + incol * ldz );
          } // 210
        }
      }
    }
  } // 220
}
//************************************************************************
LaPack2.dlaqtr = function( ltran, lreal, n, T, ldt, b, w,
scale, x, work, info, iofft, ioffb, ioffx, ioffwork ) {
  throw new Error("not programmed");
}
//************************************************************************
LaPack2.dlarfx = function( side, m, n, v, tau, C, ldc, work,
ioffv, ioffc, ioffwork ) {
  if ( tau == 0. ) return;
  if ( side.charAt(0).toUpperCase() == 'L' ) {
    LaPack1.dlarf( side, m, n, v, 1, tau, C, ldc, work,
      ioffv, ioffc, ioffwork );
//      Blas2.dgemv( 'Transpose', m, n, 1., C, ldc, v, 1, 0., work, 1,
//        ioffc, ioffv, ioffwork );
//      Blas2.dger( m, n, -tau, v, 1, work, 1, C, ldc, ioffv,
//        ioffwork, ioffc );
  } else {
    LaPack1.dlarf( side, m, n, v, 1, tau, C, ldc, work,
      ioffv, ioffc, ioffwork );
//      Blas2.dgemv( 'No transpose', m, n, 1., C, ldc, v, 1, 0., work,
//        1, ioffc, ioffv, ioffwork );
//      Blas2.dger( m, n, -tau, work, 1, v, 1, C, ldc, ioffwork,
//        ioffv, ioffc );
  }
}
LaPack2.zlarfx = function( side, m, n, v, tau, C, ldc, work ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dlarrv = function( n, vl, vu, d, l, pivmin, isplit, m,
dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers, Z,
ldz, isuppz, work, iwork, info, ioffd, ioffl, ioffisplit, ioffw,
ioffwerr, ioffwgap, ioffiblock, ioffindexw, ioffgers, ioffz, ioffisuppz,
ioffwork, ioffiwork ) {
  throw new Error("not tested");
  var maxitr = 10;
  var indld = n + 1;
  var indlld = 2 * n + 1;
  var indwrk = 3 * n + 1;
  var minwsize = 12 * n;
  for ( var i = 1; i <= minwsize; i ++ ) {
    work[ ioffwork + i - 1 ] = 0.;
  } // 5
  var iindr = 0;
  var iindc1 = n;
  var iindc2 = 2 * n;
  var iindwk = 3 * n + 1;
  var miniwsize = 7 * n;
  for ( i = 1; i <= miniwsize; i ++ ) iwork[ ioffiwork + i - 1 ] = 0.;
  var zusedl = 1;
  if ( dol > 1 ) zusedl = dol - 1;
  var zusedu = m;
  if ( dou < m ) zusedu = dou + 1;
  var zusedw = zusedu - zusedl + 1;
  LaPack0.dlaset( 'Full', n, zusedw, 0., 0., Z, ldz,
    ioffz + ( zusedl - 1 ) * ldz );
  var eps = LaPack0.dlamch( 'Precision' );
  var rqtol = 2. * eps;
  var tryrqc = true;
  if ( dol == 1 && dou == m );
  else {
    var rtol1 = 4. * eps;
    var rtol2 = 4. * eps;
  }
  var done = 0;
  var ibegin = 1;
  var wbegin = 1;
  for ( var jblk = 1; jblk <= iblock[ ioffiblock + m - 1 ];
  jblk ++ ) { // ==> 170
    var iend = isplit[ ioffisplit + jblk - 1 ];
    var sigma = l[ ioffl + iend - 1 ];
    var wend = wbegin - 1;
    while ( true ) { // 15
      if ( wend < m ) {
        if ( iblock[ ioffiblock + wend ] == jblk ) wend ++;
        else break;
      } else break;
    }
    if ( wend < wbegin) {
      ibegin = iend + 1;
      continue; // goto 170
    } else if ( wend < dol || wbegin > dou ) {
      ibegin = iend + 1;
      wbegin = wend + 1;
      continue; // goto 170
    }
    var gl = gers[ ioffgers + 2 * ibegin - 2 ];
    var gu = gers[ ioffgers + 2 * ibegin - 1 ];
    for ( i = ibegin + 1; i <= iend; i ++ ) {
      gl = Math.min( gers[ ioffgers + 2 * i - 2 ], gl );
      gu = Math.max( gers[ ioffgers + 2 * i - 1 ], gu );
    } // 20
    var spdiam = gu - gl;
    var oldien = ibegin - 1;
    var in2 = iend - ibegin + 1;
    var im = wend - wbegin + 1;
    if ( ibegin == iend ) {
      done ++;
      Z[ ioffz + ibegin - 1 + ( wbegin - 1 ) * ldz ] = 1.;
      isuppz[ ioffisuppz + 2 * wbegin - 2 ] = ibegin;
      isuppz[ ioffisuppz + 2 * wbegin - 1 ] = ibegin;
      w[ ioffw + wbegin - 1 ] += sigma;
      work[ ioffwork + wbegin - 1 ] = w[ ioffw + wbegin - 1 ];
      ibegin = iend + 1;
      wbegin ++;
      continue; // goto 170
    }
    Blas1.dcopy( im, w, 1, work, 1,
      ioffw + wbegin - 1, ioffwork + wbegin - 1 );
    for ( i = 1; i <= im; i ++ ) { // 30
      w[ ioffw + wbegin + i - 2 ] += sigma;
    }
    var ndepth = 0;
    var parity = 1;
    var nclus = 1;
    iwork[ ioffiwork + iindc1 ] = 1;
    iwork[ ioffiwork + iindc1 + 1 ] = im;
    var idone = 0;
    while ( idone < im ) { // 40
      if ( ndepth > m ) {
        info.setValue( -2 );
        return;
      }
      var oldncl = nclus;
      nclus = 0;
      parity = 1 - parity;
      if ( parity == 0 ) {
        var oldcls = iindc1;
        var newcls = iindc2;
      } else {
        oldcls = iindc2;
        newcls = iindc1;
      }
      for ( i = 1; i <= oldncl; i ++ ) { // ==> 150
        var j = oldcls + 2 * i;
        var oldfst = iwork[ ioffiwork + j - 2 ];
        var oldlst = iwork[ ioffiwork + j - 1 ];
        if ( ndepth > 0 ) {
          if ( dol == 1 && dou == m ) j = wbegin + oldfst - 1;
          else {
            if ( wbegin + oldfst - 1 < dol ) j = dol - 1;
            else if ( wbegin + oldfst - 1 > dou ) j = dou;
            else j = wbegin + oldfst - 1;
          }
          Blas1.dcopy( in2, Z, 1, d, 1,
            ioffz + ibegin - 1 + ( j - 1 ) * ldz, ioffd + ibegin - 1 );
          Blas1.dcopy( in2 - 1, Z, 1, l, 1,
            ioffz + ibegin - 1 + j * ldz, ioffl + ibegin - 1 );
          sigma = Z[ ioffz + iend - 1 + j * ldz ];
          LaPack0.dlaset( 'Full', in2, 2, 0., 0., Z, ldz,
            ioffz + ibegin - 1 + ( j - 1 ) * ldz );
        } // ndepth > 0
        for ( j = ibegin; j <= iend - 1; j ++ ) {
          var tmp = d[ ioffd + j - 1 ] * l[ ioffl + j - 1 ];
          work[ ioffwork + indld + j - 2 ] = tmp;
          work[ ioffwork + indlld + j - 2 ] = tmp * l[ ioffl + j - 1 ];
        } // 50
        if ( ndepth > 0 ) {
          var p = indexw[ ioffindexw + wbegin + oldfst - 2 ];
          var q = indexw[ ioffindexw + wbegin + oldlst - 2 ];
          var offset = indexw[ ioffindexw + wbegin - 1 ] - 1;
          var iinfo = new IntReference();
          LaPack1.dlarrb( in2, d, work, p, q, rtol1, rtol2, offset,
            work, wgap, werr, work, iwork, pivmin, spdiam, in2, iinfo,
            ioffd + ibegin - 1, ioffwork + indlld + ibegin - 2,
            ioffwork + wbegin - 1, ioffwgap + wbegin - 1,
            ioffwerr + wbegin - 1, ioffwork + indwrk - 1,
            ioffiwork + iindwk - 1 );
          if ( iinfo.getValue() != 0 ) {
            info.setValue( -1 );
            return;
          }
          if ( oldfst > 1 ) {
            wgap[ ioffwgap + wbegin + oldfst - 3 ] =
              Math.max( wgap[ ioffwgap + wbegin + oldfst - 3 ],
              w[ ioffw + wbegin + oldfst - 2 ]
              - werr[ ioffwerr + wbegin + oldfst - 2 ]
              - w[ ioffw + wbegin + oldfst - 3 ]
              - werr[ ioffwerr + wbegin + oldfst - 3 ] );
          }
          if ( wbegin + oldlst - 1 < wend ) {
            wgap[ ioffwgap + wbegin + oldlst - 2 ] =
              Math.max( wgap[ ioffwgap + wbegin +oldlst - 2 ],
              w[ ioffw + wbegin + oldlst - 1 ]
              - werr[ ioffwerr + wbegin + oldlst - 1 ]
              - w[ ioffw + wbegin + oldlst - 2 ]
              - werr[ ioffwerr + wbegin + oldlst - 2 ] );
          }
          for ( j = oldfst; j <= oldlst; j ++ ) {
            w[ ioffw + wbegin + j - 2 ] += sigma;
          } // 53
        } // ndepth > 0
        var newfst = oldfst;
        for ( j = oldfst; j <= oldlst; j ++ ) { // ==> 140
          var goto139 = false;
          if ( j == oldlst ) var newlst = j;
          else if ( wgap[ ioffwgap + wbegin + j - 2 ] >=
          minrgp * Math.abs( work[ ioffwork + wbegin + j - 2 ] ) ) {
            newlst = j;
          } else continue;
          var newsiz = newlst - newfst + 1;
          if ( dol == 1 && dou == m ) {
            var newftt = wbegin + newfst - 1;
          } else {
            if ( wbegin + newfst - 1 < dol ) newftt = dol - 1;
            else if ( wbegin + newfst - 1 > dou ) newftt = dou;
            else newftt = wbegin + newfst - 1;
          }
          if ( newsiz > 1 ) {
            if ( newfst == 1 ) {
              var lgap = Math.max( 0.,
                w[ ioffw + wbegin - 1 ] - werr[ ioffwerr + wbegin - 1 ]
                - vl );
            } else lgap = wgap[ ioffwgap + wbegin + newfst - 3 ];
            var rgap = wgap[ ioffwgap + wbegin + newlst - 2 ];
            for ( var k = 1; k <= 2; k ++ ) { // ==> 55
              p = ( k == 1 ? indexw[ ioffindexw + wbegin + newfst - 2 ]
                : indexw[ ioffindexw + wbegin + newlst - 2 ] );
              offset = indexw[ ioffindexw + wbegin - 1 ] -1;
              LaPack1.dlarrb( in2, d, work, p, p, rqtol, rqtol, offset,
                work, wgap, werr, work, iwork, pivmin, spdiam, in2,
                iinfo, ioffd + ibegin - 1,
                ioffwork + indlld + ibegin - 2,
                ioffwork + wbegin - 1, ioffwgap + wbegin - 1,
                ioffwerr + wbegin - 1, ioffwork + indwrk - 1,
                ioffiwork + iindwk - 1 );
            } // 55
            if ( wbegin + newlst -1 < dol ||
            wbegin + newfst - 1 > dou ) {
              idone += newlst - newfst + 1;
              goto139 = true;
            }
            if ( ! goto139 ) {
              var tau = new NumberReference();
              LaPack1.dlarrf( in2, d, l, work, newfst, newlst, work,
                wgap, werr, spdiam, lgap, rgap, pivmin, tau, Z, Z,
                work, iinfo, ioffd + ibegin - 1, ioffl + ibegin - 1,
                ioffwork + indld + ibegin - 2, ioffwork + wbegin - 1,
                wgap + wbegin - 1, werr + wbegin - 1,
                ioffz + ibegin - 1 + ( newftt - 1 ) * ldz,
                ioffz + ibegin - 1 + newftt * ldz,
                ioffwork + indwrk - 1 );
              if ( iinfo.getValue() == 0 ) {
                var ssigma = sigma + tau.getValue();
                Z[ ioffz + iend - 1 + newftt * ldz ] = ssigma;
                for ( k = newfst; k <= newlst; k ++ ) { // ==> 116
                  var fudge = 3. * eps
                    * Math.abs( work[ ioffwork + wbegin + k - 2 ] );
                  work[ ioffwork + wbegin + k - 2 ] -= tau.getValue();
                  fudge += 4. * eps
                    * Math.abs( work[ ioffwork + wbegin + k - 2 ] );
                  werr[ ioffwerr + wbegin + k - 2 ] += fudge;
                } // 116
                nclus ++;
                k + newcls + 2 * nclus;
                iwork[ ioffiwork + k - 2 ] = newfst;
                iwork[ ioffiwork + k - 1 ] = newlst;
              } else {
                info.setValue( -2 );
                return;
              } // iinfo == 0
            }
          } else {
            var iter = 0;
            var tol = 4. * Math.log( Number( in2 ) ) * eps;
            k = newfst;
            var windex = wbegin + k - 1;
            var windmn = Math.max( windex - 1, 1 );
            var windpl = Math.min( windex + 1, m );
            var lambda = work[ ioffwork + windex - 1 ];
            done ++;
            var goto125 = false;
            if ( windex < dol || windex > dou ) {
              var eskip = true;
              goto125 = true;
            } else eskip = false;
            if ( ! goto125 ) {
              var left = work[ ioffwork + windex - 1 ]
                - werr[ ioffwerr + windex - 1 ];
              var right = work[ ioffwork + windex - 1 ]
                + werr[ ioffwerr + windex - 1 ];
              var indeig = indexw[ ioffindexw + windex - 1 ];
              lgap = ( k == 1 ?
                eps * Math.max( Math.abs( left ), Math.abs( right ) ) :
                wgap[ ioffwgap + windmn - 1 ] );
              rgap = ( k == im ?
                eps * Math.max( Math.abs( left ), Math.abs( right ) ) :
                wgap[ ioffwgap + windex - 1 ] );
              var gap = Math.min( lgap, rgap );
              var gaptol =
                ( k == 1 || k == im  ? 0. : gap * eps );
              var isupmn = in2;
              var isupmx = 1;
              var savgap = wgap[ ioffwgap + windex - 1 ];
              wgap[ ioffwgap + windex - 1 ] = gap;
              var usedbs = false;
              var usedrq = false;
              var needbs = ! tryrqc;
              while ( true ) { // 120
                if ( needbs ) {
                  usedbs = true;
                  var itmp1 =
                    iwork[ ioffiwork + iindr + windex - 1 ];
                  offset = indexw[ ioffindexw + wbegin - 1 ] - 1;
                  LaPack1.dlarrb( in2, d, work, indeig, indeig, 0.,
                    2. * eps, offset, work, wgap, werr, work, iwork,
                    pivmin, spdiam, itmp1, iinfo, ioffd + ibegin - 1,
                    ioffwork + indlld + ibegin - 2,
                    ioffwork + wbegin - 1, ioffwgap + wbegin - 1,
                    ioffwerr + wbegin - 1, ioffwork + indwrk - 1,
                    ioffiwork + iindwk - 1 );
                  if ( iinfo.getValue() != 0 ) {
                    info.setValue( -3 );
                    return;
                  }
                  lambda = work [ ioffwork + windex - 1 ];
                  iwork[ ioffiwork + iindr + windex - 1 ] = 0;
                } // needbs
                var mingma = new NumberReference();
                var nrminv = new NumberReference();
                var resid = new NumberReference();
                var rqcorr = new NumberReference();
                var ztz = new NumberReference();
                var negcnt = new IntReference();
                var r = new IntReference();
                LaPack1.dlar1v( in2, 1, in2, lambda, d, l, work, work,
                  pivmin, gaptol, Z, ! usedbs, negcnt, ztz, mingma, r,
                  isuppz, nrminv, resid, rqcorr, work,
                  ioffd + ibegin - 1, ioffl + ibegin - 1,
                  ioffwork + indld + ibegin - 2,
                  ioffwork + indlld + ibegin - 2,
                  ioffz + ibegin - 1 + ( windex - 1 ) * ldz,
                  ioffisuppz + 2 * windex - 2, ioffwork + indwrk - 1 );
                if ( iter == 0 ) {
                  var bstres = resid.getValue();
                  var bstw = lambda;
                } else if ( resid.getValue() < bstres ) {
                  bstres = resid.getValue();
                  bstw = lambda;
                }
                iwork[ ioffiwork + iindr + windex - 1 ] = r.getValue();
                isupmn = Math.min( isupmn,
                  isuppz[ ioffisuppz + 2 * windex - 2 ] );
                isupmx = Math.max( isupmx,
                  isuppz[ ioffisuppz + 2 * windex - 1 ] );
                iter ++;
                if ( resid.getValue() > tol * gap &&
                Math.abs( rqcorr.getValue() ) > rqtol * Math.abs( lambda )
                && ! usedbs ) {
                  var sgndef =
                    ( indeig <= negcnt.getValue() ? -1. : 1. );
                  if ( rqcorr.getValue() * sgndef >= 0. &&
                  lambda + rqcorr.getValue() <= right &&
                  lambda + rqcorr.getValue() >= left ) {
                    usedrq = true;
                    if ( sgndef == 1. ) left = lambda;
                    else right = lambda;
                    work[ ioffwork + windex - 1 ] =
                      0.5 * ( right + left );
                    lambda += rqcorr.getValue();
                    werr[ ioffwerr + windex - 1 ] =
                      0.5 * ( right - left );
                  } else needbs = true;
                  if ( right - left < rqtol * Math.abs( lambda ) ) {
                    usedbs = true;
                    continue; // goto 120
                  } else if ( iter < maxitr ) continue; // goto 120
                  else if ( iter == maxitr ) {
                    needbs = true;
                    continue; // goto 120
                  } else {
                    info.setValue( 5 );
                    return;
                  }
                } else {
                  var stp2ii = false;
                  if ( usedrq && usedbs && bstres <= resid.getValue() ) {
                    lambda = bstw;
                    stp2ii = true;
                  }
                  if ( stp2ii ) {
                    LaPack1.dlar1v( in2, 1, in2, lambda, d, l, work,
                      work, pivmin, gaptol, Z, ! usedbs, negcnt, ztz,
                      mingma, r, 
                      isuppz, nrminv, resid, rqcorr, work,
                      ioffd + ibegin - 1, ioffl + ibegin - 1,
                      ioffwork + indld + ibegin - 2,
                      ioffwork + indlld + ibegin - 2,
                      ioffz + ibegin - 1 + ( windex - 1 ) * ldz,
                      ioffisuppz + 2 * windex - 2,
                      ioffwork + indwrk - 1 );
                    iwork[ ioffiwork + iindr + windex - 1 ] =
                      r.getValue();
                  }
                  work[ ioffwork + windex - 1 ] = lambda;
                  break;
                }
              } // 120
              isuppz[ ioffisuppz + 2 * windex - 2 ] += oldien;
              isuppz[ ioffisuppz + 2 * windex - 1 ] += oldien;
              var zfrom = isuppz[ ioffisuppz + 2 * windex - 2 ];
              var zto = isuppz[ ioffisuppz + 2 * windex - 1 ];
              isupmn += oldien;
              isupmx += oldien;
              if ( isupmn < zfrom ) {
                for ( var ii = isupmn; ii <= zfrom - 1; ii ++ ) {
                  Z[ ioffz + ii - 1 + ( windex - 1 ) * ldz ] = 0.;
                } // 122
              }
              if ( isupmx < zto ) {
                for ( ii = zto + 1; ii <= isupmx; ii ++ ) {
                  Z[ ioffz + ii - 1 + ( windex - 1 ) * ldz ] = 0.;
                } // 123
              }
              Blas1.dscal( zto - zfrom + 1, nrminv.getValue(), Z, 1,
                ioffz + zfrom - 1 + ( windex - 1 ) * ldz );
            } // 125
            w[ ioffw + windex - 1 ] = lambda + sigma;
            if ( ! eskip ) {
              if ( k > 1 ) {
                wgap[ ioffwgap + windmn - 1 ] = Math.max(
                  wgap[ ioffwgap + windmn - 1 ],
                  w[ ioffw + windex - 1 ]
                  - werr[ ioffwerr + windex - 1 ]
                  - w[ ioffw + windmn - 1 ]
                  - werr[ ioffwerr + windmn - 1 ] );
              }
              if ( windex < wend ) {
                wgap[ ioffwgap + windex - 1 ] = Math.max( savgap,
                  w[ ioffw + windpl - 1 ]
                  - werr[ ioffwerr + windpl - 1 ]
                  - w[ ioffw + windex - 1 ]
                  - werr[ ioffwerr + windex - 1 ] );
              }
            }
            idone ++;
          } // newsiz > 1; 139
          newfst = j + 1;
        } // 140
      } // 150
      ndepth ++;
    } // goto 40
    ibegin = iend + 1;
    wbegin = wend + 1;
  } // 170
}
//************************************************************************
LaPack2.dlartgs = function( x, y, sigma, cs, sn ) {
  var thresh = LaPack0.dlamch( 'E' );
  if ( ( sigma == 0. && Math.abs( x ) < thresh ) ||
  ( Math.abs( x ) == sigma && y == 0. ) ) {
    var z = 0.;
    var w = 0.;
  } else if ( sigma == 0. ) {
    if ( x >= 0. ) {
      z = x;
      w = y;
    } else {
      z = -x;
      w = -y;
    }
  } else if ( Math.abs( x ) < thresh ) {
    z = - sigma * sigma;
    w = 0.;
  } else {
    var s = ( x >= 0. ? 1. : -1. );
    z = s * ( Math.abs( x ) - sigma ) * ( s + sigma / x );
    w = s * y;
  }
  var r = new NumberReference();
  LaPack1.dlartgp( w, z, sn, cs, r );
}
//************************************************************************
LaPack2.dlasd4 = function( n, i, d, z, delta, rho, sigma, work,
info, ioffd, ioffz, ioffdelta, ioffwork) {
  var maxit = 64;
  var dd = new Array( 3 );
  var zz = new Array( 3 );
  info.setValue( 0 );
  if ( n == 1 ) {
    sigma.setValue( Math.sqrt( d[ ioffd ] * d[ ioffd ]
      + rho * z[ ioffz ] * z[ ioffz ] ) );
    delta[ ioffdelta ] = 1.;
    work[ ioffwork ] = 1.;
    return;
  }
  if ( n == 2 ) {
    LaPack0.dlasd5( i, d, z, delta, rho, sigma, work, ioffd, ioffz,
      ioffdelta, ioffwork );
    return;
  }
  var eps = LaPack0.dlamch( 'Epsilon' );
  var rhoinv = 1. / rho;
  if ( i == n ) {
    var ii = n - 1;
    var niter = 1;
    var temp = rho / 2.;
    var temp1 = temp / ( d[ ioffd + n - 1 ]
      + Math.sqrt( d[ ioffd + n - 1 ] * d[ ioffd + n - 1 ] + temp ) );
    for ( var j = 1; j <= n; j ++ ) {
      work[ ioffwork + j - 1 ] = d[ ioffd + j - 1 ]
        + d[ ioffd + n - 1 ] + temp1;
      delta[ ioffdelta + j - 1 ] =
        ( d[ ioffd + j - 1 ] - d[ ioffd + n - 1 ] ) - temp1;
    } // 10
    var psi = 0.;
    for ( j = 1; j <= n - 2; j ++ ) {
      psi += z[ ioffz + j - 1 ] * z[ ioffz + j - 1 ]
        / ( delta[ ioffdelta + j - 1 ] * work[ ioffwork + j - 1 ] );
    } // 20
    var c = rhoinv + psi;
    var w = c + z[ ioffz + ii - 1 ] * z[ ioffz + ii - 1 ]
      / ( delta[ ioffdelta + ii - 1 ] * work[ ioffwork + ii - 1 ] )
      + z[ ioffz + n - 1 ] * z[ ioffz + n - 1 ]
      / ( delta[ ioffdelta + n - 1 ] * work[ ioffwork + n - 1 ] );
    if ( w <= 0. ) {
      temp1 =
        Math.sqrt( d[ ioffd + n - 1 ] * d[ ioffd + n - 1 ] + rho );
      temp = z[ ioffz + n - 2 ] * z[ ioffz + n - 2 ]
        / ( d[ ioffd + n - 2 ]
        + temp1 * ( d[ ioffd + n - 1 ] - d[ ioffd + n - 2 ]
        + rho / ( d[ ioffd + n - 1 ] + temp1 ) ) )
        + z[ ioffz + n - 1 ] * z[ ioffz + n - 1 ] / rho;
      if ( c <= temp ) var tau = rho;
      else {
        var delsq = ( d[ ioffd + n - 1 ] - d[ ioffd + n - 2 ] )
          * ( d[ ioffd + n - 1 ] + d[ ioffd + n - 2 ] );
        var a = - c * delsq
          + z[ ioffz + n - 2 ] * z[ ioffz + n - 2 ]
          + z[ ioffz + n - 1 ] * z[ ioffz + n - 1 ];
        var b = z[ ioffz + n - 1 ] * z[ ioffz + n - 1 ] * delsq;
        tau = ( a < 0. ?
          2. * b / ( Math.sqrt( a * a + 4. * b * c ) - a ) :
          ( a + Math.sqrt( a * a + 4. * b * c ) ) / ( 2. * c ) );
      }
    } else {
      delsq = ( d[ ioffd + n - 1 ] - d[ ioffd + n - 2 ] )
        * ( d[ ioffd + n - 1 ] + d[ ioffd + n - 2 ] );
      a = - c * delsq
        + z[ ioffz + n - 2 ] * z[ ioffz + n - 2 ]
        + z[ ioffz + n - 1 ] * z[ ioffz + n - 1 ];
      b = z[ ioffz + n - 1 ] * z[ ioffz + n - 1 ] * delsq;
      tau = ( a < 0. ?
        2. * b / ( Math.sqrt( a * a + 4. * b * c ) - a ) :
        ( a + Math.sqrt( a * a + 4. * b * c ) ) / ( 2. * c ) );
    }
    var eta = tau / ( d[ ioffd + n - 1 ]
      + Math.sqrt( d[ ioffd + n - 1 ] * d[ ioffd + n - 1 ] + tau ) );
    sigma.setValue( d[ ioffd + n - 1 ] + eta );
    for ( j = 1; j <= n; j ++ ) {
      delta[ ioffdelta + j - 1 ] =
        ( d[ ioffd + j - 1 ] - d[ ioffd + i - 1 ] ) - eta;
      work[ ioffwork + j - 1 ] =
        d[ ioffd + j - 1 ] + d[ ioffd + i - 1 ] + eta;
    } // 30
    var dpsi = 0.;
    psi = 0.;
    var erretm = 0.;
    for ( j = 1; i <= ii; j ++ ) {
      temp = z[ ioffz + j - 1 ]
        / ( delta[ ioffdelta + j - 1 ] * work[ ioffwork + j - 1 ] );
      psi += z[ ioffz + j - 1 ] * temp;
      dpsi += temp * temp;
      erretm += psi;
    } // 40
    erretm = Math.abs( erretm );
    temp = z[ ioffz + n - 1 ]
      / ( delta[ ioffdelta + n - 1 ] * work[ ioffwork + n - 1 ] );
    var phi = z[ ioffz + n - 1 ] * temp;
    var dphi = temp * temp;
    erretm += 8. * ( - phi - psi ) - phi + rhoinv
      + Math.abs( tau ) * ( dpsi + dphi );
    w = rhoinv + phi + psi;
    if ( Math.abs( w ) <= eps * erretm ) return;
    niter ++;
    var dtnsq1 = work[ ioffwork + n - 2 ]
      * delta[ ioffdelta + n - 2 ];
    var dtnsq = work[ ioffwork + n - 1 ]
      * delta[ ioffdelta + n - 1 ];
    c = w - dtnsq1 * dpsi - dtnsq * dphi;
    a = ( dtnsq + dtnsq1 ) * w - dtnsq * dtnsq1 * ( dpsi + dphi );
    b = dtnsq * dtnsq1 * w;
    if ( c < 0. ) c = Math.abs( c );
    if ( c == 0. ) eta = rho - sigma.getValue() * sigma.getValue();
    else if ( a >= 0. ) {
      eta = ( a + Math.sqrt( Math.abs( a * a - 4. * b * c ) ) )
        / ( 2. * c );
    } else {
      eta = 2. * b
        / ( a - Math.sqrt( Math.abs( a * a - 4. * b * c ) ) );
    }
    if ( w * eta > 0. ) eta = - w / ( dpsi + dphi );
    temp = eta - dtnsq;
    if ( temp > rho ) eta = rho + dtnsq;
    tau += eta;
    eta /= sigma.getValue() + Math.sqrt( eta + sigma.getValue() * sigma.getValue() );
    for ( j = 1; j <= n; j ++ ) {
      delta[ ioffdelta + j - 1 ] -= eta;
      work[ ioffwork + j - 1 ] += eta;
    } // 50
    sigma.setValue( sigma.getValue() + eta );
    dpsi = 0.;
    psi = 0.;
    erretm = 0.;
    for ( j = 1; j <= ii; j ++ ) {
      temp = z[ ioffz + j - 1 ]
        / ( work[ ioffwork + j - 1 ] * delta[ ioffdelta + j - 1 ] );
      psi += z[ ioffz + j - 1 ] * temp;
      dpsi += temp * temp;
      erretm += psi;
    }
    erretm = Math.abs( erretm );
    temp = z[ ioffz + n - 1 ]
      / ( work[ ioffwork + n - 1 ] * delta[ ioffdelta + n - 1 ] );
    phi = z[ ioffz + n - 1 ] * temp;
    dphi = temp * temp;
    erretm += 8. * ( - phi - psi ) - phi + rhoinv
      + Math.abs( tau ) * ( dpsi + dphi );
    w = rhoinv + phi + psi;
    var iter = niter + 1;
    for ( niter = iter; niter <= maxit; niter ++ ) {
      if ( Math.abs( w ) <= eps * erretm ) return;
      dtnsq1 = work[ ioffwork + n - 2 ] * delta[ ioffdelta + n - 2 ];
      dtnsq = work[ ioffwork + n - 1 ] * delta[ ioffdelta + n - 1 ];
      c = w - dtnsq1 * dpsi - dtnsq * dphi;
      a = ( dtnsq + dtnsq1 ) * w - dtnsq1 * dtnsq * ( dpsi + dphi );
      b = dtnsq1 * dtnsq * w;
      if ( a >= 0. ) {
        eta = ( a + Math.sqrt( Math.abs( a * a - 4. * b * c ) ) )
          / ( 2. * c );
      } else {
        eta = 2. * b
          / ( a - Math.sqrt( Math.abs( a * a - 4. * b * c ) ) );
      }
      if ( w * eta > 0. ) eta = - w / ( dpsi + dphi );
      temp = eta - dtnsq;
      if ( temp <= 0. ) eta /= 2.;
      tau += eta;
      eta /= sigma.getValue()
        + Math.sqrt( eta + sigma.getValue() * sigma.getValue() );
      for ( j = 1; j <= n; j ++ ) {
        delta[ ioffdelta + j - 1 ] -= eta;
        work[ ioffwork + j - 1 ] += eta;
      } // 70
      sigma.setValue( sigma.getValue() + eta );
      dpsi = 0.;
      psi = 0.;
      erretm = 0.;
      for ( j = 1; j <= ii; j ++ ) {
        temp = z[ ioffz + j - 1 ]
          / ( work[ ioffwork + j - 1 ] * delta[ ioffdelta + j - 1 ] );
        psi += z[ ioffz + j - 1 ] * temp;
        dpsi += temp * temp;
        erretm += psi;
      }
      erretm = Math.abs( erretm );
      temp = z[ ioffz + n - 1 ]
        / ( work[ ioffwork + n - 1 ] * delta[ ioffdelta + n - 1 ] );
      phi = z[ ioffz + n - 1 ] * temp;
      dphi = temp * temp;
      erretm += 8. * ( - phi - psi ) - phi + rhoinv
        + Math.abs( tau ) * ( dpsi + dphi );
      w = rhoinv + phi + psi;
    } // 90
    info.setValue( 1 );
    return;
  } else {
    niter = 1;
    var ip1 = i + 1;
    delsq = ( d[ ioffd + ip1 - 1 ] - d[ ioffd + i - 1 ] )
      * ( d[ ioffd + ip1 - 1 ] + d[ ioffd + i - 1 ] );
    var delsq2 = delsq / 2.
    temp = delsq2 / ( d[ ioffd + i - 1 ]
      + Math.sqrt( d[ ioffd + i - 1 ] * d[ ioffd + i - 1 ]
      + delsq2 ) );
    for ( j = 1; j <= n; j ++ ) {
      work[ ioffwork + j - 1 ] =
        d[ ioffd + j - 1 ] + d[ ioffd + i - 1 ] + temp;
      delta[ ioffdelta + j - 1 ] =
        ( d[ ioffd + j - 1 ] - d[ ioffd + i - 1 ] ) - temp;
    } // 100
    psi = 0.;
    for ( j = 1; j <= i - 1; j ++ ) {
      psi += z[ ioffz + j - 1 ] * z[ ioffz + j - 1 ]
        / ( work[ ioffwork + j - 1 ] * delta[ ioffdelta + j - 1 ] );
    } // 110
    phi = 0.;
    for ( j = n; j >= i + 2; j -- ) {
      phi += z[ ioffz + j - 1 ] * z[ ioffz + j - 1 ]
        / ( work[ ioffwork + j - 1 ] * delta[ ioffdelta + j - 1 ] );
    } // 120
    c = rhoinv + psi + phi;
    w = c + z[ ioffz + i - 1 ] * z[ ioffz + i - 1 ]
      / ( work[ ioffwork + i - 1 ] * delta[ ioffdelta + i - 1 ] )
      + z[ ioffz + ip1 - 1 ] * z[ ioffz + ip1 - 1 ]
      / ( work[ ioffwork + ip1 - 1 ] * delta[ ioffdelta + ip1 - 1 ] );
    if ( w > 0. ) {
      var orgati = true;
      var sg2lb = 0.;
      var sg2ub = delsq2;
      a = c * delsq + z[ ioffz + i - 1 ] * z[ ioffz + i - 1 ]
        + z[ ioffz + ip1 - 1 ] * z[ ioffz + ip1 - 1 ];
      b = z[ ioffz + i - 1 ] * z[ ioffz + i - 1 ] * delsq;
      if ( a > 0. ) {
        tau = 2. * b
          / ( a + Math.sqrt( Math.abs( a * a - 4. * b * c ) ) );
      } else {
        tau = ( a - Math.sqrt( Math.abs( a * a - 4. * b * c ) ) )
          / ( 2. * c );
      }
      eta = tau / ( d[ ioffd + i - 1 ]
        + Math.sqrt( d[ ioffd + i - 1 ] * d[ ioffd + i - 1 ] + tau ) );
    } else {
      orgati = false;
      sg2lb = - delsq2;
      sg2ub = 0.;
      a = c * delsq - z[ ioffz + i - 1 ] * z[ ioffz + i - 1 ]
        - z[ ioffz + ip1 - 1 ] * z[ ioffz + ip1 - 1 ];
      b = z[ ioffz + ip1 - 1 ] * z[ ioffz + ip1 - 1 ] * delsq;
      if ( a < 0. ) {
        tau = 2. * b
          / ( a - Math.sqrt( Math.abs( a * a + 4. * b * c ) ) );
      } else {
        tau = - ( a + Math.sqrt( Math.abs( a * a + 4. * b * c ) ) )
          / ( 2. * c );
      }
      eta = tau / ( d[ ioffd + ip1 - 1 ]
        + Math.sqrt( d[ ioffd + ip1 - 1 ] * d[ ioffd + ip1 - 1 ]
        + tau ) );
    }
    if ( orgati ) {
      ii = i;
      sigma.setValue( d[ ioffd + i - 1 ] + eta );
      for ( j = 1; j <= n; j ++ ) {
        work[ ioffwork + j - 1 ] =
          d[ ioffd + j - 1 ] + d[ ioffd + i - 1 ] + eta;
        delta[ ioffdelta + j - 1 ] =
          ( d[ ioffd + j - 1 ] - d[ ioffd + i - 1 ] ) - eta;
      } // 130
    } else {
      ii = i + 1;
      sigma.setValue( d[ ioffd + ip1 - 1 ] + eta );
      for ( j = 1; j <= n; j ++ ) {
        work[ ioffwork + j - 1 ] =
          d[ ioffd + j - 1 ] + d[ ioffd + ip1 - 1 ] + eta;
        delta[ ioffdelta + j - 1 ] =
          ( d[ ioffd + j - 1 ] - d[ ioffd + ip1 - 1 ] ) - eta;
      } // 130
    }
    var iim1 = ii - 1;
    var iip1 = ii + 1;
    dpsi = 0.;
    psi = 0.;
    erretm = 0.;
    for ( j = 1; j <= iim1; j ++ ) {
      temp = z[ ioffz + j - 1 ]
        / ( work[ ioffwork + j - 1 ] * delta[ ioffdelta + j - 1 ] );
      psi += z[ ioffz + j - 1 ] * temp;
      dpsi += temp * temp;
      erretm += psi;
    }
    erretm = Math.abs( erretm );
    dphi = 0.;
    phi = 0.;
    for ( j = n; j >= iip1; j -- ) {
      temp = z[ ioffz = j - 1 ]
        / ( work[ ioffwork + j - 1 ] * delta[ ioffdelta + j - 1 ] );
        phi += z[ ioffz + j - 1 ] * temp;
        dphi += temp * temp;
        erretm += phi;
    } // 160
    w = rhoinv + phi + psi;
    var swtch3 = false;
    if ( orgati ) {
      if ( w < 0. ) swtch3 = true;
    } else {
      if ( w > 0. ) swtch3 = true;
    }
    if ( ii == 1 || ii == n ) swtch3 = false;
    temp = z[ ioffz + ii - 1 ]
      / ( work[ ioffwork + ii - 1 ] * delta[ ioffdelta + ii - 1 ] );
    var dw = dpsi + dphi + temp * temp;
    temp = z[ ioffz + ii - 1 ] * temp;
    w += temp;
    erretm += 8. * ( phi - psi ) + 2. * rhoinv
      + 3. * Math.abs( temp ) + Math.abs( tau ) * dw;
    if ( Math.abs( w ) <= eps * erretm ) return;
    if ( w <= 0. ) sg2lb = Math.max( sg2lb, tau );
    else sg2ub = Math.min( sg2ub, tau );
    niter ++;
    if ( ! swtch3 ) {
      var dtipsq = work[ ioffwork + ip1 - 1 ]
        * delta[ ioffdelta + ip1 - 1 ];
      var dtisq = work[ ioffwork + i - 1 ]
        * delta[ ioffdelta + i - 1 ];
      c = ( orgati ?
        w - dtipsq * dw
          + delsq * Math.pow( z[ ioffz + i - 1 ] / dtisq, 2 ) :
        w - dtisq * dw
          - delsq * Math.pow( z[ ioffz + ip1 - 1 ] / dtipsq, 2 ) );
      a = ( dtipsq + dtisq ) * w - dtipsq * dtisq * dw;
      b = dtipsq * dtisq * w;
      if ( c == 0. ) {
        if ( a == 0. ) {
          a = ( orgati ?
            z[ ioffz + i - 1 ] * z[ ioffz + i - 1 ]
              + dtipsq * dtipsq * ( dpsi + dphi ) :
            z[ ioffz + ip1 - 1 ] * z[ ioffz + ip1 - 1 ]
              + dtisq * dtisq * ( dpsi + dphi ) );
        }
        eta = b / a;
      } else if ( a <= 0. ) {
        eta = ( a - Math.sqrt( Math.abs( a * a - 4. * b * c ) ) )
          / ( 2. * c );
      } else {
        eta = 2. * b
          / ( a + Math.sqrt( Math.abs( a * a - 4. * b * c ) ) );
      }
    } else {
      var dtiim =
        work[ ioffwork + iim1 - 1 ] * delta[ ioffdelta + iim1 - 1 ];
      var dtiip =
        work[ ioffwork + iip1 - 1 ] * delta[ ioffdelta + iip1 - 1 ];
      temp = rhoinv + psi + phi;
      if ( orgati ) {
        temp1 = z[ ioffz + iim1 - 1 ] / dtiim;
        temp1 = temp1 * temp1;
        c = ( temp - dtiip * ( dpsi + dphi ) )
          - ( d[ ioffd + iim1 - 1 ] - d[ ioffd + iip1 - 1 ] )
          * ( d[ ioffd + iim1 - 1 ] + d[ ioffd + iip1 - 1 ] ) * temp1;
        zz[ 0 ] = z[ ioffz + iim1 - 1 ] * z[ ioffz + iim1 - 1 ];
        zz[ 2 ] = ( dpsi < temp1 ?  dtiip * dtiip * dphi :
          dtiip * dtiip * ( ( dpsi - temp1 ) + dphi ) );
      } else {
        temp1 = z[ ioffz + iip1 - 1 ] / dtiip;
        temp1 = temp1 * temp1;
        c = ( temp - dtiim * ( dpsi + dphi ) )
          - ( d[ ioffd + iip1 - 1 ] - d[ ioffd + iim1 - 1 ] )
          * ( d[ ioffd + iim1 - 1 ] + d[ ioffd + iip1 - 1 ] ) * temp1;
        zz[ 0 ] = ( dphi < temp1 ?  dtiim * dtiim * dpsi :
          dtiim * dtiim * ( dpsi + ( dphi - temp1 ) ) );
        zz[ 2 ] = z[ ioffz + iip1 - 1 ] * z[ ioffz + iip1 - 1 ];
      }
      zz[ 1 ] = z[ ioffz + ii - 1 ] * z[ ioffz + ii - 1 ];
      dd[ 0 ] = dtiim;
      dd[ 1 ] =
        delta[ ioffdelta + ii - 1 ] * work[ ioffwork + ii - 1 ];
      dd[ 2 ] = dtiip;
      var tauref = new NumberReference( eta );
      LaPack1.dlaed6( niter, orgati, c, dd, zz, w, tauref, info,
        0, 0 );
      eta = tauref.getValue();
      if ( info.getValue() != 0 ) return;
    }
    if ( w * eta >= 0. ) eta = - w / dw;
    if ( orgati ) {
      temp1 = work[ ioffwork + i - 1 ] * delta[ ioffdelta + i - 1 ];
      temp = eta - temp1;
    } else {
      temp1 =
        work[ ioffwork + ip1 - 1 ] * delta[ ioffdelta + ip1 - 1 ];
      temp = eta - temp1;
    }
    if ( temp > sg2ub || temp < sg2lb ) {
      eta = ( w < 0. ? ( sg2ub - tau ) / 2. : ( sg2lb - tau ) / 2. );
    }
    tau += eta;
    eta /= sigma.getValue() + Math.sqrt( sigma.getValue() * sigma.getValue() + eta );
    var prew = w;
    sigma.setValue( sigma.getValue() + eta );
    for ( j = 1; j <= n; j ++ ) {
      work[ ioffwork + j - 1 ] += eta;
      delta[ ioffdelta + j - 1 ] -= eta;
    } // 170
    dpsi = 0.;
    psi = 0.;
    erretm = 0.;
    for ( j = 1; j <= iim1; j ++ ) {
      temp = z[ ioffz + j - 1 ]
        / ( work[ ioffwork + j - 1 ] * delta[ ioffdelta + j - 1 ] );
      psi += z[ ioffz + j - 1 ] * temp;
      dpsi += temp * temp;
      erretm += psi;
    } // 180
    erretm = Math.abs( erretm );
    dphi = 0.;
    psi = 0.
    for ( j = n; j >= iip1; j -- ) {
      temp = z[ ioffz + j - 1 ]
        / ( work[ ioffwork + j - 1 ] * delta[ ioffdelta + j - 1 ] );
      phi += z[ ioffz + j - 1 ] * temp;
      dphi += temp * temp;
      erretm += phi;
    } // 190
    temp = z[ ioffz + ii - 1 ]
      / ( work[ ioffwork + ii - 1 ] * delta[ ioffdelta + ii - 1 ] );
    dw = dpsi + dphi + temp * temp;
    temp *= z[ ioffz + ii - 1 ];
    w = rhoinv + phi + psi + temp;
    erretm += 8. * ( phi - psi ) + 2. * rhoinv + 3. * Math.abs( temp )
      + Math.abs( tau ) * dw;
    if ( w <= 0. ) sg2lb = Math.max( sg2lb, tau );
    else sg2ub = Math.min( sg2ub, tau );
    var swtch = false;
    if ( orgati ) {
      if ( - w > Math.abs( prew ) / 10. ) swtch = true;
    } else {
      if ( w > Math.abs( prew ) / 10. ) swtch = true;
    }
    iter = niter + 1;
    for ( niter = iter; niter <= maxit; niter ++ ) {
      if ( Math.abs( w ) <= eps * erretm ) return;
      if ( ! swtch3 ) {
        dtipsq = work[ ioffwork + ip1 - 1 ]
          * delta[ ioffdelta + ip1 - 1 ];
        dtisq = work[ ioffwork + i - 1 ] * delta[ ioffdelta + i - 1 ];
        if ( ! swtch ) {
          c = ( orgati ?
            w - dtipsq * dw
              + delsq * Math.pow( z[ ioffz + i - 1 ] / dtisq, 2 ) :
            w - dtisq * dw
              - delsq * Math.pow( z[ ioffz + ip1 - 1 ] / dtipsq, 2 ) );
        } else {
          temp = z[ ioffz + ii - 1 ]
            / ( work[ ioffwork + ii - 1 ]
            * delta[ ioffdelta + ii - 1 ] );
          if ( orgati ) dpsi += temp * temp;
          else dphi += temp * temp;
          c = w - dtisq * dpsi - dtipsq * dphi;
        }
        a = ( dtipsq + dtisq ) * w - dtipsq * dtisq * dw;
        b = dtipsq * dtisq * w;
        if ( c == 0. ) {
          if ( a == 0. ) {
            if ( ! swtch ) {
              a = ( orgati ?
                z[ ioffz + i - 1 ] * z[ ioffz + i - 1 ]
                  + dtipsq * dtipsq * ( dpsi + dphi ) :
                z[ ioffz + ip1 - 1 ] * z[ ioffz + ip1 - 1 ]
                  + dtisq * dtisq * ( dpsi + dphi ) );
            } else {
              a = dtisq * dtisq * dpsi + dtipsq * dtipsq * dphi;
            }
          }
          eta = b / a;
        } else if ( a <= 0. ) {
          eta = ( a - Math.sqrt( Math.abs( a * a - 4. * b * c ) ) )
            / ( 2. * c );
        } else {
          eta = 2. * b
            / ( a + Math.sqrt( Math.abs( a * a - 4. * b * c ) ) );
        }
      } else {
        dtiim =
          work[ ioffwork + iim1 - 1 ] * delta[ ioffdelta + iim1 - 1 ];
        dtiip =
          work[ ioffwork + iip1 - 1 ] * delta[ ioffdelta + iip1 - 1 ];
        temp = rhoinv + psi + phi;
        if ( swtch ) {
          c = temp - dtiim * dpsi - dtiip * dphi;
          zz[ 0 ] = dtiim * dtiim * dpsi;
          zz[ 2 ] = dtiip * dtiip * dphi;
        } else {
          if ( orgati ) {
            temp1 = z[ ioffz + iim1 - 1 ] / dtiim;
            temp1 = temp1 * temp1;
            var temp2 =
              ( d[ ioffd + iim1 - 1 ] - d[ ioffd + iip1 - 1 ] )
              * ( d[ ioffd + iim1 - 1 ] + d[ ioffd + iip1 - 1 ] )
              * temp1;
            c = temp - dtiip * ( dpsi + dphi ) - temp2;
            zz[ 0 ] = z[ ioffz + iim1 - 1 ] * z[ ioffz + iim1 - 1 ];
            zz[ 2 ] = ( dpsi < temp1 ?  dtiip * dtiip * dphi :
              dtiip * dtiip * ( ( dpsi - temp1 ) + dphi ) );
          } else {
            temp1 = z[ ioffz + iip1 - 1 ] / dtiip;
            temp1 = temp1 * temp1;
            temp2 = ( d[ ioffd + iip1 - 1 ] - d[ ioffd + iim1 - 1 ] )
              * ( d[ ioffd + iim1 - 1 ] + d[ ioffd + iip1 - 1 ] )
              * temp1;
            c = temp - dtiim * ( dpsi + dphi )  - temp2;
            zz[ 0 ] = ( dphi < temp1 ?  dtiim * dtiim * dpsi :
              dtiim * dtiim * ( dpsi + ( dphi - temp1 ) ) );
            zz[ 2 ] = z[ ioffz + iip1 - 1 ] * z[ ioffz + iip1 - 1 ];
          }
        }
        dd[ 0 ] = dtiim;
        dd[ 1 ] =
          delta[ ioffdelta + ii - 1 ] * work[ ioffwork + ii - 1 ];
        dd[ 2 ] = dtiip;
        tauref.setValue( eta );
        LaPack1.dlaed6( niter, orgati, c, dd, zz, w, tauref, info,
          0, 0 );
        eta = tauref.getValue();
        if ( info.getValue() != 0 ) return;
      }
      if ( w * eta >= 0. ) eta = - w / dw;
      if ( orgati ) {
        temp1 = work[ ioffwork + i - 1 ] * delta[ ioffdelta + i - 1 ];
        temp = eta - temp1;
      } else {
        temp1 =
          work[ ioffwork + ip1 - 1 ] * delta[ ioffdelta + ip1 - 1 ];
        temp = eta - temp1;
      }
      if ( temp > sg2ub || temp < sg2lb ) {
        eta = ( w < 0. ? ( sg2ub - tau ) / 2. : ( sg2lb - tau ) / 2. );
      }
      tau += eta;
      eta /= sigma.getValue()
        + Math.sqrt( sigma.getValue() * sigma.getValue() + eta );
      sigma.setValue( sigma.getValue() + eta );
      for ( j = 1; j <= n; j ++ ) {
        work[ ioffwork + j - 1 ] += eta;
        delta[ ioffdelta + j - 1 ] -= eta;
      } // 2000
      prew = w;
      dpsi = 0.;
      psi = 0.;
      erretm = 0.;
      for ( j = 1; j <= iim1; j ++ ) {
        temp = z[ ioffz + j - 1 ]
          / ( work[ ioffwork + j - 1 ] * delta[ ioffdelta + j - 1 ] );
        psi += z[ ioffz + j - 1 ] * temp;
        dpsi += temp * temp;
        erretm += psi;
      } // 210
      erretm = Math.abs( erretm );
      dphi = 0.;
      phi = 0.
      for ( j = n; j >= iip1; j -- ) {
        temp = z[ ioffz + j - 1 ]
          / ( work[ ioffwork + j - 1 ] * delta[ ioffdelta + j - 1 ] );
        phi += z[ ioffz + j - 1 ] * temp;
        dphi += temp * temp;
        erretm += phi;
      } // 220
      temp = z[ ioffz + ii - 1 ]
        / ( work[ ioffwork + ii - 1 ] * delta[ ioffdelta + ii - 1 ] );
      dw = dpsi + dphi + temp * temp;
      temp *= z[ ioffz + ii - 1 ];
      w = rhoinv + phi + psi + temp;
      erretm += 8. * ( phi - psi ) + 2. * rhoinv
        + 3. * Math.abs( temp ) + Math.abs( tau ) * dw;
      if ( w * prew > 0. && Math.abs( w ) > Math.abs( prew ) / 10. ) {
        swtch = ! swtch;
      }
      if ( w <= 0. ) sg2lb = Math.max( sg2lb, tau );
      else sg2ub = Math.min( sg2ub, tau );
    } // 230
    info.setValue( 1 );
  }
}
//************************************************************************
LaPack2.dlasq3 = function( i0, n0, z, pp, dmin, sigma, desig,
qmax, nfail, iter, ndiv, ieee, ttype, dmin1, dmin2, dn, dn1, dn2, g, tau,
ioffz ) {
  throw new Error("not tested: complicated input");
  var cbias = 1.5;
  var n0in = n0.getValue();
  var eps = LaPack0.dlamch( 'Precision' );
  var tol = eps * 100.;
  var tol2 = tol * tol;
  while ( true ) { // 10
    if ( n0.getValue() < i0 ) return;
    var goto20 = false;
    var goto30 = false;
    var goto40 = false;
    if ( n0.getValue() == i0 ) goto20 = true;
    if ( ! goto20 ) {
      var nn = 4 * n0.getValue() + pp.getValue();
      if ( n0.getValue() == ( i0 + 1 ) ) goto40 = true;
      if ( ! goto40 ) {
        if ( z[ ioffz + nn - 6 ] >
        tol2 * ( sigma.getValue() + z[ ioffz + nn - 4 ] ) &&
        z[ ioffz + nn - 2 * pp.getValue() - 5 ] >
        tol2 * z[ ioffz + nn - 8 ] ) {
          goto30 = true;
        }
      }
    } // 20
    if ( ! goto30 && ! goto40 ) {
      z[ ioffz + 4 * n0.getValue() - 4 ] =
        z[ ioffz + 4 * n0.getValue() + pp.getValue() - 4 ] + sigma.getValue();
      n0.setValue( n0.getValue() - 1 );
      continue;
    } // 30
    if ( ! goto40 ) {
      if ( z[ ioffz + nn - 10 ] > tol2 * sigma.getValue() &&
      z[ ioffz + nn - 2 * pp.getValue() - 9 ] >
      tol2 * z[ ioffz + nn - 12 ] ) {
        break;
      }
    } // 40
    if ( z[ ioffz + nn - 4 ] > z[ ioffz + nn - 8 ] ) {
      var s = z[ ioffz + nn - 4 ];
      z[ ioffz + nn - 4 ] = z[ ioffz + nn - 8 ];
      z[ ioffz + nn - 8 ] = s;
    }
    if ( z[ ioffz + nn - 6 ] > z[ ioffz + nn - 4 ] * tol2 ) {
      var t =
        0.5 * ( ( z[ ioffz + nn - 8 ] - z[ ioffz+ nn - 4 ] )
        + z[ ioffz + nn - 6 ] );
      s = z[ ioffz + nn - 4 ] * ( z[ ioffz + nn - 6 ] / t );
      if ( s <= t ) {
        s = z[ ioffz + nn - 4 ] * ( z[ ioffz + nn - 6 ]
          / ( t * ( 1. + Math.sqrt( 1. + s / t ) ) ) );
      } else {
        s = z[ ioffz + nn - 4 ] * ( z[ ioffz + nn - 6 ]
          / ( t + Math.sqrt( t ) * Math.sqrt( t + s ) ) );
      }
      t = z[ ioffz + nn - 8 ] + ( s + z[ ioffz + nn - 6 ] );
      z[ ioffz + nn - 4 ] *= z[ ioffz + nn - 8 ] / t;
      z[ ioffz + nn - 8 ] = t;
    }
    z[ ioffz + 4 * n0.getValue() - 8 ] = z[ ioffz + nn - 8 ] * sigma.getValue();
    z[ ioffz + 4 * n0.getValue() - 4 ] = z[ ioffz + nn - 4 ] * sigma.getValue();
    n0.setValue( n0.getValue() - 2 );
  } // 50
  if ( pp.getValue() == 2 ) pp.getValue() = 0;
  if ( dmin.getValue() <= 0. || n0.getValue() < n0in ) {
    if ( cbias * z[ ioffz + 4 * i0 + pp.getValue() - 4 ] <
    z[ ioffz + 4 * n0.getValue() + pp.getValue() - 4 ] ) {
      var ipn4 = 4 * ( i0 + n0.getValue() );
      for ( var j4 = 4 * i0; j4 <= 2 * ( i0 + n0.getValue() - 1 );
      j4 += 4 ) {
        var temp = z[ ioffz + j4 - 4 ];
        z[ ioffz + j4 - 4 ] = z[ ioffz + ipn4 - j4 - 4 ];
        z[ ioffz + ipn4 - j4 - 4 ] = temp;
        temp = z[ ioffz + j4 - 3 ];
        z[ ioffz + j4 - 3 ] = z[ ioffz + ipn4 - j4 - 3 ];
        z[ ioffz + ipn4 - j4 - 3 ] = temp;
        temp = z[ ioffz + j4 - 2 ];
        z[ ioffz + j4 - 2 ] = z[ ioffz + ipn4 - j4 - 6 ];
        z[ ioffz + ipn4 - j4 - 6 ] = temp;
        temp = z[ ioffz + j4 - 1 ];
        z[ ioffz + j4 - 1 ] = z[ ioffz + ipn4 - j4 - 5 ];
        z[ ioffz + ipn4 - j4 - 5 ] = temp;
      } // 60
      if ( n0.getValue() - i0 <= 4 ) {
        z[ ioffz + 4 * n0.getValue() + pp.getValue() - 2 ] =
          z[ ioffz + 4 * i0 + pp.getValue() - 2 ];
        z[ ioffz + 4 * n0.getValue() - pp.getValue() - 1 ] =
          z[ ioffz + 4 * i0 - pp.getValue() - 1 ];
      }
      dmin2.setValue( Math.min( dmin2.getValue(),
        z[ ioffz + 4 * n0.getValue() + pp.getValue() - 2 ] ) );
      z[ ioffz + 4 * n0.getValue() + pp.getValue() - 2 ] = Math.min(
        Math.min( z[ ioffz + 4 * n0.getValue() + pp.getValue() - 2 ],
        z[ ioffz + 4 * i0 + pp.getValue() - 2 ] ),
        z[ ioffz + 4 * i0 + pp.getValue() + 2 ] );
      z[ ioffz + 4 * n0.getValue() - pp.getValue() - 1 ] = Math.min(
        Math.min( z[ ioffz + 4 * n0.getValue() - pp.getValue() - 1 ],
        z[ ioffz + 4 * i0 - pp.getValue() - 1 ] ),
        z[ ioffz + 4 * i0 - pp.getValue() + 3 ] );
      qmax.setValue( Math.max( qmax.getValue(),
        Math.max( z[ ioffz + 4 * i0 + pp.getValue() - 4 ],
        z[ ioffz + 4 * i0 + pp.getValue() ] ) ) );
      dmin.setValue( - 0. );
    }
  }
  LaPack0.dlasq4( i0, n0.getValue(), z, pp.getValue(), n0in, dmin.getValue(),
    dmin1.getValue(), dmin2.getValue(), dn.getValue(), dn1.getValue(), dn2.getValue(), tau,
    ttype, g, ioffz );
  var goto80 = false;
  var goto90 = false;
  while ( true ) { // 70
    LaPack0.dlasq5( i0, n0.getValue(), z, pp.getValue(), tau.getValue(), dmin, dmin1,
      dmin2, dn, dn1, dn2, ieee, ioffz );
    ndiv.setValue( ndiv.getValue() + ( n0.getValue() - i0 + 2 ) );
    iter.setValue( iter.getValue() + 1 );
    if ( dmin.getValue() >= 0. && dmin1.getValue() > 0. ) {
      goto90 = true;
      break;
    } else if ( dmin.getValue() < 0. && dmin1.getValue() > 0. &&
    z[ ioffz + 4 * ( n0.getValue() - 1 ) - pp.getValue() - 1 ] <
    tol * ( sigma.getValue() + dn1.getValue() ) &&
    Math.abs( dn.getValue() ) < tol * sigma.getValue() ) {
      z[ ioffz + 4 * ( n0.getValue() - 1 ) - pp.getValue() + 1 ] = 0.;
      dmin.setValue( 0. );
      goto90 = true;
      break;
    } else if ( dmin.getValue() < 0. ) {
      nfail.setValue( nfail.getValue() + 1 );
      if ( ttype.getValue() < -22 )  tau.getValue() = 0.;
      else if ( dmin1.getValue() > 0. ) {
        tau.setValue( ( tau.getValue() + dmin.getValue() ) * ( 1. - 2. * eps ) );
        ttype.setValue( ttype.getValue() - 11 );
      } else {
        tau.setValue( 0.25 * tau.getValue() );
        ttype.setValue( ttype.getValue() -12 );
      }
      continue;
    } else if ( isNaN( dmin.getValue() ) ) {
      if ( tau.getValue() == 0. ) {
        goto80 = true;
        break;
      } else {
        tau.setValue( 0. );
        continue;
      }
    } else {
      goto80 = true;
      break;
    }
  } // 80
  if ( ! goto90 ) {
    LaPack1.dlasq6( i0, n0.getValue(), z, pp.getValue(), dmin, dmin1, dmin2, dn,
      dn1, dn2, ioffz );
    ndiv.setValue( ndiv.getValue() + ( n0.getValue() - i0 + 2 ) );
    iter.setValue( iter.getValue() + 1 );
    tau.setValue( 0. );
  } // 90
  if ( tau.getValue() < sigma.getValue() ) {
    desig.setValue( desig.getValue() + tau.getValue() );
    t = sigma.getValue() + desig.getValue();
    desig.setValue( desig.getValue() - ( t - sigma.getValue() ) );
  } else {
    t = sigma.getValue() + tau.getValue();
    desig.setValue( sigma.getValue() - ( t - tau.getValue() ) + desig.getValue() );
  }
  sigma.setValue( t );
}
/* old version:
LaPack2.dlasq3 = function( n, q, e, qq, ee, sup, sigma, kend,
off, iphase, iconv, eps, tol2, small2, ioffq, ioffe, ioffqq, ioffee ) {
  var npp = 32;
  var ipp = 5;
  var iflmax = 2;
  var icnt = 0;
  var tau = 0.;
  var dm = sup.getValue();
  var tolx = sigma.getValue() * tol2;
  var tolz = Math.max( small2, sigma.getValue() ) * tol2;
  var maxit = 100 * n.getValue();
  var ic = 2;
  if ( n.getValue() > 3 ) {
    if ( iphase.getValue() == 1 ) {
      for ( var i = 1; i <= n.getValue() - 2; i ++ ) {
        if ( q[ ioffq + i - 1 ] > q[ ioffq + i ] ) ic ++;
        if ( e[ ioffe + i - 1 ] > e[ ioffe + i ] ) ic ++;
      }
      if ( q[ ioffq + n.getValue() - 2 ] > q[ ioffq + n.getValue() - 1 ] ) ic ++;
      if ( ic < n.getValue() ) {
        Blas1.dcopy( n.getValue(), q, 1, qq, -1, ioffq, ioffqq );
        Blas1.dcopy( n.getValue() - 1, e, 1, ee, -1, ioffe, ioffee );
        if ( kend.getValue() != 0 ) kend.getValue() = n.getValue() - kend.getValue() + 1;
        iphase.setValue( 2 );
      }
    } else {
      for ( i = 1; i <= n.getValue() - 2; i ++ ) {
        if ( qq[ ioffqq + i - 1 ] > qq[ ioffqq + i ] ) ic ++;
        if ( ee[ ioffee + i - 1 ] > ee[ ioffee + i ] ) ic ++;
      }
      if ( qq[ ioffqq + n.getValue() - 2 ] > qq[ ioffqq + n.getValue() - 1 ] ) {
        ic ++;
      }
      if ( ic < n.getValue() ) {
        Blas1.dcopy( n.getValue(), qq, 1, q, -1, ioffqq, ioffq );
        Blas1.dcopy( n.getValue() - 1, ee, 1, e, -1, ioffee, ioffe );
        if ( kend.getValue() != 0 ) kend.getValue() = n.getValue() - kend.getValue() + 1;
        iphase.setValue( 1 );
      }
    }
  }
  var goto80 = false;
  var goto130 = false;
  var goto140 = false;
  var goto180 = false;
  if ( iconv.getValue() == -3 ) {
    if ( iphase.getValue() == 1 ) goto180 = true;
    else goto80 = true;
  } else if ( iphase.getValue() == 2 ) goto130 = true;
  while ( true  ) { // 30
    if ( ! goto80 && ! goto130 && ! goto140 && ! goto180 ) {
      if ( kend.getValue() == 0 || sup.getValue() == 0. ) tau = 0.;
      else if ( icnt > 0 && dm <= tolz ) tau = 0.;
      else {
        var ip = Math.max( ipp, n.getValue() / npp );
        var n2 = 2 * ip + 1;
        if ( n2 >= n.getValue() ) {
          var n1 = 1;
          n2 = n.getValue();
        } else if ( kend.getValue() + ip > n.getValue() ) n1 = n.getValue() - 2 * ip;
        else if ( kend.getValue() - ip < 1 ) n1 = 1;
        else n1 = kend.getValue() - ip;
        var tauref = new NumberReference( tau );
        LaPack0.dlasq4( n2, q, e, tauref, sup, ioffq + n1 - 1,
          ioffe + n1 - 1 );
        tau = tauref.getValue();
      }
    }
    if ( ! goto130 && ! goto140 && ! goto180 ) {
      var ifl = 0;
      while ( true ) { // 40, ends at 130
        var goto120 = false;
        if ( ! goto80 ) {
          icnt ++;
          if ( icnt > maxit ) {
            sup.setValue( -1. );
            return;
          }
          if ( tau == 0. ) { // dqd algorithm in ping
            var d = q[ ioffq ];
            dm = d;
            var ke = 0;
            for ( i = 1; i <- n.getValue() - 3; i ++ ) { // 50
              qq[ ioffqq + i - 1 ] = d + e[ ioffe + i - 1 ];
              d = ( d / qq[ ioffqq + i - 1 ] ) * q[ ioffq + i ];
              if ( dm > d ) {
                dm = d;
                ke = i;
              }
            }
            ke ++;
            var k2end = ke;
            qq[ ioffqq + n.getValue() - 3 ] = d + e[ ioffe + n.getValue() - 3 ];
            d = ( d / qq[ ioffqq + n.getValue() - 3] )
              * q[ ioffq + n.getValue() - 2 ];
            if ( dm > d ) {
              dm = d;
              ke = n.getValue() - 1;
            }
            var k1end = ke;
            qq[ ioffqq + n.getValue() - 2 ] = d + e[ ioffe + n.getValue() - 2 ];
            d = ( d / qq[ ioffqq + n.getValue() - 2 ] )
              * q[ ioffq + n.getValue() - 1 ];
            if ( dm > d ) {
              dm = d;
              ke = n.getValue();
            }
            qq[ ioffqq + n.getValue() - 1 ] = d;
          } else { // dqds algorithm in ping
            d = q[ ioffq ] - tau;
            dm = d;
            ke = 0;
            if ( d < 0. ) goto120 = true;
            if ( ! goto120 ) {
              for ( i = 1; i <= n.getValue() - 3; i ++ ) { // 60
                qq[ ioffqq + i - 1 ] = d + e[ ioffe + i - 1 ];
                d = ( d / qq[ ioffqq + i - 1 ] ) * q[ ioffq + i ]
                  - tau;
                if ( dm > d ) {
                  dm = d;
                  ke = i;
                  if ( d < 0. ) {
                    goto120 = true;
                    break;
                  }
                }
              }
            }
            if ( ! goto120 ) {
              ke ++;
              k2end = ke;
              qq[ ioffqq + n.getValue() - 3 ] =
                d + e[ ioffe + n.getValue() - 3 ];
              d = ( d / qq[ ioffqq + n.getValue() - 3 ] )
                * q[ ioffq + n.getValue() - 2 ] - tau;
              if ( dm < d ) {
                dm = d;
                ke = n.getValue() - 1;
                if ( d < 0. ) goto120 = true;
              }
            }
            if ( ! goto120 ) {
              k1end = ke;
              qq[ ioffqq + n.getValue() - 2 ] =
                d + e[ ioffe + n.getValue() - 2 ];
              d = ( d / qq[ ioffqq + n.getValue() - 2 ] )
                * q[ ioffq + n.getValue() - 1 ] - tau;
              if ( dm < d ) {
                dm = d;
                ke = n.getValue();
              }
              qq[ ioffqq + n.getValue() - 1 ] = d;
            }
          }
          if ( ! goto120 ) {
            if ( Math.abs( qq[ ioffqq + n.getValue() - 1 ] )
            <= sigma.getValue() * tol2 ) {
              qq[ ioffqq + n.getValue() - 1 ] = 0.;
              dm = 0.;
              ke = n.getValue();
            }
            if ( qq[ ioffqq + n.getValue() - 1 ] < 0. ) goto120 = true;
          }
          if ( ! goto120 ) {
            for ( i = 1; i <= n.getValue() - 1; i ++ ) { // 70
              ee[ ioffee + i - 1 ] = ( ee[ ioffee + i - 1 ]
                / qq[ ioffqq + i - 1 ] ) * q[ ioffq + i ];
            }
            sigma.setValue( sigma.getValue() + tau );
            iphase.setValue( 2 );
          }
        } // 80
        goto80 = false;
        if ( ! goto120 ) {
          tolx = sigma.getValue() * tol2;
          var toly = sigma.getValue() * eps;
          tolz = Math.max( sigma.getValue(), small2 ) * tol2;
          while ( true ) { // 90
            if ( n.getValue() <= 2 ) return;
            var ldef = false;
            if ( ee[ ioffee + n.getValue() - 2 ] <= tolz ) ldef = true;
            else if ( sigma.getValue() < 0. ) {
              if ( ee[ ioffee + n.getValue() - 2 ]
              <= eps * ( sigma.getValue() + qq[ ioffqq + n.getValue() - 1 ] ) ) {
                if ( ee[ ioffee + n.getValue() - 2 ]
                * ( qq[ ioffqq + n.getValue() - 1 ]
                / ( qq[ ioffqq + n.getValue() - 1 ] + sigma.getValue() ) ) <=
                tol2 * ( qq[ ioffqq + n.getValue() - 1 ] + sigma.getValue() ) ) {
                  ldef = true;
                }
              }
            } else {
              if ( ee[ ioffee + n.getValue() - 2 ] <=
              qq[ ioffqq + n.getValue() - 1 ] * tol2 ) {
                ldef = true;
              }
            }
            if ( ldef ) {
              q[ ioffq + n.getValue() - 1 ] =
                qq[ ioffqq + n.getValue() - 1 ] + sigma.getValue();
              n.setValue( n.getValue() - 1 );
              iconv.setValue( iconv.getValue() + 1 );
              continue;
            }
            ldef = false;
            if ( ee[ ioffee + n.getValue() - 3 ] <= tolz ) ldef = true;
            else if ( sigma.getValue() > 0. ) {
              t1 = sigma.getValue() + ee[ ioffee + n.getValue() - 2 ]
               * ( sigma.getValue() / ( sigma.getValue()
               + qq[ ioffqq + n.getValue() - 1 ] ) );
              if ( ee[ ioffee + n.getValue() -3 ]
              * ( t1 / ( qq[ ioffqq + n.getValue() - 2 ] + t1 ) ) <= toly )
              {
                if ( ee[ ioffee + n.getValue() - 3 ]
                * ( qq[ ioffqq + n.getValue() - 2 ]
                / ( qq[ ioffqq + n.getValue() - 2 ] + t1 ) ) <= tolx ) {
                  ldef = true;
                }
              }
            } else {
              if ( ee[ ioffee + n.getValue() - 3 ] <=
              ( qq[ ioffqq + n.getValue() - 1 ]
              / ( qq[ ioffqq + n.getValue() - 1 ]
              + ee[ ioffee + n.getValue() - 2 ]
              + qq[ ioffqq + n.getValue() - 2 ] ) )
              * qq[ ioffqq + n.getValue() - 2 ] * tol2 ) {
                ldef = true;
              }
            }
            if ( ldef ) {
              var qemax = Math.max( Math.max(
                qq[ ioffqq + n.getValue() - 1 ] ,
                qq[ ioffqq + n.getValue() - 2 ] ),
                ee[ ioffee + n.getValue() - 2 ] );
              if ( qemax != 0. ) {
                if ( qemax == qq[ ioffqq + n.getValue() - 2 ] ) {
                  var xx = 0.5 * ( qq[ ioffqq + n.getValue() - 1 ]
                    + qq[ ioffqq + n.getValue() - 2 ]
                    + ee[ ioffee + n.getValue() - 2 ]
                    + qemax * Math.sqrt( Math.pow( (
                    qq[ ioffqq + n.getValue() - 1 ]
                    - qq[ ioffqq + n.getValue() - 2 ]
                    + ee[ ioffee + n.getValue() - 2 ] ) / qemax , 2 )
                    + 4. * ee[ ioffee + n.getValue() - 2 ] / qemax ) );
                } else if ( qemax == qq[ ioffqq + n.getValue() - 1 ] ) {
                  xx = 0.5 * ( qq[ ioffqq + n.getValue() - 1 ]
                    + qq[ ioffqq + n.getValue() - 2 ]
                    + ee[ ioffee + n.getValue() - 2 ]
                    + qemax * Math.sqrt( Math.pow(
                    ( qq[ ioffqq + n.getValue() - 2 ]
                    - qq[ ioffqq + n.getValue() - 1 ]
                    + ee[ ioffee + n.getValue() - 2 ] ) / qemax ,
                    2 ) + 4. * ee[ ioffee + n.getValue() - 2 ] / qemax ) );
                } else {
                  xx = 0.5 * ( qq[ ioffqq + n.getValue() - 1 ]
                    + qq[ ioffqq + n.getValue() - 2 ]
                    + ee[ ioffee + n.getValue() - 2 ]
                    + qemax * Math.sqrt( Math.pow(
                    ( qq[ ioffqq + n.getValue() - 1 ]
                    - qq[ ioffqq + n.getValue() - 2 ]
                    + ee[ ioffee + n.getValue() - 2 ] ) / qemax , 2 )
                    + 4. * ee[ ioffee + n.getValue() - 2 ] / qemax ) );
                }
                yy = ( Math.max( qq[ ioffqq + n.getValue() - 1 ] ,
                  qq[ ioffqq + n.getValue() - 2 ] ) / xx )
                  * Math.min( qq[ ioffqq + n.getValue() - 1 ] ,
                  qq[ ioffqq + n.getValue() - 2 ] );
              } else {
                xx = 0.;
                yy = 0.;
              }
              q[ ioffq + n.getValue() - 2 ] = sigma.getValue() + xx;
              q[ ioffq + n.getValue() - 1 ] = yy + sigma.getValue();
              n.setValue( n.getValue() -2 );
              iconv.setValue( iconv.getValue() + 2 );
            } else break;
          }
          if ( iconv.getValue() == 0 ) { // update bound before pong
            kend.setValue( ke );
            sup.setValue( Math.min( dm, sup.getValue() - tau ) );
          } else if ( iconv.getValue() > 0 ) {
            sup.setValue( Math.min( Math.min( qq[ ioffqq + n.getValue() - 1 ],
              qq[ ioffqq + n.getValue() - 2 ] ),
              qq[ ioffqq + n.getValue() - 3 ] ) );
            sup.setValue( Math.min( Math.min( qq[ ioffqq ],
              qq[ ioffqq + 1 ] ),
              Math.min( qq[ ioffqq + 2 ] , sup.getValue() ) ) );
            if ( iconv.getValue() == 1 ) kend.getValue() = k1end;
            else if ( iconv.getValue() == 2 ) kend.getValue() = k2end;
            else kend.setValue( n.getValue() );
            icnt = 0;
            maxit = 100 * n.getValue();
          }
          var lsplit = false;
          for ( var ks = n.getValue() - 3; ks >= 3; ks -- ) { // 100
            if ( ee[ ioffee + ks - 1 ] <= toly ) {
              if ( ee[ ioffee + ks - 1 ]
              * ( Math.min( qq[ ioffqq + ks ] ,
              qq[ ioffqq + ks - 1 ] ) / ( Math.min( qq[ ioffqq + ks ],
              qq [ ioffqq + ks - 1 ] ) + sigma.getValue() ) ) <= tolx ) {
                lsplit = true;
                break;
              }
            }
          }
          if ( ! lsplit ) {
            ks = 2;
            if ( ee[ ioffee + 1 ] <= tolz ) {
              lsplit = true;
            } else if ( sigma.getValue() > 0. ) {
              t1 = sigma.getValue() + ee[ ioffee ] * ( sigma.getValue()
                / ( sigma.getValue() + qq[ ioffqq ] ) );
              if ( ee[ ioffe + 1 ] * ( t1 / ( qq[ ioffqq ] + t1 ) )
              <= toly ) {
                if ( ee[ ioffe + 1 ] * ( qq[ ioffqq ]
                / ( qq[ ioffqq ] + t1 ) ) <= tolx ) {
                  lsplit = true;
                }
              }
            } else {
              if ( ee[ ioffe + 1 ] <= ( qq[ ioffqq ]
              / ( qq[ ioffqq ] + ee[ ioffee ]
              + qq[ ioffqq + 1 ] ) ) * qq[ ioffqq + 1 ] * tol2 ) {
                lsplit = true;
              }
            }
          }
          if ( ! lsplit ) {
            ks = 1;
            if ( ee[ ioffee ] <= tolz ) lsplit = true;
            else if ( sigma.getValue() > 0. ) {
              if ( ee[ ioffee ] <=
              eps * ( sigma.getValue() + qq[ ioffqq ] ) ) {
                if ( ee[ ioffee ] * ( qq[ ioffqq ]
                / ( qq[ ioffqq ] + sigma.getValue() ) )
                <= tol2 * ( qq[ ioffqq ] + sigma.getValue() ) ) {
                  lsplit = true;
                }
              }
            } else {
              if ( ee[ ioffee ] <= qq[ ioffqq ] * tol2 ) lsplit = true;
            }
          } // 110
          if ( lsplit ) {
            sup.setValue( Math.min( Math.min( qq[ ioffqq + n.getValue() - 1],
              qq[ ioffqq + n.getValue() - 2 ] ),
              qq[ ioffqq + n.getValue() - 3 ] ) );
            var isp = - ( off.getValue() + 1 );
            off.setValue( off.getValue() + ks );
            n.setValue( n.getValue() - ks );
            kend.setValue( Math.max( 1, kend.getValue() - ks ) );
            e[ ioffe + ks - 1 ] = sigma.getValue();
            ee[ ioffee + ks - 1 ] = isp;
            iconv.setValue( 0 );
            return;
          }
          if ( tau == 0. && dm <= tolz && kend.getValue() != n.getValue()
          && iconv.getValue() == 0 && icnt > 0 ) {
            Blas1.dcopy( n.getValue() - ke, e, 1, qq, 1, ioffe + ke - 1,
              ioffqq + ke - 1 );
            qq[ ioffqq + n.getValue() - 1 ] = 0.;
            Blas1.dcopy( n.getValue() - ke, q, 1, ee, 1, ioffq + ke,
              ioffee + ke - 1 );
            sup.setValue( 0. );
          }
          iconv.setValue( 0 );
          break;
        } // 120
        ifl ++;
        sup.setValue( tau );
        if ( sup.getValue() <= tolz || ifl >= iflmax ) tau = 0.;
        else {
          tau = Math.max( tau + d, 0. );
          if ( tau <= tolz ) tau = 0.;
        }
      } // 130, began at 40
    }
    if ( ! goto140 && ! goto180 ) {
      ifl = 0;
      if ( kend.getValue() == 0 && sup.getValue() == 0. ) tau = 0.;
      else if ( icnt > 0 && dm <= tolz ) tau = 0.;
      else {
        ip = Math.max( ipp, n.getValue() / npp );
        n2 = 2 * ip + 1;
        if ( n2 >= n.getValue() ) {
          n1 = 1;
          n2 = n.getValue();
        } else if ( kend.getValue() + ip > n.getValue() ) n1 = n.getValue() - 2 * ip;
        else if ( kend.getValue() - ip < 1 ) n1 = 1;
        else n1 = kend.getValue() - ip;
        tauref.setValue( tau );
        LaPack0.dlasq4( n2, qq, ee, tauref, sup, ioffqq + n1 - 1,
          ioffee + n1 - 1 );
        tau = tauref.getValue()
      }
    }
    goto140 = false;
    var goto220 = false;
    if ( ! goto180 ) { // 140
      icnt ++;
      if ( icnt > maxit ) {
        sup.setValue( - sup.getValue() );
        return;
      }
      if ( tau == 0. ) {
        d = qq[ ioffqq ];
        dm = d;
        ke = 0;
        for ( i = 1; i <= n.getValue() - 3; i ++ ) { // 150
          q[ ioffq + i - 1 ] = d + ee[ ioffee + i - 1 ];
          d = ( d / q[ ioffq + i - 1 ] ) * qq[ ioffqq + i ];
          if ( dm > d ) {
            dm = d;
            ke = i;
          }
        }
        ke ++;
        k2end = ke;
        q[ ioffq + n.getValue() - 3 ] = d + ee[ ioffee + n.getValue() - 3 ];
        d = ( d / q[ ioffq + n.getValue() - 3 ] )
          * qq[ ioffqq + n.getValue() - 2 ];
        if ( dm > d ) {
          dm = d;
          ke = n.getValue() - 1;
        }
        k1end = ke;
        q[ ioffq + n.getValue() - 2 ] = d + ee[ ioffee + n.getValue() - 2 ];
        d = ( d / q[ ioffq + n.getValue() - 2 ] )
          * qq[ ioffqq + n.getValue() - 1 ];
        if ( dm > d ) {
          dm = d;
          ke = n.getValue();
        }
        q[ ioffq + n.getValue() - 1 ] = d;
      } else {
        d = qq[ ioffqq ] - tau;
        dm = d;
        ke = 0;
        if ( d < 0. ) break;
        for ( i = 1; i <= n.getValue() - 3; i ++ ) { // 160
          q[ ioffq + i - 1 ] = d + ee[ ioffee + i - 1 ];
          d = ( d / q[ ioffq + i - 1 ] ) * qq[ ioffqq + i ] - tau;
          if ( dm > d ) {
            dm = d;
            ke = i;
            if ( d < 0. ) {
              goto220 = true;
              break;
            }
          }
        }
        if ( ! goto220 ) {
          ke ++;
          k2end = ke;
          q[ ioffq + n.getValue() - 3 ] = d + ee[ ioffee + n.getValue() - 3 ];
          d = ( d / q[ ioffq + n.getValue() - 3 ] )
            * qq[ ioffqq + n.getValue() - 2 ] - tau;
          if ( dm > d ) {
            dm = d;
            ke = n.getValue() - 1;
            if ( d < 0. ) goto220 = true;
          }
        }
        if ( ! goto220 ) {
          k1end = ke;
          q[ ioffq + n.getValue() - 2 ] = d + ee[ ioffee + n.getValue() - 2 ];
          d = ( d / q[ ioffq + n.getValue() - 2 ] )
            * qq[ ioffqq + n.getValue() - 1 ] - tau;
          if ( dm > d ) {
            dm = d;
            ke = n.getValue();
          }
          q[ ioffq + n.getValue() - 1 ] = d;
        }
      }
      if ( ! goto220 ) {
        if ( Math.abs( q[ ioffq + n.getValue() - 1 ] )
        <= sigma.getValue() * tol2 ) {
          q[ ioffq + n.getValue() - 1 ] = 0.;
          dm = 0.;
          ke = n.getValue();
        }
        if ( q[ ioffq + n.getValue() - 1 ] < 0. ) goto220 = true;
      }
      if ( ! goto220 ) {
        for ( i = 1; i <= n.getValue() - 1; i ++ ) { // 170
          e[ ioffe + i - 1 ] = ( ee[ ioffee + i - 1 ]
            / qq[ ioffqq + i - 1 ] ) * qq[ ioffqq + i ];
        }
        sigma.setValue( sigma.getValue() + tau );
      }
    } // 180
    if ( ! goto220 ) {
      iphase.setValue( 1 );
      tolx = sigma.getValue() * tol2;
      toly = sigma.getValue() * eps;
      while ( true ) { // 190
        if ( n.getValue() <= 2 ) return;
        ldef = false;
        if ( e[ ioffe + n.getValue() - 2 ] <= tolz ) ldef = true;
        else if ( sigma.getValue() > 0. ) {
          if ( e[ ioffe + n.getValue() - 2 ] <= eps * ( sigma.getValue()
          + q[ ioffq + n.getValue() - 1 ] ) ) {
            if ( e[ ioffe + n.getValue() - 2 ] * ( q[ ioffq + n.getValue() - 1 ]
            / ( q[ ioffq + n.getValue() - 1 ] + sigma.getValue() ) )
            <= tol2 * ( q[ ioffq + n.getValue() - 1 ] + sigma.getValue() ) ) {
              ldef = true;
            }
          }
        } else {
          if ( e[ ioffe + n.getValue() - 2 ]
          <= q[ ioffq + n.getValue() - 1 ] * tol2 ) {
            ldef = true;
          }
        }
        if ( ldef ) {
          q[ ioffq + n.getValue() - 1 ] += sigma.getValue();
          n.setValue( n.getValue() - 1 );
          iconv.setValue( iconv.getValue() + 1 );
          continue;
        }
        ldef = false;
        if ( e[ ioffe + n.getValue() - 3 ] <= tolz ) ldef = true;
        else if ( sigma.getValue() > 0. ) {
          var t1 =
            sigma.getValue() + ( e[ ioffe + n.getValue() - 2 ] * ( sigma.getValue()
            / ( sigma.getValue() + q[ ioffq + n.getValue() - 1 ] ) ) );
          if ( e[ ioffe + n.getValue() - 3 ] * ( t1
          / ( q[ ioffq + n.getValue() - 2 ] + t1 ) ) <= toly ) {
            if ( e[ ioffe + n.getValue() - 3 ] * ( q[ ioffq + n.getValue() - 2 ]
            / ( q[ ioffq + n.getValue() - 2 ] + t1 ) ) <= tolx ) {
              ldef = true;
            }
          }
        } else {
          if ( e[ ioffe + n.getValue() - 3 ] <= ( q[ ioffq + n.getValue() - 1 ]
          / ( q[ ioffq + n.getValue() - 1 ] + ee[ ioffee + n.getValue() - 2 ]
          + q[ ioffq + n.getValue() - 2 ] ) * q[ ioffq + n.getValue() - 2 ] )
          * tol2 ) {
            ldef = true;
          }
        }
        if ( ldef ) {
          qemax = Math.max( Math.max( q[ ioffq + n.getValue() - 1 ],
            q[ ioffq + n.getValue() - 2 ] ), e[ ioffe + n.getValue() - 2 ] );
          if ( qemax != 0. ) {
            if ( qemax == q[ ioffq + n.getValue() - 2 ] ) {
              xx = 0.5 * ( q[ ioffq + n.getValue() - 1 ]
                + q[ ioffq + n.getValue() - 2 ]
                + e[ ioffe + n.getValue() - 2 ]
                + qemax * Math.sqrt( Math.pow(
                  ( q[ ioffq + n.getValue() - 1 ] - q[ ioffq + n.getValue() - 2 ]
                  + e[ ioffe + n.getValue() - 2 ] ) / qemax, 2 )
                + 4. * e[ ioffe + n.getValue() - 2 ] / qemax ) );
            } else if ( qemax == q[ ioffq + n.getValue() - 1 ] ) {
              xx = 0.5 * ( q[ ioffq + n.getValue() - 1 ]
                + q[ ioffq + n.getValue() - 2 ]
                + e[ ioffe + n.getValue() - 2 ]
                + qemax * Math.sqrt( Math.pow(
                  ( q[ ioffq + n.getValue() - 2 ] - q[ ioffq + n.getValue() - 1 ]
                  + e[ ioffe + n.getValue() - 2 ] ) / qemax, 2 )
                + 4. * e[ ioffe + n.getValue() - 2 ] / qemax ) );
            } else {
              xx = 0.5 * ( q[ ioffq + n.getValue() - 1 ]
                + q[ ioffq + n.getValue() - 2 ]
                + e[ ioffe + n.getValue() - 2 ]
                + qemax * Math.sqrt( Math.pow(
                  ( q[ ioffq + n.getValue() - 1 ] - q[ ioffq + n.getValue() - 2 ]
                  + e[ ioffe + n.getValue() - 2 ] ) / qemax, 2 )
                + 4. * e[ ioffe + n.getValue() - 2 ] / qemax ) );
            }
            var yy = ( Math.max( q[ ioffq + n.getValue() - 1 ],
              q[ ioffq + n.getValue() - 2 ] ) / xx )
              * Math.min( q[ ioffq + n.getValue() - 1 ] ,
              q[ ioffq + n.getValue() - 2 ] );
          } else {
            xx = 0.
            yy = 0.;
          }
          q[ ioffq + n.getValue() - 2 ] = sigma.getValue() + xx;
          q[ ioffq + n.getValue() - 1 ] = yy + sigma.getValue();
          n.setValue( n.getValue() - 2 );
          iconv.setValue( iconv.getValue() + 2 );
        } else break;
      }
      if ( iconv.getValue() == 0 ) {
        kend.setValue( ke );
        sup.setValue( Math.min( dm, sup.getValue() - tau ) );
      } else if ( iconv.getValue() > 0 ) {
        sup.setValue( Math.min( Math.min( q[ ioffq + n.getValue() - 1 ],
          q[ ioffq + n.getValue() - 2 ] ), q[ ioffq + n.getValue() - 3 ] ) );
        sup.setValue( Math.min( Math.min( q[ ioffq ], q[ ioffq + 1 ] ),
          Math.min( q[ ioffq + 2 ], sup.getValue() ) ) );
        if ( iconv.getValue() == 1 ) kend.getValue() = k1end;
        else if ( iconv.getValue() == 2 ) kend.setValue( k2end );
        else kend.setValue( n.getValue() );
        icnt = 0;
        maxit = 100 * n.getValue();
      }
      lsplit = false;
      for ( ks = n.getValue() - 3; ks >= 3; ks -- ) { // 200
        if ( e[ ioffe + ks - 1 ] <= toly ) {
          if ( e[ ioffe + ks - 1 ] * ( Math.min( q[ ioffq + ks ],
          q[ ioffq + ks - 1 ] )
          / ( Math.min( q[ ioffq + ks ], q[ ioffq + ks - 1 ] )
          + sigma.getValue() ) ) <= tolx ) {
            lsplit = true;
          }
        }
      }
      if ( ! lsplit ) {
        ks = 2;
        if ( e[ ioffe + 1 ] <= tolz ) lsplit = true;
        else if ( sigma.getValue() > 0. ) {
          t1 = sigma.getValue() + e[ ioffe ] * ( sigma.getValue()
            / ( sigma.getValue() + q[ ioffq ] ) );
          if ( e[ ioffe + 1 ] * ( t1 / ( q[ ioffq ] + t1 ) ) <= toly )
          {
            if ( e[ ioffe + 1 ] * ( q[ ioffq ]
            / ( q[ ioffq ] + t1 ) ) <= tolx ) {
              lsplit = true;
            }
          }
        } else {
          if ( e[ ioffe + 1 ] <= ( q[ ioffq ]
          / ( q[ ioffq ] + e[ ioffe ] + q[ ioffq + 1 ] ) )
          * q[ ioffq + 1 ] * tol2 ) {
            lsplit = true;
          }
        }
      }
      if ( ! lsplit ) {
        ks = 1;
        if ( e[ ioffe ] <= tolz ) lsplit = true;
        else if ( sigma.getValue() > 0. ) {
          if ( e[ ioffe ] <= eps * ( sigma.getValue() + q[ ioffq ] ) ) {
            if ( e[ ioffe ] * ( q[ ioffq ]
            / ( q[ ioffq ] + sigma.getValue() ) )
            <= tol2 * ( q[ ioffq ] + sigma.getValue() ) ) {
              lsplit = true;
            }
          }
        } else {
          if ( e[ ioffe ] <= q[ ioffq ] * tol2 ) lsplit = true;
        }
      } // 210
      if ( lsplit ) {
        sup.setValue( Math.min( Math.min( q[ ioffq + n.getValue() - 1 ],
          q[ ioffq + n.getValue() - 2 ] ), q[ ioffq + n.getValue() - 3 ] ) );
        isp = off.getValue() + 1;
        off.setValue( off.getValue() + ks );
        kend.setValue( Math.max( 1, kend.getValue() - ks ) );
        n.setValue( n.getValue() - ks );
        e[ ioffe + ks - 1 ] = sigma.getValue();
        ee[ ioffee + ks - 1 ] = isp;
        iconv.setValue( 0 );
        return;
      }
      if ( tau == 0. && dm <= tolz && kend.getValue() != n.getValue()
      && iconv.getValue() == 0 && icnt > 0 ) {
        Blas1.dcopy( n.getValue() - ke, ee, 1, q, 1, ioffee + ke - 1,
          ioffq + ke - 1 );
        q[ ioffq + n.getValue() - 1 ] = 0.;
        Blas1.dcopy( n.getValue() - ke, qq, 1, e, 1, ioffqq + ke,
          ioffe + ke - 1 );
        sup.setValue( 0. );
      }
      iconv.setValue( 0. );
      continue;
    } // 220
    ifl ++;
    sup.setValue( tau );
    if ( sup.getValue() <= tolz || ifl >= iflmax ) {
      tau = 0.;
      goto140 = true;
    } else {
      tau = Math.max( tau + d, 0. )
      if ( tau <= tolz ) tau = 0.;
      goto140 = true;
    }
  }
}
*/
//************************************************************************
LaPack2.dlatrd = function( uplo, n, nb, A, lda, e, tau, W, ldw,
ioffa, ioffe, iofftau, ioffw ) {
  if ( n <= 0 ) return;
  var alpha = new NumberReference();
  var tauref = new NumberReference();
  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    for ( var i = n; i >= n - nb + 1; i -- ) {
      var iw = i - n + nb;
      if ( i < n ) {
        Blas2.dgemv( 'No transpose', i, n - i, -1., A, lda, W, ldw, 1.,
          A, 1, ioffa + i * lda, ioffw + i - 1 + iw * ldw,
          ioffa + ( i - 1 ) * lda );
        Blas2.dgemv( 'No transpose', i, n - i, -1., W, ldw, A, lda, 1.,
          A, 1, ioffw + iw * ldw, ioffa + i - 1 + i * lda,
          ioffa + ( i - 1 ) * lda );
      }
      if ( i > 1 ) {
        alpha.setValue( A[ ioffa + i - 2 + ( i - 1 ) * lda ] );
        tauref.setValue( tau[ iofftau + i - 2 ] );
        LaPack1.dlarfg( i - 1, alpha, A, 1, tauref,
          ioffa + ( i - 1 ) * lda );
        A[ ioffa + i - 2 + ( i - 1 ) * lda ] = alpha.getValue();
        tau[ iofftau + i - 2 ] = tauref.getValue();
        e[ ioffe + i - 2 ] = A[ ioffa + i - 2 + ( i - 1 ) * lda ];
        A[ ioffa + i - 2 + ( i - 1 ) * lda ] = 1.;
        Blas2.dsymv( 'Upper', i - 1, 1., A, lda, A, 1, 0., W, 1,
          ioffa, ioffa + ( i - 1 ) * lda, ioffw + ( iw - 1 ) * ldw );
        if ( i < n ) {
          Blas2.dgemv( 'Transpose', i - 1, n - i, 1., W, ldw, A, 1, 0.,
            W, 1, ioffw + iw * ldw, ioffa + ( i - 1 ) * lda,
            ioffw + i + ( iw - 1 ) * ldw );
          Blas2.dgemv( 'No transpose', i - 1, n - i, -1., A, lda,
            W, 1, 1., W, 1, ioffa + i * lda,
            ioffw + i + ( iw - 1 ) * ldw, ioffw + ( iw - 1 ) * ldw );
          Blas2.dgemv( 'Transpose', i - 1, n - i, 1., A, lda, A, 1,
            0., W, 1, ioffa + i * lda, ioffa + ( i - 1 ) * lda,
            ioffw + i + ( iw - 1 ) * ldw );
          Blas2.dgemv( 'No transpose', i - 1, n - i, -1., W, ldw,
            W, 1, 1., W, 1, ioffw + iw * ldw,
            ioffw + i + ( iw - 1 ) * ldw, ioffw + ( iw - 1 ) * ldw );
        }
        Blas1.dscal( i - 1, tau[ iofftau + i - 2 ], W, 1,
          ioffw + ( iw - 1 ) * ldw );
        alpha.setValue( - 0.5 * tau[ iofftau + i - 2 ]
          * Blas1.ddot( i - 1, W, 1, A, 1, ioffw + ( iw - 1 ) * ldw,
          ioffa + ( i - 1 ) * lda ) );
        Blas1.daxpy( i - 1, alpha.getValue(), A, 1, W, 1,
          ioffa + ( i - 1 ) * lda, ioffw + ( iw - 1 ) * ldw );
      }
    }
  } else {
    for ( i = 1; i <= nb; i ++ ) {
      Blas2.dgemv( 'No transpose', n - i + 1, i - 1, -1., A, lda,
        W, ldw, 1., A, 1, ioffa + i - 1, ioffw + i - 1,
        ioffa + i - 1 + ( i - 1 ) * lda );
      Blas2.dgemv( 'No transpose', n - i + 1, i - 1, -1., W, ldw,
        A, lda, 1., A, 1, ioffw + i - 1, ioffa + i - 1,
        ioffa + i - 1 + ( i - 1 ) * lda );
      if ( i < n ) {
        alpha.setValue( A[ ioffa + i + ( i - 1 ) * lda ] );
        tauref.setValue( tau[ iofftau + i - 1 ] );
        LaPack1.dlarfg( n - i, alpha, A, 1, tauref,
          ioffa + Math.min( i + 2, n ) - 1 + ( i - 1 ) * lda );
        A[ ioffa + i + ( i - 1 ) * lda ] = alpha.getValue();
        tau[ iofftau + i - 1 ] = tauref.getValue();
        e[ ioffe + i - 1 ] = A[ ioffa + i + ( i - 1 ) * lda ];
        A[ ioffa + i + ( i - 1 ) * lda ] = 1.;
        Blas2.dsymv( 'Lower', n - i, 1., A, lda, A, 1, 0., W, 1,
          ioffa + i + i * lda, ioffa + i + ( i - 1 ) * lda,
          ioffw + i + ( i - 1 ) * ldw );
        Blas2.dgemv( 'Transpose', n - i, i - 1, 1., W, ldw, A, 1, 0.,
          W, 1, ioffw + i, ioffa + i + ( i - 1 ) * lda,
          ioffw + ( i - 1 ) * ldw );
        Blas2.dgemv( 'No transpose', n - i, i - 1, -1., A, lda,
          W, 1, 1., W, 1, ioffa + i, ioffw + ( i - 1 ) * ldw,
          ioffw + i + ( i - 1 ) * ldw );
        Blas2.dgemv( 'Transpose', n - i, i - 1, 1., A, lda, A, 1,
          0., W, 1, ioffa + i, ioffa + i + ( i - 1 ) * lda,
          ioffw + ( i - 1 ) * ldw );
        Blas2.dgemv( 'No transpose', n - i, i - 1, -1., W, ldw,
          W, 1, 1., W, 1, ioffw + i, ioffw + ( i - 1 ) * ldw,
          ioffw + i + ( i - 1 ) * ldw );
        Blas1.dscal( n - i, tau[ iofftau + i - 1 ], W, 1,
          ioffw + i + ( i - 1 ) * ldw );
        alpha.setValue( - 0.5 * tau[ iofftau + i - 1 ]
          * Blas1.ddot( n - i, W, 1, A, 1, ioffw + i + ( i - 1 ) * ldw,
          ioffa + i + ( i - 1 ) * lda ) );
        Blas1.daxpy( n - i, alpha.getValue(), A, 1, W, 1,
          ioffa + i + ( i - 1 ) * lda, ioffw + i + ( i - 1 ) * ldw );
      }
    }
  }
}
LaPack2.zlatrd = function( uplo, n, nb, A, lda, e, tau, W,
ldw ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dlatrz = function( m, n, l, A, lda, tau, work, ioffa,
iofftau, ioffwork ) {
  if ( m == 0 ) return;
  else if ( m == n ) {
    for ( var i = 1; i <= n; i ++ ) tau[ iofftau + i - 1 ] = 0.;
    return;
  }
  for ( i = m; i >= 1; i -- ) {
    var alpha =
      new NumberReference( A[ ioffa + i - 1 + ( i - 1 ) * lda ] );
    var tauref =
      new NumberReference( tau[ iofftau + i - 1 ] );
    LaPack1.dlarfg( l + 1, alpha, A, lda, tauref,
      ioffa + i - 1 + ( n - l ) * lda );
    A[ ioffa + i - 1 + ( i - 1 ) * lda ] = alpha.getValue();
    tau[ iofftau + i - 1 ] = tauref.getValue();
    LaPack0.dlarz( 'Right', i - 1, n - i + 1, l, A, lda,
      tau[ iofftau + i - 1 ], A, lda, work,
      ioffa + i - 1 + ( n - l ) * lda, ioffa + ( i - 1 ) * lda,
      ioffwork );
  } // 20
}
//************************************************************************
LaPack2.dopmtr = function( side, uplo, trans, m, n, AP, tau, C,
ldc, work, info, ioffap, iofftau, ioffc, ioffwork) {
  throw new Error("not tested: input from packed matrix");
  info.setValue( 0 );
  var left = ( side.charAt(0).toUpperCase() == 'L' );
  var notran = ( trans.charAt(0).toUpperCase() == 'N' );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var nq = ( left ? m : n );
  if ( ! left && side.charAt(0).toUpperCase() != 'R' ) {
    info.setValue( -1 );
  } else if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -2 );
  } else if ( ! notran && trans.charAt(0).toUpperCase() != 'T' ) {
    info.setValue( -3 );
  } else if ( m < 0 ) info.setValue( -4 );
  else if ( n < 0 ) info.setValue( -5 );
  else if ( ldc < Math.max( 1, m ) ) info.setValue( -9 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dopmtr', - info.getValue() );
    return;
  }
  if ( m == 0 || n == 0 ) return;
  if ( upper ) {
    var forwrd = ( left && notran ) || ( ! left && ! notran );
    if ( forwrd) {
      var i1 = 1;
      var i2 = nq - 1;
      var i3 = 1;
      var ii = 2;
    } else {
      i1 = nq - 1;
      i2 = 1;
      i3 = -1;
      ii = ( nq * ( nq + 1 ) ) / 2 - 1;
    }
    if ( left ) var ni = n;
    else var mi = m;
    for ( var i = i1;
    ( i3 == 1 && i <= i2 || i3 == -1 && i >= i2 ); i += i3 ) {
      if ( left ) mi = i;
      else ni = i;
      var aii = AP[ ioffap + ii - 1 ];
      AP[ ioffap + ii - 1 ] = 1.;
      LaPack1.dlarf( side, mi, ni, AP, 1, tau[ iofftau + i - 1 ], C,
        ldc, work, ioffap + ii - i, ioffc, ioffwork );
      AP[ ioffap + ii - 1 ] = aii;
      if ( forwrd ) ii += i + 2;
      else ii -= i + 1;
    }
  } else {
    forwrd = ( left && ! notran ) || ( ! left && notran );
    if ( forwrd) {
      i1 = 1;
      i2 = nq - 1;
      i3 = 1;
      ii = 2;
    } else {
      i1 = nq - 1;
      i2 = 1;
      i3 = -1;
      ii = ( nq * ( nq + 1 ) ) / 2 - 1;
    }
    if ( left ) {
      ni = n;
      var jc = 1;
    } else {
      mi = m;
      var ic = 1;
    }
    for ( i = i1; ( i3 == 1 && i <= i2 || i3 == -1 && i >= i2 );
    i += i3 ) {
      aii= AP[ ioffap + ii - 1 ];
      AP[ ioffap + ii - 1 ] = 1.;
      if ( left ) {
        mi = m - i;
        ic = i + 1;
      } else {
        ni = n - i;
        jc = i + 1;
      }
      LaPack1.dlarf( side, mi, ni, AP, 1, tau[ iofftau + i - 1 ], C,
        ldc, work, ioffap + ii - 1, ioffc + ic - 1 + ( jc - 1 ) * ldc,
        ioffwork );
      AP[ ioffap + ii - 1 ] = aii;
      if ( forwrd ) ii += nq - i + 1;
      else ii -= nq - i + 2;
    }
  }
}
//************************************************************************
LaPack2.dorbdb = function( trans, signs, m, p, q, X11, ldx11,
X12, ldx12, X21, ldx21, X22, ldx22, theta, phi, taup1, taup2, tauq1,
tauq2, work, lwork, info, ioffx11, ioffx12, ioffx21, ioffx22, iofftheta,
ioffphi, iofftaup1, iofftaup2, iofftauq1, iofftauq2, ioffwork ) {
  throw new Error("not tested");
  info.setValue( 0 );
  var colmajor = trans.charAt(0).toUpperCase() != 'T';
  if ( signs.charAt(0).toUpperCase() != 'O' ) {
    var z1 = 1.;
    var z2 = 1.;
    var z3 = 1.;
    var z4 = 1.;
  } else {
    z1 = 1.;
    z2 = -1.;
    z3 = 1.;
    z4 = -1.;
  }
  var lquery = ( lwork == -1 );
  if ( m < 0 ) info.setValue( -3 );
  else if ( p < 0 || p > m ) info.setValue( -4 );
  else if ( q < 0 || q > p || q > m - p || q > m - q ) info.setValue( -5 );
  else if ( colmajor && ldx11 < Math.max( 1, p ) ) info.setValue( -7 );
  else if ( ! colmajor && ldx11 < Math.max( 1, q ) ) info.setValue( -7 );
  else if ( colmajor && ldx12 < Math.max( 1, p ) ) info.setValue( -9 );
  else if ( ! colmajor && ldx12 < Math.max( 1, m - q ) ) {
    info.setValue( -9 );
  } else if ( colmajor && ldx21 < Math.max( 1, m - p ) ) {
    info.setValue( -11 );
  } else if ( ! colmajor && ldx21 < Math.max( 1, q ) ) {
    info.setValue( -11 );
  } else if ( colmajor && ldx22 < Math.max( 1, m - p ) ) {
    info.setValue( -13 );
  } else if ( ! colmajor && ldx22 < Math.max( 1, m - q ) ) {
    info.setValue( -13 );
  }
  if ( info.getValue() == 0 ) {
    var lworkopt = m - q;
    var lworkmin = m - q;
    work[ ioffwork ] = lworkopt;
    if ( lwork < lworkmin  && ! lquery ) info.setValue( -21 );
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dorbdb', -info.getValue() );
    return;
  } else if ( lquery ) return;
  if ( colmajor ) {
    for ( i = 1; i <= q; i ++ ) {
      if ( i == 1 ) {
        Blas1.dscal( p - i + 1, z1, X11, 1,
          ioffx11 + i - 1 + ( i - 1 ) * ioffx11 );
      } else {
        Blas1.dscal( p - i + 1,
          z1 * Math.cos( phi[ ioffphi + i - 2 ] ), X11, 1,
          ioffx11 + i - 1 + ( i - 1 ) * ioffx11 );
        Blas1.daxpy( p - i + 1,
          - z1 * z3 * z4 * Math.sin( phi[ ioffphi + i - 2 ] ),
          X12, 1, X11, 1, ioffx12 + i - 1 + ( i - 2 ) * ldx12,
          ioffx11 + i - 1 + ( i - 1 ) * ldx11 );
      }
      if ( i == 1 ) {
        Blas1.dscal( m - p - i + 1, z2, X21, 1,
          ioffx21 + i - 1 + ( i - 1 ) * ioffx21 );
      } else {
        Blas1.dscal( m - p - i + 1,
          z2 * Math.cos( phi[ ioffphi + i - 2 ] ), X21, 1,
          ioffx21 + i - 1 + ( i - 1 ) * ioffx21 );
        Blas1.daxpy( m - p - i + 1,
          - z2 * z3 * z4 * Math.sin( phi[ ioffphi + i - 2 ] ),
          X22, 1, X21, 1, ioffx22 + i - 1 + ( i - 2 ) * ldx22,
          ioffx21 + i - 1 + ( i - 1 ) * ldx21 );
      }
      theta[ iofftheta + i - 1 ] = Math.atan2(
        Blas1.dnrm2( m - p - i + 1, X21, 1,
          ioffx21 + i - 1 + ( i - 1 ) * ldx21 ),
        Blas1.dnrm2( p - i + 1, X11, 1,
          ioffx11 + i - 1 + ( i - 1 ) * ldx11 ) );
      var alpha = new NumberReference(
        X11[ ioffx11 + i - 1 + ( i - 1 ) * ldx11 ] );
      var tau =
        new NumberReference( taup1[ iofftaup1 + i - 1 ] );
      LaPack1.dlarfgp( p - i + 1, alpha, X11, 1, tau,
        ioffx11 + i + ( i - 1 ) * ldx11 );
      X11[ ioffx11 + i - 1 + ( i - 1 ) * ldx11 ] = alpha.getValue();
      taup1[ iofftaup1 + i - 1 ] = tau.getValue();
      X11[ ioffx11 + i - 1 + ( i - 1 ) * ldx11 ] = 1.;
      alpha.setValue( X21[ ioffx21 + i - 1 + ( i - 1 ) * ldx21 ] );
      tau.setValue( taup2[ iofftaup2 + i - 1 ] );
      LaPack1.dlarfgp( m - p - i + 1, alpha, X21, 1, tau,
        ioffx21 + i + ( i - 1 ) * ldx21 );
      X21[ ioffx21 + i - 1 + ( i - 1 ) * ldx21 ] = alpha.getValue();
      taup2[ iofftaup2 + i - 1 ] = tau.getValue();
      X21[ ioffx21 + i - 1 + ( i - 1 ) * ldx21 ] = 1.;
      LaPack1.dlarf( 'L', p - i + 1, q - i, X11, 1,
        taup1[ iofftaup1 + i - 1 ], X11, ldx11, work,
        ioffx11 + i - 1 + ( i - 1 ) * ldx11,
        ioffx11 + i - 1 + i * ldx11, ioffwork );
      LaPack1.dlarf( 'L', p - i + 1, m - q - i + 1, X11, 1,
        taup1[ iofftaup1 + i - 1 ], X12, ldx12, work,
        ioffx11 + i - 1 + ( i - 1 ) * ldx11,
        ioffx12 + i - 1 + ( i - 1 ) * ldx12, ioffwork );
      LaPack1.dlarf( 'L', m - p - i + 1, q - i, X21, 1,
        taup2[ iofftaup2 + i - 1 ], X21, ldx21, work,
        ioffx21 + i - 1 + ( i - 1 ) * ldx21,
        ioffx21 + i - 1 + i * ldx21, ioffwork );
      LaPack1.dlarf( 'L', m - p - i + 1, m - q - i + 1, X21, 1,
        taup2[ iofftaup2 + i - 1 ], X22, ldx22, work,
        ioffx21 + i - 1 + ( i - 1 ) * ldx21,
        ioffx22 + i - 1 + ( i - 1 ) * ldx22, ioffwork );
      if ( i < q ) {
        Blas1.dscal( q - i,
          - z1 * z3 * Math.sin( theta[ iofftheta + i - 1 ] ),
          X11, ldx11, ioffx11 + i - 1 + i * ldx11 );
        Blas1.daxpy( q - i,
          z2 * z3 * Math.cos( theta[ iofftheta + i - 1 ] ),
          X21, ldx21, X11, ldx11, ioffx21 + i - 1 + i * ldx21,
          ioffx11 + i - 1 + i * ldx11 );
      }
      Blas1.dscal( m - q - i + 1,
        - z1 * z4 * Math.sin( theta[ iofftheta + i - 1 ] ),
        X12, ldx12, ioffx12 + i - 1 + ( i - 1 ) * ldx12 );
      Blas1.daxpy( m - q - i + 1,
        z2 * z4 * Math.cos( theta[ iofftheta + i - 1 ] ),
        X22, ldx22, X12, ldx12, ioffx22 + i - 1 + ( i - 1 ) * ldx22,
        ioffx12 + i - 1 + ( i - 1 ) * ldx12 );
      if ( i < q ) {
        phi[ ioffphi + i - 1 ] = Math.atan2(
          Blas1.dnrm2( q - i, X11, ldx11,
            ioffx11 + i - 1 + i * ldx11 ),
          Blas1.dnrm2( m - q - i + 1, X12, ldx12,
            ioffx12 + i - 1 + ( i - 1 ) * ldx12 ) );
      }
      if ( i < q ) {
        alpha.setValue( X11[ ioffx11 + i - 1 + i * ldx11 ] );
        tau.setValue( tauq1[ iofftauq1 + i - 1 ] );
        LaPack1.dlarfgp( q - i, alpha, X11, ldx11, tau,
          ioffx11 + i - 1 + ( i + 1 ) * ldx11 );
        X11[ ioffx11 + i - 1 + i * ldx11 ] = alpha.getValue();
        tauq1[ iofftauq1 + i - 1 ] = tau.getValue();
        X11[ ioffx11 + i - 1 + i * ldx11 ] = 1.;
      }
      alpha.setValue( X12[ ioffx12 + i - 1 + ( i - 1 ) * ldx12 ] );
      tau.setValue( tauq2[ iofftauq2 + i - 1 ] );
      LaPack1.dlarfgp( m - q - i + 1, alpha, X12, ldx12, tau,
        ioffx12 + i - 1 + i * ldx12 );
      X12[ ioffx12 + i - 1 + ( i - 1 ) * ldx12 ] = alpha.getValue();
      tauq2[ iofftauq2 + i - 1 ] = tau.getValue();
      X12[ ioffx12 + i - 1 + ( i - 1 ) * ldx12 ] = 1.;
      if ( i < q ) {
        LaPack1.dlarf( 'R', p - i, q - i, X11, ldx11,
          tauq1[ iofftauq1 + i - 1 ], X11, ldx11, work,
          ioffx11 + i - 1 + i * ldx11, ioffx11 + i + i * ldx11,
          ioffwork );
        LaPack1.dlarf( 'R', m - p - i, q - i, X11, ldx11,
          tauq1[ iofftauq1 + i - 1 ], X21, ldx21, work,
          ioffx11 + i - 1 + i * ldx11, ioffx21 + i + i * ldx21,
          ioffwork );
      }
      LaPack1.dlarf( 'R', p - i, m - q - i + 1, X12, ldx12,
        tauq2[ iofftauq2 + i - 1 ], X12, ldx12, work,
        ioffx12 + i - 1 + ( i - 1 ) * ldx12,
        ioffx12 + i + ( i - 1 ) * ldx12, ioffwork );
      LaPack1.dlarf( 'R', m - p - i, m - q - i + 1, X12, ldx12,
        tauq2[ iofftauq2 + i - 1 ], X22, ldx22, work,
        ioffx12 + i - 1 + ( i - 1 ) * ldx12,
        ioffx22 + i + ( i - 1 ) * ldx22, ioffwork );
    }
    for ( i = q + 1; i <= p; i ++ ) {
      Blas1.dscal( m - q - i + 1, - z1 * z4, X12, ldx12,
        ioffx12 + i - 1 + ( i - 1 ) * ldx12 );
      alpha.setValue( X12[ ioffx12 + i- 1 + ( i - 1 ) * ldx12 ] );
      tau.setValue( tauq2[ iofftauq2 + i - 1 ] );
      LaPack1.dlarfgp( m - q - i + 1, alpha, X12, ldx12, tau,
        ioffx12 + i - 1 + i * ldx12 );
      X12[ ioffx12 + i- 1 + ( i - 1 ) * ldx12 ] = alpha.getValue();
      tauq2[ iofftauq2 + i - 1 ] = tau.getValue();
      X12[ ioffx12 + i- 1 + ( i - 1 ) * ldx12 ] = 1.;
      LaPack1.dlarf( 'R', p - i, m - q - i + 1, X12, ldx12,
        tauq2[ iofftauq2 + i - 1 ], X12, ldx12, work,
        ioffx12 + i - 1 + ( i - 1 ) * ldx12,
        ioffx12 + i + ( i - 1 ) * ldx12, ioffwork );
      if ( m - p - q >= 1 ) {
        LaPack1.dlarf( 'R', m - p - q, m - q - i + 1, X12, ldx12,
          tauq2[ iofftauq2 + i - 1 ], X22, ldx22, work,
          ioffx12 + i - 1 + ( i - 1 ) * ldx12,
          ioffx22 + q + ( i - 1 ) * ldx22, ioffwork );
      }
    }
    for ( i = 1; i <= m - p - q; i ++ ) {
      Blas1.dscal( m - p - q - i + 1, z2 * z4, X22, ldx22,
        ioffx22 + q + i - 1 + ( p + i - 1 ) * ldx22 );
      alpha.setValue( X22[ ioffx22 + q + i - 1 + ( p + i - 1 ) * ldx22 ] );
      tau.setValue( tauq2[ iofftauq2 + p + i - 1 ] );
      LaPack1.dlarfgp( m - p - q - i + 1, alpha,X22, ldx22, tau,
        ioffx22 + q + i - 1 + ( p + i ) * ldx22 );
      X22[ ioffx22 + q + i - 1 + ( p + i - 1 ) * ldx22 ] = alpha.getValue();
      tauq2[ iofftauq2 + p + i - 1 ] = tau.getValue();
      X22[ ioffx22 + q + i - 1 + ( p + i - 1 ) * ldx22 ] = 1.;
      LaPack1.dlarf( 'R', m - p - q - i, m - p - q - i + 1, X22, ldx22,
        tauq2[ iofftauq2 + p + i - 1 ], X22, ldx22, work,
        ioffx22 + q + i - 1 + ( p + i - 1 ) * ldx22,
        ioffx22 + q + i + ( p + i - 1 ) * ldx22 );
    }
  } else {
    for ( i = 1; i <= q; i ++ ) {
      if ( i == 1 ) {
        Blas1.dscal( p - i + 1, z1, X11, ldx11,
          ioffx11 + i - 1 + ( i - 1 ) * ioffx11 );
      } else {
        Blas1.dscal( p - i + 1,
          z1 * Math.cos( phi[ ioffphi + i - 2 ] ), X11, ldx11,
          ioffx11 + i - 1 + ( i - 1 ) * ioffx11 );
        Blas1.daxpy( p - i + 1,
          - z1 * z3 * z4 * Math.sin( phi[ ioffphi + i - 2 ] ),
          X12, ldx12, X11, ldx11, ioffx12 + i - 2 + ( i - 1 ) * ldx12,
          ioffx11 + i - 1 + ( i - 1 ) * ldx11 );
      }
      if ( i == 1 ) {
        Blas1.dscal( m - p - i + 1, z2, X21, ldx21,
          ioffx21 + i - 1 + ( i - 1 ) * ioffx21 );
      } else {
        Blas1.dscal( m - p - i + 1,
          z2 * Math.cos( phi[ ioffphi + i - 2 ] ), X21, ldx21,
          ioffx21 + i - 1 + ( i - 1 ) * ioffx21 );
        Blas1.daxpy( m - p - i + 1,
          - z2 * z3 * z4 * Math.sin( phi[ ioffphi + i - 2 ] ),
          X22, ldx22, X21, ldx21, ioffx22 + i - 2 + ( i - 1 ) * ldx22,
          ioffx21 + i - 1 + ( i - 1 ) * ldx21 );
      }
      theta[ iofftheta + i - 1 ] = Math.atan2(
        Blas1.dnrm2( m - p - i + 1, X21, ldx21,
          ioffx21 + i - 1 + ( i - 1 ) * ldx21 ),
        Blas1.dnrm2( p - i + 1, X11, ldx11,
          ioffx11 + i - 1 + ( i - 1 ) * ldx11 ) );
      var alpha = new NumberReference(
        X11[ ioffx11 + i - 1 + ( i - 1 ) * ldx11 ] );
      var tau =
        new NumberReference( taup1[ iofftaup1 + i - 1 ] );
      LaPack1.dlarfgp( p - i + 1, alpha, X11, ldx11, tau,
        ioffx11 + i - 1 + i * ldx11 );
      X11[ ioffx11 + i - 1 + ( i - 1 ) * ldx11 ] = alpha.getValue();
      taup1[ iofftaup1 + i - 1 ] = tau.getValue();
      X11[ ioffx11 + i - 1 + ( i - 1 ) * ldx11 ] = 1.;
      alpha.setValue( X21[ ioffx21 + i - 1 + ( i - 1 ) * ldx21 ] );
      tau.setValue( taup2[ iofftaup2 + i - 1 ] );
      LaPack1.dlarfgp( m - p - i + 1, alpha, X21, ldx21, tau,
        ioffx21 + i - 1 + i * ldx21 );
      X21[ ioffx21 + i - 1 + ( i - 1 ) * ldx21 ] = alpha.getValue();
      taup2[ iofftaup2 + i - 1 ] = tau.getValue();
      X21[ ioffx21 + i - 1 + ( i - 1 ) * ldx21 ] = 1.;
      LaPack1.dlarf( 'R', q - i, p - i + 1, X11, ldx11,
        taup1[ iofftaup1 + i - 1 ], X11, ldx11, work,
        ioffx11 + i - 1 + ( i - 1 ) * ldx11,
        ioffx11 + i + ( i - 1 ) * ldx11, ioffwork );
      LaPack1.dlarf( 'R', m - q - i + 1, p - i + 1, X11, ldx11,
        taup1[ iofftaup1 + i - 1 ], X12, ldx12, work,
        ioffx11 + i - 1 + ( i - 1 ) * ldx11,
        ioffx12 + i - 1 + ( i - 1 ) * ldx12, ioffwork );
      LaPack1.dlarf( 'R', q - i, m - p - i + 1, X21, ldx21,
        taup2[ iofftaup2 + i - 1 ], X21, ldx21, work,
        ioffx21 + i - 1 + ( i - 1 ) * ldx21,
        ioffx21 + i + ( i - 1 ) * ldx21, ioffwork );
      LaPack1.dlarf( 'R', m - q - i + 1, m - p - i + 1, X21, ldx21,
        taup2[ iofftaup2 + i - 1 ], X22, ldx22, work,
        ioffx21 + i - 1 + ( i - 1 ) * ldx21,
        ioffx22 + i - 1 + ( i - 1 ) * ldx22, ioffwork );
      if ( i < q ) {
        Blas1.dscal( q - i,
          - z1 * z3 * Math.sin( theta[ iofftheta + i - 1 ] ),
          X11, 1, ioffx11 + i + ( i - 1 ) * ldx11 );
        Blas1.daxpy( q - i,
          z2 * z3 * Math.cos( theta[ iofftheta + i - 1 ] ),
          X21, 1, X11, 1, ioffx21 + i + ( i - 1 ) * ldx21,
          ioffx11 + i + ( i - 1 ) * ldx11 );
      }
      Blas1.dscal( m - q - i + 1,
        - z1 * z4 * Math.sin( theta[ iofftheta + i - 1 ] ),
        X12, 1, ioffx12 + i - 1 + ( i - 1 ) * ldx12 );
      Blas1.daxpy( m - q - i + 1,
        z2 * z4 * Math.cos( theta[ iofftheta + i - 1 ] ),
        X22, 1, X12, 1, ioffx22 + i - 1 + ( i - 1 ) * ldx22,
        ioffx12 + i - 1 + ( i - 1 ) * ldx12 );
      if ( i < q ) {
        phi[ ioffphi + i - 1 ] = Math.atan2(
          Blas1.dnrm2( q - i, X11, 1,
            ioffx11 + i + ( i - 1 ) * ldx11 ),
          Blas1.dnrm2( m - q - i + 1, X12, 1,
            ioffx12 + i - 1 + ( i - 1 ) * ldx12 ) );
      }
      if ( i < q ) {
        alpha.setValue( X11[ ioffx11 + i + ( i - 1 ) * ldx11 ] );
        tau.setValue( tauq1[ iofftauq1 + i - 1 ] );
        LaPack1.dlarfgp( q - i, alpha, X11, 1, tau,
          ioffx11 + i + 1 + ( i - 1 ) * ldx11 );
        X11[ ioffx11 + i + ( i - 1 ) * ldx11 ] = alpha.getValue();
        tauq1[ iofftauq1 + i - 1 ] = tau.getValue();
        X11[ ioffx11 + i + ( i - 1 ) * ldx11 ] = 1.;
      }
      alpha.setValue( X12[ ioffx12 + i - 1 + ( i - 1 ) * ldx12 ] );
      tau.setValue( tauq2[ iofftauq2 + i - 1 ] );
      LaPack1.dlarfgp( m - q - i + 1, alpha, X12, 1, tau,
        ioffx12 + i + ( i - 1 ) * ldx12 );
      X12[ ioffx12 + i - 1 + ( i - 1 ) * ldx12 ] = alpha.getValue();
      tauq2[ iofftauq2 + i - 1 ] = tau.getValue();
      X12[ ioffx12 + i - 1 + ( i - 1 ) * ldx12 ] = 1.;
      if ( i < q ) {
        LaPack1.dlarf( 'L', q - i, p - i, X11, 1,
          tauq1[ iofftauq1 + i - 1 ], X11, ldx11, work,
          ioffx11 + i + ( i - 1 ) * ldx11, ioffx11 + i + i * ldx11,
          ioffwork );
        LaPack1.dlarf( 'L', q - i, m - p - i, X11, 1,
          tauq1[ iofftauq1 + i - 1 ], X21, ldx21, work,
          ioffx11 + i + ( i - 1 ) * ldx11, ioffx21 + i + i * ldx21,
          ioffwork );
      }
      LaPack1.dlarf( 'L', m - q - i + 1, p - i, X12, 1,
        tauq2[ iofftauq2 + i - 1 ], X12, ldx12, work,
        ioffx12 + i - 1 + ( i - 1 ) * ldx12,
        ioffx12 + i - 1 + i * ldx12, ioffwork );
      LaPack1.dlarf( 'L', m - q - i + 1, m - p - i, X12, 1,
        tauq2[ iofftauq2 + i - 1 ], X22, ldx22, work,
        ioffx12 + i - 1 + ( i - 1 ) * ldx12,
        ioffx22 + i - 1 + i * ldx22, ioffwork );
    }
    for ( i = q + 1; i <= p; i ++ ) {
      Blas1.dscal( m - q - i + 1, - z1 * z4, X12, 1,
        ioffx12 + i - 1 + ( i - 1 ) * ldx12 );
      alpha.setValue( X12[ ioffx12 + i- 1 + ( i - 1 ) * ldx12 ] );
      tau.setValue( tauq2[ iofftauq2 + i - 1 ] );
      LaPack1.dlarfgp( m - q - i + 1, alpha, X12, 1, tau,
        ioffx12 + i + ( i - 1 ) * ldx12 );
      X12[ ioffx12 + i - 1 + ( i - 1 ) * ldx12 ] = alpha.getValue();
      tauq2[ iofftauq2 + i - 1 ] = tau.getValue();
      X12[ ioffx12 + i - 1 + ( i - 1 ) * ldx12 ] = 1.;
      LaPack1.dlarf( 'L', m - q - i + 1, p - i, X12, 1,
        tauq2[ iofftauq2 + i - 1 ], X12, ldx12, work,
        ioffx12 + i - 1 + ( i - 1 ) * ldx12,
        ioffx12 + i - 1 + i * ldx12, ioffwork );
      if ( m - p - q >= 1 ) {
        LaPack1.dlarf( 'L', m - q - i + 1, m - p - q, X12, 1,
          tauq2[ iofftauq2 + i - 1 ], X22, ldx22, work,
          ioffx12 + i - 1 + ( i - 1 ) * ldx12,
          ioffx22 + i - 1 + q * ldx22, ioffwork );
      }
    }
    for ( i = 1; i <= m - p - q; i ++ ) {
      Blas1.dscal( m - p - q - i + 1, z2 * z4, X22, 1,
        ioffx22 + p + i - 1 + ( q + i - 1 ) * ldx22 );
      alpha.setValue( X22[ ioffx22 + p + i - 1 + ( q + i - 1 ) * ldx22 ] );
      tau.setValue( tauq2[ iofftauq2 + p + i - 1 ] );
      LaPack1.dlarfgp( m - p - q - i + 1, alpha, X22, 1, tau,
        ioffx22 + p + i + ( q + i - 1 ) * ldx22 );
      X22[ ioffx22 + p + i - 1 + ( q + i - 1 ) * ldx22 ] = alpha.getValue();
      tauq2[ iofftauq2 + p + i - 1 ] = tau.getValue();
      X22[ ioffx22 + p + i - 1 + ( q + i - 1 ) * ldx22 ] = 1.;
      LaPack1.dlarf( 'L', m - p - q - i, m - p - q - i + 1, X22, 1,
        tauq2[ iofftauq2 + p + i - 1 ], X22, ldx22, work,
        ioffx22 + p + i - 1 + ( q + i - 1 ) * ldx22,
        ioffx22 + p + i - 1 + ( q + i ) * ldx22 );
    }
  }
}
//************************************************************************
LaPack2.dorg2l = function( m, n, k, A, lda, tau, work, info,
ioffa, iofftau, ioffwork ) {
  info.setValue( 0 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 || n > m ) info.setValue( -2 );
  else if ( k < 0 || k > n ) info.setValue( -3 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -5 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dorg2l', - info.getValue() );
    return;
  }
  if ( n <= 0 ) return;
  for ( var j = 1; j <= n - k; j ++ ) {
    for ( var l = 1; l <= m; l ++ ) {
      A[ ioffa + l - 1 + ( j - 1 ) * lda ] = 0.;
    }
    A[ ioffa + m - n + j - 1 + ( j - 1 ) * lda ] = 1.;
  }
  for ( var i = 1; i <= k; i ++ ) {
    var ii = n - k + i;
    A[ ioffa + m - n + ii - 1 + ( ii - 1 ) * lda ] = 1.;
    LaPack1.dlarf( 'Left', m - n + ii, ii - 1, A, 1,
      tau[ iofftau + i - 1 ], A, lda, work,
      ioffa + ( ii - 1 ) * lda, ioffa, ioffwork );
    Blas1.dscal( m - n + ii - 1, -tau[ iofftau + i - 1 ], A, 1,
      ioffa + ( ii - 1 ) * lda );
    A[ ioffa + m - n + ii - 1 + ( ii - 1 ) * lda ] =
      1. - tau[ iofftau + i - 1 ];
    for ( l = m - n + ii + 1; l <= m; l ++ ) {
      A[ ioffa + l - 1 + ( ii - 1 ) * lda ] = 0.;
    }
  }
}
//************************************************************************
LaPack2.dorg2r = function( m, n, k, A, lda, tau, work, info,
ioffa, iofftau, ioffwork ) {
  info.setValue( 0 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 || n > m ) info.setValue( -2 );
  else if ( k < 0 || k > n ) info.setValue( -3 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -5 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dorg2r', - info.getValue() );
    return;
  }
  if ( n <= 0 ) return;
  for ( var j = k + 1; j <= n; j ++ ) {
    for ( var l = 1; l <= m; l ++ ) {
      A[ ioffa + l - 1 + ( j - 1 ) * lda ] = 0.;
    }
    A[ ioffa + j - 1 + ( j - 1 ) * lda ] = 1.;
  }
  for ( var i = k; i >= 1; i -- ) {
    if ( i < n ) {
      A[ ioffa + i - 1 + ( i - 1 ) * lda ] = 1.;
      LaPack1.dlarf( 'Left', m - i + 1, n - i, A, 1,
        tau[ iofftau + i - 1 ], A, lda, work,
        ioffa + i - 1 + ( i - 1 ) * lda,
        ioffa + i - 1 + i * lda, ioffwork );
    }
    if ( i < m ) {
      Blas1.dscal( m - i, -tau[ iofftau + i - 1 ], A, 1,
        ioffa + i + ( i - 1 ) * lda );
    }
    A[ ioffa + i - 1 + ( i - 1 ) * lda ] = 1. - tau[ iofftau + i - 1 ];
    for ( l = 1; l <= i - 1; l ++ ) {
      A[ ioffa + l - 1 + ( i - 1 ) * lda ] = 0.;
    }
  }
}
//************************************************************************
LaPack2.dorgl2 = function( m, n, k, A, lda, tau, work, info,
ioffa, iofftau, ioffwork ) {
  info.setValue( 0 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < m ) info.setValue( -2 );
  else if ( k < 0 || k > m ) info.setValue( -3 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -5 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dorgl2', - info.getValue() );
    return;
  }
  if ( m <= 0 ) return;
  if ( k < m ) {
    for ( var j = 1; j <= n; j ++ ) {
      for ( var l = k + 1; l <= m; l ++ ) {
        A[ ioffa + l - 1 + ( j - 1 ) * lda ] = 0.;
      }
      if ( j > k && j <= m ) {
        A[ ioffa + j - 1 + ( j - 1 ) * lda ] = 1.;
      }
    }
  }
  for ( var i = k; i >= 1; i -- ) {
    if ( i < n ) {
      if ( i < m ) {
        A[ ioffa + i - 1 + ( i - 1 ) * lda ] = 1.;
        LaPack1.dlarf( 'Right', m - i, n - i + 1, A, lda,
          tau[ iofftau + i - 1 ], A, lda,
          work, ioffa + i - 1 + ( i - 1 ) * lda,
          ioffa + i + ( i  - 1 )* lda, ioffwork );
      }
      Blas1.dscal( n - i, -tau[ iofftau + i - 1 ], A, lda,
        ioffa + i - 1  + i * lda );
    }
    A[ ioffa + i - 1 + ( i - 1 ) * lda ] = 1. - tau[ iofftau + i - 1 ];
    for ( l = 1; l <= i - 1; l ++ ) {
      A[ ioffa + i - 1 + ( l - 1 ) * lda ] = 0.;
    }
  }
}
//************************************************************************
LaPack2.dorgr2 = function( m, n, k, A, lda, tau, work, info,
ioffa, iofftau, ioffwork ) {
  info.setValue( 0 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < m ) info.setValue( -2 );
  else if ( k < 0 || k > m ) info.setValue( -3 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -5 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dorgr2', - info.getValue() );
    return;
  }
  if ( m <= 0 ) return;
  if ( k < m ) {
    for ( var j = 1; j <= n; j ++ ) {
      for ( var l = 1; l <= m - k; l ++ ) {
        A[ ioffa + l - 1 + ( j - 1 ) * lda ] = 0.;
      }
      if ( j > n - m && j <= n - k ) {
        A[ ioffa + m - n + j - 1 + ( j - 1 ) * lda ] = 1.;
      }
    }
  }
  for ( var i = 1; i <= k; i ++ ) {
    var ii = m - k + i;
    A[ ioffa + ii - 1 + ( n - m + ii - 1 ) * lda ] = 1.;
    LaPack1.dlarf( 'Right', ii - 1, n - m + ii, A, lda,
      tau[ iofftau + i - 1 ], A, lda, work,
      ioffa + ii - 1, ioffa, ioffwork );
    Blas1.dscal( n - m + ii - 1, -tau[ iofftau + i - 1 ], A, lda,
      ioffa + ii - 1 );
    A[ ioffa + ii - 1 + ( n - m + ii - 1 ) * lda ] =
      1. - tau[ iofftau + i - 1 ];
    for ( l = n - m + ii + 1; l <= n; l ++ ) {
      A[ ioffa + ii - 1 + ( l - 1 ) * lda ] = 0.;
    }
  }
}
//************************************************************************
LaPack2.dorm2l = function( side, trans, m, n, k, A, lda, tau,
C, ldc, work, info, ioffa, iofftau, ioffc, ioffwork) {
  info.setValue( 0 );
  var left = ( side.charAt(0).toUpperCase() == 'L' );
  var notran = ( trans.charAt(0).toUpperCase() == 'N' );
  var nq = ( left ? m : n );
  if ( ! left && side.charAt(0).toUpperCase() != 'R' ) {
    info.setValue( -1 );
  } else if ( ! notran && trans.charAt(0).toUpperCase() != 'T' ) {
    info.setValue( -2 );
  } else if ( m < 0 ) info.setValue( -3 );
  else if ( n < 0 ) info.setValue( -4 );
  else if ( k < 0 || k > nq ) info.setValue( -5 );
  else if ( lda < Math.max( 1, nq ) ) info.setValue( -7 );
  else if ( ldc < Math.max( 1, m ) ) info.setValue( -10 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dorm2l', - info.getValue() );
    return;
  }
  if ( m == 0 || n == 0 || k == 0 ) return;
  if ( ( left && notran ) || ( ! left && ! notran) ) {
    var i1 = 1;
    var i2 = k;
    var i3 = 1;
  } else {
    i1 = k;
    i2 = 1;
    i3 = -1;
  }
  if ( left ) var ni = n;
  else var mi = m;
  for ( var i = i1; ( i3 == 1 && i <= i2 ) ||
  ( i3 == -1 && i >= i2 ); i += i3 ) {
    if ( left) mi = m - k + i;
    else ni = n - k + i;
    var aii = A[ ioffa + nq - k + i - 1 + ( i - 1 ) * lda ];
    A[ ioffa + nq - k + i - 1 + ( i - 1 ) * lda ] = 1.;
    LaPack1.dlarf( side, mi, ni, A, 1, tau[ iofftau + i - 1 ], C, ldc,
      work, ioffa + ( i - 1 ) * lda, ioffc, ioffwork );
    A[ ioffa + nq - k + i - 1 + ( i - 1 ) * lda ] = aii;
  }
}
//************************************************************************
LaPack2.dorm2r = function( side, trans, m, n, k, A, lda, tau,
C, ldc, work, info, ioffa, iofftau, ioffc, ioffwork) {
  info.setValue( 0 );
  var left = ( side.charAt(0).toUpperCase() == 'L' );
  var notran = ( trans.charAt(0).toUpperCase() == 'N' );
  var nq = ( left ? m : n );
  if ( ! left && side.charAt(0).toUpperCase() != 'R' ) {
    info.setValue( -1 );
  } else if ( ! notran && trans.charAt(0).toUpperCase() != 'T' ) {
    info.setValue( -2 );
  } else if ( m < 0 ) info.setValue( -3 );
  else if ( n < 0 ) info.setValue( -4 );
  else if ( k < 0 || k > nq ) info.setValue( -5 );
  else if ( lda < Math.max( 1, nq ) ) info.setValue( -7 );
  else if ( ldc < Math.max( 1, m ) ) info.setValue( -10 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dorm2r', - info.getValue() );
    return;
  }
  if ( m == 0 || n == 0 || k == 0 ) return;
  if ( ( left && ! notran ) || ( ! left && notran) ) {
    var i1 = 1;
    var i2 = k;
    var i3 = 1;
  } else {
    i1 = k;
    i2 = 1;
    i3 = -1;
  }
  if ( left ) {
    var ni = n;
    var jc = 1;
  } else {
    var mi = m;
    var ic = 1;
  }
  for ( var i = i1; ( i3 == 1 && i <= i2 ) ||
  ( i3 == -1 && i >= i2 ); i += i3 ) {
    if ( left) {
      mi = m - i + 1;
      ic = i;
    } else {
      ni = n - i + 1;
      jc = i;
    }
    var aii = A[ ioffa + i - 1 + ( i - 1 ) * lda ];
    A[ ioffa + i - 1 + ( i - 1 ) * lda ] = 1.;
    LaPack1.dlarf( side, mi, ni, A, 1, tau[ iofftau + i - 1 ], C, ldc,
      work, ioffa + i - 1 + ( i - 1 ) * lda,
      ioffc + ic - 1 + ( jc - 1 ) * ldc, ioffwork );
    A[ ioffa + i - 1 + ( i - 1 ) * lda ] = aii;
  }
}
//************************************************************************
LaPack2.dorml2 = function( side, trans, m, n, k, A, lda, tau,
C, ldc, work, info, ioffa, iofftau, ioffc, ioffwork) {
  info.setValue( 0 );
  var left = ( side.charAt(0).toUpperCase() == 'L' );
  var notran = ( trans.charAt(0).toUpperCase() == 'N' );
  var nq = ( left ? m : n );
  if ( ! left && side.charAt(0).toUpperCase() != 'R' ) {
    info.setValue( -1 );
  } else if ( ! notran && trans.charAt(0).toUpperCase() != 'T' ) {
    info.setValue( -2 );
  } else if ( m < 0 ) info.setValue( -3 );
  else if ( n < 0 ) info.setValue( -4 );
  else if ( k < 0 || k > nq ) info.setValue( -5 );
  else if ( lda < Math.max( 1, k ) ) info.setValue( -7 );
  else if ( ldc < Math.max( 1, m ) ) info.setValue( -10 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dorml2', - info.getValue() );
    return;
  }
  if ( m == 0 || n == 0 || k == 0 ) return;
  if ( ( left && notran ) || ( ! left && ! notran) ) {
    var i1 = 1;
    var i2 = k;
    var i3 = 1;
  } else {
    i1 = k;
    i2 = 1;
    i3 = -1;
  }
  if ( left ) {
    var ni = n;
    var jc = 1;
  } else {
    var mi = m;
    var ic = 1;
  }
  for ( var i = i1; ( i3 == 1 && i <= i2 ) ||
  ( i3 == -1 && i >= i2 ); i += i3 ) {
    if ( left) {
      mi = m - i + 1;
      ic = i;
    } else {
      ni = n - i + 1;
      jc = i;
    }
    var aii = A[ ioffa + i - 1 + ( i - 1 ) * lda ];
    A[ ioffa + i - 1 + ( i - 1 ) * lda ] = 1.;
    LaPack1.dlarf( side, mi, ni, A, lda, tau[ iofftau + i - 1 ],
      C, ldc, work, ioffa + i - 1 + ( i - 1 ) * lda,
      ioffc + ic - 1 + ( jc - 1 ) * ldc, ioffwork );
    A[ ioffa + i - 1 + ( i - 1 ) * lda ] = aii;
  }
}
//************************************************************************
LaPack2.dormlq = function( side, trans, m, n, k, A, lda, tau,
C, ldc, work, lwork, info, ioffa, iofftau, ioffc, ioffwork ) {
  throw new Error("not tested: complicated input");
  var nbmax = 64;
  var ldt = nbmax + 1;
  var T = new Array( ldt * nbmax );
  info.setValue( 0 );
  var left = ( side.charAt(0).toUpperCase() == 'L' );
  var notran = ( trans.charAt(0).toUpperCase() == 'N' );
  var lquery = ( lwork == -1 );
  if ( left ) {
    var nq = m;
    var nw = n;
  } else {
    nq = n;
    nw = m;
  }
  if ( ! left && side.charAt(0).toUpperCase() != 'R' ) info.setValue( -1 );
  else if ( ! notran && trans.charAt(0).toUpperCase() != 'T' ) {
    info.setValue( -2 );
  } else if ( m < 0 ) info.setValue( -3 );
  else if ( n < 0 ) info.setValue( -4 );
  else if ( k < 0 || k > nq ) info.setValue( -5 );
  else if ( lda < Math.max( 1, k ) ) info.setValue( -7 );
  else if ( ldc < Math.max( 1, m ) ) info.setValue( -10 );
  else if ( lwork < Math.max( 1, nw ) && ! lquery ) info.setValue( -12 );
  if ( info.getValue() == 0 ) {
    var nb = Math.min( nbmax,
      LaPack0.ilaenv( 1, 'dormlq', side + trans, m, n, k, -1 ) );
    var lwkopt = Math.max( 1, nw ) * nb;
    work[ ioffwork ] = lwkopt;
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dormlq', -info.getValue() );
    return;
  } else if ( lquery ) return;
  if ( m == 0 || n == 0 || k == 0 ) {
    work[ ioffwork ] = 1;
    return;
  }
  var nbmin = 2;
  var ldwork = nw;
  if ( nb > 1 && nb < k ) {
    var iws = nw * nb;
    if ( lwork < iws ) {
      nb = lwork / ldwork;
      nbmin = Math.max( 2, 
        LaPack0.ilaenv( 2, 'dormlq', side + trans, m, n, k, -1 ) );
    }
  } else iws = nw;
  if ( nb < nbmin || nb >= k ) {
    var iinfo = new IntReference;
    LaPack1.dorml2( side, trans, m, n, k, A, lda, tau, C, ldc, work,
      iinfo, ioffa, iofftau, ioffc, ioffwork );
  } else {
    if ( ( left && notran ) || ( ! left && ! notran ) ) {
      var i1 = 1;
      var i2 = k;
      var i3 = nb;
    } else {
      i1 = ( ( k - 1 ) / nb ) * nb + 1;
      i2 = 1;
      i3 = -nb;
    }
    if ( left ) {
      var ni = n;
      var jc = 1;
    } else {
      var mi = m;
      var ic = 1;
    }
    var transt = ( notran ? 'T' : 'N' );
    for ( var i = i1; ( i3 == nb && i <= i2 ) ||
    ( i3 == -nb && i >= i2 ); i += i3 ) {
      var ib = Math.min( nb, k - i + 1 );
      LaPack0.dlarft( 'Forward', 'Rowwise', nq - i + 1, ib,
        A, lda, tau, T, ldt, ioffa + i - 1 + ( i - 1 ) * lda,
        iofftau + i - 1, 0 );
      if ( left ) {
        mi = m - i + 1;
        ic = i;
      } else {
        ni = n - i + 1;
        jc = i;
      }
      LaPack1.dlarfb( side, transt, 'Forward', 'Rowwise', mi, ni, ib,
        A, lda, T, ldt, C, ldc, work, ldwork,
        ioffa + i - 1 + ( i - 1 ) * lda, 0,
        ioffc + ic - 1 + ( jc - 1 ) * ldc, ioffwork );
    }
  }
  work[ ioffwork ] = lwkopt;
}
//************************************************************************
LaPack2.dormql = function( side, trans, m, n, k, A, lda, tau,
C, ldc, work, lwork, info, ioffa, iofftau, ioffc, ioffwork ) {
  throw new Error("not tested: complicated input");
  var nbmax = 64;
  var ldt = nbmax + 1;
  var T = new Array( ldt * nbmax );
  info.setValue( 0 );
  var left = ( side.charAt(0).toUpperCase() == 'L' );
  var notran = ( trans.charAt(0).toUpperCase() == 'N' );
  var lquery = ( lwork == -1 );
  if ( left ) {
    var nq = m;
    var nw = Math.max( 1, n );
  } else {
    nq = n;
    nw = Math.max( 1, m );
  }
  if ( ! left && side.charAt(0).toUpperCase() != 'R' ) info.setValue( -1 );
  else if ( ! notran && trans.charAt(0).toUpperCase() != 'T' ) {
    info.setValue( -2 );
  } else if ( m < 0 ) info.setValue( -3 );
  else if ( n < 0 ) info.setValue( -4 );
  else if ( k < 0 || k > nq ) info.setValue( -5 );
  else if ( lda < Math.max( 1, nq ) ) info.setValue( -7 );
  else if ( ldc < Math.max( 1, m ) ) info.setValue( -10 );
  if ( info.getValue() == 0 ) {
    if ( m == 0 || n == 0 ) var lwkopt = 1;
    else {
      var nb = Math.min( nbmax,
        LaPack0.ilaenv( 1, 'dormql', side + trans, m, n, k, -1 ) );
      lwkopt = nw * nb;
    }
    work[ ioffwork ] = lwkopt;
    if ( lwork < nw && ! lquery) info.setValue( -12 );
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dormql', -info.getValue() );
    return;
  } else if ( lquery ) return;
  if ( m == 0 || n == 0 ) return;
  var nbmin = 2;
  var ldwork = nw;
  if ( nb > 1 && nb < k ) {
    var iws = nw * nb;
    if ( lwork < iws ) {
      nb = lwork / ldwork;
      nbmin = Math.max( 2,
        LaPack0.ilaenv( 2, 'dormql', side + trans, m, n, k, -1 ) );
    }
  } else iws = nw;
  if ( nb < nbmin || nb >= k ) {
    var iinfo = new IntReference;
    LaPack1.dorm2l( side, trans, m, n, k, A, lda, tau, C, ldc, work,
      iinfo, ioffa, iofftau, ioffc, ioffwork );
  } else {
    if ( ( left && notran ) || ( ! left && ! notran ) ) {
      var i1 = 1;
      var i2 = k;
      var i3 = nb;
    } else {
      i1 = ( ( k - 1 ) / nb ) * nb + 1;
      i2 = 1;
      i3 = -nb;
    }
    if ( left ) var ni = n;
    else var mi = m;
    for ( var i = i1; ( i3 == nb && i <= i2 ) ||
    ( i3 == -nb && i >= i2 ); i += i3 ) {
      var ib = Math.min( nb, k - i + 1 );
      LaPack0.dlarft( 'Backward', 'Columnwise', nq - k + i + ib - 1,
        ib, A, lda, tau, T, ldt, ioffa + ( i - 1 ) * lda,
        iofftau + i - 1, 0 );
      if ( left ) mi = m - k + i + ib - 1;
      else ni = n - k + i + ib - 1;
      LaPack1.dlarfb( side, trans, 'Backward', 'Columnwise', mi, ni,
        ib, A, lda, T, ldt, C, ldc, work, ldwork,
        ioffa + ( i - 1 ) * lda, 0, ioffc, ioffwork );
    }
  }
  work[ ioffwork ] = lwkopt;
}
//************************************************************************
LaPack2.dormr2 = function( side, trans, m, n, k, A, lda, tau,
C, ldc, work, info, ioffa, iofftau, ioffc, ioffwork) {
  info.setValue( 0 );
  var left = ( side.charAt(0).toUpperCase() == 'L' );
  var notran = ( trans.charAt(0).toUpperCase() == 'N' );
  var nq = ( left ? m : n );
  if ( ! left && side.charAt(0).toUpperCase() != 'R' ) info.setValue( -1 );
  else if ( ! notran && trans.charAt(0).toUpperCase() != 'T' ) {
    info.setValue( -2 );
  } else if ( m < 0 ) info.setValue( -3 );
  else if ( n < 0 ) info.setValue( -4 );
  else if ( k < 0 || k > nq ) info.setValue( -5 );
  else if ( lda < Math.max( 1, k ) ) info.setValue( -7 );
  else if ( ldc < Math.max( 1, m ) ) info.setValue( -10 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dormr2', -info.getValue() );
    return;
  }
  if ( m == 0 || n == 0 || k == 0 ) return;
  if ( ( left && ! notran ) || ( ! left && notran ) ) {
    var i1 = 1;
    var i2 = k;
    var i3 = 1;
  } else {
    i1 = k;
    i2 = 1;
    i3 = -1;
  }
  if ( left ) var ni = n;
  else var mi = m;
  for ( var i = i1; ( i3 == 1 && i <= i2 ) ||
  ( i3 == -1 && i >= i2 ); i += i3 ) {
    if ( left )  mi = m - k + i;
    else ni = n - k + i;
    var aii = A[ ioffa + i- 1 + ( nq - k + i - 1 ) * lda ];
    A[ ioffa + i- 1 + ( nq - k + i - 1 ) * lda ] = 1.;
    LaPack1.dlarf( side, mi, ni, A, lda, tau[ iofftau + i - 1 ],
      C, ldc, work, ioffa + i - 1, ioffc, ioffwork );
    A[ ioffa + i - 1 + ( nq - k + i - 1 ) * lda ] = aii;
  }
}
//************************************************************************
LaPack2.dormrq = function( side, trans, m, n, k, A, lda, tau,
C, ldc, work, lwork, info, ioffa, iofftau, ioffc, ioffwork ) {
  throw new Error("not tested: complicated input");
  var nbmax = 64;
  var ldt = nbmax + 1;
  var T = new Array( ldt * nbmax );
  info.setValue( 0 );
  var left = ( side.charAt(0).toUpperCase() == 'L' );
  var notran = ( trans.charAt(0).toUpperCase() == 'N' );
  var lquery = ( lwork == -1 );
  if ( left ) {
    var nq = m;
    var nw = Math.max( 1, n );
  } else {
    nq = n;
    nw = Math.max( 1, m );
  }
  if ( ! left && side.charAt(0).toUpperCase() != 'R' ) info.setValue( -1 );
  else if ( ! notran && trans.charAt(0).toUpperCase() != 'T' ) {
    info.setValue( -2 );
  } else if ( m < 0 ) info.setValue( -3 );
  else if ( n < 0 ) info.setValue( -4 );
  else if ( k < 0 || k > nq ) info.setValue( -5 );
  else if ( lda < Math.max( 1, k ) ) info.setValue( -7 );
  else if ( ldc < Math.max( 1, m ) ) info.setValue( -10 );
  if ( info.getValue() == 0 ) {
    if ( m == 0 || n == 0 ) var lwkopt = 1;
    else {
      var nb = Math.min( nbmax,
        LaPack0.ilaenv( 1, 'dormrq', side + trans, m, n, k, -1 ) );
      lwkopt = nw * nb;
    }
    work[ ioffwork ] = lwkopt;
    if ( lwork < nw  && ! lquery ) info.setValue( -12 );
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dormrq', -info.getValue() );
    return;
  } else if ( lquery ) return;
  if ( m == 0 || n == 0 ) return;
  var nbmin = 2;
  var ldwork = nw;
  if ( nb > 1 && nb < k ) {
    var iws = nw * nb;
    if ( lwork < iws ) {
      nb = lwork / ldwork;
      nbmin = Math.max( 2,
        LaPack0.ilaenv( 2, 'dormrq', side + trans, m, n, k, -1 ) );
    }
  } else iws = nw;
  if ( nb < nbmin || nb >= k ) {
    var iinfo = new IntReference;
    LaPack1.dormr2( side, trans, m, n, k, A, lda, tau, C, ldc, work,
      iinfo, ioffa, iofftau, ioffc, ioffwork );
  } else {
    if ( ( left && ! notran ) || ( ! left && notran ) ) {
      var i1 = 1;
      var i2 = k;
      var i3 = nb;
    } else {
      i1 = ( ( k - 1 ) / nb ) * nb + 1;
      i2 = 1;
      i3 = -nb;
    }
    if ( left ) var ni = n;
    else var mi = m;
    var trans = ( notran ? 'T' : 'N' );
    for ( var i = i1; ( i3 == nb && i <= i2 ) ||
    ( i3 == -nb && i >= i2 ); i += i3 ) {
      var ib = Math.min( nb, k - i + 1 );
      LaPack0.dlarft( 'Backward', 'Rowwise', nq - k + i + ib - 1, ib,
        A, lda, tau, T, ldt, ioffa + i - 1, iofftau + i - 1, 0 );
      if ( left ) mi = m - k + i + ib - 1;
      else ni = n - k + i + ib - 1;
      LaPack1.dlarfb( side, trans, 'Backward', 'Rowwise', mi, ni, ib,
        A, lda, T, ldt, C, ldc, work, ldwork, ioffa + i - 1, 0, ioffc,
        ioffwork );
    }
  }
  work[ 0 ] = lwkopt;
}
//************************************************************************
LaPack2.dormrz = function( side, trans, m, n, k, l, A, lda,
tau, C, ldc, work, lwork, info, ioffa, iofftau, ioffc, ioffwork ) {
  throw new Error("not tested: complicated input");
  var nbmax = 64;
  var ldt = nbmax + 1;
  var iofft = 0;
  var T = new Array( iofft + ldt * nbmax );
  info.setValue( 0 );
  var left = ( side.charAt(0).toUpperCase() == 'L' );
  var notran = ( trans.charAt(0).toUpperCase() == 'N' );
  var lquery = ( lwork == -1 );
  if ( left ) {
    var nq = m;
    var nw = Math.max( 1, n );
  } else {
    nq = n;
    nw = Math.max( 1, m );
  }
  if ( ! left && side.charAt(0).toUpperCase() != 'R' ) info.setValue( -1 );
  else if ( ! notran && trans.charAt(0).toUpperCase() != 'T' ) {
    info.getValue() == -2;
  } else if ( m < 0 ) info.setValue( -3 );
  else if ( n < 0 ) info.setValue( -4 );
  else if ( k < 0 || k > nq ) info.setValue( -5 );
  else if ( l < 0 || ( left && l > m ) || ( ! left && l > n ) ) {
    info.setValue( -6 );
  } else if ( lda < Math.max( 1, k ) ) info.setValue( -8 );
  else if ( ldc < Math.max( 1, m ) ) info.setValue( -11 );
  if ( info.getValue() == 0 ) {
    if ( m == 0 || n == 0 ) var lwkopt = 1;
    else {
      var nb = Math.min( nbmax,
        LaPack0.ilaenv( 1, 'dormrq', side + trans, m, n, k, -1 ) );
      lwkopt = nw * nb;
    }
    work[ ioffwork ] = lwkopt;
    if ( lwork < Math.max( 1, nw ) && ! lquery ) info.setValue( -13 );
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dormrz', - info.getValue() );
    return;
  } else if ( lquery ) return;
  if ( m == 0 || n == 0 ) {
    work[ ioffwork ] = 1;
    return;
  }
  var nbmin = 2;
  var ldwork = nw;
  if ( nb > 1 && nb < k ) {
    var iws = nw * nb;
    if ( lwork < iws ) {
      nb = lwork / ldwork;
      nbmin = Math.max( 2,
        LaPack0.ilaenv( 2, 'dormrq', side + trans, m, n, k, -1 ) );
    }
  } else iws = nw;
  if ( nb < nbmin || nb >= k ) {
    var iinfo = new IntReference();
    LaPack1.dormr3( side, trans, m, n, k, l, A, lda, tau, C, ldc,
      work, iinfo, ioffa, iofftau, ioffc, ioffwork );
  } else {
    if ( ( left && ! notran ) || ( ! left && notran ) ) {
      var i1 = 1;
      var i2 = k
      var i3 = nb;
    } else {
      i1 = ( ( k - 1 ) / nb ) * nb + 1;
      i2 = 1;
      i3 = - nb;
    }
    if ( left ) {
      var ni = n;
      var jc = 1;
      var ja = m - l + 1;
    } else {
      var mi = m;
      var ic = 1;
      ja = n - l + 1;
    }
    var transt = ( notran ? 'T' : 'N' );
    for ( var i = i1;
    ( i3 == nb && i <= i2 ) || ( i3 == - nb && i >= i2 ); i += i3 ) {
      var ib = Math.min( nb, k - i + 1 );
      LaPack0.dlarzt( 'Backward', 'Rowwise', l, ib, A, lda, tau,
        T, ldt, ioffa + i - 1 + ( ja - 1 ) * lda, iofftau + i - 1,
        iofft );
      if ( left ) {
        mi = m - i + 1;
        ic = i;
      } else {
        ni = n - i + 1;
        jc = i;
      }
      LaPack0.dlarzb( side, transt, 'Backward', 'Rowwise', mi, ni, ib,
        l, A, lda, T, ldt, C, ldc, work, ldwork,
        ioffa + i - 1 + ( ja - 1 ) * lda, iofft,
        ioffc + ic - 1 + ( jc - 1 ) * ldc, ioffwork );
    }
  }
  work[ ioffwork ] = lwkopt;
}
//************************************************************************
LaPack2.dpbcon = function( uplo, n, kd, AB, ldab, anorm, rcond,
work, iwork, info ) {
  throw new Error("not programmed: band matrix");
}
//************************************************************************
LaPack2.dpocon = function( uplo, n, A, lda, anorm, rcond, work,
iwork, info, ioffa, ioffwork, ioffiwork ) {
  var isave = new Array( 3 );
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -4 );
  else if ( anorm < 0. ) info.setValue( -5 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dpocon', -info.getValue() );
    return;
  }
  rcond.setValue( 0. );
  if ( n == 0 ) {
    rcond.setValue( 1. );
    return;
  } else if ( anorm == 0. ) return;
  var smlnum = LaPack0.dlamch( 'Safe minimum' );
  var kase = new IntReference(0);
  var normin = 'N';
  var ainvnm = new NumberReference();
  while ( true ) { // 10
    LaPack0.dlacn2( n, work, work, iwork, ainvnm, kase, isave,
      ioffwork + n, ioffwork, ioffiwork, 0 );
    if ( kase.getValue() != 0 ) {
      var scalel = new NumberReference();
      var scaleu = new NumberReference();
      if ( upper ) {
        LaPack1.dlatrs( 'Upper', 'Transpose', 'Non-unit', normin, n, A,
          lda, work, scalel, work, info, ioffa, ioffwork,
          ioffwork + 2 * n );
        normin = 'Y';
        LaPack1.dlatrs( 'Upper', 'No transpose', 'Non-unit', normin, n,
          A, lda, work, scaleu, work, info, ioffa, ioffwork,
          ioffwork + 2 * n );
      } else {
        LaPack1.dlatrs( 'Lower', 'No transpose', 'Non-unit', normin, n,
          A, lda, work, scalel, work, info, ioffa, ioffwork,
          ioffwork + 2 * n );
        normin = 'Y';
        LaPack1.dlatrs( 'Lower', 'Transpose', 'Non-unit', normin, n, A,
          lda, work, scaleu, work, info, ioffa, ioffwork,
          ioffwork + 2 * n );
      }
      var scale = scalel.getValue() * scaleu.getValue();
      if ( scale != 1. ) {
        var ix = Blas1.idamax( n, work, 1, ioffwork );
        if ( scale < Math.abs( work[ ioffwork + ix - 1 ] ) * smlnum ||
        scale == 0. ) {
          break;
        }
        LaPack1.drscl( n, scale, work, 1, ioffwork );
      }
      continue;
    } else break;
  }
  if ( ainvnm.getValue() != 0. ) {
    rcond.setValue( ( 1. / ainvnm.getValue() ) / anorm );
  }
}
//************************************************************************
LaPack2.dppcon = function( uplo, n, AP, anorm, rcond, work,
iwork, info ) {
  throw new Error("not programmed: packed matrix");
}
LaPack2.zppcon = function( uplo, n, AP, anorm, rcond, work,
iwork, info ) {
  throw new Error("not programmed: complex packed matrix");
}
//************************************************************************
LaPack2.dptrfs = function( n, nrhs, d, e, df, ef, B, ldb, X,
ldx, ferr, berr, work, info, ioffd, ioffe, ioffdf, ioffef, ioffb, ioffx,
ioffferr, ioffberr, ioffwork ) {
  var itmax = 5;
  info.setValue( 0 );
  if ( n < 0 ) info.setValue( -1 );
  else if ( nrhs < 0 ) info.setValue( -2 );
  else if ( ldb < Math.max( 1, n ) ) info.setValue( -8 );
  else if ( ldx < Math.max( 1, n ) ) info.setValue( -10 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dptrfs', -info.getValue() );
    return;
  }
  if ( n == 0 || nrhs == 0 ) {
    for ( var j = 1; j <= nrhs; j ++ ) {
      ferr[ ioffferr + j - 1 ] = 0.;
      berr[ ioffberr + j - 1 ] = 0.;
    }
    return;
  }
  var nz = 4;
  var eps = LaPack0.dlamch( 'Epsilon' );
  var safmin = LaPack0.dlamch( 'Safe minimum' );
  var safe1 = nz * safmin;
  var safe2 = safe1 / eps;
  for ( j = 1; j <= nrhs; j ++ ) {
    var count = 1;
    var lstres = 3.;
    while ( true ) {
      if ( n == 1 ) {
        var bi = B[ ioffb + ( j - 1 ) * ldb ];
        var dx = d[ ioffd ] * X[ ioffx + ( j - 1 ) * ldx ];
        work[ ioffwork + n ] = bi - dx;
        work[ ioffwork ] = Math.abs( bi ) + Math.abs( dx );
      } else {
        bi = B[ ioffb + ( j - 1 ) * ldb ];
        dx = d[ ioffd ] * X[ ioffx + ( j - 1 ) * ldx ];
        var ex = e[ ioffe ] * X[ ioffx + 1 + ( j - 1 ) * ldx ];
        work[ ioffwork + n ] = bi - dx - ex;
        work[ ioffwork ] =
          Math.abs( bi ) + Math.abs( dx ) + Math.abs( ex );
        for ( var i = 2; i <= n - 1; i ++ ) {
          bi = B[ ioffb + i - 1 + ( j - 1 ) * ldb ];
          var cx = e[ ioffe + i - 2 ]
            * X[ ioffx + i - 2 + ( j - 1 ) * ldx ];
          dx = d[ ioffd + i - 1 ]
            * X[ ioffx + i - 1 + ( j - 1 ) * ldx ];
          ex = e[ ioffe + i - 1 ] * X[ ioffx + i + ( j - 1 ) * ldx ];
          work[ ioffwork + n + i - 1 ] = bi - cx - dx - ex;
          work[ ioffwork + i - 1 ] = Math.abs( bi ) + Math.abs( cx )
            + Math.abs( dx ) + Math.abs( ex );
        }
        bi = B[ ioffb + n - 1 + ( j - 1 ) * ldb ];
        cx = e[ ioffe + n - 2 ] * X[ ioffx + n - 2 + ( j - 1 ) * ldx ];
        dx = d[ ioffd + n - 1 ] * X[ ioffx + n - 1 + ( j - 1 ) * ldx ];
        work[ ioffwork + n - 1 + n ] = bi - cx - dx;
        work[ ioffwork + n - 1 ] = Math.abs( bi ) + Math.abs( cx )
          + Math.abs( dx );
      }
      var s = 0.;
      for ( i = 1; i <= n; i ++ ) {
        if ( work[ ioffwork + i - 1 ] > safe2 ) {
          s = Math.max( s, Math.abs( work[ ioffwork + n + i - 1 ] )
            / work[ ioffwork + i - 1 ] );
        } else {
          s = Math.max( s, Math.abs( work[ ioffwork + n + i - 1 ]
            + safe1 ) / ( work[ ioffwork + i - 1 ] + safe1 ) );
        }
      }
      berr[ ioffberr + j - 1 ] = s;
      if ( berr[ ioffberr + j - 1 ] > eps &&
      2. * berr[ ioffberr + j - 1 ] <= lstres && count <= itmax ) {
        LaPack1.dpttrs( n, 1, df, ef, work, n, info, ioffdf, ioffef,
          ioffwork + n );
        Blas1.daxpy( n, 1., work, 1, X, 1, ioffwork + n,
          ioffx + ( j - 1 ) * ldx );
        lstres = berr[ ioffberr + j - 1 ];
        count ++;
        continue
      } else break;
    }
    for ( i = 1; i <= n; i ++ ) {
      if ( work[ ioffwork + i - 1 ] > safe2 ) {
        work[ ioffwork + i - 1 ] = Math.abs(
          work[ ioffwork + n + i - 1 ] )
          + nz * eps * work[ ioffwork + i - 1 ];
      } else {
        work[ ioffwork + i - 1 ] = Math.abs(
          work[ ioffwork + n + i - 1 ] )
          + nz * eps * work[ ioffwork + i - 1 ] + safe1;
      }
    }
    var ix = Blas1.idamax( n, work, 1, ioffwork );
    ferr[ ioffferr + j - 1 ] = work[ ioffwork + ix - 1 ];
    work[ ioffwork ] = 1.;
    for ( i = 2; i <= n; i ++ ) {
      work[ ioffwork + i - 1 ] = 1. + work[ ioffwork + i - 2 ]
        * Math.abs( ef[ ioffef + i - 2 ] );
    }
    work[ ioffwork + n - 1 ] /= df[ ioffdf + n - 1 ];
    for ( i = n - 1; i >= 1; i -- ) {
      work[ ioffwork + i - 1 ] = work[ ioffwork + i - 1 ]
        / df[ ioffdf + i - 1 ]
        + work[ ioffwork + i ] * Math.abs( ef[ ioffef + i - 1 ] );
    }
    ix = Blas1.idamax( n, work, 1, ioffwork );
    ferr[ ioffferr + j - 1 ] *= Math.abs( work[ ioffwork + ix - 1 ] );
    lstres = 0.;
    for ( i = 1; i <= n; i ++ ) {
      lstres = Math.max( lstres,
        Math.abs( X[ ioffx + i - 1 + ( j - 1 ) * ldx ] ) );
    }
    if ( lstres != 0. ) ferr[ ioffferr + j - 1 ] /= lstres;
  }
}
//************************************************************************
LaPack2.dsbgst = function( vect, uplo, n, ka, kb, AB, ldab, BB,
ldbb, X, ldx, work, info ) {
  throw new Error("not programmed: banded matrix");
}
//************************************************************************
LaPack2.dsbtrd = function( vect, uplo, n, kd, AB, ldab, d, e,
Q, ldq, work, info ) {
  throw new Error("not programmed: banded matrix");
}
//************************************************************************
LaPack2.dsptrd = function( uplo, n, AP, d, e, tau, info,
ioffap, ioffd, ioffe, iofftau ) {
  throw new Error("not programmed: packed storage");
}
//************************************************************************
LaPack2.dstein = function( n, d, e, m, w, iblock, isplit, Z,
ldz, work, iwork, ifail, info, ioffd, ioffe, ioffw, ioffiblock,
ioffisplit, ioffz, ioffwork, ioffiwork, ioffifail ) {
  throw new Error("not tested: complicated input");
  var maxits = 5;
  var extra = 2;
  var iseed = new Array( 4 );
  info.setValue( 0 );
  for ( var i = 1; i <= m; i ++ ) ifail[ ioffifail + i - 1 ] = 0;
  if ( n < 0 ) info.setValue( -1 );
  else if ( m < 0 || m > n ) info.setValue( -4 );
  else if ( ldz < Math.max( 1, n ) ) info.setValue( -9 );
  else {
    for ( var j = 2; j <= m; j ++ ) {
      if ( iblock[ ioffiblock + j - 1 ] <
      iblock[ ioffiblock + j - 2 ] ) {
        info.setValue( -6 );
        break;
      }
      if ( iblock[ ioffiblock + j - 1 ] ==
      iblock[ ioffiblock + j - 2 ] &&
      w[ ioffw + j - 1 ] < w[ ioffw + j - 2 ] ) {
        info.setValue( -5 );
        break;
      }
    } // 20
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dstein', - info.getValue() );
    return;
  }
  if ( n == 0 || m == 0 ) return;
  else if ( n == 1 ) {
    Z[ ioffz ] = 1.;
    return;
  }
  var eps = LaPack0.dlamch( 'Precision' );
  for ( i = 1; i <= 4; i ++ ) iseed[ i - 1 ] = 1;
  var indrv1 = 0;
  var indrv2 = indrv1 + n;
  var indrv3 = indrv2 + n;
  var indrv4 = indrv3 + n;
  var indrv5 = indrv4 + n;
  var j1 = 1;
  for ( var nblk = 1; nblk <= iblock[ ioffiblock + m - 1];
  nblk ++ ) {
    var b1 =
      ( nblk == 1 ? 1 : isplit[ ioffisplit + nblk - 2 ] + 1 );
    var bn = isplit[ ioffisplit + nblk - 1 ];
    var blksiz = bn - b1 + 1;
    if ( blksiz != 1 ) {
      var gpind = b1;
      var onenrm = Math.abs( d[ ioffd + b1 - 1 ] )
        + Math.abs( e[ ioffe + b1 - 1 ] );
      onenrm = Math.max( onenrm,
        Math.abs( d[ ioffd + bn - 1 ] )
        + Math.abs( e[ ioffe + bn - 2 ] ) );
      for ( i = b1 + 1; i <= bn - 1; i ++ ) {
        onenrm = Math.max( onenrm, Math.abs( d[ ioffd + i - 1 ] )
          + Math.abs( e[ ioffe + i - 2 ] )
          + Math.abs( e[ ioffe + i - 1 ] ) )
      } // 50
      var ortol = 1.e-3 * onenrm;
      var dtpcrt = Math.sqrt( 1.e-1 / blksiz );
    } // 60
    var jblk = 0;
    var xjm = Number.POSITIVE_INFINITY;
    for ( j = j1; j <= m; j ++ ) {
      if ( iblock[ ioffiblock + j - 1 ] != nblk ) {
        j1 = j;
        break;
      }
      jblk ++;
      var xj = w[ ioffw + j - 1 ];
      if ( blksiz == 1 ) {
        work[ ioffwork + indrv1 ] = 1.;
      } else {
        if ( jblk == 1 ) {
          var eps1 = Math.abs( eps * xj );
          var pertol = 10. * eps1;
          var sep = xj - xjm;
          if ( sep < pertol ) xj = xjm + pertol;
        }
        var its = 0;
        var nrmchk = 0;
        LaPack0.dlarnv( 2, iseed, blksiz, work, 0, ioffwork + indrv1 );
        Blas1.dcopy( blksiz, d, 1, work, 1, ioffd + b1 - 1,
          ioffwork + indrv4 );
        Blas1.dcopy( blksiz - 1, e, 1, work, 1, ioffe + b1 - 1,
          ioffwork + indrv2 + 1 );
        Blas1.dcopy( blksiz - 1, e, 1, work, 1, ioffe + b1 - 1,
          ioffwork + indrv3 );
        var tol = new NumberReference( 0. );
        var iinfo = new IntReference();
        LaPack1.dlagtf( blksiz, work, xj, work, work, tol.getValue(), work,
          iwork, iinfo, ioffwork + indrv4, ioffwork + indrv2 + 1,
          ioffwork + indrv3, indrv5, ioffiwork );
        var goto100 = false;
        while ( true ) { // 70
          its ++;
          if ( its > maxits ) { // 100
            goto100 = true;
          }
          if ( ! goto100 ) {
            scl = blksiz * onenrm * Math.max( eps,
              Math.abs( work[ ioffwork + indrv4 + blksiz - 1 ] ) )
              / Blas1.dasum( blksiz, work, 1, ioffwork + indrv1 );
            Blas1.dscal( blksiz, scl, work, 1, ioffwork + indrv1 );
            LaPack1.dlagts( -1, blksiz, work, work, work, work, iwork,
              work, tol, iinfo, ioffwork + indrv4, ioffwork + indrv2,
              ioffwork + indrv3, ioffwork + indrv5, ioffiwork,
              ioffwork + indrv1 );
            if ( jblk != 1 ) {
              if ( Math.abs( xj - xjm ) > ortol ) gpind = j;
              if ( gpind != j ) {
                for ( i = gpind; i <= j - 1; i ++ ) {
                  var ztr = - Blas1.ddot( blksiz, work, 1, Z,
                    1, ioffwork + indrv1,
                    ioffz + b1 - 1 + ( i - 1 ) * ldz );
                  Blas1.daxpy( blksiz, ztr, Z, 1, work, 1,
                    ioffz + b1 - 1 + ( i - 1 ) * ldz,
                    ioffwork + indrv1 );
                } // 80
              }
            } // 90
            var jmax =
              Blas1.idamax( blksiz, work, 1, ioffwork + indrv1 );
            var nrm =
              Math.abs( work[ ioffwork + indrv1 + jmax - 1 ] );
            if ( nrm < dtpcrt ) continue;
            nrmchk ++;
            if ( nrmchk < extra + 1 ) continue;
            break;
          }
        } // 100
        if ( goto100 ) {
          info.setValue( info.getValue() + 1 );
          ifail[ ioffifail + info.getValue() - 1 ] = j;
        } // 110
        var scl =
          1. / Blas1.dnrm2( blksiz, work, 1, ioffwork + indrv1 );
        jmax = Blas1.idamax( blksiz, work, 1, ioffwork + indrv1 );
        if ( work[ ioffwork + indrv1 + jmax - 1 ] < 0. ) scl = - scl;
        Blas1.dscal( blksiz, scl, work, 1, ioffwork + indrv1 );
      } // 120
      for ( i = 1; i <= n; i ++ ) {
        Z[ ioffz + i - 1 + ( j - 1 ) * ldz ] = 0.;
      }
      for ( i = 1; i <= blksiz; i ++ ) {
        Z[ ioffz + b1 + i - 2 + ( j - 1 ) * ldz ] =
          work[ ioffwork + indrv1 + i - 1 ];
      }
      xjm = xj;
    } // 150
  } // 160
}
LaPack2.zstein = function( n, d, e, m, w, iblock, isplit, Z,
ldz, work, iwork, ifail, info ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dsteqr = function( compz, n, d, e, Z, ldz, work, info,
ioffd, ioffe, ioffz, ioffwork ) {
  var maxit = 30;
  info.setValue( 0 );
  var icompz = -1;
  if ( compz.charAt(0).toUpperCase() === 'N' ) {
    icompz = 0;
  } else if ( compz.charAt(0).toUpperCase() == 'V' ) {
    icompz = 1;
  } else if ( compz.charAt(0).toUpperCase() == 'I' ) {
    icompz = 2;
  }
  if ( icompz < 0 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  else if ( ldz < 1 || ( icompz > 0 && ldz < Math.max( 1, n ) ) ) {
    info.setValue( -6 );
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dsteqr', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  if ( n == 1 ) {
    if ( icompz == 2 ) Z[ ioffz ] = 1.;
    return;
  }
  var eps = LaPack0.dlamch( 'E' );
  var eps2 = eps * eps;
  var safmin = LaPack0.dlamch( 'S' );
  var safmax = 1. / safmin;
  var ssfmax = Math.sqrt( safmax ) / 3.;
  var ssfmin = Math.sqrt( safmin ) / eps2;
  if ( icompz == 2 ) {
    LaPack0.dlaset( 'Full', n, n, 0., 1., Z, ldz, ioffz );
  }
  var nmaxit = n * maxit;
  var jtot = 0;
  var l1 = 1;
  var nm1 = n - 1;
  var goto160 = false;
  while ( true ) { // 10 --> 160
    if ( l1 > n ) {
      goto160 = true;
      break;
    }
    if ( l1 > 1 ) e[ ioffe + l1 - 2 ] = 0.;
    var goto30 = false;
    if ( l1 <= nm1 ) {
      for ( var m = l1; m <= nm1; m ++ ) {
        var tst = Math.abs( e[ ioffe + m - 1 ] );
        if ( tst == 0. ) {
          goto30 = true;
          break;
        }
        if ( tst <= ( Math.sqrt( Math.abs( d[ ioffd + m - 1 ] ) )
        * Math.sqrt( Math.abs( d[ ioffd + m ] ) ) ) * eps ) {
          e[ ioffe + m - 1 ] = 0.;
          goto30 = true;
          break;
        }
      } // 20
    }
    if ( ! goto30) m = n;
    var l = l1;
    var lsv = l
    var lend = m;
    var lendsv = lend;
    l1 = m + 1;
    if ( lend == l ) continue;
    var anorm = LaPack1.dlanst( 'M', lend - l + 1, d, e,
      ioffd + l - 1, ioffe + l - 1 );
    var iscale = 0;
    if ( anorm == 0. ) continue;
    if ( anorm > ssfmax ) {
      iscale = 1;
      LaPack1.dlascl( 'G', 0, 0, anorm, ssfmax, lend - l + 1, 1, d, n,
        info, ioffd + l - 1 );
      LaPack1.dlascl( 'G', 0, 0, anorm, ssfmax, lend - l, 1, e, n,
        info, ioffe + l - 1 );
    } else if ( anorm < ssfmin ) {
      iscale = 2;
      LaPack1.dlascl( 'G', 0, 0, anorm, ssfmin, lend - l + 1, 1, d, n,
        info, ioffd + l - 1 );
      LaPack1.dlascl( 'G', 0, 0, anorm, ssfmin, lend - l, 1, e, n,
        info, ioffe + l - 1 );
    }
    if ( Math.abs( d[ ioffd + lend - 1 ] ) <
    Math.abs( d[ ioffd + l - 1 ] ) ) {
      lend = lsv;
      l = lendsv;
    }
    var rt1 = new NumberReference();
    var rt2 = new NumberReference();
    var c = new NumberReference();
    var s = new NumberReference();
    if ( lend > l ) {
      while ( true ) { // 40
        if ( l != lend ) {
          var lendm1 = lend - 1;
          var goto60 = false;
          for ( m = l; m <= lendm1; m ++ ) {
            tst = Math.pow( Math.abs( e[ ioffe + m - 1 ] ), 2 );
            if ( tst <= ( eps2 * Math.abs( d[ ioffd + m - 1 ] ) )
            * Math.abs( d[ ioffd + m ] ) + safmin ) {
              goto60 = true;
              break;
            }
          }
        }
        if ( ! goto60 ) m = lend;
        if ( m < lend ) e[ ioffe + m - 1 ] = 0.;
        var p = d[ ioffd + l - 1 ];
        if ( m != l ) {
          if ( m == l + 1 ) {
            if ( icompz > 0 ) {
              LaPack0.dlaev2( d[ ioffd + l - 1 ], e[ ioffe + l - 1 ],
                d[ ioffd + l ], rt1, rt2, c, s );
              work[ ioffwork + l - 1 ] = c.getValue();
              work[ ioffwork + n - 2 + l ] = s.getValue();
              LaPack0.dlasr( 'R', 'V', 'B', n, 2, work, work, Z, ldz,
                ioffwork + l - 1, ioffwork + n - 2 + l,
                ioffz + ( l - 1 ) * ldz );
            } else {
              LaPack0.dlae2( d[ ioffd + l - 1 ], e[ ioffe + l - 1 ],
                d[ ioffd + l ], rt1, rt2 );
            }
            d[ ioffd + l - 1 ] = rt1.getValue();
            d[ ioffd + l ] = rt2.getValue();
            e[ ioffe + l - 1 ] = 0.;
            l += 2;
            if ( l <= lend ) continue;
            break; // ie, go to 140
          }
          if ( jtot == nmaxit ) break;
          jtot ++;
          var g =
            ( d[ ioffd + l ] - p ) / ( 2. * e[ ioffe + l - 1 ] );
          var r =
            new NumberReference( LaPack0.dlapy2( g, 1. ) );
          g = d[ ioffd + m - 1 ] - p + ( e[ ioffe + l - 1 ]
            / ( g + ( g >= 0 ? r.getValue() : -r.getValue() ) ) );
          s.setValue( 1. );
          c.setValue( 1. );
          p = 0.;
          var mm1 = m - 1;
          for ( var i = mm1; i >= l; i -- ) {
            var f = s.getValue() * e[ ioffe + i - 1 ];
            var b = c.getValue() * e[ ioffe + i - 1 ];
            LaPack1.dlartg( g, f, c, s, r );
            if ( i != m - 1 ) e[ ioffe + i ] = r.getValue();
            g = d[ ioffd + i ] - p;
            r.setValue( ( d[ ioffd + i - 1 ] - g ) * s.getValue()
              + 2. * c.getValue() * b );
            p = s.getValue() * r.getValue();
            d[ ioffd + i ] = g + p;
            g = c.getValue() * r.getValue() - b;
            if ( icompz > 0 ) {
              work[ ioffwork + i - 1 ] = c.getValue();
              work[ ioffwork + n - 2 + i ] = - s.getValue();
            }
          } // 70
          if ( icompz > 0 ) {
            var mm = m - l + 1;
            LaPack0.dlasr( 'R', 'V', 'B', n, mm, work, work, Z, ldz,
              ioffwork + l - 1, ioffwork + n - 2 + l,
              ioffz + ( l - 1 ) * ldz );
          }
          d[ ioffd + l - 1 ] -= p;
          e[ ioffe + l - 1 ] = g;
          continue;
        } // 80
        d[ ioffd + l - 1 ] = p;
        l ++;
        if ( l <= lend) continue;
        else break; // ie, go to 140
      }
    } else {
      while ( true ) { // 90 --> 140
        var goto110 = false;
        if ( l != lend) {
          var lendp1 = lend + 1;
          for ( m = l; l >= lendp1; l -- ) {
            tst = Math.pow( Math.abs( e[ ioffe + m - 2 ] ), 2 );
            if ( tst <= ( eps2 * Math.abs( d[ ioffd + m - 1 ] ) )
            * Math.abs( d[ ioffd + m - 2 ] ) + safmin ) {
              goto110 = true;
              break;
            }
          }
        }
        if ( ! goto110 ) m = lend;
        if ( m > lend ) e[ ioffe + m - 2 ] = 0.;
        p = d[ ioffd + l - 1 ];
        if ( m != l ) {
          if ( m == l - 1 ) {
            if ( icompz > 0 ) {
              LaPack0.dlaev2( d[ ioffd + l - 2], e[ ioffe + l - 2 ],
                d[ ioffd + l - 1 ], rt1, rt2, c, s );
              work[ ioffwork + m - 1 ] = c.getValue();
              work[ ioffwork + n - 2 + m ] = s.getValue();
              LaPack0.dlasr( 'R', 'V', 'F', n, 2, work, work, Z, ldz,
                ioffwork + m - 1, ioffwork + n - 2 + m,
                ioffz + ( l - 2 ) * ldz );
            } else {
              LaPack0.dlae2( d[ ioffd + l - 2 ], e[ ioffe + l - 2 ],
                d[ ioffd + l - 1 ], rt1, rt2 );
            }
            d[ ioffd + l - 2 ] = rt1.getValue();
            d[ ioffd + l - 1 ] = rt2.getValue();
            e[ ioffe + l - 2 ] = 0.;
            l -= 2;
            if ( l >= lend ) continue; // ie, go to 90
            break; // ie, go to 140
          }
          if ( jtot == nmaxit ) break;
          jtot ++;
          g = ( d[ ioffd + l - 2 ] - p ) / ( 2. * e[ ioffe + l - 2 ] );
          r.setValue( LaPack0.dlapy2( g, 1. ) );
          g = d[ ioffd + m - 1 ] - p + ( e[ ioffe + l - 2 ] / ( g
            + ( g >= 0 ? r.getValue() : -r.getValue() ) ) );
          s.setValue( 1. );
          c.setValue( 1. );
          p = 0.;
          var lm1 = l - 1;
          for ( i = m; i <= lm1; i ++ ) {
            f = s.getValue() * e[ ioffe + i - 1 ];
            b = c.getValue() * e[ ioffe + i - 1 ];
            LaPack1.dlartg( g, f, c, s, r );
            if ( i != m ) e[ ioffe + i - 2 ] = r.getValue();
            g = d[ ioffd + i  - 1 ] - p;
            r.setValue(
              ( d[ ioffd + i ] - g ) * s.getValue() + 2. * c.getValue() * b );
            p = s.getValue() * r.getValue();
            d[ ioffd + i  - 1 ] = g + p;
            g = c.getValue() * r.getValue() - b;
            if ( icompz > 0 ) {
              work[ ioffwork + i - 1 ] = c.getValue();
              work[ ioffwork + n - 2 + i ] = s.getValue();
            }
          } // 120
          if ( icompz > 0 ) {
            mm = l - m + 1;
            LaPack0.dlasr( 'R', 'V', 'F', n, mm, work, work, Z, ldz,
              ioffwork + m - 1, ioffwork + n - 2 + m,
              ioffz + ( m - 1 ) * ldz );
          }
          d[ ioffd + l - 1 ] -= p;
          e[ ioffe + lm1 - 1 ] = g;
          continue; // ie, go to 90
        } // 130
        d[ ioffd + l - 1 ] = p;
        l --;
        if ( l >- lend ) continue;
        break;
      }
    } // 140
    if ( iscale == 1 ) {
      LaPack1.dlascl( 'G', 0, 0, ssfmax, anorm, lendsv - lsv + 1, 1,
        d, n, info, ioffd + lsv - 1 );
      LaPack1.dlascl( 'G', 0, 0, ssfmax, anorm, lendsv - lsv, 1,
        e, n, info, ioffe + lsv - 1 );
    } else if ( iscale == 2 ) {
      LaPack1.dlascl( 'G', 0, 0, ssfmin, anorm, lendsv - lsv + 1, 1,
        d, n, info, ioffd + lsv - 1 );
      LaPack1.dlascl( 'G', 0, 0, ssfmin, anorm, lendsv - lsv, 1,
        e, n, info, ioffe + lsv - 1 );
    }
    if ( jtot < nmaxit ) continue; // ie, go to 10
    else break;
  }
  if ( ! goto160 ) {
    for ( i = 1; i <= n - 1; i ++ ) {
      if ( e[ ioffe + i - 1 ] != 0. ) info.getValue() = info.getValue() + 1;
    } // 150
    return;
  } // 160
  if ( icompz == 0 ) LaPack0.dlasrt( 'I', n, d, info, ioffd );
  else {
    for ( var ii = 2; ii <= n; ii ++ ) {
      i = ii - 1;
      var k = i;
      p = d[ ioffd + i - 1 ];
      for ( var j = ii; j <= n; j ++ ) {
        if ( d[ ioffd + j - 1 ] < p ) {
          k = j;
          p = d[ ioffd + j - 1 ];
        }
      } // 170
      if ( k != i ) {
        d[ ioffd + k - 1 ] = d[ ioffd + i - 1 ];
        d[ ioffd + i - 1 ] = p;
        Blas1.dswap( n, Z, 1, Z, 1, ioffz + ( i - 1 ) * ldz,
          ioffz + ( k - 1 ) * ldz );
      }
    } // 180
  } // 190
}
LaPack2.zsteqr = function( compz, n, d, e, Z, ldz, work, info )
{
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack2.dsterf = function( n, d, e, info, ioffd, ioffe ) {
  var maxit = 30;
  info.setValue( 0 );
  if ( n < 0 ) {
    info.setValue( -1 );
    Blas2.xerbla( 'dsterf', - info.getValue() );
    return;
  }
  if ( n <= 1 ) return;
  var eps = LaPack0.dlamch( 'E' );
  var eps2 = eps * eps;
  var safmin = LaPack0.dlamch( 'S' );
  var safmax = 1. / safmin;
  var ssfmax = Math.sqrt( safmax ) / 3.;
  var ssfmin = Math.sqrt( safmin ) / eps2;
  var nmaxit = n * maxit;
  var sigma = 0.;
  var jtot = 0;
  var l1 = 1;
  var nm1 = n - 1;
  while ( true ) { // 10
    if ( l1 > n ) break;
    if ( l1 > 1 ) e[ ioffe + l1 - 2 ] = 0.;
    var goto30 = false;
    if ( l1 <= nm1 ) {
      for ( var m = l1; m <= nm1; m ++ ) {
        var tst = Math.abs( e[ ioffe + m - 1 ] );
        if ( tst == 0. ) {
          goto30 = true;
          break;
        }
        if ( tst <= ( Math.sqrt( Math.abs( d[ ioffd + m - 1 ] ) )
        * Math.sqrt( Math.abs( d[ ioffd + m ] ) ) ) * eps ) {
          e[ ioffe + m - 1 ] = 0.;
          goto30 = true;
          break;
        }
      }
    }
    if ( ! goto30 ) m = n;
    var l = l1;
    var lsv = l;
    var lend = m;
    var lendsv = lend;
    l1 = m + 1;
    if ( lend == l ) continue;
    var anorm = LaPack1.dlanst( 'I', lend - l + 1, d, e,
      ioffd + l - 1, ioffe + l - 1 );
    var iscale = 0;
    if ( anorm > ssfmax ) {
      iscale = 1;
      LaPack1.dlascl( 'G', 0, 0, anorm, ssfmax, lend - l + 1, 1, d,
        n, info, ioffd + l - 1 );
      LaPack1.dlascl( 'G', 0, 0, anorm, ssfmax, lend - l, 1, e,
        n, info, ioffe + l - 1 );
    } else if ( anorm < ssfmin ) {
      iscale = 2;
      LaPack1.dlascl( 'G', 0, 0, anorm, ssfmin, lend - l + 1, 1, d,
        n, info, ioffd + l - 1 );
      LaPack1.dlascl( 'G', 0, 0, anorm, ssfmin, lend - l, 1, e,
        n, info, ioffe + l - 1 );
    }
    for ( var i = l; i <= lend - 1; i ++ ) {
      e[ ioffe + i - 1 ] = Math.pow( e[ ioffe + i - 1 ] , 2 );
    }
    if ( Math.abs( d[ ioffd + lend - 1 ] ) <
    Math.abs( d[ ioffd + l - 1 ] ) ) {
      lend = lsv;
      l = lendsv;
    }
    if ( lend >= l ) { // QL
      while ( true ) { // 50
        var goto70 = false;
        if ( l != lend ) {
          var lendm1 = lend - 1;
          for ( m = l; m <= lendm1; m ++ ) {
            tst = Math.abs( e[ ioffe + m - 1 ] );
            if ( tst <= eps2 * Math.abs( d[ ioffd + m - 1 ]
            * d[ ioffd + m ] ) ) {
              goto70 = true;
              break;
            }
          }
        }
        if ( ! goto70 ) m = lend;
        if ( m < lend ) e[ ioffe + m - 1 ] = 0.;
        var p = d[ ioffd + l - 1 ];
        if ( m != l ) {
          if ( m == l + 1 ) {
            rte = Math.sqrt( e[ ioffe + l - 1 ] );
            var rt1 = new NumberReference();
            var rt2 = new NumberReference();
            LaPack0.dlae2( d[ ioffd + l - 1 ], rte, d[ ioffd + l ],
              rt1, rt2 );
            d[ ioffd + l - 1 ] = rt1.getValue();
            d[ ioffd + l ] = rt2.getValue();
            e[ ioffe + l - 1 ] = 0.;
            l += 2;
            if ( l <= lend ) continue;
            break;
          }
          if ( jtot == nmaxit ) break;
          jtot ++;
          var rte = Math.sqrt( e[ ioffe + l - 1 ] );
          sigma = ( d[ ioffd + l ] - p ) / ( 2. * rte );
          var r = LaPack0.dlapy2( sigma, 1. );
          sigma = p - ( rte / ( sigma + ( sigma >= 0. ? r : -r ) ) );
          var c = 1.;
          var s = 0.;
          var gamma = d[ ioffd + m - 1 ] - sigma;
          p = gamma * gamma;
          var mm1 = m - 1;
          for ( i = mm1; i >= l; i -- ) {
            var bb = e[ ioffe + i - 1 ];
            r = p + bb;
            if ( i != m - 1 ) e[ ioffe + i ] = s * r;
            var oldc = c;
            c = p / r;
            s = bb / r;
            var oldgam = gamma;
            var alpha = d[ ioffd + i - 1 ];
            gamma = c * ( alpha - sigma ) - s * oldgam;
            d[ ioffd + i ] = oldgam + ( alpha - gamma );
            p = ( c != 0. ? ( gamma * gamma ) / c : oldc * bb );
          }
          e[ ioffe + l - 1 ] = s * p;
          d[ ioffd + l - 1 ] = sigma + gamma;
          continue;
        } // 90
        d[ ioffd + l - 1 ] = p;
        l ++;
        if ( l <= lend ) continue;
        break;
      }
    } else { // QR
      while ( true ) { // 100
        var goto120 = false;
        if ( l != lend ) {
          var lendp1 = lend + 1;
          for ( m = l; m >= lendp1; m -- ) {
            tst = Math.abs( e[ ioffe + m - 2 ] );
            if ( tst <= eps2 * Math.abs( d[ ioffd + m - 1 ]
            * d[ ioffd + m - 2 ] ) ) {
              goto120 = true;
              break;
            }
          }
        }
        if ( ! goto120 ) m = lend;
        if ( m > lend ) e[ ioffe + m - 2 ] = 0.;
        p = d[ ioffd + l - 1 ];
        if ( m != l ) {
          if ( m == l - 1 ) {
            rte = Math.sqrt( e[ ioffe + l - 1 ] );
            LaPack0.dlae2( d[ ioffd + l - 1 ], rte,
              d[ ioffd + l - 2 ], rt1, rt2 );
            d[ ioffd + l - 1 ] = rt1.getValue();
            d[ ioffd + l - 2 ] = rt2.getValue();
            e[ ioffe + l - 2 ] = 0.;
            l -= 2;
            if ( l <= lend ) continue;
            break;
          }
          if ( jtot == nmaxit ) break;
          jtot ++;
          rte = Math.sqrt( e[ ioffe + l - 2 ] );
          sigma = ( d[ ioffd + l - 2 ] - p ) / ( 2. * rte );
          r = LaPack0.dlapy2( sigma, 1. );
          sigma = p - ( rte / ( sigma + ( sigma >= 0. ? r : -r ) ) );
          c = 1.;
          s = 0.;
          gamma = d[ ioffd + m - 1 ] - sigma;
          p = gamma * gamma;
          var lm1 = l - 1;
          for ( i = m; i <= lm1; i ++ ) {
            bb = e[ ioffe + i - 1 ];
            r = p + bb;
            if ( i != m ) e[ ioffe + i - 2 ] = s * r;
            oldc = c;
            c = p / r;
            s = bb / r;
            oldgam = gamma;
            alpha = d[ ioffd + i ];
            gamma = c * ( alpha - sigma ) - s * oldgam;
            d[ ioffd + i - 1 ] = oldgam + ( alpha - gamma );
            p = ( c != 0. ? ( gamma * gamma ) / c : oldc * bb );
          }
          e[ ioffe + lm1 - 1 ] = s * p;
          d[ ioffd + l - 1 ] = sigma + gamma;
          continue;
        } // 140
        d[ ioffd + l - 1 ] = p;
        l --;
        if ( l >= lend ) continue;
        break;
      }
    } // 150
    if ( iscale == 1 ) {
      LaPack1.dlascl( 'G', 0, 0, ssfmax, anorm, lendsv - lsv + 1, 1,
        d, n, info, ioffd + lsv - 1 );
    }
    if ( iscale == 2 ) {
      LaPack1.dlascl( 'G', 0, 0, ssfmin, anorm, lendsv - lsv + 1, 1,
        d, n, info, ioffd + lsv - 1 );
    }
    if ( jtot == nmaxit ) {
      for ( i = 1; i <= n - 1; i ++ ) {
        if ( e[ ioffe + i - 1 ] != 0. ) info.getValue() = info.getValue() + 1;
      }
      return;
    }
  } // 170
  LaPack0.dlasrt( 'I', n, d, info, ioffd );
}
//************************************************************************
LaPack2.dsytd2 = function( uplo, n, A, lda, d, e, tau, info,
ioffa, ioffd, ioffe, iofftau ) {
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -4 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dsytd2', - info.getValue() );
    return;
  }
  if ( n <= 0 ) return;
  var alpha = new NumberReference();
  var taui = new NumberReference();
  if ( upper ) { 
    for ( var i = n - 1; i >= 1; i -- ) {
      alpha.setValue( A[ ioffa + i - 1 + i * lda ] );
      LaPack1.dlarfg( i, alpha, A, 1, taui, ioffa + i * lda );
      A[ ioffa + i - 1 + i * lda ] = alpha.getValue();
      e[ ioffe + i - 1 ] = A[ ioffa + i - 1 + i * lda ];
      if ( taui.getValue() != 0. ) {
        A[ ioffa + i - 1 + i * lda ] = 1.;
        Blas2.dsymv( uplo, i, taui.getValue(), A, lda, A, 1, 0., tau, 1,
          ioffa, ioffa + i * lda, iofftau );
        alpha.setValue( - 0.5 * taui.getValue()
          * Blas1.ddot( i, tau, 1, A, 1, iofftau, ioffa + i * lda ) );
        Blas1.daxpy( i, alpha.getValue(), A, 1, tau, 1, ioffa + i * lda,
          iofftau );
        Blas2.dsyr2( uplo, i, -1., A, 1, tau, 1, A, lda,
          ioffa + i * lda, iofftau, ioffa );
        A[ ioffa + i - 1 + i * lda ] = e[ ioffe + i - 1 ];
      }
      d[ ioffd + i ] = A[ ioffa + i + i * lda ];
      tau[ iofftau + i - 1 ] = taui.getValue();
    }
    d[ ioffd ] = A[ ioffa ];
  } else {
    for ( i = 1; i <= n - 1; i ++ ) {
      alpha.setValue( A[ ioffa + i + ( i - 1 ) * lda ] );
      LaPack1.dlarfg( n - i, alpha, A, 1, taui,
        ioffa + Math.min( i + 2, n ) - 1 + ( i - 1 ) * lda );
      A[ ioffa + i + ( i - 1 ) * lda ] = alpha.getValue();
      e[ ioffe + i - 1 ] = A[ ioffa + i + ( i - 1 ) * lda ];
      if ( taui.getValue() != 0. ) {
        A[ ioffa + i + ( i - 1 ) * lda ] = 1.;
        Blas2.dsymv( uplo, n - i, taui.getValue(), A, lda, A, 1, 0., tau, 1,
          ioffa + i + i * lda, ioffa + i + ( i - 1 ) * lda,
          iofftau + i - 1 );
        alpha.setValue( - 0.5 * taui.getValue()
          * Blas1.ddot( n - i, tau, 1, A, 1, iofftau + i - 1,
          ioffa + i + ( i - 1 ) * lda ) );
        Blas1.daxpy( n - i, alpha.getValue(), A, 1, tau, 1,
          ioffa + i + ( i - 1 ) * lda, iofftau + i - 1 );
        Blas2.dsyr2( uplo, n - i, -1., A, 1, tau, 1, A, lda,
          ioffa + i + ( i - 1 ) * lda, iofftau + i - 1,
          ioffa + i + i * lda );
        A[ ioffa + i + ( i - 1 ) * lda ] = e[ ioffe + i - 1 ];
      }
      d[ ioffd + i - 1 ] = A[ ioffa + i - 1 + ( i - 1 ) * lda ];
      tau[ iofftau + i - 1 ] = taui.getValue();
    }
    d[ ioffd + n - 1 ] = A[ ioffa + n - 1 + ( n - 1 ) * lda ];
  }
}
//************************************************************************
LaPack2.dsytri2x = function( uplo, n, A, lda, ipiv, work, nb,
info, ioffa, ioffipiv, ioffwork ) {
  throw new Error("not programmed: TRF format");
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, n ) ) info.getValue() = -4;
  if ( info.getValue() != 0 ) {
    Blas1.xerbla( 'dsytri2x', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
}
//************************************************************************
LaPack2.dtftri = function( transr, uplo, diag, n, A, info,
ioffa ) {
  throw new Error("not programmed: RFP format");
}
//************************************************************************
LaPack2.dtgevc = function( side, howmny, select, n, A, lda, B,
ldb, vl, ldvl, vr, ldvr, mm, m, work, info ) {
  throw new Error("not programmed: generalized eigenvalue problem");
}
LaPack2.ztgevc = function( side, howmny, select, n, A, lda, B,
ldb, vl, ldvl, vr, ldvr, mm, m, work, info ) {
  throw new Error("not programmed: generalized eigenvalue problem");
}
//************************************************************************
LaPack2.dtpcon = function( norm, uplo, diag, n, AP, rcond,
work, iwork, info ) {
  throw new Error("not programmed: packed matrix");
}
LaPack2.ztpcon = function( norm, uplo, diag, n, AP, rcond,
work, iwork, info ) {
  throw new Error("not programmed: complex packed matrix");
}
//************************************************************************
LaPack2.dtpqrt2 = function( m, n, l, A, lda, B, ldb, T, ldt,
info, ioffa, ioffb, iofft ) {
  throw new Error("not programmed: WY format");
}
//************************************************************************
LaPack2.dtrevc = function( side, howmny, select, n, T, ldt,
VL, ldvl, VR, ldvr, mm, m, work, info, ioffselect, iofft, ioffvl, ioffvr,
ioffwork ) {
  throw new Error("not programmed");
}
LaPack2.ztrevc = function( side, howmny, select, n, T, ldt, VL,
ldvl, VR, ldvr, mm, m, work, info ) {
  throw new Error("not programmed");
}
//************************************************************************
LaPack2.dtrsyl = function( trana, tranb, isgn, m, n, A, lda, B,
ldb, c, ldc, scale, info, ioffa, ioffb, ioffc ) {
  throw new Error("not programmed");
}
LaPack2.ztrsyl = function( trana, tranb, isgn, m, n, A, lda, B,
ldb, c, ldc, scale, info ) {
  throw new Error("not programmed");
}
//************************************************************************
/*  deprecated use dtzrzf
LaPack2.dtzrqf = function( m, n, A, lda, tau, info, ioffa,
iofftau ) {
  info.setValue( 0 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < m ) info.setValue( -2 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -4 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dtzrqf', - info.getValue() );
    return;
  }
  if ( m == 0 ) return;
  if ( m == n ) {
    for ( var i = 1; i <= n; i ++ ) tau[ iofftau + i - 1 ] = 0.;
  } else {
    var m1 = Math.min( m + 1, n );
    for ( var k = m; k >= 1; k -- ) {
      var alpha =
        new NumberReference( A[ ioffa + k - 1 + ( k - 1 ) * lda ] );
      var tauref =
        new NumberReference( tau[ iofftau + k - 1 ] );
      LaPack1.dlarfg( n - m + 1, alpha, A, lda, tauref,
        ioffa + k - 1 + ( m1 - 1 ) * lda );
      A[ ioffa + k - 1 + ( k - 1 ) * lda ] = alpha.getValue();
      tau[ iofftau + k - 1 ] = tauref.getValue();
      if ( tau[ iofftau + k - 1 ] != 0. && k > 1 ) {
        Blas1.dcopy( k - 1, A, 1, tau, 1, ioffa + ( k - 1 ) * lda,
          iofftau );
        Blas2.dgemv( 'No transpose', k - 1, n - m, 1., A, lda, A, lda,
          1., tau, 1, ioffa + ( m1 - 1 ) * lda,
          ioffa + k - 1 + ( m1 - 1 ) * lda, iofftau );
        Blas1.daxpy( k - 1, - tau[ k - 1 ], tau, 1, A, 1, iofftau,
          ioffa + ( k - 1 ) * lda );
        Blas2.dger( k - 1, n - m, -tau[ k - 1 ], tau, 1, A, lda,
          A, lda, iofftau, ioffa + k - 1 + ( m1 - 1 ) * lda
          , ioffa + ( m1 - 1 ) * lda );
      }
    }
  }
}
LaPack2.ztzrqf = function( m, n, A, lda, tau, info ) {
  throw new Error("not programmed: complex matrix");
}
*/
function testLaPack2() {
  test_dgebd2();
//test_dgecon();
//test_dgehd2();
//test_dgelq2();
//test_dgeql2();
//test_dgeqr2();
  test_dgeqr2p();
//test_dgerfs();
//test_dgerq2();
//test_dgetrf();
//test_dgtcon();
//test_dgtrfs();
//test_dlaed4();
//test_dlaein();
//test_dlaexc();
//test_dlapll();
//test_dlarfx();
//test_dlatrd();
//test_dopmtr();
//test_dorg2l();
//test_dorg2r();
//test_dorgl2();
//test_dorgr2();
//test_dorm2l();
//test_dorm2r();
//test_dorml2();
//test_dormr2();
//test_dormr3();
//test_dpocon();
//test_dptrfs();
//test_dsteqr();
//test_dsterf();
//test_dsytd2();
//test_dsytrf();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgebd2() {
  document.getElementById("debug_textarea").value +=
    "testing dgebd2 *************" + "\n";
  var m = 4;
  var n = 3;
  var lda = 5;
  var ldb = 4;
  var ioffa = 1;
  var ioffb = 2;
  var A = new Array( ioffa + lda * n );
  var B = new Array( ioffb + ldb * m );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + j );
      B[ ioffb + j + i * ldb ] = 1. / Number( 1 + i + j );
    }
  }
  var ioffd = 3;
  var ioffe = 4;
  var iofftauq = 5;
  var iofftaup = 6;
  var ioffwork = 7;
  var d = new Array( ioffd + Math.min( m, n ) );
  var e = new Array( ioffe + Math.min( m, n ) - 1 );
  var tauq = new Array( iofftauq + Math.min( m, n ) );
  var taup = new Array( iofftaup + Math.min( m, n ) );
  var work = new Array( ioffwork + Math.max( m, n ) );
  var info = new IntReference();
  LaPack2.dgebd2(m,n,A,lda,d,e,tauq,taup,work,info,
    ioffa,ioffd,ioffe,iofftauq,iofftaup,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dgebd2(m>n): info = " + info.value + "\n";
  document.getElementById("debug_textarea").value +=
    "A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 3 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 3 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 3 + 2 * lda ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "d = "
    + d[ ioffd + 0 ] + " "
    + d[ ioffd + 1 ] + " "
    + d[ ioffd + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "e = "
    + e[ ioffe + 0 ] + " "
    + e[ ioffe + 1 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "tauq = "
    + tauq[ iofftauq + 0 ] + " "
    + tauq[ iofftauq + 1 ] + " "
    + tauq[ iofftauq + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "taup = "
    + taup[ iofftaup + 0 ] + " "
    + taup[ iofftaup + 1 ] + " "
    + taup[ iofftaup + 2 ]  + "\n";

  LaPack2.dgebd2(n,m,B,ldb,d,e,tauq,taup,work,info,
    ioffb,ioffd,ioffe,iofftauq,iofftaup,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dgebd2(m<n): info = " + info.value + "\n";
  document.getElementById("debug_textarea").value +=
    "B = "
    + B[ ioffb + 0 + 0 * ldb ] + " "
    + B[ ioffb + 1 + 0 * ldb ] + " "
    + B[ ioffb + 2 + 0 * ldb ] + " "
    + B[ ioffb + 0 + 1 * ldb ] + " "
    + B[ ioffb + 1 + 1 * ldb ] + " "
    + B[ ioffb + 2 + 1 * ldb ] + " "
    + B[ ioffb + 0 + 2 * ldb ] + " "
    + B[ ioffb + 1 + 2 * ldb ] + " "
    + B[ ioffb + 2 + 2 * ldb ] + " "
    + B[ ioffb + 0 + 3 * ldb ] + " "
    + B[ ioffb + 1 + 3 * ldb ] + " "
    + B[ ioffb + 2 + 3 * ldb ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "d = "
    + d[ ioffd + 0 ] + " "
    + d[ ioffd + 1 ] + " "
    + d[ ioffd + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "e = "
    + e[ ioffe + 0 ] + " "
    + e[ ioffe + 1 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "tauq = "
    + tauq[ iofftauq + 0 ] + " "
    + tauq[ iofftauq + 1 ] + " "
    + tauq[ iofftauq + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "taup = "
    + taup[ iofftaup + 0 ] + " "
    + taup[ iofftaup + 1 ] + " "
    + taup[ iofftaup + 2 ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgecon() {
  document.getElementById("debug_textarea").value +=
    "testing dgecon *************" + "\n";
  var n = 3;
  var lda = 4;
  var anorm = 4.;
  var ioffa = 1;
  var ioffipiv = 2;
  var ioffwork = 3;
  var ioffiwork = 4;
  var A = new Array( ioffa + lda * n );
  var ipiv = new Array( ioffipiv + n );
  var iwork = new Array( ioffiwork + n );
  var work = new Array( ioffwork + 4*n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + 2 * j );
    }
  }
  var info = new IntReference();
  LaPack2.dgetrf(n,n,A,lda,ipiv,info,ioffa,ioffipiv);
  var rcond = new NumberReference();
  LaPack2.dgecon('O',n,A,lda,anorm,rcond,work,iwork,info,
    ioffa,ioffwork,ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dgecon('O',...): info,rcond = " + info.value + " "
    + rcond.value + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + 2 * j );
    }
  }
  LaPack2.dgetrf(n,n,A,lda,ipiv,info,ioffa,ioffipiv);
  LaPack2.dgecon('I',n,A,lda,anorm,rcond,work,iwork,info,
    ioffa,ioffwork,ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dgecon('I',...): info,rcond = " + info.value + " "
    + rcond.value + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgehd2() {
  document.getElementById("debug_textarea").value +=
    "testing dgehd2 *************" + "\n";
  var n = 4;
  var lda = 5;
  var ilo = 1;
  var ihi = 4;
  var ioffa = 1;
  var iofftau = 2;
  var ioffwork = 3;
  var A = new Array( ioffa + lda * n );
  var tau = new Array( iofftau + n - 1 );
  var work = new Array( ioffwork + n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + 2 * j );
    }
  }
  var info = new IntReference();
  LaPack2.dgehd2(n,ilo,ihi,A,lda,tau,work,info,ioffa,iofftau,ioffwork);
  document.getElementById("debug_textarea").value +=
    "info,tau = " + info.value + " "
    + tau[ iofftau + 0 ] + " "
    + tau[ iofftau + 1 ] + " "
    + tau[ iofftau + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 3 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 3 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 3 + 2 * lda ] + " "
    + A[ ioffa + 0 + 3 * lda ] + " "
    + A[ ioffa + 1 + 3 * lda ] + " "
    + A[ ioffa + 2 + 3 * lda ] + " "
    + A[ ioffa + 3 + 3 * lda ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgelq2() {
  document.getElementById("debug_textarea").value +=
    "testing dgelq2 *************" + "\n";
  var n = 4;
  var m = 5;
  var lda = 6;
  var ldb = 5;
  var ioffa = 1;
  var ioffb = 4;
  var iofftau = 2;
  var ioffwork = 3;
  var A = new Array( ioffa + lda * n );
  var B = new Array( ioffb + ldb * m );
  var tau = new Array( iofftau + Math.min( m, n ) );
  var work = new Array( ioffwork + Math.max( m, n ) );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + 2 * j );
      B[ ioffb + j + i * ldb ] = 1. / Number( 1 + i + 2 * j );
    }
  }
  var info = new IntReference();
  LaPack2.dgelq2(m,n,A,lda,tau,work,info,ioffa,iofftau,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dgelq2(m>n): info,tau = " + info.value + " "
    + tau[ iofftau + 0 ] + " "
    + tau[ iofftau + 1 ] + " "
    + tau[ iofftau + 2 ] + " "
    + tau[ iofftau + 3 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 3 + 0 * lda ] + " "
    + A[ ioffa + 4 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 3 + 1 * lda ] + " "
    + A[ ioffa + 4 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 3 + 2 * lda ] + " "
    + A[ ioffa + 4 + 2 * lda ] + " "
    + A[ ioffa + 0 + 3 * lda ] + " "
    + A[ ioffa + 1 + 3 * lda ] + " "
    + A[ ioffa + 2 + 3 * lda ] + " "
    + A[ ioffa + 3 + 3 * lda ] + " "
    + A[ ioffa + 4 + 3 * lda ]  + "\n";
  LaPack2.dgelq2(n,m,B,ldb,tau,work,info,ioffb,iofftau,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dgelq2(m<n): info,tau = " + info.value + " "
    + tau[ iofftau + 0 ] + " "
    + tau[ iofftau + 1 ] + " "
    + tau[ iofftau + 2 ] + " "
    + tau[ iofftau + 3 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "B = "
    + B[ ioffb + 0 + 0 * ldb ] + " "
    + B[ ioffb + 1 + 0 * ldb ] + " "
    + B[ ioffb + 2 + 0 * ldb ] + " "
    + B[ ioffb + 3 + 0 * ldb ] + " "
    + B[ ioffb + 0 + 1 * ldb ] + " "
    + B[ ioffb + 1 + 1 * ldb ] + " "
    + B[ ioffb + 2 + 1 * ldb ] + " "
    + B[ ioffb + 3 + 1 * ldb ] + " "
    + B[ ioffb + 0 + 2 * ldb ] + " "
    + B[ ioffb + 1 + 2 * ldb ] + " "
    + B[ ioffb + 2 + 2 * ldb ] + " "
    + B[ ioffb + 3 + 2 * ldb ] + " "
    + B[ ioffb + 0 + 3 * ldb ] + " "
    + B[ ioffb + 1 + 3 * ldb ] + " "
    + B[ ioffb + 2 + 3 * ldb ] + " "
    + B[ ioffb + 3 + 3 * ldb ] + " "
    + B[ ioffb + 0 + 4 * ldb ] + " "
    + B[ ioffb + 1 + 4 * ldb ] + " "
    + B[ ioffb + 2 + 4 * ldb ] + " "
    + B[ ioffb + 3 + 4 * ldb ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgeql2() {
  document.getElementById("debug_textarea").value +=
    "testing dgeql2 *************" + "\n";
  var n = 4;
  var m = 5;
  var lda = 6;
  var ldb = 5;
  var ioffa = 1;
  var ioffb = 4;
  var iofftau = 2;
  var ioffwork = 3;
  var A = new Array( ioffa + lda * n );
  var B = new Array( ioffb + ldb * m );
  var tau = new Array( iofftau + Math.min( m, n ) );
  var work = new Array( ioffwork + Math.max( m, n ) );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + 2 * j );
      B[ ioffb + j + i * ldb ] = 1. / Number( 1 + i + 2 * j );
    }
  }
  var info = new IntReference();
  LaPack2.dgeql2(m,n,A,lda,tau,work,info,ioffa,iofftau,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dgeql2(m>n): info,tau = " + info.value + " "
    + tau[ iofftau + 0 ] + " "
    + tau[ iofftau + 1 ] + " "
    + tau[ iofftau + 2 ] + " "
    + tau[ iofftau + 3 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 3 + 0 * lda ] + " "
    + A[ ioffa + 4 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 3 + 1 * lda ] + " "
    + A[ ioffa + 4 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 3 + 2 * lda ] + " "
    + A[ ioffa + 4 + 2 * lda ] + " "
    + A[ ioffa + 0 + 3 * lda ] + " "
    + A[ ioffa + 1 + 3 * lda ] + " "
    + A[ ioffa + 2 + 3 * lda ] + " "
    + A[ ioffa + 3 + 3 * lda ] + " "
    + A[ ioffa + 4 + 3 * lda ]  + "\n";
  LaPack2.dgeql2(n,m,B,ldb,tau,work,info,ioffb,iofftau,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dgeql2(m<n): info,tau = " + info.value + " "
    + tau[ iofftau + 0 ] + " "
    + tau[ iofftau + 1 ] + " "
    + tau[ iofftau + 2 ] + " "
    + tau[ iofftau + 3 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "B = "
    + B[ ioffb + 0 + 0 * ldb ] + " "
    + B[ ioffb + 1 + 0 * ldb ] + " "
    + B[ ioffb + 2 + 0 * ldb ] + " "
    + B[ ioffb + 3 + 0 * ldb ] + " "
    + B[ ioffb + 0 + 1 * ldb ] + " "
    + B[ ioffb + 1 + 1 * ldb ] + " "
    + B[ ioffb + 2 + 1 * ldb ] + " "
    + B[ ioffb + 3 + 1 * ldb ] + " "
    + B[ ioffb + 0 + 2 * ldb ] + " "
    + B[ ioffb + 1 + 2 * ldb ] + " "
    + B[ ioffb + 2 + 2 * ldb ] + " "
    + B[ ioffb + 3 + 2 * ldb ] + " "
    + B[ ioffb + 0 + 3 * ldb ] + " "
    + B[ ioffb + 1 + 3 * ldb ] + " "
    + B[ ioffb + 2 + 3 * ldb ] + " "
    + B[ ioffb + 3 + 3 * ldb ] + " "
    + B[ ioffb + 0 + 4 * ldb ] + " "
    + B[ ioffb + 1 + 4 * ldb ] + " "
    + B[ ioffb + 2 + 4 * ldb ] + " "
    + B[ ioffb + 3 + 4 * ldb ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgeqr2() {
  document.getElementById("debug_textarea").value +=
    "testing dgeqr2 *************" + "\n";
  var n = 4;
  var m = 5;
  var lda = 6;
  var ldb = 5;
  var ioffa = 1;
  var ioffb = 4;
  var iofftau = 2;
  var ioffwork = 3;
  var A = new Array( ioffa + lda * n );
  var B = new Array( ioffb + ldb * m );
  var tau = new Array( iofftau + Math.min( m, n ) );
  var work = new Array( ioffwork + Math.max( m, n ) );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + 2 * j );
      B[ ioffb + j + i * ldb ] = 1. / Number( 1 + i + 2 * j );
    }
  }
  var info = new IntReference();
  LaPack2.dgeqr2(m,n,A,lda,tau,work,info,ioffa,iofftau,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dgeqr2(m>n): info,tau = " + info.value + " "
    + tau[ iofftau + 0 ] + " "
    + tau[ iofftau + 1 ] + " "
    + tau[ iofftau + 2 ] + " "
    + tau[ iofftau + 3 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 3 + 0 * lda ] + " "
    + A[ ioffa + 4 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 3 + 1 * lda ] + " "
    + A[ ioffa + 4 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 3 + 2 * lda ] + " "
    + A[ ioffa + 4 + 2 * lda ] + " "
    + A[ ioffa + 0 + 3 * lda ] + " "
    + A[ ioffa + 1 + 3 * lda ] + " "
    + A[ ioffa + 2 + 3 * lda ] + " "
    + A[ ioffa + 3 + 3 * lda ] + " "
    + A[ ioffa + 4 + 3 * lda ]  + "\n";
  LaPack2.dgeqr2(n,m,B,ldb,tau,work,info,ioffb,iofftau,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dgeqr2(m<n): info,tau = " + info.value + " "
    + tau[ iofftau + 0 ] + " "
    + tau[ iofftau + 1 ] + " "
    + tau[ iofftau + 2 ] + " "
    + tau[ iofftau + 3 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "B = "
    + B[ ioffb + 0 + 0 * ldb ] + " "
    + B[ ioffb + 1 + 0 * ldb ] + " "
    + B[ ioffb + 2 + 0 * ldb ] + " "
    + B[ ioffb + 3 + 0 * ldb ] + " "
    + B[ ioffb + 0 + 1 * ldb ] + " "
    + B[ ioffb + 1 + 1 * ldb ] + " "
    + B[ ioffb + 2 + 1 * ldb ] + " "
    + B[ ioffb + 3 + 1 * ldb ] + " "
    + B[ ioffb + 0 + 2 * ldb ] + " "
    + B[ ioffb + 1 + 2 * ldb ] + " "
    + B[ ioffb + 2 + 2 * ldb ] + " "
    + B[ ioffb + 3 + 2 * ldb ] + " "
    + B[ ioffb + 0 + 3 * ldb ] + " "
    + B[ ioffb + 1 + 3 * ldb ] + " "
    + B[ ioffb + 2 + 3 * ldb ] + " "
    + B[ ioffb + 3 + 3 * ldb ] + " "
    + B[ ioffb + 0 + 4 * ldb ] + " "
    + B[ ioffb + 1 + 4 * ldb ] + " "
    + B[ ioffb + 2 + 4 * ldb ] + " "
    + B[ ioffb + 3 + 4 * ldb ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgeqr2p() {
  document.getElementById("debug_textarea").value +=
    "testing dgeqr2p *************" + "\n";
  var n = 4;
  var m = 5;
  var lda = 6;
  var ioffa = 1;
  var iofftau = 2;
  var ioffwork = 3;
  var A = new Array( ioffa + lda * n );
  var tau = new Array( iofftau + Math.min( m, n ) );
  var work = new Array( ioffwork + n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + 2 * j );
    }
  }
  var info = new IntReference();
  LaPack2.dgeqr2p(m,n,A,lda,tau,work,info,ioffa,iofftau,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dgeqr2p: info,tau = " + info.value + " "
    + tau[ iofftau + 0 ] + " "
    + tau[ iofftau + 1 ] + " "
    + tau[ iofftau + 2 ] + " "
    + tau[ iofftau + 3 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 3 + 0 * lda ] + " "
    + A[ ioffa + 4 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 3 + 1 * lda ] + " "
    + A[ ioffa + 4 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 3 + 2 * lda ] + " "
    + A[ ioffa + 4 + 2 * lda ] + " "
    + A[ ioffa + 0 + 3 * lda ] + " "
    + A[ ioffa + 1 + 3 * lda ] + " "
    + A[ ioffa + 2 + 3 * lda ] + " "
    + A[ ioffa + 3 + 3 * lda ] + " "
    + A[ ioffa + 4 + 3 * lda ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgerfs() {
  document.getElementById("debug_textarea").value +=
    "testing dgerfs *************" + "\n";
  var n = 3;
  var lda = 4;
  var ldaf = 5;
  var ldb = 6;
  var ldx = 7;
  var nrhs = 1;
  var ioffa = 1;
  var ioffaf = 2;
  var ioffb = 3;
  var ioffx = 4;
  var ioffipiv = 5;
  var ioffiwork = 6;
  var ioffwork = 7;
  var ioffferr = 8;
  var ioffberr = 9;
  var A = new Array( ioffa + lda * n );
  var AF = new Array( ioffaf + ldaf * n );
  var work = new Array( ioffwork + 3 * n );
  var iwork = new Array( ioffiwork + n );
  var ipiv = new Array( ioffipiv + n );
  var B = new Array( ioffb + ldb * nrhs );
  var X = new Array( ioffx + ldx * nrhs );
  var ferr = new Array( ioffferr + nrhs );
  var berr = new Array( ioffberr + nrhs );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + 2 * j );
    }
  }
  LaPack0.dlacpy('A',n,n,A,lda,AF,ldaf,ioffa,ioffaf);
  var info = new IntReference();
  LaPack2.dgetrf(n,n,AF,ldaf,ipiv,info,ioffaf,ioffipiv);
  for ( i=0; i < n; i ++ ) B[ ioffb + i + 0 * ldb ] = i;
  Blas1.dcopy(n,B,1,X,1,ioffb,ioffx);
  LaPack1.dgetrs('N',n,nrhs,AF,ldaf,ipiv,X,ldx,info,
    ioffaf,ioffipiv,ioffx);
  LaPack2.dgerfs('N',n,nrhs,A,lda,AF,ldaf,ipiv,B,ldb,X,ldx,ferr,berr,
    work,iwork,info,ioffa,ioffaf,ioffipiv,ioffb,ioffx,ioffferr,ioffberr,
    ioffwork,ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dgerfs('N',...): info,ferr,berr = " + info.value + " "
    + ferr[ ioffferr ] + " " + berr[ ioffberr ] + "\n";

  Blas1.dcopy(n,B,1,X,1,ioffb,ioffx);
  LaPack1.dgetrs('T',n,nrhs,AF,ldaf,ipiv,X,ldx,info,
    ioffaf,ioffipiv,ioffx);
  LaPack2.dgerfs('T',n,nrhs,A,lda,AF,ldaf,ipiv,B,ldb,X,ldx,ferr,berr,
    work,iwork,info,ioffa,ioffaf,ioffipiv,ioffb,ioffx,ioffferr,ioffberr,
    ioffwork,ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dgerfs('T',...): info,ferr,berr = " + info.value + " "
    + ferr[ ioffferr ] + " " + berr[ ioffberr ] + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgerq2() {
  document.getElementById("debug_textarea").value +=
    "testing dgerq2 *************" + "\n";
  var n = 4;
  var m = 5;
  var lda = 6;
  var ldb = 5;
  var ioffa = 1;
  var ioffb = 4;
  var iofftau = 2;
  var ioffwork = 3;
  var A = new Array( ioffa + lda * n );
  var B = new Array( ioffb + ldb * m );
  var tau = new Array( iofftau + Math.min( m, n ) );
  var work = new Array( ioffwork + Math.max( m, n ) );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + 2 * j );
      B[ ioffb + j + i * ldb ] = 1. / Number( 1 + i + 2 * j );
    }
  }
  var info = new IntReference();
  LaPack2.dgerq2(m,n,A,lda,tau,work,info,ioffa,iofftau,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dgerq2(m>n): info,tau = " + info.value + " "
    + tau[ iofftau + 0 ] + " "
    + tau[ iofftau + 1 ] + " "
    + tau[ iofftau + 2 ] + " "
    + tau[ iofftau + 3 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 3 + 0 * lda ] + " "
    + A[ ioffa + 4 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 3 + 1 * lda ] + " "
    + A[ ioffa + 4 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 3 + 2 * lda ] + " "
    + A[ ioffa + 4 + 2 * lda ] + " "
    + A[ ioffa + 0 + 3 * lda ] + " "
    + A[ ioffa + 1 + 3 * lda ] + " "
    + A[ ioffa + 2 + 3 * lda ] + " "
    + A[ ioffa + 3 + 3 * lda ] + " "
    + A[ ioffa + 4 + 3 * lda ]  + "\n";
  LaPack2.dgerq2(n,m,B,ldb,tau,work,info,ioffb,iofftau,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dgerq2(m<n): info,tau = " + info.value + " "
    + tau[ iofftau + 0 ] + " "
    + tau[ iofftau + 1 ] + " "
    + tau[ iofftau + 2 ] + " "
    + tau[ iofftau + 3 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "B = "
    + B[ ioffb + 0 + 0 * ldb ] + " "
    + B[ ioffb + 1 + 0 * ldb ] + " "
    + B[ ioffb + 2 + 0 * ldb ] + " "
    + B[ ioffb + 3 + 0 * ldb ] + " "
    + B[ ioffb + 0 + 1 * ldb ] + " "
    + B[ ioffb + 1 + 1 * ldb ] + " "
    + B[ ioffb + 2 + 1 * ldb ] + " "
    + B[ ioffb + 3 + 1 * ldb ] + " "
    + B[ ioffb + 0 + 2 * ldb ] + " "
    + B[ ioffb + 1 + 2 * ldb ] + " "
    + B[ ioffb + 2 + 2 * ldb ] + " "
    + B[ ioffb + 3 + 2 * ldb ] + " "
    + B[ ioffb + 0 + 3 * ldb ] + " "
    + B[ ioffb + 1 + 3 * ldb ] + " "
    + B[ ioffb + 2 + 3 * ldb ] + " "
    + B[ ioffb + 3 + 3 * ldb ] + " "
    + B[ ioffb + 0 + 4 * ldb ] + " "
    + B[ ioffb + 1 + 4 * ldb ] + " "
    + B[ ioffb + 2 + 4 * ldb ] + " "
    + B[ ioffb + 3 + 4 * ldb ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgetrf() {
  document.getElementById("debug_textarea").value +=
    "testing dgetrf *************" + "\n";
  var m = 4;
  var n = 3;
  var lda = 5;
  var ioffa = 1;
  var ioffipiv = 2;
  var A = new Array( ioffa + lda * n );
  var ipiv = new Array( ioffipiv + Math.min( m , n ) );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + 2 * i + Math.pow( i + 1, j );
    }
  }
  var info = new IntReference();
  LaPack2.dgetrf(m,n,A,lda,ipiv,info,ioffa,ioffipiv);
  document.getElementById("debug_textarea").value +=
    "dgetrf: info = " + info.value + "\n";
  document.getElementById("debug_textarea").value +=
    "ipiv = "
    + ipiv[ ioffipiv + 0 ] + " "
    + ipiv[ ioffipiv + 1 ] + " "
    + ipiv[ ioffipiv + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 3 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 3 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 3 + 2 * lda ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgtcon() {
  document.getElementById("debug_textarea").value +=
    "testing dgtcon *************" + "\n";
  var n = 4;
  var ioffd = 1;
  var ioffdl = 2;
  var ioffdu = 3;
  var ioffdu2 = 4;
  var ioffipiv = 5;
  var ioffwork = 6;
  var ioffiwork = 7;
  var d = new Array( ioffd + n );
  var dl = new Array( ioffdl + n - 1 );
  var du = new Array( ioffdu + n - 1 );
  var du2 = new Array( ioffdu2 + n - 2 );
  var work = new Array( ioffwork + 2 * n );
  var iwork = new Array( ioffiwork + n );
  for ( var i = 0; i < n; i ++ ) d[ ioffd + i ] = 2.;
  for ( i=0; i < n - 1; i ++ ) {
    dl[ ioffdl + i ] = -1.;
    du[ ioffdu + i ] = -1.;
  }
  var ipiv = new Array( ioffipiv + n );
  var info = new IntReference();
  LaPack0.dgttrf(n,dl,d,du,du2,ipiv,info,
    ioffdl,ioffd,ioffdu,ioffdu2,ioffipiv);
  var anorm = 4.;
  var rcond = new NumberReference();
  LaPack2.dgtcon('O',n,dl,d,du,du2,ipiv,anorm,rcond,work,iwork,info,
    ioffdl,ioffd,ioffdu,ioffdu2,ioffipiv,ioffwork,ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dgtcon('O',...): rcond = " + rcond.value + "\n";
  LaPack2.dgtcon('I',n,dl,d,du,du2,ipiv,anorm,rcond,work,iwork,info,
    ioffdl,ioffd,ioffdu,ioffdu2,ioffipiv,ioffwork,ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dgtcon('I',...): rcond = " + rcond.value + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgtrfs() {
  document.getElementById("debug_textarea").value +=
    "testing dgtrfs *************" + "\n";
  var n = 4;
  var ioffd = 1;
  var ioffdl = 2;
  var ioffdu = 3;
  var ioffdu2 = 4;
  var ioffipiv = 5;
  var ioffwork = 6;
  var ioffiwork = 7;
  var d = new Array( ioffd + n );
  var dl = new Array( ioffdl + n - 1 );
  var du = new Array( ioffdu + n - 1 );
  var du2 = new Array( ioffdu2 + n - 2 );
  var work = new Array( ioffwork + 2 * n );
  var iwork = new Array( ioffiwork + n );
  for ( var i = 0; i < n; i ++ ) d[ ioffd + i ] = 2.;
  for ( i=0; i < n - 1; i ++ ) {
    dl[ ioffdl + i ] = -1.;
    du[ ioffdu + i ] = -1.;
  }
  var ipiv = new Array( ioffipiv + n );

  var ioffdf = 8;
  var ioffdlf = 9;
  var ioffduf = 10;
  var df = new Array( ioffdf + n );
  var dlf = new Array( ioffdlf + n - 1 );
  var duf = new Array( ioffduf + n - 1 );
  Blas1.dcopy(n,d,1,df,1,ioffd,ioffdf);
  Blas1.dcopy(n-1,dl,1,dlf,1,ioffdl,ioffdlf);
  Blas1.dcopy(n-1,du,1,duf,1,ioffdu,ioffduf);
  var info = new IntReference();
  LaPack0.dgttrf(n,dlf,df,duf,du2,ipiv,info,
    ioffdlf,ioffdf,ioffduf,ioffdu2,ioffipiv);

  var nrhs = 1;
  var ioffb = 11;
  var ioffx = 12;
  var ldb = 5;
  var ldx = 6;
  var B = new Array( ioffb + ldb * nrhs );
  var X = new Array( ioffx + ldx * nrhs );
  B[ ioffb + 0 ] = 0.;
  B[ ioffb + 1 ] = 0.;
  B[ ioffb + 2 ] = 0.;
  B[ ioffb + 3 ] = 5.;
  Blas1.dcopy(n,B,1,X,1,ioffb,ioffx);
  LaPack1.dgttrs('N',n,nrhs,dlf,df,duf,du2,ipiv,X,ldx,info,
    ioffdlf,ioffdf,ioffduf,ioffdu2,ioffipiv,ioffx);
  var ioffferr = 13;
  var ioffberr = 14;
  var ferr = new Array( ioffferr + nrhs );
  var berr = new Array( ioffberr + nrhs );
  LaPack2.dgtrfs('N',n,nrhs,dl,d,du,dlf,df,duf,du2,ipiv,B,ldb,X,ldx,
    ferr,berr,work,iwork,info,
    ioffdl,ioffd,ioffdu,ioffdlf,ioffdf,ioffduf,ioffdu2,ioffipiv,
      ioffb,ioffx,ioffferr,ioffberr,ioffwork,ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dgtrfs('N',...): berr,ferr = " + berr[ ioffberr ] + " "
    + ferr[ ioffferr ] + "\n";
  document.getElementById("debug_textarea").value +=
    "X = "
    + X[ ioffx + 0 ] + " "
    + X[ ioffx + 1 ] + " "
    + X[ ioffx + 2 ] + " "
    + X[ ioffx + 3 ]  + "\n";

  Blas1.dcopy(n,B,1,X,1,ioffb,ioffx);
  LaPack1.dgttrs('T',n,nrhs,dlf,df,duf,du2,ipiv,X,ldx,info,
    ioffdlf,ioffdf,ioffduf,ioffdu2,ioffipiv,ioffx);
  LaPack2.dgtrfs('T',n,nrhs,dl,d,du,dlf,df,duf,du2,ipiv,B,ldb,X,ldx,
    ferr,berr,work,iwork,info,
    ioffdl,ioffd,ioffdu,ioffdlf,ioffdf,ioffduf,ioffdu2,ioffipiv,
      ioffb,ioffx,ioffferr,ioffberr,ioffwork,ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dgtrfs('T',...): berr,ferr = " + berr[ ioffberr ] + " "
    + ferr[ ioffferr ] + "\n";
  document.getElementById("debug_textarea").value +=
    "X = "
    + X[ ioffx + 0 ] + " "
    + X[ ioffx + 1 ] + " "
    + X[ ioffx + 2 ] + " "
    + X[ ioffx + 3 ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlaed4() {
  document.getElementById("debug_textarea").value +=
    "testing dlaed4 *************" + "\n";
  var ioffd = 1;
  var ioffdelta = 2;
  var ioffz = 3;
  var d = new Array( ioffd + 3 );
  var delta = new Array( ioffdelta + 3 );
  var z = new Array( ioffz + 3 );
  for ( var i = 0; i < 3; i ++ ) {
    d[ ioffd + i ] = i + 1;
  }
  z[ ioffz ] = 1.;
  var rho = 2.;
  var info = new IntReference();
  var dlam = new NumberReference();
  var zn = 0.;
  LaPack2.dlaed4(1,1,d,z,delta,rho,dlam,info,ioffd,ioffz,ioffdelta);
  document.getElementById("debug_textarea").value +=
    "dlaed4(1,...): info,dlam = " + info.value + " " + dlam.value + "\n";
  document.getElementById("debug_textarea").value +=
    "delta = " + delta[ ioffdelta ] + "\n";

  for ( i=0; i < 2; i ++ ) z[ ioffz + i ] = 3 - i;
  zn = 1. / Blas1.dnrm2( 2, z, 1, ioffz );
  Blas1.dscal( 2, zn, z, 1, ioffz );
  for ( i = 1; i <= 2; i ++ ) {
    LaPack2.dlaed4(2,i,d,z,delta,rho,dlam,info,ioffd,ioffz,ioffdelta);
    document.getElementById("debug_textarea").value +=
      "dlaed4(2," + i + ",...): info,dlam = " + info.value + " "
      + dlam.value + "\n";
    document.getElementById("debug_textarea").value +=
      "delta = "
      + delta[ ioffdelta + 0 ] + " "
      + delta[ ioffdelta + 1 ] + "\n";
  }

  for ( i=0; i < 3; i ++ ) z[ ioffz + i ] = 3 - i;
  zn = 1. / Blas1.dnrm2( 3, z, 1, ioffz );
  Blas1.dscal( 3, zn, z, 1, ioffz );
  for ( i = 1; i <= 3; i ++ ) {
    LaPack2.dlaed4(3,i,d,z,delta,rho,dlam,info,ioffd,ioffz,ioffdelta);
    document.getElementById("debug_textarea").value +=
      "dlaed4(3," + i + ",...): info,dlam = " + info.value + " "
      + dlam.value + "\n";
    document.getElementById("debug_textarea").value +=
      "delta = "
      + delta[ ioffdelta + 0 ] + " "
      + delta[ ioffdelta + 1 ] + " "
      + delta[ ioffdelta + 2 ] + "\n";
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlaein() {
  document.getElementById("debug_textarea").value +=
    "testing dlaein *************" + "\n";
  var n = 3;
  var ldh = 4;
  var ldb = 5;
  var eps3 = 1.e-15;
  var smlnum = 1.e-300;
  var bignum = 1.e+300;
  var ioffh = 1;
  var ioffb = 2;
  var ioffvi = 3;
  var ioffvr = 4;
  var ioffwork = 5;
  var H = new Array( ioffh + ldh * n );
  var B = new Array( ioffb + ldb * n );
  var vi = new Array( ioffvi + n );
  var vr = new Array( ioffvr + n );
  var work = new Array( ioffwork + n );
//    companion matrix for (x-1)*(x-2)*(x-3)
  H[ ioffh + 0 + 0 * ldh ] = 6.;
  H[ ioffh + 1 + 0 * ldh ] = 1.;
  H[ ioffh + 2 + 0 * ldh ] = 0.;
  H[ ioffh + 0 + 1 * ldh ] = -11.;
  H[ ioffh + 1 + 1 * ldh ] = 0.;
  H[ ioffh + 2 + 1 * ldh ] = 1.;
  H[ ioffh + 0 + 2 * ldh ] = 6.;
  H[ ioffh + 1 + 2 * ldh ] = 0.;
  H[ ioffh + 2 + 2 * ldh ] = 0.;
  var wr = 1.;
  var wi = 0.;
  vr[ ioffvr + 0 ] = 1.;
  vr[ ioffvr + 1 ] = 0.;
  vr[ ioffvr + 2 ] = 0.;
  var info = new IntReference();
  LaPack2.dlaein(true,false,n,H,ldh,wr,wi,vr,vi,B,ldb,work,
    eps3,smlnum,bignum,info, ioffh,ioffvr,ioffvi,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaein(T,F,...): info,w = " + info.value + " " + wr + " "
    + wi + "\n";
  document.getElementById("debug_textarea").value +=
    "vr = "
    + vr[ ioffvr + 0 ] + " "
    + vr[ ioffvr + 1 ] + " "
    + vr[ ioffvr + 2 ]  + "\n";
  vr[ ioffvr + 0 ] = 1.;
  vr[ ioffvr + 1 ] = 0.;
  vr[ ioffvr + 2 ] = 0.;
  LaPack2.dlaein(false,false,n,H,ldh,wr,wi,vr,vi,B,ldb,work,
    eps3,smlnum,bignum,info, ioffh,ioffvr,ioffvi,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaein(F,F,...): info,w = " + info.value + " " + wr + " "
    + wi + "\n";
  document.getElementById("debug_textarea").value +=
    "vr = "
    + vr[ ioffvr + 0 ] + " "
    + vr[ ioffvr + 1 ] + " "
    + vr[ ioffvr + 2 ]  + "\n";
  LaPack2.dlaein(true,true,n,H,ldh,wr,wi,vr,vi,B,ldb,work,
    eps3,smlnum,bignum,info, ioffh,ioffvr,ioffvi,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaein(T,T,...): info,w = " + info.value + " " + wr + " "
    + wi + "\n";
  document.getElementById("debug_textarea").value +=
    "vr = "
    + vr[ ioffvr + 0 ] + " "
    + vr[ ioffvr + 1 ] + " "
    + vr[ ioffvr + 2 ]  + "\n";
  LaPack2.dlaein(false,true,n,H,ldh,wr,wi,vr,vi,B,ldb,work,
    eps3,smlnum,bignum,info, ioffh,ioffvr,ioffvi,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaein(F,T,...): info,w = " + info.value + " " + wr + " "
    + wi + "\n";
  document.getElementById("debug_textarea").value +=
    "vr = "
    + vr[ ioffvr + 0 ] + " "
    + vr[ ioffvr + 1 ] + " "
    + vr[ ioffvr + 2 ]  + "\n";
  
//    companion matrix for (x-1)^3
  H[ ioffh + 0 + 0 * ldh ] = 3.;
  H[ ioffh + 1 + 0 * ldh ] = 1.;
  H[ ioffh + 2 + 0 * ldh ] = 0.;
  H[ ioffh + 0 + 1 * ldh ] = -3.;
  H[ ioffh + 1 + 1 * ldh ] = 0.;
  H[ ioffh + 2 + 1 * ldh ] = 1.;
  H[ ioffh + 0 + 2 * ldh ] = 1.;
  H[ ioffh + 1 + 2 * ldh ] = 0.;
  H[ ioffh + 2 + 2 * ldh ] = 0.;
  wr = 1.;
  wi = 0.;
  vr[ ioffvr + 0 ] = 1.;
  vr[ ioffvr + 1 ] = 0.;
  vr[ ioffvr + 2 ] = 0.;
  LaPack2.dlaein(true,false,n,H,ldh,wr,wi,vr,vi,B,ldb,work,
    eps3,smlnum,bignum,info, ioffh,ioffvr,ioffvi,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaein(T,F,...): info,w = " + info.value + " " + wr + " "
    + wi + "\n";
  document.getElementById("debug_textarea").value +=
    "vr = "
    + vr[ ioffvr + 0 ] + " "
    + vr[ ioffvr + 1 ] + " "
    + vr[ ioffvr + 2 ]  + "\n";
  vr[ ioffvr + 0 ] = 1.;
  vr[ ioffvr + 1 ] = 0.;
  vr[ ioffvr + 2 ] = 0.;
  LaPack2.dlaein(false,false,n,H,ldh,wr,wi,vr,vi,B,ldb,work,
    eps3,smlnum,bignum,info, ioffh,ioffvr,ioffvi,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaein(F,F,...): info,w = " + info.value + " " + wr + " "
    + wi + "\n";
  document.getElementById("debug_textarea").value +=
    "vr = "
    + vr[ ioffvr + 0 ] + " "
    + vr[ ioffvr + 1 ] + " "
    + vr[ ioffvr + 2 ]  + "\n";
  LaPack2.dlaein(true,true,n,H,ldh,wr,wi,vr,vi,B,ldb,work,
    eps3,smlnum,bignum,info, ioffh,ioffvr,ioffvi,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaein(T,T,...): info,w = " + info.value + " " + wr + " "
    + wi + "\n";
  document.getElementById("debug_textarea").value +=
    "vr = "
    + vr[ ioffvr + 0 ] + " "
    + vr[ ioffvr + 1 ] + " "
    + vr[ ioffvr + 2 ]  + "\n";
  LaPack2.dlaein(false,true,n,H,ldh,wr,wi,vr,vi,B,ldb,work,
    eps3,smlnum,bignum,info, ioffh,ioffvr,ioffvi,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaein(F,T,...): info,w = " + info.value + " " + wr + " "
    + wi + "\n";
  document.getElementById("debug_textarea").value +=
    "vr = "
    + vr[ ioffvr + 0 ] + " "
    + vr[ ioffvr + 1 ] + " "
    + vr[ ioffvr + 2 ]  + "\n";
  
//    companion matrix for (x-1)*(x-[1+i])*(x-[1-i]) = x^3-3x^2+4x-2
  H[ ioffh + 0 + 0 * ldh ] = 3.;
  H[ ioffh + 1 + 0 * ldh ] = 1.;
  H[ ioffh + 2 + 0 * ldh ] = 0.;
  H[ ioffh + 0 + 1 * ldh ] = -4.;
  H[ ioffh + 1 + 1 * ldh ] = 0.;
  H[ ioffh + 2 + 1 * ldh ] = 1.;
  H[ ioffh + 0 + 2 * ldh ] = 2.;
  H[ ioffh + 1 + 2 * ldh ] = 0.;
  H[ ioffh + 2 + 2 * ldh ] = 0.;
  wr = 1.;
  wi = 0.;
  vr[ ioffvr + 0 ] = 1.;
  vr[ ioffvr + 1 ] = 0.;
  vr[ ioffvr + 2 ] = 0.;
  LaPack2.dlaein(true,false,n,H,ldh,wr,wi,vr,vi,B,ldb,work,
    eps3,smlnum,bignum,info, ioffh,ioffvr,ioffvi,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaein(T,F,...): info,w = " + info.value + " " + wr + " "
    + wi + "\n";
  document.getElementById("debug_textarea").value +=
    "vr = "
    + vr[ ioffvr + 0 ] + " "
    + vr[ ioffvr + 1 ] + " "
    + vr[ ioffvr + 2 ]  + "\n";
  vr[ ioffvr + 0 ] = 1.;
  vr[ ioffvr + 1 ] = 0.;
  vr[ ioffvr + 2 ] = 0.;
  LaPack2.dlaein(false,false,n,H,ldh,wr,wi,vr,vi,B,ldb,work,
    eps3,smlnum,bignum,info, ioffh,ioffvr,ioffvi,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaein(F,F,...): info,w = " + info.value + " " + wr + " "
    + wi + "\n";
  document.getElementById("debug_textarea").value +=
    "vr = "
    + vr[ ioffvr + 0 ] + " "
    + vr[ ioffvr + 1 ] + " "
    + vr[ ioffvr + 2 ]  + "\n";
  LaPack2.dlaein(true,true,n,H,ldh,wr,wi,vr,vi,B,ldb,work,
    eps3,smlnum,bignum,info, ioffh,ioffvr,ioffvi,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaein(T,T,...): info,w = " + info.value + " " + wr + " "
    + wi + "\n";
  document.getElementById("debug_textarea").value +=
    "vr = "
    + vr[ ioffvr + 0 ] + " "
    + vr[ ioffvr + 1 ] + " "
    + vr[ ioffvr + 2 ] + "\n";
  LaPack2.dlaein(false,true,n,H,ldh,wr,wi,vr,vi,B,ldb,work,
    eps3,smlnum,bignum,info, ioffh,ioffvr,ioffvi,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaein(F,T,...): info,w = " + info.value + " " + wr + " "
    + wi + "\n";
  document.getElementById("debug_textarea").value +=
    "vr = "
    + vr[ ioffvr + 0 ] + " "
    + vr[ ioffvr + 1 ] + " "
    + vr[ ioffvr + 2 ]  + "\n";

  wr = 1.;
  wi = 1.;
  vr[ ioffvr + 0 ] = 1.;
  vr[ ioffvr + 1 ] = 0.;
  vr[ ioffvr + 2 ] = 0.;
  vi[ ioffvi + 0 ] = 0.;
  vi[ ioffvi + 1 ] = 0.;
  vi[ ioffvi + 2 ] = 1.;
  LaPack2.dlaein(true,false,n,H,ldh,wr,wi,vr,vi,B,ldb,work,
    eps3,smlnum,bignum,info, ioffh,ioffvr,ioffvi,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaein(T,F,...): info,w = " + info.value + " " + wr + " "
    + wi + "\n";
  var c = vr[ ioffvr ];
  var d = vi[ ioffvi ];
  var p = new NumberReference();
  var q = new NumberReference();
  for ( var i = 0; i < n; i ++ ) {
    var a = vr[ ioffvr + i ];
    var b = vi[ ioffvi + i ];
    LaPack0.dladiv(a,b,c,d,p,q);
    vr[ ioffvr + i ] = p.value;
    vi[ ioffvi + i ] = q.value;
  }
  document.getElementById("debug_textarea").value +=
    "vr = "
    + vr[ ioffvr + 0 ] + " "
    + vr[ ioffvr + 1 ] + " "
    + vr[ ioffvr + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "vi = "
    + vi[ ioffvi + 0 ] + " "
    + vi[ ioffvi + 1 ] + " "
    + vi[ ioffvi + 2 ]  + "\n";
  vr[ ioffvr + 0 ] = 1.;
  vr[ ioffvr + 1 ] = 0.;
  vr[ ioffvr + 2 ] = 0.;
  vi[ ioffvi + 0 ] = 0.;
  vi[ ioffvi + 1 ] = 0.;
  vi[ ioffvi + 2 ] = 1.;
  LaPack2.dlaein(false,false,n,H,ldh,wr,wi,vr,vi,B,ldb,work,
    eps3,smlnum,bignum,info, ioffh,ioffvr,ioffvi,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaein(F,F,...): info,w = " + info.value + " " + wr + " "
    + wi + "\n";
  c = vr[ ioffvr ];
  d = vi[ ioffvi ];
  for ( i = 0; i < n; i ++ ) {
    a = vr[ ioffvr + i ];
    b = vi[ ioffvi + i ];
    LaPack0.dladiv(a,b,c,d,p,q);
    vr[ ioffvr + i ] = p.value;
    vi[ ioffvi + i ] = q.value;
  }
  document.getElementById("debug_textarea").value +=
    "vr = "
    + vr[ ioffvr + 0 ] + " "
    + vr[ ioffvr + 1 ] + " "
    + vr[ ioffvr + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "vi = "
    + vi[ ioffvi + 0 ] + " "
    + vi[ ioffvi + 1 ] + " "
    + vi[ ioffvi + 2 ]  + "\n";
  LaPack2.dlaein(true,true,n,H,ldh,wr,wi,vr,vi,B,ldb,work,
    eps3,smlnum,bignum,info, ioffh,ioffvr,ioffvi,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaein(T,T,...): info,w = " + info.value + " " + wr + " "
    + wi + "\n";
  c = vr[ ioffvr ];
  d = vi[ ioffvi ];
  for ( i = 0; i < n; i ++ ) {
    a = vr[ ioffvr + i ];
    b = vi[ ioffvi + i ];
    LaPack0.dladiv(a,b,c,d,p,q);
    vr[ ioffvr + i ] = p.value;
    vi[ ioffvi + i ] = q.value;
  }
  document.getElementById("debug_textarea").value +=
    "vr = "
    + vr[ ioffvr + 0 ] + " "
    + vr[ ioffvr + 1 ] + " "
    + vr[ ioffvr + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "vi = "
    + vi[ ioffvi + 0 ] + " "
    + vi[ ioffvi + 1 ] + " "
    + vi[ ioffvi + 2 ]  + "\n";
  LaPack2.dlaein(false,true,n,H,ldh,wr,wi,vr,vi,B,ldb,work,
    eps3,smlnum,bignum,info, ioffh,ioffvr,ioffvi,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaein(F,T,...): info,w = " + info.value + " " + wr + " "
    + wi + "\n";
  c = vr[ ioffvr ];
  d = vi[ ioffvi ];
  for ( i = 0; i < n; i ++ ) {
    a = vr[ ioffvr + i ];
    b = vi[ ioffvi + i ];
    LaPack0.dladiv(a,b,c,d,p,q);
    vr[ ioffvr + i ] = p.value;
    vi[ ioffvi + i ] = q.value;
  }
  document.getElementById("debug_textarea").value +=
    "vr = "
    + vr[ ioffvr + 0 ] + " "
    + vr[ ioffvr + 1 ] + " "
    + vr[ ioffvr + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "vi = "
    + vi[ ioffvi + 0 ] + " "
    + vi[ ioffvi + 1 ] + " "
    + vi[ ioffvi + 2 ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlapll() {
  document.getElementById("debug_textarea").value +=
    "testing dlapll *************" + "\n";
  var n = 3;
  var incx = 2;
  var incy = 3;
  var ioffx = 1;
  var ioffy = 2;
  var x = new Array( ioffx + n * incx );
  var y = new Array( ioffy + n * incy );
  for ( var i = 0; i < n * incx; i ++ ) {
    x[ ioffx + i ] = i + 1;
  }
  for ( i = 0; i < n * incy; i ++ ) {
    y[ ioffy + i ] = 2 * i  - 3;
  }
  var ssmin = new NumberReference();
  LaPack2.dlapll(n,x,incx,y,incy,ssmin,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
    "dlapll: ssmin = " + ssmin.value + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ] + " "
    + x[ ioffx + 3 ] + " "
    + x[ ioffx + 4 ] + " "
    + x[ ioffx + 5 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "y = "
    + y[ ioffy + 0 ] + " "
    + y[ ioffy + 1 ] + " "
    + y[ ioffy + 2 ] + " "
    + y[ ioffy + 3 ] + " "
    + y[ ioffy + 4 ] + " "
    + y[ ioffy + 5 ] + " "
    + y[ ioffy + 6 ] + " "
    + y[ ioffy + 7 ] + " "
    + y[ ioffy + 8 ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlarfx() {
  document.getElementById("debug_textarea").value +=
    "testing dlarfx *************" + "\n";
  var ldc = 5;
  var m = 4;
  var n = 3;
  var ioffC = 1;
  var ioffv = 2;
  var ioffwork = 3;
  var C = new Array( ioffC + ldc * n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] = i + 2 * j;
    }
  }
  var work = new Array( ioffwork + Math.max(m,n) );
  var v = new Array( ioffv + m );
  for ( i = 0; i < m; i ++ ) v[ ioffv + i ] = 1 - 2 * i;
  var tau = 3.;
  LaPack2.dlarfx('L',m,n,v,tau,C,ldc,work,ioffv,ioffC,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlarfx('L',...): C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 3 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 3 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + " "
    + C[ ioffC + 3 + 2 * ldc ]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] = i + 2 * j;
    }
  }
  v = new Array( ioffv + n );
  for ( j = 0; j < n; j ++ ) v[ ioffv + j ] = 1 - 2 * j;
  LaPack2.dlarfx('R',m,n,v,tau,C,ldc,work,ioffv,ioffC,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlarfx('R',...): C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 3 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 3 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + " "
    + C[ ioffC + 3 + 2 * ldc ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlatrd() {
  document.getElementById("debug_textarea").value +=
    "testing dlatrd *************" + "\n";
  var n = 4;
  var nb = 2;
  var lda = 5;
  var ldw = 6;
  var ioffa = 1;
  var ioffe = 2;
  var iofftau = 3;
  var ioffw = 4;
  var A = new Array( ioffa + lda * n );
  var e = new Array( ioffe + n - 1 );
  var tau = new Array( iofftau + n - 1 );
  var W = new Array( ioffw + ldw * nb );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + j );
    }
  }
  LaPack2.dlatrd('U',n,nb,A,lda,e,tau,W,ldw,ioffa,ioffe,iofftau,ioffw);
  document.getElementById("debug_textarea").value +=
    "dlatrd('U',...): e = "
    + e[ ioffe + 1 ] + " "
    + e[ ioffe + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "tau = "
    + tau[ iofftau + 1 ] + " "
    + tau[ iofftau + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 3 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 3 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 3 + 2 * lda ] + " "
    + A[ ioffa + 0 + 3 * lda ] + " "
    + A[ ioffa + 1 + 3 * lda ] + " "
    + A[ ioffa + 2 + 3 * lda ] + " "
    + A[ ioffa + 3 + 3 * lda ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "W = "
    + W[ ioffw + 0 + 0 * ldw ] + " "
    + W[ ioffw + 1 + 0 * ldw ] + " "
    + W[ ioffw + 2 + 0 * ldw ] + " "
    + W[ ioffw + 3 + 0 * ldw ] + " "
    + W[ ioffw + 0 + 1 * ldw ] + " "
    + W[ ioffw + 1 + 1 * ldw ] + " "
    + W[ ioffw + 2 + 1 * ldw ] + " "
    + W[ ioffw + 3 + 1 * ldw ]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + j );
    }
  }
  LaPack2.dlatrd('L',n,nb,A,lda,e,tau,W,ldw,ioffa,ioffe,iofftau,ioffw);
  document.getElementById("debug_textarea").value +=
    "dlatrd('L',...): e = "
    + e[ ioffe + 0 ] + " "
    + e[ ioffe + 1 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "tau = "
    + tau[ iofftau + 0 ] + " "
    + tau[ iofftau + 1 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 3 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 3 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 3 + 2 * lda ] + " "
    + A[ ioffa + 0 + 3 * lda ] + " "
    + A[ ioffa + 1 + 3 * lda ] + " "
    + A[ ioffa + 2 + 3 * lda ] + " "
    + A[ ioffa + 3 + 3 * lda ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "W = "
    + W[ ioffw + 0 + 0 * ldw ] + " "
    + W[ ioffw + 1 + 0 * ldw ] + " "
    + W[ ioffw + 2 + 0 * ldw ] + " "
    + W[ ioffw + 3 + 0 * ldw ] + " "
    + W[ ioffw + 0 + 1 * ldw ] + " "
    + W[ ioffw + 1 + 1 * ldw ] + " "
    + W[ ioffw + 2 + 1 * ldw ] + " "
    + W[ ioffw + 3 + 1 * ldw ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dorg2l() {
  document.getElementById("debug_textarea").value +=
    "testing dorg2l *************" + "\n";
  var lda = 5;
  var m = 4;
  var n = 3;
  var k = 2;
  var ioffa = 1;
  var iofftau = 2;
  var ioffwork = 3;
  var A = new Array( ioffa + lda * n );
  var tau = new Array( iofftau + k );
  var work = new Array( ioffwork + n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  for ( i = 0; i < k; i ++ ) {
    tau[ iofftau + i ] = 3 * i + 2;
  }
  var info = new IntReference();
  LaPack2.dorg2l(m,n,k,A,lda,tau,work,info,ioffa,iofftau,ioffwork);
  document.getElementById("debug_textarea").value +=
    "A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 3 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 3 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 3 + 2 * lda ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dorg2r() {
  document.getElementById("debug_textarea").value +=
    "testing dorg2r *************" + "\n";
  var lda = 5;
  var m = 4;
  var n = 3;
  var k = 2;
  var ioffa = 1;
  var iofftau = 2;
  var ioffwork = 3;
  var A = new Array( ioffa + lda * n );
  var tau = new Array( iofftau + k );
  var work = new Array( ioffwork + n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  for ( i = 0; i < k; i ++ ) {
    tau[ iofftau + i ] = 3 * i + 2;
  }
  var info = new IntReference();
  LaPack2.dorg2r(m,n,k,A,lda,tau,work,info,ioffa,iofftau,ioffwork);
  document.getElementById("debug_textarea").value +=
    "A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 3 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 3 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 3 + 2 * lda ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dorgl2() {
  document.getElementById("debug_textarea").value +=
    "testing dorgl2 *************" + "\n";
  var lda = 4;
  var m = 3;
  var n = 4;
  var k = 2;
  var ioffa = 1;
  var iofftau = 2;
  var ioffwork = 3;
  var A = new Array( ioffa + lda * n );
  var tau = new Array( iofftau + k );
  var work = new Array( ioffwork + n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  for ( i = 0; i < k; i ++ ) {
    tau[ iofftau + i ] = 3 * i + 2;
  }
  var info = new IntReference();
  LaPack2.dorgl2(m,n,k,A,lda,tau,work,info,ioffa,iofftau,ioffwork);
  document.getElementById("debug_textarea").value +=
    "A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 0 + 3 * lda ] + " "
    + A[ ioffa + 1 + 3 * lda ] + " "
    + A[ ioffa + 2 + 3 * lda ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dorgr2() {
  document.getElementById("debug_textarea").value +=
    "testing dorgr2 *************" + "\n";
  var lda = 4;
  var m = 3;
  var n = 4;
  var k = 2;
  var ioffa = 1;
  var iofftau = 2;
  var ioffwork = 3;
  var A = new Array( ioffa + lda * n );
  var tau = new Array( iofftau + k );
  var work = new Array( ioffwork + n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  for ( i = 0; i < k; i ++ ) {
    tau[ iofftau + i ] = 3 * i + 2;
  }
  var info = new IntReference();
  LaPack2.dorgr2(m,n,k,A,lda,tau,work,info,ioffa,iofftau,ioffwork);
  document.getElementById("debug_textarea").value +=
    "A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 0 + 3 * lda ] + " "
    + A[ ioffa + 1 + 3 * lda ] + " "
    + A[ ioffa + 2 + 3 * lda ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dorm2l() {
  document.getElementById("debug_textarea").value +=
    "testing dorm2l *************" + "\n";
  var lda = 5;
  var ldc = 6;
  var m = 4;
  var n = 3;
  var k = 2;
  var ioffa = 1;
  var ioffc = 2;
  var iofftau = 3;
  var ioffwork = 4;
  var A = new Array( ioffa + lda * k );
  var C = new Array( ioffc + ldc * n );
  var tau = new Array( iofftau + k );
  var work = new Array( ioffwork + Math.max( m, n ) );
  for ( var j = 0; j < k; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = 2 * i + j - 1;
    }
  }
  for ( i = 0; i < k; i ++ ) tau[ iofftau + i ] = 3 * i + 2;
  var info = new IntReference();
  LaPack2.dorm2l('L','N',m,n,k,A,lda,tau,C,ldc,work,info,
    ioffa,iofftau,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dorm2l('L','N',...): C = "
    + C[ ioffc + 0 + 0 * ldc ] + " "
    + C[ ioffc + 1 + 0 * ldc ] + " "
    + C[ ioffc + 2 + 0 * ldc ] + " "
    + C[ ioffc + 3 + 0 * ldc ] + " "
    + C[ ioffc + 0 + 1 * ldc ] + " "
    + C[ ioffc + 1 + 1 * ldc ] + " "
    + C[ ioffc + 2 + 1 * ldc ] + " "
    + C[ ioffc + 3 + 1 * ldc ] + " "
    + C[ ioffc + 0 + 2 * ldc ] + " "
    + C[ ioffc + 1 + 2 * ldc ] + " "
    + C[ ioffc + 2 + 2 * ldc ] + " "
    + C[ ioffc + 3 + 2 * ldc ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = 2 * i + j - 1;
    }
  }
  LaPack2.dorm2l('L','T',m,n,k,A,lda,tau,C,ldc,work,info,
    ioffa,iofftau,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dorm2l('L','T',...): C = "
    + C[ ioffc + 0 + 0 * ldc ] + " "
    + C[ ioffc + 1 + 0 * ldc ] + " "
    + C[ ioffc + 2 + 0 * ldc ] + " "
    + C[ ioffc + 3 + 0 * ldc ] + " "
    + C[ ioffc + 0 + 1 * ldc ] + " "
    + C[ ioffc + 1 + 1 * ldc ] + " "
    + C[ ioffc + 2 + 1 * ldc ] + " "
    + C[ ioffc + 3 + 1 * ldc ] + " "
    + C[ ioffc + 0 + 2 * ldc ] + " "
    + C[ ioffc + 1 + 2 * ldc ] + " "
    + C[ ioffc + 2 + 2 * ldc ] + " "
    + C[ ioffc + 3 + 2 * ldc ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = 2 * i + j - 1;
    }
  }
  LaPack2.dorm2l('R','N',m,n,k,A,lda,tau,C,ldc,work,info,
    ioffa,iofftau,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dorm2l('R','N',...): C = "
    + C[ ioffc + 0 + 0 * ldc ] + " "
    + C[ ioffc + 1 + 0 * ldc ] + " "
    + C[ ioffc + 2 + 0 * ldc ] + " "
    + C[ ioffc + 3 + 0 * ldc ] + " "
    + C[ ioffc + 0 + 1 * ldc ] + " "
    + C[ ioffc + 1 + 1 * ldc ] + " "
    + C[ ioffc + 2 + 1 * ldc ] + " "
    + C[ ioffc + 3 + 1 * ldc ] + " "
    + C[ ioffc + 0 + 2 * ldc ] + " "
    + C[ ioffc + 1 + 2 * ldc ] + " "
    + C[ ioffc + 2 + 2 * ldc ] + " "
    + C[ ioffc + 3 + 2 * ldc ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = 2 * i + j - 1;
    }
  }
  LaPack2.dorm2l('R','T',m,n,k,A,lda,tau,C,ldc,work,info,
    ioffa,iofftau,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dorm2l('R','T',...): C = "
    + C[ ioffc + 0 + 0 * ldc ] + " "
    + C[ ioffc + 1 + 0 * ldc ] + " "
    + C[ ioffc + 2 + 0 * ldc ] + " "
    + C[ ioffc + 3 + 0 * ldc ] + " "
    + C[ ioffc + 0 + 1 * ldc ] + " "
    + C[ ioffc + 1 + 1 * ldc ] + " "
    + C[ ioffc + 2 + 1 * ldc ] + " "
    + C[ ioffc + 3 + 1 * ldc ] + " "
    + C[ ioffc + 0 + 2 * ldc ] + " "
    + C[ ioffc + 1 + 2 * ldc ] + " "
    + C[ ioffc + 2 + 2 * ldc ] + " "
    + C[ ioffc + 3 + 2 * ldc ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dorm2r() {
  document.getElementById("debug_textarea").value +=
    "testing dorm2r *************" + "\n";
  var lda = 5;
  var ldc = 6;
  var m = 4;
  var n = 3;
  var k = 2;
  var ioffa = 1;
  var ioffc = 2;
  var iofftau = 3;
  var ioffwork = 4;
  var A = new Array( ioffa + lda * k );
  var C = new Array( ioffc + ldc * n );
  var tau = new Array( iofftau + k );
  var work = new Array( ioffwork + Math.max( m, n ) );
  for ( var j = 0; j < k; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = 2 * i + j - 1;
    }
  }
  for ( i = 0; i < k; i ++ ) tau[ iofftau + i ] = 3 * i + 2;
  var info = new IntReference();
  LaPack2.dorm2r('L','N',m,n,k,A,lda,tau,C,ldc,work,info,
    ioffa,iofftau,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dorm2r('L','N',...): C = "
    + C[ ioffc + 0 + 0 * ldc ] + " "
    + C[ ioffc + 1 + 0 * ldc ] + " "
    + C[ ioffc + 2 + 0 * ldc ] + " "
    + C[ ioffc + 3 + 0 * ldc ] + " "
    + C[ ioffc + 0 + 1 * ldc ] + " "
    + C[ ioffc + 1 + 1 * ldc ] + " "
    + C[ ioffc + 2 + 1 * ldc ] + " "
    + C[ ioffc + 3 + 1 * ldc ] + " "
    + C[ ioffc + 0 + 2 * ldc ] + " "
    + C[ ioffc + 1 + 2 * ldc ] + " "
    + C[ ioffc + 2 + 2 * ldc ] + " "
    + C[ ioffc + 3 + 2 * ldc ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = 2 * i + j - 1;
    }
  }
  LaPack2.dorm2r('L','T',m,n,k,A,lda,tau,C,ldc,work,info,
    ioffa,iofftau,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dorm2r('L','T',...): C = "
    + C[ ioffc + 0 + 0 * ldc ] + " "
    + C[ ioffc + 1 + 0 * ldc ] + " "
    + C[ ioffc + 2 + 0 * ldc ] + " "
    + C[ ioffc + 3 + 0 * ldc ] + " "
    + C[ ioffc + 0 + 1 * ldc ] + " "
    + C[ ioffc + 1 + 1 * ldc ] + " "
    + C[ ioffc + 2 + 1 * ldc ] + " "
    + C[ ioffc + 3 + 1 * ldc ] + " "
    + C[ ioffc + 0 + 2 * ldc ] + " "
    + C[ ioffc + 1 + 2 * ldc ] + " "
    + C[ ioffc + 2 + 2 * ldc ] + " "
    + C[ ioffc + 3 + 2 * ldc ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = 2 * i + j - 1;
    }
  }
  LaPack2.dorm2r('R','N',m,n,k,A,lda,tau,C,ldc,work,info,
    ioffa,iofftau,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dorm2r('R','N',...): C = "
    + C[ ioffc + 0 + 0 * ldc ] + " "
    + C[ ioffc + 1 + 0 * ldc ] + " "
    + C[ ioffc + 2 + 0 * ldc ] + " "
    + C[ ioffc + 3 + 0 * ldc ] + " "
    + C[ ioffc + 0 + 1 * ldc ] + " "
    + C[ ioffc + 1 + 1 * ldc ] + " "
    + C[ ioffc + 2 + 1 * ldc ] + " "
    + C[ ioffc + 3 + 1 * ldc ] + " "
    + C[ ioffc + 0 + 2 * ldc ] + " "
    + C[ ioffc + 1 + 2 * ldc ] + " "
    + C[ ioffc + 2 + 2 * ldc ] + " "
    + C[ ioffc + 3 + 2 * ldc ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = 2 * i + j - 1;
    }
  }
  LaPack2.dorm2r('R','T',m,n,k,A,lda,tau,C,ldc,work,info,
    ioffa,iofftau,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dorm2r('R','T',...): C = "
    + C[ ioffc + 0 + 0 * ldc ] + " "
    + C[ ioffc + 1 + 0 * ldc ] + " "
    + C[ ioffc + 2 + 0 * ldc ] + " "
    + C[ ioffc + 3 + 0 * ldc ] + " "
    + C[ ioffc + 0 + 1 * ldc ] + " "
    + C[ ioffc + 1 + 1 * ldc ] + " "
    + C[ ioffc + 2 + 1 * ldc ] + " "
    + C[ ioffc + 3 + 1 * ldc ] + " "
    + C[ ioffc + 0 + 2 * ldc ] + " "
    + C[ ioffc + 1 + 2 * ldc ] + " "
    + C[ ioffc + 2 + 2 * ldc ] + " "
    + C[ ioffc + 3 + 2 * ldc ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dorml2() {
  document.getElementById("debug_textarea").value +=
    "testing dorml2 *************" + "\n";
  var lda = 5;
  var ldc = 6;
  var m = 4;
  var n = 3;
  var k = 2;
  var ioffa = 1;
  var ioffc = 2;
  var iofftau = 3;
  var ioffwork = 4;
  var A = new Array( ioffa + lda * Math.max( m, n ) );
  var C = new Array( ioffc + ldc * n );
  var tau = new Array( iofftau + k );
  var work = new Array( ioffwork + Math.max( m, n ) );
  for ( var j = 0; j < Math.max( m, n ); j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = 2 * i + j - 1;
    }
  }
  for ( i = 0; i < k; i ++ ) tau[ iofftau + i ] = 3 * i + 2;
  var info = new IntReference();
  LaPack2.dorml2('L','N',m,n,k,A,lda,tau,C,ldc,work,info,
    ioffa,iofftau,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dorml2('L','N',...): C = "
    + C[ ioffc + 0 + 0 * ldc ] + " "
    + C[ ioffc + 1 + 0 * ldc ] + " "
    + C[ ioffc + 2 + 0 * ldc ] + " "
    + C[ ioffc + 3 + 0 * ldc ] + " "
    + C[ ioffc + 0 + 1 * ldc ] + " "
    + C[ ioffc + 1 + 1 * ldc ] + " "
    + C[ ioffc + 2 + 1 * ldc ] + " "
    + C[ ioffc + 3 + 1 * ldc ] + " "
    + C[ ioffc + 0 + 2 * ldc ] + " "
    + C[ ioffc + 1 + 2 * ldc ] + " "
    + C[ ioffc + 2 + 2 * ldc ] + " "
    + C[ ioffc + 3 + 2 * ldc ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = 2 * i + j - 1;
    }
  }
  LaPack2.dorml2('L','T',m,n,k,A,lda,tau,C,ldc,work,info,
    ioffa,iofftau,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dorml2('L','T',...): C = "
    + C[ ioffc + 0 + 0 * ldc ] + " "
    + C[ ioffc + 1 + 0 * ldc ] + " "
    + C[ ioffc + 2 + 0 * ldc ] + " "
    + C[ ioffc + 3 + 0 * ldc ] + " "
    + C[ ioffc + 0 + 1 * ldc ] + " "
    + C[ ioffc + 1 + 1 * ldc ] + " "
    + C[ ioffc + 2 + 1 * ldc ] + " "
    + C[ ioffc + 3 + 1 * ldc ] + " "
    + C[ ioffc + 0 + 2 * ldc ] + " "
    + C[ ioffc + 1 + 2 * ldc ] + " "
    + C[ ioffc + 2 + 2 * ldc ] + " "
    + C[ ioffc + 3 + 2 * ldc ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = 2 * i + j - 1;
    }
  }
  LaPack2.dorml2('R','N',m,n,k,A,lda,tau,C,ldc,work,info,
    ioffa,iofftau,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dorml2('R','N',...): C = "
    + C[ ioffc + 0 + 0 * ldc ] + " "
    + C[ ioffc + 1 + 0 * ldc ] + " "
    + C[ ioffc + 2 + 0 * ldc ] + " "
    + C[ ioffc + 3 + 0 * ldc ] + " "
    + C[ ioffc + 0 + 1 * ldc ] + " "
    + C[ ioffc + 1 + 1 * ldc ] + " "
    + C[ ioffc + 2 + 1 * ldc ] + " "
    + C[ ioffc + 3 + 1 * ldc ] + " "
    + C[ ioffc + 0 + 2 * ldc ] + " "
    + C[ ioffc + 1 + 2 * ldc ] + " "
    + C[ ioffc + 2 + 2 * ldc ] + " "
    + C[ ioffc + 3 + 2 * ldc ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = 2 * i + j - 1;
    }
  }
  LaPack2.dorml2('R','T',m,n,k,A,lda,tau,C,ldc,work,info,
    ioffa,iofftau,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dorml2('R','T',...): C = "
    + C[ ioffc + 0 + 0 * ldc ] + " "
    + C[ ioffc + 1 + 0 * ldc ] + " "
    + C[ ioffc + 2 + 0 * ldc ] + " "
    + C[ ioffc + 3 + 0 * ldc ] + " "
    + C[ ioffc + 0 + 1 * ldc ] + " "
    + C[ ioffc + 1 + 1 * ldc ] + " "
    + C[ ioffc + 2 + 1 * ldc ] + " "
    + C[ ioffc + 3 + 1 * ldc ] + " "
    + C[ ioffc + 0 + 2 * ldc ] + " "
    + C[ ioffc + 1 + 2 * ldc ] + " "
    + C[ ioffc + 2 + 2 * ldc ] + " "
    + C[ ioffc + 3 + 2 * ldc ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dormr2() {
  document.getElementById("debug_textarea").value +=
    "testing dormr2 *************" + "\n";
  var lda = 5;
  var ldc = 6;
  var m = 4;
  var n = 3;
  var k = 2;
  var ioffa = 1;
  var ioffc = 2;
  var iofftau = 3;
  var ioffwork = 4;
  var A = new Array( ioffa + lda * Math.max( m, n ) );
  var C = new Array( ioffc + ldc * n );
  var tau = new Array( iofftau + k );
  var work = new Array( ioffwork + Math.max( m, n ) );
  for ( var j = 0; j < Math.max( m, n ); j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = 2 * i + j - 1;
    }
  }
  for ( i = 0; i < k; i ++ ) tau[ iofftau + i ] = 3 * i + 2;
  var info = new IntReference();
  LaPack2.dormr2('L','N',m,n,k,A,lda,tau,C,ldc,work,info,
    ioffa,iofftau,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dormr2('L','N',...): C = "
    + C[ ioffc + 0 + 0 * ldc ] + " "
    + C[ ioffc + 1 + 0 * ldc ] + " "
    + C[ ioffc + 2 + 0 * ldc ] + " "
    + C[ ioffc + 3 + 0 * ldc ] + " "
    + C[ ioffc + 0 + 1 * ldc ] + " "
    + C[ ioffc + 1 + 1 * ldc ] + " "
    + C[ ioffc + 2 + 1 * ldc ] + " "
    + C[ ioffc + 3 + 1 * ldc ] + " "
    + C[ ioffc + 0 + 2 * ldc ] + " "
    + C[ ioffc + 1 + 2 * ldc ] + " "
    + C[ ioffc + 2 + 2 * ldc ] + " "
    + C[ ioffc + 3 + 2 * ldc ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = 2 * i + j - 1;
    }
  }
  LaPack2.dormr2('L','T',m,n,k,A,lda,tau,C,ldc,work,info,
    ioffa,iofftau,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dormr2('L','T',...): C = "
    + C[ ioffc + 0 + 0 * ldc ] + " "
    + C[ ioffc + 1 + 0 * ldc ] + " "
    + C[ ioffc + 2 + 0 * ldc ] + " "
    + C[ ioffc + 3 + 0 * ldc ] + " "
    + C[ ioffc + 0 + 1 * ldc ] + " "
    + C[ ioffc + 1 + 1 * ldc ] + " "
    + C[ ioffc + 2 + 1 * ldc ] + " "
    + C[ ioffc + 3 + 1 * ldc ] + " "
    + C[ ioffc + 0 + 2 * ldc ] + " "
    + C[ ioffc + 1 + 2 * ldc ] + " "
    + C[ ioffc + 2 + 2 * ldc ] + " "
    + C[ ioffc + 3 + 2 * ldc ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = 2 * i + j - 1;
    }
  }
  LaPack2.dormr2('R','N',m,n,k,A,lda,tau,C,ldc,work,info,
    ioffa,iofftau,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dormr2('R','N',...): C = "
    + C[ ioffc + 0 + 0 * ldc ] + " "
    + C[ ioffc + 1 + 0 * ldc ] + " "
    + C[ ioffc + 2 + 0 * ldc ] + " "
    + C[ ioffc + 3 + 0 * ldc ] + " "
    + C[ ioffc + 0 + 1 * ldc ] + " "
    + C[ ioffc + 1 + 1 * ldc ] + " "
    + C[ ioffc + 2 + 1 * ldc ] + " "
    + C[ ioffc + 3 + 1 * ldc ] + " "
    + C[ ioffc + 0 + 2 * ldc ] + " "
    + C[ ioffc + 1 + 2 * ldc ] + " "
    + C[ ioffc + 2 + 2 * ldc ] + " "
    + C[ ioffc + 3 + 2 * ldc ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = 2 * i + j - 1;
    }
  }
  LaPack2.dormr2('R','T',m,n,k,A,lda,tau,C,ldc,work,info,
    ioffa,iofftau,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dormr2('R','T',...): C = "
    + C[ ioffc + 0 + 0 * ldc ] + " "
    + C[ ioffc + 1 + 0 * ldc ] + " "
    + C[ ioffc + 2 + 0 * ldc ] + " "
    + C[ ioffc + 3 + 0 * ldc ] + " "
    + C[ ioffc + 0 + 1 * ldc ] + " "
    + C[ ioffc + 1 + 1 * ldc ] + " "
    + C[ ioffc + 2 + 1 * ldc ] + " "
    + C[ ioffc + 3 + 1 * ldc ] + " "
    + C[ ioffc + 0 + 2 * ldc ] + " "
    + C[ ioffc + 1 + 2 * ldc ] + " "
    + C[ ioffc + 2 + 2 * ldc ] + " "
    + C[ ioffc + 3 + 2 * ldc ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dormr3() {
  document.getElementById("debug_textarea").value +=
    "testing dormr3 *************" + "\n";
  var lda = 5;
  var ldc = 6;
  var m = 4;
  var n = 3;
  var k = 2;
  var l = 2;
  var ioffa = 1;
  var ioffc = 2;
  var iofftau = 3;
  var ioffwork = 4;
  var A = new Array( ioffa + lda * Math.max( m, n ) );
  var C = new Array( ioffc + ldc * n );
  var tau = new Array( iofftau + k );
  var work = new Array( ioffwork + Math.max( m, n ) );
  for ( var j = 0; j < Math.max( m, n ); j ++ ) {
    for ( var i = 0; i < k; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = 2 * i + j - 1;
    }
  }
  for ( i = 0; i < k; i ++ ) tau[ iofftau + i ] = 3 * i + 2;
  var info = new IntReference();
  LaPack1.dormr3('L','N',m,n,k,l,A,lda,tau,C,ldc,work,info,
    ioffa,iofftau,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dormr3('L','N',...): C = "
    + C[ ioffc + 0 + 0 * ldc ] + " "
    + C[ ioffc + 1 + 0 * ldc ] + " "
    + C[ ioffc + 2 + 0 * ldc ] + " "
    + C[ ioffc + 3 + 0 * ldc ] + " "
    + C[ ioffc + 0 + 1 * ldc ] + " "
    + C[ ioffc + 1 + 1 * ldc ] + " "
    + C[ ioffc + 2 + 1 * ldc ] + " "
    + C[ ioffc + 3 + 1 * ldc ] + " "
    + C[ ioffc + 0 + 2 * ldc ] + " "
    + C[ ioffc + 1 + 2 * ldc ] + " "
    + C[ ioffc + 2 + 2 * ldc ] + " "
    + C[ ioffc + 3 + 2 * ldc ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = 2 * i + j - 1;
    }
  }
  LaPack1.dormr3('L','T',m,n,k,l,A,lda,tau,C,ldc,work,info,
    ioffa,iofftau,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dormr3('L','T',...): C = "
    + C[ ioffc + 0 + 0 * ldc ] + " "
    + C[ ioffc + 1 + 0 * ldc ] + " "
    + C[ ioffc + 2 + 0 * ldc ] + " "
    + C[ ioffc + 3 + 0 * ldc ] + " "
    + C[ ioffc + 0 + 1 * ldc ] + " "
    + C[ ioffc + 1 + 1 * ldc ] + " "
    + C[ ioffc + 2 + 1 * ldc ] + " "
    + C[ ioffc + 3 + 1 * ldc ] + " "
    + C[ ioffc + 0 + 2 * ldc ] + " "
    + C[ ioffc + 1 + 2 * ldc ] + " "
    + C[ ioffc + 2 + 2 * ldc ] + " "
    + C[ ioffc + 3 + 2 * ldc ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = 2 * i + j - 1;
    }
  }
  LaPack1.dormr3('R','N',m,n,k,l,A,lda,tau,C,ldc,work,info,
    ioffa,iofftau,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dormr3('R','N',...): C = "
    + C[ ioffc + 0 + 0 * ldc ] + " "
    + C[ ioffc + 1 + 0 * ldc ] + " "
    + C[ ioffc + 2 + 0 * ldc ] + " "
    + C[ ioffc + 3 + 0 * ldc ] + " "
    + C[ ioffc + 0 + 1 * ldc ] + " "
    + C[ ioffc + 1 + 1 * ldc ] + " "
    + C[ ioffc + 2 + 1 * ldc ] + " "
    + C[ ioffc + 3 + 1 * ldc ] + " "
    + C[ ioffc + 0 + 2 * ldc ] + " "
    + C[ ioffc + 1 + 2 * ldc ] + " "
    + C[ ioffc + 2 + 2 * ldc ] + " "
    + C[ ioffc + 3 + 2 * ldc ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = 2 * i + j - 1;
    }
  }
  LaPack1.dormr3('R','T',m,n,k,l,A,lda,tau,C,ldc,work,info,
    ioffa,iofftau,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dormr3('R','T',...): C = "
    + C[ ioffc + 0 + 0 * ldc ] + " "
    + C[ ioffc + 1 + 0 * ldc ] + " "
    + C[ ioffc + 2 + 0 * ldc ] + " "
    + C[ ioffc + 3 + 0 * ldc ] + " "
    + C[ ioffc + 0 + 1 * ldc ] + " "
    + C[ ioffc + 1 + 1 * ldc ] + " "
    + C[ ioffc + 2 + 1 * ldc ] + " "
    + C[ ioffc + 3 + 1 * ldc ] + " "
    + C[ ioffc + 0 + 2 * ldc ] + " "
    + C[ ioffc + 1 + 2 * ldc ] + " "
    + C[ ioffc + 2 + 2 * ldc ] + " "
    + C[ ioffc + 3 + 2 * ldc ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dpocon() {
  document.getElementById("debug_textarea").value +=
    "testing dpocon *************" + "\n";
  var n = 4;
  var lda = 5;
  var anorm = 4.;
  var ioffa = 1;
  var ioffwork = 2;
  var ioffiwork = 3;
  var A = new Array( ioffa + lda * n );
  var iwork = new Array( ioffiwork + n );
  var work = new Array( ioffwork + 3*n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + j );
    }
  }
  var info = new IntReference();
  LaPack1.dpotrf('U',n,A,lda,info,ioffa);
  var rcond = new NumberReference();
  LaPack2.dpocon('U',n,A,lda,anorm,rcond,work,iwork,info,
    ioffa,ioffwork,ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dpocon('U',...): info,rcond = " + info.value + " "
    + rcond.value + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + j );
    }
  }
  LaPack1.dpotrf('L',n,A,lda,info,ioffa);
  LaPack2.dpocon('L',n,A,lda,anorm,rcond,work,iwork,info,
    ioffa,ioffwork,ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dpocon('L',...): info,rcond = " + info.value + " "
    + rcond.value + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dptrfs() {
  document.getElementById("debug_textarea").value +=
    "testing dptrfs *************" + "\n";
  var n = 4;
  var nrhs = 1;
  var ioffb = 1;
  var ioffberr = 2;
  var ioffd = 3;
  var ioffdf = 4;
  var ioffe = 5;
  var ioffef = 6;
  var ioffferr = 7;
  var ioffwork = 8;
  var ioffx = 9;
  var d = new Array( ioffd + n );
  var df = new Array( ioffd + n );
  var e = new Array( ioffe + n - 1 );
  var ef = new Array( ioffe + n - 1 );
  var berr = new Array( ioffferr + nrhs );
  var ferr = new Array( ioffferr + nrhs );
  var work = new Array( ioffwork + 2 * n );
  var ldb = 5;
  var B = new Array( ioffb + ldb * nrhs );
  var ldx = 6;
  var X = new Array( ioffx + ldx * nrhs );
  for ( var i = 0; i < n; i ++ ) d[ ioffd + i ] = 2.;
  for ( i = 0; i < n - 1; i ++ ) e[ ioffe + i ] = -1.;
  Blas1.dcopy( n, d, 1, df, 1, ioffd, ioffdf );
  Blas1.dcopy( n - 1, e, 1, ef, 1, ioffe, ioffef );
  var info = new IntReference();
  LaPack0.dpttrf( n, df, ef, info, ioffdf, ioffef );
  B[ ioffb + 0 + 0 * ldb ] = 0.;
  B[ ioffb + 1 + 0 * ldb ] = 0.;
  B[ ioffb + 2 + 0 * ldb ] = 0.;
  B[ ioffb + 3 + 0 * ldb ] = 5.;
  Blas1.dcopy( n, B, 1, X, 1, ioffb, ioffx );
  LaPack1.dpttrs( n, nrhs, df, ef, X, ldx, info,
    ioffdf, ioffef, ioffx );
  LaPack2.dptrfs( n, nrhs, d, e, df, ef, B, ldb, X, ldx, ferr, berr,
    work, info, ioffd, ioffe, ioffdf, ioffef, ioffb, ioffx, ioffferr,
    ioffberr, ioffwork );
  document.getElementById("debug_textarea").value +=
    "dptrfs: berr, ferr = " + berr[ ioffberr ] + " "
    + ferr[ ioffferr ] + "\n";
  document.getElementById("debug_textarea").value +=
    "X = "
    + X[ ioffx + 0 + 0 * ldx ] + " "
    + X[ ioffx + 1 + 0 * ldx ] + " "
    + X[ ioffx + 2 + 0 * ldx ] + " "
    + X[ ioffx + 3 + 0 * ldx ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dsteqr() {
  document.getElementById("debug_textarea").value +=
    "testing dsteqr *************" + "\n";
  var n = 64;
  var ldz = 65;
  var ioffd = 1;
  var ioffe = 2;
  var ioffwork = 3;
  var ioffz = 4;
  var d = new Array( ioffd + n );
  var e = new Array( ioffd + n - 1 );
  var work = new Array( ioffwork + 2 * n - 2 );
  var Z = new Array( ioffz + ldz * n );
  for ( var i = 0; i < n; i ++ ) d[ ioffd + i ] = 2.;
  for ( i = 0; i < n; i ++ ) e[ ioffe + i ] = -1.;
  e[ ioffe + 32 ] = 0.;
  var info = new IntReference();
  LaPack2.dsteqr('N',n,d,e,Z,ldz,work,info,ioffd,ioffe,ioffz,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dsteqr('N',...): info = " + info.value + "\n";
  for ( i = 0; i < n; i ++ ) {
    document.getElementById("debug_textarea").value +=
      "d[" + i + "] = " + d[ ioffd + i ] + "\n";
  }

  for ( i = 0; i < n; i ++ ) d[ ioffd + i ] = 2.;
  for ( i = 0; i < n; i ++ ) e[ ioffe + i ] = -1.;
  e[ ioffe + 32 ] = 0.;
  LaPack2.dsteqr('I',n,d,e,Z,ldz,work,info,ioffd,ioffe,ioffz,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dsteqr('I',...): info = " + info.value + "\n";
  for ( i = 0; i < n; i ++ ) {
    document.getElementById("debug_textarea").value +=
      "d[" + i + "] = " + d[ ioffd + i ] + "\n";
  }
  for ( i = 0; i < n; i ++ ) {
    document.getElementById("debug_textarea").value +=
      "Z[" + i + ",0] = " + Z[ ioffz + i ] + "\n";
  }

  for ( i = 0; i < n; i ++ ) d[ ioffd + i ] = 2.;
  for ( i = 0; i < n; i ++ ) e[ ioffe + i ] = -1.;
  for ( var j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) Z[ ioffz + i + j * ldz ] = 0.;
    Z[ ioffz + j + j * ldz ] = 1.;
  }
  e[ ioffe + 32 ] = 0.;
  LaPack2.dsteqr('V',n,d,e,Z,ldz,work,info,ioffd,ioffe,ioffz,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dsteqr('V',...): info = " + info.value + "\n";
  for ( i = 0; i < n; i ++ ) {
    document.getElementById("debug_textarea").value +=
      "d[" + i + "] = " + d[ ioffd + i ] + "\n";
  }
  for ( i = 0; i < n; i ++ ) {
    document.getElementById("debug_textarea").value +=
      "Z[" + i + ",0] = " + Z[ ioffz + i ] + "\n";
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dsterf() {
  document.getElementById("debug_textarea").value +=
    "testing dsterf *************" + "\n";
  var n = 64;
  var ioffd = 1;
  var ioffe = 2;
  var d = new Array( ioffd + n );
  var e = new Array( ioffd + n - 1 );
  for ( var i = 0; i < n; i ++ ) d[ ioffd + i ] = 2.;
  for ( i = 0; i < n; i ++ ) e[ ioffe + i ] = -1.;
  e[ ioffe + 32 ] = 0.;
  var info = new IntReference();
  LaPack2.dsterf(n,d,e,info,ioffd,ioffe);
  document.getElementById("debug_textarea").value +=
    "dsterf: info = " + info.value + "\n";
  for ( i = 0; i < n; i ++ ) {
    document.getElementById("debug_textarea").value +=
      "d[" + i + "] = " + d[ ioffd + i ] + "\n";
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dsytd2() {
  document.getElementById("debug_textarea").value +=
    "testing dsytd2 *************" + "\n";
  var n = 4;
  var lda = 5;
  var ioffa = 1;
  var ioffd = 2;
  var ioffe = 3;
  var iofftau = 4;
  var A = new Array( ioffa + lda * n );
  var d = new Array( ioffd + n );
  var e = new Array( ioffe + n - 1 );
  var tau = new Array( iofftau + n - 1 );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( i + j + 1 );
    }
  }
  var info = new IntReference();
  LaPack2.dsytd2('U',n,A,lda,d,e,tau,info,ioffa,ioffd,ioffe,iofftau);
  document.getElementById("debug_textarea").value +=
    "dsytd2('U',...): info = " + info.value + "\n";
  document.getElementById("debug_textarea").value +=
    "d = "
    + d[ ioffd + 0 ] + " "
    + d[ ioffd + 1 ] + " "
    + d[ ioffd + 2 ] + " "
    + d[ ioffd + 3 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "e = "
    + e[ ioffe + 0 ] + " "
    + e[ ioffe + 1 ] + " "
    + e[ ioffe + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "tau = "
    + tau[ iofftau + 0 ] + " "
    + tau[ iofftau + 1 ] + " "
    + tau[ iofftau + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 3 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 3 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 3 + 2 * lda ] + " "
    + A[ ioffa + 0 + 3 * lda ] + " "
    + A[ ioffa + 1 + 3 * lda ] + " "
    + A[ ioffa + 2 + 3 * lda ] + " "
    + A[ ioffa + 3 + 3 * lda ]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( i + j + 1 );
    }
  }
  LaPack2.dsytd2('L',n,A,lda,d,e,tau,info,ioffa,ioffd,ioffe,iofftau);
  document.getElementById("debug_textarea").value +=
    "dsytd2('L',...): info = " + info.value + "\n";
  document.getElementById("debug_textarea").value +=
    "d = "
    + d[ ioffd + 0 ] + " "
    + d[ ioffd + 1 ] + " "
    + d[ ioffd + 2 ] + " "
    + d[ ioffd + 3 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "e = "
    + e[ ioffe + 0 ] + " "
    + e[ ioffe + 1 ] + " "
    + e[ ioffe + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "tau = "
    + tau[ iofftau + 0 ] + " "
    + tau[ iofftau + 1 ] + " "
    + tau[ iofftau + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 3 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 3 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 3 + 2 * lda ] + " "
    + A[ ioffa + 0 + 3 * lda ] + " "
    + A[ ioffa + 1 + 3 * lda ] + " "
    + A[ ioffa + 2 + 3 * lda ] + " "
    + A[ ioffa + 3 + 3 * lda ]  + "\n";
}
