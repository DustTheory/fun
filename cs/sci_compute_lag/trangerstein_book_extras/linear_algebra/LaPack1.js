function LaPack1() {
}
//************************************************************************
LaPack1.ddisna = function( job, m, n, d, sep, info, ioffd,
ioffsep ) {
  info.setValue( 0 );
  var eigen = ( job.charAt(0).toUpperCase() == 'E' );
  var left = ( job.charAt(0).toUpperCase() == 'L' );
  var right = ( job.charAt(0).toUpperCase() == 'R' );
  var sing = left || right;
  var k = 0;
  if ( eigen ) k = m;
  else if ( sing ) k = Math.min( m, n );
  if ( ! eigen && ! sing ) info.setValue( -1 );
  else if ( m < 0 ) info.setValue( -2 );
  else if ( k < 0 ) info.setValue( -3 );
  else {
    var incr = true;
    var decr = true;
    for ( var i = 1; i <= k - 1; i ++ ) {
      if ( incr ) incr = ( d[ ioffd + i - 1 ] <= d[ ioffd + i ] );
      if ( decr ) decr = ( d[ ioffd + i - 1 ] >= d[ ioffd + i ] );
    }
    if ( sing && k > 0 ) {
      if ( incr ) incr = ( 0. <= d[ ioffd ] );
      if ( decr ) decr = ( d[ ioffd + k - 1 ] >= 0. );
    }
    if ( ! ( incr || decr ) ) info.setValue( -4 );
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'ddisna', -info.getValue() );
    return;
  }
  if ( k == 0 ) return;
  if ( k == 1 ) sep[ ioffsep ] = LaPack0.dlamch( 'O' );
  else {
    var oldgap = Math.abs( d[ ioffd + 1 ] - d[ ioffd ] );
    sep[ioffsep ] = oldgap;
    for ( i = 2; i <= k - 1; i ++ ) {
      var newgap =
        Math.abs( d[ ioffd + i ] - d[ ioffd + i - 1 ] );
      sep[ ioffsep + i - 1 ] = Math.min( oldgap, newgap );
      oldgap = newgap;
    }
    sep[ ioffsep + k - 1 ] = oldgap;
  }
  if ( sing ) {
    if ( ( left && m > n ) || ( right && m < n ) ) {
      if ( incr ) {
        sep[ ioffsep ] = Math.min( sep[ ioffsep ], d[ ioffd ] );
      }
      if ( decr ) {
        sep[ ioffsep + k - 1 ] =
          Math.min( sep[ ioffsep + k - 1 ], d[ ioffd + k - 1 ] );
      }
    }
  }
  var eps = LaPack0.dlamch( 'E' );
  var safmin = LaPack0.dlamch( 'S' );
  var anorm = Math.max( Math.abs( d[ ioffd ] ) ,
    Math.abs( d[ ioffd + k - 1 ] ) );
  var thresh =
    ( anorm == 0. ? eps : Math.max( eps * anorm, safmin ) );
  for ( i = 1; i <= k; i ++ ) {
    sep[ ioffsep + i - 1 ] =
      Math.max( sep[ ioffsep + i - 1 ] , thresh );
  }
}
//************************************************************************
LaPack1.dgbequ = function( m, n, kl, ku, AB, ldab, r, c,
rowcndReference, colcndReference, amaxReference, info ) {
  throw new Error("not programmed: band matrix");
}
LaPack1.zgbequ = function( m, n, kl, ku, AB, ldab, r, c,
rowcndReference, colcndReference, amaxReference, info ) {
  throw new Error("not programmed: complex band matrix");
}
//************************************************************************
LaPack1.dgbequb = function( m, n, kl, ku, AB, ldab, r, c,
rowcndReference, colcndReference, amaxReference, info ) {
  throw new Error("not programmed: band matrix");
}
//************************************************************************
LaPack1.dgbrfs = function( trans, n, kl, ku, nrhs, AB, ldab,
AFB, ldafb, ipiv, B, ldb, X, ldx, ferr, berr, work, iwork, info ) {
  throw new Error("not programmed: band matrix");
}
LaPack1.zgbrfs = function( trans, n, kl, ku, nrhs, AB, ldab,
AFB, ldafb, ipiv, B, ldb, X, ldx, ferr, berr, work, iwork, info ) {
  throw new Error("not programmed: band matrix");
}
//************************************************************************
LaPack1.dgbtrf = function( m, n, kl, ku, AB, ldab, ipiv,
info ) {
  throw new Error("not programmed: band matrix");
}
LaPack1.zgbtrf = function( m, n, kl, ku, AB, ldab, ipiv,
info ) {
  throw new Error("not programmed: band matrix");
}
//************************************************************************
LaPack1.dgebal = function( job, n, A, lda, ilo, ihi, scale,
info, ioffa, ioffscale) {
  var sclfac = 2.;
  var factor = 0.95;
  info.setValue( 0 );
  if ( job.charAt(0).toUpperCase() != 'N'
  && job.charAt(0).toUpperCase() != 'P'
  && job.charAt(0).toUpperCase() != 'S'
  && job.charAt(0).toUpperCase() != 'B' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1 , n ) ) info.setValue( -4 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgebal', -info.getValue() );
    return;
  }
  var k = 1;
  var l = n;
  ilo.setValue( k );
  ihi.setValue( l );
  if ( n == 0 ) return;
  if ( job.charAt(0).toUpperCase() == 'N' ) {
    for ( var i = 1; i <= n; i ++ ) {
      scale[ ioffscale + i - 1 ] = 1.;
    }
    return;
  }
  var iexc = 1;
  var j = -1;
  if ( job.charAt(0).toUpperCase() != 'S' ) {
    var goto50 = true;
    var m = -1;
    while ( true ) { // 20
      var goto80 = false;
      if ( ! goto50 ) {
        scale[ ioffscale + m - 1 ] = j;
        if ( j != m ) {
          Blas1.dswap( l, A, 1, A, 1, ioffa + ( j - 1 ) * lda,
            ioffa + ( m - 1 ) * lda );
          Blas1.dswap( n - k + 1, A, lda, A, lda,
            ioffa + j - 1 + ( k - 1 ) * lda,
            ioffa + m - 1 + ( k - 1 ) * lda );
        }
        goto80 = ( iexc == 2 );
        if ( ! goto80 ) { // 40
          if ( l == 1 ) {
            ilo.setValue( k );
            ihi.setValue( l );
            return;
          }
          l --;
        }
      } // 50
      goto50 = false;
      var goto20 = false;
      var goto90 = false;
      if ( ! goto80 ) {
        for ( j = l; j >= 1; j -- ) {
          var goto70 = false;
          for ( i = 1; i <= l; i ++ ) {
            if ( i == j ) continue;
            if ( A[ ioffa + j - 1 + ( i - 1 ) * lda ] != 0. ) {
              goto70 = true;
              break;
            }
          }
          if ( ! goto70 ) {
            m = l;
            iexc = 1;
            goto20 = true;
            break;
          }
        } // 70
        if ( goto20 ) {
          continue;
        }
        goto90 = true;
      } // 80
      if ( ! goto90 ) k ++; // 90
      goto20 = false;
      for ( j = k; j <= l; j ++ ) {
        var goto110 = false;
        for ( i = k; i <= l; i ++ ) {
          if ( i == j ) continue;
          if ( A[ ioffa + i - 1 + ( j - 1 ) * lda ] != 0. ) {
            goto110 = true;
            break;
          }
        }
        if ( ! goto110 ) {
          m = k;
          iexc = 2;
          goto20 = true;
          break;
        }
      } // 110
    } // end while loop 20
  } // 120: end job != 'S'
  for ( i = k; i <= l; i ++ ) scale[ ioffscale + i - 1 ] = 1.;
  if ( job.charAt(0).toUpperCase() == 'P' ) {
    ilo.setValue( k );
    ihi.setValue( l );
    return;
  }
  var sfmin1 = LaPack0.dlamch( 'S' ) / LaPack0.dlamch( 'P' );
  var sfmax1 = 1. / sfmin1;
  var sfmin2 = sfmin1 * sclfac;
  var sfmax2 = 1. / sfmin2;
  while ( true ) { // 140
    var noconv = false;
    for ( i = k; i <= l; i ++ ) {
      var c = 0.;
      var r = 0.;
      for ( j = k; j <= l; j ++ ) {
        if ( j == i ) continue;
        c += Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ] );
        r += Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
      }
      var ica = Blas1.idamax( l, A, 1, ioffa + ( i - 1 ) * lda );
      var ca =
        Math.abs( A[ ioffa + ica - 1 + ( i - 1 ) * lda ] );
      var ira = Blas1.idamax( n - k + 1, A, lda,
        ioffa + i - 1 + ( k - 1 ) * lda );
      var ra =
        Math.abs( A[ ioffa + i - 1 + ( ira + k - 2 ) * lda ] );
      if ( c == 0. || r == 0. ) continue;
      var g = r / sclfac;
      var f = 1.;
      var s = c + r;
      while ( ! ( c >= g || Math.max( Math.max( f, c ), ca ) >= sfmax2
      || Math.min( Math.min( r, g), ra ) <= sfmin2 ) ) { // 160
        f *= sclfac;
        c *= sclfac;
        ca *= sclfac;
        r /= sclfac;
        g /= sclfac;
        ra /= sclfac;
      } // 170
      g = c / sclfac;
      while ( ! ( g < r || Math.max( r, ra ) >= sfmax2 || // 180
      Math.min( Math.min( f, c ), Math.min( g, ca ) ) <= sfmin2 ) ) {
        f /= sclfac;
        c /= sclfac;
        g /= sclfac;
        ca *= sclfac;
        r *= sclfac;
        ra *= sclfac;
      } // 190
      if ( ( c + r ) >= factor * s ) continue;
      if ( f < 1. && scale[ ioffscale + i - 1 ] < 1. ) {
        if ( f * scale[ ioffscale + i - 1 ] <= sfmin1 ) continue;
      }
      if ( f > 1. && scale[ ioffscale + i - 1 ] > 1. ) {
        if ( scale[ ioffscale + i - 1 ] >= sfmax1 / f ) continue;
      }
      g = 1. / f;
      scale[ ioffscale + i - 1 ] *= f;
      noconv = true;
      Blas1.dscal( n - k + 1, g, A, lda,
        ioffa + i - 1 + ( k - 1 ) * lda );
      Blas1.dscal( l, f, A, 1, ioffa + ( i - 1 ) * lda );
    } // 200
    if ( ! noconv ) break;
  }
  ilo.setValue( k );
  ihi.setValue( l );
}
LaPack1.zgebal = function( job, n, A, ilo, ihi, scale, info) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack1.dgeequ = function( m, n, A, lda, r, c, rowcndReference,
colcndReference, amaxReference, info, ioffa, ioffr, ioffc ) {
  info.setValue( 0 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -4 );
  if ( info.getValue() != 0) {
    Blas2.xerbla( 'dgeequ', -info.getValue() );
    return;
  }
  if ( m == 0 || n == 0 ) {
    rowcndReference.setValue( 1. );
    colcndReference.setValue( 1. );
    amaxReference.setValue( 0. );
    return;
  }
  var smlnum = LaPack0.dlamch( 'S' );
  var bignum = 1. / smlnum;
  for ( var i = 1; i <= m; i ++ ) r[ ioffr + i - 1 ] = 0.;
  for ( var j = 1; j <= n; j ++ ) {
    for ( i = 1; i <=m; i ++ ) {
      r[ ioffr + i - 1 ] = Math.max( r[ ioffr + i - 1 ],
        Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] ) );
    }
  }
  var rcmin = bignum;
  var rcmax = 0.;
  for ( i = 1; i <= m; i ++ ) {
    rcmax = Math.max( rcmax, r[ ioffr + i - 1 ] );
    rcmin = Math.min( rcmin, r[ ioffr + i - 1 ] );
  }
  amaxReference.setValue( rcmax );
  if ( rcmin == 0. ) {
    for ( i = 1; i <= m; i ++ ) {
      if ( r[ ioffr + i - 1 ] == 0. ) {
        info.setValue( i );
        return;
      }
    }
  } else {
    for ( i = 1; i <= m; i ++ ) {
      r[ ioffr + i - 1 ] = 1. / Math.min( Math.max(
        r[ ioffr + i - 1] , smlnum ), bignum );
    }
    rowcndReference.setValue(
      Math.max( rcmin, smlnum ) / Math.min( rcmax, bignum ) );
  }
  for ( j = 1; j <= n; j ++ ) c[ ioffc + j - 1 ] = 0.;
  for ( j = 1; j <= n; j ++ ) {
    for ( i = 1; i <= m; i ++ ) {
      c[ ioffc + j - 1 ] = Math.max( c[ ioffc + j - 1 ],
        Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] )
        * r [ ioffr + i - 1 ] );
    }
  }
  rcmin = bignum;
  rcmax = 0.;
  for ( j = 1; j <= n; j ++ ) {
    rcmin = Math.min( rcmin, c[ ioffc + j - 1 ] );
    rcmax = Math.max( rcmax, c[ ioffc + j - 1 ] );
  }
  if ( rcmin == 0. ) {
    for ( j = 1; j <= n; j ++ ) {
      if ( c[ ioffc + j - 1 ] == 0. ) {
        info.setValue( m + j );
        return;
      }
    }
  } else {
    for ( j = 1; j <= n; j ++ ) {
      c[ ioffc + j - 1 ] = 1. / Math.min( Math.max(
        c[ ioffc + j - 1 ] , smlnum ), bignum );
    }
    colcndReference.setValue(
      Math.max( rcmin, smlnum ) / Math.min( rcmax, bignum ) );
  }
}
LaPack1.zgeequ = function( m, n, A, lda, r, c, rowcndReference,
colcndReference, amaxReference, info ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack1.dgeequb = function( m, n, A, lda, r, c, rowcndReference,
colcndReference, amaxReference, info, ioffa, ioffr, ioffc ) {
  info.setValue( 0 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -4 );
  if ( info.getValue() != 0) {
    Blas2.xerbla( 'dgeequb', -info.getValue() );
    return;
  }
  if ( m == 0 || n == 0 ) {
    rowcndReference.setValue( 1. );
    colcndReference.setValue( 1. );
    amaxReference.setValue( 0. );
    return;
  }
  var smlnum = LaPack0.dlamch( 'S' );
  var bignum = 1. / smlnum;
  var radix = LaPack0.dlamch( 'B' );
  var logrdx = Math.log( radix );
  for ( var i = 1; i <= m; i ++ ) r[ ioffr + i - 1 ] = 0.;
  for ( var j = 1; j <= n; j ++ ) {
    for ( i = 1; i <=m; i ++ ) {
      r[ ioffr + i - 1 ] = Math.max( r[ ioffr + i - 1 ],
        Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] ) );
    }
  }
  for ( i = 1; i <= m; i ++ ) {
    if ( r[ ioffr + i - 1 ] > 0. ) {
      r[ ioffr + i - 1 ] = Math.pow( radix,
        Math.round( Math.log( r[ ioffr + i - 1 ] ) / logrdx ) );
    }
  }
  var rcmin = bignum;
  var rcmax = 0.;
  for ( i = 1; i <= m; i ++ ) {
    rcmax = Math.max( rcmax, r[ ioffr + i - 1 ] );
    rcmin = Math.min( rcmin, r[ ioffr + i - 1 ] );
  }
  amaxReference.setValue( rcmax );
  if ( rcmin == 0. ) {
    for ( i = 1; i <= m; i ++ ) {
      if ( r[ ioffr + i - 1 ] == 0. ) {
        info.setValue( i );
        return;
      }
    }
  } else {
    for ( i = 1; i <= m; i ++ ) {
      r[ ioffr + i - 1 ] = 1. / Math.min( Math.max(
        r[ ioffr + i - 1] , smlnum ), bignum );
    }
    rowcndReference.setValue(
      Math.max( rcmin, smlnum ) / Math.min( rcmax, bignum ) );
  }
  for ( j = 1; j <= n; j ++ ) c[ ioffc + j - 1 ] = 0.;
  for ( j = 1; j <= n; j ++ ) {
    for ( i = 1; i <= m; i ++ ) {
      c[ ioffc + j - 1 ] = Math.max( c[ ioffc + j - 1 ],
        Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] )
        * r [ ioffr + i - 1 ] );
    }
    if ( c[ ioffc + j - 1 ] > 0. ) {
      c[ ioffc + j - 1 ] = Math.pow( radix,
        Math.round( Math.log( c[ ioffc + j - 1 ] ) / logrdx ) );
    }
  }
  rcmin = bignum;
  rcmax = 0.;
  for ( j = 1; j <= n; j ++ ) {
    rcmin = Math.min( rcmin, c[ ioffc + j - 1 ] );
    rcmax = Math.max( rcmax, c[ ioffc + j - 1 ] );
  }
  if ( rcmin == 0. ) {
    for ( j = 1; j <= n; j ++ ) {
      if ( c[ ioffc + j - 1 ] == 0. ) {
        info.setValue( m + j );
        return;
      }
    }
  } else {
    for ( j = 1; j <= n; j ++ ) {
      c[ ioffc + j - 1 ] = 1. / Math.min( Math.max(
        c[ ioffc + j - 1 ] , smlnum ), bignum );
    }
    colcndReference.setValue(
      Math.max( rcmin, smlnum ) / Math.min( rcmax, bignum ) );
  }
}
//************************************************************************
LaPack1.dgesc2 = function( n, A, lda, rhs, ipiv, jpiv,
scaleReference, ioffa, ioffrhs, ioffipiv, ioffjpiv ) {
  var eps = LaPack0.dlamch( 'P' );
  var smlnumReference =
    new NumberReference( LaPack0.dlamch( 'S' ) / eps );
  var bignumReference =
    new NumberReference( 1. / smlnumReference.getValue() );
  LaPack0.dlabad( smlnumReference, bignumReference );
  LaPack0.dlaswp( 1, rhs, lda, 1, n - 1, ipiv, 1, ioffrhs, ioffipiv );
  for ( var i = 1; i <= n - 1; i ++ ) {
    for ( var j = i + 1; j <= n; j ++ ) {
      rhs[ ioffrhs + j - 1 ] -= A[ ioffa + j - 1 + ( i - 1 ) * lda ]
        * rhs[ ioffrhs + i - 1 ];
    }
  }
  scaleReference.setValue( 1. );
  i = Blas1.idamax( n, rhs, 1, ioffrhs );
  if ( 2. * smlnumReference.getValue()
  * Math.abs( rhs[ ioffrhs + i - 1 ] ) >
  Math.abs( A[ ioffa + n - 1 + ( n- 1 ) * lda ] ) ) {
    var temp = ( 1. / 2. ) / Math.abs( rhs[ ioffrhs + i - 1 ] );
    Blas1.dscal( n, temp, rhs, 1, ioffrhs );
    scaleReference.setValue( scaleReference.getValue() * temp );
  }
  for ( i = n; i >= 1; i -- ) {
    temp = 1. / A[ ioffa + i - 1 + ( i - 1 ) * lda ];
    rhs[ ioffrhs + i - 1 ] *= temp;
    for ( j = i + 1; j <= n; j ++ ) {
      rhs[ ioffrhs + i - 1 ] -= rhs[ ioffrhs + j - 1 ]
        * ( A[ ioffa + i - 1 + ( j - 1 ) * lda ] * temp );
    }
  }
  LaPack0.dlaswp( 1, rhs, lda, 1, n - 1, jpiv, -1, ioffrhs, ioffjpiv );
}
//************************************************************************
LaPack1.dgetc2 = function( n, A, lda, ipiv, jpiv, info, ioffa,
ioffipiv, ioffjpiv) {
  info.setValue( 0 );
  var eps = LaPack0.dlamch( 'P' );
  var smlnumReference =
    new NumberReference( LaPack0.dlamch( 'S' ) / eps );
  var bignumReference =
    new NumberReference( 1. / smlnumReference.getValue() );
  LaPack0.dlabad( smlnumReference, bignumReference );
  for ( var i = 1; i <= n - 1; i ++ ) {
    var xmax = 0.;
    for ( var jp = i; jp <= n; jp ++ ) { // loops reversed from LAPACK
      for ( var ip = i; ip <= n; ip ++ ) {
        if ( Math.abs( A[ ioffa + ip - 1 + ( jp - 1 ) * lda ] )
        > xmax ) { // LAPACK used >= not >
          xmax = Math.abs( A[ ioffa + ip - 1 + ( jp - 1 ) * lda ] );
          var ipv = ip;
          var jpv = jp;
        }
      }
    }
    if ( i == 1 ) {
      var smin = Math.max( eps * xmax, smlnumReference.getValue() );
    }
    if ( ipv != i ) {
      Blas1.dswap( n, A, lda, A, lda, ioffa + ipv - 1, ioffa + i - 1 );
    }
    ipiv[ ioffipiv + i - 1 ] = ipv;
    if ( jpv != i ) {
      Blas1.dswap( n, A, 1, A, 1, ioffa + ( jpv - 1 ) * lda,
        ioffa + ( i - 1 ) * lda );
    }
    jpiv[ ioffjpiv + i - 1 ] = jpv;
    if ( Math.abs( A[ ioffa + i - 1 + ( i - 1 ) * lda ] ) < smin ) {
      info.setValue( i );
      A[ ioffa + i - 1 + ( i - 1 ) * lda ] = smin;
    }
    for ( var j = i + 1; j <= n; j ++ ) {
      A[ ioffa + j - 1 + ( i - 1 ) * lda ] /=
        A[ ioffa + i - 1 + ( i - 1 ) * lda ];
    }
    Blas2.dger( n - i, n - i, -1., A, 1, A, lda, A, lda,
      ioffa + i + ( i - 1 ) * lda, ioffa + i - 1 + i * lda,
      ioffa + i + i * lda );
  }
  if ( Math.abs( A[ ioffa + n - 1 + ( n - 1 ) * lda ] ) < smin ) {
    info.setValue( n );
    A[ ioffa + n - 1 + ( n - 1 ) * lda ] = smin;
  }
}
//************************************************************************
LaPack1.dgetf2 = function( m, n, A, lda, ipiv, info, ioffA,
ioffipiv ) { 
  info.setValue( 0 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1 , m ) ) info.setValue( -4 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgetf2' , - info.getValue() );
    return;
  }
  if ( m == 0 || n == 0 ) return;
  var sfmin = LaPack0.dlamch('S');
  for ( var j = 1; j <= Math.min( m, n ); j ++ ) {
    var jp = j - 1 + Blas1.idamax( m - j + 1 , A, 1,
      ioffA + j - 1 + ( j - 1 ) * lda );
    ipiv[ ioffipiv + j -1 ] = jp;
    if ( A[ ioffA + jp - 1 + ( j - 1 ) * lda ] != 0. ) {
      if ( jp != j ) {
        Blas1.dswap( n, A, lda, A, lda, ioffA + j - 1,
          ioffA + jp - 1 );
      }
      if ( j < m ) {
        if ( Math.abs( A[ ioffA + j - 1 + ( j - 1 ) * lda ] )
        >= sfmin ) {
          Blas1.dscal( m - j,
            1. / A[ ioffA + j - 1 + ( j - 1 ) * lda ],
            A, 1, ioffA + j + ( j - 1 ) * lda );
        } else {
          for ( var i = 1; i <= m - j; i ++ ) {
            A[ ioffA + j + i - 1 + ( j - 1 ) * lda ] /=
              A[ ioffA + j - 1 + ( j - 1 ) * lda ];
          }
        }
      }
    } else if ( info.getValue() == 0 ) info.setValue( j );
    if ( j < Math.min( m , n ) ) {
      Blas2.dger( m - j, n - j, -1., A, 1, A, lda, A, lda, 
        ioffA + j + ( j - 1 ) * lda,
        ioffA + j - 1 + j * lda, ioffA + j + j * lda );
    }
  }
}
LaPack1.zgetf2 = function( m, n, A, lda, ipiv, info ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack1.dgetrs = function( trans, n, nrhs, A, lda, ipiv, B,
ldb, info, ioffa, ioffipiv, ioffb) {
//document.getElementById("debug_textarea").value +=
//  "entering dgetrs, n,nrhs,lda,ldb = " + n + " " + nrhs + " " + lda
//  + " " + ldb + "\n";
//document.getElementById("debug_textarea").value += "A = \n";
//for ( var i = 0; i < n; i ++ ) {
//  for ( var j = 0; j < n; j ++ ) {
//    var Aij = A[ ioffa + i + j * lda ];
//    document.getElementById("debug_textarea").value += Aij + " ";
//  }
//  document.getElementById("debug_textarea").value += "\n";
//}
//document.getElementById("debug_textarea").value += "ipiv = ";
//for ( var i = 0; i < n; i ++ ) {
//  var ipivi = ipiv[ ioffipiv + i ];
//  document.getElementById("debug_textarea").value += ipivi + " ";
//}
//document.getElementById("debug_textarea").value += "\n";
//document.getElementById("debug_textarea").value += "B = \n";
//for ( var i = 0; i < n; i ++ ) {
//  for ( var j = 0; j < nrhs; j ++ ) {
//    var Bij = B[ ioffb + i + j * ldb ];
//    document.getElementById("debug_textarea").value += Bij + " ";
//  }
//  document.getElementById("debug_textarea").value += "\n";
//}
  info.setValue( 0 );
  var notran = ( trans.charAt(0).toUpperCase() == 'N' );
  if ( ! notran && trans.charAt(0).toUpperCase() != 'T'
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( nrhs < 0 ) info.setValue( -3 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -5 );
  else if ( ldb < Math.max( 1, n ) ) info.setValue( -8 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgetrs', -info.getValue() );
    return;
  }
  if ( n == 0 || nrhs == 0 ) return;
  if ( notran ) {
    LaPack0.dlaswp( nrhs, B, ldb, 1, n, ipiv, 1, ioffb, ioffipiv );
    Blas3.dtrsm( 'Left', 'Lower', 'No transpose', 'Unit', n, nrhs, 1.,
      A, lda, B, ldb, ioffa, ioffb );
    Blas3.dtrsm( 'Left', 'Upper', 'No transpose', 'Non-unit', n, nrhs,
      1., A, lda, B, ldb, ioffa, ioffb );
  } else {
    Blas3.dtrsm( 'Left', 'Upper', 'Transpose', 'Non-unit', n, nrhs,
      1., A, lda, B, ldb, ioffa, ioffb );
    Blas3.dtrsm( 'Left', 'Lower', 'Transpose', 'Unit', n, nrhs, 1.,
      A, lda, B, ldb, ioffa, ioffb );
    LaPack0.dlaswp( nrhs, B, ldb, 1, n, ipiv, -1, ioffb, ioffipiv );
  }
//document.getElementById("debug_textarea").value +=
//  "leaving dgetrs\n";
//document.getElementById("debug_textarea").value += "B = \n";
//for ( var i = 0; i < n; i ++ ) {
//  for ( var j = 0; j < nrhs; j ++ ) {
//    var Bij = B[ ioffb + i + j * ldb ];
//    document.getElementById("debug_textarea").value += Bij + " ";
//  }
//  document.getElementById("debug_textarea").value += "\n";
//}
}
LaPack1.zgetrs = function( trans, n, nrhs, A, lda, B, ldb,
info ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack1.dggbal = function( job, n, A, lda, B, ldb, ilo, ihi,
lscale, rscale, work, info, ioffa, ioffb, iofflscale, ioffrscale,
ioffwork ) {
  throw new Error("not tested: generalized eigenvalues");
  info.setValue( 0 );
  if ( job.charAt(0).toUpperCase() != 'N'
  && job.charAt(0).toUpperCase() != 'P'
  && job.charAt(0).toUpperCase() != 'S'
  && job.charAt(0).toUpperCase() != 'B' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1 , n ) ) info.setValue( -4 );
  else if ( ldb < Math.max( 1 , n ) ) info.setValue( -5 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dggbal', -info.getValue() );
    return;
  }
  var k = 1;
  var l = n;
  if ( n == 0 ) return;
  if ( job.charAt(0).toUpperCase() == 'N' ) {
    ilo.setValue( 1 );
    ihi.setValue( n );
    for ( var i = 1; i <= n; i ++ ) {
      lscale[ iofflscale + i - 1 ] = 1.;
      rscale[ ioffrscale + i - 1 ] = 1.;
    }
    return;
  }
  if ( k == l ) {
    ilo.setValue( 1 );
    ihi.setValue( 1 );
    for ( i = 1; i <= n; i ++ ) {
      lscale[ iofflscale + i - 1 ] = 1.;
      rscale[ ioffrscale + i - 1 ] = 1.;
    }
    return;
  }
  var iexc = 1;
  var j = -1;
  ilo.setValue( k );
  ihi.setValue( l );
  if ( job.charAt(0).toUpperCase() != 'S' ) {
    var lm1 = -1;
    var iflow = 1;
    var notfirst = true;
    while ( true ) {
      var goto160 = false;
      if ( iflow == 1 ) { // 20
        if ( notfirst ) {
          l = lm1;
          if ( l == 1 ) {
            rscale[ ioffrscale ] = 1.;
            lscale[ iofflscale ] = 1.;
            break;
          }
          notfirst = true;
        } // 30
        lm1 = l - 1;
        for ( i = l; i >= 1; i -- ) { // 80
          var jp1 = -1;
          var found = false;
          for ( j = 1; j <= lm1; j ++ ) { // 40
            jp1 = j + 1;
            if ( A[ ioffa + i - 1 + ( j - 1 ) * lda ] != 0. ||
            B[ ioffb + i - 1 + ( j - 1 ) * ldb ] != 0. ) {
              found = true;
              break;
            }
          } // 40
          if ( ! found ) j = l;
          else { // 50
            found = false;
            for ( j = jp1; j <= l; j ++ ) { // 60
              if ( A[ ioffa + i - 1 + ( j - 1 ) * lda ] != 0. ||
              B[ ioffb + i - 1 + ( j - 1 ) * ldb ] != 0. ) {
                found = true;
                break;
              }
            } // 60
            if ( ! found )  j = jp1 - 1;
          }
          if ( ! found ) { // 70
            var m = l;
            iflow = 1;
            goto160 = true;
            break;
          }
        } // 80
        if ( goto160 ) break;
      }
      if ( ! goto160 ) {
        if ( iflow != 1 ) k ++; // 90
        found = false; // 100
        for ( j = k; j <= l; j ++ ) { // 150
          for ( i = l; i <= lm1; i ++ ) { // 110
            var ip1 = i + 1;
            if ( A[ ioffa + i - 1 + ( j - 1 ) * lda ] != 0. ||
            B[ ioffb + i - 1 + ( j - 1 ) * ldb ] != 0. ) {
              found = true;
              break;
            }
          } // 110
          if ( ! found ) i = l;
          else { // 120
            found = false;
            for ( i = ip1; i <= l; i ++ ) { // 130
              if ( A[ ioffa + i - 1 + ( j - 1 ) * lda ] != 0. ||
              B[ ioffb + i - 1 + ( j - 1 ) * ldb ] != 0. ) {
                found = true;
                break;
              }
            }
            if ( ! found ) i = ip1 - 1;
          } // 140
          m = k;
          iflow = 2;
          goto160 = true;
          break;
        } // 150
      }
      if ( !goto160 ) break;
      lscale[ iofflscale + m - 1 ] = i;
      if ( i != m ) {
        Blas1.dswap( n - k + 1, A, lda, A, lda,
          ioffa + i - 1 + ( k - 1 ) * lda,
          ioffa + m - 1 + ( k - 1 ) * lda );
        Blas1.dswap( n - k + 1, B, ldb, B, ldb,
          ioffb + i - 1 + ( k - 1 ) * ldb,
          ioffb + m - 1 + ( k - 1 ) * ldb );
      } // 170
      rscale[ ioffrscale + m - 1 ] = j;
      if ( j != m ) {
        Blas1.dswap( l, A, 1, A, 1, ioffa + ( j - 1 ) * lda,
          ioffa + ( m - 1 ) * lda );
        Blas1.dswap( l, B, 1, B, 1, ioffb + ( j - 1 ) * ldb,
          ioffb + ( m - 1 ) * ldb );
      } // 180
    }
  } 
  ilo.setValue( k );
  ihi.setValue( l );
  if ( ilo.getValue() == ihi.getValue() ) return;
  if ( job.charAt(0).toUpperCase() == 'P' ) return;
  var nr = ihi.getValue() - ilo.getValue() + 1;
  for ( i = ilo.getValue(); i <= ihi.getValue(); i ++ ) {
    rscale[ ioffrscale + i - 1 ] = 0.;
    lscale[ iofflscale + i - 1 ] = 0.;
    work[ i - 1         ] = 0.;
    work[ i - 1 +     n ] = 0.;
    work[ i - 1 + 2 * n ] = 0.;
    work[ i - 1 + 3 * n ] = 0.;
    work[ i - 1 + 4 * n ] = 0.;
    work[ i - 1 + 5 * n ] = 0.;
  }
  var sclfac = 1.e1;
  var basl = Math.log( sclfac ) / Math.LN10;
  for ( i = ilo.getValue(); i <= ihi.getValue(); i ++ ) {
    for ( j = ilo.getValue(); i <= ihi.getValue(); j ++ ) {
      var tb = B[ ioffb + i - 1 + ( j - 1 ) * ldb ];
      var ta = A[ ioffa + i - 1 + ( j - 1 ) * lda ];
      if ( ta != 0. ) ta =
        ( Math.log( Math.abs( ta ) ) / Math.LN10 ) / basl;
      if ( tb != 0. ) tb =
        ( Math.log( Math.abs( tb ) ) / Math.LN10 ) / basl;
      work[ i - 1 + 4 * n ] -= ta + tb;
      work[ j - 1 + 5 * n ] -= ta + tb;
    }
  }
  var coef = 1. / Number( 2 * nr );
  var coef2 = coef * coef;
  var coef5 = 0.5 * coef2;
  var nrp2 = nr + 2;
  var beta = 0.;
  var pgamma = Number.POSITIVE_INFINITY;
  var it = 1;
  do {
    var gamma =
      Blas1.ddot( nr, work, 1, work, 1, ilo.getValue() - 1 + 4 * n,
        ilo.getValue() - 1 + 4 * n )
     +Blas1.ddot( nr, work, 1, work, 1, ilo.getValue() - 1 + 5 * n,
        ilo.getValue() - 1 + 5 * n )
    var ew = 0.;
    var ewc = 0.;
    for ( i = ilo.getValue(); i <= ihi.getValue(); i ++ ) {
      ew += work[ i - 1 + 4 * n ];
      ewc += work[ i - 1 + 5 * n ];
    }
    gamma = coef * gamma - coef2 * ( ew * ew + ewc * ewc )
      - coef5 * Math.pow( ew - ewc, 2 );
    if ( gamma == 0. ) break;
    if ( it != 1 ) beta = gamma / pgamma;
    var t = coef5 * ( ewc - 3. * ew );
    var tc = coef5 * ( ew - 3. * ewc );
    Blas1.dscal( nr, beta, work, 1, ilo.getValue() - 1 );
    Blas1.dscal( nr, beta, work, 1, ilo.getValue() - 1 + n );
    Blas1.daxpy( nr, coef, work, 1, work, 1, ilo.getValue() - 1 + 4 * n,
      ilo.getValue() - 1 + n );
    Blas1.daxpy( nr, coef, work, 1, work, 1, ilo.getValue() - 1 + 5 * n,
      ilo.getValue() - 1 );
    for ( i = ilo.getValue(); i <= ihi.getValue(); i ++ ) {
      work[ i - 1 ] += tc;
      work[ i - 1 + n ] += t;
    }
    for ( i = ilo.getValue(); i <= ihi.getValue(); i ++ ) {
      var kount = 0;
      var sum = 0.;
      for ( j = ilo.getValue(); j <= ihi.getValue(); j ++ ) {
        if ( A[ ioffa + i - 1 + ( j - 1 ) * lda ] != 0. ) {
          kount ++;
          sum += work[ j - 1 ];
        }
        if ( B[ ioffb + i - 1 + ( j - 1 ) * ldb ] != 0. ) {
          kount ++;
          sum += work[ j - 1 ];
        }
        work[ i - 1 + 2 * n ] =
          Number( kount ) * work[ i - 1 + n ] + sum;
      }
    }
    for ( j = ilo.getValue(); j <= ihi.getValue(); j ++ ) {
      kount = 0;
      sum = 0.;
      for ( i = ilo.getValue(); i <= ihi.getValue(); i ++ ) {
        if ( A[ ioffa + i - 1 + ( j - 1 ) * lda ] != 0. ) {
          kount ++;
          sum += work[ i - 1 + n ];
        }
        if ( B[ ioffb + i - 1 + ( j - 1 ) * ldb ] != 0. ) {
          kount ++;
          sum += work[ i - 1 + n ];
        }
        work[ j - 1 + 3 * n ] =
          Number( kount ) * work[ j - 1 ] + sum;
      }
    }
    sum = Blas1.ddot( nr, work, 1, work, 1, ilo.getValue() - 1 + n,
        ilo.getValue() - 1 + 2 * n );
      + Blas1.ddot( nr, work, 1, work, 1, ilo.getValue() - 1,
        ilo.getValue() - 1 + 3 * n );
    var alpha = gamma / sum;
    var cmax = 0.;
    for ( i = ilo.getValue(); i <= ihi.getValue(); i ++ ) {
      var cor = alpha * work[ i - 1 + n ];
      if ( Math.abs( cor ) > cmax ) cmax = Math.abs( cor );
      lscale[ iofflscale + i - 1 ] += cor;
      cor = alpha * work[ i - 1 ];
      if ( Math.abs( cor ) > cmax ) cmax = Math.abs( cor );
      lscale[ iofflscale + i - 1 ] += cor;
    }
    if ( cmax < 0.5 ) break;
    Blas1.daxpy( nr, -alpha, work, 1, work, 1, ilo.getValue() - 1 + 2 * n,
      ilo.getValue() - 1 + 4 * n );
    Blas1.daxpy( nr, -alpha, work, 1, work, 1, ilo.getValue() - 1 + 3 * n,
      ilo.getValue() - 1 + 5 * n );
    pgamma = gamma;
    it ++;
  } while ( it <= nrp2 );
  var sfmin = LaPack0.dlamch( 'S' );
  var sfmax = 1. / sfmin;
  var lsfmin = ( Math.log( sfmin ) / Math.LN10 ) / basl + 1.;
  var lsfmax = Math.log( sfmax ) / Math.LN10 / basl;
  for ( i = ilo.getValue(); i <= ihi.getValue(); i ++ ) {
    var irab = Blas1.idamax( n - ilo.getValue() + 1, A, lda,
      ioffa + i - 1 + ( ilo.getValue() - 1 ) * lda );
    var rab = Math.abs( A[ ioffa + i - 1
                      + ( irab + ilo.getValue() - 2 ) * lda ] );
    irab = Blas1.idamax( n - ilo.getValue() + 1, B, ldb,
      ioffb + i - 1 + ( ilo.getValue() - 1 ) * ldb );
    rab = Math.max( rab, Math.abs( B[ ioffb + i - 1
                 + ( irab + ilo.getValue() - 2 ) * ldb ] ) );
    var lrab = ( Math.log( rab + sfmin ) / Math.LN10 ) / basl + 1.;
    var ir = lscale[ iofflscale + i - 1 ]
      + ( lscale[ iofflscale + i - 1 ] >= 0. ? 0.5 : -0.5 );
    ir = Math.min( Math.max( ir, lsfmin ),
                   Math.min( lsfmax, lsfmax - lrab ) );
    lscale[ iofflscale + i - 1 ] = Math.pow( sclfac, ir );
    var icab = Blas1.idamax( ihi.getValue(), A, 1,
      ioffa + ( i - 1 ) * lda );
    var cab =
      Math.abs( A[ ioffa + icab - 1 + ( i - 1 ) * lda ] );
    icab = Blas1.idamax( ihi.getValue(), B, 1, ioffb + ( i - 1 ) * ldb );
    cab = Math.max( cab,
      Math.abs( B[ ioffb + icab - 1 + ( i - 1 ) * ldb ] ) );
    var lcab = ( Math.log( cab + sfmin ) / Math.LN10 ) / basl + 1.;
    var jc = rscale[ ioffrscale + i - 1 ]
      + ( rscale[ ioffrscale + i - 1 ] >= 0. ? 0.5 : -0.5 );
    jc = Math.min( Math.max( jc, lsfmin ),
                   Math.min( lsfmax, lsfmax - lcab ) );
    rscale[ ioffrscale + i - 1 ] = Math.pow( sclfac, jc );
  }
  for ( i = ilo.getValue(); i <= ihi.getValue(); i ++ ) {
    Blas1.dscal( n - ilo.getValue() + 1, lscale[ iofflscale + i - 1 ],
      A, lda, ioffa + i - 1 + ( ilo.getValue() - 1 ) * lda );
    Blas1.dscal( n - ilo.getValue() + 1, lscale[ iofflscale + i - 1 ],
      B, ldb, ioffb + i - 1 + ( ilo.getValue() - 1 ) * ldb );
  }
  for ( i = ilo.getValue(); i <= ihi.getValue(); i ++ ) {
    Blas1.dscal( ihi.getValue(), rscale[ ioffrscale + i - 1 ], A, lda,
      ioffa + ( j - 1 ) * lda );
    Blas1.dscal( ihi.getValue(), rscale[ ioffrscale + i - 1 ], B, ldb,
      ioffb + ( j - 1 ) * ldb );
  }
}
LaPack1.zggbal = function( job, n, A, lda, B, ldb, ilo, ihi,
lscale, rscale, work, info ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack1.dgttrs = function( trans, n, nrhs, dl, d, du, du2,
ipiv, B, ldb, info, ioffdl, ioffd, ioffdu, ioffdu2, ioffipiv, ioffb ) {
  info.setValue( 0 );
  var notran = ( trans.charAt(0).toUpperCase() == 'N' );
  if ( ! notran && trans.charAt(0).toUpperCase() != 'T' &&
  trans.charAt(0).toUpperCase() != 'C' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( nrhs < 0 ) info.setValue( -3 );
  else if ( ldb < Math.max( n , 1 ) ) info.setValue( -10 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgttrs', - info.getValue() );
    return;
  }
  if ( n == 0 || nrhs == 0 ) return;
  var itrans = ( notran ? 0 : 1 );
  var nb = ( nrhs == 1 ? 1 : Math.max( 1,
    LaPack0.ilaenv( 1, 'dgttrs', trans, n, nrhs, -1, -1 ) ) );
  if ( nb >= nrhs ) {
    LaPack0.dgtts2( itrans, n, nrhs, dl, d, du, du2, ipiv, B, ldb,
      ioffdl, ioffd, ioffdu, ioffdu2, ioffipiv, ioffb );
  } else {
    for ( var j = 1; j <= nrhs; j += nb ) {
      var jb = Math.min( nrhs - j + 1, nb );
      LaPack0.dgtts2( itrans, n, jb, dl, d, du, du2, ipiv, B, ldb,
        ioffdl, ioffd, ioffdu, ioffdu2, ioffipiv,
        ioffb + ( j - 1 ) * ldb );
    }
  }
/*
  var i = -1;
  var j = -1;
  var temp = Number.POSITIVE_INFINITY;
  if ( notran ) {
    for ( j = 1; j <= nrhs; j ++ ) {
      for ( i = 1; i <= n - 1; i ++ ) {
        if ( ipiv[ ioffipiv + i - 1 ] == i ) {
          B[ ioffb + i + (j - 1 ) * ldb ] -=
            dl[ ioffdl + i - 1 ] *
            B[ ioffb + i - 1 + ( j - 1 ) * ldb ];
        } else {
          temp = B[ ioffb + i - 1 + ( j - 1 ) * ldb ];
          B[ ioffb + i - 1 + ( j - 1 ) * ldb ] =
            B[ ioffb + i + ( j - 1 ) * ldb ];
          B[ ioffb + i + ( j - 1 ) * ldb ] =
            temp - dl[ ioffdl + i - 1 ] *
            B[ ioffb + i - 1 + ( j - 1 ) * ldb ];
        }
      }
      B[ ioffb + n - 1 + ( j - 1 ) * ldb ] /= d[ ioffd + n - 1 ];
      if ( n > 1 ) {
        B[ ioffb + n - 2 + ( j - 1 ) * ldb ] =
          ( B[ ioffb + n - 2 + ( j - 1 ) * ldb ]
          - du[ ioffdu + n - 2 ]
          * B[ ioffb + n - 1 + ( j - 1 ) * ldb ] )
          / d[ ioffd + n - 2 ];
      }
      for ( i = n - 2; i >= 1; i -- ) {
        B[ ioffb + i - 1 + ( j - 1 ) * ldb ] =
          ( B[ ioffb + i - 1 + ( j - 1 ) * ldb ]
          - du[ ioffdu + i - 1 ] * B[ ioffb + i + ( j - 1 ) * ldb ]
          - du2[ ioffdu2 + i - 1 ] *
          B[ ioffb + i + 1 + ( j - 1 ) * ldb ] ) / d[ ioffd + i - 1 ];
      }
    }
  } else {
    for ( j = 1; j <= nrhs; j ++ ) {
      B[ ioffb + ( j - 1 ) * ldb ] /= d[ ioffd ];
      if ( n > 1 ) {
        B[ ioffb + 1 + ( j - 1 ) * ldb ] =
          ( B[ ioffb + 1 + ( j - 1 ) * ldb ]
          - du[ ioffdu ] * B[ ioffb + ( j - 1 ) * ldb ] )
          / d[ ioffd + 1 ];
      }
      for ( i = 3; i <= n; i ++ ) {
        B[ ioffb + i - 1 + ( j - 1 ) * ldb ] =
          ( B[ ioffb + i - 1 + ( j - 1 ) * ldb ]
          - du[ ioffdu + i - 2 ] * B[ ioffb + i - 2 + ( j - 1 ) * ldb ]
          - du2[ ioffdu2 + i - 3 ] *
          B[ ioffb + i - 3 + ( j - 1 ) * ldb ] ) / d[ ioffd + i - 1 ];
      }
      for ( i = n - 1; i >= 1; i -- ) {
        if ( ipiv[ ioffipiv + i - 1 ] == i ) {
          B[ ioffb + i - 1 + ( j - 1 ) * ldb ] -=
            dl[ ioffdl + i - 1 ] * B[ ioffb + i + ( j - 1 ) * ldb ];
        } else {
          temp = B[ ioffb + i + ( j - 1 ) * ldb ];
          B[ ioffb + i + ( j - 1 ) * ldb ] =
            B[ ioffb + i - 1 + ( j - 1 ) * ldb ]
            - dl[ ioffdl + i - 1 ] * temp;
          B[ ioffb + i - 1 + ( j - 1 ) * ldb ] = temp;
        }
      }
    }
  }
*/
}
LaPack1.zgttrs = function( trans, n, nrhs, dl, d, du, du2,
ipiv, B, ldb, info ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack1.dla_gbamv = function( trans, m, n, kl, ku, alpha, AB,
ldab, x, incx, beta, y, incy, ioffab, ioffx, ioffy) {
  throw new Error("not programmed: banded matrix");
}
//************************************************************************
LaPack1.dla_gbrcond = function( trans, n, kl, ku, AB, ldab,
AFB, ldafb, ipiv, cmode, C, info, work, iwork, ioffab, ioffafb, ioffipiv,
ioffc, ioffwork, ioffiwork ) {
  throw new Error("not programmed: banded matrix");
}
//************************************************************************
LaPack1.dla_geamv = function( trans, m, n, alpha, A, lda, x,
incx, beta, y, incy, ioffa, ioffx, ioffy ) {
  var info = new IntReference( 0 );
  if ( ! ( trans == LaPack0.ilatrans('N' ) ||
  trans == LaPack0.ilatrans('T') ||
  trans == LaPack0.ilatrans('C') ) ) {
    info.setValue( 1 );
  } else if ( m < 0 ) info.setValue( 2 );
  else if ( n < 0 ) info.setValue( 3 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( 6 );
  else if ( incx == 0 ) info.setValue( 8 );
  else if ( incy == 0 ) info.setValue( 11 );
  if ( info != 0 ) {
    Blas2.xerbla( 'dla_geamv', info );
    return;
  }
  if ( m == 0 || n == 0 || ( alpha == 0. && beta == 1. ) ) return;
  if ( trans == LaPack0.ilatrans( 'N' ) ) {
    var lenx = n;
    var leny = m;
  } else {
    lenx = m;
    leny = n;
  }
  var kx = ( incx > 0 ? 1 : 1 - ( lenx - 1 ) * incx );
  var ky = ( incy > 0 ? 1 : 1 - ( leny - 1 ) * incy );
  var safe1 = LaPack0.dlamch( 'Safe minimum' );
  safe1 *= n + 1;
  var iy = ky;
  if ( trans == LaPack0.ilatrans( 'N' ) ) {
    for ( var i = 1; i <= leny; i ++ ) {
      if ( beta == 0. ) {
        var symb_zero = true;
        y[ ioffy + iy - 1 ] = 0.;
      } else if ( y[ ioffy + iy - 1 ] == 0. ) {
        symb_zero = true;
      } else {
        symb_zero = false;
        y[ ioffy + iy - 1 ] *= beta;
      }
      if ( alpha != 0. ) {
        var jx = kx;
        for ( var j = 1; j <= lenx; j ++ ) {
          var temp =
            Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
          symb_zero = symb_zero &&
            ( x[ ioffx + jx - 1 ] == 0. || temp == 0. );
          y[ ioffy + iy - 1 ] +=
            alpha * Math.abs( x[ ioffx + jx - 1 ] ) * temp;
          jx += incx;
        }
      }
      if ( ! symb_zero ) {
        y[ ioffy + iy - 1 ] +=
          ( y[ ioffy + iy - 1 ] >= 0. ? safe1 : - safe1 );
      }
      iy += incy;
    }
  } else {
    for ( i = 1; i <= leny; i ++ ) {
      if ( beta == 0. ) {
        symb_zero = true;
        y[ ioffy + iy - 1 ] = 0.;
      } else if ( y[ ioffy + iy - 1 ] == 0. ) {
        symb_zero = true;
      } else {
        symb_zero = false;
        y[ ioffy + iy - 1 ] *= beta;
      }
      if ( alpha != 0. ) {
        jx = kx;
        for ( j = 1; j <= lenx; j ++ ) {
          temp = Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ] );
          symb_zero = symb_zero &&
            ( x[ ioffx + jx - 1 ] == 0. || temp == 0. );
          y[ ioffy + iy - 1 ] +=
            alpha * Math.abs( x[ ioffx + jx - 1 ] ) * temp;
          jx += incx;
        }
      }
      if ( ! symb_zero ) {
        y[ ioffy + iy - 1 ] +=
          ( y[ ioffy + iy - 1 ] >= 0. ? safe1 : - safe1 );
      }
      iy += incy;
    }
  }
}
//************************************************************************
LaPack1.dla_lin_berr = function( n, nz, nrhs, Res, Ayb, berr,
ioffres, ioffayb, ioffberr ) {
  var safe1 = LaPack0.dlamch( 'Safe minimum' );
  safe1 *= nz + 1;
  for ( var j = 1; j <= nrhs; j ++ ) {
    berr[ ioffberr + j - 1 ] = 0.;
    for ( var i = 1; i <= n; i ++ ) {
      if ( Ayb[ ioffayb + i - 1 + ( j - 1 ) * n ] != 0. ) {
        var tmp = ( safe1
          + Math.abs( Res[ ioffres + i - 1 + ( j - 1 ) * n ] ) )
          / Ayb[ ioffayb + i - 1 + ( j - 1 ) * n ];
        berr[ ioffberr + j - 1 ] =
          Math.max( berr[ ioffberr + j - 1 ], tmp );
      }
    }
  }
}
//************************************************************************
LaPack1.dla_porcond = function( uplo, n, A, lda, AF, ldaf,
cmode, c, info, work, iwork, ioffa, ioffaf, ioffc, ioffwork, ioffiwork ) {
  var isave = new Array( 3 );
  info.setValue( 0 );
  if ( n < 0 ) info.setValue( -2 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dla_porcond', -info.getValue() );
    return 0.;
  }
  if ( n == 0 ) return 1.
  var up = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( up ) {
    for ( var i = 1; i <= n; i ++ ) {
      var tmp = 0.;
      if ( cmode == 1 ) {
        for ( var j = 1; j <= i; j ++ ) {
          tmp += Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ] 
            * c[ ioffc + j - 1 ] );
        }
        for ( j = i + 1; j <= n; j ++ ) {
          tmp += Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] 
            * c[ ioffc + j - 1 ] );
        }
      } else if ( cmode == 0 ) {
        for ( j = 1; j <= i; j ++ ) {
          tmp += Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ] ); 
        }
        for ( j = i + 1; j <= n; j ++ ) {
          tmp += Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] ); 
        }
      } else {
        for ( j = 1; j <= i; j ++ ) {
          tmp += Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ] 
            / c[ ioffc + j - 1 ] );
        }
        for ( j = i + 1; j <= n; j ++ ) {
          tmp += Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] 
            / c[ ioffc + j - 1 ] );
        }
      }
      work[ ioffwork + i - 1 + 2 * n ] = tmp;
    }
  } else {
    for ( i = 1; i <= n; i ++ ) {
      tmp = 0.;
      if ( cmode == 1 ) {
        for ( j = 1; j <= i; j ++ ) {
          tmp += Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] 
            * c[ ioffc + j - 1 ] );
        }
        for ( j = i + 1; j <= n; j ++ ) {
          tmp += Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ] 
            * c[ ioffc + j - 1 ] );
        }
      } else if ( cmode == 0 ) {
        for ( j = 1; j <= i; j ++ ) {
          tmp += Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] ); 
        }
        for ( j = i + 1; j <= n; j ++ ) {
          tmp += Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ] ); 
        }
      } else {
        for ( j = 1; j <= i; j ++ ) {
          tmp += Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] 
            / c[ ioffc + j - 1 ] );
        }
        for ( j = i + 1; j <= n; j ++ ) {
          tmp += Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ] 
            / c[ ioffc + j - 1 ] );
        }
      }
      work[ ioffwork + i - 1 + 2 * n ] = tmp;
    }
  }
  var ainvnmReference = new NumberReference( 0. );
  var kase = new IntReference( 0 );
  while ( true ) { // 10
    LaPack0.dlacn2( n, work, work, iwork, ainvnm, kase, isave,
      ioffwork + n, ioffwork, ioffiwork, 0 );
    if ( kase.getValue() != 0 ) {
      if ( kase.getValue() == 2 ) {
        for ( i = 1; i <= n; i ++ ) {
          work[ ioffwork + i - 1 ] *= work[ ioffwork + i - 1 + 2 * n ];
        }
        if ( up ) {
          LaPack0.dpotrs( 'Upper', n, 1, AF, ldaf, work, n, info,
            ioffaf, ioffwork );
        } else {
          LaPack0.dpotrs( 'Lower', n, 1, AF, ldaf, work, n, info,
            ioffaf, ioffwork );
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
        if ( up ) {
          LaPack0.dpotrs( 'Upper', n, 1, AF, ldaf, work, n, info,
            ioffaf, ioffwork );
        } else {
          LaPack0.dpotrs( 'Lower', n, 1, AF, ldaf, work, n, info,
            ioffaf, ioffwork );
        }
        for ( i = 1; i <= n; i ++ ) {
          work[ ioffwork + i - 1 ] *= work[ ioffwork + i - 1 + 2 * n ];
        }
      }
    } else break;
  }
  if ( ainvnm.getValue() != 0. ) return 1. / ainvnm.getValue();
  return 0.;
}
//************************************************************************
//  cf Blas2.dsymv
LaPack1.dla_syamv = function( uplo, n, alpha, A, lda, x, incx,
beta, y, incy, ioffa, ioffx, ioffy ) {
  var info = new IntReference( 0 );
  if ( uplo != LaPack0.ilauplo( 'U' ) &&
  uplo != LaPack0.ilauplo( 'L' ) ) {
    info.setValue( 1 );
  } else if ( n < 0 ) info.setValue( 2 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( 5 );
  else if ( incx == 0 ) info.setValue( 7 );
  else if ( incy == 0 ) info.setValue( 10 );
  if ( info != 0 ) {
    Blas2.xerbla( 'dla_syamv', info );
    return;
  }
  if ( n == 0 || ( alpha == 0. && beta == 1. ) ) return;
  var kx = ( incx > 0 ? 1 : 1 - ( n - 1 ) * incx );
  var ky = ( incy > 0 ? 1 : 1 - ( n - 1 ) * incy );
  var safe1 = LaPack0.dlamch( 'Safe minimum' );
  safe1 *= n + 1;
  var iy = ky;
  if ( uplo == LaPack0.ilauplo( 'U' ) ) {
    for ( var i = 1; i <= n; i ++ ) {
      if ( beta == 0. ) {
        var symb_zero = true;
        y[ ioffy + iy - 1 ] = 0.;
      } else if ( y[ ioffy + iy - 1 ] == 0. ) {
        symb_zero = true;
      } else {
        symb_zero = false;
        y[ ioffy + iy - 1 ] = beta * Math.abs( y[ ioffy + iy - 1 ] );
      }
      var jx = kx;
      if ( alpha != 0. ) {
        for ( var j = 1; j <= i; j ++ ) {
          var temp =
            Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ] );
          symb_zero = symb_zero &&
            ( x[ ioffx + j - 1 ] == 0. || temp == 0. );
          y[ ioffy + iy - 1 ] += alpha *
            Math.abs( x[ ioffx + jx - 1 ] ) * temp;
          jx += incx;
        }
        for ( j = i + 1; j <= n; j ++ ) {
          temp = Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
          symb_zero = symb_zero &&
            ( x[ ioffx + j - 1 ] == 0. || temp == 0. );
          y[ ioffy + iy - 1 ] += alpha *
            Math.abs( x[ ioffx + jx - 1 ] ) * temp;
          jx += incx;
        }
      }
      if ( ! symb_zero ) {
        y[ ioffy + iy - 1 ] +=
          ( y[ ioffy + iy - 1 ] >= 0. ? safe1 : - safe1 );
      }
      iy += incy;
    }
  } else {
    for ( i = 1; i <= n; i ++ ) {
      if ( beta == 0. ) {
        symb_zero = true;
        y[ ioffy + iy - 1 ] = 0.;
      } else if ( y[ ioffy + iy - 1 ] == 0. ) {
        symb_zero = true;
      } else {
        symb_zero = false;
        y[ ioffy + iy - 1 ] = beta * Math.abs( y[ ioffy + iy - 1 ] );
      }
      jx = kx;
      if ( alpha != 0. ) {
        for ( j = 1; j <= i; j ++ ) {
          temp = Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
          symb_zero = symb_zero &&
            ( x[ ioffx + j - 1 ] == 0. || temp == 0. );
          y[ ioffy + iy - 1 ] += alpha *
            Math.abs( x[ ioffx + jx - 1 ] ) * temp;
          jx += incx;
        }
        for ( j = i + 1; j <= n; j ++ ) {
          temp = Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ] );
          symb_zero = symb_zero &&
            ( x[ ioffx + j - 1 ] == 0. || temp == 0. );
          y[ ioffy + iy - 1 ] += alpha *
            Math.abs( x[ ioffx + jx - 1 ] ) * temp;
          jx += incx;
        }
      }
      if ( ! symb_zero ) {
        y[ ioffy + iy - 1 ] +=
          ( y[ ioffy + iy - 1 ] >= 0. ? safe1 : - safe1 );
      }
      iy += incy;
    }
  }
}
//************************************************************************
LaPack1.dla_syrcond = function( uplo, n, A, lda, AF, ldaf,
ipiv, cmode, c, info, work, iwork, ioffa, ioffaf, ioffipiv, ioffc,
ioffwork, ioffiwork ) {
  throw new Error("not tested");
  var isave = new Array( 3 );
  info.setValue( 0 );
  if ( n < 0 ) info.setValue( -2 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dla_syrcond', - info.getValue() );
    return 0.;
  }
  if ( n == 0 ) return 1.;
  var up = false;
  if ( uplo.charAt(0).toUpperCase() == 'U' ) up = true;
  if ( up ) {
    for ( var i = 1; i <= n; i ++ ) {
      var tmp = 0.;
      if ( cmode == 1 ) {
        for ( var j = 1; j <= i; j ++ ) {
          tmp += Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ]
            * c[ ioffc + j - 1 ] );
        }
        for ( j = i + 1; j <= n; j ++ ) {
          tmp += Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ]
            * c[ ioffc + j - 1 ] );
        }
      } else if ( cmode == 0 ) {
        for ( j = 1; j <= i; j ++ ) {
          tmp += Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ] );
        }
        for ( j = i + 1; j <= n; j ++ ) {
          tmp += Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
        }
      } else {
        for ( j = 1; j <= i; j ++ ) {
          tmp += Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ]
            / c[ ioffc + j - 1 ] );
        }
        for ( j = i + 1; j <= n; j ++ ) {
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
        for ( j = 1; j <= i; j ++ ) {
          tmp += Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ]
            * c[ ioffc + j - 1 ] );
        }
        for ( j = i + 1; j <= n; j ++ ) {
          tmp += Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ]
            * c[ ioffc + j - 1 ] );
        }
      } else if ( cmode == 0 ) {
        for ( j = 1; j <= i; j ++ ) {
          tmp += Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
        }
        for ( j = i + 1; j <= n; j ++ ) {
          tmp += Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ] );
        }
      } else {
        for ( j = 1; j <= i; j ++ ) {
          tmp += Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ]
            / c[ ioffc + j - 1 ] );
        }
        for ( j = i + 1; j <= n; j ++ ) {
          tmp += Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ]
            / c[ ioffc + j - 1 ] );
        }
      }
      work[ ioffwork + 2 * n + i - 1 ] = tmp;
    }
  }
  var smlnum = LaPack0.dlamch( 'Safe minimum' );
  var ainvnmReference = new NumberReference( 0. );
  var normin = 'N';
  var kase = new IntReference( 0 );
  while ( true ) {
    LaPack0.dlacn2( n, work, work, iwork, ainvnm, kase, isave,
      ioffwork + n, ioffwork, ioffiwork, 0 );
    if ( kase.getValue() != 0 ) {
      if ( kase.getValue() == 2 ) {
        for ( i = 1; i <= n; i ++ ) {
          work[ ioffwork + i - 1 ] *= work[ ioffwork + 2 * n + i - 1 ];
        }
        if ( up ) {
          LaPack0.dsytrs( 'U', n, 1, AF, ldaf, ipiv, work, n, info,
            ioffaf, ioffipiv, ioffwork );
        } else {
          LaPack0.dsytrs( 'L', n, 1, AF, ldaf, ipiv, work, n, info,
            ioffaf, ioffipiv, ioffwork );
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
        if ( up ) {
          LaPack0.dsytrs( 'U', n, 1, AF, ldaf, ipiv, work, n, info,
            ioffaf, ioffipiv, ioffwork );
        } else {
          LaPack0.dsytrs( 'L', n, 1, AF, ldaf, ipiv, work, n, info,
            ioffaf, ioffipiv, ioffwork );
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
LaPack1.dlaed2 = function( k, n, n1, d, Q, ldq, indxq,
rho, z, dlamda, w, q2, indx, indxc, indxp, coltyp, info, ioffd,
ioffq, ioffindxq, ioffz, ioffdlamda, ioffw, ioffq2, ioffindx, ioffindxc,
ioffindxp, ioffcoltyp ) {
  throw new Error("not tested: complicated input");
  var ctot = new Array(4);
  var psm = new Array(4);
  info.setValue( 0 );
  if ( n < 0 ) info.setValue( -2 );
  else if ( ldq < Math.max( 1, n ) ) info.setValue( -6 );
  else if ( Math.min( 1, n / 2 ) > n1 || n / 2 < n1 ) info.setValue( -3 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlaed2', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  var n2 = n - n1;
  var n1p1 = n1 + 1;
  if ( rho.getValue() < 0. ) Blas1.dscal( n2, -1., z, 1, ioffz + n1p1 - 1 );
  var t = 1. / Math.sqrt( 2. );
  Blas1.dscal( n, t, z, 1, ioffz );
  rho.setValue( Math.abs( 2. * rho.getValue() ) );
  for ( i = n1p1; i <= n; i ++ ) {
    indxq[ ioffindxq + i - 1 ] += n1;
  }
  for ( var i = 1; i <= n; i ++ ) {
    dlamda[ ioffdlamda + i - 1 ] =
      d[ ioffd + indxq[ ioffindxq + i - 1 ] - 1 ];
  }
  LaPack0.dlamrg( n1, n2, dlamda, 1, 1, indxc, ioffdlamda, ioffindxc );
  for ( i = 1; i <= n; i ++ ) {
    indx[ ioffindx + i - 1 ] =
      indxq[ ioffindxq + indxc[ ioffindxc + i - 1 ] - 1 ];
  }
  var imax = Blas1.idamax( n, z, 1, ioffz );
  var jmax = Blas1.idamax( n, d, 1, ioffd );
  var eps = LaPack0.dlamch( 'Epsilon' );
  var tol = 8. * eps * Math.max( Math.abs(
    d[ ioffd + jmax - 1 ] ), Math.abs( z[ ioffz + imax - 1 ] ) )
  if ( rho.getValue() * Math.abs( z[ ioffz + imax - 1 ] ) <= tol ) {
    k.setValue( 0 );
    var iq2 = 1;
    for ( var j = 1; j <= n; j ++ ) {
      i = indx[ ioffindx + j - 1 ];
      Blas1.dcopy( n, Q, 1, q2, 1, ioffq + ( i - 1 ) * ldq,
        ioffq2 + iq2 - 1 );
      dlamda[ ioffdlamda + j - 1 ] = d[ ioffd + i - 1 ];
      iq2 += n;
    }
    LaPack0.dlacpy( 'A', n, n, q2, n, Q, ldq, ioffq2, ioffq );
    Blas1.dcopy( n, dlamda, 1, d, 1, ioffdlamda, ioffd );
    return;
  }
  for ( i = 1; i <= n1; i ++ ) coltyp[ ioffcoltyp + i - 1 ] = 1;
  for ( i = n1p1; i <= n; i ++ ) coltyp[ ioffcoltyp + i - 1 ] = 3;
  k.setValue( 0 );
  var k2 = n + 1;
  var goto100 = false;
  for ( j = 1; j <= n; j ++ ) {
    var nj = indx[ ioffindx + j - 1 ];
    if ( rho.getValue() * Math.abs( z[ ioffz + nj - 1 ] ) <= tol ) {
      k2 --;
      coltyp[ ioffcoltyp + nj - 1 ] = 4;
      indxp[ ioffindxp + k2 - 1 ] = nj;
      if ( j == n ) {
        goto100 = true;
        break;
      }
    } else {
      var pj = nj;
      break;
    }
  }
  if ( ! goto100 ) {
    while ( true ) { // 80
      j ++;
      nj = indx[ ioffindx + j - 1 ];
      if ( j > n ) break;
      if ( rho.getValue() * Math.abs( z[ ioffz + nj - 1 ] ) <= tol ) {
        k2 --;
        coltyp[ ioffcoltyp + nj - 1 ] = 4;
        indxp[ ioffindxp + k2 - 1 ] = nj;
      } else {
        var s = z[ ioffz + pj - 1 ];
        var c = z[ ioffz + nj - 1 ];
        var tau = LaPack0.dlapy2( c, s );
        t = d[ ioffd + nj - 1 ] - d[ ioffd + pj - 1];
        c /= tau;
        s = -s / tau;
        if ( Math.abs( t * c * s ) <= tol ) {
          z[ ioffz + nj - 1 ] = tau;
          z[ ioffz + pj - 1 ] = 0.;
          if ( coltyp[ ioffcoltyp + nj - 1 ] !=
          coltyp[ ioffcoltyp + pj - 1 ] ) {
            coltyp[ ioffcoltyp + nj - 1 ] = 2;
          }
          coltyp[ ioffcoltyp + pj - 1 ] = 4;
          Blas1.drot( n, Q, 1, Q, 1, c, s, ioffq + ( pj - 1 ) * ldq,
            ioffq + ( nj - 1 ) * ldq );
          t = d[ ioffd + pj - 1 ] * c * c + d[ ioffd + nj - 1 ]
            * s * s;
          d[ ioffd + nj - 1 ] = d[ ioffd + pj - 1 ] * s * s
            + d[ ioffd + nj - 1 ] * c * c;
          d[ ioffd + pj - 1 ] = t;
          k2 --;
          i = 1;
          while ( true ) { // 90
            if ( k2 + i <= n ) {
              if ( d[ ioffd + pj - 1 ] <
              d[ ioffd + indxp[ ioffindxp + k2 + i - 1 ] -1 ] ) {
                indxp[ ioffindxp + k2 + i - 2 ] =
                  indxp[ ioffindxp + k2 + i - 1 ];
                indxp[ ioffindxp + k2 + i - 1 ] = pj;
                i ++;
                continue;
              } else {
                indxp[ ioffindxp + k2 + i - 2 ] = pj;
                break;
              }
            } else {
              indxp[ ioffindxp + k2 + i - 2 ] = pj;
              break;
            }
          }
          pj = nj;
        } else {
          k.setValue( k.getValue() + 1 );
          dlamda[ ioffdlamda + k.getValue() - 1 ] = d[ ioffd + pj - 1 ];
          w[ ioffw + k.getValue() - 1 ] = z[ ioffz + pj - 1 ];
          indxp[ ioffindxp + k.getValue() - 1 ] = pj;
          pj = nj;
        }
      }
    }
  } // 100
  k.setValue( k.getValue() + 1 );
  dlamda[ ioffdlamda + k.getValue() - 1 ] = d[ ioffd + pj - 1 ];
  w[ ioffw + k.getValue() - 1 ] = z[ ioffz + pj - 1 ];
  indxp[ ioffindxp + k.getValue() - 1 ] = pj;
  for ( j = 1; j<= 4; j ++ ) ctot[ j - 1 ] = 0;
  for ( j = 1; j <= n; j ++ ) {
    var ct = coltyp[ ioffcoltyp + j - 1 ];
    ctot[ ct - 1 ] ++;
  }
  psm[ 0 ] = 1;
  psm[ 1 ] = 1 + ctot[ 0 ];
  psm[ 2 ] = psm[ 1 ] + ctot[ 1 ];
  psm[ 3 ] = psm[ 2 ] + ctot[ 2 ];
  k.setValue( n - ctot[ 3 ] );
  for ( j = 1; j <= n; j ++ ) {
    var js = indxp[ ioffindxp + j - 1 ];
    ct = coltyp[ ioffcoltyp + js - 1 ];
    indx[ ioffindx + psm[ ct - 1 ] - 1 ] = js;
    indxc[ ioffindxc + psm[ ct - 1 ] - 1 ] = j;
    psm[ ct - 1 ] ++;
  } // 130
  i = 1;
  var iq1 = 1;
  iq2 = 1 + ( ctot[ 0 ] + ctot[ 1 ] ) * n1;
  for ( j = 1; j <= ctot[ 0 ]; j ++ ) {
    js = indx[ ioffindx + i - 1 ];
    Blas1.dcopy( n1, Q, 1, q2, 1, ioffq + ( js - 1 ) * ldq,
      ioffq2 + iq1 - 1 );
    z[ ioffz + i - 1 ] = d[ ioffd + js - 1 ];
    i ++;
    iq1 += n1;
  } // 140
  for ( j = 1; j <= ctot[ 1 ]; j ++ ) {
    js = indx[ ioffindx + i - 1 ];
    Blas1.dcopy( n1, Q, 1, q2, 1, ioffq + ( js - 1 ) * ldq,
      ioffq2 + iq1 - 1 );
    Blas1.dcopy( n2, Q, 1, q2, 1, ioffq + n1 + ( js - 1 ) * ldq,
      ioffq2 + iq2 - 1 );
    z[ ioffz + i - 1 ] = d[ ioffd + js - 1 ];
    i ++;
    iq1 += n1;
    iq2 += n2;
  } // 150
  for ( j = 1; j <= ctot[ 2 ]; j ++ ) {
    js = indx[ ioffindx + i - 1 ];
    Blas1.dcopy( n2, Q, 1, q2, 1, ioffq + n1 + ( js - 1 ) * ldq,
      ioffq2 + iq2 - 1 );
    z[ ioffz + i - 1 ] = d[ ioffd + js - 1 ];
    i ++;
    iq2 += n2;
  } // 160
  iq1 = iq2;
  for ( j = 1; j <= ctot[ 3 ]; j ++ ) {
    js = indx[ ioffindx + i - 1 ];
    Blas1.dcopy( n, Q, 1, q2, 1, ioffq + ( js - 1 ) * ldq,
      ioffq2 + iq2 - 1 );
    iq2 += n;
    z[ ioffz + i - 1 ] = d[ ioffd + js - 1 ];
    i ++;
  } // 170
  if ( k.getValue() < n ) {
    LaPack0.dlacpy( 'A', n, ctot[ 3 ], q2, n, Q, ldq,
      ioffq2 + iq1 - 1, ioffq + k.getValue() * ldq );
    Blas1.dcopy( n - k.getValue(), z, 1, d, 1, ioffz + k.getValue(),
      ioffd + k.getValue() );
  }
  for ( j = 1; j <= 4; j ++ ) {
    coltyp[ ioffcoltyp + j - 1 ] = ctot[ j - 1 ];
  }
}
//************************************************************************
LaPack1.dlaed6 = function( kniter, orgati, rho, d, z, finit,
tauReference, info, ioffd, ioffz) {
  var maxit = 40;
  var dscale = new Array(3);
  var zscale = new Array(3);
  info.setValue( 0 );
  if ( orgati ) {
    var lbd = d[ ioffd + 1 ];
    var ubd = d[ ioffd + 2 ];
  } else {
    lbd = d[ ioffd ];
    ubd = d[ ioffd + 1 ];
  }
  if ( finit < 0. ) lbd = 0.;
  else ubd = 0.;
  var niter = 1;
  tauReference.setValue( 0. );
  if ( kniter == 2 ) {
    if ( orgati ) {
      var temp = ( d[ ioffd + 2 ] - d[ ioffd + 1 ] ) / 2.;
      var c = rho
        + z[ ioffz ] / ( ( d[ ioffd ] - d[ ioffd + 1 ] ) - temp );
      var a = c * ( d[ ioffd + 1 ] + d[ ioffd + 2 ] )
        + z[ ioffz + 1 ] + z[ ioffz + 2 ];
      var b = c * d[ ioffd + 1 ] * d[ ioffd + 2 ]
        + z[ ioffz + 1 ] * d[ ioffd + 2 ]
        + z[ ioffz + 2 ] * d[ ioffd + 1 ];
    } else {
      temp = ( d[ ioffd ] - d[ ioffd + 1 ] ) / 2.;
      c = rho + z[ ioffz + 2 ]
        / ( ( d[ ioffd + 2 ] - d[ ioffd + 1 ] ) - temp );
      a = c * ( d[ ioffd ] + d[ ioffd + 1 ] )
        + z[ ioffz ] + z[ ioffz + 1 ];
      b = c * d[ ioffd ] * d[ ioffd + 1 ]
        + z[ ioffz ] * d[ ioffd + 1 ] + z[ ioffz + 1 ] * d[ ioffd ];
    }
    temp = Math.max( Math.max( Math.abs( a ), Math.abs( b ) ),
      Math.abs( c ) );
    a /= temp;
    b /= temp;
    c /= temp;
    if ( c == 0. ) tauReference.setValue( b / a );
    else if ( a <= 0. ) {
      tauReference.setValue( ( a
        - Math.sqrt( Math.abs( a * a - 4. * b * c ) ) ) / ( 2. * c ) );
    } else {
      tauReference.setValue( 2. * b
                / ( a + Math.sqrt( Math.abs( a * a - 4 * b * c ) ) ) );
    }
    if ( tauReference.getValue() < lbd || tauReference.getValue() > ubd )
    {
      tauReference.setValue( ( lbd + ubd ) / 2. );
    }
    if ( d[ ioffd ] == tauReference.getValue() ||
    d[ ioffd + 1 ] == tauReference.getValue() ||
    d[ ioffd + 2 ] == tauReference.getValue() ) {
      tauReference.setValue( 0. );
    } else {
      temp = finit
        + tauReference.getValue() * z[ ioffz ]
          / ( d[ ioffd ] * ( d[ ioffd ] - tauReference.getValue() ) )
        + tauReference.getValue() * z[ ioffz  + 1 ]
          / ( d[ ioffd  + 1 ] * ( d[ ioffd  + 1 ]
            - tauReference.getValue() ) )
        + tauReference.getValue() * z[ ioffz  + 2 ]
          / ( d[ ioffd  + 2 ] * ( d[ ioffd  + 2 ]
            - tauReference.getValue() ) );
      if ( temp <= 0. ) lbd = tauReference.getValue();
      else ubd = tauReference.getValue();
      if ( Math.abs( finit ) <= Math.abs( temp ) ) {
        tauReference.setValue( 0. );
      }
    }
  }
  var eps = LaPack0.dlamch( 'Epsilon' );
  var base = LaPack0.dlamch( 'Base' );
  var small1 = Math.pow( base,
    Math.round( Math.log( LaPack0.dlamch( 'SafMin' ) )
    / Math.log( base ) / 3. ) );
  var sminv1 = 1. / small1;
  var small2 = small1 * small1;
  var sminv2 = sminv1 * sminv1;
  temp = ( orgati ?
    Math.min( Math.abs( d[ ioffd + 1 ] - tauReference.getValue() ),
              Math.abs( d[ ioffd + 2 ] - tauReference.getValue() ) ) :
    Math.min( Math.abs( d[ ioffd ] - tauReference.getValue() ),
              Math.abs( d[ ioffd + 1 ] - tauReference.getValue() ) ) );
  var scal = false;
  if ( temp <= small1 ) {
    scal = true;
    if ( temp <= small2 ) {
      var sclfac = sminv2;
      var sclinv = small2;
    } else {
      sclfac = sminv1;
      sclinv = small1;
    }
    for ( var i = 1; i <= 3; i ++ ) {
      dscale[ i - 1 ] = d[ ioffd + i - 1 ] * sclfac;
      zscale[ i - 1 ] = z[ ioffz + i - 1 ] * sclfac;
    }
    tauReference.setValue( tauReference.getValue() * sclfac );
    lbd *= sclfac;
    ubd *= sclfac;
  } else {
    for ( i = 1; i <= 3; i ++ ) {
      dscale[ i - 1 ] = d[ ioffd + i - 1 ];
      zscale[ i - 1 ] = z[ ioffz + i - 1 ];
    }
  }
  var fc = 0.;
  var df = 0.;
  var ddf = 0.;
  for ( i = 1; i <= 3; i ++ ) {
    temp = 1. / ( dscale[ i - 1 ] - tauReference.getValue() );
    var temp1 = zscale[ i - 1 ] * temp;
    var temp2 = temp1 * temp;
    var temp3 = temp2 * temp;
    fc += temp1 / dscale[ i - 1 ];
    df += temp2;
    ddf += temp3;
  }
  var f = finit + tauReference.getValue() * fc; 
  var goto60 = ( Math.abs( f ) <= 0. );
  if ( ! goto60 ) {
    if ( f <= 0. ) lbd = tauReference.getValue();
    else ubd = tauReference.getValue();
    var iter = niter + 1;
    for ( niter = iter; niter <= maxit; niter ++ ) {
      if ( orgati ) {
        temp1 = dscale[ 1 ] - tauReference.getValue();
        temp2 = dscale[ 2 ] - tauReference.getValue();
      } else {
        temp1 = dscale[ 0 ] - tauReference.getValue();
        temp2 = dscale[ 1 ] - tauReference.getValue();
      }
      a = ( temp1 + temp2 ) * f - temp1 * temp2 * df;
      b = temp1 * temp2 * f;
      c = f - ( temp1 + temp2 ) * df + temp1 * temp2 * ddf;
      temp = Math.max( Math.max( Math.abs( a ), Math.abs( b ) ),
        Math.abs( c ) );
      a /= temp;
      b /= temp;
      c /= temp;
      if ( c == 0. ) var eta = b / a;
      else if ( a <= 0. ) {
        eta = ( a - Math.sqrt( Math.abs( a * a - 4. * b * c ) ) )
            / ( 2. * c );
      } else {
        eta = 2. * b
            / ( a + Math.sqrt( Math.abs( a * a - 4 * b * c ) ) );
      }
      if ( f * eta >= 0. ) eta = - f / df;
      tauReference.setValue( tauReference.getValue() + eta );
      if ( tauReference.getValue() < lbd ||
      tauReference.getValue() > ubd ) {
        tauReference.setValue( ( lbd + ubd ) / 2. );
      }
      fc = 0.;
      var erretm = 0.;
      df = 0.;
      ddf = 0.;
      for ( i = 1; i <= 3; i ++ ) {
        temp = 1. / ( dscale[ i - 1 ] - tauReference.getValue() );
        temp1 = zscale[ i - 1 ] * temp;
        temp2 = temp1 * temp;
        temp3 = temp2 * temp;
        var temp4 = temp1 / dscale[ i - 1 ];
        fc += temp4;
        erretm += Math.abs( temp4 );
        df += temp2;
        ddf += temp3;
      }
      f = finit + tauReference.getValue() * fc;
      erretm = 8. * (  Math.abs( finit )
        + Math.abs( tauReference.getValue() ) * erretm )
        + Math.abs( tauReference.getValue() ) * df;
      if ( Math.abs( f ) <= eps * erretm ) {
        goto60 = true;
        break;
      }
      if ( f <= 0. ) lbd = tauReference.getValue();
      else ubd = tauReference.getValue();
    } // 50
  }
  if ( ! goto60 ) info.setValue( 1 );
  if ( scal ) tauReference.setValue( tauReference.getValue() * sclinv );
}
//************************************************************************
LaPack1.dlaed8 = function( icompq, k, n, qsiz, d, Q, ldq, indxq,
rhoReference, cutpnt, z, dlamda, Q2, ldq2, w, perm, givptr, Givcol,
Givnum, indxp, indx, info, ioffd, ioffq, ioffindxq, ioffz, ioffdlamda,
ioffq2, ioffw, ioffperm, ioffgivcol, ioffgivnum, ioffindxp, ioffindx ) {
  throw new Error("not tested: complicated input");
  info.setValue( 0 );
  if ( icompq < 0 || icompq > 1 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -3 );
  else if ( icompq ==1 && qsiz < n ) info.setValue( -4 );
  else if ( ldq < Math.max( 1, n ) ) info.setValue( -7 );
  else if ( cutpnt < Math.min( 1, n ) || cutpnt > n ) info.setValue( -10 );
  else if ( ldq2 < Math.max( 1, n ) ) info.setValue( -14 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlaed8', -info.getValue() );
    return;
  }
  givptr.setValue( 0 );
  if ( n == 0 ) return;
  var n1 = cutpnt;
  var n2 = n - n1;
  var n1p1 = n1 + 1;
  if ( rho.getValue() < 0. ) Blas1.dscal( n2, -1, z, 1, ioffz + n1p1 - 1 );
  var t = 1. / Math.sqrt( 2. );
  for ( var j = 1; j <= n; j ++ ) indx[ ioffindx + j - 1 ] = j;
  Blas1.dscal( n, t, z, 1, ioffz );
  rho.setValue( Math.abs( 2. * rho.getValue() ) );
  for ( var i = cutpnt + 1; i <= n; i ++ ) { 
    indxq[ ioffindxq + i - 1 ] += cutpnt;
  }
  for ( i = 1; i <= n; i ++ ) {
    dlamda[ ioffdlamda + i - 1 ] =
      d[ ioffd + indxq[ ioffindxq + i - 1 ] - 1 ];
    w[ ioffw + i - 1 ] = z[ ioffz + indxq[ ioffindxq + i - 1 ] - 1 ];
  }
  i = 1;
  j = cutpnt + 1;
  LaPack0.dlamrg( n1, n2, dlamda, 1, 1, indx, ioffdlamda, ioffindx );
  for ( i = 1; i <= n; i ++ ) {
    d[ ioffd + i - 1 ] =
      dlamda[ ioffdlamda + indx[ ioffindx + i - 1 ] - 1 ];
    z[ ioffz + i - 1 ] = w[ ioffw + indx[ ioffindx + i - 1 ] - 1 ];
  }
  var imax = Blas1.idamax( n, z, 1, ioffz );
  var jmax = Blas1.idamax( n, d, 1, ioffd );
  var eps = LaPack0.dlamch( 'Epsilon' );
  var tol = 8. * eps * Math.abs( d[ ioffd + jmax - 1 ] );
  if ( rho.getValue() * Math.abs( z[ ioffz + imax - 1 ] ) <= tol ) {
    k.setValue( 0 );
    if ( icompq == 0 ) {
      for ( j = 1; j <= n; j ++ ) {
        perm[ ioffperm + j - 1 ] =
          indxq[ ioffindxq + indx[ ioffindx + j - 1 ] - 1 ];
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        perm[ ioffperm + j - 1 ] =
          indxq[ ioffindxq + indx[ ioffindx + j - 1 ] - 1 ];
        Blas1.dcopy( qsiz, Q, 1, Q2, 1,
          ioffq + ( perm[ ioffperm + j - 1 ] - 1 ) * ldq,
          ioffq2 + ( j - 1 ) * ldq2 );
      }
      LaPack0.dlacpy( 'A', qsiz, n, Q2, ldq2, Q, ldq, ioffq2, ioffq );
    }
    return;
  }
  k.setValue( 0 );
  var k2 = n + 1;
  var goto110 = false;
  for ( j = 1; j <= n; j ++ ) {
    if ( rho.getValue() * Math.abs( z[ ioffz + j - 1 ] ) <= tol ) {
      k2 --;
      indxp[ ioffindxp + k2 - 1 ] = j;
      if ( j == n ) {
        goto110 = true;
        break;
      }
    } else {
      var jlam = j;
      break;
    }
  }
  if ( ! goto110 ) {
    while ( true ) { // 80
      j ++;
      if ( j > n ) break;
      if ( rho.getValue() * Math.abs( z[ ioffz + j - 1 ] ) <= tol ) {
        k2 --;
        indxp[ ioffindxp + k2 - 1 ] = j;
      } else {
        var s = z[ ioffz + jlam - 1 ];
        var c = z[ ioffz + j - 1 ];
        var tau = LaPack0.dlapy2( c, s );
        t = d[ ioffd + j  - 1 ] - d[ ioffd + jlam - 1 ];
        c /= tau;
        s = - s / tau;
        if ( Math.abs( t * c * s ) <= tol ) {
          z[ ioffz + j - 1 ] = tau;
          z[ ioffz + jlam - 1 ] = 0.;
          givptr.setValue( givptr.getValue() + 1 );
          Givcol[ ioffgivcol + ( givptr.getValue() - 1 ) * 2 ] =
            indxq[ ioffindxq + indx[ ioffindx + jlam - 1 ] -1 ];
          Givcol[ ioffgivcol + 1 + ( givptr.getValue() - 1 ) * 2 ] =
            indxq[ ioffindxq + indx[ ioffindx + j - 1 ] -1 ];
          Givnum[ ioffgivnum + ( givptr.getValue() - 1 ) * 2 ] = c;
          Givnum[ ioffgivnum + 1 + ( givptr.getValue() - 1 ) * 2 ] = s;
          if ( icompq == 1 ) {
            Blas1.drot( qsiz, Q, 1, Q, 1, c, s,
              ioffq + ( indxq[ ioffindxq
              + indx[ ioffindx + jlam - 1 ] -1 ] -1 ) * ldq,
              ioffq + ( indxq[ ioffindxq
              + indx[ ioffindx + j - 1 ] -1 ] -1 ) * ldq );
          }
          t = d[ ioffd + jlam - 1 ] * c * c
            + d[ ioffd + j - 1 ] * s * s;
          d[ ioffd + j - 1 ] = d[ ioffd + jlam - 1 ] * s * s
            + d[ ioffd + j - 1 ] * c * c;
          d[ ioffd + jlam - 1 ] = t;
          k2 --;
          i = 1;
          while ( true ) {
            if ( k2 + i <= n ) {
              if ( d[ ioffd + jlam - 1 ] <
              d[ ioffd + indxp[ ioffindxp + k2 + i - 1 ] -1 ] ) {
                indxp[ ioffindxp + k2 + i - 2 ] =
                  indxp[ ioffindxp + k2 + i - 1 ];
                indxp[ ioffindxp + k2 + i - 1 ] = jlam;
                i ++;
              } else {
                indxp[ ioffindxp + k2 + i - 2 ] = jlam;
                break;
              }
            } else {
              indxp[ ioffindxp + k2 + i - 2 ] = jlam;
              break;
            }
          }
          jlam = j;
        } else {
          k.setValue( k.getValue() + 1 );
          w[ ioffw + k.getValue() - 1 ] = z[ ioffz + jlam - 1];
          dlamda[ ioffdlamda + k.getValue() - 1 ] = d[ ioffd + jlam - 1 ];
          indxp[ ioffindxp + k.getValue() - 1 ] = jlam;
          jlam = j;
        }
      }
    } // 100
    k.setValue( k.getValue() + 1 );
    w[ ioffw + k.getValue() - 1 ] = z[ ioffz + jlam - 1];
    dlamda[ ioffdlamda + k.getValue() - 1 ] = d[ ioffd + jlam - 1 ];
    indxp[ ioffindxp + k.getValue() - 1 ] = jlam;
  } // 110
  if ( icompq == 0 ) {
    for ( j = 1; j <= n; j ++ ) {
      var jp = indxp[ ioffindxp + j - 1 ];
      dlamda[ ioffdlamda + j - 1 ] = d[ ioffd + jp - 1 ];
      perm[ ioffperm + j - 1 ] =
        indxq[ ioffindxq + indx[ ioffindx + jp - 1 ] - 1 ];
    }
  } else {
    for ( j = 1; j <= n; j ++ ) {
      jp = indxp[ ioffindxp + j - 1 ];
      dlamda[ ioffdlamda + j - 1 ] = d[ ioffd + jp - 1 ];
      perm[ ioffperm + j - 1 ] =
        indxq[ ioffindxq + indx[ ioffindx + jp - 1 ] - 1 ];
      Blas1.dcopy( qsiz, Q, 1, Q2, 1,
        ioffq + ( perm[ ioffperm + j - 1 ] -1 ) * ldq,
        ioffq2 + ( j - 1 ) * ldq2 );
    }
  }
  if ( k.getValue() < n ) {
    if ( icompq == 0 ) {
      Blas1.dcopy( n - k.getValue(), dlamda, 1, d, 1, ioffdlamda + k.getValue(),
        ioffd + k.getValue() );
    } else {
      Blas1.dcopy( n - k.getValue(), dlamda, 1, d, 1, ioffdlamda + k.getValue(),
        ioffd + k.getValue() );
      LaPack0.dlacpy( 'A', qsiz, n - k.getValue(), Q2, ldq2, Q, ldq,
        ioffq2 + k.getValue() * ldq2, ioffq + k.getValue() * ldq );
    }
  }
}
LaPack1.zlaed8 = function( icompq, k, n, qsiz, d, Q, ldq, indxq,
rhoReference, cutpnt, z, dlamda, Q2, ldq2, w, perm, givptr, Givcol,
Givnum, indxp, indx, info, ioffd, ioffq, ioffindxq, ioffz, ioffdlamda,
ioffq2, ioffw, ioffperm, ioffgivcol, ioffgivnum, ioffindxp, ioffindx ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack1.dlagtf = function( n, a, lambda, b, c, tol, d, in2,
info, ioffa, ioffb, ioffc, ioffd, ioffin2) {
  info.setValue( 0 );
  if ( n < 0 ) {
    info.setValue( -1 );
    Blas2.xerbla( 'dlagtf', -info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  a[ ioffa ] -= lambda;
  in2[ ioffin2 + n - 1 ] = 0;
  if ( n == 1 ) {
    if ( a[ ioffa ] == 0. ) in2[ ioffin2 ] = 1;
    return;
  }
  var eps = LaPack0.dlamch( 'Epsilon' );
  var tl = Math.max( tol, eps );
  var scale1 = Math.abs( a[ ioffa ] ) + Math.abs( b[ ioffb ] );
  for ( var k = 1; k <= n - 1; k ++ ) {
    a[ ioffa + k ] -= lambda;
    var scale2 = Math.abs( c[ ioffc + k - 1 ] )
      + Math.abs( a[ ioffa + k ] );
    if ( k < n - 1 ) scale2 += Math.abs( b[ ioffb + k ] );
    var piv1 =
      ( a[ ioffa + k -1 ] == 0. ? 0.
      : Math.abs( a[ ioffa + k - 1 ] ) / scale1 );
    if ( c[ ioffc + k - 1 ] == 0. ) {
      in2[ ioffin2 + k - 1 ] = 0;
      var piv2 = 0.;
      scale1 = scale2;
      if ( k < n - 1 ) d[ ioffd + k - 1 ] = 0.;
    } else {
      piv2 = Math.abs( c[ ioffc + k - 1 ] ) / scale2;
      if ( piv2 <= piv1 ) {
        in2[ ioffin2 + k - 1 ] = 0;
        scale1 = scale2;
        c[ ioffc + k - 1 ] /= a[ ioffa + k - 1 ];
        a[ ioffa + k ] -= c[ ioffc + k - 1 ] * b[ ioffb + k - 1 ];
        if ( k < n - 1 ) d[ ioffd + k - 1 ] = 0.;
      } else {
        in2[ ioffin2 + k - 1 ] = 1;
        var mult = a[ ioffa + k - 1 ] / c[ ioffc + k - 1 ];
        a[ ioffa + k - 1 ] = c[ ioffc + k - 1 ];
        var temp = a[ ioffa + k ];
        a[ ioffa + k ] = b[ ioffb + k - 1 ] - mult * temp;
        if ( k < n - 1 ) {
          d[ ioffd + k - 1 ] = b[ ioffb + k ];
          b[ ioffb + k ] = - mult * d[ ioffd + k - 1 ];
        }
        b[ ioffb + k - 1 ] = temp;
        c[ ioffc + k - 1 ] = mult;
      }
    }
    if ( ( Math.max( piv1, piv2 ) <= tl ) &&
    ( in2[ ioffin2 + n - 1 ] == 0 ) ) {
      in2[ ioffin2 + n - 1 ] = k;
    }
  }
  if ( ( Math.abs( a[ ioffa + n - 1 ] ) <= scale1 * tl ) &&
  ( in2[ ioffin2 + n - 1 ] == 0 ) ) {
    in2[ ioffin2 + n - 1 ] = n;
  }
}
//************************************************************************
LaPack1.dlagts = function( job, n, a, b, c, d, in2, y,
tolReference, info, ioffa, ioffb, ioffc, ioffd, ioffin2, ioffy ) {
  info.setValue( 0 );
  if ( ( Math.abs( job ) > 2 ) || ( job == 0 ) ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlagts', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  var eps = LaPack0.dlamch( 'Epsilon' );
  var sfmin = LaPack0.dlamch( 'Safe minimum' );
  var bignum = 1. / sfmin;
  if ( job < 0 ) {
    if ( tolReference.getValue() <= 0. ) {
      tolReference.setValue( Math.abs( a[ ioffa ] ) );
      if ( n > 1 ) {
        tolReference.setValue( Math.max( tolReference.getValue(),
                                Math.max( Math.abs( a[ ioffa + 1 ] ),
                                          Math.abs( b[ ioffb ] ) ) ));
      }
      for ( var k = 3; k <= n; k ++ ) {
        tolReference.setValue( Math.max(
          Math.max( tolReference.getValue(),
                    Math.abs( a[ ioffa + k - 1 ] ) ),
          Math.max( Math.abs( b[ ioffb + k - 2 ] ),
          Math.abs( d[ ioffd + k - 3 ] ) ) ) );
      }
      tolReference.setValue( tolReference.getValue() * eps );
      if ( tolReference.getValue() == 0. ) tolReference.setValue( eps );
    }
  }
  if ( Math.abs( job ) == 1 ) {
    for ( k = 2; k <= n; k ++ ) {
      if ( in2[ ioffin2 + k - 2 ] == 0 ) {
        y[ ioffy + k - 1 ] -= c[ ioffc + k - 2 ] * y[ ioffy + k - 2 ];
      } else {
        var temp = y[ ioffy + k - 2 ];
        y[ ioffy + k - 2 ] = y[ ioffy + k - 1 ];
        y[ ioffy + k - 1 ] = temp
          - c[ ioffc + k - 2 ] * y[ ioffy + k - 1 ];
      }
    }
    if ( job == 1 ) {
      for ( k = n; k >= 1; k -- ) {
        if ( k <= n - 2 ) {
          temp = y[ ioffy + k - 1 ]
            - b[ ioffb + k - 1 ] * y[ ioffy + k ]
            - d[ ioffd + k - 1 ] * y[ ioffy + k + 1 ];
        } else if ( k == n - 1 ) {
          temp = y[ ioffy + k - 1 ]
            - b[ ioffb + k - 1 ] * y[ ioffy + k ];
        } else temp = y[ ioffy + k - 1 ];
        var ak = a[ ioffa + k - 1 ];
        var absak = Math.abs( ak );
        if ( absak < 1. ) {
          if ( absak < sfmin ) {
            if ( absak == 0. || Math.abs( temp ) * sfmin > absak ) {
              info.setValue( k );
              return;
            } else {
              temp *= bignum;
              ak *= bignum;
            }
          } else if ( Math.abs( temp ) > absak * bignum ) {
            info.setValue( k );
            return;
          }
        }
        y[ ioffy + k - 1 ] = temp / ak;
      }
    } else {
      for ( k = n; k >= 1; k -- ) {
        if ( k <= n - 2 ) {
          temp = y[ ioffy + k - 1 ]
            - b[ ioffb + k - 1 ] * y[ ioffy + k ]
            - d[ ioffd + k - 1 ] * y[ ioffy + k + 1 ];
        } else if ( k == n - 1 ) {
          temp = y[ ioffy + k - 1 ]
            - b[ ioffb + k - 1 ] * y[ ioffy + k ];
        } else temp = y[ ioffy + k - 1 ];
        ak = a[ ioffa + k - 1 ];
        var pert = ( ak >= 0. ? tolReference.getValue() :
          -tolReference.getValue() );
        while ( true ) { // 40
          absak = Math.abs( ak );
          if ( absak < 1. ) {
            if ( absak < sfmin ) {
              if ( absak == 0. || Math.abs( temp ) * sfmin > absak ) {
                ak += pert;
                pert *= 2.;
              } else {
                temp *= bignum;
                ak *= bignum;
                break;
              }
            } else if ( Math.abs( temp ) > absak * bignum ) {
              ak += pert;
              pert *= 2.;
            }
          } else break;
        }
        y[ ioffy + k - 1 ] = temp / ak;
      } // 50
    }
  } else {
    if ( job == 2 ) {
      for ( k = 1; k <= n; k ++ ) {
        if ( k >= 3 ) {
          temp = y[ ioffy + k - 1 ]
            - b[ ioffb + k - 2 ] * y[ ioffy + k - 2 ]
            - d[ ioffd + k - 3 ] * y[ ioffy + k - 3 ];
        } else if ( k == 2 ) {
          temp = y[ ioffy + k - 1 ]
            - b[ ioffb + k - 2 ] * y[ ioffy + k - 2 ];
        } else temp = y[ ioffy + k - 1 ];
        ak = a[ ioffa + k - 1 ];
        absak = Math.abs( ak );
        if ( absak < 1. ) {
          if ( absak < sfmin ) {
            if ( absak == 0. || Math.abs( temp ) * sfmin > absak ) {
              info.setValue( k );
              return;
            } else {
              temp *= bignum;
              ak *= bignum;
            }
          } else if ( Math.abs( temp ) > absak * bignum ) {
            info.setValue( k );
            return;
          }
        }
        y[ ioffy + k - 1 ] = temp / ak;
      }
    } else {
      for ( k = 1; k <= n; k ++ ) {
        if ( k >= 3 ) {
          temp = y[ ioffy + k - 1 ]
            - b[ ioffb + k - 2 ] * y[ ioffy + k - 2 ]
            - d[ ioffd + k - 3 ] * y[ ioffy + k - 3 ];
        } else if ( k == 2 ) {
          temp = y[ ioffy + k - 1 ]
            - b[ ioffb + k - 2 ] * y[ ioffy + k - 2 ];
        } else temp = y[ ioffy + k - 1 ];
        ak = a[ ioffa + k - 1 ];
        pert = ( ak >= 0. ? tolReference.getValue() :
          - tolReference.getValue() );
        while ( true ) { // 70
          absak = Math.abs( ak );
          if ( absak < 1. ) {
            if ( absak < sfmin ) {
              if ( absak == 0. || Math.abs( temp ) * sfmin > absak ) {
                ak += pert;
                pert *= 2.;
              } else {
                temp *= bignum;
                ak *= bignum;
                break;
              }
            } else if ( Math.abs( temp ) > absak * bignum ) {
              ak += pert;
              pert *= 2;
            }
          } else break;
        }
        y[ ioffy + k - 1 ] = temp / ak;
      } // 80
    }
    for ( k = n; k >= 2; k -- ) {
      if ( in2[ ioffin2 + k - 2 ] == 0 ) {
        y[ ioffy + k - 2 ] -= c[ ioffc + k - 2 ] * y[ ioffy + k - 1 ];
      } else {
        temp = y[ ioffy + k - 2 ];
        y[ ioffy + k - 2 ] = y[ ioffy + k - 1 ];
        y[ ioffy + k - 1 ] -= c[ ioffc + k - 2 ] * y[ ioffy + k - 1 ];
      }
    }
  }
}
//************************************************************************
LaPack1.dlaic1 = function( job, j, x, sest, w, gamma,
sestprReference, sReference, cReference, ioffx, ioffw ) {
  var eps = LaPack0.dlamch( 'Epsilon' );
  var alpha = Blas1.ddot( j, x, 1, w, 1, ioffx, ioffw );
  var absalp = Math.abs( alpha );
  var absgam = Math.abs( gamma );
  var absest = Math.abs( sest );
  if ( job == 1 ) {
    if ( sest == 0. ) {
      var s1 = Math.max( absgam, absalp );
      if ( s1 == 0. ) {
        sReference.setValue( 0. );
        cReference.setValue( 1. );
        sestprReference.setValue( 0. );
      } else {
        sReference.setValue( alpha / s1 );
        cReference.setValue( gamma / s1 );
        var tmp = Math.sqrt( sReference.getValue() * sReference.getValue()
                       + cReference.getValue() * cReference.getValue() );
        sReference.setValue( sReference.getValue() / tmp );
        cReference.setValue( cReference.getValue() / tmp );
        sestprReference.setValue( s1 * tmp );
      }
      return;
    } else if ( absgam <= eps * absest ) {
      sReference.setValue( 1. );
      cReference.setValue( 0. );
      tmp = Math.max( absest, absalp );
      s1 = absest / tmp; 
      var s2 = absalp / tmp; 
      sestprReference.setValue( tmp * Math.sqrt( s1 * s1 + s2 * s2 ) );
      return;
    } else if ( absalp <= eps * absest ) {
      s1 = absgam;
      s2 = absest;
      if ( s1 <= s2 ) {
        sReference.setValue( 1. );
        cReference.setValue( 0. );
        sestprReference.setValue( s2 );
      } else {
        sReference.setValue( 0. );
        cReference.setValue( 1. );
        sestprReference.setValue( s1 );
      }
      return;
    } else if ( absest <= eps * absalp || absest <= eps * absgam ) {
      s1 = absgam;
      s2 = absalp;
      if ( s1 <= s2 ) {
        tmp = s1 / s2;
        sReference.setValue( Math.sqrt( 1. + tmp * tmp ) );
        sestprReference.setValue( s2 * sReference.getValue() );
        cReference.setValue( ( gamma / s2 ) / sReference.getValue() );
        sReference.setValue( ( alpha >= 0. ? 1. : -1. )
          / sReference.getValue() );
      } else {
        tmp = s2 / s1;
        cReference.setValue( Math.sqrt( 1. + tmp * tmp ) );
        sestprReference.setValue( s1 * cReference.getValue() );
        sReference.setValue( ( alpha / s1 ) / cReference.getValue() );
        cReference.setValue( ( gamma >= 0. ? 1. : -1. )
          / cReference.getValue() );
      }
      return;
    } else {
      var zeta1 = alpha / absest;
      var zeta2 = gamma / absest;
      var b = ( 1. - zeta1 * zeta1 - zeta2 * zeta2 ) * 0.5;
      cReference.setValue( zeta1 * zeta1 );
      var t = ( b > 0.  ? cReference.getValue()
        / ( b + Math.sqrt( b * b + cReference.getValue() ) )
        : Math.sqrt( b * b + cReference.getValue() ) - b ); 
      var sine = - zeta1 / t;
      var cosine = - zeta2 / ( 1. + t );
      tmp = Math.sqrt( sine * sine + cosine * cosine );
      sReference.setValue( sine / tmp );
      cReference.setValue( cosine / tmp );
      sestprReference.setValue( Math.sqrt( t + 1. ) * absest );
      return;
    }
  } else if ( job == 2 ) {
    if ( sest == 0. ) {
      sestprReference.setValue( 0. );
      if ( Math.max( absgam, absalp ) == 0. ) {
        sine = 1.;
        cosine = 0.;
      } else {
        sine = - gamma;
        cosine = alpha;
      }
      s1 = Math.max( Math.abs( sine ), Math.abs( cosine ) );
      sReference.setValue( sine / s1 );
      cReference.setValue( cosine / s1 );
      tmp = Math.sqrt( sReference.getValue() * sReference.getValue()
                     + cReference.getValue() * cReference.getValue() );
      sReference.setValue( sReference.getValue() / tmp );
      cReference.setValue( cReference.getValue() / tmp );
      return;
    } else if ( absgam <= eps * absest ) {
      sReference.setValue( 0. );
      cReference.setValue( 1. );
      sestprReference.setValue( absgam );
      return;
    } else if ( absalp <= eps * absest ) {
      s1 = absgam;
      s2 = absest;
      if ( s1 <= s2 ) {
        sReference.setValue( 0. );
        cReference.setValue( 1. );
        sestprReference.setValue( s1 );
      } else {
        sReference.setValue( 1. );
        cReference.setValue( 0. );
        sestprReference.setValue( s2 );
      }
      return;
    } else if ( absest <= eps * absalp || absest <= eps * absgam ) {
      s1 = absgam;
      s2 = absalp;
      if ( s1 <= s2 ) {
        tmp = s1 / s2;
        cReference.setValue( Math.sqrt( 1. + tmp * tmp ) );
        sestprReference.setValue( absest
          * ( tmp / cReference.getValue() ) );
        sReference.setValue( - ( gamma / s2 ) / cReference.getValue() );
        cReference.setValue( ( alpha >= 0. ? 1. : -1. )
          / cReference.getValue() );
      } else {
        tmp = s2 / s1;
        sReference.setValue( Math.sqrt( 1. + tmp * tmp ) );
        sestprReference.setValue( absest / sReference.getValue() );
        cReference.setValue( ( alpha / s1 ) / sReference.getValue() );
//      sReference.setValue( - ( gamma >= 0. ? 1. : -1. )
//        / sReference.getValue() );
        if ( gamma >= 0. ) {
          sReference.setValue( -1. / sReference.getValue() );
        } else sReference.setValue( 1. / sReference.getValue() );
      }
      return;
    } else {
      zeta1 = alpha / absest;
      zeta2 = gamma / absest;
      var norma =
        Math.max( 1. + zeta1 * zeta1 + Math.abs( zeta1 * zeta2 ),
                  Math.abs( zeta1 * zeta2 ) + zeta2 * zeta2 );
      var test =
        1. + 2. * ( zeta1 - zeta2 ) * ( zeta1 + zeta2 );
      if ( test >= 0. ) {
        b = ( zeta1 * zeta1 + zeta2 * zeta2 + 1. ) * 0.5;
        cReference.setValue( zeta2 * zeta2 );
        t = cReference.getValue()
          /( b + Math.sqrt( Math.abs( b * b - cReference.getValue() ) ) );
        sine = zeta1 / ( 1. - t );
        cosine = - zeta2 / t;
        sestprReference.setValue(
          Math.sqrt( t + 4. * eps * eps * norma ) * absest );
      } else {
        b = ( zeta2 * zeta2 + zeta1 * zeta1 - 1. ) * 0.5;
        cReference.setValue( zeta1 * zeta1 );
        if ( b >= 0. ) {
          t = - cReference.getValue()
           /( b + Math.sqrt( Math.abs( b*b + cReference.getValue() ) ) );
        } else {
          t = b - Math.sqrt( b * b + cReference.getValue() );
        }
        sine = - zeta1 / t;
        cosine = -zeta2 / ( 1. + t );
        sestprReference.setValue(
          Math.sqrt( 1. + t + 4. * eps * eps * norma ) * absest );
      }
      tmp = Math.sqrt( sine * sine + cosine * cosine );
      sReference.setValue( sine / tmp );
      cReference.setValue( cosine / tmp );
      return;
    }
  }
}
LaPack1.zlaic1 = function( job, j, x, sest, w, gamma,
sedtprReference, s, c ) {
  throw new Error("not programmed: complex array");
}
//************************************************************************
LaPack1.dlaln2 = function( ltrans, na, nw, smin, ca, A, lda, d1,
d2, B, ldb, wr, wi, X, ldx, scaleReference, xnormReference, info, ioffa,
ioffb, ioffx ) {
  throw new Error("not tested: too many test cases");
  var rswap = new Array( false, true, false, true );
  var zswap = new Array( false, false, true, true );
  var ipivot =
    new Array( 1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4, 3, 2, 1 );
  var ci = new Array( 4 );
  var cr = new Array( 4 );
  var smlnum = 2. * LaPack0.dlamch( 'Safe minimum' );
  var bignum = 1. / smlnum;
  var smini = Math.max( smin, smlnum );
  info.setValue( 0 );
  scale.setValue( 1. );
  if ( na == 1 ) {
    if ( nw == 1 ) {
      var csr = ca * A[ ioffa ] - wr * d1;
      var cnorm = Math.abs( csr );
      if ( cnorm < smini ) {
        csr = smini;
        cnorm = smini;
        info.setValue( 1 );
      }
      var bnorm = Math.abs( B[ ioffb ] );
      if ( cnorm < 1. && bnorm > 1. ) {
        if ( bnorm > bignum * cnorm ) scale.setValue( 1. / bnorm );
      }
      X[ ioffx ] = ( B[ ioffb ] * scale.getValue() ) / csr;
      xnorm.setValue( Math.abs( X[ ioffx ] ) );
    } else {
      csr = ca * A[ ioffa ] - wr * d1;
      var csi = - wi * d1;
      cnorm = Math.abs ( csr ) + Math.abs( csi );
      if ( cnorm < smini ) {
        csr = smini;
        csi = 0.;
        cnorm = smini;
        info.setValue( 1 );
      }
      bnorm = Math.abs( B[ ioffb ] ) + Math.abs( B[ ioffb + ldb ] );
      if ( cnorm < 1. && bnorm > 1. ) {
        if ( bnorm > bignum * cnorm ) scale.setValue( 1. / bnorm );
      }
      LaPack0.dladiv( scale.getValue() * B[ ioffb ],
        scale.getValue() * B[ ioffb + ldb ],
        csr, csi, X[ ioffx ] , X[ ioffx + ldx ] );
      xnorm.setValue( Math.abs( X[ ioffx ] )
        + Math.abs( X[ ioffx + ldx ] ) );
    }
  } else {
    cr[ 0 ] = ca * A[ ioffa ] - wr * d1;
    cr[ 3 ] = ca * A[ ioffa + 1 + lda ] - wr * d2;
    if ( ltrans ) {
      cr[ 2 ] = ca * A[ ioffa + 1 ];
      cr[ 1 ] = ca * A[ ioffa + lda ];
    } else {
      cr[ 1 ] = ca * A[ ioffa + 1 ];
      cr[ 2 ] = ca * A[ ioffa + lda ];
    }
    if ( nw == 1 ) {
      var cmax = 0.;
      var icmax = 0;
      for ( var j = 1; j <= 4; j ++ ) {
        if ( Math.abs( cr[ j - 1 ] ) > cmax ) {
          cmax = Math.abs( cr[ j - 1 ] );
          icmax = j;
        }
      }
      if ( cmax < smini ) {
        bnorm = Math.max( Math.abs( B[ ioffb ] ) ,
          Math.abs( B[ ioffb + 1 ] ) );
        if ( smini < 1. && bnorm > 1. ) {
          if ( bnorm > bignum * smini ) scale.setValue( 1. / bnorm );
        }
        var temp = scale.getValue() / smini;
        X[ ioffx ] = temp * B[ ioffb ];
        X[ ioffx + 1 ] = temp * B[ ioffb + 1 ];
        xnorm.setValue( temp * bnorm );
        info.setValue( 1 );
        return;
      }
      var ur11 = cr[ icmax - 1 ];
      var cr21 = cr[ ipivot[ 1 + ( icmax - 1 ) * 4 ] - 1 ]; 
      var ur12 = cr[ ipivot[ 2 + ( icmax - 1 ) * 4 ] - 1 ]; 
      var cr22 = cr[ ipivot[ 3 + ( icmax - 1 ) * 4 ] - 1 ]; 
      var ur11r = 1. / ur11;
      var lr21 = ur11r * cr21;
      var ur22 = cr22 - ur12 * lr21;
      if ( Math.abs( ur22 ) < smini ) {
        ur22 = smini;
        info.setValue( 1 );
      }
      if ( rswap[ icmax - 1 ] ) {
        var br1 = B[ ioffb + 1 ];
        var br2 = B[ ioffb ];
      } else {
        br1 = B[ ioffb ];
        br2 = B[ ioffb + 1 ];
      }
      br2 -= lr21 * br1;
      var bbnd = Math.max( Math.abs( br1 * ( ur22 * ur11r ) ),
        Math.abs( br2 ) );
      if ( bbnd > 1. && Math.abs( ur22 ) < 1. ) {
        if ( bbnd >= bignum * Math.abs( ur22 ) ) {
          scale.setValue( 1. / bbnd );
        }
      }
      var xr2 = ( br2 * scale.getValue() ) / ur22;
      var xr1 = ( scale.getValue() * br1 ) * ur11r
        - xr2 * ( ur11r * ur12 );
      if ( zswap[ icmax - 1 ] ) {
        X[ ioffx ] = xr2;
        X[ ioffx + 1 ] = xr1;
      } else {
        X[ ioffx ] = xr1;
        X[ ioffx + 1 ] = xr2;
      }
      xnorm.setValue( Math.max( Math.abs( xr1 ), Math.abs( xr2 ) ) );
      if ( xnorm.getValue() > 1. && cmax > 1. ) {
        if ( xnorm.getValue() > bignum / cmax ) {
          temp = cmax / bignum;
          X[ ioffx ] *= temp;
          X[ ioffx + 1 ] *= temp;
          xnorm.setValue( xnorm.getValue() * temp );
          scale.setValue( scale.getValue() * temp );
        }
      }
    } else {
      ci[ 0 ] = - wi * d1;
      ci[ 1 ] = 0.;
      ci[ 2 ] = 0.;
      ci[ 3 ] = - wi * d2;
      cmax = 0.;
      icmax = 0;
      for ( j = 1; j <= 4; j ++ ) {
        if ( Math.abs( crv[ j - 1 ] )
        + Math.abs( civ[ j - 1 ] ) > cmax ) {
          cmax = Math.abs( crv[ j - 1 ] ) + Math.abs( civ[ j - 1 ] );
          icmax = j;
        }
      }
      if ( cmax < smini ) {
        bnorm = Math.max( Math.abs( B[ ioffb ] )
          + Math.abs( B[ ioffb + ldb ] ),
          Math.abs( B[ ioffb + 1 ] )
          + Math.abs( B[ ioffb + 1 + ldb ] ) );
        if ( smini < 1. && bnorm > 1. ) {
          if ( bnorm > bignum * smini ) scale.setValue( 1. / bnorm );
        }
        temp = scale.getValue() / smini;
        X[ ioffx ] = temp * B[ ioffb ];
        X[ ioffx + 1 ] = temp * B[ ioffb + 1 ];
        X[ ioffx + ldx ] = temp * B[ ioffb + ldb ];
        X[ ioffx + 1 + ldx ] = temp * B[ ioffb + 1 + ldb ];
        xnorm.setValue( temp * bnorm );
        info.setValue( 1 );
        return;
      }
      ur11 = cr[ icmax - 1 ];
      ui11 = ci[ icmax - 1 ];
      cr21 = cr[ ipivot[ 1 + ( icmax - 1 ) * 4 ] - 1 ]; 
      var ci21 = ci[ ipivot[ 1 + ( icmax - 1 ) * 4 ] - 1 ]; 
      ur12 = cr[ ipivot[ 2 + ( icmax - 1 ) * 4 ] - 1 ]; 
      var ui12 = ci[ ipivot[ 2 + ( icmax - 1 ) * 4 ] - 1 ]; 
      cr22 = cr[ ipivot[ 3 + ( icmax - 1 ) * 4 ] - 1 ]; 
      var ci22 = ci[ ipivot[ 3 + ( icmax - 1 ) * 4 ] - 1 ]; 
      if ( icmax ==1 || icmax == 4 ) {
        if ( Math.abs( ur11 ) > Math.abs( ui11 ) ) {
          temp = ui11 / ur11;
          ur11r = 1. / ( ur11 * ( 1. + temp * temp ) );
          var ui11r = - temp * ur11r;
        } else {
          temp = ur11 / ui11;
          ui11r = -1. / ( ui11 * ( 1. + temp * temp ) );
          ur11r = - temp * ui11r;
        }
        lr21 = cr21 * ur11r;
        var li21 = cr21 * ui11r;
        var ur12s = ur12 * ur11r;
        var ui12s = ur12 * ui11r;
        ur22 = cr22 - ur12 * lr21;
        var ui22 = ci22 - ur12 * li21;
      } else {
        ur11r = 1. / ur11;
        ui11r = 0.;
        lr21 = cr21 * ur11r;
        li21 = ci21 * ur11r;
        ur12s = ur12 * ur11r;
        ui12s = ui12 * ur11r;
        ur22 = cr22 - ur12 * lr21 + ui12 * li21;
        ui22 = - ur12 * li21 - ui12 * lr21;
      }
      var u22abs = Math.abs( ur22 ) + Math.abs( ui22 );
      if ( u22abs < smini ) {
        ur22 = smini;
        ui22 = 0.;
        info.setValue( 1 );
      }
      if ( rswap[ icmax - 1 ] ) {
        br2 = B[ ioffb ];
        br1 = B[ ioffb + 1 ];
        var bi2 = B[ ioffb + ldb ];
        var bi1 = B[ ioffb + 1 + ldb ];
      } else {
        br1 = B[ ioffb ];
        br2 = B[ ioffb + 1 ];
        bi1 = B[ ioffb + ldb ];
        bi2 = B[ ioffb + 1 + ldb ];
      }
      br2 += - lr21 * br1 + li21 * bi1;
      bi2 += - li21 * br1 - lr21 * bi1;
      bbnd = Math.max( ( Math.abs( br1 ) + Math.abs( bi1 ) )
        * ( u22abs * ( Math.abs( ur11r ) + Math.abs( ui11r ) ) ),
        Math.abs( br2 ) + Math.abs( bi2 ) );
      if ( bbnd > 1. && u22abs < 1. ) {
        if ( bbnd >= bignum * u22abs ) {
          scale.setValue( 1. / bbnd );
          br1 *= scale.getValue();
          bi1 *= scale.getValue();
          br2 *= scale.getValue();
          bi2 *= scale.getValue();
        }
      }
      var xi2 = Number.POSITIVE_INFINITY;
      var xr2refReference = new NumberReference( xr2 );
      var xi2refReference = new NumberReference( xi2 );
      LaPack0.dladiv( br2, bi2, ur22, ui22, xr2ref, xi2ref );
      xr2 = xr2ref.getValue();
      xi2 = xi2ref.getValue();
      xr1 = ur11r * br1 - ui11r * bi1 - ur12s * xr2 + ui12s * xi2;
      var xi1 =
        ui11r * br1 + ur11r * bi1 - ui12s * xr2 - ur12s * xi2;
      if ( zswap[ icmax - 1 ] ) {
        X[ ioffx ] = xr2;
        X[ ioffx + 1 ] = xr1;
        X[ ioffx + ldx ] = xi2;
        X[ ioffx + 1 + ldx ] = xi1;
      } else {
        X[ ioffx ] = xr1;
        X[ ioffx + 1 ] = xr2;
        X[ ioffx + ldx ] = xi1;
        X[ ioffx + 1 + ldx ] = xi2;
      }
      xnorm.setValue( Math.max( Math.abs( xr1 ) + Math.abs( xi1 ),
        Math.abs( xr2 ) + Math.abs( xi2 ) ) );
      if ( xnorm.getValue() > 1. && cmax > 1. ) {
        if ( xnorm.getValue() > bignum / cmax ) {
          temp = cmax / bignum;
          X[ ioffx ] *= temp;
          X[ ioffx + 1 ] *= temp;
          X[ ioffx + ldx ] *= temp;
          X[ ioffx + 1 + ldx ] *= temp;
          xnorm.setValue( xnorm.getValue() * temp );
          scale.setValue( scale.getValue() * temp );
        }
      }
    }
  }
}
//************************************************************************
LaPack1.dlangb = function( norm, n, kl, ku, AB, ldab, work ) {
  throw new Error("not programmed: band matrix");
}
LaPack1.zlangb = function( norm, n, kl, ku, AB, ldab, work ) {
  throw new Error("not programmed: band matrix");
}
//************************************************************************
LaPack1.dlange = function( norm, m, n, A, lda, work, ioffa,
ioffwork ) {
  var ret = 0.;
  if ( Math.min( m, n ) == 0 ) return ret;
  if ( norm.charAt(0).toUpperCase() == 'M' ) {
    for ( var j = 1; j <= n; j ++ ) {
      for ( var i = 1; i <= m; i ++ ) {
        ret = Math.max( ret,
          Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] ) );
      }
    }
  } else if ( norm.charAt(0).toUpperCase() == 'O' ||
  norm.charAt(0) == '1' ) {
    for ( j = 1; j <= n; j ++ ) {
      var sum = 0.;
      for ( i = 1; i <= m; i ++ ) {
        sum += Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
      }
      ret = Math.max( ret, sum );
    }
  } else if ( norm.charAt(0).toUpperCase() == 'I' ) {
    for ( i = 1; i <= m; i ++ ) work[ ioffwork + i - 1 ] = 0.;
    for ( j = 1; j <= n; j ++ ) {
      for ( i = 1; i <= m; i ++ ) { 
        work[ ioffwork + i - 1 ] +=
          Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
      }
    }
    for ( i = 1; i <= m; i ++ ) {
      ret = Math.max( ret, work[ ioffwork + i - 1 ] );
    }
  } else if ( norm.charAt(0).toUpperCase() == 'F' ) {
    var scaleReference = new NumberReference( 0. );
    var sumrReference = new NumberReference( 1. );
    for ( j = 1; j <= n; j ++ ) {
      LaPack0.dlassq( m, A, 1, scaleReference, sumrReference,
        ioffa + ( j - 1 ) * lda );
    }
    ret = scaleReference.getValue()
      * Math.sqrt( sumrReference.getValue() );
  }
  return ret;
}
LaPack1.zlange = function( norm, m, n, A, lda, work, ioffa,
ioffwork ) {
  var ret = 0.;
  if ( Math.min( m, n ) == 0 ) return ret;
  if ( norm.charAt(0).toUpperCase() == 'M' ) {
    for ( var j = 1; j <= n; j ++ ) {
      for ( var i = 1; i <= m; i ++ ) {
        ret = Math.max( ret,
          ComplexMath.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] ) );
      }
    }
  } else if ( norm.charAt(0).toUpperCase() == 'O' ||
  norm.charAt(0) == '1' ) {
    for ( j = 1; j <= n; j ++ ) {
      var sum = 0.;
      for ( i = 1; i <= m; i ++ ) {
        sum +=
          ComplexMath.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
      }
      ret = Math.max( ret, sum );
    }
  } else if ( norm.charAt(0).toUpperCase() == 'I' ) {
    for ( i = 1; i <= m; i ++ ) work[ ioffwork + i - 1 ] = 0.;
    for ( j = 1; j <= n; j ++ ) {
      for ( i = 1; i <= m; i ++ ) { 
        work[ ioffwork + i - 1 ] +=
          ComplexMath.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
      }
    }
    for ( i = 1; i <= m; i ++ ) {
      ret = Math.max( ret, work[ ioffwork + i - 1 ] );
    }
  } else if ( norm.charAt(0).toUpperCase() == 'F' ) {
    var scaleReference = new NumberReference( 0. );
    var sumrReference = new NumberReference( 1. );
    for ( j = 1; j <= n; j ++ ) {
      LaPack0.zlassq( m, A, 1, scaleReference, sumrReference,
        ioffa + ( j - 1 ) * lda );
    }
    ret = scaleReference.getValue()
      * Math.sqrt( sumrReference.getValue() );
  }
  return ret;
}
//************************************************************************
LaPack1.dlangt = function( norm, n, dl, d, du, ioffdl, ioffd,
ioffdu ) {
  var anorm = Number.POSITIVE_INFINITY;
  if ( n <= 0 ) anorm = 0.;
  else if ( norm.charAt(0).toUpperCase() == 'M' ) {
    anorm = Math.abs( d[ ioffd + n - 1 ] );
    for ( var i = 1; i <= n - 1; i ++ ) {
      anorm = Math.max( anorm, Math.abs( dl[ ioffdl + i - 1 ] ) );
      anorm = Math.max( anorm, Math.abs( d[ ioffd + i - 1 ] ) );
      anorm = Math.max( anorm, Math.abs( du[ ioffdu + i - 1 ] ) );
    }
  } else if ( norm.charAt(0).toUpperCase() == 'O' ||
  norm.charAt(0).toUpperCase() == '1' ) {
    if ( n == 1 ) anorm = Math.abs( d[ ioffd ] );
    else {
      anorm = Math.max( Math.abs( d[ ioffd ] )
        + Math.abs( dl[ ioffdl ] ),
        Math.abs( d[ ioffd + n - 1 ] )
        + Math.abs( du[ ioffdu + n - 2 ] ) );
      for ( i = 2; i <= n - 1; i ++ ) {
        anorm = Math.max( anorm, Math.abs( d[ ioffd + i - 1 ] )
          + Math.abs( dl[ ioffdl + i - 1 ] )
          + Math.abs( du[ ioffdu + i - 2 ] ) );
      }
    }
  } else if ( norm.charAt(0).toUpperCase() == 'I' ) {
    if ( n == 1 ) anorm = Math.abs( d[ ioffd ] );
    else {
      anorm = Math.max( Math.abs( d[ ioffd ] )
        + Math.abs( du[ ioffdu ] ),
        Math.abs( d[ ioffd + n - 1 ] )
        + Math.abs( dl[ ioffdl + n - 2 ] ) );
      for ( i = 2; i <= n - 1; i ++ ) {
        anorm = Math.max( anorm, Math.abs( d[ ioffd + i - 1 ] )
          + Math.abs( du[ ioffdu + i - 1 ] )
          + Math.abs( dl[ ioffdl + i - 2 ] ) );
      }
    }
  } else if ( norm.charAt(0).toUpperCase() == 'F' ||
    norm.charAt(0).toUpperCase() == 'E' ) {
      var scaleReference = new NumberReference( 0. );
    var sumReference = new NumberReference( 1. );
    LaPack0.dlassq( n, d, 1, scaleReference, sumReference, ioffd );
    if ( n > 1 ) {
      LaPack0.dlassq( n - 1, dl, 1, scaleReference, sumReference,
        ioffdl );
      LaPack0.dlassq( n - 1, du, 1, scaleReference, sumReference,
        ioffdu );
    }
    anorm =
      scaleReference.getValue() * Math.sqrt( sumReference.getValue() );
  }
  return anorm;
}
LaPack1.zlangt = function( norm, n, dl, d, du ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack1.dlanhs = function( norm, n, A, lda, work, ioffa,
ioffwork ) {
  throw new Error("not programmed: hessenberg matrix");
}
LaPack1.zlanhs = function( norm, n, A, lda, work ) {
  throw new Error("not programmed: complex hessenberg matrix");
}
//************************************************************************
LaPack1.dlansb = function( norm, uplo, n, k, AB, ldab, work ) {
  throw new Error("not programmed: band matrix");
}
LaPack1.zlansb = function( norm, uplo, n, k, AB, ldab, work ) {
  throw new Error("not programmed: complex band matrix");
}
//************************************************************************
LaPack1.dlansf = function( norm, transr, uplo, n, A, work,
ioffa, ioffwork ) {
  throw new Error("not programmed: RFP format");
}
//************************************************************************
LaPack1.dlansp = function( norm, uplo, n, AP, work ) {
  throw new Error("not programmed: packed matrix");
}
LaPack1.zlansp = function( norm, uplo, n, AP, work ) {
  throw new Error("not programmed: complex packed matrix");
}
//************************************************************************
LaPack1.dlanst = function( norm, n, d, e, ioffd, ioffe) {
  var anorm = Number.POSITIVE_INFINITY;
  if ( n <= 0. ) anorm = 0.;
  else if ( norm.charAt(0).toUpperCase() == 'M' ) {
    anorm = Math.abs( d[ ioffd + n - 1 ] );
    for ( var i = 1; i <= n - 1; i ++ ) {
      anorm = Math.max( anorm, Math.abs( d[ ioffd + i - 1 ] ) );
      anorm = Math.max( anorm, Math.abs( e[ ioffe + i - 1 ] ) );
    }
  } else if ( norm.charAt(0).toUpperCase() == 'O' ||
  norm.charAt(0).toUpperCase() == '1' ||
  norm.charAt(0).toUpperCase() == 'I' ) {
    if ( n == 1 ) anorm = Math.abs( d[ ioffd ] );
    else {
      anorm = Math.max( Math.abs( d[ ioffd ] )
        + Math.abs( e[ ioffe ] ),
        Math.abs( e[ ioffe + n - 2 ] )
        + Math.abs( d[ ioffd + n - 1 ] ) );
      for ( i = 2; i <= n - 1; i ++ ) {
        anorm = Math.max( anorm, Math.abs( d[ ioffd + i - 1 ] )
          + Math.abs( e[ ioffe + i - 1 ] )
          + Math.abs( e[ ioffe + i - 2 ] ) );
      }
    }
  } else if ( norm.charAt(0).toUpperCase() == 'F' ||
  norm.charAt(0).toUpperCase() == 'E' ) {
    var scaleReference = new NumberReference( 0. );
    var sumReference = new NumberReference( 1. );
    if ( n > 1 ) {
      LaPack0.dlassq( n - 1, e, 1, scaleReference, sumReference, ioffe );
      sumReference.setValue( sumReference.getValue() * 2 );
    }
    LaPack0.dlassq( n, d, 1, scaleReference, sumReference, ioffd );
    anorm = scaleReference.getValue()
      * Math.sqrt( sumReference.getValue() );
  }
  return anorm;
}
//************************************************************************
LaPack1.dlansy = function( norm, uplo, n, A, lda, work, ioffa,
ioffwork ) {
  var value = 0.;
  if ( n == 0 ) return value;
  else if ( norm.charAt(0).toUpperCase() == 'M' ) {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      for ( var j = 1; j <= n; j ++ ) {
        for ( var i = 1; i <= j; i ++ ) {
          value = Math.max( value,
            Math.abs( A[ ioffa + i -1 + ( j - 1 ) * lda ] ) );
        }
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = j; i <= n; i ++ ) {
          value = Math.max( value,
            Math.abs( A[ ioffa + i -1 + ( j - 1 ) * lda ] ) );
        }
      }
    }
  } else if ( norm.charAt(0).toUpperCase() == 'I' ||
  norm.charAt(0).toUpperCase() == 'O' ||
  norm.charAt(0).toUpperCase() == '1' ) {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      for ( j = 1; j <= n; j ++ ) {
        var sum = 0.;
        for ( i = 1; i <= j - 1; i ++ ) {
          var absa =
            Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
          sum += absa;
          work[ ioffwork + i - 1 ] += absa;
        }
        work[ ioffwork + j - 1 ] =
          sum + Math.abs( A[ ioffa + j -1 + ( j - 1 ) * lda ] );
      }
      for ( i = 1; i <= n; i ++ ) {
        value = Math.max( value, work[ ioffwork + i - 1 ] );
      }
    } else {
      for ( i = 1; i <= n; i ++ ) {
        work[ ioffwork + i - 1 ] = 0.;
      }
      for ( j = 1; j <= n; j ++ ) {
        sum = work[ ioffwork + j - 1 ]
          + Math.abs( A[ ioffa + j - 1 + ( j - 1 ) * lda ] );
        for ( i = j + 1; i <= n; i ++ ) {
          absa = Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
          sum += absa;
          work[ ioffwork + i - 1 ] += absa;
        }
        value = Math.max( value, sum );
      }
    }
  } else if ( norm.charAt(0).toUpperCase() == 'F' ||
  norm.charAt(0).toUpperCase() == 'E' ) {
    var scaleReference = new NumberReference( 0. );
    var sumref = new NumberReference( 1. );
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      for ( j = 2; j <= n; j ++ ) {
        LaPack0.dlassq( j - 1, A, 1, scaleReference, sumref,
          ioffa + ( j - 1 ) * lda );
      }
    } else {
      for ( j = 1; j <= n - 1; j ++ ) {
        LaPack0.dlassq( n - j, A, 1, scaleReference, sumref,
          ioffa + j + ( j - 1 ) * lda );
      }
    }
    sumref.setValue( sumref.getValue() * 2. );
    LaPack0.dlassq( n, A, lda + 1, scaleReference, sumref, ioffa );
    value = scaleReference.getValue() * Math.sqrt( sumref.getValue() );
  }
  return value;
}
LaPack1.zlansy = function( norm, uplo, n, A, lda, work ) {
  throw new Error("not programmed: complex matrix");
  return 0.;
}
//************************************************************************
LaPack1.dlantb = function( norm, uplo, diag, n, k, AB, ldab,
work ) {
  throw new Error("not programmed: band matrix");
  return 0.;
}
LaPack1.zlantb = function( norm, uplo, diag, n, k, AB, ldab,
work ) {
  throw new Error("not programmed: complex band matrix");
  return 0.;
}
//************************************************************************
LaPack1.dlantp = function( norm, uplo, diag, n, AP, work ) {
  throw new Error("not programmed: packed matrix");
}
LaPack1.zlantp = function( norm, uplo, diag, n, AP, work ) {
  throw new Error("not programmed: complex packed matrix");
}
//************************************************************************
LaPack1.dlantr = function( norm, uplo, diag, m, n, A, lda, work,
ioffa, ioffwork) {
  var value = 0.;
  if ( Math.min( m, n ) == 0 ) return value;
  else if ( norm.charAt(0).toUpperCase() == 'M' ) {
    if ( diag.charAt(0).toUpperCase() == 'U' ) {
      value = 1.;
      if ( uplo.charAt(0).toUpperCase() == 'U' ) {
        for ( var j = 1; j <= n; j ++ ) {
          for ( var i = 1; i <= Math.min( m, j - 1 ); i ++ ) {
            value = Math.max( value,
              Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] ) );
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = j + 1; i <= m; i ++ ) {
            value = Math.max( value,
              Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] ) );
          }
        }
      }
    } else {
      if ( uplo.charAt(0).toUpperCase() == 'U' ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= Math.min( m, j ); i ++ ) {
            value = Math.max( value,
              Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] ) );
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = j; i <= m; i ++ ) {
            value = Math.max( value,
              Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] ) );
          }
        }
      }
    }
  } else if ( norm.charAt(0).toUpperCase() == 'O' ||
  norm.charAt(0).toUpperCase() == '1' ) {
    var udiag = ( diag.charAt(0).toUpperCase() == 'U' );
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      for ( j = 1; j <= n; j ++ ) {
        if ( udiag && j <= m ) {
          var sum = 1.;
          for ( i = 1; i <= j - 1; i ++ ) {
            sum += Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
          }
        } else {
          sum = 0.;
          for ( i = 1; i <= Math.min( m, j ); i ++ ) {
            sum += Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
          }
        }
        value = Math.max( value, sum );
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        if ( udiag ) {
          sum = 1.;
          for ( i = j + 1; i <= m; i ++ ) {
            sum += Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
          }
        } else {
          sum = 0.;
          for ( i = j; i <= m; i ++ ) {
            sum += Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
          }
        }
        value = Math.max( value, sum );
      }
    }
  } else if ( norm.charAt(0).toUpperCase() == 'I' ) {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      if ( diag.charAt(0).toUpperCase() == 'U' ) {
        for ( i = 1; i <= m; i ++ ) work[ ioffwork + i - 1 ] = 1.;
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= Math.min( m, j - 1 ); i ++ ) {
            work[ ioffwork + i - 1 ] +=
              Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
          }
        }
      } else {
        for ( i = 1; i <= m; i ++ ) work[ ioffwork + i - 1 ] = 0.;
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= Math.min( m, j ); i ++ ) {
            work[ ioffwork + i - 1 ] +=
              Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
          }
        }
      }
    } else {
      if ( diag.charAt(0).toUpperCase() == 'U' ) {
        for ( i = 1; i <= n; i ++ ) work[ ioffwork + i - 1 ] = 1.;
        for ( i = n + 1; i <= m; i ++ ) work[ ioffwork + i - 1 ] = 0.;
        for ( j = 1; j <= n; j ++ ) {
          for ( i = j + 1; i <= m; i ++ ) {
            work[ ioffwork + i - 1 ] +=
              Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
          }
        }
      } else {
        for ( i = 1; i <= m; i ++ ) work[ ioffwork + i - 1 ] = 0.;
        for ( j = 1; j <= n; j ++ ) {
          for ( i = j; i <= m; i ++ ) {
            work[ ioffwork + i - 1 ] +=
              Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
          }
        }
      }
    }
    for ( i = 1; i <= m; i ++ ) {
      value = Math.max( value, work[ ioffwork + i - 1 ] );
    }
  } else if ( norm.charAt(0).toUpperCase() == 'F' ||
  norm.charAt(0).toUpperCase() == 'E' ) {
    var scaleReference = new NumberReference();
    var sumref = new NumberReference();
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      if ( diag.charAt(0).toUpperCase() == 'U' ) {
        scaleReference.setValue( 1. );
        sumref.setValue( Math.min( m, n ) );
        for ( j = 2; j <= n; j ++ ) {
          LaPack0.dlassq( Math.min( m, j - 1 ), A, 1, scaleReference,
            sumref, ioffa + ( j - 1 ) * lda );
        }
      } else {
        scaleReference.setValue( 0. );
        sumref.setValue( 1. );
        for ( j = 1; j <= n; j ++ ) {
          LaPack0.dlassq( Math.min( m, j ), A, 1, scaleReference, sumref,
            ioffa + ( j - 1 ) * lda );
        }
      }
    } else {
      if ( diag.charAt(0).toUpperCase() == 'U' ) {
        scaleReference.setValue( 1. );
        sumref.setValue( Math.min( m, n ) );
        for ( j = 1; j <= n; j ++ ) {
          LaPack0.dlassq( m - j, A, 1, scaleReference, sumref,
            ioffa + Math.min( m - 1, j ) + ( j - 1 ) * lda );
        }
      } else {
        scaleReference.setValue( 0. );
        sumref.setValue( 1. );
        for ( j = 1; j <= n; j ++ ) {
          LaPack0.dlassq( m - j + 1, A, 1, scaleReference, sumref,
            ioffa + j - 1 + ( j - 1 ) * lda );
        }
      }
    }
    value = scaleReference.getValue() * Math.sqrt( sumref.getValue() );
  }
  return value;
}
LaPack1.zlantr = function( norm, uplo, diag, m, n, A, lda,
work ) {
  throw new Error("not programmed: complex matrix");
  return 0.;
}
//************************************************************************
LaPack1.dlanv2 = function( aReference, bReference, cReference,
dReference, rt1rReference, rt1iReference, rt2rReference, rt2iReference,
csReference, snReference ) {
  var multpl = 4.;
  var eps = LaPack0.dlamch( 'P' );
  if ( cReference.getValue() == 0. ) {
    csReference.setValue( 1. );
    snReference.setValue( 0. );
  } else if ( bReference.getValue() == 0. ) {
    csReference.setValue( 0. );
    snReference.setValue( 1. );
    var temp = dReference.getValue();
    dReference.setValue( aReference.getValue() );
    aReference.setValue( temp );
    bReference.setValue( -cReference.getValue() );
    cReference.setValue( 0. );
  } else if ( ( aReference.getValue() - dReference.getValue() ) == 0. &&
  ( bReference.getValue() >= 0. != cReference.getValue() >= 0. ) ) {
    csReference.setValue( 1. );
    snReference.setValue( 0. );
  } else {
    temp = aReference.getValue() - dReference.getValue();
    var p = 0.5 * temp;
    var bcmax =
      Math.max( Math.abs( bReference.getValue() ),
                Math.abs( cReference.getValue() ) );
    var bcmis = Math.min( Math.abs( bReference.getValue() ),
      Math.abs( cReference.getValue() ) ) *
      ( bReference.getValue() >= 0. ? 1. : -1. )
      * ( cReference.getValue() >= 0. ? 1. : -1. );
    var scale = Math.max( Math.abs( p ), bcmax );
    var z = ( p / scale ) * p + ( bcmax / scale ) * bcmis;
    if ( z >= multpl * eps ) {
      z = p + Math.sqrt( scale ) * Math.sqrt ( z ) 
        * ( p >= 0. ? 1. : -1. );
      aReference.setValue( dReference.getValue() + z );
      dReference.setValue( dReference.getValue()
        - ( bcmax / z ) * bcmis );
      var tau = LaPack0.dlapy2( cReference.getValue(), z );
      csReference.setValue( z / tau );
      snReference.setValue( cReference.getValue() / tau );
      bReference.setValue( bReference.getValue()
        - cReference.getValue() );
      cReference.setValue( 0. );
    } else {
      var sigma = bReference.getValue() + cReference.getValue();
      tau = LaPack0.dlapy2( sigma, temp );
      csReference.setValue(
        Math.sqrt( 0.5 * ( 1. + Math.abs( sigma ) / tau ) ) );
      snReference.setValue( - ( p / ( tau * csReference.getValue() ) )
        * ( sigma >= 0. ? 1. : -1. ) );
      var aa = aReference.getValue() * csReference.getValue()
        + bReference.getValue() * snReference.getValue();
      var bb = - aReference.getValue() * snReference.getValue()
        + bReference.getValue() * csReference.getValue();
      var cc = cReference.getValue() * csReference.getValue()
        + dReference.getValue() * snReference.getValue();
      var dd = - cReference.getValue() * snReference.getValue()
        + dReference.getValue() * csReference.getValue();
      aReference.setValue( aa * csReference.getValue()
        + cc * snReference.getValue() );
      bReference.setValue( bb * csReference.getValue()
        + dd * snReference.getValue() );
      cReference.setValue( - aa * snReference.getValue()
        + cc * csReference.getValue() );
      dReference.setValue( - bb * snReference.getValue()
        + dd * csReference.getValue() );
      temp = 0.5 * ( aReference.getValue() + dReference.getValue() );
      aReference.setValue( temp );
      dReference.setValue( temp );
      if ( cReference.getValue() != 0. ) {
        if ( bReference.getValue() != 0. ) {
          if ( ( bReference.getValue() >= 0. ) ==
          ( cReference.getValue() >= 0. ) ) {
            var sab = Math.sqrt( Math.abs( bReference.getValue() ) );
            var sac = Math.sqrt( Math.abs( cReference.getValue() ) );
            p = sab * sac * ( cReference.getValue() >= 0. ? 1. : - 1. );
            tau = 1. / Math.sqrt( Math.abs( bReference.getValue()
              + cReference.getValue() ) );
            aReference.setValue( temp + p );
            dReference.setValue( temp - p );
            bReference.setValue( bReference.getValue()
              - cReference.getValue() );
            cReference.setValue( 0. );
            var cs1 = sab * tau;
            var sn1 = sac * tau;
            temp = csReference.getValue() * cs1
              - snReference.getValue() * sn1;
            snReference.setValue( csReference.getValue() * sn1
              + snReference.getValue() * cs1 );
            csReference.setValue( temp );
          }
        } else {
          bReference.setValue( -cReference.getValue() );
          cReference.setValue( 0. );
          temp = csReference.getValue();
          csReference.setValue( - snReference.getValue() );
          snReference.setValue( temp );
        }
      }
    }
  }
  rt1rReference.setValue( aReference.getValue() );
  rt2rReference.setValue( dReference.getValue() );
  if ( cReference.getValue() == 0. ) {
    rt1iReference.setValue( 0. );
    rt2iReference.setValue( 0. );
  } else {
    rt1iReference.setValue( Math.sqrt( Math.abs( bReference.getValue() ) )
      * Math.sqrt( Math.abs( cReference.getValue() ) ) );
    rt2iReference.setValue( - rt1iReference.getValue() );
  }
}
//************************************************************************
LaPack1.dlaqgb = function( m, n, kl, ku, AB, ldab, r, c,
rowcndReference, colcndReference, amax, equedReference ) {
  throw new Error("not programmed: band matrix");
}
LaPack1.zlaqgb = function( m, n, kl, ku, AB, ldab, r, c,
rowcndReference, colcndReference, amax, equedReference ) {
  throw new Error("not programmed: complex band matrix");
}
//************************************************************************
LaPack1.dlaqge = function( m, n, A, lda, r, c, rowcnd, colcnd,
amax, equedReference, ioffa, ioffr, ioffc ) {
  var thresh = 0.1;
  if ( m <= 0 || n <= 0 ) {
    equedReference.setValue( 'N' );
    return;
  }
  var small = LaPack0.dlamch( 'Safe minimum' )
    / LaPack0.dlamch( 'Precision' );
  var large = 1. / small;
  if ( rowcnd >= thresh && amax >= small && amax <= large ) {
    if ( colcnd >= thresh ) equedReference.setValue( 'N' );
    else {
      for ( var j = 1; j <= n; j ++ ) {
        var cj = c[ ioffc + j - 1 ];
        for ( var i = 1; i <= m; i ++ ) {
          A[ ioffa + i - 1 + ( j - 1 ) * lda ] *= cj;
        }
      }
      equedReference.setValue( 'C' );
    }
  } else if ( colcnd >= thresh ) {
    for ( j = 1; j <= n; j ++ ) {
      for ( i = 1; i <= m; i ++ ) {
        A[ ioffa + i - 1 + ( j - 1 ) * lda ] *= r[ ioffr + i - 1 ];
      }
    }
    equedReference.setValue( 'R' );
  } else {
    for ( j = 1; j <= n; j ++ ) {
      cj = c[ ioffc + j - 1 ];
      for ( i = 1; i <= m; i ++ ) {
        A[ ioffa + i - 1 + ( j - 1 ) * lda ] *=
          cj * r[ ioffr + i - 1 ];
      }
    }
    equedReference.setValue( 'B' );
  }
}
LaPack1.zlaqge = function( m, n, A, lda, r, c, rowcnd, colcnd,
amax, equedReference ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack1.dlaqsb = function( uplo, n, kd, AB, ldab, s, scond,
amax, equedReference) {
  throw new Error("not programmed: band matrix");
}
LaPack1.zlaqsb = function( uplo, n, kd, AB, ldab, s, scond,
amax, equedReference) {
  throw new Error("not programmed: band matrix");
}
//************************************************************************
LaPack1.dlaqsp = function( uplo, n, AP, s, scond, amax,
equed ) {
  throw new Error("not programmed: packed matrix");
}
LaPack1.zlaqsp = function( uplo, n, AP, s, scond, amax,
equed ) {
  throw new Error("not programmed: complex packed matrix");
}
//************************************************************************
LaPack1.dlaqsy = function( uplo, n, A, lda, s, scond, amax,
equedReference, ioffa, ioffs ) {
  var thresh = 0.1;
  if ( n <= 0 ) {
    equedReference.setValue( 'N' );
    return;
  }
  var small = LaPack0.dlamch( 'Safe minimum' )
    / LaPack0.dlamch( 'Precision' );
  var large = 1. / small;
  if ( scond >= thresh && amax >= small && amax <= large ) {
    equedReference.setValue( 'N' );
  } else {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      for ( var j = 1; j <= n; j ++ ) {
        var cj = s[ ioffs + j - 1 ];
        for ( var i = 1; i <= j; i ++ ) {
          A[ ioffa + i - 1 + ( j - 1 ) * lda ] *=
            cj * s[ ioffs + i - 1 ];
        }
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        cj = s[ ioffs + j - 1 ];
        for ( i = j; i <= n; i ++ ) {
          A[ ioffa + i - 1 + ( j - 1 ) * lda ] *=
            cj * s[ ioffs + i - 1 ];
        }
      }
    }
    equedReference.setValue( 'Y' );
  }
}
LaPack1.zlaqsy = function( uplo, n, A, lda, s, scond, amax,
equedReference ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack1.dlar1v = function( n, b1, bn, lambda, d, l, ld, lld,
pivmin, gaptol, z, wantnc, negcnt, ztzReference, mingmaReference, r,
isuppz, nrminvReference, residReference, rqcorrReference, work, ioffd,
ioffl, ioffld, iofflld, ioffz, ioffisuppz, ioffwork ) {
  throw new Error("not tested: complicated input");
  var eps = LaPack0.dlamch( 'Precision' );
  if ( r.getValue() == 0 ) {
    var r1 = b1;
    var r2 = bn;
  } else {
    r1 = r.getValue();
    r2 = r.getValue();
  }
  var indlpl = 0;
  var indumn = n
  var inds = 2 * n + 1;
  var indp = 3 * n + 1;
  if ( b1 == 1 ) work[ ioffwork + inds - 1 ] = 0.;
  else work[ ioffwork + inds + b1 - 2 ] = lld[ iofflld + b1 - 2 ];
  var sawnan1 = false;
  var neg1 = 0;
  var s = work[ ioffwork + inds + b1 - 2 ] - lambda;
  for ( var i = b1; i <= r1 - 1; i ++ ) {
    var dplus = d[ ioffd + i - 1 ] + s;
    work[ ioffwork + indlpl + i - 1 ] = ld[ ioffld + i - 1 ] / dplus;
    if ( dplus < 0. ) neg1 ++;
    work[ ioffwork + inds + i - 1 ] =
      s * work[ ioffwork + indlpl + i - 1 ] * l[ ioffl + i - 1 ];
    s = work[ ioffwork + inds + i - 1 ] - lambda;
  }
  sawnan1 = isNaN( s );
  if ( ! sawnan1 ) {
    for ( i = r1; i <= r2 - 1; i ++ ) {
      dplus = d[ ioffd + i - 1 ] + s;
      work[ ioffwork + indlpl + i - 1 ] = ld[ ioffld + i - 1 ] / dplus;
      work[ ioffwork + inds + i - 1 ] =
        s * work[ ioffwork + indlpl + i - 1 ] * l[ ioffl + i - 1 ];
      s = work[ ioffwork + inds + i - 1 ] - lambda;
    }
    sawnan1 = isNaN( s );
  }
  if ( sawnan1 ) {
    neg1 = 0;
    s = work[ ioffwork + inds + b1 - 2 ] - lambda;
    for ( i = b1; i <= r1 - 1; i ++ ) {
      dplus = d[ ioffd + i - 1 ] + s;
      if ( Math.abs( dplus ) < pivmin ) dplus = - pivmin;
      work[ ioffwork + indlpl + i - 1 ] = ld[ ioffld + i - 1 ] / dplus;
      if ( dplus < 0. ) neg1 ++;
      work[ ioffwork + inds + i - 1 ] =
        s * work[ ioffwork + indlpl + i - 1 ] * l[ ioffl + i - 1 ];
      if ( work[ ioffwork + indlpl + i - 1 ] == 0. ) {
        work[ ioffwork + inds + i - 1 ] = lld[ iofflld + i - 1 ];
      }
      s = work[ ioffwork + inds + i - 1 ] - lambda;
    }
    for ( i = r1; i <= r2 - 1; i ++ ) {
      dplus = d[ ioffd + i - 1 ] + s;
      if ( Math.abs( dplus ) < pivmin ) dplus = - pivmin;
      work[ ioffwork + indlpl + i - 1 ] = ld[ ioffld + i - 1 ] / dplus;
      work[ ioffwork + inds + i - 1 ] =
        s * work[ indlpl + i - 1 ] * l[ ioffl + i - 1 ];
      if ( work[ ioffwork + indlpl + i - 1 ] == 0. ) {
        work[ ioffwork + inds + i - 1 ] = lld[ iofflld + i - 1 ];
      }
      s = work[ ioffwork + inds + i - 1 ] - lambda;
    }
  }
  var sawnan2 = false;
  var neg2 = 0;
  work[ ioffwork + indp + bn - 2 ] = d[ ioffd + bn - 1 ] - lambda;
  for ( i = bn - 1; i >= r1; i -- ) {
    var dminus =
      lld[ iofflld + i - 1 ] + work[ ioffwork + indp + i - 1 ];
    var tmp = d[ ioffd + i - 1 ] / dminus;
    if ( dminus < 0. ) neg2 ++;
    work[ ioffwork + indumn + i - 1 ] = l[ ioffl + i - 1 ] * tmp;
    work[ ioffwork + indp + i - 2 ] =
      work[ ioffwork + indp + i - 1 ] * tmp - lambda;
  }
  tmp = work[ ioffwork + indp + r1 - 2 ];
  sawnan2 = isNaN( tmp );
  if ( sawnan2 ) {
    neg2 = 0;
    for ( i = bn - 1; i >= r1; i -- ) {
      dminus =
        lld[ iofflld + i - 1 ] + work[ ioffwork + indp + i - 1 ];
      if ( Math.abs( dminus ) < pivmin ) dminus = - pivmin;
      tmp = d[ ioffd + i - 1 ] / dminus;
      if ( dminus < 0. ) neg2 ++;
      work[ ioffwork + indumn + i - 1 ] = l[ ioffl + i - 1 ] * tmp;
      work[ ioffwork + indp + i - 2 ] =
        work[ ioffwork + indp + i - 1 ] * tmp - lambda;
      if ( tmp == 0. ) {
        work[ ioffwork + indp + i - 2 ] = d[ ioffd + i - 1 ] - lambda;
      }
    }
  }
  mingma.setValue( work[ ioffwork + inds + r1 - 2 ]
    + work[ ioffwork + indp + r1 - 2 ] );
  if ( mingma.getValue() < 0. ) neg1 ++;
  negcnt.setValue( ( wantnc ? neg1 + neg2 : -1 ) );
  if ( Math.abs( mingma.getValue() ) == 0. ) {
    mingma.setValue( eps * work[ ioffwork + inds + r1 - 2 ] );
  }
  r.setValue( r1 );
  for ( i = r1; i <= r2 - 1; i ++ ) {
    tmp = work[ ioffwork + inds + i - 1 ]
      + work[ ioffwork + indp + i - 1 ];
    if ( tmp == 0. ) tmp = eps * work[ ioffwork + inds + i - 1 ];
    if ( Math.abs( tmp ) <= Math.abs( mingma.getValue() ) ) {
      mingma.setValue( tmp );
      r.setValue( i + 1 );
    }
  }
  isuppz[ ioffisuppz ] = b1;
  isuppz[ ioffisuppz + 1 ] = bn;
  z[ ioffz + r.getValue() - 1 ] = 1.;
  ztz.setValue( 1. );
  if ( ! sawnan1 && ! sawnan2 ) {
    for ( i = r.getValue() - 1; i >= b1; i -- ) {
      z[ ioffz + i - 1 ] =
        - ( work[ ioffwork + indlpl + i - 1 ] * z[ ioffz + i ] );
      if ( ( Math.abs( z[ ioffz + i - 1 ] )
      + Math.abs( z[ ioffz + i ] ) ) * Math.abs( ld[ ioffld + i - 1 ] )
      < gaptol ) {
        z[ ioffz + i - 1 ] = 0.;
        isuppz[ ioffisuppz ] = i + 1;
        break;
      }
      ztz.setValue( ztz.getValue() + z[ ioffz + i - 1 ] * z[ ioffz + i - 1 ] );
    }
  } else {
    for ( i = r.getValue() - 1; i >= b1; i -- ) {
      if ( z[ ioffz + i ] == 0. ) {
        z[ ioffz + i - 1 ] = - ( ld[ ioffld + i ]
          / ld[ ioffld + i - 1 ] ) * ( z[ ioffz + i + 1 ] );
      } else {
        z[ ioffz + i - 1 ] = - ( work[ ioffwork+ indlpl + i - 1 ]
          * z[ ioffz + i ] );
      }
      if ( Math.abs( z[ ioffz + i - 1 ] + Math.abs( z[ ioffz + i ] ) )
      * Math.abs( ld[ ioffld + i - 1 ] ) < gaptol ) {
        z[ ioffz + i - 1 ] = 0.;
        isuppz[ ioffisuppz ] = i + 1;
        break;
      }
      ztz.setValue( ztz.getValue() + z[ ioffz + i - 1 ] * z[ ioffz + i - 1 ] );
    }
  }
  if ( ! sawnan1 && ! sawnan2 ) {
    for ( i = r.getValue(); i <= bn - 1; i ++ ) {
      z[ ioffz + i ] = - ( work[ ioffwork + indumn + i - 1 ]
        * z[ ioffz + i - 1 ] );
      if ( ( Math.abs( z[ ioffz + i - 1 ] )
      + Math.abs( z[ ioffz + i ] ) ) * Math.abs( ld[ ioffld + i - 1 ] )
      < gaptol ) {
        z[ ioffz + i ] = 0.;
        isuppz[ ioffisuppz + 1 ] = i;
        break;
      }
      ztz.setValue( ztz.getValue() + z[ ioffz + i ] * z[ ioffz + i ] );
    }
  } else {
    for ( i = r.getValue(); i <= bn - 1; i ++ ) {
      if ( z[ ioffz + i - 1 ] == 0. ) {
        z[ ioffz + i ] = - ( ld[ ioffld + i - 2 ]
          / ld[ ioffld + i - 1 ] ) * z[ ioffz + i - 2 ];
      } else {
        z[ ioffz + i ] = - ( work[ ioffwork + indumn + i - 1 ]
          + z[ ioffz + i - 1 ] );
      }
      if ( Math.abs( z[ ioffz + i - 1 ] )
      + Math.abs( z[ ioffz + i ] ) * Math.abs( ld[ ioffld + i - 1 ] )
      < gaptol ) {
        z[ ioffz + i ] = 0.;
        isuppz[ ioffisuppz + 1 ] = i;
        break;
      }
      ztz.setValue( ztz.getValue() + z[ ioffz + i ] * z[ ioffz + i ] );
    }
  }
  tmp = 1. / ztz.getValue();
  nrminv.setValue( Math.sqrt( tmp ) );
  resid.setValue( Math.abs( mingma.getValue() ) * nrminv.getValue() );
  rqcorr.setValue( mingma.getValue() * tmp );
}
//************************************************************************
LaPack1.dlarf = function( side, m, n, v, incv, tau, C, ldc,
work, ioffv, ioffc, ioffwork ) {
  var applyleft = ( side.charAt(0).toUpperCase() == 'L' );
  var lastv = 0;
  var lastc = 0;
  if ( tau != 0. ) {
    lastv = ( applyleft ? m : n );
    var i = ( incv > 0 ? 1 + ( lastv - 1 ) * incv : 1 );
    while ( lastv > 0 && v[ ioffv + i - 1 ] == 0. ) {
      lastv --;
      i -= incv;
    }
    lastc = ( applyleft ? LaPack0.iladlc( lastv, n, C, ldc, ioffc )
      : LaPack0.iladlr( m, lastv, C, ldc, ioffc ) );
  }
  if ( applyleft ) {
    if ( lastv > 0 ) {
      Blas2.dgemv( 'Transpose', lastv, lastc, 1., C, ldc, v, incv, 0.,
        work, 1, ioffc, ioffv, ioffwork );
      Blas2.dger( lastv, lastc, -tau, v, incv, work, 1, C, ldc, ioffv,
        ioffwork, ioffc );
    }
  } else {
    if ( lastv > 0 ) {
      Blas2.dgemv( 'No transpose', lastc, lastv, 1., C, ldc, v, incv,
        0., work, 1, ioffc, ioffv, ioffwork );
      Blas2.dger( lastc, lastv, -tau, work, 1, v, incv, C, ldc,
        ioffwork, ioffv, ioffc );
    }
  }
}
LaPack1.zlarf = function( side, m, n, v, incv, tau, C,
ldc, work, ioffv, ioffc, ioffwork ) {
  var one = new Complex( 1., 0. );
  var zero = new Complex( 0., 0. );
  if ( side.charAt(0).toUpperCase() == 'L' ) {
    if ( ! tau.equals( zero ) ) {
      Blas2.zgemv( 'C', m, n, one, C, ldc, v, incv, zero, work, 1,
        ioffc, ioffv, ioffwork );
      Blas2.zgerc( m, n, tau.minus(), v, incv, work, 1, C, ldc,
        ioffv, ioffwork, ioffc );
    }
  } else {
    if ( ! tau.equals( zero ) ) {
      Blas2.zgemv( 'N', m, n, one, C, ldc, v, incv, zero, work, 1,
        ioffc, ioffv, ioffwork );
      Blas2.zgerc( m, n, tau.minus(), work, 1, v, incv, C, ldc,
        ioffwork, ioffv, ioffc );
    }
  }
}
//***********************************************************************
LaPack1.dlarfb = function( side, trans, direct, storev, m, n,
k, V, ldv, T, ldt, C, ldc, Work, ldwork, ioffv, iofft, ioffc, ioffwork ) {
  if ( m <= 0 || n <= 0 ) return;
  var transt =
    ( trans.charAt(0).toUpperCase() == 'N' ? 'T' : 'N' );
  if ( storev.charAt(0).toUpperCase() == 'C' ) {
    if ( direct.charAt(0).toUpperCase() == 'F' ) {
      if ( side.charAt(0).toUpperCase() == 'L' ) {
//          V=[V1] kxk unit lower
//            [V2] (m-k)xk
//          C=[C1] kxn
//            [C2] (m-k)xn
        var lastv =
          Math.max( k, LaPack0.iladlr( m, k, V, ldv, ioffv ) );
        var lastc = LaPack0.iladlc( lastv, n, C, ldc, ioffc );
        for ( var j = 1; j <= k; j ++ ) {
          Blas1.dcopy( lastc, C, ldc, Work, 1, ioffc + j - 1,
            ioffwork + ( j - 1 ) * ldwork );
        }
        Blas3.dtrmm( 'Right', 'Lower', 'No transpose', 'Unit', lastc,
          k, 1., V, ldv, Work, ldwork, ioffv, ioffwork );
        if ( lastv > k ) {
          Blas3.dgemm( 'Transpose', 'No transpose', lastc, k,
            lastv - k, 1., C, ldc, V, ldv, 1., Work, ldwork, ioffc + k,
            ioffv + k, ioffwork );
        }
        Blas3.dtrmm( 'Right', 'Upper', transt, 'Non-unit', lastc, k,
          1., T, ldt, Work, ldwork, iofft, ioffwork );
        if ( lastv > k ) {
          Blas3.dgemm( 'No transpose', 'Transpose', lastv - k, lastc,
            k, -1., V, ldv, Work, ldwork, 1., C, ldc, ioffv + k,
            ioffwork, ioffc + k );
        }
        Blas3.dtrmm( 'Right', 'Lower', 'Transpose', 'Unit', lastc, k,
          1., V, ldv, Work, ldwork, ioffv, ioffwork );
        for ( j = 1; j <= k; j ++ ) {
          for ( var i = 1; i <= lastc; i ++ ) {
            C[ ioffc + j - 1 + ( i - 1 ) * ldc ] -=
              Work[ ioffwork + i - 1 + ( j - 1 ) * ldwork ];
          }
        }
      } else if ( side.charAt(0).toUpperCase() == 'R' ) {
//          V=[V1] kxk unit lower
//            [V2] (n-k)xk
//          C=[C1,C2] mxk, mx(n-k)
        lastv = Math.max( k, LaPack0.iladlr( n, k, V, ldv, ioffv ) );
        lastc = LaPack0.iladlr( m, lastv, C, ldc, ioffc );
        for ( j = 1; j <= k; j ++ ) {
          Blas1.dcopy( lastc, C, 1, Work, 1, ioffc + ( j - 1 ) * ldc,
            ioffwork + ( j - 1 ) * ldwork );
        }
        Blas3.dtrmm( 'Right', 'Lower', 'No transpose', 'Unit', lastc,
          k, 1., V, ldv, Work, ldwork, ioffv, ioffwork );
        if ( lastv > k ) {
          Blas3.dgemm( 'No transpose', 'No transpose', lastc, k,
            lastv - k, 1., C, ldc, V, ldv, 1., Work, ldwork,
            ioffc + k * ldc, ioffv + k, ioffwork );
        }
        Blas3.dtrmm( 'Right', 'Upper', trans, 'Non-unit', lastc, k, 1.,
          T, ldt, Work, ldwork, iofft, ioffwork );
        if ( lastv > k ) {
          Blas3.dgemm( 'No transpose', 'Transpose', lastc, lastv - k,
            k, -1., Work, ldwork, V, ldv, 1., C, ldc, ioffwork,
            ioffv + k, ioffc + k * ldc );
        }
        Blas3.dtrmm( 'Right', 'Lower', 'Transpose', 'Unit', lastc, k,
          1., V, ldv, Work, ldwork, ioffv, ioffwork );
        for ( j = 1; j <= k; j ++ ) {
          for ( i = 1; i <= lastc; i ++ ) {
            C[ ioffc + i - 1 + ( j - 1 ) * ldc ] -=
              Work[ ioffwork + i - 1 + ( j - 1 ) * ldwork ];
          }
        }
      }
    } else {
      if ( side.charAt(0).toUpperCase() == 'L' ) {
//          V=[V1] (m-k)xk
//            [V2] kxk unit upper
//          C=[C1] (m-k)xn
//            [C2] kxn
        lastv = Math.max( k, LaPack0.iladlr( m, k, V, ldv, ioffv ) );
        lastc = LaPack0.iladlc( lastv, n, C, ldc, ioffc );
        for ( j = 1; j <= k; j ++ ) {
          Blas1.dcopy( lastc, C, ldc, Work, 1,
            ioffc + lastv - k + j - 1, ioffwork + ( j - 1 ) * ldwork );
        }
        Blas3.dtrmm( 'Right', 'Upper', 'No transpose', 'Unit', lastc,
          k, 1., V, ldv, Work, ldwork, ioffv + lastv - k, ioffwork );
        if ( lastv > k ) {
          Blas3.dgemm( 'Transpose', 'No transpose', lastc, k,
            lastv - k, 1., C, ldc, V, ldv, 1., Work, ldwork, ioffc,
            ioffv, ioffwork );
        }
        Blas3.dtrmm( 'Right', 'Lower', transt, 'Non-unit', lastc, k,
          1., T, ldt, Work, ldwork, iofft, ioffwork );
        if ( lastv > k ) {
          Blas3.dgemm( 'No transpose', 'Transpose', lastv - k, lastc,
            k, -1., V, ldv, Work, ldwork, 1., C, ldc, ioffv, ioffwork,
            ioffc );
        }
        Blas3.dtrmm( 'Right', 'Upper', 'Transpose', 'Unit', lastc, k,
          1., V, ldv, Work, ldwork, ioffv + lastv - k, ioffwork );
        for ( j = 1; j <= k; j ++ ) {
          for ( i = 1; i <= lastc; i ++ ) {
            C[ ioffc + lastv - k + j - 1 + ( i - 1 ) * ldc ] -=
              Work[ ioffwork + i - 1 + ( j - 1 ) * ldwork ];
          }
        }
      } else if ( side.charAt(0).toUpperCase() == 'R' ) {
//          V=[V1] (n-k)xk
//            [V2] kxk unit upper
//          C=[C1,C2] mx(n-k),mxk
        lastv = Math.max( k, LaPack0.iladlr( n, k, V, ldv, ioffv ) );
        lastc = LaPack0.iladlr( m, lastv, C, ldc, ioffc );
        for ( j = 1; j <= k; j ++ ) {
          Blas1.dcopy( lastc, C, 1, Work, 1,
            ioffc + ( n - k + j - 1 ) * ldc,
            ioffwork + ( j - 1 ) * ldwork );
        }
        Blas3.dtrmm( 'Right', 'Upper', 'No transpose', 'Unit', lastc,
          k, 1., V, ldv, Work, ldwork, ioffv + lastv - k, ioffwork );
        if ( lastv > k ) {
          Blas3.dgemm( 'No transpose', 'No transpose', lastc, k, 
            lastv - k, 1., C, ldc, V, ldv, 1., Work, ldwork, ioffc,
            ioffv, ioffwork );
        }
        Blas3.dtrmm( 'Right', 'Lower', trans, 'Non-unit', lastc, k, 1.,
          T, ldt, Work, ldwork, iofft, ioffwork );
        if ( lastv > k ) {
          Blas3.dgemm( 'No transpose', 'Transpose', lastc, lastv - k,
            k, -1., Work, ldwork, V, ldv, 1., C, ldc, ioffwork, ioffv,
            ioffc );
        }
        Blas3.dtrmm( 'Right', 'Upper', 'Transpose', 'Unit', lastc, k,
          1., V, ldv, Work, ldwork, ioffv + lastv - k, ioffwork );
        for ( j = 1; j <= k; j ++ ) {
          for ( i = 1; i <= lastc; i ++ ) {
            C[ ioffc + i - 1 + ( lastv - k + j - 1 ) * ldc ] -=
              Work[ ioffwork + i - 1 + ( j - 1 ) * ldwork ];
          }
        }
      }
    }
  } else if ( storev.charAt(0).toUpperCase() == 'R' ) {
    if ( direct.charAt(0).toUpperCase() == 'F' ) {
      if ( side.charAt(0).toUpperCase() == 'L' ) {
//          V=[V1,V2] kxk unit upper, kx(m-k)
//          C=[C1] kxn
//            [C2] (m-k)xn
        lastv = Math.max( k, LaPack0.iladlc( k, m, V, ldv, ioffv ) );
        lastc = LaPack0.iladlc( lastv, n, C, ldc, ioffc );
        for ( j = 1; j <= k; j ++ ) {
          Blas1.dcopy( lastc, C, ldc, Work, 1, ioffc + j - 1,
            ioffwork + ( j - 1 ) * ldwork );
        }
        Blas3.dtrmm( 'Right', 'Upper', 'Transpose', 'Unit', lastc, k,
          1., V, ldv, Work, ldwork, ioffv, ioffwork );
        if ( lastv > k ) {
          Blas3.dgemm( 'Transpose', 'Transpose', lastc, k,
            lastv - k, 1., C, ldc, V, ldv, 1., Work, ldwork, ioffc + k,
            ioffv + k * ldv, ioffwork );
        }
        Blas3.dtrmm( 'Right', 'Upper', transt, 'Non-unit', lastc, k,
          1., T, ldt, Work, ldwork, iofft, ioffwork );
        if ( lastv > k ) {
          Blas3.dgemm( 'Transpose', 'Transpose', lastv - k, lastc, k,
            -1., V, ldv, Work, ldwork, 1., C, ldc, ioffv + k * ldv,
            ioffwork, ioffc + k );
        }
        Blas3.dtrmm( 'Right', 'Upper', 'No transpose', 'Unit', lastc,
          k, 1., V, ldv, Work, ldwork, ioffv, ioffwork );
        for ( j = 1; j <= k; j ++ ) {
          for ( i = 1; i <= lastc; i ++ ) {
            C[ ioffc + j - 1 + ( i - 1 ) * ldc ] -=
              Work[ ioffwork + i - 1 + ( j - 1 ) * ldwork ];
          }
        }
      } else if ( side.charAt(0).toUpperCase() == 'R' ) {
//          V=[V1,V2] kxk unit upper, kx(n-k)
//          C=[C1,C2] mxk, mx(n-k)
        lastv = Math.max( k, LaPack0.iladlc( k, n, V, ldv, ioffv ) );
        lastc = LaPack0.iladlr( m, lastv, C, ldc, ioffc );
        for ( j = 1; j <= k; j ++ ) {
          Blas1.dcopy( lastc, C, 1, Work, 1, ioffc + (j - 1 ) * ldc,
            ioffwork + (j - 1 ) * ldwork );
        }
        Blas3.dtrmm( 'Right', 'Upper', 'Transpose', 'Unit', lastc, k,
          1., V, ldv, Work, ldwork, ioffv, ioffwork );
        if ( lastv > k ) {
          Blas3.dgemm( 'No transpose', 'Transpose', lastc, k,
            lastv - k, 1., C, ldc, V, ldv, 1., Work, ldwork,
            ioffc + k * ldc, ioffv + k * ldv, ioffwork );
        }
        Blas3.dtrmm( 'Right', 'Upper', trans, 'Non-unit', lastc, k, 1.,
          T, ldt, Work, ldwork, iofft, ioffwork );
        if ( lastv > k ) {
          Blas3.dgemm( 'No transpose', 'No transpose', lastc,
            lastv - k, k, -1., Work, ldwork, V, ldv, 1., C, ldc,
            ioffwork, ioffv + k * ldv, ioffc + k * ldc);
        }
        Blas3.dtrmm( 'Right', 'Upper', 'No transpose', 'Unit', lastc,
          k, 1., V, ldv, Work, ldwork, ioffv, ioffwork );
        for ( j = 1; j <= k; j ++ ) {
          for (i = 1; i <= lastc; i ++ ) {
            C[ ioffc + i - 1 + ( j - 1 ) * ldc ] -=
              Work[ ioffwork + i - 1 + ( j - 1 ) * ldwork ];
          }
        }
      }
    } else {
      if ( side.charAt(0).toUpperCase() == 'L' ) {
//          V=[V1,V2] kx(m-k), kxk unit lower
//          C=[C1] (m-k)xn
//            [C2] kxn
        lastv = Math.max( k, LaPack0.iladlc( k, m, V, ldv, ioffv ) );
        lastc = LaPack0.iladlc( lastv, n, C, ldc, ioffc );
        for ( j = 1; j <= k; j ++ ) {
          Blas1.dcopy( lastc, C, ldc, Work, 1,
            ioffc + lastv - k + j - 1, ioffwork + ( j - 1 ) * ldwork );
        }
        Blas3.dtrmm( 'Right', 'Lower', 'Transpose', 'Unit', lastc, k,
          1., V, ldv, Work, ldwork, ioffv + ( lastv - k ) * ldv,
          ioffwork );
        if ( lastv > k ) {
          Blas3.dgemm( 'Transpose', 'Transpose', lastc, k,
            lastv - k, 1., C, ldc, V, ldv, 1., Work, ldwork, ioffc,
            ioffv, ioffwork );
        }
        Blas3.dtrmm( 'Right', 'Lower', transt, 'Non-unit', lastc, k,
          1., T, ldt, Work, ldwork, iofft, ioffwork );
        if ( lastv > k ) {
          Blas3.dgemm( 'Transpose', 'Transpose', lastv - k, lastc, k,
            -1., V, ldv, Work, ldwork, 1., C, ldc, ioffv, ioffwork,
            ioffc );
        }
        Blas3.dtrmm( 'Right', 'Lower', 'No transpose', 'Unit', lastc,
          k, 1., V, ldv, Work, ldwork, ioffv + ( lastv - k ) * ldv,
          ioffwork );
        for ( j = 1; j <= k; j ++ ) {
          for ( i = 1; i <= lastc; i ++ ) {
            C[ ioffc + lastv - k + j - 1 + ( i - 1 ) * ldc ] -=
              Work[ ioffwork + i - 1 + ( j - 1 ) * ldwork ];
          }
        }
      } else if ( side.charAt(0).toUpperCase() == 'R' ) {
//          V=[V1,V2] kx(n-k), kxk unit lower
//          C=[C1,C2] mx(n-k), mxk
        lastv = Math.max( k, LaPack0.iladlc( k, n, V, ldv, ioffv ) );
        lastc = LaPack0.iladlr( m, lastv, C, ldc, ioffc );
        for ( j = 1; j <= k; j ++ ) {
          Blas1.dcopy( lastc, C, 1, Work, 1,
            ioffc + ( lastv - k + j - 1 ) * ldc,
            ioffwork + ( j - 1 ) * ldwork );
        }
        Blas3.dtrmm( 'Right', 'Lower', 'Transpose', 'Unit', lastc, k,
          1., V, ldv, Work, ldwork, ioffv + ( lastv - k ) * ldv,
          ioffwork );
        if ( lastv > k ) {
          Blas3.dgemm( 'No transpose', 'Transpose', lastc, k,
            lastv - k, 1., C, ldc, V, ldv, 1., Work, ldwork, ioffc,
            ioffv, ioffwork );
        }
        Blas3.dtrmm( 'Right', 'Lower', trans, 'Non-unit', lastc, k, 1.,
          T, ldt, Work, ldwork, iofft, ioffwork );
        if ( lastv > k ) {
          Blas3.dgemm( 'No transpose', 'No transpose', lastc,
            lastv - k, k, -1., Work, ldwork, V, ldv, 1., C, ldc,
            ioffwork, ioffv, ioffc );
        }
        Blas3.dtrmm( 'Right', 'Lower', 'No transpose', 'Unit', lastc,
          k, 1., V, ldv, Work, ldwork, ioffv + ( lastv - k ) * ldv,
          ioffwork );
        for ( j = 1; j <= k; j ++ ) {
          for ( i = 1; i <= lastc; i ++ ) {
            C[ ioffc + i - 1 + ( lastv - k + j - 1 ) * ldc ] -=
              Work[ ioffwork + i - 1 + ( j - 1 ) * ldwork ];
          }
        }
      }
    }
  }
}
//************************************************************************
LaPack1.dlarfg = function( n, alphaReference, x, incx,
tauReference, ioffx ) {
  if ( n <= 1 ) {
    tauReference.setValue( 0. );
    return;
  }
  var xnorm = Blas1.dnrm2( n - 1, x, incx, ioffx );
  if ( xnorm == 0. ) tauReference.setValue( 0. );
  else {
    var beta = - LaPack0.dlapy2( alphaReference.getValue(), xnorm )
     * ( alphaReference.getValue() >= 0. ? 1. : -1. );
    var safmin = LaPack0.dlamch( 'S' ) / LaPack0.dlamch( 'E' );
    var knt = 0;
    if ( Math.abs( beta ) < safmin ) {
      var rsafmn = 1. / safmin;
      do {
        knt ++;
        Blas1.dscal( n - 1, rsafmn, x, incx, ioffx );
        beta *= rsafmn;
        alphaReference.setValue( alphaReference.getValue() * rsafmn );
      } while ( Math.abs( beta ) < safmin );
      xnorm = Blas1.dnrm2( n - 1, x, incx, ioffx );
      beta = - LaPack0.dlapy2( alphaReference.getValue(), xnorm )
        * ( alphaReference.getValue() >= 0. ? 1. : -1. );
    }
    tauReference.setValue( ( beta - alphaReference.getValue() ) / beta );
    Blas1.dscal( n - 1, 1. / ( alphaReference.getValue() - beta ), x,
      incx, ioffx );
    for ( var j = 1; j <= knt; j ++ ) beta *= safmin;
    alphaReference.setValue( beta );
  }
}
LaPack1.zlarfg = function( n, alpha, x, incx, tauReference ) {
  throw new Error("not programmed: complex array");
}
//************************************************************************
LaPack1.dlarfgp = function( n, alphaReference, x, incx,
tauReference, ioffx ) {
  if ( n <= 0 ) {
    tauReference.setValue( 0. );
    return;
  }
  var xnorm = Blas1.dnrm2( n - 1, x, incx, ioffx );
  if ( xnorm == 0. ) {
    if ( alphaReference.getValue() >= 0. ) tauReference.setValue( 0. );
    else {
      tauReference.setValue( 2. );
      for ( var j = 1; j <= n - 1; j ++ ) {
        x[ ioffx + ( j - 1 ) * incx ] = 0.;
      }
      alphaReference.setValue( - alphaReference.getValue() );
    }
  } else {
    var beta = LaPack0.dlapy2( alphaReference.getValue(), xnorm )
      * ( alphaReference.getValue() >= 0. ? 1. : -1. );
    var smlnum = LaPack0.dlamch( 'S' ) / LaPack0.dlamch( 'E' );
    var knt = 0;
    if ( Math.abs( beta ) < smlnum ) {
      var bignum = 1. / smlnum;
      do {
        knt ++;
        Blas1.dscal( n - 1, bignum, x, incx, ioffx );
        beta *= bignum;
        alphaReference.setValue( alphaReference.getValue() * bignum );
      } while ( Math.abs( beta ) < smlnum );
      xnorm = Blas1.dnrm2( n - 1, x, incx, ioffx );
      beta = LaPack0.dlapy2( alphaReference.getValue(), xnorm )
        * ( alphaReference.getValue() >= 0. ? 1. : -1. );
    }
    var savealphaReference = alphaReference.getValue();
    alphaReference.setValue( alphaReference.getValue() + beta );
    if ( beta < 0. ) {
      beta = - beta;
      tauReference.setValue( - alphaReference.getValue() / beta );
    } else {
      alphaReference.setValue( xnorm *
        ( xnorm / alphaReference.getValue() ) );
      tauReference.setValue( alphaReference.getValue() / beta );
      alphaReference.setValue( - alphaReference.getValue() );
    }
    if ( Math.abs( tauReference.getValue() ) <= smlnum ) {
      if ( savealphaReference >= 0. ) tauReference.getValue() = 0.;
      else {
        tauReference.setValue( 2. );
        for ( j = 1; j <= n - 1; j ++ ) {
          x[ ioffx + ( j - 1 ) * incx ] = 0.;
        }
        beta = - savealphaReference;
      }
    } else {
      Blas1.dscal( n - 1, 1. / alphaReference.getValue(), x, incx,
        ioffx );
    }
    for ( j = 1; j <= knt; j ++ ) beta *= smlnum;
    alphaReference.setValue( beta );
  }
}
//************************************************************************
LaPack1.dlarrb = function( n, d, lld, ifirst, ilast, rtol1,
rtol2, offset, w, wgap, werr, work, iwork, pivmin, spdiam, twist, info,
ioffd, iofflld, ioffw, ioffwgap, ioffwerr, ioffwork, ioffiwork) {
  throw new Error("not tested: unsure what lld is");
  info.setValue( 0 );
  var maxitr = Math.round( ( Math.log( spdiam + pivmin )
    - Math.log( pivmin ) ) / Math.log( 2. ) ) + 2;
  var mnwdth = 2. * pivmin;
  var r = twist;
  if ( r < 1 || r > n ) r = n;
  var i1 = ifirst;
  var nint = 0;
  var prev = 0;
  var rgap = wgap( ioffwgap + i1 - offset - 1 );
  for ( var i = i1; i <= ilast; i ++ ) {
    var k = 2 * i;
    var ii = i - offset;
    var left = w[ ioffw + ii - 1 ] - werr[ ioffwerr + ii - 1 ];
    var right = w[ ioffw + ii - 1 ] + werr[ ioffwerr + ii - 1 ];
    var lgap = rgap;
    rgap = wgap[ ioffwgap + ii - 1 ];
    var gap = Math.min( lgap, rgap );
    var back = werr[ ioffwerr + ii - 1 ];
    while ( true ) { // 20
      var negcnt =
        LaPack0.dlaneg( n, d, lld, left, pivmin, r, ioffd, iofflld );
      if ( negcng > i - 1 ) {
        left -= back;
        back *= 2.;
      } else break;
    }
    back = werr[ ioffwerr + ii - 1 ];
    while ( true ) { // 50
      negcnt =
        LaPack0.dlaneg( n, d, lld, right, pivmin, r, ioffd, iofflld);
      if ( negcng < i ) {
        right += back;
        back *= 2.;
      } else break;
    }
    width = 0.5 * Math.abs( left - right );
    var tmp = Math.max( Math.abs( left), Math.abs( right ) );
    var cvrgd = Math.max( rtol1 * gap, rtol2 * tmp );
    if ( width <= cvrgd || width <= mnwdth ) {
      iwork[ ioffiwork + k - 2 ] = -1;
      if ( i == i1 && i < ilast ) i1 = i + 1;
      if ( prev >= i1 && i <= ilast ) {
        iwork[ ioffiwork + 2 * prev - 2 ] = i + 1;
      }
    } else {
      prev = i;
      nint ++;
      iwork[ ioffiwork + k - 2 ] = i + 1;
      iwork[ ioffiwork + k - 1 ] = negcnt;
    }
    work[ ioffwork + k - 2 ] = left;
    work[ ioffwork + k - 1 ] = right;
  } // 75
  var iter = 0;
  do { // 80
    prev = ii - 1
    i = i1;
    var olnint = nint;
    for ( var ip = 1; ip <= olnint; ip ++ ) {
      k = 2 * i;
      ii = i - offset;
      rgap = wgap[ ioffwgap + ii - 1 ];
      lgap = rgap;
      if ( ii > 1 ) lgap = wgap[ ioffwgap + ii - 2 ];
      gap = Math.min( lgap, rgap );
      var next = iwork[ ioffiwork + k - 2 ];
      left = work[ ioffwork + k - 2 ];
      right = work[ ioffwork + k - 1 ];
      var mid = 0.5 * ( left + right );
      var width = right - mid;
      tmp = Math.max( Math.abs( left ), Math.abs( right ) );
      cvrgd = Math.max( rtol1 * gap, rtol2 * tmp );
      if ( width <= cvrgd || width <= mnwdth || iter == maxitr ) {
        nint --;
        iwork[ ioffiwork + k - 2 ] = 0;
        if ( i1 == i ) i1 = next;
        else if ( prev >= i1 ) {
          iwork[ ioffiwork + 2 * prev - 2 ] = next;
        }
        i = next;
        continue;
      }
      prev = i;
      negcnt =
        LaPack0.dlaneg( n, d, lld, mid, pivmin, r, ioffd, iofflld);
      if ( negcnt <= i - 1 ) work[ ioffwork + k - 2 ] = mid;
      else work[ ioffwork + k - 1 ] = mid;
      i = next;
    } // 100
    iter ++;
  } while ( nint > 0 && iter <= maxitr );
  for ( i = ifirst; i <= ilast; i ++ ) {
    k = 2 * i;
    ii = i - offset;
    if ( iwork[ ioffiwork + k - 2 ] == 0 ) {
      w[ ioffw + ii - 1 ] = 0.5 * ( work[ ioffwork + k - 2 ]
        + work[ ioffwork + k - 1 ] );
      werr[ ioffwerr + ii - 1 ] = work[ ioffwork + k - 1 ]
        - w[ ioffw + ii - 1 ];
    }
  }
  for ( i = ifirst + 1; i <= ilast; i ++ ) {
    k = 2 * i;
    ii = i - offset;
    wgap[ ioffwgap + ii - 2 ] = Math.max( 0.,
      w[ ioffw + ii - 1 ] - werr[ ioffwerr + ii - 1 ]
      - w[ ioffw + ii - 2 ] - werr[ ioffwerr + ii - 2 ] );
  }
}
//************************************************************************
LaPack1.dlarrd = function( range, order, n, vl, vu, il, iu,
gers, reltol, d, e, e2, pivmin, nsplit, isplit, m, w, werr, wlReference,
wuReference, iblock, indexw, work, iwork, info, ioffgers, ioffd, ioffe,
ioffe2, ioffisplit, ioffw, ioffwerr, ioffiblock, ioffindexw, ioffwork,
ioffiwork ) {
  throw new Error("not tested: complicated input");
  var fudge = 2.;
  var allrng = 1;
  var valrng = 2;
  var indrng = 3;
  var idumma = new Array( 1 );
  info.setValue( 0 );
  var irange = 0;
  if ( range.charAt(0).toUppeCase() == 'A' ) irange = allrng;
  else if ( range.charAt(0).toUppeCase() == 'V' ) irange = valrng;
  else if ( range.charAt(0).toUppeCase() == 'I' ) irange = indrng;
  if ( irange <= 0 ) info.setValue( -1 );
  else if ( order.charAt(0).toUpperCase() != 'B' &&
  order.charAt(0).toUpperCase() != 'E' ) {
    info.setValue( -2 );
  } else if ( n < 0 ) info.setValue( -3 );
  else if ( irange == valrng ) {
    if ( vl >= vu ) info.setValue( -5 );
  } else if ( irange == indrng && ( il < 1 || il > Math.max( 1, n ) ) )
  {
    info.setValue( -6 );
  } else if ( irange == indrng &&
  ( iu < Math.min( n, il ) || iu > n ) ) {
    info.setValue( -7 );
  }
  if ( info.getValue() != 0 ) return;
  info.setValue( 0 );
  var ncnvrg = false;
  var toofew = false;
  m.setValue( 0 );
  if ( n == 0 ) return;
  if ( irange == indrng && il == 1 && iu == n ) irange = 1;
  var eps = LaPack0.dlamch( 'P' );
  var uflow = LaPack0.dlamch( 'U' );
  if ( n == 1 ) {
    if ( irange == allrng ||
    ( irange == valrng && d[ ioffd ] > vl && d[ ioffd ] <= vu ) ||
    ( irange == indrng && il == 1 && iu == 1 ) ) {
      m.setValue( 1 );
      w[ ioffw ] = d[ ioffd ];
      werr[ ioffwerr ] = 0.;
      iblock[ ioffiblock ] = 1;
      indexw[ ioffindexw ] = 1;
    }
    return;
  }
  var nb = LaPack0.ilaenv( 1, 'dstebz', ' ', n, -1, -1, -1 );
  if ( nb <= 1 ) nb = 0;
  var gl = d[ ioffd ];
  var gu = d[ ioffd ];
  for ( var i = 1; i <= n; i ++ ) {
    gl = Math.min( gl, gers[ ioffgers + 2 * i - 2 ] );
    gu = Math.max( gu, gers[ ioffgers + 2 * i - 1 ] );
  } // 5
  var tnorm = Math.max( Math.abs( gl ), Math.abs( gu ) );
  gl -= fudge * tnorm * eps * n + fudge * 2. * pivmin;
  gu += fudge * tnorm * eps * n + fudge * 2. * pivmin;
  var rtoli = reltol;
  var atoli = fudge * 2. * uflow + fudge * 2. * pivmin;
  if ( irange == indrng ) {
    var itmax = Math.round( ( Math.log( tnorm + pivmin )
      - Math.log( pivmin ) ) / Math.log( 2. ) ) + 2;
    work[ ioffwork + n ] = gl;
    work[ ioffwork + n + 1 ] = gl;
    work[ ioffwork + n + 2 ] = gu;
    work[ ioffwork + n + 3 ] = gu;
    work[ ioffwork + n + 4 ] = gl;
    work[ ioffwork + n + 5 ] = gu;
    iwork[ ioffiwork ] = -1;
    iwork[ ioffiwork + 1 ] = -1;
    iwork[ ioffiwork + 2 ] = n + 1;
    iwork[ ioffiwork + 3 ] = n + 1;
    iwork[ ioffiwork + 4 ] = il - 1;
    iwork[ ioffiwork + 5 ] = iu;
    var iout = new IntReference();
    var iinfo = new IntReference();
    LaPack0.dlaebz( 3, itmax, n, 2, 2, nb, atoli, rtoli, pivmin,
      d, e, e2, iwork, work, work, iout, iwork, w, iblock, iinfo,
      ioffd, ioffe, ioffe2, ioffwork + 4, ioffwork + n,
      ioffwork + n + 4, ioffiwork, ioffw, ioffiblock );
    if ( iinfo.getValue() != 0 ) {
      info.setValue( iinfo.getValue() );
      return;
    }
    if ( iwork[ ioffiwork + 5 ] == iu ) {
      wl.setValue( work[ ioffwork + n ] );
      var wlu = work[ ioffwork + n + 2 ];
      var nwl = iwork[ ioffiwork ];
      wu.setValue( work[ ioffwork + n + 3 ] );
      var wul = work[ ioffwork + n + 1 ];
      var nwu = iwork[ ioffiwork + 3 ];
    } else {
      wl.setValue( work[ ioffwork + n + 1 ] );
      wlu = work[ ioffwork + n + 3 ];
      nwl = iwork[ ioffiwork + 1 ];
      wu.setValue( work[ ioffwork + n + 2 ] );
      wul = work[ ioffwork + n ];
      nwu = iwork[ ioffiwork + 2 ];
    }
    if ( nwl < 0 || nwl >= n || nwu < 1 || nwu > n ) {
      info.setValue( 4 );
      return;
    }
  } else if ( irange == valrng ) {
    wl.setValue( vl );
    wu.setValue( vu );
  } else if ( irange == allrng ) {
    wl.setValue( gl );
    wu.setValue( gu );
  }
  m.setValue( 0 );
  var iend = 0;
  info.setValue( 0 );
  nwl = 0;
  nwu = 0;
  for ( var jblk = 1; jblk <= nsplit; jblk ++ ) {
    var ioff = iend;
    var ibegin = ioff + 1;
    iend = isplit[ ioffisplit + jblk - 1 ];
    var in2 = iend - ioff;
    if ( in2 == 1 ) {
      if ( wl.getValue() >= d[ ioffd + ibegin - 1 ] - pivmin ) nwl ++;
      if ( wu.getValue() >= d[ ioffd + ibegin - 1 ] - pivmin ) nwu ++;
      if ( irange == allrng ||
      ( wl.getValue() < d[ ioffd + ibegin - 1 ] - pivmin &&
      wu.getValue() >= d[ ioffd + ibegin - 1 ] - pivmin ) ) {
        m.setValue( m.getValue() + 1 );
        w[ ioffw + m.getValue() - 1 ] = d[ ioffd + ibegin - 1 ];
        werr[ ioffwerr + m.getValue() - 1 ] = 0.;
        iblock[ ioffiblock + m.getValue() - 1 ] = jblk;
        indexw[ ioffindexw + m.getValue() - 1 ] = 1;
      }
    } else {
      gu = d[ ioffd + ibegin - 1 ];
      gl = d[ ioffd + ibegin - 1 ];
      var tmp1 = 0.;
      for ( var j = ibegin; j <= iend; j ++ ) {
        gl = Math.min( gl, gers[ ioffgers + 2 * j - 2 ] );
        gu = Math.max( gu, gers[ ioffgers + 2 * j - 1 ] );
      } // 40
      gl -= fudge * tnorm * eps * in2 + fudge * pivmin;
      gu += fudge * tnorm * eps * in2 + fudge * pivmin;
      if ( irange > 1 ) {
        if ( gu < wl.getValue() ) {
          nwl += in2;
          nwu += in2;
          continue;
        }
        gl = Math.max( gl, wl.getValue() );
        gu = Math.min( gu, wu.getValue() );
        if ( gl >= gu ) continue;
      }
      work[ ioffwork + n ] = gl;
      work[ ioffwork + n + in2 ] = gu;
      var im = new IntReference();
      LaPack0.dlaebz( 1, 0, in2, in2, 1, nb, atoli, rtoli, pivmin,
        d, e, e2, idumma, work, work, im, iwork, w, iblock, iinfo,
        ioffd + ibegin - 1, ioffe + ibegin - 1, ioffe2 + ibegin - 1,
        0, ioffwork + n, ioffwork + n + 2 * in2, ioffiwork, 
        ioffw + m.getValue(), ioffiblock + m.getValue() );
      if ( iinfo.getValue() != 0 ) {
        info.setValue( iinfo.getValue() );
        return;
      }
      nwl += iwork[ ioffiwork ];
      nwu += iwork[ ioffiwork + in2 ];
      var iwoff = m.getValue() - iwork[ ioffiwork ];
      itmax = Math.round( ( Math.log( gu - gl + pivmin )
        - Math.log( pivmin ) ) / Math.log( 2. ) ) + 2;
      LaPack0.dlaebz( 2, itmax, in2, in2, 1, nb, atoli, rtoli, pivmin,
        d, e, e2, idumma, work, work, iout, iwork, w, iblock, iinfo,
        ioffd + ibegin - 1, ioffe + ibegin - 1, ioffe2 + ibegin - 1,
        0, ioffwork + n, ioffwork + n + 2 * in2, ioffiwork, 
        ioffw + m.getValue(), ioffiblock + m.getValue() );
      if ( iinfo.getValue() != 0 ) {
        info.setValue( iinfo.getValue() );
        return;
      }
      for ( j = 1; j <= iout.getValue(); j ++ ) {
        tmp1 = 0.5 * ( work[ ioffwork + j + n - 1 ]
          + work[ ioffwork + j + in2 + n - 1 ] );
        var tmp2 =
          0.5 * Math.abs( work[ ioffwork + j + n - 1 ]
          - work[ ioffwork + j + in2 + n - 1 ] );
        if ( j > iout.getValue() - iinfo.getValue() ) {
          ncnvrg = true;
          var ib = - jblk;
        } else ib = jblk;
        for ( var je = iwork[ ioffiwork + j - 1 ] + 1 + iwoff;
        je <= iwork[ ioffiwork + j + in2 - 1 ] + iwoff; je ++ ) {
          w[ ioffw + je - 1 ] = tmp1;
          werr[ ioffwerr + je - 1 ] = tmp2;
          indexw[ ioffindexw + je - 1 ] = je - iwoff;
          iblock[ ioffiblock + je - 1 ] = ib;
        } // 50
      } // 60
      m.setValue( m.getValue() + im.getValue() );
    }
  } // 70
  if ( irange == indrng ) {
    var idiscl = il - 1 - nwl;
    var idiscu = nwu - iu;
    if ( idiscl > 0 ) {
      im.setValue( 0 );
      for ( je = 1; je <= m.getValue(); je ++ ) {
        if ( w[ ioffw + je - 1 ] <= wlu && idiscl > 0 ) idiscl --;
        else {
          im.setValue( im.getValue() + 1 );
          w[ ioffw + im.getValue() - 1 ] = w[ ioffw + je - 1 ];
          werr[ ioffwerr + im.getValue() - 1 ] = werr[ ioffwerr + je - 1 ];
          indexw[ ioffindexw + im.getValue() - 1 ] =
            indexw[ ioffindexw + je - 1 ];
          iblock[ ioffiblock + im.getValue() - 1 ] =
            iblock[ ioffiblock + je - 1 ];
        }
      } // 80
      m.setValue( im.getValue() );
    }
    if ( idiscu > 0 ) {
      im.setValue( m.getValue() + 1 );
      for ( je = m.getValue(); je >= 1; je -- ) {
        if ( w[ ioffw + je - 1 ] >= wul && idiscu > 0 ) idiscu --;
        else {
          im.setValue( im.getValue() - 1 );
          w[ ioffw + im.getValue() - 1 ] = w[ ioffw + je - 1 ];
          werr[ ioffwerr + im.getValue() - 1 ] = werr[ ioffwerr + je - 1 ];
          indexw[ ioffindexw + im.getValue() - 1 ] =
            indexw[ ioffindexw + je - 1 ];
          iblock[ ioffiblock + im.getValue() - 1 ] =
            iblock[ ioffiblock + je - 1 ];
        }
      } // 81
      var jee = 0;
      for ( je = im.getValue(); je <= m.getValue(); je ++ ) {
        jee ++;
        w[ ioffw + jee - 1 ] = w[ ioffw + je - 1 ];
        werr[ ioffwerr + jee - 1 ] = werr[ ioffwerr + je - 1 ];
        indexw[ ioffindexw + jee - 1 ] =
          indexw[ ioffindexw + je - 1 ];
        iblock[ ioffiblock + jee - 1 ] =
          iblock[ ioffiblock + je - 1 ];
      } // 82
      m.setValue( m.getValue() - im.getValue() + 1 );
    }
    if ( idiscl > 0 || idiscu > 0 ) {
      if ( idiscl > 0 ) {
        var wkill = wu.getValue();
        for ( var jdisc = 1; jdisc <= idiscl; jdisc ++ ) {
          var iw = 0;
          for ( je = 1; je <= m.getValue(); je ++ ) {
            if ( iblock[ ioffiblock + je - 1 ] != 0 &&
            ( w[ ioffw + je - 1 ] < wkill || iw == 0 ) ) {
              iw = je;
              wkill = w[ ioffw + je - 1 ];
            }
          } // 90
          iblock[ ioffiblock + iw - 1 ] = 0;
        } // 100
      }
      if ( idiscu > 0 ) {
        wkill = wl.getValue();
        for ( jdisc = 1; jdisc <= idiscu; jdisc ++ ) {
          iw = 0;
          for ( je = 1; je <= m.getValue(); je ++ ) {
            if ( iblock[ ioffiblock + je - 1 ] != 0 &&
            ( w[ ioffw + je - 1 ] >= wkill || iw == 0 ) ) {
              iw = je;
              wkill = w[ ioffw + je - 1 ];
            }
          } // 110
          iblock[ ioffiblock + iw - 1 ] = 0;
        } // 120
      }
      im.setValue( 0 );
      for ( je = 1; je <= m.getValue(); je ++ ) {
        if ( iblock[ ioffiblock + je - 1 ] != 0 ) {
          im.setValue( im.getValue() + 1 );
          w[ ioffw + im.getValue() - 1 ] = w[ ioffw + je - 1 ];
          werr[ ioffwerr + im.getValue() - 1 ] = werr[ ioffwerr + je - 1 ];
          indexw[ ioffindexw + im.getValue() - 1 ] =
            indexw[ ioffindexw + je - 1 ];
          iblock[ ioffiblock + im.getValue() - 1 ] =
            iblock[ ioffiblock + je - 1 ];
        }
      } // 130
      m.setValue( im.getValue() );
    }
    if ( idiscl < 0 || idiscu < 0 ) toofew = true;
  }
  if ( ( irange == allrng && m.getValue() != n ) ||
  ( irange == indrng && m.getValue() != iu - il + 1 ) ) {
    toofew = true;
  }
  if ( order.charAt(0).toUpperCase() == 'E' && nsplit > 1 ) {
    for ( je = 1; je <= m.getValue() - 1; je ++ ) {
      var ie = 0;
      tmp1 = w[ ioffw + je - 1 ];
      for ( j = je + 1; je <= m.getValue(); je ++ ) {
        if ( w[ ioffw + j - 1 ] < tmp1 ) {
          ie = j;
          tmp1 = w[ ioffw + j - 1 ];
        }
      } // 140
      if ( ie != 0 ) {
        tmp2 = werr[ ioffwerr + ie - 1 ];
        var itmp1 = iblock[ ioffiblock + ie - 1 ];
        var itmp2 = indexw[ ioffindexw + ie - 1 ];
        w[ ioffw + ie - 1 ] = w[ ioffw + je - 1 ];
        werr[ ioffwerr + ie - 1 ] = werr[ ioffwerr + je - 1 ];
        iblock[ ioffiblock + ie - 1 ] = iblock[ ioffiblock + je - 1 ];
        indexw[ ioffindexw + ie - 1 ] = indexw[ ioffindexw + je - 1 ];
        w[ ioffw + je - 1 ] = tmp1;
        werr[ ioffwerr + je - 1 ] = tmp2;
        iblock[ ioffiblock + je - 1 ] = itmp1;
        indexw[ ioffindexw + je - 1 ] = itmp2;
      }
    } // 150
  }
  info.setValue( 0 );
  if ( ncnvrg ) info.getValue() = info.getValue() + 1;
  if ( toofew ) info.getValue() = info.getValue() + 2;
}
//************************************************************************
LaPack1.dlarrf = function( n, d, l, ld, clstrt, clend, w, wgap,
werr, spdiam, clgapl, clgapr, pivmin, sigmaReference, dplus, lplus, work,
info, ioffd, ioffl, ioffld, ioffw, ioffwgap, ioffwerr, ioffdplus,
iofflplus, ioffwork ) {
  throw new Error("not tested: complicated input");
  var maxgrowth1 = 8.;
  var maxgrowth2 = 8.;
  var ktrymax = 1;
  var sleft = 1;
  var sright = 2;
  info.setValue( 0 );
  var fact = Math.pow( 2., ktrymax );
  var eps = LaPack0.dlamch( 'Precision' );
  var shift = 0;
  var forcer = false;
  var nofail = true;
  var clwdth = Math.abs( w[ ioffw + clend -1 ]
    - w[ ioffw + clstrt - 1 ] ) + werr[ ioffwerr + clend -1 ]
    - werr[ ioffwerr + clstrt - 1 ];
  var avgap = clwdth / Number( clend - clstrt );
  var mingap = Math.min( clgapl, clgapr );
  var lsigma = Math.min( w[ ioffw + clstrt - 1 ],
    w[ ioffw + clend - 1 ] ) - werr[ ioffwerr + clstrt - 1 ];
  var rsigma = Math.max( w[ ioffw + clstrt - 1 ],
    w[ ioffw + clend - 1 ] ) + werr[ ioffwerr + clend - 1 ];
  lsigma -= Math.abs( lsigma ) * 4. * eps;
  rsigma += Math.abs( rsigma ) * 4. * eps;
  var ldmax = 0.25 * mingap + 2. * pivmin;
  var rdmax = 0.25 * mingap + 2. * pivmin;
  var ldelta = Math.max( avgap, wgap[ ioffwgap + clstrt - 1 ] )
    / fact;
  var rdelta = Math.max( avgap, wgap[ ioffwgap + clend - 2 ] )
    / fact;
  var s = LaPack0.dlamch( 'S' );
  var smlgrowth = 1. / s;
  var fail = Number( n - 1 ) * mingap / ( spdiam * eps );
  var fail2 = Number( n - 1 ) * mingap
    / ( spdiam * Math.sqrt( eps ) );
  var bestshift = lsigma;
  var ktry = 0;
  var growthbound = maxgrowth1 * spdiam;
  while ( true ) { // 5
    var goto50 = false;
    var sawnan1 = false;
    var sawnan2 = false;
    ldelta = Math.min( ldmax, ldelta );
    rdelta = Math.min( rdmax, rdelta );
    s = - lsigma;
    dplus[ ioffdplus ] = d[ ioffd ] + s;
    if ( Math.abs( dplus[ ioffdplus ] ) < pivmin ) {
      dplus[ ioffdplus ] = - pivmin;
      sawnan1 = true;
    }
    var max1 = Math.abs( dplus[ ioffdplus ] );
    for ( var i = 1; i <= n - 1; i ++ ) {
      lplus[ iofflplus + i - 1 ] = ld[ ioffld + i - 1 ]
        / dplus[ ioffdplus + i - 1 ];
      s = s * lplus[ iofflplus + i - 1 ] * l[ ioffl + i - 1 ] - lsigma;
      dplus[ ioffdplus + i ] = d[ ioffd + i ] + s;
      if ( Math.abs( dlplus[ ioffdplus + i ] ) < pivmin ) {
        dplus[ ioffdplus + i ] = - pivmin;
        sawnan1 = true;
      }
      max1 = Math.max( max1, Math.abs( dplus[ ioffdplus + i ] ) );
    } // 6
    sawnan1 = sawnan1 || isNaN( max1 );
    if ( forcer || ( max1 <= growthbound && ! sawnan1 ) ) {
      sigma.setValue( lsigma );
      shift = sleft;
      break;
    }
    s = - rsigma;
    work[ ioffwork ] = d[ ioffd ] + s;
    if ( Math.abs( work[ ioffwork ] ) < pivmin ) {
      work[ ioffwork ] = - pivmin;
      sawnan2 = true;
    }
    var max2 = Math.abs( work[ ioffwork ] );
    for ( i = 1; i <= n - 1; i ++ ) {
      work[ ioffwork + n + i - 1 ] =
        ld[ ioffld + i - 1 ] / work[ ioffwork + i - 1 ];
      s = s * work[ ioffwork + n + i - 1 ] * l[ ioffl + i - 1 ]
        - rsigma;
      work[ ioffwork + i ] = d[ ioffd + i ] + s;
      if ( Math.abs( work[ ioffwork + i ] ) < pivmin ) {
        work[ ioffwork + i ] = - pivmin;
        sawnan2 = true;
      }
      max2 = Math.max( max2, Math.abs( work[ iofffwork + i ] ) );
    } // 7
    sawnan2 = sawnan2 || isNaN( max2 );
    if ( forcer || ( max2 <= growthbound && ! sawnan2 ) ) {
      sigma.setValue( rsigma );
      shift = sright;
      break;
    }
    if ( sawnan1 && sawnan2 ) goto50 = true;
    else {
      if ( ! sawnan1 ) {
        var indx = 1;
        if ( max1 <= smlgrowth ) {
          smlgrowth = max1;
          bestshift = lsigma;
        }
      }
      if ( ! sawnan2 ) {
        if ( sawnan1 || max2 <= max1 ) indx = 2;
        if ( max2 <= smlgrowth ) {
          smlgrowth = max2;
          bestshift = rsigma;
        }
      }
    }
    if ( ! goto50 ) {
      var dorrr1 = ( clwdth < mingap / 128.
        && Math.min( max1, max2 ) < fail2 && ! sawnan1 && ! sawnan2 );
      var tryrrr1 = true;
      if ( tryrrr1 && dorrr1 ) {
        if ( indx == 1 ) {
          var tmp = Math.abs( dplus[ ioffdplus + n - 1 ] );
          var znm2 = 1.;
          var prod = 1.;
          var oldp = 1.;
          for ( i = n - 1; i >= 1; i -- ) {
            if ( prod <= eps ) {
              prod = ( ( dplus[ ioffdplus + i ]
                * work[ ioffwork + n + i ] )
                / ( dplus[ ioffdplus + i - 1 ]
                * work[ ioffwork + n + i - 1 ] ) ) * oldp;
            } else {
              prod *= Math.abs( work[ ioffwork + n + i - 1 ] );
            }
            oldp = prod;
            znm2 += prod * prod;
            tmp = Math.max( tmp,
              Math.abs( dplus[ ioffdplus + i - 1 ] * prod ) );
          } // 15
          rrr1 = tmp / ( spdiam * Math.sqrt( znm2 ) );
          if ( rrr1 <= maxgrowth2 ) {
            sigma.setValue( lsigma );
            shift = sleft;
            break;
          }
        } else if ( indx == 2 ) {
          tmp = Math.abs( work[ ioffwork + n - 1 ] );
          znm2 = 1.;
          prod = 1.;
          oldp = 1.;
          for ( i = n - 1; i >= 1; i -- ) {
            if ( prod <= eps ) {
              prod = ( ( work[ ioffwork + i ]
                * lplus[ iofflplus + i ] )
                / ( work[ ioffwork + i - 1 ]
                * lplus[ iofflplus + i - 1 ] ) ) * oldp;
            } else {
              prod *= Math.abs( lplus[ iofflplus + i - 1 ] );
            }
            oldp = prod;
            znm2 += prod * prod;
            tmp = Math.max( tmp,
              Math.abs( work[ ioffwork + i - 1 ] * prod ) );
          } // 16
          var rrr2 = tmp / ( spdiam * Math.sqrt( znm2 ) );
          if ( rrr2 <= maxgrowth2 ) {
            sigma.setValue( rsigma );
            shift = sright;
            break;
          }
        }
      }
    } // 50
    if ( ktry < ktrymax ) {
      lsigma = Math.max( lsigma - ldelta, lsigma - ldmax );
      rsigma = Math.min( rsigma + rdelta, rsigma + rdmax );
      ldelta *= 2.;
      rdelta *= 2.;
      ktry ++
    } else {
      if ( smlgrowth < fail || nofail ) {
        lsigma = bestshift;
        rsigma = bestshift;
        forcer = true;
      } else {
        info.setValue( 1 );
        return;
      }
    }
  } // 100
  if ( shift == sleft );
  else if ( shift == sright ) {
    Blas1.dcopy( n, work, 1, dplus, 1, ioffwork, ioffdplus );
    Blas1.dcopy( n - 1, work, 1, lplus, 1, ioffwork + n, iofflplus );
  }
}
//************************************************************************
LaPack1.dlarrk = function( n, iw, gl, gu, d, e2, pivmin, reltol,
wReference, werrReference, info, ioffd, ioffe2 ) {
  var fudge = 2.;
  var eps = LaPack0.dlamch( 'P' );
  var tnorm = Math.max( Math.abs( gl ), Math.abs( gu ) );
  var rtoli = reltol;
  var atoli = fudge * 2. * pivmin;
  var itmax = Math.round( ( Math.log( tnorm + pivmin )
    - Math.log( pivmin ) ) / Math.log( 2. ) ) + 2;
  info.setValue( -1 );
  var left = gl - fudge * tnorm * eps * n - fudge * 2. * pivmin;
  var right = gu + fudge * tnorm * eps * n
    + fudge * 2. * pivmin;
  var it = 0;
  while ( true ) { // 10
    var tmp1 = Math.abs( right - left );
    var tmp2 = Math.max( Math.abs( right), Math.abs( left ) );
    if ( tmp1 < Math.max( Math.max( atoli, pivmin ), rtoli * tmp2 ) ) {
      info.setValue( 0 );
      break;
    }
    if ( it > itmax ) break;
    it ++;
    var mid = 0.5 * ( left + right );
    var negcnt = 0;
    tmp1 = d[ ioffd ] - mid;
    if ( Math.abs( tmp1 ) < pivmin ) tmp1 = - pivmin;
    if ( tmp1 <= 0. ) negcnt ++;
    for ( var i = 2; i <= n; i ++ ) {
      tmp1 = d[ ioffd + i - 1 ] - e2[ ioffe2 + i - 2 ] / tmp1 - mid;
      if ( Math.abs( tmp1 ) < pivmin ) tmp1 = - pivmin;
      if ( tmp1 <= 0. ) negcnt ++;
    } // 20
    if ( negcnt >= iw ) right = mid;
    else left = mid;
  }
  wReference.setValue( 0.5 * ( left + right ) );
  werrReference.setValue( 0.5 * Math.abs( right - left ) );
}
//************************************************************************
LaPack1.dlarrr = function( n, d, e, info, ioffd, ioffe ) {
  throw new Error("not tested: hard to check");
  var relcond = 0.999;
  info.setValue( 1 );
  var safmin = LaPack0.dlamch( 'Safe minimum' );
  var eps = LaPack0.dlamch( 'Precision' );
  var smlnum = safmin / eps;
  var rmin = Math.sqrt( smlnum );
  var yesrel = true;
  var offdig = 0.;
  var tmp = Math.sqrt( Math.abs( d[ ioffd ] ) );
  if ( tmp < rmin ) yesrel = false;
  if ( yesrel ) {
    for ( var i = 2; i <= n; i ++ ) {
      var tmp2 = Math.sqrt( Math.abs( d[ ioffd + i - 1 ] ) );
      if ( tmp2 < rmin ) {
        yesrel = false;
        break;
      }
      var offdig2 = Math.abs( e[ ioffe + i - 2 ] )
        / ( tmp * tmp2 );
      if ( offdig + offdig2 >= relcond ) {
        yesrel = false;
        break;
      }
      tmp = tmp2;
      offdig = offdig2;
    }
  }
  if ( yesrel ) {
    info.setValue( 0 );
    return;
  } else;
}
//************************************************************************
LaPack1.dlartg = function( f, g, csReference, snReference, rReference ) {
  var safmin = LaPack0.dlamch( 'S' );
  var eps = LaPack0.dlamch( 'E' );
  var safmn2 = Math.pow( LaPack0.dlamch( 'B' ),
    Math.round( Math.log( safmin / eps )
       / Math.log( LaPack0.dlamch( 'B' ) ) / 2. ) );
  var safmx2 = 1. / safmn2;
  if ( g == 0. ) {
    csReference.setValue( 1. );
    snReference.setValue( 0. );
    rReference.setValue( f );
  } else if ( f == 0. ) {
    csReference.setValue( 0. );
    snReference.setValue( 1. );
    rReference.setValue( g );
  } else {
    var f1 = f;
    var g1 = g;
    var scale = Math.max( Math.abs( f1 ), Math.abs( g1 ) );
    if ( scale >= safmx2 ) {
      var count = 0;
      do { // 10
        count ++;
        f1 *= safmn2;
        g1 *= safmn2;
        scale = Math.max( Math.abs( f1 ), Math.abs( g1 ) );
      } while ( scale >= safmx2 );
      rReference.setValue( Math.sqrt( f1 * f1 + g1 * g1 ) );
      csReference.setValue( f1 / rReference.getValue() );
      snReference.setValue( g1 / rReference.getValue() );
      for ( var i = 1; i <= count; i ++ ) {
        rReference.setValue( rReference.getValue() * safmx2 );
      } // 20
    } else if ( scale <= safmn2 ) {
      count = 0;
      do { // 30
        count ++;
        f1 *= safmx2;
        g1 *= safmx2;
        scale = Math.max( Math.abs( f1 ), Math.abs( g1 ) );
      } while ( scale <= safmn2 );
      rReference.setValue( Math.sqrt( f1 * f1 + g1 * g1 ) );
      csReference.setValue( f1 / rReference.getValue() );
      snReference.setValue( g1 / rReference.getValue() );
      for ( i = 1; i <= count; i ++ ) {
        rReference.setValue( rReference.getValue() * safmn2 );
      }
    } else {
      rReference.setValue( Math.sqrt( f1 * f1 + g1 * g1 ) );
      csReference.setValue( f1 / rReference.getValue() );
      snReference.setValue( g1 / rReference.getValue() );
    }
    if ( Math.abs( f ) > Math.abs( g ) && csReference.getValue() < 0. ) {
      csReference.setValue( - csReference.getValue() );
      snReference.setValue( - snReference.getValue() );
      rReference.setValue( - rReference.getValue() );
    }
  }
}
LaPack1.zlartg = function( f, g, cs, sn, r ) {
  throw new Error("not programmed: complex numbers");
}
//************************************************************************
LaPack1.dlartgp = function( f, g, csReference, snReference,
rReference ) {
  var safmin = LaPack0.dlamch( 'S' );
  var eps = LaPack0.dlamch( 'E' );
  var safmn2 = Math.pow( LaPack0.dlamch( 'B' ),
    Math.round( Math.log( safmin / eps )
    / Math.log( LaPack0.dlamch( 'B' ) ) / 2. ) );
  var safmx2 = 1. / safmn2;
  if ( g == 0. ) {
    csReference.setValue( ( f >= 0. ? 1.: -1. ) );
    snReference.setValue( 0. );
    rReference.setValue( Math.abs( f ) );
  } else if ( f == 0. ) {
    csReference.setValue( 0. );
    snReference.setValue( ( g >= 0. ? 1. : -1. ) );
    rReference.setValue( Math.abs( g ) );
  } else {
    var f1 = f;
    var g1 = g;
    var scale = Math.max( Math.abs( f1 ), Math.abs( g1 ) );
    if ( scale >= safmx2 ) {
      var count = 0;
      do {
        count ++;
        f1 *= safmn2;
        g1 *= safmn2;
        scale = Math.max( Math.abs( f1 ), Math.abs( g1 ) );
      } while ( scale >= safmx2 );
      rReference.setValue( Math.sqrt( f1 * f1 + g1 * g1 ) );
      csReference.setValue( f1 / rReference.getValue() );
      snReference.setValue( g1 / rReference.getValue() );
      for ( var i = 1; i <= count; i ++ ) {
        rReference.setValue( rReference.getValue() * safmx2 );
      }
    } else if ( scale <= safmn2 ) {
      count = 0;
      do {
        count ++;
        f1 *= safmx2;
        g1 *= safmx2;
        scale = Math.max( Math.abs( f1 ), Math.abs( g1 ) );
      } while ( scale <= safmn2 );
      rReference.setValue( Math.sqrt( f1 * f1 + g1 * g1 ) );
      csReference.setValue( f1 / rReference.getValue() );
      snReference.setValue( g1 / rReference.getValue() );
      for ( i = 1; i <= count; i ++ ) {
        rReference.setValue( rReference.getValue() * safmn2 );
      }
    } else {
      rReference.setValue( Math.sqrt( f1 * f1 + g1 * g1 ) );
      csReference.setValue( f1 / rReference.getValue() );
      snReference.setValue( g1 / rReference.getValue() );
    }
    if ( rReference.getValue() < 0. ) {
      csReference.setValue( - csReference.getValue() );
      snReference.setValue( - snReference.getValue() );
      rReference.setValue( - rReference.getValue() );
    }
  }
}
//************************************************************************
LaPack1.dlascl = function( type, kl, ku, cfrom, cto, m, n, A,
lda, info, ioffa ) {
  info.setValue( 0 );
  if ( type.charAt(0).toUpperCase() == 'G' ) var itype = 0;
  else if ( type.charAt(0).toUpperCase() == 'L' ) itype = 1;
  else if ( type.charAt(0).toUpperCase() == 'U' ) itype = 2;
  else if ( type.charAt(0).toUpperCase() == 'H' ) itype = 3;
  else if ( type.charAt(0).toUpperCase() == 'B' ) itype = 4;
  else if ( type.charAt(0).toUpperCase() == 'Q' ) itype = 5;
  else if ( type.charAt(0).toUpperCase() == 'Z' ) itype = 6;
  else itype = -1;
  if ( itype == -1 ) info.setValue( -1 );
  else if ( cfrom == 0. || isNaN( cfrom ) ) info.setValue( -4 );
  else if ( isNaN( cto ) ) info.setValue( -5 );
  else if ( m < 0 ) info.setValue( -6 );
  else if ( n < 0 || ( itype == 4 && n != m ) ||
  ( itype == 5 && n != m ) ) {
    info.setValue( -7 );
  } else if ( itype <= 3 && lda < Math.max( 1, m ) ) info.setValue( -9 );
  else if ( itype >= 4 ) {
    if ( kl < 0 || kl > Math.max( m - 1, 0 ) ) info.setValue( -2 );
    else if ( ku < 0 || ku > Math.max( n - 1, 0 ) ||
    ( ( itype == 4 || itype == 5 ) && kl != ku ) ) {
      info.setValue( -3 );
    } else if ( ( itype == 4 && lda < kl + 1 ) ||
    ( itype == 5 && lda < ku + 1 ) ||
    ( itype == 6 && lda < 2 * kl + ku + 1 ) ) {
      info.setValue( -9 );
    }
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlascl', -info.getValue() );
    return;
  }
  if ( n == 0 || m == 0 ) return;
  var smlnum = LaPack0.dlamch( 'S' );
  var bignum = 1. / smlnum;
  var cfromc = cfrom;
  var ctoc = cto;
  var done = false;
  while ( ! done ) { // 10
    var cfrom1 = cfromc * smlnum;
    if ( cfrom1 == cfromc ) {
      var mul = ctoc / cfromc;
      done = true;
      var cto1 = ctoc;
    } else {
      cto1 = ctoc / bignum;
      if ( cto1 == ctoc ) {
        mul = ctoc;
        done = true;
        cfromc = 1.;
      } else if ( Math.abs( cfrom1 ) > Math.abs( ctoc ) && ctoc != 0. )
      {
        mul = smlnum;
        done = false;
        cfromc = cfrom1;
      } else if ( Math.abs( cto1 ) > Math.abs( cfromc ) ) {
        mul = bignum;
        done = false;
        ctoc = cto1;
      } else {
        mul = ctoc / cfromc;
        done = true;
      }
    }
    if ( itype == 0 ) {
      for ( var j = 1; j <= n; j ++ ) {
        for ( var i = 1; i <= m; i ++ ) {
          A[ ioffa + i - 1 + ( j - 1 ) * lda ] *= mul;
        }
      }
    } else if ( itype == 1 ) {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = j; i <= m; i ++ ) {
          A[ ioffa + i - 1 + ( j - 1 ) * lda ] *= mul;
        }
      }
    } else if ( itype == 2 ) {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = 1; i <= Math.min( j, m ); i ++ ) {
          A[ ioffa + i - 1 + ( j - 1 ) * lda ] *= mul;
        }
      }
    } else if ( itype == 3 ) {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = 1; i <= Math.min( j + 1, m ); i ++ ) {
          A[ ioffa + i - 1 + ( j - 1 ) * lda ] *= mul;
        }
      }
    } else if ( itype == 4 ) {
      var k3 = kl + 1;
      var k4 = n + 1;
      for ( j = 1; j <= n; j ++ ) {
        for ( i = 1; i <= Math.min( k3, k4 - j ); i ++ ) {
          A[ ioffa + i - 1 + ( j - 1 ) * lda ] *= mul;
        }
      }
    } else if ( itype == 5 ) {
      var k1 = ku + 2;
      k3 = ku + 1;
      for ( j = 1; j <= n; j ++ ) {
        for ( i = Math.max( k1 - j, 1 ); i <= k3; i ++ ) {
          A[ ioffa + i - 1 + ( j - 1 ) * lda ] *= mul;
        }
      }
    } else if ( itype == 6 ) {
      k1 = kl + ku + 2;
      var k2 = kl + 1;
      k3 = 2 * kl + ku + 1;
      k4 = kl + ku + 1 + m;
      for ( j = 1; j <= n; j ++ ) {
        for ( i = Math.max( k1 - j, k2 ); i <= Math.min( k3, k4 - j );
        i ++ ) {
          A[ ioffa + i - 1 + ( j - 1 ) * lda ] *= mul;
        }
      }
    }
  }
}
LaPack1.zlascl = function( type, kl, ku, cfrom, cto, m, n, A,
lda, info ) {
  throw new Error("not programmed: complex array");
}
//************************************************************************
LaPack1.dlasd2 = function( nl, nr, sqre, k, d, z, alpha, beta,
U, ldu, Vt, ldvt, dsigma, U2, ldu2, Vt2, ldvt2, idxp, idx, idxc, idxq,
coltyp, info, ioffd, ioffz, ioffu, ioffvt, ioffdsigma, ioffu2, ioffvt2,
ioffidxp, ioffidx, ioffidxc, ioffidxq, ioffcoltyp) {
  throw new Error("not tested: complicated input");
  var ctot = new Array( 4 );
  var psm = new Array( 4 );
  info.setValue( 0 );
  if ( nl < 1 ) info.setValue( -1 );
  else if ( nr < 1 ) info.setValue( -2 );
  else if ( sqre != 1 && sqrt != 0 ) info.setValue( -3 );
  var n = nl + nr + 1;
  var m = n + sqre;
  if ( ldu < n ) info.setValue( -10 );
  else if ( ldvt < m ) info.setValue( -12 );
  else if ( ldu2 < n ) info.setValue( -15 );
  else if ( ldvt2 < m ) info.setValue( -17 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlasd2', - info.getValue() );
    return;
  }
  var nlp1 = nl + 1;
  var nlp2 = nl + 2;
  var z1 =
    alpha * vt[ ioffvt + nlp1 - 1 + ( nlp1 - 1 ) * ldvt ];
  z[ ioffz ] = z1;
  for ( var i = nl; i >= 1; i -- ) {
    z[ ioffz + i ] =
      alpha * Vt[ ioffvt + i - 1 + ( nlp1 - 1 ) * ldvt ];
    d[ ioffd + i ] = d[ ioffd + i - 1 ];
    idxq[ ioffidxq + i ] = idxq[ ioffidxq + i - 1 ] + 1;
  } // 10
  for ( i = nlp2; i <= m; i ++ ) {
    z[ ioffz + i - 1 ] =
      beta * Vt[ ioffvt + i - 1 + ( nlp2 - 1 ) * ioffvt ];
  } // 20
  for ( i = 2; i <= nlp1; i ++ ) coltyp[ ioffcoltyp + i - 1 ] = 1;
  for ( i = nlp2; i <= n; i ++ ) coltyp[ ioffcoltyp + i - 1 ] = 2;
  for ( i = nlp2; i <= n; i ++ ) idxq[ ioffidxq + i - 1 ] += nlp1;
  for ( i = 2; i <= n; i ++ ) {
    dsigma[ ioffdsigma + i - 1 ] =
      d[ ioffd + idxq[ ioffidxq + i - 1 ] - 1 ];
    U2[ ioffu2 + i - 1 ] =
      z[ ioffz + idxq[ ioffidxq + i - 1 ] - 1 ];
    idxc[ ioffidxc + i - 1 ] =
      coltyp[ ioffcoltyp + idxq[ ioffidxq + i - 1 ] - 1 ];
  } // 60
  LaPack0.dlamrg( nl, nr, dsigma, 1, 1, idx,
    ioffdsigma + 1, ioffidx + 1 );
  for ( i = 2; i <= n; i ++ ) {
    var idxi = 1 + idx[ ioffidx + i - 1 ];
    d[ ioffd + i - 1 ] = dsigma[ ioffdsigma + idxi - 1 ];
    z[ ioffz + i - 1 ] = U2[ ioffu2 + idxi - 1 ];
    coltyp[ ioffcoltyp + i - 1 ] = idxc[ ioffidxc + idxi - 1 ];
  } // 70
  var eps = LaPack0.dlamch( 'Epsilon' );
  var tol = Math.max( Math.abs( alpha ), Math.abs( beta ) );
  tol = 8. * eps * Math.max( Math.abs( d[ ioffd + n - 1 ] ), tol );
  k.setValue( 1 );
  var k2 = n + 1;
  var goto120 = false;
  for ( var j = 2; j <= n; j ++ ) {
    if ( Math.abs( z[ ioffz + j - 1 ] ) <= tol  ) {
      k2 --;
      idxp[ ioffidxp + k2 - 1 ] = j;
      coltyp[ ioffcoltyp + j - 1 ] = 4;
      if ( j == n ) {
        goto120 = true;
        break;
      }
    } else {
      var jprev = j;
      break;
    }
  } // 80
  if ( ! goto120 ) { // 90
    j = jprev;
    while ( true ) { // 100
      j ++;
      if ( j > n ) break;
      if ( Math.abs( z[ ioffz + j - 1 ] )  <= tol ) {
        k2 --;
        idxp[ ioffidxp + k2 - 1 ] = j;
        coltyp[ ioffcoltyp + j - 1 ] = 4;
      } else {
        if ( Math.abs( d[ ioffd + j - 1 ] - d[ ioffd + jprev - 1 ] )
        <= tol ) {
          var s = z[ ioffz + jprev - 1 ];
          var c = z[ ioffz + j - 1 ];
          var tau = LaPack0.dlapy2( c, s );
          c /= tau;
          s /= -tau;
          z[ ioffz + j - 1 ] = tau;
          z[ ioffz + jprev - 1 ] = 0.;
          var idxjp =
            idxq[ ioffidxq + idx[ ioffidx + jprev - 1 ] ];
          var idxj = idxq[ ioffidxq + idx[ ioffidx + j - 1 ] ];
          if ( idxjp <= nlp1 ) idxjp --;
          if ( idxj <= nlp1 ) idxj --;
          Blas1.drot( n, U, 1, U, 1, c, s, ioffu + ( idxjp - 1 ) * ldu,
            ioffu + ( idxj - 1 ) * ldu );
          Blas1.drot( m, Vt, ldvt, Vt, ldvt, c, s,  ioffvt + idxjp - 1,
            ioffvt + idxj - 1 );
          if ( coltyp[ ioffcoltyp + j - 1 ] !=
          coltyp[ ioffcoltyp + jprev - 1 ] ) {
            coltyp[ ioffcoltyp + j - 1 ] = 3;
          }
          coltyp[ ioffcoltyp + jprev - 1 ] = 4;
          k2 --;
          idxp[ ioffidxp + k2 - 1 ] = jprev;
          jprev = j;
        } else {
          k.setValue( k.getValue() + 1 );
          U2[ ioffu2 + k - 1 ] = z[ ioffz + jprev - 1 ];
          dsigma[ ioffdsigma + k - 1 ] = d[ ioffd + jprev - 1 ];
          idxp[ ioffidxp + k - 1 ] = jprev;
          jprev = j;
        }
      }
    } // 110
    k.setValue( k.getValue() + 1 );
    U2[ ioffu2 + k - 1 ] = z[ ioffz + jprev - 1 ];
    dsigma[ ioffdsigma + k - 1 ] = d[ ioffd + jprev - 1 ];
    idxp[ ioffidxp + k - 1 ] = jprev;
  } // 120
  for ( j = 1; j <= 4; j ++ ) ctot[ j - 1 ] = 0;
  for ( j = 2; j <= n; j ++ ) {
    var ct = coltyp[ ioffcoltyp + j - 1 ]
    ctot[ ct - 1 ] ++; 
  } // 140
  psm[ 0 ] = 2;
  psm[ 1 ] = 2 + ctot[ 0 ];
  psm[ 2 ] = psm[ 1 ] + ctot[ 1 ];
  psm[ 3 ] = psm[ 2 ] + ctot[ 2 ];
  for ( j = 2; j <= n; j ++ ) {
    var jp = idxp[ ioffidxp + j - 1 ];
    ct = coltyp[ ioffcoltyp + jp - 1 ];
    idxc[ ioffidxc + psm[ ct - 1 ] - 1 ] = j;
    psm[ ct - 1 ] ++;
  } // 150
  for ( j = 2; j <= n; j ++ ) {
    jp = idxp[ ioffidxp + j - 1 ];
    dsigma[ ioffdsigma + j - 1 ] = d[ ioffd + jp - 1 ];
    idxj = idxq[ ioffidxq + idxp[ ioffidxp + idxc[ ioffidxc + j - 1 ]
      - 1 ] ];
    if ( idxj <= nlp1 ) idxj --;
    Blas1.dcopy( n, U, 1, U2, 1, ioffu + ( idxj - 1 ) * ldu,
      ioffu2 + ( j - 1 ) * ioffu2 );
    Blas1.dcopy( m, Vt, ldvt, Vt2, ldvt2, ioffvt + idxj - 1,
      ioffvt2 + j - 1 );
  }
  dsigma[ ioffdsigma ] = 0.;
  var hlftol = tol / 2.;
  if ( Math.abs( dsigma[ ioffdsigma + 1 ] ) <= hlftol ) {
    dsigma[ ioffdsigma + 1 ] = hlftol;
  }
  if ( m > n ) {
    z[ ioffz ] = LaPack0.dlapy2( z1, z[ ioffz + m - 1 ] );
    if ( z[ ioffz ] <= tol ) {
      c = 1.;
      s = 0.;
      z[ ioffz ] = tol;
    } else {
      c = z1 / z[ ioffz ];
      s = z[ ioffz + m - 1 ] / z[ ioffz ];
    }
  } else {
    z[ ioffz ] = ( Math.abs( z1 ) <= tol ? tol : z1 );
  }
  Blas1.dcopy( k.getValue() - 1, U2, 1, z, 1, ioffu2 + 1, ioffz + 1 );
  LaPack0.dlaset( 'A', n, 1, 0., 0., U2, ldu2, ioffu2 );
  U2[ ioffu2 + nlp1 - 1 ] = 1.;
  if ( m > n ) {
    for ( i = 1; i <= nlp1; i ++ ) {
      Vt[ ioffvt + m - 1 + ( i - 1 ) * ldvt ] =
        - s * Vt[ ioffvt + nlp1 - 1 + ( i - 1 ) * ldvt ];
      Vt2[ ioffvt2 + ( i - 1 ) * ldvt2 ] =
        c * Vt[ ioffvt + nlp1 - 1 + ( i - 1 ) * ldvt ];
    } // 170
    for ( i = nlp2; i <= m; i ++ ) {
      Vt2[ ioffvt2 + ( i - 1 ) * ldvt2 ] =
        s * Vt[ ioffvt + m - 1 + ( i - 1 ) * ldvt ];
      Vt[ ioffvt + m - 1 + ( i - 1 ) * ldvt ] *= c;
    } // 180
  } else {
    Blas1.dcopy( m, Vt, ldvt, Vt2, ldvt2, ioffvt + nlp1 - 1, ioffvt2 );
  }
  if ( m > n ) {
    Blas1.dcopy( m, Vt, ldvt, Vt2, ldvt2, ioffvt + m - 1,
      ioffvt2 + m - 1 );
  }
  if ( n > k.getValue() ) {
    Blas1.dcopy( n - k.getValue(), dsigma, 1, d, 1, ioffdsigma + k.getValue(),
      ioffd + k.getValue() );
    LaPack0.dlacpy( 'A', n, n - k.getValue(), U2, ldu2, U, ldu,
      ioffu2 + k.getValue() * ldu2, ioffu + k.getValue() * ldu );
    LaPack0.dlacpy( 'A', n - k.getValue(), m, Vt2, ldvt2, Vt, ldvt,
      ioffvt2 + k.getValue(), ioffvt + k.getValue() );
  }
  for ( j = 1; j <= 4; j ++ ) {
    coltyp[ ioffcoltyp + j - 1 ] = ctot[ j - 1 ];
  }
}
//************************************************************************
LaPack1.dlasd7 = function( icompq, nl, nr, sqre, k, d, z, zw,
vf, vfw, vl, vlw, alpha, beta, dsigma, idx, idxp, idxq, perm, givptr,
Givcol, ldgcol, Givnum, ldgnum, cReference, sReference, info, ioffd,
ioffz, ioffzw, ioffvf, ioffvfw, ioffvl, ioffvlw, ioffdsigma, ioffidx,
ioffidxp, ioffidxq, ioffperm, ioffgivcol, ioffgivnum ) {
  throw new Error("not tested: complicated input");
  info.setValue( 0 );
  var n = nl + nr + 1;
  var m = n + sqre;
  if ( icompq < 0 || icompq > 1 ) info.setValue( -1 );
  else if ( nl < 1 ) info.setValue( -2 );
  else if ( nr < 1 ) info.setValue( -3 );
  else if ( sqre < 0 || sqre > 1 ) info.setValue( -4 );
  else if ( ldgcol < n ) info.setValue( -22 );
  else if ( ldgnum < n ) info.setValue( -24 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlasd7', - info.getValue() );
  }
  var nlp1 = nl + 1;
  var nlp2 = nl + 2;
  if ( icompq == 1 ) givptr.setValue( 0 );
  var z1 = alpha * vl[ ioffvl + nlp1 - 1 ];
  vl[ ioffvl + nlp1 - 1 ] = 0.;
  var tau = vf[ ioffvf + nlp1 - 1 ];
  for ( var i = nl; i >= 1; i -- ) {
    z[ ioffz + i ] = alpha * vl[ ioffvl + i - 1 ];
    vl[ ioffvl + i - 1 ] = 0.;
    vf[ ioffvf + i ] = vf[ ioffvf + i - 1 ];
    d[ ioffd + i ] = d[ ioffd + i - 1 ];
    idxq[ ioffidxq + i ] = idxq[ ioffidxq + i - 1 ] + 1;
  } // 10
  vf[ ioffvf ] = tau;
  for ( i = nlp2; i <= m; i ++ ) {
    z[ ioffz + i - 1 ] = beta * vf[ ioffvf + i - 1 ];
    vf[ ioffvf + i - 1 ] = 0.;
  } // 20
  for ( i = nlp2; i <= n; i ++ ) idxq[ ioffidxq + i - 1 ] += nlp1;
  for ( i = 2; i <= n; i ++ ) {
    var idxqi = idxq[ ioffidxq + i - 1 ];
    dsigma[ ioffdsigma + i - 1 ] = d[ ioffd + idxqi - 1 ];
    zw[ ioffzw + i - 1 ] = z[ ioffz + idxqi - 1 ];
    vfw[ ioffvfw + i - 1 ] = vf[ ioffvf + idxqi - 1 ];
    vlw[ ioffvlw + i - 1 ] = vl[ ioffvl + idxqi - 1 ];
  } // 40
  LaPack0.dlamrg( nl, nr, dsigma, 1, 1, idx, ioffdsigma + 1,
    ioffidx + 1 );
  for ( i = 2; i <= n; i ++ ) {
    var idxi = 1 + idx[ ioffidx + i - 1 ];
    d[ ioffd + i - 1 ] = dsigma[ ioffdsigma + idxi - 1 ];
    z[ ioffz + i - 1 ] = zw[ ioffzw + idxi - 1 ];
    vf[ ioffvf + i - 1 ] = vfw[ ioffvfw + idxi - 1 ];
    vl[ ioffvl + i - 1 ] = vlw[ ioffvlw + idxi - 1 ];
  } // 50
  var eps = LaPack0.dlamch( 'Epsilon' );
  var tol = Math.max( Math.abs( alpha ), Math.abs( beta ) );
  tol = 8. * 8. * eps
    * Math.max( Math.abs( d[ ioffd + n - 1 ] ), tol );
  k.setValue( 1 );
  var k2 = n + 1;
  var goto100 = false;
  for ( var j = 2; j <= n; j ++ ) {
    if ( Math.abs( z[ ioffz + j - 1 ] ) <= tol ) {
      k2 --;
      idxp[ ioffidxp + k2 - 1 ] = j;
      if ( j == n ) {
        goto100 = true;
        break;
      }
    } else {
      var jprev = j;
      break;
    }
  } // 60
  if ( ! goto100 ) { // 70
    j = jprev;
    while ( true ) { // 80
      j ++;
      if ( j > n ) break;
      if ( Math.abs( z[ ioffz + j - 1 ] ) <= tol ) {
        k2 --;
        idxp[ ioffidxp + k2 - 1 ] = j;
      } else {
        if ( Math.abs( d[ ioffd + j - 1 ] - d[ ioffd + jprev - 1 ] )
        <= tol ) {
          s.setValue( z[ ioffz + jprev - 1 ] );
          c.setValue( z[ ioffz + j - 1 ] );
          tau = LaPack0.dlapy2( c.getValue(), s.getValue() );
          z[ ioffz + j - 1 ] = tau;
          z[ ioffz + jprev - 1 ] = 0.;
          c.setValue( c.getValue() / tau );
          s.setValue( - s.getValue() / tau );
          if ( icompq == 1 ) {
            givptr.setValue( givptr.getValue() + 1 );
            var idxjp =
              idxq[ ioffidxq + idx[ ioffidx + jprev - 1 ] ];
            var idxj = idxq[ ioffidxq + idx[ ioffidx + j - 1 ] ];
            if ( idxjp <= nlp1 ) idxjp --;
            if ( idxj <= nlp1 ) idxj --;
            Givcol[ ioffgivcol + givptr - 1 + ldgcol ] = idxjp;
            Givcol[ ioffgivcol + givptr - 1 ] = idxj;
            Givnum[ ioffgivnum + givptr - 1 + ldgnum ] = c.getValue();
            Givnum[ ioffgivnum + givptr - 1 ] = s.getValue();
          }
          Blas1.drot( 1, vf, 1, vf, 1, c.getValue(), s.getValue(),
            ioffvf + jprev - 1, ioffvf + j - 1 );
          Blas1.drot( 1, vl, 1, vl, 1, c.getValue(), s.getValue(),
            ioffvl + jprev - 1, ioffvl + j - 1 );
          k2 --;
          idxp[ ioffidxp + k2 - 1 ] = jprev;
          jprev = j;
        } else {
          k.setValue( k.getValue() + 1 );
          zw[ ioffzw + k.getValue() - 1 ] = z[ ioffz + jprev - 1 ];
          dsigma[ ioffdsigma + k.getValue() - 1 ] = d[ ioffd + jprev - 1 ];
          idxp[ ioffidxp + k.getValue() - 1 ] = jprev;
          jprev = j;
        }
      }
    } // 90
    k.setValue( k.getValue() + 1 );
    zw[ ioffzw + k.getValue() - 1 ] = z[ ioffz + jprev - 1 ];
    dsigma[ ioffdsigma + k.getValue() - 1 ] = d[ ioffd + jprev - 1 ];
    idxp[ ioffidxp + k.getValue() - 1 ] = jprev;
  } // 100
  for ( j = 2; j <= n; j ++ ) {
    var jp = idxp[ ioffidxp + j - 1 ];
    dsigma[ ioffdsigma + j - 1 ] = d[ ioffd + jp - 1 ];
    vfw[ ioffvfw + j - 1 ] = vf[ ioffvf + jp - 1 ];
    vlw[ ioffvlw + j - 1 ] = vl[ ioffvl + jp - 1 ];
  } // 110
  if ( icompq == 1 ) {
    for ( j = 2; j <= n; j ++ ) {
      jp = idxp[ ioffidxp + j - 1 ];
      perm[ ioffperm + j - 1 ] =
        idxq[ ioffidxq + idx[ ioffidx + jp - 1 ] ];
      if ( perm[ ioffperm + j - 1 ] <= nlp1 ) {
        perm[ ioffperm + j - 1 ] --;
      }
    } // 120
  }
  Blas1.dcopy( n - k.getValue(), dsigma, 1, d, 1, ioffdsigma + k.getValue(),
    ioffd + k.getValue() );
  dsigma[ ioffdsigma ] = 0.;
  var hlftol = tol / 2.;
  if ( Math.abs( dsigma[ ioffdsigma + 1 ] ) <= hlftol ) {
    dsigma[ ioffdsigma + 1 ] = hlftol;
  }
  if ( m > n ) {
    z[ ioffz ] = LaPack0.dlapy2( z1, z[ ioffz + m - 1 ] );
    if ( z[ ioffz ] <= tol ) {
      c.setValue( 1. );
      s.setValue( 0. );
      z[ ioffz ] = tol;
    } else {
      c.setValue( z1 / z[ ioffz ] );
      s.setValue( - z[ ioffz + m - 1 ] / z[ ioffz ] );
    }
    Blas1.drot( 1, vf, 1, vf, 1, c.getValue(), s.getValue(), ioffvf + m - 1,
      ioffvf );
    Blas1.drot( 1, vl, 1, vl, 1, c.getValue(), s.getValue(), ioffvl + m - 1,
      ioffvl );
  } else z[ ioffz ] = ( Math.abs( z1 ) <= tol ? tol : z1 );
  Blas1.dcopy( k.getValue() - 1, zw, 1, z, 1, ioffzw + 1, ioffz + 1 );
  Blas1.dcopy( n - 1, vfw, 1, vf, 1, ioffvfw + 1, ioffvf + 1 );
  Blas1.dcopy( n - 1, vlw, 1, vl, 1, ioffvlw + 1, ioffvl + 1 );
}
//************************************************************************
LaPack1.dlasq6 = function( i0, n0, z, pp, dminReference,
dmin1Reference, dmin2Reference, dnReference, dnm1Reference,
dnm2Reference, ioffz ) {
  throw new Error("not tested: complicated input");
  if ( ( n0 - i0 - 1 ) <= 0 ) return;
  var safmin = LaPack0.dlamch( 'Safe minimum' );
  var j4 = 4 * i0 + pp - 3; 
  var emin = z[ ioffz + j4 + 3 ];
  var d = z[ ioffz + j4 - 1];
  dmin.setValue( d ); 
  if ( pp == 0 ) {
    for ( j4 = 4 * i0; j4 <= 4 * ( n0 - 3 ); j4 += 4 ) {
      z[ ioffz + j4 - 3 ] = d + z[ ioffz + j4 - 2 ]; 
      if ( z[ ioffz + j4 - 3 ] == 0. ) {
        z[ ioffz + j4 - 1 ] = 0.;
        d = z[ ioffz + j4 ];
        dmin.setValue( d ); 
        emin = 0.;
      } else if ( safmin * z[ ioffz + j4 ] < z[ ioffz + j4 - 3 ] &&
      safmin * z[ ioffz + j4 - 3 ] < z[ ioffz + j4 ] ) {
        var temp = z[ ioffz + j4 ] / z[ ioffz + j4 - 3 ];
        z[ ioffz + j4  - 1 ] = z[ ioffz + j4 - 2 ] * temp;
        d *= temp;
      } else { 
        z[ ioffz + j4  - 1 ] = z[ ioffz + j4 ]
          * ( z[ ioffz + j4 - 2 ] / z[ ioffz + j4 - 3 ] ) ;
        d = z[ ioffz + j4 ] * ( d / z[ ioffz + j4 - 3 ] ); 
      }
      dmin.setValue( Math.min( dmin.getValue(), d ) );
      emin = Math.min( emin, z[ ioffz + j4  - 1 ] ); 
    }
  } else {
    for ( j4 = 4 * i0; j4 <= 4 * ( n0 - 3 ); j4 += 4 ) {
      z[ ioffz + j4 - 4 ] = d + z[ ioffz + j4  - 1 ]; 
      if ( z[ ioffz + j4 - 4 ] == 0. ) {
        z[ ioffz + j4 - 2 ] = 0.;
        d = z[ ioffz + j4 + 1 ];
        dmin.setValue( d ); 
        emin = 0.;
      } else if( safmin * z[ ioffz + j4 + 1 ] < z[ ioffz + j4 - 4 ] 
      && safmin * z[ ioffz + j4 - 4 ] < z[ ioffz + j4 + 1 ] ) {
        temp = z[ ioffz + j4 + 1 ] / z[ ioffz + j4 - 4 ];
        z[ ioffz + j4 - 2 ] = z[ ioffz + j4  - 1 ] * temp;
        d *= temp;
      } else { 
        z[ ioffz + j4 - 2 ] = z[ ioffz + j4 + 1 ]
          * ( z[ ioffz + j4  - 1 ] / z[ ioffz + j4 - 4 ] ); 
        d = z[ ioffz + j4 + 1 ] * ( d / z[ ioffz + j4 - 4 ] ); 
      }
      dmin.setValue( Math.min( dmin.getValue(), d ) );
      emin = Math.min( emin, z[ ioffz + j4 - 2 ] ); 
    }
  }
  dnm2.setValue( d );
  dmin2.setValue( dmin.getValue() );
  j4 = 4 * ( n0 - 2 ) - pp;
  var j4p2 = j4 + 2 * pp - 1;
  z[ ioffz + j4 - 3 ] = dnm2.getValue() + z[ ioffz + j4p2  - 1 ];
  if ( z[ ioffz + j4 - 3 ] == 0. ) {
    z[ ioffz + j4  - 1 ] = 0.;
    dnm1.setValue( z[ ioffz + j4p2 + 1 ] );
    dmin.setValue( dnm1.getValue() );
    emin = 0.;
  } else if ( safmin * z[ ioffz + j4p2 + 1 ] < z[ ioffz + j4 - 3 ]  &&
  safmin * z[ ioffz + j4 - 3 ] < z[ ioffz + j4p2 + 1 ] ) {
    temp = z[ ioffz + j4p2 + 1 ] / z[ ioffz + j4 - 3 ]
    z[ ioffz + j4  - 1 ] = z[ ioffz + j4p2  - 1 ] * temp;
    dnm1.setValue( dnm2.getValue()  * temp );
  } else {
     z[ ioffz + j4  - 1 ] = z[ ioffz + j4p2 + 1 ]
       * ( z[ ioffz + j4p2  - 1 ] / z[ ioffz + j4 - 3 ] );
     dnm1.setValue( z[ ioffz + j4p2 + 1 ]
       * ( dnm2.getValue() / z[ ioffz + j4 - 3 ] ) );
  }
  dmin.setValue( Math.min( dmin.getValue(), dnm1.getValue() ) );
  dmin1.setValue( dmin.getValue() );
  j4 += 4;
  j4p2 = j4 + 2 * pp - 1;
  z[ ioffz + j4 - 3 ] = dnm1.getValue() + z[ ioffz + j4p2  - 1 ];
  if( z[ ioffz + j4 - 3 ] == 0. ) {
    z[ ioffz + j4  - 1 ] = 0.;
    dn.setValue( z[ ioffz + j4p2 + 1 ] );
    dmin.setValue( dn.getValue() );
    emin = 0.;
  } else if( safmin * z[ ioffz + j4p2 + 1 ] < z[ ioffz + j4 - 3 ]  &&
  safmin * z[ ioffz + j4 - 3 ] < z[ ioffz + j4p2 + 1 ] ) {
    temp = z[ ioffz + j4p2 + 1 ] / z[ ioffz + j4 - 3 ];
    z[ ioffz + j4  - 1 ] = z[ ioffz + j4p2  - 1 ] * temp;
    dn.setValue( dnm1.getValue() * temp );
  } else {
     z[ ioffz + j4  - 1 ] = z[ ioffz + j4p2 + 1 ]
       * ( z[ ioffz + j4p2  - 1 ] / z[ ioffz + j4 - 3 ] );
     dn.setValue( z[ ioffz + j4p2 + 1 ]
       * ( dnm1.getValue() / z[ ioffz + j4 - 3 ] ) );
  }
  dmin.setValue( Math.min( dmin.getValue(), dn.getValue() ) );
  z[ ioffz + j4 + 1 ] = dn.getValue();
  z[ ioffz + 4 * n0 - pp  - 1 ] = emin;
}
//************************************************************************
LaPack1.dlasv2 = function( f, g, h, ssminReference,
ssmaxReference, snrReference, csrReference, snlReference, cslReference) {
  var ft = f;
  var fa = Math.abs( ft );
  var ht = h;
  var ha = Math.abs( h );
  var pmax = 1;
  var swap = ( ha > fa );
  if ( swap ) {
    pmax = 3;
    var temp = ft;
    ft = ht;
    ht = temp;
    temp = fa;
    fa = ha;
    ha = temp;
  }
  var gt = g;
  var ga = Math.abs( gt );
  if ( ga == 0. ) {
    ssminReference.setValue( ha );
    ssmaxReference.setValue( fa );
    var clt = 1.;
    var crt = 1.;
    var slt = 0.;
    var srt = 0.;
  } else {
    var gasmal = true;
    if ( ga > fa ) {
      pmax = 2;
      if ( ( fa / ga ) < LaPack0.dlamch( 'Eps' ) ) {
        gasmal = false;
        ssmaxReference.setValue( ga );
        ssminReference.setValue(
          ( ha > 1. ? fa / ( ga / ha ) : ( fa / ga ) * ha ) );
        clt = 1.;
        slt = ht / gt;
        srt = 1.;
        crt = ft / gt;
      }
    }
    if ( gasmal ) {
      var d = fa - ha;
      var l = ( d == fa ? 1. : d / fa );
      var m = gt / ft;
      var t = 2. - l;
      var mm = m * m;
      var tt = t * t;
      var s = Math.sqrt( tt + mm );
      var r =
        ( l == 0. ? Math.abs( m ) : Math.sqrt( l * l + mm ) );
      var a = 0.5 * ( s + r );
      ssminReference.setValue( ha / a );
      ssmaxReference.setValue( fa * a );
      if ( mm == 0. ) {
        t = ( l == 0. ?
          ( ft >= 0. ? 2. : -2. ) * ( gt >= 0. ? 1. : -1. ) :
          gt / ( ft >= 0. ? Math.abs( d ) : - Math.abs( d ) )
          + m / t );
      } else {
        t = ( m / ( s + t ) + m / ( r + l ) ) * ( 1. + a );
      }
      l = Math.sqrt( t * t + 4. );
      crt = 2. / l;
      srt = t / l;
      clt = ( crt + srt * m ) / a;
      slt = ( ht / ft ) * srt / a;
    }
  }
  if ( swap ) {
    cslReference.setValue( srt );
    snlReference.setValue( crt );
    csrReference.setValue( slt );
    snrReference.setValue( clt );
  } else {
    cslReference.setValue( clt );
    snlReference.setValue( slt );
    csrReference.setValue( crt );
    snrReference.setValue( srt );
  }
  if ( pmax == 1 ) {
    var tsign = ( csrReference.getValue() >= 0. ? 1. : -1. )
      * ( cslReference.getValue() >= 0. ? 1. : -1. )
      * ( f >= 0. ? 1. : -1. );
  }
  if ( pmax == 2 ) {
    tsign = ( snrReference.getValue() >= 0. ? 1. : -1. )
      * ( cslReference.getValue() >= 0. ? 1. : -1. )
      * ( g >= 0. ? 1. : -1. );
  }
  if ( pmax == 3 ) {
    tsign = ( snrReference.getValue() >= 0. ? 1. : -1. )
      * ( snlReference.getValue() >= 0. ? 1. : -1. )
      * ( h >= 0. ? 1. : -1. );
  }
  ssmaxReference.setValue( ( tsign >= 0. ?
    Math.abs( ssmaxReference.getValue() ) :
    - Math.abs( ssmaxReference.getValue() ) ) );
  ssminReference.setValue( ( tsign * ( f >= 0. ? 1. : -1. )
    * ( h>= 0. ? 1. : -1. ) >= 0. ?
    Math.abs( ssminReference.getValue() ) :
    - Math.abs( ssminReference.getValue() ) ) );
}
//************************************************************************
LaPack1.dlasy2 = function( ltranl, ltranr, isgn, n1, n2, Tl,
ldtl, Tr, ldtr, B, ldb, scaleReference, X, ldx, xnormReference, info,
iofftl, iofftr, ioffb, ioffx ) {
  var bswpiv = new Array( false, true, false, true );
  var xswpiv = new Array( false, false, true, true );
  var jpiv = new Array(4);
  var locl21 = new Array( 2, 1, 4, 3 );
  var locu12 = new Array( 3, 4, 1, 2 );
  var locu22 = new Array( 4, 3, 2, 1 );
  var btmp = new Array(4);
  var t16 = new Array(4*4);
  var tmp = new Array(4);
  var x2 = new Array(2);

  info.setValue( 0 );
  if ( n1 == 0 || n2 == 0 ) return;
  var eps = LaPack0.dlamch( 'P' );
  var smlnum = LaPack0.dlamch( 'S' ) / eps;
  var sgn = isgn;
  var k = n1 + n1 + n2 - 2
  if ( k <= 1 || k > 4 ) { // 1 x 1
    var tau1 = Tl[ iofftl ] + sgn * Tr[ iofftr ];
    var bet = Math.abs( tau1 );
    if ( bet <= smlnum ) {
      tau1 = smlnum;
      bet = smlnum;
      info.setValue( 1 );
    }
    scaleReference.setValue( 1. );
    var gam = Math.abs( B[ ioffb ] );
    if ( smlnum * gam > bet ) scaleReference.setValue( 1. / gam );
    X[ ioffx ] = ( B[ ioffb ] * scaleReference.getValue() ) / tau1;
    xnormReference.setValue( Math.abs( X[ ioffx ] ) );
    return;
  }
  if ( k == 2 ) { // 1 x 2
    var trmax = Math.max(
      Math.max( Math.abs( Tr[ iofftr ] ),
                Math.abs( Tr[ iofftr + ldtr ] ) ),
      Math.max( Math.abs( Tr[ iofftr + 1 ] ),
                Math.abs( Tr[ iofftr + 1 + ldtr ] ) ) );
    var smin = Math.max( eps * Math.max(
      Math.abs( Tl[ iofftl ] ), trmax ), smlnum );
    tmp[ 0 ] = Tl[ iofftl ] + sgn * Tr[ iofftr ];
    tmp[ 3 ] = Tl[ iofftl ] + sgn * Tr[ iofftr + 1 + ldtr ];
    if ( ltranr ) {
      tmp[ 1 ] = sgn * Tr[ iofftr + 1 ];
      tmp[ 2 ] = sgn * Tr[ iofftr + ldtr ];
    } else {
      tmp[ 1 ] = sgn * Tr[ iofftr + ldtr ];
      tmp[ 2 ] = sgn * Tr[ iofftr + 1 ];
    }
    btmp[ 0 ] = B[ ioffb ];
    btmp[ 1 ] = B[ ioffb + ldb ];
  }
  if ( k == 3 ) { // 2 x 1
    var tlmax = Math.max(
      Math.max( Math.abs( Tl[ iofftl ] ),
                Math.abs( Tl[ iofftl + ldtl ] ) ),
      Math.max( Math.abs( Tl[ iofftl + 1 ] ),
                Math.abs( Tl[ iofftl + 1 + ldtl ] ) ) );
    smin = Math.max( eps * Math.max(
      Math.abs( Tr[ iofftr ] ), tlmax ), smlnum );
    tmp[ 0 ] = Tl[ iofftl ] + sgn * Tr[ iofftr ];
    tmp[ 3 ] = Tl[ iofftl + 1 + ldtl ] + sgn * Tr[ iofftr ];
    if ( ltranl ) {
      tmp[ 1 ] = Tl[ iofftl + ldtl ];
      tmp[ 2 ] = Tl[ iofftl + 1 ];
    } else {
      tmp[ 1 ] = Tl[ iofftl + 1 ];
      tmp[ 2 ] = Tl[ iofftl + ldtl ];
    }
    btmp[ 0 ] = B[ ioffb ];
    btmp[ 1 ] = B[ ioffb + 1 ];
  }
  if ( k == 2 || k == 3 ) {
    var ipiv = Blas1.idamax( 4, tmp, 1, 0 );
    var u11 = tmp[ ipiv - 1 ];
    if ( Math.abs( u11 ) <= smin ) {
      info.setValue( 1 );
      u11 = smin;
    }
    var u12 = tmp[ locu12[ ipiv - 1 ] - 1 ];
    var l21 = tmp[ locl21[ ipiv - 1 ] - 1 ] / u11;
    var u22 = tmp[ locu22[ ipiv - 1 ] - 1 ] - u12 * l21;
    var xswap = xswpiv[ ipiv - 1 ];
    var bswap = bswpiv[ ipiv - 1 ];
    if ( Math.abs( u22 ) <= smin ) {
      info.setValue( 1 );
      u22 = smin;
    }
    if ( bswap ) {
      var temp = btmp[ 1 ];
      btmp[ 1 ] = btmp[ 0 ] - l21 * temp;
      btmp[ 0 ] = temp;
    } else btmp[ 1 ] -= l21 * btmp[ 0 ];
    scaleReference.setValue( 1. );
    if ( ( 2. * smlnum ) * Math.abs( btmp[ 1 ] ) > Math.abs( u22 ) ||
    ( 2. * smlnum ) * Math.abs( btmp[ 0 ] ) > Math.abs( u11 ) ) {
      scaleReference.setValue( 0.5
        / Math.max( Math.abs( btmp[ 0 ] ), Math.abs( btmp[ 1 ] ) ) );
      btmp[ 0 ] *= scaleReference.getValue();
      btmp[ 1 ] *= scaleReference.getValue();
    }
    x2[ 1 ] = btmp[ 1 ] / u22;
    x2[ 0 ] = btmp[ 0 ] / u11 - ( u12 / u11 ) * x2[ 1 ];
    if ( xswap ) {
      temp = x2[ 1 ];
      x2[ 1 ] = x2[ 0 ];
      x2[ 0 ] = temp;
    }
    X[ ioffx ] = x2[ 0 ];
    if ( n1 == 1 ) {
      X[ ioffx + ldx ] = x2[ 1 ];
      xnormReference.setValue( Math.abs( X[ ioffx ] )
                  + Math.abs( X[ ioffx + ldx ] ) );
    } else {
      X[ ioffx + 1 ] = x2[ 1 ];
      xnormReference.setValue( Math.max( Math.abs( X[ ioffx ] ),
                              Math.abs( X[ ioffx + 1 ] ) ) );
    }
    return;
  }
  trmax = Math.max( Math.max( Math.abs( Tr[ iofftr ] ),
                              Math.abs( Tr[ iofftr + ldtr ] ) ),
                    Math.max( Math.abs( Tr[ iofftr + 1 ] ),
                              Math.abs( Tr[ iofftr + 1 + ldtr ] ) ) );
  tlmax = Math.max( Math.max( Math.abs( Tl[ iofftl ] ),
                              Math.abs( Tl[ iofftl + ldtl ] ) ),
                    Math.max( Math.abs( Tl[ iofftl + 1 ] ),
                              Math.abs( Tl[ iofftl + 1 + ldtl ] ) ) );
  smin = Math.max( eps * Math.max( trmax, tlmax ), smlnum );
  btmp[ 0 ] = 0.;
  Blas1.dcopy( 16, btmp, 0, t16, 1, 0, 0 );
  t16[ 0 + 0 * 4 ] = Tl[ iofftl ] + sgn * Tr[ iofftr ];
  t16[ 1 + 1 * 4 ] = Tl[ iofftl + 1 + ldtl ] + sgn * Tr[ iofftr ];
  t16[ 2 + 2 * 4 ] = Tl[ iofftl ] + sgn * Tr[ iofftr + 1 + ldtr ];
  t16[ 3 + 3 * 4 ] = Tl[ iofftl + 1 + ldtl ]
    + sgn * Tr[ iofftr + 1 + ldtr ];
  if ( ltranl ) {
    t16[ 0 + 1 * 4 ] = Tl[ iofftl + 1 ];
    t16[ 1 + 0 * 4 ] = Tl[ iofftl + ldtl ];
    t16[ 2 + 3 * 4 ] = Tl[ iofftl + 1 ];
    t16[ 3 + 2 * 4 ] = Tl[ iofftl + ldtl ];
  } else {
    t16[ 0 + 1 * 4 ] = Tl[ iofftl + ldtl ];
    t16[ 1 + 0 * 4 ] = Tl[ iofftl + 1 ];
    t16[ 2 + 3 * 4 ] = Tl[ iofftl + ldtl ];
    t16[ 3 + 2 * 4 ] = Tl[ iofftl + 1 ];
  }
  if ( ltranr ) {
    t16[ 0 + 2 * 4 ] = sgn * Tr[ iofftr + ldtr ];
    t16[ 1 + 3 * 4 ] = sgn * Tr[ iofftr + ldtr ];
    t16[ 2 + 0 * 4 ] = sgn * Tr[ iofftr + 1 ];
    t16[ 3 + 1 * 4 ] = sgn * Tr[ iofftr + 1 ];
  } else {
    t16[ 0 + 2 * 4 ] = sgn * Tr[ iofftr + 1 ];
    t16[ 1 + 3 * 4 ] = sgn * Tr[ iofftr + 1 ];
    t16[ 2 + 0 * 4 ] = sgn * Tr[ iofftr + ldtr ];
    t16[ 3 + 1 * 4 ] = sgn * Tr[ iofftr + ldtr ];
  }
  btmp[ 0 ] = B[ ioffb ];
  btmp[ 1 ] = B[ ioffb + 1 ];
  btmp[ 2 ] = B[ ioffb + ldb ];
  btmp[ 3 ] = B[ ioffb + 1 + ldb ];
  for ( var i = 1; i <= 3; i ++ ) {
    var xmax = 0.;
    for ( var ip = i; ip <= 4; ip ++ ) {
      for ( var jp = i; jp <= 4; jp ++ ) {
        if ( Math.abs( t16[ ip - 1 + ( jp - 1 ) * 4 ] ) >= xmax ) {
          xmax = Math.abs( t16[ ip - 1 + ( jp - 1 ) * 4 ] );
          var ipsv = ip;
          var jpsv = jp;
        }
      }
    }
    if ( ipsv != i ) {
      Blas1.dswap( 4, t16, 4, t16, 4, ipsv - 1, i - 1 );
      temp = btmp[ i - 1 ];
      btmp[ i - 1 ] = btmp[ ipsv - 1 ];
      btmp[ ipsv - 1 ] = temp;
    }
    if ( jpsv != i ) {
      Blas1.dswap( 4, t16, 1, t16, 1, ( jpsv - 1 ) * 4,
        ( i - 1 ) * 4 );
    }
    jpiv[ i - 1 ] = jpsv;
    if ( Math.abs( t16[ i - 1 + ( i - 1 ) * 4 ] ) < smin ) {
      info.setValue( 1 );
      t16[ i - 1 + ( i - 1 ) * 4 ] = smin;
    }
    for ( var j = i + 1; j <= 4; j ++ ) {
      t16[ j - 1 + ( i - 1 ) * 4 ] /= t16[ i - 1 + ( i - 1 ) * 4 ];
      btmp[ j - 1 ] -=
        t16[ j - 1 + ( i - 1 ) * 4 ] * btmp[ i - 1 ];
      for ( k = i + 1; k <= 4; k ++ ) {
        t16[ j - 1 + ( k - 1 ) * 4 ] -= t16[ j - 1 + ( i - 1 ) * 4 ]
          * t16[ i - 1 + ( k - 1 ) * 4 ];
      }
    }
  }
  if ( Math.abs( t16[ 3 + 3 * 4 ] ) < smin ) t16[ 3 + 3 * 4 ] = smin;
  scaleReference.setValue( 1. );
  if ( ( 8. * smlnum ) * Math.abs( btmp[ 0 ] ) > Math.abs( t16[ 0 ] ) 
  || ( 8. * smlnum ) * Math.abs( btmp[ 1 ] ) > Math.abs( t16[ 5 ] )
  || ( 8. * smlnum ) * Math.abs( btmp[ 2 ] ) > Math.abs( t16[ 10 ] )
  || ( 8. * smlnum ) * Math.abs( btmp[ 3 ] ) > Math.abs( t16[ 15 ] ) )
  {
    scaleReference.setValue( ( 1. / 8. ) / Math.max( 
      Math.max( Math.abs( btmp[ 0 ] ), Math.abs( btmp[ 1 ] ) ),
      Math.max( Math.abs( btmp[ 2 ] ), Math.abs( btmp[ 3 ] ) ) ) );
    btmp[ 0 ] *= scaleReference.getValue();
    btmp[ 1 ] *= scaleReference.getValue();
    btmp[ 2 ] *= scaleReference.getValue();
    btmp[ 3 ] *= scaleReference.getValue();
  }
  for ( i = 1; i <= 4; i ++ ) {
    k = 5 - i;
    temp = 1. / t16[ k - 1 + ( k - 1 ) * 4 ];
    tmp[ k - 1 ] = btmp[ k - 1 ] * temp;
    for ( j = k + 1; j <= 4; j ++ ) {
      tmp[ k - 1 ] -=
        ( temp * t16[ k - 1 + ( j - 1 ) * 4 ] ) * tmp[ j - 1 ];
    }
  }
  for ( i = 1; i <= 3; i ++ ) {
    if ( jpiv[ 3 - i ] != 4 - i ) {
      temp = tmp[ 3 - i ];
      tmp[ 3 - i ] = tmp[ jpiv[ 3 - i ] - 1 ];
      tmp[ jpiv[ 3 - i ] - 1 ] = temp;
    }
  }
  X[ ioffx ] = tmp[ 0 ];
  X[ ioffx + 1 ] = tmp[ 1 ];
  X[ ioffx + ldx ] = tmp[ 2 ];
  X[ ioffx + 1 + ldx ] = tmp[ 3 ];
  xnormReference.setValue( Math.max( 
    Math.abs( tmp[ 0 ] ) + Math.abs( tmp[ 2 ] ),
    Math.abs( tmp[ 1 ] ) + Math.abs( tmp[ 3 ] ) ) );
}
//************************************************************************
LaPack1.dlatbs = function( uplo, trans, diag, normin, n, kd,
AB, ldab, x, scaleReference, cnorm, info ) {
  throw new Error("not programmed: band matrix");
}
LaPack1.zlatbs = function( uplo, trans, diag, normin, n, kd,
AB, ldab, x, scale, cnorm, info ) {
  throw new Error("not programmed: complex band matrix");
}
//************************************************************************
LaPack1.dlatps = function( uplo, trans, diag, normin, n, AP, x,
scaleReference, cnorm, info ) {
  throw new Error("not programmed: packed matrix");
}
LaPack1.zlatps = function( uplo, trans, diag, normin, n, AP, x,
scaleReference, cnorm, info ) {
  throw new Error("not programmed: complex packed matrix");
}
//************************************************************************
LaPack1.dlatrs = function( uplo, trans, diag, normin, n, A,
lda, x, scaleReference, cnorm, info, ioffa, ioffx, ioffcnorm) {
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var notran = ( trans.charAt(0).toUpperCase() == 'N' );
  var nounit = ( diag.charAt(0).toUpperCase() == 'N' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( ! notran && trans.charAt(0).toUpperCase() != 'T' &&
  trans.charAt(0).toUpperCase() != 'C' ) {
    info.setValue( -2 );
  } else if ( ! nounit && diag.charAt(0).toUpperCase() != 'U' ) {
    info.setValue( -3 );
  } else if ( normin.charAt(0).toUpperCase() != 'Y' &&
  normin.charAt(0).toUpperCase() != 'N' ) {
    info.setValue( -4 );
  } else if ( n < 0 ) info.setValue( -5 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -7 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlatrs', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  var smlnum = LaPack0.dlamch( 'Safe minimum' )
    / LaPack0.dlamch( 'Precision' );
  var bignum = 1. / smlnum;
  scaleReference.setValue( 1. );
  if ( normin.charAt(0).toUpperCase() == 'N' ) {
    if ( upper ) {
      for ( var j = 1; j <= n; j ++ ) {
        cnorm[ ioffcnorm + j - 1 ] =
          Blas1.dasum( j - 1, A, 1, ioffa + ( j - 1 ) * lda );
      }
    } else {
      for ( j = 1; j <= n - 1; j ++ ) {
        cnorm[ ioffcnorm + j - 1 ] =
          Blas1.dasum( n - j, A, 1, ioffa + j + ( j - 1 ) * lda );
      }
      cnorm[ ioffcnorm + n - 1 ] = 0.;
    }
  }
  var imax = Blas1.idamax( n, cnorm, 1, ioffcnorm );
  var tmax = cnorm[ ioffcnorm + imax - 1 ];
  if ( tmax <= bignum ) var tscal = 1.;
  else {
    tscal = 1. / ( smlnum * tmax );
    Blas1.dscal( n, tscal, cnorm, 1, ioffcnorm );
  }
  j = Blas1.idamax( n, x, 1, ioffx );
  var xmax = Math.abs( x[ ioffx + j - 1 ] );
  var xbnd = xmax;
  if ( notran ) {
    if ( upper ) {
      var jfirst = n;
      var jlast = 1;
      var jinc = -1;
    } else {
      jfirst = 1;
      jlast = n;
      jinc = 1;
    }
    if ( tscal != 1. ) {
      var grow = 0.;
    } else {
      if ( nounit ) {
        grow = 1. / Math.max( xbnd, smlnum );
        xbnd = grow;
        var goto50 = false;
        for ( j = jfirst; ( jinc == 1 && j <= jlast ) ||
        ( jinc == -1 && j >= jlast ); j += jinc ) {
          if ( grow <= smlnum ) {
            goto50 = true;
            break;
          }
          var tjj =
            Math.abs( A[ ioffa + j - 1 + ( j - 1 ) * lda ] );
          xbnd = Math.min( xbnd, Math.min( 1., tjj ) * grow );
          if ( tjj + cnorm[ ioffcnorm + j - 1 ] >= smlnum ) {
            grow *= tjj / ( tjj + cnorm[ ioffcnorm + j - 1 ] );
          } else grow = 0.;
        } // 30
        if ( ! goto50 ) grow = xbnd;
      } else {
        grow = Math.min( 1., 1. / Math.max( xbnd, smlnum ) );
        for ( j = jfirst; ( jinc == 1 && j <= jlast ) ||
        ( jinc == -1 && j >= jlast ); j += jinc ) {
          if ( grow <= smlnum ) break;
          grow *= 1. / ( 1. + cnorm[ ioffcnorm + j - 1 ] );
        } // 40
      }
    } // 50
  } else {
    if ( upper ) {
      jfirst = 1;
      jlast = n;
      jinc = 1;
    } else {
      jfirst = n;
      jlast = 1;
      jinc = -1;
    }
    if ( tscal != 1. ) {
      grow = 0.;
    } else {
      if ( nounit ) {
        grow = 1. / Math.max( xbnd, smlnum );
        xbnd = grow;
        var goto80 = false;
        for ( j = jfirst; ( jinc == 1 && j <= jlast ) ||
        ( jinc == -1 && j >= jlast ); j += jinc ) {
          if ( grow <= smlnum ) {
            goto80 = true;
            break;
          }
          var xj = 1. + cnorm[ ioffcnorm + j - 1 ];
          grow = Math.min( grow, xbnd / xj );
          tjj = Math.abs( A[ ioffa + j - 1 + ( j - 1 ) * lda ] );
          if ( xj > tjj ) xbnd *= tjj / xj;
        } // 60
        if ( ! goto80 ) grow = Math.min( grow, xbnd );
      } else {
        grow = Math.min( 1., 1. / Math.max( xbnd, smlnum ) );
        for ( j = jfirst; ( jinc == 1 && j <= jlast ) ||
        ( jinc == -1 && j >= jlast ); j += jinc ) {
          if ( grow <= smlnum ) break;
          xj = 1. + cnorm[ ioffcnorm + j - 1 ];
          grow /= xj;
        } // 70
      }
    } // 80
  }
  if ( ( grow * tscal ) > smlnum ) {
    Blas2.dtrsv( uplo, trans, diag, n, A, lda, x, 1, ioffa, ioffx );
  } else {
    if ( xmax > bignum ) {
      scaleReference.setValue( bignum / xmax );
      Blas1.dscal( n, scaleReference.getValue(), x, 1, ioffx );
      xmax = bignum;
    }
    if ( notran ) {
      for ( j = jfirst; ( jinc == 1 && j <= jlast ) ||
      ( jinc == -1 && j >= jlast ); j += jinc ) {
        xj = Math.abs( x[ ioffx + j - 1 ] );
        var goto100 = false;
        if ( nounit ) {
          var tjjs =
            A[ ioffa + j - 1 + ( j - 1 ) * lda ] * tscal;
        } else {
          tjjs = tscal;
          if ( tscal == 1. ) goto100 = true;
        }
        if ( ! goto100 ) {
          tjj = Math.abs( tjjs );
          if ( tjj > smlnum ) {
            if ( tjj < 1. ) {
              if ( xj > tjj * bignum ) {
                var rec = 1. / xj;
                Blas1.dscal( n, rec, x, 1, ioffx );
                scaleReference.setValue( scaleReference.getValue() * rec );
                xmax *= rec;
              }
            }
            x[ ioffx + j - 1 ] /= tjjs;
            xj = Math.abs( x[ ioffx + j - 1 ] );
          } else if ( tjj > 0. ) {
            if ( xj > tjj * bignum ) {
              rec = ( tjj * bignum ) / xj;
              if ( cnorm[ ioffcnorm + j - 1 ] > 1. ) {
                rec /= cnorm[ ioffcnorm + j - 1 ];
              }
              Blas1.dscal( n, rec, x, 1, ioffx );
              scaleReference.setValue( scaleReference.getValue() * rec );
              xmax *= rec;
            }
            x[ ioffx + j - 1 ] /= tjjs;
            xj = Math.abs( x[ ioffx + j - 1 ] );
          } else {
            for ( var i = 1; i <= n; i ++ ) {
              x[ ioffx + i - 1 ] = 0.;
            }
            x[ ioffx + j - 1 ] = 1.;
            xj = 1.;
            scaleReference.setValue( 0. );
            xmax = 0.;
          }
        } // 100
        if ( xj > 1. ) {
          rec = 1. / xj;
          if ( cnorm[ ioffcnorm + j - 1 ] > ( bignum - xmax ) * rec ) {
            rec *= 0.5;
            Blas1.dscal( n, rec, x, 1, ioffx );
            scaleReference.setValue( scaleReference.getValue() * rec );
          }
        } else if ( xj * cnorm[ ioffcnorm + j - 1 ]
        > ( bignum - xmax ) ) {
          Blas1.dscal( n, 0.5, x, 1, ioffx );
          scaleReference.setValue( scaleReference.getValue() * 0.5 );
        }
        if ( upper ) {
          if ( j > 1 ) {
            Blas1.daxpy( j - 1, - x[ ioffx + j - 1 ] * tscal, A, 1,
              x, 1, ioffa + ( j - 1 ) * lda, ioffx );
            i = Blas1.idamax( j - 1, x, 1, ioffx );
            xmax = Math.abs( x[ ioffx + i - 1 ] );
          }
        } else {
          if ( j < n ) {
            Blas1.daxpy( n - j, - x[ ioffx + j - 1 ] * tscal, A, 1,
              x, 1, ioffa + j + ( j - 1 ) * lda, ioffx + j );
            i = j + Blas1.idamax( n - j, x, 1, ioffx + j );
            xmax = Math.abs( x[ ioffx + i - 1 ] );
          }
        }
      }
    } else {
      for ( j = jfirst; ( jinc == 1 && j <= jlast ) ||
      ( jinc == -1 && j >= jlast ); j += jinc ) {
        xj = Math.abs( x[ ioffx + j - 1 ] );
        var uscal = tscal;
        rec = 1. / Math.max( xmax, 1. );
        if ( cnorm[ ioffcnorm + j - 1 ] > ( bignum - xj ) * rec ) {
          rec *= 0.5;
          if ( nounit ) {
            tjjs = A[ ioffa + j - 1 + ( j - 1 ) * lda ] * tscal;
          } else tjjs = tscal;
          tjj = Math.abs( tjjs );
          if ( tjj > 1. ) {
            rec = Math.min( 1., rec * tjj );
            uscal /= tjjs;
          }
          if ( rec < 1. ) {
            Blas1.dscal( n, rec, x, 1, ioffx );
            scaleReference.setValue( scaleReference.getValue() * rec );
            xmax *= rec;
          }
        }
        var sumj = 0.;
        if ( uscal == 1. ) {
          if ( upper ) {
            sumj = Blas1.ddot( j - 1, A, 1, x, 1,
              ioffa + ( j - 1 ) * lda, ioffx );
          } else if ( j < n ) {
            sumj = Blas1.ddot( n - j, A, 1, x, 1,
              ioffa + j + ( j - 1 ) * lda, ioffx + j );
          }
        } else {
          if ( upper ) {
            for ( i = 1; i <= j - 1; i ++ ) {
              sumj += ( A[ ioffa + i - 1 + ( j - 1 ) * lda ] * uscal )
                * x[ ioffx + i - 1 ];
            }
          } else if ( j < n ) {
            for ( i = j + 1; i <= n; i ++ ) {
              sumj += ( A[ ioffa + i - 1 + ( j - 1 ) * lda ] * uscal )
                * x[ ioffx + i - 1 ];
            }
          }
        }
        if ( uscal == tscal ) {
          x[ ioffx + j - 1 ] -= sumj;
          xj = Math.abs( x[ ioffx + j - 1 ] );
          var goto150 = false;
          if ( nounit ) {
            tjjs = A[ ioffa + j - 1 + ( j - 1 ) * lda ] * tscal;
          } else {
            tjjs = tscal;
            if ( tscal == 1. ) goto150 = true;
          }
          if ( ! goto150 ) {
            tjj = Math.abs( tjjs );
            if ( tjj > smlnum ) {
              if ( tjj < 1. ) {
                if ( xj > tjj * bignum ) {
                  rec = 1. / xj;
                  Blas1.dscal( n, rec, x, 1, ioffx );
                  scaleReference.setValue( scaleReference.getValue() * rec );
                  xmax *= rec;
                }
              }
              x[ ioffx + j - 1 ] /= tjjs;
            } else if ( tjj > 0. ) {
              if ( xj > tjj * bignum ) {
                rec = ( tjj * bignum ) / xj;
                Blas1.dscal( n, rec, x, 1, ioffx );
                scaleReference.setValue( scaleReference.getValue() * rec );
                xmax *= rec;
              }
              x[ ioffx + j - 1 ] /= tjjs;
            } else {
              for ( i = 1; i <= n; i ++ ) x[ ioffx + i - 1 ] = 0.;
              x[ ioffx + j - 1 ] = 1.;
              scaleReference.setValue( 0. );
              xmax = 0.;
            }
          } // 150
        } else x[ ioffx + j - 1 ] = x[ ioffx + j - 1 ] / tjjs - sumj;
      }
      xmax = Math.max( xmax, Math.abs( x[ ioffx + j - 1 ] ) );
    }
    scaleReference.setValue( scaleReference.getValue() / tscal );
  }
  if ( tscal != 1. ) Blas1.dscal( n, 1. / tscal, cnorm, 1, ioffcnorm );
}
LaPack1.zlatrs = function( uplo, trans, diag, normin, n, A,
lda, x, scaleReference, cnorm, info ) {
  throw new Error("not programmed: complex array");
}
//************************************************************************
LaPack1.dlauum = function( uplo, n, A, lda, info, ioffa ) {
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -4 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlauum', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  var nb = LaPack0.ilaenv( 1, 'dlauum', uplo, n, -1, -1, -1 );
  if ( nb <= 1 || nb >= n ) {
    LaPack0.dlauu2( uplo, n, A, lda, info, ioffa );
  } else {
    if ( upper ) {
      for ( var i = 1; i <= n; i += nb ) {
        var ib = Math.min( nb, n - i + 1 );
        Blas3.dtrmm( 'Right', 'Upper', 'Transpose', 'Non-unit', i - 1,
          ib, 1., A, lda, A, lda, ioffa + i - 1 + ( i - 1 ) * lda,
          ioffa + ( i - 1 ) * lda );
        LaPack0.dlauu2( 'Upper', ib, A, lda, info,
          ioffa + i - 1 + ( i - 1 ) * lda );
        if ( i + ib <= n ) {
          Blas3.dgemm( 'No transpose', 'Transpose', i - 1, ib,
            n - i - ib + 1, 1., A, lda, A, lda, 1., A, lda,
            ioffa + ( i + ib - 1 ) * lda,
            ioffa + i - 1 + ( i + ib - 1 ) * lda,
            ioffa + ( i - 1 ) * lda );
          Blas3.dsyrk( 'Upper', 'No Transpose', ib, n - i - ib + 1,
            1., A, lda, 1., A, lda,
            ioffa + i - 1 + ( i + ib - 1 ) * lda,
            ioffa + i - 1 + ( i - 1 ) * lda );
        }
      }
    } else {
      for ( i = 1; i <= n; i += nb ) {
        ib = Math.min( nb, n - i + 1 );
        Blas3.dtrmm( 'Left', 'Lower', 'Transpose', 'Non-unit', ib,
          i - 1, 1., A, lda, A, lda, ioffa + i - 1 + ( i - 1 ) * lda,
          ioffa + i - 1 );
        LaPack0.dlauu2( 'Lower', ib, A, lda, info,
          ioffa + i - 1 + ( i - 1 ) * lda );
        if ( i + ib <= n ) {
          Blas3.dgemm( 'Transpose', 'No transpose', ib, i - 1,
            n - i - ib + 1, 1., A, lda, A, lda, 1., A, lda,
            ioffa + i + ib - 1 + ( i - 1 ) * lda,
            ioffa + i + ib - 1, ioffa + i - 1 );
          Blas3.dsyrk( 'Lower', 'Transpose', ib, n - i - ib + 1, 1.,
            A, lda, 1., A, lda, ioffa + i + ib - 1 + ( i - 1 ) * lda,
            ioffa + i - 1 + ( i - 1 ) * lda );
        }
      }
    }
  }
}
LaPack1.zlauum = function( uplo, n, A, lda, info ) {
  throw new Error("not programmed: complex array");
}
//************************************************************************
LaPack1.dormr3 = function( side, trans, m, n, k, l, A, lda,
tau, C, ldc, work, info, ioffa, iofftau, ioffc, ioffwork) {
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
  else if ( l < 0 || ( left && l > m ) || ( ! left && l > n ) ) {
    info.setValue( -6 );
  } else if ( lda < Math.max( 1, k ) ) info.setValue( -8 );
  else if ( ldc < Math.max( 1, m ) ) info.setValue( -11 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dormr3', - info.getValue() );
    return;
  }
  if ( m == 0 || n == 0 || k == 0 ) return;
  if ( left && ! notran || ! left && notran ) {
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
    var ja = m - l + 1;
    var jc = 1;
  } else {
    var mi = m;
    ja = n - l + 1;
    var ic = 1;
  }
  for ( var i = i1;
  ( i3 == 1 && i <= i2 || i3 == -1 && i >= i2 ); i += i3 ) {
    if ( left ) {
      mi = m - i + 1;
      ic = i;
    } else {
      ni = n - i + 1;
      jc = i;
    }
    LaPack0.dlarz( side, mi, ni, l, A, lda, tau[ iofftau + i - 1 ],
      C, ldc, work, ioffa + i - 1 + ( ja - 1 ) * lda,
      ioffc + ic - 1 + ( jc - 1 ) * ldc, ioffwork );
  }
}
//************************************************************************
LaPack1.dpbrfs = function( uplo, n, kd, nrhs, AB, ldab, AFB,
ldafb, B, ldb, X, ldx, ferr, berr, work, iwork, info ) {
  throw new Error("not programmed: banded array");
}
LaPack1.zpbrfs = function( uplo, n, kd, nrhs, AB, ldab, AFB,
ldafb, B, ldb, X, ldx, ferr, berr, work, iwork, info ) {
  throw new Error("not programmed: complex banded array");
}
//************************************************************************
LaPack1.dpbtrf = function( uplo, n, kd, AB, ldab, info ) {
  throw new Error("not programmed: banded array");
}
LaPack1.zpbtrf = function( uplo, n, kd, AB, ldab, info ) {
  throw new Error("not programmed: complex banded array");
}
//************************************************************************
LaPack1.dpftrs = function( transr, uplo, n, nrhs, A, B, ldb,
info, ioffa, ioffb ) {
  throw new Error("not programmed: RFP storage");
}
//************************************************************************
LaPack1.dpoequb = function( n, A, lda, s, scondReference,
amaxReference, info, ioffa, ioffs ) {
  info.setValue( 0 );
  if ( n < 0 ) info.setValue( -1 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -3 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dpoequb', - info.getValue() );
  }
  if ( n == 0 ) {
    scondReference.setValue( 1. );
    amaxReference.setValue( 0. );
    return;
  }
  var base = LaPack0.dlamch( 'B' );
  var tmp = - 0.5 / Math.log( base );
  s[ ioffs ] = A[ ioffa ];
  var smin = s[ ioffs ];
  amaxReference.setValue( s[ ioffs ] );
  for ( var i = 2; i <= n; i ++ ) {
    s[ ioffs + i - 1 ] = A[ ioffa + i - 1 + ( i - 1 ) * lda ];
    smin = Math.min( smin, s[ ioffs + i - 1 ] );
    amaxReference.setValue( Math.max( amaxReference.getValue(),
      s[ ioffs + i - 1 ] ) );
  }
  if ( smin <= 0. ) {
    for ( i = 1; i <= n; i ++ ) {
      if ( s[ ioffs + i - 1 ] <= 0. ) {
        info.setValue( i );
        return;
      }
    }
  } else {
    for ( i = 1; i <= n; i ++ ) {
      s[ ioffs + i - 1 ] =
        Math.pow( base,
        Math.round( tmp * Math.log( s[ ioffs + i - 1 ] ) ) );
    }
    scondReference.setValue( Math.sqrt( smin )
      / Math.sqrt( amaxReference.getValue() ) );
  }
}
//************************************************************************
LaPack1.dporfs = function( uplo, n, nrhs, A, lda, AF, ldaf, B,
ldb, X, ldx, ferr, berr, work, iwork, info, ioffa, ioffaf, ioffb, ioffx,
ioffferr, ioffberr, ioffwork, ioffiwork ) {
  throw new Error("not tested");
  var itmax = 5;
  var isave = new Array( 3 );
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( nrhs < 0 ) info.setValue( -3 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -5 );
  else if ( ldaf < Math.max( 1, n ) ) info.setValue( -7 );
  else if ( ldb < Math.max( 1, n ) ) info.setValue( -9 );
  else if ( ldx < Math.max( 1, n ) ) info.setValue( -11 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dporfs', - info.getValue() );
    return;
  }
  if ( n == 0 || nrhs == 0 ) {
    for ( var j = 1; j <= nrhs; j ++ ) {
      ferr[ ioffferr + j - 1 ] = 0.;
      berr[ ioffberr + j - 1 ] = 0.;
    }
    return;
  }
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
      Blas2.dsymv( uplo, n, -1., A, lda, X, 1, 1., work, 1,
        ioffa, ioffx + ( j - 1 ) * ldx, ioffwork + n );
      for ( var i = 1; i <= n; i ++ ) {
        work[ ioffwork + i - 1 ] =
          Math.abs( B[ ioffb + i - 1 + ( j - 1 ) * ldb ] );
      }
      if ( upper ) {
        for ( var k = 1; k <= n; k ++ ) {
          var s = 0.;
          var xk =
            Math.abs( X[ ioffx + k - 1 + ( j - 1 ) * ldx ] );
          for ( i = 1; i <= k - 1; i ++ ) {
            work[ ioffwork + i - 1 ] +=
              Math.abs( A[ ioffa + i - 1 + ( k - 1 ) * lda ] ) * xk;
            s += Math.abs( A[ ioffa + i - 1 + ( k - 1 ) * lda ] )
              + Math.abs( X[ ioffx + i - 1 + ( j - 1 ) * ldx ] );
          }
          work[ ioffwork + k - 1 ] +=
            Math.abs( A[ ioffa + k - 1 + ( k - 1 ) * lda ] ) * xk + s;
        }
      } else {
        for ( k = 1; k <= n; k ++ ) {
          s = 0.;
          xk = Math.abs( X[ ioffx + k - 1 + ( j - 1 ) * ldx ] );
          work[ ioffwork + k - 1 ] +=
            Math.abs( A[ ioffa + k - 1 + ( k - 1 ) * lda ] ) * xk;
          for ( i = k + 1; i <= n; i ++ ) {
            work[ ioffwork + i - 1 ] +=
              Math.abs( A[ ioffa + i - 1 + ( k - 1 ) * lda ] ) * xk;
            s += Math.abs( A[ ioffa + i - 1 + ( k - 1 ) * lda ] )
              + Math.abs( X[ ioffx + i - 1 + ( j - 1 ) * ldx ] );
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
          s = Math.max( s, ( Math.abs( work[ ioffwork + n + i - 1 ] )
            + safe1 )/ ( work[ ioffwork + i - 1 ] + safe1 ) );
        }
      }
      berr[ ioffberr + j - 1 ] = s;
      if ( berr[ ioffberr + j - 1 ] > eps &&
      2. * berr[ ioffberr + j - 1 ] <= lstres && count <= itmax ) {
        LaPack0.dpotrs( uplo, n, 1, AF, ldaf, work, n, info,
          ioffaf, ioffwork + n );
        Blas1.daxpy( n, 1., work, 1, X, 1, ioffwork + n,
          ioffx + ( j - 1 ) * ldx );
        lstres = berr[ ioffberr + j - 1 ];
        count ++;
        continue;
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
    var estReference = new NumberReference( );
    while ( true ) { // 100
      est.setValue( ferr[ ioffferr + j - 1 ] );
      LaPack0.dlacn2( n, work, work, iwork, est, kase, isave,
        ioffwork + 2 * n, ioffwork + n, ioffiwork, 0 );
      ferr[ ioffferr + j - 1 ] = est.getValue();
      if ( kase.getValue() != 0 ) {
        if ( kase.getValue() == 1 ) {
          LaPack0.dpotrs( uplo, n, 1, AF, ldaf, work, n, info,
            ioffaf, ioffwork + n );
          for ( i = 1; i <= n; i ++ ) {
            work[ ioffwork + n + i - 1 ] *= work[ ioffwork + i - 1 ];
          }
        } else if ( kase.getValue() == 2 ) {
          for ( i = 1; i <= n; i ++ ) {
            work[ ioffwork + n + i - 1 ] *= work[ ioffwork + i - 1 ];
          }
          LaPack0.dpotrs( uplo, n, 1, AF, ldaf, work, n, info,
            ioffaf, ioffwork + n );
        }
        continue;
      } else break;
    }
    lstres = 0.;
    for ( i = 1; i <= n; i ++ ) {
      lstres = Math.max( lstres,
        Math.abs( X[ ioffx + i - 1 + ( j - 1 ) * ldx ] ) );
    }
    if ( lstres != 0. ) ferr[ ioffferr + j - 1 ] /= lstres;
  } // 140
}
LaPack1.zporfs = function( uplo, n, nrhs, A, lda, AF, ldaf, B,
ldb, X, ldx, ferr, berr, work, iwork, info) {
  throw new Error("not programmed: complex array");
}
//************************************************************************
LaPack1.dpotrf = function( uplo, n, A, lda, info, ioffa ) {
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -4 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dpotrf', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  var nb = LaPack0.ilaenv( 1, 'dpotrf', uplo, n, -1, -1, -1 );
  var goto30 = false;
  if ( nb <= 1 || nb >= n ) {
    LaPack0.dpotf2( uplo, n, A, lda, info, ioffa );
  } else {
    if ( upper ) {
      for ( var j = 1; j <= n; j += nb ) {
        var jb = Math.min( nb, n - j + 1 );
        Blas3.dsyrk( 'Upper', 'Transpose', jb, j - 1, -1., A, lda, 1.,
          A, lda, ioffa + ( j - 1 ) * lda,
          ioffa + j - 1 + ( j - 1 ) * lda );
        LaPack0.dpotf2( 'Upper', jb, A, lda, info,
          ioffa + j - 1 + ( j - 1 ) * lda );
        if ( info.getValue() != 0 ) {
          goto30 = true;
          break;
        }
        if ( j + jb <= n ) {
          Blas3.dgemm( 'Transpose', 'No transpose', jb, n - j - jb + 1,
            j - 1, -1., A, lda, A, lda, 1., A, lda,
            ioffa + ( j - 1 ) * lda, ioffa + ( j + jb - 1 ) * lda,
            ioffa + j - 1 + ( j + jb - 1 ) * lda );
          Blas3.dtrsm( 'Left', 'Upper', 'Transpose', 'Non-unit', jb,
            n - j - jb + 1, 1., A, lda, A, lda,
            ioffa + j - 1 + ( j - 1 ) * lda,
            ioffa + j - 1 + ( j + jb - 1 ) * lda );
        }
      }
    } else {
      for ( j = 1; j <= n; j += nb ) {
        jb = Math.min( nb, n - j + 1 );
        Blas3.dsyrk( 'Lower', 'No transpose', jb, j - 1, -1., A, lda,
          1., A, lda, ioffa + j - 1, ioffa + j - 1 + ( j - 1 ) * lda );
        LaPack0.dpotf2( 'Lower', jb, A, lda, info,
          ioffa + j - 1 + ( j - 1 ) * lda );
        if ( info.getValue() != 0 ) {
          goto30 = true;
          break;
        }
        if ( j + jb <= n ) {
          Blas3.dgemm( 'No transpose', 'Transpose', n - j - jb + 1,
            jb, j - 1, -1., A, lda, A, lda, 1., A, lda,
            ioffa + j + jb - 1, ioffa + j - 1,
            ioffa + j + jb - 1 + ( j - 1 ) * lda );
          Blas3.dtrsm( 'Right', 'Lower', 'Transpose', 'Non-unit', 
            n - j - jb + 1, jb, 1., A, lda, A, lda,
            ioffa + j - 1 + ( j - 1 ) * lda,
            ioffa + j + jb - 1 + ( j - 1 ) * lda );
        }
      }
    }
  }
  if ( goto30 ) info.setValue( info.getValue() + j - 1 );
}
LaPack1.zpotrf = function( uplo, n, A, lda, info, ioffa ) {
  throw new Error("not programmed: complex array");
}
//************************************************************************
LaPack1.dpprfs = function( uplo, n, nrhs, AP, AFB, B, ldb, X,
ldx, ferr, berr, Work, iwork, info ) {
  throw new Error("not programmed: packed array");
}
//************************************************************************
LaPack1.dpstf2 = function( uplo, n, A, lda, piv, rank, tol,
work, info, ioffa, ioffpiv, ioffwork ) {
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -4 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dpstf2', -info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  for ( var i = 1; i <= n; i ++ ) piv[ ioffpiv + i - 1 ] = i;
  var pvt = 1;
  var ajj = A[ ioffa + pvt - 1 + ( pvt - 1 ) * lda ];
  for ( i = 2; i <= n; i ++ ) {
    if ( A[ ioffa + i - 1 + ( i - 1 ) * lda ] > ajj ) {
      pvt = i;
      ajj = A[ ioffa + pvt - 1 + ( pvt - 1 ) * lda ];
    }
  }
  if ( ajj == 0. || isNaN( ajj ) ) {
    rank.setValue( 0 );
    info.setValue( 1 );
    return;
  }
  var dstop =
    ( tol < 0. ? n * LaPack0.dlamch( 'Epsilon' ) * ajj : tol );
  for ( i = 1; i <= n; i ++ ) work[ ioffwork + i - 1 ] = 0.;
  var goto160 = false;
  if ( upper ) {
    for ( var j = 1; j <= n; j ++ ) {
      for ( i = j; i <= n; i ++ ) {
        if ( j > 1 ) {
          work[ ioffwork + i - 1 ] +=
            Math.pow( A[ ioffa + j - 2 + ( i - 1 ) * lda ], 2 );
        }
        work[ ioffwork + n + i - 1 ] =
          A[ ioffa + i - 1 + ( i - 1 ) * lda ]
          - work[ ioffwork + i - 1 ];
      } // 120
      if ( j > 1 ) {
        var itemp =
          Blas1.maxloc( n + j, 2 * n, work, 1, ioffwork );
        pvt = itemp + j - 1;
        ajj = work[ ioffwork + n + pvt - 1 ];
        if ( ajj <= dstop || isNaN( ajj ) ) {
          A[ ioffa + j - 1 + ( j - 1 ) * lda ] = ajj;
          goto160 = true;
          break;
        }
      }
      if (j != pvt ) {
        A[ ioffa + pvt - 1 + ( pvt - 1 ) * lda ] = 
          A[ ioffa + j - 1 + ( j - 1 ) * lda ];
        Blas1.dswap( j - 1, A, 1, A, 1, ioffa + ( j - 1 ) * lda,
          ioffa + ( pvt - 1 ) * lda );
        if ( pvt < n ) {
          Blas1.dswap( n - pvt, A, lda, A, lda,
            ioffa + j - 1 + pvt * lda,
            ioffa + pvt - 1 + pvt * lda );
        }
        Blas1.dswap( pvt - j - 1, A, lda, A, 1,
          ioffa + j - 1 + j * lda, ioffa + j + ( pvt - 1 ) * lda );
        var dtemp = work[ ioffwork + j - 1 ];
        work[ ioffwork + j - 1 ] = work[ ioffwork + pvt - 1 ];
        work[ ioffwork + pvt - 1 ] = dtemp;
        itemp = piv[ ioffpiv + pvt - 1 ];
        piv[ ioffpiv + pvt - 1 ] = piv[ ioffpiv + j - 1 ];
        piv[ ioffpiv + j - 1 ] = itemp;
      }
      ajj = Math.sqrt( ajj );
      A[ ioffa + j - 1 + ( j - 1 ) * lda ] = ajj;
      if ( j < n ) {
        Blas2.dgemv( 'Trans', j - 1, n - j, -1., A, lda, A, 1,
          1., A, lda, ioffa + j * lda, ioffa + ( j - 1 ) * lda,
          ioffa + j - 1 + j * lda );
        Blas1.dscal( n - j, 1. / ajj, A, lda,
          ioffa + j - 1 + j * lda );
      }
    } // 130
  } else {
    for ( j = 1; j <= n; j ++ ) {
      for ( i = j; i <= n; i ++ ) {
        if ( j > 1 ) {
          work[ ioffwork + i - 1 ] +=
            Math.pow( A[ ioffa + i - 1 + ( j - 2 ) * lda ], 2 );
        }
        work[ ioffwork + n + i - 1 ] =
          A[ ioffa + i - 1 + ( i - 1 ) * lda ]
          - work[ ioffwork + i - 1 ];
      } // 140
      if ( j > 1 ) {
        itemp = Blas1.maxloc( n + j, 2 * n, work, 1, ioffwork );
        pvt = itemp + j - 1;
        ajj = work[ ioffwork + n + pvt - 1 ];
        if ( ajj <= dstop || isNaN( ajj ) ) {
          A[ ioffa + j - 1 + ( j - 1 ) * lda ] = ajj;
          goto160 = true;
          break;
        }
      }
      if (j != pvt ) {
        A[ ioffa + pvt - 1 + ( pvt - 1 ) * lda ] = 
          A[ ioffa + j - 1 + ( j - 1 ) * lda ];
        Blas1.dswap( j - 1, A, lda, A, lda, ioffa + j - 1,
          ioffa + pvt - 1 );
        if ( pvt < n ) {
          Blas1.dswap( n - pvt, A, 1, A, 1,
            ioffa + pvt + ( j - 1 ) * lda,
            ioffa + pvt + ( pvt - 1 ) * lda );
        }
        Blas1.dswap( pvt - j - 1, A, 1, A, lda,
          ioffa + j + ( j - 1 ) * lda, ioffa + pvt - 1 + j * lda );
        dtemp = work[ ioffwork + j - 1 ];
        work[ ioffwork + j - 1 ] = work[ ioffwork + pvt - 1 ];
        work[ ioffwork + pvt - 1 ] = dtemp;
        itemp = piv[ ioffpiv + pvt - 1 ];
        piv[ ioffpiv + pvt - 1 ] = piv[ ioffpiv + j - 1 ];
        piv[ ioffpiv + j - 1 ] = itemp;
      }
      ajj = Math.sqrt( ajj );
      A[ ioffa + j - 1 + ( j - 1 ) * lda ] = ajj;
      if ( j < n ) {
        Blas2.dgemv( 'No trans', n - j, j - 1, -1., A, lda, A, lda,
          1., A, 1, ioffa + j, ioffa + j - 1,
          ioffa + j + ( j - 1 ) * lda );
        Blas1.dscal( n - j, 1. / ajj, A, 1,
          ioffa + j + ( j - 1 ) * lda );
      }
    } // 150
  }
  if ( !goto160 ) {
    rank.setValue( n );
    return;
  }
  rank.setValue( j - 1 );
  info.setValue( 1 );
}
//************************************************************************
LaPack1.dpttrs = function( n, nrhs, d, e, B, ldb, info, ioffd,
ioffe, ioffb) {
  info.setValue( 0 );
  if ( n < 0 ) info.setValue( -1 );
  else if ( nrhs < 0 ) info.setValue( -2 );
  else if ( ldb < Math.max( 1, n ) ) info.setValue( -6 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dpttrs', -info.getValue() );
    return;
  }
  if ( n == 0 || nrhs == 0 ) return;
  var nb = ( nrhs == 1 ? 1 : Math.max( 1,
    LaPack0.ilaenv( 1, 'dpttrs', ' ', n, nrhs, -1, -1 ) ) );
  if ( nb >= nrhs ) {
    LaPack0.dptts2( n, nrhs, d, e, B, ldb, ioffd, ioffe, ioffb );
  } else {
    for ( var j = 1; j <= nrhs; j += nb ) {
      var jb = Math.min( nrhs - j + 1, nb );
      LaPack0.dptts2( n, jb, d, e, B, ldb,
        ioffd, ioffe, ioffb + ( j - 1 ) * ldb );
    }
  }
}
LaPack1.zpttrs = function( n, nrhs, d, e, B, ldb, info ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack1.drscl = function( n, sa, sx, incx, ioffsx ) {
  if ( n <= 0 ) return;
  var smlnumReference =
    new NumberReference( LaPack0.dlamch( 'S' ) );
  var bignumReference =
    new NumberReference( 1. / smlnumReference.getValue() );
  LaPack0.dlabad( smlnumReference, bignumReference );
  var cden = sa;
  var cnum = 1.;
  var done = false;
  do {
    var cden1 = cden * smlnumReference.getValue();
    var cnum1 = cnum / bignumReference.getValue();
    if ( Math.abs( cden1 ) > Math.abs( cnum ) && cnum != 0. ) {
      var mul = smlnumReference.getValue();
      cden = cden1;
    } else if ( Math.abs( cnum1 ) > Math.abs( cden ) ) {
      mul = bignumReference.getValue();
      cnum = cnum1;
    } else {
      mul = cnum / cden;
      done = true;
    }
    Blas1.dscal( n, mul, sx, incx, ioffsx );
  } while ( ! done );
}
LaPack1.zdrscl = function( n, sa, sx, incx) {
  throw new Error("not programmed: complex array");
}
//************************************************************************
LaPack1.dspcon = function( uplo, n, AP, ipiv, anorm,
rcondReference, work, iwork, info ) {
  throw new Error("not programmed: packed array");
}
LaPack1.zspcon = function( uplo, n, AP, ipiv, anorm,
rcondReference, work, iwork, info ) {
  throw new Error("not programmed: complex packed array");
}
//************************************************************************
LaPack1.dsprfs = function( uplo, n, nrhs, AP, AFP, ipiv, B,
ldb, X, ldx, ferr, berr, work, iwork, info ) {
  throw new Error("not programmed: packed array");
}
LaPack1.zsprfs = function( uplo, n, nrhs, AP, AFP, ipiv, B,
ldb, X, ldx, ferr, berr, work, iwork, info ) {
  throw new Error("not programmed: complex packed array");
}
//************************************************************************
LaPack1.dstebz = function( range, order, n, vl, vu, il, iu,
abstol, d, e, m, nsplit, w, iblock, isplit, work, iwork, info, ioffd,
ioffe, ioffw, ioffiblock, ioffisplit, ioffwork, ioffiwork) {
  var fudge =2.1;
  var relfac = 2.;
  var idumma = new Array( 1 );
  info.setValue( 0 );
  var irange = 0;
  if ( range.charAt(0).toUpperCase() == 'A' ) irange = 1;
  else if ( range.charAt(0).toUpperCase() == 'V' ) irange = 2;
  else if ( range.charAt(0).toUpperCase() == 'I' ) irange = 3;
  var iorder = 0;
  if ( order.charAt(0).toUpperCase() == 'B' ) iorder = 2;
  else if ( order.charAt(0).toUpperCase() == 'E' ) iorder = 1;
  if ( irange <= 0 ) info.setValue( -1 );
  else if ( iorder <= 0 ) info.setValue( -2 );
  else if ( n < 0 ) info.setValue( -3 );
  else if ( irange == 2 ) {
    if ( vl >= vu ) info.setValue( -5 );
  } else if ( irange == 3 && ( il < 1 || il > Math.max( 1, n ) ) ) {
    info.setValue( -6 );
  } else if ( irange == 3 && ( iu < Math.min( n, il ) || iu > n ) ) {
    info.setValue( -7 );
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dstebz', - info.getValue() );
    return;
  }
  info.setValue( 0 );
  var ncnvrg = false;
  var toofew = false;
  m.setValue( 0 );
  if ( n == 0 ) return;
  if ( irange == 3 && il == 1 && iu == n ) irange = 1;
  var safemn = LaPack0.dlamch( 'S' );
  var ulp = LaPack0.dlamch( 'P' );
  var rtoli = ulp * relfac;
  var nb = LaPack0.ilaenv( 1, 'dstebz', ' ', n, -1, -1, -1 );
  if ( nb <= 1 ) nb = 0;
  if ( n == 1 ) {
    nsplit.setValue( 1 );
    isplit[ ioffisplit ] = 1;
    if ( irange == 2 && ( vl >= d[ ioffd ] || vu < d[ ioffd ] ) ) {
      m.setValue( 0 );
    } else {
      w[ ioffw ] = d[ ioffd ];
      iblock[ ioffiblock ] = 1;
      m.setValue( 1 );
    }
    return;
  }
  nsplit.setValue( 1 );
  work[ ioffwork + n - 1 ] = 0.;
  var pivmin = 1.;
  for ( var j = 2; j <= n; j ++ ) {
    var tmp1 = Math.pow( e[ ioffe + j - 2 ], 2 );
    if ( Math.abs( d[ ioffd + j - 1 ] * d[ ioffd + j - 2 ] )
    * Math.pow( ulp, 2 ) + safemn > tmp1 ) {
      isplit[ ioffisplit + nsplit.getValue() - 1 ] = j - 1;
      nsplit.setValue( nsplit.getValue() + 1 );
      work[ ioffwork + j - 2 ] = 0.;
    } else {
      work[ ioffwork + j - 2 ] = tmp1;
      pivmin = Math.max( pivmin, tmp1 );
    }
  }
  isplit[ ioffisplit + nsplit.getValue() - 1 ] = n;
  pivmin *= safemn;
  var iinfo = new IntReference( 0 );
  var iout = new IntReference( 0 );
  if ( irange == 3 ) {
    var gu = d[ ioffd ];
    var gl = d[ ioffd ];
    tmp1 = 0.;
    for ( j = 1; j <= n - 1; j ++ ) {
      var tmp2 = Math.sqrt( work[ ioffwork + j - 1 ] );
      gu = Math.max( gu, d[ ioffd + j - 1 ] + tmp1 + tmp2 );
      gl = Math.min( gl, d[ ioffd + j - 1 ] - tmp1 - tmp2 );
      tmp1 = tmp2;
    }
    gu = Math.max( gu, d[ ioffd + n - 1 ] + tmp1 );
    gl = Math.min( gl, d[ ioffd + n - 1 ] - tmp1 );
    var tnorm = Math.max( Math.abs( gl ), Math.abs( gu ) );
    gl -= fudge * tnorm * ulp * n + fudge * 2. * pivmin;
    gu += fudge * tnorm * ulp * n + fudge * pivmin;
    var itmax = Math.round( ( Math.log( tnorm + pivmin )
      - Math.log( pivmin ) ) / Math.log( 2. ) ) + 2;
    var atoli = ( abstol <= 0. ? ulp * tnorm : abstol );
    work[ ioffwork + n ] = gl;
    work[ ioffwork + n + 1 ] = gl;
    work[ ioffwork + n + 2 ] = gu;
    work[ ioffwork + n + 3 ] = gu;
    work[ ioffwork + n + 4 ] = gl;
    work[ ioffwork + n + 5 ] = gu;
    iwork[ ioffiwork ] = -1;
    iwork[ ioffiwork + 1 ] = -1;
    iwork[ ioffiwork + 2 ] = n + 1;
    iwork[ ioffiwork + 3 ] = n + 1;
    iwork[ ioffiwork + 4 ] = il - 1;
    iwork[ ioffiwork + 5 ] = iu;
    LaPack0.dlaebz( 3, itmax, n, 2, 2, nb, atoli, rtoli, pivmin, d, e,
      work, iwork, work, work, iout, iwork, w, iblock, iinfo , ioffd,
      ioffe, ioffwork, ioffiwork + 4, ioffwork + n, ioffwork + n + 4,
      ioffiwork, ioffw, ioffiblock );
    if ( iwork[ ioffiwork + 5 ] == iu ) {
      var wl = work[ ioffwork + n ];
      var wlu = work[ ioffwork + n + 2 ];
      var nwl = iwork[ ioffiwork ];
      var wu = work[ ioffwork + n + 3 ];
      var wul = work[ ioffwork + n + 1 ];
      var nwu = iwork[ ioffiwork + 3 ];
    } else {
      wl = work[ ioffwork + n  + 1 ];
      wlu = work[ ioffwork + n + 3 ];
      nwl = iwork[ ioffiwork + 1 ];
      wu = work[ ioffwork + n + 2 ];
      wul = work[ ioffwork + n ];
      nwu = iwork[ ioffiwork + 2 ];
    }
    if ( nwl < 0 || nwl >= n || nwu < 1 || nwu > n ) {
      info.setValue( 4 );
      return;
    }
  } else {
    tnorm = Math.max( Math.abs( d[ ioffd ] ) + Math.abs( e[ ioffe ] ),
      Math.abs( d[ ioffd + n - 1 ] )
      + Math.abs( e[ ioffe + n - 2 ] ) );
    for ( j = 2; j <= n - 1; j ++ ) {
      tnorm = Math.max( tnorm, Math.abs( d[ ioffd + j - 1 ] )
        + Math.abs( e[ ioffe + j - 2 ] )
        + Math.abs( e[ ioffe + j - 1 ] ) );
    }
    atoli = ( abstol <= 0. ? ulp * tnorm : abstol );
    if ( irange == 2 ) {
      wl = vl;
      wu = vu;
    } else {
      wl = 0.;
      wu = 0.;
    }
  }
  m.setValue( 0 );
  var iend = 0;
  info.setValue( 0 );
  nwl = 0;
  nwu = 0;
  for ( var jb = 1; jb <= nsplit.getValue(); jb ++ ) {
    var ioff = iend;
    var ibegin = ioff + 1;
    iend = isplit[ ioffisplit + jb - 1 ];
    var in2 = iend - ioff;
    if ( in2 == 1 ) {
      if ( irange == 1 || wl >= d[ ioffd + ibegin - 1 ] - pivmin ) {
        nwl ++;
      }
      if ( irange == 1 || wu >= d[ ioffd + ibegin - 1 ] - pivmin ) {
        nwu ++;
      }
      if ( irange == 1 || ( wl < d[ ioffd + ibegin - 1 ] - pivmin &&
      wu >= d[ ioffd + ibegin - 1 ] - pivmin ) ) {
        m.setValue( m.getValue() + 1 );
        w[ ioffw + m.getValue() - 1 ] = d[ ioffd + ibegin - 1 ];
        iblock[ ioffiblock + m.getValue() - 1 ] = jb;
      }
    } else {
      gu = d[ ioffd + ibegin - 1 ];
      gl = d[ ioffd + ibegin - 1 ];
      tmp1 = 0.;
      for ( j = ibegin; j <= iend - 1; j ++ ) {
        tmp2 = Math.abs( e[ ioffe + j - 1 ] );
        gu = Math.max( gu, d[ ioffd + j - 1 ] + tmp1 + tmp2 );
        gl = Math.min( gl, d[ ioffd + j - 1 ] - tmp1 - tmp2 );
        tmp1 = tmp2;
      }
      gu = Math.max( gu, d[ ioffd + iend - 1 ] + tmp1 );
      gl = Math.min( gl, d[ ioffd + iend - 1 ] - tmp1 );
      var bnorm = Math.max( Math.abs( gl ), Math.abs( gu ) );
      gl -= fudge * bnorm * ulp * in2 + fudge * pivmin;
      gu += fudge * bnorm * ulp * in2 + fudge * pivmin;
      atoli = ( abstol <= 0. ?
        ulp * Math.max( Math.abs( gl ), Math.abs( gu ) ) : abstol );
      if ( irange > 1 ) {
        if ( gu < wl ) {
          nwl += in2;
          nwu += in2;
          continue;
        }
        gl = Math.max( gl, wl );
        gu = Math.min( gu, wu );
        if ( gl >= gu ) continue;
      }
      work[ ioffwork + n ] = gl;
      work[ ioffwork + n + in2 ] = gu;
      var imref = new IntReference();
      LaPack0.dlaebz( 1, 0, in2, in2, 1, nb, atoli, rtoli, pivmin,
        d, e, work, idumma, work, work, imref, iwork, w, iblock, iinfo,
        ioffd + ibegin - 1, ioffe + ibegin - 1, ioffwork + ibegin - 1,
        0, ioffwork + n, ioffwork + n + 2 * in2, ioffiwork,
        ioffw + m.getValue(), ioffiblock + m.getValue() );
      nwl += iwork[ ioffiwork ];
      nwu += iwork[ ioffiwork + in2 ];
      var iwoff = m.getValue() - iwork[ ioffiwork ];
      itmax = Math.round( ( Math.log( gu - gl + pivmin )
        - Math.log( pivmin ) ) / Math.log( 2. ) ) + 2;
      LaPack0.dlaebz( 2, itmax, in2, in2, 1, nb, atoli, rtoli, pivmin,
        d, e, work, idumma, work, work, iout, iwork, w, iblock, iinfo,
        ioffd + ibegin - 1, ioffe + ibegin - 1, ioffwork + ibegin - 1,
        0, ioffwork + n, ioffwork + n + 2 * in2, ioffiwork,
        ioffw + m.getValue(), ioffiblock + m.getValue() );
      for ( j = 1; j <= iout.getValue(); j ++ ) {
        tmp1 = 0.5 * ( work[ ioffwork + j + n - 1 ]
          + work[ ioffwork + j + in2 + n -1 ] );
        if ( j > iout.getValue() - iinfo.getValue() ) {
          ncnvrg = true;
          var ib = - jb;
        } else ib = jb;
        for ( var je = iwork[ ioffiwork + j - 1 ] + 1 + iwoff;
        je <= iwork[ ioffiwork + j + in2 - 1 ] + iwoff; je ++ ) {
          w[ ioffw + je - 1 ] = tmp1;
          iblock[ ioffiblock + je - 1 ] = ib;
        }
      } // 60
      m.setValue( m.getValue() + imref.getValue() );
    }
  } // 70
  if ( irange == 3 ) {
    var im = 0;
    var idiscl = il - 1 - nwl;
    var idiscu = nwu - iu;
    if ( idiscl > 0 || idiscu > 0 ) {
      for ( je = 1; je <= m.getValue(); je ++ ) {
        if ( w[ ioffw + je - 1 ] <= wlu && idiscl > 0 ) idiscl --;
        else if ( w[ ioffw + je - 1 ] >= wul && idiscu > 0 ) {
          idiscu --;
        } else {
          im ++;
          w[ ioffw + im - 1 ] = w[ ioffw + je - 1 ];
          iblock[ ioffiblock + im - 1 ] =
            iblock[ ioffiblock + je - 1 ];
        }
      }
      m.setValue( im );
    }
    if ( idiscl > 0 || idiscu > 0 ) {
      if ( idiscl > 0 ) {
        var wkill = wu;
        for ( var jdisc = 1; jdisc <= idiscl; jdisc ++ ) {
          var iw = 0;
          for ( je = 1; je <= m.getValue(); je ++ ) {
            if ( iblock[ ioffiblock + je - 1 ] != 0 &&
            ( w[ ioffw + je - 1 ] < wkill || iw == 0 ) ) {
              iw = je;
              wkill = w[ ioffw + je - 1 ];
            }
          }
          iblock[ ioffiblock + iw ] = 0;
        }
      }
      if ( idiscu > 0 ) {
        wkill = wl;
        for ( jdisc = 1; jdisc <= idiscu; jdisc ++ ) {
          iw = 0;
          for ( je = 1; je <= m.getValue(); je ++ ) {
            if ( iblock[ ioffiblock + je - 1 ] != 0 &&
            ( w[ ioffw + je - 1 ] > wkill || iw == 0 ) ) {
              iw = je;
              wkill = w[ ioffw + je - 1 ];
            }
          }
          iblock[ ioffiblock + iw ] = 0;
        }
      }
      im = 0;
      for ( je = 1; je <= m.getValue(); je ++ ) {
        if ( iblock[ ioffiblock + je - 1 ] != 0 ) {
          im ++;
          w[ ioffw + im - 1 ] = w[ ioffw + je - 1 ];
          iblock[ ioffiblock + im - 1 ] =
            iblock[ ioffiblock + je - 1 ];
        }
      }
      m.setValue( im );
    }
    if ( idiscl < 0 || idiscu < 0 ) toofew = true;
  }
  if ( iorder == 1 && nsplit.getValue() > 1 ) {
    for ( je = 1; je <= m.getValue() - 1; je ++ ) {
      var ie = 0;
      tmp1 = w[ ioffw + je - 1 ];
      for ( j = je + 1; j <= m.getValue(); j ++ ) {
        if ( w[ ioffw + j - 1 ] < tmp1 ) {
          ie = j;
          tmp1 = w[ ioffw + j - 1 ];
        }
      }
      if ( ie != 0 ) {
        var itmp1 = iblock[ ioffiblock + ie - 1 ];
        w[ ioffw + ie - 1 ] = w[ ioffw + je - 1 ];
        iblock[ ioffiblock + ie - 1 ] = iblock[ ioffiblock + je - 1 ];
        w[ ioffw + je - 1 ] = tmp1;
        iblock[ ioffiblock + je - 1 ] = itmp1;
      }
    }
  }
  info.setValue( 0 );
  if ( ncnvrg) info.setValue( info.getValue() + 1 );
  if ( toofew ) info.setValue( info.getValue() + 2 );
}
//************************************************************************
LaPack1.dsycon = function( uplo, n, A, lda, ipiv, anorm,
rcondReference, work, iwork, info, ioffa, ioffipiv, ioffwork,
ioffiwork ) {
  var isave = new Array( 3 );
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -4 );
  else if ( anorm < 0. ) info.setValue( -6 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dsycon', - info.getValue() );
    return;
  }
  rcondReference.setValue( 0. );
  if ( n == 0 ) {
    rcondReference.setValue( 1. );
    return;
  } else if ( anorm <= 0. ) return;
  if ( upper ) {
    for ( var i = n; i >= 1; i -- ) {
      if ( ipiv[ ioffipiv + i - 1 ] > 0 &&
      A[ ioffa + i - 1 + ( i - 1 ) * lda ] == 0. ) {
        return;
      }
    }
  } else {
    for ( i = 1; i <= n; i ++ ) {
      if ( ipiv[ ioffipiv + i - 1 ] > 0 &&
      A[ ioffa + i - 1 + ( i - 1 ) * lda ] == 0. ) {
        return;
      }
    }
  }
  var kaseReference = new IntReference(0);
  var ainvnmReference = new NumberReference(); 
  while ( true ) { // 30
    LaPack0.dlacn2( n, work, work, iwork, ainvnm, kaseReference, isave,
      ioffwork + n, ioffwork, ioffiwork, 0 );
    if ( kaseReference.getValue() != 0 ) {
      LaPack0.dsytrs( uplo, n, 1, A, lda, ipiv, work, n, info,
        ioffa, ioffipiv, ioffwork );
    } else break;
  }
  if ( ainvnm.getValue() != 0. ) {
    rcondReference.setValue( ( 1. / ainvnm.getValue() ) / anorm );
  }
}
LaPack1.zsycon = function( uplo, n, A, lda, ipiv, anorm,
rcondReference, work, iwork, info ) {
  throw new Error("not programmed: complex array");
}
//************************************************************************
LaPack1.dsyequb = function( uplo, n, A, lda, s, scondReference,
amaxReference, work, info, ioffa, ioffs, ioffwork ) {
  var max_iter = 100;
  info.setValue( 0 );
  if ( uplo.charAt(0).toUpperCase() != 'U' &&
  uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -4 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dsyequb', - info.getValue() );
    return;
  }
  var up = ( uplo.charAt(0).toUpperCase() == 'U' );
  amaxReference.setValue( 0. );
  if ( n == 0 ) {
    scondReference.setValue( 1. );
    return;
  }
  for ( var i = 1; i <= n; i ++ ) s[ ioffs + i - 1 ] = 0.;
//    amaxReference.setValue( 0. );
  if ( up ) {
    for ( var j = 1; j <= n; j ++ ) {
      for ( i = 1; i <= j - 1; i ++ ) {
        var t =
          Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
        s[ ioffs + i - 1 ] = Math.max( s[ ioffs + i - 1 ], t );
        s[ ioffs + j - 1 ] = Math.max( s[ ioffs + j - 1 ], t );
        amaxReference.setValue( Math.max( amaxReference.getValue(), t ) );
      }
      t = Math.abs( A[ ioffa + j - 1 + ( j - 1 ) * lda ] );
      s[ ioffs + j - 1 ] = Math.max( s[ ioffs + j - 1 ], t );
      amaxReference.setValue( Math.max( amaxReference.getValue(), t ) );
    }
  } else {
    for ( j = 1; j <= n; j ++ ) {
      t = Math.abs( A[ ioffa + j - 1 + ( j - 1 ) * lda ] );
      s[ ioffs + j - 1 ] = Math.max( s[ ioffs + j - 1 ], t );
      amaxReference.setValue( Math.max( amaxReference.getValue(), t ) );
      for ( i = j + 1; i <= n; i ++ ) {
        t = Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
        s[ ioffs + i - 1 ] = Math.max( s[ ioffs + i - 1 ], t );
        s[ ioffs + j - 1 ] = Math.max( s[ ioffs + j - 1 ], t );
        amaxReference.setValue( Math.max( amaxReference.getValue(), t ) );
      }
    }
  }
  for ( j = 1; j <= n; j ++ ) {
    s[ ioffs + j - 1 ] = 1. / s[ ioffs + j - 1 ];
  }
  var tol = 1. / Math.sqrt( 2. * n );
  for ( var iter = 1; iter <= max_iter; iter ++ ) {
    var scaleReference = new NumberReference( 0. );
    var sumsqReference = new NumberReference( 0. );
    for ( i = 1; i <= n; i ++ ) work[ ioffwork + i - 1 ] = 0.;
    if ( up ) {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = 1; i <= j - 1; i ++ ) {
          t = Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
          work[ ioffwork + i - 1 ] += t * s[ ioffs + j - 1 ];
          work[ ioffwork + j - 1 ] += t * s[ ioffs + i - 1 ];
        }
        work[ ioffwork + j - 1 ] +=
          Math.abs( A[ ioffa + j - 1 + ( j - 1 ) * lda ] )
          * s[ ioffs + j - 1 ];
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        work[ ioffwork + j - 1 ] +=
          Math.abs( A[ ioffa + j - 1 + ( j - 1 ) * lda ] )
          * s[ ioffs + j - 1 ];
        for ( i = j + 1; i <= n; i ++ ) {
          t = Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
          work[ ioffwork + i - 1 ] += t * s[ ioffs + j - 1 ];
          work[ ioffwork + j - 1 ] += t * s[ ioffs + i - 1 ];
        }
      }
    }
    var avg = 0.;
    for ( i = 1; i <= n; i ++ ) {
      avg += s[ ioffs + i - 1 ] * work[ ioffwork + i - 1 ];
    }
    avg /= n;
    var std = 0.;
    for ( i = 2 * n + 1; i <= 3 * n; i ++ ) {
      work[ ioffwork + i - 1 ] = s[ ioffs + i - 2 * n - 1 ]
        * work[ ioffwork + i - 2 * n - 1 ] - avg;
    }
    LaPack0.dlassq( n, work, 1, scaleReference, sumsqReference,
      ioffwork + 2 * n );
    std = scaleReference.getValue()
      * Math.sqrt( sumsqReference.getValue() / n );
    if ( std < tol * avg ) break;
    for ( i = 1; i <= n; i ++ ) {
      t = Math.abs( A[ ioffa + i - 1 + ( i - 1 ) * lda ] );
      var si = s[ ioffs + i - 1 ];
      var c2 = ( n - 1 ) * t;
      var c1 =
        ( n - 2 ) * ( work[ ioffwork + i - 1 ] - t * si );
      var c0 = - ( t * si ) * si
        + 2. * work[ ioffwork + i - 1 ] * si - n * avg;
      var d = c1 * c1 - 4. * c0 * c2;
      if ( d <= 0. ) {
        info.setValue( -1 );
        return;
      }
      si = - 2. * c0 / ( c1 + Math.sqrt( d ) );
      d = si - s[ ioffs + i - 1 ];
      var u = 0.;
      if ( up ) {
        for ( j = 1; i <= i; j ++ ) {
          t = Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ] );
          u += s[ ioffs + j - 1 ] * t;
          work[ ioffwork + j - 1 ] += d * t;
        }
        for ( j = i + 1; i <= n; j ++ ) {
          t = Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
          u += s[ ioffs + j - 1 ] * t;
          work[ ioffwork + j - 1 ] += d * t;
        }
      } else {
        for ( j = 1; i <= i; j ++ ) {
          t = Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
          u += s[ ioffs + j - 1 ] * t;
          work[ ioffwork + j - 1 ] += d * t;
        }
        for ( j = i + 1; i <= n; j ++ ) {
          t = Math.abs( A[ ioffa + j - 1 + ( i - 1 ) * lda ] );
          u += s[ ioffs + j - 1 ] * t;
          work[ ioffwork + j - 1 ] += d * t;
        }
      }
      avg += ( u + work[ ioffwork + i - 1 ] ) * d / n;
      s[ ioffs + i - 1 ] = si;
    }
  }
  var smlnum = LaPack0.dlamch( 'Safmin' );
  var bignum = 1. / smlnum;
  var smin = bignum;
  var smax = 0.;
  t = 1. / Math.sqrt( avg );
  var base = LaPack0.dlamch( 'B' );
  u = 1. / Math.log( base );
  for ( i = 1; i <= n; i ++ ) {
    s[ ioffs + i - 1 ] = Math.pow( base,
      Math.round( u * Math.log( s[ ioffs + i - 1 ] * t ) ) );
    smin = Math.min( smin, s[ ioffs + i - 1 ] );
    smax = Math.max( smax, s[ ioffs + i - 1 ] );
  }
  scondReference.setValue( Math.max( smin, smlnum )
    / Math.min( smax, bignum ) );
}
//************************************************************************
LaPack1.dsygst = function( itype, uplo, n, A, lda, B, ldb,
info, ioffa, ioffb) {
  throw new Error("not programmed: generalized eigenvalue problem");
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( itype < 1 || itype > 3 ) info.setValue( -1 );
  else if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -2 );
  } else if ( n < 0 ) info.setValue( -3 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -5 );
  else if ( ldb < Math.max( 1, n ) ) info.setValue( -6 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dsygst', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  var nb = LaPack0.ilaenv( 1, 'dsygst', uplo, n, -1, -1, -1 );
  if ( nb <= 1 || nb >= n ) {
    LaPack0.dsygs2( itype, uplo, n, A, lda, B, ldb, info, ioffa,
      ioffb );
  } else {
    if ( itype == 1 ) {
      if ( upper ) {
        for ( var k = 1; k <= n; k += nb ) {
          var kb = Math.min( n - k + 1, nb );
          LaPack0.dsygs2( itype, uplo, kb, A, lda, B, ldb, info, 
            ioffa + k - 1 + ( k - 1 ) * lda,
            ioffb + k - 1 + ( k - 1 ) * ldb );
          if ( k + kb <= n ) {
            Blas3.dtrsm( 'Left', uplo, 'Transpose', 'Non-unit', kb,
              n - k - kb + 1, 1., B, ldb, A, lda,
              ioffb + k - 1 + ( k - 1 ) * ldb,
              ioffa + k - 1 + ( k + kb - 1 ) * lda );
            Blas3.dsymm( 'Left', uplo, kb, n - k - kb + 1, -0.5, A,
              lda, B, ldb, 1., A, lda, ioffa + k - 1 + ( k - 1 ) * lda,
              ioffb + k - 1 + ( k + kb - 1 ) * ldb,
              ioffa + k - 1 + ( k + kb - 1 ) * lda );
            Blas3.dsyr2k( uplo, 'Transpose', n - k - kb + 1, kb, -1.,
              A, lda, B, ldb, 1., A, lda,
              ioffa + k - 1 + ( k + kb - 1 ) * lda,
              ioffb + k - 1 + ( k + kb - 1 ) * ldb,
              ioffa + k + kb - 1 + ( k + kb - 1 ) * lda );
            Blas3.dsymm( 'Left', uplo, kb, n - k - kb + 1, -0.5, A,
              lda, B, ldb, 1., A, lda, ioffa + k - 1 + ( k - 1 ) * lda,
              ioffb + k - 1 + ( k + kb - 1 ) * ldb,
              ioffa + k - 1 + ( k + kb - 1 ) * lda );
            Blas3.dtrsm( 'Right', uplo, 'No transpose', 'Non-unit', kb,
              n - k - kb + 1, 1., B, ldb, A, lda,
              ioffb + k + kb - 1 + ( k + kb - 1 ) * ldb,
              ioffa + k - 1 + ( k + kb - 1 ) * lda );
          }
        }
      } else {
        for ( k = 1; k <= n; k += nb ) {
          kb = Math.min( n - k + 1, nb );
          LaPack0.dsygs2( itype, uplo, kb, A, lda, B, ldb, info, 
            ioffa + k - 1 + ( k - 1 ) * lda,
            ioffb + k - 1 + ( k - 1 ) * ldb );
          if ( k + kb <= n ) {
            Blas3.dtrsm( 'Right', uplo, 'Transpose', 'Non-unit',
              n - k - kb + 1, kb, 1., B, ldb, A, lda,
              ioffb + k - 1 + ( k - 1 ) * ldb,
              ioffa + k + kb - 1 + ( k - 1 ) * lda );
            Blas3.dsymm( 'Right', uplo, n - k - kb + 1, kb, -0.5, A,
              lda, B, ldb, 1., A, lda, ioffa + k - 1 + ( k - 1 ) * lda,
              ioffb + k + kb - 1 + ( k - 1 ) * ldb,
              ioffa + k + kb - 1 + ( k - 1 ) * lda );
            Blas3.dsyr2k( uplo, 'No transpose', n - k - kb + 1, kb,
              -1., A, lda, B, ldb, 1., A, lda,
              ioffa + k + kb - 1 + ( k - 1 ) * lda,
              ioffb + k + kb - 1 + ( k - 1 ) * ldb,
              ioffa + k + kb - 1 + ( k + kb - 1 ) * lda );
            Blas3.dsymm( 'Right', uplo, n - k - kb + 1, kb, -0.5, A,
              lda, B, ldb, 1., A, lda, ioffa + k - 1 + ( k - 1 ) * lda,
              ioffb + k + kb - 1 + ( k - 1 ) * ldb,
              ioffa + k + kb - 1 + ( k - 1 ) * lda );
            Blas3.dtrsm( 'Left', uplo, 'No transpose', 'Non-unit',
              n - k - kb + 1, kb, 1., B, ldb, A, lda,
              ioffb + k + kb - 1 + ( k + kb - 1 ) * ldb,
              ioffa + k + kb - 1 + ( k - 1 ) * lda );
          }
        }
      }
    } else {
      if ( upper ) {
        for ( k = 1; k <= n; k += nb ) {
          kb = Math.min( n - k + 1, nb );
          Blas3.dtrmm( 'Left', uplo, 'No transpose', 'Non-unit', k - 1,
            kb, 1., B, ldb, A, lda, ioffb, ioffa + ( k - 1 ) * lda );
          Blas3.dsymm( 'Right', uplo, k - 1, kb, 0.5, A,
            lda, B, ldb, 1., A, lda, ioffa + k - 1 + ( k - 1 ) * lda,
            ioffb + ( k - 1 ) * ldb, ioffa + ( k - 1 ) * lda );
          Blas3.dsyr2k( uplo, 'No transpose', k - 1, kb, 1., A, lda,
            B, ldb, 1., A, lda, ioffa + ( k - 1 ) * lda,
            ioffb + ( k - 1 ) * ldb, ioffa );
          Blas3.dsymm( 'Right', uplo, k - 1, kb, 0.5, A, lda, B, ldb,
            1., A, lda, ioffa + k - 1 + ( k - 1 ) * lda,
            ioffb + ( k - 1 ) * ldb, ioffa + ( k - 1 ) * lda );
          Blas3.dtrmm( 'Right', uplo, 'Transpose', 'Non-unit', k - 1,
            kb, 1., B, ldb, A, lda, ioffb + k - 1 + ( k - 1 ) * ldb,
            ioffa + ( k - 1 ) * lda );
          LaPack0.dsygs2( itype, uplo, kb, A, lda, B, ldb, info, 
            ioffa + k - 1 + ( k - 1 ) * lda,
            ioffb + k - 1 + ( k - 1 ) * ldb );
        }
      } else {
        for ( k = 1; k <= n; k += nb ) {
          kb = Math.min( n - k + 1, nb );
          Blas3.dtrmm( 'Right', uplo, 'No transpose', 'Non-unit', kb,
            k - 1, 1., B, ldb, A, lda, ioffb, ioffa + k - 1 );
          Blas3.dsymm( 'Left', uplo, kb, k - 1, 0.5, A,
            lda, B, ldb, 1., A, lda, ioffa + k - 1 + ( k - 1 ) * lda,
            ioffb + k - 1, ioffa + k - 1 );
          Blas3.dsyr2k( uplo, 'Transpose', k - 1, kb, 1., A, lda,
            B, ldb, 1., A, lda, ioffa + k - 1, ioffb + k - 1, ioffa );
          Blas3.dsymm( 'Left', uplo, kb, k - 1, 0.5, A, lda, B, ldb,
            1., A, lda, ioffa + k - 1 + ( k - 1 ) * lda,
            ioffb + k - 1, ioffa + k - 1 );
          Blas3.dtrmm( 'Left', uplo, 'Transpose', 'Non-unit', kb,
            k - 1, 1., B, ldb, A, lda, ioffb + k - 1 + ( k - 1 ) * ldb,
            ioffa + k - 1 );
          LaPack0.dsygs2( itype, uplo, kb, A, lda, B, ldb, info, 
            ioffa + k - 1 + ( k - 1 ) * lda,
            ioffb + k - 1 + ( k - 1 ) * ldb );
        }
      }
    }
  }
}
//************************************************************************
LaPack1.dsyrfs = function( uplo, n, nrhs, A, lda, AF, ldaf,
ipiv, B, ldb, X, ldx, ferr, berr, work, iwork, info, ioffa, ioffaf,
ioffipiv, ioffb, ioffx, ioffferr, ioffberr, ioffwork, ioffiwork) {
  var itmax = 5;
  var isave = new Array( 3 );
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( nrhs < 0 ) info.setValue( -3 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -5 );
  else if ( ldaf < Math.max( 1, n ) ) info.setValue( -7 );
  else if ( ldb < Math.max( 1, n ) ) info.setValue( -10 );
  else if ( ldx < Math.max( 1, n ) ) info.setValue( -12 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dsyrfs', - info.getValue() );
    return;
  }
  if ( n == 0 || nrhs == 0 ) {
    for ( var j = 1; j <= nrhs; j ++ ) {
      ferr[ ioffferr + j - 1 ] = 0.;
      berr[ ioffberr + j - 1 ] = 0.;
    }
    return;
  }
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
      Blas2.dsymv( uplo, n, -1., A, lda, X, 1, 1., work, 1, ioffa,
        ioffx + ( j - 1 ) * ldx, ioffwork + n );
      for ( var i = 1; i <= n; i ++ ) {
        work[ ioffwork + i - 1 ] =
          Math.abs( B[ ioffb + i - 1 + ( j - 1 ) * ldb ] ); 
      }
      if ( upper ) {
        for ( var k = 1; k <= n; k ++ ) {
          var s = 0.;
          var xk =
            Math.abs( X[ ioffx + k - 1 + ( j - 1 ) * ldx ] );
          for ( i = 1; i <= k - 1; i ++ ) {
            work[ ioffwork + i - 1 ] +=
              Math.abs( A[ ioffa + i - 1 + ( k - 1 ) * lda ] ) * xk;
            s += Math.abs( A[ ioffa + i - 1 + ( k - 1 ) * lda ] )
              * Math.abs( X[ ioffx + i - 1 + ( j - 1 ) * ldx ] );
          }
          work[ ioffwork + k - 1 ] +=
            Math.abs( A[ ioffa + k - 1 + ( k - 1 ) * lda ] ) * xk + s;
        }
      } else {
        for ( k = 1; k <= n; k ++ ) {
          s = 0.;
          xk = Math.abs( X[ ioffx + k - 1 + ( j - 1 ) * ldx ] );
          for ( i = k + 1; i <= n; i ++ ) {
            work[ ioffwork + i - 1 ] +=
              Math.abs( A[ ioffa + i - 1 + ( k - 1 ) * lda ] ) * xk;
            s += Math.abs( A[ ioffa + i - 1 + ( k - 1 ) * lda ] )
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
          s = Math.max( s, ( Math.abs( work[ ioffwork + n + i - 1 ] )
            + safe1 ) / ( work[ ioffwork + i - 1 ] + safe1 ) );
        }
      }
      berr[ ioffberr + j - 1 ] = s;
      if ( berr[ ioffberr + j - 1 ] > eps &&
      2. * berr[ ioffberr + j - 1 ] <= lstres && count <= itmax ) {
        LaPack0.dsytrs( uplo, n, 1, AF, ldaf, ipiv, work, n, info,
          ioffaf, ioffipiv, ioffwork + n );
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
          + nz * eps * work[ ioffwork + i - 1 ];
      } else {
        work[ ioffwork + i - 1 ] =
          Math.abs( work[ ioffwork + n + i - 1 ] )
          + nz * eps * work[ ioffwork + i - 1 ] + safe1;
      }
    }
    var kase = new IntReference( 0 );
    while ( true ) { // 100
      var estReference =
        new NumberReference( ferr[ ioffferr + j - 1 ] );
      LaPack0.dlacn2( n, work, work, iwork, estReference, kase, isave,
        ioffwork + 2 * n, ioffwork + n, ioffiwork, 0 );
      ferr[ ioffferr + j - 1 ]= estReference.getValue();
      if ( kase.getValue() != 0 ) {
        if ( kase.getValue() == 1 ) {
          LaPack0.dsytrs( uplo, n, 1, AF, ldaf, ipiv, work, n, info,
            ioffaf, ioffipiv, ioffwork + n );
          for ( i = 1; i <= n; i ++ ) {
            work[ ioffwork + n + i - 1 ] *= work[ ioffwork + i - 1 ];
          }
        } else if ( kase.getValue() == 2 ) {
          for ( i = 1; i <= n; i ++ ) {
            work[ ioffwork + n + i - 1 ] *= work[ ioffwork + i - 1 ];
          }
          LaPack0.dsytrs( uplo, n, 1, AF, ldaf, ipiv, work, n, info,
            ioffaf, ioffipiv, ioffwork + n );
        }
      } else break;
    }
    lstres = 0.;
    for ( i = 1; i <= n; i ++ ) {
      lstres = Math.max( lstres,
        Math.abs( X[ ioffx + i - 1 + ( j - 1 ) * ldx ] ) );
    }
    if ( lstres != 0.) ferr[ ioffferr + j - 1 ] /= lstres;
  } // 140
}
LaPack1.zsyrfs = function( uplo, n, nrhs, A, lda, AF, ldaf,
ipiv, B, ldb, X, ldx, ferr, berr, work, iwork, info ) {
  throw new Error("not programmed: complex array");
}
//************************************************************************
LaPack1.dsytrf = function( uplo, n, A, lda, ipiv, work, lwork,
info, ioffa, ioffipiv, ioffwork ) {
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var lquery = ( lwork == -1 );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -4 );
  else if ( lwork < 1 && ! lquery ) info.setValue( -7 );
  if ( info.getValue() == 0 ) {
    var nb = LaPack0.ilaenv( 1, 'dsytrf', uplo, n, -1, -1, -1 );
    var lwkopt = n * nb;
    work[ ioffwork ] = lwkopt;
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dsytrf', - info.getValue() );
    return;
  } else if ( lquery ) return;
  var nbmin = 2;
  var ldwork = n;
  if ( nb > 1 && nb < n ) {
    var iws = ldwork * nb;
    if ( lwork < iws ) {
      nb = Math.max( lwork / ldwork, 1 );
      nbmin = Math.max( 2,
        LaPack0.ilaenv( 2, 'dsytrf', uplo, n, -1, -1, -1 ) );
    }
  } else iws = 1;
  if ( nb < nbmin ) nb = n;
  var kb = new IntReference();
  var iinfo = new IntReference();
  if ( upper ) {
    var k = n;
    while ( true ) { // 10
      if ( k < 1 ) break;
      if ( k > nb ) {
        LaPack0.dlasyf( uplo, k, nb, kb, A, lda, ipiv, work, ldwork,
          iinfo, ioffa, ioffipiv, ioffwork );
      } else {
        LaPack0.dsytf2( uplo, k, A, lda, ipiv, iinfo, ioffa, ioffipiv );
        kb.setValue( k );
      }
      if ( info.getValue() == 0 && iinfo.getValue() > 0 ) {
        info.setValue( iinfo.getValue() );
      }
      k -= kb.getValue();
    }
  } else {
    k = 1;
    while ( true ) { // 20
      if ( k > n ) break;
      if ( k <= n - nb ) {
        LaPack0.dlasyf( uplo, n - k + 1, nb, kb, A, lda, ipiv, work,
          ldwork, iinfo, ioffa + k - 1 + ( k - 1 ) * lda,
          ioffipiv + k - 1, ioffwork );
      } else {
        LaPack0.dsytf2( uplo, n - k + 1, A, lda, ipiv, iinfo, 
          ioffa + k - 1 + ( k - 1 ) * lda, ioffipiv + k - 1 );
        kb.setValue( n - k + 1 );
      }
      if ( info.getValue() == 0 && iinfo.getValue() > 0 ) {
        info.setValue( iinfo.getValue() + k - 1 );
      }
      for ( var j = k; j <= k + kb.getValue() - 1; j ++ ) {
        if ( ipiv[ ioffipiv + j - 1 ] > 0 ) {
          ipiv[ ioffipiv + j - 1 ] += k - 1;
        } else ipiv[ ioffipiv + j - 1 ] -= k - 1;
      }
      k += kb.getValue();
    }
  }
  work[ ioffwork ] = lwkopt;
}
LaPack1.zsytrf = function( uplo, n, A, lda, ipiv, work, lwork,
info ) {
}
//************************************************************************
LaPack1.dsytrs2 = function( uplo, n, nrhs, A, lda, ipiv, B, ldb,
work, info, ioffa, ioffipiv, ioffb, ioffwork) {
  throw new Error("not tested: TRF format");
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( nrhs < 0 ) info.setValue( -3 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -5 );
  else if ( ldb < Math.max( 1, n ) ) info.setValue( -8 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dsytrs2', - info.getValue() );
    return;
  }
  if ( n == 0 || nrhs == 0 ) return;
  var iinfo = new IntReference();
  LaPack0.dsyconv( uplo, 'C', n, A, lda, ipiv, work, iinfo,
    ioffa, ioffipiv, ioffwork );
  if ( upper ) {
    var k = n;
    while ( k >= 1 ) {
      if ( ipiv[ ioffipiv + k - 1 ] > 0 ) {
        var kp = ipiv[ ioffipiv + k - 1 ];
        if ( kp != k ) {
          Blas1.dswap( nrhs, B, ldb, B, ldb,
            ioffb + k - 1, ioffb + kp - 1 );
          k --;
        }
      } else {
        kp = - ipiv[ ioffipiv + k - 1 ];
        if ( kp == - ipiv[ ioffipiv + k - 2 ] ) {
          Blas1.dswap( nrhs, B, ldb, B, ldb,
            ioffb + k - 2, ioffb + kp - 1 );
        }
        k -= 2;
      }
    }
    Blas3.dtrsm( 'L', 'U', 'N', 'U', n, nrhs, 1., A, lda, B, ldb,
      ioffa, ioffb );
    var i = n;
    while ( i >= 1 ) {
      if ( ipiv[ ioffipiv + i - 1 ] > 0 ) {
        Blas1.dscal( nrhs, 1. / A[ ioffa + i - 1 + ( i - 1 ) * lda ],
          B, ldb, ioffb + i - 1 );
      } else if ( i > 1 ) {
        if ( ipiv[ ioffipiv + i - 2 ] == ipiv[ ioffipiv + i - 1 ] ) {
          var akm1k = work[ ioffwork + i - 1 ];
          var akm1 =
            A[ ioffa + i - 2 + ( i - 2 ) * lda ] / akm1k;
          var ak = A[ ioffa + i - 1 + ( i - 1 ) * lda ] / akm1k;
          var denom = akm1 * ak - 1.;
          for ( var j = 1; j <= nrhs; j ++ ) {
            var bkm1 =
              B[ ioffb + i - 2 + ( j - 1 ) * ldb ] / akm1k;
            var bk =
              B[ ioffb + i - 1 + ( j - 1 ) * ldb ] / akm1k;
            B[ ioffb + i - 2 + ( j - 1 ) * ldb ] =
              ( ak * bkm1 - bk ) / denom;
            B[ ioffb + i - 1 + ( j - 1 ) * ldb ] =
              ( akm1 * bk - bkm1 ) / denom;
          }
          i --;
        }
      }
      i --;
    }
    Blas3.dtrsm( 'L', 'U', 'T', 'U', n, nrhs, 1., A, lda, B, ldb,
      ioffa, ioffb );
    k = 1;
    while ( k <= n ) {
      if ( ipiv[ ioffipiv + k - 1 ] > 0 ) {
        kp = ipiv[ ioffipiv + k - 1 ];
        if ( kp != k ) {
          Blas1.dswap( nrhs, B, ldb, B, ldb,
            ioffb + k - 1, ioffb + kp - 1 );
        }
        k ++;
      } else {
        kp = - ipiv[ ioffipiv + k - 1 ];
        if ( k < n && kp == - ipiv[ ioffipiv + k ] ) {
          Blas1.dswap( nrhs, B, ldb, B, ldb,
            ioffb + k - 1, ioffb + kp - 1 );
        }
        k += 2;
      }
    }
  } else {
    k = 1;
    while ( k <= n ) {
      if ( ipiv[ ioffipiv + k - 1 ] > 0 ) {
        kp = ipiv[ ioffipiv + k - 1 ];
        if ( kp != k ) {
          Blas1.dswap( nrhs, B, ldb, B, ldb,
            ioffb + k - 1, ioffb + kp - 1 );
          k ++;
        }
      } else {
        kp = - ipiv[ ioffipiv + k ];
        if ( kp == - ipiv[ ioffipiv + k - 1 ] ) {
          Blas1.dswap( nrhs, B, ldb, B, ldb,
            ioffb + k, ioffb + kp - 1 );
        }
        k += 2;
      }
    }
    Blas3.dtrsm( 'L', 'L', 'N', 'U', n, nrhs, 1., A, lda, B, ldb,
      ioffa, ioffb );
    i = 1;
    while ( i <= n ) {
      if ( ipiv[ ioffipiv + i - 1 ] > 0 ) {
        Blas1.dscal( nrhs, 1. / A[ ioffa + i - 1 + ( i - 1 ) * lda ],
          B, ldb, ioffb + i - 1 );
      } else {
        akm1k = work[ ioffwork + i - 1 ];
        akm1 = A[ ioffa + i - 1 + ( i - 1 ) * lda ] / akm1k;
        ak = A[ ioffa + i + i * lda ] / akm1k;
        denom = akm1 * ak - 1.;
        for ( j = 1; j <= nrhs; j ++ ) {
          bkm1 = B[ ioffb + i - 1 + ( j - 1 ) * ldb ] / akm1k;
          bk = B[ ioffb + i + ( j - 1 ) * ldb ] / akm1k;
          B[ ioffb + i - 1 + ( j - 1 ) * ldb ] =
            ( ak * bkm1 - bk ) / denom;
          B[ ioffb + i + ( j - 1 ) * ldb ] =
            ( akm1 * bk - bkm1 ) / denom;
        }
        i ++;
      }
      i ++;
    }
    Blas3.dtrsm( 'L', 'L', 'T', 'U', n, nrhs, 1., A, lda, B, ldb,
      ioffa, ioffb );
    k = n;
    while ( k >= 1 ) {
      if ( ipiv[ ioffipiv + k - 1 ] > 0 ) {
        kp = ipiv[ ioffipiv + k - 1 ];
        if ( kp != k ) {
          Blas1.dswap( nrhs, B, ldb, B, ldb,
            ioffb + k - 1, ioffb + kp - 1 );
        }
        k --;
      } else {
        kp = - ipiv[ ioffipiv + k - 1 ];
        if ( k > 1 && kp == - ipiv[ ioffipiv + k - 2 ] ) {
          Blas1.dswap( nrhs, B, ldb, B, ldb,
            ioffb + k - 1, ioffb + kp - 1 );
        }
        k -= 2;
      }
    }
  }
  LaPack0.dsyconv( uplo, 'R', n, A, lda, ipiv, work, iinfo,
    ioffa, ioffipiv, ioffwork );
}
//************************************************************************
LaPack1.dtprfs = function( uplo, trans, diag, n, nrhs, AP, B,
ldb, X, ldx, ferr, berr, work, iwork, info) {
  throw new Error("not programmed: packed array");
}
LaPack1.ztprfs = function( uplo, trans, diag, n, nrhs, AP, B,
ldb, X, ldx, ferr, berr, work, iwork, info) {
  throw new Error("not programmed: complex packed array");
}
//************************************************************************
LaPack1.dtrtri = function( uplo, diag, n, A, lda, info,
ioffa ) {
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var nounit = ( diag.charAt(0).toUpperCase() == 'N' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( ! nounit && diag.charAt(0).toUpperCase() != 'U' ) {
    info.setValue( -2 );
  } else if ( n < 0 ) info.setValue( -3 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -5 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dtrtri', -info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  if ( nounit ) {
    for ( var i = 1; i <= n; i ++ ) {
      if ( A[ ioffa + i - 1 + ( i - 1 ) * lda ] == 0. ) {
        info.setValue( i );
        return;
      }
    }
  }
  var nb =
    LaPack0.ilaenv( 1, 'dtrtri', uplo + diag, n, -1, -1, -1 );
  if ( nb <= 1 || nb >= n ) {
    LaPack0.dtrti2( uplo, diag, n, A, lda, info, ioffa );
  } else {
    if ( upper ) {
      for ( var j = 1; j <= n; j += nb ) {
        var jb = Math.min( nb, n - j + 1 );
        Blas3.dtrmm( 'Left', 'Upper', 'No transpose', diag, j - 1, jb,
          1., A, lda, A, lda, ioffa, ioffa + ( j - 1 ) * lda );
        Blas3.dtrsm( 'Right', 'Upper', 'No transpose', diag, j - 1, jb,
          -1., A, lda, A, lda, ioffa + j - 1 + ( j - 1 ) * lda,
          ioffa + ( j - 1 ) * lda );
        LaPack0.dtrti2( 'Upper', diag, jb, A, lda, info,
          ioffa + j - 1 + ( j - 1 ) * lda );
      }
    } else {
      var nn = ( ( n - 1 ) / nb ) * nb + 1;
      for ( j = nn; j >= 1; j -= nb ) {
        jb = Math.min( nb, n - j + 1 );
        if ( j + jb <= n ) {
          Blas3.dtrmm( 'Left', 'Lower', 'No transpose', diag,
            n - j - jb + 1, jb, 1., A, lda, A, lda, 
            ioffa + j + jb - 1 + ( j + jb - 1 ) * lda,
            ioffa + j + jb - 1 + ( j - 1 ) * lda );
          Blas3.dtrsm( 'Right', 'Lower', 'No transpose', diag,
            n - j - jb + 1, jb, -1., A, lda, A, lda,
            ioffa + j - 1 + ( j - 1 ) * lda,
            ioffa + j + jb - 1 + ( j - 1 ) * lda );
        }
        LaPack0.dtrti2( 'Lower', diag, jb, A, lda, info,
          ioffa + j - 1 + ( j - 1 ) * lda );
      }
    }
  }
}
LaPack1.ztrtri = function( uplo, diag, n, A, lda, info ) {
  throw new Error("not programmed: complex array");
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function testLaPack1() {
  test_ddisna();
//test_dgebal();
//test_dgeequ();
//test_dgeequb();
//test_dgesc2();
//test_dgetc2();
//test_dgetf2();
//test_dgetrs();
//test_dgttrs();
//test_dla_geamv();
//test_dla_syamv();
//test_dlaed6();
//test_dlagtf();
//test_dlagts();
//test_dlaic1();
//test_dlange();
//test_dlangt();
//test_dlanst();
//test_dlansy();
//test_dlantr();
//test_dlanv2();
//test_dlaqge();
//test_dlaqsy();
//test_dlarf();
//test_dlarfb();
//test_dlarfg();
//test_dlarfgp();
//test_dlarrk();
//test_dlartg();
//test_dlartgp();
//test_dlascl();
//test_dlasv2();
//test_dlasy2();
//test_dlatrs();
//test_dlauum();
//test_dpoequb();
//test_dpotrf();
//test_dpstf2();
//test_dpttrs();
//test_drscl();
//test_dstebz();
/*test_dsycon(); // not programmed yet */
//test_dsyequb();
//test_dsyrfs();
//test_dsytrf();
//test_dsytrs2(); /* TRF format */
//test_dtrtri();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_ddisna() {
  document.getElementById("debug_textarea").value +=
    "testing ddisna *************" + "\n";
  var m = 4;
  var n = 3;
  var ioffd = 1;
  var ioffsep = 2;

  var d = new Array( ioffd + m );
  var sep = new Array( ioffsep + m );
  for ( var i = 0; i < m; i ++ ) d[ ioffd + i ] = i * i + 1;
  var info = new IntReference();
  LaPack1.ddisna('E',m,n,d,sep,info,ioffd,ioffsep);
  document.getElementById("debug_textarea").value +=
    "ddisna('E',...): sep = "
    + sep[ ioffsep + 0 ] + " "
    + sep[ ioffsep + 1 ] + " "
    + sep[ ioffsep + 2 ] + " "
    + sep[ ioffsep + 3 ]  + "\n";

  d = new Array( ioffd + Math.min( m, n ) );
  sep = new Array( ioffsep + Math.min( m, n ) );
  for ( i = 0; i < Math.min( m, n ); i ++ ) {
    d[ ioffd + i ] = i * i + 1;
  }
  LaPack1.ddisna('L',m,n,d,sep,info,ioffd,ioffsep);
  document.getElementById("debug_textarea").value +=
    "ddisna('L',...): sep = "
    + sep[ ioffsep + 0 ] + " "
    + sep[ ioffsep + 1 ] + " "
    + sep[ ioffsep + 2 ]  + "\n";
  LaPack1.ddisna('L',n,m,d,sep,info,ioffd,ioffsep);
  document.getElementById("debug_textarea").value +=
    "ddisna('L',...): sep = "
    + sep[ ioffsep + 0 ] + " "
    + sep[ ioffsep + 1 ] + " "
    + sep[ ioffsep + 2 ]  + "\n";

  d = new Array( ioffd + Math.min( m, n ) );
  sep = new Array( ioffsep + Math.min( m, n ) );
  for ( i = 0; i < Math.min( m, n ); i ++ ) {
    d[ ioffd + i ] = i * i + 1;
  }
  LaPack1.ddisna('R',m,n,d,sep,info,ioffd,ioffsep);
  document.getElementById("debug_textarea").value +=
    "ddisna('R',...): sep = "
    + sep[ ioffsep + 0 ] + " "
    + sep[ ioffsep + 1 ] + " "
    + sep[ ioffsep + 2 ]  + "\n";
  LaPack1.ddisna('R',n,m,d,sep,info,ioffd,ioffsep);
  document.getElementById("debug_textarea").value +=
    "ddisna('R',...): sep = "
    + sep[ ioffsep + 0 ] + " "
    + sep[ ioffsep + 1 ] + " "
    + sep[ ioffsep + 2 ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgebal() {
  document.getElementById("debug_textarea").value +=
    "testing dgebal *************" + "\n";
  var n = 5;
  var lda = 6;
  var ioffa = 1;
  var ioffscale = 2;
  var A = new Array( ioffa + lda * n );
  var scale = new Array( ioffscale + n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 0.;
    }
  }
  var k = 1;
  for ( j = 0; j < n; j ++ ) {
    for ( i = j; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = k ++;
    }
  }
  var ilo = new IntReference();
  var ihi = new IntReference();
  var info = new IntReference();
  LaPack1.dgebal('N',n,A,lda,ilo,ihi,scale,info,ioffa,ioffscale);
  document.getElementById("debug_textarea").value +=
    "dgebal('N',...): ilo,ihi = " + ilo + " " + ihi + "\n";
  document.getElementById("debug_textarea").value +=
    "scale = "
    + scale[ ioffscale + 0 ] + " "
    + scale[ ioffscale + 1 ] + " "
    + scale[ ioffscale + 2 ] + " "
    + scale[ ioffscale + 3 ] + " "
    + scale[ ioffscale + 4 ]  + "\n";
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
    + A[ ioffa + 4 + 3 * lda ] + " "
    + A[ ioffa + 0 + 4 * lda ] + " "
    + A[ ioffa + 1 + 4 * lda ] + " "
    + A[ ioffa + 2 + 4 * lda ] + " "
    + A[ ioffa + 3 + 4 * lda ] + " "
    + A[ ioffa + 4 + 4 * lda ]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 0.;
    }
  }
  k = 1;
  for ( j = 0; j < n; j ++ ) {
    for ( i = j; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = k ++;
    }
  }
  LaPack1.dgebal('P',n,A,lda,ilo,ihi,scale,info,ioffa,ioffscale);
  document.getElementById("debug_textarea").value +=
    "dgebal('P',...): ilo,ihi = " + ilo + " " + ihi + "\n";
  document.getElementById("debug_textarea").value +=
    "scale = "
    + scale[ ioffscale + 0 ] + " "
    + scale[ ioffscale + 1 ] + " "
    + scale[ ioffscale + 2 ] + " "
    + scale[ ioffscale + 3 ] + " "
    + scale[ ioffscale + 4 ]  + "\n";
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
    + A[ ioffa + 4 + 3 * lda ] + " "
    + A[ ioffa + 0 + 4 * lda ] + " "
    + A[ ioffa + 1 + 4 * lda ] + " "
    + A[ ioffa + 2 + 4 * lda ] + " "
    + A[ ioffa + 3 + 4 * lda ] + " "
    + A[ ioffa + 4 + 4 * lda ]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 0.;
    }
  }
  k = 1;
  for ( j = 0; j < n; j ++ ) {
    for ( i = j; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = k ++;
    }
  }
  LaPack1.dgebal('S',n,A,lda,ilo,ihi,scale,info,ioffa,ioffscale);
  document.getElementById("debug_textarea").value +=
    "dgebal('S',...): ilo,ihi = " + ilo + " " + ihi + "\n";
  document.getElementById("debug_textarea").value +=
    "scale = "
    + scale[ ioffscale + 0 ] + " "
    + scale[ ioffscale + 1 ] + " "
    + scale[ ioffscale + 2 ] + " "
    + scale[ ioffscale + 3 ] + " "
    + scale[ ioffscale + 4 ]  + "\n";
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
    + A[ ioffa + 4 + 3 * lda ] + " "
    + A[ ioffa + 0 + 4 * lda ] + " "
    + A[ ioffa + 1 + 4 * lda ] + " "
    + A[ ioffa + 2 + 4 * lda ] + " "
    + A[ ioffa + 3 + 4 * lda ] + " "
    + A[ ioffa + 4 + 4 * lda ]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 0.;
    }
  }
  k = 1;
  for ( j = 0; j < n; j ++ ) {
    for ( i = j; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = k ++;
    }
  }
  LaPack1.dgebal('B',n,A,lda,ilo,ihi,scale,info,ioffa,ioffscale);
  document.getElementById("debug_textarea").value +=
    "dgebal('B',...): ilo,ihi = " + ilo + " " + ihi + "\n";
  document.getElementById("debug_textarea").value +=
    "scale = "
    + scale[ ioffscale + 0 ] + " "
    + scale[ ioffscale + 1 ] + " "
    + scale[ ioffscale + 2 ] + " "
    + scale[ ioffscale + 3 ] + " "
    + scale[ ioffscale + 4 ]  + "\n";
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
    + A[ ioffa + 4 + 3 * lda ] + " "
    + A[ ioffa + 0 + 4 * lda ] + " "
    + A[ ioffa + 1 + 4 * lda ] + " "
    + A[ ioffa + 2 + 4 * lda ] + " "
    + A[ ioffa + 3 + 4 * lda ] + " "
    + A[ ioffa + 4 + 4 * lda ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgeequ() {
  document.getElementById("debug_textarea").value +=
    "testing dgeequ *************" + "\n";
  var m = 4;
  var n = 3;
  var lda = 5;
  var ioffa = 1;
  var ioffc = 2;
  var ioffr = 3;
  var A = new Array( ioffa + lda * n );
  var r = new Array( ioffr + m );
  var c = new Array( ioffc + n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 0.;
    }
  }
  var k = 1;
  for ( j = 0; j < n; j ++ ) {
    for ( i = j; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = k ++;
    }
  }
  var rowcndReference = new NumberReference();
  var colcndReference = new NumberReference();
  var amaxReference = new NumberReference();
  var info = new IntReference();
  LaPack1.dgeequ(m,n,A,lda,r,c,rowcndReference,colcndReference,
    amaxReference,info,ioffa,ioffr,ioffc);
  document.getElementById("debug_textarea").value +=
    "dgeequ(...): amax,colcnd,rowcnd = "
    + amaxReference + " "
    + colcndReference + " "
    + rowcndReference + "\n";
  document.getElementById("debug_textarea").value +=
    "c = "
    + c[ ioffc + 0 ] + " "
    + c[ ioffc + 1 ] + " "
    + c[ ioffc + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "r = "
    + r[ ioffr + 0 ] + " "
    + r[ ioffr + 1 ] + " "
    + r[ ioffr + 2 ] + " "
    + r[ ioffr + 3 ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgeequb() {
  document.getElementById("debug_textarea").value +=
    "testing dgeequb *************" + "\n";
  var m = 4;
  var n = 3;
  var lda = 5;
  var ioffa = 1;
  var ioffc = 2;
  var ioffr = 3;
  var A = new Array( ioffa + lda * n );
  var r = new Array( ioffr + m );
  var c = new Array( ioffc + n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 0.;
    }
  }
  var k = 1;
  for ( j = 0; j < n; j ++ ) {
    for ( i = j; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = k ++;
    }
  }
  var rowcndReference = new NumberReference();
  var colcndReference = new NumberReference();
  var amaxReference = new NumberReference();
  var info = new IntReference();
  LaPack1.dgeequb(m,n,A,lda,r,c,rowcndReference,colcndReference,
    amaxReference,info,ioffa,ioffr,ioffc);
  document.getElementById("debug_textarea").value +=
     "dgeequb(...): amax,colcnd,rowcnd = "
     + amaxReference + " "
     + colcndReference + " "
     + rowcndReference + "\n";
  document.getElementById("debug_textarea").value +=
    "c = "
    + c[ ioffc + 0 ] + " "
    + c[ ioffc + 1 ] + " "
    + c[ ioffc + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "r = "
    + r[ ioffr + 0 ] + " "
    + r[ ioffr + 1 ] + " "
    + r[ ioffr + 2 ] + " "
    + r[ ioffr + 3 ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgesc2() {
  document.getElementById("debug_textarea").value +=
    "testing dgesc2 *************" + "\n";
  var n = 4;
  var lda = 5;
  var ioffA = 1;
  var A = new Array( ioffA + lda * n );
  var ioffipiv = 2;
  var ioffjpiv = 3;
  var ipiv = new Array( ioffipiv + n );
  var jpiv = new Array( ioffjpiv + n );
  var info = new IntReference();
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < n - 1; j ++ ) {
    for ( i = 0; i < j; i ++ ) {
      A[ ioffA + i + j * lda ] = 0.;
    }
    A[ ioffA + j + j * lda ] = 1.;
    for ( i = j+1; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = -1.;
    }
  }
  for ( i = 0; i < n - 1; i ++ ) {
    A[ ioffA + i + j * lda ] = -1.;
  }
  A[ ioffA + j + (n - 1 ) * lda ] = 1.;
  LaPack1.dgetc2(n,A,lda,ipiv,jpiv,info,ioffA,ioffipiv,ioffjpiv );
  var ioffrhs = 4;
  var rhs = new Array( ioffrhs + n );
  rhs[ ioffrhs + 0 ] = -3.;
  rhs[ ioffrhs + 1 ] = -3.;
  rhs[ ioffrhs + 2 ] = -4.;
  rhs[ ioffrhs + 3 ] = -2.;
  var scaleReference = new NumberReference();
  LaPack1.dgesc2(n,A,lda,rhs,ipiv,jpiv,scaleReference,
    ioffA,ioffrhs,ioffipiv,ioffjpiv );
  document.getElementById("debug_textarea").value +=
    "scale = " + scaleReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "rhs = "
    + rhs[ ioffrhs + 0 ] + " "
    + rhs[ ioffrhs + 1 ] + " "
    + rhs[ ioffrhs + 2 ] + " "
    + rhs[ ioffrhs + 3 ] + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgetc2() {
  document.getElementById("debug_textarea").value +=
    "testing dgetc2 *************" + "\n";
  var n = 4;
  var ioffA = 1;
  var lda = 5;
  var A = new Array( ioffA + lda * n );
  var ioffipiv = 2;
  var ioffjpiv = 3;
  var ipiv = new Array( ioffipiv + n );
  var jpiv = new Array( ioffjpiv + n );
  var info = new IntReference();
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < n - 1; j ++ ) {
    for ( i = 0; i < j; i ++ ) {
      A[ ioffA + i + j * lda ] = 0.;
    }
    A[ ioffA + j + j * lda ] = 1.;
    for ( i = j+1; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = -1.;
    }
  }
  for ( i = 0; i < n - 1; i ++ ) {
    A[ ioffA + i + j * lda ] = -1.;
  }
  A[ ioffA + j + (n - 1 ) * lda ] = 1.;
  LaPack1.dgetc2(n,A,lda,ipiv,jpiv,info,ioffA,ioffipiv,ioffjpiv );
  document.getElementById("debug_textarea").value +=
    "ipiv = "
    + ipiv[ ioffipiv + 0 ] + " "
    + ipiv[ ioffipiv + 1 ] + " "
    + ipiv[ ioffipiv + 2 ] + "\n";
  document.getElementById("debug_textarea").value +=
    "jpiv = "
    + jpiv[ ioffjpiv + 0 ] + " "
    + jpiv[ ioffjpiv + 1 ] + " "
    + jpiv[ ioffjpiv + 2 ] + "\n";
  document.getElementById("debug_textarea").value +=
    "A = "
    + A[ ioffA + 0 + 0 * lda ] + " "
    + A[ ioffA + 1 + 0 * lda ] + " "
    + A[ ioffA + 2 + 0 * lda ] + " "
    + A[ ioffA + 3 + 0 * lda ] + " "
    + A[ ioffA + 0 + 1 * lda ] + " "
    + A[ ioffA + 1 + 1 * lda ] + " "
    + A[ ioffA + 2 + 1 * lda ] + " "
    + A[ ioffA + 3 + 1 * lda ] + " "
    + A[ ioffA + 0 + 2 * lda ] + " "
    + A[ ioffA + 1 + 2 * lda ] + " "
    + A[ ioffA + 2 + 2 * lda ] + " "
    + A[ ioffA + 3 + 2 * lda ] + " "
    + A[ ioffA + 0 + 3 * lda ] + " "
    + A[ ioffA + 1 + 3 * lda ] + " "
    + A[ ioffA + 2 + 3 * lda ] + " "
    + A[ ioffA + 3 + 3 * lda ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgetf2() {
  document.getElementById("debug_textarea").value +=
    "testing dgetf2 *************" + "\n";
  var m = 4;
  var n = 3;
  var ioffA = 1;
  var lda = 5;
  var A = new Array( ioffA + lda * n );
  var ioffipiv = 2;
  var ipiv = new Array( ioffipiv + Math.min( m, n ) );
  var info = new IntReference();
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
//        A[ ioffA + i + j * lda ] = 1. / ( 1. + i + j );
      A[ ioffA + i + j * lda ] = 1 + 2 * i + Math.pow( i + 1, j );
    }
  }
  LaPack1.dgetf2( m, n, A, lda, ipiv, info, ioffA, ioffipiv );
  document.getElementById("debug_textarea").value +=
    "dgetf2(m,n,A,lda,ipiv,info) , ipiv = "
    + ipiv[ ioffipiv + 0 ] + " "
    + ipiv[ ioffipiv + 1 ] + " "
    + ipiv[ ioffipiv + 2 ] + "\n";
  document.getElementById("debug_textarea").value +=
    "dgetf2(m,n,A,lda,ipiv,info) , A = "
    + A[ ioffA + 0 + 0 * lda ] + " "
    + A[ ioffA + 1 + 0 * lda ] + " "
    + A[ ioffA + 2 + 0 * lda ] + " "
    + A[ ioffA + 3 + 0 * lda ] + " "
    + A[ ioffA + 0 + 1 * lda ] + " "
    + A[ ioffA + 1 + 1 * lda ] + " "
    + A[ ioffA + 2 + 1 * lda ] + " "
    + A[ ioffA + 3 + 1 * lda ] + " "
    + A[ ioffA + 0 + 2 * lda ] + " "
    + A[ ioffA + 1 + 2 * lda ] + " "
    + A[ ioffA + 2 + 2 * lda ] + " "
    + A[ ioffA + 3 + 2 * lda ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgetrs() {
  document.getElementById("debug_textarea").value +=
    "testing dgetrs *************" + "\n";
  var n = 3;
  var lda = 5;
  var ioffa = 1;
  var ioffipiv = 2;
  var A = new Array( ioffa + lda * n );
  var ipiv = new Array( ioffipiv + n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + 2 * i + Math.pow( i + 1, j );
    }
  }
  var info = new IntReference();
  LaPack2.dgetrf(n,n,A,lda,ipiv,info,ioffa,ioffipiv);
  var nrhs = 1;
  var ldb = 6;
  var ioffb = 2;
  var B = new Array( ioffb + ldb * 1 );
  B[ ioffb + 0 + 0 * ldb ] = 12.;
  B[ ioffb + 1 + 0 * ldb ] = 35.;
  B[ ioffb + 2 + 0 * ldb ] = 64.;
  B[ ioffb + 3 + 0 * ldb ] = 99.;
  LaPack1.dgetrs('N',n,nrhs,A,lda,ipiv,B,ldb,info,
    ioffa,ioffipiv,ioffb);
  document.getElementById("debug_textarea").value +=
    "B = "
    + B[ ioffb + 0 + 0 * ldb ] + " "
    + B[ ioffb + 1 + 0 * ldb ] + " "
    + B[ ioffb + 2 + 0 * ldb ]  + "\n";
  B[ ioffb + 0 + 0 * ldb ] = 4.;
  B[ ioffb + 1 + 0 * ldb ] = 4.;
  B[ ioffb + 2 + 0 * ldb ] = 6.;
  LaPack1.dgetrs('T',n,nrhs,A,lda,ipiv,B,ldb,info,
    ioffa,ioffipiv,ioffb);
  document.getElementById("debug_textarea").value +=
    "B = "
    + B[ ioffb + 0 + 0 * ldb ] + " "
    + B[ ioffb + 1 + 0 * ldb ] + " "
    + B[ ioffb + 2 + 0 * ldb ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgttrs() {
  document.getElementById("debug_textarea").value +=
    "testing dgttrs *************" + "\n";
  var n = 3;
  var ioffdl = 1;
  var ioffd = 2;
  var ioffdu = 3;
  var ioffdu2 = 4;
  var ioffipiv = 5;
  var dl = new Array( ioffdl + n - 1 );
  var d = new Array( ioffd + n );
  var du = new Array( ioffdu + n - 1 );
  var du2 = new Array( ioffdu2 + n - 2 );
  var ipiv = new Array( ioffipiv + n );
  var info = new IntReference( 0 );
  for ( var i = 0; i < n; i ++ ) d[ ioffd + i ] = 2.;
  for ( i = 0; i < n - 1; i ++ ) {
    dl[ ioffdl + i ] = -1.;
    du[ ioffdu + i ] = -1.;
  }
  LaPack0.dgttrf( n, dl, d, du, du2, ipiv, info,
    ioffdl, ioffd, ioffdu, ioffdu2, ioffipiv );
  var nrhs = 1;
  var ldb = 3;
  var ioffB = 6;
  var B = new Array( ioffB + ldb * nrhs );

  B[ ioffB + 0 ] = 0.;
  B[ ioffB + 1 ] = 0.;
  B[ ioffB + 2 ] = 4.;
  LaPack1.dgttrs( 'N', n, nrhs, dl, d, du, du2, ipiv, B, ldb, info,
    ioffdl, ioffd, ioffdu, ioffdu2, ioffipiv, ioffB );
  document.getElementById("debug_textarea").value +=
    "dgttrf, B = "
    + B[ ioffB + 0 ] + " "
    + B[ ioffB + 1 ] + " "
    + B[ ioffB + 2 ]  + "\n";

  B[ ioffB + 0 ] = 0.;
  B[ ioffB + 1 ] = 0.;
  B[ ioffB + 2 ] = 4.;
  LaPack1.dgttrs( 'T', n, nrhs, dl, d, du, du2, ipiv, B, ldb, info,
    ioffdl, ioffd, ioffdu, ioffdu2, ioffipiv, ioffB );
  document.getElementById("debug_textarea").value +=
    "dgttrf, B = "
    + B[ ioffB + 0 ] + " "
    + B[ ioffB + 1 ] + " "
    + B[ ioffB + 2 ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dla_geamv() {
  document.getElementById("debug_textarea").value +=
    "testing dla_geamv *************" + "\n";
  var BLAS_NO_TRANS = LaPack0.ilatrans('N');
  var BLAS_TRANS = LaPack0.ilatrans('T');
  var BLAS_CONJ_TRANS = LaPack0.ilatrans('C');
  var ioffA = 1;
  var ioffx = 2;
  var ioffy = 3;
  var A = new Array( ioffA + 15 );
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 8 );
  var m = 4;
  var n = 3;
  var lda = 5;
  for (var j = 0; j < n; j ++ ) {
    for (var i = 0; i < m; i ++ ) {
      A[ ioffA + i + j * lda ] = Number( i + j );
    }
  }
  var incy = 2;
  var incx = 2;
  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  LaPack1.dla_geamv(BLAS_NO_TRANS,m,n,0.,A,lda,x,incx,1.,y,incy,
    ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
    "dla_geamv(N,m,n,0.,A,lda,x,incx,1.,y,incy): y = "
    + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4] + " "
    + y[ioffy + 6] + "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  LaPack1.dla_geamv(BLAS_NO_TRANS,m,n,0.,A,lda,x,incx,2.,y,incy,
    ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
    "dla_geamv(N,m,n,0.,A,lda,x,incx,2.,y,incy): y = "
    + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4] + " "
    + y[ioffy + 6] + "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  LaPack1.dla_geamv(BLAS_NO_TRANS,m,n,3.,A,lda,x,incx,1.,y,incy,
    ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
    "dla_geamv(N,m,n,3.,A,lda,x,incx,1.,y,incy): y = "
    + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4] + " "
    + y[ioffy + 6] + "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  LaPack1.dla_geamv(BLAS_NO_TRANS,m,n,3.,A,lda,x,incx,2.,y,incy,
    ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
    "dla_geamv(N,m,n,3.,A,lda,x,incx,2.,y,incy): y = "
    + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4] + " "
    + y[ioffy + 6] + "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  LaPack1.dla_geamv(BLAS_TRANS,m,n,0.,A,lda,y,incy,1.,x,incx,
    ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
    "dla_geamv(T,m,n,0.,A,lda,y,incy,1.,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  LaPack1.dla_geamv(BLAS_TRANS,m,n,0.,A,lda,y,incy,2.,x,incx,
    ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
    "dla_geamv(T,m,n,0.,A,lda,y,incy,2.,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  LaPack1.dla_geamv(BLAS_TRANS,m,n,3.,A,lda,y,incy,1.,x,incx,
    ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
    "dla_geamv(T,m,n,3.,A,lda,y,incy,1.,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  LaPack1.dla_geamv(BLAS_TRANS,m,n,3.,A,lda,y,incy,2.,x,incx,
    ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
    "dla_geamv(T,m,n,3.,A,lda,y,incy,2.,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dla_syamv() {
  document.getElementById("debug_textarea").value +=
    "testing dla_syamv *************" + "\n";
  var BLAS_UPPER = 121;
  var BLAS_LOWER = 122;
  var ioffA = 1;
  var ioffx = 2;
  var ioffy = 3;
  var A = new Array( ioffA + 12 );
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 6 );
  var n = 3;
  var lda = 4;
  for (var j = 0; j < n; j ++ ) {
    for (var i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = Number( i + j*2 );
    }
  }
  var incy = 2;
  var incx = 2;
  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  LaPack1.dla_syamv(BLAS_UPPER,n,0.,A,lda,x,incx,1.,y,incy,
    ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
    "dla_syamv(U,n,0.,A,lda,x,incx,1.,y,incy): y = "
    + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  LaPack1.dla_syamv(BLAS_UPPER,n,0.,A,lda,x,incx,2.,y,incy,
    ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
    "dla_syamv(U,n,0.,A,lda,x,incx,2.,y,incy): y = "
    + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  LaPack1.dla_syamv(BLAS_UPPER,n,3.,A,lda,x,incx,1.,y,incy,
    ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
    "dla_syamv(U,n,3.,A,lda,x,incx,1.,y,incy): y = "
    + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  LaPack1.dla_syamv(BLAS_UPPER,n,3.,A,lda,x,incx,2.,y,incy,
    ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
    "dla_syamv(U,n,3.,A,lda,x,incx,2.,y,incy): y = "
    + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  LaPack1.dla_syamv(BLAS_LOWER,n,0.,A,lda,y,incy,1.,x,incx,
    ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
    "dla_syamv(L,n,0.,A,lda,y,incy,1.,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  LaPack1.dla_syamv(BLAS_LOWER,n,0.,A,lda,y,incy,2.,x,incx,
    ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
    "dla_syamv(L,n,0.,A,lda,y,incy,2.,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  LaPack1.dla_syamv(BLAS_LOWER,n,3.,A,lda,y,incy,1.,x,incx,
    ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
    "dla_syamv(L,n,3.,A,lda,y,incy,1.,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  LaPack1.dla_syamv(BLAS_LOWER,n,3.,A,lda,y,incy,2.,x,incx,
    ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
    "dla_syamv(L,n,3.,A,lda,y,incy,2.,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlaed6() {
  document.getElementById("debug_textarea").value +=
    "testing dlaed6 *************" + "\n";
  var ioffd = 1;
  var ioffz = 2;
  var d = new Array(ioffd + 3);
  var z = new Array(ioffz + 3);
  d[ ioffd + 0] = -3.;
  d[ ioffd + 1] = -2.;
  d[ ioffd + 2] = 1.;
  z[ ioffz + 0] = 4.;
  z[ ioffz + 1] = 5.;
  z[ ioffz + 2] = 6.;
  var tau = -1.;
  var rho = - z[ ioffz + 0] / ( d[ ioffd + 0] - tau )
    - z[ ioffz + 1] / ( d[ ioffd + 1] - tau )
    - z[ ioffz + 2] / ( d[ ioffd + 2] - tau );
  var finit = rho + z[ ioffz + 0] / d[ ioffd + 0]
    + z[ ioffz + 1] / d[ ioffd + 1] + z[ ioffz + 2] / d[ ioffd + 2];
  var tauref = new NumberReference();
  var info = new IntReference();
  LaPack1.dlaed6(1,true,rho,d,z,finit,tauref,info,ioffd,ioffz);
  document.getElementById("debug_textarea").value +=
    "dlaed6(1,true,...): info,tau = " + info.getValue() + " "
    + tauref.getValue() + "\n";
  LaPack1.dlaed6(2,true,rho,d,z,finit,tauref,info,ioffd,ioffz);
  document.getElementById("debug_textarea").value +=
    "dlaed6(2,true,...): info,tau = " + info.getValue() + " "
    + tauref.getValue() + "\n";

  tau = -1.9;
  rho = - z[ ioffz + 0] / ( d[ ioffd + 0] - tau )
    - z[ ioffz + 1] / ( d[ ioffd + 1] - tau )
    - z[ ioffz + 2] / ( d[ ioffd + 2] - tau );
  finit = rho + z[ ioffz + 0] / d[ ioffd + 0]
    + z[ ioffz + 1] / d[ ioffd + 1] + z[ ioffz + 2] / d[ ioffd + 2];
  LaPack1.dlaed6(1,true,rho,d,z,finit,tauref,info,ioffd,ioffz);
  document.getElementById("debug_textarea").value +=
    "dlaed6(1,true,...): info,tau = " + info.getValue() + " "
    + tauref.getValue() + "\n";
  LaPack1.dlaed6(2,true,rho,d,z,finit,tauref,info,ioffd,ioffz);
  document.getElementById("debug_textarea").value +=
    "dlaed6(2,true,...): info,tau = " + info.getValue() + " "
    + tauref.getValue() + "\n";

  d[ ioffd + 0] = -1.;
  d[ ioffd + 1] = 2.;
  d[ ioffd + 2] = 3.;
  tau = 1.;
  rho = - z[ ioffz + 0] / ( d[ ioffd + 0] - tau )
    - z[ ioffz + 1] / ( d[ ioffd + 1] - tau )
    - z[ ioffz + 2] / ( d[ ioffd + 2] - tau );
  finit = rho + z[ ioffz + 0] / d[ ioffd + 0]
    + z[ ioffz + 1] / d[ ioffd + 1] + z[ ioffz + 2] / d[ ioffd + 2];
  LaPack1.dlaed6(1,false,rho,d,z,finit,tauref,info,ioffd,ioffz);
  document.getElementById("debug_textarea").value +=
    "dlaed6(1,false,...): info,tau = " + info.getValue() + " "
    + tauref.getValue() + "\n";
  LaPack1.dlaed6(2,false,rho,d,z,finit,tauref,info,ioffd,ioffz);
  document.getElementById("debug_textarea").value +=
    "dlaed6(2,false,...): info,tau = " + info.getValue() + " "
    + tauref.getValue() + "\n";

  tau = 1.9;
  rho = - z[ ioffz + 0] / ( d[ ioffd + 0] - tau )
    - z[ ioffz + 1] / ( d[ ioffd + 1] - tau )
    - z[ ioffz + 2] / ( d[ ioffd + 2] - tau );
  finit = rho + z[ ioffz + 0] / d[ ioffd + 0]
    + z[ ioffz + 1] / d[ ioffd + 1] + z[ ioffz + 2] / d[ ioffd + 2];
  LaPack1.dlaed6(1,false,rho,d,z,finit,tauref,info,ioffd,ioffz);
  document.getElementById("debug_textarea").value +=
    "dlaed6(1,false,...): info,tau = " + info.getValue() + " "
    + tauref.getValue() + "\n";
  LaPack1.dlaed6(2,false,rho,d,z,finit,tauref,info,ioffd,ioffz);
  document.getElementById("debug_textarea").value +=
    "dlaed6(2,false,...): info,tau = " + info.getValue() + " "
    + tauref.getValue() + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlagtf() {
  document.getElementById("debug_textarea").value +=
    "testing dlagtf *************" + "\n";
  var n = 4;
  var ioffa = 1;
  var ioffb = 2;
  var ioffc = 3;
  var ioffd = 4;
  var ioffin2 = 5;
  var a = new Array( ioffa + n );
  var b = new Array( ioffb + n - 1 );
  var c = new Array( ioffc + n - 1 );
  var d = new Array( ioffd + n - 2 );
  var in2 = new Array( ioffin2 + n );
  for ( var i = 0; i < n; i ++ ) a[ ioffa + i ] = 2.;
  for ( i=0; i < n - 1; i ++ ) {
    b[ ioffb + i ] = -1.;
    c[ ioffc + i ] = -1.;
  }
  var lambda = -1.;
  var tol = 0.;
  var info = new IntReference();
  LaPack1.dlagtf(n,a,lambda,b,c,tol,d,in2,info,
    ioffa,ioffb,ioffc,ioffd,ioffin2);
  document.getElementById("debug_textarea").value +=
    "a = "
    + a[ ioffa + 0 ] + " "
    + a[ ioffa + 1 ] + " "
    + a[ ioffa + 2 ] + " "
    + a[ ioffa + 3 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "b = "
    + b[ ioffb + 0 ] + " "
    + b[ ioffb + 1 ] + " "
    + b[ ioffb + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "c = "
    + c[ ioffc + 0 ] + " "
    + c[ ioffc + 1 ] + " "
    + c[ ioffc + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "d = "
    + d[ ioffd + 0 ] + " "
    + d[ ioffd + 1 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "in2 = "
    + in2[ ioffin2 + 0 ] + " "
    + in2[ ioffin2 + 1 ] + " "
    + in2[ ioffin2 + 2 ] + " "
    + in2[ ioffin2 + 3 ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlagts() {
  document.getElementById("debug_textarea").value +=
    "testing dlagts *************" + "\n";
  var n = 5;
  var ioffa = 1;
  var ioffb = 2;
  var ioffc = 3;
  var ioffd = 4;
  var ioffin2 = 5;
  var a = new Array( ioffa + n );
  var b = new Array( ioffb + n - 1 );
  var c = new Array( ioffc + n - 1 );
  var d = new Array( ioffd + n - 2 );
  var in2 = new Array( ioffin2 + n );
  for ( var i = 0; i < n; i ++ ) a[ ioffa + i ] = 2.;
  for ( i=0; i < n - 1; i ++ ) {
    b[ ioffb + i ] = -1.;
    c[ ioffc + i ] = -1.;
  }
  var lambda = -1.;
  var tol = 0.;
  var info = new IntReference();
  LaPack1.dlagtf(n,a,lambda,b,c,tol,d,in2,info,
    ioffa,ioffb,ioffc,ioffd,ioffin2);

  var ioffy = 6;
  var y = new Array( ioffy + n );
  for ( i = 0; i < n; i ++ ) y[ ioffy + i ] = i + 1;
  var tolref = new NumberReference(0.);
  LaPack1.dlagts(1,n,a,b,c,d,in2,y,tolref,info,
    ioffa,ioffb,ioffc,ioffd,ioffin2,ioffy);
  document.getElementById("debug_textarea").value +=
    "dlagts(1,...): info,tol = " + info.getValue() + " "
    + tolref.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "y = "
    + y[ ioffy + 0 ] + " "
    + y[ ioffy + 1 ] + " "
    + y[ ioffy + 2 ] + " "
    + y[ ioffy + 3 ] + " "
    + y[ ioffy + 4 ]  + "\n";
  for ( i = 0; i < n; i ++ ) y[ ioffy + i ] = i + 1;
  tolref.setValue( 0. );
  LaPack1.dlagts(-1,n,a,b,c,d,in2,y,tolref,info,
    ioffa,ioffb,ioffc,ioffd,ioffin2,ioffy);
  document.getElementById("debug_textarea").value +=
    "dlagts(-1,...): info,tol = " + info.getValue() + " "
    + tolref.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "y = "
    + y[ ioffy + 0 ] + " "
    + y[ ioffy + 1 ] + " "
    + y[ ioffy + 2 ] + " "
    + y[ ioffy + 3 ] + " "
    + y[ ioffy + 4 ]  + "\n";

  for ( i = 0; i < n; i ++ ) y[ ioffy + i ] = i + 1;
  tolref.setValue( 0. );
  LaPack1.dlagts(2,n,a,b,c,d,in2,y,tolref,info,
    ioffa,ioffb,ioffc,ioffd,ioffin2,ioffy);
  document.getElementById("debug_textarea").value +=
    "dlagts(2,...): info,tol = " + info.getValue() + " "
    + tolref.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "y = "
    + y[ ioffy + 0 ] + " "
    + y[ ioffy + 1 ] + " "
    + y[ ioffy + 2 ] + " "
    + y[ ioffy + 3 ] + " "
    + y[ ioffy + 4 ]  + "\n";

  for ( i = 0; i < n; i ++ ) y[ ioffy + i ] = i + 1;
  tolref.setValue( 0. );
  LaPack1.dlagts(-2,n,a,b,c,d,in2,y,tolref,info,
    ioffa,ioffb,ioffc,ioffd,ioffin2,ioffy);
  document.getElementById("debug_textarea").value +=
    "dlagts(-2,...): info,tol = " + info.getValue() + " "
    + tolref.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "y = "
    + y[ ioffy + 0 ] + " "
    + y[ ioffy + 1 ] + " "
    + y[ ioffy + 2 ] + " "
    + y[ ioffy + 3 ] + " "
    + y[ ioffy + 4 ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlaic1() {
  document.getElementById("debug_textarea").value +=
    "testing dlaic1 *************" + "\n";
  var j = 5;
  var ioffx = 1;
  var ioffw = 2;
  var x = new Array( ioffx + j );
  var w = new Array( ioffw + j );
  for ( var i = 0; i < j; i ++ ) {
    x[ ioffx + i ] = 1;
    w[ ioffw + i ] = 2 - i;
  }
  var sest = 0.;
  var gamma = 0.;
  var sestprReference = new NumberReference();
  var sReference = new NumberReference();
  var cReference = new NumberReference();
//    sest = 0, s1 = 0
  LaPack1.dlaic1(1,j,x,sest,w,gamma,sestprReference,sReference,
    cReference,ioffx,ioffw);
  document.getElementById("debug_textarea").value +=
    "dlaic1(1,...): sestpr,s,c = " + sestprReference.getValue() + " "
    + sReference.getValue() + " " + cReference.getValue() + "\n";
  LaPack1.dlaic1(2,j,x,sest,w,gamma,sestprReference,sReference,
    cReference,ioffx,ioffw);
  document.getElementById("debug_textarea").value +=
    "dlaic1(2,...): sestpr,s,c = " + sestprReference.getValue() + " "
    + sReference.getValue() + " " + cReference.getValue() + "\n";

//    sest = 0, s1 != 0
  gamma = 1.;
  LaPack1.dlaic1(1,j,x,sest,w,gamma,sestprReference,sReference,
    cReference,ioffx,ioffw);
  document.getElementById("debug_textarea").value +=
    "dlaic1(1,...): sestpr,s,c = " + sestprReference.getValue() + " "
    + sReference.getValue() + " " + cReference.getValue() + "\n";
  LaPack1.dlaic1(2,j,x,sest,w,gamma,sestprReference,sReference,
    cReference,ioffx,ioffw);
  document.getElementById("debug_textarea").value +=
    "dlaic1(2,...): sestpr,s,c = " + sestprReference.getValue() + " "
    + sReference.getValue() + " " + cReference.getValue() + "\n";

//    sest != 0, | gamma | <= eps * | sest |
  sest = 1.;
  gamma = 0.;
  LaPack1.dlaic1(1,j,x,sest,w,gamma,sestprReference,sReference,
    cReference,ioffx,ioffw);
  document.getElementById("debug_textarea").value +=
    "dlaic1(1,...): sestpr,s,c = " + sestprReference.getValue() + " "
    + sReference.getValue() + " " + cReference.getValue() + "\n";
  LaPack1.dlaic1(2,j,x,sest,w,gamma,sestprReference,sReference,
    cReference,ioffx,ioffw);
  document.getElementById("debug_textarea").value +=
    "dlaic1(2,...): sestpr,s,c = " + sestprReference.getValue() + " "
    + sReference.getValue() + " " + cReference.getValue() + "\n";

//    sest != 0, | gamma | > eps * | sest |, | alpha | <= eps * | sest |
//    | gamma | <= | sest |
  sest = 1.;
  gamma = 1.;
  LaPack1.dlaic1(1,j,x,sest,w,gamma,sestprReference,sReference,
    cReference,ioffx,ioffw);
  document.getElementById("debug_textarea").value +=
    "dlaic1(1,...): sestpr,s,c = " + sestprReference.getValue() + " "
    + sReference.getValue() + " " + cReference.getValue() + "\n";
  LaPack1.dlaic1(2,j,x,sest,w,gamma,sestprReference,sReference,
    cReference,ioffx,ioffw);
  document.getElementById("debug_textarea").value +=
    "dlaic1(2,...): sestpr,s,c = " + sestprReference.getValue() + " "
    + sReference.getValue() + " " + cReference.getValue() + "\n";

//    sest != 0, | gamma | > eps * | sest |, | alpha | <= eps * | sest |
//    | gamma | > | sest |
  sest = 1.;
  gamma = 2.;
  LaPack1.dlaic1(1,j,x,sest,w,gamma,sestprReference,sReference,
    cReference,ioffx,ioffw);
  document.getElementById("debug_textarea").value +=
    "dlaic1(1,...): sestpr,s,c = " + sestprReference.getValue() + " "
    + sReference.getValue() + " " + cReference.getValue() + "\n";
  LaPack1.dlaic1(2,j,x,sest,w,gamma,sestprReference,sReference,
    cReference,ioffx,ioffw);
  document.getElementById("debug_textarea").value +=
    "dlaic1(2,...): sestpr,s,c = " + sestprReference.getValue() + " "
    + sReference.getValue() + " " + cReference.getValue() + "\n";

//    sest != 0, | gamma | > eps * | sest |, | alpha | > eps * | sest |
//    | sest | <= eps * | alpha | or | sest | <= eps * | gamma |,
//    | gamma | <= | alpha |
  for ( i = 0; i < j; i ++ ) w[ ioffw + i ] = 2 + i;
  sest = 1.e-16;
  gamma = 1.;
  LaPack1.dlaic1(1,j,x,sest,w,gamma,sestprReference,sReference,
    cReference,ioffx,ioffw);
  document.getElementById("debug_textarea").value +=
    "dlaic1(1,...): sestpr,s,c = " + sestprReference.getValue() + " "
    + sReference.getValue() + " " + cReference.getValue() + "\n";
  LaPack1.dlaic1(2,j,x,sest,w,gamma,sestprReference,sReference,
    cReference,ioffx,ioffw);
  document.getElementById("debug_textarea").value +=
    "dlaic1(2,...): sestpr,s,c = " + sestprReference.getValue() + " "
    + sReference.getValue() + " " + cReference.getValue() + "\n";

//    sest != 0, | gamma | > eps * | sest |, | alpha | > eps * | sest |
//    | sest | <= eps * | alpha | or | sest | <= eps * | gamma |,
//    | gamma | > | alpha |
  gamma = 21.;
  LaPack1.dlaic1(1,j,x,sest,w,gamma,sestprReference,sReference,
    cReference,ioffx,ioffw);
  document.getElementById("debug_textarea").value +=
    "dlaic1(1,...): sestpr,s,c = " + sestprReference.getValue() + " "
    + sReference.getValue() + " " + cReference.getValue() + "\n";
  LaPack1.dlaic1(2,j,x,sest,w,gamma,sestprReference,sReference,
    cReference,ioffx,ioffw);
  document.getElementById("debug_textarea").value +=
    "dlaic1(2,...): sestpr,s,c = " + sestprReference.getValue() + " "
    + sReference.getValue() + " " + cReference.getValue() + "\n";

//    sest != 0, | gamma | > eps * | sest |, | alpha | > eps * | sest |
//    | sest | > eps * | alpha | and | sest | > eps * | gamma |
  sest = 1.;
  LaPack1.dlaic1(1,j,x,sest,w,gamma,sestprReference,sReference,cReference,
    ioffx, ioffw);
  document.getElementById("debug_textarea").value +=
    "dlaic1(1,...): sestpr,s,cReference = "
    + sestprReference.getValue() + " "
    + sReference.getValue() + " " + cReference.getValue() + "\n";
  LaPack1.dlaic1(2,j,x,sest,w,gamma,sestprReference,sReference,cReference,
    ioffx,ioffw);
  document.getElementById("debug_textarea").value +=
    "dlaic1(2,...): sestpr,s,cReference = "
    + sestprReference.getValue() + " "
    + sReference.getValue() + " " + cReference.getValue() + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlange() {
  document.getElementById("debug_textarea").value +=
    "testing dlange *************" + "\n";
  var lda = 4;
  var m = 3;
  var n = 2;
  var ioffA = 1;
  var ioffwork = 2;
  var A = new Array( lda * n );
  var work = new Array( m );
  for ( var j = 0; j <= n - 1; j ++ ) {
    for ( var i = 0; i <= m - 1; i ++ ) {
      A[ ioffA + i + j * lda ] = i + j;
    }
  }
  var ret = LaPack1.dlange('M',m,n,A,lda,work,ioffA,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlange('M',m,n,A,lda,work) = " + ret + "\n";
  ret = LaPack1.dlange('O',m,n,A,lda,work,ioffA,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlange('O',m,n,A,lda,work) = " + ret + "\n";
  ret = LaPack1.dlange('I',m,n,A,lda,work,ioffA,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlange('I',m,n,A,lda,work) = " + ret + "\n";
  ret = LaPack1.dlange('F',m,n,A,lda,work,ioffA,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlange('F',m,n,A,lda,work) = " + ret + "\n";
  ret = LaPack1.dlange('?',m,n,A,lda,work,ioffA,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlange('?',m,n,A,lda,work) = " + ret + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlangt() {
  document.getElementById("debug_textarea").value +=
    "testing dlangt *************" + "\n";
  var n = 3;
  var ioffd = 1;
  var ioffdl = 2;
  var ioffdu = 3;
  var d = new Array( ioffd + n );
  var dl = new Array( ioffdl + n - 1 );
  var du = new Array( ioffdu + n - 1 );
  for ( var i = 0; i < n; i ++ ) d[ ioffd + i ] = i;
  for ( i = 0; i < n - 1; i ++ ) {
    dl[ ioffdl + i ] = i + 4;
    du[ ioffdu + i ] = i + 6;
  }
  du[ ioffdu + 1 ] = 8.;
  var anorm = LaPack1.dlangt('M',n,dl,d,du,ioffdl,ioffd,ioffdu);
  document.getElementById("debug_textarea").value +=
    "dlangt('M',...) = " + anorm + "\n";
  anorm = LaPack1.dlangt('O',n,dl,d,du,ioffdl,ioffd,ioffdu);
  document.getElementById("debug_textarea").value +=
    "dlangt('O',...) = " + anorm + "\n";
  anorm = LaPack1.dlangt('I',n,dl,d,du,ioffdl,ioffd,ioffdu);
  document.getElementById("debug_textarea").value +=
    "dlangt('I',...) = " + anorm + "\n";
  anorm = LaPack1.dlangt('F',n,dl,d,du,ioffdl,ioffd,ioffdu);
  document.getElementById("debug_textarea").value +=
    "dlangt('F',...) = " + anorm + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlanst() {
  document.getElementById("debug_textarea").value +=
    "testing dlanst *************" + "\n";
  var n = 3;
  var ioffd = 1;
  var ioffe = 3;
  var d = new Array( ioffd + n );
  var e = new Array( ioffe + n - 1 );
  for ( var i = 0; i < n; i ++ ) d[ ioffd + i ] = i;
  for ( i = 0; i < n - 1; i ++ ) {
    e[ ioffe + i ] = i + 4;
  }
  var anorm = LaPack1.dlanst('M',n,d,e,ioffd,ioffe);
  document.getElementById("debug_textarea").value +=
    "dlanst('M',...) = " + anorm + "\n";
  anorm = LaPack1.dlanst('O',n,d,e,ioffd,ioffe);
  document.getElementById("debug_textarea").value +=
    "dlanst('O',...) = " + anorm + "\n";
  anorm = LaPack1.dlanst('I',n,d,e,ioffd,ioffe);
  document.getElementById("debug_textarea").value +=
    "dlanst('I',...) = " + anorm + "\n";
  anorm = LaPack1.dlanst('F',n,d,e,ioffd,ioffe);
  document.getElementById("debug_textarea").value +=
    "dlanst('F',...) = " + anorm + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlansy() {
  document.getElementById("debug_textarea").value +=
    "testing dlansy *************" + "\n";
  var n = 3;
  var lda = 4;
  var ioffa = 1;
  var ioffwork = 2;
  var A = new Array( ioffa + lda * n );
  var work = new Array( ioffwork + n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = i + 2 * j;
    }
  }
  var anorm =
    LaPack1.dlansy('M','U',n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlansy('M','U',...) = " + anorm + "\n";
  anorm=LaPack1.dlansy('M','L',n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlansy('M','L',...) = " + anorm + "\n";
  anorm=LaPack1.dlansy('O','U',n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlansy('O','U',...) = " + anorm + "\n";
  anorm=LaPack1.dlansy('O','L',n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlansy('O','L',...) = " + anorm + "\n";
  anorm=LaPack1.dlansy('I','U',n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlansy('I','U',...) = " + anorm + "\n";
  anorm=LaPack1.dlansy('I','L',n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlansy('I','L',...) = " + anorm + "\n";
  anorm=LaPack1.dlansy('F','U',n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlansy('F','U',...) = " + anorm + "\n";
  anorm=LaPack1.dlansy('F','L',n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlansy('F','L',...) = " + anorm + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlantr() {
  document.getElementById("debug_textarea").value +=
    "testing dlantr *************" + "\n";
  var m = 4;
  var n = 3;
  var lda = 5;
  var ldb = 4;
  var ioffa = 1;
  var ioffb = 2;
  var ioffwork = 3;
  var A = new Array( ioffa + lda * n );
  var B = new Array( ioffb + ldb * m );
  var work = new Array( ioffwork + Math.max( m, n ) );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = i + 2 * j;
      B[ ioffb + j + i * ldb ] = 2 * i - j;
    }
  }
  var anorm =
    LaPack1.dlantr('M','U','N',m,n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('M','U','N',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('M','U','U',m,n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('M','U','U',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('M','U','N',n,m,B,ldb,work,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('M','U','N',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('M','U','U',n,m,B,ldb,work,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('M','U','U',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('M','L','N',m,n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('M','L','N',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('M','L','U',m,n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('M','L','U',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('M','L','N',n,m,B,ldb,work,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('M','L','N',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('M','L','U',n,m,B,ldb,work,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('M','L','U',...) = " + anorm + "\n";

  anorm = LaPack1.dlantr('O','U','N',m,n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('O','U','N',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('O','U','U',m,n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('O','U','U',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('O','U','N',n,m,B,ldb,work,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('O','U','N',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('O','U','U',n,m,B,ldb,work,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('O','U','U',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('O','L','N',m,n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('O','L','N',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('O','L','U',m,n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('O','L','U',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('O','L','N',n,m,B,ldb,work,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('O','L','N',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('O','L','U',n,m,B,ldb,work,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('O','L','U',...) = " + anorm + "\n";

  anorm = LaPack1.dlantr('I','U','N',m,n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('I','U','N',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('I','U','U',m,n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('I','U','U',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('I','U','N',n,m,B,ldb,work,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('I','U','N',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('I','U','U',n,m,B,ldb,work,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('I','U','U',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('I','L','N',m,n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('I','L','N',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('I','L','U',m,n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('I','L','U',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('I','L','N',n,m,B,ldb,work,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('I','L','N',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('I','L','U',n,m,B,ldb,work,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('I','L','U',...) = " + anorm + "\n";

  anorm = LaPack1.dlantr('F','U','N',m,n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('F','U','N',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('F','U','U',m,n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('F','U','U',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('F','U','N',n,m,B,ldb,work,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('F','U','N',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('F','U','U',n,m,B,ldb,work,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('F','U','U',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('F','L','N',m,n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('F','L','N',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('F','L','U',m,n,A,lda,work,ioffa,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('F','L','U',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('F','L','N',n,m,B,ldb,work,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('F','L','N',...) = " + anorm + "\n";
  anorm = LaPack1.dlantr('F','L','U',n,m,B,ldb,work,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlantr('F','L','U',...) = " + anorm + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlanv2() {
  document.getElementById("debug_textarea").value +=
    "testing dlanv2 *************" + "\n";
  var aReference = new NumberReference( 1. );
  var dReference = new NumberReference( 2. );
  var bReference = new NumberReference( 0. );
  var cReference = new NumberReference( 0. );
  var rt1rReference = new NumberReference();
  var rt1iReference = new NumberReference();
  var rt2rReference = new NumberReference();
  var rt2iReference = new NumberReference();
  var csReference = new NumberReference();
  var snReference = new NumberReference();

//    c = 0
  LaPack1.dlanv2(aReference,bReference,cReference,dReference,
    rt1rReference,rt1iReference,rt2rReference,rt2iReference,csReference,
    snReference);
  document.getElementById("debug_textarea").value +=
    "dlanv2: a,b,c,d = " + aReference.getValue() + " "
    + bReference.getValue() + " "
    + cReference.getValue() + " " + dReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "rt1r,rt1i,rt2r,rt2i = " + rt1rReference.getValue() + " "
    + rt1iReference.getValue() + " "
    + rt2rReference.getValue() + " " + rt2iReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "cs,sn = " + csReference.getValue() + " "
    + snReference.getValue() + "\n";
  
//    c != 0, b = 0
  aReference.setValue( 1. );
  dReference.setValue( 2. );
  bReference.setValue( 0. );
  cReference.setValue( -1. );
  LaPack1.dlanv2(aReference,bReference,cReference,dReference,
    rt1rReference,rt1iReference,rt2rReference,rt2iReference,csReference,
    snReference);
  document.getElementById("debug_textarea").value +=
    "dlanv2: a,b,c,d = " + aReference.getValue() + " "
    + bReference.getValue() + " "
    + cReference.getValue() + " " + dReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "rt1r,rt1i,rt2r,rt2i = " + rt1rReference.getValue() + " "
    + rt1iReference.getValue() + " "
    + rt2rReference.getValue() + " " + rt2iReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "cs,sn = " + csReference.getValue() + " "
    + snReference.getValue() + "\n";
  
//    c != 0, b != 0, a = d & b>=0 == c >= 0
  aReference.setValue( 2. );
  dReference.setValue( 2. );
  bReference.setValue( 3. );
  cReference.setValue( 1. );
  LaPack1.dlanv2(aReference,bReference,cReference,dReference,
    rt1rReference,rt1iReference,rt2rReference,rt2iReference,csReference,
    snReference);
  document.getElementById("debug_textarea").value +=
    "dlanv2: a,b,c,d = " + aReference.getValue() + " "
    + bReference.getValue() + " "
    + cReference.getValue() + " " + dReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "rt1r,rt1i,rt2r,rt2i = " + rt1rReference.getValue() + " "
    + rt1iReference.getValue() + " "
    + rt2rReference.getValue() + " " + rt2iReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "cs,sn = " + csReference.getValue() + " "
    + snReference.getValue() + "\n";
  
//    c != 0, b != 0, a != d and b>=0 == c >= 0
  aReference.setValue( 1. );
  dReference.setValue( 2. );
  bReference.setValue( 3. );
  cReference.setValue( 1. );
  LaPack1.dlanv2(aReference,bReference,cReference,dReference,
    rt1rReference,rt1iReference,rt2rReference,rt2iReference,csReference,
      snReference);
  document.getElementById("debug_textarea").value +=
    "dlanv2: a,b,c,d = " + aReference.getValue() + " "
    + bReference.getValue() + " "
    + cReference.getValue() + " " + dReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "rt1r,rt1i,rt2r,rt2i = " + rt1rReference.getValue() + " "
    + rt1iReference.getValue() + " "
    + rt2rReference.getValue() + " " + rt2iReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "cs,sn = " + csReference.getValue() + " "
    + snReference.getValue() + "\n";
  
//    c != 0, b != 0, a == d and b>=0 != c >= 0
  aReference.setValue( 1. );
  dReference.setValue( 1. );
  bReference.setValue( 3. );
  cReference.setValue( -1. );
  LaPack1.dlanv2(aReference,bReference,cReference,dReference,
    rt1rReference,rt1iReference,rt2rReference,rt2iReference,csReference,
    snReference);
  document.getElementById("debug_textarea").value +=
    "dlanv2: a,b,c,d = " + aReference.getValue() + " "
    + bReference.getValue() + " "
    + cReference.getValue() + " " + dReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "rt1r,rt1i,rt2r,rt2i = " + rt1rReference.getValue() + " "
    + rt1iReference.getValue() + " "
    + rt2rReference.getValue() + " " + rt2iReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "cs,sn = " + csReference.getValue() + " "
    + snReference.getValue() + "\n";
  
//    c != 0, b != 0, a != d and b>=0 == c >= 0
  var eps = LaPack0.dlamch( 'Epsilon' );
  aReference.setValue( 1. );
  dReference.setValue( 1. + 0.5 * eps );
  bReference.setValue( 0.5 * eps );
  cReference.setValue( 0.5 * eps );
  LaPack1.dlanv2(aReference,bReference,cReference,dReference,
    rt1rReference,rt1iReference,rt2rReference,rt2iReference,csReference,
    snReference);
  document.getElementById("debug_textarea").value +=
    "dlanv2: a,b,c,d = " + aReference.getValue() + " "
    + bReference.getValue() + " "
    + cReference.getValue() + " " + dReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "rt1r,rt1i,rt2r,rt2i = " + rt1rReference.getValue() + " "
    + rt1iReference.getValue() + " "
    + rt2rReference.getValue() + " "
    + rt2iReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "cs,sn = " + csReference.getValue() + " "
    + snReference.getValue() + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlaqge() {
  document.getElementById("debug_textarea").value +=
    "testing dlaqge *************" + "\n";
  var m = 4;
  var n = 3;
  var lda = 5;
  var ioffa = 1;
  var ioffc = 2;
  var ioffr = 3;
  var A = new Array( ioffa + lda * n );
  var c = new Array( ioffc + n );
  var r = new Array( ioffr + m );
  for ( var i = 0; i < m; i ++ ) r[ ioffr + i ] = 1 - i;
  for ( var j = 0; j < n; j ++ ) c[ ioffc + j ] = 1 + j;
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = i + j;
    }
  }
  var rowcnd = 2.;
  var amax = 5.;
  var colcnd = 1.;
  var equedReference = new StringReference();
  LaPack1.dlaqge(m,n,A,lda,r,c,rowcnd,colcnd,amax,equedReference,
    ioffa,ioffr,ioffc);
  document.getElementById("debug_textarea").value +=
    "dlaqge, equed = " + equedReference.getValue() + "\n";
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

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = i + j;
    }
  }
  rowcnd = 2.;
  amax = 5.;
  colcnd = 0.01;
  LaPack1.dlaqge(m,n,A,lda,r,c,rowcnd,colcnd,amax,equedReference,
    ioffa,ioffr,ioffc);
  document.getElementById("debug_textarea").value +=
    "dlaqge, equed = " + equedReference.getValue() + "\n";
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

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = i + j;
    }
  }
  rowcnd = 0.01;
  amax = 5.;
  colcnd = 1.;
  LaPack1.dlaqge(m,n,A,lda,r,c,rowcnd,colcnd,amax,equedReference,
    ioffa,ioffr,ioffc);
  document.getElementById("debug_textarea").value +=
    "dlaqge, equed = " + equedReference.getValue() + "\n";
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

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = i + j;
    }
  }
  rowcnd = 0.01;
  amax = 5.;
  colcnd = 0.01;
  LaPack1.dlaqge(m,n,A,lda,r,c,rowcnd,colcnd,amax,equedReference,
    ioffa,ioffr,ioffc);
  document.getElementById("debug_textarea").value +=
    "dlaqge, equed = " + equedReference.getValue() + "\n";
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
function test_dlaqsy() {
  document.getElementById("debug_textarea").value +=
    "testing dlaqsy *************" + "\n";
  var n = 3;
  var lda = 4;
  var ioffa = 1;
  var ioffs = 2;
  var A = new Array( ioffa + lda * n );
  var s = new Array( ioffs + n );
  for ( var i = 0; i < n; i ++ ) s[ ioffs + i ] = i + 1;
  for ( var j = 0; j < n; j ++ ) {
    for ( i=0; i < n; i ++ ) A[ ioffa + i + j * lda ] = i + 2 * j;
  }
  var scond = 1.;
  var amax = 1.;
  var equedReference = new StringReference();
  LaPack1.dlaqsy('U',n,A,lda,s,scond,amax,equedReference,ioffa,ioffs);
  document.getElementById("debug_textarea").value +=
    "dlaqsy('U',...): equed = " + equedReference.getValue() + "\n";
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
    + A[ ioffa + 2 + 2 * lda ]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i=0; i < n; i ++ ) A[ ioffa + i + j * lda ] = i + 2 * j;
  }
  scond = 1.;
  amax = 1.;
  LaPack1.dlaqsy('L',n,A,lda,s,scond,amax,equedReference,ioffa,ioffs);
  document.getElementById("debug_textarea").value +=
    "dlaqsy('L',...): equed = " + equedReference.getValue() + "\n";
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
    + A[ ioffa + 2 + 2 * lda ]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i=0; i < n; i ++ ) A[ ioffa + i + j * lda ] = i + 2 * j;
  }
  scond = 0.01;
  amax = 1.;
  LaPack1.dlaqsy('U',n,A,lda,s,scond,amax,equedReference,ioffa,ioffs);
  document.getElementById("debug_textarea").value +=
    "dlaqsy('U',...): equed = " + equedReference.getValue() + "\n";
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
    + A[ ioffa + 2 + 2 * lda ]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i=0; i < n; i ++ ) A[ ioffa + i + j * lda ] = i + 2 * j;
  }
  scond = 0.01;
  amax = 1.;
  LaPack1.dlaqsy('L',n,A,lda,s,scond,amax,equedReference,ioffa,ioffs);
  document.getElementById("debug_textarea").value +=
    "dlaqsy('L',...): equed = " + equedReference.getValue() + "\n";
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
    + A[ ioffa + 2 + 2 * lda ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlarf() {
  document.getElementById("debug_textarea").value +=
    "testing dlarf *************" + "\n";
  var m=3;
  var n=2;
  var ioffvl = 1;
  var ioffC = 2;
  var ioffworkl = 3;
  var vl=new Array(ioffvl + m);
  var ldc = 4;
  var C = new Array( ioffC + ldc * n );
  var workl = new Array( ioffworkl + n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( i + j );
    }
  }
  for ( i = 0; i < m; i ++ ) vl[ ioffvl + i ] = Number( 2 * i - 1 );
  LaPack1.dlarf('L',m,n,vl,1,2.,C,ldc,workl,ioffvl,ioffC,ioffworkl);
  document.getElementById("debug_textarea").value +=
    "dlarf('L',...): C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ]  + "\n";
  var ioffvr = 4;
  var ioffworkr = 5;
  var vr=new Array(ioffvr + n);
  var workr = new Array( ioffworkr + m );
  for ( i = 0; i < n; i ++ ) vr[ ioffvr + i ] = Number( 1 - 2 * i );
  LaPack1.dlarf('R',m,n,vr,1,2.,C,ldc,workr,ioffvr,ioffC,ioffworkr);
  document.getElementById("debug_textarea").value +=
    "dlarf('R',...): C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlarfb() {
  document.getElementById("debug_textarea").value +=
    "testing dlarfb *************" + "\n";

  var k = 3;
  var ldt = 5;
  var ldwork = 7;
  var ioffT = 1;
  var ioffWork = 2;
  var ioffV = 3;
  var ioffC = 4;
  var T = new Array( ioffT + ldt * k );
  var Work = new Array( ldwork * k );
  for ( var j=0; j<k; j++ ) {
    for ( var i=0; i<k; i++ ) T[ ioffT + i + j * ldt ] = i + j;
  }

//    storev='C', direct='F', side='L'
//    V=[V1] kxk unit lower
//      [V2] (m-k)xk
//    C=[C1] kxn
//      [C2] (m-k)xn
  var ldc = 6;
  var ldv = 4;
  var m = 4;
  var n = 5;
  var V = new Array( ioffV + ldv * k );
  for ( j=0; j<k; j++ ) {
    for ( i=j+1; i<m; i++ ) V[ ioffV + i + j * ldv ] = 2*i + 3*j;
  }
  var C = new Array( ioffC + ldc * n ); // m x n
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('L','N','F','C',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('L','N','F','C',...): C = "
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
    + C[ ioffC + 3 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 3 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ] + " "
    + C[ ioffC + 3 + 4 * ldc ]  + "\n";
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('L','T','F','C',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('L','T','F','C',...): C = "
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
    + C[ ioffC + 3 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 3 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ] + " "
    + C[ ioffC + 3 + 4 * ldc ]  + "\n";

  ldc=6;
  ldv=4;
  m = k;
  V = new Array( ioffV + ldv * k );
  for ( j=0; j<k; j++ ) {
    for ( i=j+1; i<m; i++ ) V[ ioffV + i + j * ldv ] = 2*i + 3*j;
  }
  C = new Array( ioffC + ldc * n );
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('L','N','F','C',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('L','N','F','C',...): C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ]  + "\n";
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('L','T','F','C',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('L','N','F','C',...): C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ]  + "\n";

//    storev='C', direct='F', side='R'
//    V=[V1] kxk unit lower
//      [V2] (n-k)xk
//    C=[C1,C2] mxk, mx(n-k)
  ldc=6;
  ldv=8;
  m = 4;
  n = 5;
  V = new Array( ioffV + ldv * k );
  for ( j=0; j<k; j++ ) {
    for ( i=j+1; i<n; i++ ) V[ ioffV + i + j * ldv ] = 2*i + 3*j;
  }
  C = new Array( ioffC + ldc * n );
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('R','N','F','C',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('R','N','F','C',...): C = "
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
    + C[ ioffC + 3 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 3 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ] + " "
    + C[ ioffC + 3 + 4 * ldc ]  + "\n";
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('R','T','F','C',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('R','T','F','C',...): C = "
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
    + C[ ioffC + 3 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 3 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ] + " "
    + C[ ioffC + 3 + 4 * ldc ]  + "\n";

  ldc=6;
  ldv=8;
  n = k;
  V = new Array( ioffV + ldv * k );
  for ( j=0; j<k; j++ ) {
    for ( i=j+1; i<n; i++ ) V[ ioffV + i + j * ldv ] = 2*i + 3*j;
  }
  C = new Array( ioffC + ldc * n );
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('R','N','F','C',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('R','N','F','C',...): C = "
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
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('R','T','F','C',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('R','T','F','C',...): C = "
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

//    storev='C', direct='B', side='L'
//    V=[V1] (m-k)xk
//      [V2] kxk unit upper
//    C=[C1] (m-k)xn
//      [C2] kxn
  ldc = 6;
  ldv = 8;
  m = 4;
  n = 5;
  V = new Array( ioffV + ldv * k );
  for ( j=0; j<k; j++ ) {
    for ( i=0; i<n+j-k; i++ ) V[ ioffV + i + j * ldv ] = 2*i + 3*j;
  }
  C = new Array( ioffC + ldc * n ); // m x n
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('L','N','B','C',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('L','N','B','C',...): C = "
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
    + C[ ioffC + 3 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 3 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ] + " "
    + C[ ioffC + 3 + 4 * ldc ]  + "\n";
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('L','T','B','C',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('L','N','B','C',...): C = "
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
    + C[ ioffC + 3 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 3 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ] + " "
    + C[ ioffC + 3 + 4 * ldc ]  + "\n";

  ldc = 6;
  ldv = 8;
  m = k;
  V = new Array( ioffV + ldv * k );
  for ( j=0; j<k; j++ ) {
    for ( i=0; i<n+j-k; i++ ) V[ ioffV + i + j * ldv ] = 2*i + 3*j;
  }
  C = new Array( ioffC + ldc * n );
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('L','N','B','C',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('L','N','B','C',...): C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ]  + "\n";
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('L','T','B','C',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('L','T','B','C',...): C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ]  + "\n";

//    storev='C', direct='B', side='R'
//    V=[V1] (n-k)xk
//      [V2] kxk unit upper
//    C=[C1,C2] mx(n-k),mxk
  ldc = 6;
  ldv = 8;
  m = 4;
  n = 5;
  V = new Array( ioffV + ldv * k );
  for ( j=0; j<k; j++ ) {
    for ( i=0; i<n+j-k; i++ ) V[ ioffV + i + j * ldv ] = 2*i + 3*j;
  }
  C = new Array( ioffC + ldc * n );
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('R','N','B','C',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('R','N','B','C',...): C = "
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
    + C[ ioffC + 3 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 3 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ] + " "
    + C[ ioffC + 3 + 4 * ldc ]  + "\n";
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('R','T','B','C',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('R','T','B','C',...): C = "
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
    + C[ ioffC + 3 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 3 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ] + " "
    + C[ ioffC + 3 + 4 * ldc ]  + "\n";

  ldc = 6;
  ldv = 8;
  n = k;
  V = new Array( ioffV + ldv * k );
  for ( j=0; j<k; j++ ) {
    for ( i=0; i<n+j-k; i++ ) V[ ioffV + i + j * ldv ] = 2*i + 3*j;
  }
  C = new Array( ioffC + ldc * n );
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('R','N','B','C',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('R','N','B','C',...): C = "
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
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('R','T','B','C',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('R','N','B','C',...): C = "
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

//    storev='R', direct='F', side='L'
//    V=[V1,V2] kxk unit upper, kx(m-k)
//    C=[C1] kxn
//      [C2] (m-k)xn
  V = new Array( ldv * m );
  ldc = 6;
  ldv = 4;
  m = 4;
  n = 5;
  for ( j=0; j<m; j++ ) {
    for ( i=0; i<Math.min(j,k); i++ ) {
      V[ ioffV + i + j * ldv ] = 2*i + 3*j;
    }
  }
  C = new Array( ioffC + ldc * n ); // m x n
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('L','N','F','R',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('L','N','F','R',...): C = "
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
    + C[ ioffC + 3 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 3 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ] + " "
    + C[ ioffC + 3 + 4 * ldc ]  + "\n";
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('L','T','F','R',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('L','N','F','R',...): C = "
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
    + C[ ioffC + 3 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 3 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ] + " "
    + C[ ioffC + 3 + 4 * ldc ]  + "\n";

  ldc = 6;
  ldv = 4;
  m = k;
  V = new Array( ioffV + ldv * m );
  for ( j=0; j<m; j++ ) {
    for ( i=0; i<Math.min(j,k); i++ ) {
      V[ ioffV + i + j * ldv ] = 2*i + 3*j;
    }
  }
  C = new Array( ioffC + ldc * n );
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('L','N','F','R',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('L','N','F','R',...): C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ]  + "\n";
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('L','T','F','R',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('L','T','F','R',...): C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ]  + "\n";

//    storev='R', direct='F', side='R'
//    V=[V1,V2] kxk unit upper, kx(n-k)
//    C=[C1,C2] mxk, mx(n-k)
  ldc = 6;
  ldv = 4;
  V = new Array( ioffV + ldv * n );
  m = 4;
  n = 5;
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<Math.min(j,k); i++ ) {
      V[ ioffV + i + j * ldv ] = 2*i + 3*j;
    }
  }
  C = new Array( ioffC + ldc * n );
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('R','N','F','R',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('R','N','F','R',...): C = "
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
    + C[ ioffC + 3 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 3 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ] + " "
    + C[ ioffC + 3 + 4 * ldc ]  + "\n";
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('R','T','F','R',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('R','T','F','R',...): C = "
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
    + C[ ioffC + 3 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 3 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ] + " "
    + C[ ioffC + 3 + 4 * ldc ]  + "\n";

  ldc = 6;
  ldv = 8;
  n = k;
  V = new Array( ioffV + ldv * n );
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<Math.min(j,k); i++ ) {
      V[ ioffV + i + j * ldv ] = 2*i + 3*j;
    }
  }
  C = new Array( ioffC + ldc * n );
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('R','N','F','R',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('R','N','F','R',...): C = "
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
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('R','T','F','R',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('R','N','F','R',...): C = "
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

//    storev='R', direct='B', side='L'
//    V=[V1,V2] kx(m-k), kxk unit lower
//    C=[C1] (m-k)xn
//      [C2] kxn
  V = new Array( ioffV + ldv * k );
  ldc=6;
  ldv=4;
  m = 4;
  n = 5;
  for ( j=0; j<m; j++ ) {
    for ( i=Math.max(0,k+j-m+1); i<k; i++ ) {
      V[ ioffV + i + j * ldv ] = 2*i + 3*j;
    }
  }
  C = new Array( ioffC + ldc * n ); // m x n
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('L','N','B','R',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('L','N','B','R',...): C = "
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
    + C[ ioffC + 3 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 3 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ] + " "
    + C[ ioffC + 3 + 4 * ldc ]  + "\n";
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('L','T','B','R',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('L','N','B','R',...): C = "
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
    + C[ ioffC + 3 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 3 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ] + " "
    + C[ ioffC + 3 + 4 * ldc ]  + "\n";

  ldc=6;
  ldv=4;
  m = k;
  V = new Array( ioffV + ldv * k );
  for ( j=0; j<k; j++ ) {
    for ( i=j+1; i<m; i++ ) V[ ioffV + i + j * ldv ] = 2*i + 3*j;
  }
  C = new Array( ioffC + ldc * n );
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('L','N','B','R',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('L','N','B','R',...): C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ]  + "\n";
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('L','T','B','R',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('L','T','B','R',...): C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ]  + "\n";

//    storev='R', direct='B', side='R'
//    V=[V1,V2] kx(n-k), kxk unit lower
//    C=[C1,C2] mx(n-k), mxk
  V = new Array( ioffV + ldv * k );
  ldc=6;
  ldv=4;
  m = 4;
  n = 5;
  for ( j=0; j<n; j++ ) {
    for ( i=Math.max(0,k+j-n+1); i<k; i++ ) {
      V[ ioffV + i + j * ldv ] = 2*i + 3*j;
    }
  }
  C = new Array( ioffC + ldc * n );
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('R','N','B','R',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('R','N','B','R',...): C = "
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
    + C[ ioffC + 3 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 3 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ] + " "
    + C[ ioffC + 3 + 4 * ldc ]  + "\n";
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('R','T','B','R',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('R','T','B','R',...): C = "
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
    + C[ ioffC + 3 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + " "
    + C[ ioffC + 2 + 3 * ldc ] + " "
    + C[ ioffC + 3 + 3 * ldc ] + " "
    + C[ ioffC + 0 + 4 * ldc ] + " "
    + C[ ioffC + 1 + 4 * ldc ] + " "
    + C[ ioffC + 2 + 4 * ldc ] + " "
    + C[ ioffC + 3 + 4 * ldc ]  + "\n";

  ldc=6;
  ldv=4;
  n = k;
  V = new Array( ioffV + ldv * k );
  for ( j=0; j<k; j++ ) {
    for ( i=j+1; i<n; i++ ) V[ ioffV + i + j * ldv ] = 2*i + 3*j;
  }
  C = new Array( ioffC + ldc * n );
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('R','N','B','R',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('R','N','B','R',...): C = "
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
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<m; i++ ) C[ ioffC + i + j * ldc ] = 3*i - 2*j;
  }
  LaPack1.dlarfb('R','T','B','R',m,n,k,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffV,ioffT,ioffC,ioffWork);
  document.getElementById("debug_textarea").value +=
    "dlarfb('R','N','B','R',...): C = "
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
function test_dlarfg() {
  document.getElementById("debug_textarea").value +=
    "testing dlarfg *************" + "\n";
  var n = 3;
  var incx = 1;
  var incy = 2;
  var safmin = LaPack0.dlamch('S')/LaPack0.dlamch('E');
  var ioffx = 1;
  var ioffy = 2;
  var x = new Array( ioffx + n * incx );
  var y = new Array( ioffy + n * incy );
  for ( var i = 0; i < n * incx; i += incx ) {
    x[ ioffx + i ] = safmin / Number( 3 * i + 1 ); 
  }
  var alphaReference = new NumberReference( 0.25 * safmin );
  var tauReference = new NumberReference();
  LaPack1.dlarfg(n,alphaReference,x,incx,tauReference,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlarfg: alpha,tau = " + alphaReference.getValue() + " "
    + tauReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";

  for ( i = 0; i < n * incy; i += incy ) {
    y[ ioffy + i ] = safmin / Number( 8 * i + 1 );
  }
  alphaReference.setValue( 0.25 * safmin );
  LaPack1.dlarfg(n,alphaReference,y,incy,tauReference,ioffy);
  document.getElementById("debug_textarea").value +=
    "dlarfg: alpha,tau = " + alphaReference.getValue() + " "
    + tauReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "y = "
    + y[ ioffy + 0 ] + " "
    + y[ ioffy + 1 ] + " "
    + y[ ioffy + 2 ] + " "
    + y[ ioffy + 3 ] + " "
    + y[ ioffy + 4 ] + " "
    + y[ ioffy + 5 ]  + "\n";

  for ( i = 0; i < n * incx; i += incx ) x[ ioffx + i ] = i + 1;
  alphaReference.setValue( 3. );
  LaPack1.dlarfg(n,alphaReference,x,incx,tauReference,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlarfg: alpha,tau = " + alphaReference.getValue() + " "
    + tauReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";

  for ( i = 0; i < n * incy; i += incy ) y[ ioffy + i ] = i + 1;
  alphaReference.setValue( 3. );
  LaPack1.dlarfg(n,alphaReference,y,incy,tauReference,ioffy);
  document.getElementById("debug_textarea").value +=
    "dlarfg: alpha,tau = " + alphaReference.getValue() + " "
    + tauReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "y = "
    + y[ ioffy + 0 ] + " "
    + y[ ioffy + 1 ] + " "
    + y[ ioffy + 2 ] + " "
    + y[ ioffy + 3 ] + " "
    + y[ ioffy + 4 ] + " "
    + y[ ioffy + 5 ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlarfgp() {
  document.getElementById("debug_textarea").value +=
    "testing dlarfgp *************" + "\n";
  var n = 3;
  var incx = 1;
  var incy = 2;
  var safmin = LaPack0.dlamch('S')/LaPack0.dlamch('E');
  var ioffx = 1;
  var ioffy = 2;
  var x = new Array( ioffx + n * incx );
  var y = new Array( ioffy + n * incy );
  for ( var i = 0; i < n * incx; i += incx ) {
    x[ ioffx + i ] = safmin / Number( 3 * i + 1 ); 
  }
  var alphaReference = new NumberReference( 0.25 * safmin );
  var tauReference = new NumberReference();
  LaPack1.dlarfgp(n,alphaReference,x,incx,tauReference,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlarfgp: alpha,tau = " + alphaReference.getValue() + " "
    + tauReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";

  for ( i = 0; i < n * incy; i += incy ) {
    y[ ioffy + i ] = safmin / Number( 8 * i + 1 );
  }
  alphaReference.setValue( 0.25 * safmin );
  LaPack1.dlarfgp(n,alphaReference,y,incy,tauReference,ioffy);
  document.getElementById("debug_textarea").value +=
    "dlarfgp: alpha,tau = " + alphaReference.getValue() + " "
    + tauReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "y = "
    + y[ ioffy + 0 ] + " "
    + y[ ioffy + 1 ] + " "
    + y[ ioffy + 2 ] + " "
    + y[ ioffy + 3 ] + " "
    + y[ ioffy + 4 ] + " "
    + y[ ioffy + 5 ]  + "\n";

  for ( i = 0; i < n * incx; i += incx ) x[ ioffx + i ] = i + 1;
  alphaReference.setValue( 3. );
  LaPack1.dlarfgp(n,alphaReference,x,incx,tauReference,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlarfgp: alpha,tau = " + alphaReference.getValue() + " "
    + tauReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";

  for ( i = 0; i < n * incy; i += incy ) y[ ioffy + i ] = i + 1;
  alphaReference.setValue( 3. );
  LaPack1.dlarfgp(n,alphaReference,y,incy,tauReference,ioffy);
  document.getElementById("debug_textarea").value +=
    "dlarfgp: alpha,tau = " + alphaReference.getValue() + " "
    + tauReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "y = "
    + y[ ioffy + 0 ] + " "
    + y[ ioffy + 1 ] + " "
    + y[ ioffy + 2 ] + " "
    + y[ ioffy + 3 ] + " "
    + y[ ioffy + 4 ] + " "
    + y[ ioffy + 5 ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlarrk() {
  document.getElementById("debug_textarea").value +=
    "testing dlarrk *************" + "\n";
  var n = 64;
  var iw = 29;
  var gl = 1.8;
  var gu = 2.2;
  var pivmin = 1.e-4;
  var reltol = 1.e-14;
  var ioffd = 1;
  var ioffe2 = 2;
  var d = new Array( ioffd + n );
  var e2 = new Array( ioffe2 + n - 1 );
  for ( var i = 0; i < n; i ++ ) d[ ioffd + i ] = 2.;
  for ( i = 0; i < n - 1; i ++ ) e2[ ioffe2 + i ] = 1.;
  var wReference = new NumberReference();
  var werrReference = new NumberReference();
  var info = new IntReference();
  LaPack1.dlarrk(n,iw,gl,gu,d,e2,pivmin,reltol,wReference,werrReference,
    info,ioffd,ioffe2);
  document.getElementById("debug_textarea").value +=
    "dlarrk: w,werr,info = " + wReference.getValue() + " "
    + werrReference.getValue() + " "
    + info.getValue()  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlartg() {
  document.getElementById("debug_textarea").value +=
    "testing dlartg *************" + "\n";
  var f = 2.;
  var g = 0.;
  var csReference = new NumberReference();
  var snReference = new NumberReference();
  var rReference = new NumberReference();
  LaPack1.dlartg(f,g,csReference,snReference,rReference);
  document.getElementById("debug_textarea").value +=
    "dlartg: cs,sn,r = " + csReference.getValue() + " "
    + snReference.getValue() + " " + rReference.getValue() + "\n";

  f=0.;
  g=3.;
  LaPack1.dlartg(f,g,csReference,snReference,rReference);
  document.getElementById("debug_textarea").value +=
    "dlartg: cs,sn,r = " + csReference.getValue() + " "
    + snReference.getValue() + " " + rReference.getValue() + "\n";

  var safmin = LaPack0.dlamch( 'S' );
  var eps = LaPack0.dlamch( 'E' );
  var safmn2 = Math.pow( LaPack0.dlamch( 'B' ),
    Math.round( Math.log( safmin / eps )
       / Math.log( LaPack0.dlamch( 'B' ) ) / 2. ) );
  var safmx2 = 1. / safmn2;
  f = safmx2 * 2.;
  g = - safmx2 * 4.;
  LaPack1.dlartg(f,g,csReference,snReference,rReference);
  document.getElementById("debug_textarea").value +=
    "dlartg: cs,sn,r = " + csReference.getValue() + " "
    + snReference.getValue() + " " + rReference.getValue() + "\n";

  f = safmn2 * 0.5;
  g = - safmn2 * 0.25;
  LaPack1.dlartg(f,g,csReference,snReference,rReference);
  document.getElementById("debug_textarea").value +=
    "dlartg: cs,sn,r = " + csReference.getValue() + " "
    + snReference.getValue() + " " + rReference.getValue() + "\n";

  f = 2.;
  g = -4.;
  LaPack1.dlartg(f,g,csReference,snReference,rReference);
  document.getElementById("debug_textarea").value +=
    "dlartg: cs,sn,r = " + csReference.getValue() + " "
    + snReference.getValue() + " " + rReference.getValue() + "\n";

  f = -2.;
  g = 4.;
  LaPack1.dlartg(f,g,csReference,snReference,rReference);
  document.getElementById("debug_textarea").value +=
    "dlartg: cs,sn,r = " + csReference.getValue() + " "
    + snReference.getValue() + " " + rReference.getValue() + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlartgp() {
  document.getElementById("debug_textarea").value +=
    "testing dlartgp *************" + "\n";
  var f = 2.;
  var g = 0.;
  var csReference = new NumberReference();
  var snReference = new NumberReference();
  var rReference = new NumberReference();
  LaPack1.dlartgp(f,g,csReference,snReference,rReference);
  document.getElementById("debug_textarea").value +=
    "dlartgp: cs,sn,r = " + csReference.getValue() + " "
    + snReference.getValue() + " " + rReference.getValue() + "\n";

  f=0.;
  g=3.;
  LaPack1.dlartgp(f,g,csReference,snReference,rReference);
  document.getElementById("debug_textarea").value +=
    "dlartgp: cs,sn,r = " + csReference.getValue() + " "
    + snReference.getValue() + " " + rReference.getValue() + "\n";

  var safmin = LaPack0.dlamch( 'S' );
  var eps = LaPack0.dlamch( 'E' );
  var safmn2 = Math.pow( LaPack0.dlamch( 'B' ),
    Math.round( Math.log( safmin / eps )
       / Math.log( LaPack0.dlamch( 'B' ) ) / 2. ) );
  var safmx2 = 1. / safmn2;
  f = safmx2 * 2.;
  g = - safmx2 * 4.;
  LaPack1.dlartgp(f,g,csReference,snReference,rReference);
  document.getElementById("debug_textarea").value +=
    "dlartgp: cs,sn,r = " + csReference.getValue() + " "
    + snReference.getValue() + " " + rReference.getValue() + "\n";

  f = safmn2 * 0.5;
  g = - safmn2 * 0.25;
  LaPack1.dlartgp(f,g,csReference,snReference,rReference);
  document.getElementById("debug_textarea").value +=
    "dlartgp: cs,sn,r = " + csReference.getValue() + " "
    + snReference.getValue() + " " + rReference.getValue() + "\n";

  f = 2.;
  g = -4.;
  LaPack1.dlartgp(f,g,csReference,snReference,rReference);
  document.getElementById("debug_textarea").value +=
    "dlartgp: cs,sn,r = " + csReference.getValue() + " "
    + snReference.getValue() + " " + rReference.getValue() + "\n";

  f = -2.;
  g = 4.;
  LaPack1.dlartgp(f,g,csReference,snReference,rReference);
  document.getElementById("debug_textarea").value +=
    "dlartgp: cs,sn,r = " + csReference.getValue() + " "
    + snReference.getValue() + " " + rReference.getValue() + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlascl() {
  document.getElementById("debug_textarea").value +=
    "testing dlascl *************" + "\n";
  var m = 4;
  var n = 3;
  var lda = 7;
  var kl = 2;
  var ku = 2;
  var cto = 4.;
  var cfrom = 2.;
  var ioffa = 1;
  var A = new Array( ioffa + lda * n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = i + j + 1;
    }
  }
  var info = new IntReference();
  LaPack1.dlascl('G',kl,ku,cfrom,cto,m,n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dlascl('G',...): A = "
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

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = i + j + 1;
    }
  }
  LaPack1.dlascl('L',kl,ku,cfrom,cto,m,n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dlascl('L',...): A = "
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

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = i + j + 1;
    }
  }
  LaPack1.dlascl('U',kl,ku,cfrom,cto,m,n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dlascl('U',...): A = "
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

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = i + j + 1;
    }
  }
  LaPack1.dlascl('H',kl,ku,cfrom,cto,n,n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dlascl('H',...): A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = i + j + 1;
    }
  }
  LaPack1.dlascl('B',kl,ku,cfrom,cto,n,n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dlascl('B',...): A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = i + j + 1;
    }
  }
  LaPack1.dlascl('Q',kl,ku,cfrom,cto,n,n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dlascl('Q',...): A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = i + j + 1;
    }
  }
  LaPack1.dlascl('Z',kl,ku,cfrom,cto,m,n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dlascl('Z',...): A = "
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
function test_dlasv2() {
  document.getElementById("debug_textarea").value +=
    "testing dlasv2 *************" + "\n";
  var ssminReference = new NumberReference();
  var ssmaxReference = new NumberReference();
  var snrReference = new NumberReference();
  var csrReference = new NumberReference();
  var snlReference = new NumberReference();
  var cslReference = new NumberReference();

  LaPack1.dlasv2(1.,0.,2.,ssminReference,ssmaxReference,snrReference,
    csrReference,snlReference,cslReference);
  document.getElementById("debug_textarea").value +=
    "dlasv2(1.,0.,2.,...): ssmin,ssmax = "
    + ssminReference.getValue() + " "
    + ssmaxReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "snr,csr,snl,csl = " + snrReference.getValue() + " "
    + csrReference.getValue() + " "
    + snlReference.getValue() + " " + cslReference.getValue() + "\n";

  LaPack1.dlasv2(2.,0.,1.,ssminReference,ssmaxReference,snrReference,
    csrReference,snlReference,cslReference);
  document.getElementById("debug_textarea").value +=
    "dlasv2(2.,0.,1.,...): ssmin,ssmax = "
    + ssminReference.getValue() + " "
    + ssmaxReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "snr,csr,snl,csl = " + snrReference.getValue() + " "
    + csrReference.getValue() + " "
    + snlReference.getValue() + " " + cslReference.getValue() + "\n";

  LaPack1.dlasv2(1.,2.,3.,ssminReference,ssmaxReference,snrReference,
    csrReference,snlReference,cslReference);
  document.getElementById("debug_textarea").value +=
    "dlasv2(1.,2.,3.,...): ssmin,ssmax = "
    + ssminReference.getValue() + " "
    + ssmaxReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "snr,csr,snl,csl = " + snrReference.getValue() + " "
    + csrReference.getValue() + " "
    + snlReference.getValue() + " " + cslReference.getValue() + "\n";

  LaPack1.dlasv2(1.,4.,3.,ssminReference,ssmaxReference,snrReference,
    csrReference,snlReference,cslReference);
  document.getElementById("debug_textarea").value +=
    "dlasv2(1.,4.,3.,...): ssmin,ssmax = "
    + ssminReference.getValue() + " "
    + ssmaxReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "snr,csr,snl,csl = " + snrReference.getValue() + " "
    + csrReference.getValue() + " "
    + snlReference.getValue() + " " + cslReference.getValue() + "\n";

  LaPack1.dlasv2(2.,1.,3.,ssminReference,ssmaxReference,snrReference,
    csrReference,snlReference,cslReference);
  document.getElementById("debug_textarea").value +=
    "dlasv2(2.,1.,3.,...): ssmin,ssmax = "
    + ssminReference.getValue() + " "
    + ssmaxReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "snr,csr,snl,csl = " + snrReference.getValue() + " "
    + csrReference.getValue() + " "
    + snlReference.getValue() + " " + cslReference.getValue() + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlasy2() {
  document.getElementById("debug_textarea").value +=
    "testing dlasy2 *************" + "\n";
  var ldtl = 3;
  var ldtr = 4;
  var ldb = 5;
  var ldx = 6;
  var ioffb = 1;
  var iofftl = 2;
  var iofftr = 3;
  var ioffx = 4;
  var B = new Array( ioffb + ldb * 2 );
  var Tl = new Array( iofftl + ldtl * 2 );
  var Tr = new Array( iofftr + ldtr * 2 );
  var X = new Array( ioffx + ldx * 2 );
  var scaleReference = new NumberReference();
  var xnormReference = new NumberReference();
  var info = new IntReference();
  for ( var j = 0; j < 2; j ++ ) {
    for ( var i = 0; i < 2; i ++ ) {
      B[ ioffb + i + j * ldb ] = 1 + i + 2 * j;
    }
  }
  Tl[ iofftl ] = 3.;
  Tr[ iofftr ] = 2.;
  LaPack1.dlasy2(true,true,1,1,1,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(true,true,1,1,1,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = " + X[ ioffx ] + "\n";
  LaPack1.dlasy2(true,true,-1,1,1,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(true,true,-1,1,1,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = " + X[ ioffx ] + "\n";

  Tr[ iofftr + 0 + 0 * ldtr ] = 2.;
  Tr[ iofftr + 1 + 0 * ldtr ] = 1.;
  Tr[ iofftr + 0 + 1 * ldtr ] = -1.;
  Tr[ iofftr + 1 + 1 * ldtr ] = 2.;
  LaPack1.dlasy2(true,true,1,1,2,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(true,true,1,1,2,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = " + X[ ioffx + 0 * ldx ] + " " + X[ ioffx + 1 * ldx] + "\n";
  LaPack1.dlasy2(true,true,-1,1,2,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(true,true,-1,1,2,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = " + X[ ioffx + 0 * ldx ] + " " + X[ ioffx + 1 * ldx] + "\n";
  LaPack1.dlasy2(false,true,1,1,2,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(false,true,1,1,2,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = " + X[ ioffx + 0 * ldx ] + " " + X[ ioffx + 1 * ldx] + "\n";
  LaPack1.dlasy2(false,true,-1,1,2,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(false,true,-1,1,2,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = " + X[ ioffx + 0 * ldx ] + " " + X[ ioffx + 1 * ldx] + "\n";
  LaPack1.dlasy2(true,false,1,1,2,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(true,false,1,1,2,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = " + X[ ioffx + 0 * ldx ] + " " + X[ ioffx + 1 * ldx] + "\n";
  LaPack1.dlasy2(true,false,-1,1,2,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(true,false,-1,1,2,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = " + X[ ioffx + 0 * ldx ] + " " + X[ ioffx + 1 * ldx] + "\n";
  LaPack1.dlasy2(false,false,1,1,2,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(false,false,1,1,2,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = " + X[ ioffx + 0 * ldx ] + " " + X[ ioffx + 1 * ldx] + "\n";
  LaPack1.dlasy2(false,false,-1,1,2,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(false,false,-1,1,2,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = " + X[ ioffx + 0 * ldx ] + " " + X[ ioffx + 1 * ldx] + "\n";

  Tl[ iofftl + 0 + 0 * ldtl ] = 2.;
  Tl[ iofftl + 1 + 0 * ldtl ] = 1.;
  Tl[ iofftl + 0 + 1 * ldtl ] = -1.;
  Tl[ iofftl + 1 + 1 * ldtl ] = 2.;
  Tr[ iofftr + 0 + 0 * ldtr ] = 3.;
  LaPack1.dlasy2(true,true,1,2,1,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(true,true,1,2,1,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = " + X[ ioffx + 0 ] + " " + X[ ioffx + 1 ] + "\n";
  LaPack1.dlasy2(true,true,-1,2,1,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(true,true,-1,2,1,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = " + X[ ioffx + 0 ] + " " + X[ ioffx + 1 ] + "\n";
  LaPack1.dlasy2(false,true,1,2,1,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(false,true,1,2,1,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = " + X[ ioffx + 0 ] + " " + X[ ioffx + 1 ] + "\n";
  LaPack1.dlasy2(false,true,-1,2,1,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(false,true,-1,2,1,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = " + X[ ioffx + 0 ] + " " + X[ ioffx + 1 ] + "\n";
  LaPack1.dlasy2(true,false,1,2,1,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(true,false,1,2,1,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = " + X[ ioffx + 0 ] + " " + X[ ioffx + 1 ] + "\n";
  LaPack1.dlasy2(true,false,-1,2,1,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(true,false,-1,2,1,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = " + X[ ioffx + 0 ] + " " + X[ ioffx + 1 ] + "\n";
  LaPack1.dlasy2(false,false,1,2,1,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(false,false,1,2,1,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = " + X[ ioffx + 0 ] + " " + X[ ioffx + 1 ] + "\n";
  LaPack1.dlasy2(false,false,-1,2,1,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(false,false,-1,2,1,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = " + X[ ioffx + 0 ] + " " + X[ ioffx + 1 ] + "\n";

  Tl[ iofftl + 0 + 0 * ldtl ] = 2.;
  Tl[ iofftl + 1 + 0 * ldtl ] = 1.;
  Tl[ iofftl + 0 + 1 * ldtl ] = -1.;
  Tl[ iofftl + 1 + 1 * ldtl ] = 2.;
  Tr[ iofftr + 0 + 0 * ldtr ] = 3.;
  Tr[ iofftr + 1 + 0 * ldtr ] = -1.;
  Tr[ iofftr + 0 + 1 * ldtr ] = 1.;
  Tr[ iofftr + 1 + 1 * ldtr ] = 3.;
  LaPack1.dlasy2(true,true,1,2,2,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(true,true,1,2,2,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = "
    + X[ ioffx + 0 + 0 * ldx ] + " "
    + X[ ioffx + 1 + 0 * ldx ] + " "
    + X[ ioffx + 0 + 1 * ldx ] + " "
    + X[ ioffx + 1 + 1 * ldx ]  + "\n";
  LaPack1.dlasy2(true,true,-1,2,2,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(true,true,-1,2,2,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = "
    + X[ ioffx + 0 + 0 * ldx ] + " "
    + X[ ioffx + 1 + 0 * ldx ] + " "
    + X[ ioffx + 0 + 1 * ldx ] + " "
    + X[ ioffx + 1 + 1 * ldx ]  + "\n";
  LaPack1.dlasy2(false,true,1,2,2,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(false,true,1,2,2,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = "
    + X[ ioffx + 0 + 0 * ldx ] + " "
    + X[ ioffx + 1 + 0 * ldx ] + " "
    + X[ ioffx + 0 + 1 * ldx ] + " "
    + X[ ioffx + 1 + 1 * ldx ]  + "\n";
  LaPack1.dlasy2(false,true,-1,2,2,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(false,true,-1,2,2,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = "
    + X[ ioffx + 0 + 0 * ldx ] + " "
    + X[ ioffx + 1 + 0 * ldx ] + " "
    + X[ ioffx + 0 + 1 * ldx ] + " "
    + X[ ioffx + 1 + 1 * ldx ]  + "\n";
  LaPack1.dlasy2(true,false,1,2,2,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(true,false,1,2,2,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = "
    + X[ ioffx + 0 + 0 * ldx ] + " "
    + X[ ioffx + 1 + 0 * ldx ] + " "
    + X[ ioffx + 0 + 1 * ldx ] + " "
    + X[ ioffx + 1 + 1 * ldx ]  + "\n";
  LaPack1.dlasy2(true,false,-1,2,2,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(true,false,-1,2,2,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = "
    + X[ ioffx + 0 + 0 * ldx ] + " "
    + X[ ioffx + 1 + 0 * ldx ] + " "
    + X[ ioffx + 0 + 1 * ldx ] + " "
    + X[ ioffx + 1 + 1 * ldx ]  + "\n";
  LaPack1.dlasy2(false,false,1,2,2,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(false,false,1,2,2,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = "
    + X[ ioffx + 0 + 0 * ldx ] + " "
    + X[ ioffx + 1 + 0 * ldx ] + " "
    + X[ ioffx + 0 + 1 * ldx ] + " "
    + X[ ioffx + 1 + 1 * ldx ]  + "\n";
  LaPack1.dlasy2(false,false,-1,2,2,Tl,ldtl,Tr,ldtr,B,ldb,scaleReference,
    X,ldx,xnormReference,info,iofftl,iofftr,ioffb,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlasy2(false,false,-1,2,2,...): info,scale,xnorm = "
    + info.getValue() + " " + scaleReference.getValue() + " "
    + xnormReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "X = "
    + X[ ioffx + 0 + 0 * ldx ] + " "
    + X[ ioffx + 1 + 0 * ldx ] + " "
    + X[ ioffx + 0 + 1 * ldx ] + " "
    + X[ ioffx + 1 + 1 * ldx ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlatrs() {
  document.getElementById("debug_textarea").value +=
    "testing dlatrs *************" + "\n";
  var n = 3;
  var lda = 4;
  var ioffa = 1;
  var ioffcnorm = 2;
  var ioffx = 3;
  var A = new Array( ioffa + lda * n );
  var cnorm = new Array( ioffcnorm + n );
  var x = new Array( ioffx + n );
  A[ ioffa + 0 + 0 * lda ] = 2.;
  A[ ioffa + 1 + 0 * lda ] = 1.;
  A[ ioffa + 2 + 0 * lda ] = 0.;
  A[ ioffa + 0 + 1 * lda ] = 3.;
  A[ ioffa + 1 + 1 * lda ] = 3.;
  A[ ioffa + 2 + 1 * lda ] = 2.;
  A[ ioffa + 0 + 2 * lda ] = 4.;
  A[ ioffa + 1 + 2 * lda ] = 4.;
  A[ ioffa + 2 + 2 * lda ] = 4.;
  for ( var i = 0; i < n; i ++ ) x[ ioffx + i ] = i - 1;
  var scaleReference = new NumberReference();
  var info = new IntReference();
  LaPack1.dlatrs('U','N','U','N',n,A,lda,x,scaleReference,cnorm,info,
    ioffa,ioffx,ioffcnorm);
  document.getElementById("debug_textarea").value +=
    "dlatrs('U','N','U','N',...): info,scale = " + info.getValue()
    + " " + scaleReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "cnorm = "
    + cnorm[ ioffcnorm + 0 ] + " "
    + cnorm[ ioffcnorm + 1 ] + " "
    + cnorm[ ioffcnorm + 2 ]  + "\n";
  for ( i = 0; i < n; i ++ ) {
    x[ ioffx + i ] = i - 1;
    cnorm[ ioffcnorm + i ] = 16.;
  }
  LaPack1.dlatrs('U','N','U','Y',n,A,lda,x,scaleReference,cnorm,info,
    ioffa,ioffx,ioffcnorm);
  document.getElementById("debug_textarea").value +=
    "dlatrs('U','N','U','Y',...): info,scale = " + info.getValue()
    + " " + scaleReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "cnorm = "
    + cnorm[ ioffcnorm + 0 ] + " "
    + cnorm[ ioffcnorm + 1 ] + " "
    + cnorm[ ioffcnorm + 2 ]  + "\n";

  for ( i = 0; i < n; i ++ ) x[ ioffx + i ] = i - 1;
  LaPack1.dlatrs('U','N','N','N',n,A,lda,x,scaleReference,cnorm,info,
    ioffa,ioffx,ioffcnorm);
  document.getElementById("debug_textarea").value +=
    "dlatrs('U','N','N','N',...): info,scale = " + info.getValue()
    + " " + scaleReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "cnorm = "
    + cnorm[ ioffcnorm + 0 ] + " "
    + cnorm[ ioffcnorm + 1 ] + " "
    + cnorm[ ioffcnorm + 2 ]  + "\n";
  for ( i = 0; i < n; i ++ ) {
    x[ ioffx + i ] = i - 1;
    cnorm[ ioffcnorm + i ] = 16.;
  }
  LaPack1.dlatrs('U','N','N','Y',n,A,lda,x,scaleReference,cnorm,info,
    ioffa,ioffx,ioffcnorm);
  document.getElementById("debug_textarea").value +=
    "dlatrs('U','N','N','Y',...): info,scale = " + info.getValue()
    + " " + scaleReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "cnorm = "
    + cnorm[ ioffcnorm + 0 ] + " "
    + cnorm[ ioffcnorm + 1 ] + " "
    + cnorm[ ioffcnorm + 2 ]  + "\n";

  for ( i = 0; i < n; i ++ ) x[ ioffx + i ] = i - 1;
  LaPack1.dlatrs('U','T','U','N',n,A,lda,x,scaleReference,cnorm,info,
    ioffa,ioffx,ioffcnorm);
  document.getElementById("debug_textarea").value +=
    "dlatrs('U','T','U','N',...): info,scale = " + info.getValue()
    + " " + scaleReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "cnorm = "
    + cnorm[ ioffcnorm + 0 ] + " "
    + cnorm[ ioffcnorm + 1 ] + " "
    + cnorm[ ioffcnorm + 2 ]  + "\n";
  for ( i = 0; i < n; i ++ ) {
    x[ ioffx + i ] = i - 1;
    cnorm[ ioffcnorm + i ] = 16.;
  }
  LaPack1.dlatrs('U','T','U','Y',n,A,lda,x,scaleReference,cnorm,info,
    ioffa,ioffx,ioffcnorm);
  document.getElementById("debug_textarea").value +=
    "dlatrs('U','T','U','Y',...): info,scale = " + info.getValue()
    + " " + scaleReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "cnorm = "
    + cnorm[ ioffcnorm + 0 ] + " "
    + cnorm[ ioffcnorm + 1 ] + " "
    + cnorm[ ioffcnorm + 2 ]  + "\n";

  for ( i = 0; i < n; i ++ ) x[ ioffx + i ] = i - 1;
  LaPack1.dlatrs('U','T','N','N',n,A,lda,x,scaleReference,cnorm,info,
    ioffa,ioffx,ioffcnorm);
  document.getElementById("debug_textarea").value +=
    "dlatrs('U','T','N','N',...): info,scale = " + info.getValue()
    + " " + scaleReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "cnorm = "
    + cnorm[ ioffcnorm + 0 ] + " "
    + cnorm[ ioffcnorm + 1 ] + " "
    + cnorm[ ioffcnorm + 2 ]  + "\n";
  for ( i = 0; i < n; i ++ ) {
    x[ ioffx + i ] = i - 1;
    cnorm[ ioffcnorm + i ] = 16.;
  }
  LaPack1.dlatrs('U','T','N','Y',n,A,lda,x,scaleReference,cnorm,info,
    ioffa,ioffx,ioffcnorm);
  document.getElementById("debug_textarea").value +=
    "dlatrs('U','T','N','Y',...): info,scale = " + info.getValue()
    + " " + scaleReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "cnorm = "
    + cnorm[ ioffcnorm + 0 ] + " "
    + cnorm[ ioffcnorm + 1 ] + " "
    + cnorm[ ioffcnorm + 2 ]  + "\n";

  for ( i = 0; i < n; i ++ ) x[ ioffx + i ] = i - 1;
  LaPack1.dlatrs('L','N','U','N',n,A,lda,x,scaleReference,cnorm,info,
    ioffa,ioffx,ioffcnorm);
  document.getElementById("debug_textarea").value +=
    "dlatrs('L','N','U','N',...): info,scale = " + info.getValue()
    + " " + scaleReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "cnorm = "
    + cnorm[ ioffcnorm + 0 ] + " "
    + cnorm[ ioffcnorm + 1 ] + " "
    + cnorm[ ioffcnorm + 2 ]  + "\n";
  for ( i = 0; i < n; i ++ ) {
    x[ ioffx + i ] = i - 1;
    cnorm[ ioffcnorm + i ] = 16.;
  }
  LaPack1.dlatrs('L','N','U','Y',n,A,lda,x,scaleReference,cnorm,info,
    ioffa,ioffx,ioffcnorm);
  document.getElementById("debug_textarea").value +=
    "dlatrs('L','N','U','Y',...): info,scale = " + info.getValue()
    + " " + scaleReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "cnorm = "
    + cnorm[ ioffcnorm + 0 ] + " "
    + cnorm[ ioffcnorm + 1 ] + " "
    + cnorm[ ioffcnorm + 2 ]  + "\n";

  for ( i = 0; i < n; i ++ ) x[ ioffx + i ] = i - 1;
  LaPack1.dlatrs('L','N','N','N',n,A,lda,x,scaleReference,cnorm,info,
    ioffa,ioffx,ioffcnorm);
  document.getElementById("debug_textarea").value +=
    "dlatrs('L','N','N','N',...): info,scale = " + info.getValue()
    + " " + scaleReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "cnorm = "
    + cnorm[ ioffcnorm + 0 ] + " "
    + cnorm[ ioffcnorm + 1 ] + " "
    + cnorm[ ioffcnorm + 2 ]  + "\n";
  for ( i = 0; i < n; i ++ ) {
    x[ ioffx + i ] = i - 1;
    cnorm[ ioffcnorm + i ] = 16.;
  }
  LaPack1.dlatrs('L','N','N','Y',n,A,lda,x,scaleReference,cnorm,info,
    ioffa,ioffx,ioffcnorm);
  document.getElementById("debug_textarea").value +=
    "dlatrs('L','N','N','Y',...): info,scale = " + info.getValue()
    + " " + scaleReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "cnorm = "
    + cnorm[ ioffcnorm + 0 ] + " "
    + cnorm[ ioffcnorm + 1 ] + " "
    + cnorm[ ioffcnorm + 2 ]  + "\n";

  for ( i = 0; i < n; i ++ ) x[ ioffx + i ] = i - 1;
  LaPack1.dlatrs('L','T','U','N',n,A,lda,x,scaleReference,cnorm,info,
    ioffa,ioffx,ioffcnorm);
  document.getElementById("debug_textarea").value +=
    "dlatrs('L','T','U','N',...): info,scale = " + info.getValue()
    + " " + scaleReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "cnorm = "
    + cnorm[ ioffcnorm + 0 ] + " "
    + cnorm[ ioffcnorm + 1 ] + " "
    + cnorm[ ioffcnorm + 2 ]  + "\n";
  for ( i = 0; i < n; i ++ ) {
    x[ ioffx + i ] = i - 1;
    cnorm[ ioffcnorm + i ] = 16.;
  }
  LaPack1.dlatrs('L','T','U','Y',n,A,lda,x,scaleReference,cnorm,info,
    ioffa,ioffx,ioffcnorm);
  document.getElementById("debug_textarea").value +=
    "dlatrs('L','T','U','Y',...): info,scale = " + info.getValue()
    + " " + scaleReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "cnorm = "
    + cnorm[ ioffcnorm + 0 ] + " "
    + cnorm[ ioffcnorm + 1 ] + " "
    + cnorm[ ioffcnorm + 2 ]  + "\n";

  for ( i = 0; i < n; i ++ ) x[ ioffx + i ] = i - 1;
  LaPack1.dlatrs('L','T','N','N',n,A,lda,x,scaleReference,cnorm,info,
    ioffa,ioffx,ioffcnorm);
  document.getElementById("debug_textarea").value +=
    "dlatrs('L','T','N','N',...): info,scale = " + info.getValue()
    + " " + scaleReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "cnorm = "
    + cnorm[ ioffcnorm + 0 ] + " "
    + cnorm[ ioffcnorm + 1 ] + " "
    + cnorm[ ioffcnorm + 2 ]  + "\n";
  for ( i = 0; i < n; i ++ ) {
    x[ ioffx + i ] = i - 1;
    cnorm[ ioffcnorm + i ] = 16.;
  }
  LaPack1.dlatrs('L','T','N','Y',n,A,lda,x,scaleReference,cnorm,info,
    ioffa,ioffx,ioffcnorm);
  document.getElementById("debug_textarea").value +=
    "dlatrs('L','T','N','Y',...): info,scale = " + info.getValue()
    + " " + scaleReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "cnorm = "
    + cnorm[ ioffcnorm + 0 ] + " "
    + cnorm[ ioffcnorm + 1 ] + " "
    + cnorm[ ioffcnorm + 2 ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlauum() {
  document.getElementById("debug_textarea").value +=
    "testing dlauum *************" + "\n";
  var n = 3;
  var lda = 4;
  var ioffa = 1;
  var A = new Array( ioffa + lda * n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i <=j; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  var info = new IntReference();
  LaPack1.dlauum('U',n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dlauum('U',...): A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ]  + "\n";

  A = new Array( ioffa + lda * n );
  for ( j = 0; j < n; j ++ ) {
    for ( i = j; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  LaPack1.dlauum('L',n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dlauum('L',...): A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dpoequb() {
  document.getElementById("debug_textarea").value +=
    "testing dpoequb *************" + "\n";
  var n = 3;
  var lda = 4
  var ioffa = 1;
  var ioffs  = 2;
  var A = new Array( ioffa + lda * n );
  var s = new Array( ioffs + n );
  for ( var j = 0; j < n; j ++ ) {
    A[ ioffa + j + j * lda ] = j + 1;
  }
  var info = new IntReference();
  var scondReference = new NumberReference();
  var amaxReference = new NumberReference();
  LaPack1.dpoequb(n,A,lda,s,scondReference,amaxReference,info,ioffa,
    ioffs);
  document.getElementById("debug_textarea").value +=
    "dpoequb: scond,amax = " + scondReference.getValue() + " "
    + amaxReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "s = " + s[ ioffs + 0 ] + " " + s[ ioffs + 1 ] + " "
    + s[ ioffs + 2 ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dpotrf() {
  document.getElementById("debug_textarea").value +=
    "testing dpotrf *************" + "\n";
  var n = 4;
  var lda = 5;
  var ioffa = 1;
  var A = new Array( ioffa + lda * n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i <= j; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + j );
    }
  }
  var info = new IntReference();
  LaPack1.dpotrf('U',n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dpotrf('U',...): A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 0 + 3 * lda ] + " "
    + A[ ioffa + 1 + 3 * lda ] + " "
    + A[ ioffa + 2 + 3 * lda ] + " "
    + A[ ioffa + 3 + 3 * lda ]  + "\n";
  A = new Array( ioffa + lda * n );
  for ( j = 0; j < n; j ++ ) {
    for ( i = j; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + j );
    }
  }
  LaPack1.dpotrf('L',n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dpotrf('L',...): A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 3 + 0 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 3 + 1 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 3 + 2 * lda ] + " "
    + A[ ioffa + 3 + 3 * lda ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dpstf2() {
  document.getElementById("debug_textarea").value +=
    "testing dpstf2 *************" + "\n";
  var n = 4;
  var lda = 5;
  var tol = -1.;
  var ioffa = 1;
  var ioffpiv = 2;
  var ioffwork = 3;
  var A = new Array( ioffa + lda * n );
  var piv = new Array( ioffpiv + n );
  var work = new Array( ioffwork + 2 * n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + j );
    }
  }
  var rank = new IntReference();
  var info = new IntReference();
  LaPack1.dpstf2('U',n,A,lda,piv,rank,tol,work,info,
    ioffa,ioffpiv,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dpstf2('U',...): info,rank = " + info.getValue() + " "
    + rank.getValue()  + "\n";
  document.getElementById("debug_textarea").value +=
    "piv = "
    + piv[ ioffpiv + 0 ] + " "
    + piv[ ioffpiv + 1 ] + " "
    + piv[ ioffpiv + 2 ] + " "
    + piv[ ioffpiv + 3 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 0 + 3 * lda ] + " "
    + A[ ioffa + 1 + 3 * lda ] + " "
    + A[ ioffa + 2 + 3 * lda ] + " "
    + A[ ioffa + 3 + 3 * lda ]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + j );
    }
  }
  LaPack1.dpstf2('L',n,A,lda,piv,rank,tol,work,info,
    ioffa,ioffpiv,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dpstf2('L',...): info,rank = " + info.getValue() + " "
    + rank.getValue()  + "\n";
  document.getElementById("debug_textarea").value +=
    "piv = "
    + piv[ ioffpiv + 0 ] + " "
    + piv[ ioffpiv + 1 ] + " "
    + piv[ ioffpiv + 2 ] + " "
    + piv[ ioffpiv + 3 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 3 + 0 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 3 + 1 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 3 + 2 * lda ] + " "
    + A[ ioffa + 3 + 3 * lda ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dpttrs() {
  document.getElementById("debug_textarea").value +=
    "testing dpttrs *************" + "\n";
  var n = 5;
  var nrhs = 1;
  var ldb = 11;
  var ioffd = 1;
  var ioffe = 2;
  var ioffb = 3;
  var d = new Array( ioffd + 11 );
  var e = new Array( ioffe + 10 );
  var B = new Array( ioffb + ldb * nrhs );
  for ( var i = 0; i < 11; i ++ ) d[ ioffd + i ] = 2 * i + 3;
  for ( i = 0; i < 10; i ++ ) e[ ioffe + i ] = i;
  for ( i = 0; i < n; i ++ ) B[ ioffb + i ] = 1.;
  var info = new IntReference();
  LaPack1.dpttrs(n,nrhs,d,e,B,ldb,info,ioffd,ioffe,ioffb);
  document.getElementById("debug_textarea").value +=
    "B = "
    + B[ ioffb + 0 ] + " "
    + B[ ioffb + 1 ] + " "
    + B[ ioffb + 2 ] + " "
    + B[ ioffb + 3 ] + " "
    + B[ ioffb + 4 ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_drscl() {
  document.getElementById("debug_textarea").value +=
    "testing drscl *************" + "\n";
  var n = 4;
  var incx = 2;
  var ioffsx = 1;
  var sx = new Array( ioffsx + n * incx );
  for ( var i = 0; i < n * incx; i ++ ) sx[ ioffsx + i ] = i;
  var sa = 2.;
  LaPack1.drscl( n, sa, sx, incx, ioffsx );
  document.getElementById("debug_textarea").value +=
    "sx = "
    + sx[ ioffsx + 0 ] + " "
    + sx[ ioffsx + 1 ] + " "
    + sx[ ioffsx + 2 ] + " "
    + sx[ ioffsx + 3 ] + " "
    + sx[ ioffsx + 4 ] + " "
    + sx[ ioffsx + 5 ] + " "
    + sx[ ioffsx + 6 ] + " "
    + sx[ ioffsx + 7 ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dstebz() {
  document.getElementById("debug_textarea").value +=
    "testing dstebz *************" + "\n";
  var n = 64;
  var vl = 1.8;
  var vu = 2.2;
  var il = 10;
  var iu = 12;
  var abstol = 0.;
  var ioffd = 1;
  var ioffe = 2;
  var ioffiblock = 3;
  var ioffisplit = 4;
  var ioffiwork = 5;
  var ioffw = 6;
  var ioffwork = 7;
  var d = new Array( ioffd + n );
  var e = new Array( ioffe + n - 1 );
  var iblock = new Array( ioffiblock + n );
  var isplit = new Array( ioffisplit + n );
  var iwork = new Array( ioffiwork + 3 * n );
  var w = new Array( ioffw + n );
  var work = new Array( ioffwork + 4 * n );
  for ( var i = 0; i < n; i ++ ) d[ ioffd + i ] = 2.;
  for ( i = 0; i < n - 1; i ++ ) e[ ioffe + i ] = -1.;
  e[ ioffe + 32 ] = 0.;
  var m = new IntReference();
  var nsplit = new IntReference();
  var info = new IntReference();
  var j = -1;
  LaPack1.dstebz('A','B',n,vl,vu,il,iu,abstol,d,e,m,nsplit,w,iblock,
    isplit,work,iwork,info,ioffd,ioffe,ioffw,ioffiblock,ioffisplit,
    ioffwork,ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dstebz('A','B',...): m,nsplit,info = " + m.getValue() + " "
    + nsplit.getValue() + " " + info.getValue() + "\n";
  for ( i = 0; i < m.value; i ++ ) {
    j = ioffw + i;
    document.getElementById("debug_textarea").value +=
      "w[" + i  + "] = " + w[ j ] + "\n";
  }
  for ( i = 0; i < m.getValue(); i ++ ) {
    j = ioffiblock + i;
    document.getElementById("debug_textarea").value +=
      "iblock[" + i  + "] = " + iblock[ j ] + "\n";
  }
  for ( i = 0; i < nsplit.getValue(); i ++ ) {
    j = ioffisplit + i;
    document.getElementById("debug_textarea").value +=
      "isplit[" + i  + "] = " + isplit[ j ] + "\n";
  }
  LaPack1.dstebz('A','E',n,vl,vu,il,iu,abstol,d,e,m,nsplit,w,iblock,
    isplit,work,iwork,info,ioffd,ioffe,ioffw,ioffiblock,ioffisplit,
    ioffwork,ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dstebz('A','E',...): m,nsplit,info = " + m.getValue() + " "
    + nsplit.getValue() + " " + info.getValue() + "\n";
  for ( i = 0; i < m.getValue(); i ++ ) {
    j = ioffw + i;
    document.getElementById("debug_textarea").value +=
      "w[" + i  + "] = " + w[ j ] + "\n";
  }
  for ( i = 0; i < m.getValue(); i ++ ) {
    j = ioffiblock + i;
    document.getElementById("debug_textarea").value +=
      "iblock[" + i  + "] = " + iblock[ j ] + "\n";
  }
  for ( i = 0; i < nsplit.getValue(); i ++ ) {
    j = ioffisplit + i;
    document.getElementById("debug_textarea").value +=
      "isplit[" + i  + "] = " + isplit[ j ] + "\n";
  }

  LaPack1.dstebz('V','B',n,vl,vu,il,iu,abstol,d,e,m,nsplit,w,iblock,
    isplit,work,iwork,info,ioffd,ioffe,ioffw,ioffiblock,ioffisplit,
    ioffwork,ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dstebz('V','B',...): m,nsplit,info = " + m.getValue() + " "
    + nsplit.getValue() + " " + info.getValue() + "\n";
  for ( i = 0; i < m.getValue(); i ++ ) {
    j = ioffw + i;
    document.getElementById("debug_textarea").value +=
      "w[" + i  + "] = " + w[ j ] + "\n";
  }
  for ( i = 0; i < m.getValue(); i ++ ) {
    j = ioffiblock + i;
    document.getElementById("debug_textarea").value +=
      "iblock[" + i  + "] = " + iblock[ j ] + "\n";
  }
  for ( i = 0; i < nsplit.getValue(); i ++ ) {
    j = ioffisplit + i;
    document.getElementById("debug_textarea").value +=
      "isplit[" + i  + "] = " + isplit[ j ] + "\n";
  }
  LaPack1.dstebz('V','E',n,vl,vu,il,iu,abstol,d,e,m,nsplit,w,iblock,
    isplit,work,iwork,info,ioffd,ioffe,ioffw,ioffiblock,ioffisplit,
    ioffwork,ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dstebz('V','E',...): m,nsplit,info = " + m.getValue() + " "
    + nsplit.getValue() + " " + info.getValue() + "\n";
  for ( i = 0; i < m.getValue(); i ++ ) {
    j = ioffw + i;
    document.getElementById("debug_textarea").value +=
      "w[" + i  + "] = " + w[ j ] + "\n";
  }
  for ( i = 0; i < m.getValue(); i ++ ) {
    j = ioffiblock + i;
    document.getElementById("debug_textarea").value +=
      "iblock[" + i  + "] = " + iblock[ j ] + "\n";
  }
  for ( i = 0; i < nsplit.getValue(); i ++ ) {
    j = ioffisplit + i;
    document.getElementById("debug_textarea").value +=
      "isplit[" + i  + "] = " + isplit[ j ] + "\n";
  }

  LaPack1.dstebz('I','B',n,vl,vu,il,iu,abstol,d,e,m,nsplit,w,iblock,
    isplit,work,iwork,info,ioffd,ioffe,ioffw,ioffiblock,ioffisplit,
    ioffwork,ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dstebz('I','B',...): m,nsplit,info = " + m.getValue() + " "
    + nsplit.getValue() + " " + info.getValue() + "\n";
  for ( i = 0; i < m.getValue(); i ++ ) {
    j = ioffw + i;
    document.getElementById("debug_textarea").value +=
      "w[" + i  + "] = " + w[ j ] + "\n";
  }
  for ( i = 0; i < m.getValue(); i ++ ) {
    j = ioffiblock + i;
    document.getElementById("debug_textarea").value +=
      "iblock[" + i  + "] = " + iblock[ j ] + "\n";
  }
  for ( i = 0; i < nsplit.getValue(); i ++ ) {
    j = ioffisplit + i;
    document.getElementById("debug_textarea").value +=
      "isplit[" + i  + "] = " + isplit[ j ] + "\n";
  }
  LaPack1.dstebz('I','E',n,vl,vu,il,iu,abstol,d,e,m,nsplit,w,iblock,
    isplit,work,iwork,info,ioffd,ioffe,ioffw,ioffiblock,ioffisplit,
    ioffwork,ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dstebz('I','E',...): m,nsplit,info = " + m.getValue() + " "
    + nsplit.getValue() + " " + info.getValue() + "\n";
  for ( i = 0; i < m.getValue(); i ++ ) {
    j = ioffw + i;
    document.getElementById("debug_textarea").value +=
      "w[" + i  + "] = " + w[ j ] + "\n";
  }
  for ( i = 0; i < m.getValue(); i ++ ) {
    j = ioffiblock + i;
    document.getElementById("debug_textarea").value +=
      "iblock[" + i  + "] = " + iblock[ j ] + "\n";
  }
  for ( i = 0; i < nsplit.getValue(); i ++ ) {
    j = ioffisplit + i;
    document.getElementById("debug_textarea").value +=
      "isplit[" + i  + "] = " + isplit[ j ] + "\n";
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dsyequb() {
  document.getElementById("debug_textarea").value +=
    "testing dsyequb *************" + "\n";
  var n = 3;
  var lda = 4
  var ioffa = 1;
  var ioffs = 2;
  var ioffwork = 3;
  var A = new Array( ioffa + lda * n );
  var s = new Array( ioffs + n );
  var work = new Array( ioffwork + 3 * n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = i + j + 1;
    }
  }
  var info = new IntReference();
  var scondReference = new NumberReference();
  var amaxReference = new NumberReference();
  LaPack1.dsyequb('U',n,A,lda,s,scondReference,amaxReference,work,info,
    ioffa,ioffs,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dsyequb('U',...): scond,amax = " + scondReference.getValue()
    + " " + amaxReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "s = " + s[ ioffs + 0 ] + " " + s[ ioffs + 1 ] + " "
    + s[ ioffs + 2 ]  + "\n";
  LaPack1.dsyequb('L',n,A,lda,s,scondReference,amaxReference,work,info,
    ioffa,ioffs,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dsyequb('L',...): scond,amax = " + scondReference.getValue() + " "
    + amaxReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "s = " + s[ ioffs + 0 ] + " " + s[ ioffs + 1 ] + " "
    + s[ ioffs + 2 ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dsyrfs() {
  document.getElementById("debug_textarea").value +=
    "testing dsyrfs *************" + "\n";
  var n = 4;
  var nrhs = 1;
  var lda = 5;
  var ldaf = 6;
  var ldb = 7;
  var ldx = 8;

  var ioffa = 1;
  var ioffwork = 2;
  var ioffipiv = 3;
  var ioffiwork = 4;
  var ioffferr = 5;
  var ioffberr = 6;
  var ioffaf = 7;
  var ioffb = 8;
  var ioffx = 9;
  var A = new Array( ioffa + lda * n );
  var AF = new Array( ioffaf + lda * n );
  var ipiv = new Array( ioffipiv + n );
  var iwork = new Array( ioffiwork + n );
  var ferr = new Array( ioffferr + nrhs );
  var berr = new Array( ioffberr + nrhs );
  var work = new Array( ioffwork + 3 * n );
  var info = new IntReference();
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  LaPack0.dlacpy('U',n,n,A,lda,AF,ldaf,ioffa,ioffaf);
  LaPack0.dsytf2('U',n,AF,ldaf,ipiv,info,ioffaf,ioffipiv);
  var B = new Array( ioffb + ldb * nrhs );
  var X = new Array( ioffx + ldx * nrhs );
  B[ ioffb + 0 + 0 * ldb ] = -18.;
  B[ ioffb + 1 + 0 * ldb ] = -19.;
  B[ ioffb + 2 + 0 * ldb ] = -22.;
  B[ ioffb + 3 + 0 * ldb ] = -22.;
  Blas1.dcopy(n,B,1,X,1,ioffb,ioffx);
  LaPack0.dsytrs('U',n,nrhs,AF,ldaf,ipiv,X,ldx,info,
    ioffaf,ioffipiv,ioffx);
  LaPack1.dsyrfs('U',n,nrhs,A,lda,AF,ldaf,ipiv,B,ldb,X,ldx,
    ferr,berr,work,iwork,info,
    ioffa,ioffaf,ioffipiv,ioffb,ioffx,ioffferr,ioffberr,ioffwork,
    ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dsyrfs('U',...): berr,ferr = " + berr[ ioffberr ] + " "
    + ferr[ ioffferr ] + "\n";
  document.getElementById("debug_textarea").value +=
    "X = "
    + X[ ioffx + 0 + 0 * ldx ] + " "
    + X[ ioffx + 1 + 0 * ldx ] + " "
    + X[ ioffx + 2 + 0 * ldx ] + " "
    + X[ ioffx + 3 + 0 * ldx ]  + "\n";

  LaPack0.dlacpy('L',n,n,A,lda,AF,ldaf,ioffa,ioffaf);
  LaPack0.dsytf2('L',n,AF,ldaf,ipiv,info,ioffaf,ioffipiv);
  B[ ioffb + 0 + 0 * ldb ] = -10.;
  B[ ioffb + 1 + 0 * ldb ] = -15.;
  B[ ioffb + 2 + 0 * ldb ] = -18.;
  B[ ioffb + 3 + 0 * ldb ] = -24.;
  Blas1.dcopy(n,B,1,X,1,ioffb,ioffx);
  LaPack0.dsytrs('L',n,nrhs,AF,ldaf,ipiv,X,ldx,info,
    ioffaf,ioffipiv,ioffx);
  LaPack1.dsyrfs('L',n,nrhs,A,lda,AF,ldaf,ipiv,B,ldb,X,ldx,
    ferr,berr,work,iwork,info,
    ioffa,ioffaf,ioffipiv,ioffb,ioffx,ioffferr,ioffberr,ioffwork,
    ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dsyrfs('L',...): berr,ferr = " + berr[ ioffberr ] + " "
    + ferr[ ioffferr ] + "\n";
  document.getElementById("debug_textarea").value +=
    "X = "
    + X[ ioffx + 0 + 0 * ldx ] + " "
    + X[ ioffx + 1 + 0 * ldx ] + " "
    + X[ ioffx + 2 + 0 * ldx ] + " "
    + X[ ioffx + 3 + 0 * ldx ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dsytrf() {
  document.getElementById("debug_textarea").value +=
    "testing dsytrf *************" + "\n";
  var n = 3;
  var lda = 5;
  var lwork = 16;
  var ioffa = 1;
  var ioffipiv = 2;
  var ioffwork = 3;
  var A = new Array( ioffa + lda * n );
  var ipiv = new Array( ioffipiv + n );
  var work = new Array( ioffwork + lwork );
  var info = new IntReference();
  A[ ioffa + 0 + 0 * lda ] = 4.;
  A[ ioffa + 1 + 0 * lda ] = 10.;
  A[ ioffa + 2 + 0 * lda ] = -2.;
  A[ ioffa + 0 + 1 * lda ] = 10.;
  A[ ioffa + 1 + 1 * lda ] = 2.;
  A[ ioffa + 2 + 1 * lda ] = 10.;
  A[ ioffa + 0 + 2 * lda ] = -2.;
  A[ ioffa + 1 + 2 * lda ] = 10.;
  A[ ioffa + 2 + 2 * lda ] = 4.;
  LaPack1.dsytrf('U',n,A,lda,ipiv,work,lwork,info,
    ioffa,ioffipiv,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dsytrf('U',...): info = " + info.getValue() + "\n";
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
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ]  + "\n";

  A[ ioffa + 0 + 0 * lda ] = 4.;
  A[ ioffa + 1 + 0 * lda ] = 10.;
  A[ ioffa + 2 + 0 * lda ] = -2.;
  A[ ioffa + 0 + 1 * lda ] = 10.;
  A[ ioffa + 1 + 1 * lda ] = 2.;
  A[ ioffa + 2 + 1 * lda ] = 10.;
  A[ ioffa + 0 + 2 * lda ] = -2.;
  A[ ioffa + 1 + 2 * lda ] = 10.;
  A[ ioffa + 2 + 2 * lda ] = 4.;
  LaPack1.dsytrf('L',n,A,lda,ipiv,work,lwork,info,
    ioffa,ioffipiv,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dsytrf('L',...): info = " + info.getValue() + "\n";
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
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dsytrs2() {
  document.getElementById("debug_textarea").value +=
    "testing dsytrs2 *************" + "\n";
  var n = 3;
  var nrhs = 1;
  var lda = 4;
  var ldb = 5;
  var ioffa = 1;
  var ioffb = 2;
  var ioffipiv = 3;
  var ioffwork = 4;
  var A = new Array( ioffa + lda * n );
  var B = new Array( ioffb + ldb * nrhs );
  var ipiv = new Array( ioffipiv + n );
  var work = new Array( ioffwork + n );
  var info = new IntReference();
  A[ ioffa + 0 + 0 * lda ] = 4.;
  A[ ioffa + 1 + 0 * lda ] = 10.;
  A[ ioffa + 2 + 0 * lda ] = -2.;
  A[ ioffa + 0 + 1 * lda ] = 10.;
  A[ ioffa + 1 + 1 * lda ] = 2.;
  A[ ioffa + 2 + 1 * lda ] = 10.;
  A[ ioffa + 0 + 2 * lda ] = -2.;
  A[ ioffa + 1 + 2 * lda ] = 10.;
  A[ ioffa + 2 + 2 * lda ] = 4.;
  B[ ioffb + 0 + 0 * ldb ] = 18.;
  B[ ioffb + 1 + 0 * ldb ] = 44.;
  B[ ioffb + 2 + 0 * ldb ] = 30.;
  LaPack1.dsytrs2('U',n,nrhs,A,lda,ipiv,B,ldb,work,info,
    ioffa, ioffipiv,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dsytrs2('U',...): B = "
    + B[ ioffb + 0 ] + " "
    + B[ ioffb + 1 ] + " "
    + B[ ioffb + 2 ]  + "\n";

  A[ ioffa + 0 + 0 * lda ] = 4.;
  A[ ioffa + 1 + 0 * lda ] = 10.;
  A[ ioffa + 2 + 0 * lda ] = -2.;
  A[ ioffa + 0 + 1 * lda ] = 10.;
  A[ ioffa + 1 + 1 * lda ] = 2.;
  A[ ioffa + 2 + 1 * lda ] = 10.;
  A[ ioffa + 0 + 2 * lda ] = -2.;
  A[ ioffa + 1 + 2 * lda ] = 10.;
  A[ ioffa + 2 + 2 * lda ] = 4.;
  B[ ioffb + 0 + 0 * ldb ] = 18.;
  B[ ioffb + 1 + 0 * ldb ] = 44.;
  B[ ioffb + 2 + 0 * ldb ] = 30.;
  LaPack1.dsytrs2('L',n,nrhs,A,lda,ipiv,B,ldb,work,info,
    ioffa, ioffipiv,ioffb,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dsytrs2('L',...): B = "
    + B[ ioffb + 0 ] + " "
    + B[ ioffb + 1 ] + " "
    + B[ ioffb + 2 ]  + "\n";
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dtrtri() {
  document.getElementById("debug_textarea").value +=
    "testing dtrtri *************" + "\n";
  var n = 4;
  var lda = 5;
  var ioffa = 1;
  var i = -1;
  var j = -1;
  var info = new IntReference();
  var A = new Array( ioffa + lda * n );
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  LaPack1.dtrtri('U','N',n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dtrtri('U','N',...): A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 0 + 3 * lda ] + " "
    + A[ ioffa + 1 + 3 * lda ] + " "
    + A[ ioffa + 2 + 3 * lda ] + " "
    + A[ ioffa + 3 + 3 * lda ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  LaPack1.dtrtri('U','U',n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dtrtri('U','U',...): A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 0 + 3 * lda ] + " "
    + A[ ioffa + 1 + 3 * lda ] + " "
    + A[ ioffa + 2 + 3 * lda ] + " "
    + A[ ioffa + 3 + 3 * lda ]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  LaPack1.dtrtri('L','N',n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dtrtri('L','N',...): A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 3 + 0 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 3 + 1 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 3 + 2 * lda ] + " "
    + A[ ioffa + 3 + 3 * lda ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  LaPack1.dtrtri('L','U',n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dtrtri('L','U',...): A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 3 + 0 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 3 + 1 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + " "
    + A[ ioffa + 3 + 2 * lda ] + " "
    + A[ ioffa + 3 + 3 * lda ]  + "\n";
}
