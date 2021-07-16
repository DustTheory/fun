function LaPack0() {
  this.dlacon_iter = 0;
  this.dlacon_itmax = 5;
  this.dlacon_j = -1;
  this.dlacon_jump = 0;
  this.dlamc1_first = undefined;
  this.dlamc1_lbeta = -1;
  this.dlamc1_lieee1 = true;
  this.dlamc1_lt = -1;
  this.dlamc1_lrnd = false;
  this.dlamc2_first = undefined;
  this.dlamc2_iwarn = false;
  this.dlamc2_lbeta = undefined;
  this.dlamc2_lemax = undefined;
  this.dlamc2_lemin = -1;
  this.dlamc2_leps = Number.POSITIVE_INFINITY;
  this.dlamc2_lrmaxReference = undefined;
  this.dlamc2_lrmin = Number.POSITIVE_INFINITY;
  this.dlamc2_lt = undefined;
  this.dlamc2_lrndReference = undefined;
  this.dlamch_first = undefined;
  this.dlamch_epsReference = undefined;
  this.dlamch_sfmin = Number.POSITIVE_INFINITY;
  this.dlamch_base = Number.POSITIVE_INFINITY;
  this.dlamch_t = Number.POSITIVE_INFINITY;
  this.dlamch_rnd = Number.POSITIVE_INFINITY;
  this.dlamch_emin = Number.POSITIVE_INFINITY;
  this.dlamch_rminReference = undefined;
  this.dlamch_emax = Number.POSITIVE_INFINITY;
  this.dlamch_rmaxReference = undefined;
  this.dlamch_prec = Number.POSITIVE_INFINITY;
}
//*************************************************************************
LaPack0.chla_transtype = function( trans ) {
  var blas_no_trans = 111;
  var blas_trans = 112;
  var blas_conj_trans = 113;
  if ( trans == blas_no_trans ) return 'N';
  else if ( trans == blas_trans ) return 'T';
  else if ( trans == blas_conj_trans ) return 'C';
  else return 'X';
}
//*************************************************************************
LaPack0.dgbtf2 = function( m, n, kl, ku, AB, ldab, ipiv, info ) {
  throw new Error("not programmed: band matrix");
}
LaPack0.zgbtf2 = function( ) {
  throw new Error("not programmed: band matrix");
}
//*************************************************************************
LaPack0.dgbtrs = function( ) {
  throw new Error("not programmed: band matrix");
}
LaPack0.zgbtrs = function( ) {
  throw new Error("not programmed: band matrix");
}
//*************************************************************************
LaPack0.dgebak = function( job, side, n, ilo, ihi, scal, m, V,
ldv, info, ioffscal, ioffv ) {
  var rightv = side.charAt(0).toUpperCase() == 'R';
  var leftv = side.charAt(0).toUpperCase() == 'L';
  info.setValue( 0 );
  if ( job.charAt(0).toUpperCase() != 'N' &&
  job.charAt(0).toUpperCase() != 'P' &&
  job.charAt(0).toUpperCase() != 'S' &&
  job.charAt(0).toUpperCase() != 'B' ) {
    info.setValue( -1 );
  } else if ( ! rightv && ! leftv ) info.setValue( -2 );
  else if ( n < 0 ) info.setValue( -3 );
  else if ( ilo < 1 || ilo > Math.max( 1 , n ) ) info.setValue( -4 );
  else if ( ihi < Math.min( ilo, n ) || ihi > n ) info.setValue( - 5 );
  else if ( m < 0 ) info.setValue( -7 );
  else if ( ldv < Math.max( 1 , n ) ) info.setValue( -9 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgebak', - info.getValue() );
    return;
  }

  if (n == 0 || m == 0 || job.charAt(0).toUpperCase() == 'N') return;
  var i = -1;
  var j = -1;
  var s = Number.POSITIVE_INFINITY;
  if ( ilo != ihi ) {
    if ( job.charAt(0).toUpperCase() == 'S' ||
    job.charAt(0).toUpperCase() == 'B' ) {
      if ( rightv ) {
        for ( i = ilo; i <= ihi; i ++ ) {
          s = scal[ ioffscal + i - 1 ];
          Blas1.dscal( m, s, V, ldv, ioffv + i - 1 );
        }
      }
      if ( leftv ) {
        for ( i = ilo; i <= ihi; i ++ ) {
          s = 1. / scal[ ioffscal + i - 1 ];
          Blas1.dscal( m, s, V, ldv, ioffv + i - 1 );
        }
      }
    }
  }
  if ( job.charAt(0).toUpperCase() == 'P' ||
  job.charAt(0).toUpperCase() == 'B' ) {
    var ii = -1;
    var k = -1;
    var temp = Number.POSITIVE_INFINITY;
    if ( rightv ) {
      for ( ii = 1; ii <= n; ii ++ ) {
        i = ii;
        if ( i >= ilo && i <= ihi ) continue;
        if ( i < ilo ) i = ilo - ii;
        k = scal[ ioffscal + i - 1 ];
        if ( k == i ) continue;
        Blas1.dswap( m, V, ldv, V, ldv, ioffv + i - 1, ioffv + k - 1 );
      }
    }
    if ( leftv ) {
      for ( ii = 1; ii <= n; ii ++ ) {
        i = ii;
        if ( i >= ilo && i <= ihi ) continue;
        if ( i < ilo ) i = ilo - ii;
        k = scal[ ioffscal + i - 1 ];
        if ( k == i ) continue;
        Blas1.dswap( m, V, ldv, V, ldv, ioffv + i - 1, ioffv + k - 1 );
      }
    }
  }
}
LaPack0.zgebak = function( job, side, n, ilo, ihi, scal, m, V,
ldv, info) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dggbak = function( job, side, n, ilo, ihi, lscale,
rscale, m, V, ldv, info, iofflscale, ioffrscale, ioffv ) {
  var rightv = side.charAt(0).toUpperCase() == 'R';
  var leftv = side.charAt(0).toUpperCase() == 'L';
  info.setValue( 0 );
  if ( job.charAt(0).toUpperCase() != 'N' &&
  job.charAt(0).toUpperCase() != 'P' &&
  job.charAt(0).toUpperCase() != 'S' &&
  job.charAt(0).toUpperCase() != 'B' ) {
    info.setValue( -1 );
  } else if ( ! rightv && ! leftv ) info.setValue( -2 );
  else if ( n < 0 ) info.setValue( -3 );
  else if ( ilo < 1 ) info.setValue( -4 );
  else if ( ihi < ilo || ihi > Math.max( 1, n ) ) info.setValue( -5 );
  else if ( m < 0 ) info.setValue( -6 );
  else if ( ldv < Math.max( 1 , n ) ) info.setValue( -10 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dggbak', - info.getValue() );
    return;
  }

  if ( n == 0 || m == 0 || job.charAt(0).toUpperCase() == 'N') {
    return;
  }
  if ( ilo != ihi ) {
    if ( job.charAt(0).toUpperCase() == 'S' ||
    job.charAt(0).toUpperCase() == 'B' ) {
      var i = -1;
      var j = -1;
      var s = Number.POSITIVE_INFINITY;
      if ( rightv ) {
        for ( i = ilo; i <= ihi; i ++ ) {
          Blas1.dscal( m, rscale[ ioffrscale + i - 1 ], V, ldv,
            ioffv + i - 1 );
        }
      }
      if ( leftv ) {
        for ( i = ilo; i <= ihi; i ++ ) {
          Blas1.dscal( m, lscale[ iofflscale + i - 1 ], V, ldv,
            ioffv + i - 1 );
        }
      }
    }
  }
  if ( job.charAt(0).toUpperCase() == 'P' ||
  job.charAt(0).toUpperCase() == 'B' ) {
    var k = -1;
    var temp = Number.POSITIVE_INFINITY;
    if ( rightv ) {
      if ( ilo != 1 ) {
        for ( i = ilo - 1; i >= 1; i -- ) {
          k = rscale[ ioffrscale + i - 1 ];
          if ( k != i ) {
            Blas1.dswap( m, V, ldv, V, ldv, ioffv + i - 1,
              ioffv + k - 1 );
          }
        }
      }
      if ( ihi != n ) {
        for ( i = ihi + 1; i <= n; i ++ ) {
          k = rscale[ ioffrscale + i - 1 ];
          if ( k != i ) {
            Blas1.dswap( m, V, ldv, V, ldv, ioffv + i - 1,
              ioffv + k - 1 );
          }
        }
      }
    }
    if ( leftv ) {
      if ( ilo != 1 ) {
        for ( i = ilo - 1; i >= 1; i -- ) {
          k = lscale[ iofflscale + i - 1 ];
          if ( k != i ) {
            Blas1.dswap( m, V, ldv, V, ldv, ioffv + i - 1,
              ioffv + k - 1 );
          }
        }
      }
      if ( ihi != n ) {
        for ( i = ihi + 1; i <= n; i ++ ) {
          k = lscale[ iofflscale + i - 1 ];
          if ( k != i ) {
            Blas1.dswap( m, V, ldv, V, ldv, ioffv + i - 1,
              ioffv + k - 1 );
          }
        }
      }
    }
  }
}
LaPack0.zggbak = function( job, side, n, ilo, ihi, lscale,
rscale, m, V, ldv, info ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dgttrf = function( n, dl, d, du, du2, ipiv, info,
ioffdl, ioffd, ioffdu, ioffdu2, ioffipiv ) {
  info.setValue( 0 );
  if ( n < 0 ) {
    info.setValue( -1 );
    Blas2.xerbla( 'dgttrf', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  for ( var i = 1; i <= n; i ++ ) ipiv[ ioffipiv + i - 1 ] = i;
  for ( i = 1; i <= n - 2; i ++ ) du2[ ioffdu2 + i - 1 ] = 0.;
  for ( i = 1; i <= n - 2; i ++ ) {
    var fact = Number.POSITIVE_INFINITY;
    if ( Math.abs( d[ ioffd + i - 1 ] ) >=
    Math.abs( dl[ ioffdl + i - 1 ] ) ) {
      if ( d[ ioffd + i - 1 ] != 0. ) {
        fact = dl[ ioffdl + i - 1 ] / d[ ioffd + i - 1 ];
        dl[ ioffdl + i - 1 ] = fact;
        d[ ioffd + i ] -= fact * du[ ioffdu + i - 1 ];
      }
    } else {
      fact = d[ ioffd + i - 1 ] / dl[ ioffdl + i - 1 ];
      d[ ioffd + i - 1 ] = dl[ ioffdl + i - 1 ];
      dl[ ioffdl + i - 1 ] = fact;
      var temp = du[ ioffdu + i - 1 ];
      du[ ioffdu + i - 1 ] = d[ ioffd + i ];
      d[ ioffd + i ] = temp - fact * d[ ioffd + i ];
      du2[ ioffdu2 + i - 1 ] = du[ ioffdu + i ];
      du[ ioffdu + i ] = - fact * du[ ioffdu + i ];
      ipiv[ ioffipiv + i - 1 ] = i + 1;
    }
  }
  if ( n > 1 ) {
    i = n - 1;
    if ( Math.abs( d[ ioffd + i - 1 ] ) >=
    Math.abs( dl[ ioffdl + i - 1 ] ) ) {
      if ( d[ ioffd + i - 1 ] != 0. ) {
        fact = dl[ ioffdl + i - 1 ] / d[ ioffd + i - 1 ];
        dl[ ioffdl + i - 1 ] = fact;
        d[ ioffd + i ] -= fact * du[ ioffdu + i - 1 ];
      }
    } else {
      fact = d[ ioffd + i - 1 ] / dl[ ioffdl + i - 1 ];
      d[ ioffd + i - 1 ] = dl[ ioffdl + i - 1 ];
      dl[ ioffdl + i - 1 ] = fact;
      temp = du[ ioffdu + i - 1 ];
      du[ ioffdu + i - 1 ] = d[ ioffd + i ];
      d[ ioffd + i ] = temp - fact * d[ ioffd + i ];
      ipiv[ ioffipiv + i - 1 ] = i + 1;
    }
  }
  for ( i = 1; i <= n; i ++ ) {
    if ( d[ ioffd + i - 1 ] == 0. ) {
      info.setValue( i );
      return;
    }
  }
}
LaPack0.zgttrf = function( n, dl, d, du, ipiv, info ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dgtts2 = function( itrans, n, nrhs, dl, d, du, du2,
ipiv, B, ldb, ioffdl, ioffd, ioffdu, ioffdu2, ioffipiv, ioffb ) {
  if ( n == 0 || nrhs == 0 ) return;
  if ( itrans == 0 ) {
    if ( nrhs <= 1 ) {
      var j = 1;
      while ( j <= nrhs ) { // 10
        for ( var i = 1; i <= n - 1; i ++ ) {
          var ip = ipiv[ ioffipiv + i - 1 ];
          var temp = B[ ioffb + i - ip + i + ( j - 1 ) * ldb ]
            - dl[ ioffdl + i - 1 ]
            * B[ ioffb + ip - 1 + ( j - 1 ) * ldb ];
          B[ ioffb + i - 1 + ( j - 1 ) * ldb ] =
            B[ ioffb + ip - 1 + ( j - 1 ) * ldb ];
          B[ ioffb + i + ( j - 1 ) * ldb ] = temp;
        } // 20
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
          - du[ ioffdu + i - 1 ]
          * B[ ioffb + i + ( j - 1 ) * ldb ]
          - du2[ ioffdu2 + i - 1 ]
          * B[ ioffb + i + 1 + ( j - 1 ) * ldb ])
          / d[ ioffd + i - 1 ];
        }
        j ++;
      }
    } else {
      for ( j = 1; j <= nrhs; j ++ ) {
        for ( i = 1; i <= n - 1; i ++ ) {
          if ( ipiv[ ioffipiv + i - 1 ] == i ) {
            B[ ioffb + i + ( j - 1 ) * ldb ] -=
              dl[ ioffdl + i - 1 ]
              * B[ ioffb + i - 1 + ( j - 1 ) * ldb ];
          } else {
            temp = B[ ioffb + i - 1 + ( j - 1 ) * ldb ];
            B[ ioffb + i - 1 + ( j - 1 ) * ldb ] =
              B[ ioffb + i + ( j - 1 ) * ldb ];
            B[ ioffb + i + ( j - 1 ) * ldb ] = temp
              - dl[ ioffdl + i - 1 ]
              * B[ ioffb + i - 1 + ( j - 1 ) * ldb ]
          }
        } // 40
        B[ ioffb + n -1 + ( j - 1 ) * ldb ] /= d[ ioffd + n - 1 ];
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
            - du[ ioffdu + i - 1 ]
            * B[ ioffb + i + ( j - 1 ) * ldb ]
            - du2[ ioffdu2 + i - 1 ]
            * B[ ioffb + i + 1 + ( j - 1 ) * ldb ] )
            / d[ ioffd + i - 1 ];
        } // 50
      } // 60
    }
  } else {
    if ( nrhs <= 1 ) {
      j = 1;
      while ( j <= nrhs ) { // 70
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
          - du[ ioffdu + i - 2 ]
          * B[ ioffb + i - 2 + ( j - 1 ) * ldb ]
          - du2[ ioffdu2 + i - 3 ]
          * B[ ioffb + i - 3 + ( j - 1 ) * ldb ])
          / d[ ioffd + i - 1 ];
        } // 80
        for ( i = n - 1; i >= 1; i -- ) {
          ip = ipiv[ ioffipiv + i - 1 ];
          temp = B[ ioffb + i - 1 + ( j - 1 ) * ldb ]
            - dl[ ioffdl + i - 1 ]
            * B[ ioffb + i + ( j - 1 ) * ldb ];
          B[ ioffb + i - 1 + ( j - 1 ) * ldb ] =
            B[ ioffb + ip - 1 + ( j - 1 ) * ldb ];
          B[ ioffb + ip - 1 + ( j - 1 ) * ldb ] = temp;
        } // 90
        j ++;
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
            - du[ ioffdu + i - 2 ]
            * B[ ioffb + i - 2 + ( j - 1 ) * ldb ]
            - du2[ ioffdu2 + i - 3 ]
            * B[ ioffb + i - 3 + ( j - 1 ) * ldb ] )
            / d[ ioffd + i - 1 ];
        } // 100
        for ( i = n - 1; i >= 1; i -- ) {
          if ( ipiv[ ioffipiv + i - 1 ] == i ) {
            B[ ioffb + i + ( j - 1 ) * ldb ] -=
              dl[ ioffdl + i - 1 ]
              * B[ ioffb + i + ( j - 1 ) * ldb ];
          } else {
            temp = B[ ioffb + i + ( j - 1 ) * ldb ];
            B[ ioffb + i + ( j - 1 ) * ldb ] = 
              B[ ioffb + i - 1 + ( j - 1 ) * ldb ]
              - dl[ ioffdl + i - 1 ] * temp;
            B[ ioffb + i - 1 + ( j - 1 ) * ldb ] = temp;
          }
        } // 110
      } // 120
    }
  }
}
//*************************************************************************
LaPack0.dla_gbrpvgrw = function( n, kl, ku, ncols, AB, ldab,
AFB, ldafb, ioffab, ioffafb) {
  throw new Error("not programmed: banded matrix");
}
//*************************************************************************
LaPack0.dla_gerpvgrw = function( n, ncols, A, lda, AF, ldaf,
ioffa, ioffaf ) {
  var rpvgrw = 1.;
  for ( var j = 1; j <= ncols; j ++ ) {
    var amax = 0.;
    var umax = 0.;
    for ( var i = 1; i <= n; i ++ ) {
      amax = Math.max(
        Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] ), amax );
    }
    for ( i = 1; i <= j; i ++ ) {
      umax = Math.max(
        Math.abs( AF[ ioffaf + i - 1 + ( j - 1 ) * ldaf ] ), umax );
    }
    if ( umax != 0. ) {
      rpvgrw = Math.min( amax / umax, rpvgrw );
    }
  }
  return rpvgrw;
}
//*************************************************************************
LaPack0.dla_porpvgrw = function( uplo, ncols, A, lda, AF, ldaf,
work, ioffa, ioffaf, ioffwork ) {
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var rpvgrw = 1.;
  for ( var i = 1; i <= 2 * ncols; i ++ ) {
    work[ ioffwork + i - 1 ] = 0.;
  }
  if ( upper ) {
    for ( var j = 1; j <= ncols; j ++ ) {
      for ( i = 1; i <= j; i ++ ) {
        work[ ioffwork + ncols + j - 1 ] =
          Math.max( Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] ),
          work[ ioffwork + ncols + j - 1 ] );
      }
    }
  } else {
    for ( j = 1; j <= ncols; j ++ ) {
      for ( i = j; i <= ncols; i ++ ) {
        work[ ioffwork + ncols + j - 1 ] =
          Math.max( Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] ),
          work[ ioffwork + ncols + j - 1 ] );
      }
    }
  }
  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    for ( j = 1; j <= ncols; j ++ ) {
      for ( i = 1; i <= j; i ++ ) {
        work[ ioffwork + j - 1 ] = Math.max(
          Math.abs( AF[ ioffaf + i - 1 + ( j - 1 ) * ldaf ] ),
          work[ ioffwork + j - 1 ] );
      }
    }
  } else {
    for ( j = 1; j <= ncols; j ++ ) {
      for ( i = j; i <= ncols; i ++ ) {
        work[ ioffwork + j - 1 ] = Math.max(
          Math.abs( AF[ ioffaf + i - 1 + ( j - 1 ) * ldaf ] ),
          work[ ioffwork + j - 1 ] );
      }
    }
  }
  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    for ( i = 1; i <= ncols; i ++ ) {
      var umax = work[ ioffwork + i - 1 ];
      var amax = work[ ioffwork + ncols + i - 1 ];
      if ( umax != 0. ) rpvgrw = Math.min( amax / umax, rpvgrw );
    }
  } else {
    for ( i = 1; i <= ncols; i ++ ) {
      umax = work[ ioffwork + i - 1 ];
      amax = work[ ioffwork + ncols + i - 1 ];
      if ( umax != 0. ) rpvgrw = Math.min( amax / umax, rpvgrw );
    }
  }
  return rpvgrw;
}
//*************************************************************************
LaPack0.dla_syrpvgrw = function( uplo, n, info, A, lda, AF,
ldaf, ipiv, work, ioffa, ioffaf, ioffipiv, ioffwork ) {
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( info == 0 ) {
    var ncols = ( upper ? 1 : n );
  } else ncols = info;
  var rpvgrw = 1.;
  for ( var i = 1; i <= 2 * n; i ++ ) {
    work[ ioffwork + i - 1 ] = 0.;
  }
  if ( upper ) {
    for ( var j = 1; j <= n; j ++ ) {
      for ( i = 1; i <= j; i ++ ) {
        work[ ioffwork + n + i - 1 ] = Math.max(
          Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] ),
          work[ ioffwork + n + i - 1 ] );
        work[ ioffwork + n + j - 1 ] = Math.max(
          Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] ),
          work[ ioffwork + n + j - 1 ] );
      }
    }
  } else {
    for ( j = 1; j <= n; j ++ ) {
      for ( i = j; i <= n; i ++ ) {
        work[ ioffwork + n + i - 1 ] = Math.max(
          Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] ),
          work[ ioffwork + n + i - 1 ] );
        work[ ioffwork + n + j - 1 ] = Math.max(
          Math.abs( A[ ioffa + i - 1 + ( j - 1 ) * lda ] ),
          work[ ioffwork + n + j - 1 ] );
      }
    }
  }
  if ( upper ) {
    var k = n;
    while ( k < ncols && k > 0 ) {
      if ( ipiv[ ioffipiv + k - 1 ] > 0 ) {
        var kp = ipiv[ ioffipiv + k - 1 ];
        if ( kp != k ) {
          var tmp = work[ ioffwork + n + k - 1 ];
          work[ ioffwork + n + k - 1 ] = work[ ioffwork + n + kp - 1 ];
          work[ ioffwork + n + kp - 1 ] = tmp;
        }
        for ( i = 1; i <= k; i ++ ) {
          work[ ioffwork + k - 1 ] = Math.max(
            Math.abs( AF[ ioffaf + i - 1 + ( k - 1 ) * ldaf ] ),
            work[ ioffwork + k - 1 ] );
        }
        k --;
      } else {
        kp = - ipiv[ ioffipiv + k - 1 ];
        tmp = work[ ioffwork + n + k - 2 ];
        work[ ioffwork + n + k - 2 ] = work[ ioffwork + n + kp - 1 ];
        work[ ioffwork + n + kp - 1 ] = tmp;
        for ( i = 1; i <= k - 1; i ++ ) {
          work[ ioffwork + k - 1 ] = Math.max(
            Math.abs( AF[ ioffaf + i - 1 + ( k - 1 ) * ldaf ] ),
            work[ ioffwork + k - 1 ] );
          work[ ioffwork + k - 2 ] = Math.max(
            Math.abs( AF[ ioffaf + i - 1 + ( k - 2 ) * ldaf ] ),
            work[ ioffwork + k - 2 ] );
        }
        work[ ioffwork + k - 1 ] = Math.max(
          Math.abs( AF[ ioffaf + k - 1 + ( k - 1 ) * ldaf ] ),
          work[ ioffwork + k - 1 ] );
        k -= 2;
      }
    }
    k = ncols;
    while ( k <= n ) {
      if ( ipiv[ ioffipiv + k - 1 ] > 0 ) {
        kp = ipiv[ ioffipiv + k - 1 ];
        if ( kp != k ) {
          tmp = work[ ioffwork + n + k - 1 ];
          work[ ioffwork + n + k - 1 ] = work[ ioffwork + n + kp - 1 ];
          work[ ioffwork + n + kp - 1 ] = tmp;
        }
        k ++;
      } else {
        kp = - ipiv[ ioffipiv + k - 1 ];
        tmp = work[ ioffwork + n + k - 1 ];
        work[ ioffwork + n + k - 1 ] = work[ ioffwork + n + kp - 1 ];
        work[ ioffwork + n + kp - 1 ] = tmp;
        k += 2;
      }
    }
  } else {
    k = 1;
    while ( k <= ncols ) {
      if ( ipiv[ ioffipiv + k - 1 ] > 0 ) {
        kp = ipiv[ ioffipiv + k - 1 ];
        if ( kp != k ) {
          tmp = work[ ioffwork + n + k - 1 ];
          work[ ioffwork + n + k - 1 ] = work[ ioffwork + n + kp - 1 ];
          work[ ioffwork + n + kp - 1 ] = tmp;
        }
        for ( i = k; i <= n; i ++ ) {
          work[ ioffwork + k - 1 ] = Math.max(
            Math.abs( AF[ ioffaf + i - 1 + ( k - 1 ) * ldaf ] ),
            work[ ioffwork + k - 1 ] );
        }
        k ++;
      } else {
        kp = - ipiv[ ioffipiv + k - 1 ];
        tmp = work[ ioffwork + n + k ];
        work[ ioffwork + n + k ] = work[ ioffwork + n + kp - 1 ];
        work[ ioffwork + n + kp - 1 ] = tmp;
        for ( i = k + 1; i <= n; i ++ ) {
          work[ ioffwork + k - 1 ] = Math.max(
            Math.abs( AF[ ioffaf + i - 1 + ( k - 1 ) * ldaf ] ),
            work[ ioffwork + k - 1 ] );
          work[ ioffwork + k ] = Math.max(
            Math.abs( AF[ ioffaf + i - 1 + k * ldaf ] ),
            work[ ioffwork + k ] );
        }
        work[ ioffwork + k - 1 ] = Math.max(
          Math.abs( AF[ ioffaf + k - 1 + ( k - 1 ) * ldaf ] ),
          work[ ioffwork + k - 1 ] );
        k += 2;
      }
    }
    k = ncols;
    while ( k >= 1 ) {
      if ( ipiv[ ioffipiv + k - 1 ] > 0 ) {
        var ip = ipiv[ ioffipiv + k - 1 ];
        if ( kp != k ) {
          tmp = work[ ioffwork + n + k - 1 ];
          work[ ioffwork + n + k - 1 ] = work[ ioffwork + n + kp - 1 ];
          work[ ioffwork + n + kp - 1 ] = tmp;
        }
        k --;
      } else {
        kp = - ipiv[ ioffipiv + k - 1 ];
        tmp = work[ ioffwork + n + k - 1 ];
        work[ ioffwork + n + k - 1 ] = work[ ioffwork + n + kp - 1 ];
        work[ ioffwork + n + kp - 1 ] = tmp;
        k -= 2;
      }
    }
  }
  if ( upper ) {
    for ( i = ncols; i <= n; i ++ ) {
      var umax = work[ ioffwork + i - 1 ];
      var amax = work[ ioffwork + n + i - 1 ];
      if ( umax != 0. ) rpvgrw = Math.min( amax / umax, rpvgrw );
    }
  } else {
    for ( i = 1; i <= ncols; i ++ ) {
      umax = work[ ioffwork + i - 1 ];
      amax = work[ ioffwork + n + i - 1 ];
      if ( umax != 0. ) rpvgrw = Math.min( amax / umax, rpvgrw );
    }
  }
  return rpvgrw;
}
//*************************************************************************
LaPack0.dla_wwaddw = function( n, x, y, w, ioffx, ioffy, ioffw) {
  for ( var i = 1; i <= n; i ++ ) {
    var s = x[ ioffx + i - 1 ] + w[ ioffw + i - 1 ];
    s = ( s + s ) - s;
    y[ ioffy + i - 1 ] = ( ( x[ ioffx + i - 1 ] - s )
      + w[ ioffw + i - 1 ] ) + y[ ioffy + i - 1 ];
    x[ ioffx + i - 1 ] = s;
  }
}
//*************************************************************************
LaPack0.dlabad = function( smallReference, largeReference ) {
  if ( Math.log( largeReference.getValue() ) / Math.LN10 > 2000. ) {
    smallReference.setValue( Math.sqrt( smallReference.getValue() ) );
    largeReference.setValue( Math.sqrt( largeReference.getValue() ) );
  }
}
//*************************************************************************
LaPack0.dlacn2 = function( n, v, x, isgn, estReference, kaseReference,
isave, ioffv, ioffx, ioffisgn, ioffisave ) {
  var itmax = 5;
  if ( kaseReference.getValue() == 0 ) {
    for ( var i = 1; i <= n; i ++ ) {
      x[ ioffx + i - 1 ] = 1. / Number( n );
    }
    kaseReference.setValue( 1 );
    isave[ ioffisave ] = 1;
    return;
  }
  switch ( isave[ ioffisave ] ) {
    default:
    case 1: { // 20
      if ( n == 1 ) {
        v[ ioffv ] = x[ ioffx ];
        estReference.setValue( Math.abs( v[ ioffv ] ) );
        kaseReference.setValue( 0 );
        return;
      }
      estReference.setValue( Blas1.dasum( n, x, 1, ioffx ) );
      for ( i = 1; i <= n; i ++ ) {
        x[ ioffx + i - 1 ] = ( x[ ioffx + i - 1 ] >= 0. ? 1. : -1. );
        isgn[ ioffisgn + i - 1 ] = Math.round( x[ ioffx + i - 1 ] );
      }
      kaseReference.setValue( 2 );
      isave[ ioffisave ] = 2;
      return;
    }
    case 2: { // 40
      isave[ ioffisave + 1 ] = Blas1.idamax( n, x, 1, ioffx );
      isave[ ioffisave + 2 ] = 2;
      for ( i = 1; i <= n; i ++ ) x[ ioffx + i - 1 ] = 0.;
      x[ ioffx + isave[ ioffisave + 1 ] - 1 ] = 1.;
      kaseReference.setValue( 1 );
      isave[ ioffisave ] = 3;
      return;
    }
    case 3: { // 70
      Blas1.dcopy( n, x, 1, v, 1, ioffx, ioffv );
      var estold = estReference.getValue();
      estReference.setValue( Blas1.dasum( n, v, 1, ioffv ) );
      var goto90 = false;
      for ( i = 1; i <= n; i ++ ) {
        if ( Math.round( ( x[ ioffx + i - 1 ] >= 0. ? 1. : -1. ) )
        != isgn[ ioffisgn + i - 1 ] ) {
          goto90 = true;
          break;
        }
      }
      if ( ! goto90 || estReference.getValue() <= estold ) {
        var altsgn = 1.;
        for ( i = 1; i <= n; i ++ ) {
          x[ ioffx + i - 1 ] =
            altsgn * ( 1. + Number( i - 1 ) / Number( n - 1 ) );
          altsgn = - altsgn;
        }
        kaseReference.setValue( 1 );
        isave[ ioffisave ] = 5;
        return;
      }
      for ( i = 1; i <= n; i ++ ) {
        x[ ioffx + i - 1 ] = ( x[ ioffx + i - 1 ] >= 0. ? 1. : -1. );
        isgn[ ioffisgn + i - 1 ] = Math.round( x[ ioffx + i - 1 ] );
      }
      kaseReference.setValue( 2 );
      isave[ ioffisave ] = 4;
      return;
    }
    case 4: { // 110
      var jlast = isave[ ioffisave + 1 ];
      isave[ ioffisave + 1 ] = Blas1.idamax( n, x, 1, ioffx );
      if ( x[ ioffx + jlast - 1 ] !=
      Math.abs( x[ ioffx + isave[ ioffisave + 1 ] - 1 ] ) &&
      isave[ ioffisave + 2 ] < itmax ) {
        isave[ ioffisave + 2 ] += 1;
        for ( i = 1; i <= n; i ++ ) x[ ioffx + i - 1 ] = 0.;
        x[ ioffx + isave[ ioffisave + 1 ] - 1 ] = 1.;
        kaseReference.setValue( 1 );
        isave[ ioffisave ] = 3;
        return;
      }
      altsgn = 1.;
      for ( i = 1; i <= n; i ++ ) {
        x[ ioffx + i - 1 ] =
          altsgn * ( 1. + Number( i - 1 ) / Number( n - 1 ) );
        altsgn = - altsgn;
      }
      kaseReference.setValue( 1 );
      isave[ ioffisave ] = 5;
      return;
    }
    case 5: { // 140
      var temp =
        2. * ( Blas1.dasum( n, x, 1, ioffx ) / Number( 3 * n ) );
      if ( temp > estReference.getValue() ) {
        Blas1.dcopy( n, x, 1, v, 1, ioffx, ioffv );
        estReference.setValue( temp );
      }
      kaseReference.setValue( 0 );
    }
  }
}
//*************************************************************************
LaPack0.dlacon = function( n, v, x, isgn, estReference, kaseReference,
ioffv, ioffx, ioffisgn ) {
  if ( kaseReference.getValue() == 0 ) {
    for ( var i = 1; i <= n; i ++ ) {
      x[ ioffx + i - 1 ] = 1. / Number( n );
    }
    kaseReference.setValue( 1 );
    this.dlacon_jump = 1;
    return;
  }
  if ( this.dlacon_jump == 1 ) { // 20
    if ( n == 1 ) {
      v[ ioffv ] = x[ ioffx ];
      estReference.setValue( Math.abs( v[ ioffv ] ) );
      kaseReference.setValue( 0 ); // 150
    } else {
      estReference.setValue( Blas1.dasum( n, x, 1, ioffx ) );
      for ( i = 1; i <= n; i ++ ) {
        x[ ioffx + i - 1 ] = ( x[ ioffx + i - 1 ] >= 0. ? 1. : -1. );
        isgn[ ioffisgn + i - 1 ] = Math.round( x[ ioffx + i - 1 ] );
      }
      kaseReference.setValue( 2 );
      this.dlacon_jump = 2;
    }
    return;
  } else if ( this.dlacon_jump == 2 ) { // 40
    this.dlacon_j = Blas1.idamax( n, x, 1, ioffx );
    this.dlacon_iter = 2;
    for ( i = 1; i <= n; i ++ ) x[ ioffx + i - 1 ] = 0.; // 50
    x[ ioffx + this.dlacon_j - 1 ] = 1.;
    kaseReference.setValue( 1 );
    this.dlacon_jump = 3;
    return;
  } else if ( this.dlacon_jump == 3 ) { // 70
    Blas1.dcopy( n, x, 1, v, 1 , ioffx, ioffv );
    var estold = estReference.getValue();
    estReference.setValue( Blas1.dasum( n, v, 1, ioffv ) );
    var found = false;
    for ( i = 1; i <= n; i ++ ) {
      if ( Math.round( x[ ioffx + i - 1 ] >= 0. ? 1. : 1. ) !=
      isgn[ ioffisgn + i - 1 ] ) {
        found = true;
        break;
      }
    }
    if ( found ) { // 90
      if ( estReference.getValue() > estold ) {
        for ( i = 1; i <= n; i ++ ) {
          x[ ioffx + i - 1 ] = ( x[ ioffx + i - 1 ] >= 0. ? 1. : -1. );
          isgn[ ioffisgn + i - 1 ] = Math.round( x[ ioffx + i - 1 ] );
        }
        kaseReference.setValue( 2 );
        this.dlacon_jump = 4;
        return;
      }
    }
    var altsgn = 1.; // 120
    for ( i = 1; i <= n; i ++ ) {
      x[ ioffx + i - 1 ] =
        altsgn * ( 1. + Number( i - 1 ) / Number ( n - 1 ) );
      altsgn = - altsgn;
    }
    kaseReference.setValue( 1 );
    this.dlacon_jump = 5;
    return;
  } else if ( this.dlacon_jump == 4 ) { // 110
    var jlast = this.dlacon_j;
    this.dlacon_j = Blas1.idamax( n, x, 1, ioffx );
    if ( x[ ioffx + jlast - 1 ] !=
    Math.abs( x[ ioffx + this.dlacon_j - 1 ] ) &&
    this.dlacon_iter < this.dlacon_itmax ) {
      this.dlacon_iter ++;
      for ( i = 1; i <= n; i ++ ) x[ ioffx + i - 1 ] = 0.; // 50
      x[ ioffx + this.dlacon_j - 1 ] = 1.;
      kaseReference.setValue( 1 );
      this.dlacon_jump = 3;
      return;
    }
    altsgn = 1.; // 120
    for ( i = 1; i <= n; i ++ ) {
      x[ ioffx + i - 1 ] =
        altsgn * ( 1. + Number( i - 1 ) / Number ( n - 1 ) );
      altsgn = - altsgn;
    }
    kaseReference.setValue( 1 );
    this.dlacon_jump = 5;
    return;
  } else if ( this.dlacon_jump == 5 ) { // 140
    var temp =
      2. * ( Blas1.dasum( n, x, 1, ioffx ) / Number( 3 * n ) );
    if ( temp > estReference.getValue() ) {
      Blas1.dcopy( n, x, 1, v, 1, ioffx, ioffv );
      estReference.setValue( temp );
    }
    kaseReference.setValue( 0 ); // 150
  }
}
LaPack0.zlacon = function( n, v, x, isgn,
estReference, kaseReference ) {
  throw new Error("not programmed: complex vector");
}
//*************************************************************************
LaPack0.dlacpy = function( uplo, m, n, A, lda, B, ldb, ioffa,
ioffb ) {
  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    for ( var j = 1; j <= n; j ++ ) {
      for ( var i = 1; i <= Math.min( j, m ); i ++ ) {
        B[ ioffb + i - 1 + ( j - 1 ) * ldb ] =
          A[ ioffa + i - 1 + ( j - 1 ) * lda ];
      }
    }
  } else if ( uplo.charAt(0).toUpperCase() == 'L' ) {
    for ( j = 1; j <= n; j ++ ) {
      for ( i = j; i <= m; i ++ ) {
        B[ ioffb + i - 1 + ( j - 1 ) * ldb ] =
          A[ ioffa + i - 1 + ( j - 1 ) * lda ];
      }
    }
  } else {
    for ( j = 1; j <= n; j ++ ) {
      for ( i = 1; i <= m; i ++ ) {
        B[ ioffb + i - 1 + ( j - 1 ) * ldb ] =
          A[ ioffa + i - 1 + ( j - 1 ) * lda ];
      }
    }
  }
}
LaPack0.zlacpy = function( uplo, m, n, A, lda, B, ldb, ioffa,
ioffb ) {
  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    for ( var j = 1; j <= n; j ++ ) {
      for ( var i = 1; i <= Math.min( j, m ); i ++ ) {
        B[ ioffb + i - 1 + ( j - 1 ) * ldb ].setValue(
          A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
      }
    }
  } else if ( uplo.charAt(0).toUpperCase() == 'L' ) {
    for ( j = 1; j <= n; j ++ ) {
      for ( i = j; i <= m; i ++ ) {
        B[ ioffb + i - 1 + ( j - 1 ) * ldb ].setValue(
          A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
      }
    }
  } else {
    for ( j = 1; j <= n; j ++ ) {
      for ( i = 0; i <= m; i ++ ) {
        B[ ioffb + i - 1 + ( j - 1 ) * ldb ].setValue(
          A[ ioffa + i - 1 + ( j - 1 ) * lda ] );
      }
    }
  }
}
//*************************************************************************
LaPack0.dladiv = function( a, b, c, d, pReference, qReference ) {
  if ( Math.abs( d ) < Math.abs( c ) ) {
    var e = d / c;
    var f = c + d * e;
    pReference.setValue( ( a + b * e ) / f );
    qReference.setValue( ( b - a * e ) / f );
  } else {
    e = c / d;
    f = d + c * e;
    pReference.setValue( ( b + a * e ) / f );
    qReference.setValue( ( -a + b * e ) / f );
  }
}
LaPack0.zladiv = function( x, y ) {
  var zrReference = new NumberReference();
  var ziReference = new NumberReference();
  dladiv( x.real, x.imag, y.real, y.imag, zrReference, ziReference );
  var z = new Complex( zrReference.getValue(), ziReference.getValue() );
  return z;
}
//*************************************************************************
LaPack0.dlae2 = function( a, b, c, rt1Reference, rt2Reference ) {
  var sm = a + c;
  var df = a - c;
  var adf = Math.abs( df );
  var tb = b + b;
  var ab = Math.abs( tb );
  if ( Math.abs( a ) > Math.abs( c ) ) {
    var acmx = a;
    var acmn = c;
  } else {
    acmx = c;
    acmn = a;
  }
  if ( adf > ab ) {
    var rt = adf * Math.sqrt( 1. + Math.pow( ab / adf , 2 ) );
  } else if ( adf < ab ) {
    rt = ab * Math.sqrt( 1. + Math.pow( adf / ab , 2 ) );
  } else rt = ab * Math.sqrt( 2. );
  if ( sm < 0. ) {
    rt1Reference.setValue( 0.5 * ( sm - rt ) );
    rt2Reference.setValue( ( acmx / rt1Reference.getValue() ) * acmn
                  - ( b / rt1Reference.getValue() ) * b );
  } else if ( sm > 0. ) {
    rt1Reference.setValue( 0.5 * ( sm + rt ) );
    rt2Reference.setValue( ( acmx / rt1Reference.getValue() ) * acmn
                  - ( b / rt1Reference.getValue() ) * b );
  } else {
    rt1Reference.setValue( 0.5 * rt );
    rt2Reference.setValue( - 0.5 * rt );
  }
}
//*************************************************************************
LaPack0.dlaebz = function( ijob, nitmax, n, mmax, minp, nbmin,
abstol, reltol, pivmin, d, e, e2, nval, AB, c, mout, NAB, work, iwork,
info, ioffd, ioffe, ioffe2, ioffnval, ioffab, ioffc, ioffnab, ioffwork,
ioffiwork) {
  info.setValue( 0 );
  if ( ijob < 1 || ijob > 3 ) {
    info.setValue( -1 );
    return;
  }
  var j = -1;
  var ji = -1;
  var jp = -1;
  var tmp1 = Number.POSITIVE_INFINITY;
  if ( ijob == 1 ) {
    mout.setValue( 0 );
    for ( ji = 1; ji <= minp; ji ++ ) {
      for ( jp = 1; jp <= 2; jp ++ ) {
        tmp1 = d[ ioffd ] - AB[ ioffab + ji - 1 + ( jp - 1 ) * mmax ];
        if ( Math.abs( tmp1 ) < pivmin ) tmp1 = - pivmin;
        NAB[ ioffnab + ji - 1 + ( jp - 1 ) * mmax ] = 0;
        if ( tmp1 <= 0. ) {
          NAB[ ioffnab + ji - 1 + ( jp - 1 ) * mmax ] = 1;
        }
        for ( j = 2; j <= n; j ++ ) {
          tmp1 = d[ ioffd + j - 1 ] - e2[ ioffe2 + j - 2 ] / tmp1
               - AB[ ioffab + ji - 1 + ( jp - 1 ) * mmax ];
          if ( Math.abs( tmp1 ) < pivmin ) tmp1 = - pivmin;
          if ( tmp1 <= 0. ) {
            NAB[ ioffnab + ji - 1 + ( jp - 1 ) * mmax ] += 1;
          }
        }
      }
      mout.setValue( mout.getValue() +
        NAB[ ioffnab + ji - 1 + mmax ] - NAB[ ioffnab + ji - 1 ] );
    }
    return;
  }
  var kf = 1;
  var kl = minp;
  if ( ijob == 2 ) {
    for ( ji = 1; ji <= minp; ji ++ ) {
      c[ ioffc + ji - 1 ] = 0.5* ( AB[ ioffab + ji - 1 ] +
        AB[ ioffab + ji - 1 + mmax ] ); 
    }
  }
  for ( var jit = 1; jit <= nitmax; jit++ ) {
    var klnew = -1;
    if ( kl - kf + 1 >= nbmin && nbmin > 0 ) {
      for ( ji = kf; ji <= kl; ji ++ ) {
        work[ ioffwork + ji - 1 ] = d[ ioffd ] - c[ ioffc + ji - 1 ];
        iwork[ ioffiwork + ji - 1 ] = 0;
        if ( work[ ioffwork + ji - 1 ] <= pivmin ) {
          iwork[ ioffiwork + ji - 1 ] = 1;
          work[ ioffwork + ji - 1 ] =
            Math.min( work[ ioffwork + ji - 1 ], - pivmin );
        }
        for ( j = 2; j <= n; j ++ ) {
          work[ ioffwork + ji - 1 ] =
            d[ ioffd + j - 1 ] - e2[ ioffe2 + j - 2 ]
              / work[ ioffwork + ji - 1 ] - c[ ioffc + ji - 1 ];
          if ( work[ ioffwork + ji - 1 ] <= pivmin ) {
            iwork[ ioffiwork + ji - 1 ] += 1;
            work[ ioffwork + ji - 1 ] =
              Math.min( work[ ioffwork + ji - 1 ], - pivmin );
          }
        }
      }
      if ( ijob <= 2 ) {
        klnew = kl;
        for ( ji = kf; ji <= kl; ji ++ ) {
          iwork[ ioffiwork + ji - 1 ] =
            Math.min( NAB[ ioffnab + ji - 1 + mmax ],
            Math.max( NAB[ ioffnab + ji - 1 ],
            iwork[ ioffiwork + ji - 1 ] ) );
          if ( iwork[ ioffiwork + ji - 1 ] ==
          NAB[ ioffnab + ji - 1 + mmax ] ) {
            AB[ ioffab + ji - 1 + mmax ] = c[ ioffc + ji - 1 ];
          } else if ( iwork[ ioffiwork + ji - 1 ] ==
          NAB[ ioffnab + ji - 1 ] ) {
            AB[ ioffab + ji - 1 ] = c[ ioffc + ji - 1 ];
          } else {
            klnew ++;
            if ( klnew <= mmax ) {
              AB[ ioffab + klnew - 1 + mmax ] =
                AB[ ioffab + ji - 1 + mmax ];
              NAB[ ioffnab + klnew - 1 + mmax ] =
                NAB[ ioffnab + ji - 1 + mmax ];
              AB[ ioffab + klnew - 1 ] = c[ ioffc + ji - 1 ];
              NAB[ ioffnab + klnew - 1 ] = iwork[ ioffiwork + ji - 1 ];
              AB[ ioffab + ji - 1 + mmax ] = c[ ioffc + ji - 1 ];
              NAB[ ioffnab + ji - 1 + mmax ] =
                iwork[ ioffiwork + ji - 1 ];
            } else info.setValue( mmax + 1 );
          }
        }
        if ( info.getValue() != 0 ) return;
        kl = klnew;
      } else {
        for ( ji = kf; ji <= kl; ji ++ ) {
          if ( iwork[ ioffiwork + ji - 1 ] <=
          nval[ ioffnval + ji - 1 ] ) {
            AB[ ioffab + ji - 1 ] = c[ ioffc + ji - 1 ];
            NAB[ ioffnab + ji - 1 ] = iwork[ ioffiwork + ji - 1 ];
          }
          if ( iwork[ ioffiwork + ji - 1 ] >=
          nval[ ioffnval + ji - 1 ] ) {
            AB[ ioffab + ji - 1 + mmax ] = c[ ioffc + ji - 1 ];
            NAB[ ioffnab + ji - 1 + mmax ] =
              iwork[ ioffiwork + ji - 1 ];
          }
        }
      }
    } else {
      klnew = kl;
      for ( ji = kf; ji <= kl; ji ++ ) {
        tmp1 = c[ ioffc + ji - 1 ];
        var tmp2 = d[ ioffd ] - tmp1;
        var itmp1 = 0;
        if ( tmp2 <= pivmin ) {
          itmp1 = 1;
          tmp2 = Math.min( tmp2, -pivmin );
        }
        for ( j = 2; j <= n; j ++ ) {
          tmp2 = d[ ioffd + j - 1]
            - e2[ ioffe2 + j - 2 ] / tmp2 - tmp1;
          if ( tmp2 <= pivmin ) {
            itmp1 ++;
            tmp2 = Math.min( tmp2, - pivmin );
          }
        }
        if ( ijob <= 2 ) {
          itmp1 = Math.min( NAB[ ioffnab + ji - 1 + mmax ],
                  Math.max( NAB[ ioffnab + ji - 1 ] , itmp1 ) );
          if ( itmp1 == NAB[ ioffnab + ji - 1 + mmax ] ) {
            AB[ ioffab + ji - 1 + mmax ] = tmp1;
          } else if ( itmp1 == NAB[ ioffnab + ji - 1 ] ) {
            AB[ ioffab + ji - 1 ] = tmp1;
          } else if ( klnew < mmax ) {
            klnew ++;
            AB[ ioffab + klnew - 1 + mmax ] =
              AB[ ioffab + ji - 1 + mmax ];
            NAB[ ioffnab + klnew - 1 + mmax ] =
              NAB[ ioffnab + ji - 1 + mmax ];
            AB[ ioffab + klnew - 1 ] = tmp1;
            NAB[ ioffnab + klnew - 1 ] = itmp1;
            AB[ ioffab + ji - 1 + mmax ] = tmp1;
            NAB[ ioffnab + ji - 1 + mmax ] = itmp1;
          } else {
            info.setValue( mmax + 1 );
            return;
          }
        } else {
          if ( itmp1 <= nval[ ioffnval + ji - 1 ] ) {
            AB[ ioffab + ji - 1 ] = tmp1;
            NAB[ ioffnab + ji - 1 ] = itmp1;
          }
          if ( itmp1 >= nval[ ioffnval + ji - 1 ] ) {
            AB[ ioffab + ji - 1 + mmax ] = tmp1;
            NAB[ ioffnab + ji - 1 + mmax ] = itmp1;
          }
        }
      }
      kl = klnew;
    }
    var kfnew = kf;
    for ( ji = kf; ji <= kl; ji ++ ) {
      tmp1 = Math.abs( AB[ ioffab + ji - 1 + mmax ]
        - AB[ ioffab + ji - 1 ] );
      tmp2 = Math.max( Math.abs( AB[ ioffab + ji - 1 + mmax ] ),
                       Math.abs( AB[ ioffab + ji - 1 ] ) );
      if ( tmp1 < Math.max( Math.max( abstol, pivmin), reltol * tmp2 )
      || NAB[ ioffnab + ji - 1 ] >= NAB[ ioffnab + ji - 1 + mmax ] ) {
        if ( ji > kfnew ) {
          tmp1 = AB[ ioffab + ji - 1 ];
          tmp2 = AB[ ioffab + ji - 1 + mmax ];
          itmp1 = NAB[ ioffnab + ji - 1 ];
          var itmp2 = NAB[ ioffnab + ji - 1 + mmax ];
          AB[ ioffab + ji - 1 ] = AB[ ioffab + kfnew - 1 ];
          AB[ ioffab + ji - 1 + mmax ] =
            AB[ ioffab + kfnew - 1 + mmax ];
          NAB[ ioffnab + ji - 1 ] = NAB[ ioffnab + kfnew - 1 ];
          NAB[ ioffnab + ji - 1 + mmax ] =
            NAB[ ioffnab + kfnew - 1 + mmax ];
          AB[ ioffab + kfnew - 1 ] = tmp1;
          AB[ ioffab + kfnew - 1 + mmax ] = tmp2;
          NAB[ ioffnab + kfnew - 1 ] = itmp1;
          NAB[ ioffnab + kfnew - 1 + mmax ] = itmp2;
          if ( ijob == 3 ) {
            itmp1 = nval[ ioffnval + ji - 1 ];
            nval[ ioffnval + ji - 1 ] = nval[ ioffnval + kfnew - 1 ];
            nval[ ioffnval + kfnew - 1 ] = itmp1;
          }
        }
        kfnew ++;
      }
    }
    kf = kfnew;
    for ( ji = kf; ji <= kl; ji ++ ) {
      c[ ioffc + ji - 1 ] = 0.5 * ( AB[ ioffab + ji - 1 ]
        + AB[ ioffab + ji - 1 + mmax ] );
    }
    if ( kf > kl ) break;
  }
  info.setValue( Math.max( kl + 1 - kf, 0 ) );
  mout.setValue( kl );
}
//*************************************************************************
LaPack0.dlaed5 = function( i, d, z, delta, rho, dlamReference,
ioffd, ioffz, ioffdelta) {
  var b = Number.POSITIVE_INFINITY;
  var c = Number.POSITIVE_INFINITY;
  var del = d[ ioffd + 1 ] - d[ ioffd ];
  var tau = Number.POSITIVE_INFINITY;
  var temp = Number.POSITIVE_INFINITY;
  if ( i == 1 ) {
    var w = 1. + 2. * rho * ( z[ ioffz + 1 ] * z[ ioffz + 1 ]
      - z[ ioffz ] * z[ ioffz ] ) / del;
    if ( w > 0. ) {
      b = del + rho * ( z[ ioffz ] * z[ ioffz ]
        + z[ ioffz + 1 ] * z[ ioffz + 1 ] );
      c = rho * z[ ioffz ] * z[ ioffz ] * del;
      tau = 2. * c / ( b + Math.sqrt( Math.abs( b * b - 4. * c ) ) );
      dlamReference.setValue( d[ ioffd ] + tau );
      delta[ ioffdelta ] = - z[ ioffz ] / tau;
      delta[ ioffdelta + 1 ] = z[ ioffz + 1 ] / ( del - tau );
    } else {
      b = - del + rho * ( z[ ioffz ] * z[ ioffz ]
        + z[ ioffz + 1 ] * z[ ioffz + 1 ] );
      c = rho * z[ ioffz + 1 ] * z[ ioffz + 1 ] * del;
      if ( b > 0. ) {
        tau = - 2. * c
            / ( b + Math.sqrt( Math.abs( b * b + 4. * c ) ) );
      } else {
        tau = ( b - Math.sqrt( b * b + 4. * c ) ) / 2.;
      }
      dlamReference.setValue( d[ ioffd + 1 ] + tau );
      delta[ ioffdelta ] = - z[ ioffz ] / ( del + tau );
      delta[ ioffdelta + 1 ] = - z[ ioffz + 1 ] / tau;
    }
    temp = Math.sqrt( delta[ ioffdelta ] * delta[ ioffdelta ]
      + delta[ ioffdelta + 1 ] * delta[ ioffdelta + 1 ]);
    delta[ ioffdelta ] /= temp;
    delta[ ioffdelta + 1 ] /= temp;
  } else {
    b = - del + rho * ( z[ ioffz ] * z[ ioffz ]
      + z[ ioffz + 1 ] * z[ ioffz + 1 ] );
    c = rho * z[ ioffz + 1 ] * z[ ioffz + 1 ] * del;
    if ( b > 0. ) {
      tau = ( b + Math.sqrt( b * b + 4. * c ) ) / 2.;
    } else {
      tau = 2. * c / ( - b + Math.sqrt( Math.abs( b * b + 4. * c ) ) );
    }
    dlamReference.setValue( d[ ioffd + 1 ] + tau );
    delta[ ioffdelta ] = - z[ ioffz ] / ( del + tau );
    delta[ ioffdelta + 1 ] = - z[ ioffz + 1 ] / tau;
    temp = Math.sqrt( delta[ ioffdelta ] * delta[ ioffdelta ]
      + delta[ ioffdelta + 1 ] * delta[ ioffdelta + 1 ]);
    delta[ ioffdelta ] /= temp;
    delta[ ioffdelta + 1 ] /= temp;
  }
}
//*************************************************************************
LaPack0.dlaeda = function( n, tlvls, curlvl, curpbm, prmptr,
perm, givptr, Givcol, Givnum, Q, qptr, z, ztemp, info, ioffprmptr,
ioffperm, ioffgivptr, ioffgivcol, ioffgivnum, ioffq, ioffqptr, ioffz,
ioffztemp ) {
  info.setValue( 0 );
  if ( n < 0 ) info.setValue( -1 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlaeda', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  var mid = n / 2 + 1;
  var ptr = 1;
  var curr = ptr + curpbm * Math.pow( 2, curlvl )
    + Math.pow( 2, curlvl - 1 ) - 1;
  var bsiz1 =
    Math.round( 0.5 + Math.sqrt( Number( qptr[ ioffqptr + curr ]
    - qptr[ ioffqptr + curr - 1 ] ) ) );
  var bsiz2 =
    Math.round( 0.5 + Math.sqrt( Number( qptr[ ioffqptr + curr + 1 ]
    - qptr[ ioffqptr + curr ] ) ) );
  for ( var k = 1; k <= mid - bsiz1 - 1; k ++ ) {
    z[ ioffz + k - 1 ] = 0.;
  }
  Blas1.dcopy( bsiz1, Q, bsiz1, z, 1,
    ioffq + qptr[ ioffqptr + curr - 1 ] + bsiz1 - 2,
    ioffz + mid - bsiz1 - 1 );
  Blas1.dcopy( bsiz2, Q, bsiz2, z, 1,
    ioffq + qptr[ ioffqptr + curr ] - 1, ioffz + mid - 1 );
  for ( k = mid + bsiz2; k <= n; k ++ ) z[ ioffz + k - 1 ] = 0.;
  ptr = Math.pow( 2, tlvls ) + 1;
  for ( k = 1; k < curlvl; k ++ ) {
    curr = ptr + curpbm * Math.pow( 2, curlvl - k )
         + Math.pow( 2, curlvl - k - 1 ) - 1; 
    var psiz1 = prmptr[ ioffprmptr + curr ]
      - prmptr[ ioffprmptr + curr - 1 ];
    var psiz2 = prmptr[ ioffprmptr + curr + 1 ]
      - prmptr[ ioffprmptr + curr ];
    var zptr1 = mid - psiz1;
    for ( var i = givptr[ ioffgivptr + curr - 1 ];
    i <= givptr[ ioffgivptr + curr ] - 1; i ++ ) {
      Blas1.drot( 1, z, 1, z, 1, Givnum[ ioffgivnum + ( i - 1 ) * 2 ],
        Givnum[ ioffgivnum + 1 + ( i - 1 ) * 2 ], 
        ioffz + zptr1 + Givcol[ ioffgivcol + ( i - 1 ) * 2 ] - 2,
        ioffz + zptr1 + Givcol[ ioffgivcol + 1 + ( i - 1 ) * 2 ] - 2 );
    }
    for ( i = givptr[ ioffgivptr + curr ];
    i <= givptr[ ioffgivptr + curr + 1 ] - 1; i ++ ) {
      Blas1.drot( 1, z, 1, z, 1, Givnum[ ioffgivnum + ( i - 1 ) * 2 ],
        Givnum[ ioffgivnum + 1 + ( i - 1 ) * 2 ], 
        ioffz + mid + Givcol[ ioffgivcol + ( i - 1 ) * 2 ] - 2,
        ioffz + mid + Givcol[ ioffgivcol + 1 + ( i - 1 ) * 2 ] - 2 );
    }
    psiz1 = prmptr[ ioffprmptr + curr ]
      - prmptr[ ioffprmptr + curr - 1 ];
    psiz2 = prmptr[ ioffprmptr + curr + 1 ]
      - prmptr[ ioffprmptr + curr ];
    for ( i = 0; i <= psiz1 - 1; i ++ ) {
      ztemp[ ioffztemp + i ] =
        z[ ioffz + zptr1 + perm[ ioffperm
          + prmptr[ ioffprmptr + curr - 1 ] + i - 1 ] - 2 ];
    }
    for ( i = 0; i <= psiz2 - 1; i ++ ) {
      ztemp[ ioffztemp + psiz1 + i ] =
        z[ ioffz + mid + perm[ ioffperm
          + prmptr[ ioffprmptr + curr ] + i - 1 ] - 2 ];
    }
    bsiz1 = Math.round( 0.5
      + Math.sqrt( Number( qptr[ ioffqptr + curr ]
      - qptr[ ioffqptr + curr - 1 ] ) ) );
    bsiz2 = Math.round( 0.5
      + Math.sqrt( Number( qptr[ ioffqptr + curr + 1 ]
      - qptr[ ioffqptr + curr ] ) ) );
    if ( bsiz1 > 0 ) {
      Blas2.dgemv( 'T', bsiz1, bsiz1, 1., Q, bsiz1, ztemp, 1, 0., z,
        1, ioffq + qptr[ ioffqptr + curr - 1 ] - 1, ioffztemp,
        ioffz + zptr1 - 1 );
    }
    Blas1.dcopy( psiz1 - bsiz1, ztemp, 1, z, 1, ioffztemp + bsiz1,
      ioffz + zptr1 + bsiz1 - 1 );
    if ( bsiz2 > 0 ) {
      Blas2.dgemv( 'T', bsiz2, bsiz2, 1., Q, bsiz2, ztemp, 1, 0., z,
        1, ioffq + qptr[ ioffqptr + curr ] - 1, ioffztemp + psiz1,
        ioffz + mid - 1 );
    }
    Blas1.dcopy( psiz2 - bsiz2, ztemp, 1, z, 1,
      ioffztemp + psiz1 + bsiz2, ioffz + mid + bsiz2 - 1 );
    ptr += Math.pow( 2, tlvls - k );
  }
}
//*************************************************************************
LaPack0.dlaev2 = function( a, b, c, rt1Reference, rt2Reference,
cs1Reference, sn1Reference ) {
  var sm = a + c;
  var df = a - c;
  var adf = Math.abs( df );
  var tb = b + b;
  var ab = Math.abs( tb );
  if ( Math.abs( a ) > Math.abs( c ) ) {
    var acmx = a;
    var acmn = c;
  } else {
    acmx = c;
    acmn = a;
  }
  var rt = Number.POSITIVE_INFINITY;
  if ( adf > ab ) {
    rt = adf * Math.sqrt( 1. + Math.pow( ab / adf , 2 ) );
  } else if ( adf < ab ) {
    rt = ab * Math.sqrt( 1. + Math.pow( adf / ab , 2 ) );
  } else rt = ab * Math.sqrt( 2. );
  var sgn1 = 0;
  if ( sm < 0. ) {
    rt1Reference.setValue( 0.5 * ( sm - rt ) );
    sgn1 = -1;
    rt2Reference.setValue( ( acmx / rt1Reference.getValue() ) * acmn
                  - ( b / rt1Reference.getValue() ) * b );
  } else if ( sm > 0. ) {
    rt1Reference.setValue( 0.5 * ( sm + rt ) );
    sgn1 = 1;
    rt2Reference.setValue( ( acmx / rt1Reference.getValue() ) * acmn
                  - ( b / rt1Reference.getValue() ) * b );
  } else {
    rt1Reference.setValue( 0.5 * rt );
    rt2Reference.setValue( - 0.5 * rt );
    sgn1 = 1;
  }
  var sgn2 = 0;
  if ( df >= 0. ) {
    var cs = df + rt;
    sgn2 = 1;
  } else {
    cs = df - rt;
    sgn2 = -1;
  }
  var acs = Math.abs( cs );
  var ct = Number.POSITIVE_INFINITY;
  var tn = Number.POSITIVE_INFINITY;
  if ( acs > ab ) {
    ct = - tb / cs;
    sn1Reference.setValue( 1. / Math.sqrt( 1. + ct * ct ) );
    cs1Reference.setValue( ct * sn1Reference.getValue() );
  } else {
    if ( ab == 0. ) {
      cs1Reference.setValue( 1. );
      sn1Reference.setValue( 0. );
    } else {
      tn = - cs / tb;
      cs1Reference.setValue( 1. / Math.sqrt( 1. + tn * tn ) );
      sn1Reference.setValue( tn * cs1Reference.getValue() );
    }
  }
  if ( sgn1 == sgn2 ) {
    tn = cs1Reference.getValue();
    cs1Reference.setValue( - sn1Reference.getValue() );
    sn1Reference.setValue( tn );
  }
}
LaPack0.zlaev2 = function( a, b, c,
rt1Reference, rt2Reference, cs1Reference, sn1Reference ) {
  throw new Error("not programmed: complex numbers");
}
//*************************************************************************
LaPack0.dlag2 = function( A, lda, B, ldb, safmin, scale1Reference,
scale2Reference, wr1Reference, wr2Reference, wiReference, ioffa, ioffb ) {
  var fuzzy1 = 1. + 1.e-5;
  var rtmin = Math.sqrt( safmin );
  var rtmax = 1. / rtmin;
  var safmax = 1. / safmin;
  var anorm = Math.max(
    Math.max( Math.abs( A[ ioffa ] ) + Math.abs( A[ ioffa + 1 ] ),
      Math.abs( A[ ioffa + lda ] )
      + Math.abs( A[ ioffa + 1 + lda ] ) ), safmin );
  var ascale = 1. / anorm;
  var a11 = ascale * A[ ioffa ];
  var a21 = ascale * A[ ioffa + 1 ];
  var a12 = ascale * A[ ioffa + lda ];
  var a22 = ascale * A[ ioffa + 1 + lda ];
  var b11 = B[ ioffb ];
  var b12 = B[ ioffb + ldb ];
  var b22 = B[ ioffb + 1 + ldb ];
  var bmin = rtmin * Math.max(
    Math.max( Math.abs( b11 ) , Math.abs( b12 ) ),
    Math.max( Math.abs( b22 ) , rtmin ) );
  if ( Math.abs( b11 ) < bmin ) b11 = ( b11 >= 0. ? bmin : - bmin );
  if ( Math.abs( b22 ) < bmin ) b22 = ( b22 >= 0. ? bmin : - bmin );
  var bnorm = Math.max(
    Math.max( Math.abs( b11 ) , Math.abs( b12 ) + Math.abs( b22 ) ),
    safmin );
  var bsize = Math.max( Math.abs( b11 ), Math.abs( b22 ) );
  var bscale = 1. / bsize;
  b11 *= bscale;
  b12 *= bscale;
  b22 *= bscale;
  var binv11 = 1. / b11;
  var binv22 = 1. / b22;
  var s1 = a11 * binv11;
  var s2 = a22 * binv22;
  if ( Math.abs( s1 ) <= Math.abs( s2 ) ) {
    var as12 = a12 - s1 * b12;
    var as22 = a22 - s1 * b22;
    var ss = a21 * ( binv11 * binv22 );
    var abi22 = as22 * binv22 - ss * b12;
    var pp = 0.5 * abi22;
    var shift = s1;
  } else {
    as12 = a12 - s2 * b12;
    var as11 = a11 - s2 * b11;
    ss = a21 * ( binv11 * binv22 );
    abi22 = - ss * b12;
    pp = 0.5 * ( as11 * binv11 + abi22 );
    shift = s2;
  }
  var qq = ss * as12;
  if ( Math.abs( pp * rtmin ) >= 1. ) {
    var discr = Math.pow( rtmin * pp , 2 ) + qq * safmin;
    var r = Math.sqrt( Math.abs( discr ) ) * rtmax;
  } else {
    if ( pp * pp + Math.abs( qq ) <= safmin ) {
      discr = Math.pow( rtmax * pp, 2 ) + qq * safmax;
      r = Math.sqrt( Math.abs( discr ) ) * rtmin;
    } else {
      discr = pp * pp + qq;
      r = Math.sqrt( Math.abs( discr ) );
    }
  }
  if ( discr >= 0. || r == 0. ) {
    var sum =
      pp + ( pp >= 0. ? Math.abs( r ) : - Math.abs( r ) );
    var diff =
      pp - ( pp >= 0. ? Math.abs( r ) : - Math.abs( r ) );
    var wbig = shift + sum;
    var wsmall = shift + diff;
    if ( 0.5 * Math.abs( wbig ) >
    Math.max( Math.abs( wsmall ), safmin ) ) {
      var wdet =
        ( a11 * a22 - a12 * a21 ) * ( binv11 * binv22 );
      wsmall = wdet / wbig;
    }
    if ( pp > abi22 ) {
      wr1Reference.setValue( Math.min( wbig, wsmall ) );
      wr2Reference.setValue( Math.max( wbig, wsmall ) );
    } else {
      wr1Reference.setValue( Math.max( wbig, wsmall ) );
      wr2Reference.setValue( Math.min( wbig, wsmall ) );
    }
    wiReference.setValue( 0. );
  } else {
    wr1Reference.setValue( shift + pp );
    wr2Reference.setValue( wr1Reference.getValue() );
    wiReference.setValue( r );
  }
  var c1 = bsize * ( safmin * Math.max( 1., ascale ) );
  var c2 = safmin * Math.max( 1., bnorm );
  var c3 = bsize * safmin;
  if ( ascale <= 1. && bsize <= 1. ) {
    var c4 = Math.min( 1., ( ascale / safmin ) * bsize );
  } else  c4 = 1.;
  if ( ascale <= 1. || bsize <= 1. ) {
    var c5 = Math.min( 1., ascale * bsize );
  } else  c5 = 1.;
  var wabs = Math.abs( wr1Reference.getValue() )
           + Math.abs( wiReference.getValue() );
  var wsize = Math.max( Math.max( safmin , c1 ),
    Math.max( fuzzy1 * ( wabs * c2 + c3 ),
              Math.min( c4, 0.5 * Math.max( wabs, c5 ) ) ) );
  if ( wsize != 1. ) {
    var wscale = 1. / wsize;
    if ( wsize > 1. ) {
      scale1Reference.setValue( ( Math.max( ascale, bsize ) * wscale )
                   * Math.min( ascale, bsize ) );
    } else {
      scale1Reference.setValue( ( Math.min( ascale, bsize ) * wscale )
                   * Math.max( ascale, bsize ) );
    }
    wr1Reference.setValue( wr1Reference.getValue() * wscale );
    if ( wiReference.getValue() != 0. ) {
      wiReference.setValue( wiReference.getValue() * wscale );
      wr2Reference.setValue( wr1Reference.getValue() );
      scale2Reference.setValue( scale1Reference.getValue() );
    }
  } else {
    scale1Reference.setValue( ascale * bsize );
    scale2Reference.setValue( scale1Reference.getValue() );
  }
  if ( wiReference.getValue() == 0. ) {
    wsize = Math.max( Math.max( safmin , c1 ),
      Math.max( fuzzy1 * ( Math.abs( wr2Reference.getValue() ) * c2 + c3 ),
      Math.min( c4, 0.5 * Math.max( Math.abs( wr2Reference.getValue() ),
                                    c5 ) ) ) );
    if ( wsize != 1. ) {
      wscale = 1. / wsize;
      if ( wsize > 1. ) {
        scale2Reference.setValue( ( Math.max( ascale, bsize ) * wscale )
               * Math.min( ascale, bsize ) );
      } else {
        scale2Reference.setValue( ( Math.min( ascale, bsize ) * wscale )
               * Math.max( ascale, bsize ) );
      }
      wr2Reference.setValue( wr2Reference.getValue() * wscale );
    } else {
      scale2Reference.setValue( ascale * bsize );
    }
  }
}
//*************************************************************************
LaPack0.dlagtm = function( trans, n, nrhs, alpha, dl, d, du, X,
ldx, beta, B, ldb, ioffdl, ioffd, ioffdu, ioffx, ioffb ) {
  if ( n == 0 ) return;
  var i = -1;
  var j = -1;
  if ( beta == 0. ) {
    for ( j = 1; j <= nrhs; j ++ ) {
      for ( i = 1; i <= n; i ++ ) {
        B[ ioffb + i - 1 + ( j - 1 ) * ldb ] = 0.;
      }
    }
  } else if ( beta == -1. ) {
    for ( j = 1; j <= nrhs; j ++ ) {
      for ( i = 1; i <= n; i ++ ) {
        B[ ioffb + i - 1 + ( j - 1 ) * ldb ] = 
          -B[ ioffb + i - 1 + ( j - 1 ) * ldb ];
      }
    }
  }
  if ( alpha == 1. ) {
    if ( trans.charAt(0).toUpperCase() == 'N' ) {
      for ( j = 1; j <= nrhs; j ++ ) {
        if ( n == 1 ) {
          B[ ioffb + ( j - 1 ) * ldb ] +=
            d[ ioffd ] * X[ ioffx + ( j - 1 ) * ldx ];
        } else {
          B[ ioffb + ( j - 1 ) * ldb ] +=
            d[ ioffd ] * X[ ioffx + ( j - 1 ) * ldx ]
            + du[ ioffdu ] * X[ ioffx + 1 + ( j - 1 ) * ldx ];
          B[ ioffb + n - 1 + ( j - 1 ) * ldb ] +=
            dl[ ioffdl + n - 2 ] * X[ ioffx + n - 2 + ( j - 1 ) * ldx ]
            + d[ ioffd + n - 1 ]
            * X[ ioffx + n - 1 + ( j - 1 ) * ldx ];
          for ( i = 2; i <= n - 1; i ++ ) {
            B[ ioffb + i - 1 + ( j - 1 ) * ldb ] +=
              dl[ ioffdl + i - 2 ]
              * X[ ioffx + i - 2 + ( j - 1 ) * ldx ]
              + d[ ioffd + i - 1 ]
              * X[ ioffx + i - 1 + ( j - 1 ) * ldx ]
              + du[ ioffdu + i - 1 ]
              * X[ ioffx + i + ( j - 1 ) * ldx ];
          }
        }
      }
    } else {
      for ( j = 1; j <= nrhs; j ++ ) {
        if ( n == 1 ) {
          B[ ioffb + ( j - 1 ) * ldb ]
            += d[ ioffd ] * X[ ioffx + ( j - 1 ) * ldx ];
        } else {
          B[ ioffb + ( j - 1 ) * ldb ] +=
            d[ ioffd ] * X[ ioffx + ( j - 1 ) * ldx ]
            + dl[ ioffdl ] * X[ ioffx + 1 + ( j - 1 ) * ldx ];
          B[ ioffb + n - 1 + ( j - 1 ) * ldb ] +=
            du[ ioffdu + n - 2 ] * X[ ioffx + n - 2 + ( j - 1 ) * ldx ]
            + d[ ioffd + n - 1 ]
            * X[ ioffx + n - 1 + ( j - 1 ) * ldx ];
          for ( i = 2; i <= n - 1; i ++ ) { 
            B[ ioffb + i - 1 + ( j - 1 ) * ldb ] +=
              du[ ioffdu + i - 2 ]
              * X[ ioffx + i - 2 + ( j - 1 ) * ldx ]
              + d[ ioffd + i - 1 ]
              * X[ ioffx + i - 1 + ( j - 1 ) * ldx ]
              + dl[ ioffdl + i - 1 ]
              * X[ ioffx + i + ( j - 1 ) * ldx ];
          }
        }
      }
    }
  } else if ( alpha == -1. ) {
    if ( trans.charAt(0).toUpperCase() == 'N' ) {
      for ( j = 1; j <= nrhs; j ++ ) {
        if ( n == 1 ) {
          B[ ioffb + ( j - 1 ) * ldb ] -=
            d[ ioffd ] * X[ ioffx + ( j - 1 ) * ldx ];
        } else {
          B[ ioffb + ( j - 1 ) * ldb ] -=
            d[ ioffd ] * X[ ioffx + ( j - 1 ) * ldx ]
            + du[ ioffdu ] * X[ ioffx + 1 + ( j - 1 ) * ldx ];
          B[ ioffb + n - 1 + ( j - 1 ) * ldb ] -=
            dl[ ioffdl + n - 2 ] * X[ ioffx + n - 2 + ( j - 1 ) * ldx ]
            + d[ ioffd + n - 1 ]
            * X[ ioffx + n - 1 + ( j - 1 ) * ldx ];
          for ( i = 2; i <= n - 1; i ++ ) { 
            B[ ioffb + i - 1 + ( j - 1 ) * ldb ] -=
              dl[ ioffdl + i - 2 ]
              * X[ ioffx + i - 2 + ( j - 1 ) * ldx ]
              + d[ ioffd + i - 1 ]
              * X[ ioffx + i - 1 + ( j - 1 ) * ldx ]
              + du[ ioffdu + i - 1 ]
              * X[ ioffx + i + ( j - 1 ) * ldx ];
          }
        }
      }
    } else {
      for ( j = 1; j <= nrhs; j ++ ) {
        if ( n == 1 ) {
          B[ ioffb + ( j - 1 ) * ldb ] -=
            d[ ioffd ] * X[ ioffx + ( j - 1 ) * ldx ];
        } else {
          B[ ioffb + ( j - 1 ) * ldb ] -=
            d[ ioffd ] * X[ ioffx + ( j - 1 ) * ldx ]
            + dl[ ioffdl ] * X[ ioffx + 1 + ( j - 1 ) * ldx ];
          B[ ioffb + n - 1 + ( j - 1 ) * ldb ] -=
            du[ ioffdu + n - 2 ] * X[ ioffx + n - 2 + ( j - 1 ) * ldx ]
            + d[ ioffd + n - 1 ]
            * X[ ioffx + n - 1 + ( j - 1 ) * ldx ];
          for ( i = 2; i <= n - 1; i ++ ) { 
            B[ ioffb + i - 1 + ( j - 1 ) * ldb ] -=
              du[ ioffdu + i - 2 ]
              * X[ ioffx + i - 2 + ( j - 1 ) * ldx ]
              + d[ ioffd + i - 1 ]
              * X[ ioffx + i - 1 + ( j - 1 ) * ldx ]
              + dl[ ioffdl + i - 1 ]
              * X[ ioffx + i + ( j - 1 ) * ldx ];
          }
        }
      }
    }
  }
}
LaPack0.zlagtm = function( trans, n, nrhs, alpha, dl, d, du, X,
ldx, beta, B, ldb) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dlamch = function( cmach ) {
  if ( this.dlamch_first == undefined ) {
    this.dlamch_first = false;
    var beta = new IntReference( -1 );
    var it = new IntReference( -1 );
    var lrndReference = new BooleanReference( false );
    var imin = new IntReference( -1 );
    var imax = new IntReference( -1 );
    this.dlamch_epsReference = new NumberReference();
    this.dlamch_rminReference = new NumberReference();
    this.dlamch_rmaxReference = new NumberReference();
    this.dlamc2( beta, it, lrndReference, this.dlamch_epsReference, imin,
      this.dlamch_rminReference, imax, this.dlamch_rmaxReference );
    this.dlamch_base = beta.getValue();
    this.dlamch_t = it.getValue();
    if ( lrndReference.getValue() ) {
      this.dlamch_rnd = 1.;
      this.dlamch_epsReference.setValue(
        Math.pow( this.dlamch_base , 1 - it.getValue() ) / 2.);
    } else {
      this.dlamch_rnd = 0.;
      this.dlamch_epsReference.setValue(
        Math.pow( this.dlamch_base , 1 - it.getValue() ) );
    }
    this.dlamch_prec =
      this.dlamch_epsReference.getValue() * this.dlamch_base;
    this.dlamch_emin = imin.getValue();
    this.dlamch_emax = imax.getValue();
    this.dlamch_sfmin = this.dlamch_rminReference.getValue();
    var small = 1. / this.dlamch_rmaxReference.getValue();
    if ( small >= this.dlamch_sfmin ) {
      this.dlamch_sfmin =
        small * ( 1. + this.dlamch_epsReference.getValue() );
    }
  }
  switch ( cmach.charAt(0).toUpperCase() ) {
    case "E": return this.dlamch_epsReference.getValue();
    case "S": return this.dlamch_sfmin;
    case "B": return this.dlamch_base;
    case "P": return this.dlamch_prec;
    case "N": return this.dlamch_t;
    case "R": return this.dlamch_rnd;
    case "M": return this.dlamch_emin;
    case "U": return this.dlamch_rminReference.getValue();
    case "L": return this.dlamch_emax;
    case "O": return this.dlamch_rmaxReference.getValue();
    default: return Number.POSITIVE_INFINITY;
  }
}
LaPack0.dlamc1 = function( betaReference, tReference, rndReference,
ieee1Reference ) {
  if ( this.dlamc1_first == undefined) {
    this.dlamc1_first = false;
    var one = 1;
    var a = 1;
    var c = 1;
    while ( c == one ) {
      a *= 2;
      c = this.dlamc3( a, one );
      c = this.dlamc3( c, -a );
    }
    var b = 1;
    c = this.dlamc3( a, b );
    while ( c == a ) {
      b *= 2;
      c = this.dlamc3( a, b );
    }
    var qtr = one / 4;
    var savec = c;
    c = this.dlamc3( c, -a );
    this.dlamc1_lbeta = Math.floor( c + qtr );
    b = this.dlamc1_lbeta;
    var f = this.dlamc3( b / 2, -b / 100 );
    c = this.dlamc3( f, a );
    this.dlamc1_lrnd = ( c == a );
    f = this.dlamc3( b / 2, b / 100 );
    c = this.dlamc3( f , a );
    if ( this.dlamc1_lrnd && ( c == a ) ) this.dlamc1_lrnd = false;

    var t1 = this.dlamc3( b / 2, a );
    var t2 = this.dlamc3( b / 2, savec );
    this.dlamc1_lieee1 = ( t1 == a ) && ( t2 > savec ) && this.dlamc1_lrnd;

    this.dlamc1_lt = 0;
    a = 1;
    c = 1;
    while ( c == one ) {
      this.dlamc1_lt ++;
      a *= this.dlamc1_lbeta;
      c = this.dlamc3( a, one );
      c = this.dlamc3( c, -a );
    }
  }
  betaReference.setValue( this.dlamc1_lbeta );
  tReference.setValue( this.dlamc1_lt );
  rndReference.setValue( this.dlamc1_lrnd );
  ieee1Reference.setValue( this.dlamc1_lieee1 );
}
LaPack0.dlamc2 = function( betaReference, tReference, rndReference,
epsReference, eminReference, rminReference, emaxReference, rmaxReference )
{
  if ( this.dlamc2_first == undefined ) {
    this.dlamc2_first = false;
    var zero = 0;
    var one = 1;
    var two = 2;
    var lieee1Reference = new BooleanReference();
    this.dlamc2_lbeta = new IntReference();
    this.dlamc2_lt = new IntReference();
    this.dlamc2_lrndReference = new BooleanReference();
    this.dlamc1( this.dlamc2_lbeta, this.dlamc2_lt,
      this.dlamc2_lrndReference, lieee1Reference );
    var b = this.dlamc2_lbeta.getValue();
    var a = Math.pow( b, -this.dlamc2_lt.getValue() );
    this.dlamc2_leps = a;

    b = two / 3;
    var half = one / 2;
    var sixth = this.dlamc3( b, -half );
    var third = this.dlamc3( sixth, sixth );
    b = this.dlamc3( third, - half );
    b = this.dlamc3( b, sixth );
    b = Math.abs( b );
    if ( b < this.dlamc2_leps ) b = this.dlamc2_leps;

    this.dlamc2_leps = 1;
    while ( this.dlamc2_leps > b && b > zero ) {
      this.dlamc2_leps = b;
      var c = this.dlamc3( half * this.dlamc2_leps,
        Math.pow( two, 5 ) * Math.pow( this.dlamc2_leps , 2 ) );
      c = this.dlamc3( half, -c );
      b = this.dlamc3( half , c )
      c = this.dlamc3( half, -b );
      b = this.dlamc3( half , c )
    }
    if ( a < this.dlamc2_leps ) this.dlamc2_leps = a;

    var rbase = one / this.dlamc2_lbeta.getValue();
    var small = one;
    for ( var i = 1; i <= 3; i ++ ) {
      small = this.dlamc3( small * rbase , zero );
    }
    a = this.dlamc3( one, small );
    var ngpmin = new IntReference(-1);
    this.dlamc4( ngpmin, one, this.dlamc2_lbeta.getValue() );
    var ngnmin = new IntReference( -1 );
    this.dlamc4( ngnmin, -one, this.dlamc2_lbeta.getValue() );
    var gpmin = new IntReference( -1 );
    this.dlamc4( gpmin, a, this.dlamc2_lbeta.getValue() );
    var gnmin = new IntReference( -1 );
    this.dlamc4( gnmin, -a, this.dlamc2_lbeta.getValue() );

    var ieee = false;
    if ( ngpmin.getValue() == ngnmin.getValue()
    && gpmin.getValue() == gnmin.getValue() ) {
      if ( ngpmin.getValue() == gpmin.getValue() ) {
        this.dlamc2_lemin = ngpmin.getValue();
      } else if ( gpmin.getValue() - ngpmin.getValue() == 3 ) {
        this.dlamc2_lemin
          = ngpmin.getValue() - 1 + this.dlamc2_lt.getValue();
        ieee = true;
      } else {
        this.dlamc2_lemin
          = Math.min( ngpmin.getValue() , gpmin.getValue() );
        this.dlamc2_iwarn = true;
      }
    } else if ( ngpmin.getValue() == gpmin.getValue() &&
    ngnmin.getValue() == gnmin.getValue() ) {
      if ( Math.abs( ngpmin.getValue() - ngnmin.getValue() ) == 1 ) {
        this.dlamc2_lemin =
          Math.max( ngpmin.getValue() , ngnmin.getValue() );
      } else {
        this.dlamc2_lemin =
          Math.min( ngpmin.getValue() , ngnmin.getValue() );
        this.dlamc2_iwarn = true;
      }
    } else if ( Math.abs( ngpmin.getValue() - ngnmin.getValue() ) == 1 &&
    gpmin.getValue() == gnmin.getValue() ) {
      if ( gpmin.getValue()
      - Math.min( ngpmin.getValue(), ngnmin.getValue() ) == 3) {
        this.dlamc2_lemin =
          Math.max( ngpmin.getValue(), ngnmin.getValue() )
                   - 1 + this.dlamc2_lt.getValue();
      } else {
        this.dlamc2_lemin =
          Math.min( ngpmin.getValue() , ngnmin.getValue() );
        this.dlamc2_iwarn = true;
      }
    } else {
      this.dlamc2_lemin =
        Math.min( Math.min( ngpmin.getValue(), ngnmin.getValue() ),
        Math.min( gpmin.getValue(), gnmin.getValue() ) );
      this.dlamc2_iwarn = true;
    }
    if ( this.dlamc2_iwarn ) {
      this.dlamc2_first = true;
      document.getElementById("debug_textarea").value +=
        "WARNING. The value EMIN may be incorrect:- EMIN = "
        + this.dlamc2_lemin + "\n";
    }
    ieee = ieee || lieee1Reference.getValue();
    this.dlamc2_lrmin = 1;
    for ( i = 1; i <= 1 - this.dlamc2_lemin; i ++ ) {
      this.dlamc2_lrmin = this.dlamc3( this.dlamc2_lrmin * rbase, zero );
    }
    this.dlamc2_lrmaxReference = new NumberReference();
    this.dlamc2_lemaxReference = new IntReference();
    this.dlamc5( this.dlamc2_lbeta.getValue(), this.dlamc2_lt.getValue(),
      this.dlamc2_lemin, ieee, this.dlamc2_lemaxReference,
      this.dlamc2_lrmaxReference );
  }
  betaReference.setValue( this.dlamc2_lbeta.getValue() );
  tReference.setValue( this.dlamc2_lt.getValue() );
  rndReference.setValue( this.dlamc2_lrndReference.getValue() );
  epsReference.setValue( this.dlamc2_leps );
  eminReference.setValue( this.dlamc2_lemin );
  rminReference.setValue( this.dlamc2_lrmin );
  emaxReference.setValue( this.dlamc2_lemaxReference.getValue() );
  rmaxReference.setValue( this.dlamc2_lrmaxReference.getValue() );
}
LaPack0.dlamc3 = function( a, b ) {
  return a + b;
}
LaPack0.dlamc4 = function( eminReference, start, base ) {
  var a = start;
  var one = 1;
  var rbase = one / base;
  var zero = 0;
  eminReference.setValue( 1 );
  var b1 = this.dlamc3( a * rbase, zero );
  var c1 = a;
  var c2 = a;
  var d1 = a;
  var d2 = a;
  while ( c1 == a && c2 == a && d1 == a && d2 == a ) {
    eminReference.setValue( eminReference.getValue() - 1 );
    a = b1;
    b1 = this.dlamc3( a / base , zero );
    c1 = this.dlamc3( b1 * base , zero );
    d1 = zero;
    for ( var i = 1; i <= base; i ++ ) d1 += b1;
    var b2 = this.dlamc3( a * rbase, zero );
    c2 = this.dlamc3( b2 / rbase , zero );
    d2 = zero;
    for ( i = 1; i <= base; i ++ ) d2 += b2;
  }
}
LaPack0.dlamc5 = function( beta, p, emin, ieee, emax, rmaxReference ) {
  var zero = 0.;
  var one = 1.;
  var lexp = 1;
  var exbits = 1;
  var guess = lexp * 2;
  while ( guess <= - emin ) {
    lexp = guess;
    exbits ++;
    guess = lexp * 2;
  }
  var uexp = lexp;
  if ( lexp != -emin ) {
    uexp = guess;
    exbits ++;
  }
  var expsum =
    ( ( uexp + emin ) > ( -lexp - emin ) ? 2 * lexp : 2 * uexp ); 
  emax.setValue( expsum + emin - 1 );
  var nbits = 1 + exbits + p;
  if ( nbits % 2 == 1 && beta == 2 ) emax.setValue( emax.getValue() - 1 );
  if ( ieee ) emax.setValue( emax.getValue() - 1 );
  var recbas = one / beta;
  var z = beta - one;
  var y = zero;
  var oldy = Number.POSITIVE_INFINITY;
  for ( var i = 1; i <= p; i ++ ) {
    z *= recbas;
    if ( y < one ) oldy = y;
    y = this.dlamc3( y , z );
  }
  if ( y >= one ) y = oldy;
  for ( i = 1; i <= emax.getValue(); i ++ ) {
    y = this.dlamc3( y * beta , zero );
  }
  rmaxReference.setValue( y );
}
//*************************************************************************
LaPack0.dlamrg = function( n1, n2, a, dtrd1, dtrd2, index,
ioffa, ioffindex ) {
  var n1sv = n1;
  var n2sv = n2;
  var ind1 = ( dtrd1 > 0 ? 1 : n1 );
  var ind2 = ( dtrd2 > 0 ? 1 + n1 : n1 + n2 );
  var i = 1;
  while ( n1sv > 0 && n2sv > 0 ) {
    if ( a[ ioffa + ind1 - 1 ] <= a[ ioffa + ind2 - 1 ] ) {
      index[ ioffindex + i - 1 ] = ind1;
      i ++;
      ind1 += dtrd1;
      n1sv --;
    } else {
      index[ ioffindex + i - 1 ] = ind2;
      i ++;
      ind2 += dtrd2;
      n2sv --;
    }
  }
  if ( n1sv == 0 ) {
    for ( n1sv = 1; n1sv <= n2sv; n1sv ++ ) {
      index[ ioffindex + i - 1 ] = ind2;
      i ++;
      ind2 += dtrd2;
    }
  } else {
    for ( n2sv = 1; n2sv <= n1sv; n2sv ++ ) {
      index[ ioffindex + i - 1 ] = ind1;
      i ++;
      ind1 += dtrd1;
    }
  }
}
//*************************************************************************
LaPack0.dlaneg = function( n, d, lld, sigma, pivmin, r, ioffd,
iofflld ) {
  var blklen = 128;
  var negcnt = 0;
  var t = - sigma;
  for ( var bj = 1; bj <= r - 1; bj += blklen ) {
    var neg1 = 0;
    var bsav = t;
    for ( var j = bj; j <= Math.min( bj + blklen - 1, r - 1 );
    j ++ ) {
      var dplus = d[ ioffd + j - 1 ] + t;
      if ( dplus < 0. ) neg1 ++;
      var tmp = t / dplus;
      t = tmp * lld[ iofflld + j - 1 ] - sigma;
    }
    var sawnan = isNaN( t );
    if ( sawnan ) {
      neg1 = 0;
      t = bsav;
      for ( j = bj; j <= Math.min( bj + blklen - 1, r - 1 ); j ++ ) {
        dplus = d[ ioffd + j - 1 ] + t;
        if ( dplus < 0. ) neg1 ++;
        tmp = t / dplus;
        if ( isNaN( tmp ) ) tmp = 1.;
        t = tmp * lld[ iofflld + j - 1 ] - sigma;
      }
    }
    negcnt += neg1;
  }
  var p = d[ ioffd + n - 1 ] - sigma;
  for ( bj = n - 1; bj >= r; bj -= blklen ) {
    var neg2 = 0;
    bsav = p;
    for ( j = bj; j >= Math.max( bj - blklen + 1, r ); j -- ) {
      var dminus = lld[ iofflld + j - 1 ] + p;
      if ( dminus < 0. ) neg2 ++;
      tmp = p / dminus;
      p = tmp * d[ ioffd + j - 1 ] - sigma;
    }
    sawnan = isNaN( p );
    if ( sawnan ) {
      neg2 = 0;
      p = bsav;
      for ( j = bj; j >= Math.max( bj - blklen + 1, r ); j -- ) {
        dminus = lld[ iofflld + j - 1 ] + p;
        if ( dminus < 0. ) neg2 ++;
        tmp = p / dminus;
        if ( isNaN( tmp ) ) tmp = 1.;
        p = tmp * d[ ioffd + j - 1 ] - sigma;
      }
    }
    negcnt += neg2;
  }
  var gamma = ( t + sigma ) + p;
  if ( gamma < 0. ) negcnt ++;
  return negcnt;
}
//*************************************************************************
LaPack0.dlapmr = function( forwrd, m, n, X, ldx, k, ioffx,
ioffk ) {
  if ( m <= 1 ) return;
  for ( var i = 1; i <= m; i ++ ) {
    k[ ioffk + i - 1 ] = - k[ ioffk + i - 1 ];
  }
  if ( forwrd ) {
    for ( i = 1; i <= m; i ++ ) {
      if ( k[ ioffk + i - 1 ] > 0 ) continue;
      var j = i;
      k[ ioffk + j - 1 ] = - k[ ioffk + j - 1 ];
      var in2 = k[ ioffk + j - 1 ];
      while ( true ) { // 20
        if ( k[ ioffk + in2 - 1 ] > 0 ) break;
        for ( var jj = 1; jj <= n; jj ++ ) {
          var temp = X[ ioffx + j - 1 + ( jj - 1 ) * ldx ];
          X[ ioffx + j - 1 + ( jj - 1 ) * ldx ] =
            X[ ioffx + in2 - 1 + ( jj - 1 ) * ldx ];
          X[ ioffx + in2 - 1 + ( jj - 1 ) * ldx ] = temp;
        } // 30
        k[ ioffk + in2 - 1 ] = - k[ ioffk + in2 - 1 ];
        j = in2;
        in2 = k[ ioffk + in2 - 1 ];
      }
    } // 50
  } else {
    for ( i = 1; i <= m; i ++ ) {
      if ( k[ ioffk + i - 1 ] > 0 ) continue;
      k[ ioffk + i - 1 ] = - k[ ioffk + i - 1 ];
      j = k[ ioffk + i - 1 ];
      while ( true ) { // 60
        if ( j == i ) break;
        for ( jj = 1; jj <= n; jj ++ ) {
          temp = X[ ioffx + i - 1 + ( jj - 1 ) * ldx ];
          X[ ioffx + i - 1 + ( jj - 1 ) * ldx ] =
            X[ ioffx + j - 1 + ( jj - 1 ) * ldx ];
          X[ ioffx + j - 1 + ( jj - 1 ) * ldx ] = temp;
        } // 70
        k[ ioffk + j - 1 ] = - k[ ioffk + j - 1 ];
        j = k[ ioffk + j - 1 ];
      }
    } // 90
  }
}
//*************************************************************************
LaPack0.dlapmt = function( forwrd, m, n, X, ldx, k, ioffx,
ioffk ) {
  if ( n <= 1 ) return;
  for ( var i = 1; i <= n; i ++ ) {
    k[ ioffk + i - 1 ] = - k[ ioffk + i - 1 ];
  }
  if ( forwrd ) {
    for ( i = 1; i <= n; i ++ ) {
      if ( k[ ioffk + i - 1 ] <= 0 ) {
        var j = i;
        k[ ioffk + j - 1 ] = - k[ ioffk + j - 1 ];
        var in2 = k [ ioffk + j - 1 ];
        while ( k[ ioffk + in2 - 1 ] <= 0 ) {
          for ( var ii = 1; ii <= m; ii ++ ) {
            var temp = X[ ioffx + ii - 1 + ( j - 1 ) * ldx ];
            X[ ioffx + ii - 1 + ( j - 1 ) * ldx ] =
              X[ ioffx + ii - 1 + ( in2 - 1 ) * ldx ];
            X[ ioffx + ii - 1 + ( in2 - 1 ) * ldx ] = temp;
          }
          k[ ioffk + in2 - 1 ] = - k[ ioffk + in2 - 1 ];
          j = in2;
          in2 = k[ ioffk + in2 - 1 ];
        }
      }
    }
  } else {
    for ( i = 1; i <= n; i ++ ) {
      if ( k[ ioffk + i - 1 ] <= 0 ) {
        k[ ioffk + i - 1 ] = - k[ ioffk + i - 1 ];
        j = k[ ioffk + i - 1 ];
        while ( j != i ) {
          for ( ii = 1; ii <= m; ii ++ ) {
            temp = X[ ioffx + ii - 1 + ( i - 1 ) * ldx ];
            X[ ioffx + ii - 1 + ( i - 1 ) * ldx ] =
              X[ ioffx + ii - 1 + ( j - 1 ) * ldx ];
            X[ ioffx + ii - 1 + ( j - 1 ) * ldx ] = temp;
          }
          k[ ioffk + j - 1 ] = - k[ ioffk + j - 1 ];
          j = k[ ioffk + j - 1 ];
        }
      }
    }
  }
}
LaPack0.zlapmt = function( forwrd, m, n, X, ldx, k ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dlapy2 = function( x, y ) {
  var xabs = Math.abs( x );
  var yabs = Math.abs( y );
  var w = Math.max( xabs, yabs );
  var z = Math.min( xabs, yabs );
  return ( z == 0. ? w : w * Math.sqrt( 1. + Math.pow( z / w , 2) ) );
}
//*************************************************************************
LaPack0.dlapy3 = function( x, y, z ) {
  var xabs = Math.abs( x );
  var yabs = Math.abs( y );
  var zabs = Math.abs( z );
  var w = Math.max( Math.max( xabs, yabs ), zabs );
  if ( w == 0. ) return xabs + yabs + zabs;
  return w * Math.sqrt( Math.pow( xabs / w, 2 )
    + Math.pow( yabs / w, 2 ) + Math.pow( zabs / w , 2 ) );
}
//*************************************************************************
LaPack0.dlaqr1 = function( n, H, ldh, sr1, si1, sr2, si2, v,
ioffh, ioffv) {
  if ( n == 2 ) {
    var s = Math.abs( H[ ioffh ] - sr2 )
      + Math.abs( si2 ) + Math.abs( H[ ioffh + 1 ] );
    if ( s == 0. ) {
      v[ ioffv ] = 0.;
      v[ ioffv + 1 ] = 0.;
    } else {
      var h21s = H[ ioffh + 1 ] / s;
      v[ ioffv ] = h21s * H[ ioffh + ldh ]
        + ( H[ ioffh ] - sr1 ) * ( ( H[ ioffh ] - sr2 ) / s )
        - si1 * ( si2 / s );
      v[ ioffv + 1 ] = h21s * ( H[ ioffh ] + H[ ioffh + 1 + ldh ]
        - sr1 - sr2 );
    }
  } else {
    s = Math.abs( H[ ioffh ] - sr2 ) + Math.abs( si2 )
      + Math.abs( H[ ioffh + 1 ] ) + Math.abs( H[ ioffh + 2 ] );
    if ( s == 0. ) {
      v[ ioffv ] = 0.;
      v[ ioffv + 1 ] = 0.;
      v[ ioffv + 2 ] = 0.;
    } else {
      h21s = H[ ioffh + 1 ] / s;
      var h31s = H[ ioffh + 2 ] / s;
      v[ ioffv ] = ( H[ ioffh ] - sr1 ) * ( ( H[ ioffh ] - sr2 ) / s )
        - si1 * ( si2 / s ) + H[ ioffh + ldh ] * h21s
        + H[ ioffh + 2 * ldh ] * h31s;
      v[ ioffv + 1 ] = h21s * ( H[ ioffh ] + H[ ioffh + 1 + ldh ]
        - sr1 - sr2 ) + H[ ioffh + 1 + 2 * ldh ] * h31s;
      v[ ioffv + 2 ] = h31s * ( H[ ioffh ] + H[ ioffh + 2 + 2 * ldh ]
        - sr1 - sr2 ) + h21s * H[ ioffh + 2 + ldh ];
    }
  }
}
//*************************************************************************
LaPack0.dlar2v = function( n, x, y, z, incx, c, s, incc, ioffx,
ioffy, ioffz, ioffc, ioffs ) {
  var ix = 1;
  var ic = 1;
  for ( var i = 1; i <= n; i ++ ) {
    var xi = x[ ioffx + ix - 1 ];
    var yi = y[ ioffy + ix - 1 ];
    var zi = z[ ioffz + ix - 1 ];
    var ci = c[ ioffc + ic - 1];
    var si = s[ ioffs + ic - 1];
    var t1 = si * zi;
    var t2 = ci * zi;
    var t3 = t2 - si * xi;
    var t4 = t2 + si * yi;
    var t5 = ci * xi + t1;
    var t6 = ci * yi - t1;
    x[ ioffx + ix - 1 ] = ci * t5 + si * t4;
    y[ ioffy + ix - 1 ] = ci * t6 - si * t3;
    z[ ioffz + ix - 1 ] = ci * t4 - si * t5;
    ix += incx;
    ic += incc;
  }
}
LaPack0.zlar2v = function( n, x, y, z, incx, c, s, incc ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dlarft = function( direct, storev, n, k, V, ldv, tau,
T, ldt, ioffv, iofftau, iofft ) {
  if ( n == 0 ) return;
  if ( direct.charAt(0).toUpperCase() == 'F' ) {
    var prevlastv = n;
    for ( var i = 1; i <= k; i ++ ) {
      prevlastv = Math.max( i, prevlastv );
      if ( tau[ iofftau + i - 1 ] == 0. ) {
        for ( var j = 1; j <= i; j ++ ) {
          T[ iofft + j - 1 + ( i - 1 ) * ldt ] = 0.;
        }
      } else {
        var vii = V[ ioffv + i - 1 + ( i - 1 ) * ldv ];
        V[ ioffv + i - 1 + ( i - 1 ) * ldv ] = 1.;
        if ( storev.charAt(0).toUpperCase() == 'C' ) {
//            V = [V1] kxk unit lower triangular
//                [V2] (n-k)xk
          for ( var lastv = n; lastv >= i + 1; lastv -- ) {
            if ( V[ ioffv + lastv - 1 + ( i - 1 ) * ldv ] != 0.) break;
          }
          j = Math.min( lastv, prevlastv );
          Blas2.dgemv( 'Transpose', j - i + 1, i - 1,
            - tau[ iofftau + i - 1 ], V, ldv, V, 1, 0., T, 1,
            ioffv + i - 1, ioffv + i - 1 + ( i - 1 ) * ldv,
            iofft + ( i - 1 ) * ldt );
        } else {
//            V = [V1,V2] kxk unit upper triangular, kx(n-k)
          for ( lastv = n; lastv >= i + 1; lastv -- ) {
            if ( V[ ioffv + i - 1 + ( lastv - 1 ) * ldv ] != 0.) break;
          }
          j = Math.min( lastv, prevlastv );
          Blas2.dgemv( 'No transpose', i - 1, j - i + 1,
            - tau[ iofftau + i - 1 ], V, ldv, V, ldv, 0., T, 1,
            ioffv + ( i - 1 ) * ldv,
            ioffv + i - 1 + ( i - 1 ) * ldv, iofft + ( i - 1 ) * ldt );
        }
        V[ ioffv + i - 1 + ( i - 1 ) * ldv ] = vii;
        Blas2.dtrmv( 'Upper', 'No transpose', 'Non-unit', i - 1, T,
          ldt, T, 1, iofft, iofft + ( i - 1 ) * ldt );
        T[ iofft + i - 1 + ( i - 1 ) * ldt ] = tau[ iofftau + i - 1 ];
        if ( i > 1 ) prevlastv = Math.max( prevlastv, lastv );
        else prevlastv = lastv;
      }
    }
  } else {
    prevlastv = 1;
    for ( i = k; i >= 1; i -- ) {
      if ( tau[ iofftau + i - 1 ] == 0. ) {
        for ( j = i; j <= k; j ++ ) {
          T[ iofft + j - 1 + ( i - 1 ) * ldt ] = 0.;
        }
      } else {
        if ( i < k ) {
          if ( storev.charAt(0).toUpperCase() == 'C' ) {
//              V = [V1] (n-k)xk
//                  [V2] kxk unit upper triangular
            vii = V[ ioffv + n - k + i - 1 + ( i - 1 ) * ldv ];
            V[ ioffv + n - k + i - 1 + ( i - 1 ) * ldv ] = 1.;
            for ( lastv = 1; lastv <= i - 1; lastv ++ ) {
              if ( V[ ioffv + lastv - 1 + ( i - 1 ) * ldv ] != 0. ) {
                break;
              }
            }
            j = Math.max( lastv, prevlastv );
            Blas2.dgemv( 'Transpose', n - k + i - j + 1, k - i,
              - tau[ iofftau + i - 1 ], V, ldv, V, 1, 0., T, 1,
              ioffv + j - 1 + i * ldv, ioffv + j - 1 + ( i - 1 ) * ldv,
              iofft + i + ( i - 1 ) * ldt );
            V[ ioffv + n - k + i - 1 + ( i - 1 ) * ldv ] = vii;
          } else {
//              V = [V1,V2] kx(n-k), kxk unit lower triangular
            vii = V[ ioffv + i - 1 + ( n - k + i - 1 ) * ldv ];
            V[ ioffv + i - 1 + ( n - k + i - 1 ) * ldv ] = 1.;
            for ( lastv = 1; lastv <= i - 1; lastv ++ ) {
              if ( V[ ioffv + i - 1 + ( lastv - 1 ) * ldv ] != 0. ) {
                break;
              }
            }
            j = Math.max( lastv, prevlastv );
            Blas2.dgemv( 'No transpose', k - i, n - k + i - j + 1,
              - tau[ iofftau + i - 1 ], V, ldv, V, ldv, 0., T, 1,
              ioffv + i + ( j - 1 ) * ldv,
              ioffv + i - 1 + ( j - 1 ) * ldv,
              iofft + i + ( i - 1 ) * ldt );
            V[ ioffv + i - 1 + ( n - k + i - 1 ) * ldv ] = vii;
          }
          Blas2.dtrmv( 'Lower', 'No transpose', 'Non-unit', k - i,
            T, ldt, T, 1, iofft + i + i * ldt,
            iofft + i + ( i - 1 ) * ldt );
          if ( i > 1 ) prevlastv = Math.min( prevlastv, lastv );
          else prevlastv = lastv;
        }
        T[ iofft + i - 1 + ( i - 1 ) * ldt ] = tau[ iofftau + i - 1 ];
      }
    }
  }
}
LaPack0.zlarft = function( direct, storev, n, k, V, ldv, tau,
t, ldt ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dlargv = function( n, x, incx, y, incy, c, incc, ioffx,
ioffy, ioffc ) {
  var ix = 1;
  var iy = 1;
  var ic = 1;
  for ( var i = 1; i <= n; i ++ ) {
    var f = x[ ioffx + ix - 1 ];
    var g = y[ ioffy + iy - 1 ];
    if ( g == 0. ) {
      c[ ioffc + ic - 1 ] = 1.;
    } else if ( f == 0. ) {
      c[ ioffc + ic - 1 ] = 0.;
      y[ ioffy + iy - 1 ] = 1.;
      x[ ioffx + ix - 1 ] = g;
    } else if ( Math.abs( f ) > Math.abs( g ) ) {
      var t = g / f;
      var tt = Math.sqrt( 1. + t * t );
      c[ ioffc + ic - 1 ] = 1. / tt;
      y[ ioffy + iy - 1 ] = t * c[ ioffc + ic - 1 ];
      x[ ioffx + ix - 1 ] = f * tt;
    } else {
      t = f / g;
      tt = Math.sqrt( 1. + t * t );
      y[ ioffy + iy - 1 ] = 1. / tt;
      c[ ioffc + ic - 1 ] = t * y[ ioffy + iy - 1 ];
      x[ ioffx + ix - 1 ] = g * tt;
    }
    ic += incc;
    iy += incy;
    ix += incx;
  }
}
LaPack0.zlargv = function( n, x, incx, y, incy, c, incc ) {
  throw new Error("not programmed: complex matrix");
}
//************************************************************************
LaPack0.dlarnv = function( idist, iseed, n, x, ioffiseed,
ioffx ) {
  throw new Error("not tested: different from LaPack");
  if ( idist == 1 ) {
    for ( var i = 1; i <= n; i ++ ) {
      x[ ioffx + i - 1 ] = Math.random();
    }
  } else if ( idist == 2 ) {
    for ( i = 1; i <= n; i ++ ) {
      x[ ioffx + i - 1 ] = 2. * Math.random() - 1.;
    }
  } else if ( idist == 3 ) {
    for ( i = 1; i <= n; i ++ ) {
      x[ ioffx + i - 1 ] = Math.sqrt( -2. * Math.log( Math.random() ) )
        * Math.cos( 2. * Math.PI * Math.random() );
    }
  }
}
//*************************************************************************
LaPack0.dlarra = function( n, d, e, e2, spltol, tnrm, nsplit,
isplit, info, ioffd, ioffe, ioffe2, ioffisplit) {
  info.setValue( 0 );
  nsplit.setValue( 1 );
  if ( spltol < 0. ) {
    var tmp1 = Math.abs( spltol ) * tnrm;
    for ( var i = 1; i <= n - 1; i ++ ) {
      var eabs = Math.abs( e[ ioffe + i - 1 ] );
      if ( eabs <= tmp1 ) {
        e[ ioffe + i - 1 ] = 0.;
        e2[ ioffe2 + i - 1 ] = 0.;
        isplit[ ioffisplit + nsplit.getValue() - 1 ] = i;
        nsplit.setValue( nsplit.getValue() + 1 );
      }
    }
  } else {
    for ( i = 1; i <= n - 1; i ++ ) {
      eabs = Math.abs( e[ ioffe + i - 1 ] );
      if ( eabs <= spltol * Math.sqrt( Math.abs( d[ ioffd + i - 1 ] ) )
      * Math.sqrt( Math.abs( d[ ioffd + i ] ) ) ) {
        e[ ioffe + i - 1 ] = 0.;
        e2[ ioffe2 + i - 1 ] = 0.;
        isplit[ ioffisplit + nsplit.getValue() - 1 ] = i;
        nsplit.setValue( nsplit.getValue() + 1 );
      }
    }
  }
  isplit[ ioffisplit + nsplit.getValue() - 1 ] = n;
}
//*************************************************************************
LaPack0.dlarrc = function( jobt, n, vl, vu, d, e, pivmin,
eigcnt, lcnt, rcnt, info, ioffd, ioffe) {
  info.setValue( 0 );
  lcnt.setValue( 0 );
  rcnt.setValue( 0 );
  eigcnt.setValue( 0 );
  var matt = ( jobt.charAt(0).toUpperCase() == 'T' );
  if ( matt ) {
    var lpivot = d[ ioffd ] - vl;
    var rpivot = d[ ioffd ] - vu;
    if ( lpivot <= 0. ) lcnt.setValue( lcnt.getValue() + 1 );
    if ( rpivot <= 0. ) rcnt.setValue( rcnt.getValue() + 1 );
    for ( var i = 1; i <= n - 1; i ++ ) {
      var tmp = Math.pow( e[ ioffe  + i - 1 ], 2 );
      lpivot = ( d[ ioffd + i ] - vl ) - tmp / lpivot;
      rpivot = ( d[ ioffd + i ] - vu ) - tmp / rpivot;
      if ( lpivot <= 0. ) lcnt.setValue( lcnt.getValue() + 1 );
      if ( rpivot <= 0. ) rcnt.setValue( rcnt.getValue() + 1 );
    }
  } else {
    var sl = - vl;
    var su = - vu;
    for ( i = 1; i <= n - 1; i ++ ) {
      lpivot = d[ ioffd + i - 1 ] + sl;
      rpivot = d[ ioffd + i - 1 ] + su;
      if ( lpivot <= 0. ) lcnt.setValue( lcnt.getValue() + 1 );
      if ( rpivot <= 0. ) rcnt.setValue( rcnt.getValue() + 1 );
      tmp = e[ ioffe + i - 1 ] * d[ ioffd + i - 1 ]
        * e[ ioffe + i - 1 ];
      var tmp2 = tmp / lpivot;
      sl = ( tmp2 == 0. ? tmp - vl : sl * tmp2 - vl );
      tmp2 = tmp / rpivot;
      su = ( tmp2 == 0. ? tmp - vu : su * tmp2 - vu );
    }
    lpivot = d[ ioffd + n - 1 ] + sl;
    rpivot = d[ ioffd + n - 1 ] + su;
    if ( lpivot <= 0. ) lcnt.setValue( lcnt.getValue() + 1 );
    if ( rpivot <= 0. ) rcnt.setValue( rcnt.getValue() + 1 );
  }
  eigcnt.setValue( rcnt.getValue() - lcnt.getValue() );
}
//*************************************************************************
LaPack0.dlarrj = function( n, d, e2, ifirst, ilast, rtol,
offset, w, werr, work, iwork, pivmin, spdiam, info, ioffd, ioffe2, ioffw,
ioffwerr, ioffwork, ioffiwork ) {
  info.setValue( 0 );
  var maxitr = Math.round( ( Math.log( spdiam + pivmin )
    - Math.log( pivmin ) ) / Math.log( 2. ) ) + 2; 
  var i1 = ifirst;
  var i2 = ilast;
  var nint = 0;
  var prev = 0;
  for ( var i = i1; i <= i2; i ++ ) {
    var k = 2 * i
    var ii = i - offset;
    var left = w[ ioffw + ii - 1 ] - werr[ ioffwerr + ii - 1 ];
    var mid = w[ ioffw + ii - 1 ];
    var right = w[ ioffw + ii - 1 ] + werr[ ioffwerr + ii - 1 ];
    var width = right - mid;
    var tmp = Math.max( Math.abs( left ), Math.abs( right ) );
    if ( width < rtol * tmp ) {
      iwork[ ioffiwork + k - 2 ] = -1;
      if ( i == i1 && i < i2 ) i1 = i + 1;
      if ( prev >= i1 && i <= i2 ) {
        iwork[ ioffiwork + 2 * prev - 2 ] = i + 1;
      }
    } else {
      prev = i;
      var fac = 1.;
      while ( true ) { // 20
        var cnt = 0;
        var s = left;
        var dplus = d[ ioffd ] - s;
        if ( dplus < 0. ) cnt ++;
        for ( var j = 2; j <= n; j ++ ) {
          dplus = d[ ioffd + j - 1 ] - s
            - e2[ ioffe2 + j - 2 ] / dplus;
          if ( dplus < 0. ) cnt ++;
        }
        if ( cnt > i - 1 ) {
          left -= werr[ ioffwerr + ii - 1 ] * fac;
          fac *= 2.;
        } else break;
      }
      fac = 1.;
      while ( true ) { // 50
        cnt = 0;
        s = right;
        dplus = d[ ioffd ] - s;
        if ( dplus < 0. ) cnt ++;
        for ( j = 2; j <= n; j ++ ) {
          dplus = d[ ioffd + j - 1 ] - s
            - e2[ ioffe2 + j - 2 ] / dplus;
          if ( dplus < 0. ) cnt ++;
        } // 60
        if ( cnt < i ) {
          right += werr[ ioffwerr + ii - 1 ] * fac;
          fac *= 2.;
        } else break;
      }
      nint ++;
      iwork[ ioffiwork + k - 2 ] = i + 1;
      iwork[ ioffiwork + k - 1 ] = cnt;
    }
    work[ ioffwork + k - 2 ] = left;
    work[ ioffwork + k - 1 ] = right;
  } // 75
  var savi1 = i1;
  var iter = 0;
  while ( true ) { // 80
    prev = i1 - 1;
    i = i1;
    var olnint = nint;
    for ( var p = 1; p <= olnint; p ++ ) {
      k = 2 * i;
      ii = i - offset;
      var next = iwork[ ioffiwork + k - 2 ];
      left = work[ ioffwork + k - 2 ];
      right = work[ ioffwork + k - 1 ];
      mid = 0.5 * ( left + right );
      width = right - mid;
      tmp = Math.max( Math.abs( left ), Math.abs( right ) );
      if ( width < rtol * tmp || iter == maxitr ) {
        nint --;
        iwork[ ioffiwork + k - 2 ] = 0;
        if ( i1 == i ) i1 = next;
        else {
          if ( prev >= i1 ) iwork[ ioffiwork + 2 * prev - 2 ] = next;
        }
        i = next;
        continue;
      }
      prev = i;
      cnt = 0;
      s = mid;
      dplus = d[ ioffd ] - s;
      if ( dplus < 0. ) cnt ++;
      for ( j = 2; j <= n; j ++ ) {
        dplus = d[ ioffd + j - 1 ] - s - e2[ ioffe2 + j - 2 ] / dplus;
        if ( dplus < 0. ) cnt ++;
      }
      if ( cnt <= i - 1 ) work[ ioffwork + k - 2 ] = mid;
      else work[ ioffwork + k - 1 ] = mid;
      i = next;
    } // 100
    iter ++;
    if ( nint > 0 && iter <= maxitr ) continue;
    else break;
  } // 80
  for ( i = savi1; i <= ilast; i ++ ) {
    k = 2 * i;
    ii = i - offset;
    if ( iwork[ ioffiwork + k - 2 ] == 0 ) {
      w[ ioffw + ii - 1 ] =
        0.5 * ( work[ ioffwork + k - 2 ] + work[ ioffwork + k - 1 ] );
      werr[ ioffwerr + ii - 1 ] = work[ ioffwork + k - 1 ]
        - w[ ioffw + ii - 1 ];
    }
  } // 110
}
//*************************************************************************
//  will eventually be replaced by BLAS_dge_diag_scale
LaPack0.dlarscl2 = function( m, n, d, X, ldx, ioffd, ioffx ) {
  for ( var j = 1; j <= n; j ++ ) {
    for ( var i = 1; i <= m; i ++ ) {
      X[ ioffx + i - 1 + ( j - 1 ) * ldx ] /= d[ ioffd + i - 1 ];
    }
  }
}
//*************************************************************************
LaPack0.dlartv = function( n, x, incx, y, incy, c, s, incc,
ioffx, ioffy, ioffc, ioffs ) {
  var ix = 1;
  var iy = 1;
  var ic = 1;
  for ( var i = 1; i <= n; i ++ ) {
    var xi = x[ ioffx + ix - 1 ];
    var yi = y[ ioffy + iy - 1 ];
    x[ ioffx + ix - 1 ] = c[ ioffc + ic - 1 ] * xi
      + s[ ioffs + ic - 1 ] * yi;
    y[ ioffy + iy - 1 ] = c[ ioffc + ic - 1 ] * yi
      - s[ ioffs + ic - 1 ] * xi;
    ix += incx;
    iy += incy;
    ic += incc;
  }
}
LaPack0.zlartv = function( n, x, incx, y, incy, c, s, incc ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dlaruv = function( iseed, n, x, ioffx) {
  var lv = 128;
  for ( var i = 1; i <= Math.min( n, lv ); i ++ ) {
    x[ ioffx + i - 1 ] = Math.random();
  }
  for ( i = 1; i <= 4; i ++ ) iseed[ i - 1 ] = i;
}
//*************************************************************************
LaPack0.dlarz = function( side, m, n, l, v, incv, tau, C, ldc,
work, ioffv, ioffC, ioffwork ) {
  if ( side.charAt(0).toUpperCase() == 'L' ) {
    if ( tau != 0. ) {
      Blas1.dcopy( n, C, ldc, work, 1, ioffC, ioffwork );
      Blas2.dgemv( 'Transpose', l, n, 1., C, ldc, v, incv, 1., work, 1,
        ioffC + m - l, ioffv, ioffwork );
      Blas1.daxpy( n, -tau, work, 1, C, ldc, ioffwork, ioffC );
      Blas2.dger( l, n, -tau, v, incv, work, 1, C, ldc, ioffv,
        ioffwork, ioffC + m - l );
    }
  } else {
    if ( tau != 0. ) {
      Blas1.dcopy( m, C, 1, work, 1, ioffC, ioffwork );
      Blas2.dgemv( 'No transpose', m, l, 1., C, ldc, v, incv, 1.,
        work, 1, ioffC + ( n - l ) * ldc, ioffv, ioffwork );
      Blas1.daxpy( m, -tau, work, 1, C, 1, ioffwork, ioffC );
      Blas2.dger( m, l, -tau, work, 1, v, incv, C, ldc, ioffwork,
        ioffv, ioffC + ( n - l ) * ldc );
    }
  }
}
//*************************************************************************
LaPack0.dlarzb = function( side, trans, direct, storev, m, n, k,
l, V, ldv, T, ldt, C, ldc, Work, ldwork, ioffV, ioffT, ioffC, ioffWork) {
  if ( m <= 0 || n <= 0 ) return;
  var info = 0;
  if ( direct.charAt(0).toUpperCase() != 'B' ) info = -3;
  else if ( storev.charAt(0).toUpperCase() != 'R' ) info = -4;
  if ( info != 0 ) {
    Blas2.xerbla( 'dlarzb', - info );
    return;
  }
  var transt =
    ( trans.charAt(0).toUpperCase() == 'N' ? 'T' : 'N' );
  if ( side.charAt(0).toUpperCase() == 'L' ) {
    for ( var j = 1; j <= k; j ++ ) {
      Blas1.dcopy( n, C, ldc, Work, 1, ioffC + j - 1,
        ioffWork + ( j - 1 ) * ldwork );
    }
    if ( l > 0 ) {
      Blas3.dgemm( 'Transpose', 'Transpose', n, k, l, 1., C, ldc,
        V, ldv, 1., Work, ldwork, ioffC + m - l, ioffV, ioffWork );
    }
    Blas3.dtrmm( 'Right', 'Lower', transt, 'Non-unit', n, k, 1., 
      T, ldt, Work, ldwork, ioffT, ioffWork );
    for ( j = 1; j <= n; j ++ ) {
      for ( var i = 1; i <= k; i ++ ) {
        C[ ioffC + i - 1 + ( j - 1 ) * ldc ] -=
          Work[ ioffWork + j - 1 + ( i - 1 ) * ldwork ];
      }
    }
    if ( l > 0 ) {
      Blas3.dgemm( 'Transpose', 'Transpose', l, n, k, -1., V, ldv,
        Work, ldwork, 1, C, ldc, ioffV, ioffWork, ioffC + m - l );
    }
  } else if ( side.charAt(0).toUpperCase() == 'R' ) {
    for ( j = 1; j <= k; j ++ ) {
      Blas1.dcopy( m, C, 1, Work, 1, ioffC + ( j - 1 ) * ldc,
        ioffWork + ( j - 1 ) * ldwork );
    }
    if ( l > 0 ) {
      Blas3.dgemm( 'No transpose', 'Transpose', m, k, l, 1., C, ldc,
        V, ldv, 1., Work, ldwork, ioffC + ( n - l ) * ldc, ioffV,
        ioffWork );
    }
    Blas3.dtrmm( 'Right', 'Lower', trans, 'Non-unit', m, k, 1.,
      T, ldt, Work, ldwork, ioffT, ioffWork );
    for ( j = 1; j <= k; j ++ ) {
      for ( i = 1; i <= m; i ++ ) {
        C[ ioffC + i - 1 + ( j - 1 ) * ldc ] -=
          Work[ ioffWork + i - 1 + ( j - 1 ) * ldwork ];
      }
    }
    if ( l > 0 ) {
      Blas3.dgemm( 'No transpose', 'No transpose', m, l, k, -1.,
        Work, ldwork, V, ldv, 1., C, ldc, ioffWork, ioffV,
        ioffC + ( n - l ) * ldc );
    }
  }
}
LaPack0.zlarzb = function( side, trans, direct, storev, m, n, k,
l, V, ldv, T, ldt, C, ldc, Work, ldwork, ioffV, ioffT, ioffC, ioffWork) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dlarzt = function( direct, storev, n, k, V, ldv, tau, T,
ldt, ioffV, iofftau, ioffT ) {
  var info = 0;
  if ( direct.charAt(0).toUpperCase() != 'B' ) info = -1;
  else if ( storev.charAt(0).toUpperCase() != 'R' ) info = -2;
  if ( info != 0 ) {
    Blas2.xerbla( 'dlarzt', - info );
    return;
  }
  for ( var i = k; i >= 1; i -- ) {
    if ( tau[ iofftau + i - 1 ] == 0. ) {
      for ( var j = i; j <= k; j ++ ) {
        T[ ioffT + j - 1 + ( i - 1 ) * ldt ] = 0.;
      }
    } else {
      if ( i < k ) {
        Blas2.dgemv( 'No transpose', k - i, n,
          - tau[ iofftau + i - 1 ], V, ldv, V, ldv, 0., T, 1,
          ioffV + i, ioffV + i - 1, ioffT + i + ( i - 1 ) * ldt );
        Blas2.dtrmv( 'Lower', 'No transpose', 'Non-unit', k - i,
          T, ldt, T, 1, ioffT + i + i * ldt,
          ioffT + i + ( i - 1 ) * ldt );
      }
      T[ ioffT + i - 1 + ( i - 1 ) * ldt ] = tau[ iofftau + i - 1 ];
    }
  }
}
LaPack0.zlarzt = function( direct, storev, n, k, V, ldv, tau, T,
ldt, ioffV, iofftau, ioffT ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dlas2 = function( f, g, h, ssminReference,
ssmaxReference) {
  var fa = Math.abs( f );
  var ga = Math.abs( g );
  var ha = Math.abs( h );
  var fhmn = Math.min( fa, ha );
  var fhmx = Math.max( fa, ha );
  if ( fhmn == 0. ) {
    ssminReference.setValue( 0. );
    ssmaxReference.setValue( ( fhmx == 0. ? ga :
      Math.max( fhmx, ga ) * Math.sqrt( 1.  + Math.pow(
        Math.min( fhmx, ga ) / Math.max( fhmx, ga ) , 2 ) ) ) );
  } else {
    if ( ga < fhmx ) {
      var as2 = 1. + fhmn / fhmx ;
      var at = ( fhmx - fhmn ) / fhmx;
      var au = Math.pow( ga / fhmx , 2 );
      var c = 2. / ( Math.sqrt( as2 * as2 + au )
                          + Math.sqrt( at * at + au ) );
      ssminReference.setValue( fhmn * c );
      ssmaxReference.setValue( fhmx / c );
    } else {
      au = fhmx / ga;
      if ( au == 0. ) {
        ssminReference.setValue( ( fhmn * fhmx ) / ga );
        ssmaxReference.setValue( ga );
      } else {
        as2 = 1. + fhmn / fhmx;
        at = ( fhmx - fhmn ) / fhmx;
        c = 1. / ( Math.sqrt( 1. + Math.pow( as2 * au , 2 ) )
                 + Math.sqrt( 1. + Math.pow( at * au , 2 ) ) );
        ssminReference.setValue( ( fhmn * c ) * au );
        ssminReference.setValue( ssminReference.getValue()
                               + ssminReference.getValue() );
        ssmaxReference.setValue( ga / ( c + c ) );
      }
    }
  }
}
//*************************************************************************
LaPack0.dlascl2 = function( m, n, d, X, ldx, ioffd, ioffx ) {
  for ( var j = 1; j <= n; j ++ ) {
    for ( var i = 1; i <= m; i ++ ) {
      X[ ioffx + i - 1 + ( j - 1 ) * ldx ] *= d[ ioffd + i - 1 ];
    }
  }
}
//*************************************************************************
LaPack0.dlasd5 = function( i, d, z, delta, rho, dsigmaReference,
work, ioffd, ioffz, ioffdelta, ioffwork ) {
  var del = d[ ioffd + 1 ] - d[ ioffd ];
  var delsq = del * ( d[ ioffd + 1 ] + d[ ioffd ] );
  if ( i == 1 ) {
    var w = 1. + 4. * rho * ( z[ ioffz + 1 ] * z[ ioffz + 1 ]
      / ( d[ ioffd ] + 3. * d[ ioffd + 1 ] )
      - z[ ioffz ] * z[ ioffz ]
      / ( 3. * d[ ioffd ] + d[ ioffd + 1 ] ) ) / del;
    if ( w > 0. ) {
      var b = delsq + rho * ( z[ ioffz ] * z[ ioffz ]
        + z[ ioffz + 1 ] * z[ ioffz + 1 ] );
      var c = rho * z[ ioffz ] * z[ ioffz ] * delsq;
      var tau = 2. * c
        / ( b + Math.sqrt( Math.abs( b * b - 4. * c ) ) );
      tau /= d[ ioffd ] + Math.sqrt( d[ ioffd ] * d[ ioffd ] + tau );
      dsigmaReference.setValue( d[ ioffd ] + tau );
      delta[ ioffdelta ] = - tau;
      delta[ ioffdelta + 1 ] = del - tau;
      work[ ioffwork ] = 2. * d[ ioffd ] + tau;
      work[ ioffwork + 1 ] = ( d[ ioffd ] + tau ) + d[ ioffd + 1 ];
    } else {
      b = - delsq + rho * ( z[ ioffz ] * z[ ioffz ]
        + z[ ioffz + 1 ] * z[ ioffz + 1 ] );
      c = rho * z[ ioffz + 1 ] * z[ ioffz + 1 ] * delsq;
      tau = ( b > 0. ? -2 * c / ( b + Math.sqrt( b * b + 4. * c ) )
        : ( b - Math.sqrt( b * b + 4. * c ) ) / 2. );
      tau /= d[ ioffd + 1 ] + Math.sqrt(
        Math.abs( d[ ioffd + 1 ] * d[ ioffd + 1 ] + tau ) );
      dsigmaReference.setValue( d[ ioffd + 1 ] + tau );
      delta[ ioffdelta ] = - ( del + tau );
      delta[ ioffdelta + 1 ] = - tau;
      work[ ioffwork ] = d[ ioffd ] + tau + d[ ioffd + 1 ];
      work[ ioffwork + 1 ] = 2. * d[ ioffd + 1 ] + tau;
    }
  } else {
    b = - delsq + rho * ( z[ ioffz ] * z[ ioffz ]
      + z[ ioffz + 1 ] * z[ ioffz + 1 ] );
    c = rho * z[ ioffz + 1 ] * z[ ioffz + 1 ] * delsq;
    tau = ( b > 0. ? ( b + Math.sqrt( b * b + 4. * c ) ) / 2.
      : 2. * c / ( -b + Math.sqrt( b * b + 4. * c ) ) );
    tau /= d[ ioffd + 1 ]
      + Math.sqrt( d[ ioffd + 1 ] * d[ ioffd + 1 ] + tau );
    dsigmaReference.setValue( d[ ioffd + 1 ] + tau );
    delta[ ioffdelta ] = - ( del + tau );
    delta[ ioffdelta + 1 ] = - tau;
    work[ ioffwork ] = d[ ioffd ] + tau + d[ ioffd + 1 ];
    work[ ioffwork + 1 ] = 2. * d[ ioffd + 1 ] + tau;
  }
}
//*************************************************************************
LaPack0.dlasdt = function( n, lvl, nd, inode, ndiml, ndimr, msub,
ioffinode, ioffndiml, ioffndimr) {
  var maxn = Math.max( 1, n );
  var temp = Math.log( Number( maxn ) / Number( msub + 1 ) )
    / Math.log( 2. );
  lvl.setValue( Math.round( temp ) + 1 );
  var i = n / 2;
  inode[ ioffinode ] = i + 1;
  ndiml[ ioffndiml ] = i;
  ndimr[ ioffndimr ] = n - i - 1;
  var il = 0;
  var ir = 1;
  var llst = 1;
  for ( var nlvl = 1; nlvl <= lvl.getValue() - 1; nlvl ++ ) {
    for ( i = 0; i <= llst - 1; i ++ ) {
      il += 2;
      ir += 2;
      var ncrnt = llst + i;
      var val = ndiml[ ioffndiml + ncrnt - 1 ] / 2;
      ndiml[ ioffndiml + il - 1 ] = val;
      ndimr[ ioffndimr + il - 1 ] = ndiml[ ioffndiml + ncrnt - 1 ]
        - ndiml[ ioffndiml + il - 1 ] - 1;
      inode[ ioffinode + il - 1 ] = inode[ ioffinode + ncrnt - 1 ]
        - ndimr[ ioffndimr + il - 1 ] - 1;
      var kkl = ioffndiml + il - 1;
      var kkr = ioffndimr + il - 1;
      var kk = ioffinode + il - 1;
      val = ndimr[ ioffndimr + ncrnt - 1 ] / 2;
      ndiml[ ioffndiml + ir - 1 ] = val;
      ndimr[ ioffndimr + ir - 1 ] = ndimr[ ioffndimr + ncrnt - 1 ]
        - ndiml[ ioffndiml + ir - 1 ] - 1;
      inode[ ioffinode + ir - 1 ] = inode[ ioffinode + ncrnt - 1 ]
        + ndiml[ ioffndiml + ir - 1 ] + 1;
      kkl = ioffndiml + ir - 1;
      kkr = ioffndimr + ir - 1;
      kk = ioffinode + ir - 1;
    }
    llst *= 2;
  }
  nd.setValue( llst * 2 - 1 );
}
//*************************************************************************
LaPack0.dlaset = function( uplo, m, n, alpha, beta, A, lda,
ioffa ) {
  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    for ( var j = 2; j <= n; j ++ ) {
      for ( var i = 1; i <= Math.min( j - 1, m ); i ++ ) {
        A[ ioffa + i - 1 + ( j - 1 ) * lda ] = alpha;
      }
    }
  } else if ( uplo.charAt(0).toUpperCase() == 'L' ) {
    for ( j = 1; j <= Math.min( m, n ); j ++ ) {
      for ( i = j + 1; i <= m; i ++ ) {
        A[ ioffa + i - 1 + ( j - 1 ) * lda ] = alpha;
      }
    }
  } else {
    for ( j = 1; j <= n; j ++ ) {
      for ( i = 1; i <= m; i ++ ) {
        A[ ioffa + i - 1 + ( j - 1 ) * lda ] = alpha;
      }
    }
  }
  for ( i = 1; i <= Math.min( m, n ); i ++ ) {
    A[ ioffa + i - 1 + ( i - 1 ) * lda ] = beta;
  }
}
LaPack0.zlaset = function( uplo, m, n, alpha, beta, A, lda ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dlasq4 = function( i0, n0, z, pp, n0in, dmin, dmin1,
dmin2, dn, dn1, dn2, tauReference, ttype, gReference, ioffz ) {
  var cnst1 = 0.5630;
  var cnst2 = 1.01;
  var cnst3 = 1.05;
  var qurtr = 0.25;
  var third = 0.333;
  if ( dmin <= 0. ) {
    tau.setValue( - dmin );
    ttype.setValue( -1 );
    return;
  }
  var nn = 4 * n0 + pp;
  if ( n0in == n0 ) {
    if ( dmin == dn || dmin == dn1 ) {
      var b1 = Math.sqrt( z[ ioffz + nn - 4 ] )
        * Math.sqrt( z[ ioffz + nn - 6 ] );
      var b2 = Math.sqrt( z[ ioffz + nn - 8 ] )
        * Math.sqrt( z[ ioffz + nn - 10 ] );
      var a2 = z[ ioffz + nn - 8 ] + z[ ioffz + nn - 6 ];
      if ( dmin == dn && dmin1 == dn1 ) {
        var gap2 = dmin2 - a2 - dmin2 * qurtr;
        var gap1 = ( gap2 > 0. && gap2 > b2 ?
          a2 - dn - ( b2 / gap2 ) * b2 : a2 - dn - ( b1 + b2 ) );
        if ( gap1 > 0. && gap1 > b1 ) {
          var s =
            Math.max( dn - ( b1 / gap1 ) * b1, 0.5 * dmin );
          ttype.setValue( -2 );
        } else {
          s = 0.;
          if ( dn > b1 ) s = dn - b1;
          if ( a2 > ( b1 + b2 ) ) s = Math.min( s, a2 - ( b1 + b2 ) );
          s = Math.max( s, third * dmin );
          ttype.setValue( -3 );
        }
      } else {
        ttype.setValue( -4 );
        s = qurtr * dmin;
        if ( dmin == dn ) {
          var gam = dn;
          a2 = 0.
          if ( z[ ioffz + nn - 6 ] > z[ ioffz + nn - 8 ] ) return;
          b2 = z[ ioffz + nn - 6 ] / z[ ioffz + nn - 8 ];
          var np = nn - 9;
        } else {
          np = nn - 2 * pp;
          b2 = z[ ioffz + np - 3 ];
          gam = dn1;
          if ( z[ ioffz + np - 5 ] > z[ ioffz + np - 3 ] ) return;
          a2 = z[ ioffz + np - 5 ] / z[ ioffz + np - 3 ];
          if ( z[ ioffz + nn - 10 ] > z[ ioffz + nn - 12 ] ) return;
          b2 = z[ ioffz + nn - 10 ] / z[ ioffz + nn - 12 ];
          np = nn - 13;
        }
        a2 += b2;
        for ( var i4 = np; i4 >= 4 * i0 - 1 + pp; i4 -= 4 ) {
          if  ( b2 == 0. ) break;
          b1 = b2;
          if ( z[ ioffz + i4 - 1 ] > z[ ioffz + i4 - 3 ] ) return;
          b2 *= z[ ioffz + i4 - 1 ] / z[ ioffz + i4 -3 ];
          a2 += b2;
          if ( 100. * Math.max( b2, b1 ) < a2 || cnst1 < a2 ) break;
        }
        a2 *= cnst3;
        if ( a2 < cnst1 ) {
          s = gam * ( 1. - Math.sqrt( a2 ) ) / ( 1. + a2 );
        }
      }
    } else if ( dmin == dn2 ) {
      ttype.setValue( -5 );
      s = qurtr * dmin;
      np = nn - 2 * pp;
      b1 = z[ ioffz + np - 3 ];
      b2 = z[ ioffz + np - 7 ];
      gam = dn2;
      if ( z[ ioffz + np - 9 ] > b2 || z[ ioffz + np - 5 ] > b1 ) {
        return;
      }
      a2 = ( z[ ioffz + np - 9 ] / b2 )
        * ( 1. + z[ ioffz + np - 5 ] / b1 );
      if ( n0 - i0 > 2 ) {
        b2 = z[ ioffz + nn - 14 ] / z[ ioffz + nn - 16 ];
        a2 += b2;
        for ( i4 = nn - 17; i4 >= 4 * i0 - 1 + pp; i4 -= 4 ) {
          if ( b2 == 0. ) break;
          b1 = b2;
          if ( z[ ioffz + i4 - 1 ] > z[ ioffz + i4 - 3 ] ) return;
          b2 *= z[ ioffz + i4 - 1 ] / z[ ioffz + i4 - 3 ];
          a2 += b2;
          if ( 100. * Math.max( b2, b1 ) < a2 || cnst1 < a2 ) break;
        }
        a2 *= cnst3;
      }
      if ( a2 < cnst1 ) {
        s = gam * ( 1. - Math.sqrt( a2 ) ) / ( 1. + a2 );
      }
    } else {
      if ( ttype.getValue() == -6 ) {
        g.setValue( g.getValue() + third * ( 1. - g.getValue() ) );
      } else if ( ttype.getValue() == -18 ) g.setValue( qurtr * third );
      else g.setValue( qurtr );
      s = g.setValue( dmin );
      ttype.setValue( -6 );
    }
  } else if ( n0in == n0 + 1 ) {
    if ( dmin1 == dn1 && dmin2 == dn2 ) {
      ttype.setValue( -7 );
      s = third * dmin1;
      if ( z[ ioffz + nn - 6 ] > z[ ioffz + nn - 8 ] ) return;
      b1 = z[ ioffz + nn - 6 ] / z[ ioffz + nn - 8 ];
      b2 = b1;
      if ( b2 != 0. ) {
        for ( i4 = 4 * n0 - 9 + pp; i4 >= 4 * i0 - 1 + pp; i4 -= 4 ) {
          a2 = b1;
          if ( z[ ioffz + i4 - 1 ] > z[ ioffz + i4 - 3 ] ) return;
          b1 *= z[ ioffz + i4 - 1 ] / z[ ioffz + i4 - 3 ];
          b2 += b1;
          if ( 100. * Math.max( b1, a2 ) < b2 ) break;
        }
      }
      b2 = Math.sqrt( cnst3 * b2 );
      a2 = dmin1 / ( 1. + b2 * b2 );
      gap2 = 0.5 * dmin2 - a2;
      if ( gap2 > 0. && gap2 > b2 * a2 ) {
        s = Math.max( s,
          a2 * ( 1. - cnst2 * a2 * ( b2 / gap2 ) * b2 ) );
      } else {
        s = Math.max( s, a2 * ( 1. - cnst2 * b2 ) );
        ttype.setValue( -8 );
      }
    } else {
      s = qurtr * dmin1;
      if ( dmin1 == dn1 ) s = 0.5 * dmin1;
      ttype.setValue( -9 );
    }
  } else if ( n0in == n0 + 2 ) {
    if ( dmin2 == dn2 &&
    2. * z[ ioffz + nn - 6 ] < z[ ioffz + nn - 8 ] ) {
      ttype.setValue( -10 );
      s = third * dmin2;
      if ( z[ ioffz + nn - 6 ] > z[ ioffz + nn - 8 ] ) return;
      b1 = z[ ioffz + nn - 6 ] / z[ ioffz + nn - 8 ];
      b2 = b1;
      if ( b2 != 0. ) {
        for ( i4 = 4 * n0 - 9 + pp; i4 >= 4 * i0 -1 + pp; i4 -= 4 ) {
          if ( z[ ioffz + i4 - 1 ] > z[ ioffz + i4 - 3 ] ) return;
          b1 *= z[ ioffz + i4 - 1 ] / z[ ioffz + i4 - 3 ];
          b2 += b1;
          if ( 100. * b1 < b2 ) break;
        }
      }
      b2 = Math.sqrt( cnst3 * b2 );
      a2 = dmin2 / ( 1. + b2 * b2 );
      gap2 = z[ ioffz + nn - 8 ] + z[ ioffz + nn - 10 ]
        - Math.sqrt( z[ ioffz + nn - 12 ] )
        * Math.sqrt( z[ ioffz + nn - 10 ] ) - a2;
      s = ( gap2 > 0. && gap2 > b2 * a2 ?
        Math.max( s, a2 * ( 1. - cnst2 * a2 * ( b2 / gap2 ) * b2) ) :
        Math.max( s, a2 * ( 1. - cnst2 * b2 ) ) );
    } else {
      s = qurtr * dmin2;
      ttype.setValue( -11 );
    }
  } else if ( n0in > n0 + 2 ) {
    s = 0.;
    ttype.setValue( -12 );
  }
  tau.setValue( s );
}
//*************************************************************************
LaPack0.dlasq5 = function( i0, n0, z, pp, tau, dminReference,
dmin1Reference, dmin2Reference, dnReference, dnm1Reference, dnm2Reference,
ieee, ioffz ) {
  if ( ( n0 - i0 - 1 ) <= 0 ) return;
  var j4 = 4 * i0 + pp - 3;
  var emin = z[ ioffz + j4 + 3 ];
  var d = z[ ioffz + j4 - 1 ] - tau;
  dmin1.setValue( - z[ ioffz + j4 - 1 ] );
  if ( ieee ) {
    if ( pp == 0 ) {
      for ( j4 = 4 * i0; j4 <= 4 * ( n0 - 3 ); j4 += 4 ) {
        z[ ioffz + j4 - 3 ] = d + z[ ioffz + j4 - 2 ];
        var temp = z[ ioffz + j4 ] / z[ ioffz + j4 - 3 ];
        d = d * temp - tau;
        dmin.setValue( Math.min( dmin.getValue(), d ) );
        z[ ioffz + j4 - 1 ] = z[ ioffz + j4 - 2 ] * temp;
        emin = Math.min( z[ ioffz + j4 - 1 ], emin );
      }
    } else {
      for ( j4 = 4 * i0; j4 <= 4 * ( n0 - 3 ); j4 += 4 ) {
        z[ ioffz + j4 - 4 ] = d + z[ ioffz + j4 - 1 ];
        temp = z[ ioffz + j4 + 1 ] / z[ ioffz + j4 - 4 ];
        d = d * temp - tau;
        dmin.setValue( Math.min( dmin.getValue(), d ) );
        z[ ioffz + j4 - 2 ] = z[ ioffz + j4 - 1 ] * temp;
        emin = Math.min( z[ ioffz + j4 - 2 ], emin );
      }
    }
    dnm2.setValue( d );
    dmin2.setValue( dmin.getValue() );
    j4 = 4 * ( n0 - 2 ) - pp;
    var j4p2 = j4 + 2 * pp - 1;
    z[ ioffz + j4 - 3 ] = dnm2.getValue() + z[ ioffz + j4p2 - 1 ];
    z[ ioffz + j4 - 1 ] = z[ ioffz + j4p2 + 1 ]
      * ( dnm2.getValue() / z[ ioffz + j4 - 3 ] ) - tau;
    dmin.setValue( Math.min( dmin.getValue(), dnm1.getValue() ) );
    dmin1.setValue( dmin.getValue() );
    j4 += 4;
    j4p2 = j4 + 2 * pp - 1;
    z[ ioffz + j4 - 3 ] = dnm1.getValue() + z[ ioffz + j4p2 - 1 ];
    z[ ioffz + j4 - 1 ] = z[ ioffz + j4p2 + 1 ]
      * ( z[ ioffz + j4p2 - 1 ] / z[ ioffz + j4 - 3 ] );
    dn.setValue( z[ ioffz + j4p2 + 1 ] 
      * ( dnm1.getValue() / z[ ioffz + j4 - 3 ] ) - tau );
    dmin.setValue( Math.min( dmin.getValue(), dn.getValue() ) );
  } else {
    if ( pp == 0 ) {
      for ( j4 = 4 * i0; j4 <= 4 * ( n0 - 3 ); j4 += 4 ) {
        z[ ioffz + j4 - 3 ] = d + z[ ioffz + j4 - 2 ];
        if ( d < 0. ) return;
        else {
          z[ ioffz + j4 - 1 ] = z[ ioffz + j4 ]
            * ( z[ ioffz + j4 - 2 ] / z[ ioffz + j4 - 3 ] );
          d = z[ ioffz + j4 ] * ( d / z[ ioffz + j4 - 3 ] ) - tau;
        }
        dmin.setValue( Math.min( dmin.getValue(), d ) );
        emin = Math.min( emin, z[ ioffz + j4 - 1 ] );
      }
    } else {
      for ( j4 = 4 * i0; j4 <= 4 * ( n0 - 3 ); j4 += 4 ) {
        z[ ioffz + j4 - 4 ] = d + z[ ioffz + j4 - 1 ];
        if ( d < 0. ) return;
        else {
          z[ ioffz + j4 - 2 ] = z[ ioffz + j4 + 1 ]
            * ( z[ ioffz + j4 - 1 ] / z[ ioffz + j4 - 4 ] );
          d = z[ ioffz + j4 + 1 ] * ( d / z[ ioffz + j4 - 4 ] ) - tau;
        }
        dmin.setValue( Math.min( dmin.getValue(), d ) );
        emin = Math.min( emin, z[ ioffz + j4 - 2 ] );
      }
    }
    dnm2.setValue( d );
    dmin2.setValue( dmin.getValue() );
    j4 = 4 * ( n0 - 2 ) - pp;
    j4p2 = j4 + 2 * pp - 1;
    z[ ioffz + j4 - 3 ] = dnm2.getValue() + z[ ioffz + j4p2 - 1 ];
    if ( dnm2.getValue() < 0. ) return
    else {
      z[ ioffz + j4 - 1 ] = z[ ioffz + j4p2 + 1 ]
        * ( z[ ioffz + j4p2 - 1 ] / z[ ioffz + j4 - 3 ] );
      dnm1.setValue( z[ ioffz + j4p2 + 1 ] ) 
        * ( dnm2.getValue() / z[ ioffz + j4 - 3 ] ) - tau;
    }
    dmin.setValue( Math.min( dmin.getValue(), dnm1.getValue() ) );
    dmin1.setValue( dmin.getValue() );
    j4 += 4;
    j4p2 = j4 + 2 * pp - 1;
    z[ ioffz + j4 - 3 ] = dnm1.getValue() + z[ ioffz + j4p2 - 1 ];
    if ( dnm1.getValue() < 0. ) return;
    else {
      z[ ioffz + j4 - 1 ] = z[ ioffz + j4p2 + 1 ]
        * ( z[ ioffz + j4p2 - 1 ] / z[ ioffz + j4 - 3 ] );
      dn.setValue( z[ ioffz + j4p2 + 1 ] 
        * ( dnm1.getValue() / z[ ioffz + j4 - 3 ] ) - tau );
    }
    dmin.setValue( Math.min( dmin.getValue(), dn.getValue() ) );
  }
  z[ ioffz + j4 + 1 ] = dn.getValue();
  z[ ioffz + 4 * n0 - pp ] = emin;
}
//*************************************************************************
LaPack0.dlasr = function( side, pivot, direct, m, n, c, s, A,
lda, ioffc, ioffs, ioffa ) {
  var info = 0;
  if ( side.charAt(0).toUpperCase() != 'L' &&
  side.charAt(0).toUpperCase() != 'R' ) {
    info = 1;
  } else if ( pivot.charAt(0).toUpperCase() != 'V' &&
  pivot.charAt(0).toUpperCase() != 'T' &&
  pivot.charAt(0).toUpperCase() != 'B' ) {
    info = 2;
  } else if ( direct.charAt(0).toUpperCase() != 'F' &&
  direct.charAt(0).toUpperCase() != 'B' ) {
    info = 3;
  } else if ( m < 0 ) info = 4;
  else if ( n < 0 ) info = 5;
  else if ( lda < Math.max( 1, m ) ) info = 9;
  if ( info != 0 ) Blas2.xerbla( 'dlasr', info );
  if ( m == 0 || n == 0 ) return;
  var i = -1;
  var j = -1;
  var ctemp = Number.POSITIVE_INFINITY;
  var stemp = Number.POSITIVE_INFINITY;
  var temp = Number.POSITIVE_INFINITY;
  if ( side.charAt(0).toUpperCase() == 'L' ) {
    if ( pivot.charAt(0).toUpperCase() == 'V' ) {
      if ( direct.charAt(0).toUpperCase() == 'F' ) {
        for ( j = 1; j <= m - 1; j ++ ) {
          ctemp = c[ ioffc + j - 1 ];
          stemp = s[ ioffs + j - 1 ];
          if ( ctemp != 1. || stemp != 0. ) {
            for ( i = 1; i <= n; i ++ ) {
              temp = A[ ioffa + j + ( i - 1 ) * lda ];
              A[ ioffa + j + ( i - 1 ) * lda ] =
                ctemp * temp
                - stemp * A[ ioffa + j - 1 + ( i - 1 ) * lda ];
              A[ ioffa + j - 1 + ( i - 1 ) * lda ] =
                stemp * temp
                + ctemp * A[ ioffa + j - 1 + ( i - 1 ) * lda ];
            }
          }
        }
      } else if ( direct.charAt(0).toUpperCase() == 'B' ) {
        for ( j = m - 1; j >= 1; j -- ) {
          ctemp = c[ ioffc + j - 1 ];
          stemp = s[ ioffs + j - 1 ];
          if ( ctemp != 1. || stemp != 0. ) {
            for ( i = 1; i <= n; i ++ ) {
              temp = A[ ioffa + j + ( i - 1 ) * lda ];
              A[ ioffa + j + ( i - 1 ) * lda ] =
                ctemp * temp
                - stemp * A[ ioffa + j - 1 + ( i - 1 ) * lda ];
              A[ ioffa + j - 1 + ( i - 1 ) * lda ] =
                stemp * temp
                + ctemp * A[ ioffa + j - 1 + ( i - 1 ) * lda ];
            }
          }
        }
      }
    } else if ( pivot.charAt(0).toUpperCase() == 'T' ) {
      if ( direct.charAt(0).toUpperCase() == 'F' ) {
        for ( j = 2; j <= m; j ++ ) {
          ctemp = c[ ioffc + j - 2 ];
          stemp = s[ ioffs + j - 2 ];
          if ( ctemp != 1. || stemp != 0. ) {
            for ( i = 1; i <= n; i ++ ) {
              temp = A[ ioffa + j - 1 + ( i - 1 ) * lda ];
              A[ ioffa + j - 1 + ( i - 1 ) * lda ] =
                ctemp * temp - stemp * A[ ioffa + ( i - 1 ) * lda ];
              A[ ioffa + ( i - 1 ) * lda ] =
                stemp * temp + ctemp * A[ ioffa + ( i - 1 ) * lda ];
            }
          }
        }
      } else if ( direct.charAt(0).toUpperCase() == 'B' ) {
        for ( j = m; j >= 2; j -- ) {
          ctemp = c[ ioffc + j - 2 ];
          stemp = s[ ioffs + j - 2 ];
          if ( ctemp != 1. || stemp != 0. ) {
            for ( i = 1; i <= n; i ++ ) {
              temp = A[ ioffa + j - 1 + ( i - 1 ) * lda ];
              A[ ioffa + j - 1 + ( i - 1 ) * lda ] =
                ctemp * temp - stemp * A[ ioffa + ( i - 1 ) * lda ];
              A[ ioffa + ( i - 1 ) * lda ] =
                stemp * temp + ctemp * A[ ioffa + ( i - 1 ) * lda ];
            }
          }
        }
      }
    } else if ( pivot.charAt(0).toUpperCase() == 'B' ) {
      if ( direct.charAt(0).toUpperCase() == 'F' ) {
        for ( j = 1; j <= m - 1; j ++ ) {
          ctemp = c[ ioffc + j - 1 ];
          stemp = s[ ioffs + j - 1 ];
          if ( ctemp != 1. || stemp != 0. ) {
            for ( i = 1; i <= n; i ++ ) {
              temp = A[ ioffa + j - 1 + ( i - 1 ) * lda ];
              A[ ioffa + j - 1 + ( i - 1 ) * lda ] =
                stemp * A[ ioffa + m - 1 + ( i - 1 ) * lda ]
                + ctemp * temp;
              A[ ioffa + m - 1 + ( i - 1 ) * lda ] =
                ctemp * A[ ioffa + m - 1 + ( i - 1 ) * lda ]
                - stemp * temp;
            }
          }
        }
      } else if ( direct.charAt(0).toUpperCase() == 'B' ) {
        for ( j = m - 1; j >= 1; j -- ) {
          ctemp = c[ ioffc + j - 1 ];
          stemp = s[ ioffs + j - 1 ];
          if ( ctemp != 1. || stemp != 0. ) {
            for ( i = 1; i <= n; i ++ ) {
              temp = A[ ioffa + j - 1 + ( i - 1 ) * lda ];
              A[ ioffa + j - 1 + ( i - 1 ) * lda ] =
                stemp * A[ ioffa + m - 1 + ( i - 1 ) * lda ]
                + ctemp * temp;
              A[ ioffa + m - 1 + ( i - 1 ) * lda ] =
                ctemp * A[ ioffa + m - 1 + ( i - 1 ) * lda ]
                - stemp * temp;
            }
          }
        }
      }
    }
  } else if ( side.charAt(0).toUpperCase() == 'R' ) {
    if ( pivot.charAt(0).toUpperCase() == 'V' ) {
      if ( direct.charAt(0).toUpperCase() == 'F' ) {
        for ( j = 1; j <= n - 1; j ++ ) {
          ctemp = c[ ioffc + j - 1 ];
          stemp = s[ ioffs + j - 1 ];
          if ( ctemp != 1. || stemp != 0. ) {
            for ( i = 1; i <= m; i ++ ) {
              temp = A[ ioffa + i - 1 + j * lda ];
              A[ ioffa + i - 1 + j * lda ] =
                ctemp * temp
                - stemp * A[ ioffa + i - 1 + ( j - 1 ) * lda ];
              A[ ioffa + i - 1 + ( j - 1 ) * lda ] =
                stemp * temp
                + ctemp * A[ ioffa + i - 1 + ( j - 1 ) * lda ];
            }
          }
        }
      } else if ( direct.charAt(0).toUpperCase() == 'B' ) {
        for ( j = n - 1; j >= 1; j -- ) {
          ctemp = c[ ioffc + j - 1 ];
          stemp = s[ ioffs + j - 1 ];
          if ( ctemp != 1. || stemp != 0. ) {
            for ( i = 1; i <= m; i ++ ) {
              temp = A[ ioffa + i - 1 + j * lda ];
              A[ ioffa + i - 1 + j * lda ] =
                ctemp * temp
                - stemp * A[ ioffa + i - 1 + ( j - 1 ) * lda ];
              A[ ioffa + i - 1 + ( j - 1 ) * lda ] =
                stemp * temp
                + ctemp * A[ ioffa + i - 1 + ( j - 1 ) * lda ];
            }
          }
        }
      }
    } else if ( pivot.charAt(0).toUpperCase() == 'T' ) {
      if ( direct.charAt(0).toUpperCase() == 'F' ) {
        for ( j = 2; j <= n; j ++ ) {
          ctemp = c[ ioffc + j - 2 ];
          stemp = s[ ioffs + j - 2 ];
          if ( ctemp != 1. || stemp != 0. ) {
            for ( i = 1; i <= m; i ++ ) {
              temp = A[ ioffa + i - 1 + ( j - 1 ) * lda ];
              A[ ioffa + i - 1 + ( j - 1 ) * lda ] =
                ctemp * temp - stemp * A[ ioffa + i - 1 ];
              A[ ioffa + i - 1 ] =
                stemp * temp + ctemp * A[ ioffa + i - 1 ];
            }
          }
        }
      } else if ( direct.charAt(0).toUpperCase() == 'B' ) {
        for ( j = n; j >= 2; j -- ) {
          ctemp = c[ ioffc + j - 2 ];
          stemp = s[ ioffs + j - 2 ];
          if ( ctemp != 1. || stemp != 0. ) {
            for ( i = 1; i <= m; i ++ ) {
              temp = A[ ioffa + i - 1 + ( j - 1 ) * lda ];
              A[ ioffa + i - 1 + ( j - 1 ) * lda ] =
                ctemp * temp - stemp * A[ ioffa + i - 1 ];
              A[ ioffa + i - 1 ] = stemp * temp
                + ctemp * A[ ioffa + i - 1 ];
            }
          }
        }
      }
    } else if ( pivot.charAt(0).toUpperCase() == 'B' ) {
      if ( direct.charAt(0).toUpperCase() == 'F' ) {
        for ( j = 1; j <= n - 1; j ++ ) {
          ctemp = c[ ioffc + j - 1 ];
          stemp = s[ ioffs + j - 1 ];
          if ( ctemp != 1. || stemp != 0. ) {
            for ( i = 1; i <= m; i ++ ) {
              temp = A[ ioffa + i - 1 + ( j - 1 ) * lda ];
              A[ ioffa + i - 1 + ( j - 1 ) * lda ] =
                stemp * A[ ioffa + i - 1 + ( n - 1 ) * lda ]
                + ctemp * temp;
              A[ ioffa + i - 1 + ( n - 1 ) * lda ] =
                ctemp * A[ ioffa + i - 1 + ( n - 1 ) * lda ]
                - stemp * temp;
            }
          }
        }
      } else if ( direct.charAt(0).toUpperCase() == 'B' ) {
        for ( j = n - 1; j >= 1; j -- ) {
          ctemp = c[ ioffc + j - 1 ];
          stemp = s[ ioffs + j - 1 ];
          if ( ctemp != 1. || stemp != 0. ) {
            for ( i = 1; i <= m; i ++ ) {
              temp = A[ ioffa + i - 1 + ( j - 1 ) * lda ];
              A[ ioffa + i - 1 + ( j - 1 ) * lda ] =
                stemp * A[ ioffa + i - 1 + ( n - 1 ) * lda ]
                + ctemp * temp;
              A[ ioffa + i - 1 + ( n - 1 ) * lda ] =
                ctemp * A[ ioffa + i - 1 + ( n - 1 ) * lda ]
                - stemp * temp;
            }
          }
        }
      }
    }
  }
}
LaPack0.zlasr = function( side, direct, pivot, m, n, c, s, A,
lda ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dlasrt = function( id, n, d, info, ioffd ) {
  info.setValue( 0 );
  var dir = -1;
  if ( id.charAt(0).toUpperCase() == 'D' ) dir = 0;
  else if ( id.charAt(0).toUpperCase() == 'I' ) dir = 1;
  if ( dir == -1 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlasrt', - info.getValue() );
    return;
  }
  if ( n <= 1 ) return;
  var stkpnt = 1;
  var stack = new Array( 64 );
  stack[ 0 ] = 1;
  stack[ 1 ] = n;
  var select = 20;
  do {
    var start = stack[ ( stkpnt - 1 ) * 2 ];
    var endd = stack[ 1 + ( stkpnt - 1 ) * 2 ];
    stkpnt --;
    if ( endd - start <= select && endd - start > 0 ) {
      if ( dir == 0 ) {
        for ( var i = start + 1; i <= endd; i ++ ) {
          for ( var j = i; j >= start + 1;  j -- ) {
            if ( d[ ioffd + j - 1 ] > d[ ioffd + j - 2 ] ) {
              var dmnmx = d[ ioffd + j - 1 ];
              d[ ioffd + j - 1 ] = d[ ioffd + j - 2 ];
              d[ ioffd + j - 2 ] = dmnmx;
            } else break;
          }
        }
      } else {
        for ( i = start + 1; i <= endd; i ++ ) {
          for ( j = i; j >= start + 1;  j -- ) {
            if ( d[ ioffd + j - 1 ] < d[ ioffd + j - 2 ] ) {
              dmnmx = d[ ioffd + j - 1 ];
              d[ ioffd + j - 1 ] = d[ ioffd + j - 2 ];
              d[ ioffd + j - 2 ] = dmnmx;
            } else break;
          }
        }
      }
    } else if ( endd - start > select ) {
      var d1 = d[ ioffd + start - 1 ];
      var d2 = d[ ioffd + endd - 1 ];
      i = ( start + endd ) / 2;
      var d3 = d[ ioffd + i - 1 ];
      if ( d1 < d2 ) {
        if ( d3 < d1 ) dmnmx = d1;
        else if ( d3 < d2 ) dmnmx = d3;
        else dmnmx = d2;
      } else {
        if ( d3 < d2 ) dmnmx = d2;
        else if ( d3 < d1 ) dmnmx = d3;
        else dmnmx = d1;
      }
      if ( dir == 0 ) {
        i = start - 1;
        j = endd + 1;
        while ( true ) {
          do {
            j --;
          } while ( d[ ioffd + j - 1 ] < dmnmx );
          do {
            i ++;
          } while ( d[ ioffd + i - 1 ] > dmnmx );
          if ( i < j ) {
            var tmp = d[ ioffd + i - 1 ];
            d[ ioffd + i - 1 ] = d[ ioffd + j - 1 ];
            d[ ioffd + j - 1 ] = tmp;
          } else break;
        }
        if ( j - start > endd - j - 1 ) {
          stkpnt ++;
          stack[ ( stkpnt - 1 ) * 2 ] = start;
          stack[ 1 + ( stkpnt - 1 ) * 2 ] = j;
          stkpnt ++;
          stack[ ( stkpnt - 1 ) * 2 ] = j + 1;
          stack[ 1 + ( stkpnt - 1 ) * 2 ] = endd;
        } else {
          stkpnt ++;
          stack[ ( stkpnt - 1 ) * 2 ] = j + 1;
          stack[ 1 + ( stkpnt - 1 ) * 2 ] = endd;
          stkpnt ++;
          stack[ ( stkpnt - 1 ) * 2 ] = start;
          stack[ 1 + ( stkpnt - 1 ) * 2 ] = j;
        }
      } else {
        i = start - 1;
        j = endd + 1;
        while ( true ) {
          do {
            j --;
          } while ( d[ ioffd + j - 1 ] > dmnmx );
          do {
            i ++;
          } while ( d[ ioffd + i - 1 ] < dmnmx );
          if ( i < j ) {
            tmp = d[ ioffd + i - 1 ];
            d[ ioffd + i - 1 ] = d[ ioffd + j - 1 ];
            d[ ioffd + j - 1 ] = tmp;
          } else break;
        }
        if ( j - start > endd - j - 1 ) {
          stkpnt ++;
          stack[ ( stkpnt - 1 ) * 2 ] = start;
          stack[ 1 + ( stkpnt - 1 ) * 2 ] = j;
          stkpnt ++;
          stack[ ( stkpnt - 1 ) * 2 ] = j + 1;
          stack[ 1 + ( stkpnt - 1 ) * 2 ] = endd;
        } else {
          stkpnt ++;
          stack[ ( stkpnt - 1 ) * 2 ] = j + 1;
          stack[ 1 + ( stkpnt - 1 ) * 2 ] = endd;
          stkpnt ++;
          stack[ ( stkpnt - 1 ) * 2 ] = start;
          stack[ 1 + ( stkpnt - 1 ) * 2 ] = j;
        }
      }
    }
  } while ( stkpnt > 0 );
}
//*************************************************************************
LaPack0.dlassq = function( n, x, incx, scaleReference,
sumsqReference, ioffx ) {
  if ( n > 0 ) {
    for ( var ix = 1; ix <= 1 + ( n - 1 ) * incx; ix += incx ) {
      if ( x[ ioffx + ix - 1 ] != 0. ) {
        var absxi = Math.abs( x[ ioffx + ix - 1 ] );
        if ( scaleReference.getValue() < absxi ) {
          sumsqReference.setValue( 1. + sumsqReference.getValue()
            * Math.pow( scaleReference.getValue() / absxi , 2 ) );
          scaleReference.setValue( absxi );
        } else {
          sumsqReference.setValue( sumsqReference.getValue()
            + Math.pow( absxi / scaleReference.getValue() , 2 ) );
        }
      }
    }
  }
}
LaPack0.zlassq = function( n, x, incx, scaleReference,
sumsqReference, ioffx ) {
  var zero = 0.;
  if ( n > 0 ) {
    for ( var ix = 1; ix <= 1 + ( n - 1 ) * incx; ix += incx ) {
      if ( x[ ioffx + ix - 1 ].real != zero ) {
        var temp1 = Math.abs( x[ ioffx + ix - 1 ].real );
        if ( scaleReference.getValue() < temp1 ) {
          sumsqReference.setValue( 1. + sumsqReference.getValue()
            * Math.pow( scaleReference.getValue() / temp1 , 2 ) );
          scaleReference.setValue( temp1 );
        } else {
          sumsqReference.setValue( sumsqReference.getValue()
            + Math.pow( temp1 / scaleReference.getValue() , 2 ) );
        }
      }
      if ( x[ ioffx + ix - 1 ].imag != zero) {
        temp1 = Math.abs( x[ ioffx + ix - 1 ].imag );
        if (scaleReference.getValue() < temp1 ) {
          sumsqReference.setValue( 1. + sumsqReference.getValue()
            * Math.pow( scaleReference.getValue() / temp1 , 2 ) );
          scaleReference.setValue( temp1 );
        } else {
          sumsqReference.setValue( sumsqReference.getValue()
            + Math.pow( temp1 / scaleReference.getValue() , 2 ) );
        }
      }
    }
  }
}
//*************************************************************************
LaPack0.dlaswp = function( n, A, lda, k1, k2, ipiv, incx, ioffa,
ioffipiv ) {
  if ( incx > 0 ) {
    var ix0 = k1;
    var i1 = k1;
    var i2 = k2;
    var inc = 1;
  } else if ( incx < 0 ) {
    ix0 = 1 + ( 1 - k2 ) * incx;
    i1 = k2;
    i2 = k1;
    inc = -1;
  } else return;
  var n32 = Math.floor( n / 32 ) * 32;
  if ( n32 != 0 ) {
    for ( var j = 1; j <= n32; j += 32 ) {
      var ix = ix0;
      for ( var i = i1;
      ( inc == 1 && i <= i2 ) || ( inc == -1 && i >= i2 ); i += inc ) {
        var ip = ipiv[ ioffipiv + ix - 1 ];
        if ( ip != i ) {
          for ( var k = j; k <= j + 31; k ++ ) {
            var temp = A[ ioffa + i - 1 + ( k - 1 ) * lda ];
            A[ ioffa + i - 1 + ( k - 1 ) * lda ]
              = A[ ioffa + ip - 1 + ( k - 1 ) * lda ];
            A[ ioffa + ip - 1 + ( k - 1 ) * lda ] = temp;
          }
        }
        ix += incx;
      }
    }
  }
  if ( n32 != n ) {
    n32 ++;
    ix = ix0;
    for ( i = i1; ( inc == 1 && i <= i2 ) || ( inc == -1 && i >= i2 );
    i += inc ) {
      ip = ipiv[ ioffipiv + ix - 1];
      if ( ip != i ) {
        for ( k = n32; k <= n; k ++ ) {
          temp = A[ ioffa + i - 1 + ( k - 1 ) * lda ];
          A[ ioffa + i - 1 + ( k - 1 ) * lda ]
            = A[ ioffa + ip - 1 + ( k - 1 ) * lda ];
          A[ ioffa + ip - 1 + ( k - 1 ) * lda ] = temp;
        }
      }
      ix += incx;
    }
  }
/* old version:
  if ( incx == 0 ) return;
  var ix = ( incx > 0 ? k1 : 1 + ( 1 - k2 ) * incx ); 
  if ( incx == 1 ) {
    for ( var i = k1; i <= k2; i ++ ) {
      var ip = ipiv[ ioffipiv + i - 1 ];
      if ( ip != i ) {
        Blas1.dswap( n, A, lda, A, lda, ioffa + i - 1,
          ioffa + ip - 1 );
      }
    }
  } else if ( incx > 1 ) {
    for ( i = k1; i <= k2; i ++ ) {
      ip = ipiv[ ioffipiv + ix - 1 ];
      if ( ip != i ) {
        Blas1.dswap( n, A, lda, A, lda, ioffa + i - 1,
          ioffa + ip - 1 );
      }
      ix += incx;
    }
  } else if ( incx < 0 ) {
    for ( i = k2; i >= k1; i -- ) {
      ip = ipiv[ ioffipiv + ix - 1 ];
      if ( ip != i ) {
        Blas1.dswap( n, A, lda, A, lda, ioffa + i - 1,
          ioffa + ip - 1 );
      }
      ix += incx;
    }
  }
*/
}
LaPack0.zlaswp = function( n, A, lda, k1, k2, ipiv, incx ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dlasyf = function( uplo, n, nb, kb, A, lda, ipiv, W,
ldw, info, ioffa, ioffipiv, ioffw ) {
  info.setValue( 0 );
  var alpha = ( 1. + Math.sqrt( 17. ) ) / 8.;
  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    var k = n;
    while ( true ) { // 10
      var kw = nb + k - n;
      if ( ( k <= n - nb + 1 && nb < n ) || k < 1 ) break;
      Blas1.dcopy( k, A, 1, W, 1, ioffa + ( k - 1 ) * lda,
        ioffw + ( kw - 1 ) * ldw );
      if ( k < n ) {
        Blas2.dgemv( 'No transpose', k, n - k, -1., A, lda, W, ldw,
          1., W, 1, ioffa + k * lda, ioffw + k - 1 + kw * ldw,
          ioffw + ( kw - 1 ) * ldw );
      }
      var kstep = 1;
      var absakk =
        Math.abs( W[ ioffw + k - 1 + ( kw - 1 ) * ldw ] );
      if ( k > 1 ) {
        var imax =
          Blas1.idamax( k - 1, W, 1, ioffw + ( kw - 1 ) * ldw );
        var colmax =
          Math.abs( W[ ioffw + imax - 1 + ( kw - 1 ) * ldw ] );
      } else colmax = 0.;
      if ( Math.max( absakk, colmax ) == 0. ) {
        if ( info.getValue() == 0 ) info.setValue( k );
        var kp = k;
      } else {
        if ( absakk >= alpha * colmax ) kp = k;
        else {
          Blas1.dcopy( imax, A, 1, W, 1, ioffa + ( imax - 1 ) * lda,
            ioffw + ( kw - 2 ) * ldw );
          Blas1.dcopy( k - imax, A, lda, W, 1,
            ioffa + imax - 1 + imax * lda,
            ioffw + imax + ( kw - 2 ) * ldw );
          if ( k < n ) {
            Blas2.dgemv( 'No transpose', k, n - k, -1., A, lda, W,
              ldw, 1., W, 1, ioffa + k * lda,
              ioffw + imax - 1 + kw * ldw, ioffw + ( kw - 2 ) * ldw );
          }
          var jmax = imax
           + Blas1.idamax( k - imax, W, 1,
           ioffw + imax + ( kw - 2 ) * ldw );
          var rowmax =
            Math.abs( W[ ioffw + jmax - 1 + ( kw - 2 ) * ldw ] );
          if ( imax > 1 ) {
            jmax =
              Blas1.idamax( imax - 1, W, 1, ioffw + ( kw - 2 ) * ldw );
            rowmax = Math.max( rowmax,
              Math.abs( W[ ioffw + jmax - 1 + ( kw - 2 ) * ldw ] ) );
          }
          if ( absakk >= alpha * colmax * ( colmax / rowmax ) ) {
            kp = k;
          } else if ( Math.abs(
          W[ ioffw + imax - 1 + ( kw - 2 ) * ldw ] ) >=
          alpha * rowmax ) {
            kp = imax;
            Blas1.dcopy( k, W, 1, W, 1, ioffw + ( kw - 2 ) * ldw,
              ioffw + ( kw - 1 ) * ldw );
          } else {
            kp = imax;
            kstep = 2;
          }
        }
        var kk = k - kstep + 1;
        var kkw = nb + kk - n;
        if ( kp !== kk ) {
          A[ ioffa + k - 1 + ( k - 1 ) * lda ] =
            A[ ioffa + kk - 1 + ( k -1 ) * lda ];
          Blas1.dcopy( k - 1 - kp, A, 1, A, lda,
            ioffa + kp + ( kk - 1 ) * lda, ioffa + kp - 1 + kp * lda );
          Blas1.dcopy( kp, A, 1, A, 1, ioffa + ( kk - 1 ) * lda,
            ioffa + ( kp - 1 ) * lda );
          Blas1.dswap( n - kk + 1, A, lda, A, lda,
            ioffa + kk - 1 + ( kk - 1 ) * lda,
            ioffa + kp - 1 + ( kk - 1 ) * lda);
          Blas1.dswap( n - kk + 1, W, ldw, W, ldw,
            ioffw + kk - 1 + ( kkw - 1 ) * ldw,
            ioffw + kp - 1 + ( kkw - 1 ) * ldw);
        }
        if ( kstep == 1 ) {
          Blas1.dcopy( k, W, 1, A, 1, ioffw + ( kw - 1 ) * ldw,
            ioffa + ( k - 1 ) * lda );
          var r1 = 1. / A[ ioffa + k - 1 + ( k - 1 ) * lda ];
          Blas1.dscal( k - 1 , r1, A, 1, ioffa + ( k - 1 ) * lda );
        } else {
          if ( k > 2 ) {
            var d21 = W[ ioffw + k - 2 + ( kw - 1 ) * ldw ];
            var d11 =
              W[ ioffw + k - 1 + ( kw - 1 ) * ldw ] / d21;
            var d22 =
              W[ ioffw + k - 2 + ( kw - 2 ) * ldw ] / d21;
            var t = 1. / ( d11 * d22 - 1. );
            d21 = t / d21;
            for ( var j = 1; j <= k - 2; j ++ ) {
              A[ ioffa + j - 1 + ( k - 2 ) * lda ] = d21
                * ( d11 * W[ ioffw + j - 1 + ( kw - 2 ) * ldw ]
                  - W[ ioffw + j - 1 + ( kw - 1 ) * ldw ] );
              A[ ioffa + j - 1 + ( k - 1 ) * lda ] = d21
                * ( d22 * W[ ioffw + j - 1 + ( kw - 1 ) * ldw ]
                  - W[ ioffw + j- 1 + ( kw - 2 ) * ldw ] );
            }
          }
          A[ ioffa + k - 2 + ( k - 2 ) * lda ] =
            W[ ioffw + k - 2 + ( kw - 2 ) * ldw ];
          A[ ioffa + k - 2 + ( k - 1 ) * lda ] =
            W[ ioffw + k - 2 + ( kw - 1 ) * ldw ];
          A[ ioffa + k - 1 + ( k - 1 ) * lda ] =
            W[ ioffw + k - 1 + ( kw - 1 ) * ldw ];
        }
      }
      if ( kstep == 1 )  ipiv[ ioffipiv + k - 1 ] = kp;
      else {
        ipiv[ ioffipiv + k - 1 ] = - kp;
        ipiv[ ioffipiv + k - 2 ] = - kp;
      }
      k -= kstep;
    } // 30
    for ( j = ( ( k - 1 ) / nb ) * nb + 1; j >= 1; j -= nb ) {
      var jb = Math.min( nb, k - j + 1 );
      for ( var jj = j; jj <= j + jb - 1; jj ++ ) {
        Blas2.dgemv( 'No transpose', jj - j + 1, n - k, -1., A, lda,
          W, ldw, 1., A, 1, ioffa + j - 1 + k * lda,
          ioffw + jj - 1 + kw * ldw, ioffa + j- 1 + ( jj - 1 ) * lda );
      }
      Blas3.dgemm( 'No transpose', 'Transpose', j - 1, jb, n - k, -1.,
        A, lda, W, ldw, 1., A, lda, ioffa + k * lda,
        ioffw + j - 1 + kw * ldw, ioffa + ( j - 1 ) * lda );
    }
    j = k + 1;
    do {
      jj = j;
      var jp = ipiv[ ioffipiv + j - 1 ];
      if ( jp < 0 ) {
        jp = - jp;
        j ++;
      }
      j ++;
      if ( jp != jj && j <= n ) {
        Blas1.dswap( n - j + 1, A, lda, A, lda,
          ioffa + jp - 1 + ( j - 1 ) * lda,
          ioffa + jj - 1 + ( j - 1 ) * lda );
      }
    } while ( j <= n );
    kb.setValue( n - k );
  } else {
    k = 1;
    while ( true ) {
      if ( ( k >= nb && nb < n ) || k > n ) break;
      Blas1.dcopy( n - k + 1, A, 1, W, 1,
        ioffa + k - 1 + ( k - 1 ) * lda,
        ioffw + k - 1 + ( k - 1 ) * ldw );
      Blas2.dgemv( 'No transpose', n - k + 1, k - 1, -1., A, lda, W,
        ldw, 1., W, 1, ioffa + k - 1, ioffw + k - 1,
        ioffw + k - 1 + ( k - 1 ) * ldw );
      kstep = 1;
      absakk = Math.abs( W[ ioffw + k - 1 + ( k - 1 ) * ldw ] );
      if ( k < n ) {
        imax = k + Blas1.idamax( n - k, W, 1,
          ioffw + k + ( k - 1 ) * ldw );
        colmax = Math.abs( W[ ioffw + imax - 1 + ( k - 1 ) * ldw ] );
      } else {
        colmax = 0.;
      }
      if ( Math.max( absakk, colmax ) == 0. ) {
        if ( info.getValue() == 0 ) info.setValue( k );
        kp = k;
      } else {
        if ( absakk >= alpha * colmax )  kp = k;
        else {
          Blas1.dcopy( imax - k, A, lda, W, 1,
            ioffa + imax - 1 + ( k - 1 ) * lda,
            ioffw + k - 1 + k * ldw );
          Blas1.dcopy( n - imax + 1, A, 1, W, 1,
            ioffa + imax - 1 + ( imax - 1 ) * lda,
            ioffw + imax - 1 + k * ldw );
          Blas2.dgemv( 'No transpose', n - k + 1, k - 1, -1., A, lda,
            W, ldw, 1., W, 1, ioffa + k - 1, ioffw + imax - 1,
            ioffw + k - 1 + k * ldw );
          jmax = k - 1 + Blas1.idamax( imax - k, W, 1,
            ioffw + k - 1 + k * ldw );
          rowmax = Math.abs( W[ ioffw + jmax - 1 + k * ldw ] );
          if ( imax < n ) {
            jmax = imax
              + Blas1.idamax( n - imax, W, 1, ioffw + imax + k * ldw );
            rowmax = Math.max( rowmax,
              Math.abs( W[ ioffw + jmax - 1 + k * ldw ] ) );
          }
          if ( absakk >= alpha * colmax * ( colmax / rowmax ) ) {
            kp = k;
          } else if ( Math.abs( W[ ioffw + imax - 1 + k * ldw ] ) >=
          alpha * rowmax ) {
            kp = imax;
            Blas1.dcopy( n - k + 1, W, 1, W, 1,
              ioffw + k - 1 + k * ldw,
              ioffw + k - 1 + ( k - 1 ) * ldw );
          } else {
            kp = imax;
            kstep = 2;
          }
        }
        kk = k + kstep - 1;
        if ( kp != kk ) {
          A[ ioffa + kp - 1 + ( k - 1 ) * lda ] =
            A[ ioffa + kk - 1 + ( k - 1 ) * lda ];
          Blas1.dcopy( kp - k - 1, A, 1, A, lda,
            ioffa + k + ( kk - 1 ) * lda, ioffa + kp - 1 + k * lda );
          Blas1.dcopy( n - kp + 1, A, 1, A, 1,
            ioffa + kp - 1 * ( kk - 1 ) * lda,
            ioffa + kp - 1 + ( kp - 1 ) * lda );
          Blas1.dswap( kk, A, lda, A, lda, ioffa + kk - 1,
            ioffa + kp - 1 );
          Blas1.dswap( kk, W, ldw, W, ldw, ioffw + kk - 1,
            ioffw + kp - 1 );
        }
        if ( kstep == 1 ) {
          Blas1.dcopy( n - k + 1, W, 1, A, 1,
            ioffw + k - 1 + ( k - 1 ) * ldw,
            ioffa + k - 1 + ( k - 1 ) * lda );
          if ( k < n ) {
            r1 = 1. / A[ ioffa + k - 1 + ( k - 1 ) * lda ];
            Blas1.dscal( n - k, r1, A, 1,
              ioffa + k + ( k - 1 ) * lda );
          }
        } else {
          if ( k < n - 1 ) {
            d21 = W[ ioffw + k + ( k - 1 ) * ldw ];
            d11 = W[ ioffw + k + k * ldw ] / d21;
            d22 = W[ ioffw + k - 1 + ( k - 1 ) * ldw ] / d21;
            t = 1. / ( d11 * d22 - 1. );
            d21 = t / d21;
            for ( j = k + 2; j <= n; j ++ ) {
              A[ ioffa + j - 1 + ( k - 1 ) * lda ] = d21
                * ( d11 * W[ ioffw + j - 1 + ( k - 1 ) * ldw ]
                  - W[ ioffw + j - 1 + k * ldw ] );
              A[ ioffa + j - 1 + k * lda ] = d21
                * ( d22 * W[ ioffw + j - 1 + k * ldw ]
                  - W[ ioffw + j - 1 + ( k - 1 ) * ldw ] );
            }
          }
          A[ ioffa + k - 1 + ( k - 1 ) * lda ] =
            W[ ioffw + k - 1 + ( k - 1 ) * ldw ];
          A[ ioffa + k + ( k - 1 ) * lda ] =
            W[ ioffw + k + ( k - 1 ) * ldw ];
          A[ ioffa + k + k * lda ] = W[ ioffw + k + k * ldw ];
          }
      }
      if ( kstep == 1 ) ipiv[ ioffipiv + k - 1 ] = kp;
      else {
        ipiv[ ioffipiv + k - 1 ] = - kp;
        ipiv[ ioffipiv + k ] = - kp;
      }
      k += kstep;
    }
    for ( j = k; j <= n; j += nb ) {
      jb = Math.min( nb, n - j + 1 );
      for ( jj = j; jj <= j + jb - 1; jj ++ ) {
        Blas2.dgemv( 'No transpose', j + jb - jj, k - 1, -1., A, lda,
          W, ldw, 1., A, 1, ioffa + jj - 1, ioffw + jj - 1,
          ioffa + jj - 1 + ( jj - 1 ) * lda );
      }
      if ( j + jb <= n ) {
        Blas3.dgemm( 'No transpose', 'Transpose', n - j - jb + 1, jb,
          k - 1, -1., A, lda, W, ldw, 1., A, lda,
          ioffa + j + jb - 1, ioffw + j - 1,
          ioffa + j + jb - 1 + ( j - 1 ) * lda );
      }
    }
    j = k - 1;
    do {
      jj = j;
      jp = ipiv[ ioffipiv + j - 1 ];
      if ( jp < 0 ) {
        jp = - jp;
        j --;
      }
      j --;
      if ( jp != jj && j >= 1 ) {
        Blas1.dswap( j, A, lda, A, lda, ioffa + jp - 1,
          ioffa + jj - 1 );
      }
    } while ( j >= 1 );
    kb.setValue( k - 1 );
  }
}
LaPack0.zlasyf = function( uplo, n, nb, kb, A, lda, ipiv, W, ldw,
info ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
/* deprecated use dormrz
LaPack0.dlatzm = function( side, m, n, v, incv, tau, C1, C2,
ldc, work, ioffv, ioffc1, ioffc2, ioffwork ) {
  if ( ( Math.min( m, n ) == 0 ) || ( tau == 0. ) ) return;
  if ( side.charAt(0).toUpperCase() == 'L' ) {
    Blas1.dcopy( n, C1, ldc, work, 1, ioffc1, ioffwork );
    Blas2.dgemv( 'Transpose', m - 1, n, 1., C2, ldc, v, incv, 1.,
      work, 1, ioffc2, ioffv, ioffwork );
    Blas1.daxpy( n, -tau, work, 1, C1, ldc, ioffwork, ioffc1 );
    Blas2.dger( m - 1, n, -tau, v, incv, work, 1, C2, ldc, ioffv,
      ioffwork, ioffc2 );
  } else if ( side.charAt(0).toUpperCase() == 'R' ) {
    Blas1.dcopy( m, C1, 1, work, 1, ioffc1, ioffwork );
    Blas2.dgemv( 'No transpose', m, n - 1, 1., C2, ldc, v, incv, 1.,
      work, 1, ioffc2, ioffv, ioffwork );
    Blas1.daxpy( m, -tau, work, 1, C1, 1, ioffwork, ioffc1 );
    Blas2.dger( m, n - 1, -tau, work, 1, v, incv, C2, ldc, ioffwork,
      ioffv, ioffc2 );
  }
}
LaPack0.zlatzm = function( side, m, n, v, incv, tau, C1,
C2, ldc, work ) {
  throw new Error("not programmed: complex matrix");
}
*/
//*************************************************************************
LaPack0.dlauu2 = function( uplo, n, A, lda, info, ioffa ) {
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -4 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlauu2', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  if ( upper ) {
    for ( var i = 1; i <= n; i ++ ) {
      var aii = A[ ioffa + i - 1 + ( i - 1 ) * lda ];
      if ( i < n ) {
        A[ ioffa + i - 1 + ( i - 1 ) * lda ] =
          Blas1.ddot( n - i + 1, A, lda, A, lda,
            ioffa + i - 1 + ( i - 1 ) * lda ,
            ioffa + i - 1 + ( i - 1 ) * lda );
          Blas2.dgemv( 'No transpose', i - 1, n - i, 1., A, lda, A,
            lda, aii, A, 1, ioffa + i * lda, ioffa + i - 1 + i * lda,
            ioffa + ( i - 1 ) * lda );
      } else {
        Blas1.dscal( i, aii, A, 1, ioffa + ( i - 1 ) * lda );
      }
    }
  } else {
    for ( i = 1; i <= n; i ++ ) {
      aii = A[ ioffa + i - 1 + ( i - 1 ) * lda ];
      if ( i < n ) {
        A[ ioffa + i - 1 + ( i - 1 ) * lda ] =
          Blas1.ddot( n - i + 1, A, 1, A, 1,
            ioffa + i - 1 + ( i - 1 ) * lda,
            ioffa + i - 1 + ( i - 1 ) * lda );
          Blas2.dgemv( 'Transpose', n - i, i - 1, 1., A, lda, A,
            1, aii, A, lda, ioffa + i, ioffa + i + ( i - 1 ) * lda,
            ioffa + i - 1 );
      } else {
        Blas1.dscal( i, aii, A, lda, ioffa + i - 1 );
      }
    }
  }
}
LaPack0.zlauu2 = function( uplo, n, A, lda, info ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dpbequ = function( uplo, n, kd, AB, ldab, S,
scondReference, amaxReference, info ) {
  throw new Error("not programmed: band matrix");
}
LaPack0.zpbequ = function( uplo, n, kd, AB, ldab, S,
scondReference, amaxReference, info ) {
  throw new Error("not programmed: complex band matrix");
}
//*************************************************************************
LaPack0.dpbstf = function( uplo, n, kd, AB, ldab, info ) {
  throw new Error("not programmed: band matrix");
}
LaPack0.zpbstf = function( uplo, n, kd, AB, ldab, info ) {
  throw new Error("not programmed: complex band matrix");
}
//*************************************************************************
LaPack0.dpbtf2 = function( uplo, n, kd, AB, ldab, info ) {
  throw new Error("not programmed: band matrix");
}
LaPack0.zpbtf2 = function( uplo, n, kd, AB, ldab, info ) {
  throw new Error("not programmed: complex band matrix");
}
//*************************************************************************
LaPack0.dpbtrs = function( uplo, n, kd, nrhs, AB, ldab, B, ldb,
info ) {
  throw new Error("not programmed: band matrix");
}
LaPack0.zpbtrs = function( uplo, n, kd, nrhs, AB, ldab, B, ldb,
info ) {
  throw new Error("not programmed: complex band matrix");
}
//*************************************************************************
LaPack0.dpoequ = function( n, A, lda, s, scondReference,
amaxReference, info, ioffa, ioffs ) {
  info.setValue( 0 );
  if ( n < 0 ) info.setValue( -1 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -3 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dpoequ', - info.getValue() );
    return;
  }
  if ( n == 0 ) {
    scondReference.setValue( 1. );
    amaxReference.setValue( 0. );
    return;
  }
  s[ ioffs ] = A[ ioffa ];
  var smin = s[ ioffs ];
  amaxReference.setValue( s[ ioffs ] );
  for ( var i = 2; i <= n; i ++ ) {
    s[ ioffs + i - 1 ] = A[ ioffa + i - 1 + ( i - 1 ) * lda ];
    smin = Math.min( smin, s[ ioffs + i - 1 ] );
    amaxReference.setValue(
      Math.max( amaxReference.getValue(), s[ ioffs + i - 1 ] ) );
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
      s[ ioffs + i - 1 ] = 1. / Math.sqrt( s[ ioffs + i - 1 ] );
    }
    scondReference.setValue(
      Math.sqrt( smin ) / Math.sqrt( amaxReference.getValue() ) );
  }
}
LaPack0.zpoequ = function( n, A, lda, s, scondReference,
amaxReference, info ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dpotf2 = function( uplo, n, A, lda, info, ioffa ) {
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -4 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dpotf2', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  if ( upper ) {
    for ( var j = 1; j <= n; j ++ ) {
      var ajj = A[ ioffa + j - 1 + ( j - 1 ) * lda ]
        - Blas1.ddot( j - 1, A, 1, A, 1, ioffa + ( j - 1 ) * lda,
                      ioffa + ( j - 1 ) * lda );
      if ( ajj <= 0. ) {
        A[ ioffa + j - 1 + ( j - 1 ) * lda ] = ajj;
        info.setValue( j );
        return;
      }
      ajj = Math.sqrt( ajj );
      A[ ioffa + j - 1 + ( j - 1 ) * lda ] = ajj;
      if ( j < n ) {
        Blas2.dgemv( 'Transpose', j - 1, n - j, -1., A, lda, A, 1, 1.,
          A, lda, ioffa + j * lda, ioffa + ( j - 1 ) * lda,
          ioffa + j - 1 + j * lda );
        Blas1.dscal( n - j, 1. / ajj, A, lda,
          ioffa + j - 1 + j * lda );
      }
    }
  } else {
    for ( j = 1; j <= n; j ++ ) {
      ajj = A[ ioffa + j - 1 + ( j - 1 ) * lda ]
        - Blas1.ddot( j - 1, A, lda, A, lda, ioffa + j - 1,
        ioffa + j - 1 );
      if ( ajj <= 0. ) { 
        A[ ioffa + j - 1 + ( j - 1 ) * lda ] = ajj;
        info.setValue( j );
        return;
      }
      ajj = Math.sqrt( ajj );
      A[ ioffa + j - 1 + ( j - 1 ) * lda ] = ajj;
      if ( j < n ) { 
        Blas2.dgemv( 'No transpose', n - j, j - 1, -1., A, lda, A, lda,
          1., A, 1, ioffa + j, ioffa + j - 1,
          ioffa + j + ( j - 1 ) * lda );
        Blas1.dscal( n - j, 1. / ajj, A, 1,
          ioffa + j + ( j - 1 ) * lda );
      }
    }
  }
}
LaPack0.zpotf2 = function( uplo, n, A, lda, info ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dpotrs = function( uplo, n, nrhs, A, lda, B, ldb, info,
ioffa, ioffb) {
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( nrhs < 0 ) info.setValue( -3 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -5 );
  else if ( ldb < Math.max( 1, n ) ) info.setValue( -7 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dpotrs', - info.getValue() );
    return;
  }
  if ( n == 0 || nrhs == 0 ) return;
  if ( upper ) {
    Blas3.dtrsm( 'Left', 'Upper', 'Transpose', 'Non-unit', n, nrhs,
      1., A, lda, B, ldb, ioffa, ioffb );
    Blas3.dtrsm( 'Left', 'Upper', 'No transpose', 'Non-unit', n, nrhs,
      1., A, lda, B, ldb, ioffa, ioffb );
  } else {
    Blas3.dtrsm( 'Left', 'Lower', 'No transpose', 'Non-unit', n, nrhs,
      1., A, lda, B, ldb, ioffa, ioffb );
    Blas3.dtrsm( 'Left', 'Lower', 'Transpose', 'Non-unit', n, nrhs,
      1., A, lda, B, ldb, ioffa, ioffb );
  }
}
LaPack0.zpotrs = function( uplo, n, nrhs, A, lda, B, ldb, info )
{
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dppequ = function( uplo, n, AP, s, scondReference,
amaxReference, info ) {
  throw new Error("not programmed: packed matrix");
}
LaPack0.zppequ = function( uplo, n, AP, s, scondReference,
amaxReference, info ) {
  throw new Error("not programmed: complex packed matrix");
}
//*************************************************************************
LaPack0.dpptrf = function( uplo, n, AP, info ) {
  throw new Error("not programmed: packed matrix");
}
LaPack0.zpptrf = function( uplo, n, AP, info ) {
  throw new Error("not programmed: complex packed matrix");
}
//*************************************************************************
LaPack0.dpptrs = function( uplo, n, nrhs, AP, B, ldb, info ) {
  throw new Error("not programmed: packed matrix");
}
LaPack0.zpptrs = function( uplo, n, nrhs, AP, B, ldb, info ) {
  throw new Error("not programmed: complex packed matrix");
}
//*************************************************************************
LaPack0.dptcon = function( n, d, e, anorm, rcondReference, work,
info, ioffd, ioffe, ioffwork ) {
  info.setValue( 0 );
  if ( n < 0 ) info.setValue( -1 );
  else if ( anorm < 0. ) info.setValue( -4 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dptcon', - info.getValue() );
    return;
  }
  rcondReference.setValue( 0. );
  if ( n == 0 ) {
    rcondReference.setValue( 1. );
    return;
  } else if ( anorm == 0. ) return;
  for ( var i = 1; i <= n; i ++ ) {
    if ( d[ ioffd + i ] <= 0. ) return;
  }
  work[ ioffwork ] = 1.;
  for ( i = 2; i <= n; i ++ ) {
    work[ ioffwork + i - 1 ] =
      1. + work[ ioffwork + i - 2 ] * Math.abs( e[ ioffe + i - 2 ] ); 
  }
  work[ ioffwork + n - 1 ] /= d[ ioffd + n - 1 ];
  for ( i = n - 1; i >= 1; i -- ) {
    work[ ioffwork + i - 1 ] =
      work[ ioffwork + i - 1 ] / d [ ioffd + i - 1 ]
      + work[ ioffwork + i ] * Math.abs( e[ ioffe + i - 1 ] );
  }
  var ix = Blas1.idamax( n, work, 1, ioffwork );
  var ainvnm = Math.abs( work[ ioffwork + ix - 1 ] );
  if ( ainvnm != 0. ) rcondReference.setValue( ( 1. / ainvnm ) / anorm );
}
LaPack0.zptcon = function( n, d, e, anorm, rcondReference, work,
info ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dpttrf = function( n, d, e, info, ioffd, ioffe ) {
  info.setValue( 0 );
  if ( n < 0 ) {
    info.setValue( -1 );
    Blas2.xerbla( 'dpttrf', - info.getValue() );
    return;
  }
  if ( n == 0. ) return;
  var i4 = ( n - 1 ) % 4;
  for ( var i = 1; i <= i4; i ++ ) {
    if ( d[ ioffd + i - 1 ] <= 0. ) {
      info.setValue( i );
      return;
    }
    var ei = e[ ioffe + i - 1 ];
    e[ ioffe + i - 1 ] = ei / d[ ioffd + i - 1 ];
    d[ ioffd + i ] -= e[ ioffe + i - 1 ] * ei;
  }
  for ( i = i4 + 1; i <= n - 4; i += 4 ) {
    if ( d[ ioffd + i - 1 ] <= 0. ) {
      info.setValue( i );
      return;
    }
    ei = e[ ioffe + i - 1 ];
    e[ ioffe + i - 1 ] /= d[ ioffd + i - 1 ];
    d[ ioffd + i ] -= e[ ioffe + i - 1 ] * ei;
    if ( d[ ioffd + i ] <= 0. ) {
      info.setValue( i + 1 );
      return;
    }
    ei = e[ ioffe + i ];
    e[ ioffe + i ] = ei / d[ ioffd + i ];
    d[ ioffd + i + 1 ] -= e[ ioffe + i ] * ei;
    if ( d[ ioffd + i + 1 ] <= 0. ) {
      info.setValue( i + 2 );
      return;
    }
    ei = e[ ioffe + i + 1 ];
    e[ ioffe + i + 1 ] = ei / d[ ioffd + i + 1 ];
    d[ ioffd + i + 2 ] -= e[ ioffe + i + 1 ] * ei;
    if ( d[ ioffd + i + 2 ] <= 0. ) {
      info.setValue( i + 3 );
      return;
    }
    ei = e[ ioffe + i + 2 ];
    e[ ioffe + i + 2 ] = ei / d[ ioffd + i + 2 ];
    d[ ioffd + i + 3 ] -= e[ ioffe + i + 2 ] * ei;
  }
  if ( d[ ioffd + n - 1 ] <= 0. ) info.setValue( n );
}
LaPack0.zpttrf = function( n, d, e, info ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dptts2 = function( n, nrhs, d, e, B, ldb, ioffd, ioffe,
ioffb) {
  if ( n <= 1 ) {
    if ( n == 1 ) Blas1.dscal( nrhs, 1. / d[ ioffd ], B, ldb, ioffb );
    return;
  }
  for ( var j = 1; j <= nrhs; j ++ ) {
    for ( var i = 2; i <= n; i ++ ) {
      B[ ioffb + i - 1 + ( j - 1 ) * ldb ] -=
        B[ ioffb + i - 2 + ( j - 1 ) * ldb ] * e[ ioffe + i - 2 ];
    }
    B[ ioffb + n - 1 + ( j - 1 ) * ldb ] /= d[ ioffd + n - 1 ];
    for ( i = n - 1; i >= 1; i -- ) {
      B[ ioffb + i - 1 + ( j - 1 ) * ldb ] =
        B[ ioffb + i - 1 + ( j - 1 ) * ldb ] / d[ ioffd + i - 1 ]
        - B[ ioffb + i + ( j - 1 ) * ldb ] * e[ ioffe + i - 1 ];
    }
  }
}
//*************************************************************************
LaPack0.dsecnd = function() {
  var d = new Date();
  return d.getTime() / 1000.;
}
//*************************************************************************
LaPack0.dsfrk = function( transr, uplo, trans, n, k, alpha, A,
lda, beta, c, ioffa, ioffc ) {
  throw new Error("not tested: RFP format");
  info.setValue( 0 );
  var normaltransr = ( transr.charAt(0).toUpperCase() == 'N' );
  var lower = ( uplo.charAt(0).toUpperCase() == 'L' );
  var notrans = ( trans.charAt(0).toUpperCase() == 'N' );
  var nrowa = ( notrans ? n : k );
  if ( ! normaltransr && transr.charAt(0).toUpperCase() != 'T' ) {
    info.setValue( -1 );
  } else if ( ! lower && uplo.charAt(0).toUpperCase() != 'U' ) {
    info.setValue( -2 );
  } else if ( ! notrans && trans.charAt(0).toUpperCase() != 'T' ) {
    info.setValue( -3 );
  } else if ( n < 0 ) info.setValue( -4 );
  else if ( k < 0 ) info.setValue( -5 );
  else if ( lda < Math.max( 1, nrowa ) ) info.setValue( -8 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dsfrk', - info.getValue() );
    return;
  }
  if ( n == 0 || ( ( alpha == 0. || k == 0 ) && beta == 1. ) ) return;
  if ( alpha == 0. && beta == 0. ) {
    for ( var j = 1; j <= ( n * ( n + 1 ) ) / 2; j ++ ) {
      c[ ioffc + j - 1 ] = 0.;
    }
    return;
  }
  if ( n % 2 == 0 ) {
    var nisodd = false;
    var nk = n / 2;
  } else {
    nisodd = true;
    if ( lower ) {
      var n2 = n / 2;
      var n1 = n - n2;
    } else {
      n1 = n / 2;
      n2 = n - n1;
    }
  }
  if ( nisodd ) {
    if ( normaltransr ) {
      if ( lower ) {
        if ( notrans ) {
          Blas3.dsyrk( 'L', 'N', n1, k, alpha, A, lda, beta, c, n,
            ioffa, ioffc );
          Blas3.dsyrk( 'U', 'N', n2, k, alpha, A, lda, beta, c, n,
            ioffa + n1, ioffc + n );
          Blas3.dgemm( 'N', 'T', n2, n1, k, alpha, A, lda, A, lda,
            beta, c, n, ioffa + n1, ioffa, ioffc + n1 );
        } else {
          Blas3.dsyrk( 'L', 'T', n1, k, alpha, A, lda, beta, c, n,
            ioffa, ioffc );
          Blas3.dsyrk( 'U', 'T', n2, k, alpha, A, lda, beta, c, n,
            ioffa + n1 * lda, ioffc + n );
          Blas3.dgemm( 'T', 'N', n2, n1, k, alpha, A, lda, A, lda,
            beta, c, n, ioffa + n1 * lda, ioffa, ioffc + n1 );
        }
      } else {
        if ( notrans ) {
          Blas3.dsyrk( 'L', 'N', n1, k, alpha, A, lda, beta, c, n,
            ioffa, ioffc + n2 );
          Blas3.dsyrk( 'U', 'N', n2, k, alpha, A, lda, beta, c, n,
            ioffa + n2 - 1, ioffc + n1 );
          Blas3.dgemm( 'N', 'T', n1, n2, k, alpha, A, lda, A, lda,
            beta, c, n, ioffa, ioffa + n2 - 1, ioffc );
        } else {
          Blas3.dsyrk( 'L', 'T', n1, k, alpha, A, lda, beta, c, n,
            ioffa, ioffc + n2 );
          Blas3.dsyrk( 'U', 'T', n2, k, alpha, A, lda, beta, c, n,
            ioffa + ( n2 - 1 ) * lda, ioffc + n1 );
          Blas3.dgemm( 'T', 'N', n1, n2, k, alpha, A, lda, A, lda,
            beta, c, n, ioffa, ioffa + ( n2 - 1 ) * lda, ioffc );
        }
      }
    } else {
      if ( lower ) {
        if ( notrans ) {
          Blas3.dsyrk( 'U', 'N', n1, k, alpha, A, lda, beta, c, n1,
            ioffa, ioffc );
          Blas3.dsyrk( 'L', 'N', n2, k, alpha, A, lda, beta, c, n1,
            ioffa + n1, ioffc + 1 );
          Blas3.dgemm( 'N', 'T', n1, n2, k, alpha, A, lda, A, lda,
            beta, c, n1, ioffa, ioffa + n1, ioffc + n1 * n1 );
        } else {
          Blas3.dsyrk( 'U', 'T', n1, k, alpha, A, lda, beta, c, n1,
            ioffa, ioffc );
          Blas3.dsyrk( 'L', 'T', n2, k, alpha, A, lda, beta, c, n1,
            ioffa + n1 * lda, ioffc + 1 );
          Blas3.dgemm( 'T', 'N', n1, n2, k, alpha, A, lda, A, lda,
            beta, c, n1, ioffa, ioffa + n1 * lda, ioffc + n1 * n1 );
        }
      } else {
        if ( notrans ) {
          Blas3.dsyrk( 'U', 'N', n1, k, alpha, A, lda, beta, c, n2,
            ioffa, ioffc + n2 * n2 );
          Blas3.dsyrk( 'L', 'N', n2, k, alpha, A, lda, beta, c, n2,
            ioffa + n1, ioffc + n1 * n2 );
          Blas3.dgemm( 'N', 'T', n2, n1, k, alpha, A, lda, A, lda,
            beta, c, n2, ioffa + n1, ioffa, ioffc );
        } else {
          Blas3.dsyrk( 'U', 'T', n1, k, alpha, A, lda, beta, c, n2,
            ioffa, ioffc + n2 * n2 );
          Blas3.dsyrk( 'L', 'T', n2, k, alpha, A, lda, beta, c, n2,
            ioffa + n1 * lda, ioffc + n1 * n2 );
          Blas3.dgemm( 'T', 'N', n2, n1, k, alpha, A, lda, A, lda,
            beta, c, n2, ioffa + n1 * lda, ioffa, ioffc );
        }
      }
    }
  } else {
    if ( normaltransr ) {
      if ( lower ) {
        if ( notrans ) {
          Blas3.dsyrk( 'L', 'N', nk, k, alpha, A, lda, beta, c, n + 1,
            ioffa, ioffc + 1 );
          Blas3.dsyrk( 'U', 'N', nk, k, alpha, A, lda, beta, c, n + 1,
            ioffa + nk, ioffc );
          Blas3.dgemm( 'N', 'T', nk, nk, k, alpha, A, lda, A, lda,
            beta, c, n + 1, ioffa + nk, ioffa, ioffc + nk + 1 );
        } else {
          Blas3.dsyrk( 'L', 'T', nk, k, alpha, A, lda, beta, c, n + 1,
            ioffa, ioffc + 1 );
          Blas3.dsyrk( 'U', 'T', nk, k, alpha, A, lda, beta, c, n + 1,
            ioffa + nk * lda, ioffc );
          Blas3.dgemm( 'T', 'N', nk, nk, k, alpha, A, lda, A, lda,
            beta, c, n + 1, ioffa + nk * lda, ioffa, ioffc + nk + 1 );
        }
      } else {
        if ( notrans ) {
          Blas3.dsyrk( 'L', 'N', nk, k, alpha, A, lda, beta, c, n + 1,
            ioffa, ioffc + nk + 1 );
          Blas3.dsyrk( 'U', 'N', nk, k, alpha, A, lda, beta, c, n + 1,
            ioffa + nk, ioffc + nk );
          Blas3.dgemm( 'N', 'T', nk, nk, k, alpha, A, lda, A, lda,
            beta, c, n + 1, ioffa, ioffa + nk, ioffc );
        } else {
          Blas3.dsyrk( 'L', 'T', nk, k, alpha, A, lda, beta, c, n + 1,
            ioffa, ioffc + nk + 1 );
          Blas3.dsyrk( 'U', 'T', nk, k, alpha, A, lda, beta, c, n + 1,
            ioffa + nk * lda, ioffc + nk );
          Blas3.dgemm( 'T', 'N', nk, nk, k, alpha, A, lda, A, lda,
            beta, c, n + 1, ioffa, ioffa + nk * lda, ioffc );
        }
      }
    } else {
      if ( lower ) {
        if ( notrans ) {
          Blas3.dsyrk( 'U', 'N', nk, k, alpha, A, lda, beta, c, nk,
            ioffa, ioffc + nk );
          Blas3.dsyrk( 'L', 'N', nk, k, alpha, A, lda, beta, c, nk,
            ioffa + nk, ioffc );
          Blas3.dgemm( 'N', 'T', nk, nk, k, alpha, A, lda, A, lda,
            beta, c, nk, ioffa, ioffa + nk, ioffc + ( nk + 1 ) * nk );
        } else {
          Blas3.dsyrk( 'U', 'T', nk, k, alpha, A, lda, beta, c, nk,
            ioffa, ioffc + nk );
          Blas3.dsyrk( 'L', 'T', nk, k, alpha, A, lda, beta, c, nk,
            ioffa + nk * lda, ioffc );
          Blas3.dgemm( 'T', 'N', nk, nk, k, alpha, A, lda, A, lda,
            beta, c, nk, ioffa, ioffa + nk * lda,
            ioffc + ( nk + 1 ) * nk );
        }
      } else {
        if ( notrans ) {
          Blas3.dsyrk( 'U', 'N', nk, k, alpha, A, lda, beta, c, nk,
            ioffa, ioffc + nk * ( nk + 1 ) );
          Blas3.dsyrk( 'L', 'N', nk, k, alpha, A, lda, beta, c, nk,
            ioffa + nk, ioffc + nk * nk );
          Blas3.dgemm( 'N', 'T', nk, nk, k, alpha, A, lda, A, lda,
            beta, c, nk, ioffa + nk, ioffa, ioffc );
        } else {
          Blas3.dsyrk( 'U', 'T', nk, k, alpha, A, lda, beta, c, nk,
            ioffa, ioffc + nk * ( nk + 1 ) );
          Blas3.dsyrk( 'L', 'T', nk, k, alpha, A, lda, beta, c, nk,
            ioffa + nk * lda, ioffc + nk * nk );
          Blas3.dgemm( 'T', 'N', nk, nk, k, alpha, A, lda, A, lda,
            beta, c, nk, ioffa + nk * lda, ioffa, ioffc );
        }
      }
    }
  }
}
//*************************************************************************
LaPack0.dspgst = function( itype, uplo, n, AP, BP, info ) {
  throw new Error("not programmed: packed matrix");
}
//*************************************************************************
LaPack0.dsptrf = function( uplo, n, AP, ipiv, info ) {
  throw new Error("not programmed: packed matrix");
}
//*************************************************************************
LaPack0.dsptri = function( uplo, n, AP, ipiv, work, info ) {
  throw new Error("not programmed: packed matrix");
}
//*************************************************************************
LaPack0.dsptrs = function( uplo, n, nrhs, AP, ipiv, B, ldb,
info ) {
  throw new Error("not programmed: packed matrix");
}
LaPack0.zsptrs = function( uplo, n, nrhs, AP, ipiv, B, ldb,
info ) {
  throw new Error("not programmed: packed matrix");
}
//*************************************************************************
LaPack0.dsyconv = function( uplo, way, n, A, lda, ipiv, work,
info, ioffa, ioffipiv, ioffwork ) {
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var convert = ( way.charAt(0).toUpperCase() == 'C' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( ! convert && way.charAt(0).toUpperCase() != 'R' ) {
    info.setValue( -2 );
  } else if ( n < 0 ) info.setValue( -3 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -5 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dsyconv', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  if ( upper ) {
    if ( convert ) {
      var i = n;
      work[ ioffwork ] = 0.;
      while ( i > 1 ) {
        if ( ipiv[ ioffipiv + i - 1 ] < 0 ) {
          work[ ioffwork + i - 1 ] =
            A[ ioffa + i - 2 + ( i - 1 ) * lda ];
          A[ ioffa + i - 2 + ( i - 1 ) * lda ] = 0.;
          i --;
        } else work[ ioffwork + i - 1 ] = 0.;
        i--;
      }
      i = n;
      while ( i >= 1 ) {
        if ( ipiv[ ioffipiv + i - 1 ] > 0 ) {
          var ip = ipiv[ ioffipiv + i - 1 ]
          if ( i < n ) {
            for ( var j = i + 1; j <= n; j ++ ) {
              var temp = A[ ioffa + ip - 1 + ( j - 1 ) * lda ];
              A[ ioffa + ip - 1 + ( j - 1 ) * lda ] =
                A[ ioffa + i - 1 + ( j - 1 ) * lda ];
              A[ ioffa + i - 1 + ( j - 1 ) * lda ] = temp;
            }
          }
        } else {
          ip = - ipiv[ ioffipiv + i - 1 ]
          if ( i < n ) {
            for ( j = i + 1; j <= n; j ++ ) {
              temp = A[ ioffa + ip - 1 + ( j - 1 ) * lda ];
              A[ ioffa + ip - 1 + ( j - 1 ) * lda ] =
                A[ ioffa + i - 2 + ( j - 1 ) * lda ];
              A[ ioffa + i - 2 + ( j - 1 ) * lda ] = temp;
            }
          }
          i --;
        }
        i --;
      }
    } else {
      i = 1;
      while ( i <= n ) {
        if ( ipiv[ ioffipiv + i - 1 ] > 0 ) {
          ip = ipiv[ ioffipiv + i - 1 ]
          if ( i < n ) {
            for ( j = i + 1; j <= n; j ++ ) {
              temp = A[ ioffa + ip - 1 + ( j - 1 ) * lda ];
              A[ ioffa + ip - 1 + ( j - 1 ) * lda ] =
                A[ ioffa + i - 1 + ( j - 1 ) * lda ];
              A[ ioffa + i - 1 + ( j - 1 ) * lda ] = temp;
            }
          }
        } else {
          ip = - ipiv[ ioffipiv + i - 1 ]
          i ++;
          if ( i < n ) {
            for ( j = i + 1; j <= n; j ++ ) {
              temp = A[ ioffa + ip - 1 + ( j - 1 ) * lda ];
              A[ ioffa + ip - 1 + ( j - 1 ) * lda ] =
                A[ ioffa + i - 2 + ( j - 1 ) * lda ];
              A[ ioffa + i - 2 + ( j - 1 ) * lda ] = temp;
            }
          }
        }
        i ++;
      }
      i = n;
      while ( i > 1 ) {
        if ( ipiv[ ioffipiv + i - 1 ] > 1 ) {
          A[ ioffa + i - 2 + ( i - 1 ) * lda ] =
            work[ ioffwork + i - 1 ];
        }
        i --;
      }
    }
  } else {
    if ( convert ) {
      i = 1;
      work[ ioffwork + n - 1 ] = 0.;
      while ( i <= n ) {
        if ( i < n && ipiv[ ioffipiv + i - 1 ] < 0 ) {
          work[ ioffwork + i - 1 ] =
            A[ ioffa + i + ( i - 1 ) * lda ];
          A[ ioffa + i + ( i - 1 ) * lda ] = 0.;
          i ++;
        } else work[ ioffwork + i - 1 ] = 0.;
        i ++;
      }
      i = 1;
      while ( i <= n ) {
        if ( ipiv[ ioffipiv + i - 1 ] > 0 ) {
          ip = ipiv[ ioffipiv + i - 1 ];
          if ( i > 1 ) {
            for ( j = 1; j <= i - 1; j ++ ) {
              temp = A[ ioffa + ip - 1 + ( j - 1 ) * lda ];
              A[ ioffa + ip - 1 + ( j - 1 ) * lda ] =
                A[ ioffa + i - 1 + ( j - 1 ) * lda ];
              A[ ioffa + i - 1 + ( j - 1 ) * lda ] = temp;
            }
          }
        } else {
          ip = - ipiv[ ioffipiv + i - 1 ]
          if ( i > 1 ) {
            for ( j = 1; j <= i - 1; j ++ ) {
              temp = A[ ioffa + ip - 1 + ( j - 1 ) * lda ];
              A[ ioffa + ip - 1 + ( j - 1 ) * lda ] =
                A[ ioffa + i + ( j - 1 ) * lda ];
              A[ ioffa + i + ( j - 1 ) * lda ] = temp;
            }
          }
          i ++;
        }
        i ++;
      }
    } else {
      i = n;
      while ( i >= 1 ) {
        if ( ipiv[ ioffipiv + i - 1 ] > 0 ) {
          ip = ipiv[ ioffipiv + i - 1 ]
          if ( i > 1 ) {
            for ( j = 1; j <= i - 1; j ++ ) {
              temp = A[ ioffa + i - 1 + ( j - 1 ) * lda ];
              A[ ioffa + i - 1 + ( j - 1 ) * lda ] =
                A[ ioffa + ip - 1 + ( j - 1 ) * lda ];
              A[ ioffa + ip - 1 + ( j - 1 ) * lda ] = temp;
            }
          }
        } else {
          ip = - ipiv[ ioffipiv + i - 1 ]
          i --;
          if ( i > 1 ) {
            for ( j = 1; j <= i - 1; j ++ ) {
              temp = A[ ioffa + i + ( j - 1 ) * lda ];
              A[ ioffa + i + ( j - 1 ) * lda ] =
                A[ ioffa + ip - 1 + ( j - 1 ) * lda ];
              A[ ioffa + ip - 1 + ( j - 1 ) * lda ] = temp;
            }
          }
        }
        i --;
      }
      i = 1;
      while ( i <= n - 1 ) {
        if ( ipiv[ ioffipiv + i - 1 ] < 1 ) {
          A[ ioffa + i + ( i - 1 ) * lda ] =
            work[ ioffwork + i - 1 ];
          i ++;
        }
        i ++;
      }
    }
  }
}
//*************************************************************************
LaPack0.dsygs2 = function( itype, uplo, n, A, lda, B, ldb, info,
ioffa, ioffb ) {
  throw new Error("not tested: generalized eigenvalue");
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( itype < 1 || itype > 3 ) info.setValue( -1 );
  else if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -2 );
  } else if ( n < 0 ) info.setValue( -3 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -5 );
  else if ( ldb < Math.max( 1, n ) ) info.setValue( -7 );
  if ( info.getValue() != 0) {
    Blas2.xerbla( 'dsygs2', - info.getValue() );
    return;
  }
  if ( itype == 1 ) {
    if ( upper ) {
      for ( var k = 1; k <= n; k ++ ) {
        var akk = A[ ioffa + k - 1 + ( k - 1 ) * lda ];
        var bkk = B[ ioffb + k - 1 + ( k - 1 ) * ldb ];
        akk /= Math.pow( bkk, 2 );
        A[ ioffa + k - 1 + ( k - 1 ) * lda ] = akk;
        if ( k < n ) {
          Blas1.dscal( n - k, 1. / bkk, A, lda,
            ioffa + k - 1 + k * lda );
          var ct = - 0.5 * akk;
          Blas1.daxpy( n - k, ct, B, ldb, A, lda,
            ioffb + k - 1 + k * ldb, ioffa + k - 1 + k * lda );
          Blas2.dsyr2( uplo, n - k, -1., A, lda, B, ldb, A, lda,
            ioffa + k - 1 + k * lda, ioffb + k - 1 + k * ldb,
            ioffa + k + k * lda );
          Blas1.daxpy( n - k, ct, B, ldb, A, lda,
            ioffb + k - 1 + k * ldb, ioffa + k - 1 + k * lda );
          Blas2.dtrsv( uplo, 'Transpose', 'Non-unit', n - k, B, ldb,
            A, lda, ioffb + k + k * ldb, ioffa + k - 1 + k * lda );
        }
      }
    } else {
      for ( k = 1; k <= n; k ++ ) {
        akk = A[ ioffa + k - 1 + ( k - 1 ) * lda ];
        bkk = B[ ioffb + k - 1 + ( k - 1 ) * ldb ];
        akk /= Math.pow( bkk, 2 );
        A[ ioffa + k - 1 + ( k - 1 ) * lda ] = akk;
        if ( k < n ) {
          Blas1.dscal( n - k, 1. / bkk, A, 1,
            ioffa + k + ( k - 1 ) * lda );
          ct = - 0.5 * akk;
          Blas1.daxpy( n - k, ct, B, 1, A, 1,
            ioffb + k + ( k - 1 ) * ldb, ioffa + k + ( k - 1 ) * lda );
          Blas2.dsyr2( uplo, n - k, -1., A, 1, B, 1, A, lda,
            ioffa + k + ( k - 1 ) * lda, ioffb + k + ( k - 1 ) * ldb,
            ioffa + k + k * lda );
          Blas1.daxpy( n - k, ct, B, 1, A, 1,
            ioffb + k + ( k - 1 ) * ldb, ioffa + k + ( k - 1 ) * lda );
          Blas2.dtrsv( uplo, 'No transpose', 'Non-unit', n - k, B, ldb,
            A, 1, ioffb + k + k * ldb, ioffa + k + ( k - 1 ) * lda );
        }
      }
    }
  } else {
    if ( upper ) {
      for ( k = 1; k <= n; k ++ ) {
        akk = A[ ioffa + k - 1 + ( k - 1 ) * lda ];
        bkk = B[ ioffb + k - 1 + ( k - 1 ) * ldb ];
        Blas2.dtrmv( uplo, 'No transpose', 'Non-unit', k - 1, B, ldb,
          A, 1, ioffb, ioffa + ( k - 1 ) * lda );
        ct = 0.5 * akk;
        Blas1.daxpy( k - 1, ct, B, 1, A, 1, ioffb + ( k - 1 ) * ldb,
          ioffa + ( k - 1 ) * lda );
        Blas2.dsyr2( uplo, k - 1, 1., A, 1, B, 1, A, lda,
          ioffa + ( k - 1 ) * lda, ioffb + ( k - 1 ) * ldb, ioffa );
        Blas1.daxpy( k - 1, ct, B, 1, A, 1, ioffb + ( k - 1 ) * ldb,
          ioffa + ( k - 1 ) * lda );
        Blas1.dscal( k - 1, bkk, A, 1, ioffa + ( k - 1 ) * lda );
        A[ ioffa + k - 1 + ( k - 1 ) * lda ] =
          akk * Math.pow( bkk, 2 );
      }
    } else {
      for ( k = 1; k <= n; k ++ ) {
        akk = A[ ioffa + k - 1 + ( k - 1 ) * lda ];
        bkk = B[ ioffb + k - 1 + ( k - 1 ) * ldb ];
        Blas2.dtrmv( uplo, 'Transpose', 'Non-unit', k - 1, B, ldb,
          A, lda, ioffb, ioffa + k - 1 );
        ct = 0.5 * akk;
        Blas1.daxpy( k - 1, ct, B, ldb, A, lda,
          ioffb + k - 1, ioffa + k - 1 );
        Blas2.dsyr2( uplo, k - 1, 1., A, lda, B, ldb, A, lda,
          ioffa + k - 1, ioffb + k - 1, ioffa );
        Blas1.daxpy( k - 1, ct, B, ldb, A, lda,
          ioffb + k - 1, ioffa + k - 1 );
        Blas1.dscal( k - 1, bkk, A, lda, ioffa + k - 1 );
        A[ ioffa + k - 1 + ( k - 1 ) * lda ] =
          akk * Math.pow( bkk, 2 );
      }
    }
  }
}
//*************************************************************************
LaPack0.dsyswapr = function( uplo, n, A, lda, i1, i2, ioffa ) {
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( upper ) {
    Blas1.dswap( i1 - 1, A, 1, A, 1, ioffa + ( i1 - 1 ) * lda,
      ioffa + ( i2 - 1 ) * lda );
    var tmp = A[ ioffa + i1 - 1 + ( i1 - 1 ) * lda ];
    A[ ioffa + i1 - 1 + ( i1 - 1 ) * lda ] =
      A[ ioffa + i2 - 1 + ( i2 - 1 ) * lda ];
    A[ ioffa + i2 - 1 + ( i2 - 1 ) * lda ] = tmp;
    for ( var i = 1; i <= i2 - i1 - 1; i ++ ) {
      tmp = A[ ioffa + i1 - 1 + ( i1 + i - 1 ) * lda ];
      A[ ioffa + i1 - 1 + ( i1 + i - 1 ) * lda ] =
        A[ ioffa + i1 + i - 1 + ( i2 - 1 ) * lda ];
      A[ ioffa + i1 + i - 1 + ( i2 - 1 ) * lda ] = tmp;
    }
    for ( i = i2 + 1; i <= n; i ++ ) {
      tmp = A[ ioffa + i1 - 1 + ( i - 1 ) * lda ];
      A[ ioffa + i1 - 1 + ( i - 1 ) * lda ] =
        A[ ioffa + i2 - 1 + ( i - 1 ) * lda ];
      A[ ioffa + i2 - 1 + ( i - 1 ) * lda ] = tmp;
    }
  } else {
    Blas1.dswap( i1 - 1, A, lda, A, lda, ioffa + i1 - 1,
      ioffa + i2 - 1 );
    tmp = A[ ioffa + i1 - 1 + ( i1 - 1 ) * lda ];
    A[ ioffa + i1 - 1 + ( i1 - 1 ) * lda ] =
      A[ ioffa + i2 - 1 + ( i2 - 1 ) * lda ];
    A[ ioffa + i2 - 1 + ( i2 - 1 ) * lda ] = tmp;
    for ( i = 1; i <= i2 - i1 - 1; i ++ ) {
      tmp = A[ ioffa + i1 + i - 1 + ( i1 - 1 ) * lda ];
      A[ ioffa + i1 + i - 1 + ( i1 - 1 ) * lda ] =
        A[ ioffa + i2 - 1 + ( i1 + i - 1 ) * lda ];
      A[ ioffa + i2 - 1 + ( i1 + i - 1 ) * lda ] = tmp;
    }
    for ( i = i2 + 1; i <= n; i ++ ) {
      tmp = A[ ioffa + i - 1 + ( i1 - 1 ) * lda ];
      A[ ioffa + i - 1 + ( i1 - 1 ) * lda ] =
        A[ ioffa + i - 1 + ( i2 - 1 ) * lda ];
      A[ ioffa + i - 1 + ( i2 - 1 ) * lda ] = tmp;
    }
  }
}
//*************************************************************************
LaPack0.dsytf2 = function( uplo, n, A, lda, ipiv, info, ioffa,
ioffipiv ) {
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != "L" ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -4 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dsytf2', -info.getValue() );
    return;
  }
  var alpha = ( 1. + Math.sqrt( 17. ) ) / 8.;
  if ( upper ) {
    var k = n;
    while ( k >= 1 ) { // 10
      var kstep = 1;
      var absakk = Math.abs( A[ ioffa + k - 1 + ( k - 1 ) * lda ] );
      if ( k > 1 ) {
        var imax = Blas1.idamax( k - 1, A, 1, ioffa + ( k - 1 ) * lda );
        var colmax = Math.abs( A[ ioffa + imax - 1 + ( k - 1 ) * lda ] );
      } else colmax = 0.;
      if ( Math.max( absakk, colmax ) == 0. || isNaN( absakk ) ) {
        if ( info.getValue() == 0 ) info.setValue( k );
        var kp = k;
      } else {
        if ( absakk >= alpha * colmax ) kp = k;
        else {
          var jmax = imax + Blas1.idamax( k - imax, A, lda,
            ioffa + imax - 1 + imax * lda );
          var rowmax =
            Math.abs( A[ ioffa + imax - 1 + ( jmax - 1 ) * lda ] );
          if ( imax > 1 ) {
            jmax = Blas1.idamax( imax - 1, A, 1,
              ioffa + ( imax - 1 ) * lda );
            rowmax = Math.max( rowmax,
              Math.abs( A[ ioffa + jmax - 1 + ( imax - 1 ) * lda ] ) );
          }
          if ( absakk >= alpha * colmax * ( colmax / rowmax ) ) {
            kp = k;
          } else if ( Math.abs( A[ ioffa + imax - 1
          + ( imax - 1 ) * lda ] ) >= alpha * rowmax ) {
            kp = imax;
          } else {
            kp = imax;
            kstep = 2;
          }
        }
        var kk = k - kstep + 1;
        if ( kp != kk ) {
          Blas1.dswap( kp - 1, A, 1, A, 1, ioffa + ( kk - 1 ) * lda,
            ioffa + ( kp - 1 ) * lda );
          Blas1.dswap( kk - kp - 1, A, 1, A, lda,
            ioffa + kp + ( kk - 1 ) * lda, ioffa + kp - 1 + kp * lda );
          var t = A[ ioffa + kk - 1 + ( kk - 1 ) * lda ];
          A[ ioffa + kk - 1 + ( kk - 1 ) * lda ] =
            A[ ioffa + kp - 1 + ( kp - 1 ) * lda ];
          A[ ioffa + kp - 1 + ( kp - 1 ) * lda ] = t;
          if ( kstep == 2 ) {
            t = A[ ioffa + k - 2 + ( k - 1 ) * lda ];
            A[ ioffa + k - 2 + ( k - 1 ) * lda ] =
              A[ ioffa + kp - 1 + ( k - 1 ) * lda ];
            A[ ioffa + kp - 1 + ( k - 1 ) * lda ] = t;
          }
        }
        if ( kstep == 1 ) {
          var r1 = 1. / A[ ioffa + k - 1 + ( k - 1 ) * lda ];
          Blas2.dsyr( uplo, k - 1, - r1, A, 1, A, lda,
            ioffa + ( k - 1 ) * lda, ioffa );
          Blas1.dscal( k - 1, r1, A, 1, ioffa + ( k - 1 ) * lda );
        } else {
          if ( k > 2 ) {
            var d12 = A[ ioffa + k - 2 + ( k - 1 ) * lda ];
            var d22 = A[ ioffa + k - 2 + ( k - 2 ) * lda ] / d12;
            var d11 = A[ ioffa + k - 1 + ( k - 1 ) * lda ] / d12;
            t = 1. / ( d11 * d22 - 1. );
            d12 = t / d12;
            for ( var j = k - 2; j >= 1; j -- ) {
              var wkm1 =
                d12 * ( d11 * A[ ioffa + j - 1 + ( k - 2 ) * lda ]
                - A[ ioffa + j - 1 + ( k - 1 ) * lda ] );
              var wk =
                d12 * ( d22 * A[ ioffa + j - 1 + ( k - 1 ) * lda ]
                - A[ ioffa + j - 1 + ( k - 2 ) * lda ] );
              for ( var i = j; i >= 1; i -- ) {
                A[ ioffa + i - 1 + ( j - 1 ) * lda ] -=
                  A[ ioffa + i - 1 + ( k - 1 ) * lda ] * wk
                  + A[ ioffa + i - 1 + ( k - 2 ) * lda ] * wkm1;
              } // 20
              A[ ioffa + j - 1 + ( k - 1 ) * lda ] = wk;
              A[ ioffa + j - 1 + ( k - 2 ) * lda ] = wkm1;
            } // 30
          }
        }
      }
      if ( kstep == 1 ) ipiv[ ioffipiv + k - 1 ] = kp;
      else {
        ipiv[ ioffipiv + k - 1 ] = - kp;
        ipiv[ ioffipiv + k - 2 ] = - kp;
      }
      k -= kstep;
    }
  } else {
    k = 1;
    while ( k <= n ) { // 40
      kstep = 1;
      absakk = Math.abs( A[ ioffa + k - 1 + ( k - 1 ) * lda ] );
      if ( k < n ) {
        imax = k + Blas1.idamax( n - k, A, 1,
          ioffa + k + ( k - 1 ) * lda );
        colmax = Math.abs( A[ ioffa + imax - 1 + ( k - 1 ) * lda ] );
      } else colmax = 0.;
      if ( Math.max( absakk, colmax ) == 0. || isNaN( absakk ) ) {
        if ( info.getValue() == 0 ) info.setValue( k );
        kp = k;
      } else {
        if ( absakk >= alpha * colmax ) kp = k;
        else {
          jmax = k - 1 + Blas1.idamax( imax - k, A, lda,
              ioffa + imax - 1 + ( k - 1 ) * lda );
          rowmax = Math.abs( A[ ioffa + imax - 1 + ( jmax - 1 ) * lda ] );
          if ( imax < n ) {
            jmax = imax + Blas1.idamax( n - imax, A, 1,
              ioffa + imax + ( imax - 1 ) * lda );
            rowmax = Math.max( rowmax,
              Math.abs( A[ ioffa + jmax - 1 + ( imax - 1 ) * lda ] ) );
          }
          if ( absakk >= alpha * colmax * ( colmax / rowmax ) ) {
            kp = k;
          } else if ( Math.abs( A[ ioffa + imax - 1
          + ( imax - 1 ) * lda ] ) >= alpha * rowmax ) {
            kp = imax;
          } else {
            kp = imax;
            kstep = 2;
          }
        }
        kk = k + kstep - 1;
        if ( kp != kk ) {
          if ( kp < n ) {
            Blas1.dswap( n - kp, A, 1, A, 1,
              ioffa + kp + ( kk - 1 ) * lda,
              ioffa + kp + ( kp - 1 ) * lda );
          }
          Blas1.dswap( kp - kk - 1, A, 1, A, lda,
            ioffa + kk + ( kk - 1 ) * lda, ioffa + kp - 1 + kk * lda );
          t = A[ ioffa + kk - 1 + ( kk - 1 ) * lda ];
          A[ ioffa + kk - 1 + ( kk - 1 ) * lda ] =
            A[ ioffa + kp - 1 + ( kp - 1 ) * lda ];
          A[ ioffa + kp - 1 + ( kp - 1 ) * lda ] = t;
          if ( kstep == 2 ) {
            t = A[ ioffa + k + ( k - 1 ) * lda ];
            A[ ioffa + k + ( k - 1 ) * lda ] =
              A[ ioffa + kp - 1 + ( k - 1 ) * lda ];
            A[ ioffa + kp - 1 + ( k - 1 ) * lda ] = t;
          }
        }
        if ( kstep == 1 ) {
          if ( k < n ) {
            d11 = 1. / A[ ioffa + k - 1 + ( k - 1 ) * lda ];
            Blas2.dsyr( uplo, n - k, - d11, A, 1, A, lda,
              ioffa + k + ( k - 1 ) * lda, ioffa + k + k * lda );
            Blas1.dscal( n - k, d11, A, 1,
              ioffa + k + ( k - 1 ) * lda );
          }
        } else {
          if ( k < n - 1 ) {
            var d21 = A[ ioffa + k + ( k - 1 ) * lda ];
            d11 = A[ ioffa + k + k * lda ] / d21;
            d22 = A[ ioffa + k - 1 + ( k - 1 ) * lda ] / d21;
            t = 1. / ( d11 * d22 - 1. );
            d21 = t / d21;
            for ( j = k + 2; j <= n; j ++ ) {
              wk = d21 * ( d11 * A[ ioffa + j - 1 + ( k - 1 ) * lda ]
                - A[ ioffa + j - 1 + k * lda ] );
              var wkp1 =
                d21 * ( d22 * A[ ioffa + j - 1 + k * lda ]
                - A[ ioffa + j - 1 + ( k - 1 ) * lda ] );
              for ( i = j; i <= n; i ++ ) {
                A[ ioffa + i - 1 + ( j - 1 ) * lda ] -=
                  A[ ioffa + i - 1 + ( k - 1 ) * lda ] * wk
                  + A[ ioffa + i - 1 + k * lda ] * wkp1;
              } // 50
              A[ ioffa + j - 1 + ( k - 1 ) * lda ] = wk;
              A[ ioffa + j - 1 + k * lda ] = wkp1
            } // 60
          }
        }
      }
      if ( kstep == 1 ) ipiv[ ioffipiv + k - 1 ] = kp;
      else {
        ipiv[ ioffipiv + k - 1 ] = - kp;
        ipiv[ ioffipiv + k ] = - kp;
      }
      k += kstep;
    }
  }
}
LaPack0.zsytf2 = function( uplo, n, A, lda, ipiv, info ) {
  throw new Error("not programmed: complex array");
}
//*************************************************************************
LaPack0.dsytri = function( uplo, n, A, lda, ipiv, work, info,
ioffa, ioffipiv, ioffwork ) {
  throw new Error("not tested: testing uses dsytrf");
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -4 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dsytri', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  if ( upper ) {
    for ( var i = n; i >= 1; i -- ) {
      if ( ipiv[ ioffipiv + i - 1 ] > 0 &&
      A[ ioffa + i - 1 + ( i - 1 ) * lda ] == 0. ) {
        info.setValue( i );
        return;
      }
    }
  } else {
    for ( i = 1; i <= n; i ++ ) {
      if ( ipiv[ ioffipiv + i - 1 ] > 0 &&
      A[ ioffa + i - 1 + ( i - 1 ) * lda ] == 0. ) {
        info.setValue( i );
        return;
      }
    }
  }
  info.setValue( 0 );
  if ( upper ) {
    var k = 1;
    while ( k <= n ) { // 30
      if ( ipiv[ ioffipiv + k - 1 ] > 0 ) {
        A[ ioffa + k - 1 + ( k - 1 ) * lda ] =
          1. / A[ ioffa + k - 1 + ( k - 1 ) * lda ]
        if ( k > 1 ) {
          Blas1.dcopy( k - 1, A, 1, work, 1, ioffa + ( k - 1 ) * lda,
            ioffwork );
          Blas2.dsymv( uplo, k - 1, -1., A, lda, work, 1, 0., A, 1,
            ioffa, ioffwork, ioffa + ( k - 1 ) * lda );
          A[ ioffa + k - 1 + ( k - 1 ) * lda ] -=
            Blas1.ddot( k -1, work, 1, A, 1, ioffwork,
              ioffa + ( k - 1 ) * lda );
        }
        var kstep = 1;
      } else {
        var t = Math.abs( A[ ioffa + k - 1 + k * lda ] );
        var ak = A[ ioffa + k - 1 + ( k - 1 ) * lda ] / t;
        var akp1 = A[ ioffa + k + k * lda ] / t;
        var akkp1 = A[ ioffa + k - 1 + k * lda ] / t;
        var d = t * ( ak * akp1 - 1. );
        A[ ioffa + k - 1 + ( k - 1 ) * lda ] = akp1 / d;
        A[ ioffa + k + k * lda ] = ak / d;
        A[ ioffa + k - 1 + k * lda ] = - akkp1 / d;
        if ( k > 1 ) {
          Blas1.dcopy( k - 1, A, 1, work, 1, ioffa + ( k - 1 ) * lda,
            ioffwork );
          Blas2.dsymv( uplo, k - 1, -1., A, lda, work, 1, 0., A, 1,
            ioffa, ioffwork, ioffa + ( k - 1 ) * lda );
          A[ ioffa + k - 1 + ( k - 1 ) * lda ] -=
            Blas1.ddot( k - 1, work, 1, A, 1, ioffwork,
              ioffa + ( k - 1 ) * lda );
          A[ ioffa + k - 1 + k * lda ] -=
            Blas1.ddot( k - 1, A, 1, A, 1, ioffa + ( k - 1 ) * lda,
              ioffa + k * lda );
          Blas1.dcopy( k - 1, A, 1, work, 1, ioffa + k * lda,
            ioffwork );
          Blas2.dsymv( uplo, k - 1, -1., A, lda, work, 1, 0., A, 1,
            ioffa, ioffwork, ioffa + k * lda );
          A[ ioffa + k + k * lda ] -=
            Blas1.ddot( k - 1, work, 1, A, 1, ioffwork,
            ioffa + k * lda );
        }
        kstep = 2;
      }
      var kp = Math.abs( ipiv[ ioffipiv + k - 1 ] );
      if ( kp != k ) {
        Blas1.dswap( kp - 1, A, 1, A, 1, ioffa + ( k - 1 ) * lda,
          ioffa + ( kp - 1 ) * lda );
        Blas1.dswap( k - kp - 1, A, 1, A, lda,
          ioffa + kp + ( k - 1 ) * lda,
          ioffa + kp - 1 + kp * lda );
        var temp = A[ ioffa + k - 1 + ( k - 1 ) * lda ];
        A[ ioffa + k - 1 + ( k - 1 ) * lda ] =
          A[ ioffa + kp - 1 + ( kp - 1 ) * lda ];
        A[ ioffa + kp - 1 + ( kp - 1 ) * lda ] = temp;
        if ( kstep == 2 ) {
          temp = A[ ioffa + k - 1 + k * lda ];
          A[ ioffa + k - 1 + k * lda ] = A[ ioffa + kp - 1 + k * lda ];
          A[ ioffa + kp - 1 + k * lda ] = temp;
        }
      }
      k += kstep;
    }
  } else {
    k = n;
    while ( k >= 1 ) { // 50
      if ( ipiv[ ioffipiv + k - 1 ] > 0 ) {
        A[ ioffa + k - 1 + ( k - 1 ) * lda ] =
          1. / A[ ioffa + k - 1 + ( k - 1 ) * lda ]
        if ( k < n ) {
          Blas1.dcopy( n - k, A, 1, work, 1,
            ioffa + k + ( k - 1 ) * lda, ioffwork );
          Blas2.dsymv( uplo, n - k, -1., A, lda, work, 1, 0., A, 1,
            ioffa + k + k * lda, ioffwork,
            ioffa + k + ( k - 1 ) * lda );
          A[ ioffa + k - 1 + ( k - 1 ) * lda ] -=
            Blas1.ddot( n - k, work, 1, A, 1, ioffwork,
              ioffa + k + ( k - 1 ) * lda );
        }
        kstep = 1;
      } else {
        t = Math.abs( A[ ioffa + k - 1 + ( k - 2 ) * lda ] );
        ak = A[ ioffa + k - 2 + ( k - 2 ) * lda ] / t;
        akp1 = A[ ioffa + k - 1 + ( k - 1 ) * lda ] / t;
        akkp1 = A[ ioffa + k - 1 + ( k - 2 ) * lda ] / t;
        d = t * ( ak * akp1 - 1. );
        A[ ioffa + k - 2 + ( k - 2 ) * lda ] = akp1 / d;
        A[ ioffa + k - 1 + ( k - 1 ) * lda ] = ak / d;
        A[ ioffa + k - 1 + ( k - 2 ) * lda ] = - akkp1 / d;
        if ( k  < n ) {
          Blas1.dcopy( n - k, A, 1, work, 1,
            ioffa + k + ( k - 1 ) * lda, ioffwork );
          Blas2.dsymv( uplo, n - k, -1., A, lda, work, 1, 0., A, 1,
            ioffa + k + k * lda, ioffwork,
            ioffa + k + ( k - 1 ) * lda );
          A[ ioffa + k - 1 + ( k - 1 ) * lda ] -=
            Blas1.ddot( n - k, work, 1, A, 1, ioffwork,
              ioffa + k + ( k - 1 ) * lda );
          A[ ioffa + k - 1 + ( k - 2 ) * lda ] -=
            Blas1.ddot( n - k, A, 1, A, 1, ioffa + k + ( k - 1 ) * lda,
              ioffa + k + ( k - 2 ) * lda );
          Blas1.dcopy( n - k, A, 1, work, 1,
            ioffa + k + ( k - 2 ) * lda, ioffwork );
          Blas2.dsymv( uplo, n - k, -1., A, lda, work, 1, 0., A, 1,
            ioffa + k + k * lda, ioffwork,
            ioffa + k + ( k - 2 ) * lda );
          A[ ioffa + k - 2 + ( k - 2 ) * lda ] -=
            Blas1.ddot( n - k, work, 1, A, 1, ioffwork,
            ioffa + k + ( k - 2 ) * lda );
        }
        kstep = 2;
      }
      kp = Math.abs( ipiv[ ioffipiv + k - 1 ] );
      if ( kp != k ) {
        if ( kp < n ) {
          Blas1.dswap( n - kp, A, 1, A, 1,
            ioffa + kp + ( k - 1 ) * lda,
            ioffa + kp + ( kp - 1 ) * lda );
        }
        Blas1.dswap( kp - k - 1, A, 1, A, lda,
          ioffa + k + ( k - 1 ) * lda,
          ioffa + kp - 1 + k * lda );
        temp = A[ ioffa + k - 1 + ( k - 1 ) * lda ];
        A[ ioffa + k - 1 + ( k - 1 ) * lda ] =
          A[ ioffa + kp - 1 + ( kp - 1 ) * lda ];
        A[ ioffa + kp - 1 + ( kp - 1 ) * lda ] = temp;
        if ( kstep == 2 ) {
          temp = A[ ioffa + k - 1 + ( k - 2 ) * lda ];
          A[ ioffa + k - 1 + ( k - 2 ) * lda ] =
            A[ ioffa + kp - 1 + ( k - 2 ) * lda ];
          A[ ioffa + kp - 1 + ( k - 2 ) * lda ] = temp;
        }
      }
      k -= kstep;
    }
  }
}
//*************************************************************************
LaPack0.dsytrs = function( uplo, n, nrhs, A, lda, ipiv, B, ldb,
info, ioffa, ioffipiv, ioffb ) {
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( nrhs < 0 ) info.setValue( -3 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -5 );
  else if ( ldb < Math.max( 1, n ) ) info.setValue( -8 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dsytrs', -info.getValue() );
    return;
  }
  if ( n == 0 || nrhs == 0 ) return;
  if ( upper ) {
    var k = n;
    while ( true ) { // 10
      if ( k < 1 ) break;
      if ( ipiv[ ioffipiv + k - 1 ] > 0 ) {
        var kp = ipiv[ ioffipiv + k - 1 ];
        if ( kp != k ) {
          Blas1.dswap( nrhs, B, ldb, B, ldb, ioffb + k - 1,
            ioffb + kp - 1 );
        }
        Blas2.dger( k - 1, nrhs, -1., A, 1, B, ldb, B, ldb,
          ioffa + ( k - 1 ) * lda, ioffb + k - 1, ioffb );
        Blas1.dscal( nrhs, 1. / A[ ioffa + k - 1 + ( k - 1 ) * lda ] ,
          B, ldb, ioffb + k - 1 );
        k --;
      } else {
        kp = - ipiv[ ioffipiv + k - 1 ];
        if ( kp != k - 1 ) {
          Blas1.dswap( nrhs, B, ldb, B, ldb, ioffb + k - 2,
            ioffb + kp - 1 );
        }
        Blas2.dger( k - 2, nrhs, -1., A, 1, B, ldb, B, ldb,
          ioffa + ( k - 1 ) * lda, ioffb + k - 1, ioffb );
        Blas2.dger( k - 2, nrhs, -1., A, 1, B, ldb, B, ldb,
          ioffa + ( k - 2 ) * lda, ioffb + k - 2, ioffb );
        var akm1k = A[ ioffa + k - 2 + ( k - 1 ) * lda ];
        var akm1 = A[ ioffa + k - 2 + ( k - 2 ) * lda ] / akm1k;
        var ak = A[ ioffa + k - 1 + ( k - 1 ) * lda ] / akm1k;
        var denom = akm1 * ak - 1.;
        for ( var j = 1; j <= nrhs; j ++ ) {
          var bkm1 =
            B[ ioffb + k - 2 + ( j - 1 ) * ldb ] / akm1k;
          var bk = B[ ioffb + k - 1 + ( j - 1 ) * ldb ] / akm1k;
          B[ ioffb + k - 2 + ( j - 1 ) * ldb ] =
            ( ak * bkm1 - bk ) / denom;
          B[ ioffb + k - 1 + ( j - 1 ) * ldb ] =
            ( akm1 * bk - bkm1 ) / denom;
        } // 20
        k -= 2;
      }
    } // 30
    k = 1;
    while ( true ) { // 40
      if ( k > n ) break;
      if ( ipiv[ ioffipiv + k - 1 ] > 0 ) {
        Blas2.dgemv( 'Transpose', k - 1, nrhs, -1., B, ldb, A, 1, 1.,
          B, ldb, ioffb, ioffa + ( k - 1 ) * lda, ioffb + k - 1 );
        kp = ipiv[ ioffipiv + k - 1 ];
        if ( kp != k ) {
          Blas1.dswap( nrhs, B, ldb, B, ldb, ioffb + k - 1,
            ioffb + kp - 1 );
        }
        k ++;
      } else {
        Blas2.dgemv( 'Transpose', k - 1, nrhs, -1., B, ldb, A, 1, 1.,
          B, ldb, ioffb, ioffa + ( k - 1 ) * lda, ioffb + k - 1 );
        Blas2.dgemv( 'Transpose', k - 1, nrhs, -1., B, ldb, A, 1, 1.,
          B, ldb, ioffb, ioffa + k * lda, ioffb + k );
        kp = - ipiv[ ioffipiv + k - 1 ];
        if ( kp != k ) {
          Blas1.dswap( nrhs, B, ldb, B, ldb, ioffb + k - 1,
            ioffb + kp - 1 );
        }
        k += 2;
      }
    }
  } else {
    k = 1;
    while ( true ) { // 60
      if ( k > n ) break;
      if ( ipiv[ ioffipiv + k - 1 ] > 0 ) {
        kp = ipiv[ ioffipiv + k - 1 ];
        if ( kp != k ) {
          Blas1.dswap( nrhs, B, ldb, B, ldb, ioffb + k - 1,
            ioffb + kp - 1 );
        }
        if ( k < n ) {
          Blas2.dger( n - k, nrhs, -1., A, 1, B, ldb, B, ldb,
            ioffa + k + ( k - 1 ) * lda, ioffb + k - 1, ioffb + k );
        }
        Blas1.dscal( nrhs, 1. / A[ ioffa + k - 1 + ( k - 1 ) * lda ],
          B, ldb, ioffb + k - 1 );
        k ++;
      } else {
        kp = - ipiv[ ioffipiv + k - 1 ];
        if ( kp != k + 1 ) {
          Blas1.dswap( nrhs, B, ldb, B, ldb, ioffb + k,
            ioffb + kp - 1 );
        }
        if ( k < n - 1 ) {
          Blas2.dger( n - k - 1, nrhs, -1., A, 1, B, ldb, B, ldb,
            ioffa + k + 1 + ( k - 1 ) * lda, ioffb + k - 1,
            ioffb + k + 1 );
          Blas2.dger( n - k - 1, nrhs, -1., A, 1, B, ldb, B, ldb,
            ioffa + k + 1 + k * lda, ioffb + k, ioffb + k + 1 );
        }
        akm1k = A[ ioffa + k + ( k - 1 ) * lda ];
        akm1 = A[ ioffa + k - 1 + ( k - 1 ) * lda ] / akm1k;
        ak = A[ ioffa + k + k * lda ] / akm1k;
        denom = akm1 * ak - 1.;
        for ( j = 1; j <= nrhs; j ++ ) {
          bkm1 = B[ ioffb + k - 1 + ( j - 1 ) * ldb ] / akm1k;
          bk = B[ ioffb + k + ( j - 1 ) * ldb ] / akm1k;
          B[ ioffb + k - 1 + ( j - 1 ) * ldb ] =
            ( ak * bkm1 - bk ) / denom;
          B[ ioffb + k + ( j - 1 ) * ldb ] =
            ( akm1 * bk - bkm1 ) / denom;
        }
        k += 2;
      }
    }
    k = n;
    while ( true ) { // 90
      if ( k < 1 ) break;
      if ( ipiv[ ioffipiv + k - 1 ] > 0 ) {
        if ( k < n ) {
          Blas2.dgemv( 'Transpose', n - k, nrhs, -1., B, ldb, A, 1, 1.,
            B, ldb, ioffb + k, ioffa + k + ( k - 1 ) * lda,
            ioffb + k - 1 );
        }
        kp = ipiv[ ioffipiv + k - 1 ];
        if ( kp != k ) {
          Blas1.dswap( nrhs, B, ldb, B, ldb, ioffb + k - 1,
            ioffb + kp - 1 );
        }
        k --;
      } else {
        if ( k < n ) {
          Blas2.dgemv( 'Transpose', n - k, nrhs, -1., B, ldb, A, 1, 1.,
            B, ldb, ioffb + k, ioffa + k + ( k - 1 ) * lda,
            ioffb + k - 1 );
          Blas2.dgemv( 'Transpose', n - k, nrhs, -1., B, ldb, A, 1, 1.,
            B, ldb, ioffb + k, ioffa + k + ( k - 2 ) * lda,
            ioffb + k - 2 );
        }
        kp = - ipiv[ ioffipiv + k - 1 ];
        if ( kp != k ) {
          Blas1.dswap( nrhs, B, ldb, B, ldb, ioffb + k - 1,
            ioffb + kp - 1 );
        }
        k -= 2;
      }
    }
  }
}
LaPack0.zsytrs = function( uplo, n, nrhs, A, lda, ipiv, B, ldb,
info ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dtbtrs = function( uplo, trans, diag, n, kd, nrhs, AB,
ldab, B, ldb, info, ioffab, ioffb ) {
  throw new Error("not programmed: band matrix");
}
//*************************************************************************
LaPack0.dtfsm = function( transr, side, uplo, trans, diag, m, n,
alpha, A, B, ldb, ioffa, ioffb ) {
  throw new Error("not programmed: RFP format");
}
//*************************************************************************
LaPack0.dtfttp = function( transr, uplo, n, ARF, AP, info,
ioffarf, ioffap) {
  throw new Error("not programmed: RFP and packed format");
}
//*************************************************************************
LaPack0.dtfttr = function( transr, uplo, n, ARF, A, lda, info,
ioffarf, ioffa ) {
  throw new Error("not programmed: RFP format");
}
//*************************************************************************
LaPack0.dtprfb = function( side, trans, direct, storev, m, n, k,
l, V, ldv, T, ldt, A, lda, B, ldb, work, ldwork, ioffv, iofft, ioffa,
ioffb, ioffwork ) {
  throw new Error("not programmed: hope we wont need this long code");
}
//*************************************************************************
LaPack0.dtptri = function( uplo, diag, n, AP, info ) {
  throw new Error("not programmed: packed matrix");
}
LaPack0.ztptri = function( uplo, diag, n, AP, info ) {
  throw new Error("not programmed: complex packed matrix");
}
//*************************************************************************
LaPack0.dtptrs = function( uplo, trans, diag, n, nrhs, AP, B,
ldb, info, ioffap, ioffb ) {
  throw new Error("not programmed: packed matrix");
}
//*************************************************************************
LaPack0.dtpttf = function( transr, uplo, n, AP, ARF, info,
ioffap, ioffarf ) {
  throw new Error("not programmed: packed format");
}
//*************************************************************************
LaPack0.dtpttr = function( uplo, n, AP, A, lda, info, ioffap,
ioffa ) {
  throw new Error("not programmed: packed format");
}
//*************************************************************************
LaPack0.dtrti2 = function( uplo, diag, n, A, lda, info, ioffa ) {
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var nounit = ( diag.charAt(0).toUpperCase() == 'N' );
  if ( !upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( ! nounit && diag.charAt(0).toUpperCase() != 'U' ) {
    info.setValue( -2 );
  } else if ( n < 0 ) info.setValue( -3 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -5 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dtrti2', -info.getValue() );
    return;
  }
  if ( upper ) {
    for ( var j = 1; j <= n; j ++ ) {
      if ( nounit ) {
        A[ ioffa + j - 1 + ( j - 1 ) * lda ] =
          1. / A[ ioffa + j - 1 + ( j - 1 ) * lda ];
        var ajj = - A[ ioffa + j - 1 + ( j - 1 ) * lda ];
      } else ajj = -1.;
      Blas2.dtrmv( 'Upper', 'No transpose', diag, j - 1, A, lda, A, 1,
        ioffa, ioffa + ( j - 1 ) * lda );
      Blas1.dscal( j - 1, ajj, A, 1, ioffa + ( j - 1 ) * lda );
    }
  } else {
    for ( j = n; j >= 1; j -- ) {
      if ( nounit ) {
        A[ ioffa + j - 1 + ( j - 1 ) * lda ] =
          1. / A[ ioffa + j - 1 + ( j - 1 ) * lda ];
        ajj = - A[ ioffa + j - 1 + ( j - 1 ) * lda ];
      } else ajj = -1.;
      if ( j < n ) {
        Blas2.dtrmv( 'Lower', 'No transpose', diag, n - j, A, lda,
          A, 1, ioffa + j + j * lda, ioffa + j + ( j - 1 ) * lda );
        Blas1.dscal( n - j, ajj, A, 1, ioffa + j + ( j - 1 ) * lda );
      }
    }
  }
}
LaPack0.ztrti2 = function( uplo, diag, n, A, lda, info ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.dtrtrs = function( uplo, trans, diag, n, nrhs, A, lda,
B, ldb, info, ioffa, ioffb ) {
  info.setValue( 0 );
  var nounit = ( diag.charAt(0).toUpperCase() == 'N' );
  if ( uplo.charAt(0).toUpperCase() != 'U' &&
  uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( trans.charAt(0).toUpperCase() != 'N' &&
  trans.charAt(0).toUpperCase() != 'T' &&
  trans.charAt(0).toUpperCase() != 'C' ) {
    info.setValue( -2 );
  } else if ( ! nounit && diag.charAt(0).toUpperCase() != 'U' ) {
    info.setValue( -3 );
  } else if ( n < 0 ) info.setValue( -4 );
  else if ( nrhs < 0 ) info.setValue( -5 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -7 );
  else if ( ldb < Math.max( 1, n ) ) info.setValue( -9 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dtrtrs', - info.getValue() );
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
  Blas3.dtrsm( 'Left', uplo, trans, diag, n, nrhs, 1., A, lda, B, ldb,
    ioffa, ioffb );
}
//*************************************************************************
LaPack0.dtrttf = function( transr, uplo, n, A, lda, ARF, info,
ioffa, ioffarf) {
  throw new Error("not programmed: packed matrix");
}
//*************************************************************************
LaPack0.dtrttp = function( uplo, n, A, lda, AP, info, ioffa,
ioffap ) {
  throw new Error("not programmed: packed matrix");
}
//*************************************************************************
LaPack0.dzsum1 = function( n, cx, incx ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.icmax1 = function( n, x, incx, ioffx) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.ieeeck = function( ispec, zero, one) {
  var posinf = one / zero;
  if ( posinf <= one ) return 0;
  var neginf = - one / zero;
  if ( neginf >= zero ) return 0;
  var negzro = one / ( neginf + one );
  if ( negzro != zero ) return 0;
  neginf = one / negzro;
  if ( neginf >= zero ) return 0;
  var newzro = negzro + zero;
  if ( newzro != zero ) return 0;
  posinf = one / newzro;
  if ( posinf <= one ) return 0;
  neginf = neginf * posinf;
  if ( neginf >= zero ) return 0;
  posinf = posinf * posinf
  if ( posinf <= one ) return 0;

  if ( ispec == 0 ) return 1;
  var nan1 = posinf + neginf;
  var nan2 = posinf / neginf;
  var nan3 = posinf / posinf;
  var nan4 = posinf * zero;
  var nan5 = neginf * negzro;
  var nan6 = nan5 * zero;
  if ( nan1 == nan1 ) return 0;
  if ( nan2 == nan2 ) return 0;
  if ( nan3 == nan3 ) return 0;
  if ( nan4 == nan4 ) return 0;
  if ( nan5 == nan5 ) return 0;
  if ( nan6 == nan6 ) return 0;
  return 1;
}
//*************************************************************************
LaPack0.ilaclc = function( m, n, A, lda, ioffa ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.ilaclr = function( m, n, A, lda, ioffa ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.iladiag = function( diag ) {
  if ( diag.charAt(0).toUpperCase() == 'N' ) return 131;
  if ( diag.charAt(0).toUpperCase() == 'U' ) return 132;
  return -1;
}
//*************************************************************************
LaPack0.iladlc = function( m, n, A, lda, ioffa ) {
  if ( n == 0 ) return n;
  else if ( A[ ioffa + ( n - 1 ) * lda ] != 0. ||
  A[ ioffa + m - 1 + ( n - 1 ) * lda ] != 0. ) {
    return n;
  } else {
    for ( var j = n; j >= 1; j -- ) {
      for ( var i = 1; i <= m; i ++ ) {
        if ( A[ ioffa + i - 1 + ( j - 1 ) * lda ] != 0. ) return j;
      }
    }
    return 0;
  }
}
//*************************************************************************
LaPack0.iladlr = function( m, n, A, lda, ioffa ) {
  if ( m == 0 ) return m;
  else if ( A[ ioffa + m - 1 ] != 0. ||
  A[ ioffa + m - 1 + ( n - 1 ) * lda ] != 0. ) {
    return m;
  } else {
    var value = 0;
    for ( var j = 1; j <= 1; j ++ ) {
      var i = m;
      while ( A[ ioffa + i - 1 + ( j - 1 ) * lda ] != 0 && i >= 1 ) {
        i --;
        if ( i == 0 ) break;
      }
      value = Math.max( value, i );
    }
    return value;
  }
}
//*************************************************************************
LaPack0.ilaenv = function( ispec, name, opts, n1, n2, n3, n4) {
  switch ( ispec ) {
    case 1:
    case 2:
    case 3: {
      var nb = 1;
      var subnam = name.toUpperCase();
      var c1 = subnam.charAt(0);
      var sname = ( c1 == "S" || c1 == "D" );
      var cname = ( c1 == "C" || c1 == "Z" );
      if ( ! ( cname || sname ) ) return nb;
      var c2 = subnam.slice(1,3);
      var c3 = subnam.slice(3,6);
      var c4 = c3.slice(1,3);
      switch (ispec) {
        case 1: { // block size
          nb = 1;
          if ( c2 == "GE" ) {
            if ( c3 == "TRF" ) nb = ( sname ? 64 : 64 );
            else if ( c3 == "QRF" || c3 == "RQF" || c3 == "LQF" ||
            c3 == "QLF" ) {
              nb = ( sname ? 32 : 32 );
            } else if ( c3 == "HRD" ) nb = ( sname ? 32 : 32 );
            else if ( c3 == "BRD" ) nb = ( sname ? 32 : 32 );
            else if ( c3 == "TRI" ) nb = ( sname ? 64 : 64 );
          } else if ( c2 == "PO" ) {
            if ( c3 == "TRF" ) nb = ( sname ? 64 : 64 );
          } else if ( c2 == "SY" ) {
            if ( c3 == "TRF" ) nb = ( sname ? 64 : 64 );
            else if ( sname && c3 == "TRD" ) nb = 32;
            else if ( sname && c3 == "GST" ) nb = 64;
          } else if ( cname && c2 == "HE" ) {
            if ( c3 == "TRF" ) nb = 64;
            else if ( c3 == "TRD" ) nb = 32;
            else if ( c3 == "GST" ) nb = 64;
          } else if ( sname && c2 == "OR" ) {
            if ( c3.charAt(0) == "G" ) {
              if ( c4 == "QR" || c4 == "RQ" || c4 == "LQ" ||
              c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR" ) {
                nb = 32;
              }
            } else if (c3.charAt(0) == "M" ) {
              if ( c4 == "QR" || c4 == "RQ" || c4 == "LQ" ||
              c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR" ) {
                nb = 32;
              }
            }
          } else if ( cname && c2 == "UN" ) {
            if ( c3.charAt(0) == "G" ) {
              if ( c4 == "QR" || c4 == "RQ" || c4 == "LQ" ||
              c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR" ) {
                nb = 32;
              }
            } else if (c3.charAt(0) == "M" ) {
              if ( c4 == "QR" || c4 == "RQ" || c4 == "LQ" ||
              c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR" ) {
                nb = 32;
              }
            }
          } else if ( c2 == "GB" ) {
            if ( c3 == "TRF" ) {
              if (sname) nb = ( n4 <= 64 ? 1 : 32 );
              else nb = ( n4 <= 64 ? 1 : 32 );
            }
          } else if ( c2 == "PB" ) {
            if ( c3 == "TRF" ) {
              if (sname) nb = ( n2 <= 64 ? 1 : 32 );
              else nb = ( n2 <= 64 ? 1 : 32 );
            }
          } else if ( c2 == "TR" ) {
            if ( c3 == "TRI" ) nb = ( sname ? 64 : 64 );
          } else if ( c2 == "LA" ) {
            if ( c3 == "UUM" ) nb = ( sname ? 64 : 64 );
          } else if ( sname && c2 == "ST" ) {
            if ( c3 == "EBZ" ) nb = 1;
          }
          return nb;
        }
        case 2: { // minimum block size
          var nbmin = 2;
          if ( c2 == "GE" ) {
            if ( c3 == "QRF" || c3 == "RQF" || c3 == "LQF" ||
            c3 == "QLF" ) {
              nbmin = ( sname ? 2 : 2);
            } else if ( c3 == "HRD" ) nbmin = ( sname ? 2 : 2);
            else if ( c3 == "BRD" ) nbmin = ( sname ? 2 : 2);
            else if ( c3 == "TRI" ) nbmin = ( sname ? 2 : 2);
          } else if ( c2 == "SY" ) {
            if ( c3 == "TRF" ) nbmin = ( sname ? 8 : 8 );
            else if ( sname && c3 == "TRD" ) nbmin = 2;
          } else if ( cname && c2 == "HE" ) {
            if ( c3 == "TRD" ) nbmin = 2;
          } else if ( sname && c2 == "OR" ) {
            if ( c3.charAt(0) == "G" ) {
              if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" ||
              c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR" ) {
                nbmin = 2;
              }
            } else if ( c3.charAt(0) == "M" ) {
              if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" ||
              c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR" ) {
                nbmin = 2;
              }
            }
          } else if ( cname && c2 == "UN" ) {
            if ( c3.charAt(0) == "G" ) {
              if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" ||
              c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR" ) {
                nbmin = 2;
              }
            } else if ( c3.charAt(0) == "M" ) {
              if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" ||
              c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR" ) {
                nbmin = 2;
              }
            }
          }
          return nbmin;
        }
        case 3: { // crossover point
          var nx = 0;
          if ( c2 == "GE" ) {
            if ( c3 == "QRF" || c3 == "RQF" || c3 == "LQF" ||
            c3 == "QLF" ) {
              nx = ( sname ? 128 : 128 );
            } else if ( c3 == "HRD" ) nx = ( sname ? 128 : 128 );
            else if ( c3 == "BRD" ) nx = ( sname ? 128 : 128 );
          } else if ( c2 == "SY" ) {
            if ( sname && c3 == "TRD" ) nx = 32;
          } else if ( cname && c2 == "HE" ) {
            if ( c3 == "TRD" ) nx = 32;
          } else if ( sname && c2 == "OR" ) {
            if ( c3.charAt(0) == "G" ) {
              if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" ||
              c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR" ) {
                nx = 128;
              }
            }
          } else if ( cname && c2 == "UN" ) {
            if ( c3.charAt(0) == "G" ) {
              if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" ||
              c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR" ) {
                nx = 128;
              }
            }
          }
          return nx;
        }
      }
    }
    case 4: // number of shifts ( used by xHSEQR )
      return 6;
    case 5: // minimum column dimension ( not used )
      return 2;
    case 6: // crossover point for SVD ( used by xGELSS and xGESVD )
      return Math.round ( Number( Math.min( n1, n2 ) ) * 1.6 );
    case 7: // number of processors ( not used )
      return 1;
    case 8: // crossover point for multishift ( used by xHSEQR )
      return 50;
    case 9: // maximum size of subproblems at bottom of comp tree
      return 25;
    case 10: // IEEE NaN arithmetic can be trusted not to trap
      return 1;
    case 11: // infinity arithmetic can be trusted not to trap
      return 1;
    case 12:
    case 13:
    case 14:
    case 15:
    case 16:
      return iparmq( ispec, name, opts, n1, n2, n3, n4 );
    default:
      return -1;
  }
}
//*************************************************************************
LaPack0.ilaprec = function( prec ) {
  if ( prec.charAt(0).toUpperCase() == 'S' ) return 211;
  if ( prec.charAt(0).toUpperCase() == 'D' ) return 212;
  if ( prec.charAt(0).toUpperCase() == 'I' ) return 213;
  if ( prec.charAt(0).toUpperCase() == 'X' ||
  prec.charAt(0).toUpperCase() == 'E' ) {
    return 214;
  }
  return -1;
}
//*************************************************************************
LaPack0.ilaslc = function( m, n, A, lda, ioffa ) {
  throw new Error("not programmed: single precision matrix");
}
//*************************************************************************
LaPack0.ilaslr = function( m, n, A, lda, ioffa ) {
  throw new Error("not programmed: single precision matrix");
}
//*************************************************************************
LaPack0.ilatrans = function( trans ) {
  if ( trans.charAt(0).toUpperCase() == 'N' ) return 111;
  if ( trans.charAt(0).toUpperCase() == 'T' ) return 112;
  if ( trans.charAt(0).toUpperCase() == 'C' ) return 113;
  return -1;
}
//*************************************************************************
LaPack0.ilauplo = function( uplo ) {
  if ( uplo.charAt(0).toUpperCase() == 'U' ) return 121;
  if ( uplo.charAt(0).toUpperCase() == 'L' ) return 122;
  return -1;
}
/*
//*************************************************************************
LaPack0.ilaver = function( stuff missing ) {
  barf;
}
*/
//*************************************************************************
LaPack0.ilazlc = function( m, n, A, lda, ioffa ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.ilazlr = function( m, n, A, lda, ioffa ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.iparmq = function( ispec, name, opts, n, ilo, ihi,
lwork ) {
  var inmin = 12;
  var inwin = 13;
  var inibl = 14;
  var ishfts = 15;
  var iacc22 = 16;
  var nmin = 75;
  var k22min = 14;
  var kacmin = 14;
  var nibble = 14;
  var knwswp = 500;
  if ( ispec == ishfts || ispec == inwin || ispec == iacc22 ) {
    var nh = ihi - ilo + 1;
    var ns = 2;
    if ( nh >= 30 ) ns = 4;
    if ( nh >= 60 ) ns = 10;
    if ( nh >= 150) {
      ns = Math.max( 10,
        nh / Math.round( Math.log( Number( nh ) ) / Math.log( 2. ) ) );
    }
    if ( nh >= 590 ) ns = 64;
    if ( nh >= 3000 ) ns = 128;
    if ( nh >= 6000 ) ns = 256;
    ns = Math.max( 2, ns - ns % 2 );
  }
  if ( ispec == inmin ) return nmin;
  else if ( ispec == inibl ) return nibble;
  else if ( ispec == ishfts ) return ns;
  else if ( ispec == inwin ) {
    return ( nh <= knwswp ? ns : 3 * ns / 2 );
  } else if ( ispec == iacc22 ) {
    if ( ns >= kacmin ) return 1;
    else if ( ns >= k22min ) return 2;
    return 0;
  } else return -1;
}
//*************************************************************************
LaPack0.izmax1 = function( n, cx, incx ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack0.lsamen = function( n, a, b ) {
  return a.substr(0,n).toUpperCase() == b.substr(0,n).toUpperCase();
}
//*************************************************************************
LaPack0.xerbla_array = function( srname_array, srname_len, info,
ioffsrname_array ) {
  var srname = '';
  for ( var i = 1; i <= Math.min( srname_len, 32 ); i ++ ) {
    srname =
      srname + srname_array[ ioffsrname_array + i - 1 ].charAt(0);
  }
  Blas2.xerbla( srname, info );
}
function testLaPack0() {
  test_dgebak();
//test_dggbak();
//test_dgttrf();
//test_dgtts2();
//test_dla_gerpvgrw();
//test_dla_wwaddw();
//test_dlacon();
//test_dlacn2();
//test_dlacpy();
//test_dladiv();
//test_dlae2();
//test_dlaebz();
//test_dlaed5();
/*test_dlaeda(); */
//test_dlaev2();
//test_dlag2();
//test_dlagtm();
//test_dlamch();
//test_dlamrg();
//test_dlapmr();
//test_dlapmt();
//test_dlapy2();
//test_dlapy3();
//test_dlaqr1();
//test_dlar2v();
//test_dlarf();
//test_dlarfb();
//test_dlarft();
//test_dlargv();
//test_dlarra();
//test_dlarrc();
//test_dlarrj();
//test_dlarscl2();
//test_dlartv();
//test_dlarz();
//test_dlarzb();
//test_dlarzt();
//test_dlas2();
//test_dlascl2();
//test_dlasd5();
//test_dlasdt();
//test_dlaset();
//test_dlasr();
//test_dlasrt();
//test_dlassq();
//test_dlaswp();
//test_dlauu2();
//test_dpoequ();
//test_dpotf2();
//test_dpotrs();
//test_dptcon();
//test_dpttrf();
//test_dptts2();
//test_dsyconv();
//test_dsyswapr();
//test_dsytf2();
//test_dsytri();
//test_dsytrs();
//test_dtrti2();
//test_dtrtrs();
//test_ieeeck();
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgebak() {
  document.getElementById("debug_textarea").value +=
    "testing dgebak *************" + "\n";
  var n = 4;
  var ilo = 2;
  var ihi = 3;
  var ioffscale = 1;
  var scale = new Array( ioffscale + n );
  var m = 2;
  var ioffV = 2;
  var ldv = 5;
  var V = new Array( ioffV + ldv * m );
  var info = new IntReference();
  for ( var j = 0; j < m; j ++ ) {
    for ( var i = 0; i < ldv; i ++ ) {
      V[ ioffV + i + j * ldv ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( i = 0; i < n; i ++ ) scale[ ioffscale + i ] = n - i;

  for ( j = 0; j < m; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  LaPack0.dgebak( 'N', 'R', n, ilo, ihi, scale, m, V, ldv, info,
    ioffscale,ioffV );
  document.getElementById("debug_textarea").value +=
    "dgebak('N','R',n,ilo,ihi,scale,m,V,ldv,info) , V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 3 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 3 + 1 * ldv ]  + "\n";

  for ( j = 0; j < m; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  LaPack0.dgebak( 'P', 'R', n, ilo, ihi, scale, m, V, ldv, info,
    ioffscale,ioffV );
  document.getElementById("debug_textarea").value +=
    "dgebak('P','R',n,ilo,ihi,scale,m,V,ldv,info) , V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 3 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 3 + 1 * ldv ]  + "\n";

  for ( j = 0; j < m; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  LaPack0.dgebak( 'S', 'R', n, ilo, ihi, scale, m, V, ldv, info,
    ioffscale,ioffV );
  document.getElementById("debug_textarea").value +=
    "dgebak('S','R',n,ilo,ihi,scale,m,V,ldv,info) , V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 3 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 3 + 1 * ldv ]  + "\n";

  for ( j = 0; j < m; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  LaPack0.dgebak( 'B', 'R', n, ilo, ihi, scale, m, V, ldv, info,
    ioffscale,ioffV );
  document.getElementById("debug_textarea").value +=
    "dgebak('B','R',n,ilo,ihi,scale,m,V,ldv,info) , V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 3 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 3 + 1 * ldv ]  + "\n";

  for ( j = 0; j < m; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  LaPack0.dgebak( 'N', 'L', n, ilo, ihi, scale, m, V, ldv, info,
    ioffscale,ioffV );
  document.getElementById("debug_textarea").value +=
    "dgebak('N','L',n,ilo,ihi,scale,m,V,ldv,info) , V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 3 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 3 + 1 * ldv ]  + "\n";

  for ( j = 0; j < m; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  LaPack0.dgebak( 'P', 'L', n, ilo, ihi, scale, m, V, ldv, info,
    ioffscale,ioffV );
  document.getElementById("debug_textarea").value +=
    "dgebak('P','L',n,ilo,ihi,scale,m,V,ldv,info) , V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 3 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 3 + 1 * ldv ]  + "\n";

  for ( j = 0; j < m; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  LaPack0.dgebak( 'S', 'L', n, ilo, ihi, scale, m, V, ldv, info,
    ioffscale,ioffV );
  document.getElementById("debug_textarea").value +=
    "dgebak('S','L',n,ilo,ihi,scale,m,V,ldv,info) , V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 3 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 3 + 1 * ldv ]  + "\n";

  for ( j = 0; j < m; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  LaPack0.dgebak( 'B', 'L', n, ilo, ihi, scale, m, V, ldv, info,
    ioffscale,ioffV );
  document.getElementById("debug_textarea").value +=
    "dgebak('B','L',n,ilo,ihi,scale,m,V,ldv,info) , V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 3 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 3 + 1 * ldv ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dggbak() {
  document.getElementById("debug_textarea").value +=
    "testing dggbak *************" + "\n";
  var m = 4;
  var n = 3;
  var ldv = 5;
  var iofflscale = 1;
  var ioffrscale = 2;
  var ioffV = 3;
  var lscale = new Array( iofflscale + n );
  var rscale = new Array( ioffrscale + n );
  for ( var j = 0; j < n; j ++ ) {
    lscale[ iofflscale + j ] = Number.POSITIVE_INFINITY;
    lscale[ ioffrscale + j ] = Number.POSITIVE_INFINITY;
  }
  var V = new Array( ioffV + ldv * n );
  for ( j = 0; j < n; j ++ ) {
    for ( var i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  var ilo = 1;
  var ihi = n;
  var info = new IntReference();
  LaPack0.dggbak('N','R',n,ilo,ihi,lscale,rscale,m,V,ldv,info,
    iofflscale,ioffrscale,ioffV);
  document.getElementById("debug_textarea").value +=
    "dggbak('N','R',...), V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 0 + 2 * ldv ] + " "
    + V[ ioffV + 1 + 2 * ldv ] + " "
    + V[ ioffV + 2 + 2 * ldv ]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  for ( j = 0; j < n; j ++ ) {
    lscale[ iofflscale + j ] = 1 - j;
    rscale[ ioffrscale + j ] = 1 - 2 * j;
  }
  ilo = 1;
  ihi = n;
  LaPack0.dggbak('S','R',n,ilo,ihi,lscale,rscale,m,V,ldv,info,
    iofflscale,ioffrscale,ioffV);
  document.getElementById("debug_textarea").value +=
    "dggbak('S','R',...), V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 0 + 2 * ldv ] + " "
    + V[ ioffV + 1 + 2 * ldv ] + " "
    + V[ ioffV + 2 + 2 * ldv ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  ilo = 2;
  ihi = 2;
  LaPack0.dggbak('S','R',n,ilo,ihi,lscale,rscale,m,V,ldv,info,
    iofflscale,ioffrscale,ioffV);
  document.getElementById("debug_textarea").value +=
    "dggbak('S','R',...), V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 0 + 2 * ldv ] + " "
    + V[ ioffV + 1 + 2 * ldv ] + " "
    + V[ ioffV + 2 + 2 * ldv ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  ilo = 1;
  ihi = n;
  LaPack0.dggbak('S','L',n,ilo,ihi,lscale,rscale,m,V,ldv,info,
    iofflscale,ioffrscale,ioffV);
  document.getElementById("debug_textarea").value +=
    "dggbak('S','L',...), V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 0 + 2 * ldv ] + " "
    + V[ ioffV + 1 + 2 * ldv ] + " "
    + V[ ioffV + 2 + 2 * ldv ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  ilo = 2;
  ihi = 2;
  LaPack0.dggbak('S','L',n,ilo,ihi,lscale,rscale,m,V,ldv,info,
    iofflscale,ioffrscale,ioffV);
  document.getElementById("debug_textarea").value +=
    "dggbak('S','L',...), V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 0 + 2 * ldv ] + " "
    + V[ ioffV + 1 + 2 * ldv ] + " "
    + V[ ioffV + 2 + 2 * ldv ]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  lscale[ iofflscale + 0 ] = 1;
  lscale[ iofflscale + 1 ] = 3;
  lscale[ iofflscale + 2 ] = 2;
  rscale[ ioffrscale + 0 ] = 2;
  rscale[ ioffrscale + 1 ] = 1;
  rscale[ ioffrscale + 2 ] = 3;
  ilo = 1;
  ihi = n;
  LaPack0.dggbak('P','R',n,ilo,ihi,lscale,rscale,m,V,ldv,info,
    iofflscale,ioffrscale,ioffV);
  document.getElementById("debug_textarea").value +=
    "dggbak('P','R',...), V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 0 + 2 * ldv ] + " "
    + V[ ioffV + 1 + 2 * ldv ] + " "
    + V[ ioffV + 2 + 2 * ldv ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  ilo = 2;
  ihi = 2;
  LaPack0.dggbak('P','R',n,ilo,ihi,lscale,rscale,m,V,ldv,info,
    iofflscale,ioffrscale,ioffV);
  document.getElementById("debug_textarea").value +=
    "dggbak('P','R',...), V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 0 + 2 * ldv ] + " "
    + V[ ioffV + 1 + 2 * ldv ] + " "
    + V[ ioffV + 2 + 2 * ldv ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  ilo = 1;
  ihi = n;
  LaPack0.dggbak('P','L',n,ilo,ihi,lscale,rscale,m,V,ldv,info,
    iofflscale,ioffrscale,ioffV);
  document.getElementById("debug_textarea").value +=
    "dggbak('P','L',...), V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 0 + 2 * ldv ] + " "
    + V[ ioffV + 1 + 2 * ldv ] + " "
    + V[ ioffV + 2 + 2 * ldv ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  ilo = 2;
  ihi = 2;
  LaPack0.dggbak('P','L',n,ilo,ihi,lscale,rscale,m,V,ldv,info,
    iofflscale,ioffrscale,ioffV);
  document.getElementById("debug_textarea").value +=
    "dggbak('P','L',...), V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 0 + 2 * ldv ] + " "
    + V[ ioffV + 1 + 2 * ldv ] + " "
    + V[ ioffV + 2 + 2 * ldv ]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  lscale[ iofflscale + 0 ] = 1;
  lscale[ iofflscale + 1 ] = 3;
  lscale[ iofflscale + 2 ] = 2;
  rscale[ ioffrscale + 0 ] = 2;
  rscale[ ioffrscale + 1 ] = 1;
  rscale[ ioffrscale + 2 ] = 3;
  ilo = 1;
  ihi = n;
  LaPack0.dggbak('B','R',n,ilo,ihi,lscale,rscale,m,V,ldv,info,
    iofflscale,ioffrscale,ioffV);
  document.getElementById("debug_textarea").value +=
    "dggbak('B','R',...), V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 0 + 2 * ldv ] + " "
    + V[ ioffV + 1 + 2 * ldv ] + " "
    + V[ ioffV + 2 + 2 * ldv ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  ilo = 2;
  ihi = 2;
  LaPack0.dggbak('B','R',n,ilo,ihi,lscale,rscale,m,V,ldv,info,
    iofflscale,ioffrscale,ioffV);
  document.getElementById("debug_textarea").value +=
    "dggbak('B','R',...), V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 0 + 2 * ldv ] + " "
    + V[ ioffV + 1 + 2 * ldv ] + " "
    + V[ ioffV + 2 + 2 * ldv ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  ilo = 1;
  ihi = n;
  LaPack0.dggbak('B','L',n,ilo,ihi,lscale,rscale,m,V,ldv,info,
    iofflscale,ioffrscale,ioffV);
  document.getElementById("debug_textarea").value +=
    "dggbak('B','L',...), V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 0 + 2 * ldv ] + " "
    + V[ ioffV + 1 + 2 * ldv ] + " "
    + V[ ioffV + 2 + 2 * ldv ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      V[ ioffV + i + j * ldv ] = i + 2 * j;
    }
  }
  ilo = 2;
  ihi = 2;
  LaPack0.dggbak('B','L',n,ilo,ihi,lscale,rscale,m,V,ldv,info,
    iofflscale,ioffrscale,ioffV);
  document.getElementById("debug_textarea").value +=
    "dggbak('B','L',...), V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 0 + 2 * ldv ] + " "
    + V[ ioffV + 1 + 2 * ldv ] + " "
    + V[ ioffV + 2 + 2 * ldv ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgttrf() {
  document.getElementById("debug_textarea").value +=
    "testing dgttrf *************" + "\n";
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
  var info = new IntReference();
  for ( var i = 0; i < n; i ++ ) d[ ioffd + i ] = 2.;
  for ( i = 0; i < n - 1; i ++ ) {
    dl[ ioffdl + i ] = -1.;
    du[ ioffdu + i ] = -1.;
  }
  LaPack0.dgttrf( n, dl, d, du, du2, ipiv, info,
    ioffdl, ioffd, ioffdu, ioffdu2, ioffipiv );
  document.getElementById("debug_textarea").value +=
    "dgttrf, dl = "
    + dl[ ioffdl + 0 ] + " "
    + dl[ ioffdl + 1 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "dgttrf, d = "
    + d[ ioffd + 0 ] + " "
    + d[ ioffd + 1 ] + " "
    + d[ ioffd + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "dgttrf, du = "
    + du[ ioffdu + 0 ] + " "
    + du[ ioffdu + 1 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "dgttrf, du2 = "
    + du2[ ioffdu2 + 0 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "dgttrf, ipiv = "
    + ipiv[ ioffipiv + 0 ] + " "
    + ipiv[ ioffipiv + 1 ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgtts2() {
  document.getElementById("debug_textarea").value +=
    "testing dgtts2 *************" + "\n";
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
  var info = new IntReference();
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
  LaPack0.dgtts2( 0, n, nrhs, dl, d, du, du2, ipiv, B, ldb, 
    ioffdl, ioffd, ioffdu, ioffdu2, ioffipiv, ioffB );
  document.getElementById("debug_textarea").value +=
    "dgtts2, B = "
    + B[ ioffB + 0 ] + " "
    + B[ ioffB + 1 ] + " "
    + B[ ioffB + 2 ]  + "\n";

  B[ ioffB + 0 ] = 0.;
  B[ ioffB + 1 ] = 0.;
  B[ ioffB + 2 ] = 4.;
  LaPack0.dgtts2( 1, n, nrhs, dl, d, du, du2, ipiv, B, ldb, 
    ioffdl, ioffd, ioffdu, ioffdu2, ioffipiv, ioffB );
  document.getElementById("debug_textarea").value +=
    "dgtts2, B = "
    + B[ ioffB + 0 ] + " "
    + B[ ioffB + 1 ] + " "
    + B[ ioffB + 2 ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dla_gerpvgrw() {
  document.getElementById("debug_textarea").value +=
    "testing dla_gerpvgrw *************" + "\n";
  var n = 3;
  var ncols = 3;
  var lda = 4;
  var ldaf = 5;
  var ioffa = 1;
  var ioffaf = 2;
  var A = new Array( ioffa + lda * ncols );
  var AF = new Array( ioffaf + ldaf * ncols );
  for ( var j = 0; j < ncols; j ++ ) {
    for ( var i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
      AF[ ioffaf + i + j * ldaf ] = 3 * i - 4 * j + 2;
    }
  }
  var val =
    LaPack0.dla_gerpvgrw(n,ncols,A,lda,AF,ldaf,ioffa,ioffaf);
  document.getElementById("debug_textarea").value +=
    "dla_gerpvgrw = " + val + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dla_wwaddw() {
  document.getElementById("debug_textarea").value +=
    "testing dla_wwaddw *************" + "\n";
  var n=3;
  var ioffx = 1;
  var ioffy = 2;
  var ioffw = 3;
  var x = new Array( ioffx + n );
  var y = new Array( ioffy + n );
  var w = new Array( ioffw + n );
  for ( var i = 0; i < n; i ++ ) {
    x[ ioffx + i ] = i + 1;
    y[ ioffy + i ] = 2 * i - 1;
    w[ ioffw + i ] = 2 - 3 * i;
  }
  LaPack0.dla_wwaddw(n,x,y,w,ioffx,ioffy,ioffw);
  document.getElementById("debug_textarea").value +=
    "x = "
    + x[ ioffx + 0 ] + " "
    + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "y = "
    + y[ ioffy + 0 ] + " "
    + y[ ioffy + 1 ] + " "
    + y[ ioffy + 2 ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlacon() {
  document.getElementById("debug_textarea").value +=
    "testing dlacon *************" + "\n";
  var n = 4;
  var ioffA = 1;
  var A = new Array( ioffA + n * n );
  A[ ioffA + 0 + 0 * n ] = 1.;
  A[ ioffA + 1 + 0 * n ] = -1.;
  A[ ioffA + 2 + 0 * n ] = -1.;
  A[ ioffA + 3 + 0 * n ] = 1.;
  A[ ioffA + 0 + 1 * n ] = 1.;
  A[ ioffA + 1 + 1 * n ] = 1.;
  A[ ioffA + 2 + 1 * n ] = -1.;
  A[ ioffA + 3 + 1 * n ] = -1.;
  A[ ioffA + 0 + 2 * n ] = 1.;
  A[ ioffA + 1 + 2 * n ] = -1.;
  A[ ioffA + 2 + 2 * n ] = 1.;
  A[ ioffA + 3 + 2 * n ] = -1.;
  A[ ioffA + 0 + 3 * n ] = 1.;
  A[ ioffA + 1 + 3 * n ] = 1.;
  A[ ioffA + 2 + 3 * n ] = 1.;
  A[ ioffA + 3 + 3 * n ] = 1.;
  var ioffv = 2;
  var ioffx = 3;
  var ioffisgn = 4;
  var v = new Array( ioffv + n );
  var x = new Array( ioffx + n );
  var isgn = new Array( ioffisgn + n );
  var estReference = new NumberReference();
  var kaseReference = new IntReference(0);

  var ioffy = 5;
  var y = new Array( ioffy + n );
  while ( true ) {
    LaPack0.dlacon( n, v, x, isgn, estReference, kaseReference, ioffv,
      ioffx, ioffisgn );
    document.getElementById("debug_textarea").value +=
      "kase = " + kaseReference.getValue() + "\n";
    document.getElementById("debug_textarea").value +=
      "estReference = " + estReference.getValue() + "\n";
    document.getElementById("debug_textarea").value +=
      "x = " + x[ ioffx + 0 ] + " " + x[ ioffx + 1 ] + " "
      + x[ ioffx + 2 ] + " " + x[ ioffx + 3 ] + "\n";
    if ( kaseReference.getValue() == 0 ) break;
    else if ( kaseReference.getValue() == 1 ) {
      Blas2.dgemv( 'N', n, n, 1., A, n, x, 1, 0., y, 1,
        ioffA, ioffx, ioffy );
    } else if ( kaseReference.getValue() == 2 ) {
      Blas2.dgemv( 'T', n, n, 1., A, n, x, 1, 0., y, 1, 
        ioffA, ioffx, ioffy );
    }
    document.getElementById("debug_textarea").value +=
      "y = " + y[ ioffy + 0 ] + " " + y[ ioffy + 1 ] + " "
      + y[ ioffy + 2 ] + " " + y[ ioffy + 3 ] + "\n";
    Blas1.dcopy( n, y, 1, x, 1, ioffy, ioffx );
  }
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlacn2() {
  document.getElementById("debug_textarea").value +=
    "testing dlacn2 *************" + "\n";
  var n = 4;
  var ioffA = 1;
  var A = new Array( ioffA + n * n );
  A[ ioffA + 0 + 0 * n ] = 1.;
  A[ ioffA + 1 + 0 * n ] = -1.;
  A[ ioffA + 2 + 0 * n ] = -1.;
  A[ ioffA + 3 + 0 * n ] = 1.;
  A[ ioffA + 0 + 1 * n ] = 1.;
  A[ ioffA + 1 + 1 * n ] = 1.;
  A[ ioffA + 2 + 1 * n ] = -1.;
  A[ ioffA + 3 + 1 * n ] = -1.;
  A[ ioffA + 0 + 2 * n ] = 1.;
  A[ ioffA + 1 + 2 * n ] = -1.;
  A[ ioffA + 2 + 2 * n ] = 1.;
  A[ ioffA + 3 + 2 * n ] = -1.;
  A[ ioffA + 0 + 3 * n ] = 1.;
  A[ ioffA + 1 + 3 * n ] = 1.;
  A[ ioffA + 2 + 3 * n ] = 1.;
  A[ ioffA + 3 + 3 * n ] = 1.;
  var ioffv = 2;
  var ioffx = 3;
  var ioffisgn = 4;
  var ioffisave = 5;
  var v = new Array( ioffv + n );
  var x = new Array( ioffx + n );
  var isgn = new Array( ioffisgn + n );
  var isave = new Array( ioffisave + 3 );
  var estReference = new NumberReference();
  var kaseReference = new IntReference(0);

  var ioffy = 5;
  var y = new Array( ioffy + n );
  while ( true ) {
    LaPack0.dlacn2( n, v, x, isgn, estReference, kaseReference, isave,
      ioffv, ioffx, ioffisgn, ioffisave );
    document.getElementById("debug_textarea").value +=
      "kase = " + kaseReference.getValue() + "\n";
    document.getElementById("debug_textarea").value +=
      "estReference = " + estReference.getValue() + "\n";
    document.getElementById("debug_textarea").value +=
      "isave = "
      + isave[ ioffisave + 0 ] + " "
      + isave[ ioffisave + 1 ] + " "
      + isave[ ioffisave + 2 ]  + "\n";
    document.getElementById("debug_textarea").value +=
      "x = " + x[ ioffx + 0 ] + " " + x[ ioffx + 1 ] + " "
      + x[ ioffx + 2 ] + " " + x[ ioffx + 3 ] + "\n";
    if ( kaseReference.getValue() == 0 ) break;
    else if ( kaseReference.getValue() == 1 ) {
      Blas2.dgemv( 'N', n, n, 1., A, n, x, 1, 0., y, 1,
        ioffA, ioffx, ioffy );
    } else if ( kaseReference.getValue() == 2 ) {
      Blas2.dgemv( 'T', n, n, 1., A, n, x, 1, 0., y, 1, 
        ioffA, ioffx, ioffy );
    }
    document.getElementById("debug_textarea").value +=
      "y = " + y[ ioffy + 0 ] + " " + y[ ioffy + 1 ] + " "
      + y[ ioffy + 2 ] + " " + y[ ioffy + 3 ] + "\n";
    Blas1.dcopy( n, y, 1, x, 1, ioffy, ioffx );
  }
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlacpy() {
  document.getElementById("debug_textarea").value +=
    "testing dlacpy *************" + "\n";
  var m = 3;
  var n = 2;
  var ioffA = 1;
  var ioffB = 2;
  var lda = 4;
  var A = new Array( ioffA + lda * n );
  var ldb = 5;
  var B = new Array( ioffB + ldb * n );
  for ( var j = 0; j <= n - 1; j ++ ) {
    for ( var i = 0; i <= lda - 1; i ++ ) {
      A[ ioffA + i + j * lda ] = i + j;
    }
  }

  for ( j = 0; j <= n - 1; j ++ ) {
    for ( i = 0; i <= ldb - 1; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number.POSITIVE_INFINITY;
    }
  }
  LaPack0.dlacpy( 'U', m, n, A, lda, B, ldb, ioffA, ioffB );
  document.getElementById("debug_textarea").value +=
    "dlacpy( 'U', m, n, A, lda, B, ldb): B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ]  + "\n";

  for ( j = 0; j <= n - 1; j ++ ) {
    for ( i = 0; i <= ldb - 1; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number.POSITIVE_INFINITY;
    }
  }
  LaPack0.dlacpy( 'L', m, n, A, lda, B, ldb, ioffA, ioffB );
  document.getElementById("debug_textarea").value +=
    "dlacpy( 'L', m, n, A, lda, B, ldb): B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ]  + "\n";

  for ( j = 0; j <= n - 1; j ++ ) {
    for ( i = 0; i <= ldb - 1; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number.POSITIVE_INFINITY;
    }
  }
  LaPack0.dlacpy( 'A', m, n, A, lda, B, ldb, ioffA, ioffB );
  document.getElementById("debug_textarea").value +=
    "dlacpy( 'A', m, n, A, lda, B, ldb): B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dladiv() {
  document.getElementById("debug_textarea").value +=
    "testing dladiv *************" + "\n";
  var pReference = new NumberReference();
  var qReference = new NumberReference();
  LaPack0.dladiv( 1., 2., 3., 4., pReference, qReference );
  document.getElementById("debug_textarea").value +=
    "dladiv( 1., 2., 3., 4., p, q ): p = " + pReference.getValue()
    + " q = " + qReference.getValue()  + "\n";
  LaPack0.dladiv( 1., 2., 4., 3., pReference, qReference );
  document.getElementById("debug_textarea").value +=
    "dladiv( 1., 2., 4., 3., p, q ): p = " + pReference.getValue()
    + " q = " + qReference.getValue()  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlae2() {
  document.getElementById("debug_textarea").value +=
    "testing dlae2 *************" + "\n";
  var rt1Reference = new NumberReference();
  var rt2Reference = new NumberReference();
  LaPack0.dlae2( 1., 2., 3., rt1Reference, rt2Reference );
  document.getElementById("debug_textarea").value +=
    "dlae2(1.,2.,3.,rt1,rt2): rt1 = " + rt1Reference.getValue()
    + " rt2 = " + rt2Reference.getValue()  + "\n";
  LaPack0.dlae2( 3., 2., 1., rt1Reference, rt2Reference );
  document.getElementById("debug_textarea").value +=
    "dlae2(3.,2.,1.,rt1,rt2): rt1 = " + rt1Reference.getValue()
    + " rt2 = " + rt2Reference.getValue()  + "\n";
  LaPack0.dlae2( 1., 0.5, 3., rt1Reference, rt2Reference );
  document.getElementById("debug_textarea").value +=
    "dlae2(1.,0.5,3.,rt1,rt2): rt1 = " + rt1Reference.getValue()
    + " rt2 = " + rt2Reference.getValue()  + "\n";
  LaPack0.dlae2( 1., 2., 1., rt1Reference, rt2Reference );
  document.getElementById("debug_textarea").value +=
    "dlae2(1.,2.,1.,rt1,rt2): rt1 = " + rt1Reference.getValue()
    + " rt2 = " + rt2Reference.getValue()  + "\n";
  LaPack0.dlae2( -1., -2., -3., rt1Reference, rt2Reference );
  document.getElementById("debug_textarea").value +=
    "dlae2(-1.,-2.,-3.,rt1,rt2): rt1 = " + rt1Reference.getValue()
    + " rt2 = " + rt2Reference.getValue()  + "\n";
  LaPack0.dlae2( 1., 2., -1., rt1Reference, rt2Reference );
  document.getElementById("debug_textarea").value +=
    "dlae2(1.,2.,-1.,rt1,rt2): rt1 = " + rt1Reference.getValue()
    + " rt2 = " + rt2Reference.getValue()  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlaebz() {
  document.getElementById("debug_textarea").value +=
    "testing dlaebz *************" + "\n";
  var nitmax = 10;
  var n = 4;
  var mmax = 3;
  var minp = 2;
  var nbmin = 0;
  var abstol = 0.;
  var reltol = 1.e-5;
  var pivmin = 1.e-5;
  var ioffd = 1;
  var ioffe = 2;
  var ioffe2 = 3;
  var ioffnval = 4;
  var ioffAB = 5;
  var ioffc = 6;
  var ioffNAB = 7;
  var ioffwork = 8;
  var ioffiwork = 9;
  var d = new Array( ioffd + n );
  var e = new Array( ioffe + n );
  var e2 = new Array( ioffe2 + n );
  var nval = new Array( ioffnval + minp );
  var AB = new Array( ioffAB + mmax * 2 );
  var c = new Array( ioffc + mmax );
  var mout = new IntReference(-1);
  var NAB = new Array( ioffNAB + mmax * 2 );
  var work = new Array( ioffwork + mmax );
  var iwork = new Array( ioffiwork + mmax );
  var info = new IntReference(0);
  for ( i=0; i<=minp-1; i++ ) nval[ ioffnval + i ] = -1; 
  for ( var i=0; i<=n-1; i++ ) d[ ioffd + i ] = 1.;
  var t = 0.1;
  for ( i=0; i<=n-2; i++ ) { 
    e[ ioffe + i ] = t; 
    t *= 0.1; 
    e2[ ioffe2 + i ] = e[ ioffe + i ] * e[ ioffe + i ];
  }
  for ( i=0; i<=mmax-1; i++ ) {
    c[ ioffc + i ] = 1.;
    work[ ioffwork + i ] = Number.POSITIVE_INFINITY;
    iwork[ ioffiwork + i ] = -1;
  }
  AB[ ioffAB + 0 + 0 * mmax ] = .8;
  AB[ ioffAB + 0 + 1 * mmax ] = .9;
  AB[ ioffAB + 1 + 0 * mmax ] = .9;
  AB[ ioffAB + 1 + 1 * mmax ] = 1.1;
  AB[ ioffAB + 2 + 0 * mmax ] = 1.1;
  AB[ ioffAB + 2 + 1 * mmax ] = 1.2;
  NAB[ ioffNAB + 0 + 0 * mmax ] = -1;
  NAB[ ioffNAB + 0 + 1 * mmax ] = -1;
  NAB[ ioffNAB + 1 + 0 * mmax ] = -1;
  NAB[ ioffNAB + 1 + 1 * mmax ] = -1;
  NAB[ ioffNAB + 2 + 0 * mmax ] = -1;
  NAB[ ioffNAB + 2 + 1 * mmax ] = -1;
  LaPack0.dlaebz(1,nitmax,n,mmax,minp,nbmin,abstol,reltol,pivmin,d,
    e,e2, nval,AB,c,mout,NAB,work,iwork,info,
    ioffd,ioffe,ioffe2,ioffnval,ioffAB,ioffc,ioffNAB,ioffwork,
    ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dlaebz(1,...): mout = " + mout.getValue()
    + " info = " + info.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "NAB = "
    + NAB[ ioffNAB + 0 + 0 * mmax ] + " "
    + NAB[ ioffNAB + 1 + 0 * mmax ] + " "
    + NAB[ ioffNAB + 2 + 0 * mmax ] + " "
    + NAB[ ioffNAB + 0 + 1 * mmax ] + " "
    + NAB[ ioffNAB + 1 + 1 * mmax ] + " "
    + NAB[ ioffNAB + 2 + 1 * mmax ]  + "\n";

  LaPack0.dlaebz(2,nitmax,n,mmax,minp,nbmin,abstol,reltol,pivmin,d,
    e,e2, nval,AB,c,mout,NAB,work,iwork,info,
    ioffd,ioffe,ioffe2,ioffnval,ioffAB,ioffc,ioffNAB,ioffwork,
    ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dlaebz(2,...): mout = " + mout.getValue()
    + " info = " + info.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "AB = "
    + AB[ ioffAB + 0 + 0 * mmax ] + " "
    + AB[ ioffAB + 1 + 0 * mmax ] + " "
    + AB[ ioffAB + 2 + 0 * mmax ] + " "
    + AB[ ioffAB + 0 + 1 * mmax ] + " "
    + AB[ ioffAB + 1 + 1 * mmax ] + " "
    + AB[ ioffAB + 2 + 1 * mmax ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "NAB = "
    + NAB[ ioffNAB + 0 + 0 * mmax ] + " "
    + NAB[ ioffNAB + 1 + 0 * mmax ] + " "
    + NAB[ ioffNAB + 2 + 0 * mmax ] + " "
    + NAB[ ioffNAB + 0 + 1 * mmax ] + " "
    + NAB[ ioffNAB + 1 + 1 * mmax ] + " "
    + NAB[ ioffNAB + 2 + 1 * mmax ]  + "\n";

  for ( i=0; i<=mmax-1; i++ ) c[ ioffc + i ] = 1.;
  AB[ ioffAB + 0 + 0 * mmax ] = .8;
  AB[ ioffAB + 0 + 1 * mmax ] = .9;
  AB[ ioffAB + 1 + 0 * mmax ] = .9;
  AB[ ioffAB + 1 + 1 * mmax ] = 1.1;
  AB[ ioffAB + 2 + 0 * mmax ] = 1.1;
  AB[ ioffAB + 2 + 1 * mmax ] = 1.2;
  LaPack0.dlaebz(1,nitmax,n,mmax,minp,nbmin,abstol,reltol,pivmin,d,
    e,e2, nval,AB,c,mout,NAB,work,iwork,info,
    ioffd,ioffe,ioffe2,ioffnval,ioffAB,ioffc,ioffNAB,ioffwork,
    ioffiwork);
  LaPack0.dlaebz(3,nitmax,n,mmax,minp,nbmin,abstol,reltol,pivmin,d,
    e,e2, nval,AB,c,mout,NAB,work,iwork,info,
    ioffd,ioffe,ioffe2,ioffnval,ioffAB,ioffc,ioffNAB,ioffwork,
    ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dlaebz(3,...): mout = " + mout.getValue()
    + " info = " + info.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "AB = "
    + AB[ ioffAB + 0 + 0 * mmax ] + " "
    + AB[ ioffAB + 1 + 0 * mmax ] + " "
    + AB[ ioffAB + 2 + 0 * mmax ] + " "
    + AB[ ioffAB + 0 + 1 * mmax ] + " "
    + AB[ ioffAB + 1 + 1 * mmax ] + " "
    + AB[ ioffAB + 2 + 1 * mmax ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "NAB = "
    + NAB[ ioffNAB + 0 + 0 * mmax ] + " "
    + NAB[ ioffNAB + 1 + 0 * mmax ] + " "
    + NAB[ ioffNAB + 2 + 0 * mmax ] + " "
    + NAB[ ioffNAB + 0 + 1 * mmax ] + " "
    + NAB[ ioffNAB + 1 + 1 * mmax ] + " "
    + NAB[ ioffNAB + 2 + 1 * mmax ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "nval = " + nval[ ioffnval + 0 ] + " " + nval[ ioffnval + 1 ] + "\n";
  document.getElementById("debug_textarea").value +=
    "c = " + c[ ioffc + 0 ] + " " + c[ ioffc + 1 ] + " "
    + c[ ioffc + 2 ] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlaed5() {
  document.getElementById("debug_textarea").value +=
    "testing dlaed5 *************" + "\n";
  var ioffd = 1;
  var ioffz = 2;
  var ioffdelta = 3;
  var d = new Array(ioffd + 2);
  var z = new Array(ioffz + 2);
  var delta = new Array(ioffdelta + 2);
  var dlamReference = new NumberReference();
  d[ioffd + 0]=2.;
  d[ioffd + 1]=3.;
  z[ioffz + 0]=4.;
  z[ioffz + 1]=5.;
  LaPack0.dlaed5(1,d,z,delta,1.,dlamReference,ioffd,ioffz,ioffdelta);
  document.getElementById("debug_textarea").value +=
    "dlaed5(1,d,z,delta,1.,dlam), dlam = " + dlamReference.getValue()
    + "\n";
  document.getElementById("debug_textarea").value +=
    "delta = " + delta[ioffdelta + 0] + " "
    + delta[ioffdelta + 1] + "\n";

  z[ioffz + 0]=5.;
  z[ioffz + 1]=4.;
  LaPack0.dlaed5(1,d,z,delta,1.,dlamReference,ioffd,ioffz,ioffdelta);
  document.getElementById("debug_textarea").value +=
    "dlaed5(1,d,z,delta,1.,dlam), dlam = " + dlamReference.getValue()
    + "\n";
  document.getElementById("debug_textarea").value +=
    "delta = " + delta[ioffdelta + 0] + " "
    + delta[ioffdelta + 1] + "\n";

  z[ioffz + 0]=4.;
  z[ioffz + 1]=5.;
  LaPack0.dlaed5(2,d,z,delta,1.,dlamReference,ioffd,ioffz,ioffdelta);
  document.getElementById("debug_textarea").value +=
    "dlaed5(2,d,z,delta,1.,dlam), dlam = " + dlamReference.getValue()
    + "\n";
  document.getElementById("debug_textarea").value +=
    "delta = " + delta[ioffdelta + 0] + " "
    + delta[ioffdelta + 1] + "\n";

  LaPack0.dlaed5(2,d,z,delta,.01,dlamReference,ioffd,ioffz,ioffdelta);
  document.getElementById("debug_textarea").value +=
    "dlaed5(2,d,z,delta,.01,dlam), dlam = " + dlamReference.getValue()
    + "\n";
  document.getElementById("debug_textarea").value +=
    "delta = " + delta[ioffdelta + 0] + " "
    + delta[ioffdelta + 1] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlaeda() {
  document.getElementById("debug_textarea").value +=
    "testing dlaeda *************" + "\n";
  var lnn = 2;
  var n = Math.pow(2,lnn);
  var tlvls = lnn;
  var curlvl = 1;
  var curpbm = 1;
  var prmptr = new Array( n * lnn );
  var perm = new Array( n * lnn );
  var givptr = new Array( n * lnn );
  var Givcol = new Array( 2 * n * lnn );
  var Givnum = new Array( 2 * n * lnn );
  var Q = new Array( n * n );
  var qptr = new Array( n + 2 );
  var z = new Array( n );
  var ztemp = new Array( n );
  var info = new IntReference();
//    not sure how to set up input
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlaev2() {
  document.getElementById("debug_textarea").value +=
    "testing dlaev2 *************" + "\n";
  var rt1Reference = new NumberReference();
  var rt2Reference = new NumberReference();
  var cs1Reference = new NumberReference();
  var sn1Reference = new NumberReference();
  LaPack0.dlaev2(3.,1.,2.,rt1Reference,rt2Reference,cs1Reference,
    sn1Reference);
  document.getElementById("debug_textarea").value +=
    "dlaev2(3.,1.,2.,...): rt1,rt2,cs1,sn1 = "
    + rt1Reference.getValue() + " " + rt2Reference.getValue() + " "
    + cs1Reference.getValue() + " " + sn1Reference.getValue() + "\n";
  LaPack0.dlaev2(-3.,-1.,-2.,rt1Reference,rt2Reference,cs1Reference,
    sn1Reference);
  document.getElementById("debug_textarea").value +=
    "dlaev2(-3.,-1.,-2.,...): rt1,rt2,cs1,sn1 = "
    + rt1Reference.getValue() + " " + rt2Reference.getValue() + " "
    + cs1Reference.getValue() + " " + sn1Reference.getValue() + "\n";
  LaPack0.dlaev2(2.,1.,2.,rt1Reference,rt2Reference,cs1Reference,
    sn1Reference);
  document.getElementById("debug_textarea").value +=
    "dlaev2(2.,1.,2.,...): rt1,rt2,cs1,sn1 = "
    + rt1Reference.getValue() + " " + rt2Reference.getValue() + " "
    + cs1Reference.getValue() + " " + sn1Reference.getValue() + "\n";
  LaPack0.dlaev2(-2.,-1.,-2.,rt1Reference,rt2Reference,cs1Reference,
    sn1Reference);
  document.getElementById("debug_textarea").value +=
    "dlaev2(-2.,-1.,-2.,...): rt1,rt2,cs1,sn1 = "
    + rt1Reference.getValue() + " " + rt2Reference.getValue() + " "
    + cs1Reference.getValue() + " " + sn1Reference.getValue() + "\n";
  LaPack0.dlaev2(5.,1.,2.,rt1Reference,rt2Reference,cs1Reference,
    sn1Reference);
  document.getElementById("debug_textarea").value +=
    "dlaev2(5.,1.,2.,...): rt1,rt2,cs1,sn1 = "
    + rt1Reference.getValue() + " " + rt2Reference.getValue() + " "
    + cs1Reference.getValue() + " " + sn1Reference.getValue() + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlag2() {
  document.getElementById("debug_textarea").value +=
    "testing dlag2 *************" + "\n";
  var lda = 3;
  var ioffA = 1;
  var ioffB = 2;
  var A = new Array( ioffA + lda * 2 );
  var ldb = 4;
  var B = new Array( ioffB + ldb * 2 );
  var safmin = Number.MIN_VALUE;
  var scale1Reference = new NumberReference();
  var scale2Reference = new NumberReference();
  var wr1Reference = new NumberReference();
  var wr2Reference = new NumberReference();
  var wiReference = new NumberReference();
  A[ ioffA + 0 + 0 * lda ] = 2.;
  A[ ioffA + 1 + 0 * lda ] = 1.;
  A[ ioffA + 0 + 1 * lda ] = -1.;
  A[ ioffA + 1 + 1 * lda ] = 3.;
  B[ ioffB + 0 + 0 * ldb ] = 2.;
  B[ ioffB + 0 + 1 * ldb ] = -1.;
  B[ ioffB + 1 + 1 * ldb ] = 2.;
  LaPack0.dlag2( A, lda, B, ldb, safmin, scale1Reference, scale2Reference,
    wr1Reference, wr2Reference, wiReference, ioffA, ioffB );
  document.getElementById("debug_textarea").value +=
    "dlag2( A, lda, B, ldb ...): wr1,wr2,wi = " + wr1Reference.getValue()
    + " " + wr2Reference.getValue() + " " + wiReference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "scale1,scale2 = " + scale1Reference.getValue() + " "
    + scale2Reference.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "leaving dlag2 *************" + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlagtm() {
  document.getElementById("debug_textarea").value +=
    "testing dlagtm *************" + "\n";
  var n = 3;
  var nrhs = 1;
  var ioffdl = 1;
  var ioffd = 2;
  var ioffdu = 3;
  var ioffX = 4;
  var ioffB = 5;
  var dl = new Array( ioffdl + n - 1 );
  var d = new Array( ioffd + n );
  var du = new Array( ioffdu + n - 1 );
  var ldx = 4;
  var X = new Array( ioffX + ldx * nrhs );
  var ldb = 3;
  var B = new Array( ioffB + ldb * nrhs );

  d[ ioffd + 0 ] = 1.;
  d[ ioffd + 1 ] = 2.;
  d[ ioffd + 2 ] = 3.;
  dl[ ioffdl + 0 ] = -1.;
  dl[ ioffdl + 1 ] = -2.;
  du[ ioffdu + 0 ] = 3.;
  du[ ioffdu + 1 ] = 4.;
  X[ ioffX + 0 + 0 * ldx ] = 2.;
  X[ ioffX + 1 + 0 * ldx ] = 3.;
  X[ ioffX + 2 + 0 * ldx ] = 4.;

  B[ ioffB + 0 + 0 * ldx ] = 2.;
  B[ ioffB + 1 + 0 * ldx ] = 1.;
  B[ ioffB + 2 + 0 * ldx ] = 0.;
  LaPack0.dlagtm('N',n,nrhs,0.,dl,d,du,X,ldx,0.,B,ldb,
    ioffdl,ioffd,ioffdu,ioffX,ioffB);
  document.getElementById("debug_textarea").value +=
    "LaPack0.dlagtm('N',n,nrhs,0.,dl,d,du,X,ldx,0.,B,ldb): B = "
    + B[ ioffB + 0 + 0 * ldx ] + " "
    + B[ ioffB + 1 + 0 * ldx ] + " "
    + B[ ioffB + 2 + 0 * ldx ] + "\n"; 

  B[ ioffB + 0 + 0 * ldx ] = 2.;
  B[ ioffB + 1 + 0 * ldx ] = 1.;
  B[ ioffB + 2 + 0 * ldx ] = 0.;
  LaPack0.dlagtm('N',n,nrhs,0.,dl,d,du,X,ldx,2.,B,ldb,
    ioffdl,ioffd,ioffdu,ioffX,ioffB);
  document.getElementById("debug_textarea").value +=
    "LaPack0.dlagtm('N',n,nrhs,0.,dl,d,du,X,ldx,2.,B,ldb): B = "
    + B[ ioffB + 0 + 0 * ldx ] + " "
    + B[ ioffB + 1 + 0 * ldx ] + " "
    + B[ ioffB + 2 + 0 * ldx ] + "\n"; 

  B[ ioffB + 0 + 0 * ldx ] = 2.;
  B[ ioffB + 1 + 0 * ldx ] = 1.;
  B[ ioffB + 2 + 0 * ldx ] = 0.;
  LaPack0.dlagtm('N',n,nrhs,1.,dl,d,du,X,ldx,0.,B,ldb,
    ioffdl,ioffd,ioffdu,ioffX,ioffB);
  document.getElementById("debug_textarea").value +=
    "LaPack0.dlagtm('N',n,nrhs,1.,dl,d,du,X,ldx,0.,B,ldb): B = "
    + B[ ioffB + 0 + 0 * ldx ] + " "
    + B[ ioffB + 1 + 0 * ldx ] + " "
    + B[ ioffB + 2 + 0 * ldx ] + "\n"; 

  B[ ioffB + 0 + 0 * ldx ] = 2.;
  B[ ioffB + 1 + 0 * ldx ] = 1.;
  B[ ioffB + 2 + 0 * ldx ] = 0.;
  LaPack0.dlagtm('N',n,nrhs,1.,dl,d,du,X,ldx,2.,B,ldb,
    ioffdl,ioffd,ioffdu,ioffX,ioffB);
  document.getElementById("debug_textarea").value +=
    "LaPack0.dlagtm('N',n,nrhs,1.,dl,d,du,X,ldx,2.,B,ldb): B = "
    + B[ ioffB + 0 + 0 * ldx ] + " "
    + B[ ioffB + 1 + 0 * ldx ] + " "
    + B[ ioffB + 2 + 0 * ldx ] + "\n"; 

  B[ ioffB + 0 + 0 * ldx ] = 2.;
  B[ ioffB + 1 + 0 * ldx ] = 1.;
  B[ ioffB + 2 + 0 * ldx ] = 0.;
  LaPack0.dlagtm('N',n,nrhs,-1.,dl,d,du,X,ldx,0.,B,ldb,
    ioffdl,ioffd,ioffdu,ioffX,ioffB);
  document.getElementById("debug_textarea").value +=
    "LaPack0.dlagtm('N',n,nrhs,-1.,dl,d,du,X,ldx,0.,B,ldb): B = "
    + B[ ioffB + 0 + 0 * ldx ] + " "
    + B[ ioffB + 1 + 0 * ldx ] + " "
    + B[ ioffB + 2 + 0 * ldx ] + "\n"; 

  B[ ioffB + 0 + 0 * ldx ] = 2.;
  B[ ioffB + 1 + 0 * ldx ] = 1.;
  B[ ioffB + 2 + 0 * ldx ] = 0.;
  LaPack0.dlagtm('N',n,nrhs,-1.,dl,d,du,X,ldx,2.,B,ldb,
    ioffdl,ioffd,ioffdu,ioffX,ioffB);
  document.getElementById("debug_textarea").value +=
    "LaPack0.dlagtm('N',n,nrhs,-1.,dl,d,du,X,ldx,2.,B,ldb): B = "
    + B[ ioffB + 0 + 0 * ldx ] + " "
    + B[ ioffB + 1 + 0 * ldx ] + " "
    + B[ ioffB + 2 + 0 * ldx ] + "\n"; 

  B[ ioffB + 0 + 0 * ldx ] = 2.;
  B[ ioffB + 1 + 0 * ldx ] = 1.;
  B[ ioffB + 2 + 0 * ldx ] = 0.;
  LaPack0.dlagtm('T',n,nrhs,0.,dl,d,du,X,ldx,0.,B,ldb,
    ioffdl,ioffd,ioffdu,ioffX,ioffB);
  document.getElementById("debug_textarea").value +=
    "LaPack0.dlagtm('T',n,nrhs,0.,dl,d,du,X,ldx,0.,B,ldb): B = "
    + B[ ioffB + 0 + 0 * ldx ] + " "
    + B[ ioffB + 1 + 0 * ldx ] + " "
    + B[ ioffB + 2 + 0 * ldx ] + "\n"; 

  B[ ioffB + 0 + 0 * ldx ] = 2.;
  B[ ioffB + 1 + 0 * ldx ] = 1.;
  B[ ioffB + 2 + 0 * ldx ] = 0.;
  LaPack0.dlagtm('T',n,nrhs,0.,dl,d,du,X,ldx,2.,B,ldb,
    ioffdl,ioffd,ioffdu,ioffX,ioffB);
  document.getElementById("debug_textarea").value +=
    "LaPack0.dlagtm('T',n,nrhs,0.,dl,d,du,X,ldx,2.,B,ldb): B = "
    + B[ ioffB + 0 + 0 * ldx ] + " "
    + B[ ioffB + 1 + 0 * ldx ] + " "
    + B[ ioffB + 2 + 0 * ldx ] + "\n"; 

  B[ ioffB + 0 + 0 * ldx ] = 2.;
  B[ ioffB + 1 + 0 * ldx ] = 1.;
  B[ ioffB + 2 + 0 * ldx ] = 0.;
  LaPack0.dlagtm('T',n,nrhs,1.,dl,d,du,X,ldx,0.,B,ldb,
    ioffdl,ioffd,ioffdu,ioffX,ioffB);
  document.getElementById("debug_textarea").value +=
    "LaPack0.dlagtm('T',n,nrhs,1.,dl,d,du,X,ldx,0.,B,ldb): B = "
    + B[ ioffB + 0 + 0 * ldx ] + " "
    + B[ ioffB + 1 + 0 * ldx ] + " "
    + B[ ioffB + 2 + 0 * ldx ] + "\n"; 

  B[ ioffB + 0 + 0 * ldx ] = 2.;
  B[ ioffB + 1 + 0 * ldx ] = 1.;
  B[ ioffB + 2 + 0 * ldx ] = 0.;
  LaPack0.dlagtm('T',n,nrhs,1.,dl,d,du,X,ldx,2.,B,ldb,
    ioffdl,ioffd,ioffdu,ioffX,ioffB);
  document.getElementById("debug_textarea").value +=
    "LaPack0.dlagtm('T',n,nrhs,1.,dl,d,du,X,ldx,2.,B,ldb): B = "
    + B[ ioffB + 0 + 0 * ldx ] + " "
    + B[ ioffB + 1 + 0 * ldx ] + " "
    + B[ ioffB + 2 + 0 * ldx ] + "\n"; 

  B[ ioffB + 0 + 0 * ldx ] = 2.;
  B[ ioffB + 1 + 0 * ldx ] = 1.;
  B[ ioffB + 2 + 0 * ldx ] = 0.;
  LaPack0.dlagtm('T',n,nrhs,-1.,dl,d,du,X,ldx,0.,B,ldb,
    ioffdl,ioffd,ioffdu,ioffX,ioffB);
  document.getElementById("debug_textarea").value +=
    "LaPack0.dlagtm('T',n,nrhs,-1.,dl,d,du,X,ldx,0.,B,ldb): B = "
    + B[ ioffB + 0 + 0 * ldx ] + " "
    + B[ ioffB + 1 + 0 * ldx ] + " "
    + B[ ioffB + 2 + 0 * ldx ] + "\n"; 

  B[ ioffB + 0 + 0 * ldx ] = 2.;
  B[ ioffB + 1 + 0 * ldx ] = 1.;
  B[ ioffB + 2 + 0 * ldx ] = 0.;
  LaPack0.dlagtm('T',n,nrhs,-1.,dl,d,du,X,ldx,2.,B,ldb,
    ioffdl,ioffd,ioffdu,ioffX,ioffB);
  document.getElementById("debug_textarea").value +=
    "LaPack0.dlagtm('T',n,nrhs,-1.,dl,d,du,X,ldx,2.,B,ldb): B = "
    + B[ ioffB + 0 + 0 * ldx ] + " "
    + B[ ioffB + 1 + 0 * ldx ] + " "
    + B[ ioffB + 2 + 0 * ldx ] + "\n"; 
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlamch() {
  document.getElementById("debug_textarea").value +=
    "testing dlamch *************" + "\n";
  var ret = LaPack0.dlamch('E');
  document.getElementById("debug_textarea").value +=
    "dlamch('E') = " + ret + "\n";
  ret = LaPack0.dlamch('S');
  document.getElementById("debug_textarea").value +=
    "dlamch('S') = " + ret + "\n";
  ret = LaPack0.dlamch('B');
  document.getElementById("debug_textarea").value +=
    "dlamch('B') = " + ret + "\n";
  ret = LaPack0.dlamch('P');
  document.getElementById("debug_textarea").value +=
    "dlamch('P') = " + ret + "\n";
  ret = LaPack0.dlamch('N');
  document.getElementById("debug_textarea").value +=
    "dlamch('N') = " + ret + "\n";
  ret = LaPack0.dlamch('R');
  document.getElementById("debug_textarea").value +=
    "dlamch('R') = " + ret + "\n";
  ret = LaPack0.dlamch('M');
  document.getElementById("debug_textarea").value +=
    "dlamch('M') = " + ret + "\n";
  ret = LaPack0.dlamch('U');
  document.getElementById("debug_textarea").value +=
    "dlamch('U') = " + ret + "\n";
  ret = LaPack0.dlamch('L');
  document.getElementById("debug_textarea").value +=
    "dlamch('L') = " + ret + "\n";
  ret = LaPack0.dlamch('O');
  document.getElementById("debug_textarea").value +=
    "dlamch('O') = " + ret + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlamrg() {
  document.getElementById("debug_textarea").value +=
    "testing dlamrg *************" + "\n";
  var ioffa = 1;
  var ioffindex = 2;
  var a = new Array(ioffa + 5);
  var index = new Array(5);

  a[ioffa + 0]=1.;
  a[ioffa + 1]=2.;
  a[ioffa + 2]=3.;
  a[ioffa + 3]=-1.;
  a[ioffa + 4]=4.;
  LaPack0.dlamrg(3,2,a,1,1,index,ioffa,ioffindex);
  document.getElementById("debug_textarea").value +=
    "dlamrg(3,2,a,1,1,index): index = "
    + index[ioffindex + 0] + " "
    + index[ioffindex + 1] + " "
    + index[ioffindex + 2] + " "
    + index[ioffindex + 3] + " "
    + index[ioffindex + 4] + "\n";

  a[ioffa + 0]=1.;
  a[ioffa + 1]=2.;
  a[ioffa + 2]=-2.;
  a[ioffa + 3]=-1.;
  a[ioffa + 4]=4.;
  LaPack0.dlamrg(2,3,a,1,1,index,ioffa,ioffindex);
  document.getElementById("debug_textarea").value +=
    "dlamrg(2,3,a,1,1,index): index = "
    + index[ioffindex + 0] + " "
    + index[ioffindex + 1] + " "
    + index[ioffindex + 2] + " "
    + index[ioffindex + 3] + " "
    + index[ioffindex + 4] + "\n";

  a[ioffa + 0]=-1.;
  a[ioffa + 1]=-2.;
  a[ioffa + 2]=-3.;
  a[ioffa + 3]=1.;
  a[ioffa + 4]=-4.;
  LaPack0.dlamrg(3,2,a,-1,-1,index,ioffa,ioffindex);
  document.getElementById("debug_textarea").value +=
    "dlamrg(3,2,a,-1,-1,index): index = "
    + index[ioffindex + 0] + " "
    + index[ioffindex + 1] + " "
    + index[ioffindex + 2] + " "
    + index[ioffindex + 3] + " "
    + index[ioffindex + 4] + "\n";

  a[ioffa + 0]=-1.;
  a[ioffa + 1]=-2.;
  a[ioffa + 2]=2.;
  a[ioffa + 3]=1.;
  a[ioffa + 4]=-4.;
  LaPack0.dlamrg(2,3,a,-1,-1,index,ioffa,ioffindex);
  document.getElementById("debug_textarea").value +=
    "dlamrg(2,3,a,-1,-1,index): index = "
    + index[ioffindex + 0] + " "
    + index[ioffindex + 1] + " "
    + index[ioffindex + 2] + " "
    + index[ioffindex + 3] + " "
    + index[ioffindex + 4] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlapmr() {
  document.getElementById("debug_textarea").value +=
    "testing dlapmr *************" + "\n";
  var m = 4;
  var n = 3;
  var ldx = 5;
  var ioffX = 1;
  var ioffk = 2;
  var X = new Array( ioffX + ldx * n );
  var k = new Array( ioffk + m );
  k[ ioffk + 0 ] = 2;
  k[ ioffk + 1 ] = 3;
  k[ ioffk + 2 ] = 4;
  k[ ioffk + 3 ] = 1;

  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      X[ ioffX + i + j * ldx ] = i + j;
    }
  }
  LaPack0.dlapmr(true,m,n,X,ldx,k,ioffX,ioffk);
  document.getElementById("debug_textarea").value +=
    "dlapmr(true,m,n,X,ldx,k): X = "
    + X[ ioffX + 0 + 0 * ldx ] + " "
    + X[ ioffX + 1 + 0 * ldx ] + " "
    + X[ ioffX + 2 + 0 * ldx ] + " "
    + X[ ioffX + 3 + 0 * ldx ] + " "
    + X[ ioffX + 0 + 1 * ldx ] + " "
    + X[ ioffX + 1 + 1 * ldx ] + " "
    + X[ ioffX + 2 + 1 * ldx ] + " "
    + X[ ioffX + 3 + 1 * ldx ] + " "
    + X[ ioffX + 0 + 2 * ldx ] + " "
    + X[ ioffX + 1 + 2 * ldx ] + " "
    + X[ ioffX + 2 + 2 * ldx ] + " "
    + X[ ioffX + 3 + 2 * ldx ]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) X[ ioffX + i + j * ldx ] = i + j;
  }
  LaPack0.dlapmr(false,m,n,X,ldx,k,ioffX,ioffk);
  document.getElementById("debug_textarea").value +=
    "dlapmr(false,m,n,X,ldx,k): X = "
    + X[ ioffX + 0 + 0 * ldx ] + " "
    + X[ ioffX + 1 + 0 * ldx ] + " "
    + X[ ioffX + 2 + 0 * ldx ] + " "
    + X[ ioffX + 3 + 0 * ldx ] + " "
    + X[ ioffX + 0 + 1 * ldx ] + " "
    + X[ ioffX + 1 + 1 * ldx ] + " "
    + X[ ioffX + 2 + 1 * ldx ] + " "
    + X[ ioffX + 3 + 1 * ldx ] + " "
    + X[ ioffX + 0 + 2 * ldx ] + " "
    + X[ ioffX + 1 + 2 * ldx ] + " "
    + X[ ioffX + 2 + 2 * ldx ] + " "
    + X[ ioffX + 3 + 2 * ldx ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlapmt() {
  document.getElementById("debug_textarea").value +=
    "testing dlapmt *************" + "\n";
  var m = 2;
  var n = 3;
  var ldx = 3;
  var ioffX = 1;
  var ioffk = 2;
  var X = new Array( ioffX + ldx * n );
  var k = new Array( ioffk + n );
  k[ ioffk + 0 ] = 2;
  k[ ioffk + 1 ] = 3;
  k[ ioffk + 2 ] = 1;

  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      X[ ioffX + i + j * ldx ] = i + j;
    }
  }
  LaPack0.dlapmt(true,m,n,X,ldx,k,ioffX,ioffk);
  document.getElementById("debug_textarea").value +=
    "dlapmt(true,m,n,X,ldx,k): X = "
    + X[ ioffX + 0 + 0 * ldx ] + " "
    + X[ ioffX + 1 + 0 * ldx ] + " "
    + X[ ioffX + 0 + 1 * ldx ] + " "
    + X[ ioffX + 1 + 1 * ldx ] + " "
    + X[ ioffX + 0 + 2 * ldx ] + " "
    + X[ ioffX + 1 + 2 * ldx ]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) X[ ioffX + i + j * ldx ] = i + j;
  }
  LaPack0.dlapmt(false,m,n,X,ldx,k,ioffX,ioffk);
  document.getElementById("debug_textarea").value +=
    "dlapmt(false,m,n,X,ldx,k): X = "
    + X[ ioffX + 0 + 0 * ldx ] + " "
    + X[ ioffX + 1 + 0 * ldx ] + " "
    + X[ ioffX + 0 + 1 * ldx ] + " "
    + X[ ioffX + 1 + 1 * ldx ] + " "
    + X[ ioffX + 0 + 2 * ldx ] + " "
    + X[ ioffX + 1 + 2 * ldx ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlapy2() {
  document.getElementById("debug_textarea").value +=
    "testing dlapy2 *************" + "\n";
  var x = 3. * Math.sqrt( Number.MAX_VALUE );
  var y = 4. * Math.sqrt( Number.MAX_VALUE );
  var z = LaPack0.dlapy2( x, y );
  document.getElementById("debug_textarea").value +=
    "dlapy2( x, y ) = " + z  + "\n";

//    x = .6 * Math.sqrt( Number.MIN_VALUE ); // sub-normal
//    y = .8 * Math.sqrt( Number.MIN_VALUE );
//    z = LaPack0.dlapy2( x, y );
//    document.getElementById("debug_textarea").value +=
//      "dlapy2( x, y ) = " + z  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlapy3() {
  document.getElementById("debug_textarea").value +=
    "testing dlapy3 *************" + "\n";
  var x = 3. * Math.sqrt( Number.MAX_VALUE );
  var y = 4. * Math.sqrt( Number.MAX_VALUE );
  var z = 12. * Math.sqrt( Number.MAX_VALUE );
  var w = LaPack0.dlapy3( x, y, z );
  document.getElementById("debug_textarea").value +=
    "dlapy3( x, y, z ) = " + w  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlaqr1() {
  document.getElementById("debug_textarea").value +=
    "testing dlaqr1 *************" + "\n";
  var ldh = 4;
  var ioffh = 1;
  var ioffv = 2;
  var H = new Array( ioffh + ldh * 3 );
  var v = new Array( ioffv + 3 );

  var n = 2;
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < n; i ++ ) {
      H[ ioffh + i + j * ldh ] = 1 + i + 2 * j;
    }
  }
  LaPack0.dlaqr1(n,H,ldh,1.,2.,3.,4.,v,ioffh,ioffv);
  document.getElementById("debug_textarea").value +=
    "v = "
    + v[ ioffv + 0 ] + " "
    + v[ ioffv + 1 ]  + "\n";

  n = 3;
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      H[ ioffh + i + j * ldh ] = 1 + i + 2 * j;
    }
  }
  LaPack0.dlaqr1(n,H,ldh,1.,2.,3.,4.,v,ioffh,ioffv);
  document.getElementById("debug_textarea").value +=
    "v = "
    + v[ ioffv + 0 ] + " "
    + v[ ioffv + 1 ] + " "
    + v[ ioffv + 2 ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlar2v() {
  document.getElementById("debug_textarea").value +=
    "testing dlar2v *************" + "\n";
  var n = 2;
  var incx = 2;
  var ioffx = 1;
  var ioffy = 2;
  var ioffz = 3;
  var ioffc = 4;
  var ioffs = 5;
  var x = new Array( ioffx + 1 + ( n - 1 ) * incx );
  var y = new Array( ioffy + 1 + ( n - 1 ) * incx );
  var z = new Array( ioffz + 1 + ( n - 1 ) * incx );
  var incc = 3;
  var c = new Array( ioffc + 1 + ( n - 1 ) * incc );
  var s = new Array( ioffs + 1 + ( n - 1 ) * incc );
  for ( var i = 0; i <= ( n - 1 ) * incx; i ++ ) {
    x[ ioffx + i ] = Number( i );
    y[ ioffy + i ] = Number( 2*i+1 );
    z[ ioffz + i ] = Number( 1-2*i );
  }
  for ( i = 0; i <= ( n - 1 ) * incc; i ++ ) {
    c[ ioffc + i ] = Number( i );
    s[ ioffs + i ] = Number( 2*i+1 );
    var t = LaPack0.dlapy2( c[ioffc + i], s[ioffs + i] );
    c[ ioffc + i ] /= t;
    s[ ioffs + i ] /= t;
  }
  LaPack0.dlar2v(n,x,y,z,incx,c,s,incc,ioffx,ioffy,ioffz,ioffc,ioffs);
  document.getElementById("debug_textarea").value +=
    "dlar2v: x = " + x[ ioffx + 0 ] + " " + x[ ioffx + 1 ] + " "
    + x[ ioffx + 2 ] + "\n";
  document.getElementById("debug_textarea").value +=
    "y = " + y[ ioffy + 0 ] + " " + y[ ioffy + 1 ] + " "
    + y[ ioffy + 2 ] + "\n";
  document.getElementById("debug_textarea").value +=
    "z = " + z[ ioffz + 0 ] + " " + z[ ioffz + 1 ] + " "
    + z[ ioffz + 2 ] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlarft() {
  document.getElementById("debug_textarea").value +=
    "testing dlarft *************" + "\n";
  var k = 3;
  var ldt = 4;
  var iofftau = 1;
  var ioffT = 2;
  var ioffV = 3;
  var tau = new Array( iofftau + k );
  for ( var i = 0; i < k; i ++ ) tau[ iofftau + i ] = i + 1;
  var T = new Array( ioffT + k * ldt ); // kxk
  var n = 5;
  var ldv = 6;

//    direct='F', storev='C'
//    V = [V1] kxk unit lower triangular
//        [V2] (n-k)xk
  var V = new Array( ioffV + k * ldv );
  for ( var j = 0; j < k; j ++ ) {
    for ( i = 0; i < k; i ++ ) {
      T[ ioffT + i + j * ldt ] = 0.;
    }
  }
  for ( j = 0; j < k; j ++ ) {
    for ( i = j+1; i < n; i ++ ) V[ ioffV + i + j * ldv ] = 2*i + 3*j;
  }
  LaPack0.dlarft('F','C',n,k,V,ldv,tau,T,ldt,ioffV,iofftau,ioffT);
  document.getElementById("debug_textarea").value +=
    "dlarft('F','C',...): T = "
    + T[ ioffT + 0 + 0 * ldt ] + " "
    + T[ ioffT + 0 + 1 * ldt ] + " "
    + T[ ioffT + 1 + 1 * ldt ] + " "
    + T[ ioffT + 0 + 2 * ldt ] + " "
    + T[ ioffT + 1 + 2 * ldt ] + " "
    + T[ ioffT + 2 + 2 * ldt ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "V = "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 3 + 0 * ldv ] + " "
    + V[ ioffV + 4 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 3 + 1 * ldv ] + " "
    + V[ ioffV + 4 + 1 * ldv ] + " "
    + V[ ioffV + 3 + 2 * ldv ] + " "
    + V[ ioffV + 4 + 2 * ldv ]  + "\n";

//    direct='F', storev='R'
//    V = [V1,V2] kxk unit upper triangular, kx(n-k)
  V = new Array( k * ldv );
  for ( j = 0; j < k; j ++ ) {
    for ( i = 0; i < k; i ++ ) {
      T[ ioffT + i + j * ldt ] = 0.;
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < Math.min(j,k); i ++ ) {
      V[ ioffV + i + j * ldv ] = 2*i + 3*j;
    }
  }
  LaPack0.dlarft('F','R',n,k,V,ldv,tau,T,ldt,ioffV,iofftau,ioffT);
  document.getElementById("debug_textarea").value +=
    "dlarft('F','R',...): T = "
    + T[ ioffT + 0 + 0 * ldt ] + " "
    + T[ ioffT + 1 + 0 * ldt ] + " "
    + T[ ioffT + 2 + 0 * ldt ] + " "
    + T[ ioffT + 0 + 1 * ldt ] + " "
    + T[ ioffT + 1 + 1 * ldt ] + " "
    + T[ ioffT + 2 + 1 * ldt ] + " "
    + T[ ioffT + 0 + 2 * ldt ] + " "
    + T[ ioffT + 1 + 2 * ldt ] + " "
    + T[ ioffT + 2 + 2 * ldt ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "V = "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 0 + 2 * ldv ] + " "
    + V[ ioffV + 1 + 2 * ldv ] + " "
    + V[ ioffV + 0 + 3 * ldv ] + " "
    + V[ ioffV + 1 + 3 * ldv ] + " "
    + V[ ioffV + 2 + 3 * ldv ] + " "
    + V[ ioffV + 0 + 4 * ldv ] + " "
    + V[ ioffV + 1 + 4 * ldv ] + " "
    + V[ ioffV + 2 + 4 * ldv ]  + "\n";

//    direct='B', storev='C'
//    V = [V1] (n-k)xk 
//        [V2] kxk unit upper triangular
  V = new Array( ioffV + k * ldv );
  for ( j = 0; j < k; j ++ ) {
    for ( i = 0; i < k; i ++ ) {
      T[ ioffT + i + j * ldt ] = 0.;
    }
  }
  for ( j = 0; j < k; j ++ ) {
    for ( i = 0; i < n+j-k; i ++ ) {
      V[ ioffV + i + j * ldv ] = 2*i + 3*j;
    }
  }
  LaPack0.dlarft('B','C',n,k,V,ldv,tau,T,ldt,ioffV,iofftau,ioffT);
  document.getElementById("debug_textarea").value +=
    "dlarft('B','C',...): T = "
    + T[ ioffT + 0 + 0 * ldt ] + " "
    + T[ ioffT + 1 + 0 * ldt ] + " "
    + T[ ioffT + 2 + 0 * ldt ] + " "
    + T[ ioffT + 1 + 1 * ldt ] + " "
    + T[ ioffT + 2 + 1 * ldt ] + " "
    + T[ ioffT + 2 + 2 * ldt ]  + "\n";
  document.getElementById("debug_textarea").value +=
     "V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 0 + 2 * ldv ] + " "
    + V[ ioffV + 1 + 2 * ldv ] + " "
    + V[ ioffV + 2 + 2 * ldv ] + " "
    + V[ ioffV + 3 + 2 * ldv ] + "\n";

//    direct='B', storev='R'
//    V = [V1,V2] kx(n-k), kxk unit lower triangular
  V = new Array( ioffV + k * ldv );
  for ( j = 0; j < k; j ++ ) {
    for ( i = 0; i < k; i ++ ) {
      T[ ioffT + i + j * ldt ] = 0.;
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = Math.max(0,k+j+1-n); i < k; i ++ ) {
      V[ ioffV + i + j * ldv ] = 2*i + 3*j;
    }
  }
  LaPack0.dlarft('B','R',n,k,V,ldv,tau,T,ldt,ioffV,iofftau,ioffT);
  document.getElementById("debug_textarea").value +=
    "dlarft('B','R',...): T = "
    + T[ ioffT + 0 + 0 * ldt ] + " "
    + T[ ioffT + 1 + 0 * ldt ] + " "
    + T[ ioffT + 2 + 0 * ldt ] + " "
    + T[ ioffT + 1 + 1 * ldt ] + " "
    + T[ ioffT + 2 + 1 * ldt ] + " "
    + T[ ioffT + 2 + 2 * ldt ] + "\n";
  document.getElementById("debug_textarea").value +=
    "V = "
    + V[ ioffV + 0 + 0 * ldv ] + " "
    + V[ ioffV + 1 + 0 * ldv ] + " "
    + V[ ioffV + 2 + 0 * ldv ] + " "
    + V[ ioffV + 0 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 1 * ldv ] + " "
    + V[ ioffV + 2 + 1 * ldv ] + " "
    + V[ ioffV + 1 + 2 * ldv ] + " "
    + V[ ioffV + 2 + 2 * ldv ] + " "
    + V[ ioffV + 2 + 3 * ldv ] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlargv() {
  document.getElementById("debug_textarea").value +=
    "testing dlargv *************"+ "\n";
  var n = 4;
  var incx = 1;
  var ioffx = 1;
  var x = new Array( ioffx + 1 + ( n - 1 ) * incx );
  var incy = 2;
  var ioffy = 2;
  var y = new Array( ioffy + 1 + ( n - 1 ) * incy );
  var incc = 3;
  var ioffc = 3;
  var c = new Array( ioffc + 1 + ( n - 1 ) * incc );
  for ( var i = 0; i < 1 + ( n - 1 ) * incx; i ++ ) {
    x[ ioffx + i ] = i;
  }
  x[ ioffx + 2 ] = 3.;
  for ( i = 0; i < 1 + ( n - 1 ) * incy; i ++ ) {
    y[ ioffy + i ] = 2 * ( i - 1 );
  }
  for ( i = 0; i < 1 + ( n - 1 ) * incc; i ++ ) {
    c[ ioffc + i ] = Number.POSITIVE_INFINITY;
  }
  LaPack0.dlargv(n,x,incx,y,incy,c,incc,ioffx,ioffy,ioffc);
  document.getElementById("debug_textarea").value +=
    "dlargv: x = " + x[ ioffx + 0*incx ] + " "
    + x[ ioffx + 1*incx ] + " " + x[ ioffx + 2*incx ] + " "
    + x[ ioffx + 3*incx ] + "\n";
  document.getElementById("debug_textarea").value +=
    "y = " + y[ ioffy + 0*incy ] + " " + y[ ioffy + 1*incy ] + " "
    + y[ ioffy + 2*incy ] + " " + y[ ioffy + 3*incy ] + "\n";
  document.getElementById("debug_textarea").value +=
    "c = " + c[ ioffc + 0*incc ] + " " + c[ ioffc + 1*incc ] + " "
    + c[ ioffc + 2*incc ] + " " + c[ ioffc + 3*incc ] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlarra() {
  document.getElementById("debug_textarea").value +=
    "testing dlarra *************"+ "\n";
  var n = 3;
  var ioffd = 1;
  var ioffe = 2;
  var ioffe2 = 3;
  var ioffisplit = 4;
  var d = new Array( ioffd + n );
  var e = new Array( ioffe + n );
  var e2 = new Array( ioffe2 + n );
  var isplit = new Array( ioffisplit + n );
  for ( var i = 0; i < n; i ++ ) {
    d[ ioffd + i ] = 2.;
    e[ ioffe + i ] = -1.;
    e2[ ioffe2 + i ] = 1.;
  }
  e[ ioffe + 1 ] = -1.e-4;
  e2[ ioffe2 + 1 ] = 1.e-8;
  var nsplit = new IntReference();
  var info = new IntReference();
  LaPack0.dlarra( n, d, e, e2, 1.e-2, 4., nsplit, isplit, info,
    ioffd, ioffe, ioffe2, ioffisplit );
  document.getElementById("debug_textarea").value +=
    "dlarra: nsplit = " + nsplit.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "isplit = "
    + isplit[ ioffisplit + 0 ] + " "
    + isplit[ ioffisplit + 1 ] + " "
    + isplit[ ioffisplit + 2 ] + "\n";
  document.getElementById("debug_textarea").value +=
    "e = "
    + e[ ioffe + 0 ] + " "
    + e[ ioffe + 1 ] + " "
    + e[ ioffe + 2 ] + "\n";
  document.getElementById("debug_textarea").value +=
    "e2 = "
    + e2[ ioffe2 + 0 ] + " "
    + e2[ ioffe2 + 1 ] + " "
    + e2[ ioffe2 + 2 ] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlarrc() {
  document.getElementById("debug_textarea").value +=
    "testing dlarrc *************"+ "\n";
  var n = 64;
  var vl = 1.;
  var vu = 3.;
  var ioffd = 1;
  var ioffe = 2;
  var d = new Array( ioffd + n );
  var e = new Array( ioffe + n );
  var pivmin = Number.POSITIVE_INFINITY;
  var eigcnt = new IntReference();
  var lcnt = new IntReference();
  var rcnt = new IntReference();
  var info = new IntReference();
  for ( var i = 0; i < n; i ++ ) {
    d[ ioffd + i ] = 2.;
    e[ ioffd + i ] = -1.;
  }
  LaPack0.dlarrc('T',n,vl,vu,d,e,pivmin,eigcnt,lcnt,rcnt,info,
    ioffd,ioffe);
  document.getElementById("debug_textarea").value +=
    "dlarrc('T',...): eigcnt,lcnt,rcnt,info = "
    + eigcnt.getValue() + " " + lcnt.getValue() + " " + rcnt.getValue() + " "
    + info.getValue()+ "\n";
  LaPack0.dlarrc('L',n,vl,vu,d,e,pivmin,eigcnt,lcnt,rcnt,info,
    ioffd,ioffe);
  document.getElementById("debug_textarea").value +=
    "dlarrc('L',...): eigcnt,lcnt,rcnt,info = "
    + eigcnt.getValue() + " " + lcnt.getValue() + " " + rcnt.getValue() + " "
    + info.getValue() + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlarrj() {
  document.getElementById("debug_textarea").value +=
    "testing dlarrj *************" + "\n";
  var n = 64;
  var ioffd = 1;
  var ioffe2 = 2;
  var ioffw = 3;
  var ioffwerr = 4;
  var ioffwork = 5;
  var ioffiwork = 6;
  var d = new Array( ioffd + n );
  var e2 = new Array( ioffe2 + n );
  var w = new Array( ioffw + n );
  var werr = new Array( ioffwerr + n );
  var work = new Array( ioffwork + 2 * n );
  var iwork = new Array( ioffiwork + 2 * n );
  for ( var i = 0; i < n; i ++ ) {
    d[ ioffd + i ] = 2.;
    e2[ ioffe2 + i ] = 1.;
    w[ ioffw + i ] = 2. + Number( i - 32 ) / 64.;
    werr[ ioffwerr + i ] = 1.;
  }
  var spdiam = 4.;
  var pivmin = 1.e-14;
  var ifirst = 32;
  var ilast = 37;
  var rtol = 1.e-12;
  var offset = 1;
  var info = new IntReference();
  LaPack0.dlarrj(n,d,e2,ifirst,ilast,rtol,offset,w,werr,work,iwork,
    pivmin,spdiam,info,ioffd,ioffe2,ioffw,ioffwerr,ioffwork,ioffiwork);
  document.getElementById("debug_textarea").value +=
    "dlarrj: ,info = " + info.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "w = "
    + w[ ioffw + 28 ] + " "
    + w[ ioffw + 29 ] + " "
    + w[ ioffw + 30 ] + " "
    + w[ ioffw + 31 ] + " "
    + w[ ioffw + 32 ] + " "
    + w[ ioffw + 33 ] + " "
    + w[ ioffw + 34 ] + " "
    + w[ ioffw + 35 ] + " "
    + w[ ioffw + 36 ] + " "
    + w[ ioffw + 37 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "werr = "
    + werr[ ioffwerr + 28 ] + " "
    + werr[ ioffwerr + 29 ] + " "
    + werr[ ioffwerr + 30 ] + " "
    + werr[ ioffwerr + 31 ] + " "
    + werr[ ioffwerr + 32 ] + " "
    + werr[ ioffwerr + 33 ] + " "
    + werr[ ioffwerr + 34 ] + " "
    + werr[ ioffwerr + 35 ] + " "
    + werr[ ioffwerr + 36 ] + " "
    + werr[ ioffwerr + 37 ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlarscl2() {
  document.getElementById("debug_textarea").value +=
    "testing dlarscl2 *************" + "\n";
  var n = 3;
  var m = 4;
  var ldx = 5;
  var ioffd = 1;
  var ioffx = 2;
  var d = new Array( ioffd + m );
  var X = new Array( ioffx + ldx * n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      X[ ioffx + i + j * ldx ] = 1 + i + 2 * j;
    }
  }
  for ( i = 0; i < m; i ++ ) {
    d[ ioffd + i ] = 1 - 2 * i;
  }
  LaPack0.dlarscl2(m,n,d,X,ldx,ioffd,ioffx);
  document.getElementById("debug_textarea").value +=
    "X = "
    + X[ ioffx + 0 + 0 * ldx ] + " "
    + X[ ioffx + 1 + 0 * ldx ] + " "
    + X[ ioffx + 2 + 0 * ldx ] + " "
    + X[ ioffx + 3 + 0 * ldx ] + " "
    + X[ ioffx + 0 + 1 * ldx ] + " "
    + X[ ioffx + 1 + 1 * ldx ] + " "
    + X[ ioffx + 2 + 1 * ldx ] + " "
    + X[ ioffx + 3 + 1 * ldx ] + " "
    + X[ ioffx + 0 + 2 * ldx ] + " "
    + X[ ioffx + 1 + 2 * ldx ] + " "
    + X[ ioffx + 2 + 2 * ldx ] + " "
    + X[ ioffx + 3 + 2 * ldx ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlartv() {
  document.getElementById("debug_textarea").value +=
    "testing dlartv *************" + "\n";
  var n=3;
  var incx=1;
  var incy=2;
  var incc=3;
  var ioffx=1;
  var ioffy=2;
  var ioffc=3;
  var ioffs=4;
  var x = new Array( ioffx + 1 + ( n - 1 ) * incx );
  var y = new Array( ioffy + 1 + ( n - 1 ) * incy );
  var c = new Array( ioffc + 1 + ( n - 1 ) * incc );
  var s = new Array( ioffs + 1 + ( n - 1 ) * incc );
  for ( var i = 0; i < n; i ++ ) {
    x[ ioffx + i * incx ] = i;
    y[ ioffy + i * incy ] = 2*(i-1);
    c[ ioffc + i * incc ] = 2*i-1;
    s[ ioffs + i * incc ] = 1-3*i;
  }
  LaPack0.dlartv(n,x,incx,y,incy,c,s,incc,ioffx,ioffy,ioffc,ioffs);
  document.getElementById("debug_textarea").value +=
    "dlartv: x = "
    + x[ ioffx + 0 * incx ] + " "
    + x[ ioffx + 1 * incx ] + " "
    + x[ ioffx + 2 * incx ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "y = "
    + y[ ioffy + 0 * incy ] + " "
    + y[ ioffy + 1 * incy ] + " "
    + y[ ioffy + 2 * incy ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "c = "
    + c[ ioffc + 0 * incc ] + " "
    + c[ ioffc + 1 * incc ] + " "
    + c[ ioffc + 2 * incc ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "s = "
    + s[ ioffs + 0 * incc ] + " "
    + s[ ioffs + 1 * incc ] + " "
    + s[ ioffs + 2 * incc ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlarz() {
  document.getElementById("debug_textarea").value +=
    "testing dlarz *************" + "\n";
  var m = 4;
  var n = 3;
  var l = 2;
  var ldc = 5;
  var incv = 2;
  var ioffc = 1;
  var ioffv = 2;
  var ioffwork = 4;
  var C = new Array( ioffc + ldc * n );
  var v = new Array( ioffv + 1 + ( l - 1 ) * incv );
  var work = new Array( ioffwork + Math.max( m, n ) );
  var tau = 4.;
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = i + 2 * j;
    }
  }
  for ( i = 0; i <= 2; i ++ ) v[ ioffv + i ] = -1;
  v[ ioffv + 0 ] = 1;
  v[ ioffv + 2 ] = -5;
  LaPack0.dlarz('L',m,n,l,v,incv,tau,C,ldc,work,ioffv,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlarz: C = "
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
      C[ ioffc + i + j * ldc ] = i + 2 * j;
    }
  }
  LaPack0.dlarz('R',m,n,l,v,incv,tau,C,ldc,work,ioffv,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlarz: C = "
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
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlarzb() {
  document.getElementById("debug_textarea").value +=
    "testing dlarzb *************" + "\n";
  var ldc = 5;
  var ldt = 6;
  var ldv = 7;
  var ldwork = 8;
  var m = 4;
  var n = 3;
  var k = 2;
  var l = 2;
  var ioffc = 1;
  var iofft = 2;
  var ioffv = 3;
  var ioffwork = 4;
  var T = new Array( iofft + k * k );
  for ( var j = 0; j < k; j ++ ) {
    for ( var i = j; i < k; i ++ ) {
      T[ iofft + i + j * ldt ] = 2*i-3*j;
    }
  }
  var V = new Array( ioffv + ldv * k );
  for ( j = 0; j < k; j ++ ) {
    for ( i = 0; i < ldv; i ++ ) {
      V[ ioffv + i + j * ldv ] = 3 * i + 2 * j;
    }
  }
  var tau = 3.;
  var C = new Array( ioffc + ldc * n );
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffc + i + j * ldc ] = i + 2 * j;
    }
  }
  var Work = new Array( ioffwork + ldwork * k );
  LaPack0.dlarzb('L','N','B','R',m,n,k,l,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffv,iofft,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlarzb('L','N',...): C = "
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
      C[ ioffc + i + j * ldc ] = i + 2 * j;
    }
  }
  LaPack0.dlarzb('L','C','B','R',m,n,k,l,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffv,iofft,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlarzb('L','C',...): C = "
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
      C[ ioffc + i + j * ldc ] = i + 2 * j;
    }
  }
  LaPack0.dlarzb('R','N','B','R',m,n,k,l,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffv,iofft,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlarzb('R','N',...): C = "
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
      C[ ioffc + i + j * ldc ] = i + 2 * j;
    }
  }
  LaPack0.dlarzb('R','C','B','R',m,n,k,l,V,ldv,T,ldt,C,ldc,Work,ldwork,
    ioffv,iofft,ioffc,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlarzb('R','C',...): C = "
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
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlarzt() {
  document.getElementById("debug_textarea").value +=
    "testing dlarzt *************" + "\n";
  var n = 3;
  var k = 2;
  var ldt = 3;
  var ldv = 4;
  var iofft = 1;
  var iofftau = 2;
  var ioffv = 3;
  var V = new Array( ioffv + ldv * n );
  var tau = new Array( iofftau + k );
  var T = new Array( iofft + ldt * k );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < k; i ++ ) {
      V[ ioffv + i + j * ldv ] = 3 * i + 2 * j;
    }
  }
  for ( i = 0; i < k; i ++ ) {
    tau[ iofftau + i ] = i + 1;
  }
  for ( j = 0; j < k; j ++ ) {
    for ( i = 0; i < k; i ++ ) {
      T[ iofft + i + j * ldt ] = 0.;
    }
  }
  LaPack0.dlarzt('B','R',n,k,V,ldv,tau,T,ldt,ioffv,iofftau,iofft);
  document.getElementById("debug_textarea").value +=
    "dlarzt: T = "
    + T[ iofft + 0 + 0 * ldt ] + " "
    + T[ iofft + 1 + 0 * ldt ] + " "
    + T[ iofft + 1 + 1 * ldt ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlas2() {
  document.getElementById("debug_textarea").value +=
    "testing dlas2 *************" + "\n";
  var ssminReference = new NumberReference();
  var ssmaxReference = new NumberReference();
  LaPack0.dlas2(1.,1.,0.,ssminReference,ssmaxReference);
  document.getElementById("debug_textarea").value +=
    "dlas2(1.,1.,0.,...): ssmin,ssmax = "
    + ssminReference.getValue() + " "
    + ssmaxReference.getValue()  + "\n";
  LaPack0.dlas2(1.,0.,2.,ssminReference,ssmaxReference);
  document.getElementById("debug_textarea").value +=
    "dlas2(1.,0.,2.,...): ssmin,ssmax = "
    + ssminReference.getValue() + " "
    + ssmaxReference.getValue()  + "\n";
  var eps=1.e-200;
  LaPack0.dlas2(eps,1./eps,eps,ssminReference,ssmaxReference);
  document.getElementById("debug_textarea").value +=
    "dlas2(eps,1./eps,eps,...): ssmin,ssmax = "
    + ssminReference.getValue() + " "
    + ssmaxReference.getValue()  + "\n";
  LaPack0.dlas2(1.,3.,2.,ssminReference,ssmaxReference);
  document.getElementById("debug_textarea").value +=
    "dlas2(1.,3.,2.,...): ssmin,ssmax = "
    + ssminReference.getValue() + " "
    + ssmaxReference.getValue()  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlascl2() {
  document.getElementById("debug_textarea").value +=
    "testing dlascl2 *************" + "\n";
  var n = 3;
  var m = 4;
  var ldx = 5;
  var ioffd = 1;
  var ioffx = 2;
  var d = new Array( ioffd + m );
  var X = new Array( ioffx + ldx * n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      X[ ioffx + i + j * ldx ] = 1 + i + 2 * j;
    }
  }
  for ( i = 0; i < m; i ++ ) {
    d[ ioffd + i ] = 1 - 2 * i;
  }
  LaPack0.dlascl2(m,n,d,X,ldx,ioffd,ioffx);
  document.getElementById("debug_textarea").value +=
    "X = "
    + X[ ioffx + 0 + 0 * ldx ] + " "
    + X[ ioffx + 1 + 0 * ldx ] + " "
    + X[ ioffx + 2 + 0 * ldx ] + " "
    + X[ ioffx + 3 + 0 * ldx ] + " "
    + X[ ioffx + 0 + 1 * ldx ] + " "
    + X[ ioffx + 1 + 1 * ldx ] + " "
    + X[ ioffx + 2 + 1 * ldx ] + " "
    + X[ ioffx + 3 + 1 * ldx ] + " "
    + X[ ioffx + 0 + 2 * ldx ] + " "
    + X[ ioffx + 1 + 2 * ldx ] + " "
    + X[ ioffx + 2 + 2 * ldx ] + " "
    + X[ ioffx + 3 + 2 * ldx ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlasd5() {
  document.getElementById("debug_textarea").value +=
    "testing dlasd5 *************" + "\n";
  var ioffd = 1;
  var ioffz = 2;
  var ioffdelta = 3;
  var ioffwork = 4;
  var d = new Array( ioffd + 2 );
  var z = new Array( ioffz + 2 );
  var delta = new Array( ioffdelta + 2 );
  var work = new Array( ioffwork + 2 );
  var rho = 1.;
  d[ ioffd ] = 1.;
  d[ ioffd + 1 ] = 2.;
  z[ ioffz ] = 4.;
  z[ ioffz + 1 ] = 3.;
  var dsigmaReference = new NumberReference();
  LaPack0.dlasd5(1,d,z,delta,rho,dsigmaReference,work,
    ioffd,ioffz,ioffdelta,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlasd5(1,...): dsigma,delta,work = "
    + dsigmaReference.getValue() + " "
    + delta[ ioffdelta ] + " "
    + delta[ ioffdelta + 1 ] + " "
    + work[ ioffwork ] + " "
    + work[ ioffwork + 1 ]  + "\n";

  z[ ioffz ] = 3.;
  z[ ioffz + 1 ] = 4.;
  LaPack0.dlasd5(1,d,z,delta,rho,dsigmaReference,work,
    ioffd,ioffz,ioffdelta,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlasd5(1,...): dsigma,delta,work = "
    + dsigmaReference.getValue() + " "
    + delta[ ioffdelta ] + " "
    + delta[ ioffdelta + 1 ] + " "
    + work[ ioffwork ] + " "
    + work[ ioffwork + 1 ]  + "\n";

  LaPack0.dlasd5(2,d,z,delta,rho,dsigmaReference,work,
    ioffd,ioffz,ioffdelta,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlasd5(2,...): dsigma,delta,work = "
    + dsigmaReference.getValue() + " "
    + delta[ ioffdelta ] + " "
    + delta[ ioffdelta + 1 ] + " "
    + work[ ioffwork ] + " "
    + work[ ioffwork + 1 ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlasdt() {
  document.getElementById("debug_textarea").value +=
    "testing dlasdt *************" + "\n";
  var n = 8;
  var msub = 3;
  var ioffinode = 1;
  var ioffndiml = 2;
  var ioffndimr = 3;
  var inode = new Array( ioffinode + n );
  var ndiml = new Array( ioffndiml + n );
  var ndimr = new Array( ioffndimr + n );
  var lvl = new IntReference();
  var nd = new IntReference();
  LaPack0.dlasdt(n,lvl,nd,inode,ndiml,ndimr,msub,
    ioffinode,ioffndiml,ioffndimr);
  document.getElementById("debug_textarea").value +=
    "dlasdt: lvl,nd = " + lvl.getValue() + " " + nd.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "inode = "
    + inode[ ioffinode + 0 ] + " "
    + inode[ ioffinode + 1 ] + " "
    + inode[ ioffinode + 2 ] + " "
    + inode[ ioffinode + 3 ] + " "
    + inode[ ioffinode + 4 ] + " "
    + inode[ ioffinode + 5 ] + " "
    + inode[ ioffinode + 6 ] + " "
    + inode[ ioffinode + 7 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "ndiml = "
    + ndiml[ ioffndiml + 0 ] + " "
    + ndiml[ ioffndiml + 1 ] + " "
    + ndiml[ ioffndiml + 2 ] + " "
    + ndiml[ ioffndiml + 3 ] + " "
    + ndiml[ ioffndiml + 4 ] + " "
    + ndiml[ ioffndiml + 5 ] + " "
    + ndiml[ ioffndiml + 6 ] + " "
    + ndiml[ ioffndiml + 7 ]  + "\n";
  document.getElementById("debug_textarea").value +=
     "ndimr = "
    + ndimr[ ioffndimr + 0 ] + " "
    + ndimr[ ioffndimr + 1 ] + " "
    + ndimr[ ioffndimr + 2 ] + " "
    + ndimr[ ioffndimr + 3 ] + " "
    + ndimr[ ioffndimr + 4 ] + " "
    + ndimr[ ioffndimr + 5 ] + " "
    + ndimr[ ioffndimr + 6 ] + " "
    + ndimr[ ioffndimr + 7 ] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlaset() {
  document.getElementById("debug_textarea").value +=
    "testing dlaset *************" + "\n";
  var m=4;
  var n=3;
  var lda=5;
  var ldb=4;
  var ioffA=1;
  var ioffB=2;
  var A = new Array( ioffA + n * lda );
  var B = new Array( ioffB + m * ldb );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffA + i + j * lda ] = Number.POSITIVE_INFINITY;
    }
  }
  LaPack0.dlaset('U',m,n,0.,1.,A,lda,ioffA);
  document.getElementById("debug_textarea").value +=
    "dlaset('U',...): A = "
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

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffA + i + j * lda ] = Number.POSITIVE_INFINITY;
    }
  }
  LaPack0.dlaset('L',m,n,0.,1.,A,lda,ioffA);
  document.getElementById("debug_textarea").value +=
    "dlaset('L',...): A = "
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

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffA + i + j * lda ] = Number.POSITIVE_INFINITY;
    }
  }
  LaPack0.dlaset('F',m,n,0.,1.,A,lda,ioffA);
  document.getElementById("debug_textarea").value +=
    "dlaset('F',...): A = "
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

  for ( j = 0; j < m; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number.POSITIVE_INFINITY;
    }
  }
  LaPack0.dlaset('U',n,m,0.,1.,B,ldb,ioffB);
  document.getElementById("debug_textarea").value +=
    "dlaset('U',...): B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + " "
    + B[ ioffB + 2 + 2 * ldb ] + " "
    + B[ ioffB + 0 + 3 * ldb ] + " "
    + B[ ioffB + 1 + 3 * ldb ] + " "
    + B[ ioffB + 2 + 3 * ldb ]  + "\n";

  for ( j = 0; j < m; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number.POSITIVE_INFINITY;
    }
  }
  LaPack0.dlaset('L',n,m,0.,1.,B,ldb,ioffB);
  document.getElementById("debug_textarea").value +=
    "dlaset('L',...): B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + " "
    + B[ ioffB + 2 + 2 * ldb ] + " "
    + B[ ioffB + 0 + 3 * ldb ] + " "
    + B[ ioffB + 1 + 3 * ldb ] + " "
    + B[ ioffB + 2 + 3 * ldb ]  + "\n";

  for ( j = 0; j < m; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number.POSITIVE_INFINITY;
    }
  }
  LaPack0.dlaset('F',n,m,0.,1.,B,ldb,ioffB);
  document.getElementById("debug_textarea").value +=
    "dlaset('F',...): B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + " "
    + B[ ioffB + 2 + 2 * ldb ] + " "
    + B[ ioffB + 0 + 3 * ldb ] + " "
    + B[ ioffB + 1 + 3 * ldb ] + " "
    + B[ ioffB + 2 + 3 * ldb ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlasr() {
  document.getElementById("debug_textarea").value +=
    "testing dlasr *************" + "\n";
  var m = 4;
  var n = 3;
  var lda = 5;
  var ioffA = 1;
  var ioffc = 2;
  var ioffs = 3;
  var A = new Array( ioffA + n * lda );
  var c = new Array( ioffc + m-1 );
  var s = new Array( ioffs + m-1 );
  for ( var i = 0; i < m - 1; i ++ ) {
    c[ ioffc + i ] = 2 * i - 1;
    s[ ioffs + i ] = 2 - 3 * i;
  }
  for ( var j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) A[ ioffA + i + j * lda ] = i + 2 * j;
  }
  LaPack0.dlasr('L','V','F',m,n,c,s,A,lda,ioffc,ioffs,ioffA);
  document.getElementById("debug_textarea").value +=
    "dlasr('L','V','F',...): A = "
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

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) A[ ioffA + i + j * lda ] = i + 2 * j;
  }
  LaPack0.dlasr('L','V','B',m,n,c,s,A,lda,ioffc,ioffs,ioffA);
  document.getElementById("debug_textarea").value +=
    "dlasr('L','V','B',...): A = "
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

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) A[ ioffA + i + j * lda ] = i + 2 * j;
  }
  LaPack0.dlasr('L','T','F',m,n,c,s,A,lda,ioffc,ioffs,ioffA);
  document.getElementById("debug_textarea").value +=
    "dlasr('L','T','F',...): A = "
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

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) A[ ioffA + i + j * lda ] = i + 2 * j;
  }
  LaPack0.dlasr('L','T','B',m,n,c,s,A,lda,ioffc,ioffs,ioffA);
  document.getElementById("debug_textarea").value +=
    "dlasr('L','T','B',...): A = "
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

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) A[ ioffA + i + j * lda ] = i + 2 * j;
  }
  LaPack0.dlasr('L','B','F',m,n,c,s,A,lda,ioffc,ioffs,ioffA);
  document.getElementById("debug_textarea").value +=
    "dlasr('L','B','F',...): A = "
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

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) A[ ioffA + i + j * lda ] = i + 2 * j;
  }
  LaPack0.dlasr('L','B','B',m,n,c,s,A,lda,ioffc,ioffs,ioffA);
  document.getElementById("debug_textarea").value +=
    "dlasr('L','B','B',...): A = "
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

  c = new Array( ioffc + n - 1 );
  s = new Array( ioffs + n - 1 );
  for ( i = 0; i < n - 1; i ++ ) {
    c[ ioffc + i ] = 2 * i - 1;
    s[ ioffs + i ] = 2 - 3 * i;
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) A[ ioffA + i + j * lda ] = i + 2 * j;
  }
  LaPack0.dlasr('R','V','F',m,n,c,s,A,lda,ioffc,ioffs,ioffA);
  document.getElementById("debug_textarea").value +=
    "dlasr('R','V','F',...): A = "
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

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) A[ ioffA + i + j * lda ] = i + 2 * j;
  }
  LaPack0.dlasr('R','V','B',m,n,c,s,A,lda,ioffc,ioffs,ioffA);
  document.getElementById("debug_textarea").value +=
    "dlasr('R','V','B',...): A = "
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

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) A[ ioffA + i + j * lda ] = i + 2 * j;
  }
  LaPack0.dlasr('R','T','F',m,n,c,s,A,lda,ioffc,ioffs,ioffA);
  document.getElementById("debug_textarea").value +=
    "dlasr('R','T','F',...): A = "
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

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) A[ ioffA + i + j * lda ] = i + 2 * j;
  }
  LaPack0.dlasr('R','T','B',m,n,c,s,A,lda,ioffc,ioffs,ioffA);
  document.getElementById("debug_textarea").value +=
    "dlasr('R','T','B',...): A = "
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

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) A[ ioffA + i + j * lda ] = i + 2 * j;
  }
  LaPack0.dlasr('R','B','F',m,n,c,s,A,lda,ioffc,ioffs,ioffA);
  document.getElementById("debug_textarea").value +=
    "dlasr('R','B','F',...): A = "
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

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) A[ ioffA + i + j * lda ] = i + 2 * j;
  }
  LaPack0.dlasr('R','B','B',m,n,c,s,A,lda,ioffc,ioffs,ioffA);
  document.getElementById("debug_textarea").value +=
    "dlasr('R','B','B',...): A = "
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
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlasrt() {
  document.getElementById("debug_textarea").value +=
    "testing dlasrt *************" + "\n";

  var n = 19;
  var ioffd = 1;
  var d = new Array( ioffd + n );
  var info = new IntReference();
  for ( var i = 0; i < n; i ++ ) d[ ioffd + i ] = Math.random();
  LaPack0.dlasrt('I',n,d,info,ioffd);
  document.getElementById("debug_textarea").value +=
    "dlasrt('I'): info = " + info.getValue() + "\n";
  for ( i = 0; i < n; i ++ ) {
    document.getElementById("debug_textarea").value +=
      "d[" + i + "] = " + d[ ioffd + i ] + "\n";
  }
  for ( i = 0; i < n; i ++ ) d[ ioffd + i ] = Math.random();
  LaPack0.dlasrt('D',n,d,info,ioffd);
  document.getElementById("debug_textarea").value +=
    "dlasrt('D'): info = " + info.getValue() + "\n";
  for ( i = 0; i < n; i ++ ) {
    document.getElementById("debug_textarea").value +=
      "d[" + i + "] = " + d[ ioffd + i ] + "\n";
  }

  n = 91;
  d = new Array( ioffd + n );
  for ( i = 0; i < n; i ++ ) d[ ioffd + i ] = Math.random();
  LaPack0.dlasrt('I',n,d,info,ioffd);
  document.getElementById("debug_textarea").value +=
    "dlasrt('I'): info = " + info.getValue() + "\n";
  for ( i = 0; i < n; i ++ ) {
    document.getElementById("debug_textarea").value +=
      "d[" + i + "] = " + d[ ioffd + i ] + "\n";
  }
  for ( i = 0; i < n; i ++ ) d[ ioffd + i ] = Math.random();
  LaPack0.dlasrt('D',n,d,info,ioffd);
  document.getElementById("debug_textarea").value +=
    "dlasrt('D'): info = " + info.getValue() + "\n";
  for ( i = 0; i < n; i ++ ) {
    document.getElementById("debug_textarea").value +=
      "d[" + i + "] = " + d[ ioffd + i ] + "\n";
  }
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlassq() {
  document.getElementById("debug_textarea").value +=
    "testing dlassq *************" + "\n";
  var n = 5;
  var scalReference = new NumberReference(2.)
  var sumsqReference = new NumberReference(1.);
  var incx = 1;
  var ioffx = 1;
  var x = new Array( ioffx + n * incx );
  for ( var i = 0; i < n; i ++ ) x[ ioffx + i * incx ] = i;
  LaPack0.dlassq(n,x,incx,scalReference,sumsqReference,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlassq: scal,sumsq = " + scalReference + " " + sumsqReference + "\n";

  scalReference.setValue( 2. );
  sumsqReference.setValue( 1. );
  incx = 2;
  x = new Array( ioffx + n * incx );
  for ( i = 0; i < n; i ++ ) x[ ioffx + i * incx ] = i;
  LaPack0.dlassq(n,x,incx,scalReference,sumsqReference,ioffx);
  document.getElementById("debug_textarea").value +=
    "dlassq: scal,sumsq = " + scalReference + " " + sumsqReference + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlaswp() {
  document.getElementById("debug_textarea").value +=
    "testing dlaswp *************" + "\n";
  var lda = 5;
  var m = 4;
  var n = 3;
  var k1 = 1;
  var k2 = m;
  var ioffA = 1;
  var ioffipiv = 2;
  var A = new Array( ioffA + lda * n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffA + i + j * lda ] = i + 2 * j;
    }
  }
  var incx = 1;
  var ipiv = new Array( ioffipiv + 12 );
  for ( i=0; i < 12; i ++ ) ipiv[ ioffipiv + i ] = -1;
  ipiv[ ioffipiv + 0 ] = 1;
  ipiv[ ioffipiv + 1 ] = 3;
  ipiv[ ioffipiv + 2 ] = 4;
  ipiv[ ioffipiv + 3 ] = 4;
  LaPack0.dlaswp(n,A,lda,k1,k2,ipiv,incx,ioffA,ioffipiv);
  document.getElementById("debug_textarea").value +=
    "dlaswp: A = "
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

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffA + i + j * lda ] = i + 2 * j;
    }
  }
  k1=2;
  k2=3;
  incx=2;
  for ( i=0; i < 12; i ++ ) ipiv[ ioffipiv + i ] = -1;
  ipiv[ ioffipiv + k1 - 1 ] = k1;
  ipiv[ ioffipiv + k1 + incx - 1 ] = m;
  LaPack0.dlaswp(n,A,lda,k1,k2,ipiv,incx,ioffA,ioffipiv);
  document.getElementById("debug_textarea").value +=
    "dlaswp: A = "
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

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffA + i + j * lda ] = i + 2 * j;
    }
  }
  k1=2;
  k2=3;
  incx=-2;
  for ( i=0; i < 12; i ++ ) ipiv[ ioffipiv + i ] = -1;
  ipiv[ ioffipiv - ( k2 - 1 ) * incx ] = k1;
  ipiv[ ioffipiv - incx ] = k1;
  LaPack0.dlaswp(n,A,lda,k1,k2,ipiv,incx,ioffA,ioffipiv);
  document.getElementById("debug_textarea").value +=
    "dlaswp: A = "
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
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlauu2() {
  document.getElementById("debug_textarea").value +=
    "testing dlauu2 *************" + "\n";
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
  LaPack0.dlauu2('U',n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dlauu2('U',...): A = "
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
  LaPack0.dlauu2('L',n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dlauu2('L',...): A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dpoequ() {
  document.getElementById("debug_textarea").value +=
    "testing dpoequ *************" + "\n";
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
  LaPack0.dpoequ(n,A,lda,s,scondReference,amaxReference,info,ioffa,ioffs);
  document.getElementById("debug_textarea").value +=
    "dpoequ: scond,amax = " + scondReference + " " + amaxReference + "\n";
  document.getElementById("debug_textarea").value +=
    "s = " + s[ ioffs + 0 ] + " " + s[ ioffs + 1 ] + " "
    + s[ ioffs + 2 ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dpotf2() {
  document.getElementById("debug_textarea").value +=
    "testing dpotf2 *************" + "\n";
  var n = 3;
  var lda = 4;
  var ioffa = 1;
  var A = new Array( ioffa + lda * n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i <= j; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + j );
    }
  }
  var info = new IntReference();
  LaPack0.dpotf2('U',n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dpotf2('U',...): A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ]  + "\n";
  A = new Array( ioffa + lda * n );
  for ( j = 0; j < n; j ++ ) {
    for ( i = j; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + j );
    }
  }
  LaPack0.dpotf2('L',n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dpotf2('L',...): A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dpotrs() {
  document.getElementById("debug_textarea").value +=
    "testing dpotrs *************" + "\n";
  var n = 3;
  var nrhs = 1;
  var lda = 4;
  var ldb = 5;
  var ioffa = 1;
  var ioffb = 2;
  var A = new Array( ioffa + lda * n );
  var B = new Array( ioffb + ldb * nrhs );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i <= j; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + j );
    }
  }
  for ( i=0; i < n; i ++ ) B[ ioffb + i ] = i + 1;
  var info = new IntReference();
  LaPack0.dpotrs('U',n,nrhs,A,lda,B,ldb,info,ioffa,ioffb);
  document.getElementById("debug_textarea").value +=
    "dpotrs('U',...): B = " + B[ ioffb + 0 ] + " " + B[ ioffb + 1 ]
    + " " + B[ ioffb + 2 ] + "\n";

  A = new Array( ioffa + lda * n );
  for ( j = 0; j < n; j ++ ) {
    for ( i = j; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + j );
    }
  }
  for ( i=0; i < n; i ++ ) B[ ioffb + i ] = i + 1;
  LaPack0.dpotrs('L',n,nrhs,A,lda,B,ldb,info,ioffa,ioffb);
  document.getElementById("debug_textarea").value +=
    "dpotrs('L',...): B = " + B[ ioffb + 0 ] + " " + B[ ioffb + 1 ]
    + " " + B[ ioffb + 2 ] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dptcon() {
  document.getElementById("debug_textarea").value +=
    "testing dptcon *************" + "\n";
  var n = 3;
  var ioffd = 1;
  var ioffe = 2;
  var ioffwork = 3;
  var d = new Array( ioffd + n );
  var e = new Array( ioffe + n - 1 );
  var work = new Array( ioffwork + n );
  var anorm = 6.;
  for ( var i = 0; i < n; i ++ ) d[ ioffd + i ] = i + 1;
  for ( i = 0; i < n - 1; i ++ ) e[ ioffe + i ] = 2 * i + 1;
  var rcondReference = new NumberReference();
  var info = new IntReference();
  LaPack0.dptcon(n,d,e,anorm,rcondReference,work,info,ioffd,ioffe,
    ioffwork);
  document.getElementById("debug_textarea").value +=
    "dptcon: rcond = " + rcondReference.getValue() + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dpttrf() {
  document.getElementById("debug_textarea").value +=
    "testing dpttrf *************" + "\n";
  var ioffd = 1;
  var ioffe = 2;
  var d = new Array( ioffd + 11 );
  var e = new Array( ioffe + 10 );
  for ( var i = 0; i < 11; i ++ ) d[ ioffd + i ] = 2 * i + 3;
  for ( i = 0; i < 10; i ++ ) e[ ioffe + i ] = i;
  var n = 3;
  var info = new IntReference();
  LaPack0.dpttrf(n,d,e,info,ioffd,ioffe);
  document.getElementById("debug_textarea").value +=
    "dpttrf(3,...): d = "
    + d[ ioffd + 0 ] + " "
    + d[ ioffd + 1 ] + " "
    + d[ ioffd + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "e = "
    + e[ ioffe + 0 ] + " "
    + e[ ioffe + 1 ]  + "\n";

  for ( i = 0; i < 11; i ++ ) d[ ioffd + i ] = 2 * i + 3;
  for ( i = 0; i < 10; i ++ ) e[ ioffe + i ] = i;
  n = 7;
  LaPack0.dpttrf(n,d,e,info,ioffd,ioffe);
  document.getElementById("debug_textarea").value +=
    "dpttrf(7,...): d = "
    + d[ ioffd + 0 ] + " "
    + d[ ioffd + 1 ] + " "
    + d[ ioffd + 2 ] + " "
    + d[ ioffd + 3 ] + " "
    + d[ ioffd + 4 ] + " "
    + d[ ioffd + 5 ] + " "
    + d[ ioffd + 6 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "e = "
    + e[ ioffe + 0 ] + " "
    + e[ ioffe + 1 ] + " "
    + e[ ioffe + 2 ] + " "
    + e[ ioffe + 3 ] + " "
    + e[ ioffe + 4 ] + " "
    + e[ ioffe + 5 ]  + "\n";

  for ( i = 0; i < 11; i ++ ) d[ ioffd + i ] = 2 * i + 3;
  for ( i = 0; i < 10; i ++ ) e[ ioffe + i ] = i;
  n = 11;
  LaPack0.dpttrf(n,d,e,info,ioffd,ioffe);
  document.getElementById("debug_textarea").value +=
    "dpttrf(11,...): d = "
    + d[ ioffd + 0 ] + " "
    + d[ ioffd + 1 ] + " "
    + d[ ioffd + 2 ] + " "
    + d[ ioffd + 3 ] + " "
    + d[ ioffd + 4 ] + " "
    + d[ ioffd + 5 ] + " "
    + d[ ioffd + 6 ] + " "
    + d[ ioffd + 7 ] + " "
    + d[ ioffd + 8 ] + " "
    + d[ ioffd + 9 ] + " "
    + d[ ioffd + 10 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "e = "
    + e[ ioffe + 0 ] + " "
    + e[ ioffe + 1 ] + " "
    + e[ ioffe + 2 ] + " "
    + e[ ioffe + 3 ] + " "
    + e[ ioffe + 4 ] + " "
    + e[ ioffe + 5 ] + " "
    + e[ ioffe + 6 ] + " "
    + e[ ioffe + 7 ] + " "
    + e[ ioffe + 8 ] + " "
    + e[ ioffe + 9 ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dptts2() {
  document.getElementById("debug_textarea").value +=
    "testing dptts2 *************" + "\n";
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
  LaPack0.dptts2(n,nrhs,d,e,B,ldb,ioffd,ioffe,ioffb);
  document.getElementById("debug_textarea").value +=
    "B = "
    + B[ ioffb + 0 ] + " "
    + B[ ioffb + 1 ] + " "
    + B[ ioffb + 2 ] + " "
    + B[ ioffb + 3 ] + " "
    + B[ ioffb + 4 ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dsyconv() {
  document.getElementById("debug_textarea").value +=
    "testing dsyconv *************" + "\n";
  var n = 4;
  var lda = 5;
  var ioffa = 1;
  var ioffipiv = 2;
  var ioffwork = 3;
  var A = new Array( ioffa + lda * n );
  var ipiv = new Array( ioffipiv + n );
  var work = new Array( ioffwork + n );
  ipiv[ ioffipiv + 0 ] = -2;
  ipiv[ ioffipiv + 1 ] = -2;
  ipiv[ ioffipiv + 2 ] = 4;
  ipiv[ ioffipiv + 3 ] = 4;
  var info = new IntReference();
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  LaPack0.dsyconv('U','C',n,A,lda,ipiv,work,info,
    ioffa,ioffipiv,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dsyconv('U','C',...): A = "
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
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  LaPack0.dsyconv('U','R',n,A,lda,ipiv,work,info,
    ioffa,ioffipiv,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dsyconv('U','R',...): A = "
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
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  LaPack0.dsyconv('L','C',n,A,lda,ipiv,work,info,
    ioffa,ioffipiv,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dsyconv('L','C',...): A = "
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
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  LaPack0.dsyconv('L','R',n,A,lda,ipiv,work,info,
    ioffa,ioffipiv,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dsyconv('L','R',...): A = "
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
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dsyswapr() {
  document.getElementById("debug_textarea").value +=
    "testing dsyswapr *************" + "\n";
  var n = 5;
  var lda = 6
  var ioffa = 1;
  var A = new Array( ioffa + lda * n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  LaPack0.dsyswapr('U',n,A,lda,2,4,ioffa);
  document.getElementById("debug_textarea").value +=
    "dsyswapr('U',...): A = "
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
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  LaPack0.dsyswapr('L',n,A,lda,2,4,ioffa);
  document.getElementById("debug_textarea").value +=
    "dsyswapr('L',...): A = "
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
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dsytf2() {
  document.getElementById("debug_textarea").value +=
    "testing dsytf2 *************" + "\n";
  var n = 3;
  var lda = 5;
  var ioffa = 1;
  var ioffipiv = 2;
  var A = new Array( ioffa + lda * n );
  var ipiv = new Array( ioffipiv + n );
  A[ ioffa + 0 + 0 * lda ] = 4;
  A[ ioffa + 1 + 0 * lda ] = 10;
  A[ ioffa + 2 + 0 * lda ] = -2;
  A[ ioffa + 0 + 1 * lda ] = 10;
  A[ ioffa + 1 + 1 * lda ] = 2;
  A[ ioffa + 2 + 1 * lda ] = 10;
  A[ ioffa + 0 + 2 * lda ] = -2;
  A[ ioffa + 1 + 2 * lda ] = 10;
  A[ ioffa + 2 + 2 * lda ] = 4;
  var info = new IntReference();
  LaPack0.dsytf2('U',n,A,lda,ipiv,info,ioffa,ioffipiv);
  document.getElementById("debug_textarea").value +=
    "dsytf2('U',...): info = " + info.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "ipiv = "
    + ipiv[ ioffipiv + 0 ] + " "
    + ipiv[ ioffipiv + 1 ] + " "
    + ipiv[ ioffipiv + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ]  + "\n";

  A[ ioffa + 0 + 0 * lda ] = 4;
  A[ ioffa + 1 + 0 * lda ] = 10;
  A[ ioffa + 2 + 0 * lda ] = -2;
  A[ ioffa + 0 + 1 * lda ] = 10;
  A[ ioffa + 1 + 1 * lda ] = 2;
  A[ ioffa + 2 + 1 * lda ] = 10;
  A[ ioffa + 0 + 2 * lda ] = -2;
  A[ ioffa + 1 + 2 * lda ] = 10;
  A[ ioffa + 2 + 2 * lda ] = 4;
  LaPack0.dsytf2('L',n,A,lda,ipiv,info,ioffa,ioffipiv);
  document.getElementById("debug_textarea").value +=
    "dsytf2('L',...): info = " + info.getValue() + "\n";
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
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dsytri() {
  document.getElementById("debug_textarea").value +=
    "testing dsytri *************" + "\n";
  var n = 3;
  var lda = 5;
  var lwork = n * n;
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
//    LaPack2.dsytrf('U',n,A,lda,ipiv,work,lwork,info,
//      ioffa,ioffipiv,ioffwork);
//    LaPack0.dsytri('U',n,A,lda,ipiv,work,info,ioffa,ioffipiv,ioffwork);
//    document.getElementById("debug_textarea").value +=
//      "dsytrf('U',...): A = "
//      + A[ ioffa + 0 + 0 * lda ] + " "
//      + A[ ioffa + 1 + 0 * lda ] + " "
//      + A[ ioffa + 2 + 0 * lda ] + " "
//      + A[ ioffa + 0 + 1 * lda ] + " "
//      + A[ ioffa + 1 + 1 * lda ] + " "
//      + A[ ioffa + 2 + 1 * lda ] + " "
//      + A[ ioffa + 0 + 2 * lda ] + " "
//      + A[ ioffa + 1 + 2 * lda ] + " "
//      + A[ ioffa + 2 + 2 * lda ]  + "\n";

  A[ ioffa + 0 + 0 * lda ] = 4.;
  A[ ioffa + 1 + 0 * lda ] = 10.;
  A[ ioffa + 2 + 0 * lda ] = -2.;
  A[ ioffa + 0 + 1 * lda ] = 10.;
  A[ ioffa + 1 + 1 * lda ] = 2.;
  A[ ioffa + 2 + 1 * lda ] = 10.;
  A[ ioffa + 0 + 2 * lda ] = -2.;
  A[ ioffa + 1 + 2 * lda ] = 10.;
  A[ ioffa + 2 + 2 * lda ] = 4.;
//    LaPack2.dsytrf('L',n,A,lda,ipiv,work,lwork,info,
//      ioffa,ioffipiv,ioffwork);
//    LaPack0.dsytri('L',n,A,lda,ipiv,work,info,ioffa,ioffipiv,ioffwork);
//    document.getElementById("debug_textarea").value +=
//      "dsytrf('L',...): A = "
//      + A[ ioffa + 0 + 0 * lda ] + " "
//      + A[ ioffa + 1 + 0 * lda ] + " "
//      + A[ ioffa + 2 + 0 * lda ] + " "
//      + A[ ioffa + 0 + 1 * lda ] + " "
//      + A[ ioffa + 1 + 1 * lda ] + " "
//      + A[ ioffa + 2 + 1 * lda ] + " "
//      + A[ ioffa + 0 + 2 * lda ] + " "
//      + A[ ioffa + 1 + 2 * lda ] + " "
//      + A[ ioffa + 2 + 2 * lda ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dsytrs() {
  document.getElementById("debug_textarea").value +=
    "testing dsytrs *************" + "\n";
  var n = 4;
  var nrhs = 1;
  var lda = 5;
  var ldb = 6;
  var ioffa = 1;
  var ioffb = 2;
  var ioffipiv = 3;
  var A = new Array( ioffa + lda * n );
  var B = new Array( ioffb + ldb * nrhs );
  var ipiv = new Array( ioffipiv + n );
  A[ ioffa + 0 + 0 * lda ] = 5.;
  A[ ioffa + 1 + 0 * lda ] = 7.;
  A[ ioffa + 2 + 0 * lda ] = 7.;
  A[ ioffa + 3 + 0 * lda ] = -7.;
  A[ ioffa + 0 + 1 * lda ] = 7.;
  A[ ioffa + 1 + 1 * lda ] = -3.;
  A[ ioffa + 2 + 1 * lda ] = -7.;
  A[ ioffa + 3 + 1 * lda ] = -1.;
  A[ ioffa + 0 + 2 * lda ] = 7.;
  A[ ioffa + 1 + 2 * lda ] = -7.;
  A[ ioffa + 2 + 2 * lda ] = 5.;
  A[ ioffa + 3 + 2 * lda ] = 7.;
  A[ ioffa + 0 + 3 * lda ] = -7.;
  A[ ioffa + 1 + 3 * lda ] = -1.;
  A[ ioffa + 2 + 3 * lda ] = 7.;
  A[ ioffa + 3 + 3 * lda ] = -3.;
  B[ ioffb + 0 + 0 * ldb ] = 12.;
  B[ ioffb + 1 + 0 * ldb ] = -24.;
  B[ ioffb + 2 + 0 * ldb ] = 36.;
  B[ ioffb + 3 + 0 * ldb ] = 0.;
  var info = new IntReference();
  LaPack0.dsytf2('U',n,A,lda,ipiv,info,ioffa,ioffipiv);
//    document.getElementById("debug_textarea").value +=
//      "info = " + info.getValue() + "\n";
//    document.getElementById("debug_textarea").value +=
//      "A = "
//      + A[ ioffa + 0 + 0 * lda ] + " "
//      + A[ ioffa + 0 + 1 * lda ] + " "
//      + A[ ioffa + 1 + 1 * lda ] + " "
//      + A[ ioffa + 0 + 2 * lda ] + " "
//      + A[ ioffa + 1 + 2 * lda ] + " "
//      + A[ ioffa + 2 + 2 * lda ] + " "
//      + A[ ioffa + 0 + 3 * lda ] + " "
//      + A[ ioffa + 1 + 3 * lda ] + " "
//      + A[ ioffa + 2 + 3 * lda ] + " "
//      + A[ ioffa + 3 + 3 * lda ]  + "\n";
//    document.getElementById("debug_textarea").value +=
//      "ipiv = "
//      + ipiv[ ioffipiv + 0 ] + " "
//      + ipiv[ ioffipiv + 1 ] + " "
//      + ipiv[ ioffipiv + 2 ] + " "
//      + ipiv[ ioffipiv + 3 ]  + "\n";
  LaPack0.dsytrs('U',n,nrhs,A,lda,ipiv,B,ldb,info,
    ioffa,ioffipiv,ioffb);
  document.getElementById("debug_textarea").value +=
    "dsytrs('U',...): B = "
    + B[ ioffb + 0 + 0 * ldb ] + " "
    + B[ ioffb + 1 + 0 * ldb ] + " "
    + B[ ioffb + 2 + 0 * ldb ] + " "
    + B[ ioffb + 3 + 0 * ldb ]  + "\n";

  A[ ioffa + 0 + 0 * lda ] = 5.;
  A[ ioffa + 1 + 0 * lda ] = 7.;
  A[ ioffa + 2 + 0 * lda ] = 7.;
  A[ ioffa + 3 + 0 * lda ] = -7.;
  A[ ioffa + 0 + 1 * lda ] = 7.;
  A[ ioffa + 1 + 1 * lda ] = -3.;
  A[ ioffa + 2 + 1 * lda ] = -7.;
  A[ ioffa + 3 + 1 * lda ] = -1.;
  A[ ioffa + 0 + 2 * lda ] = 7.;
  A[ ioffa + 1 + 2 * lda ] = -7.;
  A[ ioffa + 2 + 2 * lda ] = 5.;
  A[ ioffa + 3 + 2 * lda ] = 7.;
  A[ ioffa + 0 + 3 * lda ] = -7.;
  A[ ioffa + 1 + 3 * lda ] = -1.;
  A[ ioffa + 2 + 3 * lda ] = 7.;
  A[ ioffa + 3 + 3 * lda ] = -3.;
  B[ ioffb + 0 + 0 * ldb ] = 12.;
  B[ ioffb + 1 + 0 * ldb ] = -24.;
  B[ ioffb + 2 + 0 * ldb ] = 36.;
  B[ ioffb + 3 + 0 * ldb ] = 0.;
  LaPack0.dsytf2('L',n,A,lda,ipiv,info,ioffa,ioffipiv);
  document.getElementById("debug_textarea").value +=
    "info = " + info.getValue() + "\n";
  LaPack0.dsytrs('L',n,nrhs,A,lda,ipiv,B,ldb,info,
    ioffa,ioffipiv,ioffb);
  document.getElementById("debug_textarea").value +=
    "dsytrs('L',...): B = "
    + B[ ioffb + 0 + 0 * ldb ] + " "
    + B[ ioffb + 1 + 0 * ldb ] + " "
    + B[ ioffb + 2 + 0 * ldb ] + " "
    + B[ ioffb + 3 + 0 * ldb ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dtrti2() {
  document.getElementById("debug_textarea").value +=
    "testing dtrti2 *************" + "\n";
  var n = 3;
  var lda = 4;
  var ioffa = 1;
  var A = new Array( ioffa + lda * n );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i <= j; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  var info = new IntReference();
  LaPack0.dtrti2('U','N',n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dtrti2('U','N',...): A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ]  + "\n";

  A = new Array( ioffa + lda * n );
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < j; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  LaPack0.dtrti2('U','U',n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dtrti2('U','U',...): A = "
    + A[ ioffa + 0 + 1 * lda ] + " "
    + A[ ioffa + 0 + 2 * lda ] + " "
    + A[ ioffa + 1 + 2 * lda ]  + "\n";

  A = new Array( ioffa + lda * n );
  for ( j = 0; j < n; j ++ ) {
    for ( i = j; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  LaPack0.dtrti2('L','N',n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dtrti2('L','N',...): A = "
    + A[ ioffa + 0 + 0 * lda ] + " "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 1 + 1 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ] + " "
    + A[ ioffa + 2 + 2 * lda ] + "\n";

  A = new Array( ioffa + lda * n );
  for ( j = 0; j < n; j ++ ) {
    for ( i = j + 1; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1 + i + 2 * j;
    }
  }
  LaPack0.dtrti2('L','U',n,A,lda,info,ioffa);
  document.getElementById("debug_textarea").value +=
    "dtrti2('L','U',...): A = "
    + A[ ioffa + 1 + 0 * lda ] + " "
    + A[ ioffa + 2 + 0 * lda ] + " "
    + A[ ioffa + 2 + 1 * lda ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dtrtrs() {
  document.getElementById("debug_textarea").value +=
    "testing dtrtrs *************" + "\n";
  var n = 3;
  var nrhs = 2;
  var ioffA = 1;
  var ioffB = 2;
  var lda = 4;
  var ldb = 4;
  var A = new Array( ioffA + lda * n );
  var B = new Array( ioffB + ldb * nrhs );
  var info = new IntReference();
  var i = -1;
  var j = -1;
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = Number( 1 + i + j * 2 );
    }
  }
  for ( j = 0; j < nrhs; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  LaPack0.dtrtrs('L','N','U',n,nrhs,A,lda,B,ldb,info,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrtrs('L','N','U',n,nrhs,A,lda,B,ldb,info) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < nrhs; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  LaPack0.dtrtrs('L','N','N',n,nrhs,A,lda,B,ldb,info,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrtrs('L','N','N',n,nrhs,A,lda,B,ldb,info) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < nrhs; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  LaPack0.dtrtrs('L','T','U',n,nrhs,A,lda,B,ldb,info,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrtrs('L','T','U',n,nrhs,A,lda,B,ldb,info) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < nrhs; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  LaPack0.dtrtrs('L','T','N',n,nrhs,A,lda,B,ldb,info,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrtrs('L','T','N',n,nrhs,A,lda,B,ldb,info) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";

  for ( j = 0; j < nrhs; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  LaPack0.dtrtrs('U','N','U',n,nrhs,A,lda,B,ldb,info,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrtrs('U','N','U',n,nrhs,A,lda,B,ldb,info) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < nrhs; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  LaPack0.dtrtrs('U','N','N',n,nrhs,A,lda,B,ldb,info,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrtrs('U','N','N',n,nrhs,A,lda,B,ldb,info) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < nrhs; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  LaPack0.dtrtrs('U','T','U',n,nrhs,A,lda,B,ldb,info,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrtrs('U','T','U',n,nrhs,A,lda,B,ldb,info) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < nrhs; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  LaPack0.dtrtrs('U','T','N',n,nrhs,A,lda,B,ldb,info,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrtrs('U','T','N',n,nrhs,A,lda,B,ldb,info) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_ieeeck() {
  document.getElementById("debug_textarea").value +=
    "testing ieeeck *************" + "\n";
  var val = LaPack0.ieeeck(0,0.,1.);
  document.getElementById("debug_textarea").value +=
    "ieeeck(0,...) = " + val + "\n";
  val = LaPack0.ieeeck(1,0.,1.);
  document.getElementById("debug_textarea").value +=
    "ieeeck(1,...) = " + val + "\n";
}
