function LaPack5() {
}
//*************************************************************************
LaPack5.dbdsqr = function( uplo, n, ncvt, nru, ncc, d, e, Vt,
ldvt, U, ldu, C, ldc, work, info, ioffd, ioffe, ioffvt, ioffu, ioffc,
ioffwork ) {
  throw new Error("not tested");
  var maxitr = 6;
  info.setValue( 0 );
  var lower = ( uplo.charAt(0).toUpperCase() == 'L' );
  if ( uplo.charAt(0).toUpperCase() != 'U' && ! lower ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( ncvt < 0 ) info.setValue( -3 );
  else if ( nru < 0 ) info.setValue( -4 );
  else if ( ncc < 0 ) info.setValue( -5 );
  else if ( ( ncvt == 0 && ldvt < 1 ) ||
  ( ncvt > 0 && ldvt < Math.max( 1, n ) ) ) {
    info.setValue( -9 );
  } else if ( ldu < Math.max( 1, nru ) ) info.setValue( -11 );
  else if ( ( ncc == 0 && ldc < 1 ) ||
  ( ncc > 0 && ldc < Math.max( 1, n ) ) ) {
    info.setValue( -13 );
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dbdsqr', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  var csReference = new NumberReference();
  var snReference = new NumberReference();
  var rReference = new NumberReference();
  var goto160 =( n == 1 );
  if ( ! goto160 ) {
    var rotate = ( ncvt > 0 || nru > 0 || ncc > 0 );
    if ( ! rotate ) {
      LaPack4.dlasq1( n, d, e, work, info, ioffd, ioffe, ioffwork );
      if ( info.getValue() != 2 ) return;
      info.setValue( 0 );
    }
    var nm1 = n - 1;
    var nm12 = nm1 + nm1;
    var nm13 = nm12 + nm1;
    var idr = 0;
    var eps = LaPack0.dlamch( 'Epsilon' );
    var unfl = LaPack0.dlamch( 'Safe minimum' );
    if ( lower ) {
      for ( var i = 1; i <= n - 1; i ++ ) {
        LaPack1.dlartg( d[ ioffd + i - 1], e[ ioffe + i - 1 ],
          cs, sn, r );
        d[ ioffd + i - 1 ] = r.getValue();
        e[ ioffe + i - 1 ] = sn.getValue() * d[ ioffd + i ];
        d[ ioffd + i ] *= cs.getValue();
        work[ ioffwork + i - 1 ] = cs.getValue();
        work[ ioffwork + nm1 + i - 1 ] = sn.getValue();
      } // 10
      if ( nru > 0 ) {
        LaPack0.dlasr( 'R', 'V', 'F', nru, n, work, work, U, ldu,
          ioffwork, ioffwork + n - 1, ioffu );
      }
      if ( ncc > 0 ) {
        LaPack0.dlasr( 'L', 'V', 'F', n, ncc, work, work, C, ldc,
          ioffwork, ioffwork + n - 1, ioffc );
      }
    }
    var tolmul =
      Math.max( 10., Math.min( 100., Math.pow( eps, -0.125 ) ) );
    var tol = tolmul * eps;
    var smax = 0.;
    for ( i = 1; i <= n; i ++ ) {
      smax = Math.max( smax, Math.abs( d[ ioffd + i - 1 ] ) );
    }
    for ( i = 1; i <= n - 1; i ++ ) {
      smax = Math.max( smax, Math.abs( e[ ioffe + i - 1 ] ) );
    }
    var sminl = 0.;
    if ( tol >= 0. ) {
      var sminoa = Math.abs( d[ ioffd ] );
      if ( sminoa != 0. ) {
        var mu = sminoa;
        for ( i = 2; i <= n; i ++ ) {
          mu = Math.abs( d[ ioffd + i - 1 ] )
            * ( mu / ( mu + Math.abs( e[ ioffe + i - 2 ] ) ) );
          sminoa = Math.min( sminoa, mu );
          if ( sminoa == 0. ) break;
        } // 40
      } // 50
      sminoa /= Math.sqrt( Number( n ) );
      var thresh =
        Math.max( tol * sminoa, maxitr * n * n * unfl );
    } else {
      thresh =
        Math.max( Math.abs( tol ) * smax, maxitr * n * n * unfl );
    }
    var maxit = maxitr * n * n;
    var iter = 0;
    var oldll = -1;
    var oldm = -1;
    var m = n;
    var goto200 = false;
    var sigmnReference = new NumberReference();
    var sigmxReference = new NumberReference();
    var sinrReference = new NumberReference();
    var cosrReference = new NumberReference();
    var sinlReference = new NumberReference();
    var coslReference = new NumberReference();
    while ( true ) { // 60
      var goto60 = false;
      if ( m <= 1 ) break;
      if ( iter > maxit ) {
        goto200 = true;
        break;
      }
      if ( tol < 0. && Math.abs( d[ ioffd + m - 1 ] ) <= thresh ) {
        d[ ioffd + m - 1 ] = 0.;
      }
      var smax = Math.abs( d[ ioffd + m - 1 ] );
      var smin = smax;
      var goto80 = false;
      for ( var lll = 1; lll <= m - 1; lll ++ ) {
        var ll = m - lll;
        var abss = Math.abs( d[ ioffd + ll - 1 ] );
        var abse = Math.abs( e[ ioffe + ll - 1 ] );
        if ( tol < 0. && abss <= thresh ) d[ ioffd + ll - 1 ] = 0.;
        if ( abse <= thresh ) {
          goto80 = true;
          break;
        }
        smin = Math.min( smin, abss );
        smax = Math.max( smax, Math.max( abss, abse ) );
      } // 70
      if ( ! goto80 ) ll = 0;
      else {
        e[ ioffe + ll - 1 ] = 0.;
        if ( ll == m - 1 ) {
          m --;
          continue; // goto 60
        }
      } // 90
      ll ++;
      if ( ll == m - 1 ) {
        LaPack1.dlasv2( d[ ioffd + m - 2 ], e[ ioffe + m - 2 ],
          d[ ioffd + m - 1 ], sigmn, sigmx, sinr, cosr, sinl, cosl );
        d[ ioffd + m - 2 ] = sigmx.getValue();
        e[ ioffe + m - 2 ] = 0.;
        d[ ioffd + m - 1 ] = sigmn.getValue();
        if ( ncvt > 0 ) {
          Blas1.drot( ncvt, Vt, ldvt, Vt, ldvt, cosr.getValue(), sinr.getValue(),
            ioffvt + m - 2, ioffv2 + m - 1 );
        }
        if( nru > 0 ) {
          Blas1.drot( nru, U, 1, U, 1, cosl.getValue(), sinl.getValue(),
            ioffu + ( m - 2 ) * ldu, ioffu + ( m - 1 ) * ldu );
        }
        if ( ncc > 0 ) {
          Blas1.drot( ncc, C, ldc, C, ldc, cosl.getValue(), sinl.getValue(),
            offc + m - 2, ioffc + m - 1 );
        }
        m -= 2;
        continue; // goto 60
      }
      if ( ll > oldm || m < oldll ) {
        var idr = ( Math.abs( d[ ioffd + ll - 1 ] ) >=
        Math.abs( d[ ioffd + m - 1 ] ) ? 1 : 2 );
      }
      if ( idr == 1 ) {
        if ( Math.abs( e[ ioffe + m - 2 ] ) <=
        Math.abs( tol ) * Math.abs( d[ ioffd + m - 1 ] ) ||
        ( tol < 0. && Math.abs( e[ ioffe + m - 2 ] ) <= thresh ) ) {
          e[ ioffe + m - 2 ] = 0.;
          continue; // goto 60
        }
        if ( tol >= 0. ) {
          mu = Math.abs( d[ ioffd + ll - 1 ] );
          sminl = mu;
          for ( lll = ll; lll <= m - 1; lll ++ ) {
            if ( Math.abs( e[ ioffe + lll - 1 ] ) <= tol * mu ) {
              e[ ioffe + lll - 1 ] = 0.;
              goto60 = true;
              break;
            }
            mu = Math.abs( d[ ioffd + lll ] )
              * ( mu / ( mu + Math.abs( e[ ioffe + lll - 1 ] ) ) );
            sminl = Math.min( sminl, mu );
          } // 100
          if ( goto60 ) continue;
        }
      } else {
        if ( Math.abs( e[ ioffe + ll - 1 ] ) <=
        Math.abs( tol ) * Math.abs( d[ ioffd + ll - 1 ] ) ||
        ( tol < 0. && Math.abs( e[ ioffe + ll - 1 ] ) <= thresh ) ) {
          e[ ioffe + ll - 1 ] = 0.;
          continue; // goto 60
        }
        if ( tol >= 0. ) {
          mu = Math.abs( d[ ioffd + m - 1 ] );
          sminl = mu;
          for ( lll = m - 1; lll >= ll; lll -- ) {
            if ( Math.abs( e[ ioffe + lll - 1 ] ) <= tol * mu ) {
              e[ ioffe + lll - 1 ] = 0.;
              goto60 = true;
              break;
            }
            mu = Math.abs( d[ ioffd + lll - 1 ] )
              * ( mu / ( mu + Math.abs( e[ ioffe + lll - 1 ] ) ) );
            sminl = Math.min( sminl, mu );
          } // 110;
          if ( goto60 ) continue;
        }
      }
      oldll = ll;
      oldm = m;
      var shiftReference = new NumberReference();
      if ( tol >= 0. &&
      n * tol * ( sminl / smax ) <= Math.max( eps, 0.01 * tol ) ) {
        shift.setValue( 0. );
      } else {
        if ( idr == 1 ) {
          var sll = Math.abs( d[ ioffd + ll - 1 ] );
          LaPack0.dlas2( d[ ioffd + m - 2 ], e[ ioffe + m - 2 ],
            d[ ioffd + m - 1 ], shift, r );
        } else {
          sll = Math.abs( d[ ioffd + m - 1 ] );
          LaPack0.dlas2( d[ ioffd + ll - 1 ], ioffe[ ioffe + ll - 1 ],
            d[ ioffd + ll ], shift, r );
        }
        if ( sll > 0. ) {
          if ( Math.pow( shift.getValue() / sll, 2 ) < eps ) {
            shift.setValue( 0. );
          }
        }
      }
      iter += m - ll;
      var oldcsReference = new NumberReference( );
      var oldsnReference = new NumberReference( );
      var direfReference = new NumberReference( );
      if ( shift.getValue() == 0. ) {
        if ( idr == 1 ) {
          cs.setValue( 1. );
          oldcs.setValue( 1. );
          for ( i = ll; i <= m - 1; i ++ ) {
            LaPack1.dlartg( d[ ioffd + i - 1 ] * cs.getValue(),
              e[ ioffe + i - 1 ], cs, sn, r );
            if ( i > ll ) e[ ioffe + i - 2 ] = oldsn.getValue() * r.getValue();
            diref.setValue( d[ ioffd + i - 1 ] );
            LaPack1.dlartg( oldcs.getValue() * r.getValue(),
              d[ ioffd + i ] * sn.getValue(), oldcs, oldsn, diref );
            d[ ioffd + i - 1 ] = diref.getValue();
            work[ ioffwork + i - ll ] = cs.getValue();
            work[ ioffwork + i - ll + nm1 ] = sn.getValue();
            work[ ioffwork + i - ll + nm12 ] = oldcs.getValue();
            work[ ioffwork + i - ll + nm13 ] = oldsn.getValue();
          } // 120
          var h = d[ ioffd + m - 1 ] * cs.getValue();
          d[ ioffd + m - 1 ] = h * oldcs.getValue();
          e[ ioffe + m - 2 ] = h * oldsn.getValue();
          if ( ncvt > 0 ) {
            LaPack0.dlasr( 'L', 'V', 'F', m - ll + 1, ncvt, work, work,
              Vt, ldvt, ioffwork, ioffwork + n - 1, ioffvt + ll - 1 );
          }
          if ( nru > 0 ) {
            LaPack0.dlasr( 'R', 'V', 'F', nru, m - ll + 1, work, work,
              U, ldu, ioffwork + nm12, ioffwork + nm13,
              ioffu + ( ll - 1 ) * ldu );
          }
          if ( ncc > 0 ) {
            LaPack0.dlasr( 'L', 'V', 'F', m - ll + 1, ncc, work, work,
              C, ldc, ioffwork + nm12, ioffwork + nm13,
              ioffc + ll - 1 );
          }
          if ( Math.abs( e[ ioffe + m - 2 ] ) <= thresh ) {
            e[ ioffe + m - 2 ] = 0.;
          }
        } else {
          cs.setValue( 1. );
          oldcs.setValue( 1. );
          for ( i = m; i >= ll + 1; i -- ) {
            LaPack1.dlartg( d[ ioffd + i - 1 ] * cs.getValue(),
              e[ ioffe + i - 2 ], cs, sn, r );
            if ( i < m ) e[ ioffe + i - 1 ] = oldsn.getValue() * r.getValue();
            diref.setValue( d[ ioffd + i - 1 ] );
            LaPack1.dlartg( oldcs.getValue() * r.getValue(),
              d[ ioffd + i - 2 ] * sn.getValue(), oldcs, oldsn, diref );
            d[ ioffd + i - 1 ] = diref.getValue();
            work[ ioffwork + i - ll - 1 ] = cs.getValue();
            work[ ioffwork + i - ll + nm1 - 1 ] = - sn.getValue();
            work[ ioffwork + i - ll + nm12 - 1 ] = oldcs.getValue();
            work[ ioffwork + i - ll + nm13 - 1 ] = - oldsn.getValue();
          } // 130
          h = d[ ioffd + ll - 1 ] * cs.getValue();
          d[ ioffd + ll - 1 ] = h * oldcs.getValue();
          e[ ioffe + ll - 1 ] = h * oldsn.getValue();
          if ( ncvt > 0 ) {
            LaPack0.dlasr( 'L', 'V', 'B', m - ll + 1, ncvt, work, work,
              Vt, ldvt, ioffwork + nm12, ioffwork + nm13,
              ioffvt + ll - 1 );
          }
          if ( nru > 0 ) {
            LaPack0.dlasr( 'R', 'V', 'B', nru, m - ll + 1, work, work,
              U, ldu, ioffwork, ioffwork + n - 1,
              ioffu + ( ll - 1 ) * ldu );
          }
          if ( ncc > 0 ) {
            LaPack0.dlasr( 'L', 'V', 'B', m - ll + 1, ncc, work, work,
              C, ldc, ioffwork, ioffwork + n - 1, ioffc + ll - 1 );
          }
          if ( Math.abs( e[ ioffe + ll - 1 ] ) <= thresh ) {
            e[ ioffe + ll - 1 ] = 0.;
          }
        }
      } else {
        if ( idir == 1 ) {
          var f = ( Math.abs( d[ ioffd + ll - 1 ] ) - shift )
            * ( ( d[ ioffd + ll - 1 ] >= 0. ? 1. : -1. )
            + shift / d[ ioffd + ll - 1 ] );
          var g = e[ ioffe + ll - 1 ];
          for ( i = ll; i <= m - 1; i ++ ) {
            LaPack1.dlartg( f, g, cosr, sinr, r );
            if ( i > l ) e[ ioffe + i - 2 ] = r.getValue();
            f = cosr.getValue() * d[ ioffd + i - 1 ]
              + sinr.getValue() * e[ ioffe + i - 1 ];
            e[ ioffe + i - 1 ] = cosr.getValue() * e[ ioffe + i - 1 ]
              - sinr.getValue() * d[ ioffd + i - 1 ];
            g = sinr.getValue() * d[ ioffd + i ];
            d[ ioffd + i ] *= cosr.getValue();
            LaPack1.dlartg( f, g, cosl, sinl, r );
            d[ ioffd + i - 1 ] = r.getValue();
            f = cosl.getValue() * e[ ioffe + i - 1 ]
              + sinl.getValue() * d[ ioffd + i ];
            d[ ioffd + i ] = cosl.getValue() * d[ ioffd + i ]
              - sinl.getValue() * e[ ioffe + i - 1 ];
            if ( i < m - 1 ) {
              g = sinl.getValue() * e[ ioffe + i ];
              e[ ioffe + i ] *= cosl.getValue();
            }
            work[ ioffwork + i - ll ] = cosr.getValue();
            work[ ioffwork + i - ll + nm1 ] = sinr.getValue();
            work[ ioffwork + i - ll + nm12 ] = cosl.getValue();
            work[ ioffwork + i - ll + nm13 ] = sinl.getValue();
          } // 140
          e[ ioffe + m - 2 ] = f;
          if ( ncvt > 0 ) {
            LaPack0.dlasr( 'L', 'V', 'F', m - ll + 1, ncvt, work, work,
              Vt, ldvt, ioffwork, ioffwork + n - 1, ioffvt + ll - 1 );
          }
          if ( nru > 0 ) {
            LaPack0.dlasr( 'R', 'V', 'F', nru, m - ll + 1, work, work,
              U, ldu, ioffwork + nm12, ioffwork + nm13,
              ioffu + ( ll - 1 ) * ldu );
          }
          if ( ncc > 0 ) {
            LaPack0.dlasr( 'L', 'V', 'F', m - ll + 1, ncc, work, work,
              C, ldc, ioffwork + nm12, ioffwork + nm13,
              ioffc + ll - 1 );
          }
          if ( Math.abs( e[ ioffe + m - 2 ] ) <= thresh ) {
           e[ ioffe + m - 2 ] = 0.;
          }
        } else {
          f = ( Math.abs( d[ ioffd + m - 1 ] ) - shift )
            * ( ( d[ ioffd + m - 1 ] >= 0. ? 1. : -1. )
            + shift / d[ ioffd + m - 1 ] );
          g = e[ ioffe + m - 2 ];
          for ( i = m; i >= ll + 1; i -- ) {
            LaPack1.dlartg( f, g, cosr, sinr, r );
            if ( i < m ) e[ ioffe + i - 1 ] = r.getValue();
            f = cosr.getValue() * d[ ioffd + i - 1 ]
              + sinr.getValue() * e[ ioffe + i - 2 ];
            e[ ioffe + i - 2 ] = cosr.getValue() * e[ ioffe + i - 2 ]
              - sinr.getValue() * d[ ioffd + i - 1 ];
            g = sinr.getValue() * d[ ioffd + i - 2 ];
            d[ ioffd + i - 2 ] *= cosr.getValue();
            LaPack1.dlartg( f, g, cosl, sinl, r );
            d[ ioffd + i - 1 ] = r.getValue();
            f = cosl.getValue() * e[ ioffe + i - 2 ]
              + sinl.getValue() * d[ ioffd + i - 2 ];
            d[ ioffd + i - 2 ] = cosl.getValue() * d[ ioffd + i - 2 ]
              - sinl.getValue() * e[ ioffe + i - 2 ];
            if ( i < ll + 1 ) {
              g = sinl.getValue() * e[ ioffe + i - 3 ];
              e[ ioffe + i - 3 ] *= cosl.getValue();
            }
            work[ ioffwork + i - ll - 1 ] = cosr.getValue();
            work[ ioffwork + i - ll + nm1 - 1 ] = - sinr.getValue();
            work[ ioffwork + i - ll + nm12 - 1 ] = cosl.getValue();
            work[ ioffwork + i - ll + nm13 - 1 ] = - sinl.getValue();
          } // 150
          e[ ioffe + ll - 1 ] = f;
          if ( Math.abs( e[ ioffe + ll - 1 ] ) <= thresh ) {
            e[ ioffe + ll - 1 ] = 0.;
          }
          if ( ncvt > 0 ) {
            LaPack0.dlasr( 'L', 'V', 'B', m - ll + 1, ncvt, work, work,
              Vt, ldvt, ioffwork + nm12, ioffwork + nm13,
              ioffvt + ll - 1 );
          }
          if ( nru > 0 ) {
            LaPack0.dlasr( 'R', 'V', 'B', nru, m - ll + 1, work, work,
              U, ldu, ioffwork, ioffwork + n - 1,
              ioffu + ( ll - 1 ) * ldu );
          }
          if ( ncc > 0 ) {
            LaPack0.dlasr( 'L', 'V', 'B', m - ll + 1, ncc, work, work,
              C, ldc, ioffwork, ioffwork + n - 1,
              ioffc + ll - 1 );
          }
        }
      }
    }
  } // 160
  if ( ! goto200 ) {
    for ( i = 1; i <= n; i ++ ) {
      if ( d[ ioffd + i - 1 ] < 0. ) {
        d[ ioffd + i - 1 ] = - d[ ioffd + i - 1 ];
      }
      if ( ncvt > 0 ) {
        Blas1.dscal( ncvt, -1., Vt, ldvt, ioffvt + i - 1 );
      }
    } // 170
    for ( i = 1; i <= n - 1; i ++ ) {
      var isub = 1;
      smin = d[ ioffd ];
      for ( var j = 2; j <= n + 1 - i; j ++ ) {
        if ( d[ ioffd + j - 1 ] <= smin ) {
          isub = j;
          smin = d[ ioffd + j - 1 ];
        }
      } // 180
      if ( isub != n + 1 - i ) {
        d[ ioffd + isub - 1 ] = d[ ioffd + n - i ];
        d[ ioffd + n - i ] = smin;
        if ( ncvt > 0 ) {
          Blas1.dswap( ncvt, Vt, ldvt, Vt, ldvt, ioffvt + isub - 1,
            ioffvt + n - i );
        }
        if ( nru > 0 ) {
          Blas1.dswap( nru, U, 1, U, 1, ioffu + ( isub - 1 ) * ldu,
            ioffu + ( n - i ) * ldu );
        }
        if ( ncc > 0 ) {
          Blas1.dswap( ncc, C, ldc, C, ldc, ioffc + isub - 1,
            ioffc + n - i );
        }
      }
    } // 190
    return
  } // 200
  info.setValue( 0 );
  for ( i = 1; i <= n - 1; i ++ ) {
    if ( e[ ioffe + i - 1 ] != 0. ) info.setValue( info.getValue() + 1 );
  } // 210
}
LaPack5.zbdsqr = function( uplo, n, ncvt, nru, ncc, d, e, Vt,
ldvt, U, ldu, C, ldc, work, info ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack5.dlaed0 = function( icompq, qsiz, n, d, e, Q, ldq, Qstore,
ldqs, work, iwork, info, ioffd, ioffe, ioffq, ioffqstore, ioffwork,
ioffiwork ) {
  var smlsiz = 25;
  info.setValue( 0 );
  if ( icompq < 0 || icompq > 2 ) info.setValue( -1 );
  else if ( icompq == 1 && qsiz < Math.max( 0, n ) ) info.setValue( -2 );
  else if ( n < 0 ) info.setValue( -3 );
  else if ( ldq < Math.max( 1, n ) ) info.setValue( -7 );
  else if ( ldqs < Math.max( 1, n ) ) info.setValue( -9 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlaed0', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  iwork[ ioffiwork ] = n;
  var subpbs = 1;
  var tlvls = 0;
  while ( iwork[ ioffiwork + subpbs - 1 ] > smlsiz ) {
    for ( var j = subpbs; j >= 1; j -- ) {
      iwork[ ioffiwork + 2 * j - 1 ] =
        Math.floor( ( iwork[ ioffiwork + j - 1 ] + 1 ) / 2 );
      iwork[ ioffiwork + 2 * j - 2 ] =
        Math.floor( iwork[ ioffiwork + j - 1 ] / 2 );
    }
    tlvls ++;
    subpbs *= 2;
  }
  for ( j = 2; j <= subpbs; j ++ ) {
    iwork[ ioffiwork + j - 1 ] += iwork[ ioffiwork + j - 2 ];
  }
  var spm1 = subpbs - 1;
  for ( var i = 1; i <= spm1; i ++ ) {
    var submat = iwork[ ioffiwork + i - 1 ] + 1;
    var smm1 = submat - 1;
    d[ ioffd + smm1 - 1 ] -= Math.abs( e[ ioffe + smm1 - 1 ] );
    d[ ioffd + submat - 1 ] -= Math.abs( e[ ioffe + smm1 - 1 ] );
  }
  var indxq = 4 * n + 3;
  if ( icompq != 2 ) {
    var temp = Math.log( Number( n ) ) / Math.log( 2. );
    var lgn = Math.floor( temp );
    if ( Math.pow( 2, lgn ) < n ) lgn++;
    if ( Math.pow( 2, lgn ) < n ) lgn++;
    var iprmpt = indxq + n + 1;
    var iperm = iprmpt + n * lgn;
    var iqptr = iperm + n * lgn;
    var igivpt = iqptr + n + 2;
    var igivcl = igivpt + n * lgn;
    var igivnm = 1;
    var iq = igivnm + 2 * n * lgn;
    var iwrem = iq + n * n + 1;
    for ( i = 0; i < subpbs; i ++ ) {
      iwork[ ioffiwork + iprmpt + i - 1 ] = 1;
      iwork[ ioffiwork + igivpt + i - 1 ] = 1;
    }
    iwork[ ioffiwork + iqptr - 1 ] = 1;
  }
  var curr = 0;
  var goto130 = false;
  for ( i = 0; i <= spm1; i ++ ) {
    if ( i == 0 ) {
      submat = 1;
      var matsiz = iwork[ ioffiwork ];
    } else {
      submat = iwork[ ioffiwork + i - 1 ] + 1;
      matsiz = iwork[ ioffiwork + i ] - iwork[ ioffiwork + i - 1 ];
    }
    if ( icompq == 2 ) {
      LaPack2.dsteqr( 'I', matsiz, d, e, Q, ldq, work, info,
        ioffd + submat - 1, ioffe + submat - 1,
        ioffq + submat - 1 + ( submat - 1 ) * ldq, ioffwork );
      if ( info.getValue() != 0 ) {
        goto130 = true;
        break;
      }
    } else {
      LaPack2.dsteqr( 'I', matsiz, d, e, work, matsiz, work, info,
        ioffd + submat - 1, ioffe + submat - 1,
        ioffwork + iq - 2 + iwork[ ioffiwork + iqptr + curr - 1 ],
        ioffwork );
      if ( info.getValue() != 0 ) {
        goto130 = true;
        break;
      }
      if ( icompq == 1 ) {
        Blas3.dgemm( 'N', 'N', qsiz, matsiz, matsiz, 1., Q, ldq, work,
          matsiz, 0., Qstore, ldqs, ioffq + ( submat - 1 ) * ldq,
          ioffwork + iq - 2 + iwork[ ioffiwork + iqptr + curr - 1 ],
          ioffqstore + ( submat - 1 ) * ldqs );
      }
      iwork[ ioffiwork + iqptr + curr ] =
        iwork[ ioffiwork + iqptr + curr - 1 ] + matsiz * matsiz;
      curr ++;
    }
    var k = 1;
    for ( j = submat; j <= iwork[ ioffiwork + i ]; j ++ ) {
      iwork[ ioffiwork + indxq + j - 1 ] = k;
      k ++;
    }
  }
  if ( ! goto130 ) {
    var curlvl = 1;
    while ( subpbs > 1 ) {
      var spm2 = subpbs - 2;
      for ( i = 0; i <= spm2; i += 2 ) {
        if ( i == 0 ) {
          submat = 1;
          matsiz = iwork[ ioffiwork + 1 ];
          var msd2 = iwork[ ioffiwork ];
          var curprb = 0;
        } else {
          submat = iwork[ ioffiwork + i - 1 ] + 1;
          matsiz = iwork[ ioffiwork + i + 1 ]
                 - iwork[ ioffiwork + i - 1 ];
          msd2 = Math.floor( matsiz / 2 );
          curprb ++;
        }
        if ( icompq == 2 ) {
          LaPack4.dlaed1( matsiz, d, Q, ldq, iwork,
            e[ ioffe + submat + msd2 - 2 ], msd2, work, iwork, info,
            ioffd + submat - 1,
            ioffq + submat - 1 + ( submat - 1 ) * ldq,
            ioffiwork + indxq + submat - 1, ioffwork,
            ioffiwork + subpbs );
        } else {
          LaPack4.dlaed7( icompq, matsiz, qsiz, tlvls, curlvl, curprb,
            d, Qstore, ldqs, iwork, e[ ioffe + submat + msd2 - 2 ],
            msd2, work, iwork, iwork, iwork,
            iwork, iwork, work, work, iwork, info, ioffd + submat - 1,
            ioffqstore + ( submat - 1 ) * ldqs,
            ioffiwork + indxq + submat - 1, ioffwork + iq - 1,
            ioffiwork + iqptr - 1, ioffiwork + iprmpt - 1,
            ioffiwork + iperm - 1, ioffiwork + igivpt - 1,
            ioffiwork + igivcl - 1, ioffwork + igivnm - 1,
            ioffwork + iwrem - 1, ioffiwork + subpbs );
        }
        if ( info.getValue() != 0 ) {
          goto130 = true;
          break;
        }
      }
      if ( goto130 ) break;
      subpbs /= 2;
      curlvl ++;
    }
  }
  if ( ! goto130 ) {
    if ( icompq == 1 ) {
      for ( i = 1; i <= n; i ++ ) {
        j = iwork[ ioffiwork + indxq + i - 1 ];
        work[ ioffwork + i - 1 ] = d[ ioffd + j - 1 ];
        Blas1.dcopy( qsiz, Qstore, 1, Q, 1,
          ioffqstore + ( j - 1 ) * ldqs, ioffq + ( i - 1 ) * ldq );
      }
      Blas1.dcopy( n, work, 1, d, 1, ioffwork, ioffd );
    } else if ( icompq == 2 ) {
      for ( i = 1; i <= n; i ++ ) {
        j = iwork[ ioffiwork + indxq + i - 1 ];
        work[ ioffwork + i - 1 ] = d[ ioffd + j - 1 ];
        Blas1.dcopy( n, Q, 1, work, 1, ioffq + ( j - 1 ) * ldq,
          ioffwork + n * i );
      }
      Blas1.dcopy( n, work, 1, d, 1, ioffwork, ioffd );
      LaPack0.dlacpy( 'A', n, n, work, n, Q, ldq, ioffwork + n,
        ioffq );
    } else {
      for ( i = 1; i <= n; i ++ ) {
        j = iwork[ ioffiwork + indxq + i - 1 ];
        work[ ioffwork + i - 1 ] = d[ ioffd + j - 1 ];
      }
      Blas1.dcopy( n, work, 1, d, 1, ioffwork, ioffd );
    }
    return;
  }
  info.setValue( submat * ( n + 1 ) + submat + matsiz - 1 );
}
LaPack5.zlaed0 = function( icompq, qsiz, n, d, e, Q, ldq, Qstore,
ldqs, work, iwork, info ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack5.dlaqr2 = function( wantt, wantz, n, ktop,
kbot, nw, H, ldh, iloz, ihiz, Z, ldz, ns, nd, sr, si, V, ldv, nh, T, ldt,
nv, Wv, ldwv, work, lwork, ioffh, ioffz, ioffsr, ioffsi, ioffv, iofft,
ioffwv, ioffwork) {
  var jw = Math.min( nw, kbot - ktop + 1 );
  var info = new IntReference();
  if ( jw <= 2 ) var lwkopt = 1;
  else {
    LaPack3.dgehrd( jw, 1, jw - 1, T, ldt, work, work, -1, info,
      iofft, ioffwork, ioffwork );
    var lwk1 = Math.round( work[ ioffwork ] );
    LaPack4.dormhr( 'R', 'N', jw, jw, 1, jw - 1, T, ldt, work, V, ldv,
      work, -1, info, iofft, ioffwork, ioffv, ioffwork );
    var lwk2 = Math.round( work[ ioffwork ] );
    lwkopt = jw + Math.max( lwk1, lwk2 );
  }
  if ( lwork == -1 ) {
    work[ ioffwork ] = Number( lwkopt );
    return;
  }
  ns.setValue( 0 );
  nd.setValue( 0 );
  work[ ioffwork ] = 1.;
  if ( ktop > kbot ) return;
  if ( nw < 1 ) return;
  var safminReference = 
    new NumberReference( LaPack0.dlamch( 'Safe minimum' ) );
  var safmaxReference = 
    new NumberReference( 1. / safmin.getValue() );
  LaPack0.dlabad( safmin, safmax );
  var ulp = LaPack0.dlamch( 'Precision' );
  var smlnum = safmin.getValue() * ( Number( n ) / ulp );
  jw = Math.min( nw, kbot - ktop + 1 );
  var kwtop = kbot - jw + 1;
  var s = ( kwtop == ktop ? 0. :
    H[ ioffh + kwtop - 1 + ( kwtop - 2 ) * ldh ] );
  if ( kbot == ktop ) {
    sr[ ioffsr + kwtop - 1 ] =
      H[ ioffh + kwtop - 1 + ( kwtop - 1 ) * ldh ];
    si[ ioffsi + kwtop - 1 ] = 0.;
    ns.setValue( 1 );
    nd.setValue( 0 );
    if ( Math.abs( s ) <= Math.max( smlnum, ulp *
    Math.abs( H[ ioffh + kwtop - 1 + ( kwtop - 1 ) * ldh ] ) ) ) {
      ns.setValue( 0 );
      nd.setValue( 1 );
      if ( kwtop > ktop ) {
        H[ ioffh + kwtop - 1 + ( kwtop - 2 ) * ldh ] = 0.;
      }
    }
    work[ ioffwork ] = 1.;
    return;
  }
  LaPack0.dlacpy( 'U', jw, jw, H, ldh, T, ldt,
    ioffh + kwtop - 1 + ( kwtop - 1 ) * ldh, iofft );
  Blas1.dcopy( jw - 1, H, ldh + 1, T, ldt + 1,
    ioffh + kwtop + ( kwtop - 1 ) * ldh, iofft + 1 );
  LaPack0.dlaset( 'A', jw, jw, 0., 1., V, ldv, ioffv );
  var infqr = new IntReference();
  LaPack2.dlahqr( true, true, jw, 1, jw, T, ldt, sr, si, 1, jw, V, ldv,
    infqr, iofft, ioffsr + kwtop - 1, ioffsi + kwtop - 1, ioffv );
  for ( var j = 1; j <= jw - 3; j ++ ) {
    T[ iofft + j + 1 + ( j - 1 ) * ldt ] = 0.;
    T[ iofft + j + 2 + ( j - 1 ) * ldt ] = 0.;
  } // 10
  if ( jw > 2 ) T[ iofft + jw - 1 + ( jw - 3 ) * ldt ] = 0.;
  ns.setValue( jw );
  var ilst = new IntReference( infqr.getValue() + 1 );
  var ifst = new IntReference( );
  while ( ilst.getValue() <= ns.getValue() ) { // 20
    var bulge = ( ns.getValue() == 1 ? false :
      T[ iofft + ns.getValue() - 1 + ( ns.getValue() - 2 ) * ldt ] != 0. );
    if ( ! bulge ) {
      var foo =
        Math.abs( T[ iofft + ns.getValue() - 1 + ( ns.getValue() - 1 ) * ldt ] );
      if ( foo == 0. ) foo = Math.abs( s );
      if ( Math.abs( s * V[ ioffv + ( ns.getValue() - 1 ) * ldv ] ) <=
      Math.max( smlnum, ulp * foo ) ) {
        ns.setValue( ns.getValue() - 1 );
      } else {
        ifst.setValue( ns.getValue() );
        LaPack4.dtrexc( 'V', jw, T, ldt, V, ldv, ifst, ilst, work,
          info, iofft, ioffv, ioffwork );
        ilst.setValue( ilst.getValue() + 1 );
      }
    } else {
      foo =
        Math.abs( T[ iofft + ns.getValue() - 1 + ( ns.getValue() - 1 ) * ldt ] )
        + Math.sqrt( Math.abs(
          T[ iofft + ns.getValue() - 1 + ( ns.getValue() - 2 ) * ldt ] ) )
        * Math.sqrt( Math.abs(
          T[ iofft + ns.getValue() - 2 + ( ns.getValue() - 1 ) * ldt ] ) );
      if ( foo == 0. ) foo = Math.abs( s );
      if ( Math.max( Math.abs(
      s * V[ ioffv + ( ns.getValue() - 1 ) * ldv ] ),
      Math.abs( s * V[ ioffv + ( ns.getValue() - 2 ) * ldv ] ) ) <=
      Math.max( smlnum, ulp * foo ) ) {
        ns.setValue( ns.getValue() - 2 );
      } else {
        ifst.setValue( ns.getValue() );
        LaPack4.dtrexc( 'V', jw, T, ldt, V, ldv, ifst, ilst, work,
          info, iofft, ioffv, ioffwork );
        ilst.setValue( ilst.getValue() + 2 );
      }
    }
  }
  if ( ns.getValue() == 0 ) s = 0.;
  if ( ns.getValue() < jw ) {
    var sorted = false;
    var i = ns.getValue() + 1;
    while ( ! sorted ) { // 30
      sorted = true;
      var kend = i - 1;
      i = infqr.getValue() + 1;
      if ( i == ns.getValue() ) var k = i + 1;
      else if ( T[ iofft + i + ( i - 1 ) * ldt ] == 0. ) k = i + 1;
      else k = k + 2;
      while ( k <= kend ) { // 40
        var evi = ( k == i + 1 ?
          Math.abs( T[ iofft + i - 1 + ( i - 1 ) * ldt ] ) :
          Math.abs( T[ iofft + i - 1 + ( i - 1 ) * ldt ] )
          + Math.sqrt( Math.abs( T[ iofft + i + ( i - 1 ) * ldt ] ) )
          * Math.sqrt( Math.abs( T[ iofft + i - 1 + i * ldt ] ) ) );
        if ( k == kend ) {
          var evk =
            Math.abs( T[ iofft + k - 1 + ( k - 1 ) * ldt ] );
        } else if ( T[ iofft + k + ( k - 1 ) * ldt ] == 0. ) {
          evk = Math.abs( T[ iofft + k - 1 + ( k - 1 ) * ldt ] );
        } else {
          evk = Math.abs( T[ iofft + k - 1 + ( k - 1 ) * ldt ] )
            + Math.sqrt( Math.abs(
              T[ iofft + k + ( k - 1 ) * ldt ] ) )
            * Math.sqrt( Math.abs( T[ iofft + k - 1 + k * ldt ] ) );
        }
        if ( evi >= evk ) i = k;
        else {
          sorted = false;
          ifst.setValue( i );
          ilst.setValue( k );
          LaPack4.dtrexc( 'V', jw, T, ldt, V, ldv, ifst, ilst, work,
            info, iofft, ioffv, ioffwork );
          i = ( info.getValue() == 0 ? ilst.getValue() : k );
        }
        if ( i == kend ) k = i + 1;
        else if ( T[ iofft + i + ( i - 1 ) * ldt ] == 0. ) k = i + 1;
        else k = i + 2;
      }
    }
  } // 50
  i = jw;
  while ( i >= infqr.getValue() + 1 ) { // 60
    if ( i == infqr.getValue() + 1 ) {
      sr[ ioffsr + kwtop + i - 2 ] =
        T[ iofft + i - 1 + ( i - 1 ) * ldt ];
      si[ ioffsi + kwtop + i - 2 ] = 0.;
      i --;
    } else if ( T[ iofft + i - 1 + ( i - 2 ) * ldt ] == 0. ) {
      sr[ ioffsr + kwtop + i - 2 ] =
        T[ iofft + i - 1 + ( i - 1 ) * ldt ];
      si[ ioffsi + kwtop + i - 2 ] = 0.;
      i --;
    } else {
      var aaReference = 
        new NumberReference( T[ iofft + i - 2 + ( i - 2 ) * ldt ] );
      var ccReference = 
        new NumberReference( T[ iofft + i - 1 + ( i - 2 ) * ldt ] );
      var bbReference = 
        new NumberReference( T[ iofft + i - 2 + ( i - 1 ) * ldt ] );
      var ddReference = 
        new NumberReference( T[ iofft + i - 1 + ( i - 1 ) * ldt ] );
      var rt1rReference =
        new NumberReference( sr[ ioffsr + kwtop + i - 3 ] );
      var rt1iReference =
        new NumberReference( si[ ioffsi + kwtop + i - 3 ] );
      var rt2rReference =
        new NumberReference( sr[ ioffsr + kwtop + i - 2 ] );
      var rt2iReference =
        new NumberReference( si[ ioffsi + kwtop + i - 2 ] );
      var csReference = new NumberReference();
      var snReference = new NumberReference();
      LaPack1.dlanv2( aa, bb, cc, dd, rt1r, rt1i, rt2r, rt2i,
        cs, sn );
      sr[ ioffsr + kwtop + i - 3 ] = rt1r.getValue();
      si[ ioffsi + kwtop + i - 3 ] = rt1i.getValue();
      sr[ ioffsr + kwtop + i - 2 ] = rt2r.getValue();
      si[ ioffsi + kwtop + i - 2 ] = rt2i.getValue();
      i -= 2;
    }
  }
  if ( ns.getValue() < jw || s == 0. ) {
    if ( ns.getValue() > 1 && s != 0. ) {
      Blas1.dcopy( ns.getValue(), V, ldv, work, 1, ioffv, ioffwork );
      var betaReference = 
        new NumberReference( work[ ioffwork ] );
      var tauReference = new NumberReference( );
      LaPack1.dlarfg( ns.getValue(), beta, work, 1, tau, ioffwork + 1 );
      work[ ioffwork ] = 1.;
      LaPack0.dlaset( 'L', jw - 2, jw - 2, 0., 0., T, ldt,
        iofft + 2 );
      LaPack1.dlarf( 'L', ns.getValue(), jw, work, 1, tau.getValue(), T, ldt,
        work, ioffwork, iofft, ioffwork + jw );
      LaPack1.dlarf( 'R', ns.getValue(), ns.getValue(), work, 1, tau.getValue(),
        T, ldt, work, ioffwork, iofft, ioffwork + jw );
      LaPack1.dlarf( 'R', jw, ns.getValue(), work, 1, tau.getValue(), V, ldv,
        work, ioffwork, ioffv, ioffwork + jw );
      LaPack3.dgehrd( jw, 1, ns.getValue(), T, ldt, work, work,
        lwork - jw, info, iofft, ioffwork, ioffwork + jw );
    }
    if ( kwtop > 1 ) {
      H[ ioffh + kwtop - 1 + ( kwtop - 2 ) * ldh ] = s * V[ ioffv ];
    }
    LaPack0.dlacpy( 'U', jw, jw, T, ldt, H, ldh, iofft,
      ioffh + kwtop - 1 + ( kwtop - 1 ) * ldh );
    Blas1.dcopy( jw - 1, T, ldt + 1, H, ldh + 1,
      iofft + 1, ioffh + kwtop + ( kwtop - 1 ) * ldh );
    if ( ns.getValue() > 1 && s != 0. ) {
      LaPack4.dormhr( 'R', 'N', jw, ns.getValue(), 1, ns.getValue(), T, ldt,
        work, V, ldv, work, lwork - jw, info,
        iofft, ioffwork, ioffv, ioffwork + jw );
    }
    var ltop = ( wantt ? 1 : ktop );
    for ( var krow = ltop; krow <= ktop - 1; krow += nv ) {
      var kln = Math.min( nv, kwtop - krow );
      Blas3.dgemm( 'N', 'N', kln, jw, jw, 1., H, ldh, V, ldv, 0.,
        Wv, ldwv, ioffh + krow - 1 + ( kwtop - 1 ) * ldh, ioffv,
        ioffwv );
      LaPack0.dlacpy( 'A', kln, jw, Wv, ldwv, H, ldh, ioffwv,
        ioffh + krow - 1 + ( kwtop - 1 ) * ldh );
    } // 70
    if ( wantt ) {
      for ( var kcol = kbot + 1; kcol <= n; kcol += nh ) {
        kln = Math.min( nh, n - kcol + 1 );
        Blas3.dgemm( 'C', 'N', jw, kln, jw, 1., V, ldv, H, ldh, 0.,
          T, ldt, ioffv, ioffh + kwtop - 1 + ( kcol - 1 ) * ldh,
          iofft );
        LaPack0.dlacpy( 'A', jw, kln, T, ldt, H, ldh, iofft,
          ioffh + kwtop - 1 + ( kcol - 1 ) * ldh );
      } // 80
    }
    if ( wantz ) {
      for ( krow = iloz; krow <= ihiz; krow += nv ) {
        kln = Math.min( nv, ihiz - krow + 1 );
        Blas3.dgemm( 'N', 'N', kln, jw, jw, 1., Z, ldz, V, ldv, 0.,
          Wv, ldwv, ioffz + krow - 1 + ( kwtop - 1 ) * ldz,
          ioffv, ioffwv );
        LaPack0.dlacpy( 'A', kln, jw, Wv, ldwv, Z, ldz, ioffwv,
          ioffz + krow - 1 + ( kwtop - 1 ) * ldz );
      } // 90
    }
  }
  nd.setValue( jw - ns.getValue() );
  ns.setValue( ns.getValue() - infqr.getValue() );
  work[ ioffwork ] = Number( lwkopt );
}
//*************************************************************************
LaPack5.dstedc = function( compz, n, d, e, Z, ldz, work, lwork,
iwork, liwork, info, ioffd, ioffe, ioffz, ioffwork, ioffiwork ) {
  throw new Error("not tested");
  info.setValue( 0 );
  var lquery = ( lwork == -1 || liwork == -1 );
  if ( compz.charAt(0).toUpperCase() == 'N' ) var icompz = 0;
  else if ( compz.charAt(0).toUpperCase() == 'V' ) icompz = 1;
  else if ( compz.charAt(0).toUpperCase() == 'I' ) icompz = 2;
  else icompz = -1;
  if ( icompz < 0 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  else if ( ldz < 1 || ( icompz > 0 && ldz < Math.max( 1, n ) ) ) {
    info.setValue( -6 );
  }
  if ( info.getValue() == 0 ) {
    var smlsiz = LaPack0.ilaenv( 9, 'dstedc', ' ', 0, 0, 0, 0 );
    if ( n <= 1 || icompz == 0 ) {
      var liwmin = 1;
      var lwmin = 1;
    } else if ( n <= smlsiz ) {
      liwmin = 1;
      lwmin = 2 * ( n - 1 );
    } else {
      var lgn = Math.round( Math.log( Number( n ) ) / Math.low( 2. ) );
      if ( Math.pow( 2, lgn ) < n ) lgn ++;
      if ( Math.pow( 2, lgn ) < n ) lgn ++;
      if ( icompz == 1 ) {
        lwmin = 1 + 3 * n + 2 * n * lgn + 4 * n * n;
        liwmin = 6 + 6 * n + 5 * n * lgn;
      } else {
        lwmin = 1 + 4 * n + n * n;
        liwmin = 3 + 5 * n;
      }
    }
    work[ ioffwork ] = lwmin;
    iwork[ ioffiwork ] = liwmin;
    if ( lwork < lwmin && ! lquery ) info.setValue( -8 );
    else if ( liwork < liwmin && ! lquery ) info.setValue( -10 );
  }
  if ( info.getValue() != 0 ) {
      Blas2.xerbla( 'dstedc', - info.getValue() );
      return;
  } else if ( lquery ) return;
  if ( n == 0 ) return;
  if ( n == 1 ) {
    if ( icompz != 0 ) Z[ ioffz ] = 1.;
    return;
  }
  var goto50 = false;
  if ( icompz == 0 ) {
    LaPack2.dsterf( n, d, e, info, ioffd, ioffe );
    goto50 = true;
  }
  if ( ! goto50 ) {
    if ( n <= smlsiz ) {
      LaPack2.dsteqr( compz, n, d, e, Z, ldz, work, info,
        ioffd, ioffe, ioffz, ioffwork );
    } else {
      var storez = ( icompz == 1 ? 1 + n * n : 1 );
      if ( icompz == 2 ) {
        LaPack0.dlaset( 'Full', n, n, 0., 1., Z, ldz, ioffz );
      }
      var orgnrm = LaPack1.dlanst( 'M', n, d, e, ioffd, ioffe );
      if ( orgnrm == 0. ) goto50 = true;
      if ( ! goto50 ) {
        var eps = LaPack0.dlamch( 'Epsilon' );
        var start = 1;
        while ( start <= n ) { // 10
          var finish = start;
          while ( finish < n ) { // 20
            var tiny =
              eps * Math.sqrt( Math.abs( d[ ioffd + finish - 1 ] ) )
              * Math.sqrt( Math.abs( d[ ioffd + finish ] ) );
            if ( Math.abs( e[ ioffe + finish - 1 ] ) > tiny ) finish ++;
            else break;
          }
          var m = finish - start + 1;
          if ( m == 1 ) {
            start = finish + 1;
            continue; // goto 10
          }
          if ( m > smlsiz ) {
            orgnrm = LaPack1.dlanst( 'M', m, d, e, ioffd + start - 1,
              ioffe + start - 1 );
            LaPack1.dlascl( 'G', 0, 0, orgnrm, 1., m, 1, d, m, info,
              ioffd + start - 1 );
            LaPack1.dlascl( 'G', 0, 0, orgnrm, 1., m - 1, 1, e, m - 1,
              info, ioffe + start - 1 );
            var strtrw = ( icompz == 1 ? 1 : start );
            LaPack4.dlaed0( icompz, n, m, d, e, Z, ldz, work, n, work,
              iwork, info, ioffd + start - 1, ioffe + start - 1,
              ioffz + strtrw - 1 + ( start - 1 ) * ldz, ioffwork,
              ioffwork + storez - 1, ioffiwork );
            if ( info.getValue() != 0 ) {
              info.setValue( ( info.getValue() / ( m + 1 ) + start - 1 )
                * ( n + 1 ) + ( info.getValue() % ( m + 1 ) ) + start - 1 );
              goto50 = true;
              break;
            } 
            LaPack1.dlascl( 'G', 0, 0, 1., orgnrm, m, 1, d, m, info,
              ioffd + start - 1 );
          } else {
            if ( icompz  == 1 ) {
              LaPack2.dsteqr( 'I', m, d, e, work, m, work, info,
                ioffd + start - 1, ioffe + start - 1, ioffwork,
                ioffwork + m * m );
              LaPack0.dlacpy( 'A', n, m, Z, ldz, work, n,
                ioffz + ( start - 1 ) * ldz, ioffwork + storez - 1 );
              Blas3.dgemm( 'N', 'N', n, m, m, 1., work, n, work, m,
                0., Z, ldz, ioffwork + storez - 1, ioffwork,
                ioffz + ( start - 1 ) * ldz );
            } else if ( icompz == 2 ) {
              LaPack2.dsteqr( 'I', m, d, e, Z, ldz, work, info,
                ioffd + start - 1, ioffe + start - 1,
                ioffz + start - 1 + ( start - 1 ) * ldz, ioffwork );
            } else {
              LaPack2.dsterf( m, d, e, info, ioffd + start - 1,
                ioffe + start - 1 );
            }
            if ( info.getValue() != 0 ) {
              info.setValue( start * ( n + 1 ) + finish );
              goto50 = true;
              break;
            }
          }
          start = finish + 1;
        } // end while
      }
      if ( ! goto50 ) {
        if ( m != n ) {
          if ( icompz == 0 ) LaPack0.dlasrt( 'I', n, d, info, ioffd );
          else {
            for ( var ii = 2; ii <= n; ii ++ ) {
              var i = ii - 1;
              var k = k;
              var p = d[ ioffd + i - 1 ];
              for ( var j = ii; j <= n; j ++ ) {
                if ( d[ ioffd + j - 1 ] < p ) {
                  k = j;
                  p = d[ ioffd + j - 1 ];
                }
              } // 30
              if ( k != i ) {
                d[ ioffd + k - 1 ] - d[ ioffd + i - 1 ];
                d[ ioffd + i - 1 ] = p;
                Blas1.dswap( n, Z, 1, Z, 1, ioffz + ( i - 1 ) * ldz,
                  ioffz + ( k - 1 ) * ldz );
              }
            } // 40
          }
        }
      }
    }
  }
  work[ ioffwork ] = lwmin;
  iwork[ ioffiwork ] = liwmin;
}
LaPack5.zstedc = function( compz, n, d, e, Z, ldz, work, lwork,
iwork, liwork, info ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack5.dstemr = function( jobz, range, n, d, e, vl, vu, il, iu,
m, w, Z, ldz, nzc, isuppz, tryrac, work, lwork, iwork, liwork, info, ioffd,
ioffe, ioffw, ioffz, ioffisuppz, ioffwork, ioffiwork ) {
  throw new Error("not programmed: should have been");
  var minrgp = 1.e-3;
  var wantz = ( jobz.charAt(0).toUpperCase() == 'V' );
  var alleig = ( range.charAt(0).toUpperCase() == 'A' );
  var valeig = ( range.charAt(0).toUpperCase() == 'V' );
  var indeig = ( range.charAt(0).toUpperCase() == 'I' );
  var lquery = ( lwork == -1 || liwork == -1 );
  var zquery = ( nzc == -1 );
  if ( wantz ) {
    var lwmin = 18 * n;
    var liwmin = 10 * n;
  } else {
    lwmin = 12 * n;
    liwmin = 8 * n;
  }
  var wl = 0.;
  var wu = 0.;
  var iil = 0;
  var iiu = 0;
  if ( valeig ) {
    wl = vl;
    wu = vu;
  } else if ( indeig ) {
    iil = il;
    iiu = iu;
  }
  info.setValue( 0 );
  if ( ! ( wantz || jobz.charAt(0).toUpperCase() == 'N' ) ) {
    info.setValue( -1 );
  } else if ( ! ( alleig || valeig || indeig ) ) info.setValue( -2 );
  else if ( n < 0 ) info.setValue( -3 );
  else if ( valeig && n > 0 && wu <= wl ) info.setValue( -7 );
  else if ( indeig && ( iil < 1 || iil > n ) ) info.setValue( -8 );
  else if ( indeig && ( iiu < iil || iiu > n ) ) info.setValue( -9 );
  else if ( ldz < 1 || ( wantz && ldz < n ) ) info.setValue( -13 );
  else if ( lwork < lwmin && ! lquery ) info.setValue( -17 );
  else if ( liwork < liwmin && ! lquery ) info.setValue( -19 );
  var safmin = LaPack0.dlamch( 'Safe minimum' );
  var eps = LaPack0.dlamch( 'Precision' );
  var smlnum = safmin / eps;
  var bignum = 1. / smlnum;
  var rmin = Math.sqrt( smlnum );
  var rmax = Math.min( Math.sqrt( bignum ),
    1. / Math.sqrt( Math.sqrt( safmin ) ) );
  if ( info.getValue() == 0 ) {
    work[ ioffwork ] = lwmin;
    iwork[ ioffiwork ] = liwmin;
    if ( wantz && alleig ) var nzcmin = n;
    else if ( wantz && valeig ) {
      var eigcnt = new IntReference( );
      var itmp = new IntReference( );
      var itmp2 = new IntReference( );
      LaPack0.dlarrc( 'T', n, vl, vu, d, e, safmin, eigcnt, itmp,
        itmp2, info, ioffd, ioffe );
      nzcmin = eigcnt.getValue();
    } else if ( wantz && indeig ) nzcmin = iiu - iil + 1;
    else nzcmin = 0;
    if ( zquery && info.getValue() == 0 ) z[ ioffz ] = nzcmin;
    else if ( nzc < nzcmin && ! zquery ) info.setValue( -14 );
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dstemr', - info.getValue() );
    return;
  } else if ( lquery || zquery ) return;
  m.setValue( 0 );
  if ( n == 0 ) return;
  if ( n == 1 ) {
    if ( alleig || indeig ) {
      m.setValue( 1 );
      w[ ioffw ] = d[ ioffd ];
    } else {
      if ( wl < d[ ioffd ] && wu >= d[ ioffd ] ) {
        m.setValue( 1 );
        w[ ioffw ] = d[ ioffd ];
      }
    }
    if ( wantz && ! zquery ) {
      z[ ioffz ] = 1.;
      isuppz[ ioffisuppz ] = 1;
      isuppz[ ioffisuppz + 1 ] = 1;
    }
    return;
  }
  var r1Reference = new NumberReference();
  var r2Reference = new NumberReference();
  var csReference = new NumberReference();
  var snReference = new NumberReference();
  if ( n == 2 ) {
    if ( ! wantz ) {
      LaPack0.dlae2( d[ ioffd ], e[ ioffe ], d[ ioffd + 1 ],
        r1, r2 );
    } else if ( wantz && ! zquery ) {
      LaPack0.dlaev2( d[ ioffd ], e[ ioffe ], d[ ioffd + 1 ],
        r1, r2, cs, sn );
    }
    if ( alleig || ( valeig && r2.getValue() > wl && r2.getValue() <= wu ) ||
    ( indeig && iil == 1 ) ) {
      m.setValue( m.getValue() + 1 );
      w[ ioffw + m.getValue() - 1 ] = r2.getValue();
      if ( wantz && ! zquery ) {
        z[ ioffz + ( m.getValue() - 1 ) * ldz ] = - sn.getValue();
        z[ ioffz + 1 + ( m.getValue() - 1 ) * ldz ] = cs.getValue();
        if ( sn.getValue() != 0. ) {
          if ( cs.getValue() != 0. ) {
            isuppz[ ioffisuppz + 2 * m.getValue() - 2 ] = 1;
            isuppz[ ioffisuppz + 2 * m.getValue() - 1 ] = 2;
          } else {
            isuppz[ ioffisuppz + 2 * m.getValue() - 2 ] = 1;
            isuppz[ ioffisuppz + 2 * m.getValue() - 1 ] = 1;
          }
        } else {
          isuppz[ ioffisuppz + 2 * m.getValue() - 2 ] = 2;
          isuppz[ ioffisuppz + 2 * m.getValue() - 1 ] = 2;
        }
      }
    }
    if ( alleig || ( valeig && r1.getValue() > wl && r1.getValue() <= wu ) ||
    ( indeig && iiu == 2 ) ) {
      m.setValue( m.getValue() + 1 );
      w[ ioffw + m.getValue() - 1 ] = r1.getValue();
      if ( wantz && ! zquery ) {
        z[ ioffz + ( m.getValue() - 1 ) * ldz ] = cs.getValue();
        z[ ioffz + 1 + ( m.getValue() - 1 ) * ldz ] = sn.getValue();
        if ( sn.getValue() != 0. ) {
          if ( cs.getValue() != 0. ) {
            isuppz[ ioffisuppz + 2 * m.getValue() - 2 ] = 1;
            isuppz[ ioffisuppz + 2 * m.getValue() - 1 ] = 2;
          } else {
            isuppz[ ioffisuppz + 2 * m.getValue() - 2 ] = 1;
            isuppz[ ioffisuppz + 2 * m.getValue() - 1 ] = 1;
          }
        } else {
          isuppz[ ioffisuppz + 2 * m.getValue() - 2 ] = 2;
          isuppz[ ioffisuppz + 2 * m.getValue() - 1 ] = 2;
        }
      }
    }
    return;
  }
  var indgrs = 1;
  var inderr = 2 * n + 1;
  var indgp = 3 * n + 1;
  var indd = 4 * n + 1;
  var inde2 = 5 * n + 1;
  var indwrk = 6 * n + 1;
  var iinspl = 1;
  var iindbl = n + 1;
  var iindw = 2 * n + 1;
  var iindwk = 3 * n + 1;
  var scale = 1.;
  var tnrm = LaPack1.dlanst( 'M', n, d, e, ioffd, ioffe );
  if ( tnrm > 0. && tnrm < rmin ) scale = rmin / tnrm;
  else if ( tnrm > rmax ) scale = rmax / tnrm;
  if ( scale != 1. ) {
    Blas1.dscal( n, scale, d, 1, ioffd );
    Blas1.dscal( n - 1, scale, e, 1, ioffe );
    tnrm *= scale;
    if ( valeig ) {
      wl *= scale;
      wu *= scale;
    }
  }
  var iinfo = new IntReference();
  if ( tryrac ) LaPack1.dlarrr( n, d, e, iinfo, ioffd, ioffe );
  else iinfo.setValue( -1 );
  if ( iinfo.getValue() == 0 ) var thresh = eps;
  else {
    thresh = - eps;
    tryrac = false;
  }
  if ( tryrac ) {
    Blas1.dcopy( n, d, 1, work, 1, ioffd, ioffwork + indd - 1 );
  }
  for ( j = 1; j <= n - 1; j ++ ) {
    work[ ioffwork + inde2 + j - 2 ] =
      Math.pow( e[ ioffe + j - 1 ], 2 );
  } // 5
  if ( ! wantz ) {
    var rtol1 = 4. * eps;
    var rtol2 = 4. * eps;
  } else {
    rtol1 = Math.sqrt( eps );
    rtol2 = Math.max( Math.sqrt( eps ) * 5.e-3, 4. * eps );
  }
  var wlrefReference = new NumberReference( wl );
  var wurefReference = new NumberReference( wu );
  var nsplit = new IntReference();
  var pivminReference = new NumberReference( wu );
  LaPack4.dlarre( range, n, wlref,
    wuref, iil, iiu, d, e, work,
    rtol1, rtol2, thresh, nsplit,
    iwork, m, w, work, work,
    iwork, iwork, work, pivmin,
    work, iwork, iinfo, ioffd, ioffe,
    ioffwork + inde2 - 1, ioffiwork + iinspl - 1, ioffw,
      ioffwork + inderr - 1, ioffwork + indgp - 1,
    ioffiwork + iindbl - 1, ioffiwork + iindw - 1,
      ioffwork + indgrs - 1, ioffwork + indwrk - 1,
    ioffiwork + iindwk - 1 );
  wl = wlref.getValue();
  wu = wuref.getValue();
  if ( iinfo.getValue() != 0 ) {
    info.setValue( 10 + Math.abs( iinfo.getValue() ) );
    return;
  }
  if ( wantz ) {
    LaPack2.dlarrv( n, wl, wu, d,
      e, pivmin.getValue(), iwork, m.getValue(), 1, m.getValue(),
      minrgp, rtol1, rtol2, w, work,
      work, iwork, iwork, work, Z, ldz,
      isuppz, work, iwork, iinfo, ioffd,
      ioffe, ioffiwork + iinspl - 1, ioffw, ioffwork + inderr - 1,
        ioffwork + indgp - 1,
      ioffiwork + iindbl - 1, ioffiwork + iindw - 1,
        ioffwork + indgrs - 1, ioffz,
      ioffisuppz, ioffwork + indwrk - 1, ioffiwork + iindwk - 1 );
    if ( iinfo.getValue() != 0 ) {
      info.setValue( 20 + Math.abs( iinfo.getValue() ) );
      return;
    }
  } else {
    for ( var j = 1; j <= m.getValue(); j ++ ) {
      itmp = iwork[ ioffiwork + iindbl + j - 2 ];
      w[ ioffw + j - 1 ] +=
        e[ iwork[ ioffiwork + iinspl + itmp - 2 ] ];
    } // 20
  }
  if ( tryrac ) {
    var ibegin = 1;
    var wbegin = 1;
    for ( var jblk = 1;
    jblk <= iwork[ ioffiwork + iindbl + m.getValue() - 2 ]; jblk ++ ) {
      var iend = iwork[ ioffiwork + iinspl + jblk - 2 ];
      var in2 = iend - ibegin + 1;
      var wend = wbegin - 1;
      while ( wend < m.getValue() ) { // 36
        if ( iwork[ ioffiwork + iindbl + wend - 1 ] == jblk ) {
          wend ++;
        } else break;
      }
      if ( wend < wbegin ) {
        ibegin = iend + 1;
        continue;
      }
      var offset = iwork[ ioffiwork + iindw + wbegin - 2 ] - 1;
      var ifirst = iwork[ ioffiwork + iindw + wbegin - 2 ];
      var ilast = iwork[ ioffiwork + iindw + wend - 2 ];
      rtol2 = 4. * eps;
      LaPack0.dlarrj( in2, work, work, ifirst, ilast, rtol2, offset,
        w, work, work, iwork, pivmin.getValue(), tnrm, iinfo,
        ioffwork + indd + ibegin - 2, ioffwork + inde2 + ibegin - 2,
        ioffw + wbegin - 1, ioffwork + inderr + wbegin - 2,
        ioffwork + indwrk - 1, ioffiwork + iindwk - 1 );
      ibegin = iend + 1;
      wbegin = wend + 1;
    } // 39
  }
  if ( scale != 1. ) Blas1.dscal( m.getValue(), 1. / scale, w, 1, ioffw );
  if ( nsplit.getValue() > 1 ) {
    if ( ! wantz ) {
      LaPack0.dlasrt( 'I', m.getValue(), w, iinfo, ioffw );
      if ( iinfo.getValue() != 0 ) {
        info.setValue( 3 );
        return;
      }
    } else {
      for ( j = 1; j <= m.getValue() - 1; j ++ ) {
        var i = 0;
        var tmp = w[ ioffw + j - 1 ];
        for ( var jj = j + 1; jj <= m.getValue(); jj ++ ) {
          if ( w[ ioffw + jj - 1 ] < tmp ) {
            i = jj;
            tmp = w[ ioffw + jj - 1 ];
          }
        } // 50
        if ( i != 0 ) {
          w[ ioffw + i - 1 ] = w[ ioffw + j - 1 ];
          w[ ioffw + j - 1 ] = tmp;
          if ( wantz ) {
            Blas1.dswap( n, Z, 1, Z, 1, ioffz + ( i - 1 ) * ldz,
              ioffz + ( j - 1 ) * ldz );
            itmp = isuppz[ ioffisuppz + 2 * i - 2 ];
            isuppz[ ioffisuppz + 2 * i - 2 ] =
              isuppz[ ioffisuppz + 2 * j - 2 ];
            isuppz[ ioffisuppz + 2 * j - 2 ] = itmp;
            itmp = isuppz[ ioffisuppz + 2 * i - 1 ];
            isuppz[ ioffisuppz + 2 * i - 1 ] =
              isuppz[ ioffisuppz + 2 * j - 1 ];
            isuppz[ ioffisuppz + 2 * j - 1 ] = itmp;
          }
        }
      } // 60
    }
  }
  work[ ioffwork ] = lwmin;
  iwork[ ioffiwork ] = liwmin;
}
//*************************************************************************
LaPack5.dsyev = function( jobz, uplo, n, A, lda, w, work, lwork,
info, ioffa, ioffw, ioffwork) {
  var wantz = ( jobz.charAt(0).toUpperCase() == 'V' );
  var lower = ( uplo.charAt(0).toUpperCase() == 'L' );
  var lquery = ( lwork == -1 );
  info.setValue( 0 );
  if ( ! wantz && jobz.charAt(0).toUpperCase() != 'N' ) {
    info.setValue( -1 );
  } else if ( ! lower && uplo.charAt(0).toUpperCase() != 'U' ) {
    info.setValue( -2 );
  } else if ( n < 0 ) info.setValue( -3 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -5 );
  if ( info.getValue() == 0 ) {
    var nb = LaPack0.ilaenv( 1, 'dsytrd', uplo, n, -1, -1, -1 );
    var lwkopt = Math.max( 1, ( nb + 2 ) * n );
    work[ ioffwork ] = lwkopt;
    if ( lwork < Math.max( 1, 3 * n - 1 ) && ! lquery ) info.setValue( -8 );
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dsyev', - info.getValue() );
    return;
  } else if ( lquery ) return;
  if ( n == 0 ) return;
  if ( n == 1 ) {
    w[ ioffw ] = A[ ioffa ];
    work[ ioffwork ] = 2;
    if ( wantz ) A[ ioffa ] = 1.;
    return;
  }
  var safmin = LaPack0.dlamch( 'Safe minimum' );
  var eps = LaPack0.dlamch( 'Precision' );
  var smlnum = safmin / eps;
  var bignum = 1. / smlnum;
  var rmin = Math.sqrt( smlnum );
  var rmax = Math.sqrt( bignum );
  var anrm = LaPack1.dlansy( 'M', uplo, n, A, lda, work,
    ioffa, ioffwork );
  var iscale = 0;
  if ( anrm > 0. && anrm < rmin ) {
    iscale = 1;
    var sigma = rmin / anrm;
  } else if ( anrm > rmax ) {
    iscale = 1;
    sigma = rmax / anrm;
  }
  if ( iscale == 1 ) {
    LaPack1.dlascl( uplo, 0, 0, 1., sigma, n, n, A, lda, info, ioffa );
  }
  var inde = 1;
  var indtau = inde + n;
  var indwrk = indtau + n;
  var llwork = lwork - indwrk + 1;
  var iinfo = new IntReference();
  LaPack3.dsytrd( uplo, n, A, lda, w, work, work, work, llwork, iinfo,
    ioffa, ioffw, ioffwork + inde - 1, ioffwork + indtau - 1,
    ioffwork + indwrk - 1 );
  if ( ! wantz ) {
    LaPack2.dsterf( n, w, work, info, ioffw, ioffwork + inde - 1 );
  } else {
    LaPack4.dorgtr( uplo, n, A, lda, work, work, llwork, iinfo, ioffa,
      ioffwork + indtau - 1, ioffwork + indwrk - 1 );
    LaPack2.dsteqr( jobz, n, w, work, A, lda, work, info, ioffw,
      ioffwork + inde - 1, ioffa, ioffwork + indtau - 1 );
  }
  if ( iscale == 1 ) {
    var imax = ( info.getValue() == 0 ? n : info.getValue() - 1 ); 
    Blas1.dscal( imax, 1. / sigma, w, 1, ioffw );
  }
  work[ ioffwork ] = lwkopt;
}
//*************************************************************************
LaPack5.dsyevx = function( jobz, range, uplo, n, A, lda, vl, vu,
il, iu, abstol, m, w, Z, ldz, work, lwork, iwork, ifail, info, ioffa,
ioffw, ioffz, ioffwork, ioffiwork, ioffifail) {
  var lower = ( uplo.charAt(0).toUpperCase() == 'L' );
  var wantz = ( jobz.charAt(0).toUpperCase() == 'V' );
  var alleig = ( range.charAt(0).toUpperCase() == 'A' );
  var valeig = ( range.charAt(0).toUpperCase() == 'V' );
  var indeig = ( range.charAt(0).toUpperCase() == 'I' );
  var lquery = ( lwork == -1 );
  info.setValue( 0 );
  if ( ! ( wantz || jobz.charAt(0).toUpperCase() == 'N' ) ) {
    info.setValue( -1 );
  } else if ( ! ( alleig || valeig || indeig ) ) info.setValue( -2 );
  else if ( ! ( lower || uplo.charAt(0).toUpperCase() == 'U' ) ) {
    info.setValue( -3 );
  } else if ( n < 0 ) info.setValue( -4 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -6 );
  else {
    if ( valeig ) {
      if ( n >0 && vu <= vl ) info.setValue( -8 );
    } else if ( indeig ) {
      if ( il < 1 || il > Math.max( 1, n ) ) info.setValue( -9 );
      else if ( iu < Math.min( n, il ) || iu > n ) info.setValue( -10 );
    }
  }
  if ( info.getValue() == 0 ) {
    if ( ldz < 1 || ( wantz && ldz < n ) ) info.setValue( -15 );
  }
  if ( info.getValue() == 0 ) {
    if ( n <= 1 ) {
      var lwkmin = 1;
      work[ ioffwork ] = lwkmin;
    } else {
      lwkmin = 8 * n;
      var nb = LaPack0.ilaenv( 1, 'dsytrd', uplo, n, -1, -1, -1 );
      nb = Math.max( nb,
       LaPack0.ilaenv( 1, 'dormtr', uplo, n, -1, -1, -1 ) );
      var lwkopt = Math.max( lwkmin, ( nb + 3 ) * n );
      work[ ioffwork ] = lwkopt;
    }
    if ( lwork < lwkmin && ! lquery ) info.setValue( -17 );
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dsyevx', - info.getValue() );
    return;
  } else if ( lquery ) return;
  m.setValue( 0 );
  if ( n == 0 ) {
    work[ ioffwork ] = 1;
    return;
  }
  if ( n == 1 ) {
    if ( alleig || indeig ) {
      m.setValue( 1 );
      w[ ioffw ] = A[ ioffa ];
    } else {
      if ( vl < A[ ioffa ] && vu >= A[ ioffa ] ) {
        m.setValue( 1 );
        w[ ioffw ] = A[ ioffa ];
      }
    }
    if ( wantz ) Z[ ioffz ] = 1.;
    return;
  }
  var safmin = LaPack0.dlamch( 'Safe minimum' );
  var eps = LaPack0.dlamch( 'Precision' );
  var smlnum = safmin / eps;
  var bignum = 1. / smlnum;
  var rmin = Math.sqrt( smlnum );
  var rmax = Math.min( Math.sqrt( bignum ),
    1. / Math.sqrt( Math.sqrt( safmin ) ) );;
  var iscale = 0;
  var abstll = abstol;
  if ( valeig ) {
    var vll = vl;
    var vuu = vu;
  }
  var anrm = LaPack1.dlansy( 'M', uplo, n, A, lda, work, ioffa,
    ioffwork );
  if ( anrm > 0. && anrm < rmin ) {
    iscale = 1;
    var sigma = rmin / anrm;
  } else if ( anrm > rmax ) {
    iscale = 1;
    sigma = rmax / anrm;
  }
  if ( iscale == 1 ) {
    if ( lower ) {
      for ( var j = 1; j <= n; j ++ ) {
        Blas1.dscal( n - j + 1, sigma, A, 1,
          ioffa + j - 1 + ( j - 1 ) * lda );
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        Blas1.dscal( j, sigma, A, 1, ioffa + ( j - 1 ) * lda );
      }
    }
    if ( abstol > 0. ) abstll = abstol * sigma;
    if ( valeig ) {
      vll = vl * sigma;
      vuu = vu * sigma;
    }
  }
  var indtau = 1;
  var inde = indtau + n;
  var indd = inde + n;
  var indwrk = indd + n;
  var llwork = lwork - indwrk + 1;
  var iinfo = new IntReference();
  LaPack3.dsytrd( uplo, n, A, lda, work, work, work, work, llwork,
    iinfo, ioffa, ioffwork + indd - 1, ioffwork + inde - 1,
    ioffwork + indtau - 1, ioffwork + indwrk - 1 );
  var test = false;
  if ( indeig ) {
    if ( il == 1 && iu == n ) test = true;
  }
  var lopt = 3 * n + work[ ioffwork + indwrk - 1 ];
  var goto40 = false;
  if ( ( alleig || test ) && abstol <= 0. ) {
    Blas1.dcopy( n, work, 1, w, 1, ioffwork + indd - 1, ioffw );
    var indee = indwrk + 2 * n;
    if ( ! wantz ) {
      Blas1.dcopy( n - 1, work, 1, work, 1, ioffwork + inde - 1,
        ioffwork + indee - 1 );
      LaPack2.dsterf( n, w, work, info, ioffw, ioffwork + indee - 1 );
    } else {
      LaPack0.dlacpy( 'A', n, n, A, lda, Z, ldz, ioffa, ioffz );
      LaPack4.dorgtr( uplo, n, Z, ldz, work, work, llwork, iinfo,
        ioffz, ioffwork + indtau, ioffwork + indwrk );
      Blas1.dcopy( n - 1, work, 1, work, 1, ioffwork + inde - 1,
        ioffwork + indee - 1 );
      LaPack2.dsteqr( jobz, n, w, work, Z, ldz, work, info,
        ioffw, ioffwork + indee - 1, ioffz, ioffwork + indwrk - 1 );
      if ( info.getValue() == 0 ) {
        for ( var i = 1; i <= n; i ++ ) {
          ifail[ ioffifail + i - 1 ] = 0;
        }
      }
    }
    if ( info.getValue() == 0 ) {
      m.setValue( n );
      goto40 = true;
    }
    info.setValue( 0 );
  }
  if ( ! goto40 ) { // line 320
    var order = ( wantz ? 'B' : 'E' );
    var indibl = 1;
    var indisp = indibl + n;
    var indiwo = indisp + n;
    var nsplit = new IntReference();
    LaPack1.dstebz( range, order, n, vll, vuu, il, iu, abstll, work,
      work, m, nsplit, w, iwork, iwork, work, iwork, info,
      ioffwork + indd - 1, ioffwork + inde - 1, ioffw,
      ioffiwork + indibl - 1, ioffiwork + indisp - 1,
      ioffwork + indwrk - 1, ioffiwork + indiwo - 1 );
    if ( wantz ) {
      LaPack2.dstein( n, work, work, m.getValue(), w, iwork, iwork, Z, ldz,
        work, iwork, ifail, info, ioffwork + indd - 1,
        ioffwork + inde - 1, ioffw, ioffiwork + indibl - 1,
        ioffiwork + indisp - 1, ioffz, ioffwork + indwrk - 1,
        ioffiwork + indiwo - 1, ioffifail );
      var indwkn = inde;
      var llwrkn = lwork - indwkn + 1;
      LaPack4.dormtr( 'L', uplo, 'N', n, m.getValue(), A, lda, work, Z, ldz,
        work, llwrkn, iinfo, ioffa, ioffwork + indtau - 1, ioffz,
        ioffwork + indwkn - 1 );
    }
  }
  if ( iscale == 1 ) {
    var imax = ( info.getValue() == 0 ? m.getValue() : info.getValue() - 1 );
    Blas1.dscal( imax, 1. / sigma, w, 1, ioffw );
  }
  if ( wantz ) {
    for ( j = 1; j <= m.getValue() - 1; j ++ ) {
      i = 0;
      var tmp1 = w[ ioffw + j - 1 ];
      for ( var jj = j + 1; jj <= m.getValue(); jj ++ ) {
        if ( w[ ioffw + jj - 1 ] < tmp1 ) {
          i = jj;
          tmp1 = w[ ioffw + jj - 1 ];
        }
      }
      if ( i != 0 ) {
        var itmp1 = iwork[ ioffiwork + indibl + i - 2 ];
        w[ ioffw + i - 1 ] = w[ ioffw + j - 1 ];
        iwork[ ioffiwork + indibl + i - 2 ] =
          iwork[ ioffiwork + indibl + j - 2 ];
        w[ ioffw + j - 1 ] = tmp1;
        iwork[ ioffiwork + indibl + j - 2 ] = itmp1;
        Blas1.dswap( n, Z, 1, Z, 1, ioffz + ( i - 1 ) * ldz,
          ioffz + ( j - 1 ) * ldz );
        if ( info.getValue() != 0 ) {
          itmp1 = ifail[ ioffifail + i - 1 ];
          ifail[ ioffifail + i - 1 ] = ifail[ ioffifail + j - 1 ];
          ifail[ ioffifail + j - 1 ] = itmp1;
        }
      }
    }
  }
  work[ ioffwork ] = lwkopt;
}
//*************************************************************************
LaPack5.dtgex2 = function( wantq, wantz, n, A, lda, B, ldb, Q,
ldq, Z, ldz, j1, n1, n2, work, lwork, info, ioffa, ioffb, ioffq, ioffz,
ioffwork ) {
  throw new Error("not programmed: generalized eigenvalue");
}
//*************************************************************************
LaPack5.dtgsyl = function( trans, ijob, m, n, A, lda, B, ldb, C,
ldc, D, ldd, E, lde, F, ldf, difReference, scaleReference, work, lwork,
iwork, info, ioffa, ioffb, ioffc, ioffd, ioffe, iofff, ioffwork,
ioffiwork ) {
  throw new Error("not programmed: generalized eigenvalue");
}
//*************************************************************************
LaPack5.dtrsen = function( job, compq, select, n, T, ldt, Q,
ldq, wr, wi, m, sReference, sepReference, work, lwork, iwork, liwork, info,
ioffselect, iofft, ioffq, ioffwr, ioffwi, ioffwork, ioffiwork ) {
  var isave = new Array( 3 );
  var wantbh = ( job.charAt(0).toUpperCase() == 'B' );
  var wants = ( job.charAt(0).toUpperCase() == 'E' );
  var wantsp = ( job.charAt(0).toUpperCase() == 'V' );
  var wantq = ( compq.charAt(0).toUpperCase() == 'V' );
  info.setValue( 0 );
  var lquery = ( lwork == -1 );
  if ( job.charAt(0).toUpperCase() != 'N' && ! wants && ! wantsp ) {
    info.setValue( -1 );
  } else if ( compq.charAt(0).toUpperCase() != 'N' && ! wantq ) {
    info.setValue( -2 );
  } else if ( n < 0 ) info.setValue( -4 );
  else if ( ldt < Math.max( 1, n ) ) info.setValue( -6 );
  else if ( ldq < 1 || ( wantq && ldq < n ) ) info.setValue( -8 );
  else {
    m.setValue( 0 );
    var pair = false;
    for ( var k = 1; k <= n; k ++ ) {
      if ( pair ) pair = false;
      else {
        if ( k < n ) {
          if ( T[ iofft + k + ( k - 1 ) * ldt ] == 0. ) {
            if ( select[ ioffselect + k - 1 ] ) m.setValue( m.getValue() + 1 );
          } else {
            pair = true;
            if ( select[ ioffselect + k - 1 ]
            || select[ ioffselect + k ] ) {
              m.setValue( m.getValue() + 2 );
            }
          }
        } else {
          if ( select[ ioffselect + n - 1 ] ) m.setValue( m.getValue() + 1 );
        }
      }
    } // 10
    var n1 = m.getValue();
    var n2 = n - m.getValue();
    var nn = n1 * n2;
    if ( wantsp ) {
      var lwmin = Math.max( 1, 2 * nn );
      var liwmin = Math.max( 1, nn );
    } else if ( job.charAt(0).toUpperCase() == 'N' ) {
      lwmin = Math.max( 1, n );
      liwmin = 1;
    } else if ( job.charAt(0).toUpperCase() == 'E' ) {
      lwmin = Math.max( 1, nn );
      liwmin = 1;
    }
    if ( lwork < lwmin && ! lquery ) info.setValue( -15 );
    else if ( liwork < liwmin && ! lquery ) info.setValue( -17 );
  }
  if ( info.getValue() == 0 ) {
    work[ ioffwork ] = lwmin;
    iwork[ ioffiwork ] = liwmin;
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dtrsen', - info.getValue() );
    return;
  } else if ( lquery ) return;
  var goto40 = false;
  if ( m.getValue() == n || m.getValue() == 0 ) {
    if ( wants ) s.getValue() = 1.;
    if ( wantsp ) {
      sep.setValue(
        LaPack1.dlange( '1', n, n, T, ldt, work, iofft, ioffwork ) );
    }
    goto40 = true;
  }
  if ( !goto40 ) {
    var ks = 0;
    pair = false;
    for ( k = 1; k < n; kk ++ ) {
      if ( pair ) pair = false;
      else {
        var swap = select[ ioffselect + k - 1 ];
        if ( k < n ) {
          if ( T[ iofft + k + ( k - 1 ) * ldt ] != 0. ) {
            pair = true;
            swap = swap || select[ ioffselect + k ];
          }
        }
        if ( swap ) {
          ks ++;
          var ierr = new IntReference( 0 );
          var kk = k;
          if ( k != ks ) {
            var kkref = new IntReference( kk );
            var ksref = new IntReference( ks );
            LaPack4.dtrexc( compq, n, T, ldt, Q, ldq, kkref, ksref,
              work, ierr, iofft, ioffq, ioffwork );
            kk = kkref.getValue();
            ks = ksref.getValue();
          }
          if ( ierr.getValue() == 1 || ierr.getValue() == 2 ) {
            info.setValue( 1 );
            if ( wants ) s.setValue( 0. );
            if ( wantsp ) sep.setValue( 0. );
            goto40 = true;
            break;
          }
          if ( pair ) ks ++;
        }
      }
    } // 20
  }
  if ( !goto40 ) {
    if ( wants ) {
      LaPack0.dlacpy( 'F', n1, n2, T, ldt, work, n1, iofft + n1 * ldt,
        ioffwork );
      var scaleReference = new NumberReference( ); 
      LaPack2.dtrsyl( 'N', 'N', -1, n1, n2, T, ldt, T, ldt, work, n1,
        scale, ierr, iofft, iofft + n1 + n1 * ldt, ioffwork );
      var rnorm = LaPack1.dlange( 'F', n1, n2, work, n1, work,
        ioffwork, ioffwork );
      s.setValue( ( rnorm == 0. ? 1. : scale.getValue()
        / ( Math.sqrt( scale.getValue() * scale.getValue() / rnorm + rnorm )
        * Math.sqrt( rnorm ) ) ) );
    }
    if ( wantsp ) {
      var estReference = new NumberReference ( 0. );
      var kase = new IntReference( 0 );
      while ( true ) { // 30
        LaPack0.dlacn2( nn, work, work, iwork, est, kase, isave,
          ioffwork + nn, ioffwork, ioffiwork, 0 );
        if ( kase.getValue() != 0 ) {
          if ( kase.getValue() == 1 ) {
            LaPack2.dtrsyl( 'N', 'N', -1, n1, n2, T, ldt, T, ldt, work,
              n1, scale, ierr, iofft, iofft + n1 + n1 * ldt,
              ioffwork );
          } else {
            LaPack2.dtrsyl( 'T', 'T', -1, n1, n2, T, ldt, T, ldt, work,
              n1, scale, ierr, iofft, iofft + n1 + n1 * ldt,
              ioffwork );
          }
        } else break;
      }
      sep.setValue( scale.getValue() / est.getValue() );
    }
  } // 40
  for ( k = 1; k <= n; k ++ ) {
    wr[ ioffwr + k - 1 ] = T[ iofft + k - 1 + ( k - 1 ) * ldt ];
    wi[ ioffwi + k - 1 ] = 0.;
  } // 50
  for ( k = 1; kk <= n - 1; k ++ ) {
    if ( T[ iofft + k + ( k - 1 ) * ldt ] != 0. ) {
      wi[ ioffwi + k - 1 ] =
        Math.sqrt( Math.abs( T[ iofft + k - 1 + k * ldt ] ) )
        * Math.sqrt( Math.abs( T[ iofft + k + ( k - 1 ) * ldt ] ) );
      wi[ ioffwi + k ] = - wi[ ioffwi + k - 1 ];
    }
  } // 60
  work[ ioffwork ] = lwmin;
  iwork[ ioffiwork ] = liwmin;
}
LaPack5.ztrsen = function( job, compq, select, n, T, ldt, Q,
ldq, wr, wi, m, sReference, sepReference, work, lwork, iwork, liwork,
info ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack5.dtrsna = function( job, howmny, select, n, T, ldt, Vl,
ldvl, Vr, ldvr, s, sep, mm, m, work, ldwork, iwork, info, ioffselect,
iofft, ioffvl, ioffvr, ioffs, ioffsep, ioffwork, ioffiwork) {
  var isave = new Array( 3 );
  var dummy = new Array( 1 );
  var wantbh = ( job.charAt(0).toUpperCase() == 'B' );
  var wants = ( job.charAt(0).toUpperCase() == 'E' ) || wantbh;
  var wantsp =
    ( job.charAt(0).toUpperCase() == 'V' ) || wantbh;
  var somcon = ( howmny.charAt(0).toUpperCase() == 'S' );
  info.setValue( 0 );
  if ( ! wants && ! wantsp ) info.setValue( -1 );
  else if ( howmny.charAt(0).toUpperCase() != 'A' && ! somcon ) {
    info.setValue( -2 );
  } else if ( n < 0 ) info.setValue( -4 );
  else if ( ldt < Math.max( 1, n ) ) info.setValue( -6 );
  else if ( ldvl < 1 || ( wants && ldvl < n ) ) info.setValue( -8 );
  else if ( ldvr < 1 || ( wants && ldvr < n ) ) info.setValue( -10 );
  else {
    if ( somcon ) {
      m.setValue( 0 );
      var pair = false;
      for ( var k = 1; k <= n; k ++ ) {
        if ( pair ) pair = false;
        else {
          if ( k < n ) {
            if ( T[ iofft + k + ( k - 1 ) * ldt ] == 0. ) {
              if ( select[ ioffselect + k - 1 ] ) {
                m.setValue( m.getValue() + 1 );
              }
            } else {
              pair = true;
              if ( select[ ioffselect + k - 1 ]
              || select[ ioffselect + k ] ) {
                m.setValue( m.getValue() + 2 );
              }
            }
          } else if ( select[ ioffselect + n - 1 ] ) {
            m.setValue( m.getValue() + 1 );
          }
        }
      } // 10
    } else m.setValue( n );
    if ( mm < m.getValue() ) info.setValue( -13 );
    else if ( ldwork < 1 || ( wantsp && ldwork < n ) ) {
      info.setValue( -16 );
    }
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dtrsna', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  if ( n == 1 ) {
    if ( somcon ) {
      if ( ! select[ ioffselect ] ) return;
    }
    if ( wants ) s[ ioffs ] = 1.;
    if ( wantsp ) sep[ ioffsep ] = Math.abs( T[ iofft ] );
  }
  var eps = LaPack0.dlamch( 'P' );
  var smlnumReference =
    new NumberReference( LaPack0.dlamch( 'S' ) / eps );
  var bignumReference =
    new NumberReference( 1. / smlnum.getValue() );
  LaPack0.dlabad( smlnum, bignum );
  var ks = 0;
  pair = false;
  for ( k = 1; k <= n; k ++ ) { // ==> 60
    if ( pair ) {
      pair = false;
      continue; // goto 60
    } else {
      if ( k < n ) pair = ( T[ iofft + k + ( k - 1 ) * ldt ] != 0. );
    }
    if ( somcon ) {
      if ( pair ) {
        if ( ! select[ ioffselect + k - 1 ]
        && ! select[ ioffselect + k ] ) {
          continue; // goto 60
        }
      } else if ( ! select[ ioffselect + k - 1 ] ) continue; // goto 60
    }
    ks ++;
    if ( wants ) {
      if ( ! pair ) {
        var prod = Blas1.ddot( n, Vr, 1, Vl, 1,
          ioffvr + ( ks - 1 ) * ldvr, ioffvl + ( ks - 1 ) * ldvl );
        var rnrm =
          Blas1.dnrm2( n, Vr, 1, ioffvr + ( ks - 1 ) * ldvr );
        var lnrm =
          Blas1.dnrm2( n, Vl, 1, ioffvl + ( ks - 1 ) * ldvl );
        s[ ioffs + ks - 1 ] = Math.abs( prod ) / ( rnrm * lnrm );
      } else {
        var prod1 = Blas1.ddot( n, Vr, 1, Vl, 1,
          ioffvr + ( ks - 1 ) * ldvr, ioffvl + ( ks - 1 ) * ldvl );
        prod1 += Blas1.ddot( n, Vr, 1, Vl, 1, ioffvr + ks * ldvr,
          ioffvl + ks * ldvl );
        var prod2 =
          Blas1.ddot( n, Vl, 1, Vr, 1, ioffvl + ( ks - 1 ) * ldvl,
          ioffvr + ks * ldvr );
        prod2 -= Blas1.ddot( n, Vl, 1, Vr, 1, ioffvl + ks * ldvl, 
          ioffvr + ( ks - 1 ) * ldvr );
        rnrm = LaPack0.dlapy2(
          Blas1.dnrm2( n, Vr, 1, ioffvr + ( ks - 1 ) * ldvr ),
          Blas1.dnrm2( n, Vr, 1, ioffvr + ks * ldvr ) );
        lnrm = LaPack0.dlapy2(
          Blas1.dnrm2( n, Vl, 1, ioffvl + ( ks - 1 ) * ldvl ),
          Blas1.dnrm2( n, Vl, 1, ioffvl + ks * ldvl ) );
        var cond = LaPack0.dlapy2( prod1, prod2 )
          / ( rnrm * lnrm );
        s[ ioffs + ks - 1 ] = cond;
        s[ ioffs + ks ] = cond;
      }
    }
    if ( wantsp ) {
      LaPack0.dlacpy( 'Full', n, n, T, ldt, work, ldwork, iofft,
        ioffwork );
      var ifst = new IntReference( k );
      var ilst = new IntReference( 1 );
      var ierr = new IntReference( );
      LaPack4.dtrexc( 'No Q', n, work, ldwork, dummy, 1, ifst, ilst,
        work, ierr, ioffwork, 0, ioffwork + n * ldwork );
      var scaleReference = new NumberReference( );
      var estReference = new NumberReference( );
      if ( ierr.getValue() == 1 || ierr.getValue() == 2 ) {
        scale.setValue( 1. );
        est.setValue( bignum.getValue() );
      } else {
        if ( work[ ioffwork + 1 ] == 0. ) {
          for ( var i = 2; i <= n; i ++ ) {
            work[ ioffwork + i - 1 + ( i - 1 ) * ldwork ] -=
              work[ ioffwork ];
          } // 20
          var n2 = 1;
          var nn = n - 1;
        } else {
          var mu =
            Math.sqrt( Math.abs( work[ ioffwork + ldwork ] ) )
            * Math.sqrt( Math.abs( work[ ioffwork + 1 ] ) );
          var delta =
            LaPack0.dlapy2( mu, work[ ioffwork + 1 ] );
          var cs = mu / delta;
          var sn = - work[ ioffwork + 1 ] / delta;
          for ( var j = 3; j <= n; j ++ ) {
            work[ ioffwork + 1 + ( j - 1 ) * ldwork ] *= cs;
            work[ ioffwork + j - 1 + ( j - 1 ) * ldwork ] -=
              work[ ioffwork ];
          } // 30
          work[ ioffwork + 1 + ldwork ] = 0.;
          work[ ioffwork + n * ldwork ] = 2. * mu;
          for ( i = 2; i <= n - 1; i ++ ) {
            work[ ioffwork + i - 1 + n * ldwork ] =
              sn * work[ ioffwork + i * ldwork ];
          } // 40
          n2 = 2;
          nn = 2 * ( n -1 );
        }
        est.setValue( 0. );
        var kase = new IntReference( 0 );
        while ( true ) { // 50
          LaPack0.dlacn2( nn, work, work, iwork, est, kase, isave,
            ioffwork + ( n + 1 ) * ldwork,
            ioffwork + ( n + 3 ) * ldwork, ioffiwork, 0 );
          if ( kase.getValue() != 0 ) {
            var dumm = Number.POSITIVE_INFINITY;
            if ( kase.getValue() == 1 ) {
              if ( n2 == 1 ) {
                LaPack2.dlaqtr( true, true, n - 1, work, ldwork, dummy,
                  dumm, scale, work, work, ierr,
                  ioffwork + 1 + ldwork, 0,
                  ioffwork + ( n + 3 ) * ldwork,
                  ioffwork + ( n + 5 ) * ldwork );
              } else {
                LaPack2.dlaqtr( true, false, n - 1, work, ldwork,
                  work, mu, scale, work, work, ierr,
                  ioffwork + 1 + ldwork, ioffwork + n * ldwork,
                  ioffwork + ( n + 3 ) * ldwork,
                  ioffwork + ( n + 5 ) * ldwork );
              }
            } else {
              if ( n2 == 1 ) {
                LaPack2.dlaqtr( false, true, n - 1, work, ldwork,
                  dummy, dumm, scale, work, work, ierr,
                  ioffwork + 1 + ldwork, 0,
                  ioffwork + ( n + 3 ) * ldwork,
                  ioffwork + ( n + 5 ) * ldwork );
              } else {
                LaPack2.dlaqtr( false, false, n - 1, work, ldwork,
                  work, mu, scale, work, work, ierr,
                  ioffwork + 1 + ldwork, ioffwork + n + ldwork,
                  ioffwork + ( n + 3 ) * ldwork,
                  ioffwork + ( n + 5 ) * ldwork );
              }
            }
            continue;
          } else break;
        }
      }
      sep[ ioffsep + ks - 1 ] =
        scale.getValue() / Math.max( est.getValue(), smlnum.getValue() );
      if ( pair ) sep[ ioffsep + ks ] = sep[ ioffsep + ks - 1 ];
    }
    if ( pair ) ks ++;
  } // 60
}
