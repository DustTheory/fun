function LaPack6() {
}
//**************************************************************************
LaPack6.dlaqr4 = function( wantt, wantz, n, ilo, ihi, H, ldh, wr,
wi, iloz, ihiz, Z, ldz, work, lwork, info, ioffh, ioffwr, ioffwi, ioffz,
ioffwork ) {
  var ntiny = 11;
  var kexnw = 5;
  var kexsh = 6;
  var wilk1 = 0.75;
  var wilk2 = -0.4375;
  var zdum = new Array( 1 );
  info.setValue( 0 );
  if ( n == 0 ) {
    work[ ioffwork ] = 1.;
    return;
  }
  if ( n <= ntiny ) {
    var lwkopt = 1;
    if ( lwork != -1 ) {
      LaPack2.dlahqr( wantt, wantz, n, ilo, ihi, H, ldh, wr, wi, iloz,
        ihiz, Z, ldz, info, ioffh, ioffwr, ioffwi, ioffz );
    }
  } else {
    info.setValue( 0 );
    var jbcmpz = ( wantt ? 'S' : 'E' );
    jbcmpz = jbcmpz + ( wantz ? 'V' : 'N' );
    var nwr =
      LaPack0.ilaenv( 13, 'dlaqr4', jbcmpz, n, ilo, ihi, lwork );
    nwr = Math.max( 2, nwr );
    nwr = Math.min( Math.min( ihi - ilo + 1, ( n - 1 ) / 3 ), nwr );
    var nsr = 
      LaPack0.ilaenv( 15, 'dlaqr4', jbcmpz, n, ilo, ihi, lwork );
    nsr = Math.min( nsr, Math.min( ( n + 6 ) / 9, ihi - ilo ) );
    nsr = Math.max( 2, nsr - nsr % 2 );
    var ls = new IntReference();
    var ld = new IntReference();
    LaPack5.dlaqr2( wantt, wantz, n, ilo, ihi, nwr + 1, H, ldh, iloz,
      ihiz, Z, ldz, ls, ld, wr, wi, H, ldh, n, H, ldh, n, H, ldh,
      work, -1, ioffh, ioffz, ioffwr, ioffwi, ioffh, ioffh, ioffh,
      ioffwork );
    lwkopt = Math.max( 3 * nsr / 2, Math.round( work[ ioffwork ] ) );
    if ( lwork == -1 ) {
      work[ ioffwork ] = Number( lwkopt );
      return;
    }
    var nmin =
      LaPack0.ilaenv( 12, 'dlaqr4', jbcmpz, n, ilo, ihi, lwork );
    nmin = Math.max( ntiny, nmin );
    var nibble =
      LaPack0.ilaenv( 14, 'dlaqr4', jbcmpz, n, ilo, ihi, lwork );
    nibble = Math.max( 0, nibble );
    var kacc22 =
      LaPack0.ilaenv( 16, 'dlaqr4', jbcmpz, n, ilo, ihi, lwork );
    kacc22 = Math.max( 0, kacc22 );
    kacc22 = Math.min( 2, kacc22 );
    var nwmax = Math.min( ( n - 1 ) / 3, lwork / 2 );
    var nw = nwmax;
    var nsmax = Math.min( ( n + 6 ) / 9, 2 * lwork / 3 );
    nsmax -= nsmax % 2;
    var ndfl = 1;
    var itmax = Math.max( 30, 2 * kexsh )
      * Math.max( 10, ihi - ilo + 1 );
    var kbot = ihi;
    var goto90 = false;
    for ( var it = 1; it <= itmax; it ++ ) { // ==> 80
      if ( kbot < ilo ) {
        goto90 = true;
        break;
      }
      var goto20 = false;
      for ( var k = kbot; k >= ilo + 1; k -- ) {
        if ( H[ ioffh + k - 1 + ( k - 2 ) * ldh ] == 0. ) {
          goto20 = true;
          break;
        }
      } // 10
      if ( ! goto20 ) k = ilo;
      var ktop = k;
      var nh = kbot - ktop + 1;
      var nwupbd = Math.min( nh, nwmax );
      nw = ( ndfl < kexnw ? Math.min( nwupbd, nwr ) :
        Math.min( nwupbd, 2 * nw ) );
      if ( nw < nwmax ) {
        if ( nw >= nh - 1 ) nw = nh;
        else {
          var kwtop = kbot - nw + 1;
          if ( Math.abs( H[ ioffh + kwtop - 1 + ( kwtop - 2 ) * ldh ] )
          > Math.abs( H[ ioffh + kwtop - 2 + ( kwtop - 3 ) * ldh ] ) )
          {
            nw ++;
          }
        }
      }
      if ( ndfl < kexnw ) var ndec = -1;
      else if ( ndec >= 0 || nw >= nwupbd ) {
        ndec ++;
        if ( nw - ndec < 2 ) ndec = 0;
        nw -= ndec;
      }
      var kv = n - nw + 1;
      var kt = nw + 1;
      var nho = ( n - nw - 1 ) - kt + 1;
      var kwv = nw + 2;
      var nve = ( n - nw ) - kwv + 1;
      LaPack5.dlaqr2( wantt, wantz, n, ktop, kbot, nw, H, ldh, iloz,
        ihiz, Z, ldz, ls, ld, wr, wi, H, ldh, nho, H, ldh, nve,
        H, ldh, work, lwork, ioffh, ioffz, ioffwr, ioffwi,
        ioffh + kv - 1, ioffh + kv - 1 + ( kt - 1 ) * ldh,
        ioffh + kwv - 1, ioffwork );
      kbot -= ld.getValue();
      var ks = kbot - ls.getValue() + 1;
      if ( ld.getValue() == 0 || ( 100 * ld.getValue() <= nw * nibble &&
      kbot - ktop + 1 > Math.min( nmin, nwmax ) ) ) {
        var ns = Math.min( Math.min( nsmax, nsr ),
          Math.max( 2, kbot - ktop ) );
        ns -= ns % 2;
        var aaReference = new NumberReference( );
        var bbReference = new NumberReference( );
        var ccReference = new NumberReference( );
        var ddReference = new NumberReference( );
        var rt1rReference = new NumberReference( );
        var rt1iReference = new NumberReference( );
        var rt2rReference = new NumberReference( );
        var rt2iReference = new NumberReference( );
        var csReference = new NumberReference();
        var snReference = new NumberReference();
        if ( ndfl % kexsh == 0 ) {
          ks = kbot - ns + 1;
          for ( var i = kbot; i >= Math.max( ks + 1, ktop + 2);
          i -= 2 ) {
            var ssReference = new NumberReference(
              Math.abs( H[ ioffh + i - 1 + ( i - 2 ) * ldh ] )
              + Math.abs( H[ ioffh + i - 2 + ( i - 3 ) * ldh ] ) );
            aa.setValue( wilk1 * ss.getValue()
              + H[ ioffh + i - 1 + ( i - 1 ) * ldh ] );
            bb.setValue( ss.getValue() );
            cc.setValue( wilk2 * ss.getValue() );
            dd.setValue( aa.getValue() );
            rt1r.setValue( wr[ ioffwr + i - 2 ] );
            rt1i.setValue( wi[ ioffwi + i - 2 ] );
            rt2r.setValue( wr[ ioffwr + i - 1 ] );
            rt2i.setValue( wi[ ioffwi + i - 1 ] );
            LaPack1.dlanv2( aa, bb, cc, dd, rt1r, rt1i, rt2r, rt2i,
              cs, sn );
            wr[ ioffwr + i - 2 ] = rt1r.getValue();
            wi[ ioffwi + i - 2 ] = rt1i.getValue();
            wr[ ioffwr + i - 1 ] = rt2r.getValue();
            wi[ ioffwi + i - 1 ] = rt2i.getValue();
          } // 30
          if ( ks == ktop ) {
            wr[ ioffwr + ks ] = H[ ioffh + ks + ks * ldh ];
            wi[ ioffwi + ks ] = 0.;
            wr[ ioffwr + ks - 1 ] = wr[ ioffwr + ks ];
            wi[ ioffwi + ks - 1 ] = wi[ ioffwi + ks ];
          }
        } else {
          if ( kbot - ks + 1 <= ns / 2 ) {
            ks = kbot - ns + 1;
            kt = n - ns + 1;
            LaPack0.dlacpy( 'A', ns, ns, H, ldh, H, ldh,
              ioffh + ks - 1 + ( ks - 1 ) * ldh, ioffh + kt - 1 );
            var inf = new IntReference();
            LaPack2.dlahqr( false, false, ns, 1, ns, H, ldh, wr, wi,
              1, 1, zdum, 1, inf, ioffh + kt - 1, ioffwr + ks - 1,
              ioffwi + ks - 1, 0 );
            ks += inf.getValue();
            if ( ks >= kbot ) {
              aa.setValue( H[ ioffh + kbot - 2 + ( kbot - 2 ) * ldh ] );
              cc.setValue( H[ ioffh + kbot - 1 + ( kbot - 2 ) * ldh ] );
              bb.setValue( H[ ioffh + kbot - 2 + ( kbot - 1 ) * ldh ] );
              dd.setValue( H[ ioffh + kbot - 1 + ( kbot - 1 ) * ldh ] );
              rt1r.setValue( wr[ ioffwr + kbot - 2 ] );
              rt1i.setValue( wi[ ioffwi + kbot - 2 ] );
              rt2r.setValue( wr[ ioffwr + kbot - 1 ] );
              rt2i.setValue( wi[ ioffwi + kbot - 1 ] );
              LaPack1.dlanv2( aa, bb, cc, dd, rt1r, rt1i, rt2r, rt2i,
                cs, sn );
              wr[ ioffwr + kbot - 2 ] = rt1r.getValue();
              wi[ ioffwi + kbot - 2 ] = rt1i.getValue();
              wr[ ioffwr + kbot - 1 ] = rt2r.getValue();
              wi[ ioffwi + kbot - 1 ] = rt2i.getValue();
              ks = kbot - 1;
            }
          }
          if ( kbot - ks + 1 > ns ) {
            var sorted = false;
            for ( k = kbot; k >= ks + 1; k -- ) {
              if ( sorted ) break;
              sorted = true;
              for ( i = ks; i <= k - 1; i ++ ) {
                if ( Math.abs( wr[ ioffwr + i - 1 ] )
                + Math.abs( wi[ ioffwi + i - 1 ] ) <
                Math.abs( wr[ ioffwr + i ] )
                + Math.abs( wi[ ioffwi + i ] ) ) {
                  sorted = false;
                  var swap = wr[ ioffwr + i - 1 ];
                  wr[ ioffwr + i - 1 ] = wr[ ioffwr + i ];
                  wr[ ioffwr + i ] = swap;
                  swap = wi[ ioffwi + i - 1 ];
                  wi[ ioffwi + i - 1 ] = wi[ ioffwi + i ];
                  wi[ ioffwi + i ] = swap;
                }
              } // 40
            } // 50
          }
          for ( i = kbot; i >= ks + 2; i -= 2 ) {
            if ( wi[ ioffwi + i - 1 ] != - wi[ ioffwi + i - 2 ] ) {
              swap = wr[ ioffwr + i - 1 ];
              wr[ ioffwr + i - 1 ] = wr[ ioffwr + i - 2 ];
              wr[ ioffwr + i - 2 ] = wr[ ioffwr + i - 3 ];
              wr[ ioffwr + i - 3 ] = swap;
              swap = wi[ ioffwi + i - 1 ];
              wi[ ioffwi + i - 1 ] = wi[ ioffwi + i - 2 ];
              wi[ ioffwi + i - 2 ] = wi[ ioffwi + i - 3 ];
              wi[ ioffwi + i - 3 ] = swap;
            }
          } // 70
        }
        if ( kbot - ks + 1 == 2 ) {
          if ( wi[ ioffwi + kbot - 1 ] == 0. ) {
            if ( Math.abs( wr[ ioffwr + kbot - 1 ]
            - H[ ioffh + kbot - 1 + ( kbot - 1 ) * ldh ] ) <
            Math.abs( wr[ ioffwr + kbot - 2 ]
            - H[ ioffh + kbot - 1 + ( kbot - 1 ) * ldh ] ) ) {
              wr[ ioffwr + kbot - 2 ] = wr[ ioffwr + kbot - 1 ];
            } else {
              wr[ ioffwr + kbot - 1 ] = wr[ ioffwr + kbot - 2 ];
            }
          }
        }
        ns = Math.min( ns, kbot - ks + 1 );
        ns -= ns % 2;
        ks = kbot - ns + 1;
        var kdu = 3 * ns - 3;
        var ku = n - kdu + 1;
        var kwh = kdu + 1;
        nho = ( n - kdu + 1 - 4 ) - ( kdu + 1 ) + 1;
        kwv = kdu + 4;
        nve = n - kdu - kwv + 1;
        LaPack2.dlaqr5( wantt, wantz, kacc22, n, ktop, kbot, ns,
          wr, wi, H, ldh, iloz, ihiz, Z, ldz, work, 3, H, ldh, nve,
          H, ldh, nho, H, ldh, ioffwr + ks - 1, ioffwi + ks - 1,
          ioffh, ioffz, ioffwork, ioffh + ku - 1, ioffh + kwv - 1,
          ioffh + ku - 1 + ( kwh - 1 ) * ldh );
      }
      ndfl = ( ld.getValue() > 0 ? 1 : ndfl + 1 );
    } // 80
    if ( ! goto90 ) info.setValue( kbot );
  }
  work[ ioffwork ] = Number( lwkopt );
}
//**************************************************************************
LaPack6.dlasdq = function( uplo, sqre, n, ncvt, nru, ncc, d, e,
Vt, ldvt, U, ldu, C, ldc, work, info, ioffd, ioffe, ioffvt, ioffu, ioffc,
ioffwork ) {
  info.setValue( 0 );
  var iuplo = 0;
  if ( uplo.charAt(0).toUpperCase() == 'U' ) iuplo = 1;
  if ( uplo.charAt(0).toUpperCase() == 'L' ) iuplo = 2;
  if ( iuplo == 0 ) info.setValue( -1 );
  else if ( sqre < 0 || sqre > 1 ) info.setValue( -2 );
  else if ( n < 0 ) info.setValue( -3 );
  else if ( ncvt < 0 ) info.setValue( -4 );
  else if ( nru < 0 ) info.setValue( -5 );
  else if ( ncc < 0 ) info.setValue( -6 );
  else if ( ( ncvt == 0 && ldvt < 1 ) ||
  ( ncvt > 0 && ldvt > Math.max( 1, n ) ) ) {
    info.setValue( -10 );
  } else if ( ldu < Math.max( 1, nru ) ) info.setValue( -12 );
  else if ( ( ncc == 0  && ldc < 1 ) ||
  ( ncc > 0 && ldc < Math.max( 1, n ) ) ) {
    info.setValue( - 14 );
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlasdq', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  var rotate = ( ncvt > 0 ) || ( nru > 0 ) || ( ncc > 0 );
  var np1 = n + 1;
  var sqre1 = sqre;
  var csReference = new NumberReference();
  var snReference = new NumberReference();
  var rReference = new NumberReference();
  if ( iuplo == 1 && sqre1 == 1 ) {
    for ( var i = 1; i <= n - 1; i ++ ) {
      LaPack1.dlartg( d[ ioffd + i - 1 ], e[ ioffe + i - 1 ],
        cs, sn, r );
      d[ ioffd + i - 1 ] = r.getValue();
      e[ ioffe + i - 1 ] = sn.getValue() * d[ ioffd + i ];
      d[ ioffd + i ] = cs.getValue() * d[ ioffd + i ];
      if ( rotate ) {
        work[ ioffwork + i - 1 ] = cs.getValue();
        work[ ioffwork + n + i - 1 ] = sn.getValue();
      }
    } // 10
    LaPack1.dlartg( d[ ioffd + n - 1 ], e[ ioffe + n - 1 ],
      cs, sn, r );
    d[ ioffd + n - 1 ] = r.getValue();
    e[ ioffe + n - 1 ] = 0.;
    if ( rotate ) {
      work[ ioffwork + n - 1 ] = cs.getValue();
      work[ ioffwork + n + n - 1 ] = sn.getValue();
    }
    iuplo = 2;
    sqre1 = 0;
    if ( ncvt > 0 ) {
      LaPack0.dlasr( 'L', 'V', 'F', np1, ncvt, work, work, Vt, ldvt,
        ioffwork, ioffwork + np1 - 1, ioffvt );
    }
  }
  if ( iuplo == 2 ) {
    for ( i = 1; i <= n - 1; i ++ ) {
      LaPack1.dlartg( d[ ioffd + i - 1 ], e[ ioffe + i - 1 ],
        cs, sn, r );
      d[ ioffd + i - 1 ] = r.getValue();
      e[ ioffe + i - 1 ] = sn.getValue() * d[ ioffd + i ];
      d[ ioffd + i ] = cs.getValue() * d[ ioffd + i ];
      if ( rotate ) {
        work[ ioffwork + i - 1 ] = cs.getValue();
        work[ ioffwork + n + i - 1 ] = sn.getValue();
      }
    } // 20
    if ( sqre1 == 1 ) {
      LaPack1.dlartg( d[ ioffd + n - 1 ], e[ ioffe + n - 1 ],
        cs, sn, r );
      d[ ioffd + n - 1 ] = r.getValue();
      if ( rotate ) {
        work[ ioffwork + n - 1 ] = cs.getValue();
        work[ ioffwork + n + n - 1 ] = sn.getValue();
      }
    }
    if ( nru > 0 ) {
      if ( sqre1 == 0 ) {
        LaPack0.dlasr( 'R', 'V', 'F', nru, n, work, work, U, ldu,
          ioffwork, ioffwork + np1 - 1, ioffu );
      } else {
        LaPack0.dlasr( 'R', 'V', 'F', nru, np1, work, work, U, ldu,
          ioffwork, ioffwork + np1 - 1, ioffu );
      }
    }
    if ( ncc > 0 ) {
      if ( sqre1 == 0 ) {
        LaPack0.dlasr( 'L', 'V', 'F', n, ncc, work, work, C, ldc,
          ioffwork, ioffwork + np1 - 1, ioffc );
      } else {
        LaPack0.dlasr( 'L', 'V', 'F', np1, ncc, work, work, C, ldc,
          ioffwork, ioffwork + np1 - 1, ioffc );
      }
    }
  }
  LaPack5.dbdsqr( 'U', n, ncvt, nru, ncc, d, e, Vt, ldvt, U, ldu,
    C, ldc, work, info, ioffd, ioffe, ioffvt, ioffu, ioffc,
    ioffwork );
  for ( i = 1; i <= n; i ++ ) {
    var isub = i;
    var smin = d[ ioffd + i - 1 ];
    for ( var j = i + 1; j <= n; j ++ ) {
      if ( d[ ioffd + j - 1 ] < smin ) {
        isub = j;
        smin = d[ ioffd + j - 1 ];
      }
    } // 30
    if ( isub != i ) {
      d[ ioffd + isub - 1 ] = d[ ioffd + i - 1 ];
      d[ ioffd + i - 1 ] = smin;
      if ( ncvt > 0 ) {
        Blas1.dswap( ncvt, Vt, ldvt, Vt, ldvt, ioffvt + isub - 1,
          ioffvt + i - 1 );
      }
      if ( nru > 0 ) {
        Blas1.dswap( nru, U, 1, U, 1, ioffu + ( isub - 1 ) * ldu,
          ioffu + ( i - 1 ) * ldu );
      }
      if ( ncc > 0 ) {
        Blas1.dswap( ncc, C, ldc, C, ldc, ioffc + isub - 1,
          ioffc + i - 1 );
      }
    }
  } // 40
}
//**************************************************************************
LaPack6.dtgexc = function( wantq, wantz, n, A, lda, B, ldb, Q,
ldq, Z, ldz, ifst, ilst, work, lwork, info, ioffa, ioffb, ioffq, ioffz,
ioffwork ) {
  throw new Error("not programmed: generalized eigenvalue");
}
