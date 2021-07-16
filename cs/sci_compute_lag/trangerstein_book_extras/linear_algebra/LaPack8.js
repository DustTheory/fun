function LaPack8() {
}
//**************************************************************************
LaPack8.dbdsdc = function( uplo, compq, n, d, e, U, ldu, Vt, ldvt,
q, iq, work, iwork, info, ioffd, ioffe, ioffu, ioffvt, ioffq, ioffiq,
ioffwork, ioffiwork ) {
  throw new Error("not tested");
  info.setValue( 0 );
  var iuplo = 0;
  if ( uplo.charAt(0).toUpperCase() == 'U' ) iuplo = 1;
  if ( uplo.charAt(0).toUpperCase() == 'L' ) iuplo = 2;
  if ( compq.charAt(0).toUpperCase() == 'N' ) var icompq = 0;
  else if ( compq.charAt(0).toUpperCase() == 'P' ) icompq = 1;
  else if ( compq.charAt(0).toUpperCase() == 'I' ) icompq = 2;
  else icompq = -1;
  if ( iuplo == 0 ) info.setValue( -1 );
  else if ( icompq < 0 ) info.setValue( -2 );
  else if ( n < 0 ) info.setValue( -3 );
  else if ( ldu < 1 || ( icompq == 2 && ldu < n ) ) info.setValue( -7 );
  else if ( ldvt < 1 || ( icompq == 2 && ldvt < n ) ) info.setValue( -9 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dbdsdc', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  var smlsiz = LaPack0.ilaenv( 9, 'dbdsdc', 0, 0, 0, 0 );
  if ( n == 1 ) {
    if ( icompq == 1 ) {
      q[ ioffq ] = ( d[ ioffd ] >= 0. ? 1. : -1. );
      q[ ioffq + smlsiz * n ] = 1.;
    } else if ( icompq == 2 ) {
      U[ ioffu ] = ( d[ ioffd ] >= 0. ? 1. : -1. );
      Vt[ ioffvt ] = 1.;
    }
    d[ ioffd ] = Math.abs( d[ ioffd ] );
    return;
  }
  var nm1 = n - 1;
  var wstart = 1;
  var qstart = 3;
  if ( icompq == 1 ) {
    Blas1.dcopy( n, d, 1, q, 1, ioffd, ioffq );
    Blas1.dcopy( n - 1, e, 1, q, 1, ioffe, ioffq + n );
  }
  var csReference = new NumberReference();
  var snReference = new NumberReference();
  var rReference = new NumberReference();
  if ( iuplo == 2 ) {
    qstart = 5;
    wstart = 2 * n - 1;
    for ( var i = 1; i <= n - 1; i ++ ) {
      LaPack1.dlartg( d[ ioffd + i - 1], e[ ioffe + i - 1], cs, sn, r );
      d[ ioffd + i - 1 ] = r.getValue();
      e[ ioffe + i- 1 ] = sn.getValue() * d[ ioffd + i ];
      d[ ioffd + i ] = cs.getValue() * d[ ioffd + i ];
      if ( icompq == 1 ) {
        q[ ioffq + i + 2 * n - 1 ] = cs.getValue();
        q[ ioffq + i + 3 * n - 1 ] = sn.getValue();
      } else if ( icompq == 2 ) {
        work[ ioffwork + i - 1 ] = cs.getValue();
        work[ ioffwork + nm1 + i - 1 ] = - sn.getValue();
      }
    } // 10
  }
  var goto40 = false;
  if ( icompq == 0 ) {
    LaPack6.dlasdq( 'U', 0, n, 0, 0, 0, d, e, Vt, ldvt, U, ldu, U, ldu,
      work, info, ioffd, ioffe, ioffvt, ioffu, ioffu,
      ioffwork + wstart - 1 );
    goto40 = true;
  }
  if ( ! goto40 ) {
    if ( n <= smlsiz ) {
      if ( icompq == 2 ) {
        LaPack0.dlaset( 'A', n, n, 0., 1., U, ldu, ioffu );
        LaPack0.dlaset( 'A', n, n, 0., 1., Vt, ldvt, ioffvt );
        LaPack6.dlasdq( 'U', 0, n, n, n, 0, d, e, Vt, ldvt, U, ldu,
          U, ldu, work, info, ioffd, ioffe, ioffvt, ioffu, ioffu,
          ioffwork + wstart - 1 );
      } else if ( icompq == 1 ) {
        var iu = 1;
        var ivt = iu + n;
        LaPack0.dlaset( 'A', n, n, 0., 1., q, n,
          ioffq + iu + ( qstart - 1 ) * n );
        LaPack0.dlaset( 'A', n, n, 0., 1., q, n,
          ioffq + ivt + ( qstart - 1 ) * n );
        LaPack6.dlasdq( 'U', 0, n, n, n, 0, d, e, q, n, q, n,
          q, n, work, info, ioffd, ioffe,
          ioffq + ivt - 1 + ( qstart - 1 ) * n, 
          ioffq + iu - 1 + ( qstart - 1 ) * n, 
          ioffq + iu - 1 + ( qstart - 1 ) * n, ioffwork + wstart - 1 );
      }
      goto40 = true;
    }
  }
  var ierr = new IntReference();
  if ( ! goto40 ) {
    if ( icompq == 2 ) {
      LaPack0.dlaset( 'A', n, n, 0., 1., U, ldu, ioffu );
      LaPack0.dlaset( 'A', n, n, 0., 1., Vt, ldvt, ioffvt );
    }
    var orgnrm = LaPack1.dlanst( 'M', n, d, e, ioffd, ioffe );
    if ( orgnrm == 0. ) return;
    LaPack1.dlascl( 'G', 0, 0, orgnrm, 1., n, 1, d, n, ierr, ioffd );
    LaPack1.dlascl( 'G', 0, 0, orgnrm, 1., nm1, 1, e, nm1, ierr,
      ioffe );
    var eps = 0.9 * LaPack0.dlamch( 'Epsilon' );
    var mlvl = Math.round( Math.log( Number( n ) / Number( smlsiz + 1 ) )
      / Math.log( 2. ) ) + 1;
    var smlszp = smlsiz + 1;
    if ( icompq == 1 ) {
      iu = 1;
      ivt = 1 * smlsiz;
      var difl = ivt + smlszp;
      var difr = difl + mlvl;
      var z = difr + mlvl * 2;
      var ic = z + mlvl;
      var is2 = ic + 1;
      var poles = is2 + 1;
      var givnum = poles + 2 * mlvl;
      var k = 1;
      var givptr = 2;
      var perm = 3;
      var givcol = perm + mlvl;
    }
    for ( i = 1; i <= n; i ++ ) {
      if ( Math.abs( d[ ioffd + i - 1 ] ) < eps ) {
        d[ ioffd + i - 1 ] =
          ( d[ ioffd + i - 1 ] >= 0. ? eps : -eps );
      }
    } // 20
    var start = 1;
    var sqre = 0;
    for ( i = 1; i <= nm1; i ++ ) {
      if ( Math.abs( e[ ioffe + i - 1 ] ) < eps || i == nm1 ) {
        if ( i < nm1 ) var nsize = i - start + 1;
        else if ( Math.abs( e[ ioffe + i - 1 ] ) >= eps ) {
          nsize = n - start + 1;
        } else {
          nsize = i - start + 1;
          if ( icompq == 2 ) {
            U[ ioffu + n - 1 + ( n - 1 ) * ldu ] =
              ( d[ ioffd + n - 1 ] >= 0. ? 1. : -1. );
            Vt[ ioffvt + n - 1 + ( n - 1 ) * ldvt ] = 1.;
          } else if ( icompq == 1 ) {
            q[ ioffq + n - 1 + ( qstart - 1 ) * n ] =
              ( d[ ioffd + n - 1 ] >= 0. ? 1. : -1. );
            q[ ioffq + n - 1 + ( smlsiz + qstart - 1 ) * n ] = 1.;
          }
          d[ ioffd + n - 1 ] = Math.abs( d[ ioffd + n - 1 ] );
        }
        if ( icompq == 2 ) {
          LaPack7.dlasd0( nsize, sqre, d, e, U, ldu, Vt, ldvt, smlsiz,
            iwork, work, info, ioffd + start - 1, ioffe + start - 1,
            ioffu + start - 1 + ( start - 1 ) * ldu,
            ioffvt + start - 1 + ( start - 1 ) * ldvt,
            ioffiwork, ioffwork + wstart - 1 );
        } else {
          LaPack7.dlasda( icompq, smlsiz, nsize, sqre, d, e, q, n,
            q, iq, q, q, q, q, iq, iq, n, iq, q, q, q, work, iwork,
            info, ioffd + start - 1, ioffe + start - 1,
            ioffq + start - 1 + ( iu + qstart - 2 ) * n,
            ioffq + start - 1 + ( ivt + qstart - 2 ) * n,
            ioffiq + start - 1 + k * n,
            ioffq + start - 1 + ( difl + qstart - 2 ) * n,
            ioffq + start - 1 + ( difr + qstart - 2 ) * n,
            ioffq + start - 1 + ( z + qstart - 2 ) * n,
            ioffq + start - 1 + ( poles + qstart - 2 ) * n,
            ioffiq + start - 1 + givptr * n,
            ioffiq + start - 1 + givcol * n,
            ioffiq + start - 1 + perm * n,
            ioffq + start - 1 + ( givnum + qstart - 2 ) * n,
            ioffq + start - 1 + ( ic + start - 2 ) * n,
            ioffq + start - 1 + ( is2 + qstart - 2 ) * n,
            ioffwork + wstart - 1, ioffiwork );
        }
        if ( info.getValue() != 0 ) return;
        start = i + 1;
      }
    } // 30
    LaPack1.dlascl( 'G', 0, 0, 1., orgnrm, n, 1, d, n, ierr, ioffd );
  } // 40
  for ( var ii = 2; ii <= n; ii ++ ) {
    i = ii - 1;
    var kk = i;
    var p = d[ ioffd + i - 1 ];
    for ( var j = ii; j <= n; j ++ ) {
      if ( d[ ioffd + j - 1 ] > p ) {
        kk = j;
        p = d[ ioffd + j - 1 ];
      }
    } // 50
    if ( kk != i ) {
      d[ ioffd + kk - 1 ] = d[ ioffd + i - 1 ];
      d[ ioffd + i - 1 ] = p;
      if ( icompq == 1 ) iq[ ioffiq + i - 1 ] = kk;
      else if ( icompq == 2 ) {
        Blas1.dswap( n, U, 1, U, 1, ioffu + ( i - 1 ) * ldu,
          ioffu + ( kk - 1 ) * ldu );
        Blas1.dswap( n, Vt, ldvt, Vt, ldvt, ioffvt + i - 1,
          ioffvt + kk - 1 );
      }
    } else if ( icompq == 1 ) iq[ ioffiq + i - 1 ] = i;
  } // 60
  if ( icompq == 1 ) iq[ ioffiq + n - 1 ] = ( iuplo == 1 ? 1 : 0 );
  if ( iuplo == 2 && icompq == 2 ) {
    LaPack0.dlasr( 'L', 'V', 'B', n, n, work, work, U, ldu,
      ioffwork, ioffwork + n - 1, ioffu );
  }
}
//**************************************************************************
LaPack8.dlalsd = function( uplo, smlsiz, n, nrhs, d, e, B, ldb,
rcond, rank, work, iwork, info, ioffd, ioffe, ioffb, ioffwork, ioffiwork ) {
  throw new Error("not tested");
  info.setValue( 0 );
  if ( n < 0 ) info.setValue( -3 );
  else if ( nrhs < 1 ) info.setValue( -4 );
  else if ( ldb < 1 || ldt > n ) info.setValue( -8 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlalsd', - info.getValue() );
    return;
  }
  var eps = LaPack0.dlamch( 'Epsilon' );
  var rcnd = ( rcond <= 0. || rcond >= 1. ? eps : rcond );
  rank.setValue( 0 );
  if ( n == 0 ) return;
  else if ( n == 1 ) {
    if ( d[ ioffd ] == 0. ) {
      LaPack0.dlaset( 'A', 1, nrhs, 0., 0., B, ldb, ioffb );
    } else {
      rank.setValue( 1 );
      LaPack1.dlascl( 'G', 0, 0, d[ ioffd ], 1., 1, nrhs, B, ldb, info,
        ioffb );
      d[ ioffd ] = Math.abs( d[ ioffd ] );
    }
    return;
  }
  var csReference = new NumberReference();
  var snReference = new NumberReference();
  var rReference = new NumberReference();
  if ( uplo.charAt(0).toUpperCase() == 'L' ) {
    for ( var i = 1; i <= n - 1; i ++ ) {
      LaPack1.dlartg( d[ ioffd + i - 1 ], e[ ioffe + i - 1 ],
        cs, sn, r );
      d[ ioffd + i - 1 ] = r.getValue();
      e[ ioffe + i - 1 ] = sn.getValue() * d[ ioffd + i ];
      d[ ioffd + i ] *= cs.getValue();
      if ( nrhs == 1 ) {
        Blas1.drot( 1, B, 1, B, 1, cs.getValue(), sn.getValue(),
          ioffb + i - 1, ioffb + i );
      } else {
        work[ ioffwork + i * 2 - 2 ] = cs.getValue();
        work[ ioffwork + i * 2 - 1 ] = sn.getValue();
      }
    } // 10
    if ( nrhs > 1 ) {
      for ( i = 1; i <= nrhs; i ++ ) {
        for ( var j = 1; j <= n - 1; j ++ ) {
          cs.setValue( work[ ioffwork + j * 2 - 2 ] );
          sn.setValue( work[ ioffwork + j * 2 - 1 ] );
          Blas1.drot( 1, B, 1, B, 1, cs.getValue(), sn.getValue(),
            ioffb + j - 1, ioffb + j );
        } // 20
      } // 30
    }
  }
  var nm1 = n - 1;
  var orgnrm = LaPack1.dlanst( 'M', n, d, e, ioffd, ioffe );
  if ( orgnrm == 0. ) {
    LaPack0.dlaset( 'A', n, nrhs, 0., 0., B, ldb, ioffb );
    return;
  }
  LaPack1.dlascl( 'G', 0, 0, orgnrm, 1., n, 1, d, n, info, ioffd );
  LaPack1.dlascl( 'G', 0, 0, orgnrm, 1., nm1, 1, e, nm1, info,
    ioffe );
  if ( n <= smlsiz ) {
    nwork = 1 + n * n;
    LaPack0.dlaset( 'A', n, n, 0., 1., work, n, ioffwork );
    LaPack6.dlasdq( 'U', 0, n, n, 0, nrhs, d, e, work, n, work, n,
      B, ldb, work, info, ioffd, ioffe, ioffwork, ioffwork, ioffb,
      ioffwork + nwork - 1 );
    if ( info.getValue() != 0 ) return;
    tol = rcnd
      * Math.abs( d[ ioffd + Blas1.idamax( n, d, 1, ioffd ) - 1 ] );
    for ( i = 1; i <= n; i ++ ) {
      if ( d[ ioffd + i - 1 ] <= tol ) {
        LaPack0.dlaset( 'A', 1, nrhs, 0., 0., B, ldb, ioffb + i - 1 );
      } else {
        LaPack1.dlascl( 'G', 0, 0, d[ ioffd + i - 1 ], 1., 1, nrhs,
          B, ldb, info, ioffb + i - 1 );
        rank.setValue( rank.getValue() + 1 );
      }
    } // 40
    Blas3.dgemm( 'T', 'N', n, nrhs, n, 1., work, n, B, ldb, 0., work,
      n, ioffwork, ioffb, ioffwork + nwork - 1 );
    LaPack0.dlacpy( 'A', n, nrhs, work, n, B, ldb,
      ioffwork + nwork - 1, ioffb );
    LaPack1.dlascl( 'G', 0, 0, 1., orgnrm, n, 1, d, n, info, ioffd );
    LaPack0.dlasrt( 'D', n, d, info, ioffd );
    LaPack1.dlascl( 'G', 0, 0, orgnrm, 1., n, nrhs, B, ldb, info,
      ioffb );
    return;
  }
  var nlvl = Math.round( Math.log( Number( n ) / Number( smlsiz + 1 ) )
    / Math.log( 2. ) ) + 1;
  var smlszp = smlsiz + 1;
  var u = 1;
  var vt = 1 + smlsiz * n;
  var difl = vt + smlszp * n;
  var difr = difl + nlvl * n;
  var z = difr + nlvl * n * 2;
  var c = z + nlvl * n;
  var s = c + n;
  var poles = s + n;
  var givnum = poles + 2 * nlvl * n;
  var bx = givnum + 2 * nlvl * n;
  var nwork = bx + n * nrhs;
  var sizei = 1 + n;
  var k = sizei + n;
  var givptr = k + n;
  var perm = givptr + n;
  var givcol = perm + nlvl * n;
  var iwk = givcol + nlvl * n * 2;
  var st = 1;
  var sqre = 0;
  var icmpq1 = 1;
  var icmpq2 = 0;
  var nsub = 0;
  for ( i = 1; i <= n; i ++ ) {
    if ( Math.abs( d[ ioffd + i - 1 ] ) < eps ) {
      d[ ioffd + i - 1 ] = ( d[ ioffd + i - 1 ] >= 0. ? eps : - eps );
    }
  } // 50
  for ( i = 1; i <= nm1; i ++ ) {
    if ( Math.abs( e[ ioffe + i - 1 ] ) < eps || i == nm1 ) {
      nsub ++;
      iwork[ ioffiwork + nsub - 1 ] = st;
      if ( i < nm1 ) {
        var nsize = i - st + 1;
        iwork[ ioffiwork + sizei + nsub - 2 ] = nsize;
      } else if ( Math.abs( e[ ioffe + i - 1 ] ) >= eps ) {
        nsize = n - st + 1;
        iwork[ ioffiwork + sizei + nsub - 2 ] = nsize;
      } else {
        nsize = i - st + 1;
        iwork[ ioffiwork + sizei + nsub - 2 ] = nsize;
        nsub ++;
        iwork[ ioffiwork + nsub - 1 ] = n;
        iwork[ ioffiwork + sizei + nsub - 2 ] = 1;
        Blas1.dcopy( nrhs, B, ldb, work, n, ioffb + n - 1,
          ioffwork + bx + nm1 - 1 );
      }
      var st1 = st - 1;
      if ( nsize == 1 ) {
        Blas1.dcopy( nrhs, B, ldb, work, n, ioffb + st - 1,
          ioffwork + bx + st1 -1 );
      } else if ( nsize <= smlsiz ) {
        LaPack0.dlaset( 'A', nsize, nsize, 0., 1., work, n,
          ioffwork + vt + st1 - 1 );
        LaPack6.dlasdq( 'U', 0, nsize, nsize, 0, nrhs, d, e, work, n,
          work, n, B, ldb, work, info, ioffd + st - 1, ioffe + st - 1,
          ioffwork + vt + st1 - 1, ioffwork + nwork - 1, ioffb + st - 1,
          ioffwork + nwork - 1 );
        if ( info.getValue() != 0 ) return;
        LaPack0.dlacpy( 'A', nsize, nrhs, B, ldb, work, n,
          ioffb + st - 1, ioffwork + bx + st1 - 1 );
      } else {
        LaPack7.dlasda( icmpq1, smlsiz, nsize, sqre, d,
          e, work, n, work,
          iwork, work,
          work, work,
          work, iwork,
          iwork, n, iwork,
          work, work,
          work, work, iwork,
          info,
          ioffd + st - 1,
          ioffe + st - 1, ioffwork + u + st1 - 1,
            ioffwork + vt + st1 - 1,
          ioffiwork + k + st1 - 1, ioffwork + difl + st1 - 1,
          ioffwork + difr + st1 - 1, ioffwork + z + st1 - 1,
          ioffwork + poles + st1 - 1, ioffiwork + givptr + st1 - 1,
          ioffiwork + givcol + st1 - 1, ioffiwork + perm + st1 - 1,
          ioffwork + givnum + st1 - 1, ioffwork + c + st1 - 1,
          ioffwork + s + st1 - 1, ioffwork + nwork - 1,
          ioffiwork + iwk - 1 );
        if ( info.getValue() != 0 ) return;
        var bxst = bx + st1;
        LaPack3.dlalsa( icmpq2, smlsiz, nsize, nrhs, B,
          ldb, work, n, work, n,
          work, iwork,
          work, work,
          work, work,
          iwork, iwork, n,
          iwork, work,
          work, work, work,
          iwork, info,
          ioffb + st - 1,
          ioffwork + bxst - 1, ioffwork + u + st1 - 1,
          ioffwork + vt + st1 - 1, ioffiwork + k + st1 - 1,
          ioffwork + difl + st1 - 1, ioffwork + difr + st1 - 1,
          ioffwork + z + st1 - 1, ioffwork + poles + st1 - 1,
          ioffiwork + givptr + st1 - 1, ioffiwork + givcol + st1 - 1,
          ioffiwork + perm + st1 - 1, ioffwork + givnum + st1 - 1,
          ioffwork + c + st1 - 1, ioffwork + s + st1 - 1,
            ioffwork + nwork - 1,
          ioffiwork + iwk - 1 );
        if ( info.getValue() != 0 ) return;
      }
      st = i + 1;
    }
  } // 60
  var tol =
    rcnd * Math.abs( d[ ioffd + Blas1.idamax( n, d, 1, ioffd ) - 1 ] );
  for ( i = 1; i <= n; i ++ ) {
    if ( Math.abs( d[ ioffd + i - 1 ] ) <= tol ) {
      LaPack0.dlaset( 'A', 1, nrhs, 0., 0., work, n,
        ioffwork + bx + i - 2 );
    } else {
      rank.setValue( rank.getValue() + 1 );
      LaPack1.dlascl( 'G', 0, 0, d[ ioffd + i - 1 ], 1., 1, nrhs, work,
        n, info, ioffwork + bx + i - 2 );
    }
    d[ ioffd + i - 1 ] = Math.abs( d[ ioffd + i - 1 ] );
  } // 70
  icmpq2 = 1;
  for ( i = 1; i <= nsub; i ++ ) {
    st = iwork[ ioffiwork + i - 1 ];
    st1 = st - 1;
    nsize = iwork[ ioffiwork + sizei + i - 2 ];
    bxst = bx + st1;
    if ( nsize == 1 ) {
      Blas1.dcopy( nrhs, work, n, B, ldb, ioffwork + bxst - 1,
        ioffb + st - 1 );
    } else if ( nsize <= smlsiz ) {
      Blas3.dgemm( 'T', 'N', nsize, nrhs, nsize, 1., work, n, work, n,
        0., B, ldb, ioffwork + vt + st - 1, ioffwork + bxst - 1,
        ioffb + st - 1 );
    } else {
      LaPack3.dlalsa( icmpq2, smlsiz, nsize, nrhs, work, n,
        B, ldb, work, n,
        work, iwork,
        work, work,
        work, work,
        iwork, iwork, n,
        iwork, work,
        work, work, work,
        iwork, info,
        ioffwork + bxst - 1,
        ioffb + st - 1, ioffwork + u + st1 - 1,
        ioffwork + vt + st1 - 1, ioffiwork + k + st1 - 1,
        ioffwork + difl + st1 - 1, ioffwork + difr + st1 - 1,
        ioffwork + z + st1 - 1, ioffwork + poles + st1 - 1,
        ioffiwork + givptr + st1 - 1, ioffiwork + givcol + st1 - 1,
        ioffiwork + perm + st1 - 1, ioffwork + givnum + st1 - 1,
        ioffwork + c + st1 - 1, ioffwork + s + st1 - 1,
          ioffwork + nwork - 1,
        ioffiwork + iwk - 1 );
      if ( info.getValue() != 0 ) return;
    }
  } // 80
  LaPack1.dlascl( 'G', 0, 0, 1., orgnrm, n, 1, d, n, info, ioffd );
  LaPack0.dlasrt( 'D', n, d, info, ioffd );
  LaPack1.dlascl( 'G', 0, 0, orgnrm, 1., n, nrhs, B, ldb, info,
    ioffb );
}
//**************************************************************************
LaPack8.dlaqr0 = function( wantt, wantz, n, ilo, ihi, H, ldh, wr,
wi, iloz, ihiz, Z, ldz, work, lwork, info, ioffh, ioffwr, ioffwi, ioffz,
ioffwork ) {
  throw new Error("not tested");
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
      LaPack0.ilaenv( 13, 'dlaqr0', jbcmpz, n, ilo, ihi, lwork );
    nwr = Math.max( 2, nwr );
    nwr = Math.min( Math.min( ihi - ilo + 1, ( n - 1 ) / 3 ), nwr );
    var nsr =
      LaPack0.ilaenv( 15, 'dlaqr0', jbcmpz, n, ilo, ihi, lwork );
    nsr = Math.min( nsr, Math.min( ( n + 6 ) / 9, ihi - ilo ) );
    nsr = Math.max( 2, nsr - nsr % 2 );
    var ls = new IntReference();
    var ld = new IntReference();
    LaPack7.dlaqr3( wantt, wantz, n, ilo, ihi, nwr + 1, H, ldh, iloz,
      ihiz, Z, ldz, ls, ld, wr, wi, H, ldh, n, H, ldh, n, H, ldh, work,
      -1, ioffh, ioffz, ioffwr, ioffwi, ioffh, ioffh, ioffh );
    lwkopt = Math.max( 3 * nsr / 2, Math.round( work[ ioffwork ] ) );
    if ( lwork == -1 ) {
      work[ ioffwork ] = Number( lwkopt );
      return;
    }
    var nmin =
      LaPack0.ilaenv( 12, 'dlaqr0', jbcmpz, n, ilo, ihi, lwork );
    nmin = Math.max( ntiny, nmin );
    var nibble =
      LaPack0.ilaenv( 14, 'dlaqr0', jbcmpz, n, ilo, ihi, lwork );
    nibble = Math.max( 0, nibble );
    var kacc22 =
      LaPack0.ilaenv( 16, 'dlaqr0', jbcmpz, n, ilo, ihi, lwork );
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
    var aaReference = new NumberReference();
    var bbReference = new NumberReference();
    var ccReference = new NumberReference();
    var ddReference = new NumberReference();
    var rt1rReference = new NumberReference();
    var rt1iReference = new NumberReference();
    var rt2rReference = new NumberReference();
    var rt2iReference = new NumberReference();
    var csReference = new NumberReference();
    var snReference = new NumberReference();
    for ( var it = 1; it <= itmax; it ++ ) {
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
      var nw = ( ndfl < kexnw ? Math.min( nwupbd, nwr ) :
        Math.min( nwupbd, 2 * nw ) );
      if ( nw < nwmax ) {
        if ( nw >= nh - 1 ) nw = nh;
        else {
          kwtop = kbot - nw + 1;
          if ( Math.abs( H[ ioffh + kwtop - 1 + ( kwtop - 2 ) * ldh ] )
          > Math.abs( H[ ioffh + kwtop - 2 + ( kwtop - 3 ) * ldh ] ) ) {
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
      LaPack7.dlaqr3( wantt, wantz, n, ktop, kbot, nw, H, ldh, iloz,
        ihiz, Z, ldz, ls, ld, wr, wi, H, ldh, nho, H, ldh, nve, H, ldh,
        work, lwork, ioffh, ioffz, ioffwr, ioffwi, ioffh + kv - 1,
        ioffh + kv - 1 + ( kt - 1 ) * ldh, ioffh + kwv - 1, ioffwork );
      kbot -= ld.getValue();
      ks = kbot - ls.getValue() + 1;
      if ( ld.getValue() == 0 || ( 100 * ld.getValue() <= nw * nibble &&
      kbot - ktop + 1 > Math.min( nmin, nwmax ) ) ) {
        var ns = Math.min( Math.min( nsmax, nsr ),
          Math.max( 2, kbot - ktop ) );
        ns -= ns % 2;
        if ( ndfl % kexsh == 0 ) {
          ks = kbot - ns + 1;
          for ( var i = kbot; i >= Math.max( ks + 1, ktop + 2 );
          i -= 2 ) {
            var ss =
              Math.abs( H[ ioffh + i - 1 + ( i - 2 ) * ldh ] )
              + Math.abs( H[ ioffh + i - 2 + ( i - 3 ) * ldh ] );
            aa.setValue(
              wilk1 * ss + H[ ioffh + i - 1 + ( i - 1 ) * ldh ] );
            bb.setValue( ss );
            cc.setValue( wilk2 * ss );
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
            if ( ns > nmin ) {
              LaPack6.dlaqr4( false, false, ns, 1, ns, H, ldh, wr, wi,
                1, 1, zdum, 1, work, lwork, inf, ioffh + kt - 1,
                ioffwr + ks - 1, ioffwi + ks - 1, 0, ioffwork );
            } else {
              LaPack2.dlahqr( false, false, ns, 1, ns, H, ldh, wr, wi,
                1, 1, zdum, 1, inf, ioffh + kt - 1, ioffwr + ks - 1,
                ioffwi + ks - 1, 0 );
            }
            ks += inf.getValue();
            if ( ks >= kbot ) {
              aa.setValue( H[ ioffh + kbot - 2 + ( kbot - 2 ) * ldh ] );
              cc.setValue( H[ ioffh + kbot - 1 + ( kbot - 2 ) * ldh ] );
              bb.setValue( H[ ioffh + kbot - 2 + ( kbot - 1 ) * ldh ] );
              dd.setValue( H[ ioffh + kbot - 1 + ( kbot - 1 ) * ldh ] );
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
              swap = wr[ ioffwr + i - 1 ]
              wr[ ioffwr + i - 1 ] = wr[ ioffwr + i - 2 ];
              wr[ ioffwr + i - 2 ] = wr[ ioffwr + i - 3 ];
              wr[ ioffwr + i - 3 ] = swap;
              swap = wi[ ioffwi + i - 1 ]
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
        kdu = 3 * ns - 3;
        ku = n - ku + 1;
        kwh = kdu + 1;
        nho = ( n - kdu + 1 - 4 ) - ( kdu + 1 ) + 1;
        kwv = kdu + 4;
        nve = n - ndu - kwv + 1;
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
