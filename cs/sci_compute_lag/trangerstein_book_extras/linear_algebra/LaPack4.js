function LaPack4() {
}
//*************************************************************************
LaPack4.dgeqp3 = function( m, n, A, lda, jpvt, tau, work, lwork,
info, ioffa, ioffjpvt, iofftau, ioffwork ) {
  var inb = 1;
  var inbmin = 2;
  var ixover = 3;
  info.setValue( 0 );
  var lquery = ( lwork == -1 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -4 );
  if ( info.getValue() == 0 ) {
    var minmn = Math.min( m, n );
    if ( minmn == 0 ) {
      var iws = 1;
      var lwkopt = 1;
    } else {
      iws = 3 * n + 1;
      var nb = LaPack0.ilaenv( inb, 'dgeqrf', ' ', m, n, -1, -1 );
      lwkopt = 2 * n + ( n + 1 ) * nb;
    }
    work[ ioffwork ] = lwkopt;
    if ( lwork < iws && ! lquery ) info.setValue( -8 );
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgeqp3', - info.getValue() );
    return;
  } else if ( lquery ) return;
  if ( minmn == 0 ) return;
  var nfxd = 1;
  for ( var j = 1; j <= n; j ++ ) {
    if ( jpvt[ ioffjpvt + j - 1 ] != 0 ) {
      if ( j != nfxd ) {
        Blas1.dswap( m, A, 1, A, 1, ioffa + ( j - 1 ) * lda,
          ioffa + ( nfxd - 1 ) * lda );
        jpvt[ ioffjpvt + j - 1 ] = jpvt[ ioffjpvt + nfxd - 1 ];
        jpvt[ ioffjpvt + nfxd - 1 ] = j;
      } else jpvt[ ioffjpvt + j - 1 ] = j;
      nfxd ++;
    } else jpvt[ ioffjpvt + j - 1 ] = j;
  }
  nfxd --;
  if ( nfxd > 0 ) {
    var na = Math.min( m, nfxd );
    LaPack3.dgeqrf( m, na, A, lda, tau, work, lwork, info,
      ioffa, iofftau, ioffwork );
    iws = Math.max( iws, Math.round( work[ ioffwork ] ) );
    if ( na < n ) {
      LaPack3.dormqr( 'Left', 'Transpose', m, n - na, na, A, lda, tau,
        A, lda, work, lwork, info, ioffa, iofftau, ioffa + na * lda,
        ioffwork );
      iws = Math.max( iws, Math.round( work[ ioffwork ] ) );
    }
  }
  if ( nfxd < minmn ) {
    var sm = m - nfxd;
    var sn = n - nfxd;
    var sminmn = minmn - nfxd;
    nb = LaPack0.ilaenv( inb, 'dgeqrf', ' ', sm, sn, -1, -1 );
    var nbmin = 2;
    var nx = 0;
    if ( nb > 1 && nb < sminmn ) {
      nx = Math.max( 0,
        LaPack0.ilaenv( ixover, 'dgeqrf', ' ', sm, sn, -1, -1 ) );
      if ( nx < sminmn ) {
        var minws = 2 * sn + ( sn + 1 ) * nb;
        iws = Math.max( iws, minws );
        if ( lwork < minws ) {
          nb = ( lwork - 2 * sn ) / ( sn + 1 );
          nbmin = Math.max( 2,
            LaPack0.ilaenv( inbmin, 'dgeqrf', ' ', sm, sn, -1, -1 ) );
        }
      }
    }
    for ( j = nfxd + 1; j <= n; j ++ ) {
      work[ ioffwork + j - 1 ] =
        Blas1.dnrm2( sm, A, 1, ioffa + nfxd + ( j - 1 ) * lda );
      work[ ioffwork + n + j - 1 ] = work[ ioffwork + j - 1 ];
    }
    if ( nb >= nbmin && nb < sminmn && nx < sminmn ) {
      j = nfxd + 1;
      var topbmn = minmn - nx;
      while ( j <= topbmn ) {
        var jb = Math.min( nb, topbmn - j + 1 );
        var fjb = new IntReference();
        LaPack2.dlaqps( m, n - j + 1, j - 1, jb, fjb, A, lda, jpvt,
          tau, work, work, work, work, n - j + 1,
          ioffa + ( j - 1 ) * lda, ioffjpvt + j - 1, iofftau + j - 1,
          ioffwork + j - 1, ioffwork + n + j - 1, ioffwork + 2 * n,
          ioffwork + 2 * n + jb );
        j += fjb.getValue();
      }
    } else j = nfxd + 1;
    if ( j <= minmn ) {
      LaPack2.dlaqp2( m, n - j + 1, j - 1, A, lda, jpvt, tau, work,
        work, work, ioffa + ( j - 1 ) * lda, ioffjpvt + j - 1,
        iofftau + j - 1, ioffwork + j - 1, ioffwork + n + j - 1,
        ioffwork + 2 * n );
    }
  }
  work[ ioffwork ] = iws;
}
//*************************************************************************
LaPack4.dggqrf = function( n, m, p, A, lda, taua, B, ldb, taub,
work, lwork, info, ioffa, iofftaua, ioffb, iofftaub, ioffwork ) {
  throw new Error("not programmed: generalized factorization");
}
LaPack4.zggqrf = function( n, m, P, A, lda, taua, B, ldb, taub,
work, lwork, info ) {
  throw new Error("not programmed: complex generalized factorization");
}
//*************************************************************************
LaPack4.dggrqf = function( m, p, n, A, lda, taua, B, ldb, taub,
work, lwork, info, ioffa, iofftaua, ioffb, iofftaub, ioffwork ) {
  throw new Error("not programmed: generalized factorization");
}
LaPack4.zggrqf = function( m, p, n, A, lda, taua, B, ldb, taub,
work, lwork, info ) {
  throw new Error("not programmed: complex generalized factorization");
}
//*************************************************************************
LaPack4.dggsvp = function( jobu, jobv, jobq, m, p, n, A, lda, B,
ldb, tola, tolb, k, l, U, ldu, V, ldv, Q, ldq, iwork, tau, work, info,
ioffa, ioffb, ioffu, ioffv, ioffq, ioffiwork, iofftau, ioffwork ) {
  var wantu = ( jobu.charAt(0).toUpperCase() == 'U' );
  var wantv = ( jobv.charAt(0).toUpperCase() == 'V' );
  var wantq = ( jobq.charAt(0).toUpperCase() == 'Q' );
  var forwrd = true;
  info.setValue( 0 );
  if ( ! wantu && jobu.charAt(0).toUpperCase() != 'N' ) {
    info.setValue( -1 );
  } else if ( ! wantv && jobv.charAt(0).toUpperCase() != 'N' ) {
    info.setValue( -2 );
  } else if ( ! wantq && jobq.charAt(0).toUpperCase() != 'N' ) {
    info.setValue( -3 );
  } else if ( m < 0 ) info.setValue( -4 );
  else if ( p < 0 ) info.setValue( -5 );
  else if ( n < 0 ) info.setValue( -6 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -8 );
  else if ( ldb < Math.max( 1, p ) ) info.setValue( -10 );
  else if ( ldu < 1 || ( wantu && ldu < m ) ) info.setValue( -16 );
  else if ( ldv < 1 || ( wantv && ldv < p ) ) info.setValue( -18 );
  else if ( ldq < 1 || ( wantq && ldq < n ) ) info.setValue( -20 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dggsvp', - info.getValue() );
    return;
  }
  for ( var i = 1; i <= n; i ++ ) iwork[ ioffiwork + i - 1 ] = 0.;
  LaPack3.dgeqpf( p, n, B, ldb, iwork, tau, work, info, ioffb,
    ioffiwork, iofftau, ioffwork );
  LaPack0.dlapmt( forwrd, m, n, A, lda, iwork, ioffa, ioffiwork );
  l.setValue( 0 );
  for ( i = 1; i <= Math.min( p, n ); i ++ ) {
    if ( Math.abs( B[ ioffb + i - 1 + ( i - 1 ) * ldb ] ) > tolb ) {
      l.setValue( l.getValue() + 1 );
    }
  }
  if ( wantv ) {
    LaPack0.dlaset( 'Full', p, p, 0., 0., V, ldv, ioffv );
    if ( p > 1 ) {
      LaPack0.dlacpy( 'Lower', p - 1, n, B, ldb, V, ldv, ioffb + 1,
        ioffv + 1 );
    }
    LaPack2.dorg2r( p, p, Math.min( p, n ), V, ldv, tau, work, info,
      ioffv, iofftau, ioffwork );
  }
  for ( var j = 1; j <= l.getValue() - 1; j ++ ) {
    for ( i = j + 1; i <= l.getValue(); i ++ ) {
      B[ ioffb + i - 1 + ( j - 1 ) * ldb ] = 0.;
    }
  }
  if ( p > l.getValue() ) {
    LaPack0.dlaset( 'Full', p - l.getValue(), n, 0., 0., B, ldb,
      ioffb + l.getValue() );
  }
  if ( wantq ) {
    LaPack0.dlaset( 'Full', n, n, 0., 1., Q, ldq, ioffq );
    LaPack0.dlapmt( forwrd, n, n, Q, ldq, iwork, ioffq, ioffiwork );
  }
  if ( p >= l.getValue() && n != l.getValue() ) {
    LaPack2.dgerq2( l.getValue(), n, B, ldb, tau, work, info, ioffb,
      iofftau, ioffwork );
    LaPack2.dormr2( 'Right', 'Transpose', m, n, l.getValue(), B, ldb, tau,
      A, lda, work, info, ioffb, iofftau, ioffa, ioffwork );
    if ( wantq ) {
      LaPack2.dormr2( 'Right', 'Transpose', n, n, l.getValue(), B, ldb, tau,
        Q, ldq, work, info, ioffb, iofftau, ioffq, ioffwork );
    }
    LaPack0.dlaset( 'Full', l.getValue(), n - l.getValue(), 0., 0., B, ldb,
      ioffb );
    for ( j = n - l.getValue() + 1; j <= n; j ++ ) {
      for ( i = j - n + l.getValue() + 1; i <= l.getValue(); i ++ ) {
        B[ ioffb + i - 1 + ( j - 1 ) * ldb ] = 0.;
      }
    }
  }
  for ( i = 1; i <= n - l.getValue(); i ++ ) {
    iwork[ ioffiwork + i - 1 ] = 0.;
  }
  LaPack3.dgeqpf( m, n - l.getValue(), A, lda, iwork, tau, work, info,
    ioffa, ioffiwork, iofftau, ioffwork );
  k.setValue( 0 );
  for ( i = 1; i <= Math.min( m, n - l.getValue() ); i ++ ) {
    if ( Math.abs( A[ ioffa + i - 1 + ( i - 1 ) * lda ] ) > tola ) {
      k.setValue( k.getValue() + 1 );
    }
  }
  LaPack2.dorm2r( 'Left', 'Transpose', m, l.getValue(),
    Math.min( m, n - l.getValue() ), A, lda, tau, A, lda, work, info, ioffa,
    iofftau, ioffa + ( n - l.getValue() ) * lda, ioffwork );
  if ( wantu ) {
    LaPack0.dlaset( 'Full', m, m, 0., 0., U, ldu, ioffu );
    if ( m > 1 ) {
      LaPack0.dlacpy( 'Lower', m - 1, n - l.getValue(), A, lda, U, ldu,
        ioffa + 1, ioffu + 1 );
    }
    LaPack2.dorg2r( m, m, Math.min( m, n - l.getValue() ), U, ldu, tau,
      work, info, ioffu, iofftau, ioffwork );
  }
  if ( wantq ) {
    LaPack0.dlapmt( forwrd, n, n - l.getValue(), Q, ldq, iwork,
      ioffq, ioffiwork );
  }
  for ( j = 1; j <= k.getValue() - 1; j ++ ) {
    for ( i = j + 1; i <= k.getValue(); i ++ ) {
      A[ ioffa + i - 1 + ( j - 1 ) * lda ] = 0.;
    }
  }
  if ( m > k.getValue() ) {
    LaPack0.dlaset( 'Full', m - k.getValue(), n - l.getValue(), 0., 0., A, lda,
      ioffa + k.getValue() );
  }
  if ( n - l.getValue() > k.getValue() ) {
    LaPack2.dgerq2( k.getValue(), n - l.getValue(), A, lda, tau, work, info,
      ioffa, iofftau, ioffwork );
    if ( wantq ) {
      LaPack2.dormr2( 'Right', 'Transpose', n, n - l.getValue(), k.getValue(),
        A, lda, tau, Q, ldq, work, info, ioffa, iofftau, ioffq,
        ioffwork );
    }
    LaPack0.dlaset( 'Full', k.getValue(), n - l.getValue() - k.getValue(), 0., 0.,
      A, lda, ioffa );
    for ( j = n - l.getValue() - k.getValue() + 1; j <= n - l.getValue(); j ++ ) {
      for ( i = j - n + l.getValue() + k.getValue() + 1; i <= k.getValue(); i ++ ) {
        A[ ioffa + i - 1 + ( j - 1 ) * lda ] = 0.;
      }
    }
  }
  if ( m > k.getValue() ) {
    LaPack2.dgeqr2( m - k.getValue(), l.getValue(), A, lda, tau, work, info,
      ioffa + k.getValue() + ( n - l.getValue() ) * lda , iofftau + ioffa,
      ioffwork );
    if ( wantu ) {
      LaPack2.dorm2r( 'Right', 'No transpose', m, m - k.getValue(),
        Math.min( m - k.getValue(), l.getValue() ), A, lda, tau, U, ldu, work,
        info, ioffa + k.getValue() + ( n - l.getValue() ) * lda, iofftau,
        ioffu + k.getValue() * ldu, ioffwork );
    }
    for ( j = n - l.getValue() + 1; j <= n; j ++ ) {
      for ( i = j - n + k.getValue() + l.getValue() + 1; i <= m; i ++ ) {
        A[ ioffa + i - 1 + ( j - 1 ) * lda ] = 0.;
      }
    }
  }
}
LaPack4.zggsvp = function( jobu, jobv, jobq, m, p, n, A, lda, B,
ldb, tola, tolb, k, l, U, ldu, V, ldv, Q, ldq, iwork, tau, work, info ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack4.dlaed1 = function( n, d, Q, ldq, indxq, rho, cutpnt,
work, iwork, info, ioffd, ioffq, ioffindxq, ioffwork, ioffiwork ) {
  info.setValue( 0 );
  if ( n < 0 ) info.setValue( -1 );
  else if ( ldq < Math.max( 1, n ) ) info.setValue( -4 );
  else if ( Math.min( 1, n / 2 ) > cutpnt || n / 2 < cutpnt ) {
    info.setValue( -7 );
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlaed1', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  var iz = 1;
  var idlmda = iz + n;
  var iw = idlmda + n;
  var iq2 = iw + n;
  var indx = 1;
  var indxc = indx + n;
  var coltyp = indxc + n;
  var indxp = coltyp + n;
  Blas1.dcopy( cutpnt, Q, ldq, work, 1, ioffq + cutpnt - 1,
    ioffwork + iz - 1 );
  var zpp1 = cutpnt + 1;
  Blas1.dcopy( n - cutpnt, Q, ldq, work, 1,
    ioffq + zpp1 - 1 + ( zpp1 - 1 ) * ldq,
    ioffwork + iz + cutpnt - 1 );
  var k = new IntReference();
  var rhoref = new NumberReference( rho );
  LaPack1.dlaed2( k, n, cutpnt, d, Q, ldq, indxq, rhoref, work, work,
    work, work, iwork, iwork, iwork, iwork,
    info, ioffd, ioffq, ioffindxq, ioffwork + iz - 1,
    ioffwork + idlmda - 1, ioffwork + iw - 1, ioffwork + iq2 - 1,
    ioffiwork + indx - 1, ioffiwork + indxc - 1, ioffiwork + indxp - 1,
    ioffiwork + coltyp - 1 );
  rho = rhoref.getValue();
  if ( info.getValue() != 0 ) return;
  if ( k.getValue() != 0 ) {
    var iss = ( iwork[ ioffiwork + coltyp - 1 ]
      + iwork[ ioffiwork + coltyp ] ) * cutpnt
      + ( iwork[ ioffwork + coltyp ] + iwork[ ioffwork + coltyp + 1 ] )
      * ( n - cutpnt ) + iq2;
    LaPack3.dlaed3( k.getValue(), n, cutpnt, d, Q, ldq, rho,
      work, work, iwork, iwork, work, work, info,
      ioffd, ioffq, ioffwork + idlmda - 1, ioffwork + iq2 - 1,
      ioffiwork + indxc - 1, ioffiwork + coltyp - 1, ioffwork + iw - 1,
      ioffwork + iss - 1 );
    if ( info.getValue() != 0 ) return;
    var n1 = k.getValue();
    var n2 = n - k.getValue();
    LaPack0.dlamrg( n1, n2, d, 1, -1, indxq, ioffd, ioffindxq );
  } else {
    for ( var i = 1; i <= n; i ++ ) indxq[ ioffindxq + i - 1 ] = i;
  }
}
//*************************************************************************
LaPack4.dlaed7 = function( icompq, n, qsiz, tlvls, curlvl,
curpbm, d, Q, ldq, indxq, rho, cutpnt, qstore, qptr, prmptr, perm, givptr,
givcol, givnum, work, iwork, info, ioffd, ioffq, ioffindxq, ioffqstore,
ioffqptr, ioffprmptr, ioffperm, ioffgivptr, ioffgivcol, ioffgivnum,
ioffwork, ioffiwork ) {
  throw new Error("not tested");
  info.setValue( 0 );
  if ( icompq < 0 || icompq > 1 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  else if ( icompq == 1 && qsiz < n ) info.setValue( -4 );
  else if ( ldq < Math.max( 1, n ) ) info.setValue( -9 );
  else if ( Math.min( 1, n ) > cutpnt || n < cutpnt ) info.setValue( -12 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlaed7', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  var ldq2 = ( icompq == 1 ? qsiz : n );
  var iz = 1;
  var idlmda = iz + n;
  var iw = idlmda + n;
  var iq2 = iw + n;
  var iss = iq2 + n * ldq2;
  var indx = 1;
  var indxc = indx + n;
  var coltyp = indxc + n;
  var indxp = coltyp + n;
  var ptr = 1 + Math.pow( 2, tlvls );
  for ( var i = 1; i <= curlvl - 1; i ++ ) {
    ptr += Math.pow( 2, tlvls - i );
  }
  var curr = ptr + curpbm;
  LaPack0.dlaeda( n, tlvls, curlvl, curpbm, prmptr, perm, givptr,
    givcol, givnum, qstore, qptr, work, work, info, ioffprmptr,
    ioffperm, ioffgivptr, ioffgivcol, ioffgivnum, ioffqstore,
    ioffqptr, ioffwork + iz - 1, ioffwork + iz + n - 1 );
  if ( curlvl == tlvls ) {
    qptr[ ioffqptr + curr - 1 ] = 1;
    prmptr[ ioffprmptr + curr - 1 ] = 1;
    givptr[ ioffgivptr + curr - 1 ] = 1;
  }
  var k = new IntReference();
  var rhorefReference = new NumberReference( rho );
  var givptrref =
    new IntReference( givptr[ ioffgivptr + curr ] );
  LaPack1.dlaed8( icompq, k, n, qsiz, d, Q, ldq, indxq, rhoref, cutpnt,
    work, work, work, ldq2, work, perm, givptrref, givcol, givnum,
    iwork, iwork, info, ioffd, ioffq, ioffindxq, ioffwork + iz - 1,
    ioffwork + idlmda - 1, ioffwork + iq2 - 1, ioffwork + iw - 1,
    ioffperm + prmptr[ ioffprmptr + curr - 1 ] - 1,
    ioffgivcol + ( givptr[ ioffgivptr + curr - 1 ] - 1 ) * 2,
    ioffgivnum + ( givptr[ ioffgivptr + curr - 1 ] - 1 ) * 2,
    ioffiwork + indxp - 1, ioffiwork + indx - 1 );
  rho = rhoref.getValue();
  givptr[ ioffgivptr + curr ] = givptrref.getValue();
  prmptr[ ioffprmptr + curr ] = prmptr[ ioffprmptr + curr - 1 ] + n;
  givptr[ ioffgivptr + curr ] += givptr[ ioffgivptr + curr - 1 ];
  if ( k.getValue() != 0 ) {
    LaPack3.dlaed9( k.getValue(), 1, k.getValue(), n, d, work, k.getValue(), rho,
      work, work, qstore, k.getValue(), info, ioffd, ioffwork + iss - 1,
      ioffwork + idlmda - 1, ioffwork + iw - 1,
      ioffqstore + qptr[ ioffqptr + curr - 1 ] -1 );
    if ( info.getValue() != 0 ) return;
    if ( icompq == 1 ) {
      Blas3.dgemm( 'N', 'N', qsiz, k.getValue(), k.getValue(), 1., work, ldq2,
        qstore, k.getValue(), 0., Q, ldq, ioffwork + iq2 - 1,
        ioffqstore + qptr[ ioffqptr + curr - 1 ] - 1, ioffq );
    }
    qptr[ ioffqptr + curr ] = qptr[ ioffqptr + curr - 1 ]
      + k.getValue() * k.getValue();
    var n1 = k.getValue();
    var n2 = n - k.getValue();
    LaPack0.dlamrg( n1, n2, d, 1, -1, indxq, ioffd, ioffindxq );
  } else {
    qptr[ ioffqptr + curr ] = qptr[ ioffqptr + curr - 1];
    for ( i = 1; i <= n; i ++ ) indxq[ ioffindxq + i - 1 ] = i;
  }
}
//*************************************************************************
LaPack4.dlarre = function( range, n, vlReference, vuReference,
il, iu, d, e, e2, rtol1, rtol2, spltol, nsplit, isplit, m, w, werr, wgap,
iblock, indexw, gers, pivminReference, work, iwork, info, ioffd, ioffe,
ioffe2, ioffisplit, ioffw, ioffwerr, ioffwgap, ioffiblock, ioffindexw,
ioffgers, ioffwork, ioffiwork ) {
  throw new Error("not tested");
  var fac = 0.5;
  var pert = 8.;
  var maxgrowth = 64;
  var fudge = 2.;
  var maxtry = 6;
  var allrng = 1;
  var indrng = 2;
  var valrng = 3;
  var iseed = new Array( 4 );
  info.setValue( 0 );
  if ( range.charAt(0).toUpperCase() == 'A' ) var irange = allrng;
  else if ( range.charAt(0).toUpperCase() == 'V' ) irange = valrng;
  else if ( range.charAt(0).toUpperCase() == 'I' ) irange = indrng;
  m.setValue( 0 );
  var safmin = LaPack0.dlamch( 'S' );
  var eps = LaPack0.dlamch( 'P' );
  var rtl = Math.sqrt( eps );
  var bsrtol = Math.sqrt( eps );
  if ( n == 1 ) {
    if ( irange == allrng ||
    ( irange == valrng && d[ ioffd ] > vl.getValue() &&
    d[ ioffd ] <= vu.getValue() ) ||
    ( irange == indrng && il == 1 && iu == 1 ) ) {
      m.setValue( 1 );
      w[ ioffw ] = d[ ioffd ];
      werr[ ioffwerr ] = 0.;
      wgap[ ioffwgap ] = 0.;
      iblock[ ioffiblock ] = 1;
      indexw[ ioffindexw ] = 1;
      gers[ ioffgers ] = d[ ioffd ];
      gers[ ioffgers + 1 ] = d[ ioffd ];
    }
    e[ ioffe ] = 0.;
    return;
  }
  var gl = d[ ioffd ];
  var gu = d[ ioffd ];
  var eold = 0.;
  var emax = 0.;
  e[ ioffe + n - 1 ] = 0.;
  for ( var i = 1; i <= n; i ++ ) {
    werr[ ioffwerr + i - 1 ] = 0.;
    wgap[ ioffwgap + i - 1 ] = 0.;
    var eabs = Math.abs( e[ ioffe + i - 1 ] );
    if ( eabs >= emax ) emax = eabs;
    var tmp1 = eabs + eold;
    gers[ ioffgers + 2 * i - 2 ] = d[ ioffd + i - 1 ] - tmp1;
    gl = Math.min( gl, gers[ ioffgers + 2 * i - 2 ] );
    gers[ ioffgers + 2 * i - 1 ] = d[ ioffd + i - 1 ] + tmp1;
    gu = Math.max( gu, gers[ ioffgers + 2 * i - 1 ] );
    eold = eabs;
  } // 5
  pivmin.setValue( safmin * Math.max( 1., emax * emax ) );
  var spdiam = gu - gl;
  var iinfo = new IntReference();
  LaPack0.dlarra( n, d, e, e2, spltol, spdiam, nsplit, isplit, iinfo,
    ioffd, ioffe, ioffe2, ioffisplit );
  var forceb = false;
  var usedqd = ( irange == allrng && ! forceb );
  if ( irange == allrng && ! forceb ) {
    vl.setValue( gl );
    vu.setValue( gu );
  } else {
    var mm = new IntReference();
    LaPack1.dlarrd( range, 'B', n, vl.getValue(), vu.getValue(), il, iu,
      gers, bsrtol, d, e, e2, pivmin.getValue(), nsplit.getValue(),
      isplit, mm, w, werr, vl, vu, iblock, indexw, work, iwork, iinfo,
      ioffgers, ioffd, ioffe, ioffe2, ioffisplit, ioffw, ioffwerr,
      ioffiblock, ioffindexw, ioffwork, ioffiwork );
    if ( iinfo.getValue() != 0 ) {
      info.setValue( -1 );
      return;
    }
    for ( i = mm.getValue() + 1; i <= n; i ++ ) {
      w[ ioffw + i - 1 ] = 0.;
      werr[ ioffwerr + i - 1 ] = 0.;
      iblock[ ioffiblock + i - 1 ] = 0;
      indexw[ ioffindexw + i - 1 ] = 0;
    } // 14
  }
  var ibegin = 1;
  var wbegin = 1;
  for ( var jblk = 1; jblk <= nsplit.getValue(); jblk ++ ) {
    var iend = isplit[ ioffisplit + jblk - 1 ];
    var in2 = iend - ibegin + 1;
    if ( in2 == 1 ) {
      if ( irange == allrng ||
      ( irange == valrng && d[ ioffd + ibegin - 1 ] > vl &&
      d[ ioffd + ibegin - 1 ] <= vu.getValue() ) ||
      ( irange == indrng && iblock[ ioffiblock + wbegin - 1 ] == jblk )
      ) {
        m.setValue( m.getValue() + 1 );
        w[ ioffw + m.getValue() - 1 ] = d[ ioffd + ibegin - 1 ];
        werr[ ioffwerr + m.getValue() - 1 ] = 0.;
        wgap[ ioffwgap + m.getValue() - 1 ] = 0.;
        iblock[ ioffiblock + m.getValue() - 1 ] = jblk;
        indexw[ ioffindexw + m.getValue() - 1 ] = 1;
        wbegin ++;
      }
      e[ ioffe + iend - 1 ] = 0.;
      ibegin = iend + 1;
      continue; // goto 170
    }
    e[ ioffe + iend - 1 ] = 0.;
    gl = d[ ioffd + ibegin - 1 ];
    gu = d[ ioffd + ibegin - 1 ];
    for ( i = ibegin; i <= iend; i ++ ) {
      gl = Math.min( gers[ ioffgers + 2 * i - 2 ], gl );
      gu = Math.max( gers[ ioffgers + 2 * i - 1 ], gu );
    } // 15
    spdiam = gu - gl;
    if ( ! ( irange == allrng && ! forceb ) ) {
      var mb = 0;
      for ( i = wbegin; i <= mm.getValue(); i ++ ) {
        if ( iblock[ ioffiblock + i - 1 ] == jblk ) mb ++;
        else break;
      } // 20
      if ( mb == 0 ) {
        e[ ioffe + iend - 1 ] = 0.;
        ibegin = iend + 1;
        continue; // goto 170
      } else {
        usedqd = ( mb > fac * in2 && ! forceb );
        var wend = wbegin + mb - 1;
        var sigma = 0.;
        for ( i = wbegin; i <= wend - 1; i ++ ) {
          wgap[ ioffwgap + i - 1 ] = Math.max( 0.,
            w[ ioffw + i ] - werr[ ioffwerr ]
            - ( w[ ioffw + i - 1 ] + werr[ ioffwerr + i - 1 ] ) );
        } // 30
        wgap[ ioffwgap + wend - 1 ] = Math.max( 0., vu.getValue() - sigma
          - ( w[ ioffw + wend - 1 ] + werr[ ioffwerr + wend - 1 ] ) );
        var indl = indexw[ ioffindexw + wbegin - 1 ];
        var indu = indexw[ ioffindexw + wend - 1 ];
      }
    }
    var tmp = Number.POSITIVE_INFINITY;
    if ( ( irange == allrng && ! forceb ) || usedqd ) {
      var tmprefReference = new NumberReference( tmp );
      var tmp1refReference = new NumberReference();
      LaPack1.dlarrk( in2, 1, gl, gu, d, e2, pivmin.getValue(), rtl,
        tmpref, tmp1ref, iinfo, ioffd + ibegin - 1,
        ioffe2 + ibegin - 1 );
      tmp = tmpref.getValue();
      tmp1 = tmp1ref.getValue();
      if ( iinfo.getValue() != 0 ) {
        info.setValue( -1 );
        return;
      }
      var isleft = Math.max( gl,
        tmp - tmp1 - 100. * eps * Math.abs( tmp - tmp1 ) );
      LaPack1.dlarrk( in2, in2, gl, gu, d, e2, pivmin.getValue(), rtl,
        tmpref, tmp1ref, iinfo, ioffd + ibegin - 1,
        ioffe2 + ibegin - 1 );
      tmp = tmpref.getValue();
      tmp1 = tmp1ref.getValue();
      if ( iinfo.getValue() != 0 ) {
        info.setValue( -1 );
        return;
      }
      var isrght = Math.min( gu,
        tmp + tmp1 + 100. * eps * Math.abs( tmp + tmp1 ) );
      spdiam = isrght - isleft;
    } else {
      isleft = Math.max( gl, w[ ioffw + wbegin - 1 ]
        - werr[ ioffwerr + wbegin - 1 ]
        - 100. * eps * Math.abs( w[ ioffw + wbegin - 1 ]
        - werr[ ioffwerr + wbegin - 1 ] ) );
      isrght = Math.min( gu, w[ ioffw + wend - 1 ]
        + werr[ ioffwerr + wend - 1 ]
        + 100. * eps * Math.abs( w[ ioffw + wend - 1 ]
        + werr[ ioffwerr + wend - 1 ] ) );
    }
    if ( irange == allrng && ! forceb ) {
      usedqd = true;
      indl = 1;
      indu = in2;
      mb = in2;
      wend = wbegin + mb - 1;
      var s1 = isleft + 0.25 * spdiam;
      var s2 = isrght - 0.25 * spdiam;
    } else {
      if ( usedqd ) {
        s1 = isleft + 0.25 * spdiam;
        s2 = isrght - 0.25 * spdiam;
      } else {
        tmp = Math.min( isrght, vu.getValue() )
          - Math.max( isleft, vl.getValue() );
        s1 = Math.max( isleft, vl.getValue() ) + 0.25 * tmp;
        s2 = Math.min( isrght, vu.getValue() ) - 0.25 * tmp;
      }
    }
    if ( mb > 1 ) {
      var cnt = new IntReference();
      var cnt1 = new IntReference();
      var cnt2 = new IntReference();
      LaPack0.dlarrc( 'T', in2, s1, s2, d, e, pivmin.getValue(), cnt, cnt1,
        cnt2, iinfo, ioffd + ibegin - 1, ioffe + ibegin - 1 );
    }
    if ( mb == 1 ) {
      sigma = gl;
      var sgndef = 1.;
    } else if ( cnt1.getValue() - indl >= indu - cnt2.getValue() ) {
      if ( irange == allrng && ! forceb ) {
        sigma = Math.max( isleft, gl );
      } else if ( usedqd ) sigma = isleft;
      else sigma = Math.max( isleft, vl.getValue() );
      sgndef = 1.;
    } else {
      if ( irange == allrng && ! forceb ) {
        sigma = Math.min( isrght, gu );
      } else if ( usedqd ) sigma = isrght;
      else sigma = Math.min( isrght, vu.getValue() );
      sgndef = -1.;
    }
    if ( usedqd ) {
      var tau = spdiam * eps * n + 2. * pivmin.getValue();
      tau = Math.max( tau, 2. * eps * Math.abs( sigma ) );
    } else {
      if ( mb > 1 ) {
        var clwdth = w[ ioffw + wend - 1 ]
          + werr[ ioffwerr + wend - 1 ]
          - w[ ioffw + wbegin - 1 ]
          - werr[ ioffwerr + wbegin - 1 ];
        var avgap =
          Math.abs( clwdth / Number( wend - wbegin ) );
        if ( sgndef == 1. ) {
          tau = 0.5 * Math.max( wgap[ ioffwgap + wbegin - 1 ], avgap );
          tau = Math.max( tau, werr[ ioffwerr + wbegin - 1 ] );
        } else {
          tau = 0.5 * Math.max( wgap[ ioffwgap + wend - 1 ], avgap );
          tau = Math.max( tau, werr[ ioffwerr + wend - 1 ] );
        }
      } else {
        tau = werr[ ioffwerr + wbegin - 1 ];
      }
    }
    var goto83 = false;
    for ( var idum = 1; idum <= maxtry; idum ++ ) {
      var dpivot = d[ ioffd + ibegin - 1 ] - sigma;
      work[ ioffwork ] = dpivot;
      var dmax = Math.abs( work[ ioffwork ] );
      var j = ibegin;
      for ( i = 1; i <= in2 - 1; i ++ ) {
        work[ ioffwork + 2 * in2 + i - 1 ] =
          1. /work[ ioffwork + i - 1 ];
        tmp = e[ ioffe + j - 1 ] * work[ ioffwork + 2 * in2 + i - 1 ];
        work[ ioffwork + in2 + i - 1 ] = tmp;
        dpivot = ( d[ ioffd + j ] - sigma ) - tmp * e[ ioffe + j - 1 ];
        work[ ioffwork + i ] = dpivot;
        dmax = Math.max( dmax, Math.abs( dpivot ) );
        j ++;
      } // 70
      var norep = ( dmax > maxgrowth * spdiam ? true : false );
      if ( usedqd && ! norep ) {
        for ( i = 1; i <= in2; i ++ ) {
          tmp = sgndef * work[ ioffwork + i - 1 ];
          if ( tmp < 0. ) norep = true;
        } // 71
      }
      if ( norep ) {
        if ( idum == maxtry - 1 ) {
          sigma = ( sgndef == 1. ?
            gl - fudge * spdiam * eps * n - fudge * 2. * pivmin.getValue() :
            gu + fudge * spdiam * eps * n + fudge * 2. * pivmin.getValue()
            );
        } else {
          sigma -= sgndef * tau;
          tau *= 2.;
        }
      } else {
        goto83 = true;
        break;
      }
    } // 80
    if ( ! goto83 ) {
      info.setValue( 2 );
      return;
    }
    e[ ioffe + iend - 1 ] = sigma;
    Blas1.dcopy( in2, work, 1, d, 1, ioffwork, ioffd + ibegin - 1 );
    Blas1.dcopy( in2 - 1, work, 1, e, 1,
      ioffwork + in2, ioffe + ibegin - 1 );
    if ( mb > 1 ) {
      for ( i = 1; i <= 4; i ++ ) iseed[ i - 1 ] = 1;
      LaPack0.dlarnv( 2, iseed, 2 * in2 - 1, work, 0, ioffwork );
      for ( i = 1; i <= in2 - 1; i ++ ) {
        d[ ioffd + ibegin + i - 2 ] *=
          1. + eps * pert * work[ ioffwork + i - 1 ];
        e[ ioffe + ibegin + i - 2 ] *=
          1. + eps * pert * work[ ioffwork + in2 + i - 1 ];
      }
      d[ ioffd + iend + i - 1 ] *=
        1. + eps * 4. * work[ ioffwork + in2 - 1 ];
    }
    if ( ! usedqd ) {
      for ( j = wbegin; j <= wend; j ++ ) {
        w[ ioffw + j - 1 ] -= sigma;
        werr[ ioffwerr + j - 1 ] +=
          Math.abs( w[ ioffw + j - 1 ] ) * eps;
      } // 134
      for ( i = ibegin; i <= iend - 1; i ++ ) {
        work[ ioffwork + i - 1 ] =
          d[ ioffd + i - 1 ] * Math.pow( e[ ioffe + i - 1 ], 2 );
      } // 135
      LaPack1.dlarrb( in2, d, work, indl, indu, rtol1, rtol2,
        indl - 1, w, wgap, werr, work, iwork, pivmin.getValue(), spdiam,
        in2, iinfo, ioffd + ibegin - 1, ioffwork + ibegin - 1,
        ioffw + wbegin - 1, ioffwgap + wbegin - 1,
        ioffwerr + wbegin - 1, ioffwork + 2 * n, ioffiwork );
      if ( iinfo.getValue() != 0 ) {
        info.setValue( -4 );
        return;
      }
      wgap[ ioffwgap + wend - 1 ] = Math.max( 0., ( vu.getValue() - sigma )
        - ( w[ ioffw + wend - 1 ] + werr[ ioffwerr + wend - 1 ] ) );
      for ( i = indl; i <= indu; i ++ ) {
        m.setValue( m.getValue() + 1 );
        iblock[ ioffiblock + m.getValue() - 1 ] = jblk;
        indexw[ ioffindexw + m.getValue() - 1 ] = i;
      } // 138
    } else {
      var rtol = Math.log( Number( in2 ) ) * 4. * eps;
      j = ibegin;
      for ( i = 1; i <= in2 - 1; i ++ ) {
        work[ ioffwork + 2 * i - 2 ] = Math.abs( d[ ioffd + j - 1 ] );
        work[ ioffwork + 2 * i - 1 ] =
          e[ ioffe + j - 1 ] * e[ ioffe + j - 1 ]
          * work[ ioffwork + 2 * i - 2 ];
        j ++;
      } // 140
      work[ ioffwork + 2 * in2 - 2 ] =
        Math.abs( d[ ioffd + iend - 1 ] );
      work[ ioffwork + 2 * in2 - 1 ] = 0.;
      LaPack3.dlasq2( in2, work, iinfo, ioffwork );
      if ( iinfo.getValue() != 0 ) {
        info.setValue( -5 );
        return;
      } else {
        for ( i = 1; i <= in2; i ++ ) {
          if ( work[ ioffwork + i - 1 ] < 0. ) {
            info.setValue( -6 );
            return;
          }
        } // 149
      }
      if ( sgndef > 0. ) {
        for ( i = indl; i <= indu; i ++ ) {
          m.setValue( m.getValue() + 1 );
          w[ ioffw + m.getValue() - 1 ] = work[ ioffwork + in2 - i ];
          iblock[ ioffiblock + m.getValue() - 1 ] = jblk;
          indexw[ ioffindexw + m.getValue() - 1 ] = i;
        } // 150
      } else {
        for ( i = indl; i <= indu; i ++ ) {
          m.setValue( m.getValue() + 1 );
          w[ ioffw + m.getValue() - 1 ] = - work[ ioffwork + i - 1 ];
          iblock[ ioffiblock + m.getValue() - 1 ] = jblk;
          indexw[ ioffindexw + m.getValue() - 1 ] = i;
        } // 160
      }
      for ( i = m.getValue() - mb + 1; i <= m.getValue(); i ++ ) {
        werr[ ioffwerr + i - 1 ] =
          rtol * Math.abs( w[ ioffw + i - 1 ] );
      } // 165
      for ( i = m.getValue() - mb + 1; i <= m.getValue() - 1; i ++ ) {
        wgap[ ioffwgap + i - 1 ] = Math.max( 0., w[ ioffw + i ]
          - werr[ ioffwerr + i ] - ( w[ ioffw + i - 1 ]
          + werr[ ioffwerr + i - 1 ] ) );
      } // 166
      wgap[ ioffwgap + m.getValue() - 1 ] = Math.max( 0.,
        ( vu.getValue() - sigma ) - ( w[ ioffw + m.getValue() - 1 ]
        + werr[ ioffwerr + m.getValue() - 1 ] ) );
    }
    ibegin = iend + 1;
    wbegin = wend + 1;
  } // 170
}
//*************************************************************************
LaPack4.dlasd1 = function( nl, nr, sqre, d, alphaReference,
betaReference, U, ldu, Vt, ldvt, idxq, iwork, work, info, ioffd, ioffu,
ioffvt, ioffidxq, ioffiwork, ioffwork ) {
  throw new Error("not tested");
  info.setValue( 0 );
  if ( nl < 1 ) info.setValue( -1 );
  else if ( nr < 1 ) info.setValue( -2 );
  else if ( sqre < 0 || sqre > 1 ) info.setValue( -3 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlasd1', - info.getValue() );
    return;
  }
  var n = nl + nr + 1;
  var m = n + sqre;
  var ldu2 = n;
  var ldvt2 = m;
  var iz = 1;
  var isigma = iz + m;
  var iu2 = isigma + n;
  var ivt2 = iu2 + ldu2 * n;
  var iq = ivt2 + ldvt2 * m;
  var idx = 1;
  var idxc = idx + n;
  var coltyp = idxc + n;
  var idxp = coltyp + n;
  var orgnrm =
    Math.max( Math.abs( alpha.getValue() ), Math.abs( beta.getValue() ) );
  d[ ioffd + nl ] = 0.;
  for ( var i = 1; i <= n; i ++ ) {
    if ( Math.abs( d[ ioffd + i - 1 ] ) > orgnrm ) {
      orgnrm = Math.abs( d[ ioffd + i - 1 ] );
    }
  } // 10
  LaPack1.dlascl( 'G', 0, 0, orgnrm, 1., n, 1, d, n, info, ioffd );
  alpha.setValue( alpha.getValue() / orgnrm );
  beta.setValue( beta.getValue() / orgnrm );
  var k = new IntReference();
  LaPack1.dlasd2( nl, nr, sqre,
    k, d, work, alpha.getValue(), beta.getValue(),
    U, ldu, Vt, ldvt, work, work,
    ldu2, work, ldvt2, iwork, iwork, iwork,
    idxq, iwork, info, ioffd, ioffwork + iz - 1,
    ioffu, ioffvt, ioffwork + isigma - 1, ioffwork + iu2 - 1,
      ioffwork + ivt2 - 1,
    ioffiwork + idxp - 1, ioffiwork + idx - 1, ioffiwork + idxc - 1,
      ioffidxq, ioffiwork + coltyp - 1 );
  var ldq = k.getValue();
  LaPack3.dlasd3( nl, nr, sqre, k.getValue(),
    d, work, ldq, work, U, ldu, work,
    ldu2, Vt, ldvt, work, ldvt2, iwork,
    iwork, work, info, ioffd, ioffwork + iq - 1,
    ioffwork + isigma - 1, ioffu, ioffwork + iu2 - 1, ioffvt,
      ioffwork + ivt2 - 1,
    ioffiwork + idxc - 1, ioffiwork + coltyp - 1, ioffwork + iz - 1 );
  if ( info.getValue() != 0 ) return;
  LaPack1.dlascl( 'G', 0, 0, 1., orgnrm, n, 1, d, n, info, ioffd );
  var n1 = k.getValue();
  var n2 = n - k.getValue();
  LaPack0.dlamrg( n1, n2, d, 1, -1, idxq, ioffd, ioffidxq );
}
//*************************************************************************
LaPack4.dlasd6 = function( icompq, nl, nr, sqre, d, vf, vl,
alphaReference, betaReference, idxq, perm, givptr, givcol, ldgcol, givnum,
ldgnum, poles, difl, difr, z, k, cReference, sReference, work, iwork, info,
ioffd, ioffvf, ioffvl, ioffidxq, ioffperm, ioffgivcol, ioffgivnum,
ioffpoles, ioffdifl, ioffdifr, ioffz, ioffwork, ioffiwork ) {
  throw new Error("not tested");
  info.setValue( 0 );
  var n = nl + nr + 1;
  var m = n + sqre;
  if ( icompq < 0 || icompq > 1 ) info.setValue( -1 );
  else if ( nl < 1 ) info.setValue( -2 );
  else if ( rn < 1 ) info.setValue( -3 );
  else if ( sqre < 0 || sqre > 1 ) info.setValue( -4 );
  else if ( ldgcol < n ) info.setValue( -14 );
  else if ( ldgnum < n ) info.setValue( -16 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlasd6', - info.getValue() );
    return;
  }
  var isigma = 1;
  var iw = isigma + n;
  var ivfw = iw + m;
  var ivlw = ivfw + m;
  var idx = 1;
  var idxc = idx + n;
  var idxp = idxc + n;
  var orgnrm =
    Math.max( Math.abs( alpha.getValue() ), Math.abs( beta.getValue() ) );
  d[ ioffd + nl ] = 0.;
  for ( var i = 1; i <= n; i ++ ) {
    if ( Math.abs( d[ ioffd + i - 1 ] ) > orgnrm ) {
      orgnrm = Math.abs( d[ ioffd + i - 1 ] );
    }
  } // 10
  LaPack1.dlascl( 'G', 0, 0, orgnrm, 1., n, 1, d, n, info, ioffd );
  alpha.setValue( alpha.getValue() / orgnrm );
  beta.setValue( beta.getValue() / orgnrm );
  var cReference = new NumberReference();
  var sReference = new NumberReference();
  LaPack1.dlasd7( icompq, nl, nr, sqre, k, d, z, work, vf, work, vl,
    work, alpha.getValue(), beta.getValue(), work, iwork, iwork, idxq, perm,
    givptr, givcol, ldgcol, givnum, ldgnum, c, s, info,
    ioffd, ioffz, ioffwork + iw - 1, ioffvf, ioffwork + ivfw - 1,
    ioffvl, ioffwork + ivlw - 1, ioffwork + isigma - 1,
    ioffiwork + idx - 1, ioffiwork + idxp - 1, ioffidxq, ioffperm,
    ioffgivcol, ioffgivnum );
  LaPack3.dlasd8( icompq, k.getValue(), d, z, vf, vl, difl, difr, ldgnum,
    work, work, info, ioffd, ioffz, ioffvf, ioffvl, ioffdifl, ioffdifr,
    ioffwork + isigma - 1, ioffwork + iw - 1 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlasd8', - info.getValue() );
    return;
  }
  if ( icompq == 1 ) {
   Blas1.dcopy( k.getValue(), d, 1, poles, 1, ioffd, ioffpoles );
   Blas1.dcopy( k.getValue(), work, 1, poles, 1, ioffwork + isigma - 1,
     ioffpoles + ldgnum );
  }
  LaPack1.dlascl( 'G', 0, 0, 1., orgnrm, n, 1, d, n, info, ioffd );
  var n1 = k.getValue();
  var n2 = n - k.getValue();
  LaPack0.dlamrg( n1, n2, d, 1, -1, idxq, ioffd, ioffidxq );
}
//*************************************************************************
LaPack4.dlasq1 = function( n, d, e, work, info, ioffd, ioffe,
ioffwork ) {
  var sigmxReference = new NumberReference();
  info.setValue( 0 );
  if ( n < 0 ) {
    info.setValue( -2 );
    Blas2.xerbla( 'dlasq1', - info.getValue() );
    return;
  } else if ( n == 0 ) return;
  else if ( n == 1 ) {
    d[ ioffd ] = Math.abs( d[ ioffd ] );
    return;
  } else if ( n == 2 ) {
    var sigmnReference = new NumberReference();
    LaPack0.dlas2( d[ ioffd ], e[ ioffe ], d[ ioffd + 1 ], sigmn,
      sigmx );
    d[ ioffd ] = sigmx.getValue();
    d[ ioffd + 1 ] = sigmn.getValue();
    return;
  }
  sigmx.setValue( 0. );
  for ( var i = 1; i <= n - 1; i ++ ) {
    d[ ioffd + i - 1 ] = Math.abs( d[ ioffd + i - 1 ] );
    sigmx.setValue( Math.max( sigmx.getValue(),
      Math.abs( e[ ioffe + i - 1 ] ) ) );
  }
  d[ ioffd + n - 1 ] = Math.abs( d[ ioffd + n - 1 ] );
  var iinfo = new IntReference();
  if ( sigmx.getValue() == 0. ) {
    LaPack0.dlasrt( 'D', n, d, iinfo, ioffd );
    return;
  }
  for ( i = 1; i <= n; i ++ ) {
    sigmx.setValue( Math.max( sigmx.getValue(), d[ ioffd + i - 1 ] ) );
  } // 20
  var eps = LaPack0.dlamch( 'Precision' );
  var safmin = LaPack0.dlamch( 'Safe minimum' );
  var scale = Math.sqrt( eps / safmin );
  Blas1.dcopy( n, d, 1, work, 2, ioffd, ioffwork );
  Blas1.dcopy( n - 1, e, 1, work, 2, ioffe, ioffwork + 1 );
  LaPack1.dlascl( 'G', 0, 0, sigmx.getValue(), scale, 2 * n - 1, 1, work, 
    2 * n - 1, iinfo, ioffwork );
  for ( i = 1; i <= 2 * n - 1; i ++ ) {
    work[ ioffwork + i - 1 ] =
      Math.pow( work[ ioffwork + i - 1 ], 2 );
  } // 30
  work[ ioffwork + 2 * n - 1 ] = 0.;
  LaPack3.dlasq2( n, work, info, ioffwork );
  if ( info.getValue() == 0 ) {
    for ( i = 1; i <= n; i ++ ) {
      d[ ioffd + i - 1 ] = Math.sqrt( work[ ioffwork + i - 1 ] );
    }
    LaPack1.dlascl( 'G', 0, 0, scale, sigmx.getValue(), n, 1, d, n, iinfo,
      ioffd );
  } else if ( info.getValue() == 2 ) {
    for ( i = 1; i <= n; i ++ ) {
      d[ ioffd + i - 1 ] = Math.sqrt( work[ ioffwork + 2 * i - 2 ] );
      e[ ioffd + i - 1 ] = Math.sqrt( work[ ioffwork + 2 * i - 1 ] );
    }
    LaPack1.dlascl( 'G', 0, 0, scale, sigmx.getValue(), n, 1, d, n, iinfo,
      ioffd );
    LaPack1.dlascl( 'G', 0, 0, scale, sigmx.getValue(), n, 1, e, n, iinfo,
      ioffd );
  }
}
//*************************************************************************
LaPack4.dorcsd = function( jobu1, jobu2, jobv1t, jobv2t, trans,
signs, m, p, q, X11, ldx11, X12, ldx12, X21, ldx21, X22, ldx22, theta, U1,
ldu1, U2, ldu2, V1t, ldv1t, V2t, ldv2t, work, lwork, iwork, info, ioffx11,
ioffx12, ioffx21, ioffx22, iofftheta, ioffu1, ioffu2, ioffv1t, ioffv2t,
ioffwork, ioffiwork) {
  var piover2 = 1.57079632679489662;
  info.setValue( 0 );
  var wantu1 = ( jobu1.charAt(0).toUpperCase() == 'Y' );
  var wantu2 = ( jobu2.charAt(0).toUpperCase() == 'Y' );
  var wantv1t = ( jobv1t.charAt(0).toUpperCase() == 'Y' );
  var wantv2t = ( jobv2t.charAt(0).toUpperCase() == 'Y' );
  var colmajor = ( trans.charAt(0).toUpperCase() != 'T' );
  var defaultsigns = ( signs.charAt(0).toUpperCase() != 'O' );
  var lquery = ( lwork == -1 );
  if ( m < 0 ) info.setValue( -7 );
  else if ( p < 0 || p > m ) info.setValue( -8 );
  else if ( q < 0 || q > m ) info.setValue( -9 );
  else if ( ( colmajor && ldx11 < Math.max( 1, p ) ) ||
  ( ! colmajor && ldx11 < Math.max( 1, q ) ) ) {
    info.setValue( -11 );
  } else if ( wantu1 && ldu1 < p ) info.setValue( -20 );
  else if ( wantu2 && ldu2 < m - p ) info.setValue( -22 );
  else if ( wantv1t && ldv1t < q ) info.setValue( -24 );
  else if ( wantv2t && ldv2t < m - q ) info.setValue( -26 );
  if ( info.getValue() == 0 && Math.min( p, m - p ) < Math.min( q, m - q ) )
  {
    var transt = ( colmajor ? 'T' : 'N' ); 
    var signst = ( defaultsigns ? 'O' : 'D' );
    LaPack4.dorcsd( jobv1t, jobv2t, jobu1, jobu2, transt, signst, m,
      q, p, X11, ldx11, X21, ldx21, X12, ldx12, X22, ldx22, theta,
      V1t, ldv1t, V2t, ldv2t, U1, ldu1, U2, ldu2, work, lwork, iwork,
      info, ioffx11, ioffx21, ioffx12, ioffx22, iofftheta, ioffv1t,
      ioffv2t, ioffu1, ioffu2, ioffwork, ioffiwork );
    return;
  }
  if ( info.getValue() == 0 && m - q < q ) {
    signst = ( defaultsigns ? 'O' : 'D' );
    LaPack4.dorcsd( jobu2, jobu1, jobv2t, jobv1t, trans, signst, m,
      m - p, m - q, X22, ldx22, X21, ldx21, X12, ldx12, X11, ldx11,
      theta, U2, ldu2, U1, ldu1, V2t, ldv2t, V1t, ldv1t, work, lwork,
      iwork, info, ioffx22, ioffx21, ioffx12, ioffx11, iofftheta,
      ioffu2, ioffu1, ioffv2t, ioffv1t, ioffwork, ioffiwork );
    return;
  }
  var childinfo = new IntReference();
  var nothing;
  if ( info.getValue() == 0 ) {
    var iphi = 2;
    var itaup1 = iphi + Math.max( 1, q - 1 );
    var itaup2 = itaup1 + Math.max( 1, p );
    var itauq1 = itaup2 + Math.max( 1, m - p );
    var itauq2 = itauq1 + Math.max( 1, q );
    var iorgqr = itauq2 + Math.max( 1, m - q );
    LaPack3.dorgqr( m - q, m - q, m - q, nothing, Math.max( 1, m - q ),
      nothing, work, -1, childinfo, 0, 0, ioffwork );
    var lorgqrworkopt = Math.round( work[ ioffwork ] );
    var lorgqrworkmin = Math.max( 1, m - q );
    var iorglq = itauq2 + Math.max( 1, m - q );
    LaPack3.dorglq( m - q, m - q, m - q, nothing, Math.max( 1, m - q ),
      nothing, work, -1, childinfo, 0, 0, ioffwork );
    var lorglqworkopt = Math.round( work[ ioffwork ] );
    var lorglqworkmin = Math.max( 1, m - q );
    var iorbdb = itauq2 + Math.max( 1, m - q );
    LaPack2.dorbdb( trans, signs, m, p, q, X11, ldx11, X12, ldx12, X21,
      ldx21, X22, ldx22, nothing, nothing, nothing, nothing, nothing,
      nothing, work, -1, childinfo, ioffx11, ioffx12, ioffx21, ioffx22,
      0, 0, 0, 0, 0, 0, ioffwork );
    var lorbdbworkopt = Math.round( work[ ioffwork ] );
    var lorbdbworkmin = lorbdbworkopt;
    var ib11d = itauq2 + Math.max( 1, m - q );
    var ib11e = ib11d + Math.max( 1, q );
    var ib12d = ib11e + Math.max( 1, q - 1 );
    var ib12e = ib12d + Math.max( 1, q );
    var ib21d = ib12e + Math.max( 1, q - 1 );
    var ib21e = ib21d + Math.max( 1, q );
    var ib22d = ib21e + Math.max( 1, q - 1 );
    var ib22e = ib22d + Math.max( 1, q );
    var ibbcsd = ib22e + Math.max( 1, q - 1 );
    LaPack3.dbbcsd( jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q,
      nothing, nothing, U1, ldu1, U2, ldu2, V1t, ldv1t, V2t, ldv2t,
      nothing, nothing, nothing, nothing, nothing, nothing, nothing,
      nothing, work, -1, childinfo, 0, 0, ioffu1, ioffu2, ioffv1t,
      ioffv2t, 0, 0, 0, 0, 0, 0, 0, 0, ioffwork );
    var lbbcsdworkopt = Math.round( work[ ioffwork ] );
    var lbbcsdworkmin = lbbcsdworkopt;
    var lworkopt = Math.max(
      Math.max( iorgqr + lorgqrworkopt, iorglq + lorglqworkopt ),
      Math.max( iorbdb + lorbdbworkopt, ibbcsd + lbbcsdworkopt ) ) - 1;
    var lworkmin = Math.max(
      Math.max( iorgqr + lorgqrworkmin, iorglq + lorglqworkmin ),
      Math.max( iorbdb + lorbdbworkopt, ibbcsd + lbbcsdworkmin ) ) - 1;
    work[ ioffwork ] = Math.max( lworkopt, lworkmin );
    if ( lwork < lworkmin && ! lquery ) info.setValue( -22 );
    else {
      var lorgqrwork = lwork - iorgqr + 1;
      var lorglqwork = lwork - iorglq + 1;
      var lorbdbwork = lwork - iorbdb + 1;
      var lbbcsdwork = lwork - ibbcsd + 1;
    }
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dorcsd', - info.getValue() );
    return;
  } else if ( lquery ) return;
  LaPack2.dorbdb( trans, signs, m, p, q, X11, ldx11, X12, ldx12, X21,
    ldx21, X22, ldx22, theta, work, work, work, work, work, work,
    lorbdbwork, childinfo, ioffx11, ioffx12, ioffx21, ioffx22,
    iofftheta, ioffwork + iphi - 1, ioffwork + itaup1 - 1,
    ioffwork + itaup2 - 1, ioffwork + itauq1 - 1,
    ioffwork + itauq2 - 1, ioffwork + iorbdb - 1 );
  if ( colmajor ) {
    if ( wantu1 && p > 0 ) {
      LaPack0.dlacpy( 'L', p, q, X11, ldx11, U1, ldu1,
        ioffx11, ioffu1 );
      LaPack3.dorgqr( p, p, q, U1, ldu1, work, work, lorgqrwork, info,
        ioffu1, ioffwork + itaup1 - 1, ioffwork + iorgqr - 1 );
    }
    if ( wantu2 && m - p > 0 ) {
      LaPack0.dlacpy( 'L', m - p, q, X21, ldx21, U2, ldu2,
        ioffx21, ioffu2 );
      LaPack3.dorgqr( m - p, m - p, q, U2, ldu2, work, work,
        lorgqrwork, info, ioffu2, ioffwork + itaup2 - 1,
        ioffwork + iorgqr - 1 );
    }
    if ( wantv1t && q > 0 ) {
      LaPack0.dlacpy( 'U', q - 1, q - 1, X11, ldx11, V1t, ldv1t,
        ioffx11 + ldx11, ioffv1t + 1 + ldv1t );
      V1t[ ioffv1t ] = 1.;
      for ( var j = 2; j <= q; j ++ ) {
        V1t[ ioffv1t + ( j - 1 ) * ldv1t ] = 0.;
        V1t[ ioffv1t + j - 1 ] = 0.;
      }
      LaPack3.dorglq( q - 1, q - 1, q - 1, V1t, ldv1t, work, work,
        lorglqwork, info, ioffv1t + 1 + ldv1t, ioffwork + itauq1 - 1,
        ioffwork + iorglq - 1 );
    }
    if ( wantv2t && m - q > 0 ) {
      LaPack0.dlacpy( 'U', p, m - q, X12, ldx12, V2t, ldv2t,
        ioffx12, ioffv2t );
      LaPack0.dlacpy( 'U', m - p - q, m - p - q, X22, ldx22,
        V2t, ldv2t, ioffx22 + q + p * ldx22, ioffv2t + p + p * ldv2t );
      LaPack3.dorglq( m - q, m - q, m - q, V2t, ldv2t, work, work,
        lorglqwork, info, ioffv2t, ioffwork + itauq2 - 1,
        ioffwork + iorglq - 1 );
    }
  } else {
    if ( wantu1 && p > 0 ) {
      LaPack0.dlacpy( 'U', q, p, X11, ldx11, U1, ldu1,
        ioffx11, ioffu1 );
      LaPack3.dorgqr( p, p, q, U1, ldu1, work, work, lorglqwork, info,
        ioffu1, ioffwork + itaup1 - 1, ioffwork + iorglq - 1 );
    }
    if ( wantu2 && m - p > 0 ) {
      LaPack0.dlacpy( 'U', q, m - p, X21, ldx21, U2, ldu2,
        ioffx21, ioffu2 );
      LaPack3.dorglq( m - p, m - p, q, U2, ldu2, work, work,
        lorglqwork, info, ioffu2, ioffwork + itaup2 - 1,
        ioffwork + iorglq - 1 );
    }
    if ( wantv1t && q > 0 ) {
      LaPack0.dlacpy( 'L', q - 1, q - 1, X11, ldx11, V1t, ldv1t,
        ioffx11 + 1, ioffv1t + 1 + ldv1t );
      V1t[ ioffv1t ] = 1.;
      for ( j = 2; j <= q; j ++ ) {
        V1t[ ioffv1t + ( j - 1 ) * ldv1t ] = 0.;
        V1t[ ioffv1t + j - 1 ] = 0.;
      }
      LaPack3.dorgqr( q - 1, q - 1, q - 1, V1t, ldv1t, work, work,
        lorgqrwork, info, ioffv1t + 1 + ldv1t, ioffwork + itauq1 - 1,
        ioffwork + iorgqr - 1 );
    }
    if ( wantv2t && m - q > 0 ) {
      LaPack0.dlacpy( 'L', m - q, p, X12, ldx12, V2t, ldv2t,
        ioffx12, ioffv2t );
      LaPack0.dlacpy( 'L', m - p - q, m - p - q, X22, ldx22,
        V2t, ldv2t, ioffx22 + p + q * ldx22, ioffv2t + p + p * ldv2t );
      LaPack3.dorgqr( m - q, m - q, m - q, V2t, ldv2t, work, work,
        lorglqwork, info, ioffv2t, ioffwork + itauq2 - 1,
        ioffwork + iorgqr - 1 );
    }
  }
  LaPack3.dbbcsd( jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta,
    work, U1, ldu1, U2, ldu2, V1t, ldv1t, V2t, ldv2t, work, work,
    work, work, work, work, work, work, work, lbbcsdwork, info,
    iofftheta, ioffwork + iphi - 1, ioffu1, ioffu2, ioffv1t, ioffv2t,
    ioffwork + ib11d - 1, ioffwork + ib11e - 1, ioffwork + ib12d - 1,
    ioffwork + ib12e - 1, ioffwork + ib21d - 1, ioffwork + ib21e - 1,
    ioffwork + ib22d - 1, ioffwork + ib22e - 1, ioffwork + ibbcsd - 1
    );
  if ( q > 0 && wantu2 ) {
    for ( var i = 1; i <= q; i ++ ) {
      iwork[ ioffwork + i - 1 ] = m - p - q + i;
    }
    for ( i = q + 1; i <= m - p; i ++ ) {
      iwork[ ioffwork + i - 1 ] = i - q;
    }
    if ( colmajor ) {
      LaPack0.dlapmt( false, m - p, m - p, U2, ldu2, iwork,
        ioffu2, ioffiwork );
    } else {
      LaPack0.dlapmr( false, m - p, m - p, U2, ldu2, iwork,
        ioffu2, ioffiwork );
    }
  }
  if ( m > 0 && wantv2t ) {
    for ( i = 1; i <= p; i ++ ) {
      iwork[ ioffiwork + i - 1 ] = m - p - q + i;
    }
    for ( i = p + 1; i <= m - q; i ++ ) {
      iwork[ ioffiwork + i - 1 ] = i - p;
    }
    if ( ! colmajor ) {
      LaPack0.dlapmt( false, m - q, m - q, V2t, ldv2t, iwork,
        ioffv2t, ioffiwork );
    } else {
      LaPack0.dlapmr( false, m - q, m - q, V2t, ldv2t, iwork,
        ioffv2t, ioffiwork );
    }
  }
}
//*************************************************************************
LaPack4.dorgbr = function( vect, m, n, k, A, lda, tau, work,
lwork, info, ioffa, iofftau, ioffwork ) {
  throw new Error("not tested");
  info.setValue( 0 );
  var wantq = ( vect.charAt(0).toUpperCase() == 'Q' );
  var mn = Math.min( m, n );
  var lquery = ( lwork == -1 );
  if ( ! wantq && vect.charAt(0).toUpperCase() != 'P' ) {
    info.setValue( -1 );
  } else if ( m < 0 ) info.setValue( -2 );
  else if ( n < 0 || ( wantq && ( n > m || n < Math.min( m, k ) ) ) ||
  ( ! wantq && ( m > n || m < Math.min( n, k ) ) ) ) {
    info.setValue( -3 );
  } else if ( k < 0 ) info.setValue( -4 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -6 );
  else if ( lwork < Math.max( 1, mn ) ) info.setValue( -9 );
  var iinfo = new IntReference();
  if ( info.getValue() == 0 ) {
    work[ ioffwork ] == 1;
    if ( wantq ) {
      if ( m >= k ) {
        LaPack3.dorgqr( m, n, k, A, lda, tau, work, -1, iinfo,
          ioffa, iofftau, ioffwork );
      } else if ( m > 1 ) {
        LaPack3.dorgqr( m - 1, m - 1, m - 1, A, lda, tau, work, -1,
          iinfo, ioffa + 1 + lda, iofftau, ioffwork );
      }
    } else {
      if ( k < n ) {
        LaPack3.dorglq( m, n, k, A, lda, tau, work, -1, iinfo,
          ioffa, iofftau, ioffwork );
      } else if ( n > 1 ) {
        LaPack3.dorglq( n - 1, n - 1, n - 1, A, lda, tau, work, -1,
          iinfo, ioffa + 1 + lda, iofftau, ioffwork );
      }
    }
    var lwkopt = work[ ioffwork ];
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dorbgr', - info.getValue() );
    return;
  } else if ( lquery ) return;
  if ( m == 0 || n == 0 ) {
    work[ ioffwork ] = 1;
    return;
  }
  if ( wantq ) {
    if ( m >= k ) {
      LaPack3.dorgqr( m, n, k, A, lda, tau, work, lwork, iinfo, ioffa,
        iofftau, ioffwork );
    } else {
      for ( var j = m; j >= 2; j -- ) {
        A[ ioffa + ( j - 1 ) * lda ] = 0.;
        for ( var i = j + 1; i <= m; i ++ ) {
          A[ ioffa + i - 1 + ( j - 1 ) ] =
            A[ ioffa + i - 1 + ( j - 2 ) * lda ];
        }
      }
      A[ ioffa ] = 1.;
      for ( i = 2; i <= m; i ++ ) A[ ioffa + i - 1 ] = 0.;
      if ( m > 1 ) {
        LaPack3.dorgqr( m - 1, m - 1, m - 1, A, lda, tau, work, lwork,
          iinfo, ioffa + 1 + lda, iofftau, ioffwork );
      }
    }
  } else {
    if ( k < n ) {
      LaPack3.dorglq( m, n, k, A, lda, tau, work, lwork, iinfo, ioffa,
        iofftau, ioffwork );
    } else {
      A[ ioffa ] = 1.;
      for ( i = 2; i <= n; i ++ ) A[ ioffa + i - 1 ] = 0.;
      for ( j = 2; j <= n; j ++ ) {
        for ( i = j - 1; i >= 2; i -- ) {
          A[ ioffa + i - 1 + ( j - 1 ) * lda ] =
            A[ ioffa + i - 2 + ( j - 1 ) * lda ];
        }
        A[ ioffa + ( j - 1 ) * lda ] = 0.;
      }
      if ( n > 1 ) {
        LaPack3.dorglq( n - 1, n - 1, n - 1, A, lda, tau, work, lwork,
          iinfo, ioffa + 1 + lda, iofftau, ioffwork );
      }
    }
  }
}
//*************************************************************************
LaPack4.dorghr = function( n, ilo, ihi, A, lda, tau, work, lwork,
info, ioffa, iofftau, ioffwork ) {
  throw new Error("not tested");
  info.setValue( 0 );
  var nh = ihi - ilo;
  var lquery = ( lwork == -1 );
  if ( n < 0 ) info.setValue( -1 );
  else if ( ilo < 1 || ilo > Math.max( 1, n ) ) info.setValue( -2 );
  else if ( ihi < Math.min( ilo, n ) || ihi > n ) info.setValue( -3 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -5 );
  else if ( lwork < Math.max( 1, nh ) && ! lquery ) info.setValue( -8 );
  if ( info.getValue() == 0 ) {
    var nb = LaPack0.ilaenv( 1, 'dorgqr', ' ', nh, nh, nh, -1 );
    var lwkopt = Math.max( 1, nh ) * nb;
    work[ ioffwork ] = lwkopt;
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dorghr', - info.getValue() );
    return;
  } else if ( lquery ) return;
  if ( n == 0 ) {
    work[ ioffwork ] = 1;
    return;
  }
  for ( var j = ihi; j >= ilo + 1; j -- ) {
    for ( var i = 1; i <= j - 1; i ++ ) {
      A[ ioffa + i - 1 + ( j - 1 ) * lda ] = 0.;
    }
    for ( i = j + 1; i <= ihi; i ++ ) {
      A[ ioffa + i - 1 + ( j - 1 ) * lda ] =
        A[ ioffa + i - 1 + ( j - 2 ) * lda ];
    }
    for ( i = ihi + 1; i <= n; i ++ ) {
      A[ ioffa + i - 1 + ( j - 1 ) * lda ] = 0.;
    }
  }
  for ( j = 1; j <= ilo; j ++ ) {
    for ( i = 1; i <= n; i ++ ) {
      A[ ioffa + i - 1 + ( j - 1 ) * lda ] = 0.;
    }
    A[ ioffa + j - 1 + ( j - 1 ) * lda ] = 1.;
  }
  for ( j = ihi + 1; j <= n; j ++ ) {
    for ( i = 1; i <= n; i ++ ) {
      A[ ioffa + i - 1 + ( j - 1 ) * lda ] = 0.;
    }
    A[ ioffa + j - 1 + ( j - 1 ) * lda ] = 1.;
  }
  if ( nh > 0 ) {
    var iinfo = new IntReference();
    LaPack3.dorgqr( nh, nh, nh, A, lda, tau, work, lwork, iinfo,
      ioffa + ilo + ilo * lda, iofftau + ilo, ioffwork );
  }
  work[ ioffwork ] = lwkopt;
}
//*************************************************************************
LaPack4.dorgtr = function( uplo, n, A, lda, tau, work, lwork,
info, ioffa, iofftau, ioffwork ) {
  throw new Error("not tested");
  info.getValue() = 0;
  var lquery = ( lwork == -1 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -4 );
  else if ( lwork < Math.max( 1, n - 1 ) && ! lquery ) info.setValue( -7 );
  if ( info.getValue() == 0 ) {
    var nb = ( upper ?
      LaPack0.ilaenv( 1, 'dorgql', ' ', n - 1, n - 1, n - 1, -1 ) :
      LaPack0.ilaenv( 1, 'dorgqr', ' ', n - 1, n - 1, n - 1, -1 ) );
    var lwkopt = Math.max( 1, n - 1 ) * nb;
    work[ ioffwork ] = lwkopt;
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dorgtr', - info.getValue() );
    return;
  } else if ( lquery ) return;
  if ( n == 0 ) {
    work[ ioffwork ] = 1;
    return;
  }
  var iinfo = new IntReference();
  if ( upper ) {
    for ( var j = 1; j <= n - 1; j ++ ) {
      for ( var i = 1; i <= j - 1; i ++ ) {
        A[ ioffa + i - 1 + ( j - 1 ) * lda ] =
          A[ ioffa + i - 1 + j * lda ];
      }
      A[ ioffa + n - 1 + ( j - 1 ) * lda ] = 0.;
    }
    for ( i = 1; i <= n - 1; i ++ ) {
      A[ ioffa + i - 1 + ( n - 1 ) * lda ] = 0.;
    }
    A[ ioffa + n - 1 + ( n - 1 ) * lda ] = 1.;
    LaPack3.dorgql( n - 1, n - 1, n - 1, A, lda, tau, work, lwork,
      iinfo, ioffa, iofftau, ioffwork );
  } else {
    for ( j = n; j >= 2; j -- ) {
      A[ ioffa + ( j - 1 ) * lda ] = 0.;
      for ( i = j + 1; i <= n; i ++ ) {
        A[ ioffa + i - 1 + ( j - 1 ) * lda ] =
          A[ ioffa + i - 1 + ( j - 2 ) * lda ];
      }
    }
    A[ ioffa ] = 1.;
    for ( i = 2; i <= n; i ++ ) A[ ioffa + i - 1 ] = 0.;
    if ( n > 1 ) {
      LaPack3.dorgqr( n - 1, n - 1, n - 1, A, lda, tau, work, lwork,
        iinfo, ioffa + 1 + lda, iofftau, ioffwork );
    }
  }
  work[ ioffwork ] = lwkopt;
}
//*************************************************************************
LaPack4.dormbr = function( vect, side, trans, m, n, k, A, lda,
tau, C, ldc, work, lwork, info, ioffa, iofftau, ioffc, ioffwork ) {
  throw new Error("not tested");
  info.setValue( 0 );
  var applyq = ( vect.charAt(0).toUpperCase() == 'Q' );
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
  if ( ! applyq && vect.charAt(0).toUpperCase() != 'P' ) {
    info.setValue( -1 );
  } else if ( ! left && side.charAt(0).toUpperCase() != 'R' ) {
    info.setValue( -2 );
  } else if ( ! notran && trans.charAt(0).toUpperCase() != 'T' ) {
    info.setValue( -3 );
  } else if ( m < 0 ) info.setValue( -4 );
  else if ( n < 0 ) info.setValue( -5 );
  else if ( k < 0 ) info.setValue( -6 );
  else if ( ( applyq && lda < Math.max( 1, nq ) ) ||
  ( ! applyq && lda < Math.max( 1, Math.min( nq, k ) ) ) ) {
    info.setValue( -8 );
  } else if ( ldc < Math.max( 1, m ) ) info.setValue( -11 );
  else if ( lwork < Math.max( 1, nw ) && ! lquery ) info.setValue( -13 );
  if ( info.getValue() == 0 ) {
    if ( applyq ) {
      var nb = ( left ?
        LaPack0.ilaenv(1, 'dormqr', side + trans, m - 1, n, m-1,-1) :
        LaPack0.ilaenv(1, 'dormqr', side + trans, m, n - 1, n-1,-1) );
    } else {
      nb = ( left ?
        LaPack0.ilaenv(1, 'dormlq', side + trans, m - 1, n, m-1,-1) :
        LaPack0.ilaenv(1, 'dormlq', side + trans, m, n - 1, n-1,-1) );
    }
    var lwkopt = Math.max( 1, nw ) * nb;
    work[ ioffwork ] = lwkopt;
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dormbr', - info.getValue() );
    return;
  } else if ( lquery ) return;
  work[ ioffwork ] = 1;
  if ( m == 0 || n == 0 ) return;
  var iinfo = new IntReference();
  if ( applyq ) {
    if ( nq >= k ) {
      LaPack3.dormqr( side, trans, m, n, k, A, lda, tau, C, ldc, work,
        lwork, iinfo, ioffa, iofftau, ioffc, ioffwork );
    } else if ( nq > 1 ) {
      if ( left ) {
        var mi = m - 1;
        var ni = n;
        var i1 = 2;
        var i2 = 1;
      } else {
        mi = m;
        ni = n - 1;
        i1 = 1;
        i2 = 2;
      }
      LaPack3.dormqr( side, trans, mi, ni, nq - 1, A, lda, tau, C, ldc,
        work, lwork, iinfo, ioffa + 1, iofftau,
        ioffc + i1 - 1 + ( i2 - 1 ) * ldc, ioffwork );
    }
  } else {
    var transt = ( notran ? 'T' : 'N' );
    if ( nq > k ) {
      LaPack3.dormlq( side, transt, m, n, k, A, lda, tau, C, ldc, work,
        lwork, iinfo, ioffa, iofftau, ioffc, ioffwork );
    } else if ( nq > 1 ) {
      if ( left ) {
        mi = m - 1;
        ni = n;
        i1 = 2;
        i2 = 1;
      } else {
        mi = m;
        ni = n - 1;
        i1 = 1;
        i2 = 2;
      }
      LaPack3.dormlq( side, transt, mi, ni, nq - 1, A, lda, tau, C,
        ldc, work, lwork, iinfo, ioffa + lda, iofftau,
        ioffc + i1 - 1 + ( i2 - 1 ) * ldc, ioffwork );
    }
  }
  work[ ioffwork ] = lwkopt;
}
//*************************************************************************
LaPack4.dormhr = function( side, trans, m, n, ilo, ihi, A, lda,
tau, C, ldc, work, lwork, info, ioffa, iofftau, ioffc, ioffwork ) {
  throw new Error("not tested");
  info.setValue( 0 );
  var nh = ihi - ilo;
  var left = ( side.charAt(0).toUpperCase() == 'L' );
  var lquery = ( lwork == -1 );
  if ( left) {
    var nq = m;
    var nw = n;
  } else {
    nq = n;
    nw = m;
  }
  if ( ! left && side.charAt(0).toUpperCase() != 'R' ) info.setValue( -1 );
  else if ( trans.charAt(0).toUpperCase() != 'N' &&
  trans.charAt(0).toUpperCase() != 'T' ) {
    info.setValue( -2 );
  } else if ( m < 0 ) info.setValue( -3 );
  else if ( n < 0 ) info.setValue( -4 );
  else if ( ilo < 1 || ilo > Math.max( 1, nq ) ) info.setValue( -5 );
  else if ( ihi < Math.min( ilo, nq ) || ihi > nq ) info.setValue( -6 );
  else if ( lda < Math.max( 1, nq ) ) info.setValue( -8 );
  else if ( ldc < Math.max( 1, m ) ) info.setValue( -11 );
  else if ( lwork < Math.max( 1, nw ) && ! lquery ) info.setValue( -13 );
  if ( info.getValue() == 0 ) {
    var nb = ( left ?
      LaPack0.ilaenv( 1, 'dormqr', side + trans, nh, n, nh, -1 ) :
      LaPack0.ilaenv( 1, 'dormqr', side + trans, m, nh, nh, -1 ) );
    var lwkopt = Math.max( 1, nw ) * nb;
    work[ ioffwork ] = lwkopt;
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dormhr', - info.getValue() );
    return;
  } else if ( lquery ) return;
  if ( m == 0 || n == 0 || nh == 0 ) {
    work[ ioffwork ] = 1;
    return;
  }
  if ( left ) {
    var mi = nh;
    var ni = n;
    var i1 = ilo + 1;
    var i2 = 1;
  } else {
    mi = m;
    ni = nh;
    i1 = 1;
    i2 = ilo + 1;
  }
  var iinfo = new IntReference();
  LaPack3.dormqr( side, trans, mi, ni, nh, A, lda, tau, C, ldc, work,
    lwork, iinfo, ioffa + ilo + ( ilo - 1 ) * lda, iofftau + ilo - 1,
    ioffc + i1 - 1 + ( i2 - 1 ) * ldc, ioffwork );
  work[ ioffwork ] = lwkopt;
}
//*************************************************************************
LaPack4.dormtr = function( side, uplo, trans, m, n, A, lda, tau,
C, ldc, work, lwork, info, ioffa, iofftau, ioffc, ioffwork ) {
  throw new Error("not tested");
  info.setValue( 0 );
  var left = ( side.charAt(0).toUpperCase() == 'L' );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var lquery = ( lwork == -1 );
  if ( left ) {
    var nq = m;
    var nw = n;
  } else {
    nq = n;
    nw = m;
  }
  if ( ! left && side.charAt(0).toUpperCase() != 'R' ) {
    info.setValue( -1 );
  } else if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -2 );
  } else if ( trans.charAt(0).toUpperCase() != 'N' &&
  trans.charAt(0).toUpperCase() != 'T' ) {
    info.setValue( -3 );
  } else if ( m < 0 ) info.setValue( -4 );
  else if ( n < 0 ) info.setValue( -5 );
  else if ( lda < Math.max( 1, nq ) ) info.setValue( -7 );
  else if ( ldc < Math.max( 1, m ) ) info.setValue( -10 );
  else if ( lwork < Math.max( 1, nw ) && ! lquery ) info.setValue( -12 );
  if ( info.getValue() == 0 ) {
    if ( upper ) {
      var nb = ( left ?
        LaPack0.ilaenv(1, 'dormql', side + trans, m - 1, n, m-1,-1) :
        LaPack0.ilaenv(1, 'dormql', side + trans, m, n - 1, n-1,-1) );
    } else {
      nb = ( left ?
        LaPack0.ilaenv(1, 'dormqr', side + trans, m - 1, n, m-1,-1) :
        LaPack0.ilaenv(1, 'dormqr', side + trans, m, n - 1, n-1,-1) );
    }
    var lwkopt = Math.max( 1, nw ) * nb;
    work[ ioffwork ] = lwkopt;
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dormtr', - info.getValue() );
    return;
  } else if ( lquery ) return;
  if ( m == 0 || n == 0 || nq == 1 ) {
    work[ ioffwork ] = 1;
    return;
  }
  if ( left ) {
    var mi = m - 1;
    var ni = n;
  } else {
    mi = m;
    ni = n - 1;
  }
  var iinfo = new IntReference();
  if ( upper ) {
    LaPack3.dormql( side, trans, mi, ni, nq - 1, A, lda, tau, C, ldc,
      work, lwork, iinfo, ioffa + lda, iofftau, ioffc, ioffwork );
  } else {
    if ( left ) {
      var i1 = 2;
      var i2 = 1;
    } else {
      i1 = 1;
      i2 = 2;
    }
    LaPack3.dormqr( side, trans, mi, ni, nq - 1, A, lda, tau, C, ldc,
      work, lwork, iinfo, ioffa + 1, iofftau,
      ioffc + i1 - 1 + ( i2 - 1 ) * ldc, ioffwork );
  }
  work[ ioffwork ] = lwkopt;
}
//*************************************************************************
LaPack4.dspev = function( jobz, uplo, n, AP, w, Z, ldz, work,
info) {
  throw new Error("not programmed: packed matrix");
}
//*************************************************************************
LaPack4.dspevx = function( jobz, range, uplo, n, AP, vl, vu, il,
iu, abstol, m, w, Z, ldz, work, iwork, ifail, info, ioffap, ioffw, ioffz,
ioffwork, ioffiwork, ioffifail ) {
  throw new Error("not programmed: packed matrix");
}
//*************************************************************************
LaPack4.dtgsy2 = function( trans, ijob, m, n, A, lda, B, ldb, C,
ldc, D, ldd, E, lde, F, ldf, scaleReference, rdsumReference,
rdscaleReference, iwork, pq, info, ioffa, ioffb, ioffc, ioffd, ioffe,
iofff, ioffiwork ) {
  throw new Error("not programmed: generalized eigengetValue()");
}
//*************************************************************************
LaPack4.dtrexc = function( compq, n, T, ldt, Q, ldq, ifst, ilst,
work, info, iofft, ioffq, ioffwork ) {
  throw new Error("not tested");
  info.setValue( 0 );
  var wantq = ( compq.charAt(0).toUpperCase() == 'V' );
  if ( ! wantq && compq.charAt(0).toUpperCase() != 'N' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( ldt < Math.max( 1, n ) ) info.setValue( -4 );
  else if ( ldq < 1 || ( wantq && ldq < Math.max( 1, n ) ) ) {
    info.setValue( -6 );
  } else if ( ifst.getValue() < 1 || ifst.getValue() > n ) info.setValue( -7 );
  else if ( ilst.getValue() < 1 || ilst.getValue() > n ) info.setValue( -8 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dtrexc', - info.getValue() );
    return;
  }
  if ( n <= 1 ) return;
  if ( ifst.getValue() > 1 ) {
    if ( T[ iofft + ifst.getValue() - 1 + ( ifst.getValue() - 2 ) * ldt ] != 0. )
    {
      ifst.setValue( ifst.getValue() - 1 );
    }
  }
  var nbf = 1;
  if ( ifst.getValue() < n ) {
    if ( T[ iofft + ifst.getValue() + ( ifst.getValue() - 1 ) * ldt ] != 0. ) {
      nbf = 2;
    }
  }
  if ( ilst.getValue() > 1 ) {
    if ( T[ iofft + ilst.getValue() - 1 + ( ilst.getValue() - 2 ) * ldt ] != 0. )
    {
      ilst.setValue( ilst.getValue() - 1 );
    }
  }
  var nbl = 1;
  if ( ilst.getValue() < n ) {
    if ( T[ iofft + ilst.getValue() + ( ilst.getValue() - 1 ) * ldt ] != 0. ) {
      nbl = 2;
    }
  }
  if ( ifst.getValue() == ilst.getValue() ) return;
  if ( ifst.getValue() < ilst.getValue() ) {
    if ( nbf == 2 && nbl == 1 ) ilst.getValue() = ilst.getValue() - 1;
    if ( nbf == 1 && nbl == 2 ) ilst.getValue() = ilst.getValue() + 1;
    var here = ifst.getValue();
    do {
      if ( nbf == 1 || nbf == 2 ) {
        var nbnext = 1;
        if ( here + nbf + 1 <= n ) {
          if ( T[ iofft + here + nbf + ( here + nbf - 1 ) * ldt ] !=
          0. ) {
            nbnext = 2;
          }
        }
        LaPack3.dlaexc( wantq, n, T, ldt, Q, ldq, here, nbf, nbnext,
          work, info, iofft, ioffq, ioffwork );
        if ( info.getValue() != 0 ) {
          ilst.setValue( here );
          return;
        }
        here += nbnext;
        if ( nbf == 2 ) {
          if ( T[ iofft + here + ( here - 1 ) * ldt ] == 0. ) nbf = 3;
        }
      } else {
        nbnext = 1;
        if ( here + 3 <= n ) {
          if ( T[ iofft + here + 2 + ( here + 1 ) * ldt ] != 0. ) {
            nbnext = 2;
          }
        }
        LaPack3.dlaexc( wantq, n, T, ldt, Q, ldq, here + 1, 1, nbnext,
          work, info, iofft, ioffq, ioffwork );
        if ( info.getValue() != 0 ) {
          ilst.setValue( here );
          return;
        }
        if ( nbnext == 1 ) {
          LaPack3.dlaexc( wantq, n, T, ldt, Q, ldq, here, 1, nbnext,
            work, info, iofft, ioffq, ioffwork );
          here ++;
        } else {
          if ( T[ iofft + here + 1 + here * ldt ] == 0. ) nbnext = 1;
          if ( nbnext == 2 ) {
            LaPack3.dlaexc( wantq, n, T, ldt, Q, ldq, here, 1, nbnext,
              work, info, iofft, ioffq, ioffwork );
            if ( info.getValue() != 0 ) {
              ilst.setValue( here );
              return;
            }
            here += 2;
          } else {
            LaPack3.dlaexc( wantq, n, T, ldt, Q, ldq, here, 1, 1,
              work, info, iofft, ioffq, ioffwork );
            LaPack3.dlaexc( wantq, n, T, ldt, Q, ldq, here + 1, 1, 1,
              work, info, iofft, ioffq, ioffwork );
            here += 2;
          }
        }
      }
    } while ( here < ilst.getValue() );
  } else {
    here = ifst.getValue();
    do {
      if ( nbf == 1 || nbf == 2 ) {
        nbnext= 1;
        if ( here >= 3 ) {
          if ( T[ iofft + here - 2 + ( here - 3 ) * ldt ] != 0. ) {
            nbnext = 2;
          }
        }
        LaPack3.dlaexc( wantq, n, T, ldt, Q, ldq, here - nbnext,
          nbnext, nbf, work, info, iofft, ioffq, ioffwork );
        if ( info.getValue() != 0 ) {
          ilst.setValue( here );
          return;
        }
        here -= nbnext;
        if ( nbf == 2 ) {
          if ( T[ iofft + here + ( here - 1 ) * ldt ] == 0. ) nbf = 3;
        }
      } else {
        nbnext = 1;
        if ( here >= 3 ) {
          if ( T[ iofft + here - 2 + ( here - 3 ) * ldt ] != 0. ) {
            nbnext = 2;
          }
        }
        LaPack3.dlaexc( wantq, n, T, ldt, Q, ldq, here - nbnext, 
          nbnext, 1, work, info, iofft, ioffq, ioffwork );
        if ( info.getValue() != 0 ) {
          ilst.setValue( here );
          return;
        }
        if ( nbnext == 1 ) {
          LaPack3.dlaexc( wantq, n, T, ldt, Q, ldq, here, nbnext, 1,
            work, info, iofft, ioffq, ioffwork );
          here --;
        } else {
          if ( T[ iofft + here - 1 + ( here - 2 ) * ldt ] == 0. ) {
            nbnext = 1;
          }
          if ( nbnext == 2 ) {
            LaPack3.dlaexc( wantq, n, T, ldt, Q, ldq, here - 1, 2, 1,
              work, info, iofft, ioffq, ioffwork );
            if ( info.getValue() != 0 ) {
              ilst.setValue( here );
              return;
            }
            here -= 2;
          } else {
            LaPack3.dlaexc( wantq, n, T, ldt, Q, ldq, here, 1, 1,
              work, info, iofft, ioffq, ioffwork );
            LaPack3.dlaexc( wantq, n, T, ldt, Q, ldq, here - 1, 1, 1,
              work, info, iofft, ioffq, ioffwork );
            here -= 2;
          }
        }
      }
    } while ( here > ilst.getValue() );
  }
  ilst.setValue( here );
}
LaPack4.ztrexc = function( compq, n, T, ldt, Q, ldq, ifst, ilst
, work, info ) {
  throw new Error("not programmed: complex matrix");
}
