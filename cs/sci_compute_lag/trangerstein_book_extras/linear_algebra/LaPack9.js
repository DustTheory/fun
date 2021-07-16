function LaPack9() {
}
//**************************************************************************
LaPack9.dhseqr = function( job, compz, n, ilo, ihi, H, ldh, wr,
wi, Z, ldz, work, lwork, info, ioffh, ioffwr, ioffwi, ioffz, ioffwork ) {
  var ntiny = 11;
  var nl = 49;
  var Hl = new Array( nl * nl );
  var workl = new Array ( nl );
  var wantt = ( job.charAt(0).toUpperCase() == 'S' );
  var initz = ( compz.charAt(0).toUpperCase() == 'I' );
  var wantz =
    initz || ( compz.charAt(0).toUpperCase() == 'V' );
  work[ ioffwork ] = Math.max( 1, n );
  var lquery = ( lwork == -1 );
  info.setValue( 0 );
  if ( job.charAt(0).toUpperCase() != 'E' && ! wantt ) {
    info.setValue( -1 );
  } else if ( compz.charAt(0).toUpperCase() != 'N' && ! wantz ) {
    info.setValue( -2 );
  } else if ( n < 0 ) info.setValue( -3 );
  else if ( ilo < 1 || ilo > Math.max( 1, n ) ) info.setValue( -4 );
  else if ( ihi < Math.min( ilo, n ) || ihi > n ) info.setValue( -5 );
  else if ( ldh < Math.max( 1, n ) ) info.setValue( -7 );
  else if ( ldz < 1 || ( wantz && ldz < Math.max( 1, n ) ) ) {
    info.setValue( -11 );
  } else if ( lwork < Math.max( 1, n ) && ! lquery ) info.setValue( -13 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dhseqr', - info.getValue() );
    return;
  } else if ( n == 0 ) return;
  else if ( lquery ) {
    LaPack8.dlaqr0( wantt, wantz, n, ilo, ihi, H, ldh, wr, wi, ilo,
      ihi, Z, ldz, work, lwork, info, ioffh, ioffwr, ioffwi, ioffz,
      ioffwork);
    work[ ioffwork ] = Math.max( Math.max( 1, n ), work[ ioffwork ] );
    return;
  } else {
    for ( var i = 1; i <= ilo - 1; i ++ ) {
      wr[ ioffwr + i - 1 ] = H[ ioffh + i - 1 + ( i - 1 ) * ldh ];
      wi[ ioffwi + i - 1 ] = 0.;
    }
    for ( i = ihi + 1; i <= n; i ++ ) {
      wr[ ioffwr + i - 1 ] = H[ ioffh + i - 1 + ( i - 1 ) * ldh ];
      wi[ ioffwi + i - 1 ] = 0.;
    }
    if ( initz ) {
      LaPack0.dlaset( 'A', n, n, 0., 1., Z, ldz, ioffz );
    }
    if ( ilo == ihi ) {
      wr[ ioffwr + ilo - 1 ] =
        H[ ioffh + ilo - 1 + ( ilo - 1 ) * ldh ];
      wi[ ioffwi + ilo - 1 ] = 0.;
      return;
    }
    var nmin =
      LaPack0.ilaenv( 12, 'dhseqr', job.charAt(0) + compz.charAt(0),
        n, ilo, ihi, lwork );
    nmin = Math.max( ntiny, nmin );
    if ( n > nmin ) {
      LaPack8.dlaqr0( wantt, wantz, n, ilo, ihi, H, ldh, wr, wi, ilo,
        ihi, Z, ldz, work, lwork, info, ioffh, ioffwr, ioffwi, ioffz,
        ioffwork );
    } else {
      LaPack2.dlahqr( wantt, wantz, n, ilo, ihi, H, ldh, wr, wi, ilo,
        ihi, Z, ldz, info, ioffh, ioffwr, ioffwi, ioffz );
      if ( info.getValue() > 0 ) {
        var kbot = info.getValue();
        if ( n >= nl ) {
          LaPack8.dlaqr0( wantt, wantz, n, ilo, kbot, H, ldh, wr, wi,
            ilo, ihi, Z, ldz, work, lwork, info, ioffh, ioffwr, ioffwi,
            ioffz, ioffwork );
        } else {
          LaPack0.dlacpy( 'A', n, n, H, ldh, Hl, nl, ioffh, 0 );
          Hl[ n + ( n - 1 ) * nl ] = 0.;
          LaPack0.dlaset( 'A', nl, nl - n, 0., 0., Hl, nl, n * nl );
          LaPack8.dlaqr0( wantt, wantz, nl, ilo, kbot, Hl, nl, wr, wi,
            ilo, ihi, Z, ldz, workl, nl, info, 0, ioffwr, ioffwi,
            ioffz, 0 );
          if ( wantt || info.getValue() != 0 ) {
            LaPack0.dlacpy( 'A', n, n, Hl, nl, H, ldh, 0, ioffh );
          }
        }
      }
    }
    if ( ( wantt || info.getValue() != 0 ) && n > 2 ) {
      LaPack0.dlaset( 'L', n - 2, n - 2, 0., 0., H, ldh, 2 );
    }
    work[ ioffwork ] = Math.max( Math.max( 1, n ), work[ ioffwork ] );
  }
}
LaPack9.zhseqr = function( job, compz, n, ilo, ihi, H, ldh, w, Z,
ldz, work, lwork, info ) {
  throw new Error("not programmed: complex matrix");
}
