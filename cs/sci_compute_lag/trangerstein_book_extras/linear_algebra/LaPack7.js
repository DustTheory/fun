function LaPack7() {
}
//**************************************************************************
LaPack7.dlaqr3 = function( wantt, wantz, n, ktop, kbot, nw, H,
ldh, iloz, ihiz, Z, ldz, ns, nd, sr, si, V, ldv, nh, T, ldt, nv, Wv, ldwv,
work, lwork, ioffh, ioffz, ioffsr, ioffsi, ioffv, iofft, ioffwv,
ioffwork) {
  throw new Error("not tested");
  var jw = Math.nin( nw, kbot - ktop + 1);
  if ( jw <= 2 ) var lwkopt = 1;
  else {
    LaPack3.dgehrd( jw, 1, jw - 1, T, ldt, work, work, -1, info,
      iofft );
    var lwk1 = Math.round( work[ ioffwork ] );
    LaPack4.dormhr( 'R', 'N', jw, jw, 1, jw - 1, T, ldt, work, V, ldv,
      work, -1, info, iofft, ioffwork, ioffv, ioffwork );
    var lwk2 = Math.round( work[ ioffwork ] );
    var infqr = new IntReference();
    LaPack6.dlaqr4( true, true, jw, 1, jw, T, ldt, sr, si, 1, jw,
      V, ldv, work, -1, infqr, iofft, ioffsr, ioffsi, ioffv,
      ioffwork );
    var lwk3 = Math.round( work[ ioffwork ] );
    lwkopt = Math.max( jw + Math.max( lwk1, lwk2 ), lwk3 );
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
  var safmin = 
    new NumberReference( LaPack0.dlamch( 'Safe minimum' ) );
  var safmax = 
    new NumberReference( 1. / safmin );
  LaPack0.dlabad( safmin, safmax );
  var ulp = LaPack0.dlamch( 'Precision' );
  var smlnum = safmin.getValue() * ( Number( n ) / ulp );
  jw = Math.min( nw, kbot - ktop + 1 );
  var kwtop = kbot - jw + 1;
  var s = ( kwtop == ktop ? 0. :
    H[ ioffh + kwtop - 1 + ( kwtop - 2 ) * ldh ] );
  if ( kbot == kwtop ) {
    sr[ ioffsr + kwtop - 1 ] =
      H[ ioffh + kwtop - 1 + ( kwtop - 1 ) * ldh ];
    si[ ioffsi + kwtop - 1 ] = 0.;
    ns.setValue( 1 );
    nd.setValue( 0 );
    if ( Math.abs( s ) <= Math.max( smlnum,
    ulp * Math.abs( H[ ioffh + kwtop - 1 + ( kwtop - 1 ) * ldh ] ) ) ) {
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
  var nmin = LaPack0.ilaenv( 12, 'dlaqr3', 'sv', jw, 1, jw, lwork );
  if ( jw > nmin ) {
    LaPack6.dlaqr4( true, true, jw, 1, jw, T, ldt, sr, si, 1, jw,
      V, ldv, work, lwork, infqr, iofft, ioffsr + kwtop - 1,
      ioffsi + kwtop - 1, ioffv, ioffwork );
  } else {
    LaPack2.dlahqr( true, true, jw, 1, jw, T, ldt, sr, si, 1, jw,
      V, ldv, infqr, iofft, ioffsr + kwtop - 1, ioffsi + kwtop - 1,
      ioffv );
  }
  for ( var j = 1; j <= jw - 3; j ++ ) {
    T[ iofft + j + 1 + ( j - 1 ) * ldt ] = 0.;
    T[ iofft + j + 2 + ( j - 1 ) * ldt ] = 0.;
  } // 10
  if ( jw > 2 ) T[ iofft + jw - 1 + ( jw - 3 ) * ldt ] = 0.;
  ns.setValue( jw );
  var ilst = new IntReference( infqr.getValue() + 1 );
  var info = new IntReference();
  while ( ilst.getValue() <= ns.getValue() ) { // 20
    var bulge = ( ns.getValue() == 1 ? false :
      T[ iofft + ns.getValue() - 1 + ( ns.getValue() - 2 ) * ldt ] != 0. );
    if( ! bulge ) {
      var foo =
        Math.abs( T[ iofft + ns.getValue() - 1 + ( ns.getValue() - 1 ) * ldt ] );
      if ( foo == 0. ) foo = Math.abs( s );
      if ( Math.abs( s * V[ ioffv + ( ns.getValue() - 1 ) * ldv ] )
      <= Math.max( smlnum, ulp * foo ) ) {
        ns.setValue( ns.getValue() - 1 );
      } else {
        var ifst = 
          new IntReference( ns.getValue() );
        LaPack4.dtrexc( 'V', jw, T, ldt, V, ldv, ifst, ilst, work, info,
          iofft, ioffv, ioffwork );
        ilst.setValue( ilst.getValue() + 1 );
      }
    } else {
      foo =
        Math.abs( T[ iofft + ns.getValue() - 1 + ( ns.getValue() - 1 ) * ldt ] )
        + Math.sqrt( Math.abs(
          T[ iofft + ns.getValue() - 1 + ( ns.getValue() - 2 ) * ldt ] ))
        * Math.sqrt( Math.abs(
          T[ iofft + ns.getValue() - 2 + ( ns.getValue() - 1 ) * ldt ]));
      if ( foo == 0. ) foo = Math.abs( s );
      if ( Math.max(
      Math.abs( s * V[ ioffv + ( ns.getValue() - 1 ) * ldv ] ),
      Math.abs( s * V[ ioffv + ( ns.getValue() - 2 ) * ldv ] ) ) <=
      Math.max( smlnum, ulp * foo ) ) {
        ns.setValue( ns.getValue() - 2 );
      } else {
        ifst.setValue( ns.getValue() );
        LaPack4.dtrexc( 'V', jw, T, ldt, V, ldv, ifst, ilst, work, info,
          iofft, ioffv, ioffwork );
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
      else k = i + 2;
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
            + Math.sqrt( Math.abs( T[ iofft + k + ( k - 1 ) * ldt ] ) )
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
    } // 50
  }
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
      var aa =
        new NumberReference( T[ iofft + i - 2 + ( i - 2 ) * ldt ] );
      var cc =
        new NumberReference( T[ iofft + i - 1 + ( i - 2 ) * ldt ] );
      var bb =
        new NumberReference( T[ iofft + i - 2 + ( i - 1 ) * ldt ] );
      var dd =
        new NumberReference( T[ iofft + i - 1 + ( i - 1 ) * ldt ] );
      var rt1r =
        new NumberReference( sr[ ioffsr + kwtop + i - 3 ] );
      var rt1i =
        new NumberReference( si[ ioffsi + kwtop + i - 3 ] );
      var rt2r =
        new NumberReference( sr[ ioffsr + kwtop + i - 2 ] );
      var rt2i =
        new NumberReference( si[ ioffsi + kwtop + i - 2 ] );
      var cs = new NumberReference();
      var sn = new NumberReference();
      LaPack1.dlanv2( aa, bb, cc, dd, rt1r, rt1i, rt2r, rt2i, cs, sn );
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
      var beta =
        new NumberReference( work[ ioffwork ] );
      var tau = new NumberReference();
      LaPack1.dlarfg( ns.getValue(), beta, work, 1, tau, ioffwork + 1 );
      work[ ioffwork ] = 1.;
      LaPack0.dlaset( 'L', jw - 2, jw - 2, 0., 0., T, ldt, iofft + 2 );
      LaPack1.dlarf( 'L', ns.getValue(), jw, work, 1, tau.getValue(), T, ldt,
        work, ioffwork, iofft, ioffwork + jw );
      LaPack1.dlarf( 'R', ns.getValue(), ns.getValue(), work, 1, tau.getValue(),
        T, ldt, work, ioffwork, iofft, ioffwork + jw );
      LaPack1.dlarf( 'R', ns.getValue(), ns.getValue(), work, 1, tau.getValue(),
        V, ldv, work, ioffwork, ioffv, ioffwork + jw );
      LaPack3.dgehrd( jw, 1, ns.getValue(), T, ldt, work, work, lwork - jw,
        info, iofft, ioffwork, ioffwork + jw );
    }
    if ( kwtop > 1 ) {
      H[ ioffh + kwtop - 1 + ( kwtop - 2 ) * ldh ] = s * V[ ioffv ];
    }
    LaPack0.dlacpy( 'U', jw, jw, T, ldt, H, ldh, iofft,
      ioffh + kwtop - 1 + ( kwtop - 1 ) * ldh );
    Blas1.dcopy( jw - 1, T, ldt + 1, H, ldh + 1, iofft + 1,
      ioffh + kwtop + ( kwtop - 1 ) * ldh );
    if ( ns.getValue() > 1 && s != 0. ) {
      LaPack4.dormhr( 'R', 'N', jw, ns.getValue(), 1, ns.getValue(), T, ldt,
        work, V, ldv, work, lwork - jw, info, iofft, ioffwork, ioffv,
        ioffwork + jw );
    }
    var ltop = ( wantt ? 1 : ktop );
    for ( var krow = ltop; krow <= kwtop - 1; krow += nv ) {
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
          Wv, ldwv, ioffz + krow - 1 + ( kwtop - 1 ) * ldz, ioffv,
          ioffwv );
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
LaPack7.dlasd0 = function( n, sqre, d, e, U, ldu, Vt, ldvt,
smlsiz, iwork, work, info, ioffd, ioffe, ioffu, ioffvt, ioffiwork,
ioffwork ) {
  throw new Error("not tested");
  info.setValue( 0 );
  if ( n < 0 ) info.setValue( -1 );
  else if ( sqre < 0 || sqre > 1 ) info.setValue( -2 );
  var m = n + sqre;
  if ( ldu < n ) info.setValue( -6 );
  else if ( ldvt < m ) info.setValue( -8 );
  else if ( smlsiz < 3 ) info.setValue( -9 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlasd0', - info.getValue() );
    return;
  }
  if ( n <= smlsiz ) {
    LaPack6.dlasdq( 'U', sqre, n, m, n, 0, d, e, Vt, ldvt, U, ldu,
      work, info, ioffd, ioffe, ioffvt, ioffu, ioffwork );
    return;
  }
  var inode = 1;
  var ndiml = inode + n;
  var ndimr = ndiml + n;
  var idxq = ndimr + n;
  var iwk = idxq + n;
  var nlvl = new IntReference();
  var nd = new IntReference();
  LaPack0.dlasdt( n, nlvl, nd, iwork, iwork, iwork, smlsiz,
    ioffiwork + inode - 1, ioffiwork + ndiml - 1,
    ioffiwork + ndimr - 1 )
  var ndb1 = ( nd.getValue() + 1 ) / 2;
  var ncc = 0;
  for ( var i = ndb1; i <= nd.getValue(); i ++ ) {
    var i1 = i - 1;
    var ic = iwork[ ioffiwork + inode + i1 - 1 ];
    var nl = iwork[ ioffiwork + ndiml + i1 - 1 ];
    var nlp1 = nl + 1;
    var nr = iwork[ ioffiwork + ndimr + i1 - 1 ];
    var nrp1 = nr + 1;
    var nlf = ic - nl;
    var nrf = ic + 1;
    var sqrei = 1;
    LaPack6.dlasdq( 'U', sqrei, nl, nlp1, nl, ncc, d, e, Vt, ldvt,
      U, ldu, U, ldu, work, info, ioffd + nlf - 1, ioffe + nlf - 1,
      ioffvt + nlf - 1 + ( nlf - 1 ) * ldvt,
      ioffu + nlf - 1 + ( nlf - 1 ) * ldu,
      ioffu + nlf - 1 + ( nlf - 1 ) * ldu, ioffwork );
    if ( info.getValue() != 0 ) return;
    var itemp = idxq + nlf - 2;
    for ( var j = 1; j <= nl; j ++ ) {
      iwork[ ioffiwork + itemp + j - 1 ] = j;
    }
    sqrei = ( i == nd.getValue() ? sqre : 1 );
    nrp1 = nr + sqrei;
    LaPack6.dlasdq( 'U', sqrei, nr, nrp1, nr, ncc, d, e, Vt, ldvt,
      U, ldu, U, ldu, work, info, ioffd + nrf - 1, ioffe + nrf - 1,
      ioffvt + nrf - 1 + ( nrf - 1 ) * ldvt,
      ioffu + nrf - 1 + ( nrf - 1 ) * ldu,
      ioffu + nrf - 1 + ( nrf - 1 ) * ldu, ioffwork );
    if ( info.getValue() != 0 ) return;
    itemp = idxq + ic;
    for ( j = 1; j <= nr; j ++ ) {
      iwork[ ioffiwork + itemp + j - 2 ] = j;
    }
  }
  for ( var lvl = nlvl.getValue(); lvl >= 1; lvl -- ) {
    if ( lvl == 1 ) {
      var lf = 1;
      var ll = 1;
    } else {
      lf = Math.pow( 2, lvl - 1 );
      ll = 2 * lf - 1;
    }
    for ( i = lf; i <= ll; i ++ ) {
      var im1 = i - 1;
      var ic = iwork[ ioffiwork + inode + im1 - 1 ];
      var nl = iwork[ ioffiwork + ndiml + im1 - 1 ];
      var nr = iwork[ ioffiwork + ndimr + im1 - 1 ];
      var nlf = ic - nl;
      sqrei = ( sqre == 0 && i == ll ? sqre : 1 );
      var idxqc = idxq + nlf - 1;
      var alpha =
        new NumberReference( d[ ioffd + ic - 1 ] );
      var beta =
        new NumberReference( e[ ioffe + ic - 1 ] );
      LaPack4.dlasd1( nl, nr, sqrei, d, alpha, beta, U, ldu, Vt, ldvt,
        iwork, iwork, work, info, ioffd + nlf - 1,
        ioffu + nlf - 1 + ( nlf - 1 ) * ldu,
        ioffvt + nlf - 1 + ( nlf - 1 ) * ldvt,
        ioffiwork + idxqc - 1, ioffiwork + iwk - 1, ioffwork );
      if ( info.getValue() != 0 ) return;
    } // 40
  } // 50
}
//**************************************************************************
LaPack7.dlasda = function( icompq, smlsiz, n, sqre, d, e, U, ldu,
Vt, k, Difl, Difr, Z, Poles, givptr, Givcol, ldgcol, perm, Givnum, c, s,
work, iwork, info, ioffd, ioffe, ioffu, ioffvt, ioffk, ioffdifl, ioffdifr,
ioffz, ioffpoles, ioffgivptr, ioffgivcol, ioffperm, ioffgivnum, ioffc,
ioffs, ioffwork, ioffiwork) {
  throw new Error("not tested");
  info.setValue( 0 );
  if ( icompq < 0 || icompq > 1 ) info.setValue( -1 );
  else if ( smlsiz < 3 ) info.setValue( -2 );
  else if ( n < 0 ) info.setValue( -3 );
  else if ( sqre < 0 || sqre > 1 ) info.setValue( -4 );
  else if ( ldu < n + sqre ) info.setValue( -8 );
  else if ( ldgcol < n ) info.setValue( -17 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlasda', - info.getValue() );
    return;
  }
  var m = n + sqre;
  if ( n <= smlsiz ) {
    if ( icompq == 0 ) {
      LaPack6.dlasdq( 'U', sqre, n, 0, 0, 0, d, e, Vt, ldu, U, ldu,
        U, ldu, work, info, ioffd, ioffe, ioffvt, ioffu, ioffu,
        ioffwork );
    } else {
      LaPack6.dlasdq( 'U', sqre, n, m, n, 0, d, e, Vt, ldu, U, ldu,
        U, ldu, work, info, ioffd, ioffe, ioffvt, ioffu, ioffu,
        ioffwork );
    }
    return;
  }
  var inode = 1;
  var ndiml = inode + n;
  var ndimr = idiml + n;
  var idxq = ndimr + n;
  var iwk = idxq + n;
  var ncc = 0;
  var nru = 0;
  var smlszp = smlsiz + 1;
  var vf = 1;
  var vl = vf + m;
  var nwork1 = vl + m;
  var nwork2 = nwork1 + smlszp * smlszp;
  var nlvl = new IntReference();
  var nd = new IntReference();
  LaPack0.dlasdt( n, nlvl, nd, iwork, iwork, iwork, smlsiz,
    ioffiwork + inode - 1, ioffiwork + ndiml - 1,
    ioffiwork + ndimr - 1 );
  var ndb1 = ( nd.getValue() + 1 ) / 2;
  for ( var i = ndb1; i <= nd.getValue(); i ++ ) {
    var i1 = i - 1;
    var ic = iwork[ ioffiwork + inode + i1 - 1 ];
    var nl = iwork[ ioffiwork + ndiml + i1 - 1 ];
    var nlp1 = nl + 1;
    var nr = iwork[ ioffiwork + ndimr + i1 - 1 ];
    var nlf = ic - nl;
    var nrf = ic + 1;
    var idxqi = idxq + nlf - 2;
    var vfi = vf + nlf - 1;
    var vli = vl + nlf - 1;
    var sqrei = 1;
    if ( icompq == 0 ) {
      LaPack0.dlaset( 'A', nlp1, nlp1, 0., 1., work, smlszp,
        ioffwork + nwork1 - 1 );
      LaPack6.dlasdq( 'U', sqrei, nl, nlp1, nru, ncc, d, e, work,
        smlszp, work, nl, work, nl, work, info, ioffd + nlf - 1,
        ioffe + nlf - 1, ioffwork + nwork1 - 1, ioffwork + nwork2 - 1,
        ioffwork + nwork2 - 1, ioffwork + nwork2 - 1 );
      var itemp = nwork1 + nl * smlszp;
      Blas1.dcopy( nlp1, work, 1, work, 1, ioffwork + nwork1 - 1,
        ioffwork + vfi - 1 );
      Blas1.dcopy( nlp1, work, 1, work, 1, ioffwork + itemp - 1,
        ioffwork + vli - 1 );
    } else {
      LaPack0.dlaset( 'A', nl, nl, 0., 1., U, ldu, ioffu + nlf - 1 );
      LaPack0.dlaset( 'A', nlp1, nlp1, 0., 1., Vt, ldu,
        ioffvt + nlf - 1 );
      LaPack6.dlasdq( 'U', sqrei, nl, nlp1, nl, ncc, d, e, Vt, ldu,
        U, ldu, U, ldu, work, info, ioffd + nlf - 1, ioffe + nlf - 1,
        ioffvt + nlf - 1, ioffu + nlf - 1, ioffu + nlf - 1,
        ioffwork + nwork1 - 1 );
      Blas1.dcopy( nlp1, Vt, 1, work, 1, ioffvt + nlf - 1,
        ioffwork + vfi - 1 );
      Blas1.dcopy( nlp1, Vt, 1, work, 1,
        ioffvt + nlf - 1 + ( nlp1 - 1 ) * ldvt, ioffwork + vli - 1 );
    }
    if ( info.getValue() != 0 ) return;
    for ( var j = 1; j <= nl; j ++ ) {
      iwork[ ioffiwork + idxqi + j - 1 ] = j;
    } // 10
    var sqrei = ( i == nd.getValue() && sqre == 0 ? 0 : 1 );
    idxqi += nlp1;
    vfi += nlp1;
    vli += nlp1;
    var nrp1 = nr + sqrei;
    if ( icompq == 0 ) {
      LaPack0.dlaset( 'A', nrp1, nrp1, 0., 1., work, smlszp,
        ioffwork + nwork1 - 1 );
      LaPack6.dlasdq( 'U', sqrei, nr, nrp1, nru, ncc, d, e, work,
        smlszp, work, nr, work, nr, work, info, ioffd + nrf - 1,
        ioffe + nrf - 1, ioffwork + nwork1 - 1, ioffwork + nwork2 - 1,
        ioffwork + nwork2 - 1, ioffwork + nwork2 - 1 );
      itemp = nwork1 + ( nrp1 - 1 ) * smlszp;
      Blas1.dcopy( nrp1, work, 1, work, 1, ioffwork + nwork1 - 1,
        ioffwork + vfi - 1 );
      Blas1.dcopy( nrp1, work, 1, work, 1, ioffwork + itemp - 1,
        ioffwork + vli - 1 );
    } else {
      LaPack0.dlaset( 'A', nr, nr, 0., 1., U, ldu, ioffu + nrf - 1 );
      LaPack0.dlaset( 'A', nrp1, nrp1, 0., 1., Vt, ldu,
        ioffvt + nrf - 1 );
      LaPack6.dlasdq( 'U', sqrei, nr, nrp1, nr, ncc, d, e, Vt, ldu,
        U, ldu, U, ldu, work, info, ioffd + nrf - 1, ioffe + nrf - 1,
        ioffvt + nrf - 1, ioffu + nrf - 1, ioffu + nrf - 1,
        ioffwork + nwork1 - 1 );
      Blas1.dcopy( nrp1, Vt, 1, work, 1, ioffvt + nrf - 1,
        ioffwork + vfi - 1 );
      Blas1.dcopy( nrp1, Vt, 1, work, 1,
        ioffvt + nrf - 1 + ( nrp1 - 1 ) * ldvt, ioffwork + vli - 1 );
    }
    if ( info.getValue() != 0 ) return;
    for ( j = 1; j <= nr; j ++ ) iwork[ ioffiwork + idxqi + j - 1 ] = j;
  } // 30
  j = Math.pow( 2, nlvl.getValue() );
  for ( var lvl = nlvl.getValue(); lvl >= 1; lvl -- ) {
    var lvl2 = lvl * 2 - 1;
    if ( lvl == 1 ) {
      var lf = 1;
      var ll = 1;
    } else {
      lf = Math.pow( 2, lvl - 1 );
      ll = 2 * lf - 1;
    }
    for ( i = lf; i <= ll; i ++ ) {
      var im1 = i - 1;
      var ic = iwork[ ioffiwork + inode + im1 - 1 ];
      var nl = iwork[ ioffiwork + ndiml + im1 - 1 ];
      var nr = iwork[ ioffiwork + ndimr + im1 - 1 ];
      nlf = i  - nl;
      nrf = ic + 1;
      sqrei = ( i == ll ? sqre : 1 );
      vfi = vf + nlf - 1;
      vli = vl + nlf - 1;
      idxqi = idxq + nlf - 1;
      var alpha =
        new NumberReference( d[ ioffd + ic - 1 ] );
      var beta =
        new NumberReference( e[ ioffe + ic - 1 ] );
      var kref = new IntReference();
      var cref = new NumberReference();
      var sref = new NumberReference();
      if ( icompq == 0 ) {
        kref.setValue(k[ ioffk ] );
        cref.setValue( c[ ioffc ] );
        sref.setValue( s[ ioffc ] );
        var givptrReference = new IntReference( givptr[ ioffgivptr ] );
        LaPack4.dlasd6( icompq, nl, nr, sqrei,
          d, work, work, alpha,
          beta, iwork, perm, givptr[ ioffgivptr ],
          Givcol, ldgcol, Givnum, ldu, Poles,
          Difl, Difr, Z, kref, cref,
          sref, work, iwork, info,
          ioffd + nlf - 1, ioffwork + vfi - 1, ioffwork + vli - 1,
            ioffiwork + idxqi - 1, ioffperm,
          ioffgivcol, ioffgivnum, ioffpoles, ioffdifl,
          ioffdifr, ioffz, ioffwork + nwork1 - 1, ioffwork + iwk - 1 );
        k[ ioffk ] = kref.getValue();
        c[ ioffc ] = cref.getValue();
        s[ ioffs ] = sref.getValue();
        givptr[ ioffgivptr ] = givptrReference.getValue();
      } else {
        j --;
        kref.setValue( k[ ioffk + j - 1 ] );
        cref.setValue( c[ ioffc + j - 1 ] );
        sref.setValue( s[ ioffc + j - 1 ] );
        givptrReference = new IntReference( givptr[ ioffgivptr + j - 1 ] );
        LaPack4.dlasd6( icompq, nl, nr, sqrei,
          d, work, work, alpha,
          beta, iwork, perm, givptr[ ioffgivptr + j - 1 ],
          Givcol, ldgcol, Givnum, ldu, Poles,
          Difl, Difr, Z, kref, cref,
          sref, work, iwork, info,
          ioffd + nlf - 1, ioffwork + vfi - 1, ioffwork + vli - 1,
            ioffiwork + idxqi - 1,
            ioffperm + nlf - 1 + ( lvl - 1 ) * ldgcol,
          ioffgivcol + nlf - 1 + ( lvl2 - 1 ) * ldgcol,
          ioffgivnum + nlf - 1 + ( lvl2 - 1 ) * ldu,
          ioffpoles + nlf - 1 + ( lvl2 - 1 ) * ldu,
          ioffdifl + nlf - 1 + ( lvl - 1 ) * ldu,
          ioffdifr + nlf - 1 + ( lvl2 - 1 ) * ldu,
          ioffz + nlf - 1 + ( lvl - 1 ) * ldu, ioffwork + nwork1 - 1,
          ioffwork + iwk - 1 );
        k[ ioffk + j - 1 ] = kref.getValue();
        c[ ioffc + j - 1 ] = cref.getValue();
        s[ ioffs + j - 1 ] = sref.getValue();
        givptr[ ioffgivptr +j - 1 ] = givptrReference.getValue();
      }
      if ( info.getValue() == 0 ) return;
    } // 40
  } // 50
}
//**************************************************************************
LaPack7.dspevd = function( jobz, uplo, n, AP, w, Z, ldz, work,
lwork, iwork, liwork, info, ioffap, ioffw, ioffz, ioffwork, ioffiwork ) {
  throw new Error("not programmed: packed matrix");
}
//**************************************************************************
LaPack7.dsyevd = function( jobz, uplo, n, A, lda, w, work, lwork,
iwork, liwork, info, ioffa, ioffw, ioffwork, ioffiwork ) {
  throw new Error("not tested");
  var wantz = ( jobz.charAt(0).toUpperCase() == 'V' );
  var lower = ( uplo.charAt(0).toUpperCase() == 'L' );
  var lquery = ( lwork == -1 || liwork == -1 );
  info.setValue( 0 );
  if ( ! ( wantz || jobz.charAt(0).toUpperCase() == 'N' ) ) {
    info.setValue( -1 );
  } else if ( ! ( lower || uplo.charAt(0).toUpperCase() == 'U' ) ) {
    info.setValue( -2 );
  } else if ( n < 0 ) info.setValue( -3 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -5 );
  if ( info.getValue() == 0 ) {
    if ( n <= 1 ) {
      var liwmin = 1;
      var lwmin = 1;
      var lopt = lwmin;
      var liopt = liwmin;
    } else {
      if ( wantz ) {
        liwmin = 3 + 5 * n;
        lwmin = 1 + 6 * n + 2 * n * n;
      } else {
        liwmin = 1;
        lwmin = 2 * n + 1;
      }
      lopt = Math.max( lwmin, 2 * n
        + LaPack0.ilaenv( 1, 'dsytrd', uplo, n, -1, -1, -1 ) );
      liopt = liwmin;
    }
    work[ ioffwork ] = lopt;
    iwork[ ioffiwork ] = liopt;
    if ( lwork < lwmin && ! lquery ) info.setValue( -8 );
    else if ( liwork < liwmin && ! lquery ) info.setValue( -10 );
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dsyevd', - info.getValue() );
    return;
  } else if ( lquery ) return;
  if ( n == 0 ) return;
  if ( n == 1 ) {
    w[ ioffw ] = A[ ioffa ];
    if ( wantz ) A[ ioffa ] = 1.;
    return;
  }
  var safmin = LaPack0.dlamch( 'Safe minimum' );
  var eps = LaPack0.dlamch( 'Precision' );
  var smlnum = safmin / eps;
  var bignum = 1. / smlnum;
  var rmin = Math.sqrt( smlnum );
  var rmax = Math.sqrt( bignum );
  var anrm =
    LaPack1.dlansy( 'M', uplo, n, A, lda, work, ioffa, ioffwork );
  var iscale = 0;
  if ( anrm > 0. && anrm < rmin ) {
    var iscale = 1;
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
  var indwk2 = indwrk + n * n;
  var llwrk2 = lwork - indwk2 + 1;
  var iinfo = new IntReference();
  LaPack3.dsytrd( uplo, n, A, lda, w, work, work, work, llwork, iinfo,
    ioffa, ioffw, ioffwork + inde - 1, ioffwork + indtau - 1,
    ioffwork + indwrk - 1 );
  var lopt = 2 * n + work[ ioffwork + indwrk - 1 ];
  if ( ! wantz ) {
    LaPack2.dsterf( n, w, work, info, ioffw, ioffwork + inde - 1 );
  } else {
    LaPack6.dstedc( 'I', n, w, work, work, n, work, llwrk2, iwork,
      liwork, info, ioffw, ioffwork + inde - 1, ioffwork + indwrk - 1,
      ioffwork + indwk2 - 1, ioffiwork );
    LaPack4.dormtr( 'L', uplo, 'N', n, n, A, lda, work, work, n, work,
      llwrk2, iinfo, ioffa, ioffwork + indtau - 1,
      ioffwork + indwrk - 1, ioffwork + indwk2 - 1 );
    LaPack0.dlacpy( 'A', n, n, work, n, A, lda, ioffwork + indwrk - 1,
      ioffa );
    lopt = Math.max( lopt, 1 + 6 * n + 2 * n * n );
  }
  if ( iscale == 1 ) {
    Blas1.dscal( n, 1. / sigma, w, 1, ioffw );
  }
  work[ ioffwork ] = lopt;
  iwork[ ioffiwork ] = liopt;
}
//*************************************************************************
LaPack7.dtgsen = function( ) {
  throw new Error("not programmed: generalized matrix");
}
//*************************************************************************
LaPack7.dtgsna = function( ) {
  throw new Error("not programmed: generalized matrix");
}
