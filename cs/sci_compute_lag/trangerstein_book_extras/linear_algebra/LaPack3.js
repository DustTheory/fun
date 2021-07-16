function LaPack3() {
}
//*************************************************************************
LaPack3.dbbcsd = function( jobu1, jobu2, jobv1t, jobv2t, trans,
m, p, q, theta, phi, U1, ldu1, U2, ldu2, V1t, ldv1t, V2t, ldv2t, b11d,
b11e, b12d, b12e, b21d, b21e, b22d, b22e, work, lwork, info, iofftheta,
ioffphi, ioffu1, ioffu2, ioffv1t, ioffv2t, ioffb11d, ioffb11e, ioffb12d,
ioffb12e, ioffb21d, ioffb21e, ioffb22d, ioffb22e, ioffwork ) {
  throw new Error("not tested");
  var maxitr = 6;
  var piover2 = 1.57079632679489662;
  info.setValue( 0 );
  var lquery = ( lwork == -1 );
  var wantu1 = ( jobu1.charAt(0).toUpperCase() == 'Y' );
  var wantu2 = ( jobu2.charAt(0).toUpperCase() == 'Y' );
  var wantv1t = ( jobv1t.charAt(0).toUpperCase() == 'Y' );
  var wantv2t = ( jobv2t.charAt(0).toUpperCase() == 'Y' );
  var colmajor = ( trans.charAt(0).toUpperCase() != 'T' );
  if ( m < 0 ) info.setValue( -6 );
  else if ( p < 0 || p > m ) info.setValue( -7 );
  else if ( q < 0 || q > m ) info.setValue( -8 );
  else if ( q > p || q > m - p || q > m - q ) info.setValue( -8 );
  else if ( wantu1 && ldu1 < p ) info.setValue( -12 );
  else if ( wantu2 && ldu2 < m - p ) info.setValue( -14 );
  else if ( wantv1t && ldv1t < q ) info.setValue( -16 );
  else if ( wantv2t && ldv2t < m - q ) info.setValue( -18 );
  if ( info.getValue() == 0 && q == 0 ) {
    var lworkmin = 1;
    work[ ioffwork ] = lworkmin;
    return;
  }
  if ( info.getValue() == 0 ) {
    var iu1cs = 1;
    var iu1sn = iu1cs + q;
    var iu2cs = iu1sn + q;
    var iu2sn = iu2cs + q;
    var iv1tcs = iu2sn + q;
    var iv1tsn = iv1tcs + q;
    var iv2tcs = iv1tsn + q;
    var iv2tsn = iv2tcs + q;
    var lworkopt = iv2tsn + q - 1;
    lworkmin = lworkopt;
    work[ ioffwork ] = lworkopt;
    if ( lwork < lworkmin && ! lquery ) info.setValue( -28 );
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dbbcsd', - info.getValue() );
    return;
  } else if ( lquery ) return;
  var eps = LaPack0.dlamch( 'Epsilon' );
  var unfl = LaPack0.dlamch( 'Safe minimum' );
  var tolmul =
    Math.max( 10., Math.min( 100., Math.pow( eps, -0.125 ) ) );
  var tol = tolmul * eps;
  var thresh = Math.max( tol, maxitr * q * q * unfl );
  for ( var i = 1; i <= q; i ++ ) {
    if ( theta[ iofftheta + i - 1 ] < thresh ) {
      theta[ iofftheta + i - 1 ] = 0.;
    } else if ( theta[ iofftheta + i - 1 ] > piover2 - thresh ) {
      theta[ iofftheta + i - 1 ] = piover2;
    }
  }
  for ( i = 1; i <= q - 1; i ++ ) {
    if ( phi[ ioffphi + i - 1 ] < thresh ) {
      phi[ ioffphi + i - 1 ] = 0.;
    } else if ( phi[ ioffphi + i - 1 ] > piover2 - thresh ) {
      phi[ ioffphi + i - 1 ] = piover2;
    }
  }
  var imax = q;
  while ( imax > 1 && phi[ ioffphi + imax - 2 ] == 0. ) imax --;
  var imin = imax - 1;
  if ( imin > 1 ) {
    while ( phi[ ioffphi + imin - 2 ] != 0. ) {
      imin --;
      if ( imin <= 1 ) break;
    }
  }
  var maxit = maxitr* q * q;
  var iter = 0;
  while ( imax > 1 ) {
    b11d[ ioffb11d + imin - 1 ] =
      Math.cos( theta[ iofftheta + imin - 1 ] );
    b21d[ ioffb21d + imin - 1 ] =
      - Math.sin( theta[ iofftheta + imin - 1 ] );
    for ( i = imin; i <= imax - 1; i ++ ) {
      b11e[ ioffb11e + i - 1 ] =
        - Math.sin( theta[ iofftheta + i - 1 ] )
        * Math.sin( phi[ ioffphi + i - 1 ] );
      b11d[ ioffb11d + i ] =
        Math.cos( theta[ iofftheta + i ] )
        * Math.cos( phi[ ioffphi + i - 1 ] );
      b12d[ ioffb12d + i - 1 ] =
        Math.sin( theta[ iofftheta + i - 1 ] )
        * Math.cos( phi[ ioffphi + i - 1 ] );
      b12e[ ioffb12e + i - 1 ] =
        Math.cos( theta[ iofftheta + i ] )
        * Math.sin( phi[ ioffphi + i - 1 ] );
      b21e[ ioffb21e + i - 1 ] =
        - Math.cos( theta[ iofftheta + i - 1 ] )
        * Math.sin( phi[ ioffphi + i - 1 ] );
      b21d[ ioffb21d + i ] =
        - Math.sin( theta[ iofftheta + i ] )
        * Math.cos( phi[ ioffphi + i - 1 ] );
      b22d[ ioffb22d + i - 1 ] =
        Math.cos( theta[ iofftheta + i - 1 ] )
        * Math.cos( phi[ ioffphi + i - 1 ] );
      b22e[ ioffb22e + i - 1 ] =
        - Math.sin( theta[ iofftheta + i ] )
        * Math.sin( phi[ ioffphi + i - 1 ] );
    }
    b12d[ ioffb12d + imax - 1 ] =
      Math.sin( theta[ iofftheta + imax - 1 ] );
    b22d[ ioffb22d + imax - 1 ] =
      Math.cos( theta[ iofftheta + imax - 1 ] );
    if ( iter > maxit ) {
      info.setValue( 0 );
      for ( i = 1; i <= q; i ++ ) {
        if ( phi[ ioffphi + i - 1 ] != 0. ) {
          info.setValue( info.getValue() + 1 );
        }
      }
      return;
    }
    iter += imax - imin;
    var thetamax = theta[ iofftheta + imin ];
    var thetamin = theta[ iofftheta + imin ];
    for ( i = imin + 1; i <= imax; i ++ ) {
      if ( theta[ iofftheta + i - 1 ] > thetamax ) {
        thetamax = theta[ iofftheta + i - 1 ];
      }
      if ( theta[ iofftheta + i - 1 ] < thetamin ) {
        thetamin = theta[ iofftheta + i - 1 ];
      }
    }
    if ( thetamax > piover2 - thresh ) {
      var mu = 0.;
      var nu = 1.;
    } else if ( thetamin < thresh ) {
      mu = 1.;
      nu = 0.;
    } else {
      var sigma11Reference = new NumberReference();
      var dummyReference = new NumberReference();
      LaPack0.dlas2( b11d[ ioffb11d + imax - 2 ],
        b11e[ ioffb11e + imax - 2 ], b11d[ ioffb11d + imax - 1 ], 
        sigma11, dummy );
      var sigma21Reference = new NumberReference();
      LaPack0.dlas2( b21d[ ioffb21d + imax - 2 ],
        b21e[ ioffb21e + imax - 2 ], b21d[ ioffb21d + imax - 1 ], 
        sigma21, dummy );
      if ( sigma11.getValue() <= sigma21.getValue() ) {
        mu = sigma11.getValue();
        nu = Math.sqrt( 1. - mu * mu );
        if ( mu < thresh ) {
          mu = 0.;
          nu = 1.;
        }
      } else {
        nu = sigma21.getValue();
        mu = Math.sqrt( 1. - nu * nu );
        if ( nu < thresh ) {
          mu = 1.;
          nu = 0.;
        }
      }
    }
    var csReference =
      new NumberReference( work[ ioffwork + iv1tcs + imin - 2 ] );
    var snReference =
      new NumberReference( work[ ioffwork + iv1tsn + imin - 2 ] );
    if ( mu <= nu ) {
      LaPack2.dlartgs( b11d[ ioffb11d + imin - 1 ],
        b11e[ ioffb11e + imin - 1 ], mu, cs, sn );
    } else {
      LaPack2.dlartgs( b21d[ ioffb21d + imin - 1 ],
        b21e[ ioffb21e + imin - 1 ], nu, cs, sn );
    }
    work[ ioffwork + iv1tcs + imin - 2 ] = cs.getValue();
    work[ ioffwork + iv1tsn + imin - 2 ] = sn.getValue();
    var temp = work[ ioffwork + iv1tcs + imin - 2 ]
      * b11d[ ioffb11d + imin - 1 ]
      + work[ ioffwork + iv1tsn + imin - 2 ]
      * b11e[ ioffb11e + imin - 1 ];
    b11e[ ioffb11e + imin - 1 ] = work[ ioffwork + iv1tcs + imin - 2 ]
      * b11e[ ioffb11e + imin - 1 ]
      - work[ ioffwork + iv1tsn + imin - 2 ]
      * b11d[ ioffb11d + imin - 1 ];
    b11d[ ioffb11d + imin - 1 ] = temp;
    var b11bulge = work[ ioffwork + iv1tsn + imin - 2 ]
      * b11d[ ioffb11d + imin ];
    b11d[ ioffb11d + imin ] *= work[ ioffwork + iv1tcs + imin - 2 ];
    temp = work[ ioffwork + iv1tcs + imin - 2 ]
      * b21d[ ioffb21d + imin - 1 ]
      + work[ ioffwork + iv1tsn + imin - 2 ]
      * b21e[ ioffb21e + imin - 1 ];
    b21e[ ioffb21e + imin - 1 ] = work[ ioffwork + iv1tcs + imin - 2 ]
      * b21e[ ioffb21e + imin - 1 ]
      - work[ ioffwork + iv1tsn + imin - 2 ]
      * b21d[ ioffb21d + imin - 1 ];
    b21d[ ioffb21d + imin - 1 ] = temp;
    var b21bulge = work[ ioffwork + iv1tsn + imin - 2 ]
      * b21d[ ioffb21d + imin ];
    b21d[ ioffb21d + imin ] = work[ ioffwork + iv1tcs + imin - 2 ]
      * b21d[ ioffb21d + imin ];
    theta[ iofftheta + imin - 1 ] =
      Math.atan2( Math.sqrt(
        Math.pow(  b21d[ ioffb21d + imin - 1 ], 2 )
        + b21bulge * b21bulge ), Math.sqrt(
        Math.pow(  b11d[ ioffb11d + imin - 1 ], 2 )
        + b11bulge * b11bulge ) );
    if ( Math.pow( b11d[ ioffb11d + imin - 1 ], 2 )
    + b11bulge * b11bulge > thresh * thresh ) {
      cs.setValue( work[ ioffwork + iu1sn + imin - 2 ] );
      sn.setValue( work[ ioffwork + iu1cs + imin - 2 ] );
      var rReference = new NumberReference();
      LaPack1.dlartgp( b11bulge, b11d[ ioffb11d + imin - 1 ],
        cs, sn, r );
      work[ ioffwork + iu1sn + imin - 2 ] = cs.getValue();
      work[ ioffwork + iu1cs + imin - 2 ] = sn.getValue();
    } else if ( mu <= nu ) {
      cs.setValue( work[ ioffwork + iu1cs + imin - 2 ] );
      sn.setValue( work[ ioffwork + iu1sn + imin - 2 ] );
      LaPack2.dlartgs( b11e[ ioffb11e + imin - 1 ],
        b11d[ ioffb11d + imin ], mu, cs, sn );
      work[ ioffwork + iu1cs + imin - 2 ] = cs.getValue();
      work[ ioffwork + iu1sn + imin - 2 ] = sn.getValue();
    } else {
      cs.setValue( work[ ioffwork + iu1cs + imin - 2 ] );
      sn.setValue( work[ ioffwork + iu1sn + imin - 2 ] );
      LaPack2.dlartgs( b12d[ ioffb12d + imin - 1 ],
        b12e[ ioffb12e + imin - 1 ], nu, cs, sn );
      work[ ioffwork + iu1cs + imin - 2 ] = cs.getValue();
      work[ ioffwork + iu1sn + imin - 2 ] = sn.getValue();
    }
    if ( Math.pow( b21d[ ioffb21d + imin - 1 ], 2 )
    + b21bulge * b21bulge > thresh * thresh ) {
      cs.setValue( work[ ioffwork + iu2sn + imin - 2 ] );
      sn.setValue( work[ ioffwork + iu2cs + imin - 2 ] );
      LaPack1.dlartgp( b21bulge, b21d[ ioffb21d + imin - 1 ],
        cs, sn, r );
      work[ ioffwork + iu2sn + imin - 2 ] = cs.getValue();
      work[ ioffwork + iu2cs + imin - 2 ] = sn.getValue();
    } else if ( nu < mu ) {
      cs.setValue( work[ ioffwork + iu2cs + imin - 2 ] );
      sn.setValue( work[ ioffwork + iu2sn + imin - 2 ] );
      LaPack2.dlartgs( b21e[ ioffb21e + imin - 1 ],
        b21d[ ioffb21d + imin ], nu, cs, sn );
      work[ ioffwork + iu2cs + imin - 2 ] = cs.getValue();
      work[ ioffwork + iu2sn + imin - 2 ] = sn.getValue();
    } else {
      cs.setValue( work[ ioffwork + iu2cs + imin - 2 ] );
      sn.setValue( work[ ioffwork + iu2sn + imin - 2 ] );
      LaPack2.dlartgs( b22d[ ioffb22d + imin - 1 ],
        b22e[ ioffb22e + imin - 1 ], mu, cs, sn );
      work[ ioffwork + iu2cs + imin - 2 ] = cs.getValue();
      work[ ioffwork + iu2sn + imin - 2 ] = sn.getValue();
    }
    work[ ioffwork + iu2cs + imin - 2 ] =
      - work[ ioffwork + iu2cs + imin - 2 ];
    work[ ioffwork + iu2sn + imin - 2 ] =
      - work[ ioffwork + iu2sn + imin - 2 ];
    temp = work[ ioffwork + iu1cs + imin - 2 ]
      * b11e[ ioffb11e + imin - 1 ]
      + work[ ioffwork + iu1sn + imin - 2 ]
      * b11d[ ioffb11d + imin ];
    b11d[ ioffb11d + imin ] = work[ ioffwork + iu1cs + imin - 2 ]
      * b11d[ ioffb11d + imin ]
      - work[ ioffwork + iu1sn + imin - 2 ]
      * b11e[ ioffb11e + imin - 1 ];
    b11e[ ioffb11e + imin - 1 ] = temp;
    if ( imax > imin + 1 ) {
      b11bulge = work[ ioffwork + iu1sn + imin - 2 ]
        * b11e[ ioffb11e + imin ];
      b11e[ ioffb11e + imin ] = work[ ioffwork + iu1cs + imin - 2 ]
        * b11e[ ioffb11e + imin ];
    }
    temp = work[ ioffwork + iu1cs + imin - 2 ]
      * b12d[ ioffb12d + imin - 1 ]
      + work[ ioffwork + iu1sn + imin - 2 ]
      * b12e[ ioffb12e + imin - 1 ];
    b12e[ ioffb12e + imin - 1 ] = work[ ioffwork + iu1cs + imin - 2 ]
      * b12e[ ioffb12e + imin - 1 ]
      - work[ ioffwork + iu1sn + imin - 2 ]
      * b12d[ ioffb12d + imin - 1 ];
    b12d[ ioffb12d + imin - 1 ] = temp;
    var b12bulge = work[ ioffwork + iu1sn + imin - 2 ]
      * b12d[ ioffb12d + imin ];
    b12d[ ioffb12d + imin ] = work[ ioffwork + iu1cs + imin - 2 ]
      * b12d[ ioffb12d + imin ];
    temp = work[ ioffwork + iu2cs + imin - 2 ]
      * b21e[ ioffb21e + imin - 1 ]
      + work[ ioffwork + iu2sn + imin - 2 ]
      * b21d[ ioffb21d + imin ];
    b21d[ ioffb21d + imin ] = work[ ioffwork + iu2cs + imin - 2 ]
      * b21d[ ioffb21d + imin ]
      - work[ ioffwork + iu2sn + imin - 2 ]
      * b21e[ ioffb21e + imin - 1 ];
    b21e[ ioffb21e + imin - 1 ] = temp;
    if ( imax > imin + 1 ) {
      b21bulge = work[ ioffwork + iu2sn + imin - 2 ]
        * b21e[ ioffb21e + imin ];
      b21e[ ioffb21e + imin ] = work[ ioffwork + iu2cs + imin - 2 ]
        * b21e[ ioffb21e + imin ];
    }
    temp = work[ ioffwork + iu2cs + imin - 2 ]
      * b22d[ ioffb22d + imin - 1 ]
      + work[ ioffwork + iu2sn + imin - 2 ]
      * b22e[ ioffb22e + imin - 1 ];
    b22e[ ioffb22e + imin - 1 ] = work[ ioffwork + iu2cs + imin - 2 ]
      * b22e[ ioffb22e + imin - 1 ]
      - work[ ioffwork + iu2sn + imin - 2 ]
      * b22d[ ioffb22d + imin - 1 ];
    b22d[ ioffb22d + imin - 1 ] = temp;
    var b22bulge = work[ ioffwork + iu2sn + imin - 2 ]
      * b22d[ ioffb22d + imin ];
    b22d[ ioffb22d + imin ] *= work[ ioffwork + iu2cs + imin - 2 ];
    for ( i = imin + 1; i <= imax - 1; i ++ ) {
      var x1 = Math.sin( theta[ iofftheta + i - 2 ] )
        * b11e[ ioffb11e + i - 2 ]
        + Math.cos( theta[ iofftheta + i - 2 ] )
        * b21e[ ioffb21e + i - 2 ];
      var x2 = Math.sin( theta[ iofftheta + i - 2 ] ) * b11bulge
        + Math.cos( theta[ iofftheta + i - 2 ] ) * b21bulge;
      var y1 = Math.sin( theta[ iofftheta + i - 2 ] )
        * b12d[ ioffb12d + i - 2 ]
        + Math.cos( theta[ iofftheta + i - 2 ] )
        * b22d[ ioffb22d + i - 2 ];
      var y2 = Math.sin( theta[ iofftheta + i - 2 ] ) * b12bulge
        + Math.cos( theta[ iofftheta + i - 2 ] ) * b22bulge;
      phi[ ioffphi + i - 2 ] = Math.atan2(
        Math.sqrt( x1 * x1 + x2 * x2 ),
        Math.sqrt( y1 * y1 + y2 * y2 ) );
      var restart11 = Math.pow( b11e[ ioffb11e + i - 2 ], 2 )
        + b11bulge * b11bulge <= thresh * thresh;
      var restart21 = Math.pow( b21e[ ioffb21e + i - 2 ], 2 )
        + b21bulge * b21bulge <= thresh * thresh;
      var restart12 = Math.pow( b12d[ ioffb12d + i - 2 ], 2 )
        + b12bulge * b12bulge <= thresh * thresh;
      var restart22 = Math.pow( b22d[ ioffb22d + i - 2 ], 2 )
        + b22bulge * b22bulge <= thresh * thresh;
      if ( ! restart11 && ! restart21 ) {
        cs.setValue( work[ ioffwork + iv1tsn + i - 2 ] );
        sn.setValue( work[ ioffwork + iv1tcs + i - 2 ] );
        LaPack1.dlartgp( x2, x1, cs, sn, r );
        work[ ioffwork + iv1tsn + i - 2 ] = cs.getValue();
        work[ ioffwork + iv1tcs + i - 2 ] = sn.getValue();
      } else if ( ! restart11 && restart21 ) {
        cs.setValue( work[ ioffwork + iv1tsn + i - 2 ] );
        sn.setValue( work[ ioffwork + iv1tcs + i - 2 ] );
        LaPack1.dlartgp( b11bulge, b11e[ ioffb11e + i - 2 ],
          cs, sn, r );
        work[ ioffwork + iv1tsn + i - 2 ] = cs.getValue();
        work[ ioffwork + iv1tcs + i - 2 ] = sn.getValue();
      } else if ( restart11 && ! restart21 ) {
        cs.setValue( work[ ioffwork + iv1tsn + i - 2 ] );
        sn.setValue( work[ ioffwork + iv1tcs + i - 2 ] );
        LaPack1.dlartgp( b21bulge, b21e[ ioffb21e + i - 2 ],
          cs, sn, r );
        work[ ioffwork + iv1tsn + i - 2 ] = cs.getValue();
        work[ ioffwork + iv1tcs + i - 2 ] = sn.getValue();
      } else if ( mu <= nu ) {
        cs.setValue( work[ ioffwork + iv1tcs + i - 2 ] );
        sn.setValue( work[ ioffwork + iv1tsn + i - 2 ] );
        LaPack2.dlartgs( b11d[ ioffb11d + i - 1 ],
          b11e[ ioffb11e + i - 1 ], mu, cs, sn );
        work[ ioffwork + iv1tcs + i - 2 ] = cs.getValue();
        work[ ioffwork + iv1tsn + i - 2 ] = sn.getValue();
      } else {
        cs.setValue( work[ ioffwork + iv1tcs + i - 2 ] );
        sn.setValue( work[ ioffwork + iv1tsn + i - 2 ] );
        LaPack2.dlartgs( b21d[ ioffb21d + i - 1 ],
          b21e[ ioffb21e + i - 1 ], nu, cs, sn );
        work[ ioffwork + iv1tcs + i - 2 ] = cs.getValue();
        work[ ioffwork + iv1tsn + i - 2 ] = sn.getValue();
      }
      work[ ioffwork + iv1tcs + i - 2 ] =
        - work[ ioffwork + iv1tcs + i - 2 ];
      work[ ioffwork + iv1tsn + i - 2 ] =
        - work[ ioffwork + iv1tsn + i - 2 ];
      if ( ! restart12 && ! restart22 ) {
        cs.setValue( work[ ioffwork + iv2tsn + i - 3 ] );
        sn.setValue( work[ ioffwork + iv2tcs + i - 3 ] );
        LaPack1.dlartgp( y2, y1, cs, sn, r );
        work[ ioffwork + iv2tsn + i - 3 ] = cs.getValue();
        work[ ioffwork + iv2tcs + i - 3 ] = sn.getValue();
      } else if ( ! restart12 && restart22 ) {
        cs.setValue( work[ ioffwork + iv2tsn + i - 3 ] );
        sn.setValue( work[ ioffwork + iv2tcs + i - 3 ] );
        LaPack1.dlartgp( b12bulge, b12d[ ioffb12d + i - 2 ],
          cs, sn, r );
        work[ ioffwork + iv2tsn + i - 3 ] = cs.getValue();
        work[ ioffwork + iv2tcs + i - 3 ] = sn.getValue();
      } else if ( restart12 && ! restart22 ) {
        cs.setValue( work[ ioffwork + iv2tsn + i - 3 ] );
        sn.setValue( work[ ioffwork + iv2tcs + i - 3 ] );
        LaPack1.dlartgp( b22bulge, b22d[ ioffb22d + i - 2 ],
          cs, sn, r );
        work[ ioffwork + iv2tsn + i - 3 ] = cs.getValue();
        work[ ioffwork + iv2tcs + i - 3 ] = sn.getValue();
      } else if ( nu < mu ) {
        cs.setValue( work[ ioffwork + iv2tcs + i - 3 ] );
        sn.setValue( work[ ioffwork + iv2tsn + i - 3 ] );
        LaPack2.dlartgs( b12e[ ioffb12e + i - 2 ],
          b12d[ ioffb12d + i - 1 ], nu, cs, sn );
        work[ ioffwork + iv2tcs + i - 3 ] = cs.getValue();
        work[ ioffwork + iv2tsn + i - 3 ] = sn.getValue();
      } else {
        cs.setValue( work[ ioffwork + iv2tcs + i - 3 ] );
        sn.setValue( work[ ioffwork + iv2tsn + i - 3 ] );
        LaPack2.dlartgs( b22e[ ioffb22e + i - 2 ],
          b22d[ ioffb22d + i - 1 ], mu, cs, sn );
        work[ ioffwork + iv2tcs + i - 3 ] = cs.getValue();
        work[ ioffwork + iv2tsn + i - 3 ] = sn.getValue();
      }
      temp = work[ ioffwork + iv1tcs + i - 2 ]
        * b11d[ ioffb11d + i - 1 ]
        + work[ ioffwork + iv1tsn + i - 2 ]
        * b11e[ ioffb11e + i - 1 ];
      b11e[ ioffb11e + i - 1 ] = work[ ioffwork + iv1tcs + i - 2 ]
        * b11e[ ioffb11e + i - 1 ]
        - work[ ioffwork + iv1tsn + i - 2 ]
        * b11d[ ioffb11d + i - 1 ];
      b11d[ ioffb11d + i - 1 ] = temp;
      b11bulge = work[ ioffwork + iv1tsn + i - 2 ]
        * b11d[ ioffb11d + i ];
      b11bulge = work[ ioffwork + iv1tsn + i - 2 ]
        * b11d[ ioffb11d + i ];
      b11d[ ioffb11d + i ] = work[ ioffwork + iv1tcs + i - 2 ]
        * b11d[ ioffb11d + i ];
      temp = work[ ioffwork + iv1tcs + i - 2 ]
        * b21d[ ioffb21d + i - 1 ]
        + work[ ioffwork + iv1tsn + i - 2 ]
        * b21e[ ioffb21e + i - 1 ];
      b21e[ ioffb21e + i - 1 ] = work[ ioffwork + iv1tcs + i - 2 ]
        * b21e[ ioffb21e + i - 1 ]
        - work[ ioffwork + iv1tsn + i - 2 ]
        * b21d[ ioffb21d + i - 1 ];
      b21d[ ioffb21d + i - 1 ] = temp;
      b21bulge = work[ ioffwork + iv1tsn + i - 2 ]
        * b21d[ ioffb21d + i ];
      b21d[ ioffb21d + i ] = work[ ioffwork + iv1tcs + i - 2 ]
        * b21d[ ioffb21d + i ];
      temp = work[ ioffwork + iv2tcs + i - 3 ]
        * b12e[ ioffb12e + i - 2 ]
        + work[ ioffwork + iv2tsn + i - 3 ]
        * b12d[ ioffb12d + i - 1 ];
      b12d[ ioffb12d + i - 1 ] = work[ ioffwork + iv2tcs + i - 3 ]
        * b12d[ ioffb12d + i - 1 ]
        - work[ ioffwork + iv2tsn + i - 3 ]
        * b12e[ ioffb12e + i - 2 ];
      b12e[ ioffb12e + i - 2 ] = temp;
      b12bulge = work[ ioffwork + iv2tsn + i - 3 ]
        * b12e[ ioffb12e + i - 1 ];
      b12e[ ioffb12e + i - 1 ] = work[ ioffwork + iv2tcs + i - 3 ]
        * b12e[ ioffb12e + i - 1 ];
      temp = work[ ioffwork + iv2tcs + i - 3 ]
        * b22e[ ioffb21d + i - 2 ]
        + work[ ioffwork + iv2tsn + i - 3 ]
        * b22d[ ioffb22d + i - 1 ];
      b22d[ ioffb22d + i - 1 ] = work[ ioffwork + iv2tcs + i - 3 ]
        * b22d[ ioffb22d + i - 1 ]
        - work[ ioffwork + iv2tsn + i - 3 ]
        * b22e[ ioffb22e + i - 2 ];
      b22e[ ioffb22e + i - 2 ] = temp;
      b22bulge = work[ ioffwork + iv2tsn + i - 3 ]
        * b22e[ ioffb22e + i - 1 ];
      b22e[ ioffb22e + i - 1 ] = work[ ioffwork + iv2tcs + i - 3 ]
        * b22e[ ioffb22e + i - 1 ];
      x1 = Math.cos( phi[ ioffphi + i - 2 ] )
        * b11d[ ioffb11d + i - 1 ]
        + Math.sin( phi[ ioffphi + i - 2 ] )
        * b12e[ ioffb12e + i - 2 ];
      x2 = Math.cos( phi[ ioffphi + i - 2 ] ) * b11bulge
        + Math.sin( phi[ ioffphi + i - 2 ] ) * b12bulge;
      y1 = Math.cos( phi[ ioffphi + i - 2 ] )
        * b21d[ ioffb21d + i - 1 ]
        + Math.sin( phi[ ioffphi + i - 2 ] )
        * b22e[ ioffb22e + i - 2 ];
      y2 = Math.cos( phi[ ioffphi + i - 2 ] ) * b21bulge
        + Math.sin( phi[ ioffphi + i - 2 ] ) * b22bulge;
      theta[ iofftheta + i - 1 ] = Math.atan2(
        Math.sqrt( y1 * y1 + y2 * y2 ),
        Math.sqrt( x1 * x1 + x2 * x2 ) );
      restart11 = Math.pow( b11d[ ioffb11d + i - 1 ], 2 )
        + b11bulge * b11bulge <= thresh * thresh;
      restart12 = Math.pow( b12e[ ioffb12e + i - 2 ], 2 )
        + b12bulge * b12bulge <= thresh * thresh;
      restart21 = Math.pow( b21d[ ioffb21d + i - 1 ], 2 )
        + b21bulge * b21bulge <= thresh * thresh;
      restart22 = Math.pow( b22e[ ioffb22e + i - 2 ], 2 )
        + b22bulge * b22bulge <= thresh * thresh;
      if ( ! restart11 && ! restart12 ) {
        cs.setValue( work[ ioffwork + iu1sn + i - 2 ] );
        sn.setValue( work[ ioffwork + iu1cs + i - 2 ] );
        LaPack1.dlartgp( x2, x1, cs, sn, r );
        work[ ioffwork + iu1sn + i - 2 ] = cs.getValue();
        work[ ioffwork + iu1cs + i - 2 ] = sn.getValue();
      } else if ( ! restart11 && restart12 ) {
        cs.setValue( work[ ioffwork + iu1sn + i - 2 ] );
        sn.setValue( work[ ioffwork + iu1cs + i - 2 ] );
        LaPack1.dlartgp( b11bulge, b11d[ ioffb11d + i - 1 ],
          cs, sn, r );
        work[ ioffwork + iu1sn + i - 2 ] = cs.getValue();
        work[ ioffwork + iu1cs + i - 2 ] = sn.getValue();
      } else if ( restart11 && ! restart12 ) {
        cs.setValue( work[ ioffwork + iu1sn + i - 2 ] );
        sn.setValue( work[ ioffwork + iu1cs + i - 2 ] );
        LaPack1.dlartgp( b12bulge, b12e[ ioffb12e + i - 2 ],
          cs, sn, r );
        work[ ioffwork + iu1sn + i - 2 ] = cs.getValue();
        work[ ioffwork + iu1cs + i - 2 ] = sn.getValue();
      } else if ( mu <= nu ) {
        cs.setValue( work[ ioffwork + iu1cs + i - 2 ] );
        sn.setValue( work[ ioffwork + iu1sn + i - 2 ] );
        LaPack2.dlartgs( b11e[ ioffb11e + i - 1 ],
          b11d[ ioffb11d + i ], mu, cs, sn );
        work[ ioffwork + iu1cs + i - 2 ] = cs.getValue();
        work[ ioffwork + iu1sn + i - 2 ] = sn.getValue();
      } else {
        cs.setValue( work[ ioffwork + iu1cs + i - 2 ] );
        sn.setValue( work[ ioffwork + iu1sn + i - 2 ] );
        LaPack2.dlartgs( b12d[ ioffb12d + i - 1 ],
          b12e[ ioffb12e + i - 1 ], nu, cs, sn );
        work[ ioffwork + iu1cs + i - 2 ] = cs.getValue();
        work[ ioffwork + iu1sn + i - 2 ] = sn.getValue();
      }
      if ( ! restart21 && ! restart22 ) {
        cs.setValue( work[ ioffwork + iu2sn + i - 2 ] );
        sn.setValue( work[ ioffwork + iu2cs + i - 2 ] );
        LaPack1.dlartgp( y2, y1, cs, sn, r );
        work[ ioffwork + iu2sn + i - 2 ] = cs.getValue();
        work[ ioffwork + iu2cs + i - 2 ] = sn.getValue();
      } else if ( ! restart21 && restart22 ) {
        cs.setValue( work[ ioffwork + iu2sn + i - 2 ] );
        sn.setValue( work[ ioffwork + iu2cs + i - 2 ] );
        LaPack1.dlartgp( b21bulge, b21d[ ioffb21d + i - 1 ],
          cs, sn, r );
        work[ ioffwork + iu2sn + i - 2 ] = cs.getValue();
        work[ ioffwork + iu2cs + i - 2 ] = sn.getValue();
      } else if ( restart21 && ! restart22 ) {
        cs.setValue( work[ ioffwork + iu2sn + i - 2 ] );
        sn.setValue( work[ ioffwork + iu2cs + i - 2 ] );
        LaPack1.dlartgp( b22bulge, b22e[ ioffb22e + i - 2 ],
          cs, sn, r );
        work[ ioffwork + iu2sn + i - 2 ] = cs.getValue();
        work[ ioffwork + iu2cs + i - 2 ] = sn.getValue();
      } else if ( nu < mu ) {
        cs.setValue( work[ ioffwork + iu2cs + i - 2 ] );
        sn.setValue( work[ ioffwork + iu2sn + i - 2 ] );
        LaPack2.dlartgs( b21e[ ioffb21e + i - 1 ],
          b21e[ ioffb21e + i ], nu, cs, sn );
        work[ ioffwork + iu2cs + i - 2 ] = cs.getValue();
        work[ ioffwork + iu2sn + i - 2 ] = sn.getValue();
      } else {
        cs.setValue( work[ ioffwork + iu2cs + i - 2 ] );
        sn.setValue( work[ ioffwork + iu2sn + i - 2 ] );
        LaPack2.dlartgs( b22d[ ioffb22d + i - 1 ],
          b22e[ ioffb22e + i - 1 ], mu, cs, sn );
        work[ ioffwork + iu2cs + i - 2 ] = cs.getValue();
        work[ ioffwork + iu2sn + i - 2 ] = sn.getValue();
      }
      work[ ioffwork + iu2cs + i - 2 ] =
        - work[ ioffwork + iu2cs + i - 2 ];
      work[ ioffwork + iu2sn + i - 2 ] =
        - work[ ioffwork + iu2sn + i - 2 ];
      temp = work[ ioffwork + iu1cs + i - 2 ]
        * b11e[ ioffb11e + i - 1 ]
        + work[ ioffwork + iu1sn + i - 2 ]
        * b11d[ ioffb11d + i ];
      b11d[ ioffb11d + i ] = work[ ioffwork + iu1cs + i - 2 ]
        * b11d[ ioffb11d + i ]
        - work[ ioffwork + iu1sn + i - 2 ]
        * b11e[ ioffb11e + i - 1 ];
      b11e[ ioffb11e + i - 1 ] = temp;
      if ( i < imax - 1 ) {
        b11bulge = work[ ioffwork + iu1sn + i - 2 ]
          * b11e[ ioffb11e + i ];
        b11e[ ioffb11e + i ] *= work[ ioffwork + iu1cs + i - 2 ];
      }
      temp = work[ ioffwork + iu2cs + i - 2 ]
        * b21e[ ioffb21e + i - 1 ]
        + work[ ioffwork + iu2sn + i - 2 ]
        * b21d[ ioffb21d + i ];
      b21d[ ioffb21d + i ] = work[ ioffwork + iu2cs + i - 2 ]
        * b21d[ ioffb21d + i ]
        - work[ ioffwork + iu2sn + i - 2 ]
        * b21e[ ioffb21e + i - 1 ];
      b21e[ ioffb21e + i - 1 ] = temp;
      if ( i < imax - 1 ) {
        b21bulge = work[ ioffwork + iu2sn + i - 2 ]
          * b21e[ ioffb21e + i ];
        b21e[ ioffb21e + i ] = work[ ioffwork + iu2cs + i - 2 ]
          * b21e[ ioffb21e + i ];
      }
      temp = work[ ioffwork + iu1cs + i - 2 ]
        * b12d[ ioffb12d + i - 1 ]
        + work[ ioffwork + iu1sn + i - 2 ]
        * b12e[ ioffb12e + i - 1 ];
      b12e[ ioffb12e + i - 1 ] = work[ ioffwork + iu1cs + i - 2 ]
        * b12e[ ioffb12e + i - 1 ]
        - work[ ioffwork + iu1sn + i - 2 ]
        * b12d[ ioffb12d + i - 1 ];
      b12d[ ioffb12d + i - 1 ] = temp;
      b12bulge = work[ ioffwork + iu1sn + i - 2 ]
        * b12d[ ioffb12d + i ];
      b12d[ ioffb12d + i ] = work[ ioffwork + iu1cs + i - 2 ]
        * b12d[ ioffb12d + i ];
      temp = work[ ioffwork + iu2cs + i - 2 ]
        * b22d[ ioffb21d + i - 1 ]
        + work[ ioffwork + iu2sn + i - 2 ]
        * b22e[ ioffb22e + i - 1 ];
      b22e[ ioffb22e + i - 1 ] = work[ ioffwork + iu2cs + i - 2 ]
        * b22e[ ioffb22e + i - 1 ]
        - work[ ioffwork + iu2sn + i - 2 ]
        * b22d[ ioffb22d + i - 1 ];
      b22d[ ioffb22d + i - 1 ] = temp;
      b22bulge = work[ ioffwork + iu2sn + i - 2 ]
        * b22d[ ioffb22d + i ];
      b22d[ ioffb22d + i ] = work[ ioffwork + iu2cs + i - 2 ]
        * b22d[ ioffb22d + i ];
    }
    x1 = Math.sin( theta[ iofftheta + imax - 2 ] )
      * b11e[ ioffb11e + imax - 2 ]
      + Math.cos( theta[ iofftheta + imax - 2 ] )
      * b21e[ ioffb21e + imax - 2 ];
    y1 = Math.sin( theta[ iofftheta + imax - 2 ] )
      * b12d[ ioffb12d + imax - 2 ]
      + Math.cos( theta[ iofftheta + imax - 2 ] )
      * b22d[ ioffb22d + imax - 2 ];
    y2 = Math.sin( theta[ iofftheta + imax - 2 ] ) * b12bulge
      + Math.cos( theta[ iofftheta + imax - 2 ] ) * b22bulge;
    phi[ ioffphi + imax - 2 ] = Math.atan2(
      Math.abs( x1 ), Math.sqrt( y1 * y1 + y2 * y2 ) );
    restart12 = Math.pow( b12d[ ioffb12d + imax - 2 ], 2 )
      + b12bulge * b12bulge <= thresh * thresh;
    restart22 = Math.pow( b22d[ ioffb22d + imax - 2 ], 2 )
      + b22bulge * b22bulge <= thresh * thresh;
    if ( ! restart12 && ! restart22 ) {
      cs.setValue( work[ ioffwork + iv2tsn + imax - 3 ] );
      sn.setValue( work[ ioffwork + iv2tcs + imax - 3 ] );
      LaPack1.dlartgp( y2, y1, cs, sn, r );
      work[ ioffwork + iv2tsn + imax - 3 ] = cs.getValue();
      work[ ioffwork + iv2tcs + imax - 3 ] = sn.getValue();
    } else if ( ! restart12 && restart22 ) {
      cs.setValue( work[ ioffwork + iv2tsn + imax - 3 ] );
      sn.setValue( work[ ioffwork + iv2tcs + imax - 3 ] );
      LaPack1.dlartgp( b12bulge, b12d[ ioffb12d + imax - 2 ],
        cs, sn, r );
      work[ ioffwork + iv2tsn + imax - 3 ] = cs.getValue();
      work[ ioffwork + iv2tcs + imax - 3 ] = sn.getValue();
    } else if ( restart12 && ! restart22 ) {
      cs.setValue( work[ ioffwork + iv2tsn + imax - 3 ] );
      sn.setValue( work[ ioffwork + iv2tcs + imax - 3 ] );
      LaPack1.dlartgp( b22bulge, b22d[ ioffb22d + imax - 2 ],
        cs, sn, r );
      work[ ioffwork + iv2tsn + imax - 3 ] = cs.getValue();
      work[ ioffwork + iv2tcs + imax - 3 ] = sn.getValue();
    } else if ( nu < mu ) {
      cs.setValue( work[ ioffwork + iv2tcs + imax - 3 ] );
      sn.setValue( work[ ioffwork + iv2tsn + imax - 3 ] );
      LaPack2.dlartgs( b12e[ ioffb12e + imax - 2 ],
        b12d[ ioffb12d + imax - 1 ], nu, cs, sn );
      work[ ioffwork + iv2tcs + imax - 3 ] = cs.getValue();
      work[ ioffwork + iv2tsn + imax - 3 ] = sn.getValue();
    } else {
      cs.setValue( work[ ioffwork + iv2tcs + imax - 3 ] );
      sn.setValue( work[ ioffwork + iv2tsn + imax - 3 ] );
      LaPack2.dlartgs( b22e[ ioffb22e + imax - 2 ],
        b22d[ ioffb22d + imax - 1 ], mu, cs, sn );
      work[ ioffwork + iv2tcs + imax - 3 ] = cs.getValue();
      work[ ioffwork + iv2tsn + imax - 3 ] = sn.getValue();
    }
    temp = work[ ioffwork + iv2tcs + imax - 3 ]
      * b12e[ ioffb12e + imax - 2 ]
      + work[ ioffwork + iv2tsn + imax - 3 ]
      * b12d[ ioffb12d + imax - 1 ];
    b12d[ ioffb12d + imax - 1 ] = work[ ioffwork + iv2tcs + imax - 3 ]
      * b12d[ ioffb12d + imax - 1 ]
      - work[ ioffwork + iv2tsn + imax - 3 ]
      * b12e[ ioffb12e + imax - 2 ];
    b12e[ ioffb12e + imax - 2 ] = temp;
    temp = work[ ioffwork + iv2tcs + imax - 3 ]
      * b22e[ ioffb21e + imax - 2 ]
      + work[ ioffwork + iv2tsn + imax - 3 ]
      * b22d[ ioffb22d + imax - 1 ];
    b22d[ ioffb22d + imax - 1 ] = work[ ioffwork + iv2tcs + imax - 3 ]
      * b22d[ ioffb22d + imax - 1 ]
      - work[ ioffwork + iv2tsn + imax - 3 ]
      * b22e[ ioffb22e + imax - 2 ];
    b22e[ ioffb22e + imax - 2 ] = temp;
    if ( wantu1 ) {
      if ( colmajor ) {
        LaPack0.dlasr( 'R', 'V', 'F', p, imax - imin + 1, work, work,
          U1, ldu1, ioffwork + iu1cs + imin - 2,
          ioffwork + iu1sn + imin - 2, ioffu1 + ( imin - 1 ) * ldu1 );
      } else {
        LaPack0.dlasr( 'L', 'V', 'F', imax - imin + 1, p, work, work,
          U1, ldu1, ioffwork + iu1cs + imin - 2,
          ioffwork + iu1sn + imin - 2, ioffu1 + imin - 1 );
      }
    }
    if ( wantu2 ) {
      if ( colmajor ) {
        LaPack0.dlasr( 'R', 'V', 'F', m - p, imax - imin + 1, work,
          work, U2, ldu2, ioffwork + iu2cs + imin - 2,
          ioffwork + iu2sn + imin - 2, ioffu2 + ( imin - 1 ) * ldu2 );
      } else {
        LaPack0.dlasr( 'L', 'V', 'F', imax - imin + 1, m - p, work,
          work, U2, ldu2, ioffwork + iu2cs + imin - 2,
          ioffwork + iu2sn + imin - 2, ioffu2 + imin - 1 );
      }
    }
    if ( wantv1t ) {
      if ( colmajor ) {
        LaPack0.dlasr( 'L', 'V', 'F', imax - imin + 1, q, work, work,
          V1t, ldv1t, ioffwork + iv1tcs + imin - 2,
          ioffwork + iv1tsn + imin - 2, ioffv1t + imin - 1 );
      } else {
        LaPack0.dlasr( 'R', 'V', 'F', q, imax - imin + 1,  work, work,
          V1t, ldv1t, ioffwork + iv1tcs + imin - 2,
          ioffwork + iv1tsn + imin - 2,
          ioffv1t + ( imin - 1 ) * ldv1t );
      }
    }
    if ( wantv2t ) {
      if ( colmajor ) {
        LaPack0.dlasr( 'L', 'V', 'F', imax - imin + 1, m - q, work,
          work, V2t, ldv2t, ioffwork + iv2tcs + imin - 2,
          ioffwork + iv2tsn + imin - 2, ioffv2t + imin - 1 );
      } else {
        LaPack0.dlasr( 'R', 'V', 'F', m - q, imax - imin + 1, work,
          work, V2t, ldv2t, ioffwork + iv2tcs + imin - 2,
          ioffwork + iv2tsn + imin - 2,
          ioffv2t + ( imin - 1 ) * ldv2t );
      }
    }
    if ( b11e[ ioffb11e + imax - 2 ] +
    b21e[ ioffb21e + imax - 2 ] > 0. ) {
      b11d[ ioffb11d + imax - 1 ] = - b11d[ ioffb11d + imax - 1 ];
      b21d[ ioffb21d + imax - 1 ] = - b21d[ ioffb21d + imax - 1 ];
      if ( wantv1t ) {
        if ( colmajor ) {
          Blas1.dscal( q, -1., V1t, ldv1t, ioffv1t + imax - 1 );
        } else {
          Blas1.dscal( q, -1., V1t, 1,
            ioffv1t + ( imax - 1 ) * ldv1t );
        }
      }
    }
    x1 = Math.cos( phi[ ioffphi + imax - 2 ] )
      * b11d[ ioffb11d + imax - 1 ]
      + Math.sin( phi[ ioffphi + imax - 2 ] )
      * b12e[ ioffb12e + imax - 2 ];
    y1 = Math.cos( phi[ ioffphi + imax - 2 ] )
      * b21d[ ioffb21d + imax - 1 ]
      + Math.sin( phi[ ioffphi + imax - 2 ] )
      * b22e[ ioffb22e + imax - 2 ];
    theta[ iofftheta + imax - 1 ] =
      Math.atan2( Math.abs( y1 ), Math.abs( x1 ) );
    if ( b11d[ ioffb11d + imax - 1 ] +
    b12e[ ioffb12e + imax - 2 ] < 0. ) {
      b12d[ ioffb12d + imax - 1 ] = - b12d[ ioffb12d + imax - 1 ];
      if ( wantu1 ) {
        if ( colmajor ) {
          Blas1.dscal( p, -1., U1, 1, ioffu1 + ( imax - 1 ) * ldu1 );
        } else {
          Blas1.dscal( p, -1., U1, ldu1, ioffu1 + imax - 1 );
        }
      }
    }
    if ( b21d[ ioffb21d + imax - 1 ] +
    b22e[ ioffb22e + imax - 2 ] > 0. ) {
      b22d[ ioffb22d + imax - 1 ] = - b22d[ ioffb22d + imax - 1 ];
      if ( wantu2 ) {
        if ( colmajor ) {
          Blas1.dscal( m - p, -1., U2, 1,
            ioffu2 + ( imax - 1 ) * ldu2 );
        } else {
          Blas1.dscal( m - p, -1., U2, ldu2, ioffu2 + imax - 1 );
        }
      }
    }
    if ( b12d[ ioffb12d + imax - 1 ] +
    b22d[ ioffb22d + imax - 1 ] < 0. ) {
      if ( wantv2t ) {
        if ( colmajor ) {
          Blas1.dscal( m - q, -1., V2t, ldv2t, ioffv2t + imax - 1 );
        } else {
          Blas1.dscal( m - q, -1., V2t, 1,
            ioffv2t + ( imax - 1 ) * ldv2t );
        }
      }
    }
    for ( i = imin; i <= imax; i ++ ) {
      if ( theta[ iofftheta + i - 1 ] < thresh ) {
        theta[ iofftheta + i - 1 ] = 0.;
      } else if ( theta[ iofftheta + i - 1 ] > piover2 - thresh ) {
        theta[ iofftheta + i - 1 ] = piover2;
      }
    }
    for ( i = imin; i <= imax - 1; i ++ ) {
      if ( phi[ ioffphi + i - 1 ] < thresh ) {
        phi[ ioffphi + i - 1 ] = 0.;
      } else if ( phi[ ioffphi + i - 1 ] > piover2 - thresh ) {
        phi[ ioffphi + i - 1 ] = piover2;
      }
    }
    if ( imax > 1 ) {
      while ( phi[ ioffphi + imax - 2 ] == 0. ) {
        imax --;
        if ( imax <= 1 ) break;
      }
    }
    if ( imin > imax - 1 ) imin = imax - 1;
    if ( imin > 1 ) {
      while ( phi[ ioffphi + imin - 1 ] != 0. ) {
        imin --;
        if ( imin <= 1 ) break;
      }
    }
  }
  for ( i = 1; i <= q; i ++ ) {
    var mini = i;
    thetamin = theta[ iofftheta + i - 1 ];
    for ( var j = i + 1; j <= q; j ++ ) {
      if ( theta[ iofftheta + j - 1 ] < thetamin ) {
        mini = j;
        thetamin = theta[ iofftheta + j - 1 ];
      }
    }
    if ( mini != i ) {
      theta[ iofftheta + mini - 1 ] = theta[ iofftheta + i - 1 ];
      theta[ iofftheta + i - 1 ] = thetamin;
      if ( colmajor ) {
        if ( wantu1 ) {
          Blas1.dswap( p, U1, 1, U1, 1, ioffu1 + ( i - 1 ) * ldu1,
            ioffu1 + ( mini - 1 ) * ldu1 );
        }
        if ( wantu2 ) {
          Blas1.dswap( m - p, U2, 1, U2, 1, ioffu2 + ( i - 1 ) * ldu2,
            ioffu2 + ( mini - 1 ) * ldu2 );
        }
        if ( wantv1t ) {
          Blas1.dswap( q, V1t, ldv1t, V1t, ldv1t, ioffv1t + i - 1,
            ioffv1t + mini - 1 );
        }
        if ( wantv2t ) {
          Blas1.dswap( m - q, V2t, ldv2t, V2t, ldv2t, ioffv2t + i - 1,
            ioffv2t + mini - 1 );
        }
      } else {
        if ( wantu1 ) {
          Blas1.dswap( p, U1, ldu1, U1, ldu1, ioffu1 + i - 1,
            ioffu1 + mini - 1 );
        }
        if ( wantu2 ) {
          Blas1.dswap( m - p, U2, ldu2, U2, ldu2,
            ioffu2 + i - 1, ioffu2 + mini - 1 );
        }
        if ( wantv1t ) {
          Blas1.dswap( q, V1t, 1, V1t, 1, ioffv1t + ( i - 1 ) * ldv1t,
            ioffv1t + ( mini - 1 ) * ldv1t );
        }
        if ( wantv2t ) {
          Blas1.dswap( m - q, V2t, 1, V2t, 1,
            ioffv2t + ( i - 1 ) * ldv2t,
            ioffv2t + ( mini - 1 ) * ldv2t );
        }
      }
    }
  }
}
//*************************************************************************
LaPack3.dgebrd = function( m, n, A, lda, d, e, tauq, taup, work,
lwork, info, ioffa, ioffd, ioffe, iofftauq, iofftaup, ioffwork ) {
  throw new Error("not tested");
  info.setValue( 0 );
  var nb = Math.max( 1, 
    LaPack0.ilaenv( 1, 'dgebrd', ' ', m, n, -1, -1 ) );
  var lwkopt = ( m + n ) * nb;
  work[ ioffwork ] = Number( lwkopt );
  var lquery = ( lwork == -1 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -4 );
  else if ( lwork < Math.max( 1, Math.max( m, n ) ) && ! lquery ) {
    info.setValue( -10 );
  }
  if ( info.getValue() < 0 ) {
    Blas2.xerbla( 'dgebrd', - info.getValue() );
    return;
  } else if ( lquery ) return;
  var minmn = Math.min( m, n );
  if ( minmn === 0 ) {
    work[ ioffwork ] = 1;
    return;
  }
  var ws = Math.max( m, n );
  var ldwrkx = m;
  var ldwrky = n;
  if ( nb > 1 && nb < minmn ) {
    var nx = Math.max( nb,
      LaPack0.ilaenv( 3, 'dgebrd', ' ', m, n, -1, -1 ) );
    if ( nx < minmn ) {
      ws = ( m + n ) * nb;
      if ( lwork < ws ) {
        var nbmin =
          LaPack0.ilaenv( 2, 'dgebrd', ' ', m, n, -1, -1 );
        if ( lwork >= ( m + n ) * nbmin ) nb = lwork / ( m + n );
        else {
          nb = 1;
          nx = minmn;
        }
      }
    }
  } else nx = minmn;
  for ( var i = 1; i <= minmn - nx; i += nb ) {
    LaPack2.dlabrd( m - i + 1, n - i + 1, nb, A, lda, d, e, tauq, taup,
      work, ldwrkx, work, ldwrky, ioffa + i - 1 + ( i - 1 ) * lda,
      ioffd + i - 1, ioffe + i - 1, iofftauq + i - 1, iofftaup + i - 1,
      ioffwork, ioffwork + ldwrkx * nb );
    Blas3.dgemm( 'No transpose', 'Transpose', m - i - nb + 1,
      n - i - nb + 1, nb, -1., A, lda, work, ldwrky, 1., A, lda,
      ioffa + i + nb - 1 + ( i - 1 ) * lda,
      ioffwork + ldwrkx * nb + nb,
      ioffa + i + nb - 1 + ( i + nb - 1 ) * lda );
    Blas3.dgemm( 'No transpose', 'No transpose', m - i - nb + 1,
      n - i - nb + 1, nb, -1., work, ldwrkx, A, lda, 1., A, lda,
      ioffwork + nb, ioffa + i - 1 + ( i + nb - 1 ) * lda,
      ioffa + i + nb - 1 + ( i + nb - 1 ) * lda );
    if ( m >= n ) {
      for ( var j = i; j <= i + nb - 1; j ++ ) {
        A[ ioffa + j - 1 + ( j - 1 ) * lda ] = d[ ioffd + j - 1 ];
        A[ ioffa + j - 1 + j * lda ] = e[ ioffe + j - 1 ];
      }
    } else {
      for ( j = i; j <= i + nb - 1; j ++ ) {
        A[ ioffa + j - 1 + ( j - 1 ) * lda ] = d[ ioffd + j - 1 ];
        A[ ioffa + j + ( j - 1 ) * lda ] = e[ ioffe + j - 1 ];
      }
    }
  }
  var iinfo = new IntReference();
  LaPack2.dgebd2( m - i + 1, n - i + 1, A, lda, d, e, tauq, taup,
    work, iinfo, ioffa + i - 1 + ( i - 1 ) * lda, ioffd + i - 1,
    ioffe + i - 1, iofftauq + i - 1, iofftaup + i - 1, ioffwork );
  work[ ioffwork ] = ws;
}
//*************************************************************************
LaPack3.dgehrd = function( n, ilo, ihi, A, lda, tau, work,
lwork, info, ioffa, iofftau, ioffwork ) {
//document.getElementById("debug_textarea").value +=
//  "entering LaPack3.dgehrd, n,ilo,ihi,lda,lwork = " + n + " " + ilo + " "
//  + ihi + " " + lda + " " + lwork + "\n";
//document.getElementById("debug_textarea").value += "A = \n";
//for ( var ii = 0; ii < n; ii ++ ) {
//  for ( var jj = 0; jj < n; jj ++ ) {
//    var Aij = A[ ioffa + ii + jj * lda ];
//    document.getElementById("debug_textarea").value += Aij + " ";
//  }
//  document.getElementById("debug_textarea").value += "\n";
//}
  var nbmax = 64;
  var ldt = nbmax + 1;
  var T = new Array( ldt * nbmax );
  info.setValue( 0 );
  var nb = Math.min( nbmax,
    LaPack0.ilaenv( 1, 'dgehrd', ' ', n, ilo, ihi, -1 ) );
  var lwkopt = n * nb;
  work[ ioffwork ] = lwkopt;
  var lquery = ( lwork == -1 );
  if ( n < 0 ) info.setValue( -1 );
  else if ( ilo < 1 || ilo > Math.max( 1, n ) ) info.setValue( -2 );
  else if ( ihi < Math.min( ilo, n ) || ihi > n ) info.setValue( -3 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -5 );
  else if ( lwork < Math.max( 1, n ) && ! lquery ) info.setValue( -8 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgehrd', - info.getValue() );
    return;
  } else if ( lquery ) return;
  for ( var i = 1; i <= ilo - 1; i ++ ) {
    tau[ iofftau + i - 1 ] = 0.;
  }
  for ( i = Math.max( 1, ihi ); i <= n - 1; i ++ ) {
    tau[ iofftau + i - 1 ] = 0.;
  }
  var nh = ihi - ilo + 1;
  if ( nh <= 1 ) {
    work[ ioffwork ] = 1;
    return;
  }
  nb = Math.min( nbmax,
    LaPack0.ilaenv( 1, 'dgehrd', ' ', n, ilo, ihi, -1 ) );
  var nbmin = 2;
  var iws = 1;
  if ( nb > 1 && nb < nh ) {
    var nx = Math.max( nb,
      LaPack0.ilaenv( 3, 'dgehrd', ' ', n, ilo, ihi, -1 ) );
    if ( nx < nh ) {
      iws = n * nb;
      if ( lwork < iws ) {
        nbmin = Math.max( 2,
          LaPack0.ilaenv( 2, 'dgehrd', ' ', n, ilo, ihi, -1 ) );
        nb = ( lwork >= n * nbmin ? lwork / n : 1 );
      }
    }
  }
  var ldwork = n;
  if ( nb < nbmin || nb >= nh ) i = ilo;
  else {
    for ( i = ilo; i <= ihi - 1 - nx; i += nb ) {
      var ib = Math.min( nb, ihi - i );
      LaPack2.dlahr2( ihi, i, ib, A, lda, tau, T, ldt, work, ldwork,
        ioffa + ( i - 1 ) * lda, iofftau + i - 1, 0, ioffwork );
      var ei = A[ ioffa + i + ib - 1 + ( i + ib - 2 ) * lda ];
      A[ ioffa + i + ib - 1 + ( i + ib - 2 ) * lda ] = 1.;
      Blas3.dgemm( 'No transpose', 'Transpose', ihi, ihi - i - ib + 1,
        ib, -1., work, ldwork, A, lda, 1., A, lda, ioffwork,
        ioffa + i + ib - 1 + ( i - 1 ) * lda,
        ioffa + ( i + ib - 1 ) * lda );
      A[ ioffa + i + ib - 1 + ( i + ib - 2 ) * lda ] = ei;
      Blas3.dtrmm( 'Right', 'Lower', 'Transpose', 'Unit', i, ib - 1,
        1., A, lda, work, ldwork, ioffa + i + ( i - 1 ) * lda,
        ioffwork );
      for ( var j = 0; j <= ib - 2; j ++ ) {
        Blas1.daxpy( i, -1., work, 1, A, 1, ioffwork + ldwork * j,
          ioffa + ( i + j ) * lda );
      }
      LaPack1.dlarfb( 'Left', 'Transpose', 'Forward', 'Columnwise',
        ihi - i, n - i - ib + 1, ib, A, lda, T, ldt, A, lda, work,
        ldwork, ioffa + i + ( i - 1 ) * lda, 0,
        ioffa + i + ( i + ib - 1 ) * lda, ioffwork );
    }
  }
  var iinfo = new IntReference();
  LaPack2.dgehd2( n, i, ihi, A, lda, tau, work, iinfo, ioffa, iofftau,
    ioffwork );
  work[ ioffwork ] = iws;
//document.getElementById("debug_textarea").value +=
//  "leaving LaPack3.dgehrd\n";
//document.getElementById("debug_textarea").value += "A = \n";
//for ( var ii = 0; ii < n; ii ++ ) {
//  for ( var jj = 0; jj < n; jj ++ ) {
//    var Aij = A[ ioffa + ii + jj * lda ];
//    document.getElementById("debug_textarea").value += Aij + " ";
//  }
//  document.getElementById("debug_textarea").value += "\n";
//}
}
LaPack3.zgehrd = function( n, ilo, ihi, A, lda, tau, work,
lwork, info ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack3.dgelqf = function( m, n, A, lda, tau, work, lwork, info,
ioffa, iofftau, ioffwork ) {
  info.setValue( 0 );
  var nb = LaPack0.ilaenv( 1, 'dgelqf', ' ', m, n, -1, -1 );
  var lwkopt = m * nb;
  work[ ioffwork ] = lwkopt;
  var lquery = ( lwork == -1 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -4 );
  else if ( lwork < Math.max( 1, m ) && ! lquery ) info.setValue( -7 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgelqf', - info.getValue() );
    return;
  } else if ( lquery ) return;
  var k = Math.min( m, n );
  if ( k == 0 ) {
    work[ ioffwork ] = 1;
    return;
  }
  var nbmin = 2;
  var nx = 0;
  var iws = m;
  if ( nb > 1 && nb < k ) {
    nx = Math.max( 0,
      LaPack0.ilaenv( 3, 'dgelqf', ' ', m, n, -1, -1 ) );
    if ( nx < k ) {
      var ldwork = m;
      iws = ldwork * nb;
      if ( lwork < iws ) {
        nb = lwork / ldwork ;
        nbmin = Math.max( 2, 
          LaPack0.ilaenv( 2, 'dgelqf', ' ', m, n, -1, -1 ) );
      }
    }
  }
  var iinfo = new IntReference();
  if ( nb >= nbmin && nb < k && nx < k ) {
    for ( var i = 1; i <= k - nx; i += nb ) {
      var ib = Math.min( k - i + 1, nb );
      LaPack2.dgelq2( ib, n - i + 1, A, lda, tau, work, iinfo,
        ioffa + i - 1 + ( i - 1 ) * lda, iofftau + i - 1, ioffwork );
      if ( i + ib <= m ) {
        LaPack0.dlarft( 'Forward', 'Rowwise', n - i + 1, ib, A, lda,
          tau, work, ldwork, ioffa + i - 1 + ( i - 1 ) * lda,
          iofftau + i - 1, ioffwork );
        LaPack1.dlarfb( 'Right', 'No transpose', 'Forward', 'Rowwise',
          m - i - ib + 1, n - i + 1, ib, A, lda, work, ldwork, A, lda,
          work, ldwork, ioffa + i - 1 + ( i - 1 ) * lda, ioffwork,
          ioffa + i + ib - 1 + ( i - 1 ) * lda, ioffwork + ib );
      }
    }
  } else i = 1;
  if ( i <= k ) {
    LaPack2.dgelq2( m - i + 1, n - i + 1, A, lda, tau, work, iinfo,
      ioffa + i - 1 + ( i - 1 ) * lda, iofftau + i - 1, ioffwork );
  }
  work[ ioffwork ] = iws;
}
LaPack3.zgelqf = function( m, n, A, lda, tau, work, lwork,
info ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack3.dgeqlf = function( m, n, A, lda, tau, work, lwork, info,
ioffa, iofftau, ioffwork ) {
  info.setValue( 0 );
  var lquery = ( lwork == -1 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -4 );
  else if ( lwork < Math.max( 1, m ) && ! lquery ) info.setValue( -7 );
  if ( info.getValue() == 0 ) {
    var k = Math.min( m, n );
    if ( k == 0 ) var lwkopt = 1;
    else {
      var nb = LaPack0.ilaenv( 1, 'dgeqlf', ' ', m, n, -1, -1 );
      lwkopt = n * nb;
    }
    work[ ioffwork ] = lwkopt;
    if ( lwork < Math.max( 1, n ) && ! lquery ) info.setValue( -7 );
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgeqlf', - info.getValue() );
    return;
  } else if ( lquery ) return;
  if ( k == 0 ) return;
  var nbmin = 2;
  var nx = 1;
  var iws = n;
  if ( nb > 1 && nb < k ) {
    nx = Math.max( 0,
      LaPack0.ilaenv( 3, 'dgeqlf', ' ', m, n, -1, -1 ) );
    if ( nx < k ) {
      var ldwork = n;
      iws = ldwork * nb;
      if ( lwork < iws ) {
        nb = lwork / ldwork ;
        nbmin = Math.max( 2, 
          LaPack0.ilaenv( 2, 'dgeqlf', ' ', m, n, -1, -1 ) );
      }
    }
  }
  var iinfo = new IntReference();
  if ( nb >= nbmin && nb < k && nx < k ) {
    var ki = ( ( k - nx - 1 ) / nb ) * nb;
    var kk = Math.min( k, ki + nb );
    for ( var i = k - kk + ki + 1; i >= k - kk + 1; i -= nb ) {
      var ib = Math.min( k - i + 1, nb );
      LaPack2.dgeql2( m - k + i + ib - 1, ib, A, lda, tau, work, iinfo,
        ioffa + ( n - k + i - 1 ) * lda, iofftau + i - 1, ioffwork );
      if ( n - k + i > 1 ) {
        LaPack0.dlarft( 'Backward', 'Columnwise', m - k + i + ib - 1,
          ib, A, lda, tau, work, ldwork,
          ioffa + ( n - k + i - 1 ) * lda,
          iofftau + i - 1, ioffwork );
        LaPack1.dlarfb( 'Left', 'Transpose', 'Backward', 'Columnwise',
          m - k + i + ib - 1, n - k + i - 1, ib, A, lda, work, ldwork,
          A, lda, work, ldwork, ioffa + ( n - k + i - 1 ) * lda,
          ioffwork, ioffa, ioffwork + ib );
      }
    }
    var mu = m - k + i + nb - 1;
    var nu = n - k + i + nb - 1;
  } else {
    mu = m;
    nu = n;
  }
  if ( mu > 0 && nu > 0 ) {
    LaPack2.dgeql2( mu, nu, A, lda, tau, work, iinfo,
      ioffa, iofftau, ioffwork );
  }
  work[ ioffwork ] = iws;
}
//*************************************************************************
LaPack3.dgeqpf = function( m, n, A, lda, jpvt, tau, work, info,
ioffa, ioffjpvt, iofftau, ioffwork ) {
  throw new Error("not tested");
  info.setValue( 0 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -4 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgeqpf', - info.getValue() );
    return;
  }
  var mn = Math.min( m, n );
  var itemp = 1;
  for ( var i = 1; i <= n; i ++ ) {
    if ( jpvt[ ioffjpvt + i - 1 ] != 0 ) {
      if ( i != itemp ) {
        Blas1.dswap( m, A, 1, A, 1, ioffa + ( i - 1 ) * lda,
          ioffa + ( itemp - 1 ) * lda );
        jpvt[ ioffjpvt + i - 1 ] = jpvt[ ioffjpvt + itemp - 1 ];
        jpvt[ ioffjpvt + itemp - 1 ] = i;
      } else jpvt[ ioffjpvt + i - 1 ] = i;
      itemp ++;
    } else jpvt[ ioffjpvt + i - 1 ] = i;
  }
  itemp --;
  if ( itemp > 0 ) {
    var ma = Math.min( itemp, m );
    LaPack2.dgeqr2( m, ma, A, lda, tau, work, info, ioffa, iofftau,
      ioffwork );
    if ( ma < n ) {
      LaPack2.dorm2r( 'Left', 'Transpose', m, n - ma, ma, A, lda, tau,
        A, lda, work, info, ioffa, iofftau, ioffa + ma * lda,
        ioffwork );
    }
  }
  if ( itemp < mn ) {
    for ( i = itemp + 1; i <= n; i ++ ) {
      work[ ioffwork + i - 1 ] = Blas1.dnrm2( m - itemp, A, 1,
        ioffa + itemp + ( i - 1 ) * lda );
      work[ ioffwork + n + i - 1 ] = work[ ioffwork + i - 1 ];
    }
    for ( i = itemp + 1; i <= mn; i ++ ) {
      var pvt = ( i - 1 )
        + Blas1.idamax( n - i + 1, work, 1, ioffwork + i - 1 );
      if ( pvt != i ) {
        Blas1.dswap( m, A, 1, A, 1, ioffa + ( pvt - 1 ) * lda,
          ioffa + ( i - 1 ) * lda );
        itemp = jpvt[ ioffjpvt + pvt - 1 ];
        jpvt[ ioffjpvt + pvt - 1 ] = jpvt[ ioffjpvt + i - 1 ];
        jpvt[ ioffjpvt + i - 1 ] = itemp;
        work[ ioffwork + pvt - 1 ] = work[ ioffwork + i - 1 ];
        work[ ioffwork + n + pvt - 1 ] = work[ ioffwork + n + i - 1 ];
      }
      if ( i < m ) {
        var alphaReference =
          new NumberReference( A[ ioffa + i - 1 + ( i - 1 ) * lda ] );
        var taurefReference =
          new NumberReference( tau[ iofftau + i - 1 ] );
        LaPack1.dlarfg( m - i + 1, alpha, A, 1, tauref,
          ioffa + i + ( i - 1 ) * lda );
        A[ ioffa + i - 1 + ( i - 1 ) * lda ] = alpha.getValue();
        tau[ iofftau + i - 1 ] = tauref.getValue();
      } else {
        alpha.setValue( A[ ioffa + m - 1 + ( m - 1 ) * lda ] );
        tauref.setValue( tau[ iofftau + m - 1 ] );
        LaPack1.dlarfg( 1, alpha, A, 1, tauref,
          ioffa + m - 1 + ( m - 1 ) * lda );
        A[ ioffa + m - 1 + ( m - 1 ) * lda ] = alpha.getValue();
        tau[ iofftau + m - 1 ] = tauref.getValue();
      }
      if ( i < n ) {
        var aii = A[ ioffa + i - 1 + ( i - 1 ) * lda ];
        A[ ioffa + i - 1 + ( i - 1 ) * lda ] = 1.;
        LaPack1.dlarf( 'Left', m - i + 1, n - i, A, 1, tau[i - 1 ], A,
          lda, work, ioffa + i - 1 + ( i - 1 ) * lda,
          ioffa + i - 1 + i * lda, ioffwork + 2 * n );
        A[ ioffa + i - 1 + ( i - 1 ) * lda ] = aii;
      }
      for ( var j = i + 1; j <= n; j ++ ) {
        if ( work[ ioffwork + j - 1 ] != 0. ) {
          var temp = 1.  - Math.pow( Math.abs(
            A[ ioffa + i - 1 + ( j - 1 ) * lda ] )
            / work[ ioffwork + j - 1 ], 2 )
          temp = Math.max( temp, 0. );
          var temp2 = 1. + 0.05 * temp
            * Math.pow( work[ ioffwork + j - 1 ]
            / work[ ioffwork + n + j - 1 ], 2 );
          if ( temp2 == 1. ) {
            if ( m - i > 0 ) {
              work[ ioffwork + j - 1 ] = Blas1.dnrm2( m - i, A, 1,
                ioffa + i + ( j - 1 ) * lda );
              work[ ioffwork + n + j - 1 ] = work[ ioffwork + j - 1 ];
            } else {
              work[ ioffwork + j - 1 ] = 0.;
              work[ ioffwork + n + j - 1 ] = 0.;
            }
          } else work[ ioffwork + j - 1 ] *= Math.sqrt( temp );
        }
      }
    }
  }
}
LaPack3.zgeqpf = function( m, n, A, lda, jpvt, tau, work, info )
{
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack3.dgeqrf = function( m, n, A, lda, tau, work, lwork, info,
ioffa, iofftau, ioffwork ) {
  info.setValue( 0 );
  var nb = LaPack0.ilaenv( 1, 'dgeqrf', ' ', m, n, -1, -1 );
  var lwkopt = n * nb;
  work[ ioffwork ] = lwkopt;
  var lquery = ( lwork == -1 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -4 );
  else if ( lwork < Math.max( 1, n ) && ! lquery ) info.setValue( -7 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgeqrf', - info.getValue() );
    return;
  } else if ( lquery ) return;
  var k = Math.min( m, n );
  if ( k == 0 ) {
    work[ ioffwork ] = 1;
    return;
  }
  var nbmin = 2;
  var nx = 0;
  var iws = n;
  if ( nb > 1 && nb < k ) {
    nx = Math.max( 0,
      LaPack0.ilaenv( 3, 'dgeqrf', ' ', m, n, -1, -1 ) );
    if ( nx < k ) {
      var ldwork = n;
      iws = ldwork * nb;
      if ( lwork < iws ) {
        nb = lwork / ldwork;
        nbmin = Math.max( 2,
          LaPack0.ilaenv( 2, 'dgeqrf', ' ', m, n, -1, -1 ) );
      }
    }
  }
  var iinfo = new IntReference();
  if ( nb >= nbmin && nb < k && nx < k ) {
    for ( var i = 1; i <= k - nx; i += nb ) {
      var ib = Math.min( k - i + 1, nb );
      LaPack2.dgeqr2( m - i + 1, ib, A, lda, tau, work, iinfo,
        ioffa + i - 1 + ( i - 1 ) * lda, iofftau + i - 1, ioffwork );
      if ( i + ib <= n ) {
        LaPack0.dlarft( 'Forward', 'Columnwise', m - i + 1, ib, 
          A, lda, tau, work, ldwork, ioffa + i - 1 + ( i - 1 ) * lda, 
          iofftau + i - 1, ioffwork );
        LaPack1.dlarfb( 'Left', 'Transpose', 'Forward', 'Columnwise',
          m - i + 1, n - i - ib + 1, ib, A, lda, work, ldwork, A, lda,
          work, ldwork, ioffa + i - 1 + ( i - 1 ) * lda, ioffwork,
          ioffa + i - 1 + ( i + ib - 1 ) * lda, ioffwork + ib );
      }
    }
  } else i = 1;
  if ( i <= k ) {
    LaPack2.dgeqr2( m - i + 1, n - i + 1, A, lda, tau, work, iinfo,
      ioffa + i - 1 + ( i - 1 ) * lda, iofftau + i - 1, ioffwork );
  }
  work[ ioffwork ] = iws;
}
LaPack3.zgeqrf = function( m, n, A, lda, tau, work, lwork,
info ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack3.dgeqrt = function( m, n, nb, A, lda, T, ldt, work, info,
ioffa, iofft, ioffwork ) {
  throw new Error("not programmed: WY storage");
}
//*************************************************************************
LaPack3.dgerqf = function( m, n, A, lda, tau, work, lwork, info,
ioffa, iofftau, ioffwork ) {
  info.setValue( 0 );
  var lquery = ( lwork == -1 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -4 );
  if ( info.getValue() == 0 ) {
    var k = Math.min( m, n );
    if ( k == 0 ) var lwkopt = 1;
    else {
      var nb = LaPack0.ilaenv( 1, 'dgerqf', ' ', m, n, -1, -1 );
      lwkopt = m * nb;
    }
    work[ ioffwork ] = lwkopt;
    if ( lwork < Math.max( 1, m ) && ! lquery ) info.setValue( -7 );
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgerqf', - info.getValue() );
    return;
  } else if ( lquery ) return;
  if ( k == 0 ) return;
  var nbmin = 2;
  var nx = 1;
  var iws = m;
  if ( nb > 1 && nb < k ) {
    nx = Math.max( 0,
      LaPack0.ilaenv( 3, 'dgerqf', ' ', m, n, -1, -1 ) );
    if ( nx < k ) {
      var ldwork = m;
      iws = ldwork * nb;
      if ( lwork < iws ) {
        nb = lwork / ldwork;
        nbmin = Math.max( 2,
          LaPack0.ilaenv( 2, 'dgerqf', ' ', m, n, -1, -1 ) );
      }
    }
  }
  var iinfo = new IntReference();
  if ( nb >= nbmin && nb < k && nx < k ) {
    var ki = ( ( k - nx - 1 ) / nb ) * nb;
    var kk = Math.min( k, ki + nb );
    for ( var i = k - kk + ki + 1; i >= k - kk + 1; i -= nb ) {
      var ib = Math.min( k - i + 1, nb );
      LaPack2.dgerq2( ib, n - k + i + ib - 1, A, lda, tau, work, iinfo,
        ioffa + m - k + i - 1, iofftau + i - 1, ioffwork );
      if ( m - k + i >= 1 ) {
        LaPack0.dlarft( 'Backward', 'Rowwise', n - k + i + ib - 1, ib, 
          A, lda, tau, work, ldwork,
          ioffa + m - k + i - 1, iofftau + i - 1, ioffwork );
        LaPack1.dlarfb( 'Right', 'No transpose', 'Backward', 'Rowwise',
          m - k + i - 1, n - k + i + ib - 1, ib, A, lda, work, ldwork,
          A, lda, work, ldwork, ioffa + m - k + i - 1, ioffwork,
          ioffa, ioffwork + ib );
      }
    }
    var mu = m - k + i + nb - 1;
    var nu = n - k + i + nb - 1;
  } else {
    mu = m;
    nu = n;
  }
  if ( mu > 0 && nu > 0 ) {
    LaPack2.dgerq2( mu, nu, A, lda, tau, work, iinfo, ioffa, iofftau,
      ioffwork );
  }
  work[ ioffwork ] = iws;
}
LaPack3.zgerqf = function( m, n, A, lda, tau, work, lwork,
info ) {
  throw new Error("not programmed: complex matrix");
}
//*************************************************************************
LaPack3.dgesvj = function( joba, jobu, jobv, m, n, A, lda, sva,
mv, V, ldv, work, lwork, info, ioffa, ioffsva, ioffv, ioffwork ) {
  var nsweep = 30;
  var fastr = new Array( 5 );
  var lsvec = ( jobu.charAt(0).toUpperCase() == 'U' );
  var uctol = ( jobu.charAt(0).toUpperCase() == 'C' );
  var rsvec = ( jobv.charAt(0).toUpperCase() == 'V' );
  var applv = ( jobv.charAt(0).toUpperCase() == 'A' );
  var upper = ( joba.charAt(0).toUpperCase() == 'U' );
  var lower = ( joba.charAt(0).toUpperCase() == 'L' );
  if ( ! ( upper || lower || joba.charAt(0).toUpperCase() == 'G' ) ) {
    info.setValue( -1 );
  } else if ( ! ( lsvec || uctol ||
  jobu.charAt(0).toUpperCase() == 'N' ) ) {
    info.setValue( -2 );
  } else if ( ! ( rsvec || applv ||
  jobv.charAt(0).toUpperCase() == 'N' ) ) {
    info.setValue( -3 );
  } else if ( m < 0 ) info.setValue( -4 );
  else if ( n < 0 || n > m ) info.setValue( -5 );
  else if ( lda < m ) info.setValue( 7 );
  else if ( mv < 0 ) info.setValue( -9 );
  else if ( ( rsvec && ldv < n ) || ( applv && ldv < mv ) ) {
    info.setValue( -11 );
  } else if ( uctol && work[ ioffwork ] <= 1. ) info.setValue( -12 );
  else if ( lwork < Math.max( m + n, 6 ) ) info.setValue( -13 );
  else info.getValue() = 0;
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dgesvj', - info.getValue() );
    return;
  }
  if ( m == 0 || n == 0 ) return;
  if ( uctol ) var ctol = work[ ioffwork ];
  else {
    ctol = ( lsvec || rsvec || applv ? Math.sqrt( Number( m ) )
      : Number( m ) );
  }
  var epsln = LaPack0.dlamch( 'Epsilon' );
  var rooteps = Math.sqrt( epsln );
  var sfmin = LaPack0.dlamch( 'SafeMinimum' );
  var rootsfmin = Math.sqrt( sfmin );
  var small = sfmin / epsln;
  var big = LaPack0.dlamch( 'Overflow' );
  var rootbig = 1. / rootsfmin;
  var large = big / Math.sqrt( Number( m * n ) );
  var bigtheta = 1. / rooteps;
  var tol = ctol * epsln;
  var roottol = Math.sqrt( tol );
  if ( Number( m ) * epsln >= 1. ) {
    info.setValue( -4 );
    Blas2.xerbla( 'dgesvj', - info.getValue() );
    return;
  }
  if ( rsvec ) {
    var mvl = n;
    LaPack0.dlaset( 'A', mvl, n, 0., 1., V, ldv, ioffv );
  } else if ( applv ) mvl = mv;
  rsvec = rsvec || applv;
  var skl = 1. / Math.sqrt( Number( m ) * Number( n ) );
  var noscale = true;
  var goscale = true;
  var aappReference = new NumberReference( );
  var aaqqReference = new NumberReference( );
  if ( lower ) {
    for ( var p = 1; p <= n; p ++ ) {
      aapp.setValue( 0. );
      aaqq.setValue( 1. );
      LaPack0.dlassq( m - p + 1, A, 1, aapp, aaqq,
        ioffa + p - 1 + ( p - 1 ) * lda );
      if ( aapp.getValue() > big ) {
        info.setValue( -6 );
        Blas2.xerbla( 'dgesvj', -info.getValue() );
        return;
      }
      aaqq.setValue( Math.sqrt( aaqq.getValue() ) );
      if ( aapp.getValue() < big / aaqq.getValue() && noscale ) {
        sva[ ioffsva + p - 1 ] = aapp.getValue() * aaqq.getValue();
      } else {
        noscale = false;
        sva[ ioffsva + p - 1 ] = aapp.getValue() * ( aaqq.getValue() * skl );
        if ( goscale ) {
          goscale = false;
          for ( var q = 1; q <=  p - 1; q ++ ) {
            sva[ ioffsva + q - 1 ] *= skl;
          } // 1873
        }
      }
    } // 1874
  } else if ( upper ) {
    for ( p = 1; p <= n; p ++ ) {
      aapp.setValue( 0. );
      aaqq.setValue( 1. );
      LaPack0.dlassq( p, A, 1, aapp, aaqq, ioffa + ( p - 1 ) * lda );
      if ( aapp.getValue() > big ) {
        info.setValue( -6 );
        Blas2.xerbla( 'dgesvj', -info.getValue() );
        return;
      }
      aaqq.setValue( Math.sqrt( aaqq.getValue() ) );
      if ( aapp.getValue() < big / aaqq.getValue() && noscale ) {
        sva[ ioffsva + p - 1 ] = aapp.getValue() * aaqq.getValue();
      } else {
        noscale = false;
        sva[ ioffsva + p - 1 ] = aapp.getValue() * ( aaqq.getValue() * skl );
        if ( goscale ) {
          goscale = false;
          for ( q = 1; q <=  p - 1; q ++ ) {
            sva[ ioffsva + q - 1 ] *= skl;
          } // 2873
        }
      }
    } // 2874
  } else {
    for ( p = 1; p <= n; p ++ ) {
      aapp.setValue( 0. );
      aaqq.setValue( 1. );
      LaPack0.dlassq( m, A, 1, aapp, aaqq, ioffa + ( p - 1 ) * lda );
      if ( aapp.getValue() > big ) {
        info.setValue( -6 );
        Blas2.xerbla( 'dgesvj', -info.getValue() );
        return;
      }
      aaqq.setValue( Math.sqrt( aaqq.getValue() ) );
      if ( aapp.getValue() < big / aaqq.getValue() && noscale ) {
        sva[ ioffsva + p - 1 ] = aapp.getValue() * aaqq.getValue();
      } else {
        noscale = false;
        sva[ ioffsva + p - 1 ] = aapp.getValue() * ( aaqq.getValue() * skl );
        if ( goscale ) {
          goscale = false;
          for ( q = 1; q <=  p - 1; q ++ ) {
            sva[ ioffsva + q - 1 ] *= skl;
          } // 3873
        }
      }
    } // 3874
  }
  if ( noscale ) skl = 1.;
  aapp.setValue( 0. );
  aaqq.setValue( big );
  for ( p = 1; p <= n; p ++ ) {
    if ( sva[ ioffsva + p - 1 ] != 0. ) {
      aaqq.setValue( Math.min( aaqq.getValue(), sva[ ioffsva + p - 1 ] ) );
    }
    aapp.setValue( Math.max( aapp.getValue(), sva[ ioffsva + p - 1 ] ) );
  } // 4781
  if ( aapp.getValue() == 0. ) {
    if ( lsvec ) LaPack0.dlaset( 'G', m, n, 0., 1., A, lda, ioffa );
    work[ ioffwork ] = 1.;
    work[ ioffwork + 1 ] = 0.;
    work[ ioffwork + 2 ] = 0.;
    work[ ioffwork + 3 ] = 0.;
    work[ ioffwork + 4 ] = 0.;
    work[ ioffwork + 5 ] = 0.;
    return;
  }
  var ierr = new IntReference();
  if ( n == 1 ) {
    if ( lsvec ) {
      LaPack1.dlascl( 'G', 0, 0, sva[ ioffsva ], skl, m, 1, A, lda,
        ierr, ioffa );
    }
    work[ ioffwork ] = 1. / skl;
    work[ ioffwork + 1 ] = ( sva[ ioffsva ] >= sfmin ? 1. : 0. );
    work[ ioffwork + 2 ] = 0.;
    work[ ioffwork + 3 ] = 0.;
    work[ ioffwork + 4 ] = 0.;
    work[ ioffwork + 5 ] = 0.;
    return;
  }
  var sn = Math.sqrt( sfmin / epsln );
  var temp1 = Math.sqrt( big / Number( n ) );
  if ( aapp.getValue() <= sn || aaqq.getValue() >= temp1 ||
  ( sn <= aaqq.getValue() && aapp.getValue() <= temp1 ) ) {
    temp1 = Math.min( big, temp1 / aapp.getValue() );
  } else if ( aaqq.getValue() <= sn && aapp.getValue() <= temp1 ) {
    temp1 = Math.min( sn / aaqq.getValue(),
      big / ( aapp.getValue() * Math.sqrt( Number( n ) ) ) );
  } else if ( aaqq.getValue() >= sn && aapp.getValue() >= temp1 ) {
    temp1 = Math.max( sn / aaqq.getValue() , temp1 / aapp.getValue() );
  } else if ( aaqq.getValue() <= sn && aapp.getValue() >= temp1 ) {
    temp1 = Math.min( sn / aaqq.getValue(),
      big / ( Math.sqrt( Number( n ) ) * aapp.getValue() ) );
  } else temp1 = 1.
  if ( temp1 != 1. ) {
    LaPack1.dlascl( 'G', 0, 0, 1., temp1, n, 1, sva, n, ierr,
      ioffsva );
  }
  skl *= temp1;
  if ( skl != 1. ) {
    LaPack1.dlascl( joba, 0, 0, 1., skl, m, n, A, lda, ierr, ioffa );
    skl = 1. / skl;
  }
  var emptsw = ( n * ( n - 1 ) ) / 2;
  var notrot = 0;
  fastr[ 0 ] = 0.;
  for ( q = 1; q <= n; q ++ ) work[ ioffwork + q - 1 ] = 1.;
  var swband = 3;
  var kbl = Math.min( 8, n );
  var nbl = n / kbl;
  if ( nbl * kbl != n ) nbl ++;
  var blskip = kbl * kbl;
  var rowskip = Math.min( 5, kbl );
  var lkahead = 1;
  if ( ( lower || upper ) && n > Math.max( 64, 4 * kbl) ) {
    var n4 = n / 4;
    var n2 = n / 2;
    var n34 = 3 * n4;
    q = ( applv ? 0 : 1 );
    if ( lower ) {
      LaPack2.dgsvj0( jobv, m - n34, n - n34, A, lda, work, sva, mvl,
        V, ldv, epsln, sfmin, tol, 2, work, lwork - n, ierr,
        ioffa + n34 + n34 * lda, ioffwork + n34, ioffsva + n34,
        ioffv + n34 * q + n34 * ldv, ioffwork + n );
      LaPack2.dgsvj0( jobv, m - n2, n34 - n2, A, lda, work, sva, mvl,
        V, ldv, epsln, sfmin, tol, 2, work, lwork - n, ierr,
        ioffa + n2 + n2 * lda, ioffwork + n2, ioffsva + n2,
        ioffv + n2 * q + n2 * ldv, ioffwork + n );
      LaPack2.dgsvj1( jobv, m - n2, n - n2, n4, A, lda, work, sva, mvl,
        V, ldv, epsln, sfmin, tol, 1, work, lwork - n, ierr,
        ioffa + n2 + n2 * lda, ioffwork + n2, ioffsva + n2,
        ioffv + n2 * q + n2 * ldv, ioffwork + n );
      LaPack2.dgsvj0( jobv, m - n4, n2 - n4, A, lda, work, sva, mvl,
        V, ldv, epsln, sfmin, tol, 1, work, lwork - n, ierr,
        ioffa + n4 + n4 * lda, ioffwork + n4, ioffsva + n4,
        ioffv + n4 * q + n4 * ldv, ioffwork + n );
      LaPack2.dgsvj0( jobv, m, n4, A, lda, work, sva, mvl, V, ldv,
        epsln, sfmin, tol, 1, work, lwork - n, ierr,
        ioffa, ioffwork, ioffsva, ioffv, ioffwork + n );
      LaPack2.dgsvj1( jobv, m, n2, n4, A, lda, work, sva, mvl, V, ldv,
        epsln, sfmin, tol, 1, work, lwork - n, ierr,
        ioffa, ioffwork, ioffsva, ioffv, ioffwork + n );
    } else if ( upper ) {
      LaPack2.dgsvj0( jobv, n4, n4, A, lda, work, sva, mvl, V, ldv,
        epsln, sfmin, tol, 2, work, lwork - n, ierr,
        ioffa, ioffwork, ioffsva, ioffv, ioffwork + n );
      LaPack2.dgsvj0( jobv, n2, n4, A, lda, work, sva, mvl, V, ldv,
        epsln, sfmin, tol, 1, work, lwork - n, ierr,
        ioffa + n4 * lda, ioffwork + n4, ioffsva + n4,
        ioffv + n4 * q + n4 * ldv, ioffwork + n );
      LaPack2.dgsvj1( jobv, n2, n2, n4, A, lda, work, sva, mvl, V, ldv,
        epsln, sfmin, tol, 1, work, lwork - n, ierr,
        ioffa, ioffwork, ioffsva, ioffv, ioffwork + n );
      LaPack2.dgsvj0( jobv, n2 + n4, n4, A, lda, work, sva, mvl,
        V, ldv, epsln, sfmin, tol, 1, work, lwork - n, ierr,
        ioffa + n2 * lda, ioffwork + n2, ioffsva + n2,
        ioffv + n2 * q + n2 * ldv, ioffwork + n );
    }
  }
  var goto1994 = false;
  for ( var i = 1; i <= nsweep; i ++ ) { // ==> 1993
    var mxaapq = 0.;
    var mxsinj = 0.;
    var iswrot = 0;
    notrot = 0;
    var pskipped = 0;
    for ( var ibr = 1; ibr <= nbl; ibr ++ ) { // ==> 2000
      var igl = ( ibr - 1 ) * kbl + 1;
      for ( var ir1 = 0; ir1 <= Math.min( lkahead, nbl - ibr );
      ir1 ++ ) { // ==> 1002
        igl += ir1 * kbl;
        for ( p = igl; p <= Math.min( igl + kbl - 1, n - 1 ); p ++ )
        { // ==> 2001
          q = Blas1.idamax( n - p + 1, sva, 1, ioffsva + p - 1 )
            + p - 1;
          if ( p != q ) {
            Blas1.dswap( m, A, 1, A, 1, ioffa + ( p - 1 ) * lda,
              ioffa + ( q - 1 ) * lda );
            if ( rsvec ) {
              Blas1.dswap( mvl, V, 1, V, 1, ioffv + ( p - 1 ) * ldv,
                ioffv + ( q - 1 ) * ldv );
            }
            temp1 = sva[ ioffsva + p - 1 ];
            sva[ ioffsva + p - 1 ] = sva[ ioffsva + q - 1 ];
            sva[ ioffsva + q - 1 ] = temp1;
            temp1 = work[ ioffwork + p - 1 ];
            work[ ioffwork + p - 1 ] = work[ ioffwork + q - 1 ];
            work[ ioffwork + q - 1 ] = temp1;
          }
          if ( ir1 == 0 ) {
            if ( sva[ ioffsva + p - 1 ] < rootbig &&
            sva[ ioffsva + p - 1 ] > rootsfmin ) {
              sva[ ioffsva + p - 1 ] =
                Blas1.dnrm2( m, A, 1, ioffa + ( p - 1 ) * lda )
                * work[ ioffwork + p - 1 ];
            } else {
              var temp1refReference = new NumberReference( 0. );
              aapp.setValue( 1. );
              LaPack0.dlassq( m, A, 1, temp1ref, aapp,
                ioffa + ( p - 1 ) * lda );
              temp1 = temp1ref.getValue();
              sva[ ioffsva + p - 1 ] = temp1 * Math.sqrt( aapp.getValue() )
                * work[ ioffwork + p - 1 ];
            }
            aapp.setValue( sva[ ioffsva + p - 1 ] );
          } else aapp.setValue( sva[ ioffsva + p - 1 ] );
          if ( aapp.getValue() > 0. ) {
            pskipped = 0;
            for ( q = p + 1; q <= Math.min( igl + kbl - 1, n ); q ++ )
            { // ==> 2002
              aaqq.setValue( sva[ ioffsva + q - 1 ] );
              if ( aaqq.getValue() > 0. ) {
                var aapp0 = aapp.getValue();
                if ( aaqq.getValue() >= 1. ) {
                  var rotok =
                    ( ( small * aapp.getValue() ) <= aaqq.getValue() );
                  if ( aapp.getValue() < big / aaqq.getValue() ) {
                    var aapq = ( Blas1.ddot( m, A, 1, A, 1,
                      ioffa + ( p - 1 ) * lda,
                      ioffa + ( q - 1 ) * lda )
                      * work[ ioffwork + p - 1 ]
                      * work[ ioffwork + q - 1 ] / aaqq.getValue() )
                      / aapp.getValue();
                  } else {
                    Blas1.dcopy( m, A, 1, work, 1,
                      ioffa + ( p - 1 ) * lda, ioffwork + n );
                    LaPack1.dlascl( 'G', 0, 0, aapp.getValue(),
                      work[ ioffwork + p - 1 ], m, 1, work, lda, ierr,
                      ioffwork + n );
                    aapq = Blas1.ddot( m, work, 1, A, 1, ioffwork + n,
                      ioffa + ( q - 1 ) * lda ) * work[ ioffwork + q ]
                      / aaqq.getValue();
                  }
                } else {
                  rotok = ( aapp.getValue() <= aaqq.getValue() / small );
                  if ( aapp.getValue() > small / aaqq.getValue() ) {
                    aapq = ( Blas1.ddot( m, A, 1, A, 1,
                      ioffa + ( p - 1 ) * lda,
                      ioffa + ( q - 1 ) * lda )
                      * work[ ioffwork + p - 1 ]
                      * work[ ioffwork + q - 1 ] / aaqq.getValue() )
                      / aapp.getValue();
                  } else {
                    Blas1.dcopy( m, A, 1, work, 1,
                      ioffa + ( q - 1 ) * lda, ioffwork + n );
                    LaPack1.dlascl( 'G', 0, 0, aaqq.getValue(),
                      work[ ioffwork + q - 1 ], m, 1, work, lda, ierr,
                      ioffwork + n );
                    aapq = Blas1.ddot( m, work, 1, A, 1, ioffwork + n,
                      ioffa + ( p - 1 ) * lda ) * work[ ioffwork + p ]
                      / aapp.getValue();
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
                    var aqoap = aaqq.getValue() / aapp.getValue();
                    var apoaq = aapp.getValue() / aaqq.getValue();
                    var theta =
                      - 0.5 * Math.abs( aqoap - apoaq ) / aapq;
                    if ( Math.abs( theta ) > bigtheta ) {
                      var t = 0.5 / theta;
                      fastr[ 2 ] = t * work[ ioffwork + p - 1 ]
                        / work[ ioffwork + q - 1 ];
                      fastr[ 3 ] = - t * work[ ioffwork + q - 1 ]
                        / work[ ioffwork + p - 1 ];
                      Blas1.drotm( m, A, 1, A, 1, fastr,
                        ioffa + ( p - 1 ) * lda,
                        ioffa + ( q - 1 ) * lda, 0 );
                      if ( rsvec ) {
                        Blas1.drotm( mvl, V, 1, V, 1, fastr,
                          ioffv + ( p - 1 ) * ldv,
                          ioffv + ( q - 1 ) * ldv, 0 );
                      }
                      sva[ ioffsva + q - 1 ] = aaqq.getValue() * Math.sqrt(
                        Math.max( 0., 1. + t * apoaq * aapq ) );
                      aapp.getValue() = aapp.getValue() * Math.sqrt( Math.max(
                        0., 1. - t * aqoap * aapq ) );
                      mxsinj = Math.max( mxsinj, Math.abs( t ) );
                    } else {
                      var thsign = ( aapq >= 0. ? -1. : 1. );
                      t = 1. / ( theta
                        + thsign * Math.sqrt( 1. + theta * theta ) );
                      var cs = Math.sqrt( 1. / ( 1. + t * t ) );
                      sn = t * cs;
                      mxsinj = Math.max( mxsinj, Math.abs( sn ) );
                      sva[ ioffsva + q - 1 ] = aaqq.getValue() * Math.sqrt(
                        Math.max( 0., 1. + t * apoaq * aapq ) );
                      aapp.getValue() = aapp.getValue() * Math.sqrt( Math.max(
                        0., 1. - t * aqoap * aapq ) );
                      apoaq = work[ ioffwork + p - 1 ]
                        / work[ ioffwork + q - 1 ];
                      aqoap = work[ ioffwork + q - 1 ]
                        / work[ ioffwork + p - 1 ];
                      if ( work[ ioffwork + p - 1 ] >= 1. ) {
                        if ( work[ ioffwork + q - 1 ] >= 1. ) {
                          fastr[ 2 ] = t * apoaq;
                          fastr[ 3 ] = - t * aqoap;
                          work[ ioffwork + p - 1 ] *= cs;
                          work[ ioffwork + q - 1 ] *= cs;
                          Blas1.drotm( m, A, 1, A, 1, fastr,
                            ioffa + ( p - 1 ) * lda,
                            ioffa + ( q - 1 ) * lda, 0 );
                          if ( rsvec ) {
                            Blas1.drotm( mvl, V, 1, V, 1, fastr,
                              ioffv + ( p - 1 ) * ldv,
                              ioffv + ( q - 1 ) * ldv, 0 );
                          }
                        } else {
                          Blas1.daxpy( m, -t * aqoap, A, 1, A, 1,
                            ioffa + ( q - 1 ) * lda,
                            ioffa + ( p - 1 ) * lda );
                          Blas1.daxpy( m, cs * sn * apoaq, A, 1, A, 1,
                            ioffa + ( p - 1 ) * lda,
                            ioffa + ( q - 1 ) * lda );
                          work[ ioffwork + p - 1 ] *= cs;
                          work[ ioffwork + q - 1 ] /= cs;
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
                        if ( work[ ioffwork + q - 1 ] >= 1. ) {
                          Blas1.daxpy( m, t * apoaq, A, 1, A, 1,
                            ioffa + ( p - 1 ) * lda,
                            ioffa + ( q - 1 ) * lda );
                          Blas1.daxpy( m, -cs * sn * aqoap, A, 1, A, 1,
                            ioffa + ( q - 1 ) * lda,
                            ioffa + ( p - 1 ) * lda );
                          work[ ioffwork + p - 1 ] /= cs;
                          work[ ioffwork + q - 1 ] *= cs;
                          if ( rsvec ) {
                            Blas1.daxpy( mvl, t * apoaq, V, 1, V, 1,
                              ioffv + ( p - 1 ) * ldv,
                              ioffv + ( q - 1 ) * ldv );
                            Blas1.daxpy( mvl, -cs * sn * aqoap, V, 1,
                              V, 1, ioffv + ( q - 1 ) * ldv,
                              ioffv + ( p - 1 ) * ldv );
                          }
                        } else {
                          if ( work[ ioffwork + p - 1 ] >=
                          work[ ioffwork + q - 1 ] ) {
                            Blas1.daxpy( m, -t * aqoap, A, 1, A, 1,
                              ioffa + ( q - 1 ) * lda,
                              ioffa + ( p - 1 ) * lda );
                            Blas1.daxpy( m, cs * sn * apoaq, A, 1,
                              A, 1, ioffa + ( p - 1 ) * lda,
                              ioffa + ( q - 1 ) * lda );
                            work[ ioffwork + p - 1 ] *= cs;
                            work[ ioffwork + q - 1 ] /= cs;
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
                            work[ ioffwork + p - 1 ] /= cs;
                            work[ ioffwork + q - 1 ] *= cs;
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
                    Blas1.dcopy( m, A, 1, work, 1,
                      ioffa + ( p - 1 ) * lda, ioffwork + n );
                    LaPack1.dlascl( 'G', 0, 0, aapp.getValue(), 1., m, 1,
                      work, lda, ierr, ioffwork + n );
                    LaPack1.dlascl( 'G', 0, 0, aaqq.getValue(), 1., m, 1,
                      A, lda, ierr, ioffa + ( q - 1 ) * lda );
                    temp1 = - aapq * work[ ioffwork + p - 1 ]
                      / work[ ioffwork + q - 1 ];
                    Blas1.daxpy( m, temp1, work, 1, A, 1, ioffwork + n,
                      ioffa + ( q - 1 ) * lda );
                    LaPack1.dlascl( 'G', 0, 0, 1., aaqq.getValue(), m, 1,
                      A, lda, ierr, ioffa + ( q - 1 ) * lda );
                    sva[ ioffsva + q - 1 ] = aaqq.getValue() * Math.sqrt(
                      Math.max( 0., 1. - aapq * aapq ) );
                    mxsinj = Math.max( mxsinj, sfmin );
                  }
                  if ( Math.pow( sva[ ioffsva + q - 1 ] / aaqq.getValue(),
                  2 ) <= rooteps ) {
                    if ( aaqq.getValue() < rootbig &&
                    aaqq.getValue() > rootsfmin ) {
                      sva[ ioffsva + q - 1 ] =
                        Blas1.dnrm2( m, A, 1, ioffa + ( q - 1 ) * lda )
                        * work[ ioffwork + q - 1 ];
                    } else {
                      var trefReference =
                        new NumberReference( 0. );
                      aaqq.setValue( 1. );
                      LaPack0.dlassq( m, A, 1, tref, aaqq,
                        ioffa + ( q - 1 ) * lda );
                      t = tref.getValue();
                      sva[ ioffsva + q - 1 ] = t * Math.sqrt(
                        aaqq.getValue() ) * work[ ioffwork + q - 1 ];
                    }
                  }
                  if ( aapp.getValue() / aapp0 <= rooteps ) {
                    if ( aapp.getValue() < rootbig &&
                    aapp.getValue() > rootsfmin ) {
                      aapp.getValue() =
                        Blas1.dnrm2( m, A, 1, ioffa + ( p - 1 ) * lda )
                        * work[ ioffwork + p - 1 ];
                    } else {
                      tref.setValue( 0. );
                      aapp.setValue( 1. );
                      LaPack0.dlassq( m, A, 1, tref, aapp,
                        ioffa + ( p - 1 ) * lda );
                      t = tref.getValue();
                      aapp.getValue() = t * Math.sqrt( aapp.getValue() )
                        * work[ ioffwork + p - 1 ];
                    }
                    sva[ ioffsva + p - 1 ] = aapp.getValue();
                  }
                } else {
                  if ( ir1 == 0 ) notrot ++;
                  pskipped ++;
                }
              } else {
                if ( ir1 == 0 ) notrot ++;
                pskipped ++;
              }
              if ( i <= swband && pskipped > rowskip ) {
                if ( ir1 == 0 ) aapp.setValue( - aapp.getValue() );
                notrot = 0;
                break; // goto 2103
              }
            } // 2002
            sva[ ioffsva + p - 1 ] = aapp.getValue(); // 2103
          } else {
            sva[ ioffsva + p - 1 ] = aapp.getValue();
            if ( ir1 == 0 && aapp.getValue() == 0. ) {
              notrot += Math.min( igl + kbl - 1, n ) - p;
            }
          }
        } // 2001
      } // 1002
      igl = ( ibr - 1 ) * kbl + 1;
      var goto2011 = false;
      for ( var jbc = ibr + 1; jbc <= nbl; jbc ++ ) { // ==> 2010
        var jgl = ( jbc - 1 ) * kbl + 1;
        var ijblsk = 0;
        for ( p = igl; p <= Math.min( igl + kbl - 1, n ); p ++ )
        { // ==> 2100
          aapp.setValue( sva[ ioffsva + p - 1 ] );
          if ( aapp.getValue() > 0. ) {
            pskipped = 0;
            for ( q = jgl; q <= Math.min( jgl + kbl - 1, n ); q ++ )
            { // ==> 2200
              aaqq.setValue( sva[ ioffsva + q - 1 ] );
              if ( aaqq.getValue() > 0. ) {
                aapp0 = aapp.getValue();
                if ( aaqq.getValue() >= 1. ) {
                  if ( aapp.getValue() >= aaqq.getValue() ) {
                    rotok = ( small * aapp.getValue() <= aaqq.getValue() );
                  } else {
                    rotok = ( small * aaqq.getValue() <= aapp.getValue() );
                  }
                  if ( aapp.getValue() < big / aaqq.getValue() ) {
                    aapq = ( Blas1.ddot( m, A, 1, A, 1,
                      ioffa + ( p - 1 ) * lda,
                      ioffa + ( q - 1 ) * lda )
                      * work[ ioffwork + p - 1 ]
                      * work[ ioffwork + q - 1 ] / aaqq.getValue() )
                      / aapp.getValue();
                  } else {
                    Blas1.dcopy( m, A, 1, work, 1,
                      ioffa + ( p - 1 ) * lda, ioffwork + n );
                    LaPack1.dlascl( 'G', 0, 0, aapp.getValue(),
                      work[ ioffwork + p - 1 ], m, 1, work, lda, ierr,
                      ioffwork + n );
                    aapq = Blas1.ddot( m, work, 1, A, 1, ioffwork + n,
                      ioffa + ( q - 1 ) * lda ) * work[ ioffwork + q ]
                      / aaqq.getValue();
                  }
                } else {
                  if ( aapp.getValue() >= aaqq.getValue() ) {
                    rotok = ( aapp.getValue() <= aaqq.getValue() / small );
                  } else {
                    rotok = ( aaqq.getValue() <= aapp.getValue() / small );
                  }
                  if ( aapp.getValue() > small / aaqq.getValue() ) {
                    aapq = ( Blas1.ddot( m, A, 1, A, 1,
                      ioffa + ( p - 1 ) * lda,
                      ioffa + ( q - 1 ) * lda )
                      * work[ ioffwork + p - 1 ]
                      * work[ ioffwork + q - 1 ] / aaqq.getValue() )
                      / aapp.getValue();
                  } else {
                    Blas1.dcopy( m, A, 1, work, 1,
                      ioffa + ( q - 1 ) * lda, ioffwork + n );
                    LaPack1.dlascl( 'G', 0, 0, aaqq.getValue(),
                      work[ ioffwork + q - 1 ], m, 1, work, lda, ierr,
                      ioffwork + n );
                    aapq = Blas1.ddot( m, work, 1, A, 1, ioffwork + n,
                      ioffa + ( p - 1 ) * lda ) * work[ ioffwork + p ]
                      / aapp.getValue();
                  }
                }
                mxaapq = Math.max( mxaapq, Math.abs( aapq ) );
                if ( Math.abs( aapq ) > tol ) {
                  notrot = 0;
                  pskipped = 0;
                  iswrot ++;
                  if ( rotok ) {
                    aqoap = aaqq.getValue() / aapp.getValue();
                    apoaq = aapp.getValue() / aaqq.getValue();
                    theta = - 0.5 * Math.abs( aqoap - apoaq ) / aapq;
                    if ( aaqq.getValue() > aapp0 ) theta = - theta;
                    if ( Math.abs( theta ) > bigtheta ) {
                      t = 0.5 / theta;
                      fastr[ 2 ] = t * work[ ioffwork + p - 1 ]
                        / work[ ioffwork + q - 1 ];
                      fastr[ 3 ] = - t * work[ ioffwork + q - 1 ]
                        / work[ ioffwork + p - 1 ];
                      Blas1.drotm( m, A, 1, A, 1, fastr,
                        ioffa + ( p - 1 ) * lda,
                        ioffa + ( q - 1 ) * lda, 0 );
                      if ( rsvec ) {
                        Blas1.drotm( mvl, V, 1, V, 1, fastr,
                          ioffv + ( p - 1 ) * ldv,
                          ioffv + ( q - 1 ) * ldv, 0 );
                      }
                      sva[ ioffsva + q - 1 ] = aaqq.getValue() * Math.sqrt(
                        Math.max( 0., 1. + t * apoaq * aapq ) );
                      aapp.getValue() = aapp.getValue() * Math.sqrt( Math.max(
                        0., 1. - t * aqoap * aapq ) );
                      mxsinj = Math.max( mxsinj, Math.abs( t ) );
                    } else {
                      thsign = ( aapq >= 0. ? -1. : 1. );
                      if ( aaqq.getValue() > aapp0 ) thsign = - thsign;
                      t = 1. / ( theta
                        + thsign * Math.sqrt( 1. + theta * theta ) );
                      cs = Math.sqrt( 1. / ( 1. + t * t ) );
                      sn = t * cs;
                      mxsinj = Math.max( mxsinj, Math.abs( sn ) );
                      sva[ ioffsva + q - 1 ] = aaqq.getValue() * Math.sqrt(
                        Math.max( 0., 1. + t * apoaq * aapq ) );
                      aapp.setValue( aapp.getValue() * Math.sqrt( Math.max(
                        0., 1. - t * aqoap * aapq ) ) );
                      apoaq = work[ ioffwork + p - 1 ]
                        / work[ ioffwork + q - 1 ];
                      aqoap = work[ ioffwork + q - 1 ]
                        / work[ ioffwork + p - 1 ];
                      if ( work[ ioffwork + p - 1 ] >= 1. ) {
                        if ( work[ ioffwork + q - 1 ] >= 1. ) {
                          fastr[ 2 ] = t * apoaq;
                          fastr[ 3 ] = - t * aqoap;
                          work[ ioffwork + p - 1 ] *= cs;
                          work[ ioffwork + q - 1 ] *= cs;
                          Blas1.drotm( m, A, 1, A, 1, fastr,
                            ioffa + ( p - 1 ) * lda,
                            ioffa + ( q - 1 ) * lda, 0 );
                          if ( rsvec ) {
                            Blas1.drotm( mvl, V, 1, V, 1, fastr,
                              ioffv + ( p - 1 ) * ldv,
                              ioffv + ( q - 1 ) * ldv, 0 );
                          }
                        } else {
                          Blas1.daxpy( m, -t * aqoap, A, 1, A, 1,
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
                          work[ ioffwork + p - 1 ] *= cs;
                          work[ ioffwork + q - 1 ] /= cs;
                        }
                      } else {
                        if ( work[ ioffwork + q - 1 ] >= 1. ) {
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
                          work[ ioffwork + p - 1 ] /= cs;
                          work[ ioffwork + q - 1 ] *= cs;
                        } else {
                          if ( work[ ioffwork + p - 1 ] >=
                          work[ ioffwork + q - 1 ] ) {
                            Blas1.daxpy( m, -t * aqoap, A, 1, A, 1,
                              ioffa + ( q - 1 ) * lda,
                              ioffa + ( p - 1 ) * lda );
                            Blas1.daxpy( m, cs * sn * apoaq, A, 1,
                              A, 1, ioffa + ( p - 1 ) * lda,
                              ioffa + ( q - 1 ) * lda );
                            work[ ioffwork + p - 1 ] *= cs;
                            work[ ioffwork + q - 1 ] /= cs;
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
                            work[ ioffwork + p - 1 ] /= cs;
                            work[ ioffwork + q - 1 ] *= cs;
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
                    if ( aapp.getValue() > aaqq.getValue() ) {
                      Blas1.dcopy( m, A, 1, work, 1,
                        ioffa + ( p - 1 ) * lda, ioffwork + n );
                      LaPack1.dlascl( 'G', 0, 0, aapp.getValue(), 1., m, 1,
                        work, lda, ierr, ioffwork + n );
                      LaPack1.dlascl( 'G', 0, 0, aaqq.getValue(), 1., m, 1,
                        A, lda, ierr, ioffa + ( q - 1 ) * lda );
                      temp1 = - aapq * work[ ioffwork + p - 1 ]
                        / work[ ioffwork + q - 1 ];
                      Blas1.daxpy( m, temp1, work, 1, A, 1,
                        ioffwork + n, ioffa + ( q - 1 ) * lda );
                      LaPack1.dlascl( 'G', 0, 0, 1., aaqq.getValue(), m, 1,
                        A, lda, ierr, ioffa + ( q - 1 ) * lda );
                      sva[ ioffsva + q - 1 ] = aaqq.getValue() * Math.sqrt(
                        Math.max( 0., 1. - aapq * aapq ) );
                      mxsinj = Math.max( mxsinj, sfmin )
                    } else {
                      Blas1.dcopy( m, A, 1, work, 1,
                        ioffa + ( q - 1 ) * lda, ioffwork + n );
                      LaPack1.dlascl( 'G', 0, 0, aaqq.getValue(), 1., m, 1,
                        work, lda, ierr, ioffwork + n );
                      LaPack1.dlascl( 'G', 0, 0, aapp.getValue(), 1., m, 1,
                        A, lda, ierr, ioffa + ( p - 1 ) * lda );
                      temp1 = - aapq * work[ ioffwork + q - 1 ]
                        / work[ ioffwork + p - 1 ];
                      Blas1.daxpy( m, temp1, work, 1, A, 1,
                        ioffwork + n, ioffa + ( p - 1 ) * lda );
                      LaPack1.dlascl( 'G', 0, 0, 1., aapp.getValue(), m, 1,
                        A, lda, ierr, ioffa + ( p - 1 ) * lda );
                      sva[ ioffsva + p - 1 ] = aapp.getValue() * Math.sqrt(
                        Math.max( 0., 1. - aapq * aapq ) );
                      mxsinj = Math.max( mxsinj, sfmin )
                    }
                  }
                  if ( Math.pow( sva[ ioffsva + q - 1 ] / aaqq.getValue(),
                  2 ) <= rooteps ) {
                    if ( aaqq.getValue() < rootbig &&
                    aaqq.getValue() > rootsfmin ) {
                      sva[ ioffsva + q - 1 ] =
                        Blas1.dnrm2( m, A, 1, ioffa + ( q - 1 ) * lda )
                        * work[ ioffwork + q - 1 ];
                  } else {
                    tref.setValue( 0. );
                    aaqq.setValue( 1. );
                    LaPack0.dlassq( m, A, 1, tref, aaqq,
                      ioffa + ( q - 1 ) * lda );
                    t = tref.getValue();
                    sva[ ioffsva + q - 1 ] = t * Math.sqrt(
                      aaqq.getValue() ) * work[ ioffwork + q - 1 ];
                    }
                  }
                  if ( Math.pow( aapp.getValue() / aapp0, 2 ) <= rooteps ) {
                    if ( aapp.getValue() < rootbig &&
                    aapp.getValue() > rootsfmin ) {
                      aapp.setValue(
                        Blas1.dnrm2( m, A, 1, ioffa + ( p - 1 ) * lda )
                        * work[ ioffwork + p - 1 ] );
                    } else {
                      tref.setValue( 0. );
                      aapp.setValue( 1. );
                      LaPack0.dlassq( m, A, 1, tref, aapp,
                        ioffa + ( p - 1 ) * lda );
                      t = tref.getValue();
                      aapp.setValue( t * Math.sqrt( aapp.getValue() )
                        * work[ ioffwork + p - 1 ] );
                    }
                    sva[ ioffsva + p - 1 ] = aapp.getValue();
                  }
                } else {
                  notrot ++;
                  pskipped ++;
                  ijblsk ++;
                }
              } else {
                notrot ++;
                pskipped ++;
                ijblsk ++;
              }
              if ( i <= swband && ijblsk >= blskip ) {
                sva[ ioffsva + p - 1 ] = aapp.getValue();
                notrot = 0;
                goto2011 = true;
                break;
              }
              if ( i <= swband && pskipped > rowskip ) {
                aapp.setValue( - aapp.getValue() );
                notrot = 0;
                break; // goto 2203
              }
            } // 2200
            sva[ ioffsva + p - 1 ] = aapp.getValue(); // 2203
          } else {
            if ( aapp.getValue() == 0. ) {
              notrot += Math.min( jgl + kbl - 1, n ) - jgl + 1;
            }
            if ( aapp.getValue() < 0. ) notrot = 0;
          }
        } // 2100
        if ( goto2011 ) break;
      } // 2010
      for ( p = igl; p <= Math.min( igl + kbl - 1, n ); p ++ ) { //2011
        sva[ ioffsva + p - 1 ] = Math.abs( sva[ ioffsva + p - 1 ] );
      }
    } // 2000
    if ( sva[ ioffsva + n - 1 ] < rootbig &&
    sva[ ioffsva + n - 1 ] > rootsfmin ) {
      sva[ ioffsva + n - 1 ] =
        Blas1.dnrm2( m, A, 1, ioffa + ( n - 1 ) * lda )
        * work[ ioffwork + n - 1 ];
    } else {
      tref.setValue( 0. );
      aapp.setValue( 1. );
      LaPack0.dlassq( m, A, 1, tref, aapp, ioffa + ( n - 1 ) * lda );
      t = tref.getValue();
      sva[ ioffsva + n - 1 ] = t * Math.sqrt( aapp.getValue() )
        * work[ ioffwork + n - 1 ];
    }
    if ( ( i < swband && mxaapq <= roottol ) || iswrot <= n ) {
      swband = i;
    }
    if ( i > swband + 1 &&
    mxaapq < Math.sqrt( Number( n ) ) * tol &&
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
  n2 = 0; // 1995
  n4 = 0;
  for ( p = 1; p <= n - 1; p ++ ) {
    q = Blas1.idamax( n - p + 1, sva, 1, ioffsva + p - 1 ) + p - 1;
    if ( p != q ) {
      temp1 = sva[ ioffsva + p - 1 ];
      sva[ ioffsva + p - 1 ] = sva[ ioffsva + q - 1 ];
      sva[ ioffsva + q - 1 ] = temp1;
      temp1 = work[ ioffwork + p - 1 ];
      work[ ioffwork + p - 1 ] = work[ ioffwork + q - 1 ];
      work[ ioffwork + q - 1 ] = temp1;
      Blas1.dswap( m, A, 1, A, 1, ioffa + ( p - 1 ) * lda,
        ioffa + ( q - 1 ) * lda );
      if ( rsvec ) {
        Blas1.dswap( mvl, V, 1, V, 1, ioffv + ( p - 1 ) * ldv,
          ioffv + ( q - 1 ) * ldv );
      }
    }
    if ( sva[ ioffsva + p - 1 ] != 0. ) {
      n4 ++;
      if ( sva[ ioffsva + p - 1 ] * skl > sfmin ) n2 ++;
    }
  } // 5991
  if ( sva[ ioffsva + n - 1 ] != 0. ) {
    n4 ++;
    if ( sva[ ioffsva + n - 1 ] * skl > sfmin ) n2 ++;
  }
  if ( lsvec || uctol ) {
    for ( p = 1; p <= n2; p ++ ) {
      Blas1.dscal( m,
        work[ ioffwork + p - 1 ] / sva[ ioffsva + p - 1 ], A, 1,
        ioffa + ( p - 1 ) * lda );
    } // 1998
  }
  if ( rsvec ) {
    if ( applv ) {
      for ( p = 1; p <= n; p ++ ) {
        Blas1.dscal( mvl, work[ ioffwork + p - 1 ], V, 1,
          ioffv + ( p - 1 ) * ldv );
      }
    } else {
      for ( p = 1; p <= n; p ++ ) {
        temp1 = 1. / Blas1.dnrm2( mvl, V, 1, ioffv + ( p - 1 ) * ldv );
        Blas1.dscal( mvl, temp1, V, 1, ioffv + ( p - 1 ) * ldv );
      } // 2399
    }
  }
  if ( ( skl > 1. && sva[ ioffsva ] < big / skl ) ||
  ( skl < 1. && sva[ ioffsva + n2 - 1 ] > sfmin / skl ) ) {
    for ( p = 1; p <= n; p ++ ) {
      sva[ ioffsva + p - 1 ] *= skl;
    }
    skl = 1.;
  }
  work[ ioffwork ] = skl;
  work[ ioffwork + 1 ] = Number( n4 );
  work[ ioffwork + 2 ] = Number( n2 );
  work[ ioffwork + 3 ] = Number( i );
  work[ ioffwork + 4 ] = mxaapq;
  work[ ioffwork + 5 ] = mxsinj;
}
//*************************************************************************
LaPack3.dlaed3 = function( k, n, n1, d, Q, ldq, rho, dlamda, q2,
indx, ctot, w, s, info, ioffd, ioffq, ioffdlamda, ioffq2, ioffindx,
ioffctot, ioffw, ioffs ) {
  throw new Error("not tested");
  if ( k < 0 ) info.setValue( -1 );
  else if ( n < k ) info.setValue( -2 );
  else if ( ldq < Math.max( 1, n ) ) info.setValue( -6 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlaed3', - info.getValue() );
    return;
  }
  if ( k == 0 ) return;
  for( var i = 1; i <= k; i ++ ) {
    dlamda[ ioffdlamda + i - 1 ] =
      LaPack0.dlamc3( dlamda[ ioffdlamda + i - 1 ],
      dlamda[ ioffdlamda + i - 1 ] ) - dlamda[ ioffdlamda + i - 1 ];
  }
  for ( var j = 1; j <= k; j ++ ) {
    var dlam = new NumberReference( d[ ioffd + j - 1 ] );
    LaPack2.dlaed4( k, j, dlamda, w, Q, rho, dlam, info,
      ioffdlamda, ioffw, ioffq + ( j - 1 ) * ldq );
    d[ ioffd + j - 1 ] = dlam.getValue();
    if ( info.getValue() != 0 ) return;
  }
  if ( k == 2 ) {
    for ( j = 1; j <= k; j ++ ) {
      w[ ioffw ] = Q[ ioffq + ( j - 1 ) * ldq ];
      w[ ioffw + 1 ] = Q[ ioffq + 1 + ( j - 1 ) * ldq ];
      var ii = indx[ ioffindx ];
      Q[ ioffq + ( j - 1 ) * ldq ] = w[ ioffw + ii - 1 ];
      ii = indx[ ioffindx + 1 ];
      Q[ ioffq + 1 + ( j - 1 ) * ldq ] = w[ ioffw + ii - 1 ];
    }
  }
  if ( k != 1 && k != 2 ) {
    Blas1.dcopy( k, w, 1, s, 1, ioffw, ioffs );
    Blas1.dcopy( k, Q, ldq + 1, w, 1, ioffq, ioffw );
    for ( j = 1; j <= k; j ++ ) {
      for ( i = 1; i <= j - 1; i ++ ) {
        w[ ioffw + i - 1 ] *= Q[ ioffq + i - 1 + ( j - 1 ) * ldq ]
          / ( dlamda[ ioffdlamda + i - 1 ]
          - dlamda[ ioffdlamda + j - 1 ] );
      }
      for ( i = j + 1; i <= k; i ++ ) {
        w[ ioffw + i - 1 ] *= Q[ ioffq + i - 1 + ( j - 1 ) * ldq ]
          / ( dlamda[ ioffdlamda + i - 1 ]
          - dlamda[ ioffdlamda + j - 1 ] );
      }
    }
    for ( i = 1; i <= k; i ++ ) {
      w[ ioffw + i - 1 ] = ( s[ ioffs + i - 1 ] >= 0. ?
        Math.sqrt( - w[ ioffw + i - 1 ] ) :
        - Math.sqrt( - w[ ioffw + i - 1 ] ) );
    }
    for ( j = 1; j <= k; j ++ ) {
      for ( i = 1; i <= k; i ++ ) {
        s[ ioffs + i - 1 ] = w[ ioffw + i - 1 ]
          / Q[ ioffq + i - 1 + ( j - 1 ) * ldq ];
      }
      var temp = Blas1.dnrm2( k, s, 1, ioffs );
      for ( i = 1; i <= k; i ++ ) {
        ii = indx[ ioffindx + i - 1 ];
        Q[ ioffq + i - 1 + ( j - 1 ) * ldq ] =
          s[ ioffs + ii - 1 ] / temp;
      }
    }
  }
  var n2 = n - n1;
  var n12 = ctot[ ioffctot ] + ctot[ ioffctot + 1 ];
  var n23 = ctot[ ioffctot + 1 ] + ctot[ ioffctot + 2 ];
  LaPack0.dlacpy( 'A', n23, k, Q, ldq, s, n23,
    ioffq + ctot[ ioffctot ], ioffs );
  var iq2 = n1 * n12 + 1;
  if ( n23 != 0 ) {
    Blas3.dgemm( 'N', 'N', n2, k, n23, 1., q2, n2, s, n23, 0., Q, ldq,
      ioffq2 + iq2, ioffs, ioffq + n1 );
  } else {
    LaPack0.dlaset( 'A', n2, k, 0., 0., Q, ldq, ioffq + n1 );
  }
  LaPack0.dlacpy( 'A', n12, k, Q, ldq, s, n12, ioffq, ioffs );
  if ( n12 != 0 ) {
    Blas3.dgemm( 'N', 'N', n1, k, n12, 1., q2, n1, s, n12, 0., Q, ldq,
      ioffq2, ioffs, ioffq );
  } else {
    LaPack0.dlaset( 'A', n1, k, 0., 0., Q, ldq, ioffq );
  }
}
/* old version:
LaPack3.dlaed3 = function( k, kstart, kstop, n, d, Q, ldq, rho,
cutpnt, dlamda, Q2, ldq2, indxc, ctot, w, S, lds, info, ioffd, ioffq,
ioffdlamda, ioffq2, ioffindxc, ioffctot, ioffw, ioffs ) {
  info.setValue( 0 );
  if ( k < 0 ) info.setValue( -1 );
  else if ( kstart < 1 || kstart > Math.max( 1, k ) ) info.setValue( -2 );
  else if ( Math.max( 1, kstop ) < kstart || kstop > Math.max( 1, k ) )
  {
    info.setValue( -3 );
  } else if ( n < k ) info.setValue( -4 );
  else if ( ldq < Math.max( 1, n ) ) info.setValue( -7 );
  else if ( ldq2 < Math.max( 1, n ) ) info.setValue( -12 );
  else if ( lds < Math.max( 1, k ) ) info.setValue( -17 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlaed3', - info.getValue() );
    return;
  }
  if ( k == 0 ) return;
  for( var i = 1; i <= n; i ++ ) {
    dlamda[ ioffdlamda + i - 1 ] =
      LaPack0.dlamc3( dlamda[ ioffdlamda + i - 1 ],
      dlamda[ ioffdlamda + i - 1 ] ) - dlamda[ ioffdlamda + i - 1 ];
  }
  var ktemp = kstop - kstart + 1;
  for ( var j = kstart; j <= kstop; j ++ ) {
    var dlamReference =
      new NumberReference( d[ ioffd + j - 1 ] );
    LaPack2.dlaed4( k, j, dlamda, w, Q, rho, dlam, info,
      ioffdlamda, ioffw, ioffq + ( j - 1 ) * ldq );
    if ( info.getValue() != 0 ) return;
  }
  if ( k == 1 || k == 2 ) {
    for ( i = 1; i <= k; i ++ ) {
      for ( j = 1; j <= k; j ++ ) {
        var jc = indxc[ ioffindxc + j - 1 ];
        S[ ioffs + j - 1 + ( i - 1 ) * lds ] =
          Q[ ioffq + jc - 1 + ( i - 1 ) * ldq ];
      }
    }
  } else {
    Blas1.dcopy( k, w, 1, S, 1, ioffw, ioffs );
    Blas1.dcopy( k, Q, ldq + 1, w, 1, ioffq, ioffw );
    for ( j = 1; j <= k; j ++ ) {
      for ( i = 1; i <= j - 1; i ++ ) {
        w[ ioffw + i - 1 ] *= Q[ ioffq + i - 1 + ( j - 1 ) * ldq ]
          / ( dlamda[ ioffdlamda + i - 1 ]
          - dlamda[ ioffdlamda + j - 1 ] );
      }
      for ( i = j + 1; i <= k; i ++ ) {
        w[ ioffw + i - 1 ] *= Q[ ioffq + i - 1 + ( j - 1 ) * ldq ]
          / ( dlamda[ ioffdlamda + i - 1 ]
          - dlamda[ ioffdlamda + j - 1 ] );
      }
    }
    for ( i = 1; i <= k; i ++ ) {
      w[ ioffw + i - 1 ] = ( S[ ioffs + i - 1 ] >= 0. ?
        Math.sqrt( - w[ ioffw + i - 1 ] ) :
        - Math.sqrt( - w[ ioffw + i - 1 ] ) );
    }
    for ( j = 1; j <= k; j ++ ) {
      for ( i = 1; i <= k; i ++ ) {
        Q[ ioffq + i - 1 + ( j - 1 ) * ldq ] = w[ ioffw + i - 1 ]
          / Q[ ioffq + i - 1 + ( j - 1 ) * ldq ];
      }
      var temp =
        Blas1.dnrm2( k, Q, 1, ioffq + ( j - 1 ) * ldq );
      for ( i = 1; i <= k; i ++ ) {
        jc = indxc[ ioffindxc + i - 1 ];
        S[ ioffs + i - 1 + ( j - 1 ) * lds ] =
          Q[ ioffq + jc - 1 + ( j - 1 ) * ldq ] / temp;
      }
    }
  } 120
  var parts = 0;
  if ( ctot[ ioffctot ] > 0 ) {
    parts ++;
    Blas3.dgemm( 'N', 'N', cutpnt, ktemp, ctot[ ioffctot ], 1., Q2,
      ldq2, S, lds, 0., Q, ldq, ioffq2, ioffs + ( kstart - 1 ) * lds,
      ioffq + ( kstart - 1 ) * ldq );
  }
  if ( ctot[ ioffctot + 1 ] > 0 ) {
    parts += 2;
    Blas3.dgemm( 'N', 'N', n - cutpnt, ktemp, ctot[ ioffctot + 1 ],
      1., Q2, ldq2, S, lds, 0., Q, ldq,
      cutpnt + ctot[ ioffctot ] * ldq2, 
      ioffs + ctot[ ioffctot ] + ( kstart - 1 ) * lds,
      ioffq + cutpnt + ( kstart - 1 ) * ldq );
  }
  if ( parts == 1 ) {
    LaPack0.dlaset( 'A', cutpnt, ktemp, 0., 0., Q, ldq, 
      ioffq + ( kstart - 1 ) * ldq );
  }
  if ( ctot[ ioffctot + 2 ] > 0 ) {
    if ( parts > 0 ) {
      Blas3.dgemm( 'N', 'N', n, ktemp, ctot[ ioffctot + 2 ], 1., Q2,
        ldq2, S, lds, 1., Q, ldq,
        ioffq2 + ( ctot[ ioffctot ] + ctot[ ioffctot + 1 ] ) * ldq2,
        ioffs + ctot[ ioffctot ] + ctot[ ioffctot + 1 ]
        + ( kstart - 1 ) * lds, ioffq + ( kstart - 1 ) * ldq );
    } else {
      Blas3.dgemm( 'N', 'N', n, ktemp, ctot[ ioffctot + 2 ], 1., Q2,
        ldq2, S, lds, 0., Q, ldq,
        ioffq2 + ( ctot[ ioffctot ] + ctot[ ioffctot + 1 ] ) * ldq2,
        ioffs + ctot[ ioffctot ] + ctot[ ioffctot + 1 ]
        + ( kstart - 1 ) * lds, ioffq + ( kstart - 1 ) * ldq );
    }
  }
}
*/
//*************************************************************************
LaPack3.dlaed9 = function( k, kstart, kstop, n, d, Q, ldq, rho,
dlamda, w, S, lds, info, ioffd, ioffq, ioffdlamda, ioffw, ioffs ) {
  info.setValue( 0 );
  if ( k < 0 ) info.setValue( -1 );
  else if ( kstart < 1 || kstart > Math.max( 1, k ) ) info.setValue( -2 );
  else if ( Math.max( 1, kstop ) < kstart || kstop > Math.max( 1, k ) )
  {
    info.setValue( -3 );
  } else if ( n < k ) info.setValue( -4 );
  else if ( ldq < Math.max( 1, k ) ) info.setValue( -7 );
  else if ( lds < Math.max( 1, k ) ) info.setValue( -12 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlaed9', - info.getValue() );
    return;
  }
  if ( k == 0 ) return;
  for ( var i = 1; i <= n; i ++ ) {
    dlamda[ ioffdlamda + i - 1 ] =
      LaPack0.dlamc3( dlamda[ ioffdlamda + i - 1 ] ,
      dlamda[ ioffdlamda + i - 1 ] ) - dlamda[ ioffdlamda + i - 1 ];
  }
  for ( var j = kstart; j <= kstop; j ++ ) {
    var dlamReference =
      new NumberReference( d[ ioffd + j - 1 ] );
    LaPack2.dlaed4( k, j, dlamda, w, Q, rho, dlam, info, ioffdlamda,
      ioffw, ioffq + ( j - 1 ) * ldq );
    d[ ioffd + j - 1 ] = dlam.getValue();
  }
  if ( k == 1 || k == 2 ) {
    for ( i = 1; i <= k; i ++ ) {
      for ( j = 1; j <= k; j ++ ) {
        S[ ioffs + j - 1 + ( i - 1 ) * lds ] =
          Q[ ioffq + j - 1 + ( i - 1 ) * ldq ];
      }
    }
    return;
  }
  Blas1.dcopy( k, w, 1, S, 1, ioffw, ioffs );
  Blas1.dcopy( k, Q, ldq + 1, w, 1, ioffq, ioffw );
  for ( j = 1; j <= k; j ++ ) {
    for ( i = 1; i <= j - 1; i ++ ) {
      w[ ioffw + i - 1 ] *= Q[ ioffq + i - 1 + ( j - 1 ) * ldq ]
        / ( dlamda[ ioffdlamda + i - 1 ]
        - dlamda[ ioffdlamda + j - 1 ] );
    }
    for ( i = j + 1; i <= k; i ++ ) {
      w[ ioffw + i - 1 ] *= Q[ ioffq + i - 1 + ( j - 1 ) * ldq ]
        / ( dlamda[ ioffdlamda + i - 1 ]
        - dlamda[ ioffdlamda + j - 1 ] );
    }
  }
  for ( i = 1; i <= k; i ++ ) {
    w[ ioffw + i - 1 ] =
      ( S[ ioffs + i - 1 ] >= 0 ? Math.sqrt( - w[ ioffw + i - 1 ] ) :
      - Math.sqrt( - w[ ioffw + i - 1 ] ) );
  }
  for ( j = 1; j <= k; j ++ ) {
    for ( i = 1; i <= k; i ++ ) {
      Q[ ioffq + i - 1 + ( j - 1 ) * ldq ] = w[ ioffw + i - 1 ]
        / Q[ ioffq + i - 1 + ( j - 1 ) * ldq ];
    }
    var temp = Blas1.dnrm2( k, Q, 1, ioffq + ( j - 1 ) * ldq );
    for ( i = 1; i <= k; i ++ ) {
      S[ ioffs + i - 1 + ( j - 1 ) * lds ] =
        Q[ ioffq + i - 1 + ( j - 1 ) * ldq ] / temp;
    }
  }
}
//*************************************************************************
LaPack3.dlaexc = function( wantq, n, T, ldt, Q, ldq, j1,
n1, n2, work, info, iofft, ioffq, ioffwork ) {
  var ldd = 4;
  var ldx = 2;
  var D = new Array( ldd * 4 );
  var u = new Array( 3 );
  var u1 = new Array( 3 );
  var u2 = new Array( 3 );
  var X = new Array( ldx * 2 );
  info.setValue( 0 );
  if ( n == 0 || n1 == 0 || n2 == 0 ) return;
  if ( j1 + n1 > n ) return;
  var j2 = j1 + 1;
  var j3 = j1 + 2;
  var j4 = j1 + 3;
  var csReference = new NumberReference();
  var snReference = new NumberReference();
  var tempReference = new NumberReference();
  var scaleReference = new NumberReference(); 
  var xnormReference = new NumberReference(); 
  var alphaReference = new NumberReference();
  var tauReference = new NumberReference();
  if ( n1 == 1 && n2 == 1 ) {
    var t11 = T[ iofft + j1 - 1 + ( j1 - 1 ) * ldt ];
    var t22 = T[ iofft + j2 - 1 + ( j2 - 1 ) * ldt ];
    LaPack1.dlartg( T[ iofft + j1 - 1 + ( j2 - 1 ) * ldt ],
      t22 - t11, cs, sn, temp );
    if ( j3 <= n ) {
      Blas1.drot( n - j1 - 1, T, ldt, T, ldt, cs.getValue(), sn.getValue(),
        iofft + j1 - 1 + ( j3 - 1 ) * ldt,
        iofft + j2 - 1 + ( j3 - 1 ) * ldt );
    }
    Blas1.drot( j1 - 1, T, 1, T, 1, cs.getValue(), sn.getValue(),
      iofft + ( j1 - 1 ) * ldt, iofft + ( j2 - 1 ) * ldt );
    T[ iofft + j1 - 1 + ( j1 - 1 ) * ldt ] = t22;
    T[ iofft + j2 - 1 + ( j2 - 1 ) * ldt ] = t11;
    if ( wantq ) {
      Blas1.drot( n, Q, 1, Q, 1, cs.getValue(), sn.getValue(),
        ioffq + ( j1 - 1 ) * ldq, ioffq + ( j2 - 1 ) * ldq );
    }
  } else {
    var nd = n1 + n2;
    LaPack0.dlacpy( 'Full', nd, nd, T, ldt, D, ldd,
      iofft + j1 - 1 + ( j1 - 1 ) * ldt, 0 );
    var dnorm = LaPack1.dlange( 'Max', nd, nd, D, ldd, work,
      0, ioffwork );
    var eps = LaPack0.dlamch( 'P' );
    var smlnum = LaPack0.dlamch( 'S' ) / eps;
    var thresh = Math.max( 10. * eps * dnorm, smlnum );
    var ierr = new IntReference();
    LaPack1.dlasy2( false, false, -1, n1, n2, D, ldd, D, ldd, D, ldd,
      scale, X, ldx, xnorm, ierr, 0, n1 + n1 * ldd, n1 * ldd, 0 );
    var k = n1 + n1 + n2 - 3;
    if ( k != 2 && k != 3 ) { // n1 = 1, n2 = 2
      u[ 0 ] = scale.getValue();
      u[ 1 ] = X[ 0 ];
      u[ 2 ] = X[ ldx ];
      alpha.setValue( u[ 2 ] );
      LaPack1.dlarfg( 3, alpha, u, 1, tau, 0 );
      u[2] = alpha.getValue();
      u[ 2 ] = 1.;
      t11 = T[ iofft + j1 - 1 + ( j1 - 1 ) * ldt ];
      LaPack2.dlarfx( 'L', 3, 3, u, tau.getValue(), D, ldd, work, 0, 0,
        ioffwork );
      LaPack2.dlarfx( 'R', 3, 3, u, tau.getValue(), D, ldd, work, 0, 0,
        ioffwork );
      if ( Math.max( Math.max( Math.abs( D[ 2 ] ),
      Math.abs( D[ 2 + ldd ] ) ),
      Math.abs( D[ 2 + 2 * ldd ] - t11 ) ) > thresh ) {
        info.setValue( 1 );
        return;
      }
      LaPack2.dlarfx( 'L', 3, n - j1 + 1, u, tau.getValue(), T, ldt, work,
       0, iofft + j1 - 1 + ( j1 - 1 ) * ldt, ioffwork );
      LaPack2.dlarfx( 'R', j2, 3, u, tau.getValue(), T, ldt, work, 0, 
        iofft + ( j1 - 1 ) * ldt, ioffwork );
      T[ iofft + j3 - 1 + ( j1 - 1 ) * ldt ] = 0.;
      T[ iofft + j3 - 1 + ( j2 - 1 ) * ldt ] = 0.;
      T[ iofft + j3 - 1 + ( j3 - 1 ) * ldt ] = t11;
      if ( wantq ) {
        LaPack2.dlarfx( 'R', n, 3, u, tau.getValue(), Q, ldq, work, 0, 
          ioffq + ( j1 - 1 ) * ldq, ioffwork );
      }
    } else if ( k == 2 ) { // 20
      u[ 0 ] = - X[ 0 ];
      u[ 1 ] = - X[ 1 ];
      u[ 2 ] = scale.getValue();
      alpha.setValue( u[ 0 ] );
      LaPack1.dlarfg( 3, alpha, u, 1, tau, 1 );
      u[ 0 ] = alpha.getValue();
      u[ 0 ] = 1.;
      var t33 = T[ iofft + j3 - 1 + ( j3 - 1 ) * ldt ];
      LaPack2.dlarfx( 'L', 3, 3, u, tau.getValue(), D, ldd, work, 0, 0,
        ioffwork );
      LaPack2.dlarfx( 'R', 3, 3, u, tau.getValue(), D, ldd, work, 0, 0,
        ioffwork );
      if ( Math.max( Math.max( Math.abs( D[ 1 ] ),
      Math.abs( D[ 2 ] ) ), Math.abs( D[ 0 ] - t33 ) ) > thresh ) {
        info.setValue( 1 );
        return;
      }
      LaPack2.dlarfx( 'R', j3, 3, u, tau.getValue(), T, ldt, work, 0, 
        iofft + ( j1 - 1 ) * ldt, ioffwork );
      LaPack2.dlarfx( 'L', 3, n - j1, u, tau.getValue(), T, ldt, work,
       0, iofft + j1 - 1 + ( j2 - 1 ) * ldt, ioffwork );
      T[ iofft + j1 - 1 + ( j1 - 1 ) * ldt ] = t33;
      T[ iofft + j2 - 1 + ( j1 - 1 ) * ldt ] = 0.;
      T[ iofft + j3 - 1 + ( j1 - 1 ) * ldt ] = 0.;
      if ( wantq ) {
        LaPack2.dlarfx( 'R', n, 3, u, tau.getValue(), Q, ldq, work, 0, 
          ioffq + ( j1 - 1 ) * ldq, ioffwork );
      }
    } else if ( k == 3 ) { // 30
      u1[ 0 ] = - X[ 0 ];
      u1[ 1 ] = - X[ 1 ];
      u1[ 2 ] = scale.getValue();
      alpha.setValue( u1[ 0 ] );
      var tau1Reference = new NumberReference();
      LaPack1.dlarfg( 3, alpha, u1, 1, tau1, 1 );
      u1[ 0 ] = 1.;
      temp.setValue(
        - tau1.getValue() * ( X[ ldx ] + u1[ 1 ] * X[ 1 + ldx ] ) );
      u2[ 0 ] = - temp.getValue() * u1[ 1 ] - X[ 1 + ldx ];
      u2[ 1 ] = - temp.getValue() * u1[ 2 ];
      u2[ 2 ] = scale.getValue();
      alpha.setValue( u2[ 0 ] );
      var tau2Reference = new NumberReference();
      LaPack1.dlarfg( 3, alpha, u2, 1, tau2, 1 );
      u2[ 0 ] = 1.;
      LaPack2.dlarfx( 'L', 3, 4, u1, tau1.getValue(), D, ldd, work, 0, 0,
        ioffwork );
      LaPack2.dlarfx( 'R', 4, 3, u1, tau1.getValue(), D, ldd, work, 0, 0,
        ioffwork );
      LaPack2.dlarfx( 'L', 3, 4, u2, tau2.getValue(), D, ldd, work, 0, 1,
        ioffwork );
      LaPack2.dlarfx( 'R', 4, 3, u2, tau2.getValue(), D, ldd, work, 0, ldd,
        ioffwork );
      if ( Math.max( Math.max( Math.abs( D[ 2 ] ),
      Math.abs( D[ 2 + ldd ] ) ), Math.max( Math.abs( D[ 3 ] ),
      Math.abs( D[ 3 + ldd ] ) ) ) > thresh ) {
        info.setValue( 1 );
        return;
      }
      LaPack2.dlarfx( 'L', 3, n - j1 + 1, u1, tau1.getValue(), T, ldt,
        work, 0, iofft + j1 - 1 + ( j1 - 1 ) * ldt, ioffwork );
      LaPack2.dlarfx( 'R', j4, 3, u1, tau1.getValue(), T, ldt, work, 0, 
        iofft + ( j1 - 1 ) * ldt, ioffwork );
      LaPack2.dlarfx( 'L', 3, n - j1 + 1, u2, tau2.getValue(), T, ldt,
        work, 0, iofft + j2 - 1 + ( j1 - 1 ) * ldt, ioffwork );
      LaPack2.dlarfx( 'R', j4, 3, u2, tau2.getValue(), T, ldt, work, 0, 
        iofft + ( j2 - 1 ) * ldt, ioffwork );
      T[ iofft + j3 - 1 + ( j1 - 1 ) * ldt ] = 0.;
      T[ iofft + j3 - 1 + ( j2 - 1 ) * ldt ] = 0.;
      T[ iofft + j4 - 1 + ( j1 - 1 ) * ldt ] = 0.;
      T[ iofft + j4 - 1 + ( j2 - 1 ) * ldt ] = 0.;
      if ( wantq ) {
        LaPack2.dlarfx( 'R', n, 3, u1, tau1.getValue(), Q, ldq, work, 0, 
          ioffq + ( j1 - 1 ) * ldq, ioffwork );
        LaPack2.dlarfx( 'R', n, 3, u2, tau2.getValue(), Q, ldq, work, 0, 
          ioffq + ( j2 - 1 ) * ldq, ioffwork );
      }
    } // 40
    var aReference = new NumberReference();
    var bReference = new NumberReference();
    var cReference = new NumberReference();
    var dReference = new NumberReference();
    var wr1Reference = new NumberReference();
    var wi1Reference = new NumberReference();
    var wr2Reference = new NumberReference();
    var wi2Reference = new NumberReference();
    if ( n2 == 2 ) {
      a.setValue( T[ iofft + j1 - 1 + ( j1 - 1 ) * ldt ] );
      b.setValue( T[ iofft + j1 - 1 + ( j2 - 1 ) * ldt ] );
      c.setValue( T[ iofft + j2 - 1 + ( j1 - 1 ) * ldt ] );
      d.setValue( T[ iofft + j2 - 1 + ( j2 - 1 ) * ldt ] );
      LaPack1.dlanv2( a, b, c, d, wr1, wi1, wr2, wi2, cs, sn );
      T[ iofft + j1 - 1 + ( j1 - 1 ) * ldt ] = a.getValue();
      T[ iofft + j1 - 1 + ( j2 - 1 ) * ldt ] = b.getValue();
      T[ iofft + j2 - 1 + ( j1 - 1 ) * ldt ] = c.getValue();
      T[ iofft + j2 - 1 + ( j2 - 1 ) * ldt ] = d.getValue();
      Blas1.drot( n - j1 - 1, T, ldt, T, ldt, cs.getValue(), sn.getValue(),
        iofft + j1 - 1 + ( j1 + 1 ) * ldt,
        iofft + j2 - 1 + ( j1 + 1 ) * ldt );
      Blas1.drot( j1 - 1, T, 1, T, 1, cs.getValue(), sn.getValue(),
        iofft + ( j1 - 1 ) * ldt, iofft + ( j2 - 1 ) * ldt );
      if ( wantq ) {
        Blas1.drot( n, Q, 1, Q, 1, cs.getValue(), sn.getValue(),
          ioffq + ( j1 - 1 ) * ldq, ioffq + ( j2 - 1 ) * ldq );
      }
    }
    if ( n1 == 2 ) {
      j3 = j1 + n2;
      j4 = j3 + 1;
      a.setValue( T[ iofft + j3 - 1 + ( j3 - 1 ) * ldt ] );
      b.setValue( T[ iofft + j3 - 1 + ( j4 - 1 ) * ldt ] );
      c.setValue( T[ iofft + j4 - 1 + ( j3 - 1 ) * ldt ] );
      d.setValue( T[ iofft + j4 - 1 + ( j4 - 1 ) * ldt ] );
      LaPack1.dlanv2( a, b, c, d, wr1, wi1, wr2, wi2, cs, sn );
      T[ iofft + j3 - 1 + ( j3 - 1 ) * ldt ] = a.getValue();
      T[ iofft + j3 - 1 + ( j4 - 1 ) * ldt ] = b.getValue();
      T[ iofft + j4 - 1 + ( j3 - 1 ) * ldt ] = c.getValue();
      T[ iofft + j4 - 1 + ( j4 - 1 ) * ldt ] = d.getValue();
      if ( j3 + 2 <= n ) {
        Blas1.drot( n - j3 - 1, T, ldt, T, ldt, cs.getValue(), sn.getValue(),
          iofft + j3 - 1 + ( j3 + 1 ) * ldt,
          iofft + j4 - 1 + ( j3 + 1 ) * ldt );
      }
      Blas1.drot( j3 - 1, T, 1, T, 1, cs.getValue(), sn.getValue(),
        iofft + ( j3 - 1 ) * ldt, iofft + ( j4 - 1 ) * ldt );
      if ( wantq ) {
        Blas1.drot( n, Q, 1, Q, 1, cs.getValue(), sn.getValue(),
          ioffq + ( j3 - 1 ) * ldq, ioffq + ( j4 - 1 ) * ldq );
      }
    }
  }
}
//*************************************************************************
LaPack3.dlalsa = function( icompq, smlsiz, n, nrhs, B, ldb, Bx,
ldbx, U, ldu, Vt, k, difl, difr, Z, poles, givptr, givcol, ldgcol, perm,
givnum, c, s, work, iwork, info, ioffb, ioffbx, ioffu, ioffvt, ioffk,
ioffdifl, ioffdifr, ioffz, ioffpoles, ioffgivptr, ioffgivcol, ioffperm,
ioffgivnum, ioffc, ioffs, ioffwork, ioffiwork) {
  throw new Error("not tested");
  info.setValue( 0 );
  if ( icompq < 0 || icompq > 1 ) info.setValue( -1 );
  else if ( smlsiz < 3 ) info.setValue( -2 );
  else if ( n < smlsiz ) info.setValue( -3 );
  else if ( nrhs < 1 ) info.setValue( -4 );
  else if ( ldb < n ) info.setValue( -6 );
  else if ( ldbx < n ) info.setValue( -8 );
  else if ( ldu < n ) info.setValue( -10 );
  else if ( ldgcol < n ) info.setValue( -19 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlalsa', - info.getValue() );
    return;
  }
  var inode = 1;
  var ndiml = inode + n;
  var ndimr = ndiml + n;
  LaPack0.dlasdt( n, nlvl, nd, iwork, iwork, iwork, smlsiz,
    ioffiwork + inode - 1, ioffiwork + ndiml - 1,
    ioffwork + ndimr - 1 );
  if ( icompq != 1 ) {
    var ndb1 = ( nd + 1 ) / 2;
    for ( var i = ndb1; i <= nd; i ++ ) {
      var i1 = i - 1;
      var ic = iwork[ ioffwork + inode + i1 - 1 ];
      var nl = iwork[ ioffiwork + ndiml + i1 - 1 ];
      var nr = iwork[ ioffiwork + ndimr + i1 - 1 ];
      var nlf = ic - nl;
      var nrf = ic + 1;
      Blas3.dgemm( 'T', 'N', nl, nrhs, nl, 1., U, ldu, B, ldb, 0.,
        Bx, ldbx, ioffu + nlf - 1, ioffb + nlf - 1, ioffbx + nlf - 1 );
      Blas3.dgemm( 'T', 'N', nr, nrhs, nr, 1., U, ldu, B, ldb, 0.,
        Bx, ldbx, ioffu + nrf - 1, ioffb + nrf - 1, ioffbx + nrf - 1 );
    } // 10
    for ( i = 1; i <= nd; i ++ ) {
      ic = iwork[ ioffiwork + inode + i - 2 ];
      Blas1.dcopy( nrhs, B, ldb, Bx, ldbx, ioffb + ic - 1,
        ioffbx + ic - 1 );
    } // 20
    var j = Math.pow( 2, nlvl );
    var sqre = 0;
    for ( var lvl = nlvl; lvl >= 1; lvl -- ) {
      var lvl2 = 2 * lvl - 1;
      if ( lvl == 1 ) {
        var lf = 1;
        var ll = 1;
      } else {
        lf = Math.pow( 2, lvl - 1 );
        ll = 2 * lf - 1;
      }
      for ( i = lf; i <= ll; i ++ ) {
        var im1 = i - 1;
        ic = iwork[ ioffiwork + inode + im1 - 1 ];
        nl = iwork[ ioffiwork + ndiml + im1 - 1 ];
        nr = iwork[ ioffiwork + ndimr + im1 - 1 ];
        nlf = ic - nl;
        nrf = ic + 1;
        j --;
        LaPack2.dlals0( icompq, nl, nr, sqre, nrhs, Bx, ldbx, B, ldb,
          perm, givptr[ ioffgivptr + j - 1 ], 
          givcol[ ioffgivcol + nlf - 1 + ( lvl2 - 1 ) * ldgcol ],
          ldgcol, givnum, ldu, poles, difl, difr,
          Z, k[ ioffk + j - 1 ], c[ ioffc + j - 1 ],
          s[ ioffs + j - 1 ], work, info,
          ioffbx + nlf - 1, ldb + nlf - 1,
          ioffperm + nlf - 1 + ( lvl - 1 ) * ldgcol,
          ioffgivnum + nlf - 1 + ( lvl2 - 1 ) * ldu,
          ioffpoles + nlf - 1 + ( lvl2 - 1 ) * ldu,
          ioffdifl + nlf - 1 + ( lvl - 1 ) * ldu,
          ioffdifr + nlf - 1 + ( lvl2 - 1 ) * ldu,
          ioffz + nlf - 1 + ( lvl - 1 ) * ldu, ioffwork );
      } // 30
    } // 40
  } else {
    j = 0;
    for ( lvl = 1; lvl <= nlvl; lvl ++ ) {
      lvl2 = 2 * lvl - 1;
      if ( lvl == 1 ) {
        lf = 1;
        ll = 1;
      } else {
        lf = Math.pow( 2, lvl - 1 );
        ll = 2 * lf - 1;
      }
      for ( i = ll; i >= lf; i -- ) {
        im1 = i - 1;
        ic = iwork[ ioffiwork + inode + im1 - 1 ];
        nl = iwork[ ioffiwork + ndiml + im1 - 1 ];
        nr = iwork[ ioffiwork + ndimr + im1 - 1 ];
        nlf = ic - nl;
        nrf = ic + 1;
        sqre = ( i == ll ? 0 : 1 );
        j ++;
        LaPack2.dlals0( icompq, nl, nr, sqre, nrhs, B, ldb, Bx, ldbx,
          perm, givptr[ ioffgivptr + j - 1 ], 
          givcol[ ioffgivcol + nlf - 1 + ( lvl2 - 1 ) * ldgcol ],
          ldgcol, givnum, ldu, poles, difl, difr,
          Z, k[ ioffk + j - 1 ], c[ ioffc + j - 1 ],
          s[ ioffs + j - 1 ], work, info,
          ioffb + nlf - 1, ldbx + nlf - 1,
          ioffperm + nlf - 1 + ( lvl - 1 ) * ldgcol,
          ioffgivnum + nlf - 1 + ( lvl2 - 1 ) * ldu,
          ioffpoles + nlf - 1 + ( lvl2 - 1 ) * ldu,
          ioffdifl + nlf - 1 + ( lvl - 1 ) * ldu,
          ioffdifr + nlf - 1 + ( lvl2 - 1 ) * ldu,
          ioffz + nlf - 1 + ( lvl - 1 ) * ldu, ioffwork );
      } // 60
    } // 70
    var ndb1 = ( nd + 1 ) / 2;
    for ( i = ndb1; i <= nd; i ++ ) {
      var i1 = i - 1;
      ic = iwork[ ioffiwork + inode + i1 - 1 ];
      nl = iwork[ ioffiwork + ndiml + i1 - 1 ];
      nr = iwork[ ioffiwork + ndimr + i1 - 1 ];
      nlp1 = nl + 1;
      nrp1 = ( i == nd ? nr : nr + 1 );
      nlf = ic - nl;
      nrf = ic + 1;
      Blas3.dgemm( 'T', 'N', nlp1, nrhs, nlp1, 1., Vt, ldu, B, ldb, 0.,
        Bx, ldbx, ioffvt + nlf - 1, ioffb + nlf - 1,
        ioffbx + nlf - 1 );
      Blas3.dgemm( 'T', 'N', nrp1, nrhs, nrp1, 1., Vt, ldu, B, ldb, 0.,
        Bx, ldbx, ioffvt + nrf - 1, ioffb + nrf - 1,
        ioffbx + nrf - 1 );
    }
  } // 90
}
//*************************************************************************
LaPack3.dlasd3 = function( nl, nr, sqre, k, d, Q, ldq, dsigma,
U, ldu, U2, ldu2, Vt, ldvt, Vt2, ldvt2, idxc, ctot, z, info, ioffd, ioffq,
ioffdsigma, ioffu, ioffu2, ioffvt, ioffvt2, ioffidxc, ioffctot, ioffz ) {
  throw new Error("not tested");
  info.setValue( 0 );
  if ( nl < 1 ) info.setValue( -1 );
  else if ( nr < 1 ) info.setValue( -2 );
  else if ( sqre != 1 && sqre != 0 ) info.setValue( -3 );
  var n = nl + nr + 1;
  var m = n + sqre;
  var nlp1 = nl + 1;
  var nlp2 = nl + 2;
  if ( k < 1 || k > n ) info.setValue( -4 );
  else if ( ldq < k ) info.setValue( -7 );
  else if ( ldu < n ) info.setValue( -10 );
  else if ( ldu2 < n ) info.setValue( -12 );
  else if ( ldvt < m ) info.setValue( -14 );
  else if ( ldvt2 < m ) info.setValue( -16 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlasd3', - info.getValue() );
    return;
  }
  if ( k == 1 ) {
    d[ ioffd ] = Math.abs( z[ ioffz ] );
    Blas1.dcopy( m, Vt2, ldvt2, Vt, ldvt, ioffvt2, ioffvt );
    if ( z[ ioffz ] > 0. ) {
      Blas1.dcopy( n, U2, 1, U, 1, ioffu2, ioffu );
    } else {
      for ( i = 1; i <= n; i ++ ) {
        U[ ioffu + i - 1 ] = - U2[ ioffu2 + i - 1 ];
      }
    }
    return;
  }
  for ( i = 1; i <= k; i ++ ) {
    dsigma[ ioffdsigma + i - 1 ] =
      LaPack0.dlamc3( dsigma[ ioffdsigma + i - 1 ],
      dsigma[ ioffdsigma + i - 1 ] ) - dsigma[ ioffdsigma + i - 1 ];
  }
  Blas1.dcopy( k, z, 1, Q, 1, ioffz, ioffq );
  var rho = Blas1.dnrm2( k, z, 1, ioffz );
  LaPack1.dlascl( 'G', 0, 0, rho, 1., k, 1, z, k, info, ioffz );
  rho = rho * rho;
  for ( var j = 1; j <= k; j ++ ) {
    var sigmaReference =
      new NumberReference( d[ ioffd + j - 1 ] );
    LaPack2.dlasd4( k, j, dsigma, z, U, rho, sigma, Vt, info,
      ioffdsigma, ioffz, ioffu + ( j - 1 ) * ldu,
      ioffvt + ( j - 1 ) * ldvt );
    d[ ioffd + j - 1 ] = sigma.getValue();
    if ( info.getValue() != 0 ) return;
  } // 30
  for ( var i = 1; i <= k; i ++ ) {
    z[ ioffz + i - 1 ] = U[ ioffu + i - 1 + ( k - 1 ) * ldu ]
      * Vt[ ioffvt + i - 1 + ( k - 1 ) * ldvt ];
    for ( j = 1; j <= i - 1; j ++ ) {
      z[ ioffz + i - 1 ] *= U[ ioffu + i - 1 + ( j - 1 ) * ldu ]
        * Vt[ ioffvt + i - 1 + ( j - 1 ) * ldvt ]
        / ( dsigma[ ioffdsigma + i - 1 ]
          - dsigma[ ioffdsigma + j - 1 ] )
        / ( dsigma[ ioffdsigma + i - 1 ]
          + dsigma[ ioffdsigma + j - 1 ] );
    } // 40
    for ( j = i; j <= k - 1; j ++ ) {
      z[ ioffz + i - 1 ] *= U[ ioffu + i - 1 + ( j - 1 ) * ldu ]
        * Vt[ ioffvt + i - 1 + ( j - 1 ) * ldvt ]
        / ( dsigma[ ioffdsigma + i - 1 ]
          - dsigma[ ioffdsigma + j ] )
        / ( dsigma[ ioffdsigma + i - 1 ]
          + dsigma[ ioffdsigma + j ] );
    } // 50
    z[ ioffz + i - 1 ] = ( Q[ ioffq + i - 1 ] >= 0. ?
      Math.sqrt( Math.abs( z[ ioffz + i - 1 ] ) ) :
      - Math.sqrt( Math.abs( z[ ioffz + i - 1 ] ) ) );
  } // 60
  for ( i = 1; i <= k; i ++ ) {
    Vt[ ioffvt + ( i - 1 ) * ldvt ] = z[ ioffz ]
      / U[ ioffu + ( i - 1 ) * ldu ]
      / Vt[ ioffvt + ( i - 1 ) * ldvt ];
    U[ ioffu + ( i - 1 ) * ldu ] = -1.;
    for ( j = 2; j <= k; j ++ ) {
      Vt[ ioffvt + j - 1 + ( i - 1 ) * ldvt ] = z[ ioffz + j - 1 ]
        / U[ ioffu + j - 1 + ( i - 1 ) * ldu ]
        / Vt[ ioffvt + j - 1 + ( i - 1 ) * ldvt ];
      U[ ioffu + j - 1 + ( i - 1 ) * ldu ] =
        dsigma[ ioffdsigma + j - 1 ]
        * Vt[ ioffvt + j - 1 + ( i - 1 ) * ldvt ];
    } // 70
    var temp = Blas1.dnrm2( k, U, 1, ioffu + ( i - 1 ) * ldu );
    Q[ ioffq + ( i - 1 ) * ldq ] = U[ ioffu + ( i - 1 ) * ldu ] / temp;
    for ( j = 2; j <= k; j ++ ) {
      var jc = idxc[ ioffidxc + j - 1 ];
      Q[ ioffq + j - 1 + ( i - 1 ) * ldq ] =
        U[ ioffu + jc - 1 + ( i - 1 ) * ldu ] / temp;
    } // 80
  } // 90
  if ( k == 2 ) {
    Blas3.dgemm( 'N', 'N', n, k, k, 1., U2, ldu2, Q, ldq, 0., U, ldu,
      ioffu2, ioffq, ioffu );
  } else {
    if ( ctot[ ioffctot ] > 0 ) {
      Blas3.dgemm( 'N', 'N', nl, k, ctot[ ioffctot ], 1., U2, ldu2,
        Q, ldq, 0., U, ldu, ioffu2 + ldu2, ioffq + 1, ioffu );
      if ( ctot[ ioffctot + 2 ] > 0 ) {
        var ktemp = 2 + ctot[ ioffctot ] + ctot[ ioffctot + 1 ];
        Blas3.dgemm( 'N', 'N', nl, k, ctot[ ioffctot + 2 ], 1.,
          U2, ldu2, Q, ldq, 1., U, ldu, ioffu2 + ( ktemp - 1 ) * ldu2,
          ioffq + ktemp - 1, ioffu );
      }
    } else if ( ctot[ ioffctot + 2 ] > 0 ) {
      ktemp = 2 + ctot[ ioffctot ] + ctot[ ioffctot + 1 ];
      Blas3.dgemm( 'N', 'N', nl, k, ctot[ ioffctot + 2 ], 1.,
        U2, ldu2, Q, ldq, 0., U, ldu, ioffu2 + ( ktemp - 1 ) * ldu2,
        ioffq + ktemp - 1, ioffu );
    } else {
      LaPack0.dlacpy( 'F', nl, k, U2, ldu2, U, ldu, ioffu2, ioffu );
    }
    Blas1.dcopy( k, Q, ldq, U, ldu, ioffq, ioffu + nlp1 - 1 );
    ktemp = 2 + ctot[ ioffctot ];
    var ctemp = ctot[ ioffctot + 1 ] + ctot[ ioffctot + 2 ];
    Blas3.dgemm( 'N', 'N', nr, k, ctemp, 1., U2, ldu2, Q, ldq, 0.,
      U, ldu, ioffu2 + nlp2 - 1 + ( ktemp - 1 ) * ldu2,
      ioffq + ktemp - 1, ioffu + nlp2 - 1 );
  }
  for ( i = 1; i <= k; i ++ ) {
    temp = Blas1.dnrm2( k, Vt, 1, ioffvt + ( i- 1 ) * ldvt );
    Q[ ioffq + i - 1 ] = Vt[ ioffvt + ( i - 1 ) * ldvt ] / temp;
    for ( j = 2; j <= k; j ++ ) {
      jc = idxc[ ioffidxc + j - 1 ];
      Q[ ioffq + i - 1 + ( j - 1 ) * ldq ] =
        Vt[ ioffvt + jc - 1 + ( i - 1 ) * ldvt ] / temp;
    } // 110
  } // 120
  if ( k == 2 ) {
    Blas3.dgemm( 'N', 'N', k, m, k, 1., Q, ldq, Vt2, ldvt2, 0.,
      Vt, ldvt, ioffq, ioffvt2, ioffvt );
    return;
  }
  ktemp = 1 + ctot[ ioffctot ];
  Blas3.dgemm( 'N', 'N', k, nlp1, ktemp, 1., Q, ldq, Vt2, ldvt2, 0.,
    Vt, ldvt, ioffq, ioffvt2, ioffvt );
  ktemp = 2 + ctot[ ioffctot ] + ctot[ ioffctot + 1 ];
  if ( ktemp <= ldvt2 ) {
    Blas3.dgemm( 'N', 'N', k, nlp1, ctot[ ioffctot + 2 ], 1., Q, ldq,
      Vt2, ldvt2, 1., Vt, ldvt, ioffq + ( ktemp - 1 ) * ldq,
      ioffvt2 + ktemp - 1, ioffvt );
  }
  ktemp = ctot[ ioffctot ] + 1;
  var nrp1 = nr + sqre;
  if ( ktemp > 1 ) {
    for ( i = 1; i <= k; i ++ ) {
      Q[ ioffq + i - 1 + ( ktemp - 1 ) * ldq ] = Q[ ioffq + i - 1 ];
    } // 130
    for ( i = nlp2; i <= m; i ++ ) {
      Vt2[ ioffvt2 + ktemp - 1 + ( i - 1 ) * ldvt2 ] =
        Vt2[ ioffvt2 + ( i - 1 ) * ldvt2 ];
    } // 140
  }
  ctemp = 1 + ctot[ ioffctot + 1 ] + ctot[ ioffctot + 2 ];
  Blas3.dgemm( 'N', 'N', k, nrp1, ctemp, 1., Q, ldq, Vt2, ldvt2, 0.,
    Vt, ldvt, ioffq + ( ktemp - 1 ) * ldq,
    ioffvt2 + ktemp - 1 + ( nlp2 - 1 ) * ldvt2,
    ioffvt + ( nlp2 - 1 ) * ldvt );
}
//*************************************************************************
LaPack3.dlasd8 = function( icompq, k, d, z, vf, vl, difl, Difr,
lddifr, dsigma, work, info, ioffd, ioffz, ioffvf, ioffvl, ioffdifl,
ioffdifr, ioffdsigma, ioffwork) {
  throw new Error("not tested");
  info.setValue( 0 );
  if ( icompq < 0 || icompq > 1 ) info.setValue( -1 );
  else if ( k < 1 ) info.setValue( -2 );
  else if ( lddifr < k ) info.setValue( -9 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dlasd8', - info.getValue() );
    return;
  }
  if ( k == 1 ) {
    d[ ioffd ] = Math.abs( z[ ioffz ] );
    difl[ ioffdifl ] = d[ ioffd ];
    if ( icompq == 1 ) {
      difl[ ioffdifl + 1 ] = 1.;
      Difr[ ioffdifr + lddifr ] = 1.;
    }
    return;
  }
  for ( var i = 1; i <= k; i ++ ) {
    dsigma[ ioffdsigma + i - 1 ] =
      LaPack0.dlamc3( dsigma[ ioffdsigma + i - 1 ],
        dsigma[ ioffdsigma + i - 1 ] ) - dsigma[ ioffdsigma + i - 1 ];
  } // 10  
  var iwk1 = 1;
  var iwk2 = iwk1 + k;
  var iwk3 = iwk2 + k;
  var iwk2i = iwk2 - 1;
  var iwk3i = iwk3 - 1;
  var rho = Blas1.dnrm2( k, z, 1, ioffz );
  LaPack1.dlascl( 'G', 0, 0, rho, 1., k, 1, z, k, info, ioffz );
  rho = rho * rho;
  LaPack0.dlaset( 'A', k, 1, 1., 1., work, k, ioffwork + iwk3 - 1 );
  var sigmaReference = new NumberReference();
  for ( var j = 1; j <= k; j ++ ) {
    sigma.setValue( d[ ioffd + j - 1 ] );
    LaPack2.dlasd4( k, j, dsigma, z, work, rho, sigma, work, info,
      ioffdsigma, ioffz, ioffwork + iwk1 - 1, ioffwork + iwk2 - 1 );
    d[ ioffd + j - 1 ] = sigma.getValue();
    if ( info.getValue() != 0 ) {
      Blas2.xerbla( 'dlasd4', - info.getValue() );
      return;
    }
    work[ ioffwork + iwk3i + j - 1 ] *=
      work[ ioffwork + j - 1 ] * work[ ioffwork + iwk2i + j - 1 ];
    difl[ ioffdifl + j - 1 ] = - work[ ioffwork + j - 1 ];
    Difr[ ioffdifr + j - 1 ] = - work[ ioffwork + j ];
    for ( i = 1; i <= j - 1; i ++ ) {
      work[ ioffwork + iwk3i + i - 1 ] *=
        work[ ioffwork + i - 1 ] * work[ ioffwork + iwk2i + i - 1 ]
        / ( dsigma[ ioffdsigma + i - 1 ]
          - dsigma[ ioffdsigma + j - 1 ] )
        / ( dsigma[ ioffdsigma + i - 1 ]
          + dsigma[ ioffdsigma + j - 1 ] );
    } // 20
    for ( i = j + 1; i <= k; i ++ ) {
      work[ ioffwork + iwk3i + i - 1 ] *=
        work[ ioffwork + i - 1 ] * work[ ioffwork + iwk2i + i - 1 ]
        / ( dsigma[ ioffdsigma + i - 1 ]
          - dsigma[ ioffdsigma + j - 1 ] )
        / ( dsigma[ ioffdsigma + i - 1 ]
          + dsigma[ ioffdsigma + j - 1 ] );
    } // 30
  } // 40
  for ( i = 1; i <= k; i ++ ) {
    z[ ioffz + i - 1 ] = ( z[ ioffz + i - 1 ] >= 0 ?
      Math.sqrt( Math.abs( work[ ioffwork + iwk3i + i - 1 ] ) ) :
      - Math.sqrt( Math.abs( work[ ioffwork + iwk3i + i - 1 ] ) ) );
  } // 50
  for ( j = 1; j <= k; j ++ ) {
    var diflj = difl[ ioffdifl + j - 1 ];
    var dj = d[ ioffd + j - 1 ];
    var dsigj = - dsigma[ ioffdsigma + j - 1 ];
    if ( j < k ) {
      var difrj = - Difr[ ioffdifr + j - 1 ];
      var dsigjp = - dsigma[ ioffdsigma + j ];
    }
    work[ ioffwork + j - 1 ] = - z[ ioffz + j - 1 ] / diflj
      / ( dsigma[ ioffdsigma + j - 1 ] + dj );
    for ( i = 1; i <= j - 1; i ++ ) {
      work[ ioffwork + i - 1 ] = z[ ioffz + i - 1 ]
        / ( LaPack0.dlamc3( dsigma[ ioffdsigma + i - 1 ], dsigj )
        - diflj ) / ( dsigma[ ioffdsigma + i - 1 ] + dj );
    } // 60
    for ( i = j + 1; i <= k; i ++ ) {
      work[ ioffwork + i - 1 ] = z[ ioffz + i - 1 ]
        / ( LaPack0.dlamc3( dsigma[ ioffdsigma + i - 1 ], dsigjp )
        + difrj ) / ( dsigma[ ioffdsigma + i - 1 ] + dj );
    } // 70
    var temp = Blas1.dnrm2( k, work, 1, ioffwork );
    work[ ioffwork + iwk2i + j - 1 ] =
      Blas1.ddot( k, work, 1, vf, 1, ioffwork, ioffvf ) / temp;
    work[ ioffwork + iwk3i + j - 1 ] =
      Blas1.ddot( k, work, 1, vl, 1, ioffwork, ioffvl ) / temp;
    if ( icompq == 1 ) {
      Difr[ ioffdifr + j - 1 + lddifr ] = temp;
    }
  } // 80
  Blas1.dcopy( k, work, 1, vf, 1, ioffwork + iwk2 - 1, ioffvf );
  Blas1.dcopy( k, work, 1, vl, 1, ioffwork + iwk3 - 1, ioffvl );
}
//*************************************************************************
LaPack3.dlasq2 = function( n, z, info, ioffz ) {
  throw new Error("not tested: complicated input");
  var iinfo = new IntReference();
  var cbias = 1.5;
  info.setValue( 0 );
  var eps = LaPack0.dlamch( 'Precision' );
  var safmin = LaPack0.dlamch( 'Safe minimum' );
  var tol = eps * 100.;
  var tol2 = tol * tol;
  if ( n < 0 ) {
    info.setValue( -1 );
    Blas2.xerbla( 'dlasq2', 1 );
    return;
  } else if ( n == 0 ) return;
  else if ( n == 1 ) {
    if ( z[ ioffz ] < 0. ) {
      info.setValue( -201 );
      Blas2.xerbla( 'dlasq2', 2 );
    }
    return;
  } else if ( n == 2 ) {
    if ( z[ ioffz + 1 ] < 0. || z[ ioffz + 2 ] < 0. ) {
      info.setValue( -2 );
      Blas2.xerbla( 'dlasq2', 2 );
      return;
    } else if ( z[ ioffz + 2 ] > z[ ioffz ] ) {
      var d = z[ ioffz + 2 ];
      z[ ioffz + 2 ] = z[ ioffz ];
      z[ ioffz ] = d;
    }
    z[ ioffz + 4 ] = z[ ioffz ] + z[ ioffz + 1 ] + z[ ioffz + 2 ];
    if ( z[ ioffz + 1 ] > z[ ioffz + 2 ] * tol2 ) {
      var t =
        0.5 * ( ( z[ ioffz ] - z[ ioffz + 2 ] ) + z[ ioffz + 1 ] );
      var s = z[ ioffz + 2 ] * ( z[ ioffz + 1 ] / t );
      if ( s <= t ) {
        s = z[ ioffz + 2 ] * ( z[ ioffz + 1 ] 
          / ( t * ( 1. + Math.sqrt( 1. + s / t ) ) ) );
      } else {
        s = z[ ioffz + 2 ] * ( z[ ioffz + 1 ]
          / ( t + Math.sqrt( t ) * Math.sqrt( t + s ) ) );
      }
      t = z[ ioffz ] + ( s + z[ ioffz + 1 ] );
      z[ ioffz + 2 ] *= z[ ioffz ] / t;
      z[ ioffz ] = t;
    }
    z[ ioffz + 1 ] = z[ ioffz + 2 ];
    z[ ioffz + 5 ] = z[ ioffz + 1 ] + z[ ioffz ];
    return;
  }
  z[ ioffz + 2 * n ] = 0.;
  var emin = z[ ioffz + 1 ];
  var qmax = 0.;
  var zmax = 0.;
  d = 0.;
  var e = 0.;
  for ( var k = 1; k <= 2 * ( n - 1 ); k += 2 ) {
    if ( z[ ioffz + k - 1 ] < 0. ) {
      info.setValue( - ( 200 + k ) );
      Blas2.xerbla( 'dlasq2', 2 );
      return;
    } else if ( z[ ioffz + k ] < 0. ) {
      info.setValue( - ( 200 + k  + 1 ) );
      Blas2.xerbla( 'dlasq2', 2 );
      return;
    }
    d += z[ ioffz + k - 1 ];
    e += z[ ioffz + k ];
    qmax = Math.max( qmax, z[ ioffz + k - 1 ] );
    emin = Math.min( emin, z[ ioffz + k ] );
    zmax = Math.max( Math.max( qmax, zmax ), z[ ioffz + k ] );
  } // 10
  if ( z[ ioffz + 2 * n - 2 ] < 0. ) {
    info.setValue( - ( 200 + 2 * n - 1 ) );
    Blas2.xerbla( 'dlasq2', 2 );
    return;
  }
  d += z[ ioffz + 2 * n - 2 ];
  qmax = Math.max( qmax, z[ ioffz + 2 * n - 2 ] );
  zmax = Math.max( qmax, zmax );
  if ( e == 0. ) {
    for ( k = 2; k <= n; k ++ ) {
      z[ ioffz + k - 1 ] = z[ ioffz + 2 * k - 2 ];
    } // 20
    LaPack0.dlasrt( 'D', n, z, iinfo, ioffz );
    z[ ioffz + 2 * n - 2 ] = d;
    return;
  }
  var trace = d + e;
  if ( trace == 0. ) {
    z[ ioffz + 2 * n - 2 ] = 0.;
    return;
  }
  var ieee =
    LaPack0.ilaenv( 10, 'dlasq2', 'N', 1, 2, 3, 4 ) == 1 &&
    LaPack0.ilaenv( 11, 'dlasq2', 'N', 1, 2, 3, 4 ) == 1;
  for ( k = 2 * n; k >= 2; k -= 2 ) {
    z[ ioffz + 2 * k - 1 ] = 0.;
    z[ ioffz + 2 * k - 2 ] = z[ ioffz + k - 1 ];
    z[ ioffz + 2 * k - 3 ] = 0.;
    z[ ioffz + 2 * k - 4 ] = z[ ioffz + k - 2 ];
  } // 30
  var i0 = 1;
  var n0 = n;
  if ( cbias * z[ ioffz + 4 * i0 - 4 ] < z[ ioffz + 4 * n0 - 4 ] ) {
    var ipn4 = 4 * ( i0 + n0 );
    for ( var i4 = 4 * i0; i4 <= 2 * ( i0 + n0 - 1 ); i4 += 4 ) {
      var temp = z[ ioffz + i4 - 4 ];
      z[ ioffz + i4 - 4 ] = z[ ioffz + ipn4 - i4 - 4 ];
      z[ ioffz + ipn4 - i4 - 4 ] = temp;
      temp = z[ ioffz + i4 - 2 ];
      z[ ioffz + i4 - 2 ] = z[ ioffz + ipn4 - i4 - 6 ];
      z[ ioffz + ipn4 - i4 - 6 ] = temp;
    } // 40
  }
  var pp = 0;
  for ( k = 1; k <= 2; k ++ ) {
    d = z[ ioffz + 4 * n0 + pp - 4 ];
    for ( i4 = 4 * ( n0 - 1 ) + pp; i4 >= 4 * i0 + pp; i4 -= 4 ) {
      if ( z[ ioffz + i4 - 2 ] <= tol2 * d ) {
        z[ ioffz + i4 - 2 ] = - 0.;
        d = z[ ioffz + i4 - 4 ];
      } else {
        d = z[ ioffz + i4 - 4 ] * ( d / ( d + z[ ioffz + i4 - 2 ] ) );
      }
    } // 50
    emin = z[ ioffz + 4 * i0 + pp ];
    d = z[ ioffz + 4 * i0 + pp - 4 ];
    for ( i4 = 4 * i0 + pp; i4 <= 4 * ( n0 - 1 ) + pp; i4 += 4 ) {
      z[ ioffz + i4 - 2 * pp - 3 ] = d + z[ ioffz + i4 - 2 ];
      if ( z[ ioffz + i4 - 2 ] <= tol2 * d ) {
        z[ ioffz + i4 - 2 ] = - 0.;
        z[ ioffz + i4 - 2 * pp - 3 ] = d;
        z[ ioffz - 2 * pp - 1 ] = 0.;
        d = z[ ioffz + i4 ];
      } else if ( safmin * z[ ioffz + i4 ] <
      z[ ioffz + i4 - 2 * pp - 3 ]
      && safmin * z[ ioffz + i4 - 2 * pp - 3 ] < z[ ioffz + i4 ] ) {
        temp = z[ ioffz + i4 ] / z[ ioffz + i4 - 2 * pp - 3 ];
        z[ ioffz + i4 - 2 * pp - 1 ] = z[ ioffz + i4 - 2 ] * temp;
        d *= temp;
      } else {
        z[ ioffz + i4 - 2 * pp - 1 ] = z[ ioffz + i4 ]
          * ( z[ ioffz + i4 - 2 ] / z[ ioffz + i4 - 2 * pp - 3 ] );
        d = z[ ioffz + i4 ] * ( d / z[ ioffz + i4 - 2 * pp - 3 ] );
      }
      emin = Math.min( emin, z[ ioffz + i4 - 2 * pp - 1 ] );
    } // 60
    z[ ioffz + 4 * n0 - pp - 3 ] = d;
    qmax = z[ ioffz + 4 * i0 - pp - 3 ];
    for ( i4 = 4 * i0 - pp + 2; i4 <= 4 * n0 - pp - 2; i4 += 4 ) {
      qmax = Math.max( qmax, z[ ioffz + i4 - 1 ] );
    } // 70
    pp = 1 - pp;
  } // 80
  var ttype = new IntReference( 0 );
  var dmin1Reference = new NumberReference( 0. );
  var dmin2Reference = new NumberReference( 0. );
  var dnReference = new NumberReference( 0. );
  var dn1Reference = new NumberReference( 0. );
  var dn2Reference = new NumberReference( 0. );
  var gReference = new NumberReference( 0. );
  var tauReference = new NumberReference( 0. );
  var iter = new IntReference( 2 );
  var nfail = new IntReference( 0 );
  var ndiv = new IntReference( 2 * ( n0 - i0 ) );
  var goto170 = false;
  for ( var iwhila = 1; iwhila <= n + 1; iwhila ++ ) {
    if ( n0 < 1 ) {
      goto170 = true;
      break;
    }
    var desigReference = new NumberReference( 0. );
    var sigmaReference = new NumberReference( ( n0 == n ?
      0. : -z[ ioffz + 4 * n0 - 2 ] ) );
    if ( sigma.getValue() < 0. ) {
      info.setValue( 1 );
      return;
    }
    var emax = 0.;
    emin = ( n0 > i0 ? Math.abs( z[ ioffz + 4 * n0 - 6 ] ) : 0. );
    var qmin = z[ ioffz + 4 * n0 - 4 ];
    qmax = qmin;
    var goto100 = false;
    for ( i4 = 4 * n0; i4 >= 8; i4 -= 4 ) {
      if ( z[ ioffz + i4 - 6 ] <= 0. ) {
        goto100 = true;
        break;
      }
      if ( qmin >= 4. * emax ) {
        qmin = Math.min( qmin, z[ ioffz + i4 - 4 ] );
        emax = Math.max( emax, z[ ioffz + i4 - 6 ] );
      }
      qmax =
        Math.max( qmax, z[ ioffz + i4 - 8 ] + z[ ioffz + i4 - 6 ] );
      emin = Math.min( emin, z[ ioffz + i4 - 6 ] );
    }
    if ( ! goto100 ) i4 = 4;
    i0 = i4 / 4;
    pp = 0;
    if ( n0 - i0 > 1 ) {
      var dee = z[ ioffz + 4 * i0 - 4 ];
      var deemin = dee;
      var kmin = i0;
      for ( i4 = 4 * i0 + 1; i4 <= 4 * n0 - 3; i4 += 4 ) {
        dee = z[ ioffz + i4 - 1 ]
          * ( dee / ( dee + z[ ioffz + i4 - 3 ] ) );
        if ( dee <= deemin ) {
          deemin = dee;
          kmin = ( i4 + 3 ) / 4 ;
        }
      } // 110
      if ( ( kmin - i0 ) * 2 < n0 - kmin &&
      deemin <= 0.5 * z[ ioffz + 4 * n0 - 4 ] ) {
        ipn4 = 4 * ( i0 + n0 );
        pp = 2;
        for ( i4 = 4 * i0; i4 <= 2 ( i0 + n0 - 1 ); i4 += 4 ) {
          temp = z[ ioffz + i4 - 4 ];
          z[ ioffz + i4 - 4 ] = z[ ioffz + ipn4 - i4 - 4 ];
          z[ ioffz + ipn4 - i4 - 4 ] = temp;
          temp = z[ ioffz + i4 - 3];
          z[ ioffz + i4 - 3] = z[ ioffz + ipn4 - i4 - 3 ];
          z[ ioffz + ipn4 - i4 - 3 ] = temp;
          temp = z[ ioffz + i4 - 2 ];
          z[ ioffz + i4 - 2 ] = z[ ioffz + ipn4 - i4 - 6];
          z[ ioffz + ipn4 - i4 - 6] = temp;
          temp = z[ ioffz + i4 - 1 ];
          z[ ioffz + i4 - 1 ] = z[ ioffz + ipn4 - i4 - 5];
          z[ ioffz + ipn4 - i4 - 5] = temp;
        } // 120
      }
    }
    var dminReference = new NumberReference( - Math.max( 0.,
      qmin - 2. * Math.sqrt( qmin ) * Math.sqrt( emax ) ) );
    var nbig = 100 * ( n0 - i0 + 1 );
    var goto150 = false;
    for ( var iwhilb = 1; iwhilb <= nbig; iwhilb ++ ) {
      if ( i0 > n0 ) {
        goto150 = true;
        break;
      }
      var n0ref = new IntReference( n0 );
      var ppref = new IntReference( pp );
      var qmaxrefReference = new NumberReference( qmax );
      LaPack2.dlasq3( i0, n0ref, z, ppref, dmin, sigma, desig, qmaxref,
        nfail, iter, ndiv, ieee, ttype, dmin1, dmin2, dn, dn1, dn2, g,
        tau, ioffz );
      n0 = n0ref.getValue();
      pp = ppref.getValue();
      qmax = qmaxref.getValue();
      pp = 1 - pp;
      if ( pp == 0 && n0 - i0 >= 3 ) {
        if ( z[ ioffz + 4 * n0 - 1 ] <= tol2 * qmax ||
        z[ ioffz + 4 * n0 - 2 ] <= tol2 * sigma.getValue() ) {
          var splt = i0 - 1;
          qmax = z[ ioffz + 4 * i0 - 4 ];
          emin = z[ ioffz + 4 * i0 - 2 ];
          var oldemn = z[ ioffz + 4 * i0 - 1 ];
          for ( i4 = 4 * i0; i4 <= 4 * ( n0 - 3 ); i4 += 4 ) {
            if ( z[ ioffz + i4 - 1 ] <= tol2 * z[ ioffz + i4 - 4 ] ||
            z[ ioffz + i4 - 2 ] <= tol2 * sigma.getValue() ) {
              z[ ioffz + i4 - 2 ] = - sigma.getValue();
              splt = i4 / 4;
              qmax = 0.;
              emin = z[ ioffz + i4 + 2];
              oldemn = z[ ioffz + i4 + 3 ];
            } else {
              qmax = Math.max( qmax, z[ ioffz + i4 ] );
              emin = Math.min( emin, z[ ioffz + i4 - 2 ] );
              oldemn = Math.min( oldemn, z[ ioffz + i4 - 1 ] )
            }
          } // 130
          z[ ioffz + 4 * n0 - 2 ] = emin;
          z[ ioffz + 4 * n0 - 1 ] = oldemn;
          i0 = splt + 1;
        }
      }
    } // 140
    if ( ! goto150 ) {
      info.setValue( 2 );
      var i1 = i0;
      var n1 = n0;
      while ( true ) { // 145
        var tempq = z[ ioffz + 4 * i0 - 4 ];
        z[ ioffz + 4 * i0 - 4 ] += sigma.getValue();
        for ( k = i0 + 1; k <= n0; k ++ ) {
          var tempe = z[ ioffz + 4 * k - 6 ];
          z[ ioffz + 4 * k - 6 ] *= tempq / z[ ioffz + 4 * k - 8 ];
          tempq = z[ ioffz + 4 * k - 4 ];
          z[ ioffz + 4 * k - 4 ] +=
            sigma.getValue() + tempe - z[ ioffz + 4 * k - 6 ];
        }
        if ( i1 > 1 ) {
          n1 = i1 - 1;
          while ( ( i1 >= 2 ) &&
          ( z[ ioffz + 4 * i1 - 6 ] >= 0. ) ) {
            i1 --;
          }
          sigma.setValue( - z[ ioffz + 4 * n1 - 2 ] );
          continue;
        }
        for ( k = 1; k <= n; k ++ ) {
          z[ ioffz + 2 * k - 2 ] = z[ ioffz + 4 * k - 4 ];
          z[ ioffz + 2 * k - 1 ] =
            ( k < n0 ? z[ ioffz + 4 * k - 2 ] : 0. );
        }
        return;
      }
    } // 150
  } // 160
  if ( ! goto170 ) {
    info.setValue( 3 );
    return;
  } // 170
  for ( k = 2; k <= n; k ++ ) {
    z[ ioffz + k - 1 ] = z[ ioffz + 4 * k - 4 ];
  }
  LaPack0.dlasrt( 'D', n, z, iinfo, ioffz );
  e = 0.;
  for ( k = n; k >= 1; k -- ) {
    e += z[ ioffz + k - 1 ];
  }
  z[ ioffz + 2 * n ] = trace;
  z[ ioffz + 2 * n + 1 ] = e;
  z[ ioffz + 2 * n + 2 ] = Number( iter );
  z[ ioffz + 2 * n + 3 ] = Number( ndiv ) / Number( n * n );
  z[ ioffz + 2 * n + 4 ] = 100. * nfail.getValue() / Number( iter );
}
/* old version:
LaPack3.dlasq2 = function( m, q, e, qq, ee, eps, tol2, small2,
supReference, kend, info, ioffq, ioffe, ioffqq, ioffee ) {
  var n = m;
  var off = 0;
  var off1 = off + 1;
  var sigmaReference = new NumberReference( 0. );
  var xinf = 0.;
  var iconv = new IntReference( 0 );
  var iphase = new IntReference( 2 );
  while ( n <= 2 ) {
    if ( ee[ ioffee + n - 2 ] <=
    Math.max( qq[ ioffqq + n - 1 ],
    Math.max( xinf, small2 ) ) * tol2 ) {
      q[ ioffq + n - 1 ] = qq[ ioffqq + n - 1 ];
      n = n - 1;
      if ( kend.getValue() > n ) kend.setValue( n );
      sup.setValue(
        Math.min( qq[ ioffqq + n - 1 ], qq[ ioffqq + n - 2 ] ) );
      continue;
    }
    if ( ee[ ioffee + n - 3 ] <= Math.max( Math.max( xinf, small2 ),
    ( qq[ ioffqq + n - 1 ] / ( qq[ ioffqq + n - 1 ]
    + ee[ ioffee + n - 2 ] + qq[ ioffqq + n - 2 ] ) )
    * qq[ ioffqq + n - 2 ] ) * tol2 ) {
      var qemax = Math.max(
        Math.max( qq[ ioffqq + n - 1 ], qq[ ioffqq + n - 2 ] ),
        ee[ ioffee + n - 2 ] );
      if ( qemax != 0. ) {
        if ( qemax == qq[ ioffqq + n - 2 ] ) {
          var xx =
            0.5 * ( qq[ ioffqq + n - 1 ] + qq[ ioffqq + n - 2 ]
            + ee[ ioffee + n - 2 ] + qemax * Math.sqrt(
            Math.pow( ( qq[ ioffqq + n - 1 ] - qq[ ioffqq + n - 2 ]
            + ee[ ioffee + n - 2 ] ) / qemax, 2 )
            + 4. * ee[ ioffee + n - 2 ] / qemax ) );
        } else if ( qemax == qq[ ioffqq + n - 1 ] ) {
          xx = 0.5 * ( qq[ ioffqq + n - 1 ] + qq[ ioffqq + n - 2 ]
            + ee[ ioffee + n - 2 ] + qemax * Math.sqrt( Math.pow(
            ( qq[ ioffqq + n - 2 ] - qq[ ioffqq + n - 1 ]
            + ee[ ioffee + n - 2 ] ) / qemax, 2 )
            + 4. * ee[ ioffee + n - 2 ] / qemax ) );
        } else {
          xx = 0.5 * ( qq[ ioffqq + n - 1 ] + qq[ ioffqq + n - 2 ]
            + ee[ ioffee + n - 2 ] + qemax * Math.sqrt( Math.pow(
            ( qq[ ioffqq + n - 1 ] - qq[ ioffqq + n - 2 ]
            + ee[ ioffee + n - 2 ] ) / qemax, 2 )
            + 4. * qq[ ioffqq + n - 2 ] / qemax ) );
        }
        var yy = ( Math.max( qq[ ioffqq + n - 1 ],
          qq[ ioffqq + n - 2 ] ) / xx )
          * Math.min( qq[ ioffqq + n - 1 ], qq[ ioffqq + n - 2 ] );
      } else {
       xx = 0.;
       yy = 0.;
      }
      q[ ioffq + n - 2 ] = xx;
      q[ ioffq + n - 1 ] = yy;
      n = n - 2;
      if ( kend.getValue() > n ) kend.setValue( n );
      sup.setValue( qq[ ioffqq + n - 1] );
      continue;
    } else break;
  }
  while ( true ) {
    if ( n == 0 ) {
      if ( off == 0 ) return;
      else {
        xinf = 0.;
        if ( ee[ ioffee + off - 1 ] > 0. ) {
          var isp = Math.round( ee[ ioffee + off - 1 ] );
          iphase.setValue( 1 );
        } else {
          isp = - Math.round( ee[ ioffee + off - 1 ] );
          iphase.setValue( 2 );
        }
        sigma.setValue( e[ ioffe + off - 1] );
        n = off - isp + 1;
        off1 = isp;
        off = off1 - 1;
        if ( n <= 2 ) continue;
        if ( iphase.getValue() == 1 ) {
          sup.setValue( Math.min( Math.min( q[ ioffq + n + off - 1 ],
            q[ ioffq + n - 2 + off ] ), q[ ioffq + n - 3 + off ] ) );
        } else {
          sup.setValue( Math.min( Math.min( qq[ ioffqq + n - 2 + off ],
            qq[ ioffqq + n - 3 + off ] ), qq[ ioffqq + n - 3 + off ] ) );
        }
        kend.setValue( 0 );
        iconv.setValue( -3 );
      }
      break;
    } else if ( n == 1 ) {
      if ( iphase.getValue() == 1 ) q[ ioffq + off1 - 1 ] += sigma.getValue();
      else {
        q[ ioffq + off1 - 1 ] = qq[ ioffqq + off1 - 1 ] + sigma.getValue();
      }
      n = 0;
      continue;
    } else if ( n == 2 ) {
      if ( iphase.getValue() == 2 ) {
        qemax = Math.max( Math.max( qq[ ioffqq + n + off - 1 ],
          qq[ ioffqq + n + off ] ), ee[ ioffee + n - 2 + off ] );
        if ( qemax != 0. ) {
          if ( qemax == qq[ ioffqq + n - 2 + off ] ) {
            xx = 0.5 ( qq[ ioffqq + n + off - 1 ]
              + qq[ ioffqq + n - 2 + off ]
              + ee[ ioffee + n - 2 + off ]
              + qemax * Math.sqrt( Math.pow(
              ( qq[ ioffqq + n + off - 1 ] - qq[ ioffqq + n - 2 + off ]
              + ee[ ioffee + n - 2 + off ] ) / qemax, 2 )
              + 4. * ee[ ioffee + off + n - 2 ] / qemax ) );
          } else if ( qemax == qq[ ioffqq + n + off - 1 ] ) {
            xx = 0.5 ( qq[ ioffqq + n + off - 1 ]
              + qq[ ioffqq + n - 2 + off ]
              + ee[ ioffee + n - 2 + off ]
              + qemax * Math.sqrt( Math.pow(
              ( qq[ ioffqq + n + off - 2 ] - qq[ ioffqq + n - 1 + off ]
              + ee[ ioffee + n - 2 + off ] ) / qemax, 2 )
              + 4. * ee[ ioffee + off + n - 2 ] / qemax ) );
          } else {
            xx = 0.5 ( qq[ ioffqq + n + off - 1 ]
              + qq[ ioffqq + n - 2 + off ]
              + ee[ ioffee + n - 2 + off ]
              + qemax * Math.sqrt( Math.pow(
              ( qq[ ioffqq + n + off - 1 ] - qq[ ioffqq + n - 2 + off ]
              + ee[ ioffee + n - 2 + off ] ) / qemax, 2 )
              + 4. * qq[ ioffqq + off + n - 2 ] / qemax ) );
          }
          yy = ( Math.max( qq[ ioffqq + n + off - 1 ],
            qq[ ioffqq + n - 2 + off ] ) / xx ) * Math.min(
            qq[ ioffqq + n + off - 1 ],
            qq[ ioffqq + n - 2 + off ] );
        } else {
          xx = 0.;
          yy = 0.;
        }
      } else {
        qemax = Math.max( Math.max( q[ ioffq + n + off - 1 ],
          q[ ioffq + n - 2 + off ] ), e[ ioffe + n - 2 + off ] );
        if ( qemax != 0. ) {
          if ( qemax == q[ ioffq + n - 2 + off ] ) {
            xx = 0.5 ( q[ ioffq + n + off - 1 ]
              + q[ ioffq + n - 2 + off ]
              + e[ ioffe + n - 2 + off ] + qemax * Math.sqrt( Math.pow(
              ( q[ ioffq + n + off - 1 ] - q[ ioffq + n - 2 + off ]
              + e[ ioffe + n - 2 + off ] ) / qemax, 2 )
              + 4. * e[ ioffe + off + n - 2 ] / qemax ) );
          } else if ( qemax == q[ ioffq + n + off - 1 ] ) {
            xx = 0.5 ( q[ ioffq + n + off - 1 ]
              + q[ ioffq + n - 2 + off ]
              + e[ ioffe + n - 2 + off ] + qemax * Math.sqrt( Math.pow(
              ( q[ ioffq + n + off - 2 ] - q[ ioffq + n - 1 + off ]
              + e[ ioffe + n - 2 + off ] ) / qemax, 2 )
              + 4. * e[ ioffe + off + n - 2 ] / qemax ) );
          } else {
            xx = 0.5 ( q[ ioffq + n + off - 1 ]
              + q[ ioffq + n - 2 + off ]
              + e[ ioffe + n - 2 + off ] + qemax * Math.sqrt( Math.pow(
              ( q[ ioffq + n + off - 1 ] - q[ ioffq + n - 2 + off ]
              + e[ ioffe + n - 2 + off ] ) / qemax, 2 )
              + 4. * q[ ioffq + off + n - 2 ] / qemax ) );
          }
          yy = ( Math.max( q[ ioffq + n + off - 1 ],
            q[ ioffq + n - 2 + off ] ) / xx ) * Math.min(
            q[ ioffq + n + off - 1 ], q[ ioffq + n - 2 + off ] );
        } else {
          xx = 0.;
          yy = 0.;
        }
      }
      q[ ioffq + n - 2 + off ] = sigma.getValue() + xx;
      q[ ioffq + n - 1 + off ] = yy + sigma.getValue();
      n = 0;
      continue;
    }
//      var nref = new IntReference( n );
    var offref = new IntReference( off );
    LaPack2.dlasq3( n, q, e, qq, ee, sup, sigma, kend, offref,
      iphase, iconv, eps, tol2, small2, ioffq + off1 - 1,
      ioffe + off1 - 1, ioffqq + off1 - 1, ioffee + off1 - 1 );
//      n = nref.getValue();
    off = offref.getValue();
    if ( sup.getValue() < 0. ) {
      info.setValue( n + off );
      return;
    }
    off1 = off + 1;
  }
}
*/
//*************************************************************************
LaPack3.dlatdf = function( ijob, n, Z, ldz, rhs, rdsumReference,
rdscalReference, ipiv, jpiv, ioffz, ioffipiv, ioffjpiv ) {
  throw new Error("not tested");
  var maxdim = 8;
  var iwork = new Array( maxdim );
  var work = new Array( 4 * maxdim );
  var xm = new Array( maxdim );
  var xp = new Array( maxdim );
  var rdscalReference = new NumberReference();
  var rdsumReference = new NumberReference();
  if ( ijob != 2 ) {
    LaPack0.dlaswp( 1, rhs, ldz, 1, n - 1, ipiv, 1,
      ioffrhs, ioffipiv );
    var pmone = -1.;
    for ( var j = 1; j <= n - 1; j ++ ) {
      var bp = rhs[ ioffrhs + j - 1 ] + 1.;
      var bm = rhs[ ioffrhs + j - 1 ] - 1.;
      var splus = 1.;
      splus += Blas1.ddot( n - j, Z, 1, Z, 1,
        ioffz + j + ( j - 1 ) * ldz, ioffz + j + ( j - 1 ) * ldz );
      var sminu = Blas1.ddot( n - j, Z, 1, rhs, 1,
        ioffz + j + ( j - 1 ) * ldz, ioffrhs + j );
      splus *= rhs[ ioffrhs + j - 1 ];
      if ( splus > sminu ) rhs[ ioffrhs + j - 1 ] = bp;
      else if ( sminu > splus ) rhs[ ioffrhs + j - 1 ] = bm;
      else {
        rhs[ ioffrhs + j - 1 ] += pmone;
        phone = 1.;
      }
      var temp = - rhs[ ioffrhs + j - 1 ];
      Blas1.daxpy( n - j, temp, Z, 1, rhs, 1, ioffz + j, ioffrhs + j );
    } // 10
    Blas1.dcopy( n - 1, rhs, 1, xp, 1, ioffrhs, 0 );
    xp[ n - 1 ] = rhs[ ioffrhs + n - 1 ] + 1.;
    rhs[ ioffrhs + n - 1 ] -= 1.;
    splus = 0.;
    sminu = 0.;
    for ( var i = n; i >= 1; i -- ) {
      var temp = 1. / z[ ioffz + i - 1 + ( i - 1 ) * ldz ];
      xp[ i - 1 ] *= temp;
      rhs[ ioffrhs + i - 1 ] *= temp;
      for ( var k = i + 1; k <= n; k ++ ) {
        xp[ i - 1 ] -= xp[ k - 1 ]
          * ( z[ ioffz + i - 1 + ( k - 1 ) * ldz ] * temp );
        rhs[ ioffrhs + i - 1 ] -= rhs[ ioffrhs + k - 1 ]
          * ( z[ ioffz + i - 1 + ( k - 1 ) * ldz ] * temp );
      } // 20
      splus += Math.abs( xp[ i - 1 ] );
      sminu += Math.abs( rhs[ ioffrhs + i - 1 ] );
    } // 30
    if ( splus > sminu ) {
      Blas1.dcopy( n, xp, 1, rhs, 1, 0, ioffrhs );
    }
    LaPack0.dlaswp( 1, rhs, ldz, 1, n - 1, jpiv, -1,
      ioffrhs, ioffjpiv );
    LaPack0.dlassq( n, rhs, 1, rdscal, rdsum, ioffrhs );
  } else {
    LaPack2.dgecon( 'I', n, Z, ldz, 1., temp, work, iwork, info,
      ioffz, 0, 0 );
    Blas1.dcopy( n, work, 1, xm, 1, n, 0 );
    LaPack0.dlaswp( 1, xm, ldz, 1, n - 1, ipiv, -1, 0, ioffipiv );
    temp = 1. / Math.sqrt( Blas1.ddot( n, xm, 1, xm, 1, 0, 0 ) );
    Blas1.dscal( n, temp, xm, 1, 0 );
    Blas1.dcopy( n, xm, 1, xp, 1, 0, 0 );
    Blas1.daxpy( n, 1., rhs, 1, xp, 1, ioffrhs, 0 );
    Blas1.daxpy( n, -1., xm, 1, rhs, 1, 0, ioffrhs );
    var temprefReference = new NumberReference( temp );
    LaPack1.dgesc2( n, Z, ldz, rhs, ipiv, jpiv, tempref,
      ioffz, ioffrhs, ioffipiv, ioffjpiv );
    temp = tempref.getValue();
    LaPack1.dgesc2( n, Z, ldz, xp, ipiv, jpiv, tempref,
      ioffz, 0, ioffipiv, ioffjpiv );
    temp = tempref.getValue();
    if ( Blas1.dasum( n, xp, 1, 0 ) >
    Blas1.dasum( n, rhs, 1, ioffrhs ) ) {
      Blas1.dcopy( n, xp, 1, rhs, 1, 0, ioffrhs );
    }
    LaPack0.dlassq( n, rhs, 1, rdscal, rdsum, ioffrhs );
  }
}
//*************************************************************************
LaPack3.dopgtr = function( uplo, n, ap, tau, Q, ldq, work, info,
ioffap, iofftau, ioffq, ioffwork ) {
  throw new Error("not debugged: packed matrix");
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( ldq < Math.max( 1, n ) ) info.setValue( -6 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dopgtr', - info.getValue() );
    return;
  }
  if ( n == 0 ) return;
  if ( upper ) {
    var ij = 2;
    for ( var j = 1; j <= n - 1; j ++ ) {
      for ( var i = 1; i <= j - 1; i ++ ) {
        Q[ i - 1 + ( j - 1 ) * ldq ] = ap[ ioffap + ij ];
        ij ++;
      }
      ij += 2;
      Q[ n - 1 + ( j - 1 ) * ldq ] = 0.;
    }
    for ( i = 1; i <= n - 1; i ++ ) {
      Q[ i - 1 + ( n - 1 ) * ldq ] = 0.;
    }
    Q[ n - 1 + ( n - 1 ) * ldq ] = 1.;
    var iinfo = new IntReference();
    LaPack2.dorg2l( n - 1, n - 1, n - 1, Q, ldq, tau, work, iinfo,
      0, iofftau, ioffwork );
  } else {
    Q[ 0 ] = 1.;
    for ( i = 2; i <= n; i ++ ) Q[ i - 1 ] = 0.;
    ij = 3;
    for ( j = 2; j <= n; j ++ ) {
      Q[ ( j - 1 ) * ldq ] = 0.;
      for ( i = j + 1; i <= n; i ++ ) {
        Q[ i - 1 + ( j - 1 ) * ldq ] = ap[ ioffap + ij ];
        ij ++;
      }
      ij += 2;
    }
    if ( n > 1 ) {
      LaPack2.dorg2l( n - 1, n - 1, n - 1, Q, ldq, tau, work, iinfo,
        1 + ldq, iofftau, ioffwork );
    }
  }
}
//*************************************************************************
LaPack3.dorglq = function( m, n, k, A, lda, tau, work, lwork,
info, ioffa, iofftau, ioffwork ) {
  throw new Error("not tested: complicated input");
  info.setValue( 0 );
  var nb = LaPack0.ilaenv( 1, 'dorglq', ' ', m, n, k, -1 );
  var lwkopt = Math.max( 1, m ) * nb;
  work[ ioffwork ] = lwkopt;
  var lquery = ( lwork == -1 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < m ) info.setValue( -2 );
  else if ( k < 0 || k > m ) info.setValue( -3 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -5 );
  else if ( lwork < Math.max( 1, m ) && ! lquery ) info.setValue( -8 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dorglq', - info.getValue() );
    return;
  } else if ( lquery ) return;
  if ( m <= 0 ) {
    work[ ioffwork ] = 1;
    return;
  }
  var nbmin = 2;
  var nx = 0;
  var iws = m;
  if ( nb > 1 && nb < k ) {
    nx = Math.max( 0,
      LaPack0.ilaenv( 3, 'dorglq', ' ', m, n, k, -1 ) );
    if ( nx < k ) {
      var ldwork = m;
      iws = ldwork * nb;
      if ( lwork < iws ) {
        nb = lwork / ldwork;
        nbmin = Math.max( 2,
          LaPack0.ilaenv( 2, 'dorglq', ' ', m, n, k, -1 ) );
      }
    }
  }
  if ( nb >= nbmin && nb < k && nx < k ) {
    var ki = ( ( k - nx - 1 ) / nb ) * nb;
    var kk = Math.min( k, ki + nb );
    for ( var j = 1; j <= kk; j ++ ) {
      for ( var i = kk + 1; i <= m; i ++ ) {
        A[ ioffa + i - 1 + ( j - 1 ) * lda ] = 0.;
      }
    }
  } else kk = 0;
  var iinfo = new IntReference();
  if ( kk < m ) {
    LaPack2.dorgl2( m - kk, n - kk, k - kk, A, lda, tau, work, iinfo,
      ioffa + kk + kk * lda, iofftau + kk, ioffwork );
  }
  if ( kk > 0 ) {
    for ( i = ki + 1; i >= 1; i -= nb ) {
      var ib = Math.min( nb, k - i + 1 );
      if ( i + ib <= m ) {
        LaPack0.dlarft( 'Forward', 'Rowwise', n - i + 1, ib, A, lda,
          tau, work, ldwork, ioffa + i - 1 + ( i - 1 ) * lda,
          iofftau + i - 1, ioffwork );
        LaPack1.dlarfb( 'Right', 'Transpose', 'Forward', 'Rowwise',
          m - i - ib + 1, n - i + 1, ib, A, lda, work, ldwork, A, lda,
          work, ldwork, ioffa + i - 1 + ( i - 1 ) * lda, ioffwork,
          ioffa + i + ib - 1 + ( i - 1 ) * lda, ioffwork + ib );
      }
      LaPack2.dorgl2( ib, n - i + 1, ib, A, lda, tau, work, iinfo,
        ioffa + i - 1 + ( i - 1 ) * lda, iofftau + i - 1, ioffwork );
      for( j = 1; j <= i - 1; j ++ ) {
        for ( var l = i; l <= i + ib - 1; l ++ ) {
          A[ ioffa + l - 1 + ( j - 1 ) * lda ] = 0.;
        }
      }
    }
  }
  work[ ioffwork ] = iws;
}
//*************************************************************************
LaPack3.dorgql = function( m, n, k, A, lda, tau, work, lwork,
info, ioffa, iofftau, ioffwork ) {
  throw new Error("not tested: complicated input");
  info.setValue( 0 );
  var lquery = ( lwork == -1 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 || n > m ) info.setValue( -2 );
  else if ( k < 0 || k > n ) info.setValue( -3 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -5 );
  if ( info.getValue() == 0 ) {
    if ( n == 0 ) var lwkopt = 1;
    else {
      var nb = LaPack0.ilaenv( 1, 'dorgql', ' ', m, n, k, -1 );
      lwkopt = n * nb;
    }
    work[ ioffwork ] = lwkopt;
    if ( lwork < Math.max( 1, n ) && ! lquery ) {
      info.setValue( -8 );
    }
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dorgql', - info.getValue() );
    return;
  } else if ( lquery ) return;
  if ( n <= 0 ) return;
  var nbmin = 2;
  var nx = 0;
  var iws = n;
  if ( nb > 1 && nb < k ) {
    nx = Math.max( 0,
      LaPack0.ilaenv( 3, 'dorgql', ' ', m, n, k, -1 ) );
    if ( nx < k ) {
      var ldwork = n;
      iws = ldwork * nb;
      if ( lwork < iws ) {
        nb = lwork / ldwork;
        nbmin = Math.max( 2,
          LaPack0.ilaenv( 2, 'dorgql', ' ', m, n, k, -1 ) );
      }
    }
  }
  if ( nb >= nbmin && nb < k && nx < k ) {
    var kk = Math.min( k, ( ( k - nx + nb - 1 ) / nb ) * nb );
    for ( var j = 1; j <= n - kk; j ++ ) {
      for ( var i = m - kk + 1; i <= m; i ++ ) {
        A[ ioffa + i - 1 + ( j - 1 ) * lda ] = 0.;
      }
    }
  } else kk = 0;
  var iinfo = new IntReference();
  LaPack2.dorg2l( m - kk, n - kk, k - kk, A, lda, tau, work, iinfo,
    ioffa, iofftau, ioffwork );
  if ( kk > 0 ) {
    for ( i = k - kk + 1; i <= k; i += nb ) {
      var ib = Math.min( nb, k - i + 1 );
      if ( n - k + i > 1 ) {
        LaPack0.dlarft( 'Backward', 'Columnwise', m - k + i + ib - 1,
          ib, A, lda, tau, work, ldwork,
          ioffa + ( n - k + i - 1 ) * lda,
          iofftau + i - 1, ioffwork );
        LaPack1.dlarfb( 'Left', 'No transpose', 'Backward',
          'Columnwise', m - k + i + ib - 1, n - k + i - 1, ib, A, lda,
          work, ldwork, A, lda, work, ldwork,
          ioffa + ( n - k + i - 1 ) * lda,
          ioffwork, ioffa, ioffwork + ib );
      }
      LaPack2.dorgl2( m - k + i + ib - 1, ib, ib, A, lda, tau, work,
        iinfo, ioffa + ( n - k + i - 1 ) * lda, iofftau + i - 1,
        ioffwork );
      for( j = n - k + i; j <= n - k + i + ib - 1; j ++ ) {
        for ( var l = m - k + i + ib; l <= m; l ++ ) {
          A[ ioffa + l - 1 + ( j - 1 ) * lda ] = 0.;
        }
      }
    }
  }
  work[ ioffwork ] = iws;
}
//*************************************************************************
LaPack3.dorgqr = function( m, n, k, A, lda, tau, work, lwork,
info, ioffa, iofftau, ioffwork ) {
  throw new Error("not tested: complicated input");
  info.setValue( 0 );
  var nb = LaPack0.ilaenv( 1, 'dorgqr', ' ', m, n, k, -1 );
  var lwkopt = Math.max( 1, n ) * nb;
  work[ ioffwork ] = lkwopt;
  var lquery = ( lwork == -1 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < 0 || n > m ) info.setValue( -2 );
  else if ( k < 0 || k > n ) info.setValue( -3 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -5 );
  else if ( lwork < Math.max( 1, n )  && ! lquery ) info.setValue( -8 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dorgqr', - info.getValue() );
    return;
  } else if ( lquery ) return;
  if ( n <= 0 ) {
    work[ ioffwork ] = 1;
    return;
  }
  var nbmin = 2;
  var nx = 0;
  var iws = n;
  if ( nb > 1 && nb < k ) {
    nx = Math.max( 0,
      LaPack0.ilaenv( 3, 'dorgqr', ' ', m, n, k, -1 ) );
    if ( nx < k ) {
      var ldwork = n;
      iws = ldwork * nb;
      if ( lwork < iws ) {
        nb = lwork / ldwork;
        nbmin = Math.max( 2,
          LaPack0.ilaenv( 2, 'dorgqr', ' ', m, n, k, -1 ) );
      }
    }
  }
  if ( nb >= nbmin && nb < k && nx < k ) {
    var ki = ( ( k - nx - 1 ) / nb ) * nb;
    var kk = Math.min( k, ki + nb );
    for ( var j = kk + 1; j <= n; j ++ ) {
      for ( var i = 1; i <= kk; i ++ ) {
        A[ ioffa + i - 1 + ( j - 1 ) * lda ] = 0.;
      }
    }
  } else kk = 0;
  var iinfo = new IntReference();
  if ( kk < n ) {
    LaPack2.dorg2r( m - kk, n - kk, k - kk, A, lda, tau, work, iinfo,
      ioffa + kk + kk * lda, iofftau + kk, ioffwork );
  }
  if ( kk > 0 ) {
    for ( i = ki + 1; i >= 1; i -= nb ) {
      var ib = Math.min( nb, k - i + 1 );
      if ( i + ib <= n ) {
        LaPack0.dlarft( 'Forward', 'Columnwise', m - i + 1, ib, A, lda,
          tau, work, ldwork, ioffa + i - 1 + ( i - 1 ) * lda,
          iofftau + i - 1, ioffwork );
        LaPack1.dlarfb( 'Left', 'No transpose', 'Forward',
          'Columnwise', m - i + 1, n - i - ib + 1, ib, A, lda,
          work, ldwork, A, lda, work, ldwork,
          ioffa + i - 1 + ( i - 1 ) * lda, ioffwork,
          ioffa + i - 1 + ( i + ib - 1 ) * lda, ioffwork + ib );
      }
      LaPack2.dorg2r( m - i + 1, ib, ib, A, lda, tau, work, iinfo,
        ioffa + i - 1 + ( i - 1 ) * lda, iofftau + i - 1, ioffwork );
      for( j = i; j <= i + ib - 1; j ++ ) {
        for ( var l = 1; l <= i - 1; l ++ ) {
          A[ ioffa + l - 1 + ( j - 1 ) * lda ] = 0.;
        }
      }
    }
  }
  work[ ioffwork ] = iws;
}
//*************************************************************************
LaPack3.dormqr = function( side, trans, m, n, k, A, lda, tau, C,
ldc, work, lwork, info, ioffa, iofftau, ioffc, ioffwork ) {
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
  else if ( lda < Math.max( 1, nq ) ) info.setValue( -7 );
  else if ( ldc < Math.max( 1, m ) ) info.setValue( -10 );
  else if ( lwork < Math.max( 1, nw ) && ! lquery ) info.setValue( -12 );
  if ( info.getValue() == 0 ) {
    var nb = Math.min( nbmax,
      LaPack0.ilaenv( 1, 'dormqr', side + trans, m, n, k, -1 ) );
    var lwkopt = Math.max( 1, nw ) * nb;
    work[ ioffwork ] = lwkopt;
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dormqr', -info.getValue() );
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
        LaPack0.ilaenv( 2, 'dormqr', side + trans, m, n, k, -1 ) );
    }
  } else iws = nw;
  if ( nb < nbmin || nb >= k ) {
    var iinfo = new IntReference;
    LaPack2.dorm2r( side, trans, m, n, k, A, lda, tau, C, ldc, work,
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
    if ( left ) {
      var ni = n;
      var jc = 1;
    } else {
      var mi = m;
      var ic = 1;
    }
    for ( var i = i1; ( i3 == nb && i <= i2 ) ||
    ( i3 == -nb && i >= i2 ); i += i3 ) {
      var ib = Math.min( nb, k - i + 1 );
      LaPack0.dlarft( 'Forward', 'Columnwise', nq - i + 1, ib, A, lda,
        tau, T, ldt, ioffa + i - 1 + ( i - 1 ) * lda, iofftau + i - 1,
        0 );
      if ( left ) {
        mi = m - i + 1;
        ic = i;
      } else {
        ni = n - i + 1;
        jc = i;
      }
      LaPack1.dlarfb( side, trans, 'Forward', 'Columnwise', mi, ni,
        ib, A, lda, T, ldt, C, ldc, work, ldwork,
        ioffa + i - 1 + ( i - 1 ) * lda, 0,
        ioffc + ic - 1 + ( jc - 1 ) * ldc, ioffwork );
    }
  }
  work[ ioffwork ] = lwkopt;
}
//*************************************************************************
LaPack3.dsytrd = function( uplo, n, A, lda, d, e, tau, work,
lwork, info, ioffa, ioffd, ioffe, iofftau, ioffwork ) {
  throw new Error("not tested");
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var lquery = ( lwork == -1 );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -4 );
  else if ( lwork < 1  && ! lquery ) info.setValue( -9 );
  if ( info.getValue() == 0 ) {
    var nb = LaPack0.ilaenv( 1, 'dsytrd', uplo, n, -1, -1, -1 );
    var lwkopt = n * nb;
    work[ ioffwork ] = lwkopt;
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dsytrd', - info.getValue() );
    return;
  } else if ( lquery ) return;
  if ( n == 0 ) {
    work[ ioffwork ] = 1;
    return;
  }
  var nx = n;
  var iws = 1;
  if ( nb > 1 && nb < n ) {
    nx = Math.max( nb,
      LaPack0.ilaenv( 3, 'dsytrd', uplo, n, -1, -1, -1 ) );
    if ( nx < n ) {
      var ldwork = n;
      iws = ldwork * nb;
      if ( lwork < iws ) {
        nb = Math.max( lwork / ldwork, 1 );
        var nbmin =
          LaPack0.ilaenv( 2, 'dsytrd', uplo, n, -1, -1, -1 );
        if ( nb < nbmin ) nx = n;
      }
    } else nx = n;
  } else nb = 1;
  var iinfo = new IntReference();
  if ( upper ) {
    var kk = n - ( ( n - nx + nb - 1 ) / nb ) * nb;
    for ( var i = n - nb + 1; i >= kk + 1; i -= nb ) {
      LaPack2.dlatrd( uplo, i + nb - 1, nb, A, lda, e, tau, work,
        ldwork, ioffa, ioffe, iofftau, ioffwork );
      Blas3.dsyr2k( uplo, 'No transpose', i - 1, nb, -1., A, lda,
        work, ldwork, 1., A, lda, ioffa + ( i - 1 ) * lda, ioffwork,
        ioffa );
      for ( var j = i; j <= i + nb - 1; j ++ ) {
        A[ ioffa + j - 2 + ( j - 1 ) * lda ] = e[ ioffe + j - 2 ];
        d[ ioffd + j - 1 ] = A[ ioffa + j - 1 + ( j - 1 ) * lda ];
      }
    }
    LaPack2.dsytd2( uplo, kk, A, lda, d, e, tau, iinfo,
      ioffa, ioffd, ioffe, iofftau );
  } else {
    for ( i = 1; i <= n - nx; i += nb ) {
      LaPack2.dlatrd( uplo, n - i + 1, nb, A, lda, e, tau, work,
        ldwork, ioffa + i - 1 + ( i - 1 ) * lda, ioffe + i - 1,
        iofftau + i - 1, ioffwork );
      Blas3.dsyr2k( uplo, 'No transpose', n - i - nb + 1, nb, -1.,
        A, lda, work, ldwork, 1., A, lda,
        ioffa + i + nb - 1 + ( i - 1 ) * lda,
        ioffwork + nb, ioffa + i + nb - 1 + ( i + nb - 1 ) * lda );
      for ( j = i; j <= i + nb - 1; j ++ ) {
        A[ ioffa + j + ( j -1 ) * lda ] = e[ ioffe + j - 1 ];
        d[ ioffd + j - 1 ] = A[ ioffa + j - 1 + ( j - 1 ) * lda ];
      }
    }
    LaPack2.dsytd2( uplo, n - i + 1, A, lda, d, e, tau, iinfo,
      ioffa + i - 1 + ( i - 1 ) * lda, ioffd + i - 1, ioffe + i - 1,
      iofftau + i - 1 );
  }
  work[ ioffwork ] = iws;
}
//*************************************************************************
LaPack3.dsytri2 = function( uplo, n, A, lda, ipiv, work, lwork,
info, ioffa, ioffipiv, ioffwork ) {
  throw new Error("not tested");
  info.setValue( 0 );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var lquery = ( lwork == -1 );
  var nbmax = LaPack0.ilaenv( 1, 'dsytrf', uplo, n, -1, -1, -1 );
  var minsize =
    ( nbmax >= n ? n : ( n + nbmax + 1 ) * ( nbmax + 3 ) );
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info.setValue( -1 );
  } else if ( n < 0 ) info.setValue( -2 );
  else if ( lda < Math.max( 1, n ) ) info.setValue( -4 );
  else if ( lwork < minsize && ! lquery ) info.setValue( -7 );
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dsytri2', - info.getValue() );
    return;
  } else if ( lquery ) {
    work[ ioffwork ] = minsize;
    return;
  }
  if ( n == 0 ) return;
  if ( nbmax >= n ) {
    LaPack0.dsytri( uplo, n, A, lda, ipiv, work, info,
      ioffa, ioffipiv, ioffwork );
  } else {
    LaPack2.dsytri2x( uplo, n, A, lda, ipiv, work, nbmax, info,
      ioffa, ioffipiv, ioffwork );
  }
}
//*************************************************************************
LaPack3.dtgsja = function( jobu, jobv, jobq, m, p, n, k, l, A,
lda, B, ldb, tola, tolb, alphaReference, betaReference, U, ldu, V, ldv, Q,
ldq, work, ncycle, info ) {
  throw new Error("not programmed: generalized eigenvalue problem");
}
//*************************************************************************
LaPack3.dtpqrt = function( m, n, l, nb, A, lda, B, ldb, T, ldt,
work, info, ioffa, ioffb, iofft, ioffwork) {
  throw new Error("not programmed: WY storage");
}
//*************************************************************************
LaPack3.dtzrzf = function( m, n, A, lda, tau, work, lwork, info,
ioffa, iofftau, ioffwork ) {
  info.setValue( 0 );
  var lquery = ( lwork == -1 );
  if ( m < 0 ) info.setValue( -1 );
  else if ( n < m ) info.setValue( -2 );
  else if ( lda < Math.max( 1, m ) ) info.setValue( -4 );
  if ( info.getValue() == 0 ) {
    if ( m == 0 || m == n ) {
      var lwkopt = 1;
      var lwkmin = 1;
    } else {
      var nb = LaPack0.ilaenv( 1, 'dgerqf', ' ', m, n, -1, -1 );
      lwkopt = m * nb;
      lwkmin = Math.max( 1, m );
    }
    work[ ioffwork ] = lwkopt;
    if ( lwork < lwkmin && ! lquery ) info.setValue( -7 );
  }
  if ( info.getValue() != 0 ) {
    Blas2.xerbla( 'dtzrzf', - info.getValue() );
    return;
  } else if ( lquery ) return;
  if ( m == 0 ) return;
  else if ( m == n ) {
    for ( var i = 1; i <= n; i ++ ) tau[ iofftau + i - 1 ] = 0.;
    return;
  }
  var nbmin = 2;
  var nx = 1;
  var iws = m;
  if ( nb > 1 && nb < m ) {
    nx = Math.max( 0,
      LaPack0.ilaenv( 3, 'dgerqf', ' ', m, n, -1, -1 ) );
    if ( nx < m ) {
      var ldwork = m;
      iws = ldwork * nb;
      if ( lwork < iws ) {
        nb = lwork / ldwork;
        nbmin = Math.max( 2,
          LaPack0.ilaenv( 2, 'dgerqf', ' ', m, n, -1, -1 ) );
      }
    }
  }
  if ( nb >= nbmin && nb < m && nx < m ) {
    var m1 = Math.min( m + 1, n );
    var ki = ( ( m - nx - 1 ) / nb ) * nb;
    var kk = Math.min( m, ki + nb );
    for ( i = m - kk + ki + 1; i >= m - kk + 1;  i -= nb ) {
      var ib = Math.min( m - i + 1, nb );
      LaPack2.dlatrz( ib, n - i + 1, n - m, A, lda, tau, work,
        ioffa + i - 1 + ( i - 1 ) * lda, iofftau + i - 1, ioffwork );
      if ( i > 1 ) {
        LaPack0.dlarzt( 'Backward', 'Rowwise', n - m, ib, A, lda,
          tau, work, ldwork, ioffa + i - 1 + ( m1 - 1 ) * lda,
          iofftau + i - 1, ioffwork );
        LaPack0.dlarzb( 'Right', 'No transpose', 'Backward', 'Rowwise',
          i - 1, n - i + 1, ib, n - m, A, lda, work, ldwork, A, lda,
          work, ldwork, ioffa + i - 1 + ( m1 - 1 ) * lda, ioffwork,
          ioffa + ( i - 1 ) * lda, ioffwork + ib );
      }
    }
    var mu = i + nb - 1;
  } else mu = m;
  if ( mu > 0 ) {
    LaPack2.dlatrz( mu, n, n - m, A, lda, tau, work,
      ioffa, iofftau, ioffwork );
  }
  work[ ioffwork ] = lwkopt;
}
function testLaPack3() {
//test_dgehrd();
//test_dgelqf();
//test_dgeqlf();
//test_dgeqrf();
//test_dgerqf();
//test_dgesvj();
  test_dlaexc();
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgehrd() {
  document.getElementById("debug_textarea").value +=
    "testing dgehrd *************" + "\n";
  var n = 9;
  var lda = 10;
  var lwork = 11;
  var ioffa = 1;
  var iofftau = 2;
  var ioffwork = 3;
  var A = new Array( ioffa + lda * n );
  var tau = new Array( iofftau + n - 1 );
  var work = new Array( ioffwork + lwork );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + 2 * j );
    }
  }
  var ilo = 3;
  var ihi = 6;
  for( j = 0; j < ilo; j ++ ) {
    for ( i = j + 1; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 0.;
    }
  }
  for( j = ihi; j < n; j ++ ) {
    for ( i = j + 1; i < n; i ++ ) {
      A[ ioffa + i + j * lda ] = 0.;
    }
  }
  var info = new IntReference();
  LaPack3.dgehrd(n,ilo,ihi,A,lda,tau,work,lwork,info,
    ioffa,iofftau,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dgehrd: info = " + info.value + "\n";
  document.getElementById("debug_textarea").value +=
    "tau = "
    + tau[ iofftau + 0 ] + " "
    + tau[ iofftau + 1 ] + " "
    + tau[ iofftau + 2 ] + " "
    + tau[ iofftau + 3 ] + " "
    + tau[ iofftau + 4 ] + " "
    + tau[ iofftau + 5 ] + " "
    + tau[ iofftau + 6 ] + " "
    + tau[ iofftau + 7 ]  + "\n";
  for ( j = 0; j < n; j ++ ) {
    document.getElementById("debug_textarea").value +=
      "A[*," + j + "] = "
      + A[ ioffa + 0 + j * lda ] + " "
      + A[ ioffa + 1 + j * lda ] + " "
      + A[ ioffa + 2 + j * lda ] + " "
      + A[ ioffa + 3 + j * lda ] + " "
      + A[ ioffa + 4 + j * lda ] + " "
      + A[ ioffa + 5 + j * lda ] + " "
      + A[ ioffa + 6 + j * lda ] + " "
      + A[ ioffa + 7 + j * lda ] + " "
      + A[ ioffa + 8 + j * lda ]  + "\n";
  }
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgelqf() {
  document.getElementById("debug_textarea").value +=
    "testing dgelqf *************" + "\n";
  var m = 4;
  var n = 3;
  var lda = 5;
  var lwork = 4;
  var ioffa = 1;
  var iofftau = 2;
  var ioffwork = 3;
  var A = new Array( ioffa + lda * n );
  var tau = new Array( iofftau + Math.min( m, n ) );
  var work = new Array( ioffwork + lwork );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + 2 * j );
    }
  }
  var info = new IntReference();
  LaPack3.dgelqf(m,n,A,lda,tau,work,lwork,info,ioffa,iofftau,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dgelqf: info = " + info.value + "\n";
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
    + A[ ioffa + 3 + 2 * lda ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgeqlf() {
  document.getElementById("debug_textarea").value +=
    "testing dgeqlf *************" + "\n";
  var m = 4;
  var n = 3;
  var lda = 5;
  var lwork = 4;
  var ioffa = 1;
  var iofftau = 2;
  var ioffwork = 3;
  var A = new Array( ioffa + lda * n );
  var tau = new Array( iofftau + Math.min( m, n ) );
  var work = new Array( ioffwork + lwork );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + 2 * j );
    }
  }
  var info = new IntReference();
  LaPack3.dgeqlf(m,n,A,lda,tau,work,lwork,info,ioffa,iofftau,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dgeqlf: info = " + info.value + "\n";
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
    + A[ ioffa + 3 + 2 * lda ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgeqrf() {
  document.getElementById("debug_textarea").value +=
    "testing dgeqrf *************" + "\n";
  var m = 4;
  var n = 3;
  var lda = 5;
  var lwork = 4;
  var ioffa = 1;
  var iofftau = 2;
  var ioffwork = 3;
  var A = new Array( ioffa + lda * n );
  var tau = new Array( iofftau + Math.min( m, n ) );
  var work = new Array( ioffwork + lwork );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + 2 * j );
    }
  }
  var info = new IntReference();
  LaPack3.dgeqrf(m,n,A,lda,tau,work,lwork,info,ioffa,iofftau,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dgeqrf: info = " + info.value + "\n";
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
    + A[ ioffa + 3 + 2 * lda ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgerqf() {
  document.getElementById("debug_textarea").value +=
    "testing dgerqf *************" + "\n";
  var m = 4;
  var n = 3;
  var lda = 5;
  var lwork = 4;
  var ioffa = 1;
  var iofftau = 2;
  var ioffwork = 3;
  var A = new Array( ioffa + lda * n );
  var tau = new Array( iofftau + Math.min( m, n ) );
  var work = new Array( ioffwork + lwork );
  for ( var j = 0; j < n; j ++ ) {
    for ( var i = 0; i < m; i ++ ) {
      A[ ioffa + i + j * lda ] = 1. / Number( 1 + i + 2 * j );
    }
  }
  var info = new IntReference();
  LaPack3.dgerqf(m,n,A,lda,tau,work,lwork,info,ioffa,iofftau,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dgerqf: info = " + info.value + "\n";
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
    + A[ ioffa + 3 + 2 * lda ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgesvj() {
  document.getElementById("debug_textarea").value +=
    "testing dgesvj *************" + "\n";
  var ioffl = 1;
  var ioffv = 2;
  var ioffsva = 3;
  var ioffwork = 4;
  var ldl = 4;
  var ldv = 4;
  var L = new Array( ioffl + ldl * 3 );
  var V = new Array( ioffv + ldv * 4 );
  var sva = new Array( ioffsva + 4 );
  var work = new Array( ioffwork + 8 );
  var info = new IntReference();
  L[ ioffl + 0 + 0 * ldl ] = 1.;
  L[ ioffl + 1 + 0 * ldl ] = 1.;
  L[ ioffl + 2 + 0 * ldl ] = 1.;
  L[ ioffl + 3 + 0 * ldl ] = 1.;
  L[ ioffl + 0 + 1 * ldl ] = 0.;
  L[ ioffl + 1 + 1 * ldl ] = 0.;
  L[ ioffl + 2 + 1 * ldl ] = 1.;
  L[ ioffl + 3 + 1 * ldl ] = 1.;
  L[ ioffl + 0 + 2 * ldl ] = 0.;
  L[ ioffl + 1 + 2 * ldl ] = 0.;
  L[ ioffl + 2 + 2 * ldl ] = 1.;
  L[ ioffl + 3 + 2 * ldl ] = 1.;
  LaPack3.dgesvj('L','U','V',4,3,L,4,sva,4,V,4,work,8,info,
    ioffl,ioffsva,ioffv,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dgesvj('L',...): A = "
    + L[ ioffl + 0 + 0 * ldl ] + " "
    + L[ ioffl + 1 + 0 * ldl ] + " "
    + L[ ioffl + 2 + 0 * ldl ] + " "
    + L[ ioffl + 3 + 0 * ldl ] + " "
    + L[ ioffl + 0 + 1 * ldl ] + " "
    + L[ ioffl + 1 + 1 * ldl ] + " "
    + L[ ioffl + 2 + 1 * ldl ] + " "
    + L[ ioffl + 3 + 1 * ldl ] + " "
    + L[ ioffl + 0 + 2 * ldl ] + " "
    + L[ ioffl + 1 + 2 * ldl ] + " "
    + L[ ioffl + 2 + 2 * ldl ] + " "
    + L[ ioffl + 3 + 2 * ldl ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "sva = "
    + sva[ ioffsva + 0 ] + " "
    + sva[ ioffsva + 1 ] + " "
    + sva[ ioffsva + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "V = "
    + V[ ioffv + 0 + 0 * ldv ] + " "
    + V[ ioffv + 1 + 0 * ldv ] + " "
    + V[ ioffv + 2 + 0 * ldv ] + " "
    + V[ ioffv + 0 + 1 * ldv ] + " "
    + V[ ioffv + 1 + 1 * ldv ] + " "
    + V[ ioffv + 2 + 1 * ldv ] + " "
    + V[ ioffv + 0 + 2 * ldv ] + " "
    + V[ ioffv + 1 + 2 * ldv ] + " "
    + V[ ioffv + 2 + 2 * ldv ]  + "\n";

  var ioffu = 5;
  var ldu = 3;
  var U = new Array( ioffu + ldu * 3 );
  U[ ioffu + 0 + 0 * ldu ] = 1.;
  U[ ioffu + 1 + 0 * ldu ] = 0.;
  U[ ioffu + 2 + 0 * ldu ] = 0.;
  U[ ioffu + 0 + 1 * ldu ] = 1.;
  U[ ioffu + 1 + 1 * ldu ] = 1.;
  U[ ioffu + 2 + 1 * ldu ] = 0.;
  U[ ioffu + 0 + 2 * ldu ] = 1.;
  U[ ioffu + 1 + 2 * ldu ] = 1.;
  U[ ioffu + 2 + 2 * ldu ] = 1.;
  LaPack3.dgesvj('U','U','V',3,3,U,3,sva,3,V,3,work,8,info,
    ioffu,ioffsva,ioffv,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dgesvj('L',...): A = "
    + U[ ioffu + 0 + 0 * ldu ] + " "
    + U[ ioffu + 1 + 0 * ldu ] + " "
    + U[ ioffu + 2 + 0 * ldu ] + " "
    + U[ ioffu + 0 + 1 * ldu ] + " "
    + U[ ioffu + 1 + 1 * ldu ] + " "
    + U[ ioffu + 2 + 1 * ldu ] + " "
    + U[ ioffu + 0 + 2 * ldu ] + " "
    + U[ ioffu + 1 + 2 * ldu ] + " "
    + U[ ioffu + 2 + 2 * ldu ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "sva = "
    + sva[ ioffsva + 0 ] + " "
    + sva[ ioffsva + 1 ] + " "
    + sva[ ioffsva + 2 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "V = "
    + V[ ioffv + 0 + 0 * ldv ] + " "
    + V[ ioffv + 1 + 0 * ldv ] + " "
    + V[ ioffv + 2 + 0 * ldv ] + " "
    + V[ ioffv + 0 + 1 * ldv ] + " "
    + V[ ioffv + 1 + 1 * ldv ] + " "
    + V[ ioffv + 2 + 1 * ldv ] + " "
    + V[ ioffv + 0 + 2 * ldv ] + " "
    + V[ ioffv + 1 + 2 * ldv ] + " "
    + V[ ioffv + 2 + 2 * ldv ]  + "\n";

  var ioffa = 6;
  var lda = 4;
  var A = new Array( ioffa + lda * 4 );
  for ( var j = 0; j < 4; j ++ ) {
    for ( var i = 0; i < 4; i ++ ) {
      A[ ioffa + i + j * lda ] = i + j;
    }
  }
  LaPack3.dgesvj('G','U','V',4,4,A,4,sva,4,V,4,work,8,info,
    ioffa,ioffsva,ioffv,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dgesvj('L',...): A = "
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
    "sva = "
    + sva[ ioffsva + 0 ] + " "
    + sva[ ioffsva + 1 ] + " "
    + sva[ ioffsva + 2 ] + " "
    + sva[ ioffsva + 3 ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "V = "
    + V[ ioffv + 0 + 0 * ldv ] + " "
    + V[ ioffv + 1 + 0 * ldv ] + " "
    + V[ ioffv + 2 + 0 * ldv ] + " "
    + V[ ioffv + 3 + 0 * ldv ] + " "
    + V[ ioffv + 0 + 1 * ldv ] + " "
    + V[ ioffv + 1 + 1 * ldv ] + " "
    + V[ ioffv + 2 + 1 * ldv ] + " "
    + V[ ioffv + 3 + 1 * ldv ] + " "
    + V[ ioffv + 0 + 2 * ldv ] + " "
    + V[ ioffv + 1 + 2 * ldv ] + " "
    + V[ ioffv + 2 + 2 * ldv ] + " "
    + V[ ioffv + 3 + 2 * ldv ] + " "
    + V[ ioffv + 0 + 3 * ldv ] + " "
    + V[ ioffv + 1 + 3 * ldv ] + " "
    + V[ ioffv + 2 + 3 * ldv ] + " "
    + V[ ioffv + 3 + 3 * ldv ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dlaexc() {
  document.getElementById("debug_textarea").value +=
    "testing dlaexc *************" + "\n";
  var ldt = 5;
  var ldq = 6;
  var ioffq = 1;
  var iofft = 2;
  var ioffwork = 3;
  var n = 4;
  var Q = new Array( ioffq + ldq * n );
  var T = new Array( iofft + ldt * n );
  var work = new Array( ioffwork + n );
  var info = new IntReference();
  var i = -1;
  var j = -1;

  n=2;
  T[ iofft + 0 + 0 * ldt ] = 1.;
  T[ iofft + 1 + 0 * ldt ] = 0.;
  T[ iofft + 0 + 1 * ldt ] = 3.;
  T[ iofft + 1 + 1 * ldt ] = 2.;
  LaPack3.dlaexc(false,n,T,ldt,Q,ldq,1,1,1,work,info,
    iofft,ioffq,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaexc(F,...,1,1,1,...): info = " + info.value + "\n";
  document.getElementById("debug_textarea").value +=
    "T = "
    + T[ iofft + 0 + 0 * ldt ] + " "
    + T[ iofft + 1 + 0 * ldt ] + " "
    + T[ iofft + 0 + 1 * ldt ] + " "
    + T[ iofft + 1 + 1 * ldt ]  + "\n";
  T[ iofft + 0 + 0 * ldt ] = 1.;
  T[ iofft + 1 + 0 * ldt ] = 0.;
  T[ iofft + 0 + 1 * ldt ] = 3.;
  T[ iofft + 1 + 1 * ldt ] = 2.;
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      Q[ ioffq + i + j * ldq ] = 0.;
    }
    Q[ ioffq + j + j * ldq ] = 1.;
  }
  LaPack3.dlaexc(true,n,T,ldt,Q,ldq,1,1,1,work,info,
    iofft,ioffq,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaexc(T,...,1,1,1,...): info = " + info.value + "\n";
  document.getElementById("debug_textarea").value +=
    "T = "
    + T[ iofft + 0 + 0 * ldt ] + " "
    + T[ iofft + 1 + 0 * ldt ] + " "
    + T[ iofft + 0 + 1 * ldt ] + " "
    + T[ iofft + 1 + 1 * ldt ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "Q = "
    + Q[ ioffq + 0 + 0 * ldq ] + " "
    + Q[ ioffq + 1 + 0 * ldq ] + " "
    + Q[ ioffq + 0 + 1 * ldq ] + " "
    + Q[ ioffq + 1 + 1 * ldq ]  + "\n";

  n=3;
  T[ iofft + 0 + 0 * ldt ] = 1.;
  T[ iofft + 1 + 0 * ldt ] = 0.;
  T[ iofft + 2 + 0 * ldt ] = 0.;
  T[ iofft + 0 + 1 * ldt ] = 4.;
  T[ iofft + 1 + 1 * ldt ] = 2.;
  T[ iofft + 2 + 1 * ldt ] = 3.;
  T[ iofft + 0 + 2 * ldt ] = 5.;
  T[ iofft + 1 + 2 * ldt ] = -3.;
  T[ iofft + 2 + 2 * ldt ] = 2.;
  LaPack3.dlaexc(false,n,T,ldt,Q,ldq,1,1,2,work,info,
    iofft,ioffq,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaexc(F,...,1,1,2,...): info = " + info.value + "\n";
  document.getElementById("debug_textarea").value +=
    "T = "
    + T[ iofft + 0 + 0 * ldt ] + " "
    + T[ iofft + 1 + 0 * ldt ] + " "
    + T[ iofft + 2 + 0 * ldt ] + " "
    + T[ iofft + 0 + 1 * ldt ] + " "
    + T[ iofft + 1 + 1 * ldt ] + " "
    + T[ iofft + 2 + 1 * ldt ] + " "
    + T[ iofft + 0 + 2 * ldt ] + " "
    + T[ iofft + 1 + 2 * ldt ] + " "
    + T[ iofft + 2 + 2 * ldt ]  + "\n";
  T[ iofft + 0 + 0 * ldt ] = 1.;
  T[ iofft + 1 + 0 * ldt ] = 0.;
  T[ iofft + 2 + 0 * ldt ] = 0.;
  T[ iofft + 0 + 1 * ldt ] = 4.;
  T[ iofft + 1 + 1 * ldt ] = 2.;
  T[ iofft + 2 + 1 * ldt ] = 3.;
  T[ iofft + 0 + 2 * ldt ] = 5.;
  T[ iofft + 1 + 2 * ldt ] = -3.;
  T[ iofft + 2 + 2 * ldt ] = 2.;
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      Q[ ioffq + i + j * ldq ] = 0.;
    }
    Q[ ioffq + j + j * ldq ] = 1.;
  }
  LaPack3.dlaexc(true,n,T,ldt,Q,ldq,1,1,2,work,info,
    iofft,ioffq,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaexc(T,...,1,1,2,...): info = " + info.value + "\n";
  document.getElementById("debug_textarea").value +=
    "T = "
    + T[ iofft + 0 + 0 * ldt ] + " "
    + T[ iofft + 1 + 0 * ldt ] + " "
    + T[ iofft + 2 + 0 * ldt ] + " "
    + T[ iofft + 0 + 1 * ldt ] + " "
    + T[ iofft + 1 + 1 * ldt ] + " "
    + T[ iofft + 2 + 1 * ldt ] + " "
    + T[ iofft + 0 + 2 * ldt ] + " "
    + T[ iofft + 1 + 2 * ldt ] + " "
    + T[ iofft + 2 + 2 * ldt ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "Q = "
    + Q[ ioffq + 0 + 0 * ldq ] + " "
    + Q[ ioffq + 1 + 0 * ldq ] + " "
    + Q[ ioffq + 2 + 0 * ldq ] + " "
    + Q[ ioffq + 0 + 1 * ldq ] + " "
    + Q[ ioffq + 1 + 1 * ldq ] + " "
    + Q[ ioffq + 2 + 1 * ldq ] + " "
    + Q[ ioffq + 0 + 2 * ldq ] + " "
    + Q[ ioffq + 1 + 2 * ldq ] + " "
    + Q[ ioffq + 2 + 2 * ldq ]  + "\n";

  n=3;
  T[ iofft + 0 + 0 * ldt ] = 2.;
  T[ iofft + 1 + 0 * ldt ] = 3.;
  T[ iofft + 2 + 0 * ldt ] = 0.;
  T[ iofft + 0 + 1 * ldt ] = -3.;
  T[ iofft + 1 + 1 * ldt ] = 2.;
  T[ iofft + 2 + 1 * ldt ] = 0.;
  T[ iofft + 0 + 2 * ldt ] = 4.;
  T[ iofft + 1 + 2 * ldt ] = 5.;
  T[ iofft + 2 + 2 * ldt ] = 1.;
  LaPack3.dlaexc(false,n,T,ldt,Q,ldq,1,2,1,work,info,
    iofft,ioffq,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaexc(F,...,1,2,1,...): info = " + info.value + "\n";
  document.getElementById("debug_textarea").value +=
    "T = "
    + T[ iofft + 0 + 0 * ldt ] + " "
    + T[ iofft + 1 + 0 * ldt ] + " "
    + T[ iofft + 2 + 0 * ldt ] + " "
    + T[ iofft + 0 + 1 * ldt ] + " "
    + T[ iofft + 1 + 1 * ldt ] + " "
    + T[ iofft + 2 + 1 * ldt ] + " "
    + T[ iofft + 0 + 2 * ldt ] + " "
    + T[ iofft + 1 + 2 * ldt ] + " "
    + T[ iofft + 2 + 2 * ldt ]  + "\n";
  T[ iofft + 0 + 0 * ldt ] = 2.;
  T[ iofft + 1 + 0 * ldt ] = 3.;
  T[ iofft + 2 + 0 * ldt ] = 0.;
  T[ iofft + 0 + 1 * ldt ] = -3.;
  T[ iofft + 1 + 1 * ldt ] = 2.;
  T[ iofft + 2 + 1 * ldt ] = 0.;
  T[ iofft + 0 + 2 * ldt ] = 4.;
  T[ iofft + 1 + 2 * ldt ] = 5.;
  T[ iofft + 2 + 2 * ldt ] = 1.;
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      Q[ ioffq + i + j * ldq ] = 0.;
    }
    Q[ ioffq + j + j * ldq ] = 1.;
  }
  LaPack3.dlaexc(true,n,T,ldt,Q,ldq,1,2,1,work,info,
    iofft,ioffq,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaexc(T,...,1,2,1,...): info = " + info.value + "\n";
  document.getElementById("debug_textarea").value +=
    "T = "
    + T[ iofft + 0 + 0 * ldt ] + " "
    + T[ iofft + 1 + 0 * ldt ] + " "
    + T[ iofft + 2 + 0 * ldt ] + " "
    + T[ iofft + 0 + 1 * ldt ] + " "
    + T[ iofft + 1 + 1 * ldt ] + " "
    + T[ iofft + 2 + 1 * ldt ] + " "
    + T[ iofft + 0 + 2 * ldt ] + " "
    + T[ iofft + 1 + 2 * ldt ] + " "
    + T[ iofft + 2 + 2 * ldt ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "Q = "
    + Q[ ioffq + 0 + 0 * ldq ] + " "
    + Q[ ioffq + 1 + 0 * ldq ] + " "
    + Q[ ioffq + 2 + 0 * ldq ] + " "
    + Q[ ioffq + 0 + 1 * ldq ] + " "
    + Q[ ioffq + 1 + 1 * ldq ] + " "
    + Q[ ioffq + 2 + 1 * ldq ] + " "
    + Q[ ioffq + 0 + 2 * ldq ] + " "
    + Q[ ioffq + 1 + 2 * ldq ] + " "
    + Q[ ioffq + 2 + 2 * ldq ]  + "\n";

  n=4;
  T[ iofft + 0 + 0 * ldt ] = 1.;
  T[ iofft + 1 + 0 * ldt ] = 2.;
  T[ iofft + 2 + 0 * ldt ] = 0.;
  T[ iofft + 3 + 0 * ldt ] = 0.;
  T[ iofft + 0 + 1 * ldt ] = -2.;
  T[ iofft + 1 + 1 * ldt ] = 1.;
  T[ iofft + 2 + 1 * ldt ] = 0.;
  T[ iofft + 3 + 1 * ldt ] = 0.;
  T[ iofft + 0 + 2 * ldt ] = 5.;
  T[ iofft + 1 + 2 * ldt ] = 6.;
  T[ iofft + 2 + 2 * ldt ] = 3.;
  T[ iofft + 3 + 2 * ldt ] = 4.;
  T[ iofft + 0 + 3 * ldt ] = 7.;
  T[ iofft + 1 + 3 * ldt ] = 8.;
  T[ iofft + 2 + 3 * ldt ] = -4.;
  T[ iofft + 3 + 3 * ldt ] = 3.;
  LaPack3.dlaexc(false,n,T,ldt,Q,ldq,1,2,2,work,info,
    iofft,ioffq,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaexc(F,...,1,2,2,...): info = " + info.value + "\n";
  document.getElementById("debug_textarea").value +=
    "T = "
    + T[ iofft + 0 + 0 * ldt ] + " "
    + T[ iofft + 1 + 0 * ldt ] + " "
    + T[ iofft + 2 + 0 * ldt ] + " "
    + T[ iofft + 3 + 0 * ldt ] + " "
    + T[ iofft + 0 + 1 * ldt ] + " "
    + T[ iofft + 1 + 1 * ldt ] + " "
    + T[ iofft + 2 + 1 * ldt ] + " "
    + T[ iofft + 3 + 1 * ldt ] + " "
    + T[ iofft + 0 + 2 * ldt ] + " "
    + T[ iofft + 1 + 2 * ldt ] + " "
    + T[ iofft + 2 + 2 * ldt ] + " "
    + T[ iofft + 3 + 2 * ldt ] + " "
    + T[ iofft + 0 + 3 * ldt ] + " "
    + T[ iofft + 1 + 3 * ldt ] + " "
    + T[ iofft + 2 + 3 * ldt ] + " "
    + T[ iofft + 3 + 3 * ldt ]  + "\n";
  T[ iofft + 0 + 0 * ldt ] = 1.;
  T[ iofft + 1 + 0 * ldt ] = 2.;
  T[ iofft + 2 + 0 * ldt ] = 0.;
  T[ iofft + 3 + 0 * ldt ] = 0.;
  T[ iofft + 0 + 1 * ldt ] = -2.;
  T[ iofft + 1 + 1 * ldt ] = 1.;
  T[ iofft + 2 + 1 * ldt ] = 0.;
  T[ iofft + 3 + 1 * ldt ] = 0.;
  T[ iofft + 0 + 2 * ldt ] = 5.;
  T[ iofft + 1 + 2 * ldt ] = 6.;
  T[ iofft + 2 + 2 * ldt ] = 3.;
  T[ iofft + 3 + 2 * ldt ] = 4.;
  T[ iofft + 0 + 3 * ldt ] = 7.;
  T[ iofft + 1 + 3 * ldt ] = 8.;
  T[ iofft + 2 + 3 * ldt ] = -4.;
  T[ iofft + 3 + 3 * ldt ] = 3.;
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      Q[ ioffq + i + j * ldq ] = 0.;
    }
    Q[ ioffq + j + j * ldq ] = 1.;
  }
  LaPack3.dlaexc(true,n,T,ldt,Q,ldq,1,2,2,work,info,
    iofft,ioffq,ioffwork);
  document.getElementById("debug_textarea").value +=
    "dlaexc(T,...,1,2,2,...): info = " + info.value + "\n";
  document.getElementById("debug_textarea").value +=
    "T = "
    + T[ iofft + 0 + 0 * ldt ] + " "
    + T[ iofft + 1 + 0 * ldt ] + " "
    + T[ iofft + 2 + 0 * ldt ] + " "
    + T[ iofft + 3 + 0 * ldt ] + " "
    + T[ iofft + 0 + 1 * ldt ] + " "
    + T[ iofft + 1 + 1 * ldt ] + " "
    + T[ iofft + 2 + 1 * ldt ] + " "
    + T[ iofft + 3 + 1 * ldt ] + " "
    + T[ iofft + 0 + 2 * ldt ] + " "
    + T[ iofft + 1 + 2 * ldt ] + " "
    + T[ iofft + 2 + 2 * ldt ] + " "
    + T[ iofft + 3 + 2 * ldt ] + " "
    + T[ iofft + 0 + 3 * ldt ] + " "
    + T[ iofft + 1 + 3 * ldt ] + " "
    + T[ iofft + 2 + 3 * ldt ] + " "
    + T[ iofft + 3 + 3 * ldt ]  + "\n";
  document.getElementById("debug_textarea").value +=
    "Q = "
    + Q[ ioffq + 0 + 0 * ldq ] + " "
    + Q[ ioffq + 1 + 0 * ldq ] + " "
    + Q[ ioffq + 2 + 0 * ldq ] + " "
    + Q[ ioffq + 3 + 0 * ldq ] + " "
    + Q[ ioffq + 0 + 1 * ldq ] + " "
    + Q[ ioffq + 1 + 1 * ldq ] + " "
    + Q[ ioffq + 2 + 1 * ldq ] + " "
    + Q[ ioffq + 3 + 1 * ldq ] + " "
    + Q[ ioffq + 0 + 2 * ldq ] + " "
    + Q[ ioffq + 1 + 2 * ldq ] + " "
    + Q[ ioffq + 2 + 2 * ldq ] + " "
    + Q[ ioffq + 3 + 2 * ldq ] + " "
    + Q[ ioffq + 0 + 3 * ldq ] + " "
    + Q[ ioffq + 1 + 3 * ldq ] + " "
    + Q[ ioffq + 2 + 3 * ldq ] + " "
    + Q[ ioffq + 3 + 3 * ldq ]  + "\n";
}
