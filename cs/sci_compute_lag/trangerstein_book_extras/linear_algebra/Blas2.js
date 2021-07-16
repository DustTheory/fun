function Blas2() {
}
//************************************************************************
Blas2.xerbla = function( srname, info ) {
  document.getElementById("debug_textarea").value +=
    " ** On entry to " + srname + " parameter number " + info
    + " had an illegal value" + "\n";
}
//*************************************************************************
//  x <-- A x or x <-- A^T x or x <-- A^H x with A triangular
Blas2.dtrmv = function ( uplo, trans, diag, n, A, lda, x, incx, ioffA,
ioffx) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U' &&
  uplo.charAt(0).toUpperCase() != 'L') {
    info = 1;
  } else if ( trans.charAt(0).toUpperCase() != 'N'
  && trans.charAt(0).toUpperCase() != 'T'
  && trans.charAt(0).toUpperCase() != 'C') {
    info = 2;
  } else if ( diag.charAt(0).toUpperCase() != 'U'
  && diag.charAt(0).toUpperCase() != 'N') {
    info = 3;
  } else if (n < 0 ) info = 4;
  else if ( lda < Math.max( 1, n ) ) info = 6;
  else if ( incx == 0 ) info = 8;
  if ( info != 0) {
    Blas2.xerbla( 'dtrmv' , info );
    return;
  }

  if ( n == 0 ) return;
  var nounit = ( diag.charAt(0).toUpperCase() == 'N' );
  var kx = 1;
  if ( incx <= 0 ) kx = 1 - ( n - 1 ) * incx;
  else if ( incx != 1 ) kx = 1;
  var i = -1;
  var ix = -1;
  var j = -1;
  var jx = -1;
  var temp = Number.POSITIVE_INFINITY;
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      jx = kx;
      for ( j = 1; j <= n; j++) {
        if ( x[ ioffx + jx - 1 ] != 0. ) {
          temp = x[ ioffx + jx - 1 ];
          ix = kx;
          for ( i = 1; i <= j - 1; i ++ ) {
            x[ ioffx + ix - 1 ] +=
              temp * A[ ioffA + i - 1 + ( j - 1 ) * lda ];
            ix += incx;
          }
          if (nounit) x[ ioffx + jx - 1 ] *=
            A[ ioffA + j - 1 + ( j - 1 ) * lda ];
        }
        jx += incx;
      }
    } else {
      kx += ( n - 1 ) * incx;
      jx = kx;
      for ( j = n; j >= 1; j -- ) {
        if ( x [ ioffx + jx - 1 ] != 0. ) {
          temp = x[ ioffx + jx - 1 ];
          ix = kx;
          for ( i = n; i >= j + 1; i-- ) {
            x[ ioffx + ix - 1 ] +=
              temp * A[ ioffA + i - 1 + ( j - 1 ) * lda ];
            ix -= incx;
          }
          if ( nounit ) {
            x[ ioffx + jx - 1 ] *=
              A[ ioffA + j - 1 + ( j - 1 ) * lda ];
          }
        }
        jx -= incx;
      }
    }
  } else {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      jx = kx + ( n - 1 ) * incx;
      for ( j = n; j >= 1; j -- ) {
        temp = x[ ioffx + jx - 1 ];
        ix = jx;
        if ( nounit ) temp *= A[ ioffA + j - 1 + ( j - 1 ) * lda ];
        for ( i = j-1; i >= 1; i -- ) {
          ix -= incx;
          temp += A[ ioffA + i - 1 + ( j - 1 ) * lda ]
            * x[ ioffx + ix - 1 ];
        }
        x[ ioffx + jx - 1 ] = temp;
        jx -= incx;
      }
    } else {
      jx = kx;
      for ( j = 1; j <= n; j++ ) {
        temp = x[ ioffx + jx - 1 ];
        ix = jx;
        if ( nounit ) temp *= A[ ioffA + j - 1 + ( j - 1 ) * lda ];
        for ( i = j+1; i <= n; i++ ) {
          ix += incx;
          temp += A[ ioffA + i - 1 + ( j - 1 ) * lda ]
            * x[ ioffx + ix - 1 ];
        }
        x[ ioffx + jx - 1 ] = temp;
        jx += incx;
      }
    }
  }
}
Blas2.ztrmv = function( uplo, trans, diag, n, A, lda, x, incx, ioffA,
ioffx) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U'
  && uplo.charAt(0).toUpperCase() != 'L') {
    info = 1;
  } else if ( trans.charAt(0).toUpperCase() != 'N'
  && trans.charAt(0).toUpperCase() != 'T'
  && trans.charAt(0).toUpperCase() != 'C') {
    info = 2;
  } else if ( diag.charAt(0).toUpperCase() != 'U'
  && diag.charAt(0).toUpperCase() != 'N') {
    info = 3;
  } else if (n < 0 ) info = 4;
  else if ( lda < Math.max( 1, n ) ) info = 6;
  else if ( incx == 0 ) info = 8;
  if ( info != 0) {
    Blas2.xerbla( 'ztrmv' , info );
    return;
  }

  if ( n == 0 ) return;
  var noconj = ( trans.charAt(0).toUpperCase() == 'T' );
  var nounit = ( diag.charAt(0).toUpperCase() == 'N' );
  var i = -1;
  var ix = -1;
  var j = -1;
  var jx = -1;
  var kx = 1;
  var temp = new Complex();
  var zero = new Complex( 0., 0. );
  if ( incx <= 0 ) kx = 1 - ( n - 1 ) * incx;
  else if ( incx != 1 ) kx = 1;
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      jx = kx;
      for ( j = 1; j <= n; j++) {
        if ( ! x[ ioffx + jx - 1 ].equals( zero ) ) {
          temp.setValue( x[ ioffx + jx - 1 ] );
          ix = kx;
          for ( i = 1; i <= j - 1; i ++ ) {
            x[ ioffx + ix - 1 ].plusEquals( ComplexMath.times( temp ,
              A[ ioffA + i - 1 + ( j - 1 ) * lda ] ) );
            ix += incx;
          }
          if (nounit) {
            x[ ioffx + jx - 1 ].timesEquals(
              A[ ioffA + j - 1 + ( j - 1 ) * lda ] );
          }
        }
        jx += incx;
      }
    } else {
      kx += ( n - 1 ) * incx;
      jx = kx;
      for ( j = n; j >= 1; j -- ) {
        if ( ! x [ ioffx + jx - 1 ].equals( zero ) ) {
          temp.setValue( x[ ioffx + jx - 1 ] );
          ix = kx;
          for ( i = n; i >= j + 1; i-- ) {
            x[ ioffx + ix - 1 ].plusEquals( ComplexMath.times( temp ,
              A[ ioffA + i - 1 + ( j - 1 ) * lda ]) );
            ix -= incx;
          }
          if ( nounit ) {
            x[ ioffx + jx - 1 ].timesEquals(
              A[ ioffA + j - 1 + ( j - 1 ) * lda ] );
          }
        }
        jx -= incx;
      }
    }
  } else {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      jx = kx + ( n - 1 ) * incx;
      for ( j = n; j >= 1; j -- ) {
        temp.setValue( x[ ioffx + jx - 1 ] );
        ix = jx;
        if ( noconj ) {
          if ( nounit ) {
            temp.timesEquals( A[ ioffA + j - 1 + ( j - 1 ) * lda ] );
          }
          for ( i = j-1; i >= 1; i -- ) {
            ix -= incx;
            temp.plusEquals( ComplexMath.times(
              A[ ioffA + i - 1 + ( j - 1 ) * lda ] ,
              x[ ioffx + ix - 1 ] ) );
          }
        } else {
          if ( nounit ) {
            temp.timesEquals( ComplexMath.conj(
              A[ ioffA + j - 1 + ( j - 1 ) * lda ] ) );
          }
          for ( i = j-1; i >= 1; i -- ) {
            ix -= incx;
            temp.plusEquals( ComplexMath.times( ComplexMath.conj(
              A[ ioffA + i - 1 + ( j - 1 ) * lda ] ) ,
              x[ ioffx + ix - 1 ] ) );
          }
        }
        x[ ioffx + jx - 1 ].setValue( temp );
        jx -= incx;
      }
    } else {
      jx = kx;
      for ( j = 1; j <= n; j++ ) {
        temp.setValue( x[ ioffx + jx - 1 ] );
        ix = jx;
        if ( noconj ) {
          if ( nounit ) {
            temp.timesEquals( A[ ioffA + j - 1 + ( j - 1 ) * lda ] );
          }
          for ( i = j+1; i <= n; i++ ) {
            ix += incx;
            temp.plusEquals( ComplexMath.times(
              A[ ioffA + i - 1 + ( j - 1 ) * lda ] ,
              x[ ioffx + ix - 1 ] ) );
          }
        } else {
          if ( nounit ) {
            temp.timesEquals( ComplexMath.conj(
              A[ ioffA + j - 1 + ( j - 1 ) * lda ] ) );
          }
          for ( i = j+1; i <= n; i++ ) {
            ix += incx;
            temp.plusEquals( ComplexMath.times( ComplexMath.conj(
              A[ ioffA + i - 1 + ( j - 1 ) * lda ] ) ,
              x[ ioffx + ix - 1 ] ) );
          }
        }
        x[ ioffx + jx - 1 ].setValue( temp );
        jx += incx;
      }
    }
  }
}
//*************************************************************************
//  x <-- A x or x <-- A^T x, A triangular banded
Blas2.dtbmv = function( uplo, trans, diag, n, k, A, lda, x, incx, ioffA,
ioffx) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U'
  && uplo.charAt(0).toUpperCase() != 'L') {
    info = 1;
  } else if ( trans.charAt(0).toUpperCase() != 'N'
  && trans.charAt(0).toUpperCase() != 'T'
  && trans.charAt(0).toUpperCase() != 'C') {
    info = 2;
  } else if ( diag.charAt(0).toUpperCase() != 'U'
  && diag.charAt(0).toUpperCase() != 'N') {
    info = 3;
  } else if ( n < 0 ) info = 4;
  else if ( k < 0 ) info = 5;
  else if ( lda < k+1 ) info = 7;
  else if ( incx == 0 ) info = 9;
  if ( info != 0 ) {
    Blas2.xerbla( 'dtbmv', info );
    return;
  }

  if ( n == 0 ) return;
  var nounit = ( diag.charAt(0).toUpperCase() == 'N' );
  var kx = 1;
  if ( incx <= 0 ) kx = 1 - ( n - 1 ) * incx;
  else if ( incx != 1 ) kx = 1;

  var i = -1;
  var ix = -1;
  var j = -1;
  var jx = -1;
  var l = -1;
  var temp = Number.POSITIVE_INFINITY;
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      var kplus1 = k + 1;
      jx = kx;
      for ( j = 1; j <= n; j ++ ) {
        if ( x[ ioffx + jx - 1 ] != 0. ) {
          temp = x[ ioffx + jx - 1 ];
          ix = kx;
          l = kplus1 - j;
          for ( i = Math.max( 1, j - k ); i <= j - 1; i ++ ) {
            x[ ioffx + ix - 1 ] +=
              temp * A[ ioffA + l + i - 1 + ( j - 1 ) * lda ];
            ix += incx;
          }
          if ( nounit ) {
            x[ ioffx + jx - 1 ] *=
              A[ ioffA + kplus1 - 1 + ( j - 1 ) * lda ];
          }
        }
        jx += incx;
        if ( j > k ) kx += incx;
      }
    } else {
      kx += ( n - 1 ) * incx;
      jx = kx;
      for ( j = n; j >= 1; j -- ) {
        if ( x[ ioffx + jx - 1 ] != 0. ) {
          temp = x[ ioffx + jx - 1 ];
          ix = kx;
          l = 1 - j;
          for ( i = Math.min( n , j + k ); i >= j+1; i -- ) {
            x[ ioffx + ix - 1 ] +=
              temp * A[ ioffA + l + i - 1 + ( j - 1 ) * lda ];
            ix -= incx;
          }
          if ( nounit ) {
            x[ ioffx + jx - 1 ] *= A[ ioffA + ( j - 1 ) * lda ];
          }
        }
        jx -= incx;
        if ( n - j >= k ) kx -= incx;
      }
    }
  } else {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      kplus1 = k + 1;
      kx += ( n - 1 ) * incx;
      jx = kx;
      for ( j = n; j >= 1; j -- ) {
        temp = x[ ioffx + jx - 1 ];
        kx -= incx;
        ix = kx;
        l = kplus1 - j;
        if ( nounit ) {
          temp *= A[ ioffA + kplus1 - 1 + ( j - 1 ) * lda ];
        }
        for ( i = j - 1; i >= Math.max( 1, j - k ); i -- ) {
          temp += A[ ioffA + l + i - 1 + ( j - 1 ) * lda ]
            * x[ ioffx + ix - 1 ];
          ix -= incx;
        }
        x[ ioffx + jx - 1 ] = temp;
        jx -= incx;
      }
    } else {
      jx = kx;
      for ( j = 1; j <= n; j ++ ) {
        temp = x[ ioffx + jx - 1 ];
        kx += incx;
        ix = kx;
        l = 1 - j;
        if ( nounit) temp *= A[ ioffA + ( j - 1 ) * lda ];
        for ( i = j + 1; i <= Math.min( n, j + k); i ++ ) {
          temp += A[ ioffA + l + i - 1 + ( j - 1 ) * lda ]
            * x[ ioffx + ix - 1 ];
          ix += incx;
        }
        x[ ioffx + jx - 1 ] = temp;
        jx += incx;
      }
    }
  }
}
Blas2.ztbmv = function( uplo, trans, diag, n, k, A, lda, x, incx, ioffA,
ioffx) {
  throw new Error("not tested: banded matrix");
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U'
  && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( trans.charAt(0).toUpperCase() != 'N'
  && trans.charAt(0).toUpperCase() != 'T'
  && trans.charAt(0).toUpperCase() != 'C') {
    info = 2;
  } else if ( diag.charAt(0).toUpperCase() != 'U'
  && diag.charAt(0).toUpperCase() != 'N' ) {
    info = 3;
  } else if ( n < 0 ) info = 4;
  else if ( k < 0 ) info = 5;
  else if ( lda < k+1 ) info = 7;
  else if ( incx == 0 ) info = 9;
  if ( info != 0 ) {
    Blas2.xerbla( 'ztbmv', info );
    return;
  }

  if ( n == 0 ) return;
  var noconj = ( trans.charAt(0).toUpperCase() == 'T' );
  var nounit = ( diag.charAt(0).toUpperCase() == 'N' );
  var kx = 1;
  if ( incx <= 0 ) kx = 1 - ( n - 1 ) * incx;
  else if ( incx != 1 ) kx = 1;

  var i = -1;
  var ix = -1;
  var j = -1;
  var jx = -1;
  var l = -1;
  var temp = new Complex();
  var zero = new Complex( 0., 0. );
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      var kplus1 = k + 1;
      jx = kx;
      for ( j = 1; j <= n; j ++ ) {
        if ( ! x[ ioffx + jx - 1 ].equals( zero ) ) {
          temp.setValue( x[ ioffx + jx - 1 ] );
          ix = kx;
          l = kplus1 - j;
          for ( i = Math.max( 1, j - k ); i <= j-1; i ++ ) {
            x[ ioffx + ix - 1 ].plusEquals( ComplexMath.times( temp ,
              A[ ioffA + l + i - 1 + ( j - 1 ) * lda ] ) );
            ix += incx;
          }
          if ( nounit ) {
            x[ ioffx + jx - 1 ].timesEquals(
              A[ ioffA + kplus1 - 1 + (j - 1 ) * lda ] );
          }
        }
        jx += incx;
        if ( j > k ) kx += incx;
      }
    } else {
      kx += ( n - 1 ) * incx;
      jx = kx;
      for ( j = n; j >= 1; j -- ) {
        if ( ! x[ ioffx + jx - 1 ].equals( zero ) ) {
          temp.setValue( x[ ioffx + jx - 1 ] );
          ix = kx;
          l = 1 - j;
          for ( i = Math.min( n , j + k ); i >= j+1; i -- ) {
            x[ ioffx + ix - 1 ].plusEquals( ComplexMath.times( temp ,
              A[ ioffA + l + i - 1 + ( j - 1 ) * lda ] ) );
            ix -= incx;
          }
          if ( nounit ) {
            x[ ioffx + jx - 1 ].timesEquals(
              A[ ioffA + ( j - 1 ) * lda ] );
          }
        }
        jx -= incx;
        if ( n - j >= k ) kx -= incx;
      }
    }
  } else {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      kplus1 = k + 1;
      kx += ( n - 1 ) * incx;
      jx = kx;
      for ( j = n; j >= 1; j -- ) {
        temp.setValue( x[ ioffx + jx - 1 ] );
        kx -= incx;
        ix = kx;
        l = kplus1 - j;
        if ( noconj ) {
          if ( nounit ) {
            temp.timesEquals(
              A[ ioffA + kplus1 - 1 + ( j - 1 ) * lda ] );
          }
          for ( i = j - 1; i >= Math.max( 1, j - k ); i -- ) {
            temp.plusEquals( ComplexMath.times(
              A[ ioffA + l + i - 1 + ( j - 1 ) * lda ] ,
              x[ ioffx + ix - 1 ] ) );
            ix -= incx;
          }
        } else {
          if ( nounit ) temp.timesEquals( ComplexMath.conj(
            A[ ioffA + kplus1 - 1 + ( j - 1 ) * lda ] ) );
          for ( i = j - 1; i >= Math.max( 1, j - k ); i -- ) {
            temp.plusEquals( ComplexMath.times( ComplexMath.conj(
              A[ ioffA + l + i - 1 + ( j - 1 ) * lda ] ) ,
              x[ ioffx + ix - 1 ] ) );
            ix -= incx;
          }
        }
        x[ ioffx + jx - 1 ].setValue( temp );
        jx -= incx;
      }
    } else {
      jx = kx;
      for ( j = 1; j <= n; j ++ ) {
        temp.setValue( x[ ioffx + jx - 1 ] );
        kx += incx;
        ix = kx;
        l = 1 - j;
        if ( noconj ) {
          if ( nounit) {
            temp.timesEquals( A[ ioffA + ( j - 1 ) * lda ] );
          }
          for ( i = j + 1; i <= Math.min( n, j + k); i ++ ) {
            temp.plusEquals( ComplexMath.times(
              A[ ioffA + l + i - 1 + ( j - 1 ) * lda ] ,
              x[ ioffx + ix - 1 ] ) );
            ix += incx;
          }
        } else {
          if ( nounit) {
            temp.timesEquals( ComplexMath.conj(
              A[ ioffA + ( j - 1 ) * lda ] ) );
          }
          for ( i = j + 1; i <= Math.min( n, j + k); i ++ ) {
            temp.plusEquals( ComplexMath.times( ComplexMath.conj( 
              A[ ioffA + l + i - 1 + ( j - 1 ) * lda ] ) ,
              x[ ioffx + ix - 1 ] ) );
            ix += incx;
          }
        }
        x[ ioffx + jx - 1 ].setValue( temp );
        jx += incx;
      }
    }
  }
}
//*************************************************************************
//  x <-- A x or x <-- A^T x, A triangular packed
Blas2.dtpmv = function( uplo, trans, diag, n, AP, x, incx, ioffAP,
ioffx ) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U'
  && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( trans.charAt(0).toUpperCase() != 'N'
  && trans.charAt(0).toUpperCase() != 'T'
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info = 2;
  } else if ( diag.charAt(0).toUpperCase() != 'U'
  && diag.charAt(0).toUpperCase() != 'N' ) {
    info = 3;
  } else if ( n < 0 ) info = 4;
  else if ( incx == 0 ) info = 7;
  if ( info != 0 ) {
    Blas2.xerbla( 'dtpmv', info );
    return;
  }

  if ( n == 0 ) return;
  var nounit = ( diag.charAt(0).toUpperCase() == 'N' );
  var kx = 1;
  if ( incx <= 0 ) kx = 1 - ( n - 1 ) * incx;
  else if ( incx != 1 ) kx = 1;

  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    var ix = -1;
    var j = -1;
    var jx = -1;
    var k = -1;
    var kk = -1;
    var temp = Number.POSITIVE_INFINITY;
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      kk = 1;
      jx = kx;
      for ( j = 1; j <= n; j ++ ) {
        if ( x[ ioffx + jx - 1 ] != 0. ) {
          temp = x[ ioffx + jx - 1 ];
          ix = kx;
          for ( k = kk; k <= kk + j - 2 ; k ++ ) {
            x[ ioffx + ix - 1 ] += temp * AP[ ioffAP + k - 1 ];
            ix += incx;
          }
          if ( nounit) {
            x[ ioffx + jx - 1 ] *= AP[ ioffAP + kk + j - 2 ];
          }
        }
        jx += incx;
        kk += j;
      }
    } else {
      kk = ( n * ( n + 1 ) ) / 2;
      kx += ( n - 1 ) * incx;
      jx = kx;
      for ( j = n; j >= 1; j -- ) {
        if ( x[ ioffx + jx - 1 ] != 0. ) {
          temp = x[ ioffx + jx - 1 ];
          ix = kx;
          for ( k = kk; k >= kk - n + j + 1; k -- ) {
            x[ ioffx + ix - 1 ] += temp * AP[ ioffAP + k - 1 ];
            ix -= incx;
          }
          if ( nounit ) {
            x[ ioffx + jx - 1 ] *= AP[ ioffAP + kk - n + j - 1 ];
          }
        }
        jx -= incx;
        kk -= n - j + 1;
      }
    }
  } else {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      kk = ( n * ( n + 1 ) ) / 2;
      jx = kx + ( n - 1 ) * incx;
      for ( j = n; j >= 1; j -- ) {
        temp = x[ ioffx + jx - 1 ];
        ix = jx;
        if ( nounit ) temp *= AP[ ioffAP + kk - 1 ];
        for ( k = kk - 1; k >= kk - j + 1; k -- ) {
          ix -= incx;
          temp += AP[ ioffAP + k - 1 ] * x[ ioffx + ix - 1 ];
        }
        x[ ioffx + jx - 1 ] = temp;
        jx -= incx;
        kk -= j;
      }
    } else {
      kk = 1;
      jx = kx;
      for ( j = 1; j <= n; j ++ ) {
        temp = x[ ioffx + jx - 1 ];
        ix = jx;
        if ( nounit) temp *= AP[ ioffAP + kk - 1 ];
        for ( k = kk + 1; k <= kk + n - j; k ++ ) {
          ix += incx;
          temp += AP[ ioffAP + k - 1 ] * x[ ioffx + ix - 1 ];
        }
        x[ ioffx + jx - 1 ] = temp;
        jx += incx;
        kk += n - j + 1;
      }
    }
  }
}
Blas2.ztpmv = function( uplo, trans, diag, n, AP, x, incx, ioffAP,
ioffx ) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U'
  && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( trans.charAt(0).toUpperCase() != 'N'
  && trans.charAt(0).toUpperCase() != 'T'
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info = 2;
  } else if ( diag.charAt(0).toUpperCase() != 'U'
  && diag.charAt(0).toUpperCase() != 'N' ) {
    info = 3;
  } else if ( n < 0 ) info = 4;
  else if ( incx == 0 ) info = 7;
  if ( info != 0 ) {
    Blas2.xerbla( 'ztpmv', info );
    return;
  }

  if ( n == 0 ) return;
  var noconj = ( trans.charAt(0).toUpperCase() == 'T' );
  var nounit = ( diag.charAt(0).toUpperCase() == 'N' );
  var kx = 1;
  if ( incx <= 0 ) kx = 1 - ( n - 1 ) * incx;
  else if ( incx != 1 ) kx = 1;

  var temp = new Complex( );
  var zero = new Complex( 0., 0. );
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    var ix = -1;
    var j = -1;
    var jx = -1;
    var k = -1;
    var kk = -1;
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      kk = 1;
      jx = kx;
      for ( j = 1; j <= n; j ++ ) {
        if ( ! x[ ioffx + jx - 1 ].equals( zero ) ) {
          temp.setValue( x[ ioffx + jx - 1 ] );
          ix = kx;
          for ( k = kk; k <= kk + j - 2; k ++ ) {
            x[ ioffx + ix - 1 ].plusEquals( ComplexMath.times( temp ,
              AP[ ioffAP + k - 1 ] ));
            ix += incx;
          }
          if ( nounit) {
            x[ ioffx + jx - 1 ].timesEquals(
              AP[ ioffAP + kk + j - 2 ] );
          }
        }
        jx += incx;
        kk += j;
      }
    } else {
      kk = ( n * ( n + 1 ) ) / 2;
      kx += ( n - 1 ) * incx;
      jx = kx;
      for ( j = n; j >= 1; j -- ) {
        if ( ! x[ ioffx + jx - 1 ].equals( zero ) ) {
          temp.setValue( x[ ioffx + jx - 1 ] );
          ix = kx;
          for ( k = kk; k >= kk - n + j + 1; k -- ) {
            x[ ioffx + ix - 1 ].plusEquals( ComplexMath.times( temp ,
              AP[ ioffAP + k - 1 ] ));
            ix -= incx;
          }
          if ( nounit ) {
            x[ ioffx + jx - 1 ].
              timesEquals( AP[ ioffAP + kk - n + j - 1 ] );
          }
        }
        jx -= incx;
        kk -= n - j + 1;
      }
    }
  } else {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      kk = ( n * ( n + 1 ) ) / 2;
      jx = kx + ( n - 1 ) * incx;
      for ( j = n; j >= 1; j -- ) {
        temp.setValue( x[ ioffx + jx - 1 ] );
        ix = jx;
        if ( noconj ) {
          if ( nounit ) temp.timesEquals( AP[ ioffAP + kk - 1 ] );
          for ( k = kk - 1; k >= kk - j + 1; k -- ) {
            ix -= incx;
            temp.plusEquals( ComplexMath.times( AP[ ioffAP + k - 1 ] ,
              x[ ioffx + ix - 1 ] ));
          }
        } else {
          if ( nounit ) {
            temp.timesEquals( ComplexMath.conj(
              AP[ ioffAP + kk - 1 ] ) );
          }
          for ( k = kk - 1; k >= kk - j + 1; k -- ) {
            ix -= incx;
            temp.plusEquals( ComplexMath.times( ComplexMath.conj(
              AP[ ioffAP + k - 1 ] ) , x[ ioffx + ix - 1 ] ));
          }
        }
        x[ ioffx + jx - 1 ].setValue( temp );
        jx -= incx;
        kk -= j;
      }
    } else {
      kk = 1;
      jx = kx;
      for ( j = 1; j <= n; j ++ ) {
        temp.setValue( x[ ioffx + jx - 1 ] );
        ix = jx;
        if ( noconj ) {
          if ( nounit) temp.timesEquals( AP[ ioffAP + kk - 1 ] );
          for ( k = kk + 1; k <= kk + n - j; k ++ ) {
            ix += incx;
            temp.plusEquals( ComplexMath.times( AP[ ioffAP + k - 1 ] ,
              x[ ioffx + ix - 1 ] ));
          }
        } else {
          if ( nounit) {
            temp.timesEquals(ComplexMath.conj( AP[ ioffAP + kk - 1 ]));
          }
          for ( k = kk + 1; k <= kk + n - j; k ++ ) {
            ix += incx;
            temp.plusEquals( ComplexMath.times( ComplexMath.conj(
              AP[ ioffAP + k - 1 ] ) , x[ ioffx + ix - 1 ] ));
          }
        }
        x[ ioffx + jx - 1 ].setValue( temp );
        jx += incx;
        kk += n - j + 1;
      }
    }
  }
}
//*************************************************************************
//  X <-- A^{-1} x or x <-- A^{-T} x or x <-- A^{-H} x
Blas2.dtrsv = function( uplo, trans, diag, n, A, lda, x, incx, ioffA,
ioffx) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U'
  && uplo.charAt(0).toUpperCase() != 'L') {
    info = 1;
  } else if ( trans.charAt(0).toUpperCase() != 'N'
  && trans.charAt(0).toUpperCase() != 'T'
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info = 2;
  } else if ( diag.charAt(0).toUpperCase() != 'U'
  && diag.charAt(0).toUpperCase() != 'N' ) {
    info = 3;
  } else if (n < 0 ) info = 4;
  else if ( lda < Math.max( 1, n ) ) info = 6;
  else if ( incx == 0 ) info = 8;
  if ( info != 0){
    Blas2.xerbla( 'dtrsv' , info );
    return;
  }

  if ( n == 0 ) return;
  var nounit = ( diag.charAt(0).toUpperCase() == 'N' );
  var i = -1;
  var ix = -1;
  var j = -1;
  var jx = -1;
  var kx = 1;
  var temp = Number.POSITIVE_INFINITY;
  if ( incx <= 0 ) kx = 1 - ( n - 1 ) * incx;
  else if ( incx != 1 ) kx = 1;
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      jx = kx + ( n - 1 ) * incx;
      for ( j = n; j >= 1; j -- ) {
        if ( x[ ioffx + jx - 1 ] != 0. ) {
          if ( nounit ) {
            x[ ioffx + jx - 1 ] /=
              A[ ioffA + j - 1 + ( j - 1 ) * lda ];

            var ij = ioffx + jx - 1;

          }
          temp = x [ ioffx + jx - 1 ];
          ix = jx;
          for ( i = j-1; i >= 1; i -- ) {
            ix -= incx;
            x[ ioffx + ix - 1 ] -=
              temp * A[ ioffA + i - 1 + ( j - 1 ) * lda ];
          }
        }
        jx -= incx;
      }
    } else {
      jx = kx;
      for ( j = 1; j <= n; j ++ ) {
        if ( x[ ioffx + jx - 1 ] != 0. ) {
          if ( nounit ) {
            x[ ioffx + jx - 1 ] /=
              A[ ioffA + j - 1 + ( j - 1 ) * lda ];
          }
          temp = x [ ioffx + jx - 1 ];
          ix = jx;
          for ( i = j + 1 ; i <= n; i ++ ) {
            ix += incx;
            x[ ioffx + ix - 1 ] -=
              temp * A [ ioffA + i - 1 + ( j - 1 ) * lda ];
          }
        }
        jx += incx;
      }
    }
  } else {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      jx = kx;
      for ( j = 1; j <= n; j++ ) {
        temp = x [ ioffx + jx - 1 ];
        ix = kx;
        for ( i = 1; i <= j - 1; i ++ ) {
          temp -= A[ ioffA + i - 1 + ( j - 1 ) * lda ]
            * x[ ioffx + ix - 1 ];
          ix += incx;
        }
        if ( nounit ) temp /= A[ ioffA + j - 1 + ( j - 1 ) * lda ];
        x[ ioffx + jx - 1 ] = temp;
        jx += incx;
      }
    } else {
      kx += ( n - 1 ) * incx;
      jx = kx;
      for ( j = n; j >= 1; j -- ) {
        temp = x[ ioffx + jx - 1 ];
        ix = kx;
        for ( i = n; i >= j + 1; i -- ) {
          temp -= A[ ioffA + i - 1 + ( j - 1 ) * lda ]
            * x[ ioffx + ix - 1 ];
          ix -= incx;
        }
        if ( nounit ) temp /= A[ ioffA + j - 1 + ( j - 1 ) * lda ];
        x[ ioffx + jx - 1 ] = temp;
        jx -= incx;
      }
    }
  }
}
Blas2.ztrsv = function( uplo, trans, diag, n, A, lda, x, incx, ioffA,
ioffx) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U'
  && uplo.charAt(0).toUpperCase() != 'L') {
    info = 1;
  } else if ( trans.charAt(0).toUpperCase() != 'N'
  && trans.charAt(0).toUpperCase() != 'T'
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info = 2;
  } else if ( diag.charAt(0).toUpperCase() != 'U'
  && diag.charAt(0).toUpperCase() != 'N' ) {
    info = 3;
  } else if (n < 0 ) info = 4;
  else if ( lda < Math.max( 1, n ) ) info = 6;
  else if ( incx == 0 ) info = 8;
  if ( info != 0) {
    Blas2.xerbla( 'ztrsv' , info );
    return;
  }

  if ( n == 0 ) return;
  var noconj = ( trans.charAt(0).toUpperCase() == 'T' );
  var nounit = ( diag.charAt(0).toUpperCase() == 'N' );
  var i = -1;
  var ix = -1;
  var j = -1;
  var jx = -1;
  var kx = 1;
  var temp = new Complex();
  var zero = new Complex( 0., 0. );
  if ( incx < 0 ) kx = 1 - ( n - 1 ) * incx;
  else if ( incx != 1 ) kx = 1;
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      jx = kx + ( n - 1 ) * incx;
      for ( j = n; j >= 1; j -- ) {
        if ( ! x[ ioffx + jx - 1 ].equals( zero ) ) {
          if ( nounit ) {
            x[ ioffx + jx - 1 ].divideEquals(
              A[ ioffA + j - 1 + ( j - 1 ) * lda ] );
          }
          temp.setValue( x [ ioffx + jx - 1 ] );
          ix = jx;
          for ( i = j-1; i >= 1; i -- ) {
            ix -= incx;
            x[ ioffx + ix - 1 ].minusEquals(
              ComplexMath.times( temp ,
              A[ ioffA + i - 1 + ( j - 1 ) * lda ] ) );
          }
        }
        jx -= incx;
      }
    } else {
      jx = kx;
      for ( j = 1; j <= n; j ++ ) {
        if ( ! x[ ioffx + jx - 1 ].equals( zero ) ) {
          if ( nounit ) {
            x[ ioffx + jx - 1 ].divideEquals(
              A[ ioffA + j - 1 + ( j - 1 ) * lda ] );
          }
          temp.setValue( x [ ioffx + jx - 1 ] );
          ix = jx;
          for ( i = j + 1 ; i <= n; i ++ ) {
            ix += incx;
            x[ ioffx + ix - 1 ].minusEquals( ComplexMath.times( temp ,
              A [ ioffA + i - 1 + ( j - 1 ) * lda ] ) );
          }
        }
        jx += incx;
      }
    }
  } else {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      jx = kx;
      for ( j = 1; j <= n; j++ ) {
        ix = kx;
        temp.setValue( x [ ioffx + jx - 1 ] );
        if ( noconj ) {
          for ( i = 1; i <= j - 1; i ++ ) {
            temp.minusEquals( ComplexMath.times(
              A[ ioffA + i - 1 + ( j - 1 ) * lda ] ,
              x[ ioffx + ix - 1 ] ) );
            ix += incx;
          }
          if ( nounit ) {
            temp.divideEquals( A[ ioffA + j - 1 + ( j - 1 ) * lda ] );
          }
        } else {
          for ( i = 1; i <= j - 1; i ++ ) {
            temp.minusEquals( ComplexMath.times( ComplexMath.conj(
              A[ ioffA + i - 1 + (j - 1 ) * lda ] ) ,
              x[ ioffx + ix - 1 ] ) );
            ix += incx;
          }
          if ( nounit ) {
            temp.divideEquals( ComplexMath.conj(
              A[ ioffA + j - 1 + ( j - 1 ) * lda ] ) );
          }
        }
        x[ ioffx + jx - 1 ].setValue( temp );
        jx += incx;
      }
    } else {
      kx += ( n - 1 ) * incx;
      jx = kx;
      for ( j = n; j >= 1; j -- ) {
        ix = kx;
        temp.setValue( x[ ioffx + jx - 1 ] );
        if ( noconj ) {
          for ( i = n; i >= j + 1; i -- ) {
            temp.minusEquals( ComplexMath.times(
              A[ ioffA + i - 1 + ( j - 1 ) * lda ] ,
              x[ ioffx + ix - 1 ] ) );
            ix -= incx;
          }
          if ( nounit ) {
            temp.divideEquals( A[ ioffA + j - 1 + ( j - 1 ) * lda ] );
          }
        } else {
          for ( i = n; i >= j + 1; i -- ) {
            temp.minusEquals( ComplexMath.times( ComplexMath.conj(
              A[ ioffA + i - 1 + ( j - 1 ) * lda ] ) ,
              x[ ioffx + ix - 1 ] ) );
            ix -= incx;
          }
          if ( nounit ) {
            temp.divideEquals( ComplexMath.conj(
              A[ ioffA + j - 1 + ( j - 1 ) * lda ] ) );
          }
        }
        x[ ioffx + jx - 1 ].setValue( temp );
        jx -= incx;
      }
    }
  }
}
//*************************************************************************
//  x <-- A^{-1} x or x <-- A^{-T} x, A triangular banded
Blas2.dtbsv = function( uplo, trans, diag, n, k, A, lda, x, incx, ioffA,
ioffx) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U'
  && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( trans.charAt(0).toUpperCase() != 'N'
  && trans.charAt(0).toUpperCase() != 'T'
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info = 2;
  } else if ( diag.charAt(0).toUpperCase() != 'U'
  && diag.charAt(0).toUpperCase() != 'N') {
    info = 3;
  } else if (n < 0 ) info = 4;
  else if ( lda < k + 1 ) info = 7;
  else if ( incx == 0 ) info = 9;
  if ( info != 0) {
    Blas2.xerbla( 'dtbsv' , info );
    return;
  }

  if ( n == 0 ) return;
  var nounit = ( diag.charAt(0).toUpperCase() == 'N' );
  var i = -1;
  var ix = -1;
  var j = -1;
  var jx = -1;
  var l = -1;
  var lx = -1;
  var kx = 1;
  if ( incx <= 0 ) kx = 1 - ( n - 1 ) * incx;
  else if ( incx != 1 ) kx = 1;
  var temp = Number.POSITIVE_INFINITY;
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      var kplus1 = k + 1;
      kx += ( n - 1 ) * incx;
      jx = kx;
      for ( j = n; j >= 1; j -- ) {
        kx -= incx;
        if ( x[ ioffx + jx - 1 ] != 0. ) {
          ix = kx;
          l = kplus1 - j;
          if ( nounit) {
            x[ ioffx + jx - 1 ] /=
              A[ ioffA + kplus1 - 1 + ( j - 1 ) * lda ];
          }
          temp = x[ ioffx + jx - 1 ];
          for ( i = j-1; i >= Math.max( 1, j - k ); i -- ) {
            x[ ioffx + ix - 1 ] -= temp *
              A[ ioffA + l + i - 1 + ( j - 1 ) * lda ];
            ix -= incx;
          }
        }
        jx -= incx;
      }
    } else {
      jx = kx;
      for ( j = 1; j <= n; j ++ ) {
        kx += incx;
        if ( x[ ioffx + jx - 1 ] != 0. ) {
          ix = kx;
          l = 1 - j;
          if ( nounit ) {
            x[ ioffx + jx - 1 ] /= A[ ioffA + ( j - 1 ) *lda ];
          }
          temp = x[ ioffx + jx - 1 ];
          for ( i = j + 1; i <= Math.min( n , j + k ); i ++ ) {
            x[ ioffx + ix - 1 ] -=
              temp * A[ ioffA + l + i - 1 + ( j - 1 ) * lda ];
            ix += incx;
          }
        }
        jx += incx;
      }
    }
  } else {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      kplus1 = k + 1;
      jx = kx;
      for ( j = 1; j <= n; j ++ ) {
        temp = x[ ioffx + jx - 1 ];
        ix = kx;
        l = kplus1 - j;
        for ( i = Math.max( 1, j - k ); i <= j - 1; i ++ ) {
          temp -= A[ ioffA + l + i - 1 + ( j - 1 ) * lda ]
            * x[ ioffx + ix - 1 ];
          ix += incx;
        }
        if ( nounit ) {
          temp /= A[ ioffA + kplus1 - 1 + ( j - 1 ) * lda ];
        }
        x[ ioffx + jx - 1 ] = temp;
        jx += incx;
        if ( j > k ) kx += incx;
      }
    } else {
      kx += ( n - 1 ) * incx;
      jx = kx;
      for ( j = n; j >= 1; j -- ) {
        temp = x[ ioffx + jx - 1 ];
        ix = kx;
        l = 1 - j;
        for ( i = Math.min( n , j + k ); i >= j + 1; i -- ) {
          temp -= A[ ioffA + l + i - 1 + ( j - 1 ) * lda ]
            * x[ ioffx + ix - 1 ];
          ix -= incx;
        }
        if ( nounit ) temp /= A[ ioffA + ( j - 1 ) * lda ];
        x[ ioffx + jx - 1 ] = temp;
        jx -= incx;
        if ( n - j >= k ) kx -= incx;
      }
    }
  }
}
Blas2.ztbsv = function( uplo, trans, diag, n, k, A, lda, x, incx, ioffA,
ioffx) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U'
  && uplo.charAt(0).toUpperCase() != 'L') {
    info = 1;
  } else if ( trans.charAt(0).toUpperCase() != 'N'
  && trans.charAt(0).toUpperCase() != 'T'
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info = 2;
  } else if ( diag.charAt(0).toUpperCase() != 'U'
  && diag.charAt(0).toUpperCase() != 'N') {
    info = 3;
  } else if (n < 0 ) info = 4;
  else if ( lda < k + 1 ) info = 7;
  else if ( incx == 0 ) info = 9;
  if ( info != 0) Blas2.xerbla( 'ztbsv' , info );

  if ( n == 0 ) return;
  var noconj = ( trans.charAt(0).toUpperCase() == 'T' );
  var nounit = ( diag.charAt(0).toUpperCase() == 'N' );
  var i = -1;
  var ix = -1;
  var j = -1;
  var jx = -1;
  var l = -1;
  var lx = -1;
  var kx = 1;
  if ( incx <= 0 ) kx = 1 - ( n - 1 ) * incx;
  else if ( incx != 1 ) kx = 1;
  var temp = new Complex();
  var zero = new Complex( 0. , 0. );
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      var kplus1 = k + 1;
      kx += ( n - 1 ) * incx;
      jx = kx;
      for ( j = n; j >= 1; j -- ) {
        kx -= incx;
        if ( ! x[ ioffx + jx - 1 ].equals( zero ) ) {
          ix = kx;
          l = kplus1 - j;
          if ( nounit) {
            x[ ioffx + jx - 1 ].divideEquals(
              A[ ioffA + kplus1 - 1 + ( j - 1 ) * lda ] );
          }
          temp.setValue( x[ ioffx + jx - 1 ] );
          for ( i = j-1; i >= Math.max( 1, j - k ); i -- ) {
            x[ ioffx + ix - 1 ].minusEquals( ComplexMath.times( temp ,
              A[ ioffA + l + i - 1 + ( j - 1 ) * lda ] ) );
            ix -= incx;
          }
        }
        jx -= incx;
      }
    } else {
      jx = kx;
      for ( j = 1; j <= n; j ++ ) {
        kx += incx;
        if ( ! x[ ioffx + jx - 1 ].equals( zero ) ) {
          ix = kx;
          l = 1 - j;
          if ( nounit ) {
            x[ ioffx + jx - 1 ].divideEquals(
              A[ ioffA + ( j - 1 ) *lda ] );
          }
          temp.setValue( x[ ioffx + jx - 1 ] );
          for ( i = j + 1; i <= Math.min( n , j + k ); i ++ ) {
            x[ ioffx + ix - 1 ].minusEquals( ComplexMath.times( temp ,
              A[ ioffA + l + i - 1 + ( j - 1 ) * lda ] ) );
            ix += incx;
          }
        }
        jx += incx;
      }
    }
  } else {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      kplus1 = k + 1;
      jx = kx;
      for ( j = 1; j <= n; j ++ ) {
        temp.setValue( x[ ioffx + jx - 1 ] );
        ix = kx;
        l = kplus1 - j;
        if ( noconj ) {
          for ( i = Math.max( 1, j - k ); i <= j - 1; i ++ ) {
            temp.minusEquals( ComplexMath.times(
              A[ ioffA + l + i - 1 + ( j - 1 ) * lda ] ,
              x[ ioffx + ix - 1 ] ) );
            ix += incx;
          }
          if ( nounit ) {
            temp.divideEquals(
              A[ ioffA + kplus1 - 1 + ( j - 1 ) * lda ] );
          }
        } else {
          for ( i = Math.max( 1, j - k ); i <= j - 1; i ++ ) {
            temp.minusEquals( ComplexMath.times( ComplexMath.conj(
              A[ ioffA + l + i - 1 + ( j - 1 ) * lda ] ) ,
              x[ ioffx + ix - 1 ] ) );
            ix += incx;
          }
          if ( nounit ) {
            temp.divideEquals( ComplexMath.conj(
            A[ ioffA + kplus1 - 1 + ( j - 1 ) * lda ] ) );
          }
        }
        x[ ioffx + jx - 1 ].setValue( temp );
        jx += incx;
        if ( j > k ) kx += incx;
      }
    } else {
      kx += ( n - 1 ) * incx;
      jx = kx;
      for ( j = n; j >= 1; j -- ) {
        temp.setValue( x[ ioffx + jx - 1 ] );
        ix = kx;
        l = 1 - j;
        if ( noconj ) {
          for ( i = Math.min( n , j + k ); i >= j + 1; i -- ) {
            temp.minusEquals( ComplexMath.times(
              A[ ioffA + l + i - 1 + ( j - 1 ) * lda ] ,
              x[ ioffx + ix - 1 ] ) );
            ix -= incx;
          }
          if ( nounit ) {
            temp.divideEquals( A[ ioffA + ( j - 1 ) * lda ] );
          }
        } else {
          for ( i = Math.min( n, j + k ); i >= j + 1; i -- ) {
            temp.minusEquals( ComplexMath.times( ComplexMath.conj(
              A[ ioffA + l + i - 1 + ( j - 1 ) * lda ] ) ,
              x[ ioffx + ix - 1 ] ) );
            ix -= incx;
          }
          if ( nounit ) {
            temp.divideEquals( ComplexMath.conj(
                  A[ ioffA + ( j - 1 ) * lda ] ) );
          }
        }
        x[ ioffx + jx - 1 ].setValue( temp );
        jx -= incx;
        if ( n - j >= k ) kx -= incx;
      }
    }
  }
}
//*************************************************************************
//  x <-- A^{-1} x or x <-- A^{-T} x, A triangular packed
Blas2.dtpsv = function( uplo, trans, diag, n, AP, x, incx, ioffAP,
ioffx ) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U'
  && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( trans.charAt(0).toUpperCase() != 'N'
  && trans.charAt(0).toUpperCase() != 'T'
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info = 2;
  } else if ( diag.charAt(0).toUpperCase() != 'U'
  && diag.charAt(0).toUpperCase() != 'N' ) {
    info = 3;
  } else if (n < 0 ) info = 4;
  else if ( incx == 0 ) info = 7;
  if ( info != 0) Blas2.xerbla( 'dtpsv' , info );
  
  if ( n == 0 ) return;
  var nounit = ( diag.charAt(0).toUpperCase() == 'N' );
  var kx = 1;
  if ( incx <= 0 ) kx = - ( n - 1 ) * incx;
  else if ( incx != 1 ) kx = 1;
  var i = -1;
  var ix = -1;
  var j = -1;
  var jx = -1;
  var k = -1;
  var kk = -1;
  var temp = Number.POSITIVE_INFINITY;
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      kk = ( n * ( n + 1 ) ) / 2;
      jx = kx + ( n - 1 ) * incx;
      for ( j = n; j >= 1; j -- ) {
        if ( x[ ioffx + jx - 1 ] != 0. ) {
          if ( nounit ) x[ ioffx + jx - 1 ] /= AP[ ioffAP + kk - 1 ];
          temp = x[ ioffx + jx - 1 ];
          ix = jx;
          for ( k = kk - 1; k >= kk - j + 1; k -- ) {
            ix -= incx;
            x[ ioffx + ix - 1 ] -= temp * AP[ ioffAP + k - 1 ];
          }
        }
        jx -= incx;
        kk -= j;
      }
    } else {
      kk = 1;
      jx = kx;
      for ( j = 1; j <= n; j ++ ) {
        if ( x[ ioffx + jx - 1 ] != 0. ) {
          if ( nounit ) x[ ioffx + jx - 1 ] /= AP[ ioffAP + kk - 1 ];
          temp = x[ ioffx + jx - 1 ];
          ix = jx;
          for ( k = kk + 1; k <= kk + n - j; k ++ ) {
            ix += incx;
            x[ ioffx + ix - 1 ] -= temp * AP[ ioffAP + k - 1 ];
          }
        }
        jx += incx;
        kk += n - j + 1;
      }
    }
  } else {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      kk = 1;
      jx = kx;
      for ( j = 1; j <= n; j ++ ) {
        temp = x[ ioffx + jx - 1 ];
        ix = kx;
        for ( k = kk; k <= kk + j - 2; k ++ ) {
          temp -= AP[ ioffAP + k - 1 ] * x[ ioffx + ix - 1 ];
          ix += incx;
        }
        if ( nounit ) temp /= AP[ ioffAP + kk + j - 2 ];
        x[ ioffx + jx - 1 ] = temp;
        jx += incx;
        kk += j;
      }
    } else {
      kk = ( n * ( n + 1 ) ) / 2;
      kx += ( n - 1 ) * incx;
      jx = kx;
      for ( j = n; j >= 1; j -- ) {
        temp = x[ ioffx + jx - 1 ];
        ix = kx;
        for ( k = kk; k >= kk - n + j + 1; k -- ) {
          temp -= AP[ioffAP + k - 1 ] * x[ ioffx + ix - 1 ];
          ix -= incx;
        }
        if ( nounit ) temp /= AP[ ioffAP + kk - n + j - 1 ];
        x[ ioffx + jx - 1 ] = temp;
        jx -= incx;
        kk -= n - j + 1;
      }
    }
  }
}
Blas2.ztpsv = function( uplo, trans, diag, n, AP, x, incx, ioffAP,
ioffx ) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U'
  && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( trans.charAt(0).toUpperCase() != 'N'
  && trans.charAt(0).toUpperCase() != 'T'
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info = 2;
  } else if ( diag.charAt(0).toUpperCase() != 'U'
  && diag.charAt(0).toUpperCase() != 'N' ) {
    info = 3;
  } else if (n < 0 ) info = 4;
  else if ( incx == 0 ) info = 7;
  if ( info != 0) Blas2.xerbla( 'ztpsv' , info );
  
  if ( n == 0 ) return;
  var noconj = ( trans.charAt(0).toUpperCase() == 'T' );
  var nounit = ( diag.charAt(0).toUpperCase() == 'N' );
  var kx = 1;
  if ( incx <= 0 ) kx = - ( n - 1 ) * incx;
  else if ( incx != 1 ) kx = 1;
  var i = -1;
  var ix = -1;
  var j = -1;
  var jx = -1;
  var k = -1;
  var kk = -1;
  var temp = new Complex();
  var zero = new Complex( 0., 0. );
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      kk = ( n * ( n + 1 ) ) / 2;
      jx = kx + ( n - 1 ) * incx;
      for ( j = n; j >= 1; j -- ) {
        if ( ! x[ ioffx + jx - 1 ].equals( zero ) ) {
          if ( nounit ) {
            x[ ioffx + jx - 1 ].divideEquals( AP[ ioffAP + kk - 1 ] );
          }
          temp.setValue( x[ ioffx + jx - 1 ] );
          ix = jx;
          for ( k = kk - 1; k >= kk - j + 1; k -- ) {
            ix -= incx;
            x[ ioffx + ix - 1 ].minusEquals( ComplexMath.times( temp ,
              AP[ ioffAP + k - 1 ]));
          }
        }
        jx -= incx;
        kk -= j;
      }
    } else {
      kk = 1;
      jx = kx;
      for ( j = 1; j <= n; j ++ ) {
        if ( ! x[ ioffx + jx - 1 ].equals( zero ) ) {
          if ( nounit ) {
            x[ ioffx + jx - 1 ].divideEquals( AP[ ioffAP + kk - 1 ] );
          }
          temp.setValue( x[ ioffx + jx - 1 ] );
          ix = jx;
          for ( k = kk + 1; k <= kk + n - j; k ++ ) {
            ix += incx;
            x[ ioffx + ix - 1 ].minusEquals( ComplexMath.times( temp ,
              AP[ ioffAP + k - 1 ] ));
          }
        }
        jx += incx;
        kk += n - j + 1;
      }
    }
  } else {
    if ( uplo.charAt(0).toUpperCase() == 'U' ) {
      kk = 1;
      jx = kx;
      for ( j = 1; j <= n; j ++ ) {
        temp.setValue( x[ ioffx + jx - 1 ] );
        ix = kx;
        if ( noconj ) {
          for ( k = kk; k <= kk + j - 2; k ++ ) {
            temp.minusEquals( ComplexMath.times( AP[ ioffAP + k - 1 ] ,
              x[ ioffx + ix - 1 ] ) );
            ix += incx;
          }
          if ( nounit ) temp.divideEquals( AP[ ioffAP + kk + j - 2 ] );
        } else {
          for ( k = kk; k <= kk + j - 2; k ++ ) {
            temp.minusEquals( ComplexMath.times( ComplexMath.conj (
              AP[ ioffAP + k - 1 ] ) , x[ ioffx + ix - 1 ] ) );
            ix += incx;
          }
          if ( nounit ) {
            temp.divideEquals(
              ComplexMath.conj( AP[ ioffAP + kk + j - 2 ] ) );
          }
        }
        x[ ioffx + jx - 1 ].setValue( temp );
        jx += incx;
        kk += j;
      }
    } else {
      kk = ( n * ( n + 1 ) ) / 2;
      kx += ( n - 1 ) * incx;
      jx = kx;
      for ( j = n; j >= 1; j -- ) {
        temp.setValue( x[ ioffx + jx - 1 ] );
        ix = kx;
        if ( noconj ) {
          for ( k = kk; k >= kk - n + j + 1; k -- ) {
            temp.minusEquals( ComplexMath.times( AP[ ioffAP + k - 1 ] ,
              x[ ioffx + ix - 1 ] ) );
            ix -= incx;
          }
          if ( nounit ) {
            temp.divideEquals( AP[ ioffAP + kk - n + j - 1 ] );
          }
        } else {
          for ( k = kk; k >= kk - n + j + 1; k -- ) {
            temp.minusEquals( ComplexMath.times( ComplexMath.conj (
              AP[ ioffAP + k - 1 ] ) , x[ ioffx + ix - 1 ] ) );
            ix -= incx;
          }
          if ( nounit ) {
            temp.divideEquals(
              ComplexMath.conj ( AP[ ioffAP + kk - n + j - 1 ] ) );
          }
        }
        x[ ioffx + jx - 1 ].setValue( temp );
        jx -= incx;
        kk -= n - j + 1;
      }
    }
  }
}
//*************************************************************************
//  y <--- A x alpha + y beta or y <-- A^H x alpha + y beta
Blas2.dgemv = function( trans, m, n, alpha, A, lda, x, incx, beta, y,
incy, ioffA, ioffx, ioffy ) {
  var info = 0;
  if ( trans.charAt(0).toUpperCase() != 'N'
  && trans.charAt(0).toUpperCase() != 'T'
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info = 1;
  } else if ( m < 0 ) info = 2;
  else if ( n < 0 ) info = 3;
  else if ( lda < Math.max( 1, m ) ) info = 6;
  else if ( incx == 0 ) info = 8;
  else if ( incy == 0 ) info = 11;
  if ( info != 0) {
    Blas2.xerbla( 'dgemv' , info );
    return;
  }

  if ( m == 0 || n == 0 || ( alpha == 0. && beta == 1. ) ) return;
  var lenx = m;
  var leny = n;
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    lenx = n;
    leny = m;
  }
  var kx = ( incx > 0 ? 1 : 1 - ( lenx - 1 ) * incx );
  var ky = ( incy > 0 ? 1 : 1 - ( leny - 1 ) * incy );

  var i = -1;
  var iy = -1;
  if ( beta != 1. ) {
    iy = ky;
    if ( beta == 0. ) {
      for ( i = 1; i <= leny; i ++ ) {
        y[ ioffy + iy - 1 ] = 0.;
        iy += incy;
      }
    } else {
      for ( i = 1; i <= leny; i ++ ) {
        y[ ioffy + iy - 1 ] *= beta;
        iy += incy;
      }
    }
  }
  if ( alpha == 0. ) return
  var temp = Number.POSITIVE_INFINITY;
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    var j = -1;
    var jx = kx;
    for ( j = 1; j <= n; j ++ ) {
      if ( x[ ioffx + jx - 1 ] != 0. ) {
        temp = alpha * x [ ioffx + jx - 1 ];
        iy = ky;
        for ( i = 1; i <= m; i ++ ) {
          y[ ioffy + iy - 1 ] +=
            temp * A[ ioffA + i - 1 + ( j - 1 ) * lda ];
          iy += incy;
        }
      }
      jx += incx;
    }
  } else {
    var jy = ky;
    for ( j = 1; j <= n; j ++ ) {
      temp = 0.;
      var ix = kx;
      for ( i = 1; i <= m; i ++ ) {
        temp += A[ ioffA + i - 1 + ( j - 1 ) * lda ]
          * x[ ioffx + ix - 1 ];
        ix += incx;
      }
      y[ ioffy + jy - 1 ] += alpha * temp;
      jy += incy;
    }
  }
}
Blas2.zgemv = function( trans, m, n, alpha, A, lda, x, incx, beta, y,
incy, ioffA, ioffx, ioffy ) {
  var info = 0;
  if ( trans.charAt(0).toUpperCase() != 'N'
  && trans.charAt(0).toUpperCase() != 'T'
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info = 1;
  } else if ( m < 0 ) info = 2;
  else if ( n < 0 ) info = 3;
  else if ( lda < Math.max( 1, m ) ) info = 6;
  else if ( incx == 0 ) info = 8;
  else if ( incy == 0 ) info = 11;
  if ( info != 0) {
    Blas2.xerbla( 'zgemv' , info );
    return;
  }

  var zero = new Complex( 0., 0. );
  var one = new Complex( 1., 0. );
  if ( m == 0 || n == 0 ||
  ( alpha.equals(zero) && beta.equals(one) ) ) {
    return;
  }
  var noconj = ( trans.charAt(0).toUpperCase() == 'T' );
  var lenx = m;
  var leny = n;
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    lenx = n;
    leny = m;
  }
  var kx = ( incx > 0 ? 1 : 1 - ( lenx - 1 ) * incx );
  var ky = ( incy > 0 ? 1 : 1 - ( leny - 1 ) * incy );

  var i = -1;
  var iy = -1;
  if ( ! beta.equals( one ) ) {
    iy = ky;
    if ( beta.equals( zero ) ) {
      for ( i = 1; i <= leny; i ++ ) {
        y[ ioffy + iy - 1 ].setValue( zero );
        iy += incy;
      }
    } else {
      for ( i = 1; i <= leny; i ++ ) {
        y[ ioffy + iy - 1 ].timesEquals( beta );
        iy += incy;
      }
    }
  }
  if ( alpha.equals ( zero ) ) return;
  var temp = new Complex();
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    var j = -1;
    var jx = kx;
    for ( j = 1; j <= n; j ++ ) {
      if ( x[ ioffx + jx - 1 ] != 0. ) {
        temp.setValue( ComplexMath.times( alpha , x [ ioffx + jx - 1 ] ) );
        iy = ky;
        for ( i = 1; i <= m; i ++ ) {
          y[ ioffy + iy - 1 ].plusEquals( ComplexMath.times( temp ,
            A[ ioffA + i - 1 + ( j - 1 ) * lda ] ) );
          iy += incy;
        }
      }
      jx += incx;
    }
  } else {
    var jy = ky;
    for ( j = 1; j <= n; j ++ ) {
      temp.setValue( zero );
      var ix = kx;
      if ( noconj ) {
        for ( i = 1; i <= m; i ++ ) {
          temp.plusEquals( ComplexMath.times(
            A[ ioffA + i - 1 + ( j - 1 ) * lda ] ,
            x[ ioffx + ix - 1 ] ) );
          ix += incx;
        }
      } else {
        for ( i = 1; i <= m; i ++ ) {
          temp.plusEquals( ComplexMath.times( ComplexMath.conj(
            A[ ioffA + i - 1 + ( j - 1 ) * lda ] ) ,
            x[ ioffx + ix - 1 ] ) );
          ix += incx;
        }
      }
      y[ ioffy + jy - 1 ].plusEquals(
        ComplexMath.times( alpha , temp ) );
      jy += incy;
    }
  }
}
//*************************************************************************
//  y <-- A x alpha + y beta or y <-- A^H x alpha + y beta
Blas2.dgbmv = function( trans, m, n, kl, ku, alpha, A, lda, x, incx, beta,
y, incy, ioffA, ioffx, ioffy ) {
  var info = 0;
  if ( trans.charAt(0).toUpperCase() != 'N'
  && trans.charAt(0).toUpperCase() != 'T'
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info = 1;
  } else if ( m < 0 ) info = 2;
  else if ( n < 0 ) info = 3;
  else if ( kl < 0 ) info = 4;
  else if ( ku < 0 ) info = 5;
  else if ( lda < ( kl + ku + 1 ) ) info = 8;
  else if ( incx == 0 ) info = 10;
  else if ( incy == 0 ) info = 13;
  if ( info != 0 ) {
    Blas2.xerbla( 'dgbmv', info );
    return;
  }

  if ( m == 0 || n == 0 || ( alpha == 0. && beta == 1. ) ) return;
  var lenx = m;
  var leny = n;
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    lenx = n;
    leny = m;
  }
  var kx = ( incx > 0 ? 1 : 1 - ( lenx - 1 ) * incx );
  var ky = ( incy > 0 ? 1 : 1 - ( leny - 1 ) * incy );

  var i = -1;
  var iy = -1;
  if ( beta != 1. ) {
    iy = ky;
    if ( beta == 0. ) {
      for ( i = 1; i <= leny; i ++ ) {
        y[ ioffy + iy - 1 ] = 0.;
        iy += incy;
      }
    } else {
      for ( i = 1; i <= leny; i ++ ) {
        y[ ioffy + iy - 1 ] *= beta;
        iy += incy;
      }
    }
  }
  if ( alpha == 0. ) return;
  var kup1 = ku + 1;
  var temp = Number.POSITIVE_INFINITY;
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    var jx = kx;
    for ( var j = 1; j <= n; j ++ ) {
      if ( x[ ioffx + jx - 1 ] != 0. ) {
        temp = alpha * x[ ioffx + jx - 1 ];
        iy = ky;
        var k = kup1 - j;
        for ( i = Math.max( 1 , j - ku );
        i <= Math.min( m , j + kl ); i ++ ) {
          y[ ioffy + iy - 1 ] +=
            temp * A[ ioffA + k + i - 1 + ( j - 1 ) * lda ];
          iy += incy;
        }
      }
      jx += incx;
      if ( j > ku ) ky += incy;
    }
  } else {
    var jy = ky;
    for ( j = 1; j <= n; j ++ ) {
      temp = 0.;
      var ix = kx;
      k = kup1 - j;
      for ( i = Math.max( 1, j - ku );
      i <= Math.min( m , j + kl ); i ++ ) {
        temp += A[ ioffA + k + i - 1 + ( j - 1 ) * lda ]
          * x[ ioffx + ix - 1 ];
        ix += incx;
      }
      y[ ioffy + jy - 1 ] += alpha * temp;
      jy += incy;
      if ( j > ku ) kx += incx;
    }
  }
}
Blas2.zgbmv = function( trans, m, n, kl, ku, alpha, A, lda, x, incx, beta,
y, incy, ioffA, ioffx, ioffy) {
  var info = 0;
  if ( trans.charAt(0).toUpperCase() != 'N'
  && trans.charAt(0).toUpperCase() != 'T'
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info = 1;
  } else if ( m < 0 ) info = 2;
  else if ( n < 0 ) info = 3;
  else if ( kl < 0 ) info = 4;
  else if ( ku < 0 ) info = 5;
  else if ( lda < ( kl + ku + 1 ) ) info = 8;
  else if ( incx == 0 ) info = 10;
  else if ( incy == 0 ) info = 13;
  if ( info != 0 ) {
    Blas2.xerbla( 'zgbmv', info );
    return;
  }

  var zero = new Complex( 0., 0. );
  var one = new Complex( 1., 0. );
  if ( m == 0 || n == 0 ||
  ( alpha.equals( zero ) && beta.equals( one ) ) ) {
    return;
  }
  var noconj = ( trans.charAt(0).toUpperCase()=='T' );
  var lenx = m;
  var leny = n;
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    lenx = n;
    leny = m;
  }
  var kx = ( incx > 0 ? 1 : 1 - ( lenx - 1 ) * incx );
  var ky = ( incy > 0 ? 1 : 1 - ( leny - 1 ) * incy );

  var i = -1;
  var iy = -1;
  if ( ! beta.equals( one ) ) {
    iy = ky;
    if ( beta.equals( zero ) ) {
      for ( i = 1; i <= leny; i ++ ) {
        y[ ioffy + iy - 1 ].setValue( zero );
        iy += incy;
      }
    } else {
      for ( i = 1; i <= leny; i ++ ) {
        y[ ioffy + iy - 1 ].timesEquals( beta );
        iy += incy;
      }
    }
  }
  if ( alpha.equals( zero ) ) return;
  var kup1 = ku + 1;
  var temp = new Complex();
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    var jx = kx;
    for ( var j = 1; j <= n; j ++ ) {
      if ( ! x[ ioffx + jx - 1 ].equals( zero ) ) {
        temp.setValue( ComplexMath.times( alpha , x[ ioffx + jx - 1 ] ) );
        iy = ky;
        var k = kup1 - j;
        for ( i = Math.max( 1 , j - ku );
        i <= Math.min( m , j + kl ); i ++ ) {
          y[ ioffy + iy - 1 ].plusEquals( ComplexMath.times(
            temp , A[ ioffA + k + i - 1 + ( j - 1 ) * lda ] ) );
          iy += incy;
        }
      }
      jx += incx;
      if ( j > ku ) ky += incy;
    }
  } else {
    var jy = ky;
    for ( j = 1; j <= n; j ++ ) {
      temp.setValue( zero );
      var ix = kx;
      k = kup1 - j;
      if ( noconj ) {
        for ( i = Math.max( 1, j - ku );
        i <= Math.min( m , j + kl ); i ++ ) {
          temp.plusEquals( ComplexMath.times(
            A[ ioffA + k + i - 1 + ( j - 1 ) *lda ] ,
            x[ ioffx + ix - 1 ] ) );
          ix += incx;
        }
      } else {
        for ( i = Math.max( 1, j - ku );
        i <= Math.min( m, j + kl ); i ++ ) {
          temp.plusEquals( ComplexMath.times( ComplexMath.conj(
            A[ ioffA + k + i - 1 + ( j - 1 ) *lda ] ) ,
            x[ ioffx + ix - 1 ] ) );
          ix += incx;
        }
      }
      y[ ioffy + jy - 1 ].plusEquals(
        ComplexMath.times( alpha , temp ) );
      jy += incy;
      if ( j > ku ) kx += incx;
    }
  }
}
//*************************************************************************
//  y <-- A x alpha + y beta, A symmetric
Blas2.dsymv = function( uplo, n, alpha, A, lda, x, incx, beta, y, incy,
ioffA, ioffx, ioffy ) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U'
  && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( n < 0 ) info = 2;
  else if ( lda < Math.max( 1, n ) ) info = 5;
  else if ( incx == 0 ) info = 7;
  else if ( incy == 0 ) info = 10;
  if ( info != 0 ) {
    Blas2.xerbla( 'dsymv', info );
    return;
  }

  if ( n == 0 || ( alpha == 0. && beta == 1. ) ) return;
  var kx = ( incx > 0 ? 1 : 1 - ( n - 1 ) * incx );
  var ky = ( incy > 0 ? 1 : 1 - ( n - 1 ) * incy );
  var i = -1;
  var iy = -1;
  var j = -1;
  var jx = -1;
  var jy = -1;
  var temp1 = Number.POSITIVE_INFINITY;
  var temp2 = Number.POSITIVE_INFINITY;
  if ( beta != 1. ) {
    iy = ky;
    if ( beta == 0. ) {
      for ( i = 1; i <= n; i ++ ) {
        y[ ioffy + iy - 1 ] = 0.;
        iy += incy;
      }
    } else {
      for ( i = 1; i <= n; i ++ ) {
        y[ ioffy + iy - 1 ] *= beta;
        iy += incy;
      }
    }
  }
  if ( alpha == 0. ) return;
  var ix = -1;
  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    jx = kx;
    jy = ky;
    for ( j = 1; j <= n; j ++ ) {
      temp1 = alpha * x[ ioffx + jx - 1 ];
      temp2 = 0.;
      ix = kx;
      iy = ky;
      for ( i = 1; i <= j-1; i ++ ) {
        y[ ioffy + iy - 1 ] +=
          temp1 * A[ ioffA + i - 1 + ( j - 1 ) * lda ];
        temp2 += A[ ioffA + i - 1 + ( j - 1 ) * lda ]
          * x[ ioffx + ix - 1 ];
        ix += incx;
        iy += incy;
      }
      y[ ioffy + jy - 1 ] +=
        temp1 * A[ ioffA + j - 1 + ( j - 1 ) * lda ] + alpha * temp2;
      jx += incx;
      jy += incy;
    }
  } else {
    jx = kx;
    jy = ky;
    for ( j = 1; j <= n; j ++ ) {
      temp1 = alpha * x[ ioffx + jx - 1 ];
      temp2 = 0.;
      y[ ioffy + jy - 1 ] +=
        temp1 * A[ ioffA + j - 1 + ( j - 1 ) * lda ];
      ix = jx;
      iy = jy;
      for ( i = j+1; i <= n; i ++ ) {
        ix += incx;
        iy += incy;
        y[ ioffy + iy - 1 ] +=
          temp1 * A[ ioffA + i - 1 + ( j - 1 ) * lda ];
        temp2 += A[ ioffA + i - 1 + ( j - 1 ) * lda ]
          * x[ ioffx + ix - 1 ];
      }
      y[ ioffy + jy - 1 ] += alpha * temp2;
      jx += incx;
      jy += incy;
    }
  }
}
Blas2.zhemv = function( uplo, n, alpha, A, lda, x, incx, beta, y, incy,
ioffA, ioffx, ioffy ) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U' 
  && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( n < 0 ) info = 2;
  else if ( lda < Math.max( 1, n ) ) info = 5;
  else if ( incx == 0 ) info = 7;
  else if ( incy == 0 ) info = 10;
  if ( info != 0 ) {
    Blas2.xerbla( 'zhemv', info );
    return;
  }

  var zero = new Complex( 0., 0. );
  var one = new Complex( 1., 0. );
  if ( n == 0 || ( alpha.equals( zero ) && beta.equals( one ) ) ) {
    return;
  }
  var kx = ( incx > 0 ? 1 : 1 - ( n - 1 ) * incx );
  var ky = ( incy > 0 ? 1 : 1 - ( n - 1 ) * incy );
  var i = -1;
  var iy = -1;
  var j = -1;
  var jx = -1;
  var jy = -1;
  var temp1 = new Complex();
  var temp2 = new Complex();
  if ( ! beta.equals( one ) ) {
    iy = ky;
    if ( beta.equals( zero ) ) {
      for ( i = 1; i <= n; i ++ ) {
        y[ ioffy + iy - 1 ].setValue( zero );
        iy += incy;
      }
    } else {
      for ( i = 1; i <= n; i ++ ) {
        y[ ioffy + iy - 1 ].timesEquals( beta );
        iy += incy;
      }
    }
  }
  if ( alpha.equals( zero ) ) return;
  var ix = -1;
  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    jx = kx;
    jy = ky;
    for ( j = 1; j <= n; j ++ ) {
      temp1.setValue( ComplexMath.times( alpha , x[ ioffx + jx - 1 ] ) );
      temp2.setValue( zero );
      ix = kx;
      iy = ky;
      for ( i = 1; i <= j-1; i ++ ) {
        y[ ioffy + iy - 1 ].plusEquals( ComplexMath.times( temp1 ,
          A[ ioffA + i - 1 + ( j - 1 ) * lda ] ) );
        temp2.plusEquals( ComplexMath.times( 
          ComplexMath.conj( A[ ioffA + i - 1 + ( j - 1 ) * lda ] ) ,
          x[ ioffx + ix - 1 ] ) );
        ix += incx;
        iy += incy;
      }
      y[ ioffy + jy - 1 ].plusEquals( ComplexMath.plus(
        ComplexMath.timesNumber( temp1 ,
                         A[ ioffA + j - 1 + ( j - 1 ) * lda ].getReal() ) ,
        ComplexMath.times( alpha , temp2 ) ) );
      jx += incx;
      jy += incy;
    }
  } else {
    jx = kx;
    jy = ky;
    for ( j = 1; j <= n; j ++ ) {
      temp1.setValue( ComplexMath.times( alpha , x[ ioffx + jx - 1 ] ) );
      temp2.setValue( zero );
      y[ ioffy + jy - 1 ].plusEquals( ComplexMath.timesNumber( temp1 ,
                       A[ ioffA + j - 1 + ( j - 1 ) * lda ].getReal() ) );
      ix = jx;
      iy = jy;
      for ( i = j+1; i <= n; i ++ ) {
        ix += incx;
        iy += incy;
        y[ ioffy + iy - 1 ].plusEquals( ComplexMath.times(
          temp1 , A[ ioffA + i - 1 + ( j - 1 ) * lda ] ) );
        temp2.plusEquals( ComplexMath.times(
          ComplexMath.conj( A[ ioffA + i - 1 + ( j - 1 ) * lda ] ) ,
                              x[ ioffx + ix - 1 ] ) );
      }
      y[ ioffy + jy - 1 ].plusEquals(
        ComplexMath.times( alpha , temp2 ) );
      jx += incx;
      jy += incy;
    }
  }
}
//*************************************************************************
//  y <-- A x alpha + y beta , A symmetric banded
Blas2.dsbmv = function( uplo, n, k, alpha, A, lda, x, incx, beta, y, incy,
ioffA, ioffx, ioffy ) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U' 
  && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( n < 0 ) info = 2;
  else if ( k < 0 ) info = 3;
  else if ( lda < k + 1 ) info = 6;
  else if ( incx == 0 ) info = 8;
  else if ( incy == 0 ) info = 11;
  if ( info != 0 ) {
    Blas2.xerbla( 'dsbmv', info );
    return;
  }

  if ( n == 0 || ( alpha == 0. && beta == 1. ) ) return;
  var kx = ( incx > 0 ? 1 : 1 - ( n - 1 ) * incx );
  var ky = ( incy > 0 ? 1 : 1 - ( n - 1 ) * incy );
  
  var iy = -1;
  if ( beta != 1. ) {
    iy = ky;
    if ( beta == 0. ) {
      for ( i = 1; i <= n; i ++ ) {
        y[ ioffy + iy - 1 ] = 0.;
        iy += incy;
      }
    } else {
      for ( i = 1; i <= n; i ++ ) {
        y[ ioffy + iy - 1 ] *= beta;
        iy += incy;
      }
    }
  }
  var i = -1;
  var ix = -1;
  var j = -1;
  var jx = -1;
  var jy = -1;
  var l = -1;
  var temp1 = Number.POSITIVE_INFINITY;
  var temp2 = Number.POSITIVE_INFINITY;
  if ( alpha == 0. ) return;
  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    var kplus1 = k + 1;
    jx = kx;
    jy = ky;
    for ( j = 1; j <= n; j ++ ) {
      temp1 = alpha * x[ ioffx + jx - 1 ];
      temp2 = 0.;
      ix = kx;
      iy = ky;
      l = kplus1 - j;
      for ( i = Math.max( 1 , j - k ); i <= j - 1; i ++ ) {
        y[ ioffy + iy - 1 ] +=
          temp1 * A[ ioffA + l + i - 1 + ( j - 1 ) * lda ];
        temp2 += A[ ioffA + l + i - 1 + ( j - 1 ) * lda ]
          * x[ ioffx + ix - 1 ];
        ix += incx;
        iy += incy;
      }
      y[ ioffy + jy - 1 ] +=
        temp1 * A[ ioffA + kplus1 - 1 + ( j - 1 ) * lda ]
        + alpha * temp2;
      jx += incx;
      jy += incy;
      if (j > k ) {
        kx += incx;
        ky += incy;
      }
    }
  } else {
    jx = kx;
    jy = ky;
    for ( j = 1; j <= n; j ++ ) {
      temp1 = alpha * x[ ioffx + jx - 1 ];
      temp2 = 0.;
      y[ ioffy + jy - 1 ] += temp1 * A[ ioffA + ( j - 1 ) * lda ];
      l = 1 - j;
      ix = jx;
      iy = jy;
      for ( i = j + 1; i <= Math.min( n , j + k ); i ++ ) {
        ix += incx;
        iy += incy;
        y[ ioffy + iy - 1 ] +=
          temp1 * A[ ioffA + l + i - 1 + ( j - 1 ) * lda ];
        temp2 += A[ ioffA + l + i - 1 + ( j - 1 ) * lda ]
          * x[ ioffx + ix - 1 ];
      }
      y[ ioffy + jy - 1 ] += alpha * temp2;
      jx += incx;
      jy += incy;
    }
  }
}
Blas2.zhbmv = function( uplo, n, k, alpha, A, lda, x, incx, beta, y, incy,
ioffA, ioffx, ioffy ) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U' 
  && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( n < 0 ) info = 2;
  else if ( k < 0 ) info = 3;
  else if ( lda < k + 1 ) info = 6;
  else if ( incx == 0 ) info = 8;
  else if ( incy == 0 ) info = 11;
  if ( info != 0 ) {
    Blas2.xerbla( 'zhbmv', info );
    return;
  }

  var zero = new Complex( 0., 0. );
  var one = new Complex( 1., 0. );
  if ( n == 0 || ( alpha.equals( zero ) && beta.equals(one) ) ) {
    return;
  }
  var kx = ( incx > 0 ? 1 : 1 - ( n - 1 ) * incx );
  var ky = ( incy > 0 ? 1 : 1 - ( n - 1 ) * incy );

  var iy = -1;
  if ( ! beta.equals( one ) ) {
    iy = ky;
    if ( beta.equals(zero) ) {
      for ( i = 1; i <= n; i ++ ) {
        y[ ioffy + iy - 1 ].setValue( zero );
        iy += incy;
      }
    } else {
      for ( i = 1; i <= n; i ++ ) {
        y[ ioffy + iy - 1 ].timesEquals( beta );
        iy += incy;
      }
    }
  }
  var i = -1;
  var ix = -1;
  var j = -1;
  var jx = -1;
  var jy = -1;
  var l = -1;
  var temp1 = new Complex();
  var temp2 = new Complex();
  if ( alpha.equals( zero ) ) return;
  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    var kplus1 = k + 1;
    jx = kx;
    jy = ky;
    for ( j = 1; j <= n; j ++ ) {
      temp1.setValue( ComplexMath.times( alpha , x[ ioffx + jx - 1 ] ) );
      temp2.setValue( zero );
      ix = kx;
      iy = ky;
      l = kplus1 - j;
      for ( i = Math.max( 1 , j - k ); i <= j - 1; i ++ ) {
        y[ ioffy + iy - 1 ].plusEquals( ComplexMath.times( temp1 ,
          A[ ioffA + l + i - 1 + ( j - 1 ) * lda ] ) );
        temp2.plusEquals( ComplexMath.times( ComplexMath.conj(
          A[ ioffA + l + i - 1 + ( j - 1 ) * lda ] ) ,
          x[ ioffx + ix - 1 ] ) );
        ix += incx;
        iy += incy;
      }
      y[ ioffy + jy - 1 ].plusEquals( ComplexMath.plus(
        ComplexMath.timesNumber( temp1 ,
        A[ ioffA + kplus1 - 1 + ( j - 1 ) * lda ].getReal()) ,
        ComplexMath.times( alpha , temp2 ) ) );
      jx += incx;
      jy += incy;
      if (j > k ) {
        kx += incx;
        ky += incy;
      }
    }
  } else {
    jx = kx;
    jy = ky;
    for ( j = 1; j <= n; j ++ ) {
      temp1.setValue( ComplexMath.times( alpha , x[ ioffx + jx - 1 ] ) );
      temp2.setValue( zero );
      y[ ioffy + jy - 1 ].plusEquals(
        ComplexMath.timesNumber( temp1 ,
                               A[ ioffA + ( j - 1 ) * lda ].getReal() ) );
      l = 1 - j;
      ix = jx;
      iy = jy;
      for ( i = j + 1; i <= Math.min( n, j + k ); i ++ ) {
        ix += incx;
        iy += incy;
        y[ ioffy + iy - 1 ].plusEquals( ComplexMath.times( temp1 ,
          A[ ioffA + l + i - 1 + ( j - 1 ) * lda ] ) );
        temp2.plusEquals( ComplexMath.times( ComplexMath.conj( 
          A[ ioffA + l + i - 1 + ( j - 1 ) * lda ] ) ,
          x[ ioffx + ix - 1 ] ) );
      }
      y[ ioffy + jy - 1 ].plusEquals(
        ComplexMath.times( alpha , temp2 ) );
      jx += incx;
      jy += incy;
    }
  }
}
//*************************************************************************
//  y <-- A x alpha + y beta, A symmetric packed storage
Blas2.dspmv = function( uplo, n, alpha, AP, x, incx, beta, y, incy,
ioffAP, ioffx, ioffy ) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U' 
  && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( n < 0 ) info = 2;
  else if ( incx == 0 ) info = 6;
  else if ( incy == 0 ) info = 9;
  if ( info != 0 ) {
    Blas2.xerbla( 'dspmv', info );
    return;
  }

  if ( n == 0 || ( alpha == 0. && beta == 1. ) ) return;
  var kx = ( incx > 0 ? 1 : 1 - ( n - 1 ) * incx );
  var ky = ( incy > 0 ? 1 : 1 - ( n - 1 ) * incy );

  var i = -1;
  var j = -1;
  var k = -1;
  var temp1 = Number.POSITIVE_INFINITY;
  var temp2 = Number.POSITIVE_INFINITY;
  if ( beta != 1. ) {
    var iy = ky;
    if ( beta == 0. ) {
      for ( i = 1; i <= n; i ++ ) {
        y[ ioffy + iy - 1 ] = 0.;
        iy += incy;
      }
    } else {
      for ( i = 1; i <= n; i ++ ) {
        y[ ioffy + iy - 1 ] *= beta;
        iy += incy;
      }
    }
  }
  if ( alpha == 0. ) return;
  var ix = -1;
  var kk = 1;
  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    var jx = kx;
    var jy = ky;
    for ( j = 1; j <= n; j ++ ) {
      temp1 = alpha * x[ ioffx + jx - 1 ];
      temp2 = 0.;
      ix = kx;
      iy = ky;
      for ( k = kk; k <= kk + j - 2; k ++ ) {
        y[ ioffy + iy - 1 ] += temp1 * AP[ ioffAP + k - 1 ];
        temp2 += AP[ ioffAP + k - 1 ] * x[ ioffx + ix - 1 ];
        ix += incx;
        iy += incy;
      }
      y[ ioffy + jy - 1 ] += temp1 * AP[ ioffAP + kk + j - 2 ]
                      + alpha * temp2;
      jx += incx;
      jy += incy;
      kk += j;
    }
  } else {
    jx = kx;
    jy = ky;
    for ( j = 1; j <= n; j ++ ) {
      temp1 = alpha * x[ ioffx + jx - 1 ];
      temp2 = 0.;
      y[ ioffy + jy - 1 ] += temp1 * AP[ ioffAP + kk - 1 ];
      ix = jx;
      iy = jy;
      for ( k = kk + 1; k <= kk + n - j; k ++ ) {
        ix += incx;
        iy += incy;
        y[ ioffy + iy - 1 ] += temp1 * AP[ ioffAP + k - 1 ];
        temp2 += AP[ ioffAP + k - 1 ] * x[ ioffx + ix - 1 ];
      }
      y[ ioffy + jy - 1 ] += alpha * temp2;
      jx += incx;
      jy += incy;
      kk += n - j + 1;
    }
  }
}
Blas2.zhpmv = function( uplo, n, alpha, AP, x, incx, beta, y, incy,
ioffAP, ioffx, ioffy ) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U' 
  && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( n < 0 ) info = 2;
  else if ( incx == 0 ) info = 6;
  else if ( incy == 0 ) info = 9;
  if ( info != 0 ) {
    Blas2.xerbla( 'zhpmv', info );
    return;
  }

  var zero = new Complex( 0., 0. );
  var one = new Complex( 1., 0. );
  if ( n == 0 || ( alpha.equals( zero ) && beta.equals(one) ) ) {
    return;
  }
  var kx = ( incx > 0 ? 1 : 1 - ( n - 1 ) * incx );
  var ky = ( incy > 0 ? 1 : 1 - ( n - 1 ) * incy );

  var i = -1;
  var j = -1;
  var k = -1;
  var temp1 = new Complex();
  var temp2 = new Complex();
  if ( ! beta.equals(one) ) {
    var iy = ky;
    if ( beta.equals( zero ) ) {
      for ( i = 1; i <= n; i ++ ) {
        y[ ioffy + iy - 1 ].setValue( zero );
        iy += incy;
      }
    } else {
      for ( i = 1; i <= n; i ++ ) {
        y[ ioffy + iy - 1 ].timesEquals( beta );
        iy += incy;
      }
    }
  }
  if ( alpha.equals(zero) ) return;
  var ix = -1;
  var kk = 1;
  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    var jx = kx;
    var jy = ky;
    for ( j = 1; j <= n; j ++ ) {
      temp1.setValue( ComplexMath.times( alpha , x[ ioffx + jx - 1 ] ) );
      temp2.setValue( zero );
      ix = kx;
      iy = ky;
      for ( k = kk; k <= kk + j - 2; k ++ ) {
        y[ ioffy + iy - 1 ].plusEquals(
          ComplexMath.times( temp1 , AP[ ioffAP + k - 1 ] ) );
        temp2.plusEquals( ComplexMath.times( ComplexMath.conj(
          AP[ ioffAP + k - 1 ] ) , x[ ioffx + ix - 1 ] ) );
        ix += incx;
        iy += incy;
      }
      y[ ioffy + jy - 1 ].plusEquals( ComplexMath.plus(
        ComplexMath.timesNumber( temp1 ,
                                 AP[ ioffAP + kk + j - 2 ].getReal() ) ,
        ComplexMath.times( alpha , temp2 ) ) );
      jx += incx;
      jy += incy;
      kk += j;
    }
  } else {
    jx = kx;
    jy = ky;
    for ( j = 1; j <= n; j ++ ) {
      temp1.setValue( ComplexMath.times( alpha , x[ ioffx + jx - 1 ] ) );
      temp2.setValue( zero );
      y[ ioffy + jy - 1 ].plusEquals(
        ComplexMath.timesNumber( temp1 ,
                                 AP[ ioffAP + kk - 1 ].getReal() ) );
      ix = jx;
      iy = jy;
      for ( k = kk + 1; k <= kk + n - j; k ++ ) {
        ix += incx;
        iy += incy;
        y[ ioffy + iy - 1 ].plusEquals(
          ComplexMath.times( temp1 , AP[ ioffAP + k - 1 ] ) );
        temp2.plusEquals( ComplexMath.times(
          ComplexMath.conj( AP[ ioffAP + k - 1 ] ) ,
          x[ ioffx + ix - 1 ] ) );
      }
      y[ ioffy + jy - 1 ].plusEquals(
        ComplexMath.times( alpha , temp2 ) );
      jx += incx;
      jy += incy;
      kk += n - j + 1;
    }
  }
}
//*************************************************************************
//  A += x alpha y^T or A += x alpha y^H
Blas2.dger = function( m, n, alpha, x, incx, y, incy, A, lda, ioffx,
ioffy, ioffA ) {
  var info = 0;
  if ( m < 0 ) info = 1;
  else if ( n < 0 ) info = 2;
  else if ( incx == 0 ) info = 5;
  else if ( incy == 0 ) info = 7;
  else if ( lda < Math.max( 1 , m ) ) info = 9;
  if ( info != 0 ) {
    Blas2.xerbla( 'dger', info );
    return;
  }

  if ( m == 0 || n == 0 || alpha == 0. ) return;
  var jy = ( incy > 0 ? 1 : 1 - ( n - 1 ) * incy );
  var kx = ( incx > 0 ? 1 : 1 - ( m - 1 ) * incx );
  for ( var j = 1; j <= n; j ++ ) {
    if ( y[ ioffy + jy - 1 ] != 0. ) {
      var temp = alpha * y[ ioffy + jy - 1 ];
      var ix = kx;
      for ( var i = 1; i <= m; i ++ ) {
        A[ ioffA + i - 1 +  ( j - 1 ) * lda ] +=
          x[ ioffx + ix - 1 ] * temp;
        ix += incx;
      }
    }
    jy += incy;
  }
}
Blas2.zgerc = function( m, n, alpha, x, incx, y, incy, A, lda, ioffx,
ioffy, ioffA ) {
  var info = 0;
  if ( m < 0 ) info = 1;
  else if ( n < 0 ) info = 2;
  else if ( incx == 0 ) info = 5;
  else if ( incy == 0 ) info = 7;
  else if ( lda < Math.max( 1 , m ) ) info = 9;
  if ( info != 0 ) {
    Blas2.xerbla( 'zgerc', info );
    return;
  }

  var zero = new Complex( 0., 0. );
  if ( m == 0 || n == 0 || alpha.equals( zero ) ) return;
  var jy = ( incy > 0 ? 1 : 1 - ( n - 1 ) * incy );
  var kx = ( incx > 0 ? 1 : 1 - ( m - 1 ) * incx );
  for ( var j = 1; j <= n; j ++ ) {
    if ( ! y[ ioffy + jy - 1 ].equals( zero ) ) {
      var temp
        = ComplexMath.times( alpha ,
                             ComplexMath.conj( y[ ioffy + jy - 1 ] ) );
      var ix = kx;
      for ( var i = 1; i <= m; i ++ ) {
        A[ ioffA + i - 1 + ( j - 1 ) * lda ].plusEquals( 
          ComplexMath.times( x[ ioffx + ix - 1 ] , temp ) );
        ix += incx;
      }
    }
    jy += incy;
  }
}
Blas2.zgeru = function( m, n, alpha, x, incx, y, incy, A, lda, ioffx,
ioffy, ioffA ) {
  var info = 0;
  if ( m < 0 ) info = 1;
  else if ( n < 0 ) info = 2;
  else if ( incx == 0 ) info = 5;
  else if ( incy == 0 ) info = 7;
  else if ( lda < Math.max( 1 , m ) ) info = 9;
  if ( info != 0 ) {
    Blas2.xerbla( 'zgeru', info );
    return;
  }

  var zero = new Complex( 0., 0. );
  if ( m == 0 || n == 0 || alpha.equals( zero ) ) return;
  var jy = ( incy > 0 ? 1 : 1 - ( n - 1 ) * incy );
  var kx = ( incx > 0 ? 1 : 1 - ( m - 1 ) * incx );
  for ( var j = 1; j <= n; j ++ ) {
    if ( ! y[ ioffy + jy - 1 ].equals( zero ) ) {
      var temp =
        ComplexMath.times( alpha , y[ ioffy + jy - 1 ] );
      var ix = kx;
      for ( var i = 1; i <= m; i ++ ) {
        A[ ioffA + i - 1 + ( j - 1 ) * lda ].plusEquals( 
          ComplexMath.times( x[ ioffx + ix - 1 ] , temp ) );
        ix += incx;
      }
    }
    jy += incy;
  }
}
//*************************************************************************
//  A += x alpha x^H
Blas2.dsyr = function( uplo, n, alpha, x, incx, A, lda, ioffx, ioffA ) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U' 
  && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( n < 0 ) info = 2;
  else if ( incx == 0 ) info = 5;
  else if ( lda < Math.max( 1 , n ) ) info = 7;
  if ( info != 0 ) {
    Blas2.xerbla( 'dsyr', info );
    return;
  }

  if ( n == 0 || alpha == 0. ) return;
  var i = -1;
  var ix = -1;
  var j = -1;
  var jx = -1;
  var kx = ( incx <= 0 ? 1 - ( n - 1 ) * incx : 1 );
  var temp = Number.POSITIVE_INFINITY;
  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    jx = kx;
    for ( j = 1; j <= n; j ++ ) {
      if ( x[ ioffx + jx - 1 ] != 0. ) {
        temp = alpha * x[ ioffx + jx - 1 ];
        ix = kx;
        for ( i = 1; i <= j; i ++ ) {
          A[ ioffA + i - 1 + ( j - 1 ) * lda ] +=
            x[ ioffx + ix - 1 ] * temp;
          ix += incx;
        }
      }
      jx += incx;
    }
  } else {
    jx = kx;
    for ( j = 1; j <= n; j ++ ) {
      if ( x[ ioffx + jx - 1 ] != 0. ) { 
        temp = alpha * x[ ioffx + jx - 1 ];
        ix = jx;
        for ( i = j; i <= n; i ++ ) {
          A[ ioffA + i - 1 + ( j - 1 ) * lda ] +=
            x[ ioffx + ix - 1 ] * temp;
          ix += incx;
        }
      }
      jx += incx;
    }
  }
}
Blas2.zher = function( uplo, n, alpha, x, incx, A, lda, ioffx, ioffA ) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U' 
  && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( n < 0 ) info = 2;
  else if ( incx == 0 ) info = 5;
  else if ( lda < Math.max( 1 , n ) ) info = 7;
  if ( info != 0 ) {
    Blas2.xerbla( 'zher', info );
    return;
  }

  if ( n == 0 || alpha == 0. ) return;
  var i = -1;
  var ix = -1;
  var j = -1;
  var jx = -1;
  var kx = ( incx <= 0 ? 1 - ( n - 1 ) * incx : 1 );
  var temp = new Complex();
  var zero = new Complex( 0., 0. );
  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    jx = kx;
    for ( j = 1; j <= n; j ++ ) {
      if ( ! x[ ioffx + jx - 1 ].equals( zero ) ) {
        temp.setValue( ComplexMath.timesNumber(
          ComplexMath.conj( x[ ioffx + jx - 1 ] ) , alpha ) );
        ix = kx;
        for ( i = 1; i <= j - 1; i ++ ) {
          A[ ioffA + i - 1  + ( j - 1 ) * lda ].plusEquals(
            ComplexMath.times( x[ ioffx + ix - 1 ] , temp ) );
          ix += incx;
        }
        A[ ioffA + j - 1 + ( j - 1 ) * lda ].plusEquals(
          ComplexMath.times( x[ ioffx + jx - 1 ] , temp ) );
      }
      A[ ioffA + j - 1 + ( j - 1 ) * lda ].pureReal();
      jx += incx;
    }
  } else {
    jx = kx;
    for ( j = 1; j <= n; j ++ ) {
      if ( ! x[ ioffx + jx - 1 ].equals( zero ) ) { 
        temp.setValue( ComplexMath.timesNumber(
          ComplexMath.conj( x[ ioffx + jx - 1 ] ) , alpha ) );
        A[ ioffA + j - 1 + ( j - 1 ) * lda ].plusEquals(
          ComplexMath.times( temp , x[ ioffx + jx - 1 ] ) );
        ix = jx;
        for ( i = j + 1; i <= n; i ++ ) {
          ix += incx;
          A[ ioffA + i - 1 + ( j - 1 ) * lda ].plusEquals(
            ComplexMath.times( x[ ioffx + ix - 1 ] , temp ) );
        }
      }
      A[ ioffA + j - 1 + ( j - 1 ) * lda ].pureReal();
      jx += incx;
    }
  }
}
//*************************************************************************
//  A += x alpha x^H, A symmetric packed
Blas2.dspr = function( uplo, n, alpha, x, incx, AP, ioffx, ioffAP ) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U' 
  && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( n < 0 ) info = 2;
  else if ( incx == 0 ) info = 5;
  if ( info != 0 ) {
    Blas2.xerbla( 'dspr', info );
    return;
  }

  if ( n == 0 || alpha == 0. ) return;
  var j = -1;
  var jx = -1;
  var k = -1;
  var kx = ( incx <= 0 ? 1 - ( n - 1 ) * incx : 1 );
  var temp = Number.POSITIVE_INFINITY;

  var ix = -1;
  var kk = 1;
  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    jx = kx;
    for ( j = 1; j <= n; j ++ ) {
      if ( x[ ioffx + jx - 1 ] != 0. ) {
        temp = alpha * x[ ioffx + jx - 1 ];
        ix = kx;
        for ( k = kk; k <= kk + j - 1; k ++ ) {
          AP[ ioffAP + k - 1 ] += x[ ioffx + ix - 1 ] * temp;
          ix += incx;
        }
      }
      jx += incx;
      kk += j;
    }
  } else {
    jx = kx;
    for ( j = 1; j <= n; j ++ ) {
      if ( x[ ioffx + jx - 1 ] != 0. ) {
        temp = alpha * x[ ioffx + jx - 1 ];
        ix = jx;
        for ( k = kk; k <= kk + n - j; k ++ ) {
          AP[ ioffAP + k - 1 ] += x[ ioffx + ix - 1 ] * temp;
          ix += incx;
        }
      }
      jx += incx;
      kk += n - j + 1;
    }
  }
}
Blas2.zhpr = function( uplo, n, alpha, x, incx, AP, ioffx, ioffAP ) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U' 
  && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( n < 0 ) info = 2;
  else if ( incx == 0 ) info = 5;
  if ( info != 0 ) {
    Blas2.xerbla( 'zhpr', info );
    return;
  }

  var zero = new Complex( 0., 0. );
  if ( n == 0 || alpha == 0. ) return;
  var j = -1;
  var jx = -1;
  var k = -1;
  var kx = ( incx <= 0 ? 1 - ( n - 1 ) * incx : 1 );
  var temp = new Complex();

  var ix = -1;
  var kk = 1;
  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    jx = kx;
    for ( j = 1; j <= n; j ++ ) {
      if ( ! x[ ioffx + jx - 1 ].equals( zero ) ) {
        temp.setValue( ComplexMath.timesNumber(
          ComplexMath.conj( x[ ioffx + jx - 1 ] ) , alpha ) );
        ix = kx;
        for ( k = kk; k <= kk + j - 2; k ++ ) {
          AP[ ioffAP + k - 1 ].plusEquals(
            ComplexMath.times( x[ ioffx + ix - 1 ] , temp ) );
          ix += incx;
        }
        AP[ ioffAP + kk + j - 2 ].plusEquals(
          ComplexMath.times( x[ ioffx + jx - 1 ] , temp ) );
      }
      AP[ ioffAP + kk + j  - 2 ].pureReal();
      jx += incx;
      kk += j;
    }
  } else {
    jx = kx;
    for ( j = 1; j <= n; j ++ ) {
      if ( ! x[ ioffx + jx - 1 ].equals( zero ) ) {
        temp.setValue( ComplexMath.timesNumber( 
          ComplexMath.conj( x[ ioffx + jx - 1 ] ) , alpha ) );
        AP[ ioffAP + kk - 1 ].plusEquals(
          ComplexMath.times( temp , x[ ioffx + jx - 1 ] ) );
        ix = jx;
        for ( k = kk + 1; k <= kk + n - j; k ++ ) {
          ix += incx;
          AP[ ioffAP + k - 1 ].plusEquals(
            ComplexMath.times( x[ ioffx + ix - 1 ] , temp ) );
        }
      }
      AP[ ioffAP + kk - 1 ].pureReal();
      jx += incx;
      kk += n - j + 1;
    }
  }
}
//*************************************************************************
//  A += x alpha y^T + y alpha x^T, A symmetric
Blas2.dsyr2 = function( uplo, n, alpha, x, incx, y, incy, A, lda, ioffx,
ioffy, ioffA ) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U' 
  && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( n < 0 ) info = 2;
  else if ( incx == 0 ) info = 5;
  else if ( incy == 0 ) info = 7;
  else if ( lda < Math.max( 1 , n ) ) info = 9;
  if ( info != 0 ) {
    Blas2.xerbla( 'dsyr2', info );
    return;
  }

  if ( n == 0 || alpha == 0. ) return;
  var kx = ( incx > 0 ? 1 : 1 - ( n - 1 ) * incx );
  var ky = ( incy > 0 ? 1 : 1 - ( n - 1 ) * incy );
  var i = -1;
  var ix = -1;
  var iy = -1;
  var j = -1;
  var jx = kx;
  var jy = ky;
  var temp1 = Number.POSITIVE_INFINITY;
  var temp2 = Number.POSITIVE_INFINITY;

  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    for ( j = 1; j <= n; j ++ ) {
      if ( x[ ioffx + jx - 1 ] != 0. || y[ ioffy + jy - 1 ] != 0. ) {
        temp1 = alpha * y[ ioffy + jy - 1 ];
        temp2 = alpha * x[ ioffx + jx - 1 ];
        ix = kx;
        iy = ky;
        for ( i = 1; i <= j; i ++ ) {
          A[ ioffA + i - 1 + ( j - 1 ) * lda ] +=
            x[ ioffx + ix - 1 ] * temp1
            + y[ ioffy + iy - 1 ] * temp2;
          ix += incx;
          iy += incy;
        }
      }
      jx += incx;
      jy += incy;
    }
  } else {
    for ( j = 1; j <= n; j ++ ) {
      if ( x[ ioffx + jx - 1 ] != 0. || y[ ioffy + jy - 1 ] != 0. ) {
        temp1 = alpha * y[ ioffy + jy - 1 ];
        temp2 = alpha * x[ ioffx + jx - 1 ];
        ix = jx;
        iy = jy;
        for ( i = j; i <= n; i ++ ) {
          A[ ioffA + i - 1 + ( j - 1 ) * lda ] +=
            x[ ioffx + ix - 1 ] * temp1
            + y[ ioffy + iy - 1 ] * temp2;
          ix += incx;
          iy += incy;
        }
      }
      jx += incx;
      jy += incy;
    }
  }
}
Blas2.zher2 = function( uplo, n, alpha, x, incx, y, incy, A, lda, ioffx,
ioffy, ioffA ) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U' 
  && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( n < 0 ) info = 2;
  else if ( incx == 0 ) info = 5;
  else if ( incy == 0 ) info = 7;
  else if ( lda < Math.max( 1 , n ) ) info = 9;
  if ( info != 0 ) {
    Blas2.xerbla( 'zher2', info );
    return;
  }

  var zero = new Complex( 0., 0. );
  if ( n == 0 || alpha.equals( zero ) ) return;
  var kx = ( incx > 0 ? 1 : 1 - ( n - 1 ) * incx );
  var ky = ( incy > 0 ? 1 : 1 - ( n - 1 ) * incy );
  var i = -1;
  var ix = -1;
  var iy = -1;
  var j = -1;
  var jx = kx;
  var jy = ky;
  var temp1 = new Complex();
  var temp2 = new Complex();

  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    for ( j = 1; j <= n; j ++ ) {
      if ( ! x[ ioffx + jx - 1 ].equals( zero ) ||
      ! y[ ioffy + jy - 1 ].equals( zero ) ) {
        temp1.setValue( ComplexMath.times( alpha ,
          ComplexMath.conj( y[ ioffy + jy - 1 ] ) ) );
        temp2.setValue( ComplexMath.conj(
          ComplexMath.times( alpha , x[ ioffx + jx - 1 ] ) ) );
        ix = kx;
        iy = ky;
        for ( i = 1; i <= j - 1; i ++ ) {
          A[ ioffA + i - 1 + ( j - 1 ) * lda ].plusEquals(
            ComplexMath.plus(
            ComplexMath.times( x[ ioffx + ix - 1 ] , temp1 ) ,
            ComplexMath.times( y[ ioffy + iy - 1 ] , temp2 ) ) );
          ix += incx;
          iy += incy;
        }
        A[ ioffA + j - 1 + ( j - 1 ) * lda ].plusEquals(
          ComplexMath.plus(
          ComplexMath.times( x[ ioffx + jx - 1 ] , temp1 ) ,
          ComplexMath.times( y[ ioffy + jy - 1 ] , temp2 ) ) );
      }
      A[ ioffA + j - 1 + ( j - 1 ) * lda ].pureReal();
      jx += incx;
      jy += incy;
    }
  } else {
    for ( j = 1; j <= n; j ++ ) {
      if ( ! x[ ioffx + jx - 1 ].equals( zero ) ||
      ! y[ ioffy + jy - 1 ].equals( zero ) ) {
        temp1.setValue( ComplexMath.times( alpha ,
          ComplexMath.conj( y[ ioffy + jy - 1 ] ) ) );
        temp2.setValue( ComplexMath.conj(
          ComplexMath.times( alpha , x[ ioffx + jx - 1 ] ) ) );
        A[ ioffA + j - 1 + ( j - 1 ) * lda ].plusEquals(
          ComplexMath.plus(
          ComplexMath.times( x[ ioffx + jx - 1 ] , temp1 ) ,
          ComplexMath.times( y[ ioffy + jy - 1 ] , temp2 ) ) );
        ix = jx;
        iy = jy;
        for ( i = j + 1; i <= n; i ++ ) {
          ix += incx;
          iy += incy;
          A[ ioffA + i - 1 + ( j - 1 ) * lda ].plusEquals(
            ComplexMath.plus(
            ComplexMath.times( x[ ioffx + ix - 1 ] , temp1 ) ,
            ComplexMath.times( y[ ioffy + iy - 1 ] , temp2 ) ) );
        }
      }
      A[ ioffA + j - 1 + ( j - 1 ) * lda ].pureReal();
      jx += incx;
      jy += incy;
    }
  }
}
//*************************************************************************
//  A += x alpha y^T + y alpha x^T, A symmetric packed
Blas2.dspr2 = function( uplo, n, alpha, x, incx, y, incy, AP, ioffx,
ioffy, ioffAP ) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U' 
  && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( n < 0 ) info = 2;
  else if ( incx == 0 ) info = 5;
  else if ( incy == 0 ) info = 7;
  if ( info != 0 ) {
    Blas2.xerbla( 'dspr2', info );
    return;
  }

  if ( n == 0 || alpha == 0. ) return;
  var kx = ( incx > 0 ? 1 : 1 - ( n - 1 ) * incx );
  var ky = ( incy > 0 ? 1 : 1 - ( n - 1 ) * incy );

  var ix = -1;
  var iy = -1;
  var j = -1;;
  var jx = kx;
  var jy = ky;
  var k = -1;
  var kk = 1;
  var temp1 = Number.POSITIVE_INFINITY;
  var temp2 = Number.POSITIVE_INFINITY;
  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    for ( j = 1; j <= n; j ++ ) {
      if ( x[ ioffx + jx - 1 ] != 0. || y[ ioffy + jy - 1 ] != 0. ) {
        temp1 = alpha * y[ ioffy + jy - 1 ];
        temp2 = alpha * x[ ioffx + jx - 1 ];
        ix = kx;
        iy = ky;
        for ( k = kk; k <= kk + j - 1; k ++ ) {
          AP[ ioffAP + k - 1 ] += x[ ioffx + ix - 1 ] * temp1
                           + y[ ioffy + iy - 1 ] * temp2;
          ix += incx;
          iy += incy;
        }
      }
      jx += incx;
      jy += incy;
      kk += j;
    }
  } else {
    for ( j = 1; j <= n; j ++ ) {
      if ( x[ ioffx + jx - 1 ] != 0. || y[ ioffy + jy - 1 ] != 0. ) {
        temp1 = alpha * y[ ioffy + jy - 1 ];
        temp2 = alpha * x[ ioffx + jx - 1 ];
        ix = jx;
        iy = jy;
        for ( k = kk; k <= kk + n - j; k ++ ) {
          AP[ ioffAP + k - 1 ] += x[ ioffx + ix - 1 ] * temp1
                           + y [ ioffy + iy - 1 ] * temp2;
          ix += incx;
          iy += incy;
        }
      }
      jx += incx;
      jy += incy;
      kk += n - j + 1;
    }
  }
}
Blas2.zhpr2 = function( uplo, n, alpha, x, incx, y, incy, AP, ioffx,
ioffy, ioffAP ) {
  var info = 0;
  if ( uplo.charAt(0).toUpperCase() != 'U' 
  && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( n < 0 ) info = 2;
  else if ( incx == 0 ) info = 5;
  else if ( incy == 0 ) info = 7;
  if ( info != 0 ) {
    Blas2.xerbla( 'zhpr2', info );
    return;
  }

  var zero = new Complex( 0., 0. );
  if ( n == 0 || alpha.equals(zero) ) return;
  var kx = ( incx > 0 ? 1 : 1 - ( n - 1 ) * incx );
  var ky = ( incy > 0 ? 1 : 1 - ( n - 1 ) * incy );

  var ix = -1;
  var iy = -1;
  var j = -1;;
  var jx = kx;
  var jy = ky;
  var k = -1;
  var kk = 1;
  var temp1 = new Complex();
  var temp2 = new Complex();
  if ( uplo.charAt(0).toUpperCase() == 'U' ) {
    for ( j = 1; j <= n; j ++ ) {
      if ( ! x[ ioffx + jx - 1 ].equals( zero ) ||
      ! y[ ioffy + jy - 1 ].equals( zero ) ) {
        temp1.setValue( ComplexMath.times( alpha ,
          ComplexMath.conj( y[ ioffy + jy - 1 ] ) ) );
        temp2.setValue( ComplexMath.conj(
          ComplexMath.times( alpha , x[ ioffx + jx - 1 ] ) ) );
        ix = kx;
        iy = ky;
        for ( k = kk; k <= kk + j - 2; k ++ ) {
          AP[ ioffAP + k - 1 ].plusEquals( ComplexMath.plus(
            ComplexMath.times( x[ ioffx + ix - 1 ] , temp1 ) ,
            ComplexMath.times( y[ ioffy + iy - 1 ] , temp2 ) ) );
          ix += incx;
          iy += incy;
        }
        AP[ ioffAP + kk + j - 2 ].plusEquals( ComplexMath.plus(
          ComplexMath.times( x[ ioffx + jx - 1 ] , temp1 ) ,
          ComplexMath.times( y[ ioffy + jy - 1 ] , temp2 ) ) );
      }
      AP[ ioffAP + kk + j - 2 ].pureReal();
      jx += incx;
      jy += incy;
      kk += j;
    }
  } else {
    for ( j = 1; j <= n; j ++ ) {
      if ( ! x[ ioffx + jx - 1 ].equals( zero ) ||
      ! y[ ioffy + jy - 1 ].equals( zero ) ) {
        temp1.setValue( ComplexMath.times( alpha ,
          ComplexMath.conj( y[ ioffy + jy - 1 ] ) ) );
        temp2.setValue( ComplexMath.conj(
          ComplexMath.times( alpha , x[ ioffx + jx - 1 ] ) ) );
        AP[ ioffAP + kk - 1 ].plusEquals( ComplexMath.plus(
          ComplexMath.times( x[ ioffx + jx - 1 ] , temp1 ) ,
          ComplexMath.times( y[ ioffy + jy - 1 ] , temp2 ) ) );
        ix = jx;
        iy = jy;
        for ( k = kk + 1; k <= kk + n - j; k ++ ) {
          ix += incx;
          iy += incy;
          AP[ ioffAP + k - 1 ].plusEquals( ComplexMath.plus(
            ComplexMath.times( x[ ioffx + ix - 1 ] , temp1 ) ,
            ComplexMath.times( y[ ioffy + iy - 1 ] , temp2 ) ) );
        }
      } 
      AP[ ioffAP + kk - 1 ].pureReal();
      jx += incx;
      jy += incy;
      kk += n - j + 1;
    }
  }
}

function testBlas2() {
//test_dtrmv();
//test_ztrmv();
//test_dtbmv();
/*test_ztbmv();*/
/*test_dtpmv();*/
/*test_ztpmv();*/

//test_dtrsv();
//test_ztrsv();
/*test_dtbsv();*/
/*test_ztbsv();*/
/*test_dtpsv();*/
/*test_ztpsv();*/

//test_dgemv();
//test_zgemv();
/*test_dgbmv();*/
/*test_zgbmv();*/
//test_dsymv();
//test_zhemv();
/*test_dsbmv();*/
/*test_zhbmv();*/
/*test_dspmv();*/
/*test_zhpmv();*/

//test_dger();
//test_zgerc();
//test_zgeru();

//test_dsyr();
//test_zher();
/*test_dspr();*/
/*test_zhpr();*/
  test_dsyr2();
  test_zher2();
/*test_dspr2();*/
/*test_zhpr2();*/
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dtrmv() {
  document.getElementById("debug_textarea").value +=
    "testing dtrmv *************" + "\n";
  var ioffA = 1;
  var ioffx = 2;
  var ioffy = 3;
  var A = new Array( ioffA + 12 );
  var x = new Array( ioffx + 3 );
  var y = new Array( ioffy + 6 );
  var n = 3;
  var lda = 4;
  A[ ioffA + 0 ] = 2.;
  A[ ioffA + 1 ] = 1.;
  A[ ioffA + 2 ] = 0.;
  A[ ioffA + 4 ] = 3.;
  A[ ioffA + 5 ] = 3.;
  A[ ioffA + 6 ] = 2.;
  A[ ioffA + 8 ] = 4.;
  A[ ioffA + 9 ] = 4.;
  A[ ioffA + 10 ] = 4.;

  var i = -1;
  var incx = 1;
  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  var incy = 2;
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrmv('U','N','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
    "dtrmv('U','N','U',n,A,lda,x,incx) , x = "
        + x[ioffx + 0] + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrmv('U','N','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
    "dtrmv('U','N','U',n,A,lda,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " "
        + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrmv('U','N','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
    "dtrmv('U','N','N',n,A,lda,x,incx) , x = "
        + x[ioffx + 0] + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrmv('U','N','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
    "dtrmv('U','N','N',n,A,lda,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " "
        + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrmv('U','T','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
    "dtrmv('U','T','U',n,A,lda,x,incx) , x = "
        + x[ioffx + 0] + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrmv('U','T','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
    "dtrmv('U','T','U',n,A,lda,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " "
        + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrmv('U','T','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
    "dtrmv('U','T','N',n,A,lda,x,incx) , x = "
        + x[ioffx + 0] + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrmv('U','T','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
    "dtrmv('U','T','N',n,A,lda,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " "
        + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrmv('U','C','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
    "dtrmv('U','C','U',n,A,lda,x,incx) , x = "
        + x[ioffx + 0] + " " + x[ioffx + 1]
        + " " + x[ioffx + 2] + "\n";
  Blas2.dtrmv('U','C','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
    "dtrmv('U','C','U',n,A,lda,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " "
        + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrmv('U','C','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
    "dtrmv('U','C','N',n,A,lda,x,incx) , x = "
        + x[ioffx + 0] + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrmv('U','C','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
    "dtrmv('U','C','N',n,A,lda,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " "
        + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrmv('L','N','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
    "dtrmv('L','N','U',n,A,lda,x,incx) , x = "
        + x[ioffx + 0] + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrmv('L','N','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
    "dtrmv('L','N','U',n,A,lda,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " "
        + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrmv('L','N','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
    "dtrmv('L','N','N',n,A,lda,x,incx) , x = "
        + x[ioffx + 0] + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrmv('L','N','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
    "dtrmv('L','N','N',n,A,lda,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " "
        + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrmv('L','T','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
    "dtrmv('L','T','U',n,A,lda,x,incx) , x = "
        + x[ioffx + 0] + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrmv('L','T','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
    "dtrmv('L','T','U',n,A,lda,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " "
        + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrmv('L','T','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
    "dtrmv('L','T','N',n,A,lda,x,incx) , x = "
        + x[ioffx + 0] + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrmv('L','T','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
    "dtrmv('L','T','N',n,A,lda,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " "
        + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrmv('L','C','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
    "dtrmv('L','C','U',n,A,lda,x,incx) , x = "
        + x[ioffx + 0] + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrmv('L','C','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
    "dtrmv('L','C','U',n,A,lda,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " "
        + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrmv('L','C','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
    "dtrmv('L','C','N',n,A,lda,x,incx) , x = "
        + x[ioffx + 0] + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrmv('L','C','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
    "dtrmv('L','C','N',n,A,lda,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " "
        + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_ztrmv() {
  document.getElementById("debug_textarea").value +=
   "testing ztrmv *************" + "\n";
  var ioffA = 1;
  var ioffx = 2;
  var ioffy = 3;
  var A = new Array( ioffA + 12 );
  var x = new Array( ioffx + 3 );
  var y = new Array( ioffy + 6 );
  var n = 3;
  var lda = 4;
  A[ ioffA + 0 ] = new Complex( 2. , 3. );
  A[ ioffA + 1 ] = new Complex( 1. , 2. );
  A[ ioffA + 2 ] = new Complex( 0. , 1. );
  A[ ioffA + 4 ] = new Complex( 3. , 4. );
  A[ ioffA + 5 ] = new Complex( 3. , 3. );
  A[ ioffA + 6 ] = new Complex( 2. , 2. );
  A[ ioffA + 8 ] = new Complex( 4. , 5. );
  A[ ioffA + 9 ] = new Complex( 4. , 4. );
  A[ ioffA + 10 ] = new Complex( 4. , 3. );

  var i = -1;
  var incx = 1;
  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  var incy = 2;
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrmv('U','N','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
 "ztrmv('U','N','U',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrmv('U','N','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
 "ztrmv('U','N','U',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] +"\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrmv('U','N','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
 "ztrmv('U','N','N',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrmv('U','N','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
 "ztrmv('U','N','N',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrmv('U','T','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrmv('U','T','U',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
          + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrmv('U','T','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrmv('U','T','U',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
          + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
          + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrmv('U','T','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrmv('U','T','N',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
          + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrmv('U','T','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrmv('U','T','N',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
          + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
          + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrmv('U','C','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrmv('U','C','U',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
          + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrmv('U','C','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrmv('U','C','U',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
          + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
          + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrmv('U','C','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrmv('U','C','N',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
          + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrmv('U','C','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
    "ztrmv('U','C','N',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
          + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
          + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrmv('L','N','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrmv('L','N','U',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
          + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrmv('L','N','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrmv('L','N','U',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
          + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
          + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrmv('L','N','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrmv('L','N','N',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
          + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrmv('L','N','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrmv('L','N','N',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
          + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
          + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrmv('L','T','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrmv('L','T','U',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
          + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrmv('L','T','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrmv('L','T','U',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
          + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
          + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrmv('L','T','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrmv('L','T','N',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
          + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrmv('L','T','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrmv('L','T','N',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
          + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
          + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrmv('L','C','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrmv('L','C','U',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
          + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrmv('L','C','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrmv('L','C','U',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
          + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
          + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrmv('L','C','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrmv('L','C','N',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
          + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrmv('L','C','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrmv('L','C','N',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
          + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
          + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dtbmv() {
  document.getElementById("debug_textarea").value +=
   "testing dtbmv *************"+ "\n";
  var ioffA = 1;
  var ioffB = 2;
  var ioffx = 3;
  var ioffy = 4;
  var A = new Array( ioffA + 12 );
  var B = new Array( ioffB + 12 );
  var x = new Array( ioffx + 3 );
  var y = new Array( ioffy + 6 );
  var n = 3;
  var lda = 4;
  var ldb = 4;
  A[ ioffA + 0 ] = 2.;
  A[ ioffA + 1 ] = 1.;
  A[ ioffA + 2 ] = 0.;
  A[ ioffA + 4 ] = 5.;
  A[ ioffA + 5 ] = 3.;
  A[ ioffA + 6 ] = -1.;
  A[ ioffA + 8 ] = 0.;
  A[ ioffA + 9 ] = 4.;
  A[ ioffA + 10 ] = 6.;
  var k = 1;

  var i = -1;
  var j = -1;
  var m = -1;
  for ( j = 0; j < n; j ++ ) {
    m = k - j;
    for ( i = Math.max( 0 , j - k ); i <= j; i ++ ) {
      B[ ioffB + m + i + j * ldb ] = A[ ioffA + i + j * lda ];
    }
  }
  var incx = 1;
  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  var incy = 2;
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbmv('U','N','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbmv('U','N','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
          + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbmv('U','N','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbmv('U','N','U',3,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
          + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
          + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbmv('U','N','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbmv('U','N','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
          + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbmv('U','N','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbmv('U','N','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
          + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
          + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbmv('U','T','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbmv('U','T','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbmv('U','T','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbmv('U','T','U',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbmv('U','T','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbmv('U','T','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbmv('U','T','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbmv('U','T','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbmv('U','C','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbmv('U','C','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbmv('U','C','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbmv('U','C','U',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbmv('U','C','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbmv('U','C','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbmv('U','C','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbmv('U','C','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( j = 0; j < n; j ++ ) {
    m = - j;
    for ( i = j; i <= Math.min( n - 1 , j + k ); i ++ ) {
      B[ ioffB + m + i + j * ldb ] = A[ ioffA + i + j * lda ];
    }
  }
  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbmv('L','N','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbmv('L','N','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbmv('L','N','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbmv('L','N','U',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbmv('L','N','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbmv('L','N','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbmv('L','N','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbmv('L','N','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbmv('L','T','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbmv('L','T','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbmv('L','T','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbmv('L','T','U',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbmv('L','T','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbmv('L','T','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbmv('L','T','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbmv('L','T','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbmv('L','C','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbmv('L','C','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbmv('L','C','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbmv('L','C','U',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbmv('L','C','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbmv('L','C','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbmv('L','C','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbmv('L','C','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_ztbmv() {
  document.getElementById("debug_textarea").value +=
   "testing ztbmv *************"+ "\n";
  var ioffA = 1;
  var ioffB = 2;
  var ioffx = 3;
  var ioffy = 4;
  var A = new Array( ioffA + 12 );
  var B = new Array( ioffB + 12 );
  var x = new Array( ioffx + 3 );
  var y = new Array( ioffy + 6 );
  var n = 3;
  var lda = 4;
  var ldb = 4;
  A[ ioffA + 0 ] = new Complex( 2. , 3. );
  A[ ioffA + 1 ] = new Complex( 1. , 2. );
  A[ ioffA + 2 ] = new Complex( 0. , 0. );
  A[ ioffA + 4 ] = new Complex( 5. , 4. );
  A[ ioffA + 5 ] = new Complex( 3. , 3. );
  A[ ioffA + 6 ] = new Complex( -1. , 2. );
  A[ ioffA + 8 ] = new Complex( 0. , 0. );
  A[ ioffA + 9 ] = new Complex( 4. , 4. );
  A[ ioffA + 10 ] = new Complex( 6. , 3. );
  var k = 1;

  var i = -1;
  var j = -1;
  var m = -1;
  for ( j = 0; j < n; j ++ ) {
    m = k - j;
    for ( i = Math.max( 0 , j - k ); i <= j; i ++ ) {
      B[ ioffB + m + i + j * ldb ] = new Complex();
      B[ ioffB + m + i + j * ldb ].value = A[ ioffA + i + j * lda ];
    }
  }
  var incx = 1;
  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  var incy = 2;
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbmv('U','N','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbmv('U','N','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
/*
  Blas2.ztbmv('U','N','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbmv('U','N','U',3,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbmv('U','N','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbmv('U','N','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbmv('U','N','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbmv('U','N','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbmv('U','T','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbmv('U','T','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbmv('U','T','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbmv('U','T','U',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbmv('U','T','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbmv('U','T','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbmv('U','T','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbmv('U','T','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbmv('U','C','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbmv('U','C','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbmv('U','C','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbmv('U','C','U',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbmv('U','C','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbmv('U','C','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbmv('U','C','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbmv('U','C','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( j = 0; j < n; j ++ ) {
    m = - j;
    for ( i = j; i <= Math.min( n - 1 , j + k ); i ++ ) {
      B[ ioffB + m + i + j * ldb ] = new Complex();
      B[ ioffB + m + i + j * ldb ].value = A[ ioffA + i + j * lda ];
    }
  }
  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbmv('L','N','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbmv('L','N','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbmv('L','N','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbmv('L','N','U',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbmv('L','N','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbmv('L','N','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbmv('L','N','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbmv('L','N','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbmv('L','T','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbmv('L','T','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbmv('L','T','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbmv('L','T','U',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbmv('L','T','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbmv('L','T','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbmv('L','T','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbmv('L','T','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbmv('L','C','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbmv('L','C','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbmv('L','C','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbmv('L','C','U',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbmv('L','C','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbmv('L','C','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbmv('L','C','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbmv('L','C','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";
*/
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dtpmv() {
  document.getElementById("debug_textarea").value +=
   "testing dtpmv *************"+ "\n";
  var ioffA = 1;
  var ioffP = 2;
  var ioffx = 3;
  var ioffy = 4;
  var A = new Array( ioffA + 12 );
  var P = new Array( ioffP + 6 );
  var x = new Array( ioffx + 3 );
  var y = new Array( ioffy + 6 );
  var n = 3;
  var lda = 4;
  A[ ioffA + 0 ] = 2.;
  A[ ioffA + 1 ] = 1.;
  A[ ioffA + 2 ] = 0.;
  A[ ioffA + 4 ] = 5.;
  A[ ioffA + 5 ] = 3.;
  A[ ioffA + 6 ] = -1.;
  A[ ioffA + 8 ] = 0.;
  A[ ioffA + 9 ] = 4.;
  A[ ioffA + 10 ] = 6.;

  var i = -1;
  var j = -1;
  var m = -1;
  var k = 0;
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i <= j; i ++ ) {
      P[ ioffP + k++ ] = A[ ioffA + i + j * lda ];
    }
  }
  var incx = 1;
  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  var incy = 2;
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpmv('U','N','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpmv('U','N','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";

  Blas2.dtpmv('U','N','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpmv('U','N','U',3,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpmv('U','N','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpmv('U','N','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpmv('U','N','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpmv('U','N','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpmv('U','T','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpmv('U','T','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpmv('U','T','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpmv('U','T','U',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpmv('U','T','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpmv('U','T','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpmv('U','T','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpmv('U','T','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpmv('U','C','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpmv('U','C','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpmv('U','C','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpmv('U','C','U',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpmv('U','C','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpmv('U','C','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpmv('U','C','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpmv('U','C','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  k = 0;
  for ( j = 0; j < n; j ++ ) {
    for ( i = j; i < n; i ++ ) {
      P[ ioffP + k++ ] = A[ ioffA + i + j * lda ];
    }
  }
  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpmv('L','N','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpmv('L','N','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpmv('L','N','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpmv('L','N','U',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpmv('L','N','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpmv('L','N','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpmv('L','N','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpmv('L','N','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpmv('L','T','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpmv('L','T','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpmv('L','T','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpmv('L','T','U',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpmv('L','T','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpmv('L','T','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpmv('L','T','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpmv('L','T','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpmv('L','C','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpmv('L','C','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpmv('L','C','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpmv('L','C','U',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpmv('L','C','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpmv('L','C','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpmv('L','C','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpmv('L','C','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_ztpmv() {
  document.getElementById("debug_textarea").value +=
   "testing ztpmv *************"+ "\n";
  var ioffA = 1;
  var ioffP = 2;
  var ioffx = 3;
  var ioffy = 4;
  var A = new Array( ioffA + 12 );
  var P = new Array( ioffP + 6 );
  var x = new Array( ioffx + 3 );
  var y = new Array( ioffy + 6 );
  var n = 3;
  var lda = 4;
  A[ ioffA + 0 ] = new Complex( 2. , 3. );
  A[ ioffA + 1 ] = new Complex( 1. , 2. );
  A[ ioffA + 2 ] = new Complex( 0. , 0. );
  A[ ioffA + 4 ] = new Complex( 5. , 4. );
  A[ ioffA + 5 ] = new Complex( 3. , 3. );
  A[ ioffA + 6 ] = new Complex( -1. , 2. );
  A[ ioffA + 8 ] = new Complex( 0. , 0. );
  A[ ioffA + 9 ] = new Complex( 4. , 4. );
  A[ ioffA + 10 ] = new Complex( 6. , 3. );

  var i = -1;
  var j = -1;
  var m = -1;
  var k = 0;
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i <= j; i ++ ) {
      P[ ioffP + k ] = new Complex();
      P[ ioffP + k++ ].value = A[ ioffA + i + j * lda ];
    }
  }
  var incx = 1;
  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  var incy = 2;
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpmv('U','N','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpmv('U','N','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpmv('U','N','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpmv('U','N','U',3,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpmv('U','N','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpmv('U','N','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpmv('U','N','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpmv('U','N','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpmv('U','T','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpmv('U','T','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpmv('U','T','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpmv('U','T','U',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpmv('U','T','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpmv('U','T','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpmv('U','T','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpmv('U','T','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpmv('U','C','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpmv('U','C','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpmv('U','C','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpmv('U','C','U',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpmv('U','C','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpmv('U','C','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpmv('U','C','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpmv('U','C','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  k = 0;
  for ( j = 0; j < n; j ++ ) {
    for ( i = j; i < n; i ++ ) {
      P[ ioffP + k ] = new Complex();
      P[ ioffP + k++ ].value = A[ ioffA + i + j * lda ];
    }
  }
  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpmv('L','N','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpmv('L','N','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpmv('L','N','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpmv('L','N','U',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpmv('L','N','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpmv('L','N','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpmv('L','N','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpmv('L','N','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpmv('L','T','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpmv('L','T','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpmv('L','T','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpmv('L','T','U',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpmv('L','T','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpmv('L','T','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpmv('L','T','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpmv('L','T','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpmv('L','C','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpmv('L','C','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpmv('L','C','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpmv('L','C','U',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpmv('L','C','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpmv('L','C','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpmv('L','C','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpmv('L','C','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dtrsv() {
  document.getElementById("debug_textarea").value +=
   "testing dtrsv *************"+ "\n";
  var ioffA = 1;
  var ioffx = 2;
  var ioffy = 3;
  var A = new Array( ioffA + 12 );
  var x = new Array( ioffx + 3 );
  var y = new Array( ioffy + 6 );
  var n = 3;
  var lda = 4;
  A[ ioffA + 0 ] = 2.;
  A[ ioffA + 1 ] = 1.;
  A[ ioffA + 2 ] = 0.;
  A[ ioffA + 4 ] = 3.;
  A[ ioffA + 5 ] = 3.;
  A[ ioffA + 6 ] = 2.;
  A[ ioffA + 8 ] = 4.;
  A[ ioffA + 9 ] = 4.;
  A[ ioffA + 10 ] = 4.;

  var i = -1;
  var incx = 1;
  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  var incy = 2;
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrsv('U','N','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtrsv('U','N','U',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrsv('U','N','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtrsv('U','N','U',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrsv('U','N','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtrsv('U','N','N',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrsv('U','N','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtrsv('U','N','N',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrsv('U','T','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtrsv('U','T','U',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrsv('U','T','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtrsv('U','T','U',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrsv('U','T','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtrsv('U','T','N',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrsv('U','T','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtrsv('U','T','N',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrsv('U','C','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtrsv('U','C','U',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrsv('U','C','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtrsv('U','C','U',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrsv('U','C','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtrsv('U','C','N',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrsv('U','C','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtrsv('U','C','N',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrsv('L','N','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtrsv('L','N','U',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrsv('L','N','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtrsv('L','N','U',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrsv('L','N','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtrsv('L','N','N',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrsv('L','N','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtrsv('L','N','N',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrsv('L','T','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtrsv('L','T','U',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrsv('L','T','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtrsv('L','T','U',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrsv('L','T','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
    "dtrsv('L','T','N',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrsv('L','T','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtrsv('L','T','N',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrsv('L','C','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtrsv('L','C','U',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrsv('L','C','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtrsv('L','C','U',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtrsv('L','C','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtrsv('L','C','N',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtrsv('L','C','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtrsv('L','C','N',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_ztrsv() {
  document.getElementById("debug_textarea").value +=
   "testing ztrsv *************"+ "\n";
  var ioffA = 1;
  var ioffx = 2;
  var ioffy = 3;
  var A = new Array( ioffA + 12 );
  var x = new Array( ioffx + 3 );
  var y = new Array( ioffy + 6 );
  var n = 3;
  var lda = 4;
  A[ ioffA + 0 ] = new Complex( 2. , 3. );
  A[ ioffA + 1 ] = new Complex( 1. , 2. );
  A[ ioffA + 2 ] = new Complex( 0. , 1. );
  A[ ioffA + 4 ] = new Complex( 3. , 4. );
  A[ ioffA + 5 ] = new Complex( 3. , 3. );
  A[ ioffA + 6 ] = new Complex( 2. , 2. );
  A[ ioffA + 8 ] = new Complex( 4. , 5. );
  A[ ioffA + 9 ] = new Complex( 4. , 4. );
  A[ ioffA + 10 ] = new Complex( 4. , 3. );

  var i = -1;
  var incx = 1;
  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  var incy = 2;
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrsv('U','N','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrsv('U','N','U',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrsv('U','N','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrsv('U','N','U',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrsv('U','N','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrsv('U','N','N',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrsv('U','N','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrsv('U','N','N',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrsv('U','T','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
    "ztrsv('U','T','U',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrsv('U','T','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrsv('U','T','U',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrsv('U','T','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrsv('U','T','N',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrsv('U','T','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrsv('U','T','N',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrsv('U','C','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrsv('U','C','U',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrsv('U','C','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrsv('U','C','U',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrsv('U','C','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrsv('U','C','N',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrsv('U','C','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrsv('U','C','N',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrsv('L','N','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrsv('L','N','U',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrsv('L','N','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrsv('L','N','U',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrsv('L','N','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrsv('L','N','N',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrsv('L','N','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrsv('L','N','N',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrsv('L','T','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrsv('L','T','U',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrsv('L','T','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrsv('L','T','U',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrsv('L','T','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrsv('L','T','N',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrsv('L','T','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrsv('L','T','N',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrsv('L','C','U',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrsv('L','C','U',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrsv('L','C','U',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrsv('L','C','U',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztrsv('L','C','N',n,A,lda,x,incx,ioffA,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztrsv('L','C','N',n,A,lda,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztrsv('L','C','N',n,A,lda,y,incy,ioffA,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztrsv('L','C','N',n,A,lda,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dtbsv() {
  document.getElementById("debug_textarea").value +=
   "testing dtbsv *************"+ "\n";
  var ioffA = 1;
  var ioffB = 2;
  var ioffx = 3;
  var ioffy = 4;
  var A = new Array( ioffA + 12 );
  var B = new Array( ioffB + 12 );
  var x = new Array( ioffx + 3 );
  var y = new Array( ioffy + 6 );
  var n = 3;
  var lda = 4;
  var ldb = 4;
  A[ ioffA + 0 ] = 2.;
  A[ ioffA + 1 ] = 1.;
  A[ ioffA + 2 ] = 0.;
  A[ ioffA + 4 ] = 5.;
  A[ ioffA + 5 ] = 3.;
  A[ ioffA + 6 ] = -1.;
  A[ ioffA + 8 ] = 0.;
  A[ ioffA + 9 ] = 4.;
  A[ ioffA + 10 ] = 6.;
  var k = 1;

  var i = -1;
  var j = -1;
  var m = -1;
  for ( j = 0; j < n; j ++ ) {
    m = k - j;
    for ( i = Math.max( 0 , j - k ); i <= j; i ++ ) {
      B[ ioffB + m + i + j * ldb ] = A[ ioffA + i + j * lda ];
    }
  }
  var incx = 1;
  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  var incy = 2;
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbsv('U','N','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbsv('U','N','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbsv('U','N','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbsv('U','N','U',3,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbsv('U','N','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbsv('U','N','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbsv('U','N','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbsv('U','N','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbsv('U','T','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbsv('U','T','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbsv('U','T','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbsv('U','T','U',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbsv('U','T','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbsv('U','T','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbsv('U','T','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbsv('U','T','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbsv('U','C','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbsv('U','C','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbsv('U','C','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbsv('U','C','U',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbsv('U','C','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbsv('U','C','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbsv('U','C','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbsv('U','C','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( j = 0; j < n; j ++ ) {
    m = - j;
    for ( i = j; i <= Math.min( n - 1 , j + k ); i ++ ) {
      B[ ioffB + m + i + j * ldb ] = A[ ioffA + i + j * lda ];
    }
  }
  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbsv('L','N','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbsv('L','N','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbsv('L','N','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbsv('L','N','U',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbsv('L','N','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbsv('L','N','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbsv('L','N','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbsv('L','N','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbsv('L','T','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbsv('L','T','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbsv('L','T','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbsv('L','T','U',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbsv('L','T','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbsv('L','T','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbsv('L','T','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbsv('L','T','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbsv('L','C','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbsv('L','C','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbsv('L','C','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbsv('L','C','U',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtbsv('L','C','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtbsv('L','C','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtbsv('L','C','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtbsv('L','C','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_ztbsv() {
  document.getElementById("debug_textarea").value +=
   "testing ztbsv *************"+ "\n";
  var ioffA = 1;
  var ioffB = 2;
  var ioffx = 3;
  var ioffy = 4;
  var A = new Array( ioffA + 12 );
  var B = new Array( ioffB + 12 );
  var x = new Array( ioffx + 3 );
  var y = new Array( ioffy + 6 );
  var n = 3;
  var lda = 4;
  var ldb = 4;
  A[ ioffA + 0 ] = new Complex( 2. , 3. );
  A[ ioffA + 1 ] = new Complex( 1. , 2. );
  A[ ioffA + 2 ] = new Complex( 0. , 0. );
  A[ ioffA + 4 ] = new Complex( 5. , 4. );
  A[ ioffA + 5 ] = new Complex( 3. , 3. );
  A[ ioffA + 6 ] = new Complex( -1. , 2. );
  A[ ioffA + 8 ] = new Complex( 0. , 0. );
  A[ ioffA + 9 ] = new Complex( 4. , 4. );
  A[ ioffA + 10 ] = new Complex( 6. , 3. );
  var k = 1;

  var i = -1;
  var j = -1;
  var m = -1;
  for ( j = 0; j < n; j ++ ) {
    m = k - j;
    for ( i = Math.max( 0 , j - k ); i <= j; i ++ ) {
      B[ ioffB + m + i + j * ldb ] = new Complex();
      B[ ioffB + m + i + j * ldb ].value = A[ ioffA + i + j * lda ];
    }
  }
  var incx = 1;
  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  var incy = 2;
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbsv('U','N','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbsv('U','N','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbsv('U','N','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbsv('U','N','U',3,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbsv('U','N','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbsv('U','N','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbsv('U','N','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbsv('U','N','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbsv('U','T','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbsv('U','T','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbsv('U','T','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbsv('U','T','U',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbsv('U','T','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbsv('U','T','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbsv('U','T','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbsv('U','T','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbsv('U','C','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbsv('U','C','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbsv('U','C','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbsv('U','C','U',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbsv('U','C','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbsv('U','C','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbsv('U','C','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbsv('U','C','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( j = 0; j < n; j ++ ) {
    m = - j;
    for ( i = j; i <= Math.min( n - 1 , j + k ); i ++ ) {
      B[ ioffB + m + i + j * ldb ] = new Complex();
      B[ ioffB + m + i + j * ldb ].value = A[ ioffA + i + j * lda ];
    }
  }
  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbsv('L','N','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbsv('L','N','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbsv('L','N','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbsv('L','N','U',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbsv('L','N','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbsv('L','N','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbsv('L','N','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbsv('L','N','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbsv('L','T','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbsv('L','T','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbsv('L','T','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbsv('L','T','U',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbsv('L','T','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbsv('L','T','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbsv('L','T','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbsv('L','T','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbsv('L','C','U',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbsv('L','C','U',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbsv('L','C','U',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbsv('L','C','U',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztbsv('L','C','N',n,k,B,ldb,x,incx,ioffB,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztbsv('L','C','N',n,k,B,ldb,x,incx) , x = " + x[ioffx + 0]
        + " " + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztbsv('L','C','N',n,k,B,ldb,y,incy,ioffB,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztbsv('L','C','N',n,k,B,ldb,y,incy) , y = " + y[ioffy + 0]
        + " " + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3]
        + " " + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dtpsv() {
  document.getElementById("debug_textarea").value +=
   "testing dtpsv *************"+ "\n";
  var ioffA = 1;
  var ioffP = 2;
  var ioffx = 3;
  var ioffy = 4;
  var A = new Array( ioffA + 12 );
  var P = new Array( ioffP + 6 );
  var x = new Array( ioffx + 3 );
  var y = new Array( ioffy + 6 );
  var n = 3;
  var lda = 4;
  A[ ioffA + 0 ] = 2.;
  A[ ioffA + 1 ] = 1.;
  A[ ioffA + 2 ] = 0.;
  A[ ioffA + 4 ] = 5.;
  A[ ioffA + 5 ] = 3.;
  A[ ioffA + 6 ] = -1.;
  A[ ioffA + 8 ] = 0.;
  A[ ioffA + 9 ] = 4.;
  A[ ioffA + 10 ] = 6.;

  var i = -1;
  var j = -1;
  var m = -1;
  var k = 0;
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i <= j; i ++ ) {
      P[ ioffP + k++ ] = A[ ioffA + i + j * lda ];
    }
  }
  var incx = 1;
  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  var incy = 2;
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpsv('U','N','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpsv('U','N','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpsv('U','N','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpsv('U','N','U',3,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpsv('U','N','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpsv('U','N','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpsv('U','N','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpsv('U','N','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpsv('U','T','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpsv('U','T','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpsv('U','T','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpsv('U','T','U',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpsv('U','T','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpsv('U','T','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpsv('U','T','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpsv('U','T','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpsv('U','C','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpsv('U','C','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpsv('U','C','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpsv('U','C','U',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpsv('U','C','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpsv('U','C','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpsv('U','C','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpsv('U','C','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  k = 0;
  for ( j = 0; j < n; j ++ ) {
    for ( i = j; i < n; i ++ ) {
      P[ ioffP + k++ ] = A[ ioffA + i + j * lda ];
    }
  }
  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpsv('L','N','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpsv('L','N','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpsv('L','N','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpsv('L','N','U',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpsv('L','N','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpsv('L','N','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpsv('L','N','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpsv('L','N','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpsv('L','T','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpsv('L','T','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpsv('L','T','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpsv('L','T','U',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpsv('L','T','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpsv('L','T','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpsv('L','T','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpsv('L','T','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpsv('L','C','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpsv('L','C','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpsv('L','C','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpsv('L','C','U',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) { x[ ioffx + i ] = Number( i - 1 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dtpsv('L','C','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "dtpsv('L','C','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.dtpsv('L','C','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "dtpsv('L','C','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_ztpsv() {
  document.getElementById("debug_textarea").value +=
   "testing ztpsv *************"+ "\n";
  var ioffA = 1;
  var ioffP = 2;
  var ioffx = 3;
  var ioffy = 4;
  var A = new Array( ioffA + 12 );
  var P = new Array( ioffP + 6 );
  var x = new Array( ioffx + 3 );
  var y = new Array( ioffy + 6 );
  var n = 3;
  var lda = 4;
  A[ ioffA + 0 ] = new Complex( 2. , 3. );
  A[ ioffA + 1 ] = new Complex( 1. , 2. );
  A[ ioffA + 2 ] = new Complex( 0. , 0. );
  A[ ioffA + 4 ] = new Complex( 5. , 4. );
  A[ ioffA + 5 ] = new Complex( 3. , 3. );
  A[ ioffA + 6 ] = new Complex( -1. , 2. );
  A[ ioffA + 8 ] = new Complex( 0. , 0. );
  A[ ioffA + 9 ] = new Complex( 4. , 4. );
  A[ ioffA + 10 ] = new Complex( 6. , 3. );

  var i = -1;
  var j = -1;
  var m = -1;
  var k = 0;
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i <= j; i ++ ) {
      P[ ioffP + k ] = new Complex();
      P[ ioffP + k++ ].value = A[ ioffA + i + j * lda ];
    }
  }
  var incx = 1;
  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  var incy = 2;
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpsv('U','N','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpsv('U','N','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpsv('U','N','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpsv('U','N','U',3,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpsv('U','N','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpsv('U','N','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpsv('U','N','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpsv('U','N','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpsv('U','T','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpsv('U','T','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpsv('U','T','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpsv('U','T','U',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpsv('U','T','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpsv('U','T','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpsv('U','T','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpsv('U','T','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpsv('U','C','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpsv('U','C','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpsv('U','C','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpsv('U','C','U',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpsv('U','C','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpsv('U','C','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpsv('U','C','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpsv('U','C','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  k = 0;
  for ( j = 0; j < n; j ++ ) {
    for ( i = j; i < n; i ++ ) {
      P[ ioffP + k ] = new Complex();
      P[ ioffP + k++ ].value = A[ ioffA + i + j * lda ];
    }
  }
  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpsv('L','N','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpsv('L','N','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpsv('L','N','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpsv('L','N','U',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpsv('L','N','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpsv('L','N','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpsv('L','N','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpsv('L','N','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpsv('L','T','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpsv('L','T','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpsv('L','T','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpsv('L','T','U',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpsv('L','T','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpsv('L','T','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpsv('L','T','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpsv('L','T','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpsv('L','C','U',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
    "ztpsv('L','C','U',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpsv('L','C','U',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpsv('L','C','U',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";

  for ( i = 0; i < 3; i ++ ) {
    x[ ioffx + i ] = new Complex( i - 1 , 2 * i );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1 , i - 2 );
  }
  Blas2.ztpsv('L','C','N',n,P,x,incx,ioffP,ioffx);
  document.getElementById("debug_textarea").value +=
   "ztpsv('L','C','N',n,P,x,incx) , x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + "\n";
  Blas2.ztpsv('L','C','N',n,P,y,incy,ioffP,ioffy);
  document.getElementById("debug_textarea").value +=
   "ztpsv('L','C','N',n,P,y,incy) , y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2] + " " + y[ioffy + 3] + " "
        + y[ioffy + 4] + " " + y[ioffy + 5] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgemv() {
  document.getElementById("debug_textarea").value +=
   "testing dgemv *************"+ "\n";
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
  Blas2.dgemv("N",m,n,0.,A,lda,x,incx,1.,y,incy,ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dgemv(N,m,n,0.,A,lda,x,incx,1.,y,incy): y = " + y[ioffy + 0]
    + " " + y[ioffy + 2] + " " + y[ioffy + 4] + " " + y[ioffy + 6]+ "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dgemv("N",m,n,0.,A,lda,x,incx,2.,y,incy,ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dgemv(N,m,n,0.,A,lda,x,incx,2.,y,incy): y = " + y[ioffy + 0]
    + " " + y[ioffy + 2] + " " + y[ioffy + 4] + " " + y[ioffy + 6]+ "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dgemv("N",m,n,3.,A,lda,x,incx,1.,y,incy,ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dgemv(N,m,n,3.,A,lda,x,incx,1.,y,incy): y = " + y[ioffy + 0]
    + " " + y[ioffy + 2] + " " + y[ioffy + 4] + " " + y[ioffy + 6]+ "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dgemv("N",m,n,3.,A,lda,x,incx,2.,y,incy,ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dgemv(N,m,n,3.,A,lda,x,incx,2.,y,incy): y = " + y[ioffy + 0]
    + " " + y[ioffy + 2] + " " + y[ioffy + 4] + " " + y[ioffy + 6]+ "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dgemv("T",m,n,0.,A,lda,y,incy,1.,x,incx,ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "dgemv(T,m,n,0.,A,lda,y,incy,1.,x,incx): x = " + x[ioffx + 0]
    + " " + x[ioffx + 2] + " " + x[ioffx + 4]+ "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dgemv("T",m,n,0.,A,lda,y,incy,2.,x,incx,ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "dgemv(T,m,n,0.,A,lda,y,incy,2.,x,incx): x = " + x[ioffx + 0]
    + " " + x[ioffx + 2] + " " + x[ioffx + 4]+ "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dgemv("T",m,n,3.,A,lda,y,incy,1.,x,incx,ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "dgemv(T,m,n,3.,A,lda,y,incy,1.,x,incx): x = " + x[ioffx + 0]
    + " " + x[ioffx + 2] + " " + x[ioffx + 4]+ "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dgemv("T",m,n,3.,A,lda,y,incy,2.,x,incx,ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "dgemv(T,m,n,3.,A,lda,y,incy,2.,x,incx): x = " + x[ioffx + 0]
    + " " + x[ioffx + 2] + " " + x[ioffx + 4]+ "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zgemv() {
  document.getElementById("debug_textarea").value +=
   "testing zgemv *************"+ "\n";
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
      A[ ioffA + i + j * lda ] = new Complex( i + j , i * j );
    }
  }
  var incy = 2;
  var incx = 2;
  var zero = new Complex( 0., 0. );
  var one = new Complex( 1., 0. );
  var alpha = new Complex( 3., 4. );
  var beta = new Complex( 2., 5. );
  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2 , j + 2 );
  }
  Blas2.zgemv("N",m,n,zero,A,lda,x,incx,one,y,incy,ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zgemv(N,m,n,zero,A,lda,x,incx,one,y,incy): y = "
        + y[ioffy + 0].toString() + " " + y[ioffy + 2].toString() + " "
        + y[ioffy + 4].toString() + " " + y[ioffy + 6].toString()+ "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2 , j + 2 );
  }
  Blas2.zgemv("N",m,n,zero,A,lda,x,incx,beta,y,incy,ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zgemv(N,m,n,zero,A,lda,x,incx,beta,y,incy): y = "
        + y[ioffy + 0].toString() + " " + y[ioffy + 2].toString() + " "
        + y[ioffy + 4].toString() + " " + y[ioffy + 6].toString()+ "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2 , j + 2 );
  }
  Blas2.zgemv("N",m,n,alpha,A,lda,x,incx,one,y,incy,ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zgemv(N,m,n,alpha,A,lda,x,incx,one,y,incy): y = "
    + y[ioffy + 0].toString() + " " + y[ioffy + 2].toString()
    + " " + y[ioffy + 4].toString() + " " + y[ioffy + 6].toString()+ "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2 , j + 2 );
  }
  Blas2.zgemv("N",m,n,alpha,A,lda,x,incx,beta,y,incy,
    ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zgemv(N,m,n,alpha,A,lda,x,incx,beta,y,incy): y = "
        + y[ioffy + 0].toString() + " " + y[ioffy + 2].toString() + " "
        + y[ioffy + 4].toString() + " " + y[ioffy + 6].toString()+ "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2 , j + 2 );
  }
  Blas2.zgemv("T",m,n,zero,A,lda,y,incy,one,x,incx,ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "zgemv(T,m,n,zero,A,lda,y,incy,one,x,incx): x = "
        + x[ioffx + 0].toString() + " " + x[ioffx + 2].toString() + " "
        + x[ioffx + 4].toString()+ "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2 , j + 2 );
  }
  Blas2.zgemv("T",m,n,zero,A,lda,y,incy,beta,x,incx,ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "zgemv(T,m,n,zero,A,lda,y,incy,beta,x,incx): x = "
        + x[ioffx + 0].toString() + " " + x[ioffx + 2].toString() + " "
        + x[ioffx + 4].toString()+ "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
      for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2 , j + 2 );
  }
  Blas2.zgemv("T",m,n,alpha,A,lda,y,incy,one,x,incx,ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "zgemv(T,m,n,alpha,A,lda,y,incy,one,x,incx): x = "
        + x[ioffx + 0].toString() + " " + x[ioffx + 2].toString() + " "
        + x[ioffx + 4].toString()+ "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2 , j + 2 );
  }
  Blas2.zgemv("T",m,n,alpha,A,lda,y,incy,beta,x,incx,
    ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "zgemv(T,m,n,alpha,A,lda,y,incy,beta,x,incx): x = "
        + x[ioffx + 0].toString() + " " + x[ioffx + 2].toString() + " "
        + x[ioffx + 4].toString()+ "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2 , j + 2 );
  }
  Blas2.zgemv("C",m,n,zero,A,lda,y,incy,one,x,incx,ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "zgemv(C,m,n,zero,A,lda,y,incy,one,x,incx): x = "
        + x[ioffx + 0].toString() + " " + x[ioffx + 2].toString() + " "
        + x[ioffx + 4].toString()+ "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2 , j + 2 );
  }
  Blas2.zgemv("C",m,n,zero,A,lda,y,incy,beta,x,incx,ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "zgemv(C,m,n,zero,A,lda,y,incy,beta,x,incx): x = "
        + x[ioffx + 0].toString() + " " + x[ioffx + 2].toString() + " "
        + x[ioffx + 4].toString()+ "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2 , j + 2 );
  }
  Blas2.zgemv("C",m,n,alpha,A,lda,y,incy,one,x,incx,ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "zgemv(C,m,n,alpha,A,lda,y,incy,one,x,incx): x = "
        + x[ioffx + 0].toString() + " " + x[ioffx + 2].toString() + " "
        + x[ioffx + 4].toString()+ "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2 , j + 2 );
  }
  Blas2.zgemv("C",m,n,alpha,A,lda,y,incy,beta,x,incx,
    ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "zgemv(C,m,n,alpha,A,lda,y,incy,beta,x,incx): x = "
        + x[ioffx + 0].toString() + " " + x[ioffx + 2].toString() + " "
        + x[ioffx + 4].toString()+ "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgbmv() {
  document.getElementById("debug_textarea").value +=
   "testing dgbmv *************"+ "\n";
  var ioffA = 1;
  var ioffB = 2;
  var ioffx = 3;
  var ioffy = 4;
  var A = new Array( ioffA + 15 );
  var B = new Array( ioffB + 12 );
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 8 );
  var m = 4;
  var n = 3;
  var kl = 2;
  var ku = 1;
  var lda = 5;
  for (var j = 0; j < n; j ++ ) {
    for (var i=Math.max(0,j-ku);i<=Math.min(m-1,j+kl);i++) {
      A[ ioffA + i + j * lda ] = Number( i + j );
    }
  }
  var ldb = 4;
  for (j = 0; j < n; j ++ ) {
    var k=ku-j;
    for ( i = Math.max(0,j-ku); i <= Math.min(m-1,j+kl); i++ ) {
      B[ ioffB + k+ i + j * ldb ] = A[ ioffA + i + j * lda ];
    }
  }
  var incy = 2;
  var incx = 2;
  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dgbmv("N",m,n,kl,ku,0.,B,ldb,x,incx,1.,y,incy,
    ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dgbmv(N,m,n,kl,ku,0.,B,ldb,x,incx,1.,y,incy): y = "
    + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4] + " "
    + y[ioffy + 6] + "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dgbmv("N",m,n,kl,ku,0.,B,ldb,x,incx,2.,y,incy,
    ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dgbmv(N,m,n,kl,ku,0.,B,ldb,x,incx,2.,y,incy): y = "
    + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4] + " "
    + y[ioffy + 6] + "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dgbmv("N",m,n,kl,ku,3.,B,ldb,x,incx,1.,y,incy,
    ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dgbmv(N,m,n,kl,ku,3.,B,ldb,x,incx,1.,y,incy): y = "
    + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4] + " "
    + y[ioffy + 6] + "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dgbmv("N",m,n,kl,ku,3.,B,ldb,x,incx,2.,y,incy,
    ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dgbmv(N,m,n,kl,ku,3.,B,ldb,x,incx,2.,y,incy): y = "
    + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4] + " "
    + y[ioffy + 6] + "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dgbmv("T",m,n,kl,ku,0.,B,ldb,y,incy,1.,x,incx,
    ioffB,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "dgbmv(T,m,n,kl,ku,0.,B,ldb,y,incy,1.,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dgbmv("T",m,n,kl,ku,0.,B,ldb,y,incy,2.,x,incx,
    ioffB,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "dgbmv(T,m,n,kl,ku,0.,B,ldb,y,incy,2.,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dgbmv("T",m,n,kl,ku,3.,B,ldb,y,incy,1.,x,incx,
    ioffB,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "dgbmv(T,m,n,kl,ku,3.,B,ldb,y,incy,1.,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 8; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dgbmv("T",m,n,kl,ku,3.,B,ldb,y,incy,2.,x,incx,
    ioffB,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "dgbmv(T,m,n,kl,ku,3.,B,ldb,y,incy,2.,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zgbmv() {
  document.getElementById("debug_textarea").value +=
   "testing zgbmv *************" + "\n";
  var ioffA = 1;
  var ioffB = 2;
  var ioffx = 3;
  var ioffy = 4;
  var A = new Array( ioffA + 15 );
  var B = new Array( ioffB + 12 );
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 8 );
  var m = 4;
  var n = 3;
  var kl = 2;
  var ku = 1;
  var lda = 5;
  for (var j = 0; j < n; j ++ ) {
    for (var i=Math.max(0,j-ku);i<=Math.min(m-1,j+kl);i++) {
      A[ ioffA + i + j * lda ] = new Complex(i+j,i*j);
    }
  }
  var ldb = 4;
  for (j = 0; j < n; j ++ ) {
    var k=ku-j;
    for ( i = Math.max(0,j-ku); i <= Math.min(m-1,j+kl); i++ ) {
      B[ ioffB + k+ i + j * ldb ] = new Complex();
      B[ ioffB + k+ i + j * ldb ].value = A[ ioffA + i + j * lda ];
    }
  }
  var zero = new Complex( 0., 0. );
  var one = new Complex( 1., 0. );
  var alpha = new Complex( 3. , 4. );
  var beta = new Complex( 2. , 5. );
  var incy = 2;
  var incx = 2;
  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2  , j + 2 );
  }
  Blas2.zgbmv("N",m,n,kl,ku,zero,B,ldb,x,incx,one,y,incy,
    ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zgbmv(N,m,n,kl,ku,zero,B,ldb,x,incx,one,y,incy): y = "
    + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4] + " "
    + y[ioffy + 6] + "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2  , j + 2 );
  }
  Blas2.zgbmv("N",m,n,kl,ku,zero,B,ldb,x,incx,beta,y,incy,
    ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zgbmv(N,m,n,kl,ku,zero,B,ldb,x,incx,beta,y,incy): y = "
    + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4] + " "
    + y[ioffy + 6] + "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2  , j + 2 );
  }
  Blas2.zgbmv("N",m,n,kl,ku,alpha,B,ldb,x,incx,one,y,incy,
    ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zgbmv(N,m,n,kl,ku,alpha,B,ldb,x,incx,one,y,incy): y = "
    + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4] + " "
    + y[ioffy + 6] + "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2  , j + 2 );
  }
  Blas2.zgbmv("N",m,n,kl,ku,alpha,B,ldb,x,incx,beta,y,incy,
    ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zgbmv(N,m,n,kl,ku,alpha,B,ldb,x,incx,beta,y,incy): y = "
    + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4] + " "
    + y[ioffy + 6] + "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2  , j + 2 );
  }
  Blas2.zgbmv("T",m,n,kl,ku,zero,B,ldb,y,incy,one,x,incx,
    ioffB,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "zgbmv(T,m,n,kl,ku,zero,B,ldb,y,incy,one,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2  , j + 2 );
  }
  Blas2.zgbmv("T",m,n,kl,ku,zero,B,ldb,y,incy,beta,x,incx,
    ioffB,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "zgbmv(T,m,n,kl,ku,zero,B,ldb,y,incy,beta,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2  , j + 2 );
  }
  Blas2.zgbmv("T",m,n,kl,ku,alpha,B,ldb,y,incy,one,x,incx,
    ioffB,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "zgbmv(T,m,n,kl,ku,alpha,B,ldb,y,incy,one,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2  , j + 2 );
  }
  Blas2.zgbmv("T",m,n,kl,ku,alpha,B,ldb,y,incy,beta,x,incx,
    ioffB,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "zgbmv(T,m,n,kl,ku,alpha,B,ldb,y,incy,beta,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2  , j + 2 );
  }
  Blas2.zgbmv("C",m,n,kl,ku,zero,B,ldb,y,incy,one,x,incx,
    ioffB,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "zgbmv(C,m,n,kl,ku,zero,B,ldb,y,incy,one,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2  , j + 2 );
  }
  Blas2.zgbmv("C",m,n,kl,ku,zero,B,ldb,y,incy,beta,x,incx,
    ioffB,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "zgbmv(C,m,n,kl,ku,zero,B,ldb,y,incy,beta,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2  , j + 2 );
  }
  Blas2.zgbmv("C",m,n,kl,ku,alpha,B,ldb,y,incy,one,x,incx,
    ioffB,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "zgbmv(C,m,n,kl,ku,alpha,B,ldb,y,incy,one,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 8; i ++ ) {
    y[ioffy + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2  , j + 2 );
  }
  Blas2.zgbmv("C",m,n,kl,ku,alpha,B,ldb,y,incy,beta,x,incx,
    ioffB,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "zgbmv(C,m,n,kl,ku,alpha,B,ldb,y,incy,beta,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dsymv() {
  document.getElementById("debug_textarea").value +=
   "testing dsymv *************" + "\n";
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
  Blas2.dsymv("U",n,0.,A,lda,x,incx,1.,y,incy,
    ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dsymv(U,n,0.,A,lda,x,incx,1.,y,incy): y = " + y[ioffy + 0]
    + " " + y[ioffy + 2] + " " + y[ioffy + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dsymv("U",n,0.,A,lda,x,incx,2.,y,incy,
    ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dsymv(U,n,0.,A,lda,x,incx,2.,y,incy): y = " + y[ioffy + 0]
    + " " + y[ioffy + 2] + " " + y[ioffy + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dsymv("U",n,3.,A,lda,x,incx,1.,y,incy,
    ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dsymv(U,n,3.,A,lda,x,incx,1.,y,incy): y = " + y[ioffy + 0]
    + " " + y[ioffy + 2] + " " + y[ioffy + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dsymv("U",n,3.,A,lda,x,incx,2.,y,incy,
    ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dsymv(U,n,3.,A,lda,x,incx,2.,y,incy): y = " + y[ioffy + 0]
    + " " + y[ioffy + 2] + " " + y[ioffy + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dsymv("L",n,0.,A,lda,y,incy,1.,x,incx,
    ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "dsymv(L,n,0.,A,lda,y,incy,1.,x,incx): x = " + x[ioffx + 0]
    + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dsymv("L",n,0.,A,lda,y,incy,2.,x,incx,
    ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "dsymv(L,n,0.,A,lda,y,incy,2.,x,incx): x = " + x[ioffx + 0]
    + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dsymv("L",n,3.,A,lda,y,incy,1.,x,incx,
    ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "dsymv(L,n,3.,A,lda,y,incy,1.,x,incx): x = " + x[ioffx + 0]
    + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dsymv("L",n,3.,A,lda,y,incy,2.,x,incx,
    ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "dsymv(L,n,3.,A,lda,y,incy,2.,x,incx): x = " + x[ioffx + 0]
    + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zhemv() {
  document.getElementById("debug_textarea").value +=
   "testing zhemv *************" + "\n";
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
      A[ ioffA + i + j * lda ] = new Complex( i + j*2 , i*3+j);
    }
  }
  var zero = new Complex(0.,0.);
  var one = new Complex(1.,0.);
  var alpha = new Complex(3.,4.);
  var beta = new Complex(2.,5.);
  var incy = 2;
  var incx = 2;

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = new Complex( i + 1 , 1-i);
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = new Complex( j * 2 , j+2);
  Blas2.zhemv("U",n,zero,A,lda,x,incx,one,y,incy,
    ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zhemv(U,n,zero,A,lda,x,incx,one,y,incy): y = " + y[ioffy + 0]
    + " " + y[ioffy + 2] + " " + y[ioffy + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = new Complex( i + 1 , 1-i);
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = new Complex( j * 2 , j+2);
  Blas2.zhemv("U",n,zero,A,lda,x,incx,beta,y,incy,
    ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zhemv(U,n,zero,A,lda,x,incx,beta,y,incy): y = "
    + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = new Complex( i + 1 , 1-i);
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = new Complex( j * 2 , j+2);
  Blas2.zhemv("U",n,alpha,A,lda,x,incx,one,y,incy,
    ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zhemv(U,n,alpha,A,lda,x,incx,one,y,incy): y = "
    + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = new Complex( i + 1 , 1-i);
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = new Complex( j * 2 , j+2);
  Blas2.zhemv("U",n,alpha,A,lda,x,incx,beta,y,incy,
    ioffA,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zhemv(U,n,alpha,A,lda,x,incx,beta,y,incy): y = "
    + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = new Complex( i + 1 , 1-i);
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = new Complex( j * 2 , j+2);
  Blas2.zhemv("L",n,zero,A,lda,y,incy,one,x,incx,
    ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "zhemv(L,n,zero,A,lda,y,incy,one,x,incx): x = " + x[ioffx + 0]
    + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = new Complex( i + 1 , 1-i);
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = new Complex( j * 2 , j+2);
  Blas2.zhemv("L",n,zero,A,lda,y,incy,beta,x,incx,
    ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "zhemv(L,n,zero,A,lda,y,incy,beta,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = new Complex( i + 1 , 1-i);
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = new Complex( j * 2 , j+2);
  Blas2.zhemv("L",n,alpha,A,lda,y,incy,one,x,incx,
    ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "zhemv(L,n,alpha,A,lda,y,incy,one,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";

  for ( i = 0; i < 6; i ++ ) y[ioffy + i] = new Complex( i + 1 , 1-i);
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = new Complex( j * 2 , j+2);
  Blas2.zhemv("L",n,alpha,A,lda,y,incy,beta,x,incx,
    ioffA,ioffy,ioffx);
  document.getElementById("debug_textarea").value +=
   "zhemv(L,n,alpha,A,lda,y,incy,beta,x,incx): x = "
    + x[ioffx + 0] + " " + x[ioffx + 2] + " " + x[ioffx + 4] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dsbmv() {
  document.getElementById("debug_textarea").value +=
   "testing dsbmv *************" + "\n";
  var ioffA = 1;
  var ioffB = 2;
  var ioffx = 3;
  var ioffy = 4;
  var A = new Array( ioffA + 12 );
  var B = new Array( ioffB + 9 );
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 6 );
  var n = 3;
  var lda = 4;
  var i = -1;
  var j = -1;
  for (j=0;j<n;j++) {
    for (i=0;i<n;i++) {
      A[ioffA + i+j*lda]=Number(i+j*2);
    }
  }
  var k = 1;
  var ldb = 3;
  for (j=0;j<n;j++) {
    var m=k-j
    for (i=Math.max(0,j-k);i<=j;i++) {
      B[ioffB + m+i+j*ldb]=A[ioffA + i+j*lda];
    }
  }
  var incx = 2;
  for ( j = 0; j < 6; j ++ ) { x[ ioffx + j ] = Number( j*2 ); }
  var incy = 2;
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dsbmv('U',n,k,0.,B,ldb,x,incx,1.,y,incy,ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dsbmv('U',n,k,0.,B,ldb,x,incx,1.,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dsbmv('U',n,k,0.,B,ldb,x,incx,2.,y,incy,ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dsbmv('U',n,k,0.,B,ldb,x,incx,2.,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dsbmv('U',n,k,3.,B,ldb,x,incx,1.,y,incy,ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dsbmv('U',n,k,3.,B,ldb,x,incx,1.,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dsbmv('U',n,k,3.,B,ldb,x,incx,2.,y,incy,ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dsbmv('U',n,k,3.,B,ldb,x,incx,2.,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    m = - j;
    for ( i = j; i <= Math.min( n - 1 , j + k ); i ++ ) {
      B[ ioffB + m + i + j * ldb ] = A[ ioffA + i + j * lda ];
    }
  }
  for ( j = 0; j < 6; j ++ ) { x[ ioffx + j ] = Number( j*2 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dsbmv('L',n,k,0.,B,ldb,x,incx,1.,y,incy,ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dsbmv('L',n,k,0.,B,ldb,x,incx,1.,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dsbmv('L',n,k,0.,B,ldb,x,incx,2.,y,incy,ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dsbmv('L',n,k,0.,B,ldb,x,incx,2.,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dsbmv('L',n,k,3.,B,ldb,x,incx,1.,y,incy,ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dsbmv('L',n,k,3.,B,ldb,x,incx,1.,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dsbmv('L',n,k,3.,B,ldb,x,incx,2.,y,incy,ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dsbmv('L',n,k,3.,B,ldb,x,incx,2.,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zhbmv() {
  document.getElementById("debug_textarea").value +=
   "testing zhbmv *************" + "\n";
  var ioffA = 1;
  var ioffB = 2;
  var ioffx = 3;
  var ioffy = 4;
  var A = new Array( ioffA + 12 );
  var B = new Array( ioffB + 9 );
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 6 );
  var n = 3;
  var lda = 4;
  var i = -1;
  var j = -1;
  for (j=0;j<n;j++) {
    for (i=0;i<n;i++) {
      A[ioffA + i+j*lda]=new Complex(i+j*2,i*3+j);
    }
  }
  var k = 1;
  var ldb = 3;
  for (j=0;j<n;j++) {
    var m=k-j
    for (i=Math.max(0,j-k);i<=j;i++) {
      B[ioffB + m+i+j*ldb]=new Complex();
      B[ioffB + m+i+j*ldb].value=A[ioffA + i+j*lda];
    }
  }
  var zero=new Complex(0.,0.);
  var one=new Complex(1.,0.);
  var alpha=new Complex(3.,4.);
  var beta=new Complex(2.,5.);
  var incx = 2;
  for ( j = 0; j < 6; j ++ ) {
    x[ ioffx + j ] = new Complex( j*2 ,j+2);
  }
  var incy = 2;
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1,1-3*i );
  }
  Blas2.zhbmv('U',n,k,zero,B,ldb,x,incx,one,y,incy,ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zhbmv('U',n,k,zero,B,ldb,x,incx,one,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1,1-3*i );
  }
  Blas2.zhbmv('U',n,k,zero,B,ldb,x,incx,beta,y,incy,ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zhbmv('U',n,k,zero,B,ldb,x,incx,beta,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1,1-3*i );
  }
  Blas2.zhbmv('U',n,k,alpha,B,ldb,x,incx,one,y,incy,ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zhbmv('U',n,k,alpha,B,ldb,x,incx,one,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1,1-3*i );
  }
  Blas2.zhbmv('U',n,k,alpha,B,ldb,x,incx,beta,y,incy,
    ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zhbmv('U',n,k,alpha,B,ldb,x,incx,beta,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( j = 0; j < n; j ++ ) {
    m = - j;
    for ( i = j; i <= Math.min( n - 1 , j + k ); i ++ ) {
      B[ioffB + m+i+j*ldb]=new Complex();
      B[ ioffB + m + i + j * ldb ].value = A[ ioffA + i + j * lda ];
    }
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ ioffx + j ] = new Complex( j*2 ,j+2);
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1,1-3*i );
  }
  Blas2.zhbmv('L',n,k,zero,B,ldb,x,incx,one,y,incy,ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zhbmv('L',n,k,zero,B,ldb,x,incx,one,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1,1-3*i );
  }
  Blas2.zhbmv('L',n,k,zero,B,ldb,x,incx,beta,y,incy,ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zhbmv('L',n,k,zero,B,ldb,x,incx,beta,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1,1-3*i );
  }
  Blas2.zhbmv('L',n,k,alpha,B,ldb,x,incx,one,y,incy,ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zhbmv('L',n,k,alpha,B,ldb,x,incx,one,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1,1-3*i );
  }
  Blas2.zhbmv('L',n,k,alpha,B,ldb,x,incx,beta,y,incy,
    ioffB,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zhbmv('L',n,k,alpha,B,ldb,x,incx,beta,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dspmv() {
  document.getElementById("debug_textarea").value +=
   "testing dspmv *************" + "\n";
  var ioffA = 1;
  var ioffP = 2;
  var ioffx = 3;
  var ioffy = 4;
  var A = new Array( ioffA + 12 );
  var P = new Array( ioffP + 6 );
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 6 );
  var n = 3;
  var lda = 4;
  var i = -1;
  var j = -1;
  for (j=0;j<n;j++) {
    for (i=0;i<n;i++) {
      A[ioffA + i+j*lda]=Number(i+j*2);
    }
  }
  var k = 0;
  for (j=0;j<n;j++) {
    for (i=0;i<=j;i++) {
      P[ioffP + k++]=A[ioffA + i+j*lda];
    }
  }
  var incx = 2;
  for ( j = 0; j < 6; j ++ ) { x[ ioffx + j ] = Number( j*2 ); }
  var incy = 2;
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dspmv('U',n,0.,P,x,incx,1.,y,incy,ioffP,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dspmv('U',n,0.,P,x,incx,1.,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dspmv('U',n,0.,P,x,incx,2.,y,incy,ioffP,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dspmv('U',n,0.,P,x,incx,2.,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dspmv('U',n,3.,P,x,incx,1.,y,incy,ioffP,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dspmv('U',n,3.,P,x,incx,1.,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dspmv('U',n,3.,P,x,incx,2.,y,incy,ioffP,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dspmv('U',n,3.,P,x,incx,2.,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  k=0;
  for (j=0;j<n;j++) {
    for (i=j;i<n;i++) {
      P[ioffP + k++]=A[ioffA + i+j*lda];
    }
  }
  for ( j = 0; j < 6; j ++ ) { x[ ioffx + j ] = Number( j*2 ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dspmv('L',n,0.,P,x,incx,1.,y,incy,ioffP,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dspmv('L',n,0.,P,x,incx,1.,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dspmv('L',n,0.,P,x,incx,2.,y,incy,ioffP,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dspmv('L',n,0.,P,x,incx,2.,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dspmv('L',n,3.,P,x,incx,1.,y,incy,ioffP,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dspmv('L',n,3.,P,x,incx,1.,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( i + 1 ); }
  Blas2.dspmv('L',n,3.,P,x,incx,2.,y,incy,ioffP,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "dspmv('L',n,3.,P,x,incx,2.,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zhpmv() {
  document.getElementById("debug_textarea").value +=
   "testing zhpmv *************" + "\n";
  var ioffA = 1;
  var ioffP = 2;
  var ioffx = 3;
  var ioffy = 4;
  var A = new Array( ioffA + 12 );
  var P = new Array( ioffP + 6 );
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 6 );
  var n = 3;
  var lda = 4;
  var i = -1;
  var j = -1;
  for (j=0;j<n;j++) {
    for (i=0;i<n;i++) {
      A[ioffA + i+j*lda]=new Complex(i+j*2,i*3+j);
    }
  }
  var k = 0;
  for (j=0;j<n;j++) {
    for (i=0;i<=j;i++) {
      P[ioffP + k]=new Complex();
      P[ioffP + k++].value=A[ioffA + i+j*lda];
    }
  }
  var zero=new Complex(0.,0.);
  var one=new Complex(1.,0.);
  var alpha=new Complex(3.,4.);
  var beta=new Complex(2.,5.);
  var incx = 2;
  for ( j = 0; j < 6; j ++ ) {
    x[ ioffx + j ] = new Complex( j*2 ,j+2);
  }
  var incy = 2;
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1,1-3*i );
  }
  Blas2.zhpmv('U',n,zero,P,x,incx,one,y,incy,ioffP,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zhpmv('U',n,zero,P,x,incx,one,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1,1-3*i );
  }
  Blas2.zhpmv('U',n,zero,P,x,incx,beta,y,incy,ioffP,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zhpmv('U',n,zero,P,x,incx,beta,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1,1-3*i );
  }
  Blas2.zhpmv('U',n,alpha,P,x,incx,one,y,incy,ioffP,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zhpmv('U',n,alpha,P,x,incx,one,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1,1-3*i );
  }
  Blas2.zhpmv('U',n,alpha,P,x,incx,beta,y,incy,ioffP,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zhpmv('U',n,alpha,P,x,incx,beta,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  k=0;
  for ( j = 0; j < n; j ++ ) {
    for ( i = j; i < n; i ++ ) {
      P[ioffP + k]=new Complex();
      P[ioffP + k++].value=A[ioffA + i+j*lda];
    }
  }
  for ( j = 0; j < 6; j ++ ) {
    x[ ioffx + j ] = new Complex( j*2 ,j+2);
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1,1-3*i );
  }
  Blas2.zhpmv('L',n,zero,P,x,incx,one,y,incy,ioffP,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zhpmv('L',n,zero,P,x,incx,one,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1,1-3*i );
  }
  Blas2.zhpmv('L',n,zero,P,x,incx,beta,y,incy,ioffP,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zhpmv('L',n,zero,P,x,incx,beta,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1,1-3*i );
  }
  Blas2.zhpmv('L',n,alpha,P,x,incx,one,y,incy,ioffP,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zhpmv('L',n,alpha,P,x,incx,one,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";

  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( i + 1,1-3*i );
  }
  Blas2.zhpmv('L',n,alpha,P,x,incx,beta,y,incy,ioffP,ioffx,ioffy);
  document.getElementById("debug_textarea").value +=
   "zhpmv('L',n,alpha,P,x,incx,beta,y,incy) , y = "
        + y[ioffy + 0] + " " + y[ioffy + 2] + " " + y[ioffy + 4]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dger() {
  document.getElementById("debug_textarea").value +=
   "testing dger *************" + "\n";
  var ioffA = 1;
  var ioffx = 2;
  var ioffy = 3;
  var A = new Array( ioffA + 15 );
  var x = new Array( ioffx + 8 );
  var y = new Array( ioffy + 6 );
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
  for ( i = 0; i < 8; i ++ ) x[ioffx + i] = Number( i + 1 );
  for ( j = 0; j < 6; j ++ ) y[ioffy + j] = Number( j * 2 );
  Blas2.dger(m,n,0.,x,incx,y,incy,A,lda,ioffx,ioffy,ioffA);
  document.getElementById("debug_textarea").value +=
   "dger(m,n,0.,x,incx,y,incy,A,lda): A = " + A[ioffA + 0] + " "
    + A[ioffA + 1] + " " + A[ioffA + 2] + " " + A[ioffA + 3] + " "
    + A[ioffA + 5] + " " + A[ioffA + 6] + " " + A[ioffA + 7] + " "
    + A[ioffA + 8] + " " + A[ioffA + 10] + " " + A[ioffA + 11] + " "
    + A[ioffA + 12] + " " + A[ioffA + 13] + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffA + i + j * lda ] = Number( i + j );
    }
  }
  Blas2.dger(m,n,2.,x,incx,y,incy,A,lda,ioffx,ioffy,ioffA);
  document.getElementById("debug_textarea").value +=
   "dger(m,n,2.,x,incx,y,incy,A,lda): A = " + A[ioffA + 0] + " "
    + A[ioffA + 1] + " " + A[ioffA + 2] + " " + A[ioffA + 3] + " "
    + A[ioffA + 5] + " " + A[ioffA + 6] + " " + A[ioffA + 7] + " "
    + A[ioffA + 8] + " " + A[ioffA + 10] + " " + A[ioffA + 11] + " "
    + A[ioffA + 12] + " " + A[ioffA + 13] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zgerc() {
  document.getElementById("debug_textarea").value +=
   "testing zgerc *************" + "\n";
  var ioffA = 1;
  var ioffx = 2;
  var ioffy = 3;
  var A = new Array( ioffA + 15 );
  var x = new Array( ioffx + 8 );
  var y = new Array( ioffy + 6 );
  var m = 4;
  var n = 3;
  var lda = 5;
  for (var j = 0; j < n; j ++ ) {
    for (var i = 0; i < m; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j , i * j );
    }
  }
  var incy = 2;
  var incx = 2;
  for ( i = 0; i < 8; i ++ ) {
    x[ioffx + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    y[ioffy + j] = new Complex( j * 2 , j + 2 );
  }
  var zero = new Complex(0.,0.);
  Blas2.zgerc(m,n,zero,x,incx,y,incy,A,lda,ioffx,ioffy,ioffA);
  document.getElementById("debug_textarea").value +=
   "zgerc(m,n,zero,x,incx,y,incy,A,lda): A = "
    + A[ioffA + 0].toString() + " " + A[ioffA + 1].toString() + " "
    + A[ioffA + 2].toString() + " " + A[ioffA + 3].toString() + " "
    + A[ioffA + 5].toString() + " " + A[ioffA + 6].toString() + " "
    + A[ioffA + 7].toString() + " " + A[ioffA + 8].toString() + " "
    + A[ioffA + 10].toString() + " " + A[ioffA + 11].toString() + " "
    + A[ioffA + 12].toString() + " " + A[ioffA + 13].toString() + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j , i * j );
    }
  }
  var alpha = new Complex(2.,3.);
  Blas2.zgerc(m,n,alpha,x,incx,y,incy,A,lda,ioffx,ioffy,ioffA);
  document.getElementById("debug_textarea").value +=
   "zgerc(m,n,alpha,x,incx,y,incy,A,lda): A = "
    + A[ioffA + 0].toString() + " " + A[ioffA + 1].toString() + " "
    + A[ioffA + 2].toString() + " " + A[ioffA + 3].toString() + " "
    + A[ioffA + 5].toString() + " " + A[ioffA + 6].toString() + " "
    + A[ioffA + 7].toString() + " " + A[ioffA + 8].toString() + " "
    + A[ioffA + 10].toString() + " " + A[ioffA + 11].toString() + " "
    + A[ioffA + 12].toString() + " " + A[ioffA + 13].toString() + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zgeru() {
  document.getElementById("debug_textarea").value +=
   "testing zgeru *************" + "\n";
  var ioffA = 1;
  var ioffx = 2;
  var ioffy = 3;
  var A = new Array( ioffA + 15 );
  var x = new Array( ioffx + 8 );
  var y = new Array( ioffy + 6 );
  var m = 4;
  var n = 3;
  var lda = 5;
  for (var j = 0; j < n; j ++ ) {
    for (var i = 0; i < m; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j , i * j );
    }
  }
  var incy = 2;
  var incx = 2;
  for ( i = 0; i < 8; i ++ ) {
    x[ioffx + i] = new Complex( i + 1 , 1 - i );
  }
  for ( j = 0; j < 6; j ++ ) {
    y[ioffy + j] = new Complex( j * 2 , j + 2 );
  }
  var zero = new Complex(0.,0.);
  Blas2.zgeru(m,n,zero,x,incx,y,incy,A,lda,ioffx,ioffy,ioffA);
  document.getElementById("debug_textarea").value +=
   "zgeru(m,n,zero,x,incx,y,incy,A,lda): A = "
    + A[ioffA + 0].toString() + " " + A[ioffA + 1].toString() + " "
    + A[ioffA + 2].toString() + " " + A[ioffA + 3].toString() + " "
    + A[ioffA + 5].toString() + " " + A[ioffA + 6].toString() + " "
    + A[ioffA + 7].toString() + " " + A[ioffA + 8].toString() + " "
    + A[ioffA + 10].toString() + " " + A[ioffA + 11].toString() + " "
    + A[ioffA + 12].toString() + " " + A[ioffA + 13].toString() + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j , i * j );
    }
  }
  var alpha = new Complex(2.,3.);
  Blas2.zgeru(m,n,alpha,x,incx,y,incy,A,lda,ioffx,ioffy,ioffA);
  document.getElementById("debug_textarea").value +=
   "zgeru(m,n,alpha,x,incx,y,incy,A,lda): A = "
    + A[ioffA + 0].toString() + " " + A[ioffA + 1].toString() + " "
    + A[ioffA + 2].toString() + " " + A[ioffA + 3].toString() + " "
    + A[ioffA + 5].toString() + " " + A[ioffA + 6].toString() + " "
    + A[ioffA + 7].toString() + " " + A[ioffA + 8].toString() + " "
    + A[ioffA + 10].toString() + " " + A[ioffA + 11].toString() + " "
    + A[ioffA + 12].toString() + " " + A[ioffA + 13].toString() + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dsyr() {
  document.getElementById("debug_textarea").value +=
   "testing dsyr *************" + "\n";
  var ioffA = 1;
  var ioffx = 2;
  var A = new Array( ioffA + 12 );
  var x = new Array( ioffx + 6 );
  var n = 3;
  var lda = 4;
  for (var j = 0; j < n; j ++ ) {
    for (var i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = Number( i + j );
    }
  }
  var incx = 2;
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dsyr("U",n,0.,x,incx,A,lda,ioffx,ioffA);
  document.getElementById("debug_textarea").value +=
   "dsyr(U,n,0.,x,incx,A,lda): A = " + A[ioffA + 0] + " "
    + A[ioffA + 1] + " " + A[ioffA + 2] + " " + A[ioffA + 4] + " "
    + A[ioffA + 5] + " " + A[ioffA + 6] + " " + A[ioffA + 8] + " "
    + A[ioffA + 9] + " " + A[ioffA + 10] + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = Number( i + j );
    }
  }
  Blas2.dsyr("U",n,2.,x,incx,A,lda,ioffx,ioffA);
  document.getElementById("debug_textarea").value +=
   "dsyr(U,n,2.,x,incx,A,lda): A = " + A[ioffA + 0] + " "
    + A[ioffA + 1] + " " + A[ioffA + 2] + " " + A[ioffA + 4] + " "
    + A[ioffA + 5] + " " + A[ioffA + 6] + " " + A[ioffA + 8] + " "
    + A[ioffA + 9] + " " + A[ioffA + 10] + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = Number( i + j );
    }
  }
  Blas2.dsyr("L",n,2.,x,incx,A,lda,ioffx,ioffA);
  document.getElementById("debug_textarea").value +=
   "dsyr(L,n,2.,x,incx,A,lda): A = " + A[ioffA + 0] + " "
    + A[ioffA + 1] + " " + A[ioffA + 2] + " " + A[ioffA + 4] + " "
    + A[ioffA + 5] + " " + A[ioffA + 6] + " " + A[ioffA + 8] + " "
    + A[ioffA + 9] + " " + A[ioffA + 10] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zher() {
  document.getElementById("debug_textarea").value +=
   "testing zher *************" + "\n";
  var ioffA = 1;
  var ioffx = 2;
  var A = new Array( ioffA + 12 );
  var x = new Array( ioffx + 6 );
  var n = 3;
  var lda = 4;
  for (var j = 0; j < n; j ++ ) {
    for (var i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j , i * j );
    }
  }
  var incx = 2;
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2 , j + 2 );
  }
  Blas2.zher("U",n,0.,x,incx,A,lda,ioffx,ioffA);
  document.getElementById("debug_textarea").value +=
   "zher(U,n,0.,x,incx,A,lda): A = " + A[ioffA + 0] + " "
    + A[ioffA + 1] + " " + A[ioffA + 2] + " " + A[ioffA + 4] + " "
    + A[ioffA + 5] + " " + A[ioffA + 6] + " " + A[ioffA + 8] + " "
    + A[ioffA + 9] + " " + A[ioffA + 10] + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j , i * j );
    }
  }
  Blas2.zher("U",n,2.,x,incx,A,lda,ioffx,ioffA);
  document.getElementById("debug_textarea").value +=
   "zher(U,n,2.,x,incx,A,lda): A = " + A[ioffA + 0] + " "
    + A[ioffA + 1] + " " + A[ioffA + 2] + " " + A[ioffA + 4] + " "
    + A[ioffA + 5] + " " + A[ioffA + 6] + " " + A[ioffA + 8] + " "
    + A[ioffA + 9] + " " + A[ioffA + 10] + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j , i * j );
    }
  }
  Blas2.zher("L",n,2.,x,incx,A,lda,ioffx,ioffA);
  document.getElementById("debug_textarea").value +=
   "zher(L,n,2.,x,incx,A,lda): A = " + A[ioffA + 0] + " "
    + A[ioffA + 1] + " " + A[ioffA + 2] + " " + A[ioffA + 4] + " "
    + A[ioffA + 5] + " " + A[ioffA + 6] + " " + A[ioffA + 8] + " "
    + A[ioffA + 9] + " " + A[ioffA + 10] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dspr() {
  document.getElementById("debug_textarea").value +=
   "testing dspr *************" + "\n";
  var ioffA = 1;
  var ioffP = 2;
  var ioffx = 3;
  var A = new Array( ioffA + 12 );
  var P = new Array( ioffP + 6 );
  var x = new Array( ioffx + 6 );
  var n = 3;
  var lda = 4;
  for (var j = 0; j < n; j ++ ) {
    for (var i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = Number( i + j );
    }
  }
  var k=0;
  for (j=0;j<n;j++) {
    for (i=0;i<=j;i++) {
      P[ioffP + k++]=A[ioffA + i+j*lda];
    }
  }
  var incx = 2;
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  Blas2.dspr("U",n,0.,x,incx,P,ioffx,ioffP);
  document.getElementById("debug_textarea").value +=
   "dspr(U,n,0.,x,incx,P): P = " + P[ioffP + 0] + " "
    + P[ioffP + 1] + " " + P[ioffP + 2] + " " + P[ioffP + 3] + " "
    + P[ioffP + 4] + " " + P[ioffP + 5] + "\n";

  k=0;
  for (j=0;j<n;j++) {
    for (i=0;i<=j;i++) {
      P[ioffP + k++]=A[ioffA + i+j*lda];
    }
  }
  Blas2.dspr("U",n,2.,x,incx,P,ioffx,ioffP);
  document.getElementById("debug_textarea").value +=
   "dspr(U,n,2.,x,incx,P): P = " + P[ioffP + 0] + " "
    + P[ioffP + 1] + " " + P[ioffP + 2] + " " + P[ioffP + 3] + " "
    + P[ioffP + 4] + " " + P[ioffP + 5] + "\n";

  k=0;
  for (j=0;j<n;j++) {
    for (i=j;i<n;i++) {
      P[ioffP + k++]=A[ioffA + i+j*lda];
    }
  }
  Blas2.dspr("L",n,2.,x,incx,P,ioffx,ioffP);
  document.getElementById("debug_textarea").value +=
   "dspr(L,n,2.,x,incx,P): P = " + P[ioffP + 0] + " "
    + P[ioffP + 1] + " " + P[ioffP + 2] + " " + P[ioffP + 3] + " "
    + P[ioffP + 4] + " " + P[ioffP + 5] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zhpr() {
  document.getElementById("debug_textarea").value +=
   "testing zhpr *************" + "\n";
  var ioffA = 1;
  var ioffP = 2;
  var ioffx = 3;
  var A = new Array( ioffA + 12 );
  var P = new Array( ioffP + 6 );
  var x = new Array( ioffx + 6 );
  var n = 3;
  var lda = 4;
  for (var j = 0; j < n; j ++ ) {
    for (var i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j , i * j );
    }
  }
  var k=0;
  for (j=0;j<n;j++) {
    for (i=0;i<=j;i++) {
      P[ioffP + k] = new Complex();
      P[ioffP + k++].value=A[ioffA + i+j*lda];
    }
  }
  var incx = 2;
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2 , j + 2 );
  }
  Blas2.zhpr("U",n,0.,x,incx,P,ioffx,ioffP);
  document.getElementById("debug_textarea").value +=
   "zhpr(U,n,0.,x,incx,P): P = " + P[ioffP + 0] + " "
    + P[ioffP + 1] + " " + P[ioffP + 2] + " " + P[ioffP + 3] + " "
    + P[ioffP + 4] + " " + P[ioffP + 5] + "\n";

  k=0;
  for (j=0;j<n;j++) {
    for (i=0;i<=j;i++) {
      P[ioffP + k] = new Complex();
      P[ioffP + k++].value=A[ioffA + i+j*lda];
    }
  }
  Blas2.zhpr("U",n,2.,x,incx,P,ioffx,ioffP);
  document.getElementById("debug_textarea").value +=
   "zhpr(U,n,2.,x,incx,P): P = " + P[ioffP + 0] + " "
    + P[ioffP + 1] + " " + P[ioffP + 2] + " " + P[ioffP + 3] + " "
    + P[ioffP + 4] + " " + P[ioffP + 5] + "\n";

  k=0;
  for (j=0;j<n;j++) {
    for (i=j;i<n;i++) {
      P[ioffP + k] = new Complex();
      P[ioffP + k++].value=A[ioffA + i+j*lda];
    }
  }
  Blas2.zhpr("L",n,2.,x,incx,P,ioffx,ioffP);
  document.getElementById("debug_textarea").value +=
   "zhpr(L,n,2.,x,incx,P): P = " + P[ioffP + 0] + " "
    + P[ioffP + 1] + " " + P[ioffP + 2] + " " + P[ioffP + 3] + " "
    + P[ioffP + 4] + " " + P[ioffP + 5] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dsyr2() {
  document.getElementById("debug_textarea").value +=
   "testing dsyr2 *************" + "\n";
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
      A[ ioffA + i + j * lda ] = Number( i + j );
    }
  }
  var incx = 2;
  var incy = 2;
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  for ( j = 0; j < 6; j ++ ) y[ioffy + j] = Number( j + 1 );
  Blas2.dsyr2("U",n,0.,x,incx,y,incy,A,lda,ioffx,ioffy,ioffA);
  document.getElementById("debug_textarea").value +=
   "dsyr2(U,n,0.,x,incx,y,incy,A,lda): A = " + A[ioffA + 0] + " "
    + A[ioffA + 1] + " " + A[ioffA + 2] + " " + A[ioffA + 4] + " "
    + A[ioffA + 5] + " " + A[ioffA + 6] + " " + A[ioffA + 8] + " "
    + A[ioffA + 9] + " " + A[ioffA + 10] + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = Number( i + j );
    }
  }
  Blas2.dsyr2("U",n,2.,x,incx,y,incy,A,lda,ioffx,ioffy,ioffA);
  document.getElementById("debug_textarea").value +=
   "dsyr2(U,n,2.,x,incx,y,incy,A,lda): A = " + A[ioffA + 0] + " "
    + A[ioffA + 1] + " " + A[ioffA + 2] + " " + A[ioffA + 4] + " "
    + A[ioffA + 5] + " " + A[ioffA + 6] + " " + A[ioffA + 8] + " "
    + A[ioffA + 9] + " " + A[ioffA + 10] + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = Number( i + j );
    }
  }
  Blas2.dsyr2("L",n,2.,x,incx,y,incy,A,lda,ioffx,ioffy,ioffA);
  document.getElementById("debug_textarea").value +=
   "dsyr2(L,n,2.,x,incx,y,incy,A,lda): A = " + A[ioffA + 0] + " "
    + A[ioffA + 1] + " " + A[ioffA + 2] + " " + A[ioffA + 4] + " "
    + A[ioffA + 5] + " " + A[ioffA + 6] + " " + A[ioffA + 8] + " "
    + A[ioffA + 9] + " " + A[ioffA + 10] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zher2() {
  document.getElementById("debug_textarea").value +=
   "testing zher2 *************" + "\n";
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
      A[ ioffA + i + j * lda ] = new Complex( i * j , i - j );
    }
  }
  var incx = 2;
  var incy = 2;
  var zero = new Complex(0.,0.);
  var alpha = new Complex(2.,3.);
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2 , j + 2 );
  }
  for ( j = 0; j < 6; j ++ ) {
    y[ioffy + j] = new Complex( j + 1 , 1 - j );
  }
  Blas2.zher2("U",n,zero,x,incx,y,incy,A,lda,ioffx,ioffy,ioffA);
  document.getElementById("debug_textarea").value +=
   "zher2(U,n,zero,x,incx,y,incy,A,lda): A = " + A[ioffA + 0]
    + " " + A[ioffA + 1] + " " + A[ioffA + 2] + " " + A[ioffA + 4]
    + " " + A[ioffA + 5] + " " + A[ioffA + 6] + " " + A[ioffA + 8]
    + " " + A[ioffA + 9] + " " + A[ioffA + 10] + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i * j , i - j );
    }
  }
  Blas2.zher2("U",n,alpha,x,incx,y,incy,A,lda,ioffx,ioffy,ioffA);
  document.getElementById("debug_textarea").value +=
   "zher2(U,n,alpha,x,incx,y,incy,A,lda): A = " + A[ioffA + 0]
    + " " + A[ioffA + 1] + " " + A[ioffA + 2] + " " + A[ioffA + 4]
    + " " + A[ioffA + 5] + " " + A[ioffA + 6] + " " + A[ioffA + 8]
    + " " + A[ioffA + 9] + " " + A[ioffA + 10] + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i * j , i - j );
    }
  }
  Blas2.zher2("L",n,alpha,x,incx,y,incy,A,lda,ioffx,ioffy,ioffA);
  document.getElementById("debug_textarea").value +=
   "zher2(L,n,alpha,x,incx,y,incy,A,lda): A = " + A[ioffA + 0]
    + " " + A[ioffA + 1] + " " + A[ioffA + 2] + " " + A[ioffA + 4]
    + " " + A[ioffA + 5] + " " + A[ioffA + 6] + " " + A[ioffA + 8]
    + " " + A[ioffA + 9] + " " + A[ioffA + 10] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dspr2() {
  document.getElementById("debug_textarea").value +=
   "testing dspr2 *************" + "\n";
  var ioffA = 1;
  var ioffP = 2;
  var ioffx = 3;
  var ioffy = 4;
  var A = new Array( ioffA + 12 );
  var P = new Array( ioffP + 6 );
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 6 );
  var n = 3;
  var lda = 4;
  for (var j = 0; j < n; j ++ ) {
    for (var i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = Number( i + j );
    }
  }
  var k=0;
  for (j=0;j<n;j++) {
    for (i=0;i<=j;i++) {
      P[ioffP + k++]=A[ioffA + i+j*lda];
    }
  }
  var incx = 2;
  var incy = 2;
  for ( j = 0; j < 6; j ++ ) x[ioffx + j] = Number( j * 2 );
  for ( j = 0; j < 6; j ++ ) y[ioffy + j] = Number( j + 1 );
  Blas2.dspr2("U",n,0.,x,incx,y,incy,P,ioffx,ioffy,ioffP);
  document.getElementById("debug_textarea").value +=
   "dspr2(U,n,0.,x,incx,y,incy,P): P = " + P[ioffP + 0] + " "
    + P[ioffP + 1] + " " + P[ioffP + 2] + " " + P[ioffP + 3] + " "
    + P[ioffP + 4] + " " + P[ioffP + 5] + "\n";

  k=0;
  for (j=0;j<n;j++) {
    for (i=0;i<=j;i++) {
      P[ioffP + k++]=A[ioffA + i+j*lda];
    }
  }
  Blas2.dspr2("U",n,2.,x,incx,y,incy,P,ioffx,ioffy,ioffP);
  document.getElementById("debug_textarea").value +=
   "dspr2(U,n,2.,x,incx,y,incy,P): P = " + P[ioffP + 0] + " "
    + P[ioffP + 1] + " " + P[ioffP + 2] + " " + P[ioffP + 3] + " "
    + P[ioffP + 4] + " " + P[ioffP + 5] + "\n";

  k=0;
  for (j=0;j<n;j++) {
    for (i=j;i<n;i++) {
      P[ioffP + k++]=A[ioffA + i+j*lda];
    }
  }
  Blas2.dspr2("L",n,2.,x,incx,y,incy,P,ioffx,ioffy,ioffP);
  document.getElementById("debug_textarea").value +=
   "dspr2(L,n,2.,x,incx,y,incy,P): P = " + P[ioffP + 0] + " "
    + P[ioffP + 1] + " " + P[ioffP + 2] + " " + P[ioffP + 3] + " "
    + P[ioffP + 4] + " " + P[ioffP + 5] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zhpr2() {
  document.getElementById("debug_textarea").value +=
    "testing zhpr2 *************" + "\n";
  var ioffA = 1;
  var ioffP = 2;
  var ioffx = 3;
  var ioffy = 4;
  var A = new Array( ioffA + 12 );
  var P = new Array( ioffP + 6 );
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 6 );
  var n = 3;
  var lda = 4;
  for (var j = 0; j < n; j ++ ) {
    for (var i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j , i * j );
    }
  }
  var k=0;
  for (j=0;j<n;j++) {
    for (i=0;i<=j;i++) {
      P[ioffP + k] = new Complex();
      P[ioffP + k++].value=A[ioffA + i+j*lda];
    }
  }
  var incx = 2;
  var incy = 2;
  for ( j = 0; j < 6; j ++ ) {
    x[ioffx + j] = new Complex( j * 2 , j + 2 );
  }
  for ( j = 0; j < 6; j ++ ) {
    y[ioffy + j] = new Complex( i + 1 , 1 - i );
  }
  var zero = new Complex(0.,0.);
  var alpha = new Complex(2.,3.);
  Blas2.zhpr2("U",n,zero,x,incx,y,incy,P,ioffx,ioffy,ioffP);
  document.getElementById("debug_textarea").value +=
   "zhpr2(U,n,zero,x,incx,y,incy,P): P = " + P[ioffP + 0] + " "
    + P[ioffP + 1] + " " + P[ioffP + 2] + " " + P[ioffP + 3] + " "
    + P[ioffP + 4] + " " + P[ioffP + 5] + "\n";

  k=0;
  for (j=0;j<n;j++) {
    for (i=0;i<=j;i++) {
      P[ioffP + k] = new Complex();
      P[ioffP + k++].value=A[ioffA + i+j*lda];
    }
  }
  Blas2.zhpr2("U",n,alpha,x,incx,y,incy,P,ioffx,ioffy,ioffP);
  document.getElementById("debug_textarea").value +=
   "zhpr2(U,n,alpha,x,incx,y,incy,P): P = " + P[ioffP + 0] + " "
    + P[ioffP + 1] + " " + P[ioffP + 2] + " " + P[ioffP + 3] + " "
    + P[ioffP + 4] + " " + P[ioffP + 5] + "\n";

  k=0;
  for (j=0;j<n;j++) {
    for (i=j;i<n;i++) {
      P[ioffP + k] = new Complex();
      P[ioffP + k++].value=A[ioffA + i+j*lda];
    }
  }
  Blas2.zhpr2("L",n,alpha,x,incx,y,incy,P,ioffx,ioffy,ioffP);
  document.getElementById("debug_textarea").value +=
   "zhpr2(L,n,alpha,x,incx,y,incy,P): P = " + P[ioffP + 0] + " "
    + P[ioffP + 1] + " " + P[ioffP + 2] + " " + P[ioffP + 3] + " "
    + P[ioffP + 4] + " " + P[ioffP + 5] + "\n";
}
