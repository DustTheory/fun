function Blas3() {
}
//*************************************************************************
//  C <-- op(A) alpha op(B) + C beta, op(X) = X or X^T or X^H
Blas3.dgemm = function( transa , transb , m , n , k , alpha , A , lda ,
B , ldb , beta , C , ldc, ioffA, ioffB, ioffC ) {
  var nota = ( transa.charAt(0).toUpperCase() == 'N' );
  var notb = ( transb.charAt(0).toUpperCase() == 'N' );
  var nrowa = k;
  var ncola = m;
  if ( nota ) {
    nrowa = m;
    ncola = k;
  }
  var nrowb = n;
  if ( notb ) nrowb = k;
  var info = 0;
  if ( ! nota && ! transa.charAt(0).toUpperCase() != 'C' &&
  transa.charAt(0).toUpperCase() != 'T' ) {
    info = 1;
  } else if ( ! notb && transb.charAt(0).toUpperCase() != 'C' &&
  transb.charAt(0).toUpperCase() != 'T' ) {
    info = 2;
  } else if ( m < 0 ) info = 3;
  else if ( n < 0 ) info = 4;
  else if ( k < 0 ) info = 5;
  else if ( lda < Math.max( 1 , nrowa ) ) info = 8;
  else if ( ldb < Math.max( 1 , nrowb ) ) info = 10;
  else if ( ldc < Math.max( 1 , m ) ) info = 13;
  if ( info != 0 ) {
    Blas2.xerbla( 'dgemm' , info );
    return;
  }
  
  if ( m == 0 || n == 0 || ( ( alpha == 0. || k == 0 ) && beta == 1. ) ) {
    return;
  }
  var i = -1;
  var j = -1;
  var l = -1;
  var temp = Number.POSITIVE_INFINITY;
  if ( alpha == 0. ) {
    if ( beta == 0. ) {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = 1; i <= m; i ++ ) {
          C[ ioffC + i - 1 + ( j - 1 ) * ldc ] = 0.;
        }
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = 1; i <= m; i ++ ) {
          C[ ioffC + i - 1 + ( j - 1 ) * ldc ] *= beta;
        }
      }
    }
    return;
  }
  if ( notb ) {
    if ( nota ) {
      for ( j = 1; j <= n; j ++ ) {
        if ( beta == 0. ) {
          for ( i = 1; i <= m; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ] = 0.;
          }
        } else if ( beta != 1. ) {
          for ( i = 1; i <= m; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ] *= beta;
          }
        }
        for ( l = 1; l <= k; l ++ ) {
          if ( B[ ioffB + l - 1 + ( j - 1 ) * ldb ] != 0. ) {
            temp = alpha * B[ ioffB + l - 1 + ( j - 1 ) * ldb ];
            for ( i = 1; i <= m; i ++ ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ] +=
                temp * A[ ioffA + i - 1 + ( l - 1 ) * lda ];
            }
          }
        }
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = 1; i <= m; i ++ ) {
          temp = 0.;
          for ( l = 1; l <= k; l ++ ) {
            temp += A[ ioffA + l - 1 + ( i - 1 ) * lda ]
                  * B[ ioffB + l - 1 + ( j - 1 ) * ldb ];
          }
          if ( beta == 0. ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ] = alpha * temp;
          } else {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ] = alpha * temp
              + beta * C[ ioffC + i - 1 + ( j - 1 ) * ldc ];
          }
        }
      }
    }
  } else {
    if ( nota ) {
      for ( j = 1; j <= n; j ++ ) {
        if ( beta == 0. ) {
          for ( i = 1; i <= m; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ] = 0.;
          }
        } else if ( beta != 1. ) {
          for ( i = 1; i <= m; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ] *= beta;
          }
        }
        for ( l = 1; l <= k; l ++ ) {
          if ( B[ ioffB + j - 1 + ( l - 1 ) * ldb ] != 0. ) {
            temp = alpha * B[ ioffB + j - 1 + ( l - 1 ) * ldb ];
            for ( i = 1; i <= m; i ++ ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ] +=
                temp * A[ ioffA + i - 1 + ( l - 1 ) * lda ];
            }
          }
        }
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = 1; i <= m; i ++ ) {
          temp = 0.;
          for ( l = 1; l <= k; l ++ ) {
            temp += A[ ioffA + l - 1 + ( i - 1 ) * lda ]
                  * B[ ioffB + j - 1 + ( l - 1 ) * ldb ];
          }
          if ( beta == 0. ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ] = alpha * temp;
          } else {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ] = alpha * temp
              + beta * C[ ioffC + i - 1 + ( j - 1 ) * ldc ];
          }
        }
      }
    }
  }
}
Blas3.zgemm = function( transa , transb , m , n , k , alpha , A , lda ,
B , ldb , beta , C , ldc, ioffA, ioffB, ioffC ) {
  var nota = ( transa.charAt(0).toUpperCase() == 'N' );
  var notb = ( transb.charAt(0).toUpperCase() == 'N' );
  var conja = ( transa.charAt(0).toUpperCase() == 'C' );
  var conjb = ( transb.charAt(0).toUpperCase() == 'C' );
  var nrowa = k;
  var ncola = m;
  if ( nota ) {
    nrowa = m;
    ncola = k;
  }
  var nrowb = n;
  if ( notb ) nrowb = k;
  var info = 0;
  if ( ! nota && ! conja && transa.charAt(0).toUpperCase() != 'T' ) {
    info = 1;
  } else if ( ! notb && ! conjb &&
  transb.charAt(0).toUpperCase() != 'T' ) {
    info = 2;
  } else if ( m < 0 ) info = 3;
  else if ( n < 0 ) info = 4;
  else if ( k < 0 ) info = 5;
  else if ( lda < Math.max( 1 , nrowa ) ) info = 8;
  else if ( ldb < Math.max( 1 , nrowb ) ) info = 10;
  else if ( ldc < Math.max( 1 , m ) ) info = 13;
  if ( info != 0 ) {
    Blas2.xerbla( 'zgemm' , info );
    return;
  }
  
  var zero = new Complex( 0., 0. );
  var one = new Complex( 1., 0. );
  if ( m == 0 || n == 0 ||
  ( ( alpha.equals( zero ) || k == 0 ) && beta.equals( one ) ) ) {
    return;
  }
  var i = -1;
  var j = -1;
  var l = -1;
  var temp = new Complex();
  if ( alpha.equals( zero ) ) {
    if ( beta.equals( zero ) ) {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = 1; i <= m; i ++ ) {
          C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
        }
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = 1; i <= m; i ++ ) {
          C[ ioffC + i - 1 + ( j - 1 ) * ldc ].timesEquals( beta );
        }
      }
    }
    return;
  }
  if ( notb ) {
    if ( nota ) {
      for ( j = 1; j <= n; j ++ ) {
        if ( beta.equals( zero ) ) {
          for ( i = 1; i <= m; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
          }
        } else if ( ! beta.equals( one ) ) {
          for ( i = 1; i <= m; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].timesEquals( beta );
          }
        }
        for ( l = 1; l <= k; l ++ ) {
          if ( B[ ioffB + l - 1 + ( j - 1 ) * ldb ] != 0. ) {
            temp.setValue( ComplexMath.times( alpha ,
              B[ ioffB + l - 1 + ( j - 1 ) * ldb ] ) );
            for ( i = 1; i <= m; i ++ ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].plusEquals(
                ComplexMath.times( temp ,
                A[ ioffA + i - 1 + ( l - 1 ) * lda ] ) );
            }
          }
        }
      }
    } else if ( conja ) {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = 1; i <= m; i ++ ) {
          temp.setValue( zero );
          for ( l = 1; l <= k; l ++ ) {
            temp.plusEquals( ComplexMath.times( ComplexMath.conj(
              A[ ioffA + l - 1 + ( i - 1 ) * lda ] ) ,
              B[ ioffB + l - 1 + ( j - 1 ) * ldb ] ) );
          }
          if ( beta.equals( zero ) ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
              ComplexMath.times( alpha , temp ) );
          } else {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
              ComplexMath.plus( ComplexMath.times( alpha , temp ) ,
              ComplexMath.times( beta ,
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ] ) ) );
          }
        }
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = 1; i <= m; i ++ ) {
          temp.setValue( zero );
          for ( l = 1; l <= k; l ++ ) {
            temp.plusEquals( ComplexMath.times(
              A[ ioffA + l - 1 + ( i - 1 ) * lda ] ,
              B[ ioffB + l - 1 + ( j - 1 ) * ldb ] ) );
          }
          if ( beta.equals( zero ) ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
              ComplexMath.times( alpha , temp ) );
          } else {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
              ComplexMath.plus( ComplexMath.times( alpha , temp ) ,
              ComplexMath.times( beta ,
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ] ) ) );
          }
        }
      }
    }
  } else {
    if ( nota ) {
      if ( conjb ) {
        for ( j = 1; j <= n; j ++ ) {
          if ( beta.equals( zero ) ) {
            for ( i = 1; i <= m; i ++ ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
            }
          } else if ( ! beta.equals( one ) ) {
            for ( i = 1; i <= m; i ++ ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].timesEquals(
                beta );
            }
          }
          for ( l = 1; l <= k; l ++ ) {
            if ( B[ ioffB + j - 1 + ( l - 1 ) * ldb ] != 0. ) {
              temp.setValue( ComplexMath.times( alpha ,
                ComplexMath.conj(
                B[ ioffB + j - 1 + ( l - 1 ) * ldb ] ) ) );
              for ( i = 1; i <= m; i ++ ) {
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ].plusEquals(
                  ComplexMath.times( temp ,
                    A[ ioffA + i - 1 + ( l - 1 ) * lda ] ) );
              }
            }
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          if ( beta.equals( zero ) ) {
            for ( i = 1; i <= m; i ++ ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
            }
          } else if ( ! beta.equals( one ) ) {
            for ( i = 1; i <= m; i ++ ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].timesEquals(
                beta );
            }
          }
          for ( l = 1; l <= k; l ++ ) {
            if ( B[ ioffB + j - 1 + ( l - 1 ) * ldb ] != 0. ) {
              temp.setValue( ComplexMath.times(
                alpha , B[ ioffB + j - 1 + ( l - 1 ) * ldb ] ) );
              for ( i = 1; i <= m; i ++ ) {
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ].plusEquals(
                  ComplexMath.times( temp ,
                    A[ ioffA + i - 1 + ( l - 1 ) * lda ] ) );
              }
            }
          }
        }
      }
    } else {
      if ( conja) {
        if ( conjb ) {
          for ( j = 1; j <= n; j ++ ) {
            for ( i = 1; i <= m; i ++ ) {
              temp.setValue( zero );
              for ( l = 1; l <= k; l ++ ) {
                temp.plusEquals( ComplexMath.times(
                  ComplexMath.conj(
                  A[ ioffA + l - 1 + ( i - 1 ) * lda ] ) ,
                  ComplexMath.conj(
                  B[ ioffB + j - 1 + ( l - 1 ) * ldb ] ) ) );
              }
              if ( beta.equals( zero ) ) {
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                  ComplexMath.times( alpha , temp ) );
              } else {
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                  ComplexMath.plus(
                  ComplexMath.times( alpha , temp ) ,
                  ComplexMath.times( beta ,
                    C[ ioffC + i - 1 + ( j - 1 ) * ldc ] ) ) );
              }
            }
          }
        } else {
          for ( j = 1; j <= n; j ++ ) {
            for ( i = 1; i <= m; i ++ ) {
              temp.setValue( zero );
              for ( l = 1; l <= k; l ++ ) {
                temp.plusEquals( ComplexMath.times( ComplexMath.conj(
                  A[ ioffA + l - 1 + ( i - 1 ) * lda ] ) ,
                  B[ ioffB + j - 1 + ( l - 1 ) * ldb ] ) );
              }
              if ( beta.equals( zero ) ) {
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                  ComplexMath.times( alpha , temp ) );
              } else {
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                  ComplexMath.plus(
                  ComplexMath.times( alpha , temp ) ,
                  ComplexMath.times( beta ,
                    C[ ioffC + i - 1 + ( j - 1 ) * ldc ] ) ) );
              }
            }
          }
        }
      } else {
        if ( conjb ) {
          for ( j = 1; j <= n; j ++ ) {
            for ( i = 1; i <= m; i ++ ) {
              temp.setValue( zero );
              for ( l = 1; l <= k; l ++ ) {
                temp.plusEquals( ComplexMath.times(
                  A[ ioffA + l - 1 + ( i - 1 ) * lda ] ,
                  ComplexMath.conj(
                  B[ ioffB + j - 1 + ( l - 1 ) * ldb ] ) ) );
              }
              if ( beta.equals( zero ) ) {
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                  ComplexMath.times( alpha , temp ) );
              } else {
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                  ComplexMath.plus(
                  ComplexMath.times( alpha , temp ) ,
                  ComplexMath.times( beta ,
                    C[ ioffC + i - 1 + ( j - 1 ) * ldc ] ) ) );
              }
            }
          }
        } else {
          for ( j = 1; j <= n; j ++ ) {
            for ( i = 1; i <= m; i ++ ) {
              temp.setValue( zero );
              for ( l = 1; l <= k; l ++ ) {
                temp.plusEquals( ComplexMath.times(
                  A[ ioffA + l - 1 + ( i - 1 ) * lda ] ,
                    B[ ioffB + j - 1 + ( l - 1 ) * ldb ] ) );
              }
              if ( beta.equals( zero ) ) {
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                  ComplexMath.times( alpha , temp ) );
              } else {
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                  ComplexMath.plus(
                  ComplexMath.times( alpha , temp ) ,
                  ComplexMath.times( beta ,
                    C[ ioffC + i - 1 + ( j - 1 ) * ldc ] ) ) );
              }
            }
          }
        }
      }
    }
  }
}
//*************************************************************************
//  C <-- A alpha B + C beta or C <-- B alpha A + C beta, A symmetric
Blas3.dsymm = function( side , uplo , m , n , alpha , A , lda , B , ldb ,
beta , C , ldc, ioffA, ioffB, ioffC ) {
  var nrowa = n;
  if ( side.charAt(0).toUpperCase() == 'L' ) nrowa = m;
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var info = 0;
  if ( side.charAt(0).toUpperCase() != 'L'
  && side.charAt(0).toUpperCase() != 'R' ) {
    info = 1;
  } else if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' )  {
    info = 2;
  } else if ( m < 0 ) info = 3;
  else if ( n < 0 ) info = 4;
  else if ( lda < Math.max( 1 , nrowa ) ) info = 7;
  else if ( ldb < Math.max( 1 , m ) ) info = 9;
  else if ( ldc < Math.max( 1 , m ) ) info = 12;
  if ( info != 0 ) {
    Blas2.xerbla( 'dsymm' , info );
    return;
  }

  if ( m == 0 || n == 0 || ( alpha == 0. && beta == 1. ) ) return;
  var i = -1;
  var j = -1;
  if ( alpha == 0. ) {
    if ( beta == 0. ) {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = 1; i <= m; i ++ ) {
          C[ ioffC + i - 1 + ( j - 1 ) * ldc ] = 0.;
        }
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = 1; i <= m; i ++ ) {
          C[ ioffC + i - 1 + ( j - 1 ) * ldc ] *= beta;
        }
      }
    }
  } else {
    var k = -1;
    var temp1 = Number.POSITIVE_INFINITY;
    var temp2 = Number.POSITIVE_INFINITY;
    if ( side.charAt(0).toUpperCase() == 'L' ) {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= m; i ++ ) {
            temp1 = alpha * B[ ioffB + i - 1 + ( j - 1 ) * ldb ];
            temp2 = 0.;
            for ( k = 1; k <= i - 1; k ++ ) {
              C[ ioffC + k - 1 + ( j - 1 ) * ldc ] +=
                temp1 * A[ ioffA + k - 1 + ( i - 1 ) * lda ];
              temp2 += B[ ioffB + k - 1 + ( j - 1 ) * ldb ]
                     * A[ ioffA + k - 1 + ( i - 1 ) * lda ];
            }
            if ( beta == 0. ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ] =
                temp1 * A[ ioffA + i - 1 + ( i - 1 ) * lda ]
                + alpha * temp2;
            } else {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ] =
                beta * C[ ioffC + i - 1 + ( j - 1 ) * ldc ]
                + temp1 * A[ ioffA + i - 1 + ( i - 1 ) * lda ]
                + alpha * temp2;
            }
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = m; i >= 1; i -- ) {
            temp1 = alpha * B[ ioffB + i - 1 + ( j - 1 ) * ldb ];
            temp2 = 0.;
            for ( k = i + 1; k <= m; k ++ ) {
              C[ ioffC + k - 1 + ( j - 1 ) * ldc ] +=
                temp1 * A[ ioffA + k - 1 + ( i - 1 ) * lda ];
              temp2 += B[ ioffB + k - 1 + ( j - 1 ) * ldb ]
                     * A[ ioffA + k - 1 + ( i - 1 ) * lda ];
            }
            if ( beta == 0. ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ] =
                temp1 * A[ ioffA + i - 1 + ( i - 1 ) * lda ]
                + alpha * temp2;
            } else {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ] =
                beta * C[ ioffC + i - 1 + ( j - 1 ) * ldc ]
                + temp1 * A[ ioffA + i - 1 + ( i - 1 ) * lda ]
                + alpha * temp2;
            }
          }
        }
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        temp1 = alpha * A[ ioffA + j - 1 + ( j - 1 ) * lda ];
        if ( beta == 0. ) {
          for ( i = 1; i <= m; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ] =
              temp1 * B[ ioffB + i - 1 + ( j - 1 ) * ldb ];
          }
        } else {
          for ( i = 1; i <= m; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ] =
              beta * C[ ioffC + i - 1 + ( j - 1 ) * ldc ]
              + temp1 * B[ ioffB + i - 1 + ( j - 1 ) * ldb ];
          }
        }
        for ( k = 1; k < j; k ++ ) {
          if ( upper ) {
            temp1 = alpha * A[ ioffA + k - 1 + ( j - 1 ) * lda ];
          } else {
            temp1 = alpha * A[ ioffA + j - 1 + ( k - 1 ) * lda ];
          }
          for ( i = 1; i <= m; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ] +=
              temp1 * B[ ioffB + i - 1 + ( k - 1 ) * ldb ];
          }
        }
        for ( k = j + 1; k <= n; k ++ ) {
          if ( upper ) {
            temp1 = alpha * A[ ioffA + j - 1 + ( k - 1 ) * lda ];
          } else {
            temp1 = alpha * A[ ioffA + k - 1 + ( j - 1 ) * lda ];
          }
          for ( i = 1; i <= m; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ] +=
              temp1 * B[ ioffB + i - 1 + ( k - 1 ) * ldb ];
          }
        }
      }
    }
  }
}
Blas3.zsymm = function( side , uplo , m , n , alpha , A , lda , B , ldb ,
beta , C , ldc, ioffA, ioffB, ioffC ) {
  var nrowa = n;
  if ( side.charAt(0).toUpperCase() == 'L' ) nrowa = m;
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var info = 0;
  if ( side.charAt(0).toUpperCase() != 'L' 
  && side.charAt(0).toUpperCase() != 'R' ) {
    info = 1;
  } else if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 2;
  } else if ( m < 0 ) info = 3;
  else if ( n < 0 ) info = 4;
  else if ( lda < Math.max( 1 , nrowa ) ) info = 7;
  else if ( ldb < Math.max( 1 , m ) ) info = 9;
  else if ( ldc < Math.max( 1 , m ) ) info = 12;
  if ( info != 0 ) {
    Blas2.xerbla( 'zsymm' , info );
    return;
  }

  var zero = new Complex( 0., 0. );
  var one = new Complex( 1., 0. );
  if ( m == 0 || n == 0 ||
  ( alpha.equals(zero) && beta.equals(one) )) {
    return;
  }
  var i = -1;
  var j = -1;
  if ( alpha.equals(zero) ) {
    if ( beta.equals(zero) ) {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = 1; i <= m; i ++ ) {
          C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
        }
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = 1; i <= m; i ++ ) {
          C[ ioffC + i - 1 + ( j - 1 ) * ldc ].timesEquals(beta);
        }
      }
    }
  } else {
    var k = -1;
    var temp1 = new Complex();
    var temp2 = new Complex();
    if ( side.charAt(0).toUpperCase() == 'L' ) {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= m; i ++ ) {
            temp1.setValue( ComplexMath.times( alpha ,
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ] ) );
            temp2.setValue( zero );
            for ( k = 1; k <= i - 1; k ++ ) {
              C[ ioffC + k - 1 + ( j - 1 ) * ldc ].plusEquals(
                ComplexMath.times( temp1 ,
                A[ ioffA + k - 1 + ( i - 1 ) * lda ] ) );
              temp2.plusEquals( ComplexMath.times(
                B[ ioffB + k - 1 + ( j - 1 ) * ldb ] ,
                A[ ioffA + k - 1 + ( i - 1 ) * lda ] ));
            }
            if ( beta.equals( zero ) ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                ComplexMath.plus( ComplexMath.times( temp1 ,
                A[ ioffA + i - 1 + ( i - 1 ) * lda ] ) ,
                ComplexMath.times( alpha , temp2 ) ) );
            } else {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                ComplexMath.plus( ComplexMath.times( beta , 
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ] ) ,
                ComplexMath.plus( ComplexMath.times( temp1 ,
                A[ ioffA + i - 1 + ( i - 1 ) * lda ] ) ,
                ComplexMath.times( alpha , temp2 ) ) ) );
            }
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = m; i >= 1; i -- ) {
            temp1.setValue( ComplexMath.times( alpha ,
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ] ) );
            temp2.setValue( zero );
            for ( k = i + 1; k <= m; k ++ ) {
              C[ ioffC + k - 1 + ( j - 1 ) * ldc ].plusEquals(
                ComplexMath.times( temp1 ,
                A[ ioffA + k - 1 + ( i - 1 ) * lda ] ) );
              temp2.plusEquals( ComplexMath.times(
                B[ ioffB + k - 1 + ( j - 1 ) * ldb ] ,
                A[ ioffA + k - 1 + ( i - 1 ) * lda ] ) );
            }
            if ( beta.equals( zero ) ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                ComplexMath.plus( ComplexMath.times( temp1 ,
                A[ ioffA + i - 1 + ( i - 1 ) * lda ] ) ,
                ComplexMath.times( alpha , temp2 ) ) );
            } else {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                ComplexMath.plus( ComplexMath.times( beta ,
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ] ) ,
                ComplexMath.plus( ComplexMath.times( temp1 ,
                A[ ioffA + i - 1 + ( i - 1 ) * lda ] ) ,
                ComplexMath.times( alpha , temp2 ) ) ) );
            }
          }
        }
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        temp1.setValue( ComplexMath.times( alpha ,
          A[ ioffA + j - 1 + ( j - 1 ) * lda ] ) );
        if ( beta.equals( zero ) ) {
          for ( i = 1; i <= m; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
              ComplexMath.times( temp1 ,
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ] ) );
          }
        } else {
          for ( i = 1; i <= m; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
              ComplexMath.plus( ComplexMath.times( beta ,
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ] ) ,
              ComplexMath.times( temp1 ,
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ] ) ) );
          }
        }
        for ( k = 1; k <= j - 1; k ++ ) {
          if ( upper ) {
            temp1.setValue( ComplexMath.times( alpha ,
              A[ ioffA + k - 1 + ( j - 1 ) * lda ] ) );
          } else {
            temp1.setValue( ComplexMath.times( alpha ,
              A[ ioffA + j - 1 + ( k - 1 ) * lda ] ) );
          }
          for ( i = 1; i <= m; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].plusEquals(
              ComplexMath.times( temp1 ,
              B[ ioffB + i - 1 + ( k - 1 ) * ldb ] ) );
          }
        }
        for ( k = j + 1; k <= n; k ++ ) {
          if ( upper ) {
            temp1.setValue( ComplexMath.times( alpha ,
              A[ ioffA + j - 1 + ( k - 1 ) * lda ] ) );
          } else {
            temp1.setValue( ComplexMath.times( alpha ,
              A[ ioffA + k - 1 + ( j - 1 ) * lda ] ) );
          }
          for ( i = 1; i <= m; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].plusEquals(
              ComplexMath.times( temp1 ,
              B[ ioffB + i - 1 + ( k - 1 ) * ldb ] ) );
          }
        }
      }
    }
  }
}
Blas3.zhemm = function( side , uplo , m , n , alpha , A , lda , B , ldb ,
beta , C , ldc, ioffA, ioffB, ioffC ) {
  var nrowa = n;
  if ( side.charAt(0).toUpperCase() == 'L' ) nrowa = m;
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var info = 0;
  if ( side.charAt(0).toUpperCase() != 'L' 
  && side.charAt(0).toUpperCase() != 'R' ) {
    info = 1;
  } else if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 2;
  } else if ( m < 0 ) info = 3;
  else if ( n < 0 ) info = 4;
  else if ( lda < Math.max( 1 , nrowa ) ) info = 7;
  else if ( ldb < Math.max( 1 , m ) ) info = 9;
  else if ( ldc < Math.max( 1 , m ) ) info = 12;
  if ( info != 0 ) {
    Blas2.xerbla( 'zhemm' , info );
    return;
  }

  var zero = new Complex( 0. , 0. );
  var one = new Complex( 1. , 0. );
  if ( m == 0 || n == 0 ||
  ( alpha.equals( zero ) && beta.equals( one ) ) ) {
    return;
  }
  var i = -1;
  var j = -1;
  if ( alpha.equals( zero ) ) {
    if ( beta.equals( zero ) ) {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = 1; i <= m; i ++ ) {
          C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
        }
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = 1; i <= m; i ++ ) {
          C[ ioffC + i - 1 + ( j - 1 ) * ldc ].timesEquals( beta );
        }
      }
    }
  } else {
    var k = -1;
    var temp1 = new Complex();
    var temp2 = new Complex();
    if ( side.charAt(0).toUpperCase() == 'L' ) {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= m; i ++ ) {
            temp1.setValue( ComplexMath.times( alpha ,
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ] ) );
            temp2.setValue( zero );
            for ( k = 1; k <= i - 1; k ++ ) {
              C[ ioffC + k - 1 + ( j - 1 ) * ldc ].plusEquals(
                ComplexMath.times( temp1 ,
                A[ ioffA + k - 1 + ( i - 1 ) * lda ] ) );
              temp2.plusEquals( ComplexMath.times(
                B[ ioffB + k - 1 + ( j - 1 ) * ldb ] ,
                ComplexMath.conj(
                A[ ioffA + k - 1 + ( i - 1 ) * lda ] ) ) );
            }
            if ( beta.equals( zero ) ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                ComplexMath.plus( ComplexMath.timesNumber( temp1 ,
                  A[ ioffA + i - 1 + ( i - 1 ) * lda ].getReal() ) ,
                ComplexMath.times( alpha , temp2 ) ) );
            } else {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                ComplexMath.plus( ComplexMath.times( beta ,
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ] ) ,
                ComplexMath.plus( ComplexMath.timesNumber( temp1 ,
                    A[ ioffA + i - 1 + ( i - 1 ) * lda ].getReal() ) ,
                  ComplexMath.times( alpha , temp2 ) ) ) );
            }
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = m; i >= 1; i -- ) {
            temp1.setValue( ComplexMath.times( alpha ,
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ] ) );
            temp2.setValue( zero );
            for ( k = i + 1; k <= m; k ++ ) {
              C[ ioffC + k - 1 + ( j - 1 ) * ldc ].plusEquals(
                ComplexMath.times( temp1 ,
                A[ ioffA + k - 1 + ( i - 1 ) * lda ] ) );
              temp2.plusEquals( ComplexMath.times(
                B[ ioffB + k - 1 + ( j - 1 ) * ldb ] ,
                ComplexMath.conj(
                A[ ioffA + k - 1 + ( i - 1 ) * lda ] ) ) );
            }
            if ( beta.equals( zero ) ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                ComplexMath.plus( ComplexMath.timesNumber( temp1 ,
                  A[ ioffA + i - 1 + ( i - 1 ) * lda ].getReal() ) ,
                ComplexMath.times( alpha , temp2 ) ) );
            } else {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                ComplexMath.plus( ComplexMath.times( beta ,
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ] ) ,
                ComplexMath.plus( ComplexMath.timesNumber( temp1 ,
                    A[ ioffA + i - 1 + ( i - 1 ) * lda ].getReal() ) ,
                  ComplexMath.times( alpha , temp2 ) ) ) );
            }
          }
        }
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        temp1.setValue( ComplexMath.timesNumber( alpha ,
          A[ ioffA + j - 1 + ( j - 1 ) * lda ].getReal() ) );
        if ( beta.equals( zero ) ) {
          for ( i = 1; i <= m; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
              ComplexMath.times( temp1 ,
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ] ) );
          }
        } else {
          for ( i = 1; i <= m; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
              ComplexMath.plus( ComplexMath.times( beta ,
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ] ) ,
              ComplexMath.times( temp1 ,
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ] ) ) );
          }
        }
        for ( k = 1; k <= j - 1; k ++ ) {
          if ( upper ) {
            temp1.setValue( ComplexMath.times( alpha ,
              A[ ioffA + k - 1 + ( j - 1 ) * lda ] ) );
          } else {
            temp1.setValue( ComplexMath.times( alpha , ComplexMath.conj(
              A[ ioffA + j - 1 + ( k - 1 ) * lda ] ) ) );
          }
          for ( i = 1; i <= m; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].plusEquals(
              ComplexMath.times( temp1 ,
              B[ ioffB + i - 1 + ( k - 1 ) * ldb ] ) );
          }
        }
        for ( k = j + 1; k <= n; k ++ ) {
          if ( upper ) {
            temp1.setValue( ComplexMath.times( alpha , ComplexMath.conj(
              A[ ioffA + j - 1 + ( k - 1 ) * lda ] ) ) );
          } else {
            temp1.setValue( ComplexMath.times( alpha ,
              A[ ioffA + k - 1 + ( j - 1 ) * lda ] ) );
          }
          for ( i = 1; i <= m; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].plusEquals( 
              ComplexMath.times( temp1 ,
              B[ ioffB + i - 1 + ( k - 1 ) * ldb ] ) );
          }
        }
      }
    }
  }
}
//*************************************************************************
//  C <-- A alpha A^T + C beta or C <-- A^T alpha A + C beta, C symmetric
Blas3.dsyrk = function( uplo , trans , n , k , alpha , A , lda , beta ,
C , ldc, ioffA, ioffC ) {
  var nrowa =  ( trans.charAt(0).toUpperCase() == 'N' ? n : k );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var info = 0;
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( trans.charAt(0).toUpperCase() != 'N'
  && trans.charAt(0).toUpperCase() != 'T' 
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info = 2;
  } else if ( n < 0 ) info = 3;
  else if ( k < 0 ) info = 4;
  else if ( lda < Math.max( 1 , nrowa ) ) info = 7;
  else if ( ldc < Math.max( 1 , n ) ) info = 10;
  if ( info != 0 ) {
    Blas2.xerbla( 'dsyrk' , info );
    return;
  }

  if ( n == 0 || ( ( alpha == 0. || k == 0 ) && beta == 1. ) ) return;
  var j = -1;
  var i = -1;
  if ( alpha == 0. ) {
    if ( upper ) {
      if ( beta == 0. ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= j; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ] = 0.;
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= j; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ] *= beta;
          }
        }
      }
    } else {
      if ( beta == 0. ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = j; i <= n; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ] = 0.;
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = j; i <= n; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ] *= beta;
          }
        }
      }
    }
  } else {
    var l = -1;
    var temp = Number.POSITIVE_INFINITY;
    if ( trans.charAt(0).toUpperCase() == 'N' ) {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          if ( beta == 0. ) {
            for ( i = 1; i <= j; i ++ ) {
              C[ioffC + i - 1 + ( j - 1 ) * ldc ] = 0.;
            }
          } else if ( beta != 1. ) {
            for ( i = 1; i <= j; i ++ ) {
              C[ioffC + i - 1 + ( j - 1 ) * ldc ] *= beta;
            }
          }
          for ( l = 1; l <= k; l ++ ) {
            if ( A[ ioffA + j - 1 + ( l - 1 ) * lda ] != 0. ) {
              temp = alpha * A[ ioffA + j - 1 + ( l - 1 ) * lda ];
              for ( i = 1; i <= j; i ++ ) {
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ] +=
                  temp * A[ ioffA + i - 1 + ( l - 1 ) * lda ];
              }
            }
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          if ( beta == 0. ) { 
            for ( i = j; i <= n; i ++ ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ] = 0.;
            }
          } else if ( beta != 1. ) {
            for ( i = j; i <= n; i ++ ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ] *= beta;
            }
          }
          for ( l = 1; l <= k; l ++ ) {
            if ( A[ ioffA + j - 1 + ( l - 1 ) * lda ] != 0. ) { 
              temp = alpha * A[ ioffA + j - 1 + ( l - 1 ) * lda ];
              for ( i = j; i <= n; i ++ ) {
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ] +=
                  temp * A[ ioffA + i - 1 + ( l - 1 ) * lda ];
              }
            }
          }
        }
      }
    } else {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= j; i ++ ) {
            temp = 0.;
            for ( l = 1; l <= k; l ++ ) {
              temp += A[ ioffA + l - 1 + ( i - 1 ) * lda ]
                    * A[ ioffA + l - 1 + ( j - 1 ) * lda ];
            }
            if ( beta == 0. ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ] = alpha * temp;
            } else {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ] = alpha * temp
                + beta * C[ ioffC + i - 1 + ( j - 1 ) * ldc ];
            }
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = j; i <= n; i ++ ) {
            temp = 0.;
            for ( l = 1; l <= k; l ++ ) {
              temp += A[ ioffA + l - 1 + ( i - 1 ) * lda ]
                    * A[ ioffA + l - 1 + ( j - 1 ) * lda ];
            }
            if ( beta == 0. ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ] = alpha * temp;
            } else {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ] = alpha * temp
                + beta * C[ ioffC + i - 1 + ( j - 1 ) * ldc ];
            }
          }
        }
      }
    }
  }
}
Blas3.zsyrk = function( uplo , trans , n , k , alpha , A , lda , beta ,
C , ldc, ioffA, ioffC ) {
  var nrowa = ( trans.charAt(0).toUpperCase() == 'N' ? n : k );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var info = 0;
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( trans.charAt(0).toUpperCase() != 'N' 
  && trans.charAt(0).toUpperCase() != 'T' ) {
    info = 2;
  } else if ( n < 0 ) info = 3;
  else if ( k < 0 ) info = 4;
  else if ( lda < Math.max( 1 , nrowa ) ) info = 7;
  else if ( ldc < Math.max( 1 , n ) ) info = 10;
  if ( info != 0 ) {
    Blas2.xerbla( 'zsyrk' , info );
    return;
  }

  var zero = new Complex( 0., 0. );
  var one = new Complex( 1., 0. );
  if ( n == 0 ||
  ( ( alpha.equals( zero ) || k == 0 ) && beta.equals( one ) ) ) {
    return;
  }
  var j = -1;
  var i = -1;
  if ( alpha.equals( zero ) ) {
    if ( upper ) {
      if ( beta.equals( zero ) ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= j; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= j; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].timesEquals( beta );
          }
        }
      }
    } else {
      if ( beta.equals( zero ) ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = j; i <= n; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = j; i <= n; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].timesEquals( beta );
          }
        }
      }
    }
  } else {
    var l = -1;
    var temp = new Complex();
    if ( trans.charAt(0).toUpperCase() == 'N' ) {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          if ( beta.equals( zero ) ) {
            for ( i = 1; i <= j; i ++ ) {
              C[ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
            }
          } else if ( ! beta.equals( one ) ) {
            for ( i = 1; i <= j; i ++ ) {
              C[ioffC + i - 1 + ( j - 1 ) * ldc ].timesEquals( beta );
            }
          }
          for ( l = 1; l <= k; l ++ ) {
            if ( ! A[ ioffA + j - 1 + ( l - 1 ) * lda ].equals( zero )
            ) {
              temp.setValue( ComplexMath.times( alpha ,
                A[ ioffA + j - 1 + ( l - 1 ) * lda ] ) );
              for ( i = 1; i <= j; i ++ ) {
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ].plusEquals(
                  ComplexMath.times( temp ,
                  A[ ioffA + i - 1 + ( l - 1 ) * lda ] ) );
              }
            }
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          if ( beta.equals( zero ) ) { 
            for ( i = j; i <= n; i ++ ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
            }
          } else if ( ! beta.equals( one ) ) {
            for ( i = j; i <= n; i ++ ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].timesEquals( beta );
            }
          }
          for ( l = 1; l <= k; l ++ ) {
            if ( A[ ioffA + j - 1 + ( l - 1 ) * lda ] != 0. ) { 
              temp.setValue( ComplexMath.times( alpha ,
                A[ ioffA + j - 1 + ( l - 1 ) * lda ] ) );
              for ( i = j; i <= n; i ++ ) {
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ].plusEquals(
                  ComplexMath.times( temp ,
                    A[ ioffA + i - 1 + ( l - 1 ) * lda ] ) );
              }
            }
          }
        }
      }
    } else {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= j; i ++ ) {
            temp.setValue( zero );
            for ( l = 1; l <= k; l ++ ) {
              temp.plusEquals( ComplexMath.times(
                A[ ioffA + l - 1 + ( i - 1 ) * lda ] ,
                A[ ioffA + l - 1 + ( j - 1 ) * lda ] ) );
            }
            if ( beta.equals( zero ) ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                ComplexMath.times( alpha , temp ) );
            } else {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                ComplexMath.plus( ComplexMath.times( alpha , temp ) ,
                ComplexMath.times( beta ,
                  C[ ioffC + i - 1 + ( j - 1 ) * ldc ] ) ) );
            }
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = j; i <= n; i ++ ) {
            temp.setValue( zero );
            for ( l = 1; l <= k; l ++ ) {
              temp.plusEquals( ComplexMath.times(
                A[ ioffA + l - 1 + ( i - 1 ) * lda ] ,
                A[ ioffA + l - 1 + ( j - 1 ) * lda ] ) );
            }
            if ( beta.equals( zero ) ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                ComplexMath.times( alpha , temp ) );
            } else {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                ComplexMath.plus( ComplexMath.times( alpha , temp ) ,
                ComplexMath.times( beta ,
                  C[ ioffC + i - 1 + ( j - 1 ) * ldc ] ) ) );
            }
          }
        }
      }
    }
  }
}
Blas3.zherk = function( uplo , trans , n , k , alpha , A , lda , beta ,
C , ldc, ioffA, ioffC ) {
  var nrowa =  ( trans.charAt(0).toUpperCase() == 'N' ? n : k );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var info = 0;
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( trans.charAt(0).toUpperCase() != 'N' 
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info = 2;
  } else if ( n < 0 ) info = 3;
  else if ( k < 0 ) info = 4;
  else if ( lda < Math.max( 1 , nrowa ) ) info = 7;
  else if ( ldc < Math.max( 1 , n ) ) info = 10;
  if ( info != 0 ) {
    Blas2.xerbla( 'zherk' , info );
    return;
  }

  if ( n == 0 || ( ( alpha == 0. || k == 0 ) && beta == 1. ) ) {
    return;
  }
  var j = -1;
  var i = -1;
  var zero = new Complex( 0., 0. );
  if ( alpha == 0. ) {
    if ( upper ) {
      if ( beta == 0. ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= j; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= j; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].timesNumberEquals(
              beta );
          }
          C[ ioffC + j - 1 + ( j - 1 ) * ldc ].pureReal();
        }
      }
    } else {
      if ( beta == 0. ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = j; i <= n; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = j; i <= n; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].timesNumberEquals(
              beta );
          }
          C[ ioffC + j - 1 + ( j - 1 ) * ldc ].pureReal();
        }
      }
    }
    return;
  }
  var l = -1;
  var temp = new Complex();
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    if ( upper ) {
      for ( j = 1; j <= n; j ++ ) {
        if ( beta == 0. ) {
          for ( i = 1; i <= j; i ++ ) {
            C[ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
          }
        } else if ( beta != 1. ) {
          for ( i = 1; i <= j; i ++ ) {
            C[ioffC + i - 1 + ( j - 1 ) * ldc ].timesNumberEquals(
              beta );
          }
          C[ ioffC + j - 1 + ( j - 1 ) * ldc ].pureReal();
        }
        for ( l = 1; l <= k; l ++ ) {
          if ( ! A[ ioffA + j - 1 + ( l - 1 ) * lda ].equals( zero )
          ) {
            temp.setValue( ComplexMath.timesNumber( ComplexMath.conj(
              A[ ioffA + j - 1 + ( l - 1 ) * lda ] ) , alpha ) );
            for ( i = 1; i <= j; i ++ ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].plusEquals(
                ComplexMath.times( temp ,
                A[ ioffA + i - 1 + ( l - 1 ) * lda ] ) );
            }
            C[ ioffC + j - 1 + ( j - 1 ) * ldc ].pureReal();
          }
        }
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        if ( beta == 0. ) { 
          for ( i = j; i <= n; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
          }
        } else if ( beta != 1. ) {
          C[ ioffC + j - 1 + ( j - 1 ) * ldc ].timesNumberEquals(
            beta );
          C[ ioffC + j - 1 + ( j - 1 ) * ldc ].pureReal();
          for ( i = j + 1; i <= n; i ++ ) {
            C[ioffC + i - 1 + ( j - 1 ) * ldc ].timesNumberEquals(
              beta );
          }
        }
        for ( l = 1; l <= k; l ++ ) {
          if ( A[ ioffA + j - 1  + ( l - 1 ) * lda ] != 0. ) { 
            temp.setValue( ComplexMath.timesNumber( ComplexMath.conj(
              A[ ioffA + j - 1 + ( l - 1 ) * lda ] ) , alpha ) );
            for ( i = j; i <= n; i ++ ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].plusEquals(
                ComplexMath.times( temp ,
                  A[ ioffA + i - 1 + ( l - 1 ) * lda ] ) );
            }
            C[ ioffC + j - 1 + ( j - 1 ) * ldc ].pureReal();
          }
        }
      }
    }
  } else {
    var rtemp = Number.POSITIVE_INFINITY;
    if ( upper ) {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = 1; i <= j - 1; i ++ ) {
          temp.setValue( zero );
          for ( l = 1; l <= k; l ++ ) {
            temp.plusEquals( ComplexMath.times( ComplexMath.conj(
              A[ ioffA + l - 1 + ( i - 1 ) * lda ] ) , 
              A[ ioffA + l - 1 + ( j - 1 ) * lda ] ) );
          }
          if ( beta == 0. ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
              ComplexMath.timesNumber( temp , alpha ) );
          } else {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
              ComplexMath.plus(
              ComplexMath.timesNumber( temp , alpha ) , 
              ComplexMath.timesNumber(
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ] , beta ) ) );
          }
        }
        rtemp = 0.;
        for ( l = 1; l <= k; l ++ ) {
          rtemp += ComplexMath.norm(
            A[ ioffA + l - 1 + ( j - 1 ) * lda ] );
        }
        if ( beta == 0. ) {
          C[ ioffC + j - 1 + ( j - 1 ) * ldc ].setValue(
            new Complex( alpha * rtemp , 0. ) ); 
        } else {
          C[ ioffC + j - 1 + ( j - 1 ) * ldc ].setValue(
            new Complex( alpha * rtemp + beta *
            C[ ioffC + j - 1 + ( j - 1 ) * ldc ].getReal() , 0. ) );
        }
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        rtemp = 0.;
        for ( l = 1; l <= k; l ++ ) {
          rtemp += ComplexMath.norm(
            A[ ioffA + l - 1 + ( j - 1 ) * lda ] );
        }
        if ( beta == 0. ) {
          C[ ioffC + j - 1 + ( j - 1 ) * ldc ].setValue(
            new Complex( alpha * rtemp , 0. ) );
        } else {
          C[ ioffC + j - 1 + ( j - 1 ) * ldc ].setValue(
            new Complex( alpha * rtemp + beta *
            C[ ioffC + j - 1 + ( j - 1 ) * ldc ].getReal() , 0. ) ); 
        }
        for ( i = j + 1; i <= n; i ++ ) {
          temp.setValue( zero );
          for ( l = 1; l <= k; l ++ ) {
            temp.plusEquals( ComplexMath.times(
              ComplexMath.conj(
              A[ ioffA + l - 1 + ( i - 1 ) * lda ] ) ,
              A[ ioffA + l - 1 + ( j - 1 ) * lda ] ) );
          }
          if ( beta == 0. ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
              ComplexMath.timesNumber( temp , alpha ) );
          } else {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
              ComplexMath.plus(
              ComplexMath.timesNumber( temp , alpha ) , 
              ComplexMath.timesNumber(
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ] , beta ) ) );
          }
        }
      }
    }
  }
}
//*************************************************************************
//  C <-- A alpha B^T + B alpha A^T + C beta or
//  C <-- A^T alpha B + B^T alpha A + C beta, C symmetric
Blas3.dsyr2k = function( uplo , trans , n , k , alpha , A , lda , B ,
ldb , beta , C , ldc, ioffA, ioffB, ioffC ) {
  var nrowa =  ( trans.charAt(0).toUpperCase() == 'N' ? n : k );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var info = 0;
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( trans.charAt(0).toUpperCase() != 'N' 
  && trans.charAt(0).toUpperCase() != 'T' 
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info = 2;
  } else if ( n < 0 ) info = 3;
  else if ( k < 0 ) info = 4;
  else if ( lda < Math.max( 1 , nrowa ) ) info = 7;
  else if ( ldb < Math.max( 1 , nrowa ) ) info = 9;
  else if ( ldc < Math.max( 1 , n ) ) info = 12;
  if ( info != 0 ) {
    Blas2.xerbla( 'dsyr2k' , info );
    return;
  }

  if ( n == 0 || ( ( alpha == 0. || k == 0 ) && beta == 1. ) ) return;
  var j = -1;
  var i = -1;
  if ( alpha == 0. ) {
    if ( upper ) {
      if ( beta == 0. ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= j; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ] = 0.;
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= j; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ] *= beta;
          }
        }
      }
    } else {
      if ( beta == 0. ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = j; i <= n; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ] = 0.;
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = j; i <= n; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ] *= beta;
          }
        }
      }
    }
  } else {
    var l = -1;
    var temp1 = Number.POSITIVE_INFINITY;
    var temp2 = Number.POSITIVE_INFINITY;
    if ( trans.charAt(0).toUpperCase() == 'N' ) {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          if ( beta == 0. ) {
            for ( i = 1; i <= j; i ++ ) {
              C[ioffC + i - 1 + ( j - 1 ) * ldc ] = 0.;
            }
          } else if ( beta != 1. ) {
            for ( i = 1; i <= j; i ++ ) {
              C[ioffC + i - 1 + ( j - 1 ) * ldc ] *= beta;
            }
          }
          for ( l = 1; l <= k; l ++ ) {
            if ( A[ ioffA + j - 1 + ( l - 1 ) * lda ] != 0.
            || B[ ioffB + j - 1 + ( l - 1 ) * ldb ] != 0.) {
              temp1 = alpha * B[ ioffB + j - 1 + ( l - 1 ) * lda ];
              temp2 = alpha * A[ ioffA + j - 1 + ( l - 1 ) * lda ];
              for ( i = 1; i <= j; i ++ ) {
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ] +=
                  temp1 * A[ ioffA + i - 1 + ( l - 1 ) * lda ]
                  + temp2 * B[ ioffB + i - 1 + ( l - 1 ) * ldb ];
              }
            }
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          if ( beta == 0. ) { 
            for ( i = j; i <= n; i ++ ) {
              C[ioffC + i - 1 + ( j - 1 ) * ldc ] = 0.;
            }
          } else if ( beta != 1. ) {
            for ( i = j; i <= n; i ++ ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ] *= beta;
            }
          }
          for ( l = 1; l <= k; l ++ ) {
            if ( A[ ioffA + j - 1 + ( l - 1 ) * lda ] != 0.
            || B[ ioffB + j - 1 + ( l - 1 ) * ldb ] != 0.) { 
              temp1 = alpha * B[ ioffB + j - 1 + ( l - 1 ) * lda ];
              temp2 = alpha * A[ ioffA + j - 1 + ( l - 1 ) * lda ];
              for ( i = j; i <= n; i ++ ) {
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ] +=
                  temp1 * A[ ioffA + i - 1 + ( l - 1 ) * lda ]
                  + temp2 * B[ ioffB + i - 1 + ( l - 1 ) * ldb ];
              }
            }
          }
        }
      }
    } else {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 0; i <= j; i ++ ) {
            temp1 = 0.;
            temp2 = 0.;
            for ( l = 1; l <= k; l ++ ) {
              temp1 += A[ ioffA + l - 1 + ( i - 1 ) * lda ]
                     * B[ ioffB + l - 1 + ( j - 1 ) * lda ];
              temp2 += B[ ioffB + l - 1 + ( i - 1 ) * lda ]
                     * A[ ioffA + l - 1 + ( j - 1 ) * lda ];
            }
            if ( beta == 0. ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ] =
                alpha * temp1 + alpha * temp2;
            } else {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ] =
                beta * C[ ioffC + i - 1 + ( j - 1 ) * ldc ]
                + alpha * temp1 + alpha * temp2;
            }
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = j; i <= n; i ++ ) {
            temp1 = 0.;
            temp2 = 0.;
            for ( l = 1; l <= k; l ++ ) {
              temp1 += A[ ioffA + l - 1 + ( i - 1 ) * lda ]
                * B[ ioffB + l - 1 + ( j - 1 ) * lda ];
              temp2 += B[ ioffB + l - 1 + ( i - 1 ) * lda ]
                * A[ ioffA + l - 1 + ( j - 1 ) * lda ];
            }
            if ( beta == 0. ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ] =
                alpha * temp1 + alpha * temp2;
            } else {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ] =
                beta * C[ ioffC + i - 1 + ( j - 1 ) * ldc ]
                + alpha * temp1 + alpha * temp2;
            }
          }
        }
      }
    }
  }
}
Blas3.zsyr2k = function( uplo , trans , n , k , alpha , A , lda , B ,
ldb , beta , C , ldc, ioffA, ioffB, ioffC ) {
  var nrowa =  ( trans.charAt(0).toUpperCase() == 'N' ? n : k );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var info = 0;
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( trans.charAt(0).toUpperCase() != 'N' 
  && trans.charAt(0).toUpperCase() != 'T' ) {
    info = 2;
  } else if ( n < 0 ) info = 3;
  else if ( k < 0 ) info = 4;
  else if ( lda < Math.max( 1 , nrowa ) ) info = 7;
  else if ( ldb < Math.max( 1 , nrowa ) ) info = 9;
  else if ( ldc < Math.max( 1 , n ) ) info = 12;
  if ( info != 0 ) {
    Blas2.xerbla( 'zsyr2k' , info );
    return;
  }

  var zero = new Complex( 0., 0. );
  var one = new Complex( 1., 0. );
  if ( n == 0 ||
  ( ( alpha.equals( zero ) || k != 0 ) && beta.equals( one ) ) ) {
    return;
  }
  var j = -1;
  var i = -1;
  if ( alpha.equals( zero ) ) {
    if ( upper ) {
      if ( beta.equals( zero ) ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= j; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= j; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].timesEquals( beta );
          }
        }
      }
    } else {
      if ( beta.equals( zero ) ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = j; i <= n; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = j; i <= n; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].timesEquals( beta );
          }
        }
      }
    }
  } else {
    var l = -1;
    var temp1 = new Complex();
    var temp2 = new Complex();
    if ( trans.charAt(0).toUpperCase() == 'N' ) {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          if ( beta.equals( zero ) ) {
            for ( i = 1; i <= j; i ++ ) {
              C[ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
            }
          } else if ( ! beta.equals( one ) ) {
            for ( i = 1; i <= j; i ++ ) {
              C[ioffC + i - 1 + ( j - 1 ) * ldc ].timesEquals( beta );
            }
          }
          for ( l = 1; l <= k; l ++ ) {
            if ( ! A[ ioffA + j - 1 + ( l - 1 ) * lda ].equals( zero )
            || ! B[ ioffB + j - 1 + ( l - 1 ) * ldb ].equals( zero ) )
            {
              temp1.setValue( ComplexMath.times( alpha ,
                B[ ioffB + j - 1 + ( l - 1 ) * lda ] ) );
              temp2.setValue( ComplexMath.times( alpha ,
                A[ ioffA + j - 1 + ( l - 1 ) * lda ] ) );
              for ( i = 1; i <= j; i ++ ) {
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ].plusEquals(
                  ComplexMath.plus(
                  ComplexMath.times( temp1 ,
                  A[ ioffA + i - 1 + ( l - 1 ) * lda ] ) ,
                  ComplexMath.times( temp2 ,
                  B[ ioffB + i - 1 + ( l - 1 ) * ldb ] ) ) );
              }
            }
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          if ( beta.equals( zero ) ) { 
            for ( i = j; i <= n; i ++ ) {
              C[ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
            }
          } else if ( ! beta.equals( one ) ) {
            for ( i = j; i <= n; i ++ ) {
              C[ ioffC + i - 1 + ( j - 1) * ldc ].timesEquals( beta );
            }
          }
          for ( l = 1; l <= k; l ++ ) {
            if ( ! A[ ioffA + j - 1 + ( l - 1 ) * lda ].equals( zero )
            || ! B[ ioffB + j - 1 + ( l - 1 ) * ldb ].equals( zero ) )
            { 
              temp1.setValue( ComplexMath.times( alpha ,
                B[ ioffB + j - 1 + ( l - 1 ) * lda ] ) );
              temp2.setValue( ComplexMath.times( alpha ,
                A[ ioffA + j - 1 + ( l - 1 ) * lda ] ) );
              for ( i = j; i <= n; i ++ ) {
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ].plusEquals(
                  ComplexMath.plus(
                  ComplexMath.times( temp1 ,
                  A[ ioffA + i - 1 + ( l - 1 ) * lda ] ) ,
                  ComplexMath.times( temp2 ,
                  B[ ioffB + i - 1 + ( l - 1 ) * ldb ] ) ) );
              }
            }
          }
        }
      }
    } else {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= j; i ++ ) {
            temp1.setValue( zero );
            temp2.setValue( zero );
            for ( l = 1; l <= k; l ++ ) {
              temp1.plusEquals( ComplexMath.times(
                A[ ioffA + l - 1 + ( i - 1 ) * lda ] ,
                B[ ioffB + l - 1 + ( j - 1 ) * lda ] ) );
              temp2.plusEquals( ComplexMath.times(
                B[ ioffB + l - 1 + ( i - 1 ) * lda ] ,
                A[ ioffA + l - 1 + ( j - 1 ) * lda ] ) );
            }
            if ( beta.equals( zero ) ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                ComplexMath.plus(
                ComplexMath.times( alpha , temp1 ) ,
                ComplexMath.times( alpha , temp2 ) ) );
            } else {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                ComplexMath.plus(
                ComplexMath.times( alpha , temp1 ) ,
                ComplexMath.plus( 
                  ComplexMath.times( alpha , temp2 ) ,
                  ComplexMath.times( beta ,
                  C[ ioffC + i - 1 + ( j - 1 ) * ldc ] ) ) ) );
            }
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = j; i <= n; i ++ ) {
            temp1.setValue( zero );
            temp2.setValue( zero );
            for ( l = 1; l <= k; l ++ ) {
              temp1.plusEquals( ComplexMath.times(
                A[ ioffA + l - 1 + ( i - 1 ) * lda ] ,
                B[ ioffB + l - 1 + ( j - 1 ) * lda ] ) );
              temp2.plusEquals( ComplexMath.times(
                B[ ioffB + l - 1 + ( i - 1 ) * lda ] ,
                A[ ioffA + l - 1 + ( j - 1 ) * lda ] ) );
            }
            if ( beta.equals( zero ) ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                ComplexMath.plus(
                ComplexMath.times( alpha , temp1 ) ,
                ComplexMath.times( alpha , temp2 ) ) );
            } else {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
                ComplexMath.plus( ComplexMath.times( alpha , temp1 ) ,
                ComplexMath.plus( ComplexMath.times( alpha , temp2 ) ,
                  ComplexMath.times( beta ,
                  C[ ioffC + i - 1 + ( j - 1 ) * ldc ] ) ) ) );
            }
          }
        }
      }
    }
  }
}
Blas3.zher2k = function( uplo , trans , n , k , alpha , A , lda , B ,
ldb , beta , C , ldc, ioffA, ioffB, ioffC ) {
  var nrowa =  ( trans.charAt(0).toUpperCase() == 'N' ? n : k );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var info = 0;
  if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) {
    info = 1;
  } else if ( trans.charAt(0).toUpperCase() != 'N' 
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info = 2;
  } else if ( n < 0 ) info = 3;
  else if ( k < 0 ) info = 4;
  else if ( lda < Math.max( 1 , nrowa ) ) info = 7;
  else if ( ldb < Math.max( 1 , nrowa ) ) info = 9;
  else if ( ldc < Math.max( 1 , n ) ) info = 12;
  if ( info != 0 ) {
    Blas2.xerbla( 'zher2k' , info );
    return;
  }

  var zero = new Complex( 0. , 0. );
  if ( n == 0 || ( ( alpha.equals( zero ) || k == 0 ) && beta == 1. ) )
  {
    return;
  }
  var j = -1;
  var i = -1;
  if ( alpha.equals( zero ) ) {
    if ( upper ) {
      if ( beta == 0. ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= j; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= j; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].timesNumberEquals(
              beta );
          }
          C[ ioffC + j - 1 + ( j - 1 ) * ldc ].pureReal();
        }
      }
    } else {
      if ( beta == 0. ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = j; i <= n; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = j; i <= n; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].timesNumberEquals(
              beta );
          }
          C[ ioffC + j - 1 + ( j - 1 ) * ldc ].pureReal();
        }
      }
    }
    return;
  }
  var l = -1;
  var temp1 = new Complex();
  var temp2 = new Complex();
  if ( trans.charAt(0).toUpperCase() == 'N' ) {
    if ( upper ) {
      for ( j = 1; j <= n; j ++ ) {
        if ( beta == 0. ) {
          for ( i = 1; i <= j; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
          }
        } else if ( beta != 1. ) {
          for ( i = 1; i <= j; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].timesNumberEquals(
              beta );
          }
          C[ ioffC + j - 1 + ( j - 1 ) * ldc ].pureReal();
        } else {
          C[ ioffC + j - 1 + ( j - 1 ) * ldc ].pureReal();
        }
        for ( l = 1; l <= k; l ++ ) {
          if ( ! A[ ioffA + j - 1 + ( l - 1 ) * lda ].equals( zero )
          && ! B[ ioffB + j - 1 + ( l - 1 ) * ldb ].equals( zero ) )
          {
            temp1.setValue( ComplexMath.times( alpha ,
              ComplexMath.conj(
              B[ ioffB + j - 1 + ( l - 1 ) * lda ] ) ) );
            temp2.setValue( ComplexMath.conj( ComplexMath.times(
              alpha , A[ ioffA + j - 1 + ( l - 1 ) * lda ] ) ) );
            for ( i = 1; i <= j; i ++ ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].plusEquals(
                ComplexMath.plus( ComplexMath.times( temp1 ,
                A[ ioffA + i - 1 + ( l - 1 ) * lda ] ) ,
                ComplexMath.times( temp2 ,
                B[ ioffB + i - 1 + ( l - 1 ) * lda ] ) ) );
            }
            C[ ioffC + j - 1 + ( j - 1 ) * ldc ].pureReal();
          }
        }
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        if ( beta == 0. ) { 
          for ( i = j; i <= n; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue( zero );
          }
        } else if ( beta != 1. ) {
          for ( i = j; i <= n; i ++ ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].timesNumberEquals(
              beta );
          }
          C[ ioffC + j - 1 + ( j - 1 ) * ldc ].pureReal();
        } else {
          C[ ioffC + j - 1 + ( j - 1 ) * ldc ].pureReal();
        }
        for ( l = 1; l <= k; l ++ ) {
          if ( ! A[ ioffA + j - 1 + ( l - 1 ) * lda ].equals( zero )
          || ! B[ ioffB + j - 1 + ( l - 1 ) * ldb ].equals(zero) ) { 
            temp1.setValue( ComplexMath.times( alpha , ComplexMath.conj(
              B[ ioffB + j - 1 + ( l - 1 ) * lda ] ) ) );
            temp2.setValue( ComplexMath.conj( ComplexMath.times(
              alpha , A[ ioffA + j - 1 + ( l - 1 ) * lda ] ) ) );
            for ( i = j; i <= n; i ++ ) {
              C[ ioffC + i - 1 + ( j - 1 ) * ldc ].plusEquals(
                ComplexMath.plus( ComplexMath.times( temp1 ,
                A[ ioffA + i - 1 + ( l - 1 ) * lda ] ) ,
                ComplexMath.times( temp2 ,
                B[ ioffB + i - 1 + ( l - 1 ) * lda ] ) ) );
            }
            C[ ioffC + j - 1 + ( j - 1 ) * ldc ].pureReal();
          }
        }
      }
    }
  } else {
    if ( upper ) {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = 1; i <= j; i ++ ) {
          temp1.setValue( zero );
          temp2.setValue( zero );
          for ( l = 1; l <= k; l ++ ) {
            temp1.plusEquals( ComplexMath.times( ComplexMath.conj(
              A[ ioffA + l - 1 + ( i - 1 ) * lda ] ) ,
              B[ ioffB + l - 1 + ( j - 1 ) * lda ] ) );
            temp2.plusEquals( ComplexMath.times( ComplexMath.conj(
              B[ ioffB + l - 1 + ( i - 1 ) * lda ] ) ,
              A[ ioffA + l - 1 + ( j - 1 ) * lda ] ) );
          }
          if ( beta == 0. ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
              ComplexMath.plus( ComplexMath.times( alpha , temp1 ) ,
              ComplexMath.times( ComplexMath.conj( alpha ),temp2)) );
          } else {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
              ComplexMath.plus( ComplexMath.plus(
                ComplexMath.times( alpha , temp1 ) , 
                ComplexMath.times(ComplexMath.conj(alpha),temp2)) ),
              ComplexMath.timesNumber(
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ] , beta ) );
          }
        }
        C[ ioffC + j - 1 + ( j - 1 ) * ldc ].pureReal();
      }
    } else {
      for ( j = 1; j <= n; j ++ ) {
        for ( i = j; i <= n; i ++ ) {
          temp1.setValue( zero );
          temp2.setValue( zero );
          for ( l = 1; l <= k; l ++ ) {
            temp1.plusEquals( ComplexMath.times( ComplexMath.conj(
              A[ ioffA + l - 1 + ( i - 1 ) * lda ] ) ,
              B[ ioffB + l - 1 + ( j - 1 ) * lda ] ) );
            temp2.plusEquals( ComplexMath.times( ComplexMath.conj(
              B[ ioffB + l - 1 + ( i - 1 ) * lda ] ) ,
              A[ ioffA + l - 1 + ( j - 1 ) * lda ] ) );
          }
          if ( beta == 0. ) {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
              ComplexMath.plus( ComplexMath.times( alpha , temp1 ) ,
              ComplexMath.times(ComplexMath.conj(alpha),temp2 ) ) );
              ;
          } else {
            C[ ioffC + i - 1 + ( j - 1 ) * ldc ].setValue(
              ComplexMath.plus( ComplexMath.plus(
                ComplexMath.times( alpha , temp1 ) , 
                ComplexMath.times(ComplexMath.conj(alpha),temp2)),
              ComplexMath.timesNumber(
                C[ ioffC + i - 1 + ( j - 1 ) * ldc ] , beta ) ) );
          }
        }
        C[ ioffC + j - 1 + ( j - 1 ) * ldc ].pureReal();
      }
    }
  }
}
//*************************************************************************
//  B <-- op( A ) alpha B or B <-- B alpha op( A ), op(A) = A or A^T,
//  A triangular
Blas3.dtrmm = function( side , uplo , trans , diag , m , n , alpha , A ,
lda , B , ldb, ioffA, ioffB ) {
  var lside = ( side.charAt(0).toUpperCase() == 'L' );
  var nrowa = ( lside ? m : n );
  var nounit = ( diag.charAt(0).toUpperCase() == 'N' );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var info = 0;
  if ( ! lside && side.charAt(0).toUpperCase() != 'R' ) info = 1;
  else if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) info = 2;
  else if ( trans.charAt(0).toUpperCase() != 'N' 
  && trans.charAt(0).toUpperCase() != 'T' 
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info = 3;
  } else if ( diag.charAt(0).toUpperCase() != 'U' &&
  diag.charAt(0).toUpperCase() != 'N' ) {
    info = 4;
  } else if ( m < 0 ) info = 5;
  else if ( n < 0 ) info = 6;
  else if ( lda < Math.max( 1 , nrowa ) ) info = 9;
  else if ( ldb < Math.max( 1 , m ) ) info = 11;
  if ( info != 0 ) {
    Blas2.xerbla( 'dtrmm' , info );
    return;
  }

  if ( n == 0 ) return;
  var i = -1;
  var j = -1;
  if ( alpha == 0. ) {
    for ( j = 1; j <= n; j ++ ) {
      for ( i = 1; i <= m; i ++ ) {
         B[ ioffB + i - 1 + ( j - 1 ) * ldb ] = 0.;
      }
    }
    return;
  }
  var k = -1;
  var temp = Number.POSITIVE_INFINITY;
  if ( lside ) {
    if ( trans.charAt(0).toUpperCase() == 'N' ) {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( k = 1; k <= m; k ++ ) {
            if ( B[ ioffB + k - 1 + ( j - 1 ) * ldb ] != 0. ) {
              temp = alpha * B[ ioffB + k - 1 + ( j - 1 ) * ldb ];
              for ( i = 1; i <= k - 1; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1) * ldb ] +=
                  temp * A[ ioffA + i - 1 + ( k - 1 ) * lda ];
              }
              if ( nounit ) {
                temp *= A[ ioffA + k - 1  + ( k - 1 ) * lda ];
              }
              B[ ioffB + k - 1 + ( j - 1 ) * ldb ] = temp;
            }
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( k = m; k >= 1; k -- ) {
            if ( B[ ioffB + k - 1 + ( j - 1 ) * ldb ] != 0. ) {
              temp = alpha * B[ ioffB + k - 1 + ( j - 1 ) * ldb ];
              B[ ioffB + k - 1 + ( j - 1 ) * ldb ] = temp;
              if ( nounit ) {
                B[ ioffB + k - 1 + ( j - 1 ) * ldb ] *=
                  A[ ioffA + k - 1 + ( k - 1 ) * lda ];
              }
              for ( i = k + 1; i <= m; i ++ ) {
                B[ ioffB + i - 1 +  ( j - 1 ) * ldb ] +=
                  temp * A[ ioffA + i - 1 + ( k - 1 ) * lda ];
              }
            }
          }
        }
      }
    } else {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = m; i >= 1; i -- ) {
            temp = B[ ioffB + i - 1 + ( j - 1 ) * ldb ];
            if ( nounit ) temp *= A[ ioffA + i - 1 + ( i - 1 ) * lda ];
            for ( k = 1; k <= i - 1; k ++ ) {
              temp += A[ ioffA + k - 1 +  ( i - 1 ) * lda ]
                    * B[ ioffB + k - 1 + ( j - 1 ) * ldb ];
            }
            B[ ioffB + i - 1 + ( j - 1 ) * ldb ] = alpha * temp;
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= m; i ++ ) {
            temp = B[ ioffB + i - 1 + ( j - 1 ) * ldb ];
            if ( nounit ) temp *= A[ ioffA + i - 1 + ( i - 1 ) * lda ];
            for ( k = i + 1; k <= m; k ++ ) {
              temp += A[ ioffA + k - 1 + ( i - 1 ) * lda ]
                    * B[ ioffB + k - 1 + ( j - 1 ) * ldb ];
            }
            B[ ioffB + i - 1 + ( j - 1 ) * ldb ] = alpha * temp;
          }
        }
      }
    }
  } else {
    if ( trans.charAt(0).toUpperCase() == 'N' ) {
      if ( upper ) {
        for ( j = n; j >= 1; j -- ) {
          temp = alpha;
          if ( nounit ) temp *= A[ ioffA + j - 1 + ( j - 1 ) * lda ];
          for ( i = 1; i <= m; i ++ ) {
            B[ ioffB + i - 1 + ( j - 1 ) * ldb ] *= temp;
          }
          for ( k = 1; k <= j - 1; k ++ ) {
            if ( A[ ioffA + k - 1 + ( j - 1 ) * lda ] != 0. ) {
              temp = alpha * A[ ioffA + k - 1 + ( j - 1 ) * lda ];
              for ( i = 1; i <= m; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ] +=
                  temp * B[ ioffB + i - 1 + ( k - 1 ) * ldb ];
              }
            }
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          temp = alpha;
          if ( nounit ) temp *= A[ ioffA + j - 1 + ( j - 1 ) * lda ];
          for ( i = 1; i <= m; i ++ ) {
            B[ ioffB + i - 1 + ( j - 1 ) * ldb ] *= temp;
          }
          for ( k = j + 1; k <= n; k ++ ) {
            if ( A[ ioffA + k - 1 + ( j - 1 ) * lda ] != 0. ) {
              temp = alpha * A[ ioffA + k - 1 + ( j - 1 ) * lda ];
              for ( i = 1; i <= m; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ] +=
                  temp * B[ ioffB + i - 1 + ( k - 1 ) * ldb ];
              }
            }
          }
        }
      }
    } else {
      if ( upper ) {
        for ( k = 1; k <= n; k ++ ) {
          for ( j = 1; j <= k - 1; j ++ ) {
            if ( A[ ioffA + j - 1 + ( k - 1 ) * lda ] != 0. ) {
              temp = alpha * A[ ioffA + j - 1 + ( k - 1 ) * lda ];
              for ( i = 1; i <= m; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ] +=
                  temp * B[ ioffB + i - 1 + ( k - 1 ) * ldb ];
              }
            }
          }
          temp = alpha;
          if ( nounit ) temp *= A[ ioffA + k - 1 + ( k - 1 ) * lda ];
          if ( temp != 1. ) {
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( k - 1 ) * ldb ] *= temp;
            }
          }
        }
      } else {
        for ( k = n; k >= 1; k -- ) {
          for ( j = k + 1; j <= n; j ++ ) {
            if ( A[ ioffA + j - 1 + ( k - 1 ) * lda ] != 0. ) {
              temp = alpha * A[ ioffA + j - 1 + ( k - 1 ) * lda ];
              for ( i = 1; i <= m; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ] +=
                  temp * B[ ioffB + i - 1 + ( k - 1 ) * ldb ];
              }
            }
          }
          temp = alpha;
          if ( nounit ) temp *= A[ ioffA + k - 1 + ( k - 1 ) * lda ];
          if ( temp != 1. ) {
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( k - 1 ) * ldb ] *= temp;
            }
          }
        }
      }
    }
  }
}
Blas3.ztrmm = function( side , uplo , trans , diag , m , n , alpha , A ,
lda , B , ldb, ioffA, ioffB ) {
  var lside = ( side.charAt(0).toUpperCase() == 'L' );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var nrowa = ( lside ? m : n );
  var noconj = ( trans.charAt(0).toUpperCase() == 'T' );
  var nounit = ( diag.charAt(0).toUpperCase() == 'N' );
  var info = 0;
  if ( ! lside && side.charAt(0).toUpperCase() != 'R' ) info = 1;
  else if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) info = 2;
  else if ( trans.charAt(0).toUpperCase() != 'N' 
  && trans.charAt(0).toUpperCase() != 'T' 
  && trans.charAt(0).toUpperCase() != 'C' ) {
    info = 3;
  } else if ( diag.charAt(0).toUpperCase() != 'U' &&
  diag.charAt(0).toUpperCase() != 'N' ) {
    info = 4;
  } else if ( m < 0 ) info = 5;
  else if ( n < 0 ) info = 6;
  else if ( lda < Math.max( 1 , nrowa ) ) info = 9;
  else if ( ldb < Math.max( 1 , m ) ) info = 11;
  if ( info != 0 ) {
    Blas2.xerbla( 'ztrmm' , info );
    return;
  }

  if ( n == 0 ) return;
  var i = -1;
  var j = -1;
  var zero = new Complex( 0., 0. );
  if ( alpha.equals( zero ) ) {
    for ( j = 1; j <= n; j ++ ) {
      for ( i = 1; i <= m; i ++ ) {
         B[ ioffB + i - 1 + ( j - 1 ) * ldb ].setValue( zero );
      }
    }
    return;
  }
  var k = -1;
  var temp = new Complex();
  if ( lside ) {
    if ( trans.charAt(0).toUpperCase() == 'N' ) {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( k = 1; k <= m; k ++ ) {
            if ( ! B[ ioffB + k - 1 + ( j - 1 ) * ldb ].equals( zero )
            ) {
              temp.setValue( ComplexMath.times( alpha ,
                B[ ioffB + k - 1 + ( j - 1 ) * ldb ] ) );
              for ( i = 1; i <= k - 1; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ].plusEquals(
                  ComplexMath.times( temp ,
                    A[ ioffA + i - 1 + ( k - 1 ) * lda ] ) );
              }
              if ( nounit ) {
                temp.timesEquals(
                  A[ ioffA + k - 1 + ( k - 1 ) * lda ] );
              }
              B[ ioffB + k - 1 + ( j - 1 ) * ldb ].setValue( temp );
            }
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( k = m; k >= 1; k -- ) {
            if ( ! B[ ioffB + k - 1 + ( j - 1 ) * ldb ].equals( zero )
            ) {
              temp.setValue( ComplexMath.times( alpha ,
                B[ ioffB + k - 1 + ( j - 1 ) * ldb ] ) );
              B[ ioffB + k - 1 + ( j - 1 ) * ldb ].setValue( temp );
              if ( nounit ) {
                B[ ioffB + k - 1 + ( j - 1 ) * ldb ].timesEquals(
                  A[ ioffA + k - 1 + ( k - 1 ) * lda ] );
              }
              for ( i = k + 1; i <= m; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ].plusEquals(
                  ComplexMath.times( temp ,
                    A[ ioffA + i - 1 + ( k - 1 ) * lda ] ));
              }
            }
          }
        }
      }
    } else {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = m; i >= 1; i -- ) {
            temp.setValue( B[ ioffB + i - 1 + ( j - 1 ) * ldb ] );
            if ( noconj ) {
              if ( nounit ) {
                temp.timesEquals(
                  A[ ioffA + i - 1 + ( i - 1 ) * lda ] );
              }
              for ( k = 1; k <= i - 1; k ++ ) {
                temp.plusEquals( ComplexMath.times(
                  A[ ioffA + k - 1 + ( i - 1 ) * lda ] ,
                  B[ ioffB + k - 1 + ( j - 1 ) * ldb ] ) );
              }
            } else {
              if ( nounit ) {
                temp.timesEquals( ComplexMath.conj(
                  A[ ioffA + i - 1 + ( i - 1 ) * lda ] ) );
              }
              for ( k = 1; k <= i - 1; k ++ ) {
                temp.plusEquals( ComplexMath.times( ComplexMath.conj(
                  A[ ioffA + k - 1 + ( i - 1 ) * lda ] ) ,
                  B[ ioffB + k - 1 + ( j - 1 ) * ldb ] ) );
              }
            }
            B[ ioffB + i - 1 + ( j - 1 ) * ldb ].setValue(
              ComplexMath.times( alpha , temp ) );
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= m; i ++ ) {
            temp.setValue( B[ ioffB + i - 1 + ( j - 1 ) * ldb ] );
            if ( noconj) {
              if ( nounit ) {
                temp.timesEquals(
                  A[ ioffA + i - 1 + ( i - 1 ) * lda ] );
              }
              for ( k = i + 1; k <= m; k ++ ) {
                temp.plusEquals( ComplexMath.times(
                  A[ ioffA + k - 1 + ( i - 1 ) * lda ] ,
                  B[ ioffB + k - 1 + ( j - 1 ) * ldb ] ) );
              }
            } else {
              if ( nounit ) {
                temp.timesEquals( ComplexMath.conj(
                  A[ ioffA + i - 1 + ( i - 1 ) * lda ] ) );
              }
              for ( k = i + 1; k <= m; k ++ ) {
                temp.plusEquals( ComplexMath.times( ComplexMath.conj(
                  A[ ioffA + k - 1 + ( i - 1 ) * lda ] ),
                  B[ ioffB + k - 1 + ( j - 1 ) * ldb ] ) );
              }
            }
            B[ ioffB + i - 1 + ( j - 1 ) * ldb ].setValue(
              ComplexMath.times( alpha , temp ) );
          }
        }
      }
    }
  } else {
    if ( trans.charAt(0).toUpperCase() == 'N' ) {
      if ( upper ) {
        for ( j = n; j >= 1; j -- ) {
          temp.setValue( alpha );
          if ( nounit ) {
            temp.timesEquals( A[ ioffA + j - 1 + ( j - 1 ) * lda ] );
          }
          for ( i = 1; i <= m; i ++ ) {
            B[ ioffB + i - 1 + ( j - 1 ) * ldb ].timesEquals( temp );
          }
          for ( k = 1; k <= j - 1; k ++ ) {
            if ( A[ ioffA + k - 1 + ( j - 1 ) * lda ] != 0. ) {
              temp.setValue( ComplexMath.times( alpha ,
                A[ ioffA + k - 1 + ( j - 1 ) * lda ] ) );
              for ( i = 1; i <= m; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ].plusEquals(
                  ComplexMath.times( temp ,
                    B[ ioffB + i - 1 + ( k - 1 ) * ldb ] ) );
              }
            }
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          temp.setValue( alpha );
          if ( nounit ) {
            temp.timesEquals( A[ ioffA + j - 1 + ( j - 1 ) * lda ] );
          }
          for ( i = 1; i <= m; i ++ ) {
            B[ ioffB + i - 1 + ( j - 1 ) * ldb ].timesEquals( temp );
          }
          for ( k = j + 1; k <= n; k ++ ) {
            if ( A[ ioffA + k - 1 + ( j - 1 ) * lda ] != 0. ) {
              temp.setValue( ComplexMath.times( alpha ,
                A[ ioffA + k - 1 + ( j - 1 ) * lda ] ) );
              for ( i = 1; i <= m; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ].plusEquals(
                  ComplexMath.times( temp ,
                    B[ ioffB + i - 1 + ( k - 1 ) * ldb ] ) );
              }
            }
          }
        }
      }
    } else {
      var one = new Complex( 1., 0. );
      if ( upper ) {
        for ( k = 1; k <= n; k ++ ) {
          for ( j = 1; j <= k - 1; j ++ ) {
            if ( ! A[ ioffA + j - 1 + ( k - 1 ) * lda ].equals( zero )
            ) {
              if ( noconj ) {
                temp.setValue( ComplexMath.times( alpha ,
                  A[ ioffA + j - 1 + ( k - 1 ) * lda ] ) );
              } else {
                temp.setValue( ComplexMath.times( alpha , 
                  ComplexMath.conj (
                  A[ ioffA + j - 1 + ( k - 1 ) * lda ] ) ) );
              }
              for ( i = 1; i <= m; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ].plusEquals(
                  ComplexMath.times( temp ,
                    B[ ioffB + i - 1 + ( k - 1 ) * ldb ] ) );
              }
            }
          }
          temp.setValue( alpha );
          if ( nounit ) {
            if ( noconj ) {
              temp.timesEquals( A[ ioffA + k - 1 + ( k - 1 ) * lda ] );
            } else {
              temp.timesEquals( ComplexMath.conj(
                A[ ioffA + k - 1 + ( k -1 ) * lda ] ) );
            }
          }
          if ( ! temp.equals( one ) ) {
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( k - 1 ) * ldb ].timesEquals( temp );
            }
          }
        }
      } else {
        for ( k = n; k >= 1; k -- ) {
          for ( j = k + 1; j <= n; j ++ ) {
            if ( ! A[ ioffA + j - 1 + ( k - 1 ) * lda ].equals( zero )
            ) {
              if ( noconj ) {
                temp.setValue( ComplexMath.times( alpha ,
                    A[ ioffA + j - 1 + ( k - 1 ) * lda ] ) );
              } else {
                temp.setValue( ComplexMath.times( alpha ,
                  ComplexMath.conj(
                  A[ ioffA + j - 1 + ( k - 1 ) * lda ] ) ) );
              }
              for ( i = 1; i <= m; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ].plusEquals(
                  ComplexMath.times( temp ,
                    B[ ioffB + i - 1 + ( k - 1 ) * ldb ] ) );
              }
            }
          }
          temp.setValue( alpha );
          if ( nounit ) {
            if ( noconj ) {
              temp.timesEquals( A[ ioffA + k - 1 + ( k - 1 ) * lda ] );
            } else {
              temp.timesEquals( ComplexMath.conj(
                A[ ioffA + k - 1 + ( k - 1 ) * lda ] ) );
            }
          }
          if ( ! temp.equals( one ) ) {
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( k - 1 ) * ldb ].timesEquals( temp );
            }
          }
        }
      }
    }
  }
}
//*************************************************************************
//  solve op( A ) X = alpha B or X op( A ) = B alpha , op(A) = A or A^T,
//  A triangular
Blas3.dtrsm = function( side , uplo , transa , diag , m , n , alpha , A ,
lda , B , ldb, ioffA, ioffB ) {
  var lside = ( side.charAt(0).toUpperCase() == 'L' );
  var nrowa = ( lside ? m : n );
  var nounit = ( diag.charAt(0).toUpperCase() == 'N' );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var info = 0;
  if ( ! lside && side.charAt(0).toUpperCase() != 'R' ) info = 1;
  else if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) info = 2;
  else if ( transa.charAt(0).toUpperCase() != 'N' && 
  transa.charAt(0).toUpperCase() != 'T'
  && transa.charAt(0).toUpperCase() != 'C' ) {
    info = 3;
  } else if ( diag.charAt(0).toUpperCase() != 'U' 
  && diag.charAt(0).toUpperCase() != 'N' ) {
    info = 4;
  } else if ( m < 0 ) info = 5;
  else if ( n < 0 ) info = 6;
  else if ( lda < Math.max( 1 , nrowa ) ) info = 9;
  else if ( ldb < Math.max( 1 , m ) ) info = 11;
  if ( info != 0 ) {
    Blas2.xerbla( 'dtrsm' , info );
    return;
  }

  if ( m == 0 || n == 0 ) return;
  var i = -1;
  var j = -1;
  if ( alpha == 0. ) {
    for ( j = 1; j <= n; j ++ ) {
      for ( i = 1; i <= m; i ++ ) {
         B[ ioffB + i - 1 + ( j - 1 ) * ldb ] = 0.;
      }
    }
    return;
  }
  var k = -1;
  var temp = Number.POSITIVE_INFINITY;
  if ( lside ) {
    if ( transa.charAt(0).toUpperCase() == 'N' ) {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          if ( alpha != 1.) {
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ] *= alpha;
            }
          }
          for ( k = m; k >= 1; k -- ) {
            if ( B[ ioffB + k - 1 + ( j - 1 ) * ldb ] != 0. ) {
              if ( nounit ) {
                B[ ioffB + k - 1 + ( j - 1 ) * ldb ] /=
                  A[ ioffA + k - 1 + ( k - 1 ) * lda ];
              }
              for ( i = 1; i <= k - 1; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ] -=
                  B[ ioffB + k - 1 + ( j - 1 ) * ldb]
                  * A[ ioffA + i - 1 + ( k - 1 ) * lda ];
              }
            }
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          if ( alpha != 1. ) {
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ] *= alpha;
            }
          }
          for ( k = 1; k <= m; k ++ ) {
            if ( B[ ioffB + k - 1 + ( j - 1 ) * ldb ] != 0. ) {
              if ( nounit ) {
                B[ ioffB + k - 1 + ( j - 1 ) * ldb ] /=
                  A[ ioffA + k - 1 + ( k - 1 ) * lda ];
              }
              for ( i = k + 1; i <= m; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ] -=
                  B[ ioffB + k - 1 + ( j - 1 ) * ldb ]
                  * A[ ioffA + i - 1 + ( k - 1 ) * lda ];
              }
            }
          }
        }
      }
    } else {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= m; i ++ ) {
            temp = alpha * B[ ioffB + i - 1 + ( j - 1 ) * ldb ];
            for ( k = 1; k <= i - 1; k ++ ) {
              temp -= A[ ioffA + k - 1 + ( i - 1 ) * lda ]
                * B[ ioffB + k - 1 + ( j - 1 ) * ldb ];
            }
            if ( nounit ) temp /= A[ ioffA + i - 1 + ( i - 1 ) * lda ];
            B[ ioffB + i - 1 + ( j - 1 ) * ldb ] = temp;
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = m; i >= 1; i -- ) {
            temp = alpha * B[ ioffB + i - 1  + ( j - 1 ) * ldb ];
            for ( k = i + 1; k <= m; k ++ ) {
              temp -= A[ ioffA + k - 1 + ( i - 1 ) * lda ]
                * B[ ioffB + k - 1 + ( j - 1 ) * ldb ];
            }
            if ( nounit ) temp /= A[ ioffA + i - 1 + ( i - 1 ) * lda ];
            B[ ioffB + i - 1 + ( j - 1 ) * ldb ] = temp;
          }
        }
      }
    }
  } else {
    if ( transa.charAt(0).toUpperCase() == 'N' ) {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          if ( alpha != 1. ) {
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ] *= alpha;
            }
          }
          for ( k = 1; k <= j - 1; k ++ ) {
            if ( A[ ioffA + k - 1 + ( j - 1 ) * lda ] != 0. ) {
              for ( i = 1; i <= m; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ] -=
                  A[ ioffA + k - 1 + ( j - 1 ) * lda ]
                    * B[ ioffB + i - 1 + ( k - 1 ) * ldb ];
              }
            }
          }
          if ( nounit ) {
            temp = 1. / A[ ioffA + j - 1 + ( j - 1 ) * lda ];
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ] *= temp;
            }
          }
        }
      } else {
        for ( j = n; j >= 1; j -- ) {
          if ( alpha != 1. ) {
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ] *= alpha;
            }
          }
          for ( k = j + 1; k <= n; k ++ ) {
            if ( A[ ioffA + k - 1 + ( j - 1 ) * lda ] != 0. ) {
              for ( i = 1; i <= m; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ] -=
                  A[ ioffA + k - 1 + ( j - 1 ) * lda ]
                    * B[ ioffB + i - 1 + ( k - 1 ) * ldb ];
              }
            }
          }
          if ( nounit ) {
            temp = 1. / A[ ioffA + j - 1 + ( j - 1 ) * lda ];
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ] *= temp;
            }
          }
        }
      }
    } else {
      if ( upper ) {
        for ( k = n; k >= 1; k -- ) {
          if ( nounit) {
            temp = 1. / A[ ioffA + k - 1 + ( k - 1 ) * lda ];
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( k - 1 ) * ldb ] *= temp;
            }
          }
          for ( j = 1; j <= k - 1; j ++ ) {
            if ( A[ ioffA + j - 1 + ( k - 1 ) * lda ] != 0. ) {
              temp = A[ ioffA + j - 1 + ( k - 1 ) * lda ];
              for ( i = 1; i <= m; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ] -=
                  temp * B[ ioffB + i - 1 + ( k - 1 ) * ldb ];
              }
            }
          }
          if ( alpha != 1. ) {
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( k - 1 ) * ldb ] *= alpha;
            }
          }
        }
      } else {
        for ( k = 1; k <= n; k ++ ) {
          if ( nounit ) {
            temp = 1. / A[ ioffA + k - 1 + ( k - 1 ) * lda ];
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( k - 1 ) * ldb ] *= temp;
            }
          }
          for ( j = k + 1; j <= n; j ++ ) {
            if ( A[ ioffA + j - 1 + ( k - 1 ) * lda ] != 0. ) {
              temp = A[ ioffA + j - 1 + ( k - 1 ) * lda ];
              for ( i = 1; i <= m; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ]
                  -= temp * B[ ioffB + i - 1 + ( k - 1 ) * ldb ];
              }
            }
          }
          if ( alpha != 1. ) {
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( k - 1 ) * ldb ] *= alpha;
            }
          }
        }
      }
    }
  }
}
Blas3.ztrsm = function( side , uplo , transa , diag , m , n , alpha , A ,
lda , B , ldb, ioffA, ioffB ) {
  var lside = ( side.charAt(0).toUpperCase() == 'L' );
  var nrowa = ( lside ? m : n );
  var noconj = ( transa.charAt(0).toUpperCase() == 'T' );
  var nounit = ( diag.charAt(0).toUpperCase() == 'N' );
  var upper = ( uplo.charAt(0).toUpperCase() == 'U' );
  var info = 0;
  if ( ! lside && side.charAt(0).toUpperCase() != 'R' ) info = 1;
  else if ( ! upper && uplo.charAt(0).toUpperCase() != 'L' ) info = 2;
  else if ( transa.charAt(0).toUpperCase() != 'N'
  && transa.charAt(0).toUpperCase() != 'T'
  && transa.charAt(0).toUpperCase() != 'C' ) {
    info = 3;
  } else if ( diag.charAt(0).toUpperCase() != 'U' 
  && diag.charAt(0).toUpperCase() != 'N' ) {
    info = 4;
  } else if ( m < 0 ) info = 5;
  else if ( n < 0 ) info = 6;
  else if ( lda < Math.max( 1 , nrowa ) ) info = 9;
  else if ( ldb < Math.max( 1 , m ) ) info = 11;
  if ( info != 0 ) {
    Blas2.xerbla( 'ztrsm' , info );
    return;
  }

  if ( n == 0 ) return;
  var i = -1;
  var j = -1;
  var zero = new Complex( 0., 0. );
  var one = new Complex( 1., 0. );
  if ( alpha.equals( zero ) ) {
    for ( j = 1; j <= n; j ++ ) {
      for ( i = 1; i <= m; i ++ ) {
         B[ ioffB + i - 1 + ( j - 1 ) * ldb ].setValue( zero );
      }
    }
    return;
  }
  var k = -1;
  var temp = new Complex();
  if ( lside ) {
    if ( transa.charAt(0).toUpperCase() == 'N' ) {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          if ( !alpha.equals( one ) ) {
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ].timesEquals(
                alpha );
            }
          }
          for ( k = m; k >= 1; k -- ) {
            if ( ! B[ ioffB + k - 1 + ( j - 1 ) * ldb ].equals( zero )
            ) {
              if ( nounit ) {
                B[ ioffB + k - 1 + ( j - 1 ) * ldb ].divideEquals(
                  A[ ioffA + k - 1 + ( k - 1 ) * lda ] );
              }
              for ( i = 1; i <= k - 1; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ].minusEquals(
                  ComplexMath.times(
                    B[ ioffB + k - 1 + ( j - 1 ) * ldb ] ,
                    A[ ioffA + i - 1 + ( k - 1 ) * lda ] ) );
              }
            }
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          if ( ! alpha.equals( one ) ) {
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ].timesEquals(
                alpha );
            }
          }
          for ( k = 1; k <= m; k ++ ) {
            if ( ! B[ ioffB + k - 1 + ( j - 1 ) * ldb ].equals( zero )
            ) {
              if ( nounit ) {
                B[ ioffB + k - 1 + ( j - 1 ) * ldb ].divideEquals(
                  A[ ioffA + k - 1 + ( k - 1 ) * lda ] );
              }
              for ( i = k + 1; i <= m; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ].minusEquals(
                  ComplexMath.times(
                    B[ ioffB + k - 1 + ( j - 1 ) * ldb ] ,
                    A[ ioffA + i - 1 + ( k - 1 ) * lda ] ) );
              }
            }
          }
        }
      }
    } else {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = 1; i <= m; i ++ ) {
            temp.setValue(
              ComplexMath.times( alpha,
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ] ) );
            if ( noconj ) {
              for ( k = 1; k <= i - 1; k ++ ) {
                temp.minusEquals( ComplexMath.times(
                  A[ ioffA + k - 1 + ( i - 1 ) * lda ] ,
                  B[ ioffB + k - 1 + ( j - 1 ) * ldb ] ) );
              }
              if ( nounit ) {
                temp.divideEquals(
                  A[ ioffA + i - 1 + ( i - 1 ) * lda ] );
              }
            } else {
              for ( k = 1; k <= i - 1; k ++ ) {
                temp.minusEquals( ComplexMath.times( ComplexMath.conj(
                  A[ ioffA + k - 1 + ( i - 1 ) * lda ] ) ,
                  B[ ioffB + k - 1 + ( j - 1 ) * ldb ] ) );
              }
              if ( nounit ) {
                temp.divideEquals( ComplexMath.conj(
                  A[ ioffA + i - 1 + ( i - 1 ) * lda ] ) );
              }
            }
            B[ ioffB + i - 1 + ( j - 1 ) * ldb ].setValue( temp );
          }
        }
      } else {
        for ( j = 1; j <= n; j ++ ) {
          for ( i = m; i >= 1; i -- ) {
            temp.setValue( ComplexMath.times( alpha,
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ] ) );
            if ( noconj ) {
              for ( k = i + 1; k <= m; k ++ ) {
                temp.minusEquals( ComplexMath.times(
                  A[ ioffA + k - 1 + ( i - 1 ) * lda ] ,
                  B[ ioffB + k - 1 + ( j - 1 ) * ldb ] ) );
              }
              if ( nounit ) {
                temp.divideEquals(
                  A[ ioffA + i - 1 + ( i - 1 ) * lda ] );
              }
            } else {
              for ( k = i + 1; k <= m; k ++ ) {
                temp.minusEquals( ComplexMath.times(
                  ComplexMath.conj(
                  A[ ioffA + k - 1 + ( i - 1 ) * lda ] ) ,
                  B[ ioffB + k - 1 + ( j - 1 ) * ldb ] ) );
              }
              if ( nounit ) {
                temp.divideEquals( ComplexMath.conj(
                  A[ ioffA + i - 1 + ( i - 1 ) * lda ] ) );
              }
            }
            B[ ioffB + i - 1 + ( j - 1 ) * ldb ].setValue( temp );
          }
        }
      }
    }
  } else {
    if ( transa.charAt(0).toUpperCase() == 'N' ) {
      if ( upper ) {
        for ( j = 1; j <= n; j ++ ) {
          if ( ! alpha.equals( one ) ) {
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ].timesEquals(
                alpha );
            }
          }
          for ( k = 1; k <= j - 1; k ++ ) {
            if ( ! A[ ioffA + k - 1 + ( j - 1 ) * lda ].equals( zero )
            ) {
              for ( i = 1; i <= m; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ].minusEquals(
                  ComplexMath.times(
                  A[ ioffA + k - 1 + ( j - 1 ) * lda ] ,
                  B[ ioffB + i - 1 + ( k - 1 ) * ldb ] ) );
              }
            }
          }
          if ( nounit ) {
            temp.setValue( ComplexMath.divide( one,
              A[ ioffA + j - 1 + ( j - 1 ) * lda ] ) );
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ].timesEquals( temp );
            }
          }
        }
      } else {
        for ( j = n; j >= 1; j -- ) {
          if ( ! alpha.equals( one ) ) {
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ].timesEquals(
                alpha );
            }
          }
          for ( k = j + 1; k <= n; k ++ ) {
            if ( ! A[ ioffA + k - 1 + ( j - 1 ) * lda ].equals( zero )
            ) {
              for ( i = 1; i <= m; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ].minusEquals(
                  ComplexMath.times(
                  A[ ioffA + k - 1 + ( j - 1 ) * lda ] ,
                  B[ ioffB + i - 1 + ( k - 1 ) * ldb ] ) );
              }
            }
          }
          if ( nounit ) {
            temp.setValue(
              ComplexMath.divide( one ,
              A[ ioffA + j - 1 + ( j - 1 ) * lda ] ) );
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( j - 1 ) * ldb ].timesEquals( temp );
            }
         }
        }
      }
    } else {
      if ( upper ) {
        for ( k = n; k >= 1; k -- ) {
          if ( nounit ) {
            if ( noconj ) {
              temp.setValue( ComplexMath.divide( one,
                A[ ioffA + k - 1 + ( k - 1 ) * lda ] ) );
            } else {
              temp.setValue( ComplexMath.divide( one, ComplexMath.conj(
                A[ ioffA + k - 1 + ( k - 1 ) * lda ] ) ) );
            }
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( k - 1 ) * ldb ].timesEquals( temp );
            }
          }
          for ( j = 1; j <= k - 1; j ++ ) {
            if ( ! A[ ioffA + j - 1 + ( k - 1 ) * lda ].equals( zero )
            ) {
              if ( noconj ) {
                temp.setValue( A[ ioffA + j - 1 + ( k - 1 ) * lda ] );
              } else {
                temp.setValue( ComplexMath.conj(
                  A[ ioffA + j - 1 + ( k - 1 ) * lda ] ) );
              }
              for ( i = 1; i <= m; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ].minusEquals(
                  ComplexMath.times( temp ,
                    B[ ioffB + i - 1 + ( k - 1 ) * ldb ] ) );
              }
            }
          }
          if ( ! alpha.equals( one ) ) {
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( k - 1 ) * ldb ].timesEquals(
                alpha );
            }
          }
        }
      } else {
        for ( k = 1; k <= n; k ++ ) {
          if ( nounit ) {
            if ( noconj ) {
              temp.setValue( ComplexMath.divide( one,
                A[ ioffA + k - 1 + ( k - 1 ) * lda ] ) );
            } else {
              temp.setValue( ComplexMath.divide( one, ComplexMath.conj(
                A[ ioffA + k - 1 + ( k - 1 ) * lda ] ) ) );
            }
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( k - 1 ) * ldb ].timesEquals( temp );
            }
          }
          for ( j = k + 1; j <= n; j ++ ) {
            if ( ! A[ ioffA + j - 1 + ( k - 1 ) * lda ].equals( zero )
            ) {
              if ( noconj ) {
                temp.setValue( A[ ioffA + j - 1 + ( k - 1 ) * lda ] );
              } else {
                temp.setValue( ComplexMath.conj(
                  A[ ioffA + j - 1 + ( k - 1 ) * lda ] ) );
              }
              for ( i = 1; i <= m; i ++ ) {
                B[ ioffB + i - 1 + ( j - 1 ) * ldb ].minusEquals(
                  ComplexMath.times( temp ,
                    B[ ioffB + i - 1 + ( k - 1 ) * ldb ] ) );
              }
            }
          }
          if ( ! alpha.equals( one ) ) {
            for ( i = 1; i <= m; i ++ ) {
              B[ ioffB + i - 1 + ( k - 1 ) * ldb ].timesEquals(
                alpha );
            }
          }
        }
      }
    }
  }
}
function testBlas3() {
//test_dgemm();
//test_zgemm();
//test_dsymm();
//test_zsymm();
//test_zhemm();

//test_dsyrk();
//test_zsyrk();
//test_zherk();

//test_dsyr2k();
//test_zsyr2k();
//test_zher2k();

//test_dtrmm();
//test_ztrmm();

  test_dtrsm();
  test_ztrsm();
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dgemm() {
  document.getElementById("debug_textarea").value +=
    "testing dgemm *************" + "\n";
  var ioffA = 1;
  var ioffB = 2;
  var ioffC = 3;
  var A = new Array( ioffA + 8 );
  var B = new Array( ioffB + 20 );
  var C = new Array( ioffC + 16 );
  var lda = 4;
  var ldb = 5;
  var ldc = 4;

//    C = A alpha B + C beta, A 3x2, B 2x4, C 3x4
  var k = 2;
  var m = 3;
  var n = 4;
  var i = -1;
  var j = -1;
  for ( j = 0; j < 2; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < 4; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < 4; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number.POSITIVE_INFINITY;
    }
  }
  A[ ioffA + 0 + 0 * lda ] = 2.;
  A[ ioffA + 1 + 0 * lda ] = 1.;
  A[ ioffA + 2 + 0 * lda ] = 0.;
  A[ ioffA + 0 + 1 * lda ] = 3.;
  A[ ioffA + 1 + 1 * lda ] = 2.;
  A[ ioffA + 2 + 1 * lda ] = 1.;
  B[ ioffB + 0 + 0 * ldb ] = 1.;
  B[ ioffB + 1 + 0 * ldb ] = 0.;
  B[ ioffB + 0 + 1 * ldb ] = 0.;
  B[ ioffB + 1 + 1 * ldb ] = -1.;
  B[ ioffB + 0 + 2 * ldb ] = -1.;
  B[ ioffB + 1 + 2 * ldb ] = -2.;
  B[ ioffB + 0 + 3 * ldb ] = -2.;
  B[ ioffB + 1 + 3 * ldb ] = -3.;
  C[ ioffC + 0 + 0 * ldc ] = 4.;
  C[ ioffC + 1 + 0 * ldc ] = 2.;
  C[ ioffC + 2 + 0 * ldc ] = 0.;
  C[ ioffC + 0 + 1 * ldc ] = 6.;
  C[ ioffC + 1 + 1 * ldc ] = 2.;
  C[ ioffC + 2 + 1 * ldc ] = -2.;
  C[ ioffC + 0 + 2 * ldc ] = 8.;
  C[ ioffC + 1 + 2 * ldc ] = 0.;
  C[ ioffC + 2 + 2 * ldc ] = -8.;
  C[ ioffC + 0 + 3 * ldc ] = 10.;
  C[ ioffC + 1 + 3 * ldc ] = -6.;
  C[ ioffC + 2 + 3 * ldc ] = -22.;
  Blas3.dgemm('N','N',m,n,k,0.,A,lda,B,ldb,1.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dgemm('N','N',m,n,k,0.,A,lda,B,ldb,1.,C,ldc) , C = "
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
    + C[ ioffC + 2 + 3 * ldc ] + "\n";
  C[ ioffC + 0 + 0 * ldc ] = 4.;
  C[ ioffC + 1 + 0 * ldc ] = 2.;
  C[ ioffC + 2 + 0 * ldc ] = 0.;
  C[ ioffC + 0 + 1 * ldc ] = 6.;
  C[ ioffC + 1 + 1 * ldc ] = 2.;
  C[ ioffC + 2 + 1 * ldc ] = -2.;
  C[ ioffC + 0 + 2 * ldc ] = 8.;
  C[ ioffC + 1 + 2 * ldc ] = 0.;
  C[ ioffC + 2 + 2 * ldc ] = -8.;
  C[ ioffC + 0 + 3 * ldc ] = 10.;
  C[ ioffC + 1 + 3 * ldc ] = -6.;
  C[ ioffC + 2 + 3 * ldc ] = -22.;
  Blas3.dgemm('N','N',m,n,k,0.,A,lda,B,ldb,0.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dgemm('N','N',m,n,k,0.,A,lda,B,ldb,0.,C,ldc) , C = "
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
    + C[ ioffC + 2 + 3 * ldc ] + "\n";
  C[ ioffC + 0 + 0 * ldc ] = 4.;
  C[ ioffC + 1 + 0 * ldc ] = 2.;
  C[ ioffC + 2 + 0 * ldc ] = 0.;
  C[ ioffC + 0 + 1 * ldc ] = 6.;
  C[ ioffC + 1 + 1 * ldc ] = 2.;
  C[ ioffC + 2 + 1 * ldc ] = -2.;
  C[ ioffC + 0 + 2 * ldc ] = 8.;
  C[ ioffC + 1 + 2 * ldc ] = 0.;
  C[ ioffC + 2 + 2 * ldc ] = -8.;
  C[ ioffC + 0 + 3 * ldc ] = 10.;
  C[ ioffC + 1 + 3 * ldc ] = -6.;
  C[ ioffC + 2 + 3 * ldc ] = -22.;
  Blas3.dgemm('N','N',m,n,k,0.,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dgemm('N','N',m,n,k,0.,A,lda,B,ldb,2.,C,ldc) , C = "
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
    + C[ ioffC + 2 + 3 * ldc ] + "\n";
  C[ ioffC + 0 + 0 * ldc ] = 4.;
  C[ ioffC + 1 + 0 * ldc ] = 2.;
  C[ ioffC + 2 + 0 * ldc ] = 0.;
  C[ ioffC + 0 + 1 * ldc ] = 6.;
  C[ ioffC + 1 + 1 * ldc ] = 2.;
  C[ ioffC + 2 + 1 * ldc ] = -2.;
  C[ ioffC + 0 + 2 * ldc ] = 8.;
  C[ ioffC + 1 + 2 * ldc ] = 0.;
  C[ ioffC + 2 + 2 * ldc ] = -8.;
  C[ ioffC + 0 + 3 * ldc ] = 10.;
  C[ ioffC + 1 + 3 * ldc ] = -6.;
  C[ ioffC + 2 + 3 * ldc ] = -22.;
  Blas3.dgemm('N','N',m,n,k,3.,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dgemm('N','N',m,n,k,3.,A,lda,B,ldb,2.,C,ldc) , C = "
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
    + C[ ioffC + 2 + 3 * ldc ] + "\n";

//    C = A' alpha B + C beta, A 3x2, B 3x4, C 2x4
  k = 3;
  m = 2;
  n = 4;
  for ( j = 0; j < 2; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < 4; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < 4; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number.POSITIVE_INFINITY;
    }
  }
  A[ ioffA + 0 + 0 * lda ] = 2.;
  A[ ioffA + 1 + 0 * lda ] = 1.;
  A[ ioffA + 2 + 0 * lda ] = 0.;
  A[ ioffA + 0 + 1 * lda ] = 3.;
  A[ ioffA + 1 + 1 * lda ] = 2.;
  A[ ioffA + 2 + 1 * lda ] = 1.;
  B[ ioffB + 0 + 0 * ldb ] = 1.;
  B[ ioffB + 1 + 0 * ldb ] = 0.;
  B[ ioffB + 2 + 0 * ldb ] = -1.;
  B[ ioffB + 0 + 1 * ldb ] = 0.;
  B[ ioffB + 1 + 1 * ldb ] = -1.;
  B[ ioffB + 2 + 1 * ldb ] = -2.;
  B[ ioffB + 0 + 2 * ldb ] = -1.;
  B[ ioffB + 1 + 2 * ldb ] = -2.;
  B[ ioffB + 2 + 2 * ldb ] = -3.;
  B[ ioffB + 0 + 3 * ldb ] = -2.;
  B[ ioffB + 1 + 3 * ldb ] = -3.;
  B[ ioffB + 2 + 3 * ldb ] = -4.;
  C[ ioffC + 0 + 0 * ldc ] = 4.;
  C[ ioffC + 1 + 0 * ldc ] = 2.;
  C[ ioffC + 0 + 1 * ldc ] = 6.;
  C[ ioffC + 1 + 1 * ldc ] = 2.;
  C[ ioffC + 0 + 2 * ldc ] = 8.;
  C[ ioffC + 1 + 2 * ldc ] = 0.;
  C[ ioffC + 0 + 3 * ldc ] = 10.;
  C[ ioffC + 1 + 3 * ldc ] = -6.;
  Blas3.dgemm('T','N',m,n,k,0.,A,lda,B,ldb,1.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dgemm('T','N',m,n,k,0.,A,lda,B,ldb,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + "\n";
  C[ ioffC + 0 + 0 * ldc ] = 4.;
  C[ ioffC + 1 + 0 * ldc ] = 2.;
  C[ ioffC + 0 + 1 * ldc ] = 6.;
  C[ ioffC + 1 + 1 * ldc ] = 2.;
  C[ ioffC + 0 + 2 * ldc ] = 8.;
  C[ ioffC + 1 + 2 * ldc ] = 0.;
  C[ ioffC + 0 + 3 * ldc ] = 10.;
  C[ ioffC + 1 + 3 * ldc ] = -6.;
  Blas3.dgemm('T','N',m,n,k,0.,A,lda,B,ldb,0.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dgemm('T','N',m,n,k,0.,A,lda,B,ldb,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + "\n";
  C[ ioffC + 0 + 0 * ldc ] = 4.;
  C[ ioffC + 1 + 0 * ldc ] = 2.;
  C[ ioffC + 0 + 1 * ldc ] = 6.;
  C[ ioffC + 1 + 1 * ldc ] = 2.;
  C[ ioffC + 0 + 2 * ldc ] = 8.;
  C[ ioffC + 1 + 2 * ldc ] = 0.;
  C[ ioffC + 0 + 3 * ldc ] = 10.;
  C[ ioffC + 1 + 3 * ldc ] = -6.;
  Blas3.dgemm('T','N',m,n,k,0.,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dgemm('T','N',m,n,k,0.,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + "\n";
  C[ ioffC + 0 + 0 * ldc ] = 4.;
  C[ ioffC + 1 + 0 * ldc ] = 2.;
  C[ ioffC + 0 + 1 * ldc ] = 6.;
  C[ ioffC + 1 + 1 * ldc ] = 2.;
  C[ ioffC + 0 + 2 * ldc ] = 8.;
  C[ ioffC + 1 + 2 * ldc ] = 0.;
  C[ ioffC + 0 + 3 * ldc ] = 10.;
  C[ ioffC + 1 + 3 * ldc ] = -6.;
  Blas3.dgemm('T','N',m,n,k,3.,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dgemm('T','N',m,n,k,3.,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + "\n";

//    C = A alpha B' + C beta, A 3x2, B 4x2, C 3x4
  k = 2;
  m = 3;
  n = 4;
  for ( j = 0; j < 2; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < 4; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < 4; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number.POSITIVE_INFINITY;
    }
  }
  A[ ioffA + 0 + 0 * lda ] = 2.;
  A[ ioffA + 1 + 0 * lda ] = 1.;
  A[ ioffA + 2 + 0 * lda ] = 0.;
  A[ ioffA + 0 + 1 * lda ] = 3.;
  A[ ioffA + 1 + 1 * lda ] = 2.;
  A[ ioffA + 2 + 1 * lda ] = 1.;
  B[ ioffB + 0 + 0 * ldb ] = 1.;
  B[ ioffB + 0 + 1 * ldb ] = 0.;
  B[ ioffB + 1 + 0 * ldb ] = 0.;
  B[ ioffB + 1 + 1 * ldb ] = -1.;
  B[ ioffB + 2 + 0 * ldb ] = -1.;
  B[ ioffB + 2 + 1 * ldb ] = -2.;
  B[ ioffB + 3 + 0 * ldb ] = -2.;
  B[ ioffB + 3 + 1 * ldb ] = -3.;
  C[ ioffC + 0 + 0 * ldc ] = 4.;
  C[ ioffC + 1 + 0 * ldc ] = 2.;
  C[ ioffC + 2 + 0 * ldc ] = 0.;
  C[ ioffC + 0 + 1 * ldc ] = 6.;
  C[ ioffC + 1 + 1 * ldc ] = 2.;
  C[ ioffC + 2 + 1 * ldc ] = -2.;
  C[ ioffC + 0 + 2 * ldc ] = 8.;
  C[ ioffC + 1 + 2 * ldc ] = 0.;
  C[ ioffC + 2 + 2 * ldc ] = -8.;
  C[ ioffC + 0 + 3 * ldc ] = 10.;
  C[ ioffC + 1 + 3 * ldc ] = -6.;
  C[ ioffC + 2 + 3 * ldc ] = -22.;
  Blas3.dgemm('N','T',m,n,k,0.,A,lda,B,ldb,1.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dgemm('N','T',m,n,k,0.,A,lda,B,ldb,1.,C,ldc) , C = "
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
    + C[ ioffC + 2 + 3 * ldc ] + "\n";
  C[ ioffC + 0 + 0 * ldc ] = 4.;
  C[ ioffC + 1 + 0 * ldc ] = 2.;
  C[ ioffC + 2 + 0 * ldc ] = 0.;
  C[ ioffC + 0 + 1 * ldc ] = 6.;
  C[ ioffC + 1 + 1 * ldc ] = 2.;
  C[ ioffC + 2 + 1 * ldc ] = -2.;
  C[ ioffC + 0 + 2 * ldc ] = 8.;
  C[ ioffC + 1 + 2 * ldc ] = 0.;
  C[ ioffC + 2 + 2 * ldc ] = -8.;
  C[ ioffC + 0 + 3 * ldc ] = 10.;
  C[ ioffC + 1 + 3 * ldc ] = -6.;
  C[ ioffC + 2 + 3 * ldc ] = -22.;
  Blas3.dgemm('N','T',m,n,k,0.,A,lda,B,ldb,0.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dgemm('N','T',m,n,k,0.,A,lda,B,ldb,0.,C,ldc) , C = "
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
    + C[ ioffC + 2 + 3 * ldc ] + "\n";
  C[ ioffC + 0 + 0 * ldc ] = 4.;
  C[ ioffC + 1 + 0 * ldc ] = 2.;
  C[ ioffC + 2 + 0 * ldc ] = 0.;
  C[ ioffC + 0 + 1 * ldc ] = 6.;
  C[ ioffC + 1 + 1 * ldc ] = 2.;
  C[ ioffC + 2 + 1 * ldc ] = -2.;
  C[ ioffC + 0 + 2 * ldc ] = 8.;
  C[ ioffC + 1 + 2 * ldc ] = 0.;
  C[ ioffC + 2 + 2 * ldc ] = -8.;
  C[ ioffC + 0 + 3 * ldc ] = 10.;
  C[ ioffC + 1 + 3 * ldc ] = -6.;
  C[ ioffC + 2 + 3 * ldc ] = -22.;
  Blas3.dgemm('N','T',m,n,k,0.,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dgemm('N','T',m,n,k,0.,A,lda,B,ldb,2.,C,ldc) , C = "
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
    + C[ ioffC + 2 + 3 * ldc ] + "\n";
  C[ ioffC + 0 + 0 * ldc ] = 4.;
  C[ ioffC + 1 + 0 * ldc ] = 2.;
  C[ ioffC + 2 + 0 * ldc ] = 0.;
  C[ ioffC + 0 + 1 * ldc ] = 6.;
  C[ ioffC + 1 + 1 * ldc ] = 2.;
  C[ ioffC + 2 + 1 * ldc ] = -2.;
  C[ ioffC + 0 + 2 * ldc ] = 8.;
  C[ ioffC + 1 + 2 * ldc ] = 0.;
  C[ ioffC + 2 + 2 * ldc ] = -8.;
  C[ ioffC + 0 + 3 * ldc ] = 10.;
  C[ ioffC + 1 + 3 * ldc ] = -6.;
  C[ ioffC + 2 + 3 * ldc ] = -22.;
  Blas3.dgemm('N','T',m,n,k,3.,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dgemm('N','T',m,n,k,3.,A,lda,B,ldb,2.,C,ldc) , C = "
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
    + C[ ioffC + 2 + 3 * ldc ] + "\n";

//    C = A' alpha B' + C beta, A 3x2, B 4x3, C 2x4
  k = 3;
  m = 2;
  n = 4;
  for ( j = 0; j < 2; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < 4; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < 4; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number.POSITIVE_INFINITY;
    }
  }
  A[ ioffA + 0 + 0 * lda ] = 2.;
  A[ ioffA + 1 + 0 * lda ] = 1.;
  A[ ioffA + 2 + 0 * lda ] = 0.;
  A[ ioffA + 0 + 1 * lda ] = 3.;
  A[ ioffA + 1 + 1 * lda ] = 2.;
  A[ ioffA + 2 + 1 * lda ] = 1.;
  B[ ioffB + 0 + 0 * ldb ] = 1.;
  B[ ioffB + 0 + 1 * ldb ] = 0.;
  B[ ioffB + 0 + 2 * ldb ] = -1.;
  B[ ioffB + 1 + 0 * ldb ] = 0.;
  B[ ioffB + 1 + 1 * ldb ] = -1.;
  B[ ioffB + 1 + 2 * ldb ] = -2.;
  B[ ioffB + 2 + 0 * ldb ] = -1.;
  B[ ioffB + 2 + 1 * ldb ] = -2.;
  B[ ioffB + 2 + 2 * ldb ] = -3.;
  B[ ioffB + 3 + 0 * ldb ] = -2.;
  B[ ioffB + 3 + 1 * ldb ] = -3.;
  B[ ioffB + 3 + 2 * ldb ] = -4.;
  C[ ioffC + 0 + 0 * ldc ] = 4.;
  C[ ioffC + 1 + 0 * ldc ] = 2.;
  C[ ioffC + 0 + 1 * ldc ] = 6.;
  C[ ioffC + 1 + 1 * ldc ] = 2.;
  C[ ioffC + 0 + 2 * ldc ] = 8.;
  C[ ioffC + 1 + 2 * ldc ] = 0.;
  C[ ioffC + 0 + 3 * ldc ] = 10.;
  C[ ioffC + 1 + 3 * ldc ] = -6.;
  Blas3.dgemm('T','T',m,n,k,0.,A,lda,B,ldb,1.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dgemm('T','T',m,n,k,0.,A,lda,B,ldb,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + "\n";
  C[ ioffC + 0 + 0 * ldc ] = 4.;
  C[ ioffC + 1 + 0 * ldc ] = 2.;
  C[ ioffC + 0 + 1 * ldc ] = 6.;
  C[ ioffC + 1 + 1 * ldc ] = 2.;
  C[ ioffC + 0 + 2 * ldc ] = 8.;
  C[ ioffC + 1 + 2 * ldc ] = 0.;
  C[ ioffC + 0 + 3 * ldc ] = 10.;
  C[ ioffC + 1 + 3 * ldc ] = -6.;
  Blas3.dgemm('T','T',m,n,k,0.,A,lda,B,ldb,0.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dgemm('T','T',m,n,k,0.,A,lda,B,ldb,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + "\n";
  C[ ioffC + 0 + 0 * ldc ] = 4.;
  C[ ioffC + 1 + 0 * ldc ] = 2.;
  C[ ioffC + 0 + 1 * ldc ] = 6.;
  C[ ioffC + 1 + 1 * ldc ] = 2.;
  C[ ioffC + 0 + 2 * ldc ] = 8.;
  C[ ioffC + 1 + 2 * ldc ] = 0.;
  C[ ioffC + 0 + 3 * ldc ] = 10.;
  C[ ioffC + 1 + 3 * ldc ] = -6.;
  Blas3.dgemm('T','T',m,n,k,0.,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dgemm('T','T',m,n,k,0.,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + "\n";
  C[ ioffC + 0 + 0 * ldc ] = 4.;
  C[ ioffC + 1 + 0 * ldc ] = 2.;
  C[ ioffC + 0 + 1 * ldc ] = 6.;
  C[ ioffC + 1 + 1 * ldc ] = 2.;
  C[ ioffC + 0 + 2 * ldc ] = 8.;
  C[ ioffC + 1 + 2 * ldc ] = 0.;
  C[ ioffC + 0 + 3 * ldc ] = 10.;
  C[ ioffC + 1 + 3 * ldc ] = -6.;
  Blas3.dgemm('T','T',m,n,k,3.,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dgemm('T','T',m,n,k,3.,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 0 + 3 * ldc ] + " "
    + C[ ioffC + 1 + 3 * ldc ] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zgemm() {
  document.getElementById("debug_textarea").value +=
    "testing zgemm *************" + "\n";
  var ioffA = 1;
  var ioffB = 2;
  var ioffC = 3;
  var A = new Array( ioffA + 8 );
  var B = new Array( ioffB + 20 );
  var C = new Array( ioffC + 16 );
  var lda = 4;
  var ldb = 5;
  var ldc = 4;
  var zero = new Complex(0.,0.);
  var one = new Complex(1.,0.);
  var alpha = new Complex(3.,4.);
  var beta = new Complex(2.,5.);

//    C = A alpha B + C beta, A 3x2, B 2x4, C 3x4
  var k = 2;
  var m = 3;
  var n = 4;
  var i = -1;
  var j = -1;
  for ( j = 0; j < 2; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 4; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex();
    }
  }
  for ( j = 0; j < 4; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * ldc ] = new Complex();
    }
  }
  A[ ioffA + 0 + 0 * lda ].setValue( new Complex(2.,0.) );
  A[ ioffA + 1 + 0 * lda ].setValue( new Complex(1.,1.) );
  A[ ioffA + 2 + 0 * lda ].setValue( new Complex(0.,2.) );
  A[ ioffA + 0 + 1 * lda ].setValue( new Complex(3.,1.) );
  A[ ioffA + 1 + 1 * lda ].setValue( new Complex(2.,2.) );
  A[ ioffA + 2 + 1 * lda ].setValue( new Complex(1.,3.) );
  B[ ioffB + 0 + 0 * ldb ].setValue( new Complex(1.,-1.) );
  B[ ioffB + 1 + 0 * ldb ].setValue( new Complex(0.,0.) );
  B[ ioffB + 0 + 1 * ldb ].setValue( new Complex(0.,-2.) );
  B[ ioffB + 1 + 1 * ldb ].setValue( new Complex(-1.,-1.) );
  B[ ioffB + 0 + 2 * ldb ].setValue( new Complex(-1.,-3.) );
  B[ ioffB + 1 + 2 * ldb ].setValue( new Complex(-2.,-2.) );
  B[ ioffB + 0 + 3 * ldb ].setValue( new Complex(-2.,-4.) );
  B[ ioffB + 1 + 3 * ldb ].setValue( new Complex(-3.,-3.) );
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 0 * ldc ].setValue( new Complex(0.,4.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 1 * ldc ].setValue( new Complex(-2.,6.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  2 + 2 * ldc ].setValue( new Complex(-8.,8.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22. ));
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  C[ ioffC +  2 + 3 * ldc ].setValue( new Complex(-22.,10.) );
  Blas3.zgemm('N','N',m,n,k,zero,A,lda,B,ldb,one,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('N','N',m,n,k,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  2 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  2 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  2 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + " "
    + C[ ioffC +  2 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 0 * ldc ].setValue( new Complex(0.,4.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 1 * ldc ].setValue( new Complex(-2.,6.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  2 + 2 * ldc ].setValue( new Complex(-8.,8.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22. ));
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  C[ ioffC +  2 + 3 * ldc ].setValue( new Complex(-22.,10.) );
  Blas3.zgemm('N','N',m,n,k,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('N','N',m,n,k,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  2 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  2 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  2 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + " "
    + C[ ioffC +  2 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 0 * ldc ].setValue( new Complex(0.,4.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 1 * ldc ].setValue( new Complex(-2.,6.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  2 + 2 * ldc ].setValue( new Complex(-8.,8.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22. ));
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  C[ ioffC +  2 + 3 * ldc ].setValue( new Complex(-22.,10.) );
  Blas3.zgemm('N','N',m,n,k,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('N','N',m,n,k,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  2 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  2 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  2 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + " "
    + C[ ioffC +  2 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 0 * ldc ].setValue( new Complex(0.,4.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 1 * ldc ].setValue( new Complex(-2.,6.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  2 + 2 * ldc ].setValue( new Complex(-8.,8.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22. ));
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  C[ ioffC +  2 + 3 * ldc ].setValue( new Complex(-22.,10.) );
  Blas3.zgemm('N','N',m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('N','N',m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  2 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  2 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  2 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + " "
    + C[ ioffC +  2 + 3 * ldc ].toString() + "\n";

//    C = A' alpha B + C beta, A 3x2, B 3x4, C 2x4
  k = 3;
  m = 2;
  n = 4;
  for ( j = 0; j < 2; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 4; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 4; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC +  i + j * lda ] = new Complex();
    }
  }
  A[ ioffA + 0 + 0 * lda ].setValue( new Complex(2.,0.) );
  A[ ioffA + 1 + 0 * lda ].setValue( new Complex(1.,1.) );
  A[ ioffA + 2 + 0 * lda ].setValue( new Complex(0.,2.) );
  A[ ioffA + 0 + 1 * lda ].setValue( new Complex(3.,1.) );
  A[ ioffA + 1 + 1 * lda ].setValue( new Complex(2.,2.) );
  A[ ioffA + 2 + 1 * lda ].setValue( new Complex(1.,3.) );
  B[ ioffB + 0 + 0 * ldb ].setValue( new Complex(1.,-1.) );
  B[ ioffB + 1 + 0 * ldb ].setValue( new Complex(0.,0.) );
  B[ ioffB + 2 + 0 * ldb ].setValue( new Complex(-1.,1.) );
  B[ ioffB + 0 + 1 * ldb ].setValue( new Complex(0.,-2.) );
  B[ ioffB + 1 + 1 * ldb ].setValue( new Complex(-1.,-1.) );
  B[ ioffB + 2 + 1 * ldb ].setValue( new Complex(-2.,0.) );
  B[ ioffB + 0 + 2 * ldb ].setValue( new Complex(-1.,-3.) );
  B[ ioffB + 1 + 2 * ldb ].setValue( new Complex(-2.,-2.) );
  B[ ioffB + 2 + 2 * ldb ].setValue( new Complex(-3.,-1.) );
  B[ ioffB + 0 + 3 * ldb ].setValue( new Complex(-2.,-4.) );
  B[ ioffB + 1 + 3 * ldb ].setValue( new Complex(-3.,-3.) );
  B[ ioffB + 2 + 3 * ldb ].setValue( new Complex(-4.,-2.) );
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('T','N',m,n,k,zero,A,lda,B,ldb,one,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('T','N',m,n,k,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('T','N',m,n,k,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('T','N',m,n,k,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('T','N',m,n,k,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('T','N',m,n,k,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('T','N',m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('T','N',m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('C','N',m,n,k,zero,A,lda,B,ldb,one,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('C','N',m,n,k,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('C','N',m,n,k,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('C','N',m,n,k,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('C','N',m,n,k,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('C','N',m,n,k,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('C','N',m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('C','N',m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";

//    C = A alpha B' + C beta, A 3x2, B 4x2, C 3x4
  k = 2;
  m = 3;
  n = 4;
  for ( j = 0; j < 2; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 4; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 4; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC +  i + j * lda ] = new Complex();
    }
  }
  A[ ioffA + 0 + 0 * lda ].setValue( new Complex(2.,0.) );
  A[ ioffA + 1 + 0 * lda ].setValue( new Complex(1.,1.) );
  A[ ioffA + 2 + 0 * lda ].setValue( new Complex(0.,2.) );
  A[ ioffA + 0 + 1 * lda ].setValue( new Complex(3.,1.) );
  A[ ioffA + 1 + 1 * lda ].setValue( new Complex(2.,2.) );
  A[ ioffA + 2 + 1 * lda ].setValue( new Complex(1.,3.) );
  B[ ioffB + 0 + 0 * ldb ].setValue( new Complex(1.,-1.) );
  B[ ioffB + 0 + 1 * ldb ].setValue( new Complex(0.,0.) );
  B[ ioffB + 1 + 0 * ldb ].setValue( new Complex(0.,-2.) );
  B[ ioffB + 1 + 1 * ldb ].setValue( new Complex(-1.,-1.) );
  B[ ioffB + 2 + 0 * ldb ].setValue( new Complex(-1.,-3.) );
  B[ ioffB + 2 + 1 * ldb ].setValue( new Complex(-2.,-2.) );
  B[ ioffB + 3 + 0 * ldb ].setValue( new Complex(-2.,-4.) );
  B[ ioffB + 3 + 1 * ldb ].setValue( new Complex(-3.,-3.) );
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 0 * ldc ].setValue( new Complex(0.,4.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 1 * ldc ].setValue( new Complex(-2.,6.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  2 + 2 * ldc ].setValue( new Complex(-8.,8.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22. ));
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  C[ ioffC +  2 + 3 * ldc ].setValue( new Complex(-22.,10.) );
  Blas3.zgemm('N','T',m,n,k,zero,A,lda,B,ldb,one,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('N','T',m,n,k,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  2 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  2 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  2 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + " "
    + C[ ioffC +  2 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 0 * ldc ].setValue( new Complex(0.,4.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 1 * ldc ].setValue( new Complex(-2.,6.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  2 + 2 * ldc ].setValue( new Complex(-8.,8.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22. ));
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  C[ ioffC +  2 + 3 * ldc ].setValue( new Complex(-22.,10.) );
  Blas3.zgemm('N','T',m,n,k,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('N','T',m,n,k,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  2 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  2 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  2 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + " "
    + C[ ioffC +  2 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 0 * ldc ].setValue( new Complex(0.,4.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 1 * ldc ].setValue( new Complex(-2.,6.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  2 + 2 * ldc ].setValue( new Complex(-8.,8.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22. ));
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  C[ ioffC +  2 + 3 * ldc ].setValue( new Complex(-22.,10.) );
  Blas3.zgemm('N','T',m,n,k,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('N','T',m,n,k,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  2 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  2 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  2 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + " "
    + C[ ioffC +  2 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 0 * ldc ].setValue( new Complex(0.,4.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 1 * ldc ].setValue( new Complex(-2.,6.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  2 + 2 * ldc ].setValue( new Complex(-8.,8.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22. ));
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  C[ ioffC +  2 + 3 * ldc ].setValue( new Complex(-22.,10.) );
  Blas3.zgemm('N','T',m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('N','T',m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  2 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  2 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  2 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + " "
    + C[ ioffC +  2 + 3 * ldc ].toString() + "\n";

//    C = A' alpha B' + C beta, A 3x2, B 4x3, C 2x4
  k = 3;
  m = 2;
  n = 4;
  for ( j = 0; j < 2; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 4; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 4; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC +  i + j * lda ] = new Complex();
    }
  }
  A[ ioffA + 0 + 0 * lda ].setValue( new Complex(2.,0.) );
  A[ ioffA + 1 + 0 * lda ].setValue( new Complex(1.,1.) );
  A[ ioffA + 2 + 0 * lda ].setValue( new Complex(0.,2.) );
  A[ ioffA + 0 + 1 * lda ].setValue( new Complex(3.,1.) );
  A[ ioffA + 1 + 1 * lda ].setValue( new Complex(2.,2.) );
  A[ ioffA + 2 + 1 * lda ].setValue( new Complex(1.,3.) );
  B[ ioffB + 0 + 0 * ldb ].setValue( new Complex(1.,-1.) );
  B[ ioffB + 0 + 1 * ldb ].setValue( new Complex(0.,0.) );
  B[ ioffB + 0 + 2 * ldb ].setValue( new Complex(-1.,1.) );
  B[ ioffB + 1 + 0 * ldb ].setValue( new Complex(0.,-2.) );
  B[ ioffB + 1 + 1 * ldb ].setValue( new Complex(-1.,-1.) );
  B[ ioffB + 1 + 2 * ldb ].setValue( new Complex(-2.,0.) );
  B[ ioffB + 2 + 0 * ldb ].setValue( new Complex(-1.,-3.) );
  B[ ioffB + 2 + 1 * ldb ].setValue( new Complex(-2.,-2.) );
  B[ ioffB + 2 + 2 * ldb ].setValue( new Complex(-3.,-1.) );
  B[ ioffB + 3 + 0 * ldb ].setValue( new Complex(-2.,-4.) );
  B[ ioffB + 3 + 1 * ldb ].setValue( new Complex(-3.,-3.) );
  B[ ioffB + 3 + 2 * ldb ].setValue( new Complex(-4.,-2.) );
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('T','T',m,n,k,zero,A,lda,B,ldb,one,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('T','T',m,n,k,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('T','T',m,n,k,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('T','T',m,n,k,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('T','T',m,n,k,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('T','T',m,n,k,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('T','T',m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('T','T',m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('C','T',m,n,k,zero,A,lda,B,ldb,one,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('C','T',m,n,k,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('C','T',m,n,k,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('C','T',m,n,k,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('C','T',m,n,k,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('C','T',m,n,k,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('C','T',m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('C','T',m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";

//    C = A alpha B' + C beta, A 3x2, B 4x2, C 3x4
  k = 2;
  m = 3;
  n = 4;
  for ( j = 0; j < 2; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 4; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 4; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC +  i + j * lda ] = new Complex();
    }
  }
  A[ ioffA + 0 + 0 * lda ].setValue( new Complex(2.,0.) );
  A[ ioffA + 1 + 0 * lda ].setValue( new Complex(1.,1.) );
  A[ ioffA + 2 + 0 * lda ].setValue( new Complex(0.,2.) );
  A[ ioffA + 0 + 1 * lda ].setValue( new Complex(3.,1.) );
  A[ ioffA + 1 + 1 * lda ].setValue( new Complex(2.,2.) );
  A[ ioffA + 2 + 1 * lda ].setValue( new Complex(1.,3.) );
  B[ ioffB + 0 + 0 * ldb ].setValue( new Complex(1.,-1.) );
  B[ ioffB + 0 + 1 * ldb ].setValue( new Complex(0.,0.) );
  B[ ioffB + 1 + 0 * ldb ].setValue( new Complex(0.,-2.) );
  B[ ioffB + 1 + 1 * ldb ].setValue( new Complex(-1.,-1.) );
  B[ ioffB + 2 + 0 * ldb ].setValue( new Complex(-1.,-3.) );
  B[ ioffB + 2 + 1 * ldb ].setValue( new Complex(-2.,-2.) );
  B[ ioffB + 3 + 0 * ldb ].setValue( new Complex(-2.,-4.) );
  B[ ioffB + 3 + 1 * ldb ].setValue( new Complex(-3.,-3.) );
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 0 * ldc ].setValue( new Complex(0.,4.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 1 * ldc ].setValue( new Complex(-2.,6.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  2 + 2 * ldc ].setValue( new Complex(-8.,8.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22. ));
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  C[ ioffC +  2 + 3 * ldc ].setValue( new Complex(-22.,10.) );
  Blas3.zgemm('N','C',m,n,k,zero,A,lda,B,ldb,one,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('N','C',m,n,k,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  2 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  2 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  2 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + " "
    + C[ ioffC +  2 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 0 * ldc ].setValue( new Complex(0.,4.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 1 * ldc ].setValue( new Complex(-2.,6.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  2 + 2 * ldc ].setValue( new Complex(-8.,8.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22. ));
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  C[ ioffC +  2 + 3 * ldc ].setValue( new Complex(-22.,10.) );
  Blas3.zgemm('N','C',m,n,k,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('N','C',m,n,k,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  2 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  2 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  2 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + " "
    + C[ ioffC +  2 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 0 * ldc ].setValue( new Complex(0.,4.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 1 * ldc ].setValue( new Complex(-2.,6.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  2 + 2 * ldc ].setValue( new Complex(-8.,8.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22. ));
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  C[ ioffC +  2 + 3 * ldc ].setValue( new Complex(-22.,10.) );
  Blas3.zgemm('N','C',m,n,k,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('N','C',m,n,k,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  2 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  2 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  2 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + " "
    + C[ ioffC +  2 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 0 * ldc ].setValue( new Complex(0.,4.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  2 + 1 * ldc ].setValue( new Complex(-2.,6.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  2 + 2 * ldc ].setValue( new Complex(-8.,8.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22. ));
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  C[ ioffC +  2 + 3 * ldc ].setValue( new Complex(-22.,10.) );
  Blas3.zgemm('N','C',m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('N','C',m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  2 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  2 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  2 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + " "
    + C[ ioffC +  2 + 3 * ldc ].toString() + "\n";

//    C = A' alpha B' + C beta, A 3x2, B 4x3, C 2x4
  k = 3;
  m = 2;
  n = 4;
  for ( j = 0; j < 2; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 4; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 4; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC +  i + j * lda ] = new Complex();
    }
  }
  A[ ioffA + 0 + 0 * lda ].setValue( new Complex(2.,0.) );
  A[ ioffA + 1 + 0 * lda ].setValue( new Complex(1.,1.) );
  A[ ioffA + 2 + 0 * lda ].setValue( new Complex(0.,2.) );
  A[ ioffA + 0 + 1 * lda ].setValue( new Complex(3.,1.) );
  A[ ioffA + 1 + 1 * lda ].setValue( new Complex(2.,2.) );
  A[ ioffA + 2 + 1 * lda ].setValue( new Complex(1.,3.) );
  B[ ioffB + 0 + 0 * ldb ].setValue( new Complex(1.,-1.) );
  B[ ioffB + 0 + 1 * ldb ].setValue( new Complex(0.,0.) );
  B[ ioffB + 0 + 2 * ldb ].setValue( new Complex(-1.,1.) );
  B[ ioffB + 1 + 0 * ldb ].setValue( new Complex(0.,-2.) );
  B[ ioffB + 1 + 1 * ldb ].setValue( new Complex(-1.,-1.) );
  B[ ioffB + 1 + 2 * ldb ].setValue( new Complex(-2.,0.) );
  B[ ioffB + 2 + 0 * ldb ].setValue( new Complex(-1.,-3.) );
  B[ ioffB + 2 + 1 * ldb ].setValue( new Complex(-2.,-2.) );
  B[ ioffB + 2 + 2 * ldb ].setValue( new Complex(-3.,-1.) );
  B[ ioffB + 3 + 0 * ldb ].setValue( new Complex(-2.,-4.) );
  B[ ioffB + 3 + 1 * ldb ].setValue( new Complex(-3.,-3.) );
  B[ ioffB + 3 + 2 * ldb ].setValue( new Complex(-4.,-2.) );
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('T','C',m,n,k,zero,A,lda,B,ldb,one,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('T','C',m,n,k,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('T','C',m,n,k,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('T','C',m,n,k,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.));
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('T','C',m,n,k,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('T','C',m,n,k,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.));
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('T','C',m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('T','C',m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('C','C',m,n,k,zero,A,lda,B,ldb,one,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('C','C',m,n,k,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('C','C',m,n,k,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('C','C',m,n,k,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('C','C',m,n,k,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('C','C',m,n,k,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
  C[ ioffC +  0 + 0 * ldc ].setValue( new Complex(4.,0.) );
  C[ ioffC +  1 + 0 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 1 * ldc ].setValue( new Complex(6.,-2.) );
  C[ ioffC +  1 + 1 * ldc ].setValue( new Complex(2.,2.) );
  C[ ioffC +  0 + 2 * ldc ].setValue( new Complex(8.,-8.) );
  C[ ioffC +  1 + 2 * ldc ].setValue( new Complex(0.,0.) );
  C[ ioffC +  0 + 3 * ldc ].setValue( new Complex(10.,-22.) );
  C[ ioffC +  1 + 3 * ldc ].setValue( new Complex(-6.,-6.) );
  Blas3.zgemm('C','C',m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zgemm('C','C',m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC +  0 + 0 * ldc ].toString() + " "
    + C[ ioffC +  1 + 0 * ldc ].toString() + " "
    + C[ ioffC +  0 + 1 * ldc ].toString() + " "
    + C[ ioffC +  1 + 1 * ldc ].toString() + " "
    + C[ ioffC +  0 + 2 * ldc ].toString() + " "
    + C[ ioffC +  1 + 2 * ldc ].toString() + " "
    + C[ ioffC +  0 + 3 * ldc ].toString() + " "
    + C[ ioffC +  1 + 3 * ldc ].toString() + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dsymm() {
  document.getElementById("debug_textarea").value +=
    "testing dsymm *************" + "\n";
  var ioffA = 1;
  var ioffB = 2;
  var ioffC = 3;
  var A = new Array( ioffA + 12 );
  var B = new Array( ioffB + 12 );
  var C = new Array( ioffC + 12 );
  var lda = 4;
  var ldb = 4;
  var ldc = 4;

//    C = A alpha B + C beta, A 3x3, B 3x2, C 3x2
  var m = 3;
  var n = 2;
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
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < m; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffA + i + j * lda ] = Number( i + j * 2 );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsymm('L','L',m,n,0.,A,lda,B,ldb,1.,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsymm('L','L',m,n,0.,A,lda,B,ldb,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsymm('L','L',m,n,0.,A,lda,B,ldb,0.,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsymm('L','L',m,n,0.,A,lda,B,ldb,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsymm('L','L',m,n,0.,A,lda,B,ldb,2.,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsymm('L','L',m,n,0.,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsymm('L','L',m,n,3.,A,lda,B,ldb,2.,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsymm('L','L',m,n,3.,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsymm('L','U',m,n,0.,A,lda,B,ldb,1.,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsymm('L','U',m,n,0.,A,lda,B,ldb,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsymm('L','U',m,n,0.,A,lda,B,ldb,0.,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsymm('L','U',m,n,0.,A,lda,B,ldb,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsymm('L','U',m,n,0.,A,lda,B,ldb,2.,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsymm('L','U',m,n,0.,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsymm('L','U',m,n,3.,A,lda,B,ldb,2.,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsymm('L','U',m,n,3.,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + "\n";

//    C = B alpha A + C beta, A 3x3, B 2x3, C 2x3
  m = 2;
  n = 3;
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
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = Number( i + j * 2 );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsymm('R','L',m,n,0.,A,lda,B,ldb,1.,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsymm('R','L',m,n,0.,A,lda,B,ldb,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsymm('R','L',m,n,0.,A,lda,B,ldb,0.,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsymm('R','L',m,n,0.,A,lda,B,ldb,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsymm('R','L',m,n,0.,A,lda,B,ldb,2.,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsymm('R','L',m,n,0.,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsymm('R','L',m,n,3.,A,lda,B,ldb,2.,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsymm('R','L',m,n,3.,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsymm('R','U',m,n,0.,A,lda,B,ldb,1.,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsymm('R','U',m,n,0.,A,lda,B,ldb,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsymm('R','U',m,n,0.,A,lda,B,ldb,0.,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsymm('R','U',m,n,0.,A,lda,B,ldb,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsymm('R','U',m,n,0.,A,lda,B,ldb,2.,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsymm('R','U',m,n,0.,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsymm('R','U',m,n,3.,A,lda,B,ldb,2.,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsymm('R','U',m,n,3.,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zsymm() {
  document.getElementById("debug_textarea").value +=
    "testing zsymm *************" + "\n";
  var ioffA = 1;
  var ioffB = 2;
  var ioffC = 3;
  var A = new Array( ioffA + 12 );
  var B = new Array( ioffB + 12 );
  var C = new Array( ioffC + 12 );
  var lda = 4;
  var ldb = 4;
  var ldc = 4;
  var zero = new Complex( 0. , 0. );
  var one = new Complex( 1. , 0. );
  var alpha = new Complex( 3. , 4. );
  var beta = new Complex( 2. , 5. );

//    C = A alpha B + C beta, A 3x3, B 3x2, C 3x2
  var m = 3;
  var n = 2;
  var i = -1;
  var j = -1;
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * ldc ] = new Complex();
    }
  }
  for ( j = 0; j < m; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j * 2 , 2 * i + j );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsymm('L','L',m,n,zero,A,lda,B,ldb,one,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsymm('L','L',m,n,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsymm('L','L',m,n,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsymm('L','L',m,n,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsymm('L','L',m,n,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsymm('L','L',m,n,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsymm('L','L',m,n,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsymm('L','L',m,n,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsymm('L','U',m,n,zero,A,lda,B,ldb,one,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsymm('L','U',m,n,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsymm('L','U',m,n,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsymm('L','U',m,n,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsymm('L','U',m,n,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsymm('L','U',m,n,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsymm('L','U',m,n,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsymm('L','U',m,n,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + "\n";

//    C = B alpha A + C beta, A 3x3, B 2x3, C 2x3
  m = 2;
  n = 3;
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j * 2 , 2 * i + j );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsymm('R','L',m,n,zero,A,lda,B,ldb,one,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsymm('R','L',m,n,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsymm('R','L',m,n,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsymm('R','L',m,n,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsymm('R','L',m,n,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsymm('R','L',m,n,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsymm('R','L',m,n,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsymm('R','L',m,n,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsymm('R','U',m,n,zero,A,lda,B,ldb,one,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsymm('R','U',m,n,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsymm('R','U',m,n,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsymm('R','U',m,n,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsymm('R','U',m,n,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsymm('R','U',m,n,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsymm('R','U',m,n,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsymm('R','U',m,n,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zhemm() {
  document.getElementById("debug_textarea").value +=
    "testing zhemm *************" + "\n";
  var ioffA = 1;
  var ioffB = 2;
  var ioffC = 3;
  var A = new Array( ioffA + 12 );
  var B = new Array( ioffB + 12 );
  var C = new Array( ioffC + 12 );
  var lda = 4;
  var ldb = 4;
  var ldc = 4;
  var zero = new Complex( 0. , 0. );
  var one = new Complex( 1. , 0. );
  var alpha = new Complex( 3. , 4. );
  var beta = new Complex( 2. , 5. );

//    C = A alpha B + C beta, A 3x3, B 3x2, C 3x2
  var m = 3;
  var n = 2;
  var i = -1;
  var j = -1;
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * ldc ] = new Complex();
    }
  }
  for ( j = 0; j < m; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffA + i + j * lda ] =
        new Complex( i + j * 2 , 2 * i + j );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zhemm('L','L',m,n,zero,A,lda,B,ldb,one,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zhemm('L','L',m,n,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zhemm('L','L',m,n,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
     "zhemm('L','L',m,n,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zhemm('L','L',m,n,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zhemm('L','L',m,n,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zhemm('L','L',m,n,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zhemm('L','L',m,n,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zhemm('L','U',m,n,zero,A,lda,B,ldb,one,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zhemm('L','U',m,n,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zhemm('L','U',m,n,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zhemm('L','U',m,n,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zhemm('L','U',m,n,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zhemm('L','U',m,n,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zhemm('L','U',m,n,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zhemm('L','U',m,n,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + "\n";

//    C = B alpha A + C beta, A 3x3, B 2x3, C 2x3
  m = 2;
  n = 3;
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j * 2 , 2 * i + j );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zhemm('R','L',m,n,zero,A,lda,B,ldb,one,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zhemm('R','L',m,n,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zhemm('R','L',m,n,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zhemm('R','L',m,n,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zhemm('R','L',m,n,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zhemm('R','L',m,n,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zhemm('R','L',m,n,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zhemm('R','L',m,n,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zhemm('R','U',m,n,zero,A,lda,B,ldb,one,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zhemm('R','U',m,n,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zhemm('R','U',m,n,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zhemm('R','U',m,n,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zhemm('R','U',m,n,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zhemm('R','U',m,n,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zhemm('R','U',m,n,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA, ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zhemm('R','U',m,n,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dsyrk() {
  document.getElementById("debug_textarea").value +=
    "testing dsyrk *************" + "\n";
  var ioffA = 1;
  var ioffC = 2;
  var A = new Array( ioffA + 12 );
  var C = new Array( ioffC + 12 );
  var lda = 4;
  var ldc = 4;

//    C = A alpha A' + C beta, A 3x2, C 3x3
  var k = 2;
  var n = 3;
  var i = -1;
  var j = -1;
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < k; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = Number( i + j * 2 );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyrk('L','N',n,k,0.,A,lda,1.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyrk('L','N',n,k,0.,A,lda,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyrk('L','N',n,k,0.,A,lda,0.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyrk('L','N',n,k,0.,A,lda,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyrk('L','N',n,k,0.,A,lda,2.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyrk('L','N',n,k,0.,A,lda,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyrk('L','N',n,k,3.,A,lda,2.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyrk('L','N',n,k,3.,A,lda,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyrk('U','N',n,k,0.,A,lda,1.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyrk('U','N',n,k,0.,A,lda,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyrk('U','N',n,k,0.,A,lda,0.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyrk('U','N',n,k,0.,A,lda,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyrk('U','N',n,k,0.,A,lda,2.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyrk('U','N',n,k,0.,A,lda,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyrk('U','N',n,k,3.,A,lda,2.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyrk('U','N',n,k,3.,A,lda,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";

//    C = A' alpha A + C beta, A 2x3, C 3x3
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < k; i ++ ) {
      A[ ioffA + i + j * lda ] = Number( i + j * 2 );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyrk('L','T',n,k,0.,A,lda,1.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyrk('L','T',n,k,0.,A,lda,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyrk('L','T',n,k,0.,A,lda,0.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyrk('L','T',n,k,0.,A,lda,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyrk('L','T',n,k,0.,A,lda,2.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyrk('L','T',n,k,0.,A,lda,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyrk('L','T',n,k,3.,A,lda,2.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyrk('L','T',n,k,3.,A,lda,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyrk('U','T',n,k,0.,A,lda,1.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyrk('U','T',n,k,0.,A,lda,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyrk('U','T',n,k,0.,A,lda,0.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyrk('U','T',n,k,0.,A,lda,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyrk('U','T',n,k,0.,A,lda,2.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyrk('U','T',n,k,0.,A,lda,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyrk('U','T',n,k,3.,A,lda,2.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyrk('U','T',n,k,3.,A,lda,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zsyrk() {
  document.getElementById("debug_textarea").value +=
    "testing zsyrk *************" + "\n";
  var ioffA = 1;
  var ioffC = 2;
  var A = new Array( ioffA + 12 );
  var C = new Array( ioffC + 12 );
  var lda = 4;
  var ldc = 4;
  var zero = new Complex( 0., 0. );
  var one = new Complex( 1., 0. );
  var alpha = new Complex( 3., 4. );
  var beta = new Complex( 2., 5. );

//    C = A alpha A' + C beta, A 3x2, C 3x3
  var k = 2;
  var n = 3;
  var i = -1;
  var j = -1;
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * ldc ] = new Complex();
    }
  }
  for ( j = 0; j < k; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] =
        new Complex( i + j * 2 , 2 * i + j );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyrk('L','N',n,k,zero,A,lda,one,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyrk('L','N',n,k,zero,A,lda,one,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyrk('L','N',n,k,zero,A,lda,zero,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyrk('L','N',n,k,zero,A,lda,zero,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyrk('L','N',n,k,zero,A,lda,beta,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyrk('L','N',n,k,zero,A,lda,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyrk('L','N',n,k,alpha,A,lda,beta,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyrk('L','N',n,k,alpha,A,lda,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyrk('U','N',n,k,zero,A,lda,one,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyrk('U','N',n,k,zero,A,lda,one,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyrk('U','N',n,k,zero,A,lda,zero,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyrk('U','N',n,k,zero,A,lda,zero,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] 
        = new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyrk('U','N',n,k,zero,A,lda,beta,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyrk('U','N',n,k,zero,A,lda,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyrk('U','N',n,k,alpha,A,lda,beta,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyrk('U','N',n,k,alpha,A,lda,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";

//    C = A' alpha A + C beta, A 2x3, C 3x3
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < k; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j * 2 , 2 * i + j );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyrk('L','T',n,k,zero,A,lda,one,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyrk('L','T',n,k,zero,A,lda,one,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyrk('L','T',n,k,zero,A,lda,zero,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyrk('L','T',n,k,zero,A,lda,zero,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyrk('L','T',n,k,zero,A,lda,beta,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyrk('L','T',n,k,zero,A,lda,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyrk('L','T',n,k,alpha,A,lda,beta,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyrk('L','T',n,k,alpha,A,lda,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyrk('U','T',n,k,zero,A,lda,one,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyrk('U','T',n,k,zero,A,lda,one,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyrk('U','T',n,k,zero,A,lda,zero,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyrk('U','T',n,k,zero,A,lda,zero,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyrk('U','T',n,k,zero,A,lda,beta,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyrk('U','T',n,k,zero,A,lda,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyrk('U','T',n,k,alpha,A,lda,beta,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyrk('U','T',n,k,alpha,A,lda,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zherk() {
  document.getElementById("debug_textarea").value +=
    "testing zherk *************" + "\n";
  var ioffA = 1;
  var ioffC = 2;
  var A = new Array( ioffA + 12 );
  var C = new Array( ioffC + 12 );
  var lda = 4;
  var ldc = 4;

//    C = A alpha A' + C beta, A 3x2, C 3x3
  var k = 2;
  var n = 3;
  var i = -1;
  var j = -1;
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * ldc ] = new Complex();
    }
  }
  for ( j = 0; j < k; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j * 2 , 2 * i + j );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zherk('L','N',n,k,0.,A,lda,1.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zherk('L','N',n,k,0.,A,lda,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zherk('L','N',n,k,0.,A,lda,0.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zherk('L','N',n,k,0.,A,lda,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zherk('L','N',n,k,0.,A,lda,2.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zherk('L','N',n,k,0.,A,lda,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zherk('L','N',n,k,3.,A,lda,2.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zherk('L','N',n,k,3.,A,lda,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zherk('U','N',n,k,0.,A,lda,1.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zherk('U','N',n,k,0.,A,lda,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zherk('U','N',n,k,0.,A,lda,0.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zherk('U','N',n,k,0.,A,lda,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zherk('U','N',n,k,0.,A,lda,2.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zherk('U','N',n,k,0.,A,lda,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zherk('U','N',n,k,3.,A,lda,2.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zherk('U','N',n,k,3.,A,lda,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";

//    C = A' alpha A + C beta, A 2x3, C 3x3
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < k; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j * 2 , 2 * i + j );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zherk('L','C',n,k,0.,A,lda,1.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zherk('L','C',n,k,0.,A,lda,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zherk('L','C',n,k,0.,A,lda,0.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zherk('L','C',n,k,0.,A,lda,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zherk('L','C',n,k,0.,A,lda,2.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zherk('L','C',n,k,0.,A,lda,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zherk('L','C',n,k,3.,A,lda,2.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zherk('L','C',n,k,3.,A,lda,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zherk('U','C',n,k,0.,A,lda,1.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zherk('U','C',n,k,0.,A,lda,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zherk('U','C',n,k,0.,A,lda,0.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
     "zherk('U','C',n,k,0.,A,lda,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zherk('U','C',n,k,0.,A,lda,2.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zherk('U','C',n,k,0.,A,lda,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zherk('U','C',n,k,3.,A,lda,2.,C,ldc,ioffA,ioffC);
  document.getElementById("debug_textarea").value +=
    "zherk('U','C',n,k,3.,A,lda,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dsyr2k() {
  document.getElementById("debug_textarea").value +=
    "testing dsyr2k *************" + "\n";
  var ioffA = 1;
  var ioffB = 2;
  var ioffC = 3;
  var A = new Array( ioffA + 12 );
  var B = new Array( ioffB + 12 );
  var C = new Array( ioffC + 12 );
  var lda = 4;
  var ldb = 4;
  var ldc = 4;

//    C = A alpha B' + B alpha A' + C beta, A 3x2, B 3x2, C 3x3
  var n = 3;
  var k = 2;
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
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < k; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = Number( i + j * 2 );
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyr2k('L','N',n,k,0.,A,lda,B,ldb,1.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyr2k('L','N',n,k,0.,A,lda,B,ldb,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyr2k('L','N',n,k,0.,A,lda,B,ldb,0.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyr2k('L','N',n,k,0.,A,lda,B,ldb,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyr2k('L','N',n,k,0.,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyr2k('L','N',n,k,0.,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyr2k('L','N',n,k,3.,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyr2k('L','N',n,k,3.,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyr2k('U','N',n,k,0.,A,lda,B,ldb,1.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyr2k('U','N',n,k,0.,A,lda,B,ldb,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyr2k('U','N',n,k,0.,A,lda,B,ldb,0.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyr2k('U','N',n,k,0.,A,lda,B,ldb,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyr2k('U','N',n,k,0.,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyr2k('U','N',n,k,0.,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyr2k('U','N',n,k,3.,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyr2k('U','N',n,k,3.,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";

//    C = A' alpha B + B' alpha A + C beta, A 2x3, B 2x3, C 3x3
  k = 2;
  n = 3;
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
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number.POSITIVE_INFINITY;
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < k; i ++ ) {
      A[ ioffA + i + j * lda ] = Number( i + j * 2 );
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyr2k('L','T',n,k,0.,A,lda,B,ldb,1.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyr2k('L','T',n,k,0.,A,lda,B,ldb,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyr2k('L','T',n,k,0.,A,lda,B,ldb,0.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyr2k('L','T',n,k,0.,A,lda,B,ldb,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyr2k('L','T',n,k,0.,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyr2k('L','T',n,k,0.,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyr2k('L','T',n,k,3.,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyr2k('L','T',n,k,3.,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 0 * ldc ] + " "
    + C[ ioffC + 2 + 0 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 1 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyr2k('U','T',n,k,0.,A,lda,B,ldb,1.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyr2k('U','T',n,k,0.,A,lda,B,ldb,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyr2k('U','T',n,k,0.,A,lda,B,ldb,0.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyr2k('U','T',n,k,0.,A,lda,B,ldb,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyr2k('U','T',n,k,0.,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyr2k('U','T',n,k,0.,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] = Number( 2 * i - 3 * j );
    }
  }
  Blas3.dsyr2k('U','T',n,k,3.,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "dsyr2k('U','T',n,k,3.,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ] + " "
    + C[ ioffC + 0 + 1 * ldc ] + " "
    + C[ ioffC + 1 + 1 * ldc ] + " "
    + C[ ioffC + 0 + 2 * ldc ] + " "
    + C[ ioffC + 1 + 2 * ldc ] + " "
    + C[ ioffC + 2 + 2 * ldc ] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zsyr2k() {
  document.getElementById("debug_textarea").value +=
    "testing zsyr2k *************" + "\n";
  var ioffA = 1;
  var ioffB = 2;
  var ioffC = 3;
  var A = new Array( ioffA + 12 );
  var B = new Array( ioffB + 12 );
  var C = new Array( ioffC + 12 );
  var lda = 4;
  var ldb = 4;
  var ldc = 4;
  var zero = new Complex( 0. , 0. );
  var one = new Complex( 1. , 0. );
  var alpha = new Complex( 3. , 4. );
  var beta = new Complex( 2. , 5. );

//    C = A alpha B' + B alpha A' + C beta, A 3x2, B 3x2, C 3x3
  var n = 3;
  var k = 2;
  var i = -1;
  var j = -1;
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * ldc ] = new Complex();
    }
  }
  for ( j = 0; j < k; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j * 2 , 2 * i + j );
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyr2k('L','N',n,k,zero,A,lda,B,ldb,one,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
     "zsyr2k('L','N',n,k,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyr2k('L','N',n,k,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyr2k('L','N',n,k,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyr2k('L','N',n,k,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyr2k('L','N',n,k,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyr2k('L','N',n,k,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyr2k('L','N',n,k,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyr2k('U','N',n,k,zero,A,lda,B,ldb,one,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyr2k('U','N',n,k,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyr2k('U','N',n,k,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyr2k('U','N',n,k,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyr2k('U','N',n,k,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyr2k('U','N',n,k,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyr2k('U','N',n,k,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyr2k('U','N',n,k,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";

//    C = A' alpha B + B' alpha A + C beta, A 2x3, B 2x3, C 3x3
  k = 2;
  n = 3;
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < k; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j * 2 , 2 * i + j );
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyr2k('L','T',n,k,zero,A,lda,B,ldb,one,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyr2k('L','T',n,k,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyr2k('L','T',n,k,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyr2k('L','T',n,k,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyr2k('L','T',n,k,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyr2k('L','T',n,k,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyr2k('L','T',n,k,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyr2k('L','T',n,k,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyr2k('U','T',n,k,zero,A,lda,B,ldb,one,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyr2k('U','T',n,k,zero,A,lda,B,ldb,one,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyr2k('U','T',n,k,zero,A,lda,B,ldb,zero,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyr2k('U','T',n,k,zero,A,lda,B,ldb,zero,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyr2k('U','T',n,k,zero,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyr2k('U','T',n,k,zero,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zsyr2k('U','T',n,k,alpha,A,lda,B,ldb,beta,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zsyr2k('U','T',n,k,alpha,A,lda,B,ldb,beta,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zher2k() {
  document.getElementById("debug_textarea").value +=
    "testing zher2k *************" + "\n";
  var ioffA = 1;
  var ioffB = 2;
  var ioffC = 3;
  var A = new Array( ioffA + 12 );
  var B = new Array( ioffB + 12 );
  var C = new Array( ioffC + 12 );
  var lda = 4;
  var ldb = 4;
  var ldc = 4;
  var zero = new Complex( 0. , 0. );
  var alpha = new Complex( 3. , 4. );

//    C = A alpha B' + B alpha A' + C beta, A 3x2, B 3x2, C 3x3
  var n = 3;
  var k = 2;
  var i = -1;
  var j = -1;
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * ldc ] = new Complex();
    }
  }
  for ( j = 0; j < k; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j * 2 , 2 * i + j );
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zher2k('L','N',n,k,zero,A,lda,B,ldb,1.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zher2k('L','N',n,k,zero,A,lda,B,ldb,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zher2k('L','N',n,k,zero,A,lda,B,ldb,0.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zher2k('L','N',n,k,zero,A,lda,B,ldb,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zher2k('L','N',n,k,zero,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zher2k('L','N',n,k,zero,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zher2k('L','N',n,k,alpha,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zher2k('L','N',n,k,alpha,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zher2k('U','N',n,k,zero,A,lda,B,ldb,1.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zher2k('U','N',n,k,zero,A,lda,B,ldb,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zher2k('U','N',n,k,zero,A,lda,B,ldb,0.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zher2k('U','N',n,k,zero,A,lda,B,ldb,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zher2k('U','N',n,k,zero,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zher2k('U','N',n,k,zero,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zher2k('U','N',n,k,alpha,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zher2k('U','N',n,k,alpha,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";

//    C = A' alpha B + B' alpha A + C beta, A 2x3, B 2x3, C 3x3
  k = 2;
  n = 3;
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldc; i ++ ) {
      C[ ioffC + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < k; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j * 2 , 2 * i + j );
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zher2k('L','C',n,k,zero,A,lda,B,ldb,1.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zher2k('L','C',n,k,zero,A,lda,B,ldb,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zher2k('L','C',n,k,zero,A,lda,B,ldb,0.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zher2k('L','C',n,k,zero,A,lda,B,ldb,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zher2k('L','C',n,k,zero,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zher2k('L','C',n,k,zero,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zher2k('L','C',n,k,alpha,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zher2k('L','C',n,k,alpha,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 0 * ldc ].toString() + " "
    + C[ ioffC + 2 + 0 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 1 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";

  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zher2k('U','C',n,k,zero,A,lda,B,ldb,1.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zher2k('U','C',n,k,zero,A,lda,B,ldb,1.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zher2k('U','C',n,k,zero,A,lda,B,ldb,0.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zher2k('U','C',n,k,zero,A,lda,B,ldb,0.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zher2k('U','C',n,k,zero,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zher2k('U','C',n,k,zero,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      C[ ioffC + i + j * ldc ] =
        new Complex( 2 * i - 3 * j , 3 * i - 2 * j );
    }
  }
  Blas3.zher2k('U','C',n,k,alpha,A,lda,B,ldb,2.,C,ldc,
    ioffA,ioffB,ioffC);
  document.getElementById("debug_textarea").value +=
    "zher2k('U','C',n,k,alpha,A,lda,B,ldb,2.,C,ldc) , C = "
    + C[ ioffC + 0 + 0 * ldc ].toString() + " "
    + C[ ioffC + 0 + 1 * ldc ].toString() + " "
    + C[ ioffC + 1 + 1 * ldc ].toString() + " "
    + C[ ioffC + 0 + 2 * ldc ].toString() + " "
    + C[ ioffC + 1 + 2 * ldc ].toString() + " "
    + C[ ioffC + 2 + 2 * ldc ].toString() + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dtrmm() {
  document.getElementById("debug_textarea").value +=
    "testing dtrmm *************" + "\n";
  var ioffA = 1;
  var ioffB = 2;
  var A = new Array( ioffA + 12 );
  var B = new Array( ioffB + 12 );
  var lda = 4;
  var ldb = 4;

//    B = A alpha B, A 3x3, B 3x2
  var m = 3;
  var n = 2;
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
  for ( j = 0; j < m; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffA + i + j * lda ] = Number( i + j * 2 );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrmm('L','L','N','U',m,n,0.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
     "dtrmm('L','L','N','U',m,n,0.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrmm('L','L','N','U',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrmm('L','L','N','U',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrmm('L','L','N','N',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrmm('L','L','N','N',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrmm('L','U','N','U',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrmm('L','U','N','U',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrmm('L','U','N','N',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrmm('L','U','N','N',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";

//    B = A' alpha B, A 3x3, B 3x2
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrmm('L','L','T','U',m,n,0.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrmm('L','L','T','U',m,n,0.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrmm('L','L','T','U',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrmm('L','L','T','U',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrmm('L','L','T','N',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrmm('L','L','T','N',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrmm('L','U','T','U',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrmm('L','U','T','U',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrmm('L','U','T','N',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrmm('L','U','T','N',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";

//    B = B alpha A, A 3x3, B 2x3
  m = 2;
  n = 3;
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
      A[ ioffA + i + j * lda ] = Number( i + j * 2 );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrmm('R','L','N','U',m,n,0.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrmm('R','L','N','U',m,n,0.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrmm('R','L','N','U',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrmm('R','L','N','U',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrmm('R','L','N','N',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrmm('R','L','N','N',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrmm('R','U','N','U',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrmm('R','U','N','U',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrmm('R','U','N','N',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrmm('R','U','N','N',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";

//    B = B alpha A', A 3x3, B 2x3
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrmm('R','L','T','U',m,n,0.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrmm('R','L','T','U',m,n,0.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrmm('R','L','T','U',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrmm('R','L','T','U',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrmm('R','L','T','N',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrmm('R','L','T','N',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrmm('R','U','T','U',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrmm('R','U','T','U',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrmm('R','U','T','N',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrmm('R','U','T','N',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_ztrmm() {
  document.getElementById("debug_textarea").value +=
    "testing ztrmm *************" + "\n";
  var ioffA = 1;
  var ioffB = 2;
  var A = new Array( ioffA + 12 );
  var B = new Array( ioffB + 12 );
  var lda = 4;
  var ldb = 4;
  var zero = new Complex( 0., 0. );
  var alpha = new Complex( 2., 3. );

//    B = A alpha B, A 3x3, B 3x2
  var m = 3;
  var n = 2;
  var i = -1;
  var j = -1;
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex();
    }
  }
  for ( j = 0; j < m; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j * 2 , 2 * i + j );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('L','L','N','U',m,n,zero,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('L','L','N','U',m,n,zero,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('L','L','N','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('L','L','N','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('L','L','N','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('L','L','N','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('L','U','N','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('L','U','N','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('L','U','N','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('L','U','N','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";

//    B = A' alpha B, A 3x3, B 3x2
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('L','L','T','U',m,n,zero,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('L','L','T','U',m,n,zero,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('L','L','T','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('L','L','T','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('L','L','T','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('L','L','T','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('L','U','T','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('L','U','T','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('L','U','T','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('L','U','T','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('L','L','C','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('L','L','C','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('L','L','C','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('L','L','C','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('L','U','C','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('L','U','C','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('L','U','C','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('L','U','C','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";

//    B = B alpha A, A 3x3, B 2x3
  m = 2;
  n = 3;
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex();
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex( i + j * 2 , 2 * i + j );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('R','L','N','U',m,n,zero,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('R','L','N','U',m,n,zero,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('R','L','N','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('R','L','N','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('R','L','N','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('R','L','N','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('R','U','N','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('R','U','N','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('R','U','N','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('R','U','N','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";

//    B = B alpha A', A 3x3, B 2x3
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('R','L','T','U',m,n,zero,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('R','L','T','U',m,n,zero,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('R','L','T','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('R','L','T','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('R','L','T','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('R','L','T','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('R','U','T','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('R','U','T','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('R','U','T','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('R','U','T','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('R','L','C','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('R','L','C','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('R','L','C','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('R','L','C','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('R','U','C','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('R','U','C','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrmm('R','U','C','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrmm('R','U','C','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dtrsm() {
  document.getElementById("debug_textarea").value +=
    "testing dtrsm *************" + "\n";
  var ioffA = 1;
  var ioffB = 2;
  var A = new Array( ioffA + 12 );
  var B = new Array( ioffB + 12 );
  var lda = 4;
  var ldb = 4;

//    B = A alpha B, A 3x3, B 3x2
  var m = 3;
  var n = 2;
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
  for ( j = 0; j < m; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffA + i + j * lda ] = Number( 1 + i + j * 2 );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrsm('L','L','N','U',m,n,0.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrsm('L','L','N','U',m,n,0.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrsm('L','L','N','U',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrsm('L','L','N','U',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrsm('L','L','N','N',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrsm('L','L','N','N',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrsm('L','U','N','U',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrsm('L','U','N','U',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrsm('L','U','N','N',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrsm('L','U','N','N',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";

//    B = A' alpha B, A 3x3, B 3x2
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrsm('L','L','T','U',m,n,0.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrsm('L','L','T','U',m,n,0.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrsm('L','L','T','U',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrsm('L','L','T','U',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrsm('L','L','T','N',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrsm('L','L','T','N',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrsm('L','U','T','U',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrsm('L','U','T','U',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrsm('L','U','T','N',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrsm('L','U','T','N',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";

//    B = B alpha A, A 3x3, B 2x3
  m = 2;
  n = 3;
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
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrsm('R','L','N','U',m,n,0.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrsm('R','L','N','U',m,n,0.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrsm('R','L','N','U',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrsm('R','L','N','U',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrsm('R','L','N','N',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrsm('R','L','N','N',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrsm('R','U','N','U',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrsm('R','U','N','U',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrsm('R','U','N','N',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrsm('R','U','N','N',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";

//    B = B alpha A', A 3x3, B 2x3
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrsm('R','L','T','U',m,n,0.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrsm('R','L','T','U',m,n,0.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrsm('R','L','T','U',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrsm('R','L','T','U',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrsm('R','L','T','N',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrsm('R','L','T','N',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrsm('R','U','T','U',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrsm('R','U','T','U',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = Number( 3 * i + j );
    }
  }
  Blas3.dtrsm('R','U','T','N',m,n,2.,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "dtrsm('R','U','T','N',m,n,2.,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_ztrsm() {
  document.getElementById("debug_textarea").value +=
    "testing ztrsm *************" + "\n";
  var ioffA = 1;
  var ioffB = 2;
  var A = new Array( ioffA + 12 );
  var B = new Array( ioffB + 12 );
  var lda = 4;
  var ldb = 4;
  var zero = new Complex( 0., 0. );
  var alpha = new Complex( 2., 3. );

//    B = A alpha B, A 3x3, B 3x2
  var m = 3;
  var n = 2;
  var i = -1;
  var j = -1;
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex();
    }
  }
  for ( j = 0; j < m; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      A[ ioffA + i + j * lda ] =
        new Complex( 1 + i + j * 2 , 2 * i + j );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('L','L','N','U',m,n,zero,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('L','L','N','U',m,n,zero,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('L','L','N','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('L','L','N','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('L','L','N','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('L','L','N','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('L','U','N','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('L','U','N','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('L','U','N','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('L','U','N','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";

//    B = A' alpha B, A 3x3, B 3x2
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('L','L','T','U',m,n,zero,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('L','L','T','U',m,n,zero,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('L','L','T','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('L','L','T','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('L','L','T','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('L','L','T','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('L','U','T','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('L','U','T','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('L','U','T','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('L','U','T','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('L','L','C','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('L','L','C','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('L','L','C','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('L','L','C','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('L','U','C','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('L','U','C','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('L','U','C','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('L','U','C','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 2 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 2 + 1 * ldb ] + "\n";

//    B = B alpha A, A 3x3, B 2x3
  m = 2;
  n = 3;
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < lda; i ++ ) {
      A[ ioffA + i + j * lda ] = new Complex();
    }
  }
  for ( j = 0; j < 3; j ++ ) {
    for ( i = 0; i < ldb; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex();
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < n; i ++ ) {
      A[ ioffA + i + j * lda ] =
        new Complex( 1 + i + j * 2 , 2 * i + j );
    }
  }
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('R','L','N','U',m,n,zero,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('R','L','N','U',m,n,zero,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('R','L','N','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('R','L','N','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('R','L','N','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('R','L','N','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('R','U','N','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('R','U','N','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('R','U','N','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('R','U','N','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";

//    B = B alpha A', A 3x3, B 2x3
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('R','L','T','U',m,n,zero,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('R','L','T','U',m,n,zero,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('R','L','T','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('R','L','T','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('R','L','T','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('R','L','T','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('R','U','T','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('R','U','T','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('R','U','T','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('R','U','T','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('R','L','C','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('R','L','C','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('R','L','C','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('R','L','C','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('R','U','C','U',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('R','U','C','U',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
  for ( j = 0; j < n; j ++ ) {
    for ( i = 0; i < m; i ++ ) {
      B[ ioffB + i + j * ldb ] = new Complex( 3 * i + j , i + 3 * j );
    }
  }
  Blas3.ztrsm('R','U','C','N',m,n,alpha,A,lda,B,ldb,ioffA,ioffB);
  document.getElementById("debug_textarea").value +=
    "ztrsm('R','U','C','N',m,n,alpha,A,lda,B,ldb) , B = "
    + B[ ioffB + 0 + 0 * ldb ] + " "
    + B[ ioffB + 1 + 0 * ldb ] + " "
    + B[ ioffB + 0 + 1 * ldb ] + " "
    + B[ ioffB + 1 + 1 * ldb ] + " "
    + B[ ioffB + 0 + 2 * ldb ] + " "
    + B[ ioffB + 1 + 2 * ldb ] + "\n";
}
