function dividedDifference( n, x, y, difs ) {
  for ( var j_dd = 0; j_dd <= n; j_dd ++ ) difs[ j_dd ] = y[ j_dd ];
  for ( var j_dd = 1; j_dd <= n; j_dd ++ ) {
    for ( var i_dd = n; i_dd >= j_dd; i_dd -- ) {
      difs[ i_dd ] = ( difs[ i_dd ] - difs[ i_dd - 1 ] )
                   / ( x[ i_dd ] - x[ i_dd - j_dd ] );
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
function newtonPolynomial( n, difs, x, t ) {
  var np = difs[ n ];
  for ( var j_np = n - 1; j_np >= 0; j_np -- ) {
    np = difs[ j_np ] + np * ( t - x[ j_np ] );
  }
  return np;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
function addPoint( n, x0, y0, difs, x ) {
  for ( var j = 1; j <= n + 1; j ++ ) x[ j ] = x[ j - 1 ];
  x[ 0 ] = x0;
  var temp1 = difs[ 0 ];
  difs[ 0 ] = y0;
  for ( var j = 1; j <= n + 1; j ++ ) {
    var temp2 = temp1;
    temp1 = difs[ j ];
    difs[ j ] = ( temp2 - difs[ j - 1 ] )
                    / ( x[ j ] - x0 );
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
function lagrangePolynomial( n, x, y, t ) {
  var lp = 0;
  for ( var j = 0; j <= n; j ++ ) {
    var prod = y[ j ];
    for ( var i = 0; i < j; i ++ ) {
      prod *= ( t - x[ i ] ) / ( x[ j ] - x[ i ] );
    }
    for ( var i = j + 1; i <= n; i ++ ) {
      prod *= ( t - x[ i ] ) / ( x[ j ] - x[ i ] );
    }
    lp += prod;
  }
  return lp;
}
