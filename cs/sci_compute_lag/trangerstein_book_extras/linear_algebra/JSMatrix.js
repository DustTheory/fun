function JSMatrix( m, n, scalar ) {
//document.getElementById("debug_textarea").value +=
//  "in JSMatrix, m,n = " + m + " " + n + "\n";
  if ( m == undefined || m < 0 ) m=0;
  if ( n == undefined || n < 0 ) n=0;
  this.nrows=m;
  this.ncols=n;
  this.data = new Array( m * n );
  if ( scalar == undefined ) scalar = Math.POSITIVE_INFINITY;
  this.fillWith(scalar);
}
JSMatrix.prototype.fillWith = function ( s ) {
  for ( var i = 0; i < this.data.length; i++ ) this.data[i] = s; 
}
JSMatrix.prototype.toString = function() {
  var s = new String();
  if (this.data.length > 0 ) {
    for ( var i = 0; i < this.nrows - 1; i++ ) {
      for ( var j = 0; j < this.ncols - 1; j ++ ) {
        s += this.data[ i + j * this.nrows ].toString() + " , ";
      }
      s += this.data[ i + (this.ncols - 1 ) * this.nrows ].toString()
        + " ; ";
    }
    for ( j = 0; j < this.ncols - 1; j ++ ) {
      s += this.data[ this.nrows - 1 + j * this.nrows ].toString() + " , ";
    }
    s += this.data[ this.nrows - 1 + (this.ncols - 1 ) * this.nrows
      ].toString() + " ; ";
  }
  return s;
}
JSMatrix.prototype.parse = function( t ) {
//document.getElementById("debug_textarea").value +=
//  "in JSMatrix.parse, t = " + t + "\n";
  this.nrows = 0;
  this.ncols = 0;
  var a = t.split( /;/ );
  for (var i = 0; i < a.length; i++) {
    var b = a[i].split( /,/ );
//  document.getElementById("debug_textarea").value +=
//    "a[ " + i + " ] = " + b + "\n";
    this.ncols = Math.max( this.ncols, b.length );
    this.nrows ++;
  }
//document.getElementById("debug_textarea").value +=
//  "nrows,ncols = " + this.nrows + " " + this.ncols + "\n";
  var nr = 0;
  this.data = new Array( this.nrows * this.ncols );
  this.fillWith(0);
  a = t.split( /;/ );
  for ( i = 0; i < a.length; i++) {
    b = a[i].split( /,/ );
    for ( var j = 0; j < b.length; j ++ ) {
      this.data[ nr + j * this.nrows ] = Number( b[j] );
    }
    nr ++;
  }
//for ( i = 0; i < this.nrows; i ++ ) {
//  document.getElementById("debug_textarea").value +=
//    "matrix[ " + i + " , * ] = ";
//  for ( j = 0; j <= i; j ++ ) {
//    document.getElementById("debug_textarea").value +=
//      this.getEntry( i, j ) + " ";
//  }
//  document.getElementById("debug_textarea").value += "\n";
//}
}
JSMatrix.prototype.numberRows = function() { return this.nrows; }
JSMatrix.prototype.numberCols = function() { return this.ncols; }
JSMatrix.prototype.dataArray = function() { return this.data; }
JSMatrix.prototype.length = function() { return this.data.length; }
JSMatrix.prototype.getEntry = function( i, j ) { 
  return this.data[ i + j * this.nrows ];
}
JSMatrix.prototype.setEntry = function( i, j, v ) {
  this.data[ i + j * this.nrows ] = v;
}
JSMatrix.prototype.copy = function( c ) {
  this.data.length = c.data.length;
  this.nrows = c.nrows;
  this.ncols = c.ncols;
  for ( var i = 0; i < this.data.length; i++ ) this.data[i] = c.data[i]; 
}
JSMatrix.prototype.resize = function( m, n ) {
  if (! ( m * n == this.data.length ) ) {
    this.data.length = m * n;
    this.nrows = m;
    this.ncols = n;
    this.fillWith(Number.POSITIVE_INFINITY);
  }
}
JSMatrix.prototype.plusEquals = function( c ) {
  for ( var i = 0; i < this.data.length; i++ ) this.data[i] += c.data[i];
}
JSMatrix.prototype.minusEquals = function( c ) {
  for ( var i = 0; i < this.data.length; i++ ) this.data[i] -= c.data[i];
}
JSMatrix.prototype.timesEquals = function( s ) {
  for ( var i = 0; i < this.data.length; i++ ) this.data[i] *= s;
}
JSMatrix.prototype.divideEquals = function( s ) {
  for ( var i = 0; i < this.data.length; i++ ) this.data[i] /= s;
}
JSMatrix.prototype.plus = function( c ) {
  var sum = new JSMatrix();
  sum.copy( this );
  sum.plusEquals( c );
  return sum;
}
JSMatrix.prototype.minus = function( c ) {
  var dif = new JSMatrix();
  dif.copy( this );
  dif.minusEquals( c );
  return dif;
}
JSMatrix.prototype.times = function( s ) {
  var prod = new JSMatrix();
  prod.copy( this );
  prod.timesEquals( s );
  return prod;
}
JSMatrix.prototype.divide = function( s ) {
  var quot = new JSMatrix();
  quot.copy( this );
  quot.divideEquals( s );
  return quot;
}
JSMatrix.prototype.norm = function( p ) {
  var nrm = 0;
  var i;
  if ( p <= 1. ) {
    var max_col_sum=0;
    for ( j = 0; j < this.ncols; j++ ) {
      var col_sum=0;
      for ( i = 0; i < this.nrows; i++ ) {
        col_sum += Math.abs( data[i+j*nrows] );
      }
      max_col_sum=Math.max(max_col_sum,col_sum);
    }
    return max_col_sum;
  } else { // infinity norm
    var max_row_sum=0;
    for ( i = 0; i < this.nrows; i++ ) {
      var row_sum=0;
      for ( j = 0; j < this.ncols; j++ ) {
        row_sum += Math.abs( data[i+j*nrows] );
      }
      max_row_sum=Math.max(max_row_sum,row_sum);
    }
    return max_row_sum;
  }
}
//JSMatrix.prototype.debug = function() {
//  trace("  length = " + this.data.length);
//  trace("  nrows,ncols = " + this.nrows + " " + this.ncols);
//  for (var i = 0; i < this.data.length; i++) {
//    trace(" data[" + i + "] = " + this.data[i]);
//  }
//}
