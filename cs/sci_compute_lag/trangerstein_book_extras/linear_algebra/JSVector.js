function JSVector( n, scalar ) {
  if (n<0) n=0;
  this.data = new Array(n);
  if ( scalar == undefined ) scalar = Math.POSITIVE_INFINITY;
  this.fillWith(scalar);
}
JSVector.prototype.fillWith = function ( s ) {
  for ( var i = 0; i < this.data.length; i++ ) this.data[i] = s; 
}
JSVector.prototype.toString = function() {
  var s = new String();
  if (this.data.length > 0 ) {
    for ( var i = 0; i < this.data.length-1; i++ ) {
      s += this.data[i].toString() + " ; ";
    }
    s += this.data[this.data.length-1].toString();
  }
  return s;
}
JSVector.prototype.parse = function ( s ) {
  var a = s.split( /;/ );
  for (var i = 0; i < a.length; i++) this.data[i] = Number( a[i] );
}
JSVector.prototype.dataArray = function() { return this.data; }
JSVector.prototype.length = function() { return this.data.length; }
JSVector.prototype.getEntry = function( i ) { return this.data[i]; }
JSVector.prototype.setEntry = function( i, v ) { this.data[i] = v; }
JSVector.prototype.copy = function( c ) {
  this.data.length = c.data.length;
  for ( var i = 0; i < this.data.length; i++ ) {
    this.data[i] = c.data[i]; 
  }
}
JSVector.prototype.resize = function( n ) {
  if (! ( n == this.data.length ) ) {
    this.data.length=n;
    this.fillWith(Number.POSITIVE_INFINITY);
  }
}
JSVector.prototype.plusEquals = function( c ) {
  for ( var i = 0; i < this.data.length; i++ ) {
    this.data[i] += c.data[i];
  }
}
JSVector.prototype.minusEquals = function( c ) {
  for ( var i = 0; i < this.data.length; i++ ) {
    this.data[i] -= c.data[i];
  }
}
JSVector.prototype.timesEquals = function( s ) {
  for ( var i = 0; i < this.data.length; i++ ) this.data[i] *= s;
}
JSVector.prototype.divideEquals = function( s ) {
  for ( var i = 0; i < this.data.length; i++ ) this.data[i] /= s;
}
JSVector.prototype.plus = function( c ) {
  var sum = new JSVector();
  sum.copy( this );
  sum.plusEquals( c );
  return sum;
}
JSVector.prototype.minus = function( c ) {
  var dif = new JSVector();
  dif.copy( this );
  dif.minusEquals( c );
  return dif;
}
JSVector.prototype.times = function( s ) {
  var prod = new JSVector();
  prod.copy( this );
  prod.timesEquals( s );
  return prod;
}
JSVector.prototype.divide = function( s ) {
  var quot = new JSVector();
  quot.copy( this );
  quot.divideEquals( s );
  return quot;
}
JSVector.prototype.norm = function( p ) {
  var nrm = 0;
  var i;
  if ( p <= 1 ) {
    for ( i = 0; i < this.data.length; i++ ) {
      nrm += Math.abs( this.data[i] );
    }
    return nrm;
  } else {
    var m = 0;
    for ( i = 0; i < this.data.length; i++ ) {
      m=Math.max(m,Math.abs(this.data[i]));
    }
    if ( p >= Number.POSITIVE_INFINITY || m <= 0 ) return m;
    for ( i = 0; i < this.data.length; i++ ) {
      nrm += Math.pow( Math.abs( this.data[i]/m ) , p );
    }
    return m*Math.pow( nrm , 1 / p );
  }
}
//JSVector.prototype.debug = function() { 
//  trace("  length = " + this.data.length);
//  for (var i:int = 0; i < data.length; i++) {
//    trace(" this.data[" + i + "] = " + this.data[i]);
//  }
//}
