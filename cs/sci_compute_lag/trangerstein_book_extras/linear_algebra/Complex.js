function Complex( r, c ) {
  if ( r == undefined ) r = Math.POSITIVE_INFINITY;
  if ( c == undefined ) c = Math.POSITIVE_INFINITY;
  this.real_part = r;
  this.imaginary_part = c;
}
Complex.prototype.toString = function() {
  var s = new String();
  s = "< " + this.real_part + " , " + this.imaginary_part + " >";
  return s;
}
Complex.prototype.getReal = function() { return this.real_part; }
Complex.prototype.getImag = function() { return this.imaginary_part; }
Complex.prototype.setValue = function( z ) {
  this.real_part = z.real_part;
  this.imaginary_part = z.imaginary_part;
} 
Complex.prototype.pureReal = function() { this.imaginary_part = 0; } 
Complex.prototype.pureImag = function() { this.real_part = 0; } 
Complex.prototype.equals = function( z ) {
  return this.real_part == z.real_part &&
         this.imaginary_part == z.imaginary_part;
} 
Complex.prototype.plusEquals = function( z ) {
  this.real_part += z.real_part;
  this.imaginary_part += z.imaginary_part;
  return this;
} 
Complex.prototype.minus = function() {
  return new Complex( - this.real_part, - this.imaginary_part );
} 
Complex.prototype.minusEquals = function( z ) {
  this.real_part -= z.real_part;
  this.imaginary_part -= z.imaginary_part;
  return this;
} 
Complex.prototype.timesNumberEquals = function( x ) {
  this.real_part *= x;
  this.imaginary_part *= x
  return this;
} 
Complex.prototype.timesEquals = function( z ) {
  var r = this.real_part * z.real_part
        - this.imaginary_part * z.imaginary_part;
  this.imaginary_part = this.imaginary_part * z.real_part
                      + this.real_part * z.imaginary_part;
  this.real_part = r;
  return this;
} 
Complex.prototype.divideNumberEquals = function( x ) {
  this.real_part /= x;
  this.imaginary_part /= x
  return this;
} 
Complex.prototype.divideEquals = function( z ) {
  var abs_z = Math.sqrt( z.real_part * z.real_part
                       + z.imaginary_part * z.imaginary_part );
  var x = z.real_part / abs_z;
  var y = z.imaginary_part / abs_z;
  var r = ( this.real_part * x + this.imaginary_part * y ) / abs_z;
  this.imaginary_part
    = ( this.imaginary_part * x - this.real_part * y ) / abs_z;
  this.real_part = r;
  return this;
} 

function ComplexMath() {}

ComplexMath.abs = function( z ) {
  var m = Math.max( Math.abs( z.getReal() ), Math.abs( z.getImag() ) );
  if ( m > 0 ) {
    var x = z.getReal() / m;
    var y = z.getImag() / m;
    m *= Math.sqrt( x * x + y * y );
  }
  return m;
}
ComplexMath.plusNumber = function( z, r ) {
  return new Complex( z.getReal() + r , z.getImag() );
} 
ComplexMath.plus = function( z1, z2 ) {
  return new Complex( z1.getReal() + z2.getReal(), 
                      z1.getImag() + z2.getImag() );
} 
ComplexMath.minusNumber = function( z, r ) {
  return new Complex( z.getReal() - r , z.getImag() );
} 
ComplexMath.minus = function( z1, z2 ) {
  return new Complex( z1.getReal() - z2.getReal(),
                      z1.getImag() - z2.getImag() );
} 
ComplexMath.timesNumber = function( z, r ) {
  return new Complex( z.getReal() * r , z.getImag() * r );
} 
ComplexMath.times = function( z1, z2 ) {
  return new Complex( z1.getReal() * z2.getReal()
                    - z1.getImag() * z2.getImag(),
                      z1.getImag() * z2.getReal()
                    + z1.getReal() * z2.getImag() );
} 
ComplexMath.divideNumber = function( z, r ) {
  return new Complex( z.getReal() / r , z.getImag() / r );
} 
ComplexMath.divide = function( z1, z2 ) {
  var abs_z2 = ComplexMath.abs ( z2 );
  var x = z2.getReal() / abs_z2;
  var y = z2.getImag() / abs_z2;
  return new Complex( ( z1.getReal() * x + z1.getImag() * y ) / abs_z2 ,
                      ( z1.getImag() * x - z1.getReal() * y ) / abs_z2 );
} 
ComplexMath.arg = function( z ) {
  var angle = 0.;
  if ( Math.abs( z.getReal() ) > 0. ) {
    angle = Math.atan( z.getImag() / z.getReal() );
  } else if ( z.getImag() > 0. ) {
    angle = 0.5 * Math.PI;
  } else if ( z.getImag() < 0. ) {
    angle = -0.5 * Math.PI;
  }
  return angle;
}
ComplexMath.norm = function( z ) {
  var x = z.getReal();
  var y = z.getImag();
  return x * x + y * y;
}
ComplexMath.conj = function( z ) {
  return new Complex( z.getReal(), - z.getImag() );
}
ComplexMath.polar = function( radius, angle ) {
  return new Complex( radius * Math.cos( angle ), 
                      radius * Math.sin( angle ) );
}
ComplexMath.cos = function( z ) {
  var x = z.getReal();
  var y = z.getImag();
  return new Complex(
     Math.cos( x ) * ( Math.exp( y ) + Math.exp( -y ) ) * 0.5 ,
    -Math.sin( x ) * ( Math.exp( y ) - Math.exp( -y ) ) * 0.5 );
}
ComplexMath.cosh = function( z ) {
  var x = z.getReal();
  var y = z.getImag();
  return new Complex(
    ( Math.exp( x ) + Math.exp( -x ) ) * 0.5 * Math.cos( y ),
   -( Math.exp( x ) - Math.exp( -x ) ) * 0.5 * Math.sin( y ) );
}
ComplexMath.exp = function( z ) {
  return this.polar( Math.exp( z.getReal() ), z.getImag() );
}
ComplexMath.log = function( z ) {
  return new Complex( Math.log( this.abs(z) ) , this.arg(z) );
}
ComplexMath.log10 = function( z ) {
  return new Complex( Math.log( abs(z) ) * Math.LOG10E ,
                      this.arg(z) * Math.LOG10E );
}
ComplexMath.powNumber = function( z, r ) {
  var radius = this.abs( z );
  var angle = this.arg( z );
  return this.polar( Math.pow( radius, r ), angle * r );
}
ComplexMath.powInt = function( z, n ) {
  var radius = this.abs( z );
  var angle = this.arg( z );
  return this.polar( Math.pow( radius, n ), angle * Number(n) );
}
ComplexMath.pow = function( z, c ) {
  return this.exp( this.times( c , this.log( z ) ) );
}
ComplexMath.sin = function( z ) {
  var x = z.getReal();
  var y = z.getImag();
  return new Complex(
    Math.sin( x ) * ( Math.exp( y ) + Math.exp( -y ) ) * 0.5 ,
    Math.cos( x ) * ( Math.exp( y ) - Math.exp( -y ) ) * 0.5 );
}
ComplexMath.sinh = function( z ) {
  var x = z.getReal();
  var y = z.getImag();
  return new Complex(
    ( Math.exp( x ) - Math.exp( -x ) ) * 0.5 * Math.cos( y ),
    ( Math.exp( x ) + Math.exp( -x ) ) * 0.5 * Math.sin( y ) );
}
ComplexMath.sqrt = function( z ) {
  return ComplexMath.polar( Math.sqrt( ComplexMath.abs(z) ),
    0.5*ComplexMath.arg(z) );
}
ComplexMath.sqrt2 = function( z ) {
  return ComplexMath.polar( Math.sqrt( ComplexMath.abs(z) ),
    0.5*ComplexMath.arg(z)+Math.PI );
}
ComplexMath.tan = function( z ) {
  return this.divide ( this.sin( z ) , this.cos( z ) );
}
ComplexMath.tanh = function( z ) {
  return this.divide ( this.sinh( z ) , this.cosh( z ) );
}
