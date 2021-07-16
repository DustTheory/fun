function NumberReference( r ) {
  if ( r == undefined ) r=Number.POSITIVE_INFINITY;
  this.the_number = r;
}
NumberReference.prototype.toString = function() {
  var s = new String();
  s = this.the_number.toString();
  return s;
}
NumberReference.prototype.setValue = function ( r ) { 
  this.the_number = r;
} 
NumberReference.prototype.getValue = function( ) { return this.the_number; } 
