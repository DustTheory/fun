function BooleanReference( r ) {
  if ( r == undefined ) r = false;
  this.the_number = r;
}
BooleanReference.prototype.toString = function() {
  var s = new String();
  s = this.the_number.toString();
  return s;
}
BooleanReference.prototype.setValue = function( r ) { this.the_number = r; }
BooleanReference.prototype.getValue = function( ) { 
  return this.the_number;
} 
