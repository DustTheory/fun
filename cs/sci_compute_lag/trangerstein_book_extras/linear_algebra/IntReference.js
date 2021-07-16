function IntReference( r ) {
  if ( r == undefined ) r=Number.MAX_VALUE;
  this.the_number = Math.round( r );
}
IntReference.prototype.toString = function() {
  var s = new String();
  s = this.the_number.toString();
  return s;
}
IntReference.prototype.setValue = function ( r ) { 
  this.the_number = Math.round( r );
} 
IntReference.prototype.getValue = function( ) { return this.the_number; } 
