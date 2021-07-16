function StringReference( s ) {
  this.the_string = s;
}
StringReference.prototype.toString = function() {
  return this.the_string;
}
StringReference.prototype.setValue = function( s ) { this.the_string = s; } 
StringReference.prototype.getValue = function( ) { return this.the_string; } 
