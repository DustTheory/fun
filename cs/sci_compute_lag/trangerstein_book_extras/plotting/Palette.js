function Palette( s ) {
  this.max_colors = 256;
  this.num_named_colors = 0;
  this.map_index = new Array( this.max_colors );
  this.color_value = new Array( this.max_colors );
  if ( s == undefined || s == "default" ) {
    this.insertEntry(  0, [ 255, 255, 255 ] ); // white
    this.insertEntry(  1, [   0,   0,   0 ] ); // black
    this.insertEntry(  2, [   0,   0, 255 ] ); // blue
    this.insertEntry( 19, [   0, 255, 255 ] ); // cyan
    this.insertEntry( 39, [   0, 255,   0 ] ); // green
    this.insertEntry( 59, [ 255, 255,   0 ] ); // yellow
    this.insertEntry( 79, [ 255,   0,   0 ] ); // red
    this.name = "default";
  } else this.name = s;
}
Palette.prototype.isOk = function() {
  if ( this.num_named_colors < 3 ) return false;
  if ( this.map_index[ 0 ] != 0 ) return false;
  if ( this.map_index[ 1 ] != 1 ) return false;
  for ( var pal_i = 2; pal_i < this.num_named_colors; pal_i ++ ) {
    if ( this.map_index[ pal_i - 1 ] >= this.map_index[ pal_i ] ) {
      return false;
    }
  }
  return true;
}
Palette.prototype.insertEntry = function( index, rgb ) {
  index = Math.max( 0, Math.min( this.max_colors - 1, index ) );
  var pal_i = undefined;
  if ( this.num_named_colors == 0 ||
  this.map_index[ this.num_named_colors - 1 ] < index ) {
    pal_i = this.num_named_colors;
  } else {
    for ( pal_i = 0; pal_i < this.num_named_colors; pal_i ++ ) {
      if ( this.map_index[ pal_i ] > index ) break;
    }
    for ( var pal_j = this.num_named_colors; pal_j > pal_i; pal_j -- ) {
      this.map_index[ pal_j ] = this.map_index[ pal_j - 1 ];
      this.color_value[ pal_j ] = this.color_value[ pal_j - 1 ];
    }
  }
  this.map_index[ pal_i ] = index;
  this.color_value[ pal_i ] = rgb;
  this.num_named_colors ++;
}
Palette.prototype.findEntry = function( rgb ) {
  for ( var pal_i = 0; pal_i < this.num_named_colors; pal_i ++ ) {
    if ( rgb[ 0 ] == this.color_value[ pal_i ][ 0 ] &&
    rgb[ 1 ] == this.color_value[ pal_i ][ 1 ] &&
    rgb[ 2 ] == this.color_value[ pal_i ][ 2 ] ) {
      return pal_i;
    }
  }
  return -1;
}
