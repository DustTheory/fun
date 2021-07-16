function Colormap( p ) {
//document.getElementById("debug_textarea").value +=
//  "entering Colormap constructor\n";
  if ( ! p.isOk() ) alert( "palette is not OK" );
  this.palette = p;
  this.num_cmap_colors = p.map_index[ p.num_named_colors - 1 ] + 1;
//document.getElementById("debug_textarea").value +=
//  "num_cmap_colors = " + this.num_cmap_colors + "\n";
  this.cmap_colors = new Array( this.num_cmap_colors );
  for ( var cm_i = 0; cm_i < p.num_named_colors; cm_i ++ ) {
    this.cmap_colors[ p.map_index[ cm_i ] ] = p.color_value[ cm_i ]; 
//  document.getElementById("debug_textarea").value +=
//    "cmap_colors[" + p.map_index[ cm_i ] + "] = "
//    + this.cmap_colors[ p.map_index[ cm_i ] ] + "\n";
  }
//interpolate between Palette colors:
  for ( var cm_i = 3; cm_i < p.num_named_colors; cm_i ++ ) {
    var cm_map_lo = this.palette.map_index[ cm_i - 1 ]; 
    var cm_map_hi = this.palette.map_index[ cm_i ]; 
    var cm_denom = cm_map_hi - cm_map_lo;
    var cm_low = [ this.cmap_colors[ cm_map_lo ][ 0 ],
                this.cmap_colors[ cm_map_lo ][ 1 ],
                this.cmap_colors[ cm_map_lo ][ 2 ] ];
    var cm_high = [ this.cmap_colors[ cm_map_hi ][ 0 ],
                 this.cmap_colors[ cm_map_hi ][ 1 ], 
                 this.cmap_colors[ cm_map_hi ][ 2 ] ]; 
    for ( var cm_j = cm_map_lo + 1; cm_j < cm_map_hi; cm_j ++ ) {
      var cm_lo_ratio = ( cm_map_hi - cm_j ) / cm_denom;
      var cm_hi_ratio = ( cm_j - cm_map_lo ) / cm_denom;
      this.cmap_colors[ cm_j ] = new Array( 3 );
      for ( var cm_k = 0; cm_k < 3; cm_k ++ ) {
        this.cmap_colors[ cm_j ][ cm_k ] =
          Math.round( cm_low[ cm_k ] * cm_lo_ratio
          + cm_high[ cm_k ] * cm_hi_ratio );
      }
    }
  }
//for ( var cm_i = 0; cm_i < this.num_cmap_colors; cm_i ++ ) {
//  document.getElementById("debug_textarea").value +=
//    "cmap_colors[" + cm_i + "] = " + this.cmap_colors[ cm_i ] + "\n";
//}
//document.getElementById("debug_textarea").value +=
//  "leaving Colormap constructor\n";
}
//add-on colors are not interpolated
Colormap.prototype.allocColor = function( rgb ) { 
  var cm_color_index = this.palette.findEntry( rgb );
  if ( cm_color_index < 0 ) {
    this.palette.insertEntry( this.num_cmap_colors, rgb );
    this.cmap_colors[ this.num_cmap_colors ] = 
      palette.color_value[ this.num_cmap_colors ];
    return this.cmap_colors[ this.num_cmap_colors ++ ];
  } else return cmap_colors[ cm_color_index ];
}
Colormap.prototype.setColor = function( i, rgb ) {
  if ( i >= 0 && i < this.num_cmap_colors ) {
    this.cmap_colors[ i ][ 0 ] = rgb[ 0 ];
    this.cmap_colors[ i ][ 1 ] = rgb[ 1 ];
    this.cmap_colors[ i ][ 2 ] = rgb[ 2 ];
  }
}
Colormap.prototype.getHexColor = function( i ) {
//document.getElementById("debug_textarea").value +=
//  "entering Colormap.getHexColor\n";
  if ( i < 0 || i >= this.num_cmap_colors ) alert("invalid colormap index");
//document.getElementById("debug_textarea").value +=
//  "i = " + i + "\n";
  var cm_hex_string = new String();
  var cm_numerals = "0123456789ABCDEF";
  var cm_ci = this.cmap_colors[ i ];
  cm_hex_string = "#";
  for ( var cm_j = 0; cm_j < 3; cm_j ++ ) {
    cm_hex_string +=
      cm_numerals[ Math.floor( cm_ci[ cm_j ] / 16 ) ].toString()
      + cm_numerals[ cm_ci[ cm_j ] % 16 ].toString();
  }
//document.getElementById("debug_textarea").value +=
//  "leaving Colormap.getHexColor\n";
  return cm_hex_string;
}
Colormap.prototype.getGLColor = function( i ) {
  if ( i < 0 || i >= this.num_cmap_colors ) alert("invalid colormap index");
  var cm_ci = this.cmap_colors[ i ];
  for ( var cm_j = 0; cm_j < 3; cm_j ++ ) cm_ci[ cm_j ] /= 255;
  return cm_ci;
}
Colormap.prototype.getGLColorFromRatio = function( r ) {
  r = Math.max( 0, Math.min( 1, r ) );
  var cm_i = Math.round( 2 + r * ( this.num_cmap_colors - 3 ) ); 
  var cm_ci = new Array( 3 );
  for ( var cm_j = 0; cm_j < 3; cm_j ++ ) cm_ci[ cm_j ]
    = this.cmap_colors[ cm_i ][ cm_j ] / 255;
  return cm_ci;
}
