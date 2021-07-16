function AreaGraphTool( c , s , xlab, ylab, low, high, cm ) {
  if ( c == undefined ) {
    alert("AreaGraphTool constructor: canvas is undefined");
  }
  if ( cm == undefined ) {
    p = new Palette();
    cm = new Colormap( p );
  }
  this.canvas = c;
  this.name = s;
  this.context = this.canvas.getContext("2d");
  this.colormap = cm;
  this.xlabel = xlab; // x axis label
  this.ylabel = ylab; // y axis label
  this.watch_mouse = false;
  this.user_low = [ undefined, undefined ]; // user coordinates
  this.user_high = [ undefined, undefined ];
  this.low = [ undefined, undefined ]; // expanded user coordinates
  this.high = [ undefined, undefined ];
  this.len = [ undefined, undefined ];
  this.mouse_coords = [ undefined, undefined ]; // user coordinates of mouse
  this.button = undefined;
  this.mouseDown = undefined; // function to be called when mouse down
  this.mouseMove = undefined; // function to be called when mouse moves
  this.mouseUp = undefined;   // function to be called when mouse up

  this.rescale( low, high );
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.functionVectorField = function( f, nc ) {
  var agt_len = [ ( this.user_high[ 0 ] - this.user_low[ 0 ] ) / nc[ 0 ],
          ( this.user_high[ 1 ] - this.user_low[ 1 ] ) / nc[ 1 ] ];
  var agt_farray = new Array( nc[ 0 ] * nc[ 1 ] );
  var agt_fmax = 0;
  for ( var agt_j = 0; agt_j < nc[ 1 ]; agt_j ++ ) {
    var agt_y = this.user_low[ 1 ] + agt_len[ 1 ] * ( agt_j + 0.5 );
    for ( var agt_i = 0; agt_i < nc[ 0 ]; agt_i ++ ) {
      var agt_x = this.user_low[ 0 ] + agt_len[ 0 ] * ( agt_i + 0.5 );
      var agt_fval = f( agt_x, agt_y );
      agt_farray[ agt_i + agt_j * nc[ 0 ] ] = agt_fval;
      agt_fmax = Math.max( agt_fmax,
        Math.sqrt( agt_fval[0]*agt_fval[0] + agt_fval[1]*agt_fval[1] ) );
    }
  }
  var agt_angle = 2.6;
  var agt_arrow = new Complex( Math.cos( agt_angle ) * .3,
    Math.sin( agt_angle ) * .3 );
  var agt_scale = 0.5 * Math.min( agt_len[0] , agt_len[1] ) / agt_fmax;
  for ( var agt_j = 0; agt_j < nc[ 1 ]; agt_j ++ ) {
    var agt_y = this.user_low[ 1 ] + agt_len[ 1 ] * ( agt_j + 0.5 );
    for ( var agt_i = 0; agt_i < nc[ 0 ]; agt_i ++ ) {
      var agt_x = this.user_low[ 0 ] + agt_len[ 0 ] * ( agt_i + 0.5 );
      var agt_base = new Complex( agt_x, agt_y );
      var agt_fval = agt_farray[ agt_i + agt_j * nc[ 0 ] ];
      var agt_head = new Complex( agt_fval[0] * agt_scale,
        agt_fval[1] * agt_scale );
      this.drawVector( agt_base, agt_head, agt_arrow );
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.colorFunction = function( f, nc ) {
  var agt_len = [ ( this.user_high[ 0 ] - this.user_low[ 0 ] ) / nc[ 0 ],
              ( this.user_high[ 1 ] - this.user_low[ 1 ] ) / nc[ 1 ] ];
  var agt_farray = new Array( ( nc[ 0 ] + 1 ) * ( nc[ 1 ] + 1 ) );
  var agt_fmin = Number.MAX_VALUE;
  var agt_fmax = - Number.MAX_VALUE;
  var agt_nc0 = nc[ 0 ] + 1;
  for ( var agt_j = 0; agt_j <= nc[ 1 ]; agt_j ++ ) {
    var agt_y = this.user_low[ 1 ] + agt_len[ 1 ] * agt_j; 
    for ( var agt_i = 0; agt_i <= nc[ 0 ]; agt_i ++ ) {
      var agt_x = this.user_low[ 0 ] + agt_len[ 0 ] * agt_i; 
      var agt_z = f( agt_x, agt_y );
      agt_farray[ agt_i + agt_j * agt_nc0 ] = agt_z;
      agt_fmin = Math.min( agt_fmin, agt_z );
      agt_fmax = Math.max( agt_fmax, agt_z );
    }
  }
  this.setBgColor();
  this.newPage();
  this.setFgColor();
  this.drawAxes();
  var agt_corners = new Array( 4 );
  for ( var agt_j = 0; agt_j < 4; agt_j ++ ) {
    agt_corners[ agt_j ] = new Array( 2 );
  }
  for ( var agt_j = 0; agt_j < nc[ 1 ]; agt_j ++ ) {
    agt_corners[ 1 ][ 1 ] = agt_corners[ 0 ][ 1 ] =
      this.user_low[ 1 ] + agt_len[ 1 ] * agt_j;
    agt_corners[ 3 ][ 1 ] = agt_corners[ 2 ][ 1 ] =
      agt_corners[ 1 ][ 1 ] + agt_len[ 1 ];
    for ( var agt_i = 0; agt_i < nc[ 0 ]; agt_i ++ ) {
      agt_corners[ 3 ][ 0 ] = agt_corners[ 0 ][ 0 ] =
        this.user_low[ 0 ] + agt_len[ 0 ] * agt_i;
      agt_corners[ 2 ][ 0 ] = agt_corners[ 1 ][ 0 ] =
        agt_corners[ 0 ][ 0 ] + agt_len[ 0 ];
      var agt_ratio = ( agt_farray[ agt_i + agt_j * agt_nc0 ] - agt_fmin )
                / ( agt_fmax - agt_fmin );
      agt_ratio = Math.max( 0, Math.min( 1, agt_ratio ) );
      this.colorPolygon( agt_corners, agt_ratio );
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.contourFunction = function( f, ncell, ncontour ) {
  var agt_len = [ ( this.user_high[ 0 ] - this.user_low[ 0 ] ) / ncell[ 0 ],
              ( this.user_high[ 1 ] - this.user_low[ 1 ] ) / ncell[ 1 ] ];
  var agt_farray = new Array( ( ncell[ 0 ] + 1 ) * ( ncell[ 1 ] + 1 ) );
  var agt_fmin = Number.MAX_VALUE;
  var agt_fmax = - Number.MAX_VALUE; 
  var agt_nc0 = ncell[ 0 ] + 1;
  for ( var agt_j = 0; agt_j <= ncell[ 1 ]; agt_j ++ ) {
    var agt_y = this.user_low[ 1 ] + agt_len[ 1 ] * agt_j; 
    for ( var agt_i = 0; agt_i <= ncell[ 0 ]; agt_i ++ ) {
      var agt_x = this.user_low[ 0 ] + agt_len[ 0 ] * agt_i; 
      var agt_z = f( agt_x, agt_y );
      agt_farray[ agt_i + agt_j * agt_nc0 ] = agt_z;
      agt_fmin = Math.min( agt_fmin, agt_z );
      agt_fmax = Math.max( agt_fmax, agt_z );
    }
  }
  this.setBgColor();
  this.newPage();
  this.setFgColor();
  this.drawAxes();
  var agt_dif = new Array( 5 );
  var agt_xpos = new Array( 5 );
  var agt_ypos = new Array( 5 );
  var agt_xy = new Array( 4 );
  for ( var agt_ic = 0; agt_ic < ncontour; agt_ic ++ ) {
    var agt_ratio = ( agt_ic + 1 ) / ( ncontour + 1 );
    var agt_value = agt_fmin + agt_ratio * ( agt_fmax - agt_fmin );
    this.setFgColor( agt_ratio );
    for ( var agt_j = 0; agt_j < ncell[ 1 ]; agt_j ++ ) {
      agt_ypos[ 4 ] = agt_ypos[ 1 ] = agt_ypos[ 0 ] =
        this.user_low[ 1 ] + agt_len[ 1 ] * agt_j;
      agt_ypos[ 3 ] = agt_ypos[ 2 ] = agt_ypos[ 1 ] + agt_len[ 1 ];
      for ( var agt_i = 0; agt_i < ncell[ 0 ]; agt_i ++ ) {
        agt_xpos[ 4 ] = agt_xpos[ 3 ] = agt_xpos[ 0 ] =
          this.user_low[ 0 ] + agt_len[ 0 ] * agt_i;
        agt_xpos[ 2 ] = agt_xpos[ 1 ] = agt_xpos[ 0 ] + agt_len[ 0 ];
        var agt_ij = agt_i + agt_j * agt_nc0;
        agt_dif[ 0 ] = agt_farray[ agt_ij ] - agt_value;
        agt_dif[ 1 ] = agt_farray[ agt_ij + 1 ] - agt_value;
        agt_dif[ 2 ] = agt_farray[ agt_ij + agt_nc0 + 1 ] - agt_value;
        agt_dif[ 3 ] = agt_farray[ agt_ij + agt_nc0 ] - agt_value;
        agt_dif[ 4 ] = agt_dif[ 0 ];
        var agt_count = 0;
        var agt_first = -1;
        var agt_last = -1;
        for ( var agt_k = 0; agt_k < 4; agt_k ++ ) {
          if ( agt_dif[ agt_k ] * agt_dif[ agt_k + 1 ] <= 0 &&
          Math.abs( agt_dif[ agt_k + 1 ] ) > 0 ) {
            if ( agt_count == 0 ) agt_first = agt_k;
            else agt_last = agt_k;
            agt_count ++;
            var agt_pos = agt_dif[ agt_k ]
              / ( agt_dif[ agt_k ] - agt_dif[ agt_k + 1 ] );
            agt_xy[ agt_k ] =
              [ agt_xpos[ agt_k ]
              + agt_pos * ( agt_xpos[ agt_k + 1 ] - agt_xpos[ agt_k ] ),
                agt_ypos[ agt_k ]
                + agt_pos * ( agt_ypos[ agt_k + 1 ] - agt_ypos[ agt_k ] ) ];
          }
        }
        if ( agt_count == 4 ) {
          this.beginDrawing();
            this.movePen( agt_agt_xy[ 0 ] );
            this.drawLine( agt_xy[ 2 ] );
          this.endDrawing();
          this.beginDrawing();
            this.movePen( agt_xy[ 1 ] );
            this.drawLine( agt_xy[ 3 ] );
          this.endDrawing();
        } else if ( agt_count == 2 ) {
          this.beginDrawing();
            this.movePen( agt_xy[ agt_first ] );
            this.drawLine( agt_xy[ agt_last ] );
          this.endDrawing();
        }
      }
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.listen = function() {
//see p 785 for mouse events
  this.canvas.graphtool = this;
  this.canvas.addEventListener( "mousedown", AreaGraphTool.mouseHandler,
    false );
  this.canvas.addEventListener( "mousemove", AreaGraphTool.mouseHandler,
    false );
  this.canvas.addEventListener( "mouseup", AreaGraphTool.mouseHandler,
    false );
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.newPage = function() {
  this.context.clearRect( 0, 0, this.canvas.width, this.canvas.height );
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.setLineWidth = function( lw ) { 
  this.context.lineWidth = lw;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.setFgColor = function( ratio ) { 
//document.getElementById("debug_textarea").value +=
//  "entering AreaGraphTool.setFgColor\n";
//document.getElementById("debug_textarea").value +=
//  "ratio = " + ratio + "\n";
  if ( ratio == undefined ) { // for text:
    var agt_hc = this.colormap.getHexColor( 1 );
    this.context.fillStyle = agt_hc;
    this.context.strokeStyle = agt_hc;
  } else { // for color lines or fill:
    var agt_color = 0;
    if (ratio < 0 ) agt_color = 2;
    else if ( ratio > 1 ) agt_color = this.colormap.num_cmap_colors - 1;
    else {
      agt_color = 2
        + Math.round( ratio * ( this.colormap.num_cmap_colors - 3 ) );
    }
    var agt_hc = this.colormap.getHexColor( agt_color );
    this.context.fillStyle = agt_hc;
    this.context.strokeStyle = agt_hc;
  }
//document.getElementById("debug_textarea").value +=
//  "leaving AreaGraphTool.setFgColor\n";
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.setBgColor = function( ) { // clear the canvas: 
//document.getElementById("debug_textarea").value +=
//  "entering AreaGraphTool.setBgColor\n";
  this.context.save();
  this.context.strokeStyle = this.colormap.getHexColor( 0 );
  this.context.fillRect( 0, 0, this.canvas.width, this.canvas.height );
  this.context.restore();
//document.getElementById("debug_textarea").value +=
//  "leaving AreaGraphTool.setBgColor\n";
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.setFont = function( s ) { this.context.font = s; }
AreaGraphTool.prototype.screenCoords = function( user_coords ) {
  return [
    this.canvas.width
      * ( ( user_coords[ 0 ] - this.low[ 0 ] ) / this.len[ 0 ] ),
    this.canvas.height
      * ( ( this.high[ 1 ] - user_coords[ 1 ] ) / this.len[ 1 ] )
  ];
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.userCoords = function( screen_coords) {
  return [
    this.low[ 0 ]
      + this.len[ 0 ] * ( screen_coords[ 0 ] / this.canvas.width ),
    this.high[ 1 ]
      - this.len[ 1 ] * ( screen_coords[ 1 ] / this.canvas.height )
  ];
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.putString = function( user_pos, s, align_left_right,
align_up_down, vertical ) {
  this.context.save();
  if ( align_left_right != undefined ) {
    this.context.textAlign = align_left_right; // see Flanagan p. 651
  }
  if ( align_up_down != undefined ) {
    this.context.textBaseline = align_up_down;
  }
  var agt_screen_coords = this.screenCoords( user_pos );
  if ( vertical ) {
    this.context.rotate( - Math.PI * 0.5 );
    this.context.fillText( s, - agt_screen_coords[ 1 ],
      agt_screen_coords[ 0 ] );
  } else {
    this.context.fillText( s, agt_screen_coords[ 0 ],
      agt_screen_coords[ 1 ] );
  }
  this.context.restore();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.rescale = function( low, high ) {
  this.user_low[ 0 ] = Math.min( low[ 0 ], high[ 0 ] );
  this.user_high[ 0 ] = Math.max( low[ 0 ], high[ 0 ] );
  if ( this.user_high[ 0 ] <= this.user_low[ 0 ] ) {
    this.user_high[ 0 ] += Math.abs( this.user_high[ 0 ] );
    if ( this.user_high[ 0 ] <= this.user_low[ 0 ] ) {
      this.user_high[ 0 ] += 1;
    }
  }
  this.user_low[ 1 ] = Math.min( low[ 1 ], high[ 1 ] );
  this.user_high[ 1 ] = Math.max( low[ 1 ], high[ 1 ] );
  if ( this.user_high[ 1 ] <= this.user_low[ 1 ] ) {
    this.user_high[ 1 ] += Math.abs( this.user_high[ 1 ] );
    if ( this.user_high[ 1 ] <= this.user_low[ 1 ] ) {
      this.user_high[ 1 ] += 1;
    }
  }

  var agt_dx = ( this.user_high[ 0 ] - this.user_low[ 0 ] ) * 0.05;
  this.low[ 0 ] = this.user_low[ 0 ] - agt_dx;
  this.high[ 0 ] = this.user_high[ 0 ] + agt_dx;
  this.len[ 0 ] = this.high[ 0 ] - this.low[ 0 ];
  var agt_dy = ( this.user_high[ 1 ] - this.user_low[ 1 ] ) * 0.05;
  this.low[ 1 ] = this.user_low[ 1 ] - agt_dy;
  this.high[ 1 ] = this.user_high[ 1 ] + agt_dy;
  this.len[ 1 ] = this.high[ 1 ] - this.low[ 1 ];
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.watchingMouse = function() {
  return this.watch_mouse;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.watchMouse = function() { this.watch_mouse = true; }
AreaGraphTool.prototype.ignoreMouse = function() {
  this.watch_mouse = false;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.mouseHandler = function( e ) {
  if ( ! e.currentTarget.graphtool.watch_mouse ) return;

  if ( e.type == "mousedown" || e.type == "mousemove" ||
  e.type == "mouseup") {
    var agt_screen_coords = new Array( 2 );
/*
//http://www.quirksmode.org/js/events_properties.html
    if ( e.pageX != undefined && e.pageY != undefined ) { //firefox and chrome
//    the following is incorrect for both firefox and chrome
      document.getElementById("debug_textarea").value +=
        "in AreaGraphTool.mouseHandler, e.pageX and e.pageY defined\n";
      agt_screen_coords = [ e.pageX, e.pageY ];
    } else if ( e.clientX != undefined && e.clientY != undefined ) {
      document.getElementById("debug_textarea").value +=
        "in AreaGraphTool.mouseHandler, e.clientX and e.clientY defined\n";
      agt_screen_coords = [ e.clientX + document.body.scrollLeft
                      + document.documentElement.scrollLeft,
                        e.clientY + document.body.scrollTop
                      + document.documentElement.scrollTop ];
    } else {
      alert("AreaGraphTool.mouseHandler does not know how to find screen coordinates for your browser");
    }
*/

//
//  http://miloq.blogspot.com/2011/05/coordinates-mouse-click-canvas.html
    if ( e.x != undefined && e.y != undefined ) { // chrome, microsoft
//    document.getElementById("debug_textarea").value +=
//      "in AreaGraphTool.mouseHandler, e.x and e.y defined\n";
      if ( e.pageX != undefined && e.pageY != undefined ) { // chrome
        agt_screen_coords = [ e.pageX, e.pageY ];
      } else {
        agt_screen_coords = [ e.x, e.y ];
      }
      e.currentTarget.graphtool.button = e.button; // L: 1, R:2, M: 4
    } else { // firefox
//    document.getElementById("debug_textarea").value +=
//      "in AreaGraphTool.mouseHandler, e.x or e.y undefined\n";
      agt_screen_coords = [ e.clientX + document.body.scrollLeft
                      + document.documentElement.scrollLeft,
                        e.clientY + document.body.scrollTop
                      + document.documentElement.scrollTop ];
      if ( e.button == 0 ) e.currentTarget.graphtool.button = 1;
      else if ( e.button == 1 ) e.currentTarget.graphtool.button = 4;
      else e.currentTarget.graphtool.button = 2;
    }
    agt_screen_coords = [ agt_screen_coords[ 0 ] - e.currentTarget.offsetLeft,
                          agt_screen_coords[ 1 ] - e.currentTarget.offsetTop];
//
    e.currentTarget.graphtool.mouse_coords = 
      e.currentTarget.graphtool.userCoords( agt_screen_coords );
    if ( e.type == "mousedown" ) {
      if ( e.currentTarget.graphtool.mouseDown != undefined ) {
        e.currentTarget.graphtool.mouseDown( e );
      }
    }
    if ( e.type == "mousedown" ) {
      if ( e.currentTarget.graphtool.mouseMove != undefined ) {
        e.currentTarget.graphtool.mouseMove( e );
      }
    }
    if ( e.type == "mouseup" ) {
      if ( e.currentTarget.graphtool.mouseUp != undefined ) {
        e.currentTarget.graphtool.mouseUp( e );
      }
      e.currentTarget.graphtool.watch_mouse = false;
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.getMouse = function() {
  return {
    mx : this.mouse_coords[ 0 ],
    my : this.mouse_coords[ 1 ],
    bt : this.button
  };
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.beginDrawing = function() {
  this.context.beginPath();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.endDrawing = function() { this.context.stroke(); }
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.movePen = function( user_pos ) {
  var agt_screen_coords = this.screenCoords( user_pos );
  this.context.moveTo( agt_screen_coords[ 0 ], agt_screen_coords[ 1 ] );
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.drawLine = function( user_pos ) {
  var agt_screen_coords = this.screenCoords( user_pos );
  this.context.lineTo( agt_screen_coords[ 0 ], agt_screen_coords[ 1 ] );
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.drawBoundingBox = function() {
  this.beginDrawing();
    this.context.moveTo( 0, 0 );
    this.context.lineTo( this.canvas.width, 0 );
    this.context.lineTo( this.canvas.width, this.canvas.height );
    this.context.lineTo( 0, this.canvas.height );
    this.context.lineTo( 0, 0 );
  this.endDrawing();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.drawPlusGivenCorners = function( low, high ) {
  var agt_mid = [ 0.5 * ( low[ 0 ] + high[ 0 ] ),
              0.5 * ( low[ 1 ] + high[ 1 ] ) ];
  this.beginDrawing();
    this.movePen( [ low[ 0 ], agt_mid[ 1 ] ] );
    this.drawLine( [ high[ 0 ], agt_mid[ 1 ] ] );
  this.endDrawing();
  this.beginDrawing();
    this.movePen( [ agt_mid[ 0 ], low[ 1 ] ] );
    this.drawLine( [ agt_mid[ 0 ], high[ 1 ] ] );
  this.endDrawing();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.drawPlusGivenCenter = function( center, inc ) {
  this.beginDrawing();
    this.movePen( [ center[ 0 ] - inc[ 0 ], center[ 1 ] ] );
    this.drawLine( [ center[ 0 ] + inc[ 0 ], center[ 1 ] ] );
  this.endDrawing();
  this.beginDrawing();
    this.movePen( [ center[ 0 ] , center[ 1 ] - inc[ 1 ] ] );
    this.drawLine( [ center[ 0 ] , center[ 1 ] + inc[ 1 ] ] );
  this.endDrawing();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.drawPlus = function( center, dx ) {
  if ( dx == undefined ) dx = 0.5;
  this.drawPlusGivenCenter( center,
    [ dx, dx * this.len[ 1 ] / this.len[ 0 ] ] );
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.drawCrossGivenCorners = function( low, high ) {
  this.beginDrawing();
    this.movePen( low );
    this.drawLine( high );
  this.endDrawing();
  this.beginDrawing();
    this.movePen( [ high[ 0 ], low[ 1 ] ] );
    this.drawLine( [ low[ 0 ], high[ 1 ] ] );
  this.endDrawing();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.drawCrossGivenCenter = function( center, inc ) {
  this.drawCrossGivenCorners( 
    [ center[ 0 ] - inc[ 0 ], center[ 1 ] - inc[ 1 ] ],
    [ center[ 0 ] + inc[ 0 ], center[ 1 ] + inc[ 1 ] ] );
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.drawCross = function( center, dx ) {
  if ( dx == undefined ) dx = 0.5;
  this.drawCrossGivenCenter( center, 
    [ dx, dx * this.len[ 1 ]/ this.len[ 0 ] ] );
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.drawDiamondGivenCorners = function( low, high ) {
  var agt_mid = [ 0.5 * ( low[ 0 ] + high[ 0 ] ),
              0.5 * ( low[ 1 ] + high[ 1 ] ) ];
  this.beginDrawing();
    this.movePen( [ high[ 0 ], agt_mid[ 1 ] ] );
    this.drawLine( [ agt_mid[ 0 ], high[ 1 ] ] );
    this.drawLine( [ low[ 0 ], agt_mid[ 1 ] ] );
    this.drawLine( [ agt_mid[ 0 ], low[ 1 ] ] );
    this.drawLine( [ high[ 0 ], agt_mid[ 1 ] ] );
  this.endDrawing();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.drawDiamondGivenCenter = function( center, dx ) {
  if ( dx == undefined ) dx = 0.5;
  var agt_xxhi = center[ 0 ] + dx;
  var agt_dy = dx * this.len[ 1 ] / this.len[ 0 ];
  this.beginDrawing();
    this.movePen( [ agt_xxhi, center[ 1 ] ] );
    this.drawLine( [ center[ 0 ], center[ 1 ] + agt_dy ] );
    this.drawLine( [ center[ 0 ] - dx, center[ 1 ] ] );
    this.drawLine( [ center[ 0 ], center[ 1 ] - agt_dy ] );
    this.drawLine( [ agt_xxhi, center[ 1 ] ] );
  this.endDrawing();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.drawBoxGivenCorners = function( low, high ) {
  this.beginDrawing();
    this.movePen( low );
    this.drawLine( [ high[ 0 ], low[ 1 ] ] );
    this.drawLine( high );
    this.drawLine( [ low[ 0 ], high[ 1 ] ] );
    this.drawLine( low );
  this.endDrawing();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.drawBoxGivenCenter = function( center, dx ) {
  if ( dx == undefined ) dx = 0.5;
  var agt_dy = dx * this.len[ 1 ] / this.len[ 0 ];
  this.drawBoxGivenCorners( 
    [ center[ 0 ] - dx, center[ 1 ] - agt_dy ],
    [ center[ 0 ] + dx, center[ 1 ] + agt_dy ] );
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.drawVector = function( base, head, arrow ) {
  var agt_position = ComplexMath.plus( base , head );
  this.beginDrawing();
    this.movePen( [ base.getReal(), base.getImag() ] );
    this.drawLine( [ agt_position.getReal(), agt_position.getImag() ] );
    agt_position.plusEquals( ComplexMath.times( head , arrow ) );   
    this.drawLine( [ agt_position.getReal(), agt_position.getImag() ] );
  this.endDrawing();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.colorPolygon = function( corners, ratio ) {
  var agt_color = new Number();
  if ( ratio < 0 ) agt_color = 2;
  else if ( ratio > 1 ) agt_color = this.colormap.num_cmap_colors - 1;
  else {
    agt_color =
      2 + Math.round( ratio * ( this.colormap.num_cmap_colors - 3 ) );
  }
  this.context.fillStyle = this.colormap.getHexColor( agt_color );
  var old_line_width = this.context.lineWidth;
  this.context.lineWidth = 0;
  this.context.beginPath();
    this.movePen( corners[ 0 ] );
    for ( var agt_i = 1; agt_i < corners.length; agt_i ++ ) {
      this.drawLine( corners[ agt_i ] );
    }
  this.context.closePath();
  this.context.fill();
  this.context.lineWidth = old_line_width;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//static member function
AreaGraphTool.ticSpacing = function( lo, hi ) {
  var agt_sig_digit=Math.pow(10,Math.floor(Math.LOG10E*Math.log(hi-lo)));
  var agt_digit_count=Math.ceil(hi/agt_sig_digit)-Math.ceil(lo/agt_sig_digit);
  if (agt_digit_count>10) return 2*agt_sig_digit;
  else if (agt_digit_count>=5) return agt_sig_digit;
  else if (agt_digit_count>=2) return 0.5*agt_sig_digit;
  else return 0.25*agt_sig_digit;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
AreaGraphTool.prototype.drawAxes = function() {
  var agt_xtic = AreaGraphTool.ticSpacing( this.user_low[ 0 ],
    this.user_high[ 0 ] );
  var agt_xfirst_tic = Math.floor( this.user_low[ 0 ] / agt_xtic );
  var agt_xlast_tic = Math.ceil( this.user_high[ 0 ] / agt_xtic );
  var agt_xaxis = Math.max( agt_xtic * agt_xfirst_tic,
    Math.min( 0, agt_xtic * agt_xlast_tic ));

  var agt_ytic = AreaGraphTool.ticSpacing( this.user_low[ 1 ],
    this.user_high[ 1 ] );
  var agt_yfirst_tic = Math.floor( this.user_low[ 1 ] / agt_ytic );
  var agt_ylast_tic = Math.ceil( this.user_high[ 1 ] / agt_ytic );
  var agt_yaxis = Math.max( agt_ytic * agt_yfirst_tic,
    Math.min( 0, agt_ytic * agt_ylast_tic ));

  var agt_TICK_LENGTH=0.01;
  var agt_tic_pixels = 
    agt_TICK_LENGTH*Math.max( this.canvas.width, this.canvas.height );
  var agt_tic_xlen = this.len[ 0 ] * agt_tic_pixels / this.canvas.width;
  var agt_tic_ylen = this.len[ 1 ] * agt_tic_pixels / this.canvas.height;
  this.beginDrawing();
    this.movePen( [ this.user_low[ 0 ] , agt_yaxis ] );
    this.drawLine( [ this.user_high[ 0 ], agt_yaxis ] );
  this.endDrawing();
  for ( var agt_i = agt_xfirst_tic; agt_i <= agt_xlast_tic; agt_i++ ) {
    var agt_x = agt_i * agt_xtic;
    this.beginDrawing();
      this.movePen( [ agt_x, agt_yaxis - agt_tic_ylen ] );
      this.drawLine( [ agt_x, agt_yaxis + agt_tic_ylen ] );
    this.endDrawing();
    var agt_pos = [ agt_x, agt_yaxis ];
    if ( agt_i < 0 ) {
      this.putString( agt_pos, agt_x.toString(), "right", "top" , false );
    } else if ( agt_i > 0 ) {
      this.putString( agt_pos, agt_x.toString(), "left", "top" , false );
    }
  }
  this.putString( [ agt_xlast_tic * agt_xtic, agt_yaxis ], this.xlabel, 
    ( agt_xlast_tic > 0 ? "right" : "left" ), "bottom" , false );

  this.beginDrawing();
    this.movePen( [ agt_xaxis, this.user_low[ 1 ] ] );
    this.drawLine( [ agt_xaxis, this.user_high[ 1 ] ] );
  this.endDrawing();
  for ( var agt_i = agt_yfirst_tic; agt_i <= agt_ylast_tic; agt_i++ ) {
    var agt_y = agt_i * agt_ytic;
    this.beginDrawing();
      this.movePen( [ agt_xaxis - agt_tic_xlen, agt_y ] );
      this.drawLine( [ agt_xaxis + agt_tic_xlen, agt_y ] );
    this.endDrawing();
    if ( agt_i < 0 ) {
      this.putString( [ agt_xaxis, agt_y ], agt_y.toString(), "right",
        "bottom", true );
    } else if ( agt_i > 0 ) {
      this.putString( [ agt_xaxis, agt_y ], agt_y.toString(), "right",
        "bottom" , true);
    }
  }
  this.putString( [ agt_xaxis, agt_ylast_tic * agt_ytic ], this.ylabel,
    "right", "top", true );
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//static function:
AreaGraphTool.onload = function( canvas_name ) {
}
