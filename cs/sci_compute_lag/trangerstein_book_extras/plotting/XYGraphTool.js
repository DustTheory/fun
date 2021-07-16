function XYGraphTool( c , s , xlab, ylab, xxlo, xxhi, yylo, yyhi ) {
//document.getElementById("debug_textarea").value +=
//  "entering XYGraphTool constructor, xlab,ylab = " + xlab + " "
//  + ylab + "\n";
  if ( c == undefined ) {
    alert("XYGraphTool constructor: canvas is undefined");
  }
  this.canvas = c;
  this.name = s;
  this.context = this.canvas.getContext("2d");
  this.xlabel = xlab; // x axis label
  this.ylabel = ylab; // y axis label
  this.watch_mouse = false;
  this.user_xlo = undefined; // user coordinates
  this.user_xhi = undefined;
  this.user_ylo = undefined;
  this.user_yhi = undefined;
  this.xlo = undefined; // expanded user coordinates
  this.xhi = undefined;
  this.xlen = undefined;
  this.ylo = undefined;
  this.yhi = undefined;
  this.ylen = undefined;
  this.mouseX = undefined; // user coordinates of mouse
  this.mouseY = undefined;
  this.button = undefined;
  this.mouseDown = undefined; // function to be called when mouse down
  this.mouseMove = undefined; // function to be called when mouse moves
  this.mouseUp = undefined;   // function to be called when mouse up

  this.rescale( xxlo, xxhi, yylo, yyhi );
//document.getElementById("debug_textarea").value +=
//  "leaving XYGraphTool constructor\n";
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.plotFunction = function( f, a, b, n ) {
  var xygt_xarray = new Array( n + 1 );
  var xygt_farray = new Array( n + 1 );
  var xygt_dx = ( b - a ) / n;
  var xygt_x = a;
  var xygt_y = f( a );
  xygt_xarray[ 0 ] = a;
  xygt_farray[ 0 ] = xygt_y;
  var xygt_fmin = xygt_y;
  var xygt_fmax = xygt_y;
  for ( var xygt_i = 1; xygt_i < n; xygt_i ++ ) {
    xygt_x = a + xygt_i * xygt_dx;
    xygt_y = f( xygt_x );
    xygt_xarray[ xygt_i ] = xygt_x;
    xygt_farray[ xygt_i ] = xygt_y;
    xygt_fmin = Math.min( xygt_fmin, xygt_y );
    xygt_fmax = Math.max( xygt_fmax, xygt_y );
  }
  xygt_x = b;
  xygt_y = f( b );
  xygt_xarray[ n ] = b;
  xygt_farray[ n ] = xygt_y;
  xygt_fmin = Math.min( xygt_fmin, xygt_y );
  xygt_fmax = Math.max( xygt_fmax, xygt_y );
  this.rescale( a, b, xygt_fmin, xygt_fmax );
  this.setBgColor( "white" );
  this.newPage();
  this.setFgColor( "black" );
  this.drawBoundingBox();
  this.drawAxes();
  this.setFgColor( "blue" );
  this.beginDrawing();
    this.movePen( xygt_xarray[ 0 ], xygt_farray[ 0 ] );
    for ( var xygt_i = 1; xygt_i < xygt_xarray.length; xygt_i ++ ) {
      this.drawLine( xygt_xarray[ xygt_i ], xygt_farray[ xygt_i ] );
    }
  this.endDrawing();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.listen = function() {
//document.getElementById("debug_textarea").value +=
//  "entering XYGraphTool.listen, watch_mouse = " + this.watch_mouse + "\n";
  this.canvas.graphtool = this;
  this.canvas.addEventListener( "mousedown", XYGraphTool.mouseHandler,
    false );
  this.canvas.addEventListener( "mousemove", XYGraphTool.mouseHandler,
    false );
  this.canvas.addEventListener( "mouseup", XYGraphTool.mouseHandler,
    false );
//document.getElementById("debug_textarea").value +=
//  "leaving XYGraphTool.listen\n";
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.newPage = function() {
  this.context.clearRect( 0, 0, this.canvas.width, this.canvas.height );
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.setLineWidth = function( lw ) {
  this.context.lineWidth = lw;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.setFgColor = function( css_color_string ) { 
  this.context.strokeStyle = css_color_string;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.setBgColor = function( css_color_string ) { 
  this.context.save();
  this.context.strokeStyle = css_color_string;
  this.context.fillRect( 0, 0, this.canvas.width, this.canvas.height );
  this.context.restore();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.setFont = function( s ) { this.context.font = s; }
XYGraphTool.prototype.screenCoords = function( userX, userY ) {
  return {
    screenX : this.canvas.width * ( ( userX - this.xlo ) / this.xlen ),
    screenY : this.canvas.height * ( ( this.yhi - userY ) / this.ylen )
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.userCoords = function( screenX, screenY ) {
  return {
    userX : this.xlo + this.xlen * ( screenX / this.canvas.width ),
    userY : this.yhi - this.ylen * ( screenY / this.canvas.height )
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.putString = function( x, y, s, align_left_right,
align_up_down, vertical ) {
  this.context.save();
  if ( align_left_right != undefined ) {
    this.context.textAlign = align_left_right; // see Flanagan p. 651
  }
  if ( align_up_down != undefined ) {
    this.context.textBaseline = align_up_down;
  }
  var screen_coords = this.screenCoords( x, y );
  if ( vertical ) {
    this.context.rotate( - Math.PI * 0.5 );
    this.context.fillText( s, - screen_coords.screenY,
      screen_coords.screenX );
  } else {
    this.context.fillText( s, screen_coords.screenX,
      screen_coords.screenY );
  }
  this.context.restore();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.rescale = function( xxlo, xxhi, yylo, yyhi ) {
//document.getElementById("debug_textarea").value +=
//  "entering rescale\n";
//document.getElementById("debug_textarea").value +=
//  "xxlo,xxhi = " + xxlo + " " + xxhi + "\n";
//document.getElementById("debug_textarea").value +=
//  "yylo,yyhi = " + yylo + " " + yyhi + "\n";
  this.user_xlo = Math.min( xxlo, xxhi );
  this.user_xhi = Math.max( xxlo, xxhi );
  if ( this.user_xhi <= this.user_xlo ) {
    this.user_xhi += Math.abs( this.user_xhi );
    if ( this.user_xhi <= this.user_xlo ) this.user_xhi += 1;
  }
  var xygt_xtic = XYGraphTool.ticSpacing( this.user_xlo, this.user_xhi );
//document.getElementById("debug_textarea").value +=
//  "xtic = " + xygt_xtic + "\n";
  var xygt_xfirst_tic = Math.floor( this.user_xlo / xygt_xtic );
  var xygt_xlast_tic = Math.ceil( this.user_xhi / xygt_xtic );
  if ( xygt_xlast_tic <= xygt_xfirst_tic ) {
    xygt_xlast_tic = xygt_xfirst_tic * ( 1. + 1.e-15 );
  }
  if ( xygt_xlast_tic <= xygt_xfirst_tic ) {
    xygt_xlast_tic = xygt_xfirst_tic * ( 1. - 1.e-15 );
  }
//document.getElementById("debug_textarea").value +=
//  "xygt_xfirst_tic,xygt_xlast_tic = " + xygt_xfirst_tic + " "
//  + xygt_xlast_tic + "\n";
  this.user_xlo = xygt_xtic * xygt_xfirst_tic;
  this.user_xhi = xygt_xtic * xygt_xlast_tic;
//document.getElementById("debug_textarea").value +=
//  "user_xlo,user_xhi = " + this.user_xlo + " " + this.user_xhi + "\n";

  this.user_ylo = Math.min( yylo, yyhi );
  this.user_yhi = Math.max( yylo, yyhi );
  if ( this.user_yhi <= this.user_ylo ) {
    this.user_yhi += Math.abs( this.user_yhi );
    if ( this.user_yhi <= this.user_ylo ) this.user_yhi += 1;
  }
//document.getElementById("debug_textarea").value +=
//  "user_ylo,user_yhi = " + this.user_ylo + " " + this.user_yhi + "\n";
  var xygt_ytic = XYGraphTool.ticSpacing( this.user_ylo, this.user_yhi );
  var xygt_yfirst_tic = Math.floor( this.user_ylo / xygt_ytic );
  var xygt_ylast_tic = Math.ceil( this.user_yhi / xygt_ytic );
  if ( xygt_ylast_tic <= xygt_yfirst_tic ) {
    xygt_ylast_tic = xygt_yfirst_tic * ( 1. + 1.e-15 );
  }
  if ( xygt_ylast_tic <= xygt_yfirst_tic ) {
    xygt_ylast_tic = xygt_yfirst_tic * ( 1. - 1.e-15 );
  }
  this.user_ylo = xygt_ytic * xygt_yfirst_tic;
  this.user_yhi = xygt_ytic * xygt_ylast_tic;
//document.getElementById("debug_textarea").value +=
//  "user_ylo,user_yhi = " + this.user_ylo + " " + this.user_yhi + "\n";

  var xygt_dx = ( this.user_xhi - this.user_xlo ) * 0.05;
  this.xlo = this.user_xlo - xygt_dx;
  this.xhi = this.user_xhi + xygt_dx;
  this.xlen = this.xhi - this.xlo;
  var dy = ( this.user_yhi - this.user_ylo ) * 0.05;
  this.ylo = this.user_ylo - dy;
  this.yhi = this.user_yhi + dy;
  this.ylen = this.yhi - this.ylo;
//document.getElementById("debug_textarea").value +=
//  "ylo,yhi = " + this.ylo + " " + this.yhi + "\n";
//document.getElementById("debug_textarea").value +=
//  "leaving rescale\n";
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.watchingMouse = function() {
  return this.watch_mouse;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.watchMouse = function() { 
//document.getElementById("debug_textarea").value +=
//  "entering watchMouse\n";
  this.watch_mouse = true;
//document.getElementById("debug_textarea").value +=
//  "leaving watchMouse, watch_mouse = " + this.watch_mouse + "\n";
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.ignoreMouse = function() { this.watch_mouse = false; }
XYGraphTool.mouseHandler = function( e ) {
  if ( ! e.currentTarget.graphtool.watch_mouse ) return;
  if ( e.type == "mousedown" || e.type == "mousemove" ||
  e.type == "mouseup") {
//  document.getElementById("debug_textarea").value +=
//    "in XYGraphTool.mouseHandler\n";
    var sx = new Number();
    var sy = new Number();
    if ( e.x != undefined && e.y != undefined ) { // chrome, microsoft
//    http://www.quirksmode.org/js/events_properties.html
//    http://miloq.blogspot.com/2011/05/coordinates-mouse-click-canvas.html
      if ( e.pageX != undefined && e.pageY != undefined ) { // chrome
        sx = e.pageX;
        sy = e.pageY;
      } else { // who knows?
        sx = e.x;
        sy = e.y;
      }
      e.currentTarget.graphtool.button = e.button; // L: 1, R:2, M: 4
//    document.getElementById("debug_textarea").value +=
//      "mouseHandler microsoft sx,sy = " + sx + " " + sy + "\n";
    } else { // firefox
      sx = e.clientX + document.body.scrollLeft
         + document.documentElement.scrollLeft;
      sy = e.clientY + document.body.scrollTop
         + document.documentElement.scrollTop;
      if ( e.button == 0 ) e.currentTarget.graphtool.button = 1;
      else if ( e.button == 1 ) e.currentTarget.graphtool.button = 4;
      else e.currentTarget.graphtool.button = 2;
//    document.getElementById("debug_textarea").value +=
//      "mouseHandler firefox sx,sy = " + sx + " " + sy + "\n";
    }
    sx -= e.currentTarget.offsetLeft;
    sy -= e.currentTarget.offsetTop;
//  document.getElementById("debug_textarea").value +=
//    "mouseHandler sx,sy = " + sx + " " + sy + "\n";
    var user_coords = e.currentTarget.graphtool.userCoords( sx, sy );
    e.currentTarget.graphtool.mouseX = user_coords.userX;
    e.currentTarget.graphtool.mouseY = user_coords.userY;
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
//      document.getElementById("debug_textarea").value +=
//        "mouseHandler sx,sy = " + sx + " " + sy + "\n";
//      document.getElementById("debug_textarea").value +=
//        "user_coords = " + user_coords.userX + " " + user_coords.userY
//        + "\n";
        e.currentTarget.graphtool.mouseUp( e );
      }
      e.currentTarget.graphtool.watch_mouse = false;
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.getMouse = function() {
  return {
    mx : this.mouseX,
    my : this.mouseY,
    bt : this.button
  };
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.beginDrawing = function() {
  this.context.beginPath();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.endDrawing = function() { this.context.stroke(); }
XYGraphTool.prototype.movePen = function( x, y ) {
  var screen_coords = this.screenCoords( x, y );
  this.context.moveTo( screen_coords.screenX, screen_coords.screenY );
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.drawLine = function( x, y ) {
  var screen_coords = this.screenCoords( x, y );
  this.context.lineTo( screen_coords.screenX, screen_coords.screenY );
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.drawBoundingBox = function() {
  this.beginDrawing();
    this.context.moveTo( 0, 0 );
    this.context.lineTo( this.canvas.width, 0 );
    this.context.lineTo( this.canvas.width, this.canvas.height );
    this.context.lineTo( 0, this.canvas.height );
    this.context.lineTo( 0, 0 );
  this.endDrawing();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.drawPlusGivenCorners = function( xxlo, yylo,
xxhi, yyhi ) {
  var xmid = 0.5 * ( xxlo + xxhi );
  var ymid = 0.5 * ( yylo + yyhi );
  this.beginDrawing();
    this.movePen( xxlo, ymid );
    this.drawLine( xxhi, ymid );
  this.endDrawing();
  this.beginDrawing();
    this.movePen( xmid, yylo );
    this.drawLine( xmid, yyhi );
  this.endDrawing();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.drawPlusGivenCenter = function( x, y, dx, dy ) {
  this.beginDrawing();
    this.movePen( x - dx, y );
    this.drawLine( x + dx, y );
  this.endDrawing();
  this.beginDrawing();
    this.movePen( x , y - dy );
    this.drawLine( x , y + dy );
  this.endDrawing();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.drawPlus = function( x, y, dx ) {
  if ( dx == undefined ) dx = 0.5;
  this.drawPlusGivenCenter( x, y, dx, dx * this.ylen / this.xlen );
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.drawCrossGivenCorners = function( xxlo, yylo,
xxhi, yyhi ) {
  this.beginDrawing();
    this.movePen( xxlo, yylo );
    this.drawLine( xxhi, yyhi );
  this.endDrawing();
  this.beginDrawing();
    this.movePen( xxhi, yylo );
    this.drawLine( xxlo, yyhi );
  this.endDrawing();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.drawCrossGivenCenter = function( x, y, dx, dy ) {
  this.drawCrossGivenCorners( x - dx, y - dy, x + dx, y + dy );
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.drawCross = function( x, y, dx ) {
  if ( dx == undefined ) dx = 0.5;
  this.drawCrossGivenCenter( x, y, dx, dx * this.ylen/ this.xlen );
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.drawDiamondGivenCorners = function( xxlo, yylo,
xxhi, yyhi ) {
  var xmid = 0.5 * ( xxlo + xxhi );
  var ymid = 0.5 * ( yylo + yyhi );
  this.beginDrawing();
    this.movePen( xxhi, ymid );
    this.drawLine( xmid, yyhi );
    this.drawLine( xxlo, ymid );
    this.drawLine( xmid, yylo );
    this.drawLine( xxhi, ymid );
  this.endDrawing();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.drawDiamondGivenCenter = function( x, y, dx ) {
  if ( dx == undefined ) dx = 0.5;
  var xxhi = x + dx;
  var dy = dx * this.ylen / this.xlen;
  this.beginDrawing();
    this.movePen( xxhi, y );
    this.drawLine( x, y + dy );
    this.drawLine( x - dx, y );
    this.drawline( x, y - dy );
    this.drawLine( xxhi, y );
  this.endDrawing();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.drawBoxGivenCorners = function( xxlo, yylo,
xxhi, yyhi ) {
  this.beginDrawing();
    this.movePen( xxlo, yylo );
    this.drawLine( xxhi, yylo );
    this.drawline( xxhi, yyhi );
    this.drawLine( xxlo, yyhi );
    this.drawLine( xxlo, yylo );
  this.endDrawing();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.drawBoxGivenCenter = function( x, y, dx ) {
  if ( dx == undefined ) dx = 0.5;
  var dy = dx * this.ylen / this.xlen;
  this.drawBoxGivenCorners( x - dx, y - dy, x + dx, y + dy );
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
XYGraphTool.prototype.drawVector = function( base, head, arrow ) {
  var position = ComplexMath.plus( base , head );
  this.beginDrawing();
    this.movePen( base.getReal(), base.getImag() );
    this.drawline( position.getReal(), position.getImag() );
    position.plusEquals( ComplexMath.times( head , arrow ) );   
    this.drawLine( position.getReal(), position.getImag() );
  this.endDrawing();
}
*/
/*
XYGraphTool.prototype.colorPolygon = function() {
}
*/
/*
XYGraphTool.prototype.colorQuad = function() {
}
*/
//static member function
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.ticSpacing = function( lo, hi ) {
//document.getElementById("debug_textarea").value +=
//  "entering XYGraphTool.ticSpacing, lo,hi = " + lo + " " + hi + "\n";
  var sig_digit=Math.pow(10,Math.floor(Math.LOG10E*Math.log(hi-lo)));
  var digit_count=Math.ceil(hi/sig_digit)-Math.ceil(lo/sig_digit);
//document.getElementById("debug_textarea").value +=
//  "  sig_digit,digit_count = " + sig_digit + " " + digit_count + "\n";
//document.getElementById("debug_textarea").value +=
//  "leaving XYGraphTool.ticSpacing\n";
  if (digit_count>10) return 2*sig_digit;
  else if (digit_count>=5) return sig_digit;
  else if (digit_count>=2) return 0.5*sig_digit;
  else return 0.25*sig_digit;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
XYGraphTool.prototype.drawAxes = function() {
//document.getElementById("debug_textarea").value +=
//  "entering XYGraphTool.drawAxes, xlabel,ylab = " + this.xlabel + " "
//  + this.ylabel + "\n";
  var xygt_xtic = XYGraphTool.ticSpacing( this.user_xlo, this.user_xhi );
  var xygt_ytic = XYGraphTool.ticSpacing( this.user_ylo, this.user_yhi );
//document.getElementById("debug_textarea").value +=
//  "  xlo,xhi = " + this.user_xlo + " " + this.user_xhi + "\n";
//document.getElementById("debug_textarea").value +=
//  "  ylo,yhi = " + this.user_ylo + " " + this.user_yhi + "\n";
//document.getElementById("debug_textarea").value +=
//  "  xtic,ytic = " + xygt_xtic + " " + xygt_ytic + "\n";

  if ( xygt_xtic <=0 || xygt_ytic <= 0 ) return;
  var xygt_xfirst_tic = Math.floor( this.user_xlo / xygt_xtic );
  var xygt_xlast_tic = Math.ceil( this.user_xhi / xygt_xtic );
  var xaxis = Math.max( xygt_xtic * xygt_xfirst_tic,
    Math.min( 0, xygt_xtic * xygt_xlast_tic ));
//document.getElementById("debug_textarea").value +=
//  "  xygt_xfirst_tic,xygt_xlast_tic = " + xygt_xfirst_tic + " "
//  + xygt_xlast_tic + "\n";

  var xygt_yfirst_tic = Math.floor( this.user_ylo / xygt_ytic );
  var xygt_ylast_tic = Math.ceil( this.user_yhi / xygt_ytic );
  var yaxis = Math.max( xygt_ytic * xygt_yfirst_tic,
    Math.min( 0, xygt_ytic * xygt_ylast_tic ));
//document.getElementById("debug_textarea").value +=
//  "  xygt_yfirst_tic,xygt_ylast_tic,yaxis = " + xygt_yfirst_tic + " "
//  + xygt_ylast_tic + " " //  + yaxis + "\n";

  var TICK_LENGTH=0.01;
  var tic_pixels = 
    TICK_LENGTH*Math.max( this.canvas.width, this.canvas.height );
  var tic_xlen = this.xlen * tic_pixels / this.canvas.width;
  var tic_ylen = this.ylen * tic_pixels / this.canvas.height;
  this.beginDrawing();
    this.movePen( this.user_xlo , yaxis );
    this.drawLine( this.user_xhi, yaxis );
  this.endDrawing();
  for ( var xygt_i = xygt_xfirst_tic; xygt_i <= xygt_xlast_tic; ) {
    var xygt_x = xygt_i * xygt_xtic;
    this.beginDrawing();
      this.movePen( xygt_x, yaxis - tic_ylen );
      this.drawLine( xygt_x, yaxis + tic_ylen );
    this.endDrawing();
    if ( xygt_i < 0 ) {
      this.putString( xygt_x, yaxis, xygt_x.toString(), "right", "top" ,
        false );
    } else if ( xygt_i > 0 ) {
      this.putString( xygt_x, yaxis, xygt_x.toString(), "left", "top" ,
        false );
    }
//  JavaScript uses floats for loop indices
//  if the index becomes too large, incrementing by one will round to
//    the previous index
    var ii = xygt_i + 1;
    if ( ii <= xygt_i ) ii = xygt_i + 2;
    if ( ii <= xygt_i ) ii = xygt_i + 4;
    if ( ii <= xygt_i ) ii = xygt_i + 8;
    if ( ii <= xygt_i ) ii = xygt_i + 16;
    xygt_i = ii;
  }
//document.getElementById("debug_textarea").value +=
//  "  after x loop\n";
  this.putString( xygt_xlast_tic*xygt_xtic, yaxis, this.xlabel, 
    ( xygt_xlast_tic > 0 ? "right" : "left" ), "bottom" , false );

  this.beginDrawing();
    this.movePen( xaxis, this.user_ylo );
    this.drawLine( xaxis, this.user_yhi );
  this.endDrawing();
  for ( var xygt_i = xygt_yfirst_tic; xygt_i <= xygt_ylast_tic; xygt_i++ ) {
    var xygt_y = xygt_i * xygt_ytic;
    this.beginDrawing();
      this.movePen( xaxis - tic_xlen, xygt_y );
      this.drawLine( xaxis + tic_xlen, xygt_y );
    this.endDrawing();
    if ( xygt_i < 0 ) {
      this.putString( xaxis, xygt_y, xygt_y.toString(), "right", "bottom",
        true );
    } else if ( xygt_i > 0 ) {
      this.putString( xaxis, xygt_y, xygt_y.toString(), "right", "bottom" ,
        true );
    }
    ii = xygt_i + 1;
    if ( ii <= xygt_i ) ii = xygt_i + 2;
    if ( ii <= xygt_i ) ii = xygt_i + 4;
    if ( ii <= xygt_i ) ii = xygt_i + 8;
    if ( ii <= xygt_i ) ii = xygt_i + 16;
    xygt_i = ii;
  }
//document.getElementById("debug_textarea").value +=
//  "  after y loop\n";
  this.putString( xaxis, xygt_ylast_tic * xygt_ytic, this.ylabel, "right",
    "top", true );
//document.getElementById("debug_textarea").value +=
//  "leaving XYGraphTool.drawAxes\n";
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//static function:
XYGraphTool.onload = function( canvas_name ) { // see p 785 for mouse events
/*
  document.getElementById(canvas_name).addEventListener("mousedown",
    this.mouseHandler,false); //see p. 463 for 3rd arg discussion
  document.getElementById(canvas_name).addEventListener("mousemove",
    this.mouseHandler,false); //see p. 463 for 3rd arg discussion
  document.getElementById(canvas_name).addEventListener("mouseup",
    this.mouseHandler,false); //see p. 463 for 3rd arg discussion
*/
}
