//got help from
//  https://github.com/mdn/webgl-examples/blob/gh-pages/tutorial/sample7/webgl-demo.js
glMatrixArrayType = ( typeof Float32Array !="undefined" ?  Float32Array :
  ( typeof WebGLFloatArray != "undefined" ?  WebGLFloatArray : Array ) );
// W3C: see http://www.quirksmode.org/js/events_properties.html
var LEFT_BUTTON = 0;
var MIDDLE_BUTTON = 1;
var RIGHT_BUTTON = 2;
var TRACKBALLSIZE = 1;
var GRAPH_FUDGE = 0.05;
var AMBIENT_COLOR = [ .6, .6, .6 ];
var DIRECTIONAL_COLOR = [ .4, .4, .4 ];
var LIGHTING_DIRECTION = [ -.5, -.5, -.707 ];
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
function loadShader(gl, type, source) {
  const shader = gl.createShader(type);
  gl.shaderSource(shader, source);
  gl.compileShader(shader);
  if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
    alert('An error occurred compiling the shaders: ' + gl.getShaderInfoLog(shader));
    document.getElementById("debug_textarea").value +=
      "error in loadShader\n" + gl.getShaderInfoLog(shader);
    gl.deleteShader(shader);
    return null;
  }
  return shader;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
function VolGraphTool( c , s , p ) {
//document.getElementById("debug_textarea").value +=
//  "entering VolGraphTool constructor\n";
  if ( c == undefined ) {
    alert("VolGraphTool constructor: canvas is undefined");
  }
  this.canvas = c;
  this.colormap = new Colormap( p );
//this.name = s;
//this.xlabel = "x"; // x axis label
//this.ylabel = "y"; // y axis label
//this.zlabel = "z"; // z axis label
  this.user_low = new Array( 3 ); // user coordinates
  this.user_high = new Array( 3 );
  this.user_len = new Array( 3 );
  this.user_center = new Array( 3 );
  this.new_screen_coords = new Array( 2 ); // origin at top left
  this.old_screen_coords = new Array( 2 );
  this.old_unit_sphere_coords = vec3.create();
  this.new_unit_sphere_coords = vec3.create();
  this.current_quat = new Array( 4 );
  this.old_quat = new Array( 4 );
  this.button_pressed = false;
  this.button = undefined;
  this.use_lighting = false;
/*
//this.context_lost = undefined; // function to be called when context lost
//this.context_restored = undefined; // called when context restored
*/
  this.shader_program = undefined;
  this.current_vertices = undefined;//vertices for drawing line or surface
  this.current_normals = undefined;//normals at vertices for lighting
  this.current_colors = undefined; // colors for drawing a surface
  this.current_vertex_indices = undefined; // element indices for drawing
  this.current_vertex_number = 0;
  this.mv_matrix = mat4.create(); //model view: object coords -> eye coords
  this.p_matrix = mat4.create();
  this.mv_matrix_inverse = mat4.create();
  this.p_matrix_inverse = mat4.create();

  this.surface_drawing_array = new Array();
  this.curve_drawing_array = new Array();
  this.bounding_box_array = new Array();
  this.clip_bounding_box_array = new Array();
  this.clip_direction = undefined;
  this.clip_edge_vertex = new Array( 3 );
  this.cursor_position = new Array( 4 );
  this.crosshair_array = new Array();

  this.f3d = undefined;
  this.fmin = Number.MAX_VALUE;
  this.fmax = - Number.MAX_VALUE;
  this.ncell = new Array( 3 );
  this.drawing_clip_planes = true;

  this.gl = null;
  try {
    this.gl = this.canvas.getContext("webgl") ||
              this.canvas.getContext("experimental-webgl");
  } catch(e) {
  }
  if ( ! this.gl ) {
    alert("Unable to initialize WebGL. Your browser may not support it.");
  }

//the GLSL reference manual is at
//  https://www.khronos.org/registry/OpenGL/specs/gl/GLSLangSpec.4.40.pdf
//the following was suggested by
//  https://github.com/mdn/webgl-examples/blob/gh-pages/tutorial/sample7/webgl-demo.js

//Vertex shader program
  const vsSource = ` 
    precision highp float;
    attribute vec3 aVertexPosition;
    attribute vec4 aVertexColor;
    attribute vec3 aVertexNormal;
    uniform mat4 uModelViewMatrix;
    uniform mat4 uProjectionMatrix;
    uniform mat4 uNormalMatrix;
    uniform vec3 uAmbientColor;
    uniform vec3 uDirectionalColor;
    uniform vec3 uLightingDirection;
    varying vec4 vColor;
    varying vec3 vLighting;
    void main(void) {
      gl_Position = uProjectionMatrix * uModelViewMatrix * vec4( aVertexPosition, 1. );
      vColor = aVertexColor;
      highp vec4 transformedNormal = uNormalMatrix * vec4(aVertexNormal, 1.0);
      highp float directional = max( dot(transformedNormal.xyz, uLightingDirection), 0.0 );  
      vLighting = uAmbientColor + (uDirectionalColor * directional);
    }
  `;
//Fragment shader program: vector product performed componentwise
  const fsSource = `
    precision highp float;
    varying vec4 vColor;
    varying vec3 vLighting;
    void main(void) {
      gl_FragColor = vec4( vColor.rgb * vLighting, vColor.a);
    }     
  `;

//Initialize a shader program; this is where all the lighting
//for the vertices and so forth is established.
  const vertex_shader =
    loadShader(this.gl, this.gl.VERTEX_SHADER, vsSource);
  const fragment_shader =
    loadShader(this.gl, this.gl.FRAGMENT_SHADER, fsSource);

  this.shader_program = this.gl.createProgram();
  this.gl.attachShader( this.shader_program, vertex_shader );
  this.gl.attachShader( this.shader_program, fragment_shader );
  this.gl.linkProgram( this.shader_program );
  if ( !this.gl.getProgramParameter( this.shader_program,
  this.gl.LINK_STATUS ) ) {
    alert('Unable to initialize the shader program: ' + this.gl.getProgramInfoLog(this.shader_program));
    return;
  }

  this.programInfo = {
    attribLocations: {
      vertexColor: this.gl.getAttribLocation(this.shader_program, 'aVertexColor'),
      vertexPosition: this.gl.getAttribLocation(this.shader_program, 'aVertexPosition'),
      vertexNormal: this.gl.getAttribLocation(this.shader_program, 'aVertexNormal'),
    },
    uniformLocations: {
      projectionMatrix: this.gl.getUniformLocation(this.shader_program, 'uProjectionMatrix'),
      modelViewMatrix: this.gl.getUniformLocation(this.shader_program, 'uModelViewMatrix'),
      normalMatrix: this.gl.getUniformLocation(this.shader_program, 'uNormalMatrix'),
      ambientColor: this.gl.getUniformLocation(this.shader_program, 'uAmbientColor'),
      lightingDirection: this.gl.getUniformLocation(this.shader_program, 'uLightingDirection'),
      directionalColor: this.gl.getUniformLocation(this.shader_program, 'uDirectionalColor'),
    },
  };
//this.gl.useProgram( this.shader_program );

//this.shader_program.vertexPositionAttribute =
//  this.gl.getAttribLocation( this.shader_program, "aVertexPosition" );
//this.gl.enableVertexAttribArray(
//  this.shader_program.vertexPositionAttribute );
//this.shader_program.vertexNormalAttribute =
//  this.gl.getAttribLocation( this.shader_program, "aVertexNormal" );
//this.gl.enableVertexAttribArray( this.shader_program, 
//  this.shader_program. vertexNormalAttribute );
//this.shader_program.vertexColorAttribute =
//  this.gl.getAttribLocation( this.shader_program, "aVertexColor" );
//this.gl.enableVertexAttribArray(
//  this.shader_program.vertexColorAttribute );
//this.shader_program.pMatrixUniform =
//  this.gl.getUniformLocation( this.shader_program, "uPMatrix" );
//this.shader_program.mvMatrixUniform =
//  this.gl.getUniformLocation( this.shader_program, "uMVMatrix" );
//this.shader_program.nMatrixUniform =
//  this.gl.getUniformLocation( this.shader_program, "uNMatrix" );
//this.shader_program.lightingDirectionUniform =
//  this.gl.getUniformLocation( this.shader_program,
//  "uLightingDirection" );
//this.shader_program.ambientColorUniform =
//  this.gl.getUniformLocation( this.shader_program, "uAmbientColor" );
//this.shader_program.directionalColorUniform =
//  this.gl.getUniformLocation( this.shader_program,
//  "uDirectionalColor" );

  var rt3 = Math.sqrt( 3 );
  mat4.ortho( - rt3, rt3,  -rt3, rt3,  -rt3, rt3,  this.p_matrix );

/*
//this.canvas.my_graphtool = this;
//this.canvas.addEventListener( 'webglcontextlost', this.handleContextLost,
//  false );
//this.canvas.addEventListener( 'webglcontextrestored', 
//  this.handleContextRestored, false );
*/
  this.gl.enable( this.gl.DEPTH_TEST );
//this.gl.depthFunc( this.gl.LEQUAL );

//this.gl.frontFace( this.gl.CCW );//default: front = counter-clock wind
//this.gl.enable( this.gl.CULL_FACE );// not default
//this.gl.cullFace( this.gl.BACK );//default: back-facing triangle not draw
//document.getElementById("debug_textarea").value +=
//  "leaving VolGraphTool constructor\n";
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.draw2DFunction = function( f, l, h, nc ) {
  this.use_lighting = true;
  this.drawing_clip_planes = false;

  var vgt_low = new Array( 3 );
  var vgt_high = new Array( 3 );
  for ( var vgt_d=0; vgt_d < 2; vgt_d ++ ) {
    vgt_low[ vgt_d ] = Math.min( l[ vgt_d ], h[ vgt_d ] );
    vgt_high[ vgt_d ] = Math.max( l[ vgt_d ], h[ vgt_d ] );
    if ( vgt_high[ 0 ] <= vgt_low[ 0 ] ) { 
      vgt_high[ 0 ] += Math.abs( vgt_high[ 0 ] );
      if ( vgt_high[ 0 ] <= vgt_low[ 0 ] ) {
        vgt_high[ 0 ] += 1;
      }
    }
  }
  var vgt_ncorn = [ nc[ 0 ] + 1, nc[ 1 ] + 1 ];
  var vgt_farray = new Array( vgt_ncorn[ 0 ] * vgt_ncorn[ 1 ] );
  var vgt_len = [ ( vgt_high[ 0 ] - vgt_low[ 0 ] ) / nc[ 0 ],
          ( vgt_high[ 1 ] - vgt_low[ 1 ] ) / nc[ 1 ] ];
  var vgt_zmin = Number.MAX_VALUE;
  var vgt_zmax = - Number.MAX_VALUE;
  for ( var vgt_j = 0; vgt_j <= nc[ 0 ]; vgt_j ++ ) {
    var vgt_y = vgt_low[ 1 ] + vgt_len[ 1 ] * vgt_j;
    for ( var vgt_i = 0; vgt_i <= nc[ 0 ]; vgt_i ++ ) {
      var vgt_x = vgt_low[ 0 ] + vgt_len[ 0 ] * vgt_i;
      var vgt_z = f( vgt_x, vgt_y );
      vgt_farray[ vgt_i + vgt_j * vgt_ncorn[ 0 ] ] = vgt_z;
      vgt_zmin = Math.min( vgt_zmin, vgt_z );
      vgt_zmax = Math.max( vgt_zmax, vgt_z );
    }
  }
  if ( vgt_zmin >= vgt_zmax) {
    vgt_zmax = Math.abs( vgt_zmax );
    vgt_zmin = - vgt_zmax;
  }
  if ( vgt_zmin >= vgt_zmax ) {
    vgt_zmax = 1;
    vgt_zmin = 0.
  }
  vgt_low[ 2 ] = vgt_zmin;
  vgt_high[ 2 ] = vgt_zmax;

  this.rescale( vgt_low, vgt_high );
  this.current_quat = [ 0, 0, 0, 1 ]; // z toward viewer, x to right
  quat4.toMat4( this.current_quat, this.mv_matrix );
  this.clip_low =
    [ this.user_low[ 0 ], this.user_low[ 1 ], this.user_low[ 2 ] ];
  this.clip_high =
    [ this.user_high[ 0 ], this.user_high[ 1 ], this.user_high[ 2 ] ];
  this.makeOriginBox();
  this.makeBoundingBox();

  this.beginDrawing();
    var vgt_ratios = new Array( 4 );
    var vgt_corners = new Array( 4 );
    for ( var vgt_k = 0; vgt_k < 4; vgt_k ++ ) {
      vgt_corners[ vgt_k ] = new Array( 3 );
    }
    for ( var vgt_j = 0; vgt_j < nc[ 0 ]; vgt_j ++ ) {
      vgt_corners[ 1 ][ 1 ] = vgt_corners[ 0 ][ 1 ] =
        vgt_low[ 1 ] + vgt_len[ 1 ] * vgt_j;
      vgt_corners[ 3 ][ 1 ] = vgt_corners[ 2 ][ 1 ] =
        vgt_corners[ 1 ][ 1 ] + vgt_len[ 1 ];
      for ( var vgt_i = 0; vgt_i < nc[ 0 ]; vgt_i ++ ) {
        vgt_corners[ 3 ][ 0 ] = vgt_corners[ 0 ][ 0 ] =
          vgt_low[ 0 ] + vgt_len[ 0 ] * vgt_i;
        vgt_corners[ 2 ][ 0 ] = vgt_corners[ 1 ][ 0 ] =
          vgt_corners[ 0 ][ 0 ] + vgt_len[ 0 ];
        vgt_corners[ 0 ][ 2 ] = vgt_farray[ vgt_i + vgt_j * vgt_ncorn[ 0 ] ];
        vgt_corners[ 1 ][ 2 ] = vgt_farray[ vgt_i + 1 + vgt_j * vgt_ncorn[ 0 ] ];
        vgt_corners[ 2 ][ 2 ] =
          vgt_farray[ vgt_i + 1 + ( vgt_j + 1 ) * vgt_ncorn[ 0 ] ];
        vgt_corners[ 3 ][ 2 ] =
          vgt_farray[ vgt_i + ( vgt_j + 1 ) * vgt_ncorn[ 0 ] ];
        for ( var vgt_k = 0; vgt_k < 4; vgt_k ++ ) {
          vgt_ratios[ vgt_k ] =
            ( vgt_corners[ vgt_k ][ 2 ] - vgt_zmin ) / ( vgt_zmax - vgt_zmin );
        }
        this.colorQuadrilateral(
          vgt_corners[ 0 ], vgt_corners[ 1 ], vgt_corners[ 2 ], vgt_corners[ 3 ],
          vgt_ratios[ 0 ], vgt_ratios[ 1 ], vgt_ratios[ 2 ],
          vgt_ratios[ 3 ] );
      }
    }
  this.endDrawing();
  this.newPage();
  this.expose();
  this.watchMouse();
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.draw3DIsoSurfaces = function( f, l, h, nc, fval ) {
//document.getElementById("debug_textarea").value +=
//  "entering draw3DIsoSurfaces\n";
  this.drawing_clip_planes = false;
  this.use_lighting = true;
  this.rescale( l, h );
  this.current_quat = [ .5, .5, .5, .5 ]; // z up, y to right
  quat4.toMat4( this.current_quat, this.mv_matrix );
  this.clip_low =
    [ this.user_low[ 0 ], this.user_low[ 1 ], this.user_low[ 2 ] ];
  this.clip_high =
    [ this.user_high[ 0 ], this.user_high[ 1 ], this.user_high[ 2 ] ];
  this.makeOriginBox();
  this.makeBoundingBox();
  this.f3d = f;
  this.ncell = [ nc[ 0 ], nc[ 1 ], nc[ 2 ] ];
  this.fmin = Number.MAX_VALUE;
  this.fmax = - Number.MAX_VALUE;
  var vgt_inc = new Array( 3 );
  for ( var vgt_i = 0; vgt_i < 3; vgt_i ++ ) {
    vgt_inc[ vgt_i ] = this.user_len[ vgt_i ] / this.ncell[ vgt_i ];
  }
  for ( var vgt_k = 0; vgt_k <= nc[ 2 ]; vgt_k ++ ) {
    var vgt_z = this.user_low[ 2 ] + vgt_inc[ 2 ] * vgt_k;
    for ( var vgt_j = 0; vgt_j <= nc[ 1 ]; vgt_j ++ ) {
      var vgt_y = this.user_low[ 1 ] + vgt_inc[ 1 ] * vgt_j;
      for ( var vgt_i = 0; vgt_i <= nc[ 0 ]; vgt_i ++ ) {
        var vgt_x = this.user_low[ 0 ] + vgt_inc[ 0 ] * vgt_i;
        var fv = f( vgt_x, vgt_y, vgt_z );
        this.fmin = Math.min( this.fmin, fv );
        this.fmax = Math.max( this.fmax, fv );
      }
    }
  }
  if ( this.fmin >= this.fmax ) {
    this.fmax = Math.abs( this.fmax );
    this.fmin = - this.fmax;
  }
  if ( this.fmin >= this.fmax ) {
    this.fmax = 1;
    this.fmin = 0;
  }
//fdif = ( this.fmax - this.fmin ) / ( ns + 1 );
//fval = this.fmin + fdif;
//for ( s = 0; s < ns; s ++, fval += fdif ) {
    this.drawIsoSurface( fval );
//}
  this.newPage();
  this.expose();
  this.watchMouse();
//document.getElementById("debug_textarea").value +=
//  "leaving draw3DIsoSurfaces\n";
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//need to save f3d, ncell, fmin, fmax
VolGraphTool.prototype.drawIsoSurface = function( fv ) {
//document.getElementById("debug_textarea").value +=
//  "entering drawIsoSurface\n";
  this.beginDrawing();
    var vgt_ratio = Math.max( 0, Math.min( 1,
      ( fv - this.fmin ) / ( this.fmax - this.fmin ) ) );
    var vgt_values = new Array( 2 ); // function values at cell corners
    var vgt_corners = new Array( 2 );// user coordinates of cell corners
    for ( var vgt_h0 = 0; vgt_h0 < 2; vgt_h0 ++ ) {
      vgt_corners[ vgt_h0 ] = new Array( 2 );
      vgt_values[ vgt_h0 ] = new Array( 2 );
      for ( var vgt_h1 = 0; vgt_h1 < 2; vgt_h1 ++ ) {
        vgt_corners[ vgt_h0 ][ vgt_h1 ] = new Array( 2 );
        vgt_values[ vgt_h0 ][ vgt_h1 ] = new Array( 2 );
        for ( var vgt_h2 = 0; vgt_h2 < 2; vgt_h2 ++ ) {
          vgt_corners[ vgt_h0 ][ vgt_h1 ][ vgt_h2 ] = new Array( 3 );
          vgt_values[ vgt_h0 ][ vgt_h1 ][ vgt_h2 ] = new Number();
        }
      }
    }
    var vgt_intersections = new Array( 3 ); // edge location of iso surface
    for ( var vgt_d = 0; vgt_d < 3; vgt_d ++ ) {
      vgt_intersections[ vgt_d ] = new Array( 2 );
      for ( var vgt_h1 = 0; vgt_h1 < 2; vgt_h1 ++ ) {
        vgt_intersections[ vgt_d ][ vgt_h1 ] = new Array( 2 );
        for ( var vgt_h2 = 0; vgt_h2 < 2; vgt_h2 ++ ) {
          vgt_intersections[ vgt_d ][ vgt_h1 ][ vgt_h2 ] = new Array( 3 );
        }
      }
    }
    var vgt_inc = new Array( 3 );
    for ( var vgt_d = 0; vgt_d < 3; vgt_d ++ ) {
      vgt_inc[ vgt_d ] = this.user_len[ vgt_d ] / this.ncell[ vgt_d ];
    }
    var vgt_hl = new Array( 3 );
    var vgt_hr = new Array( 3 );
    for ( var vgt_kkk = 0; vgt_kkk < this.ncell[ 2 ]; vgt_kkk ++ ) {
      vgt_corners[ 0 ][ 0 ][ 0 ][ 2 ] = vgt_corners[ 1 ][ 0 ][ 0 ][ 2 ]
        = vgt_corners[ 0 ][ 1 ][ 0 ][ 2 ] = vgt_corners[ 1 ][ 1 ][ 0 ][ 2 ]
        = this.user_low[ 2 ] + vgt_inc[ 2 ] * vgt_kkk;
      vgt_corners[ 0 ][ 0 ][ 1 ][ 2 ] = vgt_corners[ 1 ][ 0 ][ 1 ][ 2 ]
        = vgt_corners[ 0 ][ 1 ][ 1 ][ 2 ] = vgt_corners[ 1 ][ 1 ][ 1 ][ 2 ]
        = vgt_corners[ 0 ][ 0 ][ 0 ][ 2 ] + vgt_inc[ 2 ];
      for ( var vgt_jjj = 0; vgt_jjj < this.ncell[ 1 ]; vgt_jjj ++ ) {
        vgt_corners[ 0 ][ 0 ][ 0 ][ 1 ] = vgt_corners[ 1 ][ 0 ][ 0 ][ 1 ]
          = vgt_corners[ 0 ][ 0 ][ 1 ][ 1 ] = vgt_corners[ 1 ][ 0 ][ 1 ][ 1 ]
          = this.user_low[ 1 ] + vgt_inc[ 1 ] * vgt_jjj;
        vgt_corners[ 0 ][ 1 ][ 0 ][ 1 ] = vgt_corners[ 1 ][ 1 ][ 0 ][ 1 ]
          = vgt_corners[ 0 ][ 1 ][ 1 ][ 1 ] = vgt_corners[ 1 ][ 1 ][ 1 ][ 1 ]
          = vgt_corners[ 0 ][ 0 ][ 0 ][ 1 ] + vgt_inc[ 1 ];
        for ( var vgt_iii = 0; vgt_iii < this.ncell[ 0 ]; vgt_iii ++ ) {
          vgt_corners[ 0 ][ 0 ][ 0 ][ 0 ] = vgt_corners[ 0 ][ 1 ][ 0 ][ 0 ]
            = vgt_corners[ 0 ][ 0 ][ 1 ][ 0 ]
            = vgt_corners[ 0 ][ 1 ][ 1 ][ 0 ]
            = this.user_low[ 0 ] + vgt_inc[ 0 ] * vgt_iii;
          vgt_corners[ 1 ][ 0 ][ 0 ][ 0 ] = vgt_corners[ 1 ][ 1 ][ 0 ][ 0 ]
            = vgt_corners[ 1 ][ 0 ][ 1 ][ 0 ]
            = vgt_corners[ 1 ][ 1 ][ 1 ][ 0 ]
            = vgt_corners[ 0 ][ 0 ][ 0 ][ 0 ] + vgt_inc[ 0 ];
//        compute function values at cell vgt_corners
          for ( var vgt_h2 = 0; vgt_h2 < 2; vgt_h2 ++ ) {
            for ( var vgt_h1 = 0; vgt_h1 < 2; vgt_h1 ++ ) {
              for ( var vgt_h0 = 0; vgt_h0 < 2; vgt_h0 ++ ) {
                vgt_values[ vgt_h0 ][ vgt_h1 ][ vgt_h2 ] =
                  this.f3d( vgt_corners[ vgt_h0 ][ vgt_h1 ][ vgt_h2 ][ 0 ],
                            vgt_corners[ vgt_h0 ][ vgt_h1 ][ vgt_h2 ][ 1 ],
                      vgt_corners[ vgt_h0 ][ vgt_h1 ][ vgt_h2 ][ 2 ] ) - fv;
              }
            }
          }
//        compute vgt_intersections of iso-surface with cell edges
          for ( var vgt_d = 0; vgt_d < 3; vgt_d ++ ) {
            var vgt_d1 = ( vgt_d + 1 ) %3;
            var vgt_d2 = ( vgt_d + 2 ) %3;
            for ( var vgt_h1 = 0; vgt_h1 < 2; vgt_h1 ++ ) {
              vgt_hl[ vgt_d1 ] = vgt_hr[ vgt_d1 ] = vgt_h1;
              for ( var vgt_h2 = 0; vgt_h2 < 2; vgt_h2 ++ ) {
                vgt_hl[ vgt_d2 ] = vgt_hr[ vgt_d2 ] = vgt_h2;
                vgt_hl[ vgt_d ] = 0;
                vgt_hr[ vgt_d ] = 1;
                var vgt_vl =
                  vgt_values[ vgt_hl[ 0 ] ][ vgt_hl[ 1 ] ][ vgt_hl[ 2 ] ]; 
                var vgt_vr =
                  vgt_values[ vgt_hr[ 0 ] ][ vgt_hr[ 1 ] ][ vgt_hr[ 2 ] ]; 
                if ( vgt_vl * vgt_vr <= 0 ) {
                  var vgt_alpha = ( Math.abs( vgt_vl - vgt_vr ) > 0 ?
                    vgt_vl / ( vgt_vl - vgt_vr ): 0.5 );
                  for ( var vgt_dd = 0; vgt_dd < 3; vgt_dd ++ ) {
                    vgt_intersections[ vgt_d ][ vgt_h1 ][ vgt_h2 ][ vgt_dd ]
                    = vgt_corners[ vgt_hl[ 0 ] ][ vgt_hl[ 1 ] ]
                      [ vgt_hl[ 2 ] ][ vgt_dd ] * ( 1 - vgt_alpha )
                    + vgt_corners[ vgt_hr[ 0 ] ][ vgt_hr[ 1 ] ]
                      [ vgt_hr[ 2 ] ][ vgt_dd ] * vgt_alpha;
                  }
                }
              }
            }
          }
          var vgt_index = VolGraphTool.cubeIndex( vgt_values );
          var vgt_n = marching_cubes_table[ vgt_index ][ 1 ];
          var vgt_entry = new Array( );
          if ( vgt_n > 0 ) {
            vgt_entry = marching_cubes_table[ vgt_index ];
          } else if ( vgt_n < 0 ) {
            var vgt_index0 =
              VolGraphTool.subTableIndex( vgt_index, vgt_values );
            vgt_entry = marching_cubes_sub_table[ vgt_index0 ];
          }
          if ( vgt_n != 0 ) {
            var vgt_current = 6;
            var vgt_k = 1;
            var vgt_i = vgt_current;
            var vgt_size = vgt_entry[ vgt_k ];
            while ( vgt_size > 0 ) {
              var vgt_vertices = new Array( vgt_size );
              for ( var vgt_m = 0; vgt_m < vgt_size; vgt_m ++ ) {
                vgt_vertices[ vgt_m ] = new Array( 3 );
                vgt_n = vgt_entry[ vgt_i++ ];
                var vgt_h2 = vgt_n % 2;
                var vgt_nn = ( vgt_n - vgt_h2 ) / 2;
                var vgt_h1 = vgt_nn % 2;
                var vgt_dd = ( vgt_nn - vgt_h1 ) / 2;
                for ( var vgt_d = 0; vgt_d < 3; vgt_d ++ ) {
                  vgt_vertices[ vgt_m ][ vgt_d ] =
                    vgt_intersections[ vgt_dd ][ vgt_h1 ][ vgt_h2 ][ vgt_d ];
                }
              }
              var vgt_centroid = [ 0, 0, 0 ];
              for ( var vgt_m0 = 0; vgt_m0 < vgt_size; vgt_m0 ++ ) {
                for ( var vgt_d = 0; vgt_d < 3; vgt_d ++ ) {
                  vgt_centroid[ vgt_d ] += vgt_vertices[ vgt_m0 ][ vgt_d ];
                }
              }
              for ( var vgt_d = 0; vgt_d < 3; vgt_d ++ ) {
                vgt_centroid[ vgt_d ] /= vgt_size;
              }
              for ( var vgt_m = 0; vgt_m < vgt_size - 1; vgt_m ++ ) {
                this.colorTriangle( vgt_centroid, vgt_vertices[ vgt_m ],
                  vgt_vertices[ vgt_m + 1 ], vgt_ratio, vgt_ratio,
                  vgt_ratio );
              }
              this.colorTriangle( vgt_centroid, vgt_vertices[ vgt_size - 1 ],
                vgt_vertices[ 0 ], vgt_ratio, vgt_ratio, vgt_ratio );
              vgt_current += vgt_size;
              vgt_size = vgt_entry[ ++vgt_k ];
            }
          }
        }
      }
    }
  this.endDrawing();
//document.getElementById("debug_textarea").value +=
//  "leaving drawIsoSurface\n";
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.cubeIndex = function( v ) {
  var vgt_m = 1;
  var vgt_cube_index = 0;
  for ( var vgt_i0 = 0; vgt_i0 < 2; vgt_i0 ++ ) {
    for ( var vgt_i1 = 0; vgt_i1 < 2; vgt_i1 ++ ) {
      for ( var vgt_i2 = 0; vgt_i2 < 2; vgt_i2 ++ ) {
        if ( v[ vgt_i0 ][ vgt_i1 ][ vgt_i2 ] > 0 ) vgt_cube_index += vgt_m;
        vgt_m *= 2;
      }
    }
  }
  return vgt_cube_index;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.subTableIndex = function( index, values ) {
  var vgt_indx = 0;
  var vgt_indx0 = marching_cubes_table[ index ][ 2 ];
  var vgt_m = - marching_cubes_table[ index ][ 1 ];
  var vgt_val = [
    values[ 0 ][ 0 ][ 0 ],
    values[ 0 ][ 0 ][ 1 ],
    values[ 0 ][ 1 ][ 0 ],
    values[ 0 ][ 1 ][ 1 ],
    values[ 1 ][ 0 ][ 0 ],
    values[ 1 ][ 0 ][ 1 ],
    values[ 1 ][ 1 ][ 0 ],
    values[ 1 ][ 1 ][ 1 ] 
  ];
  switch ( vgt_m ) {
    case 1: {
      var vgt_k0 = marching_cubes_table[ index ][ 5 ];
      var vgt_k1 = marching_cubes_table[ index ][ 3 ];
      var vgt_k2 = marching_cubes_table[ index ][ 4 ];
      var vgt_k3 = ( vgt_k1 +vgt_k2 ) - vgt_k0;
      if ( ( vgt_val[ vgt_k0 ] * vgt_val[ vgt_k3 ] )
      > ( vgt_val[ vgt_k1 ] * vgt_val[ vgt_k2 ] ) ) {
        vgt_indx = 2 * vgt_indx0;
      } else {
        vgt_indx = 2 * vgt_indx0 + 1;
      }
      break;
    }
    case 2: {
      var vgt_s1 = 0;
      var vgt_s2 = 0;
      var vgt_k1 = marching_cubes_table[ index ][ 3 ];
      var vgt_k2 = marching_cubes_table[ index ][ 4 ];
      var vgt_k0 = marching_cubes_table[ index ][ 5 ];
      var vgt_k3 = ( vgt_k1 + vgt_k2 ) - vgt_k0;
      if ( ( vgt_val[ vgt_k0 ] * vgt_val[ vgt_k3 ] ) >
      ( vgt_val[ vgt_k1 ] * vgt_val[ vgt_k2 ] ) ) {
        vgt_s1 = 0;
      } else {
        vgt_s1 = 1;
      }
      vgt_k1 = marching_cubes_table[ index ][ 6 ];
      vgt_k2 = marching_cubes_table[ index ][ 7 ];
      vgt_k0 = marching_cubes_table[ index ][ 8 ];
      vgt_k3 = ( vgt_k1 + vgt_k2 ) - vgt_k0;
      if ( ( vgt_val[ vgt_k0 ] * vgt_val[ vgt_k3 ] ) >
      ( vgt_val[ vgt_k1 ] * vgt_val[ vgt_k2 ] ) ) {
        vgt_s2 = 0;
      } else {
        vgt_s2 = 1;
      }
      vgt_indx = 2 * vgt_indx0 + vgt_s1 * 2 + vgt_s2;
      break;
    }
    case 3: {
      var vgt_n = 3;
      var vgt_s = new Array( 3 );
      for ( var vgt_i = 0; vgt_i < 3; vgt_i++ ) {
        var vgt_k1 = marching_cubes_table[ index ][ vgt_n++ ];
        var vgt_k2 = marching_cubes_table[ index ][ vgt_n++ ];
        var vgt_k0 = marching_cubes_table[ index ][ vgt_n++ ];
        var vgt_k3 = ( vgt_k1 + vgt_k2 ) - vgt_k0;
        if( ( vgt_val[ vgt_k0 ] * vgt_val[ vgt_k3 ] ) >
        ( vgt_val[ vgt_k1 ] * vgt_val[ vgt_k2 ] ) ) {
          vgt_s[ vgt_i ] = 0;
        } else {
          vgt_s[ vgt_i ] = 1;
        }
      }
      vgt_indx = 2 * vgt_indx0 + vgt_s[ 0 ] * 4 + vgt_s[ 1 ] * 2
        + vgt_s[ 2 ];
      break;
    }
    case 6: {
      var vgt_n = 3;
      var vgt_s = new Array( 6 );
      for ( var vgt_i = 0; vgt_i < 6; vgt_i++ ) {
        var vgt_k1 = marching_cubes_table[ index ][ vgt_n++ ];
        var vgt_k2 = marching_cubes_table[ index ][ vgt_n++ ];
        var vgt_k0 = marching_cubes_table[ index ][ vgt_n++ ];
        var vgt_k3 = ( vgt_k1 + vgt_k2 ) - vgt_k0;
        if ( ( vgt_val[ vgt_k0 ] * vgt_val[ vgt_k3 ] ) >
        ( vgt_val[ vgt_k1 ] * vgt_val[ vgt_k2 ] ) ) {
          vgt_s[ vgt_i ] = 0;
        } else {
          vgt_s[ vgt_i ] = 1;
        }
      }
      vgt_indx = 2 * vgt_indx0 + vgt_s[ 0 ] * 32 + vgt_s[ 1 ] * 16
        + vgt_s[ 2 ] * 8 + vgt_s[ 3 ] * 4 + vgt_s[ 4 ] * 2 + vgt_s[ 5 ];
      break;
    }
    default: {
      alert("Error in SubTableIndex()!");
    }
  }
  return vgt_indx;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.draw3DClipPlanes = function( f, l, h, nc ) {
//document.getElementById("debug_textarea").value +=
//  "entering draw3DClipPlanes\n";
  this.drawing_clip_planes = true;
  this.use_lighting = true;
  this.rescale( l, h );
  this.current_quat = [ .5, .5, .5, .5 ]; // z up, y to right
  quat4.toMat4( this.current_quat, this.mv_matrix );
  this.clip_low =
    [ this.user_low[ 0 ], this.user_low[ 1 ], this.user_low[ 2 ] ];
  this.clip_high =
    [ this.user_high[ 0 ], this.user_high[ 1 ], this.user_high[ 2 ] ];
//document.getElementById("debug_textarea").value +=
//  "clip_low = " + this.clip_low[ 0 ] + " " + this.clip_low[ 1 ] + " "
//  + this.clip_low[ 2 ] + "\n";
//document.getElementById("debug_textarea").value +=
//  "clip_high = " + this.clip_high[ 0 ] + " " + this.clip_high[ 1 ] + " "
//  + this.clip_high[ 2 ] + "\n";
  this.makeOriginBox();
  this.makeBoundingBox();
  this.f3d = f;
  this.ncell = [ nc[ 0 ], nc[ 1 ], nc[ 2 ] ];
  this.fmin = Number.MAX_VALUE;
  this.fmax = - Number.MAX_VALUE;
  var vgt_inc = new Array( 3 );
  for ( var vgt_i = 0; vgt_i < 3; vgt_i ++ ) {
    vgt_inc[ vgt_i ] = this.user_len[ vgt_i ] / this.ncell[ vgt_i ];
  }
  for ( var vgt_k = 0; vgt_k <= nc[ 2 ]; vgt_k ++ ) {
    var vgt_z = this.user_low[ 2 ] + vgt_inc[ 2 ] * vgt_k;
    for ( var vgt_j = 0; vgt_j <= nc[ 1 ]; vgt_j ++ ) {
      var vgt_y = this.user_low[ 1 ] + vgt_inc[ 1 ] * vgt_j;
      for ( var vgt_i = 0; vgt_i <= nc[ 0 ]; vgt_i ++ ) {
        var vgt_x = this.user_low[ 0 ] + vgt_inc[ 0 ] * vgt_i;
        var vgt_fval = f( vgt_x, vgt_y, vgt_z );
        this.fmin = Math.min( this.fmin, vgt_fval );
        this.fmax = Math.max( this.fmax, vgt_fval );
      }
    }
  }
  if ( this.fmin >= this.fmax ) {
    this.fmax = Math.abs( this.fmax );
    this.fmin = - this.fmax;
  }
  if ( this.fmin >= this.fmax ) {
    this.fmax = 1;
    this.fmin = 0;
  }
  for ( var vgt_dd = 0; vgt_dd < 3; vgt_dd ++ ) {
    for ( var fgt_h = 0; fgt_h < 2; fgt_h ++ ) {
      var vgt_clip_val = ( fgt_h == 0 ? this.clip_low[ vgt_dd ] :
        this.clip_high[ vgt_dd ] );
      this.drawClipPlane( vgt_dd, fgt_h, vgt_clip_val );
    }
  }
  this.newPage();
  this.expose();
  this.watchMouse();
//document.getElementById("debug_textarea").value +=
//  "leaving draw3DClipPlanes\n";
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.drawClipPlane = function( d, h, cv ) {
//document.getElementById("debug_textarea").value +=
//  "entering drawClipPlane\n";
  if (! this.drawing_clip_planes ) return;
  if ( this.f3d == undefined ) return;
//document.getElementById("debug_textarea").value +=
//  "d,h,cv = " + d + " " + h + " " + cv + "\n";
//document.getElementById("debug_textarea").value +=
//  "fmin,fmax = " + this.fmin + " " + this.fmax + "\n";
//document.getElementById("debug_textarea").value +=
//  "user_low = " + this.user_low[ 0 ] + " " + this.user_low[ 1 ] + " "
//  + this.user_low[ 2 ] + "\n";
//document.getElementById("debug_textarea").value +=
//  "user_high = " + this.user_high[ 0 ] + " " + this.user_high[ 1 ] + " "
//  + this.user_high[ 2 ] + "\n";

  var vgt_d1 = ( d + 1 ) % 3;
  var vgt_d2 = ( d + 2 ) % 3;
  var vgt_low = [ this.user_low[ 0 ], this.user_low[ 1 ], this.user_low[ 2 ] ];
  var vgt_high = [ this.user_high[ 0 ], this.user_high[ 1 ],
    this.user_high[ 2 ] ];
  vgt_low[ d ] = cv;
  vgt_high[ d ] = cv;
  var vgt_ratios = new Array( 4 );
  var vgt_corners = new Array( 4 );
  for ( var vgt_k = 0; vgt_k < 4; vgt_k ++ ) {
    vgt_corners[ vgt_k ] = new Array( 3 );
    vgt_corners[ vgt_k ][ d ] = cv;
  }
  var vgt_inc = new Array( 3 );
  for ( var vgt_i = 0; vgt_i < 3; vgt_i ++ ) {
    vgt_inc[ vgt_i ] = this.user_len[ vgt_i ] / this.ncell[ vgt_i ];
  }
//document.getElementById("debug_textarea").value +=
//  "vgt_inc = " + vgt_inc[ 0 ] + " " + vgt_inc[ 1 ] + " " + vgt_inc[ 2 ] + "\n";
  this.beginDrawing();
    for ( var vgt_i2 = 0; vgt_i2 < this.ncell[ vgt_d2 ]; vgt_i2 ++ ) {
      vgt_corners[ 0 ][ vgt_d2 ] = vgt_corners[ 1 ][ vgt_d2 ] =
        this.user_low[ vgt_d2 ] + vgt_inc[ vgt_d2 ] * vgt_i2; 
      vgt_corners[ 2 ][ vgt_d2 ] = vgt_corners[ 3 ][ vgt_d2 ] =
        vgt_corners[ 0 ][ vgt_d2 ] + vgt_inc[ vgt_d2 ];
      for ( var vgt_i1 = 0; vgt_i1 < this.ncell[ vgt_d1 ]; vgt_i1 ++ ) {
        vgt_corners[ 0 ][ vgt_d1 ] = vgt_corners[ 3 ][ vgt_d1 ] =
          this.user_low[ vgt_d1 ] + vgt_inc[ vgt_d1 ] * vgt_i1; 
        vgt_corners[ 1 ][ vgt_d1 ] = vgt_corners[ 2 ][ vgt_d1 ] =
          vgt_corners[ 0 ][ vgt_d1 ] + vgt_inc[ vgt_d1 ];
        for ( var vgt_k = 0; vgt_k < 4; vgt_k ++ ) {
          var fv = this.f3d( vgt_corners[ vgt_k ][ 0 ],
            vgt_corners[ vgt_k ][ 1],
            vgt_corners[ vgt_k ][ 2 ] );
          vgt_ratios[ vgt_k ] = ( fv - this.fmin )
            / ( this.fmax - this.fmin );
        }
        if ( h == 0 ) {
          this.colorQuadrilateral(
            vgt_corners[ 0 ], vgt_corners[ 1 ], vgt_corners[ 2 ],
            vgt_corners[ 3 ], vgt_ratios[ 0 ], vgt_ratios[ 1 ],
            vgt_ratios[ 2 ], vgt_ratios[ 3 ] );
        } else {
          this.colorQuadrilateral(
            vgt_corners[ 3 ], vgt_corners[ 2 ], vgt_corners[ 1 ],
            vgt_corners[ 0 ], vgt_ratios[ 3 ], vgt_ratios[ 2 ],
            vgt_ratios[ 1 ], vgt_ratios[ 0 ] );
        }
      }
    }
  this.endDrawing();
//document.getElementById("debug_textarea").value +=
//  "leaving drawClipPlane\n";
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.rescale = function( l, h ) {
  for ( var vgt_d=0; vgt_d < 3; vgt_d ++ ) {
    this.user_low[ vgt_d ] = Math.min( l[ vgt_d ], h[ vgt_d ] );
    this.user_high[ vgt_d ] = Math.max( l[ vgt_d ], h[ vgt_d ] );
    if ( this.user_high[ 0 ] <= this.user_low[ 0 ] ) {
      this.user_high[ 0 ] += Math.abs( this.user_high[ 0 ] );
      if ( this.user_high[ 0 ] <= this.user_low[ 0 ] ) {
        this.user_high[ 0 ] += 1;
      }
    }
    this.user_center[ vgt_d ] =
      ( this.user_low[ vgt_d ] + this.user_high[ vgt_d ] ) * 0.5;
    this.user_len[ vgt_d ] =
      this.user_high[ vgt_d ] - this.user_low[ vgt_d ];
  }
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.getShader = function( gl, id ) {
  var vgt_shader_script = document.getElementById( id );
  if ( ! vgt_shader_script ) return null;
  var vgt_str = "";
  var vgt_k = vgt_shader_script.firstChild;
  while ( vgt_k ) {
    if ( vgt_k.nodeType == 3 ) vgt_str += vgt_k.textContent;
    vgt_k = vgt_k.nextSibling;
  }
  var vgt_shader;
  if ( vgt_shader_script.type == "x-shader/x-fragment" ) {
    vgt_shader = gl.createShader(gl.FRAGMENT_SHADER);
  } else if ( vgt_shader_script.type == "x-shader/x-vertex" ) {
    vgt_shader = gl.createShader(gl.VERTEX_SHADER);
  } else return null;
  gl.shaderSource( vgt_shader, vgt_str );
  gl.compileShader( vgt_shader );
  if ( !gl.getShaderParameter( vgt_shader, gl.COMPILE_STATUS ) ) {
    alert( gl.getShaderInfoLog( vgt_shader ) );
    return null;
  }
  return vgt_shader;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.userCoordsToObjectCoords = function( user_coords ) {
// object coords lie in (-1,1)x(-1,1)x(-1,1)
  return [
    2 * ( ( user_coords[ 0 ] - this.user_center[ 0 ] )
      / this.user_len[ 0 ] ),
    2 * ( ( user_coords[ 1 ] - this.user_center[ 1 ] )
      / this.user_len[ 1 ] ),
    2 * ( ( user_coords[ 2 ] - this.user_center[ 2 ] )
      / this.user_len[ 2 ] )
  ];
}
VolGraphTool.prototype.objectCoordsToUserCoords = function( obj_coords ) {
  return [
    this.user_center[ 0 ]
    + obj_coords[ 0 ] * this.user_len[ 0 ] * 0.5,
    this.user_center[ 1 ]
    + obj_coords[ 1 ] * this.user_len[ 1 ] * 0.5,
    this.user_center[ 2 ]
    + obj_coords[ 2 ] * this.user_len[ 2 ] * 0.5
  ];
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.watchMouse = function() {
//document.getElementById("debug_textarea").value +=
//  "entering VolGraphTool.watchMouse\n";
//see p 785 for mouse events
  this.canvas.my_graphtool = this;
  this.canvas.addEventListener( "click", VolGraphTool.mouseHandler,
    false );
  this.canvas.addEventListener( "contextmenu", VolGraphTool.mouseHandler,
    false );
  this.canvas.addEventListener( "mousedown", VolGraphTool.mouseHandler,
    false );
  this.canvas.addEventListener( "mousemove", VolGraphTool.mouseHandler,
    false );
  this.canvas.addEventListener( "mouseup", VolGraphTool.mouseHandler,
    false );
//document.getElementById("debug_textarea").value +=
//  "leaving VolGraphTool.watchMouse\n";
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.beginDrawing = function() {
//document.getElementById("debug_textarea").value +=
//  "entering VolGraphTool.beginDrawing\n";
  this.current_vertices = new Array();
  this.current_normals = new Array(); // new
  this.current_colors = new Array();
  this.current_vertex_indices = new Array();
  this.current_vertex_number = 0;
//document.getElementById("debug_textarea").value +=
//  "leaving VolGraphTool.beginDrawing\n";
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.colorTriangle = function( v0, v1, v2, r0, r1, r2 ) {
//document.getElementById("debug_textarea").value +=
//  "entering VolGraphTool.colorTriangle\n";
//document.getElementById("debug_textarea").value +=
//  "v0 = " + v0[ 0 ] + " " + v0[ 1 ] + " " + v0[ 2 ] + "\n";
//document.getElementById("debug_textarea").value +=
//  "v1 = " + v1[ 0 ] + " " + v1[ 1 ] + " " + v1[ 2 ] + "\n";
//document.getElementById("debug_textarea").value +=
//  "v2 = " + v2[ 0 ] + " " + v2[ 1 ] + " " + v2[ 2 ] + "\n";
//document.getElementById("debug_textarea").value +=
//  "r0,r1,r2 = " + r0 + " " + r1 + " " + r2 + "\n";
//v0 -> v1 -> v2 should be oriented by right-hand rule with normal
  var vgt_oc0 = this.userCoordsToObjectCoords( v0 );
  this.current_vertices.push( vgt_oc0[ 0 ], vgt_oc0[ 1 ], vgt_oc0[ 2 ] );
  var vgt_oc1 = this.userCoordsToObjectCoords( v1 );
  this.current_vertices.push( vgt_oc1[ 0 ], vgt_oc1[ 1 ], vgt_oc1[ 2 ] );
  var vgt_oc2 = this.userCoordsToObjectCoords( v2 );
  this.current_vertices.push( vgt_oc2[ 0 ], vgt_oc2[ 1 ], vgt_oc2[ 2 ] );

  var vgt_normal = VolGraphTool.vertexNormal( vgt_oc0, vgt_oc1, vgt_oc2 );
  this.current_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ],
    vgt_normal[ 2 ] );
  this.current_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ],
    vgt_normal[ 2 ] );
  this.current_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ],
    vgt_normal[ 2 ] );

  var vgt_c0=this.colormap.getGLColorFromRatio( r0 );
  var vgt_c1=this.colormap.getGLColorFromRatio( r1 );
  var vgt_c2=this.colormap.getGLColorFromRatio( r2 );
  this.current_colors.push( vgt_c0[ 0 ], vgt_c0[ 1 ], vgt_c0[ 2 ], 1);
  this.current_colors.push( vgt_c1[ 0 ], vgt_c1[ 1 ], vgt_c1[ 2 ], 1);
  this.current_colors.push( vgt_c2[ 0 ], vgt_c2[ 1 ], vgt_c2[ 2 ], 1);

  var cvn = this.current_vertex_number;
  this.current_vertex_indices.push( cvn, cvn + 1, cvn + 2);
  this.current_vertex_number += 3;
//document.getElementById("debug_textarea").value +=
//  "leaving VolGraphTool.colorTriangle\n";
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.vertexNormal = function( oc0, oc1, oc2 ) {
  var vgt_oc10 = new Array( 3 );
  vec3.subtract( oc1, oc0, vgt_oc10 );
  var vgt_oc20 = new Array( 3 );
  vec3.subtract( oc2, oc0, vgt_oc20 );
  var vgt_cp = new Array( 3 );
  vec3.cross( vgt_oc10, vgt_oc20, vgt_cp );
  vec3.normalize( vgt_cp );
  return vgt_cp;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.colorQuadrilateral = function( v0, v1, v2, v3,
r0, r1, r2, r3 ) {
//document.getElementById("debug_textarea").value +=
//  "entering VolGraphTool.colorQuadrilateral\n";
//
// v0 -> v1 -> v2 -> v3 should be oriented by right-hand rule with normal
  var vgt_center = new Array( 3 );
  for ( var vgt_d = 0; vgt_d < 3; vgt_d ++ ) {
    vgt_center[ vgt_d ] = ( v0[ vgt_d ] + v1[ vgt_d ] + v2[ vgt_d ]
      + v3[ vgt_d ] ) * 0.25;
  }
  var vgt_r = ( r0 + r1 + r2 + r3 ) * 0.25;
//document.getElementById("debug_textarea").value +=
//  "r0,r1,r2,r3,vgt_r = " + r0 + " " + r1 + " " + r2 + " " + r3 + " "
//  + vgt_r + "\n";

  var vgt_oc0 = this.userCoordsToObjectCoords( v0 );
  this.current_vertices.push( vgt_oc0[ 0 ], vgt_oc0[ 1 ], vgt_oc0[ 2 ] );
  var vgt_oc1 = this.userCoordsToObjectCoords( v1 );
  this.current_vertices.push( vgt_oc1[ 0 ], vgt_oc1[ 1 ], vgt_oc1[ 2 ] );
  var vgt_oc2 = this.userCoordsToObjectCoords( v2 );
  this.current_vertices.push( vgt_oc2[ 0 ], vgt_oc2[ 1 ], vgt_oc2[ 2 ] );
  var vgt_oc3 = this.userCoordsToObjectCoords( v3 );
  this.current_vertices.push( vgt_oc3[ 0 ], vgt_oc3[ 1 ], vgt_oc3[ 2 ] );
  var vgt_oc = this.userCoordsToObjectCoords( vgt_center );
  this.current_vertices.push( vgt_oc[ 0 ], vgt_oc[ 1 ], vgt_oc[ 2 ] );

  var vgt_normal = VolGraphTool.vertexNormal( vgt_oc0, vgt_oc1, vgt_oc3 );
  this.current_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ],
    vgt_normal[ 2 ] );
  vgt_normal = VolGraphTool.vertexNormal( vgt_oc1, vgt_oc2, vgt_oc0 );
  this.current_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ],
    vgt_normal[ 2 ] );
  vgt_normal = VolGraphTool.vertexNormal( vgt_oc2, vgt_oc3, vgt_oc1 );
  this.current_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ],
    vgt_normal[ 2 ] );
  vgt_normal = VolGraphTool.vertexNormal( vgt_oc3, vgt_oc0, vgt_oc2 );
  this.current_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ],
    vgt_normal[ 2 ] );
  var vgt_oc20 = new Array( 3 );
  vec3.subtract( vgt_oc2, vgt_oc0, vgt_oc20 );
  var vgt_oc31 = new Array( 3 );
  vec3.subtract( vgt_oc3, vgt_oc1, vgt_oc31 );
  vgt_normal = new Array( 3 );
  vec3.cross( vgt_oc20, vgt_oc31, vgt_normal );
  vec3.normalize( vgt_normal );
  this.current_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ],
    vgt_normal[ 2 ] );

  var vgt_c0=this.colormap.getGLColorFromRatio( r0 );
  var vgt_c1=this.colormap.getGLColorFromRatio( r1 );
  var vgt_c2=this.colormap.getGLColorFromRatio( r2 );
  var vgt_c3=this.colormap.getGLColorFromRatio( r3 );
  var vgt_c=this.colormap.getGLColorFromRatio( vgt_r );
  this.current_colors.push( vgt_c0[ 0 ], vgt_c0[ 1 ], vgt_c0[ 2 ], 1 );
  this.current_colors.push( vgt_c1[ 0 ], vgt_c1[ 1 ], vgt_c1[ 2 ], 1 );
  this.current_colors.push( vgt_c2[ 0 ], vgt_c2[ 1 ], vgt_c2[ 2 ], 1 );
  this.current_colors.push( vgt_c3[ 0 ], vgt_c3[ 1 ], vgt_c3[ 2 ], 1 );
  this.current_colors.push( vgt_c[ 0 ], vgt_c[ 1 ], vgt_c[ 2 ], 1 );

  var vgt_cvn = this.current_vertex_number;
  this.current_vertex_indices.push(
    vgt_cvn    , vgt_cvn + 1, vgt_cvn + 4,
    vgt_cvn + 1, vgt_cvn + 2, vgt_cvn + 4,
    vgt_cvn + 2, vgt_cvn + 3, vgt_cvn + 4,
    vgt_cvn + 3, vgt_cvn    , vgt_cvn + 4
  );
  this.current_vertex_number += 5;
//document.getElementById("debug_textarea").value +=
//  "leaving VolGraphTool.colorQuadrilateral\n";
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.endDrawing = function() {
//document.getElementById("debug_textarea").value +=
//  "entering VolGraphTool.endDrawing\n";
  var vgt_vertex_position_buffer = this.gl.createBuffer();
  this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_position_buffer );
  this.gl.bufferData( this.gl.ARRAY_BUFFER,
    new Float32Array( this.current_vertices ), this.gl.STATIC_DRAW );
  vgt_vertex_position_buffer.my_numItems = this.current_vertices.length / 3;

//begin new:
  var vgt_vertex_normal_buffer = this.gl.createBuffer();
  this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_normal_buffer );
  this.gl.bufferData( this.gl.ARRAY_BUFFER,
    new Float32Array( this.current_normals ), this.gl.STATIC_DRAW );
  vgt_vertex_normal_buffer.my_numItems = this.current_normals.length / 3;
//end new

  var vgt_vertex_color_buffer = this.gl.createBuffer();
  this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_color_buffer );
  this.gl.bufferData( this.gl.ARRAY_BUFFER,
    new Float32Array( this.current_colors ), this.gl.STATIC_DRAW );
  vgt_vertex_color_buffer.my_numItems = this.current_colors.length / 4;

  var vgt_vertex_index_buffer = this.gl.createBuffer();
  this.gl.bindBuffer( this.gl.ELEMENT_ARRAY_BUFFER,
    vgt_vertex_index_buffer );
  this.gl.bufferData( this.gl.ELEMENT_ARRAY_BUFFER,
    new Uint16Array( this.current_vertex_indices ), this.gl.STATIC_DRAW );
  vgt_vertex_index_buffer.my_numItems = this.current_vertex_indices.length;

  this.surface_drawing_array.push( {
    vpb : vgt_vertex_position_buffer,
    vnb : vgt_vertex_normal_buffer,
    vcb : vgt_vertex_color_buffer,
    vib : vgt_vertex_index_buffer
  } );
//document.getElementById("debug_textarea").value +=
//  "leaving VolGraphTool.endDrawing\n";
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.newPage = function() {
//document.getElementById("debug_textarea").value +=
//  "entering VolGraphTool.newPage\n";
  this.gl.viewport( 0, 0, this.canvas.width, this.canvas.height );
  this.gl.clear( this.gl.COLOR_BUFFER_BIT | this.gl.DEPTH_BUFFER_BIT );
//document.getElementById("debug_textarea").value +=
//  "leaving VolGraphTool.newPage\n";
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.expose = function() {
//document.getElementById("debug_textarea").value +=
//  "entering VolGraphTool.expose\n";
//document.getElementById("debug_textarea").value +=
//  "use_lighting = " + this.use_lighting + "\n";

  var vgt_normal_matrix = mat4.create();
  mat4.inverse(this.mv_matrix, vgt_normal_matrix);
  mat4.transpose( vgt_normal_matrix );

  var vgt_adjustedLD = vec3.create();
  vec3.normalize(LIGHTING_DIRECTION, vgt_adjustedLD);
  vec3.scale(vgt_adjustedLD, -1);

  for ( var vgt_s = 0; vgt_s < this.surface_drawing_array.length; vgt_s ++ )
  {
//  document.getElementById("debug_textarea").value +=
//    "drawing surface = " + vgt_s + "\n";
    var vgt_vertex_position_buffer = this.surface_drawing_array[ vgt_s ].vpb;
    this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_position_buffer );
    this.gl.vertexAttribPointer(
      this.programInfo.attribLocations.vertexPosition, 3, this.gl.FLOAT,
      false, 0, 0 );
    this.gl.enableVertexAttribArray(
      this.programInfo.attribLocations.vertexPosition);

    var vgt_vertex_normal_buffer = this.surface_drawing_array[ vgt_s ].vnb;
    this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_normal_buffer );
    this.gl.vertexAttribPointer(
      this.programInfo.attribLocations.vertexNormal, 3, this.gl.FLOAT,
      false, 0, 0);
    this.gl.enableVertexAttribArray(
      this.programInfo.attribLocations.vertexNormal);
    var vgt_vertex_color_buffer = this.surface_drawing_array[ vgt_s ].vcb;
    this.gl.bindBuffer( this.gl.ARRAY_BUFFER,
      vgt_vertex_color_buffer );
    this.gl.vertexAttribPointer(
      this.programInfo.attribLocations.vertexColor, 4, this.gl.FLOAT,
      false, 0, 0 );
    this.gl.enableVertexAttribArray(
      this.programInfo.attribLocations.vertexColor);

    var vgt_vertex_index_buffer = this.surface_drawing_array[ vgt_s ].vib;
    this.gl.bindBuffer( this.gl.ELEMENT_ARRAY_BUFFER,
      vgt_vertex_index_buffer );

    this.gl.useProgram( this.shader_program );

    this.gl.uniformMatrix4fv(
      this.programInfo.uniformLocations.projectionMatrix, false,
      this.p_matrix ); // should be p_matrix_clip
    this.gl.uniformMatrix4fv(
      this.programInfo.uniformLocations.modelViewMatrix, false,
      this.mv_matrix );
    this.gl.uniformMatrix4fv(
      this.programInfo.uniformLocations.normalMatrix, false,
      vgt_normal_matrix );
//
    if ( this.use_lighting ) {
      this.gl.uniform3f( this.programInfo.uniformLocations.ambientColor,
        AMBIENT_COLOR[ 0 ], AMBIENT_COLOR[ 1 ], AMBIENT_COLOR[ 2 ] );
      this.gl.uniform3fv(
        this.programInfo.uniformLocations.lightingDirection,
        vgt_adjustedLD);
      this.gl.uniform3f(
        this.programInfo.uniformLocations.directionalColor,
        DIRECTIONAL_COLOR[ 0 ], DIRECTIONAL_COLOR[ 1 ],
        DIRECTIONAL_COLOR[ 2 ] );
    } else {
      this.gl.uniform3f( this.programInfo.uniformLocations.ambientColor,
        1., 1., 1. );
      this.gl.uniform3fv(
        this.programInfo.uniformLocations.lightingDirection,
        vgt_adjustedLD);
      this.gl.uniform3f(
        this.programInfo.uniformLocations.directionalColor, 0., 0., 0. );
    }
//
    this.gl.drawElements( this.gl.TRIANGLES,
      vgt_vertex_index_buffer.my_numItems, this.gl.UNSIGNED_SHORT, 0);
  }
/*
  for ( var vgt_c = 0; vgt_c < this.curve_drawing_array.length; vgt_c ++ ) {
    document.getElementById("debug_textarea").value +=
      "drawing surface vgt_c = " + vgt_c + "\n";
    var vgt_vertex_position_buffer =
      this.curve_drawing_array[ vgt_c ].vpb;
    this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_position_buffer );
    this.gl.vertexAttribPointer(
      this.programInfo.attribLocations.vertexPosition, 3, this.gl.FLOAT,
      false, 0, 0 );

    var vgt_vertex_color_buffer = this.curve_drawing_arrray[ vgt_c ].vcb;
    this.gl.bindBuffer( this.gl.ARRAY_BUFFER,
      vgt_vertex_color_buffer );
    this.gl.vertexAttribPointer( this.shader_program.vertexColorAttribute,
      4, this.gl.FLOAT, false, 0, 0 );

    var vgt_vertex_index_buffer = this.curve_drawing_array[ vgt_c ].vib;
    this.gl.bindBuffer( this.gl.ELEMENT_ARRAY_BUFFER,
      vgt_vertex_index_buffer );
    this.gl.uniformMatrix4fv( this.shader_program.pMatrixUniform, false,
      this.p_matrix ); // should be p_matrix_clip
    this.gl.uniformMatrix4fv( this.shader_program.mvMatrixUniform, false,
      this.mv_matrix );
    gl.lineWidth( 1 );
    this.gl.drawElements( gl.LINES, vgt_vertex_index_buffer.my_numItems,
      gl.UNSIGNED_SHORT, 0);
  }
*/
  if ( ! this.button_pressed ) {
//  document.getElementById("debug_textarea").value +=
//    "leaving VolGraphTool.expose\n";
    return;
  }
  for ( var vgt_b = 0; vgt_b < this.bounding_box_array.length; vgt_b ++ ) {
//  document.getElementById("debug_textarea").value +=
//    "drawing bounding box face vgt_b = " + vgt_b + "\n";
    var vgt_vertex_position_buffer = this.bounding_box_array[ vgt_b ].vpb;
    this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_position_buffer );
    this.gl.vertexAttribPointer(
      this.programInfo.attribLocations.vertexPosition, 3, this.gl.FLOAT,
      false, 0, 0 );
    this.gl.enableVertexAttribArray(
      this.programInfo.attribLocations.vertexPosition);

/*
    var vgt_vertex_normal_buffer = this.bounding_box_array[ vgt_b ].vnb;
    this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_normal_buffer );
    this.gl.vertexAttribPointer(
      this.shader_program.vertexNormalAttribute, 3, this.gl.FLOAT, false,
      0, 0);

    if ( this.use_lighting ) {
      this.gl.uniform3f( this.shader_program.ambientColorUniform,
        0, 0, 0 );
      this.gl.uniform3fv( this.shader_program.lightingDirectionUniform,
        vgt_adjustedLD);
      this.gl.uniform3f( this.shader_program.directionalColorUniform,
        DIRECTIONAL_COLOR[ 0 ], DIRECTIONAL_COLOR[ 1 ],
        DIRECTIONAL_COLOR[ 2 ] );
    }
*/

    var vgt_vertex_color_buffer = this.bounding_box_array[ vgt_b ].vcb;
    this.gl.bindBuffer( this.gl.ARRAY_BUFFER,
      vgt_vertex_color_buffer );
    this.gl.vertexAttribPointer(
      this.programInfo.attribLocations.vertexColor, 4, this.gl.FLOAT,
      false, 0, 0 );
    this.gl.enableVertexAttribArray(
      this.programInfo.attribLocations.vertexColor);

    var vgt_vertex_index_buffer = this.bounding_box_array[ vgt_b ].vib;
    this.gl.bindBuffer( this.gl.ELEMENT_ARRAY_BUFFER,
      vgt_vertex_index_buffer );
    this.gl.uniformMatrix4fv(
      this.programInfo.uniformLocations.projectionMatrix, false,
      this.p_matrix ); // should be p_matrix_clip
    this.gl.uniformMatrix4fv(
      this.programInfo.uniformLocations.modelViewMatrix, false,
      this.mv_matrix );
    this.gl.uniformMatrix4fv(
      this.programInfo.uniformLocations.normalMatrix, false,
      vgt_normal_matrix );

    this.gl.lineWidth( 1 );
    this.gl.drawElements( this.gl.LINES, vgt_vertex_index_buffer.my_numItems,
      this.gl.UNSIGNED_SHORT, 0);
  }
  for ( var vgt_b = 0; vgt_b < this.clip_bounding_box_array.length;
  vgt_b ++ ) {
//  document.getElementById("debug_textarea").value +=
//    "drawing clip bounding box face vgt_b = " + vgt_b + "\n";
    var vgt_vertex_position_buffer =
      this.clip_bounding_box_array[ vgt_b ].vpb;
    this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_position_buffer );
    this.gl.vertexAttribPointer(
      this.programInfo.attribLocations.vertexPosition, 3, this.gl.FLOAT,
      false, 0, 0 );
    this.gl.enableVertexAttribArray(
      this.programInfo.attribLocations.vertexPosition);

/*
    var vgt_vertex_normal_buffer = this.bounding_box_array[ vgt_b ].vnb;
    this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_normal_buffer );
    this.gl.vertexAttribPointer(
      this.shader_program.vertexNormalAttribute, 3, this.gl.FLOAT, false,
      0, 0);

    if ( this.use_lighting ) {
      this.gl.uniform3f( this.shader_program.ambientColorUniform,
        0, 0, 0 );
      this.gl.uniform3fv( this.shader_program.lightingDirectionUniform,
        vgt_adjustedLD);
      this.gl.uniform3f( this.shader_program.directionalColorUniform,
        DIRECTIONAL_COLOR[ 0 ], DIRECTIONAL_COLOR[ 1 ],
        DIRECTIONAL_COLOR[ 2 ] );
    }
*/

    var vgt_vertex_color_buffer = this.clip_bounding_box_array[ vgt_b ].vcb;
    this.gl.bindBuffer( this.gl.ARRAY_BUFFER,
      vgt_vertex_color_buffer );
    this.gl.vertexAttribPointer(
      this.programInfo.attribLocations.vertexColor, 4, this.gl.FLOAT,
      false, 0, 0 );

    var vgt_vertex_index_buffer = this.clip_bounding_box_array[ vgt_b ].vib;
    this.gl.bindBuffer( this.gl.ELEMENT_ARRAY_BUFFER,
      vgt_vertex_index_buffer );
    this.gl.uniformMatrix4fv(
      this.programInfo.uniformLocations.projectionMatrix, false,
      this.p_matrix ); // should be p_matrix_clip
    this.gl.uniformMatrix4fv(
      this.programInfo.uniformLocations.modelViewMatrix, false,
      this.mv_matrix );
    this.gl.uniformMatrix4fv(
      this.programInfo.uniformLocations.normalMatrix, false,
      vgt_normal_matrix );

    this.gl.lineWidth( 3 );
    this.gl.drawElements( this.gl.LINES, vgt_vertex_index_buffer.my_numItems,
      this.gl.UNSIGNED_SHORT, 0);
  }
  for ( var vgt_h = 0; vgt_h < this.crosshair_array.length; vgt_h ++ ) {
//  document.getElementById("debug_textarea").value +=
//    "drawing cross hairs for vgt_h = " + vgt_h + "\n";
    var vgt_vertex_position_buffer = this.crosshair_array[ vgt_h ].vpb;
    this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_position_buffer );
    this.gl.vertexAttribPointer(
      this.programInfo.attribLocations.vertexPosition, 3, this.gl.FLOAT,
      false, 0, 0 );
    this.gl.enableVertexAttribArray(
      this.programInfo.attribLocations.vertexPosition);

    var vgt_vertex_color_buffer = this.crosshair_array[ vgt_h ].vcb;
    this.gl.bindBuffer( this.gl.ARRAY_BUFFER,
      vgt_vertex_color_buffer );
    this.gl.vertexAttribPointer(
      this.programInfo.attribLocations.vertexColor, 4, this.gl.FLOAT,
      false, 0, 0 );
    this.gl.enableVertexAttribArray(
      this.programInfo.attribLocations.vertexColor);

    var vgt_vertex_index_buffer = this.crosshair_array[ vgt_h ].vib;
    this.gl.bindBuffer( this.gl.ELEMENT_ARRAY_BUFFER,
      vgt_vertex_index_buffer );
    this.gl.uniformMatrix4fv(
      this.programInfo.uniformLocations.projectionMatrix, false,
      this.p_matrix ); // should be p_matrix_clip
    this.gl.uniformMatrix4fv(
      this.programInfo.uniformLocations.modelViewMatrix, false,
      this.mv_matrix );
    this.gl.lineWidth( 1 );
    this.gl.drawElements( this.gl.LINES, vgt_vertex_index_buffer.my_numItems,
      this.gl.UNSIGNED_SHORT, 0);
  }
//document.getElementById("debug_textarea").value +=
//  "leaving VolGraphTool.expose\n";
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
VolGraphTool.prototype.pushMVMatrix = function() {
  var copy = new glMatrixArrayType(16);
  for ( var vgt_i = 0; vgt_i < 16; vgt_i ++ ) {
    copy[ vgt_i ] = this.mv_matrix[ vgt_i ];
  }
  this.mv_matrix_stack.push(copy);
}
VolGraphTool.prototype.popMVMatrix = function() {
  if ( this.mv_matrix_stack.length == 0 ) throw "Invalid popMatrix!";
  this.mv_matrix = this.mv_matrix_stack.pop();
}
*/
/*
VolGraphTool.prototype.setLineWidth = function( lw ) { 
  this.context.lineWidth = lw;
}
VolGraphTool.prototype.setFgColor = function( ratio ) { 
  if ( ratio == undefined ) { // for text:
    this.context.strokeStyle = this.colormap.getHexColor( 1 );
  } else { // for color lines or fill:
    var color = 0;
    if (ratio < 0 ) color = 2;
    else if ( ratio > 1 ) color = this.colormap.num_cmap_colors - 1;
    else {
      color = 2
        + Math.round( ratio * ( this.colormap.num_cmap_colors - 3 ) );
    }
    var hc = this.colormap.getHexColor( color );
    this.context.fillStyle = hc;
    this.context.strokeStyle = hc;
  }
}
VolGraphTool.prototype.setFont = function( s ) { this.context.font = s; }
VolGraphTool.prototype.putString = function( user_pos, s, align_left_right,
align_up_down, vertical ) {
  this.context.save();
  if ( align_left_right != undefined ) {
    this.context.textAlign = align_left_right; // see Flanagan p. 651
  }
  if ( align_up_down != undefined ) {
    this.context.textBaseline = align_up_down;
  }
  var screen_coords = this.screenCoords( user_pos );
  if ( vertical ) {
    this.context.rotate( - Math.PI * 0.5 );
    this.context.fillText( s, - screen_coords[ 1 ], screen_coords[ 0 ] );
  } else {
    this.context.fillText( s, screen_coords[ 0 ], screen_coords[ 1 ] );
  }
  this.context.restore();
}
*/
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
VolGraphTool.handleContextLost = function( e ) {
}
VolGraphTool.handleContextRestored = function( ) {
//e.currentTarget.my_graphtool.init();
//e.currentTarget.my_graphtool.f();
}
*/
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.mouseDown = function( e ) {
//document.getElementById("debug_textarea").value +=
//  "entering VolGraphTool.mouseDown\n";
  this.old_screen_coords =
    [ this.new_screen_coords[ 0 ], this.new_screen_coords[ 1 ] ];
//document.getElementById("debug_textarea").value +=
//  "button = " + this.button + "\n";
//document.getElementById("debug_textarea").value +=
//  "old_screen_coords = " + this.old_screen_coords[ 0 ] + " "
//  + this.old_screen_coords[ 1 ] + "\n";
  switch ( this.button ) {
    case LEFT_BUTTON :
      this.old_unit_sphere_coords =
        this.screenCoordsToUnitSphereCoords(this.old_screen_coords );
      quat4.set( this.current_quat, this.old_quat );
      break;
    case MIDDLE_BUTTON :
      this.selectClipPlane();
      break;
    case RIGHT_BUTTON:
      this.markFaceLocation();
      break;
  }
//document.getElementById("debug_textarea").value +=
//  "leaving VolGraphTool.mouseDown\n";
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.mouseMove = function( e ) {
//document.getElementById("debug_textarea").value +=
//  "entering VolGraphTool.mouseMove\n";
//document.getElementById("debug_textarea").value +=
//  "button = " + this.button + "\n";
//document.getElementById("debug_textarea").value +=
//  "new_screen_coords = " + this.new_screen_coords[ 0 ] + " "
//  + this.new_screen_coords[ 1 ] + "\n";
  if ( Math.sqrt( 
    Math.pow( this.new_screen_coords[ 0 ]
            - this.old_screen_coords[ 0 ], 2 )
  + Math.pow( this.new_screen_coords[ 1 ]
            - this.old_screen_coords[ 1 ], 2 )) < 4 ) {
    return;
  }
  switch ( this.button ) {
    case LEFT_BUTTON :
      this.rotateImage();
      break;
    case MIDDLE_BUTTON :
      this.moveClipPlane();
      break;
    case RIGHT_BUTTON:
      this.moveCrossHairs();
      this.old_screen_coords =
        [ this.new_screen_coords[ 0 ], this.new_screen_coords[ 1 ] ];  
  }
//document.getElementById("debug_textarea").value +=
//  "leaving VolGraphTool.mouseMove\n";
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.mouseUp = function( e ) {
//document.getElementById("debug_textarea").value +=
//  "entering VolGraphTool.mouseUp\n";
//document.getElementById("debug_textarea").value +=
//  "button = " + this.button + "\n";
  switch ( this.button ) {
    case LEFT_BUTTON :
      break;
    case MIDDLE_BUTTON :
//    plot_obj->plot( active_clip_direction, this.active_clip_hand );
      break;
    case RIGHT_BUTTON:
      this.crosshair_array = new Array();
      break;
  }
//document.getElementById("debug_textarea").value +=
//  "leaving VolGraphTool.mouseUp\n";
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.eventToScreenCoords = function( e ) {
//http://www.quirksmode.org/js/events_properties.html
//http://miloq.blogspot.com/2011/05/coordinates-mouse-click-canvas.html
  var vgt_sc = new Array( 2 );
  if ( e.x != undefined && e.y != undefined ) { // chrome, microsoft
    if ( e.pageX != undefined && e.pageY != undefined ) { // chrome
      vgt_sc = [ e.pageX, e.pageY ];
    } else {
      vgt_sc = [ e.x, e.y ];
    }
  } else { // firefox : W3C
    vgt_sc = [ e.clientX + document.body.scrollLeft
                         + document.documentElement.scrollLeft,
               e.clientY + document.body.scrollTop
                         + document.documentElement.scrollTop ];
  }
  vgt_sc[ 0 ] -= e.currentTarget.offsetLeft;
  vgt_sc[ 1 ] -= e.currentTarget.offsetTop;
  return vgt_sc;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
//mouse position is read-only in JavaScript:
//http://stackoverflow.com/questions/1208729/jquery-set-mouse-position-not-cursor-position
VolGraphTool.screenCoordsToEvent = function( sc, e ) {
  var oc = [
    sc[ 0 ] + e.currentTarget.offsetLeft,
    sc[ 1 ] + e.currentTarget.offsetTop
  ];
  if ( e.x != undefined && e.y != undefined ) { // microsoft
    e.x = sc[ 0 ];
    e.y = sc[ 1 ];
  } else {
    e.clientX = sc[ 0 ] - document.body.scrollLeft
                        - document.documentElement.scrollLeft;
    e.clientY = sc[ 1 ] - document.body.scrollTop
                        - document.documentElement.scrollTop;
  }
}
*/
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.mouseHandler = function( e ) {
//document.getElementById("debug_textarea").value +=
//  "entering VolGraphTool.mouseHandler\n";
  e.preventDefault();
  var vgt_gt=e.currentTarget.my_graphtool;
  if ( e.type != "mousedown" && ! vgt_gt.button_pressed ) return;

  vgt_gt.new_screen_coords = new Array( 2 );
//http://www.quirksmode.org/js/events_properties.html
//http://miloq.blogspot.com/2011/05/coordinates-mouse-click-canvas.html
  if ( e.type == "mousedown" ) {
    if ( navigator.userAgent.indexOf("Chrome") != -1 ) {
      if ( e.button == 0 ) vgt_gt.button = LEFT_BUTTON;
      else if ( e.button == 2 ) vgt_gt.button = RIGHT_BUTTON;
      else vgt_gt.button = MIDDLE_BUTTON;
    } else if ( navigator.userAgent.indexOf("Firefox") != -1 ) {
      if ( e.button == 0 ) vgt_gt.button = LEFT_BUTTON;
      else if ( e.button == 2 ) vgt_gt.button = RIGHT_BUTTON;
      else vgt_gt.button = MIDDLE_BUTTON;
    } else if ( navigator.userAgent.indexOf("MSIE") != -1 ) {
// don't know if these are correct
      if ( e.button == 0 ) vgt_gt.button = LEFT_BUTTON;
      else if ( e.button == 2 ) vgt_gt.button = RIGHT_BUTTON;
      else vgt_gt.button = MIDDLE_BUTTON;
    } else if ( navigator.userAgent.indexOf("Opera") != -1 ) {
// don't know if these are correct
      if ( e.button == 0 ) vgt_gt.button = LEFT_BUTTON;
      else if ( e.button == 2 ) vgt_gt.button = RIGHT_BUTTON;
      else vgt_gt.button = MIDDLE_BUTTON;
    } else { // unknown browser
    }
  }
  vgt_gt.new_screen_coords = VolGraphTool.eventToScreenCoords( e ); 

  switch ( e.type ) {
    case "mousedown" :
      vgt_gt.button_pressed = true;
      vgt_gt.mouseDown( e );
      break;
    case "mousemove" :
      vgt_gt.mouseMove( e );
      break;
    case "mouseup" :
      vgt_gt.mouseUp( e );
      vgt_gt.button_pressed = false;
      break;
  }
  vgt_gt.newPage();
  vgt_gt.expose();
//document.getElementById("debug_textarea").value +=
//  "leaving VolGraphTool.mouseHandler\n";
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.makeOriginBox = function() {
//document.getElementById("debug_textarea").value +=
//  "entering VolGraphTool.makeOriginBox\n";
//put small box near origin
  var vgt_lo_lo = new Array( 3 );
  for ( var vgt_d = 0; vgt_d < 3; vgt_d ++ ) {
    var vgt_inc = ( this.user_high[ vgt_d ] - this.user_low[ vgt_d ] )
      * GRAPH_FUDGE;
    vgt_lo_lo[ vgt_d ] = this.user_low[ vgt_d ] - vgt_inc;
  }
  this.beginDrawing();
    this.colorQuadrilateral(
      [ vgt_lo_lo[ 0 ], vgt_lo_lo[ 1 ], vgt_lo_lo[ 2 ] ],
      [ vgt_lo_lo[ 0 ], vgt_lo_lo[ 1 ], this.user_low[ 2 ] ],
      [ vgt_lo_lo[ 0 ], this.user_low[ 1 ], this.user_low[ 2 ] ],
      [ vgt_lo_lo[ 0 ], this.user_low[ 1 ], vgt_lo_lo[ 2 ] ],
      1, 1, 1, 1 );
    this.colorQuadrilateral(
      [ this.user_low[ 0 ], vgt_lo_lo[ 1 ], vgt_lo_lo[ 2 ] ],
      [ this.user_low[ 0 ], this.user_low[ 1 ], vgt_lo_lo[ 2 ] ],
      [ this.user_low[ 0 ], this.user_low[ 1 ], this.user_low[ 2 ] ],
      [ this.user_low[ 0 ], vgt_lo_lo[ 1 ], this.user_low[ 2 ] ],
      1, 1, 1, 1 );
    this.colorQuadrilateral(
      [ vgt_lo_lo[ 0 ], vgt_lo_lo[ 1 ], vgt_lo_lo[ 2 ] ],
      [ this.user_low[ 0 ], vgt_lo_lo[ 1 ], vgt_lo_lo[ 2 ] ],
      [ this.user_low[ 0 ], vgt_lo_lo[ 1 ], this.user_low[ 2 ] ],
      [ vgt_lo_lo[ 0 ], vgt_lo_lo[ 1 ], this.user_low[ 2 ] ],
      .5, .5, .5, .5 );
    this.colorQuadrilateral(
      [ vgt_lo_lo[ 0 ], this.user_low[ 1 ], vgt_lo_lo[ 2 ] ],
      [ vgt_lo_lo[ 0 ], this.user_low[ 1 ], this.user_low[ 2 ] ],
      [ this.user_low[ 0 ], this.user_low[ 1 ], this.user_low[ 2 ] ],
      [ this.user_low[ 0 ], this.user_low[ 1 ], vgt_lo_lo[ 2 ] ],
      .5, .5, .5, .5 );
    this.colorQuadrilateral(
      [ vgt_lo_lo[ 0 ], vgt_lo_lo[ 1 ], vgt_lo_lo[ 2 ] ],
      [ vgt_lo_lo[ 0 ], this.user_low[ 1 ], vgt_lo_lo[ 2 ] ],
      [ this.user_low[ 0 ], this.user_low[ 1 ], vgt_lo_lo[ 2 ] ],
      [ this.user_low[ 0 ], vgt_lo_lo[ 1 ], vgt_lo_lo[ 2 ] ],
      0, 0, 0, 0 );
    this.colorQuadrilateral(
      [ vgt_lo_lo[ 0 ], vgt_lo_lo[ 1 ], this.user_low[ 2 ] ],
      [ this.user_low[ 0 ], vgt_lo_lo[ 1 ], this.user_low[ 2 ] ],
      [ this.user_low[ 0 ], this.user_low[ 1 ], this.user_low[ 2 ] ],
      [ vgt_lo_lo[ 0 ], this.user_low[ 1 ], this.user_low[ 2 ] ],
      0, 0, 0, 0 );
  this.endDrawing();
//document.getElementById("debug_textarea").value +=
//  "leaving VolGraphTool.makeOriginBox\n";
}; 
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.makeBoundingBox = function() {
//document.getElementById("debug_textarea").value +=
//  "entering VolGraphTool.makeBoundingBox\n";
  this.clip_bounding_box_array = new Array();
  this.bounding_box_array = new Array();

  var vgt_vertices = new Array();
  var vgt_normals = new Array();
  var vgt_colors = new Array();
  var vgt_indices = new Array();
  var vgt_n = 0;
  for ( var vgt_d = 0; vgt_d < 3; vgt_d ++ ) {
    var vgt_d1 = ( vgt_d + 1 ) % 3;
    var vgt_d2 = ( vgt_d + 2 ) % 3;
    var vgt_c = [ 0, 0, 0 ];
    vgt_c[ vgt_d ] = 1;
    var vgt_lo = new Array( 3 );
    var vgt_hi = new Array( 3 );
    for ( var vgt_i = 0; vgt_i < 3; vgt_i ++ ) {
      vgt_lo[ vgt_i ] = vgt_hi[ vgt_i ] = this.user_low[ vgt_i ];
    }
    vgt_lo[ vgt_d ] = this.clip_low[ vgt_d ];
    vgt_hi[ vgt_d ] = this.clip_high[ vgt_d ];
    var vgt_normal = [ 0, 0, 0 ];

    var vgt_oc = this.userCoordsToObjectCoords( vgt_lo );
    vgt_vertices.push( vgt_oc[ 0 ], vgt_oc[ 1 ], vgt_oc[ 2 ] );
    vgt_normal[ vgt_d1 ] = Math.SQRT1_2;
    vgt_normal[ vgt_d2 ] = Math.SQRT1_2;
    vgt_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ], vgt_normal[ 2 ] );
    vgt_colors.push( vgt_c[ 0 ], vgt_c[ 1 ], vgt_c[ 2 ], 1 );
    vgt_indices.push( vgt_n++ );
    vgt_oc = this.userCoordsToObjectCoords( vgt_hi );
    vgt_vertices.push( vgt_oc[ 0 ], vgt_oc[ 1 ], vgt_oc[ 2 ] );
    vgt_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ], vgt_normal[ 2 ] );
    vgt_colors.push( vgt_c[ 0 ], vgt_c[ 1 ], vgt_c[ 2 ], 1 );
    vgt_indices.push( vgt_n++ );

    vgt_lo[ vgt_d1 ] = vgt_hi[ vgt_d1 ] = this.user_high[ vgt_d1 ];
    var vgt_oc = this.userCoordsToObjectCoords( vgt_lo );
    vgt_vertices.push( vgt_oc[ 0 ], vgt_oc[ 1 ], vgt_oc[ 2 ] );
    vgt_normal[ vgt_d1 ] = - Math.SQRT1_2;
    vgt_normal[ vgt_d2 ] = Math.SQRT1_2;
    vgt_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ], vgt_normal[ 2 ] );
    vgt_colors.push( vgt_c[ 0 ], vgt_c[ 1 ], vgt_c[ 2 ], 1 );
    vgt_indices.push( vgt_n++ );
    vgt_oc = this.userCoordsToObjectCoords( vgt_hi );
    vgt_vertices.push( vgt_oc[ 0 ], vgt_oc[ 1 ], vgt_oc[ 2 ] );
    vgt_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ], vgt_normal[ 2 ] );
    vgt_colors.push( vgt_c[ 0 ], vgt_c[ 1 ], vgt_c[ 2 ], 1 );
    vgt_indices.push( vgt_n++ );

    vgt_lo[ vgt_d2 ] = vgt_hi[ vgt_d2 ] = this.user_high[ vgt_d2 ];
    var vgt_oc = this.userCoordsToObjectCoords( vgt_lo );
    vgt_vertices.push( vgt_oc[ 0 ], vgt_oc[ 1 ], vgt_oc[ 2 ] );
    vgt_normal[ vgt_d1 ] = - Math.SQRT1_2;
    vgt_normal[ vgt_d2 ] = - Math.SQRT1_2;
    vgt_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ], vgt_normal[ 2 ] );
    vgt_colors.push( vgt_c[ 0 ], vgt_c[ 1 ], vgt_c[ 2 ], 1 );
    vgt_indices.push( vgt_n++ );
    vgt_oc = this.userCoordsToObjectCoords( vgt_hi );
    vgt_vertices.push( vgt_oc[ 0 ], vgt_oc[ 1 ], vgt_oc[ 2 ] );
    vgt_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ], vgt_normal[ 2 ] );
    vgt_colors.push( vgt_c[ 0 ], vgt_c[ 1 ], vgt_c[ 2 ], 1 );
    vgt_indices.push( vgt_n++ );

    vgt_lo[ vgt_d1 ] = vgt_hi[ vgt_d1 ] = this.user_low[ vgt_d1 ];
    var vgt_oc = this.userCoordsToObjectCoords( vgt_lo );
    vgt_vertices.push( vgt_oc[ 0 ], vgt_oc[ 1 ], vgt_oc[ 2 ] );
    vgt_normal[ vgt_d1 ] = Math.SQRT1_2;
    vgt_normal[ vgt_d2 ] = - Math.SQRT1_2;
    vgt_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ], vgt_normal[ 2 ] );
    vgt_colors.push( vgt_c[ 0 ], vgt_c[ 1 ], vgt_c[ 2 ], 1 );
    vgt_indices.push( vgt_n++ );
    vgt_oc = this.userCoordsToObjectCoords( vgt_hi );
    vgt_vertices.push( vgt_oc[ 0 ], vgt_oc[ 1 ], vgt_oc[ 2 ] );
    vgt_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ], vgt_normal[ 2 ] );
    vgt_colors.push( vgt_c[ 0 ], vgt_c[ 1 ], vgt_c[ 2 ], 1 );
    vgt_indices.push( vgt_n++ );
  }

  var vgt_vertex_position_buffer = this.gl.createBuffer();
  this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_position_buffer );
  this.gl.bufferData( this.gl.ARRAY_BUFFER, new Float32Array(vgt_vertices),
    this.gl.STATIC_DRAW );

  var vgt_vertex_normal_buffer = this.gl.createBuffer();
  this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_normal_buffer );
  this.gl.bufferData( this.gl.ARRAY_BUFFER,
    new Float32Array( this.current_normals ), this.gl.STATIC_DRAW );
    vgt_vertex_normal_buffer.my_numItems = vgt_normals.length / 3;

  vgt_vertex_color_buffer = this.gl.createBuffer();
  this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_color_buffer );
  this.gl.bufferData( this.gl.ARRAY_BUFFER,
    new Float32Array( vgt_colors ), this.gl.STATIC_DRAW );

  var vgt_vertex_index_buffer = this.gl.createBuffer();
  this.gl.bindBuffer( this.gl.ELEMENT_ARRAY_BUFFER,
    vgt_vertex_index_buffer );
  this.gl.bufferData( this.gl.ELEMENT_ARRAY_BUFFER,
     new Uint16Array( vgt_indices ), this.gl.STATIC_DRAW );
  vgt_vertex_index_buffer.my_numItems = vgt_indices.length;

  this.clip_bounding_box_array.push( {
    vpb : vgt_vertex_position_buffer,
    vnb : vgt_vertex_normal_buffer,
    vcb : vgt_vertex_color_buffer,
    vib : vgt_vertex_index_buffer
  } );

  var vgt_vertices = new Array();
  var vgt_colors = new Array();
  var vgt_indices = new Array();
  var vgt_n = 0;
  for ( var vgt_d = 0; vgt_d < 3; vgt_d ++ ) {
    var vgt_d1 = ( vgt_d + 1 ) % 3;
    var vgt_d2 = ( vgt_d + 2 ) % 3;
    var vgt_c = [ 0, 0, 0 ];
    vgt_c[ vgt_d ] = 1;
    var vgt_lo = new Array( 3 );
    var vgt_hi = new Array( 3 );
    for ( var vgt_i = 0; vgt_i < 3; vgt_i ++ ) {
      vgt_lo[ vgt_i ] = vgt_hi[ vgt_i ] = this.user_low[ vgt_i ];
    }
    vgt_hi[ vgt_d ] = this.user_high[ vgt_d ];
    var vgt_oc = this.userCoordsToObjectCoords( vgt_lo );
    vgt_vertices.push( vgt_oc[ 0 ], vgt_oc[ 1 ], vgt_oc[ 2 ] );
    vgt_normal[ vgt_d1 ] = Math.SQRT1_2;
    vgt_normal[ vgt_d2 ] = Math.SQRT1_2;
    vgt_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ], vgt_normal[ 2 ] );
    vgt_colors.push( vgt_c[ 0 ], vgt_c[ 1 ], vgt_c[ 2 ], 1 );
    vgt_indices.push( vgt_n++ );
    var vgt_oc = this.userCoordsToObjectCoords( vgt_hi );
    vgt_vertices.push( vgt_oc[ 0 ], vgt_oc[ 1 ], vgt_oc[ 2 ] );
    vgt_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ], vgt_normal[ 2 ] );
    vgt_colors.push( vgt_c[ 0 ], vgt_c[ 1 ], vgt_c[ 2 ], 1 );
    vgt_indices.push( vgt_n++ );

    vgt_lo[ vgt_d1 ] = vgt_hi[ vgt_d1 ] = this.user_high[ vgt_d1 ];
    var vgt_oc = this.userCoordsToObjectCoords( vgt_lo );
    vgt_vertices.push( vgt_oc[ 0 ], vgt_oc[ 1 ], vgt_oc[ 2 ] );
    vgt_normal[ vgt_d1 ] = - Math.SQRT1_2;
    vgt_normal[ vgt_d2 ] = Math.SQRT1_2;
    vgt_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ], vgt_normal[ 2 ] );
    vgt_colors.push( vgt_c[ 0 ], vgt_c[ 1 ], vgt_c[ 2 ], 1 );
    vgt_indices.push( vgt_n++ );
    var vgt_oc = this.userCoordsToObjectCoords( vgt_hi );
    vgt_vertices.push( vgt_oc[ 0 ], vgt_oc[ 1 ], vgt_oc[ 2 ] );
    vgt_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ], vgt_normal[ 2 ] );
    vgt_colors.push( vgt_c[ 0 ], vgt_c[ 1 ], vgt_c[ 2 ], 1 );
    vgt_indices.push( vgt_n++ );

    vgt_lo[ vgt_d2 ] = vgt_hi[ vgt_d2 ] = this.user_high[ vgt_d2 ];
    var vgt_oc = this.userCoordsToObjectCoords( vgt_lo );
    vgt_vertices.push( vgt_oc[ 0 ], vgt_oc[ 1 ], vgt_oc[ 2 ] );
    vgt_normal[ vgt_d1 ] = - Math.SQRT1_2;
    vgt_normal[ vgt_d2 ] = - Math.SQRT1_2;
    vgt_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ], vgt_normal[ 2 ] );
    vgt_colors.push( vgt_c[ 0 ], vgt_c[ 1 ], vgt_c[ 2 ], 1 );
    vgt_indices.push( vgt_n++ );
    var vgt_oc = this.userCoordsToObjectCoords( vgt_hi );
    vgt_vertices.push( vgt_oc[ 0 ], vgt_oc[ 1 ], vgt_oc[ 2 ] );
    vgt_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ], vgt_normal[ 2 ] );
    vgt_colors.push( vgt_c[ 0 ], vgt_c[ 1 ], vgt_c[ 2 ], 1 );
    vgt_indices.push( vgt_n++ );

    vgt_lo[ vgt_d1 ] = vgt_hi[ vgt_d1 ] = this.user_low[ vgt_d1 ];
    var vgt_oc = this.userCoordsToObjectCoords( vgt_lo );
    vgt_vertices.push( vgt_oc[ 0 ], vgt_oc[ 1 ], vgt_oc[ 2 ] );
    vgt_normal[ vgt_d1 ] = Math.SQRT1_2;
    vgt_normal[ vgt_d2 ] = - Math.SQRT1_2;
    vgt_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ], vgt_normal[ 2 ] );
    vgt_colors.push( vgt_c[ 0 ], vgt_c[ 1 ], vgt_c[ 2 ], 1 );
    vgt_indices.push( vgt_n++ );
    var vgt_oc = this.userCoordsToObjectCoords( vgt_hi );
    vgt_vertices.push( vgt_oc[ 0 ], vgt_oc[ 1 ], vgt_oc[ 2 ] );
    vgt_normals.push( vgt_normal[ 0 ], vgt_normal[ 1 ], vgt_normal[ 2 ] );
    vgt_colors.push( vgt_c[ 0 ], vgt_c[ 1 ], vgt_c[ 2 ], 1 );
    vgt_indices.push( vgt_n++ );
  }

  var vgt_vertex_position_buffer = this.gl.createBuffer();
  this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_position_buffer );
  this.gl.bufferData( this.gl.ARRAY_BUFFER, new Float32Array(vgt_vertices),
    this.gl.STATIC_DRAW );

  vgt_vertex_normal_buffer = this.gl.createBuffer();
  this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_normal_buffer );
  this.gl.bufferData( this.gl.ARRAY_BUFFER,
    new Float32Array( this.current_normals ), this.gl.STATIC_DRAW );
    vgt_vertex_normal_buffer.my_numItems = vgt_normals.length / 3;

  vgt_vertex_color_buffer = this.gl.createBuffer();
  this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_color_buffer );
  this.gl.bufferData( this.gl.ARRAY_BUFFER,
    new Float32Array( vgt_colors ), this.gl.STATIC_DRAW );

  var vgt_vertex_index_buffer = this.gl.createBuffer();
  this.gl.bindBuffer( this.gl.ELEMENT_ARRAY_BUFFER,
    vgt_vertex_index_buffer );
  this.gl.bufferData( this.gl.ELEMENT_ARRAY_BUFFER,
     new Uint16Array( vgt_indices ), this.gl.STATIC_DRAW );
  vgt_vertex_index_buffer.my_numItems = vgt_indices.length;

  this.bounding_box_array.push( {
    vpb : vgt_vertex_position_buffer,
    vnb : vgt_vertex_normal_buffer,
    vcb : vgt_vertex_color_buffer,
    vib : vgt_vertex_index_buffer
  } );
//document.getElementById("debug_textarea").value +=
//  "leaving VolGraphTool.makeBoundingBox\n";
}; 
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.beginCurve = function( uc , r ) {
  this.current_vertices = new Array();
  this.current_colors = new Array();
  this.current_vertex_indices = new Array();
  this.current_vertex_number = 0;
  this.continueCurve( uc, r );
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.continueCurve = function( uc , r ) {
  var vgt_oc = this.userCoordsToObjectCoords( uc );
  this.current_vertices.push( vgt_oc[ 0 ], vgt_oc[ 1 ], vgt_oc[ 2 ] );
  if ( !r ) r = 1;
  var vgt_c = ( !r ? this.colormap.getGLColor( 1 ) : // foreground color
    this.colormap.getGLColorFromRatio( r ) );
  this.current_colors.push( vgt_c[ 0 ], vgt_c[ 1 ], vgt_c[ 2 ], 1 );
  this.current_vertex_indices.push( this.current_vertex_number );
  this.current_vertex_number ++;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.endCurve = function( uc , r ) {
// http://blogoben.wordpress.com/2011/04/16/webgl-basics-4-wireframe-3d-object/
//LINE_STRIP: 2 first vertices delimit the ends of the first segment,
//  each new vertex defines the end of a new segment starting at the end
//  of the previous one
//LINE_LOOP: same as LINE_STRIP with an additional segment between
//  the first and last vertices (closing the contiguous segments set)
//LINE: pairs of vertices delimit the ends of each individual segment,
//  giving a non-contiguous set of segments

  this.continueCurve( uc, r );
  var vgt_vertex_position_buffer = this.gl.createBuffer();
  this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_position_buffer );
  this.gl.bufferData( this.gl.ARRAY_BUFFER,
    new Float32Array(this.current_vertices),this.gl.STATIC_DRAW );
  this.gl.vertexAttribPointer(
    this.programInfo.attribLocations.vertexPosition, 3, this.gl.FLOAT,
    false, 0, 0);
  this.gl.enableVertexAttribArray(
    this.programInfo.attribLocations.vertexPosition);
  this.gl.uniformMatrix4fv( this.programInfo.projectionMatrix, 0,
    this.p_matrix );
  this.gl.drawArrays( this.gl.LINE_STRIP, 0,
    this.current_vertices.length / 3 );

  vgt_vertex_color_buffer = this.gl.createBuffer();
  this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_color_buffer );
  this.gl.bufferData( this.gl.ARRAY_BUFFER,
    new Float32Array( this.current_colors ), this.gl.STATIC_DRAW );
  vgt_vertex_color_buffer.my_numItems = this.current_colors.length / 4;

  var vgt_vertex_index_buffer = this.gl.createBuffer();
  this.gl.bindBuffer( this.gl.ELEMENT_ARRAY_BUFFER,
    vgt_vertex_index_buffer );
  this.gl.bufferData( this.gl.ELEMENT_ARRAY_BUFFER,
    new Uint16Array( this.current_vertex_indices ), this.gl.STATIC_DRAW );
  vgt_vertex_index_buffer.my_numItems = this.current_vertex_indices.length;

  this.curve_drawing_array.push( {
    vpb : vgt_vertex_position_buffer,
    vcb : vgt_vertex_color_buffer,
    vib : vgt_vertex_index_buffer
  } );
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.screenCoordsToUnitSphereCoords =
function( screen_coords ) {
  var vgt_unit_square_coords = [
    2 * screen_coords[ 0 ] / this.canvas.width - 1, 
    1 - 2 * screen_coords[ 1 ] / this.canvas.height
  ];
  var vgt_norm_squared = Math.pow( vgt_unit_square_coords[ 0 ], 2 )
                   + Math.pow( vgt_unit_square_coords[ 1 ], 2 );
  var vgt_sigma_squared =
    Math.min( TRACKBALLSIZE * TRACKBALLSIZE, vgt_norm_squared );
  unit_sphere_coords = vec3.create();
  unit_sphere_coords[ 2 ] =
    Math.sqrt( TRACKBALLSIZE * TRACKBALLSIZE - vgt_sigma_squared );
  var vgt_norm = Math.sqrt( vgt_norm_squared
                      + Math.pow( unit_sphere_coords[ 2 ], 2 ) );
  unit_sphere_coords[ 0 ] = vgt_unit_square_coords[ 0 ] / vgt_norm;
  unit_sphere_coords[ 1 ] = vgt_unit_square_coords[ 1 ] / vgt_norm;
  unit_sphere_coords[ 2 ] /= vgt_norm;
  return unit_sphere_coords;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.rotateImage = function( ) {
  this.new_unit_sphere_coords =
    this.screenCoordsToUnitSphereCoords( this.new_screen_coords );

  var vgt_axis = vec3.create();
  vec3.cross( this.old_unit_sphere_coords, this.new_unit_sphere_coords,
    vgt_axis );
  var vgt_axis_norm = vec3.length( vgt_axis );
  if ( vgt_axis_norm <= 0 ) return;
  vec3.scale( vgt_axis, 1 / vgt_axis_norm );
  var vgt_dot_product = vec3.dot( vgt_axis, this.old_unit_sphere_coords );
  var vgt_axis_dot = vec3.create();
  vec3.scale( vgt_axis, vgt_dot_product, vgt_axis_dot );
  var vgt_perp = vec3.create();
  vec3.subtract( this.old_unit_sphere_coords, vgt_axis_dot, vgt_perp ); 
  var vgt_axis_cross_perp = vec3.create();
  vec3.cross( vgt_axis, vgt_perp, vgt_axis_cross_perp );
  var vgt_axis_cross_perp_norm_squared =
    vec3.dot( vgt_axis_cross_perp, vgt_axis_cross_perp );
  if ( vgt_axis_cross_perp_norm_squared <= 0 ) return;
  var vgt_sine_angle = vec3.dot( this.new_unit_sphere_coords,
    vgt_axis_cross_perp ) / vgt_axis_cross_perp_norm_squared;
  vgt_sine_angle = Math.max( -1, Math.min( 1, - vgt_sine_angle ) );
  var vgt_half_angle = 0.5 * Math.asin( vgt_sine_angle );
  var vgt_sine_half_angle = Math.sin( vgt_half_angle );
  var vgt_rotate_quat = quat4.create();
  vgt_rotate_quat[ 0 ] = vgt_axis[ 0 ] * vgt_sine_half_angle; 
  vgt_rotate_quat[ 1 ] = vgt_axis[ 1 ] * vgt_sine_half_angle; 
  vgt_rotate_quat[ 2 ] = vgt_axis[ 2 ] * vgt_sine_half_angle; 
  vgt_rotate_quat[ 3 ] = Math.cos( vgt_half_angle );
//quat4.multiply( vgt_rotate_quat, this.old_quat, this.current_quat ); //OK
  quat4.multiply( this.old_quat, vgt_rotate_quat, this.current_quat );
  quat4.toMat4( this.current_quat, this.mv_matrix ); // good
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.nearestPointOnLine = function( v0, t0, v1, t1 ) {
//document.getElementById("debug_textarea").value +=
//  "entering nearestPointOnLine\n";
  var vgt_vdif = vec3.create( v1 );
  vec3.subtract( vgt_vdif, v0 );
  var vgt_r00 = vec3.length( t0 );
  var vgt_q0 = vec3.create( t0 );
  vec3.scale( vgt_q0, 1 / vgt_r00 );
  var vgt_r01 = vec3.dot( vgt_q0, t1 );
  var vgt_b0 = vec3.dot( vgt_q0, vgt_vdif );
  var vgt_q1 = vec3.create( t1 );
  var vgt_q0scale = vec3.create( vgt_q0 );
  vec3.subtract( vgt_q1, vec3.scale( vgt_q0scale, vgt_r01 ) );
  vgt_q0scale = vec3.create( vgt_q0 );
  vec3.subtract( vgt_vdif, vec3.scale( vgt_q0scale, vgt_b0 ) ); 
  var vgt_r11 = vec3.length( vgt_q1 );
  if ( vgt_r11 <= 0 ) { // parallel lines
    return {
      a0 : 0, // position is arbitrary
      a1 : 0, // position is arbitrary
      min_dist : vec3.length( vgt_vdif )
    };
  }
  vec3.scale( vgt_q1, 1 / vgt_r11 );
  var vgt_b1 = vec3.dot( vgt_q1, vgt_vdif );
  var vgt_q1scale = vec3.create( vgt_q1 );
  vec3.subtract( vgt_vdif, vec3.scale( vgt_q1scale, vgt_b1 ) );
  var vgt_alpha1 = - vgt_b1 / vgt_r11;
  var vgt_alpha0 = ( vgt_b0 + vgt_r01 * vgt_alpha1 ) / vgt_r00;
//document.getElementById("debug_textarea").value +=
//  "leaving nearestPointOnLine\n";
  return {
    a0 : vgt_alpha0, // v0 + t0 * a0 is nearest point on first line
    a1 : vgt_alpha1, // v1 + t1 * a1 is nearest point on second line
    min_dist : vec3.length( vgt_vdif ) // min distance between lines
  }
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.screenCoordsToClipCoords = function( sc ) {
// screen coords lie in (0,w)x(0,h)
// clip coords lie in (-1,1)x(-1,1)
  return [
    ( 2 * sc[ 0 ] / this.canvas.width - 1 ),
    ( 1 - 2 * sc[ 1 ] / this.canvas.height )
  ];
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.clipCoordsToScreenCoords = function( cc ) {
// screen coords lie in (0,w)x(0,h)
// clip coords lie in (-1,1)x(-1,1)
  return [ ( cc[ 0 ] + 1 ) * this.canvas.width * .5,
           ( 1 - cc[ 1 ] ) * this.canvas.height * .5 ];
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.selectClipPlane = function( ) {
//document.getElementById("debug_textarea").value +=
//  "entering selectClipPlane\n";
  mat4.inverse( this.p_matrix , this.p_matrix_inverse );
  mat4.inverse( this.mv_matrix, this.mv_matrix_inverse );
  var vgt_t1 = [ 0, 0, 1, 0 ]; // direction of line through point on screen
                       // in homogeneous clip coordinates
  mat4.multiplyVec4( this.p_matrix_inverse, vgt_t1 ); //direction in eye coords
  mat4.multiplyVec4( this.mv_matrix_inverse, vgt_t1 );//direction in obj coords
//document.getElementById("debug_textarea").value +=
//  "vgt_t1 obj coords = " + vgt_t1[ 0 ] + " " + vgt_t1[ 1 ] + " " + vgt_t1[ 2 ] + "\n";
//point on user-selected line in homogeneous clip coordinates
  var vgt_clip_coords =
    this.screenCoordsToClipCoords( this.new_screen_coords );
  var vgt_v1 = [ vgt_clip_coords[ 0 ], vgt_clip_coords[ 1 ], 0, 1 ];
  mat4.multiplyVec4( this.p_matrix_inverse, vgt_v1 ); // now in eye coords
  mat4.multiplyVec4( this.mv_matrix_inverse, vgt_v1 ); // now in object coords
//document.getElementById("debug_textarea").value +=
//  "vgt_v1 obj coords = " + vgt_v1[ 0 ] + " " + vgt_v1[ 1 ] + " " + vgt_v1[ 2 ] + "\n";
  var vgt_min_distance = Number.MAX_VALUE;
  this.clip_direction = 0;
  var vgt_clip_position_object_coords = undefined;
  var vgt_alpha0 = undefined;
  var vgt_alpha1 = undefined;
  for ( var vgt_d = 0; vgt_d < 3; vgt_d ++ ) {
    var vgt_t0 = [ 0, 0, 0 ];
    vgt_t0[ vgt_d ] = 1;
    var vgt_v0 = [ -1, -1, -1 ];
    vgt_v0[ vgt_d ] = 0;
    var vgt_d1 = ( vgt_d + 1 ) %3;
    var vgt_d2 = ( vgt_d + 2 ) %3;
    for ( var vgt_c1 = 0; vgt_c1 < 2; vgt_c1++, vgt_v0[ vgt_d1 ] = 1 ) {
      vgt_v0[ vgt_d2 ] = -1;
      for ( var vgt_c2 = 0; vgt_c2 < 2; vgt_c2++, vgt_v0[ vgt_d2 ] = 1 ) {
        var vgt_stuff = VolGraphTool.nearestPointOnLine( vgt_v0, vgt_t0,
          vgt_v1, vgt_t1 );
//      document.getElementById("debug_textarea").value +=
//        "vgt_d,vgt_c1,vgt_c2,distance = " + vgt_d + " " + vgt_c1 + " "
//        + vgt_c2 + " " + vgt_stuff.min_dist + "\n";
        if ( vgt_stuff.min_dist < vgt_min_distance && vgt_stuff.a0 >= -1 &&
        vgt_stuff.a0 <= 1 ) {
          vgt_min_distance = vgt_stuff.min_dist;
          this.clip_direction = vgt_d;
          vgt_clip_position_object_coords = vgt_stuff.a0;
          this.clip_edge_vertex = [ vgt_v0[ 0 ], vgt_v0[ 1 ], vgt_v0[ 2 ] ];
          vgt_alpha0 = vgt_stuff.a0;
          vgt_alpha1 = vgt_stuff.a1;
        }
      }
    }
  }
//document.getElementById("debug_textarea").value +=
//  "clip_direction = " + this.clip_direction + "\n";
//document.getElementById("debug_textarea").value +=
//  "vgt_clip_position_object_coords " + vgt_clip_position_object_coords + "\n";
  var vgt_clip_position_object_coords =
    Math.max( -1, Math.min( 1, vgt_clip_position_object_coords ) );
  var vgt_clip_position_user_coords = this.user_center[ this.clip_direction ]
    + vgt_clip_position_object_coords
    * this.user_len[ this.clip_direction ] * 0.5;
//document.getElementById("debug_textarea").value +=
//  "vgt_clip_position_user_coords" + vgt_clip_position_user_coords + "\n";
  if ( vgt_clip_position_user_coords * 2 <=
  this.clip_low[ this.clip_direction ]
  + this.clip_high[ this.clip_direction ] ) {
    this.clip_low[ this.clip_direction ] = vgt_clip_position_user_coords;
  } else {
    this.clip_high[ this.clip_direction ] = vgt_clip_position_user_coords;
  }
//document.getElementById("debug_textarea").value +=
//  "clip low, high = " + this.clip_low[ this.clip_direction ] + " "
//  + this.clip_high[ this.clip_direction ] + "\n";
//document.getElementById("debug_textarea").value +=
//  "leaving selectClipPlane\n";
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.moveClipPlane = function( ) {
//document.getElementById("debug_textarea").value +=
//  "entering moveClipPlane\n";
  var vgt_t1 = [ 0, 0, 1, 0 ];
    // direction of line normal to screen in clip coords
  mat4.multiplyVec4( this.p_matrix_inverse, vgt_t1 );//direction in eye coord
  mat4.multiplyVec4( this.mv_matrix_inverse, vgt_t1 );//dir in obj coords
//document.getElementById("debug_textarea").value +=
//  "vgt_t1 obj coords = " + vgt_t1[ 0 ] + " " + vgt_t1[ 1 ] + " " + vgt_t1[ 2 ] + "\n";
  var vgt_clip_coords =
    this.screenCoordsToClipCoords( this.new_screen_coords );
//point on line through current cursor in homogeneous clip coordinates
  var vgt_v1 = [ vgt_clip_coords[ 0 ], vgt_clip_coords[ 1 ], 0, 1 ];
  mat4.multiplyVec4( this.p_matrix_inverse, vgt_v1 ); // now in eye coords
  mat4.multiplyVec4( this.mv_matrix_inverse, vgt_v1 ); //now in object coords
//document.getElementById("debug_textarea").value +=
//  "vgt_v1 obj coords = " + vgt_v1[ 0 ] + " " + vgt_v1[ 1 ] + " "
//  + vgt_v1[ 2 ] + "\n";
  var vgt_t0 = [ 0, 0, 0 ];
  vgt_t0[ this.clip_direction ] = 1;
  var vgt_stuff =
    VolGraphTool.nearestPointOnLine( this.clip_edge_vertex, vgt_t0, vgt_v1,
      vgt_t1 );
  var vgt_clip_position_object_coords =
    Math.max( -1, Math.min( 1, vgt_stuff.a0 ) );
  vgt_clip_position_user_coords = this.user_center[ this.clip_direction ]
    + vgt_clip_position_object_coords
    * this.user_len[ this.clip_direction ] * 0.5;
//document.getElementById("debug_textarea").value +=
//  "vgt_clip_position_object_coords = " + vgt_clip_position_object_coords
//  + "\n";
  var vgt_ncp = 2 * this.clip_direction;
  if ( vgt_clip_position_user_coords * 2 <=
  this.clip_low[ this.clip_direction ]
  + this.clip_high[ this.clip_direction ] ) {
    this.clip_low[ this.clip_direction ] = vgt_clip_position_user_coords;
    this.drawClipPlane( this.clip_direction, 0,
      vgt_clip_position_user_coords );
  } else {
    this.clip_high[ this.clip_direction ] = vgt_clip_position_user_coords;
    this.drawClipPlane( this.clip_direction, 1,
      vgt_clip_position_user_coords );
    vgt_ncp ++;
  }
//document.getElementById("debug_textarea").value +=
//  "vgt_ncp = " + vgt_ncp + "\n";
//document.getElementById("debug_textarea").value +=
//  "clip low, high = " + this.clip_low[ this.clip_direction ] + " "
//  + this.clip_high[ this.clip_direction ] + "\n";
  this.makeBoundingBox();
//remember that makeOriginBox stored in the first surface
//drawClipPlane pushed recent drawing into last surface_drawing_array entry
  if (this.drawing_clip_planes ) {
    this.surface_drawing_array[ vgt_ncp + 1 ] =
      this.surface_drawing_array.pop();
  }
//document.getElementById("debug_textarea").value +=
//  "leaving moveClipPlane\n";
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.markFaceLocation = function( ) {
  mat4.inverse( this.p_matrix , this.p_matrix_inverse );
  mat4.inverse( this.mv_matrix, this.mv_matrix_inverse );

  var vgt_t = [ 0, 0, 1, 0 ]; // tangent for line through point on screen
                      // in homogeneous clip coordinates
  mat4.multiplyVec4( this.p_matrix_inverse, vgt_t );
  mat4.multiplyVec4( this.mv_matrix_inverse, vgt_t ); // now in object coords
//point on user-selected line in homogeneous clip coordinates
  var vgt_clip_coords =
    this.screenCoordsToClipCoords( this.new_screen_coords );
  var vgt_v = [ vgt_clip_coords[ 0 ], vgt_clip_coords[ 1 ], 0, 1 ];
  mat4.multiplyVec4( this.p_matrix_inverse, vgt_v );
  mat4.multiplyVec4( this.mv_matrix_inverse, vgt_v ); // now in object coords
  var vgt_max_z = - Number.MAX_VALUE;
  this.cursor_position = new Array( 4 );
  for ( var vgt_d = 0; vgt_d < 3; vgt_d ++ ) {
    for ( var vgt_vd = -1; vgt_vd < 2; vgt_vd += 2 ) {
      if ( Math.abs( vgt_t[ vgt_d ] ) > 0 ) {
        var vgt_a = ( vgt_vd - vgt_v[ vgt_d ] ) / vgt_t[ vgt_d ];
        var vgt_q = new Array( 4 ); // homogeneous object coords
        for ( var vgt_i = 0; vgt_i < 4; vgt_i ++ ) {
          vgt_q[ vgt_i ] = vgt_v[ vgt_i ] + vgt_t[ vgt_i ] * vgt_a;
        }
          //obj coord
        var vgt_c = new Array( 4 );
        for ( var vgt_i = 0; vgt_i < 4; vgt_i ++ ) {
          vgt_c[ vgt_i ] = vgt_q[ vgt_i ];
        }
        mat4.multiplyVec4( this.mv_matrix, vgt_c );
        mat4.multiplyVec4( this.mv_matrix, vgt_c ); // now in clip coords
        if ( Math.abs( vgt_q[ ( vgt_d + 1 ) % 3 ] ) <= 1 &&
        Math.abs( vgt_q[ ( vgt_d + 2 ) % 3 ] ) <= 1 &&
        vgt_c[ 2 ] > vgt_max_z ) {
          vgt_max_z = vgt_c[ 2 ];
          for ( vgt_i = 0; vgt_i < 4; vgt_i ++ ) {
            this.cursor_position[ vgt_i ] = vgt_q[ vgt_i ];
          }
        }
      }
    }
  }
  this.drawCrossHairs();
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.drawCrossHairs = function( ) {
  var vgt_vertices = new Array();
  var vgt_colors = new Array();
  var vgt_indices = new Array();
  var vgt_n = 0;
  var vgt_start = new Array( 3 );
  var vgt_finish = new Array( 3 );
  for ( var vgt_d = 0; vgt_d < 3; vgt_d ++ ) {
    vec3.set( this.cursor_position, vgt_start );
    vec3.set( this.cursor_position, vgt_finish );
    vgt_start[ vgt_d ] = -1;
    vgt_finish[ vgt_d ] = 1;
      + "\n";
    var vgt_c = [ 0, 0, 0 ];
    vgt_c[ vgt_d ] = 1;
    vgt_vertices.push( vgt_start[ 0 ], vgt_start[ 1 ], vgt_start[ 2 ] );
    vgt_colors.push( vgt_c[ 0 ], vgt_c[ 1 ], vgt_c[ 2 ], 1 );
    vgt_indices.push( vgt_n++ );
    vgt_vertices.push( vgt_finish[ 0 ], vgt_finish[ 1 ], vgt_finish[ 2 ] );
    vgt_colors.push( vgt_c[ 0 ], vgt_c[ 1 ], vgt_c[ 2 ], 1 );
    vgt_indices.push( vgt_n++ );
  }
  var vgt_vertex_position_buffer = this.gl.createBuffer();
  this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_position_buffer );
  this.gl.bufferData( this.gl.ARRAY_BUFFER, new Float32Array(vgt_vertices),
    this.gl.STATIC_DRAW );

  vgt_vertex_color_buffer = this.gl.createBuffer();
  this.gl.bindBuffer( this.gl.ARRAY_BUFFER, vgt_vertex_color_buffer );
  this.gl.bufferData( this.gl.ARRAY_BUFFER,
    new Float32Array( vgt_colors ), this.gl.STATIC_DRAW );

  var vgt_vertex_index_buffer = this.gl.createBuffer();
  this.gl.bindBuffer( this.gl.ELEMENT_ARRAY_BUFFER,
    vgt_vertex_index_buffer );
  this.gl.bufferData( this.gl.ELEMENT_ARRAY_BUFFER,
     new Uint16Array( vgt_indices ), this.gl.STATIC_DRAW );
  vgt_vertex_index_buffer.my_numItems = vgt_indices.length;

  this.crosshair_array.push( {
    vpb : vgt_vertex_position_buffer,
    vcb : vgt_vertex_color_buffer,
    vib : vgt_vertex_index_buffer
  } );

  var user_coords = this.objectCoordsToUserCoords( this.cursor_position );
  document.getElementById("crosshair_textarea").value =
    "crosshairs centered at [ " + user_coords[ 0 ] + " "
    + user_coords[ 1 ] + " " + user_coords[ 2 ] + " ]";
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool.prototype.moveCrossHairs = function( ) {
//tangent to mouse motion in clip coordinates
  var vgt_cco = this.screenCoordsToClipCoords( this.old_screen_coords );
  var vgt_ccn = this.screenCoordsToClipCoords( this.new_screen_coords );
  var vgt_t = [ vgt_ccn[ 0 ] - vgt_cco[ 0 ], vgt_ccn[ 1 ] - vgt_cco[ 1 ],
    0, 0 ];//clip coords
  mat4.multiplyVec4( this.p_matrix_inverse, vgt_t ); // tangent in eye coords
  mat4.multiplyVec4( this.mv_matrix_inverse, vgt_t );// tangent in obj coords
  var vgt_max_t = 0;
  var vgt_direction = 0;
  for ( var vgt_d = 0; vgt_d < 3; vgt_d ++ ) {
    if ( Math.abs( vgt_t[ vgt_d ] ) > vgt_max_t ) {
      vgt_max_t = Math.abs( vgt_t[ vgt_d ] ); 
      vgt_direction = vgt_d;
    }
  }
  var new_pos = this.cursor_position[ vgt_direction ]
    + vgt_t[ vgt_direction ];
  this.cursor_position[ vgt_direction ] =
    Math.max( -1, Math.min( 1, new_pos ) );
/*
//mouse position is read only, cannot be set by JavaScript:
//  http://stackoverflow.com/questions/1208729/jquery-set-mouse-position-not-cursor-position
  quat4.set( this.cursor_position, vgt_t ); // cursor position in object coords
  mat4.multiplyVec4( this.mv_matrix, vgt_t ); // cursor in eye coords
  mat4.multiplyVec4( this.p_matrix, vgt_t );  // cursor in clip coords
  sc = this.clipCoordsToScreenCoords( vgt_t ); // cursor in screen coords
*/
  this.crosshair_array = new Array();
  this.drawCrossHairs();
//return sc;
}
