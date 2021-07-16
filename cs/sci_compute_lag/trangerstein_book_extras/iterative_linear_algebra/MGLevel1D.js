var SMOOTHER =
  { GAUSS_SEIDEL : 0, GAUSS_SEIDEL_RED_BLACK : 1, RICHARDSON : 2 };
var RESTRICTION_PROLONGATION = { ALGEBRAIC_MULTIGRID: 0, FINITE_ELEMENT: 1 };
var smoother_name =
  [ "Gauss-Seidel", "Gauss-Seidel Red-Black", "Richardson" ];
var restriction_prolongation_nam = 
  [ "algebraic multigrid", "finite element" ];

// number_cells,matrix,smoother_iterations,smoother,restriction_prolongation
function Level( nc, si, s, rp ) {
//document.getElementById("debug_textarea").value +=
//  "entering Level constructor, nc,si,s,rp = " + nc + " " + si + " " + s
//  + " " + rp + "\n";
  this.coarser = null;
  this.finer = null;
  this.ncells = nc;
  this.level_number = 0;
  this.matrix_copy = null;
  this.smoother_iterations = si;
  this.smoother = s;
  this.restriction_prolongation = rp;
  this.stride = nc + 1;
  this.ilast = nc - 1;
  if ( nc % 2 == 0 && nc > 2 ) {
    this.coarser = new Level( nc / 2, si, s, rp );
    this.coarser.finer = this;
  }
//document.getElementById("debug_textarea").value +=
//  "leaving Level constructor\n"
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Level.prototype.setup = function( m ) {
  this.matrix = new Array( 3 * this.stride );
  this.Ax_b = new Array( this.stride );
  this.first_residual = new Array( this.stride );
  this.solution_increment = new Array( this.stride );
  this.residual = new Array( this.stride );
  this.richardson_residual = new Array( this.stride );
  this.prolongation_vector = new Array( this.stride );
//for plotting:
  this.solution_increment0 = new Array( this.stride ); // from preSmooth
  this.residual1 = new Array( this.stride ); // before restrict
  this.solution_increment1 = new Array( this.stride ); // after prolong

  var stride2 = 2 * this.stride;
  if ( this.finer == null ) {
    for ( var i = 0; i <= 3 * this.stride; i ++ ) this.matrix[ i ] = m[ i ];
  } else {
    this.level_number = this.finer.level_number + 1;
//  compute A_c = Prolongation^T * A_f * Prolongation
    for ( var i = 1; i <= this.finer.ilast; i ++ ) {
      this.finer.first_residual[ i ] = 0.;
    }
    for ( var START = 1; START <= Math.min( 3, this.ilast ); START ++ ) {
      for ( var i = 1; i <= this.finer.ilast; i ++ ) {
        this.finer.solution_increment[ i ] = 0;
      }
      for ( var I = 1; I <= this.ilast; I ++ ) {
        this.solution_increment[ I ] = 0.;
      }
      for ( var I = START; I <= this.ilast; I += 3 ) {
        this.solution_increment[ I ] = 1.;
      }
      this.finer.prolongation();
      this.finer.computeResidual( this.finer.solution_increment,
        this.finer.first_residual, this.finer.Ax_b );
      this.finer.restriction();
      for ( var I = START; I <= this.ilast; I += 3 ) {
        if ( I > 1 ) {
          this.matrix[ I - 1 + stride2 ] = this.first_residual[ I - 1 ];
        }
        this.matrix[ I + this.stride ] = this.first_residual[ I ];
        if ( I < this.ilast ) {
          this.matrix[ I + 1 ] = this.first_residual[ I + 1 ];
        }
      }
    }
  }
  if ( smoother == SMOOTHER.RICHARDSON ) {
//  Use Gerschgorin to estimate Richardson smoother mu
    this.mu = Math.max(
        Math.abs( this.matrix[ 1 + this.stride ] )
      + Math.abs( this.matrix[ 1 + stride2 ] ),
        Math.abs( this.matrix[ this.ilast ] )
      + Math.abs( this.matrix[ this.ilast + this.stride ] ) );
    for ( var i = 2; i < this.ilast; i ++ ) {
      var sum = Math.abs( this.matrix[ i ] )
              + Math.abs( this.matrix[ i + this.stride ] )
              + Math.abs( this.matrix[ i + stride2 ] );
      this.mu = Math.max( this.mu, sum );
    }
  }
  for ( var i = 0; i <= this.ncells; i ++ ) {
    this.solution_increment[ i ] = 0.;
  }
  if ( this.coarser == null ) {
    this.matrix_copy = new Array( 3 * this.stride ); 
    for ( var i = 0; i < 3 * this.stride; i ++ ) {
      this.matrix_copy[ i ] = this.matrix[ i ];
    }
//  Gaussian factorization on coarsest level:
    for ( var i = 2; i <= this.ilast; i ++ ) {
      this.matrix_copy[ i ] /= this.matrix_copy[ i - 1 + this.stride ];
      this.matrix_copy[ i + this.stride ] -=
        this.matrix_copy[ i ] * this.matrix_copy[ i - 1 + stride2 ];
    }
  } else this.coarser.setup( undefined );
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Level.prototype.numberLevels = function() {
  if ( coarser == null ) return this.level_number + 1;
  else return coarser.numberLevels();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Level.prototype.preSmooth = function() {
  for ( var it = 0; it < this.smoother_iterations; it ++ ) {
    var stride2 = this.stride * 2;
    switch ( this.smoother ) {
      case SMOOTHER.GAUSS_SEIDEL : {
        var i = 1;
        var residual = 
            this.matrix[ i + this.stride ] * this.solution_increment[ i ]
          + this.matrix[ i + stride2 ] * this.solution_increment[ i + 1 ]
          - this.first_residual[ i ];
        this.solution_increment[ i ] -=
          residual / this.matrix[ i + this.stride ];
        for ( var i = 2; i <= this.ilast - 1; i ++ ) {
          residual = this.matrix[ i ] * this.solution_increment[ i - 1 ]
            + this.matrix[ i + this.stride ] * this.solution_increment[ i ]
            + this.matrix[ i + stride2 ] * this.solution_increment[ i + 1 ]
            - this.first_residual[ i ];
          this.solution_increment[ i ] -=
            residual / this.matrix[ i + this.stride ];
        }
        i = this.ilast;
        residual = this.matrix[ i ] * this.solution_increment[ i - 1 ]
          + this.matrix[ i + this.stride ] * this.solution_increment[ i ]
          - this.first_residual[ i ];
        this.solution_increment[ i ] -=
          residual / this.matrix[ i + this.stride ];
        break;
      }
      case SMOOTHER.GAUSS_SEIDEL_RED_BLACK : {
//      evens:
        for ( var i = 2; i <= this.ilast - 1; i += 2 ) {
          var residual = this.matrix[ i ] * this.solution_increment[ i - 1 ]
            + this.matrix[ i + this.stride ] * this.solution_increment[ i ]
            + this.matrix[ i + stride2 ] * this.solution_increment[ i + 1 ]
            - this.first_residual[ i ];
          this.solution_increment[ i ] -=
            residual / this.matrix[ i + this.stride ];
        }
        if ( this.ilast % 2 == 0 ) {
          i = this.ilast;
          residual = this.matrix[ i ] * this.solution_increment[ i - 1 ]
            + this.matrix[ i + this.stride ] * this.solution_increment[ i ]
            - this.first_residual[ i ];
          this.solution_increment[ i ] -=
            residual / this.matrix[ i + this.stride ];
        }
//      odds:
        var i = 1;
        var residual = 
            this.matrix[ i + this.stride ] * this.solution_increment[ i ]
          + this.matrix[ i + stride2 ] * this.solution_increment[ i + 1 ]
          - this.first_residual[ i ];
        this.solution_increment[ i ] -=
          residual / this.matrix[ i + this.stride ];
        for ( var i = 3; i <= this.ilast - 1; i += 2 ) {
          var residual = this.matrix[ i ] * this.solution_increment[ i - 1 ]
            + this.matrix[ i + this.stride ] * this.solution_increment[ i ]
            + this.matrix[ i + stride2 ] * this.solution_increment[ i + 1 ]
            - this.first_residual[ i ];
          this.solution_increment[ i ] -=
            residual / this.matrix[ i + this.stride ];
        }
        if ( this.ilast % 2 == 1 ) {
          i = this.ilast;
          residual = this.matrix[ i ] * this.solution_increment[ i - 1 ]
            + this.matrix[ i + this.stride ] * this.solution_increment[ i ]
            - this.first_residual[ i ];
          this.solution_increment[ i ] -=
            residual / this.matrix[ i + this.stride ];
        }
        break;
      }
      case SMOOTHER.RICHARDSON : {
        var i = 1;
        this.richardson_residual[ i ] = 
            this.matrix[ i + this.stride ] * this.solution_increment[ i ]
          + this.matrix[ i + stride2 ] * this.solution_increment[ i + 1 ]
          - this.first_residual[ i ];
        for ( var i = 2; i <= this.ilast - 1; i ++ ) {
          this.richardson_residual[ i ] =
              this.matrix[ i ] * this.solution_increment[ i - 1 ]
            + this.matrix[ i + this.stride ] * this.solution_increment[ i ]
            + this.matrix[ i + stride2 ] * this.solution_increment[ i + 1 ]
            - this.first_residual[ i ];
        }
        i = this.ilast;
        this.richardson_residual[ i ] =
            this.matrix[ i ] * this.solution_increment[ i - 1 ]
          + this.matrix[ i + this.stride ] * this.solution_increment[ i ]
          - this.first_residual[ i ];
        for ( var i = 1; i <= this.ilast; i ++ ) {
          this.solution_increment[ i ] -=
            this.richardson_residual[ i ] / mu;
        }
        break;
      }
    }
  }
  for ( var i = 0; i < this.stride; i ++ ) {
    this.solution_increment0[ i ] = this.solution_increment[ i ];
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Level.prototype.postSmooth = function() {
  var stride2 = this.stride * 2;
  for ( var it = 0; it < this.smoother_iterations; it ++ ) {
    switch ( this.smoother ) {
      case SMOOTHER.GAUSS_SEIDEL : { // reverse of preSmooth
        var i = this.ilast;
        residual = this.matrix[ i ] * this.solution_increment[ i - 1 ]
          + this.matrix[ i + this.stride ] * this.solution_increment[ i ]
          - this.first_residual[ i ];
        this.solution_increment[ i ] -=
          residual / this.matrix[ i + this.stride ];

        for ( var i = this.ilast-1; i >= 2; i -- ) {
          residual = this.matrix[ i ] * this.solution_increment[ i - 1 ]
            + this.matrix[ i + this.stride ] * this.solution_increment[ i ]
            + this.matrix[ i + stride2 ] * this.solution_increment[ i + 1 ]
            - this.first_residual[ i ];
          this.solution_increment[ i ] -=
            residual / this.matrix[ i + this.stride ];
        }
        i = 1;
        var residual = 
            this.matrix[ i + this.stride ] * this.solution_increment[ i ]
          + this.matrix[ i + stride2 ] * this.solution_increment[ i + 1 ]
          - this.first_residual[ i ];
        this.solution_increment[ i ] -=
          residual / this.matrix[ i + this.stride ];
        break;
      }
      case SMOOTHER.GAUSS_SEIDEL_RED_BLACK : {
//      odds:
        if ( this.ilast % 2 == 1 ) {
          i = this.ilast;
          residual = this.matrix[ i ] * this.solution_increment[ i - 1 ]
            + this.matrix[ i + this.stride ] * this.solution_increment[ i ]
            - this.first_residual[ i ];
          this.solution_increment[ i ] -=
            residual / this.matrix[ i + this.stride ];
        }
        for ( var i = this.ilast - 1; i >= 3; i -- ) {
          if ( i % 2 == 1 ) {
            var residual =
              this.matrix[ i ] * this.solution_increment[ i - 1 ]
            + this.matrix[ i + this.stride ] * this.solution_increment[ i ]
            + this.matrix[ i + stride2 ] * this.solution_increment[ i + 1 ]
            - this.first_residual[ i ];
            this.solution_increment[ i ] -=
              residual / this.matrix[ i + this.stride ];
          }
        }
        var i = 1;
        var residual = 
            this.matrix[ i + this.stride ] * this.solution_increment[ i ]
          + this.matrix[ i + stride2 ] * this.solution_increment[ i + 1 ]
          - this.first_residual[ i ];
        this.solution_increment[ i ] -=
          residual / this.matrix[ i + this.stride ];
//      evens:
        if ( this.ilast % 2 == 0 ) {
          i = this.ilast;
          residual = this.matrix[ i ] * this.solution_increment[ i - 1 ]
            + this.matrix[ i + this.stride ] * this.solution_increment[ i ]
            - this.first_residual[ i ];
          this.solution_increment[ i ] -=
            residual / this.matrix[ i + this.stride ];
        }
        for ( var i = this.ilast - 1; i >= 2; i -- ) {
          if ( i % 2 == 0 ) {
            var residual =
              this.matrix[ i ] * this.solution_increment[ i - 1 ]
            + this.matrix[ i + this.stride ] * this.solution_increment[ i ]
            + this.matrix[ i + stride2 ] * this.solution_increment[ i + 1 ]
            - this.first_residual[ i ];
            this.solution_increment[ i ] -=
              residual / this.matrix[ i + this.stride ];
          }
        }
        break;
      }
      case SMOOTHER.RICHARDSON : {
        var i = 1;
        this.richardson_residual[ i ] = 
            this.matrix[ i + this.stride ] * this.solution_increment[ i ]
          + this.matrix[ i + stride2 ] * this.solution_increment[ i + 1 ]
          - this.first_residual[ i ];
        for ( var i = 2; i <= this.ilast - 1; i ++ ) {
          this.richardson_residual[ i ] =
              this.matrix[ i ] * this.solution_increment[ i - 1 ]
            + this.matrix[ i + this.stride ] * this.solution_increment[ i ]
            + this.matrix[ i + stride2 ] * this.solution_increment[ i + 1 ]
            - this.first_residual[ i ];
        }
        i = this.ilast;
        this.richardson_residual[ i ] =
            this.matrix[ i ] * this.solution_increment[ i - 1 ]
          + this.matrix[ i + this.stride ] * this.solution_increment[ i ]
          - this.first_residual[ i ];
        for ( var i = 1; i <= this.ilast; i ++ ) {
          this.solution_increment[ i ] -=
            this.richardson_residual[ i ] / mu;
        }
        break;
      }
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Level.prototype.restriction = function() {
//document.getElementById("debug_textarea").value +=
//  "entering Level.restriction\n";
  for ( var i = 0; i < this.stride; i ++ ) {
    this.residual1[ i ] = this.Ax_b[ i ];
  }
  var stride2 = this.stride * 2;
  switch( this.restriction_prolongation ) {
    case RESTRICTION_PROLONGATION.ALGEBRAIC_MULTIGRID : {
      for ( var I = 1; I <= this.coarser.ilast; I ++ ) {
        var i = 2 * I;
        this.coarser.first_residual[ I ] = this.Ax_b[ i ]
          - this.Ax_b[ i - 1 ] * this.matrix[ i - 1 + stride2 ]
                               / this.matrix[ i - 1 + this.stride ]
          - this.Ax_b[ i + 1 ] * this.matrix[ i + 1 ]
                               / this.matrix[ i + 1 + this.stride ];
      }
      break;
    }
    case RESTRICTION_PROLONGATION.FINITE_ELEMENT : {
      for ( var I = 1; I <= this.coarser.ilast; I ++ ) {
        var i = 2 * I;
        this.coarser.first_residual[ I ] = this.Ax_b[ i ]
          + 0.5 * ( this.Ax_b[ i - 1 ] + this.Ax_b[ i + 1 ] );
      }
      break;
    }
  }
//document.getElementById("debug_textarea").value +=
//  "coarser.first_residual = ";
//for ( var I = 1; I <= this.coarser.ilast; I ++ ) {
//  document.getElementById("debug_textarea").value +=
//    this.coarser.first_residual[ I ] + " ";
//}
//document.getElementById("debug_textarea").value += "\n";
//document.getElementById("debug_textarea").value +=
//  "leaving Level.restriction\n";
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Level.prototype.prolongation = function() {
  var stride2 = this.stride * 2;
  switch( this.restriction_prolongation ) {
    case RESTRICTION_PROLONGATION.ALGEBRAIC_MULTIGRID : {
      for ( var I = 1; I <= this.coarser.ilast; I ++ ) {
        var i = 2 * I;
       this.prolongation_vector[ i ] = this.coarser.solution_increment[ I ];
      }
      var i = 1;
      this.prolongation_vector[ i ] =-
        this.prolongation_vector[ i + 1 ] * this.matrix[ i + stride2 ]
                                         / this.matrix[ i + this.stride ];
      for ( var I = 1; I <= this.coarser.ilast - 1; I ++ ) {
        var i = 2 * I + 1;
        this.prolongation_vector[ i ] =-
          ( this.prolongation_vector[ i - 1 ] * this.matrix[ i ]
          + this.prolongation_vector[ i + 1 ] * this.matrix[ i + stride2 ] )
          / this.matrix[ i + this.stride ];
      }
      var I = this.coarser.ilast;
      var i = 2 * I + 1;
      this.prolongation_vector[ i ] =-
        this.prolongation_vector[ i - 1 ] * this.matrix[ i ]
                                         / this.matrix[ i + this.stride ];
      break;
    }
    case RESTRICTION_PROLONGATION.FINITE_ELEMENT : {
      for ( var I = 1; I <= this.coarser.ilast; I ++ ) {
        var i = 2 * I;
        this.prolongation_vector[ i ] = this.coarser.solution_increment[ I ];
      }
      var i = 1;
      this.prolongation_vector[ i ] = 0.5 * this.prolongation_vector[ i + 1 ];
      for ( var I = 1; I <= this.coarser.ilast - 1; I ++ ) {
        var i = 2 * I + 1;
        this.prolongation_vector[ i ] =
          0.5 * ( this.prolongation_vector[ i - 1 ]
                + this.prolongation_vector[ i + 1 ] );
      }
      var I = this.coarser.ilast;
      var i = 2 * I + 1;
      this.prolongation_vector[ i ] = 0.5 * this.prolongation_vector[ i - 1 ];
      break;
    }
  }
  for ( var i = 1; i <= this.ilast; i ++ ) {
    this.solution_increment[ i ] += this.prolongation_vector[ i ];
    this.solution_increment1[ i ] = this.solution_increment[ i ];
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Level.prototype.computeResidual = function( x, b, Ax_b ) {
  var stride2 = this.stride * 2;
  var i = 1;
  Ax_b[ i ] = this.matrix[ i + this.stride ] * x[ i ]
            + this.matrix[ i + stride2 ] * x[ i + 1 ] - b[ i ];
  for ( var i = 2; i <= this.ilast - 1; i ++ ) {
    Ax_b[ i ] = this.matrix[ i ] * x[ i - 1 ]
              + this.matrix[ i + this.stride ] * x[ i ]
              + this.matrix[ i + stride2 ] * x[ i + 1 ] - b[ i ];
  }
  i = this.ilast;
  Ax_b[ i ] = this.matrix[ i ] * x[ i - 1 ]
            + this.matrix[ i + this.stride ] * x[ i ] - b[ i ];
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Level.prototype.updateResidual = function( ) {
  var stride2 = this.stride * 2;
  var i = 1;
  this.Ax_b[ i ] =
      this.matrix[ i + this.stride ] * this.solution_increment[ i ]
    + this.matrix[ i + stride2 ] * this.solution_increment[ i + 1 ];
  for ( var i = 2; i <= this.ilast - 1; i ++ ) {
    this.Ax_b[ i ] =
        this.matrix[ i ] * this.solution_increment[i - 1 ]
      + this.matrix[ i + this.stride ] * this.solution_increment[i ]
      + this.matrix[ i + stride2 ] * this.solution_increment[i + 1 ];
  }
  i = this.ilast;
  this.Ax_b[ i ] =
      this.matrix[ i ] * this.solution_increment[ i - 1 ]
    + this.matrix[ i + this.stride ] * this.solution_increment[ i ];
  for ( var i = 1; i <= this.ilast; i ++ ) {
    this.Ax_b[ i ] = this.first_residual[ i ] - this.Ax_b[ i ];
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Level.prototype.solve = function( ) {
//document.getElementById("debug_textarea").value +=
//  "entering Level.solve\n";
//document.getElementById("debug_textarea").value +=
//  "first_residual = ";
//for ( var i = 1; i <= this.ilast; i ++ ) {
//  document.getElementById("debug_textarea").value +=
//    this.first_residual[ i ] + " ";
//}
//document.getElementById("debug_textarea").value += "\n";

//forward-solve:
  this.solution_increment[ 1 ] = this.first_residual[ 1 ];
  for ( var i = 2; i <= this.ilast; i ++ ) {
    this.solution_increment[ i ] = this.first_residual[ i ]
      - this.matrix_copy[ i ] * this.solution_increment[ i - 1 ];
  }
//back-solve:
  this.solution_increment[ this.ilast ] /=
    this.matrix_copy[ this.ilast + this.stride ];
  var stride2 = this.stride * 2;
  for ( var i = this.ilast - 1; i >= 1; i -- ) {
    this.solution_increment[ i ] = ( this.solution_increment[ i ]
      - this.matrix_copy[ i + stride2 ] * this.solution_increment[ i + 1 ] )
      / this.matrix_copy[ i + stride ];
  }
  for ( var i = 1; i <= this.ilast; i ++ ) this.Ax_b[ i ] = 0.;

//document.getElementById("debug_textarea").value +=
//  "solution_increment = ";
//for ( var i = 1; i <= this.ilast; i ++ ) {
//  document.getElementById("debug_textarea").value +=
//    this.solution_increment[ i ] + " ";
//}
//document.getElementById("debug_textarea").value += "\n";
//document.getElementById("debug_textarea").value +=
//  "leaving Level.restriction\n";
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Level.prototype.multigridStep = function( Ax_b, d ) {
//document.getElementById("debug_textarea").value +=
//  "entering Level.multigridStep, level_number = " + this.level_number
//  + "\n";
  if ( this.finer == null ) {
    for ( var i = 1; i <= this.ilast; i ++ ) {
      this.first_residual[ i ] = Ax_b[ i ];
    }
  }
  for ( var i = 1; i <= this.ilast; i ++ ) this.solution_increment[ i ] = 0.;
  if ( this.coarser == null ) this.solve();
  else {
    this.preSmooth();
    this.updateResidual();
    this.restriction();
    this.coarser.multigridStep( this.coarser.Ax_b,
      this.coarser.solution_increment );
    this.prolongation();
    this.updateResidual();
    this.postSmooth();
  }
  if ( this.finer == null ) {
    for ( var i = 1; i <= this.ilast; i ++ ) {
      d[ i ] = this.solution_increment[ i ];
    }
  }
//document.getElementById("debug_textarea").value +=
//  "leaving Level.multigridStep\n";
}
