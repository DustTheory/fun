function Minpack() {
  this.epsmch = 2.*LaPack0.dlamch( 'E' );
  this.jeval = false;
  this.inside_inner_loop = false;
  this.all_done = false;
  this.delta = 0.;
  this.fnorm = 0.;
  this.xnorm = 0.;
  this.xtol = 0.;
  this.iter = 0;
  this.maxfev = 0;
  this.nfev = 0;
  this.njev = 0;
  this.msum = 0;
  this.ncsuc = 0;
  this.ncfail = 0;
  this.nslow1 = 0;
  this.nslow2 = 0;
  this.info = 0;
  this.ldfjac = 0;
  this.lr = 0;
  this.iflag = new IntReference( 0 );
  this.sing = new BooleanReference( false );
  this.fvec = new Array( );
  this.fjac = new Array( );
  this.diag = new Array( );
  this.r = new Array( );
  this.qtf = new Array( );
  this.wa1 = new Array( ); 
  this.wa2 = new Array( ); 
  this.wa3 = new Array( ); 
  this.wa4 = new Array( ); 
  this.graph_tool = undefined;
  this.soln = new Array( );
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Minpack.prototype.chkder = function( m, n, x, fvec, fjac, ldfjac, xp, fvecp,
mode, err, ioffx, iofffvec, iofffjac, ioffxp, ioffvecp, iofferr ) {
  var eps = Math.sqrt( this.epsmch );
  if ( mode !== 2 ) {
//  first call
    for ( var j = 1; j <= n; j ++ ) { // 10
      var temp = eps * Math.abs( x[ ioffx + j - 1 ] );
      if ( temp == 0. ) temp = eps;
      xp[ ioffxp + j - 1 ] = x[ ioffx + j - 1 ] + temp;
    }
  } else {
//  not first call
    var epsf = 1.e2 * this.epsmch;
    var epslog = Math.LOG10E * Math.log( eps );
    for ( var i = 1; i <= m; i ++ ) err[ iofferr + i - 1 ] = 0.;
    for ( var j = 1; j <= n; j ++ ) { // 50
      var temp = Math.abs( x[ ioffx + j - 1 ] );
      if ( temp == 0. ) temp = 1.;
      for ( var i = 1; i <= m; i ++ ) { // 40
        err[ iofferr + i - 1 ] +=
          temp * fjac[ iofffjac + i - 1 + ( j - 1 ) * ldfjac ];
      }
    }
    for ( var i = 1; i <= m; i ++ ) { // 60
      var temp = 1.;
      if ( fvec[ iofffvec + i - 1 ] !== 0. &&
      fvecp[ iofffvecp + i - 1 ] !== 0. &&
      Math.abs( fvecp[ iofffvecp + i - 1 ] - fvec[ iofffvec + i - 1 ] ) >=
      epsf * Math.abs( fvec[ iofffvec + i - 1 ] ) ) {
        temp = eps * Math.abs( ( fvecp[ iofffvecp + i - 1 ]
                               - fvec[ iofffvec + i - 1 ] ) / eps
                               - err[ iofferr + i - 1 ] )
             / ( Math.abs( fvec[ iofffvec + i - 1 ] )
               + Math.abs( fvecp[ iofffvecp + i - 1 ] ) );
        err[ iofferr + i - 1 ] = 1.;
        if ( temp > this.epsmch && temp < eps ) {
          err[ iofferr + i - 1 ] = ( Math.LOG10E * Math.log( temp )
                                   - epslog ) / epslog;
        }
        if ( temp >= eps ) err[ iofferr + i - 1 ] = 0.;
      }
    }
  }
}
//use Blas1.dnrm2( n, x, 1, ioffx )
//Minpack.enorm( n, x, ioffx ) { }
Minpack.prototype.dogleg = function( n, x, diag, ioffx, ioffdiag ) {
//document.getElementById("debug_textarea").value +=
//  "entering dogleg, delta = " + this.delta + "\n";
  if ( this.graph_tool != undefined ) { // draw trust region
//  document.getElementById("debug_textarea").value +=
//    "drawing trust region\n";
    var np = 256;
    var dtheta = 2. * Math.PI / np;
    var theta = 0.;
    this.graph_tool.setFgColor();
    this.graph_tool.beginDrawing();
      var pos = new Array( 2 );
      pos[ 0 ] = this.soln[ 0 ] + this.delta; 
      pos[ 1 ] = this.soln[ 1 ];
//    document.getElementById("debug_textarea").value +=
//      "pos = " + pos[0] + " " + pos[1] + "\n";
      this.graph_tool.movePen( pos );
      theta += dtheta;
      for ( var i = 1; i <= np; i ++ ) {
        pos[ 0 ] = this.soln[ 0 ] + Math.cos( theta ) * this.delta; 
        pos[ 1 ] = this.soln[ 1 ] + Math.sin( theta ) * this.delta;
//      document.getElementById("debug_textarea").value +=
//        "pos = " + pos[0] + " " + pos[1] + "\n";
        this.graph_tool.drawLine( pos );
        theta += dtheta;
      }
    this.graph_tool.endDrawing();
  }
//compute Gauss-Newton step
  var jj = ( n * ( n + 1 ) ) / 2 + 1;
  for ( var k = 1; k <= n; k ++ ) { // 50
    var j = n - k + 1;
    jj -= k;
    var l = jj + 1;
    var sum = 0.;
    for ( var i = j + 1; i <= n; i ++ ) { // 10
      sum += this.r[ l - 1 ] * x[ ioffx + i - 1 ];
      l ++;
    }
    var temp = this.r[ jj - 1 ];
    if ( temp == 0. ) {
      l = j;
      for ( var i = 1; i <= j; i ++ ) { // 30
        temp = Math.max( temp, Math.abs( this.r[ l - 1 ] ) );
        l += n - i;
      }
      temp *= this.epsmch;
      if ( temp == 0. ) temp = this.epsmch;
    }
    x[ ioffx + j - 1 ] = ( this.qtf[ j - 1 ] - sum ) / temp;
  }
  if ( this.graph_tool != undefined ) { // draw Gauss-Newton step in blue
    this.graph_tool.setFgColor( 0. ); // blue
    this.graph_tool.beginDrawing();
      var pos = new Array( 2 );
      pos[ 0 ] = this.soln[ 0 ] - x[ 0 ];
      pos[ 1 ] = this.soln[ 1 ] - x[ 1 ];
      this.graph_tool.movePen( this.soln );
      this.graph_tool.drawLine( pos );
    this.graph_tool.endDrawing();
  }
  for ( var j = 1; j <= n; j ++ ) { // 60
    this.wa2[ j - 1 ] = 0.;
    this.wa3[ j - 1 ] = diag[ ioffdiag + j - 1 ] * x[ ioffx + j - 1 ];
  }
  var qnorm = Blas1.dnrm2( n, this.wa3, 1, 0 );
  if ( qnorm > this.delta ) {
//  Gauss-Newton step unacceptable
//  compute scaled gradient step
    var l = 1;
    for ( var j = 1; j <= n; j ++ ) { // 80
      var temp = this.qtf[ j - 1 ];
      for ( var i = j; i <= n; i ++ ) { // 70
        this.wa2[ i - 1 ] += this.r[ l - 1 ] * temp;
        l ++;
      }
      this.wa2[ j - 1 ] /= diag[ ioffdiag + j - 1 ];
    }
    var gnorm = Blas1.dnrm2( n, this.wa2, 1, 0 );
    var sgnorm = 0.;
    var alpha = this.delta / qnorm;
    if ( gnorm !== 0. ) {
      for ( var j = 1; j <= n; j ++ ) { // 90
        this.wa2[ j - 1 ] =
          ( this.wa2[ j - 1 ] / gnorm ) / diag[ ioffdiag + j - 1 ];
      }
      l = 1;
      for ( var j = 1; j <= n; j ++ ) { // 110
        sum = 0.;
        for ( var i = j; i <= n; i ++ ) { // 100
          sum += this.r[ l - 1 ] * this.wa2[ i - 1 ];
          l ++;
        }
        this.wa3[ j - 1 ] = sum;
      }
      temp = Blas1.dnrm2( n, this.wa3, 1, 0 );
      sgnorm = ( gnorm / temp ) / temp;
      if ( this.graph_tool != undefined ) { // scaled gradient step in red
        this.graph_tool.setFgColor( 1. ); // red
        this.graph_tool.beginDrawing();
          var pos = new Array( 2 );
          pos[ 0 ] = this.soln[ 0 ] - sgnorm * this.wa2[ 0 ];
          pos[ 1 ] = this.soln[ 1 ] - sgnorm * this.wa2[ 1 ];
          this.graph_tool.movePen( this.soln );
          this.graph_tool.drawLine( pos );
        this.graph_tool.endDrawing();
      }
      alpha = 0.;
      if ( sgnorm < this.delta ) {
//      scaled gradient step unacceptable
//      compute dogleg step
        var bnorm = Blas1.dnrm2( n, this.qtf, 1, 0 );
        var temp = ( bnorm / gnorm ) * ( bnorm / qnorm )
          * ( sgnorm / this.delta );
        temp += - ( this.delta / qnorm ) * Math.pow( sgnorm / this.delta, 2 )
          + Math.sqrt( Math.pow( temp - ( this.delta / qnorm ), 2 )
                     + ( 1. - Math.pow( this.delta / qnorm, 2 ) )
                     * ( 1. - Math.pow( sgnorm / this.delta, 2 ) ) );
        alpha = ( ( this.delta / qnorm )
          * ( 1. - Math.pow( sgnorm / this.delta, 2 ) ) ) / temp;
        if ( this.graph_tool != undefined ) { // dogleg step in green
          this.graph_tool.setFgColor( 0.5 ); // green
          this.graph_tool.beginDrawing();
            var pos = new Array( 2 );
            var mytemp = ( 1. - alpha ) * sgnorm;
            pos[ 0 ] = this.soln[ 0 ]
                     - mytemp * this.wa2[ 0 ] - alpha * x[ ioffx ];
            pos[ 1 ] = this.soln[ 1 ]
                     - mytemp * this.wa2[ 1 ] - alpha * x[ ioffx + 1 ];
            this.graph_tool.movePen( this.soln );
            this.graph_tool.drawLine( pos );
          this.graph_tool.endDrawing();
        }
      }
    }
    temp = ( 1. - alpha ) * Math.min( sgnorm, this.delta );
    for ( var j = 1; j <= n; j ++ ) { // 130
      x[ ioffx + j - 1 ] = temp * this.wa2[ j - 1 ]
        + alpha * x[ ioffx + j - 1 ];
    }
  }
//document.getElementById("debug_textarea").value +=
//  "leaving dogleg\n";
}
Minpack.prototype.fdjac1 = function( fcn, n, x, fvec, ml, mu, epsfcn,
ioffx, iofffvec ) {
  var eps = Math.sqrt( Math.max( epsfcn, this.epsmch ) );
  this.msum = ml + mu + 1;
  if ( this.msum >= n ) {
//  dense approximate jacobian
    for ( var j = 1; j <= n; j ++ ) { // 20
      var temp = x[ ioffx + j - 1 ];
      var h = eps * Math.abs( temp );
      if ( h == 0. ) h = eps;
      x[ ioffx + j - 1 ] = temp + h;
      fcn( n, x, this.wa1, this.iflag, ioffx, 0 );
      if ( this.iflag.getValue() < 0 ) break;
      x[ ioffx + j - 1 ] = temp;
      for ( var i = 1; i <= n; i ++ ) { // 10
        this.fjac[ i - 1 + ( j - 1 ) * this.ldfjac ] =
          ( this.wa1[ i - 1 ] - fvec[ iofffvec + i - 1 ] ) / h;
      }
    }
  } else { // 40
//  banded approximate jacobian
    for ( var k = 1; k <= this.msum; k ++ ) { // 90
      for ( var j = k; j <= n; j += this.msum ) { // 60
        this.wa2[ j - 1 ] = x[ ioffx + j - 1 ];
        var h = eps * Math.abs( this.wa2[ j - 1 ] );
        if ( h == 0. ) h = eps;
        x[ ioffx + j - 1 ] = this.wa2[ j - 1 ] + h;
      }
      fcn( n, x, this.wa1, this.iflag, ioffx, 0 );
      if ( this.iflag.getValue() < 0 ) break;
      for ( var j = k; j <= n; j += this.msum ) {
        x[ ioffx + j - 1 ] = this.wa2[ j - 1 ];
        h = eps * Math.abs( this.wa2[ j - 1 ] );
        if ( h == 0 ) h = eps;
        for ( var i = 1; i <= n; i ++ ) {
          this.fjac[ i - 1 + ( j - 1 ) * this.ldfjac ] = 0.;
          if ( i >= j - mu && i <= j + ml ) {
            this.fjac[ i - 1 + ( j - 1 ) * this.ldfjac ] =
              ( this.wa1[ i - 1 ] - fvec[ iofffvec + i - 1 ] ) / h;
          }
        }
      }
    }
  }
}
Minpack.prototype.fdjac2 = function( fcn, m, n, x, fvec, epsfcn, wa,
ioffx, iofffvec, ioffwa ) {
  var eps = Math.sqrt( Math.max( epsfcn, this.epsmch ) );
  for ( var j = 1; j <= n; j ++ ) {
    var temp = x[ ioffx + j - 1 ];
    var h = eps * Math.abs( temp );
    if ( h == 0. ) h = eps;
    x[ ioffx + j - 1 ] = temp + h;
    fcn( m, n, x, wa, this.iflag, ioffx, ioffwa );
    if ( this.iflag.getValue() < 0 ) break;
    x[ ioffx + j - 1 ] = temp;
    for ( var i = 1; i <= m; i ++ ) {
      this.fjac[ i - 1 + ( j - 1 ) * this.ldfjac ] =
        ( wa[ ioffwa + i - 1 ] - fvec[ iofffvec + i - 1 ] ) / h;
    }
  }
}
// always called with q = this.fjac
Minpack.prototype.qform = function(  m, n, q, ldq, ioffq ) {
  var minmn = Math.min( m, n );
  for ( var j = 2; j <= minmn; j ++ ) { // 20
    for ( var i = 1; i < j; i ++ ) { // 10
      q[ ioffq + i - 1 + ( j - 1 ) * ldq ] = 0.;
    }
  }
  for ( var j = n + 1; j <= m; j ++ ) { // 50
    for ( var i = 1; i <= m; i ++ ) { // 40
      q[ ioffq + i - 1 + ( j - 1 ) * ldq ] = 0.;
    }
    q[ ioffq + j - 1 + ( j - 1 ) * ldq ] = 1.;
  }
  for ( var l = 1; l <= minmn; l ++ ) { // 120
    var k = minmn - l + 1;
    for ( var i = k; i <= m; i ++ ) { // 70
      this.wa1[ i - 1 ] = q[ ioffq + i - 1 + ( k - 1 ) * ldq ];
      q[ ioffq + i - 1 + ( k - 1 ) * ldq ] = 0;;
    }
    q[ ioffq + k - 1 + ( k - 1 ) * ldq ] = 1.;
    if ( this.wa1[ k - 1 ] !== 0. ) {
      for ( var j = k; j <= m; j ++ ) { // 100
        var sum = 0.;
        for ( var i = k; i <= m; i ++ ) { // 80
          sum += q[ ioffq + i - 1 + ( j - 1 ) * ldq ] * this.wa1[ i - 1 ];
        }
        var temp = sum / this.wa1[ k - 1 ];
        for ( var i = k; i <= m; i ++ ) { // 90
          q[ ioffq + i - 1 + ( j - 1 ) * ldq ] -= temp * this.wa1[ i - 1 ];
        }
      }
    }
  }
}
Minpack.prototype.qrfac = function( m, n, a, lda, pivot, ipvt, lipvt, rdiag,
acnorm, ioffa, ioffipvt, ioffrdiag, ioffacnorm ) {
  for ( var j = 1; j <= n; j ++ ) { // 10
    acnorm[ ioffacnorm + j - 1 ] = Blas1.dnrm2( m, a, 1, ioffa + j * lda );
    rdiag[ ioffrdiag + j - 1 ] = acnorm[ ioffacnorm + j - 1 ];
    this.wa3[ j - 1 ] = rdiag[ ioffrdiag + j - 1 ];
    if ( pivot ) ipvt[ ioffipvt + j - 1 ] = j;
  }
  var minmn = Math.min( m, n );
  for ( var j = 1; j <= minmn; j ++ ) { // 110
    if ( pivot ) {
      var kmax = j;
      for ( var k = j; k <= n; k ++ ) {
        if ( rdiag[ ioffrdiag + k - 1 ] > rdiag[ ioffrdiag + kmax - 1 ] ) {
          kmax = k;
        }
      }
      if ( kmax !== j ) {
        for ( var i = 1; i <= m; i ++ ) {
          var temp = a[ ioffa + i - 1 + ( j - 1 ) * lda ];
          a[ ioffa + i - 1 + ( j - 1 ) * lda ] =
            a[ ioffa + i - 1 + ( kmax - 1 ) * lda ];
          a[ ioffa + i - 1 + ( kmax - 1 ) * lda ] = temp;
        }
        rdiag[ ioffrdiag + kmax - 1 ] = rdiag[ ioffrdiag + j - 1 ];
        this.wa3[ kmax - 1 ] = this.wa3[ j - 1 ];
        var k = ipvt[ ioffipvt + j - 1 ];
        ipvt[ ioffipvt + j - 1 ] = ipvt[ ioffipvt + kmax - 1 ];
        ipvt[ ioffipvt + kmax - 1 ] = k;
      }
    }
    var ajnorm =
      Blas1.dnrm2( m - j + 1, a, 1, ioffa + j - 1 + ( j - 1 ) * lda );
    if ( ajnorm !== 0. ) {
      if ( a[ ioffa + j - 1 + ( j - 1 ) * lda ] < 0. ) ajnorm = - ajnorm;
      for ( var i = j; i <= m; i ++ ) {
        a[ ioffa + i - 1 + ( j - 1 ) * lda ] /= ajnorm;
      }
      a[ ioffa + j - 1 + ( j - 1 ) * lda ] += 1.;
      var jp1 = j + 1;
      for ( var k = jp1; k <= n; k ++ ) { // 90
        var sum = 0.;
        for ( var i = j; i <= m; i ++ ) { // 60
          sum += a[ ioffa + i - 1 + ( j - 1 ) * lda ] *
                 a[ ioffa + i - 1 + ( k - 1 ) * lda ];
        }
        var temp = sum / a[ ioffa + j - 1 + ( j - 1 ) * lda ]; 
        for ( var i = j; i <= m; i ++ ) {
          a[ ioffa + i - 1 + ( k - 1 ) * lda ] -=
            temp * a[ ioffa + i - 1 + ( j - 1 ) * lda ];
        }
        if ( pivot && rdiag[ ioffrdiag + k - 1 ] !== 0. ) {
          var temp = a[ ioffa + j - 1 + ( k - 1 ) * lda ]
                   / rdiag[ ioffrdiag + k - 1 ];
          rdiag[ ioffrdiag + k - 1 ] *=
            Math.sqrt( Math.max( 0., 1. - temp * temp ) );
          if ( 0.05 * Math.pow( rdiag[ ioffrdiag + k - 1 ]
          / this.wa3[ k - 1 ] , 2 ) <= this.epsmch ) {
            rdiag[ ioffrdiag + k - 1 ] =
              Blas1.dnrm2( m - j, a, 1, ioffa + jp1 - 1 + ( k - 1 ) * lda );
            this.wa3[ k - 1 ] = rdiag[ ioffrdiag + k - 1 ];
          }
        }
      }
    }
    rdiag[ ioffrdiag + j - 1 ] = - ajnorm;
  }
}
/*
Minpack.prototype.qrsolv = function( n, ldr, ipvt, qtb, x, sdiag, 
ioffqtb, ioffx, ioffsdiag ) {
  for ( var j = 1; j <= n; j ++ ) { // 20
    for ( var i = j; i <= n; i ++ ) { // 10
      this.r[ i - 1 + ( j - 1 ) * ldr ] =
        this.r[ j - 1 + ( i - 1 ) * ldr ]
    }
    x[ ioffx + jj - 1 ] = this.r[ j - 1 + ( j - 1 ) * ldr ];
    this.wa2[ j - 1 ] = qtb[ ioffqtb + j - 1 ];
  }
  for ( var j = 1; j <= n; j ++ ) { // 100
    var l = ipvt[ ioffipvt + j - 1 ];
    if ( this.diag[ ll - 1 ] !== 0. ) {
      for ( var k = j; k <= n; k ++ ) sdiag[ ioffsdiag + k - 1 ] = 0.;
      sdiag[ ioffsdiag + j - 1 ] = this.diag[ l - 1 ];
      var qtbpj = 0.;
      for ( var k = j; k <= n; k ++ ) { // 80
        if ( sdiag[ ioffsdiag + k - 1 ] !== 0. ) {
          if ( Math.abs( this.r[ k - 1 + ( k - 1 ) * ldr ] ) <
          Math.abs( sdiag[ ioffsdiag + k - 1 ] ) ) {
            var cotan = this.r[ k - 1 + ( k - 1 ) * ldr ]
              / sdiag[ ioffsdiag + k - 1 ];
            var sin = 0.5 / Math.sqrt( 0.25 + 0.25 * cotan * cotan );
            var cos = sin * cotan;
          } else {
            var tan = sdiag[ ioffsdiag + k - 1 ]
              / this.r[ k - 1 + ( k - 1 ) * ldr ];
            var cos = 0.5 / Math.sqrt( 0.25 + 0.25 * tan * tan );
            var sin = cos * tan;
          }
          this.r[ k - 1 + ( k - 1 ) * ldr ] =
            cos * this.r[ k - 1 + ( k - 1 ) * ldr ]
            + sin * sdiag[ ioffsdiag + k - 1 ];
          var temp = cos * this.wa2[ k - 1 ] + sin * qtbpj;
          qtbpj = - sin * this.wa2[ k - 1 ] + cos * qtbpj;
          this.wa2[ k - 1 ] = temp;
          for ( var i = k + 1; i <= n; i ++ ) { // 60
            temp = cos * this.r[ i - 1 + ( k - 1 ) * ldr ]
              + sin * sdiag[ ioffsdiag + i - 1 ];
            sdiag[ ioffsdiag + i - 1 ] =
              - sin * this.r[ i - 1 + ( k - 1 ) * ldr ]
              + cos * sdiag[ ioffsdiag + i - 1 ];
            this.r[ i - 1 + ( k - 1 ) * ldr ] = temp;
          }
        }
      }
    }
    sdiag[ ioffsdiag + j - 1 ] =
      this.r[ j - 1 + ( j - 1 ) * ldr ];
    this.r[ j - 1 + ( j - 1 ) * ldr ] = x[ ioffx + j - 1 ];
  }
  var nsing = n;
  for ( var j = 1; j <= n; j ++ ) { // 110
    if ( sdiag[ ioffsdiag + j - 1 ] == 0. && nsing == n ) nsing = j - 1;
    if ( nsing < n ) this.wa2[ j - 1 ] = 0.;
  }
  for ( var k = 1; k <= nsing; k ++ ) { // 140
    var j = nsing - k + 1;
    var sum = 0.;
    for ( var i = j + 1; i <= nsing; i ++ ) { // 120
      sum += this.r[ i - 1 + ( j - 1 ) * ldr ]
        * this.wa2[ i - 1 ];
    }
    this.wa2[ j - 1 ] =
      ( this.wa2[ j - 1 ] - sum )
      / sdiag[ ioffsdiag + j - 1 ];
  }
  for ( var j = 1; j <= n; j ++ ) { // 160
    var l = ipvt[ ioffipvt + j - 1 ];
    x[ ioffx + l - 1 ] = this.wa2[ j - 1 ];
  }
}
*/
Minpack.prototype.r1mpyq = function( m, n, a, lda, v, w,
ioffa, ioffv, ioffw ) {
//document.getElementById("debug_textarea").value +=
//  "entering r1mpyq, m,n,lda = " + m + " " + n + " " + lda + "\n";
//document.getElementById("debug_textarea").value +=
//  "ioffa,ioffv,ioffw = " + ioffa + " " + ioffv + " " + ioffw + "\n";
  var nm1 = n - 1;
  for ( var nmj = 1; nmj <= nm1; nmj ++ ) { // 20
    var j = n - nmj;
    if ( Math.abs( v[ ioffv + j - 1 ] ) > 1. ) {
      var cos = 1. / v[ ioffv + j - 1 ];
      var sin = Math.sqrt( 1. - cos * cos );
    } else {
      var sin = v[ ioffv + j - 1 ];
      var cos = Math.sqrt( 1. - sin * sin );
    }
    for ( var i = 1; i <= m; i ++ ) { // 10
      var temp = cos * a[ ioffa + i - 1 + ( j - 1 ) * lda ]
        - sin * a[ ioffa + i - 1 + ( n - 1 ) * lda ];
      a[ ioffa + i - 1 + ( n - 1 ) * lda ] =
        sin * a[ ioffa + i - 1 + ( j - 1 ) * lda ]
        + cos * a[ ioffa + i - 1 + ( n - 1 ) * lda ];
      a[ ioffa + i - 1 + ( j - 1 ) * lda ] = temp;
    }
  }
  for ( var j = 1; j <= nm1; j ++ ) { // 40
    if ( Math.abs( w[ ioffw + j - 1 ] ) > 1. ) {
      var cos = 1. / w[ ioffw + j - 1 ];
      var sin = Math.sqrt( 1. - cos * cos );
    } else {
      var sin = w[ ioffw + j - 1 ];
      var cos = Math.sqrt( 1. - sin * sin );
    }
    for ( var i = 1; i <= m; i ++ ) { // 30
      var temp = cos * a[ ioffa + i - 1 + ( j - 1 ) * lda ]
        + sin * a[ ioffa + i - 1 + ( n - 1 ) * lda ];
      a[ ioffa + i - 1 + ( n - 1 ) * lda ] =
        - sin * a[ ioffa + i - 1 + ( j - 1 ) * lda ]
        + cos * a[ ioffa + i - 1 + ( n - 1 ) * lda ];
      a[ ioffa + i - 1 + ( j - 1 ) * lda ] = temp;
    }
  }
//document.getElementById("debug_textarea").value +=
//  "leaving r1mpyq\n";
}
// sing is BooleanReference
Minpack.prototype.r1updt = function( m, n, s, ls, u, v, w,
ioffs, ioffu, ioffv, ioffw ) {
//document.getElementById("debug_textarea").value +=
//  "entering r1updt\n";
  var giant = Number.MAX_VALUE;
  var jj = Math.floor( ( n * ( 2 * m - n + 1 ) ) / 2 ) - ( m - n );
  var l = jj;
  for ( var i = n; i <= m; i ++ ) { // 10
    w[ ioffw + i - 1 ] = s[ ioffs + l - 1 ];
    l ++;
  }
  var nm1 = n - 1;
  for ( var nmj = 1; nmj <= nm1; nmj ++ ) { // 60
    var j = n - nmj;
    jj -= ( m - j + 1 );
    w[ ioffw + j - 1 ] = 0.;
    if ( v[ ioffv + j - 1 ] !== 0. ) {
      if ( Math.abs( v[ ioffv + n - 1 ] ) < Math.abs( v[ ioffv + j - 1 ] ) )
      {
        var cotan = v[ ioffv + n - 1 ] / v[ ioffv + j - 1 ];
        var sin = 0.5 / Math.sqrt( 0.25 + 0.25 * cotan * cotan );
        var cos = sin * cotan;
        var tau = 1.;
        if ( Math.abs( cos ) * giant > 1. ) tau = 1. / cos;
      } else {
        var tan = v[ ioffv + j - 1 ] / v[ ioffv + n - 1 ];
        var cos = 0.5 / Math.sqrt( 0.25 + 0.25 * tan * tan );
        var sin = cos * tan;
        var tau = sin;
      }
      v[ ioffv + n - 1 ] = sin * v[ ioffv + j - 1 ]
                         + cos * v[ ioffv + n - 1 ];
      v[ ioffv + j - 1 ] = tau;
      l = jj;
      for ( var i = j; i <= m; i ++ ) { // 40
        var temp = cos * s[ ioffs + l - 1 ] - sin * w[ ioffw + i - 1 ];
        w[ ioffw + i - 1 ] = sin * s[ ioffs + l - 1 ]
                           + cos * w[ ioffw + i - 1 ];
        s[ ioffs + l - 1 ] = temp;
        l ++;
      }
    }
  }
  for ( var i = 1; i <= m; i ++ ) {
    w[ ioffw + i - 1 ] += v[ ioffv + n - 1 ] * u[ ioffu + i - 1 ];
  }
  this.sing.setValue( false );
  for ( var j = 1; j <= nm1; j ++ ) { // 130
    if ( w[ ioffw + j - 1 ] !== 0. ) {
      if ( Math.abs( s[ ioffs + jj - 1 ] ) < Math.abs( w[ ioffw + j - 1 ] ) )
      {
        var cotan = s[ ioffs + jj - 1 ] / w[ ioffw + j - 1 ];
        var sin = 0.5 / Math.sqrt( 0.25 + 0.25 * cotan * cotan );
        var cos = sin * cotan;
        var tau = 1.;
        if ( Math.abs( cos ) * giant > 1. ) tau = 1. / cos;
      } else {
        var tan = w[ ioffw + j - 1 ] / s[ ioffs + jj - 1 ];
        var cos = 0.5 / Math.sqrt( 0.25 + 0.25 * tan * tan );
        var sin = cos * tan;
        var tau = sin;
      }
      var l = jj;
      for ( var i = j; i <= m; i ++ ) { // 110
        var temp = cos * s[ ioffs + l - 1 ] + sin * w[ ioffw + i - 1 ];
        w[ ioffw + i - 1 ] = - sin * s[ ioffs + l - 1 ]
                           + cos * w[ ioffw + i - 1 ];
        s[ ioffs + l - 1 ] = temp;
        l ++;
      }
      w[ ioffw + j - 1 ] = tau;
    }
    if ( s[ ioffs + jj - 1 ] == 0. ) this.sing.setValue( true );
    jj += ( m - j + 1 );
  }
  var l = jj;
  for ( var i = n; i <= m; i ++ ) { // 150
    s[ ioffs + l - 1 ] = w[ ioffw + i - 1 ];
    l ++;
  }
  if ( s[ ioffs + jj - 1 ] == 0. ) this.sing.setValue( true );
//document.getElementById("debug_textarea").value +=
//  "leaving r1updt\n";
}
/*
// alpha is a NumberReference
Minpack.prototype.rwupdt = function( n, w, b, alpha, cos, sin,
ioffw, ioffb, ioffcos, ioffsin ) {
  for ( var j = 1; j <= n; j ++ ) { // 60
    var rowj = w[ ioffw + j - 1 ];
    for ( var i = 1; i <= j - 1; i ++ ) { // 10
      var temp = cos[ ioffcos + i - 1 ]
               * this.r[ i - 1 + ( j - 1 ) * this.ldr ]
               + sin[ ioffsin + i - 1 ] * rowj;
      var rowj = - sin[ ioffsin + i - 1 ]
               * this.r[ i - 1 + ( j - 1 ) * this.ldr ]
               + cos[ ioffcos + i - 1 ] * rowj;
      this.r[ i - 1 + ( j - 1 ) * this.ldr ] = temp;
    }
    cos[ ioffcos + j - 1 ] = 1.;
    sin[ ioffsin + j - 1 ] = 0.;
    if ( rowj !== 0. ) {
      if ( Math.abs( this.r[ j - 1 + ( j - 1 ) * this.ldr ] ) <
      Math.abs( rowj ) ) {
        var cotan = this.r[ j - 1 + ( j - 1 ) * this.ldr ] / rowj;
        sin[ ioffsin + j - 1 ] =
          0.5 / Math.sqrt( 0.25 + 0.25 * cotan * cotan );
        cos[ ioffcos + j - 1 ] = sin[ ioffsin + j - 1 ] * cotan;
      } else {
        var tan = rowj / this.r[ j - 1 + ( j - 1 ) * this.ldr ];
        cos[ ioffcos + j - 1 ] =
          0.5 / Math.sqrt( 0.25 + 0.25 * tan * tan );
        sin[ ioffsin + j - 1 ] = cos[ ioffcos + j - 1 ] * tan;
      }
      this.r[ j - 1 + ( j - 1 ) * this.ldr ] =
        cos[ ioffcos + j - 1 ] * this.r[ j - 1 + ( j - 1 ) * this.ldr ]
        + sin[ ioffsin + j - 1 ] * rowj;
      var temp = cos[ ioffcos + j - 1 ] * b[ ioffb + j - 1 ]
        + sin[ ioffsin + j - 1 ] * alpha.getValue();
      alpha.setValue( - sin[ ioffsin + j - 1 ] * b[ ioffb + j - 1 ]
        + cos[ ioffcos + j - 1 ] * alpha.getValue() );
      b[ ioffb + j - 1 ] = temp;
    }
  }
}
*/
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Minpack.prototype.hybrdTerminate = function( fcn, n, x, nprint, ioffx ) {
  if ( this.iflag.getValue() < 0 ) this.info =  this.iflag.getValue();
  this.iflag.setValue( 0 );
  if ( nprint > 0 ) fcn( n, x, this.fvec, this.iflag, ioffx, 0 );
  this.inside_inner_loop = false;
  this.all_done = true;
}
Minpack.prototype.hybrdInnerLoop = function( fcn, n, x, diag, nprint,
ioffx, ioffdiag  ) {
  if ( this.all_done ) return;
  this.inside_inner_loop = true;
  if ( nprint > 0 ) {
    this.iflag.setValue( 0 );
    if ( ( this.iter - 1 ) % nprint == 0 ) {
      fcn( n, x, this.fvec, this.iflag, ioffx, 0 );
    }
    if ( this.iflag.getValue() < 0 ) {
      this.hybrdTerminate( fcn, n, x, nprint, ioffx );
      return;
    }
  }
  this.dogleg( n, this.wa1, diag, 0, ioffdiag );
  for ( var j = 1; j <= n; j ++ ) { // 200
    this.wa1[ j - 1 ] = - this.wa1[ j - 1 ];
    this.wa2[ j - 1 ] = x[ ioffx + j - 1 ] + this.wa1[ j - 1 ];
    this.wa3[ j - 1 ] = diag[ ioffdiag + j - 1 ] * this.wa1[ j - 1 ];
  }
  var pnorm = Blas1.dnrm2( n, this.wa3, 1, 0 );
  if ( this.iter == 1 ) this.delta = Math.min( this.delta, pnorm );
  this.iflag.setValue( 1 );
  fcn( n, this.wa2, this.wa4, this.iflag, 0, 0 );
  this.nfev ++;
  if ( this.iflag.getValue() < 0 ) {
    this.hybrdTerminate( fcn, n, x, nprint, ioffx );
    return;
  }
  var fnorm1 = Blas1.dnrm2( n, this.wa4, 1, 0 );
  var actred = -1.;
  if ( fnorm1 < this.fnorm ) {
    actred = 1. - Math.pow( fnorm1 / this.fnorm, 2 );
  }
  var l = 1;
  for ( var i = 1; i <= n; i ++ ) { // 220
    var sum = 0.;
    for ( var j = i; j <= n; j ++ ) { // 210
      sum += this.r[ l - 1 ] * this.wa1[ j - 1 ];
      l ++;
    }
    this.wa3[ i - 1 ] = this.qtf[ i - 1 ] * sum;
  }
  var temp = Blas1.dnrm2( n, this.wa3, 1, 0 );
  var prered = 0.;
  if ( temp < this.fnorm ) prered = 1. - Math.pow( temp / this.fnorm, 2 );
  var ratio = 0.;
  if ( prered > 0. ) ratio = actred / prered;
  if ( ratio < 0.1 ) {
    this.ncsuc = 0;
    this.ncfail ++;
    this.delta *= 0.5;
  } else {
    this.ncfail = 0;
    this.ncsuc ++;
    if ( ratio >= 0.5 || this.ncsuc > 1 ) {
      this.delta = Math.max( this.delta, pnorm / 0.5 );
    }
    if ( Math.abs( ratio - 1. ) <= 0.1 ) this.delta = pnorm / 0.5 ;
  }
  if ( ratio >= 1.e-4 ) {
    for ( var j = 1; j <= n; j ++ ) { // 250
      x[ ioffx + j - 1 ] = this.wa2[ j - 1 ];
      this.wa2[ j - 1 ] = diag[ ioffdiag + j - 1 ] * x[ ioffx + j - 1 ];
      this.fvec[ j - 1 ] = this.wa4[ j - 1 ];
    }
    this.xnorm = Blas1.dnrm2( n, this.wa2, 1, 0 );
    this.fnorm = fnorm1;
    this.iter ++;
  }
  this.nslow1 ++;
  if ( actred >= 1.e-3 ) this.nslow1 = 0;
  if ( this.jeval ) this.nslow2 ++;
  if ( actred >= 0.1 ) this.nslow2 = 0;
  if ( this.delta <= this.xtol * this.xnorm || this.fnorm == 0. ) {
    this.info = 1;
  }
  if ( this.info != 0 ) {
    this.hybrdTerminate( fcn, n, x, nprint, ioffx );
    return;
  }
  if ( this.nfev >= this.maxfev ) this.info = 2;
  if ( 0.1 * Math.max( 0.1 * this.delta, pnorm ) <=
  this.epsmch * this.xnorm ) {
    this.info = 3;
  }
  if ( this.nslow2 == 5 ) this.info = 4;
  if ( this.nslow1 == 10 ) this.info = 5;
  if ( this.info != 0 ) {
    this.hybrdTerminate( fcn, n, x, nprint, ioffx );
    return;
  }
  if ( this.ncfail == 2 ) {
    this.inside_inner_loop = false;
    return;
  }
  for ( var j = 1; j <= n; j ++ ) { // 280
    var sum = 0.;
    for ( var i = 1; i <= n; i ++ ) { // 270
      sum += this.fjac[ i - 1 + ( j - 1 ) * this.ldfjac ]
        * this.wa4[ i - 1 ];
    }
    this.wa2[ j - 1 ] = ( sum - this.wa3[ j - 1 ] ) / pnorm;
    this.wa1[ j - 1 ] = this.diag[ j - 1 ]
      * ( ( this.diag[ j - 1 ] * this.wa1[ j - 1 ] ) / pnorm );
    if ( ratio >= 1.e-4 ) this.qtf[ j - 1 ] = sum;
  }
  this.r1updt( n, n, this.r, this.lr, this.wa1, this.wa2, this.wa3,
    0, 0, 0, 0 );
  this.r1mpyq( n, n, this.fjac, this.ldfjac, this.wa2, this.wa3, 0, 0, 0 );
  this.r1mpyq( 1, n, this.qtf, 1, this.wa2, this.wa3, 0, 0, 0 );
  this.jeval = false;
}
Minpack.prototype.hybrdOuterLoop = function( fcn, n, x, ml, mu, epsfcn,
diag, mode, factor, nprint, ioffx, ioffdiag ) {
  if ( this.all_done ) return;
  if ( this.inside_inner_loop ) return;
  this.jeval = true;
  this.iflag.setValue( 2 );
  fdjac1( fcn, n, x, this.fvec, ml, mu, epsfcn, ioffx, 0 );
  this.nfev += this.msum;
  if ( this.iflag.getValue() < 0 ) {
    this.hybrdTerminate( fcn, n, x, nprint, ioffx );
    return;
  }
  var iwa = new Array( 1 );
  this.qrfac( n, n, this.fjac, this.ldfjac, false, iwa, 1, this.wa1,
    this.wa2, 0, 0, 0, 0 );
  if ( this.iter == 1 ) {
    if ( mode !== 2 ) {
      for ( var j = 1; j <= n; j ++ ) { // 40
        diag[ ioffdiag + j - 1 ] = this.wa2[ j - 1 ];
        if ( this.wa2[ j - 1 ] == 0. ) diag[ ioffdiag + j - 1 ] = 1.;
      }
    }
    for ( var j = 1; j <= n; j ++ ) { // 60
      this.wa3[ j - 1 ] = diag[ ioffdiag + j - 1 ] * x[ ioffx + j - 1 ];
    }
    this.xnorm = Blas1.dnrm2( n, this.wa3, 1, 0 );
    this.delta = factor * this.xnorm;
    if ( this.delta == 0. ) this.delta = factor;
  }
  for ( var i = 1; i <= n; i ++ ) { // 80
    this.qtf[ i - 1 ] = this.fvec[ i - 1 ];
  }
  for ( var j = 1; j <= n; j ++ ) { // 120
    if ( this.fjac[ j - 1 + ( j - 1 ) * this.ldfjac ] !== 0 ) {
      var sum = 0.;
      for ( var i = j; i <= n; i ++ ) { // 90
        sum += this.fjac[ i - 1 + ( j - 1 ) * this.ldfjac ]
          * this.qtf[ i - 1 ];
      }
      var temp = - sum / this.fjac[ j - 1 + ( j - 1 ) * this.ldfjac ];
      for ( var i = j; i <= n; i ++ ) { // 100
        this.qtf[ i - 1 ] +=
          this.fjac[ i - 1 + ( j - 1 ) * this.ldfjac ] * temp;
      }
    }
  }
  this.sing.setValue( false );
  for ( var j = 1; j <= n; j ++ ) { // 150
    var l = j;
    for ( var i = 1; i <= j - 1; i ++ ) { // 130
      this.r[ l - 1 ] = this.fjac[ i - 1 + ( j - 1 ) * this.ldfjac ];
      l += n - i;
    }
    this.r[ l - 1 ] = this.wa1[ j - 1 ];
    if ( this.wa1[ j - 1 ] == 0. ) this.sing.setValue( true );
  }
  this.qform( n, n, this.fjac, this.ldfjac, 0 );
  if ( mode !== 2 ) {
    for ( var j = 1; j <= n; j ++ ) { // 160
      diag[ ioffdiag + j - 1 ] =
        Math.max( diag[ ioffdiag + j - 1 ], this.wa2[ j - 1 ] );
    }
  }
//while ( true ) { // 180
    this.inside_inner_loop = true;
    this.hybrdInnerLoop( fcn, n, x, diag, nprint, ioffx, ioffdiag  );
//}
}
Minpack.prototype.hybrdStep = function( fcn, n, x, ml, mu, epsfcn, diag,
mode, factor, nprint, ioffx, ioffdiag ) {
  if ( this.all_done ) return;
  if ( this.inside_inner_loop ) {
    this.hybrdInnerLoop( fcn, n, x, diag, nprint, ioffx, ioffdiag  );
  } else {
    this.hybrdOuterLoop( fcn, n, x, ml, mu, epsfcn, diag, mode, factor,
      nprint, ioffx, ioffdiag );
  }
}
Minpack.prototype.hybrd = function( fcn, n, x, xtol, maxfev, ml, mu, epsfcn,
diag, mode, factor, nprint, ioffx, ioffdiag ) {
  this.all_done = false;
  this.inside_inner_loop = false;

  this.fvec = new Array( n );
  this.fjac = new Array( n * n );
  this.ldfjac = n;
  this.info = 0;
  this.wa1 = new Array( n ); 
  this.wa2 = new Array( n ); 
  this.wa3 = new Array( n ); 
  this.wa4 = new Array( n ); 
  this.lr = ( n * ( n + 1 ) ) / 2;
  this.r = new Array( lr );
  this.qtf = new Array( n );

  this.info = 0;
  this.iflag.setValue( 0 );
  this.nfev = 0;
  if ( n <= 0 || xtol < 0. || maxfev <= 0 || ml < 0 || mu < 0 ||
  factor <= 0. ) {
    this.hybrdTerminate( fcn, n, x, nprint, ioffx );
    return;
  }
  if ( mode == 2 ) {
    for ( var j = 1; j <= n; j ++ ) {
      if ( diag[ ioffdiag + j - 1 ] <= 0. ) {
        this.hybrdTerminate( fcn, n, x, nprint, ioffx );
        return;
      }
    }
  }
  this.xtol = xtol;
  this.maxfev = maxfev;
  this.iflag.setValue( 1 );
  fcn( n, x, this.fvec, this.iflag, ioffx, 0 );
  this.nfev = 1;
  if ( this.iflag.getValue() < 0 ) {
    this.hybrdTerminate( fcn, n, x, nprint, ioffx );
    return;
  }
  this.fnorm = Blas1.dnrm2( n, this.fvec, 1, 0 );
  this.msum = Math.min( ml + mu + 1, n );
  this.iter = 1;
  this.ncsuc = 0;
  this.ncfail = 0;
  this.nslow1 = 0;
  this.nslow2 = 0;
//while ( true ) { // 30
    this.hybrdOuterLoop( fcn, n, x, ml, mu, epsfcn, mode, factor, nprint,
      ioffx );
//}
//this.hybrdTerminate( fcn, n, x, nprint, ioffx );
  return this.info;
}
Minpack.prototype.hybrd1Step = function( fcn, n, x, ioffx ) {
  if ( this.all_done ) return;
  if ( this.inside_inner_loop ) {
    this.hybrdInnerLoop( fcn, n, x, this.diag, 0, ioffx, 0  );
  } else {
    this.hybrdOuterLoop( fcn, n, x, n - 1, n - 1, 0., this.diag, 2, 1.e2, 0,
      ioffx, 0 );
  }
}
Minpack.prototype.hybrd1 = function( fcn, n, x, tol, ioffx ) {
  if ( n <= 0 || tol < 0. ) return;
  this.diag = new Array( n );
  for ( var j = 1; j <= n; j ++ ) this.diag[ j - 1 ] = 1.;
  this.hybrd( fcn, n, x, tol, 200 * ( n + 1 ), n - 1, n - 1, 0., this.diag,
    2, 1.e2, 0, ioffx, 0 );
  if ( this.info == 5 ) this.info = 4;
  return this.info;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Minpack.prototype.hybrjTerminate = function( fcn, n, x, nprint, ioffx ) {
  if ( this.iflag.getValue() < 0 ) this.info =  this.iflag.getValue();
  this.iflag.setValue( 0 );
  if ( nprint > 0 ) {
    fcn( n, x, this.fvec, this.fjac, this.ldfjac, this.iflag,
      ioffx, 0, 0 );
  }
  this.inside_inner_loop = false;
  this.all_done = true;
}
Minpack.prototype.hybrjInnerLoop = function( fcn, n, x, diag, nprint,
ioffx, ioffdiag ) {
//document.getElementById("debug_textarea").value +=
//  "entering hybrjInnerLoop, fnorm = " + this.fnorm + "\n";
  if ( this.all_done ) return;
  if ( ! this.inside_inner_loop ) return;
  if ( nprint > 0 ) {
    this.iflag.setValue( 0 );
    if ( ( this.iter - 1 ) % nprint == 0 ) {
      fcn( n, x, this.fvec, this.fjac, this.ldfjac, this.iflag,
        ioffx, 0, 0 );
    }
    if ( this.iflag.getValue() < 0 ) {
      this.hybrjTerminate( fcn, n, x, nprint, ioffx );
      return;
    }
  }
  this.dogleg( n, this.wa1, diag, 0, ioffdiag );
//document.getElementById("debug_textarea").value +=
//  "after dogleg, step = " + this.wa1[ 0 ] + " " + this.wa1[ 1 ] + "\n";
  for ( var j = 1; j <= n; j ++ ) { // 200
    this.wa1[ j - 1 ] = - this.wa1[ j - 1 ];
    this.wa2[ j - 1 ] = x[ ioffx + j - 1 ] + this.wa1[ j - 1 ];
    this.wa3[ j - 1 ] = diag[ ioffdiag + j - 1 ] * this.wa1[ j - 1 ];
  }
//document.getElementById("debug_textarea").value +=
//  "wa3 = " + this.wa3[ 0 ] + " " + this.wa3[ 1 ] + "\n";
  var pnorm = Blas1.dnrm2( n, this.wa3, 1, 0 );
  if ( this.iter == 1 ) this.delta = Math.min( this.delta, pnorm );
  this.iflag.setValue( 1 );
  fcn( n, this.wa2, this.wa4, this.fjac, this.ldfjac, this.iflag, 0, 0, 0 );
  this.nfev ++;
  if ( this.iflag.getValue() < 0 ) {
    this.hybrjTerminate( fcn, n, x, nprint, ioffx );
    return;
  }
  var fnorm1 = Blas1.dnrm2( n, this.wa4, 1, 0 );
  var actred = -1.;
  if ( fnorm1 < this.fnorm ) {
    actred = 1. - Math.pow( fnorm1 / this.fnorm, 2 );
  }
//document.getElementById("debug_textarea").value +=
//  "pnorm,fnorm,fnorm1,actred = " + pnorm + " " + this.fnorm + " " + fnorm1
//  + " " + actred + "\n";
//document.getElementById("debug_textarea").value +=
//  "r = " + this.r[ 0 ] + " " + this.r[ 1 ] + " " + this.r[ 2 ] + "\n";
//document.getElementById("debug_textarea").value +=
//  "qtf = " + this.qtf[ 0 ] + " " + this.qtf[ 1 ] + "\n";
//document.getElementById("debug_textarea").value +=
//  "wa1 = " + this.wa1[ 0 ] + " " + this.wa1[ 1 ] + "\n";
  var l = 1;
  for ( var i = 1; i <= n; i ++ ) { // 220
    var sum = 0.;
    for ( var j = i; j <= n; j ++ ) { // 210
      sum += this.r[ l - 1 ] * this.wa1[ j - 1 ];
      l ++;
    }
    this.wa3[ i - 1 ] = this.qtf[ i - 1 ] + sum;
  }
  var temp = Blas1.dnrm2( n, this.wa3, 1, 0 );
//document.getElementById("debug_textarea").value +=
//  "wa3 = " + this.wa3[ 0 ] + " " + this.wa3[ 1 ] + "\n";
//document.getElementById("debug_textarea").value +=
//  "temp = " + temp + "\n";
  var prered = 0.;
  if ( temp < this.fnorm ) prered = 1. - Math.pow( temp / this.fnorm, 2 );
  var ratio = 0.;
  if ( prered > 0. ) ratio = actred / prered;
  document.getElementById("debug_textarea").value =
    "iteration number,actual objective reduction = " + this.iter + " "
    + actred + "\n";
  document.getElementById("debug_textarea").value +=
    "predicted objective reduction,ratio of actual to predicted = "
    + prered + " " + ratio + "\n";
  if ( ratio < 0.1 ) {
    this.ncsuc = 0;
    this.ncfail ++;
    this.delta *= 0.5;
    document.getElementById("debug_textarea").value +=
      "step failed because ratio = " + ratio + " is less than 0.1\n";
    document.getElementById("debug_textarea").value +=
      "trust region radius decreased to " + this.delta + "\n";
  } else {
    this.ncfail = 0;
    this.ncsuc ++;
    if ( ratio >= 0.5 || this.ncsuc > 1 ) {
      this.delta = Math.max( this.delta, pnorm / 0.5 );
    }
    if ( Math.abs( ratio - 1. ) <= 0.1 ) this.delta = pnorm / 0.5 ;
    document.getElementById("debug_textarea").value +=
      "step succeeded because ratio = " + ratio + " is at least 0.1\n";
    document.getElementById("debug_textarea").value +=
      "trust region radius increased to " + this.delta + "\n";
  }
  if ( ratio >= 1.e-4 ) {
    for ( var j = 1; j <= n; j ++ ) { // 250
      x[ ioffx + j - 1 ] = this.wa2[ j - 1 ];
      this.wa2[ j - 1 ] = diag[ ioffdiag + j - 1 ] * x[ ioffx + j - 1 ];
      this.fvec[ j - 1 ] = this.wa4[ j - 1 ];
    }
    document.getElementById("debug_textarea").value +=
      "new position = " + x[ 0 ] + " " + x[ 1 ] + "\n";
    this.xnorm = Blas1.dnrm2( n, this.wa2, 1, 0 );
    this.fnorm = fnorm1;
    this.iter ++;
  }
  this.nslow1 ++;
  if ( actred >= 1.e-3 ) this.nslow1 = 0;
  if ( this.jeval ) this.nslow2 ++;
  if ( actred >= 0.1 ) this.nslow2 = 0;
  if ( this.delta <= this.xtol * this.xnorm || this.fnorm == 0. ) {
    this.info = 1;
  }
//document.getElementById("debug_textarea").value +=
//  "info = " + this.info + "\n";
  if ( this.info != 0 ) {
    this.hybrjTerminate( fcn, n, x, nprint, ioffx );
    return;
  }
  if ( this.nfev >= this.maxfev ) this.info = 2;
  if ( 0.1 * Math.max( 0.1 * this.delta, pnorm ) <=
  this.epsmch * this.xnorm ) {
    this.info = 3;
  }
  if ( this.nslow2 == 5 ) this.info = 4;
  if ( this.nslow1 == 10 ) this.info = 5;
//document.getElementById("debug_textarea").value +=
//  "info = " + this.info + "\n";
  if ( this.info != 0 ) {
    this.hybrjTerminate( fcn, n, x, nprint, ioffx );
    return;
  }
  if ( this.ncfail == 2 ) {
    this.inside_inner_loop = false;
    document.getElementById("debug_textarea").value +=
      "update Jacobian\n";
    return;
  }
  for ( var j = 1; j <= n; j ++ ) { // 280
    var sum = 0.;
    for ( var i = 1; i <= n; i ++ ) { // 270
      sum +=
        this.fjac[ i - 1 + ( j - 1 ) * this.ldfjac ] * this.wa4[ i - 1 ];
    }
    this.wa2[ j - 1 ] = ( sum - this.wa3[ j - 1 ] ) / pnorm;
    this.wa1[ j - 1 ] = diag[ ioffdiag + j - 1 ]
      * ( ( diag[ ioffdiag + j - 1 ] * this.wa1[ j - 1 ] ) / pnorm );
    if ( ratio >= 1.e-4 ) this.qtf[ j - 1 ] = sum;
  }
//document.getElementById("debug_textarea").value +=
//  "calling r1updt\n";
  this.r1updt( n, n, this.r, this.lr, this.wa1, this.wa2, this.wa3,
    0, 0, 0, 0 );
  this.r1mpyq( n, n, this.fjac, this.ldfjac, this.wa2, this.wa3, 0, 0, 0 );
  this.r1mpyq( 1, n, this.qtf, 1, this.wa2, this.wa3, 0, 0, 0 );
//document.getElementById("debug_textarea").value +=
//  "called r1mpyq\n";
  this.jeval = false;
//document.getElementById("debug_textarea").value +=
//  "leaving hybrjInnerLoop\n";
}
Minpack.prototype.hybrjOuterLoop = function( fcn, n, x, diag, mode, factor,
nprint, ioffx, ioffdiag ) {
//document.getElementById("debug_textarea").value +=
//  "entering hybrjOuterLoop, fnorm = " + this.fnorm + "\n";
  if ( this.all_done ) return;
  if ( this.inside_inner_loop ) return;
  this.jeval = true;
  this.iflag.setValue( 2 );
  fcn( n, x, this.fvec, this.fjac, this.ldfjac, this.iflag, ioffx, 0, 0 );
  this.njev ++;
  if ( this.iflag.getValue() < 0 ) {
    this.hybrjTerminate( fcn, n, x, nprint, ioffx );
    return;
  }
  var iwa = new Array( 1 );
  this.qrfac( n, n, this.fjac, this.ldfjac, false, iwa, 1, this.wa1,
    this.wa2, 0, 0, 0, 0 );
  if ( this.iter == 1 ) {
    if ( mode !== 2 ) {
      for ( var j = 1; j <= n; j ++ ) { // 40
        diag[ ioffdiag + j - 1 ] = this.wa2[ j - 1 ];
        if ( this.wa2[ j - 1 ] == 0. ) diag[ ioffdiag + j - 1 ] = 1.;
      }
    }
    for ( var j = 1; j <= n; j ++ ) { // 60
      this.wa3[ j - 1 ] = diag[ ioffdiag + j - 1 ] * x[ ioffx + j - 1 ];
    }
    this.xnorm = Blas1.dnrm2( n, this.wa3, 1, 0 );
    this.delta = factor * this.xnorm;
    if ( this.delta == 0. ) this.delta = factor;
  }
  for ( var i = 1; i <= n; i ++ ) { // 80
    this.qtf[ i - 1 ] = this.fvec[ i - 1 ];
  }
  for ( var j = 1; j <= n; j ++ ) { // 120
    if ( this.fjac[ j - 1 + ( j - 1 ) * this.ldfjac ] !== 0 ) {
      var sum = 0.;
      for ( var i = j; i <= n; i ++ ) { // 90
        sum += this.fjac[ i - 1 + ( j - 1 ) * this.ldfjac ]
             * this.qtf[ i - 1 ];
      }
      var temp = - sum / this.fjac[ j - 1 + ( j - 1 ) * this.ldfjac ]
      for ( var i = j; i <= n; i ++ ) { // 100
        this.qtf[ i - 1 ] +=
          this.fjac[ i - 1 + ( j - 1 ) * this.ldfjac ] * temp;
      }
    }
  }
  this.sing.setValue( false );
  for ( var j = 1; j <= n; j ++ ) { // 150
    var l = j;
    for ( var i = 1; i <= j - 1; i ++ ) { // 130
      this.r[ l - 1 ] = this.fjac[ i - 1 + ( j - 1 ) * this.ldfjac ];
      l += n - i;
    }
    this.r[ l - 1 ] = this.wa1[ j - 1 ];
    if ( this.wa1[ j - 1 ] == 0. ) this.sing.setValue( true );
  }
  this.qform( n, n, this.fjac, this.ldfjac, 0 );
  if ( mode !== 2 ) {
    for ( var j = 1; j <= n; j ++ ) { // 160
      diag[ ioffdiag + j - 1 ] =
        Math.max( diag[ ioffdiag + j - 1 ], this.wa2[ j - 1 ] );
    }
  }
//while ( true ) { // 180
    this.inside_inner_loop = true;
    this.hybrjInnerLoop( fcn, n, x, diag, nprint, ioffx, ioffdiag );
//}
//document.getElementById("debug_textarea").value +=
//  "leaving hybrjOuterLoop\n";
}
Minpack.prototype.hybrjStep = function( fcn, n, x, diag, mode, factor,
nprint, ioffx, ioffdiag ) {
//document.getElementById("debug_textarea").value +=
//  "entering hybrjStep\n";
  if ( this.all_done ) {
    document.getElementById("debug_textarea").value = "all done\n";
    return;
  }
  if ( this.inside_inner_loop ) {
    this.hybrjInnerLoop( fcn, n, x, diag, nprint, ioffx, ioffdiag  );
  } else {
    this.hybrjOuterLoop( fcn, n, x, diag, mode, factor, nprint,
      ioffx, ioffdiag );
  }
//document.getElementById("debug_textarea").value +=
//  "leaving hybrjStep\n";
}
Minpack.prototype.hybrj = function( fcn, n, x, xtol, maxfev, diag, mode,
factor, nprint, ioffx, ioffdiag ) {
//document.getElementById("debug_textarea").value +=
//  "entering hybrj, epsmch = " + this.epsmch + "\n";
  this.all_done = false;
  this.inside_inner_loop = false;

  this.fvec = new Array( n );
  this.fjac = new Array( n * n );
  this.ldfjac = n;
  this.info = 0;
  this.lr = ( n * ( n + 1 ) ) / 2;
  this.r = new Array( this.lr );
  this.qtf = new Array( n );
  this.wa1 = new Array( n ); 
  this.wa2 = new Array( n ); 
  this.wa3 = new Array( n ); 
  this.wa4 = new Array( n ); 

  this.info = 0;
  this.iflag.setValue( 0 );
  this.nfev = 0;
  this.njev = 0;

  if ( n <= 0 || xtol < 0. || maxfev <= 0 || factor <= 0. ) {
    this.hybrjTerminate( fcn, n, x, nprint, ioffx );
    return;
  }
  if ( mode == 2 ) {
    for ( var j = 1; j <= n; j ++ ) { // 10
      if ( diag[ ioffdiag + j - 1 ] <= 0. ) {
        this.hybrjTerminate( fcn, n, x, nprint, ioffx );
        return;
      }
    }
  }
  this.xtol = xtol;
  this.maxfev = maxfev;
  this.iflag.setValue( 1 );
  fcn( n, x, this.fvec, this.fjac, this.ldfjac, this.iflag, ioffx, 0, 0 );
  this.nfev = 1;
  if ( this.iflag.getValue() < 0 ) {
    this.hybrjTerminate( fcn, n, x, nprint, ioffx );
    return;
  }
  this.fnorm = Blas1.dnrm2( n, this.fvec, 1, 0 );
//document.getElementById("debug_textarea").value +=
//  "fnorm = " + this.fnorm +  "\n";
  this.iter = 1;
  this.ncsuc = 0;
  this.ncfail = 0;
  this.nslow1 = 0;
  this.nslow2 = 0;
//while ( true ) { // 30
    this.hybrjOuterLoop( fcn, n, x, diag, mode, factor, nprint,
      ioffx, ioffdiag );
//}
//document.getElementById("debug_textarea").value +=
//  "leaving hybrj\n";
  return this.info;
}
Minpack.prototype.hybrj1Step = function( fcn, n, x, ioffx ) {
//document.getElementById("debug_textarea").value +=
//  "entering hybrj1Step\n";
  if ( this.all_done ) {
    document.getElementById("debug_textarea").value = "all done\n";
    return;
  }
  if ( this.inside_inner_loop ) {
    this.hybrjInnerLoop( fcn, n, x, this.diag, 0, ioffx, 0  );
  } else {
    this.hybrjOuterLoop( fcn, n, x, this.diag, 2, 1.e2, 0, ioffx, 0 );
  }
//document.getElementById("debug_textarea").value +=
//  "leaving hybrj1Step\n";
}
Minpack.prototype.hybrj1 = function( fcn, n, x, tol, ioffx ) {
//document.getElementById("debug_textarea").value +=
//  "entering hybrj1\n";
  if ( n <= 0 || tol < 0. ) return;
  this.diag = new Array( n );
  for ( var j = 1; j <= n; j ++ ) this.diag[ j - 1 ] = 1.;
//document.getElementById("debug_textarea").value +=
//  "diag = " + this.diag[ 0 ] + " " + this.diag[ 1 ] + "\n";
  this.hybrj( fcn, n, x, tol, 100 * ( n + 1 ), this.diag, 2, 1.e2, 0,
    ioffx, 0 );
  if ( this.info == 5 ) this.info = 4;
//document.getElementById("debug_textarea").value +=
//  "leaving hybrj1\n";
  return this.info;
}
