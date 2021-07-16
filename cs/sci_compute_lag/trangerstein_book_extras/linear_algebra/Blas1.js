function Blas1() {
  var idim=[0,1,2,4];
  var ninc=[1,2,-2];
  var dx=[.6,.1,-.5,.8,.9,-.3,-.4];
  var dy=[.5,-.9,.3,.7,-.6,.2,.8];
//0., .3, .21, .62
  for (var id=0;id<idim.length;id++) {
    var val=Blas1.ddot( idim[id], dx, 1, dy, 1, 0, 0 );
  }
//0., .3, -.007, .85
  for (var id=0;id<idim.length;id++) {
    var val=Blas1.ddot( idim[id], dx, 2, dy, -2, 0, 0 );
  }
//0., .3, -.79, -.74
  for (var id=0;id<idim.length;id++) {
    var val=Blas1.ddot( idim[id], dx, -2, dy, 1, 0, 0 );
  }
//0., .3, .33, 1.27
  for (var id=0;id<idim.length;id++) {
    var val=Blas1.ddot( idim[id], dx, -1, dy, -2, 0, 0 );
  }

//[]
//[.68]
//[.68,-.87]
//[.68,-.87,.15,.94]
  for (var id=0;id<idim.length;id++) {
    var dyy=new Array(dy.length);
    for (var i=0;i<dy.length;i++) dyy[i]=dy[i];
    Blas1.daxpy( idim[id], .3, dx, 1, dyy, 1, 0, 0 );
  }
//[]
//[.68]
//[.48,-.9,.35]
//[.98,.2,-.75,.7,.57,-.9,.38]
  for (var id=0;id<idim.length;id++) {
    var dyy=new Array(dy.length);
    for (var i=0;i<dy.length;i++) dyy[i]=dy[i];
    Blas1.daxpy( idim[id], .3, dx, 2, dyy, -2, 0, 0 );
  }
//[]
//[.68]
//[.35,.72]
//[.38,-.63,.15,.88]
  for (var id=0;id<idim.length;id++) {
    var dyy=new Array(dy.length);
    for (var i=0;i<dy.length;i++) dyy[i]=dy[i];
    Blas1.daxpy( idim[id], .3, dx, -2, dyy, 1, 0, 0 );
  }
//[]
//[.68]
//[.33,-.9,.68]
//[1.04,.2,-.75,.7,.33,-.9,.68]
  for (var id=0;id<idim.length;id++) {
    var dyy=new Array(dy.length);
    for (var i=0;i<dy.length;i++) dyy[i]=dy[i];
    Blas1.daxpy( idim[id], .3, dx, -1, dyy, -2, 0, 0 );
  }

  var da=new NumberReference(.3);
  var db=new NumberReference(.4);
  var c=new NumberReference();
  var s=new NumberReference();
  Blas1.drotg(da,db,c,s); // .5,5/3,.6,.8
  da.setValue(.4);
  db.setValue(.3);
  Blas1.drotg(da,db,c,s); // .5,.6,.8,.6
  da.setValue(-.3);
  db.setValue(.4);
  Blas1.drotg(da,db,c,s); // .5,-5/3,-.6,.8
  da.setValue(-.4);
  db.setValue(.3);
  Blas1.drotg(da,db,c,s); // -.5,-.6,.8,-.6
  da.setValue(-.3);
  db.setValue(-.4);
  Blas1.drotg(da,db,c,s); // -.5,5/3,.6,.8
  da.setValue(0.);
  db.setValue(0.);
  Blas1.drotg(da,db,c,s); // 0.,0.,1.,0.
  da.setValue(0.);
  db.setValue(1.);
  Blas1.drotg(da,db,c,s); // 1.,1.,0.,1.
  da.setValue(1.);
  db.setValue(0.);
  Blas1.drotg(da,db,c,s); // 1.,0.,1.,0.

//[] & []
//[.78] & [.04]
//[.78,-.46] & [.04,-.78]
//[.78,-.46,-.22,1.06] & [.04,-.78,.54,.08]
  for (var id=0;id<idim.length;id++) {
    var dxx=new Array(dx.length);
    for (var i=0;i<dx.length;i++) dxx[i]=dx[i];
    var dyy=new Array(dy.length);
    for (var i=0;i<dy.length;i++) dyy[i]=dy[i];
    Blas1.drot( idim[id], dx, 1, dyy, 1, .8, .6, 0, 0 );
  }
//[] & []
//[.78] & [.04]
//[.66,.1,-.1] & [-.12,-.9,.7]
//[.96,.1,-.76,.8,.9,-.3,-.02] & [.28,.2,-.18,.7,-.3,-.9,.64]
  for (var id=0;id<idim.length;id++) {
    var dxx=new Array(dx.length);
    for (var i=0;i<dx.length;i++) dxx[i]=dx[i];
    var dyy=new Array(dy.length);
    for (var i=0;i<dy.length;i++) dyy[i]=dy[i];
    Blas1.drot( idim[id], dx, 2, dyy, -2, .8, .6, 0, 0 );
  }
//[] & []
//[.78] & [.04]
//[-.1,.1,-.06] & [.7,-1.08]
//[-.-2,-.3,.18,.8,-.22,.1,.9] & [.64,-1.26,.54,.2]
  for (var id=0;id<idim.length;id++) {
    var dxx=new Array(dx.length);
    for (var i=0;i<dx.length;i++) dxx[i]=dx[i];
    var dyy=new Array(dy.length);
    for (var i=0;i<dy.length;i++) dyy[i]=dy[i];
    Blas1.drot( idim[id], dx, -2, dyy, 1, .8, .6, 0, 0 );
  }
//[] & []
//[.78] & [.04]
//[.26,.78] & [.18,-.9,.04]
//[1.12,-.76,.26,.78] & [.16,.2,-.18,.7,.18,-.9,.04]
  for (var id=0;id<idim.length;id++) {
    var dxx=new Array(dx.length);
    for (var i=0;i<dx.length;i++) dxx[i]=dx[i];
    var dyy=new Array(dy.length);
    for (var i=0;i<dy.length;i++) dyy[i]=dy[i];
    Blas1.drot( idim[id], dx, -1, dyy, -2, .8, .6, 0, 0 );
  }
//[1.,2.,3.,4.,5.] & [1.,2.,3.,4.,5.]
  {
    var dxx=new Array(5);
    for (var i=0;i<5;i++) dxx[i]=i;
    var dyy=new Array(5);
    for (var i=0;i<5;i++) dyy[i]=i;
    Blas1.drot( 5, dx, 1, dyy, 1, 1., 0., 0, 0 );
  }
//[1.,2.,3.,4.,5.] & [-1.,-2.,-3.,-4.,-5.]
  {
    var dxx=new Array(5);
    for (var i=0;i<5;i++) dxx[i]=i;
    var dyy=new Array(5);
    for (var i=0;i<5;i++) dyy[i]=i;
    Blas1.drot( 5, dx, 1, dyy, 1, 0., 1., 0, 0 );
  }
//[5.,4.,3.,2.,1.] & [-1.,-2.,-3.,-4.,-5.]
  {
    var dxx=new Array(5);
    for (var i=0;i<5;i++) dxx[i]=i;
    var dyy=new Array(5);
    for (var i=0;i<5;i++) dyy[i]=5-i;
    Blas1.drot( 5, dx, 1, dyy, -1, 0., 1., 0, 0 );
  }
//[5.,4.,3.,2.,1.] & [-5.,-4.,-3.,-2.,-1.]
  {
    var dxx=new Array(5);
    for (var i=0;i<5;i++) dxx[i]=5-i;
    var dyy=new Array(5);
    for (var i=0;i<5;i++) dyy[i]=5-i;
    Blas1.drot( 5, dx, -1, dyy, -1, 0., 1., 0, 0 );
  }
//[1.,3.,5.] & [-1.,2.,-2.,4.,-3.]
  {
    var dxx=new Array(5);
    for (var i=0;i<5;i++) dxx[i]=i;
    var dyy=new Array(5);
    for (var i=0;i<5;i++) dyy[i]=i;
    Blas1.drot( 3, dx, 1, dyy, 2, 0., 1., 0, 0 );
  }
//[1.,2.,3.,4.,5.] & [-5.,-4.,-3.,-2.,-1.]
  {
    var dxx=new Array(5);
    for (var i=0;i<5;i++) dxx[i]=5-i;
    var dyy=new Array(5);
    for (var i=0;i<5;i++) dyy[i]=i;
    Blas1.drot( 5, dx, -1, dyy, 1, 0., 1., 0, 0 );
  }
//[-1.,-2.,-3.,-4.,-5.] & [1.,2.,3.,4.,5.]
  {
    var dxx=new Array(5);
    for (var i=0;i<5;i++) dxx[i]=i;
    var dyy=new Array(5);
    for (var i=0;i<5;i++) dyy[i]=i;
    Blas1.drot( 5, dx, 1, dyy, 1, 0., -1., 0, 0 );
  }
//[-5.,-4.,-3.,-2.,-1.] & [1.,2.,3.,4.,5.]
  {
    var dxx=new Array(5);
    for (var i=0;i<5;i++) dxx[i]=i;
    var dyy=new Array(5);
    for (var i=0;i<5;i++) dyy[i]=5-i;
    Blas1.drot( 5, dx, 1, dyy, -1, 0., -1., 0, 0 );
  }
//[-5.,-4.,-3.,-2.,-1.] & [ 5.,4.,3.,2.,1.]
  {
    var dxx=new Array(5);
    for (var i=0;i<5;i++) dxx[i]=5-i;
    var dyy=new Array(5);
    for (var i=0;i<5;i++) dyy[i]=5-i;
    Blas1.drot( 5, dx, -1, dyy, -1, 0., -1., 0, 0 );
  }
//[-1.,-3.,-5.] & [1.,2.,2.,4.,3.]
  {
    var dxx=new Array(5);
    for (var i=0;i<5;i++) dxx[i]=i;
    var dyy=new Array(5);
    for (var i=0;i<5;i++) dyy[i]=i;
    Blas1.drot( 3, dx, 1, dyy, 2, 0., -1., 0, 0 );
  }
//[-1.,-2.,-3.,-4.,-5.] & [5.,4.,3.,2.,1.]
  {
    var dxx=new Array(5);
    for (var i=0;i<5;i++) dxx[i]=5-i;
    var dyy=new Array(5);
    for (var i=0;i<5;i++) dyy[i]=i;
    Blas1.drot( 5, dx, -1, dyy, 1, 0., -1., 0, 0 );
  }

//[]
//[.6]
//[.6,.1]
//[.6,.1,-.5,.8]
  for (var id=0;id<idim.length;id++) {
    var dxx=new Array(dx.length);
    for (var i=0;i<dx.length;i++) dxx[i]=dx[i];
    var dyy=new Array(dy.length);
    for (var i=0;i<dy.length;i++) dyy[i]=dy[i];
    Blas1.dcopy( idim[id], dxx, 1, dyy, 1, 0, 0 );
  }
//[]
//[.6]
//[..6,-.9,-.5]
//[.6,.2,-.5,.7,.9,-.9,-.4]
  for (var id=0;id<idim.length;id++) {
    var dyy=new Array(dy.length);
    for (var i=0;i<dy.length;i++) dyy[i]=dy[i];
    Blas1.dcopy( idim[id], dxx, 2, dyy, -2, 0, 0 );
  }
//[]
//[.6]
//[-.5,.6]
//[-.4,.9,-.5,.6]
  for (var id=0;id<idim.length;id++) {
    var dyy=new Array(dy.length);
    for (var i=0;i<dy.length;i++) dyy[i]=dy[i];
    Blas1.dcopy( idim[id], dxx, -2, dyy, 1, 0, 0 );
  }
//[]
//[.6]
//[.1,-.9,.6]
//[.8,.2,-.5,.7,.1,-.9,.6]
  for (var id=0;id<idim.length;id++) {
    var dyy=new Array(dy.length);
    for (var i=0;i<dy.length;i++) dyy[i]=dy[i];
    Blas1.dcopy( idim[id], dxx, -1, dyy, -2, 0, 0 );
  }

//[] & []
//[.5] & [.6]
//[.5,-.9] & [.6,.1]
//[.5,-.9,.3,.7] & [.6,.1,-.5,.8]
  for (var id=0;id<idim.length;id++) {
    var dxx=new Array(dx.length);
    for (var i=0;i<dx.length;i++) dxx[i]=dx[i];
    var dyy=new Array(dy.length);
    for (var i=0;i<dy.length;i++) dyy[i]=dy[i];
    Blas1.dswap( idim[id], dxx, 1, dyy, 1, 0, 0 );
  }
//[] & []
//[.5] & [.6]
//[.3,.1,.5] & [.6,-.9,-.5]
//[.8,.1,-.6,.8,.3,-.3,.5] & [.6,.2,-.5,.7,.9,-.9]
  for (var id=0;id<idim.length;id++) {
    var dyy=new Array(dy.length);
    for (var i=0;i<dy.length;i++) dyy[i]=dy[i];
    Blas1.dswap( idim[id], dxx, 2, dyy, -2, 0, 0 );
  }
//[] & []
//[.5] & [.6]
//[.5,.1,-.9] & [-.5,.6]
//[.5,-.3,-.9,.8,.3,.1,.7] & [-.4,.9,-.5,.6]
  for (var id=0;id<idim.length;id++) {
    var dyy=new Array(dy.length);
    for (var i=0;i<dy.length;i++) dyy[i]=dy[i];
    Blas1.dswap( idim[id], dxx, -2, dyy, 1, 0, 0 );
  }
//[] & []
//[.5] & [.6]
//[.3,.5] & [.1,-.9,.6]
//[.8,-.6,.3,.5] & [.8,.2,-.5,.7,.1,-.9,.6]
  for (var id=0;id<idim.length;id++) {
    var dyy=new Array(dy.length);
    for (var i=0;i<dy.length;i++) dyy[i]=dy[i];
    Blas1.dswap( idim[id], dxx, -1, dyy, -2, 0, 0 );
  }

//0.
  {
    var x=[];
    var val=Blas1.dnrm2( 0, x, 1, 0 );
  }
//.3
  {
    var x=[.3];
    var val=Blas1.dnrm2( 1, x, 1, 0 );
  }
//.5
  {
    var x=[.3,-.4];
    var val=Blas1.dnrm2( 2, x, 1, 0 );
  }
//.7
  {
    var x=[.2,-.6,.3];
    var val=Blas1.dnrm2( 3, x, 1, 0 );
  }
//.6
  {
    var x=[.1,-.3,.5,-.1];
    var val=Blas1.dnrm2( 4, x, 1, 0 );
  }
//0.
  {
    var x=[];
    var val=Blas1.dnrm2( 0, x, 2, 0 );
  }
//.3
  {
    var x=[.3];
    var val=Blas1.dnrm2( 1, x, 2, 0 );
  }
//.5
  {
    var x=[.3,.2,-.4];
    var val=Blas1.dnrm2( 2, x, 2, 0 );
  }
//.7
  {
    var x=[.2,3.,-.6,5.,.3];
    var val=Blas1.dnrm2( 3, x, 2, 0 );
  }
//.6
  {
    var x=[.1,4.,-.3,6.,-.5,7.,-.1];
    var val=Blas1.dnrm2( 4, x, 2, 0 );
  }

//0.
  {
    var x=[];
    var val=Blas1.dasum( 0, x, 1, 0 );
  }
//.3
  {
    var x=[.3];
    var val=Blas1.dasum( 1, x, 1, 0 );
  }
//.7
  {
    var x=[.3,-.4];
    var val=Blas1.dasum( 2, x, 1, 0 );
  }
//1.1
  {
    var x=[.2,-.6,.3];
    var val=Blas1.dasum( 3, x, 1, 0 );
  }
//1.
  {
    var x=[.1,-.3,.5,-.1];
    var val=Blas1.dasum( 4, x, 1, 0 );
  }
//0.
  {
    var x=[];
    var val=Blas1.dasum( 0, x, 2, 0 );
  }
//.3
  {
    var x=[.3];
    var val=Blas1.dasum( 1, x, 2, 0 );
  }
//.7
  {
    var x=[.3,.2,-.4];
    var val=Blas1.dasum( 2, x, 2, 0 );
  }
//1.1
  {
    var x=[.2,3.,-.6,5.,.3];
    var val=Blas1.dasum( 3, x, 2, 0 );
  }
//1.
  {
    var x=[.1,4.,-.3,6.,-.5,7.,-.1];
    var val=Blas1.dasum( 4, x, 2, 0 );
  }

//[].
  {
    var x=[];
    Blas1.dscal( 0, x, 1, .3, 0 );
  }
//[-.3]
  {
    var x=[.3];
    Blas1.dscal( 1, x, 1, -1., 0 );
  }
//[0.,0.]
  {
    var x=[.3,-.4];
    Blas1.dscal( 2, x, 1, 0., 0 );
  }
//[.2,-.6,.3]
  {
    var x=[.2,-.6,.3];
    Blas1.dscal( 3, x, 1, 1., 0 );
  }
//[.03,-.09,.15,-.03]
  {
    var x=[.1,-.3,.5,-.1];
    Blas1.dscal( 4, x, 1, .3, 0 );
  }
//[]
  {
    var x=[];
    Blas1.dscal( 0, x, 2, .3, 0 );
  }
//[.09]
  {
    var x=[.3];
    Blas1.dscal( 1, x, 2, .3, 0 );
  }
//[.-9,2.,-.12]
  {
    var x=[.3,.2,-.4];
    Blas1.dscal( 2, x, 2, .3, 0 );
  }
//[.06,3.,-.18,5.,.09]
  {
    var x=[.2,3.,-.6,5.,.3];
    Blas1.dscal( 3, x, 2, .3, 0 );
  }
//[.03,4.,-.09,6.,-.15,7.,-.03]
  {
    var x=[.1,4.,-.3,6.,-.5,7.,-.1];
    Blas1.dscal( 4, x, 2, .3, 0 );
  }

//0
  {
    var x=[];
    var val=Blas1.idamax( 0, x, 1, .3, 0 );
  }
//1
  {
    var x=[.3];
    var val=Blas1.idamax( 1, x, 1, -1., 0 );
  }
//2
  {
    var x=[.3,-.4];
    var val=Blas1.idamax( 2, x, 1, 0., 0 );
  }
//2
  {
    var x=[.2,-.6,.3];
    var val=Blas1.idamax( 3, x, 1, 1., 0 );
  }
//3
  {
    var x=[.1,-.3,.5,-.1];
    var val=Blas1.idamax( 4, x, 1, .3, 0 );
  }
//0
  {
    var x=[];
    var val=Blas1.idamax( 0, x, 2, .3, 0 );
  }
//1
  {
    var x=[.3];
    var val=Blas1.idamax( 1, x, 2, .3, 0 );
  }
//2
  {
    var x=[.3,.2,-.4];
    var val=Blas1.idamax( 2, x, 2, .3, 0 );
  }
//2
  {
    var x=[.2,3.,-.6,5.,.3];
    var val=Blas1.idamax( 3, x, 2, .3, 0 );
  }
//3
  {
    var x=[.1,4.,-.3,6.,-.5,7.,-.1];
    var val=Blas1.idamax( 4, x, 2, .3, 0 );
  }
}
//*************************************************************************
//  [  c s ] [ a ] = [ rho ]
//  [ -s c ] [ b ] = [  0  ]
Blas1.drotg = function( a, b, c, s ) {
  var roe = b.getValue();
  if ( Math.abs(a.getValue()) > Math.abs(b.getValue()) ) {
    roe = a.getValue();
  }
  var scale = Math.abs(a.getValue()) + Math.abs(b.getValue());
  var r = Number.POSITIVE_INFINITY;
  var z = Number.POSITIVE_INFINITY;
  if ( scale <= 0 ) {
    c.setValue( 1 );
    s.setValue( 0 );
    r=0;
    z=0;
  } else {
    r = scale * Math.sqrt( Math.pow(a.getValue()/scale,2)
                         + Math.pow(b.getValue()/scale,2) );
    if ( roe < 0 ) r = -r;
    c.setValue( a.getValue() / r );
    s.setValue( b.getValue() / r );
    z = 1;
    if ( Math.abs(a.getValue()) > Math.abs(b.getValue()) ) {
      z = s.getValue();
    }
    if ( Math.abs(b.getValue()) >= Math.abs(a.getValue()) &&
    Math.abs(c.getValue()) > 0) {
      z = 1 / c.getValue();
    }
  }
  a.setValue( r );
  b.setValue( z );
}
Blas1.zrotg = function( a, b, c, s ) {
  var abs_a = ComplexMath.abs(a);
  if ( abs_a == 0 ) {
    c.setValue( 0 );
    s.setValue( new Complex( 1 , 0 ) );
    a.setValue( b );
  } else {
    var scale = abs_a + ComplexMath.abs(b);
    var norm = scale
      * Math.sqrt( Math.pow( abs_a / scale , 2 )
                 + Math.pow( ComplexMath.abs(b) / scale , 2 ) );
    var alpha = ComplexMath.divideNumber( a , abs_a );
    c.setValue( abs_a / norm );
    s.setValue( ComplexMath.divideNumber(
      ComplexMath.times( alpha, ComplexMath.conj(b) ) , norm ) );
    a.setValue( ComplexMath.timesNumber( alpha, norm ) );
  }
}
//*************************************************************************
//  [ x^T ] <-- [ h11 h12 ] [ x^T ]
//  [ y^T ] <-- [ h21 h22 ] [ y^T ]
//  where param[ 0 ] == -2 ==> H = identity
//  where param[ 0 ] == -1 ==> H given
//  where param[ 0 ] ==  0 ==>
//    [ h11 h12 ] = [  1  h12 ]
//    [ h21 h22 ] = [ h21  1  ]
//  where param[ 0 ] ==  1 ==>
//    [ h11 h12 ] = [ h11  1  ]
//    [ h21 h22 ] = [ -1  h22 ]
Blas1.drotm = function( n, x, incx, y, incy, param, ioffx, ioffy,
ioffparam ) {
  var dflag = param[ ioffparam ];
  if ( n <= 0 || ( dflag + 2 == 0. ) ) return;
  var kx = 1;
  var ky = 1;
  if ( incx < 0 ) kx = 1 + ( 1 - n ) * incx;
  if ( incy < 0 ) ky = 1 + ( 1 - n ) * incy;
  if ( dflag < 0 ) {
    var dh11 = param[ ioffparam + 1 ];
    var dh12 = param[ ioffparam + 3 ];
    var dh21 = param[ ioffparam + 2 ];
    var dh22 = param[ ioffparam + 4 ];
    for ( var i = 1; i <= n; i ++ ) {
      var w = x[ ioffx + kx - 1 ];
      var z = y[ ioffy + ky - 1 ];
      x[ ioffx + kx - 1 ] = w * dh11 + z * dh12;
      y[ ioffy + ky - 1 ] = w * dh21 + z * dh22;
      kx += incx;
      ky += incy;
    }
  } else if ( dflag == 0 ) {
    dh12 = param[ ioffparam + 3 ];
    dh21 = param[ ioffparam + 2 ];
    for ( i = 1; i <= n; i ++ ) {
      w = x[ ioffx + kx - 1 ];
      z = y[ ioffy + ky - 1 ];
      x[ ioffx + kx - 1 ] = w + z * dh12;
      y[ ioffy + ky - 1 ] = w * dh21 + z;
      kx += incx;
      ky += incy;
    }
  } else {
    dh11 = param[ ioffparam + 1 ];
    dh22 = param[ ioffparam + 4 ];
    for ( i = 1; i <= n; i ++ ) {
      w = x[ ioffx + kx - 1 ];
      z = y[ ioffy + ky - 1 ];
      x[ ioffx + kx - 1 ] = w * dh11 + z;
      y[ ioffy + ky - 1 ] = - w + z * dh22;
      kx += incx;
      ky += incy;
    }
  }
}
//*************************************************************************
//  construct modified Givens transformation
Blas1.drotmg = function( d1, d2, x1, y1, param, ioffparam) {
  var gam = 4096;
  var gamsq = 16777216;
  var rgamsq = 5.9604645e-8;
  if ( d1.getValue() < 0 ) {
    var flag = -1;
    var h11 = 0;
    var h12 = 0;
    var h21 = 0;
    var h22 = 0;
    d1.setValue( 0 );
    d2.setValue( 0 );
    x1.setValue( 0 );
  } else {
    var p2 = d2.getValue() * y1.getValue();
    if ( p2 == 0 ) {
      flag = -2;
      param[ ioffparam ] = flag;
      return;
    }
    var p1 = d1.getValue() * x1.getValue();
    var q2 = p2 * y1.getValue();
    var q1 = p1 * x1.getValue();
    if ( Math.abs( q1 ) > Math.abs( q2 ) ) {
      h21 = - y1.getValue() / x1.getValue();
      h12 = p2 / p1;
      var u = 1 - h12 * h21;
      if ( u > 0 ) {
        flag = 0;
        d1.setValue( d1.getValue() / u );
        d2.setValue( d2.getValue() / u );
        x1.setValue( x1.getValue() * u );
      }
    } else {
      if ( q2 < 0 ) {
        flag = -1;
        h11 = 0;
        h12 = 0;
        h21 = 0;
        h22 = 0;
        d1.setValue( 0 );
        d2.setValue( 0 );
        x1.setValue( 0 );
      } else {
        flag = 1;
        h11 = p1 / p2;
        h22 = x1.getValue() / y1.getValue();
        u = 1 + h11 * h22;
        var temp = d2.getValue() / u;
        d2.setValue( d1.getValue() / u );
        d1.setValue( temp );
        x1.setValue( y1.getValue() * u );
      }
    }
    if ( d1.getValue() != 0 ) {
      while ( d1.getValue() <= rgamsq || d1.getValue() >= gamsq ) {
        if ( flag == 0 ) {
          h11 = 1;
          h22 = 1;
          flag = -1;
        } else {
          h21 = -1;
          h12 = 1;
          flag = -1;
        }
        if ( d1.getValue() <= rgamsq ) {
          d1.setValue( d1.getValue() * ( gam * gam ) );
          x1.setValue( x1.getValue() / gam );
          h11 /= gam;
          h12 /= gam;
        } else {
          d1.setValue( d1.getValue() / ( gam * gam ) );
          x1.setValue( x1.getValue() * gam );
          h11 *= gam;
          h12 *= gam;
        }
      }
    }
    if ( d2.getValue() != 0 ) {
      while ( d2.getValue() <= rgamsq || d2.getValue() >= gamsq ) {
        if ( flag == 0 ) {
          h11 = 1;
          h22 = 1;
          flag = -1;
        } else {
          h21 = -1;
          h12 = 1;
          flag = -1;
        }
        if ( Math.abs( d2.getValue() ) <= rgamsq ) {
          d2.setValue( d2.getValue() * ( gam * gam ) );
          h21 /= gam;
          h22 /= gam;
        } else {
          d2.setValue( d2.getValue() / ( gam * gam ) );
          h21 *= gam;
          h22 *= gam;
        }
      }
    }
  }
  if ( flag < 0. ) {
    param[ ioffparam + 1 ] = h11;
    param[ ioffparam + 2 ] = h21;
    param[ ioffparam + 3 ] = h12;
    param[ ioffparam + 4 ] = h22;
  } else if ( flag == 0. ) {
    param[ ioffparam + 2 ] = h21;
    param[ ioffparam + 3 ] = h12;
  } else {
    param[ ioffparam + 1 ] = h11;
    param[ ioffparam + 4 ] = h22;
  }
  param[ ioffparam ] = flag;
}
//*************************************************************************
//  _nrm2 = || x ||_2
Blas1.dnrm2 = function( n, x, incx , ioffx ) {
  var norm = Number.POSITIVE_INFINITY;
  if ( n < 1 || incx < 1) {
    norm = 0;
  } else if ( n == 1 ) {
    norm = Math.abs( x[ ioffx ] );
  } else {
    var scale = 0;
    var ssq = 1;
    for (var ix = 1; ix <= 1 + ( n - 1 ) * incx; ix += incx ) {
      if ( x[ ioffx + ix - 1 ] != 0 ) {
        var absxi = Math.abs( x[ ioffx + ix - 1 ] );
        if ( scale < absxi ) {
          ssq = 1 + ssq * Math.pow( scale / absxi , 2 );
          scale = absxi;
        } else {
          ssq += Math.pow( absxi / scale , 2 );
        }
      }
    }
    norm = scale * Math.sqrt( ssq );
  }
  return norm;
}
Blas1.dznrm2 = function( n, x, incx , ioffx ) {
  var norm = Number.POSITIVE_INFINITY;
  if ( n < 1 || incx < 1) norm = 0;
  else {
    var scale = 0;
    var ssq = 1;
    for (var ix = 1; ix <= 1 + ( n - 1 ) * incx; ix += incx ) {
      if ( x[ ioffx + ix - 1 ].getReal() != 0 ) {
        var temp = Math.abs( x[ ioffx + ix - 1 ].getReal() );
        if ( scale < temp ) {
          ssq = 1 + ssq * Math.pow( scale / temp , 2 );
          scale = temp;
        } else {
          ssq += Math.pow( temp / scale , 2 );
        }
      }
      if ( x[ ioffx + ix - 1 ].getImag() != 0 ) {
        temp = Math.abs( x[ ioffx + ix - 1 ].getImag() );
        if ( scale < temp ) {
          ssq = 1 + ssq * Math.pow( scale / temp , 2 );
          scale = temp;
        } else {
          ssq += Math.pow( temp / scale , 2 );
        }
      }
    }
    norm = scale * Math.sqrt( ssq );
  }
  return norm;
}
//*************************************************************************
//  _asum = sum_i | x_i }
Blas1.dasum = function( n, x, incx , ioffx ) {
  var dtemp = 0;
  if ( n <= 0 || incx <= 0 ) return dtemp;
  var nincx = n * incx;
  for ( var i = 1; i <= nincx; i += incx ) {
    dtemp += Math.abs( x[ ioffx + i - 1 ] );
  }
  return dtemp;
}
Blas1.dzasum = function( n, x, incx , ioffx ) {
  var val = 0;
  if ( n > 0 && incx > 0 ) {
    var ix = 1;
    for ( var i = 1; i <= n; i ++ ) {
//    corresponds to netlib source, but not lapack library:
//      val += ComplexMath.abs( x[ ioffx + ix - 1 ] );
//    corresponds to lapack library:
      val += Math.abs( x[ ioffx + ix - 1 ].getReal() )
           + Math.abs( x[ ioffx + ix - 1 ].getImag() );
      ix += incx;
    }
  }
  return val;
}
//*************************************************************************
Blas1.maxloc = function( ibeg, iend, x, incx, ioffx) {
  var val = 0;
  if ( iend < ibeg || incx <= 0 ) return val;
  val = 1;
  if ( ibeg == iend ) return val;
  var xmax = x[ ioffx + ibeg - 1 ];
  for ( var i = incx; i <= iend - ibeg; i += incx ) {
    var xi = x[ ioffx + ibeg + i - 1 ];
    if ( xi > xmax ) {
      val = i + 1;
      xmax = xi;
    }
  }
  return val;
}
//*************************************************************************
//  i_amax = first i such that | x_i | = || x ||_inf
Blas1.idamax = function( n, x, incx , ioffx ) {
  var val = 0;
  if ( n < 1 || incx <= 0 ) return val;
  val = 1;
  if ( n ==1 ) return val;
  var ix = 1;
  var xmax = Math.abs( x[ ioffx ] );
  ix += incx;
  for ( var i = 2; i <= n; i ++ ) {
    if ( Math.abs( x[ ioffx + ix - 1 ] ) > xmax ) {
      val = i;
      xmax = Math.abs( x[ ioffx + ix - 1 ] );
    }
    ix += incx;
  }
  return val;
}
Blas1.izamax = function( n, x, incx , ioffx ) {
  var val = 0;
  if ( n < 1 || incx <= 0 ) return val;
  val = 1;
  if ( n ==1 ) return val;
  var ix = 1;
  var xmax = ComplexMath.abs( x[ ioffx ] );
  ix += incx;
  for ( var i = 2; i <= n; i ++ ) {
    if ( ComplexMath.abs( x[ ioffx + ix - 1 ] ) > xmax ) {
      val = i;
      xmax = ComplexMath.abs( x[ ioffx + ix - 1 ] );
    }
    ix += incx;
  }
  return val;
}
//*************************************************************************
//  x <--> y
Blas1.dswap = function( n, x, incx, y, incy, ioffx, ioffy ) {
  if ( n <= 0 ) return;
  var ix = 1;
  var iy = 1;
  if ( incx < 0 ) ix = - ( n - 1 ) * incx + 1;
  if ( incy < 0 ) iy = - ( n - 1 ) * incy + 1;
  for ( var i = 1; i <= n; i ++ ) {
//  var xi = x[ ioffx + ix ]; // NOT a reference
    var temp = x[ ioffx + ix - 1 ];
    x[ ioffx + ix - 1 ] = y[ ioffy + iy - 1 ];
    y[ ioffy + iy - 1 ] = temp;
    ix += incx;
    iy += incy;
  }
}
Blas1.zswap = function( n, x, incx, y, incy, ioffx, ioffy ) {
  if ( n <= 0 ) return;
  var ix = 1;
  var iy = 1;
  if ( incx < 0 ) ix = - ( n - 1 ) * incx + 1;
  if ( incy < 0 ) iy = - ( n - 1 ) * incy + 1;
  for ( var i = 1; i <= n; i ++ ) {
    var xi = x[ ioffx + ix - 1 ];//xi is a ref to x[ioff+ix-1]
    var yi = y[ ioffy + iy - 1 ];
    var temp = new Complex( xi.getReal(), xi.getImag() );
    xi = yi;
    yi = temp;
    ix += incx;
    iy += incy;
  }
}
//*************************************************************************
//  y <-- x
Blas1.dcopy = function( n, x, incx, y, incy, ioffx, ioffy ) {
  if ( n <= 0 ) return;
  var ix = 1;
  var iy = 1;
  if ( incx < 0 ) ix = - ( n - 1 ) * incx + 1;
  if ( incy < 0 ) iy = - ( n - 1 ) * incy + 1;
  for ( var i = 1; i <= n; i ++ ) {
    y[ioffy + iy - 1 ] = x[ ioffx + ix - 1 ];
    ix += incx;
    iy += incy;
  }
}
Blas1.zcopy = function( n, x, incx, y, incy, ioffx, ioffy ) {
  if ( n <= 0 ) return;
  var ix = 1;
  var iy = 1;
  if ( incx < 0 ) ix = - ( n - 1 ) * incx + 1;
  if ( incy < 0 ) iy = - ( n - 1 ) * incy + 1;
  for ( var i = 1; i <= n; i ++ ) {
    y[ioffy + iy - 1 ] = x[ ioffx + ix - 1 ];
    ix += incx;
    iy += incy;
  }
}
//*************************************************************************
//  x <-- x * a
Blas1.dscal = function( n, a, x, incx, ioffx ) {
  if ( n <= 0 || incx <= 0 ) return;
  var nincx = n * incx;
  for (var i = 1; i <= nincx; i += incx ) x[ ioffx + i - 1 ] *= a;
}
Blas1.zscal = function( n, a, x, incx, ioffx ) {
  if ( n <= 0 || incx <= 0 ) return;
  var ix = 1;
  for (var i = 1; i <= n; i ++ ) {
    x[ ioffx + ix - 1 ].timesEquals(a);
    ix += incx;
  }
}
//*************************************************************************
//  y += x * a
Blas1.daxpy = function( n, a, x, incx, y, incy, ioffx, ioffy ) {
  if ( n <= 0 ) return;
  if ( a == 0. ) return;
  var ix = 1;
  var iy = 1;
  if ( incx < 0 ) ix = - ( n - 1 ) * incx + 1;
  if ( incy < 0 ) iy = - ( n - 1 ) * incy + 1;
  for (var i = 1; i <= n; i ++ ) {
    y[ ioffy + iy - 1 ] += a * x[ ioffx + ix - 1 ];
    ix += incx;
    iy += incy;
  }
}
Blas1.zaxpy = function( n, a, x, incx, y, incy, ioffx, ioffy ) {
  if ( n <= 0 ) return;
  if ( ComplexMath.abs(a) == 0 ) return;
  var ix = 1;
  var iy = 1;
  if ( incx < 0 ) ix = - ( n - 1 ) * incx + 1;
  if ( incy < 0 ) iy = - ( n - 1 ) * incy + 1;
  for (var i = 1; i <= n; i ++ ) {
    var yiy = y[ ioffy + iy - 1 ];
    var xix = x[ ioffx + ix - 1 ];
    yiy.plusEquals( ComplexMath.times( xix, a ) );
    ix += incx;
    iy += incy;
  }
}
//*************************************************************************
//  _dot = x^T y  or bar{x}^T y
Blas1.ddot = function( n, x, incx, y, incy, ioffx, ioffy ) {
  var dtemp = 0;
  if ( n <= 0 ) return dtemp;
  var ix = 1;
  var iy = 1;
  if ( incx < 0 ) ix = - ( n - 1 ) * incx + 1;
  if ( incy < 0 ) iy = - ( n - 1 ) * incy + 1;
  for (var i = 1; i <= n; i ++ ) {
    dtemp += x[ ioffx + ix - 1 ] * y[ ioffy + iy - 1 ];
    ix += incx;
    iy += incy;
  }
  return dtemp;
}
Blas1.zdotu = function( n, x, incx, y, incy, ioffx, ioffy ) {
  var val = new Complex( 0 , 0 );
  if ( n <= 0 ) return val;
  var ix = 1;
  var iy = 1;
  if ( incx < 0 ) ix = - ( n - 1 ) * incx + 1;
  if ( incy < 0 ) iy = - ( n - 1 ) * incy + 1;
  for (var i = 1; i <= n; i ++ ) {
    val.plusEquals(
      ComplexMath.times( x[ ioffx + ix - 1 ] , y[ ioffy + iy - 1 ] ) );
    ix += incx;
    iy += incy;
  }
  return val;
}
Blas1.zdotc = function( n, x, incx, y, incy, ioffx, ioffy ) {
  var val = new Complex( 0, 0 );
  if ( n <= 0 ) return val;
  var ix = 1;
  var iy = 1;
  if ( incx < 0 ) ix = - ( n - 1 ) * incx + 1;
  if ( incy < 0 ) iy = - ( n - 1 ) * incy + 1;
  for (var i = 1; i <= n; i ++ ) {
    val.plusEquals( ComplexMath.times(
      ComplexMath.conj(x[ ioffx + ix - 1 ]) , y[ ioffy + iy - 1 ] ) );
    ix += incx;
    iy += incy;
  }
  return val;
}
//*************************************************************************
//  [ x^T ] <-- [  c s ] [ x^T ]
//  [ y^T ] <-- [ -s c ] [ y^T ]
Blas1.drot = function( n, x, incx, y, incy, cr, sr, ioffx, ioffy ) {
  if ( n <= 0 ) return;
  var ix = 1;
  var iy = 1;
  if ( incx < 0 ) ix = - ( n - 1 ) * incx + 1;
  if ( incy < 0 ) iy = - ( n - 1 ) * incy + 1;
  for ( var i = 1; i <= n; i ++ ) {
    var temp = cr * x[ ioffx + ix - 1 ] + sr * y[ ioffy + iy - 1 ];
    y[ ioffy + iy - 1 ] = cr*y[ ioffy + iy - 1 ] - sr*x[ ioffx + ix - 1 ];
    x[ ioffx + ix - 1 ] = temp;
    ix += incx;
    iy += incy;
  }
}

function testBlas1() {
//test_drotg();
//test_drotmg();
//test_zrotg();
//test_dnrm2();
//test_dznrm2();
//test_dasum();
//test_dzasum();
//test_idamax();
//test_izamax();
//test_dswap();
//test_zswap();
//test_dcopy();
//test_zcopy();
//test_dscal();
//test_zscal();
//test_daxpy();
//test_zaxpy();
//test_ddot();
//test_zdotu();
//test_zdotc();
  test_drot();
  test_drotm();
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_drotg() {
  document.getElementById("debug_textarea").value +=
    "testing drotg *************\n";
  var a = new NumberReference( 0. );
  var b = new NumberReference( 0. );
  var c = new NumberReference( Number.POSITIVE_INFINITY);
  var s= new NumberReference( Number.POSITIVE_INFINITY);
  Blas1.drotg( a, b, c, s);
  document.getElementById("debug_textarea").value +=
    "a,b,c,s = " + a.getValue() + " , " + b.getValue() + " , "
    + c.getValue() + " , " + s.getValue() + "\n";
  a.setValue( 1. );
  Blas1.drotg( a, b, c, s);
  document.getElementById("debug_textarea").value +=
    "a,b,c,s = " + a.getValue() + " , " + b.getValue() + " , "
    + c.getValue() + " , " + s.getValue() + "\n";

  a.setValue( 0. );
  b.setValue( 1. );
  Blas1.drotg( a, b, c, s);
  document.getElementById("debug_textarea").value +=
    "a,b,c,s = " + a.getValue() + " , " + b.getValue() + " , "
    + c.getValue() + " , " + s.getValue() + "\n";

  a.setValue( 4. );
  b.setValue( 3. );
  Blas1.drotg( a, b, c, s);
  document.getElementById("debug_textarea").value +=
    "a,b,c,s = " + a.getValue() + " , " + b.getValue() + " , "
    + c.getValue() + " , " + s.getValue() + "\n";

  a.setValue( -4. );
  b.setValue( 3. );
  Blas1.drotg( a, b, c, s);
  document.getElementById("debug_textarea").value +=
    "a,b,c,s = " + a.getValue() + " , " + b.getValue() + " , "
    + c.getValue() + " , " + s.getValue() + "\n";

  a.setValue( 3. );
  b.setValue( 4. );
  Blas1.drotg( a, b, c, s);
  document.getElementById("debug_textarea").value +=
    "a,b,c,s = " + a.getValue() + " , " + b.getValue() + " , "
    + c.getValue() + " , " + s.getValue() + "\n";

  a.setValue( 3. );
  b.setValue( -4. );
  Blas1.drotg( a, b, c, s);
  document.getElementById("debug_textarea").value +=
    "a,b,c,s = " + a.getValue() + " , " + b.getValue() + " , "
    + c.getValue() + " , " + s.getValue() + "\n";

  a.setValue( 4. );
  b.setValue( -3. );
  Blas1.drotg( a, b, c, s);
  document.getElementById("debug_textarea").value +=
    "a,b,c,s = " + a.getValue() + " , " + b.getValue() + " , "
    + c.getValue() + " , " + s.getValue() + "\n";

  a.setValue( -4. );
  b.setValue( -3. );
  Blas1.drotg( a, b, c, s);
  document.getElementById("debug_textarea").value +=
    "a,b,c,s = " + a.getValue() + " , " + b.getValue() + " , "
    + c.getValue() + " , " + s.getValue() + "\n";

  a.setValue( -3. );
  b.setValue( 4. );
  Blas1.drotg( a, b, c, s);
  document.getElementById("debug_textarea").value +=
    "a,b,c,s = " + a.getValue() + " , " + b.getValue() + " , "
    + c.getValue() + " , " + s.getValue() + "\n";

  a.setValue( -3. );
  b.setValue( -4. );
  Blas1.drotg( a, b, c, s);
  document.getElementById("debug_textarea").value +=
    "a,b,c,s = " + a.getValue() + " , " + b.getValue() + " , "
    + c.getValue() + " , " + s.getValue() + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zrotg() {
  document.getElementById("debug_textarea").value +=
    "testing zrotg *************" + "\n";
  var a = new Complex( 0., 0. );
  var b = new Complex( 0., 0. );
  var c = new NumberReference( Number.POSITIVE_INFINITY );
  var s = new Complex( );
  Blas1.zrotg( a, b, c, s );
  document.getElementById("debug_textarea").value +=
    "a,b,c,s = " + a.toString() + " , " + b.toString() + " , "
        + c.getValue() + " , " + s.toString() + "\n";

  a.setValue( new Complex( 1., 0. ) );
  Blas1.zrotg( a, b, c, s );
  document.getElementById("debug_textarea").value +=
    "a,b,c,s = " + a.toString() + " , " + b.toString() + " , "
        + c.getValue() + " , " + s.toString() + "\n";

  a.setValue( new Complex( 0., 1. ) );
  Blas1.zrotg( a, b, c, s );
  document.getElementById("debug_textarea").value +=
    "a,b,c,s = " + a.toString() + " , " + b.toString() + " , "
        + c.getValue() + " , " + s.toString() + "\n";

  a.setValue( new Complex( 0., 0. ) );
  b.setValue( new Complex( 1., 0. ) );
  Blas1.zrotg( a, b, c, s );
  document.getElementById("debug_textarea").value +=
    "a,b,c,s = " + a.toString() + " , " + b.toString() + " , "
        + c.getValue() + " , " + s.toString() + "\n";

  b.setValue( new Complex( 0., 1. ) );
  Blas1.zrotg( a, b, c, s );
  document.getElementById("debug_textarea").value +=
    "a,b,c,s = " + a.toString() + " , " + b.toString() + " , "
        + c.getValue() + " , " + s.toString() + "\n";

  a.setValue( new Complex( 4., 0. ) );
  b.setValue( new Complex( 3., 0. ) );
  Blas1.zrotg( a, b, c, s );
  document.getElementById("debug_textarea").value +=
    "a,b,c,s = " + a.toString() + " , " + b.toString() + " , "
        + c.getValue() + " , " + s.toString() + "\n";

  a.setValue( new Complex( 0., 4. ) );
  Blas1.zrotg( a, b, c, s );
  document.getElementById("debug_textarea").value +=
    "a,b,c,s = " + a.toString() + " , " + b.toString() + " , "
        + c.getValue() + " , " + s.toString() + "\n";

  a.setValue( new Complex( 4., 0. ) );
  b.setValue( new Complex( 0., 3. ) );
  Blas1.zrotg( a, b, c, s );
  document.getElementById("debug_textarea").value +=
    "a,b,c,s = " + a.toString() + " , " + b.toString() + " , "
        + c.getValue() + " , " + s.toString() + "\n";

  a.setValue( new Complex( 0., 4. ) );
  Blas1.zrotg( a, b, c, s );
  document.getElementById("debug_textarea").value +=
    "a,b,c,s = " + a.toString() + " , " + b.toString() + " , "
        + c.getValue() + " , " + s.toString() + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_drotmg() {
  document.getElementById("debug_textarea").value +=
    "testing drotmg *************" + "\n";
  var d1 = new NumberReference( -1. );
  var d2 = new NumberReference( 0. );
  var x1 = new NumberReference( 3. );
  var y1 = new NumberReference( 4. );
  var ioffparam = 1;
  var param = new Array( ioffparam + 5 );
  Blas1.drotmg( d1, d2, x1, y1, param, ioffparam );
  document.getElementById("debug_textarea").value +=
    "d1,d2,x1,y1 = " + d1.getValue() + " , " + d2.getValue() + " , "
    + x1.getValue() + " , " + y1.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "param = "
    + param[ ioffparam + 0 ] + " "
    + param[ ioffparam + 1 ] + " "
    + param[ ioffparam + 2 ] + " "
    + param[ ioffparam + 3 ] + " "
    + param[ ioffparam + 4 ]  + "\n";
  d1.setValue( 0. );
  d2.setValue( 0. );
  x1.setValue( 3. );
  y1.setValue( 4. );
  Blas1.drotmg( d1, d2, x1, y1, param, ioffparam );
  document.getElementById("debug_textarea").value +=
    "d1,d2,x1,y1 = " + d1.getValue() + " , " + d2.getValue() + " , "
    + x1.getValue() + " , " + y1.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "param = "
    + param[ ioffparam + 0 ] + " "
    + param[ ioffparam + 1 ] + " "
    + param[ ioffparam + 2 ] + " "
    + param[ ioffparam + 3 ] + " "
    + param[ ioffparam + 4 ]  + "\n";

  d1.setValue( 1. );
  d2.setValue( 1. );
  x1.setValue( 4. );
  y1.setValue( 3. );
  Blas1.drotmg( d1, d2, x1, y1, param, ioffparam );
  document.getElementById("debug_textarea").value +=
    "d1,d2,x1,y1 = " + d1.getValue() + " , " + d2.getValue() + " , "
    + x1.getValue() + " , " + y1.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "param = "
    + param[ ioffparam + 0 ] + " "
    + param[ ioffparam + 1 ] + " "
    + param[ ioffparam + 2 ] + " "
    + param[ ioffparam + 3 ] + " "
    + param[ ioffparam + 4 ]  + "\n";

  d1.setValue( 1. );
  d2.setValue( -1. );
  x1.setValue( 3. );
  y1.setValue( 4. );
  Blas1.drotmg( d1, d2, x1, y1, param, ioffparam );
  document.getElementById("debug_textarea").value +=
    "d1,d2,x1,y1 = " + d1.getValue() + " , " + d2.getValue() + " , "
    + x1.getValue() + " , " + y1.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "param = "
    + param[ ioffparam + 0 ] + " "
    + param[ ioffparam + 1 ] + " "
    + param[ ioffparam + 2 ] + " "
    + param[ ioffparam + 3 ] + " "
    + param[ ioffparam + 4 ]  + "\n";

  d1.setValue( 1. );
  d2.setValue( 1. );
  x1.setValue( 3. );
  y1.setValue( 4. );
  Blas1.drotmg( d1, d2, x1, y1, param, ioffparam );
  document.getElementById("debug_textarea").value +=
    "d1,d2,x1,y1 = " + d1.getValue() + " , " + d2.getValue() + " , "
    + x1.getValue() + " , " + y1.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "param = "
    + param[ ioffparam + 0 ] + " "
    + param[ ioffparam + 1 ] + " "
    + param[ ioffparam + 2 ] + " "
    + param[ ioffparam + 3 ] + " "
    + param[ ioffparam + 4 ]  + "\n";

  d1.setValue( 16777216. );
  d1.setValue( d1.getValue() * d1.getValue() );
  d2.setValue( 1. );
  x1.setValue( 3. );
  y1.setValue( 4. );
  Blas1.drotmg( d1, d2, x1, y1, param, ioffparam );
  document.getElementById("debug_textarea").value +=
    "d1,d2,x1,y1 = " + d1.getValue() + " , " + d2.getValue() + " , "
    + x1.getValue() + " , " + y1.getValue() + "\n";
  document.getElementById("debug_textarea").value +=
    "param = "
    + param[ ioffparam + 0 ] + " "
    + param[ ioffparam + 1 ] + " "
    + param[ ioffparam + 2 ] + " "
    + param[ ioffparam + 3 ] + " "
    + param[ ioffparam + 4 ]  + "\n";

  d1.setValue( 1. );
  d2.setValue( 16777216. );
  d2.setValue( d2.getValue() * d2.getValue() );
  x1.setValue( 3. );
  y1.setValue( 4. );
  Blas1.drotmg( d1, d2, x1, y1, param, ioffparam );
  document.getElementById("debug_textarea").value +=
        "d1,d2,x1,y1 = " + d1.getValue() + " , " + d2.getValue() + " , "
        + x1.getValue() + " , " + y1.getValue() + "\n";
      document.getElementById("debug_textarea").value +=
        "param = "
        + param[ ioffparam + 0 ] + " "
        + param[ ioffparam + 1 ] + " "
        + param[ ioffparam + 2 ] + " "
        + param[ ioffparam + 3 ] + " "
        + param[ ioffparam + 4 ]  + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dnrm2() {
  document.getElementById("debug_textarea").value +=
    "testing dnrm2 *************" + "\n";
  var ioffx = 1;
  var x = new Array( ioffx + 6 );
  for ( var i = 0; i < 6; i ++ ) { x[ ioffx + i ] = Number( i ); }
  var xnrm = Blas1.dnrm2( 0 , x , 1 , ioffx );
  document.getElementById("debug_textarea").value +=
    "dnrm2(0,x,1) = " + xnrm + "\n";
  xnrm = Blas1.dnrm2( 1 , x , 1 , ioffx );
  document.getElementById("debug_textarea").value +=
    "dnrm2(1,x,1) = " + xnrm + "\n";
  xnrm = Blas1.dnrm2( 6 , x , 1 , ioffx );
  document.getElementById("debug_textarea").value +=
    "dnrm2(6,x,1) = " + xnrm + "\n";
  xnrm = Blas1.dnrm2( 2 , x , 3 , ioffx );
  document.getElementById("debug_textarea").value +=
    "dnrm2(2,x,3) = " + xnrm + "\n";
  xnrm = Blas1.dnrm2( 3 , x , 2 , ioffx );
  document.getElementById("debug_textarea").value +=
    "dnrm2(3,x,2) = " + xnrm + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dznrm2() {
  document.getElementById("debug_textarea").value +=
    "testing dznrm2 *************" + "\n";
  var ioffx = 1;
  var x = new Array( ioffx + 6 );
  for ( var i = 0; i < 6; i ++ ) {
    x[ ioffx + i ] = new Complex( Number(i) , Number(i-1) );
  }
  var xnrm = Blas1.dznrm2( 0 , x , 1 , ioffx );
  document.getElementById("debug_textarea").value +=
    "dznrm2(0,x,1) = " + xnrm + "\n";
  xnrm = Blas1.dznrm2( 1 , x , 1 , ioffx );
  document.getElementById("debug_textarea").value +=
    "dznrm2(1,x,1) = " + xnrm + "\n";
  xnrm = Blas1.dznrm2( 6 , x , 1 , ioffx );
  document.getElementById("debug_textarea").value +=
    "dznrm2(6,x,1) = " + xnrm + "\n";
  xnrm = Blas1.dznrm2( 2 , x , 3 , ioffx );
  document.getElementById("debug_textarea").value +=
    "dznrm2(2,x,3) = " + xnrm + "\n";
  xnrm = Blas1.dznrm2( 3 , x , 2 , ioffx );
  document.getElementById("debug_textarea").value +=
    "dznrm2(3,x,2) = " + xnrm + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dasum() {
  document.getElementById("debug_textarea").value +=
    "testing dasum *************" + "\n";
  var ioffx = 1;
  var x = new Array( ioffx + 6 );
  for ( var i = 0; i < 6; i ++ ) { x[ ioffx + i ] = Number( i ); }
  var xnrm = Blas1.dasum( 0 , x , 1 , ioffx );
  document.getElementById("debug_textarea").value +=
    "dasum(0,x,1) = " + xnrm + "\n";
  xnrm = Blas1.dasum( 1 , x , 0 , ioffx );
  document.getElementById("debug_textarea").value +=
    "dasum(1,x,0) = " + xnrm + "\n";
  xnrm = Blas1.dasum( 6 , x , 1 , ioffx );
  document.getElementById("debug_textarea").value +=
    "dasum(6,x,1) = " + xnrm + "\n";
  xnrm = Blas1.dasum( 2 , x , 3 , ioffx );
  document.getElementById("debug_textarea").value +=
    "dasum(2,x,3) = " + xnrm + "\n";
  xnrm = Blas1.dasum( 3 , x , 2 , ioffx );
  document.getElementById("debug_textarea").value +=
    "dasum(3,x,2) = " + xnrm + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dzasum() {
  document.getElementById("debug_textarea").value +=
    "testing dzasum *************" + "\n";
  var ioffx = 1;
  var x = new Array( ioffx + 6 );
  for ( var i = 0; i < 6; i ++ ) {
    x[ ioffx + i ] = new Complex( Number(i) , Number(i-3) );
  }
  var xnrm = Blas1.dzasum( 0 , x , 1 , ioffx );
  document.getElementById("debug_textarea").value +=
    "dzasum(0,x,1) = " + xnrm + "\n";
  xnrm = Blas1.dzasum( 1 , x , 0 , ioffx );
  document.getElementById("debug_textarea").value +=
    "dzasum(1,x,0) = " + xnrm + "\n";
  xnrm = Blas1.dzasum( 6 , x , 1 , ioffx );
  document.getElementById("debug_textarea").value +=
    "dzasum(6,x,1) = " + xnrm + "\n";
  xnrm = Blas1.dzasum( 2 , x , 3 , ioffx );
  document.getElementById("debug_textarea").value +=
    "dzasum(2,x,3) = " + xnrm + "\n";
  xnrm = Blas1.dzasum( 3 , x , 2 , ioffx );
  document.getElementById("debug_textarea").value +=
    "dzasum(3,x,2) = " + xnrm + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_idamax() {
  document.getElementById("debug_textarea").value +=
    "testing idamax *************" + "\n";
  var ioffx = 1;
  var x = new Array( ioffx + 6 );
  x[ ioffx + 0 ] = 1.;
  x[ ioffx + 1 ] = 2.;
  x[ ioffx + 2 ] = 3.;
  x[ ioffx + 3 ] = 5.;
  x[ ioffx + 4 ] = 4.;
  x[ ioffx + 5 ] = 0.;
  var m = Blas1.idamax( 0 , x , 1 , ioffx );
  document.getElementById("debug_textarea").value +=
    "idamax(0,x,1) = " + m + "\n";
  m = Blas1.idamax( 1 , x , 0 , ioffx );
  document.getElementById("debug_textarea").value +=
    "idamax(1,x,0) = " + m + "\n";
  m = Blas1.idamax( 6 , x , 1 , ioffx );
  document.getElementById("debug_textarea").value +=
    "idamax(6,x,1) = " + m + "\n";
  m = Blas1.idamax( 2 , x , 3 , ioffx );
  document.getElementById("debug_textarea").value +=
    "idamax(2,x,3) = " + m + "\n";
  m = Blas1.idamax( 3 , x , 2 , ioffx );
  document.getElementById("debug_textarea").value +=
    "idamax(3,x,2) = " + m + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_izamax() {
  document.getElementById("debug_textarea").value +=
    "testing izamax *************" + "\n";
  var ioffx = 1;
  var x = new Array( ioffx + 6 );
  x[ ioffx + 0 ] = new Complex( 1. , -2. );
  x[ ioffx + 1 ] = new Complex( 2. , -1. );
  x[ ioffx + 2 ] = new Complex( 3. , 0. );
  x[ ioffx + 3 ] = new Complex( 5. , 2. );
  x[ ioffx + 4 ] = new Complex( 4. , 1. );
  x[ ioffx + 5 ] = new Complex( 0. , -3. );
  var m = Blas1.izamax( 0 , x , 1 , ioffx );
  document.getElementById("debug_textarea").value +=
    "izamax(0,x,1) = " + m + "\n";
  m = Blas1.izamax( 1 , x , 0 , ioffx );
  document.getElementById("debug_textarea").value +=
    "izamax(1,x,0) = " + m + "\n";
  m = Blas1.izamax( 6 , x , 1 , ioffx );
  document.getElementById("debug_textarea").value +=
    "izamax(6,x,1) = " + m + "\n";
  m = Blas1.izamax( 2 , x , 3 , ioffx );
  document.getElementById("debug_textarea").value +=
    "izamax(2,x,3) = " + m + "\n";
  m = Blas1.izamax( 3 , x , 2 , ioffx );
  document.getElementById("debug_textarea").value +=
    "izamax(3,x,2) = " + m + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dswap() {
  document.getElementById("debug_textarea").value +=
    "testing dswap *************" + "\n";
  var ioffx = 1;
  var ioffy = 2;
  var ioffz = 3;
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 3 );
  var z = new Array( ioffz + 2 );
  var i = -1;
  for ( i = 0; i < 6; i ++ ) x[ ioffx + i ] = i;
  for ( i = 0; i < 3; i ++ ) y[ ioffy + i ] = i+6;
  Blas1.dswap( 0 , x , 1 , y , 1 , ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "x = " + x[ioffx + 0] + " " + x[ioffx + 1] + " "
        + x[ioffx + 2] + " " + x[ioffx + 3] + " "
        + x[ioffx + 4] + " " + x[ioffx + 5] + "\n";
  document.getElementById("debug_textarea").value +=
    "y = " + y[ioffy + 0] + " " + y[ioffy + 1] + " "
        + y[ioffy + 2] + "\n";

  for ( i = 0; i < 6; i ++ ) x[ ioffx + i ] = i;
  for ( i = 0; i < 3; i ++ ) y[ ioffy + i ] = i+6;
  Blas1.dswap( 3 , x , 1 , y , 1 , ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "x = " + x[ioffx + 0] + " " + x[ioffx + 1] + " "
        + x[ioffx + 2] + " " + x[ioffx + 3] + " "
        + x[ioffx + 4] + " " + x[ioffx + 5] + "\n";
  document.getElementById("debug_textarea").value +=
    "y = " + y[ioffy + 0] + " " + y[ioffy + 1] + " "
        + y[ioffy + 2] + "\n";

  for ( i = 0; i < 6; i ++ ) x[ ioffx + i ] = i;
  for ( i = 0; i < 3; i ++ ) y[ ioffy + i ] = i+6;
  Blas1.dswap( 2 , x , 2 , y , 2 , ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "x = " + x[ioffx + 0] + " " + x[ioffx + 1] + " "
        + x[ioffx + 2] + " " + x[ioffx + 3] + " "
        + x[ioffx + 4] + " " + x[ioffx + 5] + "\n";
  document.getElementById("debug_textarea").value +=
    "y = " + y[ioffy + 0] + " " + y[ioffy + 1] + " "
        + y[ioffy + 2] + "\n";

  for ( i = 0; i < 6; i ++ ) x[ ioffx + i ] = i;
  for ( i = 0; i < 2; i ++ ) z[ ioffz + i ] = -1-i;
  Blas1.dswap( 2 , x , 3 , z , 1 , ioffx, ioffz );
  document.getElementById("debug_textarea").value +=
    "x = " + x[ioffx + 0] + " " + x[ioffx + 1] + " "
        + x[ioffx + 2] + " " + x[ioffx + 3] + " "
        + x[ioffx + 4] + " " + x[ioffx + 5] + "\n";
  document.getElementById("debug_textarea").value +=
    "z = " + z[ioffz + 0] + " " + z[ioffz + 1] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zswap() {
  document.getElementById("debug_textarea").value +=
    "testing zswap *************" + "\n";
  var ioffx = 1;
  var ioffy = 2;
  var ioffz = 3;
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 3 );
  var z = new Array( ioffz + 2 );
  var i = -1;
  for ( i = 0; i < 6; i ++ ) x[ ioffx + i ] = new Complex( i , i - 1 );
  for ( i = 0; i < 3; i ++ ) y[ ioffy + i ] = new Complex( i+6 , i-6 );
  Blas1.zswap( 0 , x , 1 , y , 1 , ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "x = " + x[ioffx + 0] + " " + x[ioffx + 1] + " "
        + x[ioffx + 2] + " " + x[ioffx + 3] + " "
        + x[ioffx + 4] + " " + x[ioffx + 5] + "\n";
  document.getElementById("debug_textarea").value +=
    "y = " + y[ioffy + 0] + " " + y[ioffy + 1] + " "
        + y[ioffy + 2] + "\n";

  for ( i = 0; i < 6; i ++ ) x[ ioffx + i ] = new Complex( i , i - 1 );
  for ( i = 0; i < 3; i ++ ) y[ ioffy + i ] = new Complex( i+6 , i-6 );
  Blas1.zswap( 3 , x , 1 , y , 1 , ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "x = " + x[ioffx + 0] + " " + x[ioffx + 1] + " "
        + x[ioffx + 2] + " " + x[ioffx + 3] + " "
        + x[ioffx + 4] + " " + x[ioffx + 5] + "\n";
  document.getElementById("debug_textarea").value +=
    "y = " + y[ioffy + 0] + " " + y[ioffy + 1] + " "
        + y[ioffy + 2] + "\n";

  for ( i = 0; i < 6; i ++ ) x[ ioffx + i ] = new Complex( i , i - 1 );
  for ( i = 0; i < 3; i ++ ) y[ ioffy + i ] = new Complex( i+6 , i-6 );
  Blas1.zswap( 2 , x , 2 , y , 2 , ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "x = " + x[ioffx + 0] + " " + x[ioffx + 1] + " "
        + x[ioffx + 2] + " " + x[ioffx + 3] + " "
        + x[ioffx + 4] + " " + x[ioffx + 5] + "\n";
  document.getElementById("debug_textarea").value +=
    "y = " + y[ioffy + 0] + " " + y[ioffy + 1] + " "
        + y[ioffy + 2] + "\n";

  for ( i = 0; i < 6; i ++ ) x[ ioffx + i ] = new Complex( i , i - 1 );
  for ( i = 0; i < 2; i ++ ) {
    z[ ioffz + i ] = new Complex( -1-i , 1+i );
  }
  Blas1.zswap( 2 , x , 3 , z , 1 , ioffx, ioffz );
  document.getElementById("debug_textarea").value +=
    "x = " + x[ioffx + 0] + " " + x[ioffx + 1] + " "
        + x[ioffx + 2] + " " + x[ioffx + 3] + " "
        + x[ioffx + 4] + " " + x[ioffx + 5] + "\n";
  document.getElementById("debug_textarea").value +=
    "z = " + z[ioffz + 0] + " " + z[ioffz + 1] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dcopy() {
  document.getElementById("debug_textarea").value +=
    "testing dcopy *************" + "\n";
  var ioffx = 1;
  var ioffy = 2;
  var ioffz = 3;
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 3 );
  var z = new Array( ioffz + 2 );
  var i = -1;
  for ( i = 0; i < 6; i ++ ) x[ ioffx + i ] = i;
  for ( i = 0; i < 3; i ++ ) y[ ioffy + i ] = Number.POSITIVE_INFINITY;
  Blas1.dcopy( 0 , x , 1 , y , 1 , ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "x = " + x[ioffx + 0] + " " + x[ioffx + 1] + " "
        + x[ioffx + 2] + " " + x[ioffx + 3] + " "
        + x[ioffx + 4] + " " + x[ioffx + 5] + "\n";
  document.getElementById("debug_textarea").value +=
    "y = " + y[ioffy + 0] + " " + y[ioffy + 1] + " "
        + y[ioffy + 2] + "\n";

  for ( i = 0; i < 6; i ++ ) x[ ioffx + i ] = i;
  for ( i = 0; i < 3; i ++ ) y[ ioffy + i ] = Number.POSITIVE_INFINITY;
  Blas1.dcopy( 3 , x , 1 , y , 1, ioffx, ioffy  );
  document.getElementById("debug_textarea").value +=
    "x = " + x[ioffx + 0] + " " + x[ioffx + 1] + " "
        + x[ioffx + 2] + " " + x[ioffx + 3] + " "
        + x[ioffx + 4] + " " + x[ioffx + 5] + "\n";
  document.getElementById("debug_textarea").value +=
    "y = " + y[ioffy + 0] + " " + y[ioffy + 1] + " "
        + y[ioffy + 2] + "\n";

  for ( i = 0; i < 6; i ++ ) x[ ioffx + i ] = i;
  for ( i = 0; i < 3; i ++ ) y[ ioffy + i ] = Number.POSITIVE_INFINITY;
  Blas1.dcopy( 2 , x , 2 , y , 2, ioffx, ioffy  );
  document.getElementById("debug_textarea").value +=
    "x = " + x[ioffx + 0] + " " + x[ioffx + 1] + " "
        + x[ioffx + 2] + " " + x[ioffx + 3] + " "
        + x[ioffx + 4] + " " + x[ioffx + 5] + "\n";
  document.getElementById("debug_textarea").value +=
    "y = " + y[ioffy + 0] + " " + y[ioffy + 1] + " "
        + y[ioffy + 2] + "\n";

  for ( i = 0; i < 6; i ++ ) x[ ioffx + i ] = i;
  for ( i = 0; i < 2; i ++ ) z[ ioffz + i ] = Number.POSITIVE_INFINITY;
  Blas1.dcopy( 2 , x , 3 , z , 1, ioffx, ioffz  );
  document.getElementById("debug_textarea").value +=
    "x = " + x[ioffx + 0] + " " + x[ioffx + 1] + " "
        + x[ioffx + 2] + " " + x[ioffx + 3] + " "
        + x[ioffx + 4] + " " + x[ioffx + 5] + "\n";
  document.getElementById("debug_textarea").value +=
    "z = " + z[ioffz + 0] + " " + z[ioffz + 1] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zcopy() {
  document.getElementById("debug_textarea").value +=
    "testing zcopy *************" + "\n";
  var ioffx = 1;
  var ioffy = 2;
  var ioffz = 3;
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 3 );
  var z = new Array( ioffz + 2 );
  var i = -1;
  for ( i = 0; i < 6; i ++ ) x[ ioffx + i ] = new Complex( i , i-1 );
  for ( i = 0; i < 3; i ++ ) y[ ioffy + i ] = new Complex();
  Blas1.zcopy( 0 , x , 1 , y , 1, ioffx, ioffy  );
  document.getElementById("debug_textarea").value +=
    "x = " + x[ioffx + 0] + " " + x[ioffx + 1] + " "
        + x[ioffx + 2] + " " + x[ioffx + 3] + " "
        + x[ioffx + 4] + " " + x[ioffx + 5] + "\n";
  document.getElementById("debug_textarea").value +=
    "y = " + y[ioffy + 0] + " " + y[ioffy + 1] + " "
        + y[ioffy + 2] + "\n";

  for ( i = 0; i < 6; i ++ ) x[ ioffx + i ] = new Complex( i , i-1 );
  for ( i = 0; i < 3; i ++ ) y[ ioffy + i ] = new Complex();
  Blas1.zcopy( 3 , x , 1 , y , 1, ioffx, ioffy  );
  document.getElementById("debug_textarea").value +=
    "x = " + x[ioffx + 0] + " " + x[ioffx + 1] + " "
        + x[ioffx + 2] + " " + x[ioffx + 3] + " "
        + x[ioffx + 4] + " " + x[ioffx + 5] + "\n";
  document.getElementById("debug_textarea").value +=
    "y = " + y[ioffy + 0] + " " + y[ioffy + 1] + " "
        + y[ioffy + 2] + "\n";

  for ( i = 0; i < 6; i ++ ) x[ ioffx + i ] = new Complex( i , i-1 );
  for ( i = 0; i < 3; i ++ ) y[ ioffy + i ] = new Complex();
  Blas1.zcopy( 2 , x , 2 , y , 2, ioffx, ioffy  );
  document.getElementById("debug_textarea").value +=
    "x = " + x[ioffx + 0] + " " + x[ioffx + 1] + " "
        + x[ioffx + 2] + " " + x[ioffx + 3] + " "
        + x[ioffx + 4] + " " + x[ioffx + 5] + "\n";
  document.getElementById("debug_textarea").value +=
    "y = " + y[ioffy + 0] + " " + y[ioffy + 1] + " "
        + y[ioffy + 2] + "\n";

  for ( i = 0; i < 6; i ++ ) x[ ioffx + i ] = new Complex( i , i-1 );
  for ( i = 0; i < 2; i ++ ) z[ ioffz + i ] = new Complex();
  Blas1.zcopy( 2 , x , 3 , z , 1, ioffx, ioffz  );
  document.getElementById("debug_textarea").value +=
    "x = " + x[ioffx + 0] + " " + x[ioffx + 1] + " "
        + x[ioffx + 2] + " " + x[ioffx + 3] + " "
        + x[ioffx + 4] + " " + x[ioffx + 5] + "\n";
  document.getElementById("debug_textarea").value +=
    "z = " + z[ioffz + 0] + " " + z[ioffz + 1] + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_dscal() {
  document.getElementById("debug_textarea").value +=
    "testing dscal *************" + "\n";
  var ioffx = 1;
  var x = new Array( ioffx + 6 );
  var i = -1;
  var a = 2;
  for ( i = 0; i < 6; i ++ ) { x[ ioffx + i ] = Number( i ); }
  Blas1.dscal( 0, a, x, 1, ioffx );
  document.getElementById("debug_textarea").value +=
    "dscal(0,a,x,1), x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + " "
        + x[ioffx + 3] + " " + x[ioffx + 4] + " " + x[ioffx + 5]
        + "\n";

  for ( i = 0; i < 6; i ++ ) { x[ ioffx + i ] = Number( i ); }
  Blas1.dscal( 6, a, x, 0, ioffx );
  document.getElementById("debug_textarea").value +=
    "dscal(6,a,x,0), x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + " "
        + x[ioffx + 3] + " " + x[ioffx + 4] + " " + x[ioffx + 5]
        + "\n";

  for ( i = 0; i < 6; i ++ ) { x[ ioffx + i ] = Number( i ); }
  Blas1.dscal( 6, a, x, 1, ioffx );
  document.getElementById("debug_textarea").value +=
    "dscal(6,a,x,1), x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + " "
        + x[ioffx + 3] + " " + x[ioffx + 4] + " " + x[ioffx + 5]
        + "\n";

  for ( i = 0; i < 6; i ++ ) { x[ ioffx + i ] = Number( i ); }
  Blas1.dscal( 3, a, x, 2, ioffx );
  document.getElementById("debug_textarea").value +=
    "dscal(3,a,x,2), x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + " "
        + x[ioffx + 3] + " " + x[ioffx + 4] + " " + x[ioffx + 5]
        + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zscal() {
  document.getElementById("debug_textarea").value +=
    "testing zscal *************" + "\n";
  var ioffx = 1;
  var x = new Array( ioffx + 6 );
  var i = -1;
  var a = new Complex( 2., 3. );
  for ( i = 0; i < 6; i ++ ) {
    x[ ioffx + i ] = new Complex( Number( i ), Number( i - 1 ) );
  }
  Blas1.zscal( 0, a, x, 1, ioffx );
  document.getElementById("debug_textarea").value +=
    "zscal(0,a,x,1), x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + " "
        + x[ioffx + 3] + " " + x[ioffx + 4] + " " + x[ioffx + 5]
        + "\n";

  for ( i = 0; i < 6; i ++ ) {
    x[ ioffx + i ] = new Complex( Number( i ), Number( i - 1 ) );
  }
  Blas1.zscal( 6, a, x, 0, ioffx );
  document.getElementById("debug_textarea").value +=
    "zscal(6,a,x,0), x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + " "
        + x[ioffx + 3] + " " + x[ioffx + 4] + " " + x[ioffx + 5]
        + "\n";

  for ( i = 0; i < 6; i ++ ) {
    x[ ioffx + i ] = new Complex( Number( i ), Number( i - 1 ) );
  }
  Blas1.zscal( 6, a, x, 1, ioffx );
  document.getElementById("debug_textarea").value +=
    "zscal(6,a,x,1), x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + " "
        + x[ioffx + 3] + " " + x[ioffx + 4] + " " + x[ioffx + 5]
        + "\n";

  for ( i = 0; i < 6; i ++ ) {
    x[ ioffx + i ] = new Complex( Number( i ), Number( i - 1 ) );
  }
  Blas1.zscal( 3, a, x, 2, ioffx );
  document.getElementById("debug_textarea").value +=
    "zscal(3,a,x,2), x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2] + " "
        + x[ioffx + 3] + " " + x[ioffx + 4] + " " + x[ioffx + 5]
        + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_daxpy() {
  document.getElementById("debug_textarea").value +=
    "testing daxpy *************" + "\n";
  var ioffx = 1;
  var ioffy = 2;
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 6 );
  var i = -1;
  var a = 2;
  for ( i = 0; i < 6; i ++ ) { x[ ioffx + i ] = Number( i ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( 1 - i ); }
  Blas1.daxpy( 0, a, x, 1, y, 1, ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "daxpy(0,a,x,1,y,1), y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2]
        + " " + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5]
        + "\n";

  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( 1 - i ); }
  Blas1.daxpy( 6, a, x, 1, y, 1, ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "daxpy(6,a,x,1,y,1), y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2]
        + " " + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5]
        + "\n";

  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( 1 - i ); }
  Blas1.daxpy( 2, a, x, 3, y, 2, ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "daxpy(6,a,x,3,y,2), y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2]
        + " " + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5]
        + "\n";

  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( 1 - i ); }
  Blas1.daxpy( 2, a, x, -3, y, -2, ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "daxpy(6,a,x,-3,y,-2), y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2]
        + " " + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5]
        + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zaxpy() {
  document.getElementById("debug_textarea").value +=
    "testing zaxpy *************" + "\n";
  var ioffx = 1;
  var ioffy = 2;
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 6 );
  var i = -1;
  var a = new Complex( 2., 3. );
  for ( i = 0; i < 6; i ++ ) {
    x[ ioffx + i ] = new Complex( i , i-1 );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( 1 - i , - i );
  }
  Blas1.zaxpy( 0, a, x, 1, y, 1, ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "zaxpy(0,a,x,1,y,1), y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2]
        + " " + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5]
        + "\n";

  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( 1 - i , - i );
  }
  Blas1.zaxpy( 6, a, x, 1, y, 1, ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "zaxpy(6,a,x,1,y,1), y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2]
        + " " + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5]
        + "\n";

  for ( i = 0; i < 6; i ++ ) { 
    y[ ioffy + i ] = new Complex( 1 - i , - i );
  }
  Blas1.zaxpy( 2, a, x, 3, y, 2, ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "zaxpy(6,a,x,3,y,2), y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2]
        + " " + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5]
        + "\n";

  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( 1 - i , - i );
  }
  Blas1.zaxpy( 2, a, x, -3, y, -2, ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "zaxpy(6,a,x,-3,y,-2), y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2]
        + " " + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5]
        + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_ddot() {
  document.getElementById("debug_textarea").value +=
    "testing ddot *************" + "\n";
  var ioffx = 1;
  var ioffy = 2;
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 6 );
  var i = -1;
  for ( i = 0; i < 6; i ++ ) { x[ ioffx + i ] = Number( i ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( 1 - i ); }
  var d = Blas1.ddot( 0, x, 1, y, 1, ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "ddot(0,x,1,y,1) = " + d + "\n";
  d = Blas1.ddot( 6, x, 1, y, 1, ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "ddot(6,x,1,y,1) = " + d + "\n";
  d = Blas1.ddot( 2, x, 3, y, 2, ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "ddot(2,x,3,y,2) = " + d + "\n";
  d = Blas1.ddot( 2, x, -3, y, -2, ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "ddot(2,x,-3,y,-2) = " + d + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zdotu() {
  document.getElementById("debug_textarea").value +=
    "testing zdotu *************" + "\n";
  var ioffx = 1;
  var ioffy = 2;
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 6 );
  var i = -1;
  for ( i = 0; i < 6; i ++ ) {
    x[ ioffx + i ] = new Complex( i , i-1 );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( 1 - i , -i );
  }
  var d = new Complex();
  d.setValue( Blas1.zdotu( 0, x, 1, y, 1, ioffx, ioffy ) );
  document.getElementById("debug_textarea").value +=
    "zdotu(0,x,1,y,1) = " + d.toString() + "\n";
  d.setValue( Blas1.zdotu( 6, x, 1, y, 1, ioffx, ioffy ) );
  document.getElementById("debug_textarea").value +=
    "zdotu(6,x,1,y,1) = " + d.toString() + "\n";
  d.setValue( Blas1.zdotu( 2, x, 3, y, 2, ioffx, ioffy ) );
  document.getElementById("debug_textarea").value +=
    "zdotu(2,x,3,y,2) = " + d.toString() + "\n";
  d.setValue( Blas1.zdotu( 2, x, -3, y, -2, ioffx, ioffy ) );
  document.getElementById("debug_textarea").value +=
    "zdotu(2,x,-3,y,-2) = " + d.toString() + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_zdotc() {
  document.getElementById("debug_textarea").value +=
    "testing zdotc *************" + "\n";
  var ioffx = 1;
  var ioffy = 2;
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 6 );
  var i = -1;
  for ( i = 0; i < 6; i ++ ) {
    x[ ioffx + i ] = new Complex( i , i-1 );
  }
  for ( i = 0; i < 6; i ++ ) {
    y[ ioffy + i ] = new Complex( 1 - i , -i );
  }
  var d = new Complex();
  d.setValue( Blas1.zdotc( 0, x, 1, y, 1, ioffx, ioffy ) );
  document.getElementById("debug_textarea").value +=
    "zdotc(0,x,1,y,1) = " + d.toString() + "\n";
  d.setValue( Blas1.zdotc( 6, x, 1, y, 1, ioffx, ioffy ) );
  document.getElementById("debug_textarea").value +=
    "zdotc(6,x,1,y,1) = " + d.toString() + "\n";
  d.setValue( Blas1.zdotc( 2, x, 3, y, 2, ioffx, ioffy ) );
  document.getElementById("debug_textarea").value +=
    "zdotc(2,x,3,y,2) = " + d.toString() + "\n";
  d.setValue( Blas1.zdotc( 2, x, -3, y, -2, ioffx, ioffy ) );
  document.getElementById("debug_textarea").value +=
    "zdotc(2,x,-3,y,-2) = " + d.toString() + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_drot() {
  document.getElementById("debug_textarea").value +=
    "testing drot *************" + "\n";
  var ioffx = 1;
  var ioffy = 2;
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 6 );
  var i = -1;
  for ( i = 0; i < 6; i ++ ) { x[ ioffx + i ] = Number( i ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( 1 - i ); }
  Blas1.drot( 6, x, 1, y, 1 , .8 , .6, ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "drot(6,x,1,y,1,.8,.6), x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2]
        + " " + x[ioffx + 3] + " " + x[ioffx + 4] + " " + x[ioffx + 5]
        + "\n";
  document.getElementById("debug_textarea").value +=
    "drot(6,x,1,y,1,.8,.6), y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2]
        + " " + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5]
        + "\n";

  for ( i = 0; i < 6; i ++ ) { x[ ioffx + i ] = Number( i ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( 1 - i ); }
  Blas1.drot( 2, x, 3, y, 2 , .8 , .6, ioffx, ioffy );
  document.getElementById("debug_textarea").value +=
    "drot(2,x,3,y,2,.8,.6), x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2]
        + " " + x[ioffx + 3] + " " + x[ioffx + 4] + " " + x[ioffx + 5]
        + "\n";
  document.getElementById("debug_textarea").value +=
    "drot(2,x,3,y,2,.8,.6), y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2]
        + " " + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5]
        + "\n";
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function test_drotm() {
  document.getElementById("debug_textarea").value +=
    "testing drotm *************" + "\n";
  var ioffx = 1;
  var ioffy = 2;
  var ioffparam = 3;
  var x = new Array( ioffx + 6 );
  var y = new Array( ioffy + 6 );
  var param = new Array( ioffparam + 5 );
  var i = -1;

  param[ ioffparam ] = -2.;
  for ( i = 0; i < 6; i ++ ) { x[ ioffx + i ] = Number( i ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( 1 - i ); }
  Blas1.drotm( 6, x, 1, y, 1 , param, ioffx, ioffy, ioffparam );
  document.getElementById("debug_textarea").value +=
    "drotm(6,x,1,y,1,param), x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2]
        + " " + x[ioffx + 3] + " " + x[ioffx + 4] + " " + x[ioffx + 5]
        + "\n";
  document.getElementById("debug_textarea").value +=
    "drotm(6,x,1,y,1,param), y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2]
        + " " + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5]
        + "\n";

  for ( i = 0; i < 6; i ++ ) { x[ ioffx + i ] = Number( i ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( 1 - i ); }
  Blas1.drotm( 2, x, 3, y, 2 , param, ioffx, ioffy, ioffparam );
  document.getElementById("debug_textarea").value +=
    "drotm(2,x,3,y,2,param), x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2]
        + " " + x[ioffx + 3] + " " + x[ioffx + 4] + " " + x[ioffx + 5]
        + "\n";
  document.getElementById("debug_textarea").value +=
    "drotm(2,x,3,y,2,param), y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2]
        + " " + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5]
        + "\n";

  param[ ioffparam ] = -1.;
  param[ ioffparam + 1 ] = 3.;
  param[ ioffparam + 2 ] = 2.;
  param[ ioffparam + 3 ] = 5.;
  param[ ioffparam + 4 ] = 4.;
  for ( i = 0; i < 6; i ++ ) { x[ ioffx + i ] = Number( i ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( 1 - i ); }
  Blas1.drotm( 6, x, 1, y, 1 , param, ioffx, ioffy, ioffparam );
  document.getElementById("debug_textarea").value +=
    "drotm(6,x,1,y,1,param), x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2]
        + " " + x[ioffx + 3] + " " + x[ioffx + 4] + " " + x[ioffx + 5]
        + "\n";
  document.getElementById("debug_textarea").value +=
    "drotm(6,x,1,y,1,param), y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2]
        + " " + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5]
        + "\n";

  for ( i = 0; i < 6; i ++ ) { x[ ioffx + i ] = Number( i ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( 1 - i ); }
  Blas1.drotm( 2, x, 3, y, 2 , param, ioffx, ioffy, ioffparam );
  document.getElementById("debug_textarea").value +=
    "drotm(2,x,3,y,2,param), x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2]
        + " " + x[ioffx + 3] + " " + x[ioffx + 4] + " " + x[ioffx + 5]
        + "\n";
  document.getElementById("debug_textarea").value +=
    "drotm(2,x,3,y,2,param), y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2]
        + " " + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5]
        + "\n";

  param[ ioffparam ] = 0.;
  param[ ioffparam + 1 ] = Number.POSITIVE_INFINITY;
  param[ ioffparam + 2 ] = 2.;
  param[ ioffparam + 3 ] = 5.;
  param[ ioffparam + 4 ] = Number.POSITIVE_INFINITY;
  for ( i = 0; i < 6; i ++ ) { x[ ioffx + i ] = Number( i ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( 1 - i ); }
  Blas1.drotm( 6, x, 1, y, 1 , param, ioffx, ioffy, ioffparam );
  document.getElementById("debug_textarea").value +=
    "drotm(6,x,1,y,1,param), x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2]
        + " " + x[ioffx + 3] + " " + x[ioffx + 4] + " " + x[ioffx + 5]
        + "\n";
  document.getElementById("debug_textarea").value +=
    "drotm(6,x,1,y,1,param), y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2]
        + " " + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5]
        + "\n";

  for ( i = 0; i < 6; i ++ ) { x[ ioffx + i ] = Number( i ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( 1 - i ); }
  Blas1.drotm( 2, x, 3, y, 2 , param, ioffx, ioffy, ioffparam );
  document.getElementById("debug_textarea").value +=
    "drotm(2,x,3,y,2,param), x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2]
        + " " + x[ioffx + 3] + " " + x[ioffx + 4] + " " + x[ioffx + 5]
        + "\n";
  document.getElementById("debug_textarea").value +=
    "drotm(2,x,3,y,2,param), y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2]
        + " " + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5]
        + "\n";

  param[ ioffparam ] = 1.;
  param[ ioffparam + 1 ] = 3.;
  param[ ioffparam + 2 ] = Number.POSITIVE_INFINITY;
  param[ ioffparam + 3 ] = Number.POSITIVE_INFINITY;
  param[ ioffparam + 4 ] = 4.;
  for ( i = 0; i < 6; i ++ ) { x[ ioffx + i ] = Number( i ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( 1 - i ); }
  Blas1.drotm( 6, x, 1, y, 1 , param, ioffx, ioffy, ioffparam );
  document.getElementById("debug_textarea").value +=
    "drotm(6,x,1,y,1,param), x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2]
        + " " + x[ioffx + 3] + " " + x[ioffx + 4] + " " + x[ioffx + 5]
        + "\n";
  document.getElementById("debug_textarea").value +=
    "drotm(6,x,1,y,1,param), y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2]
        + " " + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5]
        + "\n";

  for ( i = 0; i < 6; i ++ ) { x[ ioffx + i ] = Number( i ); }
  for ( i = 0; i < 6; i ++ ) { y[ ioffy + i ] = Number( 1 - i ); }
  Blas1.drotm( 2, x, 3, y, 2 , param, ioffx, ioffy, ioffparam );
  document.getElementById("debug_textarea").value +=
    "drotm(2,x,3,y,2,param), x = " + x[ioffx + 0] + " "
        + x[ioffx + 1] + " " + x[ioffx + 2]
        + " " + x[ioffx + 3] + " " + x[ioffx + 4] + " " + x[ioffx + 5]
        + "\n";
  document.getElementById("debug_textarea").value +=
    "drotm(2,x,3,y,2,param), y = " + y[ioffy + 0] + " "
        + y[ioffy + 1] + " " + y[ioffy + 2]
        + " " + y[ioffy + 3] + " " + y[ioffy + 4] + " " + y[ioffy + 5]
        + "\n";
}
