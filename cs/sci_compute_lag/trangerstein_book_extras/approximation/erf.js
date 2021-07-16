//from http://www.netlib.org/specfun/erf
//See W.J. Cody, "Rational Chebyshev Approximations for the Error Function".
//  Math. Comp. 23 (107): 631-637
function calerf(arg,jint) {
//document.getElementById("debug_textarea").value +=
//  "entering calerf\n";
//document.getElementById("debug_textarea").value +=
//  "arg,jint = " + arg + " " + jint + "\n";
  var sqrpi = 1/Math.sqrt( Math.PI );
  var thresh = .46875;
  var xbig = 26.543;
  var xhuge = 6.71e7;
  var xinf = 1.79e308;
  var xmax = 2.53e307;
  var xneg = -26.628;
  var xsmall = 1.11e-16;
  var a = [ 3.16112374387056560, 1.13864154151050156e2,
            3.77485237685302021e2, 3.20937758913846947e3,
            1.85777706184603153e-1 ];
  var b = [ 2.36012909523441209e1, 2.44024637934444173e2,
            1.28261652607737228e3, 2.84423683343917062e3 ];
  var c = [ 5.64188496988670089e-1, 8.88314979438837594,
            6.61191906371416295e1, 2.98635138197400131e2,
            8.81952221241769090e2, 1.71204761263407058e3,
            2.05107837782607147e3, 1.23033935479799725e3,
            2.15311535474403846e-8 ];
  var d = [ 1.57449261107098347e1,1.17693950891312499e2,
            5.37181101862009858e2,1.62138957456669019e3,
            3.29079923573345963e3,4.36261909014324716e3,
            3.43936767414372164e3,1.23033935480374942e3 ];
  var p = [ 3.05326634961232344e-1,3.60344899949804439e-1,
            1.25781726111229246e-1,1.60837851487422766e-2,
            6.58749161529837803e-4,1.63153871373020978e-2 ];
  var q = [ 2.56852019228982242e0,1.87295284992346047e0,
            5.27905102951428412e-1,6.05183413124413191e-2,
            2.33520497626869185e-3 ];
  var x = arg;
  var y = Math.abs(x);
  var goto300 = false;
  if (y <= thresh) { // compute erf(y)
    var ysq = 0;
    if (y > xsmall) ysq = y * y;
    var xnum = a[4]*ysq;
    var xden = ysq;
    for ( var i = 0; i < 3; i ++ ) {
      xnum = (xnum + a[i]) * ysq;
      xden = (xden + b[i]) * ysq;
    }
    var result = x * (xnum + a[3]) / (xden + b[3]);
    if (jint != 0) result = 1 - result;
    if (jint == 2) result = Math.exp(ysq) * result;
    return result;
  } else if (y <= 4) { // compute erfc(y) = 1 - erf(y)
    var xnum = c[8]*y;
    var xden = y;
    for ( var i = 0; i < 7; i ++ ) {
      xnum = (xnum + c[i]) * y;
      xden = (xden + d[i]) * y;
    }
    var result = (xnum + c[7]) / (xden + d[7]);
    if (jint != 2) {
      var ysq = Math.floor(y*16)/16;
      var del = (y-ysq)*(y+ysq);
      result = Math.exp(-ysq*ysq) * Math.exp(-del) * result;
    }
  } else { // compute erfc(y) = 1 - erf(y)
    result = 0;
    if (y >= xbig) {
      if ((jint != 2) || (y >= xmax)) goto300 = true;
      if ( ! goto300 && y >= xhuge) {
        result = sqrpi / y;
        goto300 = true;
      }
    }
    if ( ! goto300 ) {
      ysq = 1 / (y * y) ;
      var xnum = p[5]*ysq;
      var xden = ysq;
      for ( var i = 0; i < 4; i ++ ) {
        xnum = (xnum + p[i]) * ysq;
        xden = (xden + q[i]) * ysq;
      }
      var result = ysq *(xnum + p[4]) / (xden + q[4]);
      result = (sqrpi -  result) / y;
      if (jint != 2) {
        ysq = Math.floor(y*16)/16;
        del = (y-ysq)*(y+ysq);
        result = Math.exp(-ysq*ysq) * Math.exp(-del) * result;
      }
    }
  }
  if (jint == 0) { // erf
    result = (.5 - result) + .5;
    if (x < 0) result = -result;
  } else if (jint == 1) { // erfc
    if (x < 0) result = 2 - result;
  } else { // erfcx
    if (x < 0) {
      if (x < xneg) {
        result = xinf;
      } else {
        ysq = Math.ceil(x*16)/16;
        del = (x-ysq)*(x+ysq);
        y = Math.exp(ysq*ysq) * Math.exp(del);
        result = (y+y) - result;
      }
    }
  }
//document.getElementById("debug_textarea").value +=
//  "leaving calerf\n";
  return result;
}

function erf(x) { return calerf(x,0) }
function erfc(x) { return calerf(x,1); }
function erfcx(x) { return calerf(x,2); }
