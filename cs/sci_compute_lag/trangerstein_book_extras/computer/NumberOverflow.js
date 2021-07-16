function NumberOverflow() {
  var a = 1;
  while (isFinite(a)) {
//  toString first converts a to integer, then computes hex
//  document.write(a + " " + a.toString(16)+ "<br>");
    document.write(a + "<br>");
    a=2*a+1;
  }
}
