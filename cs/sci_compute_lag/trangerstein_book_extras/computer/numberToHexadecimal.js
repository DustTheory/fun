    function numberToHexadecimal(number) {
      var binary_string=new Array(64);
      for (var b=0;b<64;b++) binary_string[b]=false; // bit string for 0.
      if (number<0) {
        binary_string[0]=true;
        number=-number;
      }
      if (number>0) {
        exponent=0;
//      document.getElementById("debug_textarea").value +=
//        "number = " + number + " exponent = " + exponent + "\n";
        while (number>=2) {
          number *= 0.5;
          exponent++;
        }
        exponent+=1023;
        while (number<1 && exponent>0) {
          number *= 2;
          exponent--;
        }
        var shift=Math.pow(2,51);
        if (exponent>0) shift*=2;
        for (b=0;b<11;b++) {
          binary_string[11-b]=(exponent%2 != 0);
          exponent = Math.floor(exponent/2);
        }
        number = Math.floor(number*shift);
        b=0;
        while (number>0 && b<52) {
          binary_string[63-b]=(number%2 != 0); 
          number = Math.floor(number/2);
          b++;
        }
      }

      hex_string = new String();
      var numerals = "0123456789ABCDEF";
      var power_of_two = new Array(4);
      var p2i=1;
      for (var i=0; i<4; i++) {
        power_of_two[3-i]=p2i;
        p2i *= 2;
      }
      for (var i4=0; i4<64; i4+=4) {
        var n = 0;
        for (i=0; i<4; i++) {
          if (binary_string[i4+i]) n += power_of_two[i];
        }
        hex_string=hex_string+numerals.charAt(n);
      }
      return hex_string;
    }
