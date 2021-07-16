package computer {
  import flash.display.Sprite;
  import flash.text.*;
  import etirps.ClickSensor;
  import etirps.VScrollText;
  public class intOverflow extends Sprite {
    internal static var t:TextField = null;
    internal static var vst:VScrollText = null;
    internal static var format:TextFormat;
    public function intOverflow() {
      try {
        format = new TextFormat();
        format.size = 15;

        var int_overflow_label:TextField = new TextField();
        int_overflow_label.defaultTextFormat = format;
        int_overflow_label.htmlText = "overflow int";
        var cs:ClickSensor = new ClickSensor( this, int_overflow_label,
          myintOverflow, 0, 350 );
        addChild(cs);

        var int_negative_overflow_label:TextField = new TextField();
        int_negative_overflow_label.defaultTextFormat = format;
        int_negative_overflow_label.htmlText = "overflow neg int";
        var cs1:ClickSensor = new ClickSensor( this,
          int_negative_overflow_label, intNegativeOverflow, 100, 350 );
        addChild(cs1);

        var uint_overflow_label:TextField = new TextField();
        uint_overflow_label.defaultTextFormat = format;
        uint_overflow_label.htmlText = "overflow uint";
        var cs2:ClickSensor = new ClickSensor( this,
          uint_overflow_label, uintOverflow, 220, 350 );
        addChild(cs2);
      } catch (e:Error) {
        trace(e.message);
      }
    }
    public static function myintOverflow( s:Sprite ):void {
//    trace("myintOverflow: vst = " + vst);
      if ( t!=null ) s.removeChild( t );
      if (vst != null) vst.cleanup();

      var j:int = 1;
      t = new TextField();
      t.defaultTextFormat = format;
      t.multiline = true;
      t.htmlText = "compute j_n = 2*j_{n-1} + 1, j_0 = 1 while j_n > 0";
      while (j>0) {
        t.htmlText += j.toString() + " = hex 0x" + j.toString(16);
        j = 2*j + 1;
      }
      t.htmlText += j.toString() + " = hex 0x" + j.toString(16);

      t.autoSize = TextFieldAutoSize.NONE;
      t.width = t.textWidth+10;
      t.height = 300;
      s.addChild(t);

      vst = new VScrollText( t, s );
    }
    public static function intNegativeOverflow( s:Sprite ):void {
//    trace("intNegativeOverflow: vst = " + vst);
      if (t!=null) s.removeChild( t );
      if (vst != null) vst.cleanup();

      var j:int = -1;
      t = new TextField();
      t.defaultTextFormat = format;
      t.multiline = true;
      t.htmlText =
        "compute j_n = 2*j_{n-1} - 1, j_0 = -1 while j_n negative";
      while (j<0) {
        t.htmlText += j.toString() + " = hex 0x" + j.toString(16);
        j = 2*j - 1;
      }
      t.htmlText += j.toString() + " = hex 0x" + j.toString(16);

      t.autoSize = TextFieldAutoSize.NONE;
      t.width = t.textWidth + 10;
      t.height = 300;
      s.addChild(t);

      vst = new VScrollText( t, s );
    }
    public static function uintOverflow( s:Sprite ):void {
//    trace("uintOverflow: vst = " + vst);
      if (t!=null) s.removeChild( t );
      if (vst != null) vst.cleanup();

      var jold:uint = 0;
      var jnew:uint = 1;
      t = new TextField();
      t.defaultTextFormat = format;
      t.multiline = true;
      t.htmlText =
        "compute j_n = 2*j_{n-1} + 1, j_0 = 0 while j_n > u_{n-1}";
      while ( jnew > jold ) {
        jold = jnew;
        t.htmlText += jnew.toString() + " = hex 0x" + jnew.toString(16);
        jnew = 2*jnew + 1;
      }
      t.htmlText += jnew.toString() + " = hex 0x" + jnew.toString(16);

      t.autoSize = TextFieldAutoSize.NONE;
      t.width = t.textWidth+10;
      t.height = 300;
      s.addChild(t);

      vst = new VScrollText( t, s );
    }
  }
}
