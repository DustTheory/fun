//Flanagan p. 217
function enumeration( namesToValues ) {
  var numeration = function() { throw "Cannot instantiate enumerations"; };
  var proto = enumeration.prototype = {
    constructor: enumeration,
    toString: function() { return this.name; },
    valueOf: function() { return this.value; },
    toJSON: function() { return this.name; }
  };
  enumeration.values = [];
  for ( name in namesToValues ) {
    var e = inherit( proto);
    var e = new proto;
    e.name = name;
    e.value = namesToValues[ name ];
    enumeration[ name ] = e;
    enumeration.values.push(e);
  }
  enumeration.foreach = function(f,c) {
    for ( var i = 0; i < this.value.length; i ++ ) {
      f.call( c, this.values[i] );
    }
  };
  return enumeration;
}
