"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.compile = void 0;
function error(linen, line, errMsg) {
    throw "Error at line " + (linen + 1) + ": " + errMsg + "\n " + line.join(' ');
}
function compile(lines) {
    var memoryMap = {};
    var nextFreeMemLocation = 0;
    var assmGen = [];
    var _loop_1 = function (i) {
        var line = lines[i];
        var keyword = line[0];
        if (keyword == 'VAR') {
            var varName = line[1];
            if (memoryMap[varName])
                error(i, line, "Variable " + varName + " already defined.");
            memoryMap[varName] = nextFreeMemLocation++;
        }
        else if (keyword == 'EXPR') {
            var _ = line[0], destVarName = line[1], operator = line[2], operandVarName = line[3];
            if (memoryMap[destVarName] === undefined)
                error(i, line, "Variable " + destVarName + " not defined!");
            if (operator == '=') {
                if (!isNumeric(operandVarName))
                    error(i, line, "Invalid integer value of " + operandVarName + ": not an integer");
                var num = parseInt(operandVarName);
                if (num >= 256 || num < -127)
                    error(i, line, "Invalid integer value of " + operandVarName + ": out of range for 8 bit value");
                assmGen.push(genSet(destVarName, num));
            }
            else {
                var ALUOP = {
                    '+=': 'ADD',
                    '-=': 'SUB',
                    '*=': 'MLT',
                    '|=': 'OR',
                    '&=': 'AND',
                    '^=': 'XOR',
                    '>>=': 'SHR',
                    '<<=': 'SHL'
                }[operator];
                if (!operator)
                    error(i, line, "Non existent operator " + operator);
                if (isNumeric(operandVarName)) {
                    var num = parseInt(operandVarName);
                    if (num >= 256 || num < -127)
                        error(i, line, "Invalid integer value of " + operandVarName + ": out of range for 8 bit value");
                }
                else {
                    if (memoryMap[operandVarName] === undefined)
                        error(i, line, "Variable " + operandVarName + " not defined!");
                }
                assmGen.push(genAluOp(destVarName, ALUOP, operandVarName));
            }
        }
        else if (keyword == "ASM") {
            assmGen.push(function (_, assm) {
                assm.push(line.slice(1).join(' '));
            });
        }
        else if (keyword == "IF") {
            var _ = line[0], operand1 = line[1], operator = line[2], operand2 = line[3];
            if (isNumeric(operand1)) {
                var num = parseInt(operand1);
                if (num >= 256 || num < -127)
                    error(i, line, "Invalid integer value of " + operand1 + ": out of range for 8 bit value");
            }
            else {
                if (memoryMap[operand1] === undefined)
                    error(i, line, "Variable " + operand1 + " not defined!");
            }
            if (isNumeric(operand2)) {
                var num = parseInt(operand2);
                if (num >= 256 || num < -127)
                    error(i, line, "Invalid integer value of " + operand2 + ": out of range for 8 bit value");
            }
            else {
                if (memoryMap[operand2] === undefined)
                    error(i, line, "Variable " + operand2 + " not defined!");
            }
            if (!['==', '>', '>=', '<=', '<', '!='].includes(operator))
                error(i, line, "Unknown compare operator " + operator);
            assmGen.push(genIf(operand1, operator, operand2));
        }
        else if (keyword == 'ELSE') {
            assmGen.push(genElse());
        }
        else if (keyword == 'ENDIF') {
            assmGen.push(genEndIf());
        }
    };
    for (var i = 0; i < lines.length; i++) {
        _loop_1(i);
    }
    var assm = [];
    var ifElseArray = [];
    var helpers = {
        getMemLocation: getMemLocation,
        setIf: setIf,
        setElse: setElse,
        setEndIf: setEndIf
    };
    function getMemLocation(varname) {
        var loc16bit = memoryMap[varname];
        return [
            (loc16bit >> 8) & ((1 << 8) - 1),
            loc16bit & ((1 << 8) - 1)
        ];
    }
    function setIf(loc) {
        ifElseArray.push(['if', loc]);
    }
    function setElse(loc) {
        ifElseArray.push(['else', loc]);
    }
    function setEndIf(loc) {
        ifElseArray.push(['endif', loc]);
    }
    assmGen.forEach(function (fn) {
        fn(helpers, assm);
    });
    var stack = [];
    for (var _i = 0, ifElseArray_1 = ifElseArray; _i < ifElseArray_1.length; _i++) {
        var _a = ifElseArray_1[_i], key = _a[0], pos = _a[1];
        if (key == 'if') {
            stack.push(['if', pos]);
        }
        else if (key == 'else') {
            var _b = stack.pop(), stackTopKey = _b[0], stackTopLoc = _b[1];
            if (stackTopKey != 'if') {
                error(0, [], "ELSE IF UNMATCHED ERROR");
            }
            var elseAddr = (pos + 4) * 2;
            var elseAdrrUpper = (elseAddr >> 8) & ((1 << 8) - 1);
            var elseAddrLower = elseAddr & ((1 << 8) - 1);
            assm[stackTopLoc] = assm[stackTopLoc].replace('MATCHINGELSEUPPER', '' + elseAdrrUpper);
            assm[stackTopLoc + 1] = assm[stackTopLoc + 1].replace('MATCHINGELSELOWER', '' + elseAddrLower);
            stack.push(['else', pos]);
        }
        else if (key == 'endif') {
            var _c = stack.pop(), stackTopKey = _c[0], stackTopLoc = _c[1];
            if (stackTopKey != 'else') {
                error(0, [], "ELSE IF UNMATCHED ERROR");
            }
            var endifAddr = pos * 2;
            var endifAdrrUpper = (endifAddr >> 8) & ((1 << 8) - 1);
            var endifAddrLower = endifAddr & ((1 << 8) - 1);
            assm[stackTopLoc] = assm[stackTopLoc].replace('MATCHINGENDIFUPPER', '' + endifAdrrUpper);
            assm[stackTopLoc + 1] = assm[stackTopLoc + 1].replace('MATCHINGENDIFLOWER', '' + endifAddrLower);
        }
    }
    if (stack.length > 0) {
        error(0, [], "ELSE IF UNMATCHED ERROR");
    }
    return assm;
}
exports.compile = compile;
function isNumeric(value) {
    return /^-?\d+$/.test(value);
}
function genSet(destVarName, num) {
    return function (helpers, assm) {
        var _a = helpers.getMemLocation(destVarName), upper = _a[0], lower = _a[1];
        assm.push("SET 1 " + upper);
        assm.push("SET 2 " + lower);
        assm.push("SET 3 " + num);
        assm.push("SETADDR 1 2");
        assm.push("STORE 3");
    };
}
function genAluOp(destVarName, aluOp, operandVarName) {
    return function (_a, assm) {
        var getMemLocation = _a.getMemLocation;
        var _b = getMemLocation(destVarName), upperDest = _b[0], lowerDest = _b[1];
        if (isNumeric(operandVarName)) {
            assm.push("SET 4 " + operandVarName);
        }
        else {
            var _c = getMemLocation(operandVarName), upperOper = _c[0], lowerOper = _c[1];
            assm.push("SET 1 " + upperOper);
            assm.push("SET 2 " + lowerOper);
            assm.push("SETADDR 1 2");
            assm.push("LOAD 4");
        }
        assm.push("SET 1 " + upperDest);
        assm.push("SET 2 " + lowerDest);
        assm.push("SETADDR 1 2");
        assm.push("LOAD 3");
        assm.push(aluOp + " 3 4");
        assm.push("STORE 3");
    };
}
function genIf(operand1, operator, operand2) {
    return function (_a, assm) {
        var getMemLocation = _a.getMemLocation, setIf = _a.setIf;
        if (isNumeric(operand1)) {
            assm.push("SET 3 " + operand1);
        }
        else {
            var _b = getMemLocation(operand1), upperOper = _b[0], lowerOper = _b[1];
            assm.push("SET 1 " + upperOper);
            assm.push("SET 2 " + lowerOper);
            assm.push("SETADDR 1 2");
            assm.push("LOAD 3");
        }
        if (isNumeric(operand2)) {
            assm.push("SET 4 " + operand2);
        }
        else {
            var _c = getMemLocation(operand2), upperOper = _c[0], lowerOper = _c[1];
            assm.push("SET 1 " + upperOper);
            assm.push("SET 2 " + lowerOper);
            assm.push("SETADDR 1 2");
            assm.push("LOAD 4");
        }
        setIf(assm.length);
        assm.push("SET 1 MATCHINGELSEUPPER");
        assm.push("SET 2 MATCHINGELSELOWER");
        assm.push("SETADDR 1 2");
        if (operator == '>') {
            operator = '<';
            assm.push("SUB 4 3");
        }
        else if (operator == '<=') {
            operator = '>=';
            assm.push("SUB 4 3");
        }
        else {
            assm.push("SUB 3 4");
        }
        var j = {
            '==': 'JNZ',
            '<': 'JNC',
            '>=': 'JC',
            '!=': 'JZ',
        }[operator];
        assm.push(j);
    };
}
function genElse() {
    return function (_a, assm) {
        var getMemLocationm = _a.getMemLocationm, setElse = _a.setElse;
        var _b = ["MATCHINGENDIFUPPER", "MATCHINGENDIFLOWER"], upper = _b[0], lower = _b[1]; // nextEndIfAddr();
        setElse(assm.length);
        assm.push("SET 1 " + upper);
        assm.push("SET 2 " + lower);
        assm.push("SETADDR 1 2");
        assm.push("JMP");
    };
}
function genEndIf() {
    return function (_a, assm) {
        var setEndIf = _a.setEndIf;
        setEndIf(assm.length);
    };
}
//# sourceMappingURL=compiler.js.map