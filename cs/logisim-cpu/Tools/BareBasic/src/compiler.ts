import { type } from "os";


function error(linen: number, line: string[], errMsg: string) {
    throw `Error at line ${linen + 1}: ${errMsg}\n ${line.join(' ')}`;
}

type genFn = (getMemLocation: {[key: string]: any} , assm: string[]) => void;

function compile(lines: string[][]): (string | (() => string))[] {

    let memoryMap: { [key: string]: number } = {};
    let nextFreeMemLocation = 0;

    let assmGen: genFn[] = [];

    for (let i = 0; i < lines.length; i++) {
        let line = lines[i];
        let keyword = line[0];
        if (keyword == 'VAR') {
            let varName = line[1];
            if (memoryMap[varName])
                error(i, line, `Variable ${varName} already defined.`);
            memoryMap[varName] = nextFreeMemLocation++;
        }
        else if (keyword == 'EXPR') {
            let [_, destVarName, operator, operandVarName] = line;
            if (memoryMap[destVarName] === undefined)
                error(i, line, `Variable ${destVarName} not defined!`);

            if (operator == '=') {
                if (!isNumeric(operandVarName))
                    error(i, line, `Invalid integer value of ${operandVarName}: not an integer`);
                let num = parseInt(operandVarName);
                if (num >= 256 || num < -127)
                    error(i, line, `Invalid integer value of ${operandVarName}: out of range for 8 bit value`);
                assmGen.push(genSet(destVarName, num));
            } else {
                let ALUOP = {
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
                    error(i, line, `Non existent operator ${operator}`);

                if (isNumeric(operandVarName)) {
                    let num = parseInt(operandVarName);
                    if (num >= 256 || num < -127)
                        error(i, line, `Invalid integer value of ${operandVarName}: out of range for 8 bit value`);
                } else {
                    if (memoryMap[operandVarName] === undefined)
                        error(i, line, `Variable ${operandVarName} not defined!`);
                }
                assmGen.push(genAluOp(destVarName, ALUOP, operandVarName));
            }
        } else if (keyword == "ASM") {
            assmGen.push((_, assm) => {
                assm.push(line.slice(1).join(' '));
            });
        } else if (keyword == "IF") {
            let [_, operand1, operator, operand2] = line;

            if (isNumeric(operand1)) {
                let num = parseInt(operand1);
                if (num >= 256 || num < -127)
                    error(i, line, `Invalid integer value of ${operand1}: out of range for 8 bit value`);
            } else {
                if (memoryMap[operand1] === undefined)
                    error(i, line, `Variable ${operand1} not defined!`);
            }

            if (isNumeric(operand2)) {
                let num = parseInt(operand2);
                if (num >= 256 || num < -127)
                    error(i, line, `Invalid integer value of ${operand2}: out of range for 8 bit value`);
            } else {
                if (memoryMap[operand2] === undefined)
                    error(i, line, `Variable ${operand2} not defined!`);
            }

            if (!['==', '>', '>=', '<=', '<', '!='].includes(operator))
                error(i, line, `Unknown compare operator ${operator}`);

            assmGen.push(genIf(operand1, operator, operand2));
        } else if(keyword == 'ELSE') {
            assmGen.push(genElse());
        } else if(keyword == 'ENDIF'){
            assmGen.push(genEndIf());
        }
    }

    let assm: string[] = [];
    let ifElseArray : [('if' | 'else' | 'endif'), number][] = [];

    let helpers = {
        getMemLocation: getMemLocation,
        setIf: setIf,
        setElse: setElse,
        setEndIf: setEndIf
    };

    function getMemLocation(varname: string): [number, number] {
        let loc16bit = memoryMap[varname];
        return [
            (loc16bit >> 8) & ((1 << 8) - 1),
            loc16bit & ((1 << 8) - 1)
        ];
    }

    function setIf(loc: number){
        ifElseArray.push(['if', loc]);
    }

    function setElse(loc: number){
        ifElseArray.push(['else', loc]);
    }

    function setEndIf(loc: number){
        ifElseArray.push(['endif', loc]);
    }

    assmGen.forEach((fn: genFn) => {
        fn(helpers, assm);
    });

    let stack :  [('if' | 'else' | 'endif'), number][] = [];

    for(let [key, pos] of ifElseArray){
        if(key == 'if'){
            stack.push(['if', pos]);
        }else if(key == 'else'){
            let [stackTopKey, stackTopLoc] = stack.pop();
            if(stackTopKey != 'if'){
                error(0, [], "ELSE IF UNMATCHED ERROR");
            }
            let elseAddr = (pos + 4)*2;
            let elseAdrrUpper = (elseAddr >> 8) & ((1 << 8) - 1);
            let elseAddrLower = elseAddr & ((1 << 8) - 1);
            assm[stackTopLoc] = assm[stackTopLoc].replace('MATCHINGELSEUPPER', ''+elseAdrrUpper);
            assm[stackTopLoc+1] = assm[stackTopLoc+1].replace('MATCHINGELSELOWER', ''+elseAddrLower);
            stack.push(['else', pos]);
        }else if(key == 'endif'){
            let [stackTopKey, stackTopLoc] = stack.pop();
            if(stackTopKey != 'else'){
                error(0, [], "ELSE IF UNMATCHED ERROR");
            }
            let endifAddr = pos*2;
            let endifAdrrUpper = (endifAddr >> 8) & ((1 << 8) - 1);
            let endifAddrLower = endifAddr & ((1 << 8) - 1);
            assm[stackTopLoc] = assm[stackTopLoc].replace('MATCHINGENDIFUPPER', ''+endifAdrrUpper);
            assm[stackTopLoc+1] = assm[stackTopLoc+1].replace('MATCHINGENDIFLOWER', ''+endifAddrLower);
        }
    }

    if(stack.length > 0){
        error(0, [], "ELSE IF UNMATCHED ERROR");
    }

    return assm;
}

export {
    compile
}

function isNumeric(value: string): boolean {
    return /^-?\d+$/.test(value);
}

function genSet(destVarName: string, num: number): genFn {
    return (helpers: any, assm: string[]) => {
        let [upper, lower] = helpers.getMemLocation(destVarName);
        assm.push("SET 1 " + upper);
        assm.push("SET 2 " + lower);
        assm.push("SET 3 " + num);
        assm.push("SETADDR 1 2");
        assm.push("STORE 3");
    };
}

function genAluOp(destVarName: string, aluOp: string, operandVarName: string): genFn {
    return ({ getMemLocation }, assm: string[]) => {
        let [upperDest, lowerDest] = getMemLocation(destVarName);

        if (isNumeric(operandVarName)) {
            assm.push("SET 4 " + operandVarName);
        } else {
            let [upperOper, lowerOper] = getMemLocation(operandVarName);
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

function genIf(operand1: string, operator: string, operand2: string): genFn {
    return ({ getMemLocation, setIf }, assm: string[]) => {
        if (isNumeric(operand1)) {
            assm.push("SET 3 " + operand1);
        } else {
            let [upperOper, lowerOper] = getMemLocation(operand1);
            assm.push("SET 1 " + upperOper);
            assm.push("SET 2 " + lowerOper);
            assm.push("SETADDR 1 2");
            assm.push("LOAD 3");
        }

        if (isNumeric(operand2)) {
            assm.push("SET 4 " + operand2);
        } else {
            let [upperOper, lowerOper] = getMemLocation(operand2);
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
        } else if (operator == '<=') {
            operator = '>=';
            assm.push("SUB 4 3");
        }
        else {
            assm.push("SUB 3 4");
        }

        let j = {
            '==': 'JNZ',
            '<': 'JNC',
            '>=': 'JC',
            '!=': 'JZ',
        }[operator];

        assm.push(j);
    };
}

function genElse(): genFn {
    return ({ getMemLocationm, setElse }, assm: string[]) => {
        let [upper, lower] = ["MATCHINGENDIFUPPER", "MATCHINGENDIFLOWER"]; // nextEndIfAddr();
        setElse(assm.length);
        assm.push("SET 1 " + upper);
        assm.push("SET 2 " + lower);
        assm.push("SETADDR 1 2");
        assm.push("JMP");
    };
}


function genEndIf(): genFn {
    return ({ setEndIf }, assm: string[]) => {
        setEndIf(assm.length);
    };
}