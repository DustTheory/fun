"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.preprocess = void 0;
var fs = require("fs");
var path = require("path");
function preprocess(rawCode, pathRoot, depth) {
    if (depth === void 0) { depth = 0; }
    if (depth > 10) {
        throw "Max include depth of 10 exceded! This is done to prevent include loops so check for that.";
    }
    var lines = rawCode.split('\n').map(function (line) { return line.trim().split(' '); });
    for (var i = 0; i < lines.length; i++) {
        var line = lines[i];
        if (line[0] != 'INCLUDE')
            continue;
        var src = line[1];
        var rawCode_1 = fs.readFileSync(path.resolve(pathRoot, src), { encoding: 'utf-8', flag: 'r' });
        var includeLines = preprocess(rawCode_1, path.resolve(pathRoot, src), depth + 1);
        lines = lines.slice(0, i).concat(includeLines, lines.slice(i + 1));
        i += includeLines.length;
    }
    return lines;
}
exports.preprocess = preprocess;
//# sourceMappingURL=preprocessor.js.map