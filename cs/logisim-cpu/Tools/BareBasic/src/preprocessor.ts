import * as fs from 'fs';
import * as path from 'path';

function preprocess(rawCode: string, pathRoot: string, depth: number = 0): string[][] {
    if(depth > 10){
        throw "Max include depth of 10 exceded! This is done to prevent include loops so check for that."
    }
    let lines = rawCode.split('\n').map(line => line.trim().split(' '));
    for(let i = 0; i < lines.length; i++){
        let line = lines[i];
        if(line[0] != 'INCLUDE')
            continue;
        let src = line[1];
        let rawCode = fs.readFileSync(path.resolve(pathRoot, src), { encoding: 'utf-8', flag: 'r'});
        let includeLines = preprocess(rawCode, path.resolve(pathRoot, src), depth+1);
        lines = lines.slice(0, i).concat(includeLines, lines.slice(i+1));
        i += includeLines.length;
    }
    return lines;
}

export {
    preprocess
}
