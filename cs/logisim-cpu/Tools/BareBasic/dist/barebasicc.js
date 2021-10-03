"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
var fs = require("fs");
var path = require("path");
var commandLineArgs = require("command-line-args");
var commandLineUsage = require("command-line-usage");
var process_1 = require("process");
var preprocessor_1 = require("./preprocessor");
var compiler_1 = require("./compiler");
var optionDefinitions = [
    { name: 'src', type: String, multiple: true, defaultOption: true },
    { name: 'out', alias: 'o', type: String, defaultValue: './a.out' },
    { name: 'help', alias: 'h', type: Boolean, defaultValue: false },
];
var sections = [
    {
        header: 'BareBasic Compiler',
        content: 'A very simple compiler for a very basic language for a very fucked cpu'
    },
    {
        header: 'Synopsis',
        content: '$ barebasicc input-file-path ...'
    },
    {
        header: 'Options',
        optionList: [
            {
                name: 'src',
                typeLabel: '{underline file}',
                description: 'The input file to compile.'
            },
            {
                name: 'out',
                alias: 'o',
                typeLabel: '{underline output filepath}',
                description: 'The input file to compile.'
            },
            {
                name: 'help',
                alias: 'h',
                description: 'Print this usage guide.'
            }
        ]
    }
];
var options = commandLineArgs(optionDefinitions);
var usage = commandLineUsage(sections);
if (options.help)
    console.log(usage);
if (!options.src) {
    console.error("Source file path not specified. Use --help or -h to read usage guide.");
    (0, process_1.exit)(1);
}
var sourceFilePath = options.src[0];
var sourceFileDir = path.parse(sourceFilePath).dir;
var rawCode = fs.readFileSync(path.resolve(__dirname, sourceFilePath), { encoding: 'utf8', flag: 'r' });
var preprocessedCode = (0, preprocessor_1.preprocess)(rawCode, path.resolve(__dirname, sourceFileDir));
var assemblyCode = (0, compiler_1.compile)(preprocessedCode);
fs.writeFileSync(options.out, assemblyCode.join('\n'));
//# sourceMappingURL=barebasicc.js.map