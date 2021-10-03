import * as fs from 'fs';
import * as path from 'path';

import * as commandLineArgs from 'command-line-args';
import * as commandLineUsage from 'command-line-usage';
import { exit } from 'process';

import { preprocess } from './preprocessor';
import { compile } from './compiler';

const optionDefinitions = [
    { name: 'src', type: String, multiple: true, defaultOption: true},
    { name: 'out', alias: 'o', type: String, defaultValue: './a.out'},
    { name: 'help', alias: 'h', type: Boolean, defaultValue: false },
];
const sections = [
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

const options = commandLineArgs(optionDefinitions);
const usage = commandLineUsage(sections);

if(options.help)
  console.log(usage);

if(!options.src) {
  console.error("Source file path not specified. Use --help or -h to read usage guide.");
  exit(1);
}

let sourceFilePath = options.src[0];
let sourceFileDir = path.parse(sourceFilePath).dir;
let rawCode = fs.readFileSync(path.resolve(__dirname, sourceFilePath), {encoding:'utf8', flag:'r'});

let preprocessedCode = preprocess(rawCode, path.resolve(__dirname,sourceFileDir));
let assemblyCode = compile(preprocessedCode);

fs.writeFileSync(options.out, assemblyCode.join('\n'));