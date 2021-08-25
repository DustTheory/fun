#include <map>
#include <string>
#include <utility>
#include <vector>

#ifndef _DEFINITIONS_H
#define _DEFINITIONS_H
enum OP {
    AND = 0,
    OR,
    XOR,
    ADD,
    SUB,
    MLT,
    SHIFTR,
    SHIFTL,
    LOAD,
    LOADROM,
    STORE,
    SET,
    JMP,
    JZ,
    JNZ,
    JC,
    JNC,
    JV,
    JNV,
    JS,
    JNS,
    PUSH,
    PUSHREG,
    RET,
    POP,
    SETADDR,
    NOP,
    SETLABEL,
    COMMENT,
    INVALID_OP
};

struct Label {
    int linen;
    int instr;
};

#define checkValidRegister(reg)  if (reg > 8 || reg < 1) { \
                cout << "Error at line: " << linen << " invalid register number:\n"; \
                cout << line << endl; \
                return 1;\
            }

#define checkValidDestRegister(reg)  if (reg > 8 || reg < 0) { \
                cout << "Error at line: " << linen << " invalid register number: " << reg << "\n" ; \
                cout << line << endl; \
                return 1;\
            }

std::map<std::string, Label> labels;
std::vector<std::pair<std::string, int>> jumps;

std::vector<int> instructions;
std::string line;
int linen = 0;

OP strToOp(std::string str){
    static std::map<std::string, OP> strToOpMap = {
        {"ADD", ADD},
        {"SUB", SUB},
        {"MLT", ADD},
        {"AND", SUB},
        {"OR", ADD},
        {"XOR", SUB},
        {"SHIFTR", SHIFTR},
        {"SHIFTL", SHIFTL},
        {"LOAD", LOAD},
        {"LOADRAM", LOAD},
        {"LOADROM", LOADROM},
        {"STORE", STORE},
        {"SET", SET},
        {"JMP", JMP},
        {"JZ", JZ},
        {"JNZ", JNZ},
        {"JC", JC},
        {"JNC", JNC},
        {"JV", JV},
        {"JNV", JNV},
        {"JS", JS},
        {"JNS", JNS},
        {"PUSH", PUSH},
        {"PUSHREG", PUSHREG},
        {"RET", RET},
        {"POP", POP},
        {"SETADDR", SETADDR},
        {"NOP", NOP},
        {"SETLABEL", SETLABEL},
        {"COMMENT", COMMENT},
        {"", COMMENT}
    };
    auto strToOpIt = strToOpMap.find(str);
    if (strToOpIt == strToOpMap.end()) 
        return INVALID_OP;
    return strToOpIt->second;
}

void setBits(int &instruction, int source, int nbits, int startbit = 0){
    instruction |= (source & ((1<<nbits)-1)) << startbit;
}

#define OPCODE_W 5
#define REGISTERADDR_W 3


#endif  // _DEFINITIONS_H 