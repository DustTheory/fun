#include <iostream>
#include <fstream>
#include <map>
#include <bitset>
#include <sstream>
#include <iomanip>
#include <vector>

using namespace std;

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
    NOP,
    SETLABEL,
    COMMENT,
    EXEC
};

enum SOURCE {
    NONE = 0,
    RAMREAD,
    RAMWRITE,
    ROMREAD,
    CUADDRESSREAD,
    RAMREADPC,
    ROMREADPC,
    EXEC_INSTR
};

enum SELECT {
    FALSE = 0,
    TRUE,
    C,
    NOTC,
    Z,
    NOTZ,
    V,
    NOTV,
    S
};

void setAreg(int areg, int &instruction) {
    instruction |= areg & 0xf;
}

void setBreg(int breg, int &instruction) {
    instruction |= (breg & 0xf) << 4;
}

void setDreg(int dreg, int &instruction) {
    instruction |= (dreg & 0xf) << 8;
}

void setALUOP(int aluop, int &instruction) {
    instruction |= (aluop & 0x07) << 12;
}

void setSource(SOURCE src, int &instruction) {
    instruction |= (src & 0x07) << 15;
    if(src == EXEC_INSTR)
        cout << bitset<29>(instruction) << endl;
}

void setSelect(SELECT select, int &instruction) {
    instruction |= (select & 0x07) << 18;

}

void setAddr(int addr, int &instruction) {
    instruction |= (addr & 0xff) << 21;
}

#define checkValidRegister(reg)  if (reg > 15 || reg < 0) { \
                cout << "Error at line: " << linen+1 << " invalid register number:\n"; \
                cout << line << endl; \
                return 1;\
            }

#define checkValidDestRegister(reg)  if (reg > 15 || reg < 0) { \
                cout << "Error at line: " << linen+1 << " invalid register number: " << reg << "\n" ; \
                cout << line << endl; \
                return 1;\
            }


struct Label {
    int linen;
    int instr;
};


int main(int argc, char *argv[]) {
    map<string, Label> labels;
    vector<pair<string, int>> jumps;
    map<string, OP> strToOp {
        {"ADD", ADD},
        {"SUB", SUB},
        {"MLT", ADD},
        {"AND", SUB},
        {"OR", ADD},
        {"XOR", SUB},
        {"SHIFTR", SHIFTR},
        {"SHIFTL", SHIFTL},
        {"LOAD", LOAD},
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
        {"NOP", NOP},
        {"SETLABEL", SETLABEL},
        {"COMMENT", COMMENT},
        {"", COMMENT},
        {"EXEC", EXEC}
    };

    if (argc < 2) {
        cout << "Invalid arguments! Valid argument format" << endl;
        cout << "microcodeAssembler INPUT_FILE_PATH OUTPUT_FILE_PATH" << endl;
    }
    ifstream in(argv[1]);
    ofstream out(argv[2], ios::out);

    // ifstream in("rom.asm");
    // ofstream out("rom.bin", ios::out);
    out << "v3.0 hex words addressed\n";
    int code[256] = {};
    string line;
    int instrcnt = 0;
    int linen = 0;
    while (in.peek() != EOF) {
        if (instrcnt >= 256) {
            cout << "Program too long. Please remain under 256 instructions" << endl;
            return 1;
        }
        int instruction = 0;
        getline(in, line);
        stringstream linestream(line);

        string cmd;
        linestream >> cmd;
        auto strToOpIt = strToOp.find(cmd);
        if (strToOpIt == strToOp.end()) {
            cout << "Error at line: " << linen+1 << " invalid instruction keyword : " << cmd << endl;
            return 1;
        }
        OP op = strToOpIt->second;
        if (op >= AND && op <= SHIFTL) {
            int regA, regB, regD;
            linestream >> regA >> regB >> regD;
            checkValidRegister(regA);
            checkValidRegister(regB);
            checkValidDestRegister(regD);
            setAreg(regA, instruction);
            setBreg(regB, instruction);
            setDreg(regD, instruction);
            setALUOP(op, instruction);
        } else if(op == SET) {
            int dest, data;
            linestream >> dest;
            linestream >> data;
            checkValidDestRegister(dest);
            setDreg(dest, instruction);
            setSource(CUADDRESSREAD, instruction);
            setAddr(data, instruction);
        } else if(op == LOAD) {
            int dest;
            linestream >> dest;
            string device;
            linestream >> device;
            checkValidDestRegister(dest);
            setDreg(dest, instruction);
            if(device == "ROM")
                setSource(ROMREAD, instruction);
            else if(device == "ROMPC")
                setSource(ROMREADPC, instruction);
            else if(device == "RAMPC")
                setSource(RAMREADPC, instruction);
            else
                setSource(RAMREAD, instruction);

        } else if(op == STORE) {
            int src;
            linestream >> src;
            checkValidRegister(src);
            setAreg(src, instruction);
            setBreg(0, instruction);
            setALUOP(OR, instruction);
            setSource(RAMWRITE, instruction);
        }  else if(op == EXEC) {
            setSource(EXEC_INSTR, instruction);
        } else if(op == NOP) {
            // nothing ig
        } else if(op == SETLABEL) {
            // not a real instruction
            string label;
            linestream >> label;
            auto tmp = labels.find(label);
            if(tmp != labels.end()) {
                cout << "Error at line: " << linen+1 << " dublicate label declaration : " << label << endl;
                cout << "Same label previously set on line: " << tmp->second.linen << endl;
                return 1;
            } else {
                labels[label] = {linen, instrcnt};
            }
            linen++;
            continue;
        } else if(op == COMMENT) {
            linen++;
            continue;
        } else {
            cout << "Instruction not implemented: " << cmd << endl;
            return 1;
        }
        // conditionals and jumps after each instruction
        string conditionalType;
        string label;
        linestream >> conditionalType >> label;
        strToOpIt = strToOp.find(conditionalType);

        SELECT select = FALSE;
        if (conditionalType != "") {
            if (strToOpIt == strToOp.begin() || strToOpIt->second < JMP || strToOpIt->second > JNS) {
                cout << "Error at line: " << linen+1 << " invalid conditional type : " << conditionalType << endl;
                return 1;
            }
            select = (SELECT)((int)strToOpIt->second - (int)JMP + 1);
            jumps.push_back({label, instrcnt});
        }
        setSelect(select, instruction);

        code[instrcnt] = instruction;
        linen++;
        instrcnt++;
    }

    for(auto p : jumps) {
        auto it = labels.find(p.first);
        if(it == labels.end()) {
            cout << "Error at line: " << it->second.linen+1 << " non existent label specified : " << p.first << endl;
            return 1;
        } else {
            setAddr(it->second.instr, code[p.second]);
        }
    }

    {
        // generate setaddrpla.txt file
        auto it = labels.find("SETADDR_PLA_MARKER");
        if(it != labels.end()) {
            ofstream setaddrplafile("setaddrpla.txt", ios::out);
            setaddrplafile << "# Logisim PLA program table\n";
            int addr = it->second.instr;
            for(int i = 0; i < 8; i++) {
                for(int j = 0; j < 8; j++) {
                    setaddrplafile << bitset<3>(j)<<bitset<3>(i) << ' ' << bitset<8>(addr) << endl;
                    addr += 2;
                }
            }
            setaddrplafile.close();
        }

        // print jump instruction implementation memory location
        it= labels.find("JUMP_INSTRUCTION_IMPLEMENTATION");
        if(it == labels.end())
            cout << "JMP instruction implementation not found" << endl;
        else
            cout << "JMP instruction implementation address: " << it->second.instr << " hex: " << hex << it->second.instr << dec << endl;
    }

    for (int i = 0; i < 256; i++) {
        if(i < 20)
            cout << setfill('0') << setw(2) << hex << i << ' ' << bitset<29>(code[i]) << endl;
        if (i % 8 == 0)
            out << setfill('0') << setw(2) << hex << i << ": ";
        out << setfill('0') << setw(2) << hex << code[i] << ' ';
        if ((i + 1) % 8 == 0)
            out << endl;
    }
    in.close();
    out.close();
}