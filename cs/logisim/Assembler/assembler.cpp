#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>
#include <limits.h>
#include <iomanip>
#include <bitset>

#include "definitions.h";

using namespace std;

int main(int argc, char *argv[]) {

    if (argc < 2) {
        cout << "Invalid arguments! Valid argument format" << endl;
        cout << "microcodeAssembler INPUT_FILE_PATH OUTPUT_FILE_PATH" << endl;
    }
    ifstream in(argv[1]);
    ofstream out(argv[2], ios::out);

    //ifstream in("input.asm");
    //ofstream out("output.bin", ios::out);
    out << "v3.0 hex bytes plain big-endian\n";

    while (in.peek() != EOF) {
        int instruction = 0;
        getline(in, line);
        stringstream linestream(line);

        string opstr;
        linestream >> opstr;
        OP op = strToOp(opstr);
        setBits(instruction, op, OPCODE_W, 0);
        int reg1, reg2, reg, data;
        switch (op) {
            linen++;
        // Simple ALU OPs
        case AND:
        case OR:
        case XOR:
        case ADD:
        case SUB:
        case MLT:
        case SHIFTR:
        case SHIFTL:
            int reg1, reg2;
            linestream >> reg1 >> reg2;
            checkValidDestRegister(reg1);
            checkValidRegister(reg2);
            setBits(instruction, reg1, REGISTERADDR_W, OPCODE_W);
            setBits(instruction, reg2, REGISTERADDR_W, OPCODE_W + REGISTERADDR_W);
            break;
        case LOAD:
        case LOADROM:
            linestream >> reg;
            checkValidDestRegister(reg);
            setBits(instruction, reg, REGISTERADDR_W, OPCODE_W);
            break;
        case STORE:
            linestream >> reg;
            checkValidRegister(reg);
            setBits(instruction, reg, REGISTERADDR_W, OPCODE_W);
            break;
        case SET:
            linestream >> reg >> data;
            checkValidDestRegister(reg);
            if(data > UCHAR_MAX || data < CHAR_MIN) {
                cout << "Error at line: " << linen << " invalid data parameter:\n";
                cout << line << endl;
                return 1;
            }
            setBits(instruction, reg, REGISTERADDR_W, OPCODE_W);
            setBits(instruction, data, 8, OPCODE_W+REGISTERADDR_W);
            cout << data << endl;
            break;
        case JMP:
        case JZ:
        case JNZ:
        case JC:
        case JNC:
        case JV:
        case JNV:
        case JS:
        case JNS:
            break;
        case PUSH:
            break;
        case PUSHREG:
            linestream >> reg1 >> reg2;
            checkValidRegister(reg1);
            checkValidRegister(reg2);
            setBits(instruction, reg, REGISTERADDR_W, OPCODE_W);
            setBits(instruction, reg, REGISTERADDR_W, OPCODE_W+REGISTERADDR_W);
            break;
        case RET:
            break;
        case POP:
            linestream >> reg1 >> reg2;
            checkValidDestRegister(reg1);
            checkValidDestRegister(reg2);
            setBits(instruction, reg1, REGISTERADDR_W, OPCODE_W);
            setBits(instruction, reg2, REGISTERADDR_W, OPCODE_W+REGISTERADDR_W);
            break;
        case SETADDR:
            linestream >> reg1 >> reg2;
            checkValidDestRegister(reg1);
            checkValidDestRegister(reg2);
            setBits(instruction, reg1, REGISTERADDR_W, OPCODE_W);
            setBits(instruction, reg2, REGISTERADDR_W, OPCODE_W+REGISTERADDR_W);
            break;
        case COMMENT:
            continue;
        case INVALID_OP:
            cout << "Error at line: " << linen << " invalid instruction keyword : " << opstr << endl;
            return 1;
        case SETLABEL:
            break;

        default:
            cout << "Unknown error at line: " << linen << endl;
            cout << line << endl;
            return 1;
        }
        instructions.push_back(instruction);
    }

    for(int i = 0; i < instructions.size(); i++) {
        cout << bitset<16>(instructions[i]) << endl;
        out << setfill('0') << setw(4) << hex << (instructions[i] & 0xffff);
    }

}