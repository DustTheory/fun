var InstructionOpcode;
(function (InstructionOpcode) {
    InstructionOpcode[InstructionOpcode["AND"] = 0] = "AND";
    InstructionOpcode[InstructionOpcode["OR"] = 1] = "OR";
    InstructionOpcode[InstructionOpcode["XOR"] = 2] = "XOR";
    InstructionOpcode[InstructionOpcode["ADD"] = 3] = "ADD";
    InstructionOpcode[InstructionOpcode["SUB"] = 4] = "SUB";
    InstructionOpcode[InstructionOpcode["MLT"] = 5] = "MLT";
    InstructionOpcode[InstructionOpcode["SHL"] = 6] = "SHL";
    InstructionOpcode[InstructionOpcode["SHR"] = 7] = "SHR";
    InstructionOpcode[InstructionOpcode["LOADRAM"] = 8] = "LOADRAM";
    InstructionOpcode[InstructionOpcode["LOADROM"] = 9] = "LOADROM";
    InstructionOpcode[InstructionOpcode["STORE"] = 10] = "STORE";
    InstructionOpcode[InstructionOpcode["SET"] = 11] = "SET";
    InstructionOpcode[InstructionOpcode["JMP"] = 12] = "JMP";
    InstructionOpcode[InstructionOpcode["JZ"] = 13] = "JZ";
    InstructionOpcode[InstructionOpcode["JNZ"] = 14] = "JNZ";
    InstructionOpcode[InstructionOpcode["JC"] = 15] = "JC";
    InstructionOpcode[InstructionOpcode["JNC"] = 16] = "JNC";
    InstructionOpcode[InstructionOpcode["JV"] = 17] = "JV";
    InstructionOpcode[InstructionOpcode["JNV"] = 18] = "JNV";
    InstructionOpcode[InstructionOpcode["JS"] = 19] = "JS";
    InstructionOpcode[InstructionOpcode["JNS"] = 20] = "JNS";
    InstructionOpcode[InstructionOpcode["PUSH"] = 21] = "PUSH";
    InstructionOpcode[InstructionOpcode["PUSHADDR"] = 22] = "PUSHADDR";
    InstructionOpcode[InstructionOpcode["RET"] = 23] = "RET";
    InstructionOpcode[InstructionOpcode["POP"] = 24] = "POP";
    InstructionOpcode[InstructionOpcode["SETADDR"] = 25] = "SETADDR";
    InstructionOpcode[InstructionOpcode["DRAW"] = 26] = "DRAW";
    InstructionOpcode[InstructionOpcode["LOADKB"] = 27] = "LOADKB";
})(InstructionOpcode || (InstructionOpcode = {}));
;
class Emulator {
    constructor(drawHook) {
        this.ram = new Uint8Array(Emulator.RAMSIZE);
        this.rom = new Uint8Array(Emulator.ROMSIZE);
        this.regs = new Uint8Array(new Array(8).fill(0));
        this.cFlag = false;
        this.zFlag = false;
        this.vFlag = false;
        this.instReg = 0;
        this.addrReg = 0;
        this.prgcReg = 0;
        this.stckReg = 0;
        this.specReg = 0;
        this.instructionImplementations = [
            AND.bind(this),
            OR.bind(this),
            XOR.bind(this),
            ADD.bind(this),
            SUB.bind(this),
            MLT.bind(this),
            SHL.bind(this),
            SHR.bind(this),
            LOADRAM.bind(this),
            LOADROM.bind(this),
            STORE.bind(this),
            SET.bind(this),
            JMP.bind(this),
            JZ.bind(this),
            JNZ.bind(this),
            JC.bind(this),
            JNC.bind(this),
            JV.bind(this),
            JNV.bind(this),
            null,
            null,
            PUSH.bind(this),
            PUSHADDR.bind(this),
            RET.bind(this),
            POP.bind(this),
            SETADDR.bind(this),
            DRAW.bind(this)
        ];
        this.drawHook = drawHook;
    }
    getInstructionOpcode(instruction) {
        return getBits(instruction, 0, 5);
    }
    resetCpuState() {
        this.addrReg = 0;
        this.instReg = 0;
        this.prgcReg = 0;
        this.stckReg = 0;
        this.specReg = 0;
        this.cFlag = false;
        this.zFlag = false;
        this.vFlag = false;
        this.regs = new Uint8Array(new Array(8).fill(0));
        this.ram = new Uint8Array(new Array(Emulator.RAMSIZE).fill(0));
    }
    loadProgram(hexString) {
        let cnt = 0;
        for (let i = 0; i < hexString.length; i += 2) {
            let byte = parseInt(hexString[i] + hexString[i + 1], 16);
            this.rom[cnt++] = byte;
        }
        this.resetCpuState();
    }
    execInstr() {
        this.instReg = this.rom[this.prgcReg + 1] | (this.rom[this.prgcReg] << 8);
        this.prgcReg += 2;
        let opcode = this.getInstructionOpcode(this.instReg);
        this.instructionImplementations[opcode](this.instReg);
    }
    AluOpInstruction(instruction, op) {
        let [regA, regB] = extactParamsFromSourceDestTypeInstruction(instruction);
        this.cFlag = !!((op(this.regs[regA], this.regs[regB])) >> 8);
        this.regs[regA] = getBits(op(regA, regB), 0, 8);
        this.zFlag = this.regs[regA] == 0;
        this.vFlag = !!((this.regs[regA] * this.regs[regB]) >> 8);
    }
}
Emulator.kbBufferStartAddr = 45824;
Emulator.kbBufferEndAddr = 45567;
Emulator.stackMemStartAddr = 45568;
Emulator.stackMemEndAddr = 45823;
Emulator.vMemStartAddr = 45824;
Emulator.vMemEndAddr = 58111;
Emulator.RAMSIZE = 1 << 16;
Emulator.ROMSIZE = 1 << 16;
;
function getBits(number, begin, nBits) {
    return (number >> begin) & ((1 << nBits) - 1);
}
function extactParamsFromSourceDestTypeInstruction(instruction) {
    return [getBits(instruction, 5, 3), getBits(instruction, 8, 3)];
}
function AND(instruction) {
    this.AluOpInstruction(instruction, (a, b) => a & b);
}
function OR(instruction) {
    this.AluOpInstruction(instruction, (a, b) => a | b);
}
function XOR(instruction) {
    this.AluOpInstruction(instruction, (a, b) => a ^ b);
}
function ADD(instruction) {
    this.AluOpInstruction(instruction, (a, b) => a + b);
}
function SUB(instruction) {
    this.AluOpInstruction(instruction, (a, b) => a - b);
}
function MLT(instruction) {
    this.AluOpInstruction(instruction, (a, b) => a * b);
}
function SHL(instruction) {
    this.AluOpInstruction(instruction, (a, b) => a << b);
}
function SHR(instruction) {
    this.AluOpInstruction(instruction, (a, b) => a >> b);
}
function LOADRAM(instruction) {
    let reg = getBits(instruction, 5, 3);
    this.regs[reg] = this.ram[this.addrReg];
    this.zFlag = this.regs[reg] == 0;
}
function LOADROM(instruction) {
    let reg = getBits(instruction, 5, 3);
    this.regs[reg] = this.rom[this.addrReg];
    this.zFlag = this.regs[reg] == 0;
}
function STORE(instruction) {
    let reg = getBits(instruction, 5, 3);
    this.ram[this.addrReg] = this.regs[reg];
    this.zFlag = this.regs[reg] == 0;
}
function SET(instruction) {
    let reg = getBits(instruction, 5, 3);
    let data = getBits(instruction, 8, 8);
    this.regs[reg] = data;
    this.zFlag = this.regs[reg] == 0;
}
function JMP() {
    this.prgcReg = this.addrReg;
}
function JZ() {
    if (this.zFlag)
        this.prgcReg = this.addrReg;
}
function JNZ() {
    if (!this.zFlag)
        this.prgcReg = this.addrReg;
}
function JC() {
    if (this.cFlag)
        this.prgcReg = this.addrReg;
}
function JNC() {
    if (!this.cFlag)
        this.prgcReg = this.addrReg;
}
function JV() {
    if (this.vFlag)
        this.prgcReg = this.addrReg;
}
function JNV() {
    if (!this.vFlag)
        this.prgcReg = this.addrReg;
}
function PUSH() {
    this.ram[Emulator.stackMemStartAddr + this.stckReg] = getBits(this.prgcReg, 0, 8);
    this.ram[Emulator.stackMemStartAddr + this.stckReg + 1] = getBits(this.prgcReg, 8, 8);
    this.stckReg += 2;
}
function PUSHADDR() {
    this.ram[Emulator.stackMemStartAddr + this.stckReg] = getBits(this.addrReg, 0, 8);
    this.ram[Emulator.stackMemStartAddr + this.stckReg + 1] = getBits(this.addrReg, 8, 8);
    this.stckReg += 2;
}
function RET() {
    this.prgcReg = (this.ram[Emulator.stackMemStartAddr + this.stckReg - 1] << 8) | this.ram[Emulator.stackMemStartAddr + this.stckReg - 2];
    this.stckReg -= 2;
}
function POP() {
    this.addrReg = (this.ram[Emulator.stackMemStartAddr + this.stckReg - 1] << 8) | this.ram[Emulator.stackMemStartAddr + this.stckReg - 2];
    this.stckReg -= 2;
}
let first = 1;
function SETADDR(instruction) {
    let [regA, regB] = extactParamsFromSourceDestTypeInstruction(instruction);
    this.addrReg = (this.regs[regA] << 8) | this.regs[regB];
    this.cFlag = 0;
    this.zFlag = 0;
    this.vFlag = 0;
}
function DRAW() {
    this.drawHook(this.ram.slice(Emulator.vMemStartAddr, Emulator.vMemEndAddr));
}
//# sourceMappingURL=emulator.js.map