enum InstructionOpcode {
    AND = 0,
    OR,
    XOR,
    ADD,
    SUB,
    MLT,
    SHL,
    SHR,
    LOADRAM, // Also called just LOAD,
    LOADROM,
    STORE,
    SET,
    JMP,
    JC,
    JNC,
    JZ,
    JNZ,
    JV,
    JNV,
    JS, // JS and JNS are reserved but will probably never be implemented in hardware
    JNS, 
    PUSH,
    PUSHADDR,
    RET,
    POP,
    SETADDR,
    DRAW,
    LOADKB
};

class Emulator {

    // Memory map
    // First 45536 of the 2^16 bytes are free unreserved memory
    // Rest is divided like so:
    // 45536-45567 -> keyboard buffer
    // 45568-45823 -> stack memory
    // 45824-58111 -> video framebuffer
    // 58112-65535 -> unused for now
    static readonly kbBufferStartAddr = 45824;
    static readonly kbBufferEndAddr = 45567;
    static readonly stackMemStartAddr = 45568;
    static readonly stackMemEndAddr = 45823;
    static readonly vMemStartAddr = 45824;
    static readonly vMemEndAddr = 58111;

    static readonly RAMSIZE = 1 << 16;
    static readonly ROMSIZE = 1 << 16;

    // RAM and ROM memory
    ram: Uint8Array = new Uint8Array(Emulator.RAMSIZE); // 2^16 bytes ram memory
    rom: Uint8Array = new Uint8Array(Emulator.ROMSIZE); // 2^16 bytes rom memory

    // CPU registers
    regs: Uint8Array = new Uint8Array(new Array<number>(8).fill(0)); // general purpose registers

    // c, z and v CPU flag registers
    cFlag : boolean = false;
    zFlag : boolean = false;
    vFlag : boolean = false;

    // Hidden registers
    // Hidden registers are used internally by the cpu but can't be referenced
    // in machine code. In hardware, some are implemented by using 2 seperate 
    // 8 bit registers and some microcode logic. This will not be emulated here,
    // and instead those registers will be stored in a regular ts number and will
    // be treated as single 16 bit registers(values).

    instReg: number = 0; // 16 bit instruction register
    addrReg: number = 0; // 16 bit address register
    prgcReg: number = 0; // 16 bit program counter register
    stckReg: number = 0; // 8 bit stack pointer register
    specReg: number = 0; // 8 bit special register (used only by microcode)

    instructionImplementations: Array<(instruction: number) => void> = [
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
        JC.bind(this),
        JNC.bind(this),
        JZ.bind(this),
        JNZ.bind(this),
        JV.bind(this),
        JNV.bind(this),
        null, // JS 
        null, // JNS
        PUSH.bind(this),
        PUSHADDR.bind(this),
        RET.bind(this),
        POP.bind(this),
        SETADDR.bind(this),
        DRAW.bind(this)
        ]

    drawHook : (vMem :Uint8Array) => void;

    /**
     * Retrieves 5 bit instruction opcode from 16 bit instruction
     * @param {number} instruction 16 bit instruction as a number
     * @returns {InstructionOpcode} 5 bit instruction opcode
     */
    getInstructionOpcode(instruction: number): InstructionOpcode {
        return getBits(instruction, 0, 5);
    }

    /**
     * Reset CPU state to defaults - ready to execute new program.
     */
    resetCpuState(): void {
        // Set all registers and entire ram memory to 0
        this.addrReg = 0; 
        this.instReg = 0;
        this.prgcReg = 0;
        this.stckReg = 0;
        this.specReg = 0;
        this.cFlag = false;
        this.zFlag = false;
        this.vFlag = false;
        this.regs = new Uint8Array(new Array<number>(8).fill(0)); 
        this.ram = new Uint8Array(new Array<number>(Emulator.RAMSIZE).fill(0));
    }

    constructor(drawHook: (vMem :Uint8Array) => void){
        this.drawHook = drawHook;
    }

    /**
     * Loads program into ROM memory, and resets CPU state.
     * @param {Uint8Array} program Raw program data to be loaded into ROM memory
     */
    loadProgram(hexString: string): void {
        let cnt: number = 0;
        for(let i = 0; i < hexString.length; i+=2){
            let byte = parseInt(hexString[i] + hexString[i+1], 16);
            this.rom[cnt++] = byte;
        }

      //  this.rom.set(program.slice(Emulator.ROMSIZE)) // copy program into rom memory
        this.resetCpuState();
    }

    execInstr() {
        // Load instruction from rom by combining 2 bytes at program counter into a 16 bit instruction
        this.instReg = this.rom[this.prgcReg+1] | (this.rom[this.prgcReg] << 8);
        this.prgcReg += 2;
        // Find instruction instruction op implementation by opcode and call it.
        let opcode = this.getInstructionOpcode(this.instReg);
        this.instructionImplementations[opcode](this.instReg);
    }

    AluOpInstruction(instruction: number, op: (regA:  number, regB: number) => number) : void {
        let [regA, regB] = extactParamsFromSourceDestTypeInstruction(instruction);
        this.cFlag = !!((op(this.regs[regA], this.regs[regB])) >> 8);
        this.regs[regA] = getBits(op(this.regs[regA], this.regs[regB]), 0, 8);
        this.zFlag = this.regs[regA] == 0;
        this.vFlag = !!((this.regs[regA] * this.regs[regB]) >> 8);
    }
};

// CPU OPS implementations

/**
 * Extracts part of a binary representation of a number as a decimal value
 * @param {number} number
 * @param {number} begin Index of first bit to extract. (Starting from lsb)
 * @param {number} nBits Number of bits to include starting from begin.
 * @returns {number} Decimal representation of a portion of a binary representation of a number
 */
function getBits(number: number, begin: number, nBits: number): number {
    return (number >> begin) & ((1 << nBits) - 1);
}

/**
 * Returns source and destination register indicies from an instruction of form: opcode register register 
 * @param {number} instruction Instruction
 * @returns {[number, numbre]} Register indicies (0 - 7 inclusive) extracted from instruction
 */
function extactParamsFromSourceDestTypeInstruction(instruction: number): [number, number] {
    return [getBits(instruction, 5, 3), getBits(instruction, 8, 3)];
}
 
/**
 * AND instruction implementation
 * @param instruction 
 */
function AND(instruction: number): void {
    this.AluOpInstruction(instruction, (a: number, b: number) => a & b);
}

/**
 * OR instruction implementation
 * @param instruction 
 */
function OR(instruction: number): void {
    this.AluOpInstruction(instruction, (a: number, b: number)=> a | b);
}

/**
 * XOR instruction implementation
 * @param instruction 
 */
function XOR(instruction: number): void {
    this.AluOpInstruction(instruction, (a: number, b: number)=> a ^ b);
}

/**
 * ADD instruction implementation
 * @param instruction 
 */
function ADD(instruction: number): void {
    this.AluOpInstruction(instruction, (a: number, b: number)=> a + b);
}

/**
 * SUB instruction implementation
 * @param instruction 
 */
function SUB(instruction: number): void {
    this.AluOpInstruction(instruction, (a: number, b: number)=> a - b);
}

/**
 * MLT instruction implementation
 * @param instruction 
 */
function MLT(instruction: number): void {
    this.AluOpInstruction(instruction, (a: number, b: number)=> a * b);
}

/**
 * SHL instruction implementation
 * @param instruction 
 */
function SHL(instruction: number): void {
    this.AluOpInstruction(instruction, (a: number, b: number)=> a << b);
}

/**
 * SHR instruction implementation
 * @param instruction 
 */
function SHR(instruction: number): void {
    this.AluOpInstruction(instruction, (a: number, b: number)=> a >> b);
}

function LOADRAM(instruction: number): void {
    let reg = getBits(instruction, 5, 3);
    this.regs[reg] = this.ram[this.addrReg];
    this.zFlag = this.regs[reg] == 0;
}

function LOADROM(instruction: number): void {
    let reg = getBits(instruction, 5, 3);
    this.regs[reg] = this.rom[this.addrReg];
    this.zFlag = this.regs[reg] == 0;
}

function STORE(instruction: number): void {
    let reg = getBits(instruction, 5, 3);
    this.ram[this.addrReg] = this.regs[reg];
    this.zFlag = this.regs[reg] == 0;
}

function SET(instruction: number): void {
    let reg = getBits(instruction, 5, 3);
    let data = getBits(instruction, 8, 8);
    this.regs[reg] = data;
    this.zFlag = this.regs[reg] == 0;
}

function JMP(): void {
    this.prgcReg = this.addrReg;
}

function JZ(): void {
    if(this.zFlag)
        this.prgcReg = this.addrReg;
}

function JNZ(): void {
    if(!this.zFlag)
        this.prgcReg = this.addrReg;
}

function JC(): void {
    if(this.cFlag)
        this.prgcReg = this.addrReg;
}

function JNC(): void {
    if(!this.cFlag)
        this.prgcReg = this.addrReg;
}

function JV(): void {
    if(this.vFlag)
        this.prgcReg = this.addrReg;
}

function JNV(): void {
    if(!this.vFlag)
        this.prgcReg = this.addrReg;
}

function PUSH(): void {
    this.ram[Emulator.stackMemStartAddr + this.stckReg] = getBits(this.prgcReg, 0, 8);
    this.ram[Emulator.stackMemStartAddr + this.stckReg + 1] = getBits(this.prgcReg, 8, 8);
    this.stckReg += 2;
}

function PUSHADDR(): void {
    this.ram[Emulator.stackMemStartAddr + this.stckReg] = getBits(this.addrReg, 0, 8);
    this.ram[Emulator.stackMemStartAddr + this.stckReg + 1] = getBits(this.addrReg, 8, 8);
    this.stckReg += 2;
}

function RET(): void {
    this.prgcReg = (this.ram[Emulator.stackMemStartAddr + this.stckReg-1] << 8) |  this.ram[Emulator.stackMemStartAddr + this.stckReg-2];
    this.stckReg -= 2;
}

function POP(): void {
    this.addrReg = (this.ram[Emulator.stackMemStartAddr + this.stckReg-1] << 8) |  this.ram[Emulator.stackMemStartAddr + this.stckReg-2];
    this.stckReg -= 2;
}

let first = 1;

function SETADDR(instruction: number): void {
    let [regA, regB] = extactParamsFromSourceDestTypeInstruction(instruction);
    this.addrReg = (this.regs[regA] << 8) | this.regs[regB];
    this.cFlag = 0;
    this.zFlag = 0;
    this.vFlag = 0;
}

function DRAW() {
    this.drawHook(this.ram.slice(Emulator.vMemStartAddr, Emulator.vMemEndAddr));
}

