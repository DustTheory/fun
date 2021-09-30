const canvas = document.getElementById("display");
const ctx = canvas.getContext("2d");

let emulator = new Emulator((vMem) => {
    let cnt = 0;
    ctx.fillStyle = 'black';
    ctx.fillRect(0, 0, 32, 32);
    for(let i = 0; i < 32; i++){
        for(let j = 0; j < 32; j++){
            ctx.fillStyle = `rgb(${vMem[cnt]}, ${vMem[cnt+1]}, ${vMem[cnt+2]})`
            ctx.fillRect(j, i, 1, 1);
            cnt+=3;
        }
    }
});

document.getElementById("binaryInput").addEventListener('change', function(e) {
    let file = this.files[0];
    let reader = new FileReader();
    reader.readAsText(file);
    reader.onload = function (e) {
        var arrayBuffer = e.target.result;
        let hexString = e.target.result.split('\n')[1];
        if(!hexString)
            return;
        emulator.loadProgram(hexString);
        for(let i = 0; i < 20000; i++){
            emulator.execInstr();
        }
    }
});