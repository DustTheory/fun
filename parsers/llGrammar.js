const { checkPropertyChange } = require("json-schema")

Array.prototype.subarray = function(i, j){
    var self = this, arr = [];
    for(var n = 0;i <= j; i++, n++){
        (function(i){
            Object.defineProperty(arr, n, {       //Array is an Object
                get: function(){
                    return self[i];
                },
                set: function(value){
                    self[i] = value;
                    return value;
                }
            });   
        })(i);
    }
    return arr;
}

class Grammar{
    constructor(startSymbol){
        this.startSymbol = startSymbol
        this.productions = {}
        this.FIRST = {}
        this.FOLLOW = {}
        this.sealed = false
    }

    addProduction(head, body){
        if(this.sealed === true)
            throw "Can't add production after grammar has been sealed"
        if (!this.productions[head])
            this.productions[head] = []
        this.productions[head].push(body)
    }

    seal(){
        if(this.sealed === true)
            return console.log("Grammar has already been sealed.")
        this.sealed = true 
    }

    throwErrorIfNotSealed(){
        if(this.sealed !== true)
            throw "Can't execucte operation before grammar is sealed."
    }

    firstProductionBody(body){
        this.throwErrorIfNotSealed()
        let canBeEmpty = true
        let first = new Set()
        for(let symbol of body){
            if(!this.productions[symbol]){ // symbol is a terminal
                first.add(symbol)
                canBeEmpty = false
                break
            }
            this.firstNonTerminal(symbol)
            this.FIRST[symbol].forEach(terminal => {if(terminal != '') first.add(terminal)})
            if(!this.FIRST[symbol].has('')){
                canBeEmpty = false
                break
            }
        }
        if(canBeEmpty)
            first.add('')
        return first
    }

    firstNonTerminal(nonTerminal){
        this.throwErrorIfNotSealed()
        if(this.FIRST[nonTerminal])
            return
        
        this.FIRST[nonTerminal] = new Set()
        for(let productionBody of this.productions[nonTerminal])
            this.firstProductionBody(productionBody).forEach(terminal => this.FIRST[nonTerminal].add(terminal))
    }

    followNonTerminal(nonTerminal){
        this.throwErrorIfNotSealed()
        if(this.FOLLOW[nonTerminal])
            return 
        this.FOLLOW[nonTerminal] = new Set()
        if(nonTerminal == this.startSymbol)
            this.FOLLOW[nonTerminal].add('$')
        for(let productionHead in this.productions){
            for(let productionBody of this.productions[productionHead]){
                let indices = productionBody.reduce(function(a, e, i) {
                    if (e === nonTerminal)
                        a.push(i+1);
                    return a;
                }, []);  
                for(let startIndex of indices){
                    let first = this.firstProductionBody(productionBody.slice(startIndex))
                    first.forEach(terminal => { if(terminal != '') this.FOLLOW[nonTerminal].add(terminal) })
                    if(first.has('') || startIndex == productionBody.length){
                        this.followNonTerminal(productionHead)
                        this.FOLLOW[productionHead].forEach(terminal => this.FOLLOW[nonTerminal].add(terminal))   
                    }
                }        
            }    
        }
    }

}

g = new Grammar("E")
g.addProduction("E", ["T", "E'"])
g.addProduction("E'", ["+", "T", "E'"])
g.addProduction("E'", [''])
g.addProduction("T", ["F", "T'"])
g.addProduction("T'", ["*", "F", "T'"])
g.addProduction("T'", [''])
g.addProduction("F", ["(", "E", ")"]);
g.addProduction("F", ["id"]);

g.seal()

for (let production in g.productions){
    g.firstNonTerminal(production)
    g.followNonTerminal(production)
}

console.log("FIRST:", g.FIRST)
console.log("FOLLOW:", g.FOLLOW)
