
function areSameItemSet(I1, I2) {
    if(I1.length != I2.length)
        return false;
    let cmp = (a, b) => a.productionId == b.productionId ? a.nextElementIndex - b.nextElementIndex : a.productionId - b.productionId    
    I1 = I1.sort(cmp)
    I2 = I2.sort(cmp)
    for(let i = 0; i < I1.length; i++)
        if(I1[i].productionId != I2[i].productionId || I1[i].nextElementIndex != I2[i].nextElementIndex)
            return false
    return true;
}

class Grammar{
    constructor(startSymbol){
        this.startSymbol = startSymbol
        this.productions = {}
        this.productionIdMap = {}
        this.newProductionId = 0
        this.FIRST = {}
        this.FOLLOW = {}
        this.GOTO = {}
        this.sealed = false
        this.grammarSymbols = new Set()
    }

    addProduction(head, body){
        if(this.sealed === true)
            throw "Can't add production after grammar has been sealed"
        if (!this.productions[head])
            this.productions[head] = []
        let thisProductionId = this.newProductionId++
        let newProduction = { productionId: thisProductionId, head: head,  body: body }
        this.productions[head].push(newProduction)
        this.productionIdMap[thisProductionId] = newProduction
        body.forEach(symbol => this.grammarSymbols.add(symbol))
        return newProduction.productionId
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
        for(let productionBody of this.productions[nonTerminal]){
            this.firstProductionBody(productionBody.body).forEach(terminal => this.FIRST[nonTerminal].add(terminal))
        }
    }

    followNonTerminal(nonTerminal){
        this.throwErrorIfNotSealed()
        if(this.FOLLOW[nonTerminal])
            return 
        this.FOLLOW[nonTerminal] = new Set()
        if(nonTerminal == this.startSymbol)
            this.FOLLOW[nonTerminal].add('$')
        for(let productionHead in this.productions){
            for(let productionObj of this.productions[productionHead]){
                let productionBody = productionObj.body;
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

    closure(I){
        let closure = []
        let itemsAlreadyAdded = {}
        for(let item of I){
            if(itemsAlreadyAdded?.[item.productionId]?.[item.nextElementIndex])
                continue
            closure.push(item)
            if(!itemsAlreadyAdded[item.productionId])
                itemsAlreadyAdded[item.productionId] = {}
            itemsAlreadyAdded[item.productionId][item.nextElementIndex] = true
        }
        let addedItemToClosure
        do{
            addedItemToClosure = false
            for (let item of closure){
                let nextProductionElement = this.productionIdMap[item.productionId].body[item.nextElementIndex]
                if(!this.productions[nextProductionElement])
                    continue
                for(let production of this.productions[nextProductionElement]){
                    if(itemsAlreadyAdded?.[production.productionId]?.[0])
                        continue
                    closure.push(new Item(production.productionId, 0))
                    if(!itemsAlreadyAdded[production.productionId])
                        itemsAlreadyAdded[production.productionId] = {}
                    itemsAlreadyAdded[production.productionId][0] = true
                    addedItemToClosure = true
                }
            }
        }while(addedItemToClosure)
        return closure
    }

    goto(I, X){
        let initialClosureElements = []
        for(let item of I){
            let production = this.productionIdMap[item.productionId].body
            if(production[item.nextElementIndex] == X)
                initialClosureElements.push(new Item(item.productionId, item.nextElementIndex+1))
        }
        return this.closure(initialClosureElements)
    }
    items(){
        let nextCollectionItemId = 0
        let collection = [{id: nextCollectionItemId++, set: this.closure([new Item( 0, 0 )])}]
        let addedItemToCollection
        do{
            addedItemToCollection = false
            for(let I of collection){
                this.GOTO[I.id] = this.GOTO[I.id] || {}
                for(let X of this.grammarSymbols){
                    let goto = this.goto(I.set, X)
                    if(goto.length == 0)
                        continue
                    let alreadyExisting = collection.find(tmpI => areSameItemSet(tmpI.set, goto))
                    
                    if(alreadyExisting){
                        this.GOTO[I.id][X] = alreadyExisting.id
                        alreadyExisting.in.push({from: I.id, by: X})  
                        continue
                    }
                    let newId = nextCollectionItemId++
                    this.GOTO[I.id][X] = newId
                    collection.push({id: newId, in: [{from: I.id, by: X}], set: goto})
                }
            }
        }while(collection.length < 10)
        return collection
    }
}

class Item {
    constructor(productionId, nextElementIndex){
        this.productionId = productionId
        this.nextElementIndex = nextElementIndex
    }
}

g = new Grammar("E'")
let id = g.addProduction("E'", ["E"])
g.addProduction("E", ["E", "+", "T"])
g.addProduction("E", ["T"])
g.addProduction("T", ["T", "*", "F"])
g.addProduction("T", ["F"])
g.addProduction("F", ["(", "E", ")"])
g.addProduction("F", ["id"])
g.seal()

function printItemSetItems(itemSet) {
    itemSet.forEach( item => {
        let production = g.productionIdMap[item.productionId]
        console.log(production.head, '->', production.body.slice(0, item.nextElementIndex).join('')+ '.'+ production.body.slice(item.nextElementIndex).join(''))
    })
}
let C = g.items()
for (let production in g.productions){
    g.firstNonTerminal(production)
    g.followNonTerminal(production)
}
let ACTION = generateActionTable()

function generateActionTable(){
    let ACTION = {}
    for( let setName in C ){
        let itemSet = C[setName].set
        let setId = C[setName].id
        ACTION[setId] = {}
        for (item of itemSet) {
            production = g.productionIdMap[item.productionId]
            if(item.nextElementIndex >= production.body.length){
                if(production.head == g.startSymbol)
                    ACTION[setId]['$'] = 'ac'
                else
                    for(let x of g.FOLLOW[production.head])
                        ACTION[setId][x] = 'r'+item.productionId
            }else{
                let c = production.body[item.nextElementIndex]
                if(!g.productions[c]){ // terminal
                    ACTION[setId][c] = 's' + g.GOTO[setId][c]
                }
            }
        }
    }
    return ACTION;
}

console.log(ACTION)

function parse(terminals){
    terminals.push("$")
    let stack = [0]
    let nextTerminalIndex = 0;
    while(true){
        let currId = stack[stack.length-1]
        let curr   = C.find(itemSet => itemSet.id == currId)
        let nextTerminal = terminals[nextTerminalIndex] || 'eof'
        let action = ACTION[currId][nextTerminal] || 'error';
        if(action[0] == 'r'){ // reduce
            let production = g.productionIdMap[action.substr(1)]
            stack = stack.slice(0, stack.length-production.body.length)
            stack.push(g.GOTO[stack.slice(-1)[0]][production.head])
            console.log(production.head, '->', production.body.join(' '))
        }else if(action[0] == 's'){ // shift
            stack.push(action.substr(1))
            nextTerminalIndex++
        }else if(action == 'ac'){ // accept
            console.log("ACCEPT")
            break;
        }else{ // error
            console.log(action)
            throw "ERROR"
        }
    }
}

parse(['id', '+', '(', 'id', '+', 'id' , '*', 'id', ')'])