import functools 

class Rule:
    def __init__(self, whiteList, blackList):
      self.whiteList, self.blackList = set(whiteList), set(blackList)

    def match(self, c):
        if len(self.whiteList) > 0:
            return c in self.whiteList.difference(self.blackList)
        else:
            return not c in self.blackList and c != ''

    
    def __eq__(self, other):
        if len(self.whiteList) == 0 and len(self.blackList) == 0 :
            return self.blackList == other.blackList
        elif len(self.whiteList) == 0 or len(self.blackList) == 0:
            return False
        else:
            return self.whiteList.difference(self.blackList) == other.whiteList.difference(other.blackList)

class Transition:
    def __init__(self, rule, stateName):
        self.stateName = str(stateName)
        self.rule = rule
    
class State:
    alphabet = set()
    def __init__(self, name, accepting = False, starred = False, transitions = []):
        self.name = str(name)
        self.accepting = accepting == True
        self.starred = starred == True
        self.transitions = transitions.copy()
        for transition in self.transitions:
            self.alphabet.update(transition.rule.whiteList)

    def addTransition(self, transition):
        self.transitions.append(transition)
        self.alphabet.update(transition.rule.whiteList)
    
    def getNextStates(self, c):
        nextStates = set()
        for transition in self.transitions:
            if transition.rule.match(c):
                nextStates.add(transition.stateName)
        return nextStates

class FSM:
    states = {}
    precomputeUpToDate = False
    emptyClosureMem = {}
    alphabet = set()

    def __init__(self, states = None, startingStateName = '0'):
        if(states == None):
            states = []
        for state in states:
            self.addState(state)
        self.startingStateName = startingStateName

    def addState(self, state):
        self.precomputeUpToDate = False
        self.states[state.name] = state
        self.alphabet.update(state.alphabet)

    def getState(self, stateName):
        return self.states[stateName]

    def addTranition(self, stateName, transition):
        self.precomputeUpToDate = False
        self.states[stateName].addTranition(transition)
        self.alphabet.update(transition.rule.whiteList)
    
    def calcEmptyClosure(self):
        if self.precomputeUpToDate == True:
            return self.emptyClosureMem

        self.emptyClosureMem = {}

        def emptyClosureForState(state):
            self.emptyClosureMem[state] = set()
            for transition in self.states[state].transitions:
                if transition.rule.match(''):
                    self.emptyClosureMem[state].add(transition.stateName)
                    if not transition.stateName in self.emptyClosureMem:
                        emptyClosureForState(transition.stateName)                       
                    self.emptyClosureMem[state].update(self.emptyClosureMem[transition.stateName])
        
        for state in self.states:
            emptyClosureForState(state)

        self.precomputeUpToDate = True
        
    def move(self, T, c):
        self.calcEmptyClosure()
        def m(state, c):
            res = set()
            for transition in self.states[state].transitions:
                if(transition.rule.match(c)):
                    res.add(transition.stateName)
            return res
        return functools.reduce(lambda s, state: s.update(m(state, c) or {}) or s, T, set())

    def emptyClosure(self, T):
        self.calcEmptyClosure()
        return functools.reduce(lambda x, state: x.update(self.emptyClosureMem[state]) or x , T.copy(), T)

    def walk(self, s):          
        accepted = []
        currentStates = self.emptyClosure({self.startingStateName})
        i = 0
        while i <= len(s):
            newStates = self.emptyClosure(self.move(currentStates, s[i] if i < len(s) else 'eof'))
            currentStates = newStates
            for stateName in currentStates:
                state = self.states[stateName]
                if(state.accepting):
                    accepted.append((state.name, i))
            i+=1
        return accepted

   
def getDeterministicFSM(fsm):
    Dstates = dict()
    Dtran = dict()
    stack = set()
    newDstateId = 0

    def newDstate(T):
        nonlocal newDstateId
        stateName = str(newDstateId)
        Dstates[stateName] =  T
        Dtran[stateName] = {}
        stack.add(stateName)
        newDstateId+=1
        return stateName
    
    def getDstateName(T):
        for s in Dstates:
            if Dstates[s] == T:
                return s
        return None

    newDstate(fsm.emptyClosure({fsm.startingStateName}))
    while len(stack) > 0:
        T = stack.pop()
        for a in fsm.alphabet:
            if a == '':
                continue
            U = fsm.emptyClosure(fsm.move(Dstates[T], a))
            if(len(U) == 0):
                continue
            stateName = getDstateName(U)
            if(stateName == None):
                stateName = newDstate(U)
            
            Dtran[T].update({a: stateName})

    states = []
    for s in Dtran:
        transitions = []
        accepting = False
        for a in Dtran[s]:
            transitions.append(Transition( Rule([a], []), Dtran[s][a]))
        for ss in Dstates[s]:
            if(fsm.states[ss].accepting == True):
                accepting = True
                break
        states.append(State(s, accepting, False, transitions))
    return FSM(states, '0')
    
 
