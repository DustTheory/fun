from InitialRegexLexer import Lexer
from FSM import FSM, State, Transition, Rule, getDeterministicFSM

class Parser:
    def parse(self, reStr):
        states = [ ]
        newStateId = 1
        emptyTransitionRule = Rule([''], [])

        def makeNewState():
            nonlocal newStateId
            newState = State(newStateId)
            newStateId+=1
            states.append(newState)
            return newState

        def match(token):
            oldAttribute = lexer.attribute
            if token != lexer.lookahead:
                raise Exception("Syntax error")
            lexer.getNextToken()
            return oldAttribute

        def applyModifier(g, modifier):
            if modifier == '?':
                g[0].addTransition(Transition(emptyTransitionRule, g[1].name))
                return g
            elif modifier == '*':
                g[0].addTransition(Transition(emptyTransitionRule, g[1].name))
                g[1].addTransition(Transition(emptyTransitionRule, g[0].name))
                return g
            elif modifier == '+':
                newState = makeNewState()        
                g[1].addTransition(Transition(emptyTransitionRule, newState.name))
                newState.addTransition(Transition(emptyTransitionRule, g[0].name))
                return g
            return g

        def concatenate(g1, g2):
            if g2 == None:
                return g1
            g1[1].addTransition(Transition(emptyTransitionRule, g2[0].name))
            return g1[0], g2[1]

        def disjunction(g1, g2):
            if(g2 == None):
                return g1
            newState1 = makeNewState()
            newState2 = makeNewState()
            newState1.addTransition(Transition(emptyTransitionRule, g1[0].name))
            newState1.addTransition(Transition(emptyTransitionRule, g2[0].name))
            g1[1].addTransition(Transition(emptyTransitionRule, newState2.name))
            g2[1].addTransition(Transition(emptyTransitionRule, newState2.name))
            return newState1, newState2
        
        def createG(rule):
            newState1 = makeNewState()
            newState2 = makeNewState()
            newState1.addTransition(Transition(rule, newState2.name))
            return newState1, newState2

        def createRange(rangeStart, rangeEnd):
            if rangeEnd == None:
                return createG(Rule([rangeStart],[]))
            return createG(Rule([chr(c) for c in range(ord(rangeStart), ord(rangeEnd)+1)], []))

        def character():
            return match('C')

        def rangeRest():
            if lexer.lookahead == '-':
                match('-')
                return character()
            return None

        def rangelistRest():
            if lexer.lookahead != 'C':
                return None
            return disjunction(_range(), rangelistRest())

        def _range():
            return createRange(character(), rangeRest())

        def rangeList():
            return disjunction(_range(), rangelistRest())

        def single():
            if lexer.lookahead == '.':
                match('.')
                return createG(Rule([],[]))
            elif lexer.lookahead == '(':
                match('(')
                r = disjunctions()
                match(')')
                return r
            elif lexer.lookahead == '[':
                match('[')
                r = rangeList()
                match(']')
                return r
            else:
                return createG(Rule([character()],[]))

        def modifiedSingleModifier():
            if lexer.lookahead in ['?', '+', '*']:
                return match(lexer.lookahead)
            else:
                return ''

        def modifiedSingle():
            return applyModifier(single(), modifiedSingleModifier())         

        def concatenations():
            return concatenate(modifiedSingle(), concatenationsRest())

        def concatenationsRest():
            if not lexer.lookahead in ['C', '.', '(', '[']:
                return None
            return concatenate(modifiedSingle(), concatenationsRest())

        def disjunctions():
            return disjunction(concatenations(), disjunctionsRest())

        def disjunctionsRest():
            if(lexer.lookahead != '|'):
                return None
            match('|')
            return disjunction(concatenations(), disjunctionsRest())

        lexer = Lexer(reStr)
        
        startingState, finalState = disjunctions()
        finalState.accepting = True
        return FSM(states, startingStateName=startingState.name)
