class Lexer:
    def __init__(self, inputStr):
        self.input = inputStr
        self.pos = 0
        self.getNextToken()
        
    def getNextToken(self):
        if self.pos >= len(self.input):
            self.lookahead = 'eof'
            self.attribute = 'eof'
            return self.lookahead

        while self.pos < len(self.input):
            c = self.input[self.pos]
            self.pos+=1
            if c == '\\':
                if self.pos >= len(self.input):
                    self.syntaxError('Escape character reached eof')

                nc = self.input[self.pos]  
                if nc not in ['.', '(', ')', '+', '?', '^', '[', ']', '-', 'n', '*', '\\']:
                    self.syntaxError('Invalid escaped character '+nc)
                if nc == 'n':
                    self.attribute = '\n'
                elif nc == '\\':
                    self.attribute = '\\'
                else:
                    self.attribute = nc
                self.lookahead = 'C'
                self.pos+=1
                break
            else:
                if c in ['.', '(', ')', '+', '?', '^', '[', ']', '-', '*', '|']:
                    self.attribute = self.lookahead = c
                    break
                else:
                    self.attribute = c
                    self.lookahead = 'C'
                    break
        return self.lookahead
            
    def syntaxError(self, errMsg):
        raise Exception('Sytax error: ' + errMsg + 'at index {}'.format(self.pos))
