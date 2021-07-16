from InitialRegexParser import Parser
p = Parser()
re = p.parse('abc[0-9]*d*.*xD')

print(re.walk("abc1337dddddlolxD"))
print(re.walk("abcxD"))
print(re.walk("abc0ddxD"))


