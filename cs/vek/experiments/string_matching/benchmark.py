#Stolen from https://github.com/WojciechMula/sse4-strstr/blob/master/original/sse4_strstr-test.py

import sys, os, random

filename = "<unspecified>"

try:
	filename = sys.argv[1]
	string = open(filename, "r").read()
except:
    print("can't open '%s'" % filename)
    sys.exit(1)

try:
	random.seed(int(sys.argv[3]))
except:
	pass

def time_command(command):
	os.system('/usr/bin/time -o /tmp/measure -f "%U" ' + command)
	f = open("/tmp/measure", "r")
	t = float(f.read())
	f.close()
	return t

def time(command1, command2, command3, iters=10):
	while iters < 1000000000000000000000:
		t1 = time_command(command1.replace("__iters__", str(iters)))
		if t1 > 1:
			t2 = time_command(command2.replace("__iters__", str(iters)))
			t3 = time_command(command3.replace("__iters__", str(iters)))
			return iters, t1, t2, t3
		else:
			iters *= 10
	return 0, 0, 0, 0


def compare(filename, wordpos, word, wordlen):
	word = word.replace("%", "%%")
	cmd1 = './0x80plAlgo1 "%s" libc __iters__ "%s" > /dev/null' % (filename, word)
	cmd2 = './0x80plAlgo1 "%s" impl2 __iters__ "%s" > /dev/null' % (filename, word)
	cmd3 = './0x80plAlgo1 "%s" impl1 __iters__ "%s" > /dev/null' % (filename, word)
	# cmd4 = './0x80plAlgo1 "%s" impl2 __iters__ "%s" > /dev/null' % (filename, word)
	_, t1, t2, t3 = time(cmd1, cmd2, cmd3)
	return "[%d,%d] libc=%0.3fs naive=%0.3fs impl1=%0.2f" % (wordpos, wordlen, t1, t2, t3)


logname   = "0x80plAlgo1.log"
lognumber = 1
while True:
	if not os.path.exists(logname):
		log = open(logname, "w")
		break
	else:
		logname = "0x80plAlgo1_%d.log" % lognumber
		lognumber += 1

try:
	for n in range(4, 64):
		i1 = random.randint(   0, 64)
		i2 = random.randint(  65, 1024)
		i3 = random.randint(1024, len(string)-n)
		print ("length", n)
		for i in [i1, i2, i3]:
			word = string[i:i+n]
			for c in "\\`()<>{}\"":
				word = word.replace(c, "\\" + c)
			
			cmd = './0x80plAlgo1 "%s" verify 1 "%s"' % (filename, word)
			err = os.system(cmd)
			if err:
				print (repr(string[i:i+n]))
				sys.exit(1)
			else:
				s = compare(filename, i, word, n)
				log.write(s + "\n")
				print (s)
except:
	import traceback
	traceback.print_exc()
	log.close()