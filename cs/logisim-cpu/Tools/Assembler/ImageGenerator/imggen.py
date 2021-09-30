from PIL import Image
import numpy

image = Image.open("lena.png")
image = image.resize((32,32),Image.ANTIALIAS)
data = numpy.asarray(image)

mem = numpy.zeros(65536)

imagestart = 45824

cnt = imagestart

for i in range(32):
    for j in range(32):
        for k in range(3):
            mem[imagestart] = data[i][j][k]
            imagestart += 1


""" f = open("ramimg2.bin", "w")

f.write("v3.0 hex bytes plain big-endian\n")
for x in mem:
    f.write(f"{int(x):02x}")
 """

f = open("lenaprogram.asm", "w")


def addrTo8bit(addr):
    return [int(bin(addr)[2:10], 2), int(bin(addr)[10:], 2)]

imagestart = 45824
for i in range(imagestart, imagestart+3072):
    [reg1, reg2] = addrTo8bit(i)
    f.write("SET 1 "+str(reg1)+"\n")
    f.write("SET 2 "+str(reg2)+"\n")
    f.write("SET 3 "+str(int(mem[i]))+"\n")
    f.write("SETADDR 1 2\n");
    f.write("STORE 3\n")
f.write("DRAW");


