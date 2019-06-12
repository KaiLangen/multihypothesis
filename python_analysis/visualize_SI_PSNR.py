import re
import numpy as np
import sys
import math
from os.path import join, exists
import matplotlib.pyplot as plt

def parse_output(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    frameOrder =[]
    psnr =[]

    for line in lines:
        m = re.search('(?<=Decoding frame )\d+', line)
        if m:
            frameOrder.append(int(m.group(0)))
        m = re.search('(?<=PSNR SI: )[0-9.]+', line)
        if m:
            psnr.append(float(m.group(0)))

    return zip(frameOrder, psnr)

if __name__ == '__main__':       
    f1 = sys.argv[1]
    gop = int(sys.argv[2]);
    isBidirectional = bool(int(sys.argv[3]))
    name = sys.argv[4]
    if not exists(f1):
        print("Error: invalid file name")
        sys.exit()
    colorCycle = ['r','c','b']

    currL = []
    dist = []
    tups = parse_output(f1)
    nFrames = 50
    if isBidirectional:
        for n in range(int(math.log(gop,2)-1),-1,-1):
            dist.append(1<<n)
            currL.append(list(filter(lambda t: t[0]%(1<<n) == 0, tups)))
            tups = list(filter(lambda t: t[0]%(1<<n) != 0, tups))
    else:
        currL.append(list(filter(lambda t: t[0]%2 == 1, tups)))
        currL.append(list(filter(lambda t: t[0]%2 == 0, tups)))
    plt.figure()
    plt.ylim(20,40)
    plt.xlim(0,nFrames+7)

    proposed_labels = ["MCI", "Chroma-ME"]
    for i,c in enumerate(currL):
        order, psnr = zip(*c)
        avg = sum(psnr)/len(psnr)
        x = range(1,nFrames+8)
        y = [avg for j in x]
        if isBidirectional:
            idx=len(colorCycle)-int(math.log(gop,2))+i
            print("Frame Distance={}, Avg PSNR={}".format(order[0], avg))
            plt.bar(order, psnr, color=colorCycle[idx],
                    label="Frame distance={}".format(dist[i]))
            plt.plot(x,y, color=colorCycle[idx])
        else:
            print("{}, Avg PSNR={}".format(proposed_labels[i], avg))
            plt.bar(order, psnr, color=colorCycle[-i-1],
                    label="{}".format(proposed_labels[i]))
            plt.plot(x,y, color=colorCycle[-i-1])

        plt.text(nFrames-1, avg+.1, 'Avg PSNR\n={:.2f}'.format(avg))


    plt.legend(loc='best')
    plt.xlabel('Frame #')
    plt.ylabel('PSNR')
    if isBidirectional:
        plt.title("Original DISCOVER decoding: GOP={}".format(gop))
    else:
        plt.title("Proposed Codec: GOP={}".format(gop))
    plt.savefig("{}_si_gop{}.png".format(name,gop))

#    plt.show()

