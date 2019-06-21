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
    frameKB =[]
    psnrSI =[]
    psnrWZ =[]
    for line in lines:
        m = re.search('(?<=Decoding frame )\d+', line)
        if m:
            frameOrder.append(int(m.group(0)))
        m = re.search('(?<=PSNR SI: )[0-9.]+', line)
        if m:
            psnrSI.append(float(m.group(0)))
        m = re.search('(?<=\(Y frame\): )[0-9.]+', line)
        if m:
            frameKB.append(float(m.group(0)))
        m = re.search('(?<=PSNR WZ: )[0-9.]+', line)
        if m:
            psnrWZ.append(float(m.group(0)))

    nFrames = len(frameOrder)
    print(nFrames, len(psnrSI), len(psnrWZ))
    psnrWZ = [psnrWZ[i] - psnrSI[i] for i in range(nFrames)]
    return zip(frameOrder, psnrSI, psnrWZ, frameKB)

def parse_keystats(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        m = re.search('(?<=SNR Y\(dB\))[ |]+[0-9.]+', line)
        if m:
            avgPSNRY = float(m.group(0).lstrip(' |'))
        m = re.search('(?<=average bits/frame)[ |]+[0-9.]+', line)
        if m:
            avgBits = float(m.group(0).lstrip(' |'))
        m = re.search('(?<=Coeffs. C)[ |]+[0-9.]+', line)
        if m:
            avgChromaBits = float(m.group(0).lstrip(' |'))

    return (avgBits - avgChromaBits)/(8000), avgPSNRY 

###############################################################################

if __name__ == '__main__':       
    if len(sys.argv) < 5:
        print("Usage: python visualize_SI_PSNR.py dataFile gop type statFile")
        sys.exit()
    f1 = sys.argv[1]
    gop = int(sys.argv[2]);
    name = sys.argv[3]
    statFile = sys.argv[4]
    if not exists(f1):
        print("Error: invalid output file: {}".format(f1))
        sys.exit()
    if not exists(statFile):
        print("Error: invalid JM stats file: {}".format(statFile))
        sys.exit()
    colorCycle = ['#58508d', '#bc5090', '#f24545', '#003f5c', '#3a3a3c']
    proposed_labels = ["Chroma-ME", "MCI"]

    # parse data from file
    tup = parse_output(f1)
    avgKeyRate, avgKeyPsnr = parse_keystats(statFile)

    # sort by frame order and insert dummy vars for key frames
    tup = sorted(tup)
    cntr = 0
    idx = 0
    while idx < len(tup):
        if cntr < tup[idx][0]:
            tup.insert(idx, (cntr, 0., 0., avgKeyRate))
        idx += 1
        cntr += 1
    nFrames = len(tup)
    tup.insert(idx, (cntr, 0., 0., avgKeyRate))
    order, psnrSI, psnrWZ, nkbytes = zip(*tup)
    print("Total Rate: {} kB".format(sum(nkbytes)))

    # init figure
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Frame #')
    ax1.set_ylabel('PSNR')
    plt.title("{} Codec: GOP={}".format(name,gop))

    # get key-frame stats
    keyNo, keyPsnr = zip(*[(i, avgKeyPsnr) for i in range(0, nFrames+1, gop)])
    nFrames += len(keyNo)
    plt.bar(keyNo, keyPsnr, color=colorCycle[3],
            label='Key Frames', align='center')

    # separate lists into even and odd wz frames
    for i in range(0,2):
        avg = sum(psnrSI[i::2])/len(psnrSI[i::2])
        print("{}, Avg PSNR={}".format(proposed_labels[i], avg))
        ax1.bar(order[i::2], psnrSI[i::2], color=colorCycle[i],
                align='center', label="{}".format(proposed_labels[i]))
        if i == 0:
            ax1.bar(order[i::2], psnrWZ[i::2], color=colorCycle[2],
                    align='center', bottom=np.array(psnrSI[i::2]))
        else: 
            ax1.bar(order[i::2], psnrWZ[i::2], color=colorCycle[2],
                    align='center', label='WZ reconstruction',
                    bottom=np.array(psnrSI[i::2]))

    order = np.array(order, dtype=int)
    ax2 = ax1.twinx()
    ax2.set_ylabel('kB / frame')
    ax2.plot(order, nkbytes, color=colorCycle[4], marker='s')

    ax1.legend(loc='best')
#    ax1.set_ylim(20,50)
    ax1.set_xlim(-1,nFrames+1)
    ax1.xaxis.set_ticks(np.arange(0, nFrames, 2)) 
    fig.tight_layout()
    plt.savefig("{}_si_gop{}.png".format(name,gop))
    plt.show()
