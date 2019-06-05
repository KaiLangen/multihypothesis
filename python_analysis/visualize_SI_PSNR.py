import re
import numpy as np
import sys
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
            frameOrder.append(m.group(0))
        m = re.search('(?<=PSNR SI: )[0-9.]+', line)
        if m:
            psnr.append(m.group(0))

    return frameOrder, psnr


if __name__ == '__main__':       
    f1 = sys.argv[1]
    if not exists(f1):
        print("Error: invalid file name")
        sys.exit()
    order,psnr = parse_output(f1)

    plt.figure()
    plt.bar(order,psnr)
    plt.show()

