import matplotlib.pyplot as plt
import numpy as np
import math
with open("coeffs.txt") as f:
    it = 100

    total = 0
    totalR = 0
    for a in range(it):
        old = f.readline().split(",")
        new = f.readline().split(",")
        #print(new)
        old.pop(len(old) - 1)
        new.pop(len(new) - 1)

        for i in range(len(old)):
            #print(i)
            old[i] = abs(int(old[i]))
            new[i] = abs(int(new[i]))
        #print(old)
        maxOld = max(old)
        maxNew = max(new)
        
        ratio = maxNew // maxOld
        totalR += math.log(ratio, 2)


        GCD = np.gcd.reduce(new)
        fancyGCD = math.gcd(*new)
    

        #print(f"{GCD=}")
        #print(f"{ratio=}")
        total += GCD

        f.readline()
    average = total / it
    averageR = totalR / it
    print(f"Average GCD = {average}")
    print(f"Average log Ratio = {averageR}")



