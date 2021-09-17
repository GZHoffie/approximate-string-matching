from pymatch.algorithms import GASMAProjection, NeedlemanWunsch
import time

test_file = "/home/zhenhao/approximate-string-matching/pymatch/test/resource/sample.random.dataset.seq"
ref_file = "/home/zhenhao/approximate-string-matching/pymatch/test/resource/nw.txt"
test_items = 10000

GASMATime = 0
NWTime = 0
error = 0
correct = 0
upper_limit = 20

errorDict = {n: 0 for n in range(upper_limit)}
correctDict = {n: 0 for n in range(upper_limit)}
countDict = {n: 0 for n in range(upper_limit)}

with open(ref_file, "r") as rf:
    with open(test_file, "r") as f:
        i = 0
        while True:
            line = f.readline()
            if not line:
                break
            str1 = line[1:][:-1]
            str2 = f.readline()[1:][:-1]


            currTime = time.time()
            g1 = GASMAProjection(str1, str2, k=2, sight=7, debug=False)
            g2 = GASMAProjection(str1, str2, k=2, sight=7, maxZerosIgnored=2, debug=False)
            g3 = GASMAProjection(str1, str2, k=2, sight=7, skipHurdle=True, debug=False)

            cost1 = g1.editDistance()
            cost2 = g2.editDistance()
            cost3 = g3.editDistance()
            cost = min(cost1, cost2, cost3)
            GASMATime += time.time() - currTime

            #currTime = time.time()
            #nw = NeedlemanWunsch(str1, str2)
            #NWcost = nw.editDistance()
            #NWTime += time.time() - currTime
            NWcost = int(float(rf.readline()[:-1]))

            #error += abs(cost - NWcost)
            #correct += (cost == NWcost)
            errorDict[NWcost] += abs(cost - NWcost)
            correctDict[NWcost] += (cost == NWcost)
            countDict[NWcost] += 1

            if cost != NWcost and NWcost <= 5:
                print(str1)
                print(str2)
                print("NW Cost:", NWcost)
                print("GASMA Cost:", cost)

            i += 1
            if i >= test_items:
                break

print("GASMA Time:", GASMATime)
#print("MAE:", error / test_items)
#print("Correct rate:", correct / test_items)

MAE = [errorDict[c] / countDict[c] for c in range(upper_limit)]
correct_rate = [correctDict[c] / countDict[c] for c in range(upper_limit)]
print(countDict)

print(MAE)
print(correct_rate)

import matplotlib.pyplot as plt
plt.figure()
plt.plot(MAE[0:18])
plt.ylabel('MAE')
plt.xlabel('Optimal Edit Distance')
plt.savefig('/home/zhenhao/approximate-string-matching/pymatch/test/asset/MAE.png')

plt.figure()
plt.plot(correct_rate[0:18])
plt.ylabel('Correct Rate')
plt.xlabel('Optimal Edit Distance')
plt.savefig('/home/zhenhao/approximate-string-matching/pymatch/test/asset/CR.png')

import numpy as np
plt.figure()
plt.plot(np.array(MAE[0:18]) / np.array(range(18)))
plt.ylabel('Relative Error')
plt.xlabel('Optimal Edit Distance')
plt.savefig('/home/zhenhao/approximate-string-matching/pymatch/test/asset/CR_relative.png')