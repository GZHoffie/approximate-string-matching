from pymatch.algorithms import GASMAProjection, NeedlemanWunsch
from pymatch.metrics import Coverage, numLeaps
import time
from collections import Counter
import matplotlib.pyplot as plt

test_file = "/home/zhenhao/approximate-string-matching/pymatch/test/resource/sample.random.dataset.seq"
test_items = 1000
good_items = 0

GASMATime = 0
NWTime = 0
error = 0
correct = 0

# Counter
totalCounter = Counter()
correctCounter = Counter()
coverCounter = Counter()
leap1Counter = Counter()
leap2Counter = Counter()

with open(test_file, "r") as f:
    i = 0
    while True:
        print(i)
        line = f.readline()
        if not line:
            break
        str1 = line[1:][:-1]
        str2 = f.readline()[1:][:-1]

        print(str1)
        print(str2)
        
        currTime = time.time()
        g1 = GASMAProjection(str1, str2, k=2, sight=7, debug=False)
        #g2 = GASMAProjection(str1, str2, k=2, sight=7, maxZerosIgnored=2, debug=False)
        #g3 = GASMAProjection(str1, str2, k=2, sight=7, skipHurdle=True, debug=False)

        cost1 = g1.editDistance()
        #cost2 = g2.editDistance()
        cost = cost1#min(cost1, cost2)
        GASMATime += time.time() - currTime


        currTime = time.time()
        nw = NeedlemanWunsch(str1, str2)
        NWcost = nw.editDistance()
        nw.backtrack()
        
        NWTime += time.time() - currTime


        error += abs(cost - NWcost)
        correct += (cost == NWcost)

        # Counter
        totalCounter[NWcost] += 1
        correctCounter[NWcost] += (cost == NWcost)


        #if NWcost != cost:
            #print(nw)
            #print(g1.match["dna1"])
            #print(g1.match["dna2"])
            
            #print("NW Cost:", NWcost)
            #print("GASMA Cost:", cost)
        

        good_items += Coverage(g1.match, nw.alignment, 0, 3).compute()

        coverCounter[NWcost] += Coverage(g1.match, nw.alignment, 0, 3).compute()
        leap1Counter[NWcost] += numLeaps(g1.match).compute()
        leap2Counter[NWcost] += numLeaps(nw.alignment).compute()



        i += 1
        if i >= test_items:
            break

print("GASMA Time:", GASMATime)
print("NW Time:", NWTime)
print("MAE:", error / test_items)
print("Correct rate:", correct / test_items)
print("Better match:", good_items / test_items)

for i in totalCounter:
    correctCounter[i] /= totalCounter[i]
    coverCounter[i] /= totalCounter[i]
    leap1Counter[i] /= totalCounter[i]
    leap2Counter[i] /= totalCounter[i]

plt.bar(correctCounter.keys(), correctCounter.values())

print("Correct rate", correctCounter)
print("coverage", coverCounter)
print("num leaps in greedy", leap1Counter)
print("num leaps in optimal", leap2Counter)
