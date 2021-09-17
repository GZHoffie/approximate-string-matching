from pymatch.algorithms import GASMAProjection, NeedlemanWunsch
import time

test_file = "/home/zhenhao/approximate-string-matching/pymatch/test/resource/sample.random.dataset.seq"
test_items = 1000

GASMATime = 0
NWTime = 0
error = 0
correct = 0

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
        g2 = GASMAProjection(str1, str2, k=2, maxZerosIgnored=2, debug=False)

        cost1 = g1.editDistance()
        cost2 = g2.editDistance()
        cost = min(cost1, cost2)
        GASMATime += time.time() - currTime

        currTime = time.time()
        nw = NeedlemanWunsch(str1, str2)
        NWcost = nw.editDistance()
        NWTime += time.time() - currTime

        error += abs(cost - NWcost)
        correct += (cost == NWcost)

        if NWcost != cost:
            
            print("NW Cost:", NWcost)
            print("GASMA Cost:", cost)

        i += 1
        if i >= test_items:
            break

print("GASMA Time:", GASMATime)
print("NW Time:", NWTime)
print("MAE:", error / test_items)
print("Correct rate:", correct / test_items)