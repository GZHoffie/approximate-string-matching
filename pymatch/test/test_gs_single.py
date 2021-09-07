from pymatch.algorithms import GASMAProjection, NeedlemanWunsch
import time

test_file = "/home/zhenhao/approximate-string-matching/pymatch/test/resource/sample.random.dataset.seq"
test_items = 100

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
        g = GASMAProjection(str1, str2, k=8, sight=3, debug=False)
        cost = g.editDistance()
        GASMATime += time.time() - currTime

        currTime = time.time()
        nw = NeedlemanWunsch(str1, str2)
        NWcost = nw.editDistance()
        NWTime += time.time() - currTime

        error += abs(cost - NWcost)
        correct += (cost == NWcost)

        
        
        print("NW Cost:", NWcost)
        print("GASMA Cost:", cost)

        i += 1
        if i >= test_items:
            break

print("GASMA Time:", GASMATime)
print("NW Time:", NWTime)
print("MAE:", error / test_items)
print("Correct rate:", correct / test_items)