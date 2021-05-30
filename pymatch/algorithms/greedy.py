from pymatch.util import ApproximateStringMatching, HurdleMatrix
import time

class GASMA(ApproximateStringMatching):
    def __init__(self, dna1, dna2, k, leapCost=None, hurdleCost=1, threshold=3):
        if len(dna1) > len(dna2):
            # swap the two
            temp = dna1
            dna1 = dna2
            dna2 = temp
        
        super().__init__(dna1, dna2)
        self.hurdleMatrix = HurdleMatrix(dna1, dna2, k, mismatchCost=hurdleCost, leapCost=leapCost, threshold=threshold)
        self.k = k
        self.highways = self.hurdleMatrix.highways
        self.transformedHighways = self.transformHighways()
        self.hurdleCost = 1
        self.leapCost = 1
        #curr = time.time()
        #self.editDistance()
        #print("find route time:", time.time() - curr)

        


    def editDistance(self):
        bestHighways = self.findBestHighways()
        bestHighways_transformed = [(shift, start + length - 1, length) for shift, start, length in bestHighways]
        #print(bestHighways_transformed)
        route, hurdleCost, leapCost = self.linkHighways(bestHighways_transformed)
        #print(route)
        #print("hurdle cost:", hurdleCost)
        #print("leap cost:", leapCost)
        #print("total cost:", leapCost + hurdleCost)
        return leapCost + hurdleCost

        
    def transformHighways(self):
        transformedHighways = []
        for shift, start, length in self.highways:
            bitlist = [1] * (self.m)
            for i in range(length):
                bitlist[start + i] = 0
            transformedHighways.append((shift, bitlist))
        
        return transformedHighways

    
    def findBestHighways(self):
        def score(L):
            bitlist = [1] * (self.m + 1)
            for _, l in L:
                bitlist = [x & y for x, y in zip(bitlist, l)]
            
            numZeros = bitlist.count(0)
            numHighways = len(L)
            return numZeros - 2 * numHighways
        
        L = []
        L_info = []
        L0 = self.transformedHighways.copy()
        L0_info = self.highways.copy()
        while len(L0) > 0:
            currScore = score(L)
            scoreList = []
            for l in L0:
                scoreList.append(score(L + [l]) - currScore)
            
            maxScore = max(scoreList)
            if maxScore < 0:
                break
            else:
                bestHighway = scoreList.index(maxScore)
                L.append(L0.pop(bestHighway))
                L_info.append(L0_info.pop(bestHighway))
        
        return L_info
    
    def linkHighways(self, highways):
        def leapForwardColumn(l_, l):
            """
            Returns the number of columns the toad moves forward when leaping from lane l_ to l. 
            When l and l_ are the same lane, then the number should just be 1.
            """
            if l_ == l:
                return 0 #1 if pos < self.m else 0
            elif abs(l_) > abs(l) and l * l_ >= 0:
                return 0
            elif abs(l_) < abs(l) and l * l_ >= 0:
                return abs(l - l_)
            else:
                return abs(l - l_) - abs(l_)
        
        def leapLanePenalty(l_, l):
            """
            Returns the penalty of leaping from lane l_ to l. When l and l_ are
            the same lane, then the penalty is just the energy cost of next hurdle.
            """
            if l_ == l:
                return 0
            else:
                return 1 * abs(l_ - l)
        
        highwayDict = {}
        for h in highways:
            if h[0] not in highwayDict:
                highwayDict[h[0]] = [(h[1], h[2])]
            else:
                highwayDict[h[0]] += [(h[1], h[2])]
        
        for shift in highwayDict:
            highwayDict[shift] = sorted(highwayDict[shift], key=lambda x: x[0], reverse=True)

        currentPosition = (0, self.m - 1)
        numHighways = len(highways)
        route = [currentPosition]
        hurdleCost = 0
        leapCost = 0
        colAfterLeap = None
        while numHighways > 0:
            leastHurdleCross = float('inf')
            bestShift = None
            for shift in highwayDict:
                if len(highwayDict[shift]) == 0:
                    continue
                columnAfterLeap = leapForwardColumn(currentPosition[0], shift) + highwayDict[shift][0][0] + 1
                hurdleCross = currentPosition[1] - columnAfterLeap
                if hurdleCross < leastHurdleCross:
                    bestShift = shift
                    leastHurdleCross = hurdleCross
                    colAfterLeap = columnAfterLeap
            
            #print(leastHurdleCross)
            
            if bestShift is None:
                break
            else:
                if leastHurdleCross > 0:
                    hurdleCost += leastHurdleCross
                leapCost += leapLanePenalty(currentPosition[0], bestShift) 
                #print(leapCost)
                bestHighway = highwayDict[bestShift].pop(0)
                numHighways -= 1
                highwayEnd = (bestShift, bestHighway[0] - bestHighway[1] + 1)
                route += [(bestShift, colAfterLeap), highwayEnd]
            
            currentPosition = highwayEnd
        

        """
        if currentPosition != (0, 0):
            leapCost += leapLanePenalty(currentPosition[0], 0)
            columnAfterLeap = currentPosition[1] - leapForwardColumn(currentPosition[0], 0)
            if columnAfterLeap > 0:
                hurdleCost += 1 * columnAfterLeap
            route += [(0, 0)]
        """

        if currentPosition[1] != 0:
            hurdleCost += currentPosition[1]
        
        return route, hurdleCost, leapCost


if __name__ == "__main__":
    g = GASMA("ACTAGAACTT", "ACTTAGCACT", 2, 2)
    print(g.editDistance())
    #import time
    #a = time.time()
    #g = GASMA("GAGAACCAATCAGCACAGGGCACTCTATGTAATTCTCGAGGCGATTGACCGTCTGGTTGCGGGGCTGTGGCAATCTTTTAAGAGGGCCGTGCCATTACTG",
    #          "GAGAACCAATCAGCACAGTGCACTCTATGTAATTCTCGAGGCGATTGACCGTCTGGTTGCGGGGCTGTGGCAATCTTTTAAGAGGGCCGTGCCATTACTG",
    #          2, threshold=3)
    #print("total time:", time.time() - a)
    
    #print(g.findBestHighways())

    
    