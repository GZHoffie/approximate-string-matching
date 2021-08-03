from pymatch.util import ApproximateStringMatching, HurdleMatrix
import time

class GASMA(ApproximateStringMatching):
    def __init__(self, dna1, dna2, k, leapCost=None, hurdleCost=1, threshold=3, crossHurdleThreshold=1, reverse=False, debug=False):
        if len(dna1) > len(dna2):
            # swap the two
            temp = dna1
            dna1 = dna2
            dna2 = temp
        
        candidate = ["A", "C", "G", "T"]
        candidate.pop(candidate.index(dna1[0]))
        if dna2[0] != dna1[0]:
            candidate.pop(candidate.index(dna2[0]))
        appendix = candidate[0] * (threshold + 5)

        dna1 = appendix + dna1 + appendix
        dna2 = appendix + dna2 + appendix
        
        
        super().__init__(dna1, dna2)
        self.destinationLane = - abs(self.m - self.n)
        k = max(k, abs(self.m - self.n) + 5)
        #assert k >= abs(self.m - self.n), "k is less than the difference in length of" \
        #    " the two DNAs"
        self.hurdleMatrix = HurdleMatrix(dna1, dna2, k, mismatchCost=hurdleCost, leapCost=leapCost, threshold=threshold, crossHurdleThreshold=crossHurdleThreshold, reverse=reverse, debug=debug)
        self.k = k
        #self.highways = self.hurdleMatrix.highways
        self.debug = debug
        #print(self.highways)
        self.hurdleCost = 1
        self.leapCost = 1
        self.matrixLength = len(format(self.hurdleMatrix.hurdleMatrix[0], 'b'))
        self.threshold = threshold
        #curr = time.time()
        #self.editDistance()
        #print("find route time:", time.time() - curr)

        


    def editDistance(self):
        bestHighways = self.findBestHighways()
        bestHighways_transformed = [(shift, start + length - 1, length, hurdles) for shift, start, length, hurdles in bestHighways]
        if self.debug:
            print("best highways:", bestHighways_transformed)
        route, hurdleCost, leapCost = self.linkHighways(bestHighways_transformed)
        if self.debug:
            print(route)
            print("hurdle cost:", hurdleCost)
            print("leap cost:", leapCost)
            print("total cost:", leapCost + hurdleCost)
        return leapCost + hurdleCost

        
    def transformHighways(self):
        transformedHighways = []
        for shift, start, length, hurdles in self.highways:
            bitlist = [1] * (self.n)
            for i in range(length):
                bitlist[start + i] = 0
            transformedHighways.append((shift, bitlist, len(hurdles)))
        
        return transformedHighways

    
    def findBestHighways(self):
        def score(L):
            bitlist = [1] * (self.m + 1)
            for _, l, _ in L:
                bitlist = [x & y for x, y in zip(bitlist, l)]
            
            numZeros = bitlist.count(0)
            numHighways = len(L)
            return numZeros - 3 * numHighways# - sum([h[2] for h in L])
        
        L = []
        L_info = []
        transformedHighways = self.transformHighways()
        L0 = transformedHighways.copy()
        L0_info = self.highways.copy()
        while len(L0) > 0:
            currScore = score(L)
            maxScore = float('-inf')
            i = 0
            bestHighway = None
            for l in L0:
                lscore = score(L + [l]) - currScore
                if lscore > maxScore or (lscore == maxScore and l[0] == self.destinationLane):
                    maxScore = lscore
                    bestHighway = i
                
                i += 1

            if maxScore < 0:
                break
            else:
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
                highwayDict[h[0]] = [(h[1], h[2], h[3])]
            else:
                highwayDict[h[0]] += [(h[1], h[2], h[3])]
        
        for shift in highwayDict:
            highwayDict[shift] = sorted(highwayDict[shift], key=lambda x: x[0], reverse=True)

        currentPosition = (0, self.n)
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
                columnAfterLeap = currentPosition[1] - leapForwardColumn(currentPosition[0], shift)
                #hurdleCross = currentPosition[1] - columnAfterLeap
                hurdleCross = format(self.hurdleMatrix.hurdleMatrix[shift + self.k], 'b')[self.matrixLength-columnAfterLeap:self.matrixLength-highwayDict[shift][0][0]-1].count('1') + \
                    len([hurdle for hurdle in highwayDict[shift][0][2] if hurdle < columnAfterLeap])
                if self.debug:
                    print('shift', shift, columnAfterLeap, highwayDict[shift][0], hurdleCross)
                    print(format(self.hurdleMatrix.hurdleMatrix[shift + self.k], 'b')[self.matrixLength-columnAfterLeap:self.matrixLength-highwayDict[shift][0][0]-1])
                if hurdleCross < leastHurdleCross:
                    bestShift = shift
                    leastHurdleCross = hurdleCross
                    colAfterLeap = columnAfterLeap
            
            if self.debug:
                print('position:',currentPosition, 'col after leap:', colAfterLeap, 'best shift:', bestShift, highwayDict[bestShift][0], 'cost:', leastHurdleCross)
            
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
        


        if currentPosition != (self.destinationLane, 0):
            leapCost += leapLanePenalty(currentPosition[0], self.destinationLane)
            columnAfterLeap = currentPosition[1] - leapForwardColumn(currentPosition[0], self.destinationLane)
            if columnAfterLeap > 0:
                format(self.hurdleMatrix.hurdleMatrix[self.destinationLane + self.k], 'b')[0:columnAfterLeap].count('1')
            route += [(self.destinationLane, 0)]

        """
        if currentPosition[1] != 0:
            hurdleCost += currentPosition[1]
        """
        
        return route, hurdleCost, leapCost


if __name__ == "__main__":
    g = GASMA("CAAGTATATCAGCAGGTAGCCAGCCGCCCTTGAACTACAATCTACCGCGATACGCTATAGTACAAGACTCGGGCTAGGCCTCCATTAAGAGGCTGGCTTC", 
              "CAAGTATATCGCAGGTAGCCAGCCGCCCTTGAACTACAATCTACCGCGATACGCTATAGTACAAGACTCGGGCTAGTGCCTCCATTAAGAGGCTGGCTTC", 7, threshold=3, crossHurdleThreshold=0, debug=True)
    print(g.editDistance())
    #import time
    #a = time.time()
    #g = GASMA("GAGAACCAATCAGCACAGGGCACTCTATGTAATTCTCGAGGCGATTGACCGTCTGGTTGCGGGGCTGTGGCAATCTTTTAAGAGGGCCGTGCCATTACTG",
    #          "GAGAACCAATCAGCACAGTGCACTCTATGTAATTCTCGAGGCGATTGACCGTCTGGTTGCGGGGCTGTGGCAATCTTTTAAGAGGGCCGTGCCATTACTG",
    #          2, threshold=3)
    #print("total time:", time.time() - a)
    
    #print(g.findBestHighways())

    
    