from pymatch.util import ApproximateStringMatching, HurdleBits
from pymatch.algorithms import GASMA
import numpy as np
import gmpy


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



class GASMAProjection(ApproximateStringMatching):
    def __init__(self, dna1, dna2, k, mismatchCost=None, insertCost=None, deleteCost=None, sight=3):
        super().__init__(dna1, dna2, mismatchCost=mismatchCost, insertCost=insertCost, deleteCost=deleteCost)
        self.hurdleBits = HurdleBits(dna1=dna1, dna2=dna2, k=k, maxOnesIgnored=1)
        self.hurdleBits.shiftRight(1)
        self.length = len(format(self.hurdleBits.bits[0], 'b'))
        self.currentPosition = (0, 0)
        self.k = k
        self.nearestHighways = [0] * (2 * self.k + 1)
        self.sight = sight
        self.destinationLane = - abs(self.m - self.n)
        self.i = 0
        self.j = 0
        self.match = {"dna1": "", "dna2": ""}
        
    

    def updateNearestHighway(self):
        currentLane, currentCol = self.currentPosition
        for l in range(-self.k, self.k + 1):
            if self.nearestHighways[l + self.k] <= self.currentPosition[1] + leapForwardColumn(currentLane, l):
                self.nearestHighways[l + self.k] = self.hurdleBits.findFirstHighway(l, currentCol + leapForwardColumn(currentLane, l))
        print(self.nearestHighways)
    

    def extractHighway(self, lane, pos, length):
        """
        Extract the highway that starts from pos and at lane.
        """
        string = "1" * (self.length - length - pos) + "0" * length + "1" * pos
        return int(string, 2) | self.hurdleBits.bits[lane + self.k], int(string, 2)


    def findBestHighwayNearby(self):
        """
        Iterate through all lanes and find the best highway within sight.
        """
        currentLane, currentCol = self.currentPosition
        maxScore = float('-inf')
        bestLane = None
        bestHighway = None
        for l in range(-self.k, self.k + 1):
            if self.nearestHighways[l + self.k] - currentCol <= self.sight:
                length = self.hurdleBits.getFirstHighwayLength(l, currentCol)
                highway, processedHighway = self.extractHighway(l, self.nearestHighways[l + self.k], length)
                print(format(highway, "b"))

                hurdlesToCross = self.nearestHighways[l + self.k] - currentCol - leapForwardColumn(currentLane, l)
                score = format(highway, "b").count("0") - leapLanePenalty(currentLane, l) -\
                    leapLanePenalty(l, self.destinationLane) - hurdlesToCross
                print(l, score)
                if score > maxScore:
                    maxScore = score
                    bestHighway = (highway, processedHighway)
                    bestLane = l
        
        if bestLane is None:
            bestLane = np.argmin(self.nearestHighways) - self.k
        
        return bestLane, bestHighway[0], bestHighway[1]
    

    def print(self):
        print("bits")
        for l in range(2 * self.k + 1):
            print(format(self.hurdleBits.bits[l], 'b'))
        
        print("processed bits")
        for l in range(2 * self.k + 1):
            print(format(self.hurdleBits.processedBits[l], 'b'))


    def _updateMatch(self, l, l_, length):
        """
        Update the match between two strings according to current lane `l`, the lane on the next step `l_`
        and the length of the highway.
        """
        leapingLanes = abs(l_ - l)
        if l_ > l:
            for _ in range(leapingLanes):
                self.match["dna1"] += "-"
                self.match["dna2"] += self.dna2.string[self.j]
                self.j += 1
        else:
            for _ in range(leapingLanes):
                self.match["dna2"] += "-"
                self.match["dna1"] += self.dna1.string[self.i]
                self.i += 1
        
        for _ in range(length):
            self.match["dna1"] += self.dna1.string[self.i]
            self.match["dna2"] += self.dna2.string[self.j]
            self.j += 1
            self.i += 1



    def editDistance(self):

        # Store leap and hurdle costs
        leapCost = 0
        hurdleCost = 0

        # Find the first highway
        self.updateNearestHighway()
        lane, highway, processedHighway = self.findBestHighwayNearby()
        print(lane, highway, processedHighway)

        highwayLength = self.hurdleBits.getFirstHighwayLength(self.currentPosition[0], self.currentPosition[1])
        leapCol = leapForwardColumn(self.currentPosition[0], lane)
        colAfterLeap = self.currentPosition[1] + leapCol
        
        leapCost += leapLanePenalty(self.currentPosition[0], lane)
        self.currentPosition = (lane, colAfterLeap)
        


        if (self.hurdleBits.processedBits[lane + self.k] >> colAfterLeap) & 1 == 0:
            # We have overlap between the two lanes
            currHighway = self.hurdleBits.extractHighway(self.currentPosition[0], self.currentPosition[1])
            _highway = self.hurdleBits.extractHighway(lane, colAfterLeap)
            _highwayStartCol = gmpy.scan0(_highway)
            _overlap = format(_highway, 'b')[(self.length-_highwayStartCol):(self.length-colAfterLeap)]
            overlap = format(currHighway, 'b') [(self.length-_highwayStartCol+leapCol):(self.length-colAfterLeap+leapCol)]
            print(overlap, _overlap)
            

        




    


            

if __name__ == "__main__":
    g = GASMAProjection("AGAGCTAAACATGGCCGCACATAAATCGTTTTGAGTTGAAACTTTACCGCTGCATCTATTTTTCTCCTAGAATTATACCGTACACAGCCGACGTTCCACC", 
              "AGAGCTAAACAAGGGGCCCACATTAACGTTTTGAGCTTGAAGATCTTTACCGCGATCTATTTTTTCTCCTAGATTACCGTACACACCGACACTTCCATC", k=2)
    g.editDistance()

