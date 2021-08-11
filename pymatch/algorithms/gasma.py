from pymatch.util import ApproximateStringMatching, HurdleBits
from pymatch.algorithms import GASMA
import numpy as np


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


    def projectHighway(self, destinationLane, highway, processedHighway):
        """
        Project the highway onto the destination lane
        """
        print("destination lane", destinationLane)
        l_, c_ = self.currentPosition
        if l_ < destinationLane:
            leapingLanes = list(range(l_ + 1, destinationLane + 1))
        else:
            leapingLanes = list(range(l_ - 1, destinationLane - 1, -1))

        highwayLength = self.hurdleBits.getFirstHighwayLength(l_, c_)
        print("column", c_, "length", highwayLength)
        
        bitMask = highway
        processedBitMask = processedHighway
        newLength = 0
        for l in leapingLanes:
            #bitMask = highway#self.hurdleBits.bits[l_ + self.k] & ((((1 << highwayLength + self.currentPosition[1]) - 1) << 1) + 1)
            print(format(bitMask, 'b'))
            #processedBitMask = processedHighway
            self.hurdleBits.bits[l + self.k] = self.hurdleBits.bits[l + self.k] & (bitMask << leapForwardColumn(l_, l))
            self.hurdleBits.processedBits[l + self.k] = self.hurdleBits.processedBits[l + self.k] & (processedBitMask << leapForwardColumn(l_, l))
            c_ += leapForwardColumn(l_, l)
            #self.hurdleBits.processedBits = self.hurdleBits.removeSingleOnes(self.hurdleBits.processedBits)
            newLength = self.hurdleBits.getFirstHighwayLength(l, c_)
            bitMask, processedBitMask = self.extractHighway(l, self.nearestHighways[l + self.k], newLength)
            #highwayLength = self.hurdleBits.getFirstHighwayLength(l, c_ + leapForwardColumn(l_, l))
            #mask = 1 << newHighwayLength
            #self.currentPosition = (l, )
            l_ = l
            #self.currentPosition = (l_, c_)
            self.print()
        
        #self.currentPosition = (l_, c_ + newLength)

            


    def editDistance(self):
        self.updateNearestHighway()
        lane, highway, processedHighway = self.findBestHighwayNearby()
        self.currentPosition = (lane, self.currentPosition[1] + self.hurdleBits.getFirstHighwayLength(self.currentPosition[0], self.currentPosition[1]) + leapForwardColumn(self.currentPosition[0], lane))
        hurdleCost = 0
        leapCost = leapLanePenalty(self.currentPosition[0], lane)
        currentPos = (0, 0)
        while self.currentPosition[1] < self.length - 1 or min(self.nearestHighways) >= self.length:
            self.updateNearestHighway()
            lane_, highway_, processedHighway_ = self.findBestHighwayNearby()
            self.projectHighway(lane_, highway, processedHighway)
            
            self.currentPosition = (lane_, self.currentPosition[1] + self.hurdleBits.getFirstHighwayLength(lane, self.currentPosition[1]) + leapForwardColumn(lane, lane_))
            print("res", format(self.hurdleBits.bits[lane_ + self.k], "b")[::-1][(currentPos[1] + leapForwardColumn(lane, lane_)):self.currentPosition[1]])
            hurdleCost += format(self.hurdleBits.bits[lane_ + self.k], "b")[::-1][(currentPos[1] + leapForwardColumn(lane, lane_)):self.currentPosition[1]].count("1")
            print("position", self.currentPosition)
            leapCost += leapLanePenalty(lane, lane_)

            lane, highway, processedHighway = lane_, highway_, processedHighway_
            currentPos = self.currentPosition
        #while self.currentPosition[1] != self.length:

        print("hurdle cost:", hurdleCost)
        print("leap cost:", leapCost)



        




    


            

if __name__ == "__main__":
    g = GASMAProjection("AGAGCTAAACATGGCCGCACATAAATCGTTTTGAGTTGAAACTTTACCGCTGCATCTATTTTTCTCCTAGAATTATACCGTACACAGCCGACGTTCCACC", 
              "AGAGCTAAACAAGGGGCCCACATTAACGTTTTGAGCTTGAAGATCTTTACCGCGATCTATTTTTTCTCCTAGATTACCGTACACACCGACACTTCCATC", k=2)
    g.editDistance()

