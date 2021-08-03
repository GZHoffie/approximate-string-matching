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
        self.currentPosition = (0, 0)
        self.nearestHighways = [0] * (2 * self.k + 1)
        self.sight = sight
        self.destinationLane = - abs(self.m - self.n)
    

    def updateNearestHighway(self):
        currentLane, currentCol = self.currentPosition
        for l in range(-self.k, self.k + 1):
            if self.nearestHighways[l + self.k] < self.currentPosition + leapForwardColumn(currentLane, l):
                self.nearestHighways[l + self.k] = self.hurdleBits.findFirstHighway(l, currentCol)

    def findBestHighwayNearby(self):
        """
        Iterate through all lanes and find the best highway within sight.
        """
        currentLane, currentCol = self.currentPosition
        maxScore = float('-inf')
        bestLane = None
        for l in range(-self.k, self.k + 1):
            if self.nearestHighways[l + self.k] - self.currentPosition[1] <= self.sight:
                score = self.hurdleBits.getFirstHighwayLength(l) - leapLanePenalty(currentLane, l) -\
                    leapLanePenalty(l, self.destinationLane)
                if score > maxScore:
                    maxScore = score
                    bestLane = l
        
        if bestLane is None:
            bestLane = np.argmin(self.nearestHighways) - self.k
        
        return bestLane
    

    def projectHighway(self, destinationLane):
        """
        Project the highway onto the destination lane
        """
        l_, c_ = self.currentPosition
        if l_ < destinationLane:
            leapingLanes = list(range(l_, destinationLane + 1))
        else:
            leapingLanes = list(range(l_, destinationLane - 1, -1))

        highwayLength = self.hurdleBits.getFirstHighwayLength(l_)
        mask = 1 << highwayLength
        for l in leapingLanes:
            self.hurdleBits.bits[l] = self.hurdleBits.bits[l] & (mask << leapForwardColumn(l_, l))
            newHighwayLength = self.hurdleBits.getFirstHighwayLength(l)
            mask = 1 << newHighwayLength
            l_ = l



        




    


            

if __name__ == "__main__":
    GASMAProjection("AGAGCTAAACATGGCCGCACATAAATCGTTTTGAGTTGAAACTTTACCGCTGCATCTATTTTTCTCCTAGAATTATACCGTACACAGCCGACGTTCCACC", 
              "AGAGCTAAACAAGGGGCCCACATTAACGTTTTGAGCTTGAAGATCTTTACCGCGATCTATTTTTTCTCCTAGATTACCGTACACACCGACACTTCCATC", debug=True)

