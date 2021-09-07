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
    def __init__(self, dna1, dna2, k, mismatchCost=None, insertCost=None, deleteCost=None, sight=3, debug=False):
        super().__init__(dna1, dna2, mismatchCost=mismatchCost, insertCost=insertCost, deleteCost=deleteCost)

        k = max(k, abs(len(self.dna1.string) - len(self.dna2.string)) + 2)

        self.hurdleBits = HurdleBits(dna1=dna1, dna2=dna2, k=k, maxOnesIgnored=1, debug=debug)
        #self.hurdleBits.shiftRight(1)
        self.length = len(format(self.hurdleBits.bits[0], 'b'))
        self.currentPosition = (0, 1)
        self.k = k
        self.nearestHighways = [0] * (2 * self.k + 1)
        self.sight = sight
        self.destinationLane = - abs(self.m - self.n)
        self.i = 0
        self.j = 0
        self.match = {"dna1": "", "dna2": ""}

        # Edit distance
        self.hurdleCost = 0
        self.leapCost = 0

        self.debug = debug
        
    

    def updateNearestHighway(self, pos):
        currentLane, currentCol = pos
        for l in range(-self.k, self.k + 1):
            if self.nearestHighways[l + self.k] <= currentCol + leapForwardColumn(currentLane, l):
                #print(self.nearestHighways, currentCol + leapForwardColumn(currentLane, l))
                self.nearestHighways[l + self.k] = self.hurdleBits.findFirstHighway(l, currentCol + leapForwardColumn(currentLane, l))
        #print(self.nearestHighways)
    
    

    def extractHighway(self, lane, pos, length):
        """
        Extract the highway that starts from pos and at lane.
        """
        string = "1" * (self.length - length - pos) + "0" * length + "1" * pos
        return int(string, 2) | self.hurdleBits.bits[lane + self.k], int(string, 2)


    def findBestHighwayNearby(self, pos):
        """
        Iterate through all lanes and find the best highway within sight.
        """
        currentLane, currentCol = pos
        self.updateNearestHighway(pos)
        maxScore = float('-inf')
        bestLane = None
        bestHighway = None
        for l in range(-self.k, self.k + 1):
            if self.nearestHighways[l + self.k] - currentCol <= self.sight:
                length = self.hurdleBits.getFirstHighwayLength(l, currentCol + leapForwardColumn(currentLane, l))
            
                highway = self.hurdleBits.extractHighway(l, self.nearestHighways[l + self.k])
                if self.debug:
                    print(format(highway, "b"), length)

                hurdlesToCross = self.nearestHighways[l + self.k] - currentCol - leapForwardColumn(currentLane, l)
                score = format(highway >> (currentCol + leapForwardColumn(currentLane, l)), "b").count("0") - leapLanePenalty(currentLane, l) -\
                    (leapLanePenalty(l, self.destinationLane) - leapLanePenalty(currentLane, self.destinationLane)) - hurdlesToCross
                if self.debug:
                    print("highway on lane", l, "score:", score)
                if score > maxScore:
                    maxScore = score
                    bestHighway = (highway, length)
                    bestLane = l
        
        if bestLane is None:
            l = np.argmin(self.nearestHighways) - self.k
            length = self.hurdleBits.getFirstHighwayLength(l, currentCol + leapForwardColumn(currentLane, l))
            highway = self.hurdleBits.extractHighway(l, self.nearestHighways[l + self.k])
            return l, highway, length
        
        return bestLane, bestHighway[0], bestHighway[1]
    
    def findBestShortHighwayNearby(self, pos, targetLane, maxCol):
        """
        If leaping through multiple lanes, find if we can take short highways between the lanes we are leaping. 
        If there is such small highway, return the small highway. Otherwise return None
        """
        currentLane, currentCol = pos
        maxScore = float('-inf')
        bestLane = None
        bestHighway = None

        # Determine the range of lanes
        if pos[0] > targetLane:
            leapingLanes = range(targetLane, pos[0] + 1)

        else:
            leapingLanes = range(pos[0], targetLane + 1)

        for l in leapingLanes:
            nearestHighwayCol = gmpy.scan0(self.hurdleBits.bits[l + self.k] >> (currentCol + leapForwardColumn(currentLane, l))) + currentCol + leapForwardColumn(currentLane, l)
            #print(l, nearestHighwayCol)
            if nearestHighwayCol < maxCol - leapForwardColumn(l, targetLane):
                
                length = gmpy.scan1(self.hurdleBits.bits[l + self.k] >> nearestHighwayCol)
                if self.debug:
                    print("small highway found:", l, nearestHighwayCol, length)
                #length = self.hurdleBits.getFirstHighwayLength(l, currentCol)
                #highway, processedHighway = self.extractHighway(l, self.nearestHighways[l + self.k], length)
                #print(format(highway, "b"))

                #hurdlesToCross = self.nearestHighways[l + self.k] - currentCol - leapForwardColumn(currentLane, l)
                dist = nearestHighwayCol - currentCol - leapForwardColumn(currentLane, l) + leapLanePenalty(currentLane, l)
                if dist > 2:
                    score = - dist
                else:
                    score = length - dist
                if self.debug:
                    print(l, score, "dist", dist)
                if score > maxScore:
                    maxScore = score
                    bestHighway = (nearestHighwayCol, nearestHighwayCol + length)
                    bestLane = l
        

        
        return bestLane, bestHighway


    def decideWhereToLeap(self, highwayA, highwayB, laneA, laneB, areAdjacentLanes=True):
        """
        Given two highways in adjacent lanes, represented by bits, we decide the point where we leap.
        We can either leap at the end of first highway, or the start of the second highway.

        If the two highways have no overlap, then we simply choose the end of first highway (no difference
        in edit distance).
        If there is overlap, we look at the bit vector in the overlap of the two lanes, and choose the one
        with more zeros to go.
        
        returns the column number on which to leap from laneA to laneB.
        """
        dColumn = leapForwardColumn(laneA, laneB)
        if self.debug:
            print("dcolumn", dColumn)
        
        if areAdjacentLanes:
            choiceA = self.length - gmpy.scan0(int(format(highwayA, "b")[::-1], 2)) # End of highway A    
            choiceB = gmpy.scan0(highwayB) - dColumn # Start of highway B
        else:
            highwayBReversed = int(format(highwayB, "b")[::-1], 2)
            choiceA = self.length - (gmpy.scan1(highwayBReversed >> gmpy.scan0(highwayBReversed)) + gmpy.scan0(highwayBReversed)) - dColumn # after the last hurdle of highway B
            firstHurdle = max(self.currentPosition[1] + leapForwardColumn(self.currentPosition[0], laneA), gmpy.scan0(highwayA))
            choiceB = gmpy.scan1(highwayA >> firstHurdle) + firstHurdle # Before the First hurdle of highway A
        

        #choiceA = max(self.currentPosition[1] + leapForwardColumn(self.currentPosition[0], laneA), )
        if self.debug:
            print("A", choiceA, "B", choiceB)

        if not areAdjacentLanes:
            if choiceA < choiceB:
                return choiceA, 0, True
            else:
                return choiceB, 0, False


        if choiceA < choiceB:
            return choiceB, choiceB - choiceA, True
        else:
            overlapA = format(highwayA, "b")[::-1][choiceB:choiceA]
            overlapB = format(highwayB, "b")[::-1][(choiceB+dColumn):(choiceA+dColumn)]
            if self.debug:
                print("overlaps", overlapA, overlapB)
            if overlapA.count("0") >= overlapB.count("0"):
                if self.debug:
                    print("choosing choice A", choiceA, "over B", choiceB)
                return choiceA, overlapA.count("1"), False
            else:
                if self.debug:
                    print("choosing choice B", choiceB, "over A", choiceA)
                return choiceB, overlapB.count("1"), False



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
        if l_ < l:
            for _ in range(leapingLanes):
                if self.j >= len(self.hurdleBits.dna2.string):
                    exit(1)
                self.match["dna1"] += "-"
                self.match["dna2"] += self.hurdleBits.dna2.string[self.j]
                self.leapCost += 1
                self.j += 1
                
        else:
            for _ in range(leapingLanes):
                if self.i >= len(self.hurdleBits.dna1.string):
                    exit(1)
                self.match["dna2"] += "-"
                self.match["dna1"] += self.hurdleBits.dna1.string[self.i]
                self.leapCost += 1
                self.i += 1
                
        
        for _ in range(length):
            if self.j >= len(self.hurdleBits.dna2.string):
                exit(1)
            if self.i >= len(self.hurdleBits.dna1.string):
                exit(1)
            self.match["dna1"] += self.hurdleBits.dna1.string[self.i]
            self.match["dna2"] += self.hurdleBits.dna2.string[self.j]
            self.hurdleCost += (self.hurdleBits.dna1.string[self.i] != self.hurdleBits.dna2.string[self.j])
            self.j += 1
            self.i += 1
            
            
        if self.debug:
            print(self.match["dna1"])
            print(self.match["dna2"])



    def editDistance(self):

        # Store leap and hurdle costs
        leapCost = 0
        hurdleCost = 0

        # Find the first highway
        lane, highway, length = self.findBestHighwayNearby(self.currentPosition)
        if self.debug:
            print(lane, highway,  length)

        #highwayLength = self.hurdleBits.getFirstHighwayLength(self.currentPosition[0], self.currentPosition[1])
        leapCol = leapForwardColumn(self.currentPosition[0], lane)
        
        #FIXME
        # Update position to the end of the highway
        colAfterLeap = self.length - gmpy.scan0(int(format(highway, "b")[::-1], 2))
        #self.length - gmpy.scan0(int(format(highway, "b")[::-1], 2))#self.currentPosition[1] + leapCol + length

        while colAfterLeap < self.length - 1:
            if self.debug:
                print("new position,", lane, colAfterLeap)
            # find the next highway
            lane_, highway_, length_ = self.findBestHighwayNearby((lane, colAfterLeap))
            
            if abs(lane - lane_) <= 1:
                # In adjacent lanes
                leapColumn, hurdles, _ = self.decideWhereToLeap(highway, highway_, lane, lane_)
                if self.debug:
                    print("leaping at column", leapColumn)
                self._updateMatch(self.currentPosition[0], lane, leapColumn - self.currentPosition[1])
                self.currentPosition = (lane_, leapColumn + leapForwardColumn(lane, lane_))
                if self.debug:
                    print("current position", self.currentPosition)
                self._updateMatch(lane, lane_, 0)
            else:
                leapColumn, hurdles, hasOverlap = self.decideWhereToLeap(highway, highway_, lane, lane_, areAdjacentLanes=False)
                if self.debug:
                    print("leaping at column", leapColumn, "has overlap", hasOverlap)
                if hasOverlap:
                    self._updateMatch(self.currentPosition[0], lane, leapColumn - self.currentPosition[1])
                    self.currentPosition = (lane_, leapColumn + leapForwardColumn(lane, lane_))
                    if self.debug:
                        print("current position", self.currentPosition)
                    self._updateMatch(lane, lane_, 0)
                else:
                    #self._updateMatch(self.currentPosition[0], lane, leapColumn - self.currentPosition[1])
                    #self.currentPosition = (lane, leapColumn)
                    if self.debug:
                        print("current position", self.currentPosition)

                    highwayReversed_ = int(format(highway_, "b")[::-1], 2)
                    highwayStartCol_ = self.length - (gmpy.scan1(highwayReversed_ >> gmpy.scan0(highwayReversed_)) + gmpy.scan0(highwayReversed_))
                    if self.debug:
                        print("starting col", highwayStartCol_, gmpy.scan0(highway_))
                    while self.currentPosition[1] + leapForwardColumn(self.currentPosition[0], lane_) < highwayStartCol_:
                        shortHighwayLane, shortHighway = self.findBestShortHighwayNearby(self.currentPosition, lane_, highwayStartCol_)
                    
                        if shortHighwayLane is None:
                            break
                    
                        else:
                            if self.debug:
                                print("chosen small highway", shortHighwayLane, shortHighway)
                            self._updateMatch(self.currentPosition[0], shortHighwayLane, shortHighway[1] - (self.currentPosition[1] + leapForwardColumn(self.currentPosition[0], shortHighwayLane)))
                            self.currentPosition = (shortHighwayLane, shortHighway[1])
                        
                        if self.debug:
                            print("current position", self.currentPosition)

                    if self.currentPosition[1] + leapForwardColumn(self.currentPosition[0], lane_) < gmpy.scan0(highway_) or self.currentPosition[0] != lane_:
                        self._updateMatch(self.currentPosition[0], lane_, max(gmpy.scan0(highway_) - leapForwardColumn(self.currentPosition[1], lane_) - self.currentPosition[1], 0))
                        self.currentPosition = (lane_, max(gmpy.scan0(highway_), self.currentPosition[1] + leapForwardColumn(self.currentPosition[0], lane_)))
                        if self.debug:
                            print("current position", self.currentPosition)
                    
  

            if self.debug:
                print("hurdle cost", self.hurdleCost)
                print("leap cost", self.leapCost)
            
            
            lane, highway, length = lane_, highway_, length_
            colAfterLeap = self.length - gmpy.scan0(int(format(highway, "b")[::-1], 2))
                
                


        self._updateMatch(self.currentPosition[0], lane, min(colAfterLeap - self.currentPosition[1], self.length - 1 - self.currentPosition[1]))
        #print(colAfterLeap - self.currentPosition[1])
                # see if we can find smaller highways in between

        #leapCost += leapLanePenalty(self.currentPosition[0], lane)
        
        
        self.currentPosition = (lane, colAfterLeap)
        if lane != self.destinationLane or self.currentPosition[1] < self.length - 1:
            self._updateMatch(self.currentPosition[0], self.destinationLane, self.length - 1 - self.currentPosition[1] - leapForwardColumn(self.currentPosition[0], self.destinationLane))
        
        if self.debug:
            print("hurdle cost", self.hurdleCost)
            print("leap cost", self.leapCost)
        
        return self.hurdleCost + self.leapCost

            

        




    


            

if __name__ == "__main__":
    g = GASMAProjection("AGAGCTAAACATGGCCGCACATAAATCGTTTTGAGTTGAAACTTTACCGCTGCATCTATTTTTCTCCTAGAATTATACCGTACACAGCCGACGTTCCACC", 
              "AGAGCTAAACAAGGGGCCCACATTAACGTTTTGAGCTTGAAGATCTTTACCGCGATCTATTTTTTCTCCTAGATTACCGTACACACCGACACTTCCATC", k=3, sight=3, debug=True)
    g.editDistance()

