from pymatch.algorithms import GASMA

class GASMAShortsighted(GASMA):
    def __init__(self, dna1, dna2, k, leapCost=None, hurdleCost=1, threshold=3, sight=3, debug=False):

        appendix = "A" * (threshold + 1)
        dna1 = appendix + dna1 + appendix
        dna2 = appendix + dna2 + appendix

        super().__init__(dna1, dna2, k, leapCost, hurdleCost, threshold, debug)
        self.highways = [(shift, start + length - 1, length) for shift, start, length in self.highways]
        self.highways = sorted(self.highways, key=lambda x: x[1], reverse=True)
        if self.debug:
            print("transformed highways:",self.highways)
        self.sight = sight
        #print(self.highways)
    
    def editDistance(self):
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

        
        currentPosition = (0, self.n - 1)
        route = [currentPosition]
        hurdleCost = 0
        leapCost = 0
        while len(self.highways) > 0:
            bestHighway = 0
            bestHighwayLength = 0
            for i in range(len(self.highways)):
                if currentPosition[1] - self.highways[i][1] > self.sight:
                    break
                elif self.highways[i][2] > bestHighwayLength or (self.highways[i][2] == bestHighwayLength and self.highways[i][0] == self.destinationLane):
                    bestHighway = i
                    bestHighwayLength = self.highways[i][2]
            
            if bestHighway is None:
                print("No highway is in sight.")
            else:
                chosenHighway = self.highways.pop(bestHighway)
                if self.debug:
                    print("chosen", chosenHighway)
                leapCost += leapLanePenalty(currentPosition[0], chosenHighway[0])
                currentColumn = currentPosition[1] - leapForwardColumn(currentPosition[0], chosenHighway[0])
                if currentColumn - 1 > chosenHighway[1]:
                    #print(currentColumn, chosenHighway[1])
                    hurdleCost += currentColumn - 1 - chosenHighway[1]
                
                route.append((chosenHighway[0], chosenHighway[1]))
                currentPosition = (chosenHighway[0], chosenHighway[1] - chosenHighway[2] + 1)
                route.append(currentPosition)

            while len(self.highways) > 0:
                if self.highways[0][1] - self.highways[0][2] + 1 > currentPosition[1] - leapForwardColumn(currentPosition[0], self.highways[0][0]):
                    h = self.highways.pop(0)
                    #print("removed", h)
                else:
                    break
        
        if currentPosition != (self.destinationLane, 0):
            leapCost += leapLanePenalty(currentPosition[0], self.destinationLane)
            columnAfterLeap = currentPosition[1] - leapForwardColumn(currentPosition[0], self.destinationLane)
            if columnAfterLeap > 0:
                hurdleCost += 1 * columnAfterLeap
            route += [(self.destinationLane, 0)]
        
        #print("leap cost:", leapCost)
        #print("hurdle cost:", hurdleCost)
        return leapCost + hurdleCost, route
        




    

if __name__ == "__main__":
    g = GASMAShortsighted("CGATACCGCAGCCTATCAATAATCGAAAACGGTTCACTGGACACACAGCCCCCAGTCTTCAGAAGGAACGTCATTCCTTGCCTTGGTTGTCGGCGATTAC", 
              "CGATACCGCAGCCTATCAATAATCGAAAACGGTTCACTGGACACACAGCCCCCAGTCTTCAGAAGGAACGTCATTCCTTGCCTTGGTTGTCGGCGATTA", 2, 2, debug=True)
    cost, route = g.editDistance()
    print(cost)
    print(route)
    #print(g.editDistance())
        