from pymatch.algorithms import GASMA

class GASMAShortsighted(GASMA):
    def __init__(self, dna1, dna2, k, leapCost=None, hurdleCost=1, threshold=3, crossHurdleThreshold=1, sight=3, debug=False):

        super().__init__(dna1, dna2, k, leapCost, hurdleCost, threshold, crossHurdleThreshold, debug)
        self.debug = debug
        self.highways = [(shift, start + length, length, hurdles) for shift, start, length, hurdles in self.highways]
        self.highways = sorted(self.highways, key=lambda x: x[1], reverse=True)
        if self.debug:
            print("transformed highways:",self.highways)
        self.sight = sight
        self.match = {"dna1": "", "dna2": ""}
        self.i = 0
        self.j = 0
        #print('highways', self.highways)
    
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
        
        def _score(l, current_position):
            """
            Determine the score and the cost of a lane l from the current position.
            """
            columnAfterLeap = current_position[1] - leapForwardColumn(current_position[0], l[0])
            leapCost = leapLanePenalty(current_position[0], l[0])
            hurdleCost = len([hurdle for hurdle in l[3] if hurdle < columnAfterLeap])
            wayToHighwayCost = format(self.hurdleMatrix.hurdleMatrix[l[0] + self.k], 'b')[self.matrixLength-columnAfterLeap+1:self.matrixLength-l[1]].count('1')
            effectiveLength = min(columnAfterLeap - (l[1] - l[2]), l[2])#l[2] if columnAfterLeap >= l[1] else l[2] - (l[1] - columnAfterLeap + 1)
            score = effectiveLength - wayToHighwayCost - hurdleCost - leapCost

            if self.debug:
                print(l, score, "length", effectiveLength, "way to highway", wayToHighwayCost, "hurdle", hurdleCost, "leap", leapCost)
            return score, leapCost, hurdleCost + wayToHighwayCost, effectiveLength


        
        currentPosition = (0, self.n)
        route = [currentPosition]
        hurdleCost = 0
        leapCost = 0
        while len(self.highways) > 0:
            bestHighway = 0
            bestHighwayScore = float('-inf')
            _, lc, hc, bestLength = _score(self.highways[0], currentPosition)
            bestHighwayCost = (lc, hc)
            for i in range(len(self.highways)):
                score, this_leap_cost, this_hurdle_cost, length = _score(self.highways[i], currentPosition)
                if currentPosition[1] - leapForwardColumn(currentPosition[0], self.highways[i][0]) - self.highways[i][1] > self.sight and bestHighwayScore >= 0:
                    break
                elif score > bestHighwayScore or (score == bestHighwayScore and self.highways[i][0] == self.destinationLane):
                    bestHighway = i
                    bestHighwayScore = score
                    bestHighwayCost = (this_leap_cost, this_hurdle_cost)
                    bestLength = length
            
            chosenHighway = self.highways.pop(bestHighway)
            leapCost += bestHighwayCost[0]
            hurdleCost += bestHighwayCost[1]
                
            if self.debug:
                print("chosen", chosenHighway, "cost:", bestHighwayCost, "position:", currentPosition, "length:", bestLength)
                
                if chosenHighway[0] < currentPosition[0]:
                    self.match["dna1"] += '-' * abs(chosenHighway[0] - currentPosition[0]) 
                    for _ in range(abs(currentPosition[0] - chosenHighway[0])):
                        self.match["dna2"] += self.dna2.string[self.j] 
                        self.j += 1
                elif chosenHighway[0] > currentPosition[0]:
                    self.match["dna2"] += '-' * abs(currentPosition[0] - chosenHighway[0]) 
                    for _ in range(abs(currentPosition[0] - chosenHighway[0])):
                        self.match["dna1"] += self.dna1.string[self.i] 
                        self.i += 1
                for _ in range(bestLength):
                    try:
                        self.match["dna1"] += self.dna1.string[self.i]
                        self.match["dna2"] += self.dna2.string[self.j] 
                    except:
                        pass
                    
                    self.i += 1
                    self.j += 1
                    
                
                
                print(self.match["dna1"])
                print(self.match["dna2"])

                
                
            route.append((chosenHighway[0], chosenHighway[1]))
            currentPosition = (chosenHighway[0], chosenHighway[1] - chosenHighway[2] + 1)
            route.append(currentPosition)
            

            while len(self.highways) > 0:
                if self.highways[0][1] - self.highways[0][2] >= currentPosition[1] - leapForwardColumn(currentPosition[0], self.highways[0][0]):
                    h = self.highways.pop(0)
                    #print("removed", h)
                else:
                    break
        
        if currentPosition != (self.destinationLane, 0):
            leapCost += leapLanePenalty(currentPosition[0], self.destinationLane)
            columnAfterLeap = currentPosition[1] - leapForwardColumn(currentPosition[0], self.destinationLane)
            if columnAfterLeap > 0:
                hurdleCost += format(self.hurdleMatrix.hurdleMatrix[self.destinationLane + self.k], 'b')[self.matrixLength-columnAfterLeap+1:self.matrixLength-0].count('1')
            route += [(self.destinationLane, 0)]
        if self.debug:
            print("leap cost:", leapCost)
            print("hurdle cost:", hurdleCost)
            self.match["dna1"] += self.dna1.string[self.i:]
            self.match["dna2"] += self.dna2.string[self.j:]
            self.match["dna1"] = self.match["dna1"][(self.threshold + 7):-(self.threshold + 7)]
            self.match["dna2"] = self.match["dna2"][(self.threshold + 7):-(self.threshold + 7)]

        return leapCost + hurdleCost, route
        




    

if __name__ == "__main__":
    g = GASMAShortsighted("TAATGGTCGCCCTGCCCAAACTCCGAATTCATGCGATCCCTTTTCAAGCCTGACTTCATCCTATATATCCCACAAGCCCGGACTGATACGCTCGTGCTGG", 
              "TAATGGTCGCCCTGCCCAAACTCCGAATCTGCGATCCCTTTTCAAGCCTGACTTCATCCTATATATCCCACAAGCCCGGACTGATACGCTCGTGCGTGG", 7, threshold=2, crossHurdleThreshold=0, sight=2, debug=True)
    cost, route = g.editDistance()
    print(cost)
    print(route)
    print(g.match["dna2"])
    print(g.match["dna1"])
    print(g.destinationLane)
    #print(g.editDistance())
        