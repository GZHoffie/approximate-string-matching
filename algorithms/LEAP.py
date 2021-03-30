from util import ApproximateStringMatching
import numpy as np

class LEAP(ApproximateStringMatching):
    def __init__(self, dna1, dna2, k, E, penalty=None, forward=None, 
                 originLanes=None, destinationLanes=None, hurdleCost=1):
        if len(dna1) > len(dna2):
            # swap the two
            temp = dna1
            dna1 = dna2
            dna2 = temp
        
        super().__init__(dna1, dna2)
        self.k = k # Maximum differences allowed
        self.E = E # Maximum energy given
        
        #TODO use callable hurdleCost so that the cost is changeable 
        self.hurdleCost = hurdleCost # The cost of hurdles

        # Initialize LEAP
        self.start = np.zeros((2*k + 1, E + 1))
        self.start.fill(float('-inf'))
        self.end = np.zeros((2*k + 1, E + 1))
        self.end.fill(float('-inf'))
        self.hurdles = np.zeros((2*k + 1, self.m))
        self.initHurdleVectors()
        self.originLanes = {0: 0} if originLanes is None else originLanes 
        if destinationLanes is None:
            self.destinationLanes = {0: self.m}
        else:
            self.destinationLanes = destinationLanes

        # Define essential cost functions
        self.penalty = penalty
        self.forward = forward

        # Store result of edit distance calculation
        self.finalLane = None
        self.finalEnergy = None
    

    def leapLanePenalty(self, l_, l):
        """
        Returns the penalty of leaping from lane l_ to l. When l and l_ are
        the same lane, then the penalty is just the energy cost of next hurdle.
        """
        if self.penalty is not None:
            return self.penalty(l_, l)
        else: # Use default penalty
            if l_ == l:
                return self.hurdleCost
            else:
                return self.hurdleCost * abs(l_ - l)
        

    def leapForwardColumn(self, l_, l, pos=0):
        """
        Returns the number of columns the toad moves forward when leaping from lane l_ to l. 
        When l and l_ are the same lane, then the number should just be 1.
        """
        if self.penalty is not None:
            return self.penalty(l_, l)
        else: # Use default layout
            if l_ == l:
                return 1 if pos < self.m else 0
            elif abs(l_) > abs(l) and l * l_ >= 0:
                return 0
            elif abs(l_) < abs(l) and l * l_ >= 0:
                return abs(l - l_)
            else:
                return abs(l - l_) - abs(l_)




    def initHurdleVectors(self):
        """
        Detect hurdles in the swimming pool and encode it using bit vectors
        """
        for i in range(2 * self.k + 1):
            lane = i - self.k
            if lane <= 0:
                self.hurdles[i] = [self.match(x, x - lane) for x in range(lane + 1, self.m + 1 + lane)]
            else:
                self.hurdles[i] = [self.match(x, x - lane) for x in range(1, self.m + 1)]

        print(self.hurdles)

    
    def verticesToHurdle(self, lane, position):
        #TODO: Use De Brujin to implement this
        lane = int(lane)
        if position >= self.m - 1:
            return 0
        tempPos = int(position) if position >= 0 else 0

        print("lane", lane, "col", position)
        while self.hurdles[lane + self.k][tempPos + 1] != 0:
            print("swim!")
            tempPos += 1
            if tempPos >= self.m - 1:
                break
        print(tempPos)
        return tempPos - position

    
    def editDistance(self):
        # Initialization
        finalEnergy = float('inf')
        k = self.k
        
        for l in range(-k, k+1):
            if l in self.originLanes:
                self.start[l+k][0] = self.originLanes[l]
                length = self.verticesToHurdle(l, self.start[l+k][0])
                print("vth", length)
                self.end[l+k][0] = self.start[l+k][0] + length
        #print(self.start)
        #print(self.end)
        
        for e in range(1, self.E+1):
            for l in range(-k, k+1):

                for l_ in range(-k, k+1):
                    #print("l'", l_)
                    e_ = e - self.leapLanePenalty(l_, l)
                    #print("energy", e_)
                    if e_ >= 0:
                        candidateStart = self.end[l_+k][e_] + self.leapForwardColumn(l_, l, self.start[l_+k][e_])

                        if candidateStart > self.start[l+k][e]:
                            candidateStart = self.m if candidateStart > self.m else candidateStart
                            self.start[l+k][e] = candidateStart       
                print("start")
                print(self.start)

                length = self.verticesToHurdle(l, self.start[l+k][e])


                self.end[l+k][e] = self.start[l+k][e] + length
                print("end")
                print(self.end)
                #print("final", l, e, self.end[l+k][e])
                if l in self.destinationLanes and self.end[l+k][e] >= self.destinationLanes[l] - 1:
                    if e < finalEnergy:
                        self.finalLane = l
                        self.finalEnergy = e
                        return True
        
        return False


    def backtrack(self):
        path = []
        l = self.finalLane
        e = self.finalEnergy
        k = self.k
        pathCount = 1
        path.append({"lane": l,
                     "start": self.start[l + k][e],
                     "end": self.end[l + k][e]})
        while l not in self.originLanes or self.start[l+k][e] != self.originLanes[l]:

            for l_ in range(-k, k+1):
                e_ = e - self.leapLanePenalty(l_, l)
                print("current", l_, e_)

                if self.end[l_+k][e_] + self.leapForwardColumn(l_, l) == self.start[l+k][e]:
                    l = l_
                    e = e_
                    break

            path.insert(0, {"lane": l,
                            "start": self.start[l + k][e],
                            "end": self.end[l + k][e]})

            pathCount += 1

        return path, pathCount

        
if __name__ == "__main__":
    prob = LEAP("ACTAGAACTT", "ACTTAGCACT", 2, 10)
    prob.editDistance()
    print(prob.finalLane, prob.finalEnergy)
    path, pathCount = prob.backtrack()
    print(path, pathCount)


        
        



