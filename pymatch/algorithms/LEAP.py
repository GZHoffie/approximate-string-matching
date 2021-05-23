from pymatch.util import ApproximateStringMatching, deBrujin32Bit, bit_not
import numpy as np

class LEAP(ApproximateStringMatching):
    def __init__(self, dna1, dna2, k, E, penalty=None, forward=None, 
                 originLanes=None, destinationLanes=None, hurdleCost=1):
        if len(dna1) > len(dna2):
            # swap the two
            temp = dna1
            dna1 = dna2
            dna2 = temp

        # Enable de-Brujin acceleration
        self.deBrujin = True
        
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
        if self.deBrujin:
            self.hurdles = []
            self.lookUpTable = deBrujin32Bit()
        else:
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
        if self.forward is not None:
            return self.forward(l_, l)
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
                hurdles = [self.match(x, x - lane) for x in range(lane + 1, self.m + 1 + lane)]
            else:
                hurdles = [self.match(x, x - lane) for x in range(1, self.m + 1)]

            if self.deBrujin:
                hurdlesBits = ['0' if x or x is None else '1' for x in reversed(hurdles)]
                hurdlesInt = int("".join(hurdlesBits), 2)
                self.hurdles.append(hurdlesInt)
            else:
                self.hurdles[i] = hurdles
        
        print(self.hurdles)
        

    
    def verticesToHurdle(self, lane, position):
        lane = int(lane)
        if position >= self.m - 1:
            return 0
        tempPos = int(position) if position >= 0 else 0

        if self.deBrujin:
            shiftBitVec = int(self.hurdles[lane + self.k]) >> tempPos
            b_LSB = shiftBitVec & (~shiftBitVec + 1)
            b_LSB *= 0x6EB14F9
            b_LSB = b_LSB >> 27
            return self.lookUpTable[b_LSB]

        else:
            while self.hurdles[lane + self.k][tempPos + 1] != 0:
                tempPos += 1
                if tempPos >= self.m - 1:
                    break
            return tempPos - position

    
    def editDistance(self):
        # Initialization
        finalEnergy = float('inf')
        k = self.k
        
        for l in range(-k, k+1):
            if l in self.originLanes:
                self.start[l+k][0] = self.originLanes[l]
                length = self.verticesToHurdle(l, self.start[l+k][0])
                self.end[l+k][0] = self.start[l+k][0] + length
        
        for e in range(1, self.E+1):
            for l in range(-k, k+1):

                for l_ in range(-k, k+1):
                    e_ = e - self.leapLanePenalty(l_, l)
                    if e_ >= 0:
                        candidateStart = self.end[l_+k][e_] + self.leapForwardColumn(l_, l, self.start[l_+k][e_])
                        if candidateStart > self.start[l+k][e]:
                            candidateStart = self.m if candidateStart > self.m else candidateStart
                            self.start[l+k][e] = candidateStart       

                length = self.verticesToHurdle(l, self.start[l+k][e])
                self.end[l+k][e] = self.start[l+k][e] + length
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


        
        



