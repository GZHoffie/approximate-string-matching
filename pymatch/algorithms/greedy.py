from util import ApproximateStringMatching

def GASMA(ApproximateStringMatching):
    def __init__(self, dna1, dna2, k, E, leapPenalty=None, forward=None, 
                 originLanes=None, destinationLanes=None, hurdleCost=1):
        if len(dna1) > len(dna2):
            # swap the two
            temp = dna1
            dna1 = dna2
            dna2 = temp
        
        super().__init__(dna1, dna2)
        self.leapPenalty = lambda x: 2 * x if leapPenalty is None else leapPenalty

    
    