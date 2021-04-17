import numpy as np

class DNA:
    def __init__(self, string=None):
        self.string = string
        self.bits = self.bitEncoding()

    def bitEncoding(self):
        encoder = {'A': [0, 0], 'C': [0, 1], 'G': [1, 0], 'T': [1, 1]}
        return [encoder[char] for char in self.string]


class ApproximateStringMatching:
    def __init__(self, dna1, dna2, mismatchCost=None, insertCost=None, deleteCost=None):
        self.dna1 = DNA(dna1)
        self.dna2 = DNA(dna2)
        self.m = len(dna1) # Length of dna1
        self.n = len(dna2) # Length of dna2
        self.alignment = np.zeros((self.m, self.n))
        self.mismatchCost = mismatchCost
        self.insertCost = insertCost
        self.deleteCost = deleteCost
        self.initCostFunctions()


    def initCostFunctions(self):
        """
        initialize cost functions
        """
        if self.mismatchCost is None:
            # Use default mismatching cost
            self.mismatchCost = lambda config: 0 if self.match(config["i"], config["j"]) else 1
        if self.insertCost is None:
            # Use default insertion cost
            self.insertCost = lambda config: 1
        if self.deleteCost is None:
            # Use default deletion cost
            self.deleteCost = lambda config: 1

    
    def match(self, i, j):
        """
        check whether dna1[i] is the same as dna2[j]
        """
        if not 1 <= i <= self.m or not 1 <= j <= self.n:
            return None
        return self.dna1.string[int(i-1)] == self.dna2.string[int(j-1)]

def deBrujin32Bit():
    return [0,  1, 16,  2, 29, 17,  3, 22,
            30, 20, 18, 11, 13,  4,  7, 23,
            31, 15, 28, 21, 19, 10, 12,  6,
            14, 27,  9,  5, 26,  8, 25, 24]

def bit_not(n, numbits=32):
    return (1 << numbits) - 1 - n