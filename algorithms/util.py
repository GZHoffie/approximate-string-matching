import numpy as np
from multiprocessing import Pool
import re

class DNA:
    def __init__(self, string=None):
        self.string = string
        self.bits = self.bitEncoding()

    def bitEncoding(self):
        encoder = {'A': "00", 'C': "01", 'G': "10", 'T': "11"}
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


class HurdleMatrix():
    def __init__(self, dna1, dna2, k, mismatchCost=None, insertCost=None, deleteCost=None):
        self.dna1 = DNA(dna1)
        self.dna2 = DNA(dna2)
        self.k = k
        self.m = len(dna1) # Length of dna1
        self.n = len(dna2) # Length of dna2
        self.hurdleMatrix = self.calculateHurdleMatrix()
        self.ignoredHurdles = self.preprocessHurdleMatrix()
        print(self.hurdleMatrix)
    
    def _match_str(self, i, j):
        if not 1 <= i <= self.m or not 1 <= j <= self.n:
            return "1"
        return "0" if self.dna1.string[int(i-1)] == self.dna2.string[int(j-1)] else "1"
    
    def _get_hurdles(self, shift):
        if shift <= 0:
            hurdles = [self._match_str(x, x - shift) for x in range(shift + 1, self.m + 1 + shift)]
        else:
            hurdles = [self._match_str(x, x - shift) for x in range(1, self.m + 1)]
        return "".join(hurdles)
    
    def _remove_single_zeros(self, string):
        occurences = []
        index = [m.start() for m in re.finditer('101', string)]
        while len(index) > 0:
            occurences += index
            string = re.sub('101', '111', string)
            index = [m.start() for m in re.finditer('101', string)]
        
        #self.hurdleMatrix[shift + self.k] = string
        return string, occurences
    
    def _remove_single_ones(self, string):
        occurences = []
        index = [m.start() for m in re.finditer('010', string)]
        while len(index) > 0:
            occurences += index
            string = re.sub('010', '000', string)
            index = [m.start() for m in re.finditer('010', string)]
        
        #self.hurdleMatrix[shift + self.k] = string
        return string, occurences
    
    def _preprocess(self, shift):
        string, ones = self._remove_single_ones(self.hurdleMatrix[shift + self.k])
        string, _ = self._remove_single_zeros(string)
        return string, ones, shift + self.k


    def calculateHurdleMatrix(self):
        with Pool() as p:
            return p.map(self._get_hurdles, list(range(-self.k, self.k+1)))
    
    def preprocessHurdleMatrix(self):
        ones = []
        zeros = []

        with Pool() as p:
        #    ones = p.map(self._remove_single_ones, list(range(-self.k, self.k+1)))
            for string, occur, i in p.map(self._preprocess, list(range(-self.k, self.k+1))):
                self.hurdleMatrix[i] = string
                ones.append(occur)
        
        return ones

        #for i in range(-self.k, self.k+1):
        #    self._remove_single_zeros(i)
        #return ones, zeros

    def find


if __name__ == "__main__":
    hd = HurdleMatrix("ACTAGAACTT", "ACTTAGCACT", 2)

