import numpy as np
from multiprocessing import Pool
import re
import time

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


def deBrujin64Bit():
    return [0,  1, 48,  2, 57, 49, 28,  3,
            61, 58, 50, 42, 38, 29, 17,  4,
            62, 55, 59, 36, 53, 51, 43, 22,
            45, 39, 33, 30, 24, 18, 12,  5,
            63, 47, 56, 27, 60, 41, 37, 16,
            54, 35, 52, 21, 44, 32, 23, 11,
            46, 26, 40, 15, 34, 20, 31, 10,
            25, 14, 19,  9, 13,  8,  7,  6]

def bit_not(n, numbits=32):
    return (1 << numbits) - 1 - n


class HurdleMatrix():
    def __init__(self, dna1, dna2, k, mismatchCost=None, leapCost=None, threshold=3, crossHurdleThreshold=0, debug=False):
        self.dna1 = DNA(dna1)
        self.dna2 = DNA(dna2)
        self.k = k
        self.m = len(dna1) # Length of dna1
        self.n = len(dna2) # Length of dna2
        self.lookUpTable = deBrujin32Bit()
        self.debug = debug
        #curr = time.time()
        #print(self.calculateHurdleMatrix())
        self.hurdleMatrix = [int(s, 2) for s in self.calculateHurdleMatrix()]
        #print("Find hurdle matrix time:", time.time() - curr)
        #ignoredHurdles = self.preprocessHurdleMatrix()


        if mismatchCost is None:
            self.mismatchCost = 1
        else:
            assert isinstance(mismatchCost, int), "mismatchCost should be an int"
            self.mismatchCost = mismatchCost

        self.leapCost = leapCost
        self.threshold = threshold
        self.crossHurdleThreshold = crossHurdleThreshold
        

        #print(self.hurdleMatrix)
        #curr = time.time()
        self.highways = [h for h in self.getHighways() if h[2] >= threshold]
        #print("Find highway time:", time.time() - curr)
        #print(self.highways)
    
    def _match_str(self, i, j):
        if not 1 <= i <= self.m or not 1 <= j <= self.n:
            return "1"
        return "0" if self.dna1.string[int(i-1)] == self.dna2.string[int(j-1)] else "1"
    
    def _get_hurdles(self, shift):

        if shift <= 0:
            hurdles = [self._match_str(x, x - shift) for x in range(shift, self.n + 1 + shift)]
        else:
            hurdles = [self._match_str(x, x - shift) for x in range(0, self.n + 1)]
        
        if self.debug:
            print("".join(hurdles))
        return "".join(hurdles)
    
    def _remove_single_zeros(self, string):
        occurences = []
        index = [m.start() for m in re.finditer('101', string)]
        while len(index) > 0:
            occurences += index
            string = re.sub('101', '111', string)
            index = [m.start() + 1 for m in re.finditer('101', string)]
        
        #self.hurdleMatrix[shift + self.k] = string
        return string, occurences
    
    def _remove_single_ones(self, string):
        occurences = []
        index = [m.start() for m in re.finditer('010', string)]
        while len(index) > 0:
            occurences += index
            string = re.sub('010', '000', string)
            index = [m.start() + 1 for m in re.finditer('010', string)]
        
        #self.hurdleMatrix[shift + self.k] = string
        return string, occurences
    
    def _preprocess(self, shift):
        string = self.hurdleMatrix[shift + self.k]
        #string, ones = self._remove_single_ones(string)
        string, _ = self._remove_single_zeros(string)
        return string, shift + self.k


    def calculateHurdleMatrix(self):
        """
        with Pool() as p:
            return p.map(self._get_hurdles, list(range(-self.k, self.k+1)))
        """
        return [self._get_hurdles(i) for i in range(-self.k, self.k+1)]
    
    def preprocessHurdleMatrix(self, removeSharedHurdles=False):
        """
        Remove single zeros and ones, and convert the string to bits.
        """
        #ones = []
        #zeros = []

        #with Pool() as p:
        #    ones = p.map(self._remove_single_ones, list(range(-self.k, self.k+1)))
        #    for string, i in p.map(self._preprocess, list(range(-self.k, self.k+1))):
        #        self.hurdleMatrix[i] = int(string, 2)
                #ones.append(occur)
        
        #return ones

        #for i in range(-self.k, self.k+1):
        #    self._remove_single_zeros(i)
        #return ones, zeros
        mask = None
        if removeSharedHurdles:
            mask = int('1' * (self.m + 1), 2)
            for shift in range(-self.k, self.k+1):
                mask = mask & int(self.hurdleMatrix[shift + self.k], 2)
            
            mask = list(format(mask, 'b'))


        for shift in range(-self.k, self.k+1):
            string, i = self._preprocess(shift)
            if removeSharedHurdles:
                string = "".join([b for (b, m) in zip(string, mask) if m == '0'])
            self.hurdleMatrix[i] = int(string, 2)
    
    def _find_LSB(self, bits, findZero=True):
        if findZero:
            b_LSB = ~bits & (bits + 1)
        else:
            b_LSB = bits & (~bits + 1)
        b_LSB *= 0x6EB14F9
        b_LSB = b_LSB >> 27
        return self.lookUpTable[b_LSB % 32]
    
    def _get_highway(self, shift):
        #print(shift + self.k)
        bits = self.hurdleMatrix[shift + self.k]
        #print(bits)
        highways = []
        currentPos = 0
        tempPos = 0
        while bits > 0:
            LSZ = self._find_LSB(bits, findZero=True) # Find first zero
            #print(shift, format(bits,'b'))
            #print("LSZ", shift, LSZ)
            if LSZ is None:
                break
            currentPos += LSZ # Starting point of highway
            bits = bits >> LSZ
            if bits == 0:
                break
            nextLSB = self._find_LSB(bits, findZero=False) # Length of highway
            #print("nextLSB", shift, nextLSB)
            if nextLSB is None:
                break
            elif nextLSB == 0 and LSZ == 0:
                #tempPos += 32
                if not bits & 1:
                    highways.append((shift, currentPos, 32))
                bits = bits >> 32
                currentPos += 32

            else:
                highways.append((shift, currentPos, nextLSB))
                currentPos += nextLSB
                bits = bits >> nextLSB
            
            #print(shift, highways)
        
        return highways
    
    def _preprocess_highways(self, highways):
        processed_highways = []
        tmp = None
        #highwayHurdles = []
        for h in highways:
            if tmp is None:
                tmp = h
            elif h[1] - (tmp[1] + tmp[2]) <= self.crossHurdleThreshold:
                tmp = (tmp[0], tmp[1], tmp[2] + h[2])
            else:
                processed_highways += [tmp]
                tmp = h
            
            #print(h, tmp)
        
        if tmp is not None:
            processed_highways += [tmp]
        
        return processed_highways


    def getHighways(self):
        highways = []
        """
        with Pool() as p:
        #    ones = p.map(self._remove_single_ones, list(range(-self.k, self.k+1)))
            for h in p.map(self._get_highway, list(range(-self.k, self.k+1))):
                highways += h
        """
        for shift in range(-self.k, self.k+1):
            #print(shift, self._get_highway(shift))
            highways += self._preprocess_highways(self._get_highway(shift))
            #print(shift)

        if self.debug:
            print("highways:", highways)
        return highways


if __name__ == "__main__":
    import time
    a = time.time()
    hd = HurdleMatrix("GCAAAGTAACAACAGAGCCCCGGTATTGAGCGAGTAATCCACCCCACTCGATATGCGCTTACGATCCAAGCCTTGCTACTGAACGCAAATCCCCCTGCA", 
              "GCAAAGTAACAACAGAGCCCCGGTATTGAGCGAGTAATCCACCCCACTCGATATGCGCTTACGATCCAAGCCTTGCTACCTGAACGCAAATCCCCCTGCA",2, threshold=3)
    print(hd.highways)
    print("calculate hurdle time:", time.time() - a)

