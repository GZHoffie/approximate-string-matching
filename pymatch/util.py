import numpy as np
from multiprocessing import Pool
import re
import time
import gmpy

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
    def __init__(self, dna1, dna2, k, mismatchCost=None, leapCost=None, threshold=3, crossHurdleThreshold=1, debug=False, reverse=False):
        self.dna1 = DNA(dna1)
        self.dna2 = DNA(dna2)
        self.k = k
        self.m = len(dna1) # Length of dna1
        self.n = len(dna2) # Length of dna2
        self.lookUpTable = deBrujin32Bit()
        self.debug = debug
        #curr = time.time()
        #print(self.calculateHurdleMatrix())
        self.hurdleMatrix = [int(s, 2) for s in self.calculateHurdleMatrix(reverse=reverse)]
        #print("Find hurdle matrix time:", time.time() - curr)
        #ignoredHurdles = self.preprocessHurdleMatrix()

        """
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
        print(self.highways)
        """
        
    
    def _match_str(self, i, j):
        if not 1 <= i <= self.m or not 1 <= j <= self.n:
            return "1"
        return "0" if self.dna1.string[int(i-1)] == self.dna2.string[int(j-1)] else "1"
    
    def _get_hurdles(self, shift, reverse):
        if shift <= 0:
            if reverse:
                hurdles = [self._match_str(x, x - shift) for x in reversed(range(shift, self.n + 1 + shift))]
            else:
                hurdles = [self._match_str(x, x - shift) for x in range(shift, self.n + 1 + shift)]
        else:
            if reverse:
                hurdles = [self._match_str(x, x - shift) for x in reversed(range(0, self.n + 1))]
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


    def calculateHurdleMatrix(self, reverse=False):
        """
        with Pool() as p:
            return p.map(self._get_hurdles, list(range(-self.k, self.k+1)))
        """
        return [self._get_hurdles(i, reverse) for i in range(-self.k, self.k+1)]
    
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
        
        temp_highway = (shift, 0, 0)
        highway_transformed = []
        hurdles_in_highway = []
        for h in highways:
            #print(h, temp_highway)
            if temp_highway is None:
                temp_highway = h
            elif h[1] - (temp_highway[1] + temp_highway[2]) <= self.crossHurdleThreshold and h[2] >= 2:
                #print(shift, h[1], temp_highway[1] + temp_highway[2])
                hurdles_in_highway += list(range(temp_highway[1] + temp_highway[2], h[1]))
                #print(hurdles_in_highway)
                temp_highway = (temp_highway[0], temp_highway[1], temp_highway[2] + h[1] - (temp_highway[1] + temp_highway[2]) + h[2])
            else:
                highway_transformed += [(temp_highway[0],
                                         temp_highway[1],
                                         temp_highway[2],
                                         hurdles_in_highway)]
                temp_highway = h
                hurdles_in_highway = []
        
        if temp_highway != (shift, 0, 0):
            highway_transformed += [(temp_highway[0],
                                     temp_highway[1],
                                     temp_highway[2],
                                     hurdles_in_highway)]
        
        return highway_transformed

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
            highways += self._get_highway(shift)
            #print(shift)

        if self.debug:
            print("highways:", highways)
        return highways



class HurdleBits:
    """
    Utility for GASMAProjection. Records the position of hurdles and 
    """
    def __init__(self, dna1: DNA, dna2: DNA, k: int, maxZerosIgnored=1, maxOnesIgnored=1, reverse=True, appendices=True) -> None:
        if len(dna1) > len(dna2):
            # swap the two
            temp = dna1
            dna1 = dna2
            dna2 = temp
        
        if appendices:
            candidate = ["A", "C", "G", "T"]
            candidate.pop(candidate.index(dna1[0]))
            if dna2[0] != dna1[0]:
                candidate.pop(candidate.index(dna2[0]))
            appendix = candidate[0] * (maxZerosIgnored + 5)

            dna1 = appendix + dna1 + appendix
            dna2 = appendix + dna2 + appendix
        
        self.dna1 = DNA(dna1)
        self.dna2 = DNA(dna2)
        self.k = k
        self.m = len(dna1) # Length of dna1
        self.n = len(dna2) # Length of dna2
        self.maxZeroIgnored = maxZerosIgnored
        self.maxOnesIgnored = maxOnesIgnored

        self.hurdles = self.calculateHurdleMatrix(reverse=reverse)
        self.bits = [int(string, 2) for string in self.hurdles]
        print("hurdles")
        for l in self.hurdles:
            print(l)
        print("After removing ones")
        self.processedBits = self.removeSingleOnes(self.bits)
        self.processedBits = self.removeSingleZeros(self.processedBits)
        #self.shiftRight(1)
        self.reversedBits = [int(format(s, "b")[::-1], 2) for s in self.bits]
        self.reversedProcessedBits = [int(format(s, "b")[::-1], 2) for s in self.processedBits]
        self.length = len(format(self.bits[0], "b"))
        #self.shiftRight(1)
        for l in self.processedBits:
            print(format(l, 'b'))
        print("Reversed")
        for l in self.reversedBits:
            print(format(l, 'b'))
        
        print(format(self.extractHighway(0, 0), "b"))

    def extractHighway(self, shift, col):
        """
        extract the highway at lane `shift` and column `col`.
        """
        leftBound = gmpy.scan1(self.processedBits[shift + self.k] >> col)
        rightBound = gmpy.scan1(self.reversedProcessedBits[shift + self.k] >> (self.length - col))
        #print(leftBound, rightBound)
        if leftBound < 0 or rightBound < 0:
            return int("1" * self.length, 2)
        reconstructHighway = "1" * (self.length - col - leftBound) + "0" * (leftBound + rightBound) + "1" * (col - rightBound)
        #reconstructHighway = int(reconstructHighway, 2) | self.bits[shift + self.k]
        highway = int(reconstructHighway, 2) | self.bits[shift + self.k]
        return highway


    
    def _match_str(self, i, j):
        if not 1 <= i <= self.m or not 1 <= j <= self.n:
            return "1"
        return "0" if self.dna1.string[int(i-1)] == self.dna2.string[int(j-1)] else "1"
    
    def _get_hurdles(self, shift, reverse):
        if shift <= 0:
            if reverse:
                hurdles = [self._match_str(x, x - shift) for x in reversed(range(shift, self.n + 2 + shift))]
            else:
                hurdles = [self._match_str(x, x - shift) for x in range(shift, self.n + 2 + shift)]
        else:
            if reverse:
                hurdles = [self._match_str(x, x - shift) for x in reversed(range(0, self.n + 2))]
            else:
                hurdles = [self._match_str(x, x - shift) for x in range(0, self.n + 2)]
        
        return "".join(hurdles)
    

    def calculateHurdleMatrix(self, reverse=False):
        """
        with Pool() as p:
            return p.map(self._get_hurdles, list(range(-self.k, self.k+1)))
        """
        return [self._get_hurdles(i, reverse) for i in range(-self.k, self.k+1)]
    
    def removeSingleOnes(self, bits):
        """
        Remove the single streaks of ones that has length smaller or equal to self.maxOnesIgnored
        in self.hurdles.
        """
        #mark = 0
        bitsProcessed = bits.copy()
        for l in range(len(bitsProcessed)):
            mark = -1
            for i in range(len(self.hurdles[0])):
                #print(i)
                if (bitsProcessed[l] >> i) & 1 == 1:
                    if (i == 0 or (bitsProcessed[l] >> i-1) & 1 == 0):
                        mark = i
                elif mark > 0 and i - mark <= self.maxOnesIgnored:
                    #print(i, mark)
                    for j in range(mark, i):
                        #print(format(1 << j, 'b'))
                        bitsProcessed[l] = bitsProcessed[l] ^ (1 << j)
        
        return bitsProcessed
    
    def removeSingleZeros(self, bits):
        """
        Remove the single streaks of zeros that has length smaller or equal to self.maxZerosIgnored
        in self.hurdles.
        """
        #mark = 0
        bitsProcessed = bits.copy()
        for l in range(len(bitsProcessed)):
            mark = -1
            for i in range(len(self.hurdles[0])):
                #print(i)
                if (bitsProcessed[l] >> i) & 1 == 0:
                    if (i == 0 or (bitsProcessed[l] >> i-1) & 1 == 1):
                        mark = i
                elif mark >= 0 and i - mark <= self.maxOnesIgnored:
                    #print(i, mark)
                    for j in range(mark, i):
                        #print(format(1 << j, 'b'))
                        bitsProcessed[l] = bitsProcessed[l] ^ (1 << j)
        
        return bitsProcessed
    
    def shiftRight(self, num):
        """
        shift all self.bits to the right by `num` 
        """
        self.bits = [b >> num for b in self.bits]
        self.processedBits = [b >> num for b in self.processedBits]
    
    def findFirstHighway(self, l, col):
        """
        Return the starting point of first highway in lane l and column `col`.
        """
        l = l + self.k
        #print("looking at", col, format(self.processedBits[l] >> col, "b"))
        return gmpy.scan0(self.processedBits[l] >> (col)) + col
    
    def getFirstHighwayLength(self, l, col):
        """
        Return the length of first highway in lane l.
        """
        #print(format(self.processedBits[l + self.k], "b"))
        #print("highway", format(self.processedBits[l + self.k] >> self.findFirstHighway(l, col), "b"))
        return gmpy.scan1(self.processedBits[l + self.k] >> self.findFirstHighway(l, col))

                

    

        

if __name__ == "__main__":
    import time
    a = time.time()
    hd = HurdleBits("AGAGCTAAACATGGCCGCACATAAATCGTTTTGAGTTGAAACTTTACCGCTGCATCTATTTTTCTCCTAGAATTATACCGTACACAGCCGACGTTCCACC", 
              "AGAGCTAAACAAGGGGCCCACATTAACGTTTTGAGCTTGAAGATCTTTACCGCGATCTATTTTTTCTCCTAGATTACCGTACACACCGACACTTCCATC", k=2,  reverse=True)
    #print(hd.hurdleMatrix)
    print("calculate hurdle time:", time.time() - a)

