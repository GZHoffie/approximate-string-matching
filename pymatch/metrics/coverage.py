
def covers(str1, str2):
    """
    Return if str1 covers str2.
    """
    n = len(str1)
    m = len(str2)
    
    if n < m:
        return 0
    
    i = 0
    for j in range(m):
        if i >= n:
            return 0
        while str1[i] != str2[j]:
            print(str1[i], i, str2[j], j)
            i += 1
            if i >= n:
                return 0
        i += 1
        
    
    return 1


class Coverage:
    def __init__(self, dict1, dict2, shortestMatchAccepted=3) -> None:
        """
        A metric telling us whether the matching in dict1 'covers' the matching in dict2.
        A match `A` 'covers' another match `B` if the string formed by long matching substrings
        in `A` 'contains' that in `B`.

        Args:
            dict1, dict2 (Dict): Dicts containing two fields, "dna1" and "dna2", indicating how
                the two dnas are matched together.
        """
        self.shortestMatchAccepted = shortestMatchAccepted
        self.match1 = self.findLongMatchingSubstring(dict1)
        self.match2 = self.findLongMatchingSubstring(dict2, removeSmallMatches=True)
        print(self.match1)
        print(self.match2)
    
    def findLongMatchingSubstring(self, dnaDict, removeSmallMatches=False):
        """
        Find the long matching substring in dnaDict.
        """
        matchBits = [str(int(dnaDict["dna1"][i] != dnaDict["dna2"][i])) for i in range(len(dnaDict["dna1"]))]
        #print(matchBits)
        if removeSmallMatches:
            matchBits = self.removeSmallMatches(matchBits)
        #print(matchBits)

        return "".join([i for i, j in zip(dnaDict["dna1"], matchBits) if j == "0"])
    
    def removeSmallMatches(self, bits):
        """
        Remove the single streaks of zeros that has length smaller or equal to self.maxZerosIgnored
        in self.hurdles.
        """
        #mark = 0
        
        mark = -1
        for i in range(len(bits)):
            #print(i)
            if bits[i] == '0':
                if (i == 0 or bits[i-1] == '1'):
                    mark = i
                    #print(mark)
            elif mark >= 0 and i - mark <= self.shortestMatchAccepted:
                #print(i, mark)
                for j in range(mark, i):
                    #print(format(1 << j, 'b'))
                    bits[j] = '1'
                mark = i
        
        return bits
    
    def compute(self):
        return int(covers(self.match1, self.match2))


if __name__ == "__main__":
    dict1 = {"dna1": "AGAGCTAAACATGG-CCGCACATAAATCGTTTTGAG-TTGAA-A-CTTTACCGCTGCATCTATTTTT-CTCCTAGAATTATACCGTACACAGCCGAC-GTTCCACC",
             "dna2": "AGAGCTAAACAAGGGGCCCACATTAA-CGTTTTGAGCTTGAAGATCTTTACCGC-G-ATCTATTTTTTCTCCTAGA-TTA--CCGTACACA-CCGACACTTCCATC"}
    
    dict2 = {"dna1": "AGAGCTAAAC-ATGGCCGCACATAAATCGTTTTGAG-TTGAA-A-CTTTACCGCTGCATCTA-TTTTTCTCCTAGAATTATACCGTACACAGCCGAC-GTTCCACC",
             "dna2": "AGAGCTAAACAAGGGGCCCACATTAA-CGTTTTGAGCTTGAAGATCTTTACCGC-G-ATCTATTTTTTCTCCTAG-A-T-TACCGTACACA-CCGACACTTCCATC"}
    
    c = Coverage(dict1, dict2, 3)
    print(c.compute())

    
    
