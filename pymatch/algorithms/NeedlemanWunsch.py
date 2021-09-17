from pymatch.util import ApproximateStringMatching
import numpy as np

class NeedlemanWunsch(ApproximateStringMatching):
    def __init__(self, dna1, dna2, mismatchCost=None, insertCost=None, deleteCost=None):
        super().__init__(dna1, dna2, mismatchCost, insertCost, deleteCost)
        self.D = np.zeros((self.m + 1, self.n + 1)) # edit distance matrix
        self.alignment = {"dna1": "", "dna2": ""}


    def editDistance(self):
        for j in range(1, self.n + 1):
            self.D[0][j] = j
        for i in range(1, self.m + 1):
            self.D[i][0] = i
        for i in range(1, self.m + 1):
            for j in range(1, self.n + 1):
                delete = self.D[i-1][j] + self.deleteCost({})
                insert = self.D[i][j-1] + self.insertCost({})
                match = self.D[i-1][j-1] + self.mismatchCost({"i": i, "j": j})
                self.D[i][j] = min(insert, delete, match)
        
        return self.D[self.m][self.n]
    
    def SGEditDistance(self):
        """
        Semiglobal: end when reaching the end of matrix
        """
        self.editDistance()
        return min(np.min(self.D[:, self.n]), np.min(self.D[self.m, :]))


    def backtrack(self):
        i = self.m
        j = self.n
        while i > 0 or j > 0:
            if i > 0 and j > 0 and \
                self.D[i][j] == self.D[i-1][j-1] + self.mismatchCost({"i": i, "j": j}):
                # Match with same character
                self.alignment["dna1"] = self.dna1.string[i-1] + self.alignment["dna1"]
                self.alignment["dna2"] = self.dna2.string[j-1] + self.alignment["dna2"]
                i -= 1
                j -= 1
            elif i > 0 and self.D[i][j] == self.D[i-1][j] + self.deleteCost({}):
                self.alignment["dna1"] = self.dna1.string[i-1] + self.alignment["dna1"]
                self.alignment["dna2"] = "-" + self.alignment["dna2"]
                i -= 1
            else:
                self.alignment["dna1"] = "-" + self.alignment["dna1"]
                self.alignment["dna2"] = self.dna2.string[j-1] + self.alignment["dna2"]
                j -= 1
            #print(self.alignment)


    def __str__(self):
        self.backtrack()
        return self.alignment["dna1"] + "\n" + self.alignment["dna2"]


if __name__ == "__main__":
    prob = NeedlemanWunsch("CGCACAGGTCACCGAGAGTGCATTAGATTGTCACCGGTCGCCTTAACAGCGTGGTTATTACGGGAACTTTCCAAAAAGATTTGGGGCTCAATTGGAGATG", 
              "CGGCACAGGTCCGCGCGAGAAGTCATTAGATTGTCACCGGTCGCCCTAACAGGTGGTTATAGGGATTCCAAAAAGGATTTGGGGCGTTCAATTGAAGATG")
    print(prob.editDistance())
    print(prob.D)
    print(prob)


