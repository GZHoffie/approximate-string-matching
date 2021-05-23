from pymatch.util import ApproximateStringMatching
import numpy as np

class LandauVishkin(ApproximateStringMatching):
    def __init__(self, dna1, dna2, k, mismatchCost=None, insertCost=None, deleteCost=None):
        if len(dna1) > len(dna2):
            # swap the two
            temp = dna1
            dna1 = dna2
            dna2 = temp

        super().__init__(dna1, dna2, mismatchCost, insertCost, deleteCost)
        self.k = k # The maximum number of differences allowed
        self.L = np.zeros((2*k+3, k+3)) # Stores the farthest row reachable with limited differences


    def getL(self, diag, diff):
        """
        get from L the corresponding value, with `diag` being the diagnonal
        and `diff` being the number of differnces. The number represents the 
        maximum row i such that D[i][i+d] = e.
        """
        assert diff <= self.k, "Number of difference too large."
        assert -(self.k+1) <= diag <= self.k+1, "Diagonal \
            is not considered under {} differences.".format(self.k)
        return self.L[diag + self.k + 1][diff + 2]


    def setL(self, diag, diff, value):
        """
        Set the value in L for specified diagonal and number of differences.
        """
        assert diff <= self.k, "Number of difference too large."
        assert -(self.k+1) <= diag <= self.k+1, "Diagonal \
            is not considered under {} differences.".format(self.k)
        self.L[diag + self.k + 1][diff + 2] = value


    def editDistance(self):
        k = self.k
        for d in range(-k-1, k+2):
            self.setL(d, abs(d)-2, float('-inf'))
            if d < 0:
                self.setL(d, abs(d)-1, abs(d)-1)  # The starting row is |d|-1-th row
            else:
                self.setL(d, abs(d)-1, -1)        # The starting row is the first row
        for e in range(k+1):
            for d in range(-e, e+1):
                row = max(self.getL(d, e-1) + 1,  # The new difference is a mismatch, row += 1
                          self.getL(d-1, e-1),    # The new difference is an insertion, row unchanged
                          self.getL(d+1, e-1) + 1 # the new difference is a deletion, row += 1
                          )
                while self.match(row+1, row+1+d):
                    row += 1
                self.setL(d, e, row)
                if self.getL(d, e) == self.m:
                    return d, e # Successful match
        return None, None # Unsuccessful match

    
    def backtrack(self):
        #TODO implement the backtracking of the algorithm
        raise NotImplementedError

    
    def __str__(self):
        return str(self.L)


if __name__ == "__main__":
    prob = LandauVishkin("AGCGCTTGCTGC", "AGTCGCCGCTGCTGC", 3)
    d, e = prob.editDistance()
    if d is not None:
        print("Successful match with", e, "differences on diagonal", d)
    else:
        print("Unsuccessful match within", prob.k, "differences")
    print(prob)
