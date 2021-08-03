from pymatch.util import ApproximateStringMatching, HurdleMatrix
from pymatch.algorithms import GASMA

class GASMAProjection(ApproximateStringMatching):
    def __init__(self, dna1, dna2, k, mismatchCost=None, insertCost=None, deleteCost=None):
        super().__init__(dna1, dna2, mismatchCost=mismatchCost, insertCost=insertCost, deleteCost=deleteCost)
        self.hurdleMatrix = HurdleMatrix(dna1=dna1, dna2=dna2, k=k, threshold=0,)
    

    def findHighways(self):
        """
        Find all the available highways given the current position.
        """
        for l in range(-self.k, self.k + 1):
            for i in range(self.sight):
                pass

if __name__ == "__main__":
    GASMAProjection("AGAGCTAAACATGGCCGCACATAAATCGTTTTGAGTTGAAACTTTACCGCTGCATCTATTTTTCTCCTAGAATTATACCGTACACAGCCGACGTTCCACC", 
              "AGAGCTAAACAAGGGGCCCACATTAACGTTTTGAGCTTGAAGATCTTTACCGCGATCTATTTTTTCTCCTAGATTACCGTACACACCGACACTTCCATC", debug=True)

