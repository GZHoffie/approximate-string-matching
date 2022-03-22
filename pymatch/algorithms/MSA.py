import re
import numpy as np

def _to_PSSM(char):
    assert char in {'A', 'C', 'G', 'T', '-'}, f"Invalid character {char}."
    pssm = {
        'A': np.array([1, 0, 0, 0, 0]),
        'C': np.array([0, 1, 0, 0, 0]),
        'G': np.array([0, 0, 1, 0, 0]),
        'T': np.array([0, 0, 0, 1, 0]),
        '-': np.array([0, 0, 0, 0, 1])
    }
    return pssm.get(char)


"""
The class for profile-profile alignment
"""
class ProfileProfileAlignment:
    """
    The initializaer of the class.
    """
    def __init__(self, a1, a2) -> None:
        self.a1 = a1
        self.a2 = a2
        self.pssm1 = self._create_PSSM(self.a1)
        self.pssm2 = self._create_PSSM(self.a2)

        self.match_score = 1
        self.mismatch_score = -2

        self.score_matrix = np.zeros((5, 5))
        self.score_matrix.fill(self.mismatch_score)
        for i in range(self.score_matrix.shape[0]):
            if i != self.score_matrix.shape[0] - 1:
                self.score_matrix[i, i] = self.match_score
            else:
                self.score_matrix[i, i] = 0

        print(self.score_matrix)

        self.D = np.zeros((len(self.a1[0]) + 1, len(self.a2[0]) + 1))
        print(self._psp(1, None))
        self._dp()

    def _create_PSSM(self, alignment):
        pssm = np.zeros((len(alignment[0]), 5))

        for i in range(len(alignment[0])):
            for j in range(alignment.shape[0]):
                pssm[i, :] += _to_PSSM(alignment[j][i])
            
            pssm[i, :] /= np.sum(pssm[i, :])
        
        print(pssm)
        
        return pssm
    
    def _psp(self, i, j):
        a1 = self.pssm1[i, :] if i is not None else np.array([0, 0, 0, 0, 1])
        a2 = self.pssm2[j, :] if j is not None else np.array([0, 0, 0, 0, 1])

        # Reshape
        a1 = np.reshape(a1, (1, 5))
        a2 = np.reshape(a2, (5, 1))

        # Calculate PSP
        psp = np.dot(np.dot(a1, self.score_matrix), a2)
        return psp[0][0]
    
    def _dp(self):
        # Init
        for i in range(self.D.shape[0] - 1):
            self.D[i + 1, 0] = self.D[i, 0] + self._psp(i, None)
        for j in range(self.D.shape[0] - 1):
            self.D[0, j + 1] = self.D[0, j] + self._psp(j, None)
        
        # Dynamic programming step
        for i in range(self.D.shape[0] - 1):
            
            for j in range(self.D.shape[0] - 1):
                match = self.D[i, j] + self._psp(i, j)
                insert = self.D[i, j+1] + self._psp(i, None)
                delete = self.D[i+1, j] + self._psp(None, j)
                self.D[i+1, j+1] = max(match, insert, delete)
        
        print(self.D)




  




if __name__ == "__main__":
    ppa = ProfileProfileAlignment(np.array(["ACGT-CA", "AGGTCCA"]), np.array(["-A-CTCC", "TAGCTCC"]))

