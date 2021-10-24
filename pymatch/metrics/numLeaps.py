class numLeaps:
    def __init__(self, dict1) -> None:
        self.match = dict1
    
    def compute(self):
        numLeaps = 0
        for i in range(len(self.match["dna1"])):
            if self.match["dna1"][i] == "-" and (i == 0 or self.match["dna1"][i-1] != "-"):
                numLeaps += 1
            if self.match["dna2"][i] == "-" and (i == 0 or self.match["dna2"][i-1] != "-"):
                numLeaps += 1
        
        return numLeaps