from pymatch.algorithms import GASMA

test_file = "./resource/sample.dataset.seq"

with open(test_file, "r") as f:
    while True:
        line = f.readline()
        if not line:
            break
        str1 = line[1:][:-1]
        str2 = f.readline()[1:][:-1]
        GASMA(str1, str2, 5, threshold=3)
