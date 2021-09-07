from pymatch.algorithms import NeedlemanWunsch

test_file = "/home/zhenhao/approximate-string-matching/pymatch/test/resource/sample.random.dataset.seq"
test_items = 100000
i = 0

with open(test_file, "r") as f:
    while True:
        line = f.readline()
        if not line:
            break
        str1 = line[1:][:-1]
        str2 = f.readline()[1:][:-1]
        nw = NeedlemanWunsch(str1, str2)
        print(nw.editDistance())
        i += 1
        if i > test_items:
            break