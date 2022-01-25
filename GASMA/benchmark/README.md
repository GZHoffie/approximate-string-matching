# Benchmark for GASMA

Use the following lines to generate dataset

```c
./GASMA/bin/generate_dataset -n 1000000 -l 100 -e 0.3 -o ./pymatch/test/resource/sample.random.dataset.seq
```

/home/zhenhao/approximate-string-matching/GASMA/cmake-build-debug/hurdle-matrix-benchmark
Processed data file: simulated_5000000_100_0.050000_lt_eq.seq
...processed 100000 reads.
...processed 200000 reads.
...processed 300000 reads.
...processed 400000 reads.
...processed 500000 reads.
...processed 600000 reads.
...processed 700000 reads.
...processed 800000 reads.
...processed 900000 reads.
...complete.
===================== Benchmark Results =====================
Total number of alignments: 1000000
[Time]
=> Needleman-Wunsch | 30.560 s
=> LEAP             | 1.150 s
=> Greedy           | 1.680 s
[Accuracy] (percentage of alignments matching optimal penalty)
=> Needleman-Wunsch | 100.000 %
=> LEAP             | 99.461 %
=> Greedy           | 99.741 %
[Coverage] (percentage of alignments covering all long consecutive matches)
=> Greedy           | 99.913 %
Processed data file: simulated_5000000_100_0.100000_lt_eq.seq
...processed 100000 reads.
...processed 200000 reads.
...processed 300000 reads.
...processed 400000 reads.
...processed 500000 reads.
...processed 600000 reads.
...processed 700000 reads.
...processed 800000 reads.
...processed 900000 reads.
...complete.
===================== Benchmark Results =====================
Total number of alignments: 1000000
[Time]
=> Needleman-Wunsch | 27.340 s
=> LEAP             | 2.030 s
=> Greedy           | 2.970 s
[Accuracy] (percentage of alignments matching optimal penalty)
=> Needleman-Wunsch | 100.000 %
=> LEAP             | 97.642 %
=> Greedy           | 98.142 %
[Coverage] (percentage of alignments covering all long consecutive matches)
=> Greedy           | 99.322 %
Processed data file: simulated_5000000_100_0.150000_lt_eq.seq
...processed 100000 reads.
...processed 200000 reads.
...processed 300000 reads.
...processed 400000 reads.
...processed 500000 reads.
...processed 600000 reads.
...processed 700000 reads.
...processed 800000 reads.
...processed 900000 reads.
...complete.
===================== Benchmark Results =====================
Total number of alignments: 1000000
[Time]
=> Needleman-Wunsch | 27.270 s
=> LEAP             | 2.870 s
=> Greedy           | 3.610 s
[Accuracy] (percentage of alignments matching optimal penalty)
=> Needleman-Wunsch | 100.000 %
=> LEAP             | 94.712 %
=> Greedy           | 94.004 %
[Coverage] (percentage of alignments covering all long consecutive matches)
=> Greedy           | 97.667 %
Processed data file: simulated_5000000_100_0.200000_lt_eq.seq
...processed 100000 reads.
...processed 200000 reads.
...processed 300000 reads.
...processed 400000 reads.
...processed 500000 reads.
...processed 600000 reads.
...processed 700000 reads.
...processed 800000 reads.
...processed 900000 reads.
...complete.
===================== Benchmark Results =====================
Total number of alignments: 1000000
[Time]
=> Needleman-Wunsch | 26.220 s
=> LEAP             | 3.570 s
=> Greedy           | 4.430 s
[Accuracy] (percentage of alignments matching optimal penalty)
=> Needleman-Wunsch | 100.000 %
=> LEAP             | 92.481 %
=> Greedy           | 90.190 %
[Coverage] (percentage of alignments covering all long consecutive matches)
=> Greedy           | 96.044 %

Process finished with exit code 0





/home/zhenhao/approximate-string-matching/GASMA/cmake-build-debug/hurdle-matrix-benchmark
Processed data file: simulated_5000000_100_0.050000_lt_eq.seq
...processed 100000 reads.
...processed 200000 reads.
...processed 300000 reads.
...processed 400000 reads.
...processed 500000 reads.
...processed 600000 reads.
...processed 700000 reads.
...processed 800000 reads.
...processed 900000 reads.
...complete.
===================== Benchmark Results =====================
Total number of alignments: 1000000
[Time]
=> Needleman-Wunsch | 29.910 s
=> LEAP             | 1.490 s
=> Greedy           | 1.380 s
[Accuracy] (percentage of alignments matching optimal penalty)
=> Needleman-Wunsch | 100.000 %
=> LEAP             | 99.461 %
=> Greedy           | 99.741 %
[Coverage] (percentage of alignments covering all long consecutive matches)
=> Greedy           | 0.000 %
Processed data file: simulated_5000000_100_0.100000_lt_eq.seq
...processed 100000 reads.
...processed 200000 reads.
...processed 300000 reads.
...processed 400000 reads.
...processed 500000 reads.
...processed 600000 reads.
...processed 700000 reads.
...processed 800000 reads.
...processed 900000 reads.
...complete.
===================== Benchmark Results =====================
Total number of alignments: 1000000
[Time]
=> Needleman-Wunsch | 28.790 s
=> LEAP             | 1.940 s
=> Greedy           | 1.890 s
[Accuracy] (percentage of alignments matching optimal penalty)
=> Needleman-Wunsch | 100.000 %
=> LEAP             | 97.642 %
=> Greedy           | 98.142 %
[Coverage] (percentage of alignments covering all long consecutive matches)
=> Greedy           | 0.000 %
Processed data file: simulated_5000000_100_0.150000_lt_eq.seq
...processed 100000 reads.
...processed 200000 reads.
...processed 300000 reads.
...processed 400000 reads.
...processed 500000 reads.
...processed 600000 reads.
...processed 700000 reads.
...processed 800000 reads.
...processed 900000 reads.
...complete.
===================== Benchmark Results =====================
Total number of alignments: 1000000
[Time]
=> Needleman-Wunsch | 26.570 s
=> LEAP             | 2.820 s
=> Greedy           | 2.580 s
[Accuracy] (percentage of alignments matching optimal penalty)
=> Needleman-Wunsch | 100.000 %
=> LEAP             | 94.712 %
=> Greedy           | 94.004 %
[Coverage] (percentage of alignments covering all long consecutive matches)
=> Greedy           | 0.000 %
Processed data file: simulated_5000000_100_0.200000_lt_eq.seq
...processed 100000 reads.
...processed 200000 reads.
...processed 300000 reads.
...processed 400000 reads.
...processed 500000 reads.
...processed 600000 reads.
...processed 700000 reads.
...processed 800000 reads.
...processed 900000 reads.
...complete.
===================== Benchmark Results =====================
Total number of alignments: 1000000
[Time]
=> Needleman-Wunsch | 26.080 s
=> LEAP             | 3.570 s
=> Greedy           | 3.060 s
[Accuracy] (percentage of alignments matching optimal penalty)
=> Needleman-Wunsch | 100.000 %
=> LEAP             | 92.481 %
=> Greedy           | 90.190 %
[Coverage] (percentage of alignments covering all long consecutive matches)
=> Greedy           | 0.000 %

Process finished with exit code 0
