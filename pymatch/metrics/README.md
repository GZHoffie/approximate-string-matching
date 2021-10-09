# Metrics of Approximate String Matching

Prviously, edit distance or gap/mismatch penalty have been used to check whether two DNAs match well. However, it might not be the most suitable metric when it comes to large error rate. When error rate is large, mismatch/insertion/deletion happen very frequently, and

- The optimal match with the smallest edit distance is likely to be not unique.
- The match with the smallest edit distance might not reflect where the real indel happen.

With this, we introduce another metric called **coverage**, where we

- ignore the errornous areas, and focus on the areas where there are **long consecutive matches**, which are usually unique.
- Check whether our algorithm covers all those areas.

To do this we introduce two concepts:

1. **long consecutive matches (LCM)**: In the matching, we find the consecutive matches that are longer than a threshold. For example, we can set the threshold as 3, and look at the following matching:

   ```
   AGAGCTAAAC-ATGGCCGCACATAAATCGTTTTGAG-TTGAA-A-CTTTACCGCTGCATCTA-TTTTTCTCCTAGAATTATACCGTACACAGCCGAC-GTTCCACC
   AGAGCTAAACAAGGGGCCCACATTAA-CGTTTTGAGCTTGAAGATCTTTACCGC-G-ATCTATTTTTTCTCCTAG-A-T-TACCGTACACA-CCGACACTTCCATC
   ```

   Here, `AGAGCTAAAC` is a long consecutive match, while single matches like `A`, `GG` are ignored. We can link all long consecutive matches together, and form a **LCM substring**, like the following,

   ```
   AGAGCTAAACCACATCGTTTTGAGTTGAACTTTACCGCATCTATTTTTCTCCTAGTACCGTACACACCGACTTCCAC
   ```

   This indicates the parts of the reference and read that we are sure that are matching.

2. **String coverage**: We say that string `A` **covers** string `B` if we can construct `A` by simply inserting alphabets to `B`. For example,

   ```
   AGAGCTAAACAGGCACATAACGTTTTGAGTTGAACTTTACCGCATCTATTTTTCTCCTAGATTACCGTACACACCGACTTCCAC
   ```

   covers

   ```
   AGAGCTAAACGGCACATAACGTTTTGAGTTGAACTTTACCGCATCTATTTTTCTCCTAGTACCGTACACACCGACTTCCAC
   ```

## Ideally, ...

We want to show that the LCM substring of greedy algorithm covers the LCM substring of other traditional algorithms (Needleman-Wunsch/Landau-Vishkin).

However, I found that with higher threshold, the coverage rate is smaller. With the correct rate (rate that edit distance calculated by the greedy algorithm is optimal) being 0.81, the coverage rate is 0.82 if threshold is 0 (we are counting all the matches) and 0.72 if threshold is 1.

This is because in greedy algorithm, some matches are included in short match while the same matches are in long consecutive match when using traditional algorithms.

## To solve this, ...

I just find the LCM substring of the traditional algorithm, while keeping all the match in the result of greedy algorithm (LCM substring with threshold=0). The coverage rate would become 0.957, which is much higher than the correct rate.