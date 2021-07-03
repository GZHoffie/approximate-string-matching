# Bioinformatics Notes

## Basic Concepts

### Global/Local/Semi-global alignment

>A global alignment is defined as the *end-to-end* alignment of two strings *s* and *t*.
>
>A local alignment of string s and t is an alignment of *substrings* of *s* with substrings of *t*.
>
>A semi-global alignment of string s and t is an alignment of a substring of s with a substring of t.

### Smith-Waterman algorithm

> The **Smith–Waterman algorithm** performs local [sequence alignment](https://en.wikipedia.org/wiki/Sequence_alignment); that is, for determining similar regions between two strings of [nucleic acid sequences](https://en.wikipedia.org/wiki/Nucleic_acid_sequence) or [protein sequences](https://en.wikipedia.org/wiki/Protein_sequence). Instead of looking at the [entire](https://en.wikipedia.org/wiki/Needleman–Wunsch_algorithm) sequence, the Smith–Waterman algorithm compares segments of all possible lengths and [optimizes](https://en.wikipedia.org/wiki/Mathematical_optimization) the [similarity measure](https://en.wikipedia.org/wiki/Similarity_measure).

**Input**: Two sequences $A=a_1a_2\dots a_n$ and $B=b_1b_2\dots b_m$, e.g. `TACGGGCCCGCTAC` and `TAGCCCTATCGGTCA`,

**Output**: their similarity according to how many characters align, for example,

```
TACGGGCCCGCTA-C
||   | || ||| |
TA---G-CC-CTATC
```

Let $A=a_1a_2\dots a_n$ and $B=b_1b_2\dots b_m$ be sequences to be aligned, we use the following steps

1. Determine the substitution matrix and the gap penalty scheme.

	- ==Substitution matrix== is a matrix, which gives reward to matches and penalties to mismatriches, for example we may have a substitution matrix that looks like

		|       | A    | G    | C    | T    |
		| ----- | ---- | ---- | ---- | ---- |
		| **A** | 3    | -3   | -3   | -3   |
		| **G** | -3   | 3    | -3   | -3   |
		| **C** | -3   | -3   | 3    | -3   |
		| **T** | -3   | -3   | -3   | 3    |

		Or in compact form, 
		$$
		s(a_i,b_j)=\left\{\begin{array}{ll}+3, &a_i=b_j\\-3,&a_i\ne b_j\end{array}\right.
		$$
		different base substitutions can have different scores.

	- ==Gap penalty== $W_k$ is the penalty of a gap that has length $k$. It can be **linear** ($W_k=kW_1$), or **affine** ($W_k=uk+v\, (u>0,v>0)$). e.g. we can choose $W_k=2k$, a linear penalty.

1. Construct a scoring matrix $H$ of size $(n+1)\times (m+1)$, with the first row and first column filled with 0.
	![img](https://upload.wikimedia.org/wikipedia/commons/2/2c/Smith-Waterman-Algorithm-Example-Step1.png)

1. Fill the scoring matrix with
	$$
	H_{ij}=\max\left\{\begin{array}{l}H_{i-1,j-1}+s(a_i,b_j),\\\max_{k\geq 1}\{H_{i-k,j}-W_k\}\\ \max_{l\geq 1}\{H_{i,j-l}-W_l\}\\0\end{array}\right.
	$$
	The four options are:

	- Aligning $a_i$ and $b_j$
	- Let $a_i$ be the end of a gap of length $k$
	- Let $b_j$ be the end of a gap of length $l$
	- no alignment up to $a_i$ and $b_j$ (fresh start)

	Choosing a linear penalty reduce the calculation to
	$$
	H_{ij}=\max\left\{\begin{array}{l}H_{i-1,j-1}+s(a_i,b_j),\\H_{i-1,j}-W_1\\H_{i,j-1}-W_1\\0\end{array}\right.
	$$

1. Traceback. Find the highest score in the matrix and step backward by finding how the score comes from. Starts at the **highest score** and ends at **first zero encountered**.

- Enables local alignment: having 0 as an option gives an opportunity for the alignment to restart: previous strings has no similarity.

### Needleman-Wunsch algorithm

Also using dynamic programming to find the best alignments. We follow the following steps:

1. Choosing a ==scoring system==. Given an alignment, we can score it by counting

	- **Match**: the two letters at current index are the same
	- **Mismatch**: the two letters at current index are different
	- **Indel**: One letter aligning to a gap.

	And we can score Match by +1 and Mismatch/Indel by -1. We may also define a substitution matrix and gap penaly $d<0$.

1. We also create a matrix $F$ of size $(n+1)\times(m+1)$, and initialize it such that
	$$
	F_{i0}=di,\qquad F_{0j}=dj
	$$

1. We use a similar calculation as before to fill the whole matrix
	$$
	F_{ij}=\max(F_{i-1,j-1}+s(a_i,b_j),F_{i,j-1}+d,F_{i-1,j}+d)
	$$

1. Backtracking is also used to find the optimal alignment, starting at the **lower right of matrix** and ends at **top left cell**. 

- Considers full sequence alignment.

### Hirschberg's algorithm

Modification of Needleman-Wunsch that needs only $\mathcal{O}(\min\{m,n\})$ space. We save the space by only looking at the last two lines in the Needleman-Wunsch matrix.

**Basic Idea**: DP + Divide & Conquer

**Steps**: choose a column, find the node in the column such that the longest path passes throught this node. Then find the lognest path from top left to the node and from the node to the bottom right.

**Nice Reference**: [Youtube Video: Space-Efficient Sequence Alignment](https://www.youtube.com/watch?v=3TfDm8GpWRU)

### Method of Four Russians

Speeding up algorithms involving Boolean matrices, or matrices in which each cell may take on only a bounded number of possible values.

## LEAP

### Leaping Toad Problem (LTP)

A traversal problem in a special DAG, and find the shortest path. This graph is

* A convex 2-d vertex grid. Each row is a ==lane==, and within each row we have ==forward edges== pointing to the right.
* There are ==leap edges== connecting vertices in different lanes. There cannot be multiple edge pointing to the same lane from the same vertex. All vertices in the same lane share the same types of leap edges.
* Leap edges connecting two lanes share same, positive weight. Forward edges have zero/positive weight. If it has positive weight then it's a ==hurdle==.

**Goal**: start from the first vertex of an origin lane, either travels to the last vertes of a destination lane, or travels out of the swimming pool while exiting onto a destination lane.

### Approximate String Matching

1. Assign each element $D_{i,j}$ in the distance matrix a unique vertex.

1. Draw a directional edge $(v_{i,j},v_{i',j'})$ iff $v_{i,j}\ne v_{i',j'}$, $i'-i\leq 1$, $j'-j\leq 1$.

1. On each edge $(v_{i,j},v_{i',j'})$ we set
	$$
	w=\left\{\begin{array}{ll}0,&\text{if }i'=i+1,j'=j+1,s_i=r_j\text{ (a match)}\\1,&\text{otherwise}\end{array}\right.
	$$

1. Find minimum total edge weight from $v_{0,0}$ to $v_{L,L}$.

We can convert that to a swimming pool, and the goal becomes finding the minimum path from first vertex of center lane to the last of the center lane.
