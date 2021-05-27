# What is a sequence alignment?
One of the most problems that a DNA sequence changes, most frequently due to random errors in replication (or the copying of a DNA sequence) which causes mutations(substitutions,insertions,and deletions) .
by using Pairwise sequence alignment methods to find the best-matching piecewise (local or global) alignments of two query sequences. Pairwise alignments can only be used between two sequences with high similarity  by aligning the sequences to one another inserting gaps.

---
<br />

## the Smith-Waterman algorithm performs local sequence alignment which is to find a partially matching substring of two sequences.

---
<br />

## Step 1: Scoring matrix
 To find the local alignment of two sequence the Smith-Waterman calculates a scoring matrix first. The following code calculates this matrix for two strings with linear gap costs. For performance reasons I went for an implementation with NumPy arrays. Values for match scores and gap costs can be changed. The dimensions of the scoring matrix are 1+length of each sequence respectively. All the elements of the first row and the first column are set to 0. The extra first row and first column make it possible to align one sequence to another at any position, and setting them to 0 makes the terminal gap free from penalty.Score each element from left to right, top to bottom in the matrix, considering the outcomes of substitutions (diagonal scores) or adding gaps (horizontal and vertical scores). If none of the scores are positive, this element gets a 0. Otherwise the highest score is used and the source of that score is recorded.One of the most important distinctions is that no negative score is assigned in the scoring system of the Smith–Waterman algorithm, which enables local alignment. When any element has a score lower than zero, it means that the sequences up to this position have no similarities; this element will then be set to zero to eliminate influence from previous alignment. In this way, calculation can continue to find alignment in any position afterwards.

<br />

![smith-waterman matrix creation](https://upload.wikimedia.org/wikipedia/commons/2/28/Smith-Waterman-Algorithm-Example-Step2.png)

---
<br />

## Step 2: Backtracing

The second step is backtracing from the calculated scoring matrix to calculate the optimal alignment. we had to append with the NumPy arrays a bit because numpy.argmax() will only find the first indexes of a maximum value inside an array and backtracing involves starting with the last occurrence of the maximum .The initial scoring matrix of Smith–Waterman algorithm enables the alignment of any segment of one sequence to an arbitrary position in the other sequence. 

<br/>

![Smith-Waterman-Algorithm-Example-Step3](https://upload.wikimedia.org/wikipedia/commons/e/e6/Smith-Waterman-Algorithm-Example-Step3.png)
