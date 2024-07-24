# conflict-free-ILP

This ILP Implementation is inspired by the paper "Studying the History of Tumor Evolution from Single-Cell Sequencing Data by Exploring the Space of Binary Matrices" by Malikic et. al. 

By obtaining SCS data from tumor tissue and observing the presence or absence of certain mutations in each sample, we can reconstruct the tumor's probable evolutionary history. This evolution can be represented in a phylogenetic tree, where the root is a non-cancerous cell and successive generations of nodes are cells with malignant mutations. Since the infinite sites assumption (ISA) states that any particular mutation is gained exactly once and cannot be lost after it has been gained, all child nodes retain all the mutations of their parents \textbf{cite ISA}. 

Equivalently, we may represent the sample of cells as an $m$ by $n$ binary matrix, $B$, where there are $m$ single cell samples and $n$ possible mutations. An entry $b_{ij} = 1$ if and only if mutation $j$ is present in cell $i$, and $b_{ij} = 0$ otherwise. Furthermore, we introduce the notion of a ``conflict'' to indicate whether a binary mutation matrix can represent an evolutionary tree with a perfect phylogeny. A matrix contains a conflict if it contains the submatrix or a permutation of the rows or columns of the submatrix: [0 1 / 1 0 / 1 1]. 

Given a binary m x n matrix, D, and a positive integer, k, this ILP finds the minimum number of entries that must be flipped from 0 to 1 to produce a conflict free matrix, if k entries are allowed to be flipped from 1 to 0. 

