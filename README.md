# EssentCell 

Cancer is an evolutionary disease. As cells multiply, they can take on new mutations, and through the accumulation of those mutations, they can become cancerous. Tracking this evolution is helpful for understanding how the disease progresses and the most effective ways to treat it. We model this evolution in diagrams called phylogenetic trees. Much like a family tree, a phylogenetic tree depicts the relationships between generations of cells. These models depend on taking physical samples from cancer patients. 

But, these samples are imperfect, often leading to multiple phylogenetic trees. EssentCell finds the lineages that stay consistent among all optimal phylogenetic trees. These essential relations provide insight into the most probable ancestral relationships, painting a more precise picture of a tumor's evolutionary history. 


# Dependencies 

EssentCell is coded in Python 3.11.7
 * Gurobi 11.0.2
 * NumPy 1.26.3
 * Pandas 2.2.0
 * NetworkX 3.2.1


# Data Sources 

This ILP Implementation is inspired by \[2\]. 

Simulated data was obtained from \[\4\], which was produced using sequencing data of 77 patient with Acute Myeloid Leukemia (AML) from \[\3\]. Real data is from study of six acute lymphoblastic leukemia (ALL) patients in \[1\]. 


# References 
\[1\] GAWAD, C., KOH, W., AND QUAKE, S. R. Dissecting the clonal origins
of childhood acute lymphoblastic leukemia by single-cell genomics.
Proceedings of the National Academy of Science 111, 50 (Dec. 2014),
17947–17952.

\[2\] MALIKIC, S., MEHRABADI, F. R., AZER, E. S., EBRAHIM-ABADI,
M. H., AND SAHINALP, S. C. Studying the history of tumor evolution
from single-cell sequencing data by exploring the space of binary
matrices. J. Comput. Biol. 28, 9 (2021), 857–879

\[3\] MORITA, K., WANG, F., JAHN, K., KUIPERS, J., YAN, Y., MATTHEWS,
J., LITTLE, L., GUMBS, C., CHEN, S., ZHANG, J., SONG, X., THOMP-
SON, E., PATEL, K., BUESO-RAMOS, C., DINARDO, C. D., RAVANDI,
F., JABBOUR, E., ANDREEFF, M., CORTES, J., KONOPLEVA, M.,
BHALLA, K., GARCIA-MANERO, G., KANTARJIAN, H., BEEREN-
WINKEL, N., NAVIN, N., FUTREAL, P. A., AND TAKAHASHI, K. Clonal
evolution of acute myeloid leukemia revealed by high-throughput single-cell genomics. bioRxiv (2020).

\[4\] WEBER, L. L., AND EL-KEBIR, M. Phyolin: Identifying a linear
perfect phylogeny in single-cell DNA sequencing data of tumors. In
20th International Workshop on Algorithms in Bioinformatics, WABI
2020, September 7-9, 2020, Pisa, Italy (Virtual Conference) (2020),
C. Kingsford and N. Pisanti, Eds., vol. 172 of LIPIcs, Schloss Dagstuhl-Leibniz-Zentrum fur Informatik, pp. 5:1–5:14.

<!-- By obtaining SCS data from tumor tissue and observing the presence or absence of certain mutations in each sample, we can reconstruct the tumor's probable evolutionary history. This evolution can be represented in a phylogenetic tree, where the root is a non-cancerous cell and successive generations of nodes are cells with malignant mutations. Since the infinite sites assumption (ISA) states that any particular mutation is gained exactly once and cannot be lost after it has been gained, all child nodes retain all the mutations of their parents \textbf{cite ISA}. 

Equivalently, we may represent the sample of cells as an $m$ by $n$ binary matrix, $B$, where there are $m$ single cell samples and $n$ possible mutations. An entry $b_{ij} = 1$ if and only if mutation $j$ is present in cell $i$, and $b_{ij} = 0$ otherwise. Furthermore, we introduce the notion of a ``conflict'' to indicate whether a binary mutation matrix can represent an evolutionary tree with a perfect phylogeny. A matrix contains a conflict if it contains the submatrix or a permutation of the rows or columns of the submatrix: [0 1 / 1 0 / 1 1]. 

Given a binary m x n matrix, D, and a positive integer, k, this ILP finds the minimum number of entries that must be flipped from 0 to 1 to produce a conflict free matrix, if k entries are allowed to be flipped from 1 to 0. 
 -->
