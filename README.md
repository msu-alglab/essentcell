# conflict-free-ILP

This ILP Implementation is inspired by the paper "Studying the History of Tumor Evolution from Single-Cell Sequencing Data by Exploring the Space of Binary Matrices" by Malikic et. al. 

Given a binary m x n matrix, D, and a positive integer, k, this ILP finds the minimum number of entries that must be flipped from 0 to 1 to produce a conflict free matrix, if k entries are allowed to be flipped from 1 to 0. 

A matrix contains a conflict if it contains the submatrix [0 1 / 1 0 / 1 1]. 

The input matrix, D, corresponds to a sample of cancerous cells, where each of the m rows corresponds to a distinct cell, and each of the n columns corresponds to each of the n mutations. An entry Dij = 1 if and only if cell i has mutation j, otherwise it is 0. 
