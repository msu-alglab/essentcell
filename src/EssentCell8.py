"""
Author: Robyn Burger and Allison Shi and Adiesha Liyanage and Brendan Mumey
Date: 2024-08-07
Description: Given a positive integer k (pre-determined # of false positive(s)),
and n x m binary mutation matrix (sorted), EssentCell4 produces the 
corresponding essential order diagram. This information will be written to
results/inputFile folder.

****WITH MULTIPLICITY without y variables and with improved group testing (single u and v testing with transitive
edges being tested) ***

Optionally produced verbose detail file, containing information related to ILP calls

Example Usage:

$ python EssentCell5.py -h

$ python EssentCell5.py data.sorted.csv 3 --verbose -result_folder result4

$ python EssentCell5.py data.sorted.csv 3

$ python EssentCell5.py smallest.sorted.csv 2 -result_folder ../Results5 --verbose -random_u_sample_size_for_group_testing 2
 """
import argparse
import os
import random
import time

import networkx as nx
import numpy as np
import graphviz
import networkx
import pandas as pd
from gurobipy import *
from networkx import *
from itertools import product

parser = argparse.ArgumentParser(description="Arguments for the EssentCell program")
parser.add_argument('inputFile', type=str, help="Sorted Input file to the program")
parser.add_argument('k', type=int, default=1, help="k value for the input program")
parser.add_argument('-result_folder', type=str, default="results4", help="The result folder name")
parser.add_argument('--verbose', action='store_true', help='Increase output verbosity')
parser.add_argument('-print_trace_of_constraint', action='store_true', help='Increase the output verbosity of '
                                                                            'constraints')
parser.add_argument('-random_u_sample_size_for_group_testing', type=int, default=8, help='This is the size of the '
                                                                                         'random sample size for '
                                                                                         'initial group testing ')
args = parser.parse_args()

print(f"Input file: {args.inputFile}")
print(f"Value of k: {args.k}")
print(f"Result folder name {args.result_folder}")
print(f"Verbosity: {args.verbose}")
print(f"Print Trace Verbosity: {args.print_trace_of_constraint}")
print(f"Initial random sample size for group testing: {args.random_u_sample_size_for_group_testing}")

fileName = args.inputFile
outputName = fileName[:-4]
k = args.k
verbose = args.verbose
print_trace = args.print_trace_of_constraint
count = 0

df = pd.read_csv(fileName)
E = df.to_numpy()
df.drop_duplicates(inplace=True)
D = df.to_numpy()

n = D.shape[0]
m = D.shape[1]

# # Calculate the multiplicities of each row i in D
M = {}
for i in range(D.shape[0]):
    dups = 0
    for x in range(E.shape[0]):
        if np.all(D[i] == E[x]):
            dups += 1
    M[i] = dups

start_time = time.time()

"""
Returns sigma, the minimum number of bit flips required to make D conflict-free. 
"""


def FindOpt():
    try:
        # Silence console output
        env = Env(empty=True)
        env.setParam("OutputFlag", 0)
        env.start()
        model = Model("min_flip_model", env=env)
        model.Params.LogToConsole = 0

        # X will be constrained to be a conflict-free matrix
        X = model.addMVar((n, m), vtype=GRB.BINARY, name="X")

        # total = sum(sum(M[i]*(1 - D[i, j])*(X[i, j]) + k * M[i] * (D[i, j])*(1 - X[i, j]) for j in range(m)) for i in range(n))
        total = sum(sum(M[i] * (1 - D[i, j]) * X[i, j] for j in range(m)) for i in range(n))  # new obj function 8/5
        model.setObjective(total, GRB.MINIMIZE)

        B01 = model.addMVar((m, m), vtype=GRB.BINARY, name="B01")
        B10 = model.addMVar((m, m), vtype=GRB.BINARY, name="B10")
        B11 = model.addMVar((m, m), vtype=GRB.BINARY, name="B11")

        model.addConstrs(
            -1 * X[i, p] + X[i, q] <= B01[p, q] for p in range(m) for q in range(p + 1, m) for i in range(n))  # (1)
        model.addConstrs(
            X[i, p] - X[i, q] <= B10[p, q] for p in range(m) for q in range(p + 1, m) for i in range(n))  # (2)
        model.addConstrs(
            X[i, p] + X[i, q] - 1 <= B11[p, q] for p in range(m) for q in range(p + 1, m) for i in range(n))  # (3)
        model.addConstrs(B01[p, q] + B10[p, q] + B11[p, q] <= 2 for p in range(m) for q in range(p + 1, m))  # (4)
        # model.addConstrs(sum(sum((1 - X[i, j]) * D[i,j] == k for j in range(m)) for i in range(n))) # (5) updated on 8/5
        glb_cons = model.addConstr(sum(M[i] * (1 - X[i, j]) * D[i, j] for i in range(n) for j in range(m)) == k,
                                   name="global_constraint")  # 8/5 test
        model.update()
        if print_trace:
            print("printing global_constraint")
            print(f"{model.getRow(glb_cons)} {glb_cons.Sense} {glb_cons.RHS}")

        model.optimize()
        model.write("model.initial.lp")
        sig = model.ObjVal
        print(f"sig: {sig}")
        return sig

    except GurobiError as ex:
        print(f"*********ERROR*********\n{ex}")
        return -1


"""
Given a set of integers, V, Split(V) returns V partitioned into 2 subsets. 
"""


def Split(V):
    # V = list(V)
    i = int(len(V) / 2)
    return V[:i], V[i:]


"""
Returns the essential order relation of an input D matrix.
"""


def test_ESS(k, U, V, sig):
    # if len(U) == 0 or len(V) == 0:
    #     return True
    global D
    global count
    count += 1
    if count % 50 == 0:
        print(f"Essential relation ILP call {count}")
    try:
        # Silence console output
        env = Env(empty=True)
        env.setParam("OutputFlag", 0)
        env.start()
        model = Model("min_flip_model", env=env)
        model.Params.LogToConsole = 0

        X = model.addMVar((n, m), vtype=GRB.BINARY, name="X")
        B01 = model.addMVar((m, m), vtype=GRB.BINARY, name="B01")
        B10 = model.addMVar((m, m), vtype=GRB.BINARY, name="B10")
        B11 = model.addMVar((m, m), vtype=GRB.BINARY, name="B11")

        # total = sum(sum(M[i]*(1 - D[i, j])*(X[i, j]) + k*M[i]*(D[i, j])*(1 - X[i, j]) for j in range(m)) for i in range(n))
        # total = sum(sum(M[i] * (1 - D[i, j]) * X[i, j] for j in range(m)) for i in range(n))  # new obj function 8/5
        model.setObjective(0, GRB.MINIMIZE)

        if print_trace:
            print("printing minimizing objective")
            print(model.getObjective())

        # Numbers to the right of each constraint correspond to those in the paper
        model.addConstrs(
            X[i, q] - X[i, p] <= B01[p, q] for p in range(m) for q in range(p + 1, m) for i in range(n))  # (1)
        model.addConstrs(
            X[i, p] - X[i, q] <= B10[p, q] for p in range(m) for q in range(p + 1, m) for i in range(n))  # (2)
        model.addConstrs(
            X[i, p] + X[i, q] - 1 <= B11[p, q] for p in range(m) for q in range(p + 1, m) for i in range(n))  # (3)
        model.addConstrs(B01[p, q] + B10[p, q] + B11[p, q] <= 2 for p in range(m) for q in range(p + 1, m))  # (4)
        glb_cons = model.addConstr(sum(M[i] * (1 - X[i, j]) * D[i, j] for i in range(n) for j in range(m)) == k,
                                   name="global_constraint")  # 8/5 test

        model.update()
        if print_trace:
            print("printing global_constraint in an essential ILP call")
            print(f"{model.getRow(glb_cons)} {glb_cons.Sense} {glb_cons.RHS}")

        for u in U:
            nz = len(V)
            z = model.addMVar((m, nz), vtype=GRB.BINARY, name=f"z_{u}")
            model.update()

            # Test for EO
            for i in range(m):
                for v_index in range(nz):
                    v = V[v_index]
                    model.addConstr(X[u, i] - X[v, i] <= z[i, v_index])  # (10)
                    model.addConstr(z[i, v_index] <= (X[u, i] - X[v, i] + 1) / 2)  # (10)

            for v in range(nz):
                model.addConstr(sum(z[i, v] for i in range(m)) >= 1)  # (11)

        model.addConstr(sum(sum(M[i] * (1 - D[i, j]) * X[i, j] for i in range(n)) for j in range(m)) == sig)  # (12) 8/5
        # model.addConstr(sum(sum(M[i]*(1 - D[i, j])*(X[i, j]) + k * M[i] * (D[i, j])*(1 - X[i, j]) for j in range(m)) for i in range(n)) == sig) # (8)

        model.update()
        model.optimize()
        if print_trace:
            model.write("model.ess.lp")
        if model.Status == 3:
            print(U, V, " infeasible")
            return False
        else:
            print(U, V, " feasible *****")
            return True

    except GurobiError as ex:
        print(f"*********ERROR*********\n{ex}")


def k_ess(k, U, V, sig):
    if (set(U) & set(V)) or (not test_ESS(k, U, V, sig)):
        if len(U) * len(V) <= 16:  # this could be an optional input parameter
            result = set()
            for u, v in product(U, V):
                if u == v or (not test_ESS(k, [u], [v], sig)):
                    result.add((u, v))
            return result
        if len(U) > len(V):
            U_L, U_R = Split(U)
            return k_ess(k, U_L, V, sig).union(k_ess(k, U_R, V, sig))
        else:
            V_L, V_R = Split(V)
            return k_ess(k, U, V_L, sig).union(k_ess(k, U, V_R, sig))
    return set()


def calculate_essential_for_given_k(k):
    sig = FindOpt()

    S = list(range(D.shape[0]))

    G = networkx.DiGraph()
    G.add_nodes_from(S)
    for u in S:
        for v in S:
            if u != v:
                if nx.has_path(G, u, v) or not test_ESS(k, [u], [v], sig):
                    G.add_edge(u, v)
    return G


fileResults = outputName + "." + str(k) + ".esspairs" + ".txt"
print(f"The essential pairs will be written to {fileResults}")
pathToResult = f"{args.result_folder}/{outputName}/"
os.makedirs(pathToResult, exist_ok=True)

Graph = calculate_essential_for_given_k(k)
end_time = time.time()
networkx.write_edgelist(Graph, f"{pathToResult}{fileResults}")
if verbose:
    verboseResultFile = f"{pathToResult}{outputName}.{k}.esspairs.verbose.txt"
    print(f"Verbose result is written into {verboseResultFile}")
    f = open(verboseResultFile, "a")
    f.write(f"k value: {k}\n")
    f.write(f"n (number of samples): {n}\n")
    f.write(f"m (number of mutations): {m}\n")
    f.write(f"EssILP calls: {count}\n")
    f.write(f"Runtime: {end_time - start_time} seconds\n")
    f.write(f"Number of Nodes: {Graph.number_of_nodes()}\n")
    f.write(f"Initial sample size for u: {args.random_u_sample_size_for_group_testing}\n")
    # f.write(f"Poset Width: {Width(Graph)}\n")
    f.write(f"Essential Relation: {[edge for edge in Graph.edges]}\n\n")
    print(f"\nResults written was to file: {verboseResultFile}")
    f.close()
