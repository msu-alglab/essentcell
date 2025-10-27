"""
@author: Adiesha Liyanage and Brendan Mumey and Braeden Sopp
Date: 2025-01-07
Description: Given a positive integer k (pre-determined # of false positive(s)),
and n x m binary mutation matrix-possibly with missing values, EssentCell_With_GT.py produces the
corresponding essential order diagram. This information will be written to
results/inputFile folder.
"""
import argparse
import os
import time

import networkx
import pandas as pd
from gurobipy import *

parser = argparse.ArgumentParser(description="Arguments for the EssentCell program")
parser.add_argument('inputFile', type=str, help="Sorted Input file to the program")
parser.add_argument('k', type=int, default=1, help="k value for the input program")
parser.add_argument('-result_folder', type=str, default="results_default", help="The result folder name")
parser.add_argument('--verbose', action='store_true', help='Increase output verbosity')
parser.add_argument('-print_trace_of_constraint', action='store_true', help='Increase the output verbosity of '
                                                                            'constraints')
parser.add_argument('-timeout', type=float, default=float('inf'), help="Timeout value for ILP calls")
args = parser.parse_args()
print(f"Input file: {args.inputFile}")
print(f"Value of k: {args.k}")
print(f"Result folder name {args.result_folder}")
print(f"Verbosity: {args.verbose}")
print(f"Print Trace Verbosity: {args.print_trace_of_constraint}")
print(f"ILP timeout: {args.timeout}")

fileName = args.inputFile
outputName = fileName[:-4]  # remove the file extension
k = args.k
verbose = args.verbose
print_trace = args.print_trace_of_constraint
ilp_timeout = args.timeout
count = 0

df = pd.read_csv(fileName)
df.index = range(1, len(df) + 1)
E = df.to_numpy()
# df.drop_duplicates(inplace=True) //we do not delete duplicates anymore
D = df.to_numpy()

n = D.shape[0]
m = D.shape[1]

start_time = time.time()


def FindOpt():
    """
    Returns sigma, the minimum number of bit flips required to make D conflict-free.
    """
    try:
        # Silence console output
        # env = Env(empty=True)
        # env.setParam("OutputFlag", 0)
        # env.start()
        model = Model("min_flip_model")
        model.Params.LogToConsole = 0

        time_olp_creation_s = time.time()
        # X will be constrained to be a conflict-free matrix
        X = model.addMVar((n, m), vtype=GRB.BINARY, name="X")

        # total = sum(sum(M[i]*(1 - D[i, j])*(X[i, j]) + k * M[i] * (D[i, j])*(1 - X[i, j]) for j in range(m)) for i in range(n))
        total = sum(
            sum((1 - D[i, j]) * X[i, j] for j in range(m) if D[i, j] != -1) for i in range(n))  # new obj function 9/6
        # we check whether the original data point is not missing.
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
        glb_cons = model.addConstr(
            sum((1 - X[i, j]) * D[i, j] for i in range(n) for j in range(m) if D[i, j] == 1) == k,
            name="global_constraint")  # 9/5 test,
        # here we check whether D[i,j] == 1 to count the number changed false positives
        # model.update()
        if print_trace:
            print("printing global_constraint")
            print(f"{model.getRow(glb_cons)} {glb_cons.Sense} {glb_cons.RHS}")

        time_olp_creation_e = time.time()
        print(f"ILP creation time: {time_olp_creation_e - time_olp_creation_s}")
        time_ilp_start = time.time()
        model.optimize()
        time_ilp_end = time.time()
        print(f"ILP OPT time {time_ilp_end - time_ilp_start} seconds")
        model.write("model.initial.lp")
        sig = model.ObjVal
        print(f"sig: {sig}")
        return sig

    except GurobiError as ex:
        print(f"*********ERROR*********\n{ex}")
        return -1


def calculate_essential_for_given_k(k):
    sig = FindOpt()

    S = list(range(D.shape[0]))

    G = networkx.DiGraph()
    G.add_nodes_from(S)
    # Adding self edges
    G.add_edges_from((node, node) for node in G.nodes)
    for u in S:
        V = list(range(D.shape[0]))
        EssPairs(k, u, V, G, sig)

    return G, sig


def EssPairs(k, u, V, G, sig):
    if u in V:  # removing the self node from V
        V.remove(u)
    # we also need to remove vertices v from V, if we know that u <= v
    # basically , V - out(u)
    out_neighbors_u = list(G.successors(u))
    for out_neigh in out_neighbors_u:
        if out_neigh in V:
            V.remove(out_neigh)
    if len(V) > 0 and not test_ESS(k, [u], V, sig):
        if len(V) == 1:
            G.add_edge(u, V[0])
            # self edges are already added
            in_U = list(G.predecessors(u))
            out_V = list(G.successors(V[0]))
            for small_u in in_U:
                for small_v in out_V:
                    G.add_edge(small_u, small_v)
        else:
            (V_L, V_R) = split(V)
            EssPairs(k, u, V_L, G, sig)
            EssPairs(k, u, V_R, G, sig)


def test_ESS(k, U, V, sig):
    """
    This function takes k, U, V, sig and returns false (in feasible) if \\exists u \\in U: \\exists v \\in V: u \\leq_{e} v

    Parameters:
        k (int): Exact false positive bit flips, an integer.
        U (set): The set of cells to compare, a set.
        V (set): The set of cells to compare, a set.
        sig(float): The optimal value of the objective function.

    Returns:
        bool: returns false (infeasible) if \\exists u \\in U: \\exists v \\in V: u \\leq_{e} v, true only if \\forall u \\in U: \\forall v \\in V : u \\nleq_{e} v
    """
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
        if len(U) > 1 or len(V) > 1:
            print(f"Timeout set to {ilp_timeout} since size of V ({len(V)}) is greater than 1")
            model.Params.TimeLimit = ilp_timeout
        print(f"ILP timeout in the model {model.Params.TimeLimit}")
        time_olp_creation_s = time.time()

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
        glb_cons = model.addConstr(
            sum((1 - X[i, j]) * D[i, j] for i in range(n) for j in range(m) if D[i, j] == 1) == k,
            name="global_constraint")  # 8/5 test

        # model.update()
        if print_trace:
            print("printing global_constraint in an essential ILP call")
            print(f"{model.getRow(glb_cons)} {glb_cons.Sense} {glb_cons.RHS}")

        for u in U:
            nz = len(V)
            z = model.addMVar((m, nz), vtype=GRB.BINARY, name=f"z_{u}")
            # model.update()

            # Test for EO
            for i in range(m):
                for v_index in range(nz):
                    v = V[v_index]
                    model.addConstr(X[u, i] - X[v, i] <= z[i, v_index])  # (10)
                    model.addConstr(z[i, v_index] <= (X[u, i] - X[v, i] + 1) / 2)  # (10)

            for v in range(nz):
                model.addConstr(sum(z[i, v] for i in range(m)) >= 1)  # (11)

        model.addConstr(
            sum(sum((1 - D[i, j]) * X[i, j] for i in range(n) if D[i, j] != -1) for j in range(m)) == sig)  # (12) 8/5
        # model.addConstr(sum(sum(M[i]*(1 - D[i, j])*(X[i, j]) + k * M[i] * (D[i, j])*(1 - X[i, j]) for j in range(m)) for i in range(n)) == sig) # (8)

        # model.update()
        time_olp_creation_e = time.time()
        print(f"ILP creation time: {time_olp_creation_e - time_olp_creation_s}")

        current_time = time.strftime("%D:%H:%M:%S", time.localtime())
        print(current_time)

        time_ilp_start = time.time()
        model.optimize()
        time_ilp_end = time.time()
        current_time = time.strftime("%D:%H:%M:%S", time.localtime())
        print(f"Local time after finishing ILP Optimize: {current_time}")
        if print_trace:
            model.write("model.ess.lp")
        if model.Status == GRB.TIME_LIMIT:
            print(U, V, f" Time Limit Exceeded-----ILP time {time_ilp_end - time_ilp_start} seconds")
            return False
        if model.Status == 3:
            print(U, V, f" infeasible-----ILP time {time_ilp_end - time_ilp_start} seconds")
            return False
        else:
            print(U, V, f" feasible *****  -----ILP time {time_ilp_end - time_ilp_start} seconds")
            return True

    except GurobiError as ex:
        print(f"*********ERROR*********\n{ex}")


def test_ESS_with_mutation(k, U, V, sig, j):
    """
    This function takes k, U, V, sig, j and returns false (infeasible) if \\exists u \\in U: \\exists v \\in V: X_{uj} \\leq_{e} X_{vj}
    Parameters:
        k (int): Exact false positive bit flips, an integer.
        U (set): The set of cells to compare, a set.
        V (set): The set of cells to compare, a set.
        sig (float): The optimal value of the objective function.
        j (int): The mutation index
    Returns:
        bool: returns false (infeasible) if \\exists u \\in U: \\exists v \\in V: X_{uj} \\leq_{e} X_{vj}, true only if \\forall u \\in U: \\forall v \\in V : X_{uj} \\nleq_{e} X_{vj}
    """
    # if len(U) == 0 or len(V) == 0:
    #     return True
    # check whether the mutation index is out of bounds
    if not (0 <= j <= m):
        raise Exception("The mutation index is out of range")
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
        if len(U) > 1 or len(V) > 1:
            print(f"Timeout set to {ilp_timeout} since size of V ({len(V)}) is greater than 1")
            model.Params.TimeLimit = ilp_timeout
        print(f"ILP timeout in the model {model.Params.TimeLimit}")
        time_olp_creation_s = time.time()

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
        glb_cons = model.addConstr(
            sum((1 - X[i, j]) * D[i, j] for i in range(n) for j in range(m) if D[i, j] == 1) == k,
            name="global_constraint")  # 8/5 test

        # model.update()
        if print_trace:
            print("printing global_constraint in an essential ILP call")
            print(f"{model.getRow(glb_cons)} {glb_cons.Sense} {glb_cons.RHS}")

        for u in U:
            nz = len(V)
            z = model.addMVar((m, nz), vtype=GRB.BINARY, name=f"z_{u}")
            # model.update()

            # # Test for EO
            # for i in range(m):
            #     for v_index in range(nz):
            #         v = V[v_index]
            #         model.addConstr(X[u, i] - X[v, i] <= z[i, v_index])  # (10)
            #         model.addConstr(z[i, v_index] <= (X[u, i] - X[v, i] + 1) / 2)  # (10)

            # Add constraint for the specific mutation only.
            for v_index in range(nz):
                v = V[v_index]
                model.addConstr(X[u, j] - X[v, j] <= z[j, v_index])  # (10)
                model.addConstr(z[j, v_index] <= (X[u, j] - X[v, j] + 1) / 2)  # (10)

            # for v in range(nz):
            #     model.addConstr(sum(z[i, v] for i in range(m)) >= 1)  # (11)
            #
            # Add the constraint z_j >= 1 for each v \in V
            for v in range(nz):
                model.addConstr(z[j, v] >= 1)  # (11) for specific mutation only

        model.addConstr(
            sum(sum((1 - D[i, j]) * X[i, j] for i in range(n) if D[i, j] != -1) for j in range(m)) == sig)  # (12) 8/5
        # model.addConstr(sum(sum(M[i]*(1 - D[i, j])*(X[i, j]) + k * M[i] * (D[i, j])*(1 - X[i, j]) for j in range(m)) for i in range(n)) == sig) # (8)

        # model.update()
        time_olp_creation_e = time.time()
        print(f"ILP creation time: {time_olp_creation_e - time_olp_creation_s}")

        current_time = time.strftime("%D:%H:%M:%S", time.localtime())
        print(current_time)

        time_ilp_start = time.time()
        model.optimize()
        time_ilp_end = time.time()
        current_time = time.strftime("%D:%H:%M:%S", time.localtime())
        print(f"Local time after finishing ILP Optimize: {current_time}")
        if print_trace:
            model.write("model.ess.lp")
        if model.Status == GRB.TIME_LIMIT:
            print(U, V, f" Time Limit Exceeded-----ILP time {time_ilp_end - time_ilp_start} seconds")
            return False
        if model.Status == 3:
            print(U, V, f" infeasible-----ILP time {time_ilp_end - time_ilp_start} seconds")
            return False
        else:
            print(U, V, f" feasible *****  -----ILP time {time_ilp_end - time_ilp_start} seconds")
            return True

    except GurobiError as ex:
        print(f"*********ERROR*********\n{ex}")


def split(V):
    mid = len(V) // 2
    return V[:mid], V[mid:]


fileResults = outputName + "." + str(k) + ".esspairs" + ".txt"
print(f"The essential pairs will be written to {fileResults}")
pathToResult = f"{args.result_folder}/{outputName}/"
os.makedirs(pathToResult, exist_ok=True)

Graph, the_opt = calculate_essential_for_given_k(k)

# let's do a small test here.
# we know that 0 <= 10
for mutation in range(m):
    print(f"Mutation {mutation}: The feasibility {test_ESS_with_mutation(k, [5], [6], the_opt, mutation)}")

end_time = time.time()
networkx.write_edgelist(Graph, f"{pathToResult}{fileResults}")  # writing essential relations to file.
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
    # f.write(f"Poset Width: {Width(Graph)}\n")
    f.write(f"Essential Relation: {[edge for edge in Graph.edges]}\n\n")
    print(f"\nResults written was to file: {verboseResultFile}")
    f.close()
