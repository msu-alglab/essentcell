"""
@author: Adiesha Liyanage and Brendan Mumey and Braeden Sopp
"""
import argparse
import os
import time

import networkx
import pandas as pd
from gurobipy import *
import itertools


def increment_counter(parameters: dict, key: str):
    parameters[key] = parameters[key] + 1


def split(V):
    mid = len(V) // 2
    return V[:mid], V[mid:]


def FindOpt(parameters: dict, k: int) -> float:
    """
    Returns sigma, the minimum number of bit flips required to make D conflict-free.
    """
    n = parameters["numofrows"]
    m = parameters["numofmutations"]
    D = parameters["data"]
    print_trace = parameters["print_trace"]
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
        # verbose = parameters["verbose"]
        # if verbose:
        #     model.write("model.initial.lp")
        sig = model.ObjVal
        print(f"sig: {sig}")
        return sig

    except GurobiError as ex:
        print(f"*********ERROR*********\n{ex}")
        return -1


def test_ESS(parameters: dict, k: int, U, V, sig: float):
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
    D = parameters["data"]
    increment_counter(parameters, "count")
    ilp_timeout = parameters["ilp_timeout"]
    n = parameters["numofrows"]
    m = parameters["numofmutations"]
    print_trace = parameters["print_trace"]
    if parameters["count"] % 50 == 0:
        print(f"Essential relation ILP call {parameters['count']} times.")
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


def EssPairs(parameters: dict, k: int, u, V, G, sig: float):
    if u in V:  # removing the self node from V
        V.remove(u)
    # we also need to remove vertices v from V, if we know that u <= v
    # basically , V - out(u)
    out_neighbors_u = list(G.successors(u))
    for out_neigh in out_neighbors_u:
        if out_neigh in V:
            V.remove(out_neigh)
    if len(V) > 0 and not test_ESS(parameters, k, [u], V, sig):
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
            EssPairs(parameters, k, u, V_L, G, sig)
            EssPairs(parameters, k, u, V_R, G, sig)


def calculate_essential_for_given_k(k: int, parameters: dict):
    D = parameters["data"]  # get the data matrix
    sig = FindOpt(parameters, k)

    S = list(range(D.shape[0]))

    G = networkx.DiGraph()
    G.add_nodes_from(S)
    # Adding self edges
    G.add_edges_from((node, node) for node in G.nodes)
    disable_gt = parameters["disable_gt"]
    if not disable_gt:
        for u in S:
            V = list(range(D.shape[0]))
            EssPairs(parameters, k, u, V, G, sig)
    else:
        comb_n_c_2 = itertools.combinations(S, 2)
        for (u, v) in comb_n_c_2:
            EssPairs(parameters, k, u, [v], G, sig)
            EssPairs(parameters, k, v, [u], G, sig)

    return G, sig


def test_ESS_with_mutation(parameters: dict, k: int, U, V, sig: float, j: int):
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
    D = parameters["data"]

    ilp_timeout = parameters["ilp_timeout"]
    n = parameters["numofrows"]
    m = parameters["numofmutations"]
    print_trace = parameters["print_trace"]
    # check whether the mutation index is out of bounds
    if not (0 <= j <= m):
        raise Exception("The mutation index is out of range")
    increment_counter(parameters, "count")
    if parameters["count"] % 50 == 0:
        print(f"Essential relation ILP call {parameters['count']} times.")
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


def Width(G):
    max = 0
    for x in networkx.antichains(G):
        if len(x) > max:
            max = len(x)
    return max


def add_info_about_final_graph(g, output_name, kappa, results_folder="results"):
    path_to_result = f"{results_folder}/{output_name}/"
    verbose_result_file = f"{path_to_result}{output_name}_kappa_{kappa}.graph_info.txt"  # creates the file automatically if it is not there
    print(f"Verbose result is written into {verbose_result_file}")
    f = open(verbose_result_file, "a")
    f.write(f"Number of Nodes: {g.number_of_nodes()}\n")
    f.write(f"Poset Width: {Width(g)}\n")
    f.write(f"Essential Relation: {[edge for edge in g.edges]}\n\n")
    print(f"\nResults written was to file: {verbose_result_file}")
    f.close()
    # persisting the intersection graph with collapsed with mutation labels
    intersection_graph_file = f"{path_to_result}{output_name}_kappa_{kappa}.graph_persist.txt"
    networkx.write_edgelist(g, intersection_graph_file, delimiter="#")  # writing essential relations to file.


def find_the_mutation_labels(parameters: dict, g: networkx.DiGraph):
    kmin = parameters["kmin"]
    kmax = parameters["kmax"]
    verbose = parameters["verbose"]
    n = parameters["numofrows"]
    m = parameters["numofmutations"]

    # basically each edge in the hasse diagram we need to find mutations.
    # first find the optimal value
    for (u, v) in g.edges():  # for each edge
        # we need to pick a representative from both u and v
        u_rep = int(u.split(",")[0])
        v_rep = int(v.split(",")[0])
        print(f"u: {u_rep}, v: {v_rep}")
        labels = []
        for mutation in range(0, m):  # for each mutation
            for k in range(kmin, kmax + 1):
                sig_k = FindOpt(parameters, k)
                # if v->u is feasible for mutation, add it to labels
                feasibility = test_ESS_with_mutation(parameters, k, [v_rep], [u_rep], sig_k, mutation)
                if verbose:
                    if feasibility:
                        print(f"Mutation {mutation} is feasible for v:{v} to u:{u}")
                        labels.append(mutation)
                        break
                    else:
                        print(f"Mutation {mutation} is NOT feasible for v:{v} to u:{u}")
        g[u][v]["mutation_labels"] = str(labels)

    for u, v in g.edges():
        print(f"Edge {u}->{v} with weight {g[u][v]['mutation_labels']}")


def main():
    parser = argparse.ArgumentParser(description="Arguments for the EssentCell program")
    parser.add_argument('inputFile', type=str, help="Sorted Input file to the program")
    parser.add_argument('kmin', type=int, default=0, help="k min value for the input program")
    parser.add_argument('kmax', type=int, default=0, help="k max value for the input program")
    parser.add_argument('-result_folder', type=str, default="results_default", help="The result folder name")
    parser.add_argument('--verbose', action='store_true', help='Increase output verbosity')
    parser.add_argument('-print_trace_of_constraint', action='store_true', help='Increase the output verbosity of '
                                                                                'constraints')
    parser.add_argument('-timeout', type=float, default=float('inf'), help="Timeout value for ILP calls")
    parser.add_argument('-disable_gt', action='store_true', help="diable group testing")
    args = parser.parse_args()

    print(f"Input file: {args.inputFile}")
    print(f"Value of k min: {args.kmin}")
    print(f"Value of k max: {args.kmax}")
    print(f"Result folder name {args.result_folder}")
    print(f"Verbosity: {args.verbose}")
    print(f"Print Trace Verbosity: {args.print_trace_of_constraint}")
    print(f"Group testing disabled: {args.disable_gt}")
    print(f"ILP timeout: {args.timeout}")

    fileName = args.inputFile
    outputName = fileName[:-4]  # remove the file extension
    print(f"Input file: {fileName}")
    print(f"Output file: {outputName}")
    results_folder = args.result_folder

    # defining the necessary parameters for the whole run.
    kmin = args.kmin
    kmax = args.kmax
    verbose = args.verbose
    print_trace = args.print_trace_of_constraint
    ilp_timeout = args.timeout

    # First read the file and put the data into a dataframe
    df = pd.read_csv(fileName)
    df.index = range(1, len(df) + 1)

    D = df.to_numpy()  # convert into numpy
    n = D.shape[0]  # number of cells
    m = D.shape[1]  # number of mutations

    parameters = {"fileName": fileName, "kmin": kmin, "kmax": kmax, "verbose": verbose, "ilp_timeout": ilp_timeout,
                  "print_trace": print_trace, "data": D, "numofrows": n, "numofmutations": m, "dataframe": df,
                  "disable_gt": args.disable_gt}

    # now we need to run the program for k = kmin to kmax and generate the graphs.
    Graphs = []
    for k in range(kmin, kmax + 1):
        print(f"Value of k is : {k} and {df.head()}")
        print(f"Number of rows {parameters['numofrows']}")
        print(f"Number of rows {parameters['numofmutations']}")
        parameters['count'] = 0
        print(f"Setting the count variable to {parameters['count']}")
        fileResults = outputName + "." + str(k) + ".esspairs" + ".txt"
        print(f"The essential pairs will be written to {fileResults}")
        pathToResult = f"{args.result_folder}/{outputName}/"
        os.makedirs(pathToResult, exist_ok=True)

        start_time = time.time()  # start time of the program.
        Graph_k, the_opt = calculate_essential_for_given_k(k, parameters)
        Graphs.append(Graph_k)
        end_time = time.time()
        networkx.write_edgelist(Graph_k, f"{pathToResult}{fileResults}")  # writing essential relations to file.
        if verbose:
            verboseResultFile = f"{pathToResult}{outputName}.{k}.esspairs.verbose.txt"
            print(f"Verbose result is written into {verboseResultFile}")
            f = open(verboseResultFile, "a")
            f.write(f"k value: {k}\n")
            f.write(f"n (number of samples): {n}\n")
            f.write(f"m (number of mutations): {m}\n")
            f.write(f"EssILP calls: {parameters['count']}\n")
            f.write(f"Runtime: {end_time - start_time} seconds\n")
            f.write(f"Number of Nodes: {Graph_k.number_of_nodes()}\n")
            # f.write(f"Poset Width: {Width(Graph)}\n")
            f.write(f"Essential Relation: {[edge for edge in Graph_k.edges]}\n\n")
            print(f"\nResults written was to file: {verboseResultFile}")
            f.close()

    print("Step 2---------------------------------------------------------")
    print("Working on the strongly connected graph")
    print("Step 2.1--- Generating the intersection graph")

    intersection_of_edges = set(Graphs[0].edges)
    print(Graphs)
    print(Graphs[0])
    print(intersection_of_edges)
    i = 0
    while i <= (kmax - kmin) and len(intersection_of_edges) != 0:
        g_i = Graphs[i]
        intersection_of_edges = intersection_of_edges & set(g_i.edges)
        i += 1

    print(f"The final intersection of edges {intersection_of_edges}")
    print("Cleaning up the edge set to create the final graph")
    ess_list = list(intersection_of_edges)
    g = networkx.DiGraph(ess_list)

    print("Looking for strongly connected components")
    for scc in networkx.strongly_connected_components(g):
        l1 = [num for num in scc]
        l2 = [str(scc)[1:len(str(scc)) - 1] for i in range(len(scc))]  # removing the curly brackets
        mapping = dict(zip(l1,
                           l2))  # here we use zip to label the nodes in the same scc to have labels of all the nodes in the scc
        g = networkx.relabel_nodes(g,
                                   mapping)  # once you call the relabel_nodes, nodes with same label in the new mapping collapses
        # note that this automatically updates any old edges that we had between nodes with in the scc and going in/out to them

    # Enforce that there are no self-loop (reflexive) edges in the graph
    print("Removing self edges")
    g.remove_edges_from(networkx.selfloop_edges(g))

    # Enforce transitive property of the partial order relation
    print("Removing transitive edges")
    g = networkx.transitive_reduction(g)

    print(g.edges)

    find_the_mutation_labels(parameters, g)

    add_info_about_final_graph(g, outputName, kmax, results_folder)


if __name__ == '__main__':
    main()
