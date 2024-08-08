"""
Author: Robyn Burger and Allison Shi and Adiesha Liyanage
Date: 2024-08-07
Description: Given a inputFileName, generateIntersectionGraph.py looks at already containing data about the essential
order for input kappa value. Then it produces the intersection graph and output the graph and its info.

Example Usage:

$ python generateIntersectionGraph.py -h

$ python generateIntersectionGraph.py data.sorted.csv 2 -node_fill_color red --verbose

$ python generateIntersectionGraph.py data.sorted.csv 3
"""

import argparse
import os.path
import graphviz
import networkx


def generate_graph(g, output_name, kappa, node_border_color="black", node_fill_color="none", result_folder="results"):
    dot = graphviz.Digraph(f"{output_name}_k{kappa}",
                           node_attr={'color': node_border_color, 'style': 'filled', 'fillcolor': node_fill_color,
                                      'shape': 'circle'},
                           edge_attr={'arrowsize': '0.3'})
    # Specify the resolution of the graph
    dot.attr(dpi='1000')

    for node in g.nodes():
        dot.node(node, label=None)
    for edge in g.edges():
        dot.edge(str(edge[0]), str(edge[1]))

    dot.render(f"{result_folder}/{output_name}/{output_name}_kappa_{kappa}", format='pdf', cleanup=True)
    return g


def Width(G):
    max = 0
    for x in networkx.antichains(G):
        if len(x) > max:
            max = len(x)
    return max


def add_info_about_final_graph(g, output_name, kappa, results_folder="results"):
    path_to_result = f"{results_folder}/{output_name}/"
    verbose_result_file = f"{path_to_result}{output_name}_kappa_{kappa}.graph_info.txt"
    print(f"Verbose result is written into {verbose_result_file}")
    f = open(verbose_result_file, "a")
    f.write(f"Number of Nodes: {g.number_of_nodes()}\n")
    f.write(f"Poset Width: {Width(g)}\n")
    f.write(f"Essential Relation: {[edge for edge in g.edges]}\n\n")
    print(f"\nResults written was to file: {verbose_result_file}")
    f.close()


def main():
    parser = argparse.ArgumentParser(description="Arguments for the IntersectionGraph program")
    parser.add_argument('inputFile', type=str, help="Sorted Input file that was used to calculate the ess pairs")
    parser.add_argument('kappa', type=int, default=1, help="kappa value for the input program")
    parser.add_argument('-node_fill_color', type=str, default="none")
    parser.add_argument('-resultsFolder', type=str, default="results")
    parser.add_argument('--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    print(f"Input file: {args.inputFile}")
    print(f"Value of kappa: {args.kappa}")
    print(f"Node fill color: {args.node_fill_color}")
    print(f"Custom results folder {args.resultsFolder}")
    print(f"Verbosity: {args.verbose}")
    result_path = f"{args.resultsFolder}/{args.inputFile}/"
    print(
        f"Pattern of files to be read {result_path}{args.inputFile}.k.esspairs.txt from in which k is from {0} to {args.kappa}")
    print("First checking the availability of all files")
    for i in range(args.kappa):
        file_name = f"{result_path}{args.inputFile}.{i}.esspairs.txt"
        if os.path.exists(file_name) and os.path.isfile(file_name):
            print(f"File {file_name} exists")
        else:
            print(f"File {file_name} DOES NOT EXIST ")

    print("Generating intersection graph")
    g_0 = networkx.read_edgelist(f"{result_path}{args.inputFile}.{0}.esspairs.txt", nodetype=int,
                                 create_using=networkx.DiGraph)
    intersection_of_edges = set(g_0.edges)
    print(g_0)
    print(intersection_of_edges)
    i = 1
    while i <= args.kappa and len(intersection_of_edges) != 0:
        # read the file
        file_name = f"{result_path}{args.inputFile}.{i}.esspairs.txt"
        if os.path.exists(file_name) and os.path.isfile(file_name):
            g_i = networkx.read_edgelist(file_name, nodetype=int, create_using=networkx.DiGraph)
            intersection_of_edges = intersection_of_edges & set(g_i.edges)
        else:
            print(f"File {file_name} DOES NOT EXIST -- Exiting the program ")
            exit(1)
        i += 1

    print(f"The final intersection of edges {intersection_of_edges}")
    print("Cleaning up the edge set to create the final graph")
    ess_list = list(intersection_of_edges)
    g = networkx.DiGraph(ess_list)
    print("Looking for strongly connected components")
    for scc in networkx.strongly_connected_components(g):
        l1 = [num for num in scc]
        l2 = [str(scc)[1:len(str(scc)) - 1] for i in range(len(scc))]
        mapping = dict(zip(l1, l2))
        g = networkx.relabel_nodes(g, mapping)

    # Enforce that there are no self-loop (reflexive) edges in the graph
    print("Removing self edges")
    g.remove_edges_from(networkx.selfloop_edges(g))

    # Enforce transitive property of the partial order relation
    print("Removing transitive edges")
    g = networkx.transitive_reduction(g)

    print(g.edges)

    generate_graph(g, args.inputFile, args.kappa, "black", args.node_fill_color, args.resultsFolder)

    add_info_about_final_graph(g, args.inputFile, args.kappa, args.resultsFolder)


if __name__ == '__main__':
    main()
