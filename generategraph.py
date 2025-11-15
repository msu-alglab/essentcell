"""
@author: Adiesha Liyanage and Brendan Mumey and Braeden Sopp
"""
import argparse
import math
import os.path
import graphviz
import networkx
import matplotlib.colors as mcolors
import colorsys


def generate_graph(g, output_name, kappa, node_border_color="black", node_fill_color="none",
                   result_folder="results_default", mutations=None, cluster_id_prefix='c_'):
    dot = graphviz.Digraph(f"{output_name}_k{kappa}",
                           node_attr={'style': 'filled',
                                      'shape': 'circle'},
                           edge_attr={'arrowsize': '0.3', 'fontname': 'DejaVu Sans'})
    # Specify the resolution of the graph
    dot.attr(dpi='1000')

    # following determines the max and min size of the nodes in case we need to change the size of the nodes.
    min = math.inf
    max = -math.inf
    for node in g.nodes:
        size = len(node.split(','))
        if size < min:
            min = size
        if size > max:
            max = size

    cluster_id_dict = {}
    i = 0
    topo_sort = networkx.topological_sort(g)
    for node in topo_sort: # here I iterate over topo_sort instead of g.nodes()
        size = len(node.split(','))
        # this is code to adjust the color
        # intensity = parameter_to_intensity(size, max, min, 0, 0.5)
        # new_color = adjust_color_intensity(node_fill_color, intensity)
        cluster_id_dict[f'{cluster_id_prefix}{i}'] = node
        # this is the code for setting the size of the node as well
        # max_node_size_inch = 2.5
        # min_node_size_inch = 0.5
        # adjusted_size = (max_node_size_inch - min_node_size_inch) * (size-min_node_size_inch)/(max-min)
        # dot.node(node, label=f'cl_id_{i}', style='filled', fillcolor=new_color, width=str(adjusted_size), height=str(adjusted_size))
        # dot.node(node, label=f'cl_id_{i}', style='filled', fillcolor=new_color)
        dot.node(node, label=f'{cluster_id_prefix}{i}\n{size}', style='filled', fillcolor=node_fill_color)
        i += 1
    for u, v, data in g.edges(data=True):
        mutation_list_str = data['weight']
        mutation_list_str = mutation_list_str.replace('[', '').replace(']', '').replace(' ', '')
        mutation_list = mutation_list_str.split(',')
        new_mutation_list = [mutations[mute] for mute in mutation_list]
        new_mutation_list_str = '\u2228'.join(
            new_mutation_list)  # note that font is important if you want to display the special character
        dot.edge(str(u), str(v), label=new_mutation_list_str)

    print(cluster_id_dict)
    dot.render(f"{result_folder}/{output_name}/{output_name}_kappa_{kappa}", format='pdf', cleanup=True)
    cluster_file_name = f"{result_folder}/{output_name}/{output_name}_kappa_{kappa}.graph_cluster_id.txt"
    with open(cluster_file_name, "w") as f:
        f.write(str(cluster_id_dict))
    return g


def parameter_to_intensity(val, max, min, target_min=0.0, target_max=1.0):
    """ Convert a parameter (0 to 1) to a grayscale color. """
    value = 0 + target_min + ((val - min) / (max - min)) * (target_max - target_min)
    return value


def adjust_color_intensity(color, intensity):
    if intensity > 1 or intensity < 0:
        print(f"Intensity is not between 0 and 1..manually making it between 0 and 1")
        intensity = max(0, min(intensity, 1))
    print(mcolors.to_rgb(color))
    (h, l, s) = colorsys.rgb_to_hls(*mcolors.to_rgb(color))
    (r, g, b) = colorsys.hls_to_rgb(h, l * (1 - intensity), s)
    new_hex_color = mcolors.to_hex((r, g, b))
    return new_hex_color


def read_mutation_list(filepath):
    mutations = {}
    try:
        with open(filepath, "r") as f:
            for line in f:
                line = line.rstrip()
                split_line = line.split(" ")
                mutations[split_line[0]] = split_line[1]
            return mutations

    except FileNotFoundError:
        print(f"File not found: filepath: {filepath}")
    except Exception as e:
        print(f"Error occurred: {e}")


def main():
    parser = argparse.ArgumentParser(description="Arguments for the dot graph creation program")
    parser.add_argument('inputFile', type=str, help="Input the final graph with mutation labels")
    parser.add_argument('kappa', type=int, default=1, help="kappa value for the input program")
    parser.add_argument('-node_fill_color', type=str, default="none")
    parser.add_argument('-result_folder', type=str, default="results_default")
    parser.add_argument('--verbose', action='store_true', help='Increase output verbosity')
    parser.add_argument('-change_border_size', action='store_true', help='This would ensure that resultant graph '
                                                                         'changes intensity of the borders based on ')
    parser.add_argument('-change_node_fill_color_intensity', action='store_true', help='This would ensure that '
                                                                                       'resultant graph changes the '
                                                                                       'intensity of node fill color '
                                                                                       'based on the node size')
    parser.add_argument('-min_node_size', type=int, default=0, help="Minimum node size to be included in the graph")
    parser.add_argument('-cluster_prefix', type=str, default='c_',
                        help="Cluster prefix to be used when generating the graph")
    args = parser.parse_args()

    print(f"Input file: {args.inputFile}")
    print(f"Value of kappa: {args.kappa}")
    print(f"Node fill color: {args.node_fill_color}")
    print(f"Custom results folder {args.result_folder}")
    print(f"Verbosity: {args.verbose}")
    print(f"change_border_size: {args.change_border_size}")
    print(f"change_node_fill_color_intensity: {args.change_node_fill_color_intensity}")
    filenamewithoutextension = args.inputFile[:-4]
    k = args.kappa
    min_node_size = args.min_node_size
    print(f'Min node size: {min_node_size}')
    result_path = f"{args.result_folder}/{filenamewithoutextension}/"
    print(
        f"Pattern of files to be read {result_path}{filenamewithoutextension}_kappa_{k}.graph_persist.txt from in which k is from {0} to {args.kappa}")
    print("First checking the availability of all files")
    open_file_path = f"{result_path}{filenamewithoutextension}_kappa_{k}.graph_persist.txt"
    g = networkx.DiGraph()
    try:
        with open(open_file_path, "r") as f:
            for line in f:
                words = line.split(sep="#")
                # print(words)
                g.add_edge(words[0], words[1], weight=words[2][21:-3])
            print(f"Finished reading {open_file_path}")
    except FileNotFoundError:
        print(f"File not found: filepath: {open_file_path}")
    except Exception as e:
        print(f"Error occurred: {e}")
    print(f"Number of nodes is: {g.number_of_nodes()}")
    print(f"Number of edges is: {g.number_of_edges()}")

    mutations = read_mutation_list(f'{result_path}mutations.txt')
    # mutations = {"0": "haha0", "1": "haha1", "2": "haha2", "3": "haha3", "4": "haha4", "5": "haha5", "6": "haha6"}

    # remove nodes that is smaller than the min_node size
    # in networkx, when you remove a node, it automatically removed edges associated with it.
    nodes_to_remove = [node for node in g.nodes() if len(node.split(',')) < min_node_size]
    g.remove_nodes_from(nodes_to_remove)

    generate_graph(g, filenamewithoutextension, k, "black", args.node_fill_color, args.result_folder, mutations)


if __name__ == "__main__":
    main()
