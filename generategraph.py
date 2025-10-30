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
                   result_folder="results_default"):
    dot = graphviz.Digraph(f"{output_name}_k{kappa}",
                           node_attr={'style': 'filled',
                                      'shape': 'circle'},
                           edge_attr={'arrowsize': '0.3'})
    # Specify the resolution of the graph
    dot.attr(dpi='1000')

    min = math.inf
    max = -math.inf
    for node in g.nodes:
        size = len(node.split(','))
        if size < min:
            min = size
        if size > max:
            max = size

    for node in g.nodes():
        size = len(node.split(','))
        intensity = parameter_to_intensity(size, max, min, 0, 0.5)
        new_color = adjust_color_intensity(node_fill_color, intensity)
        dot.node(node, label=f"", style='filled', fillcolor=new_color)
    for edge in g.edges():
        dot.edge(str(edge[0]), str(edge[1]))

    dot.render(f"{result_folder}/{output_name}/{output_name}_kappa_{kappa}", format='pdf', cleanup=True)
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


def main():
    parser = argparse.ArgumentParser(description="Arguments for the dot graph creation program")
    parser.add_argument('inputFile', type=str, help="Input the final graph with mutation labels")
    parser.add_argument('kappa', type=int, default=1, help="kappa value for the input program")
    parser.add_argument('-node_fill_color', type=str, default="none")
    parser.add_argument('-resultsFolder', type=str, default="results_default")
    parser.add_argument('--verbose', action='store_true', help='Increase output verbosity')
    parser.add_argument('-change_border_size', action='store_true', help='This would ensure that resultant graph '
                                                                         'changes intensity of the borders based on ')
    parser.add_argument('-change_node_fill_color_intensity', action='store_true', help='This would ensure that '
                                                                                       'resultant graph changes the '
                                                                                       'intensity of node fill color '
                                                                                       'based on the node size')
    args = parser.parse_args()

    print(f"Input file: {args.inputFile}")
    print(f"Value of kappa: {args.kappa}")
    print(f"Node fill color: {args.node_fill_color}")
    print(f"Custom results folder {args.resultsFolder}")
    print(f"Verbosity: {args.verbose}")
    print(f"change_border_size: {args.change_border_size}")
    print(f"change_node_fill_color_intensity: {args.change_node_fill_color_intensity}")
    filenamewithoutextension = args.inputFile[:-4]
    k = args.kappa
    result_path = f"{args.resultsFolder}/{filenamewithoutextension}/"
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

    generate_graph(g, filenamewithoutextension, k, "black", args.node_fill_color, args.resultsFolder)

if __name__ == "__main__":
    main()
