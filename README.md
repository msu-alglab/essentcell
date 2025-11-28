# EssentCell

Cancer is an evolutionary disease. As cells multiply, they can take on new mutations, and through the accumulation of those mutations, they can become cancerous. Tracking this evolution is helpful for understanding how the disease progresses and the most effective ways to treat it. We model this evolution in diagrams called phylogenetic trees. Much like a family tree, a phylogenetic tree depicts the relationships between generations of cells. These models depend on taking physical samples from cancer patients.

But, these samples are imperfect, often leading to multiple phylogenetic trees. EssentCell finds the lineages that stay consistent among all optimal phylogenetic trees. These essential relations provide insight into the most probable ancestral relationships, painting a more precise picture of a tumor's evolutionary history.

## Content
  1. [Getting started](#start)
     * [Dependencies](#dep)
  2. [Usage instructions](#usage)
     * [I/O formats](#io)
     * [Examples](#example)


<a name="start"></a>
## Getting started

EssentCell is implemented using python.

<a name="dep"></a>

### Dependencies  

EssentCell has the following dependencies
 * Gurobi 11.0.2
 * NumPy 1.26.3
 * Pandas 2.2.0
 * NetworkX 3.2.1

<a name="usage"></a>
## Usage Instructions

<a name="io"></a>
### I/O formats
The input to EssentCell is a .csv file that contains n rows and m columns where n is the the number of single-cells and m is the number of mutations.  All entries in the .csv file should be either 1 if mutation j is present in cell i or 0 if mutation j is not present in cell i. Missing entries should be indicated by -1.  

```
usage: EssentCell.py [-h] [-result_folder RESULT_FOLDER] [--verbose] [-print_trace_of_constraint] [-timeout TIMEOUT] [-disable_gt] inputFile kmin kmax

Arguments for the EssentCell program

positional arguments:
  inputFile             Input file to the program; the file should be on the same folder that the EssentCell.py script resides.
  kmin                  k min value for the input program; this is the minimum fixed 1 -> 0 bit flips that will be used by the program. 
  kmax                  k max value for the input program; this is the maximum fixed 1 -> 0 bit flips that will be used by the program. Note that kmin <= kmax

optional arguments:
  -h, --help            show this help message and exit
  -result_folder RESULT_FOLDER
                        The result folder name; by default result folder is "results_default". Inside the result_folder, program will create another folder using the input file name excluding .csv
                        Output files corresponding to the input file will reside in this folder. 
  --verbose             Increase output verbosity; This will ensure that the number of ILP calls and runtime are persist in an output file.
  -print_trace_of_constraint
                        Increase the output verbosity of gurobi ILP creations and prints out the constraints created by the ILP calls
  -timeout TIMEOUT      Timeout value for ILP calls; this should be given in seconds. Note that the timeout value is only considered for ILP calls that ivolves groups larger than 2. When we run a ILP calls between two cells, the timeout value will be ignored as we need to ensure that the ILP call is either feasible or infeasible. By default timeout value is infinity.
  -disable_gt           diable group testing; if the user desires not to use group testing, then enable this flag.

```
<a name="example"></a>
### Examples

Following is an example on how to run the script using **smallest.sorted.csv** as the input file. The user must run the EssentCell.py file by going into the folder that the **EssentCell.py** resides in. The input files to this script must reside in the same folder as well.

```
python EssentCell.py smallest.sorted.csv 0 2 --verbose -timeout 300
```

This command will input the **smallest.sorted.csv** file that is residing in the root folder and run EssentCell program for k=0, k=1, and k=2 values. Then the program will create the following files.
- k = 0
  * [smallest.sorted.0.esspairs.txt](results_default/smallest.sorted/smallest.sorted.0.esspairs.txt)
  * [smallest.sorted.0.esspairs.verbose.txt](results_default/smallest.sorted/smallest.sorted.0.esspairs.verbose.txt)
- k = 1
  * [smallest.sorted.1.esspairs.txt](results_default/smallest.sorted/smallest.sorted.1.esspairs.txt)
  * [smallest.sorted.1.esspairs.verbose.txt](results_default/smallest.sorted/smallest.sorted.1.esspairs.verbose.txt)
- k = 2
  * [smallest.sorted.2.esspairs.txt](results_default/smallest.sorted/smallest.sorted.1.esspairs.txt)
  * [smallest.sorted.2.esspairs.verbose.txt](results_default/smallest.sorted/smallest.sorted.1.esspairs.verbose.txt)
- [smallest.sorted_kappa_2.graph_info.txt](results_default/smallest.sorted/smallest.sorted_kappa_2.graph_info.txt)
- [smallest.sorted_kappa_2.graph_persist.txt](results_default/smallest.sorted/smallest.sorted_kappa_2.graph_persist.txt)

For each k value, **smallest.sorted.k.esspairs.txt** file contains the essential relation graph as set of edge list (before collapsing strongly connected components). The **smallest.sorted.k.esspairs.verbose.txt** files contain the extra information such as how many ILP calls were called and the total runtime for each particular **k** value along with other information about the graph. Please take a look at the example output and actual result outputs for more details.

Next, smallest.sorted_kappa_2.graph_info.txt contains intersection graph, which contains the graph that has edges appearing in all essential relation graphs from k=min to k=kmax. Note that we further process this graph before persisting, i.e., we collapse the strongly connected components and perform transitive reduction.
Finally, smallest.sorted_kappa_2.graph_persist.txt contains the same graph but with mutation labels.

Note that if you run this command multiple times **smallest.sorted.k.esspairs.txt** files will be created from scratch while the data to the **smallest.sorted.k.esspairs.verbose.txt** files **will be appended at the end.** As you can see the verbose files contain records of multiple runs. If the user does not wish to see this behaviour, they can delete the verbose files before running the command again.

Another example:

```
python EssentCell.py Patient2.csv  0  2 -result_folder Results/Results_"$current_month"_"$current_day" --verbose -timeout 300
```

User can structure the result folder by passing variables, when calling the script inside shell script.

Example that uses an input with missing data:

```
python EssentCell.py smallest.sorted_with_missing.csv  0  2 --verbose
```

We have included a example input file called **smallest.sorted_with_missing.csv** with missing data entries. Note that this is a random file generated with missing entries. Check the corresponding output folder for the input file.

If the user wishes to not use group testing for computing the essential relation, user can pass a special flag to the script and the program will use the naive approach to compute the essential relation. Following is an example usage on how to pass that flag. Here we put the result in a new folder called **result_wo_gt**.

```
python smallest.sorted.csv 0 2 -result_folder result_wo_gt --verbose -disable_gt
```

Refer the outout folder result_wo_gt for the output files.

#### Example usage for generating the final essential relation graph

The users can use **generategraph.py** to create essential relation graph as well. In order to use this script first, use the previous script to generate the necessary output files. Then create a new file called mutations.txt that contains the mutation mutation labels. In this file each row should contain a column id and mutation label seperated by blank space. The program will read this file to output the necessary edge label values with mutation names. Please check the output folders for example **mutations.txt** file. **This file must reside in the output folder that corresponds to its input csv file.**


```
usage: generategraph.py [-h] [-node_fill_color NODE_FILL_COLOR] [-result_folder RESULT_FOLDER] [--verbose] [-min_node_size MIN_NODE_SIZE]
                        [-cluster_prefix CLUSTER_PREFIX] [-do_not_keep_connections_when_deleting] [-max_number_of_mutation_labels MAX_NUMBER_OF_MUTATION_LABELS]
                        inputFile kappa

Arguments for the dot graph creation program

positional arguments:
  inputFile             Input the final graph with mutation labels; the user should enter the input file name along with .csv extension.
  kappa                 kappa value for the input program

optional arguments:
  -h, --help            show this help message and exit
  -node_fill_color NODE_FILL_COLOR
  -result_folder RESULT_FOLDER
                        This arguement should be pointing to the result folder where your output files are residing. Inside this folder there should be another folder named using the input file name excluding the                         .csv extension. By default this folder is results_default.
  -min_node_size MIN_NODE_SIZE
                        Minimum node size to be included in the graph. This argument looks at nodes that corresponds to clusters of size less than this threshold and remove it from the graph.
                        This arguement is added to improve the clarity of the final diagram. If the user wishes to keep all nodes, then ignore this argument.
  -cluster_prefix CLUSTER_PREFIX
                        Cluster prefix to be used when generating the graph; This prefix will be used to name the names of the nodes, and cluster id will be added after cluster prefix.
  -do_not_keep_connections_when_deleting
                        This would ensure that we keep the connections of edges when we delete nodes that are too small; If the user do not want to keep the connections that goes through nodes that are removed                            using earlier arguement, then make sure this flag is enabled. By default, the program keeps the connections that goes through deleted nodes. This option is here to improve the clarity of                           the final figure.
  -max_number_of_mutation_labels MAX_NUMBER_OF_MUTATION_LABELS
                        Maximum number of mutation labels in the final graph. Default value is to include all mutation labels. Any edge label that contains more than max number of edge labels will not be                                  displayed in the final graph. This option is there to improve the clarity of the figure.

```

Example usage:

```
python generategraph.py Patient2.csv 2 --verbose -node_fill_color None -min_node_size 2 -max_number_of_mutation_labels 1
```

Note that **[mutation.txt](results_default/Patient2/mutations.txt)** file for the Patient2.csv resides in results_default/Patient2 folder. 

Two files will be created by this script.
- Patient2_kappa_2.graph_cluster_id.txt
    * This program will generate an additional file named **smallest.sorted_kappa_2.graph_cluster_id.txt** that contains the information about the cluster ids and cells allocated to each of these clusters.
    * [Refer the following file](https://github.com/msu-alglab/essentcell/blob/main/results_default/Patient2/Patient2_kappa_2.graph_cluster_id.txt)
- Patient2_kappa_2.pdf
    * contains the final output graph.
    * [Refer the following file](results_default/Patient2/Patient2_kappa_2.pdf)

If the user wishes to **customize the final graph**, then the user can edit the **generategraph.py** using graphviz attributes to generate the desired output graph.

Another example:

```
python generategraph.py smallest.sorted.csv 2 --verbose -result_folder result_wo_gt -node_fill_color None -min_node_size 0 -max_number_of_mutation_labels 6
```
In this example, we have inputted "result_wo_gt" as the input and output folder for the generategraph.py. Note that **[mutations.txt](result_wo_gt/smallest.sorted/mutations.txt)** exists in the corresponding subfolder in **result_wo_gt**.
