<p align="center">
  <img src="logo/logo.png" alt="EssentCell" width="200"/>
</p>
<h1 align="center">EssentCell</h1>


Cancer is an evolutionary disease. As cells multiply, they can take on new mutations, and through the accumulation of those mutations, they can become cancerous. Tracking this evolution is helpful for understanding how the disease progresses and the most effective ways to treat it. We model this evolution in diagrams called phylogenetic trees. Much like a family tree, a phylogenetic tree depicts the relationships between generations of cells. These models depend on taking physical samples from cancer patients.

But, these samples are imperfect, often leading to multiple phylogenetic trees. EssentCell finds the lineages that stay consistent among all optimal phylogenetic trees. These essential relations provide insight into the most probable ancestral relationships, painting a more precise picture of a tumor's evolutionary history.

## Content
  1. [Getting started](#start)
     * [Dependencies](#dep)
  2. [Usage instructions](#usage)
     * [I/O formats](#io)
     * [Examples](#example)
  3. [EssentCell Chernoff Bounds Calculation](#chernoff)
  4. [References](#references)


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

We provide two example input files for testing:
- [smallest.sorted.csv](data/smallest.sorted.csv)
- [smallest.sorted_with_missing.csv](data/smallest.sorted_with_missing.csv)

These files are not real biological datasets. They are included only to help users test and explore the program.

There are four real biological datasets that we have included in this repository, which were used for the experiments. 

- [Patient2.csv](data/Patient2.csv) and [Patient6.csv](data/Patient2.csv) contains the acute lymphoblastic leukemia (ALL) data of Patient 2 and Patient 6 in [1].
- [ER+](data/snvdata_Clonal_evo_47n_41m.csv) is an oestrogen-receptor-positive (ER+) breast cancer dataset from [2].
- [CRC1](data/CRC1.csv) is the human colorectal cancer (CRC) dataset (Patient 1) from [3].

Below is an example of how to run the script using smallest.sorted.csv as the input file.
To execute the program:
* Navigate to the directory where EssentCell.py is located.
* Ensure that the input file is in the same directory.
* Run the script using a command such as:

```
python EssentCell.py smallest.sorted.csv 0 2 --verbose -timeout 300
```

This command runs the EssentCell program using the ```smallest.sorted.csv``` file located in the root directory (copy it from the ```data``` directory to the root directory). The script will execute the analysis for ```k = 0, k = 1, and k = 2```. For groups larger than size 2, ILP calls are limited to a 300-second timeout.

Since no custom output directory was specified, the program will create a default folder named ```results_default``` in the working directory. Inside this folder, a subdirectory named ```smallest.sorted``` will be generated to store all output files associated with the ```smallest.sorted.csv``` input.

Then the program will create the following files.
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

For each value of ```k```, the file ```smallest.sorted.k.esspairs.txt``` contains the essential relation graph represented as an edge list (prior to collapsing strongly connected components).
The file ```smallest.sorted.k.esspairs.verbose.txt``` includes additional details such as the number of ILP calls made, total runtime for that specific ```k``` value, and other graph-related metadata.
For a clearer understanding of the output format and contents, please refer to the example output and the actual generated result files.

The file ```smallest.sorted_kappa_2.graph_info.txt``` contains the intersection graph, which includes edges present in all essential relation graphs from ```k = min to k = kmax```. Before saving, this graph is further processed by collapsing strongly connected components and performing a transitive reduction.

The file ```smallest.sorted_kappa_2.graph_persist.txt``` contains the same graph, but with mutation labels included.

Important: If you run the command multiple times:

```smallest.sorted.k.esspairs.txt``` files are recreated from scratch on each run.

Data in ```smallest.sorted.k.esspairs.verbose.txt``` files is appended, preserving records from previous runs.

If you prefer to start fresh without appending, you can delete the verbose files before running the command again.

#### Another example:

```
python EssentCell.py Patient2.csv  0  2 -result_folder Results/Results_"$current_month"_"$current_day" --verbose -timeout 300
```

Users can dynamically structure the output folder by passing variables when running the script from a shell script.


#### Example Using an Input with Missing Data

```
python EssentCell.py smallest.sorted_with_missing.csv  0  2 --verbose
```
We provide an example input file, ```smallest.sorted_with_missing.csv```, which contains randomly generated missing data entries. The program will handle these missing values, and the corresponding output files can be found in the output folder associated with this input.

#### Example using an Input and computing essential relations without group testing

If you want to compute essential relations without using group testing, you can pass a special flag (```-disable_gt```) to the script. This instructs the program to use the naive approach instead.
In this example, the output is stored in a new folder named ```result_wo_gt```.

```
python smallest.sorted.csv 0 2 -result_folder result_wo_gt --verbose -disable_gt
```

Refer the outout folder ```result_wo_gt``` for the output files: [result_wo_gt](result_wo_gt).

#### Example: Generating the Final Essential Relation Graph
You can use ```generategraph.py``` to create the final essential relation graph. To use this script:
1. First, run the main EssentCell script to generate the necessary output files.
2. Create a file named mutations.txt containing the mutation labels. Each row should have a column ID and mutation label, separated by a space. Start with column id 0.

Example format:
```
0 TP53
1 KRAS
2 EGFR
```
3. Place ```mutations.txt``` in the output folder corresponding to its input CSV file.

The program will read this file to label edges in the graph with the appropriate mutation names. Check the output folders for an example mutations.txt file. Example ```mutations.txt``` can be found [here](results_default/Patient2/mutations.txt).

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
                        This arguement should be pointing to the result folder where your output files are residing. Inside this folder there should be another folder named using the input file name excluding the .csv extension. By default this folder is results_default.
  -min_node_size MIN_NODE_SIZE
                        Minimum node size to be included in the graph. This argument looks at nodes that corresponds to clusters of size less than this threshold and remove it from the graph. This arguement is added to improve the clarity of the final diagram. If the user wishes to keep all nodes, then ignore this argument.
  -cluster_prefix CLUSTER_PREFIX
                        Cluster prefix to be used when generating the graph; This prefix will be used to name the names of the nodes, and cluster id will be added after cluster prefix.
  -do_not_keep_connections_when_deleting
                        This would ensure that we keep the connections of edges when we delete nodes that are too small; If the user do not want to keep the connections that goes through nodes that are removed using earlier arguement, then make sure this flag is enabled. By default, the program keeps the connections that goes through deleted nodes. This option is here to improve the clarity of the final figure.
  -max_number_of_mutation_labels MAX_NUMBER_OF_MUTATION_LABELS
                        Maximum number of mutation labels in the final graph. Default value is to include all mutation labels. Any edge label that contains more than max number of edge labels will not be displayed in the final graph. This option is there to improve the clarity of the figure.

```

#### Example: generating final graph

```
python generategraph.py Patient2.csv 2 --verbose -node_fill_color None -min_node_size 2 -max_number_of_mutation_labels 1
```

Note that **[mutation.txt](results_default/Patient2/mutations.txt)** file for the Patient2.csv is located in ```results_default/Patient2``` folder. 

Running generategraph.py will produce two output files:
- ```Patient2_kappa_2.graph_cluster_id.txt```
    * Contains cluster IDs and the cells assigned to each cluster.
    * [Refer the following file](https://github.com/msu-alglab/essentcell/blob/main/results_default/Patient2/Patient2_kappa_2.graph_cluster_id.txt)
- Patient2_kappa_2.pdf
    * Contains the final essential relation graph.
    * [Refer the following file](results_default/Patient2/Patient2_kappa_2.pdf)

If the user wishes to **customize the final graph**, you can modify ```generategraph.py``` using Graphviz attributes to adjust the output according to your preferences.

#### Another example:

```
python generategraph.py smallest.sorted.csv 2 --verbose -result_folder result_wo_gt -node_fill_color None -min_node_size 0 -max_number_of_mutation_labels 6
```
In this example, the folder ```result_wo_gt``` is used as both the input and output directory for ```generategraph.py```.
Note that ```mutations.txt``` must exist in the corresponding subfolder inside ```result_wo_gt``` for the script to label the graph nodes correctly.

<a name="chernoff"></a>
## EssentCell Chernoff Bound Calculation 
To determine the appropriate ```kappa``` value for each dataset, we used a Chernoff bound to estimate the probability that the number of false positives ```#FP``` exceeded ```k```.

Please refer the following Google doc for Chernoff bounds calculation. We have included details on how kappa values were calculated for the datasets that were used in the experiments.

[Google Doc](https://docs.google.com/spreadsheets/d/17phrXbAOQIeapo4AA4rRBmVcv8G2S6Z_5giAjfuCDZo/edit?usp=sharing)

<a name="references"></a>
## References
\[1\] Gawad, C., Koh, W. and Quake, S.R., 2014. Dissecting the clonal origins of childhood acute lymphoblastic leukemia by single-cell genomics. Proceedings of the National Academy of Sciences, 111(50), pp.17947-17952.

\[2\] Wang, Y., Waters, J., Leung, M.L., Unruh, A., Roh, W., Shi, X., Chen, K., Scheet, P., Vattathil, S., Liang, H. and Multani, A., 2014. Clonal evolution in breast cancer revealed by single nucleus genome sequencing. Nature, 512(7513), pp.155-160.
Proceedings of the National Academy of Science 111, 50 (Dec. 2014),
17947–17952.

\[3\] Leung, M.L., Davis, A., Gao, R., Casasent, A., Wang, Y., Sei, E., Vilar, E., Maru, D., Kopetz, S. and Navin, N.E., 2017. Single-cell DNA sequencing reveals a late-dissemination model in metastatic colorectal cancer. Genome research, 27(8), pp.1287-1299.
