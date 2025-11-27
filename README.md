# EssentCell

Cancer is an evolutionary disease. As cells multiply, they can take on new mutations, and through the accumulation of those mutations, they can become cancerous. Tracking this evolution is helpful for understanding how the disease progresses and the most effective ways to treat it. We model this evolution in diagrams called phylogenetic trees. Much like a family tree, a phylogenetic tree depicts the relationships between generations of cells. These models depend on taking physical samples from cancer patients.

But, these samples are imperfect, often leading to multiple phylogenetic trees. EssentCell finds the lineages that stay consistent among all optimal phylogenetic trees. These essential relations provide insight into the most probable ancestral relationships, painting a more precise picture of a tumor's evolutionary history.

## Content
  1. [Getting started](#start)
     * [Dependencies](#dep)
  2. [Usage instructions](#usage)
     * [I/O formats](#io)
     * [Example](#example)


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
  kmin                  k min value for the input program; this is the minimum fixed **1 -> 0** bit flips that will be used by the program. 
  kmax                  k max value for the input program; this is the maximum fixed **1 -> 0** bit flips that will be used by the program. Note that **kmin <= kmax**.

optional arguments:
  -h, --help            show this help message and exit
  -result_folder RESULT_FOLDER
                        The result folder name; by default result folder is "results_default". Inside the result_folder, program will create another folder using the input file name **excluding ".csv"*.
                        Output files corresponding to the input file will reside in this folder. 
  --verbose             Increase output verbosity
  -print_trace_of_constraint
                        Increase the output verbosity of gurobi ILP creations and prints out the constraints created by the ILP calls
  -timeout TIMEOUT      Timeout value for ILP calls; this should be input in seconds. Note that the timeout value is only considered for ILP calls that ivolves groups larger than 2. When we run a ILP calls between two cells, the timeout value will be ignored as we need to ensure that the ILP call is either feasible or infeasible. By default timeout value is infinity.. 
  -disable_gt           diable group testing; if the user desires not to use group testing, then enable this flag.

```

The output options are:
There will be several output files created by the program.




TODOS:
- Change code to handle the missing data.
- Change constraint 8 to handle a single z_i (new method) 
- Change input such that first row is name of mutation and column being name of cells.
- (Related) Per cell list strongly component id as an addition to output.
- Provide adjaceny matrix that list the mutations for the strongly connect components.


Figures:
- Add component/cluster id to nodes and number of cells.
- Add mutations to edges.
