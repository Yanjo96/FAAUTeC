# FAAUTeC Fully Automated Approximately Unbiased Test Comparer
* Version: 0.2

## INSTALLATION
```
python3 setup.py install  # Installation
python3 setup.py test     # Testing
```

## Requirements
* RAxML (https://github.com/stamatak/standard-RAxML)
* CONSEL (https://github.com/shimo-lab/consel)
* IQTree (http://www.iqtree.org/)
* Python3
* Biopython
* Dendropy
* Optional: ete3
* Optional: xelatex
* Optional: IQTree2 (http://www.iqtree.org/)

Don't forget to give all programs the permission to be executed

## USAGE
For a correct run I would recommend to specify the complete paths.

### Parameters
#### Mandatory
- `--alignment | -a`  
 absolute path to infile; infile in PHYLIP or FASTA format; Example: /path_to_input/test.phy

- `--constraint| -c`  
 absolute path to constraint file; infile in NEWICK format; Example: /path_to_input/tree.tre

#### Optional
- `--model`  
  Model for RAxML

- `--ml_inference`  
  Choose which program should run the ML-Tree calculation 'RAxML' or 'IQTree'

- `--path_iqtree`  
  absolute path to the iqtree executable

- `--path_raxml`  
  absolute path to the RAxML executable

- `--au_inference`  
  Choose program for AU-test calculation 'CONSEL' or 'IQTree' or 'IQTree2', multiple selection possible by ';' as delimiter, e.g. 'CONSEL;IQTree'

- `--path_consel`  
  path to consel executables

- `--path_iqtree2`  
  absolute path to the iqtree2 executable

- `--latex_format`  
  Creates a more beautiful Table with xelatex

- `--thread_number`  
  Number of maximal used threads

- `--version`  
  print version number and exit

#### On Linux
```
ALIGN=examples/input/FASTA/
CONST=examples/input/Constraint/tree_hypothesis.txt
CONSL=path/to/consel

FAAUTeC -a $ALIGN -c $CONST --ml_inference "RAxML" --au_inference "CONSEL" --path_consel $CONSL
```

#### On Windows
The program do not work on Windows currently.

## TODO
* Check requirements
* (Windows support)
