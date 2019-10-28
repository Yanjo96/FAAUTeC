# TreeTopology
* Version: 0.1

## INSTALLATION
```
python3 setup.py install  # Installation
python3 setup.py test     # Testing
```

## Requirements
* RAxML (https://github.com/amkozlov/raxml-ng)
* CONSEL (https://github.com/shimo-lab/consel)

## USAGE
#### On Linux
```
ALIGN=examples/input/TestInput1.nex
CONST=examples/input/TestConstraints.txt
RAxML=path/to/raxml
OTPUT=examples/output/TestOutput.txt

python3 scripts/treetopology_launcher_CLI.py -a $INPUT -c $METAD -r $RAxML -o $OTPUT
```
