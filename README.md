# TreeTopology
* Version: 0.1

## INSTALLATION
```
python3 setup.py install  # Installation
python3 setup.py test     # Testing
```

## Requirements
* RAxML (https://github.com/stamatak/standard-RAxML)
* CONSEL (https://github.com/shimo-lab/consel)

## USAGE
#### On Linux
```
ALIGN=examples/input/TestInput1.nex
CONST=examples/input/TestConstraints.txt
CONSL=path/to/consel
OTPUT='output_name'

python3 scripts/treetopology_launcher_CLI.py -a $INPUT -c $METAD --consel $CONSL -o $OTPUT
```
