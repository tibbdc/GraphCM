# GraphCM
A new method for processing of currency metabolites in metabolic networks based on graph theory

## About

The pipeline was written and tested with Python 3.8. The core libraries essential for the pipeline including: Cobra, Pandas, networkx, and related packages. 

## Installation

1. create bow_tie environment using conda:

```shell
$ conda create -n graph_CM python=3.8
```

2. install related packages using pip:

```shell 
$ conda activate graph_CM
$ pip install cobra
$ pip install networkx
$ pip install ipykernel
$ python -m ipykernel install --user --name graph_CM --display-name "graph_CM"
```

## Steps to reproduce the analysis in the publication

Download all data and analysis code from github (directlt download or use git clone). 

 ```shell
$ cd /file path/project save path/
$ git clone https://github.com/tibbdc/GraphCM.git
```

 All results can be reproduced by executing the Jupyter Python notebooks:

+ graph_currency_metabolites.ipynb
  + the main script of removing currency metabolites based on graph theory approach.