# locus-search
## Installation
```
$ git clone https://github.com/yksaba/locus-search.git
$ pip install ./locus-search/
```
## List of included tools
### 1. Locus Search
Tools to search for the locus of UniProt query using NCBI and Ensemble and to retrieve coordinates of genes around the query from NCBI and Ensemble.
### 2. ID Mapping
Tools to map between the identifiers used in one database, to the identifiers of another, e.g., from UniProt to Ensembl, or to PomBase, etc.  
The source code is copied from the code example provided in UniProt (https://www.uniprot.org/help/id_mapping).
### 3. UniRef Search
Tools to search for UniRef (UniRef50, 90, 100) of UniProt query.
## File structure
```
.
├── README.md
├── main.py
├── setup.py
├── notebook
│   ├── locus_search.ipynb
│   └── UniRef_search.ipynb
├── outputs
│   ├── NCBI
│   │   ├── feature_table
│   │   ├── gene_list
│   │   └── gene_table
│   ├── Ensemnl
│   │   ├── gene_list
│   │   └── gene_table
│   ├── ID_mapping
│   │   ├── from_NCBI
│   │   └── from_Ensembl
│   └── UniRef
│       ├── UniRef50
│       ├── UniRef90
│       └── UniRef100
└── src/locus_search
    ├── __init__.py
    ├── id_mapping_tools.py
    ├── locus_search_tools.py
    └── UniRef_search_tools.py
```
The repository is divided into code and outputs.  
Code contains Python implimentations of the three tools mentioned above, and the pipeline to use them in one-liner on command line, in addition to Jupyter notebooks as examples of each tool's use.  
Outputs consist of the original data obtained by API in running each tool and the data processed in Python. Each directory is briefly described below.
- `outputs/NCBI/feature_table, gene_list, gene_table`  
    An original data obtained by the API is output in `/feature table`, a formatted version of it in json format in `/gene_list`, and a table summarizing the coordinate, name, GeneID, description, and whether it is protein-coding or not for each gene in `/gene_table`, respectively.

- `outputs/Ensembl/gene_list, gene_table`  
    An original data obtained by the API is output in `/gene_list`, and a table summarizing the ID, coordinates, strand, and description of each gene is output in `/gene_table`, respectively.

- `outputs/ID_mapping/from_NCBI, from_Ensembl`  
    The results of the job to convert the gene IDs in each external database into UniProt accessions are output here.

- `outputs/UniRef`  
    UniRef search results obtained using UniRef Search for queries are output here. The output location is divided by UniRef50, UniRef90, and UniRef100.

## Computing Environment
This was originally developed using Anaconda Python 3.8.12 and the following packages and versions:
```
numpy==1.20.3
pandas==1.3.4
beautifulsoup4==4.10.0
requests==2.26.0
```
## Usage
```
$ cd locus-search
$ python main.py -h     # help
$ python main.py (UniProt accession)
```
Please refer to the notebooks for details on each tool and function.

## ChangeLog
### [1.0.1] - 2023-01-23
#### Fixed
- Fixed a problem with ignoring strands when searching for genes around a query via NCBI.
#### Changed
- Changed a few of the locus-search result outputs via NCBI to be the same as those via Ensembl.