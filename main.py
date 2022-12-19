import numpy as np
import pandas as pd
import argparse
from locus_search import *



def locus_search(query, scope=5, identity=0.5, update=False):
    
    tag_GeneID, in_NCBI, tag_Ensembl, in_Ensembl = check_external_links(query)

    df_output_NCBI, df_output_Ensembl = None, None
    
    if (in_NCBI):
        df_output_NCBI = NCBI_pipeline(tag_GeneID, scope=scope, update=update)
        query_list = df_output_NCBI['GeneID'].to_list()
        UniProtKB_accession_list = get_UniProtKB_accession(query_list, 'GeneID')
        cluster_name_list =  UniRef_pipeline(UniProtKB_accession_list, identity=identity, update=update)
        df_output_NCBI['UniProtKB accession'] = UniProtKB_accession_list
        df_output_NCBI['UniRef50 cluster'] = cluster_name_list
        
    if (in_Ensembl):
        df_output_Ensembl = Ensembl_pipeline(tag_Ensembl, scope=scope, update=update)
        query_list = df_output_Ensembl['gene_id'].to_list()
        UniProtKB_accession_list = get_UniProtKB_accession(query_list, 'Ensembl_Genomes')
        cluster_name_list =  UniRef_pipeline(UniProtKB_accession_list, identity=identity, update=update)
        df_output_Ensembl['UniProtKB accession'] = UniProtKB_accession_list
        df_output_Ensembl['UniRef50 cluster'] = cluster_name_list
    
    return df_output_NCBI, df_output_Ensembl



if (__name__ == '__main__'):
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('query', help='UniProt accession (ex. P12345, P74258)')
    parser.add_argument('--scope', type=int, default=5, help='Range of proteins around the query to be searched. (int, default = 5)')
    parser.add_argument('--identity', type=float, choices=[0.5, 0.9, 1.0], default=0.5, help='A sequence identity threshold in an UniRef cluster. (one of [0.5, 0.9, 1.0], default = 0.5)')
    parser.add_argument('--update', action='store_true', help='Whether to search again for previously searched items. (bool, default = False)')
    
    args = vars(parser.parse_args())
    
    df_output_NCBI, df_output_Ensembl = locus_search(*args.values())
    
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', None)
    
    if (df_output_NCBI is not None):
        print('---- NCBI --------')
        print(df_output_NCBI)
        print()
    
    if (df_output_Ensembl is not None):
        print('---- Ensemble ----')
        print(df_output_Ensembl)
        print()