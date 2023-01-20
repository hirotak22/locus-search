import os
import requests
import re
from bs4 import BeautifulSoup
import json
import numpy as np
import pandas as pd
import time



def check_external_links(query):
    
    UniProt_URL = f'https://rest.uniprot.org/uniprotkb/{query}.xml'
    UniProt_text = requests.get(UniProt_URL).text
    UniProt_soup = BeautifulSoup(UniProt_text, 'lxml-xml')
    
    # query_name = UniProt_soup.find('name', attrs={'type': 'primary'}).text
    
    tag_GeneID = UniProt_soup.find('dbReference', attrs={'type' : 'GeneID'})
    if (tag_GeneID is not None):
        in_NCBI = True
    else:
        in_NCBI = False
    
    tag_Ensembl = UniProt_soup.find('dbReference', attrs={'type' : ['Ensembl', 'EnsemblBacteria', 'EnsemblFungi', 'EnsemblPlants', 'EnsemblProtists']})
    if (tag_Ensembl is not None):
        in_Ensembl = True
    else:
        in_Ensembl = False
    
    return tag_GeneID, in_NCBI, tag_Ensembl, in_Ensembl



def get_query_GeneID(tag_GeneID):
    return re.sub(r'\D', '', str(tag_GeneID))       # query_GeneID



def get_query_locus_via_NCBI(query_GeneID):
    
    NCBI_URL = f'https://eutils.NCBI.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={query_GeneID}'
    NCBI_text = requests.get(NCBI_URL).text
    NCBI_soup = BeautifulSoup(NCBI_text, 'lxml-xml')
    
    tag_GenBank = NCBI_soup.find('GenomicInfo')
    query_locus = tag_GenBank.find('ChrAccVer').text
    query_location = (tag_GenBank.find('ChrStart').text,
                      tag_GenBank.find('ChrStop').text)
    
    return query_locus, query_location



def get_feature_table_via_NCBI(query_locus, update=False):
    
    file_name_ft = 'ft_' + query_locus.replace('.', '_') + '.tsv'
    if (os.path.isfile(f'outputs/NCBI/feature_table/{file_name_ft}')):
        if (not update):
            return file_name_ft
    
    NCBI_URL = f'https://eutils.NCBI.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={query_locus}&rettype=ft&retmode=text'
    NCBI_text = requests.get(NCBI_URL).text
    
    with open(f'outputs/NCBI/feature_table/{file_name_ft}', mode='w') as f:
        print(NCBI_text, end='', file=f)
    
    return file_name_ft



def process_feature_table_NCBI_into_gene_list(query_locus, file_name_ft, update=False):
    
    file_name_gl = 'gl_' + query_locus.replace('.', '_') + '.json'
    if (os.path.isfile(f'outputs/NCBI/gene_list/{file_name_gl}')):
        if (not update):
            with open(f'outputs/NCBI/gene_list/{file_name_gl}') as f:
                class_dict_list = json.load(f)
            return class_dict_list, file_name_gl
    
    df_ft = pd.read_table(f'outputs/NCBI/feature_table/{file_name_ft}',
                          skiprows=1,
                          names=['end', 'start', '_class', 'class_info', 'description'])
    
    region_list = []
    region = []
    for i in range(len(df_ft)):
        if (not pd.isna(df_ft.start[i])):
            region.append(((df_ft.start[i]), df_ft.end[i]))
        else:
            if (region != []):
                region_list.append(sorted(region))
                region = []
    
    class_dict_list = []
    i, cnt = 0, 0
    while (i < len(df_ft)):
        if (not pd.isna(df_ft._class[i])):
            if (df_ft._class[i] == 'gene'):
                if (i != 0):
                    class_dict_list.append(tmp_list)
                tmp_list = []
            class_key = df_ft._class[i]
            class_dict = {class_key : {}}
            class_dict[class_key]['region'] = region_list[cnt]
            i += 1
            cnt += 1
            while (pd.isna(df_ft._class[i])):
                if (not pd.isna(df_ft.start[i])):
                    i += 1
                    continue
                
                class_dict[class_key][df_ft.class_info[i]] = df_ft.description[i]
                i += 1
                if (i == len(df_ft)): break
            tmp_list.append(class_dict)
    class_dict_list.append(tmp_list)
    
    with open(f'outputs/NCBI/gene_list/{file_name_gl}', mode='w') as f:
        json.dump(class_dict_list, f, indent=4)
    
    return class_dict_list, file_name_gl



def process_NCBI_json_into_DataFrame(class_dict_list, query_locus, update=False):
    
    gene_table = []
    for cdl in class_dict_list:
        protein_coding = 0      # False
        for cd in cdl:
            if ('gene' in cd.keys()):
                # start, end, gene_name, GeneID
                content = [cd['gene']['region'][0][0],
                           cd['gene']['region'][0][1],
                           cd['gene']['gene'],
                           cd['gene']['db_xref'].split(':')[1],
                           cd['gene']['gene_desc']]
            if ('CDS' in cd.keys()):
                protein_coding = 1     # True
        
        content.append(protein_coding)
        gene_table.append(content)
    
    df_gt = pd.DataFrame(gene_table, columns=['start', 'end', 'gene_name', 'GeneID', 'description', 'protein_coding'])
    
    file_name_gt = 'gt_' + query_locus.replace('.', '_') + '.tsv'
    if (os.path.isfile(f'outputs/NCBI/gene_table/{file_name_gt}')):
        if (not update):
            return df_gt, file_name_gt
    
    df_gt.to_csv(f'outputs/NCBI/gene_table/{file_name_gt}', sep='\t', header=True, index=False)
    
    return df_gt, file_name_gt



def search_nearby_genes_via_NCBI(df_gt, query_GeneID, scope=5):

    df_gt_filtered = df_gt[df_gt['protein_coding'] == 1].reset_index(drop=True)
    hit_index = df_gt_filtered[(df_gt_filtered['GeneID'] == query_GeneID)].index.item()
    
    return df_gt_filtered.loc[max(hit_index-scope, 0) : min(hit_index+scope, len(df_gt_filtered))].drop(['protein_coding'], axis=1)



def get_query_Ensembl_ID(tag_Ensembl):
    
    tag_gene_ID = tag_Ensembl.find('property', attrs={'type' : 'gene ID'})
    
    return re.search(r'value=".+"', str(tag_gene_ID)).group()[7:-1]     # query_Ensembl_ID



def get_query_specy_via_Ensembl(query_Ensembl_ID):
    
    Ensembl_URL = f'https://rest.ensembl.org/lookup/id/{query_Ensembl_ID}?content-type=application/json'
    Ensembl_json = requests.get(Ensembl_URL).text
    dict_Ensembl = json.loads(Ensembl_json)
    
    query_specy = dict_Ensembl['species']
    query_seq_region_name = dict_Ensembl['seq_region_name']
    
    return query_specy, query_seq_region_name



def get_gene_list_via_Ensembl(query_specy, query_seq_region_name, update=False):
    
    file_name = f'gl_{query_specy}_{query_seq_region_name}.json'
    if (os.path.isfile(f'outputs/Ensembl/gene_list/{file_name}')):
        if (not update):
            return None
    
    Ensembl_URL = f'https://rest.ensembl.org/overlap/region/{query_specy}/{query_seq_region_name}:..:-1?content-type=application/json;biotype=protein_coding;feature=gene'

    Ensembl_json = requests.get(Ensembl_URL).text
    dict_Ensembl = json.loads(Ensembl_json)
    
    with open(f'outputs/Ensembl/gene_list/{file_name}', 'w') as f:
        json.dump(dict_Ensembl, f, indent=4)
    
    return None



def process_Ensembl_json_into_DataFrame(query_specy, query_seq_region_name, update=False):
    
    if (os.path.isfile(f'outputs/Ensembl/gene_table/gt_{query_specy}_{query_seq_region_name}.csv')):
        if (not update):
            df_gt = pd.read_csv(f'outputs/Ensembl/gene_table/gt_{query_specy}_{query_seq_region_name}.csv', header=0)
            return df_gt
    
    with open(f'outputs/Ensembl/gene_list/gl_{query_specy}_{query_seq_region_name}.json') as f:
        dict_Ensembl = json.load(f)
    
    table_id_region = []
    for de in dict_Ensembl:
        id = de['gene_id']
        start = de['start']
        end = de['end']
        strand = de['strand']
        description = de['description']
        table_id_region.append([id, start, end, strand, description])
    
    table_id_region = sorted(table_id_region, key=lambda x: x[1])
    
    df_gt = pd.DataFrame(table_id_region, columns=['gene_id', 'start', 'end', 'strand', 'description'])
    df_gt.to_csv(f'outputs/Ensembl/gene_table/gt_{query_specy}_{query_seq_region_name}.csv', header=True, index=False)
    
    return df_gt



def search_nearby_genes_via_Ensembl(df_gt, query_Ensembl_ID, scope=5):

    strand = df_gt[df_gt['gene_id'] == query_Ensembl_ID].strand.item()
    df_gt_filtered = df_gt[df_gt['strand'] == strand].reset_index()
    hit_index = df_gt_filtered[df_gt_filtered['gene_id'] == query_Ensembl_ID].index.item()
    
    return df_gt_filtered.loc[max(hit_index-scope, 0) : min(hit_index+scope, len(df_gt_filtered))]



def NCBI_pipeline(tag_GeneID, scope=5, update=False):
    
    query_GeneID = get_query_GeneID(tag_GeneID)
    query_locus, query_location = get_query_locus_via_NCBI(query_GeneID)
    file_name_ft = get_feature_table_via_NCBI(query_locus, update=update)
    NCBI_json, file_name_gl = process_feature_table_NCBI_into_gene_list(query_locus, file_name_ft, update=update)
    df_gene_table, file_name_gt = process_NCBI_json_into_DataFrame(NCBI_json, query_locus, update=update)
    
    return search_nearby_genes_via_NCBI(df_gene_table, query_GeneID, scope)



def Ensembl_pipeline(tag_Ensembl, scope=5, update=False):
    
    query_Ensembl_ID = get_query_Ensembl_ID(tag_Ensembl)
    query_specy, query_seq_region_name = get_query_specy_via_Ensembl(query_Ensembl_ID)
    get_gene_list_via_Ensembl(query_specy, query_seq_region_name, update=update)
    df_gene_table = process_Ensembl_json_into_DataFrame(query_specy, query_seq_region_name, update=update)
    
    return search_nearby_genes_via_Ensembl(df_gene_table, query_Ensembl_ID, scope)



if (__name__ == '__main__'):
    
    query = 'P12345'
    
    tag_GeneID, in_NCBI, tag_Ensembl, in_Ensembl = check_external_links(query)
    
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', None)

    if (in_NCBI):
        df_output_NCBI = NCBI_pipeline(tag_GeneID)
        print('---- NCBI --------')
        print(df_output_NCBI)
    if (in_Ensembl):
        df_output_Ensembl = Ensembl_pipeline(tag_Ensembl)
        print('---- Ensemble ----')
        print(df_output_Ensembl)