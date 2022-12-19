import os
import requests
import re
from bs4 import BeautifulSoup
import json
import numpy as np
import pandas as pd
import time

from locus_search.id_mapping_tools import *



def id_mapping_pipeline(query_list, from_db):
    
    job_id = submit_id_mapping(from_db, 'UniProtKB', query_list)
    
    time.sleep(3)
    
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_search(link)
    
    return results



def get_UniProtKB_accession(query_list, from_db):
    
    if (from_db == 'Ensembl_Genomes'): dir_name = 'ensembl'
    elif (from_db == 'GeneID'): dir_name = 'ncbi'
    else: raise ValueError(f'please set \'Ensembl_Genomes\' or \'GeneID\' to from_db')
    
    results = id_mapping_pipeline(query_list, from_db)
    
    if ('failedIds' in results.keys()):
        id_correspondence = {}
        for res in results['results']:
            accession_from = res['from']
            accession_to = res['to']['primaryAccession']
            id_correspondence[accession_from] = accession_to
            
            file_name = f'outputs/ID_mapping/{dir_name}/{accession_from}_id_mapping.json'
            json.dump(res, open(file_name, mode='w'), indent=4)
        
        UniProtKB_accession_list = [id_correspondence[query]
                                if (query in id_correspondence.keys())
                                else None
                                for query in query_list]
    
    else:
        UniProtKB_accession_list = []
        for res in results['results']:
            accession_from = res['from']
            accession_to = res['to']['primaryAccession']
            UniProtKB_accession_list.append(accession_to)
            
            file_name = f'outputs/ID_mapping/{dir_name}/{accession_from}_id_mapping.json'
            json.dump(res, open(file_name, mode='w'), indent=4)
    
    return UniProtKB_accession_list



def get_UniRef_ID(query):
    
    UniRef_URL = f'https://rest.uniprot.org/uniref/search?query=uniprot_id:{query}'
    UniRef_dict = json.loads(requests.get(UniRef_URL).text)
    
    UniRef_ID_dict = {}
    for uniref_id in [d['id'] for d in UniRef_dict['results']]:
        UniRef_ID_dict[float(re.search(r'\d+', uniref_id).group()) / 100] = uniref_id
    
    return UniRef_ID_dict



def UniRef_search(query, identity=0.5, update=False):
    
    if (identity != 0.5 and identity != 0.9 and identity != 1.0):
        raise ValueError('identitiy value is inappropriate. Set to 0.5, 0.9 or 1.0')
    
    UniRef_ID_dict = get_UniRef_ID(query)
    UniRef_ID = UniRef_ID_dict[identity]
    
    file_name = f'outputs/UniRef/UniRef{str(int(identity * 100))}/{UniRef_ID}_result.xml'
    UniRef_text = None
    if (os.path.isfile(file_name)):
        if (not update):
            with open(file_name) as f:
                UniRef_text = f.read()
    if (UniRef_text is None):
        UniRef_URL = f'https://rest.uniprot.org/uniref/{UniRef_ID}.xml'
        UniRef_text = requests.get(UniRef_URL).text
        print(UniRef_text, file=open(file_name, mode='w'))
    
    UniRef_soup = BeautifulSoup(UniRef_text, 'lxml-xml')
    
    representative = UniRef_soup.find('representativeMember')
    if (UniRef_soup.find('member') is not None):
        members = UniRef_soup.find_all('member')
    
    representative_accession = re.search(r'property type="UniProtKB accession" value=".+"', str(representative)).group()
    representative_UniParc_ID = re.search(r'property type="UniParc ID" value=".+"', str(representative)).group()
    representative_info = (re.search(r'value=".+?"', representative_accession).group()[7:-1],
                           re.search(r'value=".+?"', representative_UniParc_ID).group()[7:-1])
    
    member_accession_list = [re.search(r'property type="UniProtKB accession" value=".+"', str(member)).group()
                             if (re.search(r'property type="UniProtKB accession" value=".+"', str(member)) is not None)
                             else None
                             for member in members]
    member_UniParc_IDs = [re.search(r'property type="UniParc ID" value=".+"', str(member)).group()
                          if (re.search(r'property type="UniProtKB accession" value=".+"', str(member)) is not None)
                          else re.search(r'dbReference.+id=".+"', str(member)).group()
                          for member in members]
    members_info = [(re.search(r'value=".+?"', accession).group()[7:-1], re.search(r'value=".+?"', UniParc_ID).group()[7:-1])
                    if (accession is not None)
                    else (accession, re.search(r'id=".+?"', UniParc_ID).group()[7:-1])
                    for accession, UniParc_ID in zip(member_accession_list, member_UniParc_IDs)]
    
    cluster_name = UniRef_soup.find('name').text.split(': ')[1]
    
    return representative_info, members_info, cluster_name



def UniRef_pipeline(UniProtKB_accession_list, identity=0.5, update=False):
    
    cluster_name_list = []
    for UniProtKB_accession in UniProtKB_accession_list:
        if (UniProtKB_accession is not None):
            representative_info, members_info, cluster_name = UniRef_search(UniProtKB_accession, identity, update)
            cluster_name_list.append(cluster_name)
        else:
            cluster_name_list.append(None)
    
    return cluster_name_list