import os
import requests
import re
from bs4 import BeautifulSoup
import json
import numpy as np
import pandas as pd



def read_FASTA(path_FASTA):
    
    with open(path_FASTA) as f_fasta:
        lines = f_fasta.read().splitlines()
        header = lines[0]
        whole_sequence = ''.join(lines[1:])
    
    return header, whole_sequence



def reverse_complement(whole_sequence):
    
    return whole_sequence.translate(str.maketrans('ATGC', 'TACG'))[::-1]
    
    

def get_nucleotide_sequence_via_NCBI(query_locus, update=False):
    
    file_name_FASTA = 'nucseq_' + query_locus.replace('.', '_') + '.fasta'
    if (os.path.isfile(f'outputs/NCBI/nucleotide_sequence/whole_sequence/{file_name_FASTA}')):
        if (not update):
            return None
    
    NCBI_URL = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={query_locus}&strand=true&rettype=fasta&retmode=text'
    NCBI_text = requests.get(NCBI_URL).text
    
    with open(f'outputs/NCBI/nucleotide_sequence/whole_sequence/{file_name_FASTA}', mode='w') as f:
        f.write(NCBI_text)
    
    # get and save reverse complement
    header, whole_sequence = read_FASTA(f'outputs/NCBI/nucleotide_sequence/whole_sequence/{file_name_FASTA}')
    whole_sequence = reverse_complement(whole_sequence)
    file_name_FASTA_reverse_complement = 'nucseq_' + query_locus.replace('.', '_') + '_reverse_complement.fasta'
    idx_insert = header.find(' ')
    with open(f'outputs/NCBI/nucleotide_sequence/whole_sequence/{file_name_FASTA_reverse_complement}', mode='w') as f:
        f.write(f'{header[:idx_insert]}:c{len(whole_sequence)}-1{header[idx_insert:]}\n')    # writing header
        for i in range(0, len(whole_sequence), 70):
            f.write(whole_sequence[i:min(i+70,len(whole_sequence))] + '\n')    # writing gene sequence

    return None



def extract_gene_sequence_NCBI(query_locus, df_output_NCBI, update=False):
    
    file_name_FASTA_forward = 'nucseq_' + query_locus.replace('.', '_') + '.fasta'
    file_name_FASTA_reverse = 'nucseq_' + query_locus.replace('.', '_') + '_reverse_complement.fasta'
    
    header, whole_sequence_forward = read_FASTA(f'outputs/NCBI/nucleotide_sequence/whole_sequence/{file_name_FASTA_forward}')
    _, whole_sequence_reverse = read_FASTA(f'outputs/NCBI/nucleotide_sequence/whole_sequence/{file_name_FASTA_reverse}')
    
    idx_insert = header.find(' ')
    
    for start, end, GeneID, strand in zip(df_output_NCBI['start'], df_output_NCBI['end'], df_output_NCBI['GeneID'], df_output_NCBI['strand']):
        if (os.path.isfile(f'outputs/NCBI/nucleotide_sequence/gene_sequence/nucseq_{GeneID}.fasta')):
            if (not update):
                continue
        
        if (start[0] == '>'): start_idx = int(start[1:])
        else: start_idx = int(start)
        
        if (end[0] == '<'): end_idx = int(end[1:])
        else: end_idx = int(end)
        
        if (strand == 1):
            with open(f'outputs/NCBI/nucleotide_sequence/gene_sequence/nucseq_{GeneID}.fasta', mode='w') as f_gene:
                f_gene.write(f'{header[:idx_insert]}:{start}-{end}{header[idx_insert:]}\n')    # writing header
                gene_seq = whole_sequence_forward[(start_idx-1):end_idx]
                for i in range(0, len(gene_seq), 70):
                    f_gene.write(gene_seq[i:min(i+70,len(gene_seq))] + '\n')    # writing gene sequence
        elif (strand == -1):
            with open(f'outputs/NCBI/nucleotide_sequence/gene_sequence/nucseq_{GeneID}.fasta', mode='w') as f_gene:
                f_gene.write(f'{header[:idx_insert]}:c{start}-{end}{header[idx_insert:]}\n')    # writing header
                gene_seq = whole_sequence_reverse[-start_idx:-(end_idx-1)]
                for i in range(0, len(gene_seq), 70):
                    f_gene.write(gene_seq[i:min(i+70,len(gene_seq))] + '\n')    # writing gene sequence
        else:
            raise ValueError('strand value must be 1 or -1')
    
    return None



def get_gene_sequence_via_Ensembl(df_output_Ensembl, update=False):
    
    for gene_id in df_output_Ensembl['gene_id']:
        
        file_name_FASTA = f'nucseq_{gene_id}.fasta'
        if os.path.isfile(f'outputs/Ensembl/nucleotide_sequence/gene_sequence/{file_name_FASTA}'):
            if (not update):
                continue
        
        Ensembl_URL = f'https://rest.ensembl.org/sequence/id/{gene_id}?content-type=text/x-fasta'
        Ensembl_text = requests.get(Ensembl_URL).text
        with open(f'outputs/Ensembl/nucleotide_sequence/gene_sequence/{file_name_FASTA}', mode='w') as f:
            f.write(Ensembl_text)
    
    return None