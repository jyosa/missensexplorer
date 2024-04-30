import requests
import os
from Bio import Entrez, SeqIO
from snp_class import pipeline
import time
import json

def fetch_data_with_retry(base_url, query, max_attempts=10):
    for attempt in range(max_attempts):
        response = requests.get(base_url + query)
        if response.status_code == 200:
            try:
                data = response.json()
                return data
            except json.JSONDecodeError:
                print("Failed to parse JSON data. Retrying...")
        #print(f"Failed to retrieve data. Status code: {response.status_code}. Retrying...")
        time.sleep(5)  # Wait for 5 seconds before retrying
    return {"error": f"Failed to retrieve data after {max_attempts} attempts"}


def fetch_ucsc_track_data(genome, chromosome, position):
    """
    Fetch gene data surrounding a specific position from the UCSC Genome Browser API.
    
    Parameters:
    - genome: str, the genome assembly (e.g., 'hg38').
    - chromosome: str, the chromosome (e.g., 'chr1').
    - position: int, the position of interest.
    
    Returns:
    - The API response as JSON or an error message.
    """
    
    track = 'knownGene'
    base_url = "https://api.genome.ucsc.edu/getData/track"
    # Adjust the start and end coordinates according to UCSC rules
    start = position  # Start coordinate is 0-based
    end = position + 1  # End coordinate is 1-based, so add 1 for single nucleotide positions
    query = f"?genome={genome};track={track};chrom={chromosome};start={start};end={end}"
    return fetch_data_with_retry(base_url, query)


def parse_ucsc_output(api_output):
    """
    Parse the JSON output from the UCSC Genome Browser API to extract gene details.
    
    Parameters:
    - api_output: dict, the JSON response from the UCSC API.
    
    Returns:
    - A list of dictionaries, each containing details of a gene.
    """
    gene_details = []
    if 'knownGene' in api_output:
        for gene in api_output['knownGene']:
            details = {
                'geneName': gene['geneName'],
                'geneName2': gene['geneName2'],
                'chromosome': gene['chrom'],
                'startPosition': gene['chromStart'],
                'endPosition': gene['chromEnd'],
                'geneType': gene['geneType'],
                'strand': gene['strand']
            }
            gene_details.append(details)
    return gene_details


def fetch_mrna_sequences(transcript_ids, email='.comjuvenal.yosa@gmail.com'):
    """
    Fetch mRNA sequences for a list of transcript IDs from NCBI.

    Parameters:
    - transcript_ids: List[str], a list of RefSeq transcript IDs (e.g., ['NM_001005484', 'NM_000546']).
    - email: str, your email address, used for NCBI queries.

    Returns:
    - dict, mapping from transcript IDs to their mRNA sequences.
    """
    Entrez.email = email  # Always set your email
    sequences = {}

    for transcript_id in transcript_ids:
        handle = Entrez.efetch(db="nucleotide", id=transcript_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        sequences[transcript_id] = str(record.seq)
    
    return sequences



def classify_snp_by_rna_type(genes_at_position, snp_position):
    """
    Classify a SNP as affecting a non-coding RNA or coding RNA based on overlapping genes.

    Parameters:
    - genes_at_position: List[Dict], containing gene information from UCSC API results.
    - snp_position: int, the genomic position of the SNP.

    Returns:
    - str, classification of the SNP as 'non-coding RNA' or 'coding RNA'.
    """
    # Default classification is 'non-coding RNA'
    classification = 'non-coding RNA'
    
    # Check each gene at the SNP position for their type
    for gene in genes_at_position:
        if gene['geneType'] == 'protein_coding':
            # If any of the overlapping genes are protein-coding, classify the SNP as affecting a coding RNA
            classification = 'coding RNA'
            transcript = gene['geneName']
            if len(transcript) > 13:
                transcript = gene['geneName2']
            return transcript
            break
        # Add any other conditions here if you have more specific types to check
        else:
    
            return classification



def analyze_snp_effect(chromosome, position, ref_allele, alt_alleles, genome='hg38'):
    """
    Analyze the effect of an SNP, automating the determination of its genomic context and potential impact.
    
    Parameters:
    - chromosome: str, the chromosome number (e.g., '1').
    - position: int, the SNP position on the chromosome.
    - ref_allele: str, the reference allele.
    - alt_alleles: str, the alternate alleles.
    - genome: str, the genome assembly version (default: 'hg38').
    
    Returns:
    - A dictionary with analysis results, including genomic context and predicted effects.
    """
    # Fetch UCSC track data for the SNP position
    ucsc_response = fetch_ucsc_track_data(genome, chromosome, position)

    
    # Parse the UCSC response to get gene details
    genes_at_position = parse_ucsc_output(ucsc_response)
    
    
    if not genes_at_position:
        return 'No gene found'
    
    # Placeholder: Analysis of SNP effect based on the gene's context
    #print(genes_at_position)
    genename = classify_snp_by_rna_type(genes_at_position, position)
    return genename



def exon_check(chromosome,position,ref_allele,alt_alleles,email):                 
    result = analyze_snp_effect(chromosome, position, ref_allele, alt_alleles)
    if result == 'non-coding RNA':
        return False, None, None, None, None, None, None
    elif result == 'No gene found':
        return 'no gene found', None, None, None, None, None, None
    else:
        gene_name = result.split('.')[0]
        #print(gene_name)
        res_pipeline = pipeline(gene_name, position, alt_alleles[0])
        if res_pipeline[0] == 'nsnps':
            return 'nsnps', None, None, None, None, None, None
        elif res_pipeline[0] == 'nf':
            return 'nf', None, None, None, None, None, None
        else:
            prot_original = res_pipeline[0]
            prot_mut = res_pipeline[1]
            missense = res_pipeline[2]
            uniprot = res_pipeline[3]
            nuclotide_cds = res_pipeline[4]
            seq_cds = res_pipeline[5]
        return gene_name, prot_original, prot_mut, missense, uniprot, nuclotide_cds, seq_cds


       



# # # # # Example usage
#chromosome = '4'
#position = 2661936 # Example position, adjust as needed
#ref_allele = 'C'
#alt_alleles = 'T'
#email = 'juenal.yosa@gmail.com'

#print(exon_check(chromosome,position,ref_allele,alt_alleles,email))


