import requests
import time
from Bio.Seq import Seq
from fetch_exon import download_exon_sequence
from protein_effect import display_best_match_and_snp
import time
import json

def fetch_data_with_retry(url, headers, max_attempts=4):
    for attempt in range(max_attempts):
        response = requests.get(url, headers=headers)
        if response.ok:
            try:
                data = response.json()
                return data
            except json.JSONDecodeError:
                print("Failed to parse JSON data. Retrying...")
        print(f"Failed to retrieve data. Status code: {response.status_code}. Retrying...")
        time.sleep(5)  # Wait for 5 seconds before retrying
    return False

def fetch_protein_coding_transcripts(gene_name):
    url = f"http://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_name}?expand=1"
    headers = {"Content-Type": "application/json"}
    response = fetch_data_with_retry(url, headers)
    if response == False:
        print("Error fetching protein-coding transcripts")
        return "Error fetching protein-coding transcripts"
    else:
        data = response
        return [transcript['id'] for transcript in data['Transcript'] if transcript['biotype'] in ['protein_coding', 'protein_coding_LoF']]


def fetch_uniprot_id(transcript_id):
    url = f"https://rest.uniprot.org/uniprotkb/search?query={transcript_id}&fields=accession&format=json"
    response = requests.get(url)
    if response.ok:
        data = response.json()
        results = data.get('results', [])
        if results:
            return results[0].get('primaryAccession', "UniProt ID not found")
    else:
        print(f"Error fetching UniProt ID for transcript {transcript_id}")
        return "Error fetching UniProt ID"

def fetch_transcript_genomic_location_and_exons(transcript_id):
    url = f"http://rest.ensembl.org/lookup/id/{transcript_id}?expand=1"
    headers = {"Content-Type": "application/json"}
    response = fetch_data_with_retry(url, headers)
    if response == False:
        print(f"Error fetching genomic location and exons for transcript {transcript_id}")
        return {}, []
    else:
        data = response
        return {"start": data['start'], "end": data['end']}, [{"id": exon['id'], "start": exon['start'], "end": exon['end']} for exon in data.get('Exon', [])]
        
        

def check_snp_in_exons(snp_position, exons):
    for exon in exons:
        if exon['start'] <= snp_position <= exon['end']:
            return True, exon['id']
    return False, None

def fetch_cds_for_transcript(transcript_id):
    url = f"http://rest.ensembl.org/sequence/id/{transcript_id}?type=cds"
    headers = {"Content-Type": "text/plain"}
    response = requests.get(url, headers=headers)
    if response.ok:
        return response.text
    else:
        time.sleep(10)
        response2 = requests.get(url, headers=headers)
        if response2.ok:
            return response2.text
        else:    
            print(f"Error fetching CDS for transcript {transcript_id}")
            return None





def get_start_position_by_id(exon_list, exon_id):
    """
    Retrieves the start position of a specific exon given its ID from a list of exon information dictionaries.

    Parameters:
    - exon_list: A list of dictionaries, where each dictionary contains the 'id', 'start', and 'end' of an exon.
    - exon_id: The ID of the exon whose start position is to be retrieved.

    Returns:
    The start position of the specified exon as an integer, or None if the exon ID is not found in the list.
    """
    for exon in exon_list:
        if exon['id'] == exon_id:
            return exon['start']
    return None


def pipeline(gene_name, snp_position,alter_alle):
    transcripts = fetch_protein_coding_transcripts(gene_name)
    if transcripts  == 'Error fetching protein-coding transcripts':
        return 'nf', None, None, None, None
    else:
        snp_found_in_exons = False  # Initialize the flag here

        for transcript_id in transcripts:
            uniprot_id = fetch_uniprot_id(transcript_id)
            location, exons = fetch_transcript_genomic_location_and_exons(transcript_id)

            print(f"\nScaning transcript ID: {transcript_id}, UniProt ID: {uniprot_id}")
            print(f"Genomic Location: Start - {location['start']}, End - {location['end']}")

            in_exon, _ = check_snp_in_exons(snp_position, exons)
            if in_exon:
                print(transcript_id, alter_alle, snp_position)
                effect_snp = display_best_match_and_snp(transcript_id, alter_alle, snp_position)
                if effect_snp[0] == None:
                    return 'nsnps', None, None, None, None, None
                else:
                    prot_original = effect_snp[0]
                    prot_mut = effect_snp[1]
                    missense = effect_snp[2]
                    nucleotide_cds = effect_snp[3]
                    cds_seq = effect_snp[4]
                    return prot_original, prot_mut, missense, uniprot_id, nucleotide_cds, cds_seq
            else:
                return 'nsnps', None, None, None, None, None       

# gene_name = 'ZFYVE28'
# snp_position = 2390547
# alter_alle = 'C'


# res = pipeline(gene_name, snp_position,alter_alle)
# print(res)
