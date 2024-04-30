import requests
import time

def get_cds_annotation(transcript_id, max_attempts=5):
    """
    Fetch the genomic positions of the CDS for a given transcript from Ensembl REST API,
    ensuring they are associated with the specified transcript via the 'Parent' field.
    """
    url = f"https://rest.ensembl.org/overlap/id/{transcript_id}?feature=cds"
    headers = {"Content-Type": "application/json"}
    
    for attempt in range(max_attempts):
        response = requests.get(url, headers=headers)
        
        if response.status_code == 200:
            try:
                data = response.json()
                cds_annotations = []
                for item in data:
                    if item.get('Parent') == transcript_id:
                        cds_annotation = {
                            'start': item['start'],
                            'end': item['end'],
                            'strand': item['strand'],
                            'seq_region_name': item['seq_region_name'],
                            'parent_transcript_id': item['Parent'],
                        }
                        cds_annotations.append(cds_annotation)
                return cds_annotations
            except ValueError:
                print("Failed to parse JSON data. Retrying...")
        else:
            print(f"Error: Unable to fetch data from Ensembl REST API. Status code: {response.status_code}. Retrying...")
        time.sleep(5)  # Wait for 5 seconds before retrying
    
    print(f"Error: Unable to fetch data from Ensembl REST API. Status code: {response.status_code}")
    return []


def fetch_cds_for_transcript(transcript_id, max_attempts=5):
    """
    Fetch the CDS sequence for a given transcript.
    """
    url = f"http://rest.ensembl.org/sequence/id/{transcript_id}?type=cds"
    headers = {"Content-Type": "text/plain"}
    
    for attempt in range(max_attempts):
        response = requests.get(url, headers=headers)
        
        if response.ok:
            if response.text.strip():  # Check if response body is not empty
                return response.text
            else:
                print(f"Empty response body. Retrying...")
        else:
            print(f"Error fetching CDS for transcript {transcript_id}. Status code: {response.status_code}. Retrying...")
        
        time.sleep(5)  # Wait for 5 seconds before retrying
    
    print(f"Failed to retrieve CDS for transcript {transcript_id} after {max_attempts} attempts.")
    return None

def localize_snp_in_cds(transcript_id, snp_position):
    """
    Determine the position of a SNP within the CDS, considering the strand direction.
    """
    cds_annotations = get_cds_annotation(transcript_id)
    cds_sequence = fetch_cds_for_transcript(transcript_id)
    
    if not cds_annotations or cds_sequence is None:
        print("Error: Unable to fetch CDS annotations or sequence.")
        return None, None, None, None
    
    # Sort annotations by position; reverse if on the negative strand
    cds_annotations.sort(key=lambda x: x['start'], reverse=cds_annotations[0]['strand'] == -1)
    
    cds_position = 0
    found = False
    
    for annotation in cds_annotations:
        fragment_length = annotation['end'] - annotation['start'] + 1
        if annotation['start'] <= snp_position <= annotation['end']:
            if annotation['strand'] == 1:
                offset = snp_position - annotation['start']
            else:
                offset = annotation['end'] - snp_position
            cds_position += offset
            found = True
            break
        else:
            cds_position += fragment_length
    
    if not found:
        print("SNP does not fall within the provided CDS regions.")
        return None, None, None, None
    
    # Adjust for 0-based indexing in Python strings
    snp_cds_nucleotide = cds_sequence[cds_position] if cds_position < len(cds_sequence) else None
    
    return cds_position + 1, snp_cds_nucleotide, cds_sequence, annotation['strand']

def translate_nucleotide_to_protein(nucleotide_sequence):
    # Define the codon table
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    
    protein_sequence = ""
    
    # Ensure the sequence length is a multiple of 3


    if len(nucleotide_sequence) % 3 == 0:
        # Convert the nucleotide sequence to upper case to match the codon table keys
        nucleotide_sequence = nucleotide_sequence.upper()
        
        # Iterate over the sequence in steps of 3 to read codons
        for i in range(0, len(nucleotide_sequence), 3):
            codon = nucleotide_sequence[i:i+3]
            amino_acid = codon_table.get(codon, '?')  # Use '?' for unknown codon
            protein_sequence += amino_acid
    else:
        return "Sequence length is not a multiple of 3."

    return protein_sequence


def mutate_sequence(nucleotide_sequence, position, snp):
    """
    Mutates the nucleotide sequence at the specified position with the provided SNP.
    
    Parameters:
    - nucleotide_sequence: The original nucleotide sequence.
    - position: The position in the sequence to mutate (1-based indexing).
    - snp: The single nucleotide polymorphism (a single character) to insert at the position.
    
    Returns:
    - The mutated nucleotide sequence.
    """
    # Convert to 0-based indexing
    position -= 1
    
    if position < 0 or position >= len(nucleotide_sequence):
        return "Error: Position out of range."
    
    # Ensure the SNP is a single nucleotide
    if len(snp) != 1 or snp.upper() not in "ACGT":
        return "Error: Invalid SNP."
    
    # Perform the mutation
    mutated_sequence = nucleotide_sequence[:position] + snp + nucleotide_sequence[position+1:]
    
    return mutated_sequence

def complem(nucleotide):
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'


def find_amino_acid_change(seq1, seq2):
    """
    Finds a single amino acid change between two sequences and returns the mutation.
    
    Parameters:
    - seq1: The original amino acid sequence.
    - seq2: The mutated amino acid sequence.
    
    Returns:
    - A string describing the mutation (e.g., "G34S"), "Synonymous" if no changes,
      "Nonsense" if seq2 is shorter, or "Complex mutation" if the sequences have different lengths.
    """
    if len(seq1) == len(seq2):
        mutations = []
        for i in range(len(seq1)):
            if seq1[i] != seq2[i]:
                mutations.append(f"{seq1[i]}{i+1}{seq2[i]}")
        
        if not mutations:
            return "synonymous"
        elif len(mutations) == 1:
            return mutations[0]
        else:
            return "multiple mutations"
    elif len(seq2) < len(seq1):
        return "nonsense"
    else:
        return "complex mutation"




def display_best_match_and_snp(transcript_id, snp, snp_position):

    snp_cds_position, snp_cds_nucleotide, seq, strand = localize_snp_in_cds(transcript_id, snp_position)

    if snp_cds_position:
        # print(f"SNP position in CDS: {snp_cds_position}, Nucleotide: {snp_cds_nucleotide}")
        # print(seq)
        # print(seq)
        prot_cds = translate_nucleotide_to_protein(seq)
        if strand == -1:
            complement_nucl = complem(snp)
            mut_seq = mutate_sequence(seq, snp_cds_position, complement_nucl)
            prot_mut = translate_nucleotide_to_protein(mut_seq)
            if prot_mut != 'Sequence length is not a multiple of 3.':
                print('Sequence length is not a multiple of 3. Something wrong with the CDS sequence in data base or incomplete sequence')
                pos_aa = int(snp_cds_position/3) 
                # print(snp_cds_position)
                # print(pos_aa)
                # print(prot_cds)
                # print(prot_mut)
                aa_prot_cds = prot_cds[pos_aa]
                aa_prot_mut = prot_mut[pos_aa]
                if len(prot_cds) == len(prot_cds):
                    if aa_prot_cds == aa_prot_mut:
                        effect = "synonymous"
                    else:
                        effect = aa_prot_cds + str(pos_aa + 1) + aa_prot_mut
                else:
                    effect = "nonsense"
                return prot_cds, prot_mut, effect, complem(snp_cds_nucleotide), seq
            else:
                return None, None, None, None, None
        else:
            complement_nucl = snp
            mut_seq = mutate_sequence(seq, snp_cds_position, complement_nucl)
            prot_mut = translate_nucleotide_to_protein(mut_seq)
            if prot_mut != 'Sequence length is not a multiple of 3.':
                print('Sequence length is not a multiple of 3. Something wrong with the CDS sequence in data base or incomplete sequence')
                pos_aa = int(snp_cds_position/3) 
                # print(snp_cds_position)
                # print(pos_aa)
                # print(prot_cds)
                # print(prot_mut)
                aa_prot_cds = prot_cds[pos_aa]
                aa_prot_mut = prot_mut[pos_aa]
                if len(prot_cds) == len(prot_cds):
                    if aa_prot_cds == aa_prot_mut:
                        effect = "synonymous"
                    else:
                        effect = aa_prot_cds + str(pos_aa + 1) + aa_prot_mut
                else:
                    effect = "nonsense"
                return prot_cds, prot_mut, effect, snp_cds_nucleotide, seq
            else:
                return None, None, None, None, None

    else:
        print("SNP position could not be determined within the CDS.")
        return None, None, None, None, None

# # # Example usage
# transcript_id = "ENST00000616016"
# snp = 'T'
# snp_position = 943937  # Example SNP position

# print(display_best_match_and_snp(transcript_id, snp, snp_position))
