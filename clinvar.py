import requests

def get_amino_acid_map():
    return {
        'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 
        'C': 'Cys', 'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 
        'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys', 
        'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser', 
        'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
    }

def convert_mutation_to_hgvs(short_mutation):
    amino_acid_map = get_amino_acid_map()
    original, position, new = short_mutation[0], short_mutation[1:-1], short_mutation[-1]
    return f"p.{amino_acid_map.get(original, 'X')}{position}{amino_acid_map.get(new, 'X')}"

def get_clinvar_info(gene_name, short_mutation):
    precise_mutation = convert_mutation_to_hgvs(short_mutation)
    search_query = f"{gene_name}[Gene] AND {precise_mutation}[All Fields]"

    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "clinvar",
        "term": search_query,
        "retmode": "json",
        "retmax": 1
    }
    search_response = requests.get(search_url, params=params)
    search_result = search_response.json()

    if search_result['esearchresult']['idlist']:
        clinvar_id = search_result['esearchresult']['idlist'][0]
        summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        summary_params = {
            "db": "clinvar",
            "id": clinvar_id,
            "retmode": "json"
        }
        summary_response = requests.get(summary_url, params=summary_params)
        summary_result = summary_response.json()
        data_json = parse_clinvar_response(summary_result)
        if data_json:
            if 'result' in data_json:
                return data_json
            else:
                return {
                "Mutation Description": None,
                "Accession": None,
                "Clinical Significance": None,
                "Conditions Associated": [],
                "Molecular Consequence": [],
                "Protein Change": [],
                    "Publications and Submissions": {"SCV": [], "RCV": []}
                } 
            
    else:
        # Return a dictionary with empty values if no information is found
        return {
            "Mutation Description": None,
            "Accession": None,
            "Clinical Significance": None,
            "Conditions Associated": [],
            "Molecular Consequence": [],
            "Protein Change": [],
            "Publications and Submissions": {"SCV": [], "RCV": []}
        }

def parse_clinvar_response(clinvar_json):
    uid = clinvar_json['result']['uids'][0]
    result = clinvar_json['result'][uid]
    
    clinvar_info = {
        "Mutation Description": result.get('title'),
        "Accession": result.get('accession'),
        "Clinical Significance": result.get('germline_classification', {}).get('description'),
        "Conditions Associated": [trait_set.get('trait_name') for trait_set in result.get('germline_classification', {}).get('trait_set', [])],
        "Molecular Consequence": result.get('molecular_consequence_list', []),
        "Protein Change": result.get('protein_change', '').split(', '),
        "Publications and Submissions": {
            "SCV": result.get('supporting_submissions', {}).get('scv', []),
            "RCV": result.get('supporting_submissions', {}).get('rcv', [])
        }
    }
    return clinvar_info

# # Example usage
# info = get_clinvar_info("BRCA1", "R1720P")
# values = tuple(info.values())
# Mutation_Description = values[0]
# Accession = values[1]
# Clinical_Significance = values[2]
# Conditions_Associated = values[3]
# Molecular_Consequence = values[4]
# Protein_Change= values[5]
# Publications_Submissions = values[6]
# print(Mutation_Description,Accession,Clinical_Significance,Conditions_Associated,Molecular_Consequence,Protein_Change,Publications_Submissions)
