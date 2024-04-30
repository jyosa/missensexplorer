import requests

def get_protein_function(uniprot_id):
    """
    Fetches the function of a protein from UniProt using its UniProt ID.

    Parameters:
    uniprot_id (str): The UniProt ID of the protein.

    Returns:
    str: The function of the protein.
    """
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.txt"
    response = requests.get(url)
    
    if response.status_code == 200:
        content = response.text
        function_start = content.find("CC   -!- FUNCTION:")
        if function_start != -1:
            function_end = content.find("CC   -!-", function_start + 1)
            if function_end == -1:
                function_end = None
            function_text = content[function_start:function_end].strip()
            return function_text.replace('CC   -!- FUNCTION: ', '').replace('CC       ', '')
        else:
            return "Function not found."
    else:
        return "Failed to fetch data from UniProt."



def get_variant_info(uniprot_id, variant_description):
    """
    Fetches detailed information about a specific variant from UniProt, given its description,
    and returns selected information as a tuple.

    Parameters:
    uniprot_id (str): The UniProt ID of the protein.
    variant_description (str): The description of the variant in the format "OriginalAAPositionNewAA".
    
    Returns:
    tuple: Contains genomicLocation, codon_change, consequenceType, aggregated sources from predictions,
           and somaticStatus. Returns None if the variant is not found.
    """
    # Parse the variant description.
    original_aa, position, new_aa = variant_description[0], variant_description[1:-1], variant_description[-1]
    
    # Assuming the API URL and JSON structure based on the snippet you've shared
    url = f"https://www.ebi.ac.uk/proteins/api/variation/{uniprot_id}"
    response = requests.get(url, headers={"Accept": "application/json"})
    
    if response.status_code == 200:
        variants = response.json().get('features', [])
        for variant in variants:
            if variant.get('type') == 'VARIANT' and variant.get('wildType') == original_aa \
               and variant.get('mutatedType') == new_aa and str(variant.get('begin')) == position:
                # Aggregate sources from predictions
                sources = set()
                for prediction in variant.get('predictions', []):
                    sources.update(prediction.get('sources', []))
                
                variant_tuple = (
                    variant.get('genomicLocation'),
                    variant.get('codon'),
                    variant.get('consequenceType'),
                    list(sources),  # Convert the set to a list
                    "Somatic" if variant.get('somaticStatus') == 1 else "Germline",
                )
                return variant_tuple
        return None
    else:
        return "Failed to fetch data from UniProt."


# Test the function with the given parameters
# uniprot_id = "Q9C0B2"
# variant_description = "R218W"
# print(get_protein_function(uniprot_id))
# variant_info = get_variant_info(uniprot_id, variant_description)
# print(variant_info)


