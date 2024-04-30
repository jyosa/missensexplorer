import requests

def download_exon_sequence(transcript_id, exon_id):
    """
    Download the exon sequence given a transcript ID and exon ID using the Ensembl REST API.

    Parameters:
    - transcript_id: The ID of the transcript.
    - exon_id: The ID of the exon.

    Returns:
    The exon sequence as a string, or an error message if the sequence could not be retrieved.
    """
    # Ensembl REST API URL for fetching exon information
    url = f"http://rest.ensembl.org/sequence/id/{exon_id}"

    # Request headers to specify that we want the response in JSON format
    headers = {"Content-Type": "application/json"}

    try:
        # Make the request to the Ensembl REST API
        response = requests.get(url, headers=headers)
        response.raise_for_status()  # Raise an error for bad responses

        # Parse the JSON response
        data = response.json()

        # Check if the sequence is in the response and return it
        if 'seq' in data:
            return data['seq']
        else:
            return "Sequence not found in response."
    except requests.RequestException as e:
        return f"An error occurred while fetching the sequence: {e}"

# # Example usage
# transcript_id = "ENST00000343938"
# exon_id = "ENSE00001714593"
# print(download_exon_sequence(transcript_id, exon_id))
