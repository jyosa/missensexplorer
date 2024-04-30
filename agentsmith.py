import os
import argparse
import requests
from datetime import datetime
from dateutil.parser import parse as parse_date
from bs4 import BeautifulSoup
import PyPDF2
from openai import OpenAI
from Bio import Entrez
from dotenv import load_dotenv

def main(gene_name, mutation_type, user_email):
    load_dotenv()
    client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))

    papers = query_pubmed(gene_name, mutation_type, user_email)
    if isinstance(papers, str):
        print(papers)
    else:
        summaries = summarize_papers(papers, f"the effects of {mutation_type} mutation on protein", client)
        display_summaries(summaries)


def query_pubmed(gene, mutation, email):
    Entrez.email = email
    query = f"{gene}[Gene Name] AND {mutation}[All Fields]"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=10, sort="relevance")
    record = Entrez.read(handle)
    handle.close()
    
    ids = record['IdList']
    if not ids:
        return "There are no papers related with the gene and mutation in PubMed database."

    ids = ids[:3] if len(ids) > 3 else ids
    papers = []

    for pubmed_id in ids:
        handle = Entrez.efetch(db="pubmed", id=pubmed_id, retmode="xml")
        xml_data = handle.read()
        handle.close()
        
        soup = BeautifulSoup(xml_data, 'xml')
        article_title_tag = soup.find(['ArticleTitle', 'articletitle', 'title'])
        article_title = article_title_tag.text if article_title_tag else "Title not found"
        pub_date_tag = soup.find('PubDate')
        pub_date = parse_pubmed_date(pub_date_tag.text) if pub_date_tag else None
        abstract = soup.find('AbstractText')
        abstract_text = abstract.text if abstract else "No abstract available."
        free_full_text = soup.find('FreeFullText')
        
        if free_full_text:
            pdf_url = free_full_text.text
            pdf_response = requests.get(pdf_url)
            pdf_path = f"{pubmed_id}.pdf"
            with open(pdf_path, 'wb') as f:
                f.write(pdf_response.content)
            papers.append({'title': article_title, 'text': convert_pdf_to_text(pdf_path), 'date': pub_date})
        else:
            papers.append({'title': article_title, 'text': abstract_text, 'date': pub_date})

    return papers

def parse_pubmed_date(date_str):
    try:
        return parse_date(date_str)
    except ValueError:
        return None  # Handle date parsing issues

def convert_pdf_to_text(pdf_path):
    text = ""
    with open(pdf_path, "rb") as file:
        pdf_reader = PyPDF2.PdfReader(file)
        for page in pdf_reader.pages:
            text += page.extract_text() + "\n"
    return text

def generate_detailed_summary(text, focus_on, client):
    response = client.chat.completions.create(
        model="gpt-4-turbo-preview",
        messages=[
            {
                "role": "system",
                "content": f"You are a helpful assistant trained to generate detailed summaries of research papers focusing on {focus_on}. Please provide a comprehensive summary."
            },
            {
                "role": "user",
                "content": text
            }
        ],
        temperature=0.5,
        max_tokens=2000
    )
    return response.choices[0].message.content.strip()


def summarize_papers(papers, focus_on, client):
    summaries = {}
    for paper in papers:
        summary = generate_detailed_summary(paper['text'], focus_on, client)
        summaries[paper['title']] = summary
    return summaries


def display_summaries(summaries):
    if not summaries:
        print("No summaries to display.")
        return
    
    for title, summary in summaries.items():
        print("\n" + "="*80)
        print("Title: ", title)
        print("-" * 80)
        print("Summary:\n", summary)
        print("="*80 + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Query PubMed for papers on gene mutations and summarize them.",
        epilog="Example: python agentsmith.py TUBA1A V409A your_email@example.com"
    )
    parser.add_argument("gene_name", type=str, help="Gene name to search for in PubMed. Example: BRCA1")
    parser.add_argument("mutation_type", type=str, help="Mutation type to search for in PubMed. Example: V600E")
    parser.add_argument("user_email", type=str, help="Email address for PubMed access (required by NCBI). Example: example@example.com")

    args = parser.parse_args()

    main(args.gene_name, args.mutation_type, args.user_email)
