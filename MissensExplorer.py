import os
from multiprocessing import Pool, cpu_count
from get_annotation import analyze_snp_effect, classify_snp_by_rna_type, exon_check
from function_variant import get_protein_function, get_variant_info
import pandas as pd
from clinvar import get_clinvar_info
import time
import requests
import argparse

MAX_RETRIES = 3
RETRY_DELAY_SECONDS = 2
MAX_REQUESTS_PER_SECOND = 15



def process_vcf(vcf_file_path):
    vcf_dict_list = []
    print('\nProcessing the VCF file. Here, INDELs are not taken into account!\n')
    try:
        with open(vcf_file_path, 'r') as vcf_file:
            for line in vcf_file:
                line = line.strip()
                if line.startswith("#"):
                    continue  # Skip header lines
                fields = line.split("\t")
                chrom = fields[0]
                pos = int(fields[1])  # Convert position to integer
                ref_allele = fields[3]
                alt_alleles = fields[4].split(',')

                for alt_allele in alt_alleles:
                    # Check if ref_allele or alt_alleles contain more than one nucleotide
                    if len(ref_allele) > 1 or len(alt_allele) > 1:
                        continue  # Skip this line

                    format_ = fields[8]
                    sample_info = fields[9]
                    sample_fields = sample_info.split(":")
                    geno = sample_fields[0]  # Extract genotype

                    if geno == '0/1':
                        genotype = 'Heterozygous'
                    elif geno == '0/0':
                        genotype = 'Homozygous'
                    elif geno == '1/1':
                        genotype = 'Homozygous'
                    else:
                        genotype = 'Unknown'

                    vcf_dict = {
                        'chromosome': chrom,
                        'position': pos,
                        'ref_allele': ref_allele,
                        'alt_alleles': alt_allele,
                        'zygosity': genotype
                    }

                    vcf_dict_list.append(vcf_dict)
    except FileNotFoundError:
        print("Error: The file does not exist or the path is incorrect.")
    except Exception as e:
        print(f"An error occurred: {e}")

    return vcf_dict_list


def process_snp_vcf(info_snp, individual, email):  # Add individual as an argument
    # This function processes a single SNP
    # Extract SNP information
    chro = info_snp['chromosome']
    pos = info_snp['position']
    ref_alle = info_snp['ref_allele']
    alt_alle = info_snp['alt_alleles']
    zygosity = info_snp['zygosity']
    print(f'\nMapping SNP for chromosome: {chro}, Genome position: {pos}, Reference allele: {ref_alle}, SNP: {alt_alle[0]}, Zygosity: {zygosity}')
    effect = exon_check(chro, pos, ref_alle, alt_alle, email)
    # Process the SNP effect based on the returned value from `exon_check`
    if effect[0] == False or effect[0] == 'no gene found' or effect[0] == 'nsnps':
        print('Non-coding SNPs')
        return None
    elif effect[0] == 'nf':
        print(f'\nName of the gene is wrong in Ensembl database')
        return None
    else:
        gene_name = effect[0]
        prot_original = effect[1]
        prot_mut = effect[2]
        missense = effect[3]
        uniprot_code = effect[4]
        nucleotide_cds = effect[5]
        seq_cds = effect[6]
        if missense == 'synonymous':
            if ref_alle == nucleotide_cds:
                print(f'\nSynonymous mutations found for gene: {gene_name}, Uniprot: {uniprot_code}')
            else:
                print(f'checking: reference allele = {ref_alle} CDS nucleotide = {nucleotide_cds}')
                print('Houston we got a problem!!!!')
                print(f'{gene_name}   {pos}  {alt_alle}')
        elif missense == 'nonsense':
            if ref_alle == nucleotide_cds:
                print(f'\nNonsense mutations found for gene: {gene_name}, Uniprot: {uniprot_code}')
                prot_function = get_protein_function(uniprot_code)
                var_info = get_variant_info(uniprot_code, missense)
                genomicLocation = var_info[0]
                codon = var_info[1]
                consequenceType = var_info[2]
                sources = var_info[3]
                somaticStatus = var_info[4]
                info_mut = get_clinvar_info(gene_name, missense)
                values = tuple(info_mut.values())
                Mutation_Description = values[0]
                Accession = values[1]
                Clinical_Significance = values[2]
                Conditions_Associated = values[3]
                Molecular_Consequence = values[4]
                Protein_Change = values[5]
                Publications_Submissions = values[6]
                return {
                    'individual': individual,
                    'Chromosome': chro,
                    'Genomic position': pos,
                    'ref_alle': ref_alle,
                    'alt_alle': alt_alle,
                    'gene': gene_name,
                    'zygosity' : zygosity,
                    'Uniprot': uniprot_code,
                    'Effect': missense,
                    'seq_cds': seq_cds,
                    'prot_original': prot_original,
                    'prot_mut': prot_mut,
                    'prot_function': prot_function,
                    'genomicLocation': genomicLocation,
                    'codon': codon,
                    'consequenceType': consequenceType,
                    'sources': sources,
                    'somaticStatus': somaticStatus,
                    'Mutation_Description': Mutation_Description,
                    'Accession': Accession,
                    'Clinical_Significance': Clinical_Significance,
                    'Conditions_Associated': Conditions_Associated,
                    'Molecular_Consequence': Molecular_Consequence,
                    'Protein_Change': Protein_Change,
                    'Publications_Submissions': Publications_Submissions
                }
            else:
                print(f'checking: reference allele = {ref_alle} CDS nucleotide = {nucleotide_cds}')
                print('Houston we got a problem!!!!')
                print(f'{gene_name}   {pos}  {alt_alle}')
        else:
            if ref_alle == nucleotide_cds:
                print(f'checking: reference allele = {ref_alle} CDS nucleotide = {nucleotide_cds}')
                print(f'\nMissense mutations found for gene: {gene_name}, Uniprot: {uniprot_code}, Mutation: {missense}')
                prot_function = get_protein_function(uniprot_code)
                print(f'Protein funtion:  {prot_function}')
                var_info = get_variant_info(uniprot_code, missense)
                if var_info is None:
                    genomicLocation, codon, consequenceType, sources, somaticStatus = '', '', '', '', ''
                else:
                    genomicLocation = var_info[0]
                    codon = var_info[1]
                    consequenceType = var_info[2]
                    sources = var_info[3]
                    somaticStatus = var_info[4]
                    info_mut = get_clinvar_info(gene_name, missense)
                    values = tuple(info_mut.values())
                    Mutation_Description = values[0]
                    Accession = values[1]
                    Clinical_Significance = values[2]
                    Conditions_Associated = values[3]
                    Molecular_Consequence = values[4]
                    Protein_Change = values[5]
                    Publications_Submissions = values[6]

                    return {
                        'individual': individual,
                        'Chromosome': chro,
                        'Genomic position': pos,
                        'ref_alle': ref_alle,
                        'alt_alle': alt_alle,
                        'gene': gene_name,
                        'zygosity': zygosity,
                        'Uniprot': uniprot_code,
                        'Effect': missense,
                        'seq_cds': seq_cds,
                        'prot_original': prot_original,
                        'prot_mut': prot_mut,
                        'prot_function': prot_function,
                        'genomicLocation': genomicLocation,
                        'codon': codon,
                        'consequenceType': consequenceType,
                        'sources': sources,
                        'somaticStatus': somaticStatus,
                        'Mutation_Description': Mutation_Description,
                        'Accession': Accession,
                        'Clinical_Significance': Clinical_Significance,
                        'Conditions_Associated': Conditions_Associated,
                        'Molecular_Consequence': Molecular_Consequence,
                        'Protein_Change': Protein_Change,
                        'Publications_Submissions': Publications_Submissions
                    }

            else:
                print(f'checking: reference allele = {ref_alle} CDS nucleotide = {nucleotide_cds}')
                print('Houston we got a problem!!!!')
                print(f'{gene_name}   {pos}  {alt_alle}')


def process_snp(name_job, chro, pos, ref_alle, alt_alle, zygosity, email):  # Add individual as an argument
    # This function processes a single SNP
    # Extract SNP information
    chro = chro
    pos = pos
    ref_alle = ref_alle
    alt_alle = alt_alle
    zygosity = zygosity
    print(f'\nMapping SNP for chromosome: {chro}, Genome position: {pos}, Reference allele: {ref_alle}, SNP: {alt_alle[0]}, Zygosity: {zygosity}')
    effect = exon_check(chro, pos, ref_alle, alt_alle, email)
    # Process the SNP effect based on the returned value from `exon_check`
    if effect[0] == False or effect[0] == 'no gene found' or effect[0] == 'nsnps':
        print('Non-coding SNPs')
        return None
    elif effect[0] == 'nf':
        print(f'\nName of the gene is wrong in Ensembl database')
        return None
    else:
        gene_name = effect[0]
        prot_original = effect[1]
        prot_mut = effect[2]
        missense = effect[3]
        uniprot_code = effect[4]
        nucleotide_cds = effect[5]
        seq_cds = effect[6]
        if missense == 'synonymous':
            if ref_alle == nucleotide_cds:
                print(f'\nSynonymous mutations found for gene: {gene_name}, Uniprot: {uniprot_code}')
            else:
                print(f'checking: reference allele = {ref_alle} CDS nucleotide = {nucleotide_cds}')
                print('Houston we got a problem!!!!')
                print(f'{gene_name}   {pos}  {alt_alle}')
        elif missense == 'nonsense':
            if ref_alle == nucleotide_cds:
                print(f'\nNonsense mutations found for gene: {gene_name}, Uniprot: {uniprot_code}')
                prot_function = get_protein_function(uniprot_code)
                var_info = get_variant_info(uniprot_code, missense)
                genomicLocation = var_info[0]
                codon = var_info[1]
                consequenceType = var_info[2]
                sources = var_info[3]
                somaticStatus = var_info[4]
                info_mut = get_clinvar_info(gene_name, missense)
                values = tuple(info_mut.values())
                Mutation_Description = values[0]
                Accession = values[1]
                Clinical_Significance = values[2]
                Conditions_Associated = values[3]
                Molecular_Consequence = values[4]
                Protein_Change = values[5]
                Publications_Submissions = values[6]
                return {
                    'individual': name_job,
                    'Chromosome': chro,
                    'Genomic position': pos,
                    'ref_alle': ref_alle,
                    'alt_alle': alt_alle,
                    'gene': gene_name,
                    'zygosity' : zygosity,
                    'Uniprot': uniprot_code,
                    'Effect': missense,
                    'seq_cds': seq_cds,
                    'prot_original': prot_original,
                    'prot_mut': prot_mut,
                    'prot_function': prot_function,
                    'genomicLocation': genomicLocation,
                    'codon': codon,
                    'consequenceType': consequenceType,
                    'sources': sources,
                    'somaticStatus': somaticStatus,
                    'Mutation_Description': Mutation_Description,
                    'Accession': Accession,
                    'Clinical_Significance': Clinical_Significance,
                    'Conditions_Associated': Conditions_Associated,
                    'Molecular_Consequence': Molecular_Consequence,
                    'Protein_Change': Protein_Change,
                    'Publications_Submissions': Publications_Submissions
                }
            else:
                print(f'checking: reference allele = {ref_alle} CDS nucleotide = {nucleotide_cds}')
                print('Houston we got a problem!!!!')
                print(f'{gene_name}   {pos}  {alt_alle}')
        else:
            if ref_alle == nucleotide_cds:
                print(f'checking: reference allele = {ref_alle} CDS nucleotide = {nucleotide_cds}')
                print(f'\nMissense mutations found for gene: {gene_name}, Uniprot: {uniprot_code}, Mutation: {missense}')
                prot_function = get_protein_function(uniprot_code)
                print(f'Protein funtion:  {prot_function}')
                var_info = get_variant_info(uniprot_code, missense)
                if var_info is None:
                    genomicLocation, codon, consequenceType, sources, somaticStatus = '', '', '', '', ''
                else:
                    genomicLocation = var_info[0]
                    codon = var_info[1]
                    consequenceType = var_info[2]
                    sources = var_info[3]
                    somaticStatus = var_info[4]
                    info_mut = get_clinvar_info(gene_name, missense)
                    values = tuple(info_mut.values())
                    Mutation_Description = values[0]
                    Accession = values[1]
                    Clinical_Significance = values[2]
                    Conditions_Associated = values[3]
                    Molecular_Consequence = values[4]
                    Protein_Change = values[5]
                    Publications_Submissions = values[6]

                    return {
                        'individual': name_job,
                        'Chromosome': chro,
                        'Genomic position': pos,
                        'ref_alle': ref_alle,
                        'alt_alle': alt_alle,
                        'gene': gene_name,
                        'zygosity': zygosity,
                        'Uniprot': uniprot_code,
                        'Effect': missense,
                        'seq_cds': seq_cds,
                        'prot_original': prot_original,
                        'prot_mut': prot_mut,
                        'prot_function': prot_function,
                        'genomicLocation': genomicLocation,
                        'codon': codon,
                        'consequenceType': consequenceType,
                        'sources': sources,
                        'somaticStatus': somaticStatus,
                        'Mutation_Description': Mutation_Description,
                        'Accession': Accession,
                        'Clinical_Significance': Clinical_Significance,
                        'Conditions_Associated': Conditions_Associated,
                        'Molecular_Consequence': Molecular_Consequence,
                        'Protein_Change': Protein_Change,
                        'Publications_Submissions': Publications_Submissions
                    }

            else:
                print(f'checking: reference allele = {ref_alle} CDS nucleotide = {nucleotide_cds}')
                print('Houston we got a problem!!!!')
                print(f'{gene_name}   {pos}  {alt_alle}')

def print_dict_as_table(data):
    # Find maximum width for each field and value for alignment
    max_key_len = max(len(key) for key in data.keys()) + 2
    max_val_len = max(len(value) if isinstance(value, str) else len(str(value)) for value in data.values())

    # Setting a maximum width for value column for better readability
    max_value_display_width = 60

    # Print the header
    print(f"{'Field'.ljust(max_key_len)} | {'Value'.ljust(max_value_display_width)}")
    print('-' * (max_key_len + max_value_display_width + 3))

    # Print each key-value pair in the dictionary
    for key, value in data.items():
        value = str(value)
        if len(value) > max_value_display_width:
            # Wrap text if it exceeds the maximum display width
            wrapped_text = '\n'.join([value[i:i+max_value_display_width] for i in range(0, len(value), max_value_display_width)])
        else:
            wrapped_text = value
        print(f"{key.ljust(max_key_len)} | {wrapped_text}")

def main(vcf_file_path, individual, email, num_processors=None):
    if num_processors is None:
        num_processors = cpu_count()

    print('\nAnalyzing the VCF file.\n')
    snps = process_vcf(vcf_file_path)  # Assume process_vcf is adjusted to handle a single VCF file.
    print(f"The VCF file has {len(snps)} SNPs.")

    columns_name = ['individual', 'Chromosome', 'Genomic position', 'ref_alle', 'alt_alle', 'gene', 'zygosity', 'Uniprot',
                    'Effect', 'seq_cds', 'prot_original', 'prot_mut', 'prot_function', 'genomicLocation', 'codon',
                    'consequenceType', 'sources', 'somaticStatus', 'Mutation_Description', 'Accession',
                    'Clinical_Significance', 'Conditions_Associated', 'Molecular_Consequence', 'Protein_Change',
                    'Publications_Submissions']

    file_path_sg = 'output_snps.csv'
    failed_to_retrieve_file = 'failed_to_retrieve.csv'

    # Open the CSV file in append mode before processing any SNPs
    with open(file_path_sg, 'a') as f, open(failed_to_retrieve_file, 'a') as fail_file:
        df_snp = pd.DataFrame(columns=columns_name)
        df_snp.head(0).to_csv(f, mode='w', index=False)  # Write the header

        # Use multiprocessing to process SNPs in parallel
        with Pool(processes=num_processors) as pool:
            results = [pool.apply_async(process_snp_vcf, args=(info_snp, individual, email)) for info_snp in snps]
            for res, info_snp in zip(results, snps):  # Use zip to iterate over both results and snps
                retries = 0
                while retries < MAX_RETRIES:
                    try:
                        result = res.get(timeout=5)  # Get the result from ApplyResult object
                        if result is not None:
                            if isinstance(result, dict):  # Check if the result is a dictionary
                                # Append the result directly to the CSV file
                                new_row_df3 = pd.DataFrame([result])
                                new_row_df3.to_csv(f, mode='a', index=False, header=False)
                                f.flush()  # Flush the buffer to ensure immediate writing
                            else:
                                print("Invalid result format:", result)
                        break  # Break out of retry loop if successful
                    except requests.exceptions.HTTPError as e:
                        if e.response.status_code == 429:  # Rate limit exceeded
                            retry_after = int(e.response.headers.get('Retry-After', RETRY_DELAY_SECONDS))
                            print(f"Rate limit exceeded. Retrying after {retry_after} seconds...")
                            time.sleep(retry_after)
                        else:
                            retries += 1
                            if retries >= MAX_RETRIES:
                                print(f"Error processing SNP after {MAX_RETRIES} retries: {e}")
                                # Log failed SNP info
                                fail_info = [info_snp['chromosome'], info_snp['position'], info_snp['ref_allele'], info_snp['alt_alleles']]
                                fail_file.write(','.join(map(str, fail_info)) + '\n')
                                break
                            else:
                                print(f"Error processing SNP, retrying {retries}/{MAX_RETRIES} in {RETRY_DELAY_SECONDS} seconds.")
                                time.sleep(RETRY_DELAY_SECONDS)
                    except Exception as e:
                        retries += 1
                        if retries >= MAX_RETRIES:
                            print(f"Error processing SNP after {MAX_RETRIES} retries: {e}")
                            # Log failed SNP info
                            fail_info = [info_snp['chromosome'], info_snp['position'], info_snp['ref_allele'], info_snp['alt_alleles']]
                            fail_file.write(','.join(map(str, fail_info)) + '\n')
                            break
                        else:
                            print(f"Error processing SNP, retrying {retries}/{MAX_RETRIES} in {RETRY_DELAY_SECONDS} seconds: {e}")
                            time.sleep(RETRY_DELAY_SECONDS)

    # Optionally, close the CSV file after processing all SNPs
    # f.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="SNP Analysis Tool",
        formatter_class=argparse.RawTextHelpFormatter  # This helps in formatting the help text
    )
    parser.add_argument(
        '--mode',
        choices=['vcf', 'snp'],
        required=True,
        help="Mode of operation:\n"
             "  'vcf' - Processes a VCF file.\n"
             "  'snp' - Processes a single SNP.\n\n"
             "Usage examples:\n"
             "  Process a VCF file: --mode vcf --vcf_file path/to/file.vcf --individual JohnDoe --email user@example.com --cores 4\n"
             "  Process a single SNP: --mode snp --chro chr1 --pos 953279 --ref_alle T --alt_alle C --zygosity Homozygous --email user@example.com"
    )
    parser.add_argument('--vcf_file', help="Path to a VCF file for processing.")
    parser.add_argument('--individual', help="Identifier for the individual associated with the VCF file.")
    parser.add_argument('--email', required=True, help="Email address for notifications or results.")
    parser.add_argument('--cores', type=int, help="Number of cores to use for processing the VCF file.")
    parser.add_argument('--chro', help="Chromosome of the SNP (for SNP mode).")
    parser.add_argument('--pos', type=int, help="Position of the SNP on the chromosome (for SNP mode).")
    parser.add_argument('--ref_alle', help="Reference allele of the SNP (for SNP mode).")
    parser.add_argument('--alt_alle', help="Alternative allele of the SNP (for SNP mode).")
    parser.add_argument('--zygosity', choices=['Homozygous', 'Heterozygous', 'Unknown'], help="Zygosity of the SNP (for SNP mode).")

    args = parser.parse_args()

    if args.mode == 'vcf' and args.vcf_file and args.individual and args.email:
        main(args.vcf_file, args.individual, args.email, num_processors=args.cores)
    elif args.mode == 'snp' and all([args.chro, args.pos, args.ref_alle, args.alt_alle, args.zygosity]):
        print("SNP Processing Result:")
        result = process_snp("Single SNP Analysis", args.chro, args.pos, args.ref_alle, args.alt_alle, args.zygosity, args.email)
        if result == None:
            print('Captain, the sensors are picking up no signs of data in this sector!')
        else:
            print_dict_as_table(result)
    else:
        parser.print_help()
        print("\nMissing required arguments for the selected mode.")
