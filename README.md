# MissensExplorer


This tool analyzes single nucleotide polymorphisms (SNPs) from Variant Call Format (VCF) files, focusing on identifying and reporting potential effects of these SNPs, in human, based on various genomic databases and annotations. The tool can handle individual SNP analysis or process whole VCF files.

## Features

- Processing of VCF files to extract SNP information.
- Analysis of individual SNPs for potential effects on genes.
- Integration with external databases like ClinVar and UniProt for enriched annotation.
- Export results in CSV format with detailed information.
- Rate limiting and retries for handling API usage.
- Utilization of multiprocessing to enhance performance.

## Requirements

The script requires Python 3.x and the following Python libraries:
- `pandas`
- `requests`
- `argparse`
- `multiprocessing`


To install the required Python packages, run:
```bash
pip install pandas requests biopython
```
## Usage

This tool can be used in two different modes, allowing for the processing of either entire VCF files or individual SNPs. Below are the command line options and examples for each mode.

### Command Line Arguments

The program supports several command line arguments to specify its operation:

- `--mode`: Specifies the mode of operation. Choose `vcf` for VCF file processing or `snp` for individual SNP analysis.
- `--email`: Required for all operations. Email address to send notifications or results.

#### VCF Mode Arguments

- `--vcf_file`: Path to the VCF file to process.
- `--individual`: Identifier for the individual whose VCF file is being processed.
- `--cores`: (Optional) Number of processor cores to use for parallel processing.

#### SNP Mode Arguments

- `--chro`: Chromosome number where the SNP is located.
- `--pos`: Genomic position of the SNP.
- `--ref_alle`: Reference allele of the SNP.
- `--alt_alle`: Alternative allele of the SNP.
- `--zygosity`: Zygosity status of the SNP. Options are `Homozygous`, `Heterozygous`, or `Unknown`.

### Usage Examples

#### Process a VCF File

To process a VCF file and analyze the SNPs for a given individual:

```bash
python missenexplorer.py --mode vcf --vcf_file path/to/file.vcf --individual JohnDoe --email user@example.com --cores 4
```

For a single SNP:

```bash
python missenexplorer.py --mode snp --chro chr1 --pos 953279 --ref_alle T --alt_alle C --zygosity Homozygous --email user@example.com
```


## Output Configuration

### Files Generated

- **output_snps.csv**: This CSV file contains detailed results for each SNP processed in VCF mode. It includes information such as genomic location, gene impact, clinical significance, and other relevant details.
- **failed_to_retrieve.csv**: Records SNPs that could not be successfully processed after the maximum number of retries. This helps in tracking and re-processing these entries if necessary.

### Customizing Output

You can customize how the output data is saved by modifying the `file_path_sg` and `failed_to_retrieve_file` variables in the script:

```python
file_path_sg = 'output_snps.csv'  # Change this to your preferred output file path for successful SNP analyses
failed_to_retrieve_file = 'failed_to_retrieve.csv'  # Change this for logging failed SNP retrievals
```

[Watch our tutorial video for VCF file on YouTube!](https://www.youtube.com/watch?v=lduh43umNCc)


[Watch our tutorial video for SNP on YouTube!](https://www.youtube.com/watch?v=2y5LjD6Z92A)




### Viewing Output Data

The output CSV files can be opened with any text editor, spreadsheet software like Microsoft Excel or Google Sheets, or processed through additional scripts for further analysis.



## Agent Smith

This Python script automates the retrieval of research papers from the PubMed database based on specified gene mutations and provides detailed summaries of those papers. The summaries are generated using the OpenAI API, making this an experimental tool that also requires an OpenAI API account for operation.

## Features

- Retrieves relevant research papers from PubMed based on gene and mutation queries.
- Summarizes the papers using OpenAI's GPT models.
- Converts PDFs to text for summarization.
- Displays summarized results directly in the console.

## Requirements

- Python 3.x
- Access to the OpenAI API (requires an API key)
- Internet connection for PubMed access and API communication
- Python packages: `requests`, `bs4`, `PyPDF2`, `python-dotenv`, `biopython`

## Installation

Before running the script, ensure that you have all the required packages installed:

```bash
pip install requests beautifulsoup4 PyPDF2 python-dotenv biopython openai
```

Additionally, you need to set up an environment variable for your OpenAI API key. This can be done in a .env file in the same directory as the script:

```bash
OPENAI_API_KEY='your_openai_api_key_here'
```

### Running the Script

The script is executed from the command line and requires three arguments:

- gene_name: The name of the gene to be researched.
- mutation_type: The specific mutation of interest.
- user_email: The email address to use with PubMed queries (required by NCBI).

Command Line Syntax:

```bash
python agentsmith.py <gene_name> <mutation_type> <user_email>
```
Example:

```bash
python agentsmith.py TUBA1A V409A your_email@example.com
```

[Watch our tutorial video for SNP on YouTube!](https://www.youtube.com/watch?v=J19IQ_NZSGQ)



# Powered by

### Universidad Simón Bolívar Barranquilla 
### Escuela Nacional del Deporte

- Briyis Fontecha briyis.fontecha@unisimon.edu.co
- Jorge Leyva jorge.leyva@unisimon.edu.co
- Yecid Mina yecid.mina@endeporte.edu.co
- Juvenal Yosa juvenal.yosa@unisimon.edu.co
