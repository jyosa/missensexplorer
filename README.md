# MissensExplorer


This tool analyzes single nucleotide polymorphisms (SNPs) from Variant Call Format (VCF) files, focusing on identifying and reporting potential effects of these SNPs based on various genomic databases and annotations. The tool can handle individual SNP analysis or process whole VCF files.

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
pip install pandas requests
```
