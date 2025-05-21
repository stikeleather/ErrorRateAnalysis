# ErrorRateAnalysis

Create the environment with either conda or mamba:

```
conda create -n sera python=3.12 biopython numpy pandas openpyxl git
```
```
mamba create -n sera python=3.12 biopython numpy pandas openpyxl git
```
```
conda activate sera
```
Clone the repository
```
git clone https://github.com/stikeleather/ErrorRateAnalysis-public.git
cd ErrorRateAnalysis-public
```

1) Parse an annotation file into proteome and gene files. When selecting a --parse argument, make sure to use a feature that contains a unique ID. Locus tags usually work, but annotation files vary.

```
python parse_genbank.py --gb wt_ecoli.gbff --tt 11 --parse lt
```
```
Arguments:
--gb   :    GenBank file to be parsed (.gbff, .gbk, or.gb) [gbff preferred]
--tt   :    Translation table number <1-33> (e.g., 11)
--parse:    Parse on either locus tags with option "lt" or gene names with option "g", or proteins IDs with option "p"
```

Using Proteome Discoverer, or similar software, perform a search with the produced proteome file. Export proteins from Proteome Discoverer to Excel. (Skip for tutorial).


2) Generate Custom Database

This depends on the exported proteins tabs from Proteome Discoverer, or similar software in the .xlsx format. The requirement for --gen_xlsx is that it has a column header for either "Accession" or "Master Protein Accession".

```
python generate_custom_database.py --gen_xlsx wt_ecoli_generic.xlsx --proteome_file wt_ecoli_proteome.fasta --gene_file wt_ecoli_genes.fasta
```

This will perform the internal sequence homology check and produce the mutant database.

```
Arguments:
--gen_xlsx      : Exported proteins tab from Proteome Discoverer in .xlsx format
--proteome_file : FASTA file containing the proteome, generated from parse_genbank
--gene_file     : FASTA file containing gene sequences, generated from parse_genbank
```

Using Proteome Discoverer, or similar software, perform a search with the mutant database. Do not include any dynamic, N-terminal protein Met-loss modifications. Export peptide groups to an Excel document. (Skip for tutorial).


3) Perform Error Analysis

This depends on the exported peptide groups tab from Proteome Discoverer, or similar software in .xlsx format. The requirement for the input files is that they contain the column headers "Master Protein Accessions", "Sequence", "# PSMs", "PEP."

This will perform filtering, pre-process peptides for chemical artifacts, correct for Ile/Leu, do the analysis, and generate both amino acid and codon substitution spectra. It will also determine the likely position of codon mistranslation and generate a stats file containing data for individual codons (the stats file is only generated automatically when processing a list of files). The resultant files with “processed” in the file name are the data which have been filtered to exclude peptides that contained a chemical artifact, as well as problem homologous peptides.

For processing multiple files via a list:

```
python run_error_analysis.py --mut_fasta ih_mut_custom_wt_ecoli_proteome.fasta --file_list wt_ecoli_file_list.txt --gene_file wt_ecoli_genes.fasta --protein_file wt_ecoli_proteome.fasta --tt 11 --usage abs
```

For processing a single file:

```
python run_error_analysis.py --mut_fasta ih_mut_custom_wt_ecoli_proteome.fasta --workbook_input wt_ecoli-1.xlsx --gene_file wt_ecoli_genes.fasta --protein_file wt_ecoli_proteome.fasta --tt 11 --usage abs
```
```
Arguments:
--mut_fasta     : Mutant FASTA file, made from generate_custom_database
--workbook_input: Exported peptide groups tab from Proteome Discoverer, .xlsx format  (optional, but required if no --file_list is provided)
--file_list     : Text file containing a list of workbook_input files (optional, but required if no --workbook_input is provided)
--protein_file  : FASTA file containing protein sequences from parse_genbank
--gene_file     : FASTA file containing gene sequences from parse_genbank
--tt            : Translation table number (e.g., 11) 
--usage         : Usage type, either 'abs' (absolute) or 'rel' (relative)
--xle_corr      : Optional pickle file to load xle correction file (.pkl)
--nsp_homology  : Optional pickle file to load NSP peptide sequences that are known to be homologous between sequences with different underlying DNA; improves re-analysis run times.
```

Note, non-standard trypsin products (NSPs) searching is time-consuming. Problematic homologous NSPs are filtered during the search. NSPs that are found are saved to the custom mutant FASTA file and will improve search speeds during re-analysis of the data, or when there is an instance of an NSP appearing in multiple samples. So, try to avoid deleting the custom mutant file if you want to save time.
