import sys
import getopt
import pandas as pd
from Bio import SeqIO
import re
from tryp_and_mutate import mutate

def fasta_to_dict(fasta_file):
	sequence_dict = {}
	for record in SeqIO.parse(fasta_file, "fasta"):
		if record.id in sequence_dict.keys():
			print(f"Duplicate gene detected in fasta: {record.id}")
		else:
			sequence_dict[record.id] = str(record.seq)
	return sequence_dict

def generate_custom_database(gen_xlsx,proteome_file,gene_file):
	df = pd.read_excel(gen_xlsx)

	proteome_dic = fasta_to_dict(proteome_file)
	gene_dic = fasta_to_dict(gene_file)
	
	accessions = list()
	try:
		for accession in df["Accession"]:
			accessions.append(accession)
	except KeyError:
		for accession in df["Master Protein Accession"]:
			accessions.append(accession)
	accession_count = 0
	accessions = set(accessions)
	accession_number = len(accessions)

	error_count = 0
	unmatched_list = list()

	split_proteome = re.split('[.]',proteome_file)
	protein_output = f"custom_{split_proteome[0]}.fasta"

	split_gene = re.split('[.]',gene_file)
	gene_output = f"custom_{split_gene[0]}.fasta"

	pattern = r"^(.*?)(?=_proteome|$)"
	species = re.match(pattern, split_proteome[0])
	species = species.group(1)

	with open(protein_output, "w") as f, open(gene_output, 'w') as q:
		for accession in accessions:
			if ";" in accession:
				split_accession = re.split('[;]', accession)
				for i in range(0,len(split_accession)):
					test_key = split_accession[i].strip()
					try:
						pro_seq = proteome_dic[test_key]
						gene_seq = gene_dic[test_key]
						f.write(f">{test_key}\n")
						f.write(f"{pro_seq}\n")
						q.write(f">{test_key}\n")
						q.write(f"{gene_seq}\n")
						accession_count += 1
						error_count = 0
						break
					except KeyError:
						error_count += 1
						if error_count == len(split_accession):
							unmatched_list.append(accession)
						continue
			else:
				try:
					pro_seq = proteome_dic[accession]
					gene_seq = gene_dic[accession]
					f.write(f">{accession}\n")
					f.write(f"{pro_seq}\n")
					q.write(f">{accession}\n")
					q.write(f"{gene_seq}\n")
					accession_count += 1
				except KeyError:
					print(f"Missing ID: {accession}")

	if len(unmatched_list) > 0:
		with open("unmatched.txt", 'w') as f:
			for item in unmatched_list:
				f.write(f"{item}\n")

	print(f"Matched {accession_count} of {accession_number} entries")

	if accession_count == accession_number and len(unmatched_list) == 0:
		output_file = mutate(protein_output,2,1,gene_output)
		print("Done")
		print("For a list of files:")
		print(f"Please execute: python run_error_analysis.py --mut_fasta {output_file} --file_list {species}_file_list.txt --gene_file {species}_genes.fasta --protein_file {species}_proteome.fasta --tt <1-33> --usage abs")
		print(f"\nFor a single file:")
		print(f"Please execute: python run_error_analysis.py --mut_fasta {output_file} --workbook_input example.xlsx --gene_file {species}_genes.fasta --protein_file {species}_proteome.fasta --tt <1-33> --usage abs")
	else:
		print("Unmatched sequences detected, resolve and execute again to generate mutant database")


def print_usage():
    print("Usage:")
    print("Example: python generate_custom_database.py --gen_xlsx wt_ecoli_generic.xlsx --proteome_file wt_ecoli_proteome.fasta --gene_file wt_ecoli_genes.fasta")
    print("\nArguments:")
    print("--gen_xlsx      : Exported proteins tab from Proteome Discoverer in .xlsx format")
    print("--proteome_file : FASTA file containing the proteome, generated from parse_genbank")
    print("--gene_file     : FASTA file containing gene sequences, generated from parse_genbank")

if __name__ == '__main__':
    gen_xlsx = None
    proteome_file = None
    gene_file = None

    try:
        options, remainder = getopt.getopt(sys.argv[1:], 'h', ['gen_xlsx=', 'proteome_file=', 'gene_file=', 'help'])
    except getopt.GetoptError:
        print_usage()
        sys.exit(2)

    for opt, arg in options:
        if opt in ('-h', '--help'):
            print_usage()
            sys.exit(0)
        elif opt == '--gen_xlsx':
            gen_xlsx = arg
        elif opt == '--proteome_file':
            proteome_file = arg
        elif opt == '--gene_file':
            gene_file = arg
        else:
            print(f"Warning! Command-line argument: {opt} not recognized. Exiting...")
            sys.exit(2)

    if remainder:
        print(f"Unrecognized arguments: {remainder}. These were ignored.")

    if not gen_xlsx:
        print("Error: Missing --gen_xlsx argument (Excel file)")
        print_usage()
        sys.exit(2)
    if not proteome_file:
        print("Error: Missing --proteome_file argument (Proteome FASTA file)")
        print_usage()
        sys.exit(2)
    if not gene_file:
        print("Error: Missing --gene_file argument (Gene FASTA file)")
        print_usage()
        sys.exit(2)

    generate_custom_database(gen_xlsx, proteome_file, gene_file)