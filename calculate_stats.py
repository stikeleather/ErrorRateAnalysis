import sys
import getopt
from Bio.Data import CodonTable
import pandas as pd
import numpy as np
from generate_custom_database import fasta_to_dict
import openpyxl
from openpyxl.styles import Font, Border, Side, Alignment

def read_file_names(file_input):
	if isinstance(file_input, list):
		return file_input
	else:
		with open(file_input, 'r') as file:
			file_names = file.read().splitlines()
		return file_names

def calculate_stats(gene_file, file_list, tt, output_file):

	gene_dic = fasta_to_dict(gene_file)

	temp_table = CodonTable.unambiguous_dna_by_id[tt]
	temp_table = dict(temp_table.forward_table.items())

	codon_list = sorted(temp_table.keys())
	aa_list = [temp_table[codon] for codon in codon_list]

	df = pd.DataFrame({
		'AA': aa_list,
		'codon': codon_list})
	
	df['codon-AA'] = df['codon'] + ', ' + df['AA']
	
	df['codon_freq'] = None
	df['RSCU'] = None

	file_count = 0
	error_rate_list = []
	codon_total_list = []
	trans_error_list = []
	file_list = read_file_names(file_list)
	for file in file_list:
		file_count += 1
		temp_df = pd.read_excel(file, sheet_name='Sheet1')
		column_name = f'error_rate{file_count}'
		column_name2 = f'codon_total{file_count}'
		column_name3 = f'translation_error{file_count}'
		error_rate_list.append(column_name)
		codon_total_list.append(column_name2)
		trans_error_list.append(column_name3)

		df.at[0, column_name3] = temp_df.at[0,'translation_error']
		for _, row in temp_df.iterrows():
			codon = str(row['codon']).replace('U', 'T')
			df.loc[df['codon'] == codon, column_name] = row['error_rate']
			df.loc[df['codon'] == codon, column_name2] = row['codon_total']

	for codon in codon_list:
		codon_df = df[df['codon'] == codon]
		mean_error_rate = codon_df[error_rate_list].mean(axis=1)
		std_dev = codon_df[error_rate_list].std(axis=1, ddof=1)
		std_err = std_dev / np.sqrt(file_count)
		mean_codon_total = codon_df[codon_total_list].mean(axis=1)	

		df.loc[df['codon'] == codon, 'mean_error_rate'] = mean_error_rate
		df.loc[df['codon'] == codon, 'std_dev'] = std_dev
		df.loc[df['codon'] == codon, 'Std.Err.'] = std_err
		df.loc[df['codon'] == codon, 'mean_codon_total'] = mean_codon_total

	def calculate_rmr_and_se(row, df):
		mean_error_rate = row['mean_error_rate']
		se_codon = row['Std.Err.']
		aa_family = df[df['AA'] == row['AA']]
		mean_aa_family = aa_family['mean_error_rate'].mean()

		rep_means = []
		for col in error_rate_list:
			mean_for_replicate = aa_family[col].mean()
			rep_means.append(mean_for_replicate)

		se_aa_family = np.std(rep_means, ddof=1) / np.sqrt(len(rep_means))

		if mean_error_rate == 0 or mean_aa_family == 0:
			rmr = 0
			se_rmr = 0
		else:
			rmr = mean_error_rate / mean_aa_family
			se_rmr = rmr * np.sqrt((se_codon / mean_error_rate) ** 2 + (se_aa_family / mean_aa_family) ** 2)

		return pd.Series([rmr, se_rmr])

	df[['RMR', 'SE_RMR']] = df.apply(lambda row: calculate_rmr_and_se(row, df), axis=1)


	codon_table = {
		'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0, 
		'CTT': 0, 'CTC': 0, 'CTA': 0, 'CTG': 0, 'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 
		'ATT': 0, 'ATC': 0, 'ATA': 0, 'ATG': 0, 'AAT': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0, 
		'GTT': 0, 'GTC': 0, 'GTA': 0, 'GTG': 0, 'GAT': 0, 'GAC': 0, 'GAA': 0, 'GAG': 0, 
		'TCT': 0, 'TCC': 0, 'TCA': 0, 'TCG': 0, 'TGT': 0, 'TGC': 0, 'TGA': 0, 'TGG': 0, 
		'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0, 'CGG': 0, 
		'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0, 
		'GCT': 0, 'GCC': 0, 'GCA': 0, 'GCG': 0, 'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}


	for gene in gene_dic:
		sequence = gene_dic[gene]
		for nt in range(0, len(sequence), 3):
			codon = sequence[nt:nt+3]
			codon_table[codon] += 1


	df['codon_freq'] = df['codon'].map(codon_table)

	# Debug
	#print("Codon Frequencies:\n", df[['codon', 'codon_freq']])

	def calculate_rscu(row, df):
		aa_family = df[df['AA'] == row['AA']]
		mean_codon_freq = aa_family['codon_freq'].mean()
		rscu = row['codon_freq'] / mean_codon_freq if mean_codon_freq != 0 else 0
		return rscu

	df['RSCU'] = df.apply(lambda row: calculate_rscu(row, df), axis=1)
	df['codon'] = df['codon'].str.replace('T', 'U', regex=False)

	df['mean_translation_error'] = df[trans_error_list].mean(axis=1)
	df['trans_sdev'] = df[trans_error_list].std(axis=1,ddof=1)
	df['trans_serr'] = df['trans_sdev'] / np.sqrt(file_count)
	
	# Debug
	#print("RSCU Values:\n", df[['codon', 'RSCU']])


	df.to_excel(output_file, sheet_name='Sheet1', index=False)
	wb = openpyxl.load_workbook(output_file)
	ws = wb.active

	custom_font = Font(name='Calibri', size=11, bold=False)
	no_border = Border(left=Side(style=None),
					right=Side(style=None),
					top=Side(style=None),
					bottom=Side(style=None))
	
	alignment = Alignment(horizontal='left', vertical='center')
	
	for cell in ws[1]:
		cell.font = custom_font
		cell.border = no_border
		cell.alignment = alignment

	wb.save(output_file)


def print_usage():
	print("Usage:")
	print("File list should be codon_bias output files")
	print("Example: python calculate_stats.py --gene_file coli_genes.fasta --file_list file_list.txt --tt 11 --output_file codon_and_stats_data.xlsx")
	print("Arguments:")
	print("--gene_file   : Gene sequences in FASTA format")
	print("--file_list   : Text (.txt) file containing codon_bias output files")
	print("--tt          : Translation table (e.g., 11)")
	print("--output_file : Output Excel file for codon and stats data")

if __name__ == '__main__':
	gene_file = None
	file_list = None
	tt = None
	output_file = None

	try:
		options, remainder = getopt.getopt(sys.argv[1:], 'h', ['gene_file=', 'file_list=', 'tt=', 'output_file=', 'help'])
	except getopt.GetoptError:
		print_usage()
		sys.exit(2)

	for opt, arg in options:
		if opt in ('-h', '--help'):
			print_usage()
			sys.exit(0)
		elif opt == '--gene_file':
			gene_file = arg
		elif opt == '--file_list':
			file_list = arg
		elif opt == '--tt':
			tt = int(arg)
		elif opt == '--output_file':
			output_file = arg
		else:
			print(f"Warning! Command-line argument: {opt} not recognized. Exiting...")
			sys.exit(2)

	if not gene_file or not file_list or not tt or not output_file:
		print("Error: Missing one or more required arguments.")
		print_usage()
		sys.exit(2)

	calculate_stats(gene_file, file_list, tt, output_file)
