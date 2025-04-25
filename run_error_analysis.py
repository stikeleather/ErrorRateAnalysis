import getopt
import re
import sys
import time
import openpyxl
from Bio import SeqIO
import pickle
from matrix_aa_subs import matrix_aa_subs
from generate_custom_database import fasta_to_dict
from preprocess_peptides import preprocess_peptides
import pandas as pd
from openpyxl.styles import Font, Border, Side, Alignment
from codon_bias import codon_bias
from calculate_stats import calculate_stats
from format_matrix import format_matrix
from codon_matrix import codon_matrix


def find_match(long_pep, short_pep):
	try:
		start_index = long_pep.index(short_pep)
		end_index = start_index + len(short_pep)
		return start_index, end_index
	except ValueError:
		return -1, -1
	
def determine_sub_pos(long_pep, short_pep, mutation):
	start_range = long_pep.find(short_pep)
	pos = int(mutation.group(2))
	sub_pos = pos - start_range
	new_seq = short_pep[:sub_pos-1] + mutation.group(1) + short_pep[sub_pos:]
	return new_seq, sub_pos

def count_mismatches(test, comp):
	if len(test) != len(comp):
		raise ValueError("Test and Comp must be equal in length")
	mismatch_count = 0
	for aa1, aa2 in zip(test, comp):
		if aa1 != aa2:
			mismatch_count += 1
	return mismatch_count

def read_file_names(file_input):
	if isinstance(file_input, list):
		return file_input
	else:
		with open(file_input, 'r') as file:
			file_names = file.read().splitlines()
		return file_names

def extract_species_name(filename):
	species_pattern = r'ih_mut_custom_(.*?)_proteome\.fasta'
	match = re.search(species_pattern, filename)
	if match:
		return match.group(1)
	else:
		return None
	
def filter_pep(input_file):
	data = pd.read_excel(input_file)
	data = data[data["PEP"] < 0.01]
	file_name = re.split("[.]",input_file)
	file_name = file_name[0]
	output_file = f"{file_name}_pep99.xlsx"
	data.to_excel(output_file, index=False)
	return output_file

def run_error_analysis(mut_fasta, workbook_input, load_xle=None):
	#Removes peptides with PEP scores greater than 0.01
	workbook_out = f'analyzed_{workbook_input}'
	#Output xlsx creation
	wb = openpyxl.Workbook()
	wb.create_sheet(index=0, title='Sheet1')
	wb.remove(wb["Sheet"])
	wb.save(workbook_out)
	#Set worksheet-0 to the raw data
	#Raw data is never saved here, so any edits are temporary; integrity of raw data is secure
	wb0 = openpyxl.load_workbook(workbook_input)
	ws0 = wb0["Sheet1"]

	#Determines the number of rows; avoids calling ws0.max_row()
	row_num = 0
	for www in range(1,ws0.max_row+2):
		if ws0.cell(www,1).value is not None:
			row_num += 1
		else:
			break
	row_num = row_num+1

	#Dictionary of column names; makes it so order of columns does not matter
	col_names = dict()
	col_count = 0
	for col in range(1,ws0.max_column+1):
		if ws0.cell(1,col).value is not None:
			col_count += 1
			col_head = ws0.cell(1,col).value
			col_names.update({col_head:col_count})
		else:
			break
	
	#Sets column integars by column name
	mpa = col_names['Master Protein Accessions']
	seq_col = col_names['Sequence']
	psm = col_names['# PSMs']				

	#If empty cells are present in the accession column, add NA-count string; avoids loop terminations on the acc column due to empty cells
	na_count = 0
	total_count = 0
	for v in range(2, row_num):
		total_count += 1
		if ws0.cell(v,mpa).value is None:
			na_count += 1
			ws0.cell(v,mpa).value = f"NA-{na_count}"
	if na_count > 0:
		print("Added NA to empty cells")

	#Dictionary of peptides in the data
	search_dic = dict()
	for row in range(2,row_num):
		if ws0.cell(row,seq_col).value is not None:
			search_dic.update({row:ws0.cell(row,seq_col).value})
		else:
			break

	#Intiate Regular Expression patterns
	ItoL = re.compile(r'I\d+L')
	LtoI = re.compile(r'L\d+I')
	# The $ anchors mut_patten to the end of the string, avoiding issues with locus ids
	mut_pattern = re.compile(r'([A-Z])(\d+)([A-Z])$')
	dash_pattern = re.compile(r'[-]')
	dot_pattern = re.compile(r'[.]')
	dash_dot_pattern = re.compile(r'[-.]')
	header_pattern = re.compile(r'mutant-([^\.]+)')
	xle_correction = dict()

	#Sets homology file
	split_mut_fasta = re.split(dot_pattern,mut_fasta)
	homology_file = f"homology_{split_mut_fasta[0]}.pkl"

	#Reads the mut_fasta; creates dictionary of headers:seqs
	print("Reading input file...")
	mutant_fasta = fasta_to_dict(mut_fasta)
	wild_seqs_set = set()
	mutant_seqs_set = set()
	print("Sorting peptide database...")
	for key in mutant_fasta:
		split_key = re.split(dash_pattern,key)
		if split_key[0] == 'wild':
			wild_seqs_set.add(mutant_fasta[key])
		else:
			mutant_seqs_set.add(mutant_fasta[key])

	print("Correcting for Xle")
	#Loops through all possible Xle mutations; compares to known sequences in data; reverts to wt seq if Xle sub is found;
	#Stores raw_seq:xle_correction_seq in dict; will dump dic to pkl for future loading (makes data reruns faster); not suitable between samples
	xle_dump = False
	if not load_xle:
		for key in mutant_fasta:
			#Looks for Xle pattern in headers
			if re.search(ItoL, key) or re.search(LtoI, key):
				target_seq = mutant_fasta[key]
				for entry in search_dic:
					if search_dic[entry] == target_seq:
						header_match = re.search(header_pattern,key)
						header_match = header_match.group(1)
						wild_header = f"wild-{header_match}"
						wild_seq = mutant_fasta[wild_header]
						# print(f"{key}")
						# print(target_seq)
						# print(wild_header)
						# print(wild_seq)
						ws0.cell(entry,seq_col).value = wild_seq
						xle_correction[target_seq] = wild_seq
						xle_dump = True
		
		split_workbook = re.split(dot_pattern,workbook_input)
		workbook_name = split_workbook[0]
		xle_file_name = f'xle_{workbook_name}.pkl'

	else:
		xle_file_name = load_xle
		with open(load_xle,'rb') as f:
			xle_correction = pickle.load(f)
			print("Loaded Xle")
			for entry in search_dic:
				if search_dic[entry] in xle_correction.keys():
					raw_seq = search_dic[entry]
					ws0.cell(entry,seq_col).value = xle_correction[raw_seq]

	#Assigning output file column names; order is important
	wb = openpyxl.load_workbook(workbook_out)
	ws = wb["Sheet1"]
	column_names = ['accession', 'Sequence', 'PSM_count', 'pep_length', 'wild_aa_total', 'mutant_aa_total', 'total_aa', 'error_rate']
	for i in range(1,9): 
		ws.cell(1,i).value = column_names[i-1]

	#Dictionary where dic[header] = seq is flipped to dic[seq] = header
	print("Building Dictionaries")
	name_seq_dic = dict()
	for name, seq in mutant_fasta.items():
		name_seq_dic[seq] = name


	#First rebuild of the search_dic (post xle correction)
	print("Building Search Targets")
	search_dic = dict()
	for row in range(2,row_num):
		search_dic.update({row:ws0.cell(row,seq_col).value})	
	#Set of observed sequences
	observed_sequences = set(search_dic.values())
	final_dic_len = (len(observed_sequences))


	# name_seq_dic.keys() are peptide sequences; creates a set of peptide sequences for intersection with set of observed sequences
	#Non-standard peptide products will miss this intersection
	unique_name_seq_terms = set(name_seq_dic.keys())
	exist_set = observed_sequences.intersection(unique_name_seq_terms)


	#Compares how many sequences intersected to the number of observed sequences
	print(f"Matched Dictionary Entry {len(exist_set)} of {final_dic_len}")

	#Determines how many non-standard sequences are in the observed data
	total_unknowns = 0
	test_terms = set()
	for i in range(2,row_num):
		test = ws0.cell(i,seq_col).value
		if test not in exist_set:
			test_terms.add(test)
	total_unknowns = len(test_terms)

	#Determines nsp id starting value to ensure previous nsp ids are not overwritten
	nsp_start_number = 0
	for key in mutant_fasta:
		if "-nsp" in key:
			nsp_start_number += 1

	if total_unknowns > 0:
		print(f"{total_unknowns} Non-standard products detected")
		unknown_count = 0
		total_uknown_wilds = 0
		total_uknown_muts = 0
		#nsp_dic will be appended to mut_fasta for future data analysis with this mut_fasta (is suitable between samples)
		nsp_dic = dict()
		#First checks to see if the test seq is found in any wild-type peptide (comparison peptide)
		#Only comparison peptides used are those found to be in the data
		#If any matches are found then peptide must be wild and the test is removed from test_terms 
		#note, we are looping over a copy() of test_terms and discarding from actual test terms
		#Each NSP needs a unique key
		print("Testing for wild NSPs...")
		for test in test_terms.copy():
			for comp in wild_seqs_set:
				if test in comp and "-nsp" not in name_seq_dic[comp]:
					unknown_count += 1
					total_uknown_wilds += 1
					print(f"Matched {unknown_count} of {total_unknowns} NSPs")
					nsp_dic[f">{name_seq_dic[comp]}-nsp{nsp_start_number+unknown_count}"] = f"{test}"
					mutant_fasta[f"{name_seq_dic[comp]}-nsp{nsp_start_number+unknown_count}"] = f"{test}"
					test_terms.discard(test)
					break
		#If test terms still exist after searching for wilds, then search for mutants
		print("Testing for mutant NSPs...")
		if len(test_terms) > 0:
			for unk in test_terms.copy():
				for comp in mutant_seqs_set:
					if unk in comp and "-nsp" not in name_seq_dic[comp]:
						unknown_count += 1
						print(f"Matched {unknown_count} of {total_unknowns} NSPs")
						#If the nsp is Xle, then correct it to wild-type and determine sub position
						if re.search(ItoL, name_seq_dic[comp]) or re.search(LtoI, name_seq_dic[comp]):
							total_uknown_wilds += 1
							#If comp-unk pair is mapped to an NSP, then adjusts for mut_pattern anchoring
							if 'nsp' not in name_seq_dic[comp]:
								mutation = re.search(mut_pattern,name_seq_dic[comp])
							else:
								split_name_comp = re.split(dash_dot_pattern,name_seq_dic[comp])
								mutation = re.search(mut_pattern,split_name_comp[4])
							split_name = re.split(dash_dot_pattern, name_seq_dic[comp])
							correction = determine_sub_pos(comp,unk,mutation)
							new_seq = correction[0]
							sub_pos = correction[1]
							#Replaces an Xle sequence in the raw data with the corrected wt seq
							#The raw data workbook is never saved after editing; any edits are not permanent
							for entry in search_dic:
								if unk == search_dic[entry]:
									ws0.cell(entry,seq_col).value = new_seq					
							xle_correction[unk] = new_seq
							xle_dump = True
							nsp_dic[f">mutant-{split_name[1]}-{split_name[2]}.{split_name[3]}.{mutation.group(1)}{sub_pos}{mutation.group(3)}-nsp{nsp_start_number+unknown_count}"] = f"{unk}"		
							#nsp_dic[f">wild-{split_name[1]}-{split_name[2]}-nsp{unknown_count}"] = f"{new_seq}"
							mutant_fasta[f"wild-{split_name[1]}-{split_name[2]}-nsp{nsp_start_number+unknown_count}"] = f"{new_seq}"
							test_terms.discard(unk)
							break
						else:
							#If nsp is not xle, then determine sub position for correct naming
							#Correct naming is important for usage of NSPs in downstream codon_bias
							total_uknown_muts += 1
							if '-nsp' not in name_seq_dic[comp]:
								mutation = re.search(mut_pattern,name_seq_dic[comp])
							else:
								split_name_comp = re.split(dash_dot_pattern,name_seq_dic[comp])
								mutation = re.search(mut_pattern,split_name_comp[4])								
							correction = determine_sub_pos(comp,unk,mutation)
							sub_pos = correction[1]
							split_name = re.split(dash_dot_pattern, name_seq_dic[comp])
							nsp_dic[f">mutant-{split_name[1]}-{split_name[2]}.{split_name[3]}.{mutation.group(1)}{sub_pos}{mutation.group(3)}-nsp{nsp_start_number+unknown_count}"] = f"{unk}"
							mutant_fasta[f"mutant-{split_name[1]}-{split_name[2]}.{split_name[3]}.{mutation.group(1)}{sub_pos}{mutation.group(3)}-nsp{nsp_start_number+unknown_count}"] = f"{unk}"
							test_terms.discard(unk)
							break

		#Determines likely homology peptides that appeared in the data
		homology_discard = False
		if len(test_terms) > 0:
			discarded_homology = set()
			with open(homology_file, 'rb') as f:
				homolgy_peptides = pickle.load(f)
			for test in test_terms.copy():
				test_len = len(test)
				for comp in homolgy_peptides:
					if test_len == len(comp):
						num_mismatches = count_mismatches(test,comp)
						if num_mismatches == 1 or num_mismatches == 0:
							test_terms.discard(test)
							discarded_homology.add(test)
							homology_discard = True
							break
			if homology_discard:
				print("The following peptides were ignored due to internal homology:")
				print(discarded_homology)
		#Anything still remaining
		if len(test_terms) > 0:	
			print("Unmatched NSPs:")
			print(test_terms)

		print(f"Found {total_unknowns} non-standard products")
		print(f"{total_uknown_wilds} non-standard wilds")
		print(f"{total_uknown_muts} non-standard muts")


		if len(nsp_dic) > 0:
			print("Updating Database")
			#Stores unknowns for future runs with same mut_fasta, stored in mut_file instead of directly editing raw_data
			with open(mut_fasta, 'a') as psm_input:
				for key in nsp_dic:
					psm_input.write(f"{key}\n")
					psm_input.write(f"{nsp_dic[key]}\n")


		print("Re-Building Search Targets")
		# Collects raw peptide sequences from data input; including updated xle seqs stored in the copy of the xlsx
		search_dic = dict()
		for zz in range(2,row_num):
			search_dic.update({zz:ws0.cell(zz,seq_col).value})

		#Determines number of unique sequences, and creates set of sequences in raw data
		observed_sequences = set(search_dic.values())
		final_dic_len = (len(observed_sequences))

		#Dictionary where sequences are keys, and names are values
		name_seq_dic = {}
		for name, seq in mutant_fasta.items():
			name_seq_dic[seq] = name

		#Creates a set of all sequences to hashtable intersect observed peptide sequences
		#Results in a set of sequences
		unique_name_seq_terms = set(name_seq_dic.keys())
		exist_set = observed_sequences.intersection(unique_name_seq_terms)

		#Each sequence should be matched to a name
		psm_dic = dict()
		for term in exist_set:
			psm_dic.update({term:name_seq_dic[term]})
		print(f"Matched Dictionary Entry {len(exist_set)} of {final_dic_len}")
	else:
		psm_dic = dict()
		for term in exist_set:
			psm_dic.update({term:name_seq_dic[term]})

	print("Matching Names to PSMs")
	row_count2 = 2
	running_count = 0
	match_count = 0
	for key in psm_dic:
		running_count = 0
		for z in range(2,row_num):
			#print(ws0.cell(z,col_names['Sequence']).value)
			if key == ws0.cell(z,col_names['Sequence']).value:
				ws0.cell(z,mpa).value = psm_dic[key]
				running_count += ws0.cell(z,psm).value
				match_count += 1
				
		if running_count > 0:
			ws.cell(row_count2,1).value = psm_dic[key]
			ws.cell(row_count2,2).value = key
			ws.cell(row_count2, 3).value = running_count
			ws.cell(row_count2, 4).value = len(key)
			row_count2 += 1
			running_count = 0

	print(f"Matched PSM {match_count} of {total_count} to Name")
	
	#Math
	wb.save(workbook_out)
	wild_aa = 0
	mut_aa = 0
	total_aa = 0
	for row2 in range(2,ws.max_row+1):
		if ws.cell(row2,1).value is not None:
			acc = ws.cell(row2,1).value
			split_acc = re.split(dash_pattern, acc)
			if split_acc[0] == 'wild':
				psm_value = ws.cell(row2,3).value
				pep_len = ws.cell(row2,4).value
				wild_aa_result = psm_value * pep_len
				wild_aa += wild_aa_result
				total_aa += wild_aa_result
			else:
				psm_value = ws.cell(row2,3).value
				pep_len = ws.cell(row2,4).value
				wild_aa += (psm_value * pep_len) - psm_value
				mut_aa += psm_value
				total_aa += psm_value * pep_len
		else:
			break

	ws.cell(2,5).value = wild_aa
	ws.cell(2,6).value = mut_aa
	ws.cell(2,7).value = total_aa
	ws.cell(2,8).value = mut_aa/total_aa

	#Dumps the xle corrections made; can be used between runs of the same data input
	if xle_dump:
		print("Writing Xle correction...")
		with open(xle_file_name,'wb') as q:
			pickle.dump(xle_correction,q)
	#Final workbook save
	wb.save(workbook_out)

	print("Generating Matrix...")
	matrix_out = matrix_aa_subs(workbook_out)

	return workbook_input, workbook_out, xle_file_name, homology_file, matrix_out
	


def print_usage():
	print("Usage:")
	print("Example: python run_error_analysis.py --mut_fasta ih_mut_custom_example_proteome.fasta --workbook_input example.xlsx --load_xle xle.pkl --gene_file genes.fasta --protein_file proteome.fasta --tt 11 --usage abs")
	print("Minimum 1: python run_error_analysis.py --mut_fasta ih_mut_custom_example_proteome.fasta --workbook_input example.xlsx --gene_file genes.fasta --protein_file proteome.fasta --tt 11 --usage abs")
	print("Minimum 2: python run_error_analysis.py --mut_fasta ih_mut_custom_example_proteome.fasta --file_list example_file_list.txt --gene_file genes.fasta --protein_file proteome.fasta --tt 11 --usage abs")
	print("\nArguments:")
	print("--mut_fasta     : Mutant FASTA file, made from generate_custom_database")
	print("--workbook_input: Exported peptide groups tab from Proteome Discoverer, .xlsx format  (optional, but required if no --file_list is provided)")
	print("--load_xle      : Optional pickle file to load xle correction file (.pkl)")
	print("--file_list     : Text file containing a list of workbook_input files (optional, but required if no --workbook_input is provided)")
	print("--gene_file     : FASTA file containing gene sequences from parse_genbank")
	print("--protein_file  : FASTA file containing protein sequences from parse_genbank")
	print("--tt            : Translation table number (e.g., 11) ")
	print("--usage         : Usage type, either 'abs' (absolute) or 'rel' (relative)")

if __name__ == '__main__':
	load_xle = None
	file_list = None
	combined_data = pd.DataFrame()
	mut_fasta = None
	workbook_input = None
	gene_file = None
	protein_file = None
	usage = None
	tt = None

	try:
		options, remainder = getopt.getopt(sys.argv[1:], 'h', ['mut_fasta=', 'workbook_input=', 'load_xle=', 'file_list=', 'gene_file=', 'protein_file=', 'tt=', 'usage=', 'help'])
	except getopt.GetoptError:
		print_usage()
		sys.exit(2)

	for opt, arg in options:
		if opt in ('-h', '--help'):
			print_usage()
			sys.exit(0)
		elif opt == '--mut_fasta': 
			mut_fasta = arg
		elif opt == '--workbook_input': 
			workbook_input = str(arg)
		elif opt == '--load_xle': 
			load_xle = arg
		elif opt == '--file_list':
			file_list = arg
		elif opt == '--gene_file':
			gene_file = arg
		elif opt == '--protein_file':
			protein_file = arg
		elif opt == '--tt':
			tt = int(arg)
		elif opt == '--usage':
			usage = arg
		else: 
			print(f"Warning! Command-line argument: {opt} not recognized. Exiting...")
			sys.exit(2)

	if remainder:
		print(f"Unrecognized arguments: {remainder}. These were ignored.")

	if usage not in ['abs', 'rel']:
		print("Error: The 'usage' argument must be either 'abs' or 'rel'.")
		sys.exit(2)

	if not mut_fasta:
		print("Error: Missing mut_fasta argument")
		sys.exit(2)
	if not workbook_input and not file_list:
		print("Error: Missing workbook_input or file_list argument")
		sys.exit(2)
	if not gene_file:
		print("Error: Missing gene_file argument")
		sys.exit(2)
	if not protein_file:
		print("Error: Missing protein_file argument")
		sys.exit(2)
	if not tt:
		print("Error: Missing translation table (tt) argument")
		sys.exit(2)

	if file_list:
		workbook_inputs = read_file_names(file_list)
		filtered_file_list = []
		codon_file_list = []
	else:
		workbook_inputs = [workbook_input]

	for workbook_input in workbook_inputs:
		start_time = time.time()
		workbook_input = filter_pep(workbook_input)
		out_items = run_error_analysis(mut_fasta, workbook_input, load_xle)
		workbook_input = out_items[0]
		workbook_out = out_items[1]
		xle_file_name = out_items[2]
		homology_file = out_items[3]
		processed_file = preprocess_peptides(workbook_input, workbook_out, homology_file)
		if file_list:
			filtered_file_list.append(processed_file)
		out_items = run_error_analysis(mut_fasta, processed_file, xle_file_name)
		workbook_out = out_items[1]
		format_matrix(out_items[4], remove_artifacts=True)
		codon_output_file = codon_bias(workbook_out, gene_file, protein_file, tt, usage)
		codon_matrix(workbook_out, gene_file, protein_file, tt)
		if file_list:
			codon_file_list.append(codon_output_file)

		print("Finished")
		end_time = time.time()
		total_time = end_time - start_time
		if total_time < 60:
			print(f"Total time: {total_time:.2f} seconds")
		else:
			minutes = int(total_time // 60)
			seconds = int(total_time % 60)
			print(f"Total time: {minutes} minute(s) {seconds} second(s)")


	if file_list:
		species_name = extract_species_name(mut_fasta)
		stats_output_file = f'{species_name}_stats.xlsx'
		calculate_stats(gene_file, codon_file_list, tt, stats_output_file)
		combined_data_list = []
		for file in filtered_file_list:
			df = pd.read_excel(file)
			for index, row in df.iterrows():
				combined_data_list.append(row)
		combined_data = pd.DataFrame(combined_data_list, columns=df.columns)

		species_name = extract_species_name(mut_fasta)
		combined_file_name = f'processed_{species_name}_all_pep99.xlsx'
		combined_data.to_excel(combined_file_name, index=False)
		wb = openpyxl.load_workbook(combined_file_name)
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

		wb.save(combined_file_name)
		load_xle = None
		out_items = run_error_analysis(mut_fasta, combined_file_name, load_xle)
		workbook_out = out_items[1]
		format_matrix(out_items[4], remove_artifacts=True)
		codon_matrix(workbook_out, gene_file, protein_file, tt)
		codon_bias(workbook_out, gene_file, protein_file, tt, usage)
		
		print("Finished")
		end_time = time.time()
		total_time = end_time - start_time
		if total_time < 60:
			print(f"Total time: {total_time:.2f} seconds")
		else:
			minutes = int(total_time // 60)
			seconds = int(total_time % 60)
			print(f"Total time: {minutes} minute(s) {seconds} second(s)")