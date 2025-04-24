import getopt
import pickle
import re
import sys
import openpyxl
import pandas as pd

def mutagenize(peptide):
	AA_dict = {'A':71.04,'R':156.10,'N':114.04,'D':115.03,'C':103.00,
	'E':129.04,'Q':128.06,'G':57.02,'H':137.06,'I':113.08,'L':113.08,
	'K':128.09,'M':131.04,'F':147.07,'P':97.05,'S':87.03,'T':101.04,
	'W':186.08,'Y':163.06,'V':99.07}
	mut_pep_list = list()
	for k in range(len(peptide)):
		for j in AA_dict:
			if j != peptide[k]:
				new_pep = peptide[:k] + j + peptide[(k+1):]
				mut_pep_list.append(new_pep)
	return mut_pep_list

def preprocess_peptides(input_file, analyzed_workbook, homology_file):
	def remove_formatting(input_path):
		original_workbook = openpyxl.load_workbook(filename=input_path)

		new_workbook = openpyxl.Workbook()

		for sheet_name in original_workbook.sheetnames:
			original_sheet = original_workbook[sheet_name]
			new_sheet = new_workbook.create_sheet(title=sheet_name)

			for row in original_sheet.iter_rows(values_only=True):
				new_sheet.append(row)

		new_workbook.remove(new_workbook["Sheet"])
		new_workbook.save(filename=f"{input_path}")

	print("Removing Artifacts...")
	wb0 = openpyxl.load_workbook(input_file)
	ws0 = wb0["Sheet1"]

	df = pd.read_excel(input_file)

	col_names = dict()
	col_count = 0
	for col in range(1,ws0.max_column+1):
		if ws0.cell(1,col).value is not None:
			col_count += 1
			col_head = ws0.cell(1,col).value
			col_names.update({col_head:col_count})
		else:
			break

	seq = col_names['Sequence']
	psm = col_names['# PSMs']

	row_num0 = 0
	for www in range(1,ws0.max_row+1):
		if ws0.cell(www,psm).value is not None:
			row_num0 += 1
		else:
			break
	row_num0 = row_num0+1
	# Artifact regex patterns
	NtoD = re.compile(r'N\d+D')
	QtoE = re.compile(r'Q\d+E')
	EtoS = re.compile(r'E\d+S')
	StoD = re.compile(r'S\d+D')
	TtoE = re.compile(r'T\d+E')
	StoA = re.compile(r'S\d+A')
	YtoF = re.compile(r'Y\d+F')
	# Artifact Definitions
	ptm_definitions = {NtoD,QtoE,EtoS,StoD,TtoE,StoA,YtoF}

	wb = openpyxl.load_workbook(analyzed_workbook)
	ws = wb["Sheet1"]

	row_num = 0
	for www in range(1,ws.max_row+1):
		if ws.cell(www,1).value is not None:
			row_num += 1
		else:
			break
	row_num = row_num+1

	target_seq_set = set()
	for acc in range(2,row_num):
		if 'mutant-' in ws.cell(acc,1).value:
			accession = ws.cell(acc,1).value
			split_acc = re.split('[.]',accession)
			mutation = split_acc[2]
			for pattern in ptm_definitions:
				if pattern.search(mutation):
					target_seq = ws.cell(acc,2).value
					target_seq_set.add(target_seq)
					
	with open(homology_file, 'rb') as f:
		homolgy_peptides = pickle.load(f)

	for row in range(2,row_num):
		test = ws0.cell(row,seq).value
		if test in homolgy_peptides:
			target_seq_set.add(test)

	for item in homolgy_peptides:
		mut_peptides = mutagenize(item)
		for entry in mut_peptides:
			target_seq_set.add(entry)


	data_rows = []
	for index, row in df.iterrows():
		if pd.notnull(row['Sequence']) and row['Sequence'] != "":
			if row['Sequence'] not in target_seq_set:
				data_rows.append(row)

	processed_data = pd.DataFrame(data_rows, columns=df.columns)


	out_file = f"processed_{input_file}"
	processed_data.to_excel(out_file,index=False)
	remove_formatting(out_file)
	return out_file
	
if __name__ == '__main__':
	try:
		options, remainder = getopt.getopt(sys.argv[1:],'', ['input_file=','analyzed_workbook=', 'homology_file='])
	except getopt.GetoptError:
		print("Example: python preprocess_peptides.py --input_file wt_ecoli-1_pep99.xlsx --analyzed_workbook analyzed_wt_ecoli-1_pep99.xlsx --homology_file homology_ih_mut_custom_wt_ecoli_proteome.pkl")

	for opt, arg in options:
		if opt == '--input_file': 
			input_file=arg
		elif opt == '--analyzed_workbook': 
			analyzed_workbook=arg
		elif opt == '--homology_file': 
			homology_file = arg
		else: 
			print("Warning! Command-line argument: %s not recognized. Exiting..." % opt); sys.exit(2)

	if remainder:
		print(f"Unrecognized arguments: {remainder}. These were ignored.")

	preprocess_peptides(input_file, analyzed_workbook, homology_file)